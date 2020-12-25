#include <assert.h>
#include <math.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_dg_maxwell.h>
#include <gkyl_dg_vlasov.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_mom_calc.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_vlasov.h>
#include <gkyl_vlasov_mom.h>

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(sizeof(double)*nc, size);
  gkyl_array_clear(a, 0.0);
  return a;
}

// list of valid moment names
static const char *const valid_moment_names[] = { "M0", "M1i", "M2ij", "M2", "M3i" };

// check if name of moment is valid or not
static int
is_moment_name_valid(const char *nm)
{
  int n = sizeof(valid_moment_names)/sizeof(valid_moment_names[0]);
  for (int i=0; i<n; ++i)
    if (strcmp(valid_moment_names[i], nm) == 0)
      return 1;
  return 0;
}

// initialize local and local-ext ranges on conf-space.
static void
init_conf_ranges(int cdim, const int *cells,
  struct gkyl_range *local, struct gkyl_range *local_ext)
{
  int lower_ext[GKYL_MAX_DIM], upper_ext[GKYL_MAX_DIM];
  int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  
  for (int i=0; i<cdim; ++i) {
    lower_ext[i] = -1;
    upper_ext[i] = cells[i];

    lower[i] = 0;
    upper[i] = cells[i]-1;
  }
  gkyl_range_init(local, cdim, lower, upper);    
  gkyl_range_init(local_ext, cdim, lower_ext, upper_ext);
}

// initialize local and local-ext ranges on phase-space: each species
// can have different phase-space grid and hence ranges; note that
// this function assumes there are no ghost-cells in velocity space.
static void
init_phase_ranges(int cdim, int pdim, const int *cells,
  struct gkyl_range *local, struct gkyl_range *local_ext)
{
  int lower_ext[GKYL_MAX_DIM], upper_ext[GKYL_MAX_DIM];
  int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  
  for (int i=0; i<cdim; ++i) {
    lower_ext[i] = -1;
    upper_ext[i] = cells[i];

    lower[i] = 0;
    upper[i] = cells[i]-1;
  }
  for (int i=cdim; i<pdim; ++i) {
    lower_ext[i] = 0;
    upper_ext[i] = cells[i]-1;

    lower[i] = 0;
    upper[i] = cells[i]-1;
  }
  gkyl_range_init(local, pdim, lower, upper);
  gkyl_range_init(local_ext, pdim, lower_ext, upper_ext);
}

// data for moments
struct vm_species_moment {
    struct gkyl_mom_type *mtype;
    gkyl_mom_calc *mcalc;
    struct gkyl_array *marr;
};

// species data
struct vm_species {
    struct gkyl_vlasov_species info; // data for species
    
    struct gkyl_rect_grid grid;
    struct gkyl_range local, local_ext; // local, local-ext phase-space ranges

    // arrays for distribution function updates
    struct gkyl_array *f;

    struct vm_species_moment m1i; // for computing currents
    struct vm_species_moment *moms; // diagnostic moments
};

// field data
struct vm_field {
    struct gkyl_em_field info; // data for field
    
    // arrays for field updates
    struct gkyl_array *em;
    
    struct gkyl_dg_eqn *eqn; // Maxwell equation
    struct gkyl_hyper_dg *slvr; // solver 
};

// Vlasov object: used as opaque pointer in user code
struct gkyl_vlasov_app {
    char name[128]; // Name of app
    int cdim, vdim; // Conf, velocity space dimensions
    int poly_order; // Polynomial order
    
    struct gkyl_rect_grid grid; // config-space grid
    struct gkyl_range local, local_ext; // local, local-ext conf-space ranges
    struct gkyl_basis basis, confBasis; // phase-space, conf-space basis

    // field data
    struct vm_field field;

    // species data
    int num_species;
    struct vm_species *species;
};

// initialize field object
static void
vm_species_moment_init(struct gkyl_vlasov_app *app, struct vm_species *s,
  struct vm_species_moment *sm, const char *nm)
{
  assert(is_moment_name_valid(nm));
  
  sm->mtype = gkyl_vlasov_mom_new(&app->confBasis, &app->basis, nm);
  sm->mcalc = gkyl_mom_calc_new(&s->grid, sm->mtype);
  sm->marr = mkarr(sm->mtype->num_mom*app->confBasis.numBasis,
    app->local_ext.volume);
}

// release memory for moment data object
static void
vm_species_moment_release(struct vm_species_moment *sm)
{
  gkyl_mom_type_release(sm->mtype);
  gkyl_mom_calc_release(sm->mcalc);
  gkyl_array_release(sm->marr);
}

// initialize field object
static void
vm_field_init(struct gkyl_vm *vm, struct gkyl_vlasov_app *app, struct vm_field *f)
{
  // allocate EM arrays
  f->em = mkarr(8*app->confBasis.numBasis, app->local_ext.volume);

  // equation object
  double c = 1/sqrt(f->info.epsilon0*f->info.mu0);
  double ef = f->info.elcErrorSpeedFactor, mf = f->info.mgnErrorSpeedFactor;
  f->eqn = gkyl_dg_maxwell_new(&app->confBasis, c, ef, mf);

  int max_up_dirs[] = {0, 1, 2}, zero_flux_flags[] = {0, 0, 0};
  // Maxwell solver
  f->slvr = gkyl_hyper_dg_new(&app->grid, &app->confBasis, f->eqn,
    app->cdim, max_up_dirs, zero_flux_flags, 1);
}

// release resources for field
static void
vm_field_release(struct vm_field *f)
{
  gkyl_array_release(f->em);
  
  gkyl_dg_eqn_release(f->eqn);
  gkyl_hyper_dg_release(f->slvr);
}

// initialize species object
static void
vm_species_init(struct gkyl_vm *vm, struct gkyl_vlasov_app *app, struct vm_species *s)
{
  int cdim = app->cdim, vdim = app->vdim, poly_order = app->poly_order;
  int pdim = cdim+vdim;

  int cells[GKYL_MAX_DIM];
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];

  for (int d=0; d<cdim; ++d) {
    cells[d] = vm->cells[d];
    lower[d] = vm->lower[d];
    upper[d] = vm->upper[d];
  }
  for (int d=0; d<vdim; ++d) {
    cells[cdim+d] = s->info.cells[d];
    lower[cdim+d] = s->info.lower[d];
    upper[cdim+d] = s->info.upper[d];
  }
  gkyl_rect_grid_init(&s->grid, pdim, lower, upper, cells);
  init_phase_ranges(cdim, pdim, cells, &s->local, &s->local_ext);

  // allocate distribution function arrays
  s->f = mkarr(app->basis.numBasis, s->local_ext.volume);

  // allocate data for momentum (for use in current accumulation)
  vm_species_moment_init(app, s, &s->m1i, "M1i");

  int ndm = s->info.num_diag_moments;
  // allocate data for diagnostic moments
  s->moms = gkyl_malloc(ndm*sizeof(struct vm_species_moment));
  
  for (int m=0; m<ndm; ++m)
    vm_species_moment_init(app, s, &s->moms[m], s->info.diag_moments[m]);
}

// release resources for species
static void
vm_species_release(struct vm_species *s)
{
  // release various distribution function arrays
  gkyl_array_release(s->f);

  // release moment data
  vm_species_moment_release(&s->m1i);
  for (int i=0; i<s->info.num_diag_moments; ++i)
    vm_species_moment_release(&s->moms[i]);
  gkyl_free(s->moms);
}


gkyl_vlasov_app*
gkyl_vlasov_app_new(struct gkyl_vm vm)
{
  gkyl_vlasov_app *app = gkyl_malloc(sizeof(gkyl_vlasov_app));
  
  int cdim = app->cdim = vm.cdim;
  int vdim = app->vdim = vm.vdim;
  int pdim = cdim+vdim;
  int poly_order = app->poly_order = vm.poly_order;
  int ns = app->num_species = vm.num_species;

  strcpy(app->name, vm.name);

  // basis functions
  gkyl_cart_modal_serendip(&app->basis, pdim, poly_order);
  gkyl_cart_modal_serendip(&app->confBasis, cdim, poly_order);

  // create config grid & ranges
  gkyl_rect_grid_init(&app->grid, cdim, vm.lower, vm.upper, vm.cells);
  init_conf_ranges(cdim, vm.cells, &app->local, &app->local_ext);

  // initialize EM field
  app->field.info = vm.field;
  vm_field_init(&vm, app, &app->field);

  // allocate space to store species objects
  app->species = ns>0 ? gkyl_malloc(ns*sizeof(struct vm_species)) : 0;
  // create species grid & ranges
  for (int i=0; i<ns; ++i) {
    app->species[i].info = vm.species[i];
    vm_species_init(&vm, app, &app->species[i]);
  }
  
  return app;
}

void
gkyl_vlasov_app_init_field(gkyl_vlasov_app* app)
{
  int poly_order = app->poly_order;
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
    poly_order+1, 8, app->field.info.init, app->field.info.ctx);
  
  gkyl_proj_on_basis_advance(proj, 0.0, &app->local, app->field.em);
  gkyl_proj_on_basis_release(proj);
}

void
gkyl_vlasov_app_init_species(gkyl_vlasov_app* app, int sidx)
{
  assert(sidx < app->num_species);
  
  int poly_order = app->poly_order;
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&app->species[sidx].grid, &app->basis,
    poly_order+1, 1, app->species[sidx].info.init, app->species[sidx].info.ctx);
  
  gkyl_proj_on_basis_advance(proj, 0.0, &app->species[sidx].local, app->species[sidx].f);
  gkyl_proj_on_basis_release(proj);  
}

void
gkyl_vlasov_app_init_sim(gkyl_vlasov_app* app)
{
  gkyl_vlasov_app_init_field(app);
  for (int i=0;  i<app->num_species; ++i)
    gkyl_vlasov_app_init_species(app, i);
}

void
gkyl_vlasov_app_calc_mom(gkyl_vlasov_app* app)
{
  for (int i=0; i<app->num_species; ++i) {
    struct vm_species *s = &app->species[i];
    
    for (int m=0; m<app->species[i].info.num_diag_moments; ++m)
      gkyl_mom_calc_advance(s->moms[m].mcalc, &s->local, &app->local,
        s->f, s->moms[m].marr);
  }
}

void
gkyl_vlasov_app_write(gkyl_vlasov_app* app, double tm, int frame)
{
  do { // write EM field
    char fileNm[256];
    sprintf(fileNm, "%s-field_%d.gkyl", app->name, frame);
    
    gkyl_grid_array_write(&app->grid, &app->local, app->field.em, fileNm);
  } while(0);

  // write species distribution function
  for (int i=0; i<app->num_species; ++i) {
    char fileNm[256];
    sprintf(fileNm, "%s-%s_%d.gkyl", app->name, app->species[i].info.name, frame);
    
    gkyl_grid_array_write(&app->species[i].grid, &app->species[i].local,
      app->species[i].f, fileNm);
  }
}

void
gkyl_vlasov_app_mom_write(gkyl_vlasov_app* app, double tm, int frame)
{
  for (int i=0; i<app->num_species; ++i) {

    for (int m=0; m<app->species[i].info.num_diag_moments; ++m) {
      char fileNm[256];
      sprintf(fileNm, "%s-%s-%s_%d.gkyl", app->name, app->species[i].info.name,
        app->species[i].info.diag_moments[m], frame);
      
      gkyl_grid_array_write(&app->grid, &app->local, app->species[i].moms[m].marr, fileNm);
    }
  }
}

void
gkyl_vlasov_app_release(gkyl_vlasov_app* app)
{
  for (int i=0; i<app->num_species; ++i)
    vm_species_release(&app->species[i]);
  gkyl_free(app->species);
  
  vm_field_release(&app->field);
  
  gkyl_free(app);
}