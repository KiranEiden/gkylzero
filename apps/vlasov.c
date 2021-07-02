#include <assert.h>
#include <float.h>
#include <math.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_array_rio.h>
#include <gkyl_dg_maxwell.h>
#include <gkyl_dg_vlasov.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_mom_calc.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_vlasov.h>
#include <gkyl_vlasov_mom.h>

// ranges for use in BCs
struct skin_ghost_ranges {
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];

  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];
};

// data for moments
struct vm_species_moment {
  struct gkyl_mom_type *mtype;
  gkyl_mom_calc *mcalc;
  struct gkyl_array *marr;
  struct gkyl_array *marr_host;  
};

// species data
struct vm_species {
  struct gkyl_vlasov_species info; // data for species
    
  struct gkyl_rect_grid grid;
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  struct skin_ghost_ranges skin_ghost; // conf-space skin/ghost

  struct gkyl_array *f, *f1, *fnew; // arrays for updates
  struct gkyl_array *cflrate; // CFL rate in each cell
  struct gkyl_array *maxs_by_cell; // maximum 'speed' in each direction for global LF flux
  struct gkyl_array *bc_buffer; // buffer for BCs (used for both copy and periodic)

  struct gkyl_array *f_host; // host copy for use IO and initialization

  struct vm_species_moment m1i; // for computing currents
  struct vm_species_moment *moms; // diagnostic moments

  struct gkyl_dg_eqn *eqn; // Vlasov equation
  gkyl_hyper_dg *slvr; // solver 

  double* omegaCfl_ptr;
};

// field data
struct vm_field {
  struct gkyl_vlasov_field info; // data for field
    
  struct gkyl_array *em, *em1, *emnew; // arrays for updates
  struct gkyl_array *qmem; // array for q/m*(E,B)
  struct gkyl_array *cflrate; // CFL rate in each cell
  struct gkyl_array *maxs_by_cell; // maximum 'speed' in each direction for global LF flux
  struct gkyl_array *bc_buffer; // buffer for BCs (used for both copy and periodic)

  struct gkyl_array *em_host;  // host copy for use IO and initialization

  struct gkyl_dg_eqn *eqn; // Maxwell equation
  gkyl_hyper_dg *slvr; // solver

  double* omegaCfl_ptr;
};

// Vlasov object: used as opaque pointer in user code
struct gkyl_vlasov_app {
  char name[128]; // name of app
  int cdim, vdim; // conf, velocity space dimensions
  int poly_order; // polynomial order
  double tcurr; // current time
  double cfl; // CFL number

  bool use_gpu; // should we use GPU (if present)

  int num_periodic_dir; // number of periodic directions
  int periodic_dirs[3]; // list of periodic directions
    
  struct gkyl_rect_grid grid; // config-space grid
  struct gkyl_range local, local_ext; // local, local-ext conf-space ranges
  struct gkyl_basis basis, confBasis; // phase-space, conf-space basis

  struct skin_ghost_ranges skin_ghost; // conf-space skin/ghost

  struct vm_field field; // field data

  // species data
  int num_species;
  struct vm_species *species; // species data

  struct gkyl_vlasov_stat stat; // statistics
};

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(bool on_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (on_gpu)
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

// Compute out = c1*arr1 + c2*arr2
static inline struct gkyl_array*
array_combine(struct gkyl_array *out, double c1, const struct gkyl_array *arr1,
  double c2, const struct gkyl_array *arr2, const struct gkyl_range rng)
{
  return gkyl_array_accumulate_range(gkyl_array_set_range(out, c1, arr1, rng),
    c2, arr2, rng);
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

// Create ghost and skin sub-ranges given a parent range
static void
skin_ghost_ranges_init(struct skin_ghost_ranges *sgr,
  const struct gkyl_range *parent, const int *ghost)
{
  int ndim = parent->ndim;
  
  for (int d=0; d<ndim; ++d) {
    gkyl_skin_ghost_ranges(&sgr->lower_skin[d], &sgr->lower_ghost[d],
      d, GKYL_LOWER_EDGE, parent, ghost);
    gkyl_skin_ghost_ranges(&sgr->upper_skin[d], &sgr->upper_ghost[d],
      d, GKYL_UPPER_EDGE, parent, ghost);
  }
}

// initialize species moment object
static void
vm_species_moment_init(struct gkyl_vlasov_app *app, struct vm_species *s,
  struct vm_species_moment *sm, const char *nm)
{
  assert(is_moment_name_valid(nm));
  
  if (app->use_gpu) {
    struct gkyl_mom_type *mtype_host = gkyl_vlasov_mom_new(&app->confBasis, &app->basis, nm);
    sm->mtype = gkyl_vlasov_mom_cu_dev_new(&app->confBasis, &app->basis, nm);
    sm->mcalc = gkyl_mom_calc_cu_dev_new(&s->grid, sm->mtype);

    sm->marr = mkarr(app->use_gpu, mtype_host->num_mom*app->confBasis.numBasis,
      app->local_ext.volume);

    sm->marr_host = mkarr(false, mtype_host->num_mom*app->confBasis.numBasis,
      app->local_ext.volume);

    gkyl_mom_type_release(mtype_host);
  }
  else {
    sm->mtype = gkyl_vlasov_mom_new(&app->confBasis, &app->basis, nm);
    sm->mcalc = gkyl_mom_calc_new(&s->grid, sm->mtype);

    sm->marr = mkarr(app->use_gpu, sm->mtype->num_mom*app->confBasis.numBasis,
      app->local_ext.volume);

    sm->marr_host = sm->marr;
  }
}

// release memory for moment data object
static void
vm_species_moment_release(const struct gkyl_vlasov_app *app, const struct vm_species_moment *sm)
{
  if (app->use_gpu) {
    // TODO: release dev objects

    gkyl_array_release(sm->marr_host);
  }
  else {
    gkyl_mom_type_release(sm->mtype);
    gkyl_mom_calc_release(sm->mcalc);
  }
  
  gkyl_array_release(sm->marr);
}

// initialize field object
static void
vm_field_init(struct gkyl_vm *vm, struct gkyl_vlasov_app *app, struct vm_field *f)
{
  // allocate EM arrays
  f->em = mkarr(app->use_gpu, 8*app->confBasis.numBasis, app->local_ext.volume);
  f->em1 = mkarr(app->use_gpu, 8*app->confBasis.numBasis, app->local_ext.volume);
  f->emnew = mkarr(app->use_gpu, 8*app->confBasis.numBasis, app->local_ext.volume);
  f->qmem = mkarr(app->use_gpu, 8*app->confBasis.numBasis, app->local_ext.volume);

  f->em_host = f->em;  
  if (app->use_gpu)
    f->em_host = mkarr(false, 8*app->confBasis.numBasis, app->local_ext.volume);

  // allocate buffer for applying BCs (used for both periodic and copy BCs)
  long buff_sz = 0;
  // compute buffer size needed
  for (int d=0; d<app->cdim; ++d) {
    long vol = app->skin_ghost.lower_skin[d].volume;
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  f->bc_buffer = mkarr(app->use_gpu, 8*app->confBasis.numBasis, buff_sz);

  // allocate cflrate (scalar array)
  f->cflrate = mkarr(app->use_gpu, 1, app->local_ext.volume);
  f->maxs_by_cell = mkarr(app->use_gpu, app->cdim, app->local_ext.volume);
  if (app->use_gpu)
    f->omegaCfl_ptr = gkyl_cu_malloc_host(sizeof(double));
  else
    f->omegaCfl_ptr = gkyl_malloc(sizeof(double));

  // equation object
  double c = 1/sqrt(f->info.epsilon0*f->info.mu0);
  double ef = f->info.elcErrorSpeedFactor, mf = f->info.mgnErrorSpeedFactor;

  if (app->use_gpu)
    f->eqn = gkyl_dg_maxwell_cu_dev_new(&app->confBasis, c, ef, mf);
  else
    f->eqn = gkyl_dg_maxwell_new(&app->confBasis, c, ef, mf);

  int up_dirs[] = {0, 1, 2}, zero_flux_flags[] = {0, 0, 0};

  double maxs_init[] = {0., 0., 0.};

  // Maxwell solver
  if (app->use_gpu)
    f->slvr = gkyl_hyper_dg_cu_dev_new(&app->grid, &app->confBasis, f->eqn,
      app->cdim, up_dirs, zero_flux_flags, 1, maxs_init);
  else
    f->slvr = gkyl_hyper_dg_new(&app->grid, &app->confBasis, f->eqn,
      app->cdim, up_dirs, zero_flux_flags, 1, maxs_init);
}

// Compute the RHS for field update, returning maximum stable
// time-step.
static double
vm_field_rhs(gkyl_vlasov_app *app, struct vm_field *field,
  const struct gkyl_array *em, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();
  
  gkyl_array_clear(field->cflrate, 0.0);

  gkyl_array_clear(rhs, 0.0);
  if (app->use_gpu)
    gkyl_hyper_dg_advance_cu(field->slvr, app->local, em, field->cflrate, rhs, field->maxs_by_cell);
  else
    gkyl_hyper_dg_advance(field->slvr, app->local, em, field->cflrate, rhs, field->maxs_by_cell);

  gkyl_array_reduce_range(field->omegaCfl_ptr, field->cflrate, GKYL_MAX, app->local);
  gkyl_array_reduce_range(field->slvr->maxs, field->maxs_by_cell, GKYL_MAX, app->local);
  double omegaCfl = field->omegaCfl_ptr[0];

  app->stat.field_rhs_tm += gkyl_time_diff_now_sec(wst);
  
  return app->cfl/omegaCfl;
}

// Apply periodic BCs on EM fields
static void
vm_field_apply_periodic_bc(gkyl_vlasov_app *app, const struct vm_field *field,
  int dir, struct gkyl_array *f)
{
  gkyl_array_copy_to_buffer(field->bc_buffer->data, f, app->skin_ghost.lower_skin[dir]);
  gkyl_array_copy_from_buffer(f, field->bc_buffer->data, app->skin_ghost.upper_ghost[dir]);

  gkyl_array_copy_to_buffer(field->bc_buffer->data, f, app->skin_ghost.upper_skin[dir]);
  gkyl_array_copy_from_buffer(f, field->bc_buffer->data, app->skin_ghost.lower_ghost[dir]);
}

// Apply copy BCs on EM fields
static void
vm_field_apply_copy_bc(gkyl_vlasov_app *app, const struct vm_field *field,
  int dir, struct gkyl_array *f)
{
  gkyl_array_copy_to_buffer(field->bc_buffer->data, f, app->skin_ghost.lower_skin[dir]);
  gkyl_array_copy_from_buffer(f, field->bc_buffer->data, app->skin_ghost.lower_ghost[dir]);

  gkyl_array_copy_to_buffer(field->bc_buffer->data, f, app->skin_ghost.upper_skin[dir]);
  gkyl_array_copy_from_buffer(f, field->bc_buffer->data, app->skin_ghost.upper_ghost[dir]);
}

// Determine which directions are periodic and which directions are copy,
// and then apply boundary conditions for EM fields
static void
vm_field_apply_bc(gkyl_vlasov_app *app, const struct vm_field *field, struct gkyl_array *f)
{
  int num_periodic_dir = app->num_periodic_dir, cdim = app->cdim, is_copy[3] = {1, 1, 1};
  for (int d=0; d<num_periodic_dir; ++d) {
    vm_field_apply_periodic_bc(app, field, app->periodic_dirs[d], f);
    is_copy[app->periodic_dirs[d]] = 0;
  }
  for (int d=0; d<cdim; ++d)
    if (is_copy[d])
      vm_field_apply_copy_bc(app, field, d, f);
}

// release resources for field
static void
vm_field_release(const gkyl_vlasov_app* app, const struct vm_field *f)
{
  gkyl_array_release(f->em);
  gkyl_array_release(f->em1);
  gkyl_array_release(f->emnew);
  gkyl_array_release(f->qmem);
  gkyl_array_release(f->bc_buffer);
  gkyl_array_release(f->cflrate);
  gkyl_array_release(f->maxs_by_cell);

  if (app->use_gpu) {
    gkyl_array_release(f->em_host);
    gkyl_cu_free_host(f->omegaCfl_ptr);
  }
  else {
    gkyl_dg_eqn_release(f->eqn);
    gkyl_hyper_dg_release(f->slvr);
    gkyl_free(f->omegaCfl_ptr);
  }
}

// initialize species object
static void
vm_species_init(struct gkyl_vm *vm, struct gkyl_vlasov_app *app, struct vm_species *s)
{
  int cdim = app->cdim, vdim = app->vdim;
  int pdim = cdim+vdim;

  int cells[GKYL_MAX_DIM], ghost[GKYL_MAX_DIM];
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];

  for (int d=0; d<cdim; ++d) {
    cells[d] = vm->cells[d];
    lower[d] = vm->lower[d];
    upper[d] = vm->upper[d];
    ghost[d] = 1;
  }
  for (int d=0; d<vdim; ++d) {
    cells[cdim+d] = s->info.cells[d];
    lower[cdim+d] = s->info.lower[d];
    upper[cdim+d] = s->info.upper[d];
    ghost[cdim+d] = 0; // no ghost-cells in velocity space
  }
  gkyl_rect_grid_init(&s->grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&s->grid, ghost, &s->local_ext, &s->local);

  skin_ghost_ranges_init(&s->skin_ghost, &s->local_ext, ghost);
  
  // allocate distribution function arrays
  s->f = mkarr(app->use_gpu, app->basis.numBasis, s->local_ext.volume);
  s->f1 = mkarr(app->use_gpu, app->basis.numBasis, s->local_ext.volume);
  s->fnew = mkarr(app->use_gpu, app->basis.numBasis, s->local_ext.volume);

  s->f_host = s->f;
  if (app->use_gpu)
    s->f_host = mkarr(false, app->basis.numBasis, s->local_ext.volume);

  // allocate cflrate (scalar array) and maxs
  s->cflrate = mkarr(app->use_gpu, 1, s->local_ext.volume);
  s->maxs_by_cell = mkarr(app->use_gpu, pdim, s->local_ext.volume);
  if (app->use_gpu)
    s->omegaCfl_ptr = gkyl_cu_malloc_host(sizeof(double));
  else
    s->omegaCfl_ptr = gkyl_malloc(sizeof(double));

  // allocate buffer for applying periodic BCs
  long buff_sz = 0;
  // compute buffer size needed
  for (int d=0; d<cdim; ++d) {
    long vol = s->skin_ghost.lower_skin[d].volume;
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  s->bc_buffer = mkarr(app->use_gpu, app->basis.numBasis, buff_sz);
  
  // allocate data for momentum (for use in current accumulation)
  vm_species_moment_init(app, s, &s->m1i, "M1i");

  int ndm = s->info.num_diag_moments;
  // allocate data for diagnostic moments
  s->moms = gkyl_malloc(sizeof(struct vm_species_moment[ndm]));
  for (int m=0; m<ndm; ++m)
    vm_species_moment_init(app, s, &s->moms[m], s->info.diag_moments[m]);

  // create equation object
  if (app->use_gpu)
    s->eqn = gkyl_dg_vlasov_cu_dev_new(&app->confBasis, &app->basis, &app->local);
  else
    s->eqn = gkyl_dg_vlasov_new(&app->confBasis, &app->basis, &app->local);

  int up_dirs[GKYL_MAX_DIM], zero_flux_flags[GKYL_MAX_DIM];
  double maxs_init[GKYL_MAX_DIM];
  for (int d=0; d<cdim; ++d) {
    up_dirs[d] = d;
    zero_flux_flags[d] = 0;
    maxs_init[d] = 0.;
  }
  for (int d=cdim; d<pdim; ++d) {
    up_dirs[d] = d;
    zero_flux_flags[d] = 1; // zero-flux BCs in vel-space
    maxs_init[d] = 0.;
  }
  
  // create solver
  if (app->use_gpu)
    s->slvr = gkyl_hyper_dg_cu_dev_new(&s->grid, &app->basis, s->eqn,
      pdim, up_dirs, zero_flux_flags, 1, maxs_init);
  else 
    s->slvr = gkyl_hyper_dg_new(&s->grid, &app->basis, s->eqn,
      pdim, up_dirs, zero_flux_flags, 1, maxs_init);
}

// Compute the RHS for species update, returning maximum stable
// time-step.
static double
vm_species_rhs(gkyl_vlasov_app *app, struct vm_species *species,
  const struct gkyl_array *fin, const struct gkyl_array *qmem, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();
  
  gkyl_array_clear(species->cflrate, 0.0);
  gkyl_vlasov_set_qmem(species->eqn, qmem); // must set EM fields to use
  
  gkyl_array_clear(rhs, 0.0);
  if (app->use_gpu)
    gkyl_hyper_dg_advance_cu(species->slvr, species->local, fin, species->cflrate, rhs, species->maxs_by_cell);
  else
    gkyl_hyper_dg_advance(species->slvr, species->local, fin, species->cflrate, rhs, species->maxs_by_cell);

  gkyl_array_reduce_range(species->omegaCfl_ptr, species->cflrate, GKYL_MAX, species->local);
  gkyl_array_reduce_range(species->slvr->maxs, species->maxs_by_cell, GKYL_MAX, species->local);
  double omegaCfl = species->omegaCfl_ptr[0];

  app->stat.species_rhs_tm += gkyl_time_diff_now_sec(wst);
  
  return app->cfl/omegaCfl;
}

// Apply periodic BCs on distribution function
static void
vm_species_apply_periodic_bc(gkyl_vlasov_app *app, const struct vm_species *species,
  int dir, struct gkyl_array *f)
{
  gkyl_array_copy_to_buffer(species->bc_buffer->data, f, species->skin_ghost.lower_skin[dir]);
  gkyl_array_copy_from_buffer(f, species->bc_buffer->data, species->skin_ghost.upper_ghost[dir]);

  gkyl_array_copy_to_buffer(species->bc_buffer->data, f, species->skin_ghost.upper_skin[dir]);
  gkyl_array_copy_from_buffer(f, species->bc_buffer->data, species->skin_ghost.lower_ghost[dir]);
}

// Apply copy BCs on distribution function
static void
vm_species_apply_copy_bc(gkyl_vlasov_app *app, const struct vm_species *species,
  int dir, struct gkyl_array *f)
{
  gkyl_array_copy_to_buffer(species->bc_buffer->data, f, species->skin_ghost.lower_skin[dir]);
  gkyl_array_copy_from_buffer(f, species->bc_buffer->data, species->skin_ghost.lower_ghost[dir]);

  gkyl_array_copy_to_buffer(species->bc_buffer->data, f, species->skin_ghost.upper_skin[dir]);
  gkyl_array_copy_from_buffer(f, species->bc_buffer->data, species->skin_ghost.upper_ghost[dir]);
}

// Determine which directions are periodic and which directions are copy,
// and then apply boundary conditions for distribution function
static void
vm_species_apply_bc(gkyl_vlasov_app *app, const struct vm_species *species, struct gkyl_array *f)
{
  int num_periodic_dir = app->num_periodic_dir, cdim = app->cdim, is_copy[3] = {1, 1, 1};
  for (int d=0; d<num_periodic_dir; ++d) {
    vm_species_apply_periodic_bc(app, species, app->periodic_dirs[d], f);
    is_copy[app->periodic_dirs[d]] = 0;
  }
  for (int d=0; d<cdim; ++d)
    if (is_copy[d])
      vm_species_apply_copy_bc(app, species, d, f);
}

// release resources for species
static void
vm_species_release(const gkyl_vlasov_app* app, const struct vm_species *s)
{
  // release various arrays
  gkyl_array_release(s->f);
  gkyl_array_release(s->f1);
  gkyl_array_release(s->fnew);
  gkyl_array_release(s->cflrate);
  gkyl_array_release(s->maxs_by_cell);
  gkyl_array_release(s->bc_buffer);

  if (app->use_gpu) {
    gkyl_array_release(s->f_host);
  }

  // release moment data
  vm_species_moment_release(app, &s->m1i);
  for (int i=0; i<s->info.num_diag_moments; ++i)
    vm_species_moment_release(app, &s->moms[i]);
  gkyl_free(s->moms);

  if (app->use_gpu) {
    gkyl_cu_free_host(s->omegaCfl_ptr);
    // TODO: NOT SURE HOW TO RELEASE ON DEVICE YET
  }
  else {
    gkyl_free(s->omegaCfl_ptr);
    gkyl_dg_eqn_release(s->eqn);
    gkyl_hyper_dg_release(s->slvr);
  }
}

gkyl_vlasov_app*
gkyl_vlasov_app_new(struct gkyl_vm vm)
{
  assert(vm.num_species <= GKYL_MAX_SPECIES);
  
  gkyl_vlasov_app *app = gkyl_malloc(sizeof(gkyl_vlasov_app));
  
  int cdim = app->cdim = vm.cdim;
  int vdim = app->vdim = vm.vdim;
  int pdim = cdim+vdim;
  int poly_order = app->poly_order = vm.poly_order;
  int ns = app->num_species = vm.num_species;

  double cfl_frac = vm.cfl_frac == 0 ? 1.0 : vm.cfl_frac;
  app->cfl = cfl_frac/(2*poly_order+1);

#ifdef GKYL_HAVE_CUDA
  app->use_gpu = vm.use_gpu;
#else
  app->use_gpu = false; // can't use GPUs if we don't have them!
#endif
  
  app->num_periodic_dir = vm.num_periodic_dir;
  for (int d=0; d<cdim; ++d)
    app->periodic_dirs[d] = vm.periodic_dirs[d];

  strcpy(app->name, vm.name);
  app->tcurr = 0.0; // reset on init

  // basis functions
  gkyl_cart_modal_serendip(&app->basis, pdim, poly_order);
  gkyl_cart_modal_serendip(&app->confBasis, cdim, poly_order);

  gkyl_rect_grid_init(&app->grid, cdim, vm.lower, vm.upper, vm.cells);

  int ghost[] = { 1, 1, 1 };  
  gkyl_create_grid_ranges(&app->grid, ghost, &app->local_ext, &app->local);
  skin_ghost_ranges_init(&app->skin_ghost, &app->local_ext, ghost);

  // initialize EM field
  app->field.info = vm.field;
  vm_field_init(&vm, app, &app->field);

  // allocate space to store species objects
  app->species = ns>0 ? gkyl_malloc(sizeof(struct vm_species[ns])) : 0;
  // create species grid & ranges
  for (int i=0; i<ns; ++i) {
    app->species[i].info = vm.species[i];
    vm_species_init(&vm, app, &app->species[i]);
  }

  // initialize stat object
  app->stat = (struct gkyl_vlasov_stat) {
    .stage_2_dt_diff = { DBL_MAX, 0.0 },
    .stage_3_dt_diff = { DBL_MAX, 0.0 },
  };
  
  return app;
}

void
gkyl_vlasov_app_apply_ic(gkyl_vlasov_app* app, double t0)
{
  app->tcurr = t0;
  gkyl_vlasov_app_apply_ic_field(app, t0);
  for (int i=0;  i<app->num_species; ++i)
    gkyl_vlasov_app_apply_ic_species(app, i, t0);
}

void
gkyl_vlasov_app_apply_ic_field(gkyl_vlasov_app* app, double t0)
{
  app->tcurr = t0;
  int poly_order = app->poly_order;
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
    poly_order+1, 8, app->field.info.init, app->field.info.ctx);

  gkyl_proj_on_basis_advance(proj, t0, &app->local, app->field.em_host);
  gkyl_proj_on_basis_release(proj);
 
  if (app->use_gpu)
    gkyl_array_copy(app->field.em, app->field.em_host);
  
  vm_field_apply_bc(app, &app->field, app->field.em);  
}

void
gkyl_vlasov_app_apply_ic_species(gkyl_vlasov_app* app, int sidx, double t0)
{
  assert(sidx < app->num_species);

  app->tcurr = t0;
  int poly_order = app->poly_order;
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&app->species[sidx].grid, &app->basis,
    poly_order+1, 1, app->species[sidx].info.init, app->species[sidx].info.ctx);

  gkyl_proj_on_basis_advance(proj, t0, &app->species[sidx].local, app->species[sidx].f_host);
  gkyl_proj_on_basis_release(proj);

  if (app->use_gpu)
    gkyl_array_copy(app->species[sidx].f, app->species[sidx].f_host);

  vm_species_apply_bc(app, &app->species[sidx], app->species[sidx].f);
}

void
gkyl_vlasov_app_calc_mom(gkyl_vlasov_app* app)
{
  for (int i=0; i<app->num_species; ++i) {
    struct vm_species *s = &app->species[i];
    
    for (int m=0; m<app->species[i].info.num_diag_moments; ++m)
      if (app->use_gpu)
        gkyl_mom_calc_advance_cu(s->moms[m].mcalc, s->local, app->local,
          s->f, s->moms[m].marr);
      else
        gkyl_mom_calc_advance(s->moms[m].mcalc, s->local, app->local,
          s->f, s->moms[m].marr);
  }
}

void
gkyl_vlasov_app_write(gkyl_vlasov_app* app, double tm, int frame)
{
  gkyl_vlasov_app_write_field(app, tm, frame);
  for (int i=0; i<app->num_species; ++i)
    gkyl_vlasov_app_write_species(app, i, tm, frame);
}
 
void
gkyl_vlasov_app_write_field(gkyl_vlasov_app* app, double tm, int frame)
{
  const char *fmt = "%s-field_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, frame);

  if (app->use_gpu) {
    // copy data from device to host before writing it out    
    gkyl_array_copy(app->field.em_host, app->field.em);
    gkyl_grid_array_write(&app->grid, &app->local, app->field.em_host, fileNm);
  }
  else {
    gkyl_grid_array_write(&app->grid, &app->local, app->field.em, fileNm);
  }
}

void
gkyl_vlasov_app_write_species(gkyl_vlasov_app* app, int sidx, double tm, int frame)
{
  const char *fmt = "%s-%s_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, app->species[sidx].info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow  
  snprintf(fileNm, sizeof fileNm, fmt, app->name, app->species[sidx].info.name, frame);

  if (app->use_gpu) {
    // copy data from device to host before writing it out
    gkyl_array_copy(app->species[sidx].f_host, app->species[sidx].f);
    gkyl_grid_array_write(&app->species[sidx].grid, &app->species[sidx].local,
      app->species[sidx].f_host, fileNm);
  }
  else {
    gkyl_grid_array_write(&app->species[sidx].grid, &app->species[sidx].local,
      app->species[sidx].f, fileNm);
  }
}

void
gkyl_vlasov_app_write_mom(gkyl_vlasov_app* app, double tm, int frame)
{
  for (int i=0; i<app->num_species; ++i) {
    for (int m=0; m<app->species[i].info.num_diag_moments; ++m) {
      
      const char *fmt = "%s-%s-%s_%d.gkyl";
      int sz = gkyl_calc_strlen(fmt, app->name, app->species[i].info.name,
        app->species[i].info.diag_moments[m], frame);
      char fileNm[sz+1]; // ensures no buffer overflow  
      snprintf(fileNm, sizeof fileNm, fmt, app->name, app->species[i].info.name,
        app->species[i].info.diag_moments[m], frame);

      if (app->use_gpu)
        gkyl_array_copy(app->species[i].moms[m].marr_host, app->species[i].moms[m].marr);
      
      gkyl_grid_array_write(&app->grid, &app->local, app->species[i].moms[m].marr_host, fileNm);
    }
  }
}

// Take a forward Euler step with the suggested time-step dt. This may
// not be the actual time-step taken. However, the function will never
// take a time-step larger than dt even if it is allowed by
// stability. The actual time-step and dt_suggested are returned in
// the status object.
static void
forward_euler(gkyl_vlasov_app* app, double tcurr, double dt,
  const struct gkyl_array *fin[], const struct gkyl_array *emin,
  struct gkyl_array *fout[], struct gkyl_array *emout, struct gkyl_update_status *st)
{
  app->stat.nfeuler += 1;
  
  double dtmin = dt;

  // compute RHS of Vlasov equations
  for (int i=0; i<app->num_species; ++i) {
    double qbym = app->species[i].info.charge/app->species[i].info.mass;
    gkyl_array_set(app->field.qmem, qbym, emin);
    
    double dt1 = vm_species_rhs(app, &app->species[i], fin[i], app->field.qmem, fout[i]);
    dtmin = fmin(dtmin, dt1);
  }

  // compute RHS of Maxwell equations
  double dt1 = vm_field_rhs(app, &app->field, emin, emout);
  dtmin = fmin(dtmin, dt1);

  double dt_max_rel_diff = 0.01;
  // check if dtmin is slightly smaller than dt. Use dt if it is
  // (avoids retaking steps if dt changes are very small).
  double dt_rel_diff = (dt-dtmin)/dt;
  if (dt_rel_diff > 0 && dt_rel_diff < dt_max_rel_diff)
    dtmin = dt;

  // don't take a time-step larger that input dt
  double dta = st->dt_actual = dt < dtmin ? dt : dtmin;
  st->dt_suggested = dtmin;

  // complete update of distribution function
  for (int i=0; i<app->num_species; ++i) {
    gkyl_array_accumulate_range(gkyl_array_scale_range(fout[i], dta, app->species[i].local),
      1.0, fin[i], app->species[i].local);
    vm_species_apply_bc(app, &app->species[i], fout[i]);
  }

  struct timespec wst = gkyl_wall_clock();
  // accumulate current contribution to electric field terms
  for (int i=0; i<app->num_species; ++i) {
    struct vm_species *s = &app->species[i];    
    if (app->use_gpu)
      gkyl_mom_calc_advance_cu(s->m1i.mcalc, s->local, app->local,
        fin[i], s->m1i.marr);
    else
      gkyl_mom_calc_advance(s->m1i.mcalc, s->local, app->local,
        fin[i], s->m1i.marr);
    
    double qbyeps = s->info.charge/app->field.info.epsilon0;
    gkyl_array_accumulate_range(emout, -qbyeps, s->m1i.marr, app->local);
  }
  app->stat.current_tm += gkyl_time_diff_now_sec(wst);
  
  // complete update of field
  gkyl_array_accumulate_range(gkyl_array_scale_range(emout, dta, app->local),
    1.0, emin, app->local);
  vm_field_apply_bc(app, &app->field, emout);
}

// Take time-step using the RK3 method. Also sets the status object
// which has the actual and suggested dts used. These can be different
// from the actual time-step.
static struct gkyl_update_status
rk3(gkyl_vlasov_app* app, double dt0)
{
  const struct gkyl_array *fin[app->num_species];
  struct gkyl_array *fout[app->num_species];
  struct gkyl_update_status st = { .success = true };

  // time-stepper state
  enum { RK_STAGE_1, RK_STAGE_2, RK_STAGE_3, RK_COMPLETE } state = RK_STAGE_1;

  double tcurr = app->tcurr, dt = dt0;
  while (state != RK_COMPLETE) {
    switch (state) {
      case RK_STAGE_1:
        for (int i=0; i<app->num_species; ++i) {
          fin[i] = app->species[i].f;
          fout[i] = app->species[i].f1;
        }
        forward_euler(app, tcurr, dt, fin, app->field.em, fout, app->field.em1, &st);
        dt = st.dt_actual;
        state = RK_STAGE_2;
        break;

      case RK_STAGE_2:
        for (int i=0; i<app->num_species; ++i) {
          fin[i] = app->species[i].f1;
          fout[i] = app->species[i].fnew;
        }
        forward_euler(app, tcurr+dt, dt, fin, app->field.em1, fout, app->field.emnew, &st);
        if (st.dt_actual < dt) {

          // collect stats
          double dt_rel_diff = (dt-st.dt_actual)/st.dt_actual;
          app->stat.stage_2_dt_diff[0] = fmin(app->stat.stage_2_dt_diff[0],
            dt_rel_diff);
          app->stat.stage_2_dt_diff[1] = fmax(app->stat.stage_2_dt_diff[1],
            dt_rel_diff);
          app->stat.nstage_2_fail += 1;
            
          dt = st.dt_actual;
          state = RK_STAGE_1; // restart from stage 1

        }
        else {
          for (int i=0; i<app->num_species; ++i)
            array_combine(app->species[i].f1,
              3.0/4.0, app->species[i].f, 1.0/4.0, app->species[i].fnew, app->species[i].local_ext);
          array_combine(app->field.em1,
            3.0/4.0, app->field.em, 1.0/4.0, app->field.emnew, app->local_ext);

          state = RK_STAGE_3;
        }
        break;

      case RK_STAGE_3:
        for (int i=0; i<app->num_species; ++i) {
          fin[i] = app->species[i].f1;
          fout[i] = app->species[i].fnew;
        }
        forward_euler(app, tcurr+dt/2, dt, fin, app->field.em1, fout, app->field.emnew, &st);
        if (st.dt_actual < dt) {

          // collect stats
          double dt_rel_diff = (dt-st.dt_actual)/st.dt_actual;
          app->stat.stage_3_dt_diff[0] = fmin(app->stat.stage_3_dt_diff[0],
            dt_rel_diff);
          app->stat.stage_3_dt_diff[1] = fmax(app->stat.stage_3_dt_diff[1],
            dt_rel_diff);
          app->stat.nstage_3_fail += 1;
            
          dt = st.dt_actual;
          state = RK_STAGE_1; // restart from stage 1
            
          app->stat.nstage_2_fail += 1;
        }
        else {
          for (int i=0; i<app->num_species; ++i) {
            array_combine(app->species[i].f1,
              1.0/3.0, app->species[i].f, 2.0/3.0, app->species[i].fnew, app->species[i].local_ext);
            gkyl_array_copy_range(app->species[i].f, app->species[i].f1, app->species[i].local_ext);
          }
          array_combine(app->field.em1,
            1.0/3.0, app->field.em, 2.0/3.0, app->field.emnew, app->local_ext);
          gkyl_array_copy_range(app->field.em, app->field.em1, app->local_ext);

          state = RK_COMPLETE;
        }
        break;

      case RK_COMPLETE: // can't happen: suppresses warning
        break;
    }
  }
  
  return st;
}

struct gkyl_update_status
gkyl_vlasov_update(gkyl_vlasov_app* app, double dt)
{
  app->stat.nup += 1;
  struct timespec wst = gkyl_wall_clock();

  struct gkyl_update_status status = rk3(app, dt);
  app->tcurr += status.dt_actual;
  
  app->stat.total_tm += gkyl_time_diff_now_sec(wst);
  
  return status;
}

struct gkyl_vlasov_stat
gkyl_vlasov_app_stat(gkyl_vlasov_app* app)
{
  return app->stat;
}

void
gkyl_vlasov_app_species_ktm_rhs(gkyl_vlasov_app* app, int update_vol_term)
{
  for (int i=0; i<app->num_species; ++i) {
    
    struct vm_species *species = &app->species[i];
    
    const struct gkyl_array *qmem = app->field.qmem;
    gkyl_vlasov_set_qmem(species->eqn, qmem);

    const struct gkyl_array *fin = species->f;
    struct gkyl_array *rhs = species->f1;

    if (app->use_gpu)
      gkyl_hyper_dg_set_update_vol_cu(species->slvr, update_vol_term);
    else
      gkyl_hyper_dg_set_update_vol(species->slvr, update_vol_term);
    gkyl_array_clear_range(rhs, 0.0, species->local);
    if (app->use_gpu)
      gkyl_hyper_dg_advance_cu(species->slvr, species->local, fin,
        species->cflrate, rhs, species->maxs_by_cell);
    else
      gkyl_hyper_dg_advance(species->slvr, species->local, fin,
        species->cflrate, rhs, species->maxs_by_cell);

    // reduction to get maxs (maximum over cells of maxs_by_cell)
    gkyl_array_reduce_range(species->slvr->maxs, species->maxs_by_cell, GKYL_MAX, species->local);
  }
}

void
gkyl_vlasov_app_stat_write(const gkyl_vlasov_app* app)
{
  const char *fmt = "%s-%s";
  int sz = gkyl_calc_strlen(fmt, app->name, "stat.json");
  char fileNm[sz+1]; // ensures no buffer overflow  
  snprintf(fileNm, sizeof fileNm, fmt, app->name, "stat.json");
  
  char buff[70];
  time_t t = time(NULL);
  struct tm curr_tm = *localtime(&t);

  // append to existing file so we have a history of different runs
  FILE *fp = 0;
  with_file (fp, fileNm, "a") {
    fprintf(fp, "{\n");

    if (strftime(buff, sizeof buff, "%c", &curr_tm))
      fprintf(fp, " \"date\" : \"%s\"\n", buff);

    fprintf(fp, " \"nup\" : \"%ld\"\n", app->stat.nup);
    fprintf(fp, " \"nfeuler\" : \"%ld\"\n", app->stat.nfeuler);
    fprintf(fp, " \"nstage_2_fail\" : \"%ld\"\n", app->stat.nstage_2_fail);
    fprintf(fp, " \"nstage_3_fail\" : \"%ld\"\n", app->stat.nstage_3_fail);

    fprintf(fp, " \"stage_2_dt_diff\" : [ \"%lg\", \"%lg\" ]\n",
      app->stat.stage_2_dt_diff[0], app->stat.stage_2_dt_diff[1]);
    fprintf(fp, " \"stage_3_dt_diff\" : [ \"%lg\", \"%lg\" ]\n",
      app->stat.stage_3_dt_diff[0], app->stat.stage_3_dt_diff[1]);
    
    fprintf(fp, " \"total_tm\" : \"%lg\"\n", app->stat.total_tm);
    fprintf(fp, " \"init_species_tm\" : \"%lg\"\n", app->stat.init_species_tm);
    fprintf(fp, " \"init_field_tm\" : \"%lg\"\n", app->stat.init_field_tm);
    
    fprintf(fp, " \"species_rhs_tm\" : \"%lg\"\n", app->stat.species_rhs_tm);
    fprintf(fp, " \"field_rhs_tm\" : \"%lg\"\n", app->stat.field_rhs_tm);
    fprintf(fp, " \"current_tm\" : \"%lg\"\n", app->stat.current_tm);

    fprintf(fp, " \"nmom\" : \"%ld\"\n", app->stat.nmom);
    fprintf(fp, " \"mom_tm\" : \"%lg\"\n", app->stat.mom_tm);
  
    fprintf(fp, "}\n");
  }
}

void
gkyl_vlasov_app_release(gkyl_vlasov_app* app)
{
  for (int i=0; i<app->num_species; ++i)
    vm_species_release(app, &app->species[i]);
  gkyl_free(app->species);
  
  vm_field_release(app, &app->field);

  gkyl_free(app);
}
