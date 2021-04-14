#include <assert.h>
#include <float.h>
#include <math.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_fv_proj.h>
#include <gkyl_moment.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_wave_prop.h>
#include <gkyl_wv_maxwell.h>

// species data
struct moment_species {
    int ndim;
    char name[128]; // species name
    double charge, mass;

    int evolve; // evolve species? 1-yes, 0-no

    void *ctx; // context for initial condition init function
    // pointer to initialization function
    void (*init)(double t, const double *xn, double *fout, void *ctx);
    
    struct gkyl_array *f, *fx[3]; // arrays for updates
    struct gkyl_array *bc_buffer; // buffer for periodic BCs

    int num_equations; // number of equations in species
    gkyl_wave_prop *slvr[3]; // solver in each direction
};

struct moment_field {
    int ndim;
    double epsilon0, mu0;

    int evolve; // evolve species? 1-yes, 0-no
    
    void *ctx; // context for initial condition init function
    // pointer to initialization function
    void (*init)(double t, const double *xn, double *fout, void *ctx);    
    
    struct gkyl_array *f, *fx[3]; // arrays for updates
    struct gkyl_array *bc_buffer; // buffer for periodic BCs

    gkyl_wave_prop *slvr[3]; // solver in each direction
};

// Moment object: used as opaque pointer in user code
struct gkyl_moment_app {
    char name[128]; // name of app
    int ndim; // space dimensions
    double tcurr; // current time
    double cfl; // CFL number

    int num_periodic_dir; // number of periodic directions
    int periodic_dirs[3]; // list of periodic directions
    
    struct gkyl_rect_grid grid; // grid
    struct gkyl_range local, local_ext; // local, local-ext ranges

    int has_field; // flag to indicate if we have a field
    struct moment_field field; // field data
    
    // species data
    int num_species;
    struct moment_species *species; // species data

    struct gkyl_moment_stat stat; // statistics
};

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

// initialize species
static void
moment_species_init(const struct gkyl_moment *mom, const struct gkyl_moment_species *mom_sp,
  struct gkyl_moment_app *app, struct moment_species *sp)
{
  sp->ndim = mom->ndim;
  strcpy(sp->name, mom_sp->name);
  sp->charge = mom_sp->charge;
  sp->mass = mom_sp->mass;

  sp->ctx = mom_sp->ctx;
  sp->init = mom_sp->init;

  sp->num_equations = mom_sp->equation->num_equations;

  // choose default limiter
  enum gkyl_wave_limiter limiter =
    mom_sp->limiter == 0 ? GKYL_MONOTONIZED_CENTERED : mom_sp->limiter;

  int ndim = mom->ndim;
  // create updaters for each directional update
  for (int d=0; d<ndim; ++d)
    sp->slvr[d] = gkyl_wave_prop_new( (struct gkyl_wave_prop_inp) {
        .grid = &app->grid,
        .equation = mom_sp->equation,
        .limiter = limiter,
        .num_up_dirs = 1,
        .update_dirs = { d },
        .cfl = app->cfl
      }
    );

  int meqn = sp->num_equations;
  // allocate arrays
  sp->f = mkarr(meqn, app->local_ext.volume);
  for (int d=0; d<ndim; ++d)
    sp->fx[d] = mkarr(meqn, app->local_ext.volume);
}

// free species
static void
moment_species_release(const struct moment_species *sp)
{
  for (int d=0; d<sp->ndim; ++d)
    gkyl_wave_prop_release(sp->slvr[d]);

  gkyl_array_release(sp->f);
  for (int d=0; d<sp->ndim; ++d)
    gkyl_array_release(sp->fx[d]);
}

// initialize field
static void
moment_field_init(const struct gkyl_moment *mom, const struct gkyl_moment_field *mom_fld,
  struct gkyl_moment_app *app, struct moment_field *fld)
{
  fld->ndim = mom->ndim;
  double epsilon0 = fld->epsilon0 = mom_fld->epsilon0;
  double mu0 = fld->mu0 = mom_fld->mu0;

  fld->ctx = mom_fld->ctx;
  fld->init = mom_fld->init;

  // choose default limiter
  enum gkyl_wave_limiter limiter =
    mom_fld->limiter == 0 ? GKYL_MONOTONIZED_CENTERED : mom_fld->limiter;

  double c = 1/sqrt(epsilon0*mu0);
  struct gkyl_wv_eqn *maxwell = gkyl_wv_maxwell_new(c,
    mom_fld->elcErrorSpeedFactor, mom_fld->mgnErrorSpeedFactor);

  int ndim = mom->ndim;
  // create updaters for each directional update
  for (int d=0; d<ndim; ++d)
    fld->slvr[d] = gkyl_wave_prop_new( (struct gkyl_wave_prop_inp) {
        .grid = &app->grid,
        .equation = maxwell,
        .limiter = limiter,
        .num_up_dirs = 1,
        .update_dirs = { d },
        .cfl = app->cfl
      }
    );

  // allocate arrays
  fld->f = mkarr(8, app->local_ext.volume);
  for (int d=0; d<ndim; ++d)
    fld->fx[d] = mkarr(8, app->local_ext.volume);

  gkyl_wv_eqn_release(maxwell);
}

// free field
static void
moment_field_release(const struct moment_field *fld)
{
  for (int d=0; d<fld->ndim; ++d)
    gkyl_wave_prop_release(fld->slvr[d]);

  gkyl_array_release(fld->f);
  for (int d=0; d<fld->ndim; ++d)
    gkyl_array_release(fld->fx[d]);
}

gkyl_moment_app*
gkyl_moment_app_new(struct gkyl_moment mom)
{
  struct gkyl_moment_app *app = gkyl_malloc(sizeof(gkyl_moment_app));

  int ndim = app->ndim = mom.ndim;

  // create grid and ranges
  int ghost[3] = { 2, 2, 2 };
  gkyl_rect_grid_init(&app->grid, ndim, mom.lower, mom.upper, mom.cells);
  gkyl_create_grid_ranges(&app->grid, ghost, &app->local_ext, &app->local);

  double cfl_frac = mom.cfl_frac == 0 ? 0.95 : mom.cfl_frac;
  app->cfl = 1.0*cfl_frac;

  app->num_periodic_dir = mom.num_periodic_dir;
  for (int d=0; d<ndim; ++d)
    app->periodic_dirs[d] = mom.periodic_dirs[d];

  app->has_field = 0;
  // initialize field if we have one
  if (app->field.init) {
    app->has_field = 1;
    moment_field_init(&mom, &mom.field, app, &app->field);
  }

  int ns = app->num_species = mom.num_species;
  // allocate space to store species objects
  app->species = ns>0 ? gkyl_malloc(sizeof(struct moment_species[ns])) : 0;
  // create species grid & ranges
  for (int i=0; i<ns; ++i)
    moment_species_init(&mom, &mom.species[i], app, &app->species[i]);

  strcpy(app->name, mom.name);
  app->tcurr = 0.0; // reset on init

  return app;
}

void
gkyl_moment_app_apply_ic(gkyl_moment_app* app, double t0)
{
  app->tcurr = t0;
  gkyl_moment_app_apply_ic_field(app, t0);
  for (int i=0;  i<app->num_species; ++i)
    gkyl_moment_app_apply_ic_species(app, i, t0);
}

void
gkyl_moment_app_apply_ic_field(gkyl_moment_app* app, double t0)
{
  if (app->has_field != 1) return;

  app->tcurr = t0;
  gkyl_fv_proj *proj = gkyl_fv_proj_new(&app->grid, 2, 8, app->field.init, app->field.ctx);
  
  gkyl_fv_proj_advance(proj, t0, &app->local, app->field.f);
  gkyl_fv_proj_release(proj);
  // APPLY BCS!!!
}

void
gkyl_moment_app_apply_ic_species(gkyl_moment_app* app, int sidx, double t0)
{
  assert(sidx < app->num_species);

  app->tcurr = t0;
  gkyl_fv_proj *proj = gkyl_fv_proj_new(&app->grid, 2, app->species[sidx].num_equations,
    app->species[sidx].init, app->species[sidx].ctx);
  
  gkyl_fv_proj_advance(proj, t0, &app->local, app->species[sidx].f);
  gkyl_fv_proj_release(proj);
  // APPLY BCS!!!
}

void
gkyl_moment_app_write(gkyl_moment_app* app, double tm, int frame)
{
  gkyl_moment_app_write_field(app, tm, frame);
  for (int i=0; i<app->num_species; ++i)
    gkyl_moment_app_write_species(app, i, tm, frame);
}

void
gkyl_moment_app_write_field(gkyl_moment_app* app, double tm, int frame)
{
  if (app->has_field != 1) return;

  const char *fmt = "%s-%s_%d.gkyl";
  int sz = snprintf(NULL, 0, fmt, app->name, "field", frame);
  char fileNm[sz+1]; // ensures no buffer overflow  
  snprintf(fileNm, sizeof fileNm, fmt, app->name, "field", frame);
  
  gkyl_grid_array_write(&app->grid, &app->local, app->field.f, fileNm);
}

void
gkyl_moment_app_write_species(gkyl_moment_app* app, int sidx, double tm, int frame)
{
  const char *fmt = "%s-%s_%d.gkyl";
  int sz = snprintf(NULL, 0, fmt, app->name, app->species[sidx].name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow  
  snprintf(fileNm, sizeof fileNm, fmt, app->name, app->species[sidx].name, frame);
  
  gkyl_grid_array_write(&app->grid, &app->local, app->species[sidx].f, fileNm);
}

struct gkyl_update_status
gkyl_moment_update(gkyl_moment_app* app, double dt)
{
  return (struct gkyl_update_status) { };
}

void
gkyl_moment_app_release(gkyl_moment_app* app)
{
  for (int i=0; i<app->num_species; ++i)
    moment_species_release(&app->species[i]);
  free(app->species);

  if (app->has_field)
    moment_field_release(&app->field);
  
  gkyl_free(app);
}