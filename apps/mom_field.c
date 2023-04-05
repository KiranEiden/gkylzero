#include <gkyl_moment_priv.h>

// initialize field
void
moment_field_init(const struct gkyl_moment *mom, const struct gkyl_moment_field *mom_fld,
  struct gkyl_moment_app *app, struct moment_field *fld)
{
  fld->ndim = mom->ndim;
  double epsilon0 = fld->epsilon0 = mom_fld->epsilon0;
  double mu0 = fld->mu0 = mom_fld->mu0;

  fld->ctx = mom_fld->ctx;
  fld->init = mom_fld->init;

  fld->scheme_type = mom->scheme_type;
  
  // choose default limiter
  enum gkyl_wave_limiter limiter =
    mom_fld->limiter == 0 ? GKYL_MONOTONIZED_CENTERED : mom_fld->limiter;

  double c = 1/sqrt(epsilon0*mu0);
  struct gkyl_wv_eqn *maxwell = gkyl_wv_maxwell_new(c,
    mom_fld->elc_error_speed_fact, mom_fld->mag_error_speed_fact);

  fld->maxwell = gkyl_wv_eqn_acquire(maxwell);
  
  int ndim = mom->ndim;

  if (fld->scheme_type == GKYL_MOMENT_WAVE_PROP) {
    // create updaters for each directional update
    for (int d=0; d<ndim; ++d)
      fld->slvr[d] = gkyl_wave_prop_new( &(struct gkyl_wave_prop_inp) {
          .grid = &app->grid,
          .equation = maxwell,
          .limiter = limiter,
          .num_up_dirs = app->is_dir_skipped[d] ? 0 : 1,
          .update_dirs = { d },
          .check_inv_domain = false,
          .cfl = app->cfl,
          .geom = app->geom,
          .comm = app->comm
        }
      );

    // allocate arrays
    fld->fdup = mkarr(false, 8, app->local_ext.volume);
    for (int d=0; d<ndim+1; ++d)
      fld->f[d] = mkarr(false, 8, app->local_ext.volume);

    // set current solution so ICs and IO work properly
    fld->fcurr = fld->f[0];
  }
  else if (fld->scheme_type == GKYL_MOMENT_MP || fld->scheme_type == GKYL_MOMENT_KEP) {
    // NOTE: there is no KEP scheme for Maxwell, and we simply use MP scheme instead
    
    // determine directions to update
    int num_up_dirs = 0, update_dirs[GKYL_MAX_CDIM] = { 0 };
    for (int d=0; d<ndim; ++d)
      if (!app->is_dir_skipped[d]) {
        update_dirs[num_up_dirs] = d;
        num_up_dirs += 1;
      }

    // choose U3 as field reconstruction if using KEP
    enum gkyl_mp_recon mp_recon =
      fld->scheme_type == GKYL_MOMENT_KEP ? GKYL_MP_U3 : app->mp_recon;
    
    // single MP updater updates all directions
    fld->mp_slvr = gkyl_mp_scheme_new( &(struct gkyl_mp_scheme_inp) {
        .grid = &app->grid,
        .equation = maxwell,
        .mp_recon = mp_recon,
        .skip_mp_limiter = mom->skip_mp_limiter,
        .num_up_dirs = num_up_dirs,
        .update_dirs = { update_dirs[0], update_dirs[1], update_dirs[2] } ,
        .cfl = app->cfl,
        .geom = app->geom,
      }
    );

    // allocate arrays
    fld->f0 = mkarr(false, 8, app->local_ext.volume);
    fld->f1 = mkarr(false, 8, app->local_ext.volume);
    fld->fnew = mkarr(false, 8, app->local_ext.volume);
    fld->cflrate = mkarr(false, 1, app->local_ext.volume);

    // set current solution so ICs and IO work properly
    fld->fcurr = fld->f0;
  }

  // determine which directions are not periodic
  int num_periodic_dir = app->num_periodic_dir, is_np[3] = {1, 1, 1};
  for (int d=0; d<num_periodic_dir; ++d)
    is_np[app->periodic_dirs[d]] = 0;

  for (int i=0; i<3; ++i) {
    fld->lower_bc[i] = 0;
    fld->upper_bc[i] = 0;
  }  

  int nghost[3] = {2, 2, 2};
  for (int dir=0; dir<app->ndim; ++dir) {
    if (is_np[dir]) {
      const enum gkyl_field_bc_type *bc;
      if (dir == 0)
        bc = mom_fld->bcx;
      else if (dir == 1)
        bc = mom_fld->bcy;
      else
        bc = mom_fld->bcz;

      void (*bc_lower_func)(double t, int nc, const double *skin, double * GKYL_RESTRICT ghost, void *ctx);
      if (dir == 0)
        bc_lower_func = mom_fld->bcx_lower_func;
      else if (dir == 1)
        bc_lower_func = mom_fld->bcy_lower_func;
      else
        bc_lower_func = mom_fld->bcz_lower_func;

      void (*bc_upper_func)(double t, int nc, const double *skin, double * GKYL_RESTRICT ghost, void *ctx);
      if (dir == 0)
        bc_upper_func = mom_fld->bcx_upper_func;
      else if (dir == 1)
        bc_upper_func = mom_fld->bcy_upper_func;
      else
        bc_upper_func = mom_fld->bcz_upper_func;


      fld->lower_bct[dir] = bc[0];
      fld->upper_bct[dir] = bc[1];

      switch (bc[0]) {
        case GKYL_FIELD_PEC_WALL:
          fld->lower_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, maxwell, app->geom, dir, GKYL_LOWER_EDGE, nghost, maxwell->wall_bc_func, 0);
          break;

        case GKYL_FIELD_COPY:
        case GKYL_FIELD_WEDGE:
          fld->lower_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, maxwell, app->geom, dir, GKYL_LOWER_EDGE, nghost, bc_copy, 0);
          break;

        case GKYL_FIELD_FUNC:
          fld->lower_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, maxwell, app->geom, dir, GKYL_LOWER_EDGE, nghost,
            bc_lower_func, mom_fld->ctx);
          break;
      }

      switch (bc[1]) {
        case GKYL_FIELD_PEC_WALL:
          fld->upper_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, maxwell, app->geom, dir, GKYL_UPPER_EDGE, nghost, maxwell->wall_bc_func, 0);
          break;

        case GKYL_FIELD_COPY:
        case GKYL_FIELD_WEDGE:
          fld->upper_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, maxwell, app->geom, dir, GKYL_UPPER_EDGE, nghost, bc_copy, 0);
          break;

        case GKYL_FIELD_FUNC:
          fld->upper_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, maxwell, app->geom, dir, GKYL_UPPER_EDGE, nghost,
            bc_upper_func, mom_fld->ctx);
          break;
      }
    }
  }

  // allocate arrays for applied current/external fields
  fld->app_current = mkarr(false, 3, app->local_ext.volume);
  fld->proj_app_current = 0;
  if (mom_fld->app_current_func)
    fld->proj_app_current = gkyl_fv_proj_new(&app->grid, 2, 3, mom_fld->app_current_func, fld->ctx);

  
  fld->ext_em = mkarr(false, 6, app->local_ext.volume);
  fld->is_ext_em_static = mom_fld->is_ext_em_static;

  fld->was_ext_em_computed = false;
  fld->proj_ext_em = 0;
  if (mom_fld->ext_em_func)
    fld->proj_ext_em = gkyl_fv_proj_new(&app->grid, 2, 6, mom_fld->ext_em_func, fld->ctx);

  // allocate buffer for applying BCs (used for periodic BCs)
  long buff_sz = 0;
  // compute buffer size needed
  for (int d=0; d<app->ndim; ++d) {
    long vol = app->skin_ghost.lower_skin[d].volume;
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  fld->bc_buffer = mkarr(false, 8, buff_sz);

  gkyl_wv_eqn_release(maxwell);

  fld->integ_energy = gkyl_dynvec_new(GKYL_DOUBLE, 6);
  fld->is_first_energy_write_call = true;
}

// apply BCs to EM field
void
moment_field_apply_bc(const gkyl_moment_app *app, double tcurr,
  const struct moment_field *field, struct gkyl_array *f)
{
  int num_periodic_dir = app->num_periodic_dir, ndim = app->ndim, is_non_periodic[3] = {1, 1, 1};
  for (int d=0; d<num_periodic_dir; ++d) {
    moment_apply_periodic_bc(app, field->bc_buffer, app->periodic_dirs[d], f);
    is_non_periodic[app->periodic_dirs[d]] = 0;
  }

  for (int d=0; d<ndim; ++d)
    if (is_non_periodic[d]) {
      // handle non-wedge BCs
      if (field->lower_bct[d] != GKYL_FIELD_WEDGE)
        gkyl_wv_apply_bc_advance(field->lower_bc[d], tcurr, &app->local, f);
      if (field->upper_bct[d] != GKYL_FIELD_WEDGE)      
        gkyl_wv_apply_bc_advance(field->upper_bc[d], tcurr, &app->local, f);

      // wedge BCs for upper/lower must be handled in one shot
      if (field->lower_bct[d] == GKYL_FIELD_WEDGE)
        moment_apply_wedge_bc(app, tcurr, &app->local,
          field->bc_buffer, d, field->lower_bc[d], field->upper_bc[d], f);
    }

  gkyl_comm_gkyl_array_sync(app->comm, f);
}

double
moment_field_max_dt(const gkyl_moment_app *app, const struct moment_field *fld)
{
  double max_dt = DBL_MAX;
  if (fld->scheme_type == GKYL_MOMENT_WAVE_PROP) {  
    for (int d=0; d<app->ndim; ++d)
      max_dt = fmin(max_dt, gkyl_wave_prop_max_dt(fld->slvr[d], &app->local, fld->f[0]));
  }
  else if (fld->scheme_type == GKYL_MOMENT_MP || fld->scheme_type == GKYL_MOMENT_KEP) {
    max_dt = fmin(max_dt, gkyl_mp_scheme_max_dt(fld->mp_slvr, &app->local, fld->f0));
  }
  return max_dt;
}

// update solution: initial solution is in fld->f[0] and updated
// solution in fld->f[ndim]
struct gkyl_update_status
moment_field_update(const gkyl_moment_app *app,
  const struct moment_field *fld, double tcurr, double dt)
{
  int ndim = fld->ndim;
  struct gkyl_wave_prop_status stat = { true, DBL_MAX };

  for (int d=0; d<ndim; ++d) {
    // update solution
    stat = gkyl_wave_prop_advance(fld->slvr[d], tcurr, dt, &app->local, fld->f[d], fld->f[d+1]);

    if (!stat.success)
      return (struct gkyl_update_status) {
        .success = false,
        .dt_suggested = stat.dt_suggested
      };
    // apply BC
    moment_field_apply_bc(app, tcurr, fld, fld->f[d+1]);
  }

  return (struct gkyl_update_status) {
    .success = true,
    .dt_suggested = stat.dt_suggested
  };
}

// Compute RHS of EM equations
double
moment_field_rhs(gkyl_moment_app *app, struct moment_field *fld,
  const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec tm = gkyl_wall_clock();
  
  gkyl_array_clear(fld->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  gkyl_mp_scheme_advance(fld->mp_slvr, &app->local, fin,
    app->ql, app->qr, app->amdq, app->apdq,
    fld->cflrate, rhs);

  double omegaCfl[1];
  gkyl_array_reduce_range(omegaCfl, fld->cflrate, GKYL_MAX, app->local);

  app->stat.field_rhs_tm += gkyl_time_diff_now_sec(tm);
  
  return app->cfl/omegaCfl[0];
}

// free field
void
moment_field_release(const struct moment_field *fld)
{
  gkyl_wv_eqn_release(fld->maxwell);
  
  for (int d=0; d<fld->ndim; ++d) {
    if (fld->lower_bc[d])
      gkyl_wv_apply_bc_release(fld->lower_bc[d]);
    if (fld->upper_bc[d])    
      gkyl_wv_apply_bc_release(fld->upper_bc[d]);
  }

  if (fld->scheme_type == GKYL_MOMENT_WAVE_PROP) {
    for (int d=0; d<fld->ndim; ++d)
      gkyl_wave_prop_release(fld->slvr[d]);
    
    gkyl_array_release(fld->fdup);
    for (int d=0; d<fld->ndim+1; ++d)
      gkyl_array_release(fld->f[d]);
  }
  else if (fld->scheme_type == GKYL_MOMENT_MP) {
    gkyl_mp_scheme_release(fld->mp_slvr);
    gkyl_array_release(fld->f0);
    gkyl_array_release(fld->f1);
    gkyl_array_release(fld->fnew);
    gkyl_array_release(fld->cflrate);
  }
    
  gkyl_array_release(fld->app_current);
  if (fld->proj_app_current)
    gkyl_fv_proj_release(fld->proj_app_current);
  
  gkyl_array_release(fld->ext_em);
  if (fld->proj_ext_em)
    gkyl_fv_proj_release(fld->proj_ext_em);

  gkyl_dynvec_release(fld->integ_energy);
  gkyl_array_release(fld->bc_buffer);
}

