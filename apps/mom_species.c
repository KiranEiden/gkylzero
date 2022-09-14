#include <gkyl_moment_priv.h>

// initialize species
void
moment_species_init(const struct gkyl_moment *mom, const struct gkyl_moment_species *mom_sp,
  struct gkyl_moment_app *app, struct moment_species *sp)
{
  sp->ndim = mom->ndim;
  strcpy(sp->name, mom_sp->name);
  sp->charge = mom_sp->charge;
  sp->mass = mom_sp->mass;
  sp->ctx = mom_sp->ctx;
  sp->init = mom_sp->init;

  sp->eqn_type = mom_sp->equation->type;
  sp->num_equations = mom_sp->equation->num_equations;
  // closure parameter, used by 10 moment
  sp->k0 = mom_sp->equation->type == GKYL_EQN_TEN_MOMENT ? gkyl_wv_ten_moment_k0(mom_sp->equation) : 0.0;

  // choose default limiter
  enum gkyl_wave_limiter limiter =
    mom_sp->limiter == 0 ? GKYL_MONOTONIZED_CENTERED : mom_sp->limiter;

  int ndim = mom->ndim;
  // create updaters for each directional update
  for (int d=0; d<ndim; ++d)
    sp->slvr[d] = gkyl_wave_prop_new( &(struct gkyl_wave_prop_inp) {
        .grid = &app->grid,
        .equation = mom_sp->equation,
        .limiter = limiter,
        .num_up_dirs = app->is_dir_skipped[d] ? 0 : 1,
        .force_low_order_flux = mom_sp->force_low_order_flux,
        .check_inv_domain = true,
        .update_dirs = { d },
        .cfl = app->cfl,
        .geom = app->geom,
      }
    );

  // determine which directions are not periodic
  int num_periodic_dir = app->num_periodic_dir, is_np[3] = {1, 1, 1};
  for (int d=0; d<num_periodic_dir; ++d)
    is_np[app->periodic_dirs[d]] = 0;

  for (int i=0; i<3; ++i) {
    sp->lower_bc[i] = 0;
    sp->upper_bc[i] = 0;
  }

  int nghost[3] = {2, 2, 2};
  for (int dir=0; dir<app->ndim; ++dir) {
    if (is_np[dir]) {
      const enum gkyl_species_bc_type *bc;
      if (dir == 0)
        bc = mom_sp->bcx;
      else if (dir == 1)
        bc = mom_sp->bcy;
      else
        bc = mom_sp->bcz;

      sp->lower_bct[dir] = bc[0];
      sp->upper_bct[dir] = bc[1];

      // lower BCs
      switch (bc[0]) {
        case GKYL_SPECIES_REFLECT:
          sp->lower_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_LOWER_EDGE, nghost,
            mom_sp->equation->wall_bc_func, 0);
          break;

        case GKYL_SPECIES_NO_SLIP:
          sp->lower_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_LOWER_EDGE, nghost,
            mom_sp->equation->no_slip_bc_func, 0);
          break;

        case GKYL_SPECIES_FUNC:
          sp->lower_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_LOWER_EDGE, nghost,
            mom_sp->bc_lower_func, mom_sp->ctx);
          break;
        
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_WEDGE: // wedge also uses bc_copy
          sp->lower_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_LOWER_EDGE, nghost, bc_copy, 0);
          break;

        default:
          assert(false);
          break;
      }
      
      // upper BCs
      switch (bc[1]) {
        case GKYL_SPECIES_REFLECT:      
          sp->upper_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_UPPER_EDGE, nghost,
            mom_sp->equation->wall_bc_func, 0);
          break;

        case GKYL_SPECIES_NO_SLIP:      
          sp->upper_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_UPPER_EDGE, nghost,
            mom_sp->equation->no_slip_bc_func, 0);
          break;

        case GKYL_SPECIES_FUNC:
          sp->upper_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_UPPER_EDGE, nghost,
            mom_sp->bc_upper_func, mom_sp->ctx);
          break;
          
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_WEDGE:
          sp->upper_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_UPPER_EDGE, nghost, bc_copy, 0);
          break;

        default:
          assert(false);
          break;
      }
    }
  }

  int meqn = sp->num_equations;
  sp->fdup = mkarr(false, meqn, app->local_ext.volume);
  // allocate arrays
  for (int d=0; d<ndim+1; ++d)
    sp->f[d] = mkarr(false, meqn, app->local_ext.volume);

  // allocate array for applied acceleration/forces for each species
  sp->app_accel = mkarr(false, 3, app->local_ext.volume);
  sp->proj_app_accel = 0;
  if (mom_sp->app_accel_func)
    sp->proj_app_accel = gkyl_fv_proj_new(&app->grid, 2, 3, mom_sp->app_accel_func, sp->ctx);
  // allocate buffer for applying BCs (used for periodic BCs)
  long buff_sz = 0;
  // compute buffer size needed
  for (int d=0; d<app->ndim; ++d) {
    long vol = app->skin_ghost.lower_skin[d].volume;
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  sp->bc_buffer = mkarr(false, meqn, buff_sz);

  sp->integ_q = gkyl_dynvec_new(GKYL_DOUBLE, meqn);
  sp->is_first_q_write_call = true;
}

// apply BCs to species
void
moment_species_apply_bc(const gkyl_moment_app *app, double tcurr,
  const struct moment_species *sp, struct gkyl_array *f)
{
  int num_periodic_dir = app->num_periodic_dir, ndim = app->ndim, is_non_periodic[3] = {1, 1, 1};
  for (int d=0; d<num_periodic_dir; ++d) {
    moment_apply_periodic_bc(app, sp->bc_buffer, app->periodic_dirs[d], f);
    is_non_periodic[app->periodic_dirs[d]] = 0;
  }
  
  for (int d=0; d<ndim; ++d)
    if (is_non_periodic[d]) {
      // handle non-wedge BCs
      if (sp->lower_bct[d] != GKYL_SPECIES_WEDGE)
        gkyl_wv_apply_bc_advance(sp->lower_bc[d], tcurr, &app->local, f);
      if (sp->upper_bct[d] != GKYL_SPECIES_WEDGE)
        gkyl_wv_apply_bc_advance(sp->upper_bc[d], tcurr, &app->local, f);

      // wedge BCs for upper/lower must be handled in one shot
      if (sp->lower_bct[d] == GKYL_SPECIES_WEDGE)
        moment_apply_wedge_bc(app, tcurr, &app->local,
          sp->bc_buffer, d, sp->lower_bc[d], sp->upper_bc[d], f);
    }
}

// maximum stable time-step
double
moment_species_max_dt(const gkyl_moment_app *app, const struct moment_species *sp)
{
  double max_dt = DBL_MAX;
  for (int d=0; d<app->ndim; ++d)
    max_dt = fmin(max_dt, gkyl_wave_prop_max_dt(sp->slvr[d], &app->local, sp->f[0]));
  return max_dt;
}

// update solution: initial solution is in sp->f[0] and updated
// solution in sp->f[ndim]
struct gkyl_update_status
moment_species_update(const gkyl_moment_app *app,
  const struct moment_species *sp, double tcurr, double dt)
{
  int ndim = sp->ndim;
  double dt_suggested = DBL_MAX;
  struct gkyl_wave_prop_status stat;

  for (int d=0; d<ndim; ++d) {
    stat = gkyl_wave_prop_advance(sp->slvr[d], tcurr, dt, &app->local, sp->f[d], sp->f[d+1]);

    if (!stat.success)
      return (struct gkyl_update_status) {
        .success = false,
        .dt_suggested = stat.dt_suggested
      };
    
    dt_suggested = fmin(dt_suggested, stat.dt_suggested);
    moment_species_apply_bc(app, tcurr, sp, sp->f[d+1]);
  }

  return (struct gkyl_update_status) {
    .success = true,
    .dt_suggested = dt_suggested
  };
}

// free species
void
moment_species_release(const struct moment_species *sp)
{
  for (int d=0; d<sp->ndim; ++d)
    gkyl_wave_prop_release(sp->slvr[d]);

  for (int d=0; d<sp->ndim; ++d) {
    if (sp->lower_bc[d])
      gkyl_wv_apply_bc_release(sp->lower_bc[d]);
    if (sp->upper_bc[d])    
      gkyl_wv_apply_bc_release(sp->upper_bc[d]);
  }

  gkyl_array_release(sp->fdup);
  for (int d=0; d<sp->ndim+1; ++d)
    gkyl_array_release(sp->f[d]);

  gkyl_array_release(sp->app_accel);
  if (sp->proj_app_accel)
    gkyl_fv_proj_release(sp->proj_app_accel);

  gkyl_array_release(sp->bc_buffer);

  gkyl_dynvec_release(sp->integ_q);
}
