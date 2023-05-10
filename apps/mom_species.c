#include <gkyl_moment_priv.h>
#include <gkyl_util.h>
#include <gkyl_wv_euler.h>

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
  sp->equation = gkyl_wv_eqn_acquire(mom_sp->equation);
  
  // closure parameter, used by 10 moment
  sp->k0 = mom_sp->equation->type == GKYL_EQN_TEN_MOMENT ? gkyl_wv_ten_moment_k0(mom_sp->equation) : 0.0;
  // check if we are running with gradient-based closure
  sp->has_grad_closure = mom_sp->has_grad_closure == 0 ? 0 : mom_sp->has_grad_closure;

  sp->scheme_type = mom->scheme_type;

  // choose default limiter
  enum gkyl_wave_limiter limiter =
    mom_sp->limiter == 0 ? GKYL_MONOTONIZED_CENTERED : mom_sp->limiter;

  int ndim = mom->ndim;
  int meqn = sp->num_equations;  

  if (sp->scheme_type == GKYL_MOMENT_WAVE_PROP) {
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
          .comm = app->comm
        }
      );
      
    sp->fdup = mkarr(false, meqn, app->local_ext.volume);
    // allocate arrays
    for (int d=0; d<ndim+1; ++d)
      sp->f[d] = mkarr(false, meqn, app->local_ext.volume);
  
    // set current solution so ICs and IO work properly
    sp->fcurr = sp->f[0];
  }
  else if ( sp->scheme_type == GKYL_MOMENT_MP || sp->scheme_type == GKYL_MOMENT_KEP ) {
    // determine directions to update
    int num_up_dirs = 0, update_dirs[GKYL_MAX_CDIM] = { 0 };
    for (int d=0; d<ndim; ++d)
      if (!app->is_dir_skipped[d]) {
        update_dirs[num_up_dirs] = d;
        num_up_dirs += 1;
      }

    if (sp->scheme_type == GKYL_MOMENT_MP)
      // single MP updater updates all directions
      sp->mp_slvr = gkyl_mp_scheme_new( &(struct gkyl_mp_scheme_inp) {
          .grid = &app->grid,
          .equation = mom_sp->equation,
          .mp_recon = app->mp_recon,
          .skip_mp_limiter = mom->skip_mp_limiter,
          .num_up_dirs = num_up_dirs,
          .update_dirs = { update_dirs[0], update_dirs[1], update_dirs[2] } ,
          .cfl = app->cfl,
          .geom = app->geom,
        }
      );
    else
      // single KEP updater updates all directions
      sp->kep_slvr = gkyl_kep_scheme_new( &(struct gkyl_kep_scheme_inp) {
          .grid = &app->grid,
          .equation = mom_sp->equation,
          .use_hybrid_flux = app->use_hybrid_flux_kep,
          .num_up_dirs = num_up_dirs,
          .update_dirs = { update_dirs[0], update_dirs[1], update_dirs[2] } ,
          .cfl = app->cfl,
          .geom = app->geom,
        }
      );
    
    // allocate arrays
    sp->f0 = mkarr(false, meqn, app->local_ext.volume);
    sp->f1 = mkarr(false, meqn, app->local_ext.volume);
    sp->fnew = mkarr(false, meqn, app->local_ext.volume);
    sp->cflrate = mkarr(false, 1, app->local_ext.volume);
    sp->alpha = mkarr(false, 1, app->local_ext.volume);
    
    // set current solution so ICs and IO work properly
    sp->fcurr = sp->f0;
  }

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

      void (*bc_lower_func)(double t, int nc, const double *skin, double * GKYL_RESTRICT ghost, void *ctx);
      if (dir == 0)
        bc_lower_func = mom_sp->bcx_lower_func;
      else if (dir == 1)
        bc_lower_func = mom_sp->bcy_lower_func;
      else
        bc_lower_func = mom_sp->bcz_lower_func;

      void (*bc_upper_func)(double t, int nc, const double *skin, double * GKYL_RESTRICT ghost, void *ctx);
      if (dir == 0)
        bc_upper_func = mom_sp->bcx_upper_func;
      else if (dir == 1)
        bc_upper_func = mom_sp->bcy_upper_func;
      else
        bc_upper_func = mom_sp->bcz_upper_func;

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
            bc_lower_func, mom_sp->ctx);
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
            bc_upper_func, mom_sp->ctx);
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

  if (mom_sp->equation->type == GKYL_EQN_EULER)
    sp->integ_q = gkyl_dynvec_new(GKYL_DOUBLE, 6); // KE and PE are stored independently
  else
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
  if (ndim == 2)
    moment_apply_periodic_corner_sync_2d(app, f);
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

  gkyl_comm_array_sync(app->comm, &app->local, &app->local_ext, f);
}

// maximum stable time-step
double
moment_species_max_dt(const gkyl_moment_app *app, const struct moment_species *sp)
{
  double max_dt = DBL_MAX;
  if (sp->scheme_type == GKYL_MOMENT_WAVE_PROP) {
    for (int d=0; d<app->ndim; ++d)
      max_dt = fmin(max_dt, gkyl_wave_prop_max_dt(sp->slvr[d], &app->local, sp->f[0]));
  }
  else if (sp->scheme_type == GKYL_MOMENT_MP) {
    max_dt = fmin(max_dt, gkyl_mp_scheme_max_dt(sp->mp_slvr, &app->local, sp->f0));
  }
  else if (sp->scheme_type == GKYL_MOMENT_KEP) {
    max_dt = fmin(max_dt, gkyl_kep_scheme_max_dt(sp->kep_slvr, &app->local, sp->f0));
  }  
  return max_dt;
}

// update solution: initial solution is in sp->f[0] and updated
// solution in sp->f[ndim]
struct gkyl_update_status
moment_species_update(const gkyl_moment_app *app,
  struct moment_species *sp, double tcurr, double dt)
{
  int ndim = sp->ndim;
  double dt_suggested = DBL_MAX;
  double max_speed = 0.0;
  struct gkyl_wave_prop_status stat;

  for (int d=0; d<ndim; ++d) {
    stat = gkyl_wave_prop_advance(sp->slvr[d], tcurr, dt, &app->local, sp->f[d], sp->f[d+1]);

    double my_max_speed = stat.max_speed;
    max_speed = max_speed > my_max_speed ? max_speed : my_max_speed;

    if (!stat.success)
      return (struct gkyl_update_status) {
        .success = false,
        .dt_suggested = stat.dt_suggested
      };
    
    dt_suggested = fmin(dt_suggested, stat.dt_suggested);
    moment_species_apply_bc(app, tcurr, sp, sp->f[d+1]);
  }

  for (int d=0; d<ndim; ++d) {
    struct gkyl_wv_eqn *eqn = sp->equation;
    if (eqn->type==GKYL_EQN_MHD) {
      if (sp->eqn_type==GKYL_MHD_DIVB_GLM) {
        gkyl_wv_mhd_set_glm_ch(eqn, max_speed);
      }
    }
  }
  if (app->update_mhd_source) {
    // app->mhd_source.slvr->glm_ch = max_speed;
    // gkyl_mhd_src_set_glm_ch(app->mhd_source.slvr, max_speed);
  }

  return (struct gkyl_update_status) {
    .success = true,
    .dt_suggested = dt_suggested
  };
}

// Compute RHS of moment equations
double
moment_species_rhs(gkyl_moment_app *app, struct moment_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec tm = gkyl_wall_clock();
  
  gkyl_array_clear(species->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  if (app->scheme_type == GKYL_MOMENT_MP)
    gkyl_mp_scheme_advance(species->mp_slvr, &app->local, fin,
      app->ql, app->qr, app->amdq, app->apdq,
      species->cflrate, rhs);
  else
    gkyl_kep_scheme_advance(species->kep_slvr, &app->local, fin, species->alpha,
      species->cflrate, rhs);

  double omegaCfl[1];
  gkyl_array_reduce_range(omegaCfl, species->cflrate, GKYL_MAX, app->local);

  app->stat.species_rhs_tm += gkyl_time_diff_now_sec(tm);
  
  return app->cfl/omegaCfl[0];
}

// free species
void
moment_species_release(const struct moment_species *sp)
{
  gkyl_wv_eqn_release(sp->equation);
  
  for (int d=0; d<sp->ndim; ++d) {
    if (sp->lower_bc[d])
      gkyl_wv_apply_bc_release(sp->lower_bc[d]);
    if (sp->upper_bc[d])    
      gkyl_wv_apply_bc_release(sp->upper_bc[d]);
  }

  if (sp->scheme_type == GKYL_MOMENT_WAVE_PROP) {
    for (int d=0; d<sp->ndim; ++d)
      gkyl_wave_prop_release(sp->slvr[d]);
    
    gkyl_array_release(sp->fdup);
    for (int d=0; d<sp->ndim+1; ++d)
      gkyl_array_release(sp->f[d]);
  }
  else if (sp->scheme_type == GKYL_MOMENT_MP || sp->scheme_type == GKYL_MOMENT_KEP) {

    if (sp->scheme_type == GKYL_MOMENT_MP)
      gkyl_mp_scheme_release(sp->mp_slvr);
    else
      gkyl_kep_scheme_release(sp->kep_slvr);
    
    gkyl_array_release(sp->f0);
    gkyl_array_release(sp->f1);
    gkyl_array_release(sp->fnew);
    gkyl_array_release(sp->cflrate);
    gkyl_array_release(sp->alpha);
  }

  gkyl_array_release(sp->app_accel);
  if (sp->proj_app_accel)
    gkyl_fv_proj_release(sp->proj_app_accel);

  gkyl_array_release(sp->bc_buffer);

  gkyl_dynvec_release(sp->integ_q);
}

/** mhd_src functions */

void
mhd_src_init(const struct gkyl_moment_app *app,
             const struct gkyl_moment_species *sp,
             struct mhd_src *src)
{
  double dxyz_min = DBL_MAX;
  for (int d=0; d<app->grid.ndim; ++d) {
    double dx = app->grid.dx[d];
    dxyz_min = dx < dxyz_min ? dx : dxyz_min;
  }

  struct gkyl_mhd_src_inp src_inp = {
    .grid = &app->grid,
    .divergence_constraint = gkyl_wv_mhd_divergence_constraint(sp->equation),
    .glm_ch = gkyl_wv_mhd_glm_ch(sp->equation),
    .glm_alpha = gkyl_wv_mhd_glm_ch(sp->equation),
    .dxyz_min = dxyz_min,
  };

  src->slvr = gkyl_mhd_src_new(src_inp, &app->local_ext);
}

// update sources: 'nstrang' is 0 for the first Strang step and 1 for
// the second step
void
mhd_src_update(gkyl_moment_app *app, struct mhd_src *src, int nstrang,
               double tcurr, double dt)
{
  int sidx[] = { 0, app->ndim };
  int i = 0; // mhd has only one 'species'
  struct gkyl_array *fluid = app->species[i].f[sidx[nstrang]];

  if (app->species[i].proj_app_accel)
    gkyl_fv_proj_advance(app->species[i].proj_app_accel, tcurr, &app->local,
                         app->species[i].app_accel);

  // FIXME presently needed for computing divB etc
  moment_species_apply_bc(app, tcurr, &app->species[i], fluid);

  gkyl_mhd_src_advance(src->slvr, dt, &app->local, fluid,
                       app->species[i].app_accel);

  moment_species_apply_bc(app, tcurr, &app->species[i], fluid);

}

void
mhd_src_release(const struct mhd_src *src)
{
  gkyl_mhd_src_release(src->slvr);
}
