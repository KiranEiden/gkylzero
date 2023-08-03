#include <assert.h>
#include <float.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_priv.h>

// initialize fluid species object
void
vm_fluid_species_init(struct gkyl_vm *vm, struct gkyl_vlasov_app *app, struct vm_fluid_species *f)
{
  int cdim = app->cdim;
  int vdim = app->vdim;
  int num_eqn = f->info.num_eqn ? f->info.num_eqn : 1;
  // allocate fluid arrays
  f->fluid = mkarr(app->use_gpu, num_eqn*app->confBasis.num_basis, app->local_ext.volume);
  f->fluid1 = mkarr(app->use_gpu, num_eqn*app->confBasis.num_basis, app->local_ext.volume);
  f->fluidnew = mkarr(app->use_gpu, num_eqn*app->confBasis.num_basis, app->local_ext.volume);

  f->fluid_host = f->fluid;
  if (app->use_gpu)
    f->fluid_host = mkarr(false, num_eqn*app->confBasis.num_basis, app->local_ext.volume);

  // allocate buffer for applying BCs
  long buff_sz = 0;
  // compute buffer size needed
  for (int d=0; d<app->cdim; ++d) {
    long vol = app->skin_ghost.lower_skin[d].volume;
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  f->bc_buffer = mkarr(app->use_gpu, num_eqn*app->confBasis.num_basis, buff_sz);

  // allocate cflrate (scalar array)
  f->cflrate = mkarr(app->use_gpu, 1, app->local_ext.volume);
  if (app->use_gpu)
    f->omegaCfl_ptr = gkyl_cu_malloc(sizeof(double));
  else
    f->omegaCfl_ptr = gkyl_malloc(sizeof(double));

  int up_dirs[GKYL_MAX_DIM] = {0, 1, 2}, zero_flux_flags[GKYL_MAX_DIM] = {0, 0, 0};

  // initial pointers to primitive variables, pressure, and boolean array for if we are only using the cell average for primitive variables
  // For isothermal Euler, prim : (ux, uy, uz), p : (vth*rho)
  // For Euler, prim : (ux, uy, uz, T/m), p : (gamma - 1)*(E - 1/2 rho u^2)
  // For PKPM, prim : [ux, uy, uz, 3*Txx/m, 3*Tyy/m, 3*Tzz/m, 1/rho*div(p_par b), T_perp/m, m/T_perp]
  // p_ij : (p_par - p_perp) b_i b_j + p_perp g_ij
  f->prim = 0;
  f->p = 0;
  f->cell_avg_prim = 0;

  // initialize pointer to pkpm acceleration variables, stored in pkpm_accel: 
  // 0: div_b (divergence of magnetic field unit vector)
  // 1: bb_grad_u (bb : grad(u))
  // 2: p_force (total pressure forces in kinetic equation 1/rho div(p_parallel b_hat) - T_perp/m*div(b)
  // 3: p_perp_source (pressure source for higher Laguerre moments -> bb : grad(u) - div(u) - 2*nu)
  // 4: p_perp_div_b (p_perp/rho*div(b) = T_perp/m*div(b))
  f->pkpm_accel = 0;

  // initialize pointers for io.
  // For isothermal Euler and Euler, these are the same as the state variables
  // For PKPM we construct the 10 moment variables for ease of analysis 
  // along with an array of the various update variables, primitive and acceleration
  f->fluid_io = 0;
  f->fluid_io_host = 0;
  f->pkpm_vars_io = 0;
  f->pkpm_vars_io_host = 0;

  f->param = 0.0;
  // fluid solvers
  if (f->info.vt) {
    f->param = f->info.vt; // parameter for isothermal Euler is vt, thermal velocity
    f->eqn_id = GKYL_EQN_ISO_EULER;
    // allocate array to store primitive variables (ux, uy, uz) and pressure (vth*rho)
    f->prim = mkarr(app->use_gpu, 3*app->confBasis.num_basis, app->local_ext.volume);
    f->p = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    // boolean array for if we are only using the cell average for primitive variables
    f->cell_avg_prim = mk_int_arr(app->use_gpu, 1, app->local_ext.volume);

    f->fluid_io = f->fluid;
    f->fluid_io_host = f->fluid_host;
  }
  else if (f->info.gas_gamma) {
    f->param = f->info.gas_gamma; // parameter for Euler is gas_gamma, adiabatic index
    f->eqn_id = GKYL_EQN_EULER;
    // allocate array to store primitive variables (ux, uy, uz, T/m) and pressure (gamma - 1)*(E - 1/2 rho u^2)
    f->prim = mkarr(app->use_gpu, 4*app->confBasis.num_basis, app->local_ext.volume);
    f->p = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    // boolean array for if we are only using the cell average for primitive variables
    f->cell_avg_prim = mk_int_arr(app->use_gpu, 1, app->local_ext.volume);

    f->fluid_io = f->fluid;
    f->fluid_io_host = f->fluid_host;
  }
  else if (f->info.advection.velocity) {
    f->eqn_id = GKYL_EQN_ADVECTION;

    // setup applied advection 
    f->app_advect = mkarr(app->use_gpu, cdim*app->confBasis.num_basis, app->local_ext.volume);

    f->app_advect_host = f->app_advect;
    if (app->use_gpu)
      f->app_advect_host = mkarr(false, cdim*app->confBasis.num_basis, app->local_ext.volume);

    gkyl_proj_on_basis *app_advect_proj = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
        .grid = &app->grid,
        .basis = &app->confBasis,
        .qtype = GKYL_GAUSS_LOBATTO_QUAD,
        .num_quad = 8,
        .num_ret_vals = cdim,
        .eval = f->info.advection.velocity,
        .ctx = f->info.advection.velocity_ctx
      }
    );

    gkyl_proj_on_basis_advance(app_advect_proj, 0.0, &app->local_ext, f->app_advect_host);
    if (app->use_gpu) // note: app_advect_host is same as advect when not on GPUs
      gkyl_array_copy(f->app_advect, f->app_advect_host);
    // Free projection object
    gkyl_proj_on_basis_release(app_advect_proj);

    f->fluid_io = f->fluid;
    f->fluid_io_host = f->fluid_host;
  }
  else {
    f->eqn_id = GKYL_EQN_EULER_PKPM;
    // allocate array to store primitive moments : [ux, uy, uz, 3*Txx/m, 3*Tyy/m, 3*Tzz/m, 1/rho*div(p_par b), T_perp/m, m/T_perp]
    // and pressure p_ij : (p_par - p_perp) b_i b_j + p_perp g_ij
    f->prim = mkarr(app->use_gpu, 9*app->confBasis.num_basis, app->local_ext.volume);
    f->p = mkarr(app->use_gpu, 6*app->confBasis.num_basis, app->local_ext.volume);
    // boolean array for if we are only using the cell average for primitive variables
    f->cell_avg_prim = mk_int_arr(app->use_gpu, 1, app->local_ext.volume);

    // allocate array for pkpm acceleration variables, stored in pkpm_accel: 
    // 0: div_b (divergence of magnetic field unit vector)
    // 1: bb_grad_u (bb : grad(u))
    // 2: p_force (total pressure forces in kinetic equation 1/rho div(p_parallel b_hat) - T_perp/m*div(b)
    // 3: p_perp_source (pressure source for higher Laguerre moments -> bb : grad(u) - div(u) - nu + nu rho vth^2/p_perp)
    // 4: p_perp_div_b (p_perp/rho*div(b) = T_perp/m*div(b))
    f->pkpm_accel = mkarr(app->use_gpu, 5*app->confBasis.num_basis, app->local_ext.volume);

    f->fluid_io = mkarr(app->use_gpu, 10*app->confBasis.num_basis, app->local_ext.volume);
    f->pkpm_vars_io = mkarr(app->use_gpu, 10*app->confBasis.num_basis, app->local_ext.volume);
    f->fluid_io_host = f->fluid_io;
    f->pkpm_vars_io_host = f->pkpm_vars_io;
    if (app->use_gpu) {
      f->fluid_io_host = mkarr(false, 10*app->confBasis.num_basis, app->local_ext.volume);
      f->pkpm_vars_io_host = mkarr(false, 10*app->confBasis.num_basis, app->local_ext.volume);
    }

    // updater for computing pkpm variables 
    // pressure, primitive variables, and acceleration variables
    // also stores kernels for computing source terms, integrated variables
    // Two instances, one over extended range and one over local range for ease of handling boundary conditions
    f->calc_pkpm_vars_ext = gkyl_dg_calc_pkpm_vars_new(&app->grid, &app->confBasis, &app->local_ext, app->use_gpu);
    f->calc_pkpm_vars = gkyl_dg_calc_pkpm_vars_new(&app->grid, &app->confBasis, &app->local, app->use_gpu);

    f->pkpm_species = vm_find_species(app, f->info.pkpm_species);
    // index in fluid_species struct of fluid species kinetic species is colliding with
    f->species_index = vm_find_species_idx(app, f->info.pkpm_species);
  }

  f->advect_slvr = gkyl_dg_updater_fluid_new(&app->grid, &app->confBasis,
    &app->local, f->eqn_id, f->param, app->use_gpu);

  f->has_diffusion = false;
  f->Dij = 0;
  f->Dij_host = 0;
  f->diffusion_id = GKYL_NO_DIFFUSION;
  if (f->info.diffusion.Dij) {
    f->has_diffusion = true;
    // allocate space for full diffusion tensor 
    f->diffusion_id = GKYL_GEN_DIFFUSION;
    int szD = cdim*(cdim+1)/2;

    f->Dij = mkarr(app->use_gpu, szD*app->confBasis.num_basis, app->local_ext.volume);
    f->Dij_host = f->Dij;
    if (app->use_gpu)
      f->Dij_host = mkarr(false, szD*app->confBasis.num_basis, app->local_ext.volume);

    gkyl_proj_on_basis *diff_proj = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
        .grid = &app->grid,
        .basis = &app->confBasis,
        .qtype = GKYL_GAUSS_LOBATTO_QUAD,
        .num_quad = 8,
        .num_ret_vals = szD,
        .eval = f->info.diffusion.Dij,
        .ctx = f->info.diffusion.Dij_ctx
      }
    );
    gkyl_proj_on_basis_advance(diff_proj, 0.0, &app->local_ext, f->Dij_host);
    if (app->use_gpu) // note: Dij_host is same as Dij when not on GPUs
      gkyl_array_copy(f->Dij, f->Dij_host);
    // Free projection object
    gkyl_proj_on_basis_release(diff_proj);

    f->diff_slvr = gkyl_dg_updater_diffusion_new(&app->grid, &app->confBasis,
      &app->local, 0.0, 0, f->diffusion_id, app->use_gpu);
  }
  else if (f->info.diffusion.D) {
    f->has_diffusion = true;
    if (f->eqn_id == GKYL_EQN_ISO_EULER) 
      f->diffusion_id = GKYL_ISO_EULER_DIFFUSION;
    else if (f->eqn_id == GKYL_EQN_EULER) 
      f->diffusion_id = GKYL_EULER_DIFFUSION;
    else if (f->eqn_id == GKYL_EQN_EULER_PKPM) 
      f->diffusion_id = GKYL_PKPM_DIFFUSION;
    else 
      f->diffusion_id = GKYL_ISO_DIFFUSION;

    f->diff_slvr = gkyl_dg_updater_diffusion_new(&app->grid, &app->confBasis,
      &app->local, f->info.diffusion.D, f->info.diffusion.order, f->diffusion_id, app->use_gpu);    
  }

  // array for storing integrated moments in each cell
  f->integ_mom = mkarr(app->use_gpu, 6, app->local_ext.volume);
  if (app->use_gpu) {
    f->red_integ_diag = gkyl_cu_malloc(sizeof(double[6]));
  }
  // allocate dynamic-vector to store all-reduced integrated moments 
  f->integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, 6);
  f->is_first_integ_write_call = true;

  // set species source id
  f->source_id = f->info.source.source_id;

  // determine which directions are not periodic
  int num_periodic_dir = app->num_periodic_dir, is_np[3] = {1, 1, 1};
  for (int d=0; d<num_periodic_dir; ++d)
    is_np[app->periodic_dirs[d]] = 0;

  for (int dir=0; dir<app->cdim; ++dir) {
    f->lower_bc[dir] = f->upper_bc[dir] = GKYL_SPECIES_COPY;
    if (is_np[dir]) {
      const enum gkyl_species_bc_type *bc;
      if (dir == 0)
        bc = f->info.bcx;
      else if (dir == 1)
        bc = f->info.bcy;
      else
        bc = f->info.bcz;

      f->lower_bc[dir] = bc[0];
      f->upper_bc[dir] = bc[1];
    }
  }

  int ghost[GKYL_MAX_DIM] = {0.0};
  for (int d=0; d<app->cdim; ++d)
    ghost[d] = 1;

  // Certain operations fail if absorbing BCs used because absorbing BCs 
  // means the mass density is 0 in the ghost cells (divide by zero)
  f->bc_is_absorb = false;
  for (int d=0; d<app->cdim; ++d) {
    // Lower BC updater. Copy BCs by default.
    enum gkyl_bc_basic_type bctype = GKYL_BC_COPY;
    if (f->lower_bc[d] == GKYL_SPECIES_COPY) {
      bctype = GKYL_BC_COPY;
    }
    else if (f->lower_bc[d] == GKYL_SPECIES_ABSORB) {
      bctype = GKYL_BC_ABSORB;
      f->bc_is_absorb = true;
    }
    else if (f->lower_bc[d] == GKYL_SPECIES_REFLECT && f->eqn_id == GKYL_EQN_EULER_PKPM) {
      bctype = GKYL_BC_PKPM_MOM_REFLECT;
    }
    else if (f->lower_bc[d] == GKYL_SPECIES_NO_SLIP && f->eqn_id == GKYL_EQN_EULER_PKPM) {
      bctype = GKYL_BC_PKPM_MOM_NO_SLIP;
    }

    f->bc_lo[d] = gkyl_bc_basic_new(d, GKYL_LOWER_EDGE, &app->local_ext, ghost, bctype,
                                    app->basis_on_dev.confBasis, f->fluid->ncomp, app->cdim, app->use_gpu);

    // Upper BC updater. Copy BCs by default.
    if (f->upper_bc[d] == GKYL_SPECIES_COPY) {
      bctype = GKYL_BC_COPY;
    }
    else if (f->upper_bc[d] == GKYL_SPECIES_ABSORB) {
      bctype = GKYL_BC_ABSORB;
      f->bc_is_absorb = true;
    }
    else if (f->upper_bc[d] == GKYL_SPECIES_REFLECT && f->eqn_id == GKYL_EQN_EULER_PKPM) {
      bctype = GKYL_BC_PKPM_MOM_REFLECT;
    }
    else if (f->upper_bc[d] == GKYL_SPECIES_NO_SLIP && f->eqn_id == GKYL_EQN_EULER_PKPM) {
      bctype = GKYL_BC_PKPM_MOM_NO_SLIP;
    }

    f->bc_up[d] = gkyl_bc_basic_new(d, GKYL_UPPER_EDGE, &app->local_ext, ghost, bctype,
                                    app->basis_on_dev.confBasis, f->fluid->ncomp, app->cdim, app->use_gpu);
  }
}

void
vm_fluid_species_apply_ic(gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, double t0)
{
  int poly_order = app->poly_order;
  int num_eqn = fluid_species->info.num_eqn ? fluid_species->info.num_eqn : 1;
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
    poly_order+1, num_eqn, fluid_species->info.init, fluid_species->info.ctx);

  // run updater
  gkyl_proj_on_basis_advance(proj, t0, &app->local, fluid_species->fluid_host);
  gkyl_proj_on_basis_release(proj);

  if (app->use_gpu)
    gkyl_array_copy(fluid_species->fluid, fluid_species->fluid_host);

  // compute primitive variables at t = 0
  vm_fluid_species_prim_vars(app, fluid_species, fluid_species->fluid);

  // we are pre-computing source for now as it is time-independent
  vm_fluid_species_source_calc(app, fluid_species, t0);

}

void
vm_fluid_species_prim_vars(gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species,
  const struct gkyl_array *fluid)
{
  gkyl_array_clear(fluid_species->prim, 0.0);
  gkyl_array_clear(fluid_species->p, 0.0); 
  if (fluid_species->eqn_id == GKYL_EQN_EULER_PKPM) {
    gkyl_array_clear(fluid_species->pkpm_accel, 0.0);
    gkyl_dg_calc_pkpm_vars_pressure(fluid_species->calc_pkpm_vars, &app->local_ext, 
      app->field->bvar, fluid_species->pkpm_species->pkpm_moms.marr, 
      fluid_species->p);
    // Calculates both primitive and acceleration variables. Note acceleration variables
    // require gradients, which are computed with either averaging or recovery
    if (fluid_species->bc_is_absorb)
      gkyl_dg_calc_pkpm_vars_advance(fluid_species->calc_pkpm_vars,
        fluid_species->pkpm_species->pkpm_moms.marr, fluid, 
        fluid_species->p, fluid_species->pkpm_species->pkpm_div_ppar, 
        fluid_species->cell_avg_prim, fluid_species->prim); 
    else
      gkyl_dg_calc_pkpm_vars_advance(fluid_species->calc_pkpm_vars_ext,
        fluid_species->pkpm_species->pkpm_moms.marr, fluid, 
        fluid_species->p, fluid_species->pkpm_species->pkpm_div_ppar, 
        fluid_species->cell_avg_prim, fluid_species->prim); 
    gkyl_dg_calc_pkpm_vars_accel(fluid_species->calc_pkpm_vars, &app->local, 
      app->field->bvar, fluid_species->prim, fluid_species->pkpm_species->lbo.nu_sum, 
      fluid_species->pkpm_accel); 
  }
}

// Compute the RHS for fluid species update, returning maximum stable
// time-step.
double
vm_fluid_species_rhs(gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species,
  const struct gkyl_array *fluid, const struct gkyl_array *em, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();

  double omegaCfl = 1/DBL_MAX;

  gkyl_array_clear(fluid_species->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  if (fluid_species->eqn_id == GKYL_EQN_EULER_PKPM) {
    struct gkyl_dg_euler_pkpm_auxfields pkpm_inp = {.pkpm_prim = fluid_species->prim, 
      .p_ij = fluid_species->p, .vlasov_pkpm_moms = fluid_species->pkpm_species->pkpm_moms.marr};
    gkyl_dg_updater_fluid_advance(fluid_species->advect_slvr, 
      &app->local, &pkpm_inp, fluid, fluid_species->cflrate, rhs);
  }
  else if(fluid_species->eqn_id == GKYL_EQN_ADVECTION) {
    struct gkyl_dg_advection_auxfields adv_in = {.u_i = fluid_species->app_advect};
    gkyl_dg_updater_fluid_advance(fluid_species->advect_slvr, 
      &app->local, &adv_in, fluid, fluid_species->cflrate, rhs);
  }
  else {
    struct gkyl_dg_euler_auxfields euler_inp = {.u_i = fluid_species->prim, 
      .p_ij = fluid_species->p};
    gkyl_dg_updater_fluid_advance(fluid_species->advect_slvr, 
      &app->local, &euler_inp, fluid, fluid_species->cflrate, rhs);
  }

  if (fluid_species->has_diffusion)
    gkyl_dg_updater_diffusion_advance(fluid_species->diff_slvr, 
      &app->local, fluid_species->Dij, fluid, fluid_species->cflrate, rhs);

  // Accumulate source contribution if PKPM -> adds forces (E + u x B) to momentum equation RHS
  if (fluid_species->eqn_id == GKYL_EQN_EULER_PKPM) {
    double qbym = fluid_species->pkpm_species->info.charge/fluid_species->pkpm_species->info.mass;
    gkyl_array_set(fluid_species->pkpm_species->qmem, qbym, em);

    // Accumulate applied acceleration and/or q/m*(external electromagnetic)
    // fields onto qmem to get the total acceleration
    if (fluid_species->pkpm_species->has_accel)
      gkyl_array_accumulate(fluid_species->pkpm_species->qmem, 1.0, fluid_species->pkpm_species->accel);
    if (app->field->has_ext_em)
      gkyl_array_accumulate(fluid_species->pkpm_species->qmem, qbym, app->field->ext_em);

    gkyl_dg_calc_pkpm_vars_source(fluid_species->calc_pkpm_vars, &app->local, 
      fluid_species->pkpm_species->qmem, fluid_species->pkpm_species->pkpm_moms.marr, fluid, rhs);
  }

  gkyl_array_reduce_range(fluid_species->omegaCfl_ptr, fluid_species->cflrate, GKYL_MAX, app->local);

  double omegaCfl_ho[1];
  if (app->use_gpu)
    gkyl_cu_memcpy(omegaCfl_ho, fluid_species->omegaCfl_ptr, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    omegaCfl_ho[0] = fluid_species->omegaCfl_ptr[0];
  omegaCfl = omegaCfl_ho[0];

  app->stat.fluid_species_rhs_tm += gkyl_time_diff_now_sec(wst);

  return app->cfl/omegaCfl;
}

// Apply periodic BCs on fluid species
void
vm_fluid_species_apply_periodic_bc(gkyl_vlasov_app *app, const struct vm_fluid_species *fluid_species,
  int dir, struct gkyl_array *f)
{
  gkyl_array_copy_to_buffer(fluid_species->bc_buffer->data, f, app->skin_ghost.lower_skin[dir]);
  gkyl_array_copy_from_buffer(f, fluid_species->bc_buffer->data, app->skin_ghost.upper_ghost[dir]);

  gkyl_array_copy_to_buffer(fluid_species->bc_buffer->data, f, app->skin_ghost.upper_skin[dir]);
  gkyl_array_copy_from_buffer(f, fluid_species->bc_buffer->data, app->skin_ghost.lower_ghost[dir]);
}

// Determine which directions are periodic and which directions are copy,
// and then apply boundary conditions for fluid species
void
vm_fluid_species_apply_bc(gkyl_vlasov_app *app, const struct vm_fluid_species *fluid_species, struct gkyl_array *f)
{
  int num_periodic_dir = app->num_periodic_dir, cdim = app->cdim;
  int is_np_bc[3] = {1, 1, 1}; // flags to indicate if direction is periodic
  for (int d=0; d<num_periodic_dir; ++d) {
    vm_fluid_species_apply_periodic_bc(app, fluid_species, app->periodic_dirs[d], f);
    is_np_bc[app->periodic_dirs[d]] = 0;
  }
  for (int d=0; d<cdim; ++d) {
    if (is_np_bc[d]) {

      switch (fluid_species->lower_bc[d]) {
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_ABSORB:
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_NO_SLIP:
          gkyl_bc_basic_advance(fluid_species->bc_lo[d], fluid_species->bc_buffer, f);
          break;
        case GKYL_SPECIES_WEDGE:
        case GKYL_SPECIES_FIXED_FUNC:
          assert(false);
          break;
        default:
          break;
      }

      switch (fluid_species->upper_bc[d]) {
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_ABSORB:
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_NO_SLIP:
          gkyl_bc_basic_advance(fluid_species->bc_up[d], fluid_species->bc_buffer, f);
          break;
        case GKYL_SPECIES_WEDGE:
        case GKYL_SPECIES_FIXED_FUNC:
          assert(false);
          break;
        default:
          break;
      }
    }
  }
}

void
vm_fluid_species_calc_int_diag(gkyl_vlasov_app *app, double tm, const struct vm_fluid_species *fluid_species)
{
  if (fluid_species->eqn_id == GKYL_EQN_EULER_PKPM) {
    gkyl_array_clear(fluid_species->integ_mom, 0.0);
    gkyl_dg_calc_pkpm_integrated_vars(fluid_species->calc_pkpm_vars, &app->local, 
      fluid_species->pkpm_species->pkpm_moms.marr, fluid_species->fluid, 
      fluid_species->prim, fluid_species->integ_mom);
    gkyl_array_scale_range(fluid_species->integ_mom, app->grid.cellVolume, app->local);
    
    double int_mom[6] = { 0.0 };
    if (app->use_gpu) {
      gkyl_array_reduce_range(fluid_species->red_integ_diag, fluid_species->integ_mom, GKYL_SUM, app->local);
      gkyl_cu_memcpy(int_mom, fluid_species->red_integ_diag, sizeof(double[6]), GKYL_CU_MEMCPY_D2H);
    }
    else { 
      gkyl_array_reduce_range(int_mom, fluid_species->integ_mom, GKYL_SUM, app->local);
    }
    
    gkyl_dynvec_append(fluid_species->integ_diag, tm, int_mom);
  }
}

void
vm_fluid_species_io(gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species)
{
  gkyl_array_clear(fluid_species->fluid_io, 0.0);
  gkyl_array_clear(fluid_species->pkpm_vars_io, 0.0); 
  if (fluid_species->eqn_id == GKYL_EQN_EULER_PKPM) {
    vm_fluid_species_prim_vars(app, fluid_species, fluid_species->fluid);
    gkyl_dg_calc_pkpm_vars_io(fluid_species->calc_pkpm_vars, &app->local, 
      fluid_species->pkpm_species->pkpm_moms.marr, fluid_species->fluid, 
      fluid_species->p, fluid_species->prim, fluid_species->pkpm_accel, 
      fluid_species->fluid_io, fluid_species->pkpm_vars_io);
  }
}

// release resources for fluid species
void
vm_fluid_species_release(const gkyl_vlasov_app* app, struct vm_fluid_species *f)
{
  gkyl_array_release(f->fluid);
  gkyl_array_release(f->fluid1);
  gkyl_array_release(f->fluidnew);
  gkyl_array_release(f->bc_buffer);
  gkyl_array_release(f->cflrate);

  gkyl_dg_updater_fluid_release(f->advect_slvr);
  if (f->has_diffusion) {
    if (f->diffusion_id == GKYL_GEN_DIFFUSION) {
      gkyl_array_release(f->Dij);
      if (app->use_gpu)
        gkyl_array_release(f->Dij_host);
    }
    gkyl_dg_updater_diffusion_release(f->diff_slvr);
  }

  if (f->eqn_id == GKYL_EQN_ADVECTION) {
    gkyl_array_release(f->app_advect);
    if (app->use_gpu)
      gkyl_array_release(f->app_advect_host);
  }
  else {
    gkyl_array_release(f->prim);
    gkyl_array_release(f->p);
    gkyl_array_release(f->cell_avg_prim);
    if (f->eqn_id == GKYL_EQN_EULER_PKPM) {
      gkyl_array_release(f->fluid_io);
      gkyl_array_release(f->pkpm_vars_io);
      gkyl_array_release(f->pkpm_accel);
      gkyl_dg_calc_pkpm_vars_release(f->calc_pkpm_vars);
      gkyl_dg_calc_pkpm_vars_release(f->calc_pkpm_vars_ext);
      if (app->use_gpu) {
        gkyl_array_release(f->fluid_io_host);
        gkyl_array_release(f->pkpm_vars_io_host);        
      }
    }
  }

  gkyl_array_release(f->integ_mom);
  gkyl_dynvec_release(f->integ_diag);

  if (f->source_id) {
    vm_fluid_species_source_release(app, &f->src);
  }

  if (app->use_gpu) {
    gkyl_array_release(f->fluid_host);
    gkyl_cu_free(f->omegaCfl_ptr);
    gkyl_cu_free(f->red_integ_diag);
  }
  else {
    gkyl_free(f->omegaCfl_ptr);
  }
  // Copy BCs are allocated by default. Need to free.
  for (int d=0; d<app->cdim; ++d) {
    gkyl_bc_basic_release(f->bc_lo[d]);
    gkyl_bc_basic_release(f->bc_up[d]);
  }
}
