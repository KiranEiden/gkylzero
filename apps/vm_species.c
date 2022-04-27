#include "gkyl_alloc.h"
#include "gkyl_array_ops.h"
#include <assert.h>

#include <gkyl_app.h>
#include <gkyl_array.h>
#include <gkyl_eqn_type.h>
#include <gkyl_vlasov_priv.h>
#include <gkyl_proj_on_basis.h>

// Projection functions for p/(m*gamma) = v in special relativistic systems
// Simplifies to p/sqrt(m^2 + p^2) where c = 1
static void 
ev_p_over_gamma_1p(double t, const double *xn, double *out, void *ctx)
{
  struct gamma_ctx *gtx = (struct gamma_ctx*) ctx;
  double mass = gtx->mass;
  out[0] = xn[0]/sqrt(mass*mass + xn[0]*xn[0]);
}
static void 
ev_p_over_gamma_2p(double t, const double *xn, double *out, void *ctx)
{
  struct gamma_ctx *gtx = (struct gamma_ctx*) ctx;
  double mass = gtx->mass;
  out[0] = xn[0]/sqrt(mass*mass + xn[0]*xn[0] + xn[1]*xn[1]);
  out[1] = xn[1]/sqrt(mass*mass + xn[0]*xn[0] + xn[1]*xn[1]);
}
static void 
ev_p_over_gamma_3p(double t, const double *xn, double *out, void *ctx)
{
  struct gamma_ctx *gtx = (struct gamma_ctx*) ctx;
  double mass = gtx->mass;
  out[0] = xn[0]/sqrt(mass*mass + xn[0]*xn[0] + xn[1]*xn[1] + xn[2]*xn[2]);
  out[1] = xn[1]/sqrt(mass*mass + xn[0]*xn[0] + xn[1]*xn[1] + xn[2]*xn[2]);
  out[2] = xn[2]/sqrt(mass*mass + xn[0]*xn[0] + xn[1]*xn[1] + xn[2]*xn[2]);
}

static const evalf_t p_over_gamma_func[3] = {ev_p_over_gamma_1p, ev_p_over_gamma_2p, ev_p_over_gamma_3p};

// function to evaluate acceleration (this is needed as accel function
// provided by the user returns 3 components, while the Vlasov solver
// expects 8 components to match the EM field)
static void
eval_accel(double t, const double *xn, double *aout, void *ctx)
{
  struct vm_eval_accel_ctx *a_ctx = ctx;
  double a[3]; // output acceleration
  a_ctx->accel_func(t, xn, a, a_ctx->accel_ctx);
  
  for (int i=0; i<3; ++i) aout[i] = a[i];
  for (int i=3; i<8; ++i) aout[i] = 0.0;
}

// initialize species object
void
vm_species_init(struct gkyl_vm *vm, struct gkyl_vlasov_app *app, struct vm_species *s)
{
  int cdim = app->cdim, vdim = app->vdim;
  int pdim = cdim+vdim;

  int cells[GKYL_MAX_DIM], ghost[GKYL_MAX_DIM];
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];

  int cells_vel[GKYL_MAX_DIM], ghost_vel[GKYL_MAX_DIM];
  double lower_vel[GKYL_MAX_DIM], upper_vel[GKYL_MAX_DIM];

  for (int d=0; d<cdim; ++d) {
    cells[d] = vm->cells[d];
    lower[d] = vm->lower[d];
    upper[d] = vm->upper[d];
    ghost[d] = 1;
  }
  for (int d=0; d<vdim; ++d) {
    // full phase space grid
    cells[cdim+d] = s->info.cells[d];
    lower[cdim+d] = s->info.lower[d];
    upper[cdim+d] = s->info.upper[d];
    ghost[cdim+d] = 0; // no ghost-cells in velocity space

    // only velocity space
    cells_vel[d] = s->info.cells[d];
    lower_vel[d] = s->info.lower[d];
    upper_vel[d] = s->info.upper[d];
    ghost_vel[d] = 0; // no ghost-cells in velocity space
  }
  // full phase space grid
  gkyl_rect_grid_init(&s->grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&s->grid, ghost, &s->local_ext, &s->local);
  
  skin_ghost_ranges_init(&s->skin_ghost, &s->local_ext, ghost);
  
  // velocity space grid
  gkyl_rect_grid_init(&s->grid_vel, vdim, lower_vel, upper_vel, cells_vel);
  gkyl_create_grid_ranges(&s->grid_vel, ghost_vel, &s->local_ext_vel, &s->local_vel);

  // allocate distribution function arrays
  s->f = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
  s->f1 = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
  s->fnew = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);

  s->f_host = s->f;
  if (app->use_gpu)
    s->f_host = mkarr(false, app->basis.num_basis, s->local_ext.volume);

  // allocate cflrate (scalar array)
  s->cflrate = mkarr(app->use_gpu, 1, s->local_ext.volume);
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
  s->bc_buffer = mkarr(app->use_gpu, app->basis.num_basis, buff_sz);
  
  // determine field-type 
  s->field_id = app->has_field ? app->field->info.field_id : GKYL_FIELD_NULL;

  // allocate array to store q/m*(E,B) or potential depending on equation system
  s->qmem = 0;
  s->fac_phi = 0;
  if (s->field_id  == GKYL_FIELD_E_B || s->field_id  == GKYL_FIELD_SR_E_B)
    s->qmem = mkarr(app->use_gpu, 8*app->confBasis.num_basis, app->local_ext.volume);
  else if (s->field_id == GKYL_FIELD_PHI || s->field_id == GKYL_FIELD_PHI_A)
    s->fac_phi = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);

  // allocate array to store vector potential if present
  s->vecA = 0;
  if (s->field_id == GKYL_FIELD_PHI_A)
    s->vecA = mkarr(app->use_gpu, 3*app->confBasis.num_basis, app->local_ext.volume);

  // allocate array to store p/gamma (velocity) if present
  // Since p/gamma is a geometric quantity, can pre-compute it here
  s->p_over_gamma = 0;
  if (s->field_id  == GKYL_FIELD_SR_E_B) {
    s->p_over_gamma = mkarr(app->use_gpu, vdim*app->velBasis.num_basis, s->local_vel.volume);
    s->p_over_gamma_host = s->p_over_gamma;
    if (app->use_gpu)
      s->p_over_gamma_host = mkarr(false, vdim*app->velBasis.num_basis, s->local_vel.volume);
    struct gamma_ctx ctx = { .mass = s->info.mass };
    gkyl_proj_on_basis *p_over_gamma_proj = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
        .grid = &s->grid_vel,
        .basis = &app->velBasis,
        .qtype = GKYL_GAUSS_LOBATTO_QUAD,
        .num_quad = 8,
        .num_ret_vals = vdim,
        .eval = p_over_gamma_func[vdim-1],
        .ctx = &ctx
      }
    );  
    // run updater
    gkyl_proj_on_basis_advance(p_over_gamma_proj, 0.0, &s->local_vel, s->p_over_gamma_host);
    gkyl_proj_on_basis_release(p_over_gamma_proj);    
    if (app->use_gpu) // note: p_over_gamma_host is same as p_over_gamma when not on GPUs
      gkyl_array_copy(s->p_over_gamma, s->p_over_gamma_host);
  }

  // create solver
  s->slvr = gkyl_dg_updater_vlasov_new(&s->grid, &app->confBasis, &app->basis, &app->local, &s->local_vel, s->field_id, app->use_gpu);

  // acquire equation object
  s->eqn_vlasov = gkyl_dg_updater_vlasov_acquire_eqn(s->slvr);

  // allocate data for momentum (for use in current accumulation)
  vm_species_moment_init(app, s, &s->m1i, "M1i");
  
  int ndm = s->info.num_diag_moments;
  // allocate data for diagnostic moments
  s->moms = gkyl_malloc(sizeof(struct vm_species_moment[ndm]));
  for (int m=0; m<ndm; ++m)
    vm_species_moment_init(app, s, &s->moms[m], s->info.diag_moments[m]);

  s->has_accel = false;
  // setup applied acceleration
  if (s->info.accel) {
    s->has_accel = true;
    // we need to ensure applied acceleration has same shape as EM
    // field as it will get added to qmem
    s->accel = mkarr(app->use_gpu, 8*app->confBasis.num_basis, app->local_ext.volume);

    s->accel_host = s->accel;
    if (app->use_gpu)
      s->accel_host = mkarr(false, 8*app->confBasis.num_basis, app->local_ext.volume);

    s->accel_ctx = (struct vm_eval_accel_ctx) {
      .accel_func = s->info.accel, .accel_ctx = s->info.accel_ctx
    };
    s->accel_proj = gkyl_proj_on_basis_new(&app->grid, &app->confBasis, app->confBasis.poly_order+1,
      8, eval_accel, &s->accel_ctx);
  }

  // determine collision type to use in vlasov update
  s->collision_id = s->info.collisions.collision_id;
  s->collides_with_fluid = false;
  if (s->collision_id == GKYL_LBO_COLLISIONS) {
    if (vm_find_fluid_species(app, s->info.collisions.collide_with_fluid)) {
      s->collides_with_fluid = true;
      // index in fluid_species struct of fluid species kinetic species is colliding with
      s->fluid_index = s->info.collisions.fluid_index;
    }
    vm_species_lbo_init(app, s, &s->lbo, s->collides_with_fluid);
  }

  // determine which directions are not periodic
  int num_periodic_dir = app->num_periodic_dir, is_np[3] = {1, 1, 1};
  for (int d=0; d<num_periodic_dir; ++d)
    is_np[app->periodic_dirs[d]] = 0;

  for (int dir=0; dir<app->cdim; ++dir) {
    if (is_np[dir]) {
      const enum gkyl_species_bc_type *bc;
      if (dir == 0)
        bc = s->info.bcx;
      else if (dir == 1)
        bc = s->info.bcy;
      else
        bc = s->info.bcz;

      s->lower_bc[dir] = bc[0];
      s->upper_bc[dir] = bc[1];
    }
  }
  for (int d=0; d<3; ++d) {
    if (s->field_id == GKYL_FIELD_SR_E_B)
      s->wall_bc_func[d] = gkyl_vlasov_sr_wall_bc_create(s->eqn_vlasov, d,
        app->basis_on_dev.basis);
    else if (s->field_id == GKYL_FIELD_PHI || s->field_id == GKYL_FIELD_PHI_A)
      s->wall_bc_func[d] = gkyl_vlasov_poisson_wall_bc_create(s->eqn_vlasov, d,
        app->basis_on_dev.basis);
    else
      s->wall_bc_func[d] = gkyl_vlasov_wall_bc_create(s->eqn_vlasov, d,
        app->basis_on_dev.basis);
  }
}

void
vm_species_apply_ic(gkyl_vlasov_app *app, struct vm_species *species, double t0)
{
  int poly_order = app->poly_order;
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&species->grid, &app->basis,
    poly_order+1, 1, species->info.init, species->info.ctx);

  // run updater
  gkyl_proj_on_basis_advance(proj, t0, &species->local, species->f_host);
  gkyl_proj_on_basis_release(proj);    

  if (app->use_gpu) // note: f_host is same as f when not on GPUs
    gkyl_array_copy(species->f, species->f_host);

  // we are computing acceleration for now as it is time-independent
  vm_species_calc_accel(app, species, t0);
}

void
vm_species_calc_accel(gkyl_vlasov_app *app, struct vm_species *species, double tm)
{
  if (species->has_accel) {
    gkyl_proj_on_basis_advance(species->accel_proj, tm, &app->local_ext, species->accel_host);
    if (app->use_gpu) // note: accel_host is same as accel when not on GPUs
      gkyl_array_copy(species->accel, species->accel_host);
  }
}

// Compute the RHS for species update, returning maximum stable
// time-step.
double
vm_species_rhs(gkyl_vlasov_app *app, struct vm_species *species,
  const struct gkyl_array *fin, const struct gkyl_array *em, struct gkyl_array *rhs)
{
  gkyl_array_clear(species->cflrate, 0.0);
  if (species->field_id  == GKYL_FIELD_E_B || species->field_id  == GKYL_FIELD_SR_E_B) {
    double qbym = species->info.charge/species->info.mass;
    gkyl_array_set(species->qmem, qbym, em);

    if (species->has_accel)
      gkyl_array_accumulate(species->qmem, 1.0, species->accel);
  }
  
  gkyl_array_clear(rhs, 0.0);
 
  if (app->use_gpu)
    gkyl_dg_updater_vlasov_advance_cu(species->slvr, species->field_id, &species->local, 
      species->qmem, species->p_over_gamma, 
      fin, species->cflrate, rhs);
  else
    gkyl_dg_updater_vlasov_advance(species->slvr, species->field_id, &species->local, 
      species->qmem, species->p_over_gamma, 
      fin, species->cflrate, rhs);

  if (species->collision_id == GKYL_LBO_COLLISIONS)
    vm_species_lbo_rhs(app, species, &species->lbo, fin, rhs);

  gkyl_array_reduce_range(species->omegaCfl_ptr, species->cflrate, GKYL_MAX, species->local);
  double omegaCfl = species->omegaCfl_ptr[0];
  
  return app->cfl/omegaCfl;
}

// Apply periodic BCs on distribution function
void
vm_species_apply_periodic_bc(gkyl_vlasov_app *app, const struct vm_species *species,
  int dir, struct gkyl_array *f)
{
  gkyl_array_copy_to_buffer(species->bc_buffer->data, f, species->skin_ghost.lower_skin[dir]);
  gkyl_array_copy_from_buffer(f, species->bc_buffer->data, species->skin_ghost.upper_ghost[dir]);

  gkyl_array_copy_to_buffer(species->bc_buffer->data, f, species->skin_ghost.upper_skin[dir]);
  gkyl_array_copy_from_buffer(f, species->bc_buffer->data, species->skin_ghost.lower_ghost[dir]);
}

// Apply copy BCs on distribution function
void
vm_species_apply_copy_bc(gkyl_vlasov_app *app, const struct vm_species *species,
  int dir, enum vm_domain_edge edge, struct gkyl_array *f)
{
  if (edge == VM_EDGE_LOWER) {
    gkyl_array_copy_to_buffer(species->bc_buffer->data, f, species->skin_ghost.lower_skin[dir]);
    gkyl_array_copy_from_buffer(f, species->bc_buffer->data, species->skin_ghost.lower_ghost[dir]);
  }

  if (edge == VM_EDGE_UPPER) {
    gkyl_array_copy_to_buffer(species->bc_buffer->data, f, species->skin_ghost.upper_skin[dir]);
    gkyl_array_copy_from_buffer(f, species->bc_buffer->data, species->skin_ghost.upper_ghost[dir]);
  }
}

// Apply wall BCs on distribution function
void
vm_species_apply_wall_bc(gkyl_vlasov_app *app, const struct vm_species *species,
  int dir, enum vm_domain_edge edge, struct gkyl_array *f)
{
  int cdim = app->cdim;
  
  if (edge == VM_EDGE_LOWER) {
    gkyl_array_flip_copy_to_buffer_fn(species->bc_buffer->data, f, dir+cdim,
      species->skin_ghost.lower_skin[dir],
      species->wall_bc_func[dir]
    );
    
    gkyl_array_copy_from_buffer(f, species->bc_buffer->data, species->skin_ghost.lower_ghost[dir]);
  }

  if (edge == VM_EDGE_UPPER) {
    gkyl_array_flip_copy_to_buffer_fn(species->bc_buffer->data, f, dir+cdim,
      species->skin_ghost.upper_skin[dir],
      species->wall_bc_func[dir]
    );
    
    gkyl_array_copy_from_buffer(f, species->bc_buffer->data, species->skin_ghost.upper_ghost[dir]);
  }
}

// Determine which directions are periodic and which directions are copy,
// and then apply boundary conditions for distribution function
void
vm_species_apply_bc(gkyl_vlasov_app *app, const struct vm_species *species, struct gkyl_array *f)
{
  int num_periodic_dir = app->num_periodic_dir, cdim = app->cdim;
  int is_np_bc[3] = {1, 1, 1}; // flags to indicate if direction is periodic
  for (int d=0; d<num_periodic_dir; ++d) {
    vm_species_apply_periodic_bc(app, species, app->periodic_dirs[d], f);
    is_np_bc[app->periodic_dirs[d]] = 0;
  }
  for (int d=0; d<cdim; ++d) {
    if (is_np_bc[d]) {

      switch (species->lower_bc[d]) {
        case GKYL_SPECIES_COPY:
          vm_species_apply_copy_bc(app, species, d, VM_EDGE_LOWER, f);
          break;
        case GKYL_SPECIES_WALL:
          vm_species_apply_wall_bc(app, species, d, VM_EDGE_LOWER, f);
          break;
      }

      switch (species->upper_bc[d]) {
        case GKYL_SPECIES_COPY:
          vm_species_apply_copy_bc(app, species, d, VM_EDGE_UPPER, f);
          break;
        case GKYL_SPECIES_WALL:
          vm_species_apply_wall_bc(app, species, d, VM_EDGE_UPPER, f);
          break;
      }      
    }
  }
}

void
vm_species_coll_tm(gkyl_vlasov_app *app)
{
  for (int i=0; i<app->num_species; ++i) {
    if (app->species[i].collision_id == GKYL_LBO_COLLISIONS) {
      struct gkyl_dg_updater_lbo_vlasov_tm tm =
        gkyl_dg_updater_lbo_vlasov_get_tm(app->species[i].lbo.coll_slvr);
      app->stat.species_lbo_coll_diff_tm[i] = tm.diff_tm;
      app->stat.species_lbo_coll_drag_tm[i] = tm.drag_tm;
    }
  }
}

void
vm_species_tm(gkyl_vlasov_app *app)
{
  for (int i=0; i<app->num_species; ++i) {
    struct gkyl_dg_updater_vlasov_tm tm =
      gkyl_dg_updater_vlasov_get_tm(app->species[i].slvr);
    app->stat.species_rhs_tm += tm.vlasov_tm;
  }
}

// release resources for species
void
vm_species_release(const gkyl_vlasov_app* app, const struct vm_species *s)
{
  // release various arrays
  gkyl_array_release(s->f);
  gkyl_array_release(s->f1);
  gkyl_array_release(s->fnew);
  gkyl_array_release(s->cflrate);
  gkyl_array_release(s->bc_buffer);

  if (app->use_gpu)
    gkyl_array_release(s->f_host);
  
  // release moment data
  vm_species_moment_release(app, &s->m1i);
  for (int i=0; i<s->info.num_diag_moments; ++i)
    vm_species_moment_release(app, &s->moms[i]);
  gkyl_free(s->moms);

  gkyl_dg_eqn_release(s->eqn_vlasov);
  gkyl_dg_updater_vlasov_release(s->slvr);

  // Release arrays for different types of Vlasov equations
  if (s->field_id  == GKYL_FIELD_E_B || s->field_id  == GKYL_FIELD_SR_E_B)
    gkyl_array_release(s->qmem);
  else if (s->field_id == GKYL_FIELD_PHI || s->field_id == GKYL_FIELD_PHI_A)
    gkyl_array_release(s->fac_phi);

  if (s->field_id == GKYL_FIELD_PHI_A)
    gkyl_array_release(s->vecA);

  if (s->field_id  == GKYL_FIELD_SR_E_B) {
    gkyl_array_release(s->p_over_gamma);
    if (app->use_gpu)
      gkyl_array_release(s->p_over_gamma_host);
  }
  
  if (s->has_accel) {
    gkyl_array_release(s->accel);
    if (app->use_gpu)
      gkyl_array_release(s->accel_host);

    gkyl_proj_on_basis_release(s->accel_proj);
  }

  if (s->collision_id == GKYL_LBO_COLLISIONS)
    vm_species_lbo_release(app, &s->lbo);

  for (int d=0; d<3; ++d)
    gkyl_vlasov_wall_bc_release(s->wall_bc_func[d]);
  
  if (app->use_gpu)
    gkyl_cu_free_host(s->omegaCfl_ptr);
  else
    gkyl_free(s->omegaCfl_ptr);
}
