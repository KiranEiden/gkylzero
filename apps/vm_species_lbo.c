#include <assert.h>
#include <gkyl_vlasov_priv.h>

void 
vm_species_lbo_init(struct gkyl_vlasov_app *app, struct vm_species *s, struct vm_lbo_collisions *lbo)
{
  // TO DO: Expose nu_u and nu_vthsq arrays above species object
  //        for cross-species collisions. Just testing for now JJ 09/24/21
  int cdim = app->cdim, vdim = app->vdim;
  double v_bounds[2*GKYL_MAX_DIM];
  for (int d=0; d<vdim; ++d) {
    v_bounds[d] = s->info.lower[d];
    v_bounds[d + vdim] = s->info.upper[d];
  }

  // allocate nu and initialize it
  lbo->nu_sum = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
    app->poly_order+1, 1, s->info.nu, s->info.ctx);

  gkyl_proj_on_basis_advance(proj, 0.0, &app->local, lbo->nu_sum);
  gkyl_proj_on_basis_release(proj);

  lbo->cM = mkarr(app->use_gpu, vdim*app->confBasis.num_basis, app->local_ext.volume);
  lbo->cE = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  lbo->u_drift = mkarr(app->use_gpu, vdim*app->confBasis.num_basis, app->local_ext.volume);
  lbo->vth_sq = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  lbo->nu_u = mkarr(app->use_gpu, vdim*app->confBasis.num_basis, app->local_ext.volume);
  lbo->nu_vthsq = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);

  // allocate moments needed for LBO update
  vm_species_moment_init(app, s, &lbo->m0, "M0");
  vm_species_moment_init(app, s, &lbo->m1i, "M1i");
  vm_species_moment_init(app, s, &lbo->m2, "M2");
  
  // create collision equation object and solver
  if (app->use_gpu) {
    //lbo->coll_eqn = gkyl_dg_vlasov_lbo_cu_dev_new(&app->confBasis, &app->basis, &app->local);
    //lbo->coll_slvr = gkyl_hyper_dg_cu_dev_new(&s->grid, &app->basis, lbo->coll_eqn, num_up_dirs, up_dirs, zero_flux_flags, 1);
  }
  else {
    // edge of velocity space momentum correction 
    lbo->cM_mom = gkyl_vlasov_lbo_mom_new(&app->confBasis, &app->basis, "f", v_bounds);
    lbo->cM_bcorr = gkyl_mom_bcorr_new(&s->grid, lbo->cM_mom);

    // Edge of velocity space energy correction
    lbo->cE_mom = gkyl_vlasov_lbo_mom_new(&app->confBasis, &app->basis, "vf",  v_bounds);
    lbo->cE_bcorr = gkyl_mom_bcorr_new(&s->grid, lbo->cE_mom);
    
    // primitive moment calculators
    lbo->coll_prim = gkyl_prim_lbo_vlasov_new(&app->confBasis, &app->basis);
    lbo->coll_pcalc = gkyl_prim_lbo_calc_new(&s->grid, lbo->coll_prim);
    
    // LBO updater
    lbo->coll_slvr = gkyl_dg_lbo_updater_new(&s->grid, &app->confBasis, &app->basis, &app->local);
  }
}

double
vm_species_lbo_rhs(gkyl_vlasov_app *app, const struct vm_species *species,
  struct vm_lbo_collisions *lbo, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  // compute needed moments
  vm_species_moment_calc(&lbo->m0, species->local, app->local, fin);
  vm_species_moment_calc(&lbo->m1i, species->local, app->local, fin);
  vm_species_moment_calc(&lbo->m2, species->local, app->local, fin);

  // construct boundary corrections
  gkyl_mom_bcorr_advance(lbo->cM_bcorr, species->local, app->local, fin, lbo->cM);
  gkyl_mom_bcorr_advance(lbo->cE_bcorr, species->local, app->local, fin, lbo->cE);

  // construct primitive moments
  gkyl_prim_lbo_calc_advance(lbo->coll_pcalc, app->confBasis, app->local, 
    lbo->m0.marr, lbo->m1i.marr, lbo->m2.marr, lbo->cM, lbo->cE,
        lbo->u_drift, lbo->vth_sq);
  
  gkyl_dg_mul_op(app->confBasis, 0, lbo->nu_u, 0, lbo->u_drift, 0, lbo->nu_sum);
  gkyl_dg_mul_op(app->confBasis, 0, lbo->nu_vthsq, 0, lbo->vth_sq, 0, lbo->nu_sum);
  
  // acccumulate update due to collisions onto rhs
  gkyl_dg_lbo_updater_advance(lbo->coll_slvr, species->local,
    lbo->nu_sum, lbo->nu_u, lbo->nu_vthsq, fin, species->cflrate, rhs);

  // TODO: This needs to be set properly!
  return 0;
}

void 
vm_species_lbo_release(const struct gkyl_vlasov_app *app, const struct vm_lbo_collisions *lbo)
{
  gkyl_array_release(lbo->cM);
  gkyl_array_release(lbo->cE);
  gkyl_array_release(lbo->u_drift);
  gkyl_array_release(lbo->vth_sq);
  gkyl_array_release(lbo->nu_sum);
  gkyl_array_release(lbo->nu_u);
  gkyl_array_release(lbo->nu_vthsq);

  vm_species_moment_release(app, &lbo->m0);
  vm_species_moment_release(app, &lbo->m1i);
  vm_species_moment_release(app, &lbo->m2);
  
  if (app->use_gpu) {
  }
  else {
    gkyl_mom_type_release(lbo->cM_mom);
    gkyl_mom_type_release(lbo->cE_mom);
    gkyl_mom_bcorr_release(lbo->cM_bcorr);
    gkyl_mom_bcorr_release(lbo->cE_bcorr);
    gkyl_prim_lbo_calc_release(lbo->coll_pcalc);
    gkyl_dg_lbo_updater_release(lbo->coll_slvr);
  }
}
