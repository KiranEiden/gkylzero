#include <assert.h>
#include <gkyl_vlasov_priv.h>

// list of valid moment names
static const char *const valid_moment_names[] = {
  "M0",
  "M1i",
  "M2ij",
  "M2",
  "M3i",
  "M3ijk",
  "FiveMoments",
  "Integrated", // this is an internal flag, not for passing to moment type
  "PKPM", // internal flag for pkpm model which doesn't take a moment name
};

// check if name of moment is valid or not
static bool
is_moment_name_valid(const char *nm)
{
  int n = sizeof(valid_moment_names)/sizeof(valid_moment_names[0]);
  for (int i=0; i<n; ++i)
    if (strcmp(valid_moment_names[i], nm) == 0)
      return 1;
  return 0;
}

// initialize species moment object
void
vm_species_moment_init(struct gkyl_vlasov_app *app, struct vm_species *s,
  struct vm_species_moment *sm, const char *nm)
{
  assert(is_moment_name_valid(nm));

  // For simplicity, is_integrated flag also used by PKPM to turn on diagnostics
  bool is_integrated = strcmp(nm, "Integrated") == 0;
  sm->use_gpu = app->use_gpu;
  sm->mcalc = gkyl_dg_updater_moment_new(&s->grid, &app->confBasis, 
    &app->basis, &app->local, &s->local_vel, s->model_id, nm, is_integrated, s->info.mass, sm->use_gpu);
  int num_mom = gkyl_dg_updater_moment_num_mom(sm->mcalc);

  sm->p_over_gamma = s->p_over_gamma;
  sm->gamma = s->gamma;
  sm->gamma_inv = s->gamma_inv;
  sm->V_drift = s->V_drift;
  sm->GammaV2 = s->GammaV2;
  sm->GammaV_inv = s->GammaV_inv;

  // Using the is_integrated flag as a temporary bool for handling PKPM diagnostics so check that we are *not* PKPM
  if (is_integrated && s->model_id != GKYL_MODEL_PKPM) {
    sm->marr = mkarr(sm->use_gpu, num_mom, app->local_ext.volume);
    sm->marr_host = sm->marr;
    if (sm->use_gpu)
      sm->marr_host = mkarr(false, num_mom, app->local_ext.volume);      
  }
  else {
    sm->marr = mkarr(sm->use_gpu, num_mom*app->confBasis.num_basis,
      app->local_ext.volume);
    sm->marr_host = sm->marr;
    if (sm->use_gpu)
      sm->marr_host = mkarr(false, num_mom*app->confBasis.num_basis,
        app->local_ext.volume);
  }
}

void
vm_species_moment_calc(const struct vm_species_moment *sm,
  const struct gkyl_range phase_rng, const struct gkyl_range conf_rng,
  const struct gkyl_array *fin)
{
  struct gkyl_mom_vlasov_sr_auxfields sr_inp = {.p_over_gamma = sm->p_over_gamma, 
    .gamma = sm->gamma, .gamma_inv = sm->gamma_inv, .V_drift = sm->V_drift, 
    .GammaV2 = sm->GammaV2, .GammaV_inv = sm->GammaV_inv};
  gkyl_dg_updater_moment_advance(sm->mcalc, &phase_rng, &conf_rng, 
    &sr_inp, fin, sm->marr);
}

// release memory for moment data object
void
vm_species_moment_release(const struct gkyl_vlasov_app *app, const struct vm_species_moment *sm)
{
  if (app->use_gpu)
    gkyl_array_release(sm->marr_host);

  gkyl_dg_updater_moment_release(sm->mcalc);
  gkyl_array_release(sm->marr);
}
