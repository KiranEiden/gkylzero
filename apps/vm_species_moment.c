#include <assert.h>
#include <gkyl_vlasov_priv.h>

// list of valid moment names
static const char *const valid_moment_names[] = { "M0", "M1i", "M2ij", "M2", "M3i" };

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
  
  if (app->use_gpu) {
    struct gkyl_mom_type *mtype_host = gkyl_vlasov_mom_new(&app->confBasis, &app->basis, nm);
    sm->mtype = gkyl_vlasov_mom_cu_dev_new(&app->confBasis, &app->basis, nm);
    sm->mcalc = gkyl_mom_calc_cu_dev_new(&s->grid, sm->mtype);

    sm->marr = mkarr(app->use_gpu, mtype_host->num_mom*app->confBasis.num_basis,
      app->local_ext.volume);

    sm->marr_host = mkarr(false, mtype_host->num_mom*app->confBasis.num_basis,
      app->local_ext.volume);

    gkyl_mom_type_release(mtype_host);
  }
  else {
    sm->mtype = gkyl_vlasov_mom_new(&app->confBasis, &app->basis, nm);
    sm->mcalc = gkyl_mom_calc_new(&s->grid, sm->mtype);

    sm->marr = mkarr(app->use_gpu, sm->mtype->num_mom*app->confBasis.num_basis,
      app->local_ext.volume);

    sm->marr_host = sm->marr;
  }
}

// release memory for moment data object
void
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