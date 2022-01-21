#include "gkyl_dg_eqn.h"
#include <assert.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_dg_lbo_updater.h>
#include <gkyl_dg_lbo_updater_priv.h>
#include <gkyl_dg_vlasov_lbo_drag.h>
#include <gkyl_dg_vlasov_lbo_diff.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_util.h>

gkyl_dg_lbo_updater*
gkyl_dg_lbo_updater_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis,
  const struct gkyl_basis *pbasis, const struct gkyl_range *conf_range)
{
  gkyl_dg_lbo_updater *up = gkyl_malloc(sizeof(gkyl_dg_lbo_updater));

  up->coll_drag = gkyl_dg_vlasov_lbo_drag_new(cbasis, pbasis, conf_range);
  up->coll_diff = gkyl_dg_vlasov_lbo_diff_new(cbasis, pbasis, conf_range);

  int cdim = cbasis->ndim, pdim = pbasis->ndim;
  int vdim = pdim-cdim;
  int num_up_dirs = vdim;
  int up_dirs[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<vdim; ++d)
    up_dirs[d] = d + pbasis->ndim - vdim;

  int zero_flux_flags[GKYL_MAX_DIM] = { 0 };
  for (int d=cdim; d<pdim; ++d)
    zero_flux_flags[d] = 1;
  
  up->diff = gkyl_hyper_dg_new(grid, pbasis, up->coll_diff, num_up_dirs, up_dirs, zero_flux_flags, 1);
  up->drag = gkyl_hyper_dg_new(grid, pbasis, up->coll_drag, num_up_dirs, up_dirs, zero_flux_flags, 1);
  
  return up;
}

void
gkyl_dg_lbo_updater_advance(gkyl_dg_lbo_updater *lbo, struct gkyl_range update_rng,
  const struct gkyl_array *nu_sum, const struct gkyl_array *nu_u, const struct gkyl_array *nu_vthsq,
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs)
{
  // Set arrays needed
  gkyl_vlasov_lbo_drag_set_nuSum(lbo->coll_drag, nu_sum);
  gkyl_vlasov_lbo_drag_set_nuUSum(lbo->coll_drag, nu_u);
  gkyl_vlasov_lbo_drag_set_nuVtSqSum(lbo->coll_drag, nu_vthsq);
  gkyl_vlasov_lbo_diff_set_nuSum(lbo->coll_diff, nu_sum);
  gkyl_vlasov_lbo_diff_set_nuUSum(lbo->coll_diff, nu_u);
  gkyl_vlasov_lbo_diff_set_nuVtSqSum(lbo->coll_diff, nu_vthsq);
  
  gkyl_hyper_dg_advance(lbo->diff, update_rng, fIn, cflrate, rhs);
  gkyl_hyper_dg_advance(lbo->drag, update_rng, fIn, cflrate, rhs);
}

void
gkyl_dg_lbo_updater_release(gkyl_dg_lbo_updater* lbo)
{
  gkyl_dg_eqn_release(lbo->coll_diff);
  gkyl_dg_eqn_release(lbo->coll_drag);
  gkyl_hyper_dg_release(lbo->drag);
  gkyl_hyper_dg_release(lbo->diff);
  gkyl_free(lbo);
}

#ifdef GKYL_HAVE_CUDA

gkyl_dg_lbo_updater*
gkyl_dg_lbo_updater_cu_dev_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis,
  const struct gkyl_basis *pbasis, const struct gkyl_range *conf_range)
{
  gkyl_dg_lbo_updater *up = (gkyl_dg_lbo_updater*) gkyl_malloc(sizeof(gkyl_dg_lbo_updater));

  up->coll_drag = gkyl_dg_vlasov_lbo_drag_cu_dev_new(cbasis, pbasis, conf_range);
  up->coll_diff = gkyl_dg_vlasov_lbo_diff_cu_dev_new(cbasis, pbasis, conf_range);

  int cdim = cbasis->ndim, pdim = pbasis->ndim;
  int vdim = pdim-cdim;
  int num_up_dirs = vdim;
  int up_dirs[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<vdim; ++d)
    up_dirs[d] = d + pbasis->ndim - vdim;

  int zero_flux_flags[GKYL_MAX_DIM] = { 0 };
  for (int d=cdim; d<pdim; ++d)
    zero_flux_flags[d] = 1;
  up->diff = gkyl_hyper_dg_cu_dev_new(grid, pbasis, up->coll_diff, num_up_dirs, up_dirs, zero_flux_flags, 1);
  up->drag = gkyl_hyper_dg_cu_dev_new(grid, pbasis, up->coll_drag, num_up_dirs, up_dirs, zero_flux_flags, 1);

  return up;
}

void
gkyl_dg_lbo_updater_advance_cu(gkyl_dg_lbo_updater *lbo, struct gkyl_range update_rng,
  const struct gkyl_array *nu_sum, const struct gkyl_array *nu_u, const struct gkyl_array *nu_vthsq,
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs)
{
  // Set arrays needed
  gkyl_vlasov_lbo_drag_set_nuSum(lbo->coll_drag, nu_sum);
  gkyl_vlasov_lbo_drag_set_nuUSum(lbo->coll_drag, nu_u);
  gkyl_vlasov_lbo_drag_set_nuVtSqSum(lbo->coll_drag, nu_vthsq);
  gkyl_vlasov_lbo_diff_set_nuSum(lbo->coll_diff, nu_sum);
  gkyl_vlasov_lbo_diff_set_nuUSum(lbo->coll_diff, nu_u);
  gkyl_vlasov_lbo_diff_set_nuVtSqSum(lbo->coll_diff, nu_vthsq);

  gkyl_hyper_dg_advance_cu(lbo->diff, update_rng, fIn, cflrate, rhs);
  gkyl_hyper_dg_advance_cu(lbo->drag, update_rng, fIn, cflrate, rhs);
}

#endif

#ifndef GKYL_HAVE_CUDA

gkyl_dg_lbo_updater*
gkyl_dg_lbo_updater_cu_dev_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis,
  const struct gkyl_basis *pbasis, const struct gkyl_range *conf_range)
{
  assert(false);
}

void
gkyl_dg_lbo_updater_advance_cu(gkyl_dg_lbo_updater *lbo, struct gkyl_range update_rng,
  const struct gkyl_array *nu_sum, const struct gkyl_array *nu_u, const struct gkyl_array *nu_vthsq,                  const struct gkyl_array *fIn, struct gkyl_array *cflrate, struct gkyl_array *rhs)
{
  assert(false);
}

#endif
