#include <assert.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_diffusion_vlasov.h>
#include <gkyl_dg_updater_diffusion_vlasov.h>
#include <gkyl_dg_updater_diffusion_vlasov_priv.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_util.h>

struct gkyl_dg_eqn*
gkyl_dg_updater_diffusion_vlasov_acquire_eqn(const struct gkyl_dg_updater_diffusion_vlasov *up)
{
  return gkyl_dg_eqn_acquire(up->dgeqn);
}

struct gkyl_dg_updater_diffusion_vlasov*
gkyl_dg_updater_diffusion_vlasov_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *basis, const struct gkyl_basis *cbasis, bool is_diff_const,
  const bool *diff_in_dir, int diff_order, const struct gkyl_range *diff_range,
  const bool *is_zero_flux_dir, bool use_gpu)
{
  struct gkyl_dg_updater_diffusion_vlasov *up = gkyl_malloc(sizeof(struct gkyl_dg_updater_diffusion_vlasov));

  int cdim = cbasis->ndim;
  up->use_gpu = use_gpu;
  bool is_dir_diffusive[GKYL_MAX_CDIM];
  for (int d=0; d<cdim; d++) is_dir_diffusive[d] = diff_in_dir==NULL? true : diff_in_dir[d];

  up->dgeqn = gkyl_dg_diffusion_vlasov_new(basis, cbasis, is_diff_const, is_dir_diffusive,
                                           diff_order, diff_range, up->use_gpu);

  int num_up_dirs = 0;
  for (int d=0; d<cdim; d++) num_up_dirs += is_dir_diffusive[d]? 1 : 0;
  int up_dirs[GKYL_MAX_DIM], zero_flux_flags[GKYL_MAX_DIM];
  int linc = 0;
  for (int d=0; d<cdim; ++d) {
    if (is_dir_diffusive[d]) up_dirs[linc] = d;
    linc += 1;

    zero_flux_flags[d] = is_zero_flux_dir[d]? 1 : 0;
  }

  up->hyperdg = gkyl_hyper_dg_new(grid, basis, up->dgeqn, num_up_dirs, up_dirs, zero_flux_flags, 1, up->use_gpu);

  up->diffusion_tm = 0.0;

  return up;
}

void
gkyl_dg_updater_diffusion_vlasov_advance(struct gkyl_dg_updater_diffusion_vlasov *up,
  const struct gkyl_range *update_rng, const struct gkyl_array *coeff,
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT cflrate,
  struct gkyl_array* GKYL_RESTRICT rhs)
{
  struct timespec wst = gkyl_wall_clock();
  // Set arrays needed and call the specific advance method required
  gkyl_dg_diffusion_vlasov_set_auxfields(up->dgeqn, (struct gkyl_dg_diffusion_vlasov_auxfields) { .D = coeff });
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    gkyl_hyper_dg_advance_cu(up->hyperdg, update_rng, fIn, cflrate, rhs);
  else
    gkyl_hyper_dg_advance(up->hyperdg, update_rng, fIn, cflrate, rhs);
#else
  gkyl_hyper_dg_advance(up->hyperdg, update_rng, fIn, cflrate, rhs);
#endif
  up->diffusion_tm += gkyl_time_diff_now_sec(wst);
}

struct gkyl_dg_updater_diffusion_vlasov_tm
gkyl_dg_updater_diffusion_vlasov_get_tm(const struct gkyl_dg_updater_diffusion_vlasov *up)
{
  return (struct gkyl_dg_updater_diffusion_vlasov_tm) {
    .diffusion_tm = up->diffusion_tm,
  };
}

void
gkyl_dg_updater_diffusion_vlasov_release(struct gkyl_dg_updater_diffusion_vlasov *up)
{
  gkyl_dg_eqn_release(up->dgeqn);
  gkyl_hyper_dg_release(up->hyperdg);
  gkyl_free(up);
}
