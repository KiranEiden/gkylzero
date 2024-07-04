#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_math.h>
#include <gkyl_util.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wave_geom_priv.h>

bool
gkyl_wave_geom_is_cu_dev(const struct gkyl_wave_geom* wg)
{
  return GKYL_IS_CU_ALLOC(wg->flags);
}

void
gkyl_wave_geom_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_wave_geom *wg = container_of(ref, struct gkyl_wave_geom, ref_count);
  gkyl_array_release(wg->geom);
  if (gkyl_wave_geom_is_cu_dev(wg)) 
    gkyl_cu_free(wg->on_dev); 

  gkyl_free(wg);
}

struct gkyl_wave_geom*
gkyl_wave_geom_new(const struct gkyl_rect_grid *grid, struct gkyl_range *range,
  evalf_t mapc2p, void *ctx, bool use_gpu)
{
  struct gkyl_wave_geom_inp inp =
  {
    .mapc2p = mapc2p ? mapc2p : nomapc2p,
    .get_cov_basis = get_standard_basis,
    .get_con_basis = get_standard_basis,
    .spacetime = 0
  };
  struct gkyl_wave_geom *wg = gkyl_wave_geom_from_wg_inp(grid, range, &inp, ctx, use_gpu);
  return wg;
}

void
gkyl_wave_geom_inp_from_flag(enum gkyl_wave_coord_flag cflag, struct gkyl_wave_geom_inp *inp)
{
  switch(cflag)
  {
    case WAVE_COORD_CART:
      inp->mapc2p = nomapc2p;
      inp->get_cov_basis = get_standard_basis;
      inp->get_con_basis = get_standard_basis;
      inp->spacetime = 0;
      break;
      
    case WAVE_COORD_CYL:
      inp->mapc2p = mapc2p_cyl;
      inp->get_cov_basis = get_cyl_unit_basis;
      inp->get_con_basis = get_cyl_unit_basis;
      inp->spacetime = 0;
      break;
      
    case WAVE_COORD_SPH:
      inp->mapc2p = mapc2p_sph;
      inp->get_cov_basis = get_sph_unit_basis;
      inp->get_con_basis = get_sph_unit_basis;
      inp->spacetime = 0;
      break;
  }
}

struct gkyl_wave_geom*
gkyl_wave_geom_from_coord_flag(const struct gkyl_rect_grid *grid, struct gkyl_range *range,
  enum gkyl_wave_coord_flag cflag, void *ctx, bool use_gpu)
{
  struct gkyl_wave_geom_inp inp;
  gkyl_wave_geom_inp_from_flag(cflag, &inp);
  struct gkyl_wave_geom *wg = gkyl_wave_geom_from_wg_inp(grid, range, &inp, ctx, use_gpu);
  return wg;
}

struct gkyl_wave_geom*
gkyl_wave_geom_from_wg_inp(const struct gkyl_rect_grid *grid, struct gkyl_range *range,
  struct gkyl_wave_geom_inp *inp, void *ctx, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_wave_geom_cu_dev_from_wg_inp(grid, range, inp, ctx);
  } 
#endif 

  struct gkyl_wave_geom *wg = gkyl_malloc(sizeof(struct gkyl_wave_geom));

  wg->range = *range;
  wg->geom = gkyl_array_new(GKYL_USER, sizeof(struct gkyl_wave_cell_geom), range->volume);

  double xc[GKYL_MAX_CDIM];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {

    gkyl_rect_grid_cell_center(grid, iter.idx, xc);

    struct gkyl_wave_cell_geom *geo = gkyl_array_fetch(wg->geom, gkyl_range_idx(range, iter.idx));

    // compute geometry based on grid dimensions
    switch (grid->ndim) {
      case 1:
        calc_geom_1d(grid->dx, xc, inp, ctx, geo);
        break;

      case 2:
        calc_geom_2d(grid->dx, xc, inp, ctx, geo);
        break;

      case 3:
        calc_geom_3d(grid->dx, xc, inp, ctx, geo);
        break;
    };
  }

  wg->flags = 0;
  GKYL_CLEAR_CU_ALLOC(wg->flags);
  wg->ref_count = gkyl_ref_count_init(gkyl_wave_geom_free);
  wg->on_dev = wg; // CPU eqn obj points to itself

  return wg;
}

struct gkyl_wave_geom*
gkyl_wave_geom_acquire(const struct gkyl_wave_geom* wg)
{
  gkyl_ref_count_inc(&wg->ref_count);
  return (struct gkyl_wave_geom*) wg;
}

void
gkyl_wave_geom_release(const struct gkyl_wave_geom *wg)
{
  gkyl_ref_count_dec(&wg->ref_count);
}
