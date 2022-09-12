#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_prim_lbo_cross_calc.h>
#include <gkyl_prim_lbo_cross_calc_priv.h>
#include <gkyl_prim_lbo_kernels.h> 
#include <gkyl_prim_lbo_gyrokinetic.h>
#include <gkyl_prim_lbo_vlasov.h>
#include <gkyl_mat.h>
#include <assert.h>

struct gkyl_prim_lbo_cross_calc*
gkyl_prim_lbo_cross_calc_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *cbasis, const struct gkyl_range *conf_rng,  
  struct gkyl_prim_lbo_type *prim, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_prim_lbo_cross_calc_cu_dev_new(grid, prim);
  } 
#endif   
  gkyl_prim_lbo_cross_calc *up = gkyl_malloc(sizeof(gkyl_prim_lbo_cross_calc));
  up->grid = *grid;
  up->prim = gkyl_prim_lbo_type_acquire(prim);

  // allocate memory for use in kernels
  int nc = up->prim->num_config;
  int udim = up->prim->udim;
  int N = nc*(udim + 1);

  up->As = gkyl_nmat_new(conf_rng->volume, N, N);
  up->xs = gkyl_nmat_new(conf_rng->volume, N, 1);
  up->mem = gkyl_nmat_linsolve_lu_new(up->As->num, up->As->nr);

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host
  
  return up;
}

void
gkyl_prim_lbo_cross_calc_advance(const struct gkyl_prim_lbo_cross_calc* calc,
  const struct gkyl_range *conf_rng,
  const struct gkyl_array *greene,
  double self_m, const struct gkyl_array *self_u, const struct gkyl_array *self_vtsq,
  double cross_m, const struct gkyl_array *cross_u, const struct gkyl_array *cross_vtsq, 
  const struct gkyl_array *moms, const struct gkyl_array *boundary_corrections, 
  struct gkyl_array *u_out, struct gkyl_array *vtsq_out)
{
  int nc = calc->prim->num_config;
  int udim = calc->prim->udim;

  struct gkyl_range_iter conf_iter;
  
  // loop over configuration space cells.
  gkyl_range_iter_init(&conf_iter, conf_rng);
  long count = 0;
  while (gkyl_range_iter_next(&conf_iter)) {
    long midx = gkyl_range_idx(conf_rng, conf_iter.idx);

    struct gkyl_mat lhs = gkyl_nmat_get(calc->As, count);
    struct gkyl_mat rhs = gkyl_nmat_get(calc->xs, count);

    gkyl_mat_clear(&lhs, 0.0); gkyl_mat_clear(&rhs, 0.0);

    calc->prim->cross_prim(calc->prim, &lhs, &rhs, conf_iter.idx, gkyl_array_cfetch(greene, midx),
      self_m, gkyl_array_cfetch(self_u, midx), gkyl_array_cfetch(self_vtsq, midx),
      cross_m, gkyl_array_cfetch(cross_u, midx), gkyl_array_cfetch(cross_vtsq, midx),
      gkyl_array_cfetch(moms, midx), gkyl_array_cfetch(boundary_corrections, midx)
    );

    count += 1;
  }

  bool status = gkyl_nmat_linsolve_lu_pa(calc->mem, calc->As, calc->xs);
  gkyl_range_iter_init(&conf_iter, conf_rng);
  count = 0;
  while (gkyl_range_iter_next(&conf_iter)) {
    long midx = gkyl_range_idx(conf_rng, conf_iter.idx);

    struct gkyl_mat out = gkyl_nmat_get(calc->xs, count);
    prim_lbo_copy_sol(&out, nc, udim, 
      gkyl_array_fetch(u_out, midx), gkyl_array_fetch(vtsq_out, midx)
    );
    count += 1;

  }
}

const struct gkyl_prim_lbo_type* 
gkyl_prim_lbo_cross_calc_get_prim(struct gkyl_prim_lbo_cross_calc* calc)
{
  return calc->prim;
}

void 
gkyl_prim_lbo_cross_calc_release(struct gkyl_prim_lbo_cross_calc* up)
{
  gkyl_prim_lbo_type_release(up->prim);

  if (up->As)
    gkyl_nmat_release(up->As);
  if (up->xs)
    gkyl_nmat_release(up->xs);
  if (up->mem)
    gkyl_nmat_linsolve_lu_release(up->mem);
  
  if (GKYL_IS_CU_ALLOC(up->flags))
    gkyl_cu_free(up->on_dev);
  gkyl_free(up);
}

// "derived" class constructors
struct gkyl_prim_lbo_cross_calc* 
gkyl_prim_lbo_vlasov_cross_calc_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis,
  const struct gkyl_range *conf_rng, bool use_gpu)
{
  struct gkyl_prim_lbo_type *prim; // LBO primitive moments type
  prim = gkyl_prim_lbo_vlasov_new(cbasis, pbasis, use_gpu);
  struct gkyl_prim_lbo_cross_calc *calc = gkyl_prim_lbo_cross_calc_new(grid, cbasis, conf_rng, prim, use_gpu);
  gkyl_prim_lbo_type_release(prim);
  return calc;
}

struct gkyl_prim_lbo_cross_calc* 
gkyl_prim_lbo_gyrokinetic_cross_calc_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis,
  const struct gkyl_range *conf_rng, bool use_gpu)
{
  struct gkyl_prim_lbo_type *prim; // LBO primitive moments type
  prim = gkyl_prim_lbo_gyrokinetic_new(cbasis, pbasis, use_gpu);
  struct gkyl_prim_lbo_cross_calc *calc = gkyl_prim_lbo_cross_calc_new(grid, cbasis, conf_rng, prim, use_gpu);
  gkyl_prim_lbo_type_release(prim);
  return calc;
}

#ifndef GKYL_HAVE_CUDA

gkyl_prim_lbo_cross_calc*
gkyl_prim_lbo_cross_calc_cu_dev_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *cbasis, const struct gkyl_range *conf_rng, 
  struct gkyl_prim_lbo_type *prim)
{
  assert(false);
  return 0;
}

void
gkyl_prim_lbo_cross_calc_advance_cu(const struct gkyl_prim_lbo_cross_calc* calc,
  const struct gkyl_range *conf_rng,
  const struct gkyl_array *greene,
  double self_m, const struct gkyl_array *self_u, const struct gkyl_array *self_vtsq,
  double cross_m, const struct gkyl_array *cross_u, const struct gkyl_array *cross_vtsq, 
  const struct gkyl_array *moms, const struct gkyl_array *boundary_corrections, 
  struct gkyl_array *u_out, struct gkyl_array *vtsq_out)
{
  assert(false);
}

#endif
