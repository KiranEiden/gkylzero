/* -*- c++ -*- */
extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_mat.h>
#include <gkyl_prim_lbo_cross_calc.h>
#include <gkyl_prim_lbo_cross_calc_priv.h>
#include <gkyl_prim_lbo_kernels.h> 
#include <gkyl_prim_lbo_gyrokinetic.h>
#include <gkyl_prim_lbo_vlasov.h>
#include <gkyl_util.h>
#include <stdio.h>
}

__global__ static void
gkyl_prim_lbo_cross_calc_set_cu_ker(gkyl_prim_lbo_cross_calc* calc,
  struct gkyl_nmat *As, struct gkyl_nmat *xs,
  struct gkyl_basis cbasis, const struct gkyl_range conf_rng,
  const struct gkyl_array *greene,
  double self_m, const struct gkyl_array *self_u, const struct gkyl_array *self_vtsq,
  double cross_m, const struct gkyl_array *cross_u, const struct gkyl_array *cross_vtsq, 
  const struct gkyl_array *moms, const struct gkyl_array *boundary_corrections)
{
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_rng.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_rng, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long start = gkyl_range_idx(&conf_rng, idx);

    struct gkyl_mat lhs = gkyl_nmat_get(As, linc1);
    struct gkyl_mat rhs = gkyl_nmat_get(xs, linc1);

    const double *self_u_d = (const double*) gkyl_array_cfetch(self_u, start);
    const double *greene_d = (const double*) gkyl_array_cfetch(greene, start);
    const double *self_vtsq_d = (const double*) gkyl_array_cfetch(self_vtsq, start);
    const double *cross_u_d = (const double*) gkyl_array_cfetch(cross_u, start);
    const double *cross_vtsq_d = (const double*) gkyl_array_cfetch(cross_vtsq, start);
    const double *moms_d = (const double*) gkyl_array_cfetch(moms, start);
    const double *boundary_corrections_d = (const double*) gkyl_array_cfetch(boundary_corrections, start);

    gkyl_mat_clear(&lhs, 0.0); gkyl_mat_clear(&rhs, 0.0);

    calc->prim->cross_prim(calc->prim, &lhs, &rhs, idx, greene_d, 
      self_m, self_u_d, self_vtsq_d,
      cross_m, cross_u_d, cross_vtsq_d, 
      moms_d, boundary_corrections_d
    );
  }
}

__global__ static void
gkyl_prim_lbo_copy_sol_cu_ker(struct gkyl_nmat *xs,
  struct gkyl_basis cbasis, const struct gkyl_range conf_rng,
  int nc, int udim, 
  struct gkyl_array* u_out, struct gkyl_array* vtsq_out)
{
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_rng.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_rng, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long start = gkyl_range_idx(&conf_rng, idx);

    struct gkyl_mat out_d = gkyl_nmat_get(xs, linc1);
    double *u_d = (double*) gkyl_array_fetch(u_out, start);
    double *vtsq_d = (double*) gkyl_array_fetch(vtsq_out, start);
    
    prim_lbo_copy_sol(&out_d, nc, udim, u_d, vtsq_d);
  }
}

void
gkyl_prim_lbo_cross_calc_advance_cu(gkyl_prim_lbo_cross_calc* calc,
  struct gkyl_basis cbasis, const struct gkyl_range *conf_rng,
  const struct gkyl_array *greene,
  double self_m, const struct gkyl_array *self_u, const struct gkyl_array *self_vtsq,
  double cross_m, const struct gkyl_array *cross_u, const struct gkyl_array *cross_vtsq, 
  const struct gkyl_array *moms, const struct gkyl_array *boundary_corrections, 
  struct gkyl_array *u_out, struct gkyl_array *vtsq_out)
{
  int nc = cbasis.num_basis;
  int udim = calc->prim->udim;
  int N = nc*(udim + 1);
  
  if (calc->is_first) {
    calc->As = gkyl_nmat_cu_dev_new(conf_rng->volume, N, N);
    calc->xs = gkyl_nmat_cu_dev_new(conf_rng->volume, N, 1);
    calc->mem = gkyl_nmat_linsolve_lu_cu_dev_new(calc->As->num, calc->As->nr);
    calc->is_first = false;
  }

  gkyl_prim_lbo_cross_calc_set_cu_ker<<<conf_rng->nblocks, conf_rng->nthreads>>>(calc->on_dev,
    calc->As->on_dev, calc->xs->on_dev, 
    cbasis, *conf_rng, 
    greene->on_dev, 
    self_m, self_u->on_dev, self_vtsq->on_dev,
    cross_m, cross_u->on_dev, cross_vtsq->on_dev,
    moms->on_dev, boundary_corrections->on_dev);
  
  bool status = gkyl_nmat_linsolve_lu_pa(calc->mem, calc->As, calc->xs);
  
  gkyl_prim_lbo_copy_sol_cu_ker<<<conf_rng->nblocks, conf_rng->nthreads>>>(calc->xs->on_dev,
    cbasis, *conf_rng, nc, udim, 
    u_out->on_dev, vtsq_out->on_dev);
}

gkyl_prim_lbo_cross_calc*
gkyl_prim_lbo_cross_calc_cu_dev_new(const struct gkyl_rect_grid *grid,
  struct gkyl_prim_lbo_type *prim)
{
  gkyl_prim_lbo_cross_calc *up = (gkyl_prim_lbo_cross_calc*) gkyl_malloc(sizeof(gkyl_prim_lbo_cross_calc));
  up->grid = *grid;
  up->prim = prim;

  up->is_first = true;
  up->As = up->xs = 0;
  up->mem = 0;

  struct gkyl_prim_lbo_type *pt = gkyl_prim_lbo_type_acquire(prim);
  up->prim = pt->on_dev; // so memcpy below gets dev copy

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  gkyl_prim_lbo_cross_calc *up_cu = (gkyl_prim_lbo_cross_calc*) gkyl_cu_malloc(sizeof(gkyl_prim_lbo_cross_calc));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_prim_lbo_cross_calc), GKYL_CU_MEMCPY_H2D);

  up->prim = pt; // host portion of struct should have host copy
  up->on_dev = up_cu; // host pointer
  
  return up;
}
