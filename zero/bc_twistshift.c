#pragma once

#include <gkyl_bc_twistshift.h>
#include <gkyl_bc_twistshift_priv.h>
#include <assert.h>

struct gkyl_bc_twistshift*
gkyl_bc_twistshift_new(int dir, enum gkyl_edge_loc edge,
  const struct gkyl_range *local_range_ext, const int *num_ghosts, const struct gkyl_basis *basis,
  const struct gkyl_rect_grid *grid, int cdim,
  const struct gkyl_array *yshift, const int *ndonors, bool use_gpu) {

  // Allocate space for new updater.
  struct gkyl_bc_twistshift *up = gkyl_malloc(sizeof(struct gkyl_bc_twistshift));

  up->dir = dir;
  up->edge = edge;
  up->basis = basis;
  up->grid = grid;
  up->ndonors = ndonors.
  up->use_gpu = use_gpu;

  // Choose the kernel that does the reflection/no reflection/partial
  // reflection.
  up->kernels = gkyl_malloc(sizeof(struct gkyl_bc_twistshift_kernels));
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    up->kernels_cu = gkyl_cu_malloc(sizeof(struct gkyl_bc_twistshift_kernels));
    gkyl_bc_twistshift_choose_kernel_cu(basis, cdim, up->kernels_cu);
  } else {
    gkyl_bc_twistshift_choose_kernel(basis, cdim, up->kernels);
    up->kernels_cu = up->kernels;
  }
#else
  gkyl_bc_twistshift_choose_kernel(basis, cdim, up->kernels);
  up->kernels_cu = up->kernels;
#endif

  return up;
}

void gkyl_bc_twistshift_integral_xlimdg(struct gkyl_bc_twistshift *up,
  double sFac, const double *xLimLo, const double *xLimUp, double yLimLo, double yLimUp,
  double dyDo, double yOff, double *ySh, struct gkyl_nmat *mats, int cellidx, int doidx) {
}

void gkyl_bc_twistshift_integral_ylimdg(struct gkyl_bc_twistshift *up,
  double sFac, double xLimLo, double xLimUp, const double *yLimLo, const double *yLimUp,
  double dyDo, double yOff, const double *ySh, struct gkyl_nmat *mats, int cellidx, int doidx) {
}

void gkyl_bc_twistshift_integral_fullcelllimdg(struct gkyl_bc_twistshift *up,
  double dyDo, double yOff, const double *ySh, struct gkyl_nmat *mats, int cellidx, int doidx) {
}

void gkyl_bc_twistshift_release(struct gkyl_bc_twistshift *up) {
  // Release memory associated with this updater.
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    gkyl_cu_free(up->kernels_cu);
#endif
  gkyl_free(up->kernels);
  gkyl_free(up);
}
