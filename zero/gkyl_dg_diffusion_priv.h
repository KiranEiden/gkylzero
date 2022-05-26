#pragma once

#include <gkyl_dg_diffusion_kernels.h>
#include <gkyl_ref_count.h>
#include <gkyl_dg_eqn.h>

// private header for use in diffusion DG equation object creation
// functions

// Types for various kernels
typedef double (*diffusion_vol_t)(const double* w, const double* dx,
  const double* D, const double* q, double* GKYL_RESTRICT out);

typedef void (*diffusion_surf_t)(const double *w, const double *dx,
  const double* D, const double *ql, const double *qc, const double *qr,
  double* GKYL_RESTRICT out);

// for use in kernel tables
typedef struct { diffusion_vol_t kernels[3]; } gkyl_dg_diffusion_vol_kern_list;
typedef struct { diffusion_surf_t kernels[3]; } gkyl_dg_diffusion_surf_kern_list;

// Volume kernel list
GKYL_CU_D
static const gkyl_dg_diffusion_vol_kern_list vol_kernels[] = {
  { NULL, dg_diffusion_vol_1x_ser_p1, dg_diffusion_vol_1x_ser_p2 }, // 1
  { NULL, NULL, NULL }, // 2
  { NULL, NULL, NULL },
};

// Surface kernel list: x-direction
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list surf_x_kernels[] = {
  { NULL, dg_diffusion_surfx_1x_ser_p1, dg_diffusion_surfx_1x_ser_p2 },
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
};

// Surface kernel list: y-direction
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list surf_y_kernels[] = {
  { NULL, NULL, NULL }, // no y-direction in 1D
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
};

// Surface kernel list: z-direction
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list surf_z_kernels[] = {
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, NULL, NULL },
};

struct dg_diffusion {
  struct gkyl_dg_eqn eqn;
  diffusion_vol_t vol;
  diffusion_surf_t surf[3];
  struct gkyl_range conf_range;
  struct gkyl_dg_diffusion_auxfields auxfields;
};

/**
 * Free diffusion equation object
 *
 * @param ref Reference counter for constant diffusion equation
 */
void gkyl_diffusion_free(const struct gkyl_ref_count* ref);


GKYL_CU_D
static double
vol(const struct gkyl_dg_eqn* eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  
  return diffusion->vol(xc, dx,
    (const double*) gkyl_array_cfetch(diffusion->auxfields.D, cidx),
    qIn, qRhsOut);
}

GKYL_CU_D
static void
surf(const struct gkyl_dg_eqn* eqn, int dir,
  const double* xcL, const double* xcC, const double* xcR, 
  const double* dxL, const double* dxC, const double* dxR,
  const int* idxL, const int* idxC, const int* idxR,
  const double* qInL, const double* qInC, const double* qInR,
  double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  
  long cidx = gkyl_range_idx(&diffusion->conf_range, idxC);
  
  diffusion->surf[dir](xcC, dxC,
    (const double*) gkyl_array_cfetch(diffusion->auxfields.D, cidx), 
    qInL, qInC, qInR, qRhsOut);
}

GKYL_CU_D
static void
boundary_surf(const struct gkyl_dg_eqn* eqn, int dir,
  const double* xcEdge, const double* xcSkin,
  const double* dxEdge, const double* dxSkin,
  const int* idxEdge, const int* idxSkin, const int edge,
  const double* qInEdge, const double* qInSkin, double* GKYL_RESTRICT qRhsOut)
{
  
}
