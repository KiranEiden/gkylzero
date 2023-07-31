#pragma once

#include <gkyl_dg_gen_diffusion_kernels.h>
#include <gkyl_ref_count.h>
#include <gkyl_dg_eqn.h>

// private header for use in gen_diffusion DG equation object creation
// functions

// Types for various kernels
typedef double (*gen_diffusion_surf_t)(const double *w, const double *dx,
  const double* D, const double *q[27],
  double* GKYL_RESTRICT out);

// for use in kernel tables
typedef struct { vol_termf_t kernels[3]; } gkyl_dg_gen_diffusion_vol_kern_list;
typedef struct { gen_diffusion_surf_t kernels[3]; } gkyl_dg_gen_diffusion_surf_kern_list;

struct dg_gen_diffusion {
  struct gkyl_dg_eqn eqn;
  gen_diffusion_surf_t surf[3][3];
  struct gkyl_range conf_range;
  struct gkyl_dg_gen_diffusion_auxfields auxfields;
};

//
// Serendipity volume kernels
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static double
kernel_dg_gen_diffusion_vol_2x_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_gen_diffusion* gen_diffusion = container_of(eqn, struct dg_gen_diffusion, eqn);
  
  long cidx = gkyl_range_idx(&gen_diffusion->conf_range, idx);
  
  return dg_gen_diffusion_vol_2x_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(gen_diffusion->auxfields.Dij, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_dg_gen_diffusion_vol_2x_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_gen_diffusion* gen_diffusion = container_of(eqn, struct dg_gen_diffusion, eqn);
  
  long cidx = gkyl_range_idx(&gen_diffusion->conf_range, idx);
  
  return dg_gen_diffusion_vol_2x_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(gen_diffusion->auxfields.Dij, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_dg_gen_diffusion_vol_3x_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_gen_diffusion* gen_diffusion = container_of(eqn, struct dg_gen_diffusion, eqn);
  
  long cidx = gkyl_range_idx(&gen_diffusion->conf_range, idx);
  
  return dg_gen_diffusion_vol_3x_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(gen_diffusion->auxfields.Dij, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_dg_gen_diffusion_vol_3x_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_gen_diffusion* gen_diffusion = container_of(eqn, struct dg_gen_diffusion, eqn);
  
  long cidx = gkyl_range_idx(&gen_diffusion->conf_range, idx);
  
  return dg_gen_diffusion_vol_3x_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(gen_diffusion->auxfields.Dij, cidx),
    qIn, qRhsOut);
}

// Volume kernel list
GKYL_CU_D
static const gkyl_dg_gen_diffusion_vol_kern_list ser_vol_kernels[] = {
  { NULL, NULL, NULL }, // general diffusion is not implemented for 1D, use diffusion instead
  { NULL, kernel_dg_gen_diffusion_vol_2x_ser_p1, kernel_dg_gen_diffusion_vol_2x_ser_p2 },
  { NULL, kernel_dg_gen_diffusion_vol_3x_ser_p1, kernel_dg_gen_diffusion_vol_3x_ser_p2 },
};

// Surface kernel list: xx-direction
GKYL_CU_D
static const gkyl_dg_gen_diffusion_surf_kern_list ser_surf_xx_kernels[] = {
  { NULL, NULL, NULL }, // general diffusion is not implemented for 1D, use diffusion instead
  { NULL, dg_gen_diffusion_surfxx_2x_ser_p1, dg_gen_diffusion_surfxx_2x_ser_p2 },
  { NULL, dg_gen_diffusion_surfxx_3x_ser_p1, dg_gen_diffusion_surfxx_3x_ser_p2 },
};

// Surface kernel list: xy-direction
GKYL_CU_D
static const gkyl_dg_gen_diffusion_surf_kern_list ser_surf_xy_kernels[] = {
  { NULL, NULL, NULL },// general diffusion is not implemented for 1D, use diffusion instead
  { NULL, dg_gen_diffusion_surfxy_2x_ser_p1, dg_gen_diffusion_surfxy_2x_ser_p2 },
  { NULL, dg_gen_diffusion_surfxy_3x_ser_p1, dg_gen_diffusion_surfxy_3x_ser_p2 },
};

// Surface kernel list: xz-direction
GKYL_CU_D
static const gkyl_dg_gen_diffusion_surf_kern_list ser_surf_xz_kernels[] = {
  { NULL, NULL, NULL }, // general diffusion is not implemented for 1D, use diffusion instead
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, dg_gen_diffusion_surfxz_3x_ser_p1, dg_gen_diffusion_surfxz_3x_ser_p2 },
};

// Surface kernel list: yx-direction
GKYL_CU_D
static const gkyl_dg_gen_diffusion_surf_kern_list ser_surf_yx_kernels[] = {
  { NULL, NULL, NULL }, // general diffusion is not implemented for 1D, use diffusion instead
  { NULL, dg_gen_diffusion_surfyx_2x_ser_p1, dg_gen_diffusion_surfyx_2x_ser_p2 },
  { NULL, dg_gen_diffusion_surfyx_3x_ser_p1, dg_gen_diffusion_surfyx_3x_ser_p2 },
};

// Surface kernel list: yy-direction
GKYL_CU_D
static const gkyl_dg_gen_diffusion_surf_kern_list ser_surf_yy_kernels[] = {
  { NULL, NULL, NULL }, // general diffusion is not implemented for 1D, use diffusion instead
  { NULL, dg_gen_diffusion_surfyy_2x_ser_p1, dg_gen_diffusion_surfyy_2x_ser_p2 },
  { NULL, dg_gen_diffusion_surfyy_3x_ser_p1, dg_gen_diffusion_surfyy_3x_ser_p2 },
};

// Surface kernel list: yz-direction
GKYL_CU_D
static const gkyl_dg_gen_diffusion_surf_kern_list ser_surf_yz_kernels[] = {
  { NULL, NULL, NULL }, // general diffusion is not implemented for 1D, use diffusion instead
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, dg_gen_diffusion_surfyz_3x_ser_p1, dg_gen_diffusion_surfyz_3x_ser_p2 },
};

// Surface kernel list: zx-direction
GKYL_CU_D
static const gkyl_dg_gen_diffusion_surf_kern_list ser_surf_zx_kernels[] = {
  { NULL, NULL, NULL }, // general diffusion is not implemented for 1D, use diffusion instead
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, dg_gen_diffusion_surfzx_3x_ser_p1, dg_gen_diffusion_surfzx_3x_ser_p2 },
};

// Surface kernel list: zy-direction
GKYL_CU_D
static const gkyl_dg_gen_diffusion_surf_kern_list ser_surf_zy_kernels[] = {
  { NULL, NULL, NULL }, // general diffusion is not implemented for 1D, use diffusion instead
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, dg_gen_diffusion_surfzy_3x_ser_p1, dg_gen_diffusion_surfzy_3x_ser_p2 },
};

// Surface kernel list: zz-direction
GKYL_CU_D
static const gkyl_dg_gen_diffusion_surf_kern_list ser_surf_zz_kernels[] = {
  { NULL, NULL, NULL }, // general diffusion is not implemented for 1D, use diffusion instead
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, dg_gen_diffusion_surfzz_3x_ser_p1, dg_gen_diffusion_surfzz_3x_ser_p2 },
};

/**
 * Free gen_diffusion equation object
 *
 * @param ref Reference counter for constant gen_diffusion equation
 */
void gkyl_gen_diffusion_free(const struct gkyl_ref_count* ref);

GKYL_CU_D
static double
surf(const struct gkyl_dg_eqn* eqn, int dir1, int dir2,
  const double* xc, const double* dxc, const int* idxc,
  long sz_dim, const int idx[27][GKYL_MAX_DIM], const double* qIn[27],
  double* GKYL_RESTRICT qRhsOut)
{
  struct dg_gen_diffusion* gen_diffusion = container_of(eqn, struct dg_gen_diffusion, eqn);
  long cidx = gkyl_range_idx(&gen_diffusion->conf_range, idxc);
  
  return gen_diffusion->surf[dir1][dir2](xc, dxc,
    (const double*) gkyl_array_cfetch(gen_diffusion->auxfields.Dij, cidx), 
    qIn, qRhsOut);
}
