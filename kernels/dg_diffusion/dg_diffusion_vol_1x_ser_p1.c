#include <gkyl_dg_diffusion_kernels.h>

GKYL_CU_DH double
dg_diffusion_vol_1x_ser_p1(const double* w, const double* dx,
  const double* D, const double* q, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dxv[NDIM]: Cell spacing
  // D: Diffusion tensor
  // q: Input field
  // out: Incremented output

  const double J = 4/dx[0]/dx[0];

  out[0] += (0.0)*J;
  out[1] += (2.121320343559642*q[0]*D[1])*J;

  return D[0]*J;
}
