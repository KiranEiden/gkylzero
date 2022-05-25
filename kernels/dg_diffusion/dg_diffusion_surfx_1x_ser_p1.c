#include <gkyl_dg_diffusion_kernels.h>

GKYL_CU_DH void
dg_diffusion_surfx_1x_ser_p1(const double* w, const double* dx,
  const double* D,
  const double* ql, const double* qc, const double* qr,
  double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dxv[NDIM]: Cell spacing
  // Dl: Diffusion coefficient in the left cell
  // Dc: Diffusion coefficient in the center cell
  // Dr: Diffusion coefficient in the right cell
  // ql: Input field in the left cell
  // qc: Input field in the center cell
  // qr: Input field in the right cell
  // out: Incremented output

  const double J = 4/dx[0]/dx[0];

  out[0] += (-0.6629126073623879*D[1]*qr[1]-0.3827327723098713*D[0]*qr[1]-0.6629126073623879*D[1]*ql[1]+0.3827327723098713*D[0]*ql[1]-1.325825214724776*D[1]*qc[1]+0.6889189901577684*qr[0]*D[1]-0.6889189901577684*ql[0]*D[1]+0.3977475644174327*D[0]*qr[0]+0.3977475644174327*D[0]*ql[0]-0.7954951288348654*D[0]*qc[0])*J;
  out[1] += (-0.5358258812338196*D[1]*qr[1]-0.3093592167691142*D[0]*qr[1]+0.5358258812338196*D[1]*ql[1]-0.3093592167691142*D[0]*ql[1]-2.032931995911323*D[0]*qc[1]+0.6629126073623878*qr[0]*D[1]+0.6629126073623878*ql[0]*D[1]-3.447145558284417*qc[0]*D[1]+0.3827327723098712*D[0]*qr[0]-0.3827327723098712*D[0]*ql[0])*J;
}
