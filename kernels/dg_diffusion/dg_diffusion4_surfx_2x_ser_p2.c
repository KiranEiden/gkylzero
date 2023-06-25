#include <gkyl_dg_diffusion_kernels.h> 
GKYL_CU_DH void dg_diffusion4_surfx_2x_ser_p2(const double* w, const double* dx, double D, 
  const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: Cell-center coordinates
  // dxv[NDIM]: Cell spacing
  // D: Diffusion coefficient in the center cell
  // ql: Input field in the left cell
  // qc: Input field in the center cell
  // qr: Input field in the right cell
  // out: Incremented output

  const double dx1 = 2.0/dx[0]; 
  const double J = -1.0*pow(dx1, 4.0);

  const double *q0l = &ql[0]; 
  const double *q0c = &qc[0]; 
  const double *q0r = &qr[0]; 
  double *out0= &out[0]; 

  out0[0] += J*D*((-6.708203932499369*q0r[4])-6.708203932499369*q0l[4]+13.41640786499874*q0c[4]+8.11898816047911*q0r[1]-8.11898816047911*q0l[1]-4.6875*q0r[0]-4.6875*q0l[0]+9.375*q0c[0]); 
  out0[1] += J*D*((-9.077304717673634*q0r[4])+9.077304717673634*q0l[4]+12.65625*q0r[1]+12.65625*q0l[1]+30.9375*q0c[1]-8.11898816047911*q0r[0]+8.11898816047911*q0l[0]); 
  out0[2] += J*D*((-6.708203932499369*q0r[6])-6.708203932499369*q0l[6]+13.41640786499874*q0c[6]+8.11898816047911*q0r[3]-8.11898816047911*q0l[3]-4.6875*q0r[2]-4.6875*q0l[2]+9.375*q0c[2]); 
  out0[3] += J*D*((-9.077304717673634*q0r[6])+9.077304717673634*q0l[6]+12.65625*q0r[3]+12.65625*q0l[3]+30.9375*q0c[3]-8.11898816047911*q0r[2]+8.11898816047911*q0l[2]); 
  out0[4] += J*D*((-0.65625*q0r[4])-0.65625*q0l[4]+40.6875*q0c[4]+4.720198453190289*q0r[1]-4.720198453190289*q0l[1]-4.192627457812106*q0r[0]-4.192627457812106*q0l[0]+8.385254915624213*q0c[0]); 
  out0[5] += J*D*(8.118988160479114*q0r[7]-8.118988160479114*q0l[7]-4.6875*q0r[5]-4.6875*q0l[5]+9.375*q0c[5]); 
  out0[6] += J*D*((-0.65625*q0r[6])-0.65625*q0l[6]+40.6875*q0c[6]+4.72019845319029*q0r[3]-4.72019845319029*q0l[3]-4.192627457812105*q0r[2]-4.192627457812105*q0l[2]+8.38525491562421*q0c[2]); 
  out0[7] += J*D*(12.65625*q0r[7]+12.65625*q0l[7]+30.9375*q0c[7]-8.118988160479114*q0r[5]+8.118988160479114*q0l[5]); 

} 
