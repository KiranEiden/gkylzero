#include <gkyl_dg_diffusion_kernels.h> 
GKYL_CU_DH void dg_diffusion4_surfy_2x_tensor_p2(const double* w, const double* dx, double D, 
  const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: Cell-center coordinates
  // dxv[NDIM]: Cell spacing
  // D: Diffusion coefficient in the center cell
  // ql: Input field in the left cell
  // qc: Input field in the center cell
  // qr: Input field in the right cell
  // out: Incremented output

  const double dx1 = 2.0/dx[1]; 
  const double J = -1.0*pow(dx1, 4.0);

  const double *q0l = &ql[0]; 
  const double *q0c = &qc[0]; 
  const double *q0r = &qr[0]; 
  double *out0= &out[0]; 

  out0[0] += J*D*((-6.708203932499369*q0r[5])-6.708203932499369*q0l[5]+13.41640786499874*q0c[5]+8.11898816047911*q0r[2]-8.11898816047911*q0l[2]-4.6875*q0r[0]-4.6875*q0l[0]+9.375*q0c[0]); 
  out0[1] += J*D*((-6.708203932499369*q0r[7])-6.708203932499369*q0l[7]+13.41640786499874*q0c[7]+8.11898816047911*q0r[3]-8.11898816047911*q0l[3]-4.6875*q0r[1]-4.6875*q0l[1]+9.375*q0c[1]); 
  out0[2] += J*D*((-9.077304717673634*q0r[5])+9.077304717673634*q0l[5]+12.65625*q0r[2]+12.65625*q0l[2]+30.9375*q0c[2]-8.11898816047911*q0r[0]+8.11898816047911*q0l[0]); 
  out0[3] += J*D*((-9.077304717673634*q0r[7])+9.077304717673634*q0l[7]+12.65625*q0r[3]+12.65625*q0l[3]+30.9375*q0c[3]-8.11898816047911*q0r[1]+8.11898816047911*q0l[1]); 
  out0[4] += J*D*((-6.708203932499369*q0r[8])-6.708203932499369*q0l[8]+13.41640786499874*q0c[8]+8.118988160479114*q0r[6]-8.118988160479114*q0l[6]-4.6875*q0r[4]-4.6875*q0l[4]+9.375*q0c[4]); 
  out0[5] += J*D*((-0.65625*q0r[5])-0.65625*q0l[5]+40.6875*q0c[5]+4.720198453190289*q0r[2]-4.720198453190289*q0l[2]-4.192627457812106*q0r[0]-4.192627457812106*q0l[0]+8.385254915624213*q0c[0]); 
  out0[6] += J*D*((-9.077304717673634*q0r[8])+9.077304717673634*q0l[8]+12.65625*q0r[6]+12.65625*q0l[6]+30.9375*q0c[6]-8.118988160479114*q0r[4]+8.118988160479114*q0l[4]); 
  out0[7] += J*D*((-0.65625*q0r[7])-0.65625*q0l[7]+40.6875*q0c[7]+4.72019845319029*q0r[3]-4.72019845319029*q0l[3]-4.192627457812105*q0r[1]-4.192627457812105*q0l[1]+8.38525491562421*q0c[1]); 
  out0[8] += J*D*((-0.65625*q0r[8])-0.65625*q0l[8]+40.6875*q0c[8]+4.72019845319029*q0r[6]-4.72019845319029*q0l[6]-4.192627457812106*q0r[4]-4.192627457812106*q0l[4]+8.385254915624213*q0c[4]); 

} 
