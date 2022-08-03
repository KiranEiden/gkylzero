// Thu Jul 28 13:07:03 2022
#include <gkyl_binop_cross_mul_hyb.h>
GKYL_CU_DH
void
binop_cross_mul_2x1v_hyb_p1(const double *f, const double *g, double *fg )
{
  double tmp[12] = {0.};
  tmp[0] =  5.0000000000000000e-01*f[1]*g[1]+5.0000000000000000e-01*f[3]*g[4]+5.0000000000000000e-01*f[2]*g[2]+5.0000000000000000e-01*f[0]*g[0];
  tmp[1] =  5.0000000000000000e-01*f[1]*g[0]+5.0000000000000000e-01*f[0]*g[1]+5.0000000000000000e-01*f[3]*g[2]+5.0000000000000000e-01*f[2]*g[4];
  tmp[2] =  5.0000000000000000e-01*g[4]*f[1]+5.0000000000000000e-01*f[3]*g[1]+5.0000000000000000e-01*f[0]*g[2]+5.0000000000000000e-01*f[2]*g[0];
  tmp[3] =  5.0000000000000000e-01*f[3]*g[7]+5.0000000000000000e-01*f[1]*g[5]+5.0000000000000000e-01*f[2]*g[6]+5.0000000000000000e-01*f[0]*g[3];
  tmp[4] =  5.0000000000000000e-01*f[3]*g[0]+5.0000000000000000e-01*g[4]*f[0]+5.0000000000000000e-01*f[2]*g[1]+5.0000000000000000e-01*f[1]*g[2];
  tmp[5] =  5.0000000000000000e-01*g[5]*f[0]+5.0000000000000000e-01*f[3]*g[6]+5.0000000000000000e-01*f[1]*g[3]+5.0000000000000000e-01*f[2]*g[7];
  tmp[6] =  5.0000000000000000e-01*f[3]*g[5]+5.0000000000000000e-01*f[2]*g[3]+5.0000000000000000e-01*g[7]*f[1]+5.0000000000000000e-01*g[6]*f[0];
  tmp[7] =  5.0000000000000000e-01*g[7]*f[0]+5.0000000000000000e-01*f[2]*g[5]+5.0000000000000000e-01*g[6]*f[1]+5.0000000000000000e-01*f[3]*g[3];
  tmp[8] =  5.0000000000000000e-01*f[2]*g[10]+5.0000000000000000e-01*f[0]*g[8]+5.0000000000000000e-01*f[3]*g[11]+5.0000000000000000e-01*f[1]*g[9];
  tmp[9] =  5.0000000000000000e-01*f[2]*g[11]+5.0000000000000000e-01*f[1]*g[8]+5.0000000000000000e-01*f[0]*g[9]+5.0000000000000000e-01*g[10]*f[3];
  tmp[10] =  5.0000000000000000e-01*g[10]*f[0]+5.0000000000000000e-01*f[3]*g[9]+5.0000000000000000e-01*f[1]*g[11]+5.0000000000000000e-01*f[2]*g[8];
  tmp[11] =  5.0000000000000000e-01*f[2]*g[9]+5.0000000000000000e-01*f[3]*g[8]+5.0000000000000000e-01*g[10]*f[1]+5.0000000000000000e-01*f[0]*g[11];
 
  fg[0] = tmp[0];
  fg[1] = tmp[1];
  fg[2] = tmp[2];
  fg[3] = tmp[3];
  fg[4] = tmp[4];
  fg[5] = tmp[5];
  fg[6] = tmp[6];
  fg[7] = tmp[7];
  fg[8] = tmp[8];
  fg[9] = tmp[9];
  fg[10] = tmp[10];
  fg[11] = tmp[11];
}

