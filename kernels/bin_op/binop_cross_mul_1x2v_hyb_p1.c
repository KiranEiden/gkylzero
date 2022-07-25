// Tue Jul 12 14:38:52 2022
#include <gkyl_binop_cross_mul_hyb.h>
GKYL_CU_DH
void
binop_cross_mul_1x2v_hyb_p1(const double *f, const double *g, double *fg )
{
  fg[0] =  7.0710678118654757e-01*f[1]*g[1]+7.0710678118654757e-01*f[0]*g[0];
  fg[1] =  7.0710678118654757e-01*f[1]*g[0]+7.0710678118654757e-01*f[0]*g[1];
  fg[2] =  7.0710678118654757e-01*f[1]*g[4]+7.0710678118654757e-01*g[2]*f[0];
  fg[3] =  7.0710678118654757e-01*g[3]*f[0]+7.0710678118654757e-01*f[1]*g[5];
  fg[4] =  7.0710678118654757e-01*f[0]*g[4]+7.0710678118654757e-01*g[2]*f[1];
  fg[5] =  7.0710678118654757e-01*f[0]*g[5]+7.0710678118654757e-01*g[3]*f[1];
  fg[6] =  7.0710678118654757e-01*f[1]*g[7]+7.0710678118654757e-01*f[0]*g[6];
  fg[7] =  7.0710678118654757e-01*f[1]*g[6]+7.0710678118654757e-01*f[0]*g[7];
  fg[8] =  7.0710678118654757e-01*f[0]*g[8]+7.0710678118654757e-01*f[1]*g[9];
  fg[9] =  7.0710678118654757e-01*f[1]*g[8]+7.0710678118654757e-01*f[0]*g[9];
  fg[10] =  7.0710678118654757e-01*g[11]*f[1]+7.0710678118654757e-01*f[0]*g[10];
  fg[11] =  7.0710678118654757e-01*f[0]*g[11]+7.0710678118654757e-01*f[1]*g[10];
  fg[12] =  7.0710678118654757e-01*g[13]*f[1]+7.0710678118654757e-01*f[0]*g[12];
  fg[13] =  7.0710678118654757e-01*g[13]*f[0]+7.0710678118654757e-01*f[1]*g[12];
  fg[14] =  7.0710678118654757e-01*g[15]*f[1]+7.0710678118654757e-01*f[0]*g[14];
  fg[15] =  7.0710678118654757e-01*g[15]*f[0]+7.0710678118654757e-01*f[1]*g[14];
}

