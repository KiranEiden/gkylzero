// Thu Aug 26 15:51:57 2021
#include <gkyl_binop_cross_mul_ser.h>
GKYL_CU_DH
void
binop_cross_mul_1d_2d_ser_p2(const double *f, const double *g, double *fg )
{
  fg[0] =  7.0710678118654757e-01*f[2]*g[4]+7.0710678118654757e-01*g[1]*f[1]+7.0710678118654757e-01*g[0]*f[0];
  fg[1] =  6.3245553203367588e-01*g[1]*f[2]+6.3245553203367588e-01*f[1]*g[4]+7.0710678118654757e-01*g[1]*f[0]+7.0710678118654757e-01*g[0]*f[1];
  fg[2] =  7.0710678118654757e-01*f[2]*g[6]+7.0710678118654757e-01*g[2]*f[0]+7.0710678118654757e-01*g[3]*f[1];
  fg[3] =  6.3245553203367588e-01*g[3]*f[2]+7.0710678118654757e-01*g[2]*f[1]+6.3245553203367588e-01*g[6]*f[1]+7.0710678118654757e-01*g[3]*f[0];
  fg[4] =  7.0710678118654757e-01*g[0]*f[2]+4.5175395145262565e-01*f[2]*g[4]+6.3245553203367588e-01*g[1]*f[1]+7.0710678118654757e-01*f[0]*g[4];
  fg[5] =  7.0710678118654757e-01*f[1]*g[7]+7.0710678118654757e-01*f[0]*g[5];
  fg[6] =  4.5175395145262565e-01*f[2]*g[6]+7.0710678118654757e-01*g[6]*f[0]+7.0710678118654757e-01*g[2]*f[2]+6.3245553203367588e-01*g[3]*f[1];
  fg[7] =  6.3245553203367588e-01*f[2]*g[7]+7.0710678118654757e-01*f[0]*g[7]+7.0710678118654757e-01*f[1]*g[5];
}

