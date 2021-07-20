// Tue Jul 20 15:17:00 2021
#include <gkyl_binop_mul_ser.h>
GKYL_CU_DH
void
binop_mul_2d_ser_p1(const double *f, const double *g, double *fg )
{
  fg[0] =  5.0000000000000000e-01*g[3]*f[3]+5.0000000000000000e-01*g[1]*f[1]+5.0000000000000000e-01*g[2]*f[2]+5.0000000000000000e-01*g[0]*f[0];
  fg[1] =  5.0000000000000000e-01*g[3]*f[2]+5.0000000000000000e-01*g[1]*f[0]+5.0000000000000000e-01*g[0]*f[1]+5.0000000000000000e-01*g[2]*f[3];
  fg[2] =  5.0000000000000000e-01*g[0]*f[2]+5.0000000000000000e-01*g[1]*f[3]+5.0000000000000000e-01*g[2]*f[0]+5.0000000000000000e-01*g[3]*f[1];
  fg[3] =  5.0000000000000000e-01*g[2]*f[1]+5.0000000000000000e-01*g[1]*f[2]+5.0000000000000000e-01*g[0]*f[3]+5.0000000000000000e-01*g[3]*f[0];
}

