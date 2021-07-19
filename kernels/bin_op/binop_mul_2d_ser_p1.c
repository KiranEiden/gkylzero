// Mon Jul 19 09:31:00 2021
#include <gkyl_binop_mul_ser.h>
GKYL_CU_DH
void
binop_mul_2d_ser_p1(const double *f, const double *g, double *fg )
{
  fg[0] =  5.0000000000000000e-01*f[0]*g[0]+5.0000000000000000e-01*f[2]*g[2]+5.0000000000000000e-01*f[1]*g[1]+5.0000000000000000e-01*f[3]*g[3];
  fg[1] =  5.0000000000000000e-01*f[2]*g[3]+5.0000000000000000e-01*f[0]*g[1]+5.0000000000000000e-01*f[1]*g[0]+5.0000000000000000e-01*f[3]*g[2];
  fg[2] =  5.0000000000000000e-01*f[3]*g[1]+5.0000000000000000e-01*f[2]*g[0]+5.0000000000000000e-01*f[0]*g[2]+5.0000000000000000e-01*f[1]*g[3];
  fg[3] =  5.0000000000000000e-01*f[3]*g[0]+5.0000000000000000e-01*f[1]*g[2]+5.0000000000000000e-01*f[0]*g[3]+5.0000000000000000e-01*f[2]*g[1];
}

