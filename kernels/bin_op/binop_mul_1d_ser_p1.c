// Mon Jul 19 09:31:00 2021
#include <gkyl_binop_mul_ser.h>
GKYL_CU_DH
void
binop_mul_1d_ser_p1(const double *f, const double *g, double *fg )
{
  fg[0] =  7.0710678118654757e-01*f[0]*g[0]+7.0710678118654757e-01*f[1]*g[1];
  fg[1] =  7.0710678118654757e-01*f[0]*g[1]+7.0710678118654757e-01*f[1]*g[0];
}

