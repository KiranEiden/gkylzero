// Thu Jul 28 13:09:32 2022
#include <gkyl_binop_cross_mul_gkhyb.h>
GKYL_CU_DH
void
binop_cross_mul_1x1v_gkhyb_p1(const double *f, const double *g, double *fg )
{
  double tmp[6] = {0.};
  tmp[0] =  7.0710678118654757e-01*g[1]*f[1]+7.0710678118654757e-01*g[0]*f[0];
  tmp[1] =  7.0710678118654757e-01*g[1]*f[0]+7.0710678118654757e-01*g[0]*f[1];
  tmp[2] =  7.0710678118654757e-01*g[2]*f[0]+7.0710678118654757e-01*g[3]*f[1];
  tmp[3] =  7.0710678118654757e-01*g[2]*f[1]+7.0710678118654757e-01*g[3]*f[0];
  tmp[4] =  7.0710678118654757e-01*f[0]*g[4]+7.0710678118654757e-01*g[5]*f[1];
  tmp[5] =  7.0710678118654757e-01*f[1]*g[4]+7.0710678118654757e-01*g[5]*f[0];
 
  fg[0] = tmp[0];
  fg[1] = tmp[1];
  fg[2] = tmp[2];
  fg[3] = tmp[3];
  fg[4] = tmp[4];
  fg[5] = tmp[5];
}

