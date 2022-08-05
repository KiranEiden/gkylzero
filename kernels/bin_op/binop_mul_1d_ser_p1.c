// Thu Aug 26 15:51:37 2021
#include <gkyl_binop_mul_ser.h>
GKYL_CU_DH
void
binop_mul_1d_ser_p1(const double *f, const double *g, double *fg )
{
  fg[0] =  7.0710678118654757e-01*g[0]*f[0]+7.0710678118654757e-01*g[1]*f[1];
  fg[1] =  7.0710678118654757e-01*g[0]*f[1]+7.0710678118654757e-01*g[1]*f[0];
  // nsum = 2, nprod = 8
}

struct gkyl_kern_op_count op_count_binop_mul_1d_ser_p1(void)
{
  return (struct gkyl_kern_op_count) { .num_sum = 2, .num_prod = 8 };
}
