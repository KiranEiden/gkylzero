// Tue Jul 27 09:22:21 2021
#include <gkyl_binop_mul_ser.h>
GKYL_CU_DH
void
binop_mul_3d_ser_p1(const double *f, const double *g, double *fg )
{
  fg[0] =  3.5355339059327379e-01*f[5]*g[5]+3.5355339059327379e-01*g[3]*f[3]+3.5355339059327379e-01*f[4]*g[4]+3.5355339059327379e-01*f[6]*g[6]+3.5355339059327379e-01*f[1]*g[1]+3.5355339059327379e-01*f[0]*g[0]+3.5355339059327379e-01*g[2]*f[2]+3.5355339059327379e-01*f[7]*g[7];
  fg[1] =  3.5355339059327379e-01*f[6]*g[7]+3.5355339059327379e-01*g[4]*f[2]+3.5355339059327379e-01*g[5]*f[3]+3.5355339059327379e-01*f[1]*g[0]+3.5355339059327379e-01*f[5]*g[3]+3.5355339059327379e-01*f[7]*g[6]+3.5355339059327379e-01*f[4]*g[2]+3.5355339059327379e-01*f[0]*g[1];
  fg[2] =  3.5355339059327379e-01*f[0]*g[2]+3.5355339059327379e-01*f[6]*g[3]+3.5355339059327379e-01*f[4]*g[1]+3.5355339059327379e-01*f[5]*g[7]+3.5355339059327379e-01*f[2]*g[0]+3.5355339059327379e-01*g[6]*f[3]+3.5355339059327379e-01*g[4]*f[1]+3.5355339059327379e-01*f[7]*g[5];
  fg[3] =  3.5355339059327379e-01*g[6]*f[2]+3.5355339059327379e-01*f[4]*g[7]+3.5355339059327379e-01*g[0]*f[3]+3.5355339059327379e-01*g[5]*f[1]+3.5355339059327379e-01*f[6]*g[2]+3.5355339059327379e-01*f[5]*g[1]+3.5355339059327379e-01*f[7]*g[4]+3.5355339059327379e-01*g[3]*f[0];
  fg[4] =  3.5355339059327379e-01*f[4]*g[0]+3.5355339059327379e-01*g[7]*f[3]+3.5355339059327379e-01*g[2]*f[1]+3.5355339059327379e-01*g[4]*f[0]+3.5355339059327379e-01*f[7]*g[3]+3.5355339059327379e-01*f[2]*g[1]+3.5355339059327379e-01*f[5]*g[6]+3.5355339059327379e-01*f[6]*g[5];
  fg[5] =  3.5355339059327379e-01*g[5]*f[0]+3.5355339059327379e-01*f[7]*g[2]+3.5355339059327379e-01*f[4]*g[6]+3.5355339059327379e-01*f[5]*g[0]+3.5355339059327379e-01*g[1]*f[3]+3.5355339059327379e-01*g[7]*f[2]+3.5355339059327379e-01*g[3]*f[1]+3.5355339059327379e-01*f[6]*g[4];
  fg[6] =  3.5355339059327379e-01*f[5]*g[4]+3.5355339059327379e-01*f[4]*g[5]+3.5355339059327379e-01*g[3]*f[2]+3.5355339059327379e-01*g[2]*f[3]+3.5355339059327379e-01*g[7]*f[1]+3.5355339059327379e-01*f[6]*g[0]+3.5355339059327379e-01*f[7]*g[1]+3.5355339059327379e-01*g[6]*f[0];
  fg[7] =  3.5355339059327379e-01*f[4]*g[3]+3.5355339059327379e-01*f[5]*g[2]+3.5355339059327379e-01*g[7]*f[0]+3.5355339059327379e-01*f[6]*g[1]+3.5355339059327379e-01*g[5]*f[2]+3.5355339059327379e-01*g[4]*f[3]+3.5355339059327379e-01*f[7]*g[0]+3.5355339059327379e-01*g[6]*f[1];
  // nsum = 56, nprod = 128
}

struct gkyl_kern_op_count op_count_binop_mul_3d_ser_p1(void)
{
  return (struct gkyl_kern_op_count) { .num_sum = 56, .num_prod = 128 };
}