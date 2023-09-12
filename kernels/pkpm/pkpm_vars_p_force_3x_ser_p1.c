#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH void pkpm_vars_p_force_3x_ser_p1(const double *prim_c, const double *div_b, 
  double* GKYL_RESTRICT pkpm_accel) 
{ 
  // prim_c:     Input volume expansion of primitive variables [ux, uy, uz, 1/rho div(p_par b), T_perp/m, m/T_perp] in center cell. 
  // div_b:      Input volume expansion of div(b) in center cell. 
  // pkpm_accel: Output volume expansion of pkpm acceleration variables. 

  const double *pkpm_div_ppar = &prim_c[24]; 
  const double *T_perp_over_m = &prim_c[32]; 

  double *p_perp_div_b = &pkpm_accel[0]; 
  double *p_force = &pkpm_accel[16]; 

  binop_mul_3d_ser_p1(T_perp_over_m, div_b, p_perp_div_b); 

  p_force[0] += pkpm_div_ppar[0]-1.0*p_perp_div_b[0]; 
  p_force[1] += pkpm_div_ppar[1]-1.0*p_perp_div_b[1]; 
  p_force[2] += pkpm_div_ppar[2]-1.0*p_perp_div_b[2]; 
  p_force[3] += pkpm_div_ppar[3]-1.0*p_perp_div_b[3]; 
  p_force[4] += pkpm_div_ppar[4]-1.0*p_perp_div_b[4]; 
  p_force[5] += pkpm_div_ppar[5]-1.0*p_perp_div_b[5]; 
  p_force[6] += pkpm_div_ppar[6]-1.0*p_perp_div_b[6]; 
  p_force[7] += pkpm_div_ppar[7]-1.0*p_perp_div_b[7]; 

} 
