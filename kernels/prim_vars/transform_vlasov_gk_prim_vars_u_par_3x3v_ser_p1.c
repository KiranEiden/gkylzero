#include <gkyl_dg_prim_vars_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_3x_p1_inv.h> 
GKYL_CU_DH void transform_vlasov_gk_prim_vars_u_par_3x3v_ser_p1(const double *b_i, const double *moms, double* upar) 
{ 
  // moms: Input moments. 
  // b_i:  Contravariant components of field-aligned unit vector. 
  // upar: upar = udrift . bhat. 
 
  const double *m0 = &moms[0]; 
  const double *m1x = &moms[8]; 
  const double *m1y = &moms[16]; 
  const double *m1z = &moms[24]; 
 
  const double *b_x = &b_i[0]; 
  const double *b_y = &b_i[8]; 
  const double *b_z = &b_i[16]; 
 
  double m0_inv[8] = {0.0}; 

  double ux[8] = {0.0}; 
  double uy[8] = {0.0}; 
  double uz[8] = {0.0}; 
  ser_3x_p1_inv(m0, m0_inv); 
  binop_mul_3d_ser_p1(m1x, m0_inv, ux); 
  binop_mul_3d_ser_p1(m1y, m0_inv, uy); 
  binop_mul_3d_ser_p1(m1z, m0_inv, uz); 

  upar[0] = 0.3535533905932737*b_z[7]*uz[7]+0.3535533905932737*b_y[7]*uy[7]+0.3535533905932737*b_x[7]*ux[7]+0.3535533905932737*b_z[6]*uz[6]+0.3535533905932737*b_y[6]*uy[6]+0.3535533905932737*b_x[6]*ux[6]+0.3535533905932737*b_z[5]*uz[5]+0.3535533905932737*b_y[5]*uy[5]+0.3535533905932737*b_x[5]*ux[5]+0.3535533905932737*b_z[4]*uz[4]+0.3535533905932737*b_y[4]*uy[4]+0.3535533905932737*b_x[4]*ux[4]+0.3535533905932737*b_z[3]*uz[3]+0.3535533905932737*b_y[3]*uy[3]+0.3535533905932737*b_x[3]*ux[3]+0.3535533905932737*b_z[2]*uz[2]+0.3535533905932737*b_y[2]*uy[2]+0.3535533905932737*b_x[2]*ux[2]+0.3535533905932737*b_z[1]*uz[1]+0.3535533905932737*b_y[1]*uy[1]+0.3535533905932737*b_x[1]*ux[1]+0.3535533905932737*b_z[0]*uz[0]+0.3535533905932737*b_y[0]*uy[0]+0.3535533905932737*b_x[0]*ux[0]; 
  upar[1] = 0.3535533905932737*b_z[6]*uz[7]+0.3535533905932737*b_y[6]*uy[7]+0.3535533905932737*b_x[6]*ux[7]+0.3535533905932737*uz[6]*b_z[7]+0.3535533905932737*uy[6]*b_y[7]+0.3535533905932737*ux[6]*b_x[7]+0.3535533905932737*b_z[3]*uz[5]+0.3535533905932737*b_y[3]*uy[5]+0.3535533905932737*b_x[3]*ux[5]+0.3535533905932737*uz[3]*b_z[5]+0.3535533905932737*uy[3]*b_y[5]+0.3535533905932737*ux[3]*b_x[5]+0.3535533905932737*b_z[2]*uz[4]+0.3535533905932737*b_y[2]*uy[4]+0.3535533905932737*b_x[2]*ux[4]+0.3535533905932737*uz[2]*b_z[4]+0.3535533905932737*uy[2]*b_y[4]+0.3535533905932737*ux[2]*b_x[4]+0.3535533905932737*b_z[0]*uz[1]+0.3535533905932737*b_y[0]*uy[1]+0.3535533905932737*b_x[0]*ux[1]+0.3535533905932737*uz[0]*b_z[1]+0.3535533905932737*uy[0]*b_y[1]+0.3535533905932737*ux[0]*b_x[1]; 
  upar[2] = 0.3535533905932737*b_z[5]*uz[7]+0.3535533905932737*b_y[5]*uy[7]+0.3535533905932737*b_x[5]*ux[7]+0.3535533905932737*uz[5]*b_z[7]+0.3535533905932737*uy[5]*b_y[7]+0.3535533905932737*ux[5]*b_x[7]+0.3535533905932737*b_z[3]*uz[6]+0.3535533905932737*b_y[3]*uy[6]+0.3535533905932737*b_x[3]*ux[6]+0.3535533905932737*uz[3]*b_z[6]+0.3535533905932737*uy[3]*b_y[6]+0.3535533905932737*ux[3]*b_x[6]+0.3535533905932737*b_z[1]*uz[4]+0.3535533905932737*b_y[1]*uy[4]+0.3535533905932737*b_x[1]*ux[4]+0.3535533905932737*uz[1]*b_z[4]+0.3535533905932737*uy[1]*b_y[4]+0.3535533905932737*ux[1]*b_x[4]+0.3535533905932737*b_z[0]*uz[2]+0.3535533905932737*b_y[0]*uy[2]+0.3535533905932737*b_x[0]*ux[2]+0.3535533905932737*uz[0]*b_z[2]+0.3535533905932737*uy[0]*b_y[2]+0.3535533905932737*ux[0]*b_x[2]; 
  upar[3] = 0.3535533905932737*b_z[4]*uz[7]+0.3535533905932737*b_y[4]*uy[7]+0.3535533905932737*b_x[4]*ux[7]+0.3535533905932737*uz[4]*b_z[7]+0.3535533905932737*uy[4]*b_y[7]+0.3535533905932737*ux[4]*b_x[7]+0.3535533905932737*b_z[2]*uz[6]+0.3535533905932737*b_y[2]*uy[6]+0.3535533905932737*b_x[2]*ux[6]+0.3535533905932737*uz[2]*b_z[6]+0.3535533905932737*uy[2]*b_y[6]+0.3535533905932737*ux[2]*b_x[6]+0.3535533905932737*b_z[1]*uz[5]+0.3535533905932737*b_y[1]*uy[5]+0.3535533905932737*b_x[1]*ux[5]+0.3535533905932737*uz[1]*b_z[5]+0.3535533905932737*uy[1]*b_y[5]+0.3535533905932737*ux[1]*b_x[5]+0.3535533905932737*b_z[0]*uz[3]+0.3535533905932737*b_y[0]*uy[3]+0.3535533905932737*b_x[0]*ux[3]+0.3535533905932737*uz[0]*b_z[3]+0.3535533905932737*uy[0]*b_y[3]+0.3535533905932737*ux[0]*b_x[3]; 
  upar[4] = 0.3535533905932737*b_z[3]*uz[7]+0.3535533905932737*b_y[3]*uy[7]+0.3535533905932737*b_x[3]*ux[7]+0.3535533905932737*uz[3]*b_z[7]+0.3535533905932737*uy[3]*b_y[7]+0.3535533905932737*ux[3]*b_x[7]+0.3535533905932737*b_z[5]*uz[6]+0.3535533905932737*b_y[5]*uy[6]+0.3535533905932737*b_x[5]*ux[6]+0.3535533905932737*uz[5]*b_z[6]+0.3535533905932737*uy[5]*b_y[6]+0.3535533905932737*ux[5]*b_x[6]+0.3535533905932737*b_z[0]*uz[4]+0.3535533905932737*b_y[0]*uy[4]+0.3535533905932737*b_x[0]*ux[4]+0.3535533905932737*uz[0]*b_z[4]+0.3535533905932737*uy[0]*b_y[4]+0.3535533905932737*ux[0]*b_x[4]+0.3535533905932737*b_z[1]*uz[2]+0.3535533905932737*b_y[1]*uy[2]+0.3535533905932737*b_x[1]*ux[2]+0.3535533905932737*uz[1]*b_z[2]+0.3535533905932737*uy[1]*b_y[2]+0.3535533905932737*ux[1]*b_x[2]; 
  upar[5] = 0.3535533905932737*b_z[2]*uz[7]+0.3535533905932737*b_y[2]*uy[7]+0.3535533905932737*b_x[2]*ux[7]+0.3535533905932737*uz[2]*b_z[7]+0.3535533905932737*uy[2]*b_y[7]+0.3535533905932737*ux[2]*b_x[7]+0.3535533905932737*b_z[4]*uz[6]+0.3535533905932737*b_y[4]*uy[6]+0.3535533905932737*b_x[4]*ux[6]+0.3535533905932737*uz[4]*b_z[6]+0.3535533905932737*uy[4]*b_y[6]+0.3535533905932737*ux[4]*b_x[6]+0.3535533905932737*b_z[0]*uz[5]+0.3535533905932737*b_y[0]*uy[5]+0.3535533905932737*b_x[0]*ux[5]+0.3535533905932737*uz[0]*b_z[5]+0.3535533905932737*uy[0]*b_y[5]+0.3535533905932737*ux[0]*b_x[5]+0.3535533905932737*b_z[1]*uz[3]+0.3535533905932737*b_y[1]*uy[3]+0.3535533905932737*b_x[1]*ux[3]+0.3535533905932737*uz[1]*b_z[3]+0.3535533905932737*uy[1]*b_y[3]+0.3535533905932737*ux[1]*b_x[3]; 
  upar[6] = 0.3535533905932737*b_z[1]*uz[7]+0.3535533905932737*b_y[1]*uy[7]+0.3535533905932737*b_x[1]*ux[7]+0.3535533905932737*uz[1]*b_z[7]+0.3535533905932737*uy[1]*b_y[7]+0.3535533905932737*ux[1]*b_x[7]+0.3535533905932737*b_z[0]*uz[6]+0.3535533905932737*b_y[0]*uy[6]+0.3535533905932737*b_x[0]*ux[6]+0.3535533905932737*uz[0]*b_z[6]+0.3535533905932737*uy[0]*b_y[6]+0.3535533905932737*ux[0]*b_x[6]+0.3535533905932737*b_z[4]*uz[5]+0.3535533905932737*b_y[4]*uy[5]+0.3535533905932737*b_x[4]*ux[5]+0.3535533905932737*uz[4]*b_z[5]+0.3535533905932737*uy[4]*b_y[5]+0.3535533905932737*ux[4]*b_x[5]+0.3535533905932737*b_z[2]*uz[3]+0.3535533905932737*b_y[2]*uy[3]+0.3535533905932737*b_x[2]*ux[3]+0.3535533905932737*uz[2]*b_z[3]+0.3535533905932737*uy[2]*b_y[3]+0.3535533905932737*ux[2]*b_x[3]; 
  upar[7] = 0.3535533905932737*b_z[0]*uz[7]+0.3535533905932737*b_y[0]*uy[7]+0.3535533905932737*b_x[0]*ux[7]+0.3535533905932737*uz[0]*b_z[7]+0.3535533905932737*uy[0]*b_y[7]+0.3535533905932737*ux[0]*b_x[7]+0.3535533905932737*b_z[1]*uz[6]+0.3535533905932737*b_y[1]*uy[6]+0.3535533905932737*b_x[1]*ux[6]+0.3535533905932737*uz[1]*b_z[6]+0.3535533905932737*uy[1]*b_y[6]+0.3535533905932737*ux[1]*b_x[6]+0.3535533905932737*b_z[2]*uz[5]+0.3535533905932737*b_y[2]*uy[5]+0.3535533905932737*b_x[2]*ux[5]+0.3535533905932737*uz[2]*b_z[5]+0.3535533905932737*uy[2]*b_y[5]+0.3535533905932737*ux[2]*b_x[5]+0.3535533905932737*b_z[3]*uz[4]+0.3535533905932737*b_y[3]*uy[4]+0.3535533905932737*b_x[3]*ux[4]+0.3535533905932737*uz[3]*b_z[4]+0.3535533905932737*uy[3]*b_y[4]+0.3535533905932737*ux[3]*b_x[4]; 

} 
 
