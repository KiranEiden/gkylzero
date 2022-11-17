#include <gkyl_euler_kernels.h> 
#include <gkyl_basis_ser_2x_p1_inv.h> 
GKYL_CU_DH void euler_pkpm_p_force_2x_ser_p1(const double *bvar, const double *div_p, const double *vlasov_pkpm_moms, 
  const double *p_perp, const double *div_b, 
  double* p_force) 
{ 
  // bvar:             magnetic field unit vector (nine components; first three components, b_i, other six components, b_i b_j.) 
  // div_p:            Volume expansion of divergence of pressure tensor.
  // vlasov_pkpm_moms: [rho, p_parallel, q_parallel], Moments computed from kinetic equation in pkpm model.
  // p_perp:           Perpendicular pressure.
  // div_b:            Volume expansion of divergence of magnetic field unit vector.
  // p_force:          total pressure force = 1/rho (b . div(P) + p_perp div(b)) for Euler PKPM.

  const double *bx = &bvar[0]; 
  const double *by = &bvar[4]; 
  const double *bz = &bvar[8]; 
  const double *div_p_x = &div_p[0]; 
  const double *div_p_y = &div_p[4]; 
  const double *div_p_z = &div_p[8]; 
  const double *rho = &vlasov_pkpm_moms[0]; 

  double b_div_p[4] = {0.0}; 
  double rho_inv[4] = {0.0}; 
  ser_2x_p1_inv(rho, rho_inv); 
  b_div_p[0] = (-0.5*div_b[3]*p_perp[3])+0.5*bz[3]*div_p_z[3]+0.5*by[3]*div_p_y[3]+0.5*bx[3]*div_p_x[3]-0.5*div_b[2]*p_perp[2]+0.5*bz[2]*div_p_z[2]+0.5*by[2]*div_p_y[2]+0.5*bx[2]*div_p_x[2]-0.5*div_b[1]*p_perp[1]+0.5*bz[1]*div_p_z[1]+0.5*by[1]*div_p_y[1]+0.5*bx[1]*div_p_x[1]-0.5*div_b[0]*p_perp[0]+0.5*bz[0]*div_p_z[0]+0.5*by[0]*div_p_y[0]+0.5*bx[0]*div_p_x[0]; 
  b_div_p[1] = (-0.5*div_b[2]*p_perp[3])+0.5*bz[2]*div_p_z[3]+0.5*by[2]*div_p_y[3]+0.5*bx[2]*div_p_x[3]-0.5*p_perp[2]*div_b[3]+0.5*div_p_z[2]*bz[3]+0.5*div_p_y[2]*by[3]+0.5*div_p_x[2]*bx[3]-0.5*div_b[0]*p_perp[1]+0.5*bz[0]*div_p_z[1]+0.5*by[0]*div_p_y[1]+0.5*bx[0]*div_p_x[1]-0.5*p_perp[0]*div_b[1]+0.5*div_p_z[0]*bz[1]+0.5*div_p_y[0]*by[1]+0.5*div_p_x[0]*bx[1]; 
  b_div_p[2] = (-0.5*div_b[1]*p_perp[3])+0.5*bz[1]*div_p_z[3]+0.5*by[1]*div_p_y[3]+0.5*bx[1]*div_p_x[3]-0.5*p_perp[1]*div_b[3]+0.5*div_p_z[1]*bz[3]+0.5*div_p_y[1]*by[3]+0.5*div_p_x[1]*bx[3]-0.5*div_b[0]*p_perp[2]+0.5*bz[0]*div_p_z[2]+0.5*by[0]*div_p_y[2]+0.5*bx[0]*div_p_x[2]-0.5*p_perp[0]*div_b[2]+0.5*div_p_z[0]*bz[2]+0.5*div_p_y[0]*by[2]+0.5*div_p_x[0]*bx[2]; 
  b_div_p[3] = (-0.5*div_b[0]*p_perp[3])+0.5*bz[0]*div_p_z[3]+0.5*by[0]*div_p_y[3]+0.5*bx[0]*div_p_x[3]-0.5*p_perp[0]*div_b[3]+0.5*div_p_z[0]*bz[3]+0.5*div_p_y[0]*by[3]+0.5*div_p_x[0]*bx[3]-0.5*div_b[1]*p_perp[2]+0.5*bz[1]*div_p_z[2]+0.5*by[1]*div_p_y[2]+0.5*bx[1]*div_p_x[2]-0.5*p_perp[1]*div_b[2]+0.5*div_p_z[1]*bz[2]+0.5*div_p_y[1]*by[2]+0.5*div_p_x[1]*bx[2]; 

  p_force[0] = 0.5*b_div_p[3]*rho_inv[3]+0.5*b_div_p[2]*rho_inv[2]+0.5*b_div_p[1]*rho_inv[1]+0.5*b_div_p[0]*rho_inv[0]; 
  p_force[1] = 0.5*b_div_p[2]*rho_inv[3]+0.5*rho_inv[2]*b_div_p[3]+0.5*b_div_p[0]*rho_inv[1]+0.5*rho_inv[0]*b_div_p[1]; 
  p_force[2] = 0.5*b_div_p[1]*rho_inv[3]+0.5*rho_inv[1]*b_div_p[3]+0.5*b_div_p[0]*rho_inv[2]+0.5*rho_inv[0]*b_div_p[2]; 
  p_force[3] = 0.5*b_div_p[0]*rho_inv[3]+0.5*rho_inv[0]*b_div_p[3]+0.5*b_div_p[1]*rho_inv[2]+0.5*rho_inv[1]*b_div_p[2]; 

} 
