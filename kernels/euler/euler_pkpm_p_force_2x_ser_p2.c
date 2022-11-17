#include <gkyl_euler_kernels.h> 
#include <gkyl_basis_ser_2x_p2_inv.h> 
GKYL_CU_DH void euler_pkpm_p_force_2x_ser_p2(const double *bvar, const double *div_p, const double *vlasov_pkpm_moms, 
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
  const double *by = &bvar[8]; 
  const double *bz = &bvar[16]; 
  const double *div_p_x = &div_p[0]; 
  const double *div_p_y = &div_p[8]; 
  const double *div_p_z = &div_p[16]; 
  const double *rho = &vlasov_pkpm_moms[0]; 

  double b_div_p[8] = {0.0}; 
  double rho_inv[8] = {0.0}; 
  ser_2x_p2_inv(rho, rho_inv); 
  b_div_p[0] = (-0.5*div_b[7]*p_perp[7])+0.5*bz[7]*div_p_z[7]+0.5*by[7]*div_p_y[7]+0.5*bx[7]*div_p_x[7]-0.5*div_b[6]*p_perp[6]+0.5*bz[6]*div_p_z[6]+0.5*by[6]*div_p_y[6]+0.5*bx[6]*div_p_x[6]-0.5*div_b[5]*p_perp[5]+0.5*bz[5]*div_p_z[5]+0.5*by[5]*div_p_y[5]+0.5*bx[5]*div_p_x[5]-0.5*div_b[4]*p_perp[4]+0.5*bz[4]*div_p_z[4]+0.5*by[4]*div_p_y[4]+0.5*bx[4]*div_p_x[4]-0.5*div_b[3]*p_perp[3]+0.5*bz[3]*div_p_z[3]+0.5*by[3]*div_p_y[3]+0.5*bx[3]*div_p_x[3]-0.5*div_b[2]*p_perp[2]+0.5*bz[2]*div_p_z[2]+0.5*by[2]*div_p_y[2]+0.5*bx[2]*div_p_x[2]-0.5*div_b[1]*p_perp[1]+0.5*bz[1]*div_p_z[1]+0.5*by[1]*div_p_y[1]+0.5*bx[1]*div_p_x[1]-0.5*div_b[0]*p_perp[0]+0.5*bz[0]*div_p_z[0]+0.5*by[0]*div_p_y[0]+0.5*bx[0]*div_p_x[0]; 
  b_div_p[1] = (-0.5000000000000001*div_b[5]*p_perp[7])+0.5000000000000001*bz[5]*div_p_z[7]+0.5000000000000001*by[5]*div_p_y[7]+0.5000000000000001*bx[5]*div_p_x[7]-0.5000000000000001*p_perp[5]*div_b[7]+0.5000000000000001*div_p_z[5]*bz[7]+0.5000000000000001*div_p_y[5]*by[7]+0.5000000000000001*div_p_x[5]*bx[7]-0.447213595499958*div_b[3]*p_perp[6]+0.447213595499958*bz[3]*div_p_z[6]+0.447213595499958*by[3]*div_p_y[6]+0.447213595499958*bx[3]*div_p_x[6]-0.447213595499958*p_perp[3]*div_b[6]+0.447213595499958*div_p_z[3]*bz[6]+0.447213595499958*div_p_y[3]*by[6]+0.447213595499958*div_p_x[3]*bx[6]-0.4472135954999579*div_b[1]*p_perp[4]+0.4472135954999579*bz[1]*div_p_z[4]+0.4472135954999579*by[1]*div_p_y[4]+0.4472135954999579*bx[1]*div_p_x[4]-0.4472135954999579*p_perp[1]*div_b[4]+0.4472135954999579*div_p_z[1]*bz[4]+0.4472135954999579*div_p_y[1]*by[4]+0.4472135954999579*div_p_x[1]*bx[4]-0.5*div_b[2]*p_perp[3]+0.5*bz[2]*div_p_z[3]+0.5*by[2]*div_p_y[3]+0.5*bx[2]*div_p_x[3]-0.5*p_perp[2]*div_b[3]+0.5*div_p_z[2]*bz[3]+0.5*div_p_y[2]*by[3]+0.5*div_p_x[2]*bx[3]-0.5*div_b[0]*p_perp[1]+0.5*bz[0]*div_p_z[1]+0.5*by[0]*div_p_y[1]+0.5*bx[0]*div_p_x[1]-0.5*p_perp[0]*div_b[1]+0.5*div_p_z[0]*bz[1]+0.5*div_p_y[0]*by[1]+0.5*div_p_x[0]*bx[1]; 
  b_div_p[2] = (-0.447213595499958*div_b[3]*p_perp[7])+0.447213595499958*bz[3]*div_p_z[7]+0.447213595499958*by[3]*div_p_y[7]+0.447213595499958*bx[3]*div_p_x[7]-0.447213595499958*p_perp[3]*div_b[7]+0.447213595499958*div_p_z[3]*bz[7]+0.447213595499958*div_p_y[3]*by[7]+0.447213595499958*div_p_x[3]*bx[7]-0.5000000000000001*div_b[4]*p_perp[6]+0.5000000000000001*bz[4]*div_p_z[6]+0.5000000000000001*by[4]*div_p_y[6]+0.5000000000000001*bx[4]*div_p_x[6]-0.5000000000000001*p_perp[4]*div_b[6]+0.5000000000000001*div_p_z[4]*bz[6]+0.5000000000000001*div_p_y[4]*by[6]+0.5000000000000001*div_p_x[4]*bx[6]-0.4472135954999579*div_b[2]*p_perp[5]+0.4472135954999579*bz[2]*div_p_z[5]+0.4472135954999579*by[2]*div_p_y[5]+0.4472135954999579*bx[2]*div_p_x[5]-0.4472135954999579*p_perp[2]*div_b[5]+0.4472135954999579*div_p_z[2]*bz[5]+0.4472135954999579*div_p_y[2]*by[5]+0.4472135954999579*div_p_x[2]*bx[5]-0.5*div_b[1]*p_perp[3]+0.5*bz[1]*div_p_z[3]+0.5*by[1]*div_p_y[3]+0.5*bx[1]*div_p_x[3]-0.5*p_perp[1]*div_b[3]+0.5*div_p_z[1]*bz[3]+0.5*div_p_y[1]*by[3]+0.5*div_p_x[1]*bx[3]-0.5*div_b[0]*p_perp[2]+0.5*bz[0]*div_p_z[2]+0.5*by[0]*div_p_y[2]+0.5*bx[0]*div_p_x[2]-0.5*p_perp[0]*div_b[2]+0.5*div_p_z[0]*bz[2]+0.5*div_p_y[0]*by[2]+0.5*div_p_x[0]*bx[2]; 
  b_div_p[3] = (-0.4*div_b[6]*p_perp[7])-0.447213595499958*div_b[2]*p_perp[7]+0.4*bz[6]*div_p_z[7]+0.447213595499958*bz[2]*div_p_z[7]+0.4*by[6]*div_p_y[7]+0.447213595499958*by[2]*div_p_y[7]+0.4*bx[6]*div_p_x[7]+0.447213595499958*bx[2]*div_p_x[7]-0.4*p_perp[6]*div_b[7]-0.447213595499958*p_perp[2]*div_b[7]+0.4*div_p_z[6]*bz[7]+0.447213595499958*div_p_z[2]*bz[7]+0.4*div_p_y[6]*by[7]+0.447213595499958*div_p_y[2]*by[7]+0.4*div_p_x[6]*bx[7]+0.447213595499958*div_p_x[2]*bx[7]-0.447213595499958*div_b[1]*p_perp[6]+0.447213595499958*bz[1]*div_p_z[6]+0.447213595499958*by[1]*div_p_y[6]+0.447213595499958*bx[1]*div_p_x[6]-0.447213595499958*p_perp[1]*div_b[6]+0.447213595499958*div_p_z[1]*bz[6]+0.447213595499958*div_p_y[1]*by[6]+0.447213595499958*div_p_x[1]*bx[6]-0.4472135954999579*div_b[3]*p_perp[5]+0.4472135954999579*bz[3]*div_p_z[5]+0.4472135954999579*by[3]*div_p_y[5]+0.4472135954999579*bx[3]*div_p_x[5]-0.4472135954999579*p_perp[3]*div_b[5]+0.4472135954999579*div_p_z[3]*bz[5]+0.4472135954999579*div_p_y[3]*by[5]+0.4472135954999579*div_p_x[3]*bx[5]-0.4472135954999579*div_b[3]*p_perp[4]+0.4472135954999579*bz[3]*div_p_z[4]+0.4472135954999579*by[3]*div_p_y[4]+0.4472135954999579*bx[3]*div_p_x[4]-0.4472135954999579*p_perp[3]*div_b[4]+0.4472135954999579*div_p_z[3]*bz[4]+0.4472135954999579*div_p_y[3]*by[4]+0.4472135954999579*div_p_x[3]*bx[4]-0.5*div_b[0]*p_perp[3]+0.5*bz[0]*div_p_z[3]+0.5*by[0]*div_p_y[3]+0.5*bx[0]*div_p_x[3]-0.5*p_perp[0]*div_b[3]+0.5*div_p_z[0]*bz[3]+0.5*div_p_y[0]*by[3]+0.5*div_p_x[0]*bx[3]-0.5*div_b[1]*p_perp[2]+0.5*bz[1]*div_p_z[2]+0.5*by[1]*div_p_y[2]+0.5*bx[1]*div_p_x[2]-0.5*p_perp[1]*div_b[2]+0.5*div_p_z[1]*bz[2]+0.5*div_p_y[1]*by[2]+0.5*div_p_x[1]*bx[2]; 
  b_div_p[4] = (-0.4472135954999579*div_b[7]*p_perp[7])+0.4472135954999579*bz[7]*div_p_z[7]+0.4472135954999579*by[7]*div_p_y[7]+0.4472135954999579*bx[7]*div_p_x[7]-0.31943828249997*div_b[6]*p_perp[6]-0.5000000000000001*div_b[2]*p_perp[6]+0.31943828249997*bz[6]*div_p_z[6]+0.5000000000000001*bz[2]*div_p_z[6]+0.31943828249997*by[6]*div_p_y[6]+0.5000000000000001*by[2]*div_p_y[6]+0.31943828249997*bx[6]*div_p_x[6]+0.5000000000000001*bx[2]*div_p_x[6]-0.5000000000000001*p_perp[2]*div_b[6]+0.5000000000000001*div_p_z[2]*bz[6]+0.5000000000000001*div_p_y[2]*by[6]+0.5000000000000001*div_p_x[2]*bx[6]-0.31943828249997*div_b[4]*p_perp[4]-0.5*div_b[0]*p_perp[4]+0.31943828249997*bz[4]*div_p_z[4]+0.5*bz[0]*div_p_z[4]+0.31943828249997*by[4]*div_p_y[4]+0.5*by[0]*div_p_y[4]+0.31943828249997*bx[4]*div_p_x[4]+0.5*bx[0]*div_p_x[4]-0.5*p_perp[0]*div_b[4]+0.5*div_p_z[0]*bz[4]+0.5*div_p_y[0]*by[4]+0.5*div_p_x[0]*bx[4]-0.4472135954999579*div_b[3]*p_perp[3]+0.4472135954999579*bz[3]*div_p_z[3]+0.4472135954999579*by[3]*div_p_y[3]+0.4472135954999579*bx[3]*div_p_x[3]-0.4472135954999579*div_b[1]*p_perp[1]+0.4472135954999579*bz[1]*div_p_z[1]+0.4472135954999579*by[1]*div_p_y[1]+0.4472135954999579*bx[1]*div_p_x[1]; 
  b_div_p[5] = (-0.31943828249997*div_b[7]*p_perp[7])-0.5000000000000001*div_b[1]*p_perp[7]+0.31943828249997*bz[7]*div_p_z[7]+0.5000000000000001*bz[1]*div_p_z[7]+0.31943828249997*by[7]*div_p_y[7]+0.5000000000000001*by[1]*div_p_y[7]+0.31943828249997*bx[7]*div_p_x[7]+0.5000000000000001*bx[1]*div_p_x[7]-0.5000000000000001*p_perp[1]*div_b[7]+0.5000000000000001*div_p_z[1]*bz[7]+0.5000000000000001*div_p_y[1]*by[7]+0.5000000000000001*div_p_x[1]*bx[7]-0.4472135954999579*div_b[6]*p_perp[6]+0.4472135954999579*bz[6]*div_p_z[6]+0.4472135954999579*by[6]*div_p_y[6]+0.4472135954999579*bx[6]*div_p_x[6]-0.31943828249997*div_b[5]*p_perp[5]-0.5*div_b[0]*p_perp[5]+0.31943828249997*bz[5]*div_p_z[5]+0.5*bz[0]*div_p_z[5]+0.31943828249997*by[5]*div_p_y[5]+0.5*by[0]*div_p_y[5]+0.31943828249997*bx[5]*div_p_x[5]+0.5*bx[0]*div_p_x[5]-0.5*p_perp[0]*div_b[5]+0.5*div_p_z[0]*bz[5]+0.5*div_p_y[0]*by[5]+0.5*div_p_x[0]*bx[5]-0.4472135954999579*div_b[3]*p_perp[3]+0.4472135954999579*bz[3]*div_p_z[3]+0.4472135954999579*by[3]*div_p_y[3]+0.4472135954999579*bx[3]*div_p_x[3]-0.4472135954999579*div_b[2]*p_perp[2]+0.4472135954999579*bz[2]*div_p_z[2]+0.4472135954999579*by[2]*div_p_y[2]+0.4472135954999579*bx[2]*div_p_x[2]; 
  b_div_p[6] = (-0.4*div_b[3]*p_perp[7])+0.4*bz[3]*div_p_z[7]+0.4*by[3]*div_p_y[7]+0.4*bx[3]*div_p_x[7]-0.4*p_perp[3]*div_b[7]+0.4*div_p_z[3]*bz[7]+0.4*div_p_y[3]*by[7]+0.4*div_p_x[3]*bx[7]-0.4472135954999579*div_b[5]*p_perp[6]-0.31943828249997*div_b[4]*p_perp[6]-0.5*div_b[0]*p_perp[6]+0.4472135954999579*bz[5]*div_p_z[6]+0.31943828249997*bz[4]*div_p_z[6]+0.5*bz[0]*div_p_z[6]+0.4472135954999579*by[5]*div_p_y[6]+0.31943828249997*by[4]*div_p_y[6]+0.5*by[0]*div_p_y[6]+0.4472135954999579*bx[5]*div_p_x[6]+0.31943828249997*bx[4]*div_p_x[6]+0.5*bx[0]*div_p_x[6]-0.4472135954999579*p_perp[5]*div_b[6]-0.31943828249997*p_perp[4]*div_b[6]-0.5*p_perp[0]*div_b[6]+0.4472135954999579*div_p_z[5]*bz[6]+0.31943828249997*div_p_z[4]*bz[6]+0.5*div_p_z[0]*bz[6]+0.4472135954999579*div_p_y[5]*by[6]+0.31943828249997*div_p_y[4]*by[6]+0.5*div_p_y[0]*by[6]+0.4472135954999579*div_p_x[5]*bx[6]+0.31943828249997*div_p_x[4]*bx[6]+0.5*div_p_x[0]*bx[6]-0.5000000000000001*div_b[2]*p_perp[4]+0.5000000000000001*bz[2]*div_p_z[4]+0.5000000000000001*by[2]*div_p_y[4]+0.5000000000000001*bx[2]*div_p_x[4]-0.5000000000000001*p_perp[2]*div_b[4]+0.5000000000000001*div_p_z[2]*bz[4]+0.5000000000000001*div_p_y[2]*by[4]+0.5000000000000001*div_p_x[2]*bx[4]-0.447213595499958*div_b[1]*p_perp[3]+0.447213595499958*bz[1]*div_p_z[3]+0.447213595499958*by[1]*div_p_y[3]+0.447213595499958*bx[1]*div_p_x[3]-0.447213595499958*p_perp[1]*div_b[3]+0.447213595499958*div_p_z[1]*bz[3]+0.447213595499958*div_p_y[1]*by[3]+0.447213595499958*div_p_x[1]*bx[3]; 
  b_div_p[7] = (-0.31943828249997*div_b[5]*p_perp[7])-0.4472135954999579*div_b[4]*p_perp[7]-0.5*div_b[0]*p_perp[7]+0.31943828249997*bz[5]*div_p_z[7]+0.4472135954999579*bz[4]*div_p_z[7]+0.5*bz[0]*div_p_z[7]+0.31943828249997*by[5]*div_p_y[7]+0.4472135954999579*by[4]*div_p_y[7]+0.5*by[0]*div_p_y[7]+0.31943828249997*bx[5]*div_p_x[7]+0.4472135954999579*bx[4]*div_p_x[7]+0.5*bx[0]*div_p_x[7]-0.31943828249997*p_perp[5]*div_b[7]-0.4472135954999579*p_perp[4]*div_b[7]-0.5*p_perp[0]*div_b[7]+0.31943828249997*div_p_z[5]*bz[7]+0.4472135954999579*div_p_z[4]*bz[7]+0.5*div_p_z[0]*bz[7]+0.31943828249997*div_p_y[5]*by[7]+0.4472135954999579*div_p_y[4]*by[7]+0.5*div_p_y[0]*by[7]+0.31943828249997*div_p_x[5]*bx[7]+0.4472135954999579*div_p_x[4]*bx[7]+0.5*div_p_x[0]*bx[7]-0.4*div_b[3]*p_perp[6]+0.4*bz[3]*div_p_z[6]+0.4*by[3]*div_p_y[6]+0.4*bx[3]*div_p_x[6]-0.4*p_perp[3]*div_b[6]+0.4*div_p_z[3]*bz[6]+0.4*div_p_y[3]*by[6]+0.4*div_p_x[3]*bx[6]-0.5000000000000001*div_b[1]*p_perp[5]+0.5000000000000001*bz[1]*div_p_z[5]+0.5000000000000001*by[1]*div_p_y[5]+0.5000000000000001*bx[1]*div_p_x[5]-0.5000000000000001*p_perp[1]*div_b[5]+0.5000000000000001*div_p_z[1]*bz[5]+0.5000000000000001*div_p_y[1]*by[5]+0.5000000000000001*div_p_x[1]*bx[5]-0.447213595499958*div_b[2]*p_perp[3]+0.447213595499958*bz[2]*div_p_z[3]+0.447213595499958*by[2]*div_p_y[3]+0.447213595499958*bx[2]*div_p_x[3]-0.447213595499958*p_perp[2]*div_b[3]+0.447213595499958*div_p_z[2]*bz[3]+0.447213595499958*div_p_y[2]*by[3]+0.447213595499958*div_p_x[2]*bx[3]; 

  p_force[0] = 0.5*b_div_p[7]*rho_inv[7]+0.5*b_div_p[6]*rho_inv[6]+0.5*b_div_p[5]*rho_inv[5]+0.5*b_div_p[4]*rho_inv[4]+0.5*b_div_p[3]*rho_inv[3]+0.5*b_div_p[2]*rho_inv[2]+0.5*b_div_p[1]*rho_inv[1]+0.5*b_div_p[0]*rho_inv[0]; 
  p_force[1] = 0.5000000000000001*b_div_p[5]*rho_inv[7]+0.5000000000000001*rho_inv[5]*b_div_p[7]+0.447213595499958*b_div_p[3]*rho_inv[6]+0.447213595499958*rho_inv[3]*b_div_p[6]+0.4472135954999579*b_div_p[1]*rho_inv[4]+0.4472135954999579*rho_inv[1]*b_div_p[4]+0.5*b_div_p[2]*rho_inv[3]+0.5*rho_inv[2]*b_div_p[3]+0.5*b_div_p[0]*rho_inv[1]+0.5*rho_inv[0]*b_div_p[1]; 
  p_force[2] = 0.447213595499958*b_div_p[3]*rho_inv[7]+0.447213595499958*rho_inv[3]*b_div_p[7]+0.5000000000000001*b_div_p[4]*rho_inv[6]+0.5000000000000001*rho_inv[4]*b_div_p[6]+0.4472135954999579*b_div_p[2]*rho_inv[5]+0.4472135954999579*rho_inv[2]*b_div_p[5]+0.5*b_div_p[1]*rho_inv[3]+0.5*rho_inv[1]*b_div_p[3]+0.5*b_div_p[0]*rho_inv[2]+0.5*rho_inv[0]*b_div_p[2]; 
  p_force[3] = 0.4*b_div_p[6]*rho_inv[7]+0.447213595499958*b_div_p[2]*rho_inv[7]+0.4*rho_inv[6]*b_div_p[7]+0.447213595499958*rho_inv[2]*b_div_p[7]+0.447213595499958*b_div_p[1]*rho_inv[6]+0.447213595499958*rho_inv[1]*b_div_p[6]+0.4472135954999579*b_div_p[3]*rho_inv[5]+0.4472135954999579*rho_inv[3]*b_div_p[5]+0.4472135954999579*b_div_p[3]*rho_inv[4]+0.4472135954999579*rho_inv[3]*b_div_p[4]+0.5*b_div_p[0]*rho_inv[3]+0.5*rho_inv[0]*b_div_p[3]+0.5*b_div_p[1]*rho_inv[2]+0.5*rho_inv[1]*b_div_p[2]; 
  p_force[4] = 0.4472135954999579*b_div_p[7]*rho_inv[7]+0.31943828249997*b_div_p[6]*rho_inv[6]+0.5000000000000001*b_div_p[2]*rho_inv[6]+0.5000000000000001*rho_inv[2]*b_div_p[6]+0.31943828249997*b_div_p[4]*rho_inv[4]+0.5*b_div_p[0]*rho_inv[4]+0.5*rho_inv[0]*b_div_p[4]+0.4472135954999579*b_div_p[3]*rho_inv[3]+0.4472135954999579*b_div_p[1]*rho_inv[1]; 
  p_force[5] = 0.31943828249997*b_div_p[7]*rho_inv[7]+0.5000000000000001*b_div_p[1]*rho_inv[7]+0.5000000000000001*rho_inv[1]*b_div_p[7]+0.4472135954999579*b_div_p[6]*rho_inv[6]+0.31943828249997*b_div_p[5]*rho_inv[5]+0.5*b_div_p[0]*rho_inv[5]+0.5*rho_inv[0]*b_div_p[5]+0.4472135954999579*b_div_p[3]*rho_inv[3]+0.4472135954999579*b_div_p[2]*rho_inv[2]; 
  p_force[6] = 0.4*b_div_p[3]*rho_inv[7]+0.4*rho_inv[3]*b_div_p[7]+0.4472135954999579*b_div_p[5]*rho_inv[6]+0.31943828249997*b_div_p[4]*rho_inv[6]+0.5*b_div_p[0]*rho_inv[6]+0.4472135954999579*rho_inv[5]*b_div_p[6]+0.31943828249997*rho_inv[4]*b_div_p[6]+0.5*rho_inv[0]*b_div_p[6]+0.5000000000000001*b_div_p[2]*rho_inv[4]+0.5000000000000001*rho_inv[2]*b_div_p[4]+0.447213595499958*b_div_p[1]*rho_inv[3]+0.447213595499958*rho_inv[1]*b_div_p[3]; 
  p_force[7] = 0.31943828249997*b_div_p[5]*rho_inv[7]+0.4472135954999579*b_div_p[4]*rho_inv[7]+0.5*b_div_p[0]*rho_inv[7]+0.31943828249997*rho_inv[5]*b_div_p[7]+0.4472135954999579*rho_inv[4]*b_div_p[7]+0.5*rho_inv[0]*b_div_p[7]+0.4*b_div_p[3]*rho_inv[6]+0.4*rho_inv[3]*b_div_p[6]+0.5000000000000001*b_div_p[1]*rho_inv[5]+0.5000000000000001*rho_inv[1]*b_div_p[5]+0.447213595499958*b_div_p[2]*rho_inv[3]+0.447213595499958*rho_inv[2]*b_div_p[3]; 

} 
