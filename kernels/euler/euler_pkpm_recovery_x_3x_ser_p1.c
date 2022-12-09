#include <gkyl_euler_kernels.h> 
#include <gkyl_basis_ser_3x_p1_inv.h> 
GKYL_CU_DH void euler_pkpm_recovery_x_3x_ser_p1(const double *dxv, double nuHyp, 
  const double *bvarl, const double *bvarc, const double *bvarr, 
  const double *u_il, const double *u_ic, const double *u_ir, 
  const double *p_ijl, const double *p_ijc, const double *p_ijr, 
  const double *vlasov_pkpm_momsl, const double *vlasov_pkpm_momsc, const double *vlasov_pkpm_momsr, 
  const double *statevecl, const double *statevecc, const double *statevecr, 
  double* div_b, double* bb_grad_u, double* div_p, double* p_force, double* p_perp_source) 
{ 
  // dxv[NDIM]: Cell spacing.
  // nuHyp: Hyper-diffusion coefficient.
  // bvarl/bvarc/bvarr:  Input magnetic field unit vector in left/center/right cells.
  // u_il/u_ic/u_ir:     Input bulk velocity (ux,uy,uz) in left/center/right cells.
  // p_ijl/p_ijc/p_ijr:  Input pressure tensor in left/center/right cells.
  // vlasov_pkpm_momsl/vlasov_pkpm_momsc/vlasov_pkpm_momsr: Input pkpm moments (rho, p_parallel, q_parallel) in left/center/right cells.
  // statevecl/statevecc/statevecr:          [rho ux, rho uy, rho uz, p_perp], Fluid input state vector in center cell.
  // div_b:              Increment to volume expansion of div(b) in one direction.
  // bb_grad_u:          Increment to volume expansion of bb : grad(u) in one direction.
  // div_p:              Increment to volume expansion of div(p) in one direction.
  // p_force:            Increment to volume expansion of p_force = 1/rho * div(p_parallel b_hat) in one direction.
  // p_perp_source:      Increment to volume expansion of perpendicular pressure compression source (p_perp div(u) - p_perp bb : grad(u)).

  const double dx1 = 2.0/dxv[0]; 
  const double dx14 = dx1*dx1*dx1*dx1; 

  const double *b_l = &bvarl[0]; 
  const double *b_c = &bvarc[0]; 
  const double *b_r = &bvarr[0]; 

  const double *ux_l = &u_il[0]; 
  const double *uy_l = &u_il[8]; 
  const double *uz_l = &u_il[16]; 

  const double *ux_c = &u_ic[0]; 
  const double *uy_c = &u_ic[8]; 
  const double *uz_c = &u_ic[16]; 

  const double *ux_r = &u_ir[0]; 
  const double *uy_r = &u_ir[8]; 
  const double *uz_r = &u_ir[16]; 

  const double *bxbx = &bvarc[24]; 
  const double *bxby = &bvarc[32]; 
  const double *bxbz = &bvarc[40]; 
  const double *byby = &bvarc[48]; 
  const double *bybz = &bvarc[56]; 
  const double *bzbz = &bvarc[64]; 

  const double *Pxx_l = &p_ijl[0]; 
  const double *Pxy_l = &p_ijl[8]; 
  const double *Pxz_l = &p_ijl[16]; 
  const double *Pyy_l = &p_ijl[24]; 
  const double *Pyz_l = &p_ijl[32]; 
  const double *Pzz_l = &p_ijl[40]; 

  const double *Pxx_c = &p_ijc[0]; 
  const double *Pxy_c = &p_ijc[8]; 
  const double *Pxz_c = &p_ijc[16]; 
  const double *Pyy_c = &p_ijc[24]; 
  const double *Pyz_c = &p_ijc[32]; 
  const double *Pzz_c = &p_ijc[40]; 

  const double *Pxx_r = &p_ijr[0]; 
  const double *Pxy_r = &p_ijr[8]; 
  const double *Pxz_r = &p_ijr[16]; 
  const double *Pyy_r = &p_ijr[24]; 
  const double *Pyz_r = &p_ijr[32]; 
  const double *Pzz_r = &p_ijr[40]; 

  const double *ppar_l = &vlasov_pkpm_momsl[8]; 
  const double *ppar_c = &vlasov_pkpm_momsc[8]; 
  const double *ppar_r = &vlasov_pkpm_momsr[8]; 
  const double *rho = &vlasov_pkpm_momsc[0]; 

  const double *rhoux_l = &statevecl[0]; 
  const double *rhouy_l = &statevecl[8]; 
  const double *rhouz_l = &statevecl[16]; 
  const double *p_perp_l = &statevecl[24]; 

  const double *rhoux_c = &statevecc[0]; 
  const double *rhouy_c = &statevecc[8]; 
  const double *rhouz_c = &statevecc[16]; 
  const double *p_perp_c = &statevecc[24]; 

  const double *rhoux_r = &statevecr[0]; 
  const double *rhouy_r = &statevecr[8]; 
  const double *rhouz_r = &statevecr[16]; 
  const double *p_perp_r = &statevecr[24]; 

  const double *p_perp = &statevecc[24]; 

  double *div_p_x = &div_p[0]; 
  double *div_p_y = &div_p[8]; 
  double *div_p_z = &div_p[16]; 

  double grad_u_x[8] = {0.0}; 
  double grad_u_y[8] = {0.0}; 
  double grad_u_z[8] = {0.0}; 
  grad_u_x[0] = (-0.2886751345948129*ux_r[1])-0.2886751345948129*ux_l[1]+0.5773502691896258*ux_c[1]+0.25*ux_r[0]-0.25*ux_l[0]; 
  grad_u_x[1] = (-0.5*ux_r[1])+0.5*ux_l[1]+0.4330127018922193*ux_r[0]+0.4330127018922193*ux_l[0]-0.8660254037844386*ux_c[0]; 
  grad_u_x[2] = (-0.2886751345948129*ux_r[4])-0.2886751345948129*ux_l[4]+0.5773502691896258*ux_c[4]+0.25*ux_r[2]-0.25*ux_l[2]; 
  grad_u_x[3] = (-0.2886751345948129*ux_r[5])-0.2886751345948129*ux_l[5]+0.5773502691896258*ux_c[5]+0.25*ux_r[3]-0.25*ux_l[3]; 
  grad_u_x[4] = (-0.5*ux_r[4])+0.5*ux_l[4]+0.4330127018922193*ux_r[2]+0.4330127018922193*ux_l[2]-0.8660254037844386*ux_c[2]; 
  grad_u_x[5] = (-0.5*ux_r[5])+0.5*ux_l[5]+0.4330127018922193*ux_r[3]+0.4330127018922193*ux_l[3]-0.8660254037844386*ux_c[3]; 
  grad_u_x[6] = (-0.2886751345948129*ux_r[7])-0.2886751345948129*ux_l[7]+0.5773502691896258*ux_c[7]+0.25*ux_r[6]-0.25*ux_l[6]; 
  grad_u_x[7] = (-0.5*ux_r[7])+0.5*ux_l[7]+0.4330127018922193*ux_r[6]+0.4330127018922193*ux_l[6]-0.8660254037844386*ux_c[6]; 

  grad_u_y[0] = (-0.2886751345948129*uy_r[1])-0.2886751345948129*uy_l[1]+0.5773502691896258*uy_c[1]+0.25*uy_r[0]-0.25*uy_l[0]; 
  grad_u_y[1] = (-0.5*uy_r[1])+0.5*uy_l[1]+0.4330127018922193*uy_r[0]+0.4330127018922193*uy_l[0]-0.8660254037844386*uy_c[0]; 
  grad_u_y[2] = (-0.2886751345948129*uy_r[4])-0.2886751345948129*uy_l[4]+0.5773502691896258*uy_c[4]+0.25*uy_r[2]-0.25*uy_l[2]; 
  grad_u_y[3] = (-0.2886751345948129*uy_r[5])-0.2886751345948129*uy_l[5]+0.5773502691896258*uy_c[5]+0.25*uy_r[3]-0.25*uy_l[3]; 
  grad_u_y[4] = (-0.5*uy_r[4])+0.5*uy_l[4]+0.4330127018922193*uy_r[2]+0.4330127018922193*uy_l[2]-0.8660254037844386*uy_c[2]; 
  grad_u_y[5] = (-0.5*uy_r[5])+0.5*uy_l[5]+0.4330127018922193*uy_r[3]+0.4330127018922193*uy_l[3]-0.8660254037844386*uy_c[3]; 
  grad_u_y[6] = (-0.2886751345948129*uy_r[7])-0.2886751345948129*uy_l[7]+0.5773502691896258*uy_c[7]+0.25*uy_r[6]-0.25*uy_l[6]; 
  grad_u_y[7] = (-0.5*uy_r[7])+0.5*uy_l[7]+0.4330127018922193*uy_r[6]+0.4330127018922193*uy_l[6]-0.8660254037844386*uy_c[6]; 

  grad_u_z[0] = (-0.2886751345948129*uz_r[1])-0.2886751345948129*uz_l[1]+0.5773502691896258*uz_c[1]+0.25*uz_r[0]-0.25*uz_l[0]; 
  grad_u_z[1] = (-0.5*uz_r[1])+0.5*uz_l[1]+0.4330127018922193*uz_r[0]+0.4330127018922193*uz_l[0]-0.8660254037844386*uz_c[0]; 
  grad_u_z[2] = (-0.2886751345948129*uz_r[4])-0.2886751345948129*uz_l[4]+0.5773502691896258*uz_c[4]+0.25*uz_r[2]-0.25*uz_l[2]; 
  grad_u_z[3] = (-0.2886751345948129*uz_r[5])-0.2886751345948129*uz_l[5]+0.5773502691896258*uz_c[5]+0.25*uz_r[3]-0.25*uz_l[3]; 
  grad_u_z[4] = (-0.5*uz_r[4])+0.5*uz_l[4]+0.4330127018922193*uz_r[2]+0.4330127018922193*uz_l[2]-0.8660254037844386*uz_c[2]; 
  grad_u_z[5] = (-0.5*uz_r[5])+0.5*uz_l[5]+0.4330127018922193*uz_r[3]+0.4330127018922193*uz_l[3]-0.8660254037844386*uz_c[3]; 
  grad_u_z[6] = (-0.2886751345948129*uz_r[7])-0.2886751345948129*uz_l[7]+0.5773502691896258*uz_c[7]+0.25*uz_r[6]-0.25*uz_l[6]; 
  grad_u_z[7] = (-0.5*uz_r[7])+0.5*uz_l[7]+0.4330127018922193*uz_r[6]+0.4330127018922193*uz_l[6]-0.8660254037844386*uz_c[6]; 

  double ppar_b_l[8] = {0.0}; 
  double ppar_b_c[8] = {0.0}; 
  double ppar_b_r[8] = {0.0}; 
  ppar_b_l[0] = 0.3535533905932737*b_l[7]*ppar_l[7]+0.3535533905932737*b_l[6]*ppar_l[6]+0.3535533905932737*b_l[5]*ppar_l[5]+0.3535533905932737*b_l[4]*ppar_l[4]+0.3535533905932737*b_l[3]*ppar_l[3]+0.3535533905932737*b_l[2]*ppar_l[2]+0.3535533905932737*b_l[1]*ppar_l[1]+0.3535533905932737*b_l[0]*ppar_l[0]; 
  ppar_b_l[1] = 0.3535533905932737*b_l[6]*ppar_l[7]+0.3535533905932737*ppar_l[6]*b_l[7]+0.3535533905932737*b_l[3]*ppar_l[5]+0.3535533905932737*ppar_l[3]*b_l[5]+0.3535533905932737*b_l[2]*ppar_l[4]+0.3535533905932737*ppar_l[2]*b_l[4]+0.3535533905932737*b_l[0]*ppar_l[1]+0.3535533905932737*ppar_l[0]*b_l[1]; 
  ppar_b_l[2] = 0.3535533905932737*b_l[5]*ppar_l[7]+0.3535533905932737*ppar_l[5]*b_l[7]+0.3535533905932737*b_l[3]*ppar_l[6]+0.3535533905932737*ppar_l[3]*b_l[6]+0.3535533905932737*b_l[1]*ppar_l[4]+0.3535533905932737*ppar_l[1]*b_l[4]+0.3535533905932737*b_l[0]*ppar_l[2]+0.3535533905932737*ppar_l[0]*b_l[2]; 
  ppar_b_l[3] = 0.3535533905932737*b_l[4]*ppar_l[7]+0.3535533905932737*ppar_l[4]*b_l[7]+0.3535533905932737*b_l[2]*ppar_l[6]+0.3535533905932737*ppar_l[2]*b_l[6]+0.3535533905932737*b_l[1]*ppar_l[5]+0.3535533905932737*ppar_l[1]*b_l[5]+0.3535533905932737*b_l[0]*ppar_l[3]+0.3535533905932737*ppar_l[0]*b_l[3]; 
  ppar_b_l[4] = 0.3535533905932737*b_l[3]*ppar_l[7]+0.3535533905932737*ppar_l[3]*b_l[7]+0.3535533905932737*b_l[5]*ppar_l[6]+0.3535533905932737*ppar_l[5]*b_l[6]+0.3535533905932737*b_l[0]*ppar_l[4]+0.3535533905932737*ppar_l[0]*b_l[4]+0.3535533905932737*b_l[1]*ppar_l[2]+0.3535533905932737*ppar_l[1]*b_l[2]; 
  ppar_b_l[5] = 0.3535533905932737*b_l[2]*ppar_l[7]+0.3535533905932737*ppar_l[2]*b_l[7]+0.3535533905932737*b_l[4]*ppar_l[6]+0.3535533905932737*ppar_l[4]*b_l[6]+0.3535533905932737*b_l[0]*ppar_l[5]+0.3535533905932737*ppar_l[0]*b_l[5]+0.3535533905932737*b_l[1]*ppar_l[3]+0.3535533905932737*ppar_l[1]*b_l[3]; 
  ppar_b_l[6] = 0.3535533905932737*b_l[1]*ppar_l[7]+0.3535533905932737*ppar_l[1]*b_l[7]+0.3535533905932737*b_l[0]*ppar_l[6]+0.3535533905932737*ppar_l[0]*b_l[6]+0.3535533905932737*b_l[4]*ppar_l[5]+0.3535533905932737*ppar_l[4]*b_l[5]+0.3535533905932737*b_l[2]*ppar_l[3]+0.3535533905932737*ppar_l[2]*b_l[3]; 
  ppar_b_l[7] = 0.3535533905932737*b_l[0]*ppar_l[7]+0.3535533905932737*ppar_l[0]*b_l[7]+0.3535533905932737*b_l[1]*ppar_l[6]+0.3535533905932737*ppar_l[1]*b_l[6]+0.3535533905932737*b_l[2]*ppar_l[5]+0.3535533905932737*ppar_l[2]*b_l[5]+0.3535533905932737*b_l[3]*ppar_l[4]+0.3535533905932737*ppar_l[3]*b_l[4]; 

  ppar_b_c[0] = 0.3535533905932737*b_c[7]*ppar_c[7]+0.3535533905932737*b_c[6]*ppar_c[6]+0.3535533905932737*b_c[5]*ppar_c[5]+0.3535533905932737*b_c[4]*ppar_c[4]+0.3535533905932737*b_c[3]*ppar_c[3]+0.3535533905932737*b_c[2]*ppar_c[2]+0.3535533905932737*b_c[1]*ppar_c[1]+0.3535533905932737*b_c[0]*ppar_c[0]; 
  ppar_b_c[1] = 0.3535533905932737*b_c[6]*ppar_c[7]+0.3535533905932737*ppar_c[6]*b_c[7]+0.3535533905932737*b_c[3]*ppar_c[5]+0.3535533905932737*ppar_c[3]*b_c[5]+0.3535533905932737*b_c[2]*ppar_c[4]+0.3535533905932737*ppar_c[2]*b_c[4]+0.3535533905932737*b_c[0]*ppar_c[1]+0.3535533905932737*ppar_c[0]*b_c[1]; 
  ppar_b_c[2] = 0.3535533905932737*b_c[5]*ppar_c[7]+0.3535533905932737*ppar_c[5]*b_c[7]+0.3535533905932737*b_c[3]*ppar_c[6]+0.3535533905932737*ppar_c[3]*b_c[6]+0.3535533905932737*b_c[1]*ppar_c[4]+0.3535533905932737*ppar_c[1]*b_c[4]+0.3535533905932737*b_c[0]*ppar_c[2]+0.3535533905932737*ppar_c[0]*b_c[2]; 
  ppar_b_c[3] = 0.3535533905932737*b_c[4]*ppar_c[7]+0.3535533905932737*ppar_c[4]*b_c[7]+0.3535533905932737*b_c[2]*ppar_c[6]+0.3535533905932737*ppar_c[2]*b_c[6]+0.3535533905932737*b_c[1]*ppar_c[5]+0.3535533905932737*ppar_c[1]*b_c[5]+0.3535533905932737*b_c[0]*ppar_c[3]+0.3535533905932737*ppar_c[0]*b_c[3]; 
  ppar_b_c[4] = 0.3535533905932737*b_c[3]*ppar_c[7]+0.3535533905932737*ppar_c[3]*b_c[7]+0.3535533905932737*b_c[5]*ppar_c[6]+0.3535533905932737*ppar_c[5]*b_c[6]+0.3535533905932737*b_c[0]*ppar_c[4]+0.3535533905932737*ppar_c[0]*b_c[4]+0.3535533905932737*b_c[1]*ppar_c[2]+0.3535533905932737*ppar_c[1]*b_c[2]; 
  ppar_b_c[5] = 0.3535533905932737*b_c[2]*ppar_c[7]+0.3535533905932737*ppar_c[2]*b_c[7]+0.3535533905932737*b_c[4]*ppar_c[6]+0.3535533905932737*ppar_c[4]*b_c[6]+0.3535533905932737*b_c[0]*ppar_c[5]+0.3535533905932737*ppar_c[0]*b_c[5]+0.3535533905932737*b_c[1]*ppar_c[3]+0.3535533905932737*ppar_c[1]*b_c[3]; 
  ppar_b_c[6] = 0.3535533905932737*b_c[1]*ppar_c[7]+0.3535533905932737*ppar_c[1]*b_c[7]+0.3535533905932737*b_c[0]*ppar_c[6]+0.3535533905932737*ppar_c[0]*b_c[6]+0.3535533905932737*b_c[4]*ppar_c[5]+0.3535533905932737*ppar_c[4]*b_c[5]+0.3535533905932737*b_c[2]*ppar_c[3]+0.3535533905932737*ppar_c[2]*b_c[3]; 
  ppar_b_c[7] = 0.3535533905932737*b_c[0]*ppar_c[7]+0.3535533905932737*ppar_c[0]*b_c[7]+0.3535533905932737*b_c[1]*ppar_c[6]+0.3535533905932737*ppar_c[1]*b_c[6]+0.3535533905932737*b_c[2]*ppar_c[5]+0.3535533905932737*ppar_c[2]*b_c[5]+0.3535533905932737*b_c[3]*ppar_c[4]+0.3535533905932737*ppar_c[3]*b_c[4]; 

  ppar_b_r[0] = 0.3535533905932737*b_r[7]*ppar_r[7]+0.3535533905932737*b_r[6]*ppar_r[6]+0.3535533905932737*b_r[5]*ppar_r[5]+0.3535533905932737*b_r[4]*ppar_r[4]+0.3535533905932737*b_r[3]*ppar_r[3]+0.3535533905932737*b_r[2]*ppar_r[2]+0.3535533905932737*b_r[1]*ppar_r[1]+0.3535533905932737*b_r[0]*ppar_r[0]; 
  ppar_b_r[1] = 0.3535533905932737*b_r[6]*ppar_r[7]+0.3535533905932737*ppar_r[6]*b_r[7]+0.3535533905932737*b_r[3]*ppar_r[5]+0.3535533905932737*ppar_r[3]*b_r[5]+0.3535533905932737*b_r[2]*ppar_r[4]+0.3535533905932737*ppar_r[2]*b_r[4]+0.3535533905932737*b_r[0]*ppar_r[1]+0.3535533905932737*ppar_r[0]*b_r[1]; 
  ppar_b_r[2] = 0.3535533905932737*b_r[5]*ppar_r[7]+0.3535533905932737*ppar_r[5]*b_r[7]+0.3535533905932737*b_r[3]*ppar_r[6]+0.3535533905932737*ppar_r[3]*b_r[6]+0.3535533905932737*b_r[1]*ppar_r[4]+0.3535533905932737*ppar_r[1]*b_r[4]+0.3535533905932737*b_r[0]*ppar_r[2]+0.3535533905932737*ppar_r[0]*b_r[2]; 
  ppar_b_r[3] = 0.3535533905932737*b_r[4]*ppar_r[7]+0.3535533905932737*ppar_r[4]*b_r[7]+0.3535533905932737*b_r[2]*ppar_r[6]+0.3535533905932737*ppar_r[2]*b_r[6]+0.3535533905932737*b_r[1]*ppar_r[5]+0.3535533905932737*ppar_r[1]*b_r[5]+0.3535533905932737*b_r[0]*ppar_r[3]+0.3535533905932737*ppar_r[0]*b_r[3]; 
  ppar_b_r[4] = 0.3535533905932737*b_r[3]*ppar_r[7]+0.3535533905932737*ppar_r[3]*b_r[7]+0.3535533905932737*b_r[5]*ppar_r[6]+0.3535533905932737*ppar_r[5]*b_r[6]+0.3535533905932737*b_r[0]*ppar_r[4]+0.3535533905932737*ppar_r[0]*b_r[4]+0.3535533905932737*b_r[1]*ppar_r[2]+0.3535533905932737*ppar_r[1]*b_r[2]; 
  ppar_b_r[5] = 0.3535533905932737*b_r[2]*ppar_r[7]+0.3535533905932737*ppar_r[2]*b_r[7]+0.3535533905932737*b_r[4]*ppar_r[6]+0.3535533905932737*ppar_r[4]*b_r[6]+0.3535533905932737*b_r[0]*ppar_r[5]+0.3535533905932737*ppar_r[0]*b_r[5]+0.3535533905932737*b_r[1]*ppar_r[3]+0.3535533905932737*ppar_r[1]*b_r[3]; 
  ppar_b_r[6] = 0.3535533905932737*b_r[1]*ppar_r[7]+0.3535533905932737*ppar_r[1]*b_r[7]+0.3535533905932737*b_r[0]*ppar_r[6]+0.3535533905932737*ppar_r[0]*b_r[6]+0.3535533905932737*b_r[4]*ppar_r[5]+0.3535533905932737*ppar_r[4]*b_r[5]+0.3535533905932737*b_r[2]*ppar_r[3]+0.3535533905932737*ppar_r[2]*b_r[3]; 
  ppar_b_r[7] = 0.3535533905932737*b_r[0]*ppar_r[7]+0.3535533905932737*ppar_r[0]*b_r[7]+0.3535533905932737*b_r[1]*ppar_r[6]+0.3535533905932737*ppar_r[1]*b_r[6]+0.3535533905932737*b_r[2]*ppar_r[5]+0.3535533905932737*ppar_r[2]*b_r[5]+0.3535533905932737*b_r[3]*ppar_r[4]+0.3535533905932737*ppar_r[3]*b_r[4]; 

  div_b[0] += ((-0.2886751345948129*(b_r[1]+b_l[1]))+0.5773502691896258*b_c[1]+0.25*b_r[0]-0.25*b_l[0])*dx1; 
  div_b[1] += ((-0.5*b_r[1])+0.5*b_l[1]+0.4330127018922193*(b_r[0]+b_l[0])-0.8660254037844386*b_c[0])*dx1; 
  div_b[2] += ((-0.2886751345948129*(b_r[4]+b_l[4]))+0.5773502691896258*b_c[4]+0.25*b_r[2]-0.25*b_l[2])*dx1; 
  div_b[3] += ((-0.2886751345948129*(b_r[5]+b_l[5]))+0.5773502691896258*b_c[5]+0.25*b_r[3]-0.25*b_l[3])*dx1; 
  div_b[4] += ((-0.5*b_r[4])+0.5*b_l[4]+0.4330127018922193*(b_r[2]+b_l[2])-0.8660254037844386*b_c[2])*dx1; 
  div_b[5] += ((-0.5*b_r[5])+0.5*b_l[5]+0.4330127018922193*(b_r[3]+b_l[3])-0.8660254037844386*b_c[3])*dx1; 
  div_b[6] += ((-0.2886751345948129*(b_r[7]+b_l[7]))+0.5773502691896258*b_c[7]+0.25*b_r[6]-0.25*b_l[6])*dx1; 
  div_b[7] += ((-0.5*b_r[7])+0.5*b_l[7]+0.4330127018922193*(b_r[6]+b_l[6])-0.8660254037844386*b_c[6])*dx1; 

  bb_grad_u[0] += 0.3535533905932737*(bxbz[7]*grad_u_z[7]+bxby[7]*grad_u_y[7]+bxbx[7]*grad_u_x[7]+bxbz[6]*grad_u_z[6]+bxby[6]*grad_u_y[6]+bxbx[6]*grad_u_x[6]+bxbz[5]*grad_u_z[5]+bxby[5]*grad_u_y[5]+bxbx[5]*grad_u_x[5]+bxbz[4]*grad_u_z[4]+bxby[4]*grad_u_y[4]+bxbx[4]*grad_u_x[4]+bxbz[3]*grad_u_z[3]+bxby[3]*grad_u_y[3]+bxbx[3]*grad_u_x[3]+bxbz[2]*grad_u_z[2]+bxby[2]*grad_u_y[2]+bxbx[2]*grad_u_x[2]+bxbz[1]*grad_u_z[1]+bxby[1]*grad_u_y[1]+bxbx[1]*grad_u_x[1]+bxbz[0]*grad_u_z[0]+bxby[0]*grad_u_y[0]+bxbx[0]*grad_u_x[0])*dx1; 
  bb_grad_u[1] += 0.3535533905932737*(bxbz[6]*grad_u_z[7]+bxby[6]*grad_u_y[7]+bxbx[6]*grad_u_x[7]+grad_u_z[6]*bxbz[7]+grad_u_y[6]*bxby[7]+grad_u_x[6]*bxbx[7]+bxbz[3]*grad_u_z[5]+bxby[3]*grad_u_y[5]+bxbx[3]*grad_u_x[5]+grad_u_z[3]*bxbz[5]+grad_u_y[3]*bxby[5]+grad_u_x[3]*bxbx[5]+bxbz[2]*grad_u_z[4]+bxby[2]*grad_u_y[4]+bxbx[2]*grad_u_x[4]+grad_u_z[2]*bxbz[4]+grad_u_y[2]*bxby[4]+grad_u_x[2]*bxbx[4]+bxbz[0]*grad_u_z[1]+bxby[0]*grad_u_y[1]+bxbx[0]*grad_u_x[1]+grad_u_z[0]*bxbz[1]+grad_u_y[0]*bxby[1]+grad_u_x[0]*bxbx[1])*dx1; 
  bb_grad_u[2] += 0.3535533905932737*(bxbz[5]*grad_u_z[7]+bxby[5]*grad_u_y[7]+bxbx[5]*grad_u_x[7]+grad_u_z[5]*bxbz[7]+grad_u_y[5]*bxby[7]+grad_u_x[5]*bxbx[7]+bxbz[3]*grad_u_z[6]+bxby[3]*grad_u_y[6]+bxbx[3]*grad_u_x[6]+grad_u_z[3]*bxbz[6]+grad_u_y[3]*bxby[6]+grad_u_x[3]*bxbx[6]+bxbz[1]*grad_u_z[4]+bxby[1]*grad_u_y[4]+bxbx[1]*grad_u_x[4]+grad_u_z[1]*bxbz[4]+grad_u_y[1]*bxby[4]+grad_u_x[1]*bxbx[4]+bxbz[0]*grad_u_z[2]+bxby[0]*grad_u_y[2]+bxbx[0]*grad_u_x[2]+grad_u_z[0]*bxbz[2]+grad_u_y[0]*bxby[2]+grad_u_x[0]*bxbx[2])*dx1; 
  bb_grad_u[3] += 0.3535533905932737*(bxbz[4]*grad_u_z[7]+bxby[4]*grad_u_y[7]+bxbx[4]*grad_u_x[7]+grad_u_z[4]*bxbz[7]+grad_u_y[4]*bxby[7]+grad_u_x[4]*bxbx[7]+bxbz[2]*grad_u_z[6]+bxby[2]*grad_u_y[6]+bxbx[2]*grad_u_x[6]+grad_u_z[2]*bxbz[6]+grad_u_y[2]*bxby[6]+grad_u_x[2]*bxbx[6]+bxbz[1]*grad_u_z[5]+bxby[1]*grad_u_y[5]+bxbx[1]*grad_u_x[5]+grad_u_z[1]*bxbz[5]+grad_u_y[1]*bxby[5]+grad_u_x[1]*bxbx[5]+bxbz[0]*grad_u_z[3]+bxby[0]*grad_u_y[3]+bxbx[0]*grad_u_x[3]+grad_u_z[0]*bxbz[3]+grad_u_y[0]*bxby[3]+grad_u_x[0]*bxbx[3])*dx1; 
  bb_grad_u[4] += 0.3535533905932737*(bxbz[3]*grad_u_z[7]+bxby[3]*grad_u_y[7]+bxbx[3]*grad_u_x[7]+grad_u_z[3]*bxbz[7]+grad_u_y[3]*bxby[7]+grad_u_x[3]*bxbx[7]+bxbz[5]*grad_u_z[6]+bxby[5]*grad_u_y[6]+bxbx[5]*grad_u_x[6]+grad_u_z[5]*bxbz[6]+grad_u_y[5]*bxby[6]+grad_u_x[5]*bxbx[6]+bxbz[0]*grad_u_z[4]+bxby[0]*grad_u_y[4]+bxbx[0]*grad_u_x[4]+grad_u_z[0]*bxbz[4]+grad_u_y[0]*bxby[4]+grad_u_x[0]*bxbx[4]+bxbz[1]*grad_u_z[2]+bxby[1]*grad_u_y[2]+bxbx[1]*grad_u_x[2]+grad_u_z[1]*bxbz[2]+grad_u_y[1]*bxby[2]+grad_u_x[1]*bxbx[2])*dx1; 
  bb_grad_u[5] += 0.3535533905932737*(bxbz[2]*grad_u_z[7]+bxby[2]*grad_u_y[7]+bxbx[2]*grad_u_x[7]+grad_u_z[2]*bxbz[7]+grad_u_y[2]*bxby[7]+grad_u_x[2]*bxbx[7]+bxbz[4]*grad_u_z[6]+bxby[4]*grad_u_y[6]+bxbx[4]*grad_u_x[6]+grad_u_z[4]*bxbz[6]+grad_u_y[4]*bxby[6]+grad_u_x[4]*bxbx[6]+bxbz[0]*grad_u_z[5]+bxby[0]*grad_u_y[5]+bxbx[0]*grad_u_x[5]+grad_u_z[0]*bxbz[5]+grad_u_y[0]*bxby[5]+grad_u_x[0]*bxbx[5]+bxbz[1]*grad_u_z[3]+bxby[1]*grad_u_y[3]+bxbx[1]*grad_u_x[3]+grad_u_z[1]*bxbz[3]+grad_u_y[1]*bxby[3]+grad_u_x[1]*bxbx[3])*dx1; 
  bb_grad_u[6] += 0.3535533905932737*(bxbz[1]*grad_u_z[7]+bxby[1]*grad_u_y[7]+bxbx[1]*grad_u_x[7]+grad_u_z[1]*bxbz[7]+grad_u_y[1]*bxby[7]+grad_u_x[1]*bxbx[7]+bxbz[0]*grad_u_z[6]+bxby[0]*grad_u_y[6]+bxbx[0]*grad_u_x[6]+grad_u_z[0]*bxbz[6]+grad_u_y[0]*bxby[6]+grad_u_x[0]*bxbx[6]+bxbz[4]*grad_u_z[5]+bxby[4]*grad_u_y[5]+bxbx[4]*grad_u_x[5]+grad_u_z[4]*bxbz[5]+grad_u_y[4]*bxby[5]+grad_u_x[4]*bxbx[5]+bxbz[2]*grad_u_z[3]+bxby[2]*grad_u_y[3]+bxbx[2]*grad_u_x[3]+grad_u_z[2]*bxbz[3]+grad_u_y[2]*bxby[3]+grad_u_x[2]*bxbx[3])*dx1; 
  bb_grad_u[7] += 0.3535533905932737*(bxbz[0]*grad_u_z[7]+bxby[0]*grad_u_y[7]+bxbx[0]*grad_u_x[7]+grad_u_z[0]*bxbz[7]+grad_u_y[0]*bxby[7]+grad_u_x[0]*bxbx[7]+bxbz[1]*grad_u_z[6]+bxby[1]*grad_u_y[6]+bxbx[1]*grad_u_x[6]+grad_u_z[1]*bxbz[6]+grad_u_y[1]*bxby[6]+grad_u_x[1]*bxbx[6]+bxbz[2]*grad_u_z[5]+bxby[2]*grad_u_y[5]+bxbx[2]*grad_u_x[5]+grad_u_z[2]*bxbz[5]+grad_u_y[2]*bxby[5]+grad_u_x[2]*bxbx[5]+bxbz[3]*grad_u_z[4]+bxby[3]*grad_u_y[4]+bxbx[3]*grad_u_x[4]+grad_u_z[3]*bxbz[4]+grad_u_y[3]*bxby[4]+grad_u_x[3]*bxbx[4])*dx1; 

  div_p_x[0] += (4.871392896287466*rhoux_r[1]-4.871392896287466*rhoux_l[1]-2.8125*(rhoux_r[0]+rhoux_l[0])+5.625*rhoux_c[0])*dx14*nuHyp+((-0.2886751345948129*(Pxx_r[1]+Pxx_l[1]))+0.5773502691896258*Pxx_c[1]+0.25*Pxx_r[0]-0.25*Pxx_l[0])*dx1; 
  div_p_x[1] += (72.1875*(rhoux_r[1]+rhoux_l[1])+249.375*rhoux_c[1]-56.83291712335378*rhoux_r[0]+56.83291712335378*rhoux_l[0])*dx14*nuHyp+((-0.5*Pxx_r[1])+0.5*Pxx_l[1]+0.4330127018922193*(Pxx_r[0]+Pxx_l[0])-0.8660254037844386*Pxx_c[0])*dx1; 
  div_p_x[2] += (4.871392896287466*rhoux_r[4]-4.871392896287466*rhoux_l[4]-2.8125*(rhoux_r[2]+rhoux_l[2])+5.625*rhoux_c[2])*dx14*nuHyp+((-0.2886751345948129*(Pxx_r[4]+Pxx_l[4]))+0.5773502691896258*Pxx_c[4]+0.25*Pxx_r[2]-0.25*Pxx_l[2])*dx1; 
  div_p_x[3] += (4.871392896287466*rhoux_r[5]-4.871392896287466*rhoux_l[5]-2.8125*(rhoux_r[3]+rhoux_l[3])+5.625*rhoux_c[3])*dx14*nuHyp+((-0.2886751345948129*(Pxx_r[5]+Pxx_l[5]))+0.5773502691896258*Pxx_c[5]+0.25*Pxx_r[3]-0.25*Pxx_l[3])*dx1; 
  div_p_x[4] += (72.1875*(rhoux_r[4]+rhoux_l[4])+249.375*rhoux_c[4]-56.83291712335378*rhoux_r[2]+56.83291712335378*rhoux_l[2])*dx14*nuHyp+((-0.5*Pxx_r[4])+0.5*Pxx_l[4]+0.4330127018922193*(Pxx_r[2]+Pxx_l[2])-0.8660254037844386*Pxx_c[2])*dx1; 
  div_p_x[5] += (72.1875*(rhoux_r[5]+rhoux_l[5])+249.375*rhoux_c[5]-56.83291712335378*rhoux_r[3]+56.83291712335378*rhoux_l[3])*dx14*nuHyp+((-0.5*Pxx_r[5])+0.5*Pxx_l[5]+0.4330127018922193*(Pxx_r[3]+Pxx_l[3])-0.8660254037844386*Pxx_c[3])*dx1; 
  div_p_x[6] += (4.871392896287466*rhoux_r[7]-4.871392896287466*rhoux_l[7]-2.8125*(rhoux_r[6]+rhoux_l[6])+5.625*rhoux_c[6])*dx14*nuHyp+((-0.2886751345948129*(Pxx_r[7]+Pxx_l[7]))+0.5773502691896258*Pxx_c[7]+0.25*Pxx_r[6]-0.25*Pxx_l[6])*dx1; 
  div_p_x[7] += (72.1875*(rhoux_r[7]+rhoux_l[7])+249.375*rhoux_c[7]-56.83291712335378*rhoux_r[6]+56.83291712335378*rhoux_l[6])*dx14*nuHyp+((-0.5*Pxx_r[7])+0.5*Pxx_l[7]+0.4330127018922193*(Pxx_r[6]+Pxx_l[6])-0.8660254037844386*Pxx_c[6])*dx1; 

  div_p_y[0] += (4.871392896287466*rhouy_r[1]-4.871392896287466*rhouy_l[1]-2.8125*(rhouy_r[0]+rhouy_l[0])+5.625*rhouy_c[0])*dx14*nuHyp+((-0.2886751345948129*(Pxy_r[1]+Pxy_l[1]))+0.5773502691896258*Pxy_c[1]+0.25*Pxy_r[0]-0.25*Pxy_l[0])*dx1; 
  div_p_y[1] += (72.1875*(rhouy_r[1]+rhouy_l[1])+249.375*rhouy_c[1]-56.83291712335378*rhouy_r[0]+56.83291712335378*rhouy_l[0])*dx14*nuHyp+((-0.5*Pxy_r[1])+0.5*Pxy_l[1]+0.4330127018922193*(Pxy_r[0]+Pxy_l[0])-0.8660254037844386*Pxy_c[0])*dx1; 
  div_p_y[2] += (4.871392896287466*rhouy_r[4]-4.871392896287466*rhouy_l[4]-2.8125*(rhouy_r[2]+rhouy_l[2])+5.625*rhouy_c[2])*dx14*nuHyp+((-0.2886751345948129*(Pxy_r[4]+Pxy_l[4]))+0.5773502691896258*Pxy_c[4]+0.25*Pxy_r[2]-0.25*Pxy_l[2])*dx1; 
  div_p_y[3] += (4.871392896287466*rhouy_r[5]-4.871392896287466*rhouy_l[5]-2.8125*(rhouy_r[3]+rhouy_l[3])+5.625*rhouy_c[3])*dx14*nuHyp+((-0.2886751345948129*(Pxy_r[5]+Pxy_l[5]))+0.5773502691896258*Pxy_c[5]+0.25*Pxy_r[3]-0.25*Pxy_l[3])*dx1; 
  div_p_y[4] += (72.1875*(rhouy_r[4]+rhouy_l[4])+249.375*rhouy_c[4]-56.83291712335378*rhouy_r[2]+56.83291712335378*rhouy_l[2])*dx14*nuHyp+((-0.5*Pxy_r[4])+0.5*Pxy_l[4]+0.4330127018922193*(Pxy_r[2]+Pxy_l[2])-0.8660254037844386*Pxy_c[2])*dx1; 
  div_p_y[5] += (72.1875*(rhouy_r[5]+rhouy_l[5])+249.375*rhouy_c[5]-56.83291712335378*rhouy_r[3]+56.83291712335378*rhouy_l[3])*dx14*nuHyp+((-0.5*Pxy_r[5])+0.5*Pxy_l[5]+0.4330127018922193*(Pxy_r[3]+Pxy_l[3])-0.8660254037844386*Pxy_c[3])*dx1; 
  div_p_y[6] += (4.871392896287466*rhouy_r[7]-4.871392896287466*rhouy_l[7]-2.8125*(rhouy_r[6]+rhouy_l[6])+5.625*rhouy_c[6])*dx14*nuHyp+((-0.2886751345948129*(Pxy_r[7]+Pxy_l[7]))+0.5773502691896258*Pxy_c[7]+0.25*Pxy_r[6]-0.25*Pxy_l[6])*dx1; 
  div_p_y[7] += (72.1875*(rhouy_r[7]+rhouy_l[7])+249.375*rhouy_c[7]-56.83291712335378*rhouy_r[6]+56.83291712335378*rhouy_l[6])*dx14*nuHyp+((-0.5*Pxy_r[7])+0.5*Pxy_l[7]+0.4330127018922193*(Pxy_r[6]+Pxy_l[6])-0.8660254037844386*Pxy_c[6])*dx1; 

  div_p_z[0] += (4.871392896287466*rhouz_r[1]-4.871392896287466*rhouz_l[1]-2.8125*(rhouz_r[0]+rhouz_l[0])+5.625*rhouz_c[0])*dx14*nuHyp+((-0.2886751345948129*(Pxz_r[1]+Pxz_l[1]))+0.5773502691896258*Pxz_c[1]+0.25*Pxz_r[0]-0.25*Pxz_l[0])*dx1; 
  div_p_z[1] += (72.1875*(rhouz_r[1]+rhouz_l[1])+249.375*rhouz_c[1]-56.83291712335378*rhouz_r[0]+56.83291712335378*rhouz_l[0])*dx14*nuHyp+((-0.5*Pxz_r[1])+0.5*Pxz_l[1]+0.4330127018922193*(Pxz_r[0]+Pxz_l[0])-0.8660254037844386*Pxz_c[0])*dx1; 
  div_p_z[2] += (4.871392896287466*rhouz_r[4]-4.871392896287466*rhouz_l[4]-2.8125*(rhouz_r[2]+rhouz_l[2])+5.625*rhouz_c[2])*dx14*nuHyp+((-0.2886751345948129*(Pxz_r[4]+Pxz_l[4]))+0.5773502691896258*Pxz_c[4]+0.25*Pxz_r[2]-0.25*Pxz_l[2])*dx1; 
  div_p_z[3] += (4.871392896287466*rhouz_r[5]-4.871392896287466*rhouz_l[5]-2.8125*(rhouz_r[3]+rhouz_l[3])+5.625*rhouz_c[3])*dx14*nuHyp+((-0.2886751345948129*(Pxz_r[5]+Pxz_l[5]))+0.5773502691896258*Pxz_c[5]+0.25*Pxz_r[3]-0.25*Pxz_l[3])*dx1; 
  div_p_z[4] += (72.1875*(rhouz_r[4]+rhouz_l[4])+249.375*rhouz_c[4]-56.83291712335378*rhouz_r[2]+56.83291712335378*rhouz_l[2])*dx14*nuHyp+((-0.5*Pxz_r[4])+0.5*Pxz_l[4]+0.4330127018922193*(Pxz_r[2]+Pxz_l[2])-0.8660254037844386*Pxz_c[2])*dx1; 
  div_p_z[5] += (72.1875*(rhouz_r[5]+rhouz_l[5])+249.375*rhouz_c[5]-56.83291712335378*rhouz_r[3]+56.83291712335378*rhouz_l[3])*dx14*nuHyp+((-0.5*Pxz_r[5])+0.5*Pxz_l[5]+0.4330127018922193*(Pxz_r[3]+Pxz_l[3])-0.8660254037844386*Pxz_c[3])*dx1; 
  div_p_z[6] += (4.871392896287466*rhouz_r[7]-4.871392896287466*rhouz_l[7]-2.8125*(rhouz_r[6]+rhouz_l[6])+5.625*rhouz_c[6])*dx14*nuHyp+((-0.2886751345948129*(Pxz_r[7]+Pxz_l[7]))+0.5773502691896258*Pxz_c[7]+0.25*Pxz_r[6]-0.25*Pxz_l[6])*dx1; 
  div_p_z[7] += (72.1875*(rhouz_r[7]+rhouz_l[7])+249.375*rhouz_c[7]-56.83291712335378*rhouz_r[6]+56.83291712335378*rhouz_l[6])*dx14*nuHyp+((-0.5*Pxz_r[7])+0.5*Pxz_l[7]+0.4330127018922193*(Pxz_r[6]+Pxz_l[6])-0.8660254037844386*Pxz_c[6])*dx1; 

  double div_ppar_b[8] = {0.0}; 
  double rho_inv[8] = {0.0}; 
  ser_3x_p1_inv(rho, rho_inv); 
  div_ppar_b[0] = (-0.2886751345948129*ppar_b_r[1])-0.2886751345948129*ppar_b_l[1]+0.5773502691896258*ppar_b_c[1]+0.25*ppar_b_r[0]-0.25*ppar_b_l[0]; 
  div_ppar_b[1] = (-0.5*ppar_b_r[1])+0.5*ppar_b_l[1]+0.4330127018922193*ppar_b_r[0]+0.4330127018922193*ppar_b_l[0]-0.8660254037844386*ppar_b_c[0]; 
  div_ppar_b[2] = (-0.2886751345948129*ppar_b_r[4])-0.2886751345948129*ppar_b_l[4]+0.5773502691896258*ppar_b_c[4]+0.25*ppar_b_r[2]-0.25*ppar_b_l[2]; 
  div_ppar_b[3] = (-0.2886751345948129*ppar_b_r[5])-0.2886751345948129*ppar_b_l[5]+0.5773502691896258*ppar_b_c[5]+0.25*ppar_b_r[3]-0.25*ppar_b_l[3]; 
  div_ppar_b[4] = (-0.5*ppar_b_r[4])+0.5*ppar_b_l[4]+0.4330127018922193*ppar_b_r[2]+0.4330127018922193*ppar_b_l[2]-0.8660254037844386*ppar_b_c[2]; 
  div_ppar_b[5] = (-0.5*ppar_b_r[5])+0.5*ppar_b_l[5]+0.4330127018922193*ppar_b_r[3]+0.4330127018922193*ppar_b_l[3]-0.8660254037844386*ppar_b_c[3]; 
  div_ppar_b[6] = (-0.2886751345948129*ppar_b_r[7])-0.2886751345948129*ppar_b_l[7]+0.5773502691896258*ppar_b_c[7]+0.25*ppar_b_r[6]-0.25*ppar_b_l[6]; 
  div_ppar_b[7] = (-0.5*ppar_b_r[7])+0.5*ppar_b_l[7]+0.4330127018922193*ppar_b_r[6]+0.4330127018922193*ppar_b_l[6]-0.8660254037844386*ppar_b_c[6]; 

  p_force[0] += 0.3535533905932737*(div_ppar_b[7]*rho_inv[7]+div_ppar_b[6]*rho_inv[6]+div_ppar_b[5]*rho_inv[5]+div_ppar_b[4]*rho_inv[4]+div_ppar_b[3]*rho_inv[3]+div_ppar_b[2]*rho_inv[2]+div_ppar_b[1]*rho_inv[1]+div_ppar_b[0]*rho_inv[0])*dx1; 
  p_force[1] += 0.3535533905932737*(div_ppar_b[6]*rho_inv[7]+rho_inv[6]*div_ppar_b[7]+div_ppar_b[3]*rho_inv[5]+rho_inv[3]*div_ppar_b[5]+div_ppar_b[2]*rho_inv[4]+rho_inv[2]*div_ppar_b[4]+div_ppar_b[0]*rho_inv[1]+rho_inv[0]*div_ppar_b[1])*dx1; 
  p_force[2] += 0.3535533905932737*(div_ppar_b[5]*rho_inv[7]+rho_inv[5]*div_ppar_b[7]+div_ppar_b[3]*rho_inv[6]+rho_inv[3]*div_ppar_b[6]+div_ppar_b[1]*rho_inv[4]+rho_inv[1]*div_ppar_b[4]+div_ppar_b[0]*rho_inv[2]+rho_inv[0]*div_ppar_b[2])*dx1; 
  p_force[3] += 0.3535533905932737*(div_ppar_b[4]*rho_inv[7]+rho_inv[4]*div_ppar_b[7]+div_ppar_b[2]*rho_inv[6]+rho_inv[2]*div_ppar_b[6]+div_ppar_b[1]*rho_inv[5]+rho_inv[1]*div_ppar_b[5]+div_ppar_b[0]*rho_inv[3]+rho_inv[0]*div_ppar_b[3])*dx1; 
  p_force[4] += 0.3535533905932737*(div_ppar_b[3]*rho_inv[7]+rho_inv[3]*div_ppar_b[7]+div_ppar_b[5]*rho_inv[6]+rho_inv[5]*div_ppar_b[6]+div_ppar_b[0]*rho_inv[4]+rho_inv[0]*div_ppar_b[4]+div_ppar_b[1]*rho_inv[2]+rho_inv[1]*div_ppar_b[2])*dx1; 
  p_force[5] += 0.3535533905932737*(div_ppar_b[2]*rho_inv[7]+rho_inv[2]*div_ppar_b[7]+div_ppar_b[4]*rho_inv[6]+rho_inv[4]*div_ppar_b[6]+div_ppar_b[0]*rho_inv[5]+rho_inv[0]*div_ppar_b[5]+div_ppar_b[1]*rho_inv[3]+rho_inv[1]*div_ppar_b[3])*dx1; 
  p_force[6] += 0.3535533905932737*(div_ppar_b[1]*rho_inv[7]+rho_inv[1]*div_ppar_b[7]+div_ppar_b[0]*rho_inv[6]+rho_inv[0]*div_ppar_b[6]+div_ppar_b[4]*rho_inv[5]+rho_inv[4]*div_ppar_b[5]+div_ppar_b[2]*rho_inv[3]+rho_inv[2]*div_ppar_b[3])*dx1; 
  p_force[7] += 0.3535533905932737*(div_ppar_b[0]*rho_inv[7]+rho_inv[0]*div_ppar_b[7]+div_ppar_b[1]*rho_inv[6]+rho_inv[1]*div_ppar_b[6]+div_ppar_b[2]*rho_inv[5]+rho_inv[2]*div_ppar_b[5]+div_ppar_b[3]*rho_inv[4]+rho_inv[3]*div_ppar_b[4])*dx1; 

  p_perp_source[0] += ((-4.871392896287466*p_perp_r[1])+4.871392896287466*p_perp_l[1]+2.8125*(p_perp_r[0]+p_perp_l[0])-5.625*p_perp_c[0])*dx14*nuHyp-0.3535533905932737*(grad_u_x[7]*p_perp_c[7]+grad_u_x[6]*p_perp_c[6]+grad_u_x[5]*p_perp_c[5]+grad_u_x[4]*p_perp_c[4]+grad_u_x[3]*p_perp_c[3]+grad_u_x[2]*p_perp_c[2]+grad_u_x[1]*p_perp_c[1]+grad_u_x[0]*p_perp_c[0])*dx1+0.3535533905932737*(bb_grad_u[7]*p_perp_c[7]+bb_grad_u[6]*p_perp_c[6]+bb_grad_u[5]*p_perp_c[5]+bb_grad_u[4]*p_perp_c[4]+bb_grad_u[3]*p_perp_c[3]+bb_grad_u[2]*p_perp_c[2]+bb_grad_u[1]*p_perp_c[1]+bb_grad_u[0]*p_perp_c[0]); 
  p_perp_source[1] += ((-72.1875*(p_perp_r[1]+p_perp_l[1]))-249.375*p_perp_c[1]+56.83291712335378*p_perp_r[0]-56.83291712335378*p_perp_l[0])*dx14*nuHyp-0.3535533905932737*(grad_u_x[6]*p_perp_c[7]+p_perp_c[6]*grad_u_x[7]+grad_u_x[3]*p_perp_c[5]+p_perp_c[3]*grad_u_x[5]+grad_u_x[2]*p_perp_c[4]+p_perp_c[2]*grad_u_x[4]+grad_u_x[0]*p_perp_c[1]+p_perp_c[0]*grad_u_x[1])*dx1+0.3535533905932737*(bb_grad_u[6]*p_perp_c[7]+p_perp_c[6]*bb_grad_u[7]+bb_grad_u[3]*p_perp_c[5]+p_perp_c[3]*bb_grad_u[5]+bb_grad_u[2]*p_perp_c[4]+p_perp_c[2]*bb_grad_u[4]+bb_grad_u[0]*p_perp_c[1]+p_perp_c[0]*bb_grad_u[1]); 
  p_perp_source[2] += ((-4.871392896287466*p_perp_r[4])+4.871392896287466*p_perp_l[4]+2.8125*(p_perp_r[2]+p_perp_l[2])-5.625*p_perp_c[2])*dx14*nuHyp-0.3535533905932737*(grad_u_x[5]*p_perp_c[7]+p_perp_c[5]*grad_u_x[7]+grad_u_x[3]*p_perp_c[6]+p_perp_c[3]*grad_u_x[6]+grad_u_x[1]*p_perp_c[4]+p_perp_c[1]*grad_u_x[4]+grad_u_x[0]*p_perp_c[2]+p_perp_c[0]*grad_u_x[2])*dx1+0.3535533905932737*(bb_grad_u[5]*p_perp_c[7]+p_perp_c[5]*bb_grad_u[7]+bb_grad_u[3]*p_perp_c[6]+p_perp_c[3]*bb_grad_u[6]+bb_grad_u[1]*p_perp_c[4]+p_perp_c[1]*bb_grad_u[4]+bb_grad_u[0]*p_perp_c[2]+p_perp_c[0]*bb_grad_u[2]); 
  p_perp_source[3] += ((-4.871392896287466*p_perp_r[5])+4.871392896287466*p_perp_l[5]+2.8125*(p_perp_r[3]+p_perp_l[3])-5.625*p_perp_c[3])*dx14*nuHyp-0.3535533905932737*(grad_u_x[4]*p_perp_c[7]+p_perp_c[4]*grad_u_x[7]+grad_u_x[2]*p_perp_c[6]+p_perp_c[2]*grad_u_x[6]+grad_u_x[1]*p_perp_c[5]+p_perp_c[1]*grad_u_x[5]+grad_u_x[0]*p_perp_c[3]+p_perp_c[0]*grad_u_x[3])*dx1+0.3535533905932737*(bb_grad_u[4]*p_perp_c[7]+p_perp_c[4]*bb_grad_u[7]+bb_grad_u[2]*p_perp_c[6]+p_perp_c[2]*bb_grad_u[6]+bb_grad_u[1]*p_perp_c[5]+p_perp_c[1]*bb_grad_u[5]+bb_grad_u[0]*p_perp_c[3]+p_perp_c[0]*bb_grad_u[3]); 
  p_perp_source[4] += ((-72.1875*(p_perp_r[4]+p_perp_l[4]))-249.375*p_perp_c[4]+56.83291712335378*p_perp_r[2]-56.83291712335378*p_perp_l[2])*dx14*nuHyp-0.3535533905932737*(grad_u_x[3]*p_perp_c[7]+p_perp_c[3]*grad_u_x[7]+grad_u_x[5]*p_perp_c[6]+p_perp_c[5]*grad_u_x[6]+grad_u_x[0]*p_perp_c[4]+p_perp_c[0]*grad_u_x[4]+grad_u_x[1]*p_perp_c[2]+p_perp_c[1]*grad_u_x[2])*dx1+0.3535533905932737*(bb_grad_u[3]*p_perp_c[7]+p_perp_c[3]*bb_grad_u[7]+bb_grad_u[5]*p_perp_c[6]+p_perp_c[5]*bb_grad_u[6]+bb_grad_u[0]*p_perp_c[4]+p_perp_c[0]*bb_grad_u[4]+bb_grad_u[1]*p_perp_c[2]+p_perp_c[1]*bb_grad_u[2]); 
  p_perp_source[5] += ((-72.1875*(p_perp_r[5]+p_perp_l[5]))-249.375*p_perp_c[5]+56.83291712335378*p_perp_r[3]-56.83291712335378*p_perp_l[3])*dx14*nuHyp-0.3535533905932737*(grad_u_x[2]*p_perp_c[7]+p_perp_c[2]*grad_u_x[7]+grad_u_x[4]*p_perp_c[6]+p_perp_c[4]*grad_u_x[6]+grad_u_x[0]*p_perp_c[5]+p_perp_c[0]*grad_u_x[5]+grad_u_x[1]*p_perp_c[3]+p_perp_c[1]*grad_u_x[3])*dx1+0.3535533905932737*(bb_grad_u[2]*p_perp_c[7]+p_perp_c[2]*bb_grad_u[7]+bb_grad_u[4]*p_perp_c[6]+p_perp_c[4]*bb_grad_u[6]+bb_grad_u[0]*p_perp_c[5]+p_perp_c[0]*bb_grad_u[5]+bb_grad_u[1]*p_perp_c[3]+p_perp_c[1]*bb_grad_u[3]); 
  p_perp_source[6] += ((-4.871392896287466*p_perp_r[7])+4.871392896287466*p_perp_l[7]+2.8125*(p_perp_r[6]+p_perp_l[6])-5.625*p_perp_c[6])*dx14*nuHyp-0.3535533905932737*(grad_u_x[1]*p_perp_c[7]+p_perp_c[1]*grad_u_x[7]+grad_u_x[0]*p_perp_c[6]+p_perp_c[0]*grad_u_x[6]+grad_u_x[4]*p_perp_c[5]+p_perp_c[4]*grad_u_x[5]+grad_u_x[2]*p_perp_c[3]+p_perp_c[2]*grad_u_x[3])*dx1+0.3535533905932737*(bb_grad_u[1]*p_perp_c[7]+p_perp_c[1]*bb_grad_u[7]+bb_grad_u[0]*p_perp_c[6]+p_perp_c[0]*bb_grad_u[6]+bb_grad_u[4]*p_perp_c[5]+p_perp_c[4]*bb_grad_u[5]+bb_grad_u[2]*p_perp_c[3]+p_perp_c[2]*bb_grad_u[3]); 
  p_perp_source[7] += ((-72.1875*(p_perp_r[7]+p_perp_l[7]))-249.375*p_perp_c[7]+56.83291712335378*p_perp_r[6]-56.83291712335378*p_perp_l[6])*dx14*nuHyp-0.3535533905932737*(grad_u_x[0]*p_perp_c[7]+p_perp_c[0]*grad_u_x[7]+grad_u_x[1]*p_perp_c[6]+p_perp_c[1]*grad_u_x[6]+grad_u_x[2]*p_perp_c[5]+p_perp_c[2]*grad_u_x[5]+grad_u_x[3]*p_perp_c[4]+p_perp_c[3]*grad_u_x[4])*dx1+0.3535533905932737*(bb_grad_u[0]*p_perp_c[7]+p_perp_c[0]*bb_grad_u[7]+bb_grad_u[1]*p_perp_c[6]+p_perp_c[1]*bb_grad_u[6]+bb_grad_u[2]*p_perp_c[5]+p_perp_c[2]*bb_grad_u[5]+bb_grad_u[3]*p_perp_c[4]+p_perp_c[3]*bb_grad_u[4]); 

} 
