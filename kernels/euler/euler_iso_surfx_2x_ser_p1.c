#include <gkyl_euler_kernels.h> 
#include <gkyl_basis_ser_2x_p1_surfx1_eval_quad.h> 
#include <gkyl_basis_ser_2x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void euler_iso_surfx_2x_ser_p1(const double *w, const double *dxv, const double vth, const double *ul, const double *uc, const double *ur, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // gas_gamma: Adiabatic index.
  // ul/uc/ur: [ux, uy, uz] Fluid flow in left/center/right cells.
  // pl/pc/pr: Fluid pressure in left/center/right cells.
  // statevecl/statevecc/statevecr: [rho, rho ux, rho uy, rho uz, energy], Fluid input state vector in left/center/right cells.
  // out: Incremented output.
  const double dx1 = 2.0/dxv[0]; 
  const double *rho_l = &statevecl[0]; 
  const double *rhou0_l = &statevecl[4]; 
  const double *rhou1_l = &statevecl[8]; 
  const double *rhou2_l = &statevecl[12]; 
  const double *energy_l = &statevecl[16]; 

  const double *rho_c = &statevecc[0]; 
  const double *rhou0_c = &statevecc[4]; 
  const double *rhou1_c = &statevecc[8]; 
  const double *rhou2_c = &statevecc[12]; 
  const double *energy_c = &statevecc[16]; 

  const double *rho_r = &statevecr[0]; 
  const double *rhou0_r = &statevecr[4]; 
  const double *rhou1_r = &statevecr[8]; 
  const double *rhou2_r = &statevecr[12]; 
  const double *energy_r = &statevecr[16]; 

  const double *ul_0 = &ul[0]; 
  const double *uc_0 = &uc[0]; 
  const double *ur_0 = &ur[0]; 

  const double *ul_1 = &ul[4]; 
  const double *uc_1 = &uc[4]; 
  const double *ur_1 = &ur[4]; 

  const double *ul_2 = &ul[8]; 
  const double *uc_2 = &uc[8]; 
  const double *ur_2 = &ur[8]; 

  double *outrho = &out[0]; 
  double *outrhou0 = &out[4]; 
  double *outrhou1 = &out[8]; 
  double *outrhou2 = &out[12]; 
  double *outenergy = &out[16]; 

  double vthsq = vth*vth; 
  double uQuad_l[2] = {0.0};
  double uQuad_r[2] = {0.0};
  double uMax_l[2] = {0.0};;
  double uMax_r[2] = {0.0};
  double Ghat_rho_l[2] = {0.0}; 
  double Ghat_rho_r[2] = {0.0}; 

  double urec_0_l[2] = {0.0}; 
  double urec_0_r[2] = {0.0}; 
  double urec_1_l[2] = {0.0}; 
  double urec_1_r[2] = {0.0}; 
  double urec_2_l[2] = {0.0}; 
  double urec_2_r[2] = {0.0}; 
  double rhorec_l[2] = {0.0}; 
  double rhorec_r[2] = {0.0}; 
  double u_l_r = 0.0; 
  double u_c_l = 0.0; 
  double u_c_r = 0.0; 
  double u_r_l = 0.0; 

  u_l_r = ser_2x_p1_surfx1_eval_quad_node_0_r(ul_0); 
  u_c_l = ser_2x_p1_surfx1_eval_quad_node_0_l(uc_0); 
  u_c_r = ser_2x_p1_surfx1_eval_quad_node_0_r(uc_0); 
  u_r_l = ser_2x_p1_surfx1_eval_quad_node_0_l(ur_0); 
  uQuad_l[0] = fmax(fabs(u_l_r), fabs(u_c_l)); 
  uQuad_r[0] = fmax(fabs(u_c_r), fabs(u_r_l)); 
  u_l_r = ser_2x_p1_surfx1_eval_quad_node_1_r(ul_0); 
  u_c_l = ser_2x_p1_surfx1_eval_quad_node_1_l(uc_0); 
  u_c_r = ser_2x_p1_surfx1_eval_quad_node_1_r(uc_0); 
  u_r_l = ser_2x_p1_surfx1_eval_quad_node_1_l(ur_0); 
  uQuad_l[1] = fmax(fabs(u_l_r), fabs(u_c_l)); 
  uQuad_r[1] = fmax(fabs(u_c_r), fabs(u_r_l)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_2x_p1_upwind_quad_to_modal(uQuad_l, uMax_l); 
  ser_2x_p1_upwind_quad_to_modal(uQuad_r, uMax_r); 

  Ghat_rho_l[0] = 0.5303300858899105*rho_l[3]*ul_0[3]+0.3061862178478971*rho_l[2]*ul_0[3]+0.5303300858899105*rho_c[3]*uc_0[3]-0.3061862178478971*rho_c[2]*uc_0[3]+0.3061862178478971*ul_0[2]*rho_l[3]+0.4330127018922193*uMax_l[1]*rho_l[3]-0.3061862178478971*uc_0[2]*rho_c[3]+0.4330127018922193*uMax_l[1]*rho_c[3]+0.1767766952966368*rho_l[2]*ul_0[2]+0.1767766952966368*rho_c[2]*uc_0[2]+0.25*uMax_l[1]*rho_l[2]-0.25*uMax_l[1]*rho_c[2]+0.5303300858899105*rho_l[1]*ul_0[1]+0.3061862178478971*rho_l[0]*ul_0[1]+0.5303300858899105*rho_c[1]*uc_0[1]-0.3061862178478971*rho_c[0]*uc_0[1]+0.3061862178478971*ul_0[0]*rho_l[1]+0.4330127018922193*uMax_l[0]*rho_l[1]-0.3061862178478971*uc_0[0]*rho_c[1]+0.4330127018922193*uMax_l[0]*rho_c[1]+0.1767766952966368*rho_l[0]*ul_0[0]+0.1767766952966368*rho_c[0]*uc_0[0]+0.25*rho_l[0]*uMax_l[0]-0.25*rho_c[0]*uMax_l[0]; 
  Ghat_rho_l[1] = 0.5303300858899105*rho_l[1]*ul_0[3]+0.3061862178478971*rho_l[0]*ul_0[3]+0.5303300858899105*rho_c[1]*uc_0[3]-0.3061862178478971*rho_c[0]*uc_0[3]+0.5303300858899105*ul_0[1]*rho_l[3]+0.3061862178478971*ul_0[0]*rho_l[3]+0.4330127018922193*uMax_l[0]*rho_l[3]+0.5303300858899105*uc_0[1]*rho_c[3]-0.3061862178478971*uc_0[0]*rho_c[3]+0.4330127018922193*uMax_l[0]*rho_c[3]+0.3061862178478971*rho_l[1]*ul_0[2]+0.1767766952966368*rho_l[0]*ul_0[2]-0.3061862178478971*rho_c[1]*uc_0[2]+0.1767766952966368*rho_c[0]*uc_0[2]+0.3061862178478971*ul_0[1]*rho_l[2]+0.1767766952966368*ul_0[0]*rho_l[2]+0.25*uMax_l[0]*rho_l[2]-0.3061862178478971*uc_0[1]*rho_c[2]+0.1767766952966368*uc_0[0]*rho_c[2]-0.25*uMax_l[0]*rho_c[2]+0.4330127018922193*rho_l[1]*uMax_l[1]+0.4330127018922193*rho_c[1]*uMax_l[1]+0.25*rho_l[0]*uMax_l[1]-0.25*rho_c[0]*uMax_l[1]; 

  Ghat_rho_r[0] = 0.5303300858899105*rho_r[3]*ur_0[3]-0.3061862178478971*rho_r[2]*ur_0[3]+0.5303300858899105*rho_c[3]*uc_0[3]+0.3061862178478971*rho_c[2]*uc_0[3]-0.3061862178478971*ur_0[2]*rho_r[3]+0.4330127018922193*uMax_r[1]*rho_r[3]+0.3061862178478971*uc_0[2]*rho_c[3]+0.4330127018922193*uMax_r[1]*rho_c[3]+0.1767766952966368*rho_r[2]*ur_0[2]+0.1767766952966368*rho_c[2]*uc_0[2]-0.25*uMax_r[1]*rho_r[2]+0.25*uMax_r[1]*rho_c[2]+0.5303300858899105*rho_r[1]*ur_0[1]-0.3061862178478971*rho_r[0]*ur_0[1]+0.5303300858899105*rho_c[1]*uc_0[1]+0.3061862178478971*rho_c[0]*uc_0[1]-0.3061862178478971*ur_0[0]*rho_r[1]+0.4330127018922193*uMax_r[0]*rho_r[1]+0.3061862178478971*uc_0[0]*rho_c[1]+0.4330127018922193*uMax_r[0]*rho_c[1]+0.1767766952966368*rho_r[0]*ur_0[0]+0.1767766952966368*rho_c[0]*uc_0[0]-0.25*rho_r[0]*uMax_r[0]+0.25*rho_c[0]*uMax_r[0]; 
  Ghat_rho_r[1] = 0.5303300858899105*rho_r[1]*ur_0[3]-0.3061862178478971*rho_r[0]*ur_0[3]+0.5303300858899105*rho_c[1]*uc_0[3]+0.3061862178478971*rho_c[0]*uc_0[3]+0.5303300858899105*ur_0[1]*rho_r[3]-0.3061862178478971*ur_0[0]*rho_r[3]+0.4330127018922193*uMax_r[0]*rho_r[3]+0.5303300858899105*uc_0[1]*rho_c[3]+0.3061862178478971*uc_0[0]*rho_c[3]+0.4330127018922193*uMax_r[0]*rho_c[3]-0.3061862178478971*rho_r[1]*ur_0[2]+0.1767766952966368*rho_r[0]*ur_0[2]+0.3061862178478971*rho_c[1]*uc_0[2]+0.1767766952966368*rho_c[0]*uc_0[2]-0.3061862178478971*ur_0[1]*rho_r[2]+0.1767766952966368*ur_0[0]*rho_r[2]-0.25*uMax_r[0]*rho_r[2]+0.3061862178478971*uc_0[1]*rho_c[2]+0.1767766952966368*uc_0[0]*rho_c[2]+0.25*uMax_r[0]*rho_c[2]+0.4330127018922193*rho_r[1]*uMax_r[1]+0.4330127018922193*rho_c[1]*uMax_r[1]-0.25*rho_r[0]*uMax_r[1]+0.25*rho_c[0]*uMax_r[1]; 

  urec_0_l[0] = 0.408248290463863*ul_0[1]-0.408248290463863*uc_0[1]+0.3535533905932737*ul_0[0]+0.3535533905932737*uc_0[0]; 
  urec_0_l[1] = 0.408248290463863*ul_0[3]-0.408248290463863*uc_0[3]+0.3535533905932737*ul_0[2]+0.3535533905932737*uc_0[2]; 

  urec_1_l[0] = 0.408248290463863*ul_1[1]-0.408248290463863*uc_1[1]+0.3535533905932737*ul_1[0]+0.3535533905932737*uc_1[0]; 
  urec_1_l[1] = 0.408248290463863*ul_1[3]-0.408248290463863*uc_1[3]+0.3535533905932737*ul_1[2]+0.3535533905932737*uc_1[2]; 

  urec_2_l[0] = 0.408248290463863*ul_2[1]-0.408248290463863*uc_2[1]+0.3535533905932737*ul_2[0]+0.3535533905932737*uc_2[0]; 
  urec_2_l[1] = 0.408248290463863*ul_2[3]-0.408248290463863*uc_2[3]+0.3535533905932737*ul_2[2]+0.3535533905932737*uc_2[2]; 

  rhorec_l[0] = 0.408248290463863*rho_l[1]-0.408248290463863*rho_c[1]+0.3535533905932737*rho_l[0]+0.3535533905932737*rho_c[0]; 
  rhorec_l[1] = 0.408248290463863*rho_l[3]-0.408248290463863*rho_c[3]+0.3535533905932737*rho_l[2]+0.3535533905932737*rho_c[2]; 

  urec_0_r[0] = (-0.408248290463863*ur_0[1])+0.408248290463863*uc_0[1]+0.3535533905932737*ur_0[0]+0.3535533905932737*uc_0[0]; 
  urec_0_r[1] = (-0.408248290463863*ur_0[3])+0.408248290463863*uc_0[3]+0.3535533905932737*ur_0[2]+0.3535533905932737*uc_0[2]; 

  urec_1_r[0] = (-0.408248290463863*ur_1[1])+0.408248290463863*uc_1[1]+0.3535533905932737*ur_1[0]+0.3535533905932737*uc_1[0]; 
  urec_1_r[1] = (-0.408248290463863*ur_1[3])+0.408248290463863*uc_1[3]+0.3535533905932737*ur_1[2]+0.3535533905932737*uc_1[2]; 

  urec_2_r[0] = (-0.408248290463863*ur_2[1])+0.408248290463863*uc_2[1]+0.3535533905932737*ur_2[0]+0.3535533905932737*uc_2[0]; 
  urec_2_r[1] = (-0.408248290463863*ur_2[3])+0.408248290463863*uc_2[3]+0.3535533905932737*ur_2[2]+0.3535533905932737*uc_2[2]; 

  rhorec_r[0] = (-0.408248290463863*rho_r[1])+0.408248290463863*rho_c[1]+0.3535533905932737*rho_r[0]+0.3535533905932737*rho_c[0]; 
  rhorec_r[1] = (-0.408248290463863*rho_r[3])+0.408248290463863*rho_c[3]+0.3535533905932737*rho_r[2]+0.3535533905932737*rho_c[2]; 

  outrho[0] += 0.7071067811865475*Ghat_rho_l[0]*dx1-0.7071067811865475*Ghat_rho_r[0]*dx1; 
  outrho[1] += (-1.224744871391589*Ghat_rho_r[0]*dx1)-1.224744871391589*Ghat_rho_l[0]*dx1; 
  outrho[2] += 0.7071067811865475*Ghat_rho_l[1]*dx1-0.7071067811865475*Ghat_rho_r[1]*dx1; 
  outrho[3] += (-1.224744871391589*Ghat_rho_r[1]*dx1)-1.224744871391589*Ghat_rho_l[1]*dx1; 

  outrhou0[0] += (-0.7071067811865475*rhorec_r[0]*dx1*vthsq)+0.7071067811865475*rhorec_l[0]*dx1*vthsq-0.3061862178478971*uMax_r[1]*rhou0_r[3]*dx1+0.3061862178478971*uMax_l[1]*rhou0_l[3]*dx1-0.3061862178478971*uMax_r[1]*rhou0_c[3]*dx1+0.3061862178478971*uMax_l[1]*rhou0_c[3]*dx1+0.1767766952966368*uMax_r[1]*rhou0_r[2]*dx1+0.1767766952966368*uMax_l[1]*rhou0_l[2]*dx1-0.1767766952966368*uMax_r[1]*rhou0_c[2]*dx1-0.1767766952966368*uMax_l[1]*rhou0_c[2]*dx1-0.5*Ghat_rho_r[1]*urec_0_r[1]*dx1+0.5*Ghat_rho_l[1]*urec_0_l[1]*dx1-0.3061862178478971*uMax_r[0]*rhou0_r[1]*dx1+0.3061862178478971*uMax_l[0]*rhou0_l[1]*dx1-0.3061862178478971*uMax_r[0]*rhou0_c[1]*dx1+0.3061862178478971*uMax_l[0]*rhou0_c[1]*dx1-0.5*Ghat_rho_r[0]*urec_0_r[0]*dx1+0.5*Ghat_rho_l[0]*urec_0_l[0]*dx1+0.1767766952966368*rhou0_r[0]*uMax_r[0]*dx1-0.1767766952966368*rhou0_c[0]*uMax_r[0]*dx1+0.1767766952966368*rhou0_l[0]*uMax_l[0]*dx1-0.1767766952966368*rhou0_c[0]*uMax_l[0]*dx1; 
  outrhou0[1] += (-1.224744871391589*rhorec_r[0]*dx1*vthsq)-1.224744871391589*rhorec_l[0]*dx1*vthsq-0.5303300858899105*uMax_r[1]*rhou0_r[3]*dx1-0.5303300858899105*uMax_l[1]*rhou0_l[3]*dx1-0.5303300858899105*uMax_r[1]*rhou0_c[3]*dx1-0.5303300858899105*uMax_l[1]*rhou0_c[3]*dx1+0.3061862178478971*uMax_r[1]*rhou0_r[2]*dx1-0.3061862178478971*uMax_l[1]*rhou0_l[2]*dx1-0.3061862178478971*uMax_r[1]*rhou0_c[2]*dx1+0.3061862178478971*uMax_l[1]*rhou0_c[2]*dx1-0.8660254037844386*Ghat_rho_r[1]*urec_0_r[1]*dx1-0.8660254037844386*Ghat_rho_l[1]*urec_0_l[1]*dx1-0.5303300858899105*uMax_r[0]*rhou0_r[1]*dx1-0.5303300858899105*uMax_l[0]*rhou0_l[1]*dx1-0.5303300858899105*uMax_r[0]*rhou0_c[1]*dx1-0.5303300858899105*uMax_l[0]*rhou0_c[1]*dx1-0.8660254037844386*Ghat_rho_r[0]*urec_0_r[0]*dx1-0.8660254037844386*Ghat_rho_l[0]*urec_0_l[0]*dx1+0.3061862178478971*rhou0_r[0]*uMax_r[0]*dx1-0.3061862178478971*rhou0_c[0]*uMax_r[0]*dx1-0.3061862178478971*rhou0_l[0]*uMax_l[0]*dx1+0.3061862178478971*rhou0_c[0]*uMax_l[0]*dx1; 
  outrhou0[2] += (-0.7071067811865475*rhorec_r[1]*dx1*vthsq)+0.7071067811865475*rhorec_l[1]*dx1*vthsq-0.3061862178478971*uMax_r[0]*rhou0_r[3]*dx1+0.3061862178478971*uMax_l[0]*rhou0_l[3]*dx1-0.3061862178478971*uMax_r[0]*rhou0_c[3]*dx1+0.3061862178478971*uMax_l[0]*rhou0_c[3]*dx1+0.1767766952966368*uMax_r[0]*rhou0_r[2]*dx1+0.1767766952966368*uMax_l[0]*rhou0_l[2]*dx1-0.1767766952966368*uMax_r[0]*rhou0_c[2]*dx1-0.1767766952966368*uMax_l[0]*rhou0_c[2]*dx1-0.5*Ghat_rho_r[0]*urec_0_r[1]*dx1+0.5*Ghat_rho_l[0]*urec_0_l[1]*dx1-0.3061862178478971*rhou0_r[1]*uMax_r[1]*dx1-0.3061862178478971*rhou0_c[1]*uMax_r[1]*dx1+0.1767766952966368*rhou0_r[0]*uMax_r[1]*dx1-0.1767766952966368*rhou0_c[0]*uMax_r[1]*dx1+0.3061862178478971*rhou0_l[1]*uMax_l[1]*dx1+0.3061862178478971*rhou0_c[1]*uMax_l[1]*dx1+0.1767766952966368*rhou0_l[0]*uMax_l[1]*dx1-0.1767766952966368*rhou0_c[0]*uMax_l[1]*dx1-0.5*urec_0_r[0]*Ghat_rho_r[1]*dx1+0.5*urec_0_l[0]*Ghat_rho_l[1]*dx1; 
  outrhou0[3] += (-1.224744871391589*rhorec_r[1]*dx1*vthsq)-1.224744871391589*rhorec_l[1]*dx1*vthsq-0.5303300858899105*uMax_r[0]*rhou0_r[3]*dx1-0.5303300858899105*uMax_l[0]*rhou0_l[3]*dx1-0.5303300858899105*uMax_r[0]*rhou0_c[3]*dx1-0.5303300858899105*uMax_l[0]*rhou0_c[3]*dx1+0.3061862178478971*uMax_r[0]*rhou0_r[2]*dx1-0.3061862178478971*uMax_l[0]*rhou0_l[2]*dx1-0.3061862178478971*uMax_r[0]*rhou0_c[2]*dx1+0.3061862178478971*uMax_l[0]*rhou0_c[2]*dx1-0.8660254037844386*Ghat_rho_r[0]*urec_0_r[1]*dx1-0.8660254037844386*Ghat_rho_l[0]*urec_0_l[1]*dx1-0.5303300858899105*rhou0_r[1]*uMax_r[1]*dx1-0.5303300858899105*rhou0_c[1]*uMax_r[1]*dx1+0.3061862178478971*rhou0_r[0]*uMax_r[1]*dx1-0.3061862178478971*rhou0_c[0]*uMax_r[1]*dx1-0.5303300858899105*rhou0_l[1]*uMax_l[1]*dx1-0.5303300858899105*rhou0_c[1]*uMax_l[1]*dx1-0.3061862178478971*rhou0_l[0]*uMax_l[1]*dx1+0.3061862178478971*rhou0_c[0]*uMax_l[1]*dx1-0.8660254037844386*urec_0_r[0]*Ghat_rho_r[1]*dx1-0.8660254037844386*urec_0_l[0]*Ghat_rho_l[1]*dx1; 

  outrhou1[0] += (-0.3061862178478971*uMax_r[1]*rhou1_r[3]*dx1)+0.3061862178478971*uMax_l[1]*rhou1_l[3]*dx1-0.3061862178478971*uMax_r[1]*rhou1_c[3]*dx1+0.3061862178478971*uMax_l[1]*rhou1_c[3]*dx1+0.1767766952966368*uMax_r[1]*rhou1_r[2]*dx1+0.1767766952966368*uMax_l[1]*rhou1_l[2]*dx1-0.1767766952966368*uMax_r[1]*rhou1_c[2]*dx1-0.1767766952966368*uMax_l[1]*rhou1_c[2]*dx1-0.5*Ghat_rho_r[1]*urec_1_r[1]*dx1+0.5*Ghat_rho_l[1]*urec_1_l[1]*dx1-0.3061862178478971*uMax_r[0]*rhou1_r[1]*dx1+0.3061862178478971*uMax_l[0]*rhou1_l[1]*dx1-0.3061862178478971*uMax_r[0]*rhou1_c[1]*dx1+0.3061862178478971*uMax_l[0]*rhou1_c[1]*dx1-0.5*Ghat_rho_r[0]*urec_1_r[0]*dx1+0.5*Ghat_rho_l[0]*urec_1_l[0]*dx1+0.1767766952966368*rhou1_r[0]*uMax_r[0]*dx1-0.1767766952966368*rhou1_c[0]*uMax_r[0]*dx1+0.1767766952966368*rhou1_l[0]*uMax_l[0]*dx1-0.1767766952966368*rhou1_c[0]*uMax_l[0]*dx1; 
  outrhou1[1] += (-0.5303300858899105*uMax_r[1]*rhou1_r[3]*dx1)-0.5303300858899105*uMax_l[1]*rhou1_l[3]*dx1-0.5303300858899105*uMax_r[1]*rhou1_c[3]*dx1-0.5303300858899105*uMax_l[1]*rhou1_c[3]*dx1+0.3061862178478971*uMax_r[1]*rhou1_r[2]*dx1-0.3061862178478971*uMax_l[1]*rhou1_l[2]*dx1-0.3061862178478971*uMax_r[1]*rhou1_c[2]*dx1+0.3061862178478971*uMax_l[1]*rhou1_c[2]*dx1-0.8660254037844386*Ghat_rho_r[1]*urec_1_r[1]*dx1-0.8660254037844386*Ghat_rho_l[1]*urec_1_l[1]*dx1-0.5303300858899105*uMax_r[0]*rhou1_r[1]*dx1-0.5303300858899105*uMax_l[0]*rhou1_l[1]*dx1-0.5303300858899105*uMax_r[0]*rhou1_c[1]*dx1-0.5303300858899105*uMax_l[0]*rhou1_c[1]*dx1-0.8660254037844386*Ghat_rho_r[0]*urec_1_r[0]*dx1-0.8660254037844386*Ghat_rho_l[0]*urec_1_l[0]*dx1+0.3061862178478971*rhou1_r[0]*uMax_r[0]*dx1-0.3061862178478971*rhou1_c[0]*uMax_r[0]*dx1-0.3061862178478971*rhou1_l[0]*uMax_l[0]*dx1+0.3061862178478971*rhou1_c[0]*uMax_l[0]*dx1; 
  outrhou1[2] += (-0.3061862178478971*uMax_r[0]*rhou1_r[3]*dx1)+0.3061862178478971*uMax_l[0]*rhou1_l[3]*dx1-0.3061862178478971*uMax_r[0]*rhou1_c[3]*dx1+0.3061862178478971*uMax_l[0]*rhou1_c[3]*dx1+0.1767766952966368*uMax_r[0]*rhou1_r[2]*dx1+0.1767766952966368*uMax_l[0]*rhou1_l[2]*dx1-0.1767766952966368*uMax_r[0]*rhou1_c[2]*dx1-0.1767766952966368*uMax_l[0]*rhou1_c[2]*dx1-0.5*Ghat_rho_r[0]*urec_1_r[1]*dx1+0.5*Ghat_rho_l[0]*urec_1_l[1]*dx1-0.3061862178478971*rhou1_r[1]*uMax_r[1]*dx1-0.3061862178478971*rhou1_c[1]*uMax_r[1]*dx1+0.1767766952966368*rhou1_r[0]*uMax_r[1]*dx1-0.1767766952966368*rhou1_c[0]*uMax_r[1]*dx1+0.3061862178478971*rhou1_l[1]*uMax_l[1]*dx1+0.3061862178478971*rhou1_c[1]*uMax_l[1]*dx1+0.1767766952966368*rhou1_l[0]*uMax_l[1]*dx1-0.1767766952966368*rhou1_c[0]*uMax_l[1]*dx1-0.5*urec_1_r[0]*Ghat_rho_r[1]*dx1+0.5*urec_1_l[0]*Ghat_rho_l[1]*dx1; 
  outrhou1[3] += (-0.5303300858899105*uMax_r[0]*rhou1_r[3]*dx1)-0.5303300858899105*uMax_l[0]*rhou1_l[3]*dx1-0.5303300858899105*uMax_r[0]*rhou1_c[3]*dx1-0.5303300858899105*uMax_l[0]*rhou1_c[3]*dx1+0.3061862178478971*uMax_r[0]*rhou1_r[2]*dx1-0.3061862178478971*uMax_l[0]*rhou1_l[2]*dx1-0.3061862178478971*uMax_r[0]*rhou1_c[2]*dx1+0.3061862178478971*uMax_l[0]*rhou1_c[2]*dx1-0.8660254037844386*Ghat_rho_r[0]*urec_1_r[1]*dx1-0.8660254037844386*Ghat_rho_l[0]*urec_1_l[1]*dx1-0.5303300858899105*rhou1_r[1]*uMax_r[1]*dx1-0.5303300858899105*rhou1_c[1]*uMax_r[1]*dx1+0.3061862178478971*rhou1_r[0]*uMax_r[1]*dx1-0.3061862178478971*rhou1_c[0]*uMax_r[1]*dx1-0.5303300858899105*rhou1_l[1]*uMax_l[1]*dx1-0.5303300858899105*rhou1_c[1]*uMax_l[1]*dx1-0.3061862178478971*rhou1_l[0]*uMax_l[1]*dx1+0.3061862178478971*rhou1_c[0]*uMax_l[1]*dx1-0.8660254037844386*urec_1_r[0]*Ghat_rho_r[1]*dx1-0.8660254037844386*urec_1_l[0]*Ghat_rho_l[1]*dx1; 

  outrhou2[0] += (-0.3061862178478971*uMax_r[1]*rhou2_r[3]*dx1)+0.3061862178478971*uMax_l[1]*rhou2_l[3]*dx1-0.3061862178478971*uMax_r[1]*rhou2_c[3]*dx1+0.3061862178478971*uMax_l[1]*rhou2_c[3]*dx1+0.1767766952966368*uMax_r[1]*rhou2_r[2]*dx1+0.1767766952966368*uMax_l[1]*rhou2_l[2]*dx1-0.1767766952966368*uMax_r[1]*rhou2_c[2]*dx1-0.1767766952966368*uMax_l[1]*rhou2_c[2]*dx1-0.5*Ghat_rho_r[1]*urec_2_r[1]*dx1+0.5*Ghat_rho_l[1]*urec_2_l[1]*dx1-0.3061862178478971*uMax_r[0]*rhou2_r[1]*dx1+0.3061862178478971*uMax_l[0]*rhou2_l[1]*dx1-0.3061862178478971*uMax_r[0]*rhou2_c[1]*dx1+0.3061862178478971*uMax_l[0]*rhou2_c[1]*dx1-0.5*Ghat_rho_r[0]*urec_2_r[0]*dx1+0.5*Ghat_rho_l[0]*urec_2_l[0]*dx1+0.1767766952966368*rhou2_r[0]*uMax_r[0]*dx1-0.1767766952966368*rhou2_c[0]*uMax_r[0]*dx1+0.1767766952966368*rhou2_l[0]*uMax_l[0]*dx1-0.1767766952966368*rhou2_c[0]*uMax_l[0]*dx1; 
  outrhou2[1] += (-0.5303300858899105*uMax_r[1]*rhou2_r[3]*dx1)-0.5303300858899105*uMax_l[1]*rhou2_l[3]*dx1-0.5303300858899105*uMax_r[1]*rhou2_c[3]*dx1-0.5303300858899105*uMax_l[1]*rhou2_c[3]*dx1+0.3061862178478971*uMax_r[1]*rhou2_r[2]*dx1-0.3061862178478971*uMax_l[1]*rhou2_l[2]*dx1-0.3061862178478971*uMax_r[1]*rhou2_c[2]*dx1+0.3061862178478971*uMax_l[1]*rhou2_c[2]*dx1-0.8660254037844386*Ghat_rho_r[1]*urec_2_r[1]*dx1-0.8660254037844386*Ghat_rho_l[1]*urec_2_l[1]*dx1-0.5303300858899105*uMax_r[0]*rhou2_r[1]*dx1-0.5303300858899105*uMax_l[0]*rhou2_l[1]*dx1-0.5303300858899105*uMax_r[0]*rhou2_c[1]*dx1-0.5303300858899105*uMax_l[0]*rhou2_c[1]*dx1-0.8660254037844386*Ghat_rho_r[0]*urec_2_r[0]*dx1-0.8660254037844386*Ghat_rho_l[0]*urec_2_l[0]*dx1+0.3061862178478971*rhou2_r[0]*uMax_r[0]*dx1-0.3061862178478971*rhou2_c[0]*uMax_r[0]*dx1-0.3061862178478971*rhou2_l[0]*uMax_l[0]*dx1+0.3061862178478971*rhou2_c[0]*uMax_l[0]*dx1; 
  outrhou2[2] += (-0.3061862178478971*uMax_r[0]*rhou2_r[3]*dx1)+0.3061862178478971*uMax_l[0]*rhou2_l[3]*dx1-0.3061862178478971*uMax_r[0]*rhou2_c[3]*dx1+0.3061862178478971*uMax_l[0]*rhou2_c[3]*dx1+0.1767766952966368*uMax_r[0]*rhou2_r[2]*dx1+0.1767766952966368*uMax_l[0]*rhou2_l[2]*dx1-0.1767766952966368*uMax_r[0]*rhou2_c[2]*dx1-0.1767766952966368*uMax_l[0]*rhou2_c[2]*dx1-0.5*Ghat_rho_r[0]*urec_2_r[1]*dx1+0.5*Ghat_rho_l[0]*urec_2_l[1]*dx1-0.3061862178478971*rhou2_r[1]*uMax_r[1]*dx1-0.3061862178478971*rhou2_c[1]*uMax_r[1]*dx1+0.1767766952966368*rhou2_r[0]*uMax_r[1]*dx1-0.1767766952966368*rhou2_c[0]*uMax_r[1]*dx1+0.3061862178478971*rhou2_l[1]*uMax_l[1]*dx1+0.3061862178478971*rhou2_c[1]*uMax_l[1]*dx1+0.1767766952966368*rhou2_l[0]*uMax_l[1]*dx1-0.1767766952966368*rhou2_c[0]*uMax_l[1]*dx1-0.5*urec_2_r[0]*Ghat_rho_r[1]*dx1+0.5*urec_2_l[0]*Ghat_rho_l[1]*dx1; 
  outrhou2[3] += (-0.5303300858899105*uMax_r[0]*rhou2_r[3]*dx1)-0.5303300858899105*uMax_l[0]*rhou2_l[3]*dx1-0.5303300858899105*uMax_r[0]*rhou2_c[3]*dx1-0.5303300858899105*uMax_l[0]*rhou2_c[3]*dx1+0.3061862178478971*uMax_r[0]*rhou2_r[2]*dx1-0.3061862178478971*uMax_l[0]*rhou2_l[2]*dx1-0.3061862178478971*uMax_r[0]*rhou2_c[2]*dx1+0.3061862178478971*uMax_l[0]*rhou2_c[2]*dx1-0.8660254037844386*Ghat_rho_r[0]*urec_2_r[1]*dx1-0.8660254037844386*Ghat_rho_l[0]*urec_2_l[1]*dx1-0.5303300858899105*rhou2_r[1]*uMax_r[1]*dx1-0.5303300858899105*rhou2_c[1]*uMax_r[1]*dx1+0.3061862178478971*rhou2_r[0]*uMax_r[1]*dx1-0.3061862178478971*rhou2_c[0]*uMax_r[1]*dx1-0.5303300858899105*rhou2_l[1]*uMax_l[1]*dx1-0.5303300858899105*rhou2_c[1]*uMax_l[1]*dx1-0.3061862178478971*rhou2_l[0]*uMax_l[1]*dx1+0.3061862178478971*rhou2_c[0]*uMax_l[1]*dx1-0.8660254037844386*urec_2_r[0]*Ghat_rho_r[1]*dx1-0.8660254037844386*urec_2_l[0]*Ghat_rho_l[1]*dx1; 

} 
