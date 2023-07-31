#include <gkyl_euler_kernels.h> 
#include <gkyl_basis_ser_1x_p1_surfx1_eval_quad.h> 
GKYL_CU_DH double euler_iso_surfx_1x_ser_p1(const double *w, const double *dxv, const double vth, const double *ul, const double *uc, const double *ur, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // gas_gamma: Adiabatic index.
  // ul/uc/ur: [ux, uy, uz] Fluid flow in left/center/right cells.
  // pl/pc/pr: Fluid pressure in left/center/right cells.
  // statevecl/statevecc/statevecr: [rho, rho ux, rho uy, rho uz], Fluid input state vector in left/center/right cells.
  // out: Incremented output.
  const double dx1 = 2.0/dxv[0]; 
  const double *rho_l = &statevecl[0]; 
  const double *rhou0_l = &statevecl[2]; 
  const double *rhou1_l = &statevecl[4]; 
  const double *rhou2_l = &statevecl[6]; 

  const double *rho_c = &statevecc[0]; 
  const double *rhou0_c = &statevecc[2]; 
  const double *rhou1_c = &statevecc[4]; 
  const double *rhou2_c = &statevecc[6]; 

  const double *rho_r = &statevecr[0]; 
  const double *rhou0_r = &statevecr[2]; 
  const double *rhou1_r = &statevecr[4]; 
  const double *rhou2_r = &statevecr[6]; 

  const double *ul_0 = &ul[0]; 
  const double *uc_0 = &uc[0]; 
  const double *ur_0 = &ur[0]; 

  const double *ul_1 = &ul[2]; 
  const double *uc_1 = &uc[2]; 
  const double *ur_1 = &ur[2]; 

  const double *ul_2 = &ul[4]; 
  const double *uc_2 = &uc[4]; 
  const double *ur_2 = &ur[4]; 

  double *outrho = &out[0]; 
  double *outrhou0 = &out[2]; 
  double *outrhou1 = &out[4]; 
  double *outrhou2 = &out[6]; 

  double vthsq = vth*vth; 
  double u_l_r = ser_1x_p1_surfx1_eval_quad_node_0_r(ul_0); 
  double u_c_l = ser_1x_p1_surfx1_eval_quad_node_0_l(uc_0); 
  double u_c_r = ser_1x_p1_surfx1_eval_quad_node_0_r(uc_0); 
  double u_r_l = ser_1x_p1_surfx1_eval_quad_node_0_l(ur_0); 
  double u_max_l = fmax(fabs(u_l_r), fabs(u_c_l)); 
  double u_max_r = fmax(fabs(u_c_r), fabs(u_r_l)); 

  double Ghat_rho_l = 0.6123724356957945*rho_l[1]*vth+0.6123724356957945*rho_c[1]*vth+0.3535533905932737*rho_l[0]*vth-0.3535533905932737*rho_c[0]*vth+0.6123724356957945*rho_l[1]*u_max_l+0.6123724356957945*rho_c[1]*u_max_l+0.3535533905932737*rho_l[0]*u_max_l-0.3535533905932737*rho_c[0]*u_max_l+0.6123724356957945*rhou0_l[1]-0.6123724356957945*rhou0_c[1]+0.3535533905932737*rhou0_l[0]+0.3535533905932737*rhou0_c[0]; 
  double Ghat_rho_r = 0.6123724356957945*rho_r[1]*vth+0.6123724356957945*rho_c[1]*vth-0.3535533905932737*rho_r[0]*vth+0.3535533905932737*rho_c[0]*vth+0.6123724356957945*rho_r[1]*u_max_r+0.6123724356957945*rho_c[1]*u_max_r-0.3535533905932737*rho_r[0]*u_max_r+0.3535533905932737*rho_c[0]*u_max_r-0.6123724356957945*rhou0_r[1]+0.6123724356957945*rhou0_c[1]+0.3535533905932737*rhou0_r[0]+0.3535533905932737*rhou0_c[0]; 
  double Ghat_rhoux_l = 0.6123724356957945*rho_l[1]*vthsq-0.6123724356957945*rho_c[1]*vthsq+0.3535533905932737*rho_l[0]*vthsq+0.3535533905932737*rho_c[0]*vthsq+0.6123724356957945*rhou0_l[1]*vth+0.6123724356957945*rhou0_c[1]*vth+0.3535533905932737*rhou0_l[0]*vth-0.3535533905932737*rhou0_c[0]*vth+0.6123724356957945*rhou0_l[1]*u_max_l+0.6123724356957945*rhou0_c[1]*u_max_l+0.3535533905932737*rhou0_l[0]*u_max_l-0.3535533905932737*rhou0_c[0]*u_max_l+0.75*rhou0_l[1]*ul_2[1]+0.4330127018922193*rhou0_l[0]*ul_2[1]+0.75*rhou0_l[1]*ul_1[1]+0.4330127018922193*rhou0_l[0]*ul_1[1]+0.75*rhou0_l[1]*ul_0[1]+0.4330127018922193*rhou0_l[0]*ul_0[1]+0.75*rhou0_c[1]*uc_2[1]-0.4330127018922193*rhou0_c[0]*uc_2[1]+0.75*rhou0_c[1]*uc_1[1]-0.4330127018922193*rhou0_c[0]*uc_1[1]+0.75*rhou0_c[1]*uc_0[1]-0.4330127018922193*rhou0_c[0]*uc_0[1]+0.4330127018922193*ul_2[0]*rhou0_l[1]+0.4330127018922193*ul_1[0]*rhou0_l[1]+0.4330127018922193*ul_0[0]*rhou0_l[1]-0.4330127018922193*uc_2[0]*rhou0_c[1]-0.4330127018922193*uc_1[0]*rhou0_c[1]-0.4330127018922193*uc_0[0]*rhou0_c[1]+0.25*rhou0_l[0]*ul_2[0]+0.25*rhou0_l[0]*ul_1[0]+0.25*rhou0_l[0]*ul_0[0]+0.25*rhou0_c[0]*uc_2[0]+0.25*rhou0_c[0]*uc_1[0]+0.25*rhou0_c[0]*uc_0[0]; 
  double Ghat_rhoux_r = (-0.6123724356957945*rho_r[1]*vthsq)+0.6123724356957945*rho_c[1]*vthsq+0.3535533905932737*rho_r[0]*vthsq+0.3535533905932737*rho_c[0]*vthsq+0.6123724356957945*rhou0_r[1]*vth+0.6123724356957945*rhou0_c[1]*vth-0.3535533905932737*rhou0_r[0]*vth+0.3535533905932737*rhou0_c[0]*vth+0.6123724356957945*rhou0_r[1]*u_max_r+0.6123724356957945*rhou0_c[1]*u_max_r-0.3535533905932737*rhou0_r[0]*u_max_r+0.3535533905932737*rhou0_c[0]*u_max_r+0.75*rhou0_r[1]*ur_2[1]-0.4330127018922193*rhou0_r[0]*ur_2[1]+0.75*rhou0_r[1]*ur_1[1]-0.4330127018922193*rhou0_r[0]*ur_1[1]+0.75*rhou0_r[1]*ur_0[1]-0.4330127018922193*rhou0_r[0]*ur_0[1]+0.75*rhou0_c[1]*uc_2[1]+0.4330127018922193*rhou0_c[0]*uc_2[1]+0.75*rhou0_c[1]*uc_1[1]+0.4330127018922193*rhou0_c[0]*uc_1[1]+0.75*rhou0_c[1]*uc_0[1]+0.4330127018922193*rhou0_c[0]*uc_0[1]-0.4330127018922193*ur_2[0]*rhou0_r[1]-0.4330127018922193*ur_1[0]*rhou0_r[1]-0.4330127018922193*ur_0[0]*rhou0_r[1]+0.4330127018922193*uc_2[0]*rhou0_c[1]+0.4330127018922193*uc_1[0]*rhou0_c[1]+0.4330127018922193*uc_0[0]*rhou0_c[1]+0.25*rhou0_r[0]*ur_2[0]+0.25*rhou0_r[0]*ur_1[0]+0.25*rhou0_r[0]*ur_0[0]+0.25*rhou0_c[0]*uc_2[0]+0.25*rhou0_c[0]*uc_1[0]+0.25*rhou0_c[0]*uc_0[0]; 
  double Ghat_rhouy_l = 0.6123724356957945*rhou1_l[1]*vth+0.6123724356957945*rhou1_c[1]*vth+0.3535533905932737*rhou1_l[0]*vth-0.3535533905932737*rhou1_c[0]*vth+0.6123724356957945*rhou1_l[1]*u_max_l+0.6123724356957945*rhou1_c[1]*u_max_l+0.3535533905932737*rhou1_l[0]*u_max_l-0.3535533905932737*rhou1_c[0]*u_max_l+0.75*rhou1_l[1]*ul_2[1]+0.4330127018922193*rhou1_l[0]*ul_2[1]+0.75*rhou1_l[1]*ul_1[1]+0.4330127018922193*rhou1_l[0]*ul_1[1]+0.75*rhou1_l[1]*ul_0[1]+0.4330127018922193*rhou1_l[0]*ul_0[1]+0.75*rhou1_c[1]*uc_2[1]-0.4330127018922193*rhou1_c[0]*uc_2[1]+0.75*rhou1_c[1]*uc_1[1]-0.4330127018922193*rhou1_c[0]*uc_1[1]+0.75*rhou1_c[1]*uc_0[1]-0.4330127018922193*rhou1_c[0]*uc_0[1]+0.4330127018922193*ul_2[0]*rhou1_l[1]+0.4330127018922193*ul_1[0]*rhou1_l[1]+0.4330127018922193*ul_0[0]*rhou1_l[1]-0.4330127018922193*uc_2[0]*rhou1_c[1]-0.4330127018922193*uc_1[0]*rhou1_c[1]-0.4330127018922193*uc_0[0]*rhou1_c[1]+0.25*rhou1_l[0]*ul_2[0]+0.25*rhou1_l[0]*ul_1[0]+0.25*rhou1_l[0]*ul_0[0]+0.25*rhou1_c[0]*uc_2[0]+0.25*rhou1_c[0]*uc_1[0]+0.25*rhou1_c[0]*uc_0[0]; 
  double Ghat_rhouy_r = 0.6123724356957945*rhou1_r[1]*vth+0.6123724356957945*rhou1_c[1]*vth-0.3535533905932737*rhou1_r[0]*vth+0.3535533905932737*rhou1_c[0]*vth+0.6123724356957945*rhou1_r[1]*u_max_r+0.6123724356957945*rhou1_c[1]*u_max_r-0.3535533905932737*rhou1_r[0]*u_max_r+0.3535533905932737*rhou1_c[0]*u_max_r+0.75*rhou1_r[1]*ur_2[1]-0.4330127018922193*rhou1_r[0]*ur_2[1]+0.75*rhou1_r[1]*ur_1[1]-0.4330127018922193*rhou1_r[0]*ur_1[1]+0.75*rhou1_r[1]*ur_0[1]-0.4330127018922193*rhou1_r[0]*ur_0[1]+0.75*rhou1_c[1]*uc_2[1]+0.4330127018922193*rhou1_c[0]*uc_2[1]+0.75*rhou1_c[1]*uc_1[1]+0.4330127018922193*rhou1_c[0]*uc_1[1]+0.75*rhou1_c[1]*uc_0[1]+0.4330127018922193*rhou1_c[0]*uc_0[1]-0.4330127018922193*ur_2[0]*rhou1_r[1]-0.4330127018922193*ur_1[0]*rhou1_r[1]-0.4330127018922193*ur_0[0]*rhou1_r[1]+0.4330127018922193*uc_2[0]*rhou1_c[1]+0.4330127018922193*uc_1[0]*rhou1_c[1]+0.4330127018922193*uc_0[0]*rhou1_c[1]+0.25*rhou1_r[0]*ur_2[0]+0.25*rhou1_r[0]*ur_1[0]+0.25*rhou1_r[0]*ur_0[0]+0.25*rhou1_c[0]*uc_2[0]+0.25*rhou1_c[0]*uc_1[0]+0.25*rhou1_c[0]*uc_0[0]; 
  double Ghat_rhouz_l = 0.6123724356957945*rhou2_l[1]*vth+0.6123724356957945*rhou2_c[1]*vth+0.3535533905932737*rhou2_l[0]*vth-0.3535533905932737*rhou2_c[0]*vth+0.6123724356957945*rhou2_l[1]*u_max_l+0.6123724356957945*rhou2_c[1]*u_max_l+0.3535533905932737*rhou2_l[0]*u_max_l-0.3535533905932737*rhou2_c[0]*u_max_l+0.75*rhou0_l[1]*ul_2[1]+0.4330127018922193*rhou0_l[0]*ul_2[1]+0.75*rhou0_l[1]*ul_1[1]+0.4330127018922193*rhou0_l[0]*ul_1[1]+0.75*rhou0_l[1]*ul_0[1]+0.4330127018922193*rhou0_l[0]*ul_0[1]+0.75*rhou0_c[1]*uc_2[1]-0.4330127018922193*rhou0_c[0]*uc_2[1]+0.75*rhou0_c[1]*uc_1[1]-0.4330127018922193*rhou0_c[0]*uc_1[1]+0.75*rhou0_c[1]*uc_0[1]-0.4330127018922193*rhou0_c[0]*uc_0[1]+0.4330127018922193*ul_2[0]*rhou0_l[1]+0.4330127018922193*ul_1[0]*rhou0_l[1]+0.4330127018922193*ul_0[0]*rhou0_l[1]-0.4330127018922193*uc_2[0]*rhou0_c[1]-0.4330127018922193*uc_1[0]*rhou0_c[1]-0.4330127018922193*uc_0[0]*rhou0_c[1]+0.25*rhou0_l[0]*ul_2[0]+0.25*rhou0_l[0]*ul_1[0]+0.25*rhou0_l[0]*ul_0[0]+0.25*rhou0_c[0]*uc_2[0]+0.25*rhou0_c[0]*uc_1[0]+0.25*rhou0_c[0]*uc_0[0]; 
  double Ghat_rhouz_r = 0.6123724356957945*rhou2_r[1]*vth+0.6123724356957945*rhou2_c[1]*vth-0.3535533905932737*rhou2_r[0]*vth+0.3535533905932737*rhou2_c[0]*vth+0.6123724356957945*rhou2_r[1]*u_max_r+0.6123724356957945*rhou2_c[1]*u_max_r-0.3535533905932737*rhou2_r[0]*u_max_r+0.3535533905932737*rhou2_c[0]*u_max_r+0.75*rhou0_r[1]*ur_2[1]-0.4330127018922193*rhou0_r[0]*ur_2[1]+0.75*rhou0_r[1]*ur_1[1]-0.4330127018922193*rhou0_r[0]*ur_1[1]+0.75*rhou0_r[1]*ur_0[1]-0.4330127018922193*rhou0_r[0]*ur_0[1]+0.75*rhou0_c[1]*uc_2[1]+0.4330127018922193*rhou0_c[0]*uc_2[1]+0.75*rhou0_c[1]*uc_1[1]+0.4330127018922193*rhou0_c[0]*uc_1[1]+0.75*rhou0_c[1]*uc_0[1]+0.4330127018922193*rhou0_c[0]*uc_0[1]-0.4330127018922193*ur_2[0]*rhou0_r[1]-0.4330127018922193*ur_1[0]*rhou0_r[1]-0.4330127018922193*ur_0[0]*rhou0_r[1]+0.4330127018922193*uc_2[0]*rhou0_c[1]+0.4330127018922193*uc_1[0]*rhou0_c[1]+0.4330127018922193*uc_0[0]*rhou0_c[1]+0.25*rhou0_r[0]*ur_2[0]+0.25*rhou0_r[0]*ur_1[0]+0.25*rhou0_r[0]*ur_0[0]+0.25*rhou0_c[0]*uc_2[0]+0.25*rhou0_c[0]*uc_1[0]+0.25*rhou0_c[0]*uc_0[0]; 

  outrho[0] += 0.7071067811865475*Ghat_rho_l*dx1-0.7071067811865475*Ghat_rho_r*dx1; 
  outrho[1] += (-1.224744871391589*Ghat_rho_r*dx1)-1.224744871391589*Ghat_rho_l*dx1; 

  outrhou0[0] += 0.7071067811865475*Ghat_rhoux_l*dx1-0.7071067811865475*Ghat_rhoux_r*dx1; 
  outrhou0[1] += (-1.224744871391589*Ghat_rhoux_r*dx1)-1.224744871391589*Ghat_rhoux_l*dx1; 

  outrhou1[0] += 0.7071067811865475*Ghat_rhouy_l*dx1-0.7071067811865475*Ghat_rhouy_r*dx1; 
  outrhou1[1] += (-1.224744871391589*Ghat_rhouy_r*dx1)-1.224744871391589*Ghat_rhouy_l*dx1; 

  outrhou2[0] += 0.7071067811865475*Ghat_rhouz_l*dx1-0.7071067811865475*Ghat_rhouz_r*dx1; 
  outrhou2[1] += (-1.224744871391589*Ghat_rhouz_r*dx1)-1.224744871391589*Ghat_rhouz_l*dx1; 

  return 0.;

} 
