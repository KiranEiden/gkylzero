#include <gkyl_euler_kernels.h> 
#include <gkyl_basis_ser_1x_p2_surfx1_eval_quad.h> 
GKYL_CU_DH void euler_pkpm_surfx_1x_ser_p2(const double *w, const double *dxv, 
  const double *u_il, const double *u_ic, const double *u_ir,
  const double *p_ijl, const double *p_ijc, const double *p_ijr,
  const double *vlasov_pkpm_surf_momsl, const double *vlasov_pkpm_surf_momsr, 
  const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // u_il/u_ic/u_ir: [ux, uy, uz] Fluid flow in left/center/right cells.
  // p_ijl/p_ijc/p_ijr: Pressure tensor in left/center/right cells.
  // vlasov_pkpm_surf_momsl/vlasov_pkpm_surf_momsr: Mass flux and heat flux at left edge and right edge (computed externally) .
  // statevecl/statevecc/statevecr: [rho ux, rho uy, rho uz, energy], Fluid input state vector in left/center/right cells.
  // out: Incremented output.

  const double dx1 = 2.0/dxv[0]; 
  // Only need to fetch input energy variable, other fluxes are computed from input flow and pressure tensor.
  const double *energy_l = &statevecl[9]; 
  const double *energy_c = &statevecc[9]; 
  const double *energy_r = &statevecr[9]; 

  const double *ux_l = &u_il[0]; 
  const double *ux_c = &u_ic[0]; 
  const double *ux_r = &u_ir[0]; 

  const double *uy_l = &u_il[3]; 
  const double *uy_c = &u_ic[3]; 
  const double *uy_r = &u_ir[3]; 

  const double *uz_l = &u_il[6]; 
  const double *uz_c = &u_ic[6]; 
  const double *uz_r = &u_ir[6]; 

  const double *Pxx_l = &p_ijl[0]; 
  const double *Pxx_c = &p_ijc[0]; 
  const double *Pxx_r = &p_ijr[0]; 

  const double *Pxy_l = &p_ijl[3]; 
  const double *Pxy_c = &p_ijc[3]; 
  const double *Pxy_r = &p_ijr[3]; 

  const double *Pxz_l = &p_ijl[6]; 
  const double *Pxz_c = &p_ijc[6]; 
  const double *Pxz_r = &p_ijr[6]; 

  const double *rho_flux_l = &vlasov_pkpm_surf_momsl[0]; 
  const double *heat_flux_l = &vlasov_pkpm_surf_momsl[1]; 
  const double *rho_flux_r = &vlasov_pkpm_surf_momsr[0]; 
  const double *heat_flux_r = &vlasov_pkpm_surf_momsr[1]; 
  double *outrhou0 = &out[0]; 
  double *outrhou1 = &out[3]; 
  double *outrhou2 = &out[6]; 
  double *outenergy = &out[9]; 

  double ux_l_r = ser_1x_p2_surfx1_eval_quad_node_0_r(ux_l); 
  double ux_c_l = ser_1x_p2_surfx1_eval_quad_node_0_l(ux_c); 
  double ux_c_r = ser_1x_p2_surfx1_eval_quad_node_0_r(ux_c); 
  double ux_r_l = ser_1x_p2_surfx1_eval_quad_node_0_l(ux_r); 

  double uy_l_r = ser_1x_p2_surfx1_eval_quad_node_0_r(uy_l); 
  double uy_c_l = ser_1x_p2_surfx1_eval_quad_node_0_l(uy_c); 
  double uy_c_r = ser_1x_p2_surfx1_eval_quad_node_0_r(uy_c); 
  double uy_r_l = ser_1x_p2_surfx1_eval_quad_node_0_l(uy_r); 

  double uz_l_r = ser_1x_p2_surfx1_eval_quad_node_0_r(uz_l); 
  double uz_c_l = ser_1x_p2_surfx1_eval_quad_node_0_l(uz_c); 
  double uz_c_r = ser_1x_p2_surfx1_eval_quad_node_0_r(uz_c); 
  double uz_r_l = ser_1x_p2_surfx1_eval_quad_node_0_l(uz_r); 

  double ux_max_l = fmax(fabs(ux_l_r), fabs(ux_c_l)); 
  double ux_max_r = fmax(fabs(ux_c_r), fabs(ux_r_l)); 
  double uy_max_l = fmax(fabs(uy_l_r), fabs(uy_c_l)); 
  double uy_max_r = fmax(fabs(uy_c_r), fabs(uy_r_l)); 
  double uz_max_l = fmax(fabs(uz_l_r), fabs(uz_c_l)); 
  double uz_max_r = fmax(fabs(uz_c_r), fabs(uz_r_l)); 

  double Ghat_energy_l = (-0.7905694150420948*Pxz_l[2]*uz_max_l)+0.7905694150420948*Pxz_c[2]*uz_max_l-0.6123724356957945*Pxz_l[1]*uz_max_l-0.6123724356957945*Pxz_c[1]*uz_max_l-0.3535533905932737*Pxz_l[0]*uz_max_l+0.3535533905932737*Pxz_c[0]*uz_max_l-0.7905694150420948*Pxy_l[2]*uy_max_l+0.7905694150420948*Pxy_c[2]*uy_max_l-0.6123724356957945*Pxy_l[1]*uy_max_l-0.6123724356957945*Pxy_c[1]*uy_max_l-0.3535533905932737*Pxy_l[0]*uy_max_l+0.3535533905932737*Pxy_c[0]*uy_max_l-0.7905694150420948*energy_l[2]*ux_max_l+0.7905694150420948*energy_c[2]*ux_max_l-0.7905694150420948*Pxx_l[2]*ux_max_l+0.7905694150420948*Pxx_c[2]*ux_max_l-0.6123724356957945*energy_l[1]*ux_max_l-0.6123724356957945*energy_c[1]*ux_max_l-0.6123724356957945*Pxx_l[1]*ux_max_l-0.6123724356957945*Pxx_c[1]*ux_max_l-0.3535533905932737*energy_l[0]*ux_max_l+0.3535533905932737*energy_c[0]*ux_max_l-0.3535533905932737*Pxx_l[0]*ux_max_l+0.3535533905932737*Pxx_c[0]*ux_max_l+1.25*Pxz_l[2]*uz_l[2]+0.9682458365518543*Pxz_l[1]*uz_l[2]+0.5590169943749475*Pxz_l[0]*uz_l[2]+1.25*Pxz_c[2]*uz_c[2]-0.9682458365518543*Pxz_c[1]*uz_c[2]+0.5590169943749475*Pxz_c[0]*uz_c[2]+1.25*Pxy_l[2]*uy_l[2]+0.9682458365518543*Pxy_l[1]*uy_l[2]+0.5590169943749475*Pxy_l[0]*uy_l[2]+1.25*Pxy_c[2]*uy_c[2]-0.9682458365518543*Pxy_c[1]*uy_c[2]+0.5590169943749475*Pxy_c[0]*uy_c[2]+1.25*energy_l[2]*ux_l[2]+1.25*Pxx_l[2]*ux_l[2]+0.9682458365518543*energy_l[1]*ux_l[2]+0.9682458365518543*Pxx_l[1]*ux_l[2]+0.5590169943749475*energy_l[0]*ux_l[2]+0.5590169943749475*Pxx_l[0]*ux_l[2]+1.25*energy_c[2]*ux_c[2]+1.25*Pxx_c[2]*ux_c[2]-0.9682458365518543*energy_c[1]*ux_c[2]-0.9682458365518543*Pxx_c[1]*ux_c[2]+0.5590169943749475*energy_c[0]*ux_c[2]+0.5590169943749475*Pxx_c[0]*ux_c[2]+0.9682458365518543*ux_l[1]*energy_l[2]+0.5590169943749475*ux_l[0]*energy_l[2]-0.9682458365518543*ux_c[1]*energy_c[2]+0.5590169943749475*ux_c[0]*energy_c[2]+0.9682458365518543*uz_l[1]*Pxz_l[2]+0.5590169943749475*uz_l[0]*Pxz_l[2]-0.9682458365518543*uz_c[1]*Pxz_c[2]+0.5590169943749475*uz_c[0]*Pxz_c[2]+0.9682458365518543*uy_l[1]*Pxy_l[2]+0.5590169943749475*uy_l[0]*Pxy_l[2]-0.9682458365518543*uy_c[1]*Pxy_c[2]+0.5590169943749475*uy_c[0]*Pxy_c[2]+0.9682458365518543*ux_l[1]*Pxx_l[2]+0.5590169943749475*ux_l[0]*Pxx_l[2]-0.9682458365518543*ux_c[1]*Pxx_c[2]+0.5590169943749475*ux_c[0]*Pxx_c[2]+0.75*Pxz_l[1]*uz_l[1]+0.4330127018922193*Pxz_l[0]*uz_l[1]+0.75*Pxz_c[1]*uz_c[1]-0.4330127018922193*Pxz_c[0]*uz_c[1]+0.75*Pxy_l[1]*uy_l[1]+0.4330127018922193*Pxy_l[0]*uy_l[1]+0.75*Pxy_c[1]*uy_c[1]-0.4330127018922193*Pxy_c[0]*uy_c[1]+0.75*energy_l[1]*ux_l[1]+0.75*Pxx_l[1]*ux_l[1]+0.4330127018922193*energy_l[0]*ux_l[1]+0.4330127018922193*Pxx_l[0]*ux_l[1]+0.75*energy_c[1]*ux_c[1]+0.75*Pxx_c[1]*ux_c[1]-0.4330127018922193*energy_c[0]*ux_c[1]-0.4330127018922193*Pxx_c[0]*ux_c[1]+0.4330127018922193*ux_l[0]*energy_l[1]-0.4330127018922193*ux_c[0]*energy_c[1]+0.4330127018922193*uz_l[0]*Pxz_l[1]-0.4330127018922193*uz_c[0]*Pxz_c[1]+0.4330127018922193*uy_l[0]*Pxy_l[1]-0.4330127018922193*uy_c[0]*Pxy_c[1]+0.4330127018922193*ux_l[0]*Pxx_l[1]-0.4330127018922193*ux_c[0]*Pxx_c[1]+0.25*Pxz_l[0]*uz_l[0]+0.25*Pxz_c[0]*uz_c[0]+0.25*Pxy_l[0]*uy_l[0]+0.25*Pxy_c[0]*uy_c[0]+0.25*energy_l[0]*ux_l[0]+0.25*Pxx_l[0]*ux_l[0]+0.25*energy_c[0]*ux_c[0]+0.25*Pxx_c[0]*ux_c[0]; 
  double Ghat_energy_r = 0.7905694150420948*Pxz_r[2]*uz_max_r-0.7905694150420948*Pxz_c[2]*uz_max_r-0.6123724356957945*Pxz_r[1]*uz_max_r-0.6123724356957945*Pxz_c[1]*uz_max_r+0.3535533905932737*Pxz_r[0]*uz_max_r-0.3535533905932737*Pxz_c[0]*uz_max_r+0.7905694150420948*Pxy_r[2]*uy_max_r-0.7905694150420948*Pxy_c[2]*uy_max_r-0.6123724356957945*Pxy_r[1]*uy_max_r-0.6123724356957945*Pxy_c[1]*uy_max_r+0.3535533905932737*Pxy_r[0]*uy_max_r-0.3535533905932737*Pxy_c[0]*uy_max_r+0.7905694150420948*energy_r[2]*ux_max_r-0.7905694150420948*energy_c[2]*ux_max_r+0.7905694150420948*Pxx_r[2]*ux_max_r-0.7905694150420948*Pxx_c[2]*ux_max_r-0.6123724356957945*energy_r[1]*ux_max_r-0.6123724356957945*energy_c[1]*ux_max_r-0.6123724356957945*Pxx_r[1]*ux_max_r-0.6123724356957945*Pxx_c[1]*ux_max_r+0.3535533905932737*energy_r[0]*ux_max_r-0.3535533905932737*energy_c[0]*ux_max_r+0.3535533905932737*Pxx_r[0]*ux_max_r-0.3535533905932737*Pxx_c[0]*ux_max_r+1.25*Pxz_r[2]*uz_r[2]-0.9682458365518543*Pxz_r[1]*uz_r[2]+0.5590169943749475*Pxz_r[0]*uz_r[2]+1.25*Pxz_c[2]*uz_c[2]+0.9682458365518543*Pxz_c[1]*uz_c[2]+0.5590169943749475*Pxz_c[0]*uz_c[2]+1.25*Pxy_r[2]*uy_r[2]-0.9682458365518543*Pxy_r[1]*uy_r[2]+0.5590169943749475*Pxy_r[0]*uy_r[2]+1.25*Pxy_c[2]*uy_c[2]+0.9682458365518543*Pxy_c[1]*uy_c[2]+0.5590169943749475*Pxy_c[0]*uy_c[2]+1.25*energy_r[2]*ux_r[2]+1.25*Pxx_r[2]*ux_r[2]-0.9682458365518543*energy_r[1]*ux_r[2]-0.9682458365518543*Pxx_r[1]*ux_r[2]+0.5590169943749475*energy_r[0]*ux_r[2]+0.5590169943749475*Pxx_r[0]*ux_r[2]+1.25*energy_c[2]*ux_c[2]+1.25*Pxx_c[2]*ux_c[2]+0.9682458365518543*energy_c[1]*ux_c[2]+0.9682458365518543*Pxx_c[1]*ux_c[2]+0.5590169943749475*energy_c[0]*ux_c[2]+0.5590169943749475*Pxx_c[0]*ux_c[2]-0.9682458365518543*ux_r[1]*energy_r[2]+0.5590169943749475*ux_r[0]*energy_r[2]+0.9682458365518543*ux_c[1]*energy_c[2]+0.5590169943749475*ux_c[0]*energy_c[2]-0.9682458365518543*uz_r[1]*Pxz_r[2]+0.5590169943749475*uz_r[0]*Pxz_r[2]+0.9682458365518543*uz_c[1]*Pxz_c[2]+0.5590169943749475*uz_c[0]*Pxz_c[2]-0.9682458365518543*uy_r[1]*Pxy_r[2]+0.5590169943749475*uy_r[0]*Pxy_r[2]+0.9682458365518543*uy_c[1]*Pxy_c[2]+0.5590169943749475*uy_c[0]*Pxy_c[2]-0.9682458365518543*ux_r[1]*Pxx_r[2]+0.5590169943749475*ux_r[0]*Pxx_r[2]+0.9682458365518543*ux_c[1]*Pxx_c[2]+0.5590169943749475*ux_c[0]*Pxx_c[2]+0.75*Pxz_r[1]*uz_r[1]-0.4330127018922193*Pxz_r[0]*uz_r[1]+0.75*Pxz_c[1]*uz_c[1]+0.4330127018922193*Pxz_c[0]*uz_c[1]+0.75*Pxy_r[1]*uy_r[1]-0.4330127018922193*Pxy_r[0]*uy_r[1]+0.75*Pxy_c[1]*uy_c[1]+0.4330127018922193*Pxy_c[0]*uy_c[1]+0.75*energy_r[1]*ux_r[1]+0.75*Pxx_r[1]*ux_r[1]-0.4330127018922193*energy_r[0]*ux_r[1]-0.4330127018922193*Pxx_r[0]*ux_r[1]+0.75*energy_c[1]*ux_c[1]+0.75*Pxx_c[1]*ux_c[1]+0.4330127018922193*energy_c[0]*ux_c[1]+0.4330127018922193*Pxx_c[0]*ux_c[1]-0.4330127018922193*ux_r[0]*energy_r[1]+0.4330127018922193*ux_c[0]*energy_c[1]-0.4330127018922193*uz_r[0]*Pxz_r[1]+0.4330127018922193*uz_c[0]*Pxz_c[1]-0.4330127018922193*uy_r[0]*Pxy_r[1]+0.4330127018922193*uy_c[0]*Pxy_c[1]-0.4330127018922193*ux_r[0]*Pxx_r[1]+0.4330127018922193*ux_c[0]*Pxx_c[1]+0.25*Pxz_r[0]*uz_r[0]+0.25*Pxz_c[0]*uz_c[0]+0.25*Pxy_r[0]*uy_r[0]+0.25*Pxy_c[0]*uy_c[0]+0.25*energy_r[0]*ux_r[0]+0.25*Pxx_r[0]*ux_r[0]+0.25*energy_c[0]*ux_c[0]+0.25*Pxx_c[0]*ux_c[0]; 

  double uxrec_l = 0.3458741190809163*ux_l[2]+0.3458741190809163*ux_c[2]+0.4975526040028326*ux_l[1]-0.4975526040028326*ux_c[1]+0.3535533905932737*ux_l[0]+0.3535533905932737*ux_c[0]; 
  double uxrec_r = 0.3458741190809163*ux_r[2]+0.3458741190809163*ux_c[2]-0.4975526040028326*ux_r[1]+0.4975526040028326*ux_c[1]+0.3535533905932737*ux_r[0]+0.3535533905932737*ux_c[0]; 
  double uyrec_l = 0.3458741190809163*uy_l[2]+0.3458741190809163*uy_c[2]+0.4975526040028326*uy_l[1]-0.4975526040028326*uy_c[1]+0.3535533905932737*uy_l[0]+0.3535533905932737*uy_c[0]; 
  double uyrec_r = 0.3458741190809163*uy_r[2]+0.3458741190809163*uy_c[2]-0.4975526040028326*uy_r[1]+0.4975526040028326*uy_c[1]+0.3535533905932737*uy_r[0]+0.3535533905932737*uy_c[0]; 
  double uzrec_l = 0.3458741190809163*uz_l[2]+0.3458741190809163*uz_c[2]+0.4975526040028326*uz_l[1]-0.4975526040028326*uz_c[1]+0.3535533905932737*uz_l[0]+0.3535533905932737*uz_c[0]; 
  double uzrec_r = 0.3458741190809163*uz_r[2]+0.3458741190809163*uz_c[2]-0.4975526040028326*uz_r[1]+0.4975526040028326*uz_c[1]+0.3535533905932737*uz_r[0]+0.3535533905932737*uz_c[0]; 

  double Pxxrec_l = 0.3458741190809163*Pxx_l[2]+0.3458741190809163*Pxx_c[2]+0.4975526040028326*Pxx_l[1]-0.4975526040028326*Pxx_c[1]+0.3535533905932737*Pxx_l[0]+0.3535533905932737*Pxx_c[0]; 
  double Pxxrec_r = 0.3458741190809163*Pxx_r[2]+0.3458741190809163*Pxx_c[2]-0.4975526040028326*Pxx_r[1]+0.4975526040028326*Pxx_c[1]+0.3535533905932737*Pxx_r[0]+0.3535533905932737*Pxx_c[0]; 
  double Pxyrec_l = 0.3458741190809163*Pxy_l[2]+0.3458741190809163*Pxy_c[2]+0.4975526040028326*Pxy_l[1]-0.4975526040028326*Pxy_c[1]+0.3535533905932737*Pxy_l[0]+0.3535533905932737*Pxy_c[0]; 
  double Pxyrec_r = 0.3458741190809163*Pxy_r[2]+0.3458741190809163*Pxy_c[2]-0.4975526040028326*Pxy_r[1]+0.4975526040028326*Pxy_c[1]+0.3535533905932737*Pxy_r[0]+0.3535533905932737*Pxy_c[0]; 
  double Pxzrec_l = 0.3458741190809163*Pxz_l[2]+0.3458741190809163*Pxz_c[2]+0.4975526040028326*Pxz_l[1]-0.4975526040028326*Pxz_c[1]+0.3535533905932737*Pxz_l[0]+0.3535533905932737*Pxz_c[0]; 
  double Pxzrec_r = 0.3458741190809163*Pxz_r[2]+0.3458741190809163*Pxz_c[2]-0.4975526040028326*Pxz_r[1]+0.4975526040028326*Pxz_c[1]+0.3535533905932737*Pxz_r[0]+0.3535533905932737*Pxz_c[0]; 

  outrhou0[0] += (-0.7071067811865475*rho_flux_r[0]*dx1*uxrec_r)+0.7071067811865475*rho_flux_l[0]*dx1*uxrec_l-0.7071067811865475*Pxxrec_r*dx1+0.7071067811865475*Pxxrec_l*dx1; 
  outrhou0[1] += (-1.224744871391589*rho_flux_r[0]*dx1*uxrec_r)-1.224744871391589*rho_flux_l[0]*dx1*uxrec_l-1.224744871391589*Pxxrec_r*dx1-1.224744871391589*Pxxrec_l*dx1; 
  outrhou0[2] += (-1.58113883008419*rho_flux_r[0]*dx1*uxrec_r)+1.58113883008419*rho_flux_l[0]*dx1*uxrec_l-1.58113883008419*Pxxrec_r*dx1+1.58113883008419*Pxxrec_l*dx1; 

  outrhou1[0] += (-0.7071067811865475*rho_flux_r[0]*dx1*uyrec_r)+0.7071067811865475*rho_flux_l[0]*dx1*uyrec_l-0.7071067811865475*Pxyrec_r*dx1+0.7071067811865475*Pxyrec_l*dx1; 
  outrhou1[1] += (-1.224744871391589*rho_flux_r[0]*dx1*uyrec_r)-1.224744871391589*rho_flux_l[0]*dx1*uyrec_l-1.224744871391589*Pxyrec_r*dx1-1.224744871391589*Pxyrec_l*dx1; 
  outrhou1[2] += (-1.58113883008419*rho_flux_r[0]*dx1*uyrec_r)+1.58113883008419*rho_flux_l[0]*dx1*uyrec_l-1.58113883008419*Pxyrec_r*dx1+1.58113883008419*Pxyrec_l*dx1; 

  outrhou2[0] += (-0.7071067811865475*rho_flux_r[0]*dx1*uzrec_r)+0.7071067811865475*rho_flux_l[0]*dx1*uzrec_l-0.7071067811865475*Pxzrec_r*dx1+0.7071067811865475*Pxyrec_l*dx1; 
  outrhou2[1] += (-1.224744871391589*rho_flux_r[0]*dx1*uzrec_r)-1.224744871391589*rho_flux_l[0]*dx1*uzrec_l-1.224744871391589*Pxzrec_r*dx1-1.224744871391589*Pxyrec_l*dx1; 
  outrhou2[2] += (-1.58113883008419*rho_flux_r[0]*dx1*uzrec_r)+1.58113883008419*rho_flux_l[0]*dx1*uzrec_l-1.58113883008419*Pxzrec_r*dx1+1.58113883008419*Pxyrec_l*dx1; 

  outenergy[0] += (-0.7071067811865475*Ghat_energy_r*dx1)+0.7071067811865475*Ghat_energy_l*dx1-0.7071067811865475*heat_flux_r[0]*dx1+0.7071067811865475*heat_flux_l[0]*dx1; 
  outenergy[1] += (-1.224744871391589*Ghat_energy_r*dx1)-1.224744871391589*Ghat_energy_l*dx1-1.224744871391589*heat_flux_r[0]*dx1-1.224744871391589*heat_flux_l[0]*dx1; 
  outenergy[2] += (-1.58113883008419*Ghat_energy_r*dx1)+1.58113883008419*Ghat_energy_l*dx1-1.58113883008419*heat_flux_r[0]*dx1+1.58113883008419*heat_flux_l[0]*dx1; 

} 
