#include <gkyl_vlasov_pkpm_kernels.h> 
GKYL_CU_DH double vlasov_pkpm_surfx_1x1v_ser_p1(const double *w, const double *dxv, 
    const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
    const double *pkpm_prim_surf_l, const double *pkpm_prim_surf_c, const double *pkpm_prim_surf_r, 
    const double *fl, const double *fc, const double *fr, 
    const double *pkpm_max_b, const double *pkpm_lax, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:              Cell-center coordinates.
  // dxv[NDIM]:            Cell spacing.
  // bvar_surf_l/c/r:      Input surface magnetic field unit vector and tensor in left/center/right cells in each direction.
  // pkpm_prim_surf_l/c/r: Input surface primitive variables [u_i, 3*T_ii/m] in left/center/right cells in each direction.
  // fl/fc/fr:             Input distribution functions [F_0, T_perp/m G = T_perp/m (F_0 - F_1)] in left/center/right cells.
  // pkpm_max_b:           Surface expansion of max |b| for Lax penalization of streaming: lambda_i = |b_i|.
  // pkpm_lax:             Surface expansion of pkpm Lax penalization: lambda_i = |u_i| + sqrt(3.0*T_ii/m).
  // out:                  Incremented output distribution functions in center cell.
  const double dx1 = 2.0/dxv[0]; 
  const double dvpar = dxv[1], wvpar = w[1]; 
  const double *F_0l = &fl[0]; 
  const double *G_1l = &fl[6]; 
  const double *F_0c = &fc[0]; 
  const double *G_1c = &fc[6]; 
  const double *F_0r = &fr[0]; 
  const double *G_1r = &fr[6]; 
  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[6]; 
  double F_0_lr[3] = {0.0}; 
  double F_0_cl[3] = {0.0}; 
  double F_0_cr[3] = {0.0}; 
  double F_0_rl[3] = {0.0}; 
  double Ghat_F_0_vpar_l[3] = {0.0}; 
  double Ghat_F_0_vpar_r[3] = {0.0}; 
  double Ghat_F_0_u_l[3] = {0.0}; 
  double Ghat_F_0_u_r[3] = {0.0}; 

  double G_1_lr[3] = {0.0}; 
  double G_1_cl[3] = {0.0}; 
  double G_1_cr[3] = {0.0}; 
  double G_1_rl[3] = {0.0}; 
  double Ghat_G_1_vpar_l[3] = {0.0}; 
  double Ghat_G_1_vpar_r[3] = {0.0}; 
  double Ghat_G_1_u_l[3] = {0.0}; 
  double Ghat_G_1_u_r[3] = {0.0}; 

  F_0_lr[0] = 1.224744871391589*F_0l[1]+0.7071067811865475*F_0l[0]; 
  F_0_lr[1] = 1.224744871391589*F_0l[3]+0.7071067811865475*F_0l[2]; 
  F_0_lr[2] = 1.224744871391589*F_0l[5]+0.7071067811865475*F_0l[4]; 

  F_0_cl[0] = 0.7071067811865475*F_0c[0]-1.224744871391589*F_0c[1]; 
  F_0_cl[1] = 0.7071067811865475*F_0c[2]-1.224744871391589*F_0c[3]; 
  F_0_cl[2] = 0.7071067811865475*F_0c[4]-1.224744871391589*F_0c[5]; 

  F_0_cr[0] = 1.224744871391589*F_0c[1]+0.7071067811865475*F_0c[0]; 
  F_0_cr[1] = 1.224744871391589*F_0c[3]+0.7071067811865475*F_0c[2]; 
  F_0_cr[2] = 1.224744871391589*F_0c[5]+0.7071067811865475*F_0c[4]; 

  F_0_rl[0] = 0.7071067811865475*F_0r[0]-1.224744871391589*F_0r[1]; 
  F_0_rl[1] = 0.7071067811865475*F_0r[2]-1.224744871391589*F_0r[3]; 
  F_0_rl[2] = 0.7071067811865475*F_0r[4]-1.224744871391589*F_0r[5]; 

  G_1_lr[0] = 1.224744871391589*G_1l[1]+0.7071067811865475*G_1l[0]; 
  G_1_lr[1] = 1.224744871391589*G_1l[3]+0.7071067811865475*G_1l[2]; 
  G_1_lr[2] = 1.224744871391589*G_1l[5]+0.7071067811865475*G_1l[4]; 

  G_1_cl[0] = 0.7071067811865475*G_1c[0]-1.224744871391589*G_1c[1]; 
  G_1_cl[1] = 0.7071067811865475*G_1c[2]-1.224744871391589*G_1c[3]; 
  G_1_cl[2] = 0.7071067811865475*G_1c[4]-1.224744871391589*G_1c[5]; 

  G_1_cr[0] = 1.224744871391589*G_1c[1]+0.7071067811865475*G_1c[0]; 
  G_1_cr[1] = 1.224744871391589*G_1c[3]+0.7071067811865475*G_1c[2]; 
  G_1_cr[2] = 1.224744871391589*G_1c[5]+0.7071067811865475*G_1c[4]; 

  G_1_rl[0] = 0.7071067811865475*G_1r[0]-1.224744871391589*G_1r[1]; 
  G_1_rl[1] = 0.7071067811865475*G_1r[2]-1.224744871391589*G_1r[3]; 
  G_1_rl[2] = 0.7071067811865475*G_1r[4]-1.224744871391589*G_1r[5]; 

  const double *b_surf_lr = &bvar_surf_l[1]; 
  const double *b_surf_cl = &bvar_surf_c[0]; 
  const double *b_surf_cr = &bvar_surf_c[1]; 
  const double *b_surf_rl = &bvar_surf_r[0]; 

  const double *u_surf_lr = &pkpm_prim_surf_l[1]; 
  const double *u_surf_cl = &pkpm_prim_surf_c[0]; 
  const double *u_surf_cr = &pkpm_prim_surf_c[1]; 
  const double *u_surf_rl = &pkpm_prim_surf_r[0]; 

  const double *pkpm_max_b_l = &pkpm_max_b[0]; 
  const double *pkpm_max_b_r = &pkpm_max_b[1]; 

  const double *pkpm_lax_l = &pkpm_lax[0]; 
  const double *pkpm_lax_r = &pkpm_lax[1]; 

  double bl_r = b_surf_lr[0]; 
  double bc_l = b_surf_cl[0]; 
  double bc_r = b_surf_cr[0]; 
  double br_l = b_surf_rl[0]; 
  double avg_b_l = 0.5*(bl_r + bc_l); 
  double avg_b_r = 0.5*(bc_r + br_l); 

  double max_b_l = pkpm_max_b_l[0]; 
  double max_b_r = pkpm_max_b_r[0]; 

  if (wvpar>0) { 

  Ghat_F_0_vpar_l[0] = ((0.5*F_0_lr[0]-0.5*F_0_cl[0])*max_b_l+0.5*(F_0_lr[0]+F_0_cl[0])*avg_b_l)*wvpar+dvpar*((0.1443375672974065*F_0_lr[1]-0.1443375672974065*F_0_cl[1])*max_b_l+0.1443375672974065*(F_0_lr[1]+F_0_cl[1])*avg_b_l); 
  Ghat_F_0_vpar_l[1] = ((0.5*F_0_lr[1]-0.5*F_0_cl[1])*max_b_l+0.5*(F_0_lr[1]+F_0_cl[1])*avg_b_l)*wvpar+dvpar*((0.1290994448735806*F_0_lr[2]-0.1290994448735806*F_0_cl[2]+0.1443375672974065*F_0_lr[0]-0.1443375672974065*F_0_cl[0])*max_b_l+(0.1290994448735806*(F_0_lr[2]+F_0_cl[2])+0.1443375672974065*(F_0_lr[0]+F_0_cl[0]))*avg_b_l); 
  Ghat_F_0_vpar_l[2] = ((0.5*F_0_lr[2]-0.5*F_0_cl[2])*max_b_l+0.5*(F_0_lr[2]+F_0_cl[2])*avg_b_l)*wvpar+dvpar*((0.1290994448735806*F_0_lr[1]-0.1290994448735806*F_0_cl[1])*max_b_l+0.1290994448735806*(F_0_lr[1]+F_0_cl[1])*avg_b_l); 
  Ghat_G_1_vpar_l[0] = ((0.5*G_1_lr[0]-0.5*G_1_cl[0])*max_b_l+0.5*(G_1_lr[0]+G_1_cl[0])*avg_b_l)*wvpar+dvpar*((0.1443375672974065*G_1_lr[1]-0.1443375672974065*G_1_cl[1])*max_b_l+0.1443375672974065*(G_1_lr[1]+G_1_cl[1])*avg_b_l); 
  Ghat_G_1_vpar_l[1] = ((0.5*G_1_lr[1]-0.5*G_1_cl[1])*max_b_l+0.5*(G_1_lr[1]+G_1_cl[1])*avg_b_l)*wvpar+dvpar*((0.1290994448735806*G_1_lr[2]-0.1290994448735806*G_1_cl[2]+0.1443375672974065*G_1_lr[0]-0.1443375672974065*G_1_cl[0])*max_b_l+(0.1290994448735806*(G_1_lr[2]+G_1_cl[2])+0.1443375672974065*(G_1_lr[0]+G_1_cl[0]))*avg_b_l); 
  Ghat_G_1_vpar_l[2] = ((0.5*G_1_lr[2]-0.5*G_1_cl[2])*max_b_l+0.5*(G_1_lr[2]+G_1_cl[2])*avg_b_l)*wvpar+dvpar*((0.1290994448735806*G_1_lr[1]-0.1290994448735806*G_1_cl[1])*max_b_l+0.1290994448735806*(G_1_lr[1]+G_1_cl[1])*avg_b_l); 

  Ghat_F_0_vpar_r[0] = (0.5*(F_0_cr[0]*max_b_r+(F_0_rl[0]+F_0_cr[0])*avg_b_r)-0.5*F_0_rl[0]*max_b_r)*wvpar+dvpar*(0.1443375672974065*(F_0_cr[1]*max_b_r+(F_0_rl[1]+F_0_cr[1])*avg_b_r)-0.1443375672974065*F_0_rl[1]*max_b_r); 
  Ghat_F_0_vpar_r[1] = (0.5*(F_0_cr[1]*max_b_r+(F_0_rl[1]+F_0_cr[1])*avg_b_r)-0.5*F_0_rl[1]*max_b_r)*wvpar+dvpar*(((-0.1290994448735806*F_0_rl[2])+0.1290994448735806*F_0_cr[2]-0.1443375672974065*F_0_rl[0]+0.1443375672974065*F_0_cr[0])*max_b_r+(0.1290994448735806*(F_0_rl[2]+F_0_cr[2])+0.1443375672974065*(F_0_rl[0]+F_0_cr[0]))*avg_b_r); 
  Ghat_F_0_vpar_r[2] = (0.5*(F_0_cr[2]*max_b_r+(F_0_rl[2]+F_0_cr[2])*avg_b_r)-0.5*F_0_rl[2]*max_b_r)*wvpar+dvpar*(0.1290994448735806*(F_0_cr[1]*max_b_r+(F_0_rl[1]+F_0_cr[1])*avg_b_r)-0.1290994448735806*F_0_rl[1]*max_b_r); 
  Ghat_G_1_vpar_r[0] = (0.5*(G_1_cr[0]*max_b_r+(G_1_rl[0]+G_1_cr[0])*avg_b_r)-0.5*G_1_rl[0]*max_b_r)*wvpar+dvpar*(0.1443375672974065*(G_1_cr[1]*max_b_r+(G_1_rl[1]+G_1_cr[1])*avg_b_r)-0.1443375672974065*G_1_rl[1]*max_b_r); 
  Ghat_G_1_vpar_r[1] = (0.5*(G_1_cr[1]*max_b_r+(G_1_rl[1]+G_1_cr[1])*avg_b_r)-0.5*G_1_rl[1]*max_b_r)*wvpar+dvpar*(((-0.1290994448735806*G_1_rl[2])+0.1290994448735806*G_1_cr[2]-0.1443375672974065*G_1_rl[0]+0.1443375672974065*G_1_cr[0])*max_b_r+(0.1290994448735806*(G_1_rl[2]+G_1_cr[2])+0.1443375672974065*(G_1_rl[0]+G_1_cr[0]))*avg_b_r); 
  Ghat_G_1_vpar_r[2] = (0.5*(G_1_cr[2]*max_b_r+(G_1_rl[2]+G_1_cr[2])*avg_b_r)-0.5*G_1_rl[2]*max_b_r)*wvpar+dvpar*(0.1290994448735806*(G_1_cr[1]*max_b_r+(G_1_rl[1]+G_1_cr[1])*avg_b_r)-0.1290994448735806*G_1_rl[1]*max_b_r); 

  } else { 

  Ghat_F_0_vpar_l[0] = (0.5*(F_0_cl[0]*max_b_l+(F_0_lr[0]+F_0_cl[0])*avg_b_l)-0.5*F_0_lr[0]*max_b_l)*wvpar+dvpar*(0.1443375672974065*(F_0_cl[1]*max_b_l+(F_0_lr[1]+F_0_cl[1])*avg_b_l)-0.1443375672974065*F_0_lr[1]*max_b_l); 
  Ghat_F_0_vpar_l[1] = (0.5*(F_0_cl[1]*max_b_l+(F_0_lr[1]+F_0_cl[1])*avg_b_l)-0.5*F_0_lr[1]*max_b_l)*wvpar+dvpar*(((-0.1290994448735806*F_0_lr[2])+0.1290994448735806*F_0_cl[2]-0.1443375672974065*F_0_lr[0]+0.1443375672974065*F_0_cl[0])*max_b_l+(0.1290994448735806*(F_0_lr[2]+F_0_cl[2])+0.1443375672974065*(F_0_lr[0]+F_0_cl[0]))*avg_b_l); 
  Ghat_F_0_vpar_l[2] = (0.5*(F_0_cl[2]*max_b_l+(F_0_lr[2]+F_0_cl[2])*avg_b_l)-0.5*F_0_lr[2]*max_b_l)*wvpar+dvpar*(0.1290994448735806*(F_0_cl[1]*max_b_l+(F_0_lr[1]+F_0_cl[1])*avg_b_l)-0.1290994448735806*F_0_lr[1]*max_b_l); 
  Ghat_G_1_vpar_l[0] = (0.5*(G_1_cl[0]*max_b_l+(G_1_lr[0]+G_1_cl[0])*avg_b_l)-0.5*G_1_lr[0]*max_b_l)*wvpar+dvpar*(0.1443375672974065*(G_1_cl[1]*max_b_l+(G_1_lr[1]+G_1_cl[1])*avg_b_l)-0.1443375672974065*G_1_lr[1]*max_b_l); 
  Ghat_G_1_vpar_l[1] = (0.5*(G_1_cl[1]*max_b_l+(G_1_lr[1]+G_1_cl[1])*avg_b_l)-0.5*G_1_lr[1]*max_b_l)*wvpar+dvpar*(((-0.1290994448735806*G_1_lr[2])+0.1290994448735806*G_1_cl[2]-0.1443375672974065*G_1_lr[0]+0.1443375672974065*G_1_cl[0])*max_b_l+(0.1290994448735806*(G_1_lr[2]+G_1_cl[2])+0.1443375672974065*(G_1_lr[0]+G_1_cl[0]))*avg_b_l); 
  Ghat_G_1_vpar_l[2] = (0.5*(G_1_cl[2]*max_b_l+(G_1_lr[2]+G_1_cl[2])*avg_b_l)-0.5*G_1_lr[2]*max_b_l)*wvpar+dvpar*(0.1290994448735806*(G_1_cl[1]*max_b_l+(G_1_lr[1]+G_1_cl[1])*avg_b_l)-0.1290994448735806*G_1_lr[1]*max_b_l); 

  Ghat_F_0_vpar_r[0] = ((0.5*F_0_rl[0]-0.5*F_0_cr[0])*max_b_r+0.5*(F_0_rl[0]+F_0_cr[0])*avg_b_r)*wvpar+dvpar*((0.1443375672974065*F_0_rl[1]-0.1443375672974065*F_0_cr[1])*max_b_r+0.1443375672974065*(F_0_rl[1]+F_0_cr[1])*avg_b_r); 
  Ghat_F_0_vpar_r[1] = ((0.5*F_0_rl[1]-0.5*F_0_cr[1])*max_b_r+0.5*(F_0_rl[1]+F_0_cr[1])*avg_b_r)*wvpar+dvpar*((0.1290994448735806*F_0_rl[2]-0.1290994448735806*F_0_cr[2]+0.1443375672974065*F_0_rl[0]-0.1443375672974065*F_0_cr[0])*max_b_r+(0.1290994448735806*(F_0_rl[2]+F_0_cr[2])+0.1443375672974065*(F_0_rl[0]+F_0_cr[0]))*avg_b_r); 
  Ghat_F_0_vpar_r[2] = ((0.5*F_0_rl[2]-0.5*F_0_cr[2])*max_b_r+0.5*(F_0_rl[2]+F_0_cr[2])*avg_b_r)*wvpar+dvpar*((0.1290994448735806*F_0_rl[1]-0.1290994448735806*F_0_cr[1])*max_b_r+0.1290994448735806*(F_0_rl[1]+F_0_cr[1])*avg_b_r); 
  Ghat_G_1_vpar_r[0] = ((0.5*G_1_rl[0]-0.5*G_1_cr[0])*max_b_r+0.5*(G_1_rl[0]+G_1_cr[0])*avg_b_r)*wvpar+dvpar*((0.1443375672974065*G_1_rl[1]-0.1443375672974065*G_1_cr[1])*max_b_r+0.1443375672974065*(G_1_rl[1]+G_1_cr[1])*avg_b_r); 
  Ghat_G_1_vpar_r[1] = ((0.5*G_1_rl[1]-0.5*G_1_cr[1])*max_b_r+0.5*(G_1_rl[1]+G_1_cr[1])*avg_b_r)*wvpar+dvpar*((0.1290994448735806*G_1_rl[2]-0.1290994448735806*G_1_cr[2]+0.1443375672974065*G_1_rl[0]-0.1443375672974065*G_1_cr[0])*max_b_r+(0.1290994448735806*(G_1_rl[2]+G_1_cr[2])+0.1443375672974065*(G_1_rl[0]+G_1_cr[0]))*avg_b_r); 
  Ghat_G_1_vpar_r[2] = ((0.5*G_1_rl[2]-0.5*G_1_cr[2])*max_b_r+0.5*(G_1_rl[2]+G_1_cr[2])*avg_b_r)*wvpar+dvpar*((0.1290994448735806*G_1_rl[1]-0.1290994448735806*G_1_cr[1])*max_b_r+0.1290994448735806*(G_1_rl[1]+G_1_cr[1])*avg_b_r); 

  } 
  double ul_r = u_surf_lr[0]; 
  double uc_l = u_surf_cl[0]; 
  double uc_r = u_surf_cr[0]; 
  double ur_l = u_surf_rl[0]; 
  double avg_u_l = 0.5*(ul_r + uc_l); 
  double avg_u_r = 0.5*(uc_r + ur_l); 

  double max_speed_l = pkpm_lax_l[0]; 
  double max_speed_r = pkpm_lax_r[0]; 

  Ghat_F_0_u_l[0] = (0.5*F_0_lr[0]-0.5*F_0_cl[0])*max_speed_l+0.5*(F_0_lr[0]+F_0_cl[0])*avg_u_l; 
  Ghat_F_0_u_l[1] = (0.5*F_0_lr[1]-0.5*F_0_cl[1])*max_speed_l+0.5*(F_0_lr[1]+F_0_cl[1])*avg_u_l; 
  Ghat_F_0_u_l[2] = (0.5*F_0_lr[2]-0.5*F_0_cl[2])*max_speed_l+0.5*(F_0_lr[2]+F_0_cl[2])*avg_u_l; 
  Ghat_G_1_u_l[0] = (0.5*G_1_lr[0]-0.5*G_1_cl[0])*max_speed_l+0.5*(G_1_lr[0]+G_1_cl[0])*avg_u_l; 
  Ghat_G_1_u_l[1] = (0.5*G_1_lr[1]-0.5*G_1_cl[1])*max_speed_l+0.5*(G_1_lr[1]+G_1_cl[1])*avg_u_l; 
  Ghat_G_1_u_l[2] = (0.5*G_1_lr[2]-0.5*G_1_cl[2])*max_speed_l+0.5*(G_1_lr[2]+G_1_cl[2])*avg_u_l; 

  Ghat_F_0_u_r[0] = 0.5*(F_0_cr[0]*max_speed_r+(F_0_rl[0]+F_0_cr[0])*avg_u_r)-0.5*F_0_rl[0]*max_speed_r; 
  Ghat_F_0_u_r[1] = 0.5*(F_0_cr[1]*max_speed_r+(F_0_rl[1]+F_0_cr[1])*avg_u_r)-0.5*F_0_rl[1]*max_speed_r; 
  Ghat_F_0_u_r[2] = 0.5*(F_0_cr[2]*max_speed_r+(F_0_rl[2]+F_0_cr[2])*avg_u_r)-0.5*F_0_rl[2]*max_speed_r; 
  Ghat_G_1_u_r[0] = 0.5*(G_1_cr[0]*max_speed_r+(G_1_rl[0]+G_1_cr[0])*avg_u_r)-0.5*G_1_rl[0]*max_speed_r; 
  Ghat_G_1_u_r[1] = 0.5*(G_1_cr[1]*max_speed_r+(G_1_rl[1]+G_1_cr[1])*avg_u_r)-0.5*G_1_rl[1]*max_speed_r; 
  Ghat_G_1_u_r[2] = 0.5*(G_1_cr[2]*max_speed_r+(G_1_rl[2]+G_1_cr[2])*avg_u_r)-0.5*G_1_rl[2]*max_speed_r; 

  out_F_0[0] += ((-0.7071067811865475*Ghat_F_0_vpar_r[0])+0.7071067811865475*Ghat_F_0_vpar_l[0]-0.7071067811865475*Ghat_F_0_u_r[0]+0.7071067811865475*Ghat_F_0_u_l[0])*dx1; 
  out_F_0[1] += -1.224744871391589*(Ghat_F_0_vpar_r[0]+Ghat_F_0_vpar_l[0]+Ghat_F_0_u_r[0]+Ghat_F_0_u_l[0])*dx1; 
  out_F_0[2] += ((-0.7071067811865475*Ghat_F_0_vpar_r[1])+0.7071067811865475*Ghat_F_0_vpar_l[1]-0.7071067811865475*Ghat_F_0_u_r[1]+0.7071067811865475*Ghat_F_0_u_l[1])*dx1; 
  out_F_0[3] += -1.224744871391589*(Ghat_F_0_vpar_r[1]+Ghat_F_0_vpar_l[1]+Ghat_F_0_u_r[1]+Ghat_F_0_u_l[1])*dx1; 
  out_F_0[4] += ((-0.7071067811865475*Ghat_F_0_vpar_r[2])+0.7071067811865475*Ghat_F_0_vpar_l[2]-0.7071067811865475*Ghat_F_0_u_r[2]+0.7071067811865475*Ghat_F_0_u_l[2])*dx1; 
  out_F_0[5] += -1.224744871391589*(Ghat_F_0_vpar_r[2]+Ghat_F_0_vpar_l[2]+Ghat_F_0_u_r[2]+Ghat_F_0_u_l[2])*dx1; 
  out_G_1[0] += ((-0.7071067811865475*Ghat_G_1_vpar_r[0])+0.7071067811865475*Ghat_G_1_vpar_l[0]-0.7071067811865475*Ghat_G_1_u_r[0]+0.7071067811865475*Ghat_G_1_u_l[0])*dx1; 
  out_G_1[1] += -1.224744871391589*(Ghat_G_1_vpar_r[0]+Ghat_G_1_vpar_l[0]+Ghat_G_1_u_r[0]+Ghat_G_1_u_l[0])*dx1; 
  out_G_1[2] += ((-0.7071067811865475*Ghat_G_1_vpar_r[1])+0.7071067811865475*Ghat_G_1_vpar_l[1]-0.7071067811865475*Ghat_G_1_u_r[1]+0.7071067811865475*Ghat_G_1_u_l[1])*dx1; 
  out_G_1[3] += -1.224744871391589*(Ghat_G_1_vpar_r[1]+Ghat_G_1_vpar_l[1]+Ghat_G_1_u_r[1]+Ghat_G_1_u_l[1])*dx1; 
  out_G_1[4] += ((-0.7071067811865475*Ghat_G_1_vpar_r[2])+0.7071067811865475*Ghat_G_1_vpar_l[2]-0.7071067811865475*Ghat_G_1_u_r[2]+0.7071067811865475*Ghat_G_1_u_l[2])*dx1; 
  out_G_1[5] += -1.224744871391589*(Ghat_G_1_vpar_r[2]+Ghat_G_1_vpar_l[2]+Ghat_G_1_u_r[2]+Ghat_G_1_u_l[2])*dx1; 

  return 0.;

} 
