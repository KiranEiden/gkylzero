#include <gkyl_canonical_pb_kernels.h>
#include <gkyl_basis_ser_4x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double canonical_pb_surfy_2x2v_ser_p1(const double *w, const double *dxv, const double *hamil, 
  const double *alpha_surf_l, const double *alpha_surf_r, 
  const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
  const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // hamil: hamiltonian.
  // alpha_surf_l: Surface expansion of phase space flux on the left.
  // alpha_surf_r: Surface expansion of phase space flux on the right.
  // sgn_alpha_surf_l: sign(alpha_surf_l) at quadrature points.
  // sgn_alpha_surf_r: sign(alpha_surf_r) at quadrature points.
  // const_sgn_alpha_l: Boolean array true if sign(alpha_surf_l) is only one sign, either +1 or -1.
  // const_sgn_alpha_r: Boolean array true if sign(alpha_surf_r) is only one sign, either +1 or -1.
  // fl,fc,fr: distribution function in left, center and right cells.
  // out: output increment in center cell.

  double rdy2 = 2.0/dxv[1];

  const double *alphaL = &alpha_surf_l[8];
  const double *alphaR = &alpha_surf_r[8];
  const double *sgn_alpha_surfL = &sgn_alpha_surf_l[8];
  const double *sgn_alpha_surfR = &sgn_alpha_surf_r[8];
  const int *const_sgn_alphaL = &const_sgn_alpha_l[1];
  const int *const_sgn_alphaR = &const_sgn_alpha_r[1];

  double fUpL[8] = {0.};
  if (const_sgn_alphaL[0] == 1) {  
    if (sgn_alpha_surfL[0] == 1.0) {  
  fUpL[0] = 1.224744871391589*fl[2]+0.7071067811865475*fl[0]; 
  fUpL[1] = 1.224744871391589*fl[5]+0.7071067811865475*fl[1]; 
  fUpL[2] = 1.224744871391589*fl[7]+0.7071067811865475*fl[3]; 
  fUpL[3] = 1.224744871391589*fl[9]+0.7071067811865475*fl[4]; 
  fUpL[4] = 1.224744871391589*fl[11]+0.7071067811865475*fl[6]; 
  fUpL[5] = 1.224744871391589*fl[12]+0.7071067811865475*fl[8]; 
  fUpL[6] = 1.224744871391589*fl[14]+0.7071067811865475*fl[10]; 
  fUpL[7] = 1.224744871391589*fl[15]+0.7071067811865475*fl[13]; 
    } else { 
  fUpL[0] = 0.7071067811865475*fc[0]-1.224744871391589*fc[2]; 
  fUpL[1] = 0.7071067811865475*fc[1]-1.224744871391589*fc[5]; 
  fUpL[2] = 0.7071067811865475*fc[3]-1.224744871391589*fc[7]; 
  fUpL[3] = 0.7071067811865475*fc[4]-1.224744871391589*fc[9]; 
  fUpL[4] = 0.7071067811865475*fc[6]-1.224744871391589*fc[11]; 
  fUpL[5] = 0.7071067811865475*fc[8]-1.224744871391589*fc[12]; 
  fUpL[6] = 0.7071067811865475*fc[10]-1.224744871391589*fc[14]; 
  fUpL[7] = 0.7071067811865475*fc[13]-1.224744871391589*fc[15]; 
    } 
  } else { 
  double f_lr[8] = {0.};
  double f_cl[8] = {0.};
  double sgn_alphaUpL[8] = {0.};
  ser_4x_p1_upwind_quad_to_modal(sgn_alpha_surfL, sgn_alphaUpL); 

  f_lr[0] = 1.224744871391589*fl[2]+0.7071067811865475*fl[0]; 
  f_lr[1] = 1.224744871391589*fl[5]+0.7071067811865475*fl[1]; 
  f_lr[2] = 1.224744871391589*fl[7]+0.7071067811865475*fl[3]; 
  f_lr[3] = 1.224744871391589*fl[9]+0.7071067811865475*fl[4]; 
  f_lr[4] = 1.224744871391589*fl[11]+0.7071067811865475*fl[6]; 
  f_lr[5] = 1.224744871391589*fl[12]+0.7071067811865475*fl[8]; 
  f_lr[6] = 1.224744871391589*fl[14]+0.7071067811865475*fl[10]; 
  f_lr[7] = 1.224744871391589*fl[15]+0.7071067811865475*fl[13]; 

  f_cl[0] = 0.7071067811865475*fc[0]-1.224744871391589*fc[2]; 
  f_cl[1] = 0.7071067811865475*fc[1]-1.224744871391589*fc[5]; 
  f_cl[2] = 0.7071067811865475*fc[3]-1.224744871391589*fc[7]; 
  f_cl[3] = 0.7071067811865475*fc[4]-1.224744871391589*fc[9]; 
  f_cl[4] = 0.7071067811865475*fc[6]-1.224744871391589*fc[11]; 
  f_cl[5] = 0.7071067811865475*fc[8]-1.224744871391589*fc[12]; 
  f_cl[6] = 0.7071067811865475*fc[10]-1.224744871391589*fc[14]; 
  f_cl[7] = 0.7071067811865475*fc[13]-1.224744871391589*fc[15]; 

  fUpL[0] = (0.1767766952966368*f_lr[7]-0.1767766952966368*f_cl[7])*sgn_alphaUpL[7]+(0.1767766952966368*f_lr[6]-0.1767766952966368*f_cl[6])*sgn_alphaUpL[6]+(0.1767766952966368*f_lr[5]-0.1767766952966368*f_cl[5])*sgn_alphaUpL[5]+(0.1767766952966368*f_lr[4]-0.1767766952966368*f_cl[4])*sgn_alphaUpL[4]+(0.1767766952966368*f_lr[3]-0.1767766952966368*f_cl[3])*sgn_alphaUpL[3]+(0.1767766952966368*f_lr[2]-0.1767766952966368*f_cl[2])*sgn_alphaUpL[2]+(0.1767766952966368*f_lr[1]-0.1767766952966368*f_cl[1])*sgn_alphaUpL[1]+(0.1767766952966368*f_lr[0]-0.1767766952966368*f_cl[0])*sgn_alphaUpL[0]+0.5*(f_lr[0]+f_cl[0]); 
  fUpL[1] = (0.1767766952966368*f_lr[6]-0.1767766952966368*f_cl[6])*sgn_alphaUpL[7]+sgn_alphaUpL[6]*(0.1767766952966368*f_lr[7]-0.1767766952966368*f_cl[7])+(0.1767766952966368*f_lr[3]-0.1767766952966368*f_cl[3])*sgn_alphaUpL[5]+sgn_alphaUpL[3]*(0.1767766952966368*f_lr[5]-0.1767766952966368*f_cl[5])+(0.1767766952966368*f_lr[2]-0.1767766952966368*f_cl[2])*sgn_alphaUpL[4]+sgn_alphaUpL[2]*(0.1767766952966368*f_lr[4]-0.1767766952966368*f_cl[4])+(0.1767766952966368*f_lr[0]-0.1767766952966368*f_cl[0])*sgn_alphaUpL[1]+(0.1767766952966368*sgn_alphaUpL[0]+0.5)*f_lr[1]+(0.5-0.1767766952966368*sgn_alphaUpL[0])*f_cl[1]; 
  fUpL[2] = (0.1767766952966368*f_lr[5]-0.1767766952966368*f_cl[5])*sgn_alphaUpL[7]+sgn_alphaUpL[5]*(0.1767766952966368*f_lr[7]-0.1767766952966368*f_cl[7])+(0.1767766952966368*f_lr[3]-0.1767766952966368*f_cl[3])*sgn_alphaUpL[6]+sgn_alphaUpL[3]*(0.1767766952966368*f_lr[6]-0.1767766952966368*f_cl[6])+(0.1767766952966368*f_lr[1]-0.1767766952966368*f_cl[1])*sgn_alphaUpL[4]+sgn_alphaUpL[1]*(0.1767766952966368*f_lr[4]-0.1767766952966368*f_cl[4])+(0.1767766952966368*f_lr[0]-0.1767766952966368*f_cl[0])*sgn_alphaUpL[2]+(0.1767766952966368*sgn_alphaUpL[0]+0.5)*f_lr[2]+(0.5-0.1767766952966368*sgn_alphaUpL[0])*f_cl[2]; 
  fUpL[3] = (0.1767766952966368*f_lr[4]-0.1767766952966368*f_cl[4])*sgn_alphaUpL[7]+sgn_alphaUpL[4]*(0.1767766952966368*f_lr[7]-0.1767766952966368*f_cl[7])+(0.1767766952966368*f_lr[2]-0.1767766952966368*f_cl[2])*sgn_alphaUpL[6]+sgn_alphaUpL[2]*(0.1767766952966368*f_lr[6]-0.1767766952966368*f_cl[6])+(0.1767766952966368*f_lr[1]-0.1767766952966368*f_cl[1])*sgn_alphaUpL[5]+sgn_alphaUpL[1]*(0.1767766952966368*f_lr[5]-0.1767766952966368*f_cl[5])+(0.1767766952966368*f_lr[0]-0.1767766952966368*f_cl[0])*sgn_alphaUpL[3]+(0.1767766952966368*sgn_alphaUpL[0]+0.5)*f_lr[3]+(0.5-0.1767766952966368*sgn_alphaUpL[0])*f_cl[3]; 
  fUpL[4] = (0.1767766952966368*f_lr[3]-0.1767766952966368*f_cl[3])*sgn_alphaUpL[7]+sgn_alphaUpL[3]*(0.1767766952966368*f_lr[7]-0.1767766952966368*f_cl[7])+(0.1767766952966368*f_lr[5]-0.1767766952966368*f_cl[5])*sgn_alphaUpL[6]+sgn_alphaUpL[5]*(0.1767766952966368*f_lr[6]-0.1767766952966368*f_cl[6])+(0.1767766952966368*f_lr[0]-0.1767766952966368*f_cl[0])*sgn_alphaUpL[4]+(0.1767766952966368*sgn_alphaUpL[0]+0.5)*f_lr[4]+(0.5-0.1767766952966368*sgn_alphaUpL[0])*f_cl[4]+(0.1767766952966368*f_lr[1]-0.1767766952966368*f_cl[1])*sgn_alphaUpL[2]+sgn_alphaUpL[1]*(0.1767766952966368*f_lr[2]-0.1767766952966368*f_cl[2]); 
  fUpL[5] = (0.1767766952966368*f_lr[2]-0.1767766952966368*f_cl[2])*sgn_alphaUpL[7]+sgn_alphaUpL[2]*(0.1767766952966368*f_lr[7]-0.1767766952966368*f_cl[7])+(0.1767766952966368*f_lr[4]-0.1767766952966368*f_cl[4])*sgn_alphaUpL[6]+sgn_alphaUpL[4]*(0.1767766952966368*f_lr[6]-0.1767766952966368*f_cl[6])+(0.1767766952966368*f_lr[0]-0.1767766952966368*f_cl[0])*sgn_alphaUpL[5]+(0.1767766952966368*sgn_alphaUpL[0]+0.5)*f_lr[5]+(0.5-0.1767766952966368*sgn_alphaUpL[0])*f_cl[5]+(0.1767766952966368*f_lr[1]-0.1767766952966368*f_cl[1])*sgn_alphaUpL[3]+sgn_alphaUpL[1]*(0.1767766952966368*f_lr[3]-0.1767766952966368*f_cl[3]); 
  fUpL[6] = (0.1767766952966368*f_lr[1]-0.1767766952966368*f_cl[1])*sgn_alphaUpL[7]+sgn_alphaUpL[1]*(0.1767766952966368*f_lr[7]-0.1767766952966368*f_cl[7])+(0.1767766952966368*f_lr[0]-0.1767766952966368*f_cl[0])*sgn_alphaUpL[6]+(0.1767766952966368*sgn_alphaUpL[0]+0.5)*f_lr[6]+(0.5-0.1767766952966368*sgn_alphaUpL[0])*f_cl[6]+(0.1767766952966368*f_lr[4]-0.1767766952966368*f_cl[4])*sgn_alphaUpL[5]+sgn_alphaUpL[4]*(0.1767766952966368*f_lr[5]-0.1767766952966368*f_cl[5])+(0.1767766952966368*f_lr[2]-0.1767766952966368*f_cl[2])*sgn_alphaUpL[3]+sgn_alphaUpL[2]*(0.1767766952966368*f_lr[3]-0.1767766952966368*f_cl[3]); 
  fUpL[7] = (0.1767766952966368*f_lr[0]-0.1767766952966368*f_cl[0])*sgn_alphaUpL[7]+(0.1767766952966368*sgn_alphaUpL[0]+0.5)*f_lr[7]+(0.5-0.1767766952966368*sgn_alphaUpL[0])*f_cl[7]+(0.1767766952966368*f_lr[1]-0.1767766952966368*f_cl[1])*sgn_alphaUpL[6]+sgn_alphaUpL[1]*(0.1767766952966368*f_lr[6]-0.1767766952966368*f_cl[6])+(0.1767766952966368*f_lr[2]-0.1767766952966368*f_cl[2])*sgn_alphaUpL[5]+sgn_alphaUpL[2]*(0.1767766952966368*f_lr[5]-0.1767766952966368*f_cl[5])+(0.1767766952966368*f_lr[3]-0.1767766952966368*f_cl[3])*sgn_alphaUpL[4]+sgn_alphaUpL[3]*(0.1767766952966368*f_lr[4]-0.1767766952966368*f_cl[4]); 

  } 
  double fUpR[8] = {0.};
  if (const_sgn_alphaR[0] == 1) {  
    if (sgn_alpha_surfR[0] == 1.0) {  
  fUpR[0] = 1.224744871391589*fc[2]+0.7071067811865475*fc[0]; 
  fUpR[1] = 1.224744871391589*fc[5]+0.7071067811865475*fc[1]; 
  fUpR[2] = 1.224744871391589*fc[7]+0.7071067811865475*fc[3]; 
  fUpR[3] = 1.224744871391589*fc[9]+0.7071067811865475*fc[4]; 
  fUpR[4] = 1.224744871391589*fc[11]+0.7071067811865475*fc[6]; 
  fUpR[5] = 1.224744871391589*fc[12]+0.7071067811865475*fc[8]; 
  fUpR[6] = 1.224744871391589*fc[14]+0.7071067811865475*fc[10]; 
  fUpR[7] = 1.224744871391589*fc[15]+0.7071067811865475*fc[13]; 
    } else { 
  fUpR[0] = 0.7071067811865475*fr[0]-1.224744871391589*fr[2]; 
  fUpR[1] = 0.7071067811865475*fr[1]-1.224744871391589*fr[5]; 
  fUpR[2] = 0.7071067811865475*fr[3]-1.224744871391589*fr[7]; 
  fUpR[3] = 0.7071067811865475*fr[4]-1.224744871391589*fr[9]; 
  fUpR[4] = 0.7071067811865475*fr[6]-1.224744871391589*fr[11]; 
  fUpR[5] = 0.7071067811865475*fr[8]-1.224744871391589*fr[12]; 
  fUpR[6] = 0.7071067811865475*fr[10]-1.224744871391589*fr[14]; 
  fUpR[7] = 0.7071067811865475*fr[13]-1.224744871391589*fr[15]; 
    } 
  } else { 
  double f_cr[8] = {0.};
  double f_rl[8] = {0.};
  double sgn_alphaUpR[8] = {0.};
  ser_4x_p1_upwind_quad_to_modal(sgn_alpha_surfR, sgn_alphaUpR); 

  f_cr[0] = 1.224744871391589*fc[2]+0.7071067811865475*fc[0]; 
  f_cr[1] = 1.224744871391589*fc[5]+0.7071067811865475*fc[1]; 
  f_cr[2] = 1.224744871391589*fc[7]+0.7071067811865475*fc[3]; 
  f_cr[3] = 1.224744871391589*fc[9]+0.7071067811865475*fc[4]; 
  f_cr[4] = 1.224744871391589*fc[11]+0.7071067811865475*fc[6]; 
  f_cr[5] = 1.224744871391589*fc[12]+0.7071067811865475*fc[8]; 
  f_cr[6] = 1.224744871391589*fc[14]+0.7071067811865475*fc[10]; 
  f_cr[7] = 1.224744871391589*fc[15]+0.7071067811865475*fc[13]; 

  f_rl[0] = 0.7071067811865475*fr[0]-1.224744871391589*fr[2]; 
  f_rl[1] = 0.7071067811865475*fr[1]-1.224744871391589*fr[5]; 
  f_rl[2] = 0.7071067811865475*fr[3]-1.224744871391589*fr[7]; 
  f_rl[3] = 0.7071067811865475*fr[4]-1.224744871391589*fr[9]; 
  f_rl[4] = 0.7071067811865475*fr[6]-1.224744871391589*fr[11]; 
  f_rl[5] = 0.7071067811865475*fr[8]-1.224744871391589*fr[12]; 
  f_rl[6] = 0.7071067811865475*fr[10]-1.224744871391589*fr[14]; 
  f_rl[7] = 0.7071067811865475*fr[13]-1.224744871391589*fr[15]; 

  fUpR[0] = (0.1767766952966368*f_cr[7]-0.1767766952966368*f_rl[7])*sgn_alphaUpR[7]+(0.1767766952966368*f_cr[6]-0.1767766952966368*f_rl[6])*sgn_alphaUpR[6]+(0.1767766952966368*f_cr[5]-0.1767766952966368*f_rl[5])*sgn_alphaUpR[5]+(0.1767766952966368*f_cr[4]-0.1767766952966368*f_rl[4])*sgn_alphaUpR[4]+(0.1767766952966368*f_cr[3]-0.1767766952966368*f_rl[3])*sgn_alphaUpR[3]+(0.1767766952966368*f_cr[2]-0.1767766952966368*f_rl[2])*sgn_alphaUpR[2]+(0.1767766952966368*f_cr[1]-0.1767766952966368*f_rl[1])*sgn_alphaUpR[1]+(0.1767766952966368*f_cr[0]-0.1767766952966368*f_rl[0])*sgn_alphaUpR[0]+0.5*(f_rl[0]+f_cr[0]); 
  fUpR[1] = (0.1767766952966368*f_cr[6]-0.1767766952966368*f_rl[6])*sgn_alphaUpR[7]+sgn_alphaUpR[6]*(0.1767766952966368*f_cr[7]-0.1767766952966368*f_rl[7])+(0.1767766952966368*f_cr[3]-0.1767766952966368*f_rl[3])*sgn_alphaUpR[5]+sgn_alphaUpR[3]*(0.1767766952966368*f_cr[5]-0.1767766952966368*f_rl[5])+(0.1767766952966368*f_cr[2]-0.1767766952966368*f_rl[2])*sgn_alphaUpR[4]+sgn_alphaUpR[2]*(0.1767766952966368*f_cr[4]-0.1767766952966368*f_rl[4])+(0.1767766952966368*f_cr[0]-0.1767766952966368*f_rl[0])*sgn_alphaUpR[1]+(0.5-0.1767766952966368*sgn_alphaUpR[0])*f_rl[1]+(0.1767766952966368*sgn_alphaUpR[0]+0.5)*f_cr[1]; 
  fUpR[2] = (0.1767766952966368*f_cr[5]-0.1767766952966368*f_rl[5])*sgn_alphaUpR[7]+sgn_alphaUpR[5]*(0.1767766952966368*f_cr[7]-0.1767766952966368*f_rl[7])+(0.1767766952966368*f_cr[3]-0.1767766952966368*f_rl[3])*sgn_alphaUpR[6]+sgn_alphaUpR[3]*(0.1767766952966368*f_cr[6]-0.1767766952966368*f_rl[6])+(0.1767766952966368*f_cr[1]-0.1767766952966368*f_rl[1])*sgn_alphaUpR[4]+sgn_alphaUpR[1]*(0.1767766952966368*f_cr[4]-0.1767766952966368*f_rl[4])+(0.1767766952966368*f_cr[0]-0.1767766952966368*f_rl[0])*sgn_alphaUpR[2]+(0.5-0.1767766952966368*sgn_alphaUpR[0])*f_rl[2]+(0.1767766952966368*sgn_alphaUpR[0]+0.5)*f_cr[2]; 
  fUpR[3] = (0.1767766952966368*f_cr[4]-0.1767766952966368*f_rl[4])*sgn_alphaUpR[7]+sgn_alphaUpR[4]*(0.1767766952966368*f_cr[7]-0.1767766952966368*f_rl[7])+(0.1767766952966368*f_cr[2]-0.1767766952966368*f_rl[2])*sgn_alphaUpR[6]+sgn_alphaUpR[2]*(0.1767766952966368*f_cr[6]-0.1767766952966368*f_rl[6])+(0.1767766952966368*f_cr[1]-0.1767766952966368*f_rl[1])*sgn_alphaUpR[5]+sgn_alphaUpR[1]*(0.1767766952966368*f_cr[5]-0.1767766952966368*f_rl[5])+(0.1767766952966368*f_cr[0]-0.1767766952966368*f_rl[0])*sgn_alphaUpR[3]+(0.5-0.1767766952966368*sgn_alphaUpR[0])*f_rl[3]+(0.1767766952966368*sgn_alphaUpR[0]+0.5)*f_cr[3]; 
  fUpR[4] = (0.1767766952966368*f_cr[3]-0.1767766952966368*f_rl[3])*sgn_alphaUpR[7]+sgn_alphaUpR[3]*(0.1767766952966368*f_cr[7]-0.1767766952966368*f_rl[7])+(0.1767766952966368*f_cr[5]-0.1767766952966368*f_rl[5])*sgn_alphaUpR[6]+sgn_alphaUpR[5]*(0.1767766952966368*f_cr[6]-0.1767766952966368*f_rl[6])+(0.1767766952966368*f_cr[0]-0.1767766952966368*f_rl[0])*sgn_alphaUpR[4]+(0.5-0.1767766952966368*sgn_alphaUpR[0])*f_rl[4]+(0.1767766952966368*sgn_alphaUpR[0]+0.5)*f_cr[4]+(0.1767766952966368*f_cr[1]-0.1767766952966368*f_rl[1])*sgn_alphaUpR[2]+sgn_alphaUpR[1]*(0.1767766952966368*f_cr[2]-0.1767766952966368*f_rl[2]); 
  fUpR[5] = (0.1767766952966368*f_cr[2]-0.1767766952966368*f_rl[2])*sgn_alphaUpR[7]+sgn_alphaUpR[2]*(0.1767766952966368*f_cr[7]-0.1767766952966368*f_rl[7])+(0.1767766952966368*f_cr[4]-0.1767766952966368*f_rl[4])*sgn_alphaUpR[6]+sgn_alphaUpR[4]*(0.1767766952966368*f_cr[6]-0.1767766952966368*f_rl[6])+(0.1767766952966368*f_cr[0]-0.1767766952966368*f_rl[0])*sgn_alphaUpR[5]+(0.5-0.1767766952966368*sgn_alphaUpR[0])*f_rl[5]+(0.1767766952966368*sgn_alphaUpR[0]+0.5)*f_cr[5]+(0.1767766952966368*f_cr[1]-0.1767766952966368*f_rl[1])*sgn_alphaUpR[3]+sgn_alphaUpR[1]*(0.1767766952966368*f_cr[3]-0.1767766952966368*f_rl[3]); 
  fUpR[6] = (0.1767766952966368*f_cr[1]-0.1767766952966368*f_rl[1])*sgn_alphaUpR[7]+sgn_alphaUpR[1]*(0.1767766952966368*f_cr[7]-0.1767766952966368*f_rl[7])+(0.1767766952966368*f_cr[0]-0.1767766952966368*f_rl[0])*sgn_alphaUpR[6]+(0.5-0.1767766952966368*sgn_alphaUpR[0])*f_rl[6]+(0.1767766952966368*sgn_alphaUpR[0]+0.5)*f_cr[6]+(0.1767766952966368*f_cr[4]-0.1767766952966368*f_rl[4])*sgn_alphaUpR[5]+sgn_alphaUpR[4]*(0.1767766952966368*f_cr[5]-0.1767766952966368*f_rl[5])+(0.1767766952966368*f_cr[2]-0.1767766952966368*f_rl[2])*sgn_alphaUpR[3]+sgn_alphaUpR[2]*(0.1767766952966368*f_cr[3]-0.1767766952966368*f_rl[3]); 
  fUpR[7] = (0.1767766952966368*f_cr[0]-0.1767766952966368*f_rl[0])*sgn_alphaUpR[7]+(0.5-0.1767766952966368*sgn_alphaUpR[0])*f_rl[7]+(0.1767766952966368*sgn_alphaUpR[0]+0.5)*f_cr[7]+(0.1767766952966368*f_cr[1]-0.1767766952966368*f_rl[1])*sgn_alphaUpR[6]+sgn_alphaUpR[1]*(0.1767766952966368*f_cr[6]-0.1767766952966368*f_rl[6])+(0.1767766952966368*f_cr[2]-0.1767766952966368*f_rl[2])*sgn_alphaUpR[5]+sgn_alphaUpR[2]*(0.1767766952966368*f_cr[5]-0.1767766952966368*f_rl[5])+(0.1767766952966368*f_cr[3]-0.1767766952966368*f_rl[3])*sgn_alphaUpR[4]+sgn_alphaUpR[3]*(0.1767766952966368*f_cr[4]-0.1767766952966368*f_rl[4]); 

  } 
  double GhatL[8] = {0.};
  double GhatR[8] = {0.};
  GhatL[0] = 0.3535533905932737*(alphaL[7]*fUpL[7]+alphaL[6]*fUpL[6]+alphaL[5]*fUpL[5]+alphaL[4]*fUpL[4]+alphaL[3]*fUpL[3]+alphaL[2]*fUpL[2]+alphaL[1]*fUpL[1]+alphaL[0]*fUpL[0]); 
  GhatL[1] = 0.3535533905932737*(alphaL[6]*fUpL[7]+fUpL[6]*alphaL[7]+alphaL[3]*fUpL[5]+fUpL[3]*alphaL[5]+alphaL[2]*fUpL[4]+fUpL[2]*alphaL[4]+alphaL[0]*fUpL[1]+fUpL[0]*alphaL[1]); 
  GhatL[2] = 0.3535533905932737*(alphaL[5]*fUpL[7]+fUpL[5]*alphaL[7]+alphaL[3]*fUpL[6]+fUpL[3]*alphaL[6]+alphaL[1]*fUpL[4]+fUpL[1]*alphaL[4]+alphaL[0]*fUpL[2]+fUpL[0]*alphaL[2]); 
  GhatL[3] = 0.3535533905932737*(alphaL[4]*fUpL[7]+fUpL[4]*alphaL[7]+alphaL[2]*fUpL[6]+fUpL[2]*alphaL[6]+alphaL[1]*fUpL[5]+fUpL[1]*alphaL[5]+alphaL[0]*fUpL[3]+fUpL[0]*alphaL[3]); 
  GhatL[4] = 0.3535533905932737*(alphaL[3]*fUpL[7]+fUpL[3]*alphaL[7]+alphaL[5]*fUpL[6]+fUpL[5]*alphaL[6]+alphaL[0]*fUpL[4]+fUpL[0]*alphaL[4]+alphaL[1]*fUpL[2]+fUpL[1]*alphaL[2]); 
  GhatL[5] = 0.3535533905932737*(alphaL[2]*fUpL[7]+fUpL[2]*alphaL[7]+alphaL[4]*fUpL[6]+fUpL[4]*alphaL[6]+alphaL[0]*fUpL[5]+fUpL[0]*alphaL[5]+alphaL[1]*fUpL[3]+fUpL[1]*alphaL[3]); 
  GhatL[6] = 0.3535533905932737*(alphaL[1]*fUpL[7]+fUpL[1]*alphaL[7]+alphaL[0]*fUpL[6]+fUpL[0]*alphaL[6]+alphaL[4]*fUpL[5]+fUpL[4]*alphaL[5]+alphaL[2]*fUpL[3]+fUpL[2]*alphaL[3]); 
  GhatL[7] = 0.3535533905932737*(alphaL[0]*fUpL[7]+fUpL[0]*alphaL[7]+alphaL[1]*fUpL[6]+fUpL[1]*alphaL[6]+alphaL[2]*fUpL[5]+fUpL[2]*alphaL[5]+alphaL[3]*fUpL[4]+fUpL[3]*alphaL[4]); 

  GhatR[0] = 0.3535533905932737*(alphaR[7]*fUpR[7]+alphaR[6]*fUpR[6]+alphaR[5]*fUpR[5]+alphaR[4]*fUpR[4]+alphaR[3]*fUpR[3]+alphaR[2]*fUpR[2]+alphaR[1]*fUpR[1]+alphaR[0]*fUpR[0]); 
  GhatR[1] = 0.3535533905932737*(alphaR[6]*fUpR[7]+fUpR[6]*alphaR[7]+alphaR[3]*fUpR[5]+fUpR[3]*alphaR[5]+alphaR[2]*fUpR[4]+fUpR[2]*alphaR[4]+alphaR[0]*fUpR[1]+fUpR[0]*alphaR[1]); 
  GhatR[2] = 0.3535533905932737*(alphaR[5]*fUpR[7]+fUpR[5]*alphaR[7]+alphaR[3]*fUpR[6]+fUpR[3]*alphaR[6]+alphaR[1]*fUpR[4]+fUpR[1]*alphaR[4]+alphaR[0]*fUpR[2]+fUpR[0]*alphaR[2]); 
  GhatR[3] = 0.3535533905932737*(alphaR[4]*fUpR[7]+fUpR[4]*alphaR[7]+alphaR[2]*fUpR[6]+fUpR[2]*alphaR[6]+alphaR[1]*fUpR[5]+fUpR[1]*alphaR[5]+alphaR[0]*fUpR[3]+fUpR[0]*alphaR[3]); 
  GhatR[4] = 0.3535533905932737*(alphaR[3]*fUpR[7]+fUpR[3]*alphaR[7]+alphaR[5]*fUpR[6]+fUpR[5]*alphaR[6]+alphaR[0]*fUpR[4]+fUpR[0]*alphaR[4]+alphaR[1]*fUpR[2]+fUpR[1]*alphaR[2]); 
  GhatR[5] = 0.3535533905932737*(alphaR[2]*fUpR[7]+fUpR[2]*alphaR[7]+alphaR[4]*fUpR[6]+fUpR[4]*alphaR[6]+alphaR[0]*fUpR[5]+fUpR[0]*alphaR[5]+alphaR[1]*fUpR[3]+fUpR[1]*alphaR[3]); 
  GhatR[6] = 0.3535533905932737*(alphaR[1]*fUpR[7]+fUpR[1]*alphaR[7]+alphaR[0]*fUpR[6]+fUpR[0]*alphaR[6]+alphaR[4]*fUpR[5]+fUpR[4]*alphaR[5]+alphaR[2]*fUpR[3]+fUpR[2]*alphaR[3]); 
  GhatR[7] = 0.3535533905932737*(alphaR[0]*fUpR[7]+fUpR[0]*alphaR[7]+alphaR[1]*fUpR[6]+fUpR[1]*alphaR[6]+alphaR[2]*fUpR[5]+fUpR[2]*alphaR[5]+alphaR[3]*fUpR[4]+fUpR[3]*alphaR[4]); 

  out[0] += (0.7071067811865475*GhatL[0]-0.7071067811865475*GhatR[0])*rdy2; 
  out[1] += (0.7071067811865475*GhatL[1]-0.7071067811865475*GhatR[1])*rdy2; 
  out[2] += ((-1.224744871391589*GhatR[0])-1.224744871391589*GhatL[0])*rdy2; 
  out[3] += (0.7071067811865475*GhatL[2]-0.7071067811865475*GhatR[2])*rdy2; 
  out[4] += (0.7071067811865475*GhatL[3]-0.7071067811865475*GhatR[3])*rdy2; 
  out[5] += ((-1.224744871391589*GhatR[1])-1.224744871391589*GhatL[1])*rdy2; 
  out[6] += (0.7071067811865475*GhatL[4]-0.7071067811865475*GhatR[4])*rdy2; 
  out[7] += ((-1.224744871391589*GhatR[2])-1.224744871391589*GhatL[2])*rdy2; 
  out[8] += (0.7071067811865475*GhatL[5]-0.7071067811865475*GhatR[5])*rdy2; 
  out[9] += ((-1.224744871391589*GhatR[3])-1.224744871391589*GhatL[3])*rdy2; 
  out[10] += (0.7071067811865475*GhatL[6]-0.7071067811865475*GhatR[6])*rdy2; 
  out[11] += ((-1.224744871391589*GhatR[4])-1.224744871391589*GhatL[4])*rdy2; 
  out[12] += ((-1.224744871391589*GhatR[5])-1.224744871391589*GhatL[5])*rdy2; 
  out[13] += (0.7071067811865475*GhatL[7]-0.7071067811865475*GhatR[7])*rdy2; 
  out[14] += ((-1.224744871391589*GhatR[6])-1.224744871391589*GhatL[6])*rdy2; 
  out[15] += ((-1.224744871391589*GhatR[7])-1.224744871391589*GhatL[7])*rdy2; 

  double cflFreq = fmax(fabs(alphaL[0]), fabs(alphaR[0])); 
  return 0.5303300858899105*rdy2*cflFreq; 

} 
