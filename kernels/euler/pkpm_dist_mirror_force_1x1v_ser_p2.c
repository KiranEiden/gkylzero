#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void pkpm_dist_mirror_force_1x1v_ser_p2(const double* T_perp_over_m, const double* T_perp_over_m_inv, const double* f, const double* F_k_p_1, double* g_dist_source, double* F_k_m_1) 
{ 
  // T_perp_over_m: Input p_perp/rho = T_perp/m.
  // T_perp_over_m_inv: Input (T_perp/m)^-1.
  // f:             Input distribution function [F_0, T_perp/m G = T_perp/m (F_0 - F_1)].
  // F_k_p_1:       Input k+1 distribution function. F_2 expansion is the first NP coefficients. 
  // g_dist_source: Output 2.0*T_perp/m*(2.0*T_perp/m G + T_perp/m (F_2 - F_0)).
  // F_k_m_1:       Output k-1 distribution function. F_1 expansion is the first NP coefficients. 

  const double *F_0 = &f[0]; 
  const double *G_1 = &f[8]; 
  const double *F_2 = &F_k_p_1[0]; 

  double *out_F_1 = &F_k_m_1[0]; 

  double tmp_F_2_m_F_0[8] = {0.0}; 
  double tmp_T_perp_F_0[8] = {0.0}; 
  double tmp_T_perp_g_dist[8] = {0.0}; 
  double tmp_F_0_m_F_1[8] = {0.0}; 

  tmp_F_2_m_F_0[0] = F_2[0]-1.0*F_0[0]; 
  tmp_F_2_m_F_0[1] = F_2[1]-1.0*F_0[1]; 
  tmp_F_2_m_F_0[2] = F_2[2]-1.0*F_0[2]; 
  tmp_F_2_m_F_0[3] = F_2[3]-1.0*F_0[3]; 
  tmp_F_2_m_F_0[4] = F_2[4]-1.0*F_0[4]; 
  tmp_F_2_m_F_0[5] = F_2[5]-1.0*F_0[5]; 
  tmp_F_2_m_F_0[6] = F_2[6]-1.0*F_0[6]; 
  tmp_F_2_m_F_0[7] = F_2[7]-1.0*F_0[7]; 

  tmp_T_perp_F_0[0] = 0.7071067811865475*T_perp_over_m[2]*tmp_F_2_m_F_0[4]+0.7071067811865475*T_perp_over_m[1]*tmp_F_2_m_F_0[1]+0.7071067811865475*T_perp_over_m[0]*tmp_F_2_m_F_0[0]; 
  tmp_T_perp_F_0[1] = 0.6324555320336759*T_perp_over_m[1]*tmp_F_2_m_F_0[4]+0.6324555320336759*tmp_F_2_m_F_0[1]*T_perp_over_m[2]+0.7071067811865475*T_perp_over_m[0]*tmp_F_2_m_F_0[1]+0.7071067811865475*tmp_F_2_m_F_0[0]*T_perp_over_m[1]; 
  tmp_T_perp_F_0[2] = 0.7071067811865475*T_perp_over_m[2]*tmp_F_2_m_F_0[6]+0.7071067811865475*T_perp_over_m[1]*tmp_F_2_m_F_0[3]+0.7071067811865475*T_perp_over_m[0]*tmp_F_2_m_F_0[2]; 
  tmp_T_perp_F_0[3] = 0.632455532033676*T_perp_over_m[1]*tmp_F_2_m_F_0[6]+0.6324555320336759*T_perp_over_m[2]*tmp_F_2_m_F_0[3]+0.7071067811865475*T_perp_over_m[0]*tmp_F_2_m_F_0[3]+0.7071067811865475*T_perp_over_m[1]*tmp_F_2_m_F_0[2]; 
  tmp_T_perp_F_0[4] = 0.4517539514526256*T_perp_over_m[2]*tmp_F_2_m_F_0[4]+0.7071067811865475*T_perp_over_m[0]*tmp_F_2_m_F_0[4]+0.7071067811865475*tmp_F_2_m_F_0[0]*T_perp_over_m[2]+0.6324555320336759*T_perp_over_m[1]*tmp_F_2_m_F_0[1]; 
  tmp_T_perp_F_0[5] = 0.7071067811865475*T_perp_over_m[1]*tmp_F_2_m_F_0[7]+0.7071067811865475*T_perp_over_m[0]*tmp_F_2_m_F_0[5]; 
  tmp_T_perp_F_0[6] = 0.4517539514526256*T_perp_over_m[2]*tmp_F_2_m_F_0[6]+0.7071067811865475*T_perp_over_m[0]*tmp_F_2_m_F_0[6]+0.632455532033676*T_perp_over_m[1]*tmp_F_2_m_F_0[3]+0.7071067811865475*T_perp_over_m[2]*tmp_F_2_m_F_0[2]; 
  tmp_T_perp_F_0[7] = 0.6324555320336759*T_perp_over_m[2]*tmp_F_2_m_F_0[7]+0.7071067811865475*T_perp_over_m[0]*tmp_F_2_m_F_0[7]+0.7071067811865475*T_perp_over_m[1]*tmp_F_2_m_F_0[5]; 

  tmp_T_perp_g_dist[0] = tmp_T_perp_F_0[0]+2.0*G_1[0]; 
  tmp_T_perp_g_dist[1] = tmp_T_perp_F_0[1]+2.0*G_1[1]; 
  tmp_T_perp_g_dist[2] = tmp_T_perp_F_0[2]+2.0*G_1[2]; 
  tmp_T_perp_g_dist[3] = tmp_T_perp_F_0[3]+2.0*G_1[3]; 
  tmp_T_perp_g_dist[4] = tmp_T_perp_F_0[4]+2.0*G_1[4]; 
  tmp_T_perp_g_dist[5] = tmp_T_perp_F_0[5]+2.0*G_1[5]; 
  tmp_T_perp_g_dist[6] = tmp_T_perp_F_0[6]+2.0*G_1[6]; 
  tmp_T_perp_g_dist[7] = tmp_T_perp_F_0[7]+2.0*G_1[7]; 

  g_dist_source[0] = 1.414213562373095*T_perp_over_m[2]*tmp_T_perp_g_dist[4]+1.414213562373095*T_perp_over_m[1]*tmp_T_perp_g_dist[1]+1.414213562373095*T_perp_over_m[0]*tmp_T_perp_g_dist[0]; 
  g_dist_source[1] = 1.264911064067352*T_perp_over_m[1]*tmp_T_perp_g_dist[4]+1.264911064067352*tmp_T_perp_g_dist[1]*T_perp_over_m[2]+1.414213562373095*T_perp_over_m[0]*tmp_T_perp_g_dist[1]+1.414213562373095*tmp_T_perp_g_dist[0]*T_perp_over_m[1]; 
  g_dist_source[2] = 1.414213562373095*T_perp_over_m[2]*tmp_T_perp_g_dist[6]+1.414213562373095*T_perp_over_m[1]*tmp_T_perp_g_dist[3]+1.414213562373095*T_perp_over_m[0]*tmp_T_perp_g_dist[2]; 
  g_dist_source[3] = 1.264911064067352*T_perp_over_m[1]*tmp_T_perp_g_dist[6]+1.264911064067352*T_perp_over_m[2]*tmp_T_perp_g_dist[3]+1.414213562373095*T_perp_over_m[0]*tmp_T_perp_g_dist[3]+1.414213562373095*T_perp_over_m[1]*tmp_T_perp_g_dist[2]; 
  g_dist_source[4] = 0.9035079029052515*T_perp_over_m[2]*tmp_T_perp_g_dist[4]+1.414213562373095*T_perp_over_m[0]*tmp_T_perp_g_dist[4]+1.414213562373095*tmp_T_perp_g_dist[0]*T_perp_over_m[2]+1.264911064067352*T_perp_over_m[1]*tmp_T_perp_g_dist[1]; 
  g_dist_source[5] = 1.414213562373095*T_perp_over_m[1]*tmp_T_perp_g_dist[7]+1.414213562373095*T_perp_over_m[0]*tmp_T_perp_g_dist[5]; 
  g_dist_source[6] = 0.9035079029052515*T_perp_over_m[2]*tmp_T_perp_g_dist[6]+1.414213562373095*T_perp_over_m[0]*tmp_T_perp_g_dist[6]+1.264911064067352*T_perp_over_m[1]*tmp_T_perp_g_dist[3]+1.414213562373095*T_perp_over_m[2]*tmp_T_perp_g_dist[2]; 
  g_dist_source[7] = 1.264911064067352*T_perp_over_m[2]*tmp_T_perp_g_dist[7]+1.414213562373095*T_perp_over_m[0]*tmp_T_perp_g_dist[7]+1.414213562373095*T_perp_over_m[1]*tmp_T_perp_g_dist[5]; 

  tmp_F_0_m_F_1[0] = 0.7071067811865475*T_perp_over_m_inv[2]*G_1[4]+0.7071067811865475*G_1[1]*T_perp_over_m_inv[1]+0.7071067811865475*G_1[0]*T_perp_over_m_inv[0]; 
  tmp_F_0_m_F_1[1] = 0.6324555320336759*T_perp_over_m_inv[1]*G_1[4]+0.6324555320336759*G_1[1]*T_perp_over_m_inv[2]+0.7071067811865475*G_1[0]*T_perp_over_m_inv[1]+0.7071067811865475*T_perp_over_m_inv[0]*G_1[1]; 
  tmp_F_0_m_F_1[2] = 0.7071067811865475*T_perp_over_m_inv[2]*G_1[6]+0.7071067811865475*T_perp_over_m_inv[1]*G_1[3]+0.7071067811865475*T_perp_over_m_inv[0]*G_1[2]; 
  tmp_F_0_m_F_1[3] = 0.632455532033676*T_perp_over_m_inv[1]*G_1[6]+0.6324555320336759*T_perp_over_m_inv[2]*G_1[3]+0.7071067811865475*T_perp_over_m_inv[0]*G_1[3]+0.7071067811865475*T_perp_over_m_inv[1]*G_1[2]; 
  tmp_F_0_m_F_1[4] = 0.4517539514526256*T_perp_over_m_inv[2]*G_1[4]+0.7071067811865475*T_perp_over_m_inv[0]*G_1[4]+0.7071067811865475*G_1[0]*T_perp_over_m_inv[2]+0.6324555320336759*G_1[1]*T_perp_over_m_inv[1]; 
  tmp_F_0_m_F_1[5] = 0.7071067811865475*T_perp_over_m_inv[1]*G_1[7]+0.7071067811865475*T_perp_over_m_inv[0]*G_1[5]; 
  tmp_F_0_m_F_1[6] = 0.4517539514526256*T_perp_over_m_inv[2]*G_1[6]+0.7071067811865475*T_perp_over_m_inv[0]*G_1[6]+0.632455532033676*T_perp_over_m_inv[1]*G_1[3]+0.7071067811865475*G_1[2]*T_perp_over_m_inv[2]; 
  tmp_F_0_m_F_1[7] = 0.6324555320336759*T_perp_over_m_inv[2]*G_1[7]+0.7071067811865475*T_perp_over_m_inv[0]*G_1[7]+0.7071067811865475*T_perp_over_m_inv[1]*G_1[5]; 

  out_F_1[0] = F_0[0]-1.0*tmp_F_0_m_F_1[0]; 
  out_F_1[1] = F_0[1]-1.0*tmp_F_0_m_F_1[1]; 
  out_F_1[2] = F_0[2]-1.0*tmp_F_0_m_F_1[2]; 
  out_F_1[3] = F_0[3]-1.0*tmp_F_0_m_F_1[3]; 
  out_F_1[4] = F_0[4]-1.0*tmp_F_0_m_F_1[4]; 
  out_F_1[5] = F_0[5]-1.0*tmp_F_0_m_F_1[5]; 
  out_F_1[6] = F_0[6]-1.0*tmp_F_0_m_F_1[6]; 
  out_F_1[7] = F_0[7]-1.0*tmp_F_0_m_F_1[7]; 
} 
