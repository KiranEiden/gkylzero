#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void euler_pkpm_limiter_x_1x_ser_p2(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
  const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
  const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
  const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
  double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r)
{ 
  // limiter_fac:            Factor for relationship between cell slopes and cell average differences (by default: 1/sqrt(3)).
  // wv_eqn:                 Wave equation for computing waves for limiting characteristics.
  // geom:                   Geometry on the left (*only works with Cartesian components*).
  // prim_c:                 Input flow velocity in the center cell.
  // vlasov_pkpm_moms_l/c/r: Input pkpm moments in left/center/right cells.
  // p_ij_l/c/r:             Input volume expansion of p_ij = (p_par - p_perp) b_i b_j + p_perp g_ij in left/center/right cells.
  // euler_pkpm_l/c/r:       [rho ux, rho uy, rho uz], Fluid input and output (after limiting) state vector in left/center/right cells..

  const double *norm = geom->norm[0]; 
  const double *tau1 = geom->tau1[0]; 
  const double *tau2 = geom->tau2[0]; 

  const double *ux_c = &prim_c[0]; 
  const double *uy_c = &prim_c[3]; 
  const double *uz_c = &prim_c[6]; 

  const double *rho_l = &vlasov_pkpm_moms_l[0]; 
  const double *rho_c = &vlasov_pkpm_moms_c[0]; 
  const double *rho_r = &vlasov_pkpm_moms_r[0]; 

  const double *Pxx_l = &p_ij_l[0]; 
  const double *Pxy_l = &p_ij_l[3]; 
  const double *Pxz_l = &p_ij_l[6]; 
  const double *Pyy_l = &p_ij_l[9]; 
  const double *Pyz_l = &p_ij_l[12]; 
  const double *Pzz_l = &p_ij_l[15]; 

  const double *Pxx_c = &p_ij_c[0]; 
  const double *Pxy_c = &p_ij_c[3]; 
  const double *Pxz_c = &p_ij_c[6]; 
  const double *Pyy_c = &p_ij_c[9]; 
  const double *Pyz_c = &p_ij_c[12]; 
  const double *Pzz_c = &p_ij_c[15]; 

  const double *Pxx_r = &p_ij_r[0]; 
  const double *Pxy_r = &p_ij_r[3]; 
  const double *Pxz_r = &p_ij_r[6]; 
  const double *Pyy_r = &p_ij_r[9]; 
  const double *Pyz_r = &p_ij_r[12]; 
  const double *Pzz_r = &p_ij_r[15]; 

  double *rhoux_l = &euler_pkpm_l[0]; 
  double *rhouy_l = &euler_pkpm_l[3]; 
  double *rhouz_l = &euler_pkpm_l[6]; 

  double *rhoux_c = &euler_pkpm_c[0]; 
  double *rhouy_c = &euler_pkpm_c[3]; 
  double *rhouz_c = &euler_pkpm_c[6]; 

  double *rhoux_r = &euler_pkpm_r[0]; 
  double *rhouy_r = &euler_pkpm_r[3]; 
  double *rhouz_r = &euler_pkpm_r[6]; 

  double q_avg_l[10] = {0.0}; 
  q_avg_l[0] = rho_l[0]; 
  q_avg_l[1] = rhoux_l[0]; 
  q_avg_l[2] = rhouy_l[0]; 
  q_avg_l[3] = rhouz_l[0]; 
  q_avg_l[4] = Pxx_l[0] + q_avg_l[1]*q_avg_l[1]/q_avg_l[0]; 
  q_avg_l[5] = Pxy_l[0] + q_avg_l[1]*q_avg_l[2]/q_avg_l[0]; 
  q_avg_l[6] = Pxz_l[0] + q_avg_l[1]*q_avg_l[3]/q_avg_l[0]; 
  q_avg_l[7] = Pyy_l[0] + q_avg_l[2]*q_avg_l[2]/q_avg_l[0]; 
  q_avg_l[8] = Pyz_l[0] + q_avg_l[2]*q_avg_l[3]/q_avg_l[0]; 
  q_avg_l[9] = Pzz_l[0] + q_avg_l[3]*q_avg_l[3]/q_avg_l[0]; 
  double q_avg_c[10] = {0.0}; 
  q_avg_c[0] = rho_c[0]; 
  q_avg_c[1] = rhoux_c[0]; 
  q_avg_c[2] = rhouy_c[0]; 
  q_avg_c[3] = rhouz_c[0]; 
  q_avg_c[4] = Pxx_c[0] + q_avg_c[1]*q_avg_c[1]/q_avg_c[0]; 
  q_avg_c[5] = Pxy_c[0] + q_avg_c[1]*q_avg_c[2]/q_avg_c[0]; 
  q_avg_c[6] = Pxz_c[0] + q_avg_c[1]*q_avg_c[3]/q_avg_c[0]; 
  q_avg_c[7] = Pyy_c[0] + q_avg_c[2]*q_avg_c[2]/q_avg_c[0]; 
  q_avg_c[8] = Pyz_c[0] + q_avg_c[2]*q_avg_c[3]/q_avg_c[0]; 
  q_avg_c[9] = Pzz_c[0] + q_avg_c[3]*q_avg_c[3]/q_avg_c[0]; 
  double q_avg_r[10] = {0.0}; 
  q_avg_r[0] = rho_r[0]; 
  q_avg_r[1] = rhoux_r[0]; 
  q_avg_r[2] = rhouy_r[0]; 
  q_avg_r[3] = rhouz_r[0]; 
  q_avg_r[4] = Pxx_r[0] + q_avg_r[1]*q_avg_r[1]/q_avg_r[0]; 
  q_avg_r[5] = Pxy_r[0] + q_avg_r[1]*q_avg_r[2]/q_avg_r[0]; 
  q_avg_r[6] = Pxz_r[0] + q_avg_r[1]*q_avg_r[3]/q_avg_r[0]; 
  q_avg_r[7] = Pyy_r[0] + q_avg_r[2]*q_avg_r[2]/q_avg_r[0]; 
  q_avg_r[8] = Pyz_r[0] + q_avg_r[2]*q_avg_r[3]/q_avg_r[0]; 
  q_avg_r[9] = Pzz_r[0] + q_avg_r[3]*q_avg_r[3]/q_avg_r[0]; 
  double q_l_local[10] = {0.0}; 
  double q_c_local[10] = {0.0}; 
  double q_r_local[10] = {0.0}; 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, tau1, tau2, norm, q_avg_l, q_l_local); 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, tau1, tau2, norm, q_avg_c, q_c_local); 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, tau1, tau2, norm, q_avg_r, q_r_local); 

  double delta_l[10] = {0.0}; 
  double delta_r[10] = {0.0}; 
  delta_l[0] = limiter_fac*(q_c_local[0] - q_l_local[0]); 
  delta_l[1] = limiter_fac*(q_c_local[1] - q_l_local[1]); 
  delta_l[2] = limiter_fac*(q_c_local[2] - q_l_local[2]); 
  delta_l[3] = limiter_fac*(q_c_local[3] - q_l_local[3]); 
  delta_l[4] = limiter_fac*(q_c_local[4] - q_l_local[4]); 
  delta_l[5] = limiter_fac*(q_c_local[5] - q_l_local[5]); 
  delta_l[6] = limiter_fac*(q_c_local[6] - q_l_local[6]); 
  delta_l[7] = limiter_fac*(q_c_local[7] - q_l_local[7]); 
  delta_l[8] = limiter_fac*(q_c_local[8] - q_l_local[8]); 
  delta_l[9] = limiter_fac*(q_c_local[9] - q_l_local[9]); 
  delta_r[0] = limiter_fac*(q_r_local[0] - q_c_local[0]); 
  delta_r[1] = limiter_fac*(q_r_local[1] - q_c_local[1]); 
  delta_r[2] = limiter_fac*(q_r_local[2] - q_c_local[2]); 
  delta_r[3] = limiter_fac*(q_r_local[3] - q_c_local[3]); 
  delta_r[4] = limiter_fac*(q_r_local[4] - q_c_local[4]); 
  delta_r[5] = limiter_fac*(q_r_local[5] - q_c_local[5]); 
  delta_r[6] = limiter_fac*(q_r_local[6] - q_c_local[6]); 
  delta_r[7] = limiter_fac*(q_r_local[7] - q_c_local[7]); 
  delta_r[8] = limiter_fac*(q_r_local[8] - q_c_local[8]); 
  delta_r[9] = limiter_fac*(q_r_local[9] - q_c_local[9]); 

  double delta_c[10] = {0.0}; 
  delta_c[0] = rho_c[1]; 
  delta_c[1] = rhoux_c[1]; 
  delta_c[2] = rhouy_c[1]; 
  delta_c[3] = rhouz_c[1]; 
  delta_c[4] = 0.6324555320336759*rhoux_c[1]*ux_c[2]+0.6324555320336759*ux_c[1]*rhoux_c[2]+0.7071067811865475*rhoux_c[0]*ux_c[1]+0.7071067811865475*ux_c[0]*rhoux_c[1]+Pxx_c[1]; 
  delta_c[5] = 0.6324555320336759*rhoux_c[1]*uy_c[2]+0.6324555320336759*uy_c[1]*rhoux_c[2]+0.7071067811865475*rhoux_c[0]*uy_c[1]+0.7071067811865475*uy_c[0]*rhoux_c[1]+Pxy_c[1]; 
  delta_c[6] = 0.6324555320336759*rhoux_c[1]*uz_c[2]+0.6324555320336759*uz_c[1]*rhoux_c[2]+0.7071067811865475*rhoux_c[0]*uz_c[1]+0.7071067811865475*uz_c[0]*rhoux_c[1]+Pxz_c[1]; 
  delta_c[7] = 0.6324555320336759*rhouy_c[1]*uy_c[2]+0.6324555320336759*uy_c[1]*rhouy_c[2]+0.7071067811865475*rhouy_c[0]*uy_c[1]+0.7071067811865475*uy_c[0]*rhouy_c[1]+Pyy_c[1]; 
  delta_c[8] = 0.6324555320336759*rhouy_c[1]*uz_c[2]+0.6324555320336759*uz_c[1]*rhouy_c[2]+0.7071067811865475*rhouy_c[0]*uz_c[1]+0.7071067811865475*uz_c[0]*rhouy_c[1]+Pyz_c[1]; 
  delta_c[9] = 0.6324555320336759*rhouz_c[1]*uz_c[2]+0.6324555320336759*uz_c[1]*rhouz_c[2]+0.7071067811865475*rhouz_c[0]*uz_c[1]+0.7071067811865475*uz_c[0]*rhouz_c[1]+Pzz_c[1]; 
  double delta_c_local[10] = {0.0}; 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, tau1, tau2, norm, delta_c, delta_c_local); 

  double waves_slope_l[50] = {0.0}; 
  double waves_slope_c[50] = {0.0}; 
  double waves_slope_r[50] = {0.0}; 
  double speeds[5] = {0.0}; 
  double my_max_speed_l = gkyl_wv_eqn_waves(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta_l, q_c_local, q_c_local, waves_slope_l, speeds); 
  double my_max_speed_c = gkyl_wv_eqn_waves(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta_c_local, q_c_local, q_c_local, waves_slope_c, speeds); 
  double my_max_speed_r = gkyl_wv_eqn_waves(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta_r, q_c_local, q_c_local, waves_slope_r, speeds); 

  double mm[15] = {0.0}; 
  double slope[3] = {0.0}; 
  mm[0] = minmod(waves_slope_c[1], waves_slope_l[1], waves_slope_r[1]); 
  mm[1] = minmod(waves_slope_c[2], waves_slope_l[2], waves_slope_r[2]); 
  mm[2] = minmod(waves_slope_c[3], waves_slope_l[3], waves_slope_r[3]); 
  slope[0] += mm[0]; 
  slope[1] += mm[1]; 
  slope[2] += mm[2]; 

  mm[3] = minmod(waves_slope_c[11], waves_slope_l[11], waves_slope_r[11]); 
  mm[4] = minmod(waves_slope_c[12], waves_slope_l[12], waves_slope_r[12]); 
  mm[5] = minmod(waves_slope_c[13], waves_slope_l[13], waves_slope_r[13]); 
  slope[0] += mm[3]; 
  slope[1] += mm[4]; 
  slope[2] += mm[5]; 

  mm[6] = minmod(waves_slope_c[21], waves_slope_l[21], waves_slope_r[21]); 
  mm[7] = minmod(waves_slope_c[22], waves_slope_l[22], waves_slope_r[22]); 
  mm[8] = minmod(waves_slope_c[23], waves_slope_l[23], waves_slope_r[23]); 
  slope[0] += mm[6]; 
  slope[1] += mm[7]; 
  slope[2] += mm[8]; 

  mm[9] = minmod(waves_slope_c[31], waves_slope_l[31], waves_slope_r[31]); 
  mm[10] = minmod(waves_slope_c[32], waves_slope_l[32], waves_slope_r[32]); 
  mm[11] = minmod(waves_slope_c[33], waves_slope_l[33], waves_slope_r[33]); 
  slope[0] += mm[9]; 
  slope[1] += mm[10]; 
  slope[2] += mm[11]; 

  mm[12] = minmod(waves_slope_c[41], waves_slope_l[41], waves_slope_r[41]); 
  mm[13] = minmod(waves_slope_c[42], waves_slope_l[42], waves_slope_r[42]); 
  mm[14] = minmod(waves_slope_c[43], waves_slope_l[43], waves_slope_r[43]); 
  slope[0] += mm[12]; 
  slope[1] += mm[13]; 
  slope[2] += mm[14]; 

  // Rotate limited slope back to global coordinates 
  rhoux_c[1] = slope[0]*norm[0] + slope[1]*tau1[0] + slope[2]*tau2[0]; 
  rhouy_c[1] = slope[0]*norm[1] + slope[1]*tau1[1] + slope[2]*tau2[1]; 
  rhouz_c[1] = slope[0]*norm[2] + slope[1]*tau1[2] + slope[2]*tau2[2]; 
  for (int i = 0; i < 5; ++i) { 
    if (mm[3*i] != waves_slope_c[10*i+1]) { 
      rhoux_c[2] = 0.0; 
    } 
    if (mm[3*i+1] != waves_slope_c[10*i+2]) { 
      rhouy_c[2] = 0.0; 
    } 
    if (mm[3*i+2] != waves_slope_c[10*i+3]) { 
      rhoux_c[2] = 0.0; 
    } 
  } 
} 
