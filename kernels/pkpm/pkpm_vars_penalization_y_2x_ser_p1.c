#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_basis_ser_2x_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_ser_2x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void pkpm_vars_penalization_y_2x_ser_p1(double tol, bool force_lax, 
  const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
  const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r,
  const double *p_ij_l, const double *p_ij_r,
  const double *prim_l, const double *prim_r, 
  const double *euler_pkpm_l, const double *euler_pkpm_r,
  double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization) 
{ 
  // tol:                  Tolerance in rho^+, rho^-, and u_avg for switching to Lax fluxes.
  // force_lax:            Flag for forcing Lax fluxes to be turned on.
  // wv_eqn:               Wave equation for computing fluctuations at the interface for upwinding.
  // geom:                 Geometry for the surface update.
  // vlasov_pkpm_moms_l/r: Input pkpm moments to the left/right of the interface.
  // prim_l/r:             Input primitive variables [u_i, 3*T_ii/m] to the left/right of the interface.
  // p_ij_l/r:             Input p_ij = (p_par - p_perp) b_i b_j + p_perp g_ij to the left/right of the interface.
  // euler_pkpm_l/r:       Input [rho ux, rho uy, rho uz], Fluid input state vector to the left/right of the interface.
  // pkpm_lax:             Surface expansion of pkpm Lax penalization: lambda_i = |u_i| + sqrt(3.0*T_ii/m).
  // pkpm_penalization:    Surface expansion of the penalization term in the PKPM momentum update in each direction (cdim components).
  //                       Note: Each cell owns their *lower* edge surface evaluation.

  const double *rho_l = &vlasov_pkpm_moms_l[0]; 
  const double *rho_r = &vlasov_pkpm_moms_r[0]; 

  const double *rhoux_l = &euler_pkpm_l[0]; 
  const double *rhouy_l = &euler_pkpm_l[4]; 
  const double *rhouz_l = &euler_pkpm_l[8]; 

  const double *rhoux_r = &euler_pkpm_r[0]; 
  const double *rhouy_r = &euler_pkpm_r[4]; 
  const double *rhouz_r = &euler_pkpm_r[8]; 

  const double *u_i_l = &prim_l[4]; 
  const double *Tii_l = &prim_l[28]; 

  const double *u_i_r = &prim_r[4]; 
  const double *Tii_r = &prim_r[28]; 

  const double *Pxx_l = &p_ij_l[0]; 
  const double *Pxy_l = &p_ij_l[4]; 
  const double *Pxz_l = &p_ij_l[8]; 
  const double *Pyy_l = &p_ij_l[12]; 
  const double *Pyz_l = &p_ij_l[16]; 
  const double *Pzz_l = &p_ij_l[20]; 

  const double *Pxx_r = &p_ij_r[0]; 
  const double *Pxy_r = &p_ij_r[4]; 
  const double *Pxz_r = &p_ij_r[8]; 
  const double *Pyy_r = &p_ij_r[12]; 
  const double *Pyz_r = &p_ij_r[16]; 
  const double *Pzz_r = &p_ij_r[20]; 

  double *pkpm_lax_l = &pkpm_lax[2]; 

  double *pkpm_penalization_rhoux_l = &pkpm_penalization[6]; 
  double *pkpm_penalization_rhouy_l = &pkpm_penalization[8]; 
  double *pkpm_penalization_rhouz_l = &pkpm_penalization[10]; 

  double amdq_rhoux_quad[2] = {0.0}; 
  double apdq_rhoux_quad[2] = {0.0}; 
  double amdq_rhouy_quad[2] = {0.0}; 
  double apdq_rhouy_quad[2] = {0.0}; 
  double amdq_rhouz_quad[2] = {0.0}; 
  double apdq_rhouz_quad[2] = {0.0}; 
  double pkpm_lax_quad[2] = {0.0}; 

  double q_l[10] = {0.0}; 
  double q_r[10] = {0.0}; 
  double u_l = 0.0; 
  double u_r = 0.0; 
  double T_l = 0.0; 
  double T_r = 0.0; 
  double u_max = 0.0; 
  double vth_max = 0.0; 
  double q_l_local[10] = {0.0}; 
  double q_r_local[10] = {0.0}; 
  double delta[10] = {0.0}; 
  double my_max_speed = 0.0; 
  double lenr = 0.0; 
  double waves[50] = {0.0}; 
  double speeds[5] = {0.0}; 
  double amdq_local[10] = {0.0}; 
  double apdq_local[10] = {0.0}; 
  double amdq[10] = {0.0}; 
  double apdq[10] = {0.0}; 

  int use_lax = 0;
  q_l[0] = ser_2x_p1_surfx2_eval_quad_node_0_r(rho_l); 
  q_l[1] = ser_2x_p1_surfx2_eval_quad_node_0_r(rhoux_l); 
  q_l[2] = ser_2x_p1_surfx2_eval_quad_node_0_r(rhouy_l); 
  q_l[3] = ser_2x_p1_surfx2_eval_quad_node_0_r(rhouz_l); 
  q_l[4] = ser_2x_p1_surfx2_eval_quad_node_0_r(Pxx_l) + q_l[1]*q_l[1]/q_l[0]; 
  q_l[5] = ser_2x_p1_surfx2_eval_quad_node_0_r(Pxy_l) + q_l[1]*q_l[2]/q_l[0]; 
  q_l[6] = ser_2x_p1_surfx2_eval_quad_node_0_r(Pxz_l) + q_l[1]*q_l[3]/q_l[0]; 
  q_l[7] = ser_2x_p1_surfx2_eval_quad_node_0_r(Pyy_l) + q_l[2]*q_l[2]/q_l[0]; 
  q_l[8] = ser_2x_p1_surfx2_eval_quad_node_0_r(Pyz_l) + q_l[2]*q_l[3]/q_l[0]; 
  q_l[9] = ser_2x_p1_surfx2_eval_quad_node_0_r(Pzz_l) + q_l[3]*q_l[3]/q_l[0]; 
  q_r[0] = ser_2x_p1_surfx2_eval_quad_node_0_l(rho_r); 
  q_r[1] = ser_2x_p1_surfx2_eval_quad_node_0_l(rhoux_r); 
  q_r[2] = ser_2x_p1_surfx2_eval_quad_node_0_l(rhouy_r); 
  q_r[3] = ser_2x_p1_surfx2_eval_quad_node_0_l(rhouz_r); 
  q_r[4] = ser_2x_p1_surfx2_eval_quad_node_0_l(Pxx_r) + q_r[1]*q_r[1]/q_r[0]; 
  q_r[5] = ser_2x_p1_surfx2_eval_quad_node_0_l(Pxy_r) + q_r[1]*q_r[2]/q_r[0]; 
  q_r[6] = ser_2x_p1_surfx2_eval_quad_node_0_l(Pxz_r) + q_r[1]*q_r[3]/q_r[0]; 
  q_r[7] = ser_2x_p1_surfx2_eval_quad_node_0_l(Pyy_r) + q_r[2]*q_r[2]/q_r[0]; 
  q_r[8] = ser_2x_p1_surfx2_eval_quad_node_0_l(Pyz_r) + q_r[2]*q_r[3]/q_r[0]; 
  q_r[9] = ser_2x_p1_surfx2_eval_quad_node_0_l(Pzz_r) + q_r[3]*q_r[3]/q_r[0]; 

  T_l = ser_2x_p1_surfx2_eval_quad_node_0_r(Tii_l); 
  T_r = ser_2x_p1_surfx2_eval_quad_node_0_l(Tii_r); 
  if (T_l > 0.0 && T_r > 0.0) { 
    vth_max = fmax(sqrt(T_l), sqrt(T_r)); 
  } else if (T_l > 0.0 && T_r < 0.0) { 
    vth_max = sqrt(T_l); 
    use_lax = 1; 
  } else if (T_l < 0.0 && T_r > 0.0) { 
    vth_max = sqrt(T_r); 
    use_lax = 1; 
  } else { 
    vth_max = 0.0; 
    use_lax = 1; 
  } 
  u_l = ser_2x_p1_surfx2_eval_quad_node_0_r(u_i_l); 
  u_r = ser_2x_p1_surfx2_eval_quad_node_0_l(u_i_r); 
  u_max = fmax(fabs(u_l), fabs(u_r)); 
  pkpm_lax_quad[0] = u_max + vth_max; 

  if (q_l[0] < tol) use_lax = 1; 
  if (q_r[0] < tol) use_lax = 1; 

  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom->tau1[1], geom->tau2[1], geom->norm[1], q_l, q_l_local); 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom->tau1[1], geom->tau2[1], geom->norm[1], q_r, q_r_local); 

  delta[0] = q_r_local[0] - q_l_local[0]; 
  delta[1] = q_r_local[1] - q_l_local[1]; 
  delta[2] = q_r_local[2] - q_l_local[2]; 
  delta[3] = q_r_local[3] - q_l_local[3]; 
  delta[4] = q_r_local[4] - q_l_local[4]; 
  delta[5] = q_r_local[5] - q_l_local[5]; 
  delta[6] = q_r_local[6] - q_l_local[6]; 
  delta[7] = q_r_local[7] - q_l_local[7]; 
  delta[8] = q_r_local[8] - q_l_local[8]; 
  delta[9] = q_r_local[9] - q_l_local[9]; 
  my_max_speed = gkyl_wv_eqn_waves(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta, q_l_local, q_r_local, waves, speeds); 
  lenr = geom->lenr[1]; 
  speeds[0] *= lenr; 
  speeds[1] *= lenr; 
  speeds[2] *= lenr; 
  speeds[3] *= lenr; 
  speeds[4] *= lenr; 
  gkyl_wv_eqn_qfluct(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, q_l_local, q_r_local, waves, speeds, amdq_local, apdq_local); 
  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom->tau1[1], geom->tau2[1], geom->norm[1], amdq_local, amdq); 
  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom->tau1[1], geom->tau2[1], geom->norm[1], apdq_local, apdq); 

  amdq_rhoux_quad[0] = amdq[1]; 
  apdq_rhoux_quad[0] = apdq[1]; 
  amdq_rhouy_quad[0] = amdq[2]; 
  apdq_rhouy_quad[0] = apdq[2]; 
  amdq_rhouz_quad[0] = amdq[3]; 
  apdq_rhouz_quad[0] = apdq[3]; 

  q_l[0] = ser_2x_p1_surfx2_eval_quad_node_1_r(rho_l); 
  q_l[1] = ser_2x_p1_surfx2_eval_quad_node_1_r(rhoux_l); 
  q_l[2] = ser_2x_p1_surfx2_eval_quad_node_1_r(rhouy_l); 
  q_l[3] = ser_2x_p1_surfx2_eval_quad_node_1_r(rhouz_l); 
  q_l[4] = ser_2x_p1_surfx2_eval_quad_node_1_r(Pxx_l) + q_l[1]*q_l[1]/q_l[0]; 
  q_l[5] = ser_2x_p1_surfx2_eval_quad_node_1_r(Pxy_l) + q_l[1]*q_l[2]/q_l[0]; 
  q_l[6] = ser_2x_p1_surfx2_eval_quad_node_1_r(Pxz_l) + q_l[1]*q_l[3]/q_l[0]; 
  q_l[7] = ser_2x_p1_surfx2_eval_quad_node_1_r(Pyy_l) + q_l[2]*q_l[2]/q_l[0]; 
  q_l[8] = ser_2x_p1_surfx2_eval_quad_node_1_r(Pyz_l) + q_l[2]*q_l[3]/q_l[0]; 
  q_l[9] = ser_2x_p1_surfx2_eval_quad_node_1_r(Pzz_l) + q_l[3]*q_l[3]/q_l[0]; 
  q_r[0] = ser_2x_p1_surfx2_eval_quad_node_1_l(rho_r); 
  q_r[1] = ser_2x_p1_surfx2_eval_quad_node_1_l(rhoux_r); 
  q_r[2] = ser_2x_p1_surfx2_eval_quad_node_1_l(rhouy_r); 
  q_r[3] = ser_2x_p1_surfx2_eval_quad_node_1_l(rhouz_r); 
  q_r[4] = ser_2x_p1_surfx2_eval_quad_node_1_l(Pxx_r) + q_r[1]*q_r[1]/q_r[0]; 
  q_r[5] = ser_2x_p1_surfx2_eval_quad_node_1_l(Pxy_r) + q_r[1]*q_r[2]/q_r[0]; 
  q_r[6] = ser_2x_p1_surfx2_eval_quad_node_1_l(Pxz_r) + q_r[1]*q_r[3]/q_r[0]; 
  q_r[7] = ser_2x_p1_surfx2_eval_quad_node_1_l(Pyy_r) + q_r[2]*q_r[2]/q_r[0]; 
  q_r[8] = ser_2x_p1_surfx2_eval_quad_node_1_l(Pyz_r) + q_r[2]*q_r[3]/q_r[0]; 
  q_r[9] = ser_2x_p1_surfx2_eval_quad_node_1_l(Pzz_r) + q_r[3]*q_r[3]/q_r[0]; 

  T_l = ser_2x_p1_surfx2_eval_quad_node_1_r(Tii_l); 
  T_r = ser_2x_p1_surfx2_eval_quad_node_1_l(Tii_r); 
  if (T_l > 0.0 && T_r > 0.0) { 
    vth_max = fmax(sqrt(T_l), sqrt(T_r)); 
  } else if (T_l > 0.0 && T_r < 0.0) { 
    vth_max = sqrt(T_l); 
    use_lax = 1; 
  } else if (T_l < 0.0 && T_r > 0.0) { 
    vth_max = sqrt(T_r); 
    use_lax = 1; 
  } else { 
    vth_max = 0.0; 
    use_lax = 1; 
  } 
  u_l = ser_2x_p1_surfx2_eval_quad_node_1_r(u_i_l); 
  u_r = ser_2x_p1_surfx2_eval_quad_node_1_l(u_i_r); 
  u_max = fmax(fabs(u_l), fabs(u_r)); 
  pkpm_lax_quad[1] = u_max + vth_max; 

  if (q_l[0] < tol) use_lax = 1; 
  if (q_r[0] < tol) use_lax = 1; 

  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom->tau1[1], geom->tau2[1], geom->norm[1], q_l, q_l_local); 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom->tau1[1], geom->tau2[1], geom->norm[1], q_r, q_r_local); 

  delta[0] = q_r_local[0] - q_l_local[0]; 
  delta[1] = q_r_local[1] - q_l_local[1]; 
  delta[2] = q_r_local[2] - q_l_local[2]; 
  delta[3] = q_r_local[3] - q_l_local[3]; 
  delta[4] = q_r_local[4] - q_l_local[4]; 
  delta[5] = q_r_local[5] - q_l_local[5]; 
  delta[6] = q_r_local[6] - q_l_local[6]; 
  delta[7] = q_r_local[7] - q_l_local[7]; 
  delta[8] = q_r_local[8] - q_l_local[8]; 
  delta[9] = q_r_local[9] - q_l_local[9]; 
  my_max_speed = gkyl_wv_eqn_waves(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta, q_l_local, q_r_local, waves, speeds); 
  lenr = geom->lenr[1]; 
  speeds[0] *= lenr; 
  speeds[1] *= lenr; 
  speeds[2] *= lenr; 
  speeds[3] *= lenr; 
  speeds[4] *= lenr; 
  gkyl_wv_eqn_qfluct(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, q_l_local, q_r_local, waves, speeds, amdq_local, apdq_local); 
  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom->tau1[1], geom->tau2[1], geom->norm[1], amdq_local, amdq); 
  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom->tau1[1], geom->tau2[1], geom->norm[1], apdq_local, apdq); 

  amdq_rhoux_quad[1] = amdq[1]; 
  apdq_rhoux_quad[1] = apdq[1]; 
  amdq_rhouy_quad[1] = amdq[2]; 
  apdq_rhouy_quad[1] = apdq[2]; 
  amdq_rhouz_quad[1] = amdq[3]; 
  apdq_rhouz_quad[1] = apdq[3]; 

  ser_2x_p1_upwind_quad_to_modal(pkpm_lax_quad, pkpm_lax_l); 
  if (use_lax || force_lax) { 
  pkpm_penalization_rhoux_l[0] = (-0.4330127018922193*pkpm_lax_l[1]*rhoux_r[3])-0.4330127018922193*pkpm_lax_l[1]*rhoux_l[3]-0.4330127018922193*pkpm_lax_l[0]*rhoux_r[2]-0.4330127018922193*pkpm_lax_l[0]*rhoux_l[2]+0.25*pkpm_lax_l[1]*rhoux_r[1]-0.25*pkpm_lax_l[1]*rhoux_l[1]+0.25*pkpm_lax_l[0]*rhoux_r[0]-0.25*pkpm_lax_l[0]*rhoux_l[0]; 
  pkpm_penalization_rhoux_l[1] = (-0.4330127018922193*pkpm_lax_l[0]*rhoux_r[3])-0.4330127018922193*pkpm_lax_l[0]*rhoux_l[3]-0.4330127018922193*pkpm_lax_l[1]*rhoux_r[2]-0.4330127018922193*pkpm_lax_l[1]*rhoux_l[2]+0.25*pkpm_lax_l[0]*rhoux_r[1]-0.25*pkpm_lax_l[0]*rhoux_l[1]+0.25*rhoux_r[0]*pkpm_lax_l[1]-0.25*rhoux_l[0]*pkpm_lax_l[1]; 

  pkpm_penalization_rhouy_l[0] = (-0.4330127018922193*pkpm_lax_l[1]*rhouy_r[3])-0.4330127018922193*pkpm_lax_l[1]*rhouy_l[3]-0.4330127018922193*pkpm_lax_l[0]*rhouy_r[2]-0.4330127018922193*pkpm_lax_l[0]*rhouy_l[2]+0.25*pkpm_lax_l[1]*rhouy_r[1]-0.25*pkpm_lax_l[1]*rhouy_l[1]+0.25*pkpm_lax_l[0]*rhouy_r[0]-0.25*pkpm_lax_l[0]*rhouy_l[0]; 
  pkpm_penalization_rhouy_l[1] = (-0.4330127018922193*pkpm_lax_l[0]*rhouy_r[3])-0.4330127018922193*pkpm_lax_l[0]*rhouy_l[3]-0.4330127018922193*pkpm_lax_l[1]*rhouy_r[2]-0.4330127018922193*pkpm_lax_l[1]*rhouy_l[2]+0.25*pkpm_lax_l[0]*rhouy_r[1]-0.25*pkpm_lax_l[0]*rhouy_l[1]+0.25*rhouy_r[0]*pkpm_lax_l[1]-0.25*rhouy_l[0]*pkpm_lax_l[1]; 

  pkpm_penalization_rhouz_l[0] = (-0.4330127018922193*pkpm_lax_l[1]*rhouz_r[3])-0.4330127018922193*pkpm_lax_l[1]*rhouz_l[3]-0.4330127018922193*pkpm_lax_l[0]*rhouz_r[2]-0.4330127018922193*pkpm_lax_l[0]*rhouz_l[2]+0.25*pkpm_lax_l[1]*rhouz_r[1]-0.25*pkpm_lax_l[1]*rhouz_l[1]+0.25*pkpm_lax_l[0]*rhouz_r[0]-0.25*pkpm_lax_l[0]*rhouz_l[0]; 
  pkpm_penalization_rhouz_l[1] = (-0.4330127018922193*pkpm_lax_l[0]*rhouz_r[3])-0.4330127018922193*pkpm_lax_l[0]*rhouz_l[3]-0.4330127018922193*pkpm_lax_l[1]*rhouz_r[2]-0.4330127018922193*pkpm_lax_l[1]*rhouz_l[2]+0.25*pkpm_lax_l[0]*rhouz_r[1]-0.25*pkpm_lax_l[0]*rhouz_l[1]+0.25*rhouz_r[0]*pkpm_lax_l[1]-0.25*rhouz_l[0]*pkpm_lax_l[1]; 

  } else { 
    double amdq_rhoux[2] = {0.0}; 
    double apdq_rhoux[2] = {0.0}; 
    double amdq_rhouy[2] = {0.0}; 
    double apdq_rhouy[2] = {0.0}; 
    double amdq_rhouz[2] = {0.0}; 
    double apdq_rhouz[2] = {0.0}; 

    ser_2x_p1_upwind_quad_to_modal(amdq_rhoux_quad, amdq_rhoux); 
    ser_2x_p1_upwind_quad_to_modal(amdq_rhouy_quad, amdq_rhouy); 
    ser_2x_p1_upwind_quad_to_modal(amdq_rhouz_quad, amdq_rhouz); 

    ser_2x_p1_upwind_quad_to_modal(apdq_rhoux_quad, apdq_rhoux); 
    ser_2x_p1_upwind_quad_to_modal(apdq_rhouy_quad, apdq_rhouy); 
    ser_2x_p1_upwind_quad_to_modal(apdq_rhouz_quad, apdq_rhouz); 

    pkpm_penalization_rhoux_l[0] = 0.5*(apdq_rhoux[0] - amdq_rhoux[0]); 
    pkpm_penalization_rhouy_l[0] = 0.5*(apdq_rhouy[0] - amdq_rhouy[0]); 
    pkpm_penalization_rhouz_l[0] = 0.5*(apdq_rhouz[0] - amdq_rhouz[0]); 

    pkpm_penalization_rhoux_l[1] = 0.5*(apdq_rhoux[1] - amdq_rhoux[1]); 
    pkpm_penalization_rhouy_l[1] = 0.5*(apdq_rhouy[1] - amdq_rhouy[1]); 
    pkpm_penalization_rhouz_l[1] = 0.5*(apdq_rhouz[1] - amdq_rhouz[1]); 

  } 
} 
