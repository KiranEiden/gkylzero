#include <gkyl_euler_kernels.h> 
#include <gkyl_basis_tensor_2x_p2_surfx2_eval_quad.h> 
#include <gkyl_basis_tensor_2x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double euler_surfy_2x_tensor_p2(const double *w, const double *dxv, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r, 
    const double *u_surf_l, const double *u_surf_c, const double *u_surf_r,
    const double *p_surf_l, const double *p_surf_c, const double *p_surf_r,
    const double *fluid_l, const double *fluid_c, const double *fluid_r, 
    double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:      Cell-center coordinates.
  // dxv[NDIM]:    Cell spacing.
  // wv_eqn:       Wave equation for computing fluctuations at the interface for upwinding.
  // geom_l:       Geometry for the left surface update.
  // geom_r:       Geometry for the right surface update.
  // u_surf_l/c/r: Input surface expansion of flow velocity in left/center/right cells in each direction.
  // p_surf_l/c/r: Input surface expansion of pressure in left/center/right cells in each direction.
  //               [p_xl, p_xr, p_yl, p_yr, p_zl, p_zr] 
  // fluid_l/c/r:  [rho, rho ux, rho uy, rho uz, energy], Fluid input state vector in left/center/right cells.
  // out:          Incremented output.

  const double dx1 = 2.0/dxv[1]; 

  const double *rho_l = &fluid_l[0]; 
  const double *rhoux_l = &fluid_l[9]; 
  const double *rhouy_l = &fluid_l[18]; 
  const double *rhouz_l = &fluid_l[27]; 
  const double *energy_l = &fluid_l[36]; 

  const double *rho_c = &fluid_c[0]; 
  const double *rhoux_c = &fluid_c[9]; 
  const double *rhouy_c = &fluid_c[18]; 
  const double *rhouz_c = &fluid_c[27]; 
  const double *energy_c = &fluid_c[36]; 

  const double *rho_r = &fluid_r[0]; 
  const double *rhoux_r = &fluid_r[9]; 
  const double *rhouy_r = &fluid_r[18]; 
  const double *rhouz_r = &fluid_r[27]; 
  const double *energy_r = &fluid_r[36]; 

  const double *ux_surf_lr = &u_surf_l[21]; 
  const double *uy_surf_lr = &u_surf_l[27]; 
  const double *uz_surf_lr = &u_surf_l[33]; 

  const double *ux_surf_cl = &u_surf_c[18]; 
  const double *uy_surf_cl = &u_surf_c[24]; 
  const double *uz_surf_cl = &u_surf_c[30]; 

  const double *ux_surf_cr = &u_surf_c[21]; 
  const double *uy_surf_cr = &u_surf_c[27]; 
  const double *uz_surf_cr = &u_surf_c[33]; 

  const double *ux_surf_rl = &u_surf_r[18]; 
  const double *uy_surf_rl = &u_surf_r[24]; 
  const double *uz_surf_rl = &u_surf_r[30]; 

  const double *p_surf_lr = &p_surf_l[9]; 
  const double *p_surf_cl = &p_surf_c[6]; 
  const double *p_surf_cr = &p_surf_c[9]; 
  const double *p_surf_rl = &p_surf_r[6]; 

  double *outrho = &out[0]; 
  double *outrhoux = &out[9]; 
  double *outrhouy = &out[18]; 
  double *outrhouz = &out[27]; 
  double *outenergy = &out[36]; 

  double amdq_rho_l[3] = {0.0}; 
  double apdq_rho_l[3] = {0.0}; 
  double amdq_rhoux_l[3] = {0.0}; 
  double apdq_rhoux_l[3] = {0.0}; 
  double amdq_rhouy_l[3] = {0.0}; 
  double apdq_rhouy_l[3] = {0.0}; 
  double amdq_rhouz_l[3] = {0.0}; 
  double apdq_rhouz_l[3] = {0.0}; 
  double amdq_energy_l[3] = {0.0}; 
  double apdq_energy_l[3] = {0.0}; 

  double amdq_rho_r[3] = {0.0}; 
  double apdq_rho_r[3] = {0.0}; 
  double amdq_rhoux_r[3] = {0.0}; 
  double apdq_rhoux_r[3] = {0.0}; 
  double amdq_rhouy_r[3] = {0.0}; 
  double apdq_rhouy_r[3] = {0.0}; 
  double amdq_rhouz_r[3] = {0.0}; 
  double apdq_rhouz_r[3] = {0.0}; 
  double amdq_energy_r[3] = {0.0}; 
  double apdq_energy_r[3] = {0.0}; 

  double amdq_rho_quad_l[3] = {0.0}; 
  double apdq_rho_quad_l[3] = {0.0}; 
  double amdq_rhoux_quad_l[3] = {0.0}; 
  double apdq_rhoux_quad_l[3] = {0.0}; 
  double amdq_rhouy_quad_l[3] = {0.0}; 
  double apdq_rhouy_quad_l[3] = {0.0}; 
  double amdq_rhouz_quad_l[3] = {0.0}; 
  double apdq_rhouz_quad_l[3] = {0.0}; 
  double amdq_energy_quad_l[3] = {0.0}; 
  double apdq_energy_quad_l[3] = {0.0}; 

  double amdq_rho_quad_r[3] = {0.0}; 
  double apdq_rho_quad_r[3] = {0.0}; 
  double amdq_rhoux_quad_r[3] = {0.0}; 
  double apdq_rhoux_quad_r[3] = {0.0}; 
  double amdq_rhouy_quad_r[3] = {0.0}; 
  double apdq_rhouy_quad_r[3] = {0.0}; 
  double amdq_rhouz_quad_r[3] = {0.0}; 
  double apdq_rhouz_quad_r[3] = {0.0}; 
  double amdq_energy_quad_r[3] = {0.0}; 
  double apdq_energy_quad_r[3] = {0.0}; 

  double q_lr[5] = {0.0}; 
  double q_cl[5] = {0.0}; 
  double q_cr[5] = {0.0}; 
  double q_rl[5] = {0.0}; 
  double q_lr_local[5] = {0.0}; 
  double q_cl_local[5] = {0.0}; 
  double q_cr_local[5] = {0.0}; 
  double q_rl_local[5] = {0.0}; 
  double delta_l[5] = {0.0}; 
  double delta_r[5] = {0.0}; 
  double my_max_speed_l = 0.0; 
  double my_max_speed_r = 0.0; 
  double lenr_l = 0.0; 
  double lenr_r = 0.0; 
  double waves_l[15] = {0.0}; 
  double waves_r[15] = {0.0}; 
  double speeds_l[3] = {0.0}; 
  double speeds_r[3] = {0.0}; 
  double amdq_l_local[5] = {0.0}; 
  double apdq_l_local[5] = {0.0}; 
  double amdq_r_local[5] = {0.0}; 
  double apdq_r_local[5] = {0.0}; 
  double amdq_l[5] = {0.0}; 
  double apdq_l[5] = {0.0}; 
  double amdq_r[5] = {0.0}; 
  double apdq_r[5] = {0.0}; 

  q_lr[0] = tensor_2x_p2_surfx2_eval_quad_node_0_r(rho_l); 
  q_lr[1] = tensor_2x_p2_surfx2_eval_quad_node_0_r(rhoux_l); 
  q_lr[2] = tensor_2x_p2_surfx2_eval_quad_node_0_r(rhouy_l); 
  q_lr[3] = tensor_2x_p2_surfx2_eval_quad_node_0_r(rhouz_l); 
  q_lr[4] = tensor_2x_p2_surfx2_eval_quad_node_0_r(energy_l); 
  q_cl[0] = tensor_2x_p2_surfx2_eval_quad_node_0_l(rho_c); 
  q_cl[1] = tensor_2x_p2_surfx2_eval_quad_node_0_l(rhoux_c); 
  q_cl[2] = tensor_2x_p2_surfx2_eval_quad_node_0_l(rhouy_c); 
  q_cl[3] = tensor_2x_p2_surfx2_eval_quad_node_0_l(rhouz_c); 
  q_cl[4] = tensor_2x_p2_surfx2_eval_quad_node_0_l(energy_c); 
  q_cr[0] = tensor_2x_p2_surfx2_eval_quad_node_0_r(rho_c); 
  q_cr[1] = tensor_2x_p2_surfx2_eval_quad_node_0_r(rhoux_c); 
  q_cr[2] = tensor_2x_p2_surfx2_eval_quad_node_0_r(rhouy_c); 
  q_cr[3] = tensor_2x_p2_surfx2_eval_quad_node_0_r(rhouz_c); 
  q_cr[4] = tensor_2x_p2_surfx2_eval_quad_node_0_r(energy_c); 
  q_rl[0] = tensor_2x_p2_surfx2_eval_quad_node_0_l(rho_r); 
  q_rl[1] = tensor_2x_p2_surfx2_eval_quad_node_0_l(rhoux_r); 
  q_rl[2] = tensor_2x_p2_surfx2_eval_quad_node_0_l(rhouy_r); 
  q_rl[3] = tensor_2x_p2_surfx2_eval_quad_node_0_l(rhouz_r); 
  q_rl[4] = tensor_2x_p2_surfx2_eval_quad_node_0_l(energy_r); 

  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_l->tau1[1], geom_l->tau2[1], geom_l->norm[1], q_lr, q_lr_local); 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_l->tau1[1], geom_l->tau2[1], geom_l->norm[1], q_cl, q_cl_local); 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_r->tau1[1], geom_r->tau2[1], geom_r->norm[1], q_cr, q_cr_local); 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_r->tau1[1], geom_r->tau2[1], geom_r->norm[1], q_rl, q_rl_local); 

  delta_l[0] = q_cl_local[0] - q_lr_local[0]; 
  delta_l[1] = q_cl_local[1] - q_lr_local[1]; 
  delta_l[2] = q_cl_local[2] - q_lr_local[2]; 
  delta_l[3] = q_cl_local[3] - q_lr_local[3]; 
  delta_l[4] = q_cl_local[4] - q_lr_local[4]; 
  delta_r[0] = q_rl_local[0] - q_cr_local[0]; 
  delta_r[1] = q_rl_local[1] - q_cr_local[1]; 
  delta_r[2] = q_rl_local[2] - q_cr_local[2]; 
  delta_r[3] = q_rl_local[3] - q_cr_local[3]; 
  delta_r[4] = q_rl_local[4] - q_cr_local[4]; 

  my_max_speed_l = gkyl_wv_eqn_waves(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta_l, q_lr_local, q_cl_local, waves_l, speeds_l); 
  my_max_speed_r = gkyl_wv_eqn_waves(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta_r, q_cr_local, q_rl_local, waves_r, speeds_r); 
  lenr_l = geom_l->lenr[1]; 
  speeds_l[0] *= lenr_l; 
  speeds_l[1] *= lenr_l; 
  speeds_l[2] *= lenr_l; 
  lenr_r = geom_r->lenr[1]; 
  speeds_r[0] *= lenr_r; 
  speeds_r[1] *= lenr_r; 
  speeds_r[2] *= lenr_r; 

  gkyl_wv_eqn_qfluct(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, q_lr_local, q_cl_local, waves_l, speeds_l, amdq_l_local, apdq_l_local); 
  gkyl_wv_eqn_qfluct(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, q_cr_local, q_rl_local, waves_r, speeds_r, amdq_r_local, apdq_r_local); 

  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_l->tau1[1], geom_l->tau2[1], geom_l->norm[1], amdq_l_local, amdq_l); 
  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_l->tau1[1], geom_l->tau2[1], geom_l->norm[1], apdq_l_local, apdq_l); 
  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_r->tau1[1], geom_r->tau2[1], geom_r->norm[1], amdq_r_local, amdq_r); 
  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_r->tau1[1], geom_r->tau2[1], geom_r->norm[1], apdq_r_local, apdq_r); 

  amdq_rho_quad_l[0] = amdq_l[0]; 
  apdq_rho_quad_l[0] = apdq_l[0]; 
  amdq_rhoux_quad_l[0] = amdq_l[1]; 
  apdq_rhoux_quad_l[0] = apdq_l[1]; 
  amdq_rhouy_quad_l[0] = amdq_l[2]; 
  apdq_rhouy_quad_l[0] = apdq_l[2]; 
  amdq_rhouz_quad_l[0] = amdq_l[3]; 
  apdq_rhouz_quad_l[0] = apdq_l[3]; 
  amdq_energy_quad_l[0] = amdq_l[4]; 
  apdq_energy_quad_l[0] = apdq_l[4]; 

  amdq_rho_quad_r[0] = amdq_r[0]; 
  apdq_rho_quad_r[0] = apdq_r[0]; 
  amdq_rhoux_quad_r[0] = amdq_r[1]; 
  apdq_rhoux_quad_r[0] = apdq_r[1]; 
  amdq_rhouy_quad_r[0] = amdq_r[2]; 
  apdq_rhouy_quad_r[0] = apdq_r[2]; 
  amdq_rhouz_quad_r[0] = amdq_r[3]; 
  apdq_rhouz_quad_r[0] = apdq_r[3]; 
  amdq_energy_quad_r[0] = amdq_r[4]; 
  apdq_energy_quad_r[0] = apdq_r[4]; 

  q_lr[0] = tensor_2x_p2_surfx2_eval_quad_node_1_r(rho_l); 
  q_lr[1] = tensor_2x_p2_surfx2_eval_quad_node_1_r(rhoux_l); 
  q_lr[2] = tensor_2x_p2_surfx2_eval_quad_node_1_r(rhouy_l); 
  q_lr[3] = tensor_2x_p2_surfx2_eval_quad_node_1_r(rhouz_l); 
  q_lr[4] = tensor_2x_p2_surfx2_eval_quad_node_1_r(energy_l); 
  q_cl[0] = tensor_2x_p2_surfx2_eval_quad_node_1_l(rho_c); 
  q_cl[1] = tensor_2x_p2_surfx2_eval_quad_node_1_l(rhoux_c); 
  q_cl[2] = tensor_2x_p2_surfx2_eval_quad_node_1_l(rhouy_c); 
  q_cl[3] = tensor_2x_p2_surfx2_eval_quad_node_1_l(rhouz_c); 
  q_cl[4] = tensor_2x_p2_surfx2_eval_quad_node_1_l(energy_c); 
  q_cr[0] = tensor_2x_p2_surfx2_eval_quad_node_1_r(rho_c); 
  q_cr[1] = tensor_2x_p2_surfx2_eval_quad_node_1_r(rhoux_c); 
  q_cr[2] = tensor_2x_p2_surfx2_eval_quad_node_1_r(rhouy_c); 
  q_cr[3] = tensor_2x_p2_surfx2_eval_quad_node_1_r(rhouz_c); 
  q_cr[4] = tensor_2x_p2_surfx2_eval_quad_node_1_r(energy_c); 
  q_rl[0] = tensor_2x_p2_surfx2_eval_quad_node_1_l(rho_r); 
  q_rl[1] = tensor_2x_p2_surfx2_eval_quad_node_1_l(rhoux_r); 
  q_rl[2] = tensor_2x_p2_surfx2_eval_quad_node_1_l(rhouy_r); 
  q_rl[3] = tensor_2x_p2_surfx2_eval_quad_node_1_l(rhouz_r); 
  q_rl[4] = tensor_2x_p2_surfx2_eval_quad_node_1_l(energy_r); 

  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_l->tau1[1], geom_l->tau2[1], geom_l->norm[1], q_lr, q_lr_local); 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_l->tau1[1], geom_l->tau2[1], geom_l->norm[1], q_cl, q_cl_local); 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_r->tau1[1], geom_r->tau2[1], geom_r->norm[1], q_cr, q_cr_local); 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_r->tau1[1], geom_r->tau2[1], geom_r->norm[1], q_rl, q_rl_local); 

  delta_l[0] = q_cl_local[0] - q_lr_local[0]; 
  delta_l[1] = q_cl_local[1] - q_lr_local[1]; 
  delta_l[2] = q_cl_local[2] - q_lr_local[2]; 
  delta_l[3] = q_cl_local[3] - q_lr_local[3]; 
  delta_l[4] = q_cl_local[4] - q_lr_local[4]; 
  delta_r[0] = q_rl_local[0] - q_cr_local[0]; 
  delta_r[1] = q_rl_local[1] - q_cr_local[1]; 
  delta_r[2] = q_rl_local[2] - q_cr_local[2]; 
  delta_r[3] = q_rl_local[3] - q_cr_local[3]; 
  delta_r[4] = q_rl_local[4] - q_cr_local[4]; 

  my_max_speed_l = gkyl_wv_eqn_waves(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta_l, q_lr_local, q_cl_local, waves_l, speeds_l); 
  my_max_speed_r = gkyl_wv_eqn_waves(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta_r, q_cr_local, q_rl_local, waves_r, speeds_r); 
  lenr_l = geom_l->lenr[1]; 
  speeds_l[0] *= lenr_l; 
  speeds_l[1] *= lenr_l; 
  speeds_l[2] *= lenr_l; 
  lenr_r = geom_r->lenr[1]; 
  speeds_r[0] *= lenr_r; 
  speeds_r[1] *= lenr_r; 
  speeds_r[2] *= lenr_r; 

  gkyl_wv_eqn_qfluct(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, q_lr_local, q_cl_local, waves_l, speeds_l, amdq_l_local, apdq_l_local); 
  gkyl_wv_eqn_qfluct(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, q_cr_local, q_rl_local, waves_r, speeds_r, amdq_r_local, apdq_r_local); 

  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_l->tau1[1], geom_l->tau2[1], geom_l->norm[1], amdq_l_local, amdq_l); 
  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_l->tau1[1], geom_l->tau2[1], geom_l->norm[1], apdq_l_local, apdq_l); 
  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_r->tau1[1], geom_r->tau2[1], geom_r->norm[1], amdq_r_local, amdq_r); 
  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_r->tau1[1], geom_r->tau2[1], geom_r->norm[1], apdq_r_local, apdq_r); 

  amdq_rho_quad_l[1] = amdq_l[0]; 
  apdq_rho_quad_l[1] = apdq_l[0]; 
  amdq_rhoux_quad_l[1] = amdq_l[1]; 
  apdq_rhoux_quad_l[1] = apdq_l[1]; 
  amdq_rhouy_quad_l[1] = amdq_l[2]; 
  apdq_rhouy_quad_l[1] = apdq_l[2]; 
  amdq_rhouz_quad_l[1] = amdq_l[3]; 
  apdq_rhouz_quad_l[1] = apdq_l[3]; 
  amdq_energy_quad_l[1] = amdq_l[4]; 
  apdq_energy_quad_l[1] = apdq_l[4]; 

  amdq_rho_quad_r[1] = amdq_r[0]; 
  apdq_rho_quad_r[1] = apdq_r[0]; 
  amdq_rhoux_quad_r[1] = amdq_r[1]; 
  apdq_rhoux_quad_r[1] = apdq_r[1]; 
  amdq_rhouy_quad_r[1] = amdq_r[2]; 
  apdq_rhouy_quad_r[1] = apdq_r[2]; 
  amdq_rhouz_quad_r[1] = amdq_r[3]; 
  apdq_rhouz_quad_r[1] = apdq_r[3]; 
  amdq_energy_quad_r[1] = amdq_r[4]; 
  apdq_energy_quad_r[1] = apdq_r[4]; 

  q_lr[0] = tensor_2x_p2_surfx2_eval_quad_node_2_r(rho_l); 
  q_lr[1] = tensor_2x_p2_surfx2_eval_quad_node_2_r(rhoux_l); 
  q_lr[2] = tensor_2x_p2_surfx2_eval_quad_node_2_r(rhouy_l); 
  q_lr[3] = tensor_2x_p2_surfx2_eval_quad_node_2_r(rhouz_l); 
  q_lr[4] = tensor_2x_p2_surfx2_eval_quad_node_2_r(energy_l); 
  q_cl[0] = tensor_2x_p2_surfx2_eval_quad_node_2_l(rho_c); 
  q_cl[1] = tensor_2x_p2_surfx2_eval_quad_node_2_l(rhoux_c); 
  q_cl[2] = tensor_2x_p2_surfx2_eval_quad_node_2_l(rhouy_c); 
  q_cl[3] = tensor_2x_p2_surfx2_eval_quad_node_2_l(rhouz_c); 
  q_cl[4] = tensor_2x_p2_surfx2_eval_quad_node_2_l(energy_c); 
  q_cr[0] = tensor_2x_p2_surfx2_eval_quad_node_2_r(rho_c); 
  q_cr[1] = tensor_2x_p2_surfx2_eval_quad_node_2_r(rhoux_c); 
  q_cr[2] = tensor_2x_p2_surfx2_eval_quad_node_2_r(rhouy_c); 
  q_cr[3] = tensor_2x_p2_surfx2_eval_quad_node_2_r(rhouz_c); 
  q_cr[4] = tensor_2x_p2_surfx2_eval_quad_node_2_r(energy_c); 
  q_rl[0] = tensor_2x_p2_surfx2_eval_quad_node_2_l(rho_r); 
  q_rl[1] = tensor_2x_p2_surfx2_eval_quad_node_2_l(rhoux_r); 
  q_rl[2] = tensor_2x_p2_surfx2_eval_quad_node_2_l(rhouy_r); 
  q_rl[3] = tensor_2x_p2_surfx2_eval_quad_node_2_l(rhouz_r); 
  q_rl[4] = tensor_2x_p2_surfx2_eval_quad_node_2_l(energy_r); 

  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_l->tau1[1], geom_l->tau2[1], geom_l->norm[1], q_lr, q_lr_local); 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_l->tau1[1], geom_l->tau2[1], geom_l->norm[1], q_cl, q_cl_local); 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_r->tau1[1], geom_r->tau2[1], geom_r->norm[1], q_cr, q_cr_local); 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_r->tau1[1], geom_r->tau2[1], geom_r->norm[1], q_rl, q_rl_local); 

  delta_l[0] = q_cl_local[0] - q_lr_local[0]; 
  delta_l[1] = q_cl_local[1] - q_lr_local[1]; 
  delta_l[2] = q_cl_local[2] - q_lr_local[2]; 
  delta_l[3] = q_cl_local[3] - q_lr_local[3]; 
  delta_l[4] = q_cl_local[4] - q_lr_local[4]; 
  delta_r[0] = q_rl_local[0] - q_cr_local[0]; 
  delta_r[1] = q_rl_local[1] - q_cr_local[1]; 
  delta_r[2] = q_rl_local[2] - q_cr_local[2]; 
  delta_r[3] = q_rl_local[3] - q_cr_local[3]; 
  delta_r[4] = q_rl_local[4] - q_cr_local[4]; 

  my_max_speed_l = gkyl_wv_eqn_waves(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta_l, q_lr_local, q_cl_local, waves_l, speeds_l); 
  my_max_speed_r = gkyl_wv_eqn_waves(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta_r, q_cr_local, q_rl_local, waves_r, speeds_r); 
  lenr_l = geom_l->lenr[1]; 
  speeds_l[0] *= lenr_l; 
  speeds_l[1] *= lenr_l; 
  speeds_l[2] *= lenr_l; 
  lenr_r = geom_r->lenr[1]; 
  speeds_r[0] *= lenr_r; 
  speeds_r[1] *= lenr_r; 
  speeds_r[2] *= lenr_r; 

  gkyl_wv_eqn_qfluct(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, q_lr_local, q_cl_local, waves_l, speeds_l, amdq_l_local, apdq_l_local); 
  gkyl_wv_eqn_qfluct(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, q_cr_local, q_rl_local, waves_r, speeds_r, amdq_r_local, apdq_r_local); 

  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_l->tau1[1], geom_l->tau2[1], geom_l->norm[1], amdq_l_local, amdq_l); 
  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_l->tau1[1], geom_l->tau2[1], geom_l->norm[1], apdq_l_local, apdq_l); 
  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_r->tau1[1], geom_r->tau2[1], geom_r->norm[1], amdq_r_local, amdq_r); 
  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_r->tau1[1], geom_r->tau2[1], geom_r->norm[1], apdq_r_local, apdq_r); 

  amdq_rho_quad_l[2] = amdq_l[0]; 
  apdq_rho_quad_l[2] = apdq_l[0]; 
  amdq_rhoux_quad_l[2] = amdq_l[1]; 
  apdq_rhoux_quad_l[2] = apdq_l[1]; 
  amdq_rhouy_quad_l[2] = amdq_l[2]; 
  apdq_rhouy_quad_l[2] = apdq_l[2]; 
  amdq_rhouz_quad_l[2] = amdq_l[3]; 
  apdq_rhouz_quad_l[2] = apdq_l[3]; 
  amdq_energy_quad_l[2] = amdq_l[4]; 
  apdq_energy_quad_l[2] = apdq_l[4]; 

  amdq_rho_quad_r[2] = amdq_r[0]; 
  apdq_rho_quad_r[2] = apdq_r[0]; 
  amdq_rhoux_quad_r[2] = amdq_r[1]; 
  apdq_rhoux_quad_r[2] = apdq_r[1]; 
  amdq_rhouy_quad_r[2] = amdq_r[2]; 
  apdq_rhouy_quad_r[2] = apdq_r[2]; 
  amdq_rhouz_quad_r[2] = amdq_r[3]; 
  apdq_rhouz_quad_r[2] = apdq_r[3]; 
  amdq_energy_quad_r[2] = amdq_r[4]; 
  apdq_energy_quad_r[2] = apdq_r[4]; 

  tensor_2x_p2_upwind_quad_to_modal(amdq_rho_quad_l, amdq_rho_l); 
  tensor_2x_p2_upwind_quad_to_modal(amdq_rhoux_quad_l, amdq_rhoux_l); 
  tensor_2x_p2_upwind_quad_to_modal(amdq_rhouy_quad_l, amdq_rhouy_l); 
  tensor_2x_p2_upwind_quad_to_modal(amdq_rhouz_quad_l, amdq_rhouz_l); 
  tensor_2x_p2_upwind_quad_to_modal(amdq_energy_quad_l, amdq_energy_l); 

  tensor_2x_p2_upwind_quad_to_modal(apdq_rho_quad_l, apdq_rho_l); 
  tensor_2x_p2_upwind_quad_to_modal(apdq_rhoux_quad_l, apdq_rhoux_l); 
  tensor_2x_p2_upwind_quad_to_modal(apdq_rhouy_quad_l, apdq_rhouy_l); 
  tensor_2x_p2_upwind_quad_to_modal(apdq_rhouz_quad_l, apdq_rhouz_l); 
  tensor_2x_p2_upwind_quad_to_modal(apdq_energy_quad_l, apdq_energy_l); 

  tensor_2x_p2_upwind_quad_to_modal(amdq_rho_quad_r, amdq_rho_r); 
  tensor_2x_p2_upwind_quad_to_modal(amdq_rhoux_quad_r, amdq_rhoux_r); 
  tensor_2x_p2_upwind_quad_to_modal(amdq_rhouy_quad_r, amdq_rhouy_r); 
  tensor_2x_p2_upwind_quad_to_modal(amdq_rhouz_quad_r, amdq_rhouz_r); 
  tensor_2x_p2_upwind_quad_to_modal(amdq_energy_quad_r, amdq_energy_r); 

  tensor_2x_p2_upwind_quad_to_modal(apdq_rho_quad_r, apdq_rho_r); 
  tensor_2x_p2_upwind_quad_to_modal(apdq_rhoux_quad_r, apdq_rhoux_r); 
  tensor_2x_p2_upwind_quad_to_modal(apdq_rhouy_quad_r, apdq_rhouy_r); 
  tensor_2x_p2_upwind_quad_to_modal(apdq_rhouz_quad_r, apdq_rhouz_r); 
  tensor_2x_p2_upwind_quad_to_modal(apdq_energy_quad_r, apdq_energy_r); 

  double flux_rho_l[3] = {0.0}; 
  double flux_rho_r[3] = {0.0}; 

  flux_rho_l[0] = 0.2795084971874737*uy_surf_lr[2]*rho_l[8]+0.2795084971874737*uy_surf_cl[2]*rho_l[8]+0.2795084971874737*uy_surf_lr[2]*rho_c[8]+0.2795084971874737*uy_surf_cl[2]*rho_c[8]+0.2795084971874738*uy_surf_lr[1]*rho_l[7]+0.2795084971874738*uy_surf_cl[1]*rho_l[7]+0.2795084971874738*uy_surf_lr[1]*rho_c[7]+0.2795084971874738*uy_surf_cl[1]*rho_c[7]+0.2165063509461097*uy_surf_lr[2]*rho_l[6]+0.2165063509461097*uy_surf_cl[2]*rho_l[6]-0.2165063509461097*uy_surf_lr[2]*rho_c[6]-0.2165063509461097*uy_surf_cl[2]*rho_c[6]+0.2795084971874737*uy_surf_lr[0]*rho_l[5]+0.2795084971874737*uy_surf_cl[0]*rho_l[5]+0.2795084971874737*uy_surf_lr[0]*rho_c[5]+0.2795084971874737*uy_surf_cl[0]*rho_c[5]+0.125*uy_surf_lr[2]*rho_l[4]+0.125*uy_surf_cl[2]*rho_l[4]+0.125*uy_surf_lr[2]*rho_c[4]+0.125*uy_surf_cl[2]*rho_c[4]+0.2165063509461096*uy_surf_lr[1]*rho_l[3]+0.2165063509461096*uy_surf_cl[1]*rho_l[3]-0.2165063509461096*uy_surf_lr[1]*rho_c[3]-0.2165063509461096*uy_surf_cl[1]*rho_c[3]+0.2165063509461096*uy_surf_lr[0]*rho_l[2]+0.2165063509461096*uy_surf_cl[0]*rho_l[2]-0.2165063509461096*uy_surf_lr[0]*rho_c[2]-0.2165063509461096*uy_surf_cl[0]*rho_c[2]+0.125*rho_l[1]*uy_surf_lr[1]+0.125*rho_c[1]*uy_surf_lr[1]+0.125*rho_l[1]*uy_surf_cl[1]+0.125*rho_c[1]*uy_surf_cl[1]+0.125*rho_l[0]*uy_surf_lr[0]+0.125*rho_c[0]*uy_surf_lr[0]+0.125*rho_l[0]*uy_surf_cl[0]+0.125*rho_c[0]*uy_surf_cl[0]-0.5*apdq_rho_l[0]+0.5*amdq_rho_l[0]; 
  flux_rho_l[1] = 0.25*uy_surf_lr[1]*rho_l[8]+0.25*uy_surf_cl[1]*rho_l[8]+0.25*uy_surf_lr[1]*rho_c[8]+0.25*uy_surf_cl[1]*rho_c[8]+0.2500000000000001*uy_surf_lr[2]*rho_l[7]+0.2500000000000001*uy_surf_cl[2]*rho_l[7]+0.2795084971874738*uy_surf_lr[0]*rho_l[7]+0.2795084971874738*uy_surf_cl[0]*rho_l[7]+0.2500000000000001*uy_surf_lr[2]*rho_c[7]+0.2500000000000001*uy_surf_cl[2]*rho_c[7]+0.2795084971874738*uy_surf_lr[0]*rho_c[7]+0.2795084971874738*uy_surf_cl[0]*rho_c[7]+0.1936491673103709*uy_surf_lr[1]*rho_l[6]+0.1936491673103709*uy_surf_cl[1]*rho_l[6]-0.1936491673103709*uy_surf_lr[1]*rho_c[6]-0.1936491673103709*uy_surf_cl[1]*rho_c[6]+0.2795084971874737*uy_surf_lr[1]*rho_l[5]+0.2795084971874737*uy_surf_cl[1]*rho_l[5]+0.2795084971874737*uy_surf_lr[1]*rho_c[5]+0.2795084971874737*uy_surf_cl[1]*rho_c[5]+0.1118033988749895*uy_surf_lr[1]*rho_l[4]+0.1118033988749895*uy_surf_cl[1]*rho_l[4]+0.1118033988749895*uy_surf_lr[1]*rho_c[4]+0.1118033988749895*uy_surf_cl[1]*rho_c[4]+0.1936491673103708*uy_surf_lr[2]*rho_l[3]+0.1936491673103708*uy_surf_cl[2]*rho_l[3]+0.2165063509461096*uy_surf_lr[0]*rho_l[3]+0.2165063509461096*uy_surf_cl[0]*rho_l[3]-0.1936491673103708*uy_surf_lr[2]*rho_c[3]-0.1936491673103708*uy_surf_cl[2]*rho_c[3]-0.2165063509461096*uy_surf_lr[0]*rho_c[3]-0.2165063509461096*uy_surf_cl[0]*rho_c[3]+0.1118033988749895*rho_l[1]*uy_surf_lr[2]+0.1118033988749895*rho_c[1]*uy_surf_lr[2]+0.1118033988749895*rho_l[1]*uy_surf_cl[2]+0.1118033988749895*rho_c[1]*uy_surf_cl[2]+0.2165063509461096*uy_surf_lr[1]*rho_l[2]+0.2165063509461096*uy_surf_cl[1]*rho_l[2]-0.2165063509461096*uy_surf_lr[1]*rho_c[2]-0.2165063509461096*uy_surf_cl[1]*rho_c[2]+0.125*rho_l[0]*uy_surf_lr[1]+0.125*rho_c[0]*uy_surf_lr[1]+0.125*rho_l[0]*uy_surf_cl[1]+0.125*rho_c[0]*uy_surf_cl[1]+0.125*uy_surf_lr[0]*rho_l[1]+0.125*uy_surf_cl[0]*rho_l[1]+0.125*uy_surf_lr[0]*rho_c[1]+0.125*uy_surf_cl[0]*rho_c[1]-0.5*apdq_rho_l[1]+0.5*amdq_rho_l[1]; 
  flux_rho_l[2] = 0.1785714285714286*uy_surf_lr[2]*rho_l[8]+0.1785714285714286*uy_surf_cl[2]*rho_l[8]+0.2795084971874737*uy_surf_lr[0]*rho_l[8]+0.2795084971874737*uy_surf_cl[0]*rho_l[8]+0.1785714285714286*uy_surf_lr[2]*rho_c[8]+0.1785714285714286*uy_surf_cl[2]*rho_c[8]+0.2795084971874737*uy_surf_lr[0]*rho_c[8]+0.2795084971874737*uy_surf_cl[0]*rho_c[8]+0.2500000000000001*uy_surf_lr[1]*rho_l[7]+0.2500000000000001*uy_surf_cl[1]*rho_l[7]+0.2500000000000001*uy_surf_lr[1]*rho_c[7]+0.2500000000000001*uy_surf_cl[1]*rho_c[7]+0.138320833793122*uy_surf_lr[2]*rho_l[6]+0.138320833793122*uy_surf_cl[2]*rho_l[6]+0.2165063509461097*uy_surf_lr[0]*rho_l[6]+0.2165063509461097*uy_surf_cl[0]*rho_l[6]-0.138320833793122*uy_surf_lr[2]*rho_c[6]-0.138320833793122*uy_surf_cl[2]*rho_c[6]-0.2165063509461097*uy_surf_lr[0]*rho_c[6]-0.2165063509461097*uy_surf_cl[0]*rho_c[6]+0.2795084971874737*uy_surf_lr[2]*rho_l[5]+0.2795084971874737*uy_surf_cl[2]*rho_l[5]+0.2795084971874737*uy_surf_lr[2]*rho_c[5]+0.2795084971874737*uy_surf_cl[2]*rho_c[5]+0.07985957062499249*uy_surf_lr[2]*rho_l[4]+0.07985957062499249*uy_surf_cl[2]*rho_l[4]+0.125*uy_surf_lr[0]*rho_l[4]+0.125*uy_surf_cl[0]*rho_l[4]+0.07985957062499249*uy_surf_lr[2]*rho_c[4]+0.07985957062499249*uy_surf_cl[2]*rho_c[4]+0.125*uy_surf_lr[0]*rho_c[4]+0.125*uy_surf_cl[0]*rho_c[4]+0.1936491673103708*uy_surf_lr[1]*rho_l[3]+0.1936491673103708*uy_surf_cl[1]*rho_l[3]-0.1936491673103708*uy_surf_lr[1]*rho_c[3]-0.1936491673103708*uy_surf_cl[1]*rho_c[3]+0.2165063509461096*rho_l[2]*uy_surf_lr[2]-0.2165063509461096*rho_c[2]*uy_surf_lr[2]+0.125*rho_l[0]*uy_surf_lr[2]+0.125*rho_c[0]*uy_surf_lr[2]+0.2165063509461096*rho_l[2]*uy_surf_cl[2]-0.2165063509461096*rho_c[2]*uy_surf_cl[2]+0.125*rho_l[0]*uy_surf_cl[2]+0.125*rho_c[0]*uy_surf_cl[2]-0.5*apdq_rho_l[2]+0.5*amdq_rho_l[2]+0.1118033988749895*rho_l[1]*uy_surf_lr[1]+0.1118033988749895*rho_c[1]*uy_surf_lr[1]+0.1118033988749895*rho_l[1]*uy_surf_cl[1]+0.1118033988749895*rho_c[1]*uy_surf_cl[1]; 

  flux_rho_r[0] = 0.2795084971874737*uy_surf_rl[2]*rho_r[8]+0.2795084971874737*uy_surf_cr[2]*rho_r[8]+0.2795084971874737*uy_surf_rl[2]*rho_c[8]+0.2795084971874737*uy_surf_cr[2]*rho_c[8]+0.2795084971874738*uy_surf_rl[1]*rho_r[7]+0.2795084971874738*uy_surf_cr[1]*rho_r[7]+0.2795084971874738*uy_surf_rl[1]*rho_c[7]+0.2795084971874738*uy_surf_cr[1]*rho_c[7]-0.2165063509461097*uy_surf_rl[2]*rho_r[6]-0.2165063509461097*uy_surf_cr[2]*rho_r[6]+0.2165063509461097*uy_surf_rl[2]*rho_c[6]+0.2165063509461097*uy_surf_cr[2]*rho_c[6]+0.2795084971874737*uy_surf_rl[0]*rho_r[5]+0.2795084971874737*uy_surf_cr[0]*rho_r[5]+0.2795084971874737*uy_surf_rl[0]*rho_c[5]+0.2795084971874737*uy_surf_cr[0]*rho_c[5]+0.125*uy_surf_rl[2]*rho_r[4]+0.125*uy_surf_cr[2]*rho_r[4]+0.125*uy_surf_rl[2]*rho_c[4]+0.125*uy_surf_cr[2]*rho_c[4]-0.2165063509461096*uy_surf_rl[1]*rho_r[3]-0.2165063509461096*uy_surf_cr[1]*rho_r[3]+0.2165063509461096*uy_surf_rl[1]*rho_c[3]+0.2165063509461096*uy_surf_cr[1]*rho_c[3]-0.2165063509461096*uy_surf_rl[0]*rho_r[2]-0.2165063509461096*uy_surf_cr[0]*rho_r[2]+0.2165063509461096*uy_surf_rl[0]*rho_c[2]+0.2165063509461096*uy_surf_cr[0]*rho_c[2]+0.125*rho_r[1]*uy_surf_rl[1]+0.125*rho_c[1]*uy_surf_rl[1]+0.125*rho_r[1]*uy_surf_cr[1]+0.125*rho_c[1]*uy_surf_cr[1]+0.125*rho_r[0]*uy_surf_rl[0]+0.125*rho_c[0]*uy_surf_rl[0]+0.125*rho_r[0]*uy_surf_cr[0]+0.125*rho_c[0]*uy_surf_cr[0]-0.5*apdq_rho_r[0]+0.5*amdq_rho_r[0]; 
  flux_rho_r[1] = 0.25*uy_surf_rl[1]*rho_r[8]+0.25*uy_surf_cr[1]*rho_r[8]+0.25*uy_surf_rl[1]*rho_c[8]+0.25*uy_surf_cr[1]*rho_c[8]+0.2500000000000001*uy_surf_rl[2]*rho_r[7]+0.2500000000000001*uy_surf_cr[2]*rho_r[7]+0.2795084971874738*uy_surf_rl[0]*rho_r[7]+0.2795084971874738*uy_surf_cr[0]*rho_r[7]+0.2500000000000001*uy_surf_rl[2]*rho_c[7]+0.2500000000000001*uy_surf_cr[2]*rho_c[7]+0.2795084971874738*uy_surf_rl[0]*rho_c[7]+0.2795084971874738*uy_surf_cr[0]*rho_c[7]-0.1936491673103709*uy_surf_rl[1]*rho_r[6]-0.1936491673103709*uy_surf_cr[1]*rho_r[6]+0.1936491673103709*uy_surf_rl[1]*rho_c[6]+0.1936491673103709*uy_surf_cr[1]*rho_c[6]+0.2795084971874737*uy_surf_rl[1]*rho_r[5]+0.2795084971874737*uy_surf_cr[1]*rho_r[5]+0.2795084971874737*uy_surf_rl[1]*rho_c[5]+0.2795084971874737*uy_surf_cr[1]*rho_c[5]+0.1118033988749895*uy_surf_rl[1]*rho_r[4]+0.1118033988749895*uy_surf_cr[1]*rho_r[4]+0.1118033988749895*uy_surf_rl[1]*rho_c[4]+0.1118033988749895*uy_surf_cr[1]*rho_c[4]-0.1936491673103708*uy_surf_rl[2]*rho_r[3]-0.1936491673103708*uy_surf_cr[2]*rho_r[3]-0.2165063509461096*uy_surf_rl[0]*rho_r[3]-0.2165063509461096*uy_surf_cr[0]*rho_r[3]+0.1936491673103708*uy_surf_rl[2]*rho_c[3]+0.1936491673103708*uy_surf_cr[2]*rho_c[3]+0.2165063509461096*uy_surf_rl[0]*rho_c[3]+0.2165063509461096*uy_surf_cr[0]*rho_c[3]+0.1118033988749895*rho_r[1]*uy_surf_rl[2]+0.1118033988749895*rho_c[1]*uy_surf_rl[2]+0.1118033988749895*rho_r[1]*uy_surf_cr[2]+0.1118033988749895*rho_c[1]*uy_surf_cr[2]-0.2165063509461096*uy_surf_rl[1]*rho_r[2]-0.2165063509461096*uy_surf_cr[1]*rho_r[2]+0.2165063509461096*uy_surf_rl[1]*rho_c[2]+0.2165063509461096*uy_surf_cr[1]*rho_c[2]+0.125*rho_r[0]*uy_surf_rl[1]+0.125*rho_c[0]*uy_surf_rl[1]+0.125*rho_r[0]*uy_surf_cr[1]+0.125*rho_c[0]*uy_surf_cr[1]+0.125*uy_surf_rl[0]*rho_r[1]+0.125*uy_surf_cr[0]*rho_r[1]+0.125*uy_surf_rl[0]*rho_c[1]+0.125*uy_surf_cr[0]*rho_c[1]-0.5*apdq_rho_r[1]+0.5*amdq_rho_r[1]; 
  flux_rho_r[2] = 0.1785714285714286*uy_surf_rl[2]*rho_r[8]+0.1785714285714286*uy_surf_cr[2]*rho_r[8]+0.2795084971874737*uy_surf_rl[0]*rho_r[8]+0.2795084971874737*uy_surf_cr[0]*rho_r[8]+0.1785714285714286*uy_surf_rl[2]*rho_c[8]+0.1785714285714286*uy_surf_cr[2]*rho_c[8]+0.2795084971874737*uy_surf_rl[0]*rho_c[8]+0.2795084971874737*uy_surf_cr[0]*rho_c[8]+0.2500000000000001*uy_surf_rl[1]*rho_r[7]+0.2500000000000001*uy_surf_cr[1]*rho_r[7]+0.2500000000000001*uy_surf_rl[1]*rho_c[7]+0.2500000000000001*uy_surf_cr[1]*rho_c[7]-0.138320833793122*uy_surf_rl[2]*rho_r[6]-0.138320833793122*uy_surf_cr[2]*rho_r[6]-0.2165063509461097*uy_surf_rl[0]*rho_r[6]-0.2165063509461097*uy_surf_cr[0]*rho_r[6]+0.138320833793122*uy_surf_rl[2]*rho_c[6]+0.138320833793122*uy_surf_cr[2]*rho_c[6]+0.2165063509461097*uy_surf_rl[0]*rho_c[6]+0.2165063509461097*uy_surf_cr[0]*rho_c[6]+0.2795084971874737*uy_surf_rl[2]*rho_r[5]+0.2795084971874737*uy_surf_cr[2]*rho_r[5]+0.2795084971874737*uy_surf_rl[2]*rho_c[5]+0.2795084971874737*uy_surf_cr[2]*rho_c[5]+0.07985957062499249*uy_surf_rl[2]*rho_r[4]+0.07985957062499249*uy_surf_cr[2]*rho_r[4]+0.125*uy_surf_rl[0]*rho_r[4]+0.125*uy_surf_cr[0]*rho_r[4]+0.07985957062499249*uy_surf_rl[2]*rho_c[4]+0.07985957062499249*uy_surf_cr[2]*rho_c[4]+0.125*uy_surf_rl[0]*rho_c[4]+0.125*uy_surf_cr[0]*rho_c[4]-0.1936491673103708*uy_surf_rl[1]*rho_r[3]-0.1936491673103708*uy_surf_cr[1]*rho_r[3]+0.1936491673103708*uy_surf_rl[1]*rho_c[3]+0.1936491673103708*uy_surf_cr[1]*rho_c[3]-0.2165063509461096*rho_r[2]*uy_surf_rl[2]+0.2165063509461096*rho_c[2]*uy_surf_rl[2]+0.125*rho_r[0]*uy_surf_rl[2]+0.125*rho_c[0]*uy_surf_rl[2]-0.2165063509461096*rho_r[2]*uy_surf_cr[2]+0.2165063509461096*rho_c[2]*uy_surf_cr[2]+0.125*rho_r[0]*uy_surf_cr[2]+0.125*rho_c[0]*uy_surf_cr[2]-0.5*apdq_rho_r[2]+0.5*amdq_rho_r[2]+0.1118033988749895*rho_r[1]*uy_surf_rl[1]+0.1118033988749895*rho_c[1]*uy_surf_rl[1]+0.1118033988749895*rho_r[1]*uy_surf_cr[1]+0.1118033988749895*rho_c[1]*uy_surf_cr[1]; 

  double flux_rhoux_l[3] = {0.0}; 
  double flux_rhoux_r[3] = {0.0}; 

  flux_rhoux_l[0] = 0.3535533905932737*flux_rho_l[2]*ux_surf_lr[2]+0.3535533905932737*flux_rho_l[2]*ux_surf_cl[2]+0.3535533905932737*flux_rho_l[1]*ux_surf_lr[1]+0.3535533905932737*flux_rho_l[1]*ux_surf_cl[1]+0.3535533905932737*flux_rho_l[0]*ux_surf_lr[0]+0.3535533905932737*flux_rho_l[0]*ux_surf_cl[0]-0.5*apdq_rhoux_l[0]+0.5*amdq_rhoux_l[0]; 
  flux_rhoux_l[1] = 0.3162277660168379*flux_rho_l[1]*ux_surf_lr[2]+0.3162277660168379*flux_rho_l[1]*ux_surf_cl[2]+0.3162277660168379*ux_surf_lr[1]*flux_rho_l[2]+0.3162277660168379*ux_surf_cl[1]*flux_rho_l[2]+0.3535533905932737*flux_rho_l[0]*ux_surf_lr[1]+0.3535533905932737*flux_rho_l[0]*ux_surf_cl[1]+0.3535533905932737*ux_surf_lr[0]*flux_rho_l[1]+0.3535533905932737*ux_surf_cl[0]*flux_rho_l[1]-0.5*apdq_rhoux_l[1]+0.5*amdq_rhoux_l[1]; 
  flux_rhoux_l[2] = 0.2258769757263128*flux_rho_l[2]*ux_surf_lr[2]+0.3535533905932737*flux_rho_l[0]*ux_surf_lr[2]+0.2258769757263128*flux_rho_l[2]*ux_surf_cl[2]+0.3535533905932737*flux_rho_l[0]*ux_surf_cl[2]+0.3535533905932737*ux_surf_lr[0]*flux_rho_l[2]+0.3535533905932737*ux_surf_cl[0]*flux_rho_l[2]-0.5*apdq_rhoux_l[2]+0.5*amdq_rhoux_l[2]+0.3162277660168379*flux_rho_l[1]*ux_surf_lr[1]+0.3162277660168379*flux_rho_l[1]*ux_surf_cl[1]; 

  flux_rhoux_r[0] = 0.3535533905932737*flux_rho_r[2]*ux_surf_rl[2]+0.3535533905932737*flux_rho_r[2]*ux_surf_cr[2]+0.3535533905932737*flux_rho_r[1]*ux_surf_rl[1]+0.3535533905932737*flux_rho_r[1]*ux_surf_cr[1]+0.3535533905932737*flux_rho_r[0]*ux_surf_rl[0]+0.3535533905932737*flux_rho_r[0]*ux_surf_cr[0]-0.5*apdq_rhoux_r[0]+0.5*amdq_rhoux_r[0]; 
  flux_rhoux_r[1] = 0.3162277660168379*flux_rho_r[1]*ux_surf_rl[2]+0.3162277660168379*flux_rho_r[1]*ux_surf_cr[2]+0.3162277660168379*ux_surf_rl[1]*flux_rho_r[2]+0.3162277660168379*ux_surf_cr[1]*flux_rho_r[2]+0.3535533905932737*flux_rho_r[0]*ux_surf_rl[1]+0.3535533905932737*flux_rho_r[0]*ux_surf_cr[1]+0.3535533905932737*ux_surf_rl[0]*flux_rho_r[1]+0.3535533905932737*ux_surf_cr[0]*flux_rho_r[1]-0.5*apdq_rhoux_r[1]+0.5*amdq_rhoux_r[1]; 
  flux_rhoux_r[2] = 0.2258769757263128*flux_rho_r[2]*ux_surf_rl[2]+0.3535533905932737*flux_rho_r[0]*ux_surf_rl[2]+0.2258769757263128*flux_rho_r[2]*ux_surf_cr[2]+0.3535533905932737*flux_rho_r[0]*ux_surf_cr[2]+0.3535533905932737*ux_surf_rl[0]*flux_rho_r[2]+0.3535533905932737*ux_surf_cr[0]*flux_rho_r[2]-0.5*apdq_rhoux_r[2]+0.5*amdq_rhoux_r[2]+0.3162277660168379*flux_rho_r[1]*ux_surf_rl[1]+0.3162277660168379*flux_rho_r[1]*ux_surf_cr[1]; 

  double flux_rhouy_l[3] = {0.0}; 
  double flux_rhouy_r[3] = {0.0}; 

  flux_rhouy_l[0] = 0.3535533905932737*flux_rho_l[2]*uy_surf_lr[2]+0.3535533905932737*flux_rho_l[2]*uy_surf_cl[2]+0.3535533905932737*flux_rho_l[1]*uy_surf_lr[1]+0.3535533905932737*flux_rho_l[1]*uy_surf_cl[1]+0.3535533905932737*flux_rho_l[0]*uy_surf_lr[0]+0.3535533905932737*flux_rho_l[0]*uy_surf_cl[0]+0.5*p_surf_lr[0]+0.5*p_surf_cl[0]-0.5*apdq_rhouy_l[0]+0.5*amdq_rhouy_l[0]; 
  flux_rhouy_l[1] = 0.3162277660168379*flux_rho_l[1]*uy_surf_lr[2]+0.3162277660168379*flux_rho_l[1]*uy_surf_cl[2]+0.3162277660168379*uy_surf_lr[1]*flux_rho_l[2]+0.3162277660168379*uy_surf_cl[1]*flux_rho_l[2]+0.3535533905932737*flux_rho_l[0]*uy_surf_lr[1]+0.3535533905932737*flux_rho_l[0]*uy_surf_cl[1]+0.5*p_surf_lr[1]+0.5*p_surf_cl[1]+0.3535533905932737*uy_surf_lr[0]*flux_rho_l[1]+0.3535533905932737*uy_surf_cl[0]*flux_rho_l[1]-0.5*apdq_rhouy_l[1]+0.5*amdq_rhouy_l[1]; 
  flux_rhouy_l[2] = 0.2258769757263128*flux_rho_l[2]*uy_surf_lr[2]+0.3535533905932737*flux_rho_l[0]*uy_surf_lr[2]+0.2258769757263128*flux_rho_l[2]*uy_surf_cl[2]+0.3535533905932737*flux_rho_l[0]*uy_surf_cl[2]+0.5*p_surf_lr[2]+0.5*p_surf_cl[2]+0.3535533905932737*uy_surf_lr[0]*flux_rho_l[2]+0.3535533905932737*uy_surf_cl[0]*flux_rho_l[2]-0.5*apdq_rhouy_l[2]+0.5*amdq_rhouy_l[2]+0.3162277660168379*flux_rho_l[1]*uy_surf_lr[1]+0.3162277660168379*flux_rho_l[1]*uy_surf_cl[1]; 

  flux_rhouy_r[0] = 0.3535533905932737*flux_rho_r[2]*uy_surf_rl[2]+0.3535533905932737*flux_rho_r[2]*uy_surf_cr[2]+0.3535533905932737*flux_rho_r[1]*uy_surf_rl[1]+0.3535533905932737*flux_rho_r[1]*uy_surf_cr[1]+0.3535533905932737*flux_rho_r[0]*uy_surf_rl[0]+0.3535533905932737*flux_rho_r[0]*uy_surf_cr[0]+0.5*p_surf_rl[0]+0.5*p_surf_cr[0]-0.5*apdq_rhouy_r[0]+0.5*amdq_rhouy_r[0]; 
  flux_rhouy_r[1] = 0.3162277660168379*flux_rho_r[1]*uy_surf_rl[2]+0.3162277660168379*flux_rho_r[1]*uy_surf_cr[2]+0.3162277660168379*uy_surf_rl[1]*flux_rho_r[2]+0.3162277660168379*uy_surf_cr[1]*flux_rho_r[2]+0.3535533905932737*flux_rho_r[0]*uy_surf_rl[1]+0.3535533905932737*flux_rho_r[0]*uy_surf_cr[1]+0.5*p_surf_rl[1]+0.5*p_surf_cr[1]+0.3535533905932737*uy_surf_rl[0]*flux_rho_r[1]+0.3535533905932737*uy_surf_cr[0]*flux_rho_r[1]-0.5*apdq_rhouy_r[1]+0.5*amdq_rhouy_r[1]; 
  flux_rhouy_r[2] = 0.2258769757263128*flux_rho_r[2]*uy_surf_rl[2]+0.3535533905932737*flux_rho_r[0]*uy_surf_rl[2]+0.2258769757263128*flux_rho_r[2]*uy_surf_cr[2]+0.3535533905932737*flux_rho_r[0]*uy_surf_cr[2]+0.5*p_surf_rl[2]+0.5*p_surf_cr[2]+0.3535533905932737*uy_surf_rl[0]*flux_rho_r[2]+0.3535533905932737*uy_surf_cr[0]*flux_rho_r[2]-0.5*apdq_rhouy_r[2]+0.5*amdq_rhouy_r[2]+0.3162277660168379*flux_rho_r[1]*uy_surf_rl[1]+0.3162277660168379*flux_rho_r[1]*uy_surf_cr[1]; 

  double flux_rhouz_l[3] = {0.0}; 
  double flux_rhouz_r[3] = {0.0}; 

  flux_rhouz_l[0] = 0.3535533905932737*flux_rho_l[2]*uz_surf_lr[2]+0.3535533905932737*flux_rho_l[2]*uz_surf_cl[2]+0.3535533905932737*flux_rho_l[1]*uz_surf_lr[1]+0.3535533905932737*flux_rho_l[1]*uz_surf_cl[1]+0.3535533905932737*flux_rho_l[0]*uz_surf_lr[0]+0.3535533905932737*flux_rho_l[0]*uz_surf_cl[0]-0.5*apdq_rhouz_l[0]+0.5*amdq_rhouz_l[0]; 
  flux_rhouz_l[1] = 0.3162277660168379*flux_rho_l[1]*uz_surf_lr[2]+0.3162277660168379*flux_rho_l[1]*uz_surf_cl[2]+0.3162277660168379*uz_surf_lr[1]*flux_rho_l[2]+0.3162277660168379*uz_surf_cl[1]*flux_rho_l[2]+0.3535533905932737*flux_rho_l[0]*uz_surf_lr[1]+0.3535533905932737*flux_rho_l[0]*uz_surf_cl[1]+0.3535533905932737*uz_surf_lr[0]*flux_rho_l[1]+0.3535533905932737*uz_surf_cl[0]*flux_rho_l[1]-0.5*apdq_rhouz_l[1]+0.5*amdq_rhouz_l[1]; 
  flux_rhouz_l[2] = 0.2258769757263128*flux_rho_l[2]*uz_surf_lr[2]+0.3535533905932737*flux_rho_l[0]*uz_surf_lr[2]+0.2258769757263128*flux_rho_l[2]*uz_surf_cl[2]+0.3535533905932737*flux_rho_l[0]*uz_surf_cl[2]+0.3535533905932737*uz_surf_lr[0]*flux_rho_l[2]+0.3535533905932737*uz_surf_cl[0]*flux_rho_l[2]-0.5*apdq_rhouz_l[2]+0.5*amdq_rhouz_l[2]+0.3162277660168379*flux_rho_l[1]*uz_surf_lr[1]+0.3162277660168379*flux_rho_l[1]*uz_surf_cl[1]; 

  flux_rhouz_r[0] = 0.3535533905932737*flux_rho_r[2]*uz_surf_rl[2]+0.3535533905932737*flux_rho_r[2]*uz_surf_cr[2]+0.3535533905932737*flux_rho_r[1]*uz_surf_rl[1]+0.3535533905932737*flux_rho_r[1]*uz_surf_cr[1]+0.3535533905932737*flux_rho_r[0]*uz_surf_rl[0]+0.3535533905932737*flux_rho_r[0]*uz_surf_cr[0]-0.5*apdq_rhouz_r[0]+0.5*amdq_rhouz_r[0]; 
  flux_rhouz_r[1] = 0.3162277660168379*flux_rho_r[1]*uz_surf_rl[2]+0.3162277660168379*flux_rho_r[1]*uz_surf_cr[2]+0.3162277660168379*uz_surf_rl[1]*flux_rho_r[2]+0.3162277660168379*uz_surf_cr[1]*flux_rho_r[2]+0.3535533905932737*flux_rho_r[0]*uz_surf_rl[1]+0.3535533905932737*flux_rho_r[0]*uz_surf_cr[1]+0.3535533905932737*uz_surf_rl[0]*flux_rho_r[1]+0.3535533905932737*uz_surf_cr[0]*flux_rho_r[1]-0.5*apdq_rhouz_r[1]+0.5*amdq_rhouz_r[1]; 
  flux_rhouz_r[2] = 0.2258769757263128*flux_rho_r[2]*uz_surf_rl[2]+0.3535533905932737*flux_rho_r[0]*uz_surf_rl[2]+0.2258769757263128*flux_rho_r[2]*uz_surf_cr[2]+0.3535533905932737*flux_rho_r[0]*uz_surf_cr[2]+0.3535533905932737*uz_surf_rl[0]*flux_rho_r[2]+0.3535533905932737*uz_surf_cr[0]*flux_rho_r[2]-0.5*apdq_rhouz_r[2]+0.5*amdq_rhouz_r[2]+0.3162277660168379*flux_rho_r[1]*uz_surf_rl[1]+0.3162277660168379*flux_rho_r[1]*uz_surf_cr[1]; 

  double flux_energy_l[3] = {0.0}; 
  double flux_energy_r[3] = {0.0}; 

  flux_energy_l[0] = 0.2795084971874737*uy_surf_lr[2]*energy_l[8]+0.2795084971874737*uy_surf_cl[2]*energy_l[8]+0.2795084971874737*uy_surf_lr[2]*energy_c[8]+0.2795084971874737*uy_surf_cl[2]*energy_c[8]+0.2795084971874738*uy_surf_lr[1]*energy_l[7]+0.2795084971874738*uy_surf_cl[1]*energy_l[7]+0.2795084971874738*uy_surf_lr[1]*energy_c[7]+0.2795084971874738*uy_surf_cl[1]*energy_c[7]+0.2165063509461097*uy_surf_lr[2]*energy_l[6]+0.2165063509461097*uy_surf_cl[2]*energy_l[6]-0.2165063509461097*uy_surf_lr[2]*energy_c[6]-0.2165063509461097*uy_surf_cl[2]*energy_c[6]+0.2795084971874737*uy_surf_lr[0]*energy_l[5]+0.2795084971874737*uy_surf_cl[0]*energy_l[5]+0.2795084971874737*uy_surf_lr[0]*energy_c[5]+0.2795084971874737*uy_surf_cl[0]*energy_c[5]+0.125*uy_surf_lr[2]*energy_l[4]+0.125*uy_surf_cl[2]*energy_l[4]+0.125*uy_surf_lr[2]*energy_c[4]+0.125*uy_surf_cl[2]*energy_c[4]+0.2165063509461096*uy_surf_lr[1]*energy_l[3]+0.2165063509461096*uy_surf_cl[1]*energy_l[3]-0.2165063509461096*uy_surf_lr[1]*energy_c[3]-0.2165063509461096*uy_surf_cl[1]*energy_c[3]+0.1767766952966368*p_surf_lr[2]*uy_surf_lr[2]+0.1767766952966368*p_surf_cl[2]*uy_surf_lr[2]+0.1767766952966368*p_surf_lr[2]*uy_surf_cl[2]+0.1767766952966368*p_surf_cl[2]*uy_surf_cl[2]+0.2165063509461096*uy_surf_lr[0]*energy_l[2]+0.2165063509461096*uy_surf_cl[0]*energy_l[2]-0.2165063509461096*uy_surf_lr[0]*energy_c[2]-0.2165063509461096*uy_surf_cl[0]*energy_c[2]+0.1767766952966368*p_surf_lr[1]*uy_surf_lr[1]+0.1767766952966368*p_surf_cl[1]*uy_surf_lr[1]+0.125*energy_l[1]*uy_surf_lr[1]+0.125*energy_c[1]*uy_surf_lr[1]+0.1767766952966368*p_surf_lr[1]*uy_surf_cl[1]+0.1767766952966368*p_surf_cl[1]*uy_surf_cl[1]+0.125*energy_l[1]*uy_surf_cl[1]+0.125*energy_c[1]*uy_surf_cl[1]+0.1767766952966368*p_surf_lr[0]*uy_surf_lr[0]+0.1767766952966368*p_surf_cl[0]*uy_surf_lr[0]+0.125*energy_l[0]*uy_surf_lr[0]+0.125*energy_c[0]*uy_surf_lr[0]+0.1767766952966368*p_surf_lr[0]*uy_surf_cl[0]+0.1767766952966368*p_surf_cl[0]*uy_surf_cl[0]+0.125*energy_l[0]*uy_surf_cl[0]+0.125*energy_c[0]*uy_surf_cl[0]-0.5*apdq_energy_l[0]+0.5*amdq_energy_l[0]; 
  flux_energy_l[1] = 0.25*uy_surf_lr[1]*energy_l[8]+0.25*uy_surf_cl[1]*energy_l[8]+0.25*uy_surf_lr[1]*energy_c[8]+0.25*uy_surf_cl[1]*energy_c[8]+0.2500000000000001*uy_surf_lr[2]*energy_l[7]+0.2500000000000001*uy_surf_cl[2]*energy_l[7]+0.2795084971874738*uy_surf_lr[0]*energy_l[7]+0.2795084971874738*uy_surf_cl[0]*energy_l[7]+0.2500000000000001*uy_surf_lr[2]*energy_c[7]+0.2500000000000001*uy_surf_cl[2]*energy_c[7]+0.2795084971874738*uy_surf_lr[0]*energy_c[7]+0.2795084971874738*uy_surf_cl[0]*energy_c[7]+0.1936491673103709*uy_surf_lr[1]*energy_l[6]+0.1936491673103709*uy_surf_cl[1]*energy_l[6]-0.1936491673103709*uy_surf_lr[1]*energy_c[6]-0.1936491673103709*uy_surf_cl[1]*energy_c[6]+0.2795084971874737*uy_surf_lr[1]*energy_l[5]+0.2795084971874737*uy_surf_cl[1]*energy_l[5]+0.2795084971874737*uy_surf_lr[1]*energy_c[5]+0.2795084971874737*uy_surf_cl[1]*energy_c[5]+0.1118033988749895*uy_surf_lr[1]*energy_l[4]+0.1118033988749895*uy_surf_cl[1]*energy_l[4]+0.1118033988749895*uy_surf_lr[1]*energy_c[4]+0.1118033988749895*uy_surf_cl[1]*energy_c[4]+0.1936491673103708*uy_surf_lr[2]*energy_l[3]+0.1936491673103708*uy_surf_cl[2]*energy_l[3]+0.2165063509461096*uy_surf_lr[0]*energy_l[3]+0.2165063509461096*uy_surf_cl[0]*energy_l[3]-0.1936491673103708*uy_surf_lr[2]*energy_c[3]-0.1936491673103708*uy_surf_cl[2]*energy_c[3]-0.2165063509461096*uy_surf_lr[0]*energy_c[3]-0.2165063509461096*uy_surf_cl[0]*energy_c[3]+0.1581138830084189*p_surf_lr[1]*uy_surf_lr[2]+0.1581138830084189*p_surf_cl[1]*uy_surf_lr[2]+0.1118033988749895*energy_l[1]*uy_surf_lr[2]+0.1118033988749895*energy_c[1]*uy_surf_lr[2]+0.1581138830084189*p_surf_lr[1]*uy_surf_cl[2]+0.1581138830084189*p_surf_cl[1]*uy_surf_cl[2]+0.1118033988749895*energy_l[1]*uy_surf_cl[2]+0.1118033988749895*energy_c[1]*uy_surf_cl[2]+0.1581138830084189*uy_surf_lr[1]*p_surf_lr[2]+0.1581138830084189*uy_surf_cl[1]*p_surf_lr[2]+0.1581138830084189*uy_surf_lr[1]*p_surf_cl[2]+0.1581138830084189*uy_surf_cl[1]*p_surf_cl[2]+0.2165063509461096*uy_surf_lr[1]*energy_l[2]+0.2165063509461096*uy_surf_cl[1]*energy_l[2]-0.2165063509461096*uy_surf_lr[1]*energy_c[2]-0.2165063509461096*uy_surf_cl[1]*energy_c[2]+0.1767766952966368*p_surf_lr[0]*uy_surf_lr[1]+0.1767766952966368*p_surf_cl[0]*uy_surf_lr[1]+0.125*energy_l[0]*uy_surf_lr[1]+0.125*energy_c[0]*uy_surf_lr[1]+0.1767766952966368*p_surf_lr[0]*uy_surf_cl[1]+0.1767766952966368*p_surf_cl[0]*uy_surf_cl[1]+0.125*energy_l[0]*uy_surf_cl[1]+0.125*energy_c[0]*uy_surf_cl[1]+0.1767766952966368*uy_surf_lr[0]*p_surf_lr[1]+0.1767766952966368*uy_surf_cl[0]*p_surf_lr[1]+0.1767766952966368*uy_surf_lr[0]*p_surf_cl[1]+0.1767766952966368*uy_surf_cl[0]*p_surf_cl[1]+0.125*uy_surf_lr[0]*energy_l[1]+0.125*uy_surf_cl[0]*energy_l[1]+0.125*uy_surf_lr[0]*energy_c[1]+0.125*uy_surf_cl[0]*energy_c[1]-0.5*apdq_energy_l[1]+0.5*amdq_energy_l[1]; 
  flux_energy_l[2] = 0.1785714285714286*uy_surf_lr[2]*energy_l[8]+0.1785714285714286*uy_surf_cl[2]*energy_l[8]+0.2795084971874737*uy_surf_lr[0]*energy_l[8]+0.2795084971874737*uy_surf_cl[0]*energy_l[8]+0.1785714285714286*uy_surf_lr[2]*energy_c[8]+0.1785714285714286*uy_surf_cl[2]*energy_c[8]+0.2795084971874737*uy_surf_lr[0]*energy_c[8]+0.2795084971874737*uy_surf_cl[0]*energy_c[8]+0.2500000000000001*uy_surf_lr[1]*energy_l[7]+0.2500000000000001*uy_surf_cl[1]*energy_l[7]+0.2500000000000001*uy_surf_lr[1]*energy_c[7]+0.2500000000000001*uy_surf_cl[1]*energy_c[7]+0.138320833793122*uy_surf_lr[2]*energy_l[6]+0.138320833793122*uy_surf_cl[2]*energy_l[6]+0.2165063509461097*uy_surf_lr[0]*energy_l[6]+0.2165063509461097*uy_surf_cl[0]*energy_l[6]-0.138320833793122*uy_surf_lr[2]*energy_c[6]-0.138320833793122*uy_surf_cl[2]*energy_c[6]-0.2165063509461097*uy_surf_lr[0]*energy_c[6]-0.2165063509461097*uy_surf_cl[0]*energy_c[6]+0.2795084971874737*uy_surf_lr[2]*energy_l[5]+0.2795084971874737*uy_surf_cl[2]*energy_l[5]+0.2795084971874737*uy_surf_lr[2]*energy_c[5]+0.2795084971874737*uy_surf_cl[2]*energy_c[5]+0.07985957062499249*uy_surf_lr[2]*energy_l[4]+0.07985957062499249*uy_surf_cl[2]*energy_l[4]+0.125*uy_surf_lr[0]*energy_l[4]+0.125*uy_surf_cl[0]*energy_l[4]+0.07985957062499249*uy_surf_lr[2]*energy_c[4]+0.07985957062499249*uy_surf_cl[2]*energy_c[4]+0.125*uy_surf_lr[0]*energy_c[4]+0.125*uy_surf_cl[0]*energy_c[4]+0.1936491673103708*uy_surf_lr[1]*energy_l[3]+0.1936491673103708*uy_surf_cl[1]*energy_l[3]-0.1936491673103708*uy_surf_lr[1]*energy_c[3]-0.1936491673103708*uy_surf_cl[1]*energy_c[3]+0.1129384878631564*p_surf_lr[2]*uy_surf_lr[2]+0.1129384878631564*p_surf_cl[2]*uy_surf_lr[2]+0.2165063509461096*energy_l[2]*uy_surf_lr[2]-0.2165063509461096*energy_c[2]*uy_surf_lr[2]+0.1767766952966368*p_surf_lr[0]*uy_surf_lr[2]+0.1767766952966368*p_surf_cl[0]*uy_surf_lr[2]+0.125*energy_l[0]*uy_surf_lr[2]+0.125*energy_c[0]*uy_surf_lr[2]+0.1129384878631564*p_surf_lr[2]*uy_surf_cl[2]+0.1129384878631564*p_surf_cl[2]*uy_surf_cl[2]+0.2165063509461096*energy_l[2]*uy_surf_cl[2]-0.2165063509461096*energy_c[2]*uy_surf_cl[2]+0.1767766952966368*p_surf_lr[0]*uy_surf_cl[2]+0.1767766952966368*p_surf_cl[0]*uy_surf_cl[2]+0.125*energy_l[0]*uy_surf_cl[2]+0.125*energy_c[0]*uy_surf_cl[2]+0.1767766952966368*uy_surf_lr[0]*p_surf_lr[2]+0.1767766952966368*uy_surf_cl[0]*p_surf_lr[2]+0.1767766952966368*uy_surf_lr[0]*p_surf_cl[2]+0.1767766952966368*uy_surf_cl[0]*p_surf_cl[2]-0.5*apdq_energy_l[2]+0.5*amdq_energy_l[2]+0.1581138830084189*p_surf_lr[1]*uy_surf_lr[1]+0.1581138830084189*p_surf_cl[1]*uy_surf_lr[1]+0.1118033988749895*energy_l[1]*uy_surf_lr[1]+0.1118033988749895*energy_c[1]*uy_surf_lr[1]+0.1581138830084189*p_surf_lr[1]*uy_surf_cl[1]+0.1581138830084189*p_surf_cl[1]*uy_surf_cl[1]+0.1118033988749895*energy_l[1]*uy_surf_cl[1]+0.1118033988749895*energy_c[1]*uy_surf_cl[1]; 

  flux_energy_r[0] = 0.2795084971874737*uy_surf_rl[2]*energy_r[8]+0.2795084971874737*uy_surf_cr[2]*energy_r[8]+0.2795084971874737*uy_surf_rl[2]*energy_c[8]+0.2795084971874737*uy_surf_cr[2]*energy_c[8]+0.2795084971874738*uy_surf_rl[1]*energy_r[7]+0.2795084971874738*uy_surf_cr[1]*energy_r[7]+0.2795084971874738*uy_surf_rl[1]*energy_c[7]+0.2795084971874738*uy_surf_cr[1]*energy_c[7]-0.2165063509461097*uy_surf_rl[2]*energy_r[6]-0.2165063509461097*uy_surf_cr[2]*energy_r[6]+0.2165063509461097*uy_surf_rl[2]*energy_c[6]+0.2165063509461097*uy_surf_cr[2]*energy_c[6]+0.2795084971874737*uy_surf_rl[0]*energy_r[5]+0.2795084971874737*uy_surf_cr[0]*energy_r[5]+0.2795084971874737*uy_surf_rl[0]*energy_c[5]+0.2795084971874737*uy_surf_cr[0]*energy_c[5]+0.125*uy_surf_rl[2]*energy_r[4]+0.125*uy_surf_cr[2]*energy_r[4]+0.125*uy_surf_rl[2]*energy_c[4]+0.125*uy_surf_cr[2]*energy_c[4]-0.2165063509461096*uy_surf_rl[1]*energy_r[3]-0.2165063509461096*uy_surf_cr[1]*energy_r[3]+0.2165063509461096*uy_surf_rl[1]*energy_c[3]+0.2165063509461096*uy_surf_cr[1]*energy_c[3]+0.1767766952966368*p_surf_rl[2]*uy_surf_rl[2]+0.1767766952966368*p_surf_cr[2]*uy_surf_rl[2]+0.1767766952966368*p_surf_rl[2]*uy_surf_cr[2]+0.1767766952966368*p_surf_cr[2]*uy_surf_cr[2]-0.2165063509461096*uy_surf_rl[0]*energy_r[2]-0.2165063509461096*uy_surf_cr[0]*energy_r[2]+0.2165063509461096*uy_surf_rl[0]*energy_c[2]+0.2165063509461096*uy_surf_cr[0]*energy_c[2]+0.1767766952966368*p_surf_rl[1]*uy_surf_rl[1]+0.1767766952966368*p_surf_cr[1]*uy_surf_rl[1]+0.125*energy_r[1]*uy_surf_rl[1]+0.125*energy_c[1]*uy_surf_rl[1]+0.1767766952966368*p_surf_rl[1]*uy_surf_cr[1]+0.1767766952966368*p_surf_cr[1]*uy_surf_cr[1]+0.125*energy_r[1]*uy_surf_cr[1]+0.125*energy_c[1]*uy_surf_cr[1]+0.1767766952966368*p_surf_rl[0]*uy_surf_rl[0]+0.1767766952966368*p_surf_cr[0]*uy_surf_rl[0]+0.125*energy_r[0]*uy_surf_rl[0]+0.125*energy_c[0]*uy_surf_rl[0]+0.1767766952966368*p_surf_rl[0]*uy_surf_cr[0]+0.1767766952966368*p_surf_cr[0]*uy_surf_cr[0]+0.125*energy_r[0]*uy_surf_cr[0]+0.125*energy_c[0]*uy_surf_cr[0]-0.5*apdq_energy_r[0]+0.5*amdq_energy_r[0]; 
  flux_energy_r[1] = 0.25*uy_surf_rl[1]*energy_r[8]+0.25*uy_surf_cr[1]*energy_r[8]+0.25*uy_surf_rl[1]*energy_c[8]+0.25*uy_surf_cr[1]*energy_c[8]+0.2500000000000001*uy_surf_rl[2]*energy_r[7]+0.2500000000000001*uy_surf_cr[2]*energy_r[7]+0.2795084971874738*uy_surf_rl[0]*energy_r[7]+0.2795084971874738*uy_surf_cr[0]*energy_r[7]+0.2500000000000001*uy_surf_rl[2]*energy_c[7]+0.2500000000000001*uy_surf_cr[2]*energy_c[7]+0.2795084971874738*uy_surf_rl[0]*energy_c[7]+0.2795084971874738*uy_surf_cr[0]*energy_c[7]-0.1936491673103709*uy_surf_rl[1]*energy_r[6]-0.1936491673103709*uy_surf_cr[1]*energy_r[6]+0.1936491673103709*uy_surf_rl[1]*energy_c[6]+0.1936491673103709*uy_surf_cr[1]*energy_c[6]+0.2795084971874737*uy_surf_rl[1]*energy_r[5]+0.2795084971874737*uy_surf_cr[1]*energy_r[5]+0.2795084971874737*uy_surf_rl[1]*energy_c[5]+0.2795084971874737*uy_surf_cr[1]*energy_c[5]+0.1118033988749895*uy_surf_rl[1]*energy_r[4]+0.1118033988749895*uy_surf_cr[1]*energy_r[4]+0.1118033988749895*uy_surf_rl[1]*energy_c[4]+0.1118033988749895*uy_surf_cr[1]*energy_c[4]-0.1936491673103708*uy_surf_rl[2]*energy_r[3]-0.1936491673103708*uy_surf_cr[2]*energy_r[3]-0.2165063509461096*uy_surf_rl[0]*energy_r[3]-0.2165063509461096*uy_surf_cr[0]*energy_r[3]+0.1936491673103708*uy_surf_rl[2]*energy_c[3]+0.1936491673103708*uy_surf_cr[2]*energy_c[3]+0.2165063509461096*uy_surf_rl[0]*energy_c[3]+0.2165063509461096*uy_surf_cr[0]*energy_c[3]+0.1581138830084189*p_surf_rl[1]*uy_surf_rl[2]+0.1581138830084189*p_surf_cr[1]*uy_surf_rl[2]+0.1118033988749895*energy_r[1]*uy_surf_rl[2]+0.1118033988749895*energy_c[1]*uy_surf_rl[2]+0.1581138830084189*p_surf_rl[1]*uy_surf_cr[2]+0.1581138830084189*p_surf_cr[1]*uy_surf_cr[2]+0.1118033988749895*energy_r[1]*uy_surf_cr[2]+0.1118033988749895*energy_c[1]*uy_surf_cr[2]+0.1581138830084189*uy_surf_rl[1]*p_surf_rl[2]+0.1581138830084189*uy_surf_cr[1]*p_surf_rl[2]+0.1581138830084189*uy_surf_rl[1]*p_surf_cr[2]+0.1581138830084189*uy_surf_cr[1]*p_surf_cr[2]-0.2165063509461096*uy_surf_rl[1]*energy_r[2]-0.2165063509461096*uy_surf_cr[1]*energy_r[2]+0.2165063509461096*uy_surf_rl[1]*energy_c[2]+0.2165063509461096*uy_surf_cr[1]*energy_c[2]+0.1767766952966368*p_surf_rl[0]*uy_surf_rl[1]+0.1767766952966368*p_surf_cr[0]*uy_surf_rl[1]+0.125*energy_r[0]*uy_surf_rl[1]+0.125*energy_c[0]*uy_surf_rl[1]+0.1767766952966368*p_surf_rl[0]*uy_surf_cr[1]+0.1767766952966368*p_surf_cr[0]*uy_surf_cr[1]+0.125*energy_r[0]*uy_surf_cr[1]+0.125*energy_c[0]*uy_surf_cr[1]+0.1767766952966368*uy_surf_rl[0]*p_surf_rl[1]+0.1767766952966368*uy_surf_cr[0]*p_surf_rl[1]+0.1767766952966368*uy_surf_rl[0]*p_surf_cr[1]+0.1767766952966368*uy_surf_cr[0]*p_surf_cr[1]+0.125*uy_surf_rl[0]*energy_r[1]+0.125*uy_surf_cr[0]*energy_r[1]+0.125*uy_surf_rl[0]*energy_c[1]+0.125*uy_surf_cr[0]*energy_c[1]-0.5*apdq_energy_r[1]+0.5*amdq_energy_r[1]; 
  flux_energy_r[2] = 0.1785714285714286*uy_surf_rl[2]*energy_r[8]+0.1785714285714286*uy_surf_cr[2]*energy_r[8]+0.2795084971874737*uy_surf_rl[0]*energy_r[8]+0.2795084971874737*uy_surf_cr[0]*energy_r[8]+0.1785714285714286*uy_surf_rl[2]*energy_c[8]+0.1785714285714286*uy_surf_cr[2]*energy_c[8]+0.2795084971874737*uy_surf_rl[0]*energy_c[8]+0.2795084971874737*uy_surf_cr[0]*energy_c[8]+0.2500000000000001*uy_surf_rl[1]*energy_r[7]+0.2500000000000001*uy_surf_cr[1]*energy_r[7]+0.2500000000000001*uy_surf_rl[1]*energy_c[7]+0.2500000000000001*uy_surf_cr[1]*energy_c[7]-0.138320833793122*uy_surf_rl[2]*energy_r[6]-0.138320833793122*uy_surf_cr[2]*energy_r[6]-0.2165063509461097*uy_surf_rl[0]*energy_r[6]-0.2165063509461097*uy_surf_cr[0]*energy_r[6]+0.138320833793122*uy_surf_rl[2]*energy_c[6]+0.138320833793122*uy_surf_cr[2]*energy_c[6]+0.2165063509461097*uy_surf_rl[0]*energy_c[6]+0.2165063509461097*uy_surf_cr[0]*energy_c[6]+0.2795084971874737*uy_surf_rl[2]*energy_r[5]+0.2795084971874737*uy_surf_cr[2]*energy_r[5]+0.2795084971874737*uy_surf_rl[2]*energy_c[5]+0.2795084971874737*uy_surf_cr[2]*energy_c[5]+0.07985957062499249*uy_surf_rl[2]*energy_r[4]+0.07985957062499249*uy_surf_cr[2]*energy_r[4]+0.125*uy_surf_rl[0]*energy_r[4]+0.125*uy_surf_cr[0]*energy_r[4]+0.07985957062499249*uy_surf_rl[2]*energy_c[4]+0.07985957062499249*uy_surf_cr[2]*energy_c[4]+0.125*uy_surf_rl[0]*energy_c[4]+0.125*uy_surf_cr[0]*energy_c[4]-0.1936491673103708*uy_surf_rl[1]*energy_r[3]-0.1936491673103708*uy_surf_cr[1]*energy_r[3]+0.1936491673103708*uy_surf_rl[1]*energy_c[3]+0.1936491673103708*uy_surf_cr[1]*energy_c[3]+0.1129384878631564*p_surf_rl[2]*uy_surf_rl[2]+0.1129384878631564*p_surf_cr[2]*uy_surf_rl[2]-0.2165063509461096*energy_r[2]*uy_surf_rl[2]+0.2165063509461096*energy_c[2]*uy_surf_rl[2]+0.1767766952966368*p_surf_rl[0]*uy_surf_rl[2]+0.1767766952966368*p_surf_cr[0]*uy_surf_rl[2]+0.125*energy_r[0]*uy_surf_rl[2]+0.125*energy_c[0]*uy_surf_rl[2]+0.1129384878631564*p_surf_rl[2]*uy_surf_cr[2]+0.1129384878631564*p_surf_cr[2]*uy_surf_cr[2]-0.2165063509461096*energy_r[2]*uy_surf_cr[2]+0.2165063509461096*energy_c[2]*uy_surf_cr[2]+0.1767766952966368*p_surf_rl[0]*uy_surf_cr[2]+0.1767766952966368*p_surf_cr[0]*uy_surf_cr[2]+0.125*energy_r[0]*uy_surf_cr[2]+0.125*energy_c[0]*uy_surf_cr[2]+0.1767766952966368*uy_surf_rl[0]*p_surf_rl[2]+0.1767766952966368*uy_surf_cr[0]*p_surf_rl[2]+0.1767766952966368*uy_surf_rl[0]*p_surf_cr[2]+0.1767766952966368*uy_surf_cr[0]*p_surf_cr[2]-0.5*apdq_energy_r[2]+0.5*amdq_energy_r[2]+0.1581138830084189*p_surf_rl[1]*uy_surf_rl[1]+0.1581138830084189*p_surf_cr[1]*uy_surf_rl[1]+0.1118033988749895*energy_r[1]*uy_surf_rl[1]+0.1118033988749895*energy_c[1]*uy_surf_rl[1]+0.1581138830084189*p_surf_rl[1]*uy_surf_cr[1]+0.1581138830084189*p_surf_cr[1]*uy_surf_cr[1]+0.1118033988749895*energy_r[1]*uy_surf_cr[1]+0.1118033988749895*energy_c[1]*uy_surf_cr[1]; 

  outrho[0] += (0.7071067811865475*flux_rho_l[0]-0.7071067811865475*flux_rho_r[0])*dx1; 
  outrho[1] += (0.7071067811865475*flux_rho_l[1]-0.7071067811865475*flux_rho_r[1])*dx1; 
  outrho[2] += -1.224744871391589*(flux_rho_r[0]+flux_rho_l[0])*dx1; 
  outrho[3] += -1.224744871391589*(flux_rho_r[1]+flux_rho_l[1])*dx1; 
  outrho[4] += (0.7071067811865475*flux_rho_l[2]-0.7071067811865475*flux_rho_r[2])*dx1; 
  outrho[5] += (1.58113883008419*flux_rho_l[0]-1.58113883008419*flux_rho_r[0])*dx1; 
  outrho[6] += -1.224744871391589*(flux_rho_r[2]+flux_rho_l[2])*dx1; 
  outrho[7] += (1.58113883008419*flux_rho_l[1]-1.58113883008419*flux_rho_r[1])*dx1; 
  outrho[8] += (1.58113883008419*flux_rho_l[2]-1.58113883008419*flux_rho_r[2])*dx1; 

  outrhoux[0] += (0.7071067811865475*flux_rhoux_l[0]-0.7071067811865475*flux_rhoux_r[0])*dx1; 
  outrhoux[1] += (0.7071067811865475*flux_rhoux_l[1]-0.7071067811865475*flux_rhoux_r[1])*dx1; 
  outrhoux[2] += -1.224744871391589*(flux_rhoux_r[0]+flux_rhoux_l[0])*dx1; 
  outrhoux[3] += -1.224744871391589*(flux_rhoux_r[1]+flux_rhoux_l[1])*dx1; 
  outrhoux[4] += (0.7071067811865475*flux_rhoux_l[2]-0.7071067811865475*flux_rhoux_r[2])*dx1; 
  outrhoux[5] += (1.58113883008419*flux_rhoux_l[0]-1.58113883008419*flux_rhoux_r[0])*dx1; 
  outrhoux[6] += -1.224744871391589*(flux_rhoux_r[2]+flux_rhoux_l[2])*dx1; 
  outrhoux[7] += (1.58113883008419*flux_rhoux_l[1]-1.58113883008419*flux_rhoux_r[1])*dx1; 
  outrhoux[8] += (1.58113883008419*flux_rhoux_l[2]-1.58113883008419*flux_rhoux_r[2])*dx1; 

  outrhouy[0] += (0.7071067811865475*flux_rhouy_l[0]-0.7071067811865475*flux_rhouy_r[0])*dx1; 
  outrhouy[1] += (0.7071067811865475*flux_rhouy_l[1]-0.7071067811865475*flux_rhouy_r[1])*dx1; 
  outrhouy[2] += -1.224744871391589*(flux_rhouy_r[0]+flux_rhouy_l[0])*dx1; 
  outrhouy[3] += -1.224744871391589*(flux_rhouy_r[1]+flux_rhouy_l[1])*dx1; 
  outrhouy[4] += (0.7071067811865475*flux_rhouy_l[2]-0.7071067811865475*flux_rhouy_r[2])*dx1; 
  outrhouy[5] += (1.58113883008419*flux_rhouy_l[0]-1.58113883008419*flux_rhouy_r[0])*dx1; 
  outrhouy[6] += -1.224744871391589*(flux_rhouy_r[2]+flux_rhouy_l[2])*dx1; 
  outrhouy[7] += (1.58113883008419*flux_rhouy_l[1]-1.58113883008419*flux_rhouy_r[1])*dx1; 
  outrhouy[8] += (1.58113883008419*flux_rhouy_l[2]-1.58113883008419*flux_rhouy_r[2])*dx1; 

  outrhouz[0] += (0.7071067811865475*flux_rhouz_l[0]-0.7071067811865475*flux_rhouz_r[0])*dx1; 
  outrhouz[1] += (0.7071067811865475*flux_rhouz_l[1]-0.7071067811865475*flux_rhouz_r[1])*dx1; 
  outrhouz[2] += -1.224744871391589*(flux_rhouz_r[0]+flux_rhouz_l[0])*dx1; 
  outrhouz[3] += -1.224744871391589*(flux_rhouz_r[1]+flux_rhouz_l[1])*dx1; 
  outrhouz[4] += (0.7071067811865475*flux_rhouz_l[2]-0.7071067811865475*flux_rhouz_r[2])*dx1; 
  outrhouz[5] += (1.58113883008419*flux_rhouz_l[0]-1.58113883008419*flux_rhouz_r[0])*dx1; 
  outrhouz[6] += -1.224744871391589*(flux_rhouz_r[2]+flux_rhouz_l[2])*dx1; 
  outrhouz[7] += (1.58113883008419*flux_rhouz_l[1]-1.58113883008419*flux_rhouz_r[1])*dx1; 
  outrhouz[8] += (1.58113883008419*flux_rhouz_l[2]-1.58113883008419*flux_rhouz_r[2])*dx1; 

  outenergy[0] += (0.7071067811865475*flux_energy_l[0]-0.7071067811865475*flux_energy_r[0])*dx1; 
  outenergy[1] += (0.7071067811865475*flux_energy_l[1]-0.7071067811865475*flux_energy_r[1])*dx1; 
  outenergy[2] += -1.224744871391589*(flux_energy_r[0]+flux_energy_l[0])*dx1; 
  outenergy[3] += -1.224744871391589*(flux_energy_r[1]+flux_energy_l[1])*dx1; 
  outenergy[4] += (0.7071067811865475*flux_energy_l[2]-0.7071067811865475*flux_energy_r[2])*dx1; 
  outenergy[5] += (1.58113883008419*flux_energy_l[0]-1.58113883008419*flux_energy_r[0])*dx1; 
  outenergy[6] += -1.224744871391589*(flux_energy_r[2]+flux_energy_l[2])*dx1; 
  outenergy[7] += (1.58113883008419*flux_energy_l[1]-1.58113883008419*flux_energy_r[1])*dx1; 
  outenergy[8] += (1.58113883008419*flux_energy_l[2]-1.58113883008419*flux_energy_r[2])*dx1; 

  return 0.;

} 
