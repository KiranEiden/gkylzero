#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_basis_ser_3x_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_ser_3x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double euler_pkpm_surfy_3x_ser_p1(const double *w, const double *dxv,
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r,
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:                Cell-center coordinates.
  // dxv[NDIM]:              Cell spacing.
  // wv_eqn:                 Wave equation for computing fluctuations at the interface for upwinding.
  // geom_l:                 Geometry for the left surface update.
  // geom_r:                 Geometry for the right surface update.
  // vlasov_pkpm_moms_l/c/r: Input pkpm moments in left/center/right cells.
  // prim_surf_l/c/r:        Input surface primitive variables [u_i, 3*T_ii/m] in left/center/right cells in each direction.
  // p_ij_l/c/r:             Input volume expansion of p_ij = (p_par - p_perp) b_i b_j + p_perp g_ij in left/center/right cells.
  // euler_pkpm_l/c/r:       Input [rho ux, rho uy, rho uz], Fluid input state vector in left/center/right cells.
  // pkpm_lax:               Surface expansion of pkpm Lax penalization: lambda_i = |u_i| + sqrt(3.0*T_ii/m).
  // out: Incremented output.

  const double dx1 = 2.0/dxv[1]; 

  const double *rhoux_l = &euler_pkpm_l[0]; 
  const double *rhouy_l = &euler_pkpm_l[8]; 
  const double *rhouz_l = &euler_pkpm_l[16]; 

  const double *rhoux_c = &euler_pkpm_c[0]; 
  const double *rhouy_c = &euler_pkpm_c[8]; 
  const double *rhouz_c = &euler_pkpm_c[16]; 

  const double *rhoux_r = &euler_pkpm_r[0]; 
  const double *rhouy_r = &euler_pkpm_r[8]; 
  const double *rhouz_r = &euler_pkpm_r[16]; 

  const double *rho_l = &vlasov_pkpm_moms_l[0]; 
  const double *rho_c = &vlasov_pkpm_moms_c[0]; 
  const double *rho_r = &vlasov_pkpm_moms_r[0]; 

  const double *Pxx_l = &p_ij_l[0]; 
  const double *Pxy_l = &p_ij_l[8]; 
  const double *Pxz_l = &p_ij_l[16]; 
  const double *Pyy_l = &p_ij_l[24]; 
  const double *Pyz_l = &p_ij_l[32]; 
  const double *Pzz_l = &p_ij_l[40]; 

  const double *Pxx_c = &p_ij_c[0]; 
  const double *Pxy_c = &p_ij_c[8]; 
  const double *Pxz_c = &p_ij_c[16]; 
  const double *Pyy_c = &p_ij_c[24]; 
  const double *Pyz_c = &p_ij_c[32]; 
  const double *Pzz_c = &p_ij_c[40]; 

  const double *Pxx_r = &p_ij_r[0]; 
  const double *Pxy_r = &p_ij_r[8]; 
  const double *Pxz_r = &p_ij_r[16]; 
  const double *Pyy_r = &p_ij_r[24]; 
  const double *Pyz_r = &p_ij_r[32]; 
  const double *Pzz_r = &p_ij_r[40]; 

  const double *ux_surf_lr = &prim_surf_l[36]; 
  const double *uy_surf_lr = &prim_surf_l[44]; 
  const double *uz_surf_lr = &prim_surf_l[52]; 

  const double *ux_surf_cl = &prim_surf_c[32]; 
  const double *uy_surf_cl = &prim_surf_c[40]; 
  const double *uz_surf_cl = &prim_surf_c[48]; 

  const double *ux_surf_cr = &prim_surf_c[36]; 
  const double *uy_surf_cr = &prim_surf_c[44]; 
  const double *uz_surf_cr = &prim_surf_c[52]; 

  const double *ux_surf_rl = &prim_surf_r[32]; 
  const double *uy_surf_rl = &prim_surf_r[40]; 
  const double *uz_surf_rl = &prim_surf_r[48]; 

  const double *pkpm_lax_l = &pkpm_lax[8]; 
  const double *pkpm_lax_r = &pkpm_lax[12]; 

  double *outrhou0 = &out[0]; 
  double *outrhou1 = &out[8]; 
  double *outrhou2 = &out[16]; 

  double flux_rho_l[4] = {0.0}; 
  double flux_rho_r[4] = {0.0}; 

  double avg_p_ij_x_l[4] = {0.0}; 
  double avg_p_ij_x_r[4] = {0.0}; 
  double avg_p_ij_y_l[4] = {0.0}; 
  double avg_p_ij_y_r[4] = {0.0}; 
  double avg_p_ij_z_l[4] = {0.0}; 
  double avg_p_ij_z_r[4] = {0.0}; 

  flux_rho_l[0] = 0.1530931089239486*uy_surf_lr[3]*rho_l[7]+0.1530931089239486*uy_surf_cl[3]*rho_l[7]+0.3061862178478971*pkpm_lax_l[3]*rho_l[7]-0.1530931089239486*uy_surf_lr[3]*rho_c[7]-0.1530931089239486*uy_surf_cl[3]*rho_c[7]+0.3061862178478971*pkpm_lax_l[3]*rho_c[7]+0.1530931089239486*uy_surf_lr[2]*rho_l[6]+0.1530931089239486*uy_surf_cl[2]*rho_l[6]+0.3061862178478971*pkpm_lax_l[2]*rho_l[6]-0.1530931089239486*uy_surf_lr[2]*rho_c[6]-0.1530931089239486*uy_surf_cl[2]*rho_c[6]+0.3061862178478971*pkpm_lax_l[2]*rho_c[6]+0.0883883476483184*uy_surf_lr[3]*rho_l[5]+0.0883883476483184*uy_surf_cl[3]*rho_l[5]+0.1767766952966368*pkpm_lax_l[3]*rho_l[5]+0.0883883476483184*uy_surf_lr[3]*rho_c[5]+0.0883883476483184*uy_surf_cl[3]*rho_c[5]-0.1767766952966368*pkpm_lax_l[3]*rho_c[5]+0.1530931089239486*uy_surf_lr[1]*rho_l[4]+0.1530931089239486*uy_surf_cl[1]*rho_l[4]+0.3061862178478971*pkpm_lax_l[1]*rho_l[4]-0.1530931089239486*uy_surf_lr[1]*rho_c[4]-0.1530931089239486*uy_surf_cl[1]*rho_c[4]+0.3061862178478971*pkpm_lax_l[1]*rho_c[4]+0.0883883476483184*uy_surf_lr[2]*rho_l[3]+0.0883883476483184*uy_surf_cl[2]*rho_l[3]+0.1767766952966368*pkpm_lax_l[2]*rho_l[3]+0.0883883476483184*uy_surf_lr[2]*rho_c[3]+0.0883883476483184*uy_surf_cl[2]*rho_c[3]-0.1767766952966368*pkpm_lax_l[2]*rho_c[3]+0.1530931089239486*uy_surf_lr[0]*rho_l[2]+0.1530931089239486*uy_surf_cl[0]*rho_l[2]+0.3061862178478971*pkpm_lax_l[0]*rho_l[2]-0.1530931089239486*uy_surf_lr[0]*rho_c[2]-0.1530931089239486*uy_surf_cl[0]*rho_c[2]+0.3061862178478971*pkpm_lax_l[0]*rho_c[2]+0.0883883476483184*rho_l[1]*uy_surf_lr[1]+0.0883883476483184*rho_c[1]*uy_surf_lr[1]+0.0883883476483184*rho_l[1]*uy_surf_cl[1]+0.0883883476483184*rho_c[1]*uy_surf_cl[1]+0.1767766952966368*pkpm_lax_l[1]*rho_l[1]-0.1767766952966368*pkpm_lax_l[1]*rho_c[1]+0.0883883476483184*rho_l[0]*uy_surf_lr[0]+0.0883883476483184*rho_c[0]*uy_surf_lr[0]+0.0883883476483184*rho_l[0]*uy_surf_cl[0]+0.0883883476483184*rho_c[0]*uy_surf_cl[0]+0.1767766952966368*pkpm_lax_l[0]*rho_l[0]-0.1767766952966368*pkpm_lax_l[0]*rho_c[0]; 
  flux_rho_l[1] = 0.1530931089239486*uy_surf_lr[2]*rho_l[7]+0.1530931089239486*uy_surf_cl[2]*rho_l[7]+0.3061862178478971*pkpm_lax_l[2]*rho_l[7]-0.1530931089239486*uy_surf_lr[2]*rho_c[7]-0.1530931089239486*uy_surf_cl[2]*rho_c[7]+0.3061862178478971*pkpm_lax_l[2]*rho_c[7]+0.1530931089239486*uy_surf_lr[3]*rho_l[6]+0.1530931089239486*uy_surf_cl[3]*rho_l[6]+0.3061862178478971*pkpm_lax_l[3]*rho_l[6]-0.1530931089239486*uy_surf_lr[3]*rho_c[6]-0.1530931089239486*uy_surf_cl[3]*rho_c[6]+0.3061862178478971*pkpm_lax_l[3]*rho_c[6]+0.0883883476483184*uy_surf_lr[2]*rho_l[5]+0.0883883476483184*uy_surf_cl[2]*rho_l[5]+0.1767766952966368*pkpm_lax_l[2]*rho_l[5]+0.0883883476483184*uy_surf_lr[2]*rho_c[5]+0.0883883476483184*uy_surf_cl[2]*rho_c[5]-0.1767766952966368*pkpm_lax_l[2]*rho_c[5]+0.1530931089239486*uy_surf_lr[0]*rho_l[4]+0.1530931089239486*uy_surf_cl[0]*rho_l[4]+0.3061862178478971*pkpm_lax_l[0]*rho_l[4]-0.1530931089239486*uy_surf_lr[0]*rho_c[4]-0.1530931089239486*uy_surf_cl[0]*rho_c[4]+0.3061862178478971*pkpm_lax_l[0]*rho_c[4]+0.0883883476483184*rho_l[3]*uy_surf_lr[3]+0.0883883476483184*rho_c[3]*uy_surf_lr[3]+0.0883883476483184*rho_l[3]*uy_surf_cl[3]+0.0883883476483184*rho_c[3]*uy_surf_cl[3]+0.1767766952966368*pkpm_lax_l[3]*rho_l[3]-0.1767766952966368*pkpm_lax_l[3]*rho_c[3]+0.1530931089239486*uy_surf_lr[1]*rho_l[2]+0.1530931089239486*uy_surf_cl[1]*rho_l[2]+0.3061862178478971*pkpm_lax_l[1]*rho_l[2]-0.1530931089239486*uy_surf_lr[1]*rho_c[2]-0.1530931089239486*uy_surf_cl[1]*rho_c[2]+0.3061862178478971*pkpm_lax_l[1]*rho_c[2]+0.0883883476483184*rho_l[0]*uy_surf_lr[1]+0.0883883476483184*rho_c[0]*uy_surf_lr[1]+0.0883883476483184*rho_l[0]*uy_surf_cl[1]+0.0883883476483184*rho_c[0]*uy_surf_cl[1]+0.0883883476483184*uy_surf_lr[0]*rho_l[1]+0.0883883476483184*uy_surf_cl[0]*rho_l[1]+0.1767766952966368*pkpm_lax_l[0]*rho_l[1]+0.0883883476483184*uy_surf_lr[0]*rho_c[1]+0.0883883476483184*uy_surf_cl[0]*rho_c[1]-0.1767766952966368*pkpm_lax_l[0]*rho_c[1]+0.1767766952966368*rho_l[0]*pkpm_lax_l[1]-0.1767766952966368*rho_c[0]*pkpm_lax_l[1]; 
  flux_rho_l[2] = 0.1530931089239486*uy_surf_lr[1]*rho_l[7]+0.1530931089239486*uy_surf_cl[1]*rho_l[7]+0.3061862178478971*pkpm_lax_l[1]*rho_l[7]-0.1530931089239486*uy_surf_lr[1]*rho_c[7]-0.1530931089239486*uy_surf_cl[1]*rho_c[7]+0.3061862178478971*pkpm_lax_l[1]*rho_c[7]+0.1530931089239486*uy_surf_lr[0]*rho_l[6]+0.1530931089239486*uy_surf_cl[0]*rho_l[6]+0.3061862178478971*pkpm_lax_l[0]*rho_l[6]-0.1530931089239486*uy_surf_lr[0]*rho_c[6]-0.1530931089239486*uy_surf_cl[0]*rho_c[6]+0.3061862178478971*pkpm_lax_l[0]*rho_c[6]+0.0883883476483184*uy_surf_lr[1]*rho_l[5]+0.0883883476483184*uy_surf_cl[1]*rho_l[5]+0.1767766952966368*pkpm_lax_l[1]*rho_l[5]+0.0883883476483184*uy_surf_lr[1]*rho_c[5]+0.0883883476483184*uy_surf_cl[1]*rho_c[5]-0.1767766952966368*pkpm_lax_l[1]*rho_c[5]+0.1530931089239486*uy_surf_lr[3]*rho_l[4]+0.1530931089239486*uy_surf_cl[3]*rho_l[4]+0.3061862178478971*pkpm_lax_l[3]*rho_l[4]-0.1530931089239486*uy_surf_lr[3]*rho_c[4]-0.1530931089239486*uy_surf_cl[3]*rho_c[4]+0.3061862178478971*pkpm_lax_l[3]*rho_c[4]+0.0883883476483184*rho_l[1]*uy_surf_lr[3]+0.0883883476483184*rho_c[1]*uy_surf_lr[3]+0.0883883476483184*rho_l[1]*uy_surf_cl[3]+0.0883883476483184*rho_c[1]*uy_surf_cl[3]+0.0883883476483184*uy_surf_lr[0]*rho_l[3]+0.0883883476483184*uy_surf_cl[0]*rho_l[3]+0.1767766952966368*pkpm_lax_l[0]*rho_l[3]+0.0883883476483184*uy_surf_lr[0]*rho_c[3]+0.0883883476483184*uy_surf_cl[0]*rho_c[3]-0.1767766952966368*pkpm_lax_l[0]*rho_c[3]+0.1767766952966368*rho_l[1]*pkpm_lax_l[3]-0.1767766952966368*rho_c[1]*pkpm_lax_l[3]+0.1530931089239486*rho_l[2]*uy_surf_lr[2]-0.1530931089239486*rho_c[2]*uy_surf_lr[2]+0.0883883476483184*rho_l[0]*uy_surf_lr[2]+0.0883883476483184*rho_c[0]*uy_surf_lr[2]+0.1530931089239486*rho_l[2]*uy_surf_cl[2]-0.1530931089239486*rho_c[2]*uy_surf_cl[2]+0.0883883476483184*rho_l[0]*uy_surf_cl[2]+0.0883883476483184*rho_c[0]*uy_surf_cl[2]+0.3061862178478971*pkpm_lax_l[2]*rho_l[2]+0.3061862178478971*pkpm_lax_l[2]*rho_c[2]+0.1767766952966368*rho_l[0]*pkpm_lax_l[2]-0.1767766952966368*rho_c[0]*pkpm_lax_l[2]; 
  flux_rho_l[3] = 0.1530931089239486*uy_surf_lr[0]*rho_l[7]+0.1530931089239486*uy_surf_cl[0]*rho_l[7]+0.3061862178478971*pkpm_lax_l[0]*rho_l[7]-0.1530931089239486*uy_surf_lr[0]*rho_c[7]-0.1530931089239486*uy_surf_cl[0]*rho_c[7]+0.3061862178478971*pkpm_lax_l[0]*rho_c[7]+0.1530931089239486*uy_surf_lr[1]*rho_l[6]+0.1530931089239486*uy_surf_cl[1]*rho_l[6]+0.3061862178478971*pkpm_lax_l[1]*rho_l[6]-0.1530931089239486*uy_surf_lr[1]*rho_c[6]-0.1530931089239486*uy_surf_cl[1]*rho_c[6]+0.3061862178478971*pkpm_lax_l[1]*rho_c[6]+0.0883883476483184*uy_surf_lr[0]*rho_l[5]+0.0883883476483184*uy_surf_cl[0]*rho_l[5]+0.1767766952966368*pkpm_lax_l[0]*rho_l[5]+0.0883883476483184*uy_surf_lr[0]*rho_c[5]+0.0883883476483184*uy_surf_cl[0]*rho_c[5]-0.1767766952966368*pkpm_lax_l[0]*rho_c[5]+0.1530931089239486*uy_surf_lr[2]*rho_l[4]+0.1530931089239486*uy_surf_cl[2]*rho_l[4]+0.3061862178478971*pkpm_lax_l[2]*rho_l[4]-0.1530931089239486*uy_surf_lr[2]*rho_c[4]-0.1530931089239486*uy_surf_cl[2]*rho_c[4]+0.3061862178478971*pkpm_lax_l[2]*rho_c[4]+0.1530931089239486*rho_l[2]*uy_surf_lr[3]-0.1530931089239486*rho_c[2]*uy_surf_lr[3]+0.0883883476483184*rho_l[0]*uy_surf_lr[3]+0.0883883476483184*rho_c[0]*uy_surf_lr[3]+0.1530931089239486*rho_l[2]*uy_surf_cl[3]-0.1530931089239486*rho_c[2]*uy_surf_cl[3]+0.0883883476483184*rho_l[0]*uy_surf_cl[3]+0.0883883476483184*rho_c[0]*uy_surf_cl[3]+0.0883883476483184*uy_surf_lr[1]*rho_l[3]+0.0883883476483184*uy_surf_cl[1]*rho_l[3]+0.1767766952966368*pkpm_lax_l[1]*rho_l[3]+0.0883883476483184*uy_surf_lr[1]*rho_c[3]+0.0883883476483184*uy_surf_cl[1]*rho_c[3]-0.1767766952966368*pkpm_lax_l[1]*rho_c[3]+0.3061862178478971*rho_l[2]*pkpm_lax_l[3]+0.3061862178478971*rho_c[2]*pkpm_lax_l[3]+0.1767766952966368*rho_l[0]*pkpm_lax_l[3]-0.1767766952966368*rho_c[0]*pkpm_lax_l[3]+0.0883883476483184*rho_l[1]*uy_surf_lr[2]+0.0883883476483184*rho_c[1]*uy_surf_lr[2]+0.0883883476483184*rho_l[1]*uy_surf_cl[2]+0.0883883476483184*rho_c[1]*uy_surf_cl[2]+0.1767766952966368*rho_l[1]*pkpm_lax_l[2]-0.1767766952966368*rho_c[1]*pkpm_lax_l[2]; 

  flux_rho_r[0] = (-0.1530931089239486*uy_surf_rl[3]*rho_r[7])-0.1530931089239486*uy_surf_cr[3]*rho_r[7]+0.3061862178478971*pkpm_lax_r[3]*rho_r[7]+0.1530931089239486*uy_surf_rl[3]*rho_c[7]+0.1530931089239486*uy_surf_cr[3]*rho_c[7]+0.3061862178478971*pkpm_lax_r[3]*rho_c[7]-0.1530931089239486*uy_surf_rl[2]*rho_r[6]-0.1530931089239486*uy_surf_cr[2]*rho_r[6]+0.3061862178478971*pkpm_lax_r[2]*rho_r[6]+0.1530931089239486*uy_surf_rl[2]*rho_c[6]+0.1530931089239486*uy_surf_cr[2]*rho_c[6]+0.3061862178478971*pkpm_lax_r[2]*rho_c[6]+0.0883883476483184*uy_surf_rl[3]*rho_r[5]+0.0883883476483184*uy_surf_cr[3]*rho_r[5]-0.1767766952966368*pkpm_lax_r[3]*rho_r[5]+0.0883883476483184*uy_surf_rl[3]*rho_c[5]+0.0883883476483184*uy_surf_cr[3]*rho_c[5]+0.1767766952966368*pkpm_lax_r[3]*rho_c[5]-0.1530931089239486*uy_surf_rl[1]*rho_r[4]-0.1530931089239486*uy_surf_cr[1]*rho_r[4]+0.3061862178478971*pkpm_lax_r[1]*rho_r[4]+0.1530931089239486*uy_surf_rl[1]*rho_c[4]+0.1530931089239486*uy_surf_cr[1]*rho_c[4]+0.3061862178478971*pkpm_lax_r[1]*rho_c[4]+0.0883883476483184*uy_surf_rl[2]*rho_r[3]+0.0883883476483184*uy_surf_cr[2]*rho_r[3]-0.1767766952966368*pkpm_lax_r[2]*rho_r[3]+0.0883883476483184*uy_surf_rl[2]*rho_c[3]+0.0883883476483184*uy_surf_cr[2]*rho_c[3]+0.1767766952966368*pkpm_lax_r[2]*rho_c[3]-0.1530931089239486*uy_surf_rl[0]*rho_r[2]-0.1530931089239486*uy_surf_cr[0]*rho_r[2]+0.3061862178478971*pkpm_lax_r[0]*rho_r[2]+0.1530931089239486*uy_surf_rl[0]*rho_c[2]+0.1530931089239486*uy_surf_cr[0]*rho_c[2]+0.3061862178478971*pkpm_lax_r[0]*rho_c[2]+0.0883883476483184*rho_r[1]*uy_surf_rl[1]+0.0883883476483184*rho_c[1]*uy_surf_rl[1]+0.0883883476483184*rho_r[1]*uy_surf_cr[1]+0.0883883476483184*rho_c[1]*uy_surf_cr[1]-0.1767766952966368*pkpm_lax_r[1]*rho_r[1]+0.1767766952966368*pkpm_lax_r[1]*rho_c[1]+0.0883883476483184*rho_r[0]*uy_surf_rl[0]+0.0883883476483184*rho_c[0]*uy_surf_rl[0]+0.0883883476483184*rho_r[0]*uy_surf_cr[0]+0.0883883476483184*rho_c[0]*uy_surf_cr[0]-0.1767766952966368*pkpm_lax_r[0]*rho_r[0]+0.1767766952966368*pkpm_lax_r[0]*rho_c[0]; 
  flux_rho_r[1] = (-0.1530931089239486*uy_surf_rl[2]*rho_r[7])-0.1530931089239486*uy_surf_cr[2]*rho_r[7]+0.3061862178478971*pkpm_lax_r[2]*rho_r[7]+0.1530931089239486*uy_surf_rl[2]*rho_c[7]+0.1530931089239486*uy_surf_cr[2]*rho_c[7]+0.3061862178478971*pkpm_lax_r[2]*rho_c[7]-0.1530931089239486*uy_surf_rl[3]*rho_r[6]-0.1530931089239486*uy_surf_cr[3]*rho_r[6]+0.3061862178478971*pkpm_lax_r[3]*rho_r[6]+0.1530931089239486*uy_surf_rl[3]*rho_c[6]+0.1530931089239486*uy_surf_cr[3]*rho_c[6]+0.3061862178478971*pkpm_lax_r[3]*rho_c[6]+0.0883883476483184*uy_surf_rl[2]*rho_r[5]+0.0883883476483184*uy_surf_cr[2]*rho_r[5]-0.1767766952966368*pkpm_lax_r[2]*rho_r[5]+0.0883883476483184*uy_surf_rl[2]*rho_c[5]+0.0883883476483184*uy_surf_cr[2]*rho_c[5]+0.1767766952966368*pkpm_lax_r[2]*rho_c[5]-0.1530931089239486*uy_surf_rl[0]*rho_r[4]-0.1530931089239486*uy_surf_cr[0]*rho_r[4]+0.3061862178478971*pkpm_lax_r[0]*rho_r[4]+0.1530931089239486*uy_surf_rl[0]*rho_c[4]+0.1530931089239486*uy_surf_cr[0]*rho_c[4]+0.3061862178478971*pkpm_lax_r[0]*rho_c[4]+0.0883883476483184*rho_r[3]*uy_surf_rl[3]+0.0883883476483184*rho_c[3]*uy_surf_rl[3]+0.0883883476483184*rho_r[3]*uy_surf_cr[3]+0.0883883476483184*rho_c[3]*uy_surf_cr[3]-0.1767766952966368*pkpm_lax_r[3]*rho_r[3]+0.1767766952966368*pkpm_lax_r[3]*rho_c[3]-0.1530931089239486*uy_surf_rl[1]*rho_r[2]-0.1530931089239486*uy_surf_cr[1]*rho_r[2]+0.3061862178478971*pkpm_lax_r[1]*rho_r[2]+0.1530931089239486*uy_surf_rl[1]*rho_c[2]+0.1530931089239486*uy_surf_cr[1]*rho_c[2]+0.3061862178478971*pkpm_lax_r[1]*rho_c[2]+0.0883883476483184*rho_r[0]*uy_surf_rl[1]+0.0883883476483184*rho_c[0]*uy_surf_rl[1]+0.0883883476483184*rho_r[0]*uy_surf_cr[1]+0.0883883476483184*rho_c[0]*uy_surf_cr[1]+0.0883883476483184*uy_surf_rl[0]*rho_r[1]+0.0883883476483184*uy_surf_cr[0]*rho_r[1]-0.1767766952966368*pkpm_lax_r[0]*rho_r[1]+0.0883883476483184*uy_surf_rl[0]*rho_c[1]+0.0883883476483184*uy_surf_cr[0]*rho_c[1]+0.1767766952966368*pkpm_lax_r[0]*rho_c[1]-0.1767766952966368*rho_r[0]*pkpm_lax_r[1]+0.1767766952966368*rho_c[0]*pkpm_lax_r[1]; 
  flux_rho_r[2] = (-0.1530931089239486*uy_surf_rl[1]*rho_r[7])-0.1530931089239486*uy_surf_cr[1]*rho_r[7]+0.3061862178478971*pkpm_lax_r[1]*rho_r[7]+0.1530931089239486*uy_surf_rl[1]*rho_c[7]+0.1530931089239486*uy_surf_cr[1]*rho_c[7]+0.3061862178478971*pkpm_lax_r[1]*rho_c[7]-0.1530931089239486*uy_surf_rl[0]*rho_r[6]-0.1530931089239486*uy_surf_cr[0]*rho_r[6]+0.3061862178478971*pkpm_lax_r[0]*rho_r[6]+0.1530931089239486*uy_surf_rl[0]*rho_c[6]+0.1530931089239486*uy_surf_cr[0]*rho_c[6]+0.3061862178478971*pkpm_lax_r[0]*rho_c[6]+0.0883883476483184*uy_surf_rl[1]*rho_r[5]+0.0883883476483184*uy_surf_cr[1]*rho_r[5]-0.1767766952966368*pkpm_lax_r[1]*rho_r[5]+0.0883883476483184*uy_surf_rl[1]*rho_c[5]+0.0883883476483184*uy_surf_cr[1]*rho_c[5]+0.1767766952966368*pkpm_lax_r[1]*rho_c[5]-0.1530931089239486*uy_surf_rl[3]*rho_r[4]-0.1530931089239486*uy_surf_cr[3]*rho_r[4]+0.3061862178478971*pkpm_lax_r[3]*rho_r[4]+0.1530931089239486*uy_surf_rl[3]*rho_c[4]+0.1530931089239486*uy_surf_cr[3]*rho_c[4]+0.3061862178478971*pkpm_lax_r[3]*rho_c[4]+0.0883883476483184*rho_r[1]*uy_surf_rl[3]+0.0883883476483184*rho_c[1]*uy_surf_rl[3]+0.0883883476483184*rho_r[1]*uy_surf_cr[3]+0.0883883476483184*rho_c[1]*uy_surf_cr[3]+0.0883883476483184*uy_surf_rl[0]*rho_r[3]+0.0883883476483184*uy_surf_cr[0]*rho_r[3]-0.1767766952966368*pkpm_lax_r[0]*rho_r[3]+0.0883883476483184*uy_surf_rl[0]*rho_c[3]+0.0883883476483184*uy_surf_cr[0]*rho_c[3]+0.1767766952966368*pkpm_lax_r[0]*rho_c[3]-0.1767766952966368*rho_r[1]*pkpm_lax_r[3]+0.1767766952966368*rho_c[1]*pkpm_lax_r[3]-0.1530931089239486*rho_r[2]*uy_surf_rl[2]+0.1530931089239486*rho_c[2]*uy_surf_rl[2]+0.0883883476483184*rho_r[0]*uy_surf_rl[2]+0.0883883476483184*rho_c[0]*uy_surf_rl[2]-0.1530931089239486*rho_r[2]*uy_surf_cr[2]+0.1530931089239486*rho_c[2]*uy_surf_cr[2]+0.0883883476483184*rho_r[0]*uy_surf_cr[2]+0.0883883476483184*rho_c[0]*uy_surf_cr[2]+0.3061862178478971*pkpm_lax_r[2]*rho_r[2]+0.3061862178478971*pkpm_lax_r[2]*rho_c[2]-0.1767766952966368*rho_r[0]*pkpm_lax_r[2]+0.1767766952966368*rho_c[0]*pkpm_lax_r[2]; 
  flux_rho_r[3] = (-0.1530931089239486*uy_surf_rl[0]*rho_r[7])-0.1530931089239486*uy_surf_cr[0]*rho_r[7]+0.3061862178478971*pkpm_lax_r[0]*rho_r[7]+0.1530931089239486*uy_surf_rl[0]*rho_c[7]+0.1530931089239486*uy_surf_cr[0]*rho_c[7]+0.3061862178478971*pkpm_lax_r[0]*rho_c[7]-0.1530931089239486*uy_surf_rl[1]*rho_r[6]-0.1530931089239486*uy_surf_cr[1]*rho_r[6]+0.3061862178478971*pkpm_lax_r[1]*rho_r[6]+0.1530931089239486*uy_surf_rl[1]*rho_c[6]+0.1530931089239486*uy_surf_cr[1]*rho_c[6]+0.3061862178478971*pkpm_lax_r[1]*rho_c[6]+0.0883883476483184*uy_surf_rl[0]*rho_r[5]+0.0883883476483184*uy_surf_cr[0]*rho_r[5]-0.1767766952966368*pkpm_lax_r[0]*rho_r[5]+0.0883883476483184*uy_surf_rl[0]*rho_c[5]+0.0883883476483184*uy_surf_cr[0]*rho_c[5]+0.1767766952966368*pkpm_lax_r[0]*rho_c[5]-0.1530931089239486*uy_surf_rl[2]*rho_r[4]-0.1530931089239486*uy_surf_cr[2]*rho_r[4]+0.3061862178478971*pkpm_lax_r[2]*rho_r[4]+0.1530931089239486*uy_surf_rl[2]*rho_c[4]+0.1530931089239486*uy_surf_cr[2]*rho_c[4]+0.3061862178478971*pkpm_lax_r[2]*rho_c[4]-0.1530931089239486*rho_r[2]*uy_surf_rl[3]+0.1530931089239486*rho_c[2]*uy_surf_rl[3]+0.0883883476483184*rho_r[0]*uy_surf_rl[3]+0.0883883476483184*rho_c[0]*uy_surf_rl[3]-0.1530931089239486*rho_r[2]*uy_surf_cr[3]+0.1530931089239486*rho_c[2]*uy_surf_cr[3]+0.0883883476483184*rho_r[0]*uy_surf_cr[3]+0.0883883476483184*rho_c[0]*uy_surf_cr[3]+0.0883883476483184*uy_surf_rl[1]*rho_r[3]+0.0883883476483184*uy_surf_cr[1]*rho_r[3]-0.1767766952966368*pkpm_lax_r[1]*rho_r[3]+0.0883883476483184*uy_surf_rl[1]*rho_c[3]+0.0883883476483184*uy_surf_cr[1]*rho_c[3]+0.1767766952966368*pkpm_lax_r[1]*rho_c[3]+0.3061862178478971*rho_r[2]*pkpm_lax_r[3]+0.3061862178478971*rho_c[2]*pkpm_lax_r[3]-0.1767766952966368*rho_r[0]*pkpm_lax_r[3]+0.1767766952966368*rho_c[0]*pkpm_lax_r[3]+0.0883883476483184*rho_r[1]*uy_surf_rl[2]+0.0883883476483184*rho_c[1]*uy_surf_rl[2]+0.0883883476483184*rho_r[1]*uy_surf_cr[2]+0.0883883476483184*rho_c[1]*uy_surf_cr[2]-0.1767766952966368*rho_r[1]*pkpm_lax_r[2]+0.1767766952966368*rho_c[1]*pkpm_lax_r[2]; 

  avg_p_ij_x_l[0] = 0.6123724356957944*Pxy_l[2]-0.6123724356957944*Pxy_c[2]+0.3535533905932737*Pxy_l[0]+0.3535533905932737*Pxy_c[0]; 
  avg_p_ij_x_l[1] = 0.6123724356957944*Pxy_l[4]-0.6123724356957944*Pxy_c[4]+0.3535533905932737*Pxy_l[1]+0.3535533905932737*Pxy_c[1]; 
  avg_p_ij_x_l[2] = 0.6123724356957944*Pxy_l[6]-0.6123724356957944*Pxy_c[6]+0.3535533905932737*Pxy_l[3]+0.3535533905932737*Pxy_c[3]; 
  avg_p_ij_x_l[3] = 0.6123724356957944*Pxy_l[7]-0.6123724356957944*Pxy_c[7]+0.3535533905932737*Pxy_l[5]+0.3535533905932737*Pxy_c[5]; 

  avg_p_ij_x_r[0] = (-0.6123724356957944*Pxy_r[2])+0.6123724356957944*Pxy_c[2]+0.3535533905932737*Pxy_r[0]+0.3535533905932737*Pxy_c[0]; 
  avg_p_ij_x_r[1] = (-0.6123724356957944*Pxy_r[4])+0.6123724356957944*Pxy_c[4]+0.3535533905932737*Pxy_r[1]+0.3535533905932737*Pxy_c[1]; 
  avg_p_ij_x_r[2] = (-0.6123724356957944*Pxy_r[6])+0.6123724356957944*Pxy_c[6]+0.3535533905932737*Pxy_r[3]+0.3535533905932737*Pxy_c[3]; 
  avg_p_ij_x_r[3] = (-0.6123724356957944*Pxy_r[7])+0.6123724356957944*Pxy_c[7]+0.3535533905932737*Pxy_r[5]+0.3535533905932737*Pxy_c[5]; 

  avg_p_ij_y_l[0] = 0.6123724356957944*Pyy_l[2]-0.6123724356957944*Pyy_c[2]+0.3535533905932737*Pyy_l[0]+0.3535533905932737*Pyy_c[0]; 
  avg_p_ij_y_l[1] = 0.6123724356957944*Pyy_l[4]-0.6123724356957944*Pyy_c[4]+0.3535533905932737*Pyy_l[1]+0.3535533905932737*Pyy_c[1]; 
  avg_p_ij_y_l[2] = 0.6123724356957944*Pyy_l[6]-0.6123724356957944*Pyy_c[6]+0.3535533905932737*Pyy_l[3]+0.3535533905932737*Pyy_c[3]; 
  avg_p_ij_y_l[3] = 0.6123724356957944*Pyy_l[7]-0.6123724356957944*Pyy_c[7]+0.3535533905932737*Pyy_l[5]+0.3535533905932737*Pyy_c[5]; 

  avg_p_ij_y_r[0] = (-0.6123724356957944*Pyy_r[2])+0.6123724356957944*Pyy_c[2]+0.3535533905932737*Pyy_r[0]+0.3535533905932737*Pyy_c[0]; 
  avg_p_ij_y_r[1] = (-0.6123724356957944*Pyy_r[4])+0.6123724356957944*Pyy_c[4]+0.3535533905932737*Pyy_r[1]+0.3535533905932737*Pyy_c[1]; 
  avg_p_ij_y_r[2] = (-0.6123724356957944*Pyy_r[6])+0.6123724356957944*Pyy_c[6]+0.3535533905932737*Pyy_r[3]+0.3535533905932737*Pyy_c[3]; 
  avg_p_ij_y_r[3] = (-0.6123724356957944*Pyy_r[7])+0.6123724356957944*Pyy_c[7]+0.3535533905932737*Pyy_r[5]+0.3535533905932737*Pyy_c[5]; 

  avg_p_ij_z_l[0] = 0.6123724356957944*Pyz_l[2]-0.6123724356957944*Pyz_c[2]+0.3535533905932737*Pyz_l[0]+0.3535533905932737*Pyz_c[0]; 
  avg_p_ij_z_l[1] = 0.6123724356957944*Pyz_l[4]-0.6123724356957944*Pyz_c[4]+0.3535533905932737*Pyz_l[1]+0.3535533905932737*Pyz_c[1]; 
  avg_p_ij_z_l[2] = 0.6123724356957944*Pyz_l[6]-0.6123724356957944*Pyz_c[6]+0.3535533905932737*Pyz_l[3]+0.3535533905932737*Pyz_c[3]; 
  avg_p_ij_z_l[3] = 0.6123724356957944*Pyz_l[7]-0.6123724356957944*Pyz_c[7]+0.3535533905932737*Pyz_l[5]+0.3535533905932737*Pyz_c[5]; 

  avg_p_ij_z_r[0] = (-0.6123724356957944*Pyz_r[2])+0.6123724356957944*Pyz_c[2]+0.3535533905932737*Pyz_r[0]+0.3535533905932737*Pyz_c[0]; 
  avg_p_ij_z_r[1] = (-0.6123724356957944*Pyz_r[4])+0.6123724356957944*Pyz_c[4]+0.3535533905932737*Pyz_r[1]+0.3535533905932737*Pyz_c[1]; 
  avg_p_ij_z_r[2] = (-0.6123724356957944*Pyz_r[6])+0.6123724356957944*Pyz_c[6]+0.3535533905932737*Pyz_r[3]+0.3535533905932737*Pyz_c[3]; 
  avg_p_ij_z_r[3] = (-0.6123724356957944*Pyz_r[7])+0.6123724356957944*Pyz_c[7]+0.3535533905932737*Pyz_r[5]+0.3535533905932737*Pyz_c[5]; 

  double amdq_rhoux_l[4] = {0.0}; 
  double apdq_rhoux_l[4] = {0.0}; 
  double amdq_rhouy_l[4] = {0.0}; 
  double apdq_rhouy_l[4] = {0.0}; 
  double amdq_rhouz_l[4] = {0.0}; 
  double apdq_rhouz_l[4] = {0.0}; 

  double amdq_rhoux_r[4] = {0.0}; 
  double apdq_rhoux_r[4] = {0.0}; 
  double amdq_rhouy_r[4] = {0.0}; 
  double apdq_rhouy_r[4] = {0.0}; 
  double amdq_rhouz_r[4] = {0.0}; 
  double apdq_rhouz_r[4] = {0.0}; 

  double amdq_rhoux_quad_l[4] = {0.0}; 
  double apdq_rhoux_quad_l[4] = {0.0}; 
  double amdq_rhouy_quad_l[4] = {0.0}; 
  double apdq_rhouy_quad_l[4] = {0.0}; 
  double amdq_rhouz_quad_l[4] = {0.0}; 
  double apdq_rhouz_quad_l[4] = {0.0}; 

  double amdq_rhoux_quad_r[4] = {0.0}; 
  double apdq_rhoux_quad_r[4] = {0.0}; 
  double amdq_rhouy_quad_r[4] = {0.0}; 
  double apdq_rhouy_quad_r[4] = {0.0}; 
  double amdq_rhouz_quad_r[4] = {0.0}; 
  double apdq_rhouz_quad_r[4] = {0.0}; 

  double q_lr[10] = {0.0}; 
  double q_cl[10] = {0.0}; 
  double q_cr[10] = {0.0}; 
  double q_rl[10] = {0.0}; 
  double q_lr_local[10] = {0.0}; 
  double q_cl_local[10] = {0.0}; 
  double q_cr_local[10] = {0.0}; 
  double q_rl_local[10] = {0.0}; 
  double delta_l[10] = {0.0}; 
  double delta_r[10] = {0.0}; 
  double my_max_speed_l = 0.0; 
  double my_max_speed_r = 0.0; 
  double lenr_l = 0.0; 
  double lenr_r = 0.0; 
  double waves_l[50] = {0.0}; 
  double waves_r[50] = {0.0}; 
  double speeds_l[5] = {0.0}; 
  double speeds_r[5] = {0.0}; 
  double amdq_l_local[10] = {0.0}; 
  double apdq_l_local[10] = {0.0}; 
  double amdq_r_local[10] = {0.0}; 
  double apdq_r_local[10] = {0.0}; 
  double amdq_l[10] = {0.0}; 
  double apdq_l[10] = {0.0}; 
  double amdq_r[10] = {0.0}; 
  double apdq_r[10] = {0.0}; 

  q_lr[0] = ser_3x_p1_surfx2_eval_quad_node_0_r(rho_l); 
  q_lr[1] = q_lr[0]*(0.5*ux_surf_lr[3]-0.5*(ux_surf_lr[2]+ux_surf_lr[1])+0.5*ux_surf_lr[0]); 
  q_lr[2] = q_lr[0]*(0.5*uy_surf_lr[3]-0.5*(uy_surf_lr[2]+uy_surf_lr[1])+0.5*uy_surf_lr[0]); 
  q_lr[3] = q_lr[0]*(0.5*uz_surf_lr[3]-0.5*(uz_surf_lr[2]+uz_surf_lr[1])+0.5*uz_surf_lr[0]); 
  q_lr[4] = ser_3x_p1_surfx2_eval_quad_node_0_r(Pxx_l) + q_lr[1]*q_lr[1]/q_lr[0]; 
  q_lr[5] = ser_3x_p1_surfx2_eval_quad_node_0_r(Pxy_l) + q_lr[1]*q_lr[2]/q_lr[0]; 
  q_lr[6] = ser_3x_p1_surfx2_eval_quad_node_0_r(Pxz_l) + q_lr[1]*q_lr[3]/q_lr[0]; 
  q_lr[7] = ser_3x_p1_surfx2_eval_quad_node_0_r(Pyy_l) + q_lr[2]*q_lr[2]/q_lr[0]; 
  q_lr[8] = ser_3x_p1_surfx2_eval_quad_node_0_r(Pyz_l) + q_lr[2]*q_lr[3]/q_lr[0]; 
  q_lr[9] = ser_3x_p1_surfx2_eval_quad_node_0_r(Pzz_l) + q_lr[3]*q_lr[3]/q_lr[0]; 
  q_cl[0] = ser_3x_p1_surfx2_eval_quad_node_0_l(rho_c); 
  q_cl[1] = q_cl[0]*(0.5*ux_surf_cl[3]-0.5*(ux_surf_cl[2]+ux_surf_cl[1])+0.5*ux_surf_cl[0]); 
  q_cl[2] = q_cl[0]*(0.5*uy_surf_cl[3]-0.5*(uy_surf_cl[2]+uy_surf_cl[1])+0.5*uy_surf_cl[0]); 
  q_cl[3] = q_cl[0]*(0.5*uz_surf_cl[3]-0.5*(uz_surf_cl[2]+uz_surf_cl[1])+0.5*uz_surf_cl[0]); 
  q_cl[4] = ser_3x_p1_surfx2_eval_quad_node_0_l(Pxx_c) + q_cl[1]*q_cl[1]/q_cl[0]; 
  q_cl[5] = ser_3x_p1_surfx2_eval_quad_node_0_l(Pxy_c) + q_cl[1]*q_cl[2]/q_cl[0]; 
  q_cl[6] = ser_3x_p1_surfx2_eval_quad_node_0_l(Pxz_c) + q_cl[1]*q_cl[3]/q_cl[0]; 
  q_cl[7] = ser_3x_p1_surfx2_eval_quad_node_0_l(Pyy_c) + q_cl[2]*q_cl[2]/q_cl[0]; 
  q_cl[8] = ser_3x_p1_surfx2_eval_quad_node_0_l(Pyz_c) + q_cl[2]*q_cl[3]/q_cl[0]; 
  q_cl[9] = ser_3x_p1_surfx2_eval_quad_node_0_l(Pzz_c) + q_cl[3]*q_cl[3]/q_cl[0]; 
  q_cr[0] = ser_3x_p1_surfx2_eval_quad_node_0_r(rho_c); 
  q_cr[1] = q_cr[0]*(0.5*ux_surf_cr[3]-0.5*(ux_surf_cr[2]+ux_surf_cr[1])+0.5*ux_surf_cr[0]); 
  q_cr[2] = q_cr[0]*(0.5*uy_surf_cr[3]-0.5*(uy_surf_cr[2]+uy_surf_cr[1])+0.5*uy_surf_cr[0]); 
  q_cr[3] = q_cr[0]*(0.5*uz_surf_cr[3]-0.5*(uz_surf_cr[2]+uz_surf_cr[1])+0.5*uz_surf_cr[0]); 
  q_cr[4] = ser_3x_p1_surfx2_eval_quad_node_0_r(Pxx_c) + q_cr[1]*q_cr[1]/q_cr[0]; 
  q_cr[5] = ser_3x_p1_surfx2_eval_quad_node_0_r(Pxy_c) + q_cr[1]*q_cr[2]/q_cr[0]; 
  q_cr[6] = ser_3x_p1_surfx2_eval_quad_node_0_r(Pxz_c) + q_cr[1]*q_cr[3]/q_cr[0]; 
  q_cr[7] = ser_3x_p1_surfx2_eval_quad_node_0_r(Pyy_c) + q_cr[2]*q_cr[2]/q_cr[0]; 
  q_cr[8] = ser_3x_p1_surfx2_eval_quad_node_0_r(Pyz_c) + q_cr[2]*q_cr[3]/q_cr[0]; 
  q_cr[9] = ser_3x_p1_surfx2_eval_quad_node_0_r(Pzz_c) + q_cr[3]*q_cr[3]/q_cr[0]; 
  q_rl[0] = ser_3x_p1_surfx2_eval_quad_node_0_l(rho_r); 
  q_rl[1] = q_rl[0]*(0.5*ux_surf_rl[3]-0.5*(ux_surf_rl[2]+ux_surf_rl[1])+0.5*ux_surf_rl[0]); 
  q_rl[2] = q_rl[0]*(0.5*uy_surf_rl[3]-0.5*(uy_surf_rl[2]+uy_surf_rl[1])+0.5*uy_surf_rl[0]); 
  q_rl[3] = q_rl[0]*(0.5*uz_surf_rl[3]-0.5*(uz_surf_rl[2]+uz_surf_rl[1])+0.5*uz_surf_rl[0]); 
  q_rl[4] = ser_3x_p1_surfx2_eval_quad_node_0_l(Pxx_r) + q_rl[1]*q_rl[1]/q_rl[0]; 
  q_rl[5] = ser_3x_p1_surfx2_eval_quad_node_0_l(Pxy_r) + q_rl[1]*q_rl[2]/q_rl[0]; 
  q_rl[6] = ser_3x_p1_surfx2_eval_quad_node_0_l(Pxz_r) + q_rl[1]*q_rl[3]/q_rl[0]; 
  q_rl[7] = ser_3x_p1_surfx2_eval_quad_node_0_l(Pyy_r) + q_rl[2]*q_rl[2]/q_rl[0]; 
  q_rl[8] = ser_3x_p1_surfx2_eval_quad_node_0_l(Pyz_r) + q_rl[2]*q_rl[3]/q_rl[0]; 
  q_rl[9] = ser_3x_p1_surfx2_eval_quad_node_0_l(Pzz_r) + q_rl[3]*q_rl[3]/q_rl[0]; 

  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_l->tau1[1], geom_l->tau2[1], geom_l->norm[1], q_lr, q_lr_local); 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_l->tau1[1], geom_l->tau2[1], geom_l->norm[1], q_cl, q_cl_local); 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_r->tau1[1], geom_r->tau2[1], geom_r->norm[1], q_cr, q_cr_local); 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_r->tau1[1], geom_r->tau2[1], geom_r->norm[1], q_rl, q_rl_local); 

  delta_l[0] = q_cl_local[0] - q_lr_local[0]; 
  delta_l[1] = q_cl_local[1] - q_lr_local[1]; 
  delta_l[2] = q_cl_local[2] - q_lr_local[2]; 
  delta_l[3] = q_cl_local[3] - q_lr_local[3]; 
  delta_l[4] = q_cl_local[4] - q_lr_local[4]; 
  delta_l[5] = q_cl_local[5] - q_lr_local[5]; 
  delta_l[6] = q_cl_local[6] - q_lr_local[6]; 
  delta_l[7] = q_cl_local[7] - q_lr_local[7]; 
  delta_l[8] = q_cl_local[8] - q_lr_local[8]; 
  delta_l[9] = q_cl_local[9] - q_lr_local[9]; 
  delta_r[0] = q_rl_local[0] - q_cr_local[0]; 
  delta_r[1] = q_rl_local[1] - q_cr_local[1]; 
  delta_r[2] = q_rl_local[2] - q_cr_local[2]; 
  delta_r[3] = q_rl_local[3] - q_cr_local[3]; 
  delta_r[4] = q_rl_local[4] - q_cr_local[4]; 
  delta_r[5] = q_rl_local[5] - q_cr_local[5]; 
  delta_r[6] = q_rl_local[6] - q_cr_local[6]; 
  delta_r[7] = q_rl_local[7] - q_cr_local[7]; 
  delta_r[8] = q_rl_local[8] - q_cr_local[8]; 
  delta_r[9] = q_rl_local[9] - q_cr_local[9]; 

  my_max_speed_l = gkyl_wv_eqn_waves(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta_l, q_lr_local, q_cl_local, waves_l, speeds_l); 
  my_max_speed_r = gkyl_wv_eqn_waves(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta_r, q_cr_local, q_rl_local, waves_r, speeds_r); 
  lenr_l = geom_l->lenr[1]; 
  speeds_l[0] *= lenr_l; 
  speeds_l[1] *= lenr_l; 
  speeds_l[2] *= lenr_l; 
  speeds_l[3] *= lenr_l; 
  speeds_l[4] *= lenr_l; 
  lenr_r = geom_r->lenr[1]; 
  speeds_r[0] *= lenr_r; 
  speeds_r[1] *= lenr_r; 
  speeds_r[2] *= lenr_r; 
  speeds_r[3] *= lenr_r; 
  speeds_r[4] *= lenr_r; 

  gkyl_wv_eqn_qfluct(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, q_lr_local, q_cl_local, waves_l, speeds_l, amdq_l_local, apdq_l_local); 
  gkyl_wv_eqn_qfluct(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, q_cr_local, q_rl_local, waves_r, speeds_r, amdq_r_local, apdq_r_local); 

  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_l->tau1[1], geom_l->tau2[1], geom_l->norm[1], amdq_l_local, amdq_l); 
  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_l->tau1[1], geom_l->tau2[1], geom_l->norm[1], apdq_l_local, apdq_l); 
  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_r->tau1[1], geom_r->tau2[1], geom_r->norm[1], amdq_r_local, amdq_r); 
  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_r->tau1[1], geom_r->tau2[1], geom_r->norm[1], apdq_r_local, apdq_r); 

  amdq_rhoux_quad_l[0] = amdq_l[1]; 
  apdq_rhoux_quad_l[0] = apdq_l[1]; 
  amdq_rhouy_quad_l[0] = amdq_l[2]; 
  apdq_rhouy_quad_l[0] = apdq_l[2]; 
  amdq_rhouz_quad_l[0] = amdq_l[3]; 
  apdq_rhouz_quad_l[0] = apdq_l[3]; 

  amdq_rhoux_quad_r[0] = amdq_r[1]; 
  apdq_rhoux_quad_r[0] = apdq_r[1]; 
  amdq_rhouy_quad_r[0] = amdq_r[2]; 
  apdq_rhouy_quad_r[0] = apdq_r[2]; 
  amdq_rhouz_quad_r[0] = amdq_r[3]; 
  apdq_rhouz_quad_r[0] = apdq_r[3]; 

  q_lr[0] = ser_3x_p1_surfx2_eval_quad_node_1_r(rho_l); 
  q_lr[1] = q_lr[0]*((-0.5*ux_surf_lr[3])+0.5*ux_surf_lr[2]-0.5*ux_surf_lr[1]+0.5*ux_surf_lr[0]); 
  q_lr[2] = q_lr[0]*((-0.5*uy_surf_lr[3])+0.5*uy_surf_lr[2]-0.5*uy_surf_lr[1]+0.5*uy_surf_lr[0]); 
  q_lr[3] = q_lr[0]*((-0.5*uz_surf_lr[3])+0.5*uz_surf_lr[2]-0.5*uz_surf_lr[1]+0.5*uz_surf_lr[0]); 
  q_lr[4] = ser_3x_p1_surfx2_eval_quad_node_1_r(Pxx_l) + q_lr[1]*q_lr[1]/q_lr[0]; 
  q_lr[5] = ser_3x_p1_surfx2_eval_quad_node_1_r(Pxy_l) + q_lr[1]*q_lr[2]/q_lr[0]; 
  q_lr[6] = ser_3x_p1_surfx2_eval_quad_node_1_r(Pxz_l) + q_lr[1]*q_lr[3]/q_lr[0]; 
  q_lr[7] = ser_3x_p1_surfx2_eval_quad_node_1_r(Pyy_l) + q_lr[2]*q_lr[2]/q_lr[0]; 
  q_lr[8] = ser_3x_p1_surfx2_eval_quad_node_1_r(Pyz_l) + q_lr[2]*q_lr[3]/q_lr[0]; 
  q_lr[9] = ser_3x_p1_surfx2_eval_quad_node_1_r(Pzz_l) + q_lr[3]*q_lr[3]/q_lr[0]; 
  q_cl[0] = ser_3x_p1_surfx2_eval_quad_node_1_l(rho_c); 
  q_cl[1] = q_cl[0]*((-0.5*ux_surf_cl[3])+0.5*ux_surf_cl[2]-0.5*ux_surf_cl[1]+0.5*ux_surf_cl[0]); 
  q_cl[2] = q_cl[0]*((-0.5*uy_surf_cl[3])+0.5*uy_surf_cl[2]-0.5*uy_surf_cl[1]+0.5*uy_surf_cl[0]); 
  q_cl[3] = q_cl[0]*((-0.5*uz_surf_cl[3])+0.5*uz_surf_cl[2]-0.5*uz_surf_cl[1]+0.5*uz_surf_cl[0]); 
  q_cl[4] = ser_3x_p1_surfx2_eval_quad_node_1_l(Pxx_c) + q_cl[1]*q_cl[1]/q_cl[0]; 
  q_cl[5] = ser_3x_p1_surfx2_eval_quad_node_1_l(Pxy_c) + q_cl[1]*q_cl[2]/q_cl[0]; 
  q_cl[6] = ser_3x_p1_surfx2_eval_quad_node_1_l(Pxz_c) + q_cl[1]*q_cl[3]/q_cl[0]; 
  q_cl[7] = ser_3x_p1_surfx2_eval_quad_node_1_l(Pyy_c) + q_cl[2]*q_cl[2]/q_cl[0]; 
  q_cl[8] = ser_3x_p1_surfx2_eval_quad_node_1_l(Pyz_c) + q_cl[2]*q_cl[3]/q_cl[0]; 
  q_cl[9] = ser_3x_p1_surfx2_eval_quad_node_1_l(Pzz_c) + q_cl[3]*q_cl[3]/q_cl[0]; 
  q_cr[0] = ser_3x_p1_surfx2_eval_quad_node_1_r(rho_c); 
  q_cr[1] = q_cr[0]*((-0.5*ux_surf_cr[3])+0.5*ux_surf_cr[2]-0.5*ux_surf_cr[1]+0.5*ux_surf_cr[0]); 
  q_cr[2] = q_cr[0]*((-0.5*uy_surf_cr[3])+0.5*uy_surf_cr[2]-0.5*uy_surf_cr[1]+0.5*uy_surf_cr[0]); 
  q_cr[3] = q_cr[0]*((-0.5*uz_surf_cr[3])+0.5*uz_surf_cr[2]-0.5*uz_surf_cr[1]+0.5*uz_surf_cr[0]); 
  q_cr[4] = ser_3x_p1_surfx2_eval_quad_node_1_r(Pxx_c) + q_cr[1]*q_cr[1]/q_cr[0]; 
  q_cr[5] = ser_3x_p1_surfx2_eval_quad_node_1_r(Pxy_c) + q_cr[1]*q_cr[2]/q_cr[0]; 
  q_cr[6] = ser_3x_p1_surfx2_eval_quad_node_1_r(Pxz_c) + q_cr[1]*q_cr[3]/q_cr[0]; 
  q_cr[7] = ser_3x_p1_surfx2_eval_quad_node_1_r(Pyy_c) + q_cr[2]*q_cr[2]/q_cr[0]; 
  q_cr[8] = ser_3x_p1_surfx2_eval_quad_node_1_r(Pyz_c) + q_cr[2]*q_cr[3]/q_cr[0]; 
  q_cr[9] = ser_3x_p1_surfx2_eval_quad_node_1_r(Pzz_c) + q_cr[3]*q_cr[3]/q_cr[0]; 
  q_rl[0] = ser_3x_p1_surfx2_eval_quad_node_1_l(rho_r); 
  q_rl[1] = q_rl[0]*((-0.5*ux_surf_rl[3])+0.5*ux_surf_rl[2]-0.5*ux_surf_rl[1]+0.5*ux_surf_rl[0]); 
  q_rl[2] = q_rl[0]*((-0.5*uy_surf_rl[3])+0.5*uy_surf_rl[2]-0.5*uy_surf_rl[1]+0.5*uy_surf_rl[0]); 
  q_rl[3] = q_rl[0]*((-0.5*uz_surf_rl[3])+0.5*uz_surf_rl[2]-0.5*uz_surf_rl[1]+0.5*uz_surf_rl[0]); 
  q_rl[4] = ser_3x_p1_surfx2_eval_quad_node_1_l(Pxx_r) + q_rl[1]*q_rl[1]/q_rl[0]; 
  q_rl[5] = ser_3x_p1_surfx2_eval_quad_node_1_l(Pxy_r) + q_rl[1]*q_rl[2]/q_rl[0]; 
  q_rl[6] = ser_3x_p1_surfx2_eval_quad_node_1_l(Pxz_r) + q_rl[1]*q_rl[3]/q_rl[0]; 
  q_rl[7] = ser_3x_p1_surfx2_eval_quad_node_1_l(Pyy_r) + q_rl[2]*q_rl[2]/q_rl[0]; 
  q_rl[8] = ser_3x_p1_surfx2_eval_quad_node_1_l(Pyz_r) + q_rl[2]*q_rl[3]/q_rl[0]; 
  q_rl[9] = ser_3x_p1_surfx2_eval_quad_node_1_l(Pzz_r) + q_rl[3]*q_rl[3]/q_rl[0]; 

  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_l->tau1[1], geom_l->tau2[1], geom_l->norm[1], q_lr, q_lr_local); 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_l->tau1[1], geom_l->tau2[1], geom_l->norm[1], q_cl, q_cl_local); 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_r->tau1[1], geom_r->tau2[1], geom_r->norm[1], q_cr, q_cr_local); 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_r->tau1[1], geom_r->tau2[1], geom_r->norm[1], q_rl, q_rl_local); 

  delta_l[0] = q_cl_local[0] - q_lr_local[0]; 
  delta_l[1] = q_cl_local[1] - q_lr_local[1]; 
  delta_l[2] = q_cl_local[2] - q_lr_local[2]; 
  delta_l[3] = q_cl_local[3] - q_lr_local[3]; 
  delta_l[4] = q_cl_local[4] - q_lr_local[4]; 
  delta_l[5] = q_cl_local[5] - q_lr_local[5]; 
  delta_l[6] = q_cl_local[6] - q_lr_local[6]; 
  delta_l[7] = q_cl_local[7] - q_lr_local[7]; 
  delta_l[8] = q_cl_local[8] - q_lr_local[8]; 
  delta_l[9] = q_cl_local[9] - q_lr_local[9]; 
  delta_r[0] = q_rl_local[0] - q_cr_local[0]; 
  delta_r[1] = q_rl_local[1] - q_cr_local[1]; 
  delta_r[2] = q_rl_local[2] - q_cr_local[2]; 
  delta_r[3] = q_rl_local[3] - q_cr_local[3]; 
  delta_r[4] = q_rl_local[4] - q_cr_local[4]; 
  delta_r[5] = q_rl_local[5] - q_cr_local[5]; 
  delta_r[6] = q_rl_local[6] - q_cr_local[6]; 
  delta_r[7] = q_rl_local[7] - q_cr_local[7]; 
  delta_r[8] = q_rl_local[8] - q_cr_local[8]; 
  delta_r[9] = q_rl_local[9] - q_cr_local[9]; 

  my_max_speed_l = gkyl_wv_eqn_waves(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta_l, q_lr_local, q_cl_local, waves_l, speeds_l); 
  my_max_speed_r = gkyl_wv_eqn_waves(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta_r, q_cr_local, q_rl_local, waves_r, speeds_r); 
  lenr_l = geom_l->lenr[1]; 
  speeds_l[0] *= lenr_l; 
  speeds_l[1] *= lenr_l; 
  speeds_l[2] *= lenr_l; 
  speeds_l[3] *= lenr_l; 
  speeds_l[4] *= lenr_l; 
  lenr_r = geom_r->lenr[1]; 
  speeds_r[0] *= lenr_r; 
  speeds_r[1] *= lenr_r; 
  speeds_r[2] *= lenr_r; 
  speeds_r[3] *= lenr_r; 
  speeds_r[4] *= lenr_r; 

  gkyl_wv_eqn_qfluct(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, q_lr_local, q_cl_local, waves_l, speeds_l, amdq_l_local, apdq_l_local); 
  gkyl_wv_eqn_qfluct(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, q_cr_local, q_rl_local, waves_r, speeds_r, amdq_r_local, apdq_r_local); 

  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_l->tau1[1], geom_l->tau2[1], geom_l->norm[1], amdq_l_local, amdq_l); 
  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_l->tau1[1], geom_l->tau2[1], geom_l->norm[1], apdq_l_local, apdq_l); 
  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_r->tau1[1], geom_r->tau2[1], geom_r->norm[1], amdq_r_local, amdq_r); 
  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_r->tau1[1], geom_r->tau2[1], geom_r->norm[1], apdq_r_local, apdq_r); 

  amdq_rhoux_quad_l[1] = amdq_l[1]; 
  apdq_rhoux_quad_l[1] = apdq_l[1]; 
  amdq_rhouy_quad_l[1] = amdq_l[2]; 
  apdq_rhouy_quad_l[1] = apdq_l[2]; 
  amdq_rhouz_quad_l[1] = amdq_l[3]; 
  apdq_rhouz_quad_l[1] = apdq_l[3]; 

  amdq_rhoux_quad_r[1] = amdq_r[1]; 
  apdq_rhoux_quad_r[1] = apdq_r[1]; 
  amdq_rhouy_quad_r[1] = amdq_r[2]; 
  apdq_rhouy_quad_r[1] = apdq_r[2]; 
  amdq_rhouz_quad_r[1] = amdq_r[3]; 
  apdq_rhouz_quad_r[1] = apdq_r[3]; 

  q_lr[0] = ser_3x_p1_surfx2_eval_quad_node_2_r(rho_l); 
  q_lr[1] = q_lr[0]*(0.5*(ux_surf_lr[1]+ux_surf_lr[0])-0.5*(ux_surf_lr[3]+ux_surf_lr[2])); 
  q_lr[2] = q_lr[0]*(0.5*(uy_surf_lr[1]+uy_surf_lr[0])-0.5*(uy_surf_lr[3]+uy_surf_lr[2])); 
  q_lr[3] = q_lr[0]*(0.5*(uz_surf_lr[1]+uz_surf_lr[0])-0.5*(uz_surf_lr[3]+uz_surf_lr[2])); 
  q_lr[4] = ser_3x_p1_surfx2_eval_quad_node_2_r(Pxx_l) + q_lr[1]*q_lr[1]/q_lr[0]; 
  q_lr[5] = ser_3x_p1_surfx2_eval_quad_node_2_r(Pxy_l) + q_lr[1]*q_lr[2]/q_lr[0]; 
  q_lr[6] = ser_3x_p1_surfx2_eval_quad_node_2_r(Pxz_l) + q_lr[1]*q_lr[3]/q_lr[0]; 
  q_lr[7] = ser_3x_p1_surfx2_eval_quad_node_2_r(Pyy_l) + q_lr[2]*q_lr[2]/q_lr[0]; 
  q_lr[8] = ser_3x_p1_surfx2_eval_quad_node_2_r(Pyz_l) + q_lr[2]*q_lr[3]/q_lr[0]; 
  q_lr[9] = ser_3x_p1_surfx2_eval_quad_node_2_r(Pzz_l) + q_lr[3]*q_lr[3]/q_lr[0]; 
  q_cl[0] = ser_3x_p1_surfx2_eval_quad_node_2_l(rho_c); 
  q_cl[1] = q_cl[0]*(0.5*(ux_surf_cl[1]+ux_surf_cl[0])-0.5*(ux_surf_cl[3]+ux_surf_cl[2])); 
  q_cl[2] = q_cl[0]*(0.5*(uy_surf_cl[1]+uy_surf_cl[0])-0.5*(uy_surf_cl[3]+uy_surf_cl[2])); 
  q_cl[3] = q_cl[0]*(0.5*(uz_surf_cl[1]+uz_surf_cl[0])-0.5*(uz_surf_cl[3]+uz_surf_cl[2])); 
  q_cl[4] = ser_3x_p1_surfx2_eval_quad_node_2_l(Pxx_c) + q_cl[1]*q_cl[1]/q_cl[0]; 
  q_cl[5] = ser_3x_p1_surfx2_eval_quad_node_2_l(Pxy_c) + q_cl[1]*q_cl[2]/q_cl[0]; 
  q_cl[6] = ser_3x_p1_surfx2_eval_quad_node_2_l(Pxz_c) + q_cl[1]*q_cl[3]/q_cl[0]; 
  q_cl[7] = ser_3x_p1_surfx2_eval_quad_node_2_l(Pyy_c) + q_cl[2]*q_cl[2]/q_cl[0]; 
  q_cl[8] = ser_3x_p1_surfx2_eval_quad_node_2_l(Pyz_c) + q_cl[2]*q_cl[3]/q_cl[0]; 
  q_cl[9] = ser_3x_p1_surfx2_eval_quad_node_2_l(Pzz_c) + q_cl[3]*q_cl[3]/q_cl[0]; 
  q_cr[0] = ser_3x_p1_surfx2_eval_quad_node_2_r(rho_c); 
  q_cr[1] = q_cr[0]*(0.5*(ux_surf_cr[1]+ux_surf_cr[0])-0.5*(ux_surf_cr[3]+ux_surf_cr[2])); 
  q_cr[2] = q_cr[0]*(0.5*(uy_surf_cr[1]+uy_surf_cr[0])-0.5*(uy_surf_cr[3]+uy_surf_cr[2])); 
  q_cr[3] = q_cr[0]*(0.5*(uz_surf_cr[1]+uz_surf_cr[0])-0.5*(uz_surf_cr[3]+uz_surf_cr[2])); 
  q_cr[4] = ser_3x_p1_surfx2_eval_quad_node_2_r(Pxx_c) + q_cr[1]*q_cr[1]/q_cr[0]; 
  q_cr[5] = ser_3x_p1_surfx2_eval_quad_node_2_r(Pxy_c) + q_cr[1]*q_cr[2]/q_cr[0]; 
  q_cr[6] = ser_3x_p1_surfx2_eval_quad_node_2_r(Pxz_c) + q_cr[1]*q_cr[3]/q_cr[0]; 
  q_cr[7] = ser_3x_p1_surfx2_eval_quad_node_2_r(Pyy_c) + q_cr[2]*q_cr[2]/q_cr[0]; 
  q_cr[8] = ser_3x_p1_surfx2_eval_quad_node_2_r(Pyz_c) + q_cr[2]*q_cr[3]/q_cr[0]; 
  q_cr[9] = ser_3x_p1_surfx2_eval_quad_node_2_r(Pzz_c) + q_cr[3]*q_cr[3]/q_cr[0]; 
  q_rl[0] = ser_3x_p1_surfx2_eval_quad_node_2_l(rho_r); 
  q_rl[1] = q_rl[0]*(0.5*(ux_surf_rl[1]+ux_surf_rl[0])-0.5*(ux_surf_rl[3]+ux_surf_rl[2])); 
  q_rl[2] = q_rl[0]*(0.5*(uy_surf_rl[1]+uy_surf_rl[0])-0.5*(uy_surf_rl[3]+uy_surf_rl[2])); 
  q_rl[3] = q_rl[0]*(0.5*(uz_surf_rl[1]+uz_surf_rl[0])-0.5*(uz_surf_rl[3]+uz_surf_rl[2])); 
  q_rl[4] = ser_3x_p1_surfx2_eval_quad_node_2_l(Pxx_r) + q_rl[1]*q_rl[1]/q_rl[0]; 
  q_rl[5] = ser_3x_p1_surfx2_eval_quad_node_2_l(Pxy_r) + q_rl[1]*q_rl[2]/q_rl[0]; 
  q_rl[6] = ser_3x_p1_surfx2_eval_quad_node_2_l(Pxz_r) + q_rl[1]*q_rl[3]/q_rl[0]; 
  q_rl[7] = ser_3x_p1_surfx2_eval_quad_node_2_l(Pyy_r) + q_rl[2]*q_rl[2]/q_rl[0]; 
  q_rl[8] = ser_3x_p1_surfx2_eval_quad_node_2_l(Pyz_r) + q_rl[2]*q_rl[3]/q_rl[0]; 
  q_rl[9] = ser_3x_p1_surfx2_eval_quad_node_2_l(Pzz_r) + q_rl[3]*q_rl[3]/q_rl[0]; 

  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_l->tau1[1], geom_l->tau2[1], geom_l->norm[1], q_lr, q_lr_local); 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_l->tau1[1], geom_l->tau2[1], geom_l->norm[1], q_cl, q_cl_local); 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_r->tau1[1], geom_r->tau2[1], geom_r->norm[1], q_cr, q_cr_local); 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_r->tau1[1], geom_r->tau2[1], geom_r->norm[1], q_rl, q_rl_local); 

  delta_l[0] = q_cl_local[0] - q_lr_local[0]; 
  delta_l[1] = q_cl_local[1] - q_lr_local[1]; 
  delta_l[2] = q_cl_local[2] - q_lr_local[2]; 
  delta_l[3] = q_cl_local[3] - q_lr_local[3]; 
  delta_l[4] = q_cl_local[4] - q_lr_local[4]; 
  delta_l[5] = q_cl_local[5] - q_lr_local[5]; 
  delta_l[6] = q_cl_local[6] - q_lr_local[6]; 
  delta_l[7] = q_cl_local[7] - q_lr_local[7]; 
  delta_l[8] = q_cl_local[8] - q_lr_local[8]; 
  delta_l[9] = q_cl_local[9] - q_lr_local[9]; 
  delta_r[0] = q_rl_local[0] - q_cr_local[0]; 
  delta_r[1] = q_rl_local[1] - q_cr_local[1]; 
  delta_r[2] = q_rl_local[2] - q_cr_local[2]; 
  delta_r[3] = q_rl_local[3] - q_cr_local[3]; 
  delta_r[4] = q_rl_local[4] - q_cr_local[4]; 
  delta_r[5] = q_rl_local[5] - q_cr_local[5]; 
  delta_r[6] = q_rl_local[6] - q_cr_local[6]; 
  delta_r[7] = q_rl_local[7] - q_cr_local[7]; 
  delta_r[8] = q_rl_local[8] - q_cr_local[8]; 
  delta_r[9] = q_rl_local[9] - q_cr_local[9]; 

  my_max_speed_l = gkyl_wv_eqn_waves(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta_l, q_lr_local, q_cl_local, waves_l, speeds_l); 
  my_max_speed_r = gkyl_wv_eqn_waves(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta_r, q_cr_local, q_rl_local, waves_r, speeds_r); 
  lenr_l = geom_l->lenr[1]; 
  speeds_l[0] *= lenr_l; 
  speeds_l[1] *= lenr_l; 
  speeds_l[2] *= lenr_l; 
  speeds_l[3] *= lenr_l; 
  speeds_l[4] *= lenr_l; 
  lenr_r = geom_r->lenr[1]; 
  speeds_r[0] *= lenr_r; 
  speeds_r[1] *= lenr_r; 
  speeds_r[2] *= lenr_r; 
  speeds_r[3] *= lenr_r; 
  speeds_r[4] *= lenr_r; 

  gkyl_wv_eqn_qfluct(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, q_lr_local, q_cl_local, waves_l, speeds_l, amdq_l_local, apdq_l_local); 
  gkyl_wv_eqn_qfluct(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, q_cr_local, q_rl_local, waves_r, speeds_r, amdq_r_local, apdq_r_local); 

  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_l->tau1[1], geom_l->tau2[1], geom_l->norm[1], amdq_l_local, amdq_l); 
  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_l->tau1[1], geom_l->tau2[1], geom_l->norm[1], apdq_l_local, apdq_l); 
  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_r->tau1[1], geom_r->tau2[1], geom_r->norm[1], amdq_r_local, amdq_r); 
  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_r->tau1[1], geom_r->tau2[1], geom_r->norm[1], apdq_r_local, apdq_r); 

  amdq_rhoux_quad_l[2] = amdq_l[1]; 
  apdq_rhoux_quad_l[2] = apdq_l[1]; 
  amdq_rhouy_quad_l[2] = amdq_l[2]; 
  apdq_rhouy_quad_l[2] = apdq_l[2]; 
  amdq_rhouz_quad_l[2] = amdq_l[3]; 
  apdq_rhouz_quad_l[2] = apdq_l[3]; 

  amdq_rhoux_quad_r[2] = amdq_r[1]; 
  apdq_rhoux_quad_r[2] = apdq_r[1]; 
  amdq_rhouy_quad_r[2] = amdq_r[2]; 
  apdq_rhouy_quad_r[2] = apdq_r[2]; 
  amdq_rhouz_quad_r[2] = amdq_r[3]; 
  apdq_rhouz_quad_r[2] = apdq_r[3]; 

  q_lr[0] = ser_3x_p1_surfx2_eval_quad_node_3_r(rho_l); 
  q_lr[1] = q_lr[0]*(0.5*(ux_surf_lr[3]+ux_surf_lr[2]+ux_surf_lr[1]+ux_surf_lr[0])); 
  q_lr[2] = q_lr[0]*(0.5*(uy_surf_lr[3]+uy_surf_lr[2]+uy_surf_lr[1]+uy_surf_lr[0])); 
  q_lr[3] = q_lr[0]*(0.5*(uz_surf_lr[3]+uz_surf_lr[2]+uz_surf_lr[1]+uz_surf_lr[0])); 
  q_lr[4] = ser_3x_p1_surfx2_eval_quad_node_3_r(Pxx_l) + q_lr[1]*q_lr[1]/q_lr[0]; 
  q_lr[5] = ser_3x_p1_surfx2_eval_quad_node_3_r(Pxy_l) + q_lr[1]*q_lr[2]/q_lr[0]; 
  q_lr[6] = ser_3x_p1_surfx2_eval_quad_node_3_r(Pxz_l) + q_lr[1]*q_lr[3]/q_lr[0]; 
  q_lr[7] = ser_3x_p1_surfx2_eval_quad_node_3_r(Pyy_l) + q_lr[2]*q_lr[2]/q_lr[0]; 
  q_lr[8] = ser_3x_p1_surfx2_eval_quad_node_3_r(Pyz_l) + q_lr[2]*q_lr[3]/q_lr[0]; 
  q_lr[9] = ser_3x_p1_surfx2_eval_quad_node_3_r(Pzz_l) + q_lr[3]*q_lr[3]/q_lr[0]; 
  q_cl[0] = ser_3x_p1_surfx2_eval_quad_node_3_l(rho_c); 
  q_cl[1] = q_cl[0]*(0.5*(ux_surf_cl[3]+ux_surf_cl[2]+ux_surf_cl[1]+ux_surf_cl[0])); 
  q_cl[2] = q_cl[0]*(0.5*(uy_surf_cl[3]+uy_surf_cl[2]+uy_surf_cl[1]+uy_surf_cl[0])); 
  q_cl[3] = q_cl[0]*(0.5*(uz_surf_cl[3]+uz_surf_cl[2]+uz_surf_cl[1]+uz_surf_cl[0])); 
  q_cl[4] = ser_3x_p1_surfx2_eval_quad_node_3_l(Pxx_c) + q_cl[1]*q_cl[1]/q_cl[0]; 
  q_cl[5] = ser_3x_p1_surfx2_eval_quad_node_3_l(Pxy_c) + q_cl[1]*q_cl[2]/q_cl[0]; 
  q_cl[6] = ser_3x_p1_surfx2_eval_quad_node_3_l(Pxz_c) + q_cl[1]*q_cl[3]/q_cl[0]; 
  q_cl[7] = ser_3x_p1_surfx2_eval_quad_node_3_l(Pyy_c) + q_cl[2]*q_cl[2]/q_cl[0]; 
  q_cl[8] = ser_3x_p1_surfx2_eval_quad_node_3_l(Pyz_c) + q_cl[2]*q_cl[3]/q_cl[0]; 
  q_cl[9] = ser_3x_p1_surfx2_eval_quad_node_3_l(Pzz_c) + q_cl[3]*q_cl[3]/q_cl[0]; 
  q_cr[0] = ser_3x_p1_surfx2_eval_quad_node_3_r(rho_c); 
  q_cr[1] = q_cr[0]*(0.5*(ux_surf_cr[3]+ux_surf_cr[2]+ux_surf_cr[1]+ux_surf_cr[0])); 
  q_cr[2] = q_cr[0]*(0.5*(uy_surf_cr[3]+uy_surf_cr[2]+uy_surf_cr[1]+uy_surf_cr[0])); 
  q_cr[3] = q_cr[0]*(0.5*(uz_surf_cr[3]+uz_surf_cr[2]+uz_surf_cr[1]+uz_surf_cr[0])); 
  q_cr[4] = ser_3x_p1_surfx2_eval_quad_node_3_r(Pxx_c) + q_cr[1]*q_cr[1]/q_cr[0]; 
  q_cr[5] = ser_3x_p1_surfx2_eval_quad_node_3_r(Pxy_c) + q_cr[1]*q_cr[2]/q_cr[0]; 
  q_cr[6] = ser_3x_p1_surfx2_eval_quad_node_3_r(Pxz_c) + q_cr[1]*q_cr[3]/q_cr[0]; 
  q_cr[7] = ser_3x_p1_surfx2_eval_quad_node_3_r(Pyy_c) + q_cr[2]*q_cr[2]/q_cr[0]; 
  q_cr[8] = ser_3x_p1_surfx2_eval_quad_node_3_r(Pyz_c) + q_cr[2]*q_cr[3]/q_cr[0]; 
  q_cr[9] = ser_3x_p1_surfx2_eval_quad_node_3_r(Pzz_c) + q_cr[3]*q_cr[3]/q_cr[0]; 
  q_rl[0] = ser_3x_p1_surfx2_eval_quad_node_3_l(rho_r); 
  q_rl[1] = q_rl[0]*(0.5*(ux_surf_rl[3]+ux_surf_rl[2]+ux_surf_rl[1]+ux_surf_rl[0])); 
  q_rl[2] = q_rl[0]*(0.5*(uy_surf_rl[3]+uy_surf_rl[2]+uy_surf_rl[1]+uy_surf_rl[0])); 
  q_rl[3] = q_rl[0]*(0.5*(uz_surf_rl[3]+uz_surf_rl[2]+uz_surf_rl[1]+uz_surf_rl[0])); 
  q_rl[4] = ser_3x_p1_surfx2_eval_quad_node_3_l(Pxx_r) + q_rl[1]*q_rl[1]/q_rl[0]; 
  q_rl[5] = ser_3x_p1_surfx2_eval_quad_node_3_l(Pxy_r) + q_rl[1]*q_rl[2]/q_rl[0]; 
  q_rl[6] = ser_3x_p1_surfx2_eval_quad_node_3_l(Pxz_r) + q_rl[1]*q_rl[3]/q_rl[0]; 
  q_rl[7] = ser_3x_p1_surfx2_eval_quad_node_3_l(Pyy_r) + q_rl[2]*q_rl[2]/q_rl[0]; 
  q_rl[8] = ser_3x_p1_surfx2_eval_quad_node_3_l(Pyz_r) + q_rl[2]*q_rl[3]/q_rl[0]; 
  q_rl[9] = ser_3x_p1_surfx2_eval_quad_node_3_l(Pzz_r) + q_rl[3]*q_rl[3]/q_rl[0]; 

  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_l->tau1[1], geom_l->tau2[1], geom_l->norm[1], q_lr, q_lr_local); 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_l->tau1[1], geom_l->tau2[1], geom_l->norm[1], q_cl, q_cl_local); 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_r->tau1[1], geom_r->tau2[1], geom_r->norm[1], q_cr, q_cr_local); 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_r->tau1[1], geom_r->tau2[1], geom_r->norm[1], q_rl, q_rl_local); 

  delta_l[0] = q_cl_local[0] - q_lr_local[0]; 
  delta_l[1] = q_cl_local[1] - q_lr_local[1]; 
  delta_l[2] = q_cl_local[2] - q_lr_local[2]; 
  delta_l[3] = q_cl_local[3] - q_lr_local[3]; 
  delta_l[4] = q_cl_local[4] - q_lr_local[4]; 
  delta_l[5] = q_cl_local[5] - q_lr_local[5]; 
  delta_l[6] = q_cl_local[6] - q_lr_local[6]; 
  delta_l[7] = q_cl_local[7] - q_lr_local[7]; 
  delta_l[8] = q_cl_local[8] - q_lr_local[8]; 
  delta_l[9] = q_cl_local[9] - q_lr_local[9]; 
  delta_r[0] = q_rl_local[0] - q_cr_local[0]; 
  delta_r[1] = q_rl_local[1] - q_cr_local[1]; 
  delta_r[2] = q_rl_local[2] - q_cr_local[2]; 
  delta_r[3] = q_rl_local[3] - q_cr_local[3]; 
  delta_r[4] = q_rl_local[4] - q_cr_local[4]; 
  delta_r[5] = q_rl_local[5] - q_cr_local[5]; 
  delta_r[6] = q_rl_local[6] - q_cr_local[6]; 
  delta_r[7] = q_rl_local[7] - q_cr_local[7]; 
  delta_r[8] = q_rl_local[8] - q_cr_local[8]; 
  delta_r[9] = q_rl_local[9] - q_cr_local[9]; 

  my_max_speed_l = gkyl_wv_eqn_waves(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta_l, q_lr_local, q_cl_local, waves_l, speeds_l); 
  my_max_speed_r = gkyl_wv_eqn_waves(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta_r, q_cr_local, q_rl_local, waves_r, speeds_r); 
  lenr_l = geom_l->lenr[1]; 
  speeds_l[0] *= lenr_l; 
  speeds_l[1] *= lenr_l; 
  speeds_l[2] *= lenr_l; 
  speeds_l[3] *= lenr_l; 
  speeds_l[4] *= lenr_l; 
  lenr_r = geom_r->lenr[1]; 
  speeds_r[0] *= lenr_r; 
  speeds_r[1] *= lenr_r; 
  speeds_r[2] *= lenr_r; 
  speeds_r[3] *= lenr_r; 
  speeds_r[4] *= lenr_r; 

  gkyl_wv_eqn_qfluct(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, q_lr_local, q_cl_local, waves_l, speeds_l, amdq_l_local, apdq_l_local); 
  gkyl_wv_eqn_qfluct(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, q_cr_local, q_rl_local, waves_r, speeds_r, amdq_r_local, apdq_r_local); 

  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_l->tau1[1], geom_l->tau2[1], geom_l->norm[1], amdq_l_local, amdq_l); 
  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_l->tau1[1], geom_l->tau2[1], geom_l->norm[1], apdq_l_local, apdq_l); 
  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_r->tau1[1], geom_r->tau2[1], geom_r->norm[1], amdq_r_local, amdq_r); 
  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_r->tau1[1], geom_r->tau2[1], geom_r->norm[1], apdq_r_local, apdq_r); 

  amdq_rhoux_quad_l[3] = amdq_l[1]; 
  apdq_rhoux_quad_l[3] = apdq_l[1]; 
  amdq_rhouy_quad_l[3] = amdq_l[2]; 
  apdq_rhouy_quad_l[3] = apdq_l[2]; 
  amdq_rhouz_quad_l[3] = amdq_l[3]; 
  apdq_rhouz_quad_l[3] = apdq_l[3]; 

  amdq_rhoux_quad_r[3] = amdq_r[1]; 
  apdq_rhoux_quad_r[3] = apdq_r[1]; 
  amdq_rhouy_quad_r[3] = amdq_r[2]; 
  apdq_rhouy_quad_r[3] = apdq_r[2]; 
  amdq_rhouz_quad_r[3] = amdq_r[3]; 
  apdq_rhouz_quad_r[3] = apdq_r[3]; 

  ser_3x_p1_upwind_quad_to_modal(amdq_rhoux_quad_l, amdq_rhoux_l); 
  ser_3x_p1_upwind_quad_to_modal(amdq_rhouy_quad_l, amdq_rhouy_l); 
  ser_3x_p1_upwind_quad_to_modal(amdq_rhouz_quad_l, amdq_rhouz_l); 

  ser_3x_p1_upwind_quad_to_modal(apdq_rhoux_quad_l, apdq_rhoux_l); 
  ser_3x_p1_upwind_quad_to_modal(apdq_rhouy_quad_l, apdq_rhouy_l); 
  ser_3x_p1_upwind_quad_to_modal(apdq_rhouz_quad_l, apdq_rhouz_l); 

  ser_3x_p1_upwind_quad_to_modal(amdq_rhoux_quad_r, amdq_rhoux_r); 
  ser_3x_p1_upwind_quad_to_modal(amdq_rhouy_quad_r, amdq_rhouy_r); 
  ser_3x_p1_upwind_quad_to_modal(amdq_rhouz_quad_r, amdq_rhouz_r); 

  ser_3x_p1_upwind_quad_to_modal(apdq_rhoux_quad_r, apdq_rhoux_r); 
  ser_3x_p1_upwind_quad_to_modal(apdq_rhouy_quad_r, apdq_rhouy_r); 
  ser_3x_p1_upwind_quad_to_modal(apdq_rhouz_quad_r, apdq_rhouz_r); 

  outrhou0[0] += ((-0.1767766952966368*flux_rho_r[3]*ux_surf_rl[3])+0.1767766952966368*flux_rho_l[3]*ux_surf_lr[3]-0.1767766952966368*flux_rho_r[3]*ux_surf_cr[3]+0.1767766952966368*flux_rho_l[3]*ux_surf_cl[3]-0.1767766952966368*flux_rho_r[2]*ux_surf_rl[2]+0.1767766952966368*flux_rho_l[2]*ux_surf_lr[2]-0.1767766952966368*flux_rho_r[2]*ux_surf_cr[2]+0.1767766952966368*flux_rho_l[2]*ux_surf_cl[2]-0.1767766952966368*flux_rho_r[1]*ux_surf_rl[1]+0.1767766952966368*flux_rho_l[1]*ux_surf_lr[1]-0.1767766952966368*flux_rho_r[1]*ux_surf_cr[1]+0.1767766952966368*flux_rho_l[1]*ux_surf_cl[1]-0.1767766952966368*flux_rho_r[0]*ux_surf_rl[0]+0.1767766952966368*flux_rho_l[0]*ux_surf_lr[0]-0.1767766952966368*flux_rho_r[0]*ux_surf_cr[0]+0.1767766952966368*flux_rho_l[0]*ux_surf_cl[0]-0.7071067811865475*avg_p_ij_x_r[0]+0.7071067811865475*avg_p_ij_x_l[0]+0.3535533905932737*apdq_rhoux_r[0]-0.3535533905932737*(apdq_rhoux_l[0]+amdq_rhoux_r[0])+0.3535533905932737*amdq_rhoux_l[0])*dx1; 
  outrhou0[1] += ((-0.1767766952966368*flux_rho_r[2]*ux_surf_rl[3])+0.1767766952966368*flux_rho_l[2]*ux_surf_lr[3]-0.1767766952966368*flux_rho_r[2]*ux_surf_cr[3]+0.1767766952966368*flux_rho_l[2]*ux_surf_cl[3]-0.1767766952966368*(ux_surf_rl[2]+ux_surf_cr[2])*flux_rho_r[3]+0.1767766952966368*(ux_surf_lr[2]+ux_surf_cl[2])*flux_rho_l[3]-0.1767766952966368*flux_rho_r[0]*ux_surf_rl[1]+0.1767766952966368*flux_rho_l[0]*ux_surf_lr[1]-0.1767766952966368*flux_rho_r[0]*ux_surf_cr[1]+0.1767766952966368*flux_rho_l[0]*ux_surf_cl[1]-0.1767766952966368*(ux_surf_rl[0]+ux_surf_cr[0])*flux_rho_r[1]+0.1767766952966368*(ux_surf_lr[0]+ux_surf_cl[0])*flux_rho_l[1]-0.7071067811865475*avg_p_ij_x_r[1]+0.7071067811865475*avg_p_ij_x_l[1]+0.3535533905932737*apdq_rhoux_r[1]-0.3535533905932737*(apdq_rhoux_l[1]+amdq_rhoux_r[1])+0.3535533905932737*amdq_rhoux_l[1])*dx1; 
  outrhou0[2] += ((-0.3061862178478971*(flux_rho_r[3]*ux_surf_rl[3]+flux_rho_l[3]*ux_surf_lr[3]+flux_rho_r[3]*ux_surf_cr[3]+flux_rho_l[3]*ux_surf_cl[3]+flux_rho_r[2]*ux_surf_rl[2]+flux_rho_l[2]*ux_surf_lr[2]+flux_rho_r[2]*ux_surf_cr[2]+flux_rho_l[2]*ux_surf_cl[2]+flux_rho_r[1]*ux_surf_rl[1]+flux_rho_l[1]*ux_surf_lr[1]+flux_rho_r[1]*ux_surf_cr[1]+flux_rho_l[1]*ux_surf_cl[1]+flux_rho_r[0]*ux_surf_rl[0]+flux_rho_l[0]*ux_surf_lr[0]+flux_rho_r[0]*ux_surf_cr[0]+flux_rho_l[0]*ux_surf_cl[0]))-1.224744871391589*(avg_p_ij_x_r[0]+avg_p_ij_x_l[0])+0.6123724356957944*(apdq_rhoux_r[0]+apdq_rhoux_l[0])-0.6123724356957944*(amdq_rhoux_r[0]+amdq_rhoux_l[0]))*dx1; 
  outrhou0[3] += ((-0.1767766952966368*flux_rho_r[1]*ux_surf_rl[3])+0.1767766952966368*flux_rho_l[1]*ux_surf_lr[3]-0.1767766952966368*flux_rho_r[1]*ux_surf_cr[3]+0.1767766952966368*flux_rho_l[1]*ux_surf_cl[3]-0.1767766952966368*(ux_surf_rl[1]+ux_surf_cr[1])*flux_rho_r[3]+0.1767766952966368*(ux_surf_lr[1]+ux_surf_cl[1])*flux_rho_l[3]-0.1767766952966368*flux_rho_r[0]*ux_surf_rl[2]+0.1767766952966368*flux_rho_l[0]*ux_surf_lr[2]-0.1767766952966368*flux_rho_r[0]*ux_surf_cr[2]+0.1767766952966368*flux_rho_l[0]*ux_surf_cl[2]-0.1767766952966368*(ux_surf_rl[0]+ux_surf_cr[0])*flux_rho_r[2]+0.1767766952966368*(ux_surf_lr[0]+ux_surf_cl[0])*flux_rho_l[2]-0.7071067811865475*avg_p_ij_x_r[2]+0.7071067811865475*avg_p_ij_x_l[2]+0.3535533905932737*apdq_rhoux_r[2]-0.3535533905932737*(apdq_rhoux_l[2]+amdq_rhoux_r[2])+0.3535533905932737*amdq_rhoux_l[2])*dx1; 
  outrhou0[4] += ((-0.3061862178478971*(flux_rho_r[2]*ux_surf_rl[3]+flux_rho_l[2]*ux_surf_lr[3]+flux_rho_r[2]*ux_surf_cr[3]+flux_rho_l[2]*ux_surf_cl[3]+(ux_surf_rl[2]+ux_surf_cr[2])*flux_rho_r[3]+(ux_surf_lr[2]+ux_surf_cl[2])*flux_rho_l[3]+flux_rho_r[0]*ux_surf_rl[1]+flux_rho_l[0]*ux_surf_lr[1]+flux_rho_r[0]*ux_surf_cr[1]+flux_rho_l[0]*ux_surf_cl[1]+(ux_surf_rl[0]+ux_surf_cr[0])*flux_rho_r[1]+(ux_surf_lr[0]+ux_surf_cl[0])*flux_rho_l[1]))-1.224744871391589*(avg_p_ij_x_r[1]+avg_p_ij_x_l[1])+0.6123724356957944*(apdq_rhoux_r[1]+apdq_rhoux_l[1])-0.6123724356957944*(amdq_rhoux_r[1]+amdq_rhoux_l[1]))*dx1; 
  outrhou0[5] += ((-0.1767766952966368*flux_rho_r[0]*ux_surf_rl[3])+0.1767766952966368*flux_rho_l[0]*ux_surf_lr[3]-0.1767766952966368*flux_rho_r[0]*ux_surf_cr[3]+0.1767766952966368*flux_rho_l[0]*ux_surf_cl[3]-0.1767766952966368*(ux_surf_rl[0]+ux_surf_cr[0])*flux_rho_r[3]+0.1767766952966368*(ux_surf_lr[0]+ux_surf_cl[0])*flux_rho_l[3]-0.7071067811865475*avg_p_ij_x_r[3]+0.7071067811865475*avg_p_ij_x_l[3]+0.3535533905932737*apdq_rhoux_r[3]-0.3535533905932737*(apdq_rhoux_l[3]+amdq_rhoux_r[3])+0.3535533905932737*amdq_rhoux_l[3]-0.1767766952966368*flux_rho_r[1]*ux_surf_rl[2]+0.1767766952966368*flux_rho_l[1]*ux_surf_lr[2]-0.1767766952966368*flux_rho_r[1]*ux_surf_cr[2]+0.1767766952966368*flux_rho_l[1]*ux_surf_cl[2]-0.1767766952966368*(ux_surf_rl[1]+ux_surf_cr[1])*flux_rho_r[2]+0.1767766952966368*(ux_surf_lr[1]+ux_surf_cl[1])*flux_rho_l[2])*dx1; 
  outrhou0[6] += ((-0.3061862178478971*(flux_rho_r[1]*ux_surf_rl[3]+flux_rho_l[1]*ux_surf_lr[3]+flux_rho_r[1]*ux_surf_cr[3]+flux_rho_l[1]*ux_surf_cl[3]+(ux_surf_rl[1]+ux_surf_cr[1])*flux_rho_r[3]+(ux_surf_lr[1]+ux_surf_cl[1])*flux_rho_l[3]+flux_rho_r[0]*ux_surf_rl[2]+flux_rho_l[0]*ux_surf_lr[2]+flux_rho_r[0]*ux_surf_cr[2]+flux_rho_l[0]*ux_surf_cl[2]+(ux_surf_rl[0]+ux_surf_cr[0])*flux_rho_r[2]+(ux_surf_lr[0]+ux_surf_cl[0])*flux_rho_l[2]))-1.224744871391589*(avg_p_ij_x_r[2]+avg_p_ij_x_l[2])+0.6123724356957944*(apdq_rhoux_r[2]+apdq_rhoux_l[2])-0.6123724356957944*(amdq_rhoux_r[2]+amdq_rhoux_l[2]))*dx1; 
  outrhou0[7] += ((-0.3061862178478971*(flux_rho_r[0]*ux_surf_rl[3]+flux_rho_l[0]*ux_surf_lr[3]+flux_rho_r[0]*ux_surf_cr[3]+flux_rho_l[0]*ux_surf_cl[3]+(ux_surf_rl[0]+ux_surf_cr[0])*flux_rho_r[3]+(ux_surf_lr[0]+ux_surf_cl[0])*flux_rho_l[3]))-1.224744871391589*(avg_p_ij_x_r[3]+avg_p_ij_x_l[3])+0.6123724356957944*(apdq_rhoux_r[3]+apdq_rhoux_l[3])-0.6123724356957944*(amdq_rhoux_r[3]+amdq_rhoux_l[3])-0.3061862178478971*(flux_rho_r[1]*ux_surf_rl[2]+flux_rho_l[1]*ux_surf_lr[2]+flux_rho_r[1]*ux_surf_cr[2]+flux_rho_l[1]*ux_surf_cl[2]+(ux_surf_rl[1]+ux_surf_cr[1])*flux_rho_r[2]+(ux_surf_lr[1]+ux_surf_cl[1])*flux_rho_l[2]))*dx1; 

  outrhou1[0] += ((-0.1767766952966368*flux_rho_r[3]*uy_surf_rl[3])+0.1767766952966368*flux_rho_l[3]*uy_surf_lr[3]-0.1767766952966368*flux_rho_r[3]*uy_surf_cr[3]+0.1767766952966368*flux_rho_l[3]*uy_surf_cl[3]-0.1767766952966368*flux_rho_r[2]*uy_surf_rl[2]+0.1767766952966368*flux_rho_l[2]*uy_surf_lr[2]-0.1767766952966368*flux_rho_r[2]*uy_surf_cr[2]+0.1767766952966368*flux_rho_l[2]*uy_surf_cl[2]-0.1767766952966368*flux_rho_r[1]*uy_surf_rl[1]+0.1767766952966368*flux_rho_l[1]*uy_surf_lr[1]-0.1767766952966368*flux_rho_r[1]*uy_surf_cr[1]+0.1767766952966368*flux_rho_l[1]*uy_surf_cl[1]-0.1767766952966368*flux_rho_r[0]*uy_surf_rl[0]+0.1767766952966368*flux_rho_l[0]*uy_surf_lr[0]-0.1767766952966368*flux_rho_r[0]*uy_surf_cr[0]+0.1767766952966368*flux_rho_l[0]*uy_surf_cl[0]-0.7071067811865475*avg_p_ij_y_r[0]+0.7071067811865475*avg_p_ij_y_l[0]+0.3535533905932737*apdq_rhouy_r[0]-0.3535533905932737*(apdq_rhouy_l[0]+amdq_rhouy_r[0])+0.3535533905932737*amdq_rhouy_l[0])*dx1; 
  outrhou1[1] += ((-0.1767766952966368*flux_rho_r[2]*uy_surf_rl[3])+0.1767766952966368*flux_rho_l[2]*uy_surf_lr[3]-0.1767766952966368*flux_rho_r[2]*uy_surf_cr[3]+0.1767766952966368*flux_rho_l[2]*uy_surf_cl[3]-0.1767766952966368*(uy_surf_rl[2]+uy_surf_cr[2])*flux_rho_r[3]+0.1767766952966368*(uy_surf_lr[2]+uy_surf_cl[2])*flux_rho_l[3]-0.1767766952966368*flux_rho_r[0]*uy_surf_rl[1]+0.1767766952966368*flux_rho_l[0]*uy_surf_lr[1]-0.1767766952966368*flux_rho_r[0]*uy_surf_cr[1]+0.1767766952966368*flux_rho_l[0]*uy_surf_cl[1]-0.1767766952966368*(uy_surf_rl[0]+uy_surf_cr[0])*flux_rho_r[1]+0.1767766952966368*(uy_surf_lr[0]+uy_surf_cl[0])*flux_rho_l[1]-0.7071067811865475*avg_p_ij_y_r[1]+0.7071067811865475*avg_p_ij_y_l[1]+0.3535533905932737*apdq_rhouy_r[1]-0.3535533905932737*(apdq_rhouy_l[1]+amdq_rhouy_r[1])+0.3535533905932737*amdq_rhouy_l[1])*dx1; 
  outrhou1[2] += ((-0.3061862178478971*(flux_rho_r[3]*uy_surf_rl[3]+flux_rho_l[3]*uy_surf_lr[3]+flux_rho_r[3]*uy_surf_cr[3]+flux_rho_l[3]*uy_surf_cl[3]+flux_rho_r[2]*uy_surf_rl[2]+flux_rho_l[2]*uy_surf_lr[2]+flux_rho_r[2]*uy_surf_cr[2]+flux_rho_l[2]*uy_surf_cl[2]+flux_rho_r[1]*uy_surf_rl[1]+flux_rho_l[1]*uy_surf_lr[1]+flux_rho_r[1]*uy_surf_cr[1]+flux_rho_l[1]*uy_surf_cl[1]+flux_rho_r[0]*uy_surf_rl[0]+flux_rho_l[0]*uy_surf_lr[0]+flux_rho_r[0]*uy_surf_cr[0]+flux_rho_l[0]*uy_surf_cl[0]))-1.224744871391589*(avg_p_ij_y_r[0]+avg_p_ij_y_l[0])+0.6123724356957944*(apdq_rhouy_r[0]+apdq_rhouy_l[0])-0.6123724356957944*(amdq_rhouy_r[0]+amdq_rhouy_l[0]))*dx1; 
  outrhou1[3] += ((-0.1767766952966368*flux_rho_r[1]*uy_surf_rl[3])+0.1767766952966368*flux_rho_l[1]*uy_surf_lr[3]-0.1767766952966368*flux_rho_r[1]*uy_surf_cr[3]+0.1767766952966368*flux_rho_l[1]*uy_surf_cl[3]-0.1767766952966368*(uy_surf_rl[1]+uy_surf_cr[1])*flux_rho_r[3]+0.1767766952966368*(uy_surf_lr[1]+uy_surf_cl[1])*flux_rho_l[3]-0.1767766952966368*flux_rho_r[0]*uy_surf_rl[2]+0.1767766952966368*flux_rho_l[0]*uy_surf_lr[2]-0.1767766952966368*flux_rho_r[0]*uy_surf_cr[2]+0.1767766952966368*flux_rho_l[0]*uy_surf_cl[2]-0.1767766952966368*(uy_surf_rl[0]+uy_surf_cr[0])*flux_rho_r[2]+0.1767766952966368*(uy_surf_lr[0]+uy_surf_cl[0])*flux_rho_l[2]-0.7071067811865475*avg_p_ij_y_r[2]+0.7071067811865475*avg_p_ij_y_l[2]+0.3535533905932737*apdq_rhouy_r[2]-0.3535533905932737*(apdq_rhouy_l[2]+amdq_rhouy_r[2])+0.3535533905932737*amdq_rhouy_l[2])*dx1; 
  outrhou1[4] += ((-0.3061862178478971*(flux_rho_r[2]*uy_surf_rl[3]+flux_rho_l[2]*uy_surf_lr[3]+flux_rho_r[2]*uy_surf_cr[3]+flux_rho_l[2]*uy_surf_cl[3]+(uy_surf_rl[2]+uy_surf_cr[2])*flux_rho_r[3]+(uy_surf_lr[2]+uy_surf_cl[2])*flux_rho_l[3]+flux_rho_r[0]*uy_surf_rl[1]+flux_rho_l[0]*uy_surf_lr[1]+flux_rho_r[0]*uy_surf_cr[1]+flux_rho_l[0]*uy_surf_cl[1]+(uy_surf_rl[0]+uy_surf_cr[0])*flux_rho_r[1]+(uy_surf_lr[0]+uy_surf_cl[0])*flux_rho_l[1]))-1.224744871391589*(avg_p_ij_y_r[1]+avg_p_ij_y_l[1])+0.6123724356957944*(apdq_rhouy_r[1]+apdq_rhouy_l[1])-0.6123724356957944*(amdq_rhouy_r[1]+amdq_rhouy_l[1]))*dx1; 
  outrhou1[5] += ((-0.1767766952966368*flux_rho_r[0]*uy_surf_rl[3])+0.1767766952966368*flux_rho_l[0]*uy_surf_lr[3]-0.1767766952966368*flux_rho_r[0]*uy_surf_cr[3]+0.1767766952966368*flux_rho_l[0]*uy_surf_cl[3]-0.1767766952966368*(uy_surf_rl[0]+uy_surf_cr[0])*flux_rho_r[3]+0.1767766952966368*(uy_surf_lr[0]+uy_surf_cl[0])*flux_rho_l[3]-0.7071067811865475*avg_p_ij_y_r[3]+0.7071067811865475*avg_p_ij_y_l[3]+0.3535533905932737*apdq_rhouy_r[3]-0.3535533905932737*(apdq_rhouy_l[3]+amdq_rhouy_r[3])+0.3535533905932737*amdq_rhouy_l[3]-0.1767766952966368*flux_rho_r[1]*uy_surf_rl[2]+0.1767766952966368*flux_rho_l[1]*uy_surf_lr[2]-0.1767766952966368*flux_rho_r[1]*uy_surf_cr[2]+0.1767766952966368*flux_rho_l[1]*uy_surf_cl[2]-0.1767766952966368*(uy_surf_rl[1]+uy_surf_cr[1])*flux_rho_r[2]+0.1767766952966368*(uy_surf_lr[1]+uy_surf_cl[1])*flux_rho_l[2])*dx1; 
  outrhou1[6] += ((-0.3061862178478971*(flux_rho_r[1]*uy_surf_rl[3]+flux_rho_l[1]*uy_surf_lr[3]+flux_rho_r[1]*uy_surf_cr[3]+flux_rho_l[1]*uy_surf_cl[3]+(uy_surf_rl[1]+uy_surf_cr[1])*flux_rho_r[3]+(uy_surf_lr[1]+uy_surf_cl[1])*flux_rho_l[3]+flux_rho_r[0]*uy_surf_rl[2]+flux_rho_l[0]*uy_surf_lr[2]+flux_rho_r[0]*uy_surf_cr[2]+flux_rho_l[0]*uy_surf_cl[2]+(uy_surf_rl[0]+uy_surf_cr[0])*flux_rho_r[2]+(uy_surf_lr[0]+uy_surf_cl[0])*flux_rho_l[2]))-1.224744871391589*(avg_p_ij_y_r[2]+avg_p_ij_y_l[2])+0.6123724356957944*(apdq_rhouy_r[2]+apdq_rhouy_l[2])-0.6123724356957944*(amdq_rhouy_r[2]+amdq_rhouy_l[2]))*dx1; 
  outrhou1[7] += ((-0.3061862178478971*(flux_rho_r[0]*uy_surf_rl[3]+flux_rho_l[0]*uy_surf_lr[3]+flux_rho_r[0]*uy_surf_cr[3]+flux_rho_l[0]*uy_surf_cl[3]+(uy_surf_rl[0]+uy_surf_cr[0])*flux_rho_r[3]+(uy_surf_lr[0]+uy_surf_cl[0])*flux_rho_l[3]))-1.224744871391589*(avg_p_ij_y_r[3]+avg_p_ij_y_l[3])+0.6123724356957944*(apdq_rhouy_r[3]+apdq_rhouy_l[3])-0.6123724356957944*(amdq_rhouy_r[3]+amdq_rhouy_l[3])-0.3061862178478971*(flux_rho_r[1]*uy_surf_rl[2]+flux_rho_l[1]*uy_surf_lr[2]+flux_rho_r[1]*uy_surf_cr[2]+flux_rho_l[1]*uy_surf_cl[2]+(uy_surf_rl[1]+uy_surf_cr[1])*flux_rho_r[2]+(uy_surf_lr[1]+uy_surf_cl[1])*flux_rho_l[2]))*dx1; 

  outrhou2[0] += ((-0.1767766952966368*flux_rho_r[3]*uz_surf_rl[3])+0.1767766952966368*flux_rho_l[3]*uz_surf_lr[3]-0.1767766952966368*flux_rho_r[3]*uz_surf_cr[3]+0.1767766952966368*flux_rho_l[3]*uz_surf_cl[3]-0.1767766952966368*flux_rho_r[2]*uz_surf_rl[2]+0.1767766952966368*flux_rho_l[2]*uz_surf_lr[2]-0.1767766952966368*flux_rho_r[2]*uz_surf_cr[2]+0.1767766952966368*flux_rho_l[2]*uz_surf_cl[2]-0.1767766952966368*flux_rho_r[1]*uz_surf_rl[1]+0.1767766952966368*flux_rho_l[1]*uz_surf_lr[1]-0.1767766952966368*flux_rho_r[1]*uz_surf_cr[1]+0.1767766952966368*flux_rho_l[1]*uz_surf_cl[1]-0.1767766952966368*flux_rho_r[0]*uz_surf_rl[0]+0.1767766952966368*flux_rho_l[0]*uz_surf_lr[0]-0.1767766952966368*flux_rho_r[0]*uz_surf_cr[0]+0.1767766952966368*flux_rho_l[0]*uz_surf_cl[0]-0.7071067811865475*avg_p_ij_z_r[0]+0.7071067811865475*avg_p_ij_z_l[0]+0.3535533905932737*apdq_rhouz_r[0]-0.3535533905932737*(apdq_rhouz_l[0]+amdq_rhouz_r[0])+0.3535533905932737*amdq_rhouz_l[0])*dx1; 
  outrhou2[1] += ((-0.1767766952966368*flux_rho_r[2]*uz_surf_rl[3])+0.1767766952966368*flux_rho_l[2]*uz_surf_lr[3]-0.1767766952966368*flux_rho_r[2]*uz_surf_cr[3]+0.1767766952966368*flux_rho_l[2]*uz_surf_cl[3]-0.1767766952966368*(uz_surf_rl[2]+uz_surf_cr[2])*flux_rho_r[3]+0.1767766952966368*(uz_surf_lr[2]+uz_surf_cl[2])*flux_rho_l[3]-0.1767766952966368*flux_rho_r[0]*uz_surf_rl[1]+0.1767766952966368*flux_rho_l[0]*uz_surf_lr[1]-0.1767766952966368*flux_rho_r[0]*uz_surf_cr[1]+0.1767766952966368*flux_rho_l[0]*uz_surf_cl[1]-0.1767766952966368*(uz_surf_rl[0]+uz_surf_cr[0])*flux_rho_r[1]+0.1767766952966368*(uz_surf_lr[0]+uz_surf_cl[0])*flux_rho_l[1]-0.7071067811865475*avg_p_ij_z_r[1]+0.7071067811865475*avg_p_ij_z_l[1]+0.3535533905932737*apdq_rhouz_r[1]-0.3535533905932737*(apdq_rhouz_l[1]+amdq_rhouz_r[1])+0.3535533905932737*amdq_rhouz_l[1])*dx1; 
  outrhou2[2] += ((-0.3061862178478971*(flux_rho_r[3]*uz_surf_rl[3]+flux_rho_l[3]*uz_surf_lr[3]+flux_rho_r[3]*uz_surf_cr[3]+flux_rho_l[3]*uz_surf_cl[3]+flux_rho_r[2]*uz_surf_rl[2]+flux_rho_l[2]*uz_surf_lr[2]+flux_rho_r[2]*uz_surf_cr[2]+flux_rho_l[2]*uz_surf_cl[2]+flux_rho_r[1]*uz_surf_rl[1]+flux_rho_l[1]*uz_surf_lr[1]+flux_rho_r[1]*uz_surf_cr[1]+flux_rho_l[1]*uz_surf_cl[1]+flux_rho_r[0]*uz_surf_rl[0]+flux_rho_l[0]*uz_surf_lr[0]+flux_rho_r[0]*uz_surf_cr[0]+flux_rho_l[0]*uz_surf_cl[0]))-1.224744871391589*(avg_p_ij_z_r[0]+avg_p_ij_z_l[0])+0.6123724356957944*(apdq_rhouz_r[0]+apdq_rhouz_l[0])-0.6123724356957944*(amdq_rhouz_r[0]+amdq_rhouz_l[0]))*dx1; 
  outrhou2[3] += ((-0.1767766952966368*flux_rho_r[1]*uz_surf_rl[3])+0.1767766952966368*flux_rho_l[1]*uz_surf_lr[3]-0.1767766952966368*flux_rho_r[1]*uz_surf_cr[3]+0.1767766952966368*flux_rho_l[1]*uz_surf_cl[3]-0.1767766952966368*(uz_surf_rl[1]+uz_surf_cr[1])*flux_rho_r[3]+0.1767766952966368*(uz_surf_lr[1]+uz_surf_cl[1])*flux_rho_l[3]-0.1767766952966368*flux_rho_r[0]*uz_surf_rl[2]+0.1767766952966368*flux_rho_l[0]*uz_surf_lr[2]-0.1767766952966368*flux_rho_r[0]*uz_surf_cr[2]+0.1767766952966368*flux_rho_l[0]*uz_surf_cl[2]-0.1767766952966368*(uz_surf_rl[0]+uz_surf_cr[0])*flux_rho_r[2]+0.1767766952966368*(uz_surf_lr[0]+uz_surf_cl[0])*flux_rho_l[2]-0.7071067811865475*avg_p_ij_z_r[2]+0.7071067811865475*avg_p_ij_z_l[2]+0.3535533905932737*apdq_rhouz_r[2]-0.3535533905932737*(apdq_rhouz_l[2]+amdq_rhouz_r[2])+0.3535533905932737*amdq_rhouz_l[2])*dx1; 
  outrhou2[4] += ((-0.3061862178478971*(flux_rho_r[2]*uz_surf_rl[3]+flux_rho_l[2]*uz_surf_lr[3]+flux_rho_r[2]*uz_surf_cr[3]+flux_rho_l[2]*uz_surf_cl[3]+(uz_surf_rl[2]+uz_surf_cr[2])*flux_rho_r[3]+(uz_surf_lr[2]+uz_surf_cl[2])*flux_rho_l[3]+flux_rho_r[0]*uz_surf_rl[1]+flux_rho_l[0]*uz_surf_lr[1]+flux_rho_r[0]*uz_surf_cr[1]+flux_rho_l[0]*uz_surf_cl[1]+(uz_surf_rl[0]+uz_surf_cr[0])*flux_rho_r[1]+(uz_surf_lr[0]+uz_surf_cl[0])*flux_rho_l[1]))-1.224744871391589*(avg_p_ij_z_r[1]+avg_p_ij_z_l[1])+0.6123724356957944*(apdq_rhouz_r[1]+apdq_rhouz_l[1])-0.6123724356957944*(amdq_rhouz_r[1]+amdq_rhouz_l[1]))*dx1; 
  outrhou2[5] += ((-0.1767766952966368*flux_rho_r[0]*uz_surf_rl[3])+0.1767766952966368*flux_rho_l[0]*uz_surf_lr[3]-0.1767766952966368*flux_rho_r[0]*uz_surf_cr[3]+0.1767766952966368*flux_rho_l[0]*uz_surf_cl[3]-0.1767766952966368*(uz_surf_rl[0]+uz_surf_cr[0])*flux_rho_r[3]+0.1767766952966368*(uz_surf_lr[0]+uz_surf_cl[0])*flux_rho_l[3]-0.7071067811865475*avg_p_ij_z_r[3]+0.7071067811865475*avg_p_ij_z_l[3]+0.3535533905932737*apdq_rhouz_r[3]-0.3535533905932737*(apdq_rhouz_l[3]+amdq_rhouz_r[3])+0.3535533905932737*amdq_rhouz_l[3]-0.1767766952966368*flux_rho_r[1]*uz_surf_rl[2]+0.1767766952966368*flux_rho_l[1]*uz_surf_lr[2]-0.1767766952966368*flux_rho_r[1]*uz_surf_cr[2]+0.1767766952966368*flux_rho_l[1]*uz_surf_cl[2]-0.1767766952966368*(uz_surf_rl[1]+uz_surf_cr[1])*flux_rho_r[2]+0.1767766952966368*(uz_surf_lr[1]+uz_surf_cl[1])*flux_rho_l[2])*dx1; 
  outrhou2[6] += ((-0.3061862178478971*(flux_rho_r[1]*uz_surf_rl[3]+flux_rho_l[1]*uz_surf_lr[3]+flux_rho_r[1]*uz_surf_cr[3]+flux_rho_l[1]*uz_surf_cl[3]+(uz_surf_rl[1]+uz_surf_cr[1])*flux_rho_r[3]+(uz_surf_lr[1]+uz_surf_cl[1])*flux_rho_l[3]+flux_rho_r[0]*uz_surf_rl[2]+flux_rho_l[0]*uz_surf_lr[2]+flux_rho_r[0]*uz_surf_cr[2]+flux_rho_l[0]*uz_surf_cl[2]+(uz_surf_rl[0]+uz_surf_cr[0])*flux_rho_r[2]+(uz_surf_lr[0]+uz_surf_cl[0])*flux_rho_l[2]))-1.224744871391589*(avg_p_ij_z_r[2]+avg_p_ij_z_l[2])+0.6123724356957944*(apdq_rhouz_r[2]+apdq_rhouz_l[2])-0.6123724356957944*(amdq_rhouz_r[2]+amdq_rhouz_l[2]))*dx1; 
  outrhou2[7] += ((-0.3061862178478971*(flux_rho_r[0]*uz_surf_rl[3]+flux_rho_l[0]*uz_surf_lr[3]+flux_rho_r[0]*uz_surf_cr[3]+flux_rho_l[0]*uz_surf_cl[3]+(uz_surf_rl[0]+uz_surf_cr[0])*flux_rho_r[3]+(uz_surf_lr[0]+uz_surf_cl[0])*flux_rho_l[3]))-1.224744871391589*(avg_p_ij_z_r[3]+avg_p_ij_z_l[3])+0.6123724356957944*(apdq_rhouz_r[3]+apdq_rhouz_l[3])-0.6123724356957944*(amdq_rhouz_r[3]+amdq_rhouz_l[3])-0.3061862178478971*(flux_rho_r[1]*uz_surf_rl[2]+flux_rho_l[1]*uz_surf_lr[2]+flux_rho_r[1]*uz_surf_cr[2]+flux_rho_l[1]*uz_surf_cl[2]+(uz_surf_rl[1]+uz_surf_cr[1])*flux_rho_r[2]+(uz_surf_lr[1]+uz_surf_cl[1])*flux_rho_l[2]))*dx1; 

  return 0.;

} 
