#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>

/**
 * Compute u from state vector which contains *both* density and momentum.
 * Data always organized such that density is 0th component of statevec
 * & momentum is 1st, 2nd, and 3rd component of statevec.
 * Examples include: 
 * Isothermal Euler state vector [rho, rhoux, rhouy, rhouz],
 * Euler state vector [rho, rhoux, rhouy, rhouz, E],
 * Nonrelativistic Vlasov Five Moments array [M0, M1i, M2] (u_i = M1i/M0),
 * Special Relativistic Vlasov N_i array [Gamma*n, Gamma*n*u] where u is the bulk flow
 *
 * @param mem Pre-allocated space for use in the division 
 * @param basis Basis functions used in expansions
 * @param range Range to apply division operator
 * @param statevec Input state vector which contains *both* density and momentum
 * @param u_i Output array of bulk flow velocity
 */
void gkyl_calc_prim_vars_u_from_statevec(gkyl_dg_bin_op_mem *mem, struct gkyl_basis basis, const struct gkyl_range *range,
  const struct gkyl_array* statevec, struct gkyl_array* u_i);

/**
 * Compute u from input density and momentum vectors.
 * Examples include: 
 * Nonrelativistic Vlasov M0 & M1i (u_i = M1i/M0),
 * Special Relativistic M0 (Gamma*n) and M1i (Gamma*n*u) where u is the bulk flow
 * Parallel-kinetic-perpendicular-moment model with vlasov_pkpm_moms and euler_pkpm
 *
 * @param mem Pre-allocated space for use in the division 
 * @param basis Basis functions used in expansions
 * @param range Range to apply division operator
 * @param rho Input density 
 * @param rhou Input momentum
 * @param u_i Output array of bulk flow velocity
 */
void gkyl_calc_prim_vars_u_from_rhou(gkyl_dg_bin_op_mem *mem, struct gkyl_basis basis, const struct gkyl_range *range,
  const struct gkyl_array* rho, const struct gkyl_array* rhou, struct gkyl_array* u_i);

/**
 * Compute pressure from state vector which contains density, momentum, and energy.
 * Data always organized such that density is 0th component of statevec, momentum is 
 * 1st, 2nd, and 3rd component of statevec, and energy is 4th component of statevec.
 * Note: Pressure computation requires computation of u, bulk flow velocity.
 * Examples include: 
 * Euler state vector [rho, rhoux, rhouy, rhouz, E],
 * Nonrelativistic Vlasov Five Moments array [M0, M1i, M2]
 *
 * @param basis Basis functions used in expansions
 * @param range Range to apply division operator
 * @param p_fac Factor for obtaining pressure (gas_gamma for Euler)
 * @param u_i Input array of bulk flow velocity
 * @param statevec Input state vector which contains density, momentum, and energy
 * @param p_ij Output array of pressure (scalar pressure)
 */
void gkyl_calc_prim_vars_p_from_statevec(struct gkyl_basis basis, const struct gkyl_range *range,
  const double p_fac, const struct gkyl_array* u_i, const struct gkyl_array* statevec, 
  struct gkyl_array* p_ij);

/**
 * Compute primitive variables for parallel-kinetic-perpendicular-moment model.
 *
 * @param basis Basis functions used in expansions
 * @param range Range to apply division operator
 * @param bvar Input array of magnetic field unit vector and unit tensor
 * @param vlasov_pkpm_moms Input array parallel-kinetic-perpendicular-moment kinetic moments
 * @param euler_pkpm Input array parallel-kinetic-perpendicular-moment fluid variables
 * @param u_i Output array of flow velocity 
 * @param u_perp_i Output array of perpendicular flow velocity (u - u : bb)
 * @param rhou_perp_i Output array of perpendicular momentum density (rhou - rhou : bb)
 * @param p_perp Output array of perpendicular pressure
 * @param p_ij Output array of pressure tensor
 */
void gkyl_calc_prim_vars_pkpm(struct gkyl_basis basis, const struct gkyl_range *range,
  const struct gkyl_array* bvar, const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm, 
  struct gkyl_array* u_i, struct gkyl_array* u_perp_i, struct gkyl_array* rhou_perp_i,
  struct gkyl_array* p_perp, struct gkyl_array* p_ij);

/**
 * Compute parallel-kinetic-perpendicular-moment model source terms.
 *
 * @param basis Basis functions used in expansions
 * @param range Range to apply division operator
 * @param qmem Input array of q/m*EM fields
 * @param nu Input array of collisionality
 * @param nu_vthsq Input array of nu*vth^2, vth^2 = T/m
 * @param vlasov_pkpm_moms Input array parallel-kinetic-perpendicular-moment kinetic moments
 * @param euler_pkpm Input array parallel-kinetic-perpendicular-moment fluid variables
 * @param rhou_perp_i Input array of perpendicular momentum density (rhou - rhou : bb)
 * @param p_perp Input array of perpendicular pressure
 */
void gkyl_calc_prim_vars_pkpm_source(struct gkyl_basis basis, const struct gkyl_range *range,
  const struct gkyl_array* qmem, const struct gkyl_array* nu, const struct gkyl_array* nu_vthsq, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm,
  const struct gkyl_array* rhou_perp_i, const struct gkyl_array* p_perp, 
  struct gkyl_array* rhs);

/**
 * Compute needed gradient quantities with recovery for discretization of the 
 * parallel-kinetic-perpendicular-moment (pkpm) model. These include div(p), div(b), and bb : grad(u). 
 *
 * @param grid Grid (for getting cell spacing)
 * @param basis Basis functions used in expansions
 * @param range Range to apply division operator
 * @param bvar Input array of magnetic field unit vector and unit tensor
 * @param u_i Input array of flow velocity 
 * @param p_ij Input array of pressure tensor
 * @param div_b Output array of divergence of magnetic field unit vector
 * @param bb_grad_u Output array of bb : grad(u)
 * @param div_p Output array of divergence of pressure tensor
 */
void gkyl_calc_prim_vars_pkpm_recovery(const struct gkyl_rect_grid *grid, 
  struct gkyl_basis basis, const struct gkyl_range *range,
  const struct gkyl_array* bvar, const struct gkyl_array* u_i, const struct gkyl_array* p_ij, 
  struct gkyl_array* div_b, struct gkyl_array* bb_grad_u, struct gkyl_array* div_p);

/**
 * Compute upwinded divergence quantities for discretization of the 
 * parallel-kinetic-perpendicular-moment (pkpm) model. These include:
 * 1. div (integral v_parallel^2 b_hat F) (pressure force)
 *
 * @param phase_grid Phase space grid (for getting cell spacing and cell center)
 * @param basis Basis functions used in expansions
 * @param phase_range Phase space range 
 * @param conf_range Configuration space range
 * @param bvar Input array of magnetic field unit vector and unit tensor
 * @param f Input array of distribution function
 * @param p_force Output array of pressure force
 */
void gkyl_calc_prim_vars_pkpm_upwind_p(const struct gkyl_rect_grid *phase_grid, 
  struct gkyl_basis cbasis, const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  double mass, const struct gkyl_array* bvar, const struct gkyl_array* f, struct gkyl_array* p_force);

/**
 * Compute total pressure forces for parallel-kinetic-perpendicular-moment forces.
 * Total pressure force p_force = 1/rho (b . div(P) + p_perp div(b))
 *
 * @param basis Basis functions used in expansions
 * @param range Range to apply division operator
 * @param bvar Input array of magnetic field unit vector and unit tensor
 * @param div_p Input array of divergence of pressure tensor
 * @param vlasov_pkpm_moms Input array parallel-kinetic-perpendicular-moment kinetic moments
 * @param euler_pkpm Input array parallel-kinetic-perpendicular-moment fluid variables
 * @param div_b Input array of divergence of magnetic field unit vector
 * @param p_force Output array of total pressure forces
 */
void gkyl_calc_prim_vars_pkpm_p_force(struct gkyl_basis basis, const struct gkyl_range *range,
  const struct gkyl_array* bvar, const struct gkyl_array* div_p, const struct gkyl_array* vlasov_pkpm_moms, 
  const struct gkyl_array* euler_pkpm, const struct gkyl_array* div_b, struct gkyl_array* p_force);

/**
 * Host-side wrappers for prim vars operations on device
 */

void gkyl_calc_prim_vars_pkpm_cu(struct gkyl_basis basis, const struct gkyl_range *range,
  const struct gkyl_array* bvar, const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm, 
  struct gkyl_array* u_i, struct gkyl_array* u_perp_i, struct gkyl_array* rhou_perp_i,
  struct gkyl_array* p_perp, struct gkyl_array* p_ij);

void gkyl_calc_prim_vars_pkpm_source_cu(struct gkyl_basis basis, const struct gkyl_range *range,
  const struct gkyl_array* qmem, const struct gkyl_array* nu, const struct gkyl_array* nu_vthsq, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm,
  const struct gkyl_array* rhou_perp_i, const struct gkyl_array* p_perp, 
  struct gkyl_array* rhs);

void gkyl_calc_prim_vars_pkpm_recovery_cu(const struct gkyl_rect_grid *grid, 
  struct gkyl_basis basis, const struct gkyl_range *range,
  const struct gkyl_array* bvar, const struct gkyl_array* u_i, const struct gkyl_array* p_ij, 
  struct gkyl_array* div_b, struct gkyl_array* bb_grad_u, struct gkyl_array* div_p);

void gkyl_calc_prim_vars_pkpm_upwind_p_cu(const struct gkyl_rect_grid *phase_grid, 
  struct gkyl_basis cbasis, const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array* bvar, const struct gkyl_array* f, struct gkyl_array* p_force);

void gkyl_calc_prim_vars_pkpm_p_force_cu(struct gkyl_basis basis, const struct gkyl_range *range,
  const struct gkyl_array* bvar, const struct gkyl_array* div_p, const struct gkyl_array* vlasov_pkpm_moms, 
  const struct gkyl_array* euler_pkpm, const struct gkyl_array* div_b, struct gkyl_array* p_force);
