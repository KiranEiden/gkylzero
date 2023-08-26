#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_basis.h>

// Object type
typedef struct gkyl_dg_calc_pkpm_vars gkyl_dg_calc_pkpm_vars;

/**
 * Create new updater to compute pkpm variables needed in 
 * updates and used for diagnostics. Methods compute:
 * p_ij : (p_par - p_perp) b_i b_j + p_perp g_ij
 * prim : [ux, uy, uz, 1/rho*div(p_par b), T_perp/m, m/T_perp]
 * pkpm_accel_vars : 0 : p_perp_div_b (p_perp/rho*div(b) = T_perp/m*div(b))
                     1 : bb_grad_u (bb : grad(u))
                     2 : p_force (total pressure forces in kinetic equation 1/rho div(p_parallel b_hat) - T_perp/m*div(b)
                     3 : p_perp_source (pressure source for higher Laguerre moments -> bb : grad(u) - div(u) - 2 nu)
 * pkpm_int_vars : integrated PKPM variables (rho, rhoux, rhouy, rhouz, rhoux^2, rhouy^2, rhouz^2, p_parallel, p_perp)
 * pkpm fluid source : Explicit source terms in momentum solve q/m rho (E_i + epsilon_ijk u_j B_k)
 * 
 * Updater also stores the kernels to compute pkpm source terms and pkpm integrated moments.
 * 
 * @param conf_grid Configuration space grid (for getting cell spacing and cell center)
 * @param cbasis Configuration space basis functions
 * @param mem_range Configuration space range that sets the size of the bin_op memory
 *                  for computing primitive moments. Note range is stored so 
 *                  updater loops over consistent range for primitive moments
 * @param use_gpu bool to determine if on GPU
 * @return New updater pointer.
 */
struct gkyl_dg_calc_pkpm_vars* 
gkyl_dg_calc_pkpm_vars_new(const struct gkyl_rect_grid *conf_grid, 
  const struct gkyl_basis* cbasis, const struct gkyl_range *mem_range, 
  bool use_gpu);

/**
 * Create new updater to compute pkpm variables on
 * NV-GPU. See new() method for documentation.
 */
struct gkyl_dg_calc_pkpm_vars* 
gkyl_dg_calc_pkpm_vars_cu_dev_new(const struct gkyl_rect_grid *conf_grid, 
  const struct gkyl_basis* cbasis, const struct gkyl_range *mem_range);

/**
 * Compute pkpm primitive moments.
 *
 * @param up Updater for computing pkpm variables 
 * @param conf_range Configuration space range
 * @param vlasov_pkpm_moms Input array of pkpm kinetic moments [rho, p_parallel, p_perp]
 * @param euler_pkpm Input array of pkpm fluid variables [rho ux, rho uy, rho uz]
 * @param pkpm_div_ppar Input array of div(p_parallel b_hat) for computing pressure force
 * @param cell_avg_prim Array for storing boolean value of whether rho, p_parallel, p_perp uses *only* cell averages 
 *                      to minimize positivity violations (default: false)
 * @param prim Output array of primitive moments [ux, uy, uz, 1/rho*div(p_par b), T_perp/m, m/T_perp]
 */
void gkyl_dg_calc_pkpm_vars_advance(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm, 
  const struct gkyl_array* pkpm_div_ppar, struct gkyl_array* cell_avg_prim, 
  struct gkyl_array* prim);

/**
 * Compute pkpm surface primitive moments.
 * 2*cdim*3+2*cdim components: ux, uy, uz (3 components) at the left and right of the cell (2 components) in each dimension (cdim components)
 * Also solves for 3*Txx/m at the left and right x surfaces, 3*Tyy/m at the left and right y surfaces, and 3*Tzz/m at the left and right z surfaces
 *
 * @param up Updater for computing pkpm variables 
 * @param conf_range Configuration space range
 * @param vlasov_pkpm_moms Input array of pkpm kinetic moments [rho, p_parallel, p_perp]
 * @param euler_pkpm Input array of pkpm fluid variables [rho ux, rho uy, rho uz]
 * @param p_ij Input pressure tensor p_ij = (p_par - p_perp) b_i b_j + p_perp g_ij
 * @param cell_avg_prim Array for storing boolean value of whether rho, p_parallel, p_perp uses *only* cell averages 
 *                      to minimize positivity violations (default: false)
 * @param prim_surf Output array of surface primitive moments
 */
void gkyl_dg_calc_pkpm_vars_surf_advance(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm, 
  const struct gkyl_array* p_ij, const struct gkyl_array* cell_avg_prim, 
  struct gkyl_array* prim_surf);

/**
 * Compute pkpm pressure p_ij = (p_par - p_perp) b_i b_j + p_perp g_ij in the volume and at needed surfaces
 *
 * @param up Updater for computing pkpm variables 
 * @param conf_range Configuration space range
 * @param bvar Input array of volume expansion of magnetic field unit vector and unit tensor
 * @param bvar_surf Input array of surface expansion of magnetic field unit vector and unit tensor
 * @param vlasov_pkpm_moms Input array of pkpm kinetic moments [rho, p_parallel, p_perp]
 * @param p_ij Output array of volume expansion of pressure tensor p_ij = (p_par - p_perp) b_i b_j + p_perp g_ij
 * @param p_ij Output array of surface expansion of pressure tensor [Pxx_xl, Pxx_xr, Pxy_xl, Pxy_xr, Pxz_xl, Pxz_xr,
 *                                                                   Pyy_yl, Pyy_yr, Pxy_yl, Pxy_yr, Pyz_yl, Pyz_yr, 
 *                                                                   Pzz_zl, Pzz_zr, Pxz_zl, Pxz_zr, Pyz_zl, Pyz_zr] 
 */
void gkyl_dg_calc_pkpm_vars_pressure(struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range, 
  const struct gkyl_array* bvar, const struct gkyl_array* bvar_surf, const struct gkyl_array* vlasov_pkpm_moms, 
  struct gkyl_array* p_ij, struct gkyl_array* p_ij_surf);

/**
 * Compute pkpm acceleration variables
 *
 * @param up Updater for computing pkpm variables 
 * @param conf_range Configuration space range
 * @param prim_surf Input array of surface expansions of primitive moments [u_i, 3.0*T_ii/m]
 * @param prim Input array of volume expansion of primitive moments [ux, uy, uz, 1/rho*div(p_par b), T_perp/m, m/T_perp]
 * @param bvar Input array of magnetic field unit vector and unit tensor
 * @param div_b Input array of div(b)
 * @param nu Input array of collisionality
 * @param pkpm_lax Output array of surface expansion of Lax penalization lambda_i = |u_i| + sqrt(3.0*T_ii/m)
 * @param pkpm_accel Output arrary of pkpm acceleration variables ordered as:
 *        0: p_perp_div_b (p_perp/rho*div(b) = T_perp/m*div(b))
          1: bb_grad_u (bb : grad(u))
          2: p_force (total pressure forces in kinetic equation 1/rho div(p_parallel b_hat) - T_perp/m*div(b)
          3: p_perp_source (pressure source for higher Laguerre moments -> bb : grad(u) - div(u) - 2 nu)
 */
void gkyl_dg_calc_pkpm_vars_accel(struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range, 
  const struct gkyl_array* prim_surf, const struct gkyl_array* prim, 
  const struct gkyl_array* bvar, const struct gkyl_array* div_b, const struct gkyl_array* nu, 
  struct gkyl_array* pkpm_lax, struct gkyl_array* pkpm_accel);

/**
 * Compute integrated PKPM variables (rho, rhoux, rhouy, rhouz, rhoux^2, rhouy^2, rhouz^2, p_parallel, p_perp).
 *
 * @param up Updater for computing pkpm variables 
 * @param conf_range Configuration space range
 * @param vlasov_pkpm_moms Input array of parallel-kinetic-perpendicular-moment kinetic moments [rho, p_parallel, p_perp]
 * @param euler_pkpm Input array of parallel-kinetic-perpendicular-moment fluid variables [rho ux, rho uy, rho uz]
 * @param prim Input array of primitive moments [ux, uy, uz, 1/rho*div(p_par b), T_perp/m, m/T_perp]
 * @param int_pkpm_vars Output array of integrated variables (6 components)
 */
void gkyl_dg_calc_pkpm_integrated_vars(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_array* vlasov_pkpm_moms, 
  const struct gkyl_array* euler_pkpm, const struct gkyl_array* prim, 
  struct gkyl_array* pkpm_int_vars);

/**
 * Compute pkpm model source terms.
 *
 * @param up Updater for computing pkpm variables 
 * @param conf_range Configuration space range
 * @param qmem Input array of q/m*EM fields
 * @param vlasov_pkpm_moms Input array of parallel-kinetic-perpendicular-moment kinetic moments [rho, p_parallel, p_perp]
 * @param euler_pkpm Input array of parallel-kinetic-perpendicular-moment fluid variables [rho ux, rho uy, rho uz]
 * @param rhs Output increment to fluid variables
 */
void gkyl_dg_calc_pkpm_vars_source(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_array* qmem, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm,
  struct gkyl_array* rhs);

/**
 * Construct PKPM variables for I/O. Computes the conserved fluid variables 
 * [rho, rho ux, rho uy, rho uz, Pxx + rho ux^2, Pxy + rho ux uy, Pxz + rho ux uz, Pyy + rho uy^2, Pyz + rho uy uz, Pzz + rho uz^2]
 * And copies the pkpm primitive and acceleration variables into an array for output
 * [ux, uy, uz, T_perp/m, m/T_perp, 1/rho div(p_par b), T_perp/m div(b), bb : grad(u)]
 *
 * @param up Updater for computing pkpm variables 
 * @param conf_range Configuration space range
 * @param vlasov_pkpm_moms Input array of parallel-kinetic-perpendicular-moment kinetic moments [rho, p_parallel, p_perp]
 * @param euler_pkpm Input array of parallel-kinetic-perpendicular-moment fluid variables [rho ux, rho uy, rho uz]
 * @param p_ij Input pressure tensor p_ij = (p_par - p_perp) b_i b_j + p_perp g_ij
 * @param prim Input array of primitive moments [ux, uy, uz, 1/rho*div(p_par b), T_perp/m, m/T_perp]
 * @param pkpm_accel Input arrary of pkpm acceleration variables ordered as:
 *        0: p_perp_div_b (p_perp/rho*div(b) = T_perp/m*div(b))
          1: bb_grad_u (bb : grad(u))
          2: p_force (total pressure forces in kinetic equation 1/rho div(p_parallel b_hat) - T_perp/m*div(b)
          3: p_perp_source (pressure source for higher Laguerre moments -> bb : grad(u) - div(u) - 2 nu)
 * @param fluid_io Output array of conserved fluid variables (10 components)
 * @param pkpm_vars_io Output array of pkpm variables, primitive and acceleration (8 components)
 */
void gkyl_dg_calc_pkpm_vars_io(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_array* vlasov_pkpm_moms, 
  const struct gkyl_array* euler_pkpm, const struct gkyl_array* p_ij, 
  const struct gkyl_array* prim, const struct gkyl_array* pkpm_accel, 
  struct gkyl_array* fluid_io, struct gkyl_array* pkpm_vars_io);

/**
 * Delete pointer to updater to compute pkpm variables.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_calc_pkpm_vars_release(struct gkyl_dg_calc_pkpm_vars *up);

/**
 * Host-side wrappers for pkpm vars operations on device
 */

void gkyl_dg_calc_pkpm_vars_advance_cu(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm, 
  const struct gkyl_array* pkpm_div_ppar, 
  struct gkyl_array* cell_avg_prim, struct gkyl_array* prim);

void gkyl_dg_calc_pkpm_vars_surf_advance_cu(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm, 
  const struct gkyl_array* p_ij, const struct gkyl_array* cell_avg_prim, 
  struct gkyl_array* prim_surf);

void gkyl_dg_calc_pkpm_vars_pressure_cu(struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range, 
  const struct gkyl_array* bvar, const struct gkyl_array* vlasov_pkpm_moms, struct gkyl_array* p_ij);

void gkyl_dg_calc_pkpm_vars_accel_cu(struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range, 
  const struct gkyl_array* prim_surf, const struct gkyl_array* prim, 
  const struct gkyl_array* bvar, const struct gkyl_array* div_b, const struct gkyl_array* nu, 
  struct gkyl_array* pkpm_lax, struct gkyl_array* pkpm_accel);

void gkyl_dg_calc_pkpm_integrated_vars_cu(struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm, 
  const struct gkyl_array* prim, struct gkyl_array* pkpm_int_vars);

void gkyl_dg_calc_pkpm_vars_source_cu(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_array* qmem, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm,
  struct gkyl_array* rhs);

void gkyl_dg_calc_pkpm_vars_io_cu(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_array* vlasov_pkpm_moms, 
  const struct gkyl_array* euler_pkpm, const struct gkyl_array* p_ij, 
  const struct gkyl_array* prim, const struct gkyl_array* pkpm_accel, 
  struct gkyl_array* fluid_io, struct gkyl_array* pkpm_vars_io);
