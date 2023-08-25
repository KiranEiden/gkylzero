#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_basis.h>

// Object type
typedef struct gkyl_dg_calc_em_vars gkyl_dg_calc_em_vars;

/**
 * Create new updater to compute EM variables needed in 
 * updates and used for diagnostics. Supports:
 * 1. bvar = [b_i (3 components), b_i b_j (6 components)],
 * the magnetic field unit vector and unit tensor, b_i = B_i/|B|
 * 2. ExB = E x B/|B|^2, the E x B velocity 
 * Free using gkyl_dg_calc_em_vars_release.
 *
 * @param conf_grid Configuration space grid (for getting cell spacing and cell center)
 * @param cbasis Configuration space basis functions
 * @param mem_range Configuration space range to compute variables over
 * Note: This range sets the size of the bin_op memory and thus sets the
 * range over which the updater loops for the batched linear solves
 * @param is_ExB bool to determine if updater is for computing E x B velocity
 * @param use_gpu bool to determine if on GPU
 * @return New updater pointer.
 */
struct gkyl_dg_calc_em_vars* 
gkyl_dg_calc_em_vars_new(const struct gkyl_rect_grid *conf_grid, 
  const struct gkyl_basis* cbasis, const struct gkyl_range *mem_range, 
  bool is_ExB, bool use_gpu);

/**
 * Create new updater to compute EM variables on
 * NV-GPU. See new() method for documentation.
 */
struct gkyl_dg_calc_em_vars* 
gkyl_dg_calc_em_vars_cu_dev_new(const struct gkyl_rect_grid *conf_grid, 
  const struct gkyl_basis* cbasis, const struct gkyl_range *conf_range, 
  bool is_ExB);

/**
 * Compute either
 * 1. The magnetic field unit vector (b_i = B_i/|B|, three components) and unit tensor (b_i b_j = B_i B_j/|B|^2, 6 components) 
 * 2. The ExB velocity (E x B/|B|^2, three components)
 * Note order of operations is designed to minimize aliasing errors
 * 1. Compute B_i B_j or numerator (E x B)_i and denominator (|B|^2) using weak multiplication 
 * 2. Compute unit tensor (b_i b_j = B_i B_j/|B|^2, 6 components) or (E x B/|B|^2) using either
 *    basis_inv operator (for p=1) or weak division (p>1)
 * 3. For bvar, project diagonal components of bb onto quadrature points, evaluate square root point wise, 
 *    and project back onto modal basis using basis_sqrt to obtain b_i (see gkyl_basis_*_sqrt.h in kernels/basis/)
 *
 * @param up             Updater for computing EM variables (contains range, pre-allocated memory, and pointers to kernels)
 * @param em             Input array which contain EM fields (Ex, Ey, Ez, Bx, By, Bz)
 * @param cell_avg_magB2 Array for storing boolean value of whether |B|^2 uses *only* cell averages 
 *                       to minimize positivity violations (default: false)
 * @param out            Output array of either magnetic field unit vector and unit tensor or E x B velocity
 */
void gkyl_dg_calc_em_vars_advance(struct gkyl_dg_calc_em_vars *up, 
  const struct gkyl_array* em, struct gkyl_array* cell_avg_magB2, struct gkyl_array* out);

/**
 * Compute surface expansion of bvar
 *
 * @param up Updater for computing pkpm variables 
 * @param em                  Input array which contain EM fields (Ex, Ey, Ez, Bx, By, Bz)
 * @param cell_avg_magB2_surf Array for storing boolean value of whether |B|^2 on the surface 
 *                            uses *only* cell averages (2*cdim components: xl, xr, yl, yr, zl, & zr)
 *                            to minimize positivity violations (default: false)
 * @param bvar_surf           Output array of surface expansions of bvar
 */
void gkyl_dg_calc_em_vars_surf_advance(struct gkyl_dg_calc_em_vars *up, 
  const struct gkyl_array* em, struct gkyl_array* cell_avg_magB2_surf, struct gkyl_array* bvar_surf);

/**
 * Compute div(b) and max(|b_i|) penalization
 *
 * @param up Updater for computing pkpm variables 
 * @param conf_range Configuration space range
 * @param bvar_surf Input array of surface expansions of bvar
 * @param bvar Input array of volume expansion of bvar
 * @param max_b Output array of max(|b_i|) penalization
 * @param div_b Output array of div(b)
 */
void gkyl_dg_calc_em_vars_div_b(struct gkyl_dg_calc_em_vars *up, const struct gkyl_range *conf_range, 
  const struct gkyl_array* bvar_surf, const struct gkyl_array* bvar, 
  struct gkyl_array* max_b, struct gkyl_array* div_b);

/**
 * Delete pointer to updater to compute EM variables.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_calc_em_vars_release(struct gkyl_dg_calc_em_vars *up);

/**
 * Host-side wrappers for em vars operations on device
 */

void gkyl_dg_calc_em_vars_advance_cu(struct gkyl_dg_calc_em_vars *up, 
  const struct gkyl_array* em, struct gkyl_array* cell_avg_magB2, struct gkyl_array* out);

void gkyl_dg_calc_em_vars_surf_advance_cu(struct gkyl_dg_calc_em_vars *up, 
  const struct gkyl_array* bvar, struct gkyl_array* bvar_surf);

void gkyl_dg_calc_em_vars_div_b_cu(struct gkyl_dg_calc_em_vars *up, const struct gkyl_range *conf_range, 
  const struct gkyl_array* bvar_surf, const struct gkyl_array* bvar, 
  struct gkyl_array* max_b, struct gkyl_array* div_b);
