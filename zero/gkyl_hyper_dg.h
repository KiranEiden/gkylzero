#pragma once

#include <gkyl_util.h>
#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_hyper_dg gkyl_hyper_dg;

/**
 * Create new updater to update equations using DG algorithm.
 *
 * @param grid Grid object
 * @param basis Basis functions
 * @param equation Equation object
 * @param num_up_dirs Number of directions to update
 * @param update_dirs List of directions to update (size 'num_up_dirs')
 * @param zero_flux_flags Flags to indicate if direction has zero-flux BCs
 * @param update_vol_term Set to 0 to skip volume update
 */
gkyl_hyper_dg* gkyl_hyper_dg_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *basis, const struct gkyl_dg_eqn *equation,
  int num_up_dirs, int update_dirs[GKYL_MAX_DIM], int zero_flux_flags[GKYL_MAX_DIM],
  int update_vol_term, bool use_gpu);

/**
 * Create new updater on CUDA device to update equations using DG algorithm.
 *
 * @param grid_cu Grid object (on device)
 * @param basis Basis functions
 * @param equation Equation object
 * @param num_up_dirs Number of directions to update
 * @param update_dirs List of directions to update (size 'num_up_dirs')
 * @param zero_flux_flags Flags to indicate if direction has zero-flux BCs
 * @param update_vol_term Set to 0 to skip volume update
 */
gkyl_hyper_dg* gkyl_hyper_dg_cu_dev_new(const struct gkyl_rect_grid *grid_cu,
  const struct gkyl_basis *basis, const struct gkyl_dg_eqn *equation,
  int num_up_dirs, int update_dirs[GKYL_MAX_DIM], int zero_flux_flags[GKYL_MAX_DIM],
  int update_vol_term);

/**
 * Compute RHS of DG update. The update_rng MUST be a sub-range of the
 * range on which the array is defined. That is, it must be either the
 * same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param hdg Hyper DG updater object
 * @param update_rng Range on which to compute.
 * @param fIn Input to updater
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T])
 * @param rhs RHS output
 */
void gkyl_hyper_dg_advance(gkyl_hyper_dg *hdg, const struct gkyl_range *update_rng,
  const struct gkyl_array *fIn, struct gkyl_array *cflrate, struct gkyl_array *rhs);

// CUDA call
void gkyl_hyper_dg_advance_cu(gkyl_hyper_dg* hdg, const struct gkyl_range *update_range,
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT cflrate,
  struct gkyl_array* GKYL_RESTRICT rhs);

/**
 * Set if volume term should be computed or not.
 *
 * @param hdg Hyper DG updater object
 * @param update_vol_term Set to 1 to update vol term, 0 otherwise
 */
void gkyl_hyper_dg_set_update_vol(gkyl_hyper_dg *hdg, int update_vol_term);
// On-device version
void gkyl_hyper_dg_set_update_vol_cu(gkyl_hyper_dg *hdg, int update_vol_term);
  
/**
 * Delete updater.
 *
 * @param hdg Updater to delete.
 */
void gkyl_hyper_dg_release(gkyl_hyper_dg* hdg);
