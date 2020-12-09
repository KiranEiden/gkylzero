#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_mom_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_mom_calc gkyl_mom_calc;

/**
 * Create new updater to compute moments of distribution
 * function. Free using gkyl_mom_calc_new_release.
 *
 * @param grid Grid object
 * @param momt Pointer to moment type object
 * @return New updater pointer.
 */
gkyl_mom_calc* gkyl_mom_calc_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_mom_type *momt);

/**
 * Compute moment of distribution function. The phase_rng and conf_rng
 * MUST be a sub-ranges of the range on which the distribution
 * function and the moments are defined. These ranges must be
 * self-consistently constructed.
 *
 * @param calc Moment calculator updater to run
 * @param phase_rng Phase-space range
 * @param conf_rng Config-space rnage
 * @param fin Input distribution function array
 * @param mout Output moment array
 */
void gkyl_mom_calc_advance(const gkyl_mom_calc* calc,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng,
  const struct gkyl_array *fin, struct gkyl_array *mout);

/**
 * Delete pointer to moment calculator updater.
 *
 * @param calc Updater to delete.
 */
void gkyl_mom_calc_release(gkyl_mom_calc* calc);