#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_dg_updater_fluid gkyl_dg_updater_fluid;

// return type for drag and diffusion timers
struct gkyl_dg_updater_fluid_tm {
  double fluid_tm; // time for fluid updates
};

/**
 * Create new updater to update fluid equations using hyper dg.
 * Supports advection, isothermal Euler, and Euler
 *
 * @param grid Grid object
 * @param cbasis Configuration space basis functions
 * @param conf_range Config space range
 * @param eqn_id Enum identifier for equation type (see gkyl_eqn_type.h)
 * @param param parameter for fluid equation specification 
 * (e.g. vt for isothermal Euler, gas_gammas for Euler)
 * @param use_gpu Boolean to determine whether struct objects are on host or device
 * 
 * @return New fluid updater object
 */
gkyl_dg_updater_fluid* gkyl_dg_updater_fluid_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis *cbasis, const struct gkyl_range *conf_range, 
  enum gkyl_eqn_type eqn_id, double param, bool use_gpu);

/**
 * Acquire fluid equation object
 *
 * @param fluid fluid updater object
 * 
 * @return fluid equation object
 */
struct gkyl_dg_eqn* 
gkyl_dg_updater_fluid_acquire_eqn(const gkyl_dg_updater_fluid* fluid);

/**
 * Compute RHS of DG update. The update_rng MUST be a sub-range of the
 * range on which the array is defined. That is, it must be either the
 * same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param fluid fluid updater object
 * @param eqn_id Enum identifier for equation type (see gkyl_eqn_type.h)
 * @param update_rng Range on which to compute.
 * @param aux1 Auxiliary field 1 (always u_i)
 * @param aux2 Auxiliary field 2 (p_ij or div_p)
 * @param fIn Input to updater
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T])
 * @param rhs RHS output
 */
void gkyl_dg_updater_fluid_advance(gkyl_dg_updater_fluid *fluid,
  enum gkyl_eqn_type eqn_id, const struct gkyl_range *update_rng,
  const struct gkyl_array *aux1, const struct gkyl_array *aux2, 
  const struct gkyl_array *aux3, 
  const struct gkyl_array* GKYL_RESTRICT fIn,
  struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs);

void gkyl_dg_updater_fluid_advance_cu(gkyl_dg_updater_fluid *fluid,
  enum gkyl_eqn_type eqn_id, const struct gkyl_range *update_rng,
  const struct gkyl_array *aux1, const struct gkyl_array *aux2, 
  const struct gkyl_array *aux3, 
  const struct gkyl_array* GKYL_RESTRICT fIn,
  struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs);

/**
 * Return total time spent in drag and diffusion terms
 *
 * @param fluid Updater object
 * @return timers
 */
struct gkyl_dg_updater_fluid_tm gkyl_dg_updater_fluid_get_tm(const gkyl_dg_updater_fluid *fluid);

/**
 * Delete updater.
 *
 * @param fluid Updater to delete.
 */
void gkyl_dg_updater_fluid_release(gkyl_dg_updater_fluid* fluid);