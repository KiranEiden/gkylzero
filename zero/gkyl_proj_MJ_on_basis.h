#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_proj_MJ_on_basis gkyl_proj_MJ_on_basis;

/**
 * Create new updater to project MJ on basis functions. Free
 * using gkyl_proj_MJ_on_basis_release method.
 *
 * @param grid Grid object
 * @param conf_basis Conf-space basis functions
 * @param phase_basis Phase-space basis functions
 * @param num_quad Number of quadrature nodes
 * @param mass, mass of the species
 * @return New updater pointer.
 */
gkyl_proj_MJ_on_basis* gkyl_proj_MJ_on_basis_new(
  const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis,
  int num_quad, double mass);

/**
 * Compute projection of MJ on basis. This method takes
 * lab-frame moments to compute the projection of MJ on basis
 * functions.
 *
 * @param pob Project on basis updater to run
 * @param phase_rng Phase-space range
 * @param conf_rng Config-space range
 * @param M0 Number density moment
 * @param M1i Momentum in lab-frame
 * @param grid Grid object
 * @param conf_basis Conf-space basis functions
 * @param phase_basis Phase-space basis functions
 * @param M2 Temperature in the fluid rest frame

 * @param f_MJ Output MJ
 */
void gkyl_proj_MJ_on_basis_lab_mom(const gkyl_proj_MJ_on_basis *pob,
  const struct gkyl_range phase_range, const struct gkyl_range conf_range,
  const struct gkyl_array *M0, const struct gkyl_array *M1i, const struct gkyl_array *M2,
  const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis,
  struct gkyl_array *f_MJ);

/**
 * Compute projection of MJ on basis. This method takes
 * primitive (fluid-frame) moments to compute the projection of
 * MJ on basis functions.
 *
 * @param pob Project on basis updater to run
 * @param phase_rng Phase-space range
 * @param conf_rng Config-space range
 * @param n Number density moment
 * @param vel Velocity vector
 * @param T Temperature in the labframe
 * @param f_MJ Output MJ
 */
void gkyl_proj_MJ_on_basis_prim_mom(const gkyl_proj_MJ_on_basis *pob,
  const struct gkyl_range phase_range, const struct gkyl_range conf_range,
  const struct gkyl_array *n, const struct gkyl_array *vel, const struct gkyl_array *T,
  struct gkyl_array *f_MJ);

/**
 * Delete updater.
 *
 * @param pob Updater to delete.
 */
void gkyl_proj_MJ_on_basis_release(gkyl_proj_MJ_on_basis* pob);
