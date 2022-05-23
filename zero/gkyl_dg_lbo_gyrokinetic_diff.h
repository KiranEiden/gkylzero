#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_dg_lbo_gyrokinetic_diff_auxfields { 
  const struct gkyl_array *bmag_inv;
  const struct gkyl_array *nuSum;
  const struct gkyl_array *nuUSum;
  const struct gkyl_array *nuVtSqSum;
};

/**
 * Create a new gyrokinetic LBO diffusion term equation object.
 *
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis functions
 * @param conf_range Configuration space range for use in indexing primitive moments
 * @param mass Species mass
 * @return Pointer to LBO equation object
 */
struct gkyl_dg_eqn* gkyl_dg_lbo_gyrokinetic_diff_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, double mass);

/**
 * Create a new LBO equation object that lives on NV-GPU
 */
struct gkyl_dg_eqn* gkyl_dg_lbo_gyrokinetic_diff_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, double mass);

/**
 * Set auxiliary fields needed in updating the diffusion flux term.
 * These are bmag, nu, nu*u, and nu*vt^2.
 * 
 * @param eqn Equation pointer
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_lbo_gyrokinetic_diff_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_lbo_gyrokinetic_diff_auxfields auxin);

#ifdef GKYL_HAVE_CUDA

/**
 * CUDA device function to set auxiliary fields needed in updating the diffusion flux term.
 */
void gkyl_lbo_gyrokinetic_diff_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_lbo_gyrokinetic_diff_auxfields auxin);

#endif
