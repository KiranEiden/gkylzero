#pragma once

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_dg_diffusion_auxfields { 
  const struct gkyl_array* D;
};

/**
 * Create a new diffusion equation object.
 *
 * @param Basis functions
 * @param range Range for use in indexing diffusion tensor
 * @return Pointer to diffusion equation object
 */
struct gkyl_dg_eqn* gkyl_dg_diffusion_new(const struct gkyl_basis* basis, const struct gkyl_range* conf_range, bool use_gpu);

/**
 * Set the auxiliary fields (e.g. advection velocity u) needed in updating advection equation.
 * 
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_diffusion_set_auxfields(const struct gkyl_dg_eqn* eqn, struct gkyl_dg_diffusion_auxfields auxin);
