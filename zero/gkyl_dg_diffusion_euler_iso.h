#pragma once

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_dg_diffusion_e_i_auxfields {
  const struct gkyl_array* D;
  const struct gkyl_array* u_i;
};

/**
 * Create a new diffusion equation object.
 *
 * @param Basis functions
 * @param range Range for use in indexing diffusion tensor
 * @return Pointer to diffusion equation object
 */
struct gkyl_dg_eqn* gkyl_dg_diffusion_e_i_new(const struct gkyl_basis* cbasis, const struct gkyl_range* conf_range, bool use_gpu);

/**
 * Create a new diffusion equation object that lives on NV-GPU
 *
 * @param cbasis Configuration space basis functions
 * @param conf_range Configuration space range for use in indexing diffusion tensor
 * @return Pointer to diffusion equation object
 */
struct gkyl_dg_eqn* gkyl_dg_diffusion_e_i_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_range* conf_range);

/**
 * Set the auxiliary fields (e.g. diffusion tensor D) needed in updating diffusion equation.
 *
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_diffusion_e_i_set_auxfields(const struct gkyl_dg_eqn* eqn, struct gkyl_dg_diffusion_e_i_auxfields auxin);

#ifdef GKYL_HAVE_CUDA

/**
 * CUDA device function to set auxiliary fields (e.g. diffusion tensor D) needed in updating diffusion equation.
 *
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_diffusion_e_i_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_diffusion_e_i_auxfields auxin);

#endif
