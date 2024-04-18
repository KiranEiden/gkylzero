#pragma once

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>

// Struct containing the pointers to auxiliary fields.
// Specified hamiltonian in *Canonical* cordinates
struct gkyl_dg_canonical_pb_auxfields { 
  const struct gkyl_array *hamil;
  const struct gkyl_array *alpha_surf;
  const struct gkyl_array *sgn_alpha_surf;
  const struct gkyl_array *const_sgn_alpha;
};

/**
 * Create a new special relativistic Vlasov equation object.
 *
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis functions
 * @param conf_range Configuration space range for use in indexing EM field
 * @param vel_range Velocity space range for use in indexing p/gamma (velocity)
 * @param field_id enum to determine what type of EM fields
 * (special relativistic Vlasov-Maxwell vs. special relativistic neutrals)
 * @param use_gpu bool to determine if on GPU
 * @return Pointer to special relativistic Vlasov equation object
 */
struct gkyl_dg_eqn* gkyl_dg_canonical_pb_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* phase_range, bool use_gpu);

/**
 * Create a new special relativistic Vlasov equation object that lives on NV-GPU
 *
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis functions
 * @param phase_range Phase space range
 * (special relativistic Vlasov-Maxwell vs. special relativistic neutrals)
 * @return Pointer to special relativistic Vlasov equation object
 */
struct gkyl_dg_eqn* gkyl_dg_canonical_pb_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* phase_range);

/**
 * Set the auxiliary fields (e.g. q/m*EM) needed in updating the force terms.
 * 
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_canonical_pb_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_canonical_pb_auxfields auxin);


#ifdef GKYL_HAVE_CUDA
/**
 * CUDA device function to set auxiliary fields (e.g. q/m*EM) needed in updating the force terms.
 * 
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_canonical_pb_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_canonical_pb_auxfields auxin);

#endif
