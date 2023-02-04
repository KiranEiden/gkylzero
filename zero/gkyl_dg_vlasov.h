#pragma once

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_dg_vlasov_auxfields { 
  const struct gkyl_array *field; // q/m*(E,B) for Maxwell's, q/m*phi for Poisson's (gradient calculated in kernel)
  const struct gkyl_array *ext_field; // constant q/m*A for Poisson's (curl calculated in kernel)
  const struct gkyl_array *cot_vec; // cotangent vectors (e^i) used in volume term if general geometry enabled
  const struct gkyl_array *alpha_geo; // alpha^i (e^i . alpha) used in surface term if general geometry enabled
};

/**
 * Create a new Vlasov equation object.
 *
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis functions
 * @param conf_range Configuration space range for use in indexing EM field
 * @param phase_range Phase space range for indexing geometry quantities (alpha_geo)
 * @param model_id enum to determine what type of Vlasov model (Cartesian vs. general geometry)
 * @param field_id enum to determine what type of EM fields (Vlasov-Maxwell vs. Vlasov-Poisson vs. neutrals)
 * @param use_gpu bool to determine if on GPU
 * @return Pointer to Vlasov equation object
 */
struct gkyl_dg_eqn* gkyl_dg_vlasov_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range, const struct gkyl_range* phase_range,
  enum gkyl_model_id model_id, enum gkyl_field_id field_id, bool use_gpu);

/**
 * Create a new Vlasov equation object that lives on NV-GPU
 *
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis functions
 * @param conf_range Configuration space range for use in indexing EM field
 * @param phase_range Phase space range for indexing geometry quantities (alpha_geo)
 * @param model_id enum to determine what type of Vlasov model (Cartesian vs. general geometry)
 * @param field_id enum to determine what type of EM fields (Vlasov-Maxwell vs. Vlasov-Poisson vs. neutrals)
 * @return Pointer to Vlasov equation object
 */
struct gkyl_dg_eqn* gkyl_dg_vlasov_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range, const struct gkyl_range* phase_range,
  enum gkyl_model_id model_id, enum gkyl_field_id field_id);

/**
 * Set the auxiliary fields (e.g. q/m*EM) needed in updating the force terms.
 * 
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_vlasov_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_vlasov_auxfields auxin);

#ifdef GKYL_HAVE_CUDA
/**
 * CUDA device function to set auxiliary fields (e.g. q/m*EM) needed in updating the force terms.
 * 
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_vlasov_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_vlasov_auxfields auxin);


#endif
