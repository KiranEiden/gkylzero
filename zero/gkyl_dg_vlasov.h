#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>

/**
 * Create a new Vlasov equation object.
 *
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis functions
 * @param conf_range Configuration space range for use in indexing EM field
 * @return Pointer to Vlasov equation object
 */
struct gkyl_dg_eqn* gkyl_dg_vlasov_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range);

/**
 * Create a new Vlasov equation object that lives on NV-GPU
 *
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis functions
 * @param conf_range Configuration space range for use in indexing EM field
 * @return Pointer to Vlasov equation object
 */
struct gkyl_dg_eqn* gkyl_dg_vlasov_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range);

/**
 * Set the q/m*EM field needed in updating the force terms.
 * 
 * @param eqn Equation pointer
 * @param qmem Pointer to EM field scaled by q/m,
 */
void gkyl_vlasov_set_qmem(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *qmem);

#ifdef GKYL_HAVE_CUDA
/**
 * CUDA device functions to set the q/m*EM field needed in updating the force terms.
 * 
 * @param eqn Equation pointer
 * @param qmem Pointer to EM field scaled by q/m,
 */
void gkyl_vlasov_set_qmem_cu(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *qmem);
__global__ void gkyl_vlasov_set_qmem_cu_kernel(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *qmem);

#endif
