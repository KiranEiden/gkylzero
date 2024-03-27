// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_vlasov_lte_moments.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h> 
#include <gkyl_util.h>
#include <assert.h>

GKYL_CU_DH
static inline void
comp_to_phys(int ndim, const double *eta,
  const double * GKYL_RESTRICT dx, const double * GKYL_RESTRICT xc,
  double* GKYL_RESTRICT xout)
{
  for (int d=0; d<ndim; ++d) xout[d] = 0.5*dx[d]*eta[d]+xc[d];
}

static inline void
copy_idx_arrays(int cdim, int pdim, const int *cidx, const int *vidx, int *out)
{
  for (int i=0; i<cdim; ++i)
    out[i] = cidx[i];
  for (int i=cdim; i<pdim; ++i)
    out[i] = vidx[i-cdim];
}

struct gkyl_vlasov_lte_proj_on_basis {
  struct gkyl_rect_grid phase_grid;
  int cdim; // Configuration-space dimension
  int pdim; // Phase-space dimension
  int vdim; // Velocity-space dimension

  struct gkyl_basis conf_basis; // Configuration-space basis
  struct gkyl_basis phase_basis; // Phase-space basis
  struct gkyl_basis* phase_basis_on_dev; // Device-side Phase-space basis for storing
                                         // nodal-to-modal kernels.

  int num_conf_basis; // number of Configuration-space basis functions
  int num_phase_basis; // number of Phase-space basis functions

  bool is_relativistic; // Boolean for if we are projecting the relativistic LTE
  bool use_gpu; // Boolean if we are performing projection on device.

  bool use_quad2m; // Boolean for if we are using special nodal to modal kernels
                   // These kernels are if num_quad = p+1 and there is a one-to-one
                   // transformation of the Gauss-Legendre quadrature nodal basis
                   // and our modal basis.

  struct gkyl_range conf_qrange; // Range of Configuration-space ordinates.
  struct gkyl_range phase_qrange; // Range of Phase-space ordinates.

  // for quadrature in Phase-space
  int tot_quad; // total number of Phase-space quadrature points
  struct gkyl_array *ordinates; // Phase-space ordinates for quadrature
  struct gkyl_array *weights; // Phase-space weights for quadrature
  struct gkyl_array *basis_at_ords; // Phase-space basis functions at ordinates

  // for quadrature in Configuration-space
  int tot_conf_quad; // total number of Configuration-space quadrature points
  struct gkyl_array *conf_ordinates; // Configuration-space ordinates for quadrature
  struct gkyl_array *conf_weights; // weights for Configuration-space quadrature  
  struct gkyl_array *conf_basis_at_ords; // Configuration-space basis functions at ordinates

  struct gkyl_array *fun_at_ords; // function LTE distribution evaluated at
                                  // ordinates in a cell.
  struct gkyl_array *fun_at_ords_on_dev; // device-side function at ordinates

  int *p2c_qidx;  // Mapping between Configuration-space and Phase-space ordinates.

  struct gkyl_vlasov_lte_moments *moments_up; // LTE moment calculation routine for computing density
  struct gkyl_array *num_ratio; // Number density ratio: num_ratio = n_target/n0
  struct gkyl_dg_bin_op_mem *mem; // bin_op memory to compute ratio and rescale distribution function
};
