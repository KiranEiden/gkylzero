#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_alloc.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_array_ops.h>
#include <gkyl_mat.h>
#include <gkyl_mat_triples.h>
#include <gkyl_superlu_ops.h>

// Object type
typedef struct gkyl_fem_parproj gkyl_fem_parproj;

// Boundary condition types.
enum gkyl_fem_parproj_bc_type {
  GKYL_FEM_PARPROJ_PERIODIC = 0,
  GKYL_FEM_PARPROJ_DIRICHLET, // sets the value.
  GKYL_FEM_PARPROJ_NONE,      // does not enforce a BC.
};

/**
 * Create new updater to project a DG field onto the FEM (nodal) basis
 * in order to make the field continuous or, thanks to the option to pass
 * a multiplicative weight, solve 1D algebraic equations in which the output
 * field is continuous (but the input may not be). That is, we solve
 *    wgt*phi_{fem} \doteq rho_{dg}
 * where wgt is the weight field, phi_{fem} is the (continuous field)
 * we wish to compute, rho_{dg} is the (discontinuous) input source field,
 * and \doteq implies weak equality with respect to the FEM basis.
 * Free using gkyl_fem_parproj_release method.
 *
 * @param grid Grid object
 * @param basis Basis functions of the DG field.
 * @param bctype Type of boundary condition (see gkyl_fem_parproj_bc_type).
 * @param weight multiplicative weight on left-side of the operator.
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
struct gkyl_fem_parproj* gkyl_fem_parproj_new(
  const struct gkyl_rect_grid *grid, const struct gkyl_basis *basis,
  enum gkyl_fem_parproj_bc_type bctype, const struct gkyl_array *weight,
  bool use_gpu);

/**
 * Set the multiplicative weight.
 * Note that this triggers a re-building and decomposition
 * of the LHS matrix.
 *
 * @param up FEM project updater to run.
 * @param tm Time at which projection must be computed
 * @param inw Input array weight.
 */
void gkyl_fem_parproj_set_weight(const struct gkyl_fem_parproj *up,
  const struct gkyl_array *inw);

/**
 * Begin assembling the right-side source vector and, if necessary,
 * the stiffness matrix.
 * In parallel simulations there will be an option to use
 * non-blocking MPI here and later checking that assembly is complete.
 *
 * @param up FEM project updater to run.
 * @param tm Time at which projection must be computed
 * @param src Input source field.
 */
void gkyl_fem_parproj_begin_assembly(const struct gkyl_fem_parproj *up,
  double tm, const struct gkyl_array *src);

/**
 * Assign the right-side vector with the discontinuous (DG) source field.
 *
 * @param up FEM project updater to run.
 * @param rhsin DG field to set as RHS source.
 * @param phibc Potential to use for Dirichlet BCs (only use ghost cells).
 */
void gkyl_fem_parproj_set_rhs(struct gkyl_fem_parproj* up, const struct gkyl_array *rhsin, const struct gkyl_array *phibc);

void gkyl_fem_parproj_set_rhs_cu(struct gkyl_fem_parproj *up, const struct gkyl_array *rhsin, const struct gkyl_array *phibc);

/**
 * Solve the linear problem.
 *
 * @param up FEM project updater to run.
 */
void gkyl_fem_parproj_solve(struct gkyl_fem_parproj* up, struct gkyl_array *phiout);

void gkyl_fem_parproj_solve_cu(struct gkyl_fem_parproj* up, struct gkyl_array *phiout);

/**
 * Compute the projection onto the FEM basis. 
 *
 * @param up FEM project updater to run
 * @param tm Time at which projection must be computed
 * @param out Output array
 */
void gkyl_fem_parproj_advance(struct gkyl_fem_parproj *up,
  double tm, struct gkyl_array *out);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_fem_parproj_release(struct gkyl_fem_parproj *up);
