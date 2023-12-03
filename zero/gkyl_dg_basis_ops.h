#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>

// Type for storing preallocating memory needed in various operations
typedef struct gkyl_dg_basis_op_mem gkyl_dg_basis_op_mem;

/**
 * Given values and gradients at the corner of a 1D cell, compute the
 * expansion cofficients of modal cubic (p=3) basis functions. The
 * output @a coeff must be preallocated and have at least 4 elements.
 *
 * @param val val[0] and val[1] are the values at left/right edges
 * @param grad grad[0] and grad[1] are the gradients at left/right edges
 * @param coeff On output, the DG expansion coefficients for p=3 basis.
 */
void gkyl_dg_calc_cubic_1d(const double val[2], const double grad[2], double *coeff);

/**
 * Given values and gradients at the corners of a 2D cell, compute the
 * expansion cofficients of modal cubic (p=3) tensor basis functions. The
 * output @a coeff must be preallocated and have at least 16 elements. The nodes
 * are numbered as follows.
 *
 *  2       4
 *  *-------*
 *  |       |
 *  |       |
 *  |       |
 *  *-------*
 *  1       3
 *
 * @param val val[i] is the value at node i
 * @param gradx gradx[i] is the x-derivative at node i
 * @param grady grady[i] is the y-derivative at node i
 * @param gradxy gradxy[i] is the xy- (cross) derivative at node i
 * @param coeff On output, the DG expansion coefficients for p=3 tensor basis.
 */
void gkyl_dg_calc_cubic_2d(const double val[4],
  const double gradx[4], const double grady[4], const double gradxy[4],
  double *coeff);

/**
 * Allocate memory for use in the computing 1D cubic reconstruction
 * from nodal data. Free memory using gkyl_dg_basis_op_mem_release method.
 *
 * @param cells Number of cells in grid
 * @return Newly alloicated memory.
 */
gkyl_dg_basis_op_mem *gkyl_dg_alloc_cubic_1d(int cells);

/**
 * Allocate memory for use in the computing 2D cubic reconstruction
 * from nodal data. Free memory using gkyl_dg_basis_op_mem_release method.
 *
 * @param cells Number of cells in each direction in 2D grid
 * @return Newly alloicated memory.
 */
gkyl_dg_basis_op_mem *gkyl_dg_alloc_cubic_2d(int cells[2]);

/**
 * Release memory allocated for use in DG basis ops.
 *
 * @param mem Memory to release
 */
void gkyl_dg_basis_op_mem_release(gkyl_dg_basis_op_mem *mem);

/**
 * Compute cubic expansion from 1D nodal values. Note that the
 * nodal_vals and cubic arrays must not have any ghost cells and must
 * be properly allocated by the caller. Further, mem must also be
 * constructed from a call to gkyl_dg_alloc_cubic_1d.
 *
 * @param mem Memory space to use
 * @param cells Number of cells in grid
 * @param dx Cell spacing
 * @param nodal_vals Array holding nodal values
 * @param cubic On output, DG expansions of cubic
 */
void gkyl_dg_calc_cubic_1d_from_nodal_vals(gkyl_dg_basis_op_mem *mem, int cells, double dx,
  const struct gkyl_array *nodal_vals, struct gkyl_array *cubic);

/**
 * Compute cubic expansion from 2D nodal values. Note that the
 * nodal_vals and cubic arrays must not have any ghost cells and must
 * be properly allocated by the caller. Further, mem must also be
 * constructed from a call to gkyl_dg_alloc_cubic_1d.
 *
 * @param mem Memory space to use
 * @param cells Number of cells along each direction
 * @param dx Cell spacing in each direction
 * @param nodal_vals Array holding nodal values
 * @param cubic On output, DG expansions of cubic
 */
void gkyl_dg_calc_cubic_2d_from_nodal_vals(gkyl_dg_basis_op_mem *mem, int cells[2], double dx[2],
  const struct gkyl_array *nodal_vals, struct gkyl_array *cubic);
