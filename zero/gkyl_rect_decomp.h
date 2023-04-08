#pragma once

#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_ref_count.h>

struct gkyl_rect_decomp {
  int ndim; // dimension of decomposition
  int ndecomp; // number of sub-domains
  struct gkyl_range parent_range; // range that was decomposed
  struct gkyl_range *ranges; // decomposed ranges

  struct gkyl_ref_count ref_count;
};

/**
 * Create a new decomposition of @a range, given @a cuts in each
 * direction. The total number of decomposed ranges are product of all
 * cuts. The decomposed ranges are sub-ranges of the parent range.
 *
 * @param ndim Number of dimensions
 * @param cuts Cuts in each direction.
 * @param range Range to decompose
 * @return Decomposition of @a range
 */
struct gkyl_rect_decomp *gkyl_rect_decomp_new_from_cuts(int ndim,
  const int cuts[], const struct gkyl_range *range);

/**
 * Acquire a pointer to the decomposition.
 *
 * @param decomp Decom to acquire pointer to
 * @return New decomposition
 */
struct gkyl_rect_decomp* gkyl_rect_decomp_acquire(const struct gkyl_rect_decomp *decomp);

/**
 * Check if decomposition is  a valid covering of the range.
 *
 * NOTE: This function internally allocates memory over the complete
 * parent range. This can be a problem if the parent range is huge.
 *
 * @param decomp Demposition to check
 * @return true if this is a valid covering
 */
bool gkyl_rect_decomp_check_covering(const struct gkyl_rect_decomp *decomp);

/**
 * Free decomposition.
 *
 * @param decomp Decomposition to free
 */
void gkyl_rect_decomp_release(struct gkyl_rect_decomp *decomp);

// The functions below are utility functions to construct properly
// nested ranges that extend over the grid or over local ranges, given
// ghost cells.

/**
 * Create range and extended ranges from grid and ghost-cell data. The
 * range is a sub-range of the extended range.
 *
 * @param grid Grid to compute ranges for
 * @param nghost Number of ghost-cells in each direction
 * @param ext_range On output, extended range spanning grid+ghost-cells
 * @param range On output, range spanning grid. Sub-range of ext_range.
 */
void gkyl_create_grid_ranges(const struct gkyl_rect_grid *grid,
  const int *nghost, struct gkyl_range *ext_range,
  struct gkyl_range *range);

/**
 * Create range and extended ranges from give range and ghost-cell
 * data. The range is a sub-range of the extended range.
 *
 * @param inrange Input range to use
 * @param nghost Number of ghost-cells in each direction
 * @param ext_range On output, extended range spanning inrange+ghost-cells
 * @param range On output, range same as inrange, but sub-range of ext_range.
 */
void gkyl_create_ranges(const struct gkyl_range *inrange,
  const int *nghost, struct gkyl_range *ext_range, struct gkyl_range *range);
