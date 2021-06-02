#pragma once

#include <stdio.h>
#include <gkyl_util.h>

/**
 * Rectangular grid object.
 */
struct gkyl_rect_grid {
  int ndim; // number of dimensions
  double lower[GKYL_MAX_DIM]; // lower-left corner
  double upper[GKYL_MAX_DIM]; // upper-right corner
  int cells[GKYL_MAX_DIM]; // number of cells    
  double dx[GKYL_MAX_DIM]; // cell spacing
  double cellVolume; // cell volume
};

/**
 * Create new grid object.
 *
 * @param grid Grid object to initialize.
 * @param ndim Dimension of grid
 * @param lower Coordinates of lower-left corner of grid
 * @param upper Coordinates of upper-right corner of grid
 * @param cells Number of cells in each direction
 */
void gkyl_rect_grid_init(struct gkyl_rect_grid *grid, int ndim,
  const double *lower, const double *upper, const int *cells);

/**
 * Clone grid object on NV-GPU.
 *
 * @param grid Grid object on device to clone
 * @return Clone valid on device.
 */

struct gkyl_rect_grid* gkyl_rect_grid_clone_on_cu_dev(struct gkyl_rect_grid* grid);

/**
 * Get cell-center coordinates. idx is the zero-based cell index.
 *
 * @param grid Grid object
 * @param idx Index of cell (lower-left corner has all index (0,0,...) )
 * @param xc On output, cell-center coordinates of cell 'idx'
 */
GKYL_CU_DH
static inline void gkyl_rect_grid_cell_center(const struct gkyl_rect_grid *grid,
  const int *idx, double *xc)
{
  for (int i=0; i<grid->ndim; ++i)
    xc[i] = grid->lower[i]+(idx[i]+0.5)*grid->dx[i];
}

/**
 * Write grid data to file. File must be opened by caller of this
 * function. Data is written in binary format.
 *
 * @param grid Grid object to write
 * @param fp File handle to write to.
 */
void gkyl_rect_grid_write(const struct gkyl_rect_grid *grid, FILE *fp);
