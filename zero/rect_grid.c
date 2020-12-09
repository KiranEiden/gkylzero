#include <stdint.h>

#include <gkyl_rect_grid.h>

void
gkyl_rect_grid_init(struct gkyl_rect_grid *grid, int ndim,
  const double *lower, const double *upper, const int *cells)
{
  grid->ndim = ndim;  
  grid->cellVolume = 1.0;
  for (unsigned i=0; i<ndim; ++i) {
    grid->lower[i] = lower[i];
    grid->upper[i] = upper[i];
    grid->cells[i] = cells[i];
    grid->dx[i] = (upper[i]-lower[i])/cells[i];
    grid->cellVolume *= grid->dx[i];
  }
}

void
gkyl_rect_grid_cell_center(const struct gkyl_rect_grid *grid, const int *idx, double *xc)
{
  for (unsigned i=0; i<grid->ndim; ++i)
    xc[i] = grid->lower[i]+(idx[i]+0.5)*grid->dx[i];
}

void
gkyl_rect_grid_write(const struct gkyl_rect_grid *grid, FILE *fp)
{
  // dimension and shape are written as 64 bit integers
  uint64_t ndim = grid->ndim;
  uint64_t cells[GKYL_MAX_DIM];
  for (unsigned d=0; d<grid->ndim; ++d)
    cells[d] = grid->cells[d];
  
  fwrite(&ndim, sizeof(uint64_t), 1, fp);
  fwrite(cells, sizeof(uint64_t), grid->ndim, fp);
  fwrite(grid->lower, sizeof(double), grid->ndim, fp);
  fwrite(grid->upper, sizeof(double), grid->ndim, fp);
}