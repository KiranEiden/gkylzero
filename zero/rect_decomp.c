#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>

#include <string.h>

static void
rect_decomp_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_rect_decomp *decomp = container_of(ref, struct gkyl_rect_decomp, ref_count);
  gkyl_free(decomp->ranges);
  gkyl_free(decomp);
}    

struct gkyl_rect_decomp*
gkyl_rect_decomp_new_from_cuts(int ndim, const int cuts[], const struct gkyl_range *range)
{
  struct gkyl_rect_decomp *decomp = gkyl_malloc(sizeof(*decomp));

  int ndecomp = 1;
  decomp->ndim = ndim;  

  for (int d=0; d<ndim; ++d) ndecomp *= cuts[d];  
  decomp->ndecomp = ndecomp;
  decomp->ranges = gkyl_malloc(sizeof(struct gkyl_range[ndecomp]));

  memcpy(&decomp->parent_range, range, sizeof(struct gkyl_range));

  div_t qr[GKYL_MAX_DIM];
  for (int d=0; d<ndim; ++d)
    qr[d] = div(gkyl_range_shape(range, d), cuts[d]);

  int *sidx[GKYL_MAX_DIM], *eidx[GKYL_MAX_DIM];
  for (int d=0; d<ndim; ++d) {
    sidx[d] = gkyl_malloc(sizeof(int[cuts[d]]));
    eidx[d] = gkyl_malloc(sizeof(int[cuts[d]]));

    int *shape = gkyl_malloc(sizeof(int[cuts[d]]));

    // compute shape in direction 'd'
    for (int i=0; i<cuts[d]; ++i)
      shape[i] = i<qr[d].rem ? qr[d].quot+1 : qr[d].quot;

    sidx[d][0] = range->lower[d];
    eidx[d][0] = sidx[d][0]+shape[0]-1;
    for (int i=1; i<cuts[d]; ++i) {
      sidx[d][i] = eidx[d][i-1]+1;
      eidx[d][i] = sidx[d][i]+shape[i]-1;
    }

    gkyl_free(shape);
  }

  struct gkyl_range rcuts;
  gkyl_range_init_from_shape(&rcuts, ndim, cuts);
  struct gkyl_range_iter citer;
  gkyl_range_iter_init(&citer, &rcuts);

  int dnum = 0;
  // loop over cuts range, constructing each of the sub-ranges in the
  // decomposition
  while( gkyl_range_iter_next(&citer) ) {
    int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];

    for (int d=0; d<ndim; ++d) {
      lower[d] = sidx[d][citer.idx[d]];
      upper[d] = eidx[d][citer.idx[d]];
    }

    gkyl_sub_range_init(&decomp->ranges[dnum++], range, lower, upper);
  }

  for (int d=0; d<ndim; ++d) {
    gkyl_free(sidx[d]);
    gkyl_free(eidx[d]);
  }

  return decomp;
}

bool
gkyl_rect_decomp_check_covering(const struct gkyl_rect_decomp *decomp)
{
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 1, decomp->parent_range.volume);
  gkyl_array_clear(arr, 0.0);

  // following loops over each sub-range and increments the region it
  // indexes in 'arr'. Each index should be visited exactly once.
  for (int i=0; i<decomp->ndecomp; ++i) {
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &decomp->ranges[i]);

    while (gkyl_range_iter_next(&iter)) {
      double *d = gkyl_array_fetch(arr, gkyl_range_idx(&decomp->ranges[i], iter.idx));
      d[0] += 1.0;
    }
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &decomp->parent_range);
  while (gkyl_range_iter_next(&iter)) {
    const double *d = gkyl_array_cfetch(arr, gkyl_range_idx(&decomp->parent_range, iter.idx));
    if (d[0] != 1.0)
      return false;
  }

  gkyl_array_release(arr);
  
  return true;
}

void
gkyl_rect_decomp_release(struct gkyl_rect_decomp *decomp)
{
  gkyl_free(decomp->ranges);
  gkyl_free(decomp);
}

void
gkyl_create_grid_ranges(const struct gkyl_rect_grid *grid,
  const int *nghost, struct gkyl_range *ext_range, struct gkyl_range *range)
{
  int lower_ext[GKYL_MAX_DIM], upper_ext[GKYL_MAX_DIM];
  int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  
  for (int i=0; i<grid->ndim; ++i) {
    lower_ext[i] = 1-nghost[i];
    upper_ext[i] = grid->cells[i]+nghost[i];

    lower[i] = 1;
    upper[i] = grid->cells[i];
  }
  gkyl_range_init(ext_range, grid->ndim, lower_ext, upper_ext);
  gkyl_sub_range_init(range, ext_range, lower, upper);
}

void gkyl_create_ranges(const struct gkyl_range *inrange,
  const int *nghost, struct gkyl_range *ext_range, struct gkyl_range *range)
{
  int lower_ext[GKYL_MAX_DIM], upper_ext[GKYL_MAX_DIM];
  int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  
  for (int i=0; i<inrange->ndim; ++i) {
    lower_ext[i] = inrange->lower[i]-nghost[i];
    upper_ext[i] = inrange->upper[i]+nghost[i];

    lower[i] = inrange->lower[i];
    upper[i] = inrange->upper[i];
  }
  gkyl_range_init(ext_range, inrange->ndim, lower_ext, upper_ext);
  gkyl_sub_range_init(range, ext_range, lower, upper);  
}
