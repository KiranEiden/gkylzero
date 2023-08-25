/* -*- c++ -*- */

#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_calc_em_vars.h>
#include <gkyl_dg_calc_em_vars_priv.h>
#include <gkyl_util.h>
}

__global__ static void
gkyl_dg_calc_em_vars_set_cu_kernel(gkyl_dg_calc_em_vars* up,
  struct gkyl_nmat *As, struct gkyl_nmat *xs, struct gkyl_range conf_range,
  const struct gkyl_array* em, struct gkyl_array* cell_avg_magB2, struct gkyl_array* temp_var)
{
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&conf_range, idx);
    // fetch the correct count in the matrix (since we solve Ncomp systems in each cell)
    long count = linc1*up->Ncomp;

    const double *em_d = (const double*) gkyl_array_cfetch(em, loc);
    int *cell_avg_magB2_d = (int*) gkyl_array_fetch(cell_avg_magB2, loc);

    up->em_calc_temp(em_d, (double*) gkyl_array_fetch(temp_var, loc));
    cell_avg_magB2_d[0] = up->em_set(count, As, xs, (const double*) gkyl_array_cfetch(temp_var, loc));
  }
}

__global__ static void
gkyl_dg_calc_em_vars_copy_cu_kernel(gkyl_dg_calc_em_vars* up, 
  struct gkyl_nmat *xs, struct gkyl_range conf_range,
  const struct gkyl_array* em, struct gkyl_array* cell_avg_magB2, struct gkyl_array* out)
{
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&conf_range, idx);
    // fetch the correct count in the matrix (since we solve Ncomp systems in each cell)
    long count = linc1*up->Ncomp;

    const double *em_d = (const double*) gkyl_array_cfetch(em, loc);
    int *cell_avg_magB2_d = (int*) gkyl_array_fetch(cell_avg_magB2, loc);
    double *out_d = (double*) gkyl_array_fetch(out, loc);

    up->em_copy(count, xs, em_d, cell_avg_magB2_d, out_d);
  }
}

void gkyl_dg_calc_em_vars_advance_cu(struct gkyl_dg_calc_em_vars *up, 
  const struct gkyl_array* em, struct gkyl_array* cell_avg_magB2, struct gkyl_array* out)
{
  gkyl_array_clear(up->temp_var, 0.0);
  struct gkyl_range conf_range = up->mem_range;
  
  gkyl_dg_calc_em_vars_set_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->As->on_dev, up->xs->on_dev, conf_range,
    em->on_dev, cell_avg_magB2->on_dev, up->temp_var->on_dev);

  if (up->poly_order > 1) {
    bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
    assert(status);
  }

  gkyl_dg_calc_em_vars_copy_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->xs->on_dev, conf_range,
    em->on_dev, cell_avg_magB2->on_dev, out->on_dev);
}

__global__ static void
gkyl_dg_calc_em_vars_surf_set_cu_kernel(gkyl_dg_calc_em_vars* up,
  struct gkyl_nmat *As, struct gkyl_nmat *xs, struct gkyl_range conf_range,
  const struct gkyl_array* em, struct gkyl_array* cell_avg_magB2_surf, struct gkyl_array* temp_var_surf) 
{
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&conf_range, idx);
    // fetch the correct count in the matrix (since we solve Ncomp systems in each cell)
    long count = linc1*up->Ncomp_surf;

    const double *em_d = (const double*) gkyl_array_cfetch(em, loc);
    int *cell_avg_magB2_surf_d = (int*) gkyl_array_fetch(cell_avg_magB2_surf, loc);

    up->em_calc_surf_temp(em_d, (double*) gkyl_array_fetch(temp_var_surf, loc));
    up->em_surf_set(count, As, xs, 
      (const double*) gkyl_array_cfetch(temp_var_surf, loc), cell_avg_magB2_surf_d);
  }
}

__global__ static void
gkyl_dg_calc_em_vars_copy_cu_kernel(gkyl_dg_calc_em_vars* up, 
  struct gkyl_nmat *xs, struct gkyl_range conf_range,
  const struct gkyl_array* em, struct gkyl_array* cell_avg_magB2_surf, struct gkyl_array* bvar_surf)
{
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&conf_range, idx);
    // fetch the correct count in the matrix (since we solve Ncomp_surf systems in each cell)
    long count = linc1*up->Ncomp_surf;

    const double *em_d = (const double*) gkyl_array_cfetch(em, loc);
    int *cell_avg_magB2_surf_d = (int*) gkyl_array_fetch(cell_avg_magB2_surf, loc);
    double *bvar_surf_d = (double*) gkyl_array_fetch(bvar_surf, loc);

    up->em_surf_copy(count, xs, em_d, cell_avg_magB2_surf_d, bvar_surf_d);
  }
}

void gkyl_dg_calc_em_vars_surf_advance_cu(struct gkyl_dg_calc_em_vars *up, 
  const struct gkyl_array* em, struct gkyl_array* cell_avg_magB2_surf, struct gkyl_array* bvar_surf)
{
  struct gkyl_range conf_range = up->mem_range;
  
  gkyl_dg_calc_em_vars_surf_set_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->As_surf->on_dev, up->xs_surf->on_dev, conf_range,
    em->on_dev, cell_avg_magB2_surf->on_dev, up->temp_var_surf->on_dev);

  if (up->poly_order > 1 && up->cdim > 1) {
    bool status = gkyl_nmat_linsolve_lu_pa(up->mem_surf, up->As_surf, up->xs_surf);
    assert(status);
  }

  gkyl_dg_calc_em_vars_surf_copy_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->xs_surf->on_dev, conf_range, 
    em->on_dev, cell_avg_magB2_surf->on_dev, bvar_surf->on_dev);
}

__global__ void
gkyl_dg_calc_em_vars_div_b_cu_kernel(struct gkyl_dg_calc_em_vars *up, struct gkyl_range conf_range, 
  const struct gkyl_array* bvar_surf, const struct gkyl_array* bvar, 
  struct gkyl_array* max_b, struct gkyl_array* div_b)
{
  int cdim = up->cdim;
  int idxl[GKYL_MAX_DIM], idxc[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_range, linc1, idxc);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long linc = gkyl_range_idx(&conf_range, idxc);

    const double *bvar_surf_c = (const double*) gkyl_array_cfetch(bvar_surf, linc);
    const double *bvar_d = (const double*) gkyl_array_cfetch(bvar, linc);

    double *max_b_d = (double*) gkyl_array_fetch(max_b, linc);
    double *div_b_d = (double*) gkyl_array_fetch(div_b, linc);

    for (int dir=0; dir<cdim; ++dir) {
      gkyl_copy_int_arr(cdim, idxc, idxl);
      gkyl_copy_int_arr(cdim, idxc, idxr);

      idxl[dir] = idxl[dir]-1; idxr[dir] = idxr[dir]+1;

      long linl = gkyl_range_idx(&conf_range, idxl); 
      long linr = gkyl_range_idx(&conf_range, idxr);

      const double *bvar_surf_l = (const double*) gkyl_array_cfetch(bvar_surf, linl);
      const double *bvar_surf_r = (const double*) gkyl_array_cfetch(bvar_surf, linr);
      
      up->em_div_b[dir](up->conf_grid.dx, 
        bvar_surf_l, bvar_surf_c, bvar_surf_r, 
        bvar_d, max_b_d, div_b_d);
    }
  }
}

// Host-side wrapper for div(b) and max(|b_i|) variable calculations 
void
gkyl_dg_calc_em_vars_div_b_cu(struct gkyl_dg_calc_em_vars *up, const struct gkyl_range *conf_range, 
  const struct gkyl_array* bvar_surf, const struct gkyl_array* bvar, 
  struct gkyl_array* max_b, struct gkyl_array* div_b)
{
  int nblocks = conf_range->nblocks;
  int nthreads = conf_range->nthreads;
  gkyl_dg_calc_em_vars_div_b_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, *conf_range, 
    bvar_surf->on_dev, bvar->on_dev, 
    max_b->on_dev, div_b->on_dev);
}

// CUDA kernel to set device pointers to em vars kernel functions
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
dg_calc_em_vars_set_cu_dev_ptrs(struct gkyl_dg_calc_em_vars *up, enum gkyl_basis_type b_type,
  int cdim,int poly_order, bool is_ExB)
{
  if (is_ExB) {
    up->em_calc_temp = choose_em_calc_num_ExB_kern(b_type, cdim, poly_order);
    up->em_set = choose_em_set_ExB_kern(b_type, cdim, poly_order);
    up->em_copy = choose_em_copy_ExB_kern(b_type, cdim, poly_order);
  }
  else {
    up->em_calc_temp = choose_em_calc_BB_kern(b_type, cdim, poly_order);
    up->em_calc_surf_temp = choose_em_calc_surf_BB_kern(b_type, cdim, poly_order);
    up->em_set = choose_em_set_bvar_kern(b_type, cdim, poly_order);
    up->em_surf_set = choose_em_surf_set_bvar_kern(b_type, cdim, poly_order);
    up->em_copy = choose_em_copy_bvar_kern(b_type, cdim, poly_order);    
    up->em_surf_copy = choose_em_surf_copy_bvar_kern(b_type, cdim, poly_order);    
    // Fetch the kernels in each direction
    for (int d=0; d<cdim; ++d) 
      up->em_div_b[d] = choose_em_div_b_kern(d, b_type, cdim, poly_order);       
  }
}

gkyl_dg_calc_em_vars*
gkyl_dg_calc_em_vars_cu_dev_new(const struct gkyl_rect_grid *conf_grid, 
  const struct gkyl_basis* cbasis, const struct gkyl_range *mem_range, 
  bool is_ExB)
{
  struct gkyl_dg_calc_em_vars *up = (struct gkyl_dg_calc_em_vars*) gkyl_malloc(sizeof(gkyl_dg_calc_em_vars));

  up->conf_grid = *conf_grid;
  int nc = cbasis->num_basis;
  int cdim = cbasis->ndim;
  int poly_order = cbasis->poly_order;
  int nc_surf = cbasis->num_basis/(poly_order+1); // *only valid for tensor bases for cdim > 1*
  enum gkyl_basis_type b_type = cbasis->b_type;
  up->cdim = cdim;
  up->poly_order = poly_order;
  up->mem_range = *mem_range;

  if (is_ExB) {
    up->Ncomp = 3;
    up->Ncomp_surf = 2*cdim;
  }
  else {
    up->Ncomp = 6;
    up->Ncomp_surf = 2*cdim*4;
  }

  // There are Ncomp more linear systems to be solved 
  // 6 components of bb and 3 components of E x B
  up->As = gkyl_nmat_cu_dev_new(up->Ncomp*mem_range->volume, nc, nc);
  up->xs = gkyl_nmat_cu_dev_new(up->Ncomp*mem_range->volume, nc, 1);
  up->mem = gkyl_nmat_linsolve_lu_cu_dev_new(up->As->num, up->As->nr);
  // 6 component temporary variable for either storing B_i B_j (for computing bb) 
  // or (E x B)_i and B_i^2 (for computing E x B/|B|^2)
  up->temp_var = gkyl_array_cu_dev_new(GKYL_DOUBLE, 6*nc, mem_range->volume);

  // There are Ncomp_surf*range->volume linear systems to be solved 
  // Each linear system is nc_surf x nc_surf (only solved over the surface basis and only when poly_order and cdim > 1)
  // 2*cdim*4: bx, bxbx, bxby, bxbz (xl and xr), by, byby, bxby, bybz (yl and yr), bz, bzbz, bxbz, bybz (zl and zr)  
  up->As_surf = gkyl_nmat_cu_dev_new(up->Ncomp_surf*mem_range->volume, nc_surf, nc_surf);
  up->xs_surf = gkyl_nmat_cu_dev_new(up->Ncomp_surf*mem_range->volume, nc_surf, 1);
  up->mem_surf = gkyl_nmat_linsolve_lu_cu_dev_new(up->As_surf->num, up->As_surf->nr);
  // 5*cdim component temporary variable for storing Bx^2, By^2, Bz^2 and two components of B_i B_j at the surface
  // Temporary variables are computed at the left and right in each dimension (xl, xr, yl, yr, zl, & zr)
  up->temp_var_surf = gkyl_array_new(GKYL_DOUBLE, 10*cdim*nc_surf, mem_range->volume); 

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_dg_calc_em_vars *up_cu = (struct gkyl_dg_calc_em_vars*) gkyl_cu_malloc(sizeof(gkyl_dg_calc_em_vars));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_dg_calc_em_vars), GKYL_CU_MEMCPY_H2D);

  dg_calc_em_vars_set_cu_dev_ptrs<<<1,1>>>(up_cu, b_type, cdim, poly_order, is_ExB);

  // set parent on_dev pointer
  up->on_dev = up_cu;
  
  return up;
}
