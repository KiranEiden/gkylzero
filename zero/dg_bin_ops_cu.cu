/* -*- c++ -*- */

#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_bin_ops_priv.h>
#include <gkyl_mat.h>
#include <gkyl_util.h>
}

// start ID for use in various loops
#define START_ID (threadIdx.x + blockIdx.x*blockDim.x)

void
gkyl_get_array_range_kernel_launch_dims(dim3* dimGrid, dim3* dimBlock, gkyl_range range, int ncomp)
{
  int volume = range.volume;
  int ndim = range.ndim;
  // ac1 = size of last dimension of range (fastest moving dimension)
  int ac1 = range.iac[ndim-1] > 0 ? range.iac[ndim-1] : 1;
  dimBlock->x = min(ncomp*ac1, GKYL_DEFAULT_NUM_THREADS);
  dimGrid->x = gkyl_int_div_up(ncomp*ac1, dimBlock->x);

  dimBlock->y = gkyl_int_div_up(GKYL_DEFAULT_NUM_THREADS, ncomp*ac1);
  dimGrid->y = gkyl_int_div_up(volume, ac1*dimBlock->y);
}

__global__ void
gkyl_dg_mul_op_cu_kernel(struct gkyl_basis basis,
  int c_oop, struct gkyl_array* out,
  int c_lop, const struct gkyl_array* lop,
  int c_rop, const struct gkyl_array* rop)
{
  int num_basis = basis.num_basis;
  int ndim = basis.ndim;
  int poly_order = basis.poly_order;
  mul_op_t mul_op = choose_ser_mul_kern(ndim, poly_order);

  for (unsigned long linc = START_ID; linc < NSIZE(out); linc += blockDim.x*gridDim.x) {
    
    const double *lop_d = (const double*) gkyl_array_cfetch(lop, linc);
    const double *rop_d = (const double*) gkyl_array_cfetch(rop, linc);
    double *out_d = (double*) gkyl_array_fetch(out, linc);

    mul_op(lop_d+c_lop*num_basis, rop_d+c_rop*num_basis, out_d+c_oop*num_basis);
  }  
}

// Host-side wrapper for dg multiplication operation
void
gkyl_dg_mul_op_cu(struct gkyl_basis basis,
  int c_oop, struct gkyl_array* out,
  int c_lop, const struct gkyl_array* lop,
  int c_rop, const struct gkyl_array* rop)
{
  gkyl_dg_mul_op_cu_kernel<<<out->nblocks, out->nthreads>>>(basis, c_oop, out->on_dev, c_lop, lop->on_dev, c_rop, rop->on_dev);
}

__global__ void
gkyl_dg_mul_op_range_cu_kernel(struct gkyl_basis basis,
  int c_oop, struct gkyl_array* out,
  int c_lop, const struct gkyl_array* lop,
  int c_rop, const struct gkyl_array* rop, struct gkyl_range range)
{
  int num_basis = basis.num_basis;
  int ndim = basis.ndim;
  int poly_order = basis.poly_order;
  mul_op_t mul_op = choose_ser_mul_kern(ndim, poly_order);

  long n = NCOM(out);
  int idx[GKYL_MAX_DIM];
  int ndim = range.ndim;
  // ac1 = size of last dimension of range (fastest moving dimension)
  long ac1 = range.iac[ndim-1] > 0 ? range.iac[ndim-1] : 1;

  // 2D thread grid
  // linc1 = c + n*idx1 (contiguous data, including component index c, with idx1 = 0,.., ac1-1)
  long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
  // linc2 = idx2 + ac2*idx3 + ...
  for (unsigned long linc2 = threadIdx.y + blockIdx.y*blockDim.y;
      linc2 < range.volume/ac1;
      linc2 += gridDim.y*blockDim.y)
  {
    // full linear cell index (not including components) is 
    // idx1 + ac1*idx2 + ac1*ac2*idx3 + ... = idx1 + ac1*linc2.
    // we want to find the start linear index of each contiguous data block, 
    // which corresponds to idx1 = 0. 
    // so linear index of start of contiguous block is ac1*linc2.
    gkyl_sub_range_inv_idx(&range, ac1*linc2, idx);
    long start = gkyl_range_idx(&range, idx);

    const double *lop_d = (const double*) gkyl_array_cfetch(lop, start);
    const double *rop_d = (const double*) gkyl_array_cfetch(rop, start);
    double *out_d = (double*) gkyl_array_fetch(out, start);

    mul_op(lop_d+c_lop*num_basis, rop_d+c_rop*num_basis, out_d+c_oop*num_basis);  
  }
}

// Host-side wrapper for range-based dg multiplication operation
void
gkyl_dg_mul_op_range_cu(struct gkyl_basis basis,
  int c_oop, struct gkyl_array* out,
  int c_lop, const struct gkyl_array* lop,
  int c_rop, const struct gkyl_array* rop, struct gkyl_range range)
{
  dim3 dimGrid, dimBlock;
  gkyl_get_array_range_kernel_launch_dims(&dimGrid, &dimBlock, range, out->ncomp);

  gkyl_array_clear_range_cu_kernel<<<dimGrid, dimBlock>>>(basis, c_oop, out->on_dev, c_lop, lop->on_dev, c_rop, rop->on_dev, range);
}

__global__ void
gkyl_dg_div_set_op_cu_kernel(struct gkyl_nmat *As, struct gkyl_nmat *xs,
  struct gkyl_basis basis, struct gkyl_array* out,
  int c_lop, const struct gkyl_array* lop,
  int c_rop, const struct gkyl_array* rop)
{
  int num_basis = basis.num_basis;
  int ndim = basis.ndim;
  int poly_order = basis.poly_order;
  div_set_op_t div_set_op = choose_ser_div_set_kern(ndim, poly_order);

  for (unsigned long linc = START_ID; linc < NSIZE(out); linc += blockDim.x*gridDim.x) {
    
    const double *lop_d = (const double*) gkyl_array_cfetch(lop, linc);
    const double *rop_d = (const double*) gkyl_array_cfetch(rop, linc);

    struct gkyl_mat A = gkyl_nmat_get(As, linc);
    struct gkyl_mat x = gkyl_nmat_get(xs, linc);
    gkyl_mat_clear(&A, 0.0); gkyl_mat_clear(&x, 0.0);  
    div_set_op(&A, &x, lop_d+c_lop*num_basis, rop_d+c_rop*num_basis);
  }  
}

__global__ void
gkyl_dg_div_copy_sol_op_cu_kernel(struct gkyl_nmat *xs,
  struct gkyl_basis basis,
  int c_oop, struct gkyl_array* out)
{
  int num_basis = basis.num_basis;
  for (unsigned long linc = START_ID; linc < NSIZE(out); linc += blockDim.x*gridDim.x) {
    double *out_d = (double*) gkyl_array_fetch(out, linc);
    struct gkyl_mat x = gkyl_nmat_get(xs, linc);
    binop_div_copy_sol(&x, out_d+c_oop*num_basis);
  }  
}

// Host-side wrapper for dg division operation
void
gkyl_dg_div_op_cu(struct gkyl_basis basis,
  int c_oop, struct gkyl_array* out,
  int c_lop, const struct gkyl_array* lop,
  int c_rop, const struct gkyl_array* rop)
{
  // allocate memory for use in kernels
  struct gkyl_nmat *A_d = gkyl_nmat_cu_dev_new(out->size, num_basis, num_basis);
  struct gkyl_nmat *x_d = gkyl_nmat_cu_dev_new(out->size, num_basis, 1);

  gkyl_dg_div_set_op_cu_kernel<<<out->nblocks, out->nthreads>>>(A_d, x_d, basis, out->on_dev, c_lop, lop->on_dev, c_rop, rop->on_dev);

  bool status = gkyl_nmat_linsolve_lu(A_d, x_d);

  gkyl_dg_div_copy_sol_op_cu_kernel<<<out->nblocks, out->nthreads>>>(x_d, basis, c_oop, out->on_dev)

  gkyl_nmat_release(As);
  gkyl_nmat_release(xs);  
}

__global__ void
gkyl_dg_div_set_op_range_cu_kernel(struct gkyl_nmat *As, struct gkyl_nmat *xs,
  struct gkyl_basis basis, struct gkyl_array* out,
  int c_lop, const struct gkyl_array* lop,
  int c_rop, const struct gkyl_array* rop, struct gkyl_range range)
{
  int num_basis = basis.num_basis;
  int ndim = basis.ndim;
  int poly_order = basis.poly_order;
  div_set_op_t div_set_op = choose_ser_div_set_kern(ndim, poly_order);

  long n = NCOM(out);
  int idx[GKYL_MAX_DIM];
  int ndim = range.ndim;
  // ac1 = size of last dimension of range (fastest moving dimension)
  long ac1 = range.iac[ndim-1] > 0 ? range.iac[ndim-1] : 1;

  // 2D thread grid
  // linc1 = c + n*idx1 (contiguous data, including component index c, with idx1 = 0,.., ac1-1)
  long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
  // linc2 = idx2 + ac2*idx3 + ...
  for (unsigned long linc2 = threadIdx.y + blockIdx.y*blockDim.y;
      linc2 < range.volume/ac1;
      linc2 += gridDim.y*blockDim.y)
  {
    // full linear cell index (not including components) is 
    // idx1 + ac1*idx2 + ac1*ac2*idx3 + ... = idx1 + ac1*linc2.
    // we want to find the start linear index of each contiguous data block, 
    // which corresponds to idx1 = 0. 
    // so linear index of start of contiguous block is ac1*linc2.
    gkyl_sub_range_inv_idx(&range, ac1*linc2, idx);
    long start = gkyl_range_idx(&range, idx);
    
    const double *lop_d = (const double*) gkyl_array_cfetch(lop, start);
    const double *rop_d = (const double*) gkyl_array_cfetch(rop, start);

    struct gkyl_mat A = gkyl_nmat_get(As, start);
    struct gkyl_mat x = gkyl_nmat_get(xs, start);
    gkyl_mat_clear(&A, 0.0); gkyl_mat_clear(&x, 0.0);  
    div_set_op(&A, &x, lop_d+c_lop*num_basis, rop_d+c_rop*num_basis);
  }  
}

__global__ void
gkyl_dg_div_copy_sol_op_range_cu_kernel(struct gkyl_nmat *xs,
  struct gkyl_basis basis,
  int c_oop, struct gkyl_array* out, struct gkyl_range range)
{
  int num_basis = basis.num_basis;
  long n = NCOM(out);
  int idx[GKYL_MAX_DIM];
  int ndim = range.ndim;
  // ac1 = size of last dimension of range (fastest moving dimension)
  long ac1 = range.iac[ndim-1] > 0 ? range.iac[ndim-1] : 1;

  // 2D thread grid
  // linc1 = c + n*idx1 (contiguous data, including component index c, with idx1 = 0,.., ac1-1)
  long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
  // linc2 = idx2 + ac2*idx3 + ...
  for (unsigned long linc2 = threadIdx.y + blockIdx.y*blockDim.y;
      linc2 < range.volume/ac1;
      linc2 += gridDim.y*blockDim.y)
  {
    // full linear cell index (not including components) is 
    // idx1 + ac1*idx2 + ac1*ac2*idx3 + ... = idx1 + ac1*linc2.
    // we want to find the start linear index of each contiguous data block, 
    // which corresponds to idx1 = 0. 
    // so linear index of start of contiguous block is ac1*linc2.
    gkyl_sub_range_inv_idx(&range, ac1*linc2, idx);
    long start = gkyl_range_idx(&range, idx);

    double *out_d = (double*) gkyl_array_fetch(out, start);
    struct gkyl_mat x = gkyl_nmat_get(xs, start);
    binop_div_copy_sol(&x, out_d+c_oop*num_basis);
  }  
}

// Host-side wrapper for range-based dg division operation
void
gkyl_dg_div_op_range_cu(struct gkyl_basis basis,
  int c_oop, struct gkyl_array* out,
  int c_lop, const struct gkyl_array* lop,
  int c_rop, const struct gkyl_array* rop, struct gkyl_range range)
{
  dim3 dimGrid, dimBlock;
  gkyl_get_array_range_kernel_launch_dims(&dimGrid, &dimBlock, range, out->ncomp);

  // allocate memory for use in kernels
  struct gkyl_nmat *A_d = gkyl_nmat_cu_dev_new(out->size, num_basis, num_basis);
  struct gkyl_nmat *x_d = gkyl_nmat_cu_dev_new(out->size, num_basis, 1);

  gkyl_dg_div_set_op_cu_kernel<<<dimGrid, dimBlock>>>(A_d, x_d, basis, out->on_dev, c_lop, lop->on_dev, c_rop, rop->on_dev, range);

  bool status = gkyl_nmat_linsolve_lu(A_d, x_d);

  gkyl_dg_div_copy_sol_op_cu_kernel<<<dimGrid, dimBlock>>>(x_d, basis, c_oop, out->on_dev, range)

  gkyl_nmat_release(As);
  gkyl_nmat_release(xs);  
}