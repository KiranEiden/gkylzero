/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_fem_poisson.h>
#include <gkyl_fem_poisson_priv.h>
}

// CUDA kernel to set device pointers to l2g kernel function.
// Doing function pointer stuff in here avoids troublesome
// cudaMemcpyFromSymbol.
__global__ static void
fem_poisson_set_cu_l2gker_ptrs(struct gkyl_fem_poisson_kernels* kers, enum gkyl_basis_type b_type,
  int dim, int poly_order, const int *bckey)
{

  // Set l2g kernels.
  const local2global_kern_bcx_list_1x *local2global_1x_kernels;
  const local2global_kern_bcx_list_2x *local2global_2x_kernels;

  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      local2global_1x_kernels = ser_loc2glob_list_1x;
      local2global_2x_kernels = ser_loc2glob_list_2x;
      break;
    default:
      assert(false);
      break;
  }

  for (int k=0; k<(int)(pow(2,dim)+0.5); k++) {
    if (dim == 1) {
      kers->l2g[k] = CK1(local2global_1x_kernels, poly_order, k, bckey[0]);
    } else if ( dim == 2) {
      kers->l2g[k] = CK2(local2global_2x_kernels, poly_order, k, bckey[0], bckey[1]);
//    } else if (dim == 3) {
//      kers->l2g[k] = CK3(ser_loc2glob_list_3x, poly_order, k, bckey[0], bckey[1], bckey[2]);
    }
  }

}

__global__ static void
fem_poisson_set_cu_ker_ptrs(struct gkyl_fem_poisson_kernels* kers, enum gkyl_basis_type b_type,
  int dim, int poly_order, const int *bckey)
{

//  // Set LHS stencil kernels.
//  const lhsstencil_kern_bcx_list_1x *lhsstencil_1x_kernels;
//  const lhsstencil_kern_bcx_list_2x *lhsstencil_2x_kernels;
//
//  switch (b_type) {
//    case GKYL_BASIS_MODAL_SERENDIPITY:
//        lhsstencil_1x_kernels = ser_lhsstencil_list_1x;
//        lhsstencil_2x_kernels = ser_lhsstencil_list_2x;
//      break;
////    case GKYL_BASIS_MODAL_TENSOR:
////      break;
//    default:
//      assert(false);
//  }
//
//  for (int k=0; k<(int)(pow(3,dim)+0.5); k++) {
//    if (dim == 1) {
//      kers->lhsker[k] = CK1(lhsstencil_1x_kernels, poly_order, k, bckey[0]);
//    } else if (dim == 2) {
//      kers->lhsker[k] = CK2(lhsstencil_2x_kernels, poly_order, k, bckey[0], bckey[1]);
////  } else if (dim == 3) {
////    kers->lhsker[k] = CK3(lhsstencil_3x_kernels, poly_order, k, bckey[0], bckey[1], bckey[2]);
//    }
//  }

  // Set RHS stencil kernels.
  const srcstencil_kern_bcx_list_1x *srcstencil_1x_kernels;
  const srcstencil_kern_bcx_list_2x *srcstencil_2x_kernels;

  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
        srcstencil_1x_kernels = ser_srcstencil_list_1x;
        srcstencil_2x_kernels = ser_srcstencil_list_2x;
      break;
//    case GKYL_BASIS_MODAL_TENSOR:
//      break;
    default:
      assert(false);
  }

  for (int k=0; k<(int)(pow(3,dim)+0.5); k++) {
    if (dim == 1) {
      kers->srcker[k] = CK1(srcstencil_1x_kernels, poly_order, k, bckey[0]);
    } else if (dim == 2) {
      kers->srcker[k] = CK2(srcstencil_2x_kernels, poly_order, k, bckey[0], bckey[1]);
//  } else if (dim == 3) {
//    kers->srcker[k] = CK3(srcstencil_3x_kernels, poly_order, k, bckey[0], bckey[1], bckey[2]);
    }
  }

  // Set the get solution stencil kernel.
  const solstencil_kern_list *solstencil_kernels;

  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
        solstencil_kernels = ser_solstencil_list;
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      break;
    default:
      assert(false);
  }

  kers->solker = solstencil_kernels[dim].kernels[poly_order];

}

void
choose_kernels_cu(const struct gkyl_basis* basis, const struct gkyl_poisson_bc bcs, const bool *isdirperiodic, struct gkyl_fem_poisson_kernels *kers)
{

  int dim = basis->ndim;
  int poly_order = basis->poly_order;

  int bckey[POISSON_MAX_DIM] = {-1};
  for (int d=0; d<basis->ndim; d++) bckey[d] = isdirperiodic[d] ? 0 : 1;
  int *bckey_d = (int *) gkyl_cu_malloc(sizeof(int[POISSON_MAX_DIM]));
  gkyl_cu_memcpy(bckey_d, bckey, sizeof(int[POISSON_MAX_DIM]), GKYL_CU_MEMCPY_H2D);

  fem_poisson_set_cu_l2gker_ptrs<<<1,1>>>(kers, basis->b_type, dim, poly_order, bckey_d);
  
  for (int d=0; d<basis->ndim; d++) {
    if (bcs.lo_type[d]==GKYL_POISSON_PERIODIC && bcs.up_type[d]==GKYL_POISSON_PERIODIC) { bckey[d] = 0; }
    else if (bcs.lo_type[d]==GKYL_POISSON_DIRICHLET && bcs.up_type[d]==GKYL_POISSON_DIRICHLET) { bckey[d] = 1; }
    else if (bcs.lo_type[d]==GKYL_POISSON_DIRICHLET && bcs.up_type[d]==GKYL_POISSON_NEUMANN) { bckey[d] = 2; }
    else if (bcs.lo_type[d]==GKYL_POISSON_NEUMANN && bcs.up_type[d]==GKYL_POISSON_DIRICHLET) { bckey[d] = 3; }
    else if (bcs.lo_type[d]==GKYL_POISSON_DIRICHLET && bcs.up_type[d]==GKYL_POISSON_ROBIN) { bckey[d] = 4; }
    else if (bcs.lo_type[d]==GKYL_POISSON_ROBIN && bcs.up_type[d]==GKYL_POISSON_DIRICHLET) { bckey[d] = 5; }
    else { assert(false); }
  };
  gkyl_cu_memcpy(bckey_d, bckey, sizeof(int[POISSON_MAX_DIM]), GKYL_CU_MEMCPY_H2D);

  fem_poisson_set_cu_ker_ptrs<<<1,1>>>(kers, basis->b_type, dim, poly_order, bckey_d);

  gkyl_cu_free(bckey_d);
}

//void 
//gkyl_fem_poisson_set_rhs_cu(gkyl_fem_poisson *up, struct gkyl_array *rhsin)
//{
//  gkyl_fem_poisson_set_rhs_kernel<<<dG, dB>>>(up->rhs, rhsin->on_dev, &up->solve_range, up->globalidx_cu, up->bcvals); 
//}	
//
//void
//gkyl_fem_poisson_solve_cu(gkyl_fem_poisson *up, struct gkyl_array *phiin)
//{
//  // do linear solve with cusolver
//  gkyl_cusolver_solve(up->prob_cu);
//
//  gkyl_fem_poisson_get_sol_kernel<<<dG, dB>>>(up->rhs, rhsin->on_dev, &up->solve_range, up->globalidx_cu, up->bcvals); 
//}

//__global__ void
//gkyl_fem_poisson_set_rhs_kernel(double *rhs_global, struct gkyl_array *rhs_local, struct gkyl_range range, long *globalidx, double *bcvals, struct gkyl_fem_poisson_kernels *kers)
//{
//  int idx[GKYL_MAX_DIM];
//  int idx0[GKYL_MAX_DIM];
//  int num_cells[POISSON_MAX_DIM];
//  for (int d=0; d<POISSON_MAX_DIM; d++) num_cells[d] = range.upper[d]-range.lower[d]+1;
//
//  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
//      linc1 < range.volume;
//      linc1 += gridDim.x*blockDim.x)
//  {
//    // inverse index from linc1 to idx
//    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
//    // since update_range is a subrange
//    gkyl_sub_range_inv_idx(&range, linc1, idx);
//    
//    // convert back to a linear index on the super-range (with ghost cells)
//    // linc will have jumps in it to jump over ghost cells
//    long start = gkyl_range_idx(&range, idx);
//    
//    const double *local_d = (const double*) gkyl_array_cfetch(rhs_local, start);
//
//    int keri = idx_to_inup_ker(range->ndim, num_cells, idx);
//
//    for (size_t d=0; d<range->ndim; d++) idx0[d] = idx[d]-1;
//
//    kers->l2g[keri](num_cells, idx0, globalidx);
//
//    // Apply the RHS source stencil. It's mostly the mass matrix times a
//    // modal-to-nodal operator times the source, modified by BCs in skin cells.
//    keri = idx_to_inloup_ker(range->ndim, num_cells, idx);
//    kers->srcker[keri](local_d, bcvals, globalidx, rhs_global);
//  }
//}
