// Private header for fem_poisson updater.
#pragma once
#include <gkyl_fem_poisson_kernels.h>
#include <gkyl_basis.h>
#include <gkyl_superlu_ops.h>
#include <gkyl_cusolver_ops.h>

#ifndef GKYL_IPOW
# define GKYL_IPOW(a,e) (int)(pow(a,e)+0.5)
#endif

// Function pointer type for local-to-global mapping.
typedef void (*local2global_t)(const int *numCells, const int *idx,
  long *globalIdxs);

// For use in kernel tables.
typedef struct { local2global_t kernels[2]; } local2global_kern_loc_list_1x;
typedef struct { local2global_kern_loc_list_1x list[3]; } local2global_kern_bcx_list_1x;

typedef struct { local2global_t kernels[4]; } local2global_kern_loc_list_2x;
typedef struct { local2global_kern_loc_list_2x list[3]; } local2global_kern_bcy_list_2x;
typedef struct { local2global_kern_bcy_list_2x list[2]; } local2global_kern_bcx_list_2x;

// Serendipity local-to-global kernels.
GKYL_CU_D
static const local2global_kern_bcx_list_1x ser_loc2glob_list_1x[] = {
  // periodicx
  { .list = {{NULL, NULL},
             {fem_poisson_local_to_global_1x_ser_p1_inx_periodicx, fem_poisson_local_to_global_1x_ser_p1_upx_periodicx},
             {fem_poisson_local_to_global_1x_ser_p2_inx_periodicx, fem_poisson_local_to_global_1x_ser_p2_upx_periodicx}}, },
  // nonperiodicx
  { .list = {{NULL, NULL},
            {fem_poisson_local_to_global_1x_ser_p1_inx_nonperiodicx, fem_poisson_local_to_global_1x_ser_p1_upx_nonperiodicx},
            {fem_poisson_local_to_global_1x_ser_p2_inx_nonperiodicx, fem_poisson_local_to_global_1x_ser_p2_upx_nonperiodicx}}, }
};

GKYL_CU_D
static const local2global_kern_bcx_list_2x ser_loc2glob_list_2x[] = {
  // periodicx
  { .list = {
    // periodicy
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_local_to_global_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_local_to_global_2x_ser_p1_upx_periodicx_iny_periodicy, fem_poisson_local_to_global_2x_ser_p1_inx_periodicx_upy_periodicy, fem_poisson_local_to_global_2x_ser_p1_upx_periodicx_upy_periodicy},
               {fem_poisson_local_to_global_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_local_to_global_2x_ser_p2_upx_periodicx_iny_periodicy, fem_poisson_local_to_global_2x_ser_p2_inx_periodicx_upy_periodicy, fem_poisson_local_to_global_2x_ser_p2_upx_periodicx_upy_periodicy}},
    },
    // nonperiodicy
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_local_to_global_2x_ser_p1_inx_periodicx_iny_nonperiodicy, fem_poisson_local_to_global_2x_ser_p1_upx_periodicx_iny_nonperiodicy, fem_poisson_local_to_global_2x_ser_p1_inx_periodicx_upy_nonperiodicy, fem_poisson_local_to_global_2x_ser_p1_upx_periodicx_upy_nonperiodicy},
               {fem_poisson_local_to_global_2x_ser_p2_inx_periodicx_iny_nonperiodicy, fem_poisson_local_to_global_2x_ser_p2_upx_periodicx_iny_nonperiodicy, fem_poisson_local_to_global_2x_ser_p2_inx_periodicx_upy_nonperiodicy, fem_poisson_local_to_global_2x_ser_p2_upx_periodicx_upy_nonperiodicy},}
    }}
  },
  // nonperiodicx
  { .list = {
    // periodicy
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_local_to_global_2x_ser_p1_inx_nonperiodicx_iny_periodicy, fem_poisson_local_to_global_2x_ser_p1_upx_nonperiodicx_iny_periodicy, fem_poisson_local_to_global_2x_ser_p1_inx_nonperiodicx_upy_periodicy, fem_poisson_local_to_global_2x_ser_p1_upx_nonperiodicx_upy_periodicy},
               {fem_poisson_local_to_global_2x_ser_p2_inx_nonperiodicx_iny_periodicy, fem_poisson_local_to_global_2x_ser_p2_upx_nonperiodicx_iny_periodicy, fem_poisson_local_to_global_2x_ser_p2_inx_nonperiodicx_upy_periodicy, fem_poisson_local_to_global_2x_ser_p2_upx_nonperiodicx_upy_periodicy}},
    },
    // nonperiodicy
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_local_to_global_2x_ser_p1_inx_nonperiodicx_iny_nonperiodicy, fem_poisson_local_to_global_2x_ser_p1_upx_nonperiodicx_iny_nonperiodicy, fem_poisson_local_to_global_2x_ser_p1_inx_nonperiodicx_upy_nonperiodicy, fem_poisson_local_to_global_2x_ser_p1_upx_nonperiodicx_upy_nonperiodicy},
               {fem_poisson_local_to_global_2x_ser_p2_inx_nonperiodicx_iny_nonperiodicy, fem_poisson_local_to_global_2x_ser_p2_upx_nonperiodicx_iny_nonperiodicy, fem_poisson_local_to_global_2x_ser_p2_inx_nonperiodicx_upy_nonperiodicy, fem_poisson_local_to_global_2x_ser_p2_upx_nonperiodicx_upy_nonperiodicy},}
    }}
  }
};

// Function pointer type for lhs kernels.
typedef void (*lhsstencil_t)(const double *epsilon, const double *kSq, const double *dx, const double *bcVals,
  const long *globalIdxs, gkyl_mat_triples *tri);

// For use in kernel tables.
typedef struct { lhsstencil_t kernels[3]; } lhsstencil_kern_loc_list_1x;
typedef struct { lhsstencil_kern_loc_list_1x list[3]; } lhsstencil_kern_bcx_list_1x;

typedef struct { lhsstencil_t kernels[9]; } lhsstencil_kern_loc_list_2x;
typedef struct { lhsstencil_kern_loc_list_2x list[3]; } lhsstencil_kern_bcy_list_2x;
typedef struct { lhsstencil_kern_bcy_list_2x list[9]; } lhsstencil_kern_bcx_list_2x;

// Serendipity lhs kernels.
static const lhsstencil_kern_bcx_list_1x ser_lhsstencil_consteps_list_1x[] = {
  // periodicx
  { .list = {{NULL, NULL},
             {fem_poisson_lhs_stencil_consteps_1x_ser_p1_inx_periodicx, fem_poisson_lhs_stencil_consteps_1x_ser_p1_lox_periodicx, fem_poisson_lhs_stencil_consteps_1x_ser_p1_upx_periodicx},
             {fem_poisson_lhs_stencil_consteps_1x_ser_p2_inx_periodicx, fem_poisson_lhs_stencil_consteps_1x_ser_p2_lox_periodicx, fem_poisson_lhs_stencil_consteps_1x_ser_p2_upx_periodicx}}, },
  // dirichletx-dirichletx
  { .list = {{NULL, NULL},
             {fem_poisson_lhs_stencil_consteps_1x_ser_p1_inx_periodicx, fem_poisson_lhs_stencil_consteps_1x_ser_p1_lox_dirichletx, fem_poisson_lhs_stencil_consteps_1x_ser_p1_upx_dirichletx},
             {fem_poisson_lhs_stencil_consteps_1x_ser_p2_inx_periodicx, fem_poisson_lhs_stencil_consteps_1x_ser_p2_lox_dirichletx, fem_poisson_lhs_stencil_consteps_1x_ser_p2_upx_dirichletx}}, },
  // dirichletx-neumannx
  { .list = {{NULL, NULL},
             {fem_poisson_lhs_stencil_consteps_1x_ser_p1_inx_periodicx, fem_poisson_lhs_stencil_consteps_1x_ser_p1_lox_dirichletx, fem_poisson_lhs_stencil_consteps_1x_ser_p1_upx_neumannx},
             {fem_poisson_lhs_stencil_consteps_1x_ser_p2_inx_periodicx, fem_poisson_lhs_stencil_consteps_1x_ser_p2_lox_dirichletx, fem_poisson_lhs_stencil_consteps_1x_ser_p2_upx_neumannx}}, },
  // neumannx-dirichletx
  { .list = {{NULL, NULL},
             {fem_poisson_lhs_stencil_consteps_1x_ser_p1_inx_periodicx, fem_poisson_lhs_stencil_consteps_1x_ser_p1_lox_neumannx, fem_poisson_lhs_stencil_consteps_1x_ser_p1_upx_dirichletx},
             {fem_poisson_lhs_stencil_consteps_1x_ser_p2_inx_periodicx, fem_poisson_lhs_stencil_consteps_1x_ser_p2_lox_neumannx, fem_poisson_lhs_stencil_consteps_1x_ser_p2_upx_dirichletx}}, },
  // dirichletx-robinx
  { .list = {{NULL, NULL},
             {fem_poisson_lhs_stencil_consteps_1x_ser_p1_inx_periodicx, fem_poisson_lhs_stencil_consteps_1x_ser_p1_lox_dirichletx, fem_poisson_lhs_stencil_consteps_1x_ser_p1_upx_robinx},
             {fem_poisson_lhs_stencil_consteps_1x_ser_p2_inx_periodicx, fem_poisson_lhs_stencil_consteps_1x_ser_p2_lox_dirichletx, fem_poisson_lhs_stencil_consteps_1x_ser_p2_upx_robinx}}, },
  // robinx-dirichletx
  { .list = {{NULL, NULL},
             {fem_poisson_lhs_stencil_consteps_1x_ser_p1_inx_periodicx, fem_poisson_lhs_stencil_consteps_1x_ser_p1_lox_robinx, fem_poisson_lhs_stencil_consteps_1x_ser_p1_upx_dirichletx},
             {fem_poisson_lhs_stencil_consteps_1x_ser_p2_inx_periodicx, fem_poisson_lhs_stencil_consteps_1x_ser_p2_lox_robinx, fem_poisson_lhs_stencil_consteps_1x_ser_p2_upx_dirichletx}}, },
};

static const lhsstencil_kern_bcx_list_2x ser_lhsstencil_consteps_list_2x[] = {
  // periodicx
  { .list = {
    // periodicy
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_loy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_upy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_periodicx_loy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_periodicx_upy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_periodicx_loy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_periodicx_upy_periodicy},
               {fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_loy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_upy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_periodicx_loy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_periodicx_upy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_periodicx_loy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_periodicx_upy_periodicy},},
    },
    // dirichlety-dirichlety
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_periodicx_upy_dirichlety},
               {fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_periodicx_upy_dirichlety},},
    },
    // dirichlety-neumanny
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_upy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_periodicx_upy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_periodicx_upy_neumanny},
               {fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_upy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_periodicx_upy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_periodicx_upy_neumanny},},
    },
    // neumanny-dirichlety
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_loy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_periodicx_loy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_periodicx_loy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_periodicx_upy_dirichlety},
               {fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_loy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_periodicx_loy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_periodicx_loy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_periodicx_upy_dirichlety},}
    },
    }
  },
  // dirichletx-dirichletx
  { .list = {
    // periodicy
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_loy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_upy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_dirichletx_loy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_dirichletx_upy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_dirichletx_loy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_dirichletx_upy_periodicy},
               {fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_loy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_upy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_dirichletx_loy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_dirichletx_upy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_dirichletx_loy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_dirichletx_upy_periodicy},},
    },
    // dirichlety-dirichlety
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_dirichletx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_dirichletx_upy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_dirichletx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_dirichletx_upy_dirichlety},
               {fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_dirichletx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_dirichletx_upy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_dirichletx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_dirichletx_upy_dirichlety},},
    },
    // dirichlety-neumanny
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_upy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_dirichletx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_dirichletx_upy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_dirichletx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_dirichletx_upy_neumanny},
               {fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_upy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_dirichletx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_dirichletx_upy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_dirichletx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_dirichletx_upy_neumanny},},
    },
    // neumanny-dirichlety
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_loy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_dirichletx_loy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_dirichletx_upy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_dirichletx_loy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_dirichletx_upy_dirichlety},
               {fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_loy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_dirichletx_loy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_dirichletx_upy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_dirichletx_loy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_dirichletx_upy_dirichlety},}
    },
    }
  },
  // dirichletx-neumannx
  { .list = {
    // periodicy
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_neumannx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_loy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_upy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_dirichletx_loy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_dirichletx_upy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_neumannx_loy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_neumannx_upy_periodicy},
               {fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_neumannx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_loy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_upy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_dirichletx_loy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_dirichletx_upy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_neumannx_loy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_neumannx_upy_periodicy},},
    },
    // dirichlety-dirichlety
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_neumannx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_dirichletx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_dirichletx_upy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_neumannx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_neumannx_upy_dirichlety},
               {fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_neumannx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_dirichletx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_dirichletx_upy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_neumannx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_neumannx_upy_dirichlety},},
    },
    // dirichlety-neumanny
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_neumannx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_upy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_dirichletx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_dirichletx_upy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_neumannx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_neumannx_upy_neumanny},
               {fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_neumannx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_upy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_dirichletx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_dirichletx_upy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_neumannx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_neumannx_upy_neumanny},},
    },
    // neumanny-dirichlety
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_neumannx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_loy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_dirichletx_loy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_dirichletx_upy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_neumannx_loy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_neumannx_upy_dirichlety},
               {fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_neumannx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_loy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_dirichletx_loy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_dirichletx_upy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_neumannx_loy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_neumannx_upy_dirichlety},}
    },
    }
  },
  // neumannx-dirichletx
  { .list = {
    // periodicy
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_neumannx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_loy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_upy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_neumannx_loy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_neumannx_upy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_dirichletx_loy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_dirichletx_upy_periodicy},
               {fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_neumannx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_loy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_upy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_neumannx_loy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_neumannx_upy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_dirichletx_loy_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_dirichletx_upy_periodicy},},
    },
    // dirichlety-dirichlety
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_neumannx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_neumannx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_neumannx_upy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_dirichletx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_dirichletx_upy_dirichlety},
               {fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_neumannx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_neumannx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_neumannx_upy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_dirichletx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_dirichletx_upy_dirichlety},},
    },
    // dirichlety-neumanny
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_neumannx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_upy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_neumannx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_neumannx_upy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_dirichletx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_dirichletx_upy_neumanny},
               {fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_neumannx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_upy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_neumannx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_neumannx_upy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_dirichletx_loy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_dirichletx_upy_neumanny},},
    },
    // neumanny-dirichlety
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_neumannx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_loy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_neumannx_loy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p1_lox_neumannx_upy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_dirichletx_loy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p1_upx_dirichletx_upy_dirichlety},
               {fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_neumannx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_loy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p2_inx_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_neumannx_loy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p2_lox_neumannx_upy_dirichlety, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_dirichletx_loy_neumanny, fem_poisson_lhs_stencil_consteps_2x_ser_p2_upx_dirichletx_upy_dirichlety},}
    },
    }
  },
};

static const lhsstencil_kern_bcx_list_1x ser_lhsstencil_vareps_list_1x[] = {
  // periodicx
  { .list = {{NULL, NULL},
             {fem_poisson_lhs_stencil_vareps_1x_ser_p1_inx_periodicx, fem_poisson_lhs_stencil_vareps_1x_ser_p1_lox_periodicx, fem_poisson_lhs_stencil_vareps_1x_ser_p1_upx_periodicx},
             {fem_poisson_lhs_stencil_vareps_1x_ser_p2_inx_periodicx, fem_poisson_lhs_stencil_vareps_1x_ser_p2_lox_periodicx, fem_poisson_lhs_stencil_vareps_1x_ser_p2_upx_periodicx}}, },
  // dirichletx-dirichletx
  { .list = {{NULL, NULL},
             {fem_poisson_lhs_stencil_vareps_1x_ser_p1_inx_periodicx, fem_poisson_lhs_stencil_vareps_1x_ser_p1_lox_dirichletx, fem_poisson_lhs_stencil_vareps_1x_ser_p1_upx_dirichletx},
             {fem_poisson_lhs_stencil_vareps_1x_ser_p2_inx_periodicx, fem_poisson_lhs_stencil_vareps_1x_ser_p2_lox_dirichletx, fem_poisson_lhs_stencil_vareps_1x_ser_p2_upx_dirichletx}}, },
  // dirichletx-neumannx
  { .list = {{NULL, NULL},
             {fem_poisson_lhs_stencil_vareps_1x_ser_p1_inx_periodicx, fem_poisson_lhs_stencil_vareps_1x_ser_p1_lox_dirichletx, fem_poisson_lhs_stencil_vareps_1x_ser_p1_upx_neumannx},
             {fem_poisson_lhs_stencil_vareps_1x_ser_p2_inx_periodicx, fem_poisson_lhs_stencil_vareps_1x_ser_p2_lox_dirichletx, fem_poisson_lhs_stencil_vareps_1x_ser_p2_upx_neumannx}}, },
  // neumannx-dirichletx
  { .list = {{NULL, NULL},
             {fem_poisson_lhs_stencil_vareps_1x_ser_p1_inx_periodicx, fem_poisson_lhs_stencil_vareps_1x_ser_p1_lox_neumannx, fem_poisson_lhs_stencil_vareps_1x_ser_p1_upx_dirichletx},
             {fem_poisson_lhs_stencil_vareps_1x_ser_p2_inx_periodicx, fem_poisson_lhs_stencil_vareps_1x_ser_p2_lox_neumannx, fem_poisson_lhs_stencil_vareps_1x_ser_p2_upx_dirichletx}}, },
  // dirichletx-robinx
  { .list = {{NULL, NULL},
             {fem_poisson_lhs_stencil_vareps_1x_ser_p1_inx_periodicx, fem_poisson_lhs_stencil_vareps_1x_ser_p1_lox_dirichletx, fem_poisson_lhs_stencil_vareps_1x_ser_p1_upx_robinx},
             {fem_poisson_lhs_stencil_vareps_1x_ser_p2_inx_periodicx, fem_poisson_lhs_stencil_vareps_1x_ser_p2_lox_dirichletx, fem_poisson_lhs_stencil_vareps_1x_ser_p2_upx_robinx}}, },
  // robinx-dirichletx
  { .list = {{NULL, NULL},
             {fem_poisson_lhs_stencil_vareps_1x_ser_p1_inx_periodicx, fem_poisson_lhs_stencil_vareps_1x_ser_p1_lox_robinx, fem_poisson_lhs_stencil_vareps_1x_ser_p1_upx_dirichletx},
             {fem_poisson_lhs_stencil_vareps_1x_ser_p2_inx_periodicx, fem_poisson_lhs_stencil_vareps_1x_ser_p2_lox_robinx, fem_poisson_lhs_stencil_vareps_1x_ser_p2_upx_dirichletx}}, },
};

static const lhsstencil_kern_bcx_list_2x ser_lhsstencil_vareps_list_2x[] = {
  // periodicx
  { .list = {
    // periodicy
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_loy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_upy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_periodicx_loy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_periodicx_upy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_periodicx_loy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_periodicx_upy_periodicy},
               {fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_loy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_upy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_periodicx_loy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_periodicx_upy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_periodicx_loy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_periodicx_upy_periodicy},},
    },
    // dirichlety-dirichlety
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_periodicx_upy_dirichlety},
               {fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_periodicx_upy_dirichlety},},
    },
    // dirichlety-neumanny
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_upy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_periodicx_upy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_periodicx_upy_neumanny},
               {fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_upy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_periodicx_upy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_periodicx_upy_neumanny},},
    },
    // neumanny-dirichlety
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_loy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_periodicx_loy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_periodicx_loy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_periodicx_upy_dirichlety},
               {fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_loy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_periodicx_loy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_periodicx_loy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_periodicx_upy_dirichlety},}
    },
    }
  },
  // dirichletx-dirichletx
  { .list = {
    // periodicy
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_loy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_upy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_dirichletx_loy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_dirichletx_upy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_dirichletx_loy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_dirichletx_upy_periodicy},
               {fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_loy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_upy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_dirichletx_loy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_dirichletx_upy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_dirichletx_loy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_dirichletx_upy_periodicy},},
    },
    // dirichlety-dirichlety
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_dirichletx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_dirichletx_upy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_dirichletx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_dirichletx_upy_dirichlety},
               {fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_dirichletx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_dirichletx_upy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_dirichletx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_dirichletx_upy_dirichlety},},
    },
    // dirichlety-neumanny
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_upy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_dirichletx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_dirichletx_upy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_dirichletx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_dirichletx_upy_neumanny},
               {fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_upy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_dirichletx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_dirichletx_upy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_dirichletx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_dirichletx_upy_neumanny},},
    },
    // neumanny-dirichlety
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_loy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_dirichletx_loy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_dirichletx_upy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_dirichletx_loy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_dirichletx_upy_dirichlety},
               {fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_loy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_dirichletx_loy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_dirichletx_upy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_dirichletx_loy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_dirichletx_upy_dirichlety},}
    },
    }
  },
  // dirichletx-neumannx
  { .list = {
    // periodicy
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_neumannx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_loy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_upy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_dirichletx_loy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_dirichletx_upy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_neumannx_loy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_neumannx_upy_periodicy},
               {fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_neumannx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_loy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_upy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_dirichletx_loy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_dirichletx_upy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_neumannx_loy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_neumannx_upy_periodicy},},
    },
    // dirichlety-dirichlety
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_neumannx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_dirichletx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_dirichletx_upy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_neumannx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_neumannx_upy_dirichlety},
               {fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_neumannx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_dirichletx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_dirichletx_upy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_neumannx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_neumannx_upy_dirichlety},},
    },
    // dirichlety-neumanny
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_neumannx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_upy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_dirichletx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_dirichletx_upy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_neumannx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_neumannx_upy_neumanny},
               {fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_neumannx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_upy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_dirichletx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_dirichletx_upy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_neumannx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_neumannx_upy_neumanny},},
    },
    // neumanny-dirichlety
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_neumannx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_loy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_dirichletx_loy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_dirichletx_upy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_neumannx_loy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_neumannx_upy_dirichlety},
               {fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_neumannx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_loy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_dirichletx_loy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_dirichletx_upy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_neumannx_loy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_neumannx_upy_dirichlety},}
    },
    }
  },
  // neumannx-dirichletx
  { .list = {
    // periodicy
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_neumannx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_loy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_upy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_neumannx_loy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_neumannx_upy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_dirichletx_loy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_dirichletx_upy_periodicy},
               {fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_neumannx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_loy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_upy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_neumannx_loy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_neumannx_upy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_dirichletx_loy_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_dirichletx_upy_periodicy},},
    },
    // dirichlety-dirichlety
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_neumannx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_neumannx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_neumannx_upy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_dirichletx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_dirichletx_upy_dirichlety},
               {fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_neumannx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_neumannx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_neumannx_upy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_dirichletx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_dirichletx_upy_dirichlety},},
    },
    // dirichlety-neumanny
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_neumannx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_upy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_neumannx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_neumannx_upy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_dirichletx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_dirichletx_upy_neumanny},
               {fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_neumannx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_upy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_neumannx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_neumannx_upy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_dirichletx_loy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_dirichletx_upy_neumanny},},
    },
    // neumanny-dirichlety
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_neumannx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_loy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_neumannx_loy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p1_lox_neumannx_upy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_dirichletx_loy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p1_upx_dirichletx_upy_dirichlety},
               {fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_neumannx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_dirichletx_iny_periodicy, fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_loy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p2_inx_periodicx_upy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_neumannx_loy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p2_lox_neumannx_upy_dirichlety, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_dirichletx_loy_neumanny, fem_poisson_lhs_stencil_vareps_2x_ser_p2_upx_dirichletx_upy_dirichlety},}
    },
    }
  },
};

// Function pointer type for rhs source kernels.
typedef void (*srcstencil_t)(const double *epsilon, const double *dx, const double *rho, const double *bcVals, const long *globalIdxs,
  double *bsrc);

// For use in kernel tables.
typedef struct { srcstencil_t kernels[3]; } srcstencil_kern_loc_list_1x;
typedef struct { srcstencil_kern_loc_list_1x list[3]; } srcstencil_kern_bcx_list_1x;

typedef struct { srcstencil_t kernels[9]; } srcstencil_kern_loc_list_2x;
typedef struct { srcstencil_kern_loc_list_2x list[3]; } srcstencil_kern_bcy_list_2x;
typedef struct { srcstencil_kern_bcy_list_2x list[9]; } srcstencil_kern_bcx_list_2x;

// Serendipity src kernels.
GKYL_CU_D
static const srcstencil_kern_bcx_list_1x ser_srcstencil_consteps_list_1x[] = {
  // periodicx
  { .list = {{NULL, NULL},
             {fem_poisson_src_stencil_consteps_1x_ser_p1_inx_periodicx, fem_poisson_src_stencil_consteps_1x_ser_p1_lox_periodicx, fem_poisson_src_stencil_consteps_1x_ser_p1_upx_periodicx},
             {fem_poisson_src_stencil_consteps_1x_ser_p2_inx_periodicx, fem_poisson_src_stencil_consteps_1x_ser_p2_lox_periodicx, fem_poisson_src_stencil_consteps_1x_ser_p2_upx_periodicx}}, },
  // dirichletx-dirichletx
  { .list = {{NULL, NULL},
             {fem_poisson_src_stencil_consteps_1x_ser_p1_inx_periodicx, fem_poisson_src_stencil_consteps_1x_ser_p1_lox_dirichletx, fem_poisson_src_stencil_consteps_1x_ser_p1_upx_dirichletx},
             {fem_poisson_src_stencil_consteps_1x_ser_p2_inx_periodicx, fem_poisson_src_stencil_consteps_1x_ser_p2_lox_dirichletx, fem_poisson_src_stencil_consteps_1x_ser_p2_upx_dirichletx}}, },
  // dirichletx-neumannx
  { .list = {{NULL, NULL},
             {fem_poisson_src_stencil_consteps_1x_ser_p1_inx_periodicx, fem_poisson_src_stencil_consteps_1x_ser_p1_lox_dirichletx, fem_poisson_src_stencil_consteps_1x_ser_p1_upx_neumannx},
             {fem_poisson_src_stencil_consteps_1x_ser_p2_inx_periodicx, fem_poisson_src_stencil_consteps_1x_ser_p2_lox_dirichletx, fem_poisson_src_stencil_consteps_1x_ser_p2_upx_neumannx}}, },
  // neumannx-dirichletx
  { .list = {{NULL, NULL},
             {fem_poisson_src_stencil_consteps_1x_ser_p1_inx_periodicx, fem_poisson_src_stencil_consteps_1x_ser_p1_lox_neumannx, fem_poisson_src_stencil_consteps_1x_ser_p1_upx_dirichletx},
             {fem_poisson_src_stencil_consteps_1x_ser_p2_inx_periodicx, fem_poisson_src_stencil_consteps_1x_ser_p2_lox_neumannx, fem_poisson_src_stencil_consteps_1x_ser_p2_upx_dirichletx}}, },
  // dirichletx-robinx
  { .list = {{NULL, NULL},
             {fem_poisson_src_stencil_consteps_1x_ser_p1_inx_periodicx, fem_poisson_src_stencil_consteps_1x_ser_p1_lox_dirichletx, fem_poisson_src_stencil_consteps_1x_ser_p1_upx_robinx},
             {fem_poisson_src_stencil_consteps_1x_ser_p2_inx_periodicx, fem_poisson_src_stencil_consteps_1x_ser_p2_lox_dirichletx, fem_poisson_src_stencil_consteps_1x_ser_p2_upx_robinx}}, },
  // robinx-dirichletx
  { .list = {{NULL, NULL},
             {fem_poisson_src_stencil_consteps_1x_ser_p1_inx_periodicx, fem_poisson_src_stencil_consteps_1x_ser_p1_lox_robinx, fem_poisson_src_stencil_consteps_1x_ser_p1_upx_dirichletx},
             {fem_poisson_src_stencil_consteps_1x_ser_p2_inx_periodicx, fem_poisson_src_stencil_consteps_1x_ser_p2_lox_robinx, fem_poisson_src_stencil_consteps_1x_ser_p2_upx_dirichletx}}, },
};

GKYL_CU_D
static const srcstencil_kern_bcx_list_2x ser_srcstencil_consteps_list_2x[] = {
  // periodicx
  { .list = {
    // periodicy
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_loy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_upy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_periodicx_loy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_periodicx_upy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_periodicx_loy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_periodicx_upy_periodicy},
               {fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_loy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_upy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_periodicx_loy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_periodicx_upy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_periodicx_loy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_periodicx_upy_periodicy},},
    },
    // dirichlety-dirichlety
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_periodicx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_periodicx_upy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_periodicx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_periodicx_upy_dirichlety},
               {fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_upy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_periodicx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_periodicx_upy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_periodicx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_periodicx_upy_dirichlety},},
    },
    // dirichlety-neumanny
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_upy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_periodicx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_periodicx_upy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_periodicx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_periodicx_upy_neumanny},
               {fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_upy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_periodicx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_periodicx_upy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_periodicx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_periodicx_upy_neumanny},},
    },
    // neumanny-dirichlety
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_loy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_periodicx_loy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_periodicx_upy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_periodicx_loy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_periodicx_upy_dirichlety},
               {fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_loy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_upy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_periodicx_loy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_periodicx_upy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_periodicx_loy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_periodicx_upy_dirichlety},}
    },
    }
  },
  // dirichletx-dirichletx
  { .list = {
    // periodicy
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_loy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_upy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_dirichletx_loy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_dirichletx_upy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_dirichletx_loy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_dirichletx_upy_periodicy},
               {fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_dirichletx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_dirichletx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_loy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_upy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_dirichletx_loy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_dirichletx_upy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_dirichletx_loy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_dirichletx_upy_periodicy},},
    },
    // dirichlety-dirichlety
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_dirichletx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_dirichletx_upy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_dirichletx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_dirichletx_upy_dirichlety},
               {fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_dirichletx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_dirichletx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_upy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_dirichletx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_dirichletx_upy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_dirichletx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_dirichletx_upy_dirichlety},},
    },
    // dirichlety-neumanny
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_upy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_dirichletx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_dirichletx_upy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_dirichletx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_dirichletx_upy_neumanny},
               {fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_dirichletx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_dirichletx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_upy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_dirichletx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_dirichletx_upy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_dirichletx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_dirichletx_upy_neumanny},},
    },
    // neumanny-dirichlety
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_loy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_dirichletx_loy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_dirichletx_upy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_dirichletx_loy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_dirichletx_upy_dirichlety},
               {fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_dirichletx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_dirichletx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_loy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_upy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_dirichletx_loy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_dirichletx_upy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_dirichletx_loy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_dirichletx_upy_dirichlety},}
    },
    }
  },
  // dirichletx-neumannx
  { .list = {
    // periodicy
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_neumannx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_loy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_upy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_dirichletx_loy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_dirichletx_upy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_neumannx_loy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_neumannx_upy_periodicy},
               {fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_dirichletx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_neumannx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_loy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_upy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_dirichletx_loy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_dirichletx_upy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_neumannx_loy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_neumannx_upy_periodicy},},
    },
    // dirichlety-dirichlety
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_neumannx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_dirichletx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_dirichletx_upy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_neumannx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_neumannx_upy_dirichlety},
               {fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_dirichletx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_neumannx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_upy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_dirichletx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_dirichletx_upy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_neumannx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_neumannx_upy_dirichlety},},
    },
    // dirichlety-neumanny
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_neumannx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_upy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_dirichletx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_dirichletx_upy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_neumannx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_neumannx_upy_neumanny},
               {fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_dirichletx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_neumannx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_upy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_dirichletx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_dirichletx_upy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_neumannx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_neumannx_upy_neumanny},},
    },
    // neumanny-dirichlety
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_neumannx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_loy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_dirichletx_loy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_dirichletx_upy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_neumannx_loy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_neumannx_upy_dirichlety},
               {fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_dirichletx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_neumannx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_loy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_upy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_dirichletx_loy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_dirichletx_upy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_neumannx_loy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_neumannx_upy_dirichlety},}
    },
    }
  },
  // neumannx-dirichletx
  { .list = {
    // periodicy
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_neumannx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_loy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_upy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_neumannx_loy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_neumannx_upy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_dirichletx_loy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_dirichletx_upy_periodicy},
               {fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_neumannx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_dirichletx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_loy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_upy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_neumannx_loy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_neumannx_upy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_dirichletx_loy_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_dirichletx_upy_periodicy},},
    },
    // dirichlety-dirichlety
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_neumannx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_neumannx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_neumannx_upy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_dirichletx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_dirichletx_upy_dirichlety},
               {fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_neumannx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_dirichletx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_upy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_neumannx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_neumannx_upy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_dirichletx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_dirichletx_upy_dirichlety},},
    },
    // dirichlety-neumanny
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_neumannx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_upy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_neumannx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_neumannx_upy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_dirichletx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_dirichletx_upy_neumanny},
               {fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_neumannx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_dirichletx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_upy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_neumannx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_neumannx_upy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_dirichletx_loy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_dirichletx_upy_neumanny},},
    },
    // neumanny-dirichlety
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_neumannx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_loy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_neumannx_loy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p1_lox_neumannx_upy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_dirichletx_loy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p1_upx_dirichletx_upy_dirichlety},
               {fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_neumannx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_dirichletx_iny_periodicy, fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_loy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p2_inx_periodicx_upy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_neumannx_loy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p2_lox_neumannx_upy_dirichlety, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_dirichletx_loy_neumanny, fem_poisson_src_stencil_consteps_2x_ser_p2_upx_dirichletx_upy_dirichlety},}
    },
    }
  },
};

GKYL_CU_D
static const srcstencil_kern_bcx_list_1x ser_srcstencil_vareps_list_1x[] = {
  // periodicx
  { .list = {{NULL, NULL},
             {fem_poisson_src_stencil_vareps_1x_ser_p1_inx_periodicx, fem_poisson_src_stencil_vareps_1x_ser_p1_lox_periodicx, fem_poisson_src_stencil_vareps_1x_ser_p1_upx_periodicx},
             {fem_poisson_src_stencil_vareps_1x_ser_p2_inx_periodicx, fem_poisson_src_stencil_vareps_1x_ser_p2_lox_periodicx, fem_poisson_src_stencil_vareps_1x_ser_p2_upx_periodicx}}, },
  // dirichletx-dirichletx
  { .list = {{NULL, NULL},
             {fem_poisson_src_stencil_vareps_1x_ser_p1_inx_periodicx, fem_poisson_src_stencil_vareps_1x_ser_p1_lox_dirichletx, fem_poisson_src_stencil_vareps_1x_ser_p1_upx_dirichletx},
             {fem_poisson_src_stencil_vareps_1x_ser_p2_inx_periodicx, fem_poisson_src_stencil_vareps_1x_ser_p2_lox_dirichletx, fem_poisson_src_stencil_vareps_1x_ser_p2_upx_dirichletx}}, },
  // dirichletx-neumannx
  { .list = {{NULL, NULL},
             {fem_poisson_src_stencil_vareps_1x_ser_p1_inx_periodicx, fem_poisson_src_stencil_vareps_1x_ser_p1_lox_dirichletx, fem_poisson_src_stencil_vareps_1x_ser_p1_upx_neumannx},
             {fem_poisson_src_stencil_vareps_1x_ser_p2_inx_periodicx, fem_poisson_src_stencil_vareps_1x_ser_p2_lox_dirichletx, fem_poisson_src_stencil_vareps_1x_ser_p2_upx_neumannx}}, },
  // neumannx-dirichletx
  { .list = {{NULL, NULL},
             {fem_poisson_src_stencil_vareps_1x_ser_p1_inx_periodicx, fem_poisson_src_stencil_vareps_1x_ser_p1_lox_neumannx, fem_poisson_src_stencil_vareps_1x_ser_p1_upx_dirichletx},
             {fem_poisson_src_stencil_vareps_1x_ser_p2_inx_periodicx, fem_poisson_src_stencil_vareps_1x_ser_p2_lox_neumannx, fem_poisson_src_stencil_vareps_1x_ser_p2_upx_dirichletx}}, },
  // dirichletx-robinx
  { .list = {{NULL, NULL},
             {fem_poisson_src_stencil_vareps_1x_ser_p1_inx_periodicx, fem_poisson_src_stencil_vareps_1x_ser_p1_lox_dirichletx, fem_poisson_src_stencil_vareps_1x_ser_p1_upx_robinx},
             {fem_poisson_src_stencil_vareps_1x_ser_p2_inx_periodicx, fem_poisson_src_stencil_vareps_1x_ser_p2_lox_dirichletx, fem_poisson_src_stencil_vareps_1x_ser_p2_upx_robinx}}, },
  // robinx-dirichletx
  { .list = {{NULL, NULL},
             {fem_poisson_src_stencil_vareps_1x_ser_p1_inx_periodicx, fem_poisson_src_stencil_vareps_1x_ser_p1_lox_robinx, fem_poisson_src_stencil_vareps_1x_ser_p1_upx_dirichletx},
             {fem_poisson_src_stencil_vareps_1x_ser_p2_inx_periodicx, fem_poisson_src_stencil_vareps_1x_ser_p2_lox_robinx, fem_poisson_src_stencil_vareps_1x_ser_p2_upx_dirichletx}}, },
};

GKYL_CU_D
static const srcstencil_kern_bcx_list_2x ser_srcstencil_vareps_list_2x[] = {
  // periodicx
  { .list = {
    // periodicy
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_loy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_upy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_periodicx_loy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_periodicx_upy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_periodicx_loy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_periodicx_upy_periodicy},
               {fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_loy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_upy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_periodicx_loy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_periodicx_upy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_periodicx_loy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_periodicx_upy_periodicy},},
    },
    // dirichlety-dirichlety
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_periodicx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_periodicx_upy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_periodicx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_periodicx_upy_dirichlety},
               {fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_upy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_periodicx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_periodicx_upy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_periodicx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_periodicx_upy_dirichlety},},
    },
    // dirichlety-neumanny
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_upy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_periodicx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_periodicx_upy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_periodicx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_periodicx_upy_neumanny},
               {fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_upy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_periodicx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_periodicx_upy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_periodicx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_periodicx_upy_neumanny},},
    },
    // neumanny-dirichlety
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_loy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_periodicx_loy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_periodicx_upy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_periodicx_loy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_periodicx_upy_dirichlety},
               {fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_loy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_upy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_periodicx_loy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_periodicx_upy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_periodicx_loy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_periodicx_upy_dirichlety},}
    },
    }
  },
  // dirichletx-dirichletx
  { .list = {
    // periodicy
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_loy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_upy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_dirichletx_loy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_dirichletx_upy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_dirichletx_loy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_dirichletx_upy_periodicy},
               {fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_dirichletx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_dirichletx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_loy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_upy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_dirichletx_loy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_dirichletx_upy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_dirichletx_loy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_dirichletx_upy_periodicy},},
    },
    // dirichlety-dirichlety
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_dirichletx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_dirichletx_upy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_dirichletx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_dirichletx_upy_dirichlety},
               {fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_dirichletx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_dirichletx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_upy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_dirichletx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_dirichletx_upy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_dirichletx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_dirichletx_upy_dirichlety},},
    },
    // dirichlety-neumanny
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_upy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_dirichletx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_dirichletx_upy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_dirichletx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_dirichletx_upy_neumanny},
               {fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_dirichletx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_dirichletx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_upy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_dirichletx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_dirichletx_upy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_dirichletx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_dirichletx_upy_neumanny},},
    },
    // neumanny-dirichlety
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_loy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_dirichletx_loy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_dirichletx_upy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_dirichletx_loy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_dirichletx_upy_dirichlety},
               {fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_dirichletx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_dirichletx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_loy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_upy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_dirichletx_loy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_dirichletx_upy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_dirichletx_loy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_dirichletx_upy_dirichlety},}
    },
    }
  },
  // dirichletx-neumannx
  { .list = {
    // periodicy
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_neumannx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_loy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_upy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_dirichletx_loy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_dirichletx_upy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_neumannx_loy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_neumannx_upy_periodicy},
               {fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_dirichletx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_neumannx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_loy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_upy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_dirichletx_loy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_dirichletx_upy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_neumannx_loy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_neumannx_upy_periodicy},},
    },
    // dirichlety-dirichlety
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_neumannx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_dirichletx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_dirichletx_upy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_neumannx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_neumannx_upy_dirichlety},
               {fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_dirichletx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_neumannx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_upy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_dirichletx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_dirichletx_upy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_neumannx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_neumannx_upy_dirichlety},},
    },
    // dirichlety-neumanny
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_neumannx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_upy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_dirichletx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_dirichletx_upy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_neumannx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_neumannx_upy_neumanny},
               {fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_dirichletx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_neumannx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_upy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_dirichletx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_dirichletx_upy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_neumannx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_neumannx_upy_neumanny},},
    },
    // neumanny-dirichlety
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_neumannx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_loy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_dirichletx_loy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_dirichletx_upy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_neumannx_loy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_neumannx_upy_dirichlety},
               {fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_dirichletx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_neumannx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_loy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_upy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_dirichletx_loy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_dirichletx_upy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_neumannx_loy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_neumannx_upy_dirichlety},}
    },
    }
  },
  // neumannx-dirichletx
  { .list = {
    // periodicy
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_neumannx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_loy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_upy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_neumannx_loy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_neumannx_upy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_dirichletx_loy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_dirichletx_upy_periodicy},
               {fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_neumannx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_dirichletx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_loy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_upy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_neumannx_loy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_neumannx_upy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_dirichletx_loy_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_dirichletx_upy_periodicy},},
    },
    // dirichlety-dirichlety
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_neumannx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_neumannx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_neumannx_upy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_dirichletx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_dirichletx_upy_dirichlety},
               {fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_neumannx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_dirichletx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_upy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_neumannx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_neumannx_upy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_dirichletx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_dirichletx_upy_dirichlety},},
    },
    // dirichlety-neumanny
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_neumannx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_upy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_neumannx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_neumannx_upy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_dirichletx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_dirichletx_upy_neumanny},
               {fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_neumannx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_dirichletx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_upy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_neumannx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_neumannx_upy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_dirichletx_loy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_dirichletx_upy_neumanny},},
    },
    // neumanny-dirichlety
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_neumannx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_loy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_neumannx_loy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p1_lox_neumannx_upy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_dirichletx_loy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p1_upx_dirichletx_upy_dirichlety},
               {fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_neumannx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_dirichletx_iny_periodicy, fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_loy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p2_inx_periodicx_upy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_neumannx_loy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p2_lox_neumannx_upy_dirichlety, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_dirichletx_loy_neumanny, fem_poisson_src_stencil_vareps_2x_ser_p2_upx_dirichletx_upy_dirichlety},}
    },
    }
  },
};

// Function pointer type for sol kernels.
typedef void (*solstencil_t)(const double *sol_nodal_global, const long *globalIdxs, double *sol_modal_local);

typedef struct { solstencil_t kernels[3]; } solstencil_kern_list;

GKYL_CU_D
static const solstencil_kern_list ser_solstencil_list[] = {
  { NULL, NULL, NULL },
  // 1x kernels
  { NULL, fem_poisson_sol_stencil_1x_ser_p1, fem_poisson_sol_stencil_1x_ser_p2 }, // 0
  // 2x kernels
  { NULL, fem_poisson_sol_stencil_2x_ser_p1, fem_poisson_sol_stencil_2x_ser_p2 }, // 1
};

// "Choose Kernel" based on cdim, vdim and polyorder
#define CK1(lst,poly_order,loc,bcx) lst[bcx].list[poly_order].kernels[loc]
#define CK2(lst,poly_order,loc,bcx,bcy) lst[bcx].list[bcy].list[poly_order].kernels[loc]
#define CK3(lst,poly_order,loc,bcx,bcy,bcz) lst[bcx].list[bcy].list[bcz].list[poly_order].kernels[loc]

// Struct containing pointers to the various kernels. Needed to create a similar struct on the GPU.
struct gkyl_fem_poisson_kernels { 
  // Pointer to local-to-global kernels. 2^3, 2 (interior and upper) in each direction.
  local2global_t l2g[8];

  // LHS kernels (one for each position in the domain, up to 3^3).
  lhsstencil_t lhsker[27];

  // RHS source kernels (one for each position in the domain, up to 3^3).
  srcstencil_t srcker[27];

  solstencil_t solker;
};

// Updater type
struct gkyl_fem_poisson {
  void *ctx; // evaluation context.
  struct gkyl_rect_grid grid;
  int ndim; // grid's number of dimensions.
  int num_basis; // number of basis functions.
  enum gkyl_basis_type basis_type;
  int poly_order;
  struct gkyl_basis basis;
  int num_cells[GKYL_MAX_CDIM];
  double dx[GKYL_MAX_CDIM];
#ifdef GKYL_HAVE_CUDA
  double *dx_cu;
#endif
  bool isdirperiodic[GKYL_MAX_CDIM]; // =true if direction is periodic.

  struct gkyl_array *epsilon; // permittivity.
  bool isvareps; // indicate if permittivity is a tensor/varies in space.
  bool ishelmholtz; // if solving Helmholtz equation (kSq is not zero/NULL).

  bool isdomperiodic; // =true if all directions are periodic.
  struct gkyl_array *rhs_cellavg;
  double *rhs_avg, mavgfac;
  double *rhs_avg_cu;

  double bcvals[GKYL_MAX_CDIM*2*3]; // BC values, bc[0]*phi+bc[1]*d(phi)/dx=phi[3] at each boundary.
  double *bcvals_cu; // BC values, bc[0]*phi+bc[1]*d(phi)/dx=phi[3] at each boundary.

  struct gkyl_range local_range, local_range_ext;
  struct gkyl_range solve_range;
  struct gkyl_range_iter solve_iter;

  int numnodes_local;
  long numnodes_global;

  struct gkyl_superlu_prob* prob;
  struct gkyl_array *brhs;

#ifdef GKYL_HAVE_CUDA
  struct gkyl_cusolver_prob *prob_cu;
  struct gkyl_array *brhs_cu;
#endif

  long *globalidx;

  struct gkyl_fem_poisson_kernels *kernels;
  struct gkyl_fem_poisson_kernels *kernels_cu;
  bool use_gpu;
};

void
fem_poisson_choose_kernels_cu(const struct gkyl_basis* basis, const struct gkyl_poisson_bc* bcs, bool isvareps, const bool *isdirperiodic, struct gkyl_fem_poisson_kernels *kers);

static long
gkyl_fem_poisson_global_num_nodes(const int dim, const int poly_order, const int basis_type, const int *num_cells, bool *isdirperiodic)
{
  if (dim==1) {
    if (poly_order == 1) {
      if (isdirperiodic[0]) {
        return fem_poisson_num_nodes_global_1x_ser_p1_periodicx(num_cells);
      } else {
        return fem_poisson_num_nodes_global_1x_ser_p1_nonperiodicx(num_cells);
      }
    } else if (poly_order == 2) {
      if (isdirperiodic[0]) {
        return fem_poisson_num_nodes_global_1x_ser_p2_periodicx(num_cells);
      } else {
        return fem_poisson_num_nodes_global_1x_ser_p2_nonperiodicx(num_cells);
      }
    }
  } else if (dim==2) {
    if (poly_order == 1) {
      if (isdirperiodic[0] && isdirperiodic[1]) {
        return fem_poisson_num_nodes_global_2x_ser_p1_periodicx_periodicy(num_cells);
      } else if (!isdirperiodic[0] && isdirperiodic[1]) {
        return fem_poisson_num_nodes_global_2x_ser_p1_nonperiodicx_periodicy(num_cells);
      } else if (isdirperiodic[0] && !isdirperiodic[1]) {
        return fem_poisson_num_nodes_global_2x_ser_p1_periodicx_nonperiodicy(num_cells);
      } else {
        return fem_poisson_num_nodes_global_2x_ser_p1_nonperiodicx_nonperiodicy(num_cells);
      }
    } else if (poly_order == 2) {
      if (isdirperiodic[0] && isdirperiodic[1]) {
        return fem_poisson_num_nodes_global_2x_ser_p2_periodicx_periodicy(num_cells);
      } else if (!isdirperiodic[0] && isdirperiodic[1]) {
        return fem_poisson_num_nodes_global_2x_ser_p2_nonperiodicx_periodicy(num_cells);
      } else if (isdirperiodic[0] && !isdirperiodic[1]) {
        return fem_poisson_num_nodes_global_2x_ser_p2_periodicx_nonperiodicy(num_cells);
      } else {
        return fem_poisson_num_nodes_global_2x_ser_p2_nonperiodicx_nonperiodicy(num_cells);
      }
    }
  } else if (dim==3) {
    assert(false);  // Other dimensionalities not supported.
  }
  assert(false);  // Other dimensionalities not supported.
  return -1;
}

GKYL_CU_D
static void
fem_poisson_choose_local2global_kernels(const struct gkyl_basis* basis, const bool *isdirperiodic, local2global_t *l2gout)
{
  int dim = basis->ndim;
  int poly_order = basis->poly_order;

  int bckey[GKYL_MAX_CDIM] = {-1};
  for (int d=0; d<basis->ndim; d++) bckey[d] = isdirperiodic[d] ? 0 : 1;

  switch (basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      for (int k=0; k<(int)(pow(2,dim)+0.5); k++) {
        if (dim == 1) {
          l2gout[k] = CK1(ser_loc2glob_list_1x, poly_order, k, bckey[0]);
        } else if (dim == 2) {
          l2gout[k] = CK2(ser_loc2glob_list_2x, poly_order, k, bckey[0], bckey[1]);
//        } else if (dim == 3) {
//          l2gout[k] = CK3(ser_loc2glob_list_3x, poly_order, k, bckey[0], bckey[1], bckey[2]);
        }
      }
      break;
//    case GKYL_BASIS_MODAL_TENSOR:
//      break;
    default:
      assert(false);
      break;
  }
}

GKYL_CU_D
static void
fem_poisson_choose_lhs_kernels(const struct gkyl_basis* basis, const struct gkyl_poisson_bc *bcs, bool isvareps, lhsstencil_t *lhsout)
{
  int dim = basis->ndim;
  int poly_order = basis->poly_order;

  int bckey[GKYL_MAX_CDIM] = {-1};
  for (int d=0; d<basis->ndim; d++) {
    if (bcs->lo_type[d]==GKYL_POISSON_PERIODIC && bcs->up_type[d]==GKYL_POISSON_PERIODIC) { bckey[d] = 0; }
    else if (bcs->lo_type[d]==GKYL_POISSON_DIRICHLET && bcs->up_type[d]==GKYL_POISSON_DIRICHLET) { bckey[d] = 1; }
    else if (bcs->lo_type[d]==GKYL_POISSON_DIRICHLET && bcs->up_type[d]==GKYL_POISSON_NEUMANN) { bckey[d] = 2; }
    else if (bcs->lo_type[d]==GKYL_POISSON_NEUMANN && bcs->up_type[d]==GKYL_POISSON_DIRICHLET) { bckey[d] = 3; }
    else if (bcs->lo_type[d]==GKYL_POISSON_DIRICHLET && bcs->up_type[d]==GKYL_POISSON_ROBIN) { bckey[d] = 4; }
    else if (bcs->lo_type[d]==GKYL_POISSON_ROBIN && bcs->up_type[d]==GKYL_POISSON_DIRICHLET) { bckey[d] = 5; }
    else { assert(false); }
  };

  switch (basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      for (int k=0; k<(int)(pow(3,dim)+0.5); k++) {
        if (dim == 1) {
          lhsout[k] = isvareps?  CK1(ser_lhsstencil_vareps_list_1x, poly_order, k, bckey[0])
                               : CK1(ser_lhsstencil_consteps_list_1x, poly_order, k, bckey[0]);
        } else if (dim == 2) {                              
          lhsout[k] = isvareps?  CK2(ser_lhsstencil_vareps_list_2x, poly_order, k, bckey[0], bckey[1])
                               : CK2(ser_lhsstencil_consteps_list_2x, poly_order, k, bckey[0], bckey[1]);
//      } else if (dim == 3) {                              
//        lhsout[k] = isvareps?  CK3(ser_lhsstencil_vareps_list_3x, poly_order, k, bckey[0], bckey[1], bckey[2])
//                             : CK3(ser_lhsstencil_consteps_list_3x, poly_order, k, bckey[0], bckey[1], bckey[2]);
        }
      }
      break;
//    case GKYL_BASIS_MODAL_TENSOR:
//      break;
    default:
      assert(false);
      break;
  }
}

GKYL_CU_D
static void
fem_poisson_choose_src_kernels(const struct gkyl_basis* basis, const struct gkyl_poisson_bc *bcs, bool isvareps, srcstencil_t *srcout)
{
  int dim = basis->ndim;
  int poly_order = basis->poly_order;

  int bckey[GKYL_MAX_CDIM] = {-1};
  for (int d=0; d<basis->ndim; d++) {
    if (bcs->lo_type[d]==GKYL_POISSON_PERIODIC && bcs->up_type[d]==GKYL_POISSON_PERIODIC) { bckey[d] = 0; }
    else if (bcs->lo_type[d]==GKYL_POISSON_DIRICHLET && bcs->up_type[d]==GKYL_POISSON_DIRICHLET) { bckey[d] = 1; }
    else if (bcs->lo_type[d]==GKYL_POISSON_DIRICHLET && bcs->up_type[d]==GKYL_POISSON_NEUMANN) { bckey[d] = 2; }
    else if (bcs->lo_type[d]==GKYL_POISSON_NEUMANN && bcs->up_type[d]==GKYL_POISSON_DIRICHLET) { bckey[d] = 3; }
    else if (bcs->lo_type[d]==GKYL_POISSON_DIRICHLET && bcs->up_type[d]==GKYL_POISSON_ROBIN) { bckey[d] = 4; }
    else if (bcs->lo_type[d]==GKYL_POISSON_ROBIN && bcs->up_type[d]==GKYL_POISSON_DIRICHLET) { bckey[d] = 5; }
    else { assert(false); }
  };

  switch (basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      for (int k=0; k<(int)(pow(3,dim)+0.5); k++) {
        if (dim == 1) {
          srcout[k] = isvareps?  CK1(ser_srcstencil_vareps_list_1x, poly_order, k, bckey[0])
                               : CK1(ser_srcstencil_consteps_list_1x, poly_order, k, bckey[0]);
        } else if (dim == 2) {                              
          srcout[k] = isvareps?  CK2(ser_srcstencil_vareps_list_2x, poly_order, k, bckey[0], bckey[1])
                               : CK2(ser_srcstencil_consteps_list_2x, poly_order, k, bckey[0], bckey[1]);
//        } else if (dim == 3) {                            
//          srcout[k] = isvareps?  CK3(ser_srcstencil_vareps_list_3x, poly_order, k, bckey[0], bckey[1], bckey[2])
//                               : CK3(ser_srcstencil_consteps_list_3x, poly_order, k, bckey[0], bckey[1], bckey[2]);
        }
      }
      break;
//    case GKYL_BASIS_MODAL_TENSOR:
//      break;
    default:
      assert(false);
      break;
  }
}

GKYL_CU_D
static solstencil_t
fem_poisson_choose_sol_kernels(const struct gkyl_basis* basis)
{
  int dim = basis->ndim;
  int poly_order = basis->poly_order;

  switch (basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_solstencil_list[dim].kernels[poly_order];
      
      break;
//    case GKYL_BASIS_MODAL_TENSOR:
//      break;
    default:
      assert(false);
      break;
  }
}
