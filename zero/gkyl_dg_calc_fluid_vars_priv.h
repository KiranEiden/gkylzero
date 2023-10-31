// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_euler_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <assert.h>

typedef int (*fluid_set_t)(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *fluid);

typedef void (*fluid_copy_t)(int count, struct gkyl_nmat *x, 
  double* GKYL_RESTRICT u, double* GKYL_RESTRICT u_surf);

typedef void (*fluid_pressure_t)(double param, const double *fluid, const double *u, 
  double* GKYL_RESTRICT p, double* GKYL_RESTRICT p_surf);

typedef void (*fluid_limiter_t)(double gas_gamma, double *p_c, 
  double *fluid_l, double *fluid_c, double *fluid_r);

// for use in kernel tables
typedef struct { fluid_set_t kernels[3]; } gkyl_dg_fluid_set_kern_list;
typedef struct { fluid_copy_t kernels[3]; } gkyl_dg_fluid_copy_kern_list;
typedef struct { fluid_pressure_t kernels[3]; } gkyl_dg_fluid_pressure_kern_list;
typedef struct { fluid_limiter_t kernels[3]; } gkyl_dg_fluid_limiter_kern_list;

struct gkyl_dg_calc_fluid_vars {
  struct gkyl_rect_grid conf_grid; // Configuration space grid for cell spacing and cell center
  int cdim; // Configuration space dimensionality
  int poly_order; // polynomial order (determines whether we solve linear system or use basis_inv method)
  struct gkyl_range mem_range; // Configuration space range for linear solve

  struct gkyl_nmat *As, *xs; // matrices for LHS and RHS
  gkyl_nmat_mem *mem; // memory for use in batched linear solve
  int Ncomp; // number of components in the linear solve (6 variables being solved for)

  fluid_set_t fluid_set;  // kernel for setting matrices for linear solve
  fluid_copy_t fluid_copy; // kernel for copying solution to output; also computed needed surface expansions
  fluid_pressure_t fluid_pressure; // kernel for computing pressure (Volume and surface expansion)
  fluid_limiter_t fluid_limiter; // kernel for limiting slopes of fluid variables

  uint32_t flags;
  struct gkyl_dg_calc_fluid_vars *on_dev; // pointer to itself or device data
};

// Set matrices for computing fluid flow velocity (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_fluid_set_kern_list ser_fluid_set_kernels[] = {
  { NULL, fluid_vars_u_set_1x_ser_p1, fluid_vars_u_set_1x_ser_p2 }, // 0
  { NULL, fluid_vars_u_set_2x_ser_p1, NULL }, // 1
  { NULL, fluid_vars_u_set_3x_ser_p1, NULL }, // 2
};

// Set matrices for computing fluid flow velocity (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_fluid_set_kern_list ten_fluid_set_kernels[] = {
  { NULL, fluid_vars_u_set_1x_ser_p1, fluid_vars_u_set_1x_ser_p2 }, // 0
  { NULL, fluid_vars_u_set_2x_ser_p1, fluid_vars_u_set_2x_tensor_p2 }, // 1
  { NULL, fluid_vars_u_set_3x_ser_p1, NULL }, // 2
};

// Copy solution for fluid flow velocity (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_fluid_copy_kern_list ser_fluid_copy_kernels[] = {
  { NULL, fluid_vars_u_copy_1x_ser_p1, fluid_vars_u_copy_1x_ser_p2 }, // 0
  { NULL, fluid_vars_u_copy_2x_ser_p1, NULL }, // 1
  { NULL, fluid_vars_u_copy_3x_ser_p1, NULL }, // 2
};

// Copy solution for fluid flow velocity (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_fluid_copy_kern_list ten_fluid_copy_kernels[] = {
  { NULL, fluid_vars_u_copy_1x_ser_p1, fluid_vars_u_copy_1x_ser_p2 }, // 0
  { NULL, fluid_vars_u_copy_2x_ser_p1, fluid_vars_u_copy_2x_tensor_p2 }, // 1
  { NULL, fluid_vars_u_copy_3x_ser_p1, NULL }, // 2
};

// Scalar pressure Isothermal Euler -> p = vth*rho; Euler -> p = (gas_gamma - 1)*(E - 1/2 rho u^2) (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_fluid_pressure_kern_list ser_fluid_pressure_kernels[] = {
  { NULL, fluid_vars_pressure_1x_ser_p1, fluid_vars_pressure_1x_ser_p2 }, // 0
  { NULL, fluid_vars_pressure_2x_ser_p1, NULL }, // 1
  { NULL, fluid_vars_pressure_3x_ser_p1, NULL }, // 2
};

// Scalar pressure Isothermal Euler -> p = vth*rho; Euler -> p = (gas_gamma - 1)*(E - 1/2 rho u^2) (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_fluid_pressure_kern_list ten_fluid_pressure_kernels[] = {
  { NULL, fluid_vars_pressure_1x_ser_p1, fluid_vars_pressure_1x_ser_p2 }, // 0
  { NULL, fluid_vars_pressure_2x_ser_p1, fluid_vars_pressure_2x_tensor_p2 }, // 1
  { NULL, fluid_vars_pressure_3x_ser_p1, NULL }, // 2
};

// Scalar pressure Isothermal Euler -> p = vth*rho; Euler -> p = (gas_gamma - 1)*(E - 1/2 rho u^2) (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_fluid_limiter_kern_list ser_fluid_limiter_kernels[] = {
  { NULL, fluid_vars_limiter_1x_ser_p1, fluid_vars_limiter_1x_ser_p2 }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
};

// Scalar pressure Isothermal Euler -> p = vth*rho; Euler -> p = (gas_gamma - 1)*(E - 1/2 rho u^2) (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_fluid_limiter_kern_list ten_fluid_limiter_kernels[] = {
  { NULL, fluid_vars_limiter_1x_ser_p1, fluid_vars_limiter_1x_ser_p2 }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
};

GKYL_CU_D
static fluid_set_t
choose_fluid_set_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_fluid_set_kernels[cdim-1].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return ten_fluid_set_kernels[cdim-1].kernels[poly_order];
      break;
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static fluid_copy_t
choose_fluid_copy_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_fluid_copy_kernels[cdim-1].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return ten_fluid_copy_kernels[cdim-1].kernels[poly_order];
      break;
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static fluid_pressure_t
choose_fluid_pressure_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_fluid_pressure_kernels[cdim-1].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return ten_fluid_pressure_kernels[cdim-1].kernels[poly_order];
      break;
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static fluid_limiter_t
choose_fluid_limiter_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_fluid_limiter_kernels[cdim-1].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return ten_fluid_limiter_kernels[cdim-1].kernels[poly_order];
      break;
    default:
      assert(false);
      break;  
  }
}
