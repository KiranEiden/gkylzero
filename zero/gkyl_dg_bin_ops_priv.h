// Private header: not for direct use
#pragma once

#include <gkyl_binop_div_ser.h>
#include <gkyl_binop_mul_ser.h>
#include <gkyl_mat.h>
#include <gkyl_util.h>

// Function pointer type for multiplication
typedef void (*mul_op_t)(const double *f, const double *g, double *fg);
typedef struct gkyl_kern_op_count (*mul_op_count_t)(void);

// Function pointer type for setting matrices for division
typedef void (*div_set_op_t)(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *f, const double *g);

// Serendipity multiplication kernels
GKYL_CU_D
static struct { mul_op_t mul[4]; } ser_mul_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { binop_mul_1d_ser_p0, binop_mul_1d_ser_p1, binop_mul_1d_ser_p2, binop_mul_1d_ser_p3 },
  { binop_mul_2d_ser_p0, binop_mul_2d_ser_p1, binop_mul_2d_ser_p2, binop_mul_2d_ser_p3 },
  { binop_mul_3d_ser_p0, binop_mul_3d_ser_p1, binop_mul_3d_ser_p2, binop_mul_3d_ser_p3 }
};

static struct { mul_op_count_t mul[4]; } ser_mul_op_count_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { op_count_binop_mul_1d_ser_p0, op_count_binop_mul_1d_ser_p1, op_count_binop_mul_1d_ser_p2, op_count_binop_mul_1d_ser_p3 },
  { op_count_binop_mul_2d_ser_p0, op_count_binop_mul_2d_ser_p1, op_count_binop_mul_2d_ser_p2, op_count_binop_mul_2d_ser_p3 },
  { op_count_binop_mul_3d_ser_p0, op_count_binop_mul_3d_ser_p1, op_count_binop_mul_3d_ser_p2, op_count_binop_mul_3d_ser_p3 }
};

// Serendipity division kernels
GKYL_CU_D
static struct { div_set_op_t div[4]; } ser_div_set_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { binop_div_set_1d_ser_p0, binop_div_set_1d_ser_p1, binop_div_set_1d_ser_p2, binop_div_set_1d_ser_p3 },
  { binop_div_set_2d_ser_p0, binop_div_set_2d_ser_p1, binop_div_set_2d_ser_p2, binop_div_set_2d_ser_p3 },
  { binop_div_set_3d_ser_p0, binop_div_set_3d_ser_p1, binop_div_set_3d_ser_p2, binop_div_set_3d_ser_p3 } 
};

static mul_op_t
choose_ser_mul_kern(int dim, int poly_order)
{
  return ser_mul_list[dim].mul[poly_order];
}

static mul_op_count_t
choose_ser_mul_op_count_kern(int dim, int poly_order)
{
  return ser_mul_op_count_list[dim].mul[poly_order];
}

static div_set_op_t
choose_ser_div_set_kern(int dim, int poly_order)
{
  return ser_div_set_list[dim].div[poly_order];
}
