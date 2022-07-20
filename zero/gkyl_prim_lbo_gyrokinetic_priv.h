// Private header: not for direct use
#pragma once

#include <gkyl_prim_lbo_type.h>
#include <gkyl_prim_lbo_gyrokinetic_kernels.h>
#include <gkyl_mat.h>
#include <gkyl_util.h>

typedef void (*gyrokinetic_self_prim_t)(struct gkyl_mat *A, struct gkyl_mat *rhs, 
  const double *moms, const double *boundary_corrections);

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
static struct { int vdim[3]; } cv_index[] = {
  {-1, -1, -1}, // 0x makes no sense
  {-1,  0,  1}, // 1x kernel indices
  {-1, -1,  2}, // 2x kernel indices
  {-1, -1,  3}, // 3x kernel indices  
};

// for use in kernel tables
typedef struct { gyrokinetic_self_prim_t kernels[3]; } gkyl_prim_lbo_gyrokinetic_kern_list;

//
// Serendipity basis kernels
//

// self primitive moment kernel list
GKYL_CU_D
static const gkyl_prim_lbo_gyrokinetic_kern_list ser_self_prim_kernels[] = {
  // 1x kernels
  { NULL, gyrokinetic_self_prim_moments_1x1v_ser_p1, gyrokinetic_self_prim_moments_1x1v_ser_p2 }, // 0
  { NULL, gyrokinetic_self_prim_moments_1x2v_ser_p1, gyrokinetic_self_prim_moments_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, NULL, gyrokinetic_self_prim_moments_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, NULL, gyrokinetic_self_prim_moments_3x2v_ser_p2 }, // 3
};

struct prim_lbo_type_gyrokinetic {
  struct gkyl_prim_lbo_type prim; // Base object
  gyrokinetic_self_prim_t self_prim; // Self-primitive moments kernel
};

/**
 * Free primitive moment object.
 *
 * @param ref Reference counter for primitive moment to free
 */
void prim_lbo_gyrokinetic_free(const struct gkyl_ref_count *ref);

GKYL_CU_D
static void
self_prim(const struct gkyl_prim_lbo_type *prim, struct gkyl_mat *A, struct gkyl_mat *rhs, 
  const int* idx, const double *moms, const double *boundary_corrections)
{
  struct prim_lbo_type_gyrokinetic *prim_gyrokinetic = container_of(prim, struct prim_lbo_type_gyrokinetic, prim);

  return prim_gyrokinetic->self_prim(A, rhs, moms, boundary_corrections);
}
