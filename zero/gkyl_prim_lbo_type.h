#pragma once

#include <gkyl_mat.h>
#include <gkyl_ref_count.h>

// Forward declare for use in function pointers
struct gkyl_prim_lbo_type;

// Self-primitive moment kernel pointer type
typedef void (*self_prim_t)(const struct gkyl_prim_lbo_type *prim, struct gkyl_mat *A,
  struct gkyl_mat *rhs, const double *moms, const double *boundary_corrections);

// Cross-primitive moment kernel pointer type
typedef void (*cross_prim_t)(const struct gkyl_prim_lbo_type *prim, struct gkyl_mat *A,
  struct gkyl_mat *rhs, const double betaGreenep1, const double m_self,
  const double *u_self, const double *vtsq_self, const double m_other,
  const double *u_other, const double *vtsq_other, const double *moms, const double *boundary_corrections);

struct gkyl_prim_lbo_type {
  int cdim; // config-space dim
  int pdim; // phase-space dim
  int poly_order; // polynomal order
  int num_config; // number of basis functions in config-space
  int num_phase; // number of basis functions in phase-space
  self_prim_t self_prim; // self-primitive moment calculation kernel
  cross_prim_t cross_prim; // cross-primitive moment calculation kernels
  struct gkyl_ref_count ref_count; // reference count

  uint32_t flag;
  struct gkyl_prim_lbo_type *on_dev; // pointer to itself or device data
};

/**
 * Acquire primitive moment object pointer. Delete using the release()
 * method
 *
 * @param prim Primitive moment object.
 */
struct gkyl_prim_lbo_type* gkyl_prim_lbo_type_acquire(const struct gkyl_prim_lbo_type* prim);

/**
 * Delete primitive moment object
 *
 * @param prim Primitive moment object to delete.
 */
void gkyl_prim_lbo_type_release(const struct gkyl_prim_lbo_type* prim);
