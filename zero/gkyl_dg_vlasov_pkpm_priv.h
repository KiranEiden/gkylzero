#pragma once

// Private header, not for direct use in user code

#include <gkyl_array.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_kernels.h>
#include <gkyl_vlasov.h>

// Types for various kernels
typedef void (*vlasov_pkpm_stream_surf_t)(const double *w, const double *dxv,
  const double *bvarl, const double *bvarc, const double *bvarr, 
  const double *u_il, const double *u_ic, const double *u_ir, 
  const double *vth_sql, const double *vth_sqc, const double *vth_sqr,
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);

typedef void (*vlasov_pkpm_accel_surf_t)(const double *w, const double *dxv,
  const double *bb_grad_u, const double *p_force, const double *div_b, const double *p_perp_div_b, 
  const double *g_dist_sourcel, const double *g_dist_sourcec, const double *g_dist_sourcer, 
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

typedef void (*vlasov_pkpm_accel_boundary_surf_t)(const double *w, const double *dxv,
  const double *bb_grad_u, const double *p_force, const double *div_b, const double *p_perp_div_b, 
  const double *g_dist_sourceEdge, const double *g_dist_sourceSkin, 
  const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

// for use in kernel tables
typedef struct { vol_termf_t kernels[3]; } gkyl_dg_vlasov_pkpm_vol_kern_list;
typedef struct { vlasov_pkpm_stream_surf_t kernels[3]; } gkyl_dg_vlasov_pkpm_stream_surf_kern_list;
typedef struct { vlasov_pkpm_accel_surf_t kernels[3]; } gkyl_dg_vlasov_pkpm_accel_surf_kern_list;
typedef struct { vlasov_pkpm_accel_boundary_surf_t kernels[3]; } gkyl_dg_vlasov_pkpm_accel_boundary_surf_kern_list;

struct dg_vlasov_pkpm {
  struct gkyl_dg_eqn eqn; // Base object
  int cdim; // Config-space dimensions
  int pdim; // Phase-space dimensions
  vlasov_pkpm_stream_surf_t stream_surf[3]; // Surface terms for streaming
  vlasov_pkpm_accel_surf_t accel_surf; // Surface terms for acceleration
  vlasov_pkpm_accel_boundary_surf_t accel_boundary_surf; // Surface terms for acceleration
  struct gkyl_range conf_range; // configuration space range
  struct gkyl_range phase_range; // phase space range
  struct gkyl_dg_vlasov_pkpm_auxfields auxfields; // Auxiliary fields.
};

//
// Serendipity volume kernels
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static double
kernel_vlasov_pkpm_vol_1x1v_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_pkpm *vlasov_pkpm = container_of(eqn, struct dg_vlasov_pkpm, eqn);

  long cidx = gkyl_range_idx(&vlasov_pkpm->conf_range, idx);
  long pidx = gkyl_range_idx(&vlasov_pkpm->phase_range, idx);
  return vlasov_pkpm_vol_1x1v_ser_p1(xc, dx, 
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.bvar, cidx),
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.u_i, cidx),
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.bb_grad_u, cidx),
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.p_force, cidx),
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.div_b, cidx),
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.p_perp_source, cidx),
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.p_perp_div_b, cidx),
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.g_dist_source, pidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_pkpm_vol_1x1v_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_pkpm *vlasov_pkpm = container_of(eqn, struct dg_vlasov_pkpm, eqn);

  long cidx = gkyl_range_idx(&vlasov_pkpm->conf_range, idx);
  long pidx = gkyl_range_idx(&vlasov_pkpm->phase_range, idx);
  return vlasov_pkpm_vol_1x1v_ser_p2(xc, dx, 
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.bvar, cidx),
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.u_i, cidx),
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.bb_grad_u, cidx),
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.p_force, cidx),
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.div_b, cidx),
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.p_perp_source, cidx),
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.p_perp_div_b, cidx),
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.g_dist_source, pidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_pkpm_vol_2x1v_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_pkpm *vlasov_pkpm = container_of(eqn, struct dg_vlasov_pkpm, eqn);

  long cidx = gkyl_range_idx(&vlasov_pkpm->conf_range, idx);
  long pidx = gkyl_range_idx(&vlasov_pkpm->phase_range, idx);
  return vlasov_pkpm_vol_2x1v_ser_p1(xc, dx, 
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.bvar, cidx),
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.u_i, cidx),
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.bb_grad_u, cidx),
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.p_force, cidx),
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.div_b, cidx),
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.p_perp_source, cidx),
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.p_perp_div_b, cidx),
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.g_dist_source, pidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_pkpm_vol_2x1v_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_pkpm *vlasov_pkpm = container_of(eqn, struct dg_vlasov_pkpm, eqn);

  long cidx = gkyl_range_idx(&vlasov_pkpm->conf_range, idx);
  long pidx = gkyl_range_idx(&vlasov_pkpm->phase_range, idx);
  return vlasov_pkpm_vol_2x1v_ser_p2(xc, dx, 
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.bvar, cidx),
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.u_i, cidx),
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.bb_grad_u, cidx),
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.p_force, cidx),
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.div_b, cidx),
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.p_perp_source, cidx),
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.p_perp_div_b, cidx),
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.g_dist_source, pidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_pkpm_vol_3x1v_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_pkpm *vlasov_pkpm = container_of(eqn, struct dg_vlasov_pkpm, eqn);

  long cidx = gkyl_range_idx(&vlasov_pkpm->conf_range, idx);
  long pidx = gkyl_range_idx(&vlasov_pkpm->phase_range, idx);
  return vlasov_pkpm_vol_3x1v_ser_p1(xc, dx, 
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.bvar, cidx),
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.u_i, cidx),
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.bb_grad_u, cidx),
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.p_force, cidx),
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.div_b, cidx),
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.p_perp_source, cidx),
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.p_perp_div_b, cidx),
    (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.g_dist_source, pidx),
    qIn, qRhsOut);
}

// Volume kernel list
GKYL_CU_D
static const gkyl_dg_vlasov_pkpm_vol_kern_list ser_vol_kernels[] = {
  // 1x kernels
  { NULL, kernel_vlasov_pkpm_vol_1x1v_ser_p1, kernel_vlasov_pkpm_vol_1x1v_ser_p2 }, // 0
  // 2x kernels
  { NULL, kernel_vlasov_pkpm_vol_2x1v_ser_p1, kernel_vlasov_pkpm_vol_2x1v_ser_p2 }, // 1
  // 3x kernels
  { NULL, kernel_vlasov_pkpm_vol_3x1v_ser_p1, NULL }, // 2
};

// Streaming surface kernel list: x-direction
GKYL_CU_D
static const gkyl_dg_vlasov_pkpm_stream_surf_kern_list ser_stream_surf_x_kernels[] = {
  // 1x kernels
  { NULL, vlasov_pkpm_surfx_1x1v_ser_p1, vlasov_pkpm_surfx_1x1v_ser_p2 }, // 0
  // 2x kernels
  { NULL, vlasov_pkpm_surfx_2x1v_ser_p1, vlasov_pkpm_surfx_2x1v_ser_p2 }, // 1
  // 3x kernels
  { NULL, vlasov_pkpm_surfx_3x1v_ser_p1, NULL }, // 2
};

// Streaming surface kernel list: y-direction
GKYL_CU_D
static const gkyl_dg_vlasov_pkpm_stream_surf_kern_list ser_stream_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, vlasov_pkpm_surfy_2x1v_ser_p1, vlasov_pkpm_surfy_2x1v_ser_p2 }, // 1
  // 3x kernels
  { NULL, vlasov_pkpm_surfy_3x1v_ser_p1, NULL }, // 2
};

// Streaming surface kernel list: z-direction
GKYL_CU_D
static const gkyl_dg_vlasov_pkpm_stream_surf_kern_list ser_stream_surf_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, NULL, NULL }, // 1
  // 3x kernels
  { NULL, vlasov_pkpm_surfz_3x1v_ser_p1, NULL }, // 2
};

// Acceleration surface kernel list: vpar-direction
GKYL_CU_D
static const gkyl_dg_vlasov_pkpm_accel_surf_kern_list ser_accel_surf_vpar_kernels[] = {
  // 1x kernels
  { NULL, vlasov_pkpm_surfvpar_1x1v_ser_p1, vlasov_pkpm_surfvpar_1x1v_ser_p2 }, // 0
  // 2x kernels
  { NULL, vlasov_pkpm_surfvpar_2x1v_ser_p1, vlasov_pkpm_surfvpar_2x1v_ser_p2 }, // 1
  // 3x kernels
  { NULL, vlasov_pkpm_surfvpar_3x1v_ser_p1, NULL }, // 2
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vpar-direction
GKYL_CU_D
static const gkyl_dg_vlasov_pkpm_accel_boundary_surf_kern_list ser_accel_boundary_surf_vpar_kernels[] = {
  // 1x kernels
  { NULL, vlasov_pkpm_boundary_surfvpar_1x1v_ser_p1, vlasov_pkpm_boundary_surfvpar_1x1v_ser_p2 }, // 0
  // 2x kernels
  { NULL, vlasov_pkpm_boundary_surfvpar_2x1v_ser_p1, vlasov_pkpm_boundary_surfvpar_2x1v_ser_p2 }, // 1
  // 3x kernels
  { NULL, vlasov_pkpm_boundary_surfvpar_3x1v_ser_p1, NULL }, // 2
};

/**
 * Free vlasov_pkpm eqn object.
 *
 * @param ref Reference counter for vlasov_pkpm eqn
 */
void gkyl_vlasov_pkpm_free(const struct gkyl_ref_count *ref);

GKYL_CU_D
static void
surf(const struct gkyl_dg_eqn *eqn, 
  int dir,
  const double*  xcL, const double*  xcC, const double*  xcR, 
  const double*  dxL, const double* dxC, const double* dxR,
  const int*  idxL, const int*  idxC, const int*  idxR,
  const double* qInL, const double*  qInC, const double*  qInR, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_pkpm *vlasov_pkpm = container_of(eqn, struct dg_vlasov_pkpm, eqn);
  long cidx = gkyl_range_idx(&vlasov_pkpm->conf_range, idxC);
  long pidx = gkyl_range_idx(&vlasov_pkpm->phase_range, idxC);
  if (dir < vlasov_pkpm->cdim) {
    long cidx_l = gkyl_range_idx(&vlasov_pkpm->conf_range, idxL);
    long cidx_r = gkyl_range_idx(&vlasov_pkpm->conf_range, idxR);
    vlasov_pkpm->stream_surf[dir]
      (xcC, dxC, 
      (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.bvar, cidx_l), 
      (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.bvar, cidx), 
      (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.bvar, cidx_r), 
      (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.u_i, cidx_l), 
      (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.u_i, cidx), 
      (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.u_i, cidx_r), 
      (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.vth_sq, cidx_l), 
      (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.vth_sq, cidx), 
      (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.vth_sq, cidx_r), 
      qInL, qInC, qInR, qRhsOut);
  }
  else {
    long pidx_l = gkyl_range_idx(&vlasov_pkpm->phase_range, idxL);
    long pidx_r = gkyl_range_idx(&vlasov_pkpm->phase_range, idxR);
    vlasov_pkpm->accel_surf(xcC, dxC,
      (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.bb_grad_u, cidx),
      (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.p_force, cidx),
      (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.div_b, cidx),
      (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.p_perp_div_b, cidx),
      (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.g_dist_source, pidx_l),
      (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.g_dist_source, pidx),
      (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.g_dist_source, pidx_r),
      qInL, qInC, qInR, qRhsOut);
  }
}

GKYL_CU_D
static void
boundary_surf(const struct gkyl_dg_eqn *eqn,
  int dir,
  const double*  xcEdge, const double*  xcSkin,
  const double*  dxEdge, const double* dxSkin,
  const int* idxEdge, const int* idxSkin, const int edge,
  const double* qInEdge, const double* qInSkin, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_pkpm *vlasov_pkpm = container_of(eqn, struct dg_vlasov_pkpm, eqn);

  // only in vpar direction
  if (dir == vlasov_pkpm->cdim) {
    long cidx = gkyl_range_idx(&vlasov_pkpm->conf_range, idxSkin);
    long pidxSkin = gkyl_range_idx(&vlasov_pkpm->phase_range, idxSkin);
    long pidxEdge = gkyl_range_idx(&vlasov_pkpm->phase_range, idxEdge);
    vlasov_pkpm->accel_boundary_surf(xcSkin, dxSkin,
      (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.bb_grad_u, cidx),
      (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.p_force, cidx),
      (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.div_b, cidx),
      (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.p_perp_div_b, cidx),
      (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.g_dist_source, pidxEdge),
      (const double*) gkyl_array_cfetch(vlasov_pkpm->auxfields.g_dist_source, pidxSkin),
      edge, qInEdge, qInSkin, qRhsOut);
  }
}
