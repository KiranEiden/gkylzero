#pragma once

// Private header, not for direct use in user code

#include <math.h>
#include <gkyl_array.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

#define GRM_D1 0
#define GRM_D2 1
#define GRM_D3 2
#define GRM_B1 3
#define GRM_B2 4
#define GRM_B3 5
#define GRM_PHI 6
#define GRM_PSI 7
#define GRM_LAPSE 8
#define GRM_S1 9
#define GRM_S2 10
#define GRM_S3 11
#define GRM_DET 12
#define GRM_GAM11 13
#define GRM_GAM12 14
#define GRM_GAM13 15
#define GRM_GAM21 16
#define GRM_GAM22 17
#define GRM_GAM23 18
#define GRM_GAM31 19
#define GRM_GAM32 20
#define GRM_GAM33 21
#define GRM_GINV11 22
#define GRM_GINV12 23
#define GRM_GINV13 24
#define GRM_GINV21 25
#define GRM_GINV22 26
#define GRM_GINV23 27
#define GRM_GINV31 28
#define GRM_GINV32 29
#define GRM_GINV33 30
#define GRM_IER 31

#define GRM_NEQ 32

// Note: Can likely solve for both the covariant and contravariant components.
// This would get rid of the metric components in the state vector, but requires coding up additional fluxes.

struct wv_gr_maxwell {
  struct gkyl_wv_eqn eqn; // base object
  double c; // light speed
  double e_fact, b_fact; // electric and magnetic correction factors
};

/**
 * Free Maxwell eqn object.
 *
 * @param ref Reference counter for Maxwell eqn
 */
void gkyl_wv_gr_maxwell_free(const struct gkyl_ref_count *ref);

static inline void
gkyl_gr_maxwell_rm_rt_det(const double *vin, double *vout, const double rt_gam_det)
{
  vout[0] = vin[0] / rt_gam_det;
  vout[1] = vin[1] / rt_gam_det;
  vout[2] = vin[2] / rt_gam_det;
}

static inline void
gkyl_gr_maxwell_lower_ind(const double *vin, double *vout, const double *gam)
{
  vout[0] = gam[0]*vin[0] + gam[1]*vin[1] + gam[2]*vin[2];
  vout[1] = gam[3]*vin[0] + gam[4]*vin[1] + gam[5]*vin[2];
  vout[2] = gam[6]*vin[0] + gam[7]*vin[1] + gam[8]*vin[2];
}

/**
 * Compute maximum absolute speed.
 * 
 * @param c Speed of light
 * @param e_fact Correction speed for div(E) correction
 * @param b_fact Correction speed for div(B) correction
 * @param q Conserved variables
 * @return Maximum absolute speed for given q
 */
static inline double
gkyl_gr_maxwell_max_abs_speed(double c, double e_fact, double b_fact, const double q[GRM_NEQ])
{
  if (q[GRM_IER] < 0.0)
  {
    return 0.0;
  }
  
  // Speed for case of Euclidean metric + gauge variables
  //double spd = 0.0;
  //double spd_d;
  //for(int d = 0; d < 3; d++) {
    //spd_d = fmax(fabs(q[6] - q[7+d]), fabs(q[6] + q[7+d]));
    //spd = spd < spd_d ? spd_d : spd;
  //}
  
  return c/fabs(q[GRM_LAPSE]*q[GRM_DET]);
}

/**
 * Compute flux.
 * 
 * @param c Speed of light
 * @param e_fact Correction speed for div(E) correction
 * @param b_fact Correction speed for div(B) correction
 * @param Conserved variables
 * @param flux On output, the flux in direction 'dir'
 */
static inline void
gkyl_gr_maxwell_flux(double c, double e_fact, double b_fact, const double q[GRM_NEQ], double flux[GRM_NEQ])
{
  if (q[GRM_IER] < 0.0)
  {
    // Fluxes are 0 in excision region
    for (int i = 0; i < GRM_NEQ; ++i) flux[i] = 0.0;
    return;
  }
  
  const double c2 = c*c;
  
  double D[3];
  gkyl_gr_maxwell_rm_rt_det(&q[GRM_D1], D, q[GRM_DET]);
  double B[3];
  gkyl_gr_maxwell_rm_rt_det(&q[GRM_B1], B, q[GRM_DET]);

  //flux[0] = e_fact*c2*q[6]; // e_fact*c^2*phi
  flux[GRM_D1] = 0.0;
  flux[GRM_D2] = c2*B[2]; // c^2*Bz
  flux[GRM_D3] = -c2*B[1]; // -c^2*By
  //flux[3] = b_fact*q[7]; // b_fact*psi
  flux[GRM_B1] = 0.0;
  flux[GRM_B2] = -D[2]; // -Ez
  flux[GRM_B3] = D[1]; // Ey
  //flux[6] = e_fact*q[0]; // e_fact*Ex
  //flux[7] = b_fact*c2*q[3]; // b_fact*c^2*Bx
  
  for (int i = GRM_LAPSE; i < GRM_NEQ; ++i) flux[i] = 0.0;
}

static void
gkyl_gr_maxwell_extrapolate_flux(double c, double e_fact, double b_fact,
    const double q[GRM_NEQ], const double flux_in[GRM_NEQ], double flux_extrap[GRM_NEQ])
{
  if (q[GRM_IER] < 0.0)
  {
    // Just copy fluxes if in excision region
    for (int i = 0; i < GRM_NEQ; ++i) flux_extrap[i] = flux_in[i];
    return;
  }
  
  const double c2 = c*c;
  const double lapse = q[GRM_LAPSE];
  
  double D[3], D_cov[3];
  double B[3], B_cov[3];
  double shift_cov[3];
  
  D[0] = q[GRM_D1] / q[GRM_DET];
  D[1] = flux_in[5];
  D[2] = -flux_in[4];
  
  B[0] = q[GRM_B1] / q[GRM_DET];
  B[1] = -flux_in[2] / c2;
  B[2] = flux_in[1] / c2;
  
  gkyl_gr_maxwell_lower_ind(D, D_cov, &q[GRM_GAM11]);
  gkyl_gr_maxwell_lower_ind(B, B_cov, &q[GRM_GAM11]);
  gkyl_gr_maxwell_lower_ind(&q[GRM_S1], shift_cov, &q[GRM_GAM11]);
  
  flux_extrap[GRM_D1] = flux_in[0]; // 0 for Maxwell curl equations
  flux_extrap[GRM_D2] = c2 * (lapse*B_cov[2] + (shift_cov[1]*D_cov[0] - shift_cov[0]*D_cov[1]) / q[GRM_DET]);
  flux_extrap[GRM_D3] = -c2 * (lapse*B_cov[1] + (shift_cov[0]*D_cov[2] - shift_cov[2]*D_cov[0]) / q[GRM_DET]);
  flux_extrap[GRM_B1] = flux_in[3]; // 0 for Maxwell curl equations
  flux_extrap[GRM_B2] = -(lapse*D_cov[2] + (shift_cov[0]*B_cov[1] - shift_cov[1]*B_cov[0]) / q[GRM_DET]);
  flux_extrap[GRM_B3] = lapse*D_cov[1] + (shift_cov[2]*B_cov[0] - shift_cov[0]*B_cov[2] / q[GRM_DET]);
  
  // Include cell volume correction here (can be handled in geometry code)?
  for (int i = 0; i < GRM_LAPSE; ++i) {
    flux_extrap[i] /= q[GRM_DET];
  }
  
  for (int i = GRM_LAPSE; i < GRM_NEQ; ++i) flux_extrap[i] = 0.0;
}

static inline void
cons_to_riem(const struct gkyl_wv_eqn *eqn,
  const double *qstate, const double *qin, double *wout)
{
  // TODO: this should use proper L matrix
  for (int i=0; i<GRM_NEQ; ++i)
    wout[i] = qin[i];
}

static inline void
riem_to_cons(const struct gkyl_wv_eqn *eqn,
  const double *qstate, const double *win, double *qout)
{
  // TODO: this should use proper L matrix
  for (int i=0; i<GRM_NEQ; ++i)
    qout[i] = win[i];
}

static void
gr_maxwell_wall(double t, int nc, const double *skin, double * GKYL_RESTRICT ghost, void *ctx)
{
  // zero-tangent for E field
  ghost[GRM_D1] = skin[GRM_D1];
  ghost[GRM_D2] = -skin[GRM_D2];
  ghost[GRM_D3] = -skin[GRM_D3];

  // zero-normal for B field
  ghost[GRM_B1] = -skin[GRM_B1];
  ghost[GRM_B2] = skin[GRM_B2];
  ghost[GRM_B3] = skin[GRM_B3];
  
  for (int i = GRM_LAPSE; i < GRM_NEQ; ++i) ghost[i] = skin[i];

  // correction potential
  //ghost[6] = -skin[6];
  //ghost[7] = skin[7];
}

static inline void
rot_to_local(const double *tau1, const double *tau2, const double *norm,
  const double *GKYL_RESTRICT qglobal, double *GKYL_RESTRICT qlocal)
{
  // Rotate D to local coordinates
  qlocal[GRM_D1] = (qglobal[GRM_D1] * norm[0]) + (qglobal[GRM_D2] * norm[1]) + (qglobal[GRM_D3] * norm[2]);
  qlocal[GRM_D2] = (qglobal[GRM_D1] * tau1[0]) + (qglobal[GRM_D2] * tau1[1]) + (qglobal[GRM_D3] * tau1[2]);
  qlocal[GRM_D3] = (qglobal[GRM_D1] * tau2[0]) + (qglobal[GRM_D2] * tau2[1]) + (qglobal[GRM_D3] * tau2[2]);
  
  // Rotate B to local coordinates
  qlocal[GRM_B1] = (qglobal[GRM_B1] * norm[0]) + (qglobal[GRM_B2] * norm[1]) + (qglobal[GRM_B3] * norm[2]);
  qlocal[GRM_B2] = (qglobal[GRM_B1] * tau1[0]) + (qglobal[GRM_B2] * tau1[1]) + (qglobal[GRM_B3] * tau1[2]);
  qlocal[GRM_B3] = (qglobal[GRM_B1] * tau2[0]) + (qglobal[GRM_B2] * tau2[1]) + (qglobal[GRM_B3] * tau2[2]);
  
  // Copy scalar lapse function
  qlocal[GRM_LAPSE] = qglobal[GRM_LAPSE];
  
  // Rotate shift vector (beta) to local coordinates
  qlocal[GRM_S1] = (qglobal[GRM_S1] * norm[0]) + (qglobal[GRM_S2] * norm[1]) + (qglobal[GRM_S3] * norm[2]);
  qlocal[GRM_S2] = (qglobal[GRM_S1] * tau1[0]) + (qglobal[GRM_S2] * tau1[1]) + (qglobal[GRM_S3] * tau1[2]);
  qlocal[GRM_S3] = (qglobal[GRM_S1] * tau2[0]) + (qglobal[GRM_S2] * tau2[1]) + (qglobal[GRM_S3] * tau2[2]);
  
  // Copy the sqrt of the spatial metric determinant
  qlocal[GRM_DET] = qglobal[GRM_DET];
  
  // Transform the spatial metric (stored in row-major order)
  double x[9]; // temporary array
  x[0] = (qglobal[GRM_GAM11] * norm[0]) + (qglobal[GRM_GAM12] * norm[1]) + (qglobal[GRM_GAM13] * norm[2]);
  x[1] = (qglobal[GRM_GAM11] * tau1[0]) + (qglobal[GRM_GAM12] * tau1[1]) + (qglobal[GRM_GAM13] * tau1[2]);
  x[2] = (qglobal[GRM_GAM11] * tau2[0]) + (qglobal[GRM_GAM12] * tau2[1]) + (qglobal[GRM_GAM13] * tau2[2]);
  x[3] = (qglobal[GRM_GAM21] * norm[0]) + (qglobal[GRM_GAM22] * norm[1]) + (qglobal[GRM_GAM23] * norm[2]);
  x[4] = (qglobal[GRM_GAM21] * tau1[0]) + (qglobal[GRM_GAM22] * tau1[1]) + (qglobal[GRM_GAM23] * tau1[2]);
  x[5] = (qglobal[GRM_GAM21] * tau2[0]) + (qglobal[GRM_GAM22] * tau2[1]) + (qglobal[GRM_GAM23] * tau2[2]);
  x[6] = (qglobal[GRM_GAM31] * norm[0]) + (qglobal[GRM_GAM32] * norm[1]) + (qglobal[GRM_GAM33] * norm[2]);
  x[7] = (qglobal[GRM_GAM31] * tau1[0]) + (qglobal[GRM_GAM32] * tau1[1]) + (qglobal[GRM_GAM33] * tau1[2]);
  x[8] = (qglobal[GRM_GAM31] * tau2[0]) + (qglobal[GRM_GAM32] * tau2[1]) + (qglobal[GRM_GAM33] * tau2[2]);
  
  qlocal[GRM_GAM11] = (x[0] * norm[0]) + (x[3] * norm[1]) + (x[6] * norm[2]);
  qlocal[GRM_GAM12] = (x[0] * tau1[0]) + (x[3] * tau1[1]) + (x[6] * tau1[2]);
  qlocal[GRM_GAM13] = (x[0] * tau2[0]) + (x[3] * tau2[1]) + (x[6] * tau2[2]);
  qlocal[GRM_GAM21] = (x[1] * norm[0]) + (x[4] * norm[1]) + (x[7] * norm[2]);
  qlocal[GRM_GAM22] = (x[1] * tau1[0]) + (x[4] * tau1[1]) + (x[7] * tau1[2]);
  qlocal[GRM_GAM23] = (x[1] * tau2[0]) + (x[4] * tau2[1]) + (x[7] * tau2[2]);
  qlocal[GRM_GAM31] = (x[2] * norm[0]) + (x[5] * norm[1]) + (x[8] * norm[2]);
  qlocal[GRM_GAM32] = (x[2] * tau1[0]) + (x[5] * tau1[1]) + (x[8] * tau1[2]);
  qlocal[GRM_GAM33] = (x[2] * tau2[0]) + (x[5] * tau2[1]) + (x[8] * tau2[2]);
  
  // Transform the spatial metric inverse (stored in row-major order)
  x[0] = (qglobal[GRM_GINV11] * norm[0]) + (qglobal[GRM_GINV12] * norm[1]) + (qglobal[GRM_GINV13] * norm[2]);
  x[1] = (qglobal[GRM_GINV11] * tau1[0]) + (qglobal[GRM_GINV12] * tau1[1]) + (qglobal[GRM_GINV13] * tau1[2]);
  x[2] = (qglobal[GRM_GINV11] * tau2[0]) + (qglobal[GRM_GINV12] * tau2[1]) + (qglobal[GRM_GINV13] * tau2[2]);
  x[3] = (qglobal[GRM_GINV21] * norm[0]) + (qglobal[GRM_GINV22] * norm[1]) + (qglobal[GRM_GINV23] * norm[2]);
  x[4] = (qglobal[GRM_GINV21] * tau1[0]) + (qglobal[GRM_GINV22] * tau1[1]) + (qglobal[GRM_GINV23] * tau1[2]);
  x[5] = (qglobal[GRM_GINV21] * tau2[0]) + (qglobal[GRM_GINV22] * tau2[1]) + (qglobal[GRM_GINV23] * tau2[2]);
  x[6] = (qglobal[GRM_GINV31] * norm[0]) + (qglobal[GRM_GINV32] * norm[1]) + (qglobal[GRM_GINV33] * norm[2]);
  x[7] = (qglobal[GRM_GINV31] * tau1[0]) + (qglobal[GRM_GINV32] * tau1[1]) + (qglobal[GRM_GINV33] * tau1[2]);
  x[8] = (qglobal[GRM_GINV31] * tau2[0]) + (qglobal[GRM_GINV32] * tau2[1]) + (qglobal[GRM_GINV33] * tau2[2]);
  
  qlocal[GRM_GINV11] = (x[0] * norm[0]) + (x[3] * norm[1]) + (x[6] * norm[2]);
  qlocal[GRM_GINV12] = (x[0] * tau1[0]) + (x[3] * tau1[1]) + (x[6] * tau1[2]);
  qlocal[GRM_GINV13] = (x[0] * tau2[0]) + (x[3] * tau2[1]) + (x[6] * tau2[2]);
  qlocal[GRM_GINV21] = (x[1] * norm[0]) + (x[4] * norm[1]) + (x[7] * norm[2]);
  qlocal[GRM_GINV22] = (x[1] * tau1[0]) + (x[4] * tau1[1]) + (x[7] * tau1[2]);
  qlocal[GRM_GINV23] = (x[1] * tau2[0]) + (x[4] * tau2[1]) + (x[7] * tau2[2]);
  qlocal[GRM_GINV31] = (x[2] * norm[0]) + (x[5] * norm[1]) + (x[8] * norm[2]);
  qlocal[GRM_GINV32] = (x[2] * tau1[0]) + (x[5] * tau1[1]) + (x[8] * tau1[2]);
  qlocal[GRM_GINV33] = (x[2] * tau2[0]) + (x[5] * tau2[1]) + (x[8] * tau2[2]);
  
  // Copy in_excision_region flag
  qlocal[GRM_IER] = qglobal[GRM_IER];
  
  // Correction potentials are scalars and unchanged
  //qlocal[6] = qglobal[6];
  //qlocal[7] = qglobal[7];
}

static inline void
rot_to_global(const double *tau1, const double *tau2, const double *norm,
  const double *GKYL_RESTRICT qlocal, double *GKYL_RESTRICT qglobal)
{
  // Rotate D back to global coordinates
  qglobal[GRM_D1] = (qlocal[GRM_D1] * norm[0]) + (qlocal[GRM_D2] * tau1[0]) + (qlocal[GRM_D3] * tau2[0]);
  qglobal[GRM_D2] = (qlocal[GRM_D1] * norm[1]) + (qlocal[GRM_D2] * tau1[1]) + (qlocal[GRM_D3] * tau2[1]);
  qglobal[GRM_D3] = (qlocal[GRM_D1] * norm[2]) + (qlocal[GRM_D2] * tau1[2]) + (qlocal[GRM_D3] * tau2[2]);
  
  // Rotate B back to global coordinates
  qglobal[GRM_B1] = (qlocal[GRM_B1] * norm[0]) + (qlocal[GRM_B2] * tau1[0]) + (qlocal[GRM_B3] * tau2[0]);
  qglobal[GRM_B2] = (qlocal[GRM_B1] * norm[1]) + (qlocal[GRM_B2] * tau1[1]) + (qlocal[GRM_B3] * tau2[1]);
  qglobal[GRM_B3] = (qlocal[GRM_B1] * norm[2]) + (qlocal[GRM_B2] * tau1[2]) + (qlocal[GRM_B3] * tau2[2]);
  
  // Again just copy scalar lapse function
  qglobal[GRM_LAPSE] = qlocal[GRM_LAPSE];
  
  // Rotate shift vector (beta) to back to global coordinates
  qglobal[GRM_S1] = (qlocal[GRM_S1] * norm[0]) + (qlocal[GRM_S2] * tau1[0]) + (qlocal[GRM_S3] * tau2[0]);
  qglobal[GRM_S2] = (qlocal[GRM_S1] * norm[1]) + (qlocal[GRM_S2] * tau1[1]) + (qlocal[GRM_S3] * tau2[1]);
  qglobal[GRM_S3] = (qlocal[GRM_S1] * norm[2]) + (qlocal[GRM_S2] * tau1[2]) + (qlocal[GRM_S3] * tau2[2]);
  
  // Copy sqrt of spatial metric determinant
  qglobal[GRM_DET] = qlocal[GRM_DET];
  
  // Transform the spatial metric
  double x[9]; // again a temporary array
  x[0] = (qlocal[GRM_GAM11] * norm[0]) + (qlocal[GRM_GAM12] * tau1[0]) + (qlocal[GRM_GAM13] * tau2[0]);
  x[1] = (qlocal[GRM_GAM11] * norm[1]) + (qlocal[GRM_GAM12] * tau1[1]) + (qlocal[GRM_GAM13] * tau2[1]);
  x[2] = (qlocal[GRM_GAM11] * norm[2]) + (qlocal[GRM_GAM12] * tau1[2]) + (qlocal[GRM_GAM13] * tau2[2]);
  x[3] = (qlocal[GRM_GAM21] * norm[0]) + (qlocal[GRM_GAM22] * tau1[0]) + (qlocal[GRM_GAM23] * tau2[0]);
  x[4] = (qlocal[GRM_GAM21] * norm[1]) + (qlocal[GRM_GAM22] * tau1[1]) + (qlocal[GRM_GAM23] * tau2[1]);
  x[5] = (qlocal[GRM_GAM21] * norm[2]) + (qlocal[GRM_GAM22] * tau1[2]) + (qlocal[GRM_GAM23] * tau2[2]);
  x[6] = (qlocal[GRM_GAM31] * norm[0]) + (qlocal[GRM_GAM32] * tau1[0]) + (qlocal[GRM_GAM33] * tau2[0]);
  x[7] = (qlocal[GRM_GAM31] * norm[1]) + (qlocal[GRM_GAM32] * tau1[1]) + (qlocal[GRM_GAM33] * tau2[1]);
  x[8] = (qlocal[GRM_GAM31] * norm[2]) + (qlocal[GRM_GAM32] * tau1[2]) + (qlocal[GRM_GAM33] * tau2[2]);
  
  qglobal[GRM_GAM11] = (x[0] * norm[0]) + (x[3] * tau1[0]) + (x[6] * tau2[0]);
  qglobal[GRM_GAM12] = (x[0] * norm[1]) + (x[3] * tau1[1]) + (x[6] * tau2[1]);
  qglobal[GRM_GAM13] = (x[0] * norm[2]) + (x[3] * tau1[2]) + (x[6] * tau2[2]);
  qglobal[GRM_GAM21] = (x[1] * norm[0]) + (x[4] * tau1[0]) + (x[7] * tau2[0]);
  qglobal[GRM_GAM22] = (x[1] * norm[1]) + (x[4] * tau1[1]) + (x[7] * tau2[1]);
  qglobal[GRM_GAM23] = (x[1] * norm[2]) + (x[4] * tau1[2]) + (x[7] * tau2[2]);
  qglobal[GRM_GAM31] = (x[2] * norm[0]) + (x[5] * tau1[0]) + (x[8] * tau2[0]);
  qglobal[GRM_GAM32] = (x[2] * norm[1]) + (x[5] * tau1[1]) + (x[8] * tau2[1]);
  qglobal[GRM_GAM33] = (x[2] * norm[2]) + (x[5] * tau1[2]) + (x[8] * tau2[2]);
  
  // Transform the spatial metric inverse
  x[0] = (qlocal[GRM_GINV11] * norm[0]) + (qlocal[GRM_GINV12] * tau1[0]) + (qlocal[GRM_GINV13] * tau2[0]);
  x[1] = (qlocal[GRM_GINV11] * norm[1]) + (qlocal[GRM_GINV12] * tau1[1]) + (qlocal[GRM_GINV13] * tau2[1]);
  x[2] = (qlocal[GRM_GINV11] * norm[2]) + (qlocal[GRM_GINV12] * tau1[2]) + (qlocal[GRM_GINV13] * tau2[2]);
  x[3] = (qlocal[GRM_GINV21] * norm[0]) + (qlocal[GRM_GINV22] * tau1[0]) + (qlocal[GRM_GINV23] * tau2[0]);
  x[4] = (qlocal[GRM_GINV21] * norm[1]) + (qlocal[GRM_GINV22] * tau1[1]) + (qlocal[GRM_GINV23] * tau2[1]);
  x[5] = (qlocal[GRM_GINV21] * norm[2]) + (qlocal[GRM_GINV22] * tau1[2]) + (qlocal[GRM_GINV23] * tau2[2]);
  x[6] = (qlocal[GRM_GINV31] * norm[0]) + (qlocal[GRM_GINV32] * tau1[0]) + (qlocal[GRM_GINV33] * tau2[0]);
  x[7] = (qlocal[GRM_GINV31] * norm[1]) + (qlocal[GRM_GINV32] * tau1[1]) + (qlocal[GRM_GINV33] * tau2[1]);
  x[8] = (qlocal[GRM_GINV31] * norm[2]) + (qlocal[GRM_GINV32] * tau1[2]) + (qlocal[GRM_GINV33] * tau2[2]);
  
  qglobal[GRM_GINV11] = (x[0] * norm[0]) + (x[3] * tau1[0]) + (x[6] * tau2[0]);
  qglobal[GRM_GINV12] = (x[0] * norm[1]) + (x[3] * tau1[1]) + (x[6] * tau2[1]);
  qglobal[GRM_GINV13] = (x[0] * norm[2]) + (x[3] * tau1[2]) + (x[6] * tau2[2]);
  qglobal[GRM_GINV21] = (x[1] * norm[0]) + (x[4] * tau1[0]) + (x[7] * tau2[0]);
  qglobal[GRM_GINV22] = (x[1] * norm[1]) + (x[4] * tau1[1]) + (x[7] * tau2[1]);
  qglobal[GRM_GINV23] = (x[1] * norm[2]) + (x[4] * tau1[2]) + (x[7] * tau2[2]);
  qglobal[GRM_GINV31] = (x[2] * norm[0]) + (x[5] * tau1[0]) + (x[8] * tau2[0]);
  qglobal[GRM_GINV32] = (x[2] * norm[1]) + (x[5] * tau1[1]) + (x[8] * tau2[1]);
  qglobal[GRM_GINV33] = (x[2] * norm[2]) + (x[5] * tau1[2]) + (x[8] * tau2[2]);
  
  // Copy the in_excision_region flag
  qglobal[GRM_IER] = qlocal[GRM_IER];
  
  // Correction potentials are scalars and unchanged
  //qglobal[6] = qlocal[6];
  //qglobal[7] = qlocal[7];
}

// Flux difference for special relativistic Maxwell
static double
flux_jump(const struct gkyl_wv_eqn *eqn, const double *ql, const double *qr, double *flux_jump)
{
  bool in_excision_region = (ql[GRM_IER] < 0.0) || (qr[GRM_IER] < 0.0);
  if (in_excision_region) {
    for (int m=0; m<GRM_NEQ; ++m) flux_jump[m] = 0.0;
    return 0.0;
  }
  
  const struct wv_gr_maxwell *maxwell = container_of(eqn, struct wv_gr_maxwell, eqn);

  double c = maxwell->c;
  double e_fact = maxwell->e_fact, b_fact = maxwell->b_fact;

  double fr[GRM_NEQ], fl[GRM_NEQ];
  gkyl_gr_maxwell_flux(c, e_fact, b_fact, ql, fl);
  gkyl_gr_maxwell_flux(c, e_fact, b_fact, qr, fr);

  for (int m=0; m<GRM_NEQ; ++m) flux_jump[m] = fr[m]-fl[m];
  
  return c;
}

// Flux difference for general relativistic Maxwell
static double
flux_jump_gr(const struct gkyl_wv_eqn *eqn, const double *ql, const double *qr, double *flux_jump)
{
  bool in_excision_region = (ql[GRM_IER] < 0.0) || (qr[GRM_IER] < 0.0);
  if (in_excision_region) {
    for (int m=0; m<GRM_NEQ; ++m) flux_jump[m] = 0.0;
    return 0.0;
  }
  
  const struct wv_gr_maxwell *maxwell = container_of(eqn, struct wv_gr_maxwell, eqn);

  const double c = maxwell->c;
  const double e_fact = maxwell->e_fact, b_fact = maxwell->b_fact;

  double fl[GRM_NEQ], fr[GRM_NEQ];
  gkyl_gr_maxwell_flux(c, e_fact, b_fact, ql, fl);
  gkyl_gr_maxwell_flux(c, e_fact, b_fact, qr, fr);
  
  double fl_extrap[GRM_NEQ], fr_extrap[GRM_NEQ];
  gkyl_gr_maxwell_extrapolate_flux(c, e_fact, b_fact, ql, fl, fl_extrap);
  gkyl_gr_maxwell_extrapolate_flux(c, e_fact, b_fact, qr, fr, fr_extrap);

  for (int m=0; m<GRM_NEQ; ++m) flux_jump[m] = fr_extrap[m]-fl_extrap[m];
  
  const double msl = gkyl_gr_maxwell_max_abs_speed(c, e_fact, b_fact, ql);
  const double msr = gkyl_gr_maxwell_max_abs_speed(c, e_fact, b_fact, qr);
  return fmax(msl, msr);
}

static double
wave_lax(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr, 
  double *waves, double *s)
{
  bool in_excision_region = (ql[GRM_IER] < 0.0) || (qr[GRM_IER] < 0.0);
  if (in_excision_region) {
    for (int i=0; i<GRM_NEQ*2; ++i) waves[i] = 0.0;
    s[0] = s[1] = 0.0;
    return 0.0;
  }
  
  double df[GRM_NEQ];
  const double max_spd = flux_jump_gr(eqn, ql, qr, df);
  
  double *w0 = &waves[0]; double *w1 = &waves[GRM_NEQ];
  
  for (int i = 0; i < GRM_LAPSE; ++i) {
    w0[i] = 0.5*((qr[i] - ql[i]) - df[i]/max_spd);
    w1[i] = 0.5*((qr[i] - ql[i]) + df[i]/max_spd);
  }
  
  for (int i = GRM_LAPSE; i < GRM_NEQ; ++i) {
    w0[i] = w1[i] = 0.0;
  }
  
  // Speed rescaling due to non-Euclidean geometry (can be handled in geometry code?)
  const double gam_22 = 0.5*(ql[GRM_GAM22] + qr[GRM_GAM22]);
  const double gam_33 = 0.5*(ql[GRM_GAM33] + qr[GRM_GAM33]);
  const double gam_23 = 0.5*(ql[GRM_GAM23] + qr[GRM_GAM23]);
  const double lenr = sqrt(gam_22*gam_33 - gam_23*gam_23);
  
  s[0] = -max_spd*lenr;
  s[1] = max_spd*lenr;

  return max_spd;
}

static void
qfluct_lax(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  const double *w0 = &waves[0], *w1 = &waves[GRM_NEQ];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i = 0; i < GRM_LAPSE; ++i) {
    amdq[i] = s0m*w0[i] + s1m*w1[i];
    apdq[i] = s0p*w0[i] + s1p*w1[i];
  }
  for (int i = GRM_LAPSE; i < GRM_NEQ; ++i) {
    amdq[i] = 0.0;
    apdq[i] = 0.0;
  }
}

//static double
//wave(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  //const double *delta, const double *ql, const double *qr, 
  //double *waves, double *s)
//{
  //const struct wv_gr_maxwell *maxwell = container_of(eqn, struct wv_gr_maxwell, eqn);

  //double c = maxwell->c, c1 = 1/c;
  //double e_fact = maxwell->e_fact, b_fact = maxwell->b_fact;
    
  //// compute projections of jump
  //double a1 = 0.5*(delta[3]-delta[7]*c1);
  //double a2 = 0.5*(delta[3]+delta[7]*c1);
  //double a3 = 0.5*(delta[0]-delta[6]*c);
  //double a4 = 0.5*(delta[0]+delta[6]*c);
  //double a5 = 0.5*(delta[1]-delta[5]*c);
  //double a6 = 0.5*(delta[4]*c+delta[2]);
  //double a7 = 0.5*(delta[5]*c+delta[1]);
  //double a8 = 0.5*(delta[2]-delta[4]*c);

  //// set waves to 0.0 as most entries vanish
  //for (int i=0; i<23*6; ++i) waves[i] = 0.0;

  //double *w = 0;

  //// wave 1:
  //w = &waves[0*23];
  //w[3] = a1;
  //w[7] = -a1*c;
  //s[0] = -c*b_fact;

  //// wave 2:
  //w = &waves[1*23];
  //w[3] = a2;
  //w[7] = a2*c;
  //s[1] = c*b_fact;

  //// wave 3:
  //w = &waves[2*23];
  //w[0] = a3;
  //w[6] = -a3*c1;
  //s[2] = -c*e_fact;

  //// wave 4:
  //w = &waves[3*23];
  //w[0] = a4;
  //w[6] = a4*c1;
  //s[3] = c*e_fact;

  //// wave 5: (two waves with EV -c, -c lumped into one)
  //w = &waves[4*23];
  //w[1] = a5;
  //w[2] = a6;
  //w[4] = a6*c1;
  //w[5] = -a5*c1;
  //s[4] = -c;

  //// wave 6: (two waves with EV c, c lumped into one)
  //w = &waves[5*23];
  //w[1] = a7;
  //w[2] = a8;
  //w[4] = -a8*c1;
  //w[5] = a7*c1;
  //s[5] = c;
  
  //return c;
//}

//static void
//qfluct(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  //const double *ql, const double *qr, const double *waves, const double *s,
  //double *amdq, double *apdq)
//{
  //const double *w0 = &waves[0*23], *w1 = &waves[1*23], *w2 = &waves[2*23];
  //const double *w3 = &waves[3*23], *w4 = &waves[4*23], *w5 = &waves[5*23];
  
  //double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]), s2m = fmin(0.0, s[2]);
  //double s3m = fmin(0.0, s[3]), s4m = fmin(0.0, s[4]), s5m = fmin(0.0, s[5]);
  
  //double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]), s2p = fmax(0.0, s[2]);
  //double s3p = fmax(0.0, s[3]), s4p = fmax(0.0, s[4]), s5p = fmax(0.0, s[5]);

  //for (int i=0; i<23; ++i) {
    //amdq[i] = s0m*w0[i] + s1m*w1[i] + s2m*w2[i] + s3m*w3[i] + s4m*w4[i] + s5m*w5[i];
    //apdq[i] = s0p*w0[i] + s1p*w1[i] + s2p*w2[i] + s3p*w3[i] + s4p*w4[i] + s5p*w5[i];
  //}
//}

static bool
check_inv(const struct gkyl_wv_eqn *eqn, const double *q)
{
  return true;
}

static double
max_speed(const struct gkyl_wv_eqn *eqn, const double *q)
{
  const struct wv_gr_maxwell *maxwell = container_of(eqn, struct wv_gr_maxwell, eqn);
  if (q[GRM_IER] < 0.0) {
    return 0.0;
  }
  return maxwell->c / fabs(q[GRM_LAPSE] * q[GRM_DET]);
}

static inline void
gr_maxwell_cons_to_diag(const struct gkyl_wv_eqn *eqn,
  const double *qin, double *diag)
{
  // components of EM energy
  // TODO: update to account for modified components (sqrt(gam)*D, sqrt(gam)*B)
  for (int i=0; i<GRM_LAPSE; ++i) diag[i] = qin[i]*qin[i];
}
