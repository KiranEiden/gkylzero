#pragma once

// Private header, not for direct use in user code

#include <math.h>
#include <gkyl_array.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <stdio.h>

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
gkyl_gr_maxwell_max_abs_speed(double c, double e_fact, double b_fact, const double q[21])
{
  if (q[20] < 0.0)
  {
    return 0.0;
  }
  
  double spd = 0.0;
  double spd_d;
  for(int d = 0; d < 3; d++)
  {
    spd_d = fmax(fabs(q[6] - q[7+d]), fabs(q[6] + q[7+d]));
    spd = spd < spd_d ? spd_d : spd;
  }
  
  return c*spd;
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
gkyl_gr_maxwell_flux(double c, double e_fact, double b_fact, const double q[21], double flux[21])
{
  if (q[20] < 0.0)
  {
    // Fluxes are 0 in excision region
    for (int i = 0; i < 21; ++i) flux[i] = 0.0;
    return;
  }
  
  const double c2 = c*c;
  
  double D[3];
  gkyl_gr_maxwell_rm_rt_det(&q[0], D, q[10]);
  double B[3];
  gkyl_gr_maxwell_rm_rt_det(&q[3], B, q[10]);

  //flux[0] = e_fact*c2*q[6]; // e_fact*c^2*phi
  flux[0] = 0.0;
  flux[1] = c2*B[2]; // c^2*Bz
  flux[2] = -c2*B[1]; // -c^2*By
  //flux[3] = b_fact*q[7]; // b_fact*psi
  flux[3] = 0.0;
  flux[4] = -D[2]; // -Ez
  flux[5] = D[1]; // Ey
  //flux[6] = e_fact*q[0]; // e_fact*Ex
  //flux[7] = b_fact*c2*q[3]; // b_fact*c^2*Bx
  
  for (int i = 6; i < 21; ++i) flux[i] = 0.0;
}

static void
gkyl_gr_maxwell_extrapolate_flux(double c, double e_fact, double b_fact,
    const double q[21], const double flux_in[21], double flux_extrap[21])
{
  if (q[20] < 0.0)
  {
    // Just copy fluxes if in excision region
    for (int i = 0; i < 21; ++i) flux_extrap[i] = flux_in[i];
    return;
  }
  
  const double c2 = c*c;
  const double lapse = q[6];
  const double rt_gam_det_inv = 1. / q[10];
  
  double D[3], D_cov[3];
  double B[3], B_cov[3];
  double shift_cov[3];
  
  D[0] = q[0] * rt_gam_det_inv;
  D[1] = flux_in[5];
  D[2] = -flux_in[4];
  
  B[0] = q[3] * rt_gam_det_inv;
  B[1] = flux_in[1] / c2;
  B[2] = -flux_in[2] / c2;
  
  gkyl_gr_maxwell_lower_ind(D, D_cov, &q[11]);
  gkyl_gr_maxwell_lower_ind(B, B_cov, &q[11]);
  gkyl_gr_maxwell_lower_ind(&q[7], shift_cov, &q[11]);
  
  flux_extrap[0] = flux_in[0]; // 0 for Maxwell curl equations
  flux_extrap[1] = c2 * (lapse*B_cov[2] + rt_gam_det_inv*(shift_cov[1]*D_cov[0] - shift_cov[0]*D_cov[1]));
  flux_extrap[2] = c2 * (lapse*B_cov[1] + rt_gam_det_inv*(shift_cov[0]*D_cov[2] - shift_cov[2]*D_cov[0]));
  flux_extrap[3] = flux_in[3]; // 0 for Maxwell curl equations
  flux_extrap[4] = lapse*D_cov[2] + rt_gam_det_inv*(shift_cov[0]*B_cov[1] - shift_cov[1]*B_cov[0]);
  flux_extrap[5] = lapse*D_cov[1] + rt_gam_det_inv*(shift_cov[2]*B_cov[0] - shift_cov[0]*B_cov[2]);
  
  // Include cell volume correction here (can be handled in geometry code)?
  for (int i = 0; i < 6; ++i) {
    flux_extrap[i] *= rt_gam_det_inv;
  }
  
  for (int i = 6; i < 21; ++i) flux_extrap[i] = 0.0;
}

static inline void
cons_to_riem(const struct gkyl_wv_eqn *eqn,
  const double *qstate, const double *qin, double *wout)
{
  // TODO: this should use proper L matrix
  for (int i=0; i<21; ++i)
    wout[i] = qin[i];
}

static inline void
riem_to_cons(const struct gkyl_wv_eqn *eqn,
  const double *qstate, const double *win, double *qout)
{
  // TODO: this should use proper L matrix
  for (int i=0; i<21; ++i)
    qout[i] = win[i];
}

static void
gr_maxwell_wall(double t, int nc, const double *skin, double * GKYL_RESTRICT ghost, void *ctx)
{
  // zero-tangent for E field
  ghost[0] = skin[0];
  ghost[1] = -skin[1];
  ghost[2] = -skin[2];

  // zero-normal for B field
  ghost[3] = -skin[3];
  ghost[4] = skin[4];
  ghost[5] = skin[5];
  
  for (int i = 6; i < 21; ++i) ghost[i] = skin[i];

  // correction potential
  //ghost[6] = -skin[6];
  //ghost[7] = skin[7];
}

static inline void
rot_to_local(const double *tau1, const double *tau2, const double *norm,
  const double *GKYL_RESTRICT qglobal, double *GKYL_RESTRICT qlocal)
{
  // Rotate D to local coordinates
  qlocal[0] = (qglobal[0] * norm[0]) + (qglobal[1] * norm[1]) + (qglobal[2] * norm[2]);
  qlocal[1] = (qglobal[0] * tau1[0]) + (qglobal[1] * tau1[1]) + (qglobal[2] * tau1[2]);
  qlocal[2] = (qglobal[0] * tau2[0]) + (qglobal[1] * tau2[1]) + (qglobal[2] * tau2[2]);
  
  // Rotate B to local coordinates
  qlocal[3] = (qglobal[3] * norm[0]) + (qglobal[4] * norm[1]) + (qglobal[5] * norm[2]);
  qlocal[4] = (qglobal[3] * tau1[0]) + (qglobal[4] * tau1[1]) + (qglobal[5] * tau1[2]);
  qlocal[5] = (qglobal[3] * tau2[0]) + (qglobal[4] * tau2[1]) + (qglobal[5] * tau2[2]);
  
  // Copy scalar lapse function
  qlocal[6] = qglobal[6];
  
  // Rotate shift vector (beta) to local coordinates
  qlocal[7] = (qglobal[7] * norm[0]) + (qglobal[8] * norm[1]) + (qglobal[9] * norm[2]);
  qlocal[8] = (qglobal[7] * tau1[0]) + (qglobal[8] * tau1[1]) + (qglobal[9] * tau1[2]);
  qlocal[9] = (qglobal[7] * tau2[0]) + (qglobal[8] * tau2[1]) + (qglobal[9] * tau2[2]);
  
  // Copy the sqrt of the spatial metric determinant
  qlocal[10] = qglobal[10];
  
  // Transform the spatial metric (stored in row-major order)
  double x[9]; // temporary array
  x[0] = (qglobal[11] * norm[0]) + (qglobal[12] * norm[1]) + (qglobal[13] * norm[2]);
  x[1] = (qglobal[11] * tau1[0]) + (qglobal[12] * tau1[1]) + (qglobal[13] * tau1[2]);
  x[2] = (qglobal[11] * tau2[0]) + (qglobal[12] * tau2[1]) + (qglobal[13] * tau2[2]);
  x[3] = (qglobal[14] * norm[0]) + (qglobal[15] * norm[1]) + (qglobal[16] * norm[2]);
  x[4] = (qglobal[14] * tau1[0]) + (qglobal[15] * tau1[1]) + (qglobal[16] * tau1[2]);
  x[5] = (qglobal[14] * tau2[0]) + (qglobal[15] * tau2[1]) + (qglobal[16] * tau2[2]);
  x[6] = (qglobal[17] * norm[0]) + (qglobal[18] * norm[1]) + (qglobal[19] * norm[2]);
  x[7] = (qglobal[17] * tau1[0]) + (qglobal[18] * tau1[1]) + (qglobal[19] * tau1[2]);
  x[8] = (qglobal[17] * tau2[0]) + (qglobal[18] * tau2[1]) + (qglobal[19] * tau2[2]);
  
  qlocal[11] = (x[0] * norm[0]) + (x[3] * norm[1]) + (x[6] * norm[2]);
  qlocal[12] = (x[0] * tau1[0]) + (x[3] * tau1[1]) + (x[6] * tau1[2]);
  qlocal[13] = (x[0] * tau2[0]) + (x[3] * tau2[1]) + (x[6] * tau2[2]);
  qlocal[14] = (x[1] * norm[0]) + (x[4] * norm[1]) + (x[7] * norm[2]);
  qlocal[15] = (x[1] * tau1[0]) + (x[4] * tau1[1]) + (x[7] * tau1[2]);
  qlocal[16] = (x[1] * tau2[0]) + (x[4] * tau2[1]) + (x[7] * tau2[2]);
  qlocal[17] = (x[2] * norm[0]) + (x[5] * norm[1]) + (x[8] * norm[2]);
  qlocal[18] = (x[2] * tau1[0]) + (x[5] * tau1[1]) + (x[8] * tau1[2]);
  qlocal[19] = (x[2] * tau2[0]) + (x[5] * tau2[1]) + (x[8] * tau2[2]);
  
  // Copy in_excision_region flag
  qlocal[20] = qglobal[20];
  
  // Correction potentials are scalars and unchanged
  //qlocal[6] = qglobal[6];
  //qlocal[7] = qglobal[7];
}

static inline void
rot_to_global(const double *tau1, const double *tau2, const double *norm,
  const double *GKYL_RESTRICT qlocal, double *GKYL_RESTRICT qglobal)
{
  // Rotate D back to global coordinates
  qglobal[0] = (qlocal[0] * norm[0]) + (qlocal[1] * tau1[0]) + (qlocal[2] * tau2[0]);
  qglobal[1] = (qlocal[0] * norm[1]) + (qlocal[1] * tau1[1]) + (qlocal[2] * tau2[1]);
  qglobal[2] = (qlocal[0] * norm[2]) + (qlocal[1] * tau1[2]) + (qlocal[2] * tau2[2]);
  
  // Rotate B back to global coordinates
  qglobal[3] = (qlocal[3] * norm[0]) + (qlocal[4] * tau1[0]) + (qlocal[5] * tau2[0]);
  qglobal[4] = (qlocal[3] * norm[1]) + (qlocal[4] * tau1[1]) + (qlocal[5] * tau2[1]);
  qglobal[5] = (qlocal[3] * norm[2]) + (qlocal[4] * tau1[2]) + (qlocal[5] * tau2[2]);
  
  // Again just copy scalar lapse function
  qglobal[6] = qlocal[6];
  
  // Rotate shift vector (beta) to back to global coordinates
  qglobal[7] = (qlocal[7] * norm[0]) + (qlocal[8] * tau1[0]) + (qlocal[9] * tau2[0]);
  qglobal[8] = (qlocal[7] * norm[1]) + (qlocal[8] * tau1[1]) + (qlocal[9] * tau2[1]);
  qglobal[9] = (qlocal[7] * norm[2]) + (qlocal[8] * tau1[2]) + (qlocal[9] * tau2[2]);
  
  // Copy sqrt of spatial metric determinant
  qglobal[10] = qlocal[10];
  
  // Transform the spatial metric
  double x[9]; // again a temporary array
  x[0] = (qglobal[11] * norm[0]) + (qglobal[12] * tau1[0]) + (qglobal[13] * tau2[0]);
  x[1] = (qglobal[11] * norm[1]) + (qglobal[12] * tau1[1]) + (qglobal[13] * tau2[1]);
  x[2] = (qglobal[11] * norm[2]) + (qglobal[12] * tau1[2]) + (qglobal[13] * tau2[2]);
  x[3] = (qglobal[14] * norm[0]) + (qglobal[15] * tau1[0]) + (qglobal[16] * tau2[0]);
  x[4] = (qglobal[14] * norm[1]) + (qglobal[15] * tau1[1]) + (qglobal[16] * tau2[1]);
  x[5] = (qglobal[14] * norm[2]) + (qglobal[15] * tau1[2]) + (qglobal[16] * tau2[2]);
  x[6] = (qglobal[17] * norm[0]) + (qglobal[18] * tau1[0]) + (qglobal[19] * tau2[0]);
  x[7] = (qglobal[17] * norm[1]) + (qglobal[18] * tau1[1]) + (qglobal[19] * tau2[1]);
  x[8] = (qglobal[17] * norm[2]) + (qglobal[18] * tau1[2]) + (qglobal[19] * tau2[2]);
  
  qglobal[11] = (x[0] * norm[0]) + (x[3] * tau1[0]) + (x[6] * tau2[0]);
  qglobal[12] = (x[0] * norm[1]) + (x[3] * tau1[1]) + (x[6] * tau2[1]);
  qglobal[13] = (x[0] * norm[2]) + (x[3] * tau1[2]) + (x[6] * tau2[2]);
  qglobal[14] = (x[1] * norm[0]) + (x[4] * tau1[0]) + (x[7] * tau2[0]);
  qglobal[15] = (x[1] * norm[1]) + (x[4] * tau1[1]) + (x[7] * tau2[1]);
  qglobal[16] = (x[1] * norm[2]) + (x[4] * tau1[2]) + (x[7] * tau2[2]);
  qglobal[17] = (x[2] * norm[0]) + (x[5] * tau1[0]) + (x[8] * tau2[0]);
  qglobal[18] = (x[2] * norm[1]) + (x[5] * tau1[1]) + (x[8] * tau2[1]);
  qglobal[19] = (x[2] * norm[2]) + (x[5] * tau1[2]) + (x[8] * tau2[2]);
  
  // Copy the in_excision_region flag
  qglobal[20] = qlocal[20];
  
  // Correction potentials are scalars and unchanged
  //qglobal[6] = qlocal[6];
  //qglobal[7] = qlocal[7];
}

// Flux difference for special relativistic Maxwell
static double
flux_jump(const struct gkyl_wv_eqn *eqn, const double *ql, const double *qr, double *flux_jump)
{
  bool in_excision_region = (ql[20] < 0.0) || (qr[20] < 0.0);
  if (in_excision_region) {
    for (int m=0; m<21; ++m) flux_jump[m] = 0.0;
    return 0.0;
  }
  
  const struct wv_gr_maxwell *maxwell = container_of(eqn, struct wv_gr_maxwell, eqn);

  double c = maxwell->c;
  double e_fact = maxwell->e_fact, b_fact = maxwell->b_fact;

  double fr[21], fl[21];
  gkyl_gr_maxwell_flux(c, e_fact, b_fact, ql, fl);
  gkyl_gr_maxwell_flux(c, e_fact, b_fact, qr, fr);

  for (int m=0; m<21; ++m) flux_jump[m] = fr[m]-fl[m];
  
  return c;
}

// Flux difference for general relativistic Maxwell
static double
flux_jump_gr(const struct gkyl_wv_eqn *eqn, const double *ql, const double *qr, double *flux_jump)
{
  bool in_excision_region = (ql[20] < 0.0) || (qr[20] < 0.0);
  if (in_excision_region) {
    for (int m=0; m<21; ++m) flux_jump[m] = 0.0;
    return 0.0;
  }
  
  const struct wv_gr_maxwell *maxwell = container_of(eqn, struct wv_gr_maxwell, eqn);

  const double c = maxwell->c;
  const double e_fact = maxwell->e_fact, b_fact = maxwell->b_fact;

  double fl[21], fr[21];
  gkyl_gr_maxwell_flux(c, e_fact, b_fact, ql, fl);
  gkyl_gr_maxwell_flux(c, e_fact, b_fact, qr, fr);
  
  double fl_extrap[21], fr_extrap[21];
  gkyl_gr_maxwell_extrapolate_flux(c, e_fact, b_fact, ql, fl, fl_extrap);
  gkyl_gr_maxwell_extrapolate_flux(c, e_fact, b_fact, qr, fr, fr_extrap);

  for (int m=0; m<21; ++m) flux_jump[m] = fr_extrap[m]-fl_extrap[m];
  
  const double msl = gkyl_gr_maxwell_max_abs_speed(c, e_fact, b_fact, ql);
  const double msr = gkyl_gr_maxwell_max_abs_speed(c, e_fact, b_fact, qr);
  return fmax(msl, msr);
}

static double
wave_lax(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr, 
  double *waves, double *s)
{
  bool in_excision_region = (ql[20] < 0.0) || (qr[20] < 0.0);
  if (in_excision_region) {
    for (int i=0; i<21*2; ++i) waves[i] = 0.0;
    s[0] = s[1] = 0.0;
    return 0.0;
  }
  
  double df[21];
  const double max_spd = flux_jump_gr(eqn, ql, qr, df);
  
  double *w0 = &waves[0]; double *w1 = &waves[21];
  
  for (int i = 0; i < 6; ++i) {
    w0[i] = 0.5*((qr[i] - ql[i]) - df[i]/max_spd);
    w1[i] = 0.5*((qr[i] - ql[i]) + df[i]/max_spd);
  }
  
  for (int i = 6; i < 21; ++i) {
    w0[i] = w1[i] = 0.0;
  }
  
  // Speed rescaling due to non-Euclidean geometry (can be handled in geometry code?)
  const double gam_22 = 0.5*(ql[15] + qr[15]);
  const double gam_33 = 0.5*(ql[19] + qr[19]);
  const double gam_23 = 0.5*(ql[16] + qr[16]);
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
  const double *w0 = &waves[0], *w1 = &waves[21];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i = 0; i < 6; ++i) {
    amdq[i] = s0m*w0[i] + s1m*w1[i];
    apdq[i] = s0p*w0[i] + s1p*w1[i];
  }
  for (int i = 6; i < 21; ++i) {
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
  //for (int i=0; i<21*6; ++i) waves[i] = 0.0;

  //double *w = 0;

  //// wave 1:
  //w = &waves[0*21];
  //w[3] = a1;
  //w[7] = -a1*c;
  //s[0] = -c*b_fact;

  //// wave 2:
  //w = &waves[1*21];
  //w[3] = a2;
  //w[7] = a2*c;
  //s[1] = c*b_fact;

  //// wave 3:
  //w = &waves[2*21];
  //w[0] = a3;
  //w[6] = -a3*c1;
  //s[2] = -c*e_fact;

  //// wave 4:
  //w = &waves[3*21];
  //w[0] = a4;
  //w[6] = a4*c1;
  //s[3] = c*e_fact;

  //// wave 5: (two waves with EV -c, -c lumped into one)
  //w = &waves[4*21];
  //w[1] = a5;
  //w[2] = a6;
  //w[4] = a6*c1;
  //w[5] = -a5*c1;
  //s[4] = -c;

  //// wave 6: (two waves with EV c, c lumped into one)
  //w = &waves[5*21];
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
  //const double *w0 = &waves[0*21], *w1 = &waves[1*21], *w2 = &waves[2*21];
  //const double *w3 = &waves[3*21], *w4 = &waves[4*21], *w5 = &waves[5*21];
  
  //double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]), s2m = fmin(0.0, s[2]);
  //double s3m = fmin(0.0, s[3]), s4m = fmin(0.0, s[4]), s5m = fmin(0.0, s[5]);
  
  //double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]), s2p = fmax(0.0, s[2]);
  //double s3p = fmax(0.0, s[3]), s4p = fmax(0.0, s[4]), s5p = fmax(0.0, s[5]);

  //for (int i=0; i<21; ++i) {
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
  return maxwell->c;
}

static inline void
gr_maxwell_cons_to_diag(const struct gkyl_wv_eqn *eqn,
  const double *qin, double *diag)
{
  // components of EM energy
  // TODO: update to account for modified components (sqrt(gam)*D, sqrt(gam)*B)
  for (int i=0; i<6; ++i) diag[i] = qin[i]*qin[i];
}
