#include <stdbool.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_moment_prim_gr_maxwell.h>
#include <gkyl_wv_gr_maxwell.h>

#ifndef GR_MAXWELL_EPS
#define GR_MAXWELL_EPS 1e-16
#endif

struct wv_gr_maxwell {
  struct gkyl_wv_eqn eqn; // base object
  //double c; // light speed
  //double e_fact, b_fact; // electric and magnetic correction factors
};

static void
gr_maxwell_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_wv_eqn *base = container_of(ref, struct gkyl_wv_eqn, ref_count);
  struct wv_gr_maxwell *gr_maxwell = container_of(base, struct wv_gr_maxwell, eqn);
  gkyl_free(gr_maxwell);
}

static inline void
cons_to_riem(const struct gkyl_wv_eqn *eqn,
  const double *qstate, const double *qin, double *wout)
{
  // TODO: this should use proper L matrix
  for (int i=0; i<11; ++i)
    wout[i] = qin[i];
}
static inline void
riem_to_cons(const struct gkyl_wv_eqn *eqn,
  const double *qstate, const double *win, double *qout)
{
  // TODO: this should use proper L matrix
  for (int i=0; i<11; ++i)
    qout[i] = win[i];
}

static void
gr_maxwell_wall(double t, int nc, const double *skin, double * GKYL_RESTRICT ghost, void *ctx)
{
  // TODO: double-check for D-field and different geometry
  
  // zero-tangent for D field
  ghost[0] = skin[0];
  ghost[1] = -skin[1];
  ghost[2] = -skin[2];

  // zero-normal for B field
  ghost[3] = -skin[3];
  ghost[4] = skin[4];
  ghost[5] = skin[5];

  // correction potential
  ghost[6] = -skin[6];
  ghost[7] = skin[7];
}

static inline void
rot_to_local(const double *tau1, const double *tau2, const double *norm,
  const double *GKYL_RESTRICT qglobal, double *GKYL_RESTRICT qlocal)
{
  // TODO: Calculate actual norm, tau1, tau2 using covariant basis vector components
  
  // Rotate D to local coordinates
  qlocal[0] = qglobal[0]*norm[0] + qglobal[1]*norm[1] + qglobal[2]*norm[2];
  qlocal[1] = qglobal[0]*tau1[0] + qglobal[1]*tau1[1] + qglobal[2]*tau1[2];
  qlocal[2] = qglobal[0]*tau2[0] + qglobal[1]*tau2[1] + qglobal[2]*tau2[2];
  // Rotate B to local coordinates
  qlocal[3] = qglobal[3]*norm[0] + qglobal[4]*norm[1] + qglobal[5]*norm[2];
  qlocal[4] = qglobal[3]*tau1[0] + qglobal[4]*tau1[1] + qglobal[5]*tau1[2];
  qlocal[5] = qglobal[3]*tau2[0] + qglobal[4]*tau2[1] + qglobal[5]*tau2[2];
  // Just copy scalar lapse function
  qlocal[6] = qglobal[6];
  // Rotate shift vector (beta) to local coordinates
  qlocal[7] = qglobal[7]*norm[0] + qglobal[8]*norm[1] + qglobal[9]*norm[2];
  qlocal[8] = qglobal[7]*tau1[0] + qglobal[8]*tau1[1] + qglobal[9]*tau1[2];
  qlocal[9] = qglobal[7]*tau2[0] + qglobal[8]*tau2[1] + qglobal[9]*tau2[2];
  // Copy the in_excision_region variable
  qlocal[10] = qglobal[10];
  
  // Correction potentials are scalars and unchanged
  //qlocal[6] = qglobal[6];
  //qlocal[7] = qglobal[7];
}

static inline void
rot_to_global(const double *tau1, const double *tau2, const double *norm,
  const double *GKYL_RESTRICT qlocal, double *GKYL_RESTRICT qglobal)
{
  // Rotate D back to global coordinates
  qglobal[0] = qlocal[0]*norm[0] + qlocal[1]*tau1[0] + qlocal[2]*tau2[0];
  qglobal[1] = qlocal[0]*norm[1] + qlocal[1]*tau1[1] + qlocal[2]*tau2[1];
  qglobal[2] = qlocal[0]*norm[2] + qlocal[1]*tau1[2] + qlocal[2]*tau2[2];
  // Rotate B back to global coordinates
  qglobal[3] = qlocal[3]*norm[0] + qlocal[4]*tau1[0] + qlocal[5]*tau2[0];
  qglobal[4] = qlocal[3]*norm[1] + qlocal[4]*tau1[1] + qlocal[5]*tau2[1];
  qglobal[5] = qlocal[3]*norm[2] + qlocal[4]*tau1[2] + qlocal[5]*tau2[2];
  // Again just copy scalar lapse function
  qglobal[6] = qlocal[6];
  // Rotate shift vector (beta) to back to global coordinates
  qglobal[7] = qlocal[7]*norm[0] + qlocal[8]*tau1[0] + qlocal[9]*tau2[0];
  qglobal[8] = qlocal[7]*norm[1] + qlocal[8]*tau1[1] + qlocal[9]*tau2[1];
  qglobal[9] = qlocal[7]*norm[2] + qlocal[8]*tau1[2] + qlocal[9]*tau2[2];
  // Copy the in_excision_region variable
  qglobal[10] = qlocal[10];
  
  // Correction potentials are scalars and unchanged
  //qglobal[6] = qlocal[6];
  //qglobal[7] = qlocal[7];
}

static double
flux_jump(const struct gkyl_wv_eqn *eqn, const double *ql, const double *qr, double *flux_jump)
{
  bool in_excision_region = (ql[10] < GR_MAXWELL_EPS) || (qr[10] < GR_MAXWELL_EPS);
  if (in_excision_region)
  {
    for (int m=0; m<11; ++m) flux_jump[m] = 0.0;
    return GR_MAXWELL_EPS;
  }
  
  double fr[11], fl[11];
  gkyl_gr_maxwell_flux(ql, fl);
  gkyl_gr_maxwell_flux(qr, fr);

  for (int m=0; m<11; ++m) flux_jump[m] = fr[m]-fl[m];

  return fmax(gkyl_gr_maxwell_max_abs_speed(ql), gkyl_gr_maxwell_max_abs_speed(qr));
}

static double
wave_lax(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr, 
  double *waves, double *s)
{
  bool in_excision_region = (ql[10] < GR_MAXWELL_EPS) || (qr[10] < GR_MAXWELL_EPS);
  if (in_excision_region)
  {
    for (int i=0; i<11*2; ++i) waves[i] = 0.0;
    s[0] = s[1] = 0.0;
    return GR_MAXWELL_EPS;
  }
  
  double df[11];
  flux_jump(eqn, ql, qr, df);
  double max_spd_l = fmax(fabs(ql[6] - ql[7]), fabs(ql[6] + ql[7]));
  double max_spd_r = fmax(fabs(qr[6] - qr[7]), fabs(qr[6] + qr[7]));
  double max_spd = fmax(max_spd_l, max_spd_r);
  
  double *w0 = &waves[0]; double *w1 = &waves[11];
  
  for (int i = 0; i < 6; ++i)
  {
    w0[i] = 0.5*((qr[i] - ql[i]) - df[i]/max_spd);
    w1[i] = 0.5*((qr[i] - ql[i]) + df[i]/max_spd);
  }
  
  for (int i = 6; i < 11; ++i)
  {
    w0[i] = w1[i] = 0.0;
  }
  
  s[0] = -max_spd;
  s[1] = max_spd;

  return max_spd;
}

static void
qfluct_lax(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  const double *w0 = &waves[0], *w1 = &waves[11];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i = 0; i < 6; ++i)
  {
    amdq[i] = s0m*w0[i] + s1m*w1[i];
    apdq[i] = s0p*w0[i] + s1p*w1[i];
  }
  for (int i = 6; i < 11; ++i)
  {
    amdq[i] = 0.0;
    apdq[i] = 0.0;
  }
}

// Case where last two shift vector entries are 0 (and the first may or may not be)
static double
wave_zero_s2s3(const double *delta, const double lapse, const double shift1, double *waves, double *s)
{
  // compute projections of jump
  const double a1 = delta[0];
  const double a2 = delta[3];
  const double a3 = 0.5*(delta[2] + delta[4]);
  const double a4 = 0.5*(-delta[1] + delta[5]);
  const double a5 = 0.5*(-delta[2] + delta[4]);
  const double a6 = 0.5*(delta[1] + delta[5]);
  
  // set waves to 0.0 as most entries vanish
  for (int i=0; i<11*3; ++i) waves[i] = 0.0;

  double *w = 0;

  // Two waves with EV 0
  w = &waves[0*11];
  w[0] = a1;
  w[3] = a2;
  s[0] = 0.0;

  // Two waves with EV -lapse - shift1
  w = &waves[1*11];
  w[1] = -a4;
  w[2] = a3;
  w[4] = a3;
  w[5] = a4;
  s[1] = -lapse - shift1;

  // Two waves with EV lapse - shift1
  w = &waves[2*11];
  w[1] = a6;
  w[2] = -a5;
  w[4] = a5;
  w[5] = a6;
  s[2] = lapse - shift1;
  
  return fmax(fabs(lapse - shift1), fabs(lapse + shift1));
}

// Case where the first shift vector entry is 0, and potentially one (but not both) of the other two
static double
wave_zero_s1(const double *delta, const double lapse, const double shift2, const double shift3, 
  double *waves, double *s)
{
  bool no_s2 = fabs(shift2) < GR_MAXWELL_EPS;
  
  // compute projections of jump
  const double a1 = delta[3];
  const double a2 = no_s2 ? (shift3/lapse)*delta[0] : (-shift2/lapse)*delta[0];
  const double a3 = 0.5*((-shift3/lapse)*delta[0] + delta[2] + (-shift2/lapse)*delta[3] + delta[4]);
  const double a4 = 0.5*((shift2/lapse)*delta[0] - delta[1] + (-shift3/lapse)*delta[3] + delta[5]);
  const double a5 = 0.5*((-shift3/lapse)*delta[0] - delta[2] + (shift2/lapse)*delta[3] + delta[4]);
  const double a6 = 0.5*((shift2/lapse)*delta[0] + delta[1] + (shift3/lapse)*delta[3] + delta[5]);

  // set waves to 0.0 as most entries vanish
  for (int i=0; i<11*3; ++i) waves[i] = 0.0;

  double *w = 0;

  // Two waves with EV 0
  w = &waves[0*11];
  w[0] = no_s2 ? (lapse / shift3 * a2) : (-lapse / shift2 * a2);
  w[1] = -shift3 / lapse * a1;
  w[2] = no_s2 ? 0.0 : shift2 / lapse * a1;
  w[3] = a1;
  w[4] = no_s2 ? a2 : (-shift3 / shift2 * a2);
  w[5] = no_s2 ? 0.0 : a2;
  s[0] = 0.0;

  // Two waves with EV -lapse
  w = &waves[1*11];
  w[1] = -a4;
  w[2] = a3;
  w[4] = a3;
  w[5] = a4;
  s[1] = -lapse;

  // Two waves with EV lapse
  w = &waves[2*11];
  w[1] = a6;
  w[2] = -a5;
  w[4] = a5;
  w[5] = a6;
  s[2] = lapse;
  
  return fabs(lapse);
}

static double
wave(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr, 
  double *waves, double *s)
{
  const struct wv_gr_maxwell *gr_maxwell = container_of(eqn, struct wv_gr_maxwell, eqn);
  
  bool in_excision_region = (ql[10] < GR_MAXWELL_EPS) || (qr[10] < GR_MAXWELL_EPS);
  if (in_excision_region)
  {
    for (int i=0; i<11*3; ++i) waves[i] = 0.0;
    s[0] = s[1] = s[2] = 0.0;
    return GR_MAXWELL_EPS;
  }

  const double lapse = 0.5*(ql[6] + qr[6]);
  const double shift1 = 0.5*(ql[7] + qr[7]);
  const double shift2 = 0.5*(ql[8] + qr[8]);
  const double shift3 = 0.5*(ql[9] + qr[9]);
  const double ev2 = -(lapse + shift1); // ev1 is 0
  const double ev3 = lapse - shift1;
  const double ev_prod = ev2*ev3;
  const double kap1 = shift2*shift2 + shift3*shift3;
  const double kap2 = lapse*lapse * shift3*shift3 + shift1*shift1 * shift2*shift2;
  const double kap3 = lapse*lapse * shift2*shift2 + shift1*shift1 * shift3*shift3;
  
  // KE: should we also handle cases where ev2 = 0, ev3 = 0?
  
  if (fabs(kap1) < GR_MAXWELL_EPS)
  {
    // This one should be first, since it also accounts for zero shift1 if needed
    return wave_zero_s2s3(delta, lapse, shift1, waves, s);
  }
  else if (fabs(shift1) < GR_MAXWELL_EPS)
  {
    // Accounts also for either zero shift2 or shift3 being zero, but not both
    return wave_zero_s1(delta, lapse, shift2, shift3, waves, s);
  }
    
  // compute projections of jump
  const double a1 = ((-shift3/(shift1*ev_prod*kap1))*delta[0] + (shift2/(lapse*ev_prod*kap1))*delta[3]);
  const double a2 = ((shift2/(shift1*ev_prod*kap1))*delta[0] + (shift3/(lapse*ev_prod*kap1))*delta[3]);
  const double a3 = 0.5*((shift3/ev2)*delta[0] + delta[2] + (shift2/ev2)*delta[3] + delta[4]);
  const double a4 = 0.5*((-shift2/ev2)*delta[0] - delta[1] + (shift3/ev2)*delta[3] + delta[5]);
  const double a5 = 0.5*((-shift3/ev3)*delta[0] - delta[2] + (shift2/ev3)*delta[3] + delta[4]);
  const double a6 = 0.5*((shift2/ev3)*delta[0] + delta[1] + (shift3/ev3)*delta[3] + delta[5]);

  // set waves to 0.0 as most entries vanish
  for (int i=0; i<11*3; ++i) waves[i] = 0.0;

  double *w = 0;

  // Two waves with EV 0
  w = &waves[0*11];
  w[0] = (-shift1*shift3*ev_prod)*a1 + (shift1*shift2*ev_prod)*a2;
  w[1] = (-shift2*shift3*ev_prod)*a1 + kap2*a2;
  w[2] = -kap3*a1 + (shift2*shift3*ev_prod)*a2;
  w[3] = lapse * ev_prod * (shift2*a1 + shift3*a2);
  w[4] = lapse * shift1 * kap1 * a1;
  w[5] = lapse * shift1 * kap1 * a2;
  s[0] = 0.0;

  // Two waves with EV -lapse - shift1
  w = &waves[1*11];
  w[1] = -a4;
  w[2] = a3;
  w[4] = a3;
  w[5] = a4;
  s[1] = ev2;

  // Two waves with EV lapse - shift1
  w = &waves[2*11];
  w[1] = a6;
  w[2] = -a5;
  w[4] = a5;
  w[5] = a6;
  s[2] = ev3;
  
  return fmax(fabs(ev2), fabs(ev3));
}

static void
qfluct(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  const double *w0 = &waves[0*11], *w1 = &waves[1*11], *w2 = &waves[2*11];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]), s2m = fmin(0.0, s[2]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]), s2p = fmax(0.0, s[2]);

  for (int i=0; i<6; ++i)
  {
    amdq[i] = s0m*w0[i] + s1m*w1[i] + s2m*w2[i];
    apdq[i] = s0p*w0[i] + s1p*w1[i] + s2p*w2[i];
  }
  
  for (int i=6; i<11; ++i)
  {
    amdq[i] = 0.0;
    apdq[i] = 0.0;
  }
}

static void
ffluct(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  for (int i=0; i<11; ++i)
  {
    amdq[i] = 0.0;
    apdq[i] = 0.0;
  }
  
  for (int i=0; i<3; ++i)
  {
    for (int j=0; j<6; ++j)
    {
      if (s[i] < 0.0)
      {
        amdq[j] += waves[i*11 + j];
      }
      else if (s[i] > 0.0)
      {
        apdq[j] += waves[i*11 + j];
      }
      else
      {
        const double cor = 0.5*waves[i*11 + j];
        amdq[j] += cor;
        apdq[j] += cor;
      }
    }
  }
}

static bool
check_inv(const struct gkyl_wv_eqn *eqn, const double *q)
{
  return true; // TODO: double-check that this still works
}

static double
max_speed(const struct gkyl_wv_eqn *eqn, const double *q)
{
  return gkyl_gr_maxwell_max_abs_speed(q);
}

static inline void
gr_maxwell_cons_to_diag(const struct gkyl_wv_eqn *eqn,
  const double *qin, double *diag)
{
  // components of EM energy
  // TODO: verify what this does; convert to handle D instead of E if needed
  for (int i=0; i<6; ++i) diag[i] = qin[i]*qin[i];
}

struct gkyl_wv_eqn*
gkyl_wv_gr_maxwell_new(enum gkyl_wv_gr_maxwell_rp rp_type)
{
  struct wv_gr_maxwell *gr_maxwell = gkyl_malloc(sizeof(struct wv_gr_maxwell));

  gr_maxwell->eqn.type = GKYL_EQN_GR_MAXWELL;
  gr_maxwell->eqn.num_equations = 11; 
  gr_maxwell->eqn.num_diag = 6; // Dx^2, Dy^2, Dz^2, Bx^2, By^2, Bz^2
  
  if(rp_type == WV_GR_MAXWELL_RP_EIG)
  {
    gr_maxwell->eqn.num_waves = 3;
    gr_maxwell->eqn.waves_func = wave;
    gr_maxwell->eqn.qfluct_func = qfluct;
    gr_maxwell->eqn.ffluct_func = ffluct;
  }
  else if(rp_type == WV_GR_MAXWELL_RP_LAX)
  {
    gr_maxwell->eqn.num_waves = 2;
    gr_maxwell->eqn.waves_func = wave_lax;
    gr_maxwell->eqn.qfluct_func = qfluct_lax;
  }
  
  //gr_maxwell->c = c;
  //gr_maxwell->e_fact = e_fact;
  //gr_maxwell->b_fact = b_fact;

  gr_maxwell->eqn.flux_jump = flux_jump;
  gr_maxwell->eqn.check_inv_func = check_inv;
  gr_maxwell->eqn.max_speed_func = max_speed;
  gr_maxwell->eqn.rotate_to_local_func = rot_to_local;
  gr_maxwell->eqn.rotate_to_global_func = rot_to_global;

  gr_maxwell->eqn.cons_to_riem = cons_to_riem;
  gr_maxwell->eqn.riem_to_cons = riem_to_cons;

  gr_maxwell->eqn.wall_bc_func = gr_maxwell_wall;

  gr_maxwell->eqn.cons_to_diag = gr_maxwell_cons_to_diag;

  gr_maxwell->eqn.ref_count = gkyl_ref_count_init(gr_maxwell_free);

  return &gr_maxwell->eqn;
}
