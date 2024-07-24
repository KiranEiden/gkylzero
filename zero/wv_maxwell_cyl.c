#include <stdbool.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_moment_prim_maxwell_cyl.h>
#include <gkyl_wv_maxwell_cyl.h>

struct wv_maxwell_cyl {
  struct gkyl_wv_eqn eqn; // base object
  double c; // light speed
};

static void
maxwell_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_wv_eqn *base = container_of(ref, struct gkyl_wv_eqn, ref_count);
  struct wv_maxwell_cyl *maxwell = container_of(base, struct wv_maxwell_cyl, eqn);
  gkyl_free(maxwell);
}

static double
flux_jump(const struct gkyl_wv_eqn *eqn, const double *ql, const double *qr, double *flux_jump)
{
  const struct wv_maxwell_cyl *maxwell = container_of(eqn, struct wv_maxwell_cyl, eqn);

  double c = maxwell->c;

  double fr[9], fl[9];
  gkyl_maxwell_cyl_flux(c, ql, fl);
  gkyl_maxwell_cyl_flux(c, qr, fr);

  for (int m=0; m<9; ++m) flux_jump[m] = fr[m]-fl[m];

  double msl = gkyl_maxwell_cyl_max_abs_speed(c, ql);
  double msr = gkyl_maxwell_cyl_max_abs_speed(c, qr);
  return fmax(msl, msr);
}

static inline void
cons_to_riem(const struct gkyl_wv_eqn *eqn,
  const double *qstate, const double *qin, double *wout)
{
  // TODO: this should use proper L matrix
  for (int i=0; i<9; ++i)
    wout[i] = qin[i];
}
static inline void
riem_to_cons(const struct gkyl_wv_eqn *eqn,
  const double *qstate, const double *win, double *qout)
{
  // TODO: this should use proper L matrix
  for (int i=0; i<9; ++i)
    qout[i] = win[i];
}

static void
maxwell_wall(double t, int nc, const double *skin, double * GKYL_RESTRICT ghost, void *ctx)
{
  // zero-tangent for E field
  ghost[0] = skin[0];
  ghost[1] = -skin[1];
  ghost[2] = -skin[2];

  // zero-normal for B field
  ghost[3] = -skin[3];
  ghost[4] = skin[4];
  ghost[5] = skin[5];
  
  // TODO: Double-check what this function does
  ghost[6] = skin[6];
  ghost[7] = skin[7];
  ghost[8] = skin[8];
}

static inline void
rot_to_local(const double *tau1, const double *tau2, const double *norm,
  const double *GKYL_RESTRICT qglobal, double *GKYL_RESTRICT qlocal)
{
  // Rotate E to local coordinates
  qlocal[0] = qglobal[0]*norm[0] + qglobal[1]*norm[1] + qglobal[2]*norm[2];
  qlocal[1] = qglobal[0]*tau1[0] + qglobal[1]*tau1[1] + qglobal[2]*tau1[2];
  qlocal[2] = qglobal[0]*tau2[0] + qglobal[1]*tau2[1] + qglobal[2]*tau2[2];
  // Rotate B to local coordinates
  qlocal[3] = qglobal[3]*norm[0] + qglobal[4]*norm[1] + qglobal[5]*norm[2];
  qlocal[4] = qglobal[3]*tau1[0] + qglobal[4]*tau1[1] + qglobal[5]*tau1[2];
  qlocal[5] = qglobal[3]*tau2[0] + qglobal[4]*tau2[1] + qglobal[5]*tau2[2];
  // Rotate sqrts of metric diagonal (we want q[6] = r for azimuthal direction, 1 otherwise)
  qlocal[6] = qglobal[6]*norm[0] + qglobal[7]*norm[1] + qglobal[8]*norm[2];
  qlocal[7] = qglobal[6]*tau1[0] + qglobal[7]*tau1[1] + qglobal[8]*tau1[2];
  qlocal[8] = qglobal[6]*tau2[0] + qglobal[7]*tau2[1] + qglobal[8]*tau2[2];
}

static inline void
rot_to_global(const double *tau1, const double *tau2, const double *norm,
  const double *GKYL_RESTRICT qlocal, double *GKYL_RESTRICT qglobal)
{
  // E-field waves/fluctuations
  qglobal[0] = qlocal[0]*norm[0] + qlocal[1]*tau1[0] + qlocal[2]*tau2[0];
  qglobal[1] = qlocal[0]*norm[1] + qlocal[1]*tau1[1] + qlocal[2]*tau2[1];
  qglobal[2] = qlocal[0]*norm[2] + qlocal[1]*tau1[2] + qlocal[2]*tau2[2];
  // B-field waves/fluctuations
  qglobal[3] = qlocal[3]*norm[0] + qlocal[4]*tau1[0] + qlocal[5]*tau2[0];
  qglobal[4] = qlocal[3]*norm[1] + qlocal[4]*tau1[1] + qlocal[5]*tau2[1];
  qglobal[5] = qlocal[3]*norm[2] + qlocal[4]*tau1[2] + qlocal[5]*tau2[2];
  // Waves/fluctuations in sqrts of metric diagonal
  qglobal[6] = qlocal[6]*norm[0] + qlocal[7]*tau1[0] + qlocal[8]*tau2[0];
  qglobal[7] = qlocal[6]*norm[1] + qlocal[7]*tau1[1] + qlocal[8]*tau2[1];
  qglobal[8] = qlocal[6]*norm[2] + qlocal[7]*tau1[2] + qlocal[8]*tau2[2];
}

static double
wave_lax(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr, 
  double *waves, double *s)
{
  double df[9];
  double max_spd = flux_jump(eqn, ql, qr, df);
  
  double rl = fabs(ql[6]*ql[7]*ql[8]);
  double rr = fabs(qr[6]*qr[7]*qr[8]);
  double r = 0.5*(rl + rr);
  
  double *w0 = &waves[0]; double *w1 = &waves[9];
  
  for (int i = 0; i < 6; ++i)
  {
    w0[i] = 0.5*((qr[i] - ql[i]) - df[i]/r/max_spd);
    w1[i] = 0.5*((qr[i] - ql[i]) + df[i]/r/max_spd);
  }
  
  for (int i = 6; i < 9; ++i)
  {
    w0[i] = w1[i] = 0.0;
  }
  
  s[0] = -max_spd;
  s[1] = max_spd;

  return max_spd;
}

static double
wave(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr, 
  double *waves, double *s)
{
  const struct wv_maxwell_cyl *maxwell = container_of(eqn, struct wv_maxwell_cyl, eqn);

  double c = maxwell->c, c1 = 1.0 / c;
  //double rt_g00 = 0.5*(fabs(ql[6]) + fabs(qr[6]));
  double rfac_l = fabs(ql[7]*ql[8]), rfac_r = fabs(ql[7]*ql[8]);
  double rfac = 0.5*(rfac_l + rfac_r), rfac1 = 1.0 / rfac;
    
  // compute projections of jump
  double a1 = 0.5*(delta[1]-delta[5]*c*rfac1);
  double a2 = 0.5*(delta[2]+delta[4]*c*rfac);
  double a3 = 0.5*(delta[1]+delta[5]*c*rfac1);
  double a4 = 0.5*(delta[2]-delta[4]*c*rfac);
  double a5 = delta[0];
  double a6 = delta[3];

  // set waves to 0.0 as most entries vanish
  for (int i=0; i<9*3; ++i) waves[i] = 0.0;

  double *w = 0;

  // wave 1: waves with EV -c, -c
  w = &waves[0*9];
  w[1] = a1;
  w[2] = a2;
  w[4] = a2*c1*rfac1;
  w[5] = -a1*c1*rfac;
  s[0] = -c * rfac;

  // wave 2: waves with EV c, c
  w = &waves[1*9];
  w[1] = a3;
  w[2] = a4;
  w[4] = -a4*c1*rfac1;
  w[5] = a3*c1*rfac;
  s[1] = c * rfac;

  // wave 3: waves with EV 0, 0
  w = &waves[2*9];
  w[0] = a5;
  w[3] = a6;
  s[2] = 0.0;
  
  return c * rfac;
}

static void
qfluct_lax(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  const double *w0 = &waves[0*9], *w1 = &waves[1*9];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i=0; i<6; ++i)
  {
    amdq[i] = s0m*w0[i] + s1m*w1[i];
    apdq[i] = s0p*w0[i] + s1p*w1[i];
  }
  for (int i=6; i<9; ++i)
  {
    amdq[i] = 0.0;
    apdq[i] = 0.0;
  }
}

static void
qfluct(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  const double *w0 = &waves[0*9], *w1 = &waves[1*9], *w2 = &waves[2*9];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]), s2m = fmin(0.0, s[2]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]), s2p = fmax(0.0, s[2]);

  for (int i=0; i<6; ++i)
  {
    amdq[i] = s0m*w0[i] + s1m*w1[i] + s2m*w2[i];
    apdq[i] = s0p*w0[i] + s1p*w1[i] + s2p*w2[i];
  }
  for (int i=6; i<9; ++i)
  {
    amdq[i] = 0.0;
    apdq[i] = 0.0;
  }
}

static bool
check_inv(const struct gkyl_wv_eqn *eqn, const double *q)
{
  return true; // TODO: Double-check what this does
}

static double
max_speed(const struct gkyl_wv_eqn *eqn, const double *q)
{
  const struct wv_maxwell_cyl *maxwell = container_of(eqn, struct wv_maxwell_cyl, eqn);
  return gkyl_maxwell_cyl_max_abs_speed(maxwell->c, q);
}

static inline void
maxwell_cons_to_diag(const struct gkyl_wv_eqn *eqn,
  const double *qin, double *diag)
{
  // Radius, which our field components are weighted by
  double r = fabs(qin[6]*qin[7]*qin[8]);
  double r2 = r*r;
  
  // components of EM energy
  for (int i=0; i<6; ++i) diag[i] = qin[i]*qin[i] / r2;
  
  // The azimuthal components already had implicit 1/r factor due to the basis vector
  diag[1] *= r2; diag[4] *= r2;
}

struct gkyl_wv_eqn*
gkyl_wv_maxwell_cyl_new(double c, enum gkyl_wv_maxwell_cyl_rp rp_type)
{
  struct wv_maxwell_cyl *maxwell = gkyl_malloc(sizeof(struct wv_maxwell_cyl));

  maxwell->eqn.type = GKYL_EQN_MAXWELL_CYL;
  maxwell->eqn.num_equations = 9;  
  maxwell->eqn.num_diag = 6; // Er^2, Ephi^2, Ez^2, Br^2, Bphi^2, Bz^2
  
  maxwell->c = c;
  
  if (rp_type == WV_MAXWELL_CYL_RP_LAX)
  {
    maxwell->eqn.num_waves = 2;
    maxwell->eqn.waves_func = wave_lax;
    maxwell->eqn.qfluct_func = qfluct_lax;
  }
  else if (rp_type == WV_MAXWELL_CYL_RP_EIG)
  {
    maxwell->eqn.num_waves = 3;
    maxwell->eqn.waves_func = wave;
    maxwell->eqn.qfluct_func = qfluct;
  }

  maxwell->eqn.flux_jump = flux_jump;
  maxwell->eqn.check_inv_func = check_inv;
  maxwell->eqn.max_speed_func = max_speed;
  maxwell->eqn.rotate_to_local_func = rot_to_local;
  maxwell->eqn.rotate_to_global_func = rot_to_global;

  maxwell->eqn.cons_to_riem = cons_to_riem;
  maxwell->eqn.riem_to_cons = riem_to_cons;

  maxwell->eqn.wall_bc_func = maxwell_wall;

  maxwell->eqn.cons_to_diag = maxwell_cons_to_diag;

  maxwell->eqn.ref_count = gkyl_ref_count_init(maxwell_free);

  return &maxwell->eqn;
}
