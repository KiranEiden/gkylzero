#include <acutest.h>
#include <gkyl_moment_prim_maxwell_cyl.h>
#include <gkyl_wv_maxwell_cyl.h>

// Make indexing cleaner with the dir_shuffle
#define ER 0
#define EP 1
#define EZ 2
#define BR 3
#define BP 4
#define BZ 5

void
test_maxwell_basic()
{
  // speed of light in SI units so tests are non-trivial
  double c = 299792458.0;
  double c2 = c*c;
  struct gkyl_wv_eqn *maxwell = gkyl_wv_maxwell_cyl_new(c, WV_MAXWELL_CYL_RP_LAX);

  TEST_CHECK( maxwell->num_equations == 9 );
  TEST_CHECK( maxwell->num_waves == 2 );

  double Er = 1.0, Ephi = 0.1, Ez = 0.2;
  double Br = 10.0, Bphi = 10.1, Bz = 10.2;
  double r = 0.5;
  double q[9] = { Er, Ephi, Ez, Br, Bphi, Bz, 1.0, r, 1.0};
  
  double norm[3][3] =
  {
    { 1.0, 0.0, 0.0 },
    { 0.0, 1.0, 0.0 },
    { 0.0, 0.0, 1.0 }
  };

  double tau1[3][3] =
  {
    { 0.0, 1.0, 0.0 },
    { 1.0, 0.0, 0.0 },
    { 0.0, -1.0, 0.0 }
  };

  double tau2[3][3] =
  {
    { 0.0, 0.0, 1.0 },
    { 0.0, 0.0, -1.0 },
    { 1.0, 0.0, 0.0 }
  };

  double fluxes[3][9] =
  {
    {
      0.0,
      c2/r*r*q[BZ],
      -c2*r*r*q[BP],
      
      0.0,
      -q[EZ],
      q[EP]*r*r,
      
      0.0,
      0.0,
      0.0
    },
    {
      -c2*q[BZ],
      0.0,
      c2*q[BR],
      
      q[EZ],
      0.0,
      -q[ER],
      
      0.0,
      0.0,
      0.0
    },
    {
      c2*r*r*q[BP],
      -c2*q[BR],      
      0.0,
      
      -q[EP]*r*r,
      q[ER],
      0.0,
      
      0.0,
      0.0,
      0.0
    },

  };


  double q_local[9], flux_local[9], flux[9];

  for(int d=0; d<3; ++d)
  {
    maxwell->rotate_to_local_func(tau1[d], tau2[d], norm[d], q, q_local);
    gkyl_maxwell_cyl_flux(c, q_local, flux_local);
    maxwell->rotate_to_global_func(tau1[d], tau2[d], norm[d], flux_local, flux);

    for (int m=0; m<9; ++m)
      TEST_CHECK( gkyl_compare(flux[m], fluxes[d][m], 1e-15) );

    // check Riemann transform
    double w1[9], q1[9];
    maxwell->cons_to_riem(maxwell, q_local, q_local, w1);
    maxwell->riem_to_cons(maxwell, q_local, w1, q1);
    
    for(int m=0; m<9; ++m)
      TEST_CHECK( gkyl_compare_double(q_local[m], q1[m], 1e-14) );    
  }

  gkyl_wv_eqn_release(maxwell);
}

void
test_maxwell_waves()
{
  // speed of light in SI units so tests are non-trivial
  double c = 299792458.0;
  double c2 = c*c;
  struct gkyl_wv_eqn *maxwell = gkyl_wv_maxwell_cyl_new(c, WV_MAXWELL_CYL_RP_EIG);

  double r = 0.5;
  double ql[9] = { 0.0, 1.0, 0.0, 1.0, -0.75, 0.0, 1.0, r, 1.0};
  double qr[9] = { 0.0, -1.0, 0.0, 1.0, 0.75, 0.0, 1.0, r, 1.0};

  double norm[3][3] =
  {
    { 1.0, 0.0, 0.0 },
    { 0.0, 1.0, 0.0 },
    { 0.0, 0.0, 1.0 }
  };

  double tau1[3][3] =
  {
    { 0.0, 1.0, 0.0 },
    { 1.0, 0.0, 0.0 },
    { 1.0, 0.0, 0.0 }
  };

  double tau2[3][3] =
  {
    { 0.0, 0.0, 1.0 },
    { 0.0, 0.0, -1.0 },
    { 0.0, 1.0, 0.0 }
  };  

  for(int d=0; d<3; ++d)
  {
    double speeds[3], waves[3*9], waves_local[3*9];
    double ql_local[9], qr_local[9];
    
    // rotate to local tangent-normal frame
    gkyl_wv_eqn_rotate_to_local(maxwell, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(maxwell, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[9];
    for (int i=0; i<9; ++i) delta[i] = qr_local[i]-ql_local[i];
    
    gkyl_wv_eqn_waves(maxwell, GKYL_WV_HIGH_ORDER_FLUX, delta, ql_local, qr_local, waves_local, speeds);

    double apdq_local[9], amdq_local[9];
    gkyl_wv_eqn_qfluct(maxwell, GKYL_WV_HIGH_ORDER_FLUX, ql_local, qr_local, waves_local, speeds,
      amdq_local, apdq_local);

    // rotate waves back to global frame
    for (int mw=0; mw<3; ++mw)
      gkyl_wv_eqn_rotate_to_global(maxwell, tau1[d], tau2[d], norm[d], &waves_local[mw*9], &waves[mw*9]);

    double apdq[9], amdq[9];
    gkyl_wv_eqn_rotate_to_global(maxwell, tau1[d], tau2[d], norm[d], apdq_local, apdq);
    gkyl_wv_eqn_rotate_to_global(maxwell, tau1[d], tau2[d], norm[d], amdq_local, amdq);
    
    // check if sum of left/right going fluctuations sum to jump in flux
    double fl_local[9], fr_local[9];
    gkyl_maxwell_cyl_flux(c, ql_local, fl_local);
    gkyl_maxwell_cyl_flux(c, qr_local, fr_local);
    
    double fl[9], fr[9];
    gkyl_wv_eqn_rotate_to_global(maxwell, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(maxwell, tau1[d], tau2[d], norm[d], fr_local, fr);
    
    for (int i=0; i<9; ++i) {
      // printf("%d: %g %g (%g)\n", i, fr[i]-fl[i], amdq[i]+apdq[i], (fr[i]-fl[i])-(amdq[i]+apdq[i]));
      TEST_CHECK( gkyl_compare(fr[i]-fl[i], amdq[i]+apdq[i], 1e-14) );
    }
  }
    
  gkyl_wv_eqn_release(maxwell);
}

TEST_LIST =
{
  { "maxwell_basic", test_maxwell_basic },
  { "maxwell_waves", test_maxwell_waves },
  { NULL, NULL },
};
