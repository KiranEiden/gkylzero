#include <acutest.h>
#include <gkyl_prim_mhd.h>
#include <gkyl_wv_mhd.h>

void
calcq(double gas_gamma, const double pv[8], double q[8])
{
  double rho = pv[0], u = pv[1], v = pv[2], w = pv[3], pr = pv[4];
  q[0] = rho;
  q[1] = rho*u; q[2] = rho*v; q[3] = rho*w;
  q[5] = pv[5]; q[6] = pv[6]; q[7] = pv[7];  // B field
  double pb = 0.5*(pv[5]*pv[5]+pv[6]*pv[6]+pv[7]*pv[7]); // magnetic pressure
  q[4] = pr/(gas_gamma-1) + 0.5*rho*(u*u+v*v+w*w) + pb;
}

/* Check flux function implementation */
void
test_mhd_basic()
{
  double gas_gamma = 1.4;
  struct gkyl_wv_eqn *mhd = gkyl_wv_mhd_new(gas_gamma, "none");

  TEST_CHECK( mhd->num_equations == 8 );
  TEST_CHECK( mhd->num_waves == 7 );

  double rho = 1.0, u = 0.1, v = 0.2, w = 0.3, pr = 1.5;
  double bx = 0.11, by = 0.22, bz = 0.33;
  double q[8], pv[8] = { rho, u, v, w, pr, bx, by, bz };
  calcq(gas_gamma, pv, q);
  double pb = 0.5*(bx*bx+by*by+bz*bz);
  double u_dot_b = u*bx+v*by+w*bz;
  double E = q[4];

  double fluxes[3][8] = {
    { rho*u, rho*u*u - bx*bx + pr + pb, rho*u*v - bx*by, rho*u*w - bx*bz, (E+pr+pb)*u-bx*u_dot_b, 0.0, u*by - v*bx, u*bz - w*bx },
    { rho*v, rho*v*u - bx*by, rho*v*v - by*by + pr + pb, rho*v*w - by*bz, (E+pr+pb)*v-by*u_dot_b, v*bx - u*by, 0.0, v*bz - w*by },
    { rho*w, rho*w*u - bx*bz, rho*w*v - by*bz, rho*w*w - bz*bz + pr + pb, (E+pr+pb)*w-bz*u_dot_b, w*bx - u*bz, w*by - v*bz, 0.0 },
  };

  double norm[3][3] = {
    { 1.0, 0.0, 0.0 },
    { 0.0, 1.0, 0.0 },
    { 0.0, 0.0, 1.0 }
  };

  double tau1[3][3] = {
    { 0.0, 1.0, 0.0 },
    { 1.0, 0.0, 0.0 },
    { 1.0, 0.0, 0.0 }
  };

  double tau2[3][3] = {
    { 0.0, 0.0, 1.0 },
    { 0.0, 0.0, -1.0 },
    { 0.0, 1.0, 0.0 }
  };


  TEST_CHECK ( gkyl_compare(pr, gkyl_mhd_pressure(gas_gamma, q), 1e-14) );
 
  double q_local[8], flux_local[8], flux[8];
  for (int d=1; d<2; ++d) {
    mhd->rotate_to_local_func(d, tau1[d], tau2[d], norm[d], q, q_local);
    gkyl_mhd_flux(gas_gamma, q_local, flux_local);
    mhd->rotate_to_global_func(d, tau1[d], tau2[d], norm[d], flux_local, flux);
    
    for (int m=0; m<8; ++m)
      TEST_CHECK( gkyl_compare(flux[m], fluxes[d][m], 1e-14) );
  }

  double q_l[8], q_g[8];
  for (int d=1; d<3; ++d) {
    gkyl_wv_eqn_rotate_to_local(mhd, d, tau1[d], tau2[d], norm[d], q, q_l);
    gkyl_wv_eqn_rotate_to_global(mhd, d, tau1[d], tau2[d], norm[d], q_l, q_g);

    for (int m=0; m<8; ++m) TEST_CHECK( gkyl_compare(q[m], q_g[m], 1e-14) );
  }
  
  gkyl_wv_eqn_release(mhd);
}

/* check if sum of left/right going fluctuations sum to jump in flux */
void
test_mhd_waves()
{
  double gas_gamma = 5.0/3.0;
  struct gkyl_wv_eqn *mhd = gkyl_wv_mhd_new(gas_gamma, "none");

  double ql[8], qr[8];
  double ql_local[8], qr_local[8];

  // Bx must be the same since presently the wv_mhd does not have the divB wave
  double vl[8] = { 1.0,  0.1,  0.2,  0.3,  1.5, 0.4, 0.4,  0.3};
  double vr[8] = { 1.1, 0.13, 0.25, 0.34, 1.54, 0.4, 0.44, 0.34};

  calcq(gas_gamma, vl, ql); calcq(gas_gamma, vr, qr);

  double norm[3][3] = {
    { 1.0, 0.0, 0.0 },
    { 0.0, -1.0, 0.0 },
    { 0.0, 0.0, 1.0 }
  };

  double tau1[3][3] = {
    { 0.0, 1.0, 0.0 },
    { 1.0, 0.0, 0.0 },
    { 1.0, 0.0, 0.0 }
  };

  double tau2[3][3] = {
    { 0.0, 0.0, 1.0 },
    { 0.0, 0.0, 1.0 },
    { 0.0, 1.0, 0.0 }
  };  

  for (int d=0; d<1; ++d) {
    double speeds[7], waves[7*8], waves_local[7*8];
    // rotate to local tangent-normal frame
    gkyl_wv_eqn_rotate_to_local(mhd, d, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(mhd, d, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[8];
    for (int i=0; i<8; ++i) delta[i] = qr_local[i]-ql_local[i];
    
    gkyl_wv_eqn_waves(mhd, d, delta, ql_local, qr_local, waves_local, speeds);

    // rotate waves back to global frame
    for (int mw=0; mw<7; ++mw)
      gkyl_wv_eqn_rotate_to_global(mhd, d, tau1[d], tau2[d], norm[d], &waves_local[mw*8], &waves[mw*8]);
    
    double apdq[8], amdq[8];
    gkyl_wv_eqn_qfluct(mhd, d, ql, qr, waves, speeds, amdq, apdq);
    
    // check if sum of left/right going fluctuations sum to jump in flux
    double fl_local[8], fr_local[8];
    gkyl_mhd_flux(gas_gamma, ql_local, fl_local);
    gkyl_mhd_flux(gas_gamma, qr_local, fr_local);

    double fl[8], fr[8];
    gkyl_wv_eqn_rotate_to_global(mhd, d, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(mhd, d, tau1[d], tau2[d], norm[d], fr_local, fr);
    
    for (int i=0; i<8; ++i)
      TEST_CHECK( gkyl_compare(fr[i]-fl[i], amdq[i]+apdq[i], 1e-13) );
  }
    
  gkyl_wv_eqn_release(mhd);
}

void
test_mhd_waves_2()
{
  double gas_gamma = 5.0/3.0;
  struct gkyl_wv_eqn *mhd = gkyl_wv_mhd_new(gas_gamma, "none");

  double ql[8], qr[8];
  double ql_local[8], qr_local[8];

  // Bx must be the same since presently the wv_mhd does not have the divB wave
  double vl[8] = { 1.0,  0.1,  0.2,  0.3,  1.5, 0.4, 0.4, 0.3};
  double vr[8] = { 1.1, 0.13, 0.25, 0.34, 15.0, 0.4, 0.4, -0.3};

  calcq(gas_gamma, vl, ql); calcq(gas_gamma, vr, qr);

  double norm[3][3] = {
    { 1.0, 0.0, 0.0 },
    { 0.0, -1.0, 0.0 },
    { 0.0, 0.0, 1.0 }
  };

  double tau1[3][3] = {
    { 0.0, 1.0, 0.0 },
    { 1.0, 0.0, 0.0 },
    { 1.0, 0.0, 0.0 }
  };

  double tau2[3][3] = {
    { 0.0, 0.0, 1.0 },
    { 0.0, 0.0, 1.0 },
    { 0.0, 1.0, 0.0 }
  };  

  for (int d=0; d<2; ++d) {
    double speeds[7], waves[7*8], waves_local[7*8];
    // rotate to local tangent-normal frame
    gkyl_wv_eqn_rotate_to_local(mhd, d, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(mhd, d, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[8];
    for (int i=0; i<8; ++i) delta[i] = qr_local[i]-ql_local[i];
    
    gkyl_wv_eqn_waves(mhd, d, delta, ql_local, qr_local, waves_local, speeds);

    // rotate waves back to global frame
    for (int mw=0; mw<7; ++mw)
      gkyl_wv_eqn_rotate_to_global(mhd, d, tau1[d], tau2[d], norm[d], &waves_local[mw*8], &waves[mw*8]);
    
    double apdq[8], amdq[8];
    gkyl_wv_eqn_qfluct(mhd, d, ql, qr, waves, speeds, amdq, apdq);
    
    // check if sum of left/right going fluctuations sum to jump in flux
    double fl_local[8], fr_local[8];
    gkyl_mhd_flux(gas_gamma, ql_local, fl_local);
    gkyl_mhd_flux(gas_gamma, qr_local, fr_local);

    double fl[8], fr[8];
    gkyl_wv_eqn_rotate_to_global(mhd, d, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(mhd, d, tau1[d], tau2[d], norm[d], fr_local, fr);
    
    for (int i=0; i<8; ++i)
      TEST_CHECK( gkyl_compare(fr[i]-fl[i], amdq[i]+apdq[i], 1e-13) );
  }
    
  gkyl_wv_eqn_release(mhd);
}

void
test_mhd_waves_3()
{
  double gas_gamma = 5.0/3.0;
  struct gkyl_wv_eqn *mhd = gkyl_wv_mhd_new(gas_gamma, "glm");

  double ql[9], qr[9];
  double ql_local[9], qr_local[9];

  // Bx must be the same since presently the wv_mhd does not have the divB wave
  double vl[9] = { 1.0,  0.1,  0.2,  0.3,  1.5, 0.4, 0.4, 0.3, 0.1};
  double vr[9] = { 1.1, 0.13, 0.25, 0.34, 15.0, 0.4, 0.4, -0.3, 0.1};

  calcq(gas_gamma, vl, ql); calcq(gas_gamma, vr, qr);

  double norm[3][3] = {
    { 1.0, 0.0, 0.0 },
    { 0.0, -1.0, 0.0 },
    { 0.0, 0.0, 1.0 }
  };

  double tau1[3][3] = {
    { 0.0, 1.0, 0.0 },
    { 1.0, 0.0, 0.0 },
    { 1.0, 0.0, 0.0 }
  };

  double tau2[3][3] = {
    { 0.0, 0.0, 1.0 },
    { 0.0, 0.0, 1.0 },
    { 0.0, 1.0, 0.0 }
  };  

  for (int d=0; d<2; ++d) {
    double speeds[9], waves[9*9], waves_local[9*9];
    // rotate to local tangent-normal frame
    gkyl_wv_eqn_rotate_to_local(mhd, d, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(mhd, d, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[9];
    for (int i=0; i<9; ++i) delta[i] = qr_local[i]-ql_local[i];
    
    gkyl_wv_eqn_waves(mhd, d, delta, ql_local, qr_local, waves_local, speeds);

    // rotate waves back to global frame
    for (int mw=0; mw<9; ++mw)
      gkyl_wv_eqn_rotate_to_global(mhd, d, tau1[d], tau2[d], norm[d], &waves_local[mw*9], &waves[mw*9]);
    
    double apdq[9], amdq[9];
    gkyl_wv_eqn_qfluct(mhd, d, ql, qr, waves, speeds, amdq, apdq);
    
    // check if sum of left/right going fluctuations sum to jump in flux
    double fl_local[9], fr_local[9];
    gkyl_mhd_flux(gas_gamma, ql_local, fl_local);
    gkyl_mhd_flux(gas_gamma, qr_local, fr_local);

    double fl[9], fr[9];
    gkyl_wv_eqn_rotate_to_global(mhd, d, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(mhd, d, tau1[d], tau2[d], norm[d], fr_local, fr);
    
    for (int i=0; i<9; ++i)
      TEST_CHECK( gkyl_compare(fr[i]-fl[i], amdq[i]+apdq[i], 1e-13) );
  }
    
  gkyl_wv_eqn_release(mhd);
}

TEST_LIST = {
  { "mhd_basic", test_mhd_basic },
  { "mhd_waves", test_mhd_waves },
  { "mhd_waves_2", test_mhd_waves_2 },
  { "mhd_waves_3", test_mhd_waves_3 },
  { NULL, NULL },
};
