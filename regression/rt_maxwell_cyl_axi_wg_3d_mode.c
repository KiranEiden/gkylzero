#include <math.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_maxwell_cyl.h>
#include <rt_arg_parse.h>

// map (r,theta) -> (x,y)
void
mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  double r = xc[0], ph = xc[1], z = xc[2];
  xp[0] = r*cos(ph);
  xp[1] = r*sin(ph);
  xp[2] = z;
}

void
evalFieldInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double r = xn[0], phi = xn[1], z = xn[2];
  int m = 4, n = 1;  
  double kn = 2*M_PI*n/5.0;

  double w = 2.73598725136604; // frequency of mode
  double wkn = sqrt(w*w-kn*kn);
  double a = 1.0, b = -a*jn(m,2.0*wkn)/yn(m,2.0*wkn);

  double Ez_r = a*jn(m,r*wkn) + b*yn(m,r*wkn);
  double Ez = Ez_r*cos(m*phi)*cos(kn*z);

  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = Ez;
  fout[3] = 0.0;
  fout[4] = 0.0;
  fout[5] = 0.0;
  fout[6] = 1.0;
  fout[7] = r;
  fout[8] = 1.0;
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 48);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], 96);
  int NZ = APP_ARGS_CHOOSE(app_args.xcells[2], 96);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
  
  struct gkyl_wv_eqn *maxwell_cyl = gkyl_wv_maxwell_cyl_new(299792458.0);

  // VM app
  struct gkyl_moment app_inp = {
    .name = "maxwell_cyl_axi_wg_3d_mode",

    .ndim = 3,
    // grid in computational space
    .lower = { 2.0, 0.0, 0.0 },
    .upper = { 5.0, 2*M_PI, 5.0 },
    .cells = { NX, NY, NZ },

    .mapc2p = mapc2p, // mapping of computational to physical space

    .num_periodic_dir = 2,
    .periodic_dirs = { 2, 3 },

    .cfl_frac = 0.9,

    .field = {
      .epsilon0 = 1.0, .mu0 = 1.0,
      .evolve = 1,
      .maxwell = maxwell_cyl,
      .limiter = GKYL_NO_LIMITER,
      .init = evalFieldInit,

      .bcx = { GKYL_FIELD_PEC_WALL, GKYL_FIELD_PEC_WALL },
    }
  };

  // create app object
  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  double w = 2.73598725136604; // frequency of mode
  double tperiod = 2*M_PI/w;

  // start, end and initial time-step
  double tcurr = 0.0, tend = 2*tperiod;

  // initialize simulation
  gkyl_moment_app_apply_ic(app, tcurr);
  gkyl_moment_app_write(app, tcurr, 0);

  // compute estimate of maximum stable time-step
  double dt = gkyl_moment_app_max_dt(app);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    printf("Taking time-step %ld at t = %g ...", step, tcurr);
    struct gkyl_update_status status = gkyl_moment_update(app, dt);
    printf(" dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      printf("** Update method failed! Aborting simulation ....\n");
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;

    step += 1;
  }

  gkyl_moment_app_write(app, tcurr, 1);
  gkyl_moment_app_stat_write(app);

  struct gkyl_moment_stat stat = gkyl_moment_app_stat(app);

  // simulation complete, free resources
  gkyl_moment_app_release(app);

  printf("\n");
  printf("Number of update calls %ld\n", stat.nup);
  printf("Number of failed time-steps %ld\n", stat.nfail);
  printf("Species updates took %g secs\n", stat.species_tm);
  printf("Field updates took %g secs\n", stat.field_tm);
  printf("Total updates took %g secs\n", stat.total_tm);
  
  return 0;
}
