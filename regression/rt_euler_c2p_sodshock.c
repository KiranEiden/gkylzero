#include <math.h>
#include <stdio.h>

#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_euler.h>
#include <rt_arg_parse.h>

struct euler_ctx {
  double gas_gamma; // gas constant
};

static inline void
mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  double z0 = xc[0];
  xp[0] = 0.5*xc[0]+0.5;
}

void
evalEulerInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct euler_ctx *app = ctx;
  double gas_gamma = app->gas_gamma;

  double xp[GKYL_MAX_CDIM];
  mapc2p(0.0, xn, xp, 0);
  double x = xp[0]; // x-coordinate in physical space

  double rhol = 3.0, ul = 0.0, pl = 3.0;
  double rhor = 1.0, ur = 0.0, pr = 1.0;

  double rho = rhor, u = ur, p = pr;
  if (x<0.5) {
    rho = rhol;
    u = ul;
    p = pl;
  }

  fout[0] = rho;
  fout[1] = rho*u; fout[2] = 0.0; fout[3] = 0.0;
  fout[4] = p/(gas_gamma-1) + 0.5*rho*u*u;
}

struct euler_ctx
euler_ctx(void)
{
  return (struct euler_ctx) { .gas_gamma = 1.4 };
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);
  struct euler_ctx ctx = euler_ctx(); // context for init functions

  // equation object
  struct gkyl_wv_eqn *euler = gkyl_wv_euler_new(ctx.gas_gamma);

  struct gkyl_moment_species fluid = {
    .name = "euler",

    .equation = euler,
    .evolve = 1,
    .ctx = &ctx,
    .init = evalEulerInit,

    .bcx = { GKYL_MOMENT_COPY, GKYL_MOMENT_COPY },
  };

  // VM app
  struct gkyl_moment app_inp = {
    .name = "euler_c2p_sodshock",

    .ndim = 1,
    .lower = { -1.0 },
    .upper = { 1.0 }, 
    .cells = { 512 },

    .mapc2p = mapc2p, // mapping of computational to physical space    

    .cfl_frac = 0.9,

    .num_species = 1,
    .species = { fluid },
  };

  // create app object
  gkyl_moment_app *app = gkyl_moment_app_new(app_inp);

  // start, end and initial time-step
  double tcurr = 0.0, tend = 0.1;

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
  gkyl_wv_eqn_release(euler);
  gkyl_moment_app_release(app);

  printf("\n");
  printf("Number of update calls %ld\n", stat.nup);
  printf("Number of failed time-steps %ld\n", stat.nfail);
  printf("Species updates took %g secs\n", stat.species_tm);
  printf("Field updates took %g secs\n", stat.field_tm);
  printf("Total updates took %g secs\n", stat.total_tm);
  
  return 0;
}