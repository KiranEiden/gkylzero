#include <math.h>
#include <stdio.h>

#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_iso_euler.h>
#include <app_arg_parse.h>

struct iso_euler_ctx {
  double cs; // sound speed
};

void
evalEulerInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0];

  double rhol = 3.0, ul = 0.5;
  double rhor = 1.0, ur = 0.0;

  double rho = rhor, u = ur;
  if (x<0.5) {
    rho = rhol;
    u = ul;
  }

  fout[0] = rho;
  fout[1] = rho*u; fout[2] = 0.0; fout[3] = 0.0;
}

struct iso_euler_ctx
iso_euler_ctx(void)
{
  return (struct iso_euler_ctx) { .cs = 1.0 };
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = get_parse_app_args(argc, argv);  
  struct iso_euler_ctx ctx = iso_euler_ctx(); // context for init functions

  // equation object
  struct gkyl_wv_eqn *iso_euler = gkyl_wv_iso_euler_new(ctx.cs);

  struct gkyl_moment_species fluid = {
    .name = "euler",

    .equation = iso_euler,
    .evolve = 1,
    .ctx = &ctx,
    .init = evalEulerInit,

    .bcx = { GKYL_MOMENT_COPY, GKYL_MOMENT_COPY },
  };

  // VM app
  struct gkyl_moment app_inp = {
    .name = "iso_euler_sodshock",

    .ndim = 1,
    .lower = { 0.0 },
    .upper = { 1.0 }, 
    .cells = { 512 },

    .cfl_frac = 0.9,

    .num_species = 1,
    .species = { fluid },
  };

  // create app object
  gkyl_moment_app *app = gkyl_moment_app_new(app_inp);

  // start, end
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

  struct gkyl_moment_stat stat = gkyl_moment_app_stat(app);

  // simulation complete, free resources
  gkyl_wv_eqn_release(iso_euler);
  gkyl_moment_app_release(app);

  printf("\n");
  printf("Number of update calls %ld\n", stat.nup);
  printf("Number of failed time-steps %ld\n", stat.nfail);
  printf("Species updates took %g secs\n", stat.species_tm);
  printf("Field updates took %g secs\n", stat.field_tm);
  printf("Total updates took %g secs\n", stat.total_tm);
  
  return 0;
}
