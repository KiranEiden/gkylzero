#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>
#include <rt_arg_parse.h>

struct free_stream_ctx {
  double charge; // charge
  double mass; // mass
  double vt; // thermal velocity
  double Lx; // size of the box
};

static inline double sq(double x) { return x*x; }

static inline double
maxwellian(double n, double v, double u, double vth)
{
  double v2 = (v - u)*(v - u);
  return n/sqrt(2*M_PI*vth*vth)*exp(-v2/(2*vth*vth));
}

static inline double
maxwelljuttner1D(double n, double px, double ux, double T)
{

  // Set the normalization  
  double K1;
  if (T == 1.0){
    K1 = 0.601907230197235;
  } else {
    K1 = 0.495079105512939;
  }

  // All constants = 1 (c, m0, kb)
  double gamma = 1.0/sqrt(1.0 - ux*ux);
  return n/(2*K1)*exp(-(gamma/T)*(sqrt(1 + px*px) - ux*px ));
}

void
evalDistFunc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct free_stream_ctx *app = ctx;
  double x = xn[0], v = xn[1];
  if (x<0.5)
    fout[0] = maxwelljuttner1D(1.0, v, 0.0, 1.0); //maxwellian(1.0, v, 0.0, 1.0);
  else
    fout[0] = maxwelljuttner1D(0.125, v, 0.0, sqrt(0.1/0.125)); //maxwellian(0.125, v, 0.0, sqrt(0.1/0.125));
}

void
evalNu(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct free_stream_ctx *app = ctx;
  double x = xn[0], v = xn[1];
  fout[0] = 100.0;
}

struct free_stream_ctx
create_ctx(void)
{
  struct free_stream_ctx ctx = {
    .mass = 1.0,
    .charge = 1.0,
    .vt = 1.0,
    .Lx = 1.0,
  };
  return ctx;
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 128);
  int NV = APP_ARGS_CHOOSE(app_args.vcells[0], 32); //16

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
  struct free_stream_ctx ctx = create_ctx(); // context for init functions

  // electrons
  struct gkyl_vlasov_species neut = {
    .name = "neut",
    .model_id = GKYL_MODEL_SR,
    .charge = ctx.charge, .mass = ctx.mass,
    .lower = { -10.0*ctx.vt},
    .upper = { 10.0*ctx.vt}, 
    .cells = { NV },

    .ctx = &ctx,
    .init = evalDistFunc,

    .collisions =  {
      .collision_id = GKYL_BGK_COLLISIONS,

      .ctx = &ctx,
      .self_nu = evalNu,
    },

    //.num_diag_moments = 3,
    //.diag_moments = { "M0", "M1i", "M2" },
  };

  // VM app
  struct gkyl_vm vm = {
    .name = "neut_sr_bgk_sod_shock_1x1v_p2",

    .cdim = 1, .vdim = 1,
    .lower = { 0.0 },
    .upper = { ctx.Lx },
    .cells = { NX },
    .poly_order = 2,
    .basis_type = app_args.basis_type,

    .num_periodic_dir = 0,
    .periodic_dirs = {  },

    .num_species = 1,
    .species = { neut },
    .skip_field = true,

    .use_gpu = app_args.use_gpu,
  };

  // create app object
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&vm);

  // start, end and initial time-step
  double tcurr = 0.0, tend = 0.1;
  double dt = tend-tcurr;

  // initialize simulation
  gkyl_vlasov_app_apply_ic(app, tcurr);
  
  gkyl_vlasov_app_write(app, tcurr, 0);
  gkyl_vlasov_app_calc_mom(app); gkyl_vlasov_app_write_mom(app, tcurr, 0);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    printf("Taking time-step at t = %g ...", tcurr);
    struct gkyl_update_status status = gkyl_vlasov_update(app, dt);
    printf(" dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      printf("** Update method failed! Aborting simulation ....\n");
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;
    step += 1;
  }

  gkyl_vlasov_app_write(app, tcurr, 1);
  gkyl_vlasov_app_calc_mom(app); gkyl_vlasov_app_write_mom(app, tcurr, 1);
  gkyl_vlasov_app_stat_write(app);

  // fetch simulation statistics
  struct gkyl_vlasov_stat stat = gkyl_vlasov_app_stat(app);

  // simulation complete, free app
  gkyl_vlasov_app_release(app);

  printf("\n");
  printf("Number of update calls %ld\n", stat.nup);
  printf("Number of forward-Euler calls %ld\n", stat.nfeuler);
  printf("Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    printf("Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    printf("Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }  
  printf("Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  printf("Species RHS calc took %g secs\n", stat.species_rhs_tm);
  printf("Updates took %g secs\n", stat.total_tm);
  
  return 0;
}
