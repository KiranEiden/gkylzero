#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>
#include <rt_arg_parse.h>

struct pkpm_sheath_ctx {
  double chargeElc; // electron charge
  double massElc; // electron mass
  double chargeIon; // ion charge
  double massIon; // ion mass
  double Te_Ti; // electron to ion temperature ratio
  double vte; // electron thermal velocity
  double vti; // ion thermal velocity
  double Lx; // size of the box
  double n0; // initial number density
  double tend;
  double min_dt;
};

static inline double sq(double x) { return x*x; }

void
evalDistFuncElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_sheath_ctx *app = ctx;
  double x = xn[0], v = xn[1];
  double vt = app->vte, n0 = app->n0;
  fout[0] = n0/sqrt(2.0*M_PI*sq(vt))*(exp(-sq(v)/(2*sq(vt))));
  fout[1] = vt*vt*n0/sqrt(2.0*M_PI*sq(vt))*(exp(-sq(v)/(2*sq(vt))));
}

void
evalDistFuncIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_sheath_ctx *app = ctx;
  double x = xn[0], v = xn[1];
  double vt = app->vti, n0 = app->n0;
  fout[0] = n0/sqrt(2.0*M_PI*sq(vt))*(exp(-sq(v)/(2*sq(vt))));
  fout[1] = vt*vt*n0/sqrt(2.0*M_PI*sq(vt))*(exp(-sq(v)/(2*sq(vt))));
}

void
evalFluidElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_sheath_ctx *app = ctx;
  double x = xn[0]; 
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
}

void
evalFluidIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_sheath_ctx *app = ctx;
  double x = xn[0]; 
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
}

void
evalFieldFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_sheath_ctx *app = ctx;
  double x = xn[0];
  
  fout[0] = 0.0; fout[1] = 0.0, fout[2] = 0.0;
  fout[3] = 0.0; fout[4] = 0.0; fout[5] = 0.0;
  fout[6] = 0.0; fout[7] = 0.0;
}

void
evalExtEmFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_sheath_ctx *app = ctx;
  double x = xn[0];
  double B_x = 1.0;
  
  fout[0] = 0.0; fout[1] = 0.0, fout[2] = 0.0;
  fout[3] = B_x; fout[4] = 0.0; fout[5] = 0.0;
}

void
evalNuElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_sheath_ctx *app = ctx;
  fout[0] = 1.0e-5;
}

void
evalNuIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_sheath_ctx *app = ctx;
  fout[0] = 1.0e-5/sqrt(app->massIon)*(app->Te_Ti*sqrt(app->Te_Ti));
}

struct pkpm_sheath_ctx
create_ctx(void)
{
  struct pkpm_sheath_ctx ctx = {
    .chargeElc = -1.0,
    .massElc = 1.0,
    .chargeIon = 1.0,
    .massIon = 1836.153,
    .Te_Ti = 1.0,
    .vte = 1.0,
    .vti = ctx.vte/sqrt(ctx.Te_Ti*ctx.massIon),
    .Lx = 128.0,
    .n0 = 1.0,
    .tend = 100.0,
    .min_dt = 1.0e-4, 
  };
  return ctx;
}

void
write_data(struct gkyl_tm_trigger *iot, gkyl_vlasov_app *app, double tcurr)
{
  if (gkyl_tm_trigger_check_and_bump(iot, tcurr)) {
    gkyl_vlasov_app_write(app, tcurr, iot->curr-1);
    gkyl_vlasov_app_calc_mom(app); gkyl_vlasov_app_write_mom(app, tcurr, iot->curr-1);
  }
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 128);
  int VX = APP_ARGS_CHOOSE(app_args.vcells[0], 32);  

  struct pkpm_sheath_ctx ctx = create_ctx(); // context for init functions

  // electron momentum
  struct gkyl_vlasov_fluid_species fluid_elc = {
    .name = "fluid_elc",
    .num_eqn = 3,
    .pkpm_species = "elc",
    .ctx = &ctx,
    .init = evalFluidElc,
    .bcx = { GKYL_SPECIES_ABSORB, GKYL_SPECIES_ABSORB },
    //.diffusion = {.D = 1.0e-5, .order=4},
  };  
  
  // electrons
  struct gkyl_vlasov_species elc = {
    .name = "elc",
    .model_id = GKYL_MODEL_PKPM,
    .pkpm_fluid_species = "fluid_elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { -6.0 * ctx.vte},
    .upper = { 6.0 * ctx.vte}, 
    .cells = { VX },

    .ctx = &ctx,
    .init = evalDistFuncElc,

    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,

      .ctx = &ctx,
      .self_nu = evalNuElc,
      .normNu = true,
    },    

    .num_diag_moments = 0,
    .bcx = { GKYL_SPECIES_ABSORB, GKYL_SPECIES_ABSORB },
  };

  // ion momentum
  struct gkyl_vlasov_fluid_species fluid_ion = {
    .name = "fluid_ion",
    .num_eqn = 3,
    .pkpm_species = "ion",
    .ctx = &ctx,
    .init = evalFluidIon,
    .bcx = { GKYL_SPECIES_ABSORB, GKYL_SPECIES_ABSORB },
    //.diffusion = {.D = 1.0e-5, .order=4},
  };  
  
  // ions
  struct gkyl_vlasov_species ion = {
    .name = "ion",
    .model_id = GKYL_MODEL_PKPM,
    .pkpm_fluid_species = "fluid_ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .lower = { -6.0 * ctx.vti},
    .upper = { 6.0 * ctx.vti}, 
    .cells = { VX },

    .ctx = &ctx,
    .init = evalDistFuncIon,

    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,

      .ctx = &ctx,
      .self_nu = evalNuIon,
      .normNu = true,
    },    

    .num_diag_moments = 0,
    .bcx = { GKYL_SPECIES_ABSORB, GKYL_SPECIES_ABSORB },
  };

  // field
  struct gkyl_vlasov_field field = {
    .epsilon0 = 1.0, .mu0 = 1.0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,

    .ctx = &ctx,
    .init = evalFieldFunc,
    // Plasma EM field BCs are PEC, external field goes into conducting wall
    .bcx = { GKYL_FIELD_PEC_WALL, GKYL_FIELD_PEC_WALL }, 
    .ext_em = evalExtEmFunc,
    .ext_em_ctx = &ctx,
  };

  // VM app
  struct gkyl_vm vm = {
    .name = "pkpm_sheath_p1",

    .cdim = 1, .vdim = 1,
    .lower = { -ctx.Lx },
    .upper = { ctx.Lx },
    .cells = { NX },
    .poly_order = 1,
    .basis_type = app_args.basis_type,

    .num_periodic_dir = 0,
    .periodic_dirs = { },

    .num_species = 2,
    .species = { elc, ion },
    .num_fluid_species = 2,
    .fluid_species = { fluid_elc, fluid_ion },
    .field = field,

    .use_gpu = app_args.use_gpu,
  };

  // create app object
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&vm);

  // start, end and initial time-step
  double tcurr = 0.0, tend = ctx.tend;
  double dt = tend-tcurr;
  int nframe = 1;
  // create trigger for IO
  struct gkyl_tm_trigger io_trig = { .dt = tend/nframe };

  // initialize simulation
  gkyl_vlasov_app_apply_ic(app, tcurr);
  write_data(&io_trig, app, tcurr);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    printf("Taking time-step at t = %g ...", tcurr);
    struct gkyl_update_status status = gkyl_vlasov_update(app, dt);
    printf(" dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      printf("** Update method failed! Aborting simulation ....\n");
      break;
    }
    if (status.dt_actual < ctx.min_dt) {
      printf("** Time step crashing! Aborting simulation and writing out last output ....\n");
      gkyl_vlasov_app_write(app, tcurr, 1000);
      gkyl_vlasov_app_calc_mom(app); gkyl_vlasov_app_write_mom(app, tcurr, 1000);
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;

    write_data(&io_trig, app, tcurr);

    step += 1;
  }

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
  printf("Field RHS calc took %g secs\n", stat.field_rhs_tm);
  printf("Current evaluation and accumulate took %g secs\n", stat.current_tm);
  printf("Updates took %g secs\n", stat.total_tm);
  
  return 0;
}
