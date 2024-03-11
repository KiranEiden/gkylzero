// 2D Sod-type vacuum shock tube test in axial symmetry, in polar coordinates, for the 5-moment (Euler) equations.
// Input parameters are an axisymmetric modification of those in Section 2.6.2, with the right-hand fluid set close to vacuum and the contact discontinuity placed at r = 0.75, from the thesis:
// A. Hakim (2006), "High Resolution Wave Propagation Schemes for Two-Fluid Plasma Simulations",
// PhD Thesis, University of Washington.
// https://www.aa.washington.edu/sites/aa/files/research/cpdlab/docs/PhDthesis_hakim.pdf

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_euler.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

#include <rt_arg_parse.h>

struct axi_vac_riem_ctx
{
  // Mathematical constants (dimensionless).
  double pi;

  // Physical constants (using normalized code units).
  double gas_gamma; // Adiabatic index.

  double rhol; // Left/inner fluid mass density.
  double ul; // Left/inner fluid velocity.
  double pl; // Left/inner fluid pressure.

  double rhor; // Right/outer fluid mass density.
  double ur; // Right/outer fluid velocity.
  double pr; // Right/outer fluid pressure.

  // Simulation parameters.
  int Nr; // Cell count (radial direction).
  int Ntheta; // Cell count (angular direction).
  double Lr; // Domain size (radial direction).
  double Ltheta; // Domain size (angular direction).
  double cfl_frac; // CFL coefficient.
  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.

  double rloc; // Fluid boundary (radial coordinate).
};

struct axi_vac_riem_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using normalized code units).
  double gas_gamma = 1.4; // Adiabatic index.

  double rhol = 3.0; // Left/inner fluid mass density.
  double ul = 0.0; // Left/inner fluid velocity.
  double pl = 3.0; // Left/inner fluid pressure.

  double rhor = 1.0e-6; // Right/outer fluid mass density.
  double ur = 0.0; // Right/outer fluid velocity.
  double pr = 1.0e-10; // Right/outer fluid pressure.

  // Simulation parameters.
  int Nr = 64; // Cell count (radial direction).
  int Ntheta = 64 * 6; // Cell count (angular direction).
  double Lr = 1.0; // Domain size (radial direction).
  double Ltheta = 2.0 * pi; // Domain size (angular direction).
  double cfl_frac = 0.9; // CFL coefficient.
  double t_end = 0.1; // Final simulation time.
  int num_frames = 1; // Number of output frames.

  double rloc = 0.5 * (0.25 + 1.25); // Fluid boundary (radial coordinate).

  struct axi_vac_riem_ctx ctx = {
    .pi = pi,
    .gas_gamma = gas_gamma,
    .rhol = rhol,
    .ul = ul,
    .pl = pl,
    .rhor = rhor,
    .ur = ur,
    .pr = pr,
    .Nr = Nr,
    .Ntheta = Ntheta,
    .Lr = Lr,
    .Ltheta = Ltheta,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
    .rloc = rloc,
    .num_frames = num_frames,
  };

  return ctx;
}

void
evalEulerInit(double t, const double* GKYL_RESTRICT zc, double* GKYL_RESTRICT fout, void* ctx)
{
  double r = zc[0];
  struct axi_vac_riem_ctx *app = ctx;

  double gas_gamma = app -> gas_gamma;

  double rhol = app -> rhol;
  double ul = app -> ul;
  double pl = app -> pl;

  double rhor = app -> rhor;
  double ur = app -> ur;
  double pr = app -> pr;

  double rloc = app -> rloc;

  double rho = 0.0;
  double u = 0.0;
  double p = 0.0;

  if (r < rloc) {
    rho = rhol; // Fluid mass density (left/inner).
    u = ul; // Fluid velocity (left/inner).
    p = pl; // Fluid pressure (left/inner).
  }
  else {
    rho = rhor; // Fluid mass density (right/outer).
    u = ur; // Fluid velocity (right/outer).
    p = pr; // Fluid pressure (right/outer).
  }
  
  // Set fluid mass density.
  fout[0] = rho;
  // Set fluid momentum density.
  fout[1] = rho * u; fout[2] = 0.0; fout[3] = 0.0;
  // Set fluid total energy density.
  fout[4] = p / (gas_gamma - 1.0) + 0.5 * rho * u * u;
}

static inline void
mapc2p(double t, const double* GKYL_RESTRICT zc, double* GKYL_RESTRICT xp, void* ctx)
{
  double r = zc[0], theta = zc[1];

  // Set physical coordinates (x, y) from computational coordinates (r, theta).
  xp[0] = r * cos(theta);
  xp[1] = r * sin(theta);
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_moment_app* app, double t_curr)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr)) {
    gkyl_moment_app_write(app, t_curr, iot -> curr - 1);
  }
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Init(&argc, &argv);
  }
#endif

  if (app_args.trace_mem)  {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct axi_vac_riem_ctx ctx = create_ctx(); // Context for initialization functions.

  int Nr = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nr);
  int Ntheta = APP_ARGS_CHOOSE(app_args.xcells[1], ctx.Ntheta);

  // Fluid equations.
  struct gkyl_wv_eqn *euler = gkyl_wv_euler_new(ctx.gas_gamma, app_args.use_gpu);

  struct gkyl_moment_species fluid = {
    .name = "euler",
    .equation = euler,
    .evolve = true,
    .init = evalEulerInit,
    .ctx = &ctx,

    .bcx = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },
  };

  int nrank = 1; // Number of processes in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif

  // Create global range.
  int cells[] = { Nr, Ntheta };
  int dim = sizeof(cells) / sizeof(cells[0]);
  struct gkyl_range global_r;
  gkyl_create_global_range(dim, cells, &global_r);

  // Create decomposition.
  int cuts[dim];
#ifdef GKYL_HAVE_MPI
  for (int d = 0; d < dim; d++) {
    if (app_args.use_mpi) {
      cuts[d] = app_args.cuts[d];
    }
    else {
      cuts[d] = 1;
    }
  }
#else
  for (int d = 0; d < dim; d++) {
    cuts[d] = 1;
  }
#endif

  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(dim, cuts, &global_r);

  // Construct communicator for use in app.
  struct gkyl_comm *comm;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
        .decomp = decomp
      }
    );
  }
  else {
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
        .decomp = decomp,
        .use_gpu = app_args.use_gpu
      }
    );
  }
#else
  comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .decomp = decomp,
      .use_gpu = app_args.use_gpu
    }
  );
#endif

  int my_rank;
  gkyl_comm_get_rank(comm, &my_rank);
  int comm_size;
  gkyl_comm_get_size(comm, &comm_size);

  int ncuts = 1;
  for (int d = 0; d < dim; d++) {
    ncuts *= cuts[d];
  }

  if (ncuts != comm_size) {
    if (my_rank == 0) {
      fprintf(stderr, "*** Number of ranks, %d, does not match total cuts, %d!\n", comm_size, ncuts);
    }
    goto mpifinalize;
  }

  // Moment app.
  struct gkyl_moment app_inp = {
    .name = "euler_axi_vac_riem",

    .ndim = 2,
    .lower = { 0.25, 0.0 },
    .upper = { 0.25 + ctx.Lr, 0.0 + ctx.Ltheta },
    .cells = { Nr, Ntheta },

    .mapc2p = mapc2p,

    .num_periodic_dir = 1,
    .periodic_dirs = { 1 },

    .cfl_frac = ctx.cfl_frac,

    .num_species = 1,
    .species = { fluid },

    .has_low_inp = true,
    .low_inp = {
      .local_range = decomp -> ranges[my_rank],
      .comm = comm
    }
  };

  // Create app object.
  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  // Initial and final simulation times.
  double t_curr = 0.0, t_end = ctx.t_end;

  // Create trigger for IO.
  int num_frames = ctx.num_frames;
  struct gkyl_tm_trigger io_trig = { .dt = t_end / num_frames };

  // Initialize simulation.
  gkyl_moment_app_apply_ic(app, t_curr);
  write_data(&io_trig, app, t_curr);

  // Compute estimate of maximum stable time-step.
  double dt = gkyl_moment_app_max_dt(app);

  long step = 1;
  while ((t_curr < t_end) && (step <= app_args.num_steps)) {
    gkyl_moment_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    struct gkyl_update_status status = gkyl_moment_update(app, dt);
    gkyl_moment_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      gkyl_moment_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    write_data(&io_trig, app, t_curr);

    step += 1;
  }

  write_data(&io_trig, app, t_curr);
  gkyl_moment_app_stat_write(app);

  struct gkyl_moment_stat stat = gkyl_moment_app_stat(app);

  gkyl_moment_app_cout(app, stdout, "\n");
  gkyl_moment_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_moment_app_cout(app, stdout, "Number of failed time-steps %ld\n", stat.nfail);
  gkyl_moment_app_cout(app, stdout, "Species updates took %g secs\n", stat.species_tm);
  gkyl_moment_app_cout(app, stdout, "Field updates took %g secs\n", stat.field_tm);
  gkyl_moment_app_cout(app, stdout, "Source updates took %g secs\n", stat.sources_tm);
  gkyl_moment_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

  // Free resources after simulation completion.
  gkyl_wv_eqn_release(euler);
  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
  gkyl_moment_app_release(app);  
  
mpifinalize:
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif
  
  return 0;
}
