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

struct rt_ctx
{
  // Mathematical constants (dimensionless).
  double pi;

  // Physical constants (using normalized code units).
  double gas_gamma; // Adiabatic index.
  double epsilon0; // Permittivity of free space.
  double mu0; // Permeability of free space.
  double mass_ion; // Proton mass.
  double charge_ion; // Proton charge.
  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.
  double atwood; // Fluid atwood number.
  double g_hat; // Gravitational field strength.
  double B0; // Reference magnetic field strength.
  double beta; // Ion plasma beta.
  double Te_over_Ti; // Electron temperature / ion temperature.
  double pert_max; // Maximum magnitude of initial perturbation.

  // Derived physical quantities (using normalized code units).
  double omega_ci; // Ion cyclotron frequency.
  double nl; // Left number density.
  double nr; // Right number density.
  double Ti; // Ion temperature.
  double Te; // Electron temperature.
  double T; // Total temperature.
  double vAi; // Ion Alfven velocity.
  double grav; // Gravitational acceleration.

  // Simulation parameters.
  double Nx; // Cell count (x-direction).
  double Ny; // Cell count (y-direction).
  double Lx; // Domain size (x-direction).
  double Ly; // Domain size (y-direction).
  double t_end; // Final simulation time.

  double xloc; // Fluid boundary (x-coordinate).
};

struct rt_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = 3.141592653589793238462643383279502884;

  // Physical constants (using normalized code units).
  double gas_gamma = 5.0 / 3.0; // Adiabatic index.
  double epsilon0 = 1.0; // Permittivity of free space.
  double mu0 = 1.0; // Permeability of free space.
  double mass_ion = 25.0; // Proton mass.
  double charge_ion = 1.0; // Proton charge.
  double mass_elc = 1.0; // Electron mass.
  double charge_elc = -1.0; // Electron charge.
  double atwood = 0.17; // Fluid atwood number.
  double g_hat = 0.09; // Gravitational field strength.
  double B0 = 0.1; // Reference magnetic field strength.
  double beta = 0.071; // Ion plasma beta.
  double Te_over_Ti = 0.1; // Electron temperature / ion temperature.
  double pert_max = 0.01; // Maximum magnitude of initial perturbation.

  // Derived physical quantities (using normalized code units).
  double omega_ci = charge_ion * B0 / mass_ion; // Ion cyclotron frequency.
  double nl = 1.0; // Left number density.
  double nr = nl * (1.0 - atwood) / (1.0 + atwood); // Right number density.
  double Ti = (beta / nl) * B0 * B0 / (2.0 * mu0); // Ion temperature.
  double Te = Te_over_Ti * Ti; // Electron temperature.
  double T = Te + Ti; // Total temperature.
  double vAi = B0 / sqrt(mu0 * nl * mass_ion); // Ion Alfven velocity.
  double grav = g_hat * omega_ci * vAi; // Gravitational acceleration.

  // Simulation parameters.
  double Nx = 64; // Cell count (x-direction).
  double Ny = 64; // Cell count (y-direction).
  double Lx = 3.0; // Domain size (x-direction).
  double Ly = 3.75; // Domain size (y-direction).
  double t_end = 250.0; // Final simulation time.

  double xloc = 0.5 * Lx; // Fluid boundary (x-coordinate).
  
  struct rt_ctx ctx = {
    .pi = pi,
    .gas_gamma = gas_gamma,
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .mass_ion = mass_ion,
    .charge_ion = charge_ion,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .atwood = atwood,
    .g_hat = g_hat,
    .B0 = B0,
    .beta = beta,
    .Te_over_Ti = Te_over_Ti,
    .pert_max = pert_max,
    .omega_ci = omega_ci,
    .nl = nl,
    .nr = nr,
    .Ti = Ti,
    .Te = Te,
    .T = T,
    .vAi = vAi,
    .grav = grav,
    .Nx = Nx,
    .Ny = Ny,
    .Lx = Lx,
    .Ly = Ly,
    .t_end = t_end,
    .xloc = xloc
  };

  return ctx;
}

void
evalElcInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1];
  struct rt_ctx *app = ctx;

  double pi = app -> pi;

  double gas_gamma = app -> gas_gamma;
  double mass_elc = app -> mass_elc;
  double pert_max = app -> pert_max;

  double nl = app -> nl;
  double nr = app -> nr;
  double Te = app -> Te;
  double vAi = app -> vAi;

  double Lx = app -> Lx;
  double Ly = app -> Ly;

  double xloc = app -> xloc;

  double n = 0.0;

  if (x < xloc)
  {
    n = nl; // Electron number density (left).
  }
  else
  {
    n = nr; // Electron number density (right).
  }

  pcg64_random_t rng = gkyl_pcg64_init(0); // Random number generator;

  double rhoe = n * mass_elc; // Electron mass density.
  double momxe = 0.0;
  double kx = 2.0 * pi / Lx; // Wave number (x-direction).
  double ky = 2.0 * pi / Ly; // Wave number (y-direction).
  for (int i = 0; i < 32; i++)
  {
    for (int j = 0; j < 32; j++)
    {
      momxe += rhoe * pert_max * vAi * gkyl_pcg64_rand_double(&rng) * sin(i * kx * x + 2.0 * pi * gkyl_pcg64_rand_double(&rng))
        * sin(j * ky * y + 2.0 * pi * gkyl_pcg64_rand_double(&rng)); // Electron momentum density (x-direction).
    }
  }
  double Ee_tot = n * Te / (gas_gamma - 1.0) + 0.5 * momxe * momxe / rhoe; // Electron total energy density.

  // Set electron mass density.
  fout[0] = rhoe;
  // Set electron momentum density.
  fout[1] = momxe; fout[2] = 0.0; fout[3] = 0.0;
  // Set electron total energy density.
  fout[4] = Ee_tot;
}

void
evalIonInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1];
  struct rt_ctx *app = ctx;

  double pi = app -> pi;

  double gas_gamma = app -> gas_gamma;
  double mass_ion = app -> mass_ion;
  double pert_max = app -> pert_max;

  double nl = app -> nl;
  double nr = app -> nr;
  double Ti = app -> Ti;
  double vAi = app -> vAi;

  double Lx = app -> Lx;
  double Ly = app -> Ly;

  double xloc = app -> xloc;

  double n = 0.0;

  if (x < xloc)
  {
    n = nl; // Ion number density (left).
  }
  else
  {
    n = nr; // Ion number density (right).
  }

  pcg64_random_t rng = gkyl_pcg64_init(0); // Random number generator.

  double rhoi = n * mass_ion; // Ion mass density.
  double momxi = 0.0;
  double kx = 2.0 * pi / Lx; // Wave number (x-direction).
  double ky = 2.0 * pi / Ly; // Wave number (y-direction).
  for (int i = 0; i < 32; i++)
  {
    for (int j = 0; j < 32; j++)
    {
      momxi += rhoi * pert_max * vAi * gkyl_pcg64_rand_double(&rng) * sin(i * kx * x + 2.0 * pi * gkyl_pcg64_rand_double(&rng))
        * sin(j * ky * y + 2.0 * pi * gkyl_pcg64_rand_double(&rng)); // Ion momentum density (x-direction).
    }
  }
  double Ei_tot = n * Ti / (gas_gamma - 1.0) + 0.5 * momxi * momxi / rhoi; // Ion total energy density.

  // Set ion mass density.
  fout[0] = rhoi;
  // Set ion momentum density.
  fout[1] = momxi; fout[2] = 0.0; fout[3] = 0.0;
  // Set ion total energy density.
  fout[4] = Ei_tot;
}

void
evalFieldInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1];
  struct rt_ctx *app = ctx;

  double mu0 = app -> mu0;
  double mass_ion = app -> mass_ion;
  double B0 = app -> B0;

  double nl = app -> nl;
  double nr = app -> nr;
  double T = app -> T;
  double grav = app -> grav;

  double xloc = app -> xloc;

  double Bz = 0.0;

  if (x < xloc)
  { 
    Bz = sqrt(B0 * B0 + 2 * mu0 * (mass_ion * grav * nl * x)); // Total magnetic field (z- direction, left).
  }
  else
  {
    Bz = sqrt(B0 * B0 + 2 * mu0 * ((nl - nr) * T + mass_ion * grav * nl * xloc + mass_ion * grav * nr * (x - xloc))); // Total magnetic field (z-direction, left).
  }

  // Set electric field.
  fout[0] = 0.0, fout[1] = 0.0; fout[2] = 0.0;
  // Set magnetic field.
  fout[3] = 0.0, fout[4] = 0.0; fout[5] = Bz;
  // Set correction potentials.
  fout[6] = 0.0; fout[7] = 0.0;
}

void
evalAppAccel(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct rt_ctx *app = ctx;

  double grav = app -> grav; // Gravitational acceleration.

  // Set applied acceleration.
  fout[0] = grav; fout[1] = 0.0; fout[2] = 0.0;
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
  {
    MPI_Init(&argc, &argv);
  }
#endif

  if (app_args.trace_mem) 
  {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct rt_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], ctx.Ny);

  // Electron/ion equations.
  struct gkyl_wv_eqn *elc_euler = gkyl_wv_euler_new(ctx.gas_gamma, app_args.use_gpu);
  struct gkyl_wv_eqn *ion_euler = gkyl_wv_euler_new(ctx.gas_gamma, app_args.use_gpu);

  struct gkyl_moment_species elc = {
    .name = "elc",
    .charge = ctx.charge_elc, .mass = ctx.mass_elc,
    .equation = elc_euler,
    .evolve = 1,
    .init = evalElcInit,
    .app_accel_func = evalAppAccel,
    .ctx = &ctx,

    .bcx = { GKYL_SPECIES_REFLECT, GKYL_SPECIES_REFLECT },
  };
  struct gkyl_moment_species ion = {
    .name = "ion",
    .charge = ctx.charge_ion, .mass = ctx.mass_ion,
    .equation = ion_euler,
    .evolve = 1,
    .init = evalIonInit,
    .app_accel_func = evalAppAccel,
    .ctx = &ctx,

    .bcx = { GKYL_SPECIES_REFLECT, GKYL_SPECIES_REFLECT },    
  };

  int nrank = 1; // Number of processes in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
  {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif

  // Create global range.
  int cells[] = { NX, NY };
  struct gkyl_range globalr;
  gkyl_create_global_range(2, cells, &globalr);

  // Create decomposition.
  int cuts[] = { 1, 1 };
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
  {
    cuts[0] = app_args.cuts[0];
    cuts[1] = app_args.cuts[1];
  }
#endif

  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(2, cuts, &globalr);

  // Construct communicator for use in app.
  struct gkyl_comm *comm;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
  {
    comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp)
      {
        .mpi_comm = MPI_COMM_WORLD,
        .decomp = decomp
      }
    );
  }
  else
  {
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp)
      {
        .decomp = decomp,
        .use_gpu = app_args.use_gpu
      }
    );
  }
#else
  comm = gklyl_null_comm_inew( &(struct gkyl_null_comp_inp)
    {
      .decomp = decomp,
      .use_gpu = app_args.use_gpu
    }
  );
#endif

  int my_rank;
  gkyl_comm_get_rank(comm, &my_rank);
  int comm_size;
  gkyl_comm_get_size(comm, &comm_size);

  int ncuts = cuts[0] * cuts[1];
  if (ncuts != comm_size)
  {
    if (my_rank == 0)
    {
      fprintf(stderr, "*** Number of ranks, %d, does not match total cuts, %d!\n", comm_size, ncuts);
    }
    goto mpifinalize;
  }

  // Moment app.
  struct gkyl_moment app_inp = {
    .name = "5m_rt",

    .ndim = 2,
    .lower = { 0.0, 0.0 },
    .upper = { ctx.Lx, ctx.Ly }, 
    .cells = { NX, NY },

    .num_periodic_dir = 1,
    .periodic_dirs = { 1 },
    .cfl_frac = 1.0,

    .num_species = 2,
    .species = { elc, ion },

    .field = {
      .epsilon0 = ctx.epsilon0, .mu0 = ctx.mu0,
      .mag_error_speed_fact = 1.0,
      
      .evolve = 1,
      .init = evalFieldInit,
      .ctx = &ctx,
      
      .bcx = { GKYL_FIELD_PEC_WALL, GKYL_FIELD_PEC_WALL },
    },

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

  // Initialize simulation.
  gkyl_moment_app_apply_ic(app, t_curr);
  gkyl_moment_app_write(app, t_curr, 0);

  // Compute estimate of maximum stable time-step.
  double dt = gkyl_moment_app_max_dt(app);

  long step = 1;
  while ((t_curr < t_end) && (step <= app_args.num_steps))
  {
    gkyl_moment_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    struct gkyl_update_status status = gkyl_moment_update(app, dt);
    gkyl_moment_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    
    if (!status.success)
    {
      gkyl_moment_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    step += 1;
  }

  gkyl_moment_app_write(app, t_curr, 1);
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
  gkyl_wv_eqn_release(elc_euler);
  gkyl_wv_eqn_release(ion_euler);
  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
  gkyl_moment_app_release(app);
  
mpifinalize:
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
  {
    MPI_Finalize();
  }
#endif
  
  return 0;
}
