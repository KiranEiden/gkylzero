// Burch et al. magnetic reconnection test for the 5-moment equations.
// Input parameters match the initial conditions in Section 3, from the article:
// J. M. TenBarge, J. Ng, J. Juno, L. Wang, A. M. Hakim and A. Bhattacharjee (2019), "An Extended MHD Study of the 16 October 2015 MMS Diffusion Region Crossing",
// Journal of Geophysical Research: Space Physics, Volume 124 (11): 8474-8487.
// https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019JA026731

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

struct burch_ctx
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
  double vAe; // Electron Alfven velocity.
  double n0; // Reference number density.
  double beta2; // Magnetosheath plasma beta.
  double Ti1_over_Ti2; // Magnetospheric ion temperature / magnetosheath ion temperature.
  double Te1_over_Te2; // Magnetospheric electron temperature / magnetosheath ion temperature.

  // Derived physical quantities (using normalized code units).
  double light_speed; // Speed of light.
  double B0; // Reference magnetic field strength.
  double omega_pi; // Ion plasma frequency.
  double di; // Ion skin depth.
  double omega0; // Reference frequency.
  double psi0; // Reference magnetic scalar potential.
  double guide1; // Magnetospheric guide field strength.
  double guide2; // Magnetosheath guide field strength.
  double b1; // Magnetospheric magnetic field strength.
  double b2; // Magnetosheath magnetic field strength.
  double n1; // Magnetospheric number density.
  double n2; // Magnetosheath number density.
  double Ti2; // Magnetosheath ion temperature.

  double Te1; // Magnetospheric electron temperature.
  double Ti1; // Magnetospheric ion temperature.

  double Te2; // Magnetosheath electron temperature (so that the system is in force balance).

  // Simulation parameters.
  int Nx; // Cell count (x-direction);
  int Ny; // Cell count (y-direction);
  double Lx; // Domain size (x-direction).
  double Ly; // Domain size (y-direction).
  double cfl_frac; // CFL coefficient.
  double t_end; // Final simulation time.
};

struct burch_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using normalized code units).
  double gas_gamma = 5.0 / 3.0; // Adiabatic index.
  double epsilon0 = 1.0; // Permittivity of free space.
  double mu0 = 1.0; // Permeability of free space.
  double mass_ion = 1.0; // Proton mass.
  double charge_ion = 1.0; // Proton charge.
  double mass_elc = mass_ion / 100.0; // Electron mass.
  double charge_elc = -1.0; // Electron charge.
  double vAe = 0.2; // Electron Alfven velocity.
  double n0 = 1.0; // Reference number density.
  double beta2 = 2.748; // Magnetosheath plasma beta.
  double Ti1_over_Ti2 = 7.73 / 1.374; // Magnetospheric ion temperature / magnetosheath ion temperature.
  double Te1_over_Te2 = 1.288 / 1.374; // Magnetospheric electron temperature / magnetosheath ion temperature.

  // Derived physical quantities (using normalized code units).
  double light_speed = 1.0 / sqrt(mu0 * epsilon0); // Speed of light.
  double B0 = vAe * sqrt(n0 * mass_elc); // Reference magnetic field strength.
  double omega_pi = sqrt(n0 * charge_ion * charge_ion / (epsilon0 * mass_ion)); // Ion plasma frequency.
  double di = light_speed / omega_pi; // Ion skin depth.
  double omega0 = 1.0 * di; // Reference frequency.
  double psi0 = 0.1 * B0 * di; // Reference magnetic scalar potential.
  double guide1 = 0.099 * B0; // Magnetospheric guide field strength.
  double guide2 = guide1; // Magnetosheath guide field strength.
  double b1 = 1.696 * B0; // Magnetospheric magnetic field strength.
  double b2 = 1.00 * B0; // Magnetosheath magnetic field strength.
  double n1 = 0.062 * n0; // Magnetospheric number density.
  double n2 = 1.0 * n0; // Magnetosheath number density.
  double Ti2 = beta2 * (b2 * b2) / (2.0 * n2 * mu0); // Magnetosheath ion temperature.

  double Te1 = Ti2 * Te1_over_Te2; // Magnetospheric electron temperature.
  double Ti1 = Ti2 * Ti1_over_Ti2; // Magnetospheric ion temperature.

  double Te2 = (0.5 * (b1 * b1 - b2 * b2) + 0.5 * (guide1 * guide1 - guide2 * guide2)
    + n1 * (Ti1 + Te1) - n2 * Ti2) / n2; // Magnetosheath electron temperature (so that the system is in force balance).

  // Simulation parameters.
  double Nx = 256; // Cell count (x-direction).
  double Ny = 128; // Cell count (y-direction).
  double Lx = 40.96 * di; // Domain size (x-direction).
  double Ly = 20.48 * di; // Domain size (y-direction).
  double cfl_frac = 1.0; // CFL coefficient.
  double t_end = 250.0; // Final simulation time.
  
  struct burch_ctx ctx = {
    .pi = pi,
    .gas_gamma = gas_gamma,
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .mass_ion = mass_ion,
    .charge_ion = charge_ion,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .vAe = vAe,
    .n0 = n0,
    .beta2 = beta2,
    .Ti1_over_Ti2 = Ti1_over_Ti2,
    .Te1_over_Te2 = Te1_over_Te2,
    .light_speed = light_speed,
    .B0 = B0,
    .omega_pi = omega_pi,
    .di = di,
    .omega0 = omega0,
    .psi0 = psi0,
    .guide1 = guide1,
    .guide2 = guide2,
    .b1 = b1,
    .b2 = b2,
    .n1 = n1,
    .n2 = n2,
    .Ti2 = Ti2,
    .Te1 = Te1,
    .Ti1 = Ti1,
    .Te2 = Te2,
    .Nx = Nx,
    .Ny = Ny,
    .Lx = Lx,
    .Ly = Ly,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
  };

  return ctx;
}

void
evalElcInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1];
  struct burch_ctx *app = ctx;

  double pi = app -> pi;

  double gas_gamma = app -> gas_gamma;
  double mass_elc = app -> mass_elc;
  double charge_elc = app -> charge_elc;

  double omega0 = app -> omega0;
  double psi0 = app -> psi0;
  double guide1 = app -> guide1;
  double guide2 = app -> guide2;
  double b1 = app -> b1;
  double b2 = app -> b2;
  double n1 = app -> n1;
  double Ti2 = app -> Ti2;

  double Te1 = app -> Te1;
  double Ti1 = app -> Ti1;

  double Te2 = app -> Te2;

  double Lx = app -> Lx;
  double Ly = app -> Ly;

  double b1x = 0.5 * (b2 + b1) * (tanh((y - Ly * 0.25) / omega0) - tanh((y - Ly * 0.75) / omega0)
    + tanh((y - Ly * 1.25) / omega0) - tanh((y + Ly * 0.25) / omega0) + 1.0) + 0.5 * (b2 - b1); // Magnetospheric magnetic field (x-direction).
  double b1y = 0.0; // Magnetospheric magnetic field (y-direction).
  double b1z = 0.5 * (guide2 - guide1) * (tanh((y - Ly * 0.25) / omega0) - tanh((y - Ly * 0.75) / omega0)
    + tanh((y - Ly * 1.25) / omega0) - tanh((y + Ly * 0.25) / omega0) + 1.0) + 0.5 * (guide2 + guide1); // Magnetospheric magnetic field (z-direction).

  double Ti_tot = 0.5 * (Ti2 - Ti1) * (tanh((y - Ly * 0.25) / omega0) - tanh((y - Ly * 0.75) / omega0)
    + tanh((y - Ly * 1.25) / omega0) - tanh((y + Ly * 0.25) / omega0) + 1.0) + 0.5 * (Ti2 + Ti1); // Total ion temperature.
  double Te_tot = 0.5 * (Te2 - Te1) * (tanh((y - Ly * 0.25) / omega0) - tanh((y - Ly * 0.75) / omega0)
    + tanh((y - Ly * 1.25) / omega0) - tanh((y + Ly * 0.25) / omega0) + 1.0) + 0.5 * (Te2 + Te1); // Total electron temperature.
  double n = (0.5 * (b1 * b1 - b1x * b1x) + 0.5 * (guide1 * guide1 - b1z * b1z) + n1 * (Ti1 + Te1)) / (Ti_tot + Te_tot); // Total number density.

  double Bx = b1x - psi0 * 4.0 * pi / Ly * sin(2.0 * pi * x / Lx) * sin(4.0 * pi * y / Ly); // Total magnetic field (x-direction).
  double By = b1y + psi0 * 2.0 * pi / Lx * cos(2.0 * pi * x / Lx) * (1.0 - cos(4.0 * pi * y / Ly)); // Total magnetic field (y-direction).
  double Bz = b1z; // Total magnetic field (z-direction).

  double Te_frac = Te_tot / (Te_tot + Ti_tot); // Fraction of total temperature from electrons.
  double Ti_frac = Ti_tot / (Te_tot + Ti_tot); // Fraction of total temperature from ions;

  double Jx = 0.5 * (guide2 - guide1) / omega0 * ((1.0 / cosh((y - Ly * 0.25) / omega0)) * (1.0 / cosh((y - Ly * 0.25) / omega0)) 
    - (1.0 / cosh((y - Ly * 0.75) / omega0)) * (1.0 / cosh((y - Ly * 0.75) / omega0)) 
    + (1.0 / cosh((y - Ly * 1.25) / omega0)) * (1.0 / cosh((y - Ly * 1.25) / omega0)) 
    - (1.0 / cosh((y + Ly * 0.25) / omega0)) * (1.0 / cosh((y + Ly * 0.25) / omega0))); // Total current density (x-direction).
  double Jy = 0.0; // Total current density (y-direction).
  double Jz  = -0.5 * (b2 + b1) / omega0 * ((1.0 / cosh((y - Ly * 0.25) / omega0)) * (1.0 / cosh((y - Ly * 0.25) / omega0)) 
    - (1.0 / cosh((y - Ly * 0.75) / omega0)) * (1.0 / cosh((y - Ly * 0.75) / omega0)) 
    + (1.0 / cosh((y - Ly * 1.25) / omega0)) * (1.0 / cosh((y - Ly * 1.25) / omega0)) 
    - (1.0 / cosh((y + Ly * 0.25) / omega0)) * (1.0 / cosh((y + Ly * 0.25) / omega0))) 
    - psi0 * sin(2.0 * pi * x / Lx) * ((2.0 * pi / Lx) * (2.0 * pi / Lx) * (1.0 - cos(4.0 * pi * y / Ly)) +
      (4.0 * pi / Ly) * (4.0 * pi / Ly) * cos(4.0 * pi * y / Ly)); // Total current density (z-direction).
   
  double Jxe = Jx * Te_frac; // Electron current density (x-direction).
  double Jye = Jy * Te_frac; // Electron current density (y-direction).
  double Jze = Jz * Te_frac; // Electron current density (z-direction).

  double rhoe = mass_elc * n; // Electron mass density.
  double momxe = (mass_elc / charge_elc) * Jxe; // Electron momentum density (x-direction).
  double momye = (mass_elc / charge_elc) * Jye; // Electron momentum density (y-direction).
  double momze = (mass_elc / charge_elc) * Jze; // Electron momentum density (z-direction).
  double Ee_tot = n * Te_tot / (gas_gamma - 1.0) + 0.5 * (momxe * momxe + momye * momye + momze * momze) / rhoe; // Electron total energy density.

  // Set electron mass density.
  fout[0] = rhoe;
  // Set electron momentum density.
  fout[1] = momxe; fout[2] = momye; fout[3] = momze;
  // Set electron total energy density.
  fout[4] = Ee_tot;
}

void
evalIonInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1];
  struct burch_ctx *app = ctx;

  double pi = app -> pi;

  double gas_gamma = app -> gas_gamma;
  double mass_ion = app -> mass_ion;
  double charge_ion = app -> charge_ion;

  double omega0 = app -> omega0;
  double psi0 = app -> psi0;
  double guide1 = app -> guide1;
  double guide2 = app -> guide2;
  double b1 = app -> b1;
  double b2 = app -> b2;
  double n1 = app -> n1;
  double Ti2 = app -> Ti2;

  double Te1 = app -> Te1;
  double Ti1 = app -> Ti1;

  double Te2 = app -> Te2;

  double Lx = app -> Lx;
  double Ly = app -> Ly;

  double b1x = 0.5 * (b2 + b1) * (tanh((y - Ly * 0.25) / omega0) - tanh((y - Ly * 0.75) / omega0)
    + tanh((y - Ly * 1.25) / omega0) - tanh((y + Ly * 0.25) / omega0) + 1.0) + 0.5 * (b2 - b1); // Magnetospheric magnetic field (x-direction).
  double b1y = 0.0; // Magnetospheric magnetic field (y-direction).
  double b1z = 0.5 * (guide2 - guide1) * (tanh((y - Ly * 0.25) / omega0) - tanh((y - Ly * 0.75) / omega0)
    + tanh((y - Ly * 1.25) / omega0) - tanh((y + Ly * 0.25) / omega0) + 1.0) + 0.5 * (guide2 + guide1); // Magnetospheric magnetic field (z-direction).

  double Ti_tot = 0.5 * (Ti2 - Ti1) * (tanh((y - Ly * 0.25) / omega0) - tanh((y - Ly * 0.75) / omega0)
    + tanh((y - Ly * 1.25) / omega0) - tanh((y + Ly * 0.25) / omega0) + 1.0) + 0.5 * (Ti2 + Ti1); // Total ion temperature.
  double Te_tot = 0.5 * (Te2 - Te1) * (tanh((y - Ly * 0.25) / omega0) - tanh((y - Ly * 0.75) / omega0)
    + tanh((y - Ly * 1.25) / omega0) - tanh((y + Ly * 0.25) / omega0) + 1.0) + 0.5 * (Te2 + Te1); // Total electron temperature.
  double n = (0.5 * (b1 * b1 - b1x * b1x) + 0.5 * (guide1 * guide1 - b1z * b1z) + n1 * (Ti1 + Te1)) / (Ti_tot + Te_tot); // Total number density.

  double Bx = b1x - psi0 * 4.0 * pi / Ly * sin(2.0 * pi * x / Lx) * sin(4.0 * pi * y / Ly); // Total magnetic field (x-direction).
  double By = b1y + psi0 * 2.0 * pi / Lx * cos(2.0 * pi * x / Lx) * (1.0 - cos(4.0 * pi * y / Ly)); // Total magnetic field (y-direction).
  double Bz = b1z; // Total magnetic field (z-direction).

  double Te_frac = Te_tot / (Te_tot + Ti_tot); // Fraction of total temperature from electrons.
  double Ti_frac = Ti_tot / (Te_tot + Ti_tot); // Fraction of total temperature from ions;

  double Jx = 0.5 * (guide2 - guide1) / omega0 * ((1.0 / cosh((y - Ly * 0.25) / omega0)) * (1.0 / cosh((y - Ly * 0.25) / omega0)) 
    - (1.0 / cosh((y - Ly * 0.75) / omega0)) * (1.0 / cosh((y - Ly * 0.75) / omega0)) 
    + (1.0 / cosh((y - Ly * 1.25) / omega0)) * (1.0 / cosh((y - Ly * 1.25) / omega0)) 
    - (1.0 / cosh((y + Ly * 0.25) / omega0)) * (1.0 / cosh((y + Ly * 0.25) / omega0))); // Total current density (x-direction).
  double Jy = 0.0; // Total current density (y-direction).
  double Jz  = -0.5 * (b2 + b1) / omega0 * ((1.0 / cosh((y - Ly * 0.25) / omega0)) * (1.0 / cosh((y - Ly * 0.25) / omega0)) 
    - (1.0 / cosh((y - Ly * 0.75) / omega0)) * (1.0 / cosh((y - Ly * 0.75) / omega0)) 
    + (1.0 / cosh((y - Ly * 1.25) / omega0)) * (1.0 / cosh((y - Ly * 1.25) / omega0)) 
    - (1.0 / cosh((y + Ly * 0.25) / omega0)) * (1.0 / cosh((y + Ly * 0.25) / omega0))) 
    - psi0 * sin(2.0 * pi * x / Lx) * ((2.0 * pi / Lx) * (2.0 * pi / Lx) * (1.0 - cos(4.0 * pi * y / Ly)) +
      (4.0 * pi / Ly) * (4.0 * pi / Ly) * cos(4.0 * pi * y / Ly)); // Total current density (z-direction).
   
  double Jxi = Jx * Ti_frac; // Ion current density (x-direction).
  double Jyi = Jy * Ti_frac; // Ion current density (y-direction).
  double Jzi = Jz * Ti_frac; // Ion current density (z-direction).

  double rhoi = mass_ion * n; // Ion mass density.
  double momxi = (mass_ion / charge_ion) * Jxi; // Ion momentum density (x-direction).
  double momyi = (mass_ion / charge_ion) * Jyi; // Ion momentum density (y-direction).
  double momzi = (mass_ion / charge_ion) * Jzi; // Ion momentum density (z-direction).
  double Ei_tot = n * Ti_tot / (gas_gamma - 1.0) + 0.5 * (momxi * momxi + momyi * momyi + momzi * momzi) / rhoi; // Ion total energy density.

  // Set ion mass density,
  fout[0] = rhoi;
  // Set ion momentum density.
  fout[1] = momxi; fout[2] = momyi; fout[3] = momzi;
  // Set ion total energy density.
  fout[4] = Ei_tot;
}

void
evalFieldInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1];
  struct burch_ctx *app = ctx;

  double pi = app -> pi;

  double omega0 = app -> omega0;
  double psi0 = app -> psi0;
  double guide1 = app -> guide1;
  double guide2 = app -> guide2;
  double b1 = app -> b1;
  double b2 = app -> b2;
  double n1 = app -> n1;
  double Ti2 = app -> Ti2;

  double Te1 = app -> Te1;
  double Ti1 = app -> Ti1;

  double Te2 = app -> Te2;

  double Lx = app -> Lx;
  double Ly = app -> Ly;

  double b1x = 0.5 * (b2 + b1) * (tanh((y - Ly * 0.25) / omega0) - tanh((y - Ly * 0.75) / omega0)
    + tanh((y - Ly * 1.25) / omega0) - tanh((y + Ly * 0.25) / omega0) + 1.0) + 0.5 * (b2 - b1); // Magnetospheric magnetic field (x-direction).
  double b1y = 0.0; // Magnetospheric magnetic field (y-direction).
  double b1z = 0.5 * (guide2 - guide1) * (tanh((y - Ly * 0.25) / omega0) - tanh((y - Ly * 0.75) / omega0)
    + tanh((y - Ly * 1.25) / omega0) - tanh((y + Ly * 0.25) / omega0) + 1.0) + 0.5 * (guide2 + guide1); // Magnetospheric magnetic field (z-direction).

  double Ti_tot = 0.5 * (Ti2 - Ti1) * (tanh((y - Ly * 0.25) / omega0) - tanh((y - Ly * 0.75) / omega0)
    + tanh((y - Ly * 1.25) / omega0) - tanh((y + Ly * 0.25) / omega0) + 1.0) + 0.5 * (Ti2 + Ti1); // Total ion temperature.
  double Te_tot = 0.5 * (Te2 - Te1) * (tanh((y - Ly * 0.25) / omega0) - tanh((y - Ly * 0.75) / omega0)
    + tanh((y - Ly * 1.25) / omega0) - tanh((y + Ly * 0.25) / omega0) + 1.0) + 0.5 * (Te2 + Te1); // Total electron temperature.
  double n = (0.5 * (b1 * b1 - b1x * b1x) + 0.5 * (guide1 * guide1 - b1z * b1z) + n1 * (Ti1 + Te1)) / (Ti_tot + Te_tot); // Total number density.

  double Bx = b1x - psi0 * 4.0 * pi / Ly * sin(2.0 * pi * x / Lx) * sin(4.0 * pi * y / Ly); // Total magnetic field (x-direction).
  double By = b1y + psi0 * 2.0 * pi / Lx * cos(2.0 * pi * x / Lx) * (1.0 - cos(4.0 * pi * y / Ly)); // Total magnetic field (y-direction).
  double Bz = b1z; // Total magnetic field (z-direction).

  // Set electric field.
  fout[0] = 0.0, fout[1] = 0.0; fout[2] = 0.0;
  // Set magnetic field.
  fout[3] = Bx, fout[4] = By; fout[5] = Bz;
  // Set correction potentials.
  fout[6] = 0.0; fout[7] = 0.0;
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

  struct burch_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], ctx.Ny);

  // Electron/ion equations.
  struct gkyl_wv_eqn *elc_euler = gkyl_wv_euler_new(ctx.gas_gamma, app_args.use_gpu);
  struct gkyl_wv_eqn *ion_euler = gkyl_wv_euler_new(ctx.gas_gamma, app_args.use_gpu);

  struct gkyl_moment_species elc = {
    .name = "elc",
    .charge = ctx.charge_elc, .mass = ctx.mass_elc,
    .equation = elc_euler,
    .evolve = true,
    .init = evalElcInit,
    .ctx = &ctx,
  };
  struct gkyl_moment_species ion = {
    .name = "ion",
    .charge = ctx.charge_ion, .mass = ctx.mass_ion,
    .equation = ion_euler,
    .evolve = true,
    .init = evalIonInit,
    .ctx = &ctx,
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
  comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp)
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
    .name = "5m_burch",

    .ndim = 2,
    .lower = { 0.0, 0.0 },
    .upper = { ctx.Lx / ctx.di, ctx.Ly / ctx.di }, 
    .cells = { NX, NY },

    .num_periodic_dir = 2,
    .periodic_dirs = { 0, 1 },
    .cfl_frac = ctx.cfl_frac,

    .num_species = 2,
    .species = { elc, ion },

    .field = {
      .epsilon0 = ctx.epsilon0, .mu0 = ctx.mu0,
      .mag_error_speed_fact = 1.0,
      
      .evolve = true,
      .init = evalFieldInit,
      .ctx = &ctx,
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
