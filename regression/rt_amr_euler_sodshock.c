#include <gkyl_amr_core.h>

struct amr_euler_sodshock_ctx
{
  // Physical constants (using normalized code units).
  double gas_gamma; // Adiabatic index.

  double rhol; // Left fluid mass density.
  double ul; // Left fluid velocity.
  double pl; // Left fluid pressure.

  double rhor; // Right fluid mass density.
  double ur; // Right fluid velocity.
  double pr; // Right fluid pressure.

  // Simulation parameters.
  int Nx; // Coarse cell count (x-direction).
  int ref_factor; // Refinement factor.
  double Lx; // Coarse domain size (x-direction).
  double fine_Lx; // Fine domain size (x-direction).
  double cfl_frac; // CFL coefficient.
  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
};

struct amr_euler_sodshock_ctx
create_ctx(void)
{
  // Physical constants (using normalized code units).
  double gas_gamma = 1.4; // Adiabatic index.

  double rhol = 3.0; // Left fluid mass density.
  double ul = 0.0; // Left fluid velocity.
  double pl = 3.0; // Left fluid pressure.

  double rhor = 1.0; // Right fluid mass density.
  double ur = 0.0; // Right fluid velocity.
  double pr = 1.0; // Right fluid pressure.

  // Simulation parameters.
  int Nx = 32; // Coarse cell count (x-direction).
  int ref_factor = 4; // Refinement factor.
  double Lx = 1.0; // Coarse domain size (x-direction).
  double fine_Lx = 0.2; // Fine domain size (x-direction).
  double cfl_frac = 0.95; // CFL coefficient.
  double t_end = 0.1; // Final simulation time.
  int num_frames = 1; // Number of output frames.

  struct amr_euler_sodshock_ctx ctx = {
    .gas_gamma = gas_gamma,
    .rhol = rhol,
    .ul = ul,
    .pl = pl,
    .rhor = rhor,
    .ur = ur,
    .pr = pr,
    .Nx = Nx,
    .ref_factor = ref_factor,
    .Lx = Lx,
    .fine_Lx = fine_Lx,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
    .num_frames = num_frames,
  };

  return ctx;
}

void
evalEulerInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0];
  struct amr_euler_sodshock_ctx new_ctx = create_ctx(); // Context for initialization functions.
  struct amr_euler_sodshock_ctx *app = &new_ctx;

  double gas_gamma = app->gas_gamma;

  double rhol = app->rhol;
  double ul = app->ul;
  double pl = app->pl;

  double rhor = app->rhor;
  double ur = app->ur;
  double pr = app->pr;

  double rho = 0.0;
  double u = 0.0;
  double p = 0.0;

  if (x < 0.75) {
    rho = rhol; // Fluid mass density (left).
    u = ul; // Fluid velocity (left).
    p = pl; // Fluid pressure (left).
  }
  else {
    rho = rhor; // Fluid mass density (right).
    u = ur; // Fluid velocity (right).
    p = pr; // Fluid pressure (right).
  }
  
  // Set fluid mass density.
  fout[0] = rho;
  // Set fluid momentum density.
  fout[1] = rho * u; fout[2] = 0.0; fout[3] = 0.0;
  // Set fluid total energy density.
  fout[4] = p / (gas_gamma - 1.0) + 0.5 * rho * u * u;
}

int main(int argc, char **argv)
{
  struct amr_euler_sodshock_ctx ctx = create_ctx(); // Context for initialization functions.

  struct euler1d_single_init init = {
    .base_Nx = ctx.Nx,
    .ref_factor = ctx.ref_factor,

    .coarse_x1 = 0.25,
    .coarse_x2 = 0.25 + ctx.Lx,

    .refined_x1 = (0.25 + (0.5 * ctx.Lx)) - (0.5 * ctx.fine_Lx),
    .refined_x2 = (0.25 + (0.5 * ctx.Lx)) + (0.5 * ctx.fine_Lx),

    .eval = evalEulerInit,
    .gas_gamma = ctx.gas_gamma,

    .low_order_flux = false,
    .cfl_frac = ctx.cfl_frac,
    .t_end = ctx.t_end,
    .num_frames = ctx.num_frames,
  };

  euler1d_run_single(argc, argv, &init);
}