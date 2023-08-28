#pragma once

#include <gkyl_hyper_dg.h>
#include <gkyl_dg_eqn.h>

struct gkyl_dg_updater_diffusion {
  struct gkyl_dg_eqn *dgeqn; // Equation object.
  struct gkyl_hyper_dg *hyperdg; // solvers for specific diffusion equation.
  enum gkyl_diffusion_id diffid; // Type of diffusion model.
  bool use_gpu;

  double diffusion_tm; // total time spent in computing diffusion equation.
};
