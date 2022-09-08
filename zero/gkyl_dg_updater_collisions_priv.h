#pragma once

#include <gkyl_alloc.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_hyper_dg.h>

struct gkyl_dg_updater_collisions {
  enum gkyl_model_id model_id; // Type of LBO (e.g., Vlasov vs. PKPM model)
  
  struct gkyl_dg_eqn *coll_drag; // Collision drag equation
  struct gkyl_dg_eqn *coll_diff; // Collision diffusion equation
  struct gkyl_hyper_dg *drag; // solvers for drag terms
  struct gkyl_hyper_dg *diff; // solvers for diffusion terms

  double drag_tm; // total time spent in computing drag terms
  double diff_tm; // totat time spent in computing diffusion terms
};
