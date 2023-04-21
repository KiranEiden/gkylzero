#pragma once

#include <gkyl_eqn_type.h>
#include <gkyl_mom_type.h>
#include <gkyl_mom_calc.h>

struct gkyl_dg_updater_moment {
  enum gkyl_model_id model_id; // Identifier for model (e.g., SR, PKPM, see gkyl_eqn_type.h)
  struct gkyl_mom_type *type; // Moment type
  struct gkyl_mom_calc *up_moment; // Updater for computing moment

  double moment_tm; // total time spent in computing moment
};