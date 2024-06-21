#pragma once

#include <gkyl_wv_eqn.h>

// Type of Riemann solver to use:
enum gkyl_wv_maxwell_cyl_rp {
  WV_MAXWELL_CYL_RP_LAX = 0, // Default (Lax fluxes).
  WV_MAXWELL_CYL_RP_EIG, // Solve eigensystem to find fluctuations
};

/**
 * Create a new maxwell equation object.
 * 
 * @param c Speed of light
 * @param rp_type Type of Riemann solver
 * @return Pointer to Maxwell equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_maxwell_cyl_new(double c, enum gkyl_wv_maxwell_cyl_rp rp_type);
