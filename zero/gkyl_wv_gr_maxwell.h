#pragma once

#include <gkyl_wv_eqn.h>

// Type of Riemann solver to use:
enum gkyl_wv_gr_maxwell_rp {
  WV_GR_MAXWELL_RP_EIG = 0, // Solve eigensystem to find fluctuations (default)
  WV_GR_MAXWELL_RP_LAX // Lax fluxes
};

/**
 * Create a new GR Maxwell equation object.
 * 
 * @param rp_type Type of Riemann solver to use 
 * @return Pointer to GR Maxwell equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_gr_maxwell_new(enum gkyl_wv_gr_maxwell_rp rp_type);
