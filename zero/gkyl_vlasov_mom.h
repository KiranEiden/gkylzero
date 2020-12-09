#pragma once

#include <gkyl_basis.h>
#include <gkyl_mom_type.h>

/**
 * Create new Vlasov moment type object. Valid 'mom' strings are "M0",
 * "M1i", "M2", "M2ij", "M3i"
 *
 * @param cbasis Configuration-space basis-functions
 * @param pbasis Phase-space basis-functions
 * @param mom Name of moment to compute.
 */
struct gkyl_mom_type* gkyl_vlasov_mom_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const char *mom);