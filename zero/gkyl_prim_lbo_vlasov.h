#pragma once

#include <gkyl_basis.h>
#include <gkyl_prim_lbo.h>

/**
 * Create a new Vlasov primitive moment object.
 *
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis functions
 * @return Pointer to Vlasov primitive moment object
 */
struct gkyl_prim_lbo* gkyl_prim_lbo_vlasov_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis);