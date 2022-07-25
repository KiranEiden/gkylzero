#pragma once

#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_array.h>

// BC types in this updater.
enum gkyl_bc_basic_type { BC_ABSORB, BC_REFLECT };

// Object type
typedef struct gkyl_bc_basic gkyl_bc_basic;

/**
 * Create new updater to apply basic BCs to a field
 * in a gkyl_array. Basic BCs are those in which the
 * ghost cell depends solely on the skin cell next to it
 * via a function of type array_copy_func_t (e.g. absorb, reflect).
 *
 * @param dir Direction in which to apply BC.
 * @param edge Lower or upper edge at which to apply BC (see gkyl_edge_loc).
 * @param local_range_ext Local extended range.
 * @param num_ghosts Number of ghosts in each dimension.
 * @param bctype BC type (see gkyl_bc_basic_type).
 * @param basis Basis on which coefficients in array are expanded.
 * @param cdim Configuration space dimensions.
 * @param use_gpu Boolean to indicate whether to use the GPU.
 * @return New updater pointer.
 */
struct gkyl_bc_basic* gkyl_bc_basic_new(int dir, enum gkyl_edge_loc edge, const struct gkyl_range* local_range_ext,
  const int *num_ghosts, enum gkyl_bc_basic_type bctype, const struct gkyl_basis *basis, int cdim, bool use_gpu);

/**
 * Create new updater to apply basic BCs to a field
 * in a gkyl_array. Basic BCs are those in which the
 * ghost cell depends solely on the skin cell next to it
 * via a function of type array_copy_func_t (e.g. absorb, reflect).
 *
 * @param up BC updater.
 * @param buff_arr Buffer array, big enough for ghost cells at this boundary.
 * @param f_arr Field array to apply BC to.
 */
void gkyl_bc_basic_advance(const struct gkyl_bc_basic *up, struct gkyl_array *buff_arr, struct gkyl_array *f_arr);

/**
 * Free memory associated with bc_basic updater.
 *
 * @param up BC updater.
 */
void gkyl_bc_basic_release(struct gkyl_bc_basic *up);
