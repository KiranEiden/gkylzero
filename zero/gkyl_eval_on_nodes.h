#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_evalf_def.h>

// Object type
typedef struct gkyl_eval_on_nodes gkyl_eval_on_nodes;

/**
 * Create new updater to compute function on nodes and calculate its
 * expansion on basis functions. Free using gkyl_eval_on_nodes_release
 * method.
 *
 * @param grid Grid object
 * @param basis Basis functions to project on
 * @param num_ret_vals Number of values 'eval' sets
 * @param eval Function to project.
 * @param ctx Context for function evaluation. Can be NULL.
 * @return New updater pointer.
 */
struct gkyl_eval_on_nodes* gkyl_eval_on_nodes_new(
  const struct gkyl_rect_grid *grid, const struct gkyl_basis *basis,
  int num_ret_vals, evalf_t eval, void *ctx);

/**
 * Compute evaluation on nodes and corresponding expansion
 * coefficients. The update_rng MUST be a sub-range of the range on
 * which the array is defined. That is, it must be either the same
 * range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param up Eval on nodes updater to run
 * @param tm Time at which eval must be computed
 * @param update_rng Range on which to run eval.
 * @param out Output array
 */
void gkyl_eval_on_nodes_advance(const struct gkyl_eval_on_nodes *up,
  double tm, const struct gkyl_range *update_rng, struct gkyl_array *out);

/**
 * Perform the nodal to modal transformation.
 *
 * @param up Project on basis updater.
 * @param fun_at_nodes Function evaluated at nodes in one cell.
 * @param f Modal coefficients of the function in one cell.
 */
void gkyl_eval_on_nodes_nod2mod(const struct gkyl_eval_on_nodes *up, const struct gkyl_array *fun_at_nodes, double *f);

/**
 * Get the coordinates of a given node.
 *
 * @param up Project on basis updater.
 * @param node Index indicate the desired node.
 * @return Node coordinates.
 */
double* gkyl_eval_on_nodes_fetch_node(const struct gkyl_eval_on_nodes *up, long node);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_eval_on_nodes_release(struct gkyl_eval_on_nodes *up);
