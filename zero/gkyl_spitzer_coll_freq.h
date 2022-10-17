#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>

// Object type
typedef struct gkyl_spitzer_coll_freq gkyl_spitzer_coll_freq;

/**
 * Create new updater to either compute the Spitzer collision frequency from
 * scratch based on local parameters, or scale a normalized collision frequency
 * by the local n_r/(v_ts^2+v_tr^2)^(3/2).
 *
 * @param basis Basis object (configuration space).
 * @param num_quad Number of quadrature nodes.
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_spitzer_coll_freq* gkyl_spitzer_coll_freq_new(
  const struct gkyl_basis *basis, int num_quad, bool use_gpu);

/**
 * Scale the normalized collision frequency, normNu, by
 * n_r/(v_ts^2+v_tr^2)^(3/2) and project it on to the basis.
 *
 * @param up Spizer collision frequency updater object.
 * @param range Config-space range
 * @param vtSqSelf Thermal speed squared of this species. 
 * @param m0Other Thermal speed squared of the other species. 
 * @param vtSqOther Thermal speed squared of the other species. 
 * @param normNu Normalized collision frequency to scale.
 * @param nuOut Output collision frequency.
 */
void gkyl_spitzer_coll_freq_advance_normnu(const gkyl_spitzer_coll_freq *up,
  const struct gkyl_range *range, const struct gkyl_array *vtSqSelf,
  const struct gkyl_array *m0Other, const struct gkyl_array *vtSqOther,
  double normNu, struct gkyl_array *nuOut);

/**
 * Compute the Spitzer collision frequency from scratch. Coulomb Logarithm
 * is computed using cell averaged values.
 *
 * @param up Spizer collision frequency updater object.
 * @param range Config-space range
 * @param qSelf Charge of this species.
 * @param mSelf Mass of this species.
 * @param m0Self Thermal speed squared of the other species. 
 * @param vtSqSelf Thermal speed squared of this species. 
 * @param qOther Charge of this species.
 * @param mOther Mass of this species.
 * @param m0Other Thermal speed squared of the other species. 
 * @param vtSqOther Thermal speed squared of the other species. 
 * @param nuOut Output collision frequency.
 */
void gkyl_spitzer_coll_freq_advance(const gkyl_spitzer_coll_freq *up,
  const struct gkyl_range *range,
  double qSelf, double mSelf, const struct gkyl_array *m0Self, const struct gkyl_array *vtSqSelf,
  double qOther, double mOther, const struct gkyl_array *m0Other, const struct gkyl_array *vtSqOther,
  struct gkyl_array *nuOut);

/**
 * Delete updater.
 *
 * @param pob Updater to delete.
 */
void gkyl_spitzer_coll_freq_release(gkyl_spitzer_coll_freq* up);
