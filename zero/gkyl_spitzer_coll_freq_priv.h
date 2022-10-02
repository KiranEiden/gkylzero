#include <gkyl_spitzer_coll_freq.h>

struct gkyl_spitzer_coll_freq {
  int num_quad; // number of quadrature points to use in each direction
  int ndim; // Configuration-space dimension

  int num_basis; // number of conf-space basis functions

  bool use_gpu;

  // for quadrature in conf space
  int tot_quad; // total number of quadrature points
  struct gkyl_array *ordinates; // conf-space ordinates for quadrature
  struct gkyl_array *weights; // weights for conf-space quadrature
  struct gkyl_array *basis_at_ords; // conf-space basis functions at ordinates

  struct gkyl_array *fun_at_ords; // function (Maxwellian) evaluated at
                                  // ordinates in a cell.
};

void
gkyl_spitzer_coll_freq_advance_normnu_cu(const gkyl_spitzer_coll_freq *up,
  const struct gkyl_range *range, const struct gkyl_array *vtSqSelf,
  const struct gkyl_array *m0Other, const struct gkyl_array *vtSqOther,
  double normNu, struct gkyl_array *nuOut);

void
gkyl_spitzer_coll_freq_advance_cu(const gkyl_spitzer_coll_freq *up,
  const struct gkyl_range *range,
  double qSelf, double mSelf, const struct gkyl_array *m0Self, const struct gkyl_array *vtSqSelf,
  double qOther, double mOther, const struct gkyl_array *m0Other, const struct gkyl_array *vtSqOther,
  struct gkyl_array *nuOut);
