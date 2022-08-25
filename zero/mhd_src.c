#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_mhd_src.h>


// Makes indexing cleaner
#define DN (0)
#define MX (1)
#define MY (2)
#define MZ (3)
#define ER (4)
#define BX (5)
#define BY (6)
#define BZ (7)
#define PSI_GLM (8)

struct gkyl_mhd_src {
  struct gkyl_rect_grid grid;
  int ndim;
  enum gkyl_wv_mhd_div_constraint divergence_constraint;
  double glm_ch;
  double glm_cp;
  double dxyz_min;
  double cfl; // cfl used in the hyperbolic part
};

gkyl_mhd_src *gkyl_mhd_src_new(struct gkyl_mhd_src_inp inp) {
  gkyl_mhd_src *up = gkyl_malloc(sizeof(gkyl_mhd_src));

  up->grid = *(inp.grid);
  up->ndim = up->grid.ndim;
  up->divergence_constraint = inp.divergence_constraint;
  up->glm_ch = inp.glm_ch;
  up->glm_cp = inp.glm_cp;

  return up;
}

void gkyl_mhd_src_advance(const gkyl_mhd_src *up, double dt,
                          const struct gkyl_range *update_range,
                          struct gkyl_array *q_array,
                          const struct gkyl_array *acc_array) {
  int div_type = up->divergence_constraint;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);
  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(update_range, iter.idx);

    double *q = gkyl_array_fetch(q_array, lidx);
    const double *acc = gkyl_array_cfetch(acc_array, lidx);

    if (div_type == GKYL_MHD_DIVB_EIGHT_WAVES) {
      // TODO
    } else if (div_type == GKYL_MHD_DIVB_GLM) {
      // Dedner et al. JCP (2002) 175, 645, Equation (19).
      // Mignone & Tzeferacos, JCP (2010) 229, 2117, Equation (27).
      double ch = up->glm_ch;
      double cp = up->glm_cp;
      q[PSI_GLM] *= exp(- ch * ch / cp / cp * dt);
    }
  }
}

void gkyl_mhd_src_release(gkyl_mhd_src *up) { gkyl_free(up); }
