#include <assert.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_cx.h>
#include <gkyl_dg_cx_priv.h>
#include <gkyl_util.h>

// "Choose Kernel" based on cdim and polyorder
#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order]

gkyl_dg_cx*
gkyl_dg_cx_new(const struct gkyl_rect_grid *grid,
  struct gkyl_basis *conf_basis,
  struct gkyl_basis *phase_basis, double mass_ion,
  enum gkyl_dg_cx_type type_ion, 
  bool use_gpu)
{
  gkyl_dg_cx *up = gkyl_malloc(sizeof(gkyl_dg_cx));

  int num_basis = phase_basis->num_basis;
  int cdim = conf_basis->ndim;
  int vdim = phase_basis->ndim - cdim;
  int poly_order = conf_basis->poly_order;
  up->grid = *grid; 
  up->cdim = cdim;
  up->vdim = vdim; 
  up->poly_order = poly_order;

  /* up->elem_charge = elem_charge; */
  up->mass_ion = mass_ion;
  /* up->mass_neut = mass_neut; */
  if (type_ion == GKYL_H) {
    up->a = 1.12e-18;
    up->b = 7.15e-20;
  }
  else if (type_ion == GKYL_D) {
    up->a = 1.09e-18;
    up->b = 7.15e-20;
  }
  else if (type_ion == GKYL_NE) {
    up->a = 7.95e-19;
    up->b = 5.65e-20;
  }
  else if (type_ion == GKYL_HE) {
    up->a = 6.484e-19;
    up->b = 4.350e-20;
  }  

  up->react_rate = CK(ser_cx_react_rate_kernels, cdim, poly_order);
  
  return up;
}

// need m0Neut, uIon, uNeut, vtSqIon, vtSqNeut
void gkyl_dg_cx_react_rate(const struct gkyl_dg_cx *cx,
  const struct gkyl_range *conf_rng, const struct gkyl_range *phase_rng, 
  const struct gkyl_array *n_neut,  const struct gkyl_array *u_neut, const struct gkyl_array *vth_sq_neut,
  const struct gkyl_array *u_ion, const struct gkyl_array *vth_sq_ion,
  struct gkyl_array *cflrate, struct gkyl_array *coef_cx)
{
  struct gkyl_range vel_rng;
  struct gkyl_range_iter conf_iter, vel_iter;

  int rem_dir[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<conf_rng->ndim; ++d) rem_dir[d] = 1;

  gkyl_range_iter_init(&conf_iter, conf_rng);
  while (gkyl_range_iter_next(&conf_iter)) {
    long loc = gkyl_range_idx(conf_rng, conf_iter.idx);

    const double *n_neut_d = gkyl_array_cfetch(n_neut, loc);
    const double *u_neut_d = gkyl_array_cfetch(u_neut, loc);
    const double *u_ion_d = gkyl_array_cfetch(u_ion, loc);
    const double *vth_sq_neut_d = gkyl_array_cfetch(vth_sq_neut, loc);
    const double *vth_sq_ion_d = gkyl_array_cfetch(vth_sq_ion, loc);
    double *coef_cx_d = gkyl_array_fetch(coef_cx, loc);

    // Calculate vt_sq min for ion, neut (use same for now to test 1x1v)
    double vth_sq_ion_min;
    double vth_sq_neut_min;
    double TempMin = 0.0; 
    for (int d=0; d<cx->vdim; d++) {
      TempMin = TempMin + (1./3.)*(cx->mass_ion/6.)*cx->grid.dx[cx->cdim+d];
    }
    vth_sq_ion_min = TempMin/cx->mass_ion;
    vth_sq_neut_min = TempMin/cx->mass_ion;
    
    double cflr = cx->react_rate(cx->a, cx->b,
      n_neut_d, u_ion_d, u_neut_d, vth_sq_ion_d,
      vth_sq_ion_min, vth_sq_neut_d, vth_sq_neut_min,
      coef_cx_d);

    gkyl_range_deflate(&vel_rng, phase_rng, rem_dir, conf_iter.idx);
    gkyl_range_iter_no_split_init(&vel_iter, &vel_rng);
    // cfl associated with reaction is a *phase space* cfl
    // Need to loop over velocity space for each configuration space cell
    // to get total cfl rate in each phase space cell
    while (gkyl_range_iter_next(&vel_iter)) {
      long cfl_idx = gkyl_range_idx(&vel_rng, vel_iter.idx);
      double *cflrate_d = gkyl_array_fetch(cflrate, cfl_idx);
      cflrate_d[0] += cflr; // frequencies are additive
      }
  }
}

void
gkyl_dg_cx_release(gkyl_dg_cx* cx)
{
  free(cx);
}
