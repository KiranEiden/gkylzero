#include <string.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_const.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_proj_maxwellian_on_basis.h>
#include <gkyl_range.h>

struct gkyl_proj_maxwellian_on_basis {
  struct gkyl_rect_grid grid;
  int num_quad; // number of quadrature points to use in each direction
  int cdim; // Configuration-space dimension
  int pdim; // Phase-space dimension

  int num_conf_basis; // number of conf-space basis functions
  int num_phase_basis; // number of phase-space basis functions

  // for quadrature in phase-space
  int tot_quad; // total number of quadrature points
  struct gkyl_array *ordinates; // ordinates for quadrature
  struct gkyl_array *weights; // weights for quadrature
  struct gkyl_array *basis_at_ords; // basis functions at ordinates

  // for quadrature in conf space
  int tot_conf_quad; // total number of quadrature points
  struct gkyl_array *conf_ordinates; // conf-space ordinates for quadrature
  struct gkyl_array *conf_weights; // weights for conf-space quadrature  
  struct gkyl_array *conf_basis_at_ords; // conf-space basis functions at ordinates
};

// Sets ordinates, weights and basis functions at ords. Returns total
// number of quadrature nodes
static int
init_quad_values(const struct gkyl_basis *basis, int num_quad,
  struct gkyl_array **ordinates, struct gkyl_array **weights, struct gkyl_array **basis_at_ords)
{
  int ndim = basis->ndim;
  double ordinates1[num_quad], weights1[num_quad];

  if (num_quad <= gkyl_gauss_max) {
    // use pre-computed values if possible (these are more accurate
    // than computing them on the fly)
    memcpy(ordinates1, gkyl_gauss_ordinates[num_quad], sizeof(double[num_quad]));
    memcpy(weights1, gkyl_gauss_weights[num_quad], sizeof(double[num_quad]));
  }
  else {
    gkyl_gauleg(-1, 1, ordinates1, weights1, num_quad);
  }

  // create range to loop over quadrature points
  int qshape[GKYL_MAX_DIM];
  for (int i=0; i<ndim; ++i) qshape[i] = num_quad;
  struct gkyl_range qrange;
  gkyl_range_init_from_shape(&qrange, ndim, qshape);

  int tot_quad = qrange.volume;

  // create ordinates and weights for multi-D quadrature
  *ordinates = gkyl_array_new(GKYL_DOUBLE, ndim, tot_quad);
  *weights = gkyl_array_new(GKYL_DOUBLE, 1, tot_quad);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &qrange);

  while (gkyl_range_iter_next(&iter)) {
    long node = gkyl_range_idx(&qrange, iter.idx);
    
    // set ordinates
    double *ord = gkyl_array_fetch(*ordinates, node);
    for (int i=0; i<ndim; ++i)
      ord[i] = ordinates1[iter.idx[i]-qrange.lower[i]];
    
    // set weights
    double *wgt = gkyl_array_fetch(*weights, node);
    wgt[0] = 1.0;
    for (int i=0; i<qrange.ndim; ++i)
      wgt[0] *= weights1[iter.idx[i]-qrange.lower[i]];
  }

  // pre-compute basis functions at ordinates
  *basis_at_ords = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, tot_quad);
  for (int n=0; n<tot_quad; ++n)
    basis->eval(gkyl_array_fetch(*ordinates, n), gkyl_array_fetch(*basis_at_ords, n));

  return tot_quad;
}

gkyl_proj_maxwellian_on_basis*
gkyl_proj_maxwellian_on_basis_new(
  const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis,
  int num_quad)
{
  gkyl_proj_maxwellian_on_basis *up = gkyl_malloc(sizeof(gkyl_proj_maxwellian_on_basis));

  up->grid = *grid;
  up->num_quad = num_quad;
  up->cdim = conf_basis->ndim;
  up->pdim = phase_basis->ndim;
  up->num_conf_basis = conf_basis->num_basis;
  up->num_phase_basis = phase_basis->num_basis;

  // initialize data needed for phase-space quadrature 
  up->tot_quad = init_quad_values(phase_basis, num_quad,
    &up->ordinates, &up->weights, &up->basis_at_ords);

  // initialize data needed for conf-space quadrature 
  up->tot_conf_quad = init_quad_values(conf_basis, num_quad,
    &up->conf_ordinates, &up->conf_weights, &up->conf_basis_at_ords);

  return up;
}

static inline void
comp_to_phys(int ndim, const double *eta,
  const double * GKYL_RESTRICT dx, const double * GKYL_RESTRICT xc,
  double* GKYL_RESTRICT xout)
{
  for (int d=0; d<ndim; ++d) xout[d] = 0.5*dx[d]*eta[d]+xc[d];
}

static void
proj_on_basis(const gkyl_proj_maxwellian_on_basis *up, const struct gkyl_array *fun_at_ords, double* f)
{
  int num_basis = up->num_phase_basis;
  int tot_quad = up->tot_quad;

  const double* GKYL_RESTRICT weights = up->weights->data;
  const double* GKYL_RESTRICT basis_at_ords = up->basis_at_ords->data;
  const double* GKYL_RESTRICT func_at_ords = fun_at_ords->data;

  for (int k=0; k<num_basis; ++k) f[k] = 0.0;
  
  for (int imu=0; imu<tot_quad; ++imu) {
    double tmp = weights[imu]*func_at_ords[imu];
    for (int k=0; k<num_basis; ++k)
      f[k] += tmp*basis_at_ords[k+num_basis*imu];
  }
}

static inline void
copy_idx_arrays(int cdim, int pdim, const int *cidx, const int *vidx, int *out)
{
  for (int i=0; i<cdim; ++i)
    out[i] = cidx[i];
  for (int i=cdim; i<pdim; ++i)
    out[i] = vidx[i-cdim];
}

void
gkyl_proj_maxwellian_on_basis_lab_mom(const gkyl_proj_maxwellian_on_basis *up,
  const struct gkyl_range phase_rng, const struct gkyl_range conf_rng,
  const struct gkyl_array *M0, const struct gkyl_array *M1i, const struct gkyl_array *M2,
  struct gkyl_array *fmax)
{
  int cdim = up->cdim, pdim = up->pdim;
  int vdim = pdim-cdim;
  int num_quad = up->num_quad;
  int tot_quad = up->tot_quad;
  int num_phase_basis = up->num_phase_basis;  

  int tot_conf_quad = up->tot_conf_quad;
  int num_conf_basis = up->num_conf_basis;

  int qshape[GKYL_MAX_DIM];

  // create range to loop over config-space and phase-space quadrature points
  for (int i=0; i<cdim; ++i) qshape[i] = num_quad;
  struct gkyl_range conf_qrange;
  gkyl_range_init_from_shape(&conf_qrange, cdim, qshape);

  for (int i=0; i<pdim; ++i) qshape[i] = num_quad;
  struct gkyl_range phase_qrange;
  gkyl_range_init_from_shape(&phase_qrange, pdim, qshape);

  struct gkyl_array *fun_at_ords = gkyl_array_new(GKYL_DOUBLE, 1, tot_quad);

  struct gkyl_range vel_rng;
  struct gkyl_range_iter conf_iter, vel_iter;
  
  int pidx[GKYL_MAX_DIM], rem_dir[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<conf_rng.ndim; ++d) rem_dir[d] = 1;

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];
  double num[tot_conf_quad], vel[tot_conf_quad][vdim], vth2[tot_conf_quad];
  
  // outer loop over configuration space cells; for each
  // config-space cell inner loop walks over velocity space
  gkyl_range_iter_init(&conf_iter, &conf_rng);
  while (gkyl_range_iter_next(&conf_iter)) {
    long midx = gkyl_range_idx(&conf_rng, conf_iter.idx);

    const double *M0_d = gkyl_array_cfetch(M0, midx);
    const double *M1i_d = gkyl_array_cfetch(M1i, midx);
    const double *M2_d = gkyl_array_cfetch(M2, midx);
    
    // compute primitive moments at quadrature nodes
    for (int n=0; n<tot_conf_quad; ++n) {
      const double *b_ord = gkyl_array_cfetch(up->conf_basis_at_ords, n);

      // number density
      num[n] = 0.0;
      for (int k=0; k<num_conf_basis; ++k)
        num[n] += M0_d[k]*b_ord[k];

      // velocity vector
      for (int d=0; d<vdim; ++d) {
        double M1i_n = 0.0;
        for (int k=0; k<num_conf_basis; ++k)
          M1i_n += M1i_d[num_conf_basis*d+k]*b_ord[k];
        vel[n][d] = M1i_n/num[n];
      }

      // vth2
      double M2_n = 0.0;
      for (int k=0; k<num_conf_basis; ++k)
        M2_n += M2_d[k]*b_ord[k];

      double v2 = 0.0; // vel^2
      for (int d=0; d<vdim; ++d) v2 += vel[n][d]*vel[n][d];
      vth2[n] = (M2_n - num[n]*v2)/(num[n]*vdim);
    }

    // inner loop over velocity space
    gkyl_range_deflate(&vel_rng, &phase_rng, rem_dir, conf_iter.idx);
    gkyl_range_iter_no_split_init(&vel_iter, &vel_rng);
    while (gkyl_range_iter_next(&vel_iter)) {
      
      copy_idx_arrays(conf_rng.ndim, phase_rng.ndim, conf_iter.idx, vel_iter.idx, pidx);
      gkyl_rect_grid_cell_center(&up->grid, pidx, xc);

      struct gkyl_range_iter qiter;
      // compute Maxwellian at phase-space quadrature nodes
      gkyl_range_iter_init(&qiter, &phase_qrange);
      while (gkyl_range_iter_next(&qiter)) {

        long cqidx = gkyl_range_idx(&conf_qrange, qiter.idx);
        double nvth2_q = num[cqidx]/pow(2*GKYL_PI*vth2[cqidx], vdim/2.0);

        long pqidx = gkyl_range_idx(&phase_qrange, qiter.idx);

        comp_to_phys(pdim, gkyl_array_cfetch(up->ordinates, pqidx),
          up->grid.dx, xc, xmu);

        double efact = 0.0;        
        for (int d=0; d<vdim; ++d)
          efact += (vel[cqidx][d]-xmu[cdim+d])*(vel[cqidx][d]-xmu[cdim+d]);

        double *fq = gkyl_array_fetch(fun_at_ords, pqidx);
        fq[0] = nvth2_q*exp(-efact/(2*vth2[cqidx]));
      }

      // compute expansion coefficients of Maxwellian on basis
      long lidx = gkyl_range_idx(&vel_rng, vel_iter.idx);
      proj_on_basis(up, fun_at_ords, gkyl_array_fetch(fmax, lidx));
    }
  }
  
  gkyl_array_release(fun_at_ords);
}

void
gkyl_proj_maxwellian_on_basis_release(gkyl_proj_maxwellian_on_basis* up)
{
  gkyl_array_release(up->ordinates);
  gkyl_array_release(up->weights);
  gkyl_array_release(up->basis_at_ords);
  gkyl_array_release(up->conf_ordinates);
  gkyl_array_release(up->conf_weights);
  gkyl_array_release(up->conf_basis_at_ords);
  gkyl_free(up);
}