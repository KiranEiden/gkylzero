#include <string.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_const.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_vlasov.h>
#include <gkyl_mom_vlasov_sr.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_proj_MJ_on_basis.h>
#include <gkyl_range.h>

struct gkyl_proj_MJ_on_basis {
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

  struct gkyl_dg_bin_op_mem *mem_for_div_op; // memory for weak division
  double mass; // used for init routine
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

gkyl_proj_MJ_on_basis*
gkyl_proj_MJ_on_basis_new(
  const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis,
  int num_quad, double mass)
{
  gkyl_proj_MJ_on_basis *up = gkyl_malloc(sizeof(gkyl_proj_MJ_on_basis));

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
proj_on_basis(const gkyl_proj_MJ_on_basis *up, const struct gkyl_array *fun_at_ords, double* f)
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
gkyl_proj_MJ_on_basis_fluid_stationary_frame_mom(const gkyl_proj_MJ_on_basis *up,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng,
  const struct gkyl_array *num_fluid_frame, const struct gkyl_array *vel_fluid_frame, const struct gkyl_array *T_fluid_frame,
  struct gkyl_array *f_MJ)
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
  for (int d=0; d<conf_rng->ndim; ++d) rem_dir[d] = 1;

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];
  double num[tot_conf_quad], vel[tot_conf_quad][vdim], T[tot_conf_quad];

  // outer loop over configuration space cells; for each
  // config-space cell inner loop walks over velocity space
  gkyl_range_iter_init(&conf_iter, conf_rng);
  while (gkyl_range_iter_next(&conf_iter)) {
    long midx = gkyl_range_idx(conf_rng, conf_iter.idx);

    const double *num_fluid_frame_d = gkyl_array_cfetch(num_fluid_frame, midx);
    const double *vel_fluid_frame_d = gkyl_array_cfetch(vel_fluid_frame, midx);
    const double *T_fluid_frame_d = gkyl_array_cfetch(T_fluid_frame, midx);

    // Sum over basis for given primative moments n,vector(v),T in the flow frame
    for (int n=0; n<tot_conf_quad; ++n) {
      const double *b_ord = gkyl_array_cfetch(up->conf_basis_at_ords, n);

      // number density
      num[n] = 0.0;
      for (int k=0; k<num_conf_basis; ++k)
        num[n] += num_fluid_frame_d[k]*b_ord[k];

      // velocity vector
      for (int d=0; d<vdim; ++d) {
        double vel_fluid_frame_n = 0.0;
        for (int k=0; k<num_conf_basis; ++k)
          vel_fluid_frame_n += vel_fluid_frame_d[num_conf_basis*d+k]*b_ord[k];
        vel[n][d] = vel_fluid_frame_n;
        //vel[n][d] = vel_fluid_frame_n/num[n];
      }

      // vth2
      double T_fluid_frame_n = 0.0;
      for (int k=0; k<num_conf_basis; ++k)
        T_fluid_frame_n += T_fluid_frame_d[k]*b_ord[k];


      //Originally: Not needed for vth2 calc, we assume (for MJ) that T_fluid_frame_n def is kT/mc^2
      //double v2 = 0.0; // vel^2
      //for (int d=0; d<vdim; ++d) v2 += vel[n][d]*vel[n][d];
      //vth2[n] = (T_fluid_frame_n - num[n]*v2)/(num[n]*vdim);

      // Using new def*** vth -> T (MJ change)
      T[n] = T_fluid_frame_n; // Change to P = <nT> moment
    }

    // inner loop over velocity space
    gkyl_range_deflate(&vel_rng, phase_rng, rem_dir, conf_iter.idx);
    gkyl_range_iter_no_split_init(&vel_iter, &vel_rng);
    while (gkyl_range_iter_next(&vel_iter)) {

      copy_idx_arrays(conf_rng->ndim, phase_rng->ndim, conf_iter.idx, vel_iter.idx, pidx);
      gkyl_rect_grid_cell_center(&up->grid, pidx, xc);

      struct gkyl_range_iter qiter;
      // compute MJ at phase-space quadrature nodes
      gkyl_range_iter_init(&qiter, &phase_qrange);
      while (gkyl_range_iter_next(&qiter)) {

        long cqidx = gkyl_range_idx(&conf_qrange, qiter.idx);
        //double nvth2_q = num[cqidx]/pow(2*GKYL_PI*vth2[cqidx], vdim/2.0);
	double kb = 1.0; // Update these later using the input data
	double mass_rest_frame = 1.0;
	double c = 1.0;
	double Theta = kb*T[cqidx]/(mass_rest_frame*c*c); // T = vth2[cqidx]; (?) - Need to re-write the moments
	double norm = num[cqidx] * (1.0/(4.0*GKYL_PI*mass_rest_frame*mass_rest_frame*mass_rest_frame*c*c*c*Theta)) * (sqrt(2*Theta/GKYL_PI));
  //Normalizaton needs to be handled here *

        long pqidx = gkyl_range_idx(&phase_qrange, qiter.idx);

        comp_to_phys(pdim, gkyl_array_cfetch(up->ordinates, pqidx),
          up->grid.dx, xc, xmu);

        //double efact = 0.0;
        //for (int d=0; d<vdim; ++d)
          //efact += (vel[cqidx][d]-xmu[cdim+d])*(vel[cqidx][d]-xmu[cdim+d]);

	double uu = 0.0;
	double vu = 0.0;
	double vv = 0.0;
	for (int d=0; d<vdim; ++d){
	   vv += (vel[cqidx][d]*vel[cqidx][d])/(c*c);
	   vu += (vel[cqidx][d]*xmu[cdim+d])/(c*c);
	   uu += (xmu[cdim+d]*xmu[cdim+d])/(c*c);
	}
	double gamma_shifted = 0.0;
	gamma_shifted = 1/sqrt(1-vv); // vv cannot be greater than 1, do we want to check for this?

        double *fq = gkyl_array_fetch(fun_at_ords, pqidx);
        fq[0] = norm*exp( (1.0/Theta) - (gamma_shifted/Theta)*(sqrt(1+uu) - vu) );
	//Takes the expanded K2 Bessel function to avoid calculating it for the normalization
	//fq[0] = nvth2_q*exp(-efact/(2*vth2[cqidx]));
      }

      // compute expansion coefficients of MJ on basis
      long lidx = gkyl_range_idx(&vel_rng, vel_iter.idx);
      proj_on_basis(up, fun_at_ords, gkyl_array_fetch(f_MJ, lidx));
    }
  }

      // To do list:
      // *** Need to compute: the integral N_Unnorm = int(f*dv)
        // Can we call the sr_mom_M0 with f from here to get N_Unnorm?
        //^^ Yes, see steps below
      // *** Need to compute: normaliztion*f where normalization = num_fluid_frame/N_Unnorm
      // *** Need to pull constants from species/field structure

      //Outside loop, but same function, calculate the normalization:
      //1. Call mom_calc advance (call on f_MJ) f_MJ is the basis expansion (unnormalized)
        // Returns unnormalied density
        // a. Start with: gkyl_mom_calc_new

      //struct gkyl_basis confbasis = *conf_basis;
      //struct gkyl_basis phasebasis = *phase_basis;
      //struct gkyl_mom_type *vmM0_t = gkyl_mom_vlasov_sr_new(&confbasis, &phasebasis, &vel_rng, "M0");

      struct gkyl_basis basis, confBasis;
      int poly_order = up->num_quad - 1;
      printf("(From line 311, Proj_MJ): poly_order %d, pdim %d, cdim %d\n",poly_order,up->pdim,up->cdim);
      gkyl_cart_modal_serendip(&basis, up->pdim, poly_order);
      gkyl_cart_modal_serendip(&confBasis, up->cdim, poly_order);

      // --> Uncomment this block for my attempt at normalization
      //Make a temporary basis for phase_basis
      //struct gkyl_mom_type *vmM0_t = gkyl_mom_vlasov_sr_new(up->conf_basis_at_ords, up->basis_at_ords, &vel_rng, "M0");
      //struct gkyl_mom_type *vmM0_t = gkyl_mom_vlasov_sr_new(&confBasis, &basis, &vel_rng, "M0");
      struct gkyl_mom_type *vmM0_t = gkyl_mom_vlasov_new(&confBasis, &basis, "M0");

      gkyl_mom_calc *m0calc = gkyl_mom_calc_new(&up->grid, vmM0_t);
        // b. Then call the advance method: gkyl_mom_calc_advance
      //struct gkyl_array* m0 = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, conf_rng->volume);
      //struct gkyl_array* norm_factor = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, conf_rng->volume);
      struct gkyl_array* m0 = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, num_fluid_frame->size);
      struct gkyl_array* norm_factor = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, num_fluid_frame->size);
      gkyl_mom_calc_advance(m0calc, phase_rng, conf_rng, f_MJ, m0);
        // c. Release the memory with: gkyl_mom_calc_new_release ( see the end)

      //2. Call gkyl_dg_div_op, we want to find the factor to multiply by: norm_factor = num_fluid_frame/m0 (actual over calculated)
      struct gkyl_dg_bin_op_mem *mem_for_div_op = gkyl_dg_bin_op_mem_new(f_MJ->size, up->num_phase_basis);
      //gkyl_dg_div_op(mem_for_div_op, confBasis, 0, num_fluid_frame, 0, m0, 0, norm_factor); //
      // g_bar = h/f = g
      // gkyl_dg_div_op(mem, basis, 0, g_bar_cu, 0, h_cu, 0, distf_cu);
      gkyl_dg_div_op_range(mem_for_div_op, confBasis, 0, norm_factor, 0, num_fluid_frame, 0, m0, *conf_rng);



      //3. Call dg_bin_op_comp_phase_multi (takes factor from step 2. and f_MJ)
      //How do I select part of the velocity space at once to do this multiplication??
      // gkyl_dg_mul_op(confBasis, 0, f_MJ, 0, norm_factor, 0, f_MJ); // &num[midx] changed to --> local_num
      // w = cfield*f
      //gkyl_dg_mul_conf_phase_op_range(&cbasis, &basis, w_bar, cfield, distf, &arr_crange, &arr_range);
      gkyl_dg_mul_conf_phase_op_range(&confBasis, &basis, f_MJ, norm_factor, f_MJ, conf_rng, phase_rng);
      //printf("(From line 341, Proj_MJ): norm_factor element 0 %e\n",&norm_factor->data[0]);
      //Print some contents of norm_factor


      // Extra: don't forget to free memory
      gkyl_array_release(m0);
      gkyl_array_release(norm_factor);
      gkyl_mom_type_release(vmM0_t);
      gkyl_mom_calc_release(m0calc);
      gkyl_dg_bin_op_mem_release(mem_for_div_op);


      //Outstanding questions
      //Does gkyl_dg_mult_op need to be moved up or down a loop (is it all calc. simult)?
      //Correct use of confbasis and phasebasis?
      //Does assert() call abort()?
        // Why would there be an assert for:
        // Assertion failed: (cbasis->poly_order == pbasis->poly_order), function gkyl_mom_vlasov_sr_new, file zero/mom_vlasov_sr.c, line 39.

  gkyl_array_release(fun_at_ords);
}

void
gkyl_proj_MJ_on_basis_release(gkyl_proj_MJ_on_basis* up)
{
  gkyl_array_release(up->ordinates);
  gkyl_array_release(up->weights);
  gkyl_array_release(up->basis_at_ords);
  gkyl_array_release(up->conf_ordinates);
  gkyl_array_release(up->conf_weights);
  gkyl_array_release(up->conf_basis_at_ords);
  gkyl_free(up);
}
