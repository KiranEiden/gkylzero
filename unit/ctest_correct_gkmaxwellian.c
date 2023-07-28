#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_correct_gkmaxwellian.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_gyrokinetic.h>
#include <gkyl_proj_maxwellian_on_basis.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <math.h>

// Allocate cu_dev array
static struct gkyl_array*
mkarr_cu(long nc, long size, bool use_gpu)
{
  if (use_gpu) {
    struct gkyl_array* a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
    return a;
  } else {
    struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
    return a;
  }
}

// Create ghost ranges
struct skin_ghost_ranges
{
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];

  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];
};

// Create ghost and skin sub-ranges given a parent range
static void
skin_ghost_ranges_init(struct skin_ghost_ranges *sgr,
                       const struct gkyl_range *parent, const int *ghost)
{
  int ndim = parent->ndim;

  for (int d = 0; d < ndim; ++d){
    gkyl_skin_ghost_ranges(&sgr->lower_skin[d], &sgr->lower_ghost[d],
                           d, GKYL_LOWER_EDGE, parent, ghost);
    gkyl_skin_ghost_ranges(&sgr->upper_skin[d], &sgr->upper_ghost[d],
                           d, GKYL_UPPER_EDGE, parent, ghost);
  }
}

void eval_jacob_tot(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0;
}
void eval_M0(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0;
}
void eval_M1(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.5;
}
void eval_M2(double t, const double *xn, double *restrict fout, void *ctx)
{
  double T = 2.0;
  double x = xn[0];
  fout[0] = T;
}
void eval_vtsq(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double vtsq = 1.0;
  fout[0] = vtsq;
}
void eval_udrift_2v_gk(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.5;
}
void eval_M1_2v_gk(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double den[1], udrift[1];
  eval_M0(t, xn, den, ctx);
  eval_udrift_2v_gk(t, xn, udrift, ctx);
  fout[0] = den[0]*udrift[0];
}
void eval_M2_2v_gk(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double den[1], udrift[1], vtsq[1];
  eval_M0(t, xn, den, ctx);
  eval_udrift_2v_gk(t, xn, udrift, ctx);
  eval_vtsq(t, xn, vtsq, ctx);
  fout[0] = 3.*den[0]*vtsq[0] + den[0]*(udrift[0]*udrift[0]);
}
void eval_bmag_1x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = cos((2.*M_PI/(2.*2.*M_PI))*x);
}

void test_1x1v(int poly_order, bool use_gpu)
{
  double mass = 1.0;
  double err_max = 1e-14, iter_max = 50;
  double lower[] = {-0.5, -5.0}, upper[] = {0.5, 5.0};
  int cells[] = {2, 8};
  int vdim = 1, cdim = 1;
  int ndim = cdim + vdim;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  int confCells[] = {cells[0]};

  // Grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  // Basis functions
  struct gkyl_basis basis, confBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  // Configuration space range
  int confGhost[] = {1};
  struct gkyl_range confLocal, confLocal_ext;
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost;
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  // Phase space range
  int ghost[] = {confGhost[0], 0};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost;
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  // Create bmag and jacob_tot arrays
  struct gkyl_array *bmag_ho, *jacob_tot_ho;
  bmag_ho = mkarr_cu(confBasis.num_basis, confLocal_ext.volume, false);
  jacob_tot_ho = mkarr_cu(confBasis.num_basis, confLocal_ext.volume, false);
  gkyl_proj_on_basis *proj_bmag = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_bmag_1x, NULL);
  gkyl_proj_on_basis *proj_jac = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_jacob_tot, NULL);
  gkyl_proj_on_basis_advance(proj_bmag, 0.0, &confLocal, bmag_ho);
  gkyl_proj_on_basis_advance(proj_jac, 0.0, &confLocal, jacob_tot_ho);
  struct gkyl_array *bmag, *jacob_tot;
  if (use_gpu) { 
    // create device copies
    bmag = mkarr_cu(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    jacob_tot = mkarr_cu(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    // copy host array to device
    gkyl_array_copy(bmag, bmag_ho);
    gkyl_array_copy(jacob_tot, jacob_tot_ho);
  } else {
    bmag = bmag_ho;
    jacob_tot = jacob_tot_ho;
  }

  // Create correct moment arrays
  struct gkyl_array *m0_corr_ho, *m1_corr_ho, *m2_corr_ho;
  m0_corr_ho = mkarr_cu(confBasis.num_basis, confLocal_ext.volume, false);
  m1_corr_ho = mkarr_cu(confBasis.num_basis, confLocal_ext.volume, false);
  m2_corr_ho = mkarr_cu(confBasis.num_basis, confLocal_ext.volume, false);
  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_M0, NULL);
  gkyl_proj_on_basis *proj_m1 = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_M1, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_M2, NULL);
  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0_corr_ho);
  gkyl_proj_on_basis_advance(proj_m1, 0.0, &confLocal, m1_corr_ho);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2_corr_ho);
  struct gkyl_array *m0_corr, *m1_corr, *m2_corr;
  if (use_gpu) {
    m0_corr = mkarr_cu(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    m1_corr = mkarr_cu(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    m2_corr = mkarr_cu(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    gkyl_array_copy(m0_corr, m0_corr_ho);
    gkyl_array_copy(m1_corr, m1_corr_ho);
    gkyl_array_copy(m2_corr, m2_corr_ho);
  } else {
    m0_corr = m0_corr_ho;
    m1_corr = m1_corr_ho;
    m2_corr = m2_corr_ho;
  }
  // Write the input moments
  char fname_m0_in[1024], fname_m1_in[1024], fname_m2_in[1024];
  sprintf(fname_m0_in, "ctest_correct_maxwellian_%dx%dv_p%d_m0_in.gkyl", cdim, vdim, poly_order);
  sprintf(fname_m1_in, "ctest_correct_maxwellian_%dx%dv_p%d_m1_in.gkyl", cdim, vdim, poly_order);
  sprintf(fname_m2_in, "ctest_correct_maxwellian_%dx%dv_p%d_m2_in.gkyl", cdim, vdim, poly_order);
  gkyl_grid_sub_array_write(&confGrid, &confLocal, m0_corr, fname_m0_in);
  gkyl_grid_sub_array_write(&confGrid, &confLocal, m1_corr, fname_m1_in);
  gkyl_grid_sub_array_write(&confGrid, &confLocal, m2_corr, fname_m2_in);

  // Project the Maxwellian on basis
  // (1) proj_maxwellian expects the moments as a single array
  struct gkyl_array *moms_ho = mkarr_cu(3*confBasis.num_basis, confLocal_ext.volume, false);
  gkyl_array_set_offset(moms_ho, 1., m0_corr, 0*confBasis.num_basis);
  gkyl_array_set_offset(moms_ho, 1., m1_corr, 1*confBasis.num_basis);
  gkyl_array_set_offset(moms_ho, 1., m2_corr, 2*confBasis.num_basis);
  struct gkyl_array *moms;
  if (use_gpu)
  { // copy host array to device
    moms = mkarr_cu(3*confBasis.num_basis, confLocal_ext.volume, use_gpu);
    gkyl_array_copy(moms, moms_ho);
  }
  else
  {
    moms = moms_ho;
  }
  // (2) create distribution function array
  gkyl_proj_maxwellian_on_basis *proj_maxwellian = gkyl_proj_maxwellian_on_basis_new(&grid, &confBasis, &basis, poly_order+1, use_gpu);
  struct gkyl_array *fM = use_gpu? mkarr_cu(basis.num_basis, local_ext.volume, use_gpu)
                                 : mkarr_cu(basis.num_basis, local_ext.volume, false); 
  gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_maxwellian, &local, &confLocal, moms, bmag, jacob_tot, mass, fM);
  // Calculate the uncorrected moments
  // (1) create the calculators
  struct gkyl_mom_type *M0_t = gkyl_mom_gyrokinetic_new(&confBasis, &basis, &confLocal, mass, "M0", use_gpu);
  struct gkyl_mom_type *M1_t = gkyl_mom_gyrokinetic_new(&confBasis, &basis, &confLocal, mass, "M1", use_gpu);
  struct gkyl_mom_type *M2_t = gkyl_mom_gyrokinetic_new(&confBasis, &basis, &confLocal, mass, "M2", use_gpu);
  gkyl_gyrokinetic_set_bmag(M0_t, bmag);
  gkyl_gyrokinetic_set_bmag(M1_t, bmag);
  gkyl_gyrokinetic_set_bmag(M2_t, bmag);
  gkyl_mom_calc *m0calc = gkyl_mom_calc_new(&grid, M0_t, use_gpu);
  gkyl_mom_calc *m1calc = gkyl_mom_calc_new(&grid, M1_t, use_gpu);
  gkyl_mom_calc *m2calc = gkyl_mom_calc_new(&grid, M2_t, use_gpu);
  // (2) create moment arrays
  struct gkyl_array *m0_ho, *m1_ho, *m2_ho; 
  m0_ho = mkarr_cu(confBasis.num_basis, confLocal_ext.volume, false);
  m1_ho = mkarr_cu(confBasis.num_basis, confLocal_ext.volume, false);
  m2_ho = mkarr_cu(confBasis.num_basis, confLocal_ext.volume, false);
  struct gkyl_array *m0, *m1, *m2;
  if (use_gpu) {
    m0 = mkarr_cu(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    m1 = mkarr_cu(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    m2 = mkarr_cu(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    gkyl_mom_calc_advance(m0calc, &local, &confLocal, fM, m0);
    gkyl_mom_calc_advance(m1calc, &local, &confLocal, fM, m1);
    gkyl_mom_calc_advance(m2calc, &local, &confLocal, fM, m2);
  } else {
    gkyl_mom_calc_advance(m0calc, &local, &confLocal, fM, m0_ho);
    gkyl_mom_calc_advance(m1calc, &local, &confLocal, fM, m1_ho);
    gkyl_mom_calc_advance(m2calc, &local, &confLocal, fM, m2_ho);
    m0 = m0_ho;
    m1 = m1_ho;
    m2 = m2_ho;
  }    
  // Write the uncorrected moments  
  char fname_fM_ic[1024];
  sprintf(fname_fM_ic, "ctest_correct_maxwellian_%dx%dv_p%d.gkyl", cdim, vdim, poly_order);
  gkyl_grid_sub_array_write(&grid, &local, fM, fname_fM_ic);
  char fname_m0_ic[1024], fname_m1_ic[1024], fname_m2_ic[1024];
  sprintf(fname_m0_ic, "ctest_correct_maxwellian_%dx%dv_p%d_m0_ic.gkyl", cdim, vdim, poly_order);
  sprintf(fname_m1_ic, "ctest_correct_maxwellian_%dx%dv_p%d_m1_ic.gkyl", cdim, vdim, poly_order);
  sprintf(fname_m2_ic, "ctest_correct_maxwellian_%dx%dv_p%d_m2_ic.gkyl", cdim, vdim, poly_order);
  gkyl_grid_sub_array_write(&confGrid, &confLocal, m0, fname_m0_ic);
  gkyl_grid_sub_array_write(&confGrid, &confLocal, m1, fname_m1_ic);
  gkyl_grid_sub_array_write(&confGrid, &confLocal, m2, fname_m2_ic);

  // Create a Maxwellian with corrected moments
  gkyl_correct_gkmaxwellian *corr_max = gkyl_correct_gkmaxwellian_new(&grid, &confBasis, &basis, &confLocal, &confLocal_ext, bmag, mass, use_gpu);
  gkyl_correct_gkmaxwellian_fix(corr_max, fM, m0_corr, m1_corr, m2_corr, jacob_tot, bmag, mass, err_max, iter_max, &confLocal, &local, &confLocal_ext, poly_order, use_gpu);
  gkyl_correct_gkmaxwellian_release(corr_max);

  // Calculate the corrected moments
  gkyl_array_clear(m0, 0.0);
  gkyl_array_clear(m1, 0.0);
  gkyl_array_clear(m2, 0.0);
  gkyl_mom_calc_advance(m0calc, &local, &confLocal, fM, m0);
  gkyl_mom_calc_advance(m1calc, &local, &confLocal, fM, m1);
  gkyl_mom_calc_advance(m2calc, &local, &confLocal, fM, m2);
  // Write the output
  char fname_fM_corr[1024];
  sprintf(fname_fM_corr, "ctest_correct_maxwellian_%dx%dv_p%d_corr.gkyl", cdim, vdim, poly_order);
  gkyl_grid_sub_array_write(&grid, &local, fM, fname_fM_corr);
  char fname_m0_corr[1024], fname_m1_corr[1024], fname_m2_corr[1024];
  sprintf(fname_m0_corr, "ctest_correct_maxwellian_%dx%dv_p%d_m0_corr.gkyl", cdim, vdim, poly_order);
  sprintf(fname_m1_corr, "ctest_correct_maxwellian_%dx%dv_p%d_m1_corr.gkyl", cdim, vdim, poly_order);
  sprintf(fname_m2_corr, "ctest_correct_maxwellian_%dx%dv_p%d_m2_corr.gkyl", cdim, vdim, poly_order);
  gkyl_grid_sub_array_write(&confGrid, &confLocal, m0, fname_m0_corr);
  gkyl_grid_sub_array_write(&confGrid, &confLocal, m1, fname_m1_corr);
  gkyl_grid_sub_array_write(&confGrid, &confLocal, m2, fname_m2_corr);

  // For higher polynomial order basis functions
  /*
  if (poly_order == 2)
    for (int i = 0; i < basis.num_basis; ++i)
      TEST_CHECK(gkyl_compare_double(p2_vals[i], fv[i], 1e-12)) 
  */

  // Release memory for moment data object
  gkyl_array_release(m0_corr);
  gkyl_array_release(m1_corr);
  gkyl_array_release(m2_corr);
  gkyl_array_release(m0_corr_ho);
  gkyl_array_release(m1_corr_ho);
  gkyl_array_release(m2_corr_ho);
  gkyl_array_release(m0);
  gkyl_array_release(m1);
  gkyl_array_release(m2);
  gkyl_array_release(m0_ho);
  gkyl_array_release(m1_ho);
  gkyl_array_release(m2_ho);
  gkyl_array_release(moms_ho);
  gkyl_array_release(fM);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1);
  gkyl_proj_on_basis_release(proj_m2);
  gkyl_mom_calc_release(m0calc);
  gkyl_mom_calc_release(m1calc);
  gkyl_mom_calc_release(m2calc);
}

void test_1x2v(int poly_order, bool use_gpu)
{
  double mass = 1.0;
  double err_max = 1e-14, iter_max = 50;
  double lower[] = {-0.5, -5.0, 0.0}, upper[] = {0.5, 5.0, 5.0};
  int cells[] = {2, 8, 8};
  int vdim = 2, cdim = 1;
  int ndim = cdim + vdim;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  int confCells[] = {cells[0]};

  // Grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  // Basis functions
  struct gkyl_basis basis, confBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  // Configuration space range
  int confGhost[] = {1};
  struct gkyl_range confLocal, confLocal_ext;
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost;
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  // Phase space range
  int ghost[] = {confGhost[0], 0, 0};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost;
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  // Create bmag and jacob_tot arrays
  struct gkyl_array *bmag_ho, *jacob_tot_ho;
  bmag_ho = mkarr_cu(confBasis.num_basis, confLocal_ext.volume, false);
  jacob_tot_ho = mkarr_cu(confBasis.num_basis, confLocal_ext.volume, false);
  gkyl_proj_on_basis *proj_bmag = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_bmag_1x, NULL);
  gkyl_proj_on_basis *proj_jac = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_jacob_tot, NULL);
  gkyl_proj_on_basis_advance(proj_bmag, 0.0, &confLocal, bmag_ho);
  gkyl_proj_on_basis_advance(proj_jac, 0.0, &confLocal, jacob_tot_ho);
  struct gkyl_array *bmag, *jacob_tot;
  if (use_gpu) { 
    // create device copies
    bmag = mkarr_cu(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    jacob_tot = mkarr_cu(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    // copy host array to device
    gkyl_array_copy(bmag, bmag_ho);
    gkyl_array_copy(jacob_tot, jacob_tot_ho);
  } else {
    bmag = bmag_ho;
    jacob_tot = jacob_tot_ho;
  }

  // Create correct moment arrays
  struct gkyl_array *m0_corr_ho, *m1_corr_ho, *m2_corr_ho;
  m0_corr_ho = mkarr_cu(confBasis.num_basis, confLocal_ext.volume, false);
  m1_corr_ho = mkarr_cu(confBasis.num_basis, confLocal_ext.volume, false);
  m2_corr_ho = mkarr_cu(confBasis.num_basis, confLocal_ext.volume, false);
  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_M0, NULL);
  gkyl_proj_on_basis *proj_m1 = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_M1, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_M2, NULL);
  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0_corr_ho);
  gkyl_proj_on_basis_advance(proj_m1, 0.0, &confLocal, m1_corr_ho);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2_corr_ho);
  struct gkyl_array *m0_corr, *m1_corr, *m2_corr;
  if (use_gpu) {
    m0_corr = mkarr_cu(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    m1_corr = mkarr_cu(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    m2_corr = mkarr_cu(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    gkyl_array_copy(m0_corr, m0_corr_ho);
    gkyl_array_copy(m1_corr, m1_corr_ho);
    gkyl_array_copy(m2_corr, m2_corr_ho);
  } else {
    m0_corr = m0_corr_ho;
    m1_corr = m1_corr_ho;
    m2_corr = m2_corr_ho;
  }
  // Write the input moments
  char fname_m0_in[1024], fname_m1_in[1024], fname_m2_in[1024];
  sprintf(fname_m0_in, "ctest_correct_maxwellian_%dx%dv_p%d_m0_in.gkyl", cdim, vdim, poly_order);
  sprintf(fname_m1_in, "ctest_correct_maxwellian_%dx%dv_p%d_m1_in.gkyl", cdim, vdim, poly_order);
  sprintf(fname_m2_in, "ctest_correct_maxwellian_%dx%dv_p%d_m2_in.gkyl", cdim, vdim, poly_order);
  gkyl_grid_sub_array_write(&confGrid, &confLocal, m0_corr, fname_m0_in);
  gkyl_grid_sub_array_write(&confGrid, &confLocal, m1_corr, fname_m1_in);
  gkyl_grid_sub_array_write(&confGrid, &confLocal, m2_corr, fname_m2_in);

  // Project the Maxwellian on basis
  // (1) proj_maxwellian expects the moments as a single array
  struct gkyl_array *moms_ho = mkarr_cu(3*confBasis.num_basis, confLocal_ext.volume, false);
  gkyl_array_set_offset(moms_ho, 1., m0_corr, 0*confBasis.num_basis);
  gkyl_array_set_offset(moms_ho, 1., m1_corr, 1*confBasis.num_basis);
  gkyl_array_set_offset(moms_ho, 1., m2_corr, 2*confBasis.num_basis);
  struct gkyl_array *moms;
  if (use_gpu)
  { // copy host array to device
    moms = mkarr_cu(3*confBasis.num_basis, confLocal_ext.volume, use_gpu);
    gkyl_array_copy(moms, moms_ho);
  }
  else
  {
    moms = moms_ho;
  }
  // (2) create distribution function array
  gkyl_proj_maxwellian_on_basis *proj_maxwellian = gkyl_proj_maxwellian_on_basis_new(&grid, &confBasis, &basis, poly_order+1, use_gpu);
  struct gkyl_array *fM = use_gpu? mkarr_cu(basis.num_basis, local_ext.volume, use_gpu)
                                 : mkarr_cu(basis.num_basis, local_ext.volume, false); 
  gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_maxwellian, &local, &confLocal, moms, bmag, jacob_tot, mass, fM);
  // Calculate the uncorrected moments
  // (1) create the calculators
  struct gkyl_mom_type *M0_t = gkyl_mom_gyrokinetic_new(&confBasis, &basis, &confLocal, mass, "M0", use_gpu);
  struct gkyl_mom_type *M1_t = gkyl_mom_gyrokinetic_new(&confBasis, &basis, &confLocal, mass, "M1", use_gpu);
  struct gkyl_mom_type *M2_t = gkyl_mom_gyrokinetic_new(&confBasis, &basis, &confLocal, mass, "M2", use_gpu);
  gkyl_gyrokinetic_set_bmag(M0_t, bmag);
  gkyl_gyrokinetic_set_bmag(M1_t, bmag);
  gkyl_gyrokinetic_set_bmag(M2_t, bmag);
  gkyl_mom_calc *m0calc = gkyl_mom_calc_new(&grid, M0_t, use_gpu);
  gkyl_mom_calc *m1calc = gkyl_mom_calc_new(&grid, M1_t, use_gpu);
  gkyl_mom_calc *m2calc = gkyl_mom_calc_new(&grid, M2_t, use_gpu);
  // (2) create moment arrays
  struct gkyl_array *m0_ho, *m1_ho, *m2_ho; 
  m0_ho = mkarr_cu(confBasis.num_basis, confLocal_ext.volume, false);
  m1_ho = mkarr_cu(confBasis.num_basis, confLocal_ext.volume, false);
  m2_ho = mkarr_cu(confBasis.num_basis, confLocal_ext.volume, false);
  struct gkyl_array *m0, *m1, *m2;
  if (use_gpu) {
    m0 = mkarr_cu(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    m1 = mkarr_cu(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    m2 = mkarr_cu(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    gkyl_mom_calc_advance(m0calc, &local, &confLocal, fM, m0);
    gkyl_mom_calc_advance(m1calc, &local, &confLocal, fM, m1);
    gkyl_mom_calc_advance(m2calc, &local, &confLocal, fM, m2);
  } else {
    gkyl_mom_calc_advance(m0calc, &local, &confLocal, fM, m0_ho);
    gkyl_mom_calc_advance(m1calc, &local, &confLocal, fM, m1_ho);
    gkyl_mom_calc_advance(m2calc, &local, &confLocal, fM, m2_ho);
    m0 = m0_ho;
    m1 = m1_ho;
    m2 = m2_ho;
  }    
  // Write the uncorrected moments  
  char fname_fM_ic[1024];
  sprintf(fname_fM_ic, "ctest_correct_maxwellian_%dx%dv_p%d.gkyl", cdim, vdim, poly_order);
  gkyl_grid_sub_array_write(&grid, &local, fM, fname_fM_ic);
  char fname_m0_ic[1024], fname_m1_ic[1024], fname_m2_ic[1024];
  sprintf(fname_m0_ic, "ctest_correct_maxwellian_%dx%dv_p%d_m0_ic.gkyl", cdim, vdim, poly_order);
  sprintf(fname_m1_ic, "ctest_correct_maxwellian_%dx%dv_p%d_m1_ic.gkyl", cdim, vdim, poly_order);
  sprintf(fname_m2_ic, "ctest_correct_maxwellian_%dx%dv_p%d_m2_ic.gkyl", cdim, vdim, poly_order);
  gkyl_grid_sub_array_write(&confGrid, &confLocal, m0, fname_m0_ic);
  gkyl_grid_sub_array_write(&confGrid, &confLocal, m1, fname_m1_ic);
  gkyl_grid_sub_array_write(&confGrid, &confLocal, m2, fname_m2_ic);

  // Create a Maxwellian with corrected moments
  gkyl_correct_gkmaxwellian *corr_max = gkyl_correct_gkmaxwellian_new(&grid, &confBasis, &basis, &confLocal, &confLocal_ext, bmag, mass, use_gpu);
  gkyl_correct_gkmaxwellian_fix(corr_max, fM, m0_corr, m1_corr, m2_corr, jacob_tot, bmag, mass, err_max, iter_max, &confLocal, &local, &confLocal_ext, poly_order, use_gpu);
  gkyl_correct_gkmaxwellian_release(corr_max);

  // Calculate the corrected moments
  gkyl_array_clear(m0, 0.0);
  gkyl_array_clear(m1, 0.0);
  gkyl_array_clear(m2, 0.0);
  gkyl_mom_calc_advance(m0calc, &local, &confLocal, fM, m0);
  gkyl_mom_calc_advance(m1calc, &local, &confLocal, fM, m1);
  gkyl_mom_calc_advance(m2calc, &local, &confLocal, fM, m2);
  // Write the output
  char fname_fM_corr[1024];
  sprintf(fname_fM_corr, "ctest_correct_maxwellian_%dx%dv_p%d_corr.gkyl", cdim, vdim, poly_order);
  gkyl_grid_sub_array_write(&grid, &local, fM, fname_fM_corr);
  char fname_m0_corr[1024], fname_m1_corr[1024], fname_m2_corr[1024];
  sprintf(fname_m0_corr, "ctest_correct_maxwellian_%dx%dv_p%d_m0_corr.gkyl", cdim, vdim, poly_order);
  sprintf(fname_m1_corr, "ctest_correct_maxwellian_%dx%dv_p%d_m1_corr.gkyl", cdim, vdim, poly_order);
  sprintf(fname_m2_corr, "ctest_correct_maxwellian_%dx%dv_p%d_m2_corr.gkyl", cdim, vdim, poly_order);
  gkyl_grid_sub_array_write(&confGrid, &confLocal, m0, fname_m0_corr);
  gkyl_grid_sub_array_write(&confGrid, &confLocal, m1, fname_m1_corr);
  gkyl_grid_sub_array_write(&confGrid, &confLocal, m2, fname_m2_corr);

  // For higher polynomial order basis functions
  /*
  if (poly_order == 2)
    for (int i = 0; i < basis.num_basis; ++i)
      TEST_CHECK(gkyl_compare_double(p2_vals[i], fv[i], 1e-12)) 
  */

  // Release memory for moment data object
  gkyl_array_release(m0_corr);
  gkyl_array_release(m1_corr);
  gkyl_array_release(m2_corr);
  gkyl_array_release(m0_corr_ho);
  gkyl_array_release(m1_corr_ho);
  gkyl_array_release(m2_corr_ho);
  gkyl_array_release(m0);
  gkyl_array_release(m1);
  gkyl_array_release(m2);
  gkyl_array_release(m0_ho);
  gkyl_array_release(m1_ho);
  gkyl_array_release(m2_ho);
  gkyl_array_release(moms_ho);
  gkyl_array_release(fM);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1);
  gkyl_proj_on_basis_release(proj_m2);
  gkyl_mom_calc_release(m0calc);
  gkyl_mom_calc_release(m1calc);
  gkyl_mom_calc_release(m2calc);
}

// Run the test
void test_1x1v_p1() {test_1x1v(1, false);}
//void test_1x1v_p2() {test_1x1v(2, false);}
void test_1x2v_p1() {test_1x2v(1, false);}

TEST_LIST = {
  {"test_1x1v_p1", test_1x1v_p1},
//  {"test_1x1v_p2", test_1x1v_p2},
  {"test_1x2v_p1", test_1x2v_p1},
  {NULL, NULL},
};
