#include <acutest.h>

#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wave_geom_priv.h>

#include <math.h>

static void
my_nomapc2p(double t, const double *xc, double *xp, void *ctx)
{
  int *ndim = ctx;
  for (int i=0; i<(*ndim); ++i) xp[i] = xc[i];
}

void
test_wv_geom_1d_1()
{
  int ndim = 1;
  double lower[] = {0.0}, upper[] = {1.0};
  int cells[] = {10};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // create range
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  struct gkyl_wave_geom *wg = gkyl_wave_geom_new(&grid, &arr_range, my_nomapc2p, &ndim, false);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &arr_range);
  
  while (gkyl_range_iter_next(&iter)) {
    const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wg, iter.idx);

    TEST_CHECK( gkyl_compare_double( cg->kappa, 1.0, 1e-15) );
    TEST_CHECK( cg->lenr[0] == 1.0 );

    TEST_CHECK( cg->norm_cov[0][0] == 1.0 );
    TEST_CHECK( cg->tau1_cov[0][1] == 1.0 );
    TEST_CHECK( cg->tau2_cov[0][2] == 1.0 );
  }

  gkyl_wave_geom_release(wg);
}

static void
mapc2p(double t, const double *xc, double *xp, void *ctx)
{
  // quadratic mapping
  int *ndim = ctx;
  for (int i=0; i<(*ndim); ++i) xp[i] = xc[i]*xc[i];
}

void
test_wv_geom_1d_2()
{
  int ndim = 1;
  double lower[] = {0.0}, upper[] = {1.0};
  int cells[] = {2};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // create range
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range range, ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  struct gkyl_wave_geom *wg = gkyl_wave_geom_new(&grid, &range, mapc2p, &ndim, false);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);

  // cell 1
  do {
    const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wg, (int[]) { range.lower[0]+0 } );
    TEST_CHECK( gkyl_compare_double( cg->kappa, 0.5*0.5/0.5, 1e-15) );
    TEST_CHECK( cg->lenr[0] == 1.0 );
    TEST_CHECK( cg->norm_cov[0][0] == 1.0 );
    TEST_CHECK( cg->tau1_cov[0][1] == 1.0 );
    TEST_CHECK( cg->tau2_cov[0][2] == 1.0 );
  } while (0);

  // cell 2
  do {
    const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wg, (int[]) { range.lower[0]+1 } );
    TEST_CHECK( gkyl_compare_double( cg->kappa, (1-0.5*0.5)/0.5, 1e-15) );
    TEST_CHECK( cg->lenr[0] == 1.0 );
    TEST_CHECK( cg->norm_cov[0][0] == 1.0 );
    TEST_CHECK( cg->tau1_cov[0][1] == 1.0 );
    TEST_CHECK( cg->tau2_cov[0][2] == 1.0 );
  } while (0);  

  gkyl_wave_geom_release(wg);
}

void
test_wv_geom_2d_1()
{
  int ndim = 2;
  double lower[] = {0.0, 0.0}, upper[] = {1.0, 1.0};
  int cells[] = {2, 2};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // create range
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range range, ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  struct gkyl_wave_geom *wg = gkyl_wave_geom_new(&grid, &range, my_nomapc2p, &ndim, false);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  
  while (gkyl_range_iter_next(&iter)) {
    const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wg, iter.idx);

    TEST_CHECK( gkyl_compare_double( cg->kappa, 1.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->lenr[0], 1.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->lenr[1], 1.0, 1e-15) );

    // normal to left face is ex
    TEST_CHECK( gkyl_compare_double( cg->norm_cov[0][0], 1.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->norm_cov[0][1], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->norm_cov[0][2], 0.0, 1e-15) );

    // tangent1 to left face is ey
    TEST_CHECK( gkyl_compare_double( cg->tau1_cov[0][0], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau1_cov[0][1], 1.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau1_cov[0][2], 0.0, 1e-15) );

    // tangent2 to left face is ez
    TEST_CHECK( gkyl_compare_double( cg->tau2_cov[0][0], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau2_cov[0][1], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau2_cov[0][2], 1.0, 1e-15) );

    // normal to bottom face is ey
    TEST_CHECK( gkyl_compare_double( cg->norm_cov[1][0], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->norm_cov[1][1], 1.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->norm_cov[1][0], 0.0, 1e-15) );

    // tangent1 to bottom face is ex
    TEST_CHECK( gkyl_compare_double( cg->tau1_cov[1][0], 1.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau1_cov[1][1], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau1_cov[1][2], 0.0, 1e-15) );

    // tangent2 to bottom face is ez
    TEST_CHECK( gkyl_compare_double( cg->tau2_cov[1][0], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau2_cov[1][1], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau2_cov[1][2], -1.0, 1e-15) );
  }

  gkyl_wave_geom_release(wg);
}

static void
mapc2p_2d(double t, const double *xc, double *xp, void *ctx)
{
  double x = xc[0], y = xc[1];

  xp[0] = 0.375*((x+1.0)*y+x+1.0)-0.25*((x+1.0)*y-1.0*x-1.0)-0.125*((x-1.0)*y+x-1.0);
  xp[1] = 0.25*((x+1.0)*y+x+1.0)-0.25*((x-1.0)*y+x-1.0);
}

void
test_wv_geom_2d_2()
{
  int ndim = 2;
  double lower[] = {-1.0, -1.0}, upper[] = {1.0, 1.0};
  int cells[] = {1, 1};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // create range
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range range, ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  struct gkyl_wave_geom *wg = gkyl_wave_geom_new(&grid, &range, mapc2p_2d, &ndim, false);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  
  while (gkyl_range_iter_next(&iter)) {
    const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wg, iter.idx);

    TEST_CHECK( gkyl_compare_double( cg->kappa, 1.0/4.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->lenr[0], 1.118033988749895/2.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->lenr[1], 1.0/2.0, 1e-15) );

    // normal to left face is ex
    TEST_CHECK( gkyl_compare_double( cg->norm_cov[0][0], 0.8944271909999159, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->norm_cov[0][1], -0.4472135954999579, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->norm_cov[0][2], 0.0, 1e-15) );

    // tangent1 to left face
    TEST_CHECK( gkyl_compare_double( cg->tau1_cov[0][0], 0.4472135954999579, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau1_cov[0][1], 0.8944271909999159, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau1_cov[0][2], 0.0, 1e-15) );

    // tangent2 to left face is ez
    TEST_CHECK( gkyl_compare_double( cg->tau2_cov[0][0], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau2_cov[0][1], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau2_cov[0][2], 1.0, 1e-15) );

    // normal to bottom face is ey
    TEST_CHECK( gkyl_compare_double( cg->norm_cov[1][0], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->norm_cov[1][1], 1.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->norm_cov[1][0], 0.0, 1e-15) );

    // tangent1 to bottom face is ex
    TEST_CHECK( gkyl_compare_double( cg->tau1_cov[1][0], 1.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau1_cov[1][1], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau1_cov[1][2], 0.0, 1e-15) );

    // tangent2 to bottom face is ez
    TEST_CHECK( gkyl_compare_double( cg->tau2_cov[1][0], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau2_cov[1][1], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau2_cov[1][2], -1.0, 1e-15) );
  }

  gkyl_wave_geom_release(wg);
}

// map (r,theta) -> (x,y)
void
mapc2p_pol_local(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  double r = xc[0], th = xc[1];
  xp[0] = r*cos(th); xp[1] = r*sin(th);
}

void
test_wv_geom_2d_3()
{
  int ndim = 2;
  double r_inn = 0.25, r_out = 1.25;
  double phi_max = M_PI / 2.;
  double lower[] = {r_inn, 0}, upper[] = {r_out, phi_max};
  int cells[] = {1, 1};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  double area = 0.5*(r_out*r_out-r_inn*r_inn);
  double edge_inn = sqrt(2) * r_inn;
  double area_c = (r_out - r_inn) * phi_max;

  // create range
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range range, ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  struct gkyl_wave_geom *wg = gkyl_wave_geom_new(&grid, &range, mapc2p_pol_local, &ndim, false);
  // Also check that cylindrical coordinate flag works for 2D polar
  struct gkyl_wave_geom *wg_fl = gkyl_wave_geom_from_coord_flag(&grid, &range, WAVE_COORD_CYL, &ndim, false);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  
  while (gkyl_range_iter_next(&iter)) {
    const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wg, iter.idx);
    const struct gkyl_wave_cell_geom *cg_fl = gkyl_wave_geom_get(wg_fl, iter.idx);

    TEST_CHECK( gkyl_compare_double( cg->kappa, area / area_c, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->lenr[0], edge_inn / phi_max, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->lenr[1], 1, 1e-15) );

    // normal to left face has phi angle 45 deg
    TEST_CHECK( gkyl_compare_double( cg->norm_cov[0][0], 1/sqrt(2), 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->norm_cov[0][1], 1/sqrt(2), 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->norm_cov[0][2], 0.0, 1e-15) );

    // tangent1 to left face has phi angle 135 deg
    TEST_CHECK( gkyl_compare_double( cg->tau1_cov[0][0], -1/sqrt(2), 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau1_cov[0][1], 1/sqrt(2), 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau1_cov[0][2], 0.0, 1e-15) );

    // tangent2 to left face is ez
    TEST_CHECK( gkyl_compare_double( cg->tau2_cov[0][0], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau2_cov[0][1], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau2_cov[0][2], 1.0, 1e-15) );

    // normal to bottom face is ey
    TEST_CHECK( gkyl_compare_double( cg->norm_cov[1][0], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->norm_cov[1][1], 1.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->norm_cov[1][2], 0.0, 1e-15) );

    // tangent1 to bottom face is ex
    TEST_CHECK( gkyl_compare_double( cg->tau1_cov[1][0], 1.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau1_cov[1][1], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau1_cov[1][2], 0.0, 1e-15) );

    // tangent2 to bottom face is -ez
    TEST_CHECK( gkyl_compare_double( cg->tau2_cov[1][0], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau2_cov[1][1], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau2_cov[1][2], -1.0, 1e-15) );
    
    /********************************
     * Checks with cylindrical flag *
     ********************************/
    
    TEST_CHECK( gkyl_compare_double( cg_fl->kappa, area / area_c, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg_fl->lenr[0], edge_inn / phi_max, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg_fl->lenr[1], 1, 1e-15) );

    // normal to left face is in radial direction
    TEST_CHECK( gkyl_compare_double( cg_fl->norm_cov[0][0], 1.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg_fl->norm_cov[0][1], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg_fl->norm_cov[0][2], 0.0, 1e-15) );

    // tangent1 to left face is in azimuthal direction
    TEST_CHECK( gkyl_compare_double( cg_fl->tau1_cov[0][0], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg_fl->tau1_cov[0][1], 1.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg_fl->tau1_cov[0][2], 0.0, 1e-15) );

    // tangent2 to left face is ez
    TEST_CHECK( gkyl_compare_double( cg_fl->tau2_cov[0][0], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg_fl->tau2_cov[0][1], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg_fl->tau2_cov[0][2], 1.0, 1e-15) );
    
    // normal to bottom face is ey
    TEST_CHECK( gkyl_compare_double( cg_fl->norm_cov[1][0], 1/sqrt(2), 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg_fl->norm_cov[1][1], 1/sqrt(2), 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg_fl->norm_cov[1][2], 0.0, 1e-15) );

    // tangent1 to bottom face is ex
    TEST_CHECK( gkyl_compare_double( cg_fl->tau1_cov[1][0], 1/sqrt(2), 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg_fl->tau1_cov[1][1], -1/sqrt(2), 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg_fl->tau1_cov[1][2], 0.0, 1e-15) );

    // tangent2 to bottom face is -ez
    TEST_CHECK( gkyl_compare_double( cg_fl->tau2_cov[1][0], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg_fl->tau2_cov[1][1], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg_fl->tau2_cov[1][2], -1.0, 1e-15) );
  }

  gkyl_wave_geom_release(wg);
}

void
test_wv_geom_3d_1()
{
  int ndim = 3;
  double lower[] = {0.0, 0.0, 0.0}, upper[] = {1.0, 1.0, 1.0};
  int cells[] = {2, 2, 2};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // create range
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range range, ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  struct gkyl_wave_geom *wg = gkyl_wave_geom_new(&grid, &range, my_nomapc2p, &ndim, false);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  
  while (gkyl_range_iter_next(&iter)) {
    const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wg, iter.idx);

    TEST_CHECK( gkyl_compare_double( cg->kappa, 1.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->lenr[0], 1.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->lenr[1], 1.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->lenr[2], 1.0, 1e-15) );

    // normal to lower-x face is ex
    TEST_CHECK( gkyl_compare_double( cg->norm_cov[0][0], 1.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->norm_cov[0][1], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->norm_cov[0][2], 0.0, 1e-15) );

    // tangent1 to lower-x face is ey
    TEST_CHECK( gkyl_compare_double( cg->tau1_cov[0][0], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau1_cov[0][1], 1.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau1_cov[0][2], 0.0, 1e-15) );

    // tangent2 to lower-x face is ez
    TEST_CHECK( gkyl_compare_double( cg->tau2_cov[0][0], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau2_cov[0][1], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau2_cov[0][2], 1.0, 1e-15) );

    // normal to lower-y face is ey
    TEST_CHECK( gkyl_compare_double( cg->norm_cov[1][0], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->norm_cov[1][1], 1.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->norm_cov[1][2], 0.0, 1e-15) );

    // tangent1 to lower-y face is ez
    TEST_CHECK( gkyl_compare_double( cg->tau1_cov[1][0], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau1_cov[1][1], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau1_cov[1][2], 1.0, 1e-15) );

    // tangent2 to lower-y face is ex
    TEST_CHECK( gkyl_compare_double( cg->tau2_cov[1][0], 1.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau2_cov[1][1], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau2_cov[1][2], 0.0, 1e-15) );

    // normal to lower-z face is ez
    TEST_CHECK( gkyl_compare_double( cg->norm_cov[2][0], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->norm_cov[2][1], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->norm_cov[2][2], 1.0, 1e-15) );

    // tangent1 to lower-z face is ex
    TEST_CHECK( gkyl_compare_double( cg->tau1_cov[2][0], 1.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau1_cov[2][1], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau1_cov[2][2], 0.0, 1e-15) );

    // tangent2 to lower-z face is ey
    TEST_CHECK( gkyl_compare_double( cg->tau2_cov[2][0], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau2_cov[2][1], 1.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau2_cov[2][2], 0.0, 1e-15) );
  }

  gkyl_wave_geom_release(wg);
}

// map (r,theta) -> (x,y)
void
mapc2p_cyl_local(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  double r = xc[0], th = xc[1], z = xc[2];
  xp[0] = r*cos(th); xp[1] = r*sin(th); xp[2] = z;
}

void
test_wv_geom_3d_2()
{
  int ndim = 3;
  double z_min = 0, z_max = 1;
  double r_inn = 0.25, r_out = 1.25;
  double phi_max = M_PI / 2.;
  double lower[] = {r_inn, 0, z_min}, upper[] = {r_out, phi_max, z_max};
  int cells[] = {1, 1, 1};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  double area = 0.5*(r_out*r_out-r_inn*r_inn);
  double edge_inn = sqrt(2) * r_inn;
  double area_c = (r_out - r_inn) * phi_max;

  // create range
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range range, ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  struct gkyl_wave_geom *wg = gkyl_wave_geom_new(&grid, &range, mapc2p_cyl_local, &ndim, false);
  // Also test cylindrical flag (which gives vectors in cylindrical unit basis)
  struct gkyl_wave_geom *wg_fl = gkyl_wave_geom_from_coord_flag(&grid, &range, WAVE_COORD_CYL, &ndim, false);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  
  while (gkyl_range_iter_next(&iter)) {
    const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wg, iter.idx);
    const struct gkyl_wave_cell_geom *cg_fl = gkyl_wave_geom_get(wg_fl, iter.idx);

    TEST_CHECK( gkyl_compare_double( cg->kappa, area / area_c, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->lenr[0], edge_inn / phi_max, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->lenr[1], 1, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->lenr[2], area / area_c, 1e-15) );

    // normal to left face has phi angle 45 deg
    TEST_CHECK( gkyl_compare_double( cg->norm_cov[0][0], 1/sqrt(2), 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->norm_cov[0][1], 1/sqrt(2), 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->norm_cov[0][2], 0.0, 1e-15) );

    // tangent1 to left face has phi angle 135 deg
    TEST_CHECK( gkyl_compare_double( cg->tau1_cov[0][0], -1/sqrt(2), 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau1_cov[0][1], 1/sqrt(2), 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau1_cov[0][2], 0.0, 1e-15) );

    // tangent2 to left face is ez
    TEST_CHECK( gkyl_compare_double( cg->tau2_cov[0][0], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau2_cov[0][1], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau2_cov[0][2], 1.0, 1e-15) );
    
    TEST_CHECK( gkyl_compare_double( cg_fl->kappa, area / area_c, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg_fl->lenr[0], edge_inn / phi_max, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg_fl->lenr[1], 1, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg_fl->lenr[2], area / area_c, 1e-15) );
    
    // normal to left face is in radial direction
    TEST_CHECK( gkyl_compare_double( cg_fl->norm_cov[0][0], 1.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg_fl->norm_cov[0][1], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg_fl->norm_cov[0][2], 0.0, 1e-15) );

    // tangent1 to left face is in azimuthal direction
    TEST_CHECK( gkyl_compare_double( cg_fl->tau1_cov[0][0], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg_fl->tau1_cov[0][1], 1.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg_fl->tau1_cov[0][2], 0.0, 1e-15) );

    // tangent2 to left face is vertical direction
    TEST_CHECK( gkyl_compare_double( cg_fl->tau2_cov[0][0], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg_fl->tau2_cov[0][1], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg_fl->tau2_cov[0][2], 1.0, 1e-15) );

    /* // normal to bottom face is ey */
    /* TEST_CHECK( gkyl_compare_double( cg->norm_cov[1][0], 0.0, 1e-15) ); */
    /* TEST_CHECK( gkyl_compare_double( cg->norm_cov[1][1], 1.0, 1e-15) ); */
    /* TEST_CHECK( gkyl_compare_double( cg->norm_cov[1][0], 0.0, 1e-15) ); */
    /*  */
    /* // tangent1 to bottom face is ex */
    /* TEST_CHECK( gkyl_compare_double( cg->tau1_cov[1][0], 1.0, 1e-15) ); */
    /* TEST_CHECK( gkyl_compare_double( cg->tau1_cov[1][1], 0.0, 1e-15) ); */
    /* TEST_CHECK( gkyl_compare_double( cg->tau1_cov[1][2], 0.0, 1e-15) ); */
    /*  */
    /* // tangent2 to bottom face is ez */
    /* TEST_CHECK( gkyl_compare_double( cg->tau2_cov[1][0], 0.0, 1e-15) ); */
    /* TEST_CHECK( gkyl_compare_double( cg->tau2_cov[1][1], 0.0, 1e-15) ); */
    /* TEST_CHECK( gkyl_compare_double( cg->tau2_cov[1][2], -1.0, 1e-15) ); */
  }

  gkyl_wave_geom_release(wg);
}

void
mapc2p_sph_local(double t, const double *xc, double *xp, void *ctx)
{
  const double r = xc[0], th = xc[1], ph = xc[2];
  xp[0] = r*sin(th)*cos(ph);
  xp[1] = r*sin(th)*sin(ph);
  xp[2] = r*cos(th);
}

void
test_wv_geom_3d_3()
{
  int ndim = 3;
  double r_inn = 0.25, r_out = 1.25;
  double theta_min = M_PI / 4.;
  double phi_min = 0.0;
  double theta_max = 3. * M_PI / 4.;
  double phi_max = M_PI / 2.;
  double lower[] = {r_inn, theta_min, phi_min}, upper[] = {r_out, theta_max, phi_max};
  int cells[] = {1, 1, 1};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  double vol = sqrt(2)/6.*(r_out*r_out*r_out-r_inn*r_inn*r_inn);
  double area_1 = sqrt(2) * r_inn*r_inn;
  double area_1_c = (phi_max - phi_min) * (theta_max - theta_min);
  double area_2 = sqrt(3) / 4. * (r_out - r_inn) * (r_out + r_inn);
  double area_2_c = (r_out - r_inn) * (phi_max - phi_min);
  double area_3 = (r_out - r_inn) * (r_inn + 0.5);
  double area_3_c = (r_out - r_inn) * (theta_max - theta_min);
  double vol_c = (r_out - r_inn) * (theta_max - theta_min) * (phi_max - phi_min);

  // create range
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range range, ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  struct gkyl_wave_geom *wg = gkyl_wave_geom_new(&grid, &range, mapc2p_sph_local, &ndim, false);
  // Also test spherical flag (which gives vectors in spherical unit basis)
  struct gkyl_wave_geom *wg_fl = gkyl_wave_geom_from_coord_flag(&grid, &range, WAVE_COORD_SPH, &ndim, false);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  
  while (gkyl_range_iter_next(&iter)) {
    const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wg, iter.idx);
    const struct gkyl_wave_cell_geom *cg_fl = gkyl_wave_geom_get(wg_fl, iter.idx);

    TEST_CHECK( gkyl_compare_double( cg->kappa, vol / vol_c, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->lenr[0], area_1 / area_1_c, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->lenr[1], area_2 / area_2_c, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->lenr[2], area_3 / area_3_c, 1e-15) );

    // normal to left face has phi angle 45 deg
    TEST_CHECK( gkyl_compare_double( cg->norm_cov[0][0], 1/sqrt(2), 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->norm_cov[0][1], 1/sqrt(2), 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->norm_cov[0][2], 0.0, 1e-15) );
    
    // tangent1 to left face is in -z direction (theta = 180 deg)
    TEST_CHECK( gkyl_compare_double( cg->tau1_cov[0][0], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau1_cov[0][1], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau1_cov[0][2], -1.0, 1e-15) );

    // tangent2 to left face has phi angle 135 deg
    TEST_CHECK( gkyl_compare_double( cg->tau2_cov[0][0], -1/sqrt(2), 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau2_cov[0][1], 1/sqrt(2), 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau2_cov[0][2], 0.0, 1e-15) );
    
    TEST_CHECK( gkyl_compare_double( cg_fl->kappa, vol / vol_c, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg_fl->lenr[0], area_1 / area_1_c, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg_fl->lenr[1], area_2 / area_2_c, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg_fl->lenr[2], area_3 / area_3_c, 1e-15) );
    
    // normal to left face is in radial direction
    TEST_CHECK( gkyl_compare_double( cg_fl->norm_cov[0][0], 1.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg_fl->norm_cov[0][1], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg_fl->norm_cov[0][2], 0.0, 1e-15) );

    // tangent1 to left face is in polar direction
    TEST_CHECK( gkyl_compare_double( cg_fl->tau1_cov[0][0], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg_fl->tau1_cov[0][1], 1.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg_fl->tau1_cov[0][2], 0.0, 1e-15) );

    // tangent2 to left face is azimuthal direction
    TEST_CHECK( gkyl_compare_double( cg_fl->tau2_cov[0][0], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg_fl->tau2_cov[0][1], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg_fl->tau2_cov[0][2], 1.0, 1e-15) );
  }

  gkyl_wave_geom_release(wg);
}

static void
get_cyl_cov_basis_local(double t, const double *xc, double *e1, double *e2, double *e3, void *ctx)
{
  const double r = xc[0], ph = xc[1];
  const double cosph = cos(ph), sinph = sin(ph);
  
  e1[0] = cosph; e1[1] = sinph; e1[2] = 0.0;
  e2[0] = -r*sinph; e2[1] = r*cosph; e2[2] = 0.0;
  e3[0] = 0.0; e3[1] = 0.0; e3[2] = 1.0;
}

static void
get_cyl_con_basis_local(double t, const double *xc, double *e1, double *e2, double *e3, void *ctx)
{
  const double r = xc[0], ph = xc[1];
  const double cosph = cos(ph), sinph = sin(ph);
  
  e1[0] = cosph; e1[1] = sinph; e1[2] = 0.0;
  e2[0] = -sinph/r; e2[1] = cosph/r; e2[2] = 0.0;
  e3[0] = 0.0; e3[1] = 0.0; e3[2] = 1.0;
}

void
test_wv_geom_3d_4()
{
  int ndim = 3;
  double z_min = 0, z_max = 1;
  double r_inn = 0.25, r_out = 1.25;
  double phi_max = M_PI / 2.;
  double lower[] = {r_inn, 0, z_min}, upper[] = {r_out, phi_max, z_max};
  int cells[] = {1, 1, 1};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  double area = 0.5*(r_out*r_out-r_inn*r_inn);
  double edge_inn = sqrt(2) * r_inn;
  double area_c = (r_out - r_inn) * phi_max;

  // create range
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range range, ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  struct gkyl_wave_geom_inp inp =
  {
    .mapc2p = mapc2p_cyl_local,
    .get_cov_basis = get_cyl_cov_basis_local,
    .get_con_basis = get_cyl_con_basis_local,
    .spacetime = 0
  };
  struct gkyl_wave_geom *wg = gkyl_wave_geom_from_wg_inp(&grid, &range, &inp, &ndim, false);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  
  while (gkyl_range_iter_next(&iter)) {
    const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wg, iter.idx);
    
    TEST_CHECK( gkyl_compare_double( cg->kappa, area / area_c, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->lenr[0], edge_inn / phi_max, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->lenr[1], 1, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->lenr[2], area / area_c, 1e-15) );
    
    // Check covariant components
    
    TEST_CHECK( gkyl_compare_double( cg->norm_cov[0][0], 1.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->norm_cov[0][1], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->norm_cov[0][2], 0.0, 1e-15) );

    TEST_CHECK( gkyl_compare_double( cg->tau1_cov[0][0], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau1_cov[0][1], 0.75, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau1_cov[0][2], 0.0, 1e-15) );

    TEST_CHECK( gkyl_compare_double( cg->tau2_cov[0][0], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau2_cov[0][1], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau2_cov[0][2], 1.0, 1e-15) );
    
    // Check contravariant components
    
    TEST_CHECK( gkyl_compare_double( cg->norm_con[0][0], 1.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->norm_con[0][1], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->norm_con[0][2], 0.0, 1e-15) );

    TEST_CHECK( gkyl_compare_double( cg->tau1_con[0][0], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau1_con[0][1], 4./3., 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau1_con[0][2], 0.0, 1e-15) );

    TEST_CHECK( gkyl_compare_double( cg->tau2_con[0][0], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau2_con[0][1], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau2_con[0][2], 1.0, 1e-15) );
  }

  gkyl_wave_geom_release(wg);
}

#ifdef GKYL_HAVE_CUDA

int cu_wave_geom_test(const struct gkyl_wave_geom *wg);

void
test_wv_geom_3d_cu()
{
  int ndim = 3;
  double z_min = 0, z_max = 1;
  double r_inn = 0.25, r_out = 1.25;
  double phi_max = M_PI / 2.;
  double lower[] = {r_inn, 0, z_min}, upper[] = {r_out, phi_max, z_max};
  int cells[] = {1, 1, 1};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  double area = 0.5*(r_out*r_out-r_inn*r_inn);
  double edge_inn = sqrt(2) * r_inn;
  double area_c = (r_out - r_inn) * phi_max;

  // create range
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range range, ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  struct gkyl_wave_geom *wg = gkyl_wave_geom_new(&grid, &range, mapc2p_cylind, &ndim, true);

  cu_wave_geom_test(wg->on_dev);

  gkyl_wave_geom_release(wg);
}
#endif

TEST_LIST = {
  { "wv_geom_1d_1", test_wv_geom_1d_1 },
  { "wv_geom_1d_2", test_wv_geom_1d_2 },
  { "wv_geom_2d_1", test_wv_geom_2d_1 },
  { "wv_geom_2d_2", test_wv_geom_2d_2 },
  { "wv_geom_2d_3", test_wv_geom_2d_3 },
  { "wv_geom_3d_1", test_wv_geom_3d_1 },
  { "wv_geom_3d_2", test_wv_geom_3d_2 },
  { "wv_geom_3d_3", test_wv_geom_3d_3 },
  { "wv_geom_3d_4", test_wv_geom_3d_4 },
#ifdef GKYL_HAVE_CUDA
  { "wv_geom_3d_cu", test_wv_geom_3d_cu },
#endif
  { NULL, NULL },
};
