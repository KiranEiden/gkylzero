#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_math.h>
#include <gkyl_util.h>
#include <gkyl_wave_geom.h>

#include <math.h>

static void
nomapc2p(double t, const double *xc, double *xp, void *ctx)
{
  for (int i=0; i<3; ++i) xp[i] = xc[i];
}

// Computes 1D geometry
static void
calc_geom_1d(const double *dx, const double *xc, evalf_t mapc2p, void *ctx, struct gkyl_wave_cell_geom *geo)
{
  double xlc[GKYL_MAX_CDIM], xrc[GKYL_MAX_CDIM];
  double xlp[GKYL_MAX_CDIM], xrp[GKYL_MAX_CDIM];

  xlc[0] = xc[0]-0.5*dx[0]; // left node
  xrc[0] = xc[0]+0.5*dx[0]; // right node

  // compute coordinates of left/right nodes
  mapc2p(0.0, xlc, xlp, ctx);
  mapc2p(0.0, xrc, xrp, ctx);

  geo->kappa = fabs(xrp[0]-xlp[0])/dx[0];
  geo->lenr[0] = 1.0;

  geo->norm[0][0] = 1.0; geo->norm[0][1] = 0.0; geo->norm[0][2] = 0.0;
  geo->tau1[0][0] = 0.0; geo->tau1[0][1] = 1.0; geo->tau1[0][2] = 0.0;
  geo->tau2[0][0] = 0.0; geo->tau2[0][1] = 0.0; geo->tau2[0][2] = 1.0;
}

static void
calc_geom_2d(const double *dx, const double *xc, evalf_t mapc2p, void *ctx, struct gkyl_wave_cell_geom *geo)
{
  // ll: lower-left; lr: lower-right
  // ul: upper-left; ur: upper-right
  
  struct gkyl_vec3 xll_p = gkyl_vec3_zeros();
  struct gkyl_vec3 xlr_p = gkyl_vec3_zeros();
  struct gkyl_vec3 xul_p = gkyl_vec3_zeros();
  struct gkyl_vec3 xur_p = gkyl_vec3_zeros();

  struct gkyl_vec3 xll_c = gkyl_vec3_new(xc[0] - 0.5*dx[0],  xc[1] - 0.5*dx[1], 0.0);
  struct gkyl_vec3 xlr_c = gkyl_vec3_new(xc[0] + 0.5*dx[0],  xc[1] - 0.5*dx[1], 0.0);

  struct gkyl_vec3 xul_c = gkyl_vec3_new(xc[0] - 0.5*dx[0],  xc[1] + 0.5*dx[1], 0.0);
  struct gkyl_vec3 xur_c = gkyl_vec3_new(xc[0] + 0.5*dx[0],  xc[1] + 0.5*dx[1], 0.0);
  
  mapc2p(0.0, xll_c.x, xll_p.x, ctx);
  mapc2p(0.0, xlr_c.x, xlr_p.x, ctx);
  mapc2p(0.0, xul_c.x, xul_p.x, ctx);
  mapc2p(0.0, xur_c.x, xur_p.x, ctx);

  // need to set the final coordinate to 0.0
  xll_p.x[2] = xlr_p.x[2] = xul_p.x[2] = xur_p.x[2] = 0.0;

  // volume factor
  double area = 0.5*gkyl_vec3_len( gkyl_vec3_cross(gkyl_vec3_sub(xlr_p,xll_p), gkyl_vec3_sub(xul_p,xll_p)) )
    + 0.5*gkyl_vec3_len( gkyl_vec3_cross(gkyl_vec3_sub(xlr_p,xur_p), gkyl_vec3_sub(xul_p,xur_p)) );

  geo->kappa = area/(dx[0]*dx[1]);

  // face-area ratios for faces (a face is an edge in 2D)
  geo->lenr[0] = gkyl_vec3_len(gkyl_vec3_sub(xul_p, xll_p))/dx[1];
  geo->lenr[1] = gkyl_vec3_len(gkyl_vec3_sub(xlr_p, xll_p))/dx[0];

  // normal-tangent to left face
  struct gkyl_vec3 tau1_l = gkyl_vec3_norm(gkyl_vec3_sub(xul_p, xll_p));
  struct gkyl_vec3 tau2_l = gkyl_vec3_new(0.0, 0.0, 1.0); // ez
  struct gkyl_vec3 norm_l = gkyl_vec3_cross(tau1_l, tau2_l);
  
  for (int d=0; d<3; ++d) {
    geo->norm[0][d] = norm_l.x[d];
    geo->tau1[0][d] = tau1_l.x[d];
    geo->tau2[0][d] = tau2_l.x[d];
  }

  // normal-tangent to bottom face
  struct gkyl_vec3 tau1_b = gkyl_vec3_norm(gkyl_vec3_sub(xlr_p, xll_p));
  struct gkyl_vec3 tau2_b = gkyl_vec3_new(0.0, 0.0, -1.0); // -ez (ensures normal points into cell)
  struct gkyl_vec3 norm_b = gkyl_vec3_cross(tau1_b, tau2_b);

  for (int d=0; d<3; ++d) {
    geo->norm[1][d] = norm_b.x[d];
    geo->tau1[1][d] = tau1_b.x[d];
    geo->tau2[1][d] = tau2_b.x[d];
  }
}

static void
calc_geom_3d(const double *dx, const double *xc, evalf_t mapc2p, void *ctx, struct gkyl_wave_cell_geom *geo)
{
}

struct gkyl_wave_geom*
gkyl_wave_geom_new(const struct gkyl_rect_grid *grid, struct gkyl_range *range,
  evalf_t mapc2p, void *ctx)
{
  struct gkyl_wave_geom *wg = gkyl_malloc(sizeof(struct gkyl_wave_geom));

  wg->range = *range;
  wg->geom = gkyl_array_new(GKYL_USER, sizeof(struct gkyl_wave_cell_geom), range->volume);

  double xc[GKYL_MAX_CDIM];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {

    gkyl_rect_grid_cell_center(grid, iter.idx, xc);

    struct gkyl_wave_cell_geom *geo = gkyl_array_fetch(wg->geom, gkyl_range_idx(range, iter.idx));
    
    // compute geometry based on grid dimensions
    switch (grid->ndim) {
      case 1:
        calc_geom_1d(grid->dx, xc, mapc2p ? mapc2p : nomapc2p, ctx, geo);
        break;

      case 2:
        calc_geom_2d(grid->dx, xc, mapc2p ? mapc2p : nomapc2p, ctx, geo);
        break;

      case 3:
        calc_geom_3d(grid->dx, xc, mapc2p ? mapc2p : nomapc2p, ctx, geo);
        break;
    };   
  }

  return wg;
}


const struct gkyl_wave_cell_geom*
gkyl_wave_geom_get(const struct gkyl_wave_geom *wg, const int *idx)
{
  return gkyl_array_cfetch(wg->geom, gkyl_range_idx(&wg->range, idx));
}

void
gkyl_wave_geom_release(struct gkyl_wave_geom *wg)
{
  gkyl_array_release(wg->geom);
  gkyl_free(wg);
}