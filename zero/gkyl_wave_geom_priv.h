#pragma once

// Private header, not for direct use in user code

#include <math.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_math.h>
#include <gkyl_evalf_def.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <gkyl_gr_spacetime.h>

#define GKYL_WG_DOT(v1, v2) v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]

/**
 * Free wave geometry object.
 *
 * @param ref Reference counter for wave geometry
 */
void gkyl_wave_geom_free(const struct gkyl_ref_count *ref);

static void
nomapc2p(double t, const double *xc, double *xp, void *ctx)
{
  for (int i=0; i<3; ++i) xp[i] = xc[i];
}

static void
get_standard_basis(double t, const double *xc, double *e1, double *e2, double *e3, void *ctx)
{
  e1[0] = 1.0; e1[1] = 0.0; e1[2] = 0.0;
  e2[0] = 0.0; e2[1] = 1.0; e2[2] = 0.0;
  e3[0] = 0.0; e3[1] = 0.0; e3[2] = 1.0;
}

static void
mapc2p_cyl(double t, const double *xc, double *xp, void *ctx)
{
  const double r = xc[0], ph = xc[1], z = xc[2];
  xp[0] = r*cos(ph);
  xp[1] = r*sin(ph);
  xp[2] = z;
}

static void
get_cyl_unit_basis(double t, const double *xc, double *e1, double *e2, double *e3, void *ctx)
{
  const double ph = xc[1];
  const double cosph = cos(ph), sinph = sin(ph);
  
  e1[0] = cosph; e1[1] = sinph; e1[2] = 0.0;
  e2[0] = -sinph; e2[1] = cosph; e2[2] = 0.0;
  e3[0] = 0.0; e3[1] = 0.0; e3[2] = 1.0;
}

static void
mapc2p_sph(double t, const double *xc, double *xp, void *ctx)
{
  const double r = xc[0], th = xc[1], ph = xc[2];
  xp[0] = r*sin(th)*cos(ph);
  xp[1] = r*sin(th)*sin(ph);
  xp[2] = r*cos(th);
}

static void
get_sph_unit_basis(double t, const double *xc, double *e1, double *e2, double *e3, void *ctx)
{
  const double th = xc[1], ph = xc[2];
  const double costh = cos(th), sinth = sin(th);
  const double cosph = cos(ph), sinph = sin(ph);
  
  e1[0] = sinth*cosph; e1[1] = sinth*sinph; e1[2] = costh;
  e2[0] = costh*cosph; e2[1] = costh*sinph; e2[2] = -sinth;
  e3[0] = -sinph; e3[1] = cosph; e3[2] = 0.0;
}

// Computes 1D geometry
static void
calc_geom_1d_from_nodes(const double *dx, const double *xlp, const double *xrp,
    evalf_t mapc2p, struct gkyl_vec3 e_cov[3], struct gkyl_vec3 e_con[3],
    void *ctx, struct gkyl_wave_cell_geom *geo)
{
  geo->kappa = fabs(xrp[0]-xlp[0])/dx[0];
  geo->lenr[0] = 1.0;

  get_standard_basis(0.0, 0, geo->norm_cov[0], geo->tau1_cov[0], geo->tau2_cov[0], ctx);
  get_standard_basis(0.0, 0, geo->norm_con[0], geo->tau1_con[0], geo->tau2_con[0], ctx);
  
  geo->norm_cov[0][0] = e_cov[0].x[0];
  geo->norm_con[0][0] = e_con[0].x[0];
}

static void
calc_geom_1d_from_spacetime(const double *dx, const double* xc,
  struct gkyl_gr_spacetime *spacetime, struct gkyl_wave_cell_geom *geo)
{
  get_standard_basis(0.0, 0, geo->norm_cov[0], geo->tau1_cov[0], geo->tau2_cov[0], 0);
  get_standard_basis(0.0, 0, geo->norm_con[0], geo->tau1_con[0], geo->tau2_con[0], 0);
  
  bool in_excision_region;
  spacetime->excision_region_func(spacetime, 0.0, xc[0], 0.0, 0.0, &in_excision_region);
  if (in_excision_region)
  {
    geo->kappa = 1.0;
    geo->lenr[0] = 1.0;
    return;
  }
  
  // Covariant metric components and determinant
  double gam_det;
  double **gam = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; ++i)
  {
    gam[i] = gkyl_malloc(sizeof(double[3]));
  }

  // Contravariant metric components
  double **gam_inv = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; ++i)
  {
    gam_inv[i] = gkyl_malloc(sizeof(double[3]));
  }
  
  spacetime->spatial_metric_det_func(spacetime, 0.0, xc[0], 0.0, 0.0, &gam_det);
  spacetime->spatial_metric_tensor_func(spacetime, 0.0, xc[0]-0.5*dx[0], 0.0, 0.0, &gam);
  spacetime->spatial_inv_metric_tensor_func(spacetime, 0.0, xc[0]-0.5*dx[0], 0.0, 0.0, &gam_inv);
    
  geo->kappa = sqrt(gam_det);
  geo->lenr[0] = 1.0;
    
  geo->norm_con[0][0] = sqrt(gam_inv[0][0]);
  geo->norm_cov[0][0] = geo->norm_con[0][0] * gam[0][0];
  
  for (int i = 0; i < 3; ++i)
  {
    gkyl_free(gam[i]);
    gkyl_free(gam_inv[i]);
  }
  gkyl_free(gam);
  gkyl_free(gam_inv);
}

static void
calc_geom_1d(const double *dx, const double *xc, struct gkyl_wave_geom_inp *inp,
    void *ctx, struct gkyl_wave_cell_geom *geo)
{
  if (inp->spacetime)
  {
    // Setup wave cell geom based on spacetime and ignore other inputs
    calc_geom_1d_from_spacetime(dx, xc, inp->spacetime, geo);
    return;
  }
  
  double xlc[GKYL_MAX_CDIM], xrc[GKYL_MAX_CDIM];
  double xlp[GKYL_MAX_CDIM], xrp[GKYL_MAX_CDIM];

  xlc[0] = xc[0]-0.5*dx[0]; // left node
  xrc[0] = xc[0]+0.5*dx[0]; // right node

  // compute coordinates of left/right nodes
  inp->mapc2p(0.0, xlc, xlp, ctx);
  inp->mapc2p(0.0, xrc, xrp, ctx);
  
  struct gkyl_vec3 e_cov[3] = {gkyl_vec3_zeros(), gkyl_vec3_zeros(), gkyl_vec3_zeros()};
  inp->get_cov_basis(0.0, xc, e_cov[0].x, e_cov[1].x, e_cov[2].x, ctx);
  struct gkyl_vec3 e_con[3] = {gkyl_vec3_zeros(), gkyl_vec3_zeros(), gkyl_vec3_zeros()};
  inp->get_con_basis(0.0, xc, e_con[0].x, e_con[1].x, e_con[2].x, ctx);

  calc_geom_1d_from_nodes(dx, xlp, xrp, inp->mapc2p, e_cov, e_con, ctx, geo);
}

// Computes 2D geometry
static void
calc_geom_2d_from_nodes(const double *dx, const struct gkyl_vec3 xll_p,
  const struct gkyl_vec3 xlr_p, const struct gkyl_vec3 xul_p,
  const struct gkyl_vec3 xur_p, evalf_t mapc2p,
  struct gkyl_vec3 e_cov[3], struct gkyl_vec3 e_con[3],
  void *ctx, struct gkyl_wave_cell_geom *geo)
{
  // ll: lower-left; lr: lower-right
  // ul: upper-left; ur: upper-right

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
  
  // convert norm, tau1; bases may have 3 nonzero vectors, so we zero out last component
  struct gkyl_vec3 norm_l_cov = gkyl_vec3_change_basis(norm_l, e_cov);
  struct gkyl_vec3 tau1_l_cov = gkyl_vec3_change_basis(tau1_l, e_cov);
  struct gkyl_vec3 norm_l_con = gkyl_vec3_change_basis(norm_l, e_con);
  struct gkyl_vec3 tau1_l_con = gkyl_vec3_change_basis(tau1_l, e_con);
  norm_l_cov.x[2] = tau1_l_cov.x[2] = norm_l_con.x[2] = tau1_l_cov.x[2] = 0.0;
  
  for (int d=0; d<3; ++d) {
    geo->norm_cov[0][d] = norm_l_cov.x[d];
    geo->tau1_cov[0][d] = tau1_l_cov.x[d];
    geo->tau2_cov[0][d] = tau2_l.x[d];
    
    geo->norm_con[0][d] = norm_l_con.x[d];
    geo->tau1_con[0][d] = tau1_l_con.x[d];
    geo->tau2_con[0][d] = tau2_l.x[d];
  }

  // normal-tangent to bottom face
  struct gkyl_vec3 tau1_b = gkyl_vec3_norm(gkyl_vec3_sub(xlr_p, xll_p));
  struct gkyl_vec3 tau2_b = gkyl_vec3_new(0.0, 0.0, -1.0); // -ez (ensures normal points into cell)
  struct gkyl_vec3 norm_b = gkyl_vec3_cross(tau1_b, tau2_b);
  
  // do same changes of basis as were done with left face
  struct gkyl_vec3 norm_b_cov = gkyl_vec3_change_basis(norm_b, e_cov);
  struct gkyl_vec3 tau1_b_cov = gkyl_vec3_change_basis(tau1_b, e_cov);
  struct gkyl_vec3 norm_b_con = gkyl_vec3_change_basis(norm_b, e_con);
  struct gkyl_vec3 tau1_b_con = gkyl_vec3_change_basis(tau1_b, e_con);
  norm_b_cov.x[2] = tau1_b_cov.x[2] = norm_b_con.x[2] = tau1_b_cov.x[2] = 0.0;

  for (int d=0; d<3; ++d) {
    geo->norm_cov[1][d] = norm_b_cov.x[d];
    geo->tau1_cov[1][d] = tau1_b_cov.x[d];
    geo->tau2_cov[1][d] = tau2_b.x[d];
    
    geo->norm_con[1][d] = norm_b_con.x[d];
    geo->tau1_con[1][d] = tau1_b_con.x[d];
    geo->tau2_con[1][d] = tau2_b.x[d];
  }
}

static void
calc_geom_2d_from_spacetime(const double *dx, const double* xc,
  struct gkyl_gr_spacetime *spacetime, struct gkyl_wave_cell_geom *geo)
{
  bool in_excision_region;
  spacetime->excision_region_func(spacetime, 0.0, xc[0], xc[1], 0.0, &in_excision_region);
  if (in_excision_region)
  {
    geo->kappa = 1.0;
    for (int e = 0; e < 2; ++e)
    {
      geo->lenr[e] = 1.0;
      get_standard_basis(0.0, 0, geo->norm_cov[e], geo->tau1_cov[e], geo->tau2_cov[e], 0);
      get_standard_basis(0.0, 0, geo->norm_con[e], geo->tau1_con[e], geo->tau2_con[e], 0);
    }
    return;
  }
  
  // Covariant metric components and determinant
  double gam_det;
  double **gam = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; ++i)
  {
    gam[i] = gkyl_malloc(sizeof(double[3]));
  }

  // Contravariant metric components
  double **gam_inv = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; ++i)
  {
    gam_inv[i] = gkyl_malloc(sizeof(double[3]));
  }
  
  spacetime->spatial_metric_det_func(spacetime, 0.0, xc[0], xc[1], xc[2], &gam_det);
  geo->kappa = sqrt(gam_det);
  
  // Permutations of indices to use when calculating triads (biads?)
  int ind[2][2] =
  {
    {0, 1},
    {1, 0}
  };
  
  // Coordinates to evaluate metric at (will be center of edge)
  double eval_x[2];
  
  for (int edge_idx = 0; edge_idx < 2; ++edge_idx)
  {
    int i = ind[edge_idx][0], j = ind[edge_idx][1];
    
    eval_x[0] = xc[0]; eval_x[1] = xc[1];
    eval_x[edge_idx] -= 0.5*dx[edge_idx];
    
    spacetime->spatial_metric_tensor_func(spacetime, 0.0, eval_x[0], eval_x[1], 0.0, &gam);
    spacetime->spatial_inv_metric_tensor_func(spacetime, 0.0, eval_x[0], eval_x[1], 0.0, &gam_inv);
    
    // Sqrt of determinant of restricted spatial metric tensor
    geo->lenr[edge_idx] = sqrt(gam[j][j]);
    
    // Contravariant components of orthonormal basis vectors
    double norm_con[3];
    double tau1_con[3];
    double tau2[3];
    
    norm_con[0] = gam_inv[i][0] / sqrt(gam_inv[i][i]);
    norm_con[1] = gam_inv[i][1] / sqrt(gam_inv[i][i]);
    norm_con[2] = 0.0;
    
    tau1_con[i] = tau1_con[2] = 0.0;
    tau1_con[j] = 1. / sqrt(gam[j][j]);
    
    tau2[i] = 0.0;
    tau2[j] = 0.0;
    tau2[2] = 1.0 - 2.0*i; // Matches other geometry function
    
    for (int d = 0; d < 3; ++d)
    {
      geo->norm_cov[edge_idx][d] = GKYL_WG_DOT(norm_con, gam[d]);
      geo->tau1_cov[edge_idx][d] = GKYL_WG_DOT(tau1_con, gam[d]);
      geo->tau2_cov[edge_idx][d] = tau2[d];
      
      geo->norm_con[edge_idx][d] = norm_con[d];
      geo->tau1_con[edge_idx][d] = tau1_con[d];
      geo->tau2_con[edge_idx][d] = tau2[d];
    }
  }
  
  for (int i = 0; i < 3; ++i)
  {
    gkyl_free(gam[i]);
    gkyl_free(gam_inv[i]);
  }
  gkyl_free(gam);
  gkyl_free(gam_inv);
}

static void
calc_geom_2d(const double *dx, const double *xc, struct gkyl_wave_geom_inp *inp,
    void *ctx, struct gkyl_wave_cell_geom *geo)
{
  if (inp->spacetime)
  {
    // Setup wave cell geom based on spacetime and ignore other inputs
    calc_geom_2d_from_spacetime(dx, xc, inp->spacetime, geo);
    return;
  }
  
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
  
  inp->mapc2p(0.0, xll_c.x, xll_p.x, ctx);
  inp->mapc2p(0.0, xlr_c.x, xlr_p.x, ctx);
  inp->mapc2p(0.0, xul_c.x, xul_p.x, ctx);
  inp->mapc2p(0.0, xur_c.x, xur_p.x, ctx);
  
  struct gkyl_vec3 e_cov[3] = {gkyl_vec3_zeros(), gkyl_vec3_zeros(), gkyl_vec3_zeros()};
  inp->get_cov_basis(0.0, xc, e_cov[0].x, e_cov[1].x, e_cov[2].x, ctx);
  struct gkyl_vec3 e_con[3] = {gkyl_vec3_zeros(), gkyl_vec3_zeros(), gkyl_vec3_zeros()};
  inp->get_con_basis(0.0, xc, e_con[0].x, e_con[1].x, e_con[2].x, ctx);

  // need to set the final coordinate to 0.0
  xll_p.x[2] = xlr_p.x[2] = xul_p.x[2] = xur_p.x[2] = 0.0;

  calc_geom_2d_from_nodes(dx, xll_p, xlr_p, xul_p, xur_p, inp->mapc2p, e_cov, e_con, ctx, geo);
}

// Computes 3D geometry
static double
vol_tetra(const struct gkyl_vec3 p1,
          const struct gkyl_vec3 p2,
          const struct gkyl_vec3 p3,
          const struct gkyl_vec3 p4)
{
  struct gkyl_vec3 a = gkyl_vec3_sub(p1, p4);
  struct gkyl_vec3 b = gkyl_vec3_sub(p2, p4);
  struct gkyl_vec3 c = gkyl_vec3_sub(p3, p4);
  double vol = gkyl_vec3_dot(a, gkyl_vec3_cross(b, c)) / 6.;
  return vol > 0? vol : -vol;
}

/* The order of vertices is assumed to be (l: lower, u: upper):
[0], lll: (i,   j,   k  ) or (i-0.5, j-0.5, k-0.5)
[1], ull: (i+1, j,   k  ) or (i+0.5, j-0.5, k-0.5)
[2], uul: (i+1, j+1, k  ) or (i+0.5, j+0.5, k-0.5)
[3], lul: (i,   j+1, k  ) or (i-0.5, j+0.5, k-0.5)
[4], llu: (i,   j,   k+1) or (i-0.5, j-0.5, k+0.5)
[5], ulu: (i+1, j,   k+1) or (i+0.5, j-0.5, k+0.5)
[6], uuu: (i+1, j+1, k+1) or (i+0.5, j+0.5, k+0.5)
[7], luu: (i,   j+1, k+1) or (i-0.5, j+0.5, k+0.5)

       y
3----------2
|\     ^   |\
| \    |   | \
|  \   |   |  \
|   7------+---6
|   |  *-- |-- | -> x
0---+---\--1   |
 \  |    \  \  |
  \ |     \  \ |
   \|      z  \|
    4----------5
*/
static double
vol_hexa(const struct gkyl_vec3 *verts)
{
  // split the hexahedron into five tetrahedrons and add up their volumes
  // FIXME does this handle bad hexahedrons?
  double vol = 0;
  vol += vol_tetra(verts[1], verts[5], verts[4], verts[6]);
  vol += vol_tetra(verts[7], verts[4], verts[6], verts[3]);
  vol += vol_tetra(verts[4], verts[0], verts[1], verts[3]);
  vol += vol_tetra(verts[6], verts[1], verts[3], verts[4]);
  vol += vol_tetra(verts[1], verts[2], verts[6], verts[3]);

  return vol;
}

static inline double
triangle_area(const struct gkyl_vec3 p1,
              const struct gkyl_vec3 p2,
              const struct gkyl_vec3 p3)
{
  return 0.5*gkyl_vec3_len(
      gkyl_vec3_cross(gkyl_vec3_sub(p1, p2), gkyl_vec3_sub(p2, p3)));
}

// Points are in anti-clockwise order, i.e., p1-p3 and p2-p4 are diagonals.
static double
planar_quad_area_norm(const struct gkyl_vec3 p1,
                        const struct gkyl_vec3 p2,
                        const struct gkyl_vec3 p3,
                        const struct gkyl_vec3 p4,
                        struct gkyl_vec3 *norm)
{
  struct gkyl_vec3 v13 = gkyl_vec3_sub(p3, p1);
  struct gkyl_vec3 v24 = gkyl_vec3_sub(p4, p2);

  double area = 0.5 * gkyl_vec3_len(gkyl_vec3_cross(v13, v24));

  if (norm) {
    struct gkyl_vec3 norm_ = gkyl_vec3_norm(gkyl_vec3_cross(v13, v24));

    for (int d=0; d<3; ++d)
    {
      norm->x[d] = norm_.x[d];
    }
  }

  return area;
}

// ca * a + cb * b
static inline struct gkyl_vec3
gkyl_vec3_add_coeff(
  const double ca, struct gkyl_vec3 a, const double cb, struct gkyl_vec3 b)
{
  return (struct gkyl_vec3) { .x = {
    ca*a.x[0] + cb*b.x[0],
    ca*a.x[1] + cb*b.x[1],
    ca*a.x[2] + cb*b.x[2] } };
}

static double
quad_area_norm_tang(const struct gkyl_vec3 p1,
                    const struct gkyl_vec3 p2,
                    const struct gkyl_vec3 p3,
                    const struct gkyl_vec3 p4,
                    struct gkyl_vec3 *norm,
                    struct gkyl_vec3 *tau1,
                    struct gkyl_vec3 *tau2)
{
  // find Varignon parallelogram and use its normal as the quad's normal
  struct gkyl_vec3 pp1 = gkyl_vec3_add_coeff(0.5, p1, 0.5, p2);
  struct gkyl_vec3 pp2 = gkyl_vec3_add_coeff(0.5, p2, 0.5, p3);
  struct gkyl_vec3 pp3 = gkyl_vec3_add_coeff(0.5, p3, 0.5, p4);
  struct gkyl_vec3 pp4 = gkyl_vec3_add_coeff(0.5, p4, 0.5, p1);
  double area = planar_quad_area_norm(pp1, pp2, pp3, pp4, norm);

  // tau1 is v(p1->p2) projected on the normal plane, (tau1, tau2, norm)
  // completes a high-hand system
  struct gkyl_vec3 p12 = gkyl_vec3_sub(p2, p1);
  double d = gkyl_vec3_dot(p12, *norm);
  struct gkyl_vec3 t1 = gkyl_vec3_norm(gkyl_vec3_add_coeff(1, p12, -d, *norm));
  struct gkyl_vec3 t2 = gkyl_vec3_cross(*norm, t1);
  for (int d=0; d<3; ++d)
  {
    tau1->x[d] = t1.x[d];
    tau2->x[d] = t2.x[d];
  }

#if 1
  // project corner points onto the normal plane and use the projected planar
  // quad's area as the original quad's area
  d = gkyl_vec3_dot(*norm, p1);
  struct gkyl_vec3 proj1 = gkyl_vec3_add_coeff(1, p1, -d, *norm);
  d = gkyl_vec3_dot(*norm, p2);
  struct gkyl_vec3 proj2 = gkyl_vec3_add_coeff(1, p2, -d, *norm);
  d = gkyl_vec3_dot(*norm, p3);
  struct gkyl_vec3 proj3 = gkyl_vec3_add_coeff(1, p3, -d, *norm);
  d = gkyl_vec3_dot(*norm, p4);
  struct gkyl_vec3 proj4 = gkyl_vec3_add_coeff(1, p4, -d, *norm);
  
  area = planar_quad_area_norm(proj1, proj2, proj3, proj4, norm);
#else
  // accumulate triangle areas
  area += triangle_area(p1, pp1, pp4);
  area += triangle_area(p2, pp2, pp1);
  area += triangle_area(p3, pp3, pp2);
  area += triangle_area(p4, pp4, pp3);
#endif

  return area;
}

static void
calc_geom_3d_from_nodes(const double *dx, struct gkyl_vec3 verts[8],
    evalf_t mapc2p, struct gkyl_vec3 e_cov[3], struct gkyl_vec3 e_con[3],
    void *ctx, struct gkyl_wave_cell_geom *geo)
{
  // compute cell volume and kappa
  double vol = vol_hexa(verts);
  double cell_vol_c = dx[0] * dx[1] * dx[2];
  geo->kappa = vol / cell_vol_c;

  // for each of the thee lower quad faces 'owned' by the present cell, compute
  // area, norm, tan1, tan2, and lenr
  int pt_idx[3][4] = { // indices of verices of each face in verts; see vol_hexa
    // lower-x face, lll, lul, luu, llu; v(p1->p2)=ey
    {0, 3, 7, 4},
    // lower-y face, lll, llu, ulu, ull; v(p1->p2)=ez
    {0, 4, 5, 1},
    // lower-z face, lll, ull, uul, lul; v(p1->p2)=ex
    {0, 1, 2, 3}
  };
  for (int face_idx = 0; face_idx < 3; ++face_idx)
  {
    int ip1 = pt_idx[face_idx][0];
    int ip2 = pt_idx[face_idx][1];
    int ip3 = pt_idx[face_idx][2];
    int ip4 = pt_idx[face_idx][3];

    struct gkyl_vec3 norm, tau1, tau2;
    double area= quad_area_norm_tang(
        verts[ip1], verts[ip2], verts[ip3], verts[ip4], &norm, &tau1, &tau2);
    double cell_area_c = cell_vol_c / dx[face_idx];
    geo->lenr[face_idx] = area / cell_area_c;
    
    // Convert norm, tau1, tau2 to computational covariant and contravariant bases
    struct gkyl_vec3 norm_cov = gkyl_vec3_change_basis(norm, e_cov);
    struct gkyl_vec3 tau1_cov = gkyl_vec3_change_basis(tau1, e_cov);
    struct gkyl_vec3 tau2_cov = gkyl_vec3_change_basis(tau2, e_cov);
    
    struct gkyl_vec3 norm_con = gkyl_vec3_change_basis(norm, e_con);
    struct gkyl_vec3 tau1_con = gkyl_vec3_change_basis(tau1, e_con);
    struct gkyl_vec3 tau2_con = gkyl_vec3_change_basis(tau2, e_con);
    
    for (int d=0; d<3; ++d)
    {
      geo->norm_cov[face_idx][d] = norm_cov.x[d];
      geo->tau1_cov[face_idx][d] = tau1_cov.x[d];
      geo->tau2_cov[face_idx][d] = tau2_cov.x[d];
      
      geo->norm_con[face_idx][d] = norm_con.x[d];
      geo->tau1_con[face_idx][d] = tau1_con.x[d];
      geo->tau2_con[face_idx][d] = tau2_con.x[d];
    }
  }
}

static void
calc_geom_3d_from_spacetime(const double *dx, const double* xc,
  struct gkyl_gr_spacetime *spacetime, struct gkyl_wave_cell_geom *geo)
{
  bool in_excision_region;
  spacetime->excision_region_func(spacetime, 0.0, xc[0], xc[1], xc[2], &in_excision_region);
  if (in_excision_region)
  {
    geo->kappa = 1.0;
    for (int e = 0; e < 3; ++e)
    {
      geo->lenr[e] = 1.0;
      for (int d = 0; d < 3; ++d)
      {
        get_standard_basis(0.0, 0, geo->norm_cov[e], geo->tau1_cov[e], geo->tau2_cov[e], 0);
        get_standard_basis(0.0, 0, geo->norm_con[e], geo->tau1_con[e], geo->tau2_con[e], 0);
      }
    }
    return;
  }
  
  // Covariant metric components and determinant
  double gam_det;
  double **gam = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; ++i)
  {
    gam[i] = gkyl_malloc(sizeof(double[3]));
  }

  // Contravariant metric components
  double **gam_inv = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; ++i)
  {
    gam_inv[i] = gkyl_malloc(sizeof(double[3]));
  }
  
  spacetime->spatial_metric_det_func(spacetime, 0.0, xc[0], xc[1], xc[2], &gam_det);
  geo->kappa = sqrt(gam_det);
  
  // Permutations of indices to use when calculating triads
  int ind[3][3] =
  {
    {0, 1, 2},
    {1, 0, 2},
    {2, 0, 1}
  };
  
  // Coordinates to evaluate metric at (will be center of face)
  double eval_x[3];
  
  for (int face_idx = 0; face_idx < 3; ++face_idx)
  {
    int i = ind[face_idx][0], j = ind[face_idx][1], k = ind[face_idx][2];
    
    eval_x[0] = xc[0]; eval_x[1] = xc[1]; eval_x[2] = xc[2];
    eval_x[face_idx] -= 0.5*dx[face_idx];
    
    spacetime->spatial_metric_tensor_func(spacetime, 0.0, eval_x[0], eval_x[1], eval_x[2], &gam);
    spacetime->spatial_inv_metric_tensor_func(spacetime, 0.0, eval_x[0], eval_x[1], eval_x[2], &gam_inv);
    
    // Sqrt of determinant of restricted spatial metric tensor
    geo->lenr[face_idx] = sqrt(gam[j][j]*gam[k][k] - gam[j][k]*gam[k][j]);
    
    // Contravariant components of orthonormal basis vectors
    double norm_con[3];
    double tau1_con[3];
    double tau2_con[3];
    
    for (int d = 0; d < 3; ++d)
    {
      norm_con[d] = gam_inv[i][d] / sqrt(gam_inv[i][i]);
    }
    
    tau1_con[i] = tau1_con[k] = 0.0;
    tau1_con[j] = 1. / sqrt(gam[j][j]);
    
    tau2_con[i] = 0.0;
    tau2_con[j] = -gam[j][k] / (sqrt(gam[j][j]) * geo->lenr[face_idx]);
    tau2_con[k] = sqrt(gam[j][j]) / geo->lenr[face_idx];
    
    for (int d = 0; d < 3; ++d)
    {
      geo->norm_cov[face_idx][d] = GKYL_WG_DOT(norm_con, gam[d]);
      geo->tau1_cov[face_idx][d] = GKYL_WG_DOT(tau1_con, gam[d]);
      geo->tau2_cov[face_idx][d] = GKYL_WG_DOT(tau2_con, gam[d]);
      
      geo->norm_con[face_idx][d] = norm_con[d];
      geo->tau1_con[face_idx][d] = tau1_con[d];
      geo->tau2_con[face_idx][d] = tau2_con[d];
    }
  }
  
  for (int i = 0; i < 3; ++i)
  {
    gkyl_free(gam[i]);
    gkyl_free(gam_inv[i]);
  }
  gkyl_free(gam);
  gkyl_free(gam_inv);
}

static void
calc_geom_3d(const double *dx, const double *xc, struct gkyl_wave_geom_inp *inp,
    void *ctx, struct gkyl_wave_cell_geom *geo)
{
  if (inp->spacetime)
  {
    // Setup wave cell geom based on spacetime and ignore other inputs
    calc_geom_3d_from_spacetime(dx, xc, inp->spacetime, geo);
    return;
  }
  
  // get all vertices of the hexahedron; see vol_hexa for their order
  struct gkyl_vec3 verts_c[8] = {
    gkyl_vec3_new(xc[0]-0.5*dx[0], xc[1]-0.5*dx[1], xc[2]-0.5*dx[2]),
    gkyl_vec3_new(xc[0]+0.5*dx[0], xc[1]-0.5*dx[1], xc[2]-0.5*dx[2]),
    gkyl_vec3_new(xc[0]+0.5*dx[0], xc[1]+0.5*dx[1], xc[2]-0.5*dx[2]),
    gkyl_vec3_new(xc[0]-0.5*dx[0], xc[1]+0.5*dx[1], xc[2]-0.5*dx[2]),
    gkyl_vec3_new(xc[0]-0.5*dx[0], xc[1]-0.5*dx[1], xc[2]+0.5*dx[2]),
    gkyl_vec3_new(xc[0]+0.5*dx[0], xc[1]-0.5*dx[1], xc[2]+0.5*dx[2]),
    gkyl_vec3_new(xc[0]+0.5*dx[0], xc[1]+0.5*dx[1], xc[2]+0.5*dx[2]),
    gkyl_vec3_new(xc[0]-0.5*dx[0], xc[1]+0.5*dx[1], xc[2]+0.5*dx[2])
  };

  struct gkyl_vec3 verts[8]; // physical coordinate nodes
  for (int i=0; i<8; ++i)
    inp->mapc2p(0.0, verts_c[i].x, verts[i].x, ctx);
  
  struct gkyl_vec3 e_cov[3] = {gkyl_vec3_zeros(), gkyl_vec3_zeros(), gkyl_vec3_zeros()};
  inp->get_cov_basis(0.0, xc, e_cov[0].x, e_cov[1].x, e_cov[2].x, ctx);
  struct gkyl_vec3 e_con[3] = {gkyl_vec3_zeros(), gkyl_vec3_zeros(), gkyl_vec3_zeros()};
  inp->get_con_basis(0.0, xc, e_con[0].x, e_con[1].x, e_con[2].x, ctx);
  
  calc_geom_3d_from_nodes(dx, verts, inp->mapc2p, e_cov, e_con, ctx, geo);
}

