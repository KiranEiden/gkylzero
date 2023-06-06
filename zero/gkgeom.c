#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_gkgeom.h>
#include <gkyl_rect_grid.h>
#include <gkyl_range.h>

#include <math.h>
#include <string.h>

struct gkyl_gkgeom {
  struct gkyl_rect_grid rzgrid; // RZ grid on which psi(R,Z) is defined
  const struct gkyl_array *psiRZ; // psi(R,Z) DG representation
  struct gkyl_range rzlocal; // local range over which psiRZ is defined
  int num_rzbasis; // number of basis functions in RZ

  struct { int max_iter; double eps; } root_param;
  struct { double eps; } quad_param;

  // pointer to root finder (depends on polyorder)
  struct RdRdZ_sol (*calc_roots)(const double *psi, double psi0, double Z,
    double xc[2], double dx[2]);
};

// some helper functions
static inline double
choose_closest(int ref, double R[2], double out[2])
{
  return fabs(R[0]-ref) < fabs(R[1]-ref) ? out[0] : out[1];
}

static inline double SQ(double x) { return x*x; }

static inline int
get_idx(int dir, double x, const struct gkyl_rect_grid *grid, const struct gkyl_range *range)
{
  double xlower = grid->lower[dir], dx = grid->dx[dir];
  int idx = range->lower[dir] + (int) floor((x-xlower)/dx);
  return idx <= range->upper[dir] ? idx : range->upper[dir];
}

// struct for solutions to roots
struct RdRdZ_sol {
  int nsol;
  double R[2], dRdZ[2];
};

// Compute roots R(psi,Z) and dR/dZ(psi,Z) in a p=1 DG cell
static inline struct RdRdZ_sol
calc_RdR_p1(const double *psi, double psi0, double Z, double xc[2], double dx[2])
{
  struct RdRdZ_sol sol = { .nsol = 0 };

  double y = (Z-xc[1])/(dx[1]*0.5);
  
  double rnorm = (-(1.732050807568877*psi[2]*y)/(3.0*psi[3]*y+1.732050807568877*psi[1]))+(2.0*psi0)/(3.0*psi[3]*y+1.732050807568877*psi[1])-(1.0*psi[0])/(3.0*psi[3]*y+1.732050807568877*psi[1]) ;

  if ((-1<=rnorm) && (rnorm < 1)) {
    double drdznorm = -(3.0*(2.0*psi[3]*psi0-1.0*psi[0]*psi[3]+psi[1]*psi[2]))/SQ(3.0*psi[3]*y+1.732050807568877*psi[1]) ;
    
    sol.nsol = 1;
    sol.R[0] = rnorm*dx[0]*0.5 + xc[0];
    sol.dRdZ[0] = drdznorm*dx[0]/dx[1];
  }
  return sol;
}

// Compute roots R(psi,Z) and dR/dZ(psi,Z) in a p=2 DG cell
static inline struct RdRdZ_sol
calc_RdR_p2(const double *psi, double psi0, double Z, double xc[2], double dx[2])
{
  struct RdRdZ_sol sol = { .nsol = 0 };
  double y = (Z-xc[1])/(dx[1]*0.5);

  double aq = 2.904737509655563*psi[6]*y+1.677050983124842*psi[4]; 
  double bq = 2.904737509655563*psi[7]*SQ(y)+1.5*psi[3]*y-0.9682458365518543*psi[7]+0.8660254037844386*psi[1]; 
  double cq = 1.677050983124842*psi[5]*SQ(y)-0.9682458365518543*psi[6]*y+0.8660254037844386*psi[2]*y-1.0*psi0-0.5590169943749475*psi[5]-0.5590169943749475*psi[4]+0.5*psi[0]; 
  double delta2 = bq*bq - 4*aq*cq;

  if (delta2 > 0) {
    double r1, r2;
    double delta = sqrt(delta2);
    // compute both roots
    if (bq>=0) {
      r1 = (-bq-delta)/(2*aq);
      r2 = 2*cq/(-bq-delta);
    }
    else {
      r1 = 2*cq/(-bq+delta);
      r2 = (-bq+delta)/(2*aq);
    }

    int sidx = 0;
    if ((-1<=r1) && (r1 < 1)) {
      sol.nsol += 1;
      sol.R[sidx] = r1*dx[0]*0.5 + xc[0];

      double x = r1;
      double C = 5.809475019311126*psi[7]*x*y+3.354101966249685*psi[5]*y+2.904737509655563*psi[6]*SQ(x)+1.5*psi[3]*x-0.9682458365518543*psi[6]+0.8660254037844386*psi[2]; 
      double A = 2.904737509655563*psi[7]*SQ(y)+5.809475019311126*psi[6]*x*y+1.5*psi[3]*y+3.354101966249685*psi[4]*x-0.9682458365518543*psi[7]+0.8660254037844386*psi[1];
      sol.dRdZ[sidx] = -C/A*dx[0]/dx[1];
      
      sidx += 1;
    }
    if ((-1<=r2) && (r2 < 1)) {
      sol.nsol += 1;
      sol.R[sidx] = r2*dx[0]*0.5 + xc[0];

      double x = r2;
      double C = 5.809475019311126*psi[7]*x*y+3.354101966249685*psi[5]*y+2.904737509655563*psi[6]*SQ(x)+1.5*psi[3]*x-0.9682458365518543*psi[6]+0.8660254037844386*psi[2]; 
      double A = 2.904737509655563*psi[7]*SQ(y)+5.809475019311126*psi[6]*x*y+1.5*psi[3]*y+3.354101966249685*psi[4]*x-0.9682458365518543*psi[7]+0.8660254037844386*psi[1];
      sol.dRdZ[sidx] = -C/A*dx[0]/dx[1];
      
      sidx += 1;
    }
  }
  return sol;
}

// Compute R(psi,Z) given a psi and Z. Can return multiple solutions
// or no solutions. The number of roots found is returned and are
// copied in the array R and dR. The calling function must ensure that
// these arrays are big enough to hold all roots required
int
R_psiZ(const gkyl_gkgeom *geo, double psi, double Z, int nmaxroots,
  double *R, double *dR)
{
  int zcell = get_idx(1, Z, &geo->rzgrid, &geo->rzlocal);

  int sidx = 0;
  int idx[2] = { 0, zcell };
  double dx[2] = { geo->rzgrid.dx[0], geo->rzgrid.dx[1] };
  
  struct gkyl_range rangeR;
  gkyl_range_deflate(&rangeR, &geo->rzlocal, (int[]) { 0, 1 }, (int[]) { 0, zcell });

  struct gkyl_range_iter riter;
  gkyl_range_iter_init(&riter, &rangeR);
  
  // loop over all R cells to find psi crossing
  while (gkyl_range_iter_next(&riter) && sidx<=nmaxroots) {
    long loc = gkyl_range_idx(&rangeR, riter.idx);
    const double *psih = gkyl_array_cfetch(geo->psiRZ, loc);

    double xc[2];
    idx[0] = riter.idx[0];
    gkyl_rect_grid_cell_center(&geo->rzgrid, idx, xc);

    struct RdRdZ_sol sol = geo->calc_roots(psih, psi, Z, xc, dx);
    
    if (sol.nsol > 0)
      for (int s=0; s<sol.nsol; ++s) {
        R[sidx] = sol.R[s];
        dR[sidx] = sol.dRdZ[s];
      }
    sidx += sol.nsol;
  }
  return sidx;
}

gkyl_gkgeom*
gkyl_gkgeom_new(const struct gkyl_gkgeom_inp *inp)
{
  struct gkyl_gkgeom *geo = gkyl_malloc(sizeof(*geo));

  geo->rzgrid = *inp->rzgrid;
  geo->psiRZ = gkyl_array_acquire(inp->psiRZ);
  geo->num_rzbasis = inp->rzbasis->num_basis;
  memcpy(&geo->rzlocal, inp->rzlocal, sizeof(struct gkyl_range));

  geo->root_param.eps = inp->root_param.eps > 0 ? inp->root_param.eps : 1e-10;
  geo->root_param.max_iter = inp->root_param.max_iter > 0 ? inp->root_param.max_iter : 100;

  geo->quad_param.eps = inp->quad_param.eps > 0 ? inp->quad_param.eps : 1e-10;

  if (inp->rzbasis->poly_order == 1)
    geo->calc_roots = calc_RdR_p1;
  else if (inp->rzbasis->poly_order == 2)
    geo->calc_roots = calc_RdR_p2;
  
  return geo;
}

double
gkyl_gkgeom_integrate_psi_contour(const gkyl_gkgeom *geo, double psi,
  double zmin, double zmax, double rclose)
{
  return 0;
}

void
gkyl_gkgeom_release(gkyl_gkgeom *geo)
{
  gkyl_array_release(geo->psiRZ);
  gkyl_free(geo);
}
