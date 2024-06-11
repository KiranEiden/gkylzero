#include <math.h>

#include <gkyl_moment_prim_maxwell_cyl.h>

// Make indexing cleaner with the dir_shuffle
#define ER 0
#define EP 1
#define EZ 2
#define BR 3
#define BP 4
#define BZ 5

double
gkyl_maxwell_cyl_max_abs_speed(double c, const double q[9])
{
  return c / fabs(q[6]);
}

void
gkyl_maxwell_cyl_flux(double c, const double q[9], double flux[9])
{
  double c2 = c*c;
  double r = fabs(q[6]*q[7]*q[8]);
  double r2 = r*r;
  double g_ii = q[6]*q[6];

  flux[ER] = 0.0;
  flux[EP] = c2*q[BZ];
  flux[EZ] = -c2*r2*q[BP]/g_ii;
  flux[BR] = 0.0;
  flux[BP] = -q[EZ];
  flux[BZ] = q[EP]*r2/g_ii;
  flux[6] = flux[7] = flux[8] = 0.0;
}
