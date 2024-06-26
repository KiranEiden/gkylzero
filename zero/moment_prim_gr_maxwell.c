#include <math.h>

#include <gkyl_moment_prim_maxwell.h>

// Make indexing cleaner with the dir_shuffle
#define D1 0
#define D2 1
#define D3 2
#define B1 3
#define B2 4
#define B3 5

double
gkyl_gr_maxwell_max_abs_speed(const double q[10])
{
  return fmax(fabs(q[6] - q[7]), fabs(q[6] + q[7]));
}

void
gkyl_gr_maxwell_flux(const double q[10], double flux[10])
{
  const double lapse = q[6];
  const double shift1 = q[7];
  const double shift2 = q[8];
  const double shift3 = q[9];
  
  flux[D1] = 0.0;
  flux[D2] = lapse*q[B3] + (shift2*q[D1] - shift1*q[D2]);
  flux[D3] = -(lapse*q[B2] + (shift1*q[D3] - shift3*q[D1]));
  flux[B1] = 0.0;
  flux[B2] = -(lapse*q[D3] + (shift1*q[B2] - shift2*q[B1]));
  flux[B3] = lapse*q[D2] + (shift3*q[B1] - shift1*q[B3]);
  flux[6] = flux[7] = flux[8] = flux[9] = 0.0;
}
