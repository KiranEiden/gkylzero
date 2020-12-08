#include <gkyl_mom_type.h>

struct gkyl_mom_type*
gkyl_mom_type_aquire(const struct gkyl_mom_type* momt)
{
  gkyl_ref_count_inc(&momt->ref_count);
  return (struct gkyl_mom_type*) momt;
}

void
gkyl_mom_type_release(struct gkyl_mom_type* momt)
{
  gkyl_ref_count_dec(&momt->ref_count);
}

void
gkyl_mom_type_calc(const struct gkyl_mom_type* momt,
  const double *xc, const double *dx, const long *idx,
  const double *f, double* out)
{
  momt->kernel(xc, dx, idx, f, out);
}
