#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_calc_prim_vars.h>
#include <gkyl_dg_calc_prim_vars_priv.h>
#include <gkyl_util.h>

void gkyl_calc_prim_vars_u_from_statevec(gkyl_dg_bin_op_mem *mem, struct gkyl_basis basis, const struct gkyl_range *range,
  const struct gkyl_array* statevec, struct gkyl_array* u_i)
{
  // Find number of components of flow vector
  int num_comp = u_i->ncomp/basis.num_basis;
  for (int i = 0; i<num_comp; ++i)
    gkyl_dg_div_op_range(mem, basis, 
      i, u_i, i+1, statevec, 0, statevec, range);  
}

void gkyl_calc_prim_vars_u_from_rhou(gkyl_dg_bin_op_mem *mem, struct gkyl_basis basis, const struct gkyl_range *range,
  const struct gkyl_array* rho, const struct gkyl_array* rhou, struct gkyl_array* u_i)
{
  // Find number of components of flow vector
  int num_comp = u_i->ncomp/basis.num_basis;
  for (int i = 0; i<num_comp; ++i)
    gkyl_dg_div_op_range(mem, basis, 
      i, u_i, i, rhou, 0, rho, range);  
}

void gkyl_calc_prim_vars_p_from_statevec(struct gkyl_basis basis, const struct gkyl_range *range,
  const double p_fac,  const struct gkyl_array* u_i, const struct gkyl_array* statevec,
  struct gkyl_array* p_ij)
{
  int cdim = basis.ndim;
  int poly_order = basis.poly_order;
  euler_pressure_t pressure = choose_ser_euler_pressure_kern(cdim, poly_order);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);

    const double *statevec_d = gkyl_array_cfetch(statevec, loc);
    const double *u_i_d = gkyl_array_cfetch(u_i, loc);
    double *p_ij_d = gkyl_array_fetch(p_ij, loc);
    pressure(p_fac, u_i_d, statevec_d, p_ij_d);
  }
}