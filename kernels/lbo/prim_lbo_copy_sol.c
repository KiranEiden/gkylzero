#include <gkyl_mat.h> 
#include <gkyl_prim_lbo_vlasov_kernels.h>
 
GKYL_CU_DH void prim_lbo_copy_sol(const struct gkyl_mat *rhs, const int nc, const int vdim, double* GKYL_RESTRICT u, double* GKYL_RESTRICT vtSq)
{
  for(size_t i=0; i<nc; i++) { 
    u[i] = gkyl_mat_get(rhs, i, 0);
    vtSq[i] = gkyl_mat_get(rhs, i+nc*vdim, 0);
  }
}