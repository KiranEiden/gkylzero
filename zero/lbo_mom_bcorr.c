#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_util.h>
#include <gkyl_lbo_mom_bcorr.h>
#include <gkyl_lbo_mom_bcorr_priv.h>

static void
mom_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_mom_type *base = container_of(ref, struct gkyl_mom_type, ref_count);
  struct lbo_mom_type *mom_bcorr = container_of(base, struct lbo_mom_type, momt);
  gkyl_free(mom_bcorr);
}

void
gkyl_lbo_mom_set_vBoundary(const struct gkyl_mom_type *momt, const double vBoundary)
{
  struct lbo_mom_type *mom_bcorr = container_of(momt, struct lbo_mom_type, momt);
  mom_bcorr->vBoundary = vBoundary;
}

void
gkyl_lbo_mom_set_atLower(const struct gkyl_mom_type *momt, const bool atLower)
{
  struct lbo_mom_type *mom_bcorr = container_of(momt, struct lbo_mom_type, momt);
  mom_bcorr->atLower = atLower;
}

struct gkyl_mom_type*
gkyl_vlasov_lbo_mom_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const char *mom)
{
  assert(cbasis->poly_order == pbasis->poly_order);
  
  struct lbo_mom_type *mom_bcorr = gkyl_malloc(sizeof(struct lbo_mom_type));
  int cdim = mom_bcorr->cdim = cbasis->ndim;
  int pdim = mom_bcorr->pdim = pbasis->ndim;
  int vdim = pdim-cdim;
  int poly_order = mom_bcorr->poly_order = cbasis->poly_order;
  mom_bcorr->num_config = cbasis->num_basis;
  mom_bcorr->num_phase = pbasis->num_basis;

  mom_bcorr->momt.kernel = kernel;

  // choose kernel tables based on basis-function type
  const gkyl_lbo_mom_kern_list *boundary_integral_f_kernels, *boundary_integral_vf_kernels;

  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      boundary_integral_f_kernels = ser_boundary_integral_f_kernels;
      boundary_integral_vf_kernels = ser_boundary_integral_vf_kernels;
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      break;

    default:
      assert(false);
      break;    
  }
  if (strcmp(mom, "f") == 0) { // density
    printf("Gimme a v!\n");
    fflush(stdout);
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != boundary_integral_f_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_bcorr->kernel = boundary_integral_f_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
  }
  else if (strcmp(mom, "vf") == 0) { // momentum
    printf("Gimme a vf!\n");
    fflush(stdout);
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != boundary_integral_vf_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_bcorr->kernel = boundary_integral_vf_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
  }
  mom_bcorr->num_mom = vdim;
  if (mom_bcorr->kernel) {
    printf("%s: Success!\n", mom);
    fflush(stdout);
  }

  // set reference counter
  mom_bcorr->ref_count = (struct gkyl_ref_count) { mom_free, 1 };
    
  return &mom_bcorr->momt;
}
