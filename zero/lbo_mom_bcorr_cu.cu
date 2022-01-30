/* -*- c++ -*- */

#include <assert.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_util.h>
#include <gkyl_lbo_mom_bcorr.h>
#include <gkyl_lbo_mom_bcorr_priv.h>
}

enum { f, vf, BAD };

static int
get_mom_id(const char *mom)
{
  int mom_idx = BAD;

  if (strcmp(mom, "f") == 0) { // density
    mom_idx = f;
  }
  else if (strcmp(mom, "vf") == 0) { // momentum
    mom_idx = vf;
  }
  else {
    mom_idx = BAD;
  }

  return mom_idx;
}

__global__
static void
lbo_mom_set_cu_dev_ptrs(struct lbo_mom_type* mom_bcorr, int mom_id, enum gkyl_basis_type b_type, int vdim, int poly_order, int tblidx)
{
  mom_bcorr->momt.kernel = kernel;

  // choose kernel tables based on basis-function type
  const gkyl_lbo_mom_kern_list *boundary_integral_f_kernels, *boundary_integral_vf_kernels;

  switch (b_type) {
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

  switch (mom_id) {
  case f:
    mom_bcorr->kernel = boundary_integral_f_kernels[tblidx].kernels[poly_order];
    mom_bcorr->momt.num_mom = vdim;
    break;

  case vf:
    mom_bcorr->kernel = boundary_integral_vf_kernels[tblidx].kernels[poly_order];
    mom_bcorr->momt.num_mom = vdim;
    break;

  default: // can't happen
    break;
  }
}

struct gkyl_mom_type*
gkyl_vlasov_lbo_mom_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const char *mom, const double *vBoundary)
{
  assert(cbasis->poly_order == pbasis->poly_order);

  struct lbo_mom_type *mom_bcorr = (struct lbo_mom_type*) gkyl_malloc(sizeof(struct lbo_mom_type));

  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  mom_bcorr->momt.cdim = cdim;
  mom_bcorr->momt.pdim = pdim;
  mom_bcorr->momt.poly_order = poly_order;
  mom_bcorr->momt.num_config = cbasis->num_basis;
  mom_bcorr->momt.num_phase = pbasis->num_basis;
  for (int d=0; d<vdim; ++d) {
    mom_bcorr->vBoundary[d] = vBoundary[d];
    mom_bcorr->vBoundary[d + vdim] = vBoundary[d + vdim];
  }

  int mom_id = get_mom_id(mom);
  assert(mom_id != BAD);
  mom_bcorr->momt.num_mom = vdim; // number of moments

  mom_bcorr->momt.flag = 0;
  GKYL_SET_CU_ALLOC(mom_bcorr->momt.flag);
  mom_bcorr->momt.ref_count = gkyl_ref_count_init(gkyl_mom_free);

  // copy struct to device
  struct lbo_mom_type *mom_bcorr_cu = (struct lbo_mom_type*)
    gkyl_cu_malloc(sizeof(struct lbo_mom_type));
  gkyl_cu_memcpy(mom_bcorr_cu, mom_bcorr, sizeof(struct lbo_mom_type), GKYL_CU_MEMCPY_H2D);

  assert(cv_index[cdim].vdim[vdim] != -1);


  lbo_mom_set_cu_dev_ptrs<<<1,1>>>(mom_bcorr_cu, mom_id, cbasis->b_type,
    vdim, poly_order, cv_index[cdim].vdim[vdim]);

  mom_bcorr->momt.on_dev = &mom_bcorr_cu->momt;

  return &mom_bcorr->momt;
}
