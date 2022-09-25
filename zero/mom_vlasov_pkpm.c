#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_mom_vlasov_pkpm.h>
#include <gkyl_mom_vlasov_pkpm_priv.h>
#include <gkyl_util.h>

// "Choose Kernel" based on cdim and polyorder
#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order]

void
gkyl_vlasov_pkpm_mom_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_mom_type *momt = container_of(ref, struct gkyl_mom_type, ref_count);
  if (GKYL_IS_CU_ALLOC(momt->flags))
    gkyl_cu_free(momt->on_dev);
  gkyl_free(momt);
}

struct gkyl_mom_type*
gkyl_mom_vlasov_pkpm_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  double mass, bool use_gpu)
{
  assert(cbasis->poly_order == pbasis->poly_order);

#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_mom_vlasov_pkpm_cu_dev_new(cbasis, pbasis, conf_range, mass);
  } 
#endif    
  struct mom_type_vlasov_pkpm *mom_vlasov_pkpm = gkyl_malloc(sizeof(struct mom_type_vlasov_pkpm));
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  mom_vlasov_pkpm->momt.cdim = cdim;
  mom_vlasov_pkpm->momt.pdim = pdim;
  mom_vlasov_pkpm->momt.poly_order = poly_order;
  mom_vlasov_pkpm->momt.num_config = cbasis->num_basis;
  mom_vlasov_pkpm->momt.num_phase = pbasis->num_basis;

  // choose kernel tables based on basis-function type
  const gkyl_mom_vlasov_pkpm_kern_list *mom_vlasov_pkpm_kernels;

  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      mom_vlasov_pkpm_kernels = ser_mom_vlasov_pkpm_kernels;
      break;

    default:
      assert(false);
      break;    
  }

  mom_vlasov_pkpm->momt.kernel = CK(mom_vlasov_pkpm_kernels, cdim, poly_order);
  mom_vlasov_pkpm->momt.num_mom = 3; // rho, p_parallel, q_parallel

  mom_vlasov_pkpm->mass = mass;
    
  mom_vlasov_pkpm->momt.flags = 0;
  GKYL_CLEAR_CU_ALLOC(mom_vlasov_pkpm->momt.flags);
  mom_vlasov_pkpm->momt.ref_count = gkyl_ref_count_init(gkyl_vlasov_pkpm_mom_free);
  
  mom_vlasov_pkpm->momt.on_dev = &mom_vlasov_pkpm->momt; // on host, self-reference
    
  return &mom_vlasov_pkpm->momt;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_mom_type*
gkyl_mom_vlasov_pkpm_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  double mass)
{
  assert(false);
  return 0;
}

#endif
