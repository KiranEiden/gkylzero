/* -*- c++ -*- */

#include <assert.h>
#include <stdlib.h>
#include <string.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_mom.h>
#include <gkyl_vlasov_mom_priv.h>
}

enum { M0, M1i, M2, M2ij, M3i, M3ijk, BAD };

static int
get_mom_id(const char *mom)
{
  int mom_idx = BAD;

  if (strcmp(mom, "M0") == 0) { // density
    mom_idx = M0;
  }
  else if (strcmp(mom, "M1i") == 0) { // momentum
    mom_idx = M1i;
  }
  else if (strcmp(mom, "M2") == 0) { // energy
    mom_idx = M2;
  }
  else if (strcmp(mom, "M2ij") == 0) { // pressure tensor in lab-frame
    mom_idx = M2ij;    
  }
  else if (strcmp(mom, "M3i") == 0) { // heat-flux vector in lab-frame
    mom_idx = M3i;
  }
  else if (strcmp(mom, "M3ijk") == 0) { // heat-flux tensor in lab-frame
    mom_idx = M3ijk;
  }
  else {
    mom_idx = BAD;
  }    

  return mom_idx;
}

__global__
void vlasov_mom_set_cu_dev_ptrs(struct gkyl_mom_type* momt, int mom_id, int vdim,
  int polyOrder, int tblidx)
{
  int m3ijk_count[] = { 1, 4, 10 };
  
  switch (mom_id) {
    case M0:
      momt->kernel = m0_kernels[tblidx].kernels[polyOrder];
      momt->num_mom = 1;
      break;

    case M1i:
      momt->kernel = m1i_kernels[tblidx].kernels[polyOrder];
      momt->num_mom = vdim;
      break;

    case M2:
      momt->kernel = m2_kernels[tblidx].kernels[polyOrder];
      momt->num_mom = 1;
      break;

    case M2ij:
      momt->kernel = m2ij_kernels[tblidx].kernels[polyOrder];
      momt->num_mom = vdim*(vdim+1)/2;
      break;

    case M3i:
      momt->kernel = m3i_kernels[tblidx].kernels[polyOrder];
      momt->num_mom = vdim;
      break;

    case M3ijk:
      momt->kernel = m3ijk_kernels[tblidx].kernels[polyOrder];
      momt->num_mom = m3ijk_count[vdim-1];
      break;
      
    default: // can't happen
      break;
  }
}

struct gkyl_mom_type*
gkyl_vlasov_mom_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const char *mom)
{
  assert(cbasis->polyOrder == pbasis->polyOrder);
  
  struct gkyl_mom_type *momt = (struct gkyl_mom_type*) gkyl_malloc(sizeof(struct gkyl_mom_type));
  int cdim = momt->cdim = cbasis->ndim;
  int pdim = momt->pdim = pbasis->ndim;
  int vdim = pdim-cdim;
  int polyOrder = momt->polyOrder = cbasis->polyOrder;
  momt->num_config = cbasis->numBasis;
  momt->num_phase = pbasis->numBasis;

  // copy struct to device
  struct gkyl_mom_type *momt_cu = (struct gkyl_mom_type*) gkyl_cu_malloc(sizeof(struct gkyl_mom_type));
  gkyl_cu_memcpy(momt_cu, momt, sizeof(struct gkyl_mom_type), GKYL_CU_MEMCPY_H2D);

  assert(cv_index[cdim].vdim[vdim] != -1);
  int mom_id = get_mom_id(mom);
  assert(mom_id != BAD);

  vlasov_mom_set_cu_dev_ptrs<<<1,1>>>(momt_cu, mom_id, vdim, polyOrder, cv_index[cdim].vdim[vdim]);
  
  gkyl_free(momt);
    
  return momt_cu;
}
