/* -*- c++ -*- */

#include "gkyl_util.h"
extern "C" {
#include <assert.h>
#include <string.h>    
    
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_maxwell.h>    
#include <gkyl_dg_maxwell_priv.h>

#include "cart_modal_serendip_priv.h"
}

__global__ void static
gkyl_cart_modal_serendip_cu_dev_kern(struct gkyl_basis *basis, int ndim, int poly_order)
{
  assert(ev_list[ndim].ev[poly_order]);

  basis->ndim = ndim;
  basis->poly_order = poly_order;
  basis->num_basis = num_basis_list[ndim].count[poly_order];
  basis->b_type = GKYL_BASIS_MODAL_SERENDIPITY;  
  
  // function pointers
  basis->eval = ev_list[ndim].ev[poly_order];
  basis->eval_expand = eve_list[ndim].ev[poly_order];
  basis->eval_grad_expand = eveg_list[ndim].ev[poly_order];
  basis->flip_odd_sign = fos_list[ndim].fs[poly_order];
  basis->flip_even_sign = fes_list[ndim].fs[poly_order];
  basis->node_list = nl_list[ndim].nl[poly_order];
  basis->nodal_to_modal = n2m_list[ndim].n2m[poly_order];
}

void
gkyl_cart_modal_serendip_cu_dev(struct gkyl_basis *basis, int ndim, int poly_order)
{
  assert(ndim>0 && ndim<=6);

  struct gkyl_basis ho_basis;

  strcpy(ho_basis.id, "serendipity");
  // this copy needs to be done here as the strcpy needed in the
  // "type" field can't be done on the device
  gkyl_cu_memcpy(basis, &ho_basis, sizeof(struct gkyl_basis),
    GKYL_CU_MEMCPY_H2D);
  
  gkyl_cart_modal_serendip_cu_dev_kern<<<1,1>>>(basis, ndim, poly_order);
}
