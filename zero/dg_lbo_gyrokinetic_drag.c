#include <assert.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_dg_lbo_gyrokinetic_drag.h>
#include <gkyl_dg_lbo_gyrokinetic_drag_priv.h>
#include <gkyl_util.h>

void
gkyl_lbo_gyrokinetic_drag_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_dg_eqn* base = container_of(ref, struct gkyl_dg_eqn, ref_count);
  struct dg_lbo_gyrokinetic_drag *lbo_gyrokinetic_drag  = container_of(base, struct dg_lbo_gyrokinetic_drag, eqn);

  if (GKYL_IS_CU_ALLOC(lbo_gyrokinetic_drag->eqn.flags))
    gkyl_cu_free(lbo_gyrokinetic_drag->eqn.on_dev);
  
  gkyl_free(lbo_gyrokinetic_drag);
}

void
gkyl_lbo_gyrokinetic_drag_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_lbo_gyrokinetic_drag_auxfields auxin)
{

#ifdef GKYL_HAVE_CUDA
 if (gkyl_array_is_cu_dev(auxin.bmag_inv) && gkyl_array_is_cu_dev(auxin.nuSum) &&
     gkyl_array_is_cu_dev(auxin.nuUSum) && gkyl_array_is_cu_dev(auxin.nuVtSqSum)) {
   gkyl_lbo_gyrokinetic_drag_set_auxfields_cu(eqn->on_dev, auxin);
   return;
 }
#endif

  struct dg_lbo_gyrokinetic_drag *lbo_gyrokinetic_drag = container_of(eqn, struct dg_lbo_gyrokinetic_drag, eqn);
  lbo_gyrokinetic_drag->auxfields.bmag_inv = auxin.bmag_inv;
  lbo_gyrokinetic_drag->auxfields.nuSum = auxin.nuSum;
  lbo_gyrokinetic_drag->auxfields.nuUSum = auxin.nuUSum;
  lbo_gyrokinetic_drag->auxfields.nuVtSqSum = auxin.nuVtSqSum;
}


struct gkyl_dg_eqn*
gkyl_dg_lbo_gyrokinetic_drag_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, double mass)
{
  struct dg_lbo_gyrokinetic_drag* lbo_gyrokinetic_drag = gkyl_malloc(sizeof(struct dg_lbo_gyrokinetic_drag));

  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  lbo_gyrokinetic_drag->cdim = cdim;
  lbo_gyrokinetic_drag->pdim = pdim;

  lbo_gyrokinetic_drag->eqn.num_equations = 1;
  lbo_gyrokinetic_drag->eqn.vol_term = vol;
  lbo_gyrokinetic_drag->eqn.surf_term = surf;
  lbo_gyrokinetic_drag->eqn.boundary_surf_term = boundary_surf;

  const gkyl_dg_lbo_gyrokinetic_drag_vol_kern_list *vol_kernels;
  const gkyl_dg_lbo_gyrokinetic_drag_surf_kern_list *surf_vpar_kernels, *surf_mu_kernels;
  const gkyl_dg_lbo_gyrokinetic_drag_boundary_surf_kern_list *boundary_surf_vpar_kernels, *boundary_surf_mu_kernels;
  
  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      vol_kernels = ser_vol_kernels;
      surf_vpar_kernels = ser_surf_vpar_kernels;
      surf_mu_kernels = ser_surf_mu_kernels;
      boundary_surf_vpar_kernels = ser_boundary_surf_vpar_kernels;
      boundary_surf_mu_kernels = ser_boundary_surf_mu_kernels;
      break;

    default:
      assert(false);
      break;    
  }  

  lbo_gyrokinetic_drag->vol = CK(vol_kernels, cdim, vdim, poly_order);

  lbo_gyrokinetic_drag->surf[0] = CK(surf_vpar_kernels, cdim, vdim, poly_order);
  if (vdim>1)
    lbo_gyrokinetic_drag->surf[1] = CK(surf_mu_kernels, cdim, vdim, poly_order);

  lbo_gyrokinetic_drag->boundary_surf[0] = CK(boundary_surf_vpar_kernels, cdim, vdim, poly_order);
  if (vdim>1)
    lbo_gyrokinetic_drag->boundary_surf[1] = CK(boundary_surf_mu_kernels, cdim, vdim, poly_order);

  // ensure non-NULL pointers
  assert(lbo_gyrokinetic_drag->vol);
  for (int i=0; i<vdim; ++i) assert(lbo_gyrokinetic_drag->surf[i]);
  for (int i=0; i<vdim; ++i) assert(lbo_gyrokinetic_drag->boundary_surf[i]);

  lbo_gyrokinetic_drag->mass = mass;
  lbo_gyrokinetic_drag->auxfields.bmag_inv = 0;
  lbo_gyrokinetic_drag->auxfields.nuSum = 0;
  lbo_gyrokinetic_drag->auxfields.nuUSum = 0;
  lbo_gyrokinetic_drag->auxfields.nuVtSqSum = 0;
  lbo_gyrokinetic_drag->conf_range = *conf_range;

  lbo_gyrokinetic_drag->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(lbo_gyrokinetic_drag->eqn.flags);
  lbo_gyrokinetic_drag->eqn.ref_count = gkyl_ref_count_init(gkyl_lbo_gyrokinetic_drag_free);
  lbo_gyrokinetic_drag->eqn.on_dev = &lbo_gyrokinetic_drag->eqn;
  
  return &lbo_gyrokinetic_drag->eqn;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_eqn*
gkyl_dg_lbo_gyrokinetic_drag_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, double mass)
{
  assert(false);
  return 0;
}

#endif