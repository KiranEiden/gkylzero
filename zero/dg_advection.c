#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_advection.h>
#include <gkyl_dg_advection_priv.h>
#include <gkyl_util.h>

// "Choose Kernel" based on cdim and polyorder
#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order]

void 
gkyl_advection_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_eqn *base = container_of(ref, struct gkyl_dg_eqn, ref_count);

  if (gkyl_dg_eqn_is_cu_dev(base)) {
    // free inner on_dev object
    struct dg_advection *advection = container_of(base->on_dev, struct dg_advection, eqn);
    gkyl_cu_free(advection);
  }  
  
  struct dg_advection *advection = container_of(base, struct dg_advection, eqn);
  gkyl_free(advection);
}

struct gkyl_array_copy_func*
gkyl_advection_absorb_bc_create(const struct gkyl_dg_eqn *eqn, int dir, const struct gkyl_basis* cbasis)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_dg_eqn_is_cu_dev(eqn)) {
    return gkyl_advection_absorb_bc_create_cu(eqn->on_dev, dir, cbasis);
  }
#endif

  struct dg_advection *advection = container_of(eqn, struct dg_advection, eqn);

  struct dg_bc_ctx *ctx = (struct dg_bc_ctx*) gkyl_malloc(sizeof(struct dg_bc_ctx));
  ctx->basis = cbasis;

  struct gkyl_array_copy_func *bc = (struct gkyl_array_copy_func*) gkyl_malloc(sizeof(struct gkyl_array_copy_func));
  bc->func = advection->absorb_bc;
  bc->ctx = ctx;
  bc->ctx_on_dev = bc->ctx;

  bc->flags = 0;
  GKYL_CLEAR_CU_ALLOC(bc->flags);
  bc->on_dev = bc; // CPU eqn obj points to itself
  return bc;
}

void
gkyl_advection_bc_release(struct gkyl_array_copy_func* bc)
{
  if (gkyl_array_copy_func_is_cu_dev(bc)) {
    gkyl_cu_free(bc->ctx_on_dev);
    gkyl_cu_free(bc->on_dev);
  }
  gkyl_free(bc->ctx);
  gkyl_free(bc);
}

void
gkyl_advection_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_advection_auxfields auxin)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(auxin.u)) {
    gkyl_advection_set_auxfields_cu(eqn->on_dev, auxin);
    return;
  }
#endif

  struct dg_advection *advection = container_of(eqn, struct dg_advection, eqn);
  advection->auxfields.u = auxin.u;
}

struct gkyl_dg_eqn*
gkyl_dg_advection_new(const struct gkyl_basis* cbasis, const struct gkyl_range* conf_range, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_dg_advection_cu_dev_new(cbasis, conf_range);
  } 
#endif
  struct dg_advection *advection = gkyl_malloc(sizeof(struct dg_advection));

  int cdim = cbasis->ndim;
  int poly_order = cbasis->poly_order;

  const gkyl_dg_advection_vol_kern_list *vol_kernels;
  const gkyl_dg_advection_surf_kern_list *surf_x_kernels, *surf_y_kernels, *surf_z_kernels;

  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      vol_kernels = ser_vol_kernels;
      surf_x_kernels = ser_surf_x_kernels;
      surf_y_kernels = ser_surf_y_kernels;
      surf_z_kernels = ser_surf_z_kernels;
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      vol_kernels = ten_vol_kernels;
      surf_x_kernels = ten_surf_x_kernels;
      surf_y_kernels = ten_surf_y_kernels;
      surf_z_kernels = ten_surf_z_kernels;
      break;

    default:
      assert(false);
      break;    
  }  
    
  advection->eqn.num_equations = 1;
  advection->eqn.vol_term = vol;
  advection->eqn.surf_term = surf;
  advection->eqn.boundary_surf_term = boundary_surf;

  advection->vol =  CK(vol_kernels, cdim, poly_order);
  assert(advection->vol);

  advection->surf[0] = CK(surf_x_kernels, cdim, poly_order);
  if (cdim>1)
    advection->surf[1] = CK(surf_y_kernels, cdim, poly_order);
  if (cdim>2)
    advection->surf[2] = CK(surf_z_kernels, cdim, poly_order);

  // setup pointer for absorbing BC function
  advection->absorb_bc = advection_absorb_bc;

  // ensure non-NULL pointers 
  for (int i=0; i<cdim; ++i) assert(advection->surf[i]);

  advection->auxfields.u = 0;  
  advection->conf_range = *conf_range;

  advection->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(advection->eqn.flags);
  advection->eqn.ref_count = gkyl_ref_count_init(gkyl_advection_free);
  advection->eqn.on_dev = &advection->eqn; // CPU eqn obj points to itself
  
  return &advection->eqn;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_eqn*
gkyl_dg_advection_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_range* conf_range)
{
  assert(false);
  return 0;
}

#endif