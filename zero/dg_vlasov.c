#include "gkyl_dg_eqn.h"
#include <assert.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_dg_vlasov.h>
#include <gkyl_dg_vlasov_priv.h>
#include <gkyl_util.h>

void
gkyl_vlasov_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_eqn *base = container_of(ref, struct gkyl_dg_eqn, ref_count);
  
  if (gkyl_dg_eqn_is_cu_dev(base)) {
    // free inner on_dev object
    struct dg_vlasov *vlasov = container_of(base->on_dev, struct dg_vlasov, eqn);
    gkyl_cu_free(vlasov);
  }
  
  struct dg_vlasov *vlasov = container_of(base, struct dg_vlasov, eqn);
  gkyl_free(vlasov);
}

void
gkyl_vlasov_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_vlasov_auxfields auxin)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(auxin.qmem)) {
    gkyl_vlasov_set_auxfields_cu(eqn->on_dev, auxin);
    return;
  }
#endif

  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);
  vlasov->auxfields.qmem = auxin.qmem;
}

struct gkyl_dg_eqn*
gkyl_dg_vlasov_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  const struct gkyl_range* conf_range, enum gkyl_field_id field_id, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_dg_vlasov_cu_dev_new(cbasis, pbasis, conf_range, field_id);
  } 
#endif
  struct dg_vlasov *vlasov = gkyl_malloc(sizeof(struct dg_vlasov));


  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  vlasov->cdim = cdim;
  vlasov->pdim = pdim;

  vlasov->eqn.num_equations = 1;
  vlasov->eqn.vol_term = vol;
  vlasov->eqn.surf_term = surf;
  vlasov->eqn.boundary_surf_term = boundary_surf;

  const gkyl_dg_vlasov_stream_vol_kern_list *stream_vol_kernels;
  const gkyl_dg_vlasov_vol_kern_list *vol_kernels;
  const gkyl_dg_vlasov_stream_gen_geo_vol_kern_list *stream_gen_geo_vol_kernels; 
  const gkyl_dg_vlasov_stream_surf_kern_list *stream_surf_x_kernels, *stream_surf_y_kernels, *stream_surf_z_kernels;
  const gkyl_dg_vlasov_accel_surf_kern_list *accel_surf_vx_kernels, *accel_surf_vy_kernels, *accel_surf_vz_kernels;
  const gkyl_dg_vlasov_accel_boundary_surf_kern_list *accel_boundary_surf_vx_kernels, *accel_boundary_surf_vy_kernels,
    *accel_boundary_surf_vz_kernels;
  
  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      stream_vol_kernels = ser_stream_vol_kernels;
      stream_gen_geo_vol_kernels = ser_stream_gen_geo_vol_kernels;
      vol_kernels = ser_vol_kernels;
      stream_surf_x_kernels = ser_stream_surf_x_kernels;
      stream_surf_y_kernels = ser_stream_surf_y_kernels;
      stream_surf_z_kernels = ser_stream_surf_z_kernels;
      accel_surf_vx_kernels = ser_accel_surf_vx_kernels;
      accel_surf_vy_kernels = ser_accel_surf_vy_kernels;
      accel_surf_vz_kernels = ser_accel_surf_vz_kernels;
      accel_boundary_surf_vx_kernels = ser_accel_boundary_surf_vx_kernels;
      accel_boundary_surf_vy_kernels = ser_accel_boundary_surf_vy_kernels;
      accel_boundary_surf_vz_kernels = ser_accel_boundary_surf_vz_kernels;
      
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      stream_vol_kernels = ten_stream_vol_kernels;
      vol_kernels = ten_vol_kernels;
      stream_surf_x_kernels = ten_stream_surf_x_kernels;
      stream_surf_y_kernels = ten_stream_surf_y_kernels;
      stream_surf_z_kernels = ten_stream_surf_z_kernels;
      accel_surf_vx_kernels = ten_accel_surf_vx_kernels;
      accel_surf_vy_kernels = ten_accel_surf_vy_kernels;
      accel_surf_vz_kernels = ten_accel_surf_vz_kernels;
      accel_boundary_surf_vx_kernels = ten_accel_boundary_surf_vx_kernels;
      accel_boundary_surf_vy_kernels = ten_accel_boundary_surf_vy_kernels;
      accel_boundary_surf_vz_kernels = ten_accel_boundary_surf_vz_kernels;
      break;

    default:
      assert(false);
      break;    
  }  
  if (field_id == GKYL_FIELD_NULL)
    vlasov->vol = CK(stream_vol_kernels,cdim,vdim,poly_order);
  else if (field_id == GKYL_FIELD_NULL_GEN_GEO)
    vlasov->vol = CK(stream_gen_geo_vol_kernels,cdim,vdim,poly_order);
  else
    vlasov->vol = CK(vol_kernels,cdim,vdim,poly_order);

  vlasov->stream_surf[0] = CK(stream_surf_x_kernels,cdim,vdim,poly_order);
  if (cdim>1)
    vlasov->stream_surf[1] = CK(stream_surf_y_kernels,cdim,vdim,poly_order);
  if (cdim>2)
    vlasov->stream_surf[2] = CK(stream_surf_z_kernels,cdim,vdim,poly_order);

  vlasov->accel_surf[0] = CK(accel_surf_vx_kernels,cdim,vdim,poly_order);
  if (vdim>1)
    vlasov->accel_surf[1] = CK(accel_surf_vy_kernels,cdim,vdim,poly_order);
  if (vdim>2)
    vlasov->accel_surf[2] = CK(accel_surf_vz_kernels,cdim,vdim,poly_order);

  vlasov->accel_boundary_surf[0] = CK(accel_boundary_surf_vx_kernels,cdim,vdim,poly_order);
  if (vdim>1)
    vlasov->accel_boundary_surf[1] = CK(accel_boundary_surf_vy_kernels,cdim,vdim,poly_order);
  if (vdim>2)
    vlasov->accel_boundary_surf[2] = CK(accel_boundary_surf_vz_kernels,cdim,vdim,poly_order);

  // ensure non-NULL pointers
  assert(vlasov->vol);
  for (int i=0; i<cdim; ++i) assert(vlasov->stream_surf[i]);
  for (int i=0; i<vdim; ++i) assert(vlasov->accel_surf[i]);
  for (int i=0; i<vdim; ++i) assert(vlasov->accel_boundary_surf[i]);

  vlasov->auxfields.qmem = 0;  
  vlasov->conf_range = *conf_range;

  vlasov->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(vlasov->eqn.flags);

  vlasov->eqn.ref_count = gkyl_ref_count_init(gkyl_vlasov_free);
  vlasov->eqn.on_dev = &vlasov->eqn; // CPU eqn obj points to itself
  
  return &vlasov->eqn;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_eqn*
gkyl_dg_vlasov_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range, enum gkyl_field_id field_id)
{
  assert(false);
  return 0;
}

#endif
