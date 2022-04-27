/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_vlasov_sr.h>    
#include <gkyl_dg_vlasov_sr_priv.h>
}

#include <cassert>

__global__ static void
gkyl_vlasov_sr_wall_bc_create_set_cu_dev_ptrs(const struct gkyl_dg_eqn *eqn, int dir,
  const struct gkyl_basis* pbasis, struct species_wall_bc_ctx *ctx, struct gkyl_array_copy_func *bc)
{
  struct dg_vlasov_sr *vlasov_sr = container_of(eqn, struct dg_vlasov_sr, eqn);

  ctx->dir = dir;
  ctx->cdim = vlasov_sr->cdim;
  ctx->basis = pbasis;

  bc->func = vlasov_sr->wall_bc;
  bc->ctx = ctx;
}

struct gkyl_array_copy_func*
gkyl_vlasov_sr_wall_bc_create_cu(const struct gkyl_dg_eqn *eqn, int dir, const struct gkyl_basis* pbasis)
{
  // create host context and bc func structs
  struct species_wall_bc_ctx *ctx = (struct species_wall_bc_ctx*) gkyl_malloc(sizeof(struct species_wall_bc_ctx));
  struct gkyl_array_copy_func *bc = (struct gkyl_array_copy_func*) gkyl_malloc(sizeof(struct gkyl_array_copy_func));
  bc->ctx = ctx;

  bc->flags = 0;
  GKYL_SET_CU_ALLOC(bc->flags);

  // create device context and bc func structs
  struct species_wall_bc_ctx *ctx_cu = (struct species_wall_bc_ctx*) gkyl_cu_malloc(sizeof(struct species_wall_bc_ctx));
  struct gkyl_array_copy_func *bc_cu = (struct gkyl_array_copy_func*) gkyl_cu_malloc(sizeof(struct gkyl_array_copy_func));

  gkyl_cu_memcpy(ctx_cu, ctx, sizeof(struct species_wall_bc_ctx), GKYL_CU_MEMCPY_H2D);
  gkyl_cu_memcpy(bc_cu, bc, sizeof(struct gkyl_array_copy_func), GKYL_CU_MEMCPY_H2D);

  gkyl_vlasov_sr_wall_bc_create_set_cu_dev_ptrs<<<1,1>>>(eqn, dir, pbasis, ctx_cu, bc_cu);

  // set parent on_dev pointer 
  bc->on_dev = bc_cu;  
  return bc;
}

// CUDA kernel to set pointer to auxiliary fields.
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_vlasov_sr_set_auxfields_cu_kernel(const struct gkyl_dg_eqn *eqn, 
  const struct gkyl_array *qmem, const struct gkyl_array *p_over_gamma)
{
  struct dg_vlasov_sr *vlasov_sr = container_of(eqn, struct dg_vlasov_sr, eqn);
  vlasov_sr->auxfields.qmem = qmem;
  vlasov_sr->auxfields.p_over_gamma = p_over_gamma;
}

// Host-side wrapper for set_auxfields_cu_kernel
void
gkyl_vlasov_sr_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_vlasov_sr_auxfields auxin)
{
  gkyl_vlasov_sr_set_auxfields_cu_kernel<<<1,1>>>(eqn, auxin.qmem->on_dev, auxin.p_over_gamma->on_dev);
}

// CUDA kernel to set device pointers to range object and vlasov kernel function
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
dg_vlasov_sr_set_cu_dev_ptrs(struct dg_vlasov_sr *vlasov_sr, enum gkyl_basis_type b_type,
  int cv_index, int cdim, int vdim, int poly_order, enum gkyl_field_id field_id)
{
  vlasov_sr->auxfields.qmem = 0; 
  vlasov_sr->auxfields.p_over_gamma= 0; 

  vlasov_sr->eqn.surf_term = surf;
  vlasov_sr->eqn.boundary_surf_term = boundary_surf;

  const gkyl_dg_vlasov_sr_stream_vol_kern_list *stream_vol_kernels;
  const gkyl_dg_vlasov_sr_vol_kern_list *vol_kernels;
  const gkyl_dg_vlasov_sr_stream_surf_kern_list *stream_surf_x_kernels, *stream_surf_y_kernels, *stream_surf_z_kernels;
  const gkyl_dg_vlasov_sr_accel_surf_kern_list *accel_surf_vx_kernels, *accel_surf_vy_kernels, *accel_surf_vz_kernels;
  const gkyl_dg_vlasov_sr_accel_boundary_surf_kern_list *accel_boundary_surf_vx_kernels, *accel_boundary_surf_vy_kernels,
    *accel_boundary_surf_vz_kernels;
  
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      stream_vol_kernels = ser_stream_vol_kernels;
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

    // case GKYL_BASIS_MODAL_TENSOR:
    //   stream_vol_kernels = ten_stream_vol_kernels;
    //   vol_kernels = ten_vol_kernels;
    //   stream_surf_x_kernels = ten_stream_surf_x_kernels;
    //   stream_surf_y_kernels = ten_stream_surf_y_kernels;
    //   stream_surf_z_kernels = ten_stream_surf_z_kernels;
    //   accel_surf_vx_kernels = ten_accel_surf_vx_kernels;
    //   accel_surf_vy_kernels = ten_accel_surf_vy_kernels;
    //   accel_surf_vz_kernels = ten_accel_surf_vz_kernels;
    //   accel_boundary_surf_vx_kernels = ten_accel_boundary_surf_vx_kernels;
    //   accel_boundary_surf_vy_kernels = ten_accel_boundary_surf_vy_kernels;
    //   accel_boundary_surf_vz_kernels = ten_accel_boundary_surf_vz_kernels;
    //   break;

    default:
      assert(false);
      break;    
  }  
  if (field_id == GKYL_FIELD_SR_NULL)
    vlasov_sr->eqn.vol_term = stream_vol_kernels[cv_index].kernels[poly_order];
  else
    vlasov_sr->eqn.vol_term = vol_kernels[cv_index].kernels[poly_order];

  vlasov_sr->stream_surf[0] = stream_surf_x_kernels[cv_index].kernels[poly_order];
  if (cdim>1)
    vlasov_sr->stream_surf[1] = stream_surf_y_kernels[cv_index].kernels[poly_order];
  if (cdim>2)
    vlasov_sr->stream_surf[2] = stream_surf_z_kernels[cv_index].kernels[poly_order];

  vlasov_sr->accel_surf[0] = accel_surf_vx_kernels[cv_index].kernels[poly_order];
  if (vdim>1)
    vlasov_sr->accel_surf[1] = accel_surf_vy_kernels[cv_index].kernels[poly_order];
  if (vdim>2)
    vlasov_sr->accel_surf[2] = accel_surf_vz_kernels[cv_index].kernels[poly_order];

  vlasov_sr->accel_boundary_surf[0] = accel_boundary_surf_vx_kernels[cv_index].kernels[poly_order];
  if (vdim>1)
    vlasov_sr->accel_boundary_surf[1] = accel_boundary_surf_vy_kernels[cv_index].kernels[poly_order];
  if (vdim>2)
    vlasov_sr->accel_boundary_surf[2] = accel_boundary_surf_vz_kernels[cv_index].kernels[poly_order];

  // setup pointer for wall BC function
  vlasov_sr->wall_bc = species_wall_bc;
}

struct gkyl_dg_eqn*
gkyl_dg_vlasov_sr_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  const struct gkyl_range* conf_range, const struct gkyl_range* vel_range, 
  enum gkyl_field_id field_id)
{
  struct dg_vlasov_sr *vlasov_sr = (struct dg_vlasov_sr*) gkyl_malloc(sizeof(struct dg_vlasov_sr));

  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  vlasov_sr->cdim = cdim;
  vlasov_sr->pdim = pdim;

  vlasov_sr->eqn.num_equations = 1;
  vlasov_sr->conf_range = *conf_range;
  vlasov_sr->vel_range = *vel_range;

  vlasov_sr->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(vlasov_sr->eqn.flags);
  vlasov_sr->eqn.ref_count = gkyl_ref_count_init(gkyl_vlasov_sr_free);

  // copy the host struct to device struct
  struct dg_vlasov_sr *vlasov_sr_cu = (struct dg_vlasov_sr*) gkyl_cu_malloc(sizeof(struct dg_vlasov_sr));
  gkyl_cu_memcpy(vlasov_sr_cu, vlasov_sr, sizeof(struct dg_vlasov_sr), GKYL_CU_MEMCPY_H2D);

  dg_vlasov_sr_set_cu_dev_ptrs<<<1,1>>>(vlasov_sr_cu, cbasis->b_type, cv_index[cdim].vdim[vdim],
    cdim, vdim, poly_order, field_id);

  // set parent on_dev pointer
  vlasov_sr->eqn.on_dev = &vlasov_sr_cu->eqn;
  
  return &vlasov_sr->eqn;
}
