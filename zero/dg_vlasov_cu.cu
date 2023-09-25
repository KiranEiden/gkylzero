/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_vlasov.h>    
#include <gkyl_dg_vlasov_priv.h>
}

#include <cassert>

// CUDA kernel to set pointer to auxiliary fields.
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_vlasov_set_auxfields_cu_kernel(const struct gkyl_dg_eqn *eqn, 
  const struct gkyl_array *field, const struct gkyl_array *cot_vec, const struct gkyl_array *alpha_geo)
{
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);
  vlasov->auxfields.field = field; // q/m*(E,B) for Maxwell's, q/m*phi for Poisson's (gradient calculated in kernel)
  vlasov->auxfields.cot_vec = cot_vec; // cotangent vectors (e^i) used in volume term if general geometry enabled
  vlasov->auxfields.alpha_geo = alpha_geo; // alpha^i (e^i . alpha) used in surface term if general geometry enabled
}

// Host-side wrapper for set_auxfields_cu_kernel
void
gkyl_vlasov_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_vlasov_auxfields auxin)
{
  gkyl_vlasov_set_auxfields_cu_kernel<<<1,1>>>(eqn, auxin.field->on_dev,
    auxin.cot_vec ? auxin.cot_vec->on_dev : 0,
    auxin.alpha_geo ? auxin.alpha_geo->on_dev : 0);
}

// CUDA kernel to set device pointers to range object and vlasov kernel function
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
dg_vlasov_set_cu_dev_ptrs(struct dg_vlasov *vlasov, enum gkyl_basis_type b_type,
  int cv_index, int cdim, int vdim, int poly_order, 
  enum gkyl_model_id model_id, enum gkyl_field_id field_id)
{
  vlasov->auxfields.field = 0;
  vlasov->auxfields.cot_vec = 0;
  vlasov->auxfields.alpha_geo = 0;

  vlasov->eqn.surf_term = surf;
  vlasov->eqn.boundary_surf_term = boundary_surf;

  const gkyl_dg_vlasov_stream_vol_kern_list *stream_vol_kernels;
  const gkyl_dg_vlasov_stream_gen_geo_vol_kern_list *stream_gen_geo_vol_kernels; 
  const gkyl_dg_vlasov_poisson_vol_kern_list *poisson_vol_kernels;
  const gkyl_dg_vlasov_poisson_extem_vol_kern_list *poisson_extem_vol_kernels;
  const gkyl_dg_vlasov_vol_kern_list *vol_kernels;

  const gkyl_dg_vlasov_stream_surf_kern_list *stream_surf_x_kernels, 
    *stream_surf_y_kernels, 
    *stream_surf_z_kernels;
  const gkyl_dg_vlasov_stream_gen_geo_surf_kern_list *stream_gen_geo_surf_x_kernels, 
    *stream_gen_geo_surf_y_kernels, 
    *stream_gen_geo_surf_z_kernels;

  const gkyl_dg_vlasov_poisson_accel_surf_kern_list *poisson_accel_surf_vx_kernels, 
    *poisson_accel_surf_vy_kernels, 
    *poisson_accel_surf_vz_kernels;
  const gkyl_dg_vlasov_poisson_extem_accel_surf_kern_list *poisson_extem_accel_surf_vx_kernels, 
    *poisson_extem_accel_surf_vy_kernels, 
    *poisson_extem_accel_surf_vz_kernels;
  const gkyl_dg_vlasov_accel_surf_kern_list *accel_surf_vx_kernels, 
    *accel_surf_vy_kernels, 
    *accel_surf_vz_kernels;

  const gkyl_dg_vlasov_stream_boundary_surf_kern_list *stream_boundary_surf_x_kernels, 
    *stream_boundary_surf_y_kernels,
    *stream_boundary_surf_z_kernels;
  
  const gkyl_dg_vlasov_poisson_accel_boundary_surf_kern_list *poisson_accel_boundary_surf_vx_kernels, 
    *poisson_accel_boundary_surf_vy_kernels,
    *poisson_accel_boundary_surf_vz_kernels;
  const gkyl_dg_vlasov_poisson_extem_accel_boundary_surf_kern_list *poisson_extem_accel_boundary_surf_vx_kernels, 
    *poisson_extem_accel_boundary_surf_vy_kernels,
    *poisson_extem_accel_boundary_surf_vz_kernels;
  const gkyl_dg_vlasov_accel_boundary_surf_kern_list *accel_boundary_surf_vx_kernels, 
    *accel_boundary_surf_vy_kernels,
    *accel_boundary_surf_vz_kernels;
  
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      stream_vol_kernels = ser_stream_vol_kernels;
      stream_gen_geo_vol_kernels = ser_stream_gen_geo_vol_kernels;
      poisson_vol_kernels = ser_poisson_vol_kernels;
      poisson_extem_vol_kernels = ser_poisson_extem_vol_kernels;
      vol_kernels = ser_vol_kernels;

      stream_surf_x_kernels = ser_stream_surf_x_kernels;
      stream_surf_y_kernels = ser_stream_surf_y_kernels;
      stream_surf_z_kernels = ser_stream_surf_z_kernels;
      stream_gen_geo_surf_x_kernels = ser_stream_gen_geo_surf_x_kernels;
      stream_gen_geo_surf_y_kernels = ser_stream_gen_geo_surf_y_kernels;
      stream_gen_geo_surf_z_kernels = ser_stream_gen_geo_surf_z_kernels;

      poisson_accel_surf_vx_kernels = ser_poisson_accel_surf_vx_kernels;
      poisson_accel_surf_vy_kernels = ser_poisson_accel_surf_vy_kernels;
      poisson_accel_surf_vz_kernels = ser_poisson_accel_surf_vz_kernels;
      poisson_extem_accel_surf_vx_kernels = ser_poisson_extem_accel_surf_vx_kernels;
      poisson_extem_accel_surf_vy_kernels = ser_poisson_extem_accel_surf_vy_kernels;
      poisson_extem_accel_surf_vz_kernels = ser_poisson_extem_accel_surf_vz_kernels;
      accel_surf_vx_kernels = ser_accel_surf_vx_kernels;
      accel_surf_vy_kernels = ser_accel_surf_vy_kernels;
      accel_surf_vz_kernels = ser_accel_surf_vz_kernels;

      stream_boundary_surf_x_kernels = ser_stream_boundary_surf_x_kernels;
      stream_boundary_surf_y_kernels = ser_stream_boundary_surf_y_kernels;
      stream_boundary_surf_z_kernels = ser_stream_boundary_surf_z_kernels;
      
      poisson_accel_boundary_surf_vx_kernels = ser_poisson_accel_boundary_surf_vx_kernels;
      poisson_accel_boundary_surf_vy_kernels = ser_poisson_accel_boundary_surf_vy_kernels;
      poisson_accel_boundary_surf_vz_kernels = ser_poisson_accel_boundary_surf_vz_kernels;
      poisson_extem_accel_boundary_surf_vx_kernels = ser_poisson_extem_accel_boundary_surf_vx_kernels;
      poisson_extem_accel_boundary_surf_vy_kernels = ser_poisson_extem_accel_boundary_surf_vy_kernels;
      poisson_extem_accel_boundary_surf_vz_kernels = ser_poisson_extem_accel_boundary_surf_vz_kernels;
      accel_boundary_surf_vx_kernels = ser_accel_boundary_surf_vx_kernels;
      accel_boundary_surf_vy_kernels = ser_accel_boundary_surf_vy_kernels;
      accel_boundary_surf_vz_kernels = ser_accel_boundary_surf_vz_kernels;
      
      break;

    default:
      assert(false);
      break;    
  }
  if (model_id == GKYL_MODEL_GEN_GEO) {
    if (field_id == GKYL_FIELD_NULL) 
      vlasov->eqn.vol_term = stream_gen_geo_vol_kernels[cv_index].kernels[poly_order];
    else 
      assert(false); // general geometry not yet working for non-neutral

    // Neutral gen geo will always be 3x3v.
    vlasov->stream_surf[0] = stream_gen_geo_surf_x_kernels[cv_index].kernels[poly_order];
    vlasov->stream_surf[1] = stream_gen_geo_surf_y_kernels[cv_index].kernels[poly_order];
    vlasov->stream_surf[2] = stream_gen_geo_surf_z_kernels[cv_index].kernels[poly_order];

    // Set acceleration kernels even if they are not used
    vlasov->accel_surf[0] = accel_surf_vx_kernels[cv_index].kernels[poly_order];
    vlasov->accel_surf[1] = accel_surf_vy_kernels[cv_index].kernels[poly_order];
    vlasov->accel_surf[2] = accel_surf_vz_kernels[cv_index].kernels[poly_order];

    vlasov->accel_boundary_surf[0] = accel_boundary_surf_vx_kernels[cv_index].kernels[poly_order];
    vlasov->accel_boundary_surf[1] = accel_boundary_surf_vy_kernels[cv_index].kernels[poly_order];
    vlasov->accel_boundary_surf[2] = accel_boundary_surf_vz_kernels[cv_index].kernels[poly_order];    
  }
  else {
    if (field_id == GKYL_FIELD_NULL) {
      vlasov->eqn.vol_term = stream_vol_kernels[cv_index].kernels[poly_order];
    }
    else if (field_id == GKYL_FIELD_PHI) {
      vlasov->eqn.vol_term = poisson_vol_kernels[cv_index].kernels[poly_order];

      vlasov->accel_surf[0] = poisson_accel_surf_vx_kernels[cv_index].kernels[poly_order];
      if (vdim>1)
        vlasov->accel_surf[1] = poisson_accel_surf_vy_kernels[cv_index].kernels[poly_order];
      if (vdim>2)
        vlasov->accel_surf[2] = poisson_accel_surf_vz_kernels[cv_index].kernels[poly_order];

      vlasov->accel_boundary_surf[0] = poisson_accel_boundary_surf_vx_kernels[cv_index].kernels[poly_order];
      if (vdim>1)
        vlasov->accel_boundary_surf[1] = poisson_accel_boundary_surf_vy_kernels[cv_index].kernels[poly_order];
      if (vdim>2)
        vlasov->accel_boundary_surf[2] = poisson_accel_boundary_surf_vz_kernels[cv_index].kernels[poly_order];
    }
    else if (field_id == GKYL_FIELD_PHI_A) {
      vlasov->eqn.vol_term = poisson_extem_vol_kernels[cv_index].kernels[poly_order];

      vlasov->accel_surf[0] = poisson_extem_accel_surf_vx_kernels[cv_index].kernels[poly_order];
      if (vdim>1)
        vlasov->accel_surf[1] = poisson_extem_accel_surf_vy_kernels[cv_index].kernels[poly_order];
      if (vdim>2)
        vlasov->accel_surf[2] = poisson_extem_accel_surf_vz_kernels[cv_index].kernels[poly_order];

      vlasov->accel_boundary_surf[0] = poisson_extem_accel_boundary_surf_vx_kernels[cv_index].kernels[poly_order];
      if (vdim>1)
        vlasov->accel_boundary_surf[1] = poisson_extem_accel_boundary_surf_vy_kernels[cv_index].kernels[poly_order];
      if (vdim>2)
        vlasov->accel_boundary_surf[2] = poisson_extem_accel_boundary_surf_vz_kernels[cv_index].kernels[poly_order];
    }
    else {
      vlasov->eqn.vol_term = vol_kernels[cv_index].kernels[poly_order];

      vlasov->accel_surf[0] = accel_surf_vx_kernels[cv_index].kernels[poly_order];
      if (vdim>1)
        vlasov->accel_surf[1] = accel_surf_vy_kernels[cv_index].kernels[poly_order];
      if (vdim>2)
        vlasov->accel_surf[2] = accel_surf_vz_kernels[cv_index].kernels[poly_order];

      vlasov->accel_boundary_surf[0] = accel_boundary_surf_vx_kernels[cv_index].kernels[poly_order];
      if (vdim>1)
        vlasov->accel_boundary_surf[1] = accel_boundary_surf_vy_kernels[cv_index].kernels[poly_order];
      if (vdim>2)
        vlasov->accel_boundary_surf[2] = accel_boundary_surf_vz_kernels[cv_index].kernels[poly_order];
    }
    // Streaming kernels are the same for each field_id
    vlasov->stream_surf[0] = stream_surf_x_kernels[cv_index].kernels[poly_order];
    if (cdim>1)
      vlasov->stream_surf[1] = stream_surf_y_kernels[cv_index].kernels[poly_order];
    if (cdim>2)
      vlasov->stream_surf[2] = stream_surf_z_kernels[cv_index].kernels[poly_order];

    vlasov->stream_boundary_surf[0] = stream_boundary_surf_x_kernels[cv_index].kernels[poly_order];
    if (cdim>1)
      vlasov->stream_boundary_surf[1] = stream_boundary_surf_y_kernels[cv_index].kernels[poly_order];
    if (cdim>2)
      vlasov->stream_boundary_surf[2] = stream_boundary_surf_z_kernels[cv_index].kernels[poly_order];
  }    
}

struct gkyl_dg_eqn*
gkyl_dg_vlasov_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  const struct gkyl_range* conf_range, const struct gkyl_range* phase_range,
  enum gkyl_model_id model_id, enum gkyl_field_id field_id)
{
  struct dg_vlasov *vlasov = (struct dg_vlasov*) gkyl_malloc(sizeof(struct dg_vlasov));

  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  vlasov->cdim = cdim;
  vlasov->pdim = pdim;

  vlasov->eqn.num_equations = 1;
  vlasov->conf_range = *conf_range;
  vlasov->phase_range = *phase_range;

  vlasov->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(vlasov->eqn.flags);
  vlasov->eqn.ref_count = gkyl_ref_count_init(gkyl_vlasov_free);

  // copy the host struct to device struct
  struct dg_vlasov *vlasov_cu = (struct dg_vlasov*) gkyl_cu_malloc(sizeof(struct dg_vlasov));
  gkyl_cu_memcpy(vlasov_cu, vlasov, sizeof(struct dg_vlasov), GKYL_CU_MEMCPY_H2D);

  dg_vlasov_set_cu_dev_ptrs<<<1,1>>>(vlasov_cu, cbasis->b_type, cv_index[cdim].vdim[vdim],
    cdim, vdim, poly_order, model_id, field_id);

  // set parent on_dev pointer
  vlasov->eqn.on_dev = &vlasov_cu->eqn;
  
  return &vlasov->eqn;
}
