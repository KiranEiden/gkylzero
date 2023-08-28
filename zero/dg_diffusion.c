#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_diffusion.h>
#include <gkyl_dg_diffusion_priv.h>
#include <gkyl_util.h>

void
gkyl_diffusion_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_eqn *base = container_of(ref, struct gkyl_dg_eqn, ref_count);

  if (gkyl_dg_eqn_is_cu_dev(base)) {
    // free inner on_dev object
    struct dg_diffusion *diffusion = container_of(base->on_dev, struct dg_diffusion, eqn);
    gkyl_cu_free(diffusion);
  }
  
  struct dg_diffusion *diffusion = container_of(base, struct dg_diffusion, eqn);
  gkyl_free(diffusion);
}

void
gkyl_diffusion_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_diffusion_auxfields auxin)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(auxin.D)) {
    gkyl_diffusion_set_auxfields_cu(eqn->on_dev, auxin);
    return;
  }
#endif
  
  struct dg_diffusion *diffusion = container_of(eqn, struct dg_diffusion, eqn);
  diffusion->auxfields.D = auxin.D;
}

struct gkyl_dg_eqn*
gkyl_dg_diffusion_new(const struct gkyl_basis *basis,
  const struct gkyl_basis *cbasis, enum gkyl_diffusion_id diffusion_id, bool *diff_in_dir,
  int diff_order, const struct gkyl_range *conf_range, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    return gkyl_dg_diffusion_cu_dev_new(diffusion_id, basis, cbasis, diff_in_dir, diff_order, conf_range);
#endif
  
  struct dg_diffusion *diffusion = gkyl_malloc(sizeof(struct dg_diffusion));

  int cdim = cbasis->ndim;
  int vdim = basis->ndim - cdim;
  int poly_order = cbasis->poly_order;

  int num_equations = 1;
  if ((diffusion_id == GKYL_DIFFUSION_DIAGONAL_CONST_EULER) || (diffusion_id == GKYL_DIFFUSION_DIAGONAL_VAR_EULER))
    num_equations = 4;
  else if ((diffusion_id == GKYL_DIFFUSION_DIAGONAL_CONST_EULER_ISO) || (diffusion_id == GKYL_DIFFUSION_DIAGONAL_VAR_EULER_ISO))
    num_equations = 5;
  if ((diffusion_id == GKYL_DIFFUSION_DIAGONAL_CONST_PKPM) || (diffusion_id == GKYL_DIFFUSION_DIAGONAL_VAR_PKPM))
    num_equations = 3;

  diffusion->const_coeff = (   (diffusion_id == GKYL_DIFFUSION_DIAGONAL_CONST)
                            || (diffusion_id == GKYL_DIFFUSION_DIAGONAL_CONST_EULER)
                            || (diffusion_id == GKYL_DIFFUSION_DIAGONAL_CONST_EULER_ISO)
                            || (diffusion_id == GKYL_DIFFUSION_DIAGONAL_CONST_VLASOV)
                            || (diffusion_id == GKYL_DIFFUSION_DIAGONAL_CONST_GYROKINETIC)
                            || (diffusion_id == GKYL_DIFFUSION_DIAGONAL_CONST_PKPM) );

  const gkyl_dg_diffusion_vol_kern_list *vol_kernels;
  const gkyl_dg_diffusion_surf_kern_list *surfx_kernels;
  const gkyl_dg_diffusion_surf_kern_list *surfy_kernels;
  const gkyl_dg_diffusion_surf_kern_list *surfz_kernels; 
  const gkyl_dg_diffusion_boundary_surf_kern_list *boundary_surfx_kernels;
  const gkyl_dg_diffusion_boundary_surf_kern_list *boundary_surfy_kernels;
  const gkyl_dg_diffusion_boundary_surf_kern_list *boundary_surfz_kernels; 

  if ((diffusion_id == GKYL_DIFFUSION_DIAGONAL_CONST_VLASOV) || (diffusion_id == GKYL_DIFFUSION_DIAGONAL_VAR_VLASOV)) {
    switch (cbasis->b_type) {
      case GKYL_BASIS_MODAL_SERENDIPITY:
        vol_kernels            = diffusion->const_coeff? ser_vol_kernels_constcoeff                   : ser_vol_kernels_varcoeff                  ;
        surfx_kernels          = diffusion->const_coeff? ser_vlasov_surfx_kernels_constcoeff          : ser_vlasov_surfx_kernels_varcoeff         ;
        surfy_kernels          = diffusion->const_coeff? ser_vlasov_surfy_kernels_constcoeff          : ser_vlasov_surfy_kernels_varcoeff         ;
        surfz_kernels          = diffusion->const_coeff? ser_vlasov_surfz_kernels_constcoeff          : ser_vlasov_surfz_kernels_varcoeff         ;
        boundary_surfx_kernels = diffusion->const_coeff? ser_vlasov_boundary_surfx_kernels_constcoeff : ser_vlasov_boundary_surfx_kernels_varcoeff;
        boundary_surfy_kernels = diffusion->const_coeff? ser_vlasov_boundary_surfy_kernels_constcoeff : ser_vlasov_boundary_surfy_kernels_varcoeff;
        boundary_surfz_kernels = diffusion->const_coeff? ser_vlasov_boundary_surfz_kernels_constcoeff : ser_vlasov_boundary_surfz_kernels_varcoeff;
        break;
  
      default:
        assert(false);
        break;    
    } 
  } else {
    switch (cbasis->b_type) {
      case GKYL_BASIS_MODAL_SERENDIPITY:
        vol_kernels            = diffusion->const_coeff? ser_vol_kernels_constcoeff            : ser_vol_kernels_varcoeff           ;
        surfx_kernels          = diffusion->const_coeff? ser_surfx_kernels_constcoeff          : ser_surfx_kernels_varcoeff         ;
        surfy_kernels          = diffusion->const_coeff? ser_surfy_kernels_constcoeff          : ser_surfy_kernels_varcoeff         ;
        surfz_kernels          = diffusion->const_coeff? ser_surfz_kernels_constcoeff          : ser_surfz_kernels_varcoeff         ;
        boundary_surfx_kernels = diffusion->const_coeff? ser_boundary_surfx_kernels_constcoeff : ser_boundary_surfx_kernels_varcoeff;
        boundary_surfy_kernels = diffusion->const_coeff? ser_boundary_surfy_kernels_constcoeff : ser_boundary_surfy_kernels_varcoeff;
        boundary_surfz_kernels = diffusion->const_coeff? ser_boundary_surfz_kernels_constcoeff : ser_boundary_surfz_kernels_varcoeff;
        break;
  
      default:
        assert(false);
        break;    
    } 
  } 

  diffusion->num_equations = num_equations;
  diffusion->num_basis = basis->num_basis;
  if (diff_in_dir)
    for (size_t d=0; d<cdim; d++) diffusion->diff_in_dir[d] = diff_in_dir[d];
  else
    for (size_t d=0; d<cdim; d++) diffusion->diff_in_dir[d] = true;

  // Linear index into list of volume kernels.
  int dirs_bin_key[] = {1,2,4,8,16,32}; // Binary: 000001, 000010, 000100, 001000, 010000, 100000.
  int dirs_linidx = 0; // Binary 000000.
  for (int d=0; d<cdim; d++) {
     if (diffusion->diff_in_dir[d]) dirs_linidx = dirs_linidx | dirs_bin_key[d];
  }
  dirs_linidx -= 1;

  diffusion->eqn.num_equations = num_equations;
  diffusion->eqn.surf_term = surf;
  diffusion->eqn.boundary_surf_term = boundary_surf;

  diffusion->eqn.vol_term = CKVOL(vol_kernels, cdim, diff_order, poly_order, dirs_linidx);

  diffusion->surf[0] = CKSURF(surfx_kernels, diff_order, cdim, vdim, poly_order, diffusion_id);
  if (cdim>1)
    diffusion->surf[1] = CKSURF(surfy_kernels, diff_order, cdim, vdim, poly_order, diffusion_id);
  if (cdim>2)
    diffusion->surf[2] = CKSURF(surfz_kernels, diff_order, cdim, vdim, poly_order, diffusion_id);

  diffusion->boundary_surf[0] = CKSURF(boundary_surfx_kernels, diff_order, cdim, vdim, poly_order, diffusion_id);
  if (cdim>1)
    diffusion->boundary_surf[1] = CKSURF(boundary_surfy_kernels, diff_order, cdim, vdim, poly_order, diffusion_id);
  if (cdim>2)
    diffusion->boundary_surf[2] = CKSURF(boundary_surfz_kernels, diff_order, cdim, vdim, poly_order, diffusion_id);

  // Ensure non-NULL pointers.
  for (int i=0; i<cdim; ++i) assert(diffusion->surf[i]);

  diffusion->auxfields.D = 0;
  diffusion->conf_range = *conf_range;

  diffusion->eqn.flags = 0;
  diffusion->eqn.ref_count = gkyl_ref_count_init(gkyl_diffusion_free);
  diffusion->eqn.on_dev = &diffusion->eqn;
  
  return &diffusion->eqn;
}
