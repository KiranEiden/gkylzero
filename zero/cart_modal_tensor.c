#include <assert.h>
#include <math.h>
#include <string.h>

#include <gkyl_basis.h>
#include <gkyl_util.h>
#include <gkyl_basis_ser_kernels.h>
#include <gkyl_basis_tensor_kernels.h>

// Basis function eval for each dimension: ev_list[ndim].ev[poly_order]
static struct { void (*ev[4])(const double *z, double *b); } ev_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { eval_1d_ser_p0, eval_1d_ser_p1, eval_1d_ser_p2, eval_1d_ser_p3 },
  { eval_2d_ser_p0, eval_2d_ser_p1, eval_2d_tensor_p2, NULL },
  { eval_3d_ser_p0, eval_3d_ser_p1, eval_3d_tensor_p2, NULL },
  { eval_4d_ser_p0, eval_4d_ser_p1, eval_4d_tensor_p2, NULL },
  { eval_5d_ser_p0, eval_5d_ser_p1, eval_5d_tensor_p2, NULL },
  { eval_6d_ser_p0, eval_6d_ser_p1, NULL, NULL },
};

// Flip-sign functions: ev_list[ndim].ev[poly_order]
static struct { void (*fs[4])(int dir, const double *f, double *fout); } fs_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { flip_sign_1d_ser_p0, flip_sign_1d_ser_p1, flip_sign_1d_ser_p2, flip_sign_1d_ser_p3 },
  { flip_sign_2d_ser_p0, flip_sign_2d_ser_p1, flip_sign_2d_tensor_p2, NULL },
  { flip_sign_3d_ser_p0, flip_sign_3d_ser_p1, flip_sign_3d_tensor_p2, NULL },
  { flip_sign_4d_ser_p0, flip_sign_4d_ser_p1, flip_sign_4d_tensor_p2, NULL },
  { flip_sign_5d_ser_p0, flip_sign_5d_ser_p1, flip_sign_5d_tensor_p2, NULL },
  { flip_sign_6d_ser_p0, flip_sign_6d_ser_p1, NULL, NULL },
};

void
gkyl_cart_modal_tensor(struct gkyl_basis *basis, int ndim, int poly_order)
{
  assert(ndim>0 && ndim<=6);
  assert(ev_list[ndim].ev[poly_order]);
  
  basis->ndim = ndim;
  basis->poly_order = poly_order;
  basis->num_basis = pow(poly_order+1, ndim);
  strcpy(basis->id, "tensor");
  basis->b_type = GKYL_BASIS_MODAL_TENSOR;
  basis->eval = ev_list[ndim].ev[poly_order];
  basis->flip_sign = fs_list[ndim].fs[poly_order];
}