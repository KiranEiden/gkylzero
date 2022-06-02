#include <gkyl_dg_diffusion_kernels.h>

GKYL_CU_DH void
dg_diffusion_surfx_2x_ser_p1(const double* w, const double* dx,
  const double* D,
  const double* ql, const double* qc, const double* qr,
  double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dxv[NDIM]: Cell spacing
  // D: Diffusion coefficient in the center cell
  // ql: Input field in the left cell
  // qc: Input field in the center cell
  // qr: Input field in the right cell
  // out: Incremented output

  const double J = 4/dx[0]/dx[0];

  out[0] += J*((-0.46875*D[3]*qr[3])-0.270632938682637*D[2]*qr[3]-0.46875*D[3]*ql[3]+0.270632938682637*D[2]*ql[3]-0.9375*D[3]*qc[3]+0.4871392896287466*qr[2]*D[3]-0.4871392896287466*ql[2]*D[3]+0.28125*D[2]*qr[2]+0.28125*D[2]*ql[2]-0.5625*D[2]*qc[2]-0.46875*D[1]*qr[1]-0.270632938682637*D[0]*qr[1]-0.46875*D[1]*ql[1]+0.270632938682637*D[0]*ql[1]-0.9375*D[1]*qc[1]+0.4871392896287466*qr[0]*D[1]-0.4871392896287466*ql[0]*D[1]+0.28125*D[0]*qr[0]+0.28125*D[0]*ql[0]-0.5625*D[0]*qc[0]);
  out[1] += J*((-0.3788861141556919*D[3]*qr[3])-0.21875*D[2]*qr[3]+0.3788861141556919*D[3]*ql[3]-0.21875*D[2]*ql[3]-1.4375*D[2]*qc[3]+0.46875*qr[2]*D[3]+0.46875*ql[2]*D[3]-0.9375*qc[2]*D[3]+0.270632938682637*D[2]*qr[2]-0.270632938682637*D[2]*ql[2]-0.3788861141556919*D[1]*qr[1]-0.21875*D[0]*qr[1]+0.3788861141556919*D[1]*ql[1]-0.21875*D[0]*ql[1]-1.4375*D[0]*qc[1]+0.46875*qr[0]*D[1]+0.46875*ql[0]*D[1]-0.9375*qc[0]*D[1]+0.270632938682637*D[0]*qr[0]-0.270632938682637*D[0]*ql[0]);
  out[2] += J*((-0.46875*D[1]*qr[3])-0.270632938682637*D[0]*qr[3]-0.46875*D[1]*ql[3]+0.270632938682637*D[0]*ql[3]-0.9375*D[1]*qc[3]-0.46875*qr[1]*D[3]-0.46875*ql[1]*D[3]-0.9375*qc[1]*D[3]+0.4871392896287466*qr[0]*D[3]-0.4871392896287466*ql[0]*D[3]+0.4871392896287466*D[1]*qr[2]+0.28125*D[0]*qr[2]-0.4871392896287466*D[1]*ql[2]+0.28125*D[0]*ql[2]-0.5625*D[0]*qc[2]-0.270632938682637*qr[1]*D[2]+0.270632938682637*ql[1]*D[2]+0.28125*qr[0]*D[2]+0.28125*ql[0]*D[2]-0.5625*qc[0]*D[2]);
  out[3] += J*((-0.3788861141556919*D[1]*qr[3])-0.21875*D[0]*qr[3]+0.3788861141556919*D[1]*ql[3]-0.21875*D[0]*ql[3]-1.4375*D[0]*qc[3]-0.3788861141556919*qr[1]*D[3]+0.3788861141556919*ql[1]*D[3]+0.46875*qr[0]*D[3]+0.46875*ql[0]*D[3]-0.9375*qc[0]*D[3]+0.46875*D[1]*qr[2]+0.270632938682637*D[0]*qr[2]+0.46875*D[1]*ql[2]-0.270632938682637*D[0]*ql[2]-0.9375*D[1]*qc[2]-0.21875*qr[1]*D[2]-0.21875*ql[1]*D[2]-1.4375*qc[1]*D[2]+0.270632938682637*qr[0]*D[2]-0.270632938682637*ql[0]*D[2]);
}