#include <gkyl_dg_gen_diffusion_kernels.h>

GKYL_CU_DH void
dg_gen_diffusion_surfx_2x_ser_p1(const double* w, const double* dx,
  const double* Dij, const double* q[], double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dxv[NDIM]: Cell spacing
  // Dij: Diffusion coefficient in the center cell
  // q: Input field in the left cell
  // out: Incremented output

  const double Jxx = 4/dx[0]/dx[0];
  const double Jxy = 4/dx[0]/dx[1];

  const double* Dxx = &Dij[0];
  const double* Dxy = &Dij[4];
  const double* qll = q[0];
  const double* qlc = q[1];
  const double* qlu = q[2];
  const double* qcl = q[3];
  const double* qcc = q[4];
  const double* qcu = q[5];
  const double* qul = q[6];
  const double* quc = q[7];
  const double* quu = q[8];

  out[0] += Jxy*(0.125*Dxy[3]*quu[3]+0.07216878364870323*Dxy[2]*quu[3]+0.07216878364870323*Dxy[1]*quu[3]+0.04166666666666666*Dxy[0]*quu[3]-0.125*Dxy[3]*qul[3]-0.07216878364870323*Dxy[2]*qul[3]+0.07216878364870323*Dxy[1]*qul[3]+0.04166666666666666*Dxy[0]*qul[3]-0.1443375672974065*Dxy[1]*quc[3]-0.08333333333333333*Dxy[0]*quc[3]-0.125*Dxy[3]*qlu[3]+0.07216878364870323*Dxy[2]*qlu[3]-0.07216878364870323*Dxy[1]*qlu[3]+0.04166666666666666*Dxy[0]*qlu[3]+0.125*Dxy[3]*qll[3]-0.07216878364870323*Dxy[2]*qll[3]-0.07216878364870323*Dxy[1]*qll[3]+0.04166666666666666*Dxy[0]*qll[3]+0.1443375672974065*Dxy[1]*qlc[3]-0.08333333333333333*Dxy[0]*qlc[3]-0.1443375672974065*Dxy[2]*qcu[3]-0.08333333333333333*Dxy[0]*qcu[3]+0.1443375672974065*Dxy[2]*qcl[3]-0.08333333333333333*Dxy[0]*qcl[3]+0.1666666666666667*Dxy[0]*qcc[3]-0.1082531754730548*quu[2]*Dxy[3]+0.1082531754730548*qul[2]*Dxy[3]-0.1082531754730548*qlu[2]*Dxy[3]+0.1082531754730548*qll[2]*Dxy[3]-0.2165063509461096*qcu[2]*Dxy[3]+0.2165063509461096*qcl[2]*Dxy[3]-0.1082531754730548*quu[1]*Dxy[3]-0.1082531754730548*qul[1]*Dxy[3]+0.2165063509461096*quc[1]*Dxy[3]+0.1082531754730548*qlu[1]*Dxy[3]+0.1082531754730548*qll[1]*Dxy[3]-0.2165063509461096*qlc[1]*Dxy[3]+0.09375*quu[0]*Dxy[3]+0.09375*qul[0]*Dxy[3]-0.1875*quc[0]*Dxy[3]+0.09375*qlu[0]*Dxy[3]+0.09375*qll[0]*Dxy[3]-0.1875*qlc[0]*Dxy[3]+0.1875*qcu[0]*Dxy[3]+0.1875*qcl[0]*Dxy[3]-0.375*qcc[0]*Dxy[3]-0.0625*Dxy[2]*quu[2]-0.0625*Dxy[1]*quu[2]-0.03608439182435162*Dxy[0]*quu[2]+0.0625*Dxy[2]*qul[2]-0.0625*Dxy[1]*qul[2]-0.03608439182435162*Dxy[0]*qul[2]+0.125*Dxy[1]*quc[2]+0.07216878364870323*Dxy[0]*quc[2]+0.0625*Dxy[2]*qlu[2]-0.0625*Dxy[1]*qlu[2]+0.03608439182435162*Dxy[0]*qlu[2]-0.0625*Dxy[2]*qll[2]-0.0625*Dxy[1]*qll[2]+0.03608439182435162*Dxy[0]*qll[2]+0.125*Dxy[1]*qlc[2]-0.07216878364870323*Dxy[0]*qlc[2]-0.125*Dxy[1]*qcu[2]-0.125*Dxy[1]*qcl[2]+0.25*Dxy[1]*qcc[2]-0.0625*quu[1]*Dxy[2]-0.0625*qul[1]*Dxy[2]+0.125*quc[1]*Dxy[2]-0.0625*qlu[1]*Dxy[2]-0.0625*qll[1]*Dxy[2]+0.125*qlc[1]*Dxy[2]+0.125*qcu[1]*Dxy[2]+0.125*qcl[1]*Dxy[2]-0.25*qcc[1]*Dxy[2]+0.05412658773652741*quu[0]*Dxy[2]+0.05412658773652741*qul[0]*Dxy[2]-0.1082531754730548*quc[0]*Dxy[2]-0.05412658773652741*qlu[0]*Dxy[2]-0.05412658773652741*qll[0]*Dxy[2]+0.1082531754730548*qlc[0]*Dxy[2]-0.0625*Dxy[1]*quu[1]-0.03608439182435162*Dxy[0]*quu[1]+0.0625*Dxy[1]*qul[1]+0.03608439182435162*Dxy[0]*qul[1]+0.0625*Dxy[1]*qlu[1]-0.03608439182435162*Dxy[0]*qlu[1]-0.0625*Dxy[1]*qll[1]+0.03608439182435162*Dxy[0]*qll[1]+0.07216878364870323*Dxy[0]*qcu[1]-0.07216878364870323*Dxy[0]*qcl[1]+0.05412658773652741*quu[0]*Dxy[1]-0.05412658773652741*qul[0]*Dxy[1]+0.05412658773652741*qlu[0]*Dxy[1]-0.05412658773652741*qll[0]*Dxy[1]+0.1082531754730548*qcu[0]*Dxy[1]-0.1082531754730548*qcl[0]*Dxy[1]+0.03125*Dxy[0]*quu[0]-0.03125*Dxy[0]*qul[0]-0.03125*Dxy[0]*qlu[0]+0.03125*Dxy[0]*qll[0]) + Jxx*((-0.46875*Dxx[3]*quc[3])-0.270632938682637*Dxx[2]*quc[3]-0.46875*Dxx[3]*qlc[3]+0.270632938682637*Dxx[2]*qlc[3]-0.9375*Dxx[3]*qcc[3]+0.4871392896287466*quc[2]*Dxx[3]-0.4871392896287466*qlc[2]*Dxx[3]+0.28125*Dxx[2]*quc[2]+0.28125*Dxx[2]*qlc[2]-0.5625*Dxx[2]*qcc[2]-0.46875*Dxx[1]*quc[1]-0.270632938682637*Dxx[0]*quc[1]-0.46875*Dxx[1]*qlc[1]+0.270632938682637*Dxx[0]*qlc[1]-0.9375*Dxx[1]*qcc[1]+0.4871392896287466*quc[0]*Dxx[1]-0.4871392896287466*qlc[0]*Dxx[1]+0.28125*Dxx[0]*quc[0]+0.28125*Dxx[0]*qlc[0]-0.5625*Dxx[0]*qcc[0]);
  out[1] += Jxy*(0.2165063509461096*Dxy[3]*quu[3]+0.125*Dxy[2]*quu[3]+0.125*Dxy[1]*quu[3]+0.07216878364870323*Dxy[0]*quu[3]-0.2165063509461096*Dxy[3]*qul[3]-0.125*Dxy[2]*qul[3]+0.125*Dxy[1]*qul[3]+0.07216878364870323*Dxy[0]*qul[3]-0.25*Dxy[1]*quc[3]-0.1443375672974065*Dxy[0]*quc[3]+0.2165063509461096*Dxy[3]*qlu[3]-0.125*Dxy[2]*qlu[3]+0.125*Dxy[1]*qlu[3]-0.07216878364870323*Dxy[0]*qlu[3]-0.2165063509461096*Dxy[3]*qll[3]+0.125*Dxy[2]*qll[3]+0.125*Dxy[1]*qll[3]-0.07216878364870323*Dxy[0]*qll[3]-0.25*Dxy[1]*qlc[3]+0.1443375672974065*Dxy[0]*qlc[3]-0.1875*quu[2]*Dxy[3]+0.1875*qul[2]*Dxy[3]+0.1875*qlu[2]*Dxy[3]-0.1875*qll[2]*Dxy[3]-0.1875*quu[1]*Dxy[3]-0.1875*qul[1]*Dxy[3]+0.375*quc[1]*Dxy[3]-0.1875*qlu[1]*Dxy[3]-0.1875*qll[1]*Dxy[3]+0.375*qlc[1]*Dxy[3]+0.1623797632095822*quu[0]*Dxy[3]+0.1623797632095822*qul[0]*Dxy[3]-0.3247595264191644*quc[0]*Dxy[3]-0.1623797632095822*qlu[0]*Dxy[3]-0.1623797632095822*qll[0]*Dxy[3]+0.3247595264191644*qlc[0]*Dxy[3]-0.1082531754730548*Dxy[2]*quu[2]-0.1082531754730548*Dxy[1]*quu[2]-0.0625*Dxy[0]*quu[2]+0.1082531754730548*Dxy[2]*qul[2]-0.1082531754730548*Dxy[1]*qul[2]-0.0625*Dxy[0]*qul[2]+0.2165063509461096*Dxy[1]*quc[2]+0.125*Dxy[0]*quc[2]-0.1082531754730548*Dxy[2]*qlu[2]+0.1082531754730548*Dxy[1]*qlu[2]-0.0625*Dxy[0]*qlu[2]+0.1082531754730548*Dxy[2]*qll[2]+0.1082531754730548*Dxy[1]*qll[2]-0.0625*Dxy[0]*qll[2]-0.2165063509461096*Dxy[1]*qlc[2]+0.125*Dxy[0]*qlc[2]+0.2165063509461096*Dxy[2]*qcu[2]+0.125*Dxy[0]*qcu[2]-0.2165063509461096*Dxy[2]*qcl[2]+0.125*Dxy[0]*qcl[2]-0.25*Dxy[0]*qcc[2]-0.1082531754730548*quu[1]*Dxy[2]-0.1082531754730548*qul[1]*Dxy[2]+0.2165063509461096*quc[1]*Dxy[2]+0.1082531754730548*qlu[1]*Dxy[2]+0.1082531754730548*qll[1]*Dxy[2]-0.2165063509461096*qlc[1]*Dxy[2]+0.09375*quu[0]*Dxy[2]+0.09375*qul[0]*Dxy[2]-0.1875*quc[0]*Dxy[2]+0.09375*qlu[0]*Dxy[2]+0.09375*qll[0]*Dxy[2]-0.1875*qlc[0]*Dxy[2]-0.1875*qcu[0]*Dxy[2]-0.1875*qcl[0]*Dxy[2]+0.375*qcc[0]*Dxy[2]-0.1082531754730548*Dxy[1]*quu[1]-0.0625*Dxy[0]*quu[1]+0.1082531754730548*Dxy[1]*qul[1]+0.0625*Dxy[0]*qul[1]-0.1082531754730548*Dxy[1]*qlu[1]+0.0625*Dxy[0]*qlu[1]+0.1082531754730548*Dxy[1]*qll[1]-0.0625*Dxy[0]*qll[1]+0.09375*quu[0]*Dxy[1]-0.09375*qul[0]*Dxy[1]-0.09375*qlu[0]*Dxy[1]+0.09375*qll[0]*Dxy[1]+0.05412658773652741*Dxy[0]*quu[0]-0.05412658773652741*Dxy[0]*qul[0]+0.05412658773652741*Dxy[0]*qlu[0]-0.05412658773652741*Dxy[0]*qll[0]-0.1082531754730548*Dxy[0]*qcu[0]+0.1082531754730548*Dxy[0]*qcl[0]) + Jxx*((-0.3788861141556919*Dxx[3]*quc[3])-0.21875*Dxx[2]*quc[3]+0.3788861141556919*Dxx[3]*qlc[3]-0.21875*Dxx[2]*qlc[3]-1.4375*Dxx[2]*qcc[3]+0.46875*quc[2]*Dxx[3]+0.46875*qlc[2]*Dxx[3]-0.9375*qcc[2]*Dxx[3]+0.270632938682637*Dxx[2]*quc[2]-0.270632938682637*Dxx[2]*qlc[2]-0.3788861141556919*Dxx[1]*quc[1]-0.21875*Dxx[0]*quc[1]+0.3788861141556919*Dxx[1]*qlc[1]-0.21875*Dxx[0]*qlc[1]-1.4375*Dxx[0]*qcc[1]+0.46875*quc[0]*Dxx[1]+0.46875*qlc[0]*Dxx[1]-0.9375*qcc[0]*Dxx[1]+0.270632938682637*Dxx[0]*quc[0]-0.270632938682637*Dxx[0]*qlc[0]);
  out[2] += Jxy*(0.2165063509461096*Dxy[3]*quu[3]+0.125*Dxy[2]*quu[3]+0.125*Dxy[1]*quu[3]+0.07216878364870323*Dxy[0]*quu[3]+0.2165063509461096*Dxy[3]*qul[3]+0.125*Dxy[2]*qul[3]-0.125*Dxy[1]*qul[3]-0.07216878364870323*Dxy[0]*qul[3]+0.4330127018922193*Dxy[3]*quc[3]+0.25*Dxy[2]*quc[3]-0.2165063509461096*Dxy[3]*qlu[3]+0.125*Dxy[2]*qlu[3]-0.125*Dxy[1]*qlu[3]+0.07216878364870323*Dxy[0]*qlu[3]-0.2165063509461096*Dxy[3]*qll[3]+0.125*Dxy[2]*qll[3]+0.125*Dxy[1]*qll[3]-0.07216878364870323*Dxy[0]*qll[3]-0.4330127018922193*Dxy[3]*qlc[3]+0.25*Dxy[2]*qlc[3]-0.25*Dxy[2]*qcu[3]-0.1443375672974065*Dxy[0]*qcu[3]-0.25*Dxy[2]*qcl[3]+0.1443375672974065*Dxy[0]*qcl[3]-0.5*Dxy[2]*qcc[3]-0.1875*quu[2]*Dxy[3]-0.1875*qul[2]*Dxy[3]-0.375*quc[2]*Dxy[3]-0.1875*qlu[2]*Dxy[3]-0.1875*qll[2]*Dxy[3]-0.375*qlc[2]*Dxy[3]-0.375*qcu[2]*Dxy[3]-0.375*qcl[2]*Dxy[3]-0.75*qcc[2]*Dxy[3]-0.1875*quu[1]*Dxy[3]+0.1875*qul[1]*Dxy[3]+0.1875*qlu[1]*Dxy[3]-0.1875*qll[1]*Dxy[3]+0.1623797632095822*quu[0]*Dxy[3]-0.1623797632095822*qul[0]*Dxy[3]+0.1623797632095822*qlu[0]*Dxy[3]-0.1623797632095822*qll[0]*Dxy[3]+0.3247595264191644*qcu[0]*Dxy[3]-0.3247595264191644*qcl[0]*Dxy[3]-0.1082531754730548*Dxy[2]*quu[2]-0.1082531754730548*Dxy[1]*quu[2]-0.0625*Dxy[0]*quu[2]-0.1082531754730548*Dxy[2]*qul[2]+0.1082531754730548*Dxy[1]*qul[2]+0.0625*Dxy[0]*qul[2]-0.2165063509461096*Dxy[2]*quc[2]+0.1082531754730548*Dxy[2]*qlu[2]-0.1082531754730548*Dxy[1]*qlu[2]+0.0625*Dxy[0]*qlu[2]+0.1082531754730548*Dxy[2]*qll[2]+0.1082531754730548*Dxy[1]*qll[2]-0.0625*Dxy[0]*qll[2]+0.2165063509461096*Dxy[2]*qlc[2]-0.2165063509461096*Dxy[1]*qcu[2]+0.2165063509461096*Dxy[1]*qcl[2]-0.1082531754730548*quu[1]*Dxy[2]+0.1082531754730548*qul[1]*Dxy[2]-0.1082531754730548*qlu[1]*Dxy[2]+0.1082531754730548*qll[1]*Dxy[2]+0.2165063509461096*qcu[1]*Dxy[2]-0.2165063509461096*qcl[1]*Dxy[2]+0.09375*quu[0]*Dxy[2]-0.09375*qul[0]*Dxy[2]-0.09375*qlu[0]*Dxy[2]+0.09375*qll[0]*Dxy[2]-0.1082531754730548*Dxy[1]*quu[1]-0.0625*Dxy[0]*quu[1]-0.1082531754730548*Dxy[1]*qul[1]-0.0625*Dxy[0]*qul[1]+0.2165063509461096*Dxy[1]*quc[1]+0.125*Dxy[0]*quc[1]+0.1082531754730548*Dxy[1]*qlu[1]-0.0625*Dxy[0]*qlu[1]+0.1082531754730548*Dxy[1]*qll[1]-0.0625*Dxy[0]*qll[1]-0.2165063509461096*Dxy[1]*qlc[1]+0.125*Dxy[0]*qlc[1]+0.125*Dxy[0]*qcu[1]+0.125*Dxy[0]*qcl[1]-0.25*Dxy[0]*qcc[1]+0.09375*quu[0]*Dxy[1]+0.09375*qul[0]*Dxy[1]-0.1875*quc[0]*Dxy[1]+0.09375*qlu[0]*Dxy[1]+0.09375*qll[0]*Dxy[1]-0.1875*qlc[0]*Dxy[1]+0.1875*qcu[0]*Dxy[1]+0.1875*qcl[0]*Dxy[1]-0.375*qcc[0]*Dxy[1]+0.05412658773652741*Dxy[0]*quu[0]+0.05412658773652741*Dxy[0]*qul[0]-0.1082531754730548*Dxy[0]*quc[0]-0.05412658773652741*Dxy[0]*qlu[0]-0.05412658773652741*Dxy[0]*qll[0]+0.1082531754730548*Dxy[0]*qlc[0]) + Jxx*((-0.46875*Dxx[1]*quc[3])-0.270632938682637*Dxx[0]*quc[3]-0.46875*Dxx[1]*qlc[3]+0.270632938682637*Dxx[0]*qlc[3]-0.9375*Dxx[1]*qcc[3]-0.46875*quc[1]*Dxx[3]-0.46875*qlc[1]*Dxx[3]-0.9375*qcc[1]*Dxx[3]+0.4871392896287466*quc[0]*Dxx[3]-0.4871392896287466*qlc[0]*Dxx[3]+0.4871392896287466*Dxx[1]*quc[2]+0.28125*Dxx[0]*quc[2]-0.4871392896287466*Dxx[1]*qlc[2]+0.28125*Dxx[0]*qlc[2]-0.5625*Dxx[0]*qcc[2]-0.270632938682637*quc[1]*Dxx[2]+0.270632938682637*qlc[1]*Dxx[2]+0.28125*quc[0]*Dxx[2]+0.28125*qlc[0]*Dxx[2]-0.5625*qcc[0]*Dxx[2]);
  out[3] += Jxy*(0.375*Dxy[3]*quu[3]+0.2165063509461096*Dxy[2]*quu[3]+0.2165063509461096*Dxy[1]*quu[3]+0.125*Dxy[0]*quu[3]+0.375*Dxy[3]*qul[3]+0.2165063509461096*Dxy[2]*qul[3]-0.2165063509461096*Dxy[1]*qul[3]-0.125*Dxy[0]*qul[3]+0.75*Dxy[3]*quc[3]+0.4330127018922193*Dxy[2]*quc[3]+0.375*Dxy[3]*qlu[3]-0.2165063509461096*Dxy[2]*qlu[3]+0.2165063509461096*Dxy[1]*qlu[3]-0.125*Dxy[0]*qlu[3]+0.375*Dxy[3]*qll[3]-0.2165063509461096*Dxy[2]*qll[3]-0.2165063509461096*Dxy[1]*qll[3]+0.125*Dxy[0]*qll[3]+0.75*Dxy[3]*qlc[3]-0.4330127018922193*Dxy[2]*qlc[3]-0.3247595264191644*quu[2]*Dxy[3]-0.3247595264191644*qul[2]*Dxy[3]-0.6495190528383289*quc[2]*Dxy[3]+0.3247595264191644*qlu[2]*Dxy[3]+0.3247595264191644*qll[2]*Dxy[3]+0.6495190528383289*qlc[2]*Dxy[3]-0.3247595264191644*quu[1]*Dxy[3]+0.3247595264191644*qul[1]*Dxy[3]-0.3247595264191644*qlu[1]*Dxy[3]+0.3247595264191644*qll[1]*Dxy[3]+0.28125*quu[0]*Dxy[3]-0.28125*qul[0]*Dxy[3]-0.28125*qlu[0]*Dxy[3]+0.28125*qll[0]*Dxy[3]-0.1875*Dxy[2]*quu[2]-0.1875*Dxy[1]*quu[2]-0.1082531754730548*Dxy[0]*quu[2]-0.1875*Dxy[2]*qul[2]+0.1875*Dxy[1]*qul[2]+0.1082531754730548*Dxy[0]*qul[2]-0.375*Dxy[2]*quc[2]-0.1875*Dxy[2]*qlu[2]+0.1875*Dxy[1]*qlu[2]-0.1082531754730548*Dxy[0]*qlu[2]-0.1875*Dxy[2]*qll[2]-0.1875*Dxy[1]*qll[2]+0.1082531754730548*Dxy[0]*qll[2]-0.375*Dxy[2]*qlc[2]+0.375*Dxy[2]*qcu[2]+0.2165063509461096*Dxy[0]*qcu[2]+0.375*Dxy[2]*qcl[2]-0.2165063509461096*Dxy[0]*qcl[2]+0.75*Dxy[2]*qcc[2]-0.1875*quu[1]*Dxy[2]+0.1875*qul[1]*Dxy[2]+0.1875*qlu[1]*Dxy[2]-0.1875*qll[1]*Dxy[2]+0.1623797632095822*quu[0]*Dxy[2]-0.1623797632095822*qul[0]*Dxy[2]+0.1623797632095822*qlu[0]*Dxy[2]-0.1623797632095822*qll[0]*Dxy[2]-0.3247595264191644*qcu[0]*Dxy[2]+0.3247595264191644*qcl[0]*Dxy[2]-0.1875*Dxy[1]*quu[1]-0.1082531754730548*Dxy[0]*quu[1]-0.1875*Dxy[1]*qul[1]-0.1082531754730548*Dxy[0]*qul[1]+0.375*Dxy[1]*quc[1]+0.2165063509461096*Dxy[0]*quc[1]-0.1875*Dxy[1]*qlu[1]+0.1082531754730548*Dxy[0]*qlu[1]-0.1875*Dxy[1]*qll[1]+0.1082531754730548*Dxy[0]*qll[1]+0.375*Dxy[1]*qlc[1]-0.2165063509461096*Dxy[0]*qlc[1]+0.1623797632095822*quu[0]*Dxy[1]+0.1623797632095822*qul[0]*Dxy[1]-0.3247595264191644*quc[0]*Dxy[1]-0.1623797632095822*qlu[0]*Dxy[1]-0.1623797632095822*qll[0]*Dxy[1]+0.3247595264191644*qlc[0]*Dxy[1]+0.09375*Dxy[0]*quu[0]+0.09375*Dxy[0]*qul[0]-0.1875*Dxy[0]*quc[0]+0.09375*Dxy[0]*qlu[0]+0.09375*Dxy[0]*qll[0]-0.1875*Dxy[0]*qlc[0]-0.1875*Dxy[0]*qcu[0]-0.1875*Dxy[0]*qcl[0]+0.375*Dxy[0]*qcc[0]) + Jxx*((-0.3788861141556919*Dxx[1]*quc[3])-0.21875*Dxx[0]*quc[3]+0.3788861141556919*Dxx[1]*qlc[3]-0.21875*Dxx[0]*qlc[3]-1.4375*Dxx[0]*qcc[3]-0.3788861141556919*quc[1]*Dxx[3]+0.3788861141556919*qlc[1]*Dxx[3]+0.46875*quc[0]*Dxx[3]+0.46875*qlc[0]*Dxx[3]-0.9375*qcc[0]*Dxx[3]+0.46875*Dxx[1]*quc[2]+0.270632938682637*Dxx[0]*quc[2]+0.46875*Dxx[1]*qlc[2]-0.270632938682637*Dxx[0]*qlc[2]-0.9375*Dxx[1]*qcc[2]-0.21875*quc[1]*Dxx[2]-0.21875*qlc[1]*Dxx[2]-1.4375*qcc[1]*Dxx[2]+0.270632938682637*quc[0]*Dxx[2]-0.270632938682637*qlc[0]*Dxx[2]);
}
