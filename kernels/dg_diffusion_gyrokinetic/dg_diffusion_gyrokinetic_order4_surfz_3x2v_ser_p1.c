#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_surfz_3x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[2],4.);

  out[0] += -0.0625*(25.98076211353316*coeff[2]*qr[3]-25.98076211353316*coeff[2]*ql[3]+((-15.0*qr[0])-15.0*ql[0]+30.0*qc[0])*coeff[2])*Jfac; 
  out[1] += -0.0625*(25.98076211353316*coeff[2]*qr[7]-25.98076211353316*coeff[2]*ql[7]+((-15.0*qr[1])-15.0*ql[1]+30.0*qc[1])*coeff[2])*Jfac; 
  out[2] += -0.0625*(25.98076211353316*coeff[2]*qr[8]-25.98076211353316*coeff[2]*ql[8]-15.0*coeff[2]*qr[2]-15.0*coeff[2]*ql[2]+30.0*coeff[2]*qc[2])*Jfac; 
  out[3] += -0.0625*(33.0*coeff[2]*qr[3]+33.0*coeff[2]*ql[3]+114.0*coeff[2]*qc[3]+(25.98076211353316*ql[0]-25.98076211353316*qr[0])*coeff[2])*Jfac; 
  out[4] += -0.0625*(25.98076211353316*coeff[2]*qr[11]-25.98076211353316*coeff[2]*ql[11]-15.0*coeff[2]*qr[4]-15.0*coeff[2]*ql[4]+30.0*coeff[2]*qc[4])*Jfac; 
  out[5] += -0.0625*(25.98076211353316*coeff[2]*qr[14]-25.98076211353316*coeff[2]*ql[14]-15.0*coeff[2]*qr[5]-15.0*coeff[2]*ql[5]+30.0*coeff[2]*qc[5])*Jfac; 
  out[6] += -0.0625*(25.98076211353316*coeff[2]*qr[16]-25.98076211353316*coeff[2]*ql[16]-15.0*coeff[2]*qr[6]-15.0*coeff[2]*ql[6]+30.0*coeff[2]*qc[6])*Jfac; 
  out[7] += -0.0625*(33.0*coeff[2]*qr[7]+33.0*coeff[2]*ql[7]+114.0*coeff[2]*qc[7]+(25.98076211353316*ql[1]-25.98076211353316*qr[1])*coeff[2])*Jfac; 
  out[8] += -0.0625*(33.0*coeff[2]*qr[8]+33.0*coeff[2]*ql[8]+114.0*coeff[2]*qc[8]-25.98076211353316*coeff[2]*qr[2]+25.98076211353316*coeff[2]*ql[2])*Jfac; 
  out[9] += -0.0625*(25.98076211353316*coeff[2]*qr[18]-25.98076211353316*coeff[2]*ql[18]-15.0*coeff[2]*qr[9]-15.0*coeff[2]*ql[9]+30.0*coeff[2]*qc[9])*Jfac; 
  out[10] += -0.0625*(25.98076211353316*coeff[2]*qr[19]-25.98076211353316*coeff[2]*ql[19]-15.0*coeff[2]*qr[10]-15.0*coeff[2]*ql[10]+30.0*coeff[2]*qc[10])*Jfac; 
  out[11] += -0.0625*(33.0*coeff[2]*qr[11]+33.0*coeff[2]*ql[11]+114.0*coeff[2]*qc[11]-25.98076211353316*coeff[2]*qr[4]+25.98076211353316*coeff[2]*ql[4])*Jfac; 
  out[12] += -0.0625*(25.98076211353316*coeff[2]*qr[21]-25.98076211353316*coeff[2]*ql[21]-15.0*coeff[2]*qr[12]-15.0*coeff[2]*ql[12]+30.0*coeff[2]*qc[12])*Jfac; 
  out[13] += -0.0625*(25.98076211353316*coeff[2]*qr[22]-25.98076211353316*coeff[2]*ql[22]-15.0*coeff[2]*qr[13]-15.0*coeff[2]*ql[13]+30.0*coeff[2]*qc[13])*Jfac; 
  out[14] += -0.0625*(33.0*coeff[2]*qr[14]+33.0*coeff[2]*ql[14]+114.0*coeff[2]*qc[14]-25.98076211353316*coeff[2]*qr[5]+25.98076211353316*coeff[2]*ql[5])*Jfac; 
  out[15] += -0.0625*(25.98076211353316*coeff[2]*qr[25]-25.98076211353316*coeff[2]*ql[25]-15.0*coeff[2]*qr[15]-15.0*coeff[2]*ql[15]+30.0*coeff[2]*qc[15])*Jfac; 
  out[16] += -0.0625*(33.0*coeff[2]*qr[16]+33.0*coeff[2]*ql[16]+114.0*coeff[2]*qc[16]-25.98076211353316*coeff[2]*qr[6]+25.98076211353316*coeff[2]*ql[6])*Jfac; 
  out[17] += -0.0625*(25.98076211353316*coeff[2]*qr[26]-25.98076211353316*coeff[2]*ql[26]-15.0*coeff[2]*qr[17]-15.0*coeff[2]*ql[17]+30.0*coeff[2]*qc[17])*Jfac; 
  out[18] += -0.0625*(33.0*coeff[2]*qr[18]+33.0*coeff[2]*ql[18]+114.0*coeff[2]*qc[18]-25.98076211353316*coeff[2]*qr[9]+25.98076211353316*coeff[2]*ql[9])*Jfac; 
  out[19] += -0.0625*(33.0*coeff[2]*qr[19]+33.0*coeff[2]*ql[19]+114.0*coeff[2]*qc[19]-25.98076211353316*coeff[2]*qr[10]+25.98076211353316*coeff[2]*ql[10])*Jfac; 
  out[20] += -0.0625*(25.98076211353316*coeff[2]*qr[27]-25.98076211353316*coeff[2]*ql[27]-15.0*coeff[2]*qr[20]-15.0*coeff[2]*ql[20]+30.0*coeff[2]*qc[20])*Jfac; 
  out[21] += -0.0625*(33.0*coeff[2]*qr[21]+33.0*coeff[2]*ql[21]+114.0*coeff[2]*qc[21]-25.98076211353316*coeff[2]*qr[12]+25.98076211353316*coeff[2]*ql[12])*Jfac; 
  out[22] += -0.0625*(33.0*coeff[2]*qr[22]+33.0*coeff[2]*ql[22]+114.0*coeff[2]*qc[22]-25.98076211353316*coeff[2]*qr[13]+25.98076211353316*coeff[2]*ql[13])*Jfac; 
  out[23] += -0.0625*(25.98076211353316*coeff[2]*qr[29]-25.98076211353316*coeff[2]*ql[29]-15.0*coeff[2]*qr[23]-15.0*coeff[2]*ql[23]+30.0*coeff[2]*qc[23])*Jfac; 
  out[24] += -0.0625*(25.98076211353316*coeff[2]*qr[30]-25.98076211353316*coeff[2]*ql[30]-15.0*coeff[2]*qr[24]-15.0*coeff[2]*ql[24]+30.0*coeff[2]*qc[24])*Jfac; 
  out[25] += -0.0625*(33.0*coeff[2]*qr[25]+33.0*coeff[2]*ql[25]+114.0*coeff[2]*qc[25]-25.98076211353316*coeff[2]*qr[15]+25.98076211353316*coeff[2]*ql[15])*Jfac; 
  out[26] += -0.0625*(33.0*coeff[2]*qr[26]+33.0*coeff[2]*ql[26]+114.0*coeff[2]*qc[26]-25.98076211353316*coeff[2]*qr[17]+25.98076211353316*coeff[2]*ql[17])*Jfac; 
  out[27] += -0.0625*(33.0*coeff[2]*qr[27]+33.0*coeff[2]*ql[27]+114.0*coeff[2]*qc[27]-25.98076211353316*coeff[2]*qr[20]+25.98076211353316*coeff[2]*ql[20])*Jfac; 
  out[28] += -0.0625*(25.98076211353316*coeff[2]*qr[31]-25.98076211353316*coeff[2]*ql[31]-15.0*coeff[2]*qr[28]-15.0*coeff[2]*ql[28]+30.0*coeff[2]*qc[28])*Jfac; 
  out[29] += -0.0625*(33.0*coeff[2]*qr[29]+33.0*coeff[2]*ql[29]+114.0*coeff[2]*qc[29]-25.98076211353316*coeff[2]*qr[23]+25.98076211353316*coeff[2]*ql[23])*Jfac; 
  out[30] += -0.0625*(33.0*coeff[2]*qr[30]+33.0*coeff[2]*ql[30]+114.0*coeff[2]*qc[30]-25.98076211353316*coeff[2]*qr[24]+25.98076211353316*coeff[2]*ql[24])*Jfac; 
  out[31] += -0.0625*(33.0*coeff[2]*qr[31]+33.0*coeff[2]*ql[31]+114.0*coeff[2]*qc[31]-25.98076211353316*coeff[2]*qr[28]+25.98076211353316*coeff[2]*ql[28])*Jfac; 
  out[32] += -0.0625*(25.98076211353316*coeff[2]*qr[35]-25.98076211353316*coeff[2]*ql[35]-15.0*coeff[2]*qr[32]-15.0*coeff[2]*ql[32]+30.0*coeff[2]*qc[32])*Jfac; 
  out[33] += -0.0625*(25.98076211353316*coeff[2]*qr[38]-25.98076211353316*coeff[2]*ql[38]-15.0*coeff[2]*qr[33]-15.0*coeff[2]*ql[33]+30.0*coeff[2]*qc[33])*Jfac; 
  out[34] += -0.0625*(25.98076211353316*coeff[2]*qr[39]-25.98076211353316*coeff[2]*ql[39]-15.0*coeff[2]*qr[34]-15.0*coeff[2]*ql[34]+30.0*coeff[2]*qc[34])*Jfac; 
  out[35] += -0.0625*(33.0*coeff[2]*qr[35]+33.0*coeff[2]*ql[35]+114.0*coeff[2]*qc[35]-25.98076211353316*coeff[2]*qr[32]+25.98076211353316*coeff[2]*ql[32])*Jfac; 
  out[36] += -0.0625*(25.98076211353316*coeff[2]*qr[42]-25.98076211353316*coeff[2]*ql[42]-15.0*coeff[2]*qr[36]-15.0*coeff[2]*ql[36]+30.0*coeff[2]*qc[36])*Jfac; 
  out[37] += -0.0625*(25.98076211353316*coeff[2]*qr[43]-25.98076211353316*coeff[2]*ql[43]-15.0*coeff[2]*qr[37]-15.0*coeff[2]*ql[37]+30.0*coeff[2]*qc[37])*Jfac; 
  out[38] += -0.0625*(33.0*coeff[2]*qr[38]+33.0*coeff[2]*ql[38]+114.0*coeff[2]*qc[38]-25.98076211353316*coeff[2]*qr[33]+25.98076211353316*coeff[2]*ql[33])*Jfac; 
  out[39] += -0.0625*(33.0*coeff[2]*qr[39]+33.0*coeff[2]*ql[39]+114.0*coeff[2]*qc[39]-25.98076211353316*coeff[2]*qr[34]+25.98076211353316*coeff[2]*ql[34])*Jfac; 
  out[40] += -0.0625*(25.98076211353316*coeff[2]*qr[45]-25.98076211353316*coeff[2]*ql[45]-15.0*coeff[2]*qr[40]-15.0*coeff[2]*ql[40]+30.0*coeff[2]*qc[40])*Jfac; 
  out[41] += -0.0625*(25.98076211353316*coeff[2]*qr[46]-25.98076211353316*coeff[2]*ql[46]-15.0*coeff[2]*qr[41]-15.0*coeff[2]*ql[41]+30.0*coeff[2]*qc[41])*Jfac; 
  out[42] += -0.0625*(33.0*coeff[2]*qr[42]+33.0*coeff[2]*ql[42]+114.0*coeff[2]*qc[42]-25.98076211353316*coeff[2]*qr[36]+25.98076211353316*coeff[2]*ql[36])*Jfac; 
  out[43] += -0.0625*(33.0*coeff[2]*qr[43]+33.0*coeff[2]*ql[43]+114.0*coeff[2]*qc[43]-25.98076211353316*coeff[2]*qr[37]+25.98076211353316*coeff[2]*ql[37])*Jfac; 
  out[44] += -0.0625*(25.98076211353316*coeff[2]*qr[47]-25.98076211353316*coeff[2]*ql[47]-15.0*coeff[2]*qr[44]-15.0*coeff[2]*ql[44]+30.0*coeff[2]*qc[44])*Jfac; 
  out[45] += -0.0625*(33.0*coeff[2]*qr[45]+33.0*coeff[2]*ql[45]+114.0*coeff[2]*qc[45]-25.98076211353316*coeff[2]*qr[40]+25.98076211353316*coeff[2]*ql[40])*Jfac; 
  out[46] += -0.0625*(33.0*coeff[2]*qr[46]+33.0*coeff[2]*ql[46]+114.0*coeff[2]*qc[46]-25.98076211353316*coeff[2]*qr[41]+25.98076211353316*coeff[2]*ql[41])*Jfac; 
  out[47] += -0.0625*(33.0*coeff[2]*qr[47]+33.0*coeff[2]*ql[47]+114.0*coeff[2]*qc[47]-25.98076211353316*coeff[2]*qr[44]+25.98076211353316*coeff[2]*ql[44])*Jfac; 

  return 0.;

}

