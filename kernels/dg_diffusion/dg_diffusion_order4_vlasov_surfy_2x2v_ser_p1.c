#include <gkyl_dg_diffusion_kernels.h>

GKYL_CU_DH double dg_diffusion_order4_vlasov_surfy_2x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[1],4.);

  out[0] += -0.0625*(25.98076211353316*coeff[1]*qr[2]-25.98076211353316*coeff[1]*ql[2]+((-15.0*qr[0])-15.0*ql[0]+30.0*qc[0])*coeff[1])*Jfac; 
  out[1] += -0.0625*(25.98076211353316*coeff[1]*qr[5]-25.98076211353316*coeff[1]*ql[5]-15.0*coeff[1]*qr[1]-15.0*coeff[1]*ql[1]+30.0*coeff[1]*qc[1])*Jfac; 
  out[2] += -0.0625*(33.0*coeff[1]*qr[2]+33.0*coeff[1]*ql[2]+114.0*coeff[1]*qc[2]+(25.98076211353316*ql[0]-25.98076211353316*qr[0])*coeff[1])*Jfac; 
  out[3] += -0.0625*(25.98076211353316*coeff[1]*qr[7]-25.98076211353316*coeff[1]*ql[7]-15.0*coeff[1]*qr[3]-15.0*coeff[1]*ql[3]+30.0*coeff[1]*qc[3])*Jfac; 
  out[4] += -0.0625*(25.98076211353316*coeff[1]*qr[9]-25.98076211353316*coeff[1]*ql[9]-15.0*coeff[1]*qr[4]-15.0*coeff[1]*ql[4]+30.0*coeff[1]*qc[4])*Jfac; 
  out[5] += -0.0625*(33.0*coeff[1]*qr[5]+33.0*coeff[1]*ql[5]+114.0*coeff[1]*qc[5]-25.98076211353316*coeff[1]*qr[1]+25.98076211353316*coeff[1]*ql[1])*Jfac; 
  out[6] += -0.0625*(25.98076211353316*coeff[1]*qr[11]-25.98076211353316*coeff[1]*ql[11]-15.0*coeff[1]*qr[6]-15.0*coeff[1]*ql[6]+30.0*coeff[1]*qc[6])*Jfac; 
  out[7] += -0.0625*(33.0*coeff[1]*qr[7]+33.0*coeff[1]*ql[7]+114.0*coeff[1]*qc[7]-25.98076211353316*coeff[1]*qr[3]+25.98076211353316*coeff[1]*ql[3])*Jfac; 
  out[8] += -0.0625*(25.98076211353316*coeff[1]*qr[12]-25.98076211353316*coeff[1]*ql[12]-15.0*coeff[1]*qr[8]-15.0*coeff[1]*ql[8]+30.0*coeff[1]*qc[8])*Jfac; 
  out[9] += -0.0625*(33.0*coeff[1]*qr[9]+33.0*coeff[1]*ql[9]+114.0*coeff[1]*qc[9]-25.98076211353316*coeff[1]*qr[4]+25.98076211353316*coeff[1]*ql[4])*Jfac; 
  out[10] += -0.0625*(25.98076211353316*coeff[1]*qr[14]-25.98076211353316*coeff[1]*ql[14]-15.0*coeff[1]*qr[10]-15.0*coeff[1]*ql[10]+30.0*coeff[1]*qc[10])*Jfac; 
  out[11] += -0.0625*(33.0*coeff[1]*qr[11]+33.0*coeff[1]*ql[11]+114.0*coeff[1]*qc[11]-25.98076211353316*coeff[1]*qr[6]+25.98076211353316*coeff[1]*ql[6])*Jfac; 
  out[12] += -0.0625*(33.0*coeff[1]*qr[12]+33.0*coeff[1]*ql[12]+114.0*coeff[1]*qc[12]-25.98076211353316*coeff[1]*qr[8]+25.98076211353316*coeff[1]*ql[8])*Jfac; 
  out[13] += -0.0625*(25.98076211353316*coeff[1]*qr[15]-25.98076211353316*coeff[1]*ql[15]-15.0*coeff[1]*qr[13]-15.0*coeff[1]*ql[13]+30.0*coeff[1]*qc[13])*Jfac; 
  out[14] += -0.0625*(33.0*coeff[1]*qr[14]+33.0*coeff[1]*ql[14]+114.0*coeff[1]*qc[14]-25.98076211353316*coeff[1]*qr[10]+25.98076211353316*coeff[1]*ql[10])*Jfac; 
  out[15] += -0.0625*(33.0*coeff[1]*qr[15]+33.0*coeff[1]*ql[15]+114.0*coeff[1]*qc[15]-25.98076211353316*coeff[1]*qr[13]+25.98076211353316*coeff[1]*ql[13])*Jfac; 
  out[16] += -0.0625*(25.98076211353316*coeff[1]*qr[18]-25.98076211353316*coeff[1]*ql[18]-15.0*coeff[1]*qr[16]-15.0*coeff[1]*ql[16]+30.0*coeff[1]*qc[16])*Jfac; 
  out[17] += -0.0625*(25.98076211353316*coeff[1]*qr[20]-25.98076211353316*coeff[1]*ql[20]-15.0*coeff[1]*qr[17]-15.0*coeff[1]*ql[17]+30.0*coeff[1]*qc[17])*Jfac; 
  out[18] += -0.0625*(33.0*coeff[1]*qr[18]+33.0*coeff[1]*ql[18]+114.0*coeff[1]*qc[18]-25.98076211353316*coeff[1]*qr[16]+25.98076211353316*coeff[1]*ql[16])*Jfac; 
  out[19] += -0.0625*(25.98076211353316*coeff[1]*qr[22]-25.98076211353316*coeff[1]*ql[22]-15.0*coeff[1]*qr[19]-15.0*coeff[1]*ql[19]+30.0*coeff[1]*qc[19])*Jfac; 
  out[20] += -0.0625*(33.0*coeff[1]*qr[20]+33.0*coeff[1]*ql[20]+114.0*coeff[1]*qc[20]-25.98076211353316*coeff[1]*qr[17]+25.98076211353316*coeff[1]*ql[17])*Jfac; 
  out[21] += -0.0625*(25.98076211353316*coeff[1]*qr[23]-25.98076211353316*coeff[1]*ql[23]-15.0*coeff[1]*qr[21]-15.0*coeff[1]*ql[21]+30.0*coeff[1]*qc[21])*Jfac; 
  out[22] += -0.0625*(33.0*coeff[1]*qr[22]+33.0*coeff[1]*ql[22]+114.0*coeff[1]*qc[22]-25.98076211353316*coeff[1]*qr[19]+25.98076211353316*coeff[1]*ql[19])*Jfac; 
  out[23] += -0.0625*(33.0*coeff[1]*qr[23]+33.0*coeff[1]*ql[23]+114.0*coeff[1]*qc[23]-25.98076211353316*coeff[1]*qr[21]+25.98076211353316*coeff[1]*ql[21])*Jfac; 
  out[24] += -0.0625*(25.98076211353316*coeff[1]*qr[26]-25.98076211353316*coeff[1]*ql[26]-15.0*coeff[1]*qr[24]-15.0*coeff[1]*ql[24]+30.0*coeff[1]*qc[24])*Jfac; 
  out[25] += -0.0625*(25.98076211353316*coeff[1]*qr[28]-25.98076211353316*coeff[1]*ql[28]-15.0*coeff[1]*qr[25]-15.0*coeff[1]*ql[25]+30.0*coeff[1]*qc[25])*Jfac; 
  out[26] += -0.0625*(33.0*coeff[1]*qr[26]+33.0*coeff[1]*ql[26]+114.0*coeff[1]*qc[26]-25.98076211353316*coeff[1]*qr[24]+25.98076211353316*coeff[1]*ql[24])*Jfac; 
  out[27] += -0.0625*(25.98076211353316*coeff[1]*qr[30]-25.98076211353316*coeff[1]*ql[30]-15.0*coeff[1]*qr[27]-15.0*coeff[1]*ql[27]+30.0*coeff[1]*qc[27])*Jfac; 
  out[28] += -0.0625*(33.0*coeff[1]*qr[28]+33.0*coeff[1]*ql[28]+114.0*coeff[1]*qc[28]-25.98076211353316*coeff[1]*qr[25]+25.98076211353316*coeff[1]*ql[25])*Jfac; 
  out[29] += -0.0625*(25.98076211353316*coeff[1]*qr[31]-25.98076211353316*coeff[1]*ql[31]-15.0*coeff[1]*qr[29]-15.0*coeff[1]*ql[29]+30.0*coeff[1]*qc[29])*Jfac; 
  out[30] += -0.0625*(33.0*coeff[1]*qr[30]+33.0*coeff[1]*ql[30]+114.0*coeff[1]*qc[30]-25.98076211353316*coeff[1]*qr[27]+25.98076211353316*coeff[1]*ql[27])*Jfac; 
  out[31] += -0.0625*(33.0*coeff[1]*qr[31]+33.0*coeff[1]*ql[31]+114.0*coeff[1]*qc[31]-25.98076211353316*coeff[1]*qr[29]+25.98076211353316*coeff[1]*ql[29])*Jfac; 

  return 0.;

}

GKYL_CU_DH double dg_diffusion_order4_vlasov_surfy_2x2v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[1],4.);

  out[0] += -0.03125*((57.0*qr[5]+57.0*ql[5]+66.0*qc[5]-25.98076211353316*qr[1]+25.98076211353316*ql[1])*coeff[7]+(57.0*qr[2]+57.0*ql[2]+66.0*qc[2]-25.98076211353316*qr[0]+25.98076211353316*ql[0])*coeff[6]+25.98076211353316*coeff[5]*qr[5]-25.98076211353316*coeff[5]*ql[5]+((-15.0*qr[1])-15.0*ql[1]+30.0*qc[1])*coeff[5]+(25.98076211353316*qr[2]-25.98076211353316*ql[2]-15.0*qr[0]-15.0*ql[0]+30.0*qc[0])*coeff[4])*Jfac; 
  out[1] += -0.03125*((57.0*qr[2]+57.0*ql[2]+66.0*qc[2]-25.98076211353316*qr[0]+25.98076211353316*ql[0])*coeff[7]+(57.0*qr[5]+57.0*ql[5]+66.0*qc[5]-25.98076211353316*qr[1]+25.98076211353316*ql[1])*coeff[6]+25.98076211353316*coeff[4]*qr[5]-25.98076211353316*coeff[4]*ql[5]+(25.98076211353316*qr[2]-25.98076211353316*ql[2]-15.0*qr[0]-15.0*ql[0]+30.0*qc[0])*coeff[5]+((-15.0*qr[1])-15.0*ql[1]+30.0*qc[1])*coeff[4])*Jfac; 
  out[2] += -0.03125*((77.94228634059945*qr[5]-77.94228634059945*ql[5]-45.0*qr[1]-45.0*ql[1]+90.0*qc[1])*coeff[7]+(77.94228634059945*qr[2]-77.94228634059945*ql[2]-45.0*qr[0]-45.0*ql[0]+90.0*qc[0])*coeff[6]+33.0*coeff[5]*qr[5]+33.0*coeff[5]*ql[5]+114.0*coeff[5]*qc[5]+(25.98076211353316*ql[1]-25.98076211353316*qr[1])*coeff[5]+(33.0*qr[2]+33.0*ql[2]+114.0*qc[2]-25.98076211353316*qr[0]+25.98076211353316*ql[0])*coeff[4])*Jfac; 
  out[3] += -0.03125*((57.0*coeff[7]+25.98076211353316*coeff[5])*qr[11]+(57.0*coeff[7]-25.98076211353316*coeff[5])*ql[11]+66.0*coeff[7]*qc[11]+(57.0*coeff[6]+25.98076211353316*coeff[4])*qr[7]+(57.0*coeff[6]-25.98076211353316*coeff[4])*ql[7]+66.0*coeff[6]*qc[7]+(25.98076211353316*ql[6]-25.98076211353316*qr[6])*coeff[7]-15.0*coeff[5]*qr[6]-15.0*coeff[5]*ql[6]+30.0*coeff[5]*qc[6]+(25.98076211353316*ql[3]-25.98076211353316*qr[3])*coeff[6]+((-15.0*qr[3])-15.0*ql[3]+30.0*qc[3])*coeff[4])*Jfac; 
  out[4] += -0.03125*((57.0*coeff[7]+25.98076211353316*coeff[5])*qr[12]+(57.0*coeff[7]-25.98076211353316*coeff[5])*ql[12]+66.0*coeff[7]*qc[12]+(57.0*coeff[6]+25.98076211353316*coeff[4])*qr[9]+(57.0*coeff[6]-25.98076211353316*coeff[4])*ql[9]+66.0*coeff[6]*qc[9]+((-25.98076211353316*coeff[7])-15.0*coeff[5])*qr[8]+(25.98076211353316*coeff[7]-15.0*coeff[5])*ql[8]+30.0*coeff[5]*qc[8]+(25.98076211353316*ql[4]-25.98076211353316*qr[4])*coeff[6]-15.0*coeff[4]*qr[4]-15.0*coeff[4]*ql[4]+30.0*coeff[4]*qc[4])*Jfac; 
  out[5] += -0.03125*((77.94228634059945*qr[2]-77.94228634059945*ql[2]-45.0*qr[0]-45.0*ql[0]+90.0*qc[0])*coeff[7]+(77.94228634059945*qr[5]-77.94228634059945*ql[5]-45.0*qr[1]-45.0*ql[1]+90.0*qc[1])*coeff[6]+33.0*coeff[4]*qr[5]+33.0*coeff[4]*ql[5]+114.0*coeff[4]*qc[5]+(33.0*qr[2]+33.0*ql[2]+114.0*qc[2]-25.98076211353316*qr[0]+25.98076211353316*ql[0])*coeff[5]+(25.98076211353316*ql[1]-25.98076211353316*qr[1])*coeff[4])*Jfac; 
  out[6] += -0.03125*((57.0*coeff[6]+25.98076211353316*coeff[4])*qr[11]+(57.0*coeff[6]-25.98076211353316*coeff[4])*ql[11]+66.0*coeff[6]*qc[11]+(57.0*coeff[7]+25.98076211353316*coeff[5])*qr[7]+(57.0*coeff[7]-25.98076211353316*coeff[5])*ql[7]+66.0*coeff[7]*qc[7]+(25.98076211353316*ql[3]-25.98076211353316*qr[3])*coeff[7]+((-25.98076211353316*coeff[6])-15.0*coeff[4])*qr[6]+(25.98076211353316*coeff[6]-15.0*coeff[4])*ql[6]+30.0*coeff[4]*qc[6]+((-15.0*qr[3])-15.0*ql[3]+30.0*qc[3])*coeff[5])*Jfac; 
  out[7] += -0.03125*((77.94228634059945*coeff[7]+33.0*coeff[5])*qr[11]+(33.0*coeff[5]-77.94228634059945*coeff[7])*ql[11]+114.0*coeff[5]*qc[11]+(77.94228634059945*coeff[6]+33.0*coeff[4])*qr[7]+(33.0*coeff[4]-77.94228634059945*coeff[6])*ql[7]+114.0*coeff[4]*qc[7]+((-45.0*qr[6])-45.0*ql[6]+90.0*qc[6])*coeff[7]-25.98076211353316*coeff[5]*qr[6]+25.98076211353316*coeff[5]*ql[6]+((-45.0*qr[3])-45.0*ql[3]+90.0*qc[3])*coeff[6]+(25.98076211353316*ql[3]-25.98076211353316*qr[3])*coeff[4])*Jfac; 
  out[8] += -0.03125*((57.0*coeff[6]+25.98076211353316*coeff[4])*qr[12]+(57.0*coeff[6]-25.98076211353316*coeff[4])*ql[12]+66.0*coeff[6]*qc[12]+(57.0*coeff[7]+25.98076211353316*coeff[5])*qr[9]+(57.0*coeff[7]-25.98076211353316*coeff[5])*ql[9]+66.0*coeff[7]*qc[9]+((-25.98076211353316*coeff[6])-15.0*coeff[4])*qr[8]+(25.98076211353316*coeff[6]-15.0*coeff[4])*ql[8]+30.0*coeff[4]*qc[8]+(25.98076211353316*ql[4]-25.98076211353316*qr[4])*coeff[7]+((-15.0*qr[4])-15.0*ql[4]+30.0*qc[4])*coeff[5])*Jfac; 
  out[9] += -0.03125*((77.94228634059945*coeff[7]+33.0*coeff[5])*qr[12]+(33.0*coeff[5]-77.94228634059945*coeff[7])*ql[12]+114.0*coeff[5]*qc[12]+(77.94228634059945*coeff[6]+33.0*coeff[4])*qr[9]+(33.0*coeff[4]-77.94228634059945*coeff[6])*ql[9]+114.0*coeff[4]*qc[9]+((-45.0*coeff[7])-25.98076211353316*coeff[5])*qr[8]+(25.98076211353316*coeff[5]-45.0*coeff[7])*ql[8]+90.0*coeff[7]*qc[8]+((-45.0*qr[4])-45.0*ql[4]+90.0*qc[4])*coeff[6]-25.98076211353316*coeff[4]*qr[4]+25.98076211353316*coeff[4]*ql[4])*Jfac; 
  out[10] += -0.03125*((57.0*coeff[7]+25.98076211353316*coeff[5])*qr[15]+(57.0*coeff[7]-25.98076211353316*coeff[5])*ql[15]+66.0*coeff[7]*qc[15]+(57.0*coeff[6]+25.98076211353316*coeff[4])*qr[14]+(57.0*coeff[6]-25.98076211353316*coeff[4])*ql[14]+66.0*coeff[6]*qc[14]+((-25.98076211353316*coeff[7])-15.0*coeff[5])*qr[13]+(25.98076211353316*coeff[7]-15.0*coeff[5])*ql[13]+30.0*coeff[5]*qc[13]+((-25.98076211353316*coeff[6])-15.0*coeff[4])*qr[10]+(25.98076211353316*coeff[6]-15.0*coeff[4])*ql[10]+30.0*coeff[4]*qc[10])*Jfac; 
  out[11] += -0.03125*((77.94228634059945*coeff[6]+33.0*coeff[4])*qr[11]+(33.0*coeff[4]-77.94228634059945*coeff[6])*ql[11]+114.0*coeff[4]*qc[11]+(77.94228634059945*coeff[7]+33.0*coeff[5])*qr[7]+(33.0*coeff[5]-77.94228634059945*coeff[7])*ql[7]+114.0*coeff[5]*qc[7]+((-45.0*qr[3])-45.0*ql[3]+90.0*qc[3])*coeff[7]+((-45.0*coeff[6])-25.98076211353316*coeff[4])*qr[6]+(25.98076211353316*coeff[4]-45.0*coeff[6])*ql[6]+90.0*coeff[6]*qc[6]+(25.98076211353316*ql[3]-25.98076211353316*qr[3])*coeff[5])*Jfac; 
  out[12] += -0.03125*((77.94228634059945*coeff[6]+33.0*coeff[4])*qr[12]+(33.0*coeff[4]-77.94228634059945*coeff[6])*ql[12]+114.0*coeff[4]*qc[12]+(77.94228634059945*coeff[7]+33.0*coeff[5])*qr[9]+(33.0*coeff[5]-77.94228634059945*coeff[7])*ql[9]+114.0*coeff[5]*qc[9]+((-45.0*coeff[6])-25.98076211353316*coeff[4])*qr[8]+(25.98076211353316*coeff[4]-45.0*coeff[6])*ql[8]+90.0*coeff[6]*qc[8]+((-45.0*qr[4])-45.0*ql[4]+90.0*qc[4])*coeff[7]+(25.98076211353316*ql[4]-25.98076211353316*qr[4])*coeff[5])*Jfac; 
  out[13] += -0.03125*((57.0*coeff[6]+25.98076211353316*coeff[4])*qr[15]+(57.0*coeff[6]-25.98076211353316*coeff[4])*ql[15]+66.0*coeff[6]*qc[15]+(57.0*coeff[7]+25.98076211353316*coeff[5])*qr[14]+(57.0*coeff[7]-25.98076211353316*coeff[5])*ql[14]+66.0*coeff[7]*qc[14]+((-25.98076211353316*coeff[6])-15.0*coeff[4])*qr[13]+(25.98076211353316*coeff[6]-15.0*coeff[4])*ql[13]+30.0*coeff[4]*qc[13]+((-25.98076211353316*coeff[7])-15.0*coeff[5])*qr[10]+(25.98076211353316*coeff[7]-15.0*coeff[5])*ql[10]+30.0*coeff[5]*qc[10])*Jfac; 
  out[14] += -0.03125*((77.94228634059945*coeff[7]+33.0*coeff[5])*qr[15]+(33.0*coeff[5]-77.94228634059945*coeff[7])*ql[15]+114.0*coeff[5]*qc[15]+(77.94228634059945*coeff[6]+33.0*coeff[4])*qr[14]+(33.0*coeff[4]-77.94228634059945*coeff[6])*ql[14]+114.0*coeff[4]*qc[14]+((-45.0*coeff[7])-25.98076211353316*coeff[5])*qr[13]+(25.98076211353316*coeff[5]-45.0*coeff[7])*ql[13]+90.0*coeff[7]*qc[13]+((-45.0*coeff[6])-25.98076211353316*coeff[4])*qr[10]+(25.98076211353316*coeff[4]-45.0*coeff[6])*ql[10]+90.0*coeff[6]*qc[10])*Jfac; 
  out[15] += -0.03125*((77.94228634059945*coeff[6]+33.0*coeff[4])*qr[15]+(33.0*coeff[4]-77.94228634059945*coeff[6])*ql[15]+114.0*coeff[4]*qc[15]+(77.94228634059945*coeff[7]+33.0*coeff[5])*qr[14]+(33.0*coeff[5]-77.94228634059945*coeff[7])*ql[14]+114.0*coeff[5]*qc[14]+((-45.0*coeff[6])-25.98076211353316*coeff[4])*qr[13]+(25.98076211353316*coeff[4]-45.0*coeff[6])*ql[13]+90.0*coeff[6]*qc[13]+((-45.0*coeff[7])-25.98076211353316*coeff[5])*qr[10]+(25.98076211353316*coeff[5]-45.0*coeff[7])*ql[10]+90.0*coeff[7]*qc[10])*Jfac; 
  out[16] += -0.00625*((285.0*coeff[7]+129.9038105676658*coeff[5])*qr[20]+(285.0*coeff[7]-129.9038105676658*coeff[5])*ql[20]+330.0*coeff[7]*qc[20]+(285.0*coeff[6]+129.9038105676658*coeff[4])*qr[18]+(285.0*coeff[6]-129.9038105676658*coeff[4])*ql[18]+330.0000000000001*coeff[6]*qc[18]+((-129.9038105676658*coeff[7])-75.00000000000001*coeff[5])*qr[17]+(129.9038105676658*coeff[7]-75.00000000000001*coeff[5])*ql[17]+150.0*coeff[5]*qc[17]+((-129.9038105676658*coeff[6])-75.0*coeff[4])*qr[16]+(129.9038105676658*coeff[6]-75.0*coeff[4])*ql[16]+150.0*coeff[4]*qc[16])*Jfac; 
  out[17] += -0.00625*((285.0*coeff[6]+129.9038105676658*coeff[4])*qr[20]+(285.0*coeff[6]-129.9038105676658*coeff[4])*ql[20]+330.0000000000001*coeff[6]*qc[20]+(285.0*coeff[7]+129.9038105676658*coeff[5])*qr[18]+(285.0*coeff[7]-129.9038105676658*coeff[5])*ql[18]+330.0*coeff[7]*qc[18]+((-129.9038105676658*coeff[6])-75.0*coeff[4])*qr[17]+(129.9038105676658*coeff[6]-75.0*coeff[4])*ql[17]+150.0*coeff[4]*qc[17]+((-129.9038105676658*coeff[7])-75.00000000000001*coeff[5])*qr[16]+(129.9038105676658*coeff[7]-75.00000000000001*coeff[5])*ql[16]+150.0*coeff[5]*qc[16])*Jfac; 
  out[18] += -0.00625*((389.7114317029975*coeff[7]+165.0*coeff[5])*qr[20]+(165.0*coeff[5]-389.7114317029975*coeff[7])*ql[20]+570.0*coeff[5]*qc[20]+(389.7114317029973*coeff[6]+165.0*coeff[4])*qr[18]+(165.0*coeff[4]-389.7114317029973*coeff[6])*ql[18]+570.0*coeff[4]*qc[18]+((-225.0*coeff[7])-129.9038105676658*coeff[5])*qr[17]+(129.9038105676658*coeff[5]-225.0*coeff[7])*ql[17]+450.0*coeff[7]*qc[17]+((-225.0*coeff[6])-129.9038105676658*coeff[4])*qr[16]+(129.9038105676658*coeff[4]-225.0*coeff[6])*ql[16]+450.0000000000001*coeff[6]*qc[16])*Jfac; 
  out[19] += -0.00625*((285.0*coeff[7]+129.9038105676658*coeff[5])*qr[23]+(285.0*coeff[7]-129.9038105676658*coeff[5])*ql[23]+330.0*coeff[7]*qc[23]+(285.0*coeff[6]+129.9038105676658*coeff[4])*qr[22]+(285.0*coeff[6]-129.9038105676658*coeff[4])*ql[22]+330.0000000000001*coeff[6]*qc[22]+((-129.9038105676658*coeff[7])-75.00000000000001*coeff[5])*qr[21]+(129.9038105676658*coeff[7]-75.00000000000001*coeff[5])*ql[21]+150.0*coeff[5]*qc[21]+((-129.9038105676658*coeff[6])-75.0*coeff[4])*qr[19]+(129.9038105676658*coeff[6]-75.0*coeff[4])*ql[19]+150.0*coeff[4]*qc[19])*Jfac; 
  out[20] += -0.00625*((389.7114317029973*coeff[6]+165.0*coeff[4])*qr[20]+(165.0*coeff[4]-389.7114317029973*coeff[6])*ql[20]+570.0*coeff[4]*qc[20]+(389.7114317029975*coeff[7]+165.0*coeff[5])*qr[18]+(165.0*coeff[5]-389.7114317029975*coeff[7])*ql[18]+570.0*coeff[5]*qc[18]+((-225.0*coeff[6])-129.9038105676658*coeff[4])*qr[17]+(129.9038105676658*coeff[4]-225.0*coeff[6])*ql[17]+450.0000000000001*coeff[6]*qc[17]+((-225.0*coeff[7])-129.9038105676658*coeff[5])*qr[16]+(129.9038105676658*coeff[5]-225.0*coeff[7])*ql[16]+450.0*coeff[7]*qc[16])*Jfac; 
  out[21] += -0.00625*((285.0*coeff[6]+129.9038105676658*coeff[4])*qr[23]+(285.0*coeff[6]-129.9038105676658*coeff[4])*ql[23]+330.0000000000001*coeff[6]*qc[23]+(285.0*coeff[7]+129.9038105676658*coeff[5])*qr[22]+(285.0*coeff[7]-129.9038105676658*coeff[5])*ql[22]+330.0*coeff[7]*qc[22]+((-129.9038105676658*coeff[6])-75.0*coeff[4])*qr[21]+(129.9038105676658*coeff[6]-75.0*coeff[4])*ql[21]+150.0*coeff[4]*qc[21]+((-129.9038105676658*coeff[7])-75.00000000000001*coeff[5])*qr[19]+(129.9038105676658*coeff[7]-75.00000000000001*coeff[5])*ql[19]+150.0*coeff[5]*qc[19])*Jfac; 
  out[22] += -0.00625*((389.7114317029975*coeff[7]+165.0*coeff[5])*qr[23]+(165.0*coeff[5]-389.7114317029975*coeff[7])*ql[23]+570.0*coeff[5]*qc[23]+(389.7114317029973*coeff[6]+165.0*coeff[4])*qr[22]+(165.0*coeff[4]-389.7114317029973*coeff[6])*ql[22]+570.0*coeff[4]*qc[22]+((-225.0*coeff[7])-129.9038105676658*coeff[5])*qr[21]+(129.9038105676658*coeff[5]-225.0*coeff[7])*ql[21]+450.0*coeff[7]*qc[21]+((-225.0*coeff[6])-129.9038105676658*coeff[4])*qr[19]+(129.9038105676658*coeff[4]-225.0*coeff[6])*ql[19]+450.0000000000001*coeff[6]*qc[19])*Jfac; 
  out[23] += -0.00625*((389.7114317029973*coeff[6]+165.0*coeff[4])*qr[23]+(165.0*coeff[4]-389.7114317029973*coeff[6])*ql[23]+570.0*coeff[4]*qc[23]+(389.7114317029975*coeff[7]+165.0*coeff[5])*qr[22]+(165.0*coeff[5]-389.7114317029975*coeff[7])*ql[22]+570.0*coeff[5]*qc[22]+((-225.0*coeff[6])-129.9038105676658*coeff[4])*qr[21]+(129.9038105676658*coeff[4]-225.0*coeff[6])*ql[21]+450.0000000000001*coeff[6]*qc[21]+((-225.0*coeff[7])-129.9038105676658*coeff[5])*qr[19]+(129.9038105676658*coeff[5]-225.0*coeff[7])*ql[19]+450.0*coeff[7]*qc[19])*Jfac; 
  out[24] += -0.00625*((285.0*coeff[7]+129.9038105676658*coeff[5])*qr[28]+(285.0*coeff[7]-129.9038105676658*coeff[5])*ql[28]+330.0*coeff[7]*qc[28]+(285.0*coeff[6]+129.9038105676658*coeff[4])*qr[26]+(285.0*coeff[6]-129.9038105676658*coeff[4])*ql[26]+330.0000000000001*coeff[6]*qc[26]+((-129.9038105676658*coeff[7])-75.00000000000001*coeff[5])*qr[25]+(129.9038105676658*coeff[7]-75.00000000000001*coeff[5])*ql[25]+150.0*coeff[5]*qc[25]+((-129.9038105676658*coeff[6])-75.0*coeff[4])*qr[24]+(129.9038105676658*coeff[6]-75.0*coeff[4])*ql[24]+150.0*coeff[4]*qc[24])*Jfac; 
  out[25] += -0.00625*((285.0*coeff[6]+129.9038105676658*coeff[4])*qr[28]+(285.0*coeff[6]-129.9038105676658*coeff[4])*ql[28]+330.0000000000001*coeff[6]*qc[28]+(285.0*coeff[7]+129.9038105676658*coeff[5])*qr[26]+(285.0*coeff[7]-129.9038105676658*coeff[5])*ql[26]+330.0*coeff[7]*qc[26]+((-129.9038105676658*coeff[6])-75.0*coeff[4])*qr[25]+(129.9038105676658*coeff[6]-75.0*coeff[4])*ql[25]+150.0*coeff[4]*qc[25]+((-129.9038105676658*coeff[7])-75.00000000000001*coeff[5])*qr[24]+(129.9038105676658*coeff[7]-75.00000000000001*coeff[5])*ql[24]+150.0*coeff[5]*qc[24])*Jfac; 
  out[26] += -0.00625*((389.7114317029975*coeff[7]+165.0*coeff[5])*qr[28]+(165.0*coeff[5]-389.7114317029975*coeff[7])*ql[28]+570.0*coeff[5]*qc[28]+(389.7114317029973*coeff[6]+165.0*coeff[4])*qr[26]+(165.0*coeff[4]-389.7114317029973*coeff[6])*ql[26]+570.0*coeff[4]*qc[26]+((-225.0*coeff[7])-129.9038105676658*coeff[5])*qr[25]+(129.9038105676658*coeff[5]-225.0*coeff[7])*ql[25]+450.0*coeff[7]*qc[25]+((-225.0*coeff[6])-129.9038105676658*coeff[4])*qr[24]+(129.9038105676658*coeff[4]-225.0*coeff[6])*ql[24]+450.0000000000001*coeff[6]*qc[24])*Jfac; 
  out[27] += -0.00625*((285.0*coeff[7]+129.9038105676658*coeff[5])*qr[31]+(285.0*coeff[7]-129.9038105676658*coeff[5])*ql[31]+330.0*coeff[7]*qc[31]+(285.0*coeff[6]+129.9038105676658*coeff[4])*qr[30]+(285.0*coeff[6]-129.9038105676658*coeff[4])*ql[30]+330.0000000000001*coeff[6]*qc[30]+((-129.9038105676658*coeff[7])-75.00000000000001*coeff[5])*qr[29]+(129.9038105676658*coeff[7]-75.00000000000001*coeff[5])*ql[29]+150.0*coeff[5]*qc[29]+((-129.9038105676658*coeff[6])-75.0*coeff[4])*qr[27]+(129.9038105676658*coeff[6]-75.0*coeff[4])*ql[27]+150.0*coeff[4]*qc[27])*Jfac; 
  out[28] += -0.00625*((389.7114317029973*coeff[6]+165.0*coeff[4])*qr[28]+(165.0*coeff[4]-389.7114317029973*coeff[6])*ql[28]+570.0*coeff[4]*qc[28]+(389.7114317029975*coeff[7]+165.0*coeff[5])*qr[26]+(165.0*coeff[5]-389.7114317029975*coeff[7])*ql[26]+570.0*coeff[5]*qc[26]+((-225.0*coeff[6])-129.9038105676658*coeff[4])*qr[25]+(129.9038105676658*coeff[4]-225.0*coeff[6])*ql[25]+450.0000000000001*coeff[6]*qc[25]+((-225.0*coeff[7])-129.9038105676658*coeff[5])*qr[24]+(129.9038105676658*coeff[5]-225.0*coeff[7])*ql[24]+450.0*coeff[7]*qc[24])*Jfac; 
  out[29] += -0.00625*((285.0*coeff[6]+129.9038105676658*coeff[4])*qr[31]+(285.0*coeff[6]-129.9038105676658*coeff[4])*ql[31]+330.0000000000001*coeff[6]*qc[31]+(285.0*coeff[7]+129.9038105676658*coeff[5])*qr[30]+(285.0*coeff[7]-129.9038105676658*coeff[5])*ql[30]+330.0*coeff[7]*qc[30]+((-129.9038105676658*coeff[6])-75.0*coeff[4])*qr[29]+(129.9038105676658*coeff[6]-75.0*coeff[4])*ql[29]+150.0*coeff[4]*qc[29]+((-129.9038105676658*coeff[7])-75.00000000000001*coeff[5])*qr[27]+(129.9038105676658*coeff[7]-75.00000000000001*coeff[5])*ql[27]+150.0*coeff[5]*qc[27])*Jfac; 
  out[30] += -0.00625*((389.7114317029975*coeff[7]+165.0*coeff[5])*qr[31]+(165.0*coeff[5]-389.7114317029975*coeff[7])*ql[31]+570.0*coeff[5]*qc[31]+(389.7114317029973*coeff[6]+165.0*coeff[4])*qr[30]+(165.0*coeff[4]-389.7114317029973*coeff[6])*ql[30]+570.0*coeff[4]*qc[30]+((-225.0*coeff[7])-129.9038105676658*coeff[5])*qr[29]+(129.9038105676658*coeff[5]-225.0*coeff[7])*ql[29]+450.0*coeff[7]*qc[29]+((-225.0*coeff[6])-129.9038105676658*coeff[4])*qr[27]+(129.9038105676658*coeff[4]-225.0*coeff[6])*ql[27]+450.0000000000001*coeff[6]*qc[27])*Jfac; 
  out[31] += -0.00625*((389.7114317029973*coeff[6]+165.0*coeff[4])*qr[31]+(165.0*coeff[4]-389.7114317029973*coeff[6])*ql[31]+570.0*coeff[4]*qc[31]+(389.7114317029975*coeff[7]+165.0*coeff[5])*qr[30]+(165.0*coeff[5]-389.7114317029975*coeff[7])*ql[30]+570.0*coeff[5]*qc[30]+((-225.0*coeff[6])-129.9038105676658*coeff[4])*qr[29]+(129.9038105676658*coeff[4]-225.0*coeff[6])*ql[29]+450.0000000000001*coeff[6]*qc[29]+((-225.0*coeff[7])-129.9038105676658*coeff[5])*qr[27]+(129.9038105676658*coeff[5]-225.0*coeff[7])*ql[27]+450.0*coeff[7]*qc[27])*Jfac; 

  return 0.;

}

