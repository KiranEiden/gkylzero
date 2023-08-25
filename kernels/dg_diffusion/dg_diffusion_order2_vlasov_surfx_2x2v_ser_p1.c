#include <gkyl_dg_diffusion_kernels.h>

GKYL_CU_DH double dg_diffusion_order2_vlasov_surfx_2x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],2.);

  out[0] += -0.0625*(8.660254037844386*coeff[0]*qr[1]-8.660254037844386*coeff[0]*ql[1]-9.0*coeff[0]*qr[0]-9.0*coeff[0]*ql[0]+18.0*coeff[0]*qc[0])*Jfac; 
  out[1] += -0.0625*(7.0*coeff[0]*qr[1]+7.0*coeff[0]*ql[1]+46.0*coeff[0]*qc[1]-8.660254037844386*coeff[0]*qr[0]+8.660254037844386*coeff[0]*ql[0])*Jfac; 
  out[2] += -0.0625*(8.660254037844386*coeff[0]*qr[5]-8.660254037844386*coeff[0]*ql[5]-9.0*coeff[0]*qr[2]-9.0*coeff[0]*ql[2]+18.0*coeff[0]*qc[2])*Jfac; 
  out[3] += -0.0625*(8.660254037844386*coeff[0]*qr[6]-8.660254037844386*coeff[0]*ql[6]-9.0*coeff[0]*qr[3]-9.0*coeff[0]*ql[3]+18.0*coeff[0]*qc[3])*Jfac; 
  out[4] += -0.0625*(8.660254037844386*coeff[0]*qr[8]-8.660254037844386*coeff[0]*ql[8]-9.0*coeff[0]*qr[4]-9.0*coeff[0]*ql[4]+18.0*coeff[0]*qc[4])*Jfac; 
  out[5] += -0.0625*(7.0*coeff[0]*qr[5]+7.0*coeff[0]*ql[5]+46.0*coeff[0]*qc[5]-8.660254037844386*coeff[0]*qr[2]+8.660254037844386*coeff[0]*ql[2])*Jfac; 
  out[6] += -0.0625*(7.0*coeff[0]*qr[6]+7.0*coeff[0]*ql[6]+46.0*coeff[0]*qc[6]-8.660254037844386*coeff[0]*qr[3]+8.660254037844386*coeff[0]*ql[3])*Jfac; 
  out[7] += -0.0625*(8.660254037844386*coeff[0]*qr[11]-8.660254037844386*coeff[0]*ql[11]-9.0*coeff[0]*qr[7]-9.0*coeff[0]*ql[7]+18.0*coeff[0]*qc[7])*Jfac; 
  out[8] += -0.0625*(7.0*coeff[0]*qr[8]+7.0*coeff[0]*ql[8]+46.0*coeff[0]*qc[8]-8.660254037844386*coeff[0]*qr[4]+8.660254037844386*coeff[0]*ql[4])*Jfac; 
  out[9] += -0.0625*(8.660254037844386*coeff[0]*qr[12]-8.660254037844386*coeff[0]*ql[12]-9.0*coeff[0]*qr[9]-9.0*coeff[0]*ql[9]+18.0*coeff[0]*qc[9])*Jfac; 
  out[10] += -0.0625*(8.660254037844386*coeff[0]*qr[13]-8.660254037844386*coeff[0]*ql[13]-9.0*coeff[0]*qr[10]-9.0*coeff[0]*ql[10]+18.0*coeff[0]*qc[10])*Jfac; 
  out[11] += -0.0625*(7.0*coeff[0]*qr[11]+7.0*coeff[0]*ql[11]+46.0*coeff[0]*qc[11]-8.660254037844386*coeff[0]*qr[7]+8.660254037844386*coeff[0]*ql[7])*Jfac; 
  out[12] += -0.0625*(7.0*coeff[0]*qr[12]+7.0*coeff[0]*ql[12]+46.0*coeff[0]*qc[12]-8.660254037844386*coeff[0]*qr[9]+8.660254037844386*coeff[0]*ql[9])*Jfac; 
  out[13] += -0.0625*(7.0*coeff[0]*qr[13]+7.0*coeff[0]*ql[13]+46.0*coeff[0]*qc[13]-8.660254037844386*coeff[0]*qr[10]+8.660254037844386*coeff[0]*ql[10])*Jfac; 
  out[14] += -0.0625*(8.660254037844386*coeff[0]*qr[15]-8.660254037844386*coeff[0]*ql[15]-9.0*coeff[0]*qr[14]-9.0*coeff[0]*ql[14]+18.0*coeff[0]*qc[14])*Jfac; 
  out[15] += -0.0625*(7.0*coeff[0]*qr[15]+7.0*coeff[0]*ql[15]+46.0*coeff[0]*qc[15]-8.660254037844386*coeff[0]*qr[14]+8.660254037844386*coeff[0]*ql[14])*Jfac; 
  out[16] += -0.0625*(8.660254037844387*coeff[0]*qr[17]-8.660254037844387*coeff[0]*ql[17]-9.0*coeff[0]*qr[16]-9.0*coeff[0]*ql[16]+18.0*coeff[0]*qc[16])*Jfac; 
  out[17] += -0.0625*(7.0*coeff[0]*qr[17]+7.0*coeff[0]*ql[17]+46.0*coeff[0]*qc[17]-8.660254037844387*coeff[0]*qr[16]+8.660254037844387*coeff[0]*ql[16])*Jfac; 
  out[18] += -0.0625*(8.660254037844387*coeff[0]*qr[20]-8.660254037844387*coeff[0]*ql[20]-9.0*coeff[0]*qr[18]-9.0*coeff[0]*ql[18]+18.0*coeff[0]*qc[18])*Jfac; 
  out[19] += -0.0625*(8.660254037844387*coeff[0]*qr[21]-8.660254037844387*coeff[0]*ql[21]-9.0*coeff[0]*qr[19]-9.0*coeff[0]*ql[19]+18.0*coeff[0]*qc[19])*Jfac; 
  out[20] += -0.0625*(7.0*coeff[0]*qr[20]+7.0*coeff[0]*ql[20]+46.0*coeff[0]*qc[20]-8.660254037844387*coeff[0]*qr[18]+8.660254037844387*coeff[0]*ql[18])*Jfac; 
  out[21] += -0.0625*(7.0*coeff[0]*qr[21]+7.0*coeff[0]*ql[21]+46.0*coeff[0]*qc[21]-8.660254037844387*coeff[0]*qr[19]+8.660254037844387*coeff[0]*ql[19])*Jfac; 
  out[22] += -0.0625*(8.660254037844387*coeff[0]*qr[23]-8.660254037844387*coeff[0]*ql[23]-9.0*coeff[0]*qr[22]-9.0*coeff[0]*ql[22]+18.0*coeff[0]*qc[22])*Jfac; 
  out[23] += -0.0625*(7.0*coeff[0]*qr[23]+7.0*coeff[0]*ql[23]+46.0*coeff[0]*qc[23]-8.660254037844387*coeff[0]*qr[22]+8.660254037844387*coeff[0]*ql[22])*Jfac; 
  out[24] += -0.0625*(8.660254037844387*coeff[0]*qr[25]-8.660254037844387*coeff[0]*ql[25]-9.0*coeff[0]*qr[24]-9.0*coeff[0]*ql[24]+18.0*coeff[0]*qc[24])*Jfac; 
  out[25] += -0.0625*(7.0*coeff[0]*qr[25]+7.0*coeff[0]*ql[25]+46.0*coeff[0]*qc[25]-8.660254037844387*coeff[0]*qr[24]+8.660254037844387*coeff[0]*ql[24])*Jfac; 
  out[26] += -0.0625*(8.660254037844387*coeff[0]*qr[28]-8.660254037844387*coeff[0]*ql[28]-9.0*coeff[0]*qr[26]-9.0*coeff[0]*ql[26]+18.0*coeff[0]*qc[26])*Jfac; 
  out[27] += -0.0625*(8.660254037844387*coeff[0]*qr[29]-8.660254037844387*coeff[0]*ql[29]-9.0*coeff[0]*qr[27]-9.0*coeff[0]*ql[27]+18.0*coeff[0]*qc[27])*Jfac; 
  out[28] += -0.0625*(7.0*coeff[0]*qr[28]+7.0*coeff[0]*ql[28]+46.0*coeff[0]*qc[28]-8.660254037844387*coeff[0]*qr[26]+8.660254037844387*coeff[0]*ql[26])*Jfac; 
  out[29] += -0.0625*(7.0*coeff[0]*qr[29]+7.0*coeff[0]*ql[29]+46.0*coeff[0]*qc[29]-8.660254037844387*coeff[0]*qr[27]+8.660254037844387*coeff[0]*ql[27])*Jfac; 
  out[30] += -0.0625*(8.660254037844387*coeff[0]*qr[31]-8.660254037844387*coeff[0]*ql[31]-9.0*coeff[0]*qr[30]-9.0*coeff[0]*ql[30]+18.0*coeff[0]*qc[30])*Jfac; 
  out[31] += -0.0625*(7.0*coeff[0]*qr[31]+7.0*coeff[0]*ql[31]+46.0*coeff[0]*qc[31]-8.660254037844387*coeff[0]*qr[30]+8.660254037844387*coeff[0]*ql[30])*Jfac; 

  return 0.;

}

GKYL_CU_DH double dg_diffusion_order2_vlasov_surfx_2x2v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],2.);

  out[0] += -0.03125*((15.0*coeff[3]+8.660254037844386*coeff[2])*qr[5]+(15.0*coeff[3]-8.660254037844386*coeff[2])*ql[5]+30.0*coeff[3]*qc[5]+(15.58845726811989*ql[2]-15.58845726811989*qr[2])*coeff[3]-9.0*coeff[2]*qr[2]-9.0*coeff[2]*ql[2]+18.0*coeff[2]*qc[2]+(15.0*coeff[1]+8.660254037844386*coeff[0])*qr[1]+(15.0*coeff[1]-8.660254037844386*coeff[0])*ql[1]+30.0*coeff[1]*qc[1]+(15.58845726811989*ql[0]-15.58845726811989*qr[0])*coeff[1]-9.0*coeff[0]*qr[0]-9.0*coeff[0]*ql[0]+18.0*coeff[0]*qc[0])*Jfac; 
  out[1] += -0.03125*((12.12435565298214*coeff[3]+7.0*coeff[2])*qr[5]+(7.0*coeff[2]-12.12435565298214*coeff[3])*ql[5]+46.0*coeff[2]*qc[5]+((-15.0*qr[2])-15.0*ql[2]+30.0*qc[2])*coeff[3]-8.660254037844386*coeff[2]*qr[2]+8.660254037844386*coeff[2]*ql[2]+(12.12435565298214*coeff[1]+7.0*coeff[0])*qr[1]+(7.0*coeff[0]-12.12435565298214*coeff[1])*ql[1]+46.0*coeff[0]*qc[1]+((-15.0*qr[0])-15.0*ql[0]+30.0*qc[0])*coeff[1]-8.660254037844386*coeff[0]*qr[0]+8.660254037844386*coeff[0]*ql[0])*Jfac; 
  out[2] += -0.03125*((15.0*coeff[1]+8.660254037844386*coeff[0])*qr[5]+(15.0*coeff[1]-8.660254037844386*coeff[0])*ql[5]+30.0*coeff[1]*qc[5]+(15.0*qr[1]+15.0*ql[1]+30.0*qc[1]-15.58845726811989*qr[0]+15.58845726811989*ql[0])*coeff[3]+((-15.58845726811989*coeff[1])-9.0*coeff[0])*qr[2]+(15.58845726811989*coeff[1]-9.0*coeff[0])*ql[2]+18.0*coeff[0]*qc[2]+(8.660254037844386*qr[1]-8.660254037844386*ql[1]-9.0*qr[0]-9.0*ql[0]+18.0*qc[0])*coeff[2])*Jfac; 
  out[3] += -0.03125*((15.0*coeff[3]+8.660254037844386*coeff[2])*qr[11]+(15.0*coeff[3]-8.660254037844386*coeff[2])*ql[11]+30.0*coeff[3]*qc[11]+((-15.58845726811989*coeff[3])-9.0*coeff[2])*qr[7]+(15.58845726811989*coeff[3]-9.0*coeff[2])*ql[7]+18.0*coeff[2]*qc[7]+(15.0*coeff[1]+8.660254037844386*coeff[0])*qr[6]+(15.0*coeff[1]-8.660254037844386*coeff[0])*ql[6]+30.0*coeff[1]*qc[6]+((-15.58845726811989*coeff[1])-9.0*coeff[0])*qr[3]+(15.58845726811989*coeff[1]-9.0*coeff[0])*ql[3]+18.0*coeff[0]*qc[3])*Jfac; 
  out[4] += -0.03125*((15.0*coeff[3]+8.660254037844386*coeff[2])*qr[12]+(15.0*coeff[3]-8.660254037844386*coeff[2])*ql[12]+30.0*coeff[3]*qc[12]+((-15.58845726811989*coeff[3])-9.0*coeff[2])*qr[9]+(15.58845726811989*coeff[3]-9.0*coeff[2])*ql[9]+18.0*coeff[2]*qc[9]+(15.0*coeff[1]+8.660254037844386*coeff[0])*qr[8]+(15.0*coeff[1]-8.660254037844386*coeff[0])*ql[8]+30.0*coeff[1]*qc[8]+((-15.58845726811989*coeff[1])-9.0*coeff[0])*qr[4]+(15.58845726811989*coeff[1]-9.0*coeff[0])*ql[4]+18.0*coeff[0]*qc[4])*Jfac; 
  out[5] += -0.03125*((12.12435565298214*coeff[1]+7.0*coeff[0])*qr[5]+(7.0*coeff[0]-12.12435565298214*coeff[1])*ql[5]+46.0*coeff[0]*qc[5]+(12.12435565298214*qr[1]-12.12435565298214*ql[1]-15.0*qr[0]-15.0*ql[0]+30.0*qc[0])*coeff[3]+((-15.0*coeff[1])-8.660254037844386*coeff[0])*qr[2]+(8.660254037844386*coeff[0]-15.0*coeff[1])*ql[2]+30.0*coeff[1]*qc[2]+(7.0*qr[1]+7.0*ql[1]+46.0*qc[1]-8.660254037844386*qr[0]+8.660254037844386*ql[0])*coeff[2])*Jfac; 
  out[6] += -0.03125*((12.12435565298214*coeff[3]+7.0*coeff[2])*qr[11]+(7.0*coeff[2]-12.12435565298214*coeff[3])*ql[11]+46.0*coeff[2]*qc[11]+((-15.0*coeff[3])-8.660254037844386*coeff[2])*qr[7]+(8.660254037844386*coeff[2]-15.0*coeff[3])*ql[7]+30.0*coeff[3]*qc[7]+(12.12435565298214*coeff[1]+7.0*coeff[0])*qr[6]+(7.0*coeff[0]-12.12435565298214*coeff[1])*ql[6]+46.0*coeff[0]*qc[6]+((-15.0*coeff[1])-8.660254037844386*coeff[0])*qr[3]+(8.660254037844386*coeff[0]-15.0*coeff[1])*ql[3]+30.0*coeff[1]*qc[3])*Jfac; 
  out[7] += -0.03125*((15.0*coeff[1]+8.660254037844386*coeff[0])*qr[11]+(15.0*coeff[1]-8.660254037844386*coeff[0])*ql[11]+30.0*coeff[1]*qc[11]+((-15.58845726811989*coeff[1])-9.0*coeff[0])*qr[7]+(15.58845726811989*coeff[1]-9.0*coeff[0])*ql[7]+18.0*coeff[0]*qc[7]+(15.0*coeff[3]+8.660254037844386*coeff[2])*qr[6]+(15.0*coeff[3]-8.660254037844386*coeff[2])*ql[6]+30.0*coeff[3]*qc[6]+((-15.58845726811989*coeff[3])-9.0*coeff[2])*qr[3]+(15.58845726811989*coeff[3]-9.0*coeff[2])*ql[3]+18.0*coeff[2]*qc[3])*Jfac; 
  out[8] += -0.03125*((12.12435565298214*coeff[3]+7.0*coeff[2])*qr[12]+(7.0*coeff[2]-12.12435565298214*coeff[3])*ql[12]+46.0*coeff[2]*qc[12]+((-15.0*coeff[3])-8.660254037844386*coeff[2])*qr[9]+(8.660254037844386*coeff[2]-15.0*coeff[3])*ql[9]+30.0*coeff[3]*qc[9]+(12.12435565298214*coeff[1]+7.0*coeff[0])*qr[8]+(7.0*coeff[0]-12.12435565298214*coeff[1])*ql[8]+46.0*coeff[0]*qc[8]+((-15.0*coeff[1])-8.660254037844386*coeff[0])*qr[4]+(8.660254037844386*coeff[0]-15.0*coeff[1])*ql[4]+30.0*coeff[1]*qc[4])*Jfac; 
  out[9] += -0.03125*((15.0*coeff[1]+8.660254037844386*coeff[0])*qr[12]+(15.0*coeff[1]-8.660254037844386*coeff[0])*ql[12]+30.0*coeff[1]*qc[12]+((-15.58845726811989*coeff[1])-9.0*coeff[0])*qr[9]+(15.58845726811989*coeff[1]-9.0*coeff[0])*ql[9]+18.0*coeff[0]*qc[9]+(15.0*coeff[3]+8.660254037844386*coeff[2])*qr[8]+(15.0*coeff[3]-8.660254037844386*coeff[2])*ql[8]+30.0*coeff[3]*qc[8]+((-15.58845726811989*coeff[3])-9.0*coeff[2])*qr[4]+(15.58845726811989*coeff[3]-9.0*coeff[2])*ql[4]+18.0*coeff[2]*qc[4])*Jfac; 
  out[10] += -0.03125*((15.0*coeff[3]+8.660254037844386*coeff[2])*qr[15]+(15.0*coeff[3]-8.660254037844386*coeff[2])*ql[15]+30.0*coeff[3]*qc[15]+((-15.58845726811989*coeff[3])-9.0*coeff[2])*qr[14]+(15.58845726811989*coeff[3]-9.0*coeff[2])*ql[14]+18.0*coeff[2]*qc[14]+(15.0*coeff[1]+8.660254037844386*coeff[0])*qr[13]+(15.0*coeff[1]-8.660254037844386*coeff[0])*ql[13]+30.0*coeff[1]*qc[13]+((-15.58845726811989*coeff[1])-9.0*coeff[0])*qr[10]+(15.58845726811989*coeff[1]-9.0*coeff[0])*ql[10]+18.0*coeff[0]*qc[10])*Jfac; 
  out[11] += -0.03125*((12.12435565298214*coeff[1]+7.0*coeff[0])*qr[11]+(7.0*coeff[0]-12.12435565298214*coeff[1])*ql[11]+46.0*coeff[0]*qc[11]+((-15.0*coeff[1])-8.660254037844386*coeff[0])*qr[7]+(8.660254037844386*coeff[0]-15.0*coeff[1])*ql[7]+30.0*coeff[1]*qc[7]+(12.12435565298214*coeff[3]+7.0*coeff[2])*qr[6]+(7.0*coeff[2]-12.12435565298214*coeff[3])*ql[6]+46.0*coeff[2]*qc[6]+((-15.0*coeff[3])-8.660254037844386*coeff[2])*qr[3]+(8.660254037844386*coeff[2]-15.0*coeff[3])*ql[3]+30.0*coeff[3]*qc[3])*Jfac; 
  out[12] += -0.03125*((12.12435565298214*coeff[1]+7.0*coeff[0])*qr[12]+(7.0*coeff[0]-12.12435565298214*coeff[1])*ql[12]+46.0*coeff[0]*qc[12]+((-15.0*coeff[1])-8.660254037844386*coeff[0])*qr[9]+(8.660254037844386*coeff[0]-15.0*coeff[1])*ql[9]+30.0*coeff[1]*qc[9]+(12.12435565298214*coeff[3]+7.0*coeff[2])*qr[8]+(7.0*coeff[2]-12.12435565298214*coeff[3])*ql[8]+46.0*coeff[2]*qc[8]+((-15.0*coeff[3])-8.660254037844386*coeff[2])*qr[4]+(8.660254037844386*coeff[2]-15.0*coeff[3])*ql[4]+30.0*coeff[3]*qc[4])*Jfac; 
  out[13] += -0.03125*((12.12435565298214*coeff[3]+7.0*coeff[2])*qr[15]+(7.0*coeff[2]-12.12435565298214*coeff[3])*ql[15]+46.0*coeff[2]*qc[15]+((-15.0*coeff[3])-8.660254037844386*coeff[2])*qr[14]+(8.660254037844386*coeff[2]-15.0*coeff[3])*ql[14]+30.0*coeff[3]*qc[14]+(12.12435565298214*coeff[1]+7.0*coeff[0])*qr[13]+(7.0*coeff[0]-12.12435565298214*coeff[1])*ql[13]+46.0*coeff[0]*qc[13]+((-15.0*coeff[1])-8.660254037844386*coeff[0])*qr[10]+(8.660254037844386*coeff[0]-15.0*coeff[1])*ql[10]+30.0*coeff[1]*qc[10])*Jfac; 
  out[14] += -0.03125*((15.0*coeff[1]+8.660254037844386*coeff[0])*qr[15]+(15.0*coeff[1]-8.660254037844386*coeff[0])*ql[15]+30.0*coeff[1]*qc[15]+((-15.58845726811989*coeff[1])-9.0*coeff[0])*qr[14]+(15.58845726811989*coeff[1]-9.0*coeff[0])*ql[14]+18.0*coeff[0]*qc[14]+(15.0*coeff[3]+8.660254037844386*coeff[2])*qr[13]+(15.0*coeff[3]-8.660254037844386*coeff[2])*ql[13]+30.0*coeff[3]*qc[13]+((-15.58845726811989*coeff[3])-9.0*coeff[2])*qr[10]+(15.58845726811989*coeff[3]-9.0*coeff[2])*ql[10]+18.0*coeff[2]*qc[10])*Jfac; 
  out[15] += -0.03125*((12.12435565298214*coeff[1]+7.0*coeff[0])*qr[15]+(7.0*coeff[0]-12.12435565298214*coeff[1])*ql[15]+46.0*coeff[0]*qc[15]+((-15.0*coeff[1])-8.660254037844386*coeff[0])*qr[14]+(8.660254037844386*coeff[0]-15.0*coeff[1])*ql[14]+30.0*coeff[1]*qc[14]+(12.12435565298214*coeff[3]+7.0*coeff[2])*qr[13]+(7.0*coeff[2]-12.12435565298214*coeff[3])*ql[13]+46.0*coeff[2]*qc[13]+((-15.0*coeff[3])-8.660254037844386*coeff[2])*qr[10]+(8.660254037844386*coeff[2]-15.0*coeff[3])*ql[10]+30.0*coeff[3]*qc[10])*Jfac; 
  out[16] += -0.00625*((75.0*coeff[3]+43.30127018922193*coeff[2])*qr[20]+(75.0*coeff[3]-43.30127018922193*coeff[2])*ql[20]+150.0*coeff[3]*qc[20]+((-77.94228634059948*coeff[3])-45.0*coeff[2])*qr[18]+(77.94228634059948*coeff[3]-45.0*coeff[2])*ql[18]+90.0*coeff[2]*qc[18]+(75.00000000000001*coeff[1]+43.30127018922195*coeff[0])*qr[17]+(75.00000000000001*coeff[1]-43.30127018922195*coeff[0])*ql[17]+150.0*coeff[1]*qc[17]+((-77.94228634059945*coeff[1])-45.0*coeff[0])*qr[16]+(77.94228634059945*coeff[1]-45.0*coeff[0])*ql[16]+90.0*coeff[0]*qc[16])*Jfac; 
  out[17] += -0.002083333333333333*((181.8653347947321*coeff[3]+105.0*coeff[2])*qr[20]+(105.0*coeff[2]-181.8653347947321*coeff[3])*ql[20]+690.0*coeff[2]*qc[20]+((-225.0*coeff[3])-129.9038105676658*coeff[2])*qr[18]+(129.9038105676658*coeff[2]-225.0*coeff[3])*ql[18]+450.0*coeff[3]*qc[18]+(181.8653347947321*coeff[1]+105.0*coeff[0])*qr[17]+(105.0*coeff[0]-181.8653347947321*coeff[1])*ql[17]+690.0*coeff[0]*qc[17]+((-225.0*coeff[1])-129.9038105676658*coeff[0])*qr[16]+(129.9038105676658*coeff[0]-225.0*coeff[1])*ql[16]+450.0000000000001*coeff[1]*qc[16])*Jfac; 
  out[18] += -0.00625*((75.00000000000001*coeff[1]+43.30127018922195*coeff[0])*qr[20]+(75.00000000000001*coeff[1]-43.30127018922195*coeff[0])*ql[20]+150.0*coeff[1]*qc[20]+((-77.94228634059945*coeff[1])-45.0*coeff[0])*qr[18]+(77.94228634059945*coeff[1]-45.0*coeff[0])*ql[18]+90.0*coeff[0]*qc[18]+(75.0*coeff[3]+43.30127018922193*coeff[2])*qr[17]+(75.0*coeff[3]-43.30127018922193*coeff[2])*ql[17]+150.0*coeff[3]*qc[17]+((-77.94228634059948*coeff[3])-45.0*coeff[2])*qr[16]+(77.94228634059948*coeff[3]-45.0*coeff[2])*ql[16]+90.0*coeff[2]*qc[16])*Jfac; 
  out[19] += -0.00625*((75.0*coeff[3]+43.30127018922193*coeff[2])*qr[23]+(75.0*coeff[3]-43.30127018922193*coeff[2])*ql[23]+150.0*coeff[3]*qc[23]+((-77.94228634059948*coeff[3])-45.0*coeff[2])*qr[22]+(77.94228634059948*coeff[3]-45.0*coeff[2])*ql[22]+90.0*coeff[2]*qc[22]+(75.00000000000001*coeff[1]+43.30127018922195*coeff[0])*qr[21]+(75.00000000000001*coeff[1]-43.30127018922195*coeff[0])*ql[21]+150.0*coeff[1]*qc[21]+((-77.94228634059945*coeff[1])-45.0*coeff[0])*qr[19]+(77.94228634059945*coeff[1]-45.0*coeff[0])*ql[19]+90.0*coeff[0]*qc[19])*Jfac; 
  out[20] += -0.002083333333333333*((181.8653347947321*coeff[1]+105.0*coeff[0])*qr[20]+(105.0*coeff[0]-181.8653347947321*coeff[1])*ql[20]+690.0*coeff[0]*qc[20]+((-225.0*coeff[1])-129.9038105676658*coeff[0])*qr[18]+(129.9038105676658*coeff[0]-225.0*coeff[1])*ql[18]+450.0000000000001*coeff[1]*qc[18]+(181.8653347947321*coeff[3]+105.0*coeff[2])*qr[17]+(105.0*coeff[2]-181.8653347947321*coeff[3])*ql[17]+690.0*coeff[2]*qc[17]+((-225.0*coeff[3])-129.9038105676658*coeff[2])*qr[16]+(129.9038105676658*coeff[2]-225.0*coeff[3])*ql[16]+450.0*coeff[3]*qc[16])*Jfac; 
  out[21] += -0.002083333333333333*((181.8653347947321*coeff[3]+105.0*coeff[2])*qr[23]+(105.0*coeff[2]-181.8653347947321*coeff[3])*ql[23]+690.0*coeff[2]*qc[23]+((-225.0*coeff[3])-129.9038105676658*coeff[2])*qr[22]+(129.9038105676658*coeff[2]-225.0*coeff[3])*ql[22]+450.0*coeff[3]*qc[22]+(181.8653347947321*coeff[1]+105.0*coeff[0])*qr[21]+(105.0*coeff[0]-181.8653347947321*coeff[1])*ql[21]+690.0*coeff[0]*qc[21]+((-225.0*coeff[1])-129.9038105676658*coeff[0])*qr[19]+(129.9038105676658*coeff[0]-225.0*coeff[1])*ql[19]+450.0000000000001*coeff[1]*qc[19])*Jfac; 
  out[22] += -0.00625*((75.00000000000001*coeff[1]+43.30127018922195*coeff[0])*qr[23]+(75.00000000000001*coeff[1]-43.30127018922195*coeff[0])*ql[23]+150.0*coeff[1]*qc[23]+((-77.94228634059945*coeff[1])-45.0*coeff[0])*qr[22]+(77.94228634059945*coeff[1]-45.0*coeff[0])*ql[22]+90.0*coeff[0]*qc[22]+(75.0*coeff[3]+43.30127018922193*coeff[2])*qr[21]+(75.0*coeff[3]-43.30127018922193*coeff[2])*ql[21]+150.0*coeff[3]*qc[21]+((-77.94228634059948*coeff[3])-45.0*coeff[2])*qr[19]+(77.94228634059948*coeff[3]-45.0*coeff[2])*ql[19]+90.0*coeff[2]*qc[19])*Jfac; 
  out[23] += -0.002083333333333333*((181.8653347947321*coeff[1]+105.0*coeff[0])*qr[23]+(105.0*coeff[0]-181.8653347947321*coeff[1])*ql[23]+690.0*coeff[0]*qc[23]+((-225.0*coeff[1])-129.9038105676658*coeff[0])*qr[22]+(129.9038105676658*coeff[0]-225.0*coeff[1])*ql[22]+450.0000000000001*coeff[1]*qc[22]+(181.8653347947321*coeff[3]+105.0*coeff[2])*qr[21]+(105.0*coeff[2]-181.8653347947321*coeff[3])*ql[21]+690.0*coeff[2]*qc[21]+((-225.0*coeff[3])-129.9038105676658*coeff[2])*qr[19]+(129.9038105676658*coeff[2]-225.0*coeff[3])*ql[19]+450.0*coeff[3]*qc[19])*Jfac; 
  out[24] += -0.00625*((75.0*coeff[3]+43.30127018922193*coeff[2])*qr[28]+(75.0*coeff[3]-43.30127018922193*coeff[2])*ql[28]+150.0*coeff[3]*qc[28]+((-77.94228634059948*coeff[3])-45.0*coeff[2])*qr[26]+(77.94228634059948*coeff[3]-45.0*coeff[2])*ql[26]+90.0*coeff[2]*qc[26]+(75.00000000000001*coeff[1]+43.30127018922195*coeff[0])*qr[25]+(75.00000000000001*coeff[1]-43.30127018922195*coeff[0])*ql[25]+150.0*coeff[1]*qc[25]+((-77.94228634059945*coeff[1])-45.0*coeff[0])*qr[24]+(77.94228634059945*coeff[1]-45.0*coeff[0])*ql[24]+90.0*coeff[0]*qc[24])*Jfac; 
  out[25] += -0.002083333333333333*((181.8653347947321*coeff[3]+105.0*coeff[2])*qr[28]+(105.0*coeff[2]-181.8653347947321*coeff[3])*ql[28]+690.0*coeff[2]*qc[28]+((-225.0*coeff[3])-129.9038105676658*coeff[2])*qr[26]+(129.9038105676658*coeff[2]-225.0*coeff[3])*ql[26]+450.0*coeff[3]*qc[26]+(181.8653347947321*coeff[1]+105.0*coeff[0])*qr[25]+(105.0*coeff[0]-181.8653347947321*coeff[1])*ql[25]+690.0*coeff[0]*qc[25]+((-225.0*coeff[1])-129.9038105676658*coeff[0])*qr[24]+(129.9038105676658*coeff[0]-225.0*coeff[1])*ql[24]+450.0000000000001*coeff[1]*qc[24])*Jfac; 
  out[26] += -0.00625*((75.00000000000001*coeff[1]+43.30127018922195*coeff[0])*qr[28]+(75.00000000000001*coeff[1]-43.30127018922195*coeff[0])*ql[28]+150.0*coeff[1]*qc[28]+((-77.94228634059945*coeff[1])-45.0*coeff[0])*qr[26]+(77.94228634059945*coeff[1]-45.0*coeff[0])*ql[26]+90.0*coeff[0]*qc[26]+(75.0*coeff[3]+43.30127018922193*coeff[2])*qr[25]+(75.0*coeff[3]-43.30127018922193*coeff[2])*ql[25]+150.0*coeff[3]*qc[25]+((-77.94228634059948*coeff[3])-45.0*coeff[2])*qr[24]+(77.94228634059948*coeff[3]-45.0*coeff[2])*ql[24]+90.0*coeff[2]*qc[24])*Jfac; 
  out[27] += -0.00625*((75.0*coeff[3]+43.30127018922193*coeff[2])*qr[31]+(75.0*coeff[3]-43.30127018922193*coeff[2])*ql[31]+150.0*coeff[3]*qc[31]+((-77.94228634059948*coeff[3])-45.0*coeff[2])*qr[30]+(77.94228634059948*coeff[3]-45.0*coeff[2])*ql[30]+90.0*coeff[2]*qc[30]+(75.00000000000001*coeff[1]+43.30127018922195*coeff[0])*qr[29]+(75.00000000000001*coeff[1]-43.30127018922195*coeff[0])*ql[29]+150.0*coeff[1]*qc[29]+((-77.94228634059945*coeff[1])-45.0*coeff[0])*qr[27]+(77.94228634059945*coeff[1]-45.0*coeff[0])*ql[27]+90.0*coeff[0]*qc[27])*Jfac; 
  out[28] += -0.002083333333333333*((181.8653347947321*coeff[1]+105.0*coeff[0])*qr[28]+(105.0*coeff[0]-181.8653347947321*coeff[1])*ql[28]+690.0*coeff[0]*qc[28]+((-225.0*coeff[1])-129.9038105676658*coeff[0])*qr[26]+(129.9038105676658*coeff[0]-225.0*coeff[1])*ql[26]+450.0000000000001*coeff[1]*qc[26]+(181.8653347947321*coeff[3]+105.0*coeff[2])*qr[25]+(105.0*coeff[2]-181.8653347947321*coeff[3])*ql[25]+690.0*coeff[2]*qc[25]+((-225.0*coeff[3])-129.9038105676658*coeff[2])*qr[24]+(129.9038105676658*coeff[2]-225.0*coeff[3])*ql[24]+450.0*coeff[3]*qc[24])*Jfac; 
  out[29] += -0.002083333333333333*((181.8653347947321*coeff[3]+105.0*coeff[2])*qr[31]+(105.0*coeff[2]-181.8653347947321*coeff[3])*ql[31]+690.0*coeff[2]*qc[31]+((-225.0*coeff[3])-129.9038105676658*coeff[2])*qr[30]+(129.9038105676658*coeff[2]-225.0*coeff[3])*ql[30]+450.0*coeff[3]*qc[30]+(181.8653347947321*coeff[1]+105.0*coeff[0])*qr[29]+(105.0*coeff[0]-181.8653347947321*coeff[1])*ql[29]+690.0*coeff[0]*qc[29]+((-225.0*coeff[1])-129.9038105676658*coeff[0])*qr[27]+(129.9038105676658*coeff[0]-225.0*coeff[1])*ql[27]+450.0000000000001*coeff[1]*qc[27])*Jfac; 
  out[30] += -0.00625*((75.00000000000001*coeff[1]+43.30127018922195*coeff[0])*qr[31]+(75.00000000000001*coeff[1]-43.30127018922195*coeff[0])*ql[31]+150.0*coeff[1]*qc[31]+((-77.94228634059945*coeff[1])-45.0*coeff[0])*qr[30]+(77.94228634059945*coeff[1]-45.0*coeff[0])*ql[30]+90.0*coeff[0]*qc[30]+(75.0*coeff[3]+43.30127018922193*coeff[2])*qr[29]+(75.0*coeff[3]-43.30127018922193*coeff[2])*ql[29]+150.0*coeff[3]*qc[29]+((-77.94228634059948*coeff[3])-45.0*coeff[2])*qr[27]+(77.94228634059948*coeff[3]-45.0*coeff[2])*ql[27]+90.0*coeff[2]*qc[27])*Jfac; 
  out[31] += -0.002083333333333333*((181.8653347947321*coeff[1]+105.0*coeff[0])*qr[31]+(105.0*coeff[0]-181.8653347947321*coeff[1])*ql[31]+690.0*coeff[0]*qc[31]+((-225.0*coeff[1])-129.9038105676658*coeff[0])*qr[30]+(129.9038105676658*coeff[0]-225.0*coeff[1])*ql[30]+450.0000000000001*coeff[1]*qc[30]+(181.8653347947321*coeff[3]+105.0*coeff[2])*qr[29]+(105.0*coeff[2]-181.8653347947321*coeff[3])*ql[29]+690.0*coeff[2]*qc[29]+((-225.0*coeff[3])-129.9038105676658*coeff[2])*qr[27]+(129.9038105676658*coeff[2]-225.0*coeff[3])*ql[27]+450.0*coeff[3]*qc[27])*Jfac; 

  return 0.;

}

