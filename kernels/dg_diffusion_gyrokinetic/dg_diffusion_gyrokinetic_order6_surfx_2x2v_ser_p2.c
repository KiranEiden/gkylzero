#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order6_surfx_2x2v_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],6.);

  out[0] += 0.0625*(563.489130329947*coeff[0]*qr[11]+563.489130329947*coeff[0]*ql[11]-1126.978260659894*coeff[0]*qc[11]-545.5960043841961*coeff[0]*qr[1]+545.5960043841961*coeff[0]*ql[1]+315.0*coeff[0]*qr[0]+315.0*coeff[0]*ql[0]-630.0*coeff[0]*qc[0])*Jfac; 
  out[1] += 0.0078125*(6587.944671898812*coeff[0]*qr[11]-6587.944671898812*coeff[0]*ql[11]-7245.0*coeff[0]*qr[1]-7245.0*coeff[0]*ql[1]-15750.0*coeff[0]*qc[1]+4364.768035073569*coeff[0]*qr[0]-4364.768035073569*coeff[0]*ql[0])*Jfac; 
  out[2] += 0.0625*(563.4891303299469*coeff[0]*qr[19]+563.4891303299469*coeff[0]*ql[19]-1126.978260659894*coeff[0]*qc[19]-545.5960043841961*coeff[0]*qr[5]+545.5960043841961*coeff[0]*ql[5]+315.0*coeff[0]*qr[2]+315.0*coeff[0]*ql[2]-630.0*coeff[0]*qc[2])*Jfac; 
  out[3] += 0.0625*(563.4891303299469*coeff[0]*qr[21]+563.4891303299469*coeff[0]*ql[21]-1126.978260659894*coeff[0]*qc[21]-545.5960043841961*coeff[0]*qr[6]+545.5960043841961*coeff[0]*ql[6]+315.0*coeff[0]*qr[3]+315.0*coeff[0]*ql[3]-630.0*coeff[0]*qc[3])*Jfac; 
  out[4] += 0.0625*(563.4891303299469*coeff[0]*qr[25]+563.4891303299469*coeff[0]*ql[25]-1126.978260659894*coeff[0]*qc[25]-545.5960043841961*coeff[0]*qr[8]+545.5960043841961*coeff[0]*ql[8]+315.0*coeff[0]*qr[4]+315.0*coeff[0]*ql[4]-630.0*coeff[0]*qc[4])*Jfac; 
  out[5] += 0.0078125*(6587.944671898817*coeff[0]*qr[19]-6587.944671898817*coeff[0]*ql[19]-7245.0*coeff[0]*qr[5]-7245.0*coeff[0]*ql[5]-15750.0*coeff[0]*qc[5]+4364.768035073569*coeff[0]*qr[2]-4364.768035073569*coeff[0]*ql[2])*Jfac; 
  out[6] += 0.0078125*(6587.944671898817*coeff[0]*qr[21]-6587.944671898817*coeff[0]*ql[21]-7245.0*coeff[0]*qr[6]-7245.0*coeff[0]*ql[6]-15750.0*coeff[0]*qc[6]+4364.768035073569*coeff[0]*qr[3]-4364.768035073569*coeff[0]*ql[3])*Jfac; 
  out[7] += 0.0625*(563.489130329947*coeff[0]*qr[32]+563.489130329947*coeff[0]*ql[32]-1126.978260659894*coeff[0]*qc[32]-545.5960043841961*coeff[0]*qr[15]+545.5960043841961*coeff[0]*ql[15]+315.0*coeff[0]*qr[7]+315.0*coeff[0]*ql[7]-630.0*coeff[0]*qc[7])*Jfac; 
  out[8] += 0.0078125*(6587.944671898817*coeff[0]*qr[25]-6587.944671898817*coeff[0]*ql[25]-7245.0*coeff[0]*qr[8]-7245.0*coeff[0]*ql[8]-15750.0*coeff[0]*qc[8]+4364.768035073569*coeff[0]*qr[4]-4364.768035073569*coeff[0]*ql[4])*Jfac; 
  out[9] += 0.0625*(563.489130329947*coeff[0]*qr[35]+563.489130329947*coeff[0]*ql[35]-1126.978260659894*coeff[0]*qc[35]-545.5960043841961*coeff[0]*qr[16]+545.5960043841961*coeff[0]*ql[16]+315.0*coeff[0]*qr[9]+315.0*coeff[0]*ql[9]-630.0*coeff[0]*qc[9])*Jfac; 
  out[10] += 0.0625*(563.489130329947*coeff[0]*qr[37]+563.489130329947*coeff[0]*ql[37]-1126.978260659894*coeff[0]*qc[37]-545.5960043841961*coeff[0]*qr[17]+545.5960043841961*coeff[0]*ql[17]+315.0*coeff[0]*qr[10]+315.0*coeff[0]*ql[10]-630.0*coeff[0]*qc[10])*Jfac; 
  out[11] += -0.0078125*(405.0*coeff[0]*qr[11]+405.0*coeff[0]*ql[11]+18090.0*coeff[0]*qc[11]+1568.558255214003*coeff[0]*qr[1]-1568.558255214003*coeff[0]*ql[1]-1609.968943799849*coeff[0]*qr[0]-1609.968943799849*coeff[0]*ql[0]+3219.937887599698*coeff[0]*qc[0])*Jfac; 
  out[12] += -0.0625*(545.5960043841964*coeff[0]*qr[20]-545.5960043841964*coeff[0]*ql[20]-315.0*coeff[0]*qr[12]-315.0*coeff[0]*ql[12]+630.0*coeff[0]*qc[12])*Jfac; 
  out[13] += -0.0625*(545.5960043841964*coeff[0]*qr[23]-545.5960043841964*coeff[0]*ql[23]-315.0*coeff[0]*qr[13]-315.0*coeff[0]*ql[13]+630.0*coeff[0]*qc[13])*Jfac; 
  out[14] += -0.0625*(545.5960043841964*coeff[0]*qr[28]-545.5960043841964*coeff[0]*ql[28]-315.0*coeff[0]*qr[14]-315.0*coeff[0]*ql[14]+630.0*coeff[0]*qc[14])*Jfac; 
  out[15] += 0.0078125*(6587.944671898812*coeff[0]*qr[32]-6587.944671898812*coeff[0]*ql[32]-7245.0*coeff[0]*qr[15]-7245.0*coeff[0]*ql[15]-15750.0*coeff[0]*qc[15]+4364.768035073569*coeff[0]*qr[7]-4364.768035073569*coeff[0]*ql[7])*Jfac; 
  out[16] += 0.0078125*(6587.944671898812*coeff[0]*qr[35]-6587.944671898812*coeff[0]*ql[35]-7245.0*coeff[0]*qr[16]-7245.0*coeff[0]*ql[16]-15750.0*coeff[0]*qc[16]+4364.768035073569*coeff[0]*qr[9]-4364.768035073569*coeff[0]*ql[9])*Jfac; 
  out[17] += 0.0078125*(6587.944671898812*coeff[0]*qr[37]-6587.944671898812*coeff[0]*ql[37]-7245.0*coeff[0]*qr[17]-7245.0*coeff[0]*ql[17]-15750.0*coeff[0]*qc[17]+4364.768035073569*coeff[0]*qr[10]-4364.768035073569*coeff[0]*ql[10])*Jfac; 
  out[18] += 0.0625*(563.4891303299469*coeff[0]*qr[44]+563.4891303299469*coeff[0]*ql[44]-1126.978260659894*coeff[0]*qc[44]-545.5960043841961*coeff[0]*qr[31]+545.5960043841961*coeff[0]*ql[31]+315.0*coeff[0]*qr[18]+315.0*coeff[0]*ql[18]-630.0*coeff[0]*qc[18])*Jfac; 
  out[19] += -0.0078125*(405.0*coeff[0]*qr[19]+405.0*coeff[0]*ql[19]+18090.0*coeff[0]*qc[19]+1568.558255214004*coeff[0]*qr[5]-1568.558255214004*coeff[0]*ql[5]-1609.968943799848*coeff[0]*qr[2]-1609.968943799848*coeff[0]*ql[2]+3219.937887599697*coeff[0]*qc[2])*Jfac; 
  out[20] += -0.0078125*(7245.0*coeff[0]*qr[20]+7245.0*coeff[0]*ql[20]+15750.0*coeff[0]*qc[20]-4364.768035073571*coeff[0]*qr[12]+4364.768035073571*coeff[0]*ql[12])*Jfac; 
  out[21] += -0.0078125*(405.0*coeff[0]*qr[21]+405.0*coeff[0]*ql[21]+18090.0*coeff[0]*qc[21]+1568.558255214004*coeff[0]*qr[6]-1568.558255214004*coeff[0]*ql[6]-1609.968943799848*coeff[0]*qr[3]-1609.968943799848*coeff[0]*ql[3]+3219.937887599697*coeff[0]*qc[3])*Jfac; 
  out[22] += -0.0625*(545.5960043841964*coeff[0]*qr[33]-545.5960043841964*coeff[0]*ql[33]-315.0*coeff[0]*qr[22]-315.0*coeff[0]*ql[22]+630.0*coeff[0]*qc[22])*Jfac; 
  out[23] += -0.0078125*(7245.0*coeff[0]*qr[23]+7245.0*coeff[0]*ql[23]+15750.0*coeff[0]*qc[23]-4364.768035073571*coeff[0]*qr[13]+4364.768035073571*coeff[0]*ql[13])*Jfac; 
  out[24] += -0.0625*(545.5960043841964*coeff[0]*qr[34]-545.5960043841964*coeff[0]*ql[34]-315.0*coeff[0]*qr[24]-315.0*coeff[0]*ql[24]+630.0*coeff[0]*qc[24])*Jfac; 
  out[25] += -0.0078125*(405.0*coeff[0]*qr[25]+405.0*coeff[0]*ql[25]+18090.0*coeff[0]*qc[25]+1568.558255214004*coeff[0]*qr[8]-1568.558255214004*coeff[0]*ql[8]-1609.968943799848*coeff[0]*qr[4]-1609.968943799848*coeff[0]*ql[4]+3219.937887599697*coeff[0]*qc[4])*Jfac; 
  out[26] += -0.0625*(545.5960043841964*coeff[0]*qr[36]-545.5960043841964*coeff[0]*ql[36]-315.0*coeff[0]*qr[26]-315.0*coeff[0]*ql[26]+630.0*coeff[0]*qc[26])*Jfac; 
  out[27] += -0.0625*(545.5960043841964*coeff[0]*qr[39]-545.5960043841964*coeff[0]*ql[39]-315.0*coeff[0]*qr[27]-315.0*coeff[0]*ql[27]+630.0*coeff[0]*qc[27])*Jfac; 
  out[28] += -0.0078125*(7245.0*coeff[0]*qr[28]+7245.0*coeff[0]*ql[28]+15750.0*coeff[0]*qc[28]-4364.768035073571*coeff[0]*qr[14]+4364.768035073571*coeff[0]*ql[14])*Jfac; 
  out[29] += -0.0625*(545.5960043841964*coeff[0]*qr[41]-545.5960043841964*coeff[0]*ql[41]-315.0*coeff[0]*qr[29]-315.0*coeff[0]*ql[29]+630.0*coeff[0]*qc[29])*Jfac; 
  out[30] += -0.0625*(545.5960043841964*coeff[0]*qr[42]-545.5960043841964*coeff[0]*ql[42]-315.0*coeff[0]*qr[30]-315.0*coeff[0]*ql[30]+630.0*coeff[0]*qc[30])*Jfac; 
  out[31] += 0.0078125*(6587.944671898817*coeff[0]*qr[44]-6587.944671898817*coeff[0]*ql[44]-7245.0*coeff[0]*qr[31]-7245.0*coeff[0]*ql[31]-15750.0*coeff[0]*qc[31]+4364.768035073569*coeff[0]*qr[18]-4364.768035073569*coeff[0]*ql[18])*Jfac; 
  out[32] += -0.0078125*(405.0*coeff[0]*qr[32]+405.0*coeff[0]*ql[32]+18090.0*coeff[0]*qc[32]+1568.558255214003*coeff[0]*qr[15]-1568.558255214003*coeff[0]*ql[15]-1609.968943799849*coeff[0]*qr[7]-1609.968943799849*coeff[0]*ql[7]+3219.937887599698*coeff[0]*qc[7])*Jfac; 
  out[33] += -0.0078125*(7245.0*coeff[0]*qr[33]+7245.0*coeff[0]*ql[33]+15750.0*coeff[0]*qc[33]-4364.768035073571*coeff[0]*qr[22]+4364.768035073571*coeff[0]*ql[22])*Jfac; 
  out[34] += -0.0078125*(7245.0*coeff[0]*qr[34]+7245.0*coeff[0]*ql[34]+15750.0*coeff[0]*qc[34]-4364.768035073571*coeff[0]*qr[24]+4364.768035073571*coeff[0]*ql[24])*Jfac; 
  out[35] += -0.0078125*(405.0*coeff[0]*qr[35]+405.0*coeff[0]*ql[35]+18090.0*coeff[0]*qc[35]+1568.558255214003*coeff[0]*qr[16]-1568.558255214003*coeff[0]*ql[16]-1609.968943799849*coeff[0]*qr[9]-1609.968943799849*coeff[0]*ql[9]+3219.937887599698*coeff[0]*qc[9])*Jfac; 
  out[36] += -0.0078125*(7245.0*coeff[0]*qr[36]+7245.0*coeff[0]*ql[36]+15750.0*coeff[0]*qc[36]-4364.768035073571*coeff[0]*qr[26]+4364.768035073571*coeff[0]*ql[26])*Jfac; 
  out[37] += -0.0078125*(405.0*coeff[0]*qr[37]+405.0*coeff[0]*ql[37]+18090.0*coeff[0]*qc[37]+1568.558255214003*coeff[0]*qr[17]-1568.558255214003*coeff[0]*ql[17]-1609.968943799849*coeff[0]*qr[10]-1609.968943799849*coeff[0]*ql[10]+3219.937887599698*coeff[0]*qc[10])*Jfac; 
  out[38] += -0.0625*(545.5960043841964*coeff[0]*qr[45]-545.5960043841964*coeff[0]*ql[45]-315.0*coeff[0]*qr[38]-315.0*coeff[0]*ql[38]+630.0*coeff[0]*qc[38])*Jfac; 
  out[39] += -0.0078125*(7245.0*coeff[0]*qr[39]+7245.0*coeff[0]*ql[39]+15750.0*coeff[0]*qc[39]-4364.768035073571*coeff[0]*qr[27]+4364.768035073571*coeff[0]*ql[27])*Jfac; 
  out[40] += -0.0625*(545.5960043841964*coeff[0]*qr[46]-545.5960043841964*coeff[0]*ql[46]-315.0*coeff[0]*qr[40]-315.0*coeff[0]*ql[40]+630.0*coeff[0]*qc[40])*Jfac; 
  out[41] += -0.0078125*(7245.0*coeff[0]*qr[41]+7245.0*coeff[0]*ql[41]+15750.0*coeff[0]*qc[41]-4364.768035073571*coeff[0]*qr[29]+4364.768035073571*coeff[0]*ql[29])*Jfac; 
  out[42] += -0.0078125*(7245.0*coeff[0]*qr[42]+7245.0*coeff[0]*ql[42]+15750.0*coeff[0]*qc[42]-4364.768035073571*coeff[0]*qr[30]+4364.768035073571*coeff[0]*ql[30])*Jfac; 
  out[43] += -0.0625*(545.5960043841964*coeff[0]*qr[47]-545.5960043841964*coeff[0]*ql[47]-315.0*coeff[0]*qr[43]-315.0*coeff[0]*ql[43]+630.0*coeff[0]*qc[43])*Jfac; 
  out[44] += -0.0078125*(405.0*coeff[0]*qr[44]+405.0*coeff[0]*ql[44]+18090.0*coeff[0]*qc[44]+1568.558255214004*coeff[0]*qr[31]-1568.558255214004*coeff[0]*ql[31]-1609.968943799848*coeff[0]*qr[18]-1609.968943799848*coeff[0]*ql[18]+3219.937887599697*coeff[0]*qc[18])*Jfac; 
  out[45] += -0.0078125*(7245.0*coeff[0]*qr[45]+7245.0*coeff[0]*ql[45]+15750.0*coeff[0]*qc[45]-4364.768035073571*coeff[0]*qr[38]+4364.768035073571*coeff[0]*ql[38])*Jfac; 
  out[46] += -0.0078125*(7245.0*coeff[0]*qr[46]+7245.0*coeff[0]*ql[46]+15750.0*coeff[0]*qc[46]-4364.768035073571*coeff[0]*qr[40]+4364.768035073571*coeff[0]*ql[40])*Jfac; 
  out[47] += -0.0078125*(7245.0*coeff[0]*qr[47]+7245.0*coeff[0]*ql[47]+15750.0*coeff[0]*qc[47]-4364.768035073571*coeff[0]*qr[43]+4364.768035073571*coeff[0]*ql[43])*Jfac; 

  return 0.;

}

