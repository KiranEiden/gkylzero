#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order6_surfy_3x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[1],6.);

  out[0] += 0.0625*(563.489130329947*coeff[1]*qr[8]+563.489130329947*coeff[1]*ql[8]-1126.978260659894*coeff[1]*qc[8]-545.5960043841961*coeff[1]*qr[2]+545.5960043841961*coeff[1]*ql[2]+(315.0*qr[0]+315.0*ql[0]-630.0*qc[0])*coeff[1])*Jfac; 
  out[1] += 0.0625*(563.4891303299469*coeff[1]*qr[12]+563.4891303299469*coeff[1]*ql[12]-1126.978260659894*coeff[1]*qc[12]-545.5960043841961*coeff[1]*qr[4]+545.5960043841961*coeff[1]*ql[4]+315.0*coeff[1]*qr[1]+315.0*coeff[1]*ql[1]-630.0*coeff[1]*qc[1])*Jfac; 
  out[2] += 0.0078125*(6587.944671898812*coeff[1]*qr[8]-6587.944671898812*coeff[1]*ql[8]-7245.0*coeff[1]*qr[2]-7245.0*coeff[1]*ql[2]-15750.0*coeff[1]*qc[2]+(4364.768035073569*qr[0]-4364.768035073569*ql[0])*coeff[1])*Jfac; 
  out[3] += 0.0625*(563.4891303299469*coeff[1]*qr[14]+563.4891303299469*coeff[1]*ql[14]-1126.978260659894*coeff[1]*qc[14]-545.5960043841961*coeff[1]*qr[6]+545.5960043841961*coeff[1]*ql[6]+315.0*coeff[1]*qr[3]+315.0*coeff[1]*ql[3]-630.0*coeff[1]*qc[3])*Jfac; 
  out[4] += 0.0078125*(6587.944671898817*coeff[1]*qr[12]-6587.944671898817*coeff[1]*ql[12]-7245.0*coeff[1]*qr[4]-7245.0*coeff[1]*ql[4]-15750.0*coeff[1]*qc[4]+4364.768035073569*coeff[1]*qr[1]-4364.768035073569*coeff[1]*ql[1])*Jfac; 
  out[5] += 0.0625*(563.489130329947*coeff[1]*qr[18]+563.489130329947*coeff[1]*ql[18]-1126.978260659894*coeff[1]*qc[18]-545.5960043841961*coeff[1]*qr[10]+545.5960043841961*coeff[1]*ql[10]+315.0*coeff[1]*qr[5]+315.0*coeff[1]*ql[5]-630.0*coeff[1]*qc[5])*Jfac; 
  out[6] += 0.0078125*(6587.944671898817*coeff[1]*qr[14]-6587.944671898817*coeff[1]*ql[14]-7245.0*coeff[1]*qr[6]-7245.0*coeff[1]*ql[6]-15750.0*coeff[1]*qc[6]+4364.768035073569*coeff[1]*qr[3]-4364.768035073569*coeff[1]*ql[3])*Jfac; 
  out[7] += -0.0625*(545.5960043841964*coeff[1]*qr[11]-545.5960043841964*coeff[1]*ql[11]-315.0*coeff[1]*qr[7]-315.0*coeff[1]*ql[7]+630.0*coeff[1]*qc[7])*Jfac; 
  out[8] += -0.0078125*(405.0*coeff[1]*qr[8]+405.0*coeff[1]*ql[8]+18090.0*coeff[1]*qc[8]+1568.558255214003*coeff[1]*qr[2]-1568.558255214003*coeff[1]*ql[2]+((-1609.968943799849*qr[0])-1609.968943799849*ql[0]+3219.937887599698*qc[0])*coeff[1])*Jfac; 
  out[9] += -0.0625*(545.5960043841964*coeff[1]*qr[16]-545.5960043841964*coeff[1]*ql[16]-315.0*coeff[1]*qr[9]-315.0*coeff[1]*ql[9]+630.0*coeff[1]*qc[9])*Jfac; 
  out[10] += 0.0078125*(6587.944671898812*coeff[1]*qr[18]-6587.944671898812*coeff[1]*ql[18]-7245.0*coeff[1]*qr[10]-7245.0*coeff[1]*ql[10]-15750.0*coeff[1]*qc[10]+4364.768035073569*coeff[1]*qr[5]-4364.768035073569*coeff[1]*ql[5])*Jfac; 
  out[11] += -0.0078125*(7245.0*coeff[1]*qr[11]+7245.0*coeff[1]*ql[11]+15750.0*coeff[1]*qc[11]-4364.768035073571*coeff[1]*qr[7]+4364.768035073571*coeff[1]*ql[7])*Jfac; 
  out[12] += -0.0078125*(405.0*coeff[1]*qr[12]+405.0*coeff[1]*ql[12]+18090.0*coeff[1]*qc[12]+1568.558255214004*coeff[1]*qr[4]-1568.558255214004*coeff[1]*ql[4]-1609.968943799848*coeff[1]*qr[1]-1609.968943799848*coeff[1]*ql[1]+3219.937887599697*coeff[1]*qc[1])*Jfac; 
  out[13] += -0.0625*(545.5960043841964*coeff[1]*qr[17]-545.5960043841964*coeff[1]*ql[17]-315.0*coeff[1]*qr[13]-315.0*coeff[1]*ql[13]+630.0*coeff[1]*qc[13])*Jfac; 
  out[14] += -0.0078125*(405.0*coeff[1]*qr[14]+405.0*coeff[1]*ql[14]+18090.0*coeff[1]*qc[14]+1568.558255214004*coeff[1]*qr[6]-1568.558255214004*coeff[1]*ql[6]-1609.968943799848*coeff[1]*qr[3]-1609.968943799848*coeff[1]*ql[3]+3219.937887599697*coeff[1]*qc[3])*Jfac; 
  out[15] += -0.0625*(545.5960043841964*coeff[1]*qr[19]-545.5960043841964*coeff[1]*ql[19]-315.0*coeff[1]*qr[15]-315.0*coeff[1]*ql[15]+630.0*coeff[1]*qc[15])*Jfac; 
  out[16] += -0.0078125*(7245.0*coeff[1]*qr[16]+7245.0*coeff[1]*ql[16]+15750.0*coeff[1]*qc[16]-4364.768035073571*coeff[1]*qr[9]+4364.768035073571*coeff[1]*ql[9])*Jfac; 
  out[17] += -0.0078125*(7245.0*coeff[1]*qr[17]+7245.0*coeff[1]*ql[17]+15750.0*coeff[1]*qc[17]-4364.768035073571*coeff[1]*qr[13]+4364.768035073571*coeff[1]*ql[13])*Jfac; 
  out[18] += -0.0078125*(405.0*coeff[1]*qr[18]+405.0*coeff[1]*ql[18]+18090.0*coeff[1]*qc[18]+1568.558255214003*coeff[1]*qr[10]-1568.558255214003*coeff[1]*ql[10]-1609.968943799849*coeff[1]*qr[5]-1609.968943799849*coeff[1]*ql[5]+3219.937887599698*coeff[1]*qc[5])*Jfac; 
  out[19] += -0.0078125*(7245.0*coeff[1]*qr[19]+7245.0*coeff[1]*ql[19]+15750.0*coeff[1]*qc[19]-4364.768035073571*coeff[1]*qr[15]+4364.768035073571*coeff[1]*ql[15])*Jfac; 

  return 0.;

}

