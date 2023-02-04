#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_stream_vol_2x2v_ser_p2(const double *w, const double *dxv, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // f:         Input distribution function.
  // out:       Incremented output.
  double w2Ddx0  = w[2]/dxv[0]; 
  double dv2Ddx0 = dxv[2]/dxv[0]; 
  double w3Ddx1  = w[3]/dxv[1]; 
  double dv3Ddx1 = dxv[3]/dxv[1]; 

  out[1] += 3.464101615137754*f[0]*w2Ddx0+f[3]*dv2Ddx0; 
  out[2] += 3.464101615137754*f[0]*w3Ddx1+f[4]*dv3Ddx1; 
  out[5] += 3.464101615137754*f[1]*w3Ddx1+3.464101615137754*f[2]*w2Ddx0+f[8]*dv3Ddx1+f[7]*dv2Ddx0; 
  out[6] += 3.464101615137754*f[3]*w2Ddx0+(0.8944271909999159*f[13]+f[0])*dv2Ddx0; 
  out[7] += 3.464101615137754*f[3]*w3Ddx1+f[10]*dv3Ddx1; 
  out[8] += 3.464101615137754*f[4]*w2Ddx0+f[10]*dv2Ddx0; 
  out[9] += 3.464101615137754*f[4]*w3Ddx1+(0.8944271909999159*f[14]+f[0])*dv3Ddx1; 
  out[11] += 7.745966692414834*f[1]*w2Ddx0+2.23606797749979*f[6]*dv2Ddx0; 
  out[12] += 7.745966692414834*f[2]*w3Ddx1+2.23606797749979*f[9]*dv3Ddx1; 
  out[15] += 3.464101615137754*f[6]*w3Ddx1+3.464101615137754*f[7]*w2Ddx0+f[17]*dv3Ddx1+(0.8944271909999161*f[24]+f[2])*dv2Ddx0; 
  out[16] += 3.464101615137754*f[8]*w3Ddx1+3.464101615137754*f[9]*w2Ddx0+(0.8944271909999161*f[28]+f[1])*dv3Ddx1+f[18]*dv2Ddx0; 
  out[17] += 3.464101615137754*f[10]*w2Ddx0+(0.8944271909999161*f[27]+f[4])*dv2Ddx0; 
  out[18] += 3.464101615137754*f[10]*w3Ddx1+(0.8944271909999161*f[30]+f[3])*dv3Ddx1; 
  out[19] += 3.464101615137755*f[11]*w3Ddx1+7.745966692414834*f[5]*w2Ddx0+f[25]*dv3Ddx1+2.23606797749979*f[15]*dv2Ddx0; 
  out[20] += 7.745966692414834*f[5]*w3Ddx1+3.464101615137755*f[12]*w2Ddx0+2.23606797749979*f[16]*dv3Ddx1+f[22]*dv2Ddx0; 
  out[21] += 7.745966692414834*f[6]*w2Ddx0+(2.0*f[23]+2.23606797749979*f[1])*dv2Ddx0; 
  out[22] += 7.745966692414834*f[7]*w3Ddx1+2.23606797749979*f[18]*dv3Ddx1; 
  out[23] += 3.464101615137755*f[13]*w2Ddx0+0.8944271909999161*f[3]*dv2Ddx0; 
  out[24] += 3.464101615137755*f[13]*w3Ddx1+f[27]*dv3Ddx1; 
  out[25] += 7.745966692414834*f[8]*w2Ddx0+2.23606797749979*f[17]*dv2Ddx0; 
  out[26] += 7.745966692414834*f[9]*w3Ddx1+(2.0*f[29]+2.23606797749979*f[2])*dv3Ddx1; 
  out[28] += 3.464101615137755*f[14]*w2Ddx0+f[30]*dv2Ddx0; 
  out[29] += 3.464101615137755*f[14]*w3Ddx1+0.8944271909999161*f[4]*dv3Ddx1; 
  out[31] += 3.464101615137754*f[17]*w3Ddx1+3.464101615137754*f[18]*w2Ddx0+(0.8944271909999159*f[42]+f[6])*dv3Ddx1+(0.8944271909999159*f[40]+f[9])*dv2Ddx0; 
  out[32] += 3.464101615137755*f[21]*w3Ddx1+7.745966692414834*f[15]*w2Ddx0+f[37]*dv3Ddx1+(2.0*f[34]+2.23606797749979*f[5])*dv2Ddx0; 
  out[33] += 7.745966692414834*f[15]*w3Ddx1+3.464101615137755*f[22]*w2Ddx0+2.23606797749979*f[31]*dv3Ddx1+f[12]*dv2Ddx0; 
  out[34] += 3.464101615137755*f[23]*w3Ddx1+3.464101615137755*f[24]*w2Ddx0+f[39]*dv3Ddx1+0.8944271909999159*f[7]*dv2Ddx0; 
  out[35] += 3.464101615137755*f[25]*w3Ddx1+7.745966692414834*f[16]*w2Ddx0+f[11]*dv3Ddx1+2.23606797749979*f[31]*dv2Ddx0; 
  out[36] += 7.745966692414834*f[16]*w3Ddx1+3.464101615137755*f[26]*w2Ddx0+(2.0*f[41]+2.23606797749979*f[5])*dv3Ddx1+f[38]*dv2Ddx0; 
  out[37] += 7.745966692414834*f[17]*w2Ddx0+(2.0*f[39]+2.23606797749979*f[8])*dv2Ddx0; 
  out[38] += 7.745966692414834*f[18]*w3Ddx1+(2.0*f[43]+2.23606797749979*f[7])*dv3Ddx1; 
  out[39] += 3.464101615137755*f[27]*w2Ddx0+0.8944271909999159*f[10]*dv2Ddx0; 
  out[40] += 3.464101615137755*f[27]*w3Ddx1+f[13]*dv3Ddx1; 
  out[41] += 3.464101615137755*f[28]*w3Ddx1+3.464101615137755*f[29]*w2Ddx0+0.8944271909999159*f[8]*dv3Ddx1+f[43]*dv2Ddx0; 
  out[42] += 3.464101615137755*f[30]*w2Ddx0+f[14]*dv2Ddx0; 
  out[43] += 3.464101615137755*f[30]*w3Ddx1+0.8944271909999159*f[10]*dv3Ddx1; 
  out[44] += 3.464101615137755*f[37]*w3Ddx1+7.745966692414834*f[31]*w2Ddx0+f[21]*dv3Ddx1+(2.0*f[46]+2.23606797749979*f[16])*dv2Ddx0; 
  out[45] += 7.745966692414834*f[31]*w3Ddx1+3.464101615137755*f[38]*w2Ddx0+(2.0*f[47]+2.23606797749979*f[15])*dv3Ddx1+f[26]*dv2Ddx0; 
  out[46] += 3.464101615137755*f[39]*w3Ddx1+3.464101615137755*f[40]*w2Ddx0+f[23]*dv3Ddx1+0.8944271909999161*f[18]*dv2Ddx0; 
  out[47] += 3.464101615137755*f[42]*w3Ddx1+3.464101615137755*f[43]*w2Ddx0+0.8944271909999161*f[17]*dv3Ddx1+f[29]*dv2Ddx0; 

  return 5.0*(fabs(w2Ddx0)+0.5*dv2Ddx0+fabs(w3Ddx1)+0.5*dv3Ddx1);
} 
