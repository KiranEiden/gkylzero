#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_stream_vol_2x3v_ser_p1(const double *w, const double *dxv, const double *qmem, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // qmem:      q/m*EM fields (unused in streaming-only update).
  // f:         Input distribution function.
  // out:       Incremented output.
  double w2Ddx0  = w[2]/dxv[0]; 
  double dv2Ddx0 = dxv[2]/dxv[0]; 
  double w3Ddx1  = w[3]/dxv[1]; 
  double dv3Ddx1 = dxv[3]/dxv[1]; 

  out[1] += 3.464101615137754*f[0]*w2Ddx0+f[3]*dv2Ddx0; 
  out[2] += 3.464101615137754*f[0]*w3Ddx1+f[4]*dv3Ddx1; 
  out[6] += 3.464101615137754*f[1]*w3Ddx1+3.464101615137754*f[2]*w2Ddx0+f[9]*dv3Ddx1+f[8]*dv2Ddx0; 
  out[7] += 3.464101615137754*f[3]*w2Ddx0+(0.8944271909999159*f[32]+f[0])*dv2Ddx0; 
  out[8] += 3.464101615137754*f[3]*w3Ddx1+f[11]*dv3Ddx1; 
  out[9] += 3.464101615137754*f[4]*w2Ddx0+f[11]*dv2Ddx0; 
  out[10] += 3.464101615137754*f[4]*w3Ddx1+(0.8944271909999159*f[48]+f[0])*dv3Ddx1; 
  out[12] += 3.464101615137754*f[5]*w2Ddx0+f[14]*dv2Ddx0; 
  out[13] += 3.464101615137754*f[5]*w3Ddx1+f[15]*dv3Ddx1; 
  out[16] += 3.464101615137754*f[7]*w3Ddx1+3.464101615137754*f[8]*w2Ddx0+f[18]*dv3Ddx1+(0.8944271909999161*f[34]+f[2])*dv2Ddx0; 
  out[17] += 3.464101615137754*f[9]*w3Ddx1+3.464101615137754*f[10]*w2Ddx0+(0.8944271909999161*f[49]+f[1])*dv3Ddx1+f[19]*dv2Ddx0; 
  out[18] += 3.464101615137754*f[11]*w2Ddx0+(0.8944271909999161*f[35]+f[4])*dv2Ddx0; 
  out[19] += 3.464101615137754*f[11]*w3Ddx1+(0.8944271909999161*f[51]+f[3])*dv3Ddx1; 
  out[20] += 3.464101615137754*f[12]*w3Ddx1+3.464101615137754*f[13]*w2Ddx0+f[23]*dv3Ddx1+f[22]*dv2Ddx0; 
  out[21] += 3.464101615137754*f[14]*w2Ddx0+(0.8944271909999161*f[36]+f[5])*dv2Ddx0; 
  out[22] += 3.464101615137754*f[14]*w3Ddx1+f[25]*dv3Ddx1; 
  out[23] += 3.464101615137754*f[15]*w2Ddx0+f[25]*dv2Ddx0; 
  out[24] += 3.464101615137754*f[15]*w3Ddx1+(0.8944271909999161*f[52]+f[5])*dv3Ddx1; 
  out[26] += 3.464101615137754*f[18]*w3Ddx1+3.464101615137754*f[19]*w2Ddx0+(0.8944271909999159*f[54]+f[7])*dv3Ddx1+(0.8944271909999159*f[39]+f[10])*dv2Ddx0; 
  out[27] += 3.464101615137754*f[21]*w3Ddx1+3.464101615137754*f[22]*w2Ddx0+f[29]*dv3Ddx1+(0.8944271909999159*f[41]+f[13])*dv2Ddx0; 
  out[28] += 3.464101615137754*f[23]*w3Ddx1+3.464101615137754*f[24]*w2Ddx0+(0.8944271909999159*f[56]+f[12])*dv3Ddx1+f[30]*dv2Ddx0; 
  out[29] += 3.464101615137754*f[25]*w2Ddx0+(0.8944271909999159*f[42]+f[15])*dv2Ddx0; 
  out[30] += 3.464101615137754*f[25]*w3Ddx1+(0.8944271909999159*f[58]+f[14])*dv3Ddx1; 
  out[31] += 3.464101615137754*f[29]*w3Ddx1+3.464101615137754*f[30]*w2Ddx0+(0.8944271909999161*f[61]+f[21])*dv3Ddx1+(0.8944271909999161*f[46]+f[24])*dv2Ddx0; 
  out[33] += 3.464101615137755*f[32]*w2Ddx0+0.8944271909999161*f[3]*dv2Ddx0; 
  out[34] += 3.464101615137755*f[32]*w3Ddx1+f[35]*dv3Ddx1; 
  out[37] += 3.464101615137755*f[33]*w3Ddx1+3.464101615137755*f[34]*w2Ddx0+f[38]*dv3Ddx1+0.8944271909999159*f[8]*dv2Ddx0; 
  out[38] += 3.464101615137755*f[35]*w2Ddx0+0.8944271909999159*f[11]*dv2Ddx0; 
  out[39] += 3.464101615137755*f[35]*w3Ddx1+f[32]*dv3Ddx1; 
  out[40] += 3.464101615137755*f[36]*w2Ddx0+0.8944271909999159*f[14]*dv2Ddx0; 
  out[41] += 3.464101615137755*f[36]*w3Ddx1+f[42]*dv3Ddx1; 
  out[43] += 3.464101615137755*f[38]*w3Ddx1+3.464101615137755*f[39]*w2Ddx0+f[33]*dv3Ddx1+0.8944271909999161*f[19]*dv2Ddx0; 
  out[44] += 3.464101615137755*f[40]*w3Ddx1+3.464101615137755*f[41]*w2Ddx0+f[45]*dv3Ddx1+0.8944271909999161*f[22]*dv2Ddx0; 
  out[45] += 3.464101615137755*f[42]*w2Ddx0+0.8944271909999161*f[25]*dv2Ddx0; 
  out[46] += 3.464101615137755*f[42]*w3Ddx1+f[36]*dv3Ddx1; 
  out[47] += 3.464101615137755*f[45]*w3Ddx1+3.464101615137755*f[46]*w2Ddx0+f[40]*dv3Ddx1+0.8944271909999159*f[30]*dv2Ddx0; 
  out[49] += 3.464101615137755*f[48]*w2Ddx0+f[51]*dv2Ddx0; 
  out[50] += 3.464101615137755*f[48]*w3Ddx1+0.8944271909999161*f[4]*dv3Ddx1; 
  out[53] += 3.464101615137755*f[49]*w3Ddx1+3.464101615137755*f[50]*w2Ddx0+0.8944271909999159*f[9]*dv3Ddx1+f[55]*dv2Ddx0; 
  out[54] += 3.464101615137755*f[51]*w2Ddx0+f[48]*dv2Ddx0; 
  out[55] += 3.464101615137755*f[51]*w3Ddx1+0.8944271909999159*f[11]*dv3Ddx1; 
  out[56] += 3.464101615137755*f[52]*w2Ddx0+f[58]*dv2Ddx0; 
  out[57] += 3.464101615137755*f[52]*w3Ddx1+0.8944271909999159*f[15]*dv3Ddx1; 
  out[59] += 3.464101615137755*f[54]*w3Ddx1+3.464101615137755*f[55]*w2Ddx0+0.8944271909999161*f[18]*dv3Ddx1+f[50]*dv2Ddx0; 
  out[60] += 3.464101615137755*f[56]*w3Ddx1+3.464101615137755*f[57]*w2Ddx0+0.8944271909999161*f[23]*dv3Ddx1+f[62]*dv2Ddx0; 
  out[61] += 3.464101615137755*f[58]*w2Ddx0+f[52]*dv2Ddx0; 
  out[62] += 3.464101615137755*f[58]*w3Ddx1+0.8944271909999161*f[25]*dv3Ddx1; 
  out[63] += 3.464101615137755*f[61]*w3Ddx1+3.464101615137755*f[62]*w2Ddx0+0.8944271909999159*f[29]*dv3Ddx1+f[57]*dv2Ddx0; 
  out[65] += 3.464101615137755*f[64]*w2Ddx0+f[67]*dv2Ddx0; 
  out[66] += 3.464101615137755*f[64]*w3Ddx1+f[68]*dv3Ddx1; 
  out[69] += 3.464101615137755*f[65]*w3Ddx1+3.464101615137755*f[66]*w2Ddx0+f[72]*dv3Ddx1+f[71]*dv2Ddx0; 
  out[70] += 3.464101615137755*f[67]*w2Ddx0+f[64]*dv2Ddx0; 
  out[71] += 3.464101615137755*f[67]*w3Ddx1+f[74]*dv3Ddx1; 
  out[72] += 3.464101615137755*f[68]*w2Ddx0+f[74]*dv2Ddx0; 
  out[73] += 3.464101615137755*f[68]*w3Ddx1+f[64]*dv3Ddx1; 
  out[75] += 3.464101615137755*f[70]*w3Ddx1+3.464101615137755*f[71]*w2Ddx0+f[77]*dv3Ddx1+f[66]*dv2Ddx0; 
  out[76] += 3.464101615137755*f[72]*w3Ddx1+3.464101615137755*f[73]*w2Ddx0+f[65]*dv3Ddx1+f[78]*dv2Ddx0; 
  out[77] += 3.464101615137755*f[74]*w2Ddx0+f[68]*dv2Ddx0; 
  out[78] += 3.464101615137755*f[74]*w3Ddx1+f[67]*dv3Ddx1; 
  out[79] += 3.464101615137755*f[77]*w3Ddx1+3.464101615137755*f[78]*w2Ddx0+f[70]*dv3Ddx1+f[73]*dv2Ddx0; 

  return 3.0*(fabs(w2Ddx0)+0.5*dv2Ddx0+fabs(w3Ddx1)+0.5*dv3Ddx1);
} 
