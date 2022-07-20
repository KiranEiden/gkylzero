#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_stream_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *qmem, const double *f, double* GKYL_RESTRICT out) 
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
  out[5] += 3.464101615137754*f[1]*w3Ddx1+3.464101615137754*f[2]*w2Ddx0+f[8]*dv3Ddx1+f[7]*dv2Ddx0; 
  out[6] += 3.464101615137754*f[3]*w2Ddx0+(0.8944271909999159*f[16]+f[0])*dv2Ddx0; 
  out[7] += 3.464101615137754*f[3]*w3Ddx1+f[10]*dv3Ddx1; 
  out[8] += 3.464101615137754*f[4]*w2Ddx0+f[10]*dv2Ddx0; 
  out[9] += 3.464101615137754*f[4]*w3Ddx1+(0.8944271909999159*f[24]+f[0])*dv3Ddx1; 
  out[11] += 3.464101615137754*f[6]*w3Ddx1+3.464101615137754*f[7]*w2Ddx0+f[13]*dv3Ddx1+(0.8944271909999161*f[18]+f[2])*dv2Ddx0; 
  out[12] += 3.464101615137754*f[8]*w3Ddx1+3.464101615137754*f[9]*w2Ddx0+(0.8944271909999161*f[25]+f[1])*dv3Ddx1+f[14]*dv2Ddx0; 
  out[13] += 3.464101615137754*f[10]*w2Ddx0+(0.8944271909999161*f[19]+f[4])*dv2Ddx0; 
  out[14] += 3.464101615137754*f[10]*w3Ddx1+(0.8944271909999161*f[27]+f[3])*dv3Ddx1; 
  out[15] += 3.464101615137754*f[13]*w3Ddx1+3.464101615137754*f[14]*w2Ddx0+(0.8944271909999159*f[29]+f[6])*dv3Ddx1+(0.8944271909999159*f[22]+f[9])*dv2Ddx0; 
  out[17] += 3.464101615137755*f[16]*w2Ddx0+0.8944271909999161*f[3]*dv2Ddx0; 
  out[18] += 3.464101615137755*f[16]*w3Ddx1+f[19]*dv3Ddx1; 
  out[20] += 3.464101615137755*f[17]*w3Ddx1+3.464101615137755*f[18]*w2Ddx0+f[21]*dv3Ddx1+0.8944271909999159*f[7]*dv2Ddx0; 
  out[21] += 3.464101615137755*f[19]*w2Ddx0+0.8944271909999159*f[10]*dv2Ddx0; 
  out[22] += 3.464101615137755*f[19]*w3Ddx1+f[16]*dv3Ddx1; 
  out[23] += 3.464101615137755*f[21]*w3Ddx1+3.464101615137755*f[22]*w2Ddx0+f[17]*dv3Ddx1+0.8944271909999161*f[14]*dv2Ddx0; 
  out[25] += 3.464101615137755*f[24]*w2Ddx0+f[27]*dv2Ddx0; 
  out[26] += 3.464101615137755*f[24]*w3Ddx1+0.8944271909999161*f[4]*dv3Ddx1; 
  out[28] += 3.464101615137755*f[25]*w3Ddx1+3.464101615137755*f[26]*w2Ddx0+0.8944271909999159*f[8]*dv3Ddx1+f[30]*dv2Ddx0; 
  out[29] += 3.464101615137755*f[27]*w2Ddx0+f[24]*dv2Ddx0; 
  out[30] += 3.464101615137755*f[27]*w3Ddx1+0.8944271909999159*f[10]*dv3Ddx1; 
  out[31] += 3.464101615137755*f[29]*w3Ddx1+3.464101615137755*f[30]*w2Ddx0+0.8944271909999161*f[13]*dv3Ddx1+f[26]*dv2Ddx0; 

  return fabs(w2Ddx0)+0.5*dv2Ddx0+fabs(w3Ddx1)+0.5*dv3Ddx1;
} 
