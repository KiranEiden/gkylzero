#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_stream_vol_1x3v_ser_p2(const double *w, const double *dxv, const double *qmem, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // qmem:      q/m*EM fields (unused in streaming-only update).
  // f:         Input distribution function.
  // out:       Incremented output.
  double w1Ddx0  = w[1]/dxv[0]; 
  double dv1Ddx0 = dxv[1]/dxv[0]; 

  out[1] += 3.464101615137754*f[0]*w1Ddx0+f[2]*dv1Ddx0; 
  out[5] += 3.464101615137754*f[2]*w1Ddx0+(0.8944271909999159*f[12]+f[0])*dv1Ddx0; 
  out[6] += 3.464101615137754*f[3]*w1Ddx0+f[7]*dv1Ddx0; 
  out[8] += 3.464101615137754*f[4]*w1Ddx0+f[9]*dv1Ddx0; 
  out[11] += 7.745966692414834*f[1]*w1Ddx0+2.23606797749979*f[5]*dv1Ddx0; 
  out[15] += 3.464101615137754*f[7]*w1Ddx0+(0.8944271909999161*f[22]+f[3])*dv1Ddx0; 
  out[16] += 3.464101615137754*f[9]*w1Ddx0+(0.8944271909999161*f[26]+f[4])*dv1Ddx0; 
  out[17] += 3.464101615137754*f[10]*w1Ddx0+f[18]*dv1Ddx0; 
  out[19] += 7.745966692414834*f[5]*w1Ddx0+(2.0*f[20]+2.23606797749979*f[1])*dv1Ddx0; 
  out[20] += 3.464101615137755*f[12]*w1Ddx0+0.8944271909999161*f[2]*dv1Ddx0; 
  out[21] += 7.745966692414834*f[6]*w1Ddx0+2.23606797749979*f[15]*dv1Ddx0; 
  out[23] += 3.464101615137755*f[13]*w1Ddx0+f[24]*dv1Ddx0; 
  out[25] += 7.745966692414834*f[8]*w1Ddx0+2.23606797749979*f[16]*dv1Ddx0; 
  out[28] += 3.464101615137755*f[14]*w1Ddx0+f[29]*dv1Ddx0; 
  out[31] += 3.464101615137754*f[18]*w1Ddx0+(0.8944271909999159*f[38]+f[10])*dv1Ddx0; 
  out[32] += 7.745966692414834*f[15]*w1Ddx0+(2.0*f[33]+2.23606797749979*f[6])*dv1Ddx0; 
  out[33] += 3.464101615137755*f[22]*w1Ddx0+0.8944271909999159*f[7]*dv1Ddx0; 
  out[34] += 3.464101615137755*f[24]*w1Ddx0+f[13]*dv1Ddx0; 
  out[35] += 7.745966692414834*f[16]*w1Ddx0+(2.0*f[36]+2.23606797749979*f[8])*dv1Ddx0; 
  out[36] += 3.464101615137755*f[26]*w1Ddx0+0.8944271909999159*f[9]*dv1Ddx0; 
  out[37] += 7.745966692414834*f[17]*w1Ddx0+2.23606797749979*f[31]*dv1Ddx0; 
  out[39] += 3.464101615137755*f[27]*w1Ddx0+f[40]*dv1Ddx0; 
  out[41] += 3.464101615137755*f[29]*w1Ddx0+f[14]*dv1Ddx0; 
  out[42] += 3.464101615137755*f[30]*w1Ddx0+f[43]*dv1Ddx0; 
  out[44] += 7.745966692414834*f[31]*w1Ddx0+(2.0*f[45]+2.23606797749979*f[17])*dv1Ddx0; 
  out[45] += 3.464101615137755*f[38]*w1Ddx0+0.8944271909999161*f[18]*dv1Ddx0; 
  out[46] += 3.464101615137755*f[40]*w1Ddx0+f[27]*dv1Ddx0; 
  out[47] += 3.464101615137755*f[43]*w1Ddx0+f[30]*dv1Ddx0; 

  return 5.0*(fabs(w1Ddx0)+0.5*dv1Ddx0);
} 
