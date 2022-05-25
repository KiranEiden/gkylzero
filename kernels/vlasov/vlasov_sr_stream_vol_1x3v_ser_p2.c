#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_sr_stream_vol_1x3v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // p_over_gamma: p/gamma (velocity).
  // qmem:      q/m*EM fields (unused in streaming-only update).
  // f:         Input distribution function.
  // out:       Incremented output.
  const double dx10 = 2/dxv[0]; 
  const double *p0_over_gamma = &p_over_gamma[0]; 
  const double *p1_over_gamma = &p_over_gamma[20]; 
  const double *p2_over_gamma = &p_over_gamma[40]; 

  double alpha_mid = 0.0; 
  double alpha_cdim[48] = {0.0}; 
  alpha_cdim[0] = 1.414213562373095*p0_over_gamma[0]*dx10; 
  alpha_cdim[2] = 1.414213562373095*p0_over_gamma[1]*dx10; 
  alpha_cdim[3] = 1.414213562373095*p0_over_gamma[2]*dx10; 
  alpha_cdim[4] = 1.414213562373095*p0_over_gamma[3]*dx10; 
  alpha_cdim[7] = 1.414213562373095*p0_over_gamma[4]*dx10; 
  alpha_cdim[9] = 1.414213562373095*p0_over_gamma[5]*dx10; 
  alpha_cdim[10] = 1.414213562373095*p0_over_gamma[6]*dx10; 
  alpha_cdim[12] = 1.414213562373095*p0_over_gamma[7]*dx10; 
  alpha_cdim[13] = 1.414213562373095*p0_over_gamma[8]*dx10; 
  alpha_cdim[14] = 1.414213562373095*p0_over_gamma[9]*dx10; 
  alpha_cdim[18] = 1.414213562373095*p0_over_gamma[10]*dx10; 
  alpha_cdim[22] = 1.414213562373095*p0_over_gamma[11]*dx10; 
  alpha_cdim[24] = 1.414213562373095*p0_over_gamma[12]*dx10; 
  alpha_cdim[26] = 1.414213562373095*p0_over_gamma[13]*dx10; 
  alpha_cdim[27] = 1.414213562373095*p0_over_gamma[14]*dx10; 
  alpha_cdim[29] = 1.414213562373095*p0_over_gamma[15]*dx10; 
  alpha_cdim[30] = 1.414213562373095*p0_over_gamma[16]*dx10; 
  alpha_cdim[38] = 1.414213562373095*p0_over_gamma[17]*dx10; 
  alpha_cdim[40] = 1.414213562373095*p0_over_gamma[18]*dx10; 
  alpha_cdim[43] = 1.414213562373095*p0_over_gamma[19]*dx10; 
  alpha_mid += fabs(0.125*alpha_cdim[0]-0.1397542485937369*(alpha_cdim[14]+alpha_cdim[13]+alpha_cdim[12])); 

  out[1] += 0.4330127018922193*(alpha_cdim[43]*f[43]+alpha_cdim[40]*f[40]+alpha_cdim[38]*f[38]+alpha_cdim[30]*f[30]+alpha_cdim[29]*f[29]+alpha_cdim[27]*f[27]+alpha_cdim[26]*f[26]+alpha_cdim[24]*f[24]+alpha_cdim[22]*f[22]+alpha_cdim[18]*f[18]+alpha_cdim[14]*f[14]+alpha_cdim[13]*f[13]+alpha_cdim[12]*f[12]+alpha_cdim[10]*f[10]+alpha_cdim[9]*f[9]+alpha_cdim[7]*f[7]+alpha_cdim[4]*f[4]+alpha_cdim[3]*f[3]+alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[5] += 0.4330127018922193*(alpha_cdim[30]*f[43]+f[30]*alpha_cdim[43]+alpha_cdim[27]*f[40]+f[27]*alpha_cdim[40])+0.3872983346207416*(alpha_cdim[18]*f[38]+f[18]*alpha_cdim[38])+0.4330127018922193*(alpha_cdim[14]*f[29]+f[14]*alpha_cdim[29])+0.3872983346207416*(alpha_cdim[9]*f[26]+f[9]*alpha_cdim[26])+0.4330127018922193*(alpha_cdim[13]*f[24]+f[13]*alpha_cdim[24])+0.3872983346207416*(alpha_cdim[7]*f[22]+f[7]*alpha_cdim[22])+0.4330127018922193*(alpha_cdim[10]*f[18]+f[10]*alpha_cdim[18])+0.3872983346207416*(alpha_cdim[2]*f[12]+f[2]*alpha_cdim[12])+0.4330127018922193*(alpha_cdim[4]*f[9]+f[4]*alpha_cdim[9]+alpha_cdim[3]*f[7]+f[3]*alpha_cdim[7]+alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]); 
  out[6] += 0.4330127018922193*(alpha_cdim[29]*f[43]+f[29]*alpha_cdim[43])+0.3872983346207416*(alpha_cdim[18]*f[40]+f[18]*alpha_cdim[40])+0.4330127018922193*(alpha_cdim[26]*f[38]+f[26]*alpha_cdim[38]+alpha_cdim[14]*f[30]+f[14]*alpha_cdim[30])+0.3872983346207416*(alpha_cdim[10]*f[27]+f[10]*alpha_cdim[27]+alpha_cdim[7]*f[24]+f[7]*alpha_cdim[24])+0.4330127018922193*(alpha_cdim[12]*f[22]+f[12]*alpha_cdim[22]+alpha_cdim[9]*f[18]+f[9]*alpha_cdim[18])+0.3872983346207416*(alpha_cdim[3]*f[13]+f[3]*alpha_cdim[13])+0.4330127018922193*(alpha_cdim[4]*f[10]+f[4]*alpha_cdim[10]+alpha_cdim[2]*f[7]+f[2]*alpha_cdim[7]+alpha_cdim[0]*f[3]+f[0]*alpha_cdim[3]); 
  out[8] += 0.3872983346207416*(alpha_cdim[18]*f[43]+f[18]*alpha_cdim[43])+0.4330127018922193*(alpha_cdim[24]*f[40]+f[24]*alpha_cdim[40]+alpha_cdim[22]*f[38]+f[22]*alpha_cdim[38])+0.3872983346207416*(alpha_cdim[10]*f[30]+f[10]*alpha_cdim[30]+alpha_cdim[9]*f[29]+f[9]*alpha_cdim[29])+0.4330127018922193*(alpha_cdim[13]*f[27]+f[13]*alpha_cdim[27]+alpha_cdim[12]*f[26]+f[12]*alpha_cdim[26]+alpha_cdim[7]*f[18]+f[7]*alpha_cdim[18])+0.3872983346207416*(alpha_cdim[4]*f[14]+f[4]*alpha_cdim[14])+0.4330127018922193*(alpha_cdim[3]*f[10]+f[3]*alpha_cdim[10]+alpha_cdim[2]*f[9]+f[2]*alpha_cdim[9]+alpha_cdim[0]*f[4]+f[0]*alpha_cdim[4]); 
  out[11] += 0.9682458365518543*(alpha_cdim[43]*f[47]+alpha_cdim[40]*f[46]+alpha_cdim[38]*f[45]+alpha_cdim[30]*f[42]+alpha_cdim[29]*f[41]+alpha_cdim[27]*f[39]+alpha_cdim[26]*f[36]+alpha_cdim[24]*f[34]+alpha_cdim[22]*f[33]+alpha_cdim[18]*f[31]+alpha_cdim[14]*f[28]+alpha_cdim[13]*f[23]+alpha_cdim[12]*f[20]+alpha_cdim[10]*f[17]+alpha_cdim[9]*f[16]+alpha_cdim[7]*f[15]+alpha_cdim[4]*f[8]+alpha_cdim[3]*f[6]+alpha_cdim[2]*f[5]+alpha_cdim[0]*f[1]); 
  out[15] += 0.4330127018922193*(alpha_cdim[14]*f[43]+f[14]*alpha_cdim[43])+(0.3464101615137755*alpha_cdim[38]+0.3872983346207416*alpha_cdim[10])*f[40]+0.3464101615137755*f[38]*alpha_cdim[40]+0.3872983346207416*(f[10]*alpha_cdim[40]+alpha_cdim[9]*f[38]+f[9]*alpha_cdim[38])+0.4330127018922193*(alpha_cdim[29]*f[30]+f[29]*alpha_cdim[30])+0.3872983346207416*(alpha_cdim[18]*f[27]+f[18]*alpha_cdim[27]+alpha_cdim[18]*f[26]+f[18]*alpha_cdim[26])+(0.3464101615137755*alpha_cdim[22]+0.3872983346207416*alpha_cdim[3])*f[24]+0.3464101615137755*f[22]*alpha_cdim[24]+0.3872983346207416*(f[3]*alpha_cdim[24]+alpha_cdim[2]*f[22]+f[2]*alpha_cdim[22])+0.4330127018922193*(alpha_cdim[4]*f[18]+f[4]*alpha_cdim[18])+0.3872983346207416*(alpha_cdim[7]*f[13]+f[7]*alpha_cdim[13]+alpha_cdim[7]*f[12]+f[7]*alpha_cdim[12])+0.4330127018922193*(alpha_cdim[9]*f[10]+f[9]*alpha_cdim[10]+alpha_cdim[0]*f[7]+f[0]*alpha_cdim[7]+alpha_cdim[2]*f[3]+f[2]*alpha_cdim[3]); 
  out[16] += (0.3464101615137755*alpha_cdim[38]+0.3872983346207416*alpha_cdim[10])*f[43]+(0.3464101615137755*f[38]+0.3872983346207416*f[10])*alpha_cdim[43]+0.4330127018922193*(alpha_cdim[13]*f[40]+f[13]*alpha_cdim[40])+0.3872983346207416*(alpha_cdim[7]*f[38]+f[7]*alpha_cdim[38]+alpha_cdim[18]*f[30]+f[18]*alpha_cdim[30])+(0.3464101615137755*alpha_cdim[26]+0.3872983346207416*alpha_cdim[4])*f[29]+(0.3464101615137755*f[26]+0.3872983346207416*f[4])*alpha_cdim[29]+0.4330127018922193*(alpha_cdim[24]*f[27]+f[24]*alpha_cdim[27])+0.3872983346207416*(alpha_cdim[2]*f[26]+f[2]*alpha_cdim[26]+alpha_cdim[18]*f[22]+f[18]*alpha_cdim[22])+0.4330127018922193*(alpha_cdim[3]*f[18]+f[3]*alpha_cdim[18])+0.3872983346207416*(alpha_cdim[9]*f[14]+f[9]*alpha_cdim[14]+alpha_cdim[9]*f[12]+f[9]*alpha_cdim[12])+0.4330127018922193*(alpha_cdim[7]*f[10]+f[7]*alpha_cdim[10]+alpha_cdim[0]*f[9]+f[0]*alpha_cdim[9]+alpha_cdim[2]*f[4]+f[2]*alpha_cdim[4]); 
  out[17] += (0.3464101615137755*alpha_cdim[40]+0.3872983346207416*alpha_cdim[9])*f[43]+0.3464101615137755*f[40]*alpha_cdim[43]+0.3872983346207416*(f[9]*alpha_cdim[43]+alpha_cdim[7]*f[40]+f[7]*alpha_cdim[40])+0.4330127018922193*(alpha_cdim[12]*f[38]+f[12]*alpha_cdim[38])+(0.3464101615137755*alpha_cdim[27]+0.3872983346207416*alpha_cdim[4])*f[30]+0.3464101615137755*f[27]*alpha_cdim[30]+0.3872983346207416*(f[4]*alpha_cdim[30]+alpha_cdim[18]*f[29]+f[18]*alpha_cdim[29]+alpha_cdim[3]*f[27]+f[3]*alpha_cdim[27])+0.4330127018922193*(alpha_cdim[22]*f[26]+f[22]*alpha_cdim[26])+0.3872983346207416*(alpha_cdim[18]*f[24]+f[18]*alpha_cdim[24])+0.4330127018922193*(alpha_cdim[2]*f[18]+f[2]*alpha_cdim[18])+0.3872983346207416*(alpha_cdim[10]*f[14]+f[10]*alpha_cdim[14]+alpha_cdim[10]*f[13]+f[10]*alpha_cdim[13])+0.4330127018922193*(alpha_cdim[0]*f[10]+f[0]*alpha_cdim[10]+alpha_cdim[7]*f[9]+f[7]*alpha_cdim[9]+alpha_cdim[3]*f[4]+f[3]*alpha_cdim[4]); 
  out[19] += 0.9682458365518543*(alpha_cdim[30]*f[47]+alpha_cdim[27]*f[46])+0.8660254037844386*alpha_cdim[18]*f[45]+0.9682458365518543*(f[42]*alpha_cdim[43]+alpha_cdim[14]*f[41]+f[39]*alpha_cdim[40])+0.8660254037844386*(f[31]*alpha_cdim[38]+alpha_cdim[9]*f[36])+0.9682458365518543*alpha_cdim[13]*f[34]+0.8660254037844386*alpha_cdim[7]*f[33]+0.9682458365518543*(alpha_cdim[10]*f[31]+f[28]*alpha_cdim[29])+0.8660254037844386*f[16]*alpha_cdim[26]+0.9682458365518543*f[23]*alpha_cdim[24]+0.8660254037844386*(f[15]*alpha_cdim[22]+alpha_cdim[2]*f[20])+0.9682458365518543*(f[17]*alpha_cdim[18]+alpha_cdim[4]*f[16]+alpha_cdim[3]*f[15])+0.8660254037844386*f[5]*alpha_cdim[12]+0.9682458365518543*(f[8]*alpha_cdim[9]+f[6]*alpha_cdim[7]+alpha_cdim[0]*f[5]+f[1]*alpha_cdim[2]); 
  out[20] += 0.3872983346207416*(alpha_cdim[43]*f[43]+alpha_cdim[40]*f[40])+0.276641667586244*alpha_cdim[38]*f[38]+0.4330127018922193*(alpha_cdim[10]*f[38]+f[10]*alpha_cdim[38])+0.3872983346207416*alpha_cdim[29]*f[29]+0.276641667586244*alpha_cdim[26]*f[26]+0.4330127018922193*(alpha_cdim[4]*f[26]+f[4]*alpha_cdim[26])+0.3872983346207416*alpha_cdim[24]*f[24]+0.276641667586244*alpha_cdim[22]*f[22]+0.4330127018922193*(alpha_cdim[3]*f[22]+f[3]*alpha_cdim[22])+0.3872983346207416*alpha_cdim[18]*f[18]+0.276641667586244*alpha_cdim[12]*f[12]+0.4330127018922193*(alpha_cdim[0]*f[12]+f[0]*alpha_cdim[12])+0.3872983346207416*(alpha_cdim[9]*f[9]+alpha_cdim[7]*f[7]+alpha_cdim[2]*f[2]); 
  out[21] += 0.9682458365518543*alpha_cdim[29]*f[47]+0.8660254037844386*alpha_cdim[18]*f[46]+0.9682458365518543*(alpha_cdim[26]*f[45]+f[41]*alpha_cdim[43]+alpha_cdim[14]*f[42])+0.8660254037844386*(f[31]*alpha_cdim[40]+alpha_cdim[10]*f[39])+0.9682458365518543*f[36]*alpha_cdim[38]+0.8660254037844386*alpha_cdim[7]*f[34]+0.9682458365518543*(alpha_cdim[12]*f[33]+alpha_cdim[9]*f[31]+f[28]*alpha_cdim[30])+0.8660254037844386*(f[17]*alpha_cdim[27]+f[15]*alpha_cdim[24]+alpha_cdim[3]*f[23])+0.9682458365518543*(f[20]*alpha_cdim[22]+f[16]*alpha_cdim[18]+alpha_cdim[4]*f[17]+alpha_cdim[2]*f[15])+0.8660254037844386*f[6]*alpha_cdim[13]+0.9682458365518543*(f[8]*alpha_cdim[10]+f[5]*alpha_cdim[7]+alpha_cdim[0]*f[6]+f[1]*alpha_cdim[3]); 
  out[23] += 0.3872983346207416*alpha_cdim[43]*f[43]+0.276641667586244*alpha_cdim[40]*f[40]+0.4330127018922193*(alpha_cdim[9]*f[40]+f[9]*alpha_cdim[40])+0.3872983346207416*(alpha_cdim[38]*f[38]+alpha_cdim[30]*f[30])+0.276641667586244*alpha_cdim[27]*f[27]+0.4330127018922193*(alpha_cdim[4]*f[27]+f[4]*alpha_cdim[27])+0.276641667586244*alpha_cdim[24]*f[24]+0.4330127018922193*(alpha_cdim[2]*f[24]+f[2]*alpha_cdim[24])+0.3872983346207416*(alpha_cdim[22]*f[22]+alpha_cdim[18]*f[18])+0.276641667586244*alpha_cdim[13]*f[13]+0.4330127018922193*(alpha_cdim[0]*f[13]+f[0]*alpha_cdim[13])+0.3872983346207416*(alpha_cdim[10]*f[10]+alpha_cdim[7]*f[7]+alpha_cdim[3]*f[3]); 
  out[25] += 0.8660254037844386*alpha_cdim[18]*f[47]+0.9682458365518543*(alpha_cdim[24]*f[46]+alpha_cdim[22]*f[45])+0.8660254037844386*(f[31]*alpha_cdim[43]+alpha_cdim[10]*f[42]+alpha_cdim[9]*f[41])+0.9682458365518543*(f[34]*alpha_cdim[40]+alpha_cdim[13]*f[39]+f[33]*alpha_cdim[38]+alpha_cdim[12]*f[36]+alpha_cdim[7]*f[31])+0.8660254037844386*(f[17]*alpha_cdim[30]+f[16]*alpha_cdim[29]+alpha_cdim[4]*f[28])+0.9682458365518543*(f[23]*alpha_cdim[27]+f[20]*alpha_cdim[26]+f[15]*alpha_cdim[18]+alpha_cdim[3]*f[17]+alpha_cdim[2]*f[16])+0.8660254037844386*f[8]*alpha_cdim[14]+0.9682458365518543*(f[6]*alpha_cdim[10]+f[5]*alpha_cdim[9]+alpha_cdim[0]*f[8]+f[1]*alpha_cdim[4]); 
  out[28] += 0.276641667586244*alpha_cdim[43]*f[43]+0.4330127018922193*(alpha_cdim[7]*f[43]+f[7]*alpha_cdim[43])+0.3872983346207416*(alpha_cdim[40]*f[40]+alpha_cdim[38]*f[38])+0.276641667586244*alpha_cdim[30]*f[30]+0.4330127018922193*(alpha_cdim[3]*f[30]+f[3]*alpha_cdim[30])+0.276641667586244*alpha_cdim[29]*f[29]+0.4330127018922193*(alpha_cdim[2]*f[29]+f[2]*alpha_cdim[29])+0.3872983346207416*(alpha_cdim[27]*f[27]+alpha_cdim[26]*f[26]+alpha_cdim[18]*f[18])+0.276641667586244*alpha_cdim[14]*f[14]+0.4330127018922193*(alpha_cdim[0]*f[14]+f[0]*alpha_cdim[14])+0.3872983346207416*(alpha_cdim[10]*f[10]+alpha_cdim[9]*f[9]+alpha_cdim[4]*f[4]); 
  out[31] += (0.3464101615137755*(alpha_cdim[27]+alpha_cdim[26])+0.3872983346207416*alpha_cdim[4])*f[43]+(0.3464101615137755*(f[27]+f[26])+0.3872983346207416*f[4])*alpha_cdim[43]+(0.3464101615137755*(alpha_cdim[30]+alpha_cdim[22])+0.3872983346207416*alpha_cdim[3])*f[40]+(0.3464101615137755*(f[30]+f[22])+0.3872983346207416*f[3])*alpha_cdim[40]+(0.3464101615137755*(alpha_cdim[29]+alpha_cdim[24])+0.3872983346207416*alpha_cdim[2])*f[38]+0.3464101615137755*(f[29]+f[24])*alpha_cdim[38]+0.3872983346207416*(f[2]*alpha_cdim[38]+alpha_cdim[9]*f[30]+f[9]*alpha_cdim[30]+alpha_cdim[10]*f[29]+f[10]*alpha_cdim[29]+alpha_cdim[7]*f[27]+f[7]*alpha_cdim[27]+alpha_cdim[7]*f[26]+f[7]*alpha_cdim[26]+alpha_cdim[10]*f[24]+f[10]*alpha_cdim[24]+alpha_cdim[9]*f[22]+f[9]*alpha_cdim[22])+(0.3872983346207416*(alpha_cdim[14]+alpha_cdim[13]+alpha_cdim[12])+0.4330127018922193*alpha_cdim[0])*f[18]+0.3872983346207416*(f[14]+f[13]+f[12])*alpha_cdim[18]+0.4330127018922193*(f[0]*alpha_cdim[18]+alpha_cdim[2]*f[10]+f[2]*alpha_cdim[10]+alpha_cdim[3]*f[9]+f[3]*alpha_cdim[9]+alpha_cdim[4]*f[7]+f[4]*alpha_cdim[7]); 
  out[32] += 0.9682458365518543*alpha_cdim[14]*f[47]+(0.7745966692414833*alpha_cdim[38]+0.8660254037844386*alpha_cdim[10])*f[46]+(0.7745966692414833*alpha_cdim[40]+0.8660254037844386*alpha_cdim[9])*f[45]+0.9682458365518543*(f[28]*alpha_cdim[43]+alpha_cdim[29]*f[42]+alpha_cdim[30]*f[41])+0.8660254037844386*(f[17]*alpha_cdim[40]+alpha_cdim[18]*f[39]+f[16]*alpha_cdim[38]+alpha_cdim[18]*f[36])+(0.7745966692414833*alpha_cdim[22]+0.8660254037844386*alpha_cdim[3])*f[34]+(0.7745966692414833*alpha_cdim[24]+0.8660254037844386*alpha_cdim[2])*f[33]+(0.8660254037844386*(alpha_cdim[27]+alpha_cdim[26])+0.9682458365518543*alpha_cdim[4])*f[31]+0.8660254037844386*(f[6]*alpha_cdim[24]+alpha_cdim[7]*f[23]+f[5]*alpha_cdim[22]+alpha_cdim[7]*f[20])+0.9682458365518543*(f[8]*alpha_cdim[18]+alpha_cdim[9]*f[17]+alpha_cdim[10]*f[16])+0.8660254037844386*(alpha_cdim[13]+alpha_cdim[12])*f[15]+0.9682458365518543*(alpha_cdim[0]*f[15]+f[1]*alpha_cdim[7]+alpha_cdim[2]*f[6]+alpha_cdim[3]*f[5]); 
  out[33] += 0.3872983346207416*(alpha_cdim[29]*f[43]+f[29]*alpha_cdim[43])+0.3464101615137755*(alpha_cdim[18]*f[40]+f[18]*alpha_cdim[40])+(0.3872983346207416*alpha_cdim[27]+0.276641667586244*alpha_cdim[26]+0.4330127018922193*alpha_cdim[4])*f[38]+(0.3872983346207416*f[27]+0.276641667586244*f[26])*alpha_cdim[38]+0.4330127018922193*(f[4]*alpha_cdim[38]+alpha_cdim[10]*f[26]+f[10]*alpha_cdim[26])+0.3464101615137755*(alpha_cdim[7]*f[24]+f[7]*alpha_cdim[24])+(0.3872983346207416*alpha_cdim[13]+0.276641667586244*alpha_cdim[12]+0.4330127018922193*alpha_cdim[0])*f[22]+(0.3872983346207416*f[13]+0.276641667586244*f[12]+0.4330127018922193*f[0])*alpha_cdim[22]+0.3872983346207416*(alpha_cdim[9]*f[18]+f[9]*alpha_cdim[18])+0.4330127018922193*(alpha_cdim[3]*f[12]+f[3]*alpha_cdim[12])+0.3872983346207416*(alpha_cdim[2]*f[7]+f[2]*alpha_cdim[7]); 
  out[34] += 0.3872983346207416*(alpha_cdim[30]*f[43]+f[30]*alpha_cdim[43])+(0.276641667586244*alpha_cdim[27]+0.3872983346207416*alpha_cdim[26]+0.4330127018922193*alpha_cdim[4])*f[40]+(0.276641667586244*f[27]+0.3872983346207416*f[26]+0.4330127018922193*f[4])*alpha_cdim[40]+0.3464101615137755*(alpha_cdim[18]*f[38]+f[18]*alpha_cdim[38])+0.4330127018922193*(alpha_cdim[9]*f[27]+f[9]*alpha_cdim[27])+(0.276641667586244*alpha_cdim[13]+0.3872983346207416*alpha_cdim[12]+0.4330127018922193*alpha_cdim[0])*f[24]+(0.276641667586244*f[13]+0.3872983346207416*f[12]+0.4330127018922193*f[0])*alpha_cdim[24]+0.3464101615137755*(alpha_cdim[7]*f[22]+f[7]*alpha_cdim[22])+0.3872983346207416*(alpha_cdim[10]*f[18]+f[10]*alpha_cdim[18])+0.4330127018922193*(alpha_cdim[2]*f[13]+f[2]*alpha_cdim[13])+0.3872983346207416*(alpha_cdim[3]*f[7]+f[3]*alpha_cdim[7]); 
  out[35] += (0.7745966692414833*alpha_cdim[38]+0.8660254037844386*alpha_cdim[10])*f[47]+0.9682458365518543*alpha_cdim[13]*f[46]+0.7745966692414833*alpha_cdim[43]*f[45]+0.8660254037844386*(alpha_cdim[7]*f[45]+f[17]*alpha_cdim[43]+alpha_cdim[18]*f[42])+(0.7745966692414833*alpha_cdim[26]+0.8660254037844386*alpha_cdim[4])*f[41]+0.9682458365518543*(f[23]*alpha_cdim[40]+alpha_cdim[24]*f[39])+0.8660254037844386*f[15]*alpha_cdim[38]+(0.7745966692414833*alpha_cdim[29]+0.8660254037844386*alpha_cdim[2])*f[36]+0.9682458365518543*alpha_cdim[27]*f[34]+0.8660254037844386*alpha_cdim[18]*f[33]+(0.8660254037844386*(alpha_cdim[30]+alpha_cdim[22])+0.9682458365518543*alpha_cdim[3])*f[31]+0.8660254037844386*(f[8]*alpha_cdim[29]+alpha_cdim[9]*f[28]+f[5]*alpha_cdim[26]+alpha_cdim[9]*f[20])+0.9682458365518543*(f[6]*alpha_cdim[18]+alpha_cdim[7]*f[17])+0.8660254037844386*(alpha_cdim[14]+alpha_cdim[12])*f[16]+0.9682458365518543*(alpha_cdim[0]*f[16]+alpha_cdim[10]*f[15]+f[1]*alpha_cdim[9]+alpha_cdim[2]*f[8]+alpha_cdim[4]*f[5]); 
  out[36] += 0.3464101615137755*(alpha_cdim[18]*f[43]+f[18]*alpha_cdim[43])+0.3872983346207416*(alpha_cdim[24]*f[40]+f[24]*alpha_cdim[40])+(0.3872983346207416*alpha_cdim[30]+0.276641667586244*alpha_cdim[22]+0.4330127018922193*alpha_cdim[3])*f[38]+(0.3872983346207416*f[30]+0.276641667586244*f[22]+0.4330127018922193*f[3])*alpha_cdim[38]+0.3464101615137755*(alpha_cdim[9]*f[29]+f[9]*alpha_cdim[29])+(0.3872983346207416*alpha_cdim[14]+0.276641667586244*alpha_cdim[12]+0.4330127018922193*alpha_cdim[0])*f[26]+(0.3872983346207416*f[14]+0.276641667586244*f[12])*alpha_cdim[26]+0.4330127018922193*(f[0]*alpha_cdim[26]+alpha_cdim[10]*f[22]+f[10]*alpha_cdim[22])+0.3872983346207416*(alpha_cdim[7]*f[18]+f[7]*alpha_cdim[18])+0.4330127018922193*(alpha_cdim[4]*f[12]+f[4]*alpha_cdim[12])+0.3872983346207416*(alpha_cdim[2]*f[9]+f[2]*alpha_cdim[9]); 
  out[37] += (0.7745966692414833*alpha_cdim[40]+0.8660254037844386*alpha_cdim[9])*f[47]+(0.7745966692414833*alpha_cdim[43]+0.8660254037844386*alpha_cdim[7])*f[46]+0.9682458365518543*alpha_cdim[12]*f[45]+0.8660254037844386*f[16]*alpha_cdim[43]+0.7745966692414833*alpha_cdim[27]*f[42]+0.8660254037844386*(alpha_cdim[4]*f[42]+alpha_cdim[18]*f[41]+f[15]*alpha_cdim[40])+(0.7745966692414833*alpha_cdim[30]+0.8660254037844386*alpha_cdim[3])*f[39]+0.9682458365518543*(f[20]*alpha_cdim[38]+alpha_cdim[22]*f[36])+0.8660254037844386*alpha_cdim[18]*f[34]+0.9682458365518543*alpha_cdim[26]*f[33]+(0.8660254037844386*(alpha_cdim[29]+alpha_cdim[24])+0.9682458365518543*alpha_cdim[2])*f[31]+0.8660254037844386*(f[8]*alpha_cdim[30]+alpha_cdim[10]*f[28]+f[6]*alpha_cdim[27]+alpha_cdim[10]*f[23])+0.9682458365518543*f[5]*alpha_cdim[18]+0.8660254037844386*(alpha_cdim[14]+alpha_cdim[13])*f[17]+0.9682458365518543*(alpha_cdim[0]*f[17]+alpha_cdim[7]*f[16]+alpha_cdim[9]*f[15]+f[1]*alpha_cdim[10]+alpha_cdim[3]*f[8]+alpha_cdim[4]*f[6]); 
  out[39] += 0.3464101615137755*(alpha_cdim[18]*f[43]+f[18]*alpha_cdim[43])+(0.3872983346207416*alpha_cdim[29]+0.276641667586244*alpha_cdim[24]+0.4330127018922193*alpha_cdim[2])*f[40]+(0.3872983346207416*f[29]+0.276641667586244*f[24]+0.4330127018922193*f[2])*alpha_cdim[40]+0.3872983346207416*(alpha_cdim[22]*f[38]+f[22]*alpha_cdim[38])+0.3464101615137755*(alpha_cdim[10]*f[30]+f[10]*alpha_cdim[30])+(0.3872983346207416*alpha_cdim[14]+0.276641667586244*alpha_cdim[13]+0.4330127018922193*alpha_cdim[0])*f[27]+(0.3872983346207416*f[14]+0.276641667586244*f[13])*alpha_cdim[27]+0.4330127018922193*(f[0]*alpha_cdim[27]+alpha_cdim[9]*f[24]+f[9]*alpha_cdim[24])+0.3872983346207416*(alpha_cdim[7]*f[18]+f[7]*alpha_cdim[18])+0.4330127018922193*(alpha_cdim[4]*f[13]+f[4]*alpha_cdim[13])+0.3872983346207416*(alpha_cdim[3]*f[10]+f[3]*alpha_cdim[10]); 
  out[41] += (0.276641667586244*alpha_cdim[30]+0.3872983346207416*alpha_cdim[22]+0.4330127018922193*alpha_cdim[3])*f[43]+(0.276641667586244*f[30]+0.3872983346207416*f[22]+0.4330127018922193*f[3])*alpha_cdim[43]+0.3872983346207416*(alpha_cdim[27]*f[40]+f[27]*alpha_cdim[40])+0.3464101615137755*(alpha_cdim[18]*f[38]+f[18]*alpha_cdim[38])+0.4330127018922193*(alpha_cdim[7]*f[30]+f[7]*alpha_cdim[30])+(0.276641667586244*alpha_cdim[14]+0.3872983346207416*alpha_cdim[12]+0.4330127018922193*alpha_cdim[0])*f[29]+(0.276641667586244*f[14]+0.3872983346207416*f[12]+0.4330127018922193*f[0])*alpha_cdim[29]+0.3464101615137755*(alpha_cdim[9]*f[26]+f[9]*alpha_cdim[26])+0.3872983346207416*(alpha_cdim[10]*f[18]+f[10]*alpha_cdim[18])+0.4330127018922193*(alpha_cdim[2]*f[14]+f[2]*alpha_cdim[14])+0.3872983346207416*(alpha_cdim[4]*f[9]+f[4]*alpha_cdim[9]); 
  out[42] += (0.276641667586244*alpha_cdim[29]+0.3872983346207416*alpha_cdim[24]+0.4330127018922193*alpha_cdim[2])*f[43]+(0.276641667586244*f[29]+0.3872983346207416*f[24]+0.4330127018922193*f[2])*alpha_cdim[43]+0.3464101615137755*(alpha_cdim[18]*f[40]+f[18]*alpha_cdim[40])+0.3872983346207416*(alpha_cdim[26]*f[38]+f[26]*alpha_cdim[38])+(0.276641667586244*alpha_cdim[14]+0.3872983346207416*alpha_cdim[13]+0.4330127018922193*alpha_cdim[0])*f[30]+(0.276641667586244*f[14]+0.3872983346207416*f[13])*alpha_cdim[30]+0.4330127018922193*(f[0]*alpha_cdim[30]+alpha_cdim[7]*f[29]+f[7]*alpha_cdim[29])+0.3464101615137755*(alpha_cdim[10]*f[27]+f[10]*alpha_cdim[27])+0.3872983346207416*(alpha_cdim[9]*f[18]+f[9]*alpha_cdim[18])+0.4330127018922193*(alpha_cdim[3]*f[14]+f[3]*alpha_cdim[14])+0.3872983346207416*(alpha_cdim[4]*f[10]+f[4]*alpha_cdim[10]); 
  out[44] += (0.7745966692414833*(alpha_cdim[27]+alpha_cdim[26])+0.8660254037844386*alpha_cdim[4])*f[47]+(0.7745966692414833*(alpha_cdim[30]+alpha_cdim[22])+0.8660254037844386*alpha_cdim[3])*f[46]+(0.7745966692414833*(alpha_cdim[29]+alpha_cdim[24])+0.8660254037844386*alpha_cdim[2])*f[45]+(0.7745966692414833*(f[39]+f[36])+0.8660254037844386*f[8])*alpha_cdim[43]+(0.7745966692414833*alpha_cdim[40]+0.8660254037844386*alpha_cdim[9])*f[42]+(0.7745966692414833*alpha_cdim[38]+0.8660254037844386*alpha_cdim[10])*f[41]+0.7745966692414833*f[33]*alpha_cdim[40]+0.8660254037844386*(f[6]*alpha_cdim[40]+alpha_cdim[7]*f[39])+0.7745966692414833*f[34]*alpha_cdim[38]+0.8660254037844386*(f[5]*alpha_cdim[38]+alpha_cdim[7]*f[36]+alpha_cdim[10]*f[34]+alpha_cdim[9]*f[33])+(0.8660254037844386*(alpha_cdim[14]+alpha_cdim[13]+alpha_cdim[12])+0.9682458365518543*alpha_cdim[0])*f[31]+0.8660254037844386*(f[16]*alpha_cdim[30]+f[17]*alpha_cdim[29]+alpha_cdim[18]*f[28]+f[15]*(alpha_cdim[27]+alpha_cdim[26])+f[17]*alpha_cdim[24]+alpha_cdim[18]*f[23]+f[16]*alpha_cdim[22]+alpha_cdim[18]*f[20])+0.9682458365518543*(f[1]*alpha_cdim[18]+alpha_cdim[2]*f[17]+alpha_cdim[3]*f[16]+alpha_cdim[4]*f[15]+f[5]*alpha_cdim[10]+f[6]*alpha_cdim[9]+alpha_cdim[7]*f[8]); 
  out[45] += (0.3098386676965933*alpha_cdim[40]+0.3464101615137755*alpha_cdim[9])*f[43]+0.3098386676965933*f[40]*alpha_cdim[43]+0.3464101615137755*(f[9]*alpha_cdim[43]+alpha_cdim[7]*f[40]+f[7]*alpha_cdim[40])+(0.3872983346207416*(alpha_cdim[14]+alpha_cdim[13])+0.276641667586244*alpha_cdim[12]+0.4330127018922193*alpha_cdim[0])*f[38]+(0.3872983346207416*(f[14]+f[13])+0.276641667586244*f[12]+0.4330127018922193*f[0])*alpha_cdim[38]+0.3872983346207416*(alpha_cdim[26]*f[30]+f[26]*alpha_cdim[30])+0.3464101615137755*(alpha_cdim[18]*f[29]+f[18]*alpha_cdim[29])+0.3872983346207416*(alpha_cdim[22]*f[27]+f[22]*alpha_cdim[27])+(0.276641667586244*alpha_cdim[22]+0.4330127018922193*alpha_cdim[3])*f[26]+(0.276641667586244*f[22]+0.4330127018922193*f[3])*alpha_cdim[26]+0.3464101615137755*(alpha_cdim[18]*f[24]+f[18]*alpha_cdim[24])+0.4330127018922193*(alpha_cdim[4]*f[22]+f[4]*alpha_cdim[22])+0.3872983346207416*(alpha_cdim[2]*f[18]+f[2]*alpha_cdim[18])+0.4330127018922193*(alpha_cdim[10]*f[12]+f[10]*alpha_cdim[12])+0.3872983346207416*(alpha_cdim[7]*f[9]+f[7]*alpha_cdim[9]); 
  out[46] += (0.3098386676965933*alpha_cdim[38]+0.3464101615137755*alpha_cdim[10])*f[43]+(0.3098386676965933*f[38]+0.3464101615137755*f[10])*alpha_cdim[43]+(0.3872983346207416*alpha_cdim[14]+0.276641667586244*alpha_cdim[13]+0.3872983346207416*alpha_cdim[12]+0.4330127018922193*alpha_cdim[0])*f[40]+(0.3872983346207416*f[14]+0.276641667586244*f[13]+0.3872983346207416*f[12]+0.4330127018922193*f[0])*alpha_cdim[40]+0.3464101615137755*(alpha_cdim[7]*f[38]+f[7]*alpha_cdim[38]+alpha_cdim[18]*f[30]+f[18]*alpha_cdim[30])+0.3872983346207416*alpha_cdim[27]*f[29]+f[27]*(0.3872983346207416*alpha_cdim[29]+0.276641667586244*alpha_cdim[24]+0.4330127018922193*alpha_cdim[2])+(0.276641667586244*f[24]+0.4330127018922193*f[2])*alpha_cdim[27]+0.3872983346207416*(alpha_cdim[24]*f[26]+f[24]*alpha_cdim[26])+0.4330127018922193*(alpha_cdim[4]*f[24]+f[4]*alpha_cdim[24])+0.3464101615137755*(alpha_cdim[18]*f[22]+f[18]*alpha_cdim[22])+0.3872983346207416*(alpha_cdim[3]*f[18]+f[3]*alpha_cdim[18])+0.4330127018922193*(alpha_cdim[9]*f[13]+f[9]*alpha_cdim[13])+0.3872983346207416*(alpha_cdim[7]*f[10]+f[7]*alpha_cdim[10]); 
  out[47] += (0.276641667586244*alpha_cdim[14]+0.3872983346207416*(alpha_cdim[13]+alpha_cdim[12])+0.4330127018922193*alpha_cdim[0])*f[43]+(0.276641667586244*f[14]+0.3872983346207416*(f[13]+f[12])+0.4330127018922193*f[0])*alpha_cdim[43]+(0.3098386676965933*alpha_cdim[38]+0.3464101615137755*alpha_cdim[10])*f[40]+0.3098386676965933*f[38]*alpha_cdim[40]+0.3464101615137755*(f[10]*alpha_cdim[40]+alpha_cdim[9]*f[38]+f[9]*alpha_cdim[38])+(0.276641667586244*alpha_cdim[29]+0.3872983346207416*alpha_cdim[24]+0.4330127018922193*alpha_cdim[2])*f[30]+(0.276641667586244*f[29]+0.3872983346207416*f[24]+0.4330127018922193*f[2])*alpha_cdim[30]+(0.3872983346207416*alpha_cdim[22]+0.4330127018922193*alpha_cdim[3])*f[29]+(0.3872983346207416*f[22]+0.4330127018922193*f[3])*alpha_cdim[29]+0.3464101615137755*(alpha_cdim[18]*f[27]+f[18]*alpha_cdim[27]+alpha_cdim[18]*f[26]+f[18]*alpha_cdim[26])+0.3872983346207416*(alpha_cdim[4]*f[18]+f[4]*alpha_cdim[18])+0.4330127018922193*(alpha_cdim[7]*f[14]+f[7]*alpha_cdim[14])+0.3872983346207416*(alpha_cdim[9]*f[10]+f[9]*alpha_cdim[10]); 

  return alpha_mid; 
} 
