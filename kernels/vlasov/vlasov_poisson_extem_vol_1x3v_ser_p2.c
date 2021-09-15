#include <gkyl_vlasov_kernels.h> 

GKYL_CU_DH double vlasov_poisson_extem_vol_1x3v_ser_p2(const double *w, const double *dxv, const double *fac_phi, const double *vecA, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fac_phi:   potential (scaled by appropriate factors).
  // vecA:      vector potential (scaled by appropriate factors). Unused in pure Vlasov-Poisson. 
  // f:         Input distribution function.
  // out:       Incremented output.
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  const double dx10 = 2/dxv[0]; 
  const double *phi = &fac_phi[0]; 
  const double dv10 = 2/dxv[1]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv11 = 2/dxv[2]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double dv12 = 2/dxv[3]; 
  const double dv3 = dxv[3], wv3 = w[3]; 

  const double *A0 = &vecA[0]; 
  const double *A1 = &vecA[3]; 
  const double *A2 = &vecA[6]; 
  double alpha_mid = 0.0; 
  double alpha_cdim[48]; 
  double alpha_vdim[144]; 

  alpha_cdim[0] = 8.0*w0dx0; 
  alpha_cdim[2] = 2.309401076758503*dv0dx0; 
  alpha_mid += fabs(w0dx0)+0.5*dv0dx0; 

  alpha_vdim[0] = dv10*dx10*(4.898979485566357*(A2[1]*wv3+A1[1]*wv2)-4.898979485566357*phi[1]); 
  alpha_vdim[1] = dv10*dx10*(10.95445115010332*(A2[2]*wv3+A1[2]*wv2)-10.95445115010332*phi[2]); 
  alpha_vdim[3] = 1.414213562373095*A1[1]*dv10*dv2*dx10; 
  alpha_vdim[4] = 1.414213562373095*A2[1]*dv10*dv3*dx10; 
  alpha_vdim[6] = 3.16227766016838*A1[2]*dv10*dv2*dx10; 
  alpha_vdim[8] = 3.16227766016838*A2[2]*dv10*dv3*dx10; 
  alpha_mid += fabs(0.125*alpha_vdim[0]); 

  alpha_vdim[48] = -4.898979485566357*A1[1]*dv11*dx10*wv1; 
  alpha_vdim[49] = -10.95445115010332*A1[2]*dv11*dx10*wv1; 
  alpha_vdim[50] = -1.414213562373095*A1[1]*dv1*dv11*dx10; 
  alpha_vdim[53] = -3.16227766016838*A1[2]*dv1*dv11*dx10; 
  alpha_mid += fabs(0.125*alpha_vdim[48]); 

  alpha_vdim[96] = -4.898979485566357*A2[1]*dv12*dx10*wv1; 
  alpha_vdim[97] = -10.95445115010332*A2[2]*dv12*dx10*wv1; 
  alpha_vdim[98] = -1.414213562373095*A2[1]*dv1*dv12*dx10; 
  alpha_vdim[101] = -3.16227766016838*A2[2]*dv1*dv12*dx10; 
  alpha_mid += fabs(0.125*alpha_vdim[96]); 

  out[1] += 0.4330127018922193*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.4330127018922193*(alpha_vdim[8]*f[8]+alpha_vdim[6]*f[6]+alpha_vdim[4]*f[4]+alpha_vdim[3]*f[3]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.4330127018922193*(f[5]*alpha_vdim[53]+f[2]*alpha_vdim[50]+f[1]*alpha_vdim[49]+f[0]*alpha_vdim[48]); 
  out[4] += 0.4330127018922193*(f[5]*alpha_vdim[101]+f[2]*alpha_vdim[98]+f[1]*alpha_vdim[97]+f[0]*alpha_vdim[96]); 
  out[5] += 0.3872983346207416*(alpha_vdim[8]*f[25]+alpha_vdim[6]*f[21]+alpha_cdim[2]*f[12]+alpha_vdim[1]*f[11])+0.4330127018922193*(alpha_vdim[4]*f[8]+f[4]*alpha_vdim[8]+alpha_vdim[3]*f[6]+f[3]*alpha_vdim[6]+alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[6] += 0.3872983346207416*f[19]*alpha_vdim[53]+0.4330127018922193*(f[2]*alpha_vdim[53]+f[5]*alpha_vdim[50])+0.3872983346207416*f[11]*alpha_vdim[49]+0.4330127018922193*(f[0]*alpha_vdim[49]+f[1]*alpha_vdim[48]+alpha_cdim[2]*f[7]+alpha_cdim[0]*f[3]); 
  out[7] += (0.3872983346207416*f[20]+0.4330127018922193*f[1])*alpha_vdim[53]+0.3872983346207416*f[12]*alpha_vdim[50]+0.4330127018922193*(f[0]*alpha_vdim[50]+f[5]*alpha_vdim[49]+f[2]*alpha_vdim[48])+0.3872983346207416*alpha_vdim[6]*f[23]+0.4330127018922193*alpha_vdim[8]*f[17]+0.3872983346207416*alpha_vdim[3]*f[13]+0.4330127018922193*(alpha_vdim[4]*f[10]+alpha_vdim[1]*f[6]+f[1]*alpha_vdim[6]+alpha_vdim[0]*f[3]+f[0]*alpha_vdim[3]); 
  out[8] += 0.3872983346207416*f[19]*alpha_vdim[101]+0.4330127018922193*(f[2]*alpha_vdim[101]+f[5]*alpha_vdim[98])+0.3872983346207416*f[11]*alpha_vdim[97]+0.4330127018922193*(f[0]*alpha_vdim[97]+f[1]*alpha_vdim[96]+alpha_cdim[2]*f[9]+alpha_cdim[0]*f[4]); 
  out[9] += (0.3872983346207416*f[20]+0.4330127018922193*f[1])*alpha_vdim[101]+0.3872983346207416*f[12]*alpha_vdim[98]+0.4330127018922193*(f[0]*alpha_vdim[98]+f[5]*alpha_vdim[97]+f[2]*alpha_vdim[96])+0.3872983346207416*alpha_vdim[8]*f[28]+0.4330127018922193*alpha_vdim[6]*f[17]+0.3872983346207416*alpha_vdim[4]*f[14]+0.4330127018922193*(alpha_vdim[3]*f[10]+alpha_vdim[1]*f[8]+f[1]*alpha_vdim[8]+alpha_vdim[0]*f[4]+f[0]*alpha_vdim[4]); 
  out[10] += 0.4330127018922193*(f[15]*alpha_vdim[101]+f[7]*alpha_vdim[98]+f[6]*alpha_vdim[97]+f[3]*alpha_vdim[96]+f[16]*alpha_vdim[53]+f[9]*alpha_vdim[50]+f[8]*alpha_vdim[49]+f[4]*alpha_vdim[48]); 
  out[11] += 0.9682458365518543*(alpha_cdim[2]*f[5]+alpha_cdim[0]*f[1]); 
  out[12] += 0.9682458365518543*(alpha_vdim[8]*f[16]+alpha_vdim[6]*f[15]+alpha_vdim[4]*f[9]+alpha_vdim[3]*f[7]+alpha_vdim[1]*f[5]+alpha_vdim[0]*f[2]); 
  out[13] += 0.9682458365518543*(f[15]*alpha_vdim[53]+f[7]*alpha_vdim[50]+f[6]*alpha_vdim[49]+f[3]*alpha_vdim[48]); 
  out[14] += 0.9682458365518543*(f[16]*alpha_vdim[101]+f[9]*alpha_vdim[98]+f[8]*alpha_vdim[97]+f[4]*alpha_vdim[96]); 
  out[15] += (0.3872983346207416*(f[12]+f[11])+0.4330127018922193*f[0])*alpha_vdim[53]+(0.3872983346207416*f[20]+0.4330127018922193*f[1])*alpha_vdim[50]+0.3872983346207416*f[19]*alpha_vdim[49]+0.4330127018922193*(f[2]*alpha_vdim[49]+f[5]*alpha_vdim[48])+0.3872983346207416*(alpha_vdim[8]*f[37]+alpha_vdim[3]*f[23]+alpha_cdim[2]*f[22]+alpha_vdim[1]*f[21])+0.4330127018922193*alpha_vdim[4]*f[17]+0.3872983346207416*alpha_vdim[6]*(f[13]+f[11])+0.4330127018922193*(alpha_vdim[8]*f[10]+alpha_cdim[0]*f[7]+alpha_vdim[0]*f[6]+f[0]*alpha_vdim[6]+(alpha_cdim[2]+alpha_vdim[1])*f[3]+f[1]*alpha_vdim[3]); 
  out[16] += (0.3872983346207416*(f[12]+f[11])+0.4330127018922193*f[0])*alpha_vdim[101]+(0.3872983346207416*f[20]+0.4330127018922193*f[1])*alpha_vdim[98]+0.3872983346207416*f[19]*alpha_vdim[97]+0.4330127018922193*(f[2]*alpha_vdim[97]+f[5]*alpha_vdim[96])+0.3872983346207416*(alpha_vdim[6]*f[37]+alpha_vdim[4]*f[28]+alpha_cdim[2]*f[26]+alpha_vdim[1]*f[25])+0.4330127018922193*alpha_vdim[3]*f[17]+0.3872983346207416*alpha_vdim[8]*(f[14]+f[11])+0.4330127018922193*(alpha_vdim[6]*f[10]+alpha_cdim[0]*f[9]+alpha_vdim[0]*f[8]+f[0]*alpha_vdim[8]+(alpha_cdim[2]+alpha_vdim[1])*f[4]+f[1]*alpha_vdim[4]); 
  out[17] += 0.3872983346207416*f[32]*alpha_vdim[101]+0.4330127018922193*(f[7]*alpha_vdim[101]+f[15]*alpha_vdim[98])+0.3872983346207416*f[21]*alpha_vdim[97]+0.4330127018922193*(f[3]*alpha_vdim[97]+f[6]*alpha_vdim[96])+0.3872983346207416*f[35]*alpha_vdim[53]+0.4330127018922193*(f[9]*alpha_vdim[53]+f[16]*alpha_vdim[50])+0.3872983346207416*f[25]*alpha_vdim[49]+0.4330127018922193*(f[4]*alpha_vdim[49]+f[8]*alpha_vdim[48]+alpha_cdim[2]*f[18]+alpha_cdim[0]*f[10]); 
  out[18] += (0.3872983346207416*f[33]+0.4330127018922193*f[6])*alpha_vdim[101]+0.3872983346207416*f[22]*alpha_vdim[98]+0.4330127018922193*(f[3]*alpha_vdim[98]+f[15]*alpha_vdim[97]+f[7]*alpha_vdim[96])+(0.3872983346207416*f[36]+0.4330127018922193*f[8])*alpha_vdim[53]+0.3872983346207416*f[26]*alpha_vdim[50]+0.4330127018922193*(f[4]*alpha_vdim[50]+f[16]*alpha_vdim[49]+f[9]*alpha_vdim[48])+0.3872983346207416*(alpha_vdim[8]*f[42]+alpha_vdim[6]*f[39]+alpha_vdim[4]*f[30]+alpha_vdim[3]*f[27])+0.4330127018922193*(alpha_vdim[1]*f[17]+alpha_vdim[0]*f[10]+alpha_vdim[6]*f[8]+f[6]*alpha_vdim[8]+alpha_vdim[3]*f[4]+f[3]*alpha_vdim[4]); 
  out[19] += 0.4330127018922193*(alpha_vdim[4]*f[25]+alpha_vdim[3]*f[21])+0.8660254037844386*alpha_cdim[2]*f[20]+0.4330127018922193*alpha_vdim[0]*f[11]+0.3872983346207416*(alpha_vdim[8]*f[8]+alpha_vdim[6]*f[6])+0.9682458365518543*alpha_cdim[0]*f[5]+f[1]*(0.9682458365518543*alpha_cdim[2]+0.3872983346207416*alpha_vdim[1]); 
  out[20] += 0.8660254037844386*(alpha_vdim[8]*f[35]+alpha_vdim[6]*f[32]+alpha_vdim[1]*f[19])+0.9682458365518543*(alpha_vdim[4]*f[16]+alpha_vdim[3]*f[15])+0.4330127018922193*alpha_cdim[0]*f[12]+0.9682458365518543*(alpha_vdim[8]*f[9]+alpha_vdim[6]*f[7]+alpha_vdim[0]*f[5])+(0.3872983346207416*alpha_cdim[2]+0.9682458365518543*alpha_vdim[1])*f[2]; 
  out[21] += 0.3872983346207416*f[5]*alpha_vdim[53]+0.4330127018922193*f[19]*alpha_vdim[50]+0.3872983346207416*f[1]*alpha_vdim[49]+0.4330127018922193*f[11]*alpha_vdim[48]+0.9682458365518543*(alpha_cdim[2]*f[15]+alpha_cdim[0]*f[6]); 
  out[22] += 0.3872983346207416*(f[5]*alpha_vdim[53]+f[2]*alpha_vdim[50])+0.4330127018922193*(f[20]*alpha_vdim[49]+f[12]*alpha_vdim[48])+0.8660254037844386*alpha_vdim[6]*f[34]+0.9682458365518543*alpha_vdim[8]*f[31]+0.8660254037844386*alpha_vdim[3]*f[24]+0.9682458365518543*(alpha_vdim[4]*f[18]+alpha_vdim[1]*f[15]+alpha_vdim[0]*f[7]+f[5]*alpha_vdim[6]+f[2]*alpha_vdim[3]); 
  out[23] += 0.8660254037844386*f[32]*alpha_vdim[53]+0.9682458365518543*(f[7]*alpha_vdim[53]+f[15]*alpha_vdim[50])+0.8660254037844386*f[21]*alpha_vdim[49]+0.9682458365518543*(f[3]*alpha_vdim[49]+f[6]*alpha_vdim[48])+0.4330127018922193*(alpha_cdim[2]*f[24]+alpha_cdim[0]*f[13]); 
  out[24] += (0.8660254037844386*f[33]+0.9682458365518543*f[6])*alpha_vdim[53]+0.8660254037844386*f[22]*alpha_vdim[50]+0.9682458365518543*(f[3]*alpha_vdim[50]+f[15]*alpha_vdim[49]+f[7]*alpha_vdim[48])+0.4330127018922193*(alpha_vdim[8]*f[39]+alpha_vdim[4]*f[27]+alpha_vdim[1]*f[23]+alpha_vdim[0]*f[13])+0.3872983346207416*(alpha_vdim[6]*f[6]+alpha_vdim[3]*f[3]); 
  out[25] += 0.3872983346207416*f[5]*alpha_vdim[101]+0.4330127018922193*f[19]*alpha_vdim[98]+0.3872983346207416*f[1]*alpha_vdim[97]+0.4330127018922193*f[11]*alpha_vdim[96]+0.9682458365518543*(alpha_cdim[2]*f[16]+alpha_cdim[0]*f[8]); 
  out[26] += 0.3872983346207416*(f[5]*alpha_vdim[101]+f[2]*alpha_vdim[98])+0.4330127018922193*(f[20]*alpha_vdim[97]+f[12]*alpha_vdim[96])+0.8660254037844386*alpha_vdim[8]*f[41]+0.9682458365518543*alpha_vdim[6]*f[31]+0.8660254037844386*alpha_vdim[4]*f[29]+0.9682458365518543*(alpha_vdim[3]*f[18]+alpha_vdim[1]*f[16]+alpha_vdim[0]*f[9]+f[5]*alpha_vdim[8]+f[2]*alpha_vdim[4]); 
  out[27] += 0.4330127018922193*(f[34]*alpha_vdim[101]+f[24]*alpha_vdim[98]+f[23]*alpha_vdim[97]+f[13]*alpha_vdim[96])+0.9682458365518543*(f[31]*alpha_vdim[53]+f[18]*alpha_vdim[50]+f[17]*alpha_vdim[49]+f[10]*alpha_vdim[48]); 
  out[28] += 0.8660254037844386*f[35]*alpha_vdim[101]+0.9682458365518543*(f[9]*alpha_vdim[101]+f[16]*alpha_vdim[98])+0.8660254037844386*f[25]*alpha_vdim[97]+0.9682458365518543*(f[4]*alpha_vdim[97]+f[8]*alpha_vdim[96])+0.4330127018922193*(alpha_cdim[2]*f[29]+alpha_cdim[0]*f[14]); 
  out[29] += (0.8660254037844386*f[36]+0.9682458365518543*f[8])*alpha_vdim[101]+0.8660254037844386*f[26]*alpha_vdim[98]+0.9682458365518543*(f[4]*alpha_vdim[98]+f[16]*alpha_vdim[97]+f[9]*alpha_vdim[96])+0.4330127018922193*(alpha_vdim[6]*f[42]+alpha_vdim[3]*f[30]+alpha_vdim[1]*f[28]+alpha_vdim[0]*f[14])+0.3872983346207416*(alpha_vdim[8]*f[8]+alpha_vdim[4]*f[4]); 
  out[30] += 0.9682458365518543*(f[31]*alpha_vdim[101]+f[18]*alpha_vdim[98]+f[17]*alpha_vdim[97]+f[10]*alpha_vdim[96])+0.4330127018922193*(f[41]*alpha_vdim[53]+f[29]*alpha_vdim[50]+f[28]*alpha_vdim[49]+f[14]*alpha_vdim[48]); 
  out[31] += (0.3872983346207416*(f[22]+f[21])+0.4330127018922193*f[3])*alpha_vdim[101]+(0.3872983346207416*f[33]+0.4330127018922193*f[6])*alpha_vdim[98]+0.3872983346207416*f[32]*alpha_vdim[97]+0.4330127018922193*(f[7]*alpha_vdim[97]+f[15]*alpha_vdim[96])+(0.3872983346207416*(f[26]+f[25])+0.4330127018922193*f[4])*alpha_vdim[53]+(0.3872983346207416*f[36]+0.4330127018922193*f[8])*alpha_vdim[50]+0.3872983346207416*f[35]*alpha_vdim[49]+0.4330127018922193*(f[9]*alpha_vdim[49]+f[16]*alpha_vdim[48])+0.3872983346207416*(alpha_vdim[4]*f[42]+alpha_vdim[3]*f[39]+alpha_cdim[2]*f[38]+alpha_vdim[1]*f[37]+alpha_vdim[8]*f[30]+alpha_vdim[6]*(f[27]+f[25])+alpha_vdim[8]*f[21])+0.4330127018922193*(alpha_cdim[0]*f[18]+alpha_vdim[0]*f[17]+(alpha_cdim[2]+alpha_vdim[1])*f[10]+alpha_vdim[3]*f[8]+f[3]*alpha_vdim[8]+alpha_vdim[4]*f[6]+f[4]*alpha_vdim[6]); 
  out[32] += (0.3464101615137755*f[20]+0.3872983346207416*f[1])*alpha_vdim[53]+0.4330127018922193*f[11]*alpha_vdim[50]+0.3872983346207416*f[5]*alpha_vdim[49]+0.4330127018922193*(f[19]*alpha_vdim[48]+alpha_vdim[4]*f[37])+0.8660254037844386*alpha_cdim[2]*f[33]+0.3464101615137755*alpha_vdim[6]*f[23]+0.4330127018922193*alpha_vdim[0]*f[21]+0.3872983346207416*alpha_vdim[8]*f[17]+0.9682458365518543*alpha_cdim[0]*f[15]+0.4330127018922193*alpha_vdim[3]*f[11]+0.9682458365518543*alpha_cdim[2]*f[6]+0.3872983346207416*(alpha_vdim[1]*f[6]+f[1]*alpha_vdim[6]); 
  out[33] += 0.3464101615137755*f[19]*alpha_vdim[53]+0.3872983346207416*(f[2]*alpha_vdim[53]+f[5]*alpha_vdim[50])+0.4330127018922193*(f[12]*alpha_vdim[49]+f[20]*alpha_vdim[48])+0.8660254037844386*(alpha_vdim[8]*f[44]+alpha_vdim[3]*f[34]+alpha_vdim[1]*f[32])+0.9682458365518543*alpha_vdim[4]*f[31]+0.8660254037844386*alpha_vdim[6]*f[24]+0.4330127018922193*alpha_cdim[0]*f[22]+0.8660254037844386*alpha_vdim[6]*f[19]+0.9682458365518543*(alpha_vdim[8]*f[18]+alpha_vdim[0]*f[15])+0.3872983346207416*alpha_cdim[2]*f[7]+0.9682458365518543*(alpha_vdim[1]*f[7]+f[2]*alpha_vdim[6]+alpha_vdim[3]*f[5]); 
  out[34] += (0.8660254037844386*(f[22]+f[21])+0.9682458365518543*f[3])*alpha_vdim[53]+(0.8660254037844386*f[33]+0.9682458365518543*f[6])*alpha_vdim[50]+0.8660254037844386*f[32]*alpha_vdim[49]+0.9682458365518543*(f[7]*alpha_vdim[49]+f[15]*alpha_vdim[48])+0.4330127018922193*(alpha_vdim[4]*f[39]+alpha_vdim[8]*f[27]+alpha_cdim[0]*f[24]+alpha_vdim[0]*f[23])+0.3464101615137755*alpha_vdim[6]*f[21]+0.4330127018922193*(alpha_cdim[2]+alpha_vdim[1])*f[13]+0.3872983346207416*(alpha_vdim[3]*f[6]+f[3]*alpha_vdim[6]); 
  out[35] += (0.3464101615137755*f[20]+0.3872983346207416*f[1])*alpha_vdim[101]+0.4330127018922193*f[11]*alpha_vdim[98]+0.3872983346207416*f[5]*alpha_vdim[97]+0.4330127018922193*(f[19]*alpha_vdim[96]+alpha_vdim[3]*f[37])+0.8660254037844386*alpha_cdim[2]*f[36]+0.3464101615137755*alpha_vdim[8]*f[28]+0.4330127018922193*alpha_vdim[0]*f[25]+0.3872983346207416*alpha_vdim[6]*f[17]+0.9682458365518543*alpha_cdim[0]*f[16]+0.4330127018922193*alpha_vdim[4]*f[11]+0.9682458365518543*alpha_cdim[2]*f[8]+0.3872983346207416*(alpha_vdim[1]*f[8]+f[1]*alpha_vdim[8]); 
  out[36] += 0.3464101615137755*f[19]*alpha_vdim[101]+0.3872983346207416*(f[2]*alpha_vdim[101]+f[5]*alpha_vdim[98])+0.4330127018922193*(f[12]*alpha_vdim[97]+f[20]*alpha_vdim[96])+0.8660254037844386*(alpha_vdim[6]*f[44]+alpha_vdim[4]*f[41]+alpha_vdim[1]*f[35])+0.9682458365518543*alpha_vdim[3]*f[31]+0.8660254037844386*alpha_vdim[8]*f[29]+0.4330127018922193*alpha_cdim[0]*f[26]+0.8660254037844386*alpha_vdim[8]*f[19]+0.9682458365518543*(alpha_vdim[6]*f[18]+alpha_vdim[0]*f[16])+0.3872983346207416*alpha_cdim[2]*f[9]+0.9682458365518543*(alpha_vdim[1]*f[9]+f[2]*alpha_vdim[8]+alpha_vdim[4]*f[5]); 
  out[37] += 0.3872983346207416*f[15]*alpha_vdim[101]+0.4330127018922193*f[32]*alpha_vdim[98]+0.3872983346207416*f[6]*alpha_vdim[97]+0.4330127018922193*f[21]*alpha_vdim[96]+0.3872983346207416*f[16]*alpha_vdim[53]+0.4330127018922193*f[35]*alpha_vdim[50]+0.3872983346207416*f[8]*alpha_vdim[49]+0.4330127018922193*f[25]*alpha_vdim[48]+0.9682458365518543*(alpha_cdim[2]*f[31]+alpha_cdim[0]*f[17]); 
  out[38] += 0.3872983346207416*(f[15]*alpha_vdim[101]+f[7]*alpha_vdim[98])+0.4330127018922193*(f[33]*alpha_vdim[97]+f[22]*alpha_vdim[96])+0.3872983346207416*(f[16]*alpha_vdim[53]+f[9]*alpha_vdim[50])+0.4330127018922193*(f[36]*alpha_vdim[49]+f[26]*alpha_vdim[48])+0.8660254037844386*(alpha_vdim[8]*f[47]+alpha_vdim[6]*f[46]+alpha_vdim[4]*f[43]+alpha_vdim[3]*f[40])+0.9682458365518543*(alpha_vdim[1]*f[31]+alpha_vdim[0]*f[18]+alpha_vdim[6]*f[16]+alpha_vdim[8]*f[15]+alpha_vdim[3]*f[9]+alpha_vdim[4]*f[7]); 
  out[39] += 0.4330127018922193*(f[24]*alpha_vdim[101]+f[34]*alpha_vdim[98]+f[13]*alpha_vdim[97]+f[23]*alpha_vdim[96])+0.8660254037844386*f[44]*alpha_vdim[53]+0.9682458365518543*(f[18]*alpha_vdim[53]+f[31]*alpha_vdim[50])+0.8660254037844386*f[37]*alpha_vdim[49]+0.9682458365518543*(f[10]*alpha_vdim[49]+f[17]*alpha_vdim[48])+0.4330127018922193*(alpha_cdim[2]*f[40]+alpha_cdim[0]*f[27]); 
  out[40] += 0.4330127018922193*(f[23]*alpha_vdim[101]+f[13]*alpha_vdim[98]+f[34]*alpha_vdim[97]+f[24]*alpha_vdim[96])+(0.8660254037844386*f[45]+0.9682458365518543*f[17])*alpha_vdim[53]+0.8660254037844386*f[38]*alpha_vdim[50]+0.9682458365518543*(f[10]*alpha_vdim[50]+f[31]*alpha_vdim[49]+f[18]*alpha_vdim[48])+0.4330127018922193*(alpha_vdim[1]*f[39]+alpha_vdim[0]*f[27]+alpha_vdim[8]*f[23])+0.3872983346207416*alpha_vdim[6]*f[17]+0.4330127018922193*alpha_vdim[4]*f[13]+0.3872983346207416*alpha_vdim[3]*f[10]; 
  out[41] += (0.8660254037844386*(f[26]+f[25])+0.9682458365518543*f[4])*alpha_vdim[101]+(0.8660254037844386*f[36]+0.9682458365518543*f[8])*alpha_vdim[98]+0.8660254037844386*f[35]*alpha_vdim[97]+0.9682458365518543*(f[9]*alpha_vdim[97]+f[16]*alpha_vdim[96])+0.4330127018922193*(alpha_vdim[3]*f[42]+alpha_vdim[6]*f[30]+alpha_cdim[0]*f[29]+alpha_vdim[0]*f[28])+0.3464101615137755*alpha_vdim[8]*f[25]+0.4330127018922193*(alpha_cdim[2]+alpha_vdim[1])*f[14]+0.3872983346207416*(alpha_vdim[4]*f[8]+f[4]*alpha_vdim[8]); 
  out[42] += 0.8660254037844386*f[44]*alpha_vdim[101]+0.9682458365518543*(f[18]*alpha_vdim[101]+f[31]*alpha_vdim[98])+0.8660254037844386*f[37]*alpha_vdim[97]+0.9682458365518543*(f[10]*alpha_vdim[97]+f[17]*alpha_vdim[96])+0.4330127018922193*(f[29]*alpha_vdim[53]+f[41]*alpha_vdim[50]+f[14]*alpha_vdim[49]+f[28]*alpha_vdim[48]+alpha_cdim[2]*f[43]+alpha_cdim[0]*f[30]); 
  out[43] += (0.8660254037844386*f[45]+0.9682458365518543*f[17])*alpha_vdim[101]+0.8660254037844386*f[38]*alpha_vdim[98]+0.9682458365518543*(f[10]*alpha_vdim[98]+f[31]*alpha_vdim[97]+f[18]*alpha_vdim[96])+0.4330127018922193*(f[28]*alpha_vdim[53]+f[14]*alpha_vdim[50]+f[41]*alpha_vdim[49]+f[29]*alpha_vdim[48]+alpha_vdim[1]*f[42]+alpha_vdim[0]*f[30]+alpha_vdim[6]*f[28])+0.3872983346207416*alpha_vdim[8]*f[17]+0.4330127018922193*alpha_vdim[3]*f[14]+0.3872983346207416*alpha_vdim[4]*f[10]; 
  out[44] += (0.3464101615137755*f[33]+0.3872983346207416*f[6])*alpha_vdim[101]+0.4330127018922193*f[21]*alpha_vdim[98]+0.3872983346207416*f[15]*alpha_vdim[97]+0.4330127018922193*f[32]*alpha_vdim[96]+(0.3464101615137755*f[36]+0.3872983346207416*f[8])*alpha_vdim[53]+0.4330127018922193*f[25]*alpha_vdim[50]+0.3872983346207416*f[16]*alpha_vdim[49]+0.4330127018922193*f[35]*alpha_vdim[48]+0.8660254037844386*alpha_cdim[2]*f[45]+0.3464101615137755*(alpha_vdim[8]*f[42]+alpha_vdim[6]*f[39])+0.4330127018922193*alpha_vdim[0]*f[37]+0.9682458365518543*alpha_cdim[0]*f[31]+0.4330127018922193*(alpha_vdim[3]*f[25]+alpha_vdim[4]*f[21])+0.9682458365518543*alpha_cdim[2]*f[17]+0.3872983346207416*(alpha_vdim[1]*f[17]+alpha_vdim[6]*f[8]+f[6]*alpha_vdim[8]); 
  out[45] += 0.3464101615137755*f[32]*alpha_vdim[101]+0.3872983346207416*(f[7]*alpha_vdim[101]+f[15]*alpha_vdim[98])+0.4330127018922193*(f[22]*alpha_vdim[97]+f[33]*alpha_vdim[96])+0.3464101615137755*f[35]*alpha_vdim[53]+0.3872983346207416*(f[9]*alpha_vdim[53]+f[16]*alpha_vdim[50])+0.4330127018922193*(f[26]*alpha_vdim[49]+f[36]*alpha_vdim[48])+0.8660254037844386*(alpha_vdim[4]*f[47]+alpha_vdim[3]*f[46]+alpha_vdim[1]*f[44]+alpha_vdim[8]*f[43]+alpha_vdim[6]*f[40])+0.4330127018922193*alpha_cdim[0]*f[38]+0.8660254037844386*(alpha_vdim[6]*f[35]+alpha_vdim[8]*f[32])+0.9682458365518543*alpha_vdim[0]*f[31]+0.3872983346207416*alpha_cdim[2]*f[18]+0.9682458365518543*(alpha_vdim[1]*f[18]+alpha_vdim[3]*f[16]+alpha_vdim[4]*f[15]+alpha_vdim[6]*f[9]+f[7]*alpha_vdim[8]); 
  out[46] += 0.4330127018922193*(f[13]*alpha_vdim[101]+f[23]*alpha_vdim[98]+f[24]*alpha_vdim[97]+f[34]*alpha_vdim[96])+(0.8660254037844386*(f[38]+f[37])+0.9682458365518543*f[10])*alpha_vdim[53]+(0.8660254037844386*f[45]+0.9682458365518543*f[17])*alpha_vdim[50]+0.8660254037844386*f[44]*alpha_vdim[49]+0.9682458365518543*(f[18]*alpha_vdim[49]+f[31]*alpha_vdim[48])+0.4330127018922193*(alpha_cdim[0]*f[40]+alpha_vdim[0]*f[39])+0.3464101615137755*alpha_vdim[6]*f[37]+0.4330127018922193*((alpha_cdim[2]+alpha_vdim[1])*f[27]+alpha_vdim[4]*f[23])+0.3872983346207416*alpha_vdim[3]*f[17]+0.4330127018922193*alpha_vdim[8]*f[13]+0.3872983346207416*alpha_vdim[6]*f[10]; 
  out[47] += (0.8660254037844386*(f[38]+f[37])+0.9682458365518543*f[10])*alpha_vdim[101]+(0.8660254037844386*f[45]+0.9682458365518543*f[17])*alpha_vdim[98]+0.8660254037844386*f[44]*alpha_vdim[97]+0.9682458365518543*(f[18]*alpha_vdim[97]+f[31]*alpha_vdim[96])+0.4330127018922193*(f[14]*alpha_vdim[53]+f[28]*alpha_vdim[50]+f[29]*alpha_vdim[49]+f[41]*alpha_vdim[48]+alpha_cdim[0]*f[43]+alpha_vdim[0]*f[42])+0.3464101615137755*alpha_vdim[8]*f[37]+0.4330127018922193*((alpha_cdim[2]+alpha_vdim[1])*f[30]+alpha_vdim[3]*f[28])+0.3872983346207416*alpha_vdim[4]*f[17]+0.4330127018922193*alpha_vdim[6]*f[14]+0.3872983346207416*alpha_vdim[8]*f[10]; 

  return alpha_mid; 
} 

