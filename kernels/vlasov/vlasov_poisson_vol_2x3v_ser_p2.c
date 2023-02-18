#include <gkyl_vlasov_kernels.h> 

GKYL_CU_DH double vlasov_poisson_vol_2x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // field:     potential (scaled by appropriate factors).
  // f:         Input distribution function.
  // out:       Incremented output.
  double dv0dx0 = dxv[2]/dxv[0]; 
  double w0dx0 = w[2]/dxv[0]; 
  const double dx10 = 2/dxv[0]; 
  double dv1dx1 = dxv[3]/dxv[1]; 
  double w1dx1 = w[3]/dxv[1]; 
  const double dx11 = 2/dxv[1]; 
  const double *phi = &field[0]; 
  const double dv10 = 2/dxv[2]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv11 = 2/dxv[3]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double dv12 = 2/dxv[4]; 
  const double dv3 = dxv[4], wv3 = w[4]; 

  double cflFreq_mid = 0.0; 
  double alpha_cdim[224]; 
  double alpha_vdim[336]; 

  alpha_cdim[0] = 11.31370849898477*w0dx0; 
  alpha_cdim[3] = 3.265986323710906*dv0dx0; 
  cflFreq_mid += 5.0*(fabs(w0dx0)+0.5*dv0dx0); 

  alpha_cdim[112] = 11.31370849898477*w1dx1; 
  alpha_cdim[116] = 3.265986323710906*dv1dx1; 
  cflFreq_mid += 5.0*(fabs(w1dx1)+0.5*dv1dx1); 

  alpha_vdim[0] = -4.898979485566357*phi[1]*dv10*dx10; 
  alpha_vdim[1] = -10.95445115010332*phi[4]*dv10*dx10; 
  alpha_vdim[2] = -4.898979485566357*phi[3]*dv10*dx10; 
  alpha_vdim[6] = -10.95445115010332*phi[6]*dv10*dx10; 
  alpha_vdim[17] = -4.898979485566357*phi[7]*dv10*dx10; 
  cflFreq_mid += 5.0*fabs(0.0883883476483184*alpha_vdim[0]-0.09882117688026182*alpha_vdim[17]); 

  alpha_vdim[112] = -4.898979485566357*phi[2]*dv11*dx11; 
  alpha_vdim[113] = -4.898979485566357*phi[3]*dv11*dx11; 
  alpha_vdim[114] = -10.95445115010332*phi[5]*dv11*dx11; 
  alpha_vdim[118] = -10.95445115010332*phi[7]*dv11*dx11; 
  alpha_vdim[128] = -4.898979485566357*phi[6]*dv11*dx11; 
  cflFreq_mid += 5.0*fabs(0.0883883476483184*alpha_vdim[112]-0.09882117688026182*alpha_vdim[128]); 

  cflFreq_mid += 5.0*fabs(0.0); 

  out[1] += 0.3061862178478971*(alpha_cdim[3]*f[3]+alpha_cdim[0]*f[0]); 
  out[2] += 0.3061862178478971*(f[4]*alpha_cdim[116]+f[0]*alpha_cdim[112]); 
  out[3] += 0.3061862178478971*(alpha_vdim[17]*f[17]+alpha_vdim[6]*f[6]+alpha_vdim[2]*f[2]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[4] += 0.3061862178478971*(f[16]*alpha_vdim[128]+f[6]*alpha_vdim[118]+f[2]*alpha_vdim[114]+f[1]*alpha_vdim[113]+f[0]*alpha_vdim[112]); 
  out[6] += 0.3061862178478971*(f[9]*alpha_cdim[116]+f[1]*alpha_cdim[112]+alpha_cdim[3]*f[8]+alpha_cdim[0]*f[2]); 
  out[7] += 0.3061862178478971*alpha_vdim[17]*f[32]+0.273861278752583*(alpha_vdim[6]*f[31]+alpha_cdim[3]*f[18]+alpha_vdim[1]*f[16])+0.3061862178478971*(alpha_vdim[2]*f[6]+f[2]*alpha_vdim[6]+alpha_cdim[0]*f[3]+f[0]*alpha_cdim[3]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[8] += 0.3061862178478971*(f[11]*alpha_cdim[116]+f[3]*alpha_cdim[112])+0.273861278752583*(alpha_vdim[6]*f[32]+alpha_vdim[2]*f[17]+f[2]*alpha_vdim[17])+0.3061862178478971*(alpha_vdim[1]*f[6]+f[1]*alpha_vdim[6]+alpha_vdim[0]*f[2]+f[0]*alpha_vdim[2]); 
  out[9] += 0.273861278752583*(f[1]*alpha_vdim[128]+f[31]*alpha_vdim[118])+0.3061862178478971*(f[2]*alpha_vdim[118]+f[6]*alpha_vdim[114])+0.273861278752583*f[16]*alpha_vdim[113]+0.3061862178478971*(f[0]*alpha_vdim[113]+f[1]*alpha_vdim[112]+alpha_cdim[3]*f[11]+alpha_cdim[0]*f[4]); 
  out[10] += 0.3061862178478971*f[31]*alpha_vdim[128]+(0.273861278752583*f[32]+0.3061862178478971*f[1])*alpha_vdim[118]+(0.273861278752583*f[19]+0.3061862178478971*f[0])*alpha_cdim[116]+0.273861278752583*f[17]*alpha_vdim[114]+0.3061862178478971*(f[0]*alpha_vdim[114]+f[6]*alpha_vdim[113]+f[2]*alpha_vdim[112]+f[4]*alpha_cdim[112]); 
  out[11] += 0.3061862178478971*(f[33]*alpha_vdim[128]+f[21]*alpha_vdim[118]+f[8]*alpha_vdim[114]+f[7]*alpha_vdim[113]+f[3]*alpha_vdim[112]+alpha_vdim[17]*f[38]+alpha_vdim[6]*f[22]+alpha_vdim[2]*f[10]+alpha_vdim[1]*f[9]+alpha_vdim[0]*f[4]); 
  out[12] += 0.3061862178478971*(alpha_cdim[3]*f[14]+alpha_cdim[0]*f[5]); 
  out[13] += 0.3061862178478971*(f[15]*alpha_cdim[116]+f[5]*alpha_cdim[112]); 
  out[14] += 0.3061862178478971*(alpha_vdim[17]*f[44]+alpha_vdim[6]*f[25]+alpha_vdim[2]*f[13]+alpha_vdim[1]*f[12]+alpha_vdim[0]*f[5]); 
  out[15] += 0.3061862178478971*(f[43]*alpha_vdim[128]+f[25]*alpha_vdim[118]+f[13]*alpha_vdim[114]+f[12]*alpha_vdim[113]+f[5]*alpha_vdim[112]); 
  out[16] += 0.6846531968814573*(alpha_cdim[3]*f[7]+alpha_cdim[0]*f[1]); 
  out[17] += 0.6846531968814573*(f[10]*alpha_cdim[116]+f[2]*alpha_cdim[112]); 
  out[18] += 0.6846531968814573*(alpha_vdim[17]*f[34]+alpha_vdim[6]*f[21]+alpha_vdim[2]*f[8]+alpha_vdim[1]*f[7]+alpha_vdim[0]*f[3]); 
  out[19] += 0.6846531968814573*(f[37]*alpha_vdim[128]+f[22]*alpha_vdim[118]+f[10]*alpha_vdim[114]+f[9]*alpha_vdim[113]+f[4]*alpha_vdim[112]); 
  out[21] += 0.3061862178478971*(f[23]*alpha_cdim[116]+f[7]*alpha_cdim[112])+0.273861278752583*(alpha_cdim[3]*f[36]+alpha_vdim[2]*f[32]+alpha_vdim[1]*f[31]+alpha_vdim[6]*f[17]+f[6]*alpha_vdim[17]+alpha_vdim[6]*f[16])+0.3061862178478971*(alpha_cdim[0]*f[8]+alpha_vdim[0]*f[6]+f[0]*alpha_vdim[6]+f[2]*(alpha_cdim[3]+alpha_vdim[1])+f[1]*alpha_vdim[2]); 
  out[22] += 0.273861278752583*f[6]*alpha_vdim[128]+(0.273861278752583*(f[17]+f[16])+0.3061862178478971*f[0])*alpha_vdim[118]+(0.273861278752583*f[40]+0.3061862178478971*f[1])*alpha_cdim[116]+(0.273861278752583*f[32]+0.3061862178478971*f[1])*alpha_vdim[114]+0.273861278752583*f[31]*alpha_vdim[113]+0.3061862178478971*(f[2]*alpha_vdim[113]+f[6]*alpha_vdim[112]+f[9]*alpha_cdim[112]+alpha_cdim[3]*f[24]+alpha_cdim[0]*f[10]); 
  out[23] += 0.273861278752583*(f[7]*alpha_vdim[128]+f[56]*alpha_vdim[118])+0.3061862178478971*(f[8]*alpha_vdim[118]+f[21]*alpha_vdim[114])+0.273861278752583*f[33]*alpha_vdim[113]+0.3061862178478971*(f[3]*alpha_vdim[113]+f[7]*alpha_vdim[112]+alpha_vdim[17]*f[60])+0.273861278752583*(alpha_vdim[6]*f[59]+alpha_cdim[3]*f[39]+alpha_vdim[1]*f[37])+0.3061862178478971*(alpha_vdim[2]*f[22]+alpha_cdim[0]*f[11]+alpha_vdim[6]*f[10]+alpha_vdim[0]*f[9]+(alpha_cdim[3]+alpha_vdim[1])*f[4]); 
  out[24] += 0.3061862178478971*f[56]*alpha_vdim[128]+(0.273861278752583*f[57]+0.3061862178478971*f[7])*alpha_vdim[118]+(0.273861278752583*f[42]+0.3061862178478971*f[3])*alpha_cdim[116]+0.273861278752583*f[34]*alpha_vdim[114]+0.3061862178478971*(f[3]*alpha_vdim[114]+f[21]*alpha_vdim[113]+f[8]*alpha_vdim[112]+f[11]*alpha_cdim[112])+0.273861278752583*(alpha_vdim[6]*f[60]+alpha_vdim[2]*f[38])+0.3061862178478971*alpha_vdim[1]*f[22]+0.273861278752583*f[10]*alpha_vdim[17]+0.3061862178478971*(alpha_vdim[0]*f[10]+alpha_vdim[6]*f[9]+alpha_vdim[2]*f[4]); 
  out[25] += 0.3061862178478971*(f[28]*alpha_cdim[116]+f[12]*alpha_cdim[112]+alpha_cdim[3]*f[27]+alpha_cdim[0]*f[13]); 
  out[26] += 0.3061862178478971*alpha_vdim[17]*f[69]+0.273861278752583*(alpha_vdim[6]*f[68]+alpha_cdim[3]*f[45]+alpha_vdim[1]*f[43])+0.3061862178478971*(alpha_vdim[2]*f[25]+alpha_cdim[0]*f[14]+alpha_vdim[6]*f[13]+alpha_vdim[0]*f[12]+(alpha_cdim[3]+alpha_vdim[1])*f[5]); 
  out[27] += 0.3061862178478971*(f[30]*alpha_cdim[116]+f[14]*alpha_cdim[112])+0.273861278752583*(alpha_vdim[6]*f[69]+alpha_vdim[2]*f[44])+0.3061862178478971*alpha_vdim[1]*f[25]+0.273861278752583*f[13]*alpha_vdim[17]+0.3061862178478971*(alpha_vdim[0]*f[13]+alpha_vdim[6]*f[12]+alpha_vdim[2]*f[5]); 
  out[28] += 0.273861278752583*(f[12]*alpha_vdim[128]+f[68]*alpha_vdim[118])+0.3061862178478971*(f[13]*alpha_vdim[118]+f[25]*alpha_vdim[114])+0.273861278752583*f[43]*alpha_vdim[113]+0.3061862178478971*(f[5]*alpha_vdim[113]+f[12]*alpha_vdim[112]+alpha_cdim[3]*f[30]+alpha_cdim[0]*f[15]); 
  out[29] += 0.3061862178478971*f[68]*alpha_vdim[128]+(0.273861278752583*f[69]+0.3061862178478971*f[12])*alpha_vdim[118]+(0.273861278752583*f[46]+0.3061862178478971*f[5])*alpha_cdim[116]+0.273861278752583*f[44]*alpha_vdim[114]+0.3061862178478971*(f[5]*alpha_vdim[114]+f[25]*alpha_vdim[113]+f[13]*alpha_vdim[112]+f[15]*alpha_cdim[112]); 
  out[30] += 0.3061862178478971*(f[70]*alpha_vdim[128]+f[52]*alpha_vdim[118]+f[27]*alpha_vdim[114]+f[26]*alpha_vdim[113]+f[14]*alpha_vdim[112]+alpha_vdim[17]*f[75]+alpha_vdim[6]*f[53]+alpha_vdim[2]*f[29]+alpha_vdim[1]*f[28]+alpha_vdim[0]*f[15]); 
  out[31] += 0.3061862178478971*(f[37]*alpha_cdim[116]+f[16]*alpha_cdim[112])+0.6846531968814573*(alpha_cdim[3]*f[21]+alpha_cdim[0]*f[6]); 
  out[32] += 0.6846531968814573*(f[22]*alpha_cdim[116]+f[6]*alpha_cdim[112])+0.3061862178478971*(alpha_cdim[3]*f[34]+alpha_cdim[0]*f[17]); 
  out[33] += 0.6123724356957944*alpha_cdim[3]*f[35]+0.3061862178478971*(alpha_vdim[2]*f[31]+alpha_vdim[0]*f[16])+0.6846531968814573*alpha_cdim[0]*f[7]+0.273861278752583*alpha_vdim[6]*f[6]+f[1]*(0.6846531968814573*alpha_cdim[3]+0.273861278752583*alpha_vdim[1]); 
  out[34] += 0.6846531968814573*(f[24]*alpha_cdim[116]+f[8]*alpha_cdim[112])+0.3061862178478971*alpha_vdim[1]*f[32]+0.1956151991089878*alpha_vdim[17]*f[17]+0.3061862178478971*(alpha_vdim[0]*f[17]+f[0]*alpha_vdim[17])+0.273861278752583*(alpha_vdim[6]*f[6]+alpha_vdim[2]*f[2]); 
  out[35] += 0.6846531968814573*alpha_vdim[17]*f[57]+0.6123724356957944*(alpha_vdim[6]*f[56]+alpha_vdim[1]*f[33])+0.6846531968814573*alpha_vdim[2]*f[21]+0.3061862178478971*alpha_cdim[0]*f[18]+0.6846531968814573*(alpha_vdim[6]*f[8]+alpha_vdim[0]*f[7])+(0.273861278752583*alpha_cdim[3]+0.6846531968814573*alpha_vdim[1])*f[3]; 
  out[36] += 0.3061862178478971*(f[39]*alpha_cdim[116]+f[18]*alpha_cdim[112])+0.6123724356957944*(alpha_vdim[6]*f[57]+alpha_vdim[2]*f[34])+0.6846531968814573*alpha_vdim[1]*f[21]+0.6123724356957944*f[8]*alpha_vdim[17]+0.6846531968814573*(alpha_vdim[0]*f[8]+alpha_vdim[6]*f[7]+alpha_vdim[2]*f[3]); 
  out[37] += (0.1956151991089878*f[16]+0.3061862178478971*f[0])*alpha_vdim[128]+0.273861278752583*f[6]*alpha_vdim[118]+0.3061862178478971*f[31]*alpha_vdim[114]+0.273861278752583*f[1]*alpha_vdim[113]+0.3061862178478971*f[16]*alpha_vdim[112]+0.6846531968814573*(alpha_cdim[3]*f[23]+alpha_cdim[0]*f[9]); 
  out[38] += 0.273861278752583*f[6]*alpha_vdim[118]+0.6123724356957944*f[41]*alpha_cdim[116]+f[2]*(0.6846531968814573*alpha_cdim[116]+0.273861278752583*alpha_vdim[114])+0.3061862178478971*(f[32]*alpha_vdim[113]+f[17]*alpha_vdim[112])+0.6846531968814573*f[10]*alpha_cdim[112]; 
  out[39] += 0.3061862178478971*(f[58]*alpha_vdim[118]+f[36]*alpha_vdim[114]+f[35]*alpha_vdim[113]+f[18]*alpha_vdim[112])+0.6846531968814573*(alpha_vdim[17]*f[62]+alpha_vdim[6]*f[51]+alpha_vdim[2]*f[24]+alpha_vdim[1]*f[23]+alpha_vdim[0]*f[11]); 
  out[40] += 0.6123724356957944*(f[9]*alpha_vdim[128]+f[59]*alpha_vdim[118])+0.6846531968814573*(f[10]*alpha_vdim[118]+f[22]*alpha_vdim[114])+0.6123724356957944*f[37]*alpha_vdim[113]+0.6846531968814573*(f[4]*alpha_vdim[113]+f[9]*alpha_vdim[112])+0.3061862178478971*(alpha_cdim[3]*f[42]+alpha_cdim[0]*f[19]); 
  out[41] += 0.6846531968814573*f[59]*alpha_vdim[128]+(0.6123724356957944*f[60]+0.6846531968814573*f[9])*alpha_vdim[118]+0.273861278752583*f[4]*alpha_cdim[116]+0.6123724356957944*f[38]*alpha_vdim[114]+0.6846531968814573*(f[4]*alpha_vdim[114]+f[22]*alpha_vdim[113]+f[10]*alpha_vdim[112])+0.3061862178478971*f[19]*alpha_cdim[112]; 
  out[42] += 0.6846531968814573*(f[61]*alpha_vdim[128]+f[51]*alpha_vdim[118]+f[24]*alpha_vdim[114]+f[23]*alpha_vdim[113]+f[11]*alpha_vdim[112])+0.3061862178478971*(alpha_vdim[6]*f[65]+alpha_vdim[2]*f[41]+alpha_vdim[1]*f[40]+alpha_vdim[0]*f[19]); 
  out[43] += 0.6846531968814573*(alpha_cdim[3]*f[26]+alpha_cdim[0]*f[12]); 
  out[44] += 0.6846531968814573*(f[29]*alpha_cdim[116]+f[13]*alpha_cdim[112]); 
  out[45] += 0.6846531968814573*(alpha_vdim[17]*f[71]+alpha_vdim[6]*f[52]+alpha_vdim[2]*f[27]+alpha_vdim[1]*f[26]+alpha_vdim[0]*f[14]); 
  out[46] += 0.6846531968814573*(f[74]*alpha_vdim[128]+f[53]*alpha_vdim[118]+f[29]*alpha_vdim[114]+f[28]*alpha_vdim[113]+f[15]*alpha_vdim[112]); 
  out[47] += 0.3061862178478971*(alpha_cdim[3]*f[49]+alpha_cdim[0]*f[20]); 
  out[48] += 0.3061862178478971*(f[50]*alpha_cdim[116]+f[20]*alpha_cdim[112]); 
  out[49] += 0.3061862178478971*(alpha_vdim[6]*f[80]+alpha_vdim[2]*f[48]+alpha_vdim[1]*f[47]+alpha_vdim[0]*f[20]); 
  out[50] += 0.3061862178478971*(f[80]*alpha_vdim[118]+f[48]*alpha_vdim[114]+f[47]*alpha_vdim[113]+f[20]*alpha_vdim[112]); 
  out[51] += 0.273861278752583*f[21]*alpha_vdim[128]+(0.273861278752583*(f[34]+f[33])+0.3061862178478971*f[3])*alpha_vdim[118]+(0.273861278752583*f[66]+0.3061862178478971*f[7])*alpha_cdim[116]+(0.273861278752583*f[57]+0.3061862178478971*f[7])*alpha_vdim[114]+0.273861278752583*f[56]*alpha_vdim[113]+0.3061862178478971*(f[8]*alpha_vdim[113]+f[21]*alpha_vdim[112]+f[23]*alpha_cdim[112])+0.273861278752583*(alpha_cdim[3]*f[64]+alpha_vdim[2]*f[60]+alpha_vdim[1]*f[59]+alpha_vdim[6]*(f[38]+f[37]))+0.3061862178478971*alpha_cdim[0]*f[24]+0.273861278752583*alpha_vdim[17]*f[22]+0.3061862178478971*(alpha_vdim[0]*f[22]+(alpha_cdim[3]+alpha_vdim[1])*f[10]+alpha_vdim[2]*f[9]+f[4]*alpha_vdim[6]); 
  out[52] += 0.3061862178478971*(f[54]*alpha_cdim[116]+f[26]*alpha_cdim[112])+0.273861278752583*(alpha_cdim[3]*f[73]+alpha_vdim[2]*f[69]+alpha_vdim[1]*f[68]+alpha_vdim[6]*(f[44]+f[43]))+0.3061862178478971*alpha_cdim[0]*f[27]+0.273861278752583*alpha_vdim[17]*f[25]+0.3061862178478971*(alpha_vdim[0]*f[25]+(alpha_cdim[3]+alpha_vdim[1])*f[13]+alpha_vdim[2]*f[12]+f[5]*alpha_vdim[6]); 
  out[53] += 0.273861278752583*f[25]*alpha_vdim[128]+(0.273861278752583*(f[44]+f[43])+0.3061862178478971*f[5])*alpha_vdim[118]+(0.273861278752583*f[77]+0.3061862178478971*f[12])*alpha_cdim[116]+(0.273861278752583*f[69]+0.3061862178478971*f[12])*alpha_vdim[114]+0.273861278752583*f[68]*alpha_vdim[113]+0.3061862178478971*(f[13]*alpha_vdim[113]+f[25]*alpha_vdim[112]+f[28]*alpha_cdim[112]+alpha_cdim[3]*f[55]+alpha_cdim[0]*f[29]); 
  out[54] += 0.273861278752583*(f[26]*alpha_vdim[128]+f[91]*alpha_vdim[118])+0.3061862178478971*(f[27]*alpha_vdim[118]+f[52]*alpha_vdim[114])+0.273861278752583*f[70]*alpha_vdim[113]+0.3061862178478971*(f[14]*alpha_vdim[113]+f[26]*alpha_vdim[112]+alpha_vdim[17]*f[95])+0.273861278752583*(alpha_vdim[6]*f[94]+alpha_cdim[3]*f[76]+alpha_vdim[1]*f[74])+0.3061862178478971*(alpha_vdim[2]*f[53]+alpha_cdim[0]*f[30]+alpha_vdim[6]*f[29]+alpha_vdim[0]*f[28]+(alpha_cdim[3]+alpha_vdim[1])*f[15]); 
  out[55] += 0.3061862178478971*f[91]*alpha_vdim[128]+(0.273861278752583*f[92]+0.3061862178478971*f[26])*alpha_vdim[118]+(0.273861278752583*f[79]+0.3061862178478971*f[14])*alpha_cdim[116]+0.273861278752583*f[71]*alpha_vdim[114]+0.3061862178478971*(f[14]*alpha_vdim[114]+f[52]*alpha_vdim[113]+f[27]*alpha_vdim[112]+f[30]*alpha_cdim[112])+0.273861278752583*(alpha_vdim[6]*f[95]+alpha_vdim[2]*f[75])+0.3061862178478971*alpha_vdim[1]*f[53]+0.273861278752583*alpha_vdim[17]*f[29]+0.3061862178478971*(alpha_vdim[0]*f[29]+alpha_vdim[6]*f[28]+alpha_vdim[2]*f[15]); 
  out[56] += 0.3061862178478971*(f[61]*alpha_cdim[116]+f[33]*alpha_cdim[112])+0.6123724356957944*alpha_cdim[3]*f[58]+0.2449489742783178*alpha_vdim[6]*f[32]+(0.273861278752583*alpha_vdim[17]+0.3061862178478971*alpha_vdim[0])*f[31]+0.6846531968814573*alpha_cdim[0]*f[21]+0.3061862178478971*alpha_vdim[2]*f[16]+0.6846531968814573*alpha_cdim[3]*f[6]+0.273861278752583*(alpha_vdim[1]*f[6]+f[1]*alpha_vdim[6]); 
  out[57] += 0.6846531968814573*(f[51]*alpha_cdim[116]+f[21]*alpha_cdim[112])+0.3061862178478971*alpha_cdim[0]*f[34]+(0.1956151991089878*alpha_vdim[17]+0.3061862178478971*alpha_vdim[0])*f[32]+0.2449489742783178*alpha_vdim[6]*f[31]+0.3061862178478971*((alpha_cdim[3]+alpha_vdim[1])*f[17]+f[1]*alpha_vdim[17])+0.273861278752583*(alpha_vdim[2]*f[6]+f[2]*alpha_vdim[6]); 
  out[58] += 0.3061862178478971*(f[63]*alpha_cdim[116]+f[35]*alpha_cdim[112])+0.6123724356957944*(alpha_vdim[2]*f[57]+alpha_vdim[1]*f[56])+0.3061862178478971*alpha_cdim[0]*f[36]+0.6123724356957944*alpha_vdim[6]*(f[34]+f[33])+(0.6123724356957944*alpha_vdim[17]+0.6846531968814573*alpha_vdim[0])*f[21]+0.273861278752583*alpha_cdim[3]*f[8]+0.6846531968814573*(alpha_vdim[1]*f[8]+alpha_vdim[2]*f[7]+f[3]*alpha_vdim[6]); 
  out[59] += (0.1956151991089878*f[31]+0.3061862178478971*f[2])*alpha_vdim[128]+(0.2449489742783178*f[32]+0.273861278752583*f[1])*alpha_vdim[118]+0.3061862178478971*f[16]*(alpha_cdim[116]+alpha_vdim[114])+0.273861278752583*f[6]*alpha_vdim[113]+0.3061862178478971*(f[31]*alpha_vdim[112]+f[37]*alpha_cdim[112])+0.6846531968814573*(alpha_cdim[3]*f[51]+alpha_cdim[0]*f[22]); 
  out[60] += 0.273861278752583*f[32]*alpha_vdim[128]+(0.2449489742783178*f[31]+0.273861278752583*f[2])*alpha_vdim[118]+0.6123724356957944*f[65]*alpha_cdim[116]+f[6]*(0.6846531968814573*alpha_cdim[116]+0.273861278752583*alpha_vdim[114])+0.3061862178478971*(f[17]*alpha_vdim[113]+f[32]*alpha_vdim[112])+0.6846531968814573*f[22]*alpha_cdim[112]+0.3061862178478971*(alpha_cdim[3]*f[62]+alpha_cdim[0]*f[38]); 
  out[61] += (0.1956151991089878*f[33]+0.3061862178478971*f[3])*alpha_vdim[128]+0.273861278752583*f[21]*alpha_vdim[118]+0.3061862178478971*f[56]*alpha_vdim[114]+0.273861278752583*f[7]*alpha_vdim[113]+0.3061862178478971*f[33]*alpha_vdim[112]+0.6123724356957944*alpha_cdim[3]*f[63]+0.3061862178478971*(alpha_vdim[2]*f[59]+alpha_vdim[0]*f[37])+0.6846531968814573*alpha_cdim[0]*f[23]+0.273861278752583*alpha_vdim[6]*f[22]+(0.6846531968814573*alpha_cdim[3]+0.273861278752583*alpha_vdim[1])*f[9]; 
  out[62] += 0.273861278752583*f[21]*alpha_vdim[118]+0.6123724356957944*f[67]*alpha_cdim[116]+f[8]*(0.6846531968814573*alpha_cdim[116]+0.273861278752583*alpha_vdim[114])+0.3061862178478971*(f[57]*alpha_vdim[113]+f[34]*alpha_vdim[112])+0.6846531968814573*f[24]*alpha_cdim[112]+0.3061862178478971*alpha_vdim[1]*f[60]+(0.1956151991089878*alpha_vdim[17]+0.3061862178478971*alpha_vdim[0])*f[38]+0.273861278752583*alpha_vdim[6]*f[22]+0.3061862178478971*f[4]*alpha_vdim[17]+0.273861278752583*alpha_vdim[2]*f[10]; 
  out[63] += 0.273861278752583*f[35]*alpha_vdim[128]+0.3061862178478971*(f[36]*alpha_vdim[118]+f[58]*alpha_vdim[114]+f[18]*alpha_vdim[113]+f[35]*alpha_vdim[112])+0.6846531968814573*alpha_vdim[17]*f[88]+0.6123724356957944*(alpha_vdim[6]*f[87]+alpha_vdim[1]*f[61])+0.6846531968814573*alpha_vdim[2]*f[51]+0.3061862178478971*alpha_cdim[0]*f[39]+0.6846531968814573*(alpha_vdim[6]*f[24]+alpha_vdim[0]*f[23])+(0.273861278752583*alpha_cdim[3]+0.6846531968814573*alpha_vdim[1])*f[11]; 
  out[64] += 0.3061862178478971*(f[35]*alpha_vdim[118]+f[18]*(alpha_cdim[116]+alpha_vdim[114])+f[58]*alpha_vdim[113]+f[36]*alpha_vdim[112]+f[39]*alpha_cdim[112])+0.6123724356957944*(alpha_vdim[6]*f[88]+alpha_vdim[2]*f[62])+0.6846531968814573*alpha_vdim[1]*f[51]+0.6123724356957944*alpha_vdim[17]*f[24]+0.6846531968814573*(alpha_vdim[0]*f[24]+alpha_vdim[6]*f[23]+alpha_vdim[2]*f[11]); 
  out[65] += 0.6123724356957944*f[22]*alpha_vdim[128]+(0.6123724356957944*(f[38]+f[37])+0.6846531968814573*f[4])*alpha_vdim[118]+0.273861278752583*f[9]*alpha_cdim[116]+(0.6123724356957944*f[60]+0.6846531968814573*f[9])*alpha_vdim[114]+0.6123724356957944*f[59]*alpha_vdim[113]+0.6846531968814573*(f[10]*alpha_vdim[113]+f[22]*alpha_vdim[112])+0.3061862178478971*(f[40]*alpha_cdim[112]+alpha_cdim[3]*f[67]+alpha_cdim[0]*f[41]); 
  out[66] += 0.6123724356957944*(f[23]*alpha_vdim[128]+f[87]*alpha_vdim[118])+0.6846531968814573*(f[24]*alpha_vdim[118]+f[51]*alpha_vdim[114])+0.6123724356957944*f[61]*alpha_vdim[113]+0.6846531968814573*(f[11]*alpha_vdim[113]+f[23]*alpha_vdim[112])+0.3061862178478971*(alpha_vdim[2]*f[65]+alpha_cdim[0]*f[42]+alpha_vdim[6]*f[41]+alpha_vdim[0]*f[40]+(alpha_cdim[3]+alpha_vdim[1])*f[19]); 
  out[67] += 0.6846531968814573*f[87]*alpha_vdim[128]+(0.6123724356957944*f[88]+0.6846531968814573*f[23])*alpha_vdim[118]+0.273861278752583*f[11]*alpha_cdim[116]+0.6123724356957944*f[62]*alpha_vdim[114]+0.6846531968814573*(f[11]*alpha_vdim[114]+f[51]*alpha_vdim[113]+f[24]*alpha_vdim[112])+0.3061862178478971*(f[42]*alpha_cdim[112]+alpha_vdim[1]*f[65])+0.273861278752583*alpha_vdim[17]*f[41]+0.3061862178478971*(alpha_vdim[0]*f[41]+alpha_vdim[6]*f[40]+alpha_vdim[2]*f[19]); 
  out[68] += 0.3061862178478971*(f[74]*alpha_cdim[116]+f[43]*alpha_cdim[112])+0.6846531968814573*(alpha_cdim[3]*f[52]+alpha_cdim[0]*f[25]); 
  out[69] += 0.6846531968814573*(f[53]*alpha_cdim[116]+f[25]*alpha_cdim[112])+0.3061862178478971*(alpha_cdim[3]*f[71]+alpha_cdim[0]*f[44]); 
  out[70] += 0.6123724356957944*alpha_cdim[3]*f[72]+0.3061862178478971*(alpha_vdim[2]*f[68]+alpha_vdim[0]*f[43])+0.6846531968814573*alpha_cdim[0]*f[26]+0.273861278752583*alpha_vdim[6]*f[25]+(0.6846531968814573*alpha_cdim[3]+0.273861278752583*alpha_vdim[1])*f[12]; 
  out[71] += 0.6846531968814573*(f[55]*alpha_cdim[116]+f[27]*alpha_cdim[112])+0.3061862178478971*alpha_vdim[1]*f[69]+(0.1956151991089878*alpha_vdim[17]+0.3061862178478971*alpha_vdim[0])*f[44]+0.273861278752583*alpha_vdim[6]*f[25]+0.3061862178478971*f[5]*alpha_vdim[17]+0.273861278752583*alpha_vdim[2]*f[13]; 
  out[72] += 0.6846531968814573*alpha_vdim[17]*f[92]+0.6123724356957944*(alpha_vdim[6]*f[91]+alpha_vdim[1]*f[70])+0.6846531968814573*alpha_vdim[2]*f[52]+0.3061862178478971*alpha_cdim[0]*f[45]+0.6846531968814573*(alpha_vdim[6]*f[27]+alpha_vdim[0]*f[26])+(0.273861278752583*alpha_cdim[3]+0.6846531968814573*alpha_vdim[1])*f[14]; 
  out[73] += 0.3061862178478971*(f[76]*alpha_cdim[116]+f[45]*alpha_cdim[112])+0.6123724356957944*(alpha_vdim[6]*f[92]+alpha_vdim[2]*f[71])+0.6846531968814573*alpha_vdim[1]*f[52]+0.6123724356957944*alpha_vdim[17]*f[27]+0.6846531968814573*(alpha_vdim[0]*f[27]+alpha_vdim[6]*f[26]+alpha_vdim[2]*f[14]); 
  out[74] += (0.1956151991089878*f[43]+0.3061862178478971*f[5])*alpha_vdim[128]+0.273861278752583*f[25]*alpha_vdim[118]+0.3061862178478971*f[68]*alpha_vdim[114]+0.273861278752583*f[12]*alpha_vdim[113]+0.3061862178478971*f[43]*alpha_vdim[112]+0.6846531968814573*(alpha_cdim[3]*f[54]+alpha_cdim[0]*f[28]); 
  out[75] += 0.273861278752583*f[25]*alpha_vdim[118]+0.6123724356957944*f[78]*alpha_cdim[116]+f[13]*(0.6846531968814573*alpha_cdim[116]+0.273861278752583*alpha_vdim[114])+0.3061862178478971*(f[69]*alpha_vdim[113]+f[44]*alpha_vdim[112])+0.6846531968814573*f[29]*alpha_cdim[112]; 
  out[76] += 0.3061862178478971*(f[93]*alpha_vdim[118]+f[73]*alpha_vdim[114]+f[72]*alpha_vdim[113]+f[45]*alpha_vdim[112])+0.6846531968814573*(alpha_vdim[17]*f[97]+alpha_vdim[6]*f[86]+alpha_vdim[2]*f[55]+alpha_vdim[1]*f[54]+alpha_vdim[0]*f[30]); 
  out[77] += 0.6123724356957944*(f[28]*alpha_vdim[128]+f[94]*alpha_vdim[118])+0.6846531968814573*(f[29]*alpha_vdim[118]+f[53]*alpha_vdim[114])+0.6123724356957944*f[74]*alpha_vdim[113]+0.6846531968814573*(f[15]*alpha_vdim[113]+f[28]*alpha_vdim[112])+0.3061862178478971*(alpha_cdim[3]*f[79]+alpha_cdim[0]*f[46]); 
  out[78] += 0.6846531968814573*f[94]*alpha_vdim[128]+(0.6123724356957944*f[95]+0.6846531968814573*f[28])*alpha_vdim[118]+0.273861278752583*f[15]*alpha_cdim[116]+0.6123724356957944*f[75]*alpha_vdim[114]+0.6846531968814573*(f[15]*alpha_vdim[114]+f[53]*alpha_vdim[113]+f[29]*alpha_vdim[112])+0.3061862178478971*f[46]*alpha_cdim[112]; 
  out[79] += 0.6846531968814573*(f[96]*alpha_vdim[128]+f[86]*alpha_vdim[118]+f[55]*alpha_vdim[114]+f[54]*alpha_vdim[113]+f[30]*alpha_vdim[112])+0.3061862178478971*(alpha_vdim[6]*f[100]+alpha_vdim[2]*f[78]+alpha_vdim[1]*f[77]+alpha_vdim[0]*f[46]); 
  out[80] += 0.3061862178478971*(f[83]*alpha_cdim[116]+f[47]*alpha_cdim[112]+alpha_cdim[3]*f[82]+alpha_cdim[0]*f[48]); 
  out[81] += 0.3061862178478971*(alpha_vdim[2]*f[80]+alpha_cdim[0]*f[49]+alpha_vdim[6]*f[48]+alpha_vdim[0]*f[47]+(alpha_cdim[3]+alpha_vdim[1])*f[20]); 
  out[82] += 0.3061862178478971*(f[85]*alpha_cdim[116]+f[49]*alpha_cdim[112]+alpha_vdim[1]*f[80])+0.273861278752583*alpha_vdim[17]*f[48]+0.3061862178478971*(alpha_vdim[0]*f[48]+alpha_vdim[6]*f[47]+alpha_vdim[2]*f[20]); 
  out[83] += 0.273861278752583*f[47]*alpha_vdim[128]+0.3061862178478971*(f[48]*alpha_vdim[118]+f[80]*alpha_vdim[114]+f[20]*alpha_vdim[113]+f[47]*alpha_vdim[112]+alpha_cdim[3]*f[85]+alpha_cdim[0]*f[50]); 
  out[84] += 0.3061862178478971*(f[47]*alpha_vdim[118]+f[20]*(alpha_cdim[116]+alpha_vdim[114])+f[80]*alpha_vdim[113]+f[48]*alpha_vdim[112]+f[50]*alpha_cdim[112]); 
  out[85] += 0.3061862178478971*(f[103]*alpha_vdim[118]+f[82]*alpha_vdim[114]+f[81]*alpha_vdim[113]+f[49]*alpha_vdim[112]+alpha_vdim[6]*f[104]+alpha_vdim[2]*f[84]+alpha_vdim[1]*f[83]+alpha_vdim[0]*f[50]); 
  out[86] += 0.273861278752583*f[52]*alpha_vdim[128]+(0.273861278752583*(f[71]+f[70])+0.3061862178478971*f[14])*alpha_vdim[118]+(0.273861278752583*f[101]+0.3061862178478971*f[26])*alpha_cdim[116]+(0.273861278752583*f[92]+0.3061862178478971*f[26])*alpha_vdim[114]+0.273861278752583*f[91]*alpha_vdim[113]+0.3061862178478971*(f[27]*alpha_vdim[113]+f[52]*alpha_vdim[112]+f[54]*alpha_cdim[112])+0.273861278752583*(alpha_cdim[3]*f[99]+alpha_vdim[2]*f[95]+alpha_vdim[1]*f[94]+alpha_vdim[6]*(f[75]+f[74]))+0.3061862178478971*alpha_cdim[0]*f[55]+0.273861278752583*alpha_vdim[17]*f[53]+0.3061862178478971*(alpha_vdim[0]*f[53]+(alpha_cdim[3]+alpha_vdim[1])*f[29]+alpha_vdim[2]*f[28]+alpha_vdim[6]*f[15]); 
  out[87] += (0.1956151991089878*f[56]+0.3061862178478971*f[8])*alpha_vdim[128]+(0.2449489742783178*f[57]+0.273861278752583*f[7])*alpha_vdim[118]+0.3061862178478971*f[33]*(alpha_cdim[116]+alpha_vdim[114])+0.273861278752583*f[21]*alpha_vdim[113]+0.3061862178478971*(f[56]*alpha_vdim[112]+f[61]*alpha_cdim[112])+0.6123724356957944*alpha_cdim[3]*f[89]+0.2449489742783178*alpha_vdim[6]*f[60]+(0.273861278752583*alpha_vdim[17]+0.3061862178478971*alpha_vdim[0])*f[59]+0.6846531968814573*alpha_cdim[0]*f[51]+0.3061862178478971*alpha_vdim[2]*f[37]+0.6846531968814573*alpha_cdim[3]*f[22]+0.273861278752583*(alpha_vdim[1]*f[22]+alpha_vdim[6]*f[9]); 
  out[88] += 0.273861278752583*f[57]*alpha_vdim[128]+(0.2449489742783178*f[56]+0.273861278752583*f[8])*alpha_vdim[118]+0.6123724356957944*f[90]*alpha_cdim[116]+f[21]*(0.6846531968814573*alpha_cdim[116]+0.273861278752583*alpha_vdim[114])+0.3061862178478971*(f[34]*alpha_vdim[113]+f[57]*alpha_vdim[112])+0.6846531968814573*f[51]*alpha_cdim[112]+0.3061862178478971*alpha_cdim[0]*f[62]+(0.1956151991089878*alpha_vdim[17]+0.3061862178478971*alpha_vdim[0])*f[60]+0.2449489742783178*alpha_vdim[6]*f[59]+0.3061862178478971*(alpha_cdim[3]+alpha_vdim[1])*f[38]+0.273861278752583*alpha_vdim[2]*f[22]+0.3061862178478971*f[9]*alpha_vdim[17]+0.273861278752583*alpha_vdim[6]*f[10]; 
  out[89] += 0.273861278752583*f[58]*alpha_vdim[128]+0.3061862178478971*(f[18]*alpha_vdim[118]+f[35]*(alpha_cdim[116]+alpha_vdim[114])+f[36]*alpha_vdim[113]+f[58]*alpha_vdim[112]+f[63]*alpha_cdim[112])+0.6123724356957944*(alpha_vdim[2]*f[88]+alpha_vdim[1]*f[87])+0.3061862178478971*alpha_cdim[0]*f[64]+0.6123724356957944*alpha_vdim[6]*(f[62]+f[61])+(0.6123724356957944*alpha_vdim[17]+0.6846531968814573*alpha_vdim[0])*f[51]+0.273861278752583*alpha_cdim[3]*f[24]+0.6846531968814573*(alpha_vdim[1]*f[24]+alpha_vdim[2]*f[23]+alpha_vdim[6]*f[11]); 
  out[90] += 0.6123724356957944*f[51]*alpha_vdim[128]+(0.6123724356957944*(f[62]+f[61])+0.6846531968814573*f[11])*alpha_vdim[118]+0.273861278752583*f[23]*alpha_cdim[116]+(0.6123724356957944*f[88]+0.6846531968814573*f[23])*alpha_vdim[114]+0.6123724356957944*f[87]*alpha_vdim[113]+0.6846531968814573*(f[24]*alpha_vdim[113]+f[51]*alpha_vdim[112])+0.3061862178478971*(f[66]*alpha_cdim[112]+alpha_cdim[0]*f[67])+0.273861278752583*alpha_vdim[17]*f[65]+0.3061862178478971*(alpha_vdim[0]*f[65]+(alpha_cdim[3]+alpha_vdim[1])*f[41]+alpha_vdim[2]*f[40]+alpha_vdim[6]*f[19]); 
  out[91] += 0.3061862178478971*(f[96]*alpha_cdim[116]+f[70]*alpha_cdim[112])+0.6123724356957944*alpha_cdim[3]*f[93]+0.2449489742783178*alpha_vdim[6]*f[69]+(0.273861278752583*alpha_vdim[17]+0.3061862178478971*alpha_vdim[0])*f[68]+0.6846531968814573*alpha_cdim[0]*f[52]+0.3061862178478971*alpha_vdim[2]*f[43]+0.6846531968814573*alpha_cdim[3]*f[25]+0.273861278752583*(alpha_vdim[1]*f[25]+alpha_vdim[6]*f[12]); 
  out[92] += 0.6846531968814573*(f[86]*alpha_cdim[116]+f[52]*alpha_cdim[112])+0.3061862178478971*alpha_cdim[0]*f[71]+(0.1956151991089878*alpha_vdim[17]+0.3061862178478971*alpha_vdim[0])*f[69]+0.2449489742783178*alpha_vdim[6]*f[68]+0.3061862178478971*(alpha_cdim[3]+alpha_vdim[1])*f[44]+0.273861278752583*alpha_vdim[2]*f[25]+0.3061862178478971*f[12]*alpha_vdim[17]+0.273861278752583*alpha_vdim[6]*f[13]; 
  out[93] += 0.3061862178478971*(f[98]*alpha_cdim[116]+f[72]*alpha_cdim[112])+0.6123724356957944*(alpha_vdim[2]*f[92]+alpha_vdim[1]*f[91])+0.3061862178478971*alpha_cdim[0]*f[73]+0.6123724356957944*alpha_vdim[6]*(f[71]+f[70])+(0.6123724356957944*alpha_vdim[17]+0.6846531968814573*alpha_vdim[0])*f[52]+0.273861278752583*alpha_cdim[3]*f[27]+0.6846531968814573*(alpha_vdim[1]*f[27]+alpha_vdim[2]*f[26]+alpha_vdim[6]*f[14]); 
  out[94] += (0.1956151991089878*f[68]+0.3061862178478971*f[13])*alpha_vdim[128]+(0.2449489742783178*f[69]+0.273861278752583*f[12])*alpha_vdim[118]+0.3061862178478971*f[43]*(alpha_cdim[116]+alpha_vdim[114])+0.273861278752583*f[25]*alpha_vdim[113]+0.3061862178478971*(f[68]*alpha_vdim[112]+f[74]*alpha_cdim[112])+0.6846531968814573*(alpha_cdim[3]*f[86]+alpha_cdim[0]*f[53]); 
  out[95] += 0.273861278752583*f[69]*alpha_vdim[128]+(0.2449489742783178*f[68]+0.273861278752583*f[13])*alpha_vdim[118]+0.6123724356957944*f[100]*alpha_cdim[116]+f[25]*(0.6846531968814573*alpha_cdim[116]+0.273861278752583*alpha_vdim[114])+0.3061862178478971*(f[44]*alpha_vdim[113]+f[69]*alpha_vdim[112])+0.6846531968814573*f[53]*alpha_cdim[112]+0.3061862178478971*(alpha_cdim[3]*f[97]+alpha_cdim[0]*f[75]); 
  out[96] += (0.1956151991089878*f[70]+0.3061862178478971*f[14])*alpha_vdim[128]+0.273861278752583*f[52]*alpha_vdim[118]+0.3061862178478971*f[91]*alpha_vdim[114]+0.273861278752583*f[26]*alpha_vdim[113]+0.3061862178478971*f[70]*alpha_vdim[112]+0.6123724356957944*alpha_cdim[3]*f[98]+0.3061862178478971*(alpha_vdim[2]*f[94]+alpha_vdim[0]*f[74])+0.6846531968814573*alpha_cdim[0]*f[54]+0.273861278752583*alpha_vdim[6]*f[53]+(0.6846531968814573*alpha_cdim[3]+0.273861278752583*alpha_vdim[1])*f[28]; 
  out[97] += 0.273861278752583*f[52]*alpha_vdim[118]+0.6123724356957944*f[102]*alpha_cdim[116]+f[27]*(0.6846531968814573*alpha_cdim[116]+0.273861278752583*alpha_vdim[114])+0.3061862178478971*(f[92]*alpha_vdim[113]+f[71]*alpha_vdim[112])+0.6846531968814573*f[55]*alpha_cdim[112]+0.3061862178478971*alpha_vdim[1]*f[95]+(0.1956151991089878*alpha_vdim[17]+0.3061862178478971*alpha_vdim[0])*f[75]+0.273861278752583*(alpha_vdim[6]*f[53]+alpha_vdim[2]*f[29])+0.3061862178478971*f[15]*alpha_vdim[17]; 
  out[98] += 0.273861278752583*f[72]*alpha_vdim[128]+0.3061862178478971*(f[73]*alpha_vdim[118]+f[93]*alpha_vdim[114]+f[45]*alpha_vdim[113]+f[72]*alpha_vdim[112])+0.6846531968814573*alpha_vdim[17]*f[108]+0.6123724356957944*(alpha_vdim[6]*f[107]+alpha_vdim[1]*f[96])+0.6846531968814573*alpha_vdim[2]*f[86]+0.3061862178478971*alpha_cdim[0]*f[76]+0.6846531968814573*(alpha_vdim[6]*f[55]+alpha_vdim[0]*f[54])+(0.273861278752583*alpha_cdim[3]+0.6846531968814573*alpha_vdim[1])*f[30]; 
  out[99] += 0.3061862178478971*(f[72]*alpha_vdim[118]+f[45]*(alpha_cdim[116]+alpha_vdim[114])+f[93]*alpha_vdim[113]+f[73]*alpha_vdim[112]+f[76]*alpha_cdim[112])+0.6123724356957944*(alpha_vdim[6]*f[108]+alpha_vdim[2]*f[97])+0.6846531968814573*alpha_vdim[1]*f[86]+0.6123724356957944*alpha_vdim[17]*f[55]+0.6846531968814573*(alpha_vdim[0]*f[55]+alpha_vdim[6]*f[54]+alpha_vdim[2]*f[30]); 
  out[100] += 0.6123724356957944*f[53]*alpha_vdim[128]+(0.6123724356957944*(f[75]+f[74])+0.6846531968814573*f[15])*alpha_vdim[118]+0.273861278752583*f[28]*alpha_cdim[116]+(0.6123724356957944*f[95]+0.6846531968814573*f[28])*alpha_vdim[114]+0.6123724356957944*f[94]*alpha_vdim[113]+0.6846531968814573*(f[29]*alpha_vdim[113]+f[53]*alpha_vdim[112])+0.3061862178478971*(f[77]*alpha_cdim[112]+alpha_cdim[3]*f[102]+alpha_cdim[0]*f[78]); 
  out[101] += 0.6123724356957944*(f[54]*alpha_vdim[128]+f[107]*alpha_vdim[118])+0.6846531968814573*(f[55]*alpha_vdim[118]+f[86]*alpha_vdim[114])+0.6123724356957944*f[96]*alpha_vdim[113]+0.6846531968814573*(f[30]*alpha_vdim[113]+f[54]*alpha_vdim[112])+0.3061862178478971*(alpha_vdim[2]*f[100]+alpha_cdim[0]*f[79]+alpha_vdim[6]*f[78]+alpha_vdim[0]*f[77]+(alpha_cdim[3]+alpha_vdim[1])*f[46]); 
  out[102] += 0.6846531968814573*f[107]*alpha_vdim[128]+(0.6123724356957944*f[108]+0.6846531968814573*f[54])*alpha_vdim[118]+0.273861278752583*f[30]*alpha_cdim[116]+0.6123724356957944*f[97]*alpha_vdim[114]+0.6846531968814573*(f[30]*alpha_vdim[114]+f[86]*alpha_vdim[113]+f[55]*alpha_vdim[112])+0.3061862178478971*(f[79]*alpha_cdim[112]+alpha_vdim[1]*f[100])+0.273861278752583*alpha_vdim[17]*f[78]+0.3061862178478971*(alpha_vdim[0]*f[78]+alpha_vdim[6]*f[77]+alpha_vdim[2]*f[46]); 
  out[103] += 0.3061862178478971*(f[105]*alpha_cdim[116]+f[81]*alpha_cdim[112]+alpha_cdim[0]*f[82])+0.273861278752583*alpha_vdim[17]*f[80]+0.3061862178478971*(alpha_vdim[0]*f[80]+(alpha_cdim[3]+alpha_vdim[1])*f[48]+alpha_vdim[2]*f[47]+alpha_vdim[6]*f[20]); 
  out[104] += 0.273861278752583*f[80]*alpha_vdim[128]+0.3061862178478971*(f[20]*alpha_vdim[118]+f[47]*(alpha_cdim[116]+alpha_vdim[114])+f[48]*alpha_vdim[113]+f[80]*alpha_vdim[112]+f[83]*alpha_cdim[112]+alpha_cdim[3]*f[106]+alpha_cdim[0]*f[84]); 
  out[105] += 0.273861278752583*f[81]*alpha_vdim[128]+0.3061862178478971*(f[82]*alpha_vdim[118]+f[103]*alpha_vdim[114]+f[49]*alpha_vdim[113]+f[81]*alpha_vdim[112]+alpha_vdim[2]*f[104]+alpha_cdim[0]*f[85]+alpha_vdim[6]*f[84]+alpha_vdim[0]*f[83]+(alpha_cdim[3]+alpha_vdim[1])*f[50]); 
  out[106] += 0.3061862178478971*(f[81]*alpha_vdim[118]+f[49]*(alpha_cdim[116]+alpha_vdim[114])+f[103]*alpha_vdim[113]+f[82]*alpha_vdim[112]+f[85]*alpha_cdim[112]+alpha_vdim[1]*f[104])+0.273861278752583*alpha_vdim[17]*f[84]+0.3061862178478971*(alpha_vdim[0]*f[84]+alpha_vdim[6]*f[83]+alpha_vdim[2]*f[50]); 
  out[107] += (0.1956151991089878*f[91]+0.3061862178478971*f[27])*alpha_vdim[128]+(0.2449489742783178*f[92]+0.273861278752583*f[26])*alpha_vdim[118]+0.3061862178478971*f[70]*(alpha_cdim[116]+alpha_vdim[114])+0.273861278752583*f[52]*alpha_vdim[113]+0.3061862178478971*(f[91]*alpha_vdim[112]+f[96]*alpha_cdim[112])+0.6123724356957944*alpha_cdim[3]*f[109]+0.2449489742783178*alpha_vdim[6]*f[95]+(0.273861278752583*alpha_vdim[17]+0.3061862178478971*alpha_vdim[0])*f[94]+0.6846531968814573*alpha_cdim[0]*f[86]+0.3061862178478971*alpha_vdim[2]*f[74]+0.6846531968814573*alpha_cdim[3]*f[53]+0.273861278752583*(alpha_vdim[1]*f[53]+alpha_vdim[6]*f[28]); 
  out[108] += 0.273861278752583*f[92]*alpha_vdim[128]+(0.2449489742783178*f[91]+0.273861278752583*f[27])*alpha_vdim[118]+0.6123724356957944*f[110]*alpha_cdim[116]+f[52]*(0.6846531968814573*alpha_cdim[116]+0.273861278752583*alpha_vdim[114])+0.3061862178478971*(f[71]*alpha_vdim[113]+f[92]*alpha_vdim[112])+0.6846531968814573*f[86]*alpha_cdim[112]+0.3061862178478971*alpha_cdim[0]*f[97]+(0.1956151991089878*alpha_vdim[17]+0.3061862178478971*alpha_vdim[0])*f[95]+0.2449489742783178*alpha_vdim[6]*f[94]+0.3061862178478971*(alpha_cdim[3]+alpha_vdim[1])*f[75]+0.273861278752583*(alpha_vdim[2]*f[53]+alpha_vdim[6]*f[29])+0.3061862178478971*alpha_vdim[17]*f[28]; 
  out[109] += 0.273861278752583*f[93]*alpha_vdim[128]+0.3061862178478971*(f[45]*alpha_vdim[118]+f[72]*(alpha_cdim[116]+alpha_vdim[114])+f[73]*alpha_vdim[113]+f[93]*alpha_vdim[112]+f[98]*alpha_cdim[112])+0.6123724356957944*(alpha_vdim[2]*f[108]+alpha_vdim[1]*f[107])+0.3061862178478971*alpha_cdim[0]*f[99]+0.6123724356957944*alpha_vdim[6]*(f[97]+f[96])+(0.6123724356957944*alpha_vdim[17]+0.6846531968814573*alpha_vdim[0])*f[86]+0.273861278752583*alpha_cdim[3]*f[55]+0.6846531968814573*(alpha_vdim[1]*f[55]+alpha_vdim[2]*f[54]+alpha_vdim[6]*f[30]); 
  out[110] += 0.6123724356957944*f[86]*alpha_vdim[128]+(0.6123724356957944*(f[97]+f[96])+0.6846531968814573*f[30])*alpha_vdim[118]+0.273861278752583*f[54]*alpha_cdim[116]+(0.6123724356957944*f[108]+0.6846531968814573*f[54])*alpha_vdim[114]+0.6123724356957944*f[107]*alpha_vdim[113]+0.6846531968814573*(f[55]*alpha_vdim[113]+f[86]*alpha_vdim[112])+0.3061862178478971*(f[101]*alpha_cdim[112]+alpha_cdim[0]*f[102])+0.273861278752583*alpha_vdim[17]*f[100]+0.3061862178478971*(alpha_vdim[0]*f[100]+(alpha_cdim[3]+alpha_vdim[1])*f[78]+alpha_vdim[2]*f[77]+alpha_vdim[6]*f[46]); 
  out[111] += 0.273861278752583*f[103]*alpha_vdim[128]+0.3061862178478971*(f[49]*alpha_vdim[118]+f[81]*(alpha_cdim[116]+alpha_vdim[114])+f[82]*alpha_vdim[113]+f[103]*alpha_vdim[112]+f[105]*alpha_cdim[112]+alpha_cdim[0]*f[106])+0.273861278752583*alpha_vdim[17]*f[104]+0.3061862178478971*(alpha_vdim[0]*f[104]+(alpha_cdim[3]+alpha_vdim[1])*f[84]+alpha_vdim[2]*f[83]+alpha_vdim[6]*f[50]); 

  return cflFreq_mid; 
} 

