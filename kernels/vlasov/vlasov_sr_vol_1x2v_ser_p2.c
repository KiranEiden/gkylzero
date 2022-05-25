#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_sr_vol_1x2v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // p_over_gamma: p/gamma (velocity).
  // qmem:      q/m*EM fields.
  // f:         Input distribution function.
  // out:       Incremented output.
  const double dx10 = 2/dxv[0]; 
  const double dv10 = 2/dxv[1]; 
  const double *E0 = &qmem[0]; 
  const double *p0_over_gamma = &p_over_gamma[0]; 
  const double dv11 = 2/dxv[2]; 
  const double *E1 = &qmem[3]; 
  const double *p1_over_gamma = &p_over_gamma[8]; 

  const double *B2 = &qmem[15]; 
  double alpha_mid = 0.0; 
  double alpha_cdim[20] = {0.0}; 
  double alpha_vdim[40] = {0.0}; 

  alpha_cdim[0] = 1.414213562373095*p0_over_gamma[0]*dx10; 
  alpha_cdim[2] = 1.414213562373095*p0_over_gamma[1]*dx10; 
  alpha_cdim[3] = 1.414213562373095*p0_over_gamma[2]*dx10; 
  alpha_cdim[6] = 1.414213562373095*p0_over_gamma[3]*dx10; 
  alpha_cdim[8] = 1.414213562373095*p0_over_gamma[4]*dx10; 
  alpha_cdim[9] = 1.414213562373095*p0_over_gamma[5]*dx10; 
  alpha_cdim[14] = 1.414213562373095*p0_over_gamma[6]*dx10; 
  alpha_cdim[16] = 1.414213562373095*p0_over_gamma[7]*dx10; 
  alpha_mid += fabs(0.1767766952966368*alpha_cdim[0]-0.1976423537605236*(alpha_cdim[9]+alpha_cdim[8])); 

  alpha_vdim[0] = (B2[0]*p1_over_gamma[0]+2.0*E0[0])*dv10; 
  alpha_vdim[1] = (2.0*E0[1]+p1_over_gamma[0]*B2[1])*dv10; 
  alpha_vdim[2] = B2[0]*p1_over_gamma[1]*dv10; 
  alpha_vdim[3] = B2[0]*p1_over_gamma[2]*dv10; 
  alpha_vdim[4] = B2[1]*p1_over_gamma[1]*dv10; 
  alpha_vdim[5] = B2[1]*p1_over_gamma[2]*dv10; 
  alpha_vdim[6] = B2[0]*p1_over_gamma[3]*dv10; 
  alpha_vdim[7] = (2.0*E0[2]+p1_over_gamma[0]*B2[2])*dv10; 
  alpha_vdim[8] = B2[0]*p1_over_gamma[4]*dv10; 
  alpha_vdim[9] = B2[0]*p1_over_gamma[5]*dv10; 
  alpha_vdim[10] = B2[1]*p1_over_gamma[3]*dv10; 
  alpha_vdim[11] = p1_over_gamma[1]*B2[2]*dv10; 
  alpha_vdim[12] = B2[1]*p1_over_gamma[4]*dv10; 
  alpha_vdim[13] = B2[2]*p1_over_gamma[2]*dv10; 
  alpha_vdim[14] = B2[0]*p1_over_gamma[6]*dv10; 
  alpha_vdim[15] = B2[1]*p1_over_gamma[5]*dv10; 
  alpha_vdim[16] = B2[0]*p1_over_gamma[7]*dv10; 
  alpha_vdim[17] = B2[2]*p1_over_gamma[3]*dv10; 
  alpha_vdim[18] = B2[1]*p1_over_gamma[6]*dv10; 
  alpha_vdim[19] = B2[1]*p1_over_gamma[7]*dv10; 
  alpha_mid += fabs(0.1767766952966368*alpha_vdim[0]-0.1976423537605236*(alpha_vdim[9]+alpha_vdim[8]+alpha_vdim[7])); 

  alpha_vdim[20] = (2.0*E1[0]-1.0*B2[0]*p0_over_gamma[0])*dv11; 
  alpha_vdim[21] = (2.0*E1[1]-1.0*p0_over_gamma[0]*B2[1])*dv11; 
  alpha_vdim[22] = -1.0*B2[0]*p0_over_gamma[1]*dv11; 
  alpha_vdim[23] = -1.0*B2[0]*p0_over_gamma[2]*dv11; 
  alpha_vdim[24] = -1.0*B2[1]*p0_over_gamma[1]*dv11; 
  alpha_vdim[25] = -1.0*B2[1]*p0_over_gamma[2]*dv11; 
  alpha_vdim[26] = -1.0*B2[0]*p0_over_gamma[3]*dv11; 
  alpha_vdim[27] = (2.0*E1[2]-1.0*p0_over_gamma[0]*B2[2])*dv11; 
  alpha_vdim[28] = -1.0*B2[0]*p0_over_gamma[4]*dv11; 
  alpha_vdim[29] = -1.0*B2[0]*p0_over_gamma[5]*dv11; 
  alpha_vdim[30] = -1.0*B2[1]*p0_over_gamma[3]*dv11; 
  alpha_vdim[31] = -1.0*p0_over_gamma[1]*B2[2]*dv11; 
  alpha_vdim[32] = -1.0*B2[1]*p0_over_gamma[4]*dv11; 
  alpha_vdim[33] = -1.0*B2[2]*p0_over_gamma[2]*dv11; 
  alpha_vdim[34] = -1.0*B2[0]*p0_over_gamma[6]*dv11; 
  alpha_vdim[35] = -1.0*B2[1]*p0_over_gamma[5]*dv11; 
  alpha_vdim[36] = -1.0*B2[0]*p0_over_gamma[7]*dv11; 
  alpha_vdim[37] = -1.0*B2[2]*p0_over_gamma[3]*dv11; 
  alpha_vdim[38] = -1.0*B2[1]*p0_over_gamma[6]*dv11; 
  alpha_vdim[39] = -1.0*B2[1]*p0_over_gamma[7]*dv11; 
  alpha_mid += fabs(0.1767766952966368*alpha_vdim[20]-0.1976423537605236*(alpha_vdim[29]+alpha_vdim[28]+alpha_vdim[27])); 

  out[1] += 0.6123724356957944*(alpha_cdim[16]*f[16]+alpha_cdim[14]*f[14]+alpha_cdim[9]*f[9]+alpha_cdim[8]*f[8]+alpha_cdim[6]*f[6]+alpha_cdim[3]*f[3]+alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.6123724356957944*(alpha_vdim[19]*f[19]+alpha_vdim[18]*f[18]+alpha_vdim[17]*f[17]+alpha_vdim[16]*f[16]+alpha_vdim[15]*f[15]+alpha_vdim[14]*f[14]+alpha_vdim[13]*f[13]+alpha_vdim[12]*f[12]+alpha_vdim[11]*f[11]+alpha_vdim[10]*f[10]+alpha_vdim[9]*f[9]+alpha_vdim[8]*f[8]+alpha_vdim[7]*f[7]+alpha_vdim[6]*f[6]+alpha_vdim[5]*f[5]+alpha_vdim[4]*f[4]+alpha_vdim[3]*f[3]+alpha_vdim[2]*f[2]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[19]*alpha_vdim[39]+f[18]*alpha_vdim[38]+f[17]*alpha_vdim[37]+f[16]*alpha_vdim[36]+f[15]*alpha_vdim[35]+f[14]*alpha_vdim[34]+f[13]*alpha_vdim[33]+f[12]*alpha_vdim[32]+f[11]*alpha_vdim[31]+f[10]*alpha_vdim[30]+f[9]*alpha_vdim[29]+f[8]*alpha_vdim[28]+f[7]*alpha_vdim[27]+f[6]*alpha_vdim[26]+f[5]*alpha_vdim[25]+f[4]*alpha_vdim[24]+f[3]*alpha_vdim[23]+f[2]*alpha_vdim[22]+f[1]*alpha_vdim[21]+f[0]*alpha_vdim[20]); 
  out[4] += 0.6123724356957944*(alpha_vdim[16]*f[19]+f[16]*alpha_vdim[19]+alpha_vdim[14]*f[18]+f[14]*alpha_vdim[18])+0.5477225575051661*(alpha_vdim[10]*f[17]+f[10]*alpha_vdim[17])+0.6123724356957944*(alpha_cdim[9]*f[16]+f[9]*alpha_cdim[16]+alpha_vdim[9]*f[15]+f[9]*alpha_vdim[15])+0.5477225575051661*(alpha_cdim[6]*f[14]+f[6]*alpha_cdim[14]+alpha_vdim[5]*f[13]+f[5]*alpha_vdim[13])+0.6123724356957944*(alpha_vdim[8]*f[12]+f[8]*alpha_vdim[12])+0.5477225575051661*(alpha_vdim[4]*f[11]+f[4]*alpha_vdim[11])+0.6123724356957944*(alpha_vdim[6]*f[10]+f[6]*alpha_vdim[10])+0.5477225575051661*(alpha_cdim[2]*f[8]+f[2]*alpha_cdim[8]+alpha_vdim[1]*f[7]+f[1]*alpha_vdim[7])+0.6123724356957944*(alpha_cdim[3]*f[6]+f[3]*alpha_cdim[6]+alpha_vdim[3]*f[5]+f[3]*alpha_vdim[5]+alpha_vdim[2]*f[4]+f[2]*(alpha_vdim[4]+alpha_cdim[0])+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[5] += 0.6123724356957944*(f[16]*alpha_vdim[39]+f[14]*alpha_vdim[38])+0.5477225575051661*f[10]*alpha_vdim[37]+0.6123724356957944*(f[19]*alpha_vdim[36]+f[9]*alpha_vdim[35]+f[18]*alpha_vdim[34])+0.5477225575051661*f[5]*alpha_vdim[33]+0.6123724356957944*f[8]*alpha_vdim[32]+0.5477225575051661*(f[4]*alpha_vdim[31]+f[17]*alpha_vdim[30])+0.6123724356957944*(f[6]*alpha_vdim[30]+f[15]*alpha_vdim[29]+f[12]*alpha_vdim[28])+0.5477225575051661*f[1]*alpha_vdim[27]+0.6123724356957944*f[10]*alpha_vdim[26]+(0.5477225575051661*f[13]+0.6123724356957944*f[3])*alpha_vdim[25]+0.5477225575051661*f[11]*alpha_vdim[24]+0.6123724356957944*(f[2]*alpha_vdim[24]+f[5]*alpha_vdim[23]+f[4]*alpha_vdim[22])+0.5477225575051661*f[7]*alpha_vdim[21]+0.6123724356957944*(f[0]*alpha_vdim[21]+f[1]*alpha_vdim[20])+0.5477225575051661*(alpha_cdim[6]*f[16]+f[6]*alpha_cdim[16])+0.6123724356957944*(alpha_cdim[8]*f[14]+f[8]*alpha_cdim[14])+0.5477225575051661*(alpha_cdim[3]*f[9]+f[3]*alpha_cdim[9])+0.6123724356957944*(alpha_cdim[2]*f[6]+f[2]*alpha_cdim[6]+alpha_cdim[0]*f[3]+f[0]*alpha_cdim[3]); 
  out[6] += 0.6123724356957944*f[15]*alpha_vdim[39]+0.5477225575051661*f[10]*alpha_vdim[38]+0.6123724356957944*(f[13]*alpha_vdim[37]+f[9]*alpha_vdim[36]+f[19]*alpha_vdim[35])+0.5477225575051661*f[6]*alpha_vdim[34]+0.6123724356957944*f[17]*alpha_vdim[33]+0.5477225575051661*f[4]*alpha_vdim[32]+0.6123724356957944*f[7]*alpha_vdim[31]+0.5477225575051661*f[18]*alpha_vdim[30]+0.6123724356957944*(f[5]*alpha_vdim[30]+f[16]*alpha_vdim[29])+0.5477225575051661*f[2]*alpha_vdim[28]+0.6123724356957944*f[11]*alpha_vdim[27]+0.5477225575051661*f[14]*alpha_vdim[26]+0.6123724356957944*(f[3]*alpha_vdim[26]+f[10]*alpha_vdim[25])+0.5477225575051661*f[12]*alpha_vdim[24]+0.6123724356957944*(f[1]*alpha_vdim[24]+f[6]*alpha_vdim[23])+0.5477225575051661*f[8]*alpha_vdim[22]+0.6123724356957944*(f[0]*alpha_vdim[22]+f[4]*alpha_vdim[21]+f[2]*alpha_vdim[20])+0.5477225575051661*(alpha_vdim[10]*f[19]+f[10]*alpha_vdim[19])+0.6123724356957944*(alpha_vdim[12]*f[18]+f[12]*alpha_vdim[18]+alpha_vdim[11]*f[17]+f[11]*alpha_vdim[17])+0.5477225575051661*(alpha_vdim[6]*f[16]+f[6]*alpha_vdim[16]+alpha_vdim[5]*f[15]+f[5]*alpha_vdim[15])+0.6123724356957944*(alpha_vdim[8]*f[14]+f[8]*alpha_vdim[14]+alpha_vdim[7]*f[13]+f[7]*alpha_vdim[13]+alpha_vdim[4]*f[10]+f[4]*alpha_vdim[10])+0.5477225575051661*(alpha_vdim[3]*f[9]+f[3]*alpha_vdim[9])+0.6123724356957944*(alpha_vdim[2]*f[6]+f[2]*alpha_vdim[6]+alpha_vdim[1]*f[5]+f[1]*alpha_vdim[5]+alpha_vdim[0]*f[3]+f[0]*alpha_vdim[3]); 
  out[7] += 1.369306393762915*(alpha_cdim[16]*f[19]+alpha_cdim[14]*f[18]+alpha_cdim[9]*f[15]+alpha_cdim[8]*f[12]+alpha_cdim[6]*f[10]+alpha_cdim[3]*f[5]+alpha_cdim[2]*f[4]+alpha_cdim[0]*f[1]); 
  out[8] += 1.369306393762915*(alpha_vdim[15]*f[19]+f[15]*alpha_vdim[19])+1.224744871391589*(alpha_vdim[10]*f[18]+f[10]*alpha_vdim[18])+1.369306393762915*(alpha_vdim[13]*f[17]+f[13]*alpha_vdim[17]+alpha_vdim[9]*f[16]+f[9]*alpha_vdim[16])+1.224744871391589*(alpha_vdim[6]*f[14]+f[6]*alpha_vdim[14]+alpha_vdim[4]*f[12]+f[4]*alpha_vdim[12])+1.369306393762915*(alpha_vdim[7]*f[11]+f[7]*alpha_vdim[11]+alpha_vdim[5]*f[10]+f[5]*alpha_vdim[10])+1.224744871391589*(alpha_vdim[2]*f[8]+f[2]*alpha_vdim[8])+1.369306393762915*(alpha_vdim[3]*f[6]+f[3]*alpha_vdim[6]+alpha_vdim[1]*f[4]+f[1]*alpha_vdim[4]+alpha_vdim[0]*f[2]+f[0]*alpha_vdim[2]); 
  out[9] += 1.224744871391589*f[10]*alpha_vdim[39]+1.369306393762915*(f[12]*alpha_vdim[38]+f[11]*alpha_vdim[37])+1.224744871391589*(f[6]*alpha_vdim[36]+f[5]*alpha_vdim[35])+1.369306393762915*(f[8]*alpha_vdim[34]+f[7]*alpha_vdim[33]+f[18]*alpha_vdim[32]+f[17]*alpha_vdim[31])+(1.224744871391589*f[19]+1.369306393762915*f[4])*alpha_vdim[30]+1.224744871391589*f[3]*alpha_vdim[29]+1.369306393762915*(f[14]*alpha_vdim[28]+f[13]*alpha_vdim[27])+(1.224744871391589*f[16]+1.369306393762915*f[2])*alpha_vdim[26]+1.224744871391589*f[15]*alpha_vdim[25]+1.369306393762915*(f[1]*alpha_vdim[25]+f[10]*alpha_vdim[24])+1.224744871391589*f[9]*alpha_vdim[23]+1.369306393762915*(f[0]*alpha_vdim[23]+f[6]*alpha_vdim[22]+f[5]*alpha_vdim[21]+f[3]*alpha_vdim[20]); 
  out[10] += 0.6123724356957944*f[9]*alpha_vdim[39]+(0.4898979485566357*f[17]+0.5477225575051661*f[6])*alpha_vdim[38]+(0.4898979485566357*f[18]+0.5477225575051661*f[5])*alpha_vdim[37]+0.6123724356957944*(f[15]*alpha_vdim[36]+f[16]*alpha_vdim[35])+0.5477225575051661*f[10]*(alpha_vdim[34]+alpha_vdim[33])+(0.4898979485566357*f[11]+0.5477225575051661*f[2])*alpha_vdim[32]+0.4898979485566357*f[12]*alpha_vdim[31]+0.5477225575051661*(f[1]*alpha_vdim[31]+(f[14]+f[13])*alpha_vdim[30])+0.6123724356957944*(f[3]*alpha_vdim[30]+f[19]*alpha_vdim[29])+0.5477225575051661*f[4]*(alpha_vdim[28]+alpha_vdim[27])+(0.5477225575051661*f[18]+0.6123724356957944*f[5])*alpha_vdim[26]+(0.5477225575051661*f[17]+0.6123724356957944*f[6])*alpha_vdim[25]+0.5477225575051661*(f[8]+f[7])*alpha_vdim[24]+0.6123724356957944*(f[0]*alpha_vdim[24]+f[10]*alpha_vdim[23])+(0.5477225575051661*f[12]+0.6123724356957944*f[1])*alpha_vdim[22]+0.5477225575051661*f[11]*alpha_vdim[21]+0.6123724356957944*(f[2]*alpha_vdim[21]+f[4]*alpha_vdim[20])+(0.4898979485566357*alpha_vdim[17]+0.5477225575051661*alpha_vdim[6])*f[19]+(0.4898979485566357*f[17]+0.5477225575051661*f[6])*alpha_vdim[19]+0.6123724356957944*(alpha_vdim[8]*f[18]+f[8]*alpha_vdim[18])+0.5477225575051661*(alpha_vdim[4]*f[17]+f[4]*alpha_vdim[17])+0.4898979485566357*alpha_cdim[14]*f[16]+0.5477225575051661*((alpha_vdim[10]+alpha_cdim[3])*f[16]+f[10]*alpha_vdim[16])+(0.4898979485566357*f[14]+0.5477225575051661*f[3])*alpha_cdim[16]+(0.4898979485566357*alpha_vdim[13]+0.5477225575051661*alpha_vdim[3])*f[15]+(0.4898979485566357*f[13]+0.5477225575051661*f[3])*alpha_vdim[15]+(0.6123724356957944*alpha_vdim[12]+0.5477225575051661*alpha_cdim[2])*f[14]+0.6123724356957944*f[12]*alpha_vdim[14]+0.5477225575051661*(f[2]*alpha_cdim[14]+alpha_vdim[1]*f[13]+f[1]*alpha_vdim[13]+alpha_vdim[10]*f[11]+f[10]*alpha_vdim[11])+0.6123724356957944*(alpha_vdim[2]*f[10]+f[2]*alpha_vdim[10])+0.5477225575051661*((alpha_cdim[6]+alpha_vdim[5])*f[9]+f[5]*alpha_vdim[9]+f[6]*alpha_cdim[9]+alpha_cdim[6]*f[8]+f[6]*alpha_cdim[8]+alpha_vdim[5]*f[7]+f[5]*alpha_vdim[7])+0.6123724356957944*((alpha_vdim[4]+alpha_cdim[0])*f[6]+f[4]*alpha_vdim[6]+f[0]*alpha_cdim[6]+alpha_vdim[0]*f[5]+f[0]*alpha_vdim[5]+(alpha_cdim[2]+alpha_vdim[1])*f[3]+f[1]*alpha_vdim[3]+f[2]*alpha_cdim[3]); 
  out[11] += (0.5477225575051661*alpha_vdim[19]+1.369306393762915*alpha_cdim[9])*f[19]+(0.5477225575051661*alpha_vdim[18]+1.224744871391589*alpha_cdim[6])*f[18]+0.3912303982179757*alpha_vdim[17]*f[17]+0.6123724356957944*(alpha_vdim[6]*f[17]+f[6]*alpha_vdim[17])+f[15]*(1.369306393762915*alpha_cdim[16]+0.5477225575051661*alpha_vdim[15])+1.224744871391589*f[10]*alpha_cdim[14]+0.3912303982179757*alpha_vdim[13]*f[13]+0.6123724356957944*(alpha_vdim[3]*f[13]+f[3]*alpha_vdim[13])+(0.5477225575051661*alpha_vdim[12]+1.224744871391589*alpha_cdim[2])*f[12]+0.3912303982179757*alpha_vdim[11]*f[11]+0.6123724356957944*(alpha_vdim[2]*f[11]+f[2]*alpha_vdim[11])+(0.5477225575051661*alpha_vdim[10]+1.369306393762915*alpha_cdim[3])*f[10]+1.224744871391589*f[4]*alpha_cdim[8]+0.3912303982179757*alpha_vdim[7]*f[7]+0.6123724356957944*(alpha_vdim[0]*f[7]+f[0]*alpha_vdim[7])+f[5]*(1.369306393762915*alpha_cdim[6]+0.5477225575051661*alpha_vdim[5])+(0.5477225575051661*alpha_vdim[4]+1.369306393762915*alpha_cdim[0])*f[4]+f[1]*(1.369306393762915*alpha_cdim[2]+0.5477225575051661*alpha_vdim[1]); 
  out[12] += 1.369306393762915*(alpha_vdim[9]*f[19]+f[9]*alpha_vdim[19])+(1.095445115010332*alpha_vdim[17]+1.224744871391589*alpha_vdim[6])*f[18]+1.095445115010332*f[17]*alpha_vdim[18]+1.224744871391589*(f[6]*alpha_vdim[18]+alpha_vdim[5]*f[17]+f[5]*alpha_vdim[17])+0.5477225575051661*alpha_cdim[16]*f[16]+1.369306393762915*(alpha_vdim[15]*f[16]+f[15]*alpha_vdim[16])+(0.3912303982179757*alpha_cdim[14]+1.224744871391589*alpha_vdim[10]+0.6123724356957944*alpha_cdim[3])*f[14]+1.224744871391589*f[10]*alpha_vdim[14]+0.6123724356957944*f[3]*alpha_cdim[14]+1.224744871391589*(alpha_vdim[10]*f[13]+f[10]*alpha_vdim[13])+(1.095445115010332*alpha_vdim[11]+1.224744871391589*alpha_vdim[2])*f[12]+1.095445115010332*f[11]*alpha_vdim[12]+1.224744871391589*(f[2]*alpha_vdim[12]+alpha_vdim[1]*f[11]+f[1]*alpha_vdim[11])+1.369306393762915*(alpha_vdim[3]*f[10]+f[3]*alpha_vdim[10])+(0.3912303982179757*alpha_cdim[8]+1.224744871391589*alpha_vdim[4]+0.6123724356957944*alpha_cdim[0])*f[8]+1.224744871391589*f[4]*alpha_vdim[8]+0.6123724356957944*f[0]*alpha_cdim[8]+1.224744871391589*(alpha_vdim[4]*f[7]+f[4]*alpha_vdim[7])+0.5477225575051661*alpha_cdim[6]*f[6]+1.369306393762915*(alpha_vdim[5]*f[6]+f[5]*alpha_vdim[6]+alpha_vdim[0]*f[4]+f[0]*alpha_vdim[4])+0.5477225575051661*alpha_cdim[2]*f[2]+1.369306393762915*(alpha_vdim[1]*f[2]+f[1]*alpha_vdim[2]); 
  out[13] += 0.5477225575051661*(f[19]*alpha_vdim[39]+f[18]*alpha_vdim[38])+(0.3912303982179757*f[17]+0.6123724356957944*f[6])*alpha_vdim[37]+0.5477225575051661*f[15]*alpha_vdim[35]+(0.3912303982179757*f[13]+0.6123724356957944*f[3])*alpha_vdim[33]+0.5477225575051661*f[12]*alpha_vdim[32]+(0.3912303982179757*f[11]+0.6123724356957944*f[2])*alpha_vdim[31]+0.5477225575051661*f[10]*alpha_vdim[30]+0.3912303982179757*f[7]*alpha_vdim[27]+0.6123724356957944*(f[0]*alpha_vdim[27]+f[17]*alpha_vdim[26])+0.5477225575051661*(f[5]*alpha_vdim[25]+f[4]*alpha_vdim[24])+0.6123724356957944*(f[13]*alpha_vdim[23]+f[11]*alpha_vdim[22])+0.5477225575051661*f[1]*alpha_vdim[21]+0.6123724356957944*f[7]*alpha_vdim[20]+1.224744871391589*alpha_cdim[6]*f[19]+1.369306393762915*alpha_cdim[8]*f[18]+1.224744871391589*(f[10]*alpha_cdim[16]+alpha_cdim[3]*f[15])+1.369306393762915*(f[12]*alpha_cdim[14]+alpha_cdim[2]*f[10])+1.224744871391589*f[5]*alpha_cdim[9]+1.369306393762915*(f[4]*alpha_cdim[6]+alpha_cdim[0]*f[5]+f[1]*alpha_cdim[3]); 
  out[14] += 0.5477225575051661*f[19]*alpha_vdim[39]+(0.3912303982179757*f[18]+0.6123724356957944*f[5])*alpha_vdim[38]+0.5477225575051661*(f[17]*alpha_vdim[37]+f[16]*alpha_vdim[36])+(0.3912303982179757*f[14]+0.6123724356957944*f[3])*alpha_vdim[34]+(0.3912303982179757*f[12]+0.6123724356957944*f[1])*alpha_vdim[32]+0.5477225575051661*(f[11]*alpha_vdim[31]+f[10]*alpha_vdim[30])+(0.3912303982179757*f[8]+0.6123724356957944*f[0])*alpha_vdim[28]+0.5477225575051661*f[6]*alpha_vdim[26]+0.6123724356957944*f[18]*alpha_vdim[25]+0.5477225575051661*f[4]*alpha_vdim[24]+0.6123724356957944*f[14]*alpha_vdim[23]+0.5477225575051661*f[2]*alpha_vdim[22]+0.6123724356957944*(f[12]*alpha_vdim[21]+f[8]*alpha_vdim[20])+(1.095445115010332*alpha_vdim[18]+1.224744871391589*alpha_vdim[5])*f[19]+1.095445115010332*f[18]*alpha_vdim[19]+1.224744871391589*(f[5]*alpha_vdim[19]+alpha_vdim[4]*f[18]+f[4]*alpha_vdim[18])+1.369306393762915*(alpha_vdim[7]*f[17]+f[7]*alpha_vdim[17])+(1.095445115010332*alpha_vdim[14]+1.224744871391589*alpha_vdim[3])*f[16]+1.095445115010332*f[14]*alpha_vdim[16]+1.224744871391589*(f[3]*alpha_vdim[16]+alpha_vdim[10]*f[15]+f[10]*alpha_vdim[15]+alpha_vdim[2]*f[14]+f[2]*alpha_vdim[14])+1.369306393762915*(alpha_vdim[11]*f[13]+f[11]*alpha_vdim[13])+1.224744871391589*(alpha_vdim[10]*f[12]+f[10]*alpha_vdim[12])+1.369306393762915*(alpha_vdim[1]*f[10]+f[1]*alpha_vdim[10])+1.224744871391589*(alpha_vdim[6]*f[9]+f[6]*alpha_vdim[9]+alpha_vdim[6]*f[8]+f[6]*alpha_vdim[8])+1.369306393762915*(alpha_vdim[0]*f[6]+f[0]*alpha_vdim[6]+alpha_vdim[4]*f[5]+f[4]*alpha_vdim[5]+alpha_vdim[2]*f[3]+f[2]*alpha_vdim[3]); 
  out[15] += (1.095445115010332*f[17]+1.224744871391589*f[6])*alpha_vdim[39]+1.369306393762915*f[8]*alpha_vdim[38]+1.095445115010332*f[19]*alpha_vdim[37]+1.224744871391589*(f[4]*alpha_vdim[37]+f[10]*alpha_vdim[36])+(1.095445115010332*f[13]+1.224744871391589*f[3])*alpha_vdim[35]+1.369306393762915*f[12]*alpha_vdim[34]+(1.095445115010332*f[15]+1.224744871391589*f[1])*alpha_vdim[33]+1.369306393762915*f[14]*alpha_vdim[32]+1.224744871391589*f[10]*alpha_vdim[31]+(1.224744871391589*(f[16]+f[11])+1.369306393762915*f[2])*alpha_vdim[30]+1.224744871391589*f[5]*alpha_vdim[29]+1.369306393762915*f[18]*alpha_vdim[28]+1.224744871391589*f[5]*alpha_vdim[27]+(1.224744871391589*f[19]+1.369306393762915*f[4])*alpha_vdim[26]+(1.224744871391589*(f[9]+f[7])+1.369306393762915*f[0])*alpha_vdim[25]+(1.224744871391589*f[17]+1.369306393762915*f[6])*alpha_vdim[24]+1.224744871391589*f[15]*alpha_vdim[23]+1.369306393762915*(f[1]*alpha_vdim[23]+f[10]*alpha_vdim[22])+1.224744871391589*f[13]*alpha_vdim[21]+1.369306393762915*(f[3]*alpha_vdim[21]+f[5]*alpha_vdim[20])+0.3912303982179757*alpha_cdim[16]*f[16]+0.6123724356957944*(alpha_cdim[2]*f[16]+f[2]*alpha_cdim[16])+0.5477225575051661*alpha_cdim[14]*f[14]+0.3912303982179757*alpha_cdim[9]*f[9]+0.6123724356957944*(alpha_cdim[0]*f[9]+f[0]*alpha_cdim[9])+0.5477225575051661*(alpha_cdim[6]*f[6]+alpha_cdim[3]*f[3]); 
  out[16] += (1.095445115010332*f[18]+1.224744871391589*f[5])*alpha_vdim[39]+(1.095445115010332*f[19]+1.224744871391589*f[4])*alpha_vdim[38]+1.369306393762915*f[7]*alpha_vdim[37]+1.095445115010332*f[14]*alpha_vdim[36]+1.224744871391589*(f[3]*alpha_vdim[36]+f[10]*alpha_vdim[35])+(1.095445115010332*f[16]+1.224744871391589*f[2])*alpha_vdim[34]+1.369306393762915*f[11]*alpha_vdim[33]+1.224744871391589*f[10]*alpha_vdim[32]+1.369306393762915*f[13]*alpha_vdim[31]+(1.224744871391589*(f[15]+f[12])+1.369306393762915*f[1])*alpha_vdim[30]+1.224744871391589*f[6]*(alpha_vdim[29]+alpha_vdim[28])+1.369306393762915*f[17]*alpha_vdim[27]+(1.224744871391589*(f[9]+f[8])+1.369306393762915*f[0])*alpha_vdim[26]+(1.224744871391589*f[19]+1.369306393762915*f[4])*alpha_vdim[25]+(1.224744871391589*f[18]+1.369306393762915*f[5])*alpha_vdim[24]+(1.224744871391589*f[16]+1.369306393762915*f[2])*alpha_vdim[23]+1.224744871391589*f[14]*alpha_vdim[22]+1.369306393762915*(f[3]*alpha_vdim[22]+f[10]*alpha_vdim[21]+f[6]*alpha_vdim[20])+0.3912303982179757*alpha_vdim[19]*f[19]+0.6123724356957944*(alpha_vdim[4]*f[19]+f[4]*alpha_vdim[19])+0.5477225575051661*(alpha_vdim[18]*f[18]+alpha_vdim[17]*f[17])+0.3912303982179757*alpha_vdim[16]*f[16]+0.6123724356957944*(alpha_vdim[2]*f[16]+f[2]*alpha_vdim[16])+0.3912303982179757*alpha_vdim[15]*f[15]+0.6123724356957944*(alpha_vdim[1]*f[15]+f[1]*alpha_vdim[15])+0.5477225575051661*(alpha_vdim[14]*f[14]+alpha_vdim[13]*f[13]+alpha_vdim[10]*f[10])+0.3912303982179757*alpha_vdim[9]*f[9]+0.6123724356957944*(alpha_vdim[0]*f[9]+f[0]*alpha_vdim[9])+0.5477225575051661*(alpha_vdim[6]*f[6]+alpha_vdim[5]*f[5]+alpha_vdim[3]*f[3]); 
  out[17] += 0.5477225575051661*f[15]*alpha_vdim[39]+0.4898979485566357*f[10]*alpha_vdim[38]+(0.5477225575051661*f[14]+0.3912303982179757*f[13]+0.6123724356957944*f[3])*alpha_vdim[37]+0.5477225575051661*(f[19]*alpha_vdim[35]+f[17]*alpha_vdim[34])+(0.3912303982179757*f[17]+0.6123724356957944*f[6])*alpha_vdim[33]+0.4898979485566357*f[4]*alpha_vdim[32]+(0.5477225575051661*f[8]+0.3912303982179757*f[7]+0.6123724356957944*f[0])*alpha_vdim[31]+(0.4898979485566357*f[18]+0.5477225575051661*f[5])*alpha_vdim[30]+f[11]*(0.5477225575051661*alpha_vdim[28]+0.3912303982179757*alpha_vdim[27])+0.6123724356957944*(f[2]*alpha_vdim[27]+f[13]*alpha_vdim[26])+0.5477225575051661*f[10]*alpha_vdim[25]+(0.4898979485566357*f[12]+0.5477225575051661*f[1])*alpha_vdim[24]+0.6123724356957944*(f[17]*alpha_vdim[23]+f[7]*alpha_vdim[22])+0.5477225575051661*f[4]*alpha_vdim[21]+0.6123724356957944*f[11]*alpha_vdim[20]+(1.095445115010332*alpha_cdim[14]+0.4898979485566357*alpha_vdim[10]+1.224744871391589*alpha_cdim[3])*f[19]+0.4898979485566357*f[10]*alpha_vdim[19]+(1.095445115010332*alpha_cdim[16]+0.5477225575051661*alpha_vdim[12]+1.224744871391589*alpha_cdim[2])*f[18]+0.5477225575051661*f[12]*alpha_vdim[18]+(0.5477225575051661*alpha_vdim[16]+0.3912303982179757*alpha_vdim[11]+0.6123724356957944*alpha_vdim[2])*f[17]+(0.5477225575051661*f[16]+0.3912303982179757*f[11]+0.6123724356957944*f[2])*alpha_vdim[17]+1.224744871391589*(f[5]*alpha_cdim[16]+alpha_cdim[6]*f[15])+0.4898979485566357*(alpha_vdim[5]*f[15]+f[5]*alpha_vdim[15])+1.224744871391589*f[4]*alpha_cdim[14]+(0.5477225575051661*alpha_vdim[9]+0.3912303982179757*alpha_vdim[7]+0.6123724356957944*alpha_vdim[0])*f[13]+(0.5477225575051661*f[9]+0.3912303982179757*f[7]+0.6123724356957944*f[0])*alpha_vdim[13]+1.224744871391589*alpha_cdim[6]*f[12]+0.6123724356957944*(alpha_vdim[6]*f[11]+f[6]*alpha_vdim[11])+(1.224744871391589*(alpha_cdim[9]+alpha_cdim[8])+0.5477225575051661*alpha_vdim[4]+1.369306393762915*alpha_cdim[0])*f[10]+0.5477225575051661*f[4]*alpha_vdim[10]+0.6123724356957944*(alpha_vdim[3]*f[7]+f[3]*alpha_vdim[7])+1.369306393762915*(f[1]*alpha_cdim[6]+alpha_cdim[2]*f[5])+0.5477225575051661*(alpha_vdim[1]*f[5]+f[1]*alpha_vdim[5])+1.369306393762915*alpha_cdim[3]*f[4]; 
  out[18] += 0.5477225575051661*f[16]*alpha_vdim[39]+(0.3912303982179757*f[14]+0.5477225575051661*f[13]+0.6123724356957944*f[3])*alpha_vdim[38]+0.4898979485566357*f[10]*alpha_vdim[37]+0.5477225575051661*f[19]*alpha_vdim[36]+(0.3912303982179757*f[18]+0.6123724356957944*f[5])*alpha_vdim[34]+0.5477225575051661*f[18]*alpha_vdim[33]+(0.3912303982179757*f[8]+0.5477225575051661*f[7]+0.6123724356957944*f[0])*alpha_vdim[32]+0.4898979485566357*f[4]*alpha_vdim[31]+(0.4898979485566357*f[17]+0.5477225575051661*f[6])*alpha_vdim[30]+(0.3912303982179757*f[12]+0.6123724356957944*f[1])*alpha_vdim[28]+0.5477225575051661*(f[12]*alpha_vdim[27]+f[10]*alpha_vdim[26])+0.6123724356957944*f[14]*alpha_vdim[25]+(0.4898979485566357*f[11]+0.5477225575051661*f[2])*alpha_vdim[24]+0.6123724356957944*f[18]*alpha_vdim[23]+0.5477225575051661*f[4]*alpha_vdim[22]+0.6123724356957944*(f[8]*alpha_vdim[21]+f[12]*alpha_vdim[20])+(1.095445115010332*(alpha_vdim[14]+alpha_vdim[13])+1.224744871391589*alpha_vdim[3])*f[19]+(1.095445115010332*(f[14]+f[13])+1.224744871391589*f[3])*alpha_vdim[19]+(1.095445115010332*(alpha_vdim[16]+alpha_vdim[11])+1.224744871391589*alpha_vdim[2])*f[18]+(1.095445115010332*(f[16]+f[11])+1.224744871391589*f[2])*alpha_vdim[18]+(1.095445115010332*(alpha_vdim[15]+alpha_vdim[12])+1.224744871391589*alpha_vdim[1])*f[17]+(1.095445115010332*(f[15]+f[12])+1.224744871391589*f[1])*alpha_vdim[17]+0.4898979485566357*alpha_cdim[6]*f[16]+1.224744871391589*(alpha_vdim[5]*f[16]+f[5]*alpha_vdim[16])+0.4898979485566357*f[6]*alpha_cdim[16]+1.224744871391589*(alpha_vdim[6]*f[15]+f[6]*alpha_vdim[15])+(0.5477225575051661*alpha_cdim[9]+0.3912303982179757*alpha_cdim[8]+1.224744871391589*alpha_vdim[4]+0.6123724356957944*alpha_cdim[0])*f[14]+1.224744871391589*f[4]*alpha_vdim[14]+(0.5477225575051661*f[9]+0.3912303982179757*f[8]+0.6123724356957944*f[0])*alpha_cdim[14]+1.224744871391589*(alpha_vdim[4]*f[13]+f[4]*alpha_vdim[13]+alpha_vdim[6]*f[12]+f[6]*alpha_vdim[12]+alpha_vdim[5]*f[11]+f[5]*alpha_vdim[11])+(1.224744871391589*(alpha_vdim[9]+alpha_vdim[8]+alpha_vdim[7])+1.369306393762915*alpha_vdim[0])*f[10]+(1.224744871391589*(f[9]+f[8]+f[7])+1.369306393762915*f[0])*alpha_vdim[10]+0.6123724356957944*(alpha_cdim[3]*f[8]+f[3]*alpha_cdim[8])+0.5477225575051661*alpha_cdim[2]*f[6]+1.369306393762915*(alpha_vdim[1]*f[6]+f[1]*alpha_vdim[6])+0.5477225575051661*f[2]*alpha_cdim[6]+1.369306393762915*(alpha_vdim[2]*f[5]+f[2]*alpha_vdim[5]+alpha_vdim[3]*f[4]+f[3]*alpha_vdim[4]); 
  out[19] += (1.095445115010332*(f[14]+f[13])+1.224744871391589*f[3])*alpha_vdim[39]+(1.095445115010332*(f[16]+f[11])+1.224744871391589*f[2])*alpha_vdim[38]+(1.095445115010332*(f[15]+f[12])+1.224744871391589*f[1])*alpha_vdim[37]+(1.095445115010332*f[18]+1.224744871391589*f[5])*alpha_vdim[36]+(1.095445115010332*f[17]+1.224744871391589*f[6])*alpha_vdim[35]+(1.095445115010332*f[19]+1.224744871391589*f[4])*alpha_vdim[34]+(1.095445115010332*f[19]+1.224744871391589*f[4])*alpha_vdim[33]+(1.095445115010332*f[17]+1.224744871391589*f[6])*alpha_vdim[32]+(1.095445115010332*f[18]+1.224744871391589*f[5])*alpha_vdim[31]+(1.224744871391589*(f[9]+f[8]+f[7])+1.369306393762915*f[0])*alpha_vdim[30]+1.224744871391589*f[10]*(alpha_vdim[29]+alpha_vdim[28]+alpha_vdim[27])+(1.224744871391589*(f[15]+f[12])+1.369306393762915*f[1])*alpha_vdim[26]+(1.224744871391589*(f[16]+f[11])+1.369306393762915*f[2])*alpha_vdim[25]+(1.224744871391589*(f[14]+f[13])+1.369306393762915*f[3])*alpha_vdim[24]+(1.224744871391589*f[19]+1.369306393762915*f[4])*alpha_vdim[23]+(1.224744871391589*f[18]+1.369306393762915*f[5])*alpha_vdim[22]+1.224744871391589*f[17]*alpha_vdim[21]+1.369306393762915*(f[6]*alpha_vdim[21]+f[10]*alpha_vdim[20])+(0.3912303982179757*alpha_vdim[16]+0.5477225575051661*alpha_vdim[11]+0.6123724356957944*alpha_vdim[2])*f[19]+(0.3912303982179757*f[16]+0.5477225575051661*f[11]+0.6123724356957944*f[2])*alpha_vdim[19]+0.5477225575051661*(alpha_vdim[14]*f[18]+f[14]*alpha_vdim[18])+0.4898979485566357*(alpha_vdim[10]*f[17]+f[10]*alpha_vdim[17])+(0.3912303982179757*alpha_cdim[9]+0.5477225575051661*alpha_cdim[8])*f[16]+0.6123724356957944*((alpha_vdim[4]+alpha_cdim[0])*f[16]+f[4]*alpha_vdim[16])+(0.3912303982179757*f[9]+0.5477225575051661*f[8]+0.6123724356957944*f[0])*alpha_cdim[16]+(0.3912303982179757*alpha_vdim[9]+0.5477225575051661*alpha_vdim[7]+0.6123724356957944*alpha_vdim[0])*f[15]+(0.3912303982179757*f[9]+0.5477225575051661*f[7]+0.6123724356957944*f[0])*alpha_vdim[15]+0.4898979485566357*(alpha_cdim[6]*f[14]+f[6]*alpha_cdim[14]+alpha_vdim[5]*f[13]+f[5]*alpha_vdim[13])+0.5477225575051661*(alpha_vdim[6]*f[10]+f[6]*alpha_vdim[10])+0.6123724356957944*((alpha_cdim[2]+alpha_vdim[1])*f[9]+f[1]*alpha_vdim[9]+f[2]*alpha_cdim[9])+0.5477225575051661*(alpha_cdim[3]*f[6]+f[3]*alpha_cdim[6]+alpha_vdim[3]*f[5]+f[3]*alpha_vdim[5]); 

  return alpha_mid; 
} 
