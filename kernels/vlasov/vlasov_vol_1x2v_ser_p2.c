#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_vol_1x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // field:      q/m*EM fields.
  // cot_vec:   Only used in gen geo.
  // f:         Input distribution function.
  // out:       Incremented output.
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  const double dv10 = 2/dxv[1]; 
  const double *E0 = &field[0]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv11 = 2/dxv[2]; 
  const double *E1 = &field[3]; 
  const double dv2 = dxv[2], wv2 = w[2]; 

  const double *B2 = &field[15]; 
  double cflFreq_mid = 0.0; 
  double alpha_cdim[20] = {0.0}; 
  double alpha_vdim[40] = {0.0}; 

  alpha_cdim[0] = 5.656854249492382*w0dx0; 
  alpha_cdim[2] = 1.632993161855453*dv0dx0; 
  cflFreq_mid += 5.0*(fabs(w0dx0)+0.5*dv0dx0); 

  alpha_vdim[0] = 2.0*dv10*(B2[0]*wv2+E0[0]); 
  alpha_vdim[1] = 2.0*dv10*(B2[1]*wv2+E0[1]); 
  alpha_vdim[3] = 0.5773502691896258*B2[0]*dv10*dv2; 
  alpha_vdim[5] = 0.5773502691896258*B2[1]*dv10*dv2; 
  alpha_vdim[7] = 2.0*dv10*(B2[2]*wv2+E0[2]); 
  alpha_vdim[13] = 0.5773502691896258*B2[2]*dv10*dv2; 
  cflFreq_mid += 5.0*fabs(0.1767766952966368*alpha_vdim[0]-0.1976423537605236*alpha_vdim[7]); 

  alpha_vdim[20] = dv11*(2.0*E1[0]-2.0*B2[0]*wv1); 
  alpha_vdim[21] = dv11*(2.0*E1[1]-2.0*B2[1]*wv1); 
  alpha_vdim[22] = -0.5773502691896258*B2[0]*dv1*dv11; 
  alpha_vdim[24] = -0.5773502691896258*B2[1]*dv1*dv11; 
  alpha_vdim[27] = dv11*(2.0*E1[2]-2.0*B2[2]*wv1); 
  alpha_vdim[31] = -0.5773502691896258*B2[2]*dv1*dv11; 
  cflFreq_mid += 5.0*fabs(0.1767766952966368*alpha_vdim[20]-0.1976423537605236*alpha_vdim[27]); 

  out[1] += 0.6123724356957944*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.6123724356957944*(alpha_vdim[13]*f[13]+alpha_vdim[7]*f[7]+alpha_vdim[5]*f[5]+alpha_vdim[3]*f[3]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[11]*alpha_vdim[31]+f[7]*alpha_vdim[27]+f[4]*alpha_vdim[24]+f[2]*alpha_vdim[22]+f[1]*alpha_vdim[21]+f[0]*alpha_vdim[20]); 
  out[4] += 0.5477225575051661*(alpha_vdim[5]*f[13]+f[5]*alpha_vdim[13]+alpha_cdim[2]*f[8]+alpha_vdim[1]*f[7]+f[1]*alpha_vdim[7])+0.6123724356957944*(alpha_vdim[3]*f[5]+f[3]*alpha_vdim[5]+alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[5] += 0.5477225575051661*(f[4]*alpha_vdim[31]+f[1]*alpha_vdim[27]+f[11]*alpha_vdim[24])+0.6123724356957944*(f[2]*alpha_vdim[24]+f[4]*alpha_vdim[22])+0.5477225575051661*f[7]*alpha_vdim[21]+0.6123724356957944*(f[0]*alpha_vdim[21]+f[1]*alpha_vdim[20]+alpha_cdim[2]*f[6]+alpha_cdim[0]*f[3]); 
  out[6] += 0.6123724356957944*(f[7]*alpha_vdim[31]+f[11]*alpha_vdim[27])+(0.5477225575051661*f[12]+0.6123724356957944*f[1])*alpha_vdim[24]+0.5477225575051661*f[8]*alpha_vdim[22]+0.6123724356957944*(f[0]*alpha_vdim[22]+f[4]*alpha_vdim[21]+f[2]*alpha_vdim[20])+0.5477225575051661*alpha_vdim[5]*f[15]+0.6123724356957944*(alpha_vdim[7]*f[13]+f[7]*alpha_vdim[13])+0.5477225575051661*alpha_vdim[3]*f[9]+0.6123724356957944*(alpha_vdim[1]*f[5]+f[1]*alpha_vdim[5]+alpha_vdim[0]*f[3]+f[0]*alpha_vdim[3]); 
  out[7] += 1.369306393762915*(alpha_cdim[2]*f[4]+alpha_cdim[0]*f[1]); 
  out[8] += 1.369306393762915*(alpha_vdim[13]*f[17]+alpha_vdim[7]*f[11]+alpha_vdim[5]*f[10]+alpha_vdim[3]*f[6]+alpha_vdim[1]*f[4]+alpha_vdim[0]*f[2]); 
  out[9] += 1.369306393762915*(f[17]*alpha_vdim[31]+f[13]*alpha_vdim[27]+f[10]*alpha_vdim[24]+f[6]*alpha_vdim[22]+f[5]*alpha_vdim[21]+f[3]*alpha_vdim[20]); 
  out[10] += 0.4898979485566357*f[12]*alpha_vdim[31]+0.5477225575051661*(f[1]*alpha_vdim[31]+f[4]*alpha_vdim[27])+(0.5477225575051661*(f[8]+f[7])+0.6123724356957944*f[0])*alpha_vdim[24]+(0.5477225575051661*f[12]+0.6123724356957944*f[1])*alpha_vdim[22]+0.5477225575051661*f[11]*alpha_vdim[21]+0.6123724356957944*(f[2]*alpha_vdim[21]+f[4]*alpha_vdim[20])+0.4898979485566357*alpha_vdim[13]*f[15]+0.5477225575051661*(alpha_vdim[3]*f[15]+alpha_cdim[2]*f[14]+alpha_vdim[1]*f[13]+f[1]*alpha_vdim[13]+alpha_vdim[5]*(f[9]+f[7])+f[5]*alpha_vdim[7])+0.6123724356957944*(alpha_cdim[0]*f[6]+alpha_vdim[0]*f[5]+f[0]*alpha_vdim[5]+(alpha_cdim[2]+alpha_vdim[1])*f[3]+f[1]*alpha_vdim[3]); 
  out[11] += 0.3912303982179757*alpha_vdim[13]*f[13]+0.6123724356957944*(alpha_vdim[3]*f[13]+f[3]*alpha_vdim[13])+1.224744871391589*alpha_cdim[2]*f[12]+0.3912303982179757*alpha_vdim[7]*f[7]+0.6123724356957944*(alpha_vdim[0]*f[7]+f[0]*alpha_vdim[7])+0.5477225575051661*alpha_vdim[5]*f[5]+1.369306393762915*alpha_cdim[0]*f[4]+f[1]*(1.369306393762915*alpha_cdim[2]+0.5477225575051661*alpha_vdim[1]); 
  out[12] += 1.224744871391589*(alpha_vdim[5]*f[17]+f[10]*alpha_vdim[13]+alpha_vdim[1]*f[11])+1.369306393762915*alpha_vdim[3]*f[10]+0.6123724356957944*alpha_cdim[0]*f[8]+1.224744871391589*f[4]*alpha_vdim[7]+1.369306393762915*(alpha_vdim[5]*f[6]+alpha_vdim[0]*f[4])+(0.5477225575051661*alpha_cdim[2]+1.369306393762915*alpha_vdim[1])*f[2]; 
  out[13] += (0.3912303982179757*f[11]+0.6123724356957944*f[2])*alpha_vdim[31]+(0.3912303982179757*f[7]+0.6123724356957944*f[0])*alpha_vdim[27]+0.5477225575051661*f[4]*alpha_vdim[24]+0.6123724356957944*f[11]*alpha_vdim[22]+0.5477225575051661*f[1]*alpha_vdim[21]+0.6123724356957944*f[7]*alpha_vdim[20]+1.369306393762915*(alpha_cdim[2]*f[10]+alpha_cdim[0]*f[5]); 
  out[14] += 0.5477225575051661*(f[11]*alpha_vdim[31]+f[4]*alpha_vdim[24]+f[2]*alpha_vdim[22])+0.6123724356957944*(f[12]*alpha_vdim[21]+f[8]*alpha_vdim[20])+1.224744871391589*alpha_vdim[5]*f[19]+1.369306393762915*alpha_vdim[7]*f[17]+1.224744871391589*alpha_vdim[3]*f[16]+1.369306393762915*(f[11]*alpha_vdim[13]+alpha_vdim[1]*f[10]+alpha_vdim[0]*f[6]+f[4]*alpha_vdim[5]+f[2]*alpha_vdim[3]); 
  out[15] += 1.224744871391589*(f[10]*alpha_vdim[31]+f[5]*alpha_vdim[27]+f[17]*alpha_vdim[24])+1.369306393762915*(f[6]*alpha_vdim[24]+f[10]*alpha_vdim[22])+1.224744871391589*f[13]*alpha_vdim[21]+1.369306393762915*(f[3]*alpha_vdim[21]+f[5]*alpha_vdim[20])+0.6123724356957944*(alpha_cdim[2]*f[16]+alpha_cdim[0]*f[9]); 
  out[16] += 1.369306393762915*(f[13]*alpha_vdim[31]+f[17]*alpha_vdim[27])+(1.224744871391589*f[18]+1.369306393762915*f[5])*alpha_vdim[24]+1.224744871391589*f[14]*alpha_vdim[22]+1.369306393762915*(f[3]*alpha_vdim[22]+f[10]*alpha_vdim[21]+f[6]*alpha_vdim[20])+0.6123724356957944*alpha_vdim[1]*f[15]+0.5477225575051661*alpha_vdim[13]*f[13]+0.6123724356957944*alpha_vdim[0]*f[9]+0.5477225575051661*(alpha_vdim[5]*f[5]+alpha_vdim[3]*f[3]); 
  out[17] += (0.5477225575051661*f[8]+0.3912303982179757*f[7]+0.6123724356957944*f[0])*alpha_vdim[31]+(0.3912303982179757*f[11]+0.6123724356957944*f[2])*alpha_vdim[27]+(0.4898979485566357*f[12]+0.5477225575051661*f[1])*alpha_vdim[24]+0.6123724356957944*f[7]*alpha_vdim[22]+0.5477225575051661*f[4]*alpha_vdim[21]+0.6123724356957944*f[11]*alpha_vdim[20]+1.224744871391589*alpha_cdim[2]*f[18]+0.4898979485566357*alpha_vdim[5]*f[15]+(0.3912303982179757*alpha_vdim[7]+0.6123724356957944*alpha_vdim[0])*f[13]+(0.5477225575051661*f[9]+0.3912303982179757*f[7]+0.6123724356957944*f[0])*alpha_vdim[13]+1.369306393762915*alpha_cdim[0]*f[10]+0.6123724356957944*(alpha_vdim[3]*f[7]+f[3]*alpha_vdim[7])+1.369306393762915*alpha_cdim[2]*f[5]+0.5477225575051661*(alpha_vdim[1]*f[5]+f[1]*alpha_vdim[5]); 
  out[18] += 0.4898979485566357*f[4]*alpha_vdim[31]+0.5477225575051661*f[12]*alpha_vdim[27]+0.4898979485566357*f[11]*alpha_vdim[24]+0.5477225575051661*(f[2]*alpha_vdim[24]+f[4]*alpha_vdim[22])+0.6123724356957944*(f[8]*alpha_vdim[21]+f[12]*alpha_vdim[20])+1.095445115010332*alpha_vdim[13]*f[19]+1.224744871391589*(alpha_vdim[3]*f[19]+alpha_vdim[1]*f[17]+alpha_vdim[5]*f[16])+0.6123724356957944*alpha_cdim[0]*f[14]+1.224744871391589*(f[4]*alpha_vdim[13]+alpha_vdim[5]*f[11])+(1.224744871391589*alpha_vdim[7]+1.369306393762915*alpha_vdim[0])*f[10]+0.5477225575051661*alpha_cdim[2]*f[6]+1.369306393762915*(alpha_vdim[1]*f[6]+f[2]*alpha_vdim[5]+alpha_vdim[3]*f[4]); 
  out[19] += 1.095445115010332*f[18]*alpha_vdim[31]+1.224744871391589*(f[5]*alpha_vdim[31]+f[10]*alpha_vdim[27])+(1.224744871391589*(f[14]+f[13])+1.369306393762915*f[3])*alpha_vdim[24]+(1.224744871391589*f[18]+1.369306393762915*f[5])*alpha_vdim[22]+1.224744871391589*f[17]*alpha_vdim[21]+1.369306393762915*(f[6]*alpha_vdim[21]+f[10]*alpha_vdim[20])+0.6123724356957944*alpha_cdim[0]*f[16]+(0.5477225575051661*alpha_vdim[7]+0.6123724356957944*alpha_vdim[0])*f[15]+0.4898979485566357*(alpha_vdim[5]*f[13]+f[5]*alpha_vdim[13])+0.6123724356957944*(alpha_cdim[2]+alpha_vdim[1])*f[9]+0.5477225575051661*(alpha_vdim[3]*f[5]+f[3]*alpha_vdim[5]); 

  return cflFreq_mid; 
} 
