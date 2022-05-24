#include <gkyl_advection_kernels.h> 
GKYL_CU_DH double advection_vol_3x_ser_p2(const double *w, const double *dxv, const double *u, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // u[NDIM]:   Advection velocity.
  // f:         Input function.
  // out:       Incremented output.
  const double dx2 = 2.0/dxv[0]; 
  const double dy2 = 2.0/dxv[1]; 
  const double dz2 = 2.0/dxv[2]; 
  double alpha_mid = 0.0; 
  alpha_mid += fabs(0.1767766952966368*u[0]-0.1976423537605236*(u[9]+u[8]+u[7])); 
  alpha_mid += fabs(0.1767766952966368*u[20]-0.1976423537605236*(u[29]+u[28]+u[27])); 
  alpha_mid += fabs(0.1767766952966368*u[40]-0.1976423537605236*(u[49]+u[48]+u[47])); 

  out[1] += 0.6123724356957944*(f[19]*u[19]+f[18]*u[18]+f[17]*u[17]+f[16]*u[16]+f[15]*u[15]+f[14]*u[14]+f[13]*u[13]+f[12]*u[12]+f[11]*u[11]+f[10]*u[10]+f[9]*u[9]+f[8]*u[8]+f[7]*u[7]+f[6]*u[6]+f[5]*u[5]+f[4]*u[4]+f[3]*u[3]+f[2]*u[2]+f[1]*u[1]+f[0]*u[0])*dx2; 
  out[2] += 0.6123724356957944*(f[19]*u[39]+f[18]*u[38]+f[17]*u[37]+f[16]*u[36]+f[15]*u[35]+f[14]*u[34]+f[13]*u[33]+f[12]*u[32]+f[11]*u[31]+f[10]*u[30]+f[9]*u[29]+f[8]*u[28]+f[7]*u[27]+f[6]*u[26]+f[5]*u[25]+f[4]*u[24]+f[3]*u[23]+f[2]*u[22]+f[1]*u[21]+f[0]*u[20])*dy2; 
  out[3] += 0.6123724356957944*(f[19]*u[59]+f[18]*u[58]+f[17]*u[57]+f[16]*u[56]+f[15]*u[55]+f[14]*u[54]+f[13]*u[53]+f[12]*u[52]+f[11]*u[51]+f[10]*u[50]+f[9]*u[49]+f[8]*u[48]+f[7]*u[47]+f[6]*u[46]+f[5]*u[45]+f[4]*u[44]+f[3]*u[43]+f[2]*u[42]+f[1]*u[41]+f[0]*u[40])*dz2; 
  out[4] += (0.6123724356957944*(f[16]*u[39]+f[14]*u[38])+0.5477225575051661*f[10]*u[37]+0.6123724356957944*(f[19]*u[36]+f[9]*u[35]+f[18]*u[34])+0.5477225575051661*f[5]*u[33]+0.6123724356957944*f[8]*u[32]+0.5477225575051661*(f[4]*u[31]+f[17]*u[30]+f[1]*u[27])+0.6123724356957944*(f[6]*u[30]+f[15]*u[29]+f[12]*u[28]+f[10]*u[26])+(0.5477225575051661*f[13]+0.6123724356957944*f[3])*u[25]+0.5477225575051661*(f[11]*u[24]+f[7]*u[21])+0.6123724356957944*(f[2]*u[24]+f[5]*u[23]+f[4]*u[22]+f[0]*u[21]+f[1]*u[20]))*dy2+(0.6123724356957944*(f[15]*u[19]+u[15]*f[19])+0.5477225575051661*(f[10]*u[18]+u[10]*f[18])+0.6123724356957944*(f[13]*u[17]+u[13]*f[17]+f[9]*u[16]+u[9]*f[16])+0.5477225575051661*(f[6]*u[14]+u[6]*f[14]+f[4]*u[12]+u[4]*f[12])+0.6123724356957944*(f[7]*u[11]+u[7]*f[11]+f[5]*u[10]+u[5]*f[10])+0.5477225575051661*(f[2]*u[8]+u[2]*f[8])+0.6123724356957944*(f[3]*u[6]+u[3]*f[6]+f[1]*u[4]+u[1]*f[4]+f[0]*u[2]+u[0]*f[2]))*dx2; 
  out[5] += (0.6123724356957944*(f[16]*u[59]+f[14]*u[58])+0.5477225575051661*f[10]*u[57]+0.6123724356957944*(f[19]*u[56]+f[9]*u[55]+f[18]*u[54])+0.5477225575051661*f[5]*u[53]+0.6123724356957944*f[8]*u[52]+0.5477225575051661*(f[4]*u[51]+f[17]*u[50]+f[1]*u[47])+0.6123724356957944*(f[6]*u[50]+f[15]*u[49]+f[12]*u[48]+f[10]*u[46])+(0.5477225575051661*f[13]+0.6123724356957944*f[3])*u[45]+0.5477225575051661*(f[11]*u[44]+f[7]*u[41])+0.6123724356957944*(f[2]*u[44]+f[5]*u[43]+f[4]*u[42]+f[0]*u[41]+f[1]*u[40]))*dz2+(0.5477225575051661*(f[10]*u[19]+u[10]*f[19])+0.6123724356957944*(f[12]*u[18]+u[12]*f[18]+f[11]*u[17]+u[11]*f[17])+0.5477225575051661*(f[6]*u[16]+u[6]*f[16]+f[5]*u[15]+u[5]*f[15])+0.6123724356957944*(f[8]*u[14]+u[8]*f[14]+f[7]*u[13]+u[7]*f[13]+f[4]*u[10]+u[4]*f[10])+0.5477225575051661*(f[3]*u[9]+u[3]*f[9])+0.6123724356957944*(f[2]*u[6]+u[2]*f[6]+f[1]*u[5]+u[1]*f[5]+f[0]*u[3]+u[0]*f[3]))*dx2; 
  out[6] += (0.6123724356957944*f[15]*u[59]+0.5477225575051661*f[10]*u[58]+0.6123724356957944*(f[13]*u[57]+f[9]*u[56]+f[19]*u[55])+0.5477225575051661*f[6]*u[54]+0.6123724356957944*f[17]*u[53]+0.5477225575051661*f[4]*u[52]+0.6123724356957944*f[7]*u[51]+0.5477225575051661*(f[18]*u[50]+f[2]*u[48])+0.6123724356957944*(f[5]*u[50]+f[16]*u[49]+f[11]*u[47])+0.5477225575051661*(f[14]*u[46]+f[12]*u[44]+f[8]*u[42])+0.6123724356957944*(f[3]*u[46]+f[10]*u[45]+f[1]*u[44]+f[6]*u[43]+f[0]*u[42]+f[4]*u[41]+f[2]*u[40]))*dz2+(0.5477225575051661*f[10]*u[39]+0.6123724356957944*(f[12]*u[38]+f[11]*u[37])+0.5477225575051661*(f[6]*u[36]+f[5]*u[35])+0.6123724356957944*(f[8]*u[34]+f[7]*u[33]+f[18]*u[32]+f[17]*u[31])+(0.5477225575051661*f[19]+0.6123724356957944*f[4])*u[30]+0.5477225575051661*f[3]*u[29]+0.6123724356957944*(f[14]*u[28]+f[13]*u[27])+(0.5477225575051661*f[16]+0.6123724356957944*f[2])*u[26]+0.5477225575051661*(f[15]*u[25]+f[9]*u[23])+0.6123724356957944*(f[1]*u[25]+f[10]*u[24]+f[0]*u[23]+f[6]*u[22]+f[5]*u[21]+f[3]*u[20]))*dy2; 
  out[7] += (1.369306393762915*(f[16]*u[19]+u[16]*f[19]+f[14]*u[18]+u[14]*f[18])+1.224744871391589*(f[10]*u[17]+u[10]*f[17])+1.369306393762915*(f[9]*u[15]+u[9]*f[15])+1.224744871391589*(f[5]*u[13]+u[5]*f[13])+1.369306393762915*(f[8]*u[12]+u[8]*f[12])+1.224744871391589*(f[4]*u[11]+u[4]*f[11])+1.369306393762915*(f[6]*u[10]+u[6]*f[10])+1.224744871391589*(f[1]*u[7]+u[1]*f[7])+1.369306393762915*(f[3]*u[5]+u[3]*f[5]+f[2]*u[4]+u[2]*f[4]+f[0]*u[1]+u[0]*f[1]))*dx2; 
  out[8] += (1.369306393762915*f[15]*u[39]+1.224744871391589*f[10]*u[38]+1.369306393762915*(f[13]*u[37]+f[9]*u[36]+f[19]*u[35])+1.224744871391589*f[6]*u[34]+1.369306393762915*f[17]*u[33]+1.224744871391589*f[4]*u[32]+1.369306393762915*f[7]*u[31]+1.224744871391589*f[18]*u[30]+1.369306393762915*(f[5]*u[30]+f[16]*u[29])+1.224744871391589*f[2]*u[28]+1.369306393762915*f[11]*u[27]+1.224744871391589*f[14]*u[26]+1.369306393762915*(f[3]*u[26]+f[10]*u[25])+1.224744871391589*f[12]*u[24]+1.369306393762915*(f[1]*u[24]+f[6]*u[23])+1.224744871391589*f[8]*u[22]+1.369306393762915*(f[0]*u[22]+f[4]*u[21]+f[2]*u[20]))*dy2; 
  out[9] += (1.224744871391589*f[10]*u[59]+1.369306393762915*(f[12]*u[58]+f[11]*u[57])+1.224744871391589*(f[6]*u[56]+f[5]*u[55])+1.369306393762915*(f[8]*u[54]+f[7]*u[53]+f[18]*u[52]+f[17]*u[51])+(1.224744871391589*f[19]+1.369306393762915*f[4])*u[50]+1.224744871391589*f[3]*u[49]+1.369306393762915*(f[14]*u[48]+f[13]*u[47])+(1.224744871391589*f[16]+1.369306393762915*f[2])*u[46]+1.224744871391589*f[15]*u[45]+1.369306393762915*(f[1]*u[45]+f[10]*u[44])+1.224744871391589*f[9]*u[43]+1.369306393762915*(f[0]*u[43]+f[6]*u[42]+f[5]*u[41]+f[3]*u[40]))*dz2; 
  out[10] += (0.6123724356957944*f[9]*u[59]+(0.4898979485566357*f[17]+0.5477225575051661*f[6])*u[58]+(0.4898979485566357*f[18]+0.5477225575051661*f[5])*u[57]+0.6123724356957944*(f[15]*u[56]+f[16]*u[55])+0.5477225575051661*f[10]*(u[54]+u[53])+0.4898979485566357*(f[11]*u[52]+f[12]*u[51])+0.5477225575051661*(f[2]*u[52]+f[1]*u[51]+f[13]*u[50]+f[4]*(u[48]+u[47]))+0.6123724356957944*(f[3]*u[50]+f[19]*u[49])+0.5477225575051661*f[14]*u[50]+(0.5477225575051661*f[18]+0.6123724356957944*f[5])*u[46]+(0.5477225575051661*f[17]+0.6123724356957944*f[6])*u[45]+0.5477225575051661*(f[7]*u[44]+f[12]*u[42]+f[11]*u[41])+0.6123724356957944*(f[0]*u[44]+f[10]*u[43])+0.5477225575051661*f[8]*u[44]+0.6123724356957944*(f[1]*u[42]+f[2]*u[41]+f[4]*u[40]))*dz2+((0.4898979485566357*f[17]+0.5477225575051661*f[6])*u[39]+0.6123724356957944*f[8]*u[38]+0.5477225575051661*(f[4]*u[37]+f[10]*u[36])+0.4898979485566357*f[19]*u[37]+(0.4898979485566357*f[13]+0.5477225575051661*f[3])*u[35]+0.6123724356957944*f[12]*u[34]+(0.4898979485566357*f[15]+0.5477225575051661*f[1])*u[33]+0.6123724356957944*f[14]*u[32]+0.5477225575051661*f[10]*u[31]+(0.5477225575051661*(f[16]+f[11])+0.6123724356957944*f[2])*u[30]+0.5477225575051661*f[5]*u[29]+0.6123724356957944*f[18]*u[28]+0.5477225575051661*f[5]*u[27]+(0.5477225575051661*f[19]+0.6123724356957944*f[4])*u[26]+(0.5477225575051661*(f[9]+f[7])+0.6123724356957944*f[0])*u[25]+(0.5477225575051661*f[17]+0.6123724356957944*f[6])*u[24]+0.5477225575051661*(f[15]*u[23]+f[13]*u[21])+0.6123724356957944*(f[1]*u[23]+f[10]*u[22]+f[3]*u[21]+f[5]*u[20]))*dy2+(0.5477225575051661*(f[5]*u[19]+u[5]*f[19]+f[4]*u[18]+u[4]*f[18])+0.4898979485566357*(f[18]*u[19]+u[18]*f[19])+0.6123724356957944*(f[7]*u[17]+u[7]*f[17])+0.4898979485566357*(f[14]*u[16]+u[14]*f[16])+0.5477225575051661*(f[3]*u[16]+u[3]*f[16]+f[10]*u[15]+u[10]*f[15]+f[2]*u[14]+u[2]*f[14])+0.6123724356957944*(f[11]*u[13]+u[11]*f[13])+0.5477225575051661*(f[10]*u[12]+u[10]*f[12])+0.6123724356957944*(f[1]*u[10]+u[1]*f[10])+0.5477225575051661*(f[6]*u[9]+u[6]*f[9]+f[6]*u[8]+u[6]*f[8])+0.6123724356957944*(f[0]*u[6]+u[0]*f[6]+f[4]*u[5]+u[4]*f[5]+f[2]*u[3]+u[2]*f[3]))*dx2; 
  out[11] += (0.5477225575051661*(f[19]*u[39]+f[18]*u[38])+(0.3912303982179757*f[17]+0.6123724356957944*f[6])*u[37]+0.5477225575051661*f[15]*u[35]+(0.3912303982179757*f[13]+0.6123724356957944*f[3])*u[33]+0.5477225575051661*f[12]*u[32]+(0.3912303982179757*f[11]+0.6123724356957944*f[2])*u[31]+0.5477225575051661*f[10]*u[30]+0.6123724356957944*(f[0]*u[27]+f[17]*u[26])+0.3912303982179757*f[7]*u[27]+0.5477225575051661*(f[5]*u[25]+f[4]*u[24])+0.6123724356957944*(f[13]*u[23]+f[11]*u[22])+0.5477225575051661*f[1]*u[21]+0.6123724356957944*f[7]*u[20])*dy2+(1.369306393762915*(f[9]*u[19]+u[9]*f[19])+1.095445115010332*(f[17]*u[18]+u[17]*f[18])+1.224744871391589*(f[6]*u[18]+u[6]*f[18]+f[5]*u[17]+u[5]*f[17])+1.369306393762915*(f[15]*u[16]+u[15]*f[16])+1.224744871391589*(f[10]*u[14]+u[10]*f[14]+f[10]*u[13]+u[10]*f[13])+1.095445115010332*(f[11]*u[12]+u[11]*f[12])+1.224744871391589*(f[2]*u[12]+u[2]*f[12]+f[1]*u[11]+u[1]*f[11])+1.369306393762915*(f[3]*u[10]+u[3]*f[10])+1.224744871391589*(f[4]*u[8]+u[4]*f[8]+f[4]*u[7]+u[4]*f[7])+1.369306393762915*(f[5]*u[6]+u[5]*f[6]+f[0]*u[4]+u[0]*f[4]+f[1]*u[2]+u[1]*f[2]))*dx2; 
  out[12] += (1.369306393762915*f[9]*u[39]+(1.095445115010332*f[17]+1.224744871391589*f[6])*u[38]+(1.095445115010332*f[18]+1.224744871391589*f[5])*u[37]+1.369306393762915*(f[15]*u[36]+f[16]*u[35])+1.224744871391589*f[10]*(u[34]+u[33])+1.095445115010332*(f[11]*u[32]+f[12]*u[31])+1.224744871391589*(f[2]*u[32]+f[1]*u[31]+f[13]*u[30]+f[4]*(u[28]+u[27]))+1.369306393762915*(f[3]*u[30]+f[19]*u[29])+1.224744871391589*f[14]*u[30]+(1.224744871391589*f[18]+1.369306393762915*f[5])*u[26]+(1.224744871391589*f[17]+1.369306393762915*f[6])*u[25]+1.224744871391589*(f[7]*u[24]+f[12]*u[22]+f[11]*u[21])+1.369306393762915*(f[0]*u[24]+f[10]*u[23])+1.224744871391589*f[8]*u[24]+1.369306393762915*(f[1]*u[22]+f[2]*u[21]+f[4]*u[20]))*dy2+(0.5477225575051661*f[19]*u[19]+0.6123724356957944*(f[5]*u[18]+u[5]*f[18])+0.3912303982179757*f[18]*u[18]+0.5477225575051661*(f[17]*u[17]+f[16]*u[16])+0.3912303982179757*(f[14]*u[14]+f[12]*u[12])+0.6123724356957944*(f[3]*u[14]+u[3]*f[14]+f[1]*u[12]+u[1]*f[12])+0.5477225575051661*(f[11]*u[11]+f[10]*u[10])+0.6123724356957944*(f[0]*u[8]+u[0]*f[8])+0.3912303982179757*f[8]*u[8]+0.5477225575051661*(f[6]*u[6]+f[4]*u[4]+f[2]*u[2]))*dx2; 
  out[13] += (0.5477225575051661*(f[19]*u[59]+f[18]*u[58])+(0.3912303982179757*f[17]+0.6123724356957944*f[6])*u[57]+0.5477225575051661*f[15]*u[55]+(0.3912303982179757*f[13]+0.6123724356957944*f[3])*u[53]+0.5477225575051661*f[12]*u[52]+(0.3912303982179757*f[11]+0.6123724356957944*f[2])*u[51]+0.5477225575051661*f[10]*u[50]+0.6123724356957944*(f[0]*u[47]+f[17]*u[46])+0.3912303982179757*f[7]*u[47]+0.5477225575051661*(f[5]*u[45]+f[4]*u[44])+0.6123724356957944*(f[13]*u[43]+f[11]*u[42])+0.5477225575051661*f[1]*u[41]+0.6123724356957944*f[7]*u[40])*dz2+((1.095445115010332*f[17]+1.224744871391589*f[6])*u[19]+(1.095445115010332*u[17]+1.224744871391589*u[6])*f[19]+1.369306393762915*(f[8]*u[18]+u[8]*f[18])+1.224744871391589*(f[4]*u[17]+u[4]*f[17]+f[10]*u[16]+u[10]*f[16])+(1.095445115010332*f[13]+1.224744871391589*f[3])*u[15]+(1.095445115010332*u[13]+1.224744871391589*u[3])*f[15]+1.369306393762915*(f[12]*u[14]+u[12]*f[14])+1.224744871391589*(f[1]*u[13]+u[1]*f[13]+f[10]*u[11]+u[10]*f[11])+1.369306393762915*(f[2]*u[10]+u[2]*f[10])+1.224744871391589*(f[5]*u[9]+u[5]*f[9]+f[5]*u[7]+u[5]*f[7])+1.369306393762915*(f[4]*u[6]+u[4]*f[6]+f[0]*u[5]+u[0]*f[5]+f[1]*u[3]+u[1]*f[3]))*dx2; 
  out[14] += (0.5477225575051661*f[19]*u[59]+(0.3912303982179757*f[18]+0.6123724356957944*f[5])*u[58]+0.5477225575051661*(f[17]*u[57]+f[16]*u[56])+(0.3912303982179757*f[14]+0.6123724356957944*f[3])*u[54]+(0.3912303982179757*f[12]+0.6123724356957944*f[1])*u[52]+0.5477225575051661*(f[11]*u[51]+f[10]*u[50])+(0.3912303982179757*f[8]+0.6123724356957944*f[0])*u[48]+0.5477225575051661*f[6]*u[46]+0.6123724356957944*f[18]*u[45]+0.5477225575051661*f[4]*u[44]+0.6123724356957944*f[14]*u[43]+0.5477225575051661*f[2]*u[42]+0.6123724356957944*(f[12]*u[41]+f[8]*u[40]))*dz2+((1.095445115010332*f[18]+1.224744871391589*f[5])*u[39]+(1.095445115010332*f[19]+1.224744871391589*f[4])*u[38]+1.369306393762915*f[7]*u[37]+1.224744871391589*(f[3]*u[36]+f[10]*u[35])+1.095445115010332*f[14]*u[36]+(1.095445115010332*f[16]+1.224744871391589*f[2])*u[34]+1.369306393762915*f[11]*u[33]+1.224744871391589*f[10]*u[32]+1.369306393762915*f[13]*u[31]+(1.224744871391589*(f[15]+f[12])+1.369306393762915*f[1])*u[30]+1.224744871391589*f[6]*(u[29]+u[28])+1.369306393762915*f[17]*u[27]+(1.224744871391589*(f[9]+f[8])+1.369306393762915*f[0])*u[26]+(1.224744871391589*f[19]+1.369306393762915*f[4])*u[25]+(1.224744871391589*f[18]+1.369306393762915*f[5])*u[24]+1.224744871391589*(f[16]*u[23]+f[14]*u[22])+1.369306393762915*(f[2]*u[23]+f[3]*u[22]+f[10]*u[21]+f[6]*u[20]))*dy2; 
  out[15] += ((1.095445115010332*f[17]+1.224744871391589*f[6])*u[59]+1.369306393762915*f[8]*u[58]+1.224744871391589*(f[4]*u[57]+f[10]*u[56])+1.095445115010332*f[19]*u[57]+(1.095445115010332*f[13]+1.224744871391589*f[3])*u[55]+1.369306393762915*f[12]*u[54]+(1.095445115010332*f[15]+1.224744871391589*f[1])*u[53]+1.369306393762915*f[14]*u[52]+1.224744871391589*f[10]*u[51]+(1.224744871391589*(f[16]+f[11])+1.369306393762915*f[2])*u[50]+1.224744871391589*f[5]*u[49]+1.369306393762915*f[18]*u[48]+1.224744871391589*f[5]*u[47]+(1.224744871391589*f[19]+1.369306393762915*f[4])*u[46]+(1.224744871391589*(f[9]+f[7])+1.369306393762915*f[0])*u[45]+(1.224744871391589*f[17]+1.369306393762915*f[6])*u[44]+1.224744871391589*(f[15]*u[43]+f[13]*u[41])+1.369306393762915*(f[1]*u[43]+f[10]*u[42]+f[3]*u[41]+f[5]*u[40]))*dz2+(0.6123724356957944*(f[4]*u[19]+u[4]*f[19])+0.3912303982179757*f[19]*u[19]+0.5477225575051661*(f[18]*u[18]+f[17]*u[17])+0.3912303982179757*(f[16]*u[16]+f[15]*u[15])+0.6123724356957944*(f[2]*u[16]+u[2]*f[16]+f[1]*u[15]+u[1]*f[15])+0.5477225575051661*(f[14]*u[14]+f[13]*u[13]+f[10]*u[10])+0.6123724356957944*(f[0]*u[9]+u[0]*f[9])+0.3912303982179757*f[9]*u[9]+0.5477225575051661*(f[6]*u[6]+f[5]*u[5]+f[3]*u[3]))*dx2; 
  out[16] += ((1.095445115010332*f[18]+1.224744871391589*f[5])*u[59]+(1.095445115010332*f[19]+1.224744871391589*f[4])*u[58]+1.369306393762915*f[7]*u[57]+1.224744871391589*(f[3]*u[56]+f[10]*u[55])+1.095445115010332*f[14]*u[56]+(1.095445115010332*f[16]+1.224744871391589*f[2])*u[54]+1.369306393762915*f[11]*u[53]+1.224744871391589*f[10]*u[52]+1.369306393762915*f[13]*u[51]+(1.224744871391589*(f[15]+f[12])+1.369306393762915*f[1])*u[50]+1.224744871391589*f[6]*(u[49]+u[48])+1.369306393762915*f[17]*u[47]+(1.224744871391589*(f[9]+f[8])+1.369306393762915*f[0])*u[46]+(1.224744871391589*f[19]+1.369306393762915*f[4])*u[45]+(1.224744871391589*f[18]+1.369306393762915*f[5])*u[44]+1.224744871391589*(f[16]*u[43]+f[14]*u[42])+1.369306393762915*(f[2]*u[43]+f[3]*u[42]+f[10]*u[41]+f[6]*u[40]))*dz2+((0.3912303982179757*f[19]+0.6123724356957944*f[4])*u[39]+0.5477225575051661*(f[18]*u[38]+f[17]*u[37])+(0.3912303982179757*f[16]+0.6123724356957944*f[2])*u[36]+(0.3912303982179757*f[15]+0.6123724356957944*f[1])*u[35]+0.5477225575051661*(f[14]*u[34]+f[13]*u[33]+f[10]*u[30])+(0.3912303982179757*f[9]+0.6123724356957944*f[0])*u[29]+0.5477225575051661*(f[6]*u[26]+f[5]*u[25])+0.6123724356957944*f[19]*u[24]+0.5477225575051661*f[3]*u[23]+0.6123724356957944*(f[16]*u[22]+f[15]*u[21]+f[9]*u[20]))*dy2; 
  out[17] += (0.5477225575051661*f[15]*u[59]+0.4898979485566357*f[10]*u[58]+(0.5477225575051661*f[14]+0.3912303982179757*f[13]+0.6123724356957944*f[3])*u[57]+0.5477225575051661*(f[19]*u[55]+f[17]*u[54])+(0.3912303982179757*f[17]+0.6123724356957944*f[6])*u[53]+0.4898979485566356*f[4]*u[52]+(0.5477225575051661*f[8]+0.3912303982179757*f[7]+0.6123724356957944*f[0])*u[51]+0.5477225575051661*(f[5]*u[50]+f[11]*u[48])+0.4898979485566357*f[18]*u[50]+0.6123724356957944*(f[2]*u[47]+f[13]*u[46])+0.3912303982179757*f[11]*u[47]+0.5477225575051661*f[10]*u[45]+(0.4898979485566356*f[12]+0.5477225575051661*f[1])*u[44]+0.6123724356957944*(f[17]*u[43]+f[7]*u[42])+0.5477225575051661*f[4]*u[41]+0.6123724356957944*f[11]*u[40])*dz2+(0.4898979485566357*f[10]*u[39]+0.5477225575051661*f[12]*u[38]+(0.5477225575051661*f[16]+0.3912303982179757*f[11]+0.6123724356957944*f[2])*u[37]+0.5477225575051661*f[17]*u[36]+0.4898979485566356*f[5]*u[35]+(0.5477225575051661*f[9]+0.3912303982179757*f[7]+0.6123724356957944*f[0])*u[33]+0.5477225575051661*f[18]*u[32]+(0.3912303982179757*f[17]+0.6123724356957944*f[6])*u[31]+0.5477225575051661*(f[4]*u[30]+f[13]*u[29])+0.4898979485566357*f[19]*u[30]+0.6123724356957944*(f[3]*u[27]+f[11]*u[26])+0.3912303982179757*f[13]*u[27]+0.5477225575051661*(f[1]*u[25]+f[10]*u[24])+0.4898979485566356*f[15]*u[25]+0.6123724356957944*(f[7]*u[23]+f[17]*u[22])+0.5477225575051661*f[5]*u[21]+0.6123724356957944*f[13]*u[20])*dy2+((1.095445115010332*(f[14]+f[13])+1.224744871391589*f[3])*u[19]+(1.095445115010332*(u[14]+u[13])+1.224744871391589*u[3])*f[19]+(1.095445115010332*(f[16]+f[11])+1.224744871391589*f[2])*u[18]+(1.095445115010332*(u[16]+u[11])+1.224744871391589*u[2])*f[18]+1.095445115010332*(f[12]*u[17]+(u[15]+u[12])*f[17])+1.224744871391589*(f[1]*u[17]+u[1]*f[17]+f[5]*u[16]+u[5]*f[16]+f[6]*u[15]+u[6]*f[15]+f[4]*u[14]+u[4]*f[14]+f[4]*u[13]+u[4]*f[13]+f[6]*u[12]+u[6]*f[12]+f[5]*u[11]+u[5]*f[11]+f[7]*u[10]+(u[9]+u[8]+u[7])*f[10])+1.095445115010332*f[15]*u[17]+1.369306393762915*(f[0]*u[10]+u[0]*f[10]+f[1]*u[6]+u[1]*f[6]+f[2]*u[5]+u[2]*f[5]+f[3]*u[4]+u[3]*f[4])+1.224744871391589*(f[9]+f[8])*u[10])*dx2; 
  out[18] += (0.5477225575051661*f[16]*u[59]+(0.3912303982179757*f[14]+0.5477225575051661*f[13]+0.6123724356957944*f[3])*u[58]+0.4898979485566357*f[10]*u[57]+0.5477225575051661*f[19]*u[56]+(0.3912303982179757*f[18]+0.6123724356957944*f[5])*u[54]+0.5477225575051661*f[18]*u[53]+(0.3912303982179757*f[8]+0.5477225575051661*f[7]+0.6123724356957944*f[0])*u[52]+0.4898979485566356*f[4]*u[51]+(0.4898979485566357*f[17]+0.5477225575051661*f[6])*u[50]+(0.3912303982179757*f[12]+0.6123724356957944*f[1])*u[48]+0.5477225575051661*(f[12]*u[47]+f[10]*u[46])+0.6123724356957944*f[14]*u[45]+(0.4898979485566356*f[11]+0.5477225575051661*f[2])*u[44]+0.6123724356957944*f[18]*u[43]+0.5477225575051661*f[4]*u[42]+0.6123724356957944*(f[8]*u[41]+f[12]*u[40]))*dz2+((1.095445115010332*(f[14]+f[13])+1.224744871391589*f[3])*u[39]+(1.095445115010332*(f[16]+f[11])+1.224744871391589*f[2])*u[38]+(1.095445115010332*(f[15]+f[12])+1.224744871391589*f[1])*u[37]+(1.095445115010332*f[18]+1.224744871391589*f[5])*u[36]+(1.095445115010332*f[17]+1.224744871391589*f[6])*u[35]+(1.095445115010332*f[19]+1.224744871391589*f[4])*(u[34]+u[33])+(1.095445115010332*f[17]+1.224744871391589*f[6])*u[32]+(1.095445115010332*f[18]+1.224744871391589*f[5])*u[31]+(1.224744871391589*(f[9]+f[8]+f[7])+1.369306393762915*f[0])*u[30]+1.224744871391589*f[10]*(u[29]+u[28]+u[27])+(1.224744871391589*(f[15]+f[12])+1.369306393762915*f[1])*u[26]+(1.224744871391589*(f[16]+f[11])+1.369306393762915*f[2])*u[25]+(1.224744871391589*(f[14]+f[13])+1.369306393762915*f[3])*u[24]+(1.224744871391589*f[19]+1.369306393762915*f[4])*u[23]+1.224744871391589*(f[18]*u[22]+f[17]*u[21])+1.369306393762915*(f[5]*u[22]+f[6]*u[21]+f[10]*u[20]))*dy2+(0.4898979485566357*(f[10]*u[19]+u[10]*f[19])+(0.5477225575051661*f[15]+0.3912303982179757*f[12]+0.6123724356957944*f[1])*u[18]+(0.5477225575051661*u[15]+0.3912303982179757*u[12]+0.6123724356957944*u[1])*f[18]+0.5477225575051661*(f[11]*u[17]+u[11]*f[17])+0.4898979485566356*(f[6]*u[16]+u[6]*f[16])+0.6123724356957944*(f[0]*u[14]+u[0]*f[14]+f[5]*u[12]+u[5]*f[12])+(0.5477225575051661*f[9]+0.3912303982179757*f[8])*u[14]+(0.5477225575051661*u[9]+0.3912303982179757*u[8])*f[14]+0.5477225575051661*(f[4]*u[10]+u[4]*f[10])+0.6123724356957944*(f[3]*u[8]+u[3]*f[8])+0.5477225575051661*(f[2]*u[6]+u[2]*f[6]))*dx2; 
  out[19] += ((1.095445115010332*(f[14]+f[13])+1.224744871391589*f[3])*u[59]+(1.095445115010332*(f[16]+f[11])+1.224744871391589*f[2])*u[58]+(1.095445115010332*(f[15]+f[12])+1.224744871391589*f[1])*u[57]+(1.095445115010332*f[18]+1.224744871391589*f[5])*u[56]+(1.095445115010332*f[17]+1.224744871391589*f[6])*u[55]+(1.095445115010332*f[19]+1.224744871391589*f[4])*(u[54]+u[53])+(1.095445115010332*f[17]+1.224744871391589*f[6])*u[52]+(1.095445115010332*f[18]+1.224744871391589*f[5])*u[51]+(1.224744871391589*(f[9]+f[8]+f[7])+1.369306393762915*f[0])*u[50]+1.224744871391589*f[10]*(u[49]+u[48]+u[47])+(1.224744871391589*(f[15]+f[12])+1.369306393762915*f[1])*u[46]+(1.224744871391589*(f[16]+f[11])+1.369306393762915*f[2])*u[45]+(1.224744871391589*(f[14]+f[13])+1.369306393762915*f[3])*u[44]+(1.224744871391589*f[19]+1.369306393762915*f[4])*u[43]+1.224744871391589*(f[18]*u[42]+f[17]*u[41])+1.369306393762915*(f[5]*u[42]+f[6]*u[41]+f[10]*u[40]))*dz2+((0.3912303982179757*f[16]+0.5477225575051661*f[11]+0.6123724356957944*f[2])*u[39]+0.5477225575051661*f[14]*u[38]+0.4898979485566357*f[10]*u[37]+(0.3912303982179757*f[19]+0.6123724356957944*f[4])*u[36]+(0.3912303982179757*f[9]+0.5477225575051661*f[7]+0.6123724356957944*f[0])*u[35]+0.5477225575051661*f[18]*u[34]+0.4898979485566356*f[5]*u[33]+0.5477225575051661*f[19]*u[31]+(0.4898979485566357*f[17]+0.5477225575051661*f[6])*u[30]+(0.3912303982179757*f[15]+0.6123724356957944*f[1])*u[29]+0.5477225575051661*(f[15]*u[27]+f[10]*u[26])+(0.4898979485566356*f[13]+0.5477225575051661*f[3])*u[25]+0.6123724356957944*f[16]*u[24]+0.5477225575051661*f[5]*u[23]+0.6123724356957944*(f[19]*u[22]+f[9]*u[21]+f[15]*u[20]))*dy2+((0.3912303982179757*f[15]+0.5477225575051661*f[12]+0.6123724356957944*f[1])*u[19]+(0.3912303982179757*u[15]+0.5477225575051661*u[12]+0.6123724356957944*u[1])*f[19]+0.4898979485566357*(f[10]*u[18]+u[10]*f[18])+0.5477225575051661*(f[13]*u[17]+u[13]*f[17])+0.6123724356957944*(f[0]*u[16]+u[0]*f[16]+f[4]*u[15]+u[4]*f[15])+(0.3912303982179757*f[9]+0.5477225575051661*f[8])*u[16]+(0.3912303982179757*u[9]+0.5477225575051661*u[8])*f[16]+0.4898979485566356*(f[6]*u[14]+u[6]*f[14])+0.5477225575051661*(f[5]*u[10]+u[5]*f[10])+0.6123724356957944*(f[2]*u[9]+u[2]*f[9])+0.5477225575051661*(f[3]*u[6]+u[3]*f[6]))*dx2; 

  return alpha_mid; 
} 
