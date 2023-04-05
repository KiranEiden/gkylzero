#include <gkyl_lbo_vlasov_pkpm_kernels.h> 
GKYL_CU_DH double lbo_vlasov_pkpm_drag_vol_2x1v_ser_p2(const double *w, const double *dxv, const double *nu, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[3]:   Cell-center coordinates. 
  // dxv[3]: Cell spacing. 
  // nu:      collisionality. 
  // f:       Input distribution function.
  // out:     Incremented output 
  const double rdvpar = 2.0/dxv[2]; 
  const double dvpar = dxv[2], wvpar = w[2]; 
  const double *F_0 = &f[0]; 
  const double *G_1 = &f[20]; 
  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[20]; 
  double alphaDrag[20]; 
  // Expand rdvpar*(nu*vx) in phase basis.
  alphaDrag[0] = -1.414213562373095*nu[0]*rdvpar*wvpar; 
  alphaDrag[1] = -1.414213562373095*nu[1]*rdvpar*wvpar; 
  alphaDrag[2] = -1.414213562373095*nu[2]*rdvpar*wvpar; 
  alphaDrag[3] = -0.408248290463863*nu[0]*dvpar*rdvpar; 
  alphaDrag[4] = -1.414213562373095*nu[3]*rdvpar*wvpar; 
  alphaDrag[5] = -0.408248290463863*nu[1]*dvpar*rdvpar; 
  alphaDrag[6] = -0.408248290463863*nu[2]*dvpar*rdvpar; 
  alphaDrag[7] = -1.414213562373095*nu[4]*rdvpar*wvpar; 
  alphaDrag[8] = -1.414213562373095*nu[5]*rdvpar*wvpar; 
  alphaDrag[10] = -0.408248290463863*nu[3]*dvpar*rdvpar; 
  alphaDrag[11] = -1.414213562373095*nu[6]*rdvpar*wvpar; 
  alphaDrag[12] = -1.414213562373095*nu[7]*rdvpar*wvpar; 
  alphaDrag[13] = -0.408248290463863*nu[4]*dvpar*rdvpar; 
  alphaDrag[14] = -0.408248290463863*nu[5]*dvpar*rdvpar; 
  alphaDrag[17] = -0.408248290463863*nu[6]*dvpar*rdvpar; 
  alphaDrag[18] = -0.408248290463863*nu[7]*dvpar*rdvpar; 

  out_F_0[3] += 0.6123724356957944*F_0[18]*alphaDrag[18]+0.6123724356957944*F_0[17]*alphaDrag[17]+0.6123724356957944*F_0[14]*alphaDrag[14]+0.6123724356957944*F_0[13]*alphaDrag[13]+0.6123724356957944*F_0[12]*alphaDrag[12]+0.6123724356957944*F_0[11]*alphaDrag[11]+0.6123724356957944*F_0[10]*alphaDrag[10]+0.6123724356957944*F_0[8]*alphaDrag[8]+0.6123724356957944*F_0[7]*alphaDrag[7]+0.6123724356957944*F_0[6]*alphaDrag[6]+0.6123724356957944*F_0[5]*alphaDrag[5]+0.6123724356957944*F_0[4]*alphaDrag[4]+0.6123724356957944*F_0[3]*alphaDrag[3]+0.6123724356957944*F_0[2]*alphaDrag[2]+0.6123724356957944*F_0[1]*alphaDrag[1]+0.6123724356957944*F_0[0]*alphaDrag[0]; 
  out_F_0[5] += 0.6123724356957944*F_0[14]*alphaDrag[18]+0.6123724356957944*alphaDrag[14]*F_0[18]+0.5477225575051661*F_0[10]*alphaDrag[17]+0.5477225575051661*alphaDrag[10]*F_0[17]+0.5477225575051661*F_0[5]*alphaDrag[13]+0.5477225575051661*alphaDrag[5]*F_0[13]+0.6123724356957944*F_0[8]*alphaDrag[12]+0.6123724356957944*alphaDrag[8]*F_0[12]+0.5477225575051661*F_0[4]*alphaDrag[11]+0.5477225575051661*alphaDrag[4]*F_0[11]+0.6123724356957944*F_0[6]*alphaDrag[10]+0.6123724356957944*alphaDrag[6]*F_0[10]+0.5477225575051661*F_0[1]*alphaDrag[7]+0.5477225575051661*alphaDrag[1]*F_0[7]+0.6123724356957944*F_0[3]*alphaDrag[5]+0.6123724356957944*alphaDrag[3]*F_0[5]+0.6123724356957944*F_0[2]*alphaDrag[4]+0.6123724356957944*alphaDrag[2]*F_0[4]+0.6123724356957944*F_0[0]*alphaDrag[1]+0.6123724356957944*alphaDrag[0]*F_0[1]; 
  out_F_0[6] += 0.5477225575051661*F_0[10]*alphaDrag[18]+0.5477225575051661*alphaDrag[10]*F_0[18]+0.6123724356957944*F_0[13]*alphaDrag[17]+0.6123724356957944*alphaDrag[13]*F_0[17]+0.5477225575051661*F_0[6]*alphaDrag[14]+0.5477225575051661*alphaDrag[6]*F_0[14]+0.5477225575051661*F_0[4]*alphaDrag[12]+0.5477225575051661*alphaDrag[4]*F_0[12]+0.6123724356957944*F_0[7]*alphaDrag[11]+0.6123724356957944*alphaDrag[7]*F_0[11]+0.6123724356957944*F_0[5]*alphaDrag[10]+0.6123724356957944*alphaDrag[5]*F_0[10]+0.5477225575051661*F_0[2]*alphaDrag[8]+0.5477225575051661*alphaDrag[2]*F_0[8]+0.6123724356957944*F_0[3]*alphaDrag[6]+0.6123724356957944*alphaDrag[3]*F_0[6]+0.6123724356957944*F_0[1]*alphaDrag[4]+0.6123724356957944*alphaDrag[1]*F_0[4]+0.6123724356957944*F_0[0]*alphaDrag[2]+0.6123724356957944*alphaDrag[0]*F_0[2]; 
  out_F_0[9] += 1.224744871391589*alphaDrag[10]*F_0[19]+1.369306393762915*F_0[12]*alphaDrag[18]+1.369306393762915*alphaDrag[12]*F_0[18]+1.369306393762915*F_0[11]*alphaDrag[17]+1.369306393762915*alphaDrag[11]*F_0[17]+1.224744871391589*alphaDrag[6]*F_0[16]+1.224744871391589*alphaDrag[5]*F_0[15]+1.369306393762915*F_0[8]*alphaDrag[14]+1.369306393762915*alphaDrag[8]*F_0[14]+1.369306393762915*F_0[7]*alphaDrag[13]+1.369306393762915*alphaDrag[7]*F_0[13]+1.369306393762915*F_0[4]*alphaDrag[10]+1.369306393762915*alphaDrag[4]*F_0[10]+1.224744871391589*alphaDrag[3]*F_0[9]+1.369306393762915*F_0[2]*alphaDrag[6]+1.369306393762915*alphaDrag[2]*F_0[6]+1.369306393762915*F_0[1]*alphaDrag[5]+1.369306393762915*alphaDrag[1]*F_0[5]+1.369306393762915*F_0[0]*alphaDrag[3]+1.369306393762915*alphaDrag[0]*F_0[3]; 
  out_F_0[10] += 0.4898979485566357*F_0[17]*alphaDrag[18]+0.5477225575051661*F_0[6]*alphaDrag[18]+0.4898979485566357*alphaDrag[17]*F_0[18]+0.5477225575051661*alphaDrag[6]*F_0[18]+0.5477225575051661*F_0[5]*alphaDrag[17]+0.5477225575051661*alphaDrag[5]*F_0[17]+0.5477225575051661*F_0[10]*alphaDrag[14]+0.5477225575051661*alphaDrag[10]*F_0[14]+0.5477225575051661*F_0[10]*alphaDrag[13]+0.5477225575051661*alphaDrag[10]*F_0[13]+0.4898979485566357*F_0[11]*alphaDrag[12]+0.5477225575051661*F_0[2]*alphaDrag[12]+0.4898979485566357*alphaDrag[11]*F_0[12]+0.5477225575051661*alphaDrag[2]*F_0[12]+0.5477225575051661*F_0[1]*alphaDrag[11]+0.5477225575051661*alphaDrag[1]*F_0[11]+0.6123724356957944*F_0[3]*alphaDrag[10]+0.6123724356957944*alphaDrag[3]*F_0[10]+0.5477225575051661*F_0[4]*alphaDrag[8]+0.5477225575051661*alphaDrag[4]*F_0[8]+0.5477225575051661*F_0[4]*alphaDrag[7]+0.5477225575051661*alphaDrag[4]*F_0[7]+0.6123724356957944*F_0[5]*alphaDrag[6]+0.6123724356957944*alphaDrag[5]*F_0[6]+0.6123724356957944*F_0[0]*alphaDrag[4]+0.6123724356957944*alphaDrag[0]*F_0[4]+0.6123724356957944*F_0[1]*alphaDrag[2]+0.6123724356957944*alphaDrag[1]*F_0[2]; 
  out_F_0[13] += 0.5477225575051661*F_0[18]*alphaDrag[18]+0.3912303982179757*F_0[17]*alphaDrag[17]+0.6123724356957944*F_0[6]*alphaDrag[17]+0.6123724356957944*alphaDrag[6]*F_0[17]+0.3912303982179757*F_0[13]*alphaDrag[13]+0.6123724356957944*F_0[3]*alphaDrag[13]+0.6123724356957944*alphaDrag[3]*F_0[13]+0.5477225575051661*F_0[12]*alphaDrag[12]+0.3912303982179757*F_0[11]*alphaDrag[11]+0.6123724356957944*F_0[2]*alphaDrag[11]+0.6123724356957944*alphaDrag[2]*F_0[11]+0.5477225575051661*F_0[10]*alphaDrag[10]+0.3912303982179757*F_0[7]*alphaDrag[7]+0.6123724356957944*F_0[0]*alphaDrag[7]+0.6123724356957944*alphaDrag[0]*F_0[7]+0.5477225575051661*F_0[5]*alphaDrag[5]+0.5477225575051661*F_0[4]*alphaDrag[4]+0.5477225575051661*F_0[1]*alphaDrag[1]; 
  out_F_0[14] += 0.3912303982179757*F_0[18]*alphaDrag[18]+0.6123724356957944*F_0[5]*alphaDrag[18]+0.6123724356957944*alphaDrag[5]*F_0[18]+0.5477225575051661*F_0[17]*alphaDrag[17]+0.3912303982179757*F_0[14]*alphaDrag[14]+0.6123724356957944*F_0[3]*alphaDrag[14]+0.6123724356957944*alphaDrag[3]*F_0[14]+0.3912303982179757*F_0[12]*alphaDrag[12]+0.6123724356957944*F_0[1]*alphaDrag[12]+0.6123724356957944*alphaDrag[1]*F_0[12]+0.5477225575051661*F_0[11]*alphaDrag[11]+0.5477225575051661*F_0[10]*alphaDrag[10]+0.3912303982179757*F_0[8]*alphaDrag[8]+0.6123724356957944*F_0[0]*alphaDrag[8]+0.6123724356957944*alphaDrag[0]*F_0[8]+0.5477225575051661*F_0[6]*alphaDrag[6]+0.5477225575051661*F_0[4]*alphaDrag[4]+0.5477225575051661*F_0[2]*alphaDrag[2]; 
  out_F_0[15] += 1.095445115010332*alphaDrag[17]*F_0[19]+1.224744871391589*alphaDrag[6]*F_0[19]+1.369306393762915*F_0[8]*alphaDrag[18]+1.369306393762915*alphaDrag[8]*F_0[18]+1.224744871391589*F_0[4]*alphaDrag[17]+1.224744871391589*alphaDrag[4]*F_0[17]+1.224744871391589*alphaDrag[10]*F_0[16]+1.095445115010332*alphaDrag[13]*F_0[15]+1.224744871391589*alphaDrag[3]*F_0[15]+1.369306393762915*F_0[12]*alphaDrag[14]+1.369306393762915*alphaDrag[12]*F_0[14]+1.224744871391589*F_0[1]*alphaDrag[13]+1.224744871391589*alphaDrag[1]*F_0[13]+1.224744871391589*F_0[10]*alphaDrag[11]+1.224744871391589*alphaDrag[10]*F_0[11]+1.369306393762915*F_0[2]*alphaDrag[10]+1.369306393762915*alphaDrag[2]*F_0[10]+1.224744871391589*alphaDrag[5]*F_0[9]+1.224744871391589*F_0[5]*alphaDrag[7]+1.224744871391589*alphaDrag[5]*F_0[7]+1.369306393762915*F_0[4]*alphaDrag[6]+1.369306393762915*alphaDrag[4]*F_0[6]+1.369306393762915*F_0[0]*alphaDrag[5]+1.369306393762915*alphaDrag[0]*F_0[5]+1.369306393762915*F_0[1]*alphaDrag[3]+1.369306393762915*alphaDrag[1]*F_0[3]; 
  out_F_0[16] += 1.095445115010332*alphaDrag[18]*F_0[19]+1.224744871391589*alphaDrag[5]*F_0[19]+1.224744871391589*F_0[4]*alphaDrag[18]+1.224744871391589*alphaDrag[4]*F_0[18]+1.369306393762915*F_0[7]*alphaDrag[17]+1.369306393762915*alphaDrag[7]*F_0[17]+1.095445115010332*alphaDrag[14]*F_0[16]+1.224744871391589*alphaDrag[3]*F_0[16]+1.224744871391589*alphaDrag[10]*F_0[15]+1.224744871391589*F_0[2]*alphaDrag[14]+1.224744871391589*alphaDrag[2]*F_0[14]+1.369306393762915*F_0[11]*alphaDrag[13]+1.369306393762915*alphaDrag[11]*F_0[13]+1.224744871391589*F_0[10]*alphaDrag[12]+1.224744871391589*alphaDrag[10]*F_0[12]+1.369306393762915*F_0[1]*alphaDrag[10]+1.369306393762915*alphaDrag[1]*F_0[10]+1.224744871391589*alphaDrag[6]*F_0[9]+1.224744871391589*F_0[6]*alphaDrag[8]+1.224744871391589*alphaDrag[6]*F_0[8]+1.369306393762915*F_0[0]*alphaDrag[6]+1.369306393762915*alphaDrag[0]*F_0[6]+1.369306393762915*F_0[4]*alphaDrag[5]+1.369306393762915*alphaDrag[4]*F_0[5]+1.369306393762915*F_0[2]*alphaDrag[3]+1.369306393762915*alphaDrag[2]*F_0[3]; 
  out_F_0[17] += 0.4898979485566357*F_0[10]*alphaDrag[18]+0.4898979485566357*alphaDrag[10]*F_0[18]+0.5477225575051661*F_0[14]*alphaDrag[17]+0.3912303982179757*F_0[13]*alphaDrag[17]+0.6123724356957944*F_0[3]*alphaDrag[17]+0.5477225575051661*alphaDrag[14]*F_0[17]+0.3912303982179757*alphaDrag[13]*F_0[17]+0.6123724356957944*alphaDrag[3]*F_0[17]+0.6123724356957944*F_0[6]*alphaDrag[13]+0.6123724356957944*alphaDrag[6]*F_0[13]+0.4898979485566356*F_0[4]*alphaDrag[12]+0.4898979485566356*alphaDrag[4]*F_0[12]+0.5477225575051661*F_0[8]*alphaDrag[11]+0.3912303982179757*F_0[7]*alphaDrag[11]+0.6123724356957944*F_0[0]*alphaDrag[11]+0.5477225575051661*alphaDrag[8]*F_0[11]+0.3912303982179757*alphaDrag[7]*F_0[11]+0.6123724356957944*alphaDrag[0]*F_0[11]+0.5477225575051661*F_0[5]*alphaDrag[10]+0.5477225575051661*alphaDrag[5]*F_0[10]+0.6123724356957944*F_0[2]*alphaDrag[7]+0.6123724356957944*alphaDrag[2]*F_0[7]+0.5477225575051661*F_0[1]*alphaDrag[4]+0.5477225575051661*alphaDrag[1]*F_0[4]; 
  out_F_0[18] += 0.3912303982179757*F_0[14]*alphaDrag[18]+0.5477225575051661*F_0[13]*alphaDrag[18]+0.6123724356957944*F_0[3]*alphaDrag[18]+0.3912303982179757*alphaDrag[14]*F_0[18]+0.5477225575051661*alphaDrag[13]*F_0[18]+0.6123724356957944*alphaDrag[3]*F_0[18]+0.4898979485566357*F_0[10]*alphaDrag[17]+0.4898979485566357*alphaDrag[10]*F_0[17]+0.6123724356957944*F_0[5]*alphaDrag[14]+0.6123724356957944*alphaDrag[5]*F_0[14]+0.3912303982179757*F_0[8]*alphaDrag[12]+0.5477225575051661*F_0[7]*alphaDrag[12]+0.6123724356957944*F_0[0]*alphaDrag[12]+0.3912303982179757*alphaDrag[8]*F_0[12]+0.5477225575051661*alphaDrag[7]*F_0[12]+0.6123724356957944*alphaDrag[0]*F_0[12]+0.4898979485566356*F_0[4]*alphaDrag[11]+0.4898979485566356*alphaDrag[4]*F_0[11]+0.5477225575051661*F_0[6]*alphaDrag[10]+0.5477225575051661*alphaDrag[6]*F_0[10]+0.6123724356957944*F_0[1]*alphaDrag[8]+0.6123724356957944*alphaDrag[1]*F_0[8]+0.5477225575051661*F_0[2]*alphaDrag[4]+0.5477225575051661*alphaDrag[2]*F_0[4]; 
  out_F_0[19] += 1.095445115010332*alphaDrag[14]*F_0[19]+1.095445115010332*alphaDrag[13]*F_0[19]+1.224744871391589*alphaDrag[3]*F_0[19]+1.095445115010332*F_0[16]*alphaDrag[18]+1.095445115010332*F_0[11]*alphaDrag[18]+1.224744871391589*F_0[2]*alphaDrag[18]+1.095445115010332*alphaDrag[11]*F_0[18]+1.224744871391589*alphaDrag[2]*F_0[18]+1.095445115010332*F_0[15]*alphaDrag[17]+1.095445115010332*F_0[12]*alphaDrag[17]+1.224744871391589*F_0[1]*alphaDrag[17]+1.095445115010332*alphaDrag[12]*F_0[17]+1.224744871391589*alphaDrag[1]*F_0[17]+1.224744871391589*alphaDrag[5]*F_0[16]+1.224744871391589*alphaDrag[6]*F_0[15]+1.224744871391589*F_0[4]*alphaDrag[14]+1.224744871391589*alphaDrag[4]*F_0[14]+1.224744871391589*F_0[4]*alphaDrag[13]+1.224744871391589*alphaDrag[4]*F_0[13]+1.224744871391589*F_0[6]*alphaDrag[12]+1.224744871391589*alphaDrag[6]*F_0[12]+1.224744871391589*F_0[5]*alphaDrag[11]+1.224744871391589*alphaDrag[5]*F_0[11]+1.224744871391589*F_0[9]*alphaDrag[10]+1.224744871391589*F_0[8]*alphaDrag[10]+1.224744871391589*F_0[7]*alphaDrag[10]+1.369306393762915*F_0[0]*alphaDrag[10]+1.224744871391589*alphaDrag[8]*F_0[10]+1.224744871391589*alphaDrag[7]*F_0[10]+1.369306393762915*alphaDrag[0]*F_0[10]+1.369306393762915*F_0[1]*alphaDrag[6]+1.369306393762915*alphaDrag[1]*F_0[6]+1.369306393762915*F_0[2]*alphaDrag[5]+1.369306393762915*alphaDrag[2]*F_0[5]+1.369306393762915*F_0[3]*alphaDrag[4]+1.369306393762915*alphaDrag[3]*F_0[4]; 
  out_G_1[3] += 0.6123724356957944*G_1[18]*alphaDrag[18]+0.6123724356957944*G_1[17]*alphaDrag[17]+0.6123724356957944*G_1[14]*alphaDrag[14]+0.6123724356957944*G_1[13]*alphaDrag[13]+0.6123724356957944*G_1[12]*alphaDrag[12]+0.6123724356957944*G_1[11]*alphaDrag[11]+0.6123724356957944*G_1[10]*alphaDrag[10]+0.6123724356957944*G_1[8]*alphaDrag[8]+0.6123724356957944*G_1[7]*alphaDrag[7]+0.6123724356957944*G_1[6]*alphaDrag[6]+0.6123724356957944*G_1[5]*alphaDrag[5]+0.6123724356957944*G_1[4]*alphaDrag[4]+0.6123724356957944*G_1[3]*alphaDrag[3]+0.6123724356957944*G_1[2]*alphaDrag[2]+0.6123724356957944*G_1[1]*alphaDrag[1]+0.6123724356957944*G_1[0]*alphaDrag[0]; 
  out_G_1[5] += 0.6123724356957944*G_1[14]*alphaDrag[18]+0.6123724356957944*alphaDrag[14]*G_1[18]+0.5477225575051661*G_1[10]*alphaDrag[17]+0.5477225575051661*alphaDrag[10]*G_1[17]+0.5477225575051661*G_1[5]*alphaDrag[13]+0.5477225575051661*alphaDrag[5]*G_1[13]+0.6123724356957944*G_1[8]*alphaDrag[12]+0.6123724356957944*alphaDrag[8]*G_1[12]+0.5477225575051661*G_1[4]*alphaDrag[11]+0.5477225575051661*alphaDrag[4]*G_1[11]+0.6123724356957944*G_1[6]*alphaDrag[10]+0.6123724356957944*alphaDrag[6]*G_1[10]+0.5477225575051661*G_1[1]*alphaDrag[7]+0.5477225575051661*alphaDrag[1]*G_1[7]+0.6123724356957944*G_1[3]*alphaDrag[5]+0.6123724356957944*alphaDrag[3]*G_1[5]+0.6123724356957944*G_1[2]*alphaDrag[4]+0.6123724356957944*alphaDrag[2]*G_1[4]+0.6123724356957944*G_1[0]*alphaDrag[1]+0.6123724356957944*alphaDrag[0]*G_1[1]; 
  out_G_1[6] += 0.5477225575051661*G_1[10]*alphaDrag[18]+0.5477225575051661*alphaDrag[10]*G_1[18]+0.6123724356957944*G_1[13]*alphaDrag[17]+0.6123724356957944*alphaDrag[13]*G_1[17]+0.5477225575051661*G_1[6]*alphaDrag[14]+0.5477225575051661*alphaDrag[6]*G_1[14]+0.5477225575051661*G_1[4]*alphaDrag[12]+0.5477225575051661*alphaDrag[4]*G_1[12]+0.6123724356957944*G_1[7]*alphaDrag[11]+0.6123724356957944*alphaDrag[7]*G_1[11]+0.6123724356957944*G_1[5]*alphaDrag[10]+0.6123724356957944*alphaDrag[5]*G_1[10]+0.5477225575051661*G_1[2]*alphaDrag[8]+0.5477225575051661*alphaDrag[2]*G_1[8]+0.6123724356957944*G_1[3]*alphaDrag[6]+0.6123724356957944*alphaDrag[3]*G_1[6]+0.6123724356957944*G_1[1]*alphaDrag[4]+0.6123724356957944*alphaDrag[1]*G_1[4]+0.6123724356957944*G_1[0]*alphaDrag[2]+0.6123724356957944*alphaDrag[0]*G_1[2]; 
  out_G_1[9] += 1.224744871391589*alphaDrag[10]*G_1[19]+1.369306393762915*G_1[12]*alphaDrag[18]+1.369306393762915*alphaDrag[12]*G_1[18]+1.369306393762915*G_1[11]*alphaDrag[17]+1.369306393762915*alphaDrag[11]*G_1[17]+1.224744871391589*alphaDrag[6]*G_1[16]+1.224744871391589*alphaDrag[5]*G_1[15]+1.369306393762915*G_1[8]*alphaDrag[14]+1.369306393762915*alphaDrag[8]*G_1[14]+1.369306393762915*G_1[7]*alphaDrag[13]+1.369306393762915*alphaDrag[7]*G_1[13]+1.369306393762915*G_1[4]*alphaDrag[10]+1.369306393762915*alphaDrag[4]*G_1[10]+1.224744871391589*alphaDrag[3]*G_1[9]+1.369306393762915*G_1[2]*alphaDrag[6]+1.369306393762915*alphaDrag[2]*G_1[6]+1.369306393762915*G_1[1]*alphaDrag[5]+1.369306393762915*alphaDrag[1]*G_1[5]+1.369306393762915*G_1[0]*alphaDrag[3]+1.369306393762915*alphaDrag[0]*G_1[3]; 
  out_G_1[10] += 0.4898979485566357*G_1[17]*alphaDrag[18]+0.5477225575051661*G_1[6]*alphaDrag[18]+0.4898979485566357*alphaDrag[17]*G_1[18]+0.5477225575051661*alphaDrag[6]*G_1[18]+0.5477225575051661*G_1[5]*alphaDrag[17]+0.5477225575051661*alphaDrag[5]*G_1[17]+0.5477225575051661*G_1[10]*alphaDrag[14]+0.5477225575051661*alphaDrag[10]*G_1[14]+0.5477225575051661*G_1[10]*alphaDrag[13]+0.5477225575051661*alphaDrag[10]*G_1[13]+0.4898979485566357*G_1[11]*alphaDrag[12]+0.5477225575051661*G_1[2]*alphaDrag[12]+0.4898979485566357*alphaDrag[11]*G_1[12]+0.5477225575051661*alphaDrag[2]*G_1[12]+0.5477225575051661*G_1[1]*alphaDrag[11]+0.5477225575051661*alphaDrag[1]*G_1[11]+0.6123724356957944*G_1[3]*alphaDrag[10]+0.6123724356957944*alphaDrag[3]*G_1[10]+0.5477225575051661*G_1[4]*alphaDrag[8]+0.5477225575051661*alphaDrag[4]*G_1[8]+0.5477225575051661*G_1[4]*alphaDrag[7]+0.5477225575051661*alphaDrag[4]*G_1[7]+0.6123724356957944*G_1[5]*alphaDrag[6]+0.6123724356957944*alphaDrag[5]*G_1[6]+0.6123724356957944*G_1[0]*alphaDrag[4]+0.6123724356957944*alphaDrag[0]*G_1[4]+0.6123724356957944*G_1[1]*alphaDrag[2]+0.6123724356957944*alphaDrag[1]*G_1[2]; 
  out_G_1[13] += 0.5477225575051661*G_1[18]*alphaDrag[18]+0.3912303982179757*G_1[17]*alphaDrag[17]+0.6123724356957944*G_1[6]*alphaDrag[17]+0.6123724356957944*alphaDrag[6]*G_1[17]+0.3912303982179757*G_1[13]*alphaDrag[13]+0.6123724356957944*G_1[3]*alphaDrag[13]+0.6123724356957944*alphaDrag[3]*G_1[13]+0.5477225575051661*G_1[12]*alphaDrag[12]+0.3912303982179757*G_1[11]*alphaDrag[11]+0.6123724356957944*G_1[2]*alphaDrag[11]+0.6123724356957944*alphaDrag[2]*G_1[11]+0.5477225575051661*G_1[10]*alphaDrag[10]+0.3912303982179757*G_1[7]*alphaDrag[7]+0.6123724356957944*G_1[0]*alphaDrag[7]+0.6123724356957944*alphaDrag[0]*G_1[7]+0.5477225575051661*G_1[5]*alphaDrag[5]+0.5477225575051661*G_1[4]*alphaDrag[4]+0.5477225575051661*G_1[1]*alphaDrag[1]; 
  out_G_1[14] += 0.3912303982179757*G_1[18]*alphaDrag[18]+0.6123724356957944*G_1[5]*alphaDrag[18]+0.6123724356957944*alphaDrag[5]*G_1[18]+0.5477225575051661*G_1[17]*alphaDrag[17]+0.3912303982179757*G_1[14]*alphaDrag[14]+0.6123724356957944*G_1[3]*alphaDrag[14]+0.6123724356957944*alphaDrag[3]*G_1[14]+0.3912303982179757*G_1[12]*alphaDrag[12]+0.6123724356957944*G_1[1]*alphaDrag[12]+0.6123724356957944*alphaDrag[1]*G_1[12]+0.5477225575051661*G_1[11]*alphaDrag[11]+0.5477225575051661*G_1[10]*alphaDrag[10]+0.3912303982179757*G_1[8]*alphaDrag[8]+0.6123724356957944*G_1[0]*alphaDrag[8]+0.6123724356957944*alphaDrag[0]*G_1[8]+0.5477225575051661*G_1[6]*alphaDrag[6]+0.5477225575051661*G_1[4]*alphaDrag[4]+0.5477225575051661*G_1[2]*alphaDrag[2]; 
  out_G_1[15] += 1.095445115010332*alphaDrag[17]*G_1[19]+1.224744871391589*alphaDrag[6]*G_1[19]+1.369306393762915*G_1[8]*alphaDrag[18]+1.369306393762915*alphaDrag[8]*G_1[18]+1.224744871391589*G_1[4]*alphaDrag[17]+1.224744871391589*alphaDrag[4]*G_1[17]+1.224744871391589*alphaDrag[10]*G_1[16]+1.095445115010332*alphaDrag[13]*G_1[15]+1.224744871391589*alphaDrag[3]*G_1[15]+1.369306393762915*G_1[12]*alphaDrag[14]+1.369306393762915*alphaDrag[12]*G_1[14]+1.224744871391589*G_1[1]*alphaDrag[13]+1.224744871391589*alphaDrag[1]*G_1[13]+1.224744871391589*G_1[10]*alphaDrag[11]+1.224744871391589*alphaDrag[10]*G_1[11]+1.369306393762915*G_1[2]*alphaDrag[10]+1.369306393762915*alphaDrag[2]*G_1[10]+1.224744871391589*alphaDrag[5]*G_1[9]+1.224744871391589*G_1[5]*alphaDrag[7]+1.224744871391589*alphaDrag[5]*G_1[7]+1.369306393762915*G_1[4]*alphaDrag[6]+1.369306393762915*alphaDrag[4]*G_1[6]+1.369306393762915*G_1[0]*alphaDrag[5]+1.369306393762915*alphaDrag[0]*G_1[5]+1.369306393762915*G_1[1]*alphaDrag[3]+1.369306393762915*alphaDrag[1]*G_1[3]; 
  out_G_1[16] += 1.095445115010332*alphaDrag[18]*G_1[19]+1.224744871391589*alphaDrag[5]*G_1[19]+1.224744871391589*G_1[4]*alphaDrag[18]+1.224744871391589*alphaDrag[4]*G_1[18]+1.369306393762915*G_1[7]*alphaDrag[17]+1.369306393762915*alphaDrag[7]*G_1[17]+1.095445115010332*alphaDrag[14]*G_1[16]+1.224744871391589*alphaDrag[3]*G_1[16]+1.224744871391589*alphaDrag[10]*G_1[15]+1.224744871391589*G_1[2]*alphaDrag[14]+1.224744871391589*alphaDrag[2]*G_1[14]+1.369306393762915*G_1[11]*alphaDrag[13]+1.369306393762915*alphaDrag[11]*G_1[13]+1.224744871391589*G_1[10]*alphaDrag[12]+1.224744871391589*alphaDrag[10]*G_1[12]+1.369306393762915*G_1[1]*alphaDrag[10]+1.369306393762915*alphaDrag[1]*G_1[10]+1.224744871391589*alphaDrag[6]*G_1[9]+1.224744871391589*G_1[6]*alphaDrag[8]+1.224744871391589*alphaDrag[6]*G_1[8]+1.369306393762915*G_1[0]*alphaDrag[6]+1.369306393762915*alphaDrag[0]*G_1[6]+1.369306393762915*G_1[4]*alphaDrag[5]+1.369306393762915*alphaDrag[4]*G_1[5]+1.369306393762915*G_1[2]*alphaDrag[3]+1.369306393762915*alphaDrag[2]*G_1[3]; 
  out_G_1[17] += 0.4898979485566357*G_1[10]*alphaDrag[18]+0.4898979485566357*alphaDrag[10]*G_1[18]+0.5477225575051661*G_1[14]*alphaDrag[17]+0.3912303982179757*G_1[13]*alphaDrag[17]+0.6123724356957944*G_1[3]*alphaDrag[17]+0.5477225575051661*alphaDrag[14]*G_1[17]+0.3912303982179757*alphaDrag[13]*G_1[17]+0.6123724356957944*alphaDrag[3]*G_1[17]+0.6123724356957944*G_1[6]*alphaDrag[13]+0.6123724356957944*alphaDrag[6]*G_1[13]+0.4898979485566356*G_1[4]*alphaDrag[12]+0.4898979485566356*alphaDrag[4]*G_1[12]+0.5477225575051661*G_1[8]*alphaDrag[11]+0.3912303982179757*G_1[7]*alphaDrag[11]+0.6123724356957944*G_1[0]*alphaDrag[11]+0.5477225575051661*alphaDrag[8]*G_1[11]+0.3912303982179757*alphaDrag[7]*G_1[11]+0.6123724356957944*alphaDrag[0]*G_1[11]+0.5477225575051661*G_1[5]*alphaDrag[10]+0.5477225575051661*alphaDrag[5]*G_1[10]+0.6123724356957944*G_1[2]*alphaDrag[7]+0.6123724356957944*alphaDrag[2]*G_1[7]+0.5477225575051661*G_1[1]*alphaDrag[4]+0.5477225575051661*alphaDrag[1]*G_1[4]; 
  out_G_1[18] += 0.3912303982179757*G_1[14]*alphaDrag[18]+0.5477225575051661*G_1[13]*alphaDrag[18]+0.6123724356957944*G_1[3]*alphaDrag[18]+0.3912303982179757*alphaDrag[14]*G_1[18]+0.5477225575051661*alphaDrag[13]*G_1[18]+0.6123724356957944*alphaDrag[3]*G_1[18]+0.4898979485566357*G_1[10]*alphaDrag[17]+0.4898979485566357*alphaDrag[10]*G_1[17]+0.6123724356957944*G_1[5]*alphaDrag[14]+0.6123724356957944*alphaDrag[5]*G_1[14]+0.3912303982179757*G_1[8]*alphaDrag[12]+0.5477225575051661*G_1[7]*alphaDrag[12]+0.6123724356957944*G_1[0]*alphaDrag[12]+0.3912303982179757*alphaDrag[8]*G_1[12]+0.5477225575051661*alphaDrag[7]*G_1[12]+0.6123724356957944*alphaDrag[0]*G_1[12]+0.4898979485566356*G_1[4]*alphaDrag[11]+0.4898979485566356*alphaDrag[4]*G_1[11]+0.5477225575051661*G_1[6]*alphaDrag[10]+0.5477225575051661*alphaDrag[6]*G_1[10]+0.6123724356957944*G_1[1]*alphaDrag[8]+0.6123724356957944*alphaDrag[1]*G_1[8]+0.5477225575051661*G_1[2]*alphaDrag[4]+0.5477225575051661*alphaDrag[2]*G_1[4]; 
  out_G_1[19] += 1.095445115010332*alphaDrag[14]*G_1[19]+1.095445115010332*alphaDrag[13]*G_1[19]+1.224744871391589*alphaDrag[3]*G_1[19]+1.095445115010332*G_1[16]*alphaDrag[18]+1.095445115010332*G_1[11]*alphaDrag[18]+1.224744871391589*G_1[2]*alphaDrag[18]+1.095445115010332*alphaDrag[11]*G_1[18]+1.224744871391589*alphaDrag[2]*G_1[18]+1.095445115010332*G_1[15]*alphaDrag[17]+1.095445115010332*G_1[12]*alphaDrag[17]+1.224744871391589*G_1[1]*alphaDrag[17]+1.095445115010332*alphaDrag[12]*G_1[17]+1.224744871391589*alphaDrag[1]*G_1[17]+1.224744871391589*alphaDrag[5]*G_1[16]+1.224744871391589*alphaDrag[6]*G_1[15]+1.224744871391589*G_1[4]*alphaDrag[14]+1.224744871391589*alphaDrag[4]*G_1[14]+1.224744871391589*G_1[4]*alphaDrag[13]+1.224744871391589*alphaDrag[4]*G_1[13]+1.224744871391589*G_1[6]*alphaDrag[12]+1.224744871391589*alphaDrag[6]*G_1[12]+1.224744871391589*G_1[5]*alphaDrag[11]+1.224744871391589*alphaDrag[5]*G_1[11]+1.224744871391589*G_1[9]*alphaDrag[10]+1.224744871391589*G_1[8]*alphaDrag[10]+1.224744871391589*G_1[7]*alphaDrag[10]+1.369306393762915*G_1[0]*alphaDrag[10]+1.224744871391589*alphaDrag[8]*G_1[10]+1.224744871391589*alphaDrag[7]*G_1[10]+1.369306393762915*alphaDrag[0]*G_1[10]+1.369306393762915*G_1[1]*alphaDrag[6]+1.369306393762915*alphaDrag[1]*G_1[6]+1.369306393762915*G_1[2]*alphaDrag[5]+1.369306393762915*alphaDrag[2]*G_1[5]+1.369306393762915*G_1[3]*alphaDrag[4]+1.369306393762915*alphaDrag[3]*G_1[4]; 

  return fabs(0.8838834764831842*alphaDrag[0]-0.9882117688026182*(alphaDrag[8]+alphaDrag[7])); 

} 
