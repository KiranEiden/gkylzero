#include <gkyl_lbo_vlasov_pkpm_kernels.h> 
GKYL_CU_DH double lbo_vlasov_pkpm_drag_vol_3x1v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:       Cell-center coordinates. 
  // dxv[NDIM]:     Cell spacing. 
  // nuSum:         Collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: Sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // f:             Input distribution functions [F_0, T_perp/m G_1 = T_perp/m (F_0 - F_1)].
  // out:           Incremented output distribution functions. 
  const double rdvpar = 2.0/dxv[3]; 
  const double dvpar = dxv[3], wvpar = w[3]; 
  const double *F_0 = &f[0]; 
  const double *G_1 = &f[24]; 
  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[24]; 

  const double *sumNuUPar = &nuPrimMomsSum[0]; 

  double alphaDrag[24]; 
  alphaDrag[0] = rdvpar*(1.414213562373095*sumNuUPar[0]-1.414213562373095*nuSum[0]*wvpar); 
  alphaDrag[1] = rdvpar*(1.414213562373095*sumNuUPar[1]-1.414213562373095*nuSum[1]*wvpar); 
  alphaDrag[2] = rdvpar*(1.414213562373095*sumNuUPar[2]-1.414213562373095*nuSum[2]*wvpar); 
  alphaDrag[3] = rdvpar*(1.414213562373095*sumNuUPar[3]-1.414213562373095*nuSum[3]*wvpar); 
  alphaDrag[4] = -0.408248290463863*nuSum[0]*dvpar*rdvpar; 
  alphaDrag[5] = rdvpar*(1.414213562373095*sumNuUPar[4]-1.414213562373095*nuSum[4]*wvpar); 
  alphaDrag[6] = rdvpar*(1.414213562373095*sumNuUPar[5]-1.414213562373095*nuSum[5]*wvpar); 
  alphaDrag[7] = rdvpar*(1.414213562373095*sumNuUPar[6]-1.414213562373095*nuSum[6]*wvpar); 
  alphaDrag[8] = -0.408248290463863*nuSum[1]*dvpar*rdvpar; 
  alphaDrag[9] = -0.408248290463863*nuSum[2]*dvpar*rdvpar; 
  alphaDrag[10] = -0.408248290463863*nuSum[3]*dvpar*rdvpar; 
  alphaDrag[11] = rdvpar*(1.414213562373095*sumNuUPar[7]-1.414213562373095*nuSum[7]*wvpar); 
  alphaDrag[12] = -0.408248290463863*nuSum[4]*dvpar*rdvpar; 
  alphaDrag[13] = -0.408248290463863*nuSum[5]*dvpar*rdvpar; 
  alphaDrag[14] = -0.408248290463863*nuSum[6]*dvpar*rdvpar; 
  alphaDrag[15] = -0.408248290463863*nuSum[7]*dvpar*rdvpar; 

  out_F_0[4] += 0.4330127018922193*F_0[15]*alphaDrag[15]+0.4330127018922193*F_0[14]*alphaDrag[14]+0.4330127018922193*F_0[13]*alphaDrag[13]+0.4330127018922193*F_0[12]*alphaDrag[12]+0.4330127018922193*F_0[11]*alphaDrag[11]+0.4330127018922193*F_0[10]*alphaDrag[10]+0.4330127018922193*F_0[9]*alphaDrag[9]+0.4330127018922193*F_0[8]*alphaDrag[8]+0.4330127018922193*F_0[7]*alphaDrag[7]+0.4330127018922193*F_0[6]*alphaDrag[6]+0.4330127018922193*F_0[5]*alphaDrag[5]+0.4330127018922193*F_0[4]*alphaDrag[4]+0.4330127018922193*F_0[3]*alphaDrag[3]+0.4330127018922193*F_0[2]*alphaDrag[2]+0.4330127018922193*F_0[1]*alphaDrag[1]+0.4330127018922193*F_0[0]*alphaDrag[0]; 
  out_F_0[8] += 0.4330127018922193*F_0[14]*alphaDrag[15]+0.4330127018922193*alphaDrag[14]*F_0[15]+0.4330127018922193*F_0[10]*alphaDrag[13]+0.4330127018922193*alphaDrag[10]*F_0[13]+0.4330127018922193*F_0[9]*alphaDrag[12]+0.4330127018922193*alphaDrag[9]*F_0[12]+0.4330127018922193*F_0[7]*alphaDrag[11]+0.4330127018922193*alphaDrag[7]*F_0[11]+0.4330127018922193*F_0[4]*alphaDrag[8]+0.4330127018922193*alphaDrag[4]*F_0[8]+0.4330127018922193*F_0[3]*alphaDrag[6]+0.4330127018922193*alphaDrag[3]*F_0[6]+0.4330127018922193*F_0[2]*alphaDrag[5]+0.4330127018922193*alphaDrag[2]*F_0[5]+0.4330127018922193*F_0[0]*alphaDrag[1]+0.4330127018922193*alphaDrag[0]*F_0[1]; 
  out_F_0[9] += 0.4330127018922193*F_0[13]*alphaDrag[15]+0.4330127018922193*alphaDrag[13]*F_0[15]+0.4330127018922193*F_0[10]*alphaDrag[14]+0.4330127018922193*alphaDrag[10]*F_0[14]+0.4330127018922193*F_0[8]*alphaDrag[12]+0.4330127018922193*alphaDrag[8]*F_0[12]+0.4330127018922193*F_0[6]*alphaDrag[11]+0.4330127018922193*alphaDrag[6]*F_0[11]+0.4330127018922193*F_0[4]*alphaDrag[9]+0.4330127018922193*alphaDrag[4]*F_0[9]+0.4330127018922193*F_0[3]*alphaDrag[7]+0.4330127018922193*alphaDrag[3]*F_0[7]+0.4330127018922193*F_0[1]*alphaDrag[5]+0.4330127018922193*alphaDrag[1]*F_0[5]+0.4330127018922193*F_0[0]*alphaDrag[2]+0.4330127018922193*alphaDrag[0]*F_0[2]; 
  out_F_0[10] += 0.4330127018922193*F_0[12]*alphaDrag[15]+0.4330127018922193*alphaDrag[12]*F_0[15]+0.4330127018922193*F_0[9]*alphaDrag[14]+0.4330127018922193*alphaDrag[9]*F_0[14]+0.4330127018922193*F_0[8]*alphaDrag[13]+0.4330127018922193*alphaDrag[8]*F_0[13]+0.4330127018922193*F_0[5]*alphaDrag[11]+0.4330127018922193*alphaDrag[5]*F_0[11]+0.4330127018922193*F_0[4]*alphaDrag[10]+0.4330127018922193*alphaDrag[4]*F_0[10]+0.4330127018922193*F_0[2]*alphaDrag[7]+0.4330127018922193*alphaDrag[2]*F_0[7]+0.4330127018922193*F_0[1]*alphaDrag[6]+0.4330127018922193*alphaDrag[1]*F_0[6]+0.4330127018922193*F_0[0]*alphaDrag[3]+0.4330127018922193*alphaDrag[0]*F_0[3]; 
  out_F_0[12] += 0.4330127018922193*F_0[10]*alphaDrag[15]+0.4330127018922193*alphaDrag[10]*F_0[15]+0.4330127018922193*F_0[13]*alphaDrag[14]+0.4330127018922193*alphaDrag[13]*F_0[14]+0.4330127018922193*F_0[4]*alphaDrag[12]+0.4330127018922193*alphaDrag[4]*F_0[12]+0.4330127018922193*F_0[3]*alphaDrag[11]+0.4330127018922193*alphaDrag[3]*F_0[11]+0.4330127018922193*F_0[8]*alphaDrag[9]+0.4330127018922193*alphaDrag[8]*F_0[9]+0.4330127018922193*F_0[6]*alphaDrag[7]+0.4330127018922193*alphaDrag[6]*F_0[7]+0.4330127018922193*F_0[0]*alphaDrag[5]+0.4330127018922193*alphaDrag[0]*F_0[5]+0.4330127018922193*F_0[1]*alphaDrag[2]+0.4330127018922193*alphaDrag[1]*F_0[2]; 
  out_F_0[13] += 0.4330127018922193*F_0[9]*alphaDrag[15]+0.4330127018922193*alphaDrag[9]*F_0[15]+0.4330127018922193*F_0[12]*alphaDrag[14]+0.4330127018922193*alphaDrag[12]*F_0[14]+0.4330127018922193*F_0[4]*alphaDrag[13]+0.4330127018922193*alphaDrag[4]*F_0[13]+0.4330127018922193*F_0[2]*alphaDrag[11]+0.4330127018922193*alphaDrag[2]*F_0[11]+0.4330127018922193*F_0[8]*alphaDrag[10]+0.4330127018922193*alphaDrag[8]*F_0[10]+0.4330127018922193*F_0[5]*alphaDrag[7]+0.4330127018922193*alphaDrag[5]*F_0[7]+0.4330127018922193*F_0[0]*alphaDrag[6]+0.4330127018922193*alphaDrag[0]*F_0[6]+0.4330127018922193*F_0[1]*alphaDrag[3]+0.4330127018922193*alphaDrag[1]*F_0[3]; 
  out_F_0[14] += 0.4330127018922193*F_0[8]*alphaDrag[15]+0.4330127018922193*alphaDrag[8]*F_0[15]+0.4330127018922193*F_0[4]*alphaDrag[14]+0.4330127018922193*alphaDrag[4]*F_0[14]+0.4330127018922193*F_0[12]*alphaDrag[13]+0.4330127018922193*alphaDrag[12]*F_0[13]+0.4330127018922193*F_0[1]*alphaDrag[11]+0.4330127018922193*alphaDrag[1]*F_0[11]+0.4330127018922193*F_0[9]*alphaDrag[10]+0.4330127018922193*alphaDrag[9]*F_0[10]+0.4330127018922193*F_0[0]*alphaDrag[7]+0.4330127018922193*alphaDrag[0]*F_0[7]+0.4330127018922193*F_0[5]*alphaDrag[6]+0.4330127018922193*alphaDrag[5]*F_0[6]+0.4330127018922193*F_0[2]*alphaDrag[3]+0.4330127018922193*alphaDrag[2]*F_0[3]; 
  out_F_0[15] += 0.4330127018922193*F_0[4]*alphaDrag[15]+0.4330127018922193*alphaDrag[4]*F_0[15]+0.4330127018922193*F_0[8]*alphaDrag[14]+0.4330127018922193*alphaDrag[8]*F_0[14]+0.4330127018922193*F_0[9]*alphaDrag[13]+0.4330127018922193*alphaDrag[9]*F_0[13]+0.4330127018922193*F_0[10]*alphaDrag[12]+0.4330127018922193*alphaDrag[10]*F_0[12]+0.4330127018922193*F_0[0]*alphaDrag[11]+0.4330127018922193*alphaDrag[0]*F_0[11]+0.4330127018922193*F_0[1]*alphaDrag[7]+0.4330127018922193*alphaDrag[1]*F_0[7]+0.4330127018922193*F_0[2]*alphaDrag[6]+0.4330127018922193*alphaDrag[2]*F_0[6]+0.4330127018922193*F_0[3]*alphaDrag[5]+0.4330127018922193*alphaDrag[3]*F_0[5]; 
  out_F_0[16] += 0.8660254037844387*alphaDrag[15]*F_0[23]+0.8660254037844386*alphaDrag[14]*F_0[22]+0.8660254037844386*alphaDrag[13]*F_0[21]+0.8660254037844386*alphaDrag[12]*F_0[20]+0.8660254037844387*alphaDrag[10]*F_0[19]+0.8660254037844387*alphaDrag[9]*F_0[18]+0.8660254037844387*alphaDrag[8]*F_0[17]+0.8660254037844386*alphaDrag[4]*F_0[16]+0.9682458365518543*F_0[11]*alphaDrag[15]+0.9682458365518543*alphaDrag[11]*F_0[15]+0.9682458365518543*F_0[7]*alphaDrag[14]+0.9682458365518543*alphaDrag[7]*F_0[14]+0.9682458365518543*F_0[6]*alphaDrag[13]+0.9682458365518543*alphaDrag[6]*F_0[13]+0.9682458365518543*F_0[5]*alphaDrag[12]+0.9682458365518543*alphaDrag[5]*F_0[12]+0.9682458365518543*F_0[3]*alphaDrag[10]+0.9682458365518543*alphaDrag[3]*F_0[10]+0.9682458365518543*F_0[2]*alphaDrag[9]+0.9682458365518543*alphaDrag[2]*F_0[9]+0.9682458365518543*F_0[1]*alphaDrag[8]+0.9682458365518543*alphaDrag[1]*F_0[8]+0.9682458365518543*F_0[0]*alphaDrag[4]+0.9682458365518543*alphaDrag[0]*F_0[4]; 
  out_F_0[17] += 0.8660254037844386*alphaDrag[14]*F_0[23]+0.8660254037844387*alphaDrag[15]*F_0[22]+0.8660254037844387*alphaDrag[10]*F_0[21]+0.8660254037844387*alphaDrag[9]*F_0[20]+0.8660254037844386*alphaDrag[13]*F_0[19]+0.8660254037844386*alphaDrag[12]*F_0[18]+0.8660254037844386*alphaDrag[4]*F_0[17]+0.8660254037844387*alphaDrag[8]*F_0[16]+0.9682458365518543*F_0[7]*alphaDrag[15]+0.9682458365518543*alphaDrag[7]*F_0[15]+0.9682458365518543*F_0[11]*alphaDrag[14]+0.9682458365518543*alphaDrag[11]*F_0[14]+0.9682458365518543*F_0[3]*alphaDrag[13]+0.9682458365518543*alphaDrag[3]*F_0[13]+0.9682458365518543*F_0[2]*alphaDrag[12]+0.9682458365518543*alphaDrag[2]*F_0[12]+0.9682458365518543*F_0[6]*alphaDrag[10]+0.9682458365518543*alphaDrag[6]*F_0[10]+0.9682458365518543*F_0[5]*alphaDrag[9]+0.9682458365518543*alphaDrag[5]*F_0[9]+0.9682458365518543*F_0[0]*alphaDrag[8]+0.9682458365518543*alphaDrag[0]*F_0[8]+0.9682458365518543*F_0[1]*alphaDrag[4]+0.9682458365518543*alphaDrag[1]*F_0[4]; 
  out_F_0[18] += 0.8660254037844386*alphaDrag[13]*F_0[23]+0.8660254037844387*alphaDrag[10]*F_0[22]+0.8660254037844387*alphaDrag[15]*F_0[21]+0.8660254037844387*alphaDrag[8]*F_0[20]+0.8660254037844386*alphaDrag[14]*F_0[19]+0.8660254037844386*alphaDrag[4]*F_0[18]+0.8660254037844386*alphaDrag[12]*F_0[17]+0.8660254037844387*alphaDrag[9]*F_0[16]+0.9682458365518543*F_0[6]*alphaDrag[15]+0.9682458365518543*alphaDrag[6]*F_0[15]+0.9682458365518543*F_0[3]*alphaDrag[14]+0.9682458365518543*alphaDrag[3]*F_0[14]+0.9682458365518543*F_0[11]*alphaDrag[13]+0.9682458365518543*alphaDrag[11]*F_0[13]+0.9682458365518543*F_0[1]*alphaDrag[12]+0.9682458365518543*alphaDrag[1]*F_0[12]+0.9682458365518543*F_0[7]*alphaDrag[10]+0.9682458365518543*alphaDrag[7]*F_0[10]+0.9682458365518543*F_0[0]*alphaDrag[9]+0.9682458365518543*alphaDrag[0]*F_0[9]+0.9682458365518543*F_0[5]*alphaDrag[8]+0.9682458365518543*alphaDrag[5]*F_0[8]+0.9682458365518543*F_0[2]*alphaDrag[4]+0.9682458365518543*alphaDrag[2]*F_0[4]; 
  out_F_0[19] += 0.8660254037844386*alphaDrag[12]*F_0[23]+0.8660254037844387*alphaDrag[9]*F_0[22]+0.8660254037844387*alphaDrag[8]*F_0[21]+0.8660254037844387*alphaDrag[15]*F_0[20]+0.8660254037844386*alphaDrag[4]*F_0[19]+0.8660254037844386*alphaDrag[14]*F_0[18]+0.8660254037844386*alphaDrag[13]*F_0[17]+0.8660254037844387*alphaDrag[10]*F_0[16]+0.9682458365518543*F_0[5]*alphaDrag[15]+0.9682458365518543*alphaDrag[5]*F_0[15]+0.9682458365518543*F_0[2]*alphaDrag[14]+0.9682458365518543*alphaDrag[2]*F_0[14]+0.9682458365518543*F_0[1]*alphaDrag[13]+0.9682458365518543*alphaDrag[1]*F_0[13]+0.9682458365518543*F_0[11]*alphaDrag[12]+0.9682458365518543*alphaDrag[11]*F_0[12]+0.9682458365518543*F_0[0]*alphaDrag[10]+0.9682458365518543*alphaDrag[0]*F_0[10]+0.9682458365518543*F_0[7]*alphaDrag[9]+0.9682458365518543*alphaDrag[7]*F_0[9]+0.9682458365518543*F_0[6]*alphaDrag[8]+0.9682458365518543*alphaDrag[6]*F_0[8]+0.9682458365518543*F_0[3]*alphaDrag[4]+0.9682458365518543*alphaDrag[3]*F_0[4]; 
  out_F_0[20] += 0.8660254037844387*alphaDrag[10]*F_0[23]+0.8660254037844386*alphaDrag[13]*F_0[22]+0.8660254037844386*alphaDrag[14]*F_0[21]+0.8660254037844386*alphaDrag[4]*F_0[20]+0.8660254037844387*alphaDrag[15]*F_0[19]+0.8660254037844387*alphaDrag[8]*F_0[18]+0.8660254037844387*alphaDrag[9]*F_0[17]+0.8660254037844386*alphaDrag[12]*F_0[16]+0.9682458365518543*F_0[3]*alphaDrag[15]+0.9682458365518543*alphaDrag[3]*F_0[15]+0.9682458365518543*F_0[6]*alphaDrag[14]+0.9682458365518543*alphaDrag[6]*F_0[14]+0.9682458365518543*F_0[7]*alphaDrag[13]+0.9682458365518543*alphaDrag[7]*F_0[13]+0.9682458365518543*F_0[0]*alphaDrag[12]+0.9682458365518543*alphaDrag[0]*F_0[12]+0.9682458365518543*F_0[10]*alphaDrag[11]+0.9682458365518543*alphaDrag[10]*F_0[11]+0.9682458365518543*F_0[1]*alphaDrag[9]+0.9682458365518543*alphaDrag[1]*F_0[9]+0.9682458365518543*F_0[2]*alphaDrag[8]+0.9682458365518543*alphaDrag[2]*F_0[8]+0.9682458365518543*F_0[4]*alphaDrag[5]+0.9682458365518543*alphaDrag[4]*F_0[5]; 
  out_F_0[21] += 0.8660254037844387*alphaDrag[9]*F_0[23]+0.8660254037844386*alphaDrag[12]*F_0[22]+0.8660254037844386*alphaDrag[4]*F_0[21]+0.8660254037844386*alphaDrag[14]*F_0[20]+0.8660254037844387*alphaDrag[8]*F_0[19]+0.8660254037844387*alphaDrag[15]*F_0[18]+0.8660254037844387*alphaDrag[10]*F_0[17]+0.8660254037844386*alphaDrag[13]*F_0[16]+0.9682458365518543*F_0[2]*alphaDrag[15]+0.9682458365518543*alphaDrag[2]*F_0[15]+0.9682458365518543*F_0[5]*alphaDrag[14]+0.9682458365518543*alphaDrag[5]*F_0[14]+0.9682458365518543*F_0[0]*alphaDrag[13]+0.9682458365518543*alphaDrag[0]*F_0[13]+0.9682458365518543*F_0[7]*alphaDrag[12]+0.9682458365518543*alphaDrag[7]*F_0[12]+0.9682458365518543*F_0[9]*alphaDrag[11]+0.9682458365518543*alphaDrag[9]*F_0[11]+0.9682458365518543*F_0[1]*alphaDrag[10]+0.9682458365518543*alphaDrag[1]*F_0[10]+0.9682458365518543*F_0[3]*alphaDrag[8]+0.9682458365518543*alphaDrag[3]*F_0[8]+0.9682458365518543*F_0[4]*alphaDrag[6]+0.9682458365518543*alphaDrag[4]*F_0[6]; 
  out_F_0[22] += 0.8660254037844387*alphaDrag[8]*F_0[23]+0.8660254037844386*alphaDrag[4]*F_0[22]+0.8660254037844386*alphaDrag[12]*F_0[21]+0.8660254037844386*alphaDrag[13]*F_0[20]+0.8660254037844387*alphaDrag[9]*F_0[19]+0.8660254037844387*alphaDrag[10]*F_0[18]+0.8660254037844387*alphaDrag[15]*F_0[17]+0.8660254037844386*alphaDrag[14]*F_0[16]+0.9682458365518543*F_0[1]*alphaDrag[15]+0.9682458365518543*alphaDrag[1]*F_0[15]+0.9682458365518543*F_0[0]*alphaDrag[14]+0.9682458365518543*alphaDrag[0]*F_0[14]+0.9682458365518543*F_0[5]*alphaDrag[13]+0.9682458365518543*alphaDrag[5]*F_0[13]+0.9682458365518543*F_0[6]*alphaDrag[12]+0.9682458365518543*alphaDrag[6]*F_0[12]+0.9682458365518543*F_0[8]*alphaDrag[11]+0.9682458365518543*alphaDrag[8]*F_0[11]+0.9682458365518543*F_0[2]*alphaDrag[10]+0.9682458365518543*alphaDrag[2]*F_0[10]+0.9682458365518543*F_0[3]*alphaDrag[9]+0.9682458365518543*alphaDrag[3]*F_0[9]+0.9682458365518543*F_0[4]*alphaDrag[7]+0.9682458365518543*alphaDrag[4]*F_0[7]; 
  out_F_0[23] += 0.8660254037844386*alphaDrag[4]*F_0[23]+0.8660254037844387*alphaDrag[8]*F_0[22]+0.8660254037844387*alphaDrag[9]*F_0[21]+0.8660254037844387*alphaDrag[10]*F_0[20]+0.8660254037844386*alphaDrag[12]*F_0[19]+0.8660254037844386*alphaDrag[13]*F_0[18]+0.8660254037844386*alphaDrag[14]*F_0[17]+0.8660254037844387*alphaDrag[15]*F_0[16]+0.9682458365518543*F_0[0]*alphaDrag[15]+0.9682458365518543*alphaDrag[0]*F_0[15]+0.9682458365518543*F_0[1]*alphaDrag[14]+0.9682458365518543*alphaDrag[1]*F_0[14]+0.9682458365518543*F_0[2]*alphaDrag[13]+0.9682458365518543*alphaDrag[2]*F_0[13]+0.9682458365518543*F_0[3]*alphaDrag[12]+0.9682458365518543*alphaDrag[3]*F_0[12]+0.9682458365518543*F_0[4]*alphaDrag[11]+0.9682458365518543*alphaDrag[4]*F_0[11]+0.9682458365518543*F_0[5]*alphaDrag[10]+0.9682458365518543*alphaDrag[5]*F_0[10]+0.9682458365518543*F_0[6]*alphaDrag[9]+0.9682458365518543*alphaDrag[6]*F_0[9]+0.9682458365518543*F_0[7]*alphaDrag[8]+0.9682458365518543*alphaDrag[7]*F_0[8]; 
  out_G_1[4] += 0.4330127018922193*G_1[15]*alphaDrag[15]+0.4330127018922193*G_1[14]*alphaDrag[14]+0.4330127018922193*G_1[13]*alphaDrag[13]+0.4330127018922193*G_1[12]*alphaDrag[12]+0.4330127018922193*G_1[11]*alphaDrag[11]+0.4330127018922193*G_1[10]*alphaDrag[10]+0.4330127018922193*G_1[9]*alphaDrag[9]+0.4330127018922193*G_1[8]*alphaDrag[8]+0.4330127018922193*G_1[7]*alphaDrag[7]+0.4330127018922193*G_1[6]*alphaDrag[6]+0.4330127018922193*G_1[5]*alphaDrag[5]+0.4330127018922193*G_1[4]*alphaDrag[4]+0.4330127018922193*G_1[3]*alphaDrag[3]+0.4330127018922193*G_1[2]*alphaDrag[2]+0.4330127018922193*G_1[1]*alphaDrag[1]+0.4330127018922193*G_1[0]*alphaDrag[0]; 
  out_G_1[8] += 0.4330127018922193*G_1[14]*alphaDrag[15]+0.4330127018922193*alphaDrag[14]*G_1[15]+0.4330127018922193*G_1[10]*alphaDrag[13]+0.4330127018922193*alphaDrag[10]*G_1[13]+0.4330127018922193*G_1[9]*alphaDrag[12]+0.4330127018922193*alphaDrag[9]*G_1[12]+0.4330127018922193*G_1[7]*alphaDrag[11]+0.4330127018922193*alphaDrag[7]*G_1[11]+0.4330127018922193*G_1[4]*alphaDrag[8]+0.4330127018922193*alphaDrag[4]*G_1[8]+0.4330127018922193*G_1[3]*alphaDrag[6]+0.4330127018922193*alphaDrag[3]*G_1[6]+0.4330127018922193*G_1[2]*alphaDrag[5]+0.4330127018922193*alphaDrag[2]*G_1[5]+0.4330127018922193*G_1[0]*alphaDrag[1]+0.4330127018922193*alphaDrag[0]*G_1[1]; 
  out_G_1[9] += 0.4330127018922193*G_1[13]*alphaDrag[15]+0.4330127018922193*alphaDrag[13]*G_1[15]+0.4330127018922193*G_1[10]*alphaDrag[14]+0.4330127018922193*alphaDrag[10]*G_1[14]+0.4330127018922193*G_1[8]*alphaDrag[12]+0.4330127018922193*alphaDrag[8]*G_1[12]+0.4330127018922193*G_1[6]*alphaDrag[11]+0.4330127018922193*alphaDrag[6]*G_1[11]+0.4330127018922193*G_1[4]*alphaDrag[9]+0.4330127018922193*alphaDrag[4]*G_1[9]+0.4330127018922193*G_1[3]*alphaDrag[7]+0.4330127018922193*alphaDrag[3]*G_1[7]+0.4330127018922193*G_1[1]*alphaDrag[5]+0.4330127018922193*alphaDrag[1]*G_1[5]+0.4330127018922193*G_1[0]*alphaDrag[2]+0.4330127018922193*alphaDrag[0]*G_1[2]; 
  out_G_1[10] += 0.4330127018922193*G_1[12]*alphaDrag[15]+0.4330127018922193*alphaDrag[12]*G_1[15]+0.4330127018922193*G_1[9]*alphaDrag[14]+0.4330127018922193*alphaDrag[9]*G_1[14]+0.4330127018922193*G_1[8]*alphaDrag[13]+0.4330127018922193*alphaDrag[8]*G_1[13]+0.4330127018922193*G_1[5]*alphaDrag[11]+0.4330127018922193*alphaDrag[5]*G_1[11]+0.4330127018922193*G_1[4]*alphaDrag[10]+0.4330127018922193*alphaDrag[4]*G_1[10]+0.4330127018922193*G_1[2]*alphaDrag[7]+0.4330127018922193*alphaDrag[2]*G_1[7]+0.4330127018922193*G_1[1]*alphaDrag[6]+0.4330127018922193*alphaDrag[1]*G_1[6]+0.4330127018922193*G_1[0]*alphaDrag[3]+0.4330127018922193*alphaDrag[0]*G_1[3]; 
  out_G_1[12] += 0.4330127018922193*G_1[10]*alphaDrag[15]+0.4330127018922193*alphaDrag[10]*G_1[15]+0.4330127018922193*G_1[13]*alphaDrag[14]+0.4330127018922193*alphaDrag[13]*G_1[14]+0.4330127018922193*G_1[4]*alphaDrag[12]+0.4330127018922193*alphaDrag[4]*G_1[12]+0.4330127018922193*G_1[3]*alphaDrag[11]+0.4330127018922193*alphaDrag[3]*G_1[11]+0.4330127018922193*G_1[8]*alphaDrag[9]+0.4330127018922193*alphaDrag[8]*G_1[9]+0.4330127018922193*G_1[6]*alphaDrag[7]+0.4330127018922193*alphaDrag[6]*G_1[7]+0.4330127018922193*G_1[0]*alphaDrag[5]+0.4330127018922193*alphaDrag[0]*G_1[5]+0.4330127018922193*G_1[1]*alphaDrag[2]+0.4330127018922193*alphaDrag[1]*G_1[2]; 
  out_G_1[13] += 0.4330127018922193*G_1[9]*alphaDrag[15]+0.4330127018922193*alphaDrag[9]*G_1[15]+0.4330127018922193*G_1[12]*alphaDrag[14]+0.4330127018922193*alphaDrag[12]*G_1[14]+0.4330127018922193*G_1[4]*alphaDrag[13]+0.4330127018922193*alphaDrag[4]*G_1[13]+0.4330127018922193*G_1[2]*alphaDrag[11]+0.4330127018922193*alphaDrag[2]*G_1[11]+0.4330127018922193*G_1[8]*alphaDrag[10]+0.4330127018922193*alphaDrag[8]*G_1[10]+0.4330127018922193*G_1[5]*alphaDrag[7]+0.4330127018922193*alphaDrag[5]*G_1[7]+0.4330127018922193*G_1[0]*alphaDrag[6]+0.4330127018922193*alphaDrag[0]*G_1[6]+0.4330127018922193*G_1[1]*alphaDrag[3]+0.4330127018922193*alphaDrag[1]*G_1[3]; 
  out_G_1[14] += 0.4330127018922193*G_1[8]*alphaDrag[15]+0.4330127018922193*alphaDrag[8]*G_1[15]+0.4330127018922193*G_1[4]*alphaDrag[14]+0.4330127018922193*alphaDrag[4]*G_1[14]+0.4330127018922193*G_1[12]*alphaDrag[13]+0.4330127018922193*alphaDrag[12]*G_1[13]+0.4330127018922193*G_1[1]*alphaDrag[11]+0.4330127018922193*alphaDrag[1]*G_1[11]+0.4330127018922193*G_1[9]*alphaDrag[10]+0.4330127018922193*alphaDrag[9]*G_1[10]+0.4330127018922193*G_1[0]*alphaDrag[7]+0.4330127018922193*alphaDrag[0]*G_1[7]+0.4330127018922193*G_1[5]*alphaDrag[6]+0.4330127018922193*alphaDrag[5]*G_1[6]+0.4330127018922193*G_1[2]*alphaDrag[3]+0.4330127018922193*alphaDrag[2]*G_1[3]; 
  out_G_1[15] += 0.4330127018922193*G_1[4]*alphaDrag[15]+0.4330127018922193*alphaDrag[4]*G_1[15]+0.4330127018922193*G_1[8]*alphaDrag[14]+0.4330127018922193*alphaDrag[8]*G_1[14]+0.4330127018922193*G_1[9]*alphaDrag[13]+0.4330127018922193*alphaDrag[9]*G_1[13]+0.4330127018922193*G_1[10]*alphaDrag[12]+0.4330127018922193*alphaDrag[10]*G_1[12]+0.4330127018922193*G_1[0]*alphaDrag[11]+0.4330127018922193*alphaDrag[0]*G_1[11]+0.4330127018922193*G_1[1]*alphaDrag[7]+0.4330127018922193*alphaDrag[1]*G_1[7]+0.4330127018922193*G_1[2]*alphaDrag[6]+0.4330127018922193*alphaDrag[2]*G_1[6]+0.4330127018922193*G_1[3]*alphaDrag[5]+0.4330127018922193*alphaDrag[3]*G_1[5]; 
  out_G_1[16] += 0.8660254037844387*alphaDrag[15]*G_1[23]+0.8660254037844386*alphaDrag[14]*G_1[22]+0.8660254037844386*alphaDrag[13]*G_1[21]+0.8660254037844386*alphaDrag[12]*G_1[20]+0.8660254037844387*alphaDrag[10]*G_1[19]+0.8660254037844387*alphaDrag[9]*G_1[18]+0.8660254037844387*alphaDrag[8]*G_1[17]+0.8660254037844386*alphaDrag[4]*G_1[16]+0.9682458365518543*G_1[11]*alphaDrag[15]+0.9682458365518543*alphaDrag[11]*G_1[15]+0.9682458365518543*G_1[7]*alphaDrag[14]+0.9682458365518543*alphaDrag[7]*G_1[14]+0.9682458365518543*G_1[6]*alphaDrag[13]+0.9682458365518543*alphaDrag[6]*G_1[13]+0.9682458365518543*G_1[5]*alphaDrag[12]+0.9682458365518543*alphaDrag[5]*G_1[12]+0.9682458365518543*G_1[3]*alphaDrag[10]+0.9682458365518543*alphaDrag[3]*G_1[10]+0.9682458365518543*G_1[2]*alphaDrag[9]+0.9682458365518543*alphaDrag[2]*G_1[9]+0.9682458365518543*G_1[1]*alphaDrag[8]+0.9682458365518543*alphaDrag[1]*G_1[8]+0.9682458365518543*G_1[0]*alphaDrag[4]+0.9682458365518543*alphaDrag[0]*G_1[4]; 
  out_G_1[17] += 0.8660254037844386*alphaDrag[14]*G_1[23]+0.8660254037844387*alphaDrag[15]*G_1[22]+0.8660254037844387*alphaDrag[10]*G_1[21]+0.8660254037844387*alphaDrag[9]*G_1[20]+0.8660254037844386*alphaDrag[13]*G_1[19]+0.8660254037844386*alphaDrag[12]*G_1[18]+0.8660254037844386*alphaDrag[4]*G_1[17]+0.8660254037844387*alphaDrag[8]*G_1[16]+0.9682458365518543*G_1[7]*alphaDrag[15]+0.9682458365518543*alphaDrag[7]*G_1[15]+0.9682458365518543*G_1[11]*alphaDrag[14]+0.9682458365518543*alphaDrag[11]*G_1[14]+0.9682458365518543*G_1[3]*alphaDrag[13]+0.9682458365518543*alphaDrag[3]*G_1[13]+0.9682458365518543*G_1[2]*alphaDrag[12]+0.9682458365518543*alphaDrag[2]*G_1[12]+0.9682458365518543*G_1[6]*alphaDrag[10]+0.9682458365518543*alphaDrag[6]*G_1[10]+0.9682458365518543*G_1[5]*alphaDrag[9]+0.9682458365518543*alphaDrag[5]*G_1[9]+0.9682458365518543*G_1[0]*alphaDrag[8]+0.9682458365518543*alphaDrag[0]*G_1[8]+0.9682458365518543*G_1[1]*alphaDrag[4]+0.9682458365518543*alphaDrag[1]*G_1[4]; 
  out_G_1[18] += 0.8660254037844386*alphaDrag[13]*G_1[23]+0.8660254037844387*alphaDrag[10]*G_1[22]+0.8660254037844387*alphaDrag[15]*G_1[21]+0.8660254037844387*alphaDrag[8]*G_1[20]+0.8660254037844386*alphaDrag[14]*G_1[19]+0.8660254037844386*alphaDrag[4]*G_1[18]+0.8660254037844386*alphaDrag[12]*G_1[17]+0.8660254037844387*alphaDrag[9]*G_1[16]+0.9682458365518543*G_1[6]*alphaDrag[15]+0.9682458365518543*alphaDrag[6]*G_1[15]+0.9682458365518543*G_1[3]*alphaDrag[14]+0.9682458365518543*alphaDrag[3]*G_1[14]+0.9682458365518543*G_1[11]*alphaDrag[13]+0.9682458365518543*alphaDrag[11]*G_1[13]+0.9682458365518543*G_1[1]*alphaDrag[12]+0.9682458365518543*alphaDrag[1]*G_1[12]+0.9682458365518543*G_1[7]*alphaDrag[10]+0.9682458365518543*alphaDrag[7]*G_1[10]+0.9682458365518543*G_1[0]*alphaDrag[9]+0.9682458365518543*alphaDrag[0]*G_1[9]+0.9682458365518543*G_1[5]*alphaDrag[8]+0.9682458365518543*alphaDrag[5]*G_1[8]+0.9682458365518543*G_1[2]*alphaDrag[4]+0.9682458365518543*alphaDrag[2]*G_1[4]; 
  out_G_1[19] += 0.8660254037844386*alphaDrag[12]*G_1[23]+0.8660254037844387*alphaDrag[9]*G_1[22]+0.8660254037844387*alphaDrag[8]*G_1[21]+0.8660254037844387*alphaDrag[15]*G_1[20]+0.8660254037844386*alphaDrag[4]*G_1[19]+0.8660254037844386*alphaDrag[14]*G_1[18]+0.8660254037844386*alphaDrag[13]*G_1[17]+0.8660254037844387*alphaDrag[10]*G_1[16]+0.9682458365518543*G_1[5]*alphaDrag[15]+0.9682458365518543*alphaDrag[5]*G_1[15]+0.9682458365518543*G_1[2]*alphaDrag[14]+0.9682458365518543*alphaDrag[2]*G_1[14]+0.9682458365518543*G_1[1]*alphaDrag[13]+0.9682458365518543*alphaDrag[1]*G_1[13]+0.9682458365518543*G_1[11]*alphaDrag[12]+0.9682458365518543*alphaDrag[11]*G_1[12]+0.9682458365518543*G_1[0]*alphaDrag[10]+0.9682458365518543*alphaDrag[0]*G_1[10]+0.9682458365518543*G_1[7]*alphaDrag[9]+0.9682458365518543*alphaDrag[7]*G_1[9]+0.9682458365518543*G_1[6]*alphaDrag[8]+0.9682458365518543*alphaDrag[6]*G_1[8]+0.9682458365518543*G_1[3]*alphaDrag[4]+0.9682458365518543*alphaDrag[3]*G_1[4]; 
  out_G_1[20] += 0.8660254037844387*alphaDrag[10]*G_1[23]+0.8660254037844386*alphaDrag[13]*G_1[22]+0.8660254037844386*alphaDrag[14]*G_1[21]+0.8660254037844386*alphaDrag[4]*G_1[20]+0.8660254037844387*alphaDrag[15]*G_1[19]+0.8660254037844387*alphaDrag[8]*G_1[18]+0.8660254037844387*alphaDrag[9]*G_1[17]+0.8660254037844386*alphaDrag[12]*G_1[16]+0.9682458365518543*G_1[3]*alphaDrag[15]+0.9682458365518543*alphaDrag[3]*G_1[15]+0.9682458365518543*G_1[6]*alphaDrag[14]+0.9682458365518543*alphaDrag[6]*G_1[14]+0.9682458365518543*G_1[7]*alphaDrag[13]+0.9682458365518543*alphaDrag[7]*G_1[13]+0.9682458365518543*G_1[0]*alphaDrag[12]+0.9682458365518543*alphaDrag[0]*G_1[12]+0.9682458365518543*G_1[10]*alphaDrag[11]+0.9682458365518543*alphaDrag[10]*G_1[11]+0.9682458365518543*G_1[1]*alphaDrag[9]+0.9682458365518543*alphaDrag[1]*G_1[9]+0.9682458365518543*G_1[2]*alphaDrag[8]+0.9682458365518543*alphaDrag[2]*G_1[8]+0.9682458365518543*G_1[4]*alphaDrag[5]+0.9682458365518543*alphaDrag[4]*G_1[5]; 
  out_G_1[21] += 0.8660254037844387*alphaDrag[9]*G_1[23]+0.8660254037844386*alphaDrag[12]*G_1[22]+0.8660254037844386*alphaDrag[4]*G_1[21]+0.8660254037844386*alphaDrag[14]*G_1[20]+0.8660254037844387*alphaDrag[8]*G_1[19]+0.8660254037844387*alphaDrag[15]*G_1[18]+0.8660254037844387*alphaDrag[10]*G_1[17]+0.8660254037844386*alphaDrag[13]*G_1[16]+0.9682458365518543*G_1[2]*alphaDrag[15]+0.9682458365518543*alphaDrag[2]*G_1[15]+0.9682458365518543*G_1[5]*alphaDrag[14]+0.9682458365518543*alphaDrag[5]*G_1[14]+0.9682458365518543*G_1[0]*alphaDrag[13]+0.9682458365518543*alphaDrag[0]*G_1[13]+0.9682458365518543*G_1[7]*alphaDrag[12]+0.9682458365518543*alphaDrag[7]*G_1[12]+0.9682458365518543*G_1[9]*alphaDrag[11]+0.9682458365518543*alphaDrag[9]*G_1[11]+0.9682458365518543*G_1[1]*alphaDrag[10]+0.9682458365518543*alphaDrag[1]*G_1[10]+0.9682458365518543*G_1[3]*alphaDrag[8]+0.9682458365518543*alphaDrag[3]*G_1[8]+0.9682458365518543*G_1[4]*alphaDrag[6]+0.9682458365518543*alphaDrag[4]*G_1[6]; 
  out_G_1[22] += 0.8660254037844387*alphaDrag[8]*G_1[23]+0.8660254037844386*alphaDrag[4]*G_1[22]+0.8660254037844386*alphaDrag[12]*G_1[21]+0.8660254037844386*alphaDrag[13]*G_1[20]+0.8660254037844387*alphaDrag[9]*G_1[19]+0.8660254037844387*alphaDrag[10]*G_1[18]+0.8660254037844387*alphaDrag[15]*G_1[17]+0.8660254037844386*alphaDrag[14]*G_1[16]+0.9682458365518543*G_1[1]*alphaDrag[15]+0.9682458365518543*alphaDrag[1]*G_1[15]+0.9682458365518543*G_1[0]*alphaDrag[14]+0.9682458365518543*alphaDrag[0]*G_1[14]+0.9682458365518543*G_1[5]*alphaDrag[13]+0.9682458365518543*alphaDrag[5]*G_1[13]+0.9682458365518543*G_1[6]*alphaDrag[12]+0.9682458365518543*alphaDrag[6]*G_1[12]+0.9682458365518543*G_1[8]*alphaDrag[11]+0.9682458365518543*alphaDrag[8]*G_1[11]+0.9682458365518543*G_1[2]*alphaDrag[10]+0.9682458365518543*alphaDrag[2]*G_1[10]+0.9682458365518543*G_1[3]*alphaDrag[9]+0.9682458365518543*alphaDrag[3]*G_1[9]+0.9682458365518543*G_1[4]*alphaDrag[7]+0.9682458365518543*alphaDrag[4]*G_1[7]; 
  out_G_1[23] += 0.8660254037844386*alphaDrag[4]*G_1[23]+0.8660254037844387*alphaDrag[8]*G_1[22]+0.8660254037844387*alphaDrag[9]*G_1[21]+0.8660254037844387*alphaDrag[10]*G_1[20]+0.8660254037844386*alphaDrag[12]*G_1[19]+0.8660254037844386*alphaDrag[13]*G_1[18]+0.8660254037844386*alphaDrag[14]*G_1[17]+0.8660254037844387*alphaDrag[15]*G_1[16]+0.9682458365518543*G_1[0]*alphaDrag[15]+0.9682458365518543*alphaDrag[0]*G_1[15]+0.9682458365518543*G_1[1]*alphaDrag[14]+0.9682458365518543*alphaDrag[1]*G_1[14]+0.9682458365518543*G_1[2]*alphaDrag[13]+0.9682458365518543*alphaDrag[2]*G_1[13]+0.9682458365518543*G_1[3]*alphaDrag[12]+0.9682458365518543*alphaDrag[3]*G_1[12]+0.9682458365518543*G_1[4]*alphaDrag[11]+0.9682458365518543*alphaDrag[4]*G_1[11]+0.9682458365518543*G_1[5]*alphaDrag[10]+0.9682458365518543*alphaDrag[5]*G_1[10]+0.9682458365518543*G_1[6]*alphaDrag[9]+0.9682458365518543*alphaDrag[6]*G_1[9]+0.9682458365518543*G_1[7]*alphaDrag[8]+0.9682458365518543*alphaDrag[7]*G_1[8]; 

  return fabs(0.625*alphaDrag[0]); 

} 
