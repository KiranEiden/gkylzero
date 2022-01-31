#include <gkyl_vlasov_lbo_kernels.h> 
GKYL_CU_DH double vlasov_lbo_drag_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[4]:      Cell-center coordinates. 
  // dxv[4]:    Cell spacing. 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // f:         Input distribution function.
  // out:       Incremented output 
  const double rdvx2 = 2.0/dxv[2]; 
  const double rdvy2 = 2.0/dxv[3]; 

  double alphaDrag[32]; 
  // Expand rdv2*(nu*vx-nuUSumx) in phase basis.
  alphaDrag[0] = (2.0*nuUSum[0]-2.0*nuSum[0]*w[2])*rdvx2; 
  alphaDrag[1] = (2.0*nuUSum[1]-2.0*nuSum[1]*w[2])*rdvx2; 
  alphaDrag[2] = (2.0*nuUSum[2]-2.0*nuSum[2]*w[2])*rdvx2; 
  alphaDrag[3] = -0.5773502691896258*nuSum[0]*dxv[2]*rdvx2; 
  alphaDrag[5] = (2.0*nuUSum[3]-2.0*w[2]*nuSum[3])*rdvx2; 
  alphaDrag[6] = -0.5773502691896258*nuSum[1]*dxv[2]*rdvx2; 
  alphaDrag[7] = -0.5773502691896258*dxv[2]*nuSum[2]*rdvx2; 
  alphaDrag[11] = -0.5773502691896258*dxv[2]*nuSum[3]*rdvx2; 

  // Expand rdv2*(nu*vy-nuUSumy) in phase basis.
  alphaDrag[16] = (2.0*nuUSum[4]-2.0*nuSum[0]*w[3])*rdvy2; 
  alphaDrag[17] = (2.0*nuUSum[5]-2.0*nuSum[1]*w[3])*rdvy2; 
  alphaDrag[18] = (2.0*nuUSum[6]-2.0*nuSum[2]*w[3])*rdvy2; 
  alphaDrag[20] = -0.5773502691896258*nuSum[0]*dxv[3]*rdvy2; 
  alphaDrag[21] = (2.0*nuUSum[7]-2.0*nuSum[3]*w[3])*rdvy2; 
  alphaDrag[24] = -0.5773502691896258*nuSum[1]*dxv[3]*rdvy2; 
  alphaDrag[25] = -0.5773502691896258*nuSum[2]*dxv[3]*rdvy2; 
  alphaDrag[28] = -0.5773502691896258*dxv[3]*nuSum[3]*rdvy2; 

  out[3] += 0.4330127018922193*(alphaDrag[11]*f[11]+alphaDrag[7]*f[7]+alphaDrag[6]*f[6]+alphaDrag[5]*f[5]+alphaDrag[3]*f[3]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[4] += 0.4330127018922193*(f[12]*alphaDrag[28]+f[9]*alphaDrag[25]+f[8]*alphaDrag[24]+f[5]*alphaDrag[21]+f[4]*alphaDrag[20]+f[2]*alphaDrag[18]+f[1]*alphaDrag[17]+f[0]*alphaDrag[16]); 
  out[6] += 0.4330127018922193*(alphaDrag[7]*f[11]+f[7]*alphaDrag[11]+alphaDrag[3]*f[6]+f[3]*alphaDrag[6]+alphaDrag[2]*f[5]+f[2]*alphaDrag[5]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[7] += 0.4330127018922193*(alphaDrag[6]*f[11]+f[6]*alphaDrag[11]+alphaDrag[3]*f[7]+f[3]*alphaDrag[7]+alphaDrag[1]*f[5]+f[1]*alphaDrag[5]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[8] += 0.4330127018922193*(f[9]*alphaDrag[28]+f[12]*alphaDrag[25]+f[4]*alphaDrag[24]+f[2]*alphaDrag[21]+f[8]*alphaDrag[20]+f[5]*alphaDrag[18]+f[0]*alphaDrag[17]+f[1]*alphaDrag[16]); 
  out[9] += 0.4330127018922193*(f[8]*alphaDrag[28]+f[4]*alphaDrag[25]+f[12]*alphaDrag[24]+f[1]*alphaDrag[21]+f[9]*alphaDrag[20]+f[0]*alphaDrag[18]+f[5]*alphaDrag[17]+f[2]*alphaDrag[16]); 
  out[10] += 0.4330127018922193*(f[15]*alphaDrag[28]+f[14]*alphaDrag[25]+f[13]*alphaDrag[24]+f[11]*alphaDrag[21]+f[10]*alphaDrag[20]+f[7]*alphaDrag[18]+f[6]*alphaDrag[17]+f[3]*alphaDrag[16]+alphaDrag[11]*f[15]+alphaDrag[7]*f[14]+alphaDrag[6]*f[13]+alphaDrag[5]*f[12]+alphaDrag[3]*f[10]+alphaDrag[2]*f[9]+alphaDrag[1]*f[8]+alphaDrag[0]*f[4]); 
  out[11] += 0.4330127018922193*(alphaDrag[3]*f[11]+f[3]*alphaDrag[11]+alphaDrag[6]*f[7]+f[6]*alphaDrag[7]+alphaDrag[0]*f[5]+f[0]*alphaDrag[5]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[12] += 0.4330127018922193*(f[4]*alphaDrag[28]+f[8]*alphaDrag[25]+f[9]*alphaDrag[24]+f[0]*alphaDrag[21]+f[12]*alphaDrag[20]+f[1]*alphaDrag[18]+f[2]*alphaDrag[17]+f[5]*alphaDrag[16]); 
  out[13] += 0.4330127018922193*(f[14]*alphaDrag[28]+f[15]*alphaDrag[25]+f[10]*alphaDrag[24]+f[7]*alphaDrag[21]+f[13]*alphaDrag[20]+f[11]*alphaDrag[18]+f[3]*alphaDrag[17]+f[6]*alphaDrag[16]+alphaDrag[7]*f[15]+alphaDrag[11]*f[14]+alphaDrag[3]*f[13]+alphaDrag[2]*f[12]+alphaDrag[6]*f[10]+alphaDrag[5]*f[9]+alphaDrag[0]*f[8]+alphaDrag[1]*f[4]); 
  out[14] += 0.4330127018922193*(f[13]*alphaDrag[28]+f[10]*alphaDrag[25]+f[15]*alphaDrag[24]+f[6]*alphaDrag[21]+f[14]*alphaDrag[20]+f[3]*alphaDrag[18]+f[11]*alphaDrag[17]+f[7]*alphaDrag[16]+alphaDrag[6]*f[15]+alphaDrag[3]*f[14]+alphaDrag[11]*f[13]+alphaDrag[1]*f[12]+alphaDrag[7]*f[10]+alphaDrag[0]*f[9]+alphaDrag[5]*f[8]+alphaDrag[2]*f[4]); 
  out[15] += 0.4330127018922193*(f[10]*alphaDrag[28]+f[13]*alphaDrag[25]+f[14]*alphaDrag[24]+f[3]*alphaDrag[21]+f[15]*alphaDrag[20]+f[6]*alphaDrag[18]+f[7]*alphaDrag[17]+f[11]*alphaDrag[16]+alphaDrag[3]*f[15]+alphaDrag[6]*f[14]+alphaDrag[7]*f[13]+alphaDrag[0]*f[12]+f[10]*alphaDrag[11]+alphaDrag[1]*f[9]+alphaDrag[2]*f[8]+f[4]*alphaDrag[5]); 

  return fabs(0.125*alphaDrag[0])+fabs(0.125*alphaDrag[16]); 

} 