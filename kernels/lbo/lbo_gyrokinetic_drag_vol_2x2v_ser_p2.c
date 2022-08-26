#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH double lbo_gyrokinetic_drag_vol_2x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[4]:      cell-center coordinates. 
  // dxv[4]:    cell spacing. 
  // m_:        species mass.
  // bmag_inv:  1/(magnetic field magnitude). 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // f:         input distribution function.
  // out:       incremented output 
  double rdv2[2]; 
  rdv2[0]   = 2.0/dxv[2]; 
  rdv2[1]   = 2.0/dxv[3]; 

  double alphaDrag[96]; 
  // Expand rdv2*(nu*vpar-nuUparSum) in phase basis.
  alphaDrag[0] = rdv2[0]*(2.0*nuUSum[0]-2.0*nuSum[0]*w[2]); 
  alphaDrag[1] = rdv2[0]*(2.0*nuUSum[1]-2.0*nuSum[1]*w[2]); 
  alphaDrag[2] = rdv2[0]*(2.0*nuUSum[2]-2.0*nuSum[2]*w[2]); 
  alphaDrag[3] = -0.5773502691896258*nuSum[0]*rdv2[0]*dxv[2]; 
  alphaDrag[5] = rdv2[0]*(2.0*nuUSum[3]-2.0*w[2]*nuSum[3]); 
  alphaDrag[6] = -0.5773502691896258*rdv2[0]*nuSum[1]*dxv[2]; 
  alphaDrag[7] = -0.5773502691896258*rdv2[0]*dxv[2]*nuSum[2]; 
  alphaDrag[11] = rdv2[0]*(2.0*nuUSum[4]-2.0*w[2]*nuSum[4]); 
  alphaDrag[12] = rdv2[0]*(2.0*nuUSum[5]-2.0*w[2]*nuSum[5]); 
  alphaDrag[15] = -0.5773502691896258*rdv2[0]*dxv[2]*nuSum[3]; 
  alphaDrag[19] = rdv2[0]*(2.0*nuUSum[6]-2.0*w[2]*nuSum[6]); 
  alphaDrag[20] = rdv2[0]*(2.0*nuUSum[7]-2.0*w[2]*nuSum[7]); 
  alphaDrag[21] = -0.5773502691896258*rdv2[0]*dxv[2]*nuSum[4]; 
  alphaDrag[22] = -0.5773502691896258*rdv2[0]*dxv[2]*nuSum[5]; 
  alphaDrag[32] = -0.5773502691896258*rdv2[0]*dxv[2]*nuSum[6]; 
  alphaDrag[33] = -0.5773502691896258*rdv2[0]*dxv[2]*nuSum[7]; 

  // Expand rdv2*nu*2*mu in phase basis.
  alphaDrag[48] = -4.0*nuSum[0]*rdv2[1]*w[3]; 
  alphaDrag[49] = -4.0*nuSum[1]*rdv2[1]*w[3]; 
  alphaDrag[50] = -4.0*rdv2[1]*nuSum[2]*w[3]; 
  alphaDrag[52] = -1.154700538379252*nuSum[0]*rdv2[1]*dxv[3]; 
  alphaDrag[53] = -4.0*rdv2[1]*nuSum[3]*w[3]; 
  alphaDrag[56] = -1.154700538379252*nuSum[1]*rdv2[1]*dxv[3]; 
  alphaDrag[57] = -1.154700538379252*rdv2[1]*nuSum[2]*dxv[3]; 
  alphaDrag[59] = -4.0*rdv2[1]*w[3]*nuSum[4]; 
  alphaDrag[60] = -4.0*rdv2[1]*w[3]*nuSum[5]; 
  alphaDrag[64] = -1.154700538379252*rdv2[1]*dxv[3]*nuSum[3]; 
  alphaDrag[67] = -4.0*rdv2[1]*w[3]*nuSum[6]; 
  alphaDrag[68] = -4.0*rdv2[1]*w[3]*nuSum[7]; 
  alphaDrag[73] = -1.154700538379252*rdv2[1]*dxv[3]*nuSum[4]; 
  alphaDrag[74] = -1.154700538379252*rdv2[1]*dxv[3]*nuSum[5]; 
  alphaDrag[83] = -1.154700538379252*rdv2[1]*dxv[3]*nuSum[6]; 
  alphaDrag[84] = -1.154700538379252*rdv2[1]*dxv[3]*nuSum[7]; 

  out[3] += 0.4330127018922193*(alphaDrag[33]*f[33]+alphaDrag[32]*f[32]+alphaDrag[22]*f[22]+alphaDrag[21]*f[21]+alphaDrag[20]*f[20]+alphaDrag[19]*f[19]+alphaDrag[15]*f[15]+alphaDrag[12]*f[12]+alphaDrag[11]*f[11]+alphaDrag[7]*f[7]+alphaDrag[6]*f[6]+alphaDrag[5]*f[5]+alphaDrag[3]*f[3]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[4] += 0.4330127018922193*(f[36]*alphaDrag[84]+f[35]*alphaDrag[83]+f[26]*alphaDrag[74]+f[25]*alphaDrag[73]+f[20]*alphaDrag[68]+f[19]*alphaDrag[67]+f[16]*alphaDrag[64]+f[12]*alphaDrag[60]+f[11]*alphaDrag[59]+f[9]*alphaDrag[57]+f[8]*alphaDrag[56]+f[5]*alphaDrag[53]+f[4]*alphaDrag[52]+f[2]*alphaDrag[50]+f[1]*alphaDrag[49]+f[0]*alphaDrag[48]); 
  out[6] += 0.4330127018922193*(alphaDrag[22]*f[33]+f[22]*alphaDrag[33])+0.3872983346207416*(alphaDrag[15]*f[32]+f[15]*alphaDrag[32]+alphaDrag[6]*f[21]+f[6]*alphaDrag[21])+0.4330127018922193*(alphaDrag[12]*f[20]+f[12]*alphaDrag[20])+0.3872983346207416*(alphaDrag[5]*f[19]+f[5]*alphaDrag[19])+0.4330127018922193*(alphaDrag[7]*f[15]+f[7]*alphaDrag[15])+0.3872983346207416*(alphaDrag[1]*f[11]+f[1]*alphaDrag[11])+0.4330127018922193*(alphaDrag[3]*f[6]+f[3]*alphaDrag[6]+alphaDrag[2]*f[5]+f[2]*alphaDrag[5]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[7] += 0.3872983346207416*(alphaDrag[15]*f[33]+f[15]*alphaDrag[33])+0.4330127018922193*(alphaDrag[21]*f[32]+f[21]*alphaDrag[32])+0.3872983346207416*(alphaDrag[7]*f[22]+f[7]*alphaDrag[22]+alphaDrag[5]*f[20]+f[5]*alphaDrag[20])+0.4330127018922193*(alphaDrag[11]*f[19]+f[11]*alphaDrag[19]+alphaDrag[6]*f[15]+f[6]*alphaDrag[15])+0.3872983346207416*(alphaDrag[2]*f[12]+f[2]*alphaDrag[12])+0.4330127018922193*(alphaDrag[3]*f[7]+f[3]*alphaDrag[7]+alphaDrag[1]*f[5]+f[1]*alphaDrag[5]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[8] += 0.4330127018922193*f[26]*alphaDrag[84]+0.3872983346207416*f[16]*alphaDrag[83]+0.4330127018922193*f[36]*alphaDrag[74]+0.3872983346207416*f[8]*alphaDrag[73]+0.4330127018922193*f[12]*alphaDrag[68]+0.3872983346207416*(f[5]*alphaDrag[67]+f[35]*alphaDrag[64])+0.4330127018922193*(f[9]*alphaDrag[64]+f[20]*alphaDrag[60])+0.3872983346207416*f[1]*alphaDrag[59]+0.4330127018922193*f[16]*alphaDrag[57]+(0.3872983346207416*f[25]+0.4330127018922193*f[4])*alphaDrag[56]+0.3872983346207416*f[19]*alphaDrag[53]+0.4330127018922193*(f[2]*alphaDrag[53]+f[8]*alphaDrag[52]+f[5]*alphaDrag[50])+0.3872983346207416*f[11]*alphaDrag[49]+0.4330127018922193*(f[0]*alphaDrag[49]+f[1]*alphaDrag[48]); 
  out[9] += 0.3872983346207416*f[16]*alphaDrag[84]+0.4330127018922193*f[25]*alphaDrag[83]+0.3872983346207416*f[9]*alphaDrag[74]+0.4330127018922193*f[35]*alphaDrag[73]+0.3872983346207416*f[5]*alphaDrag[68]+0.4330127018922193*f[11]*alphaDrag[67]+(0.3872983346207416*f[36]+0.4330127018922193*f[8])*alphaDrag[64]+0.3872983346207416*f[2]*alphaDrag[60]+0.4330127018922193*f[19]*alphaDrag[59]+0.3872983346207416*f[26]*alphaDrag[57]+0.4330127018922193*(f[4]*alphaDrag[57]+f[16]*alphaDrag[56])+0.3872983346207416*f[20]*alphaDrag[53]+0.4330127018922193*(f[1]*alphaDrag[53]+f[9]*alphaDrag[52])+0.3872983346207416*f[12]*alphaDrag[50]+0.4330127018922193*(f[0]*alphaDrag[50]+f[5]*alphaDrag[49]+f[2]*alphaDrag[48]); 
  out[10] += 0.4330127018922193*(f[45]*alphaDrag[84]+f[44]*alphaDrag[83]+f[38]*alphaDrag[74]+f[37]*alphaDrag[73]+f[33]*alphaDrag[68]+f[32]*alphaDrag[67]+f[31]*alphaDrag[64]+f[22]*alphaDrag[60]+f[21]*alphaDrag[59]+f[18]*alphaDrag[57]+f[17]*alphaDrag[56]+f[15]*alphaDrag[53]+f[10]*alphaDrag[52]+f[7]*alphaDrag[50]+f[6]*alphaDrag[49]+f[3]*alphaDrag[48]+alphaDrag[33]*f[45]+alphaDrag[32]*f[44]+alphaDrag[22]*f[38]+alphaDrag[21]*f[37]+alphaDrag[20]*f[36]+alphaDrag[19]*f[35]+alphaDrag[15]*f[31]+alphaDrag[12]*f[26]+alphaDrag[11]*f[25]+alphaDrag[7]*f[18]+alphaDrag[6]*f[17]+alphaDrag[5]*f[16]+alphaDrag[3]*f[10]+alphaDrag[2]*f[9]+alphaDrag[1]*f[8]+alphaDrag[0]*f[4]); 
  out[13] += 0.8660254037844386*alphaDrag[15]*f[34]+0.9682458365518543*(alphaDrag[20]*f[33]+f[20]*alphaDrag[33]+alphaDrag[19]*f[32]+f[19]*alphaDrag[32])+0.8660254037844386*(alphaDrag[7]*f[24]+alphaDrag[6]*f[23])+0.9682458365518543*(alphaDrag[12]*f[22]+f[12]*alphaDrag[22]+alphaDrag[11]*f[21]+f[11]*alphaDrag[21]+alphaDrag[5]*f[15]+f[5]*alphaDrag[15])+0.8660254037844386*alphaDrag[3]*f[13]+0.9682458365518543*(alphaDrag[2]*f[7]+f[2]*alphaDrag[7]+alphaDrag[1]*f[6]+f[1]*alphaDrag[6]+alphaDrag[0]*f[3]+f[0]*alphaDrag[3]); 
  out[14] += 0.9682458365518543*(f[20]*alphaDrag[84]+f[19]*alphaDrag[83]+f[12]*alphaDrag[74]+f[11]*alphaDrag[73]+f[36]*alphaDrag[68]+f[35]*alphaDrag[67])+0.8660254037844386*f[41]*alphaDrag[64]+0.9682458365518543*(f[5]*alphaDrag[64]+f[26]*alphaDrag[60]+f[25]*alphaDrag[59])+(0.8660254037844386*f[29]+0.9682458365518543*f[2])*alphaDrag[57]+0.8660254037844386*f[28]*alphaDrag[56]+0.9682458365518543*(f[1]*alphaDrag[56]+f[16]*alphaDrag[53])+0.8660254037844386*f[14]*alphaDrag[52]+0.9682458365518543*(f[0]*alphaDrag[52]+f[9]*alphaDrag[50]+f[8]*alphaDrag[49]+f[4]*alphaDrag[48]); 
  out[15] += (0.3464101615137755*alphaDrag[32]+0.3872983346207416*alphaDrag[7])*f[33]+0.3464101615137755*f[32]*alphaDrag[33]+0.3872983346207416*(f[7]*alphaDrag[33]+alphaDrag[6]*f[32]+f[6]*alphaDrag[32]+alphaDrag[15]*f[22]+f[15]*alphaDrag[22]+alphaDrag[15]*f[21]+f[15]*alphaDrag[21])+(0.3464101615137755*alphaDrag[19]+0.3872983346207416*alphaDrag[2])*f[20]+0.3464101615137755*f[19]*alphaDrag[20]+0.3872983346207416*(f[2]*alphaDrag[20]+alphaDrag[1]*f[19]+f[1]*alphaDrag[19])+0.4330127018922193*(alphaDrag[3]*f[15]+f[3]*alphaDrag[15])+0.3872983346207416*(alphaDrag[5]*f[12]+f[5]*alphaDrag[12]+alphaDrag[5]*f[11]+f[5]*alphaDrag[11])+0.4330127018922193*(alphaDrag[6]*f[7]+f[6]*alphaDrag[7]+alphaDrag[0]*f[5]+f[0]*alphaDrag[5]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[16] += (0.3464101615137755*f[35]+0.3872983346207416*f[9])*alphaDrag[84]+0.3464101615137755*f[36]*alphaDrag[83]+0.3872983346207416*(f[8]*alphaDrag[83]+f[16]*(alphaDrag[74]+alphaDrag[73]))+(0.3464101615137755*f[19]+0.3872983346207416*f[2])*alphaDrag[68]+(0.3464101615137755*f[20]+0.3872983346207416*f[1])*alphaDrag[67]+(0.3872983346207416*(f[26]+f[25])+0.4330127018922193*f[4])*alphaDrag[64]+0.3872983346207416*f[5]*(alphaDrag[60]+alphaDrag[59])+(0.3872983346207416*f[36]+0.4330127018922193*f[8])*alphaDrag[57]+(0.3872983346207416*f[35]+0.4330127018922193*f[9])*alphaDrag[56]+0.3872983346207416*(f[12]+f[11])*alphaDrag[53]+0.4330127018922193*(f[0]*alphaDrag[53]+f[16]*alphaDrag[52])+(0.3872983346207416*f[20]+0.4330127018922193*f[1])*alphaDrag[50]+0.3872983346207416*f[19]*alphaDrag[49]+0.4330127018922193*(f[2]*alphaDrag[49]+f[5]*alphaDrag[48]); 
  out[17] += 0.4330127018922193*f[38]*alphaDrag[84]+0.3872983346207416*f[31]*alphaDrag[83]+0.4330127018922193*f[45]*alphaDrag[74]+0.3872983346207416*f[17]*alphaDrag[73]+0.4330127018922193*f[22]*alphaDrag[68]+0.3872983346207416*(f[15]*alphaDrag[67]+f[44]*alphaDrag[64])+0.4330127018922193*(f[18]*alphaDrag[64]+f[33]*alphaDrag[60])+0.3872983346207416*f[6]*alphaDrag[59]+0.4330127018922193*f[31]*alphaDrag[57]+(0.3872983346207416*f[37]+0.4330127018922193*f[10])*alphaDrag[56]+0.3872983346207416*f[32]*alphaDrag[53]+0.4330127018922193*(f[7]*alphaDrag[53]+f[17]*alphaDrag[52]+f[15]*alphaDrag[50])+0.3872983346207416*f[21]*alphaDrag[49]+0.4330127018922193*(f[3]*alphaDrag[49]+f[6]*alphaDrag[48]+alphaDrag[22]*f[45])+0.3872983346207416*alphaDrag[15]*f[44]+0.4330127018922193*alphaDrag[33]*f[38]+0.3872983346207416*alphaDrag[6]*f[37]+0.4330127018922193*alphaDrag[12]*f[36]+0.3872983346207416*(alphaDrag[5]*f[35]+f[31]*alphaDrag[32])+0.4330127018922193*(alphaDrag[7]*f[31]+alphaDrag[20]*f[26])+0.3872983346207416*(alphaDrag[1]*f[25]+f[17]*alphaDrag[21]+f[16]*alphaDrag[19])+0.4330127018922193*(alphaDrag[15]*f[18]+alphaDrag[3]*f[17]+alphaDrag[2]*f[16])+0.3872983346207416*f[8]*alphaDrag[11]+0.4330127018922193*(alphaDrag[6]*f[10]+alphaDrag[5]*f[9]+alphaDrag[0]*f[8]+alphaDrag[1]*f[4]); 
  out[18] += 0.3872983346207416*f[31]*alphaDrag[84]+0.4330127018922193*f[37]*alphaDrag[83]+0.3872983346207416*f[18]*alphaDrag[74]+0.4330127018922193*f[44]*alphaDrag[73]+0.3872983346207416*f[15]*alphaDrag[68]+0.4330127018922193*f[21]*alphaDrag[67]+(0.3872983346207416*f[45]+0.4330127018922193*f[17])*alphaDrag[64]+0.3872983346207416*f[7]*alphaDrag[60]+0.4330127018922193*f[32]*alphaDrag[59]+0.3872983346207416*f[38]*alphaDrag[57]+0.4330127018922193*(f[10]*alphaDrag[57]+f[31]*alphaDrag[56])+0.3872983346207416*f[33]*alphaDrag[53]+0.4330127018922193*(f[6]*alphaDrag[53]+f[18]*alphaDrag[52])+0.3872983346207416*f[22]*alphaDrag[50]+0.4330127018922193*(f[3]*alphaDrag[50]+f[15]*alphaDrag[49]+f[7]*alphaDrag[48])+0.3872983346207416*alphaDrag[15]*f[45]+0.4330127018922193*alphaDrag[21]*f[44]+0.3872983346207416*alphaDrag[7]*f[38]+0.4330127018922193*alphaDrag[32]*f[37]+0.3872983346207416*alphaDrag[5]*f[36]+0.4330127018922193*alphaDrag[11]*f[35]+f[31]*(0.3872983346207416*alphaDrag[33]+0.4330127018922193*alphaDrag[6])+0.3872983346207416*alphaDrag[2]*f[26]+0.4330127018922193*alphaDrag[19]*f[25]+0.3872983346207416*(f[18]*alphaDrag[22]+f[16]*alphaDrag[20])+0.4330127018922193*(alphaDrag[3]*f[18]+alphaDrag[15]*f[17]+alphaDrag[1]*f[16])+0.3872983346207416*f[9]*alphaDrag[12]+0.4330127018922193*(alphaDrag[7]*f[10]+alphaDrag[0]*f[9]+alphaDrag[5]*f[8]+alphaDrag[2]*f[4]); 
  out[21] += 0.3872983346207416*alphaDrag[33]*f[33]+0.276641667586244*alphaDrag[32]*f[32]+0.4330127018922193*(alphaDrag[7]*f[32]+f[7]*alphaDrag[32])+0.276641667586244*alphaDrag[21]*f[21]+0.4330127018922193*(alphaDrag[3]*f[21]+f[3]*alphaDrag[21])+0.3872983346207416*alphaDrag[20]*f[20]+0.276641667586244*alphaDrag[19]*f[19]+0.4330127018922193*(alphaDrag[2]*f[19]+f[2]*alphaDrag[19])+0.3872983346207416*alphaDrag[15]*f[15]+0.276641667586244*alphaDrag[11]*f[11]+0.4330127018922193*(alphaDrag[0]*f[11]+f[0]*alphaDrag[11])+0.3872983346207416*(alphaDrag[6]*f[6]+alphaDrag[5]*f[5]+alphaDrag[1]*f[1]); 
  out[22] += 0.276641667586244*alphaDrag[33]*f[33]+0.4330127018922193*(alphaDrag[6]*f[33]+f[6]*alphaDrag[33])+0.3872983346207416*alphaDrag[32]*f[32]+0.276641667586244*alphaDrag[22]*f[22]+0.4330127018922193*(alphaDrag[3]*f[22]+f[3]*alphaDrag[22])+0.276641667586244*alphaDrag[20]*f[20]+0.4330127018922193*(alphaDrag[1]*f[20]+f[1]*alphaDrag[20])+0.3872983346207416*(alphaDrag[19]*f[19]+alphaDrag[15]*f[15])+0.276641667586244*alphaDrag[12]*f[12]+0.4330127018922193*(alphaDrag[0]*f[12]+f[0]*alphaDrag[12])+0.3872983346207416*(alphaDrag[7]*f[7]+alphaDrag[5]*f[5]+alphaDrag[2]*f[2]); 
  out[23] += (0.7745966692414833*alphaDrag[32]+0.8660254037844386*alphaDrag[7])*f[34]+0.9682458365518543*(alphaDrag[12]*f[33]+f[12]*alphaDrag[33])+0.8660254037844386*(alphaDrag[5]*f[32]+f[5]*alphaDrag[32]+alphaDrag[15]*f[24])+(0.7745966692414833*alphaDrag[21]+0.8660254037844386*alphaDrag[3])*f[23]+0.9682458365518543*(alphaDrag[20]*f[22]+f[20]*alphaDrag[22])+0.8660254037844386*(alphaDrag[1]*f[21]+f[1]*alphaDrag[21]+alphaDrag[15]*f[19]+f[15]*alphaDrag[19])+0.9682458365518543*(alphaDrag[2]*f[15]+f[2]*alphaDrag[15])+0.8660254037844386*(alphaDrag[6]*(f[13]+f[11])+f[6]*alphaDrag[11])+0.9682458365518543*(alphaDrag[5]*f[7]+f[5]*alphaDrag[7]+alphaDrag[0]*f[6]+f[0]*alphaDrag[6]+alphaDrag[1]*f[3]+f[1]*alphaDrag[3]); 
  out[24] += 0.7745966692414833*alphaDrag[33]*f[34]+0.8660254037844386*(alphaDrag[6]*f[34]+alphaDrag[5]*f[33]+f[5]*alphaDrag[33])+0.9682458365518543*(alphaDrag[11]*f[32]+f[11]*alphaDrag[32])+0.7745966692414833*alphaDrag[22]*f[24]+0.8660254037844386*(alphaDrag[3]*f[24]+alphaDrag[15]*f[23]+alphaDrag[2]*f[22]+f[2]*alphaDrag[22])+0.9682458365518543*(alphaDrag[19]*f[21]+f[19]*alphaDrag[21])+0.8660254037844386*(alphaDrag[15]*f[20]+f[15]*alphaDrag[20])+0.9682458365518543*(alphaDrag[1]*f[15]+f[1]*alphaDrag[15])+0.8660254037844386*(alphaDrag[7]*(f[13]+f[12])+f[7]*alphaDrag[12])+0.9682458365518543*(alphaDrag[0]*f[7]+f[0]*alphaDrag[7]+alphaDrag[5]*f[6]+f[5]*alphaDrag[6]+alphaDrag[2]*f[3]+f[2]*alphaDrag[3]); 
  out[25] += 0.3872983346207416*f[36]*alphaDrag[84]+(0.276641667586244*f[35]+0.4330127018922193*f[9])*alphaDrag[83]+(0.276641667586244*f[25]+0.4330127018922193*f[4])*alphaDrag[73]+0.3872983346207416*f[20]*alphaDrag[68]+(0.276641667586244*f[19]+0.4330127018922193*f[2])*alphaDrag[67]+0.3872983346207416*f[16]*alphaDrag[64]+0.276641667586244*f[11]*alphaDrag[59]+0.4330127018922193*(f[0]*alphaDrag[59]+f[35]*alphaDrag[57])+0.3872983346207416*(f[8]*alphaDrag[56]+f[5]*alphaDrag[53])+0.4330127018922193*(f[25]*alphaDrag[52]+f[19]*alphaDrag[50])+0.3872983346207416*f[1]*alphaDrag[49]+0.4330127018922193*f[11]*alphaDrag[48]; 
  out[26] += (0.276641667586244*f[36]+0.4330127018922193*f[8])*alphaDrag[84]+0.3872983346207416*f[35]*alphaDrag[83]+(0.276641667586244*f[26]+0.4330127018922193*f[4])*alphaDrag[74]+(0.276641667586244*f[20]+0.4330127018922193*f[1])*alphaDrag[68]+0.3872983346207416*(f[19]*alphaDrag[67]+f[16]*alphaDrag[64])+(0.276641667586244*f[12]+0.4330127018922193*f[0])*alphaDrag[60]+0.3872983346207416*f[9]*alphaDrag[57]+0.4330127018922193*f[36]*alphaDrag[56]+0.3872983346207416*f[5]*alphaDrag[53]+0.4330127018922193*f[26]*alphaDrag[52]+0.3872983346207416*f[2]*alphaDrag[50]+0.4330127018922193*(f[20]*alphaDrag[49]+f[12]*alphaDrag[48]); 
  out[27] += 0.4330127018922193*(f[46]*alphaDrag[64]+f[40]*alphaDrag[57]+f[39]*alphaDrag[56]+f[34]*alphaDrag[53]+f[27]*alphaDrag[52]+f[24]*alphaDrag[50]+f[23]*alphaDrag[49]+f[13]*alphaDrag[48])+0.8660254037844386*alphaDrag[15]*f[46]+0.9682458365518543*(alphaDrag[20]*f[45]+alphaDrag[19]*f[44])+0.8660254037844386*(alphaDrag[7]*f[40]+alphaDrag[6]*f[39])+0.9682458365518543*(alphaDrag[12]*f[38]+alphaDrag[11]*f[37]+alphaDrag[33]*f[36]+alphaDrag[32]*f[35]+alphaDrag[5]*f[31])+0.8660254037844386*alphaDrag[3]*f[27]+0.9682458365518543*(alphaDrag[22]*f[26]+alphaDrag[21]*f[25]+alphaDrag[2]*f[18]+alphaDrag[1]*f[17]+alphaDrag[15]*f[16]+alphaDrag[0]*f[10]+alphaDrag[7]*f[9]+alphaDrag[6]*f[8]+alphaDrag[3]*f[4]); 
  out[28] += 0.9682458365518543*f[12]*alphaDrag[84]+(0.7745966692414833*f[41]+0.8660254037844386*f[5])*alphaDrag[83]+0.9682458365518543*f[20]*alphaDrag[74]+(0.7745966692414833*f[28]+0.8660254037844386*f[1])*alphaDrag[73]+0.9682458365518543*f[26]*alphaDrag[68]+0.8660254037844386*(f[16]*alphaDrag[67]+(f[29]+f[19])*alphaDrag[64])+0.9682458365518543*(f[2]*alphaDrag[64]+f[36]*alphaDrag[60])+0.8660254037844386*f[8]*alphaDrag[59]+(0.8660254037844386*f[41]+0.9682458365518543*f[5])*alphaDrag[57]+(0.8660254037844386*(f[14]+f[11])+0.9682458365518543*f[0])*alphaDrag[56]+(0.8660254037844386*f[35]+0.9682458365518543*f[9])*alphaDrag[53]+0.8660254037844386*f[28]*alphaDrag[52]+0.9682458365518543*(f[1]*alphaDrag[52]+f[16]*alphaDrag[50])+0.8660254037844386*f[25]*alphaDrag[49]+0.9682458365518543*(f[4]*alphaDrag[49]+f[8]*alphaDrag[48]); 
  out[29] += (0.7745966692414833*f[41]+0.8660254037844386*f[5])*alphaDrag[84]+0.9682458365518543*f[11]*alphaDrag[83]+(0.7745966692414833*f[29]+0.8660254037844386*f[2])*alphaDrag[74]+0.9682458365518543*f[19]*alphaDrag[73]+0.8660254037844386*f[16]*alphaDrag[68]+0.9682458365518543*f[25]*alphaDrag[67]+(0.8660254037844386*(f[28]+f[20])+0.9682458365518543*f[1])*alphaDrag[64]+0.8660254037844386*f[9]*alphaDrag[60]+0.9682458365518543*f[35]*alphaDrag[59]+(0.8660254037844386*(f[14]+f[12])+0.9682458365518543*f[0])*alphaDrag[57]+(0.8660254037844386*f[41]+0.9682458365518543*f[5])*alphaDrag[56]+(0.8660254037844386*f[36]+0.9682458365518543*f[8])*alphaDrag[53]+(0.8660254037844386*f[29]+0.9682458365518543*f[2])*alphaDrag[52]+0.8660254037844386*f[26]*alphaDrag[50]+0.9682458365518543*(f[4]*alphaDrag[50]+f[16]*alphaDrag[49]+f[9]*alphaDrag[48]); 
  out[30] += 0.9682458365518543*(f[33]*alphaDrag[84]+f[32]*alphaDrag[83]+f[22]*alphaDrag[74]+f[21]*alphaDrag[73]+f[45]*alphaDrag[68]+f[44]*alphaDrag[67])+0.8660254037844386*f[47]*alphaDrag[64]+0.9682458365518543*(f[15]*alphaDrag[64]+f[38]*alphaDrag[60]+f[37]*alphaDrag[59])+(0.8660254037844386*f[43]+0.9682458365518543*f[7])*alphaDrag[57]+0.8660254037844386*f[42]*alphaDrag[56]+0.9682458365518543*(f[6]*alphaDrag[56]+f[31]*alphaDrag[53])+0.8660254037844386*f[30]*alphaDrag[52]+0.9682458365518543*(f[3]*alphaDrag[52]+f[18]*alphaDrag[50]+f[17]*alphaDrag[49]+f[10]*alphaDrag[48])+0.4330127018922193*(alphaDrag[15]*f[47]+alphaDrag[7]*f[43]+alphaDrag[6]*f[42]+alphaDrag[5]*f[41]+alphaDrag[3]*f[30]+alphaDrag[2]*f[29]+alphaDrag[1]*f[28]+alphaDrag[0]*f[14]); 
  out[31] += (0.3464101615137755*f[44]+0.3872983346207416*f[18])*alphaDrag[84]+0.3464101615137755*f[45]*alphaDrag[83]+0.3872983346207416*(f[17]*alphaDrag[83]+f[31]*(alphaDrag[74]+alphaDrag[73]))+(0.3464101615137755*f[32]+0.3872983346207416*f[7])*alphaDrag[68]+(0.3464101615137755*f[33]+0.3872983346207416*f[6])*alphaDrag[67]+(0.3872983346207416*(f[38]+f[37])+0.4330127018922193*f[10])*alphaDrag[64]+0.3872983346207416*f[15]*(alphaDrag[60]+alphaDrag[59])+(0.3872983346207416*f[45]+0.4330127018922193*f[17])*alphaDrag[57]+(0.3872983346207416*f[44]+0.4330127018922193*f[18])*alphaDrag[56]+0.3872983346207416*(f[22]+f[21])*alphaDrag[53]+0.4330127018922193*(f[3]*alphaDrag[53]+f[31]*alphaDrag[52])+(0.3872983346207416*f[33]+0.4330127018922193*f[6])*alphaDrag[50]+0.3872983346207416*f[32]*alphaDrag[49]+0.4330127018922193*(f[7]*alphaDrag[49]+f[15]*alphaDrag[48])+(0.3464101615137755*alphaDrag[32]+0.3872983346207416*alphaDrag[7])*f[45]+0.3464101615137755*alphaDrag[33]*f[44]+0.3872983346207416*(alphaDrag[6]*f[44]+alphaDrag[15]*(f[38]+f[37]))+(0.3464101615137755*alphaDrag[19]+0.3872983346207416*alphaDrag[2])*f[36]+0.3464101615137755*alphaDrag[20]*f[35]+0.3872983346207416*(alphaDrag[1]*f[35]+f[18]*alphaDrag[33]+f[17]*alphaDrag[32])+(0.3872983346207416*(alphaDrag[22]+alphaDrag[21])+0.4330127018922193*alphaDrag[3])*f[31]+0.3872983346207416*(alphaDrag[5]*(f[26]+f[25])+f[9]*alphaDrag[20]+f[8]*alphaDrag[19])+0.4330127018922193*(alphaDrag[6]*f[18]+alphaDrag[7]*f[17])+0.3872983346207416*(alphaDrag[12]+alphaDrag[11])*f[16]+0.4330127018922193*(alphaDrag[0]*f[16]+f[10]*alphaDrag[15]+alphaDrag[1]*f[9]+alphaDrag[2]*f[8]+f[4]*alphaDrag[5]); 
  out[32] += 0.3464101615137755*(alphaDrag[15]*f[33]+f[15]*alphaDrag[33])+(0.3872983346207416*alphaDrag[22]+0.276641667586244*alphaDrag[21]+0.4330127018922193*alphaDrag[3])*f[32]+(0.3872983346207416*f[22]+0.276641667586244*f[21])*alphaDrag[32]+0.4330127018922193*(f[3]*alphaDrag[32]+alphaDrag[7]*f[21]+f[7]*alphaDrag[21])+0.3464101615137755*(alphaDrag[5]*f[20]+f[5]*alphaDrag[20])+(0.3872983346207416*alphaDrag[12]+0.276641667586244*alphaDrag[11]+0.4330127018922193*alphaDrag[0])*f[19]+(0.3872983346207416*f[12]+0.276641667586244*f[11]+0.4330127018922193*f[0])*alphaDrag[19]+0.3872983346207416*(alphaDrag[6]*f[15]+f[6]*alphaDrag[15])+0.4330127018922193*(alphaDrag[2]*f[11]+f[2]*alphaDrag[11])+0.3872983346207416*(alphaDrag[1]*f[5]+f[1]*alphaDrag[5]); 
  out[33] += (0.276641667586244*alphaDrag[22]+0.3872983346207416*alphaDrag[21]+0.4330127018922193*alphaDrag[3])*f[33]+(0.276641667586244*f[22]+0.3872983346207416*f[21]+0.4330127018922193*f[3])*alphaDrag[33]+0.3464101615137755*(alphaDrag[15]*f[32]+f[15]*alphaDrag[32])+0.4330127018922193*(alphaDrag[6]*f[22]+f[6]*alphaDrag[22])+(0.276641667586244*alphaDrag[12]+0.3872983346207416*alphaDrag[11]+0.4330127018922193*alphaDrag[0])*f[20]+(0.276641667586244*f[12]+0.3872983346207416*f[11]+0.4330127018922193*f[0])*alphaDrag[20]+0.3464101615137755*(alphaDrag[5]*f[19]+f[5]*alphaDrag[19])+0.3872983346207416*(alphaDrag[7]*f[15]+f[7]*alphaDrag[15])+0.4330127018922193*(alphaDrag[1]*f[12]+f[1]*alphaDrag[12])+0.3872983346207416*(alphaDrag[2]*f[5]+f[2]*alphaDrag[5]); 
  out[34] += (0.7745966692414833*(alphaDrag[22]+alphaDrag[21])+0.8660254037844386*alphaDrag[3])*f[34]+(0.7745966692414833*alphaDrag[19]+0.8660254037844386*alphaDrag[2])*f[33]+(0.7745966692414833*(f[24]+f[19])+0.8660254037844386*f[2])*alphaDrag[33]+(0.7745966692414833*alphaDrag[20]+0.8660254037844386*alphaDrag[1])*f[32]+0.7745966692414833*(f[23]+f[20])*alphaDrag[32]+0.8660254037844386*(f[1]*alphaDrag[32]+alphaDrag[6]*f[24]+alphaDrag[7]*f[23]+alphaDrag[5]*f[22]+f[5]*alphaDrag[22]+alphaDrag[5]*f[21]+f[5]*alphaDrag[21]+alphaDrag[7]*f[20]+f[7]*alphaDrag[20]+alphaDrag[6]*f[19]+f[6]*alphaDrag[19])+(0.8660254037844386*(alphaDrag[12]+alphaDrag[11])+0.9682458365518543*alphaDrag[0])*f[15]+0.8660254037844386*(f[13]+f[12]+f[11])*alphaDrag[15]+0.9682458365518543*(f[0]*alphaDrag[15]+alphaDrag[1]*f[7]+f[1]*alphaDrag[7]+alphaDrag[2]*f[6]+f[2]*alphaDrag[6]+alphaDrag[3]*f[5]+f[3]*alphaDrag[5]); 
  out[35] += 0.3464101615137755*f[16]*alphaDrag[84]+(0.3872983346207416*f[26]+0.276641667586244*f[25]+0.4330127018922193*f[4])*alphaDrag[83]+0.3872983346207416*f[35]*alphaDrag[74]+(0.276641667586244*f[35]+0.4330127018922193*f[9])*alphaDrag[73]+0.3464101615137755*f[5]*alphaDrag[68]+(0.3872983346207416*f[12]+0.276641667586244*f[11]+0.4330127018922193*f[0])*alphaDrag[67]+(0.3464101615137755*f[36]+0.3872983346207416*f[8])*alphaDrag[64]+f[19]*(0.3872983346207416*alphaDrag[60]+0.276641667586244*alphaDrag[59])+0.4330127018922193*(f[2]*alphaDrag[59]+f[25]*alphaDrag[57])+0.3872983346207416*f[16]*alphaDrag[56]+(0.3464101615137755*f[20]+0.3872983346207416*f[1])*alphaDrag[53]+0.4330127018922193*(f[35]*alphaDrag[52]+f[11]*alphaDrag[50])+0.3872983346207416*f[5]*alphaDrag[49]+0.4330127018922193*f[19]*alphaDrag[48]; 
  out[36] += (0.276641667586244*f[26]+0.3872983346207416*f[25]+0.4330127018922193*f[4])*alphaDrag[84]+0.3464101615137755*f[16]*alphaDrag[83]+(0.276641667586244*f[36]+0.4330127018922193*f[8])*alphaDrag[74]+0.3872983346207416*f[36]*alphaDrag[73]+(0.276641667586244*f[12]+0.3872983346207416*f[11]+0.4330127018922193*f[0])*alphaDrag[68]+0.3464101615137755*f[5]*alphaDrag[67]+(0.3464101615137755*f[35]+0.3872983346207416*f[9])*alphaDrag[64]+(0.276641667586244*f[20]+0.4330127018922193*f[1])*alphaDrag[60]+0.3872983346207416*(f[20]*alphaDrag[59]+f[16]*alphaDrag[57])+0.4330127018922193*f[26]*alphaDrag[56]+(0.3464101615137755*f[19]+0.3872983346207416*f[2])*alphaDrag[53]+0.4330127018922193*f[36]*alphaDrag[52]+0.3872983346207416*f[5]*alphaDrag[50]+0.4330127018922193*(f[12]*alphaDrag[49]+f[20]*alphaDrag[48]); 
  out[37] += 0.3872983346207416*f[45]*alphaDrag[84]+(0.276641667586244*f[44]+0.4330127018922193*f[18])*alphaDrag[83]+(0.276641667586244*f[37]+0.4330127018922193*f[10])*alphaDrag[73]+0.3872983346207416*f[33]*alphaDrag[68]+(0.276641667586244*f[32]+0.4330127018922193*f[7])*alphaDrag[67]+0.3872983346207416*f[31]*alphaDrag[64]+0.276641667586244*f[21]*alphaDrag[59]+0.4330127018922193*(f[3]*alphaDrag[59]+f[44]*alphaDrag[57])+0.3872983346207416*(f[17]*alphaDrag[56]+f[15]*alphaDrag[53])+0.4330127018922193*(f[37]*alphaDrag[52]+f[32]*alphaDrag[50])+0.3872983346207416*f[6]*alphaDrag[49]+0.4330127018922193*f[21]*alphaDrag[48]+0.3872983346207416*alphaDrag[33]*f[45]+(0.276641667586244*alphaDrag[32]+0.4330127018922193*alphaDrag[7])*f[44]+(0.276641667586244*alphaDrag[21]+0.4330127018922193*alphaDrag[3])*f[37]+0.3872983346207416*alphaDrag[20]*f[36]+0.276641667586244*alphaDrag[19]*f[35]+0.4330127018922193*(alphaDrag[2]*f[35]+f[18]*alphaDrag[32])+0.3872983346207416*alphaDrag[15]*f[31]+0.276641667586244*alphaDrag[11]*f[25]+0.4330127018922193*(alphaDrag[0]*f[25]+f[10]*alphaDrag[21]+f[9]*alphaDrag[19])+0.3872983346207416*(alphaDrag[6]*f[17]+alphaDrag[5]*f[16])+0.4330127018922193*f[4]*alphaDrag[11]+0.3872983346207416*alphaDrag[1]*f[8]; 
  out[38] += (0.276641667586244*f[45]+0.4330127018922193*f[17])*alphaDrag[84]+0.3872983346207416*f[44]*alphaDrag[83]+(0.276641667586244*f[38]+0.4330127018922193*f[10])*alphaDrag[74]+(0.276641667586244*f[33]+0.4330127018922193*f[6])*alphaDrag[68]+0.3872983346207416*(f[32]*alphaDrag[67]+f[31]*alphaDrag[64])+(0.276641667586244*f[22]+0.4330127018922193*f[3])*alphaDrag[60]+0.3872983346207416*f[18]*alphaDrag[57]+0.4330127018922193*f[45]*alphaDrag[56]+0.3872983346207416*f[15]*alphaDrag[53]+0.4330127018922193*f[38]*alphaDrag[52]+0.3872983346207416*f[7]*alphaDrag[50]+0.4330127018922193*(f[33]*alphaDrag[49]+f[22]*alphaDrag[48])+(0.276641667586244*alphaDrag[33]+0.4330127018922193*alphaDrag[6])*f[45]+0.3872983346207416*alphaDrag[32]*f[44]+(0.276641667586244*alphaDrag[22]+0.4330127018922193*alphaDrag[3])*f[38]+(0.276641667586244*alphaDrag[20]+0.4330127018922193*alphaDrag[1])*f[36]+0.3872983346207416*alphaDrag[19]*f[35]+0.4330127018922193*f[17]*alphaDrag[33]+0.3872983346207416*alphaDrag[15]*f[31]+0.276641667586244*alphaDrag[12]*f[26]+0.4330127018922193*(alphaDrag[0]*f[26]+f[10]*alphaDrag[22]+f[8]*alphaDrag[20])+0.3872983346207416*(alphaDrag[7]*f[18]+alphaDrag[5]*f[16])+0.4330127018922193*f[4]*alphaDrag[12]+0.3872983346207416*alphaDrag[2]*f[9]; 
  out[39] += 0.3872983346207416*(f[46]*alphaDrag[83]+f[39]*alphaDrag[73]+f[34]*alphaDrag[67])+0.4330127018922193*f[40]*alphaDrag[64]+0.3872983346207416*f[23]*alphaDrag[59]+0.4330127018922193*(f[46]*alphaDrag[57]+f[27]*alphaDrag[56]+f[24]*alphaDrag[53]+f[39]*alphaDrag[52]+f[34]*alphaDrag[50]+f[13]*alphaDrag[49]+f[23]*alphaDrag[48])+(0.7745966692414833*alphaDrag[32]+0.8660254037844386*alphaDrag[7])*f[46]+0.9682458365518543*alphaDrag[12]*f[45]+0.8660254037844386*(alphaDrag[5]*f[44]+alphaDrag[15]*f[40])+(0.7745966692414833*alphaDrag[21]+0.8660254037844386*alphaDrag[3])*f[39]+0.9682458365518543*alphaDrag[20]*f[38]+0.8660254037844386*alphaDrag[1]*f[37]+0.9682458365518543*alphaDrag[22]*f[36]+0.8660254037844386*alphaDrag[15]*f[35]+0.9682458365518543*f[26]*alphaDrag[33]+0.8660254037844386*f[16]*alphaDrag[32]+(0.8660254037844386*alphaDrag[19]+0.9682458365518543*alphaDrag[2])*f[31]+0.8660254037844386*(alphaDrag[6]*(f[27]+f[25])+f[8]*alphaDrag[21])+0.9682458365518543*alphaDrag[5]*f[18]+0.8660254037844386*alphaDrag[11]*f[17]+0.9682458365518543*(alphaDrag[0]*f[17]+alphaDrag[7]*f[16]+f[9]*alphaDrag[15]+alphaDrag[1]*f[10]+alphaDrag[3]*f[8]+f[4]*alphaDrag[6]); 
  out[40] += 0.3872983346207416*(f[46]*alphaDrag[84]+f[40]*alphaDrag[74]+f[34]*alphaDrag[68])+0.4330127018922193*f[39]*alphaDrag[64]+0.3872983346207416*f[24]*alphaDrag[60]+0.4330127018922193*(f[27]*alphaDrag[57]+f[46]*alphaDrag[56]+f[23]*alphaDrag[53]+f[40]*alphaDrag[52]+f[13]*alphaDrag[50]+f[34]*alphaDrag[49]+f[24]*alphaDrag[48])+0.7745966692414833*alphaDrag[33]*f[46]+0.8660254037844386*(alphaDrag[6]*f[46]+alphaDrag[5]*f[45])+0.9682458365518543*alphaDrag[11]*f[44]+0.7745966692414833*alphaDrag[22]*f[40]+0.8660254037844386*(alphaDrag[3]*f[40]+alphaDrag[15]*f[39]+alphaDrag[2]*f[38])+0.9682458365518543*alphaDrag[19]*f[37]+0.8660254037844386*alphaDrag[15]*f[36]+0.9682458365518543*alphaDrag[21]*f[35]+0.8660254037844386*f[16]*alphaDrag[33]+0.9682458365518543*f[25]*alphaDrag[32]+(0.8660254037844386*alphaDrag[20]+0.9682458365518543*alphaDrag[1])*f[31]+0.8660254037844386*(alphaDrag[7]*(f[27]+f[26])+f[9]*alphaDrag[22]+alphaDrag[12]*f[18])+0.9682458365518543*(alphaDrag[0]*f[18]+alphaDrag[5]*f[17]+alphaDrag[6]*f[16]+f[8]*alphaDrag[15]+alphaDrag[2]*f[10]+alphaDrag[3]*f[9]+f[4]*alphaDrag[7]); 
  out[41] += (0.7745966692414833*(f[29]+f[19])+0.8660254037844386*f[2])*alphaDrag[84]+(0.7745966692414833*(f[28]+f[20])+0.8660254037844386*f[1])*alphaDrag[83]+(0.7745966692414833*f[41]+0.8660254037844386*f[5])*alphaDrag[74]+(0.7745966692414833*f[41]+0.8660254037844386*f[5])*alphaDrag[73]+(0.7745966692414833*f[35]+0.8660254037844386*f[9])*alphaDrag[68]+(0.7745966692414833*f[36]+0.8660254037844386*f[8])*alphaDrag[67]+(0.8660254037844386*(f[14]+f[12]+f[11])+0.9682458365518543*f[0])*alphaDrag[64]+0.8660254037844386*f[16]*(alphaDrag[60]+alphaDrag[59])+(0.8660254037844386*(f[28]+f[20])+0.9682458365518543*f[1])*alphaDrag[57]+(0.8660254037844386*(f[29]+f[19])+0.9682458365518543*f[2])*alphaDrag[56]+(0.8660254037844386*(f[26]+f[25])+0.9682458365518543*f[4])*alphaDrag[53]+(0.8660254037844386*f[41]+0.9682458365518543*f[5])*alphaDrag[52]+(0.8660254037844386*f[36]+0.9682458365518543*f[8])*alphaDrag[50]+0.8660254037844386*f[35]*alphaDrag[49]+0.9682458365518543*(f[9]*alphaDrag[49]+f[16]*alphaDrag[48]); 
  out[42] += 0.9682458365518543*f[22]*alphaDrag[84]+(0.7745966692414833*f[47]+0.8660254037844386*f[15])*alphaDrag[83]+0.9682458365518543*f[33]*alphaDrag[74]+(0.7745966692414833*f[42]+0.8660254037844386*f[6])*alphaDrag[73]+0.9682458365518543*f[38]*alphaDrag[68]+0.8660254037844386*(f[31]*alphaDrag[67]+(f[43]+f[32])*alphaDrag[64])+0.9682458365518543*(f[7]*alphaDrag[64]+f[45]*alphaDrag[60])+0.8660254037844386*f[17]*alphaDrag[59]+(0.8660254037844386*f[47]+0.9682458365518543*f[15])*alphaDrag[57]+(0.8660254037844386*(f[30]+f[21])+0.9682458365518543*f[3])*alphaDrag[56]+(0.8660254037844386*f[44]+0.9682458365518543*f[18])*alphaDrag[53]+0.8660254037844386*f[42]*alphaDrag[52]+0.9682458365518543*(f[6]*alphaDrag[52]+f[31]*alphaDrag[50])+0.8660254037844386*f[37]*alphaDrag[49]+0.9682458365518543*(f[10]*alphaDrag[49]+f[17]*alphaDrag[48])+0.3872983346207416*alphaDrag[32]*f[47]+0.4330127018922193*(alphaDrag[7]*f[47]+alphaDrag[15]*f[43])+(0.3872983346207416*alphaDrag[21]+0.4330127018922193*alphaDrag[3])*f[42]+0.3872983346207416*alphaDrag[19]*f[41]+0.4330127018922193*(alphaDrag[2]*f[41]+alphaDrag[6]*f[30]+alphaDrag[5]*f[29])+0.3872983346207416*alphaDrag[11]*f[28]+0.4330127018922193*(alphaDrag[0]*f[28]+alphaDrag[1]*f[14]); 
  out[43] += (0.7745966692414833*f[47]+0.8660254037844386*f[15])*alphaDrag[84]+0.9682458365518543*f[21]*alphaDrag[83]+(0.7745966692414833*f[43]+0.8660254037844386*f[7])*alphaDrag[74]+0.9682458365518543*f[32]*alphaDrag[73]+0.8660254037844386*f[31]*alphaDrag[68]+0.9682458365518543*f[37]*alphaDrag[67]+(0.8660254037844386*(f[42]+f[33])+0.9682458365518543*f[6])*alphaDrag[64]+0.8660254037844386*f[18]*alphaDrag[60]+0.9682458365518543*f[44]*alphaDrag[59]+(0.8660254037844386*(f[30]+f[22])+0.9682458365518543*f[3])*alphaDrag[57]+(0.8660254037844386*f[47]+0.9682458365518543*f[15])*alphaDrag[56]+(0.8660254037844386*f[45]+0.9682458365518543*f[17])*alphaDrag[53]+(0.8660254037844386*f[43]+0.9682458365518543*f[7])*alphaDrag[52]+0.8660254037844386*f[38]*alphaDrag[50]+0.9682458365518543*(f[10]*alphaDrag[50]+f[31]*alphaDrag[49]+f[18]*alphaDrag[48])+(0.3872983346207416*alphaDrag[33]+0.4330127018922193*alphaDrag[6])*f[47]+0.3872983346207416*alphaDrag[22]*f[43]+0.4330127018922193*(alphaDrag[3]*f[43]+alphaDrag[15]*f[42])+0.3872983346207416*alphaDrag[20]*f[41]+0.4330127018922193*(alphaDrag[1]*f[41]+alphaDrag[7]*f[30])+0.3872983346207416*alphaDrag[12]*f[29]+0.4330127018922193*(alphaDrag[0]*f[29]+alphaDrag[5]*f[28]+alphaDrag[2]*f[14]); 
  out[44] += 0.3464101615137755*f[31]*alphaDrag[84]+(0.3872983346207416*f[38]+0.276641667586244*f[37]+0.4330127018922193*f[10])*alphaDrag[83]+0.3872983346207416*f[44]*alphaDrag[74]+(0.276641667586244*f[44]+0.4330127018922193*f[18])*alphaDrag[73]+0.3464101615137755*f[15]*alphaDrag[68]+(0.3872983346207416*f[22]+0.276641667586244*f[21]+0.4330127018922193*f[3])*alphaDrag[67]+(0.3464101615137755*f[45]+0.3872983346207416*f[17])*alphaDrag[64]+f[32]*(0.3872983346207416*alphaDrag[60]+0.276641667586244*alphaDrag[59])+0.4330127018922193*(f[7]*alphaDrag[59]+f[37]*alphaDrag[57])+0.3872983346207416*f[31]*alphaDrag[56]+(0.3464101615137755*f[33]+0.3872983346207416*f[6])*alphaDrag[53]+0.4330127018922193*(f[44]*alphaDrag[52]+f[21]*alphaDrag[50])+0.3872983346207416*f[15]*alphaDrag[49]+0.4330127018922193*f[32]*alphaDrag[48]+0.3464101615137755*alphaDrag[15]*f[45]+(0.3872983346207416*alphaDrag[22]+0.276641667586244*alphaDrag[21]+0.4330127018922193*alphaDrag[3])*f[44]+0.3872983346207416*alphaDrag[32]*f[38]+(0.276641667586244*alphaDrag[32]+0.4330127018922193*alphaDrag[7])*f[37]+0.3464101615137755*alphaDrag[5]*f[36]+(0.3872983346207416*alphaDrag[12]+0.276641667586244*alphaDrag[11]+0.4330127018922193*alphaDrag[0])*f[35]+0.3464101615137755*f[31]*alphaDrag[33]+0.4330127018922193*f[10]*alphaDrag[32]+0.3872983346207416*alphaDrag[6]*f[31]+alphaDrag[19]*(0.3872983346207416*f[26]+0.276641667586244*f[25])+0.4330127018922193*(alphaDrag[2]*f[25]+f[18]*alphaDrag[21])+0.3464101615137755*f[16]*alphaDrag[20]+0.4330127018922193*f[4]*alphaDrag[19]+0.3872983346207416*(alphaDrag[15]*f[17]+alphaDrag[1]*f[16])+0.4330127018922193*f[9]*alphaDrag[11]+0.3872983346207416*alphaDrag[5]*f[8]; 
  out[45] += (0.276641667586244*f[38]+0.3872983346207416*f[37]+0.4330127018922193*f[10])*alphaDrag[84]+0.3464101615137755*f[31]*alphaDrag[83]+(0.276641667586244*f[45]+0.4330127018922193*f[17])*alphaDrag[74]+0.3872983346207416*f[45]*alphaDrag[73]+(0.276641667586244*f[22]+0.3872983346207416*f[21]+0.4330127018922193*f[3])*alphaDrag[68]+0.3464101615137755*f[15]*alphaDrag[67]+(0.3464101615137755*f[44]+0.3872983346207416*f[18])*alphaDrag[64]+(0.276641667586244*f[33]+0.4330127018922193*f[6])*alphaDrag[60]+0.3872983346207416*(f[33]*alphaDrag[59]+f[31]*alphaDrag[57])+0.4330127018922193*f[38]*alphaDrag[56]+(0.3464101615137755*f[32]+0.3872983346207416*f[7])*alphaDrag[53]+0.4330127018922193*f[45]*alphaDrag[52]+0.3872983346207416*f[15]*alphaDrag[50]+0.4330127018922193*(f[22]*alphaDrag[49]+f[33]*alphaDrag[48])+(0.276641667586244*alphaDrag[22]+0.3872983346207416*alphaDrag[21]+0.4330127018922193*alphaDrag[3])*f[45]+0.3464101615137755*alphaDrag[15]*f[44]+(0.276641667586244*alphaDrag[33]+0.4330127018922193*alphaDrag[6])*f[38]+0.3872983346207416*alphaDrag[33]*f[37]+(0.276641667586244*alphaDrag[12]+0.3872983346207416*alphaDrag[11]+0.4330127018922193*alphaDrag[0])*f[36]+0.3464101615137755*alphaDrag[5]*f[35]+0.4330127018922193*f[10]*alphaDrag[33]+f[31]*(0.3464101615137755*alphaDrag[32]+0.3872983346207416*alphaDrag[7])+(0.276641667586244*alphaDrag[20]+0.4330127018922193*alphaDrag[1])*f[26]+0.3872983346207416*alphaDrag[20]*f[25]+0.4330127018922193*(f[17]*alphaDrag[22]+f[4]*alphaDrag[20])+0.3464101615137755*f[16]*alphaDrag[19]+0.3872983346207416*(alphaDrag[15]*f[18]+alphaDrag[2]*f[16])+0.4330127018922193*f[8]*alphaDrag[12]+0.3872983346207416*alphaDrag[5]*f[9]; 
  out[46] += 0.3872983346207416*(f[40]*alphaDrag[84]+f[39]*alphaDrag[83]+f[46]*(alphaDrag[74]+alphaDrag[73])+f[24]*alphaDrag[68]+f[23]*alphaDrag[67])+0.4330127018922193*f[27]*alphaDrag[64]+0.3872983346207416*f[34]*(alphaDrag[60]+alphaDrag[59])+0.4330127018922193*(f[39]*alphaDrag[57]+f[40]*alphaDrag[56]+f[13]*alphaDrag[53]+f[46]*alphaDrag[52]+f[23]*alphaDrag[50]+f[24]*alphaDrag[49]+f[34]*alphaDrag[48])+(0.7745966692414833*(alphaDrag[22]+alphaDrag[21])+0.8660254037844386*alphaDrag[3])*f[46]+(0.7745966692414833*alphaDrag[19]+0.8660254037844386*alphaDrag[2])*f[45]+(0.7745966692414833*alphaDrag[20]+0.8660254037844386*alphaDrag[1])*f[44]+(0.7745966692414833*alphaDrag[33]+0.8660254037844386*alphaDrag[6])*f[40]+0.7745966692414833*alphaDrag[32]*f[39]+0.8660254037844386*(alphaDrag[7]*f[39]+alphaDrag[5]*(f[38]+f[37]))+(0.7745966692414833*alphaDrag[32]+0.8660254037844386*alphaDrag[7])*f[36]+0.7745966692414833*alphaDrag[33]*f[35]+0.8660254037844386*(alphaDrag[6]*f[35]+f[9]*alphaDrag[33]+f[8]*alphaDrag[32])+(0.8660254037844386*(alphaDrag[12]+alphaDrag[11])+0.9682458365518543*alphaDrag[0])*f[31]+0.8660254037844386*(alphaDrag[15]*(f[27]+f[26]+f[25])+f[16]*(alphaDrag[22]+alphaDrag[21])+f[18]*alphaDrag[20]+f[17]*alphaDrag[19])+0.9682458365518543*(alphaDrag[1]*f[18]+alphaDrag[2]*f[17]+alphaDrag[3]*f[16]+f[4]*alphaDrag[15]+alphaDrag[5]*f[10]+alphaDrag[6]*f[9]+alphaDrag[7]*f[8]); 
  out[47] += (0.7745966692414833*(f[43]+f[32])+0.8660254037844386*f[7])*alphaDrag[84]+(0.7745966692414833*(f[42]+f[33])+0.8660254037844386*f[6])*alphaDrag[83]+(0.7745966692414833*f[47]+0.8660254037844386*f[15])*alphaDrag[74]+(0.7745966692414833*f[47]+0.8660254037844386*f[15])*alphaDrag[73]+(0.7745966692414833*f[44]+0.8660254037844386*f[18])*alphaDrag[68]+(0.7745966692414833*f[45]+0.8660254037844386*f[17])*alphaDrag[67]+(0.8660254037844386*(f[30]+f[22]+f[21])+0.9682458365518543*f[3])*alphaDrag[64]+0.8660254037844386*f[31]*(alphaDrag[60]+alphaDrag[59])+(0.8660254037844386*(f[42]+f[33])+0.9682458365518543*f[6])*alphaDrag[57]+(0.8660254037844386*(f[43]+f[32])+0.9682458365518543*f[7])*alphaDrag[56]+(0.8660254037844386*(f[38]+f[37])+0.9682458365518543*f[10])*alphaDrag[53]+(0.8660254037844386*f[47]+0.9682458365518543*f[15])*alphaDrag[52]+(0.8660254037844386*f[45]+0.9682458365518543*f[17])*alphaDrag[50]+0.8660254037844386*f[44]*alphaDrag[49]+0.9682458365518543*(f[18]*alphaDrag[49]+f[31]*alphaDrag[48])+(0.3872983346207416*(alphaDrag[22]+alphaDrag[21])+0.4330127018922193*alphaDrag[3])*f[47]+(0.3872983346207416*alphaDrag[33]+0.4330127018922193*alphaDrag[6])*f[43]+(0.3872983346207416*alphaDrag[32]+0.4330127018922193*alphaDrag[7])*f[42]+0.3872983346207416*(alphaDrag[12]+alphaDrag[11])*f[41]+0.4330127018922193*(alphaDrag[0]*f[41]+alphaDrag[15]*f[30])+(0.3872983346207416*alphaDrag[20]+0.4330127018922193*alphaDrag[1])*f[29]+0.3872983346207416*alphaDrag[19]*f[28]+0.4330127018922193*(alphaDrag[2]*f[28]+alphaDrag[5]*f[14]); 

  return fabs(0.625*alphaDrag[0]-0.6987712429686843*(alphaDrag[12]+alphaDrag[11]))+fabs(0.625*alphaDrag[48]-0.6987712429686843*(alphaDrag[60]+alphaDrag[59])); 

} 
