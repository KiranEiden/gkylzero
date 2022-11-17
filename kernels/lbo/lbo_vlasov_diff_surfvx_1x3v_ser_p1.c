#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH void lbo_vlasov_diff_surfvx_1x3v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[4]: cell-center coordinates. 
  // dxv[4]: cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[8]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fl/fc/fr: distribution function in cells 
  // out: incremented distribution function in cell 

  const double *nuVtSqSum = &nuPrimMomsSum[6];

  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 
  double incr[40]; 

  double f_xx[40] = {0.0}; 
  f_xx[0] = 0.6708203932499369*fr[16]+0.6708203932499369*fl[16]-1.341640786499874*fc[16]-1.190784930203603*fr[2]+1.190784930203603*fl[2]+0.9375*fr[0]+0.9375*fl[0]-1.875*fc[0]; 
  f_xx[1] = 0.6708203932499369*fr[17]+0.6708203932499369*fl[17]-1.341640786499874*fc[17]-1.190784930203603*fr[5]+1.190784930203603*fl[5]+0.9375*fr[1]+0.9375*fl[1]-1.875*fc[1]; 
  f_xx[2] = 0.7382874503707888*fr[16]-0.7382874503707888*fl[16]-1.453125*fr[2]-1.453125*fl[2]-5.34375*fc[2]+1.190784930203603*fr[0]-1.190784930203603*fl[0]; 
  f_xx[3] = 0.6708203932499369*fr[18]+0.6708203932499369*fl[18]-1.341640786499874*fc[18]-1.190784930203603*fr[7]+1.190784930203603*fl[7]+0.9375*fr[3]+0.9375*fl[3]-1.875*fc[3]; 
  f_xx[4] = 0.6708203932499369*fr[19]+0.6708203932499369*fl[19]-1.341640786499874*fc[19]-1.190784930203603*fr[9]+1.190784930203603*fl[9]+0.9375*fr[4]+0.9375*fl[4]-1.875*fc[4]; 
  f_xx[5] = 0.7382874503707888*fr[17]-0.7382874503707888*fl[17]-1.453125*fr[5]-1.453125*fl[5]-5.34375*fc[5]+1.190784930203603*fr[1]-1.190784930203603*fl[1]; 
  f_xx[6] = 0.6708203932499369*fr[20]+0.6708203932499369*fl[20]-1.341640786499874*fc[20]-1.190784930203603*fr[11]+1.190784930203603*fl[11]+0.9375*fr[6]+0.9375*fl[6]-1.875*fc[6]; 
  f_xx[7] = 0.7382874503707888*fr[18]-0.7382874503707888*fl[18]-1.453125*fr[7]-1.453125*fl[7]-5.34375*fc[7]+1.190784930203603*fr[3]-1.190784930203603*fl[3]; 
  f_xx[8] = 0.6708203932499369*fr[21]+0.6708203932499369*fl[21]-1.341640786499874*fc[21]-1.190784930203603*fr[12]+1.190784930203603*fl[12]+0.9375*fr[8]+0.9375*fl[8]-1.875*fc[8]; 
  f_xx[9] = 0.7382874503707888*fr[19]-0.7382874503707888*fl[19]-1.453125*fr[9]-1.453125*fl[9]-5.34375*fc[9]+1.190784930203603*fr[4]-1.190784930203603*fl[4]; 
  f_xx[10] = 0.6708203932499369*fr[22]+0.6708203932499369*fl[22]-1.341640786499874*fc[22]-1.190784930203603*fr[14]+1.190784930203603*fl[14]+0.9375*fr[10]+0.9375*fl[10]-1.875*fc[10]; 
  f_xx[11] = 0.7382874503707888*fr[20]-0.7382874503707888*fl[20]-1.453125*fr[11]-1.453125*fl[11]-5.34375*fc[11]+1.190784930203603*fr[6]-1.190784930203603*fl[6]; 
  f_xx[12] = 0.7382874503707888*fr[21]-0.7382874503707888*fl[21]-1.453125*fr[12]-1.453125*fl[12]-5.34375*fc[12]+1.190784930203603*fr[8]-1.190784930203603*fl[8]; 
  f_xx[13] = 0.6708203932499369*fr[23]+0.6708203932499369*fl[23]-1.341640786499874*fc[23]-1.190784930203603*fr[15]+1.190784930203603*fl[15]+0.9375*fr[13]+0.9375*fl[13]-1.875*fc[13]; 
  f_xx[14] = 0.7382874503707888*fr[22]-0.7382874503707888*fl[22]-1.453125*fr[14]-1.453125*fl[14]-5.34375*fc[14]+1.190784930203603*fr[10]-1.190784930203603*fl[10]; 
  f_xx[15] = 0.7382874503707888*fr[23]-0.7382874503707888*fl[23]-1.453125*fr[15]-1.453125*fl[15]-5.34375*fc[15]+1.190784930203603*fr[13]-1.190784930203603*fl[13]; 
  f_xx[16] = (-0.140625*fr[16])-0.140625*fl[16]-6.28125*fc[16]-0.3025768239224545*fr[2]+0.3025768239224545*fl[2]+0.4192627457812106*fr[0]+0.4192627457812106*fl[0]-0.8385254915624212*fc[0]; 
  f_xx[17] = (-0.140625*fr[17])-0.140625*fl[17]-6.28125*fc[17]-0.3025768239224544*fr[5]+0.3025768239224544*fl[5]+0.4192627457812105*fr[1]+0.4192627457812105*fl[1]-0.8385254915624211*fc[1]; 
  f_xx[18] = (-0.140625*fr[18])-0.140625*fl[18]-6.28125*fc[18]-0.3025768239224544*fr[7]+0.3025768239224544*fl[7]+0.4192627457812105*fr[3]+0.4192627457812105*fl[3]-0.8385254915624211*fc[3]; 
  f_xx[19] = (-0.140625*fr[19])-0.140625*fl[19]-6.28125*fc[19]-0.3025768239224544*fr[9]+0.3025768239224544*fl[9]+0.4192627457812105*fr[4]+0.4192627457812105*fl[4]-0.8385254915624211*fc[4]; 
  f_xx[20] = (-0.140625*fr[20])-0.140625*fl[20]-6.28125*fc[20]-0.3025768239224545*fr[11]+0.3025768239224545*fl[11]+0.4192627457812106*fr[6]+0.4192627457812106*fl[6]-0.8385254915624212*fc[6]; 
  f_xx[21] = (-0.140625*fr[21])-0.140625*fl[21]-6.28125*fc[21]-0.3025768239224545*fr[12]+0.3025768239224545*fl[12]+0.4192627457812106*fr[8]+0.4192627457812106*fl[8]-0.8385254915624212*fc[8]; 
  f_xx[22] = (-0.140625*fr[22])-0.140625*fl[22]-6.28125*fc[22]-0.3025768239224545*fr[14]+0.3025768239224545*fl[14]+0.4192627457812106*fr[10]+0.4192627457812106*fl[10]-0.8385254915624212*fc[10]; 
  f_xx[23] = (-0.140625*fr[23])-0.140625*fl[23]-6.28125*fc[23]-0.3025768239224544*fr[15]+0.3025768239224544*fl[15]+0.4192627457812105*fr[13]+0.4192627457812105*fl[13]-0.8385254915624211*fc[13]; 
  f_xx[24] = (-1.190784930203603*fr[26])+1.190784930203603*fl[26]+0.9375*fr[24]+0.9375*fl[24]-1.875*fc[24]; 
  f_xx[25] = (-1.190784930203603*fr[28])+1.190784930203603*fl[28]+0.9375*fr[25]+0.9375*fl[25]-1.875*fc[25]; 
  f_xx[26] = (-1.453125*fr[26])-1.453125*fl[26]-5.34375*fc[26]+1.190784930203603*fr[24]-1.190784930203603*fl[24]; 
  f_xx[27] = (-1.190784930203603*fr[30])+1.190784930203603*fl[30]+0.9375*fr[27]+0.9375*fl[27]-1.875*fc[27]; 
  f_xx[28] = (-1.453125*fr[28])-1.453125*fl[28]-5.34375*fc[28]+1.190784930203603*fr[25]-1.190784930203603*fl[25]; 
  f_xx[29] = (-1.190784930203603*fr[31])+1.190784930203603*fl[31]+0.9375*fr[29]+0.9375*fl[29]-1.875*fc[29]; 
  f_xx[30] = (-1.453125*fr[30])-1.453125*fl[30]-5.34375*fc[30]+1.190784930203603*fr[27]-1.190784930203603*fl[27]; 
  f_xx[31] = (-1.453125*fr[31])-1.453125*fl[31]-5.34375*fc[31]+1.190784930203603*fr[29]-1.190784930203603*fl[29]; 
  f_xx[32] = (-1.190784930203603*fr[34])+1.190784930203603*fl[34]+0.9375*fr[32]+0.9375*fl[32]-1.875*fc[32]; 
  f_xx[33] = (-1.190784930203603*fr[36])+1.190784930203603*fl[36]+0.9375*fr[33]+0.9375*fl[33]-1.875*fc[33]; 
  f_xx[34] = (-1.453125*fr[34])-1.453125*fl[34]-5.34375*fc[34]+1.190784930203603*fr[32]-1.190784930203603*fl[32]; 
  f_xx[35] = (-1.190784930203603*fr[38])+1.190784930203603*fl[38]+0.9375*fr[35]+0.9375*fl[35]-1.875*fc[35]; 
  f_xx[36] = (-1.453125*fr[36])-1.453125*fl[36]-5.34375*fc[36]+1.190784930203603*fr[33]-1.190784930203603*fl[33]; 
  f_xx[37] = (-1.190784930203603*fr[39])+1.190784930203603*fl[39]+0.9375*fr[37]+0.9375*fl[37]-1.875*fc[37]; 
  f_xx[38] = (-1.453125*fr[38])-1.453125*fl[38]-5.34375*fc[38]+1.190784930203603*fr[35]-1.190784930203603*fl[35]; 
  f_xx[39] = (-1.453125*fr[39])-1.453125*fl[39]-5.34375*fc[39]+1.190784930203603*fr[37]-1.190784930203603*fl[37]; 

  incr[0] = 0.7071067811865475*f_xx[1]*nuVtSqSum[1]+0.7071067811865475*f_xx[0]*nuVtSqSum[0]; 
  incr[1] = 0.7071067811865475*f_xx[0]*nuVtSqSum[1]+0.7071067811865475*nuVtSqSum[0]*f_xx[1]; 
  incr[2] = 0.7071067811865475*nuVtSqSum[1]*f_xx[5]+0.7071067811865475*nuVtSqSum[0]*f_xx[2]; 
  incr[3] = 0.7071067811865475*nuVtSqSum[1]*f_xx[6]+0.7071067811865475*nuVtSqSum[0]*f_xx[3]; 
  incr[4] = 0.7071067811865475*nuVtSqSum[1]*f_xx[8]+0.7071067811865475*nuVtSqSum[0]*f_xx[4]; 
  incr[5] = 0.7071067811865475*nuVtSqSum[0]*f_xx[5]+0.7071067811865475*nuVtSqSum[1]*f_xx[2]; 
  incr[6] = 0.7071067811865475*nuVtSqSum[0]*f_xx[6]+0.7071067811865475*nuVtSqSum[1]*f_xx[3]; 
  incr[7] = 0.7071067811865475*nuVtSqSum[1]*f_xx[11]+0.7071067811865475*nuVtSqSum[0]*f_xx[7]; 
  incr[8] = 0.7071067811865475*nuVtSqSum[0]*f_xx[8]+0.7071067811865475*nuVtSqSum[1]*f_xx[4]; 
  incr[9] = 0.7071067811865475*nuVtSqSum[1]*f_xx[12]+0.7071067811865475*nuVtSqSum[0]*f_xx[9]; 
  incr[10] = 0.7071067811865475*nuVtSqSum[1]*f_xx[13]+0.7071067811865475*nuVtSqSum[0]*f_xx[10]; 
  incr[11] = 0.7071067811865475*nuVtSqSum[0]*f_xx[11]+0.7071067811865475*nuVtSqSum[1]*f_xx[7]; 
  incr[12] = 0.7071067811865475*nuVtSqSum[0]*f_xx[12]+0.7071067811865475*nuVtSqSum[1]*f_xx[9]; 
  incr[13] = 0.7071067811865475*nuVtSqSum[0]*f_xx[13]+0.7071067811865475*nuVtSqSum[1]*f_xx[10]; 
  incr[14] = 0.7071067811865475*nuVtSqSum[1]*f_xx[15]+0.7071067811865475*nuVtSqSum[0]*f_xx[14]; 
  incr[15] = 0.7071067811865475*nuVtSqSum[0]*f_xx[15]+0.7071067811865475*nuVtSqSum[1]*f_xx[14]; 
  incr[16] = 0.7071067811865475*nuVtSqSum[1]*f_xx[17]+0.7071067811865475*nuVtSqSum[0]*f_xx[16]; 
  incr[17] = 0.7071067811865475*nuVtSqSum[0]*f_xx[17]+0.7071067811865475*nuVtSqSum[1]*f_xx[16]; 
  incr[18] = 0.7071067811865475*nuVtSqSum[1]*f_xx[20]+0.7071067811865475*nuVtSqSum[0]*f_xx[18]; 
  incr[19] = 0.7071067811865475*nuVtSqSum[1]*f_xx[21]+0.7071067811865475*nuVtSqSum[0]*f_xx[19]; 
  incr[20] = 0.7071067811865475*nuVtSqSum[0]*f_xx[20]+0.7071067811865475*nuVtSqSum[1]*f_xx[18]; 
  incr[21] = 0.7071067811865475*nuVtSqSum[0]*f_xx[21]+0.7071067811865475*nuVtSqSum[1]*f_xx[19]; 
  incr[22] = 0.7071067811865475*nuVtSqSum[1]*f_xx[23]+0.7071067811865475*nuVtSqSum[0]*f_xx[22]; 
  incr[23] = 0.7071067811865475*nuVtSqSum[0]*f_xx[23]+0.7071067811865475*nuVtSqSum[1]*f_xx[22]; 
  incr[24] = 0.7071067811865475*nuVtSqSum[1]*f_xx[25]+0.7071067811865475*nuVtSqSum[0]*f_xx[24]; 
  incr[25] = 0.7071067811865475*nuVtSqSum[0]*f_xx[25]+0.7071067811865475*nuVtSqSum[1]*f_xx[24]; 
  incr[26] = 0.7071067811865475*nuVtSqSum[1]*f_xx[28]+0.7071067811865475*nuVtSqSum[0]*f_xx[26]; 
  incr[27] = 0.7071067811865475*nuVtSqSum[1]*f_xx[29]+0.7071067811865475*nuVtSqSum[0]*f_xx[27]; 
  incr[28] = 0.7071067811865475*nuVtSqSum[0]*f_xx[28]+0.7071067811865475*nuVtSqSum[1]*f_xx[26]; 
  incr[29] = 0.7071067811865475*nuVtSqSum[0]*f_xx[29]+0.7071067811865475*nuVtSqSum[1]*f_xx[27]; 
  incr[30] = 0.7071067811865475*nuVtSqSum[1]*f_xx[31]+0.7071067811865475*nuVtSqSum[0]*f_xx[30]; 
  incr[31] = 0.7071067811865475*nuVtSqSum[0]*f_xx[31]+0.7071067811865475*nuVtSqSum[1]*f_xx[30]; 
  incr[32] = 0.7071067811865475*nuVtSqSum[1]*f_xx[33]+0.7071067811865475*nuVtSqSum[0]*f_xx[32]; 
  incr[33] = 0.7071067811865475*nuVtSqSum[0]*f_xx[33]+0.7071067811865475*nuVtSqSum[1]*f_xx[32]; 
  incr[34] = 0.7071067811865475*nuVtSqSum[1]*f_xx[36]+0.7071067811865475*nuVtSqSum[0]*f_xx[34]; 
  incr[35] = 0.7071067811865475*nuVtSqSum[1]*f_xx[37]+0.7071067811865475*nuVtSqSum[0]*f_xx[35]; 
  incr[36] = 0.7071067811865475*nuVtSqSum[0]*f_xx[36]+0.7071067811865475*nuVtSqSum[1]*f_xx[34]; 
  incr[37] = 0.7071067811865475*nuVtSqSum[0]*f_xx[37]+0.7071067811865475*nuVtSqSum[1]*f_xx[35]; 
  incr[38] = 0.7071067811865475*nuVtSqSum[1]*f_xx[39]+0.7071067811865475*nuVtSqSum[0]*f_xx[38]; 
  incr[39] = 0.7071067811865475*nuVtSqSum[0]*f_xx[39]+0.7071067811865475*nuVtSqSum[1]*f_xx[38]; 

  out[0] += incr[0]*rdvSq4; 
  out[1] += incr[1]*rdvSq4; 
  out[2] += incr[2]*rdvSq4; 
  out[3] += incr[3]*rdvSq4; 
  out[4] += incr[4]*rdvSq4; 
  out[5] += incr[5]*rdvSq4; 
  out[6] += incr[6]*rdvSq4; 
  out[7] += incr[7]*rdvSq4; 
  out[8] += incr[8]*rdvSq4; 
  out[9] += incr[9]*rdvSq4; 
  out[10] += incr[10]*rdvSq4; 
  out[11] += incr[11]*rdvSq4; 
  out[12] += incr[12]*rdvSq4; 
  out[13] += incr[13]*rdvSq4; 
  out[14] += incr[14]*rdvSq4; 
  out[15] += incr[15]*rdvSq4; 
  out[16] += incr[16]*rdvSq4; 
  out[17] += incr[17]*rdvSq4; 
  out[18] += incr[18]*rdvSq4; 
  out[19] += incr[19]*rdvSq4; 
  out[20] += incr[20]*rdvSq4; 
  out[21] += incr[21]*rdvSq4; 
  out[22] += incr[22]*rdvSq4; 
  out[23] += incr[23]*rdvSq4; 
  out[24] += incr[24]*rdvSq4; 
  out[25] += incr[25]*rdvSq4; 
  out[26] += incr[26]*rdvSq4; 
  out[27] += incr[27]*rdvSq4; 
  out[28] += incr[28]*rdvSq4; 
  out[29] += incr[29]*rdvSq4; 
  out[30] += incr[30]*rdvSq4; 
  out[31] += incr[31]*rdvSq4; 
  out[32] += incr[32]*rdvSq4; 
  out[33] += incr[33]*rdvSq4; 
  out[34] += incr[34]*rdvSq4; 
  out[35] += incr[35]*rdvSq4; 
  out[36] += incr[36]*rdvSq4; 
  out[37] += incr[37]*rdvSq4; 
  out[38] += incr[38]*rdvSq4; 
  out[39] += incr[39]*rdvSq4; 
} 
