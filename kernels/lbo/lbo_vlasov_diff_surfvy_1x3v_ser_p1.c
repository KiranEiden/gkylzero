#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH void lbo_vlasov_diff_surfvy_1x3v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[4]: cell-center coordinates. 
  // dxv[4]: cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[8]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fl/fc/fr: distribution function in cells 
  // out: incremented distribution function in cell 

  const double *nuVtSqSum = &nuPrimMomsSum[6];

  double rdvSq4 = 4.0/(dxv[2]*dxv[2]); 
  double temp_diff[40] = {0.0}; 
  double diff_incr[40] = {0.0}; 

  temp_diff[0] = 0.6708203932499369*fr[24]+0.6708203932499369*fl[24]-1.341640786499874*fc[24]-1.190784930203603*fr[3]+1.190784930203603*fl[3]+0.9375*fr[0]+0.9375*fl[0]-1.875*fc[0]; 
  temp_diff[1] = 0.6708203932499369*fr[25]+0.6708203932499369*fl[25]-1.341640786499874*fc[25]-1.190784930203603*fr[6]+1.190784930203603*fl[6]+0.9375*fr[1]+0.9375*fl[1]-1.875*fc[1]; 
  temp_diff[2] = 0.6708203932499369*fr[26]+0.6708203932499369*fl[26]-1.341640786499874*fc[26]-1.190784930203603*fr[7]+1.190784930203603*fl[7]+0.9375*fr[2]+0.9375*fl[2]-1.875*fc[2]; 
  temp_diff[3] = 0.7382874503707888*fr[24]-0.7382874503707888*fl[24]-1.453125*fr[3]-1.453125*fl[3]-5.34375*fc[3]+1.190784930203603*fr[0]-1.190784930203603*fl[0]; 
  temp_diff[4] = 0.6708203932499369*fr[27]+0.6708203932499369*fl[27]-1.341640786499874*fc[27]-1.190784930203603*fr[10]+1.190784930203603*fl[10]+0.9375*fr[4]+0.9375*fl[4]-1.875*fc[4]; 
  temp_diff[5] = 0.6708203932499369*fr[28]+0.6708203932499369*fl[28]-1.341640786499874*fc[28]-1.190784930203603*fr[11]+1.190784930203603*fl[11]+0.9375*fr[5]+0.9375*fl[5]-1.875*fc[5]; 
  temp_diff[6] = 0.7382874503707888*fr[25]-0.7382874503707888*fl[25]-1.453125*fr[6]-1.453125*fl[6]-5.34375*fc[6]+1.190784930203603*fr[1]-1.190784930203603*fl[1]; 
  temp_diff[7] = 0.7382874503707888*fr[26]-0.7382874503707888*fl[26]-1.453125*fr[7]-1.453125*fl[7]-5.34375*fc[7]+1.190784930203603*fr[2]-1.190784930203603*fl[2]; 
  temp_diff[8] = 0.6708203932499369*fr[29]+0.6708203932499369*fl[29]-1.341640786499874*fc[29]-1.190784930203603*fr[13]+1.190784930203603*fl[13]+0.9375*fr[8]+0.9375*fl[8]-1.875*fc[8]; 
  temp_diff[9] = 0.6708203932499369*fr[30]+0.6708203932499369*fl[30]-1.341640786499874*fc[30]-1.190784930203603*fr[14]+1.190784930203603*fl[14]+0.9375*fr[9]+0.9375*fl[9]-1.875*fc[9]; 
  temp_diff[10] = 0.7382874503707888*fr[27]-0.7382874503707888*fl[27]-1.453125*fr[10]-1.453125*fl[10]-5.34375*fc[10]+1.190784930203603*fr[4]-1.190784930203603*fl[4]; 
  temp_diff[11] = 0.7382874503707888*fr[28]-0.7382874503707888*fl[28]-1.453125*fr[11]-1.453125*fl[11]-5.34375*fc[11]+1.190784930203603*fr[5]-1.190784930203603*fl[5]; 
  temp_diff[12] = 0.6708203932499369*fr[31]+0.6708203932499369*fl[31]-1.341640786499874*fc[31]-1.190784930203603*fr[15]+1.190784930203603*fl[15]+0.9375*fr[12]+0.9375*fl[12]-1.875*fc[12]; 
  temp_diff[13] = 0.7382874503707888*fr[29]-0.7382874503707888*fl[29]-1.453125*fr[13]-1.453125*fl[13]-5.34375*fc[13]+1.190784930203603*fr[8]-1.190784930203603*fl[8]; 
  temp_diff[14] = 0.7382874503707888*fr[30]-0.7382874503707888*fl[30]-1.453125*fr[14]-1.453125*fl[14]-5.34375*fc[14]+1.190784930203603*fr[9]-1.190784930203603*fl[9]; 
  temp_diff[15] = 0.7382874503707888*fr[31]-0.7382874503707888*fl[31]-1.453125*fr[15]-1.453125*fl[15]-5.34375*fc[15]+1.190784930203603*fr[12]-1.190784930203603*fl[12]; 
  temp_diff[16] = (-1.190784930203603*fr[18])+1.190784930203603*fl[18]+0.9375*fr[16]+0.9375*fl[16]-1.875*fc[16]; 
  temp_diff[17] = (-1.190784930203603*fr[20])+1.190784930203603*fl[20]+0.9375*fr[17]+0.9375*fl[17]-1.875*fc[17]; 
  temp_diff[18] = (-1.453125*fr[18])-1.453125*fl[18]-5.34375*fc[18]+1.190784930203603*fr[16]-1.190784930203603*fl[16]; 
  temp_diff[19] = (-1.190784930203603*fr[22])+1.190784930203603*fl[22]+0.9375*fr[19]+0.9375*fl[19]-1.875*fc[19]; 
  temp_diff[20] = (-1.453125*fr[20])-1.453125*fl[20]-5.34375*fc[20]+1.190784930203603*fr[17]-1.190784930203603*fl[17]; 
  temp_diff[21] = (-1.190784930203603*fr[23])+1.190784930203603*fl[23]+0.9375*fr[21]+0.9375*fl[21]-1.875*fc[21]; 
  temp_diff[22] = (-1.453125*fr[22])-1.453125*fl[22]-5.34375*fc[22]+1.190784930203603*fr[19]-1.190784930203603*fl[19]; 
  temp_diff[23] = (-1.453125*fr[23])-1.453125*fl[23]-5.34375*fc[23]+1.190784930203603*fr[21]-1.190784930203603*fl[21]; 
  temp_diff[24] = (-0.140625*fr[24])-0.140625*fl[24]-6.28125*fc[24]-0.3025768239224545*fr[3]+0.3025768239224545*fl[3]+0.4192627457812106*fr[0]+0.4192627457812106*fl[0]-0.8385254915624212*fc[0]; 
  temp_diff[25] = (-0.140625*fr[25])-0.140625*fl[25]-6.28125*fc[25]-0.3025768239224544*fr[6]+0.3025768239224544*fl[6]+0.4192627457812105*fr[1]+0.4192627457812105*fl[1]-0.8385254915624211*fc[1]; 
  temp_diff[26] = (-0.140625*fr[26])-0.140625*fl[26]-6.28125*fc[26]-0.3025768239224544*fr[7]+0.3025768239224544*fl[7]+0.4192627457812105*fr[2]+0.4192627457812105*fl[2]-0.8385254915624211*fc[2]; 
  temp_diff[27] = (-0.140625*fr[27])-0.140625*fl[27]-6.28125*fc[27]-0.3025768239224544*fr[10]+0.3025768239224544*fl[10]+0.4192627457812105*fr[4]+0.4192627457812105*fl[4]-0.8385254915624211*fc[4]; 
  temp_diff[28] = (-0.140625*fr[28])-0.140625*fl[28]-6.28125*fc[28]-0.3025768239224545*fr[11]+0.3025768239224545*fl[11]+0.4192627457812106*fr[5]+0.4192627457812106*fl[5]-0.8385254915624212*fc[5]; 
  temp_diff[29] = (-0.140625*fr[29])-0.140625*fl[29]-6.28125*fc[29]-0.3025768239224545*fr[13]+0.3025768239224545*fl[13]+0.4192627457812106*fr[8]+0.4192627457812106*fl[8]-0.8385254915624212*fc[8]; 
  temp_diff[30] = (-0.140625*fr[30])-0.140625*fl[30]-6.28125*fc[30]-0.3025768239224545*fr[14]+0.3025768239224545*fl[14]+0.4192627457812106*fr[9]+0.4192627457812106*fl[9]-0.8385254915624212*fc[9]; 
  temp_diff[31] = (-0.140625*fr[31])-0.140625*fl[31]-6.28125*fc[31]-0.3025768239224544*fr[15]+0.3025768239224544*fl[15]+0.4192627457812105*fr[12]+0.4192627457812105*fl[12]-0.8385254915624211*fc[12]; 
  temp_diff[32] = (-1.190784930203603*fr[35])+1.190784930203603*fl[35]+0.9375*fr[32]+0.9375*fl[32]-1.875*fc[32]; 
  temp_diff[33] = (-1.190784930203603*fr[37])+1.190784930203603*fl[37]+0.9375*fr[33]+0.9375*fl[33]-1.875*fc[33]; 
  temp_diff[34] = (-1.190784930203603*fr[38])+1.190784930203603*fl[38]+0.9375*fr[34]+0.9375*fl[34]-1.875*fc[34]; 
  temp_diff[35] = (-1.453125*fr[35])-1.453125*fl[35]-5.34375*fc[35]+1.190784930203603*fr[32]-1.190784930203603*fl[32]; 
  temp_diff[36] = (-1.190784930203603*fr[39])+1.190784930203603*fl[39]+0.9375*fr[36]+0.9375*fl[36]-1.875*fc[36]; 
  temp_diff[37] = (-1.453125*fr[37])-1.453125*fl[37]-5.34375*fc[37]+1.190784930203603*fr[33]-1.190784930203603*fl[33]; 
  temp_diff[38] = (-1.453125*fr[38])-1.453125*fl[38]-5.34375*fc[38]+1.190784930203603*fr[34]-1.190784930203603*fl[34]; 
  temp_diff[39] = (-1.453125*fr[39])-1.453125*fl[39]-5.34375*fc[39]+1.190784930203603*fr[36]-1.190784930203603*fl[36]; 

  diff_incr[0] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[1]+0.7071067811865475*nuVtSqSum[0]*temp_diff[0]; 
  diff_incr[1] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[1]+0.7071067811865475*temp_diff[0]*nuVtSqSum[1]; 
  diff_incr[2] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[5]+0.7071067811865475*nuVtSqSum[0]*temp_diff[2]; 
  diff_incr[3] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[6]+0.7071067811865475*nuVtSqSum[0]*temp_diff[3]; 
  diff_incr[4] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[8]+0.7071067811865475*nuVtSqSum[0]*temp_diff[4]; 
  diff_incr[5] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[5]+0.7071067811865475*nuVtSqSum[1]*temp_diff[2]; 
  diff_incr[6] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[6]+0.7071067811865475*nuVtSqSum[1]*temp_diff[3]; 
  diff_incr[7] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[11]+0.7071067811865475*nuVtSqSum[0]*temp_diff[7]; 
  diff_incr[8] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[8]+0.7071067811865475*nuVtSqSum[1]*temp_diff[4]; 
  diff_incr[9] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[12]+0.7071067811865475*nuVtSqSum[0]*temp_diff[9]; 
  diff_incr[10] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[13]+0.7071067811865475*nuVtSqSum[0]*temp_diff[10]; 
  diff_incr[11] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[11]+0.7071067811865475*nuVtSqSum[1]*temp_diff[7]; 
  diff_incr[12] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[12]+0.7071067811865475*nuVtSqSum[1]*temp_diff[9]; 
  diff_incr[13] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[13]+0.7071067811865475*nuVtSqSum[1]*temp_diff[10]; 
  diff_incr[14] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[15]+0.7071067811865475*nuVtSqSum[0]*temp_diff[14]; 
  diff_incr[15] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[15]+0.7071067811865475*nuVtSqSum[1]*temp_diff[14]; 
  diff_incr[16] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[17]+0.7071067811865475*nuVtSqSum[0]*temp_diff[16]; 
  diff_incr[17] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[17]+0.7071067811865475*nuVtSqSum[1]*temp_diff[16]; 
  diff_incr[18] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[20]+0.7071067811865475*nuVtSqSum[0]*temp_diff[18]; 
  diff_incr[19] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[21]+0.7071067811865475*nuVtSqSum[0]*temp_diff[19]; 
  diff_incr[20] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[20]+0.7071067811865475*nuVtSqSum[1]*temp_diff[18]; 
  diff_incr[21] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[21]+0.7071067811865475*nuVtSqSum[1]*temp_diff[19]; 
  diff_incr[22] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[23]+0.7071067811865475*nuVtSqSum[0]*temp_diff[22]; 
  diff_incr[23] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[23]+0.7071067811865475*nuVtSqSum[1]*temp_diff[22]; 
  diff_incr[24] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[25]+0.7071067811865475*nuVtSqSum[0]*temp_diff[24]; 
  diff_incr[25] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[25]+0.7071067811865475*nuVtSqSum[1]*temp_diff[24]; 
  diff_incr[26] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[28]+0.7071067811865475*nuVtSqSum[0]*temp_diff[26]; 
  diff_incr[27] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[29]+0.7071067811865475*nuVtSqSum[0]*temp_diff[27]; 
  diff_incr[28] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[28]+0.7071067811865475*nuVtSqSum[1]*temp_diff[26]; 
  diff_incr[29] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[29]+0.7071067811865475*nuVtSqSum[1]*temp_diff[27]; 
  diff_incr[30] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[31]+0.7071067811865475*nuVtSqSum[0]*temp_diff[30]; 
  diff_incr[31] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[31]+0.7071067811865475*nuVtSqSum[1]*temp_diff[30]; 
  diff_incr[32] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[33]+0.7071067811865475*nuVtSqSum[0]*temp_diff[32]; 
  diff_incr[33] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[33]+0.7071067811865475*nuVtSqSum[1]*temp_diff[32]; 
  diff_incr[34] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[36]+0.7071067811865475*nuVtSqSum[0]*temp_diff[34]; 
  diff_incr[35] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[37]+0.7071067811865475*nuVtSqSum[0]*temp_diff[35]; 
  diff_incr[36] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[36]+0.7071067811865475*nuVtSqSum[1]*temp_diff[34]; 
  diff_incr[37] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[37]+0.7071067811865475*nuVtSqSum[1]*temp_diff[35]; 
  diff_incr[38] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[39]+0.7071067811865475*nuVtSqSum[0]*temp_diff[38]; 
  diff_incr[39] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[39]+0.7071067811865475*nuVtSqSum[1]*temp_diff[38]; 

  out[0] += diff_incr[0]*rdvSq4; 
  out[1] += diff_incr[1]*rdvSq4; 
  out[2] += diff_incr[2]*rdvSq4; 
  out[3] += diff_incr[3]*rdvSq4; 
  out[4] += diff_incr[4]*rdvSq4; 
  out[5] += diff_incr[5]*rdvSq4; 
  out[6] += diff_incr[6]*rdvSq4; 
  out[7] += diff_incr[7]*rdvSq4; 
  out[8] += diff_incr[8]*rdvSq4; 
  out[9] += diff_incr[9]*rdvSq4; 
  out[10] += diff_incr[10]*rdvSq4; 
  out[11] += diff_incr[11]*rdvSq4; 
  out[12] += diff_incr[12]*rdvSq4; 
  out[13] += diff_incr[13]*rdvSq4; 
  out[14] += diff_incr[14]*rdvSq4; 
  out[15] += diff_incr[15]*rdvSq4; 
  out[16] += diff_incr[16]*rdvSq4; 
  out[17] += diff_incr[17]*rdvSq4; 
  out[18] += diff_incr[18]*rdvSq4; 
  out[19] += diff_incr[19]*rdvSq4; 
  out[20] += diff_incr[20]*rdvSq4; 
  out[21] += diff_incr[21]*rdvSq4; 
  out[22] += diff_incr[22]*rdvSq4; 
  out[23] += diff_incr[23]*rdvSq4; 
  out[24] += diff_incr[24]*rdvSq4; 
  out[25] += diff_incr[25]*rdvSq4; 
  out[26] += diff_incr[26]*rdvSq4; 
  out[27] += diff_incr[27]*rdvSq4; 
  out[28] += diff_incr[28]*rdvSq4; 
  out[29] += diff_incr[29]*rdvSq4; 
  out[30] += diff_incr[30]*rdvSq4; 
  out[31] += diff_incr[31]*rdvSq4; 
  out[32] += diff_incr[32]*rdvSq4; 
  out[33] += diff_incr[33]*rdvSq4; 
  out[34] += diff_incr[34]*rdvSq4; 
  out[35] += diff_incr[35]*rdvSq4; 
  out[36] += diff_incr[36]*rdvSq4; 
  out[37] += diff_incr[37]*rdvSq4; 
  out[38] += diff_incr[38]*rdvSq4; 
  out[39] += diff_incr[39]*rdvSq4; 
} 
