#include <gkyl_vlasov_lbo_kernels.h> 
GKYL_CU_DH void vlasov_lbo_diff_surfvz_1x3v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[4]:         cell-center coordinates. 
  // dxv[4]:       cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[9]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      distribution function in cells 
  // out:           incremented distribution function in cell 
  double rdvSq4 = 4.0/(dxv[3]*dxv[3]); 
  double temp_diff[48] = {0.0}; 
  double diff_incr[48] = {0.0}; 

  temp_diff[0] = 0.6708203932499369*fr[14]+0.6708203932499369*fl[14]-1.341640786499874*fc[14]-1.190784930203603*fr[4]+1.190784930203603*fl[4]+0.9375*fr[0]+0.9375*fl[0]-1.875*fc[0]; 
  temp_diff[1] = 0.6708203932499369*fr[28]+0.6708203932499369*fl[28]-1.341640786499874*fc[28]-1.190784930203603*fr[8]+1.190784930203603*fl[8]+0.9375*fr[1]+0.9375*fl[1]-1.875*fc[1]; 
  temp_diff[2] = 0.6708203932499369*fr[29]+0.6708203932499369*fl[29]-1.341640786499874*fc[29]-1.190784930203603*fr[9]+1.190784930203603*fl[9]+0.9375*fr[2]+0.9375*fl[2]-1.875*fc[2]; 
  temp_diff[3] = 0.6708203932499369*fr[30]+0.6708203932499369*fl[30]-1.341640786499874*fc[30]-1.190784930203603*fr[10]+1.190784930203603*fl[10]+0.9375*fr[3]+0.9375*fl[3]-1.875*fc[3]; 
  temp_diff[4] = 0.7382874503707888*fr[14]-0.7382874503707888*fl[14]-1.453125*fr[4]-1.453125*fl[4]-5.34375*fc[4]+1.190784930203603*fr[0]-1.190784930203603*fl[0]; 
  temp_diff[5] = 0.6708203932499369*fr[41]+0.6708203932499369*fl[41]-1.341640786499874*fc[41]-1.190784930203603*fr[16]+1.190784930203603*fl[16]+0.9375*fr[5]+0.9375*fl[5]-1.875*fc[5]; 
  temp_diff[6] = 0.6708203932499369*fr[42]+0.6708203932499369*fl[42]-1.341640786499874*fc[42]-1.190784930203603*fr[17]+1.190784930203603*fl[17]+0.9375*fr[6]+0.9375*fl[6]-1.875*fc[6]; 
  temp_diff[7] = 0.6708203932499369*fr[43]+0.6708203932499369*fl[43]-1.341640786499874*fc[43]-1.190784930203603*fr[18]+1.190784930203603*fl[18]+0.9375*fr[7]+0.9375*fl[7]-1.875*fc[7]; 
  temp_diff[8] = 0.7382874503707888*fr[28]-0.7382874503707888*fl[28]-1.453125*fr[8]-1.453125*fl[8]-5.34375*fc[8]+1.190784930203603*fr[1]-1.190784930203603*fl[1]; 
  temp_diff[9] = 0.7382874503707888*fr[29]-0.7382874503707888*fl[29]-1.453125*fr[9]-1.453125*fl[9]-5.34375*fc[9]+1.190784930203603*fr[2]-1.190784930203603*fl[2]; 
  temp_diff[10] = 0.7382874503707888*fr[30]-0.7382874503707888*fl[30]-1.453125*fr[10]-1.453125*fl[10]-5.34375*fc[10]+1.190784930203603*fr[3]-1.190784930203603*fl[3]; 
  temp_diff[11] = (-1.190784930203603*fr[25])+1.190784930203603*fl[25]+0.9375*fr[11]+0.9375*fl[11]-1.875*fc[11]; 
  temp_diff[12] = (-1.190784930203603*fr[26])+1.190784930203603*fl[26]+0.9375*fr[12]+0.9375*fl[12]-1.875*fc[12]; 
  temp_diff[13] = (-1.190784930203603*fr[27])+1.190784930203603*fl[27]+0.9375*fr[13]+0.9375*fl[13]-1.875*fc[13]; 
  temp_diff[14] = (-0.140625*fr[14])-0.140625*fl[14]-6.28125*fc[14]-0.3025768239224545*fr[4]+0.3025768239224545*fl[4]+0.4192627457812106*fr[0]+0.4192627457812106*fl[0]-0.8385254915624212*fc[0]; 
  temp_diff[15] = 0.6708203932499369*fr[47]+0.6708203932499369*fl[47]-1.341640786499874*fc[47]-1.190784930203603*fr[31]+1.190784930203603*fl[31]+0.9375*fr[15]+0.9375*fl[15]-1.875*fc[15]; 
  temp_diff[16] = 0.7382874503707888*fr[41]-0.7382874503707888*fl[41]-1.453125*fr[16]-1.453125*fl[16]-5.34375*fc[16]+1.190784930203603*fr[5]-1.190784930203603*fl[5]; 
  temp_diff[17] = 0.7382874503707888*fr[42]-0.7382874503707888*fl[42]-1.453125*fr[17]-1.453125*fl[17]-5.34375*fc[17]+1.190784930203603*fr[6]-1.190784930203603*fl[6]; 
  temp_diff[18] = 0.7382874503707888*fr[43]-0.7382874503707888*fl[43]-1.453125*fr[18]-1.453125*fl[18]-5.34375*fc[18]+1.190784930203603*fr[7]-1.190784930203603*fl[7]; 
  temp_diff[19] = (-1.190784930203603*fr[35])+1.190784930203603*fl[35]+0.9375*fr[19]+0.9375*fl[19]-1.875*fc[19]; 
  temp_diff[20] = (-1.190784930203603*fr[36])+1.190784930203603*fl[36]+0.9375*fr[20]+0.9375*fl[20]-1.875*fc[20]; 
  temp_diff[21] = (-1.190784930203603*fr[37])+1.190784930203603*fl[37]+0.9375*fr[21]+0.9375*fl[21]-1.875*fc[21]; 
  temp_diff[22] = (-1.190784930203603*fr[38])+1.190784930203603*fl[38]+0.9375*fr[22]+0.9375*fl[22]-1.875*fc[22]; 
  temp_diff[23] = (-1.190784930203603*fr[39])+1.190784930203603*fl[39]+0.9375*fr[23]+0.9375*fl[23]-1.875*fc[23]; 
  temp_diff[24] = (-1.190784930203603*fr[40])+1.190784930203603*fl[40]+0.9375*fr[24]+0.9375*fl[24]-1.875*fc[24]; 
  temp_diff[25] = (-1.453125*fr[25])-1.453125*fl[25]-5.34375*fc[25]+1.190784930203603*fr[11]-1.190784930203603*fl[11]; 
  temp_diff[26] = (-1.453125*fr[26])-1.453125*fl[26]-5.34375*fc[26]+1.190784930203603*fr[12]-1.190784930203603*fl[12]; 
  temp_diff[27] = (-1.453125*fr[27])-1.453125*fl[27]-5.34375*fc[27]+1.190784930203603*fr[13]-1.190784930203603*fl[13]; 
  temp_diff[28] = (-0.140625*fr[28])-0.140625*fl[28]-6.28125*fc[28]-0.3025768239224544*fr[8]+0.3025768239224544*fl[8]+0.4192627457812105*fr[1]+0.4192627457812105*fl[1]-0.8385254915624211*fc[1]; 
  temp_diff[29] = (-0.140625*fr[29])-0.140625*fl[29]-6.28125*fc[29]-0.3025768239224544*fr[9]+0.3025768239224544*fl[9]+0.4192627457812105*fr[2]+0.4192627457812105*fl[2]-0.8385254915624211*fc[2]; 
  temp_diff[30] = (-0.140625*fr[30])-0.140625*fl[30]-6.28125*fc[30]-0.3025768239224544*fr[10]+0.3025768239224544*fl[10]+0.4192627457812105*fr[3]+0.4192627457812105*fl[3]-0.8385254915624211*fc[3]; 
  temp_diff[31] = 0.7382874503707888*fr[47]-0.7382874503707888*fl[47]-1.453125*fr[31]-1.453125*fl[31]-5.34375*fc[31]+1.190784930203603*fr[15]-1.190784930203603*fl[15]; 
  temp_diff[32] = (-1.190784930203603*fr[44])+1.190784930203603*fl[44]+0.9375*fr[32]+0.9375*fl[32]-1.875*fc[32]; 
  temp_diff[33] = (-1.190784930203603*fr[45])+1.190784930203603*fl[45]+0.9375*fr[33]+0.9375*fl[33]-1.875*fc[33]; 
  temp_diff[34] = (-1.190784930203603*fr[46])+1.190784930203603*fl[46]+0.9375*fr[34]+0.9375*fl[34]-1.875*fc[34]; 
  temp_diff[35] = (-1.453125*fr[35])-1.453125*fl[35]-5.34375*fc[35]+1.190784930203603*fr[19]-1.190784930203603*fl[19]; 
  temp_diff[36] = (-1.453125*fr[36])-1.453125*fl[36]-5.34375*fc[36]+1.190784930203603*fr[20]-1.190784930203603*fl[20]; 
  temp_diff[37] = (-1.453125*fr[37])-1.453125*fl[37]-5.34375*fc[37]+1.190784930203603*fr[21]-1.190784930203603*fl[21]; 
  temp_diff[38] = (-1.453125*fr[38])-1.453125*fl[38]-5.34375*fc[38]+1.190784930203603*fr[22]-1.190784930203603*fl[22]; 
  temp_diff[39] = (-1.453125*fr[39])-1.453125*fl[39]-5.34375*fc[39]+1.190784930203603*fr[23]-1.190784930203603*fl[23]; 
  temp_diff[40] = (-1.453125*fr[40])-1.453125*fl[40]-5.34375*fc[40]+1.190784930203603*fr[24]-1.190784930203603*fl[24]; 
  temp_diff[41] = (-0.140625*fr[41])-0.140625*fl[41]-6.28125*fc[41]-0.3025768239224545*fr[16]+0.3025768239224545*fl[16]+0.4192627457812106*fr[5]+0.4192627457812106*fl[5]-0.8385254915624212*fc[5]; 
  temp_diff[42] = (-0.140625*fr[42])-0.140625*fl[42]-6.28125*fc[42]-0.3025768239224545*fr[17]+0.3025768239224545*fl[17]+0.4192627457812106*fr[6]+0.4192627457812106*fl[6]-0.8385254915624212*fc[6]; 
  temp_diff[43] = (-0.140625*fr[43])-0.140625*fl[43]-6.28125*fc[43]-0.3025768239224545*fr[18]+0.3025768239224545*fl[18]+0.4192627457812106*fr[7]+0.4192627457812106*fl[7]-0.8385254915624212*fc[7]; 
  temp_diff[44] = (-1.453125*fr[44])-1.453125*fl[44]-5.34375*fc[44]+1.190784930203603*fr[32]-1.190784930203603*fl[32]; 
  temp_diff[45] = (-1.453125*fr[45])-1.453125*fl[45]-5.34375*fc[45]+1.190784930203603*fr[33]-1.190784930203603*fl[33]; 
  temp_diff[46] = (-1.453125*fr[46])-1.453125*fl[46]-5.34375*fc[46]+1.190784930203603*fr[34]-1.190784930203603*fl[34]; 
  temp_diff[47] = (-0.140625*fr[47])-0.140625*fl[47]-6.28125*fc[47]-0.3025768239224544*fr[31]+0.3025768239224544*fl[31]+0.4192627457812105*fr[15]+0.4192627457812105*fl[15]-0.8385254915624211*fc[15]; 

  diff_incr[0] = 0.7071067811865475*nuVtSqSum[2]*temp_diff[11]+0.7071067811865475*nuVtSqSum[1]*temp_diff[1]+0.7071067811865475*nuVtSqSum[0]*temp_diff[0]; 
  diff_incr[1] = 0.6324555320336759*nuVtSqSum[1]*temp_diff[11]+0.6324555320336759*temp_diff[1]*nuVtSqSum[2]+0.7071067811865475*nuVtSqSum[0]*temp_diff[1]+0.7071067811865475*temp_diff[0]*nuVtSqSum[1]; 
  diff_incr[2] = 0.7071067811865475*nuVtSqSum[2]*temp_diff[19]+0.7071067811865475*nuVtSqSum[1]*temp_diff[5]+0.7071067811865475*nuVtSqSum[0]*temp_diff[2]; 
  diff_incr[3] = 0.7071067811865475*nuVtSqSum[2]*temp_diff[21]+0.7071067811865475*nuVtSqSum[1]*temp_diff[6]+0.7071067811865475*nuVtSqSum[0]*temp_diff[3]; 
  diff_incr[4] = 0.7071067811865475*nuVtSqSum[2]*temp_diff[25]+0.7071067811865475*nuVtSqSum[1]*temp_diff[8]+0.7071067811865475*nuVtSqSum[0]*temp_diff[4]; 
  diff_incr[5] = 0.632455532033676*nuVtSqSum[1]*temp_diff[19]+0.6324555320336759*nuVtSqSum[2]*temp_diff[5]+0.7071067811865475*nuVtSqSum[0]*temp_diff[5]+0.7071067811865475*nuVtSqSum[1]*temp_diff[2]; 
  diff_incr[6] = 0.632455532033676*nuVtSqSum[1]*temp_diff[21]+0.6324555320336759*nuVtSqSum[2]*temp_diff[6]+0.7071067811865475*nuVtSqSum[0]*temp_diff[6]+0.7071067811865475*nuVtSqSum[1]*temp_diff[3]; 
  diff_incr[7] = 0.7071067811865475*nuVtSqSum[2]*temp_diff[32]+0.7071067811865475*nuVtSqSum[1]*temp_diff[15]+0.7071067811865475*nuVtSqSum[0]*temp_diff[7]; 
  diff_incr[8] = 0.632455532033676*nuVtSqSum[1]*temp_diff[25]+0.6324555320336759*nuVtSqSum[2]*temp_diff[8]+0.7071067811865475*nuVtSqSum[0]*temp_diff[8]+0.7071067811865475*nuVtSqSum[1]*temp_diff[4]; 
  diff_incr[9] = 0.7071067811865475*nuVtSqSum[2]*temp_diff[35]+0.7071067811865475*nuVtSqSum[1]*temp_diff[16]+0.7071067811865475*nuVtSqSum[0]*temp_diff[9]; 
  diff_incr[10] = 0.7071067811865475*nuVtSqSum[2]*temp_diff[37]+0.7071067811865475*nuVtSqSum[1]*temp_diff[17]+0.7071067811865475*nuVtSqSum[0]*temp_diff[10]; 
  diff_incr[11] = 0.4517539514526256*nuVtSqSum[2]*temp_diff[11]+0.7071067811865475*nuVtSqSum[0]*temp_diff[11]+0.7071067811865475*temp_diff[0]*nuVtSqSum[2]+0.6324555320336759*nuVtSqSum[1]*temp_diff[1]; 
  diff_incr[12] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[20]+0.7071067811865475*nuVtSqSum[0]*temp_diff[12]; 
  diff_incr[13] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[23]+0.7071067811865475*nuVtSqSum[0]*temp_diff[13]; 
  diff_incr[14] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[28]+0.7071067811865475*nuVtSqSum[0]*temp_diff[14]; 
  diff_incr[15] = 0.6324555320336759*nuVtSqSum[1]*temp_diff[32]+0.6324555320336759*nuVtSqSum[2]*temp_diff[15]+0.7071067811865475*nuVtSqSum[0]*temp_diff[15]+0.7071067811865475*nuVtSqSum[1]*temp_diff[7]; 
  diff_incr[16] = 0.6324555320336759*nuVtSqSum[1]*temp_diff[35]+0.6324555320336759*nuVtSqSum[2]*temp_diff[16]+0.7071067811865475*nuVtSqSum[0]*temp_diff[16]+0.7071067811865475*nuVtSqSum[1]*temp_diff[9]; 
  diff_incr[17] = 0.6324555320336759*nuVtSqSum[1]*temp_diff[37]+0.6324555320336759*nuVtSqSum[2]*temp_diff[17]+0.7071067811865475*nuVtSqSum[0]*temp_diff[17]+0.7071067811865475*nuVtSqSum[1]*temp_diff[10]; 
  diff_incr[18] = 0.7071067811865475*nuVtSqSum[2]*temp_diff[44]+0.7071067811865475*nuVtSqSum[1]*temp_diff[31]+0.7071067811865475*nuVtSqSum[0]*temp_diff[18]; 
  diff_incr[19] = 0.4517539514526256*nuVtSqSum[2]*temp_diff[19]+0.7071067811865475*nuVtSqSum[0]*temp_diff[19]+0.632455532033676*nuVtSqSum[1]*temp_diff[5]+0.7071067811865475*nuVtSqSum[2]*temp_diff[2]; 
  diff_incr[20] = 0.6324555320336759*nuVtSqSum[2]*temp_diff[20]+0.7071067811865475*nuVtSqSum[0]*temp_diff[20]+0.7071067811865475*nuVtSqSum[1]*temp_diff[12]; 
  diff_incr[21] = 0.4517539514526256*nuVtSqSum[2]*temp_diff[21]+0.7071067811865475*nuVtSqSum[0]*temp_diff[21]+0.632455532033676*nuVtSqSum[1]*temp_diff[6]+0.7071067811865475*nuVtSqSum[2]*temp_diff[3]; 
  diff_incr[22] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[33]+0.7071067811865475*nuVtSqSum[0]*temp_diff[22]; 
  diff_incr[23] = 0.6324555320336759*nuVtSqSum[2]*temp_diff[23]+0.7071067811865475*nuVtSqSum[0]*temp_diff[23]+0.7071067811865475*nuVtSqSum[1]*temp_diff[13]; 
  diff_incr[24] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[34]+0.7071067811865475*nuVtSqSum[0]*temp_diff[24]; 
  diff_incr[25] = 0.4517539514526256*nuVtSqSum[2]*temp_diff[25]+0.7071067811865475*nuVtSqSum[0]*temp_diff[25]+0.632455532033676*nuVtSqSum[1]*temp_diff[8]+0.7071067811865475*nuVtSqSum[2]*temp_diff[4]; 
  diff_incr[26] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[36]+0.7071067811865475*nuVtSqSum[0]*temp_diff[26]; 
  diff_incr[27] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[39]+0.7071067811865475*nuVtSqSum[0]*temp_diff[27]; 
  diff_incr[28] = 0.6324555320336759*nuVtSqSum[2]*temp_diff[28]+0.7071067811865475*nuVtSqSum[0]*temp_diff[28]+0.7071067811865475*nuVtSqSum[1]*temp_diff[14]; 
  diff_incr[29] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[41]+0.7071067811865475*nuVtSqSum[0]*temp_diff[29]; 
  diff_incr[30] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[42]+0.7071067811865475*nuVtSqSum[0]*temp_diff[30]; 
  diff_incr[31] = 0.632455532033676*nuVtSqSum[1]*temp_diff[44]+0.6324555320336759*nuVtSqSum[2]*temp_diff[31]+0.7071067811865475*nuVtSqSum[0]*temp_diff[31]+0.7071067811865475*nuVtSqSum[1]*temp_diff[18]; 
  diff_incr[32] = 0.4517539514526256*nuVtSqSum[2]*temp_diff[32]+0.7071067811865475*nuVtSqSum[0]*temp_diff[32]+0.6324555320336759*nuVtSqSum[1]*temp_diff[15]+0.7071067811865475*nuVtSqSum[2]*temp_diff[7]; 
  diff_incr[33] = 0.6324555320336759*nuVtSqSum[2]*temp_diff[33]+0.7071067811865475*nuVtSqSum[0]*temp_diff[33]+0.7071067811865475*nuVtSqSum[1]*temp_diff[22]; 
  diff_incr[34] = 0.6324555320336759*nuVtSqSum[2]*temp_diff[34]+0.7071067811865475*nuVtSqSum[0]*temp_diff[34]+0.7071067811865475*nuVtSqSum[1]*temp_diff[24]; 
  diff_incr[35] = 0.4517539514526256*nuVtSqSum[2]*temp_diff[35]+0.7071067811865475*nuVtSqSum[0]*temp_diff[35]+0.6324555320336759*nuVtSqSum[1]*temp_diff[16]+0.7071067811865475*nuVtSqSum[2]*temp_diff[9]; 
  diff_incr[36] = 0.6324555320336759*nuVtSqSum[2]*temp_diff[36]+0.7071067811865475*nuVtSqSum[0]*temp_diff[36]+0.7071067811865475*nuVtSqSum[1]*temp_diff[26]; 
  diff_incr[37] = 0.4517539514526256*nuVtSqSum[2]*temp_diff[37]+0.7071067811865475*nuVtSqSum[0]*temp_diff[37]+0.6324555320336759*nuVtSqSum[1]*temp_diff[17]+0.7071067811865475*nuVtSqSum[2]*temp_diff[10]; 
  diff_incr[38] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[45]+0.7071067811865475*nuVtSqSum[0]*temp_diff[38]; 
  diff_incr[39] = 0.6324555320336759*nuVtSqSum[2]*temp_diff[39]+0.7071067811865475*nuVtSqSum[0]*temp_diff[39]+0.7071067811865475*nuVtSqSum[1]*temp_diff[27]; 
  diff_incr[40] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[46]+0.7071067811865475*nuVtSqSum[0]*temp_diff[40]; 
  diff_incr[41] = 0.6324555320336759*nuVtSqSum[2]*temp_diff[41]+0.7071067811865475*nuVtSqSum[0]*temp_diff[41]+0.7071067811865475*nuVtSqSum[1]*temp_diff[29]; 
  diff_incr[42] = 0.6324555320336759*nuVtSqSum[2]*temp_diff[42]+0.7071067811865475*nuVtSqSum[0]*temp_diff[42]+0.7071067811865475*nuVtSqSum[1]*temp_diff[30]; 
  diff_incr[43] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[47]+0.7071067811865475*nuVtSqSum[0]*temp_diff[43]; 
  diff_incr[44] = 0.4517539514526256*nuVtSqSum[2]*temp_diff[44]+0.7071067811865475*nuVtSqSum[0]*temp_diff[44]+0.632455532033676*nuVtSqSum[1]*temp_diff[31]+0.7071067811865475*nuVtSqSum[2]*temp_diff[18]; 
  diff_incr[45] = 0.6324555320336759*nuVtSqSum[2]*temp_diff[45]+0.7071067811865475*nuVtSqSum[0]*temp_diff[45]+0.7071067811865475*nuVtSqSum[1]*temp_diff[38]; 
  diff_incr[46] = 0.6324555320336759*nuVtSqSum[2]*temp_diff[46]+0.7071067811865475*nuVtSqSum[0]*temp_diff[46]+0.7071067811865475*nuVtSqSum[1]*temp_diff[40]; 
  diff_incr[47] = 0.6324555320336759*nuVtSqSum[2]*temp_diff[47]+0.7071067811865475*nuVtSqSum[0]*temp_diff[47]+0.7071067811865475*nuVtSqSum[1]*temp_diff[43]; 

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
  out[40] += diff_incr[40]*rdvSq4; 
  out[41] += diff_incr[41]*rdvSq4; 
  out[42] += diff_incr[42]*rdvSq4; 
  out[43] += diff_incr[43]*rdvSq4; 
  out[44] += diff_incr[44]*rdvSq4; 
  out[45] += diff_incr[45]*rdvSq4; 
  out[46] += diff_incr[46]*rdvSq4; 
  out[47] += diff_incr[47]*rdvSq4; 
} 