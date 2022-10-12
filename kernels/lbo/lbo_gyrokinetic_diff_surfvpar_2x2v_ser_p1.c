#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH void lbo_gyrokinetic_diff_surfvpar_2x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // w[4]: cell-center coordinates. 
  // dxv[4]: cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[8]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fl/fc/fr: distribution function in cells 
  // out: incremented distribution function in cell 

  const double *nuVtSqSum = &nuPrimMomsSum[4];

  double rdvSq4 = 4.0/(dxv[2]*dxv[2]); 
  double temp_diff[24] = {0.0}; 
  double diff_incr[24] = {0.0}; 

  temp_diff[0] = 0.6708203932499369*fr[16]+0.6708203932499369*fl[16]-1.341640786499874*fc[16]-1.190784930203603*fr[3]+1.190784930203603*fl[3]+0.9375*fr[0]+0.9375*fl[0]-1.875*fc[0]; 
  temp_diff[1] = 0.6708203932499369*fr[17]+0.6708203932499369*fl[17]-1.341640786499874*fc[17]-1.190784930203603*fr[6]+1.190784930203603*fl[6]+0.9375*fr[1]+0.9375*fl[1]-1.875*fc[1]; 
  temp_diff[2] = 0.6708203932499369*fr[18]+0.6708203932499369*fl[18]-1.341640786499874*fc[18]-1.190784930203603*fr[7]+1.190784930203603*fl[7]+0.9375*fr[2]+0.9375*fl[2]-1.875*fc[2]; 
  temp_diff[3] = 0.7382874503707888*fr[16]-0.7382874503707888*fl[16]-1.453125*fr[3]-1.453125*fl[3]-5.34375*fc[3]+1.190784930203603*fr[0]-1.190784930203603*fl[0]; 
  temp_diff[4] = 0.6708203932499369*fr[19]+0.6708203932499369*fl[19]-1.341640786499874*fc[19]-1.190784930203603*fr[10]+1.190784930203603*fl[10]+0.9375*fr[4]+0.9375*fl[4]-1.875*fc[4]; 
  temp_diff[5] = 0.6708203932499369*fr[20]+0.6708203932499369*fl[20]-1.341640786499874*fc[20]-1.190784930203603*fr[11]+1.190784930203603*fl[11]+0.9375*fr[5]+0.9375*fl[5]-1.875*fc[5]; 
  temp_diff[6] = 0.7382874503707888*fr[17]-0.7382874503707888*fl[17]-1.453125*fr[6]-1.453125*fl[6]-5.34375*fc[6]+1.190784930203603*fr[1]-1.190784930203603*fl[1]; 
  temp_diff[7] = 0.7382874503707888*fr[18]-0.7382874503707888*fl[18]-1.453125*fr[7]-1.453125*fl[7]-5.34375*fc[7]+1.190784930203603*fr[2]-1.190784930203603*fl[2]; 
  temp_diff[8] = 0.6708203932499369*fr[21]+0.6708203932499369*fl[21]-1.341640786499874*fc[21]-1.190784930203603*fr[13]+1.190784930203603*fl[13]+0.9375*fr[8]+0.9375*fl[8]-1.875*fc[8]; 
  temp_diff[9] = 0.6708203932499369*fr[22]+0.6708203932499369*fl[22]-1.341640786499874*fc[22]-1.190784930203603*fr[14]+1.190784930203603*fl[14]+0.9375*fr[9]+0.9375*fl[9]-1.875*fc[9]; 
  temp_diff[10] = 0.7382874503707888*fr[19]-0.7382874503707888*fl[19]-1.453125*fr[10]-1.453125*fl[10]-5.34375*fc[10]+1.190784930203603*fr[4]-1.190784930203603*fl[4]; 
  temp_diff[11] = 0.7382874503707888*fr[20]-0.7382874503707888*fl[20]-1.453125*fr[11]-1.453125*fl[11]-5.34375*fc[11]+1.190784930203603*fr[5]-1.190784930203603*fl[5]; 
  temp_diff[12] = 0.6708203932499369*fr[23]+0.6708203932499369*fl[23]-1.341640786499874*fc[23]-1.190784930203603*fr[15]+1.190784930203603*fl[15]+0.9375*fr[12]+0.9375*fl[12]-1.875*fc[12]; 
  temp_diff[13] = 0.7382874503707888*fr[21]-0.7382874503707888*fl[21]-1.453125*fr[13]-1.453125*fl[13]-5.34375*fc[13]+1.190784930203603*fr[8]-1.190784930203603*fl[8]; 
  temp_diff[14] = 0.7382874503707888*fr[22]-0.7382874503707888*fl[22]-1.453125*fr[14]-1.453125*fl[14]-5.34375*fc[14]+1.190784930203603*fr[9]-1.190784930203603*fl[9]; 
  temp_diff[15] = 0.7382874503707888*fr[23]-0.7382874503707888*fl[23]-1.453125*fr[15]-1.453125*fl[15]-5.34375*fc[15]+1.190784930203603*fr[12]-1.190784930203603*fl[12]; 
  temp_diff[16] = (-0.140625*fr[16])-0.140625*fl[16]-6.28125*fc[16]-0.3025768239224545*fr[3]+0.3025768239224545*fl[3]+0.4192627457812106*fr[0]+0.4192627457812106*fl[0]-0.8385254915624212*fc[0]; 
  temp_diff[17] = (-0.140625*fr[17])-0.140625*fl[17]-6.28125*fc[17]-0.3025768239224544*fr[6]+0.3025768239224544*fl[6]+0.4192627457812105*fr[1]+0.4192627457812105*fl[1]-0.8385254915624211*fc[1]; 
  temp_diff[18] = (-0.140625*fr[18])-0.140625*fl[18]-6.28125*fc[18]-0.3025768239224544*fr[7]+0.3025768239224544*fl[7]+0.4192627457812105*fr[2]+0.4192627457812105*fl[2]-0.8385254915624211*fc[2]; 
  temp_diff[19] = (-0.140625*fr[19])-0.140625*fl[19]-6.28125*fc[19]-0.3025768239224544*fr[10]+0.3025768239224544*fl[10]+0.4192627457812105*fr[4]+0.4192627457812105*fl[4]-0.8385254915624211*fc[4]; 
  temp_diff[20] = (-0.140625*fr[20])-0.140625*fl[20]-6.28125*fc[20]-0.3025768239224545*fr[11]+0.3025768239224545*fl[11]+0.4192627457812106*fr[5]+0.4192627457812106*fl[5]-0.8385254915624212*fc[5]; 
  temp_diff[21] = (-0.140625*fr[21])-0.140625*fl[21]-6.28125*fc[21]-0.3025768239224545*fr[13]+0.3025768239224545*fl[13]+0.4192627457812106*fr[8]+0.4192627457812106*fl[8]-0.8385254915624212*fc[8]; 
  temp_diff[22] = (-0.140625*fr[22])-0.140625*fl[22]-6.28125*fc[22]-0.3025768239224545*fr[14]+0.3025768239224545*fl[14]+0.4192627457812106*fr[9]+0.4192627457812106*fl[9]-0.8385254915624212*fc[9]; 
  temp_diff[23] = (-0.140625*fr[23])-0.140625*fl[23]-6.28125*fc[23]-0.3025768239224544*fr[15]+0.3025768239224544*fl[15]+0.4192627457812105*fr[12]+0.4192627457812105*fl[12]-0.8385254915624211*fc[12]; 

  diff_incr[0] = 0.5*nuVtSqSum[3]*temp_diff[5]+0.5*nuVtSqSum[2]*temp_diff[2]+0.5*nuVtSqSum[1]*temp_diff[1]+0.5*nuVtSqSum[0]*temp_diff[0]; 
  diff_incr[1] = 0.5*nuVtSqSum[2]*temp_diff[5]+0.5*temp_diff[2]*nuVtSqSum[3]+0.5*nuVtSqSum[0]*temp_diff[1]+0.5*temp_diff[0]*nuVtSqSum[1]; 
  diff_incr[2] = 0.5*nuVtSqSum[1]*temp_diff[5]+0.5*temp_diff[1]*nuVtSqSum[3]+0.5*nuVtSqSum[0]*temp_diff[2]+0.5*temp_diff[0]*nuVtSqSum[2]; 
  diff_incr[3] = 0.5*nuVtSqSum[3]*temp_diff[11]+0.5*nuVtSqSum[2]*temp_diff[7]+0.5*nuVtSqSum[1]*temp_diff[6]+0.5*nuVtSqSum[0]*temp_diff[3]; 
  diff_incr[4] = 0.5*nuVtSqSum[3]*temp_diff[12]+0.5*nuVtSqSum[2]*temp_diff[9]+0.5*nuVtSqSum[1]*temp_diff[8]+0.5*nuVtSqSum[0]*temp_diff[4]; 
  diff_incr[5] = 0.5*nuVtSqSum[0]*temp_diff[5]+0.5*temp_diff[0]*nuVtSqSum[3]+0.5*nuVtSqSum[1]*temp_diff[2]+0.5*temp_diff[1]*nuVtSqSum[2]; 
  diff_incr[6] = 0.5*nuVtSqSum[2]*temp_diff[11]+0.5*nuVtSqSum[3]*temp_diff[7]+0.5*nuVtSqSum[0]*temp_diff[6]+0.5*nuVtSqSum[1]*temp_diff[3]; 
  diff_incr[7] = 0.5*nuVtSqSum[1]*temp_diff[11]+0.5*nuVtSqSum[0]*temp_diff[7]+0.5*nuVtSqSum[3]*temp_diff[6]+0.5*nuVtSqSum[2]*temp_diff[3]; 
  diff_incr[8] = 0.5*nuVtSqSum[2]*temp_diff[12]+0.5*nuVtSqSum[3]*temp_diff[9]+0.5*nuVtSqSum[0]*temp_diff[8]+0.5*nuVtSqSum[1]*temp_diff[4]; 
  diff_incr[9] = 0.5*nuVtSqSum[1]*temp_diff[12]+0.5*nuVtSqSum[0]*temp_diff[9]+0.5*nuVtSqSum[3]*temp_diff[8]+0.5*nuVtSqSum[2]*temp_diff[4]; 
  diff_incr[10] = 0.5*nuVtSqSum[3]*temp_diff[15]+0.5*nuVtSqSum[2]*temp_diff[14]+0.5*nuVtSqSum[1]*temp_diff[13]+0.5*nuVtSqSum[0]*temp_diff[10]; 
  diff_incr[11] = 0.5*nuVtSqSum[0]*temp_diff[11]+0.5*nuVtSqSum[1]*temp_diff[7]+0.5*nuVtSqSum[2]*temp_diff[6]+0.5*nuVtSqSum[3]*temp_diff[3]; 
  diff_incr[12] = 0.5*nuVtSqSum[0]*temp_diff[12]+0.5*nuVtSqSum[1]*temp_diff[9]+0.5*nuVtSqSum[2]*temp_diff[8]+0.5*nuVtSqSum[3]*temp_diff[4]; 
  diff_incr[13] = 0.5*nuVtSqSum[2]*temp_diff[15]+0.5*nuVtSqSum[3]*temp_diff[14]+0.5*nuVtSqSum[0]*temp_diff[13]+0.5*nuVtSqSum[1]*temp_diff[10]; 
  diff_incr[14] = 0.5*nuVtSqSum[1]*temp_diff[15]+0.5*nuVtSqSum[0]*temp_diff[14]+0.5*nuVtSqSum[3]*temp_diff[13]+0.5*nuVtSqSum[2]*temp_diff[10]; 
  diff_incr[15] = 0.5*nuVtSqSum[0]*temp_diff[15]+0.5*nuVtSqSum[1]*temp_diff[14]+0.5*nuVtSqSum[2]*temp_diff[13]+0.5*nuVtSqSum[3]*temp_diff[10]; 
  diff_incr[16] = 0.5*nuVtSqSum[3]*temp_diff[20]+0.5000000000000001*nuVtSqSum[2]*temp_diff[18]+0.5000000000000001*nuVtSqSum[1]*temp_diff[17]+0.5*nuVtSqSum[0]*temp_diff[16]; 
  diff_incr[17] = 0.5000000000000001*nuVtSqSum[2]*temp_diff[20]+0.5*nuVtSqSum[3]*temp_diff[18]+0.5*nuVtSqSum[0]*temp_diff[17]+0.5000000000000001*nuVtSqSum[1]*temp_diff[16]; 
  diff_incr[18] = 0.5000000000000001*nuVtSqSum[1]*temp_diff[20]+0.5*nuVtSqSum[0]*temp_diff[18]+0.5*nuVtSqSum[3]*temp_diff[17]+0.5000000000000001*nuVtSqSum[2]*temp_diff[16]; 
  diff_incr[19] = 0.5*nuVtSqSum[3]*temp_diff[23]+0.5000000000000001*nuVtSqSum[2]*temp_diff[22]+0.5000000000000001*nuVtSqSum[1]*temp_diff[21]+0.5*nuVtSqSum[0]*temp_diff[19]; 
  diff_incr[20] = 0.5*nuVtSqSum[0]*temp_diff[20]+0.5000000000000001*nuVtSqSum[1]*temp_diff[18]+0.5000000000000001*nuVtSqSum[2]*temp_diff[17]+0.5*nuVtSqSum[3]*temp_diff[16]; 
  diff_incr[21] = 0.5000000000000001*nuVtSqSum[2]*temp_diff[23]+0.5*nuVtSqSum[3]*temp_diff[22]+0.5*nuVtSqSum[0]*temp_diff[21]+0.5000000000000001*nuVtSqSum[1]*temp_diff[19]; 
  diff_incr[22] = 0.5000000000000001*nuVtSqSum[1]*temp_diff[23]+0.5*nuVtSqSum[0]*temp_diff[22]+0.5*nuVtSqSum[3]*temp_diff[21]+0.5000000000000001*nuVtSqSum[2]*temp_diff[19]; 
  diff_incr[23] = 0.5*nuVtSqSum[0]*temp_diff[23]+0.5000000000000001*nuVtSqSum[1]*temp_diff[22]+0.5000000000000001*nuVtSqSum[2]*temp_diff[21]+0.5*nuVtSqSum[3]*temp_diff[19]; 

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
} 
