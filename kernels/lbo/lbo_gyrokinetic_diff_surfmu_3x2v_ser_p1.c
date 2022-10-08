#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH void lbo_gyrokinetic_diff_surfmu_3x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // w[5]: cell-center coordinates. 
  // dxv[5]: cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[16]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fl/fc/fr: distribution function in cells 
  // out: incremented distribution function in cell 

  const double *nuVtSqSum = &nuPrimMomsSum[8];

  double rdvSq4 = 4.0/(dxv[4]*dxv[4]); 
  double temp_diff[48] = {0.0}; 
  double diff_incr[48] = {0.0}; 

  temp_diff[0] = (-0.5412658773652741*w[4]*fr[5])-0.270632938682637*dxv[4]*fr[5]+0.5412658773652741*w[4]*fl[5]-0.270632938682637*dxv[4]*fl[5]-0.5412658773652741*dxv[4]*fc[5]+0.5625*fr[0]*w[4]+0.5625*fl[0]*w[4]-1.125*fc[0]*w[4]+0.28125*fr[0]*dxv[4]-0.28125*fl[0]*dxv[4]; 
  temp_diff[1] = (-0.5412658773652741*w[4]*fr[12])-0.270632938682637*dxv[4]*fr[12]+0.5412658773652741*w[4]*fl[12]-0.270632938682637*dxv[4]*fl[12]-0.5412658773652741*dxv[4]*fc[12]+0.5625*fr[1]*w[4]+0.5625*fl[1]*w[4]-1.125*fc[1]*w[4]+0.28125*fr[1]*dxv[4]-0.28125*fl[1]*dxv[4]; 
  temp_diff[2] = (-0.5412658773652741*w[4]*fr[13])-0.270632938682637*dxv[4]*fr[13]+0.5412658773652741*w[4]*fl[13]-0.270632938682637*dxv[4]*fl[13]-0.5412658773652741*dxv[4]*fc[13]+0.5625*fr[2]*w[4]+0.5625*fl[2]*w[4]-1.125*fc[2]*w[4]+0.28125*fr[2]*dxv[4]-0.28125*fl[2]*dxv[4]; 
  temp_diff[3] = (-0.5412658773652741*w[4]*fr[14])-0.270632938682637*dxv[4]*fr[14]+0.5412658773652741*w[4]*fl[14]-0.270632938682637*dxv[4]*fl[14]-0.5412658773652741*dxv[4]*fc[14]+0.5625*fr[3]*w[4]+0.5625*fl[3]*w[4]-1.125*fc[3]*w[4]+0.28125*fr[3]*dxv[4]-0.28125*fl[3]*dxv[4]; 
  temp_diff[4] = (-0.5412658773652741*w[4]*fr[15])-0.270632938682637*dxv[4]*fr[15]+0.5412658773652741*w[4]*fl[15]-0.270632938682637*dxv[4]*fl[15]-0.5412658773652741*dxv[4]*fc[15]+0.5625*fr[4]*w[4]+0.5625*fl[4]*w[4]-1.125*fc[4]*w[4]+0.28125*dxv[4]*fr[4]-0.28125*dxv[4]*fl[4]; 
  temp_diff[5] = (-0.4375*w[4]*fr[5])-0.21875*dxv[4]*fr[5]-0.4375*w[4]*fl[5]+0.21875*dxv[4]*fl[5]-2.875*w[4]*fc[5]+0.5412658773652741*fr[0]*w[4]-0.5412658773652741*fl[0]*w[4]+0.270632938682637*fr[0]*dxv[4]+0.270632938682637*fl[0]*dxv[4]-0.5412658773652741*fc[0]*dxv[4]; 
  temp_diff[6] = (-0.5412658773652741*w[4]*fr[20])-0.270632938682637*dxv[4]*fr[20]+0.5412658773652741*w[4]*fl[20]-0.270632938682637*dxv[4]*fl[20]-0.5412658773652741*dxv[4]*fc[20]+0.5625*w[4]*fr[6]+0.28125*dxv[4]*fr[6]+0.5625*w[4]*fl[6]-0.28125*dxv[4]*fl[6]-1.125*w[4]*fc[6]; 
  temp_diff[7] = (-0.5412658773652741*w[4]*fr[21])-0.270632938682637*dxv[4]*fr[21]+0.5412658773652741*w[4]*fl[21]-0.270632938682637*dxv[4]*fl[21]-0.5412658773652741*dxv[4]*fc[21]+0.5625*w[4]*fr[7]+0.28125*dxv[4]*fr[7]+0.5625*w[4]*fl[7]-0.28125*dxv[4]*fl[7]-1.125*w[4]*fc[7]; 
  temp_diff[8] = (-0.5412658773652741*w[4]*fr[22])-0.270632938682637*dxv[4]*fr[22]+0.5412658773652741*w[4]*fl[22]-0.270632938682637*dxv[4]*fl[22]-0.5412658773652741*dxv[4]*fc[22]+0.5625*w[4]*fr[8]+0.28125*dxv[4]*fr[8]+0.5625*w[4]*fl[8]-0.28125*dxv[4]*fl[8]-1.125*w[4]*fc[8]; 
  temp_diff[9] = (-0.5412658773652741*w[4]*fr[23])-0.270632938682637*dxv[4]*fr[23]+0.5412658773652741*w[4]*fl[23]-0.270632938682637*dxv[4]*fl[23]-0.5412658773652741*dxv[4]*fc[23]+0.5625*w[4]*fr[9]+0.28125*dxv[4]*fr[9]+0.5625*w[4]*fl[9]-0.28125*dxv[4]*fl[9]-1.125*w[4]*fc[9]; 
  temp_diff[10] = (-0.5412658773652741*w[4]*fr[24])-0.270632938682637*dxv[4]*fr[24]+0.5412658773652741*w[4]*fl[24]-0.270632938682637*dxv[4]*fl[24]-0.5412658773652741*dxv[4]*fc[24]+0.5625*w[4]*fr[10]+0.28125*dxv[4]*fr[10]+0.5625*w[4]*fl[10]-0.28125*dxv[4]*fl[10]-1.125*w[4]*fc[10]; 
  temp_diff[11] = (-0.5412658773652741*w[4]*fr[25])-0.270632938682637*dxv[4]*fr[25]+0.5412658773652741*w[4]*fl[25]-0.270632938682637*dxv[4]*fl[25]-0.5412658773652741*dxv[4]*fc[25]+0.5625*w[4]*fr[11]+0.28125*dxv[4]*fr[11]+0.5625*w[4]*fl[11]-0.28125*dxv[4]*fl[11]-1.125*w[4]*fc[11]; 
  temp_diff[12] = (-0.4375*w[4]*fr[12])-0.21875*dxv[4]*fr[12]-0.4375*w[4]*fl[12]+0.21875*dxv[4]*fl[12]-2.875*w[4]*fc[12]+0.5412658773652741*fr[1]*w[4]-0.5412658773652741*fl[1]*w[4]+0.270632938682637*fr[1]*dxv[4]+0.270632938682637*fl[1]*dxv[4]-0.5412658773652741*fc[1]*dxv[4]; 
  temp_diff[13] = (-0.4375*w[4]*fr[13])-0.21875*dxv[4]*fr[13]-0.4375*w[4]*fl[13]+0.21875*dxv[4]*fl[13]-2.875*w[4]*fc[13]+0.5412658773652741*fr[2]*w[4]-0.5412658773652741*fl[2]*w[4]+0.270632938682637*fr[2]*dxv[4]+0.270632938682637*fl[2]*dxv[4]-0.5412658773652741*fc[2]*dxv[4]; 
  temp_diff[14] = (-0.4375*w[4]*fr[14])-0.21875*dxv[4]*fr[14]-0.4375*w[4]*fl[14]+0.21875*dxv[4]*fl[14]-2.875*w[4]*fc[14]+0.5412658773652741*fr[3]*w[4]-0.5412658773652741*fl[3]*w[4]+0.270632938682637*fr[3]*dxv[4]+0.270632938682637*fl[3]*dxv[4]-0.5412658773652741*fc[3]*dxv[4]; 
  temp_diff[15] = (-0.4375*w[4]*fr[15])-0.21875*dxv[4]*fr[15]-0.4375*w[4]*fl[15]+0.21875*dxv[4]*fl[15]-2.875*w[4]*fc[15]+0.5412658773652741*fr[4]*w[4]-0.5412658773652741*fl[4]*w[4]+0.270632938682637*dxv[4]*fr[4]+0.270632938682637*dxv[4]*fl[4]-0.5412658773652741*dxv[4]*fc[4]; 
  temp_diff[16] = (-0.5412658773652741*w[4]*fr[27])-0.270632938682637*dxv[4]*fr[27]+0.5412658773652741*w[4]*fl[27]-0.270632938682637*dxv[4]*fl[27]-0.5412658773652741*dxv[4]*fc[27]+0.5625*w[4]*fr[16]+0.28125*dxv[4]*fr[16]+0.5625*w[4]*fl[16]-0.28125*dxv[4]*fl[16]-1.125*w[4]*fc[16]; 
  temp_diff[17] = (-0.5412658773652741*w[4]*fr[28])-0.270632938682637*dxv[4]*fr[28]+0.5412658773652741*w[4]*fl[28]-0.270632938682637*dxv[4]*fl[28]-0.5412658773652741*dxv[4]*fc[28]+0.5625*w[4]*fr[17]+0.28125*dxv[4]*fr[17]+0.5625*w[4]*fl[17]-0.28125*dxv[4]*fl[17]-1.125*w[4]*fc[17]; 
  temp_diff[18] = (-0.5412658773652741*w[4]*fr[29])-0.270632938682637*dxv[4]*fr[29]+0.5412658773652741*w[4]*fl[29]-0.270632938682637*dxv[4]*fl[29]-0.5412658773652741*dxv[4]*fc[29]+0.5625*w[4]*fr[18]+0.28125*dxv[4]*fr[18]+0.5625*w[4]*fl[18]-0.28125*dxv[4]*fl[18]-1.125*w[4]*fc[18]; 
  temp_diff[19] = (-0.5412658773652741*w[4]*fr[30])-0.270632938682637*dxv[4]*fr[30]+0.5412658773652741*w[4]*fl[30]-0.270632938682637*dxv[4]*fl[30]-0.5412658773652741*dxv[4]*fc[30]+0.5625*w[4]*fr[19]+0.28125*dxv[4]*fr[19]+0.5625*w[4]*fl[19]-0.28125*dxv[4]*fl[19]-1.125*w[4]*fc[19]; 
  temp_diff[20] = (-0.4375*w[4]*fr[20])-0.21875*dxv[4]*fr[20]-0.4375*w[4]*fl[20]+0.21875*dxv[4]*fl[20]-2.875*w[4]*fc[20]+0.5412658773652741*w[4]*fr[6]+0.270632938682637*dxv[4]*fr[6]-0.5412658773652741*w[4]*fl[6]+0.270632938682637*dxv[4]*fl[6]-0.5412658773652741*dxv[4]*fc[6]; 
  temp_diff[21] = (-0.4375*w[4]*fr[21])-0.21875*dxv[4]*fr[21]-0.4375*w[4]*fl[21]+0.21875*dxv[4]*fl[21]-2.875*w[4]*fc[21]+0.5412658773652741*w[4]*fr[7]+0.270632938682637*dxv[4]*fr[7]-0.5412658773652741*w[4]*fl[7]+0.270632938682637*dxv[4]*fl[7]-0.5412658773652741*dxv[4]*fc[7]; 
  temp_diff[22] = (-0.4375*w[4]*fr[22])-0.21875*dxv[4]*fr[22]-0.4375*w[4]*fl[22]+0.21875*dxv[4]*fl[22]-2.875*w[4]*fc[22]+0.5412658773652741*w[4]*fr[8]+0.270632938682637*dxv[4]*fr[8]-0.5412658773652741*w[4]*fl[8]+0.270632938682637*dxv[4]*fl[8]-0.5412658773652741*dxv[4]*fc[8]; 
  temp_diff[23] = (-0.4375*w[4]*fr[23])-0.21875*dxv[4]*fr[23]-0.4375*w[4]*fl[23]+0.21875*dxv[4]*fl[23]-2.875*w[4]*fc[23]+0.5412658773652741*w[4]*fr[9]+0.270632938682637*dxv[4]*fr[9]-0.5412658773652741*w[4]*fl[9]+0.270632938682637*dxv[4]*fl[9]-0.5412658773652741*dxv[4]*fc[9]; 
  temp_diff[24] = (-0.4375*w[4]*fr[24])-0.21875*dxv[4]*fr[24]-0.4375*w[4]*fl[24]+0.21875*dxv[4]*fl[24]-2.875*w[4]*fc[24]+0.5412658773652741*w[4]*fr[10]+0.270632938682637*dxv[4]*fr[10]-0.5412658773652741*w[4]*fl[10]+0.270632938682637*dxv[4]*fl[10]-0.5412658773652741*dxv[4]*fc[10]; 
  temp_diff[25] = (-0.4375*w[4]*fr[25])-0.21875*dxv[4]*fr[25]-0.4375*w[4]*fl[25]+0.21875*dxv[4]*fl[25]-2.875*w[4]*fc[25]+0.5412658773652741*w[4]*fr[11]+0.270632938682637*dxv[4]*fr[11]-0.5412658773652741*w[4]*fl[11]+0.270632938682637*dxv[4]*fl[11]-0.5412658773652741*dxv[4]*fc[11]; 
  temp_diff[26] = (-0.5412658773652741*w[4]*fr[31])-0.270632938682637*dxv[4]*fr[31]+0.5412658773652741*w[4]*fl[31]-0.270632938682637*dxv[4]*fl[31]-0.5412658773652741*dxv[4]*fc[31]+0.5625*w[4]*fr[26]+0.28125*dxv[4]*fr[26]+0.5625*w[4]*fl[26]-0.28125*dxv[4]*fl[26]-1.125*w[4]*fc[26]; 
  temp_diff[27] = (-0.4375*w[4]*fr[27])-0.21875*dxv[4]*fr[27]-0.4375*w[4]*fl[27]+0.21875*dxv[4]*fl[27]-2.875*w[4]*fc[27]+0.5412658773652741*w[4]*fr[16]+0.270632938682637*dxv[4]*fr[16]-0.5412658773652741*w[4]*fl[16]+0.270632938682637*dxv[4]*fl[16]-0.5412658773652741*dxv[4]*fc[16]; 
  temp_diff[28] = (-0.4375*w[4]*fr[28])-0.21875*dxv[4]*fr[28]-0.4375*w[4]*fl[28]+0.21875*dxv[4]*fl[28]-2.875*w[4]*fc[28]+0.5412658773652741*w[4]*fr[17]+0.270632938682637*dxv[4]*fr[17]-0.5412658773652741*w[4]*fl[17]+0.270632938682637*dxv[4]*fl[17]-0.5412658773652741*dxv[4]*fc[17]; 
  temp_diff[29] = (-0.4375*w[4]*fr[29])-0.21875*dxv[4]*fr[29]-0.4375*w[4]*fl[29]+0.21875*dxv[4]*fl[29]-2.875*w[4]*fc[29]+0.5412658773652741*w[4]*fr[18]+0.270632938682637*dxv[4]*fr[18]-0.5412658773652741*w[4]*fl[18]+0.270632938682637*dxv[4]*fl[18]-0.5412658773652741*dxv[4]*fc[18]; 
  temp_diff[30] = (-0.4375*w[4]*fr[30])-0.21875*dxv[4]*fr[30]-0.4375*w[4]*fl[30]+0.21875*dxv[4]*fl[30]-2.875*w[4]*fc[30]+0.5412658773652741*w[4]*fr[19]+0.270632938682637*dxv[4]*fr[19]-0.5412658773652741*w[4]*fl[19]+0.270632938682637*dxv[4]*fl[19]-0.5412658773652741*dxv[4]*fc[19]; 
  temp_diff[31] = (-0.4375*w[4]*fr[31])-0.21875*dxv[4]*fr[31]-0.4375*w[4]*fl[31]+0.21875*dxv[4]*fl[31]-2.875*w[4]*fc[31]+0.5412658773652741*w[4]*fr[26]+0.270632938682637*dxv[4]*fr[26]-0.5412658773652741*w[4]*fl[26]+0.270632938682637*dxv[4]*fl[26]-0.5412658773652741*dxv[4]*fc[26]; 
  temp_diff[32] = (-0.5412658773652742*w[4]*fr[36])-0.2706329386826371*dxv[4]*fr[36]+0.5412658773652742*w[4]*fl[36]-0.2706329386826371*dxv[4]*fl[36]-0.5412658773652742*dxv[4]*fc[36]+0.5625*w[4]*fr[32]+0.28125*dxv[4]*fr[32]+0.5625*w[4]*fl[32]-0.28125*dxv[4]*fl[32]-1.125*w[4]*fc[32]; 
  temp_diff[33] = (-0.5412658773652742*w[4]*fr[40])-0.2706329386826371*dxv[4]*fr[40]+0.5412658773652742*w[4]*fl[40]-0.2706329386826371*dxv[4]*fl[40]-0.5412658773652742*dxv[4]*fc[40]+0.5625*w[4]*fr[33]+0.28125*dxv[4]*fr[33]+0.5625*w[4]*fl[33]-0.28125*dxv[4]*fl[33]-1.125*w[4]*fc[33]; 
  temp_diff[34] = (-0.5412658773652742*w[4]*fr[41])-0.2706329386826371*dxv[4]*fr[41]+0.5412658773652742*w[4]*fl[41]-0.2706329386826371*dxv[4]*fl[41]-0.5412658773652742*dxv[4]*fc[41]+0.5625*w[4]*fr[34]+0.28125*dxv[4]*fr[34]+0.5625*w[4]*fl[34]-0.28125*dxv[4]*fl[34]-1.125*w[4]*fc[34]; 
  temp_diff[35] = (-0.5412658773652742*w[4]*fr[42])-0.2706329386826371*dxv[4]*fr[42]+0.5412658773652742*w[4]*fl[42]-0.2706329386826371*dxv[4]*fl[42]-0.5412658773652742*dxv[4]*fc[42]+0.5625*w[4]*fr[35]+0.28125*dxv[4]*fr[35]+0.5625*w[4]*fl[35]-0.28125*dxv[4]*fl[35]-1.125*w[4]*fc[35]; 
  temp_diff[36] = (-0.4375*w[4]*fr[36])-0.21875*dxv[4]*fr[36]-0.4375*w[4]*fl[36]+0.21875*dxv[4]*fl[36]-2.875*w[4]*fc[36]+0.5412658773652742*w[4]*fr[32]+0.2706329386826371*dxv[4]*fr[32]-0.5412658773652742*w[4]*fl[32]+0.2706329386826371*dxv[4]*fl[32]-0.5412658773652742*dxv[4]*fc[32]; 
  temp_diff[37] = (-0.5412658773652742*w[4]*fr[44])-0.2706329386826371*dxv[4]*fr[44]+0.5412658773652742*w[4]*fl[44]-0.2706329386826371*dxv[4]*fl[44]-0.5412658773652742*dxv[4]*fc[44]+0.5625*w[4]*fr[37]+0.28125*dxv[4]*fr[37]+0.5625*w[4]*fl[37]-0.28125*dxv[4]*fl[37]-1.125*w[4]*fc[37]; 
  temp_diff[38] = (-0.5412658773652742*w[4]*fr[45])-0.2706329386826371*dxv[4]*fr[45]+0.5412658773652742*w[4]*fl[45]-0.2706329386826371*dxv[4]*fl[45]-0.5412658773652742*dxv[4]*fc[45]+0.5625*w[4]*fr[38]+0.28125*dxv[4]*fr[38]+0.5625*w[4]*fl[38]-0.28125*dxv[4]*fl[38]-1.125*w[4]*fc[38]; 
  temp_diff[39] = (-0.5412658773652742*w[4]*fr[46])-0.2706329386826371*dxv[4]*fr[46]+0.5412658773652742*w[4]*fl[46]-0.2706329386826371*dxv[4]*fl[46]-0.5412658773652742*dxv[4]*fc[46]+0.5625*w[4]*fr[39]+0.28125*dxv[4]*fr[39]+0.5625*w[4]*fl[39]-0.28125*dxv[4]*fl[39]-1.125*w[4]*fc[39]; 
  temp_diff[40] = (-0.4375*w[4]*fr[40])-0.21875*dxv[4]*fr[40]-0.4375*w[4]*fl[40]+0.21875*dxv[4]*fl[40]-2.875*w[4]*fc[40]+0.5412658773652742*w[4]*fr[33]+0.2706329386826371*dxv[4]*fr[33]-0.5412658773652742*w[4]*fl[33]+0.2706329386826371*dxv[4]*fl[33]-0.5412658773652742*dxv[4]*fc[33]; 
  temp_diff[41] = (-0.4375*w[4]*fr[41])-0.21875*dxv[4]*fr[41]-0.4375*w[4]*fl[41]+0.21875*dxv[4]*fl[41]-2.875*w[4]*fc[41]+0.5412658773652742*w[4]*fr[34]+0.2706329386826371*dxv[4]*fr[34]-0.5412658773652742*w[4]*fl[34]+0.2706329386826371*dxv[4]*fl[34]-0.5412658773652742*dxv[4]*fc[34]; 
  temp_diff[42] = (-0.4375*w[4]*fr[42])-0.21875*dxv[4]*fr[42]-0.4375*w[4]*fl[42]+0.21875*dxv[4]*fl[42]-2.875*w[4]*fc[42]+0.5412658773652742*w[4]*fr[35]+0.2706329386826371*dxv[4]*fr[35]-0.5412658773652742*w[4]*fl[35]+0.2706329386826371*dxv[4]*fl[35]-0.5412658773652742*dxv[4]*fc[35]; 
  temp_diff[43] = (-0.5412658773652742*w[4]*fr[47])-0.2706329386826371*dxv[4]*fr[47]+0.5412658773652742*w[4]*fl[47]-0.2706329386826371*dxv[4]*fl[47]-0.5412658773652742*dxv[4]*fc[47]+0.5625*w[4]*fr[43]+0.28125*dxv[4]*fr[43]+0.5625*w[4]*fl[43]-0.28125*dxv[4]*fl[43]-1.125*w[4]*fc[43]; 
  temp_diff[44] = (-0.4375*w[4]*fr[44])-0.21875*dxv[4]*fr[44]-0.4375*w[4]*fl[44]+0.21875*dxv[4]*fl[44]-2.875*w[4]*fc[44]+0.5412658773652742*w[4]*fr[37]+0.2706329386826371*dxv[4]*fr[37]-0.5412658773652742*w[4]*fl[37]+0.2706329386826371*dxv[4]*fl[37]-0.5412658773652742*dxv[4]*fc[37]; 
  temp_diff[45] = (-0.4375*w[4]*fr[45])-0.21875*dxv[4]*fr[45]-0.4375*w[4]*fl[45]+0.21875*dxv[4]*fl[45]-2.875*w[4]*fc[45]+0.5412658773652742*w[4]*fr[38]+0.2706329386826371*dxv[4]*fr[38]-0.5412658773652742*w[4]*fl[38]+0.2706329386826371*dxv[4]*fl[38]-0.5412658773652742*dxv[4]*fc[38]; 
  temp_diff[46] = (-0.4375*w[4]*fr[46])-0.21875*dxv[4]*fr[46]-0.4375*w[4]*fl[46]+0.21875*dxv[4]*fl[46]-2.875*w[4]*fc[46]+0.5412658773652742*w[4]*fr[39]+0.2706329386826371*dxv[4]*fr[39]-0.5412658773652742*w[4]*fl[39]+0.2706329386826371*dxv[4]*fl[39]-0.5412658773652742*dxv[4]*fc[39]; 
  temp_diff[47] = (-0.4375*w[4]*fr[47])-0.21875*dxv[4]*fr[47]-0.4375*w[4]*fl[47]+0.21875*dxv[4]*fl[47]-2.875*w[4]*fc[47]+0.5412658773652742*w[4]*fr[43]+0.2706329386826371*dxv[4]*fr[43]-0.5412658773652742*w[4]*fl[43]+0.2706329386826371*dxv[4]*fl[43]-0.5412658773652742*dxv[4]*fc[43]; 

  double diffFac[8] = {0.}; 
  diffFac[0] = 0.7071067811865475*bmag_inv[7]*nuVtSqSum[7]*m_+0.7071067811865475*bmag_inv[6]*nuVtSqSum[6]*m_+0.7071067811865475*bmag_inv[5]*nuVtSqSum[5]*m_+0.7071067811865475*bmag_inv[4]*nuVtSqSum[4]*m_+0.7071067811865475*bmag_inv[3]*nuVtSqSum[3]*m_+0.7071067811865475*bmag_inv[2]*nuVtSqSum[2]*m_+0.7071067811865475*bmag_inv[1]*nuVtSqSum[1]*m_+0.7071067811865475*bmag_inv[0]*nuVtSqSum[0]*m_; 
  diffFac[1] = 0.7071067811865475*bmag_inv[6]*nuVtSqSum[7]*m_+0.7071067811865475*nuVtSqSum[6]*bmag_inv[7]*m_+0.7071067811865475*bmag_inv[3]*nuVtSqSum[5]*m_+0.7071067811865475*nuVtSqSum[3]*bmag_inv[5]*m_+0.7071067811865475*bmag_inv[2]*nuVtSqSum[4]*m_+0.7071067811865475*nuVtSqSum[2]*bmag_inv[4]*m_+0.7071067811865475*bmag_inv[0]*nuVtSqSum[1]*m_+0.7071067811865475*nuVtSqSum[0]*bmag_inv[1]*m_; 
  diffFac[2] = 0.7071067811865475*bmag_inv[5]*nuVtSqSum[7]*m_+0.7071067811865475*nuVtSqSum[5]*bmag_inv[7]*m_+0.7071067811865475*bmag_inv[3]*nuVtSqSum[6]*m_+0.7071067811865475*nuVtSqSum[3]*bmag_inv[6]*m_+0.7071067811865475*bmag_inv[1]*nuVtSqSum[4]*m_+0.7071067811865475*nuVtSqSum[1]*bmag_inv[4]*m_+0.7071067811865475*bmag_inv[0]*nuVtSqSum[2]*m_+0.7071067811865475*nuVtSqSum[0]*bmag_inv[2]*m_; 
  diffFac[3] = 0.7071067811865475*bmag_inv[4]*nuVtSqSum[7]*m_+0.7071067811865475*nuVtSqSum[4]*bmag_inv[7]*m_+0.7071067811865475*bmag_inv[2]*nuVtSqSum[6]*m_+0.7071067811865475*nuVtSqSum[2]*bmag_inv[6]*m_+0.7071067811865475*bmag_inv[1]*nuVtSqSum[5]*m_+0.7071067811865475*nuVtSqSum[1]*bmag_inv[5]*m_+0.7071067811865475*bmag_inv[0]*nuVtSqSum[3]*m_+0.7071067811865475*nuVtSqSum[0]*bmag_inv[3]*m_; 
  diffFac[4] = 0.7071067811865475*bmag_inv[3]*nuVtSqSum[7]*m_+0.7071067811865475*nuVtSqSum[3]*bmag_inv[7]*m_+0.7071067811865475*bmag_inv[5]*nuVtSqSum[6]*m_+0.7071067811865475*nuVtSqSum[5]*bmag_inv[6]*m_+0.7071067811865475*bmag_inv[0]*nuVtSqSum[4]*m_+0.7071067811865475*nuVtSqSum[0]*bmag_inv[4]*m_+0.7071067811865475*bmag_inv[1]*nuVtSqSum[2]*m_+0.7071067811865475*nuVtSqSum[1]*bmag_inv[2]*m_; 
  diffFac[5] = 0.7071067811865475*bmag_inv[2]*nuVtSqSum[7]*m_+0.7071067811865475*nuVtSqSum[2]*bmag_inv[7]*m_+0.7071067811865475*bmag_inv[4]*nuVtSqSum[6]*m_+0.7071067811865475*nuVtSqSum[4]*bmag_inv[6]*m_+0.7071067811865475*bmag_inv[0]*nuVtSqSum[5]*m_+0.7071067811865475*nuVtSqSum[0]*bmag_inv[5]*m_+0.7071067811865475*bmag_inv[1]*nuVtSqSum[3]*m_+0.7071067811865475*nuVtSqSum[1]*bmag_inv[3]*m_; 
  diffFac[6] = 0.7071067811865475*bmag_inv[1]*nuVtSqSum[7]*m_+0.7071067811865475*nuVtSqSum[1]*bmag_inv[7]*m_+0.7071067811865475*bmag_inv[0]*nuVtSqSum[6]*m_+0.7071067811865475*nuVtSqSum[0]*bmag_inv[6]*m_+0.7071067811865475*bmag_inv[4]*nuVtSqSum[5]*m_+0.7071067811865475*nuVtSqSum[4]*bmag_inv[5]*m_+0.7071067811865475*bmag_inv[2]*nuVtSqSum[3]*m_+0.7071067811865475*nuVtSqSum[2]*bmag_inv[3]*m_; 
  diffFac[7] = 0.7071067811865475*bmag_inv[0]*nuVtSqSum[7]*m_+0.7071067811865475*nuVtSqSum[0]*bmag_inv[7]*m_+0.7071067811865475*bmag_inv[1]*nuVtSqSum[6]*m_+0.7071067811865475*nuVtSqSum[1]*bmag_inv[6]*m_+0.7071067811865475*bmag_inv[2]*nuVtSqSum[5]*m_+0.7071067811865475*nuVtSqSum[2]*bmag_inv[5]*m_+0.7071067811865475*bmag_inv[3]*nuVtSqSum[4]*m_+0.7071067811865475*nuVtSqSum[3]*bmag_inv[4]*m_; 

  diff_incr[0] = 0.3535533905932737*diffFac[7]*temp_diff[16]+0.3535533905932737*diffFac[6]*temp_diff[8]+0.3535533905932737*diffFac[5]*temp_diff[7]+0.3535533905932737*diffFac[4]*temp_diff[6]+0.3535533905932737*diffFac[3]*temp_diff[3]+0.3535533905932737*diffFac[2]*temp_diff[2]+0.3535533905932737*diffFac[1]*temp_diff[1]+0.3535533905932737*diffFac[0]*temp_diff[0]; 
  diff_incr[1] = 0.3535533905932737*diffFac[6]*temp_diff[16]+0.3535533905932737*diffFac[7]*temp_diff[8]+0.3535533905932737*diffFac[3]*temp_diff[7]+0.3535533905932737*diffFac[2]*temp_diff[6]+0.3535533905932737*temp_diff[3]*diffFac[5]+0.3535533905932737*temp_diff[2]*diffFac[4]+0.3535533905932737*diffFac[0]*temp_diff[1]+0.3535533905932737*temp_diff[0]*diffFac[1]; 
  diff_incr[2] = 0.3535533905932737*diffFac[5]*temp_diff[16]+0.3535533905932737*diffFac[3]*temp_diff[8]+0.3535533905932737*diffFac[7]*temp_diff[7]+0.3535533905932737*diffFac[1]*temp_diff[6]+0.3535533905932737*temp_diff[3]*diffFac[6]+0.3535533905932737*temp_diff[1]*diffFac[4]+0.3535533905932737*diffFac[0]*temp_diff[2]+0.3535533905932737*temp_diff[0]*diffFac[2]; 
  diff_incr[3] = 0.3535533905932737*diffFac[4]*temp_diff[16]+0.3535533905932737*diffFac[2]*temp_diff[8]+0.3535533905932737*diffFac[1]*temp_diff[7]+0.3535533905932737*temp_diff[6]*diffFac[7]+0.3535533905932737*temp_diff[2]*diffFac[6]+0.3535533905932737*temp_diff[1]*diffFac[5]+0.3535533905932737*diffFac[0]*temp_diff[3]+0.3535533905932737*temp_diff[0]*diffFac[3]; 
  diff_incr[4] = 0.3535533905932737*diffFac[7]*temp_diff[26]+0.3535533905932737*diffFac[6]*temp_diff[19]+0.3535533905932737*diffFac[5]*temp_diff[18]+0.3535533905932737*diffFac[4]*temp_diff[17]+0.3535533905932737*diffFac[3]*temp_diff[11]+0.3535533905932737*diffFac[2]*temp_diff[10]+0.3535533905932737*diffFac[1]*temp_diff[9]+0.3535533905932737*diffFac[0]*temp_diff[4]; 
  diff_incr[5] = 0.3535533905932737*diffFac[7]*temp_diff[27]+0.3535533905932737*diffFac[6]*temp_diff[22]+0.3535533905932737*diffFac[5]*temp_diff[21]+0.3535533905932737*diffFac[4]*temp_diff[20]+0.3535533905932737*diffFac[3]*temp_diff[14]+0.3535533905932737*diffFac[2]*temp_diff[13]+0.3535533905932737*diffFac[1]*temp_diff[12]+0.3535533905932737*diffFac[0]*temp_diff[5]; 
  diff_incr[6] = 0.3535533905932737*diffFac[3]*temp_diff[16]+0.3535533905932737*diffFac[5]*temp_diff[8]+0.3535533905932737*diffFac[6]*temp_diff[7]+0.3535533905932737*temp_diff[3]*diffFac[7]+0.3535533905932737*diffFac[0]*temp_diff[6]+0.3535533905932737*temp_diff[0]*diffFac[4]+0.3535533905932737*diffFac[1]*temp_diff[2]+0.3535533905932737*temp_diff[1]*diffFac[2]; 
  diff_incr[7] = 0.3535533905932737*diffFac[2]*temp_diff[16]+0.3535533905932737*diffFac[4]*temp_diff[8]+0.3535533905932737*diffFac[0]*temp_diff[7]+0.3535533905932737*temp_diff[2]*diffFac[7]+0.3535533905932737*diffFac[6]*temp_diff[6]+0.3535533905932737*temp_diff[0]*diffFac[5]+0.3535533905932737*diffFac[1]*temp_diff[3]+0.3535533905932737*temp_diff[1]*diffFac[3]; 
  diff_incr[8] = 0.3535533905932737*diffFac[1]*temp_diff[16]+0.3535533905932737*diffFac[0]*temp_diff[8]+0.3535533905932737*diffFac[4]*temp_diff[7]+0.3535533905932737*temp_diff[1]*diffFac[7]+0.3535533905932737*diffFac[5]*temp_diff[6]+0.3535533905932737*temp_diff[0]*diffFac[6]+0.3535533905932737*diffFac[2]*temp_diff[3]+0.3535533905932737*temp_diff[2]*diffFac[3]; 
  diff_incr[9] = 0.3535533905932737*diffFac[6]*temp_diff[26]+0.3535533905932737*diffFac[7]*temp_diff[19]+0.3535533905932737*diffFac[3]*temp_diff[18]+0.3535533905932737*diffFac[2]*temp_diff[17]+0.3535533905932737*diffFac[5]*temp_diff[11]+0.3535533905932737*diffFac[4]*temp_diff[10]+0.3535533905932737*diffFac[0]*temp_diff[9]+0.3535533905932737*diffFac[1]*temp_diff[4]; 
  diff_incr[10] = 0.3535533905932737*diffFac[5]*temp_diff[26]+0.3535533905932737*diffFac[3]*temp_diff[19]+0.3535533905932737*diffFac[7]*temp_diff[18]+0.3535533905932737*diffFac[1]*temp_diff[17]+0.3535533905932737*diffFac[6]*temp_diff[11]+0.3535533905932737*diffFac[0]*temp_diff[10]+0.3535533905932737*diffFac[4]*temp_diff[9]+0.3535533905932737*diffFac[2]*temp_diff[4]; 
  diff_incr[11] = 0.3535533905932737*diffFac[4]*temp_diff[26]+0.3535533905932737*diffFac[2]*temp_diff[19]+0.3535533905932737*diffFac[1]*temp_diff[18]+0.3535533905932737*diffFac[7]*temp_diff[17]+0.3535533905932737*diffFac[0]*temp_diff[11]+0.3535533905932737*diffFac[6]*temp_diff[10]+0.3535533905932737*diffFac[5]*temp_diff[9]+0.3535533905932737*diffFac[3]*temp_diff[4]; 
  diff_incr[12] = 0.3535533905932737*diffFac[6]*temp_diff[27]+0.3535533905932737*diffFac[7]*temp_diff[22]+0.3535533905932737*diffFac[3]*temp_diff[21]+0.3535533905932737*diffFac[2]*temp_diff[20]+0.3535533905932737*diffFac[5]*temp_diff[14]+0.3535533905932737*diffFac[4]*temp_diff[13]+0.3535533905932737*diffFac[0]*temp_diff[12]+0.3535533905932737*diffFac[1]*temp_diff[5]; 
  diff_incr[13] = 0.3535533905932737*diffFac[5]*temp_diff[27]+0.3535533905932737*diffFac[3]*temp_diff[22]+0.3535533905932737*diffFac[7]*temp_diff[21]+0.3535533905932737*diffFac[1]*temp_diff[20]+0.3535533905932737*diffFac[6]*temp_diff[14]+0.3535533905932737*diffFac[0]*temp_diff[13]+0.3535533905932737*diffFac[4]*temp_diff[12]+0.3535533905932737*diffFac[2]*temp_diff[5]; 
  diff_incr[14] = 0.3535533905932737*diffFac[4]*temp_diff[27]+0.3535533905932737*diffFac[2]*temp_diff[22]+0.3535533905932737*diffFac[1]*temp_diff[21]+0.3535533905932737*diffFac[7]*temp_diff[20]+0.3535533905932737*diffFac[0]*temp_diff[14]+0.3535533905932737*diffFac[6]*temp_diff[13]+0.3535533905932737*diffFac[5]*temp_diff[12]+0.3535533905932737*diffFac[3]*temp_diff[5]; 
  diff_incr[15] = 0.3535533905932737*diffFac[7]*temp_diff[31]+0.3535533905932737*diffFac[6]*temp_diff[30]+0.3535533905932737*diffFac[5]*temp_diff[29]+0.3535533905932737*diffFac[4]*temp_diff[28]+0.3535533905932737*diffFac[3]*temp_diff[25]+0.3535533905932737*diffFac[2]*temp_diff[24]+0.3535533905932737*diffFac[1]*temp_diff[23]+0.3535533905932737*diffFac[0]*temp_diff[15]; 
  diff_incr[16] = 0.3535533905932737*diffFac[0]*temp_diff[16]+0.3535533905932737*diffFac[1]*temp_diff[8]+0.3535533905932737*diffFac[2]*temp_diff[7]+0.3535533905932737*temp_diff[0]*diffFac[7]+0.3535533905932737*diffFac[3]*temp_diff[6]+0.3535533905932737*temp_diff[1]*diffFac[6]+0.3535533905932737*temp_diff[2]*diffFac[5]+0.3535533905932737*temp_diff[3]*diffFac[4]; 
  diff_incr[17] = 0.3535533905932737*diffFac[3]*temp_diff[26]+0.3535533905932737*diffFac[5]*temp_diff[19]+0.3535533905932737*diffFac[6]*temp_diff[18]+0.3535533905932737*diffFac[0]*temp_diff[17]+0.3535533905932737*diffFac[7]*temp_diff[11]+0.3535533905932737*diffFac[1]*temp_diff[10]+0.3535533905932737*diffFac[2]*temp_diff[9]+0.3535533905932737*diffFac[4]*temp_diff[4]; 
  diff_incr[18] = 0.3535533905932737*diffFac[2]*temp_diff[26]+0.3535533905932737*diffFac[4]*temp_diff[19]+0.3535533905932737*diffFac[0]*temp_diff[18]+0.3535533905932737*diffFac[6]*temp_diff[17]+0.3535533905932737*diffFac[1]*temp_diff[11]+0.3535533905932737*diffFac[7]*temp_diff[10]+0.3535533905932737*diffFac[3]*temp_diff[9]+0.3535533905932737*temp_diff[4]*diffFac[5]; 
  diff_incr[19] = 0.3535533905932737*diffFac[1]*temp_diff[26]+0.3535533905932737*diffFac[0]*temp_diff[19]+0.3535533905932737*diffFac[4]*temp_diff[18]+0.3535533905932737*diffFac[5]*temp_diff[17]+0.3535533905932737*diffFac[2]*temp_diff[11]+0.3535533905932737*diffFac[3]*temp_diff[10]+0.3535533905932737*diffFac[7]*temp_diff[9]+0.3535533905932737*temp_diff[4]*diffFac[6]; 
  diff_incr[20] = 0.3535533905932737*diffFac[3]*temp_diff[27]+0.3535533905932737*diffFac[5]*temp_diff[22]+0.3535533905932737*diffFac[6]*temp_diff[21]+0.3535533905932737*diffFac[0]*temp_diff[20]+0.3535533905932737*diffFac[7]*temp_diff[14]+0.3535533905932737*diffFac[1]*temp_diff[13]+0.3535533905932737*diffFac[2]*temp_diff[12]+0.3535533905932737*diffFac[4]*temp_diff[5]; 
  diff_incr[21] = 0.3535533905932737*diffFac[2]*temp_diff[27]+0.3535533905932737*diffFac[4]*temp_diff[22]+0.3535533905932737*diffFac[0]*temp_diff[21]+0.3535533905932737*diffFac[6]*temp_diff[20]+0.3535533905932737*diffFac[1]*temp_diff[14]+0.3535533905932737*diffFac[7]*temp_diff[13]+0.3535533905932737*diffFac[3]*temp_diff[12]+0.3535533905932737*diffFac[5]*temp_diff[5]; 
  diff_incr[22] = 0.3535533905932737*diffFac[1]*temp_diff[27]+0.3535533905932737*diffFac[0]*temp_diff[22]+0.3535533905932737*diffFac[4]*temp_diff[21]+0.3535533905932737*diffFac[5]*temp_diff[20]+0.3535533905932737*diffFac[2]*temp_diff[14]+0.3535533905932737*diffFac[3]*temp_diff[13]+0.3535533905932737*diffFac[7]*temp_diff[12]+0.3535533905932737*temp_diff[5]*diffFac[6]; 
  diff_incr[23] = 0.3535533905932737*diffFac[6]*temp_diff[31]+0.3535533905932737*diffFac[7]*temp_diff[30]+0.3535533905932737*diffFac[3]*temp_diff[29]+0.3535533905932737*diffFac[2]*temp_diff[28]+0.3535533905932737*diffFac[5]*temp_diff[25]+0.3535533905932737*diffFac[4]*temp_diff[24]+0.3535533905932737*diffFac[0]*temp_diff[23]+0.3535533905932737*diffFac[1]*temp_diff[15]; 
  diff_incr[24] = 0.3535533905932737*diffFac[5]*temp_diff[31]+0.3535533905932737*diffFac[3]*temp_diff[30]+0.3535533905932737*diffFac[7]*temp_diff[29]+0.3535533905932737*diffFac[1]*temp_diff[28]+0.3535533905932737*diffFac[6]*temp_diff[25]+0.3535533905932737*diffFac[0]*temp_diff[24]+0.3535533905932737*diffFac[4]*temp_diff[23]+0.3535533905932737*diffFac[2]*temp_diff[15]; 
  diff_incr[25] = 0.3535533905932737*diffFac[4]*temp_diff[31]+0.3535533905932737*diffFac[2]*temp_diff[30]+0.3535533905932737*diffFac[1]*temp_diff[29]+0.3535533905932737*diffFac[7]*temp_diff[28]+0.3535533905932737*diffFac[0]*temp_diff[25]+0.3535533905932737*diffFac[6]*temp_diff[24]+0.3535533905932737*diffFac[5]*temp_diff[23]+0.3535533905932737*diffFac[3]*temp_diff[15]; 
  diff_incr[26] = 0.3535533905932737*diffFac[0]*temp_diff[26]+0.3535533905932737*diffFac[1]*temp_diff[19]+0.3535533905932737*diffFac[2]*temp_diff[18]+0.3535533905932737*diffFac[3]*temp_diff[17]+0.3535533905932737*diffFac[4]*temp_diff[11]+0.3535533905932737*diffFac[5]*temp_diff[10]+0.3535533905932737*diffFac[6]*temp_diff[9]+0.3535533905932737*temp_diff[4]*diffFac[7]; 
  diff_incr[27] = 0.3535533905932737*diffFac[0]*temp_diff[27]+0.3535533905932737*diffFac[1]*temp_diff[22]+0.3535533905932737*diffFac[2]*temp_diff[21]+0.3535533905932737*diffFac[3]*temp_diff[20]+0.3535533905932737*diffFac[4]*temp_diff[14]+0.3535533905932737*diffFac[5]*temp_diff[13]+0.3535533905932737*diffFac[6]*temp_diff[12]+0.3535533905932737*temp_diff[5]*diffFac[7]; 
  diff_incr[28] = 0.3535533905932737*diffFac[3]*temp_diff[31]+0.3535533905932737*diffFac[5]*temp_diff[30]+0.3535533905932737*diffFac[6]*temp_diff[29]+0.3535533905932737*diffFac[0]*temp_diff[28]+0.3535533905932737*diffFac[7]*temp_diff[25]+0.3535533905932737*diffFac[1]*temp_diff[24]+0.3535533905932737*diffFac[2]*temp_diff[23]+0.3535533905932737*diffFac[4]*temp_diff[15]; 
  diff_incr[29] = 0.3535533905932737*diffFac[2]*temp_diff[31]+0.3535533905932737*diffFac[4]*temp_diff[30]+0.3535533905932737*diffFac[0]*temp_diff[29]+0.3535533905932737*diffFac[6]*temp_diff[28]+0.3535533905932737*diffFac[1]*temp_diff[25]+0.3535533905932737*diffFac[7]*temp_diff[24]+0.3535533905932737*diffFac[3]*temp_diff[23]+0.3535533905932737*diffFac[5]*temp_diff[15]; 
  diff_incr[30] = 0.3535533905932737*diffFac[1]*temp_diff[31]+0.3535533905932737*diffFac[0]*temp_diff[30]+0.3535533905932737*diffFac[4]*temp_diff[29]+0.3535533905932737*diffFac[5]*temp_diff[28]+0.3535533905932737*diffFac[2]*temp_diff[25]+0.3535533905932737*diffFac[3]*temp_diff[24]+0.3535533905932737*diffFac[7]*temp_diff[23]+0.3535533905932737*diffFac[6]*temp_diff[15]; 
  diff_incr[31] = 0.3535533905932737*diffFac[0]*temp_diff[31]+0.3535533905932737*diffFac[1]*temp_diff[30]+0.3535533905932737*diffFac[2]*temp_diff[29]+0.3535533905932737*diffFac[3]*temp_diff[28]+0.3535533905932737*diffFac[4]*temp_diff[25]+0.3535533905932737*diffFac[5]*temp_diff[24]+0.3535533905932737*diffFac[6]*temp_diff[23]+0.3535533905932737*diffFac[7]*temp_diff[15]; 
  diff_incr[32] = 0.3535533905932737*diffFac[7]*temp_diff[43]+0.3535533905932737*diffFac[6]*temp_diff[39]+0.3535533905932737*diffFac[5]*temp_diff[38]+0.3535533905932737*diffFac[4]*temp_diff[37]+0.3535533905932737*diffFac[3]*temp_diff[35]+0.3535533905932737*diffFac[2]*temp_diff[34]+0.3535533905932737*diffFac[1]*temp_diff[33]+0.3535533905932737*diffFac[0]*temp_diff[32]; 
  diff_incr[33] = 0.3535533905932737*diffFac[6]*temp_diff[43]+0.3535533905932737*diffFac[7]*temp_diff[39]+0.3535533905932737*diffFac[3]*temp_diff[38]+0.3535533905932737*diffFac[2]*temp_diff[37]+0.3535533905932737*diffFac[5]*temp_diff[35]+0.3535533905932737*diffFac[4]*temp_diff[34]+0.3535533905932737*diffFac[0]*temp_diff[33]+0.3535533905932737*diffFac[1]*temp_diff[32]; 
  diff_incr[34] = 0.3535533905932737*diffFac[5]*temp_diff[43]+0.3535533905932737*diffFac[3]*temp_diff[39]+0.3535533905932737*diffFac[7]*temp_diff[38]+0.3535533905932737*diffFac[1]*temp_diff[37]+0.3535533905932737*diffFac[6]*temp_diff[35]+0.3535533905932737*diffFac[0]*temp_diff[34]+0.3535533905932737*diffFac[4]*temp_diff[33]+0.3535533905932737*diffFac[2]*temp_diff[32]; 
  diff_incr[35] = 0.3535533905932737*diffFac[4]*temp_diff[43]+0.3535533905932737*diffFac[2]*temp_diff[39]+0.3535533905932737*diffFac[1]*temp_diff[38]+0.3535533905932737*diffFac[7]*temp_diff[37]+0.3535533905932737*diffFac[0]*temp_diff[35]+0.3535533905932737*diffFac[6]*temp_diff[34]+0.3535533905932737*diffFac[5]*temp_diff[33]+0.3535533905932737*diffFac[3]*temp_diff[32]; 
  diff_incr[36] = 0.3535533905932737*diffFac[7]*temp_diff[47]+0.3535533905932737*diffFac[6]*temp_diff[46]+0.3535533905932737*diffFac[5]*temp_diff[45]+0.3535533905932737*diffFac[4]*temp_diff[44]+0.3535533905932737*diffFac[3]*temp_diff[42]+0.3535533905932737*diffFac[2]*temp_diff[41]+0.3535533905932737*diffFac[1]*temp_diff[40]+0.3535533905932737*diffFac[0]*temp_diff[36]; 
  diff_incr[37] = 0.3535533905932737*diffFac[3]*temp_diff[43]+0.3535533905932737*diffFac[5]*temp_diff[39]+0.3535533905932737*diffFac[6]*temp_diff[38]+0.3535533905932737*diffFac[0]*temp_diff[37]+0.3535533905932737*diffFac[7]*temp_diff[35]+0.3535533905932737*diffFac[1]*temp_diff[34]+0.3535533905932737*diffFac[2]*temp_diff[33]+0.3535533905932737*diffFac[4]*temp_diff[32]; 
  diff_incr[38] = 0.3535533905932737*diffFac[2]*temp_diff[43]+0.3535533905932737*diffFac[4]*temp_diff[39]+0.3535533905932737*diffFac[0]*temp_diff[38]+0.3535533905932737*diffFac[6]*temp_diff[37]+0.3535533905932737*diffFac[1]*temp_diff[35]+0.3535533905932737*diffFac[7]*temp_diff[34]+0.3535533905932737*diffFac[3]*temp_diff[33]+0.3535533905932737*diffFac[5]*temp_diff[32]; 
  diff_incr[39] = 0.3535533905932737*diffFac[1]*temp_diff[43]+0.3535533905932737*diffFac[0]*temp_diff[39]+0.3535533905932737*diffFac[4]*temp_diff[38]+0.3535533905932737*diffFac[5]*temp_diff[37]+0.3535533905932737*diffFac[2]*temp_diff[35]+0.3535533905932737*diffFac[3]*temp_diff[34]+0.3535533905932737*diffFac[7]*temp_diff[33]+0.3535533905932737*diffFac[6]*temp_diff[32]; 
  diff_incr[40] = 0.3535533905932737*diffFac[6]*temp_diff[47]+0.3535533905932737*diffFac[7]*temp_diff[46]+0.3535533905932737*diffFac[3]*temp_diff[45]+0.3535533905932737*diffFac[2]*temp_diff[44]+0.3535533905932737*diffFac[5]*temp_diff[42]+0.3535533905932737*diffFac[4]*temp_diff[41]+0.3535533905932737*diffFac[0]*temp_diff[40]+0.3535533905932737*diffFac[1]*temp_diff[36]; 
  diff_incr[41] = 0.3535533905932737*diffFac[5]*temp_diff[47]+0.3535533905932737*diffFac[3]*temp_diff[46]+0.3535533905932737*diffFac[7]*temp_diff[45]+0.3535533905932737*diffFac[1]*temp_diff[44]+0.3535533905932737*diffFac[6]*temp_diff[42]+0.3535533905932737*diffFac[0]*temp_diff[41]+0.3535533905932737*diffFac[4]*temp_diff[40]+0.3535533905932737*diffFac[2]*temp_diff[36]; 
  diff_incr[42] = 0.3535533905932737*diffFac[4]*temp_diff[47]+0.3535533905932737*diffFac[2]*temp_diff[46]+0.3535533905932737*diffFac[1]*temp_diff[45]+0.3535533905932737*diffFac[7]*temp_diff[44]+0.3535533905932737*diffFac[0]*temp_diff[42]+0.3535533905932737*diffFac[6]*temp_diff[41]+0.3535533905932737*diffFac[5]*temp_diff[40]+0.3535533905932737*diffFac[3]*temp_diff[36]; 
  diff_incr[43] = 0.3535533905932737*diffFac[0]*temp_diff[43]+0.3535533905932737*diffFac[1]*temp_diff[39]+0.3535533905932737*diffFac[2]*temp_diff[38]+0.3535533905932737*diffFac[3]*temp_diff[37]+0.3535533905932737*diffFac[4]*temp_diff[35]+0.3535533905932737*diffFac[5]*temp_diff[34]+0.3535533905932737*diffFac[6]*temp_diff[33]+0.3535533905932737*diffFac[7]*temp_diff[32]; 
  diff_incr[44] = 0.3535533905932737*diffFac[3]*temp_diff[47]+0.3535533905932737*diffFac[5]*temp_diff[46]+0.3535533905932737*diffFac[6]*temp_diff[45]+0.3535533905932737*diffFac[0]*temp_diff[44]+0.3535533905932737*diffFac[7]*temp_diff[42]+0.3535533905932737*diffFac[1]*temp_diff[41]+0.3535533905932737*diffFac[2]*temp_diff[40]+0.3535533905932737*diffFac[4]*temp_diff[36]; 
  diff_incr[45] = 0.3535533905932737*diffFac[2]*temp_diff[47]+0.3535533905932737*diffFac[4]*temp_diff[46]+0.3535533905932737*diffFac[0]*temp_diff[45]+0.3535533905932737*diffFac[6]*temp_diff[44]+0.3535533905932737*diffFac[1]*temp_diff[42]+0.3535533905932737*diffFac[7]*temp_diff[41]+0.3535533905932737*diffFac[3]*temp_diff[40]+0.3535533905932737*diffFac[5]*temp_diff[36]; 
  diff_incr[46] = 0.3535533905932737*diffFac[1]*temp_diff[47]+0.3535533905932737*diffFac[0]*temp_diff[46]+0.3535533905932737*diffFac[4]*temp_diff[45]+0.3535533905932737*diffFac[5]*temp_diff[44]+0.3535533905932737*diffFac[2]*temp_diff[42]+0.3535533905932737*diffFac[3]*temp_diff[41]+0.3535533905932737*diffFac[7]*temp_diff[40]+0.3535533905932737*diffFac[6]*temp_diff[36]; 
  diff_incr[47] = 0.3535533905932737*diffFac[0]*temp_diff[47]+0.3535533905932737*diffFac[1]*temp_diff[46]+0.3535533905932737*diffFac[2]*temp_diff[45]+0.3535533905932737*diffFac[3]*temp_diff[44]+0.3535533905932737*diffFac[4]*temp_diff[42]+0.3535533905932737*diffFac[5]*temp_diff[41]+0.3535533905932737*diffFac[6]*temp_diff[40]+0.3535533905932737*diffFac[7]*temp_diff[36]; 

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
