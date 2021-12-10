#include <gkyl_vlasov_lbo_kernels.h> 
GKYL_CU_DH void vlasov_lbo_diff_surfvx_2x3v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[5]:         cell-center coordinates. 
  // dxv[5]:       cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[12]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[4]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      distribution function in cells 
  // out:           incremented distribution function in cell 
  double rdvSq4 = 4.0/(dxv[2]*dxv[2]); 
  double diff_incr[32] = {0.0}; 

  diff_incr[0] = (-0.270632938682637*nuVtSqSum[3]*fr[16])+0.270632938682637*nuVtSqSum[3]*fl[16]-0.270632938682637*nuVtSqSum[2]*fr[8]+0.270632938682637*nuVtSqSum[2]*fl[8]-0.270632938682637*nuVtSqSum[1]*fr[7]+0.270632938682637*nuVtSqSum[1]*fl[7]+0.28125*nuVtSqSum[3]*fr[6]+0.28125*nuVtSqSum[3]*fl[6]-0.5625*nuVtSqSum[3]*fc[6]-0.270632938682637*nuVtSqSum[0]*fr[3]+0.270632938682637*nuVtSqSum[0]*fl[3]+0.28125*fr[2]*nuVtSqSum[2]+0.28125*fl[2]*nuVtSqSum[2]-0.5625*fc[2]*nuVtSqSum[2]+0.28125*fr[1]*nuVtSqSum[1]+0.28125*fl[1]*nuVtSqSum[1]-0.5625*fc[1]*nuVtSqSum[1]+0.28125*fr[0]*nuVtSqSum[0]+0.28125*fl[0]*nuVtSqSum[0]-0.5625*fc[0]*nuVtSqSum[0]; 
  diff_incr[1] = (-0.270632938682637*nuVtSqSum[2]*fr[16])+0.270632938682637*nuVtSqSum[2]*fl[16]-0.270632938682637*nuVtSqSum[3]*fr[8]+0.270632938682637*nuVtSqSum[3]*fl[8]-0.270632938682637*nuVtSqSum[0]*fr[7]+0.270632938682637*nuVtSqSum[0]*fl[7]+0.28125*nuVtSqSum[2]*fr[6]+0.28125*nuVtSqSum[2]*fl[6]-0.5625*nuVtSqSum[2]*fc[6]+0.28125*fr[2]*nuVtSqSum[3]+0.28125*fl[2]*nuVtSqSum[3]-0.5625*fc[2]*nuVtSqSum[3]-0.270632938682637*nuVtSqSum[1]*fr[3]+0.270632938682637*nuVtSqSum[1]*fl[3]+0.28125*fr[0]*nuVtSqSum[1]+0.28125*fl[0]*nuVtSqSum[1]-0.5625*fc[0]*nuVtSqSum[1]+0.28125*nuVtSqSum[0]*fr[1]+0.28125*nuVtSqSum[0]*fl[1]-0.5625*nuVtSqSum[0]*fc[1]; 
  diff_incr[2] = (-0.270632938682637*nuVtSqSum[1]*fr[16])+0.270632938682637*nuVtSqSum[1]*fl[16]-0.270632938682637*nuVtSqSum[0]*fr[8]+0.270632938682637*nuVtSqSum[0]*fl[8]-0.270632938682637*nuVtSqSum[3]*fr[7]+0.270632938682637*nuVtSqSum[3]*fl[7]+0.28125*nuVtSqSum[1]*fr[6]+0.28125*nuVtSqSum[1]*fl[6]-0.5625*nuVtSqSum[1]*fc[6]+0.28125*fr[1]*nuVtSqSum[3]+0.28125*fl[1]*nuVtSqSum[3]-0.5625*fc[1]*nuVtSqSum[3]-0.270632938682637*nuVtSqSum[2]*fr[3]+0.270632938682637*nuVtSqSum[2]*fl[3]+0.28125*fr[0]*nuVtSqSum[2]+0.28125*fl[0]*nuVtSqSum[2]-0.5625*fc[0]*nuVtSqSum[2]+0.28125*nuVtSqSum[0]*fr[2]+0.28125*nuVtSqSum[0]*fl[2]-0.5625*nuVtSqSum[0]*fc[2]; 
  diff_incr[3] = (-0.21875*nuVtSqSum[3]*fr[16])-0.21875*nuVtSqSum[3]*fl[16]-1.4375*nuVtSqSum[3]*fc[16]-0.21875*nuVtSqSum[2]*fr[8]-0.21875*nuVtSqSum[2]*fl[8]-1.4375*nuVtSqSum[2]*fc[8]-0.21875*nuVtSqSum[1]*fr[7]-0.21875*nuVtSqSum[1]*fl[7]-1.4375*nuVtSqSum[1]*fc[7]+0.270632938682637*nuVtSqSum[3]*fr[6]-0.270632938682637*nuVtSqSum[3]*fl[6]-0.21875*nuVtSqSum[0]*fr[3]-0.21875*nuVtSqSum[0]*fl[3]-1.4375*nuVtSqSum[0]*fc[3]+0.270632938682637*fr[2]*nuVtSqSum[2]-0.270632938682637*fl[2]*nuVtSqSum[2]+0.270632938682637*fr[1]*nuVtSqSum[1]-0.270632938682637*fl[1]*nuVtSqSum[1]+0.270632938682637*fr[0]*nuVtSqSum[0]-0.270632938682637*fl[0]*nuVtSqSum[0]; 
  diff_incr[4] = (-0.270632938682637*nuVtSqSum[3]*fr[26])+0.270632938682637*nuVtSqSum[3]*fl[26]-0.270632938682637*nuVtSqSum[2]*fr[19]+0.270632938682637*nuVtSqSum[2]*fl[19]-0.270632938682637*nuVtSqSum[1]*fr[18]+0.270632938682637*nuVtSqSum[1]*fl[18]+0.28125*nuVtSqSum[3]*fr[17]+0.28125*nuVtSqSum[3]*fl[17]-0.5625*nuVtSqSum[3]*fc[17]-0.270632938682637*nuVtSqSum[0]*fr[11]+0.270632938682637*nuVtSqSum[0]*fl[11]+0.28125*nuVtSqSum[2]*fr[10]+0.28125*nuVtSqSum[2]*fl[10]-0.5625*nuVtSqSum[2]*fc[10]+0.28125*nuVtSqSum[1]*fr[9]+0.28125*nuVtSqSum[1]*fl[9]-0.5625*nuVtSqSum[1]*fc[9]+0.28125*nuVtSqSum[0]*fr[4]+0.28125*nuVtSqSum[0]*fl[4]-0.5625*nuVtSqSum[0]*fc[4]; 
  diff_incr[5] = (-0.270632938682637*nuVtSqSum[3]*fr[27])+0.270632938682637*nuVtSqSum[3]*fl[27]-0.270632938682637*nuVtSqSum[2]*fr[22]+0.270632938682637*nuVtSqSum[2]*fl[22]-0.270632938682637*nuVtSqSum[1]*fr[21]+0.270632938682637*nuVtSqSum[1]*fl[21]+0.28125*nuVtSqSum[3]*fr[20]+0.28125*nuVtSqSum[3]*fl[20]-0.5625*nuVtSqSum[3]*fc[20]-0.270632938682637*nuVtSqSum[0]*fr[14]+0.270632938682637*nuVtSqSum[0]*fl[14]+0.28125*nuVtSqSum[2]*fr[13]+0.28125*nuVtSqSum[2]*fl[13]-0.5625*nuVtSqSum[2]*fc[13]+0.28125*nuVtSqSum[1]*fr[12]+0.28125*nuVtSqSum[1]*fl[12]-0.5625*nuVtSqSum[1]*fc[12]+0.28125*nuVtSqSum[0]*fr[5]+0.28125*nuVtSqSum[0]*fl[5]-0.5625*nuVtSqSum[0]*fc[5]; 
  diff_incr[6] = (-0.270632938682637*nuVtSqSum[0]*fr[16])+0.270632938682637*nuVtSqSum[0]*fl[16]-0.270632938682637*nuVtSqSum[1]*fr[8]+0.270632938682637*nuVtSqSum[1]*fl[8]-0.270632938682637*nuVtSqSum[2]*fr[7]+0.270632938682637*nuVtSqSum[2]*fl[7]+0.28125*nuVtSqSum[0]*fr[6]+0.28125*nuVtSqSum[0]*fl[6]-0.5625*nuVtSqSum[0]*fc[6]-0.270632938682637*fr[3]*nuVtSqSum[3]+0.270632938682637*fl[3]*nuVtSqSum[3]+0.28125*fr[0]*nuVtSqSum[3]+0.28125*fl[0]*nuVtSqSum[3]-0.5625*fc[0]*nuVtSqSum[3]+0.28125*fr[1]*nuVtSqSum[2]+0.28125*fl[1]*nuVtSqSum[2]-0.5625*fc[1]*nuVtSqSum[2]+0.28125*nuVtSqSum[1]*fr[2]+0.28125*nuVtSqSum[1]*fl[2]-0.5625*nuVtSqSum[1]*fc[2]; 
  diff_incr[7] = (-0.21875*nuVtSqSum[2]*fr[16])-0.21875*nuVtSqSum[2]*fl[16]-1.4375*nuVtSqSum[2]*fc[16]-0.21875*nuVtSqSum[3]*fr[8]-0.21875*nuVtSqSum[3]*fl[8]-1.4375*nuVtSqSum[3]*fc[8]-0.21875*nuVtSqSum[0]*fr[7]-0.21875*nuVtSqSum[0]*fl[7]-1.4375*nuVtSqSum[0]*fc[7]+0.270632938682637*nuVtSqSum[2]*fr[6]-0.270632938682637*nuVtSqSum[2]*fl[6]+0.270632938682637*fr[2]*nuVtSqSum[3]-0.270632938682637*fl[2]*nuVtSqSum[3]-0.21875*nuVtSqSum[1]*fr[3]-0.21875*nuVtSqSum[1]*fl[3]-1.4375*nuVtSqSum[1]*fc[3]+0.270632938682637*fr[0]*nuVtSqSum[1]-0.270632938682637*fl[0]*nuVtSqSum[1]+0.270632938682637*nuVtSqSum[0]*fr[1]-0.270632938682637*nuVtSqSum[0]*fl[1]; 
  diff_incr[8] = (-0.21875*nuVtSqSum[1]*fr[16])-0.21875*nuVtSqSum[1]*fl[16]-1.4375*nuVtSqSum[1]*fc[16]-0.21875*nuVtSqSum[0]*fr[8]-0.21875*nuVtSqSum[0]*fl[8]-1.4375*nuVtSqSum[0]*fc[8]-0.21875*nuVtSqSum[3]*fr[7]-0.21875*nuVtSqSum[3]*fl[7]-1.4375*nuVtSqSum[3]*fc[7]+0.270632938682637*nuVtSqSum[1]*fr[6]-0.270632938682637*nuVtSqSum[1]*fl[6]+0.270632938682637*fr[1]*nuVtSqSum[3]-0.270632938682637*fl[1]*nuVtSqSum[3]-0.21875*nuVtSqSum[2]*fr[3]-0.21875*nuVtSqSum[2]*fl[3]-1.4375*nuVtSqSum[2]*fc[3]+0.270632938682637*fr[0]*nuVtSqSum[2]-0.270632938682637*fl[0]*nuVtSqSum[2]+0.270632938682637*nuVtSqSum[0]*fr[2]-0.270632938682637*nuVtSqSum[0]*fl[2]; 
  diff_incr[9] = (-0.270632938682637*nuVtSqSum[2]*fr[26])+0.270632938682637*nuVtSqSum[2]*fl[26]-0.270632938682637*nuVtSqSum[3]*fr[19]+0.270632938682637*nuVtSqSum[3]*fl[19]-0.270632938682637*nuVtSqSum[0]*fr[18]+0.270632938682637*nuVtSqSum[0]*fl[18]+0.28125*nuVtSqSum[2]*fr[17]+0.28125*nuVtSqSum[2]*fl[17]-0.5625*nuVtSqSum[2]*fc[17]-0.270632938682637*nuVtSqSum[1]*fr[11]+0.270632938682637*nuVtSqSum[1]*fl[11]+0.28125*nuVtSqSum[3]*fr[10]+0.28125*nuVtSqSum[3]*fl[10]-0.5625*nuVtSqSum[3]*fc[10]+0.28125*nuVtSqSum[0]*fr[9]+0.28125*nuVtSqSum[0]*fl[9]-0.5625*nuVtSqSum[0]*fc[9]+0.28125*nuVtSqSum[1]*fr[4]+0.28125*nuVtSqSum[1]*fl[4]-0.5625*nuVtSqSum[1]*fc[4]; 
  diff_incr[10] = (-0.270632938682637*nuVtSqSum[1]*fr[26])+0.270632938682637*nuVtSqSum[1]*fl[26]-0.270632938682637*nuVtSqSum[0]*fr[19]+0.270632938682637*nuVtSqSum[0]*fl[19]-0.270632938682637*nuVtSqSum[3]*fr[18]+0.270632938682637*nuVtSqSum[3]*fl[18]+0.28125*nuVtSqSum[1]*fr[17]+0.28125*nuVtSqSum[1]*fl[17]-0.5625*nuVtSqSum[1]*fc[17]-0.270632938682637*nuVtSqSum[2]*fr[11]+0.270632938682637*nuVtSqSum[2]*fl[11]+0.28125*nuVtSqSum[0]*fr[10]+0.28125*nuVtSqSum[0]*fl[10]-0.5625*nuVtSqSum[0]*fc[10]+0.28125*nuVtSqSum[3]*fr[9]+0.28125*nuVtSqSum[3]*fl[9]-0.5625*nuVtSqSum[3]*fc[9]+0.28125*nuVtSqSum[2]*fr[4]+0.28125*nuVtSqSum[2]*fl[4]-0.5625*nuVtSqSum[2]*fc[4]; 
  diff_incr[11] = (-0.21875*nuVtSqSum[3]*fr[26])-0.21875*nuVtSqSum[3]*fl[26]-1.4375*nuVtSqSum[3]*fc[26]-0.21875*nuVtSqSum[2]*fr[19]-0.21875*nuVtSqSum[2]*fl[19]-1.4375*nuVtSqSum[2]*fc[19]-0.21875*nuVtSqSum[1]*fr[18]-0.21875*nuVtSqSum[1]*fl[18]-1.4375*nuVtSqSum[1]*fc[18]+0.270632938682637*nuVtSqSum[3]*fr[17]-0.270632938682637*nuVtSqSum[3]*fl[17]-0.21875*nuVtSqSum[0]*fr[11]-0.21875*nuVtSqSum[0]*fl[11]-1.4375*nuVtSqSum[0]*fc[11]+0.270632938682637*nuVtSqSum[2]*fr[10]-0.270632938682637*nuVtSqSum[2]*fl[10]+0.270632938682637*nuVtSqSum[1]*fr[9]-0.270632938682637*nuVtSqSum[1]*fl[9]+0.270632938682637*nuVtSqSum[0]*fr[4]-0.270632938682637*nuVtSqSum[0]*fl[4]; 
  diff_incr[12] = (-0.270632938682637*nuVtSqSum[2]*fr[27])+0.270632938682637*nuVtSqSum[2]*fl[27]-0.270632938682637*nuVtSqSum[3]*fr[22]+0.270632938682637*nuVtSqSum[3]*fl[22]-0.270632938682637*nuVtSqSum[0]*fr[21]+0.270632938682637*nuVtSqSum[0]*fl[21]+0.28125*nuVtSqSum[2]*fr[20]+0.28125*nuVtSqSum[2]*fl[20]-0.5625*nuVtSqSum[2]*fc[20]-0.270632938682637*nuVtSqSum[1]*fr[14]+0.270632938682637*nuVtSqSum[1]*fl[14]+0.28125*nuVtSqSum[3]*fr[13]+0.28125*nuVtSqSum[3]*fl[13]-0.5625*nuVtSqSum[3]*fc[13]+0.28125*nuVtSqSum[0]*fr[12]+0.28125*nuVtSqSum[0]*fl[12]-0.5625*nuVtSqSum[0]*fc[12]+0.28125*nuVtSqSum[1]*fr[5]+0.28125*nuVtSqSum[1]*fl[5]-0.5625*nuVtSqSum[1]*fc[5]; 
  diff_incr[13] = (-0.270632938682637*nuVtSqSum[1]*fr[27])+0.270632938682637*nuVtSqSum[1]*fl[27]-0.270632938682637*nuVtSqSum[0]*fr[22]+0.270632938682637*nuVtSqSum[0]*fl[22]-0.270632938682637*nuVtSqSum[3]*fr[21]+0.270632938682637*nuVtSqSum[3]*fl[21]+0.28125*nuVtSqSum[1]*fr[20]+0.28125*nuVtSqSum[1]*fl[20]-0.5625*nuVtSqSum[1]*fc[20]-0.270632938682637*nuVtSqSum[2]*fr[14]+0.270632938682637*nuVtSqSum[2]*fl[14]+0.28125*nuVtSqSum[0]*fr[13]+0.28125*nuVtSqSum[0]*fl[13]-0.5625*nuVtSqSum[0]*fc[13]+0.28125*nuVtSqSum[3]*fr[12]+0.28125*nuVtSqSum[3]*fl[12]-0.5625*nuVtSqSum[3]*fc[12]+0.28125*nuVtSqSum[2]*fr[5]+0.28125*nuVtSqSum[2]*fl[5]-0.5625*nuVtSqSum[2]*fc[5]; 
  diff_incr[14] = (-0.21875*nuVtSqSum[3]*fr[27])-0.21875*nuVtSqSum[3]*fl[27]-1.4375*nuVtSqSum[3]*fc[27]-0.21875*nuVtSqSum[2]*fr[22]-0.21875*nuVtSqSum[2]*fl[22]-1.4375*nuVtSqSum[2]*fc[22]-0.21875*nuVtSqSum[1]*fr[21]-0.21875*nuVtSqSum[1]*fl[21]-1.4375*nuVtSqSum[1]*fc[21]+0.270632938682637*nuVtSqSum[3]*fr[20]-0.270632938682637*nuVtSqSum[3]*fl[20]-0.21875*nuVtSqSum[0]*fr[14]-0.21875*nuVtSqSum[0]*fl[14]-1.4375*nuVtSqSum[0]*fc[14]+0.270632938682637*nuVtSqSum[2]*fr[13]-0.270632938682637*nuVtSqSum[2]*fl[13]+0.270632938682637*nuVtSqSum[1]*fr[12]-0.270632938682637*nuVtSqSum[1]*fl[12]+0.270632938682637*nuVtSqSum[0]*fr[5]-0.270632938682637*nuVtSqSum[0]*fl[5]; 
  diff_incr[15] = (-0.270632938682637*nuVtSqSum[3]*fr[31])+0.270632938682637*nuVtSqSum[3]*fl[31]-0.270632938682637*nuVtSqSum[2]*fr[30]+0.270632938682637*nuVtSqSum[2]*fl[30]-0.270632938682637*nuVtSqSum[1]*fr[29]+0.270632938682637*nuVtSqSum[1]*fl[29]+0.28125*nuVtSqSum[3]*fr[28]+0.28125*nuVtSqSum[3]*fl[28]-0.5625*nuVtSqSum[3]*fc[28]-0.270632938682637*nuVtSqSum[0]*fr[25]+0.270632938682637*nuVtSqSum[0]*fl[25]+0.28125*nuVtSqSum[2]*fr[24]+0.28125*nuVtSqSum[2]*fl[24]-0.5625*nuVtSqSum[2]*fc[24]+0.28125*nuVtSqSum[1]*fr[23]+0.28125*nuVtSqSum[1]*fl[23]-0.5625*nuVtSqSum[1]*fc[23]+0.28125*nuVtSqSum[0]*fr[15]+0.28125*nuVtSqSum[0]*fl[15]-0.5625*nuVtSqSum[0]*fc[15]; 
  diff_incr[16] = (-0.21875*nuVtSqSum[0]*fr[16])-0.21875*nuVtSqSum[0]*fl[16]-1.4375*nuVtSqSum[0]*fc[16]-0.21875*nuVtSqSum[1]*fr[8]-0.21875*nuVtSqSum[1]*fl[8]-1.4375*nuVtSqSum[1]*fc[8]-0.21875*nuVtSqSum[2]*fr[7]-0.21875*nuVtSqSum[2]*fl[7]-1.4375*nuVtSqSum[2]*fc[7]+0.270632938682637*nuVtSqSum[0]*fr[6]-0.270632938682637*nuVtSqSum[0]*fl[6]-0.21875*fr[3]*nuVtSqSum[3]-0.21875*fl[3]*nuVtSqSum[3]-1.4375*fc[3]*nuVtSqSum[3]+0.270632938682637*fr[0]*nuVtSqSum[3]-0.270632938682637*fl[0]*nuVtSqSum[3]+0.270632938682637*fr[1]*nuVtSqSum[2]-0.270632938682637*fl[1]*nuVtSqSum[2]+0.270632938682637*nuVtSqSum[1]*fr[2]-0.270632938682637*nuVtSqSum[1]*fl[2]; 
  diff_incr[17] = (-0.270632938682637*nuVtSqSum[0]*fr[26])+0.270632938682637*nuVtSqSum[0]*fl[26]-0.270632938682637*nuVtSqSum[1]*fr[19]+0.270632938682637*nuVtSqSum[1]*fl[19]-0.270632938682637*nuVtSqSum[2]*fr[18]+0.270632938682637*nuVtSqSum[2]*fl[18]+0.28125*nuVtSqSum[0]*fr[17]+0.28125*nuVtSqSum[0]*fl[17]-0.5625*nuVtSqSum[0]*fc[17]-0.270632938682637*nuVtSqSum[3]*fr[11]+0.270632938682637*nuVtSqSum[3]*fl[11]+0.28125*nuVtSqSum[1]*fr[10]+0.28125*nuVtSqSum[1]*fl[10]-0.5625*nuVtSqSum[1]*fc[10]+0.28125*nuVtSqSum[2]*fr[9]+0.28125*nuVtSqSum[2]*fl[9]-0.5625*nuVtSqSum[2]*fc[9]+0.28125*nuVtSqSum[3]*fr[4]+0.28125*nuVtSqSum[3]*fl[4]-0.5625*nuVtSqSum[3]*fc[4]; 
  diff_incr[18] = (-0.21875*nuVtSqSum[2]*fr[26])-0.21875*nuVtSqSum[2]*fl[26]-1.4375*nuVtSqSum[2]*fc[26]-0.21875*nuVtSqSum[3]*fr[19]-0.21875*nuVtSqSum[3]*fl[19]-1.4375*nuVtSqSum[3]*fc[19]-0.21875*nuVtSqSum[0]*fr[18]-0.21875*nuVtSqSum[0]*fl[18]-1.4375*nuVtSqSum[0]*fc[18]+0.270632938682637*nuVtSqSum[2]*fr[17]-0.270632938682637*nuVtSqSum[2]*fl[17]-0.21875*nuVtSqSum[1]*fr[11]-0.21875*nuVtSqSum[1]*fl[11]-1.4375*nuVtSqSum[1]*fc[11]+0.270632938682637*nuVtSqSum[3]*fr[10]-0.270632938682637*nuVtSqSum[3]*fl[10]+0.270632938682637*nuVtSqSum[0]*fr[9]-0.270632938682637*nuVtSqSum[0]*fl[9]+0.270632938682637*nuVtSqSum[1]*fr[4]-0.270632938682637*nuVtSqSum[1]*fl[4]; 
  diff_incr[19] = (-0.21875*nuVtSqSum[1]*fr[26])-0.21875*nuVtSqSum[1]*fl[26]-1.4375*nuVtSqSum[1]*fc[26]-0.21875*nuVtSqSum[0]*fr[19]-0.21875*nuVtSqSum[0]*fl[19]-1.4375*nuVtSqSum[0]*fc[19]-0.21875*nuVtSqSum[3]*fr[18]-0.21875*nuVtSqSum[3]*fl[18]-1.4375*nuVtSqSum[3]*fc[18]+0.270632938682637*nuVtSqSum[1]*fr[17]-0.270632938682637*nuVtSqSum[1]*fl[17]-0.21875*nuVtSqSum[2]*fr[11]-0.21875*nuVtSqSum[2]*fl[11]-1.4375*nuVtSqSum[2]*fc[11]+0.270632938682637*nuVtSqSum[0]*fr[10]-0.270632938682637*nuVtSqSum[0]*fl[10]+0.270632938682637*nuVtSqSum[3]*fr[9]-0.270632938682637*nuVtSqSum[3]*fl[9]+0.270632938682637*nuVtSqSum[2]*fr[4]-0.270632938682637*nuVtSqSum[2]*fl[4]; 
  diff_incr[20] = (-0.270632938682637*nuVtSqSum[0]*fr[27])+0.270632938682637*nuVtSqSum[0]*fl[27]-0.270632938682637*nuVtSqSum[1]*fr[22]+0.270632938682637*nuVtSqSum[1]*fl[22]-0.270632938682637*nuVtSqSum[2]*fr[21]+0.270632938682637*nuVtSqSum[2]*fl[21]+0.28125*nuVtSqSum[0]*fr[20]+0.28125*nuVtSqSum[0]*fl[20]-0.5625*nuVtSqSum[0]*fc[20]-0.270632938682637*nuVtSqSum[3]*fr[14]+0.270632938682637*nuVtSqSum[3]*fl[14]+0.28125*nuVtSqSum[1]*fr[13]+0.28125*nuVtSqSum[1]*fl[13]-0.5625*nuVtSqSum[1]*fc[13]+0.28125*nuVtSqSum[2]*fr[12]+0.28125*nuVtSqSum[2]*fl[12]-0.5625*nuVtSqSum[2]*fc[12]+0.28125*nuVtSqSum[3]*fr[5]+0.28125*nuVtSqSum[3]*fl[5]-0.5625*nuVtSqSum[3]*fc[5]; 
  diff_incr[21] = (-0.21875*nuVtSqSum[2]*fr[27])-0.21875*nuVtSqSum[2]*fl[27]-1.4375*nuVtSqSum[2]*fc[27]-0.21875*nuVtSqSum[3]*fr[22]-0.21875*nuVtSqSum[3]*fl[22]-1.4375*nuVtSqSum[3]*fc[22]-0.21875*nuVtSqSum[0]*fr[21]-0.21875*nuVtSqSum[0]*fl[21]-1.4375*nuVtSqSum[0]*fc[21]+0.270632938682637*nuVtSqSum[2]*fr[20]-0.270632938682637*nuVtSqSum[2]*fl[20]-0.21875*nuVtSqSum[1]*fr[14]-0.21875*nuVtSqSum[1]*fl[14]-1.4375*nuVtSqSum[1]*fc[14]+0.270632938682637*nuVtSqSum[3]*fr[13]-0.270632938682637*nuVtSqSum[3]*fl[13]+0.270632938682637*nuVtSqSum[0]*fr[12]-0.270632938682637*nuVtSqSum[0]*fl[12]+0.270632938682637*nuVtSqSum[1]*fr[5]-0.270632938682637*nuVtSqSum[1]*fl[5]; 
  diff_incr[22] = (-0.21875*nuVtSqSum[1]*fr[27])-0.21875*nuVtSqSum[1]*fl[27]-1.4375*nuVtSqSum[1]*fc[27]-0.21875*nuVtSqSum[0]*fr[22]-0.21875*nuVtSqSum[0]*fl[22]-1.4375*nuVtSqSum[0]*fc[22]-0.21875*nuVtSqSum[3]*fr[21]-0.21875*nuVtSqSum[3]*fl[21]-1.4375*nuVtSqSum[3]*fc[21]+0.270632938682637*nuVtSqSum[1]*fr[20]-0.270632938682637*nuVtSqSum[1]*fl[20]-0.21875*nuVtSqSum[2]*fr[14]-0.21875*nuVtSqSum[2]*fl[14]-1.4375*nuVtSqSum[2]*fc[14]+0.270632938682637*nuVtSqSum[0]*fr[13]-0.270632938682637*nuVtSqSum[0]*fl[13]+0.270632938682637*nuVtSqSum[3]*fr[12]-0.270632938682637*nuVtSqSum[3]*fl[12]+0.270632938682637*nuVtSqSum[2]*fr[5]-0.270632938682637*nuVtSqSum[2]*fl[5]; 
  diff_incr[23] = (-0.270632938682637*nuVtSqSum[2]*fr[31])+0.270632938682637*nuVtSqSum[2]*fl[31]-0.270632938682637*nuVtSqSum[3]*fr[30]+0.270632938682637*nuVtSqSum[3]*fl[30]-0.270632938682637*nuVtSqSum[0]*fr[29]+0.270632938682637*nuVtSqSum[0]*fl[29]+0.28125*nuVtSqSum[2]*fr[28]+0.28125*nuVtSqSum[2]*fl[28]-0.5625*nuVtSqSum[2]*fc[28]-0.270632938682637*nuVtSqSum[1]*fr[25]+0.270632938682637*nuVtSqSum[1]*fl[25]+0.28125*nuVtSqSum[3]*fr[24]+0.28125*nuVtSqSum[3]*fl[24]-0.5625*nuVtSqSum[3]*fc[24]+0.28125*nuVtSqSum[0]*fr[23]+0.28125*nuVtSqSum[0]*fl[23]-0.5625*nuVtSqSum[0]*fc[23]+0.28125*nuVtSqSum[1]*fr[15]+0.28125*nuVtSqSum[1]*fl[15]-0.5625*nuVtSqSum[1]*fc[15]; 
  diff_incr[24] = (-0.270632938682637*nuVtSqSum[1]*fr[31])+0.270632938682637*nuVtSqSum[1]*fl[31]-0.270632938682637*nuVtSqSum[0]*fr[30]+0.270632938682637*nuVtSqSum[0]*fl[30]-0.270632938682637*nuVtSqSum[3]*fr[29]+0.270632938682637*nuVtSqSum[3]*fl[29]+0.28125*nuVtSqSum[1]*fr[28]+0.28125*nuVtSqSum[1]*fl[28]-0.5625*nuVtSqSum[1]*fc[28]-0.270632938682637*nuVtSqSum[2]*fr[25]+0.270632938682637*nuVtSqSum[2]*fl[25]+0.28125*nuVtSqSum[0]*fr[24]+0.28125*nuVtSqSum[0]*fl[24]-0.5625*nuVtSqSum[0]*fc[24]+0.28125*nuVtSqSum[3]*fr[23]+0.28125*nuVtSqSum[3]*fl[23]-0.5625*nuVtSqSum[3]*fc[23]+0.28125*nuVtSqSum[2]*fr[15]+0.28125*nuVtSqSum[2]*fl[15]-0.5625*nuVtSqSum[2]*fc[15]; 
  diff_incr[25] = (-0.21875*nuVtSqSum[3]*fr[31])-0.21875*nuVtSqSum[3]*fl[31]-1.4375*nuVtSqSum[3]*fc[31]-0.21875*nuVtSqSum[2]*fr[30]-0.21875*nuVtSqSum[2]*fl[30]-1.4375*nuVtSqSum[2]*fc[30]-0.21875*nuVtSqSum[1]*fr[29]-0.21875*nuVtSqSum[1]*fl[29]-1.4375*nuVtSqSum[1]*fc[29]+0.270632938682637*nuVtSqSum[3]*fr[28]-0.270632938682637*nuVtSqSum[3]*fl[28]-0.21875*nuVtSqSum[0]*fr[25]-0.21875*nuVtSqSum[0]*fl[25]-1.4375*nuVtSqSum[0]*fc[25]+0.270632938682637*nuVtSqSum[2]*fr[24]-0.270632938682637*nuVtSqSum[2]*fl[24]+0.270632938682637*nuVtSqSum[1]*fr[23]-0.270632938682637*nuVtSqSum[1]*fl[23]+0.270632938682637*nuVtSqSum[0]*fr[15]-0.270632938682637*nuVtSqSum[0]*fl[15]; 
  diff_incr[26] = (-0.21875*nuVtSqSum[0]*fr[26])-0.21875*nuVtSqSum[0]*fl[26]-1.4375*nuVtSqSum[0]*fc[26]-0.21875*nuVtSqSum[1]*fr[19]-0.21875*nuVtSqSum[1]*fl[19]-1.4375*nuVtSqSum[1]*fc[19]-0.21875*nuVtSqSum[2]*fr[18]-0.21875*nuVtSqSum[2]*fl[18]-1.4375*nuVtSqSum[2]*fc[18]+0.270632938682637*nuVtSqSum[0]*fr[17]-0.270632938682637*nuVtSqSum[0]*fl[17]-0.21875*nuVtSqSum[3]*fr[11]-0.21875*nuVtSqSum[3]*fl[11]-1.4375*nuVtSqSum[3]*fc[11]+0.270632938682637*nuVtSqSum[1]*fr[10]-0.270632938682637*nuVtSqSum[1]*fl[10]+0.270632938682637*nuVtSqSum[2]*fr[9]-0.270632938682637*nuVtSqSum[2]*fl[9]+0.270632938682637*nuVtSqSum[3]*fr[4]-0.270632938682637*nuVtSqSum[3]*fl[4]; 
  diff_incr[27] = (-0.21875*nuVtSqSum[0]*fr[27])-0.21875*nuVtSqSum[0]*fl[27]-1.4375*nuVtSqSum[0]*fc[27]-0.21875*nuVtSqSum[1]*fr[22]-0.21875*nuVtSqSum[1]*fl[22]-1.4375*nuVtSqSum[1]*fc[22]-0.21875*nuVtSqSum[2]*fr[21]-0.21875*nuVtSqSum[2]*fl[21]-1.4375*nuVtSqSum[2]*fc[21]+0.270632938682637*nuVtSqSum[0]*fr[20]-0.270632938682637*nuVtSqSum[0]*fl[20]-0.21875*nuVtSqSum[3]*fr[14]-0.21875*nuVtSqSum[3]*fl[14]-1.4375*nuVtSqSum[3]*fc[14]+0.270632938682637*nuVtSqSum[1]*fr[13]-0.270632938682637*nuVtSqSum[1]*fl[13]+0.270632938682637*nuVtSqSum[2]*fr[12]-0.270632938682637*nuVtSqSum[2]*fl[12]+0.270632938682637*nuVtSqSum[3]*fr[5]-0.270632938682637*nuVtSqSum[3]*fl[5]; 
  diff_incr[28] = (-0.270632938682637*nuVtSqSum[0]*fr[31])+0.270632938682637*nuVtSqSum[0]*fl[31]-0.270632938682637*nuVtSqSum[1]*fr[30]+0.270632938682637*nuVtSqSum[1]*fl[30]-0.270632938682637*nuVtSqSum[2]*fr[29]+0.270632938682637*nuVtSqSum[2]*fl[29]+0.28125*nuVtSqSum[0]*fr[28]+0.28125*nuVtSqSum[0]*fl[28]-0.5625*nuVtSqSum[0]*fc[28]-0.270632938682637*nuVtSqSum[3]*fr[25]+0.270632938682637*nuVtSqSum[3]*fl[25]+0.28125*nuVtSqSum[1]*fr[24]+0.28125*nuVtSqSum[1]*fl[24]-0.5625*nuVtSqSum[1]*fc[24]+0.28125*nuVtSqSum[2]*fr[23]+0.28125*nuVtSqSum[2]*fl[23]-0.5625*nuVtSqSum[2]*fc[23]+0.28125*nuVtSqSum[3]*fr[15]+0.28125*nuVtSqSum[3]*fl[15]-0.5625*nuVtSqSum[3]*fc[15]; 
  diff_incr[29] = (-0.21875*nuVtSqSum[2]*fr[31])-0.21875*nuVtSqSum[2]*fl[31]-1.4375*nuVtSqSum[2]*fc[31]-0.21875*nuVtSqSum[3]*fr[30]-0.21875*nuVtSqSum[3]*fl[30]-1.4375*nuVtSqSum[3]*fc[30]-0.21875*nuVtSqSum[0]*fr[29]-0.21875*nuVtSqSum[0]*fl[29]-1.4375*nuVtSqSum[0]*fc[29]+0.270632938682637*nuVtSqSum[2]*fr[28]-0.270632938682637*nuVtSqSum[2]*fl[28]-0.21875*nuVtSqSum[1]*fr[25]-0.21875*nuVtSqSum[1]*fl[25]-1.4375*nuVtSqSum[1]*fc[25]+0.270632938682637*nuVtSqSum[3]*fr[24]-0.270632938682637*nuVtSqSum[3]*fl[24]+0.270632938682637*nuVtSqSum[0]*fr[23]-0.270632938682637*nuVtSqSum[0]*fl[23]+0.270632938682637*nuVtSqSum[1]*fr[15]-0.270632938682637*nuVtSqSum[1]*fl[15]; 
  diff_incr[30] = (-0.21875*nuVtSqSum[1]*fr[31])-0.21875*nuVtSqSum[1]*fl[31]-1.4375*nuVtSqSum[1]*fc[31]-0.21875*nuVtSqSum[0]*fr[30]-0.21875*nuVtSqSum[0]*fl[30]-1.4375*nuVtSqSum[0]*fc[30]-0.21875*nuVtSqSum[3]*fr[29]-0.21875*nuVtSqSum[3]*fl[29]-1.4375*nuVtSqSum[3]*fc[29]+0.270632938682637*nuVtSqSum[1]*fr[28]-0.270632938682637*nuVtSqSum[1]*fl[28]-0.21875*nuVtSqSum[2]*fr[25]-0.21875*nuVtSqSum[2]*fl[25]-1.4375*nuVtSqSum[2]*fc[25]+0.270632938682637*nuVtSqSum[0]*fr[24]-0.270632938682637*nuVtSqSum[0]*fl[24]+0.270632938682637*nuVtSqSum[3]*fr[23]-0.270632938682637*nuVtSqSum[3]*fl[23]+0.270632938682637*nuVtSqSum[2]*fr[15]-0.270632938682637*nuVtSqSum[2]*fl[15]; 
  diff_incr[31] = (-0.21875*nuVtSqSum[0]*fr[31])-0.21875*nuVtSqSum[0]*fl[31]-1.4375*nuVtSqSum[0]*fc[31]-0.21875*nuVtSqSum[1]*fr[30]-0.21875*nuVtSqSum[1]*fl[30]-1.4375*nuVtSqSum[1]*fc[30]-0.21875*nuVtSqSum[2]*fr[29]-0.21875*nuVtSqSum[2]*fl[29]-1.4375*nuVtSqSum[2]*fc[29]+0.270632938682637*nuVtSqSum[0]*fr[28]-0.270632938682637*nuVtSqSum[0]*fl[28]-0.21875*nuVtSqSum[3]*fr[25]-0.21875*nuVtSqSum[3]*fl[25]-1.4375*nuVtSqSum[3]*fc[25]+0.270632938682637*nuVtSqSum[1]*fr[24]-0.270632938682637*nuVtSqSum[1]*fl[24]+0.270632938682637*nuVtSqSum[2]*fr[23]-0.270632938682637*nuVtSqSum[2]*fl[23]+0.270632938682637*nuVtSqSum[3]*fr[15]-0.270632938682637*nuVtSqSum[3]*fl[15]; 

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
} 
