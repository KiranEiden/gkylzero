#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH void lbo_gyrokinetic_diff_surfmu_1x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // m_:            species mass.
  // bmag_inv:      1/(magnetic field magnitude). 
  // w[3]:         cell-center coordinates. 
  // dxv[3]:       cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[6]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      distribution function in cells 
  // out:           incremented distribution function in cell 
  double rdvSq4 = 4.0/(dxv[2]*dxv[2]); 
  double temp_diff[20] = {0.0}; 
  double diff_incr[20] = {0.0}; 

  temp_diff[0] = 0.6708203932499369*w[2]*fr[9]+0.3354101966249685*dxv[2]*fr[9]+0.6708203932499369*w[2]*fl[9]-0.3354101966249685*dxv[2]*fl[9]-1.341640786499874*w[2]*fc[9]-1.190784930203603*w[2]*fr[3]-0.5953924651018015*dxv[2]*fr[3]+1.190784930203603*w[2]*fl[3]-0.5953924651018015*dxv[2]*fl[3]-1.190784930203603*dxv[2]*fc[3]+0.9375*fr[0]*w[2]+0.9375*fl[0]*w[2]-1.875*fc[0]*w[2]+0.46875*fr[0]*dxv[2]-0.46875*fl[0]*dxv[2]; 
  temp_diff[1] = 0.6708203932499369*w[2]*fr[15]+0.3354101966249685*dxv[2]*fr[15]+0.6708203932499369*w[2]*fl[15]-0.3354101966249685*dxv[2]*fl[15]-1.341640786499874*w[2]*fc[15]-1.190784930203603*w[2]*fr[5]-0.5953924651018015*dxv[2]*fr[5]+1.190784930203603*w[2]*fl[5]-0.5953924651018015*dxv[2]*fl[5]-1.190784930203603*dxv[2]*fc[5]+0.9375*fr[1]*w[2]+0.9375*fl[1]*w[2]-1.875*fc[1]*w[2]+0.46875*fr[1]*dxv[2]-0.46875*fl[1]*dxv[2]; 
  temp_diff[2] = 0.6708203932499369*w[2]*fr[16]+0.3354101966249685*dxv[2]*fr[16]+0.6708203932499369*w[2]*fl[16]-0.3354101966249685*dxv[2]*fl[16]-1.341640786499874*w[2]*fc[16]-1.190784930203603*w[2]*fr[6]-0.5953924651018015*dxv[2]*fr[6]+1.190784930203603*w[2]*fl[6]-0.5953924651018015*dxv[2]*fl[6]-1.190784930203603*dxv[2]*fc[6]+0.9375*fr[2]*w[2]+0.9375*fl[2]*w[2]-1.875*fc[2]*w[2]+0.46875*dxv[2]*fr[2]-0.46875*dxv[2]*fl[2]; 
  temp_diff[3] = 0.7382874503707888*w[2]*fr[9]+0.3691437251853944*dxv[2]*fr[9]-0.7382874503707888*w[2]*fl[9]+0.3691437251853944*dxv[2]*fl[9]-1.585502557353661*dxv[2]*fc[9]-1.453125*w[2]*fr[3]-0.7265625*dxv[2]*fr[3]-1.453125*w[2]*fl[3]+0.7265625*dxv[2]*fl[3]-5.34375*w[2]*fc[3]+1.190784930203603*fr[0]*w[2]-1.190784930203603*fl[0]*w[2]+0.5953924651018015*fr[0]*dxv[2]+0.5953924651018015*fl[0]*dxv[2]-1.190784930203603*fc[0]*dxv[2]; 
  temp_diff[4] = 0.6708203932499369*w[2]*fr[19]+0.3354101966249685*dxv[2]*fr[19]+0.6708203932499369*w[2]*fl[19]-0.3354101966249685*dxv[2]*fl[19]-1.341640786499874*w[2]*fc[19]-1.190784930203603*w[2]*fr[10]-0.5953924651018015*dxv[2]*fr[10]+1.190784930203603*w[2]*fl[10]-0.5953924651018015*dxv[2]*fl[10]-1.190784930203603*dxv[2]*fc[10]+0.9375*w[2]*fr[4]+0.46875*dxv[2]*fr[4]+0.9375*w[2]*fl[4]-0.46875*dxv[2]*fl[4]-1.875*w[2]*fc[4]; 
  temp_diff[5] = 0.7382874503707888*w[2]*fr[15]+0.3691437251853944*dxv[2]*fr[15]-0.7382874503707888*w[2]*fl[15]+0.3691437251853944*dxv[2]*fl[15]-1.585502557353661*dxv[2]*fc[15]-1.453125*w[2]*fr[5]-0.7265625*dxv[2]*fr[5]-1.453125*w[2]*fl[5]+0.7265625*dxv[2]*fl[5]-5.34375*w[2]*fc[5]+1.190784930203603*fr[1]*w[2]-1.190784930203603*fl[1]*w[2]+0.5953924651018015*fr[1]*dxv[2]+0.5953924651018015*fl[1]*dxv[2]-1.190784930203603*fc[1]*dxv[2]; 
  temp_diff[6] = 0.7382874503707888*w[2]*fr[16]+0.3691437251853944*dxv[2]*fr[16]-0.7382874503707888*w[2]*fl[16]+0.3691437251853944*dxv[2]*fl[16]-1.585502557353661*dxv[2]*fc[16]-1.453125*w[2]*fr[6]-0.7265625*dxv[2]*fr[6]-1.453125*w[2]*fl[6]+0.7265625*dxv[2]*fl[6]-5.34375*w[2]*fc[6]+1.190784930203603*fr[2]*w[2]-1.190784930203603*fl[2]*w[2]+0.5953924651018015*dxv[2]*fr[2]+0.5953924651018015*dxv[2]*fl[2]-1.190784930203603*dxv[2]*fc[2]; 
  temp_diff[7] = (-1.190784930203603*w[2]*fr[13])-0.5953924651018015*dxv[2]*fr[13]+1.190784930203603*w[2]*fl[13]-0.5953924651018015*dxv[2]*fl[13]-1.190784930203603*dxv[2]*fc[13]+0.9375*w[2]*fr[7]+0.46875*dxv[2]*fr[7]+0.9375*w[2]*fl[7]-0.46875*dxv[2]*fl[7]-1.875*w[2]*fc[7]; 
  temp_diff[8] = (-1.190784930203603*w[2]*fr[14])-0.5953924651018015*dxv[2]*fr[14]+1.190784930203603*w[2]*fl[14]-0.5953924651018015*dxv[2]*fl[14]-1.190784930203603*dxv[2]*fc[14]+0.9375*w[2]*fr[8]+0.46875*dxv[2]*fr[8]+0.9375*w[2]*fl[8]-0.46875*dxv[2]*fl[8]-1.875*w[2]*fc[8]; 
  temp_diff[9] = (-0.140625*w[2]*fr[9])-0.0703125*dxv[2]*fr[9]-0.140625*w[2]*fl[9]+0.0703125*dxv[2]*fl[9]-6.28125*w[2]*fc[9]-0.3025768239224545*w[2]*fr[3]-0.1512884119612272*dxv[2]*fr[3]+0.3025768239224545*w[2]*fl[3]-0.1512884119612272*dxv[2]*fl[3]-1.149791930905327*dxv[2]*fc[3]+0.4192627457812106*fr[0]*w[2]+0.4192627457812106*fl[0]*w[2]-0.8385254915624212*fc[0]*w[2]+0.2096313728906053*fr[0]*dxv[2]-0.2096313728906053*fl[0]*dxv[2]; 
  temp_diff[10] = 0.7382874503707888*w[2]*fr[19]+0.3691437251853944*dxv[2]*fr[19]-0.7382874503707888*w[2]*fl[19]+0.3691437251853944*dxv[2]*fl[19]-1.585502557353661*dxv[2]*fc[19]-1.453125*w[2]*fr[10]-0.7265625*dxv[2]*fr[10]-1.453125*w[2]*fl[10]+0.7265625*dxv[2]*fl[10]-5.34375*w[2]*fc[10]+1.190784930203603*w[2]*fr[4]+0.5953924651018015*dxv[2]*fr[4]-1.190784930203603*w[2]*fl[4]+0.5953924651018015*dxv[2]*fl[4]-1.190784930203603*dxv[2]*fc[4]; 
  temp_diff[11] = (-1.190784930203603*w[2]*fr[17])-0.5953924651018015*dxv[2]*fr[17]+1.190784930203603*w[2]*fl[17]-0.5953924651018015*dxv[2]*fl[17]-1.190784930203603*dxv[2]*fc[17]+0.9375*w[2]*fr[11]+0.46875*dxv[2]*fr[11]+0.9375*w[2]*fl[11]-0.46875*dxv[2]*fl[11]-1.875*w[2]*fc[11]; 
  temp_diff[12] = (-1.190784930203603*w[2]*fr[18])-0.5953924651018015*dxv[2]*fr[18]+1.190784930203603*w[2]*fl[18]-0.5953924651018015*dxv[2]*fl[18]-1.190784930203603*dxv[2]*fc[18]+0.9375*w[2]*fr[12]+0.46875*dxv[2]*fr[12]+0.9375*w[2]*fl[12]-0.46875*dxv[2]*fl[12]-1.875*w[2]*fc[12]; 
  temp_diff[13] = (-1.453125*w[2]*fr[13])-0.7265625*dxv[2]*fr[13]-1.453125*w[2]*fl[13]+0.7265625*dxv[2]*fl[13]-5.34375*w[2]*fc[13]+1.190784930203603*w[2]*fr[7]+0.5953924651018015*dxv[2]*fr[7]-1.190784930203603*w[2]*fl[7]+0.5953924651018015*dxv[2]*fl[7]-1.190784930203603*dxv[2]*fc[7]; 
  temp_diff[14] = (-1.453125*w[2]*fr[14])-0.7265625*dxv[2]*fr[14]-1.453125*w[2]*fl[14]+0.7265625*dxv[2]*fl[14]-5.34375*w[2]*fc[14]+1.190784930203603*w[2]*fr[8]+0.5953924651018015*dxv[2]*fr[8]-1.190784930203603*w[2]*fl[8]+0.5953924651018015*dxv[2]*fl[8]-1.190784930203603*dxv[2]*fc[8]; 
  temp_diff[15] = (-0.140625*w[2]*fr[15])-0.0703125*dxv[2]*fr[15]-0.140625*w[2]*fl[15]+0.0703125*dxv[2]*fl[15]-6.28125*w[2]*fc[15]-0.3025768239224544*w[2]*fr[5]-0.1512884119612272*dxv[2]*fr[5]+0.3025768239224544*w[2]*fl[5]-0.1512884119612272*dxv[2]*fl[5]-1.149791930905327*dxv[2]*fc[5]+0.4192627457812105*fr[1]*w[2]+0.4192627457812105*fl[1]*w[2]-0.8385254915624211*fc[1]*w[2]+0.2096313728906053*fr[1]*dxv[2]-0.2096313728906053*fl[1]*dxv[2]; 
  temp_diff[16] = (-0.140625*w[2]*fr[16])-0.0703125*dxv[2]*fr[16]-0.140625*w[2]*fl[16]+0.0703125*dxv[2]*fl[16]-6.28125*w[2]*fc[16]-0.3025768239224544*w[2]*fr[6]-0.1512884119612272*dxv[2]*fr[6]+0.3025768239224544*w[2]*fl[6]-0.1512884119612272*dxv[2]*fl[6]-1.149791930905327*dxv[2]*fc[6]+0.4192627457812105*fr[2]*w[2]+0.4192627457812105*fl[2]*w[2]-0.8385254915624211*fc[2]*w[2]+0.2096313728906053*dxv[2]*fr[2]-0.2096313728906053*dxv[2]*fl[2]; 
  temp_diff[17] = (-1.453125*w[2]*fr[17])-0.7265625*dxv[2]*fr[17]-1.453125*w[2]*fl[17]+0.7265625*dxv[2]*fl[17]-5.34375*w[2]*fc[17]+1.190784930203603*w[2]*fr[11]+0.5953924651018015*dxv[2]*fr[11]-1.190784930203603*w[2]*fl[11]+0.5953924651018015*dxv[2]*fl[11]-1.190784930203603*dxv[2]*fc[11]; 
  temp_diff[18] = (-1.453125*w[2]*fr[18])-0.7265625*dxv[2]*fr[18]-1.453125*w[2]*fl[18]+0.7265625*dxv[2]*fl[18]-5.34375*w[2]*fc[18]+1.190784930203603*w[2]*fr[12]+0.5953924651018015*dxv[2]*fr[12]-1.190784930203603*w[2]*fl[12]+0.5953924651018015*dxv[2]*fl[12]-1.190784930203603*dxv[2]*fc[12]; 
  temp_diff[19] = (-0.140625*w[2]*fr[19])-0.0703125*dxv[2]*fr[19]-0.140625*w[2]*fl[19]+0.0703125*dxv[2]*fl[19]-6.28125*w[2]*fc[19]-0.3025768239224545*w[2]*fr[10]-0.1512884119612272*dxv[2]*fr[10]+0.3025768239224545*w[2]*fl[10]-0.1512884119612272*dxv[2]*fl[10]-1.149791930905327*dxv[2]*fc[10]+0.4192627457812106*w[2]*fr[4]+0.2096313728906053*dxv[2]*fr[4]+0.4192627457812106*w[2]*fl[4]-0.2096313728906053*dxv[2]*fl[4]-0.8385254915624212*w[2]*fc[4]; 

  diff_incr[0] = 0.6388765649999399*bmag_inv[2]*nuVtSqSum[2]*temp_diff[7]*m_+bmag_inv[0]*nuVtSqSum[2]*temp_diff[7]*m_+nuVtSqSum[0]*bmag_inv[2]*temp_diff[7]*m_+0.8944271909999159*bmag_inv[1]*nuVtSqSum[1]*temp_diff[7]*m_+temp_diff[0]*bmag_inv[2]*nuVtSqSum[2]*m_+0.8944271909999159*bmag_inv[1]*temp_diff[1]*nuVtSqSum[2]*m_+0.8944271909999159*nuVtSqSum[1]*temp_diff[1]*bmag_inv[2]*m_+bmag_inv[0]*nuVtSqSum[1]*temp_diff[1]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_diff[1]*m_+temp_diff[0]*bmag_inv[1]*nuVtSqSum[1]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[0]*m_; 
  diff_incr[1] = 1.571428571428571*bmag_inv[1]*nuVtSqSum[2]*temp_diff[7]*m_+1.571428571428571*nuVtSqSum[1]*bmag_inv[2]*temp_diff[7]*m_+0.8944271909999159*bmag_inv[0]*nuVtSqSum[1]*temp_diff[7]*m_+0.8944271909999159*nuVtSqSum[0]*bmag_inv[1]*temp_diff[7]*m_+1.571428571428571*temp_diff[1]*bmag_inv[2]*nuVtSqSum[2]*m_+0.8944271909999159*bmag_inv[0]*temp_diff[1]*nuVtSqSum[2]*m_+0.8944271909999159*temp_diff[0]*bmag_inv[1]*nuVtSqSum[2]*m_+0.8944271909999159*nuVtSqSum[0]*temp_diff[1]*bmag_inv[2]*m_+0.8944271909999159*temp_diff[0]*nuVtSqSum[1]*bmag_inv[2]*m_+1.8*bmag_inv[1]*nuVtSqSum[1]*temp_diff[1]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[1]*m_+bmag_inv[0]*temp_diff[0]*nuVtSqSum[1]*m_+nuVtSqSum[0]*temp_diff[0]*bmag_inv[1]*m_; 
  diff_incr[2] = 0.63887656499994*bmag_inv[2]*nuVtSqSum[2]*temp_diff[11]*m_+1.0*bmag_inv[0]*nuVtSqSum[2]*temp_diff[11]*m_+1.0*nuVtSqSum[0]*bmag_inv[2]*temp_diff[11]*m_+0.8944271909999161*bmag_inv[1]*nuVtSqSum[1]*temp_diff[11]*m_+0.8944271909999159*bmag_inv[1]*nuVtSqSum[2]*temp_diff[4]*m_+0.8944271909999159*nuVtSqSum[1]*bmag_inv[2]*temp_diff[4]*m_+bmag_inv[0]*nuVtSqSum[1]*temp_diff[4]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_diff[4]*m_+bmag_inv[2]*nuVtSqSum[2]*temp_diff[2]*m_+bmag_inv[1]*nuVtSqSum[1]*temp_diff[2]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[2]*m_; 
  diff_incr[3] = 0.63887656499994*bmag_inv[2]*nuVtSqSum[2]*temp_diff[13]*m_+1.0*bmag_inv[0]*nuVtSqSum[2]*temp_diff[13]*m_+1.0*nuVtSqSum[0]*bmag_inv[2]*temp_diff[13]*m_+0.8944271909999161*bmag_inv[1]*nuVtSqSum[1]*temp_diff[13]*m_+0.8944271909999159*bmag_inv[1]*nuVtSqSum[2]*temp_diff[5]*m_+0.8944271909999159*nuVtSqSum[1]*bmag_inv[2]*temp_diff[5]*m_+bmag_inv[0]*nuVtSqSum[1]*temp_diff[5]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_diff[5]*m_+bmag_inv[2]*nuVtSqSum[2]*temp_diff[3]*m_+bmag_inv[1]*nuVtSqSum[1]*temp_diff[3]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[3]*m_; 
  diff_incr[4] = 1.571428571428572*bmag_inv[1]*nuVtSqSum[2]*temp_diff[11]*m_+1.571428571428572*nuVtSqSum[1]*bmag_inv[2]*temp_diff[11]*m_+0.8944271909999161*bmag_inv[0]*nuVtSqSum[1]*temp_diff[11]*m_+0.8944271909999161*nuVtSqSum[0]*bmag_inv[1]*temp_diff[11]*m_+1.571428571428571*bmag_inv[2]*nuVtSqSum[2]*temp_diff[4]*m_+0.8944271909999159*bmag_inv[0]*nuVtSqSum[2]*temp_diff[4]*m_+0.8944271909999159*nuVtSqSum[0]*bmag_inv[2]*temp_diff[4]*m_+1.8*bmag_inv[1]*nuVtSqSum[1]*temp_diff[4]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[4]*m_+0.8944271909999159*bmag_inv[1]*nuVtSqSum[2]*temp_diff[2]*m_+0.8944271909999159*nuVtSqSum[1]*bmag_inv[2]*temp_diff[2]*m_+bmag_inv[0]*nuVtSqSum[1]*temp_diff[2]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_diff[2]*m_; 
  diff_incr[5] = 1.571428571428572*bmag_inv[1]*nuVtSqSum[2]*temp_diff[13]*m_+1.571428571428572*nuVtSqSum[1]*bmag_inv[2]*temp_diff[13]*m_+0.8944271909999161*bmag_inv[0]*nuVtSqSum[1]*temp_diff[13]*m_+0.8944271909999161*nuVtSqSum[0]*bmag_inv[1]*temp_diff[13]*m_+1.571428571428571*bmag_inv[2]*nuVtSqSum[2]*temp_diff[5]*m_+0.8944271909999159*bmag_inv[0]*nuVtSqSum[2]*temp_diff[5]*m_+0.8944271909999159*nuVtSqSum[0]*bmag_inv[2]*temp_diff[5]*m_+1.8*bmag_inv[1]*nuVtSqSum[1]*temp_diff[5]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[5]*m_+0.8944271909999159*bmag_inv[1]*nuVtSqSum[2]*temp_diff[3]*m_+0.8944271909999159*nuVtSqSum[1]*bmag_inv[2]*temp_diff[3]*m_+bmag_inv[0]*nuVtSqSum[1]*temp_diff[3]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_diff[3]*m_; 
  diff_incr[6] = 0.6388765649999399*bmag_inv[2]*nuVtSqSum[2]*temp_diff[17]*m_+bmag_inv[0]*nuVtSqSum[2]*temp_diff[17]*m_+nuVtSqSum[0]*bmag_inv[2]*temp_diff[17]*m_+0.8944271909999159*bmag_inv[1]*nuVtSqSum[1]*temp_diff[17]*m_+0.8944271909999159*bmag_inv[1]*nuVtSqSum[2]*temp_diff[10]*m_+0.8944271909999159*nuVtSqSum[1]*bmag_inv[2]*temp_diff[10]*m_+bmag_inv[0]*nuVtSqSum[1]*temp_diff[10]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_diff[10]*m_+bmag_inv[2]*nuVtSqSum[2]*temp_diff[6]*m_+bmag_inv[1]*nuVtSqSum[1]*temp_diff[6]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[6]*m_; 
  diff_incr[7] = 2.142857142857143*bmag_inv[2]*nuVtSqSum[2]*temp_diff[7]*m_+0.6388765649999399*bmag_inv[0]*nuVtSqSum[2]*temp_diff[7]*m_+0.6388765649999399*nuVtSqSum[0]*bmag_inv[2]*temp_diff[7]*m_+1.571428571428571*bmag_inv[1]*nuVtSqSum[1]*temp_diff[7]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[7]*m_+0.6388765649999399*temp_diff[0]*bmag_inv[2]*nuVtSqSum[2]*m_+1.571428571428571*bmag_inv[1]*temp_diff[1]*nuVtSqSum[2]*m_+bmag_inv[0]*temp_diff[0]*nuVtSqSum[2]*m_+1.571428571428571*nuVtSqSum[1]*temp_diff[1]*bmag_inv[2]*m_+nuVtSqSum[0]*temp_diff[0]*bmag_inv[2]*m_+0.8944271909999159*bmag_inv[0]*nuVtSqSum[1]*temp_diff[1]*m_+0.8944271909999159*nuVtSqSum[0]*bmag_inv[1]*temp_diff[1]*m_+0.8944271909999159*temp_diff[0]*bmag_inv[1]*nuVtSqSum[1]*m_; 
  diff_incr[8] = 0.8944271909999161*bmag_inv[1]*nuVtSqSum[2]*temp_diff[12]*m_+0.8944271909999161*nuVtSqSum[1]*bmag_inv[2]*temp_diff[12]*m_+1.0*bmag_inv[0]*nuVtSqSum[1]*temp_diff[12]*m_+1.0*nuVtSqSum[0]*bmag_inv[1]*temp_diff[12]*m_+bmag_inv[2]*nuVtSqSum[2]*temp_diff[8]*m_+bmag_inv[1]*nuVtSqSum[1]*temp_diff[8]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[8]*m_; 
  diff_incr[9] = 0.8944271909999161*bmag_inv[1]*nuVtSqSum[2]*temp_diff[15]*m_+0.8944271909999161*nuVtSqSum[1]*bmag_inv[2]*temp_diff[15]*m_+1.0*bmag_inv[0]*nuVtSqSum[1]*temp_diff[15]*m_+1.0*nuVtSqSum[0]*bmag_inv[1]*temp_diff[15]*m_+bmag_inv[2]*nuVtSqSum[2]*temp_diff[9]*m_+bmag_inv[1]*nuVtSqSum[1]*temp_diff[9]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[9]*m_; 
  diff_incr[10] = 1.571428571428571*bmag_inv[1]*nuVtSqSum[2]*temp_diff[17]*m_+1.571428571428571*nuVtSqSum[1]*bmag_inv[2]*temp_diff[17]*m_+0.8944271909999159*bmag_inv[0]*nuVtSqSum[1]*temp_diff[17]*m_+0.8944271909999159*nuVtSqSum[0]*bmag_inv[1]*temp_diff[17]*m_+1.571428571428571*bmag_inv[2]*nuVtSqSum[2]*temp_diff[10]*m_+0.8944271909999159*bmag_inv[0]*nuVtSqSum[2]*temp_diff[10]*m_+0.8944271909999159*nuVtSqSum[0]*bmag_inv[2]*temp_diff[10]*m_+1.8*bmag_inv[1]*nuVtSqSum[1]*temp_diff[10]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[10]*m_+0.8944271909999159*bmag_inv[1]*nuVtSqSum[2]*temp_diff[6]*m_+0.8944271909999159*nuVtSqSum[1]*bmag_inv[2]*temp_diff[6]*m_+bmag_inv[0]*nuVtSqSum[1]*temp_diff[6]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_diff[6]*m_; 
  diff_incr[11] = 2.142857142857143*bmag_inv[2]*nuVtSqSum[2]*temp_diff[11]*m_+0.6388765649999399*bmag_inv[0]*nuVtSqSum[2]*temp_diff[11]*m_+0.6388765649999399*nuVtSqSum[0]*bmag_inv[2]*temp_diff[11]*m_+1.571428571428571*bmag_inv[1]*nuVtSqSum[1]*temp_diff[11]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[11]*m_+1.571428571428572*bmag_inv[1]*nuVtSqSum[2]*temp_diff[4]*m_+1.571428571428572*nuVtSqSum[1]*bmag_inv[2]*temp_diff[4]*m_+0.8944271909999161*bmag_inv[0]*nuVtSqSum[1]*temp_diff[4]*m_+0.8944271909999161*nuVtSqSum[0]*bmag_inv[1]*temp_diff[4]*m_+0.63887656499994*bmag_inv[2]*nuVtSqSum[2]*temp_diff[2]*m_+1.0*bmag_inv[0]*nuVtSqSum[2]*temp_diff[2]*m_+1.0*nuVtSqSum[0]*bmag_inv[2]*temp_diff[2]*m_+0.8944271909999161*bmag_inv[1]*nuVtSqSum[1]*temp_diff[2]*m_; 
  diff_incr[12] = 1.571428571428571*bmag_inv[2]*nuVtSqSum[2]*temp_diff[12]*m_+0.8944271909999159*bmag_inv[0]*nuVtSqSum[2]*temp_diff[12]*m_+0.8944271909999159*nuVtSqSum[0]*bmag_inv[2]*temp_diff[12]*m_+1.8*bmag_inv[1]*nuVtSqSum[1]*temp_diff[12]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[12]*m_+0.8944271909999161*bmag_inv[1]*nuVtSqSum[2]*temp_diff[8]*m_+0.8944271909999161*nuVtSqSum[1]*bmag_inv[2]*temp_diff[8]*m_+1.0*bmag_inv[0]*nuVtSqSum[1]*temp_diff[8]*m_+1.0*nuVtSqSum[0]*bmag_inv[1]*temp_diff[8]*m_; 
  diff_incr[13] = 2.142857142857143*bmag_inv[2]*nuVtSqSum[2]*temp_diff[13]*m_+0.6388765649999399*bmag_inv[0]*nuVtSqSum[2]*temp_diff[13]*m_+0.6388765649999399*nuVtSqSum[0]*bmag_inv[2]*temp_diff[13]*m_+1.571428571428571*bmag_inv[1]*nuVtSqSum[1]*temp_diff[13]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[13]*m_+1.571428571428572*bmag_inv[1]*nuVtSqSum[2]*temp_diff[5]*m_+1.571428571428572*nuVtSqSum[1]*bmag_inv[2]*temp_diff[5]*m_+0.8944271909999161*bmag_inv[0]*nuVtSqSum[1]*temp_diff[5]*m_+0.8944271909999161*nuVtSqSum[0]*bmag_inv[1]*temp_diff[5]*m_+0.63887656499994*bmag_inv[2]*nuVtSqSum[2]*temp_diff[3]*m_+1.0*bmag_inv[0]*nuVtSqSum[2]*temp_diff[3]*m_+1.0*nuVtSqSum[0]*bmag_inv[2]*temp_diff[3]*m_+0.8944271909999161*bmag_inv[1]*nuVtSqSum[1]*temp_diff[3]*m_; 
  diff_incr[14] = 0.8944271909999161*bmag_inv[1]*nuVtSqSum[2]*temp_diff[18]*m_+0.8944271909999161*nuVtSqSum[1]*bmag_inv[2]*temp_diff[18]*m_+1.0*bmag_inv[0]*nuVtSqSum[1]*temp_diff[18]*m_+1.0*nuVtSqSum[0]*bmag_inv[1]*temp_diff[18]*m_+bmag_inv[2]*nuVtSqSum[2]*temp_diff[14]*m_+bmag_inv[1]*nuVtSqSum[1]*temp_diff[14]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[14]*m_; 
  diff_incr[15] = 1.571428571428571*bmag_inv[2]*nuVtSqSum[2]*temp_diff[15]*m_+0.8944271909999159*bmag_inv[0]*nuVtSqSum[2]*temp_diff[15]*m_+0.8944271909999159*nuVtSqSum[0]*bmag_inv[2]*temp_diff[15]*m_+1.8*bmag_inv[1]*nuVtSqSum[1]*temp_diff[15]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[15]*m_+0.8944271909999161*bmag_inv[1]*nuVtSqSum[2]*temp_diff[9]*m_+0.8944271909999161*nuVtSqSum[1]*bmag_inv[2]*temp_diff[9]*m_+1.0*bmag_inv[0]*nuVtSqSum[1]*temp_diff[9]*m_+1.0*nuVtSqSum[0]*bmag_inv[1]*temp_diff[9]*m_; 
  diff_incr[16] = 0.8944271909999161*bmag_inv[1]*nuVtSqSum[2]*temp_diff[19]*m_+0.8944271909999161*nuVtSqSum[1]*bmag_inv[2]*temp_diff[19]*m_+1.0*bmag_inv[0]*nuVtSqSum[1]*temp_diff[19]*m_+1.0*nuVtSqSum[0]*bmag_inv[1]*temp_diff[19]*m_+bmag_inv[2]*nuVtSqSum[2]*temp_diff[16]*m_+bmag_inv[1]*nuVtSqSum[1]*temp_diff[16]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[16]*m_; 
  diff_incr[17] = 2.142857142857143*bmag_inv[2]*nuVtSqSum[2]*temp_diff[17]*m_+0.6388765649999399*bmag_inv[0]*nuVtSqSum[2]*temp_diff[17]*m_+0.6388765649999399*nuVtSqSum[0]*bmag_inv[2]*temp_diff[17]*m_+1.571428571428571*bmag_inv[1]*nuVtSqSum[1]*temp_diff[17]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[17]*m_+1.571428571428571*bmag_inv[1]*nuVtSqSum[2]*temp_diff[10]*m_+1.571428571428571*nuVtSqSum[1]*bmag_inv[2]*temp_diff[10]*m_+0.8944271909999159*bmag_inv[0]*nuVtSqSum[1]*temp_diff[10]*m_+0.8944271909999159*nuVtSqSum[0]*bmag_inv[1]*temp_diff[10]*m_+0.6388765649999399*bmag_inv[2]*nuVtSqSum[2]*temp_diff[6]*m_+bmag_inv[0]*nuVtSqSum[2]*temp_diff[6]*m_+nuVtSqSum[0]*bmag_inv[2]*temp_diff[6]*m_+0.8944271909999159*bmag_inv[1]*nuVtSqSum[1]*temp_diff[6]*m_; 
  diff_incr[18] = 1.571428571428571*bmag_inv[2]*nuVtSqSum[2]*temp_diff[18]*m_+0.8944271909999159*bmag_inv[0]*nuVtSqSum[2]*temp_diff[18]*m_+0.8944271909999159*nuVtSqSum[0]*bmag_inv[2]*temp_diff[18]*m_+1.8*bmag_inv[1]*nuVtSqSum[1]*temp_diff[18]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[18]*m_+0.8944271909999161*bmag_inv[1]*nuVtSqSum[2]*temp_diff[14]*m_+0.8944271909999161*nuVtSqSum[1]*bmag_inv[2]*temp_diff[14]*m_+1.0*bmag_inv[0]*nuVtSqSum[1]*temp_diff[14]*m_+1.0*nuVtSqSum[0]*bmag_inv[1]*temp_diff[14]*m_; 
  diff_incr[19] = 1.571428571428571*bmag_inv[2]*nuVtSqSum[2]*temp_diff[19]*m_+0.8944271909999159*bmag_inv[0]*nuVtSqSum[2]*temp_diff[19]*m_+0.8944271909999159*nuVtSqSum[0]*bmag_inv[2]*temp_diff[19]*m_+1.8*bmag_inv[1]*nuVtSqSum[1]*temp_diff[19]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[19]*m_+0.8944271909999161*bmag_inv[1]*nuVtSqSum[2]*temp_diff[16]*m_+0.8944271909999161*nuVtSqSum[1]*bmag_inv[2]*temp_diff[16]*m_+1.0*bmag_inv[0]*nuVtSqSum[1]*temp_diff[16]*m_+1.0*nuVtSqSum[0]*bmag_inv[1]*temp_diff[16]*m_; 

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
} 
