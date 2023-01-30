#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH void lbo_gyrokinetic_diff_surfmu_1x2v_tensor_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // w[3]: cell-center coordinates. 
  // dxv[3]: cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[6]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fl/fc/fr: distribution function in cells 
  // out: incremented distribution function in cell 

  const double *nuVtSqSum = &nuPrimMomsSum[3];

  double rdvSq4 = 4.0/(dxv[2]*dxv[2]); 
  double incr[27] = {0.0}; 

  double diffFac[3] = {0.}; 
  diffFac[0] = 1.414213562373095*bmag_inv[2]*nuVtSqSum[2]*m_+1.414213562373095*bmag_inv[1]*nuVtSqSum[1]*m_+1.414213562373095*bmag_inv[0]*nuVtSqSum[0]*m_; 
  diffFac[1] = 1.264911064067352*bmag_inv[1]*nuVtSqSum[2]*m_+1.264911064067352*nuVtSqSum[1]*bmag_inv[2]*m_+1.414213562373095*bmag_inv[0]*nuVtSqSum[1]*m_+1.414213562373095*nuVtSqSum[0]*bmag_inv[1]*m_; 
  diffFac[2] = 0.9035079029052515*bmag_inv[2]*nuVtSqSum[2]*m_+1.414213562373095*bmag_inv[0]*nuVtSqSum[2]*m_+1.414213562373095*nuVtSqSum[0]*bmag_inv[2]*m_+1.264911064067352*bmag_inv[1]*nuVtSqSum[1]*m_; 

  double f_xx[27] = {0.0}; 
  f_xx[0] = 0.6708203932499369*w[2]*fr[9]+0.3354101966249685*dxv[2]*fr[9]+0.6708203932499369*w[2]*fl[9]-0.3354101966249685*dxv[2]*fl[9]-1.341640786499874*w[2]*fc[9]-1.190784930203603*w[2]*fr[3]-0.5953924651018015*dxv[2]*fr[3]+1.190784930203603*w[2]*fl[3]-0.5953924651018015*dxv[2]*fl[3]-1.190784930203603*dxv[2]*fc[3]+0.9375*fr[0]*w[2]+0.9375*fl[0]*w[2]-1.875*fc[0]*w[2]+0.46875*fr[0]*dxv[2]-0.46875*fl[0]*dxv[2]; 
  f_xx[1] = 0.6708203932499369*w[2]*fr[15]+0.3354101966249685*dxv[2]*fr[15]+0.6708203932499369*w[2]*fl[15]-0.3354101966249685*dxv[2]*fl[15]-1.341640786499874*w[2]*fc[15]-1.190784930203603*w[2]*fr[5]-0.5953924651018015*dxv[2]*fr[5]+1.190784930203603*w[2]*fl[5]-0.5953924651018015*dxv[2]*fl[5]-1.190784930203603*dxv[2]*fc[5]+0.9375*fr[1]*w[2]+0.9375*fl[1]*w[2]-1.875*fc[1]*w[2]+0.46875*fr[1]*dxv[2]-0.46875*fl[1]*dxv[2]; 
  f_xx[2] = 0.6708203932499369*w[2]*fr[16]+0.3354101966249685*dxv[2]*fr[16]+0.6708203932499369*w[2]*fl[16]-0.3354101966249685*dxv[2]*fl[16]-1.341640786499874*w[2]*fc[16]-1.190784930203603*w[2]*fr[6]-0.5953924651018015*dxv[2]*fr[6]+1.190784930203603*w[2]*fl[6]-0.5953924651018015*dxv[2]*fl[6]-1.190784930203603*dxv[2]*fc[6]+0.9375*fr[2]*w[2]+0.9375*fl[2]*w[2]-1.875*fc[2]*w[2]+0.46875*dxv[2]*fr[2]-0.46875*dxv[2]*fl[2]; 
  f_xx[3] = 0.7382874503707888*w[2]*fr[9]+0.3691437251853944*dxv[2]*fr[9]-0.7382874503707888*w[2]*fl[9]+0.3691437251853944*dxv[2]*fl[9]-1.585502557353661*dxv[2]*fc[9]-1.453125*w[2]*fr[3]-0.7265625*dxv[2]*fr[3]-1.453125*w[2]*fl[3]+0.7265625*dxv[2]*fl[3]-5.34375*w[2]*fc[3]+1.190784930203603*fr[0]*w[2]-1.190784930203603*fl[0]*w[2]+0.5953924651018015*fr[0]*dxv[2]+0.5953924651018015*fl[0]*dxv[2]-1.190784930203603*fc[0]*dxv[2]; 
  f_xx[4] = 0.6708203932499369*w[2]*fr[19]+0.3354101966249685*dxv[2]*fr[19]+0.6708203932499369*w[2]*fl[19]-0.3354101966249685*dxv[2]*fl[19]-1.341640786499874*w[2]*fc[19]-1.190784930203603*w[2]*fr[10]-0.5953924651018015*dxv[2]*fr[10]+1.190784930203603*w[2]*fl[10]-0.5953924651018015*dxv[2]*fl[10]-1.190784930203603*dxv[2]*fc[10]+0.9375*w[2]*fr[4]+0.46875*dxv[2]*fr[4]+0.9375*w[2]*fl[4]-0.46875*dxv[2]*fl[4]-1.875*w[2]*fc[4]; 
  f_xx[5] = 0.7382874503707888*w[2]*fr[15]+0.3691437251853944*dxv[2]*fr[15]-0.7382874503707888*w[2]*fl[15]+0.3691437251853944*dxv[2]*fl[15]-1.585502557353661*dxv[2]*fc[15]-1.453125*w[2]*fr[5]-0.7265625*dxv[2]*fr[5]-1.453125*w[2]*fl[5]+0.7265625*dxv[2]*fl[5]-5.34375*w[2]*fc[5]+1.190784930203603*fr[1]*w[2]-1.190784930203603*fl[1]*w[2]+0.5953924651018015*fr[1]*dxv[2]+0.5953924651018015*fl[1]*dxv[2]-1.190784930203603*fc[1]*dxv[2]; 
  f_xx[6] = 0.7382874503707888*w[2]*fr[16]+0.3691437251853944*dxv[2]*fr[16]-0.7382874503707888*w[2]*fl[16]+0.3691437251853944*dxv[2]*fl[16]-1.585502557353661*dxv[2]*fc[16]-1.453125*w[2]*fr[6]-0.7265625*dxv[2]*fr[6]-1.453125*w[2]*fl[6]+0.7265625*dxv[2]*fl[6]-5.34375*w[2]*fc[6]+1.190784930203603*fr[2]*w[2]-1.190784930203603*fl[2]*w[2]+0.5953924651018015*dxv[2]*fr[2]+0.5953924651018015*dxv[2]*fl[2]-1.190784930203603*dxv[2]*fc[2]; 
  f_xx[7] = 0.6708203932499369*w[2]*fr[21]+0.3354101966249685*dxv[2]*fr[21]+0.6708203932499369*w[2]*fl[21]-0.3354101966249685*dxv[2]*fl[21]-1.341640786499874*w[2]*fc[21]-1.190784930203603*w[2]*fr[13]-0.5953924651018015*dxv[2]*fr[13]+1.190784930203603*w[2]*fl[13]-0.5953924651018015*dxv[2]*fl[13]-1.190784930203603*dxv[2]*fc[13]+0.9375*w[2]*fr[7]+0.46875*dxv[2]*fr[7]+0.9375*w[2]*fl[7]-0.46875*dxv[2]*fl[7]-1.875*w[2]*fc[7]; 
  f_xx[8] = 0.6708203932499369*w[2]*fr[22]+0.3354101966249685*dxv[2]*fr[22]+0.6708203932499369*w[2]*fl[22]-0.3354101966249685*dxv[2]*fl[22]-1.341640786499874*w[2]*fc[22]-1.190784930203603*w[2]*fr[14]-0.5953924651018015*dxv[2]*fr[14]+1.190784930203603*w[2]*fl[14]-0.5953924651018015*dxv[2]*fl[14]-1.190784930203603*dxv[2]*fc[14]+0.9375*w[2]*fr[8]+0.46875*dxv[2]*fr[8]+0.9375*w[2]*fl[8]-0.46875*dxv[2]*fl[8]-1.875*w[2]*fc[8]; 
  f_xx[9] = (-0.140625*w[2]*fr[9])-0.0703125*dxv[2]*fr[9]-0.140625*w[2]*fl[9]+0.0703125*dxv[2]*fl[9]-6.28125*w[2]*fc[9]-0.3025768239224545*w[2]*fr[3]-0.1512884119612272*dxv[2]*fr[3]+0.3025768239224545*w[2]*fl[3]-0.1512884119612272*dxv[2]*fl[3]-1.149791930905327*dxv[2]*fc[3]+0.4192627457812106*fr[0]*w[2]+0.4192627457812106*fl[0]*w[2]-0.8385254915624212*fc[0]*w[2]+0.2096313728906053*fr[0]*dxv[2]-0.2096313728906053*fl[0]*dxv[2]; 
  f_xx[10] = 0.7382874503707888*w[2]*fr[19]+0.3691437251853944*dxv[2]*fr[19]-0.7382874503707888*w[2]*fl[19]+0.3691437251853944*dxv[2]*fl[19]-1.585502557353661*dxv[2]*fc[19]-1.453125*w[2]*fr[10]-0.7265625*dxv[2]*fr[10]-1.453125*w[2]*fl[10]+0.7265625*dxv[2]*fl[10]-5.34375*w[2]*fc[10]+1.190784930203603*w[2]*fr[4]+0.5953924651018015*dxv[2]*fr[4]-1.190784930203603*w[2]*fl[4]+0.5953924651018015*dxv[2]*fl[4]-1.190784930203603*dxv[2]*fc[4]; 
  f_xx[11] = 0.6708203932499369*w[2]*fr[24]+0.3354101966249685*dxv[2]*fr[24]+0.6708203932499369*w[2]*fl[24]-0.3354101966249685*dxv[2]*fl[24]-1.341640786499874*w[2]*fc[24]-1.190784930203603*w[2]*fr[17]-0.5953924651018015*dxv[2]*fr[17]+1.190784930203603*w[2]*fl[17]-0.5953924651018015*dxv[2]*fl[17]-1.190784930203603*dxv[2]*fc[17]+0.9375*w[2]*fr[11]+0.46875*dxv[2]*fr[11]+0.9375*w[2]*fl[11]-0.46875*dxv[2]*fl[11]-1.875*w[2]*fc[11]; 
  f_xx[12] = 0.6708203932499369*w[2]*fr[25]+0.3354101966249685*dxv[2]*fr[25]+0.6708203932499369*w[2]*fl[25]-0.3354101966249685*dxv[2]*fl[25]-1.341640786499874*w[2]*fc[25]-1.190784930203603*w[2]*fr[18]-0.5953924651018015*dxv[2]*fr[18]+1.190784930203603*w[2]*fl[18]-0.5953924651018015*dxv[2]*fl[18]-1.190784930203603*dxv[2]*fc[18]+0.9375*w[2]*fr[12]+0.46875*dxv[2]*fr[12]+0.9375*w[2]*fl[12]-0.46875*dxv[2]*fl[12]-1.875*w[2]*fc[12]; 
  f_xx[13] = 0.7382874503707888*w[2]*fr[21]+0.3691437251853944*dxv[2]*fr[21]-0.7382874503707888*w[2]*fl[21]+0.3691437251853944*dxv[2]*fl[21]-1.585502557353661*dxv[2]*fc[21]-1.453125*w[2]*fr[13]-0.7265625*dxv[2]*fr[13]-1.453125*w[2]*fl[13]+0.7265625*dxv[2]*fl[13]-5.34375*w[2]*fc[13]+1.190784930203603*w[2]*fr[7]+0.5953924651018015*dxv[2]*fr[7]-1.190784930203603*w[2]*fl[7]+0.5953924651018015*dxv[2]*fl[7]-1.190784930203603*dxv[2]*fc[7]; 
  f_xx[14] = 0.7382874503707888*w[2]*fr[22]+0.3691437251853944*dxv[2]*fr[22]-0.7382874503707888*w[2]*fl[22]+0.3691437251853944*dxv[2]*fl[22]-1.585502557353661*dxv[2]*fc[22]-1.453125*w[2]*fr[14]-0.7265625*dxv[2]*fr[14]-1.453125*w[2]*fl[14]+0.7265625*dxv[2]*fl[14]-5.34375*w[2]*fc[14]+1.190784930203603*w[2]*fr[8]+0.5953924651018015*dxv[2]*fr[8]-1.190784930203603*w[2]*fl[8]+0.5953924651018015*dxv[2]*fl[8]-1.190784930203603*dxv[2]*fc[8]; 
  f_xx[15] = (-0.140625*w[2]*fr[15])-0.0703125*dxv[2]*fr[15]-0.140625*w[2]*fl[15]+0.0703125*dxv[2]*fl[15]-6.28125*w[2]*fc[15]-0.3025768239224544*w[2]*fr[5]-0.1512884119612272*dxv[2]*fr[5]+0.3025768239224544*w[2]*fl[5]-0.1512884119612272*dxv[2]*fl[5]-1.149791930905327*dxv[2]*fc[5]+0.4192627457812105*fr[1]*w[2]+0.4192627457812105*fl[1]*w[2]-0.8385254915624211*fc[1]*w[2]+0.2096313728906053*fr[1]*dxv[2]-0.2096313728906053*fl[1]*dxv[2]; 
  f_xx[16] = (-0.140625*w[2]*fr[16])-0.0703125*dxv[2]*fr[16]-0.140625*w[2]*fl[16]+0.0703125*dxv[2]*fl[16]-6.28125*w[2]*fc[16]-0.3025768239224544*w[2]*fr[6]-0.1512884119612272*dxv[2]*fr[6]+0.3025768239224544*w[2]*fl[6]-0.1512884119612272*dxv[2]*fl[6]-1.149791930905327*dxv[2]*fc[6]+0.4192627457812105*fr[2]*w[2]+0.4192627457812105*fl[2]*w[2]-0.8385254915624211*fc[2]*w[2]+0.2096313728906053*dxv[2]*fr[2]-0.2096313728906053*dxv[2]*fl[2]; 
  f_xx[17] = 0.7382874503707888*w[2]*fr[24]+0.3691437251853944*dxv[2]*fr[24]-0.7382874503707888*w[2]*fl[24]+0.3691437251853944*dxv[2]*fl[24]-1.585502557353661*dxv[2]*fc[24]-1.453125*w[2]*fr[17]-0.7265625*dxv[2]*fr[17]-1.453125*w[2]*fl[17]+0.7265625*dxv[2]*fl[17]-5.34375*w[2]*fc[17]+1.190784930203603*w[2]*fr[11]+0.5953924651018015*dxv[2]*fr[11]-1.190784930203603*w[2]*fl[11]+0.5953924651018015*dxv[2]*fl[11]-1.190784930203603*dxv[2]*fc[11]; 
  f_xx[18] = 0.7382874503707888*w[2]*fr[25]+0.3691437251853944*dxv[2]*fr[25]-0.7382874503707888*w[2]*fl[25]+0.3691437251853944*dxv[2]*fl[25]-1.585502557353661*dxv[2]*fc[25]-1.453125*w[2]*fr[18]-0.7265625*dxv[2]*fr[18]-1.453125*w[2]*fl[18]+0.7265625*dxv[2]*fl[18]-5.34375*w[2]*fc[18]+1.190784930203603*w[2]*fr[12]+0.5953924651018015*dxv[2]*fr[12]-1.190784930203603*w[2]*fl[12]+0.5953924651018015*dxv[2]*fl[12]-1.190784930203603*dxv[2]*fc[12]; 
  f_xx[19] = (-0.140625*w[2]*fr[19])-0.0703125*dxv[2]*fr[19]-0.140625*w[2]*fl[19]+0.0703125*dxv[2]*fl[19]-6.28125*w[2]*fc[19]-0.3025768239224545*w[2]*fr[10]-0.1512884119612272*dxv[2]*fr[10]+0.3025768239224545*w[2]*fl[10]-0.1512884119612272*dxv[2]*fl[10]-1.149791930905327*dxv[2]*fc[10]+0.4192627457812106*w[2]*fr[4]+0.2096313728906053*dxv[2]*fr[4]+0.4192627457812106*w[2]*fl[4]-0.2096313728906053*dxv[2]*fl[4]-0.8385254915624212*w[2]*fc[4]; 
  f_xx[20] = 0.6708203932499369*w[2]*fr[26]+0.3354101966249685*dxv[2]*fr[26]+0.6708203932499369*w[2]*fl[26]-0.3354101966249685*dxv[2]*fl[26]-1.341640786499874*w[2]*fc[26]-1.190784930203603*w[2]*fr[23]-0.5953924651018015*dxv[2]*fr[23]+1.190784930203603*w[2]*fl[23]-0.5953924651018015*dxv[2]*fl[23]-1.190784930203603*dxv[2]*fc[23]+0.9375*w[2]*fr[20]+0.46875*dxv[2]*fr[20]+0.9375*w[2]*fl[20]-0.46875*dxv[2]*fl[20]-1.875*w[2]*fc[20]; 
  f_xx[21] = (-0.140625*w[2]*fr[21])-0.0703125*dxv[2]*fr[21]-0.140625*w[2]*fl[21]+0.0703125*dxv[2]*fl[21]-6.28125*w[2]*fc[21]-0.3025768239224544*w[2]*fr[13]-0.1512884119612272*dxv[2]*fr[13]+0.3025768239224544*w[2]*fl[13]-0.1512884119612272*dxv[2]*fl[13]-1.149791930905327*dxv[2]*fc[13]+0.4192627457812106*w[2]*fr[7]+0.2096313728906053*dxv[2]*fr[7]+0.4192627457812106*w[2]*fl[7]-0.2096313728906053*dxv[2]*fl[7]-0.8385254915624212*w[2]*fc[7]; 
  f_xx[22] = (-0.140625*w[2]*fr[22])-0.0703125*dxv[2]*fr[22]-0.140625*w[2]*fl[22]+0.0703125*dxv[2]*fl[22]-6.28125*w[2]*fc[22]-0.3025768239224544*w[2]*fr[14]-0.1512884119612272*dxv[2]*fr[14]+0.3025768239224544*w[2]*fl[14]-0.1512884119612272*dxv[2]*fl[14]-1.149791930905327*dxv[2]*fc[14]+0.4192627457812106*w[2]*fr[8]+0.2096313728906053*dxv[2]*fr[8]+0.4192627457812106*w[2]*fl[8]-0.2096313728906053*dxv[2]*fl[8]-0.8385254915624212*w[2]*fc[8]; 
  f_xx[23] = 0.7382874503707888*w[2]*fr[26]+0.3691437251853944*dxv[2]*fr[26]-0.7382874503707888*w[2]*fl[26]+0.3691437251853944*dxv[2]*fl[26]-1.585502557353661*dxv[2]*fc[26]-1.453125*w[2]*fr[23]-0.7265625*dxv[2]*fr[23]-1.453125*w[2]*fl[23]+0.7265625*dxv[2]*fl[23]-5.34375*w[2]*fc[23]+1.190784930203603*w[2]*fr[20]+0.5953924651018015*dxv[2]*fr[20]-1.190784930203603*w[2]*fl[20]+0.5953924651018015*dxv[2]*fl[20]-1.190784930203603*dxv[2]*fc[20]; 
  f_xx[24] = (-0.140625*w[2]*fr[24])-0.0703125*dxv[2]*fr[24]-0.140625*w[2]*fl[24]+0.0703125*dxv[2]*fl[24]-6.28125*w[2]*fc[24]-0.3025768239224545*w[2]*fr[17]-0.1512884119612272*dxv[2]*fr[17]+0.3025768239224545*w[2]*fl[17]-0.1512884119612272*dxv[2]*fl[17]-1.149791930905327*dxv[2]*fc[17]+0.4192627457812105*w[2]*fr[11]+0.2096313728906053*dxv[2]*fr[11]+0.4192627457812105*w[2]*fl[11]-0.2096313728906053*dxv[2]*fl[11]-0.8385254915624211*w[2]*fc[11]; 
  f_xx[25] = (-0.140625*w[2]*fr[25])-0.0703125*dxv[2]*fr[25]-0.140625*w[2]*fl[25]+0.0703125*dxv[2]*fl[25]-6.28125*w[2]*fc[25]-0.3025768239224545*w[2]*fr[18]-0.1512884119612272*dxv[2]*fr[18]+0.3025768239224545*w[2]*fl[18]-0.1512884119612272*dxv[2]*fl[18]-1.149791930905327*dxv[2]*fc[18]+0.4192627457812105*w[2]*fr[12]+0.2096313728906053*dxv[2]*fr[12]+0.4192627457812105*w[2]*fl[12]-0.2096313728906053*dxv[2]*fl[12]-0.8385254915624211*w[2]*fc[12]; 
  f_xx[26] = (-0.140625*w[2]*fr[26])-0.0703125*dxv[2]*fr[26]-0.140625*w[2]*fl[26]+0.0703125*dxv[2]*fl[26]-6.28125*w[2]*fc[26]-0.3025768239224545*w[2]*fr[23]-0.1512884119612272*dxv[2]*fr[23]+0.3025768239224545*w[2]*fl[23]-0.1512884119612272*dxv[2]*fl[23]-1.149791930905327*dxv[2]*fc[23]+0.4192627457812106*w[2]*fr[20]+0.2096313728906053*dxv[2]*fr[20]+0.4192627457812106*w[2]*fl[20]-0.2096313728906053*dxv[2]*fl[20]-0.8385254915624212*w[2]*fc[20]; 

  incr[0] = 0.7071067811865475*diffFac[2]*f_xx[7]+0.7071067811865475*diffFac[1]*f_xx[1]+0.7071067811865475*diffFac[0]*f_xx[0]; 
  incr[1] = 0.6324555320336759*diffFac[1]*f_xx[7]+0.6324555320336759*f_xx[1]*diffFac[2]+0.7071067811865475*diffFac[0]*f_xx[1]+0.7071067811865475*f_xx[0]*diffFac[1]; 
  incr[2] = 0.7071067811865475*diffFac[2]*f_xx[11]+0.7071067811865475*diffFac[1]*f_xx[4]+0.7071067811865475*diffFac[0]*f_xx[2]; 
  incr[3] = 0.7071067811865475*diffFac[2]*f_xx[13]+0.7071067811865475*diffFac[1]*f_xx[5]+0.7071067811865475*diffFac[0]*f_xx[3]; 
  incr[4] = 0.632455532033676*diffFac[1]*f_xx[11]+0.6324555320336759*diffFac[2]*f_xx[4]+0.7071067811865475*diffFac[0]*f_xx[4]+0.7071067811865475*diffFac[1]*f_xx[2]; 
  incr[5] = 0.632455532033676*diffFac[1]*f_xx[13]+0.6324555320336759*diffFac[2]*f_xx[5]+0.7071067811865475*diffFac[0]*f_xx[5]+0.7071067811865475*diffFac[1]*f_xx[3]; 
  incr[6] = 0.7071067811865475*diffFac[2]*f_xx[17]+0.7071067811865475*diffFac[1]*f_xx[10]+0.7071067811865475*diffFac[0]*f_xx[6]; 
  incr[7] = 0.4517539514526256*diffFac[2]*f_xx[7]+0.7071067811865475*diffFac[0]*f_xx[7]+0.7071067811865475*f_xx[0]*diffFac[2]+0.6324555320336759*diffFac[1]*f_xx[1]; 
  incr[8] = 0.7071067811865475*diffFac[2]*f_xx[20]+0.7071067811865475*diffFac[1]*f_xx[12]+0.7071067811865475*diffFac[0]*f_xx[8]; 
  incr[9] = 0.7071067811865475*diffFac[2]*f_xx[21]+0.7071067811865475*diffFac[1]*f_xx[15]+0.7071067811865475*diffFac[0]*f_xx[9]; 
  incr[10] = 0.6324555320336759*diffFac[1]*f_xx[17]+0.6324555320336759*diffFac[2]*f_xx[10]+0.7071067811865475*diffFac[0]*f_xx[10]+0.7071067811865475*diffFac[1]*f_xx[6]; 
  incr[11] = 0.4517539514526256*diffFac[2]*f_xx[11]+0.7071067811865475*diffFac[0]*f_xx[11]+0.632455532033676*diffFac[1]*f_xx[4]+0.7071067811865475*diffFac[2]*f_xx[2]; 
  incr[12] = 0.632455532033676*diffFac[1]*f_xx[20]+0.6324555320336759*diffFac[2]*f_xx[12]+0.7071067811865475*diffFac[0]*f_xx[12]+0.7071067811865475*diffFac[1]*f_xx[8]; 
  incr[13] = 0.4517539514526256*diffFac[2]*f_xx[13]+0.7071067811865475*diffFac[0]*f_xx[13]+0.632455532033676*diffFac[1]*f_xx[5]+0.7071067811865475*diffFac[2]*f_xx[3]; 
  incr[14] = 0.7071067811865475*diffFac[2]*f_xx[23]+0.7071067811865475*diffFac[1]*f_xx[18]+0.7071067811865475*diffFac[0]*f_xx[14]; 
  incr[15] = 0.632455532033676*diffFac[1]*f_xx[21]+0.6324555320336759*diffFac[2]*f_xx[15]+0.7071067811865475*diffFac[0]*f_xx[15]+0.7071067811865475*diffFac[1]*f_xx[9]; 
  incr[16] = 0.7071067811865475*diffFac[2]*f_xx[24]+0.7071067811865475*diffFac[1]*f_xx[19]+0.7071067811865475*diffFac[0]*f_xx[16]; 
  incr[17] = 0.4517539514526256*diffFac[2]*f_xx[17]+0.7071067811865475*diffFac[0]*f_xx[17]+0.6324555320336759*diffFac[1]*f_xx[10]+0.7071067811865475*diffFac[2]*f_xx[6]; 
  incr[18] = 0.6324555320336759*diffFac[1]*f_xx[23]+0.6324555320336759*diffFac[2]*f_xx[18]+0.7071067811865475*diffFac[0]*f_xx[18]+0.7071067811865475*diffFac[1]*f_xx[14]; 
  incr[19] = 0.6324555320336759*diffFac[1]*f_xx[24]+0.6324555320336759*diffFac[2]*f_xx[19]+0.7071067811865475*diffFac[0]*f_xx[19]+0.7071067811865475*diffFac[1]*f_xx[16]; 
  incr[20] = 0.4517539514526256*diffFac[2]*f_xx[20]+0.7071067811865475*diffFac[0]*f_xx[20]+0.632455532033676*diffFac[1]*f_xx[12]+0.7071067811865475*diffFac[2]*f_xx[8]; 
  incr[21] = 0.4517539514526256*diffFac[2]*f_xx[21]+0.7071067811865475*diffFac[0]*f_xx[21]+0.632455532033676*diffFac[1]*f_xx[15]+0.7071067811865475*diffFac[2]*f_xx[9]; 
  incr[22] = 0.7071067811865475*diffFac[2]*f_xx[26]+0.7071067811865475*diffFac[1]*f_xx[25]+0.7071067811865475*diffFac[0]*f_xx[22]; 
  incr[23] = 0.4517539514526256*diffFac[2]*f_xx[23]+0.7071067811865475*diffFac[0]*f_xx[23]+0.6324555320336759*diffFac[1]*f_xx[18]+0.7071067811865475*diffFac[2]*f_xx[14]; 
  incr[24] = 0.4517539514526256*diffFac[2]*f_xx[24]+0.7071067811865475*diffFac[0]*f_xx[24]+0.6324555320336759*diffFac[1]*f_xx[19]+0.7071067811865475*diffFac[2]*f_xx[16]; 
  incr[25] = 0.6324555320336759*diffFac[1]*f_xx[26]+0.6324555320336759*diffFac[2]*f_xx[25]+0.7071067811865475*diffFac[0]*f_xx[25]+0.7071067811865475*diffFac[1]*f_xx[22]; 
  incr[26] = 0.4517539514526256*diffFac[2]*f_xx[26]+0.7071067811865475*diffFac[0]*f_xx[26]+0.6324555320336759*diffFac[1]*f_xx[25]+0.7071067811865475*diffFac[2]*f_xx[22]; 

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
} 
