#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH double lbo_gyrokinetic_diff_vol_3x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fin, double* GKYL_RESTRICT out) 
{ 
  // m_:        species mass.
  // bmag_inv:  1/(magnetic field magnitude). 
  // w[5]:      cell-center coordinates. 
  // dxv[5]:    cell spacing. 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // f:         input distribution function.
  // out:       incremented output 
  double rdv2[2]; 
  double rdvSq4[2]; 
  rdv2[0]   = 2.0/dxv[3]; 
  rdvSq4[0] = rdv2[0]*rdv2[0]; 
  rdv2[1]   = 2.0/dxv[4]; 
  rdvSq4[1] = rdv2[1]*rdv2[1]; 

  return fabs(0.3333333333333333*rdvSq4[1]*w[4]*(bmag_inv[7]*nuVtSqSum[7]+bmag_inv[6]*nuVtSqSum[6]+bmag_inv[5]*nuVtSqSum[5]+bmag_inv[4]*nuVtSqSum[4]+bmag_inv[3]*nuVtSqSum[3]+bmag_inv[2]*nuVtSqSum[2]+bmag_inv[1]*nuVtSqSum[1])*m_+nuVtSqSum[0]*(0.3333333333333333*bmag_inv[0]*rdvSq4[1]*w[4]*m_+0.4714045207910317*rdvSq4[0])); 

} 