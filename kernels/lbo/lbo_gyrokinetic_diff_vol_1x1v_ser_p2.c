#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH double lbo_gyrokinetic_diff_vol_1x1v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fin, double* GKYL_RESTRICT out) 
{ 
  // m_:        species mass.
  // bmag_inv:  1/(magnetic field magnitude). 
  // w[2]:      cell-center coordinates. 
  // dxv[2]:    cell spacing. 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // f:         input distribution function.
  // out:       incremented output 
  double rdv2[1]; 
  double rdvSq4[1]; 
  rdv2[0]   = 2.0/dxv[1]; 
  rdvSq4[0] = rdv2[0]*rdv2[0]; 

  double facDiffVpar[3];
  // Expand nuVtSqSum in phase basis.
  facDiffVpar[0] = nuVtSqSum[0]; 
  facDiffVpar[2] = nuVtSqSum[2]; 

  return fabs(rdvSq4[0]*(1.272792206135785*facDiffVpar[0]-1.42302494707577*facDiffVpar[2])); 

} 