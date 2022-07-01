#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH void lbo_gyrokinetic_diff_boundary_surfmu_1x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  // w[3]:         Cell-center coordinates. 
  // dxv[3]:       Cell spacing. 
  // m_:           species mass.
  // bmag_inv:     1/(magnetic field magnitude). 
  // nuSum:        collisionalities added (self and cross species collisionalities). 
  // nuUSum[4]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]: sum of thermal speeds squared time their respective collisionalities. 
  // fskin/edge:   Distribution function in cells 
  // out:          Incremented distribution function in cell 
  double rdvSq4 = 4.0/(dxv[2]*dxv[2]); 

  double facDiff[2]; 
  // Expand diffusion coefficient in conf basis.
  facDiff[0] = 1.414213562373095*(bmag_inv[1]*nuVtSqSum[1]+bmag_inv[0]*nuVtSqSum[0])*m_; 
  facDiff[1] = 1.414213562373095*(bmag_inv[0]*nuVtSqSum[1]+nuVtSqSum[0]*bmag_inv[1])*m_; 

  double vol_incr[8] = {0.0};
  vol_incr[3] = 0.6123724356957944*facDiff[1]*fskin[1]*dxv[2]*rdvSq4+0.6123724356957944*facDiff[0]*fskin[0]*dxv[2]*rdvSq4; 
  vol_incr[5] = 0.6123724356957944*facDiff[0]*fskin[1]*dxv[2]*rdvSq4+0.6123724356957944*fskin[0]*facDiff[1]*dxv[2]*rdvSq4; 
  vol_incr[6] = 0.6123724356957944*facDiff[1]*dxv[2]*fskin[4]*rdvSq4+0.6123724356957944*facDiff[0]*dxv[2]*fskin[2]*rdvSq4; 
  vol_incr[7] = 0.6123724356957944*facDiff[0]*dxv[2]*fskin[4]*rdvSq4+0.6123724356957944*facDiff[1]*dxv[2]*fskin[2]*rdvSq4; 

  double temp_diff[8] = {0.0}; 
  double temp_edge[8] = {0.0}; 
  double diff_incr[8] = {0.0}; 
  double edge_incr[8] = {0.0}; 

  if (edge == -1) { 

  temp_diff[0] = (-0.5412658773652741*w[2]*fskin[3])-0.270632938682637*dxv[2]*fskin[3]-0.5412658773652741*w[2]*fedge[3]-0.270632938682637*dxv[2]*fedge[3]-0.5625*fskin[0]*w[2]+0.5625*fedge[0]*w[2]-0.28125*fskin[0]*dxv[2]+0.28125*fedge[0]*dxv[2]; 
  temp_diff[1] = (-0.5412658773652741*w[2]*fskin[5])-0.270632938682637*dxv[2]*fskin[5]-0.5412658773652741*w[2]*fedge[5]-0.270632938682637*dxv[2]*fedge[5]-0.5625*fskin[1]*w[2]+0.5625*fedge[1]*w[2]-0.28125*fskin[1]*dxv[2]+0.28125*fedge[1]*dxv[2]; 
  temp_diff[2] = (-0.5412658773652741*w[2]*fskin[6])-0.270632938682637*dxv[2]*fskin[6]-0.5412658773652741*w[2]*fedge[6]-0.270632938682637*dxv[2]*fedge[6]-0.5625*fskin[2]*w[2]+0.5625*fedge[2]*w[2]-0.28125*dxv[2]*fskin[2]+0.28125*dxv[2]*fedge[2]; 
  temp_diff[3] = (-1.4375*w[2]*fskin[3])-0.71875*dxv[2]*fskin[3]-0.4375*w[2]*fedge[3]-0.21875*dxv[2]*fedge[3]-1.407291281149713*fskin[0]*w[2]+0.5412658773652741*fedge[0]*w[2]-0.7036456405748563*fskin[0]*dxv[2]+0.270632938682637*fedge[0]*dxv[2]; 
  temp_diff[4] = (-0.5412658773652741*w[2]*fskin[7])-0.270632938682637*dxv[2]*fskin[7]-0.5412658773652741*w[2]*fedge[7]-0.270632938682637*dxv[2]*fedge[7]-0.5625*w[2]*fskin[4]-0.28125*dxv[2]*fskin[4]+0.5625*w[2]*fedge[4]+0.28125*dxv[2]*fedge[4]; 
  temp_diff[5] = (-1.4375*w[2]*fskin[5])-0.71875*dxv[2]*fskin[5]-0.4375*w[2]*fedge[5]-0.21875*dxv[2]*fedge[5]-1.407291281149713*fskin[1]*w[2]+0.5412658773652741*fedge[1]*w[2]-0.7036456405748563*fskin[1]*dxv[2]+0.270632938682637*fedge[1]*dxv[2]; 
  temp_diff[6] = (-1.4375*w[2]*fskin[6])-0.71875*dxv[2]*fskin[6]-0.4375*w[2]*fedge[6]-0.21875*dxv[2]*fedge[6]-1.407291281149713*fskin[2]*w[2]+0.5412658773652741*fedge[2]*w[2]-0.7036456405748563*dxv[2]*fskin[2]+0.270632938682637*dxv[2]*fedge[2]; 
  temp_diff[7] = (-1.4375*w[2]*fskin[7])-0.71875*dxv[2]*fskin[7]-0.4375*w[2]*fedge[7]-0.21875*dxv[2]*fedge[7]-1.407291281149713*w[2]*fskin[4]-0.7036456405748563*dxv[2]*fskin[4]+0.5412658773652741*w[2]*fedge[4]+0.270632938682637*dxv[2]*fedge[4]; 

  temp_edge[3] = (-1.5*w[2]*fskin[3])+0.75*dxv[2]*fskin[3]+0.8660254037844386*fskin[0]*w[2]-0.4330127018922193*fskin[0]*dxv[2]; 
  temp_edge[5] = (-1.5*w[2]*fskin[5])+0.75*dxv[2]*fskin[5]+0.8660254037844386*fskin[1]*w[2]-0.4330127018922193*fskin[1]*dxv[2]; 
  temp_edge[6] = (-1.5*w[2]*fskin[6])+0.75*dxv[2]*fskin[6]+0.8660254037844386*fskin[2]*w[2]-0.4330127018922193*dxv[2]*fskin[2]; 
  temp_edge[7] = (-1.5*w[2]*fskin[7])+0.75*dxv[2]*fskin[7]+0.8660254037844386*w[2]*fskin[4]-0.4330127018922193*dxv[2]*fskin[4]; 

  diff_incr[0] = bmag_inv[0]*nuVtSqSum[1]*temp_diff[1]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_diff[1]*m_+temp_diff[0]*bmag_inv[1]*nuVtSqSum[1]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[0]*m_; 
  diff_incr[1] = 1.8*bmag_inv[1]*nuVtSqSum[1]*temp_diff[1]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[1]*m_+bmag_inv[0]*temp_diff[0]*nuVtSqSum[1]*m_+nuVtSqSum[0]*temp_diff[0]*bmag_inv[1]*m_; 
  diff_incr[2] = bmag_inv[0]*nuVtSqSum[1]*temp_diff[4]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_diff[4]*m_+bmag_inv[1]*nuVtSqSum[1]*temp_diff[2]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[2]*m_; 
  diff_incr[3] = bmag_inv[0]*nuVtSqSum[1]*temp_diff[5]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_diff[5]*m_+bmag_inv[1]*nuVtSqSum[1]*temp_diff[3]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[3]*m_; 
  diff_incr[4] = 1.8*bmag_inv[1]*nuVtSqSum[1]*temp_diff[4]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[4]*m_+bmag_inv[0]*nuVtSqSum[1]*temp_diff[2]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_diff[2]*m_; 
  diff_incr[5] = 1.8*bmag_inv[1]*nuVtSqSum[1]*temp_diff[5]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[5]*m_+bmag_inv[0]*nuVtSqSum[1]*temp_diff[3]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_diff[3]*m_; 
  diff_incr[6] = bmag_inv[0]*nuVtSqSum[1]*temp_diff[7]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_diff[7]*m_+bmag_inv[1]*nuVtSqSum[1]*temp_diff[6]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[6]*m_; 
  diff_incr[7] = 1.8*bmag_inv[1]*nuVtSqSum[1]*temp_diff[7]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[7]*m_+bmag_inv[0]*nuVtSqSum[1]*temp_diff[6]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_diff[6]*m_; 

  edge_incr[0] = bmag_inv[0]*nuVtSqSum[1]*temp_edge[1]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_edge[1]*m_+temp_edge[0]*bmag_inv[1]*nuVtSqSum[1]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_edge[0]*m_; 
  edge_incr[1] = 1.8*bmag_inv[1]*nuVtSqSum[1]*temp_edge[1]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_edge[1]*m_+bmag_inv[0]*temp_edge[0]*nuVtSqSum[1]*m_+nuVtSqSum[0]*temp_edge[0]*bmag_inv[1]*m_; 
  edge_incr[2] = bmag_inv[0]*nuVtSqSum[1]*temp_edge[4]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_edge[4]*m_+bmag_inv[1]*nuVtSqSum[1]*temp_edge[2]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_edge[2]*m_; 
  edge_incr[3] = bmag_inv[0]*nuVtSqSum[1]*temp_edge[5]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_edge[5]*m_+bmag_inv[1]*nuVtSqSum[1]*temp_edge[3]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_edge[3]*m_; 
  edge_incr[4] = 1.8*bmag_inv[1]*nuVtSqSum[1]*temp_edge[4]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_edge[4]*m_+bmag_inv[0]*nuVtSqSum[1]*temp_edge[2]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_edge[2]*m_; 
  edge_incr[5] = 1.8*bmag_inv[1]*nuVtSqSum[1]*temp_edge[5]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_edge[5]*m_+bmag_inv[0]*nuVtSqSum[1]*temp_edge[3]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_edge[3]*m_; 
  edge_incr[6] = bmag_inv[0]*nuVtSqSum[1]*temp_edge[7]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_edge[7]*m_+bmag_inv[1]*nuVtSqSum[1]*temp_edge[6]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_edge[6]*m_; 
  edge_incr[7] = 1.8*bmag_inv[1]*nuVtSqSum[1]*temp_edge[7]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_edge[7]*m_+bmag_inv[0]*nuVtSqSum[1]*temp_edge[6]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_edge[6]*m_; 


  } else { 

  temp_diff[0] = 0.5412658773652741*w[2]*fskin[3]-0.270632938682637*dxv[2]*fskin[3]+0.5412658773652741*w[2]*fedge[3]-0.270632938682637*dxv[2]*fedge[3]-0.5625*fskin[0]*w[2]+0.5625*fedge[0]*w[2]+0.28125*fskin[0]*dxv[2]-0.28125*fedge[0]*dxv[2]; 
  temp_diff[1] = 0.5412658773652741*w[2]*fskin[5]-0.270632938682637*dxv[2]*fskin[5]+0.5412658773652741*w[2]*fedge[5]-0.270632938682637*dxv[2]*fedge[5]-0.5625*fskin[1]*w[2]+0.5625*fedge[1]*w[2]+0.28125*fskin[1]*dxv[2]-0.28125*fedge[1]*dxv[2]; 
  temp_diff[2] = 0.5412658773652741*w[2]*fskin[6]-0.270632938682637*dxv[2]*fskin[6]+0.5412658773652741*w[2]*fedge[6]-0.270632938682637*dxv[2]*fedge[6]-0.5625*fskin[2]*w[2]+0.5625*fedge[2]*w[2]+0.28125*dxv[2]*fskin[2]-0.28125*dxv[2]*fedge[2]; 
  temp_diff[3] = (-1.4375*w[2]*fskin[3])+0.71875*dxv[2]*fskin[3]-0.4375*w[2]*fedge[3]+0.21875*dxv[2]*fedge[3]+1.407291281149713*fskin[0]*w[2]-0.5412658773652741*fedge[0]*w[2]-0.7036456405748563*fskin[0]*dxv[2]+0.270632938682637*fedge[0]*dxv[2]; 
  temp_diff[4] = 0.5412658773652741*w[2]*fskin[7]-0.270632938682637*dxv[2]*fskin[7]+0.5412658773652741*w[2]*fedge[7]-0.270632938682637*dxv[2]*fedge[7]-0.5625*w[2]*fskin[4]+0.28125*dxv[2]*fskin[4]+0.5625*w[2]*fedge[4]-0.28125*dxv[2]*fedge[4]; 
  temp_diff[5] = (-1.4375*w[2]*fskin[5])+0.71875*dxv[2]*fskin[5]-0.4375*w[2]*fedge[5]+0.21875*dxv[2]*fedge[5]+1.407291281149713*fskin[1]*w[2]-0.5412658773652741*fedge[1]*w[2]-0.7036456405748563*fskin[1]*dxv[2]+0.270632938682637*fedge[1]*dxv[2]; 
  temp_diff[6] = (-1.4375*w[2]*fskin[6])+0.71875*dxv[2]*fskin[6]-0.4375*w[2]*fedge[6]+0.21875*dxv[2]*fedge[6]+1.407291281149713*fskin[2]*w[2]-0.5412658773652741*fedge[2]*w[2]-0.7036456405748563*dxv[2]*fskin[2]+0.270632938682637*dxv[2]*fedge[2]; 
  temp_diff[7] = (-1.4375*w[2]*fskin[7])+0.71875*dxv[2]*fskin[7]-0.4375*w[2]*fedge[7]+0.21875*dxv[2]*fedge[7]+1.407291281149713*w[2]*fskin[4]-0.7036456405748563*dxv[2]*fskin[4]-0.5412658773652741*w[2]*fedge[4]+0.270632938682637*dxv[2]*fedge[4]; 

  temp_edge[3] = (-1.5*w[2]*fskin[3])-0.75*dxv[2]*fskin[3]-0.8660254037844386*fskin[0]*w[2]-0.4330127018922193*fskin[0]*dxv[2]; 
  temp_edge[5] = (-1.5*w[2]*fskin[5])-0.75*dxv[2]*fskin[5]-0.8660254037844386*fskin[1]*w[2]-0.4330127018922193*fskin[1]*dxv[2]; 
  temp_edge[6] = (-1.5*w[2]*fskin[6])-0.75*dxv[2]*fskin[6]-0.8660254037844386*fskin[2]*w[2]-0.4330127018922193*dxv[2]*fskin[2]; 
  temp_edge[7] = (-1.5*w[2]*fskin[7])-0.75*dxv[2]*fskin[7]-0.8660254037844386*w[2]*fskin[4]-0.4330127018922193*dxv[2]*fskin[4]; 

  diff_incr[0] = bmag_inv[0]*nuVtSqSum[1]*temp_diff[1]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_diff[1]*m_+temp_diff[0]*bmag_inv[1]*nuVtSqSum[1]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[0]*m_; 
  diff_incr[1] = 1.8*bmag_inv[1]*nuVtSqSum[1]*temp_diff[1]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[1]*m_+bmag_inv[0]*temp_diff[0]*nuVtSqSum[1]*m_+nuVtSqSum[0]*temp_diff[0]*bmag_inv[1]*m_; 
  diff_incr[2] = bmag_inv[0]*nuVtSqSum[1]*temp_diff[4]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_diff[4]*m_+bmag_inv[1]*nuVtSqSum[1]*temp_diff[2]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[2]*m_; 
  diff_incr[3] = bmag_inv[0]*nuVtSqSum[1]*temp_diff[5]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_diff[5]*m_+bmag_inv[1]*nuVtSqSum[1]*temp_diff[3]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[3]*m_; 
  diff_incr[4] = 1.8*bmag_inv[1]*nuVtSqSum[1]*temp_diff[4]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[4]*m_+bmag_inv[0]*nuVtSqSum[1]*temp_diff[2]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_diff[2]*m_; 
  diff_incr[5] = 1.8*bmag_inv[1]*nuVtSqSum[1]*temp_diff[5]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[5]*m_+bmag_inv[0]*nuVtSqSum[1]*temp_diff[3]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_diff[3]*m_; 
  diff_incr[6] = bmag_inv[0]*nuVtSqSum[1]*temp_diff[7]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_diff[7]*m_+bmag_inv[1]*nuVtSqSum[1]*temp_diff[6]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[6]*m_; 
  diff_incr[7] = 1.8*bmag_inv[1]*nuVtSqSum[1]*temp_diff[7]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[7]*m_+bmag_inv[0]*nuVtSqSum[1]*temp_diff[6]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_diff[6]*m_; 

  edge_incr[0] = bmag_inv[0]*nuVtSqSum[1]*temp_edge[1]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_edge[1]*m_+temp_edge[0]*bmag_inv[1]*nuVtSqSum[1]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_edge[0]*m_; 
  edge_incr[1] = 1.8*bmag_inv[1]*nuVtSqSum[1]*temp_edge[1]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_edge[1]*m_+bmag_inv[0]*temp_edge[0]*nuVtSqSum[1]*m_+nuVtSqSum[0]*temp_edge[0]*bmag_inv[1]*m_; 
  edge_incr[2] = bmag_inv[0]*nuVtSqSum[1]*temp_edge[4]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_edge[4]*m_+bmag_inv[1]*nuVtSqSum[1]*temp_edge[2]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_edge[2]*m_; 
  edge_incr[3] = bmag_inv[0]*nuVtSqSum[1]*temp_edge[5]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_edge[5]*m_+bmag_inv[1]*nuVtSqSum[1]*temp_edge[3]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_edge[3]*m_; 
  edge_incr[4] = 1.8*bmag_inv[1]*nuVtSqSum[1]*temp_edge[4]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_edge[4]*m_+bmag_inv[0]*nuVtSqSum[1]*temp_edge[2]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_edge[2]*m_; 
  edge_incr[5] = 1.8*bmag_inv[1]*nuVtSqSum[1]*temp_edge[5]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_edge[5]*m_+bmag_inv[0]*nuVtSqSum[1]*temp_edge[3]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_edge[3]*m_; 
  edge_incr[6] = bmag_inv[0]*nuVtSqSum[1]*temp_edge[7]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_edge[7]*m_+bmag_inv[1]*nuVtSqSum[1]*temp_edge[6]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_edge[6]*m_; 
  edge_incr[7] = 1.8*bmag_inv[1]*nuVtSqSum[1]*temp_edge[7]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_edge[7]*m_+bmag_inv[0]*nuVtSqSum[1]*temp_edge[6]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_edge[6]*m_; 

  } 

  out[0] += edge_incr[0]*rdvSq4+diff_incr[0]*rdvSq4+vol_incr[0]; 
  out[1] += edge_incr[1]*rdvSq4+diff_incr[1]*rdvSq4+vol_incr[1]; 
  out[2] += edge_incr[2]*rdvSq4+diff_incr[2]*rdvSq4+vol_incr[2]; 
  out[3] += edge_incr[3]*rdvSq4+diff_incr[3]*rdvSq4+vol_incr[3]; 
  out[4] += edge_incr[4]*rdvSq4+diff_incr[4]*rdvSq4+vol_incr[4]; 
  out[5] += edge_incr[5]*rdvSq4+diff_incr[5]*rdvSq4+vol_incr[5]; 
  out[6] += edge_incr[6]*rdvSq4+diff_incr[6]*rdvSq4+vol_incr[6]; 
  out[7] += edge_incr[7]*rdvSq4+diff_incr[7]*rdvSq4+vol_incr[7]; 
} 
