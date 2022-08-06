#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH void lbo_gyrokinetic_drag_boundary_surfmu_1x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[3]:     cell-center coordinates. 
  // dxv[3]:   cell spacing. 
  // m_:        species mass.
  // bmag_inv:  1/(magnetic field magnitude). 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum[4]:sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[2]; 

  double alphaDrSurf[6] = {0.0}; 
  double Ghat[6] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = nuSum[0]*(2.828427124746191*w[2]+1.414213562373095*dxv[2]); 
  alphaDrSurf[1] = nuSum[1]*(2.828427124746191*w[2]+1.414213562373095*dxv[2]); 

  Ghat[0] = -0.25*(2.449489742783178*(alphaDrSurf[1]*fEdge[5]+alphaDrSurf[0]*fEdge[3])-1.414213562373095*(alphaDrSurf[1]*fEdge[1]+alphaDrSurf[0]*fEdge[0])); 
  Ghat[1] = -0.25*(2.449489742783178*(alphaDrSurf[0]*fEdge[5]+alphaDrSurf[1]*fEdge[3])-1.414213562373095*(alphaDrSurf[0]*fEdge[1]+fEdge[0]*alphaDrSurf[1])); 
  Ghat[2] = -0.25*(2.449489742783178*(alphaDrSurf[1]*fEdge[7]+alphaDrSurf[0]*fEdge[6])-1.414213562373095*(alphaDrSurf[1]*fEdge[4]+alphaDrSurf[0]*fEdge[2])); 
  Ghat[3] = -0.25*(2.449489742783178*(alphaDrSurf[0]*fEdge[7]+alphaDrSurf[1]*fEdge[6])-1.414213562373095*(alphaDrSurf[0]*fEdge[4]+alphaDrSurf[1]*fEdge[2])); 
  Ghat[4] = -0.01666666666666667*(36.74234614174767*alphaDrSurf[1]*fEdge[11]+36.74234614174768*alphaDrSurf[0]*fEdge[10]-21.21320343559643*alphaDrSurf[1]*fEdge[9]-21.21320343559643*alphaDrSurf[0]*fEdge[8]); 
  Ghat[5] = -0.01666666666666667*(36.74234614174768*alphaDrSurf[0]*fEdge[11]+36.74234614174767*alphaDrSurf[1]*fEdge[10]-21.21320343559643*alphaDrSurf[0]*fEdge[9]-21.21320343559643*alphaDrSurf[1]*fEdge[8]); 

  out[0] += 0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += 1.224744871391589*Ghat[0]*rdv2; 
  out[4] += 0.7071067811865475*Ghat[3]*rdv2; 
  out[5] += 1.224744871391589*Ghat[1]*rdv2; 
  out[6] += 1.224744871391589*Ghat[2]*rdv2; 
  out[7] += 1.224744871391589*Ghat[3]*rdv2; 
  out[8] += 0.7071067811865475*Ghat[4]*rdv2; 
  out[9] += 0.7071067811865475*Ghat[5]*rdv2; 
  out[10] += 1.224744871391589*Ghat[4]*rdv2; 
  out[11] += 1.224744871391589*Ghat[5]*rdv2; 

  } else { 

  alphaDrSurf[0] = nuSum[0]*(2.828427124746191*w[2]-1.414213562373095*dxv[2]); 
  alphaDrSurf[1] = nuSum[1]*(2.828427124746191*w[2]-1.414213562373095*dxv[2]); 

  Ghat[0] = -0.25*(2.449489742783178*(alphaDrSurf[1]*fSkin[5]+alphaDrSurf[0]*fSkin[3])-1.414213562373095*(alphaDrSurf[1]*fSkin[1]+alphaDrSurf[0]*fSkin[0])); 
  Ghat[1] = -0.25*(2.449489742783178*(alphaDrSurf[0]*fSkin[5]+alphaDrSurf[1]*fSkin[3])-1.414213562373095*(alphaDrSurf[0]*fSkin[1]+fSkin[0]*alphaDrSurf[1])); 
  Ghat[2] = -0.25*(2.449489742783178*(alphaDrSurf[1]*fSkin[7]+alphaDrSurf[0]*fSkin[6])-1.414213562373095*(alphaDrSurf[1]*fSkin[4]+alphaDrSurf[0]*fSkin[2])); 
  Ghat[3] = -0.25*(2.449489742783178*(alphaDrSurf[0]*fSkin[7]+alphaDrSurf[1]*fSkin[6])-1.414213562373095*(alphaDrSurf[0]*fSkin[4]+alphaDrSurf[1]*fSkin[2])); 
  Ghat[4] = -0.01666666666666667*(36.74234614174767*alphaDrSurf[1]*fSkin[11]+36.74234614174768*alphaDrSurf[0]*fSkin[10]-21.21320343559643*alphaDrSurf[1]*fSkin[9]-21.21320343559643*alphaDrSurf[0]*fSkin[8]); 
  Ghat[5] = -0.01666666666666667*(36.74234614174768*alphaDrSurf[0]*fSkin[11]+36.74234614174767*alphaDrSurf[1]*fSkin[10]-21.21320343559643*alphaDrSurf[0]*fSkin[9]-21.21320343559643*alphaDrSurf[1]*fSkin[8]); 

  out[0] += -0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += -0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += -0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += 1.224744871391589*Ghat[0]*rdv2; 
  out[4] += -0.7071067811865475*Ghat[3]*rdv2; 
  out[5] += 1.224744871391589*Ghat[1]*rdv2; 
  out[6] += 1.224744871391589*Ghat[2]*rdv2; 
  out[7] += 1.224744871391589*Ghat[3]*rdv2; 
  out[8] += -0.7071067811865475*Ghat[4]*rdv2; 
  out[9] += -0.7071067811865475*Ghat[5]*rdv2; 
  out[10] += 1.224744871391589*Ghat[4]*rdv2; 
  out[11] += 1.224744871391589*Ghat[5]*rdv2; 

  } 
} 
