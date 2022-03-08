#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_ser_4x_p1_surfx4_quad.h> 
#include <gkyl_basis_ser_4x_p1_upwind.h> 
GKYL_CU_DH void lbo_vlasov_drag_boundary_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[4]:         Cell-center coordinates. 
  // dxv[4]:       Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[8]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[4]: sum of thermal speeds squared time their respective collisionalities. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[3]; 

  const double *sumNuUy = &nuUSum[4]; 

  double alphaDrSurf[8] = {0.0}; 
  double fUpwindQuad[8] = {0.0};
  double fUpwind[8] = {0.0};;
  double drag_incr[8] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = 0.7071067811865475*(nuSum[0]*(2.0*w[3]+dxv[3])-2.0*sumNuUy[0]); 
  alphaDrSurf[1] = 0.7071067811865475*(nuSum[1]*(2.0*w[3]+dxv[3])-2.0*sumNuUy[1]); 
  alphaDrSurf[2] = 0.7071067811865475*(nuSum[2]*(2.0*w[3]+dxv[3])-2.0*sumNuUy[2]); 
  alphaDrSurf[4] = 0.7071067811865475*(2.0*nuSum[3]*w[3]-2.0*sumNuUy[3]+dxv[3]*nuSum[3]); 

  if (alphaDrSurf[4]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = ser_4x_p1_surfx4_quad_0_r(fSkin); 
    fUpwindQuad[4] = ser_4x_p1_surfx4_quad_4_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_4x_p1_surfx4_quad_0_l(fEdge); 
    fUpwindQuad[4] = ser_4x_p1_surfx4_quad_4_l(fEdge); 
  } 
  if ((-alphaDrSurf[4])-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[1] = ser_4x_p1_surfx4_quad_1_r(fSkin); 
    fUpwindQuad[5] = ser_4x_p1_surfx4_quad_5_r(fSkin); 
  } else { 
    fUpwindQuad[1] = ser_4x_p1_surfx4_quad_1_l(fEdge); 
    fUpwindQuad[5] = ser_4x_p1_surfx4_quad_5_l(fEdge); 
  } 
  if ((-alphaDrSurf[4])+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[2] = ser_4x_p1_surfx4_quad_2_r(fSkin); 
    fUpwindQuad[6] = ser_4x_p1_surfx4_quad_6_r(fSkin); 
  } else { 
    fUpwindQuad[2] = ser_4x_p1_surfx4_quad_2_l(fEdge); 
    fUpwindQuad[6] = ser_4x_p1_surfx4_quad_6_l(fEdge); 
  } 
  if (alphaDrSurf[4]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = ser_4x_p1_surfx4_quad_3_r(fSkin); 
    fUpwindQuad[7] = ser_4x_p1_surfx4_quad_7_r(fSkin); 
  } else { 
    fUpwindQuad[3] = ser_4x_p1_surfx4_quad_3_l(fEdge); 
    fUpwindQuad[7] = ser_4x_p1_surfx4_quad_7_l(fEdge); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_4x_p1_upwind(fUpwindQuad, fUpwind); 

  drag_incr[0] = 0.3535533905932737*alphaDrSurf[4]*fUpwind[4]+0.3535533905932737*alphaDrSurf[2]*fUpwind[2]+0.3535533905932737*alphaDrSurf[1]*fUpwind[1]+0.3535533905932737*alphaDrSurf[0]*fUpwind[0]; 
  drag_incr[1] = 0.3535533905932737*alphaDrSurf[2]*fUpwind[4]+0.3535533905932737*fUpwind[2]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[0]*fUpwind[1]+0.3535533905932737*fUpwind[0]*alphaDrSurf[1]; 
  drag_incr[2] = 0.3535533905932737*alphaDrSurf[1]*fUpwind[4]+0.3535533905932737*fUpwind[1]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[0]*fUpwind[2]+0.3535533905932737*fUpwind[0]*alphaDrSurf[2]; 
  drag_incr[3] = 0.3535533905932737*alphaDrSurf[4]*fUpwind[7]+0.3535533905932737*alphaDrSurf[2]*fUpwind[6]+0.3535533905932737*alphaDrSurf[1]*fUpwind[5]+0.3535533905932737*alphaDrSurf[0]*fUpwind[3]; 
  drag_incr[4] = 0.3535533905932737*alphaDrSurf[0]*fUpwind[4]+0.3535533905932737*fUpwind[0]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[1]*fUpwind[2]+0.3535533905932737*fUpwind[1]*alphaDrSurf[2]; 
  drag_incr[5] = 0.3535533905932737*alphaDrSurf[2]*fUpwind[7]+0.3535533905932737*alphaDrSurf[4]*fUpwind[6]+0.3535533905932737*alphaDrSurf[0]*fUpwind[5]+0.3535533905932737*alphaDrSurf[1]*fUpwind[3]; 
  drag_incr[6] = 0.3535533905932737*alphaDrSurf[1]*fUpwind[7]+0.3535533905932737*alphaDrSurf[0]*fUpwind[6]+0.3535533905932737*alphaDrSurf[4]*fUpwind[5]+0.3535533905932737*alphaDrSurf[2]*fUpwind[3]; 
  drag_incr[7] = 0.3535533905932737*alphaDrSurf[0]*fUpwind[7]+0.3535533905932737*alphaDrSurf[1]*fUpwind[6]+0.3535533905932737*alphaDrSurf[2]*fUpwind[5]+0.3535533905932737*fUpwind[3]*alphaDrSurf[4]; 

  out[0] += 0.7071067811865475*drag_incr[0]*rdv2; 
  out[1] += 0.7071067811865475*drag_incr[1]*rdv2; 
  out[2] += 0.7071067811865475*drag_incr[2]*rdv2; 
  out[3] += 0.7071067811865475*drag_incr[3]*rdv2; 
  out[4] += 1.224744871391589*drag_incr[0]*rdv2; 
  out[5] += 0.7071067811865475*drag_incr[4]*rdv2; 
  out[6] += 0.7071067811865475*drag_incr[5]*rdv2; 
  out[7] += 0.7071067811865475*drag_incr[6]*rdv2; 
  out[8] += 1.224744871391589*drag_incr[1]*rdv2; 
  out[9] += 1.224744871391589*drag_incr[2]*rdv2; 
  out[10] += 1.224744871391589*drag_incr[3]*rdv2; 
  out[11] += 0.7071067811865475*drag_incr[7]*rdv2; 
  out[12] += 1.224744871391589*drag_incr[4]*rdv2; 
  out[13] += 1.224744871391589*drag_incr[5]*rdv2; 
  out[14] += 1.224744871391589*drag_incr[6]*rdv2; 
  out[15] += 1.224744871391589*drag_incr[7]*rdv2; 

  } else { 

  alphaDrSurf[0] = 0.7071067811865475*(nuSum[0]*(2.0*w[3]-1.0*dxv[3])-2.0*sumNuUy[0]); 
  alphaDrSurf[1] = 0.7071067811865475*(nuSum[1]*(2.0*w[3]-1.0*dxv[3])-2.0*sumNuUy[1]); 
  alphaDrSurf[2] = 0.7071067811865475*(nuSum[2]*(2.0*w[3]-1.0*dxv[3])-2.0*sumNuUy[2]); 
  alphaDrSurf[4] = 0.7071067811865475*(2.0*nuSum[3]*w[3]-2.0*sumNuUy[3]-1.0*dxv[3]*nuSum[3]); 

  if (alphaDrSurf[4]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = ser_4x_p1_surfx4_quad_0_r(fEdge); 
    fUpwindQuad[4] = ser_4x_p1_surfx4_quad_4_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_4x_p1_surfx4_quad_0_l(fSkin); 
    fUpwindQuad[4] = ser_4x_p1_surfx4_quad_4_l(fSkin); 
  } 
  if ((-alphaDrSurf[4])-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[1] = ser_4x_p1_surfx4_quad_1_r(fEdge); 
    fUpwindQuad[5] = ser_4x_p1_surfx4_quad_5_r(fEdge); 
  } else { 
    fUpwindQuad[1] = ser_4x_p1_surfx4_quad_1_l(fSkin); 
    fUpwindQuad[5] = ser_4x_p1_surfx4_quad_5_l(fSkin); 
  } 
  if ((-alphaDrSurf[4])+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[2] = ser_4x_p1_surfx4_quad_2_r(fEdge); 
    fUpwindQuad[6] = ser_4x_p1_surfx4_quad_6_r(fEdge); 
  } else { 
    fUpwindQuad[2] = ser_4x_p1_surfx4_quad_2_l(fSkin); 
    fUpwindQuad[6] = ser_4x_p1_surfx4_quad_6_l(fSkin); 
  } 
  if (alphaDrSurf[4]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = ser_4x_p1_surfx4_quad_3_r(fEdge); 
    fUpwindQuad[7] = ser_4x_p1_surfx4_quad_7_r(fEdge); 
  } else { 
    fUpwindQuad[3] = ser_4x_p1_surfx4_quad_3_l(fSkin); 
    fUpwindQuad[7] = ser_4x_p1_surfx4_quad_7_l(fSkin); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_4x_p1_upwind(fUpwindQuad, fUpwind); 

  drag_incr[0] = 0.3535533905932737*alphaDrSurf[4]*fUpwind[4]+0.3535533905932737*alphaDrSurf[2]*fUpwind[2]+0.3535533905932737*alphaDrSurf[1]*fUpwind[1]+0.3535533905932737*alphaDrSurf[0]*fUpwind[0]; 
  drag_incr[1] = 0.3535533905932737*alphaDrSurf[2]*fUpwind[4]+0.3535533905932737*fUpwind[2]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[0]*fUpwind[1]+0.3535533905932737*fUpwind[0]*alphaDrSurf[1]; 
  drag_incr[2] = 0.3535533905932737*alphaDrSurf[1]*fUpwind[4]+0.3535533905932737*fUpwind[1]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[0]*fUpwind[2]+0.3535533905932737*fUpwind[0]*alphaDrSurf[2]; 
  drag_incr[3] = 0.3535533905932737*alphaDrSurf[4]*fUpwind[7]+0.3535533905932737*alphaDrSurf[2]*fUpwind[6]+0.3535533905932737*alphaDrSurf[1]*fUpwind[5]+0.3535533905932737*alphaDrSurf[0]*fUpwind[3]; 
  drag_incr[4] = 0.3535533905932737*alphaDrSurf[0]*fUpwind[4]+0.3535533905932737*fUpwind[0]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[1]*fUpwind[2]+0.3535533905932737*fUpwind[1]*alphaDrSurf[2]; 
  drag_incr[5] = 0.3535533905932737*alphaDrSurf[2]*fUpwind[7]+0.3535533905932737*alphaDrSurf[4]*fUpwind[6]+0.3535533905932737*alphaDrSurf[0]*fUpwind[5]+0.3535533905932737*alphaDrSurf[1]*fUpwind[3]; 
  drag_incr[6] = 0.3535533905932737*alphaDrSurf[1]*fUpwind[7]+0.3535533905932737*alphaDrSurf[0]*fUpwind[6]+0.3535533905932737*alphaDrSurf[4]*fUpwind[5]+0.3535533905932737*alphaDrSurf[2]*fUpwind[3]; 
  drag_incr[7] = 0.3535533905932737*alphaDrSurf[0]*fUpwind[7]+0.3535533905932737*alphaDrSurf[1]*fUpwind[6]+0.3535533905932737*alphaDrSurf[2]*fUpwind[5]+0.3535533905932737*fUpwind[3]*alphaDrSurf[4]; 

  out[0] += -0.7071067811865475*drag_incr[0]*rdv2; 
  out[1] += -0.7071067811865475*drag_incr[1]*rdv2; 
  out[2] += -0.7071067811865475*drag_incr[2]*rdv2; 
  out[3] += -0.7071067811865475*drag_incr[3]*rdv2; 
  out[4] += 1.224744871391589*drag_incr[0]*rdv2; 
  out[5] += -0.7071067811865475*drag_incr[4]*rdv2; 
  out[6] += -0.7071067811865475*drag_incr[5]*rdv2; 
  out[7] += -0.7071067811865475*drag_incr[6]*rdv2; 
  out[8] += 1.224744871391589*drag_incr[1]*rdv2; 
  out[9] += 1.224744871391589*drag_incr[2]*rdv2; 
  out[10] += 1.224744871391589*drag_incr[3]*rdv2; 
  out[11] += -0.7071067811865475*drag_incr[7]*rdv2; 
  out[12] += 1.224744871391589*drag_incr[4]*rdv2; 
  out[13] += 1.224744871391589*drag_incr[5]*rdv2; 
  out[14] += 1.224744871391589*drag_incr[6]*rdv2; 
  out[15] += 1.224744871391589*drag_incr[7]*rdv2; 

  } 
} 
