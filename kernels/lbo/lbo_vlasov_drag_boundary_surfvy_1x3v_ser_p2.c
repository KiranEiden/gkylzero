#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_ser_4x_p2_surfx3_quad.h> 
#include <gkyl_basis_ser_4x_p2_upwind.h> 
GKYL_CU_DH void lbo_vlasov_drag_boundary_surfvy_1x3v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[4]:         Cell-center coordinates. 
  // dxv[4]:       Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[9]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]: sum of thermal speeds squared time their respective collisionalities. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[2]; 

  const double *sumNuUy = &nuUSum[3]; 

  double alphaDrSurf[20] = {0.0}; 
  double fUpwindQuad[27] = {0.0};
  double fUpwind[20] = {0.0};;
  double drag_incr[20] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = nuSum[0]*(2.0*w[2]+dxv[2])-2.0*sumNuUy[0]; 
  alphaDrSurf[1] = nuSum[1]*(2.0*w[2]+dxv[2])-2.0*sumNuUy[1]; 
  alphaDrSurf[7] = 2.0*nuSum[2]*w[2]-2.0*sumNuUy[2]+dxv[2]*nuSum[2]; 

  if (0.3162277660168379*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = ser_4x_p2_surfx3_quad_0_r(fSkin); 
    fUpwindQuad[3] = ser_4x_p2_surfx3_quad_3_r(fSkin); 
    fUpwindQuad[6] = ser_4x_p2_surfx3_quad_6_r(fSkin); 
    fUpwindQuad[9] = ser_4x_p2_surfx3_quad_9_r(fSkin); 
    fUpwindQuad[12] = ser_4x_p2_surfx3_quad_12_r(fSkin); 
    fUpwindQuad[15] = ser_4x_p2_surfx3_quad_15_r(fSkin); 
    fUpwindQuad[18] = ser_4x_p2_surfx3_quad_18_r(fSkin); 
    fUpwindQuad[21] = ser_4x_p2_surfx3_quad_21_r(fSkin); 
    fUpwindQuad[24] = ser_4x_p2_surfx3_quad_24_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_4x_p2_surfx3_quad_0_l(fEdge); 
    fUpwindQuad[3] = ser_4x_p2_surfx3_quad_3_l(fEdge); 
    fUpwindQuad[6] = ser_4x_p2_surfx3_quad_6_l(fEdge); 
    fUpwindQuad[9] = ser_4x_p2_surfx3_quad_9_l(fEdge); 
    fUpwindQuad[12] = ser_4x_p2_surfx3_quad_12_l(fEdge); 
    fUpwindQuad[15] = ser_4x_p2_surfx3_quad_15_l(fEdge); 
    fUpwindQuad[18] = ser_4x_p2_surfx3_quad_18_l(fEdge); 
    fUpwindQuad[21] = ser_4x_p2_surfx3_quad_21_l(fEdge); 
    fUpwindQuad[24] = ser_4x_p2_surfx3_quad_24_l(fEdge); 
  } 
  if (0.3535533905932737*alphaDrSurf[0]-0.3952847075210473*alphaDrSurf[7] < 0) { 
    fUpwindQuad[1] = ser_4x_p2_surfx3_quad_1_r(fSkin); 
    fUpwindQuad[4] = ser_4x_p2_surfx3_quad_4_r(fSkin); 
    fUpwindQuad[7] = ser_4x_p2_surfx3_quad_7_r(fSkin); 
    fUpwindQuad[10] = ser_4x_p2_surfx3_quad_10_r(fSkin); 
    fUpwindQuad[13] = ser_4x_p2_surfx3_quad_13_r(fSkin); 
    fUpwindQuad[16] = ser_4x_p2_surfx3_quad_16_r(fSkin); 
    fUpwindQuad[19] = ser_4x_p2_surfx3_quad_19_r(fSkin); 
    fUpwindQuad[22] = ser_4x_p2_surfx3_quad_22_r(fSkin); 
    fUpwindQuad[25] = ser_4x_p2_surfx3_quad_25_r(fSkin); 
  } else { 
    fUpwindQuad[1] = ser_4x_p2_surfx3_quad_1_l(fEdge); 
    fUpwindQuad[4] = ser_4x_p2_surfx3_quad_4_l(fEdge); 
    fUpwindQuad[7] = ser_4x_p2_surfx3_quad_7_l(fEdge); 
    fUpwindQuad[10] = ser_4x_p2_surfx3_quad_10_l(fEdge); 
    fUpwindQuad[13] = ser_4x_p2_surfx3_quad_13_l(fEdge); 
    fUpwindQuad[16] = ser_4x_p2_surfx3_quad_16_l(fEdge); 
    fUpwindQuad[19] = ser_4x_p2_surfx3_quad_19_l(fEdge); 
    fUpwindQuad[22] = ser_4x_p2_surfx3_quad_22_l(fEdge); 
    fUpwindQuad[25] = ser_4x_p2_surfx3_quad_25_l(fEdge); 
  } 
  if (0.3162277660168379*alphaDrSurf[7]+0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[2] = ser_4x_p2_surfx3_quad_2_r(fSkin); 
    fUpwindQuad[5] = ser_4x_p2_surfx3_quad_5_r(fSkin); 
    fUpwindQuad[8] = ser_4x_p2_surfx3_quad_8_r(fSkin); 
    fUpwindQuad[11] = ser_4x_p2_surfx3_quad_11_r(fSkin); 
    fUpwindQuad[14] = ser_4x_p2_surfx3_quad_14_r(fSkin); 
    fUpwindQuad[17] = ser_4x_p2_surfx3_quad_17_r(fSkin); 
    fUpwindQuad[20] = ser_4x_p2_surfx3_quad_20_r(fSkin); 
    fUpwindQuad[23] = ser_4x_p2_surfx3_quad_23_r(fSkin); 
    fUpwindQuad[26] = ser_4x_p2_surfx3_quad_26_r(fSkin); 
  } else { 
    fUpwindQuad[2] = ser_4x_p2_surfx3_quad_2_l(fEdge); 
    fUpwindQuad[5] = ser_4x_p2_surfx3_quad_5_l(fEdge); 
    fUpwindQuad[8] = ser_4x_p2_surfx3_quad_8_l(fEdge); 
    fUpwindQuad[11] = ser_4x_p2_surfx3_quad_11_l(fEdge); 
    fUpwindQuad[14] = ser_4x_p2_surfx3_quad_14_l(fEdge); 
    fUpwindQuad[17] = ser_4x_p2_surfx3_quad_17_l(fEdge); 
    fUpwindQuad[20] = ser_4x_p2_surfx3_quad_20_l(fEdge); 
    fUpwindQuad[23] = ser_4x_p2_surfx3_quad_23_l(fEdge); 
    fUpwindQuad[26] = ser_4x_p2_surfx3_quad_26_l(fEdge); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_4x_p2_upwind(fUpwindQuad, fUpwind); 

  drag_incr[0] = 0.3535533905932737*alphaDrSurf[7]*fUpwind[7]+0.3535533905932737*alphaDrSurf[1]*fUpwind[1]+0.3535533905932737*alphaDrSurf[0]*fUpwind[0]; 
  drag_incr[1] = 0.3162277660168379*alphaDrSurf[1]*fUpwind[7]+0.3162277660168379*fUpwind[1]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[0]*fUpwind[1]+0.3535533905932737*fUpwind[0]*alphaDrSurf[1]; 
  drag_incr[2] = 0.3535533905932737*alphaDrSurf[7]*fUpwind[11]+0.3535533905932737*alphaDrSurf[1]*fUpwind[4]+0.3535533905932737*alphaDrSurf[0]*fUpwind[2]; 
  drag_incr[3] = 0.3535533905932737*alphaDrSurf[7]*fUpwind[13]+0.3535533905932737*alphaDrSurf[1]*fUpwind[5]+0.3535533905932737*alphaDrSurf[0]*fUpwind[3]; 
  drag_incr[4] = 0.3162277660168379*alphaDrSurf[1]*fUpwind[11]+0.3162277660168379*fUpwind[4]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[0]*fUpwind[4]+0.3535533905932737*alphaDrSurf[1]*fUpwind[2]; 
  drag_incr[5] = 0.3162277660168379*alphaDrSurf[1]*fUpwind[13]+0.3162277660168379*fUpwind[5]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[0]*fUpwind[5]+0.3535533905932737*alphaDrSurf[1]*fUpwind[3]; 
  drag_incr[6] = 0.3535533905932737*alphaDrSurf[7]*fUpwind[17]+0.3535533905932737*alphaDrSurf[1]*fUpwind[10]+0.3535533905932737*alphaDrSurf[0]*fUpwind[6]; 
  drag_incr[7] = 0.2258769757263128*alphaDrSurf[7]*fUpwind[7]+0.3535533905932737*alphaDrSurf[0]*fUpwind[7]+0.3535533905932737*fUpwind[0]*alphaDrSurf[7]+0.3162277660168379*alphaDrSurf[1]*fUpwind[1]; 
  drag_incr[8] = 0.3535533905932737*alphaDrSurf[1]*fUpwind[12]+0.3535533905932737*alphaDrSurf[0]*fUpwind[8]; 
  drag_incr[9] = 0.3535533905932737*alphaDrSurf[1]*fUpwind[15]+0.3535533905932737*alphaDrSurf[0]*fUpwind[9]; 
  drag_incr[10] = 0.3162277660168379*alphaDrSurf[1]*fUpwind[17]+0.3162277660168379*alphaDrSurf[7]*fUpwind[10]+0.3535533905932737*alphaDrSurf[0]*fUpwind[10]+0.3535533905932737*alphaDrSurf[1]*fUpwind[6]; 
  drag_incr[11] = 0.2258769757263128*alphaDrSurf[7]*fUpwind[11]+0.3535533905932737*alphaDrSurf[0]*fUpwind[11]+0.3535533905932737*fUpwind[2]*alphaDrSurf[7]+0.3162277660168379*alphaDrSurf[1]*fUpwind[4]; 
  drag_incr[12] = 0.3162277660168379*alphaDrSurf[7]*fUpwind[12]+0.3535533905932737*alphaDrSurf[0]*fUpwind[12]+0.3535533905932737*alphaDrSurf[1]*fUpwind[8]; 
  drag_incr[13] = 0.2258769757263128*alphaDrSurf[7]*fUpwind[13]+0.3535533905932737*alphaDrSurf[0]*fUpwind[13]+0.3535533905932737*fUpwind[3]*alphaDrSurf[7]+0.3162277660168379*alphaDrSurf[1]*fUpwind[5]; 
  drag_incr[14] = 0.3535533905932737*alphaDrSurf[1]*fUpwind[18]+0.3535533905932737*alphaDrSurf[0]*fUpwind[14]; 
  drag_incr[15] = 0.3162277660168379*alphaDrSurf[7]*fUpwind[15]+0.3535533905932737*alphaDrSurf[0]*fUpwind[15]+0.3535533905932737*alphaDrSurf[1]*fUpwind[9]; 
  drag_incr[16] = 0.3535533905932737*alphaDrSurf[1]*fUpwind[19]+0.3535533905932737*alphaDrSurf[0]*fUpwind[16]; 
  drag_incr[17] = 0.2258769757263128*alphaDrSurf[7]*fUpwind[17]+0.3535533905932737*alphaDrSurf[0]*fUpwind[17]+0.3162277660168379*alphaDrSurf[1]*fUpwind[10]+0.3535533905932737*fUpwind[6]*alphaDrSurf[7]; 
  drag_incr[18] = 0.3162277660168379*alphaDrSurf[7]*fUpwind[18]+0.3535533905932737*alphaDrSurf[0]*fUpwind[18]+0.3535533905932737*alphaDrSurf[1]*fUpwind[14]; 
  drag_incr[19] = 0.3162277660168379*alphaDrSurf[7]*fUpwind[19]+0.3535533905932737*alphaDrSurf[0]*fUpwind[19]+0.3535533905932737*alphaDrSurf[1]*fUpwind[16]; 

  out[0] += 0.7071067811865475*drag_incr[0]*rdv2; 
  out[1] += 0.7071067811865475*drag_incr[1]*rdv2; 
  out[2] += 0.7071067811865475*drag_incr[2]*rdv2; 
  out[3] += 1.224744871391589*drag_incr[0]*rdv2; 
  out[4] += 0.7071067811865475*drag_incr[3]*rdv2; 
  out[5] += 0.7071067811865475*drag_incr[4]*rdv2; 
  out[6] += 1.224744871391589*drag_incr[1]*rdv2; 
  out[7] += 1.224744871391589*drag_incr[2]*rdv2; 
  out[8] += 0.7071067811865475*drag_incr[5]*rdv2; 
  out[9] += 0.7071067811865475*drag_incr[6]*rdv2; 
  out[10] += 1.224744871391589*drag_incr[3]*rdv2; 
  out[11] += 0.7071067811865475*drag_incr[7]*rdv2; 
  out[12] += 0.7071067811865475*drag_incr[8]*rdv2; 
  out[13] += 1.58113883008419*drag_incr[0]*rdv2; 
  out[14] += 0.7071067811865475*drag_incr[9]*rdv2; 
  out[15] += 1.224744871391589*drag_incr[4]*rdv2; 
  out[16] += 0.7071067811865475*drag_incr[10]*rdv2; 
  out[17] += 1.224744871391589*drag_incr[5]*rdv2; 
  out[18] += 1.224744871391589*drag_incr[6]*rdv2; 
  out[19] += 0.7071067811865475*drag_incr[11]*rdv2; 
  out[20] += 0.7071067811865475*drag_incr[12]*rdv2; 
  out[21] += 1.224744871391589*drag_incr[7]*rdv2; 
  out[22] += 1.224744871391589*drag_incr[8]*rdv2; 
  out[23] += 1.58113883008419*drag_incr[1]*rdv2; 
  out[24] += 1.58113883008419*drag_incr[2]*rdv2; 
  out[25] += 0.7071067811865475*drag_incr[13]*rdv2; 
  out[26] += 0.7071067811865475*drag_incr[14]*rdv2; 
  out[27] += 1.58113883008419*drag_incr[3]*rdv2; 
  out[28] += 0.7071067811865475*drag_incr[15]*rdv2; 
  out[29] += 0.7071067811865475*drag_incr[16]*rdv2; 
  out[30] += 1.224744871391589*drag_incr[9]*rdv2; 
  out[31] += 1.224744871391589*drag_incr[10]*rdv2; 
  out[32] += 1.224744871391589*drag_incr[11]*rdv2; 
  out[33] += 1.224744871391589*drag_incr[12]*rdv2; 
  out[34] += 1.58113883008419*drag_incr[4]*rdv2; 
  out[35] += 0.7071067811865475*drag_incr[17]*rdv2; 
  out[36] += 0.7071067811865475*drag_incr[18]*rdv2; 
  out[37] += 1.224744871391589*drag_incr[13]*rdv2; 
  out[38] += 1.224744871391589*drag_incr[14]*rdv2; 
  out[39] += 1.58113883008419*drag_incr[5]*rdv2; 
  out[40] += 1.58113883008419*drag_incr[6]*rdv2; 
  out[41] += 0.7071067811865475*drag_incr[19]*rdv2; 
  out[42] += 1.224744871391589*drag_incr[15]*rdv2; 
  out[43] += 1.224744871391589*drag_incr[16]*rdv2; 
  out[44] += 1.224744871391589*drag_incr[17]*rdv2; 
  out[45] += 1.224744871391589*drag_incr[18]*rdv2; 
  out[46] += 1.58113883008419*drag_incr[10]*rdv2; 
  out[47] += 1.224744871391589*drag_incr[19]*rdv2; 

  } else { 

  alphaDrSurf[0] = nuSum[0]*(2.0*w[2]-1.0*dxv[2])-2.0*sumNuUy[0]; 
  alphaDrSurf[1] = nuSum[1]*(2.0*w[2]-1.0*dxv[2])-2.0*sumNuUy[1]; 
  alphaDrSurf[7] = 2.0*nuSum[2]*w[2]-2.0*sumNuUy[2]-1.0*dxv[2]*nuSum[2]; 

  if (0.3162277660168379*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = ser_4x_p2_surfx3_quad_0_r(fEdge); 
    fUpwindQuad[3] = ser_4x_p2_surfx3_quad_3_r(fEdge); 
    fUpwindQuad[6] = ser_4x_p2_surfx3_quad_6_r(fEdge); 
    fUpwindQuad[9] = ser_4x_p2_surfx3_quad_9_r(fEdge); 
    fUpwindQuad[12] = ser_4x_p2_surfx3_quad_12_r(fEdge); 
    fUpwindQuad[15] = ser_4x_p2_surfx3_quad_15_r(fEdge); 
    fUpwindQuad[18] = ser_4x_p2_surfx3_quad_18_r(fEdge); 
    fUpwindQuad[21] = ser_4x_p2_surfx3_quad_21_r(fEdge); 
    fUpwindQuad[24] = ser_4x_p2_surfx3_quad_24_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_4x_p2_surfx3_quad_0_l(fSkin); 
    fUpwindQuad[3] = ser_4x_p2_surfx3_quad_3_l(fSkin); 
    fUpwindQuad[6] = ser_4x_p2_surfx3_quad_6_l(fSkin); 
    fUpwindQuad[9] = ser_4x_p2_surfx3_quad_9_l(fSkin); 
    fUpwindQuad[12] = ser_4x_p2_surfx3_quad_12_l(fSkin); 
    fUpwindQuad[15] = ser_4x_p2_surfx3_quad_15_l(fSkin); 
    fUpwindQuad[18] = ser_4x_p2_surfx3_quad_18_l(fSkin); 
    fUpwindQuad[21] = ser_4x_p2_surfx3_quad_21_l(fSkin); 
    fUpwindQuad[24] = ser_4x_p2_surfx3_quad_24_l(fSkin); 
  } 
  if (0.3535533905932737*alphaDrSurf[0]-0.3952847075210473*alphaDrSurf[7] < 0) { 
    fUpwindQuad[1] = ser_4x_p2_surfx3_quad_1_r(fEdge); 
    fUpwindQuad[4] = ser_4x_p2_surfx3_quad_4_r(fEdge); 
    fUpwindQuad[7] = ser_4x_p2_surfx3_quad_7_r(fEdge); 
    fUpwindQuad[10] = ser_4x_p2_surfx3_quad_10_r(fEdge); 
    fUpwindQuad[13] = ser_4x_p2_surfx3_quad_13_r(fEdge); 
    fUpwindQuad[16] = ser_4x_p2_surfx3_quad_16_r(fEdge); 
    fUpwindQuad[19] = ser_4x_p2_surfx3_quad_19_r(fEdge); 
    fUpwindQuad[22] = ser_4x_p2_surfx3_quad_22_r(fEdge); 
    fUpwindQuad[25] = ser_4x_p2_surfx3_quad_25_r(fEdge); 
  } else { 
    fUpwindQuad[1] = ser_4x_p2_surfx3_quad_1_l(fSkin); 
    fUpwindQuad[4] = ser_4x_p2_surfx3_quad_4_l(fSkin); 
    fUpwindQuad[7] = ser_4x_p2_surfx3_quad_7_l(fSkin); 
    fUpwindQuad[10] = ser_4x_p2_surfx3_quad_10_l(fSkin); 
    fUpwindQuad[13] = ser_4x_p2_surfx3_quad_13_l(fSkin); 
    fUpwindQuad[16] = ser_4x_p2_surfx3_quad_16_l(fSkin); 
    fUpwindQuad[19] = ser_4x_p2_surfx3_quad_19_l(fSkin); 
    fUpwindQuad[22] = ser_4x_p2_surfx3_quad_22_l(fSkin); 
    fUpwindQuad[25] = ser_4x_p2_surfx3_quad_25_l(fSkin); 
  } 
  if (0.3162277660168379*alphaDrSurf[7]+0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[2] = ser_4x_p2_surfx3_quad_2_r(fEdge); 
    fUpwindQuad[5] = ser_4x_p2_surfx3_quad_5_r(fEdge); 
    fUpwindQuad[8] = ser_4x_p2_surfx3_quad_8_r(fEdge); 
    fUpwindQuad[11] = ser_4x_p2_surfx3_quad_11_r(fEdge); 
    fUpwindQuad[14] = ser_4x_p2_surfx3_quad_14_r(fEdge); 
    fUpwindQuad[17] = ser_4x_p2_surfx3_quad_17_r(fEdge); 
    fUpwindQuad[20] = ser_4x_p2_surfx3_quad_20_r(fEdge); 
    fUpwindQuad[23] = ser_4x_p2_surfx3_quad_23_r(fEdge); 
    fUpwindQuad[26] = ser_4x_p2_surfx3_quad_26_r(fEdge); 
  } else { 
    fUpwindQuad[2] = ser_4x_p2_surfx3_quad_2_l(fSkin); 
    fUpwindQuad[5] = ser_4x_p2_surfx3_quad_5_l(fSkin); 
    fUpwindQuad[8] = ser_4x_p2_surfx3_quad_8_l(fSkin); 
    fUpwindQuad[11] = ser_4x_p2_surfx3_quad_11_l(fSkin); 
    fUpwindQuad[14] = ser_4x_p2_surfx3_quad_14_l(fSkin); 
    fUpwindQuad[17] = ser_4x_p2_surfx3_quad_17_l(fSkin); 
    fUpwindQuad[20] = ser_4x_p2_surfx3_quad_20_l(fSkin); 
    fUpwindQuad[23] = ser_4x_p2_surfx3_quad_23_l(fSkin); 
    fUpwindQuad[26] = ser_4x_p2_surfx3_quad_26_l(fSkin); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_4x_p2_upwind(fUpwindQuad, fUpwind); 

  drag_incr[0] = 0.3535533905932737*alphaDrSurf[7]*fUpwind[7]+0.3535533905932737*alphaDrSurf[1]*fUpwind[1]+0.3535533905932737*alphaDrSurf[0]*fUpwind[0]; 
  drag_incr[1] = 0.3162277660168379*alphaDrSurf[1]*fUpwind[7]+0.3162277660168379*fUpwind[1]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[0]*fUpwind[1]+0.3535533905932737*fUpwind[0]*alphaDrSurf[1]; 
  drag_incr[2] = 0.3535533905932737*alphaDrSurf[7]*fUpwind[11]+0.3535533905932737*alphaDrSurf[1]*fUpwind[4]+0.3535533905932737*alphaDrSurf[0]*fUpwind[2]; 
  drag_incr[3] = 0.3535533905932737*alphaDrSurf[7]*fUpwind[13]+0.3535533905932737*alphaDrSurf[1]*fUpwind[5]+0.3535533905932737*alphaDrSurf[0]*fUpwind[3]; 
  drag_incr[4] = 0.3162277660168379*alphaDrSurf[1]*fUpwind[11]+0.3162277660168379*fUpwind[4]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[0]*fUpwind[4]+0.3535533905932737*alphaDrSurf[1]*fUpwind[2]; 
  drag_incr[5] = 0.3162277660168379*alphaDrSurf[1]*fUpwind[13]+0.3162277660168379*fUpwind[5]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[0]*fUpwind[5]+0.3535533905932737*alphaDrSurf[1]*fUpwind[3]; 
  drag_incr[6] = 0.3535533905932737*alphaDrSurf[7]*fUpwind[17]+0.3535533905932737*alphaDrSurf[1]*fUpwind[10]+0.3535533905932737*alphaDrSurf[0]*fUpwind[6]; 
  drag_incr[7] = 0.2258769757263128*alphaDrSurf[7]*fUpwind[7]+0.3535533905932737*alphaDrSurf[0]*fUpwind[7]+0.3535533905932737*fUpwind[0]*alphaDrSurf[7]+0.3162277660168379*alphaDrSurf[1]*fUpwind[1]; 
  drag_incr[8] = 0.3535533905932737*alphaDrSurf[1]*fUpwind[12]+0.3535533905932737*alphaDrSurf[0]*fUpwind[8]; 
  drag_incr[9] = 0.3535533905932737*alphaDrSurf[1]*fUpwind[15]+0.3535533905932737*alphaDrSurf[0]*fUpwind[9]; 
  drag_incr[10] = 0.3162277660168379*alphaDrSurf[1]*fUpwind[17]+0.3162277660168379*alphaDrSurf[7]*fUpwind[10]+0.3535533905932737*alphaDrSurf[0]*fUpwind[10]+0.3535533905932737*alphaDrSurf[1]*fUpwind[6]; 
  drag_incr[11] = 0.2258769757263128*alphaDrSurf[7]*fUpwind[11]+0.3535533905932737*alphaDrSurf[0]*fUpwind[11]+0.3535533905932737*fUpwind[2]*alphaDrSurf[7]+0.3162277660168379*alphaDrSurf[1]*fUpwind[4]; 
  drag_incr[12] = 0.3162277660168379*alphaDrSurf[7]*fUpwind[12]+0.3535533905932737*alphaDrSurf[0]*fUpwind[12]+0.3535533905932737*alphaDrSurf[1]*fUpwind[8]; 
  drag_incr[13] = 0.2258769757263128*alphaDrSurf[7]*fUpwind[13]+0.3535533905932737*alphaDrSurf[0]*fUpwind[13]+0.3535533905932737*fUpwind[3]*alphaDrSurf[7]+0.3162277660168379*alphaDrSurf[1]*fUpwind[5]; 
  drag_incr[14] = 0.3535533905932737*alphaDrSurf[1]*fUpwind[18]+0.3535533905932737*alphaDrSurf[0]*fUpwind[14]; 
  drag_incr[15] = 0.3162277660168379*alphaDrSurf[7]*fUpwind[15]+0.3535533905932737*alphaDrSurf[0]*fUpwind[15]+0.3535533905932737*alphaDrSurf[1]*fUpwind[9]; 
  drag_incr[16] = 0.3535533905932737*alphaDrSurf[1]*fUpwind[19]+0.3535533905932737*alphaDrSurf[0]*fUpwind[16]; 
  drag_incr[17] = 0.2258769757263128*alphaDrSurf[7]*fUpwind[17]+0.3535533905932737*alphaDrSurf[0]*fUpwind[17]+0.3162277660168379*alphaDrSurf[1]*fUpwind[10]+0.3535533905932737*fUpwind[6]*alphaDrSurf[7]; 
  drag_incr[18] = 0.3162277660168379*alphaDrSurf[7]*fUpwind[18]+0.3535533905932737*alphaDrSurf[0]*fUpwind[18]+0.3535533905932737*alphaDrSurf[1]*fUpwind[14]; 
  drag_incr[19] = 0.3162277660168379*alphaDrSurf[7]*fUpwind[19]+0.3535533905932737*alphaDrSurf[0]*fUpwind[19]+0.3535533905932737*alphaDrSurf[1]*fUpwind[16]; 

  out[0] += -0.7071067811865475*drag_incr[0]*rdv2; 
  out[1] += -0.7071067811865475*drag_incr[1]*rdv2; 
  out[2] += -0.7071067811865475*drag_incr[2]*rdv2; 
  out[3] += 1.224744871391589*drag_incr[0]*rdv2; 
  out[4] += -0.7071067811865475*drag_incr[3]*rdv2; 
  out[5] += -0.7071067811865475*drag_incr[4]*rdv2; 
  out[6] += 1.224744871391589*drag_incr[1]*rdv2; 
  out[7] += 1.224744871391589*drag_incr[2]*rdv2; 
  out[8] += -0.7071067811865475*drag_incr[5]*rdv2; 
  out[9] += -0.7071067811865475*drag_incr[6]*rdv2; 
  out[10] += 1.224744871391589*drag_incr[3]*rdv2; 
  out[11] += -0.7071067811865475*drag_incr[7]*rdv2; 
  out[12] += -0.7071067811865475*drag_incr[8]*rdv2; 
  out[13] += -1.58113883008419*drag_incr[0]*rdv2; 
  out[14] += -0.7071067811865475*drag_incr[9]*rdv2; 
  out[15] += 1.224744871391589*drag_incr[4]*rdv2; 
  out[16] += -0.7071067811865475*drag_incr[10]*rdv2; 
  out[17] += 1.224744871391589*drag_incr[5]*rdv2; 
  out[18] += 1.224744871391589*drag_incr[6]*rdv2; 
  out[19] += -0.7071067811865475*drag_incr[11]*rdv2; 
  out[20] += -0.7071067811865475*drag_incr[12]*rdv2; 
  out[21] += 1.224744871391589*drag_incr[7]*rdv2; 
  out[22] += 1.224744871391589*drag_incr[8]*rdv2; 
  out[23] += -1.58113883008419*drag_incr[1]*rdv2; 
  out[24] += -1.58113883008419*drag_incr[2]*rdv2; 
  out[25] += -0.7071067811865475*drag_incr[13]*rdv2; 
  out[26] += -0.7071067811865475*drag_incr[14]*rdv2; 
  out[27] += -1.58113883008419*drag_incr[3]*rdv2; 
  out[28] += -0.7071067811865475*drag_incr[15]*rdv2; 
  out[29] += -0.7071067811865475*drag_incr[16]*rdv2; 
  out[30] += 1.224744871391589*drag_incr[9]*rdv2; 
  out[31] += 1.224744871391589*drag_incr[10]*rdv2; 
  out[32] += 1.224744871391589*drag_incr[11]*rdv2; 
  out[33] += 1.224744871391589*drag_incr[12]*rdv2; 
  out[34] += -1.58113883008419*drag_incr[4]*rdv2; 
  out[35] += -0.7071067811865475*drag_incr[17]*rdv2; 
  out[36] += -0.7071067811865475*drag_incr[18]*rdv2; 
  out[37] += 1.224744871391589*drag_incr[13]*rdv2; 
  out[38] += 1.224744871391589*drag_incr[14]*rdv2; 
  out[39] += -1.58113883008419*drag_incr[5]*rdv2; 
  out[40] += -1.58113883008419*drag_incr[6]*rdv2; 
  out[41] += -0.7071067811865475*drag_incr[19]*rdv2; 
  out[42] += 1.224744871391589*drag_incr[15]*rdv2; 
  out[43] += 1.224744871391589*drag_incr[16]*rdv2; 
  out[44] += 1.224744871391589*drag_incr[17]*rdv2; 
  out[45] += 1.224744871391589*drag_incr[18]*rdv2; 
  out[46] += -1.58113883008419*drag_incr[10]*rdv2; 
  out[47] += 1.224744871391589*drag_incr[19]*rdv2; 

  } 
} 
