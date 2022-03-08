#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_4x_p2_surfx4_quad.h> 
#include <gkyl_basis_ser_4x_p2_upwind.h> 
GKYL_CU_DH void vlasov_poisson_extem_boundary_surfvz_1x3v_ser_p2(const double *w, const double *dxv, const double *fac_phi, const double *vecA, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // fac_phi:     potential (scaled by appropriate factors).
  // vecA:        vector potential (scaled by appropriate factors). Unused in pure Vlasov-Poisson. 
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv12 = 2/dxv[3]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double dv3 = dxv[3], wv3 = w[3]; 
  const double *phi = &fac_phi[0]; 
  const double dx10 = 2/dxv[0]; 
  const double *A0 = &vecA[0]; 
  const double *A1 = &vecA[3]; 
  const double *A2 = &vecA[6]; 
  double alpha[20] = {0.0}; 

  alpha[0] = -3.464101615137754*A2[1]*dx10*wv1; 
  alpha[1] = -7.745966692414834*A2[2]*dx10*wv1; 
  alpha[2] = -1.0*A2[1]*dv1*dx10; 
  alpha[4] = -2.23606797749979*A2[2]*dv1*dx10; 

  double fUpwindQuad[27] = {0.0};
  double fUpwind[20] = {0.0};
  double Ghat[20] = {0.0}; 

  if (edge == -1) { 

  if (0.6363961030678926*alpha[4]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[0] = ser_4x_p2_surfx4_quad_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_4x_p2_surfx4_quad_0_l(fEdge); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[2] > 0) { 
    fUpwindQuad[1] = ser_4x_p2_surfx4_quad_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = ser_4x_p2_surfx4_quad_1_l(fEdge); 
  } 
  if ((-0.6363961030678926*alpha[4])-0.4743416490252568*alpha[2]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[2] = ser_4x_p2_surfx4_quad_2_r(fSkin); 
  } else { 
    fUpwindQuad[2] = ser_4x_p2_surfx4_quad_2_l(fEdge); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[1] > 0) { 
    fUpwindQuad[3] = ser_4x_p2_surfx4_quad_3_r(fSkin); 
  } else { 
    fUpwindQuad[3] = ser_4x_p2_surfx4_quad_3_l(fEdge); 
  } 
  if (0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[4] = ser_4x_p2_surfx4_quad_4_r(fSkin); 
  } else { 
    fUpwindQuad[4] = ser_4x_p2_surfx4_quad_4_l(fEdge); 
  } 
  if (0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[5] = ser_4x_p2_surfx4_quad_5_r(fSkin); 
  } else { 
    fUpwindQuad[5] = ser_4x_p2_surfx4_quad_5_l(fEdge); 
  } 
  if ((-0.6363961030678926*alpha[4])+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[6] = ser_4x_p2_surfx4_quad_6_r(fSkin); 
  } else { 
    fUpwindQuad[6] = ser_4x_p2_surfx4_quad_6_l(fEdge); 
  } 
  if (0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[7] = ser_4x_p2_surfx4_quad_7_r(fSkin); 
  } else { 
    fUpwindQuad[7] = ser_4x_p2_surfx4_quad_7_l(fEdge); 
  } 
  if (0.6363961030678926*alpha[4]+0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[8] = ser_4x_p2_surfx4_quad_8_r(fSkin); 
  } else { 
    fUpwindQuad[8] = ser_4x_p2_surfx4_quad_8_l(fEdge); 
  } 
  if (0.6363961030678926*alpha[4]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[9] = ser_4x_p2_surfx4_quad_9_r(fSkin); 
  } else { 
    fUpwindQuad[9] = ser_4x_p2_surfx4_quad_9_l(fEdge); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[2] > 0) { 
    fUpwindQuad[10] = ser_4x_p2_surfx4_quad_10_r(fSkin); 
  } else { 
    fUpwindQuad[10] = ser_4x_p2_surfx4_quad_10_l(fEdge); 
  } 
  if ((-0.6363961030678926*alpha[4])-0.4743416490252568*alpha[2]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[11] = ser_4x_p2_surfx4_quad_11_r(fSkin); 
  } else { 
    fUpwindQuad[11] = ser_4x_p2_surfx4_quad_11_l(fEdge); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[1] > 0) { 
    fUpwindQuad[12] = ser_4x_p2_surfx4_quad_12_r(fSkin); 
  } else { 
    fUpwindQuad[12] = ser_4x_p2_surfx4_quad_12_l(fEdge); 
  } 
  if (0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[13] = ser_4x_p2_surfx4_quad_13_r(fSkin); 
  } else { 
    fUpwindQuad[13] = ser_4x_p2_surfx4_quad_13_l(fEdge); 
  } 
  if (0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[14] = ser_4x_p2_surfx4_quad_14_r(fSkin); 
  } else { 
    fUpwindQuad[14] = ser_4x_p2_surfx4_quad_14_l(fEdge); 
  } 
  if ((-0.6363961030678926*alpha[4])+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[15] = ser_4x_p2_surfx4_quad_15_r(fSkin); 
  } else { 
    fUpwindQuad[15] = ser_4x_p2_surfx4_quad_15_l(fEdge); 
  } 
  if (0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[16] = ser_4x_p2_surfx4_quad_16_r(fSkin); 
  } else { 
    fUpwindQuad[16] = ser_4x_p2_surfx4_quad_16_l(fEdge); 
  } 
  if (0.6363961030678926*alpha[4]+0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[17] = ser_4x_p2_surfx4_quad_17_r(fSkin); 
  } else { 
    fUpwindQuad[17] = ser_4x_p2_surfx4_quad_17_l(fEdge); 
  } 
  if (0.6363961030678926*alpha[4]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[18] = ser_4x_p2_surfx4_quad_18_r(fSkin); 
  } else { 
    fUpwindQuad[18] = ser_4x_p2_surfx4_quad_18_l(fEdge); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[2] > 0) { 
    fUpwindQuad[19] = ser_4x_p2_surfx4_quad_19_r(fSkin); 
  } else { 
    fUpwindQuad[19] = ser_4x_p2_surfx4_quad_19_l(fEdge); 
  } 
  if ((-0.6363961030678926*alpha[4])-0.4743416490252568*alpha[2]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[20] = ser_4x_p2_surfx4_quad_20_r(fSkin); 
  } else { 
    fUpwindQuad[20] = ser_4x_p2_surfx4_quad_20_l(fEdge); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[1] > 0) { 
    fUpwindQuad[21] = ser_4x_p2_surfx4_quad_21_r(fSkin); 
  } else { 
    fUpwindQuad[21] = ser_4x_p2_surfx4_quad_21_l(fEdge); 
  } 
  if (0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[22] = ser_4x_p2_surfx4_quad_22_r(fSkin); 
  } else { 
    fUpwindQuad[22] = ser_4x_p2_surfx4_quad_22_l(fEdge); 
  } 
  if (0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[23] = ser_4x_p2_surfx4_quad_23_r(fSkin); 
  } else { 
    fUpwindQuad[23] = ser_4x_p2_surfx4_quad_23_l(fEdge); 
  } 
  if ((-0.6363961030678926*alpha[4])+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[24] = ser_4x_p2_surfx4_quad_24_r(fSkin); 
  } else { 
    fUpwindQuad[24] = ser_4x_p2_surfx4_quad_24_l(fEdge); 
  } 
  if (0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[25] = ser_4x_p2_surfx4_quad_25_r(fSkin); 
  } else { 
    fUpwindQuad[25] = ser_4x_p2_surfx4_quad_25_l(fEdge); 
  } 
  if (0.6363961030678926*alpha[4]+0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[26] = ser_4x_p2_surfx4_quad_26_r(fSkin); 
  } else { 
    fUpwindQuad[26] = ser_4x_p2_surfx4_quad_26_l(fEdge); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_4x_p2_upwind(fUpwindQuad, fUpwind); 

  Ghat[0] += 0.3535533905932737*(alpha[4]*fUpwind[4]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] += 0.3162277660168379*alpha[4]*fUpwind[11]+0.3162277660168379*alpha[1]*fUpwind[7]+0.3535533905932737*(alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] += 0.3162277660168379*alpha[4]*fUpwind[12]+0.3162277660168379*alpha[2]*fUpwind[8]+0.3535533905932737*(alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] += 0.3535533905932737*(alpha[4]*fUpwind[10]+alpha[2]*fUpwind[6]+alpha[1]*fUpwind[5]+alpha[0]*fUpwind[3]); 
  Ghat[4] += 0.3162277660168379*(alpha[2]*fUpwind[12]+alpha[1]*fUpwind[11])+0.3162277660168379*alpha[4]*(fUpwind[8]+fUpwind[7])+0.3535533905932737*(alpha[0]*fUpwind[4]+fUpwind[0]*alpha[4]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[5] += 0.3162277660168379*alpha[4]*fUpwind[17]+0.3162277660168379*alpha[1]*fUpwind[13]+0.3535533905932737*(alpha[2]*fUpwind[10]+alpha[4]*fUpwind[6]+alpha[0]*fUpwind[5]+alpha[1]*fUpwind[3]); 
  Ghat[6] += 0.3162277660168379*alpha[4]*fUpwind[18]+0.3162277660168379*alpha[2]*fUpwind[14]+0.3535533905932737*(alpha[1]*fUpwind[10]+alpha[0]*fUpwind[6]+alpha[4]*fUpwind[5]+alpha[2]*fUpwind[3]); 
  Ghat[7] += 0.3535533905932737*(alpha[2]*fUpwind[11]+alpha[0]*fUpwind[7])+0.3162277660168379*(alpha[4]*fUpwind[4]+alpha[1]*fUpwind[1]); 
  Ghat[8] += 0.3535533905932737*(alpha[1]*fUpwind[12]+alpha[0]*fUpwind[8])+0.3162277660168379*(alpha[4]*fUpwind[4]+alpha[2]*fUpwind[2]); 
  Ghat[9] += 0.3535533905932737*(alpha[4]*fUpwind[19]+alpha[2]*fUpwind[16]+alpha[1]*fUpwind[15]+alpha[0]*fUpwind[9]); 
  Ghat[10] += 0.3162277660168379*(alpha[2]*fUpwind[18]+alpha[1]*fUpwind[17])+0.3162277660168379*alpha[4]*(fUpwind[14]+fUpwind[13])+0.3535533905932737*(alpha[0]*fUpwind[10]+alpha[1]*fUpwind[6]+alpha[2]*fUpwind[5]+fUpwind[3]*alpha[4]); 
  Ghat[11] += 0.2828427124746191*alpha[4]*fUpwind[12]+0.3535533905932737*(alpha[0]*fUpwind[11]+alpha[2]*fUpwind[7])+0.3162277660168379*(alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]); 
  Ghat[12] += 0.3535533905932737*alpha[0]*fUpwind[12]+0.2828427124746191*alpha[4]*fUpwind[11]+0.3535533905932737*alpha[1]*fUpwind[8]+0.3162277660168379*(alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]); 
  Ghat[13] += 0.3535533905932737*(alpha[2]*fUpwind[17]+alpha[0]*fUpwind[13])+0.3162277660168379*(alpha[4]*fUpwind[10]+alpha[1]*fUpwind[5]); 
  Ghat[14] += 0.3535533905932737*(alpha[1]*fUpwind[18]+alpha[0]*fUpwind[14])+0.3162277660168379*(alpha[4]*fUpwind[10]+alpha[2]*fUpwind[6]); 
  Ghat[15] += 0.3535533905932737*(alpha[2]*fUpwind[19]+alpha[4]*fUpwind[16]+alpha[0]*fUpwind[15]+alpha[1]*fUpwind[9]); 
  Ghat[16] += 0.3535533905932737*(alpha[1]*fUpwind[19]+alpha[0]*fUpwind[16]+alpha[4]*fUpwind[15]+alpha[2]*fUpwind[9]); 
  Ghat[17] += 0.2828427124746191*alpha[4]*fUpwind[18]+0.3535533905932737*(alpha[0]*fUpwind[17]+alpha[2]*fUpwind[13])+0.3162277660168379*(alpha[1]*fUpwind[10]+alpha[4]*fUpwind[5]); 
  Ghat[18] += 0.3535533905932737*alpha[0]*fUpwind[18]+0.2828427124746191*alpha[4]*fUpwind[17]+0.3535533905932737*alpha[1]*fUpwind[14]+0.3162277660168379*(alpha[2]*fUpwind[10]+alpha[4]*fUpwind[6]); 
  Ghat[19] += 0.3535533905932737*(alpha[0]*fUpwind[19]+alpha[1]*fUpwind[16]+alpha[2]*fUpwind[15]+alpha[4]*fUpwind[9]); 

  out[0] += -0.7071067811865475*Ghat[0]*dv12; 
  out[1] += -0.7071067811865475*Ghat[1]*dv12; 
  out[2] += -0.7071067811865475*Ghat[2]*dv12; 
  out[3] += -0.7071067811865475*Ghat[3]*dv12; 
  out[4] += -1.224744871391589*Ghat[0]*dv12; 
  out[5] += -0.7071067811865475*Ghat[4]*dv12; 
  out[6] += -0.7071067811865475*Ghat[5]*dv12; 
  out[7] += -0.7071067811865475*Ghat[6]*dv12; 
  out[8] += -1.224744871391589*Ghat[1]*dv12; 
  out[9] += -1.224744871391589*Ghat[2]*dv12; 
  out[10] += -1.224744871391589*Ghat[3]*dv12; 
  out[11] += -0.7071067811865475*Ghat[7]*dv12; 
  out[12] += -0.7071067811865475*Ghat[8]*dv12; 
  out[13] += -0.7071067811865475*Ghat[9]*dv12; 
  out[14] += -1.58113883008419*Ghat[0]*dv12; 
  out[15] += -0.7071067811865475*Ghat[10]*dv12; 
  out[16] += -1.224744871391589*Ghat[4]*dv12; 
  out[17] += -1.224744871391589*Ghat[5]*dv12; 
  out[18] += -1.224744871391589*Ghat[6]*dv12; 
  out[19] += -0.7071067811865475*Ghat[11]*dv12; 
  out[20] += -0.7071067811865475*Ghat[12]*dv12; 
  out[21] += -0.7071067811865475*Ghat[13]*dv12; 
  out[22] += -0.7071067811865475*Ghat[14]*dv12; 
  out[23] += -0.7071067811865475*Ghat[15]*dv12; 
  out[24] += -0.7071067811865475*Ghat[16]*dv12; 
  out[25] += -1.224744871391589*Ghat[7]*dv12; 
  out[26] += -1.224744871391589*Ghat[8]*dv12; 
  out[27] += -1.224744871391589*Ghat[9]*dv12; 
  out[28] += -1.58113883008419*Ghat[1]*dv12; 
  out[29] += -1.58113883008419*Ghat[2]*dv12; 
  out[30] += -1.58113883008419*Ghat[3]*dv12; 
  out[31] += -1.224744871391589*Ghat[10]*dv12; 
  out[32] += -0.7071067811865475*Ghat[17]*dv12; 
  out[33] += -0.7071067811865475*Ghat[18]*dv12; 
  out[34] += -0.7071067811865475*Ghat[19]*dv12; 
  out[35] += -1.224744871391589*Ghat[11]*dv12; 
  out[36] += -1.224744871391589*Ghat[12]*dv12; 
  out[37] += -1.224744871391589*Ghat[13]*dv12; 
  out[38] += -1.224744871391589*Ghat[14]*dv12; 
  out[39] += -1.224744871391589*Ghat[15]*dv12; 
  out[40] += -1.224744871391589*Ghat[16]*dv12; 
  out[41] += -1.58113883008419*Ghat[4]*dv12; 
  out[42] += -1.58113883008419*Ghat[5]*dv12; 
  out[43] += -1.58113883008419*Ghat[6]*dv12; 
  out[44] += -1.224744871391589*Ghat[17]*dv12; 
  out[45] += -1.224744871391589*Ghat[18]*dv12; 
  out[46] += -1.224744871391589*Ghat[19]*dv12; 
  out[47] += -1.58113883008419*Ghat[10]*dv12; 

  } else { 

  if (0.6363961030678926*alpha[4]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[0] = ser_4x_p2_surfx4_quad_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_4x_p2_surfx4_quad_0_l(fSkin); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[2] > 0) { 
    fUpwindQuad[1] = ser_4x_p2_surfx4_quad_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = ser_4x_p2_surfx4_quad_1_l(fSkin); 
  } 
  if ((-0.6363961030678926*alpha[4])-0.4743416490252568*alpha[2]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[2] = ser_4x_p2_surfx4_quad_2_r(fEdge); 
  } else { 
    fUpwindQuad[2] = ser_4x_p2_surfx4_quad_2_l(fSkin); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[1] > 0) { 
    fUpwindQuad[3] = ser_4x_p2_surfx4_quad_3_r(fEdge); 
  } else { 
    fUpwindQuad[3] = ser_4x_p2_surfx4_quad_3_l(fSkin); 
  } 
  if (0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[4] = ser_4x_p2_surfx4_quad_4_r(fEdge); 
  } else { 
    fUpwindQuad[4] = ser_4x_p2_surfx4_quad_4_l(fSkin); 
  } 
  if (0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[5] = ser_4x_p2_surfx4_quad_5_r(fEdge); 
  } else { 
    fUpwindQuad[5] = ser_4x_p2_surfx4_quad_5_l(fSkin); 
  } 
  if ((-0.6363961030678926*alpha[4])+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[6] = ser_4x_p2_surfx4_quad_6_r(fEdge); 
  } else { 
    fUpwindQuad[6] = ser_4x_p2_surfx4_quad_6_l(fSkin); 
  } 
  if (0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[7] = ser_4x_p2_surfx4_quad_7_r(fEdge); 
  } else { 
    fUpwindQuad[7] = ser_4x_p2_surfx4_quad_7_l(fSkin); 
  } 
  if (0.6363961030678926*alpha[4]+0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[8] = ser_4x_p2_surfx4_quad_8_r(fEdge); 
  } else { 
    fUpwindQuad[8] = ser_4x_p2_surfx4_quad_8_l(fSkin); 
  } 
  if (0.6363961030678926*alpha[4]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[9] = ser_4x_p2_surfx4_quad_9_r(fEdge); 
  } else { 
    fUpwindQuad[9] = ser_4x_p2_surfx4_quad_9_l(fSkin); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[2] > 0) { 
    fUpwindQuad[10] = ser_4x_p2_surfx4_quad_10_r(fEdge); 
  } else { 
    fUpwindQuad[10] = ser_4x_p2_surfx4_quad_10_l(fSkin); 
  } 
  if ((-0.6363961030678926*alpha[4])-0.4743416490252568*alpha[2]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[11] = ser_4x_p2_surfx4_quad_11_r(fEdge); 
  } else { 
    fUpwindQuad[11] = ser_4x_p2_surfx4_quad_11_l(fSkin); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[1] > 0) { 
    fUpwindQuad[12] = ser_4x_p2_surfx4_quad_12_r(fEdge); 
  } else { 
    fUpwindQuad[12] = ser_4x_p2_surfx4_quad_12_l(fSkin); 
  } 
  if (0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[13] = ser_4x_p2_surfx4_quad_13_r(fEdge); 
  } else { 
    fUpwindQuad[13] = ser_4x_p2_surfx4_quad_13_l(fSkin); 
  } 
  if (0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[14] = ser_4x_p2_surfx4_quad_14_r(fEdge); 
  } else { 
    fUpwindQuad[14] = ser_4x_p2_surfx4_quad_14_l(fSkin); 
  } 
  if ((-0.6363961030678926*alpha[4])+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[15] = ser_4x_p2_surfx4_quad_15_r(fEdge); 
  } else { 
    fUpwindQuad[15] = ser_4x_p2_surfx4_quad_15_l(fSkin); 
  } 
  if (0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[16] = ser_4x_p2_surfx4_quad_16_r(fEdge); 
  } else { 
    fUpwindQuad[16] = ser_4x_p2_surfx4_quad_16_l(fSkin); 
  } 
  if (0.6363961030678926*alpha[4]+0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[17] = ser_4x_p2_surfx4_quad_17_r(fEdge); 
  } else { 
    fUpwindQuad[17] = ser_4x_p2_surfx4_quad_17_l(fSkin); 
  } 
  if (0.6363961030678926*alpha[4]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[18] = ser_4x_p2_surfx4_quad_18_r(fEdge); 
  } else { 
    fUpwindQuad[18] = ser_4x_p2_surfx4_quad_18_l(fSkin); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[2] > 0) { 
    fUpwindQuad[19] = ser_4x_p2_surfx4_quad_19_r(fEdge); 
  } else { 
    fUpwindQuad[19] = ser_4x_p2_surfx4_quad_19_l(fSkin); 
  } 
  if ((-0.6363961030678926*alpha[4])-0.4743416490252568*alpha[2]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[20] = ser_4x_p2_surfx4_quad_20_r(fEdge); 
  } else { 
    fUpwindQuad[20] = ser_4x_p2_surfx4_quad_20_l(fSkin); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[1] > 0) { 
    fUpwindQuad[21] = ser_4x_p2_surfx4_quad_21_r(fEdge); 
  } else { 
    fUpwindQuad[21] = ser_4x_p2_surfx4_quad_21_l(fSkin); 
  } 
  if (0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[22] = ser_4x_p2_surfx4_quad_22_r(fEdge); 
  } else { 
    fUpwindQuad[22] = ser_4x_p2_surfx4_quad_22_l(fSkin); 
  } 
  if (0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[23] = ser_4x_p2_surfx4_quad_23_r(fEdge); 
  } else { 
    fUpwindQuad[23] = ser_4x_p2_surfx4_quad_23_l(fSkin); 
  } 
  if ((-0.6363961030678926*alpha[4])+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[24] = ser_4x_p2_surfx4_quad_24_r(fEdge); 
  } else { 
    fUpwindQuad[24] = ser_4x_p2_surfx4_quad_24_l(fSkin); 
  } 
  if (0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[25] = ser_4x_p2_surfx4_quad_25_r(fEdge); 
  } else { 
    fUpwindQuad[25] = ser_4x_p2_surfx4_quad_25_l(fSkin); 
  } 
  if (0.6363961030678926*alpha[4]+0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[26] = ser_4x_p2_surfx4_quad_26_r(fEdge); 
  } else { 
    fUpwindQuad[26] = ser_4x_p2_surfx4_quad_26_l(fSkin); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_4x_p2_upwind(fUpwindQuad, fUpwind); 

  Ghat[0] += 0.3535533905932737*(alpha[4]*fUpwind[4]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] += 0.3162277660168379*alpha[4]*fUpwind[11]+0.3162277660168379*alpha[1]*fUpwind[7]+0.3535533905932737*(alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] += 0.3162277660168379*alpha[4]*fUpwind[12]+0.3162277660168379*alpha[2]*fUpwind[8]+0.3535533905932737*(alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] += 0.3535533905932737*(alpha[4]*fUpwind[10]+alpha[2]*fUpwind[6]+alpha[1]*fUpwind[5]+alpha[0]*fUpwind[3]); 
  Ghat[4] += 0.3162277660168379*(alpha[2]*fUpwind[12]+alpha[1]*fUpwind[11])+0.3162277660168379*alpha[4]*(fUpwind[8]+fUpwind[7])+0.3535533905932737*(alpha[0]*fUpwind[4]+fUpwind[0]*alpha[4]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[5] += 0.3162277660168379*alpha[4]*fUpwind[17]+0.3162277660168379*alpha[1]*fUpwind[13]+0.3535533905932737*(alpha[2]*fUpwind[10]+alpha[4]*fUpwind[6]+alpha[0]*fUpwind[5]+alpha[1]*fUpwind[3]); 
  Ghat[6] += 0.3162277660168379*alpha[4]*fUpwind[18]+0.3162277660168379*alpha[2]*fUpwind[14]+0.3535533905932737*(alpha[1]*fUpwind[10]+alpha[0]*fUpwind[6]+alpha[4]*fUpwind[5]+alpha[2]*fUpwind[3]); 
  Ghat[7] += 0.3535533905932737*(alpha[2]*fUpwind[11]+alpha[0]*fUpwind[7])+0.3162277660168379*(alpha[4]*fUpwind[4]+alpha[1]*fUpwind[1]); 
  Ghat[8] += 0.3535533905932737*(alpha[1]*fUpwind[12]+alpha[0]*fUpwind[8])+0.3162277660168379*(alpha[4]*fUpwind[4]+alpha[2]*fUpwind[2]); 
  Ghat[9] += 0.3535533905932737*(alpha[4]*fUpwind[19]+alpha[2]*fUpwind[16]+alpha[1]*fUpwind[15]+alpha[0]*fUpwind[9]); 
  Ghat[10] += 0.3162277660168379*(alpha[2]*fUpwind[18]+alpha[1]*fUpwind[17])+0.3162277660168379*alpha[4]*(fUpwind[14]+fUpwind[13])+0.3535533905932737*(alpha[0]*fUpwind[10]+alpha[1]*fUpwind[6]+alpha[2]*fUpwind[5]+fUpwind[3]*alpha[4]); 
  Ghat[11] += 0.2828427124746191*alpha[4]*fUpwind[12]+0.3535533905932737*(alpha[0]*fUpwind[11]+alpha[2]*fUpwind[7])+0.3162277660168379*(alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]); 
  Ghat[12] += 0.3535533905932737*alpha[0]*fUpwind[12]+0.2828427124746191*alpha[4]*fUpwind[11]+0.3535533905932737*alpha[1]*fUpwind[8]+0.3162277660168379*(alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]); 
  Ghat[13] += 0.3535533905932737*(alpha[2]*fUpwind[17]+alpha[0]*fUpwind[13])+0.3162277660168379*(alpha[4]*fUpwind[10]+alpha[1]*fUpwind[5]); 
  Ghat[14] += 0.3535533905932737*(alpha[1]*fUpwind[18]+alpha[0]*fUpwind[14])+0.3162277660168379*(alpha[4]*fUpwind[10]+alpha[2]*fUpwind[6]); 
  Ghat[15] += 0.3535533905932737*(alpha[2]*fUpwind[19]+alpha[4]*fUpwind[16]+alpha[0]*fUpwind[15]+alpha[1]*fUpwind[9]); 
  Ghat[16] += 0.3535533905932737*(alpha[1]*fUpwind[19]+alpha[0]*fUpwind[16]+alpha[4]*fUpwind[15]+alpha[2]*fUpwind[9]); 
  Ghat[17] += 0.2828427124746191*alpha[4]*fUpwind[18]+0.3535533905932737*(alpha[0]*fUpwind[17]+alpha[2]*fUpwind[13])+0.3162277660168379*(alpha[1]*fUpwind[10]+alpha[4]*fUpwind[5]); 
  Ghat[18] += 0.3535533905932737*alpha[0]*fUpwind[18]+0.2828427124746191*alpha[4]*fUpwind[17]+0.3535533905932737*alpha[1]*fUpwind[14]+0.3162277660168379*(alpha[2]*fUpwind[10]+alpha[4]*fUpwind[6]); 
  Ghat[19] += 0.3535533905932737*(alpha[0]*fUpwind[19]+alpha[1]*fUpwind[16]+alpha[2]*fUpwind[15]+alpha[4]*fUpwind[9]); 

  out[0] += 0.7071067811865475*Ghat[0]*dv12; 
  out[1] += 0.7071067811865475*Ghat[1]*dv12; 
  out[2] += 0.7071067811865475*Ghat[2]*dv12; 
  out[3] += 0.7071067811865475*Ghat[3]*dv12; 
  out[4] += -1.224744871391589*Ghat[0]*dv12; 
  out[5] += 0.7071067811865475*Ghat[4]*dv12; 
  out[6] += 0.7071067811865475*Ghat[5]*dv12; 
  out[7] += 0.7071067811865475*Ghat[6]*dv12; 
  out[8] += -1.224744871391589*Ghat[1]*dv12; 
  out[9] += -1.224744871391589*Ghat[2]*dv12; 
  out[10] += -1.224744871391589*Ghat[3]*dv12; 
  out[11] += 0.7071067811865475*Ghat[7]*dv12; 
  out[12] += 0.7071067811865475*Ghat[8]*dv12; 
  out[13] += 0.7071067811865475*Ghat[9]*dv12; 
  out[14] += 1.58113883008419*Ghat[0]*dv12; 
  out[15] += 0.7071067811865475*Ghat[10]*dv12; 
  out[16] += -1.224744871391589*Ghat[4]*dv12; 
  out[17] += -1.224744871391589*Ghat[5]*dv12; 
  out[18] += -1.224744871391589*Ghat[6]*dv12; 
  out[19] += 0.7071067811865475*Ghat[11]*dv12; 
  out[20] += 0.7071067811865475*Ghat[12]*dv12; 
  out[21] += 0.7071067811865475*Ghat[13]*dv12; 
  out[22] += 0.7071067811865475*Ghat[14]*dv12; 
  out[23] += 0.7071067811865475*Ghat[15]*dv12; 
  out[24] += 0.7071067811865475*Ghat[16]*dv12; 
  out[25] += -1.224744871391589*Ghat[7]*dv12; 
  out[26] += -1.224744871391589*Ghat[8]*dv12; 
  out[27] += -1.224744871391589*Ghat[9]*dv12; 
  out[28] += 1.58113883008419*Ghat[1]*dv12; 
  out[29] += 1.58113883008419*Ghat[2]*dv12; 
  out[30] += 1.58113883008419*Ghat[3]*dv12; 
  out[31] += -1.224744871391589*Ghat[10]*dv12; 
  out[32] += 0.7071067811865475*Ghat[17]*dv12; 
  out[33] += 0.7071067811865475*Ghat[18]*dv12; 
  out[34] += 0.7071067811865475*Ghat[19]*dv12; 
  out[35] += -1.224744871391589*Ghat[11]*dv12; 
  out[36] += -1.224744871391589*Ghat[12]*dv12; 
  out[37] += -1.224744871391589*Ghat[13]*dv12; 
  out[38] += -1.224744871391589*Ghat[14]*dv12; 
  out[39] += -1.224744871391589*Ghat[15]*dv12; 
  out[40] += -1.224744871391589*Ghat[16]*dv12; 
  out[41] += 1.58113883008419*Ghat[4]*dv12; 
  out[42] += 1.58113883008419*Ghat[5]*dv12; 
  out[43] += 1.58113883008419*Ghat[6]*dv12; 
  out[44] += -1.224744871391589*Ghat[17]*dv12; 
  out[45] += -1.224744871391589*Ghat[18]*dv12; 
  out[46] += -1.224744871391589*Ghat[19]*dv12; 
  out[47] += 1.58113883008419*Ghat[10]*dv12; 

  } 
} 
