#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_5x_p2_surfx3_quad.h> 
#include <gkyl_basis_ser_5x_p2_upwind.h> 
GKYL_CU_DH void vlasov_boundary_surfvx_2x3v_ser_p2(const double *w, const double *dxv, const double *qmem, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // qmem:        q/m*EM fields.
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv10 = 2/dxv[2]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double dv3 = dxv[4], wv3 = w[4]; 
  const double *E0 = &qmem[0]; 
  const double *B0 = &qmem[24]; 
  const double *B1 = &qmem[32]; 
  const double *B2 = &qmem[40]; 

  double alpha[48] = {0.0}; 

  alpha[0] = 2.0*(B2[0]*wv2+E0[0])-2.0*B1[0]*wv3; 
  alpha[1] = 2.0*(B2[1]*wv2+E0[1])-2.0*B1[1]*wv3; 
  alpha[2] = 2.0*(B2[2]*wv2+E0[2])-2.0*B1[2]*wv3; 
  alpha[3] = 0.5773502691896258*B2[0]*dv2; 
  alpha[4] = -0.5773502691896258*B1[0]*dv3; 
  alpha[5] = 2.0*(B2[3]*wv2+E0[3])-2.0*B1[3]*wv3; 
  alpha[6] = 0.5773502691896258*B2[1]*dv2; 
  alpha[7] = 0.5773502691896258*B2[2]*dv2; 
  alpha[8] = -0.5773502691896258*B1[1]*dv3; 
  alpha[9] = -0.5773502691896258*B1[2]*dv3; 
  alpha[11] = 2.0*(B2[4]*wv2+E0[4])-2.0*B1[4]*wv3; 
  alpha[12] = 2.0*(B2[5]*wv2+E0[5])-2.0*B1[5]*wv3; 
  alpha[15] = 0.5773502691896258*B2[3]*dv2; 
  alpha[16] = -0.5773502691896258*B1[3]*dv3; 
  alpha[19] = 2.0*(B2[6]*wv2+E0[6])-2.0*B1[6]*wv3; 
  alpha[20] = 2.0*(B2[7]*wv2+E0[7])-2.0*B1[7]*wv3; 
  alpha[21] = 0.5773502691896258*B2[4]*dv2; 
  alpha[22] = 0.5773502691896258*B2[5]*dv2; 
  alpha[25] = -0.5773502691896258*B1[4]*dv3; 
  alpha[26] = -0.5773502691896258*B1[5]*dv3; 
  alpha[32] = 0.5773502691896258*B2[6]*dv2; 
  alpha[33] = 0.5773502691896258*B2[7]*dv2; 
  alpha[35] = -0.5773502691896258*B1[6]*dv3; 
  alpha[36] = -0.5773502691896258*B1[7]*dv3; 

  double fUpwindQuad[81] = {0.0};
  double fUpwind[48] = {0.0};
  double Ghat[48] = {0.0}; 

  if (edge == -1) { 

  if (0.4024922359499621*(alpha[36]+alpha[35]+alpha[33]+alpha[32])-0.2999999999999999*(alpha[26]+alpha[25]+alpha[22]+alpha[21])-0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]-0.6037383539249431*(alpha[16]+alpha[15])+0.223606797749979*(alpha[12]+alpha[11])+0.45*(alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5])-0.3354101966249685*(alpha[4]+alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[0] = ser_5x_p2_surfx3_quad_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_5x_p2_surfx3_quad_0_l(fEdge); 
  } 
  if ((-0.5031152949374527*(alpha[35]+alpha[32]))-0.2999999999999999*alpha[26]+0.375*alpha[25]-0.2999999999999999*alpha[22]+0.375*(alpha[21]+alpha[19])+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]+0.45*(alpha[9]+alpha[7])-0.3354101966249685*(alpha[4]+alpha[3]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad[1] = ser_5x_p2_surfx3_quad_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = ser_5x_p2_surfx3_quad_1_l(fEdge); 
  } 
  if ((-0.4024922359499621*alpha[36])+0.4024922359499621*alpha[35]-0.4024922359499621*alpha[33]+0.4024922359499621*alpha[32]-0.2999999999999999*(alpha[26]+alpha[25]+alpha[22]+alpha[21])+0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]+0.6037383539249431*(alpha[16]+alpha[15])+0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[9]-0.45*alpha[8]+0.45*alpha[7]-0.45*(alpha[6]+alpha[5])-0.3354101966249685*(alpha[4]+alpha[3]+alpha[2])+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[2] = ser_5x_p2_surfx3_quad_2_r(fSkin); 
  } else { 
    fUpwindQuad[2] = ser_5x_p2_surfx3_quad_2_l(fEdge); 
  } 
  if ((-0.5031152949374527*(alpha[36]+alpha[33]))+0.375*alpha[26]-0.2999999999999999*alpha[25]+0.375*alpha[22]-0.2999999999999999*alpha[21]+0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]+0.45*(alpha[8]+alpha[6])-0.3354101966249685*(alpha[4]+alpha[3]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[3] = ser_5x_p2_surfx3_quad_3_r(fSkin); 
  } else { 
    fUpwindQuad[3] = ser_5x_p2_surfx3_quad_3_l(fEdge); 
  } 
  if (0.375*(alpha[26]+alpha[25]+alpha[22]+alpha[21])-0.2795084971874737*(alpha[12]+alpha[11])-0.3354101966249685*(alpha[4]+alpha[3])+0.25*alpha[0] > 0) { 
    fUpwindQuad[4] = ser_5x_p2_surfx3_quad_4_r(fSkin); 
  } else { 
    fUpwindQuad[4] = ser_5x_p2_surfx3_quad_4_l(fEdge); 
  } 
  if (0.5031152949374527*(alpha[36]+alpha[33])+0.375*alpha[26]-0.2999999999999999*alpha[25]+0.375*alpha[22]-0.2999999999999999*alpha[21]-0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]-0.45*(alpha[8]+alpha[6])-0.3354101966249685*(alpha[4]+alpha[3])+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[5] = ser_5x_p2_surfx3_quad_5_r(fSkin); 
  } else { 
    fUpwindQuad[5] = ser_5x_p2_surfx3_quad_5_l(fEdge); 
  } 
  if (0.4024922359499621*alpha[36]-0.4024922359499621*alpha[35]+0.4024922359499621*alpha[33]-0.4024922359499621*alpha[32]-0.2999999999999999*(alpha[26]+alpha[25]+alpha[22]+alpha[21])-0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]+0.6037383539249431*(alpha[16]+alpha[15])+0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[9]+0.45*alpha[8]-0.45*alpha[7]+0.45*alpha[6]-0.45*alpha[5]-0.3354101966249685*(alpha[4]+alpha[3])+0.3354101966249685*alpha[2]-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[6] = ser_5x_p2_surfx3_quad_6_r(fSkin); 
  } else { 
    fUpwindQuad[6] = ser_5x_p2_surfx3_quad_6_l(fEdge); 
  } 
  if (0.5031152949374527*(alpha[35]+alpha[32])-0.2999999999999999*alpha[26]+0.375*alpha[25]-0.2999999999999999*alpha[22]+0.375*alpha[21]-0.375*alpha[19]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]-0.45*(alpha[9]+alpha[7])-0.3354101966249685*(alpha[4]+alpha[3])+0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[7] = ser_5x_p2_surfx3_quad_7_r(fSkin); 
  } else { 
    fUpwindQuad[7] = ser_5x_p2_surfx3_quad_7_l(fEdge); 
  } 
  if ((-0.4024922359499621*(alpha[36]+alpha[35]+alpha[33]+alpha[32]))-0.2999999999999999*(alpha[26]+alpha[25]+alpha[22]+alpha[21])+0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]-0.6037383539249431*(alpha[16]+alpha[15])+0.223606797749979*(alpha[12]+alpha[11])-0.45*(alpha[9]+alpha[8]+alpha[7]+alpha[6])+0.45*alpha[5]-0.3354101966249685*(alpha[4]+alpha[3])+0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[8] = ser_5x_p2_surfx3_quad_8_r(fSkin); 
  } else { 
    fUpwindQuad[8] = ser_5x_p2_surfx3_quad_8_l(fEdge); 
  } 
  if (0.4024922359499621*(alpha[36]+alpha[35])-0.2999999999999999*(alpha[26]+alpha[25])-0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]-0.6037383539249431*alpha[16]+0.223606797749979*(alpha[12]+alpha[11])+0.45*(alpha[9]+alpha[8]+alpha[5])-0.3354101966249685*(alpha[4]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[9] = ser_5x_p2_surfx3_quad_9_r(fSkin); 
  } else { 
    fUpwindQuad[9] = ser_5x_p2_surfx3_quad_9_l(fEdge); 
  } 
  if ((-0.5031152949374527*alpha[35])-0.2999999999999999*alpha[26]+0.375*(alpha[25]+alpha[19])+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]+0.45*alpha[9]-0.3354101966249685*(alpha[4]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad[10] = ser_5x_p2_surfx3_quad_10_r(fSkin); 
  } else { 
    fUpwindQuad[10] = ser_5x_p2_surfx3_quad_10_l(fEdge); 
  } 
  if ((-0.4024922359499621*alpha[36])+0.4024922359499621*alpha[35]-0.2999999999999999*(alpha[26]+alpha[25])+0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]+0.6037383539249431*alpha[16]+0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[9]-0.45*(alpha[8]+alpha[5])-0.3354101966249685*(alpha[4]+alpha[2])+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[11] = ser_5x_p2_surfx3_quad_11_r(fSkin); 
  } else { 
    fUpwindQuad[11] = ser_5x_p2_surfx3_quad_11_l(fEdge); 
  } 
  if ((-0.5031152949374527*alpha[36])+0.375*alpha[26]-0.2999999999999999*alpha[25]+0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]+0.45*alpha[8]-0.3354101966249685*(alpha[4]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[12] = ser_5x_p2_surfx3_quad_12_r(fSkin); 
  } else { 
    fUpwindQuad[12] = ser_5x_p2_surfx3_quad_12_l(fEdge); 
  } 
  if (0.375*(alpha[26]+alpha[25])-0.2795084971874737*(alpha[12]+alpha[11])-0.3354101966249685*alpha[4]+0.25*alpha[0] > 0) { 
    fUpwindQuad[13] = ser_5x_p2_surfx3_quad_13_r(fSkin); 
  } else { 
    fUpwindQuad[13] = ser_5x_p2_surfx3_quad_13_l(fEdge); 
  } 
  if (0.5031152949374527*alpha[36]+0.375*alpha[26]-0.2999999999999999*alpha[25]-0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]-0.45*alpha[8]-0.3354101966249685*alpha[4]+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[14] = ser_5x_p2_surfx3_quad_14_r(fSkin); 
  } else { 
    fUpwindQuad[14] = ser_5x_p2_surfx3_quad_14_l(fEdge); 
  } 
  if (0.4024922359499621*alpha[36]-0.4024922359499621*alpha[35]-0.2999999999999999*(alpha[26]+alpha[25])-0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]+0.6037383539249431*alpha[16]+0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[9]+0.45*alpha[8]-0.45*alpha[5]-0.3354101966249685*alpha[4]+0.3354101966249685*alpha[2]-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[15] = ser_5x_p2_surfx3_quad_15_r(fSkin); 
  } else { 
    fUpwindQuad[15] = ser_5x_p2_surfx3_quad_15_l(fEdge); 
  } 
  if (0.5031152949374527*alpha[35]-0.2999999999999999*alpha[26]+0.375*alpha[25]-0.375*alpha[19]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]-0.45*alpha[9]-0.3354101966249685*alpha[4]+0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[16] = ser_5x_p2_surfx3_quad_16_r(fSkin); 
  } else { 
    fUpwindQuad[16] = ser_5x_p2_surfx3_quad_16_l(fEdge); 
  } 
  if ((-0.4024922359499621*(alpha[36]+alpha[35]))-0.2999999999999999*(alpha[26]+alpha[25])+0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]-0.6037383539249431*alpha[16]+0.223606797749979*(alpha[12]+alpha[11])-0.45*(alpha[9]+alpha[8])+0.45*alpha[5]-0.3354101966249685*alpha[4]+0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[17] = ser_5x_p2_surfx3_quad_17_r(fSkin); 
  } else { 
    fUpwindQuad[17] = ser_5x_p2_surfx3_quad_17_l(fEdge); 
  } 
  if (0.4024922359499621*(alpha[36]+alpha[35])-0.4024922359499621*(alpha[33]+alpha[32])-0.2999999999999999*(alpha[26]+alpha[25])+0.2999999999999999*(alpha[22]+alpha[21])-0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]-0.6037383539249431*alpha[16]+0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])+0.45*(alpha[9]+alpha[8])-0.45*(alpha[7]+alpha[6])+0.45*alpha[5]-0.3354101966249685*alpha[4]+0.3354101966249685*alpha[3]-0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[18] = ser_5x_p2_surfx3_quad_18_r(fSkin); 
  } else { 
    fUpwindQuad[18] = ser_5x_p2_surfx3_quad_18_l(fEdge); 
  } 
  if ((-0.5031152949374527*alpha[35])+0.5031152949374527*alpha[32]-0.2999999999999999*alpha[26]+0.375*alpha[25]+0.2999999999999999*alpha[22]-0.375*alpha[21]+0.375*alpha[19]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]+0.45*alpha[9]-0.45*alpha[7]-0.3354101966249685*alpha[4]+0.3354101966249685*alpha[3]-0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[19] = ser_5x_p2_surfx3_quad_19_r(fSkin); 
  } else { 
    fUpwindQuad[19] = ser_5x_p2_surfx3_quad_19_l(fEdge); 
  } 
  if ((-0.4024922359499621*alpha[36])+0.4024922359499621*(alpha[35]+alpha[33])-0.4024922359499621*alpha[32]-0.2999999999999999*(alpha[26]+alpha[25])+0.2999999999999999*(alpha[22]+alpha[21])+0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]+0.6037383539249431*alpha[16]-0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[9]-0.45*(alpha[8]+alpha[7])+0.45*alpha[6]-0.45*alpha[5]-0.3354101966249685*alpha[4]+0.3354101966249685*alpha[3]-0.3354101966249685*alpha[2]+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[20] = ser_5x_p2_surfx3_quad_20_r(fSkin); 
  } else { 
    fUpwindQuad[20] = ser_5x_p2_surfx3_quad_20_l(fEdge); 
  } 
  if ((-0.5031152949374527*alpha[36])+0.5031152949374527*alpha[33]+0.375*alpha[26]-0.2999999999999999*alpha[25]-0.375*alpha[22]+0.2999999999999999*alpha[21]+0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]+0.45*alpha[8]-0.45*alpha[6]-0.3354101966249685*alpha[4]+0.3354101966249685*alpha[3]-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[21] = ser_5x_p2_surfx3_quad_21_r(fSkin); 
  } else { 
    fUpwindQuad[21] = ser_5x_p2_surfx3_quad_21_l(fEdge); 
  } 
  if (0.375*(alpha[26]+alpha[25])-0.375*(alpha[22]+alpha[21])-0.2795084971874737*(alpha[12]+alpha[11])-0.3354101966249685*alpha[4]+0.3354101966249685*alpha[3]+0.25*alpha[0] > 0) { 
    fUpwindQuad[22] = ser_5x_p2_surfx3_quad_22_r(fSkin); 
  } else { 
    fUpwindQuad[22] = ser_5x_p2_surfx3_quad_22_l(fEdge); 
  } 
  if (0.5031152949374527*alpha[36]-0.5031152949374527*alpha[33]+0.375*alpha[26]-0.2999999999999999*alpha[25]-0.375*alpha[22]+0.2999999999999999*alpha[21]-0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]-0.45*alpha[8]+0.45*alpha[6]-0.3354101966249685*alpha[4]+0.3354101966249685*(alpha[3]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[23] = ser_5x_p2_surfx3_quad_23_r(fSkin); 
  } else { 
    fUpwindQuad[23] = ser_5x_p2_surfx3_quad_23_l(fEdge); 
  } 
  if (0.4024922359499621*alpha[36]-0.4024922359499621*(alpha[35]+alpha[33])+0.4024922359499621*alpha[32]-0.2999999999999999*(alpha[26]+alpha[25])+0.2999999999999999*(alpha[22]+alpha[21])-0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]+0.6037383539249431*alpha[16]-0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[9]+0.45*(alpha[8]+alpha[7])-0.45*(alpha[6]+alpha[5])-0.3354101966249685*alpha[4]+0.3354101966249685*(alpha[3]+alpha[2])-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[24] = ser_5x_p2_surfx3_quad_24_r(fSkin); 
  } else { 
    fUpwindQuad[24] = ser_5x_p2_surfx3_quad_24_l(fEdge); 
  } 
  if (0.5031152949374527*alpha[35]-0.5031152949374527*alpha[32]-0.2999999999999999*alpha[26]+0.375*alpha[25]+0.2999999999999999*alpha[22]-0.375*(alpha[21]+alpha[19])+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]-0.45*alpha[9]+0.45*alpha[7]-0.3354101966249685*alpha[4]+0.3354101966249685*(alpha[3]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad[25] = ser_5x_p2_surfx3_quad_25_r(fSkin); 
  } else { 
    fUpwindQuad[25] = ser_5x_p2_surfx3_quad_25_l(fEdge); 
  } 
  if ((-0.4024922359499621*(alpha[36]+alpha[35]))+0.4024922359499621*(alpha[33]+alpha[32])-0.2999999999999999*(alpha[26]+alpha[25])+0.2999999999999999*(alpha[22]+alpha[21])+0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]-0.6037383539249431*alpha[16]+0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])-0.45*(alpha[9]+alpha[8])+0.45*(alpha[7]+alpha[6]+alpha[5])-0.3354101966249685*alpha[4]+0.3354101966249685*(alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[26] = ser_5x_p2_surfx3_quad_26_r(fSkin); 
  } else { 
    fUpwindQuad[26] = ser_5x_p2_surfx3_quad_26_l(fEdge); 
  } 
  if (0.4024922359499621*(alpha[33]+alpha[32])-0.2999999999999999*(alpha[22]+alpha[21])-0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]-0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])+0.45*(alpha[7]+alpha[6]+alpha[5])-0.3354101966249685*(alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[27] = ser_5x_p2_surfx3_quad_27_r(fSkin); 
  } else { 
    fUpwindQuad[27] = ser_5x_p2_surfx3_quad_27_l(fEdge); 
  } 
  if ((-0.5031152949374527*alpha[32])-0.2999999999999999*alpha[22]+0.375*(alpha[21]+alpha[19])+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]+0.45*alpha[7]-0.3354101966249685*(alpha[3]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad[28] = ser_5x_p2_surfx3_quad_28_r(fSkin); 
  } else { 
    fUpwindQuad[28] = ser_5x_p2_surfx3_quad_28_l(fEdge); 
  } 
  if ((-0.4024922359499621*alpha[33])+0.4024922359499621*alpha[32]-0.2999999999999999*(alpha[22]+alpha[21])+0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]+0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[7]-0.45*(alpha[6]+alpha[5])-0.3354101966249685*(alpha[3]+alpha[2])+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[29] = ser_5x_p2_surfx3_quad_29_r(fSkin); 
  } else { 
    fUpwindQuad[29] = ser_5x_p2_surfx3_quad_29_l(fEdge); 
  } 
  if ((-0.5031152949374527*alpha[33])+0.375*alpha[22]-0.2999999999999999*alpha[21]+0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]+0.45*alpha[6]-0.3354101966249685*(alpha[3]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[30] = ser_5x_p2_surfx3_quad_30_r(fSkin); 
  } else { 
    fUpwindQuad[30] = ser_5x_p2_surfx3_quad_30_l(fEdge); 
  } 
  if (0.375*(alpha[22]+alpha[21])-0.2795084971874737*(alpha[12]+alpha[11])-0.3354101966249685*alpha[3]+0.25*alpha[0] > 0) { 
    fUpwindQuad[31] = ser_5x_p2_surfx3_quad_31_r(fSkin); 
  } else { 
    fUpwindQuad[31] = ser_5x_p2_surfx3_quad_31_l(fEdge); 
  } 
  if (0.5031152949374527*alpha[33]+0.375*alpha[22]-0.2999999999999999*alpha[21]-0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]-0.45*alpha[6]-0.3354101966249685*alpha[3]+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[32] = ser_5x_p2_surfx3_quad_32_r(fSkin); 
  } else { 
    fUpwindQuad[32] = ser_5x_p2_surfx3_quad_32_l(fEdge); 
  } 
  if (0.4024922359499621*alpha[33]-0.4024922359499621*alpha[32]-0.2999999999999999*(alpha[22]+alpha[21])-0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]+0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[7]+0.45*alpha[6]-0.45*alpha[5]-0.3354101966249685*alpha[3]+0.3354101966249685*alpha[2]-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[33] = ser_5x_p2_surfx3_quad_33_r(fSkin); 
  } else { 
    fUpwindQuad[33] = ser_5x_p2_surfx3_quad_33_l(fEdge); 
  } 
  if (0.5031152949374527*alpha[32]-0.2999999999999999*alpha[22]+0.375*alpha[21]-0.375*alpha[19]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]-0.45*alpha[7]-0.3354101966249685*alpha[3]+0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[34] = ser_5x_p2_surfx3_quad_34_r(fSkin); 
  } else { 
    fUpwindQuad[34] = ser_5x_p2_surfx3_quad_34_l(fEdge); 
  } 
  if ((-0.4024922359499621*(alpha[33]+alpha[32]))-0.2999999999999999*(alpha[22]+alpha[21])+0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]-0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])-0.45*(alpha[7]+alpha[6])+0.45*alpha[5]-0.3354101966249685*alpha[3]+0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[35] = ser_5x_p2_surfx3_quad_35_r(fSkin); 
  } else { 
    fUpwindQuad[35] = ser_5x_p2_surfx3_quad_35_l(fEdge); 
  } 
  if ((-0.2999999999999998*alpha[20])-0.2999999999999999*alpha[19]+0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[5]-0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[36] = ser_5x_p2_surfx3_quad_36_r(fSkin); 
  } else { 
    fUpwindQuad[36] = ser_5x_p2_surfx3_quad_36_l(fEdge); 
  } 
  if (0.375*alpha[19]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]-0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[37] = ser_5x_p2_surfx3_quad_37_r(fSkin); 
  } else { 
    fUpwindQuad[37] = ser_5x_p2_surfx3_quad_37_l(fEdge); 
  } 
  if (0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]+0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[5]-0.3354101966249685*alpha[2]+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[38] = ser_5x_p2_surfx3_quad_38_r(fSkin); 
  } else { 
    fUpwindQuad[38] = ser_5x_p2_surfx3_quad_38_l(fEdge); 
  } 
  if (0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[39] = ser_5x_p2_surfx3_quad_39_r(fSkin); 
  } else { 
    fUpwindQuad[39] = ser_5x_p2_surfx3_quad_39_l(fEdge); 
  } 
  if (0.25*alpha[0]-0.2795084971874737*(alpha[12]+alpha[11]) > 0) { 
    fUpwindQuad[40] = ser_5x_p2_surfx3_quad_40_r(fSkin); 
  } else { 
    fUpwindQuad[40] = ser_5x_p2_surfx3_quad_40_l(fEdge); 
  } 
  if ((-0.375*alpha[20])-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[41] = ser_5x_p2_surfx3_quad_41_r(fSkin); 
  } else { 
    fUpwindQuad[41] = ser_5x_p2_surfx3_quad_41_l(fEdge); 
  } 
  if ((-0.2999999999999998*alpha[20])+0.2999999999999999*alpha[19]+0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[5]+0.3354101966249685*alpha[2]-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[42] = ser_5x_p2_surfx3_quad_42_r(fSkin); 
  } else { 
    fUpwindQuad[42] = ser_5x_p2_surfx3_quad_42_l(fEdge); 
  } 
  if ((-0.375*alpha[19])+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]+0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[43] = ser_5x_p2_surfx3_quad_43_r(fSkin); 
  } else { 
    fUpwindQuad[43] = ser_5x_p2_surfx3_quad_43_l(fEdge); 
  } 
  if (0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]+0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[5]+0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[44] = ser_5x_p2_surfx3_quad_44_r(fSkin); 
  } else { 
    fUpwindQuad[44] = ser_5x_p2_surfx3_quad_44_l(fEdge); 
  } 
  if ((-0.4024922359499621*(alpha[33]+alpha[32]))+0.2999999999999999*(alpha[22]+alpha[21])-0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]+0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])-0.45*(alpha[7]+alpha[6])+0.45*alpha[5]+0.3354101966249685*alpha[3]-0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[45] = ser_5x_p2_surfx3_quad_45_r(fSkin); 
  } else { 
    fUpwindQuad[45] = ser_5x_p2_surfx3_quad_45_l(fEdge); 
  } 
  if (0.5031152949374527*alpha[32]+0.2999999999999999*alpha[22]-0.375*alpha[21]+0.375*alpha[19]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]-0.45*alpha[7]+0.3354101966249685*alpha[3]-0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[46] = ser_5x_p2_surfx3_quad_46_r(fSkin); 
  } else { 
    fUpwindQuad[46] = ser_5x_p2_surfx3_quad_46_l(fEdge); 
  } 
  if (0.4024922359499621*alpha[33]-0.4024922359499621*alpha[32]+0.2999999999999999*(alpha[22]+alpha[21])+0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]-0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[7]+0.45*alpha[6]-0.45*alpha[5]+0.3354101966249685*alpha[3]-0.3354101966249685*alpha[2]+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[47] = ser_5x_p2_surfx3_quad_47_r(fSkin); 
  } else { 
    fUpwindQuad[47] = ser_5x_p2_surfx3_quad_47_l(fEdge); 
  } 
  if (0.5031152949374527*alpha[33]-0.375*alpha[22]+0.2999999999999999*alpha[21]+0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]-0.45*alpha[6]+0.3354101966249685*alpha[3]-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[48] = ser_5x_p2_surfx3_quad_48_r(fSkin); 
  } else { 
    fUpwindQuad[48] = ser_5x_p2_surfx3_quad_48_l(fEdge); 
  } 
  if ((-0.375*(alpha[22]+alpha[21]))-0.2795084971874737*(alpha[12]+alpha[11])+0.3354101966249685*alpha[3]+0.25*alpha[0] > 0) { 
    fUpwindQuad[49] = ser_5x_p2_surfx3_quad_49_r(fSkin); 
  } else { 
    fUpwindQuad[49] = ser_5x_p2_surfx3_quad_49_l(fEdge); 
  } 
  if ((-0.5031152949374527*alpha[33])-0.375*alpha[22]+0.2999999999999999*alpha[21]-0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]+0.45*alpha[6]+0.3354101966249685*(alpha[3]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[50] = ser_5x_p2_surfx3_quad_50_r(fSkin); 
  } else { 
    fUpwindQuad[50] = ser_5x_p2_surfx3_quad_50_l(fEdge); 
  } 
  if ((-0.4024922359499621*alpha[33])+0.4024922359499621*alpha[32]+0.2999999999999999*(alpha[22]+alpha[21])-0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]-0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[7]-0.45*(alpha[6]+alpha[5])+0.3354101966249685*(alpha[3]+alpha[2])-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[51] = ser_5x_p2_surfx3_quad_51_r(fSkin); 
  } else { 
    fUpwindQuad[51] = ser_5x_p2_surfx3_quad_51_l(fEdge); 
  } 
  if ((-0.5031152949374527*alpha[32])+0.2999999999999999*alpha[22]-0.375*(alpha[21]+alpha[19])+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]+0.45*alpha[7]+0.3354101966249685*(alpha[3]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad[52] = ser_5x_p2_surfx3_quad_52_r(fSkin); 
  } else { 
    fUpwindQuad[52] = ser_5x_p2_surfx3_quad_52_l(fEdge); 
  } 
  if (0.4024922359499621*(alpha[33]+alpha[32])+0.2999999999999999*(alpha[22]+alpha[21])+0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]+0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])+0.45*(alpha[7]+alpha[6]+alpha[5])+0.3354101966249685*(alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[53] = ser_5x_p2_surfx3_quad_53_r(fSkin); 
  } else { 
    fUpwindQuad[53] = ser_5x_p2_surfx3_quad_53_l(fEdge); 
  } 
  if ((-0.4024922359499621*(alpha[36]+alpha[35]))+0.4024922359499621*(alpha[33]+alpha[32])+0.2999999999999999*(alpha[26]+alpha[25])-0.2999999999999999*(alpha[22]+alpha[21])-0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]+0.6037383539249431*alpha[16]-0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])-0.45*(alpha[9]+alpha[8])+0.45*(alpha[7]+alpha[6]+alpha[5])+0.3354101966249685*alpha[4]-0.3354101966249685*(alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[54] = ser_5x_p2_surfx3_quad_54_r(fSkin); 
  } else { 
    fUpwindQuad[54] = ser_5x_p2_surfx3_quad_54_l(fEdge); 
  } 
  if (0.5031152949374527*alpha[35]-0.5031152949374527*alpha[32]+0.2999999999999999*alpha[26]-0.375*alpha[25]-0.2999999999999999*alpha[22]+0.375*(alpha[21]+alpha[19])+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]-0.45*alpha[9]+0.45*alpha[7]+0.3354101966249685*alpha[4]-0.3354101966249685*(alpha[3]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad[55] = ser_5x_p2_surfx3_quad_55_r(fSkin); 
  } else { 
    fUpwindQuad[55] = ser_5x_p2_surfx3_quad_55_l(fEdge); 
  } 
  if (0.4024922359499621*alpha[36]-0.4024922359499621*(alpha[35]+alpha[33])+0.4024922359499621*alpha[32]+0.2999999999999999*(alpha[26]+alpha[25])-0.2999999999999999*(alpha[22]+alpha[21])+0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]-0.6037383539249431*alpha[16]+0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[9]+0.45*(alpha[8]+alpha[7])-0.45*(alpha[6]+alpha[5])+0.3354101966249685*alpha[4]-0.3354101966249685*(alpha[3]+alpha[2])+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[56] = ser_5x_p2_surfx3_quad_56_r(fSkin); 
  } else { 
    fUpwindQuad[56] = ser_5x_p2_surfx3_quad_56_l(fEdge); 
  } 
  if (0.5031152949374527*alpha[36]-0.5031152949374527*alpha[33]-0.375*alpha[26]+0.2999999999999999*alpha[25]+0.375*alpha[22]-0.2999999999999999*alpha[21]+0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]-0.45*alpha[8]+0.45*alpha[6]+0.3354101966249685*alpha[4]-0.3354101966249685*(alpha[3]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[57] = ser_5x_p2_surfx3_quad_57_r(fSkin); 
  } else { 
    fUpwindQuad[57] = ser_5x_p2_surfx3_quad_57_l(fEdge); 
  } 
  if ((-0.375*(alpha[26]+alpha[25]))+0.375*(alpha[22]+alpha[21])-0.2795084971874737*(alpha[12]+alpha[11])+0.3354101966249685*alpha[4]-0.3354101966249685*alpha[3]+0.25*alpha[0] > 0) { 
    fUpwindQuad[58] = ser_5x_p2_surfx3_quad_58_r(fSkin); 
  } else { 
    fUpwindQuad[58] = ser_5x_p2_surfx3_quad_58_l(fEdge); 
  } 
  if ((-0.5031152949374527*alpha[36])+0.5031152949374527*alpha[33]-0.375*alpha[26]+0.2999999999999999*alpha[25]+0.375*alpha[22]-0.2999999999999999*alpha[21]-0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]+0.45*alpha[8]-0.45*alpha[6]+0.3354101966249685*alpha[4]-0.3354101966249685*alpha[3]+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[59] = ser_5x_p2_surfx3_quad_59_r(fSkin); 
  } else { 
    fUpwindQuad[59] = ser_5x_p2_surfx3_quad_59_l(fEdge); 
  } 
  if ((-0.4024922359499621*alpha[36])+0.4024922359499621*(alpha[35]+alpha[33])-0.4024922359499621*alpha[32]+0.2999999999999999*(alpha[26]+alpha[25])-0.2999999999999999*(alpha[22]+alpha[21])-0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]-0.6037383539249431*alpha[16]+0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[9]-0.45*(alpha[8]+alpha[7])+0.45*alpha[6]-0.45*alpha[5]+0.3354101966249685*alpha[4]-0.3354101966249685*alpha[3]+0.3354101966249685*alpha[2]-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[60] = ser_5x_p2_surfx3_quad_60_r(fSkin); 
  } else { 
    fUpwindQuad[60] = ser_5x_p2_surfx3_quad_60_l(fEdge); 
  } 
  if ((-0.5031152949374527*alpha[35])+0.5031152949374527*alpha[32]+0.2999999999999999*alpha[26]-0.375*alpha[25]-0.2999999999999999*alpha[22]+0.375*alpha[21]-0.375*alpha[19]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]+0.45*alpha[9]-0.45*alpha[7]+0.3354101966249685*alpha[4]-0.3354101966249685*alpha[3]+0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[61] = ser_5x_p2_surfx3_quad_61_r(fSkin); 
  } else { 
    fUpwindQuad[61] = ser_5x_p2_surfx3_quad_61_l(fEdge); 
  } 
  if (0.4024922359499621*(alpha[36]+alpha[35])-0.4024922359499621*(alpha[33]+alpha[32])+0.2999999999999999*(alpha[26]+alpha[25])-0.2999999999999999*(alpha[22]+alpha[21])+0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]+0.6037383539249431*alpha[16]-0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])+0.45*(alpha[9]+alpha[8])-0.45*(alpha[7]+alpha[6])+0.45*alpha[5]+0.3354101966249685*alpha[4]-0.3354101966249685*alpha[3]+0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[62] = ser_5x_p2_surfx3_quad_62_r(fSkin); 
  } else { 
    fUpwindQuad[62] = ser_5x_p2_surfx3_quad_62_l(fEdge); 
  } 
  if ((-0.4024922359499621*(alpha[36]+alpha[35]))+0.2999999999999999*(alpha[26]+alpha[25])-0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]+0.6037383539249431*alpha[16]+0.223606797749979*(alpha[12]+alpha[11])-0.45*(alpha[9]+alpha[8])+0.45*alpha[5]+0.3354101966249685*alpha[4]-0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[63] = ser_5x_p2_surfx3_quad_63_r(fSkin); 
  } else { 
    fUpwindQuad[63] = ser_5x_p2_surfx3_quad_63_l(fEdge); 
  } 
  if (0.5031152949374527*alpha[35]+0.2999999999999999*alpha[26]-0.375*alpha[25]+0.375*alpha[19]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]-0.45*alpha[9]+0.3354101966249685*alpha[4]-0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[64] = ser_5x_p2_surfx3_quad_64_r(fSkin); 
  } else { 
    fUpwindQuad[64] = ser_5x_p2_surfx3_quad_64_l(fEdge); 
  } 
  if (0.4024922359499621*alpha[36]-0.4024922359499621*alpha[35]+0.2999999999999999*(alpha[26]+alpha[25])+0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]-0.6037383539249431*alpha[16]+0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[9]+0.45*alpha[8]-0.45*alpha[5]+0.3354101966249685*alpha[4]-0.3354101966249685*alpha[2]+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[65] = ser_5x_p2_surfx3_quad_65_r(fSkin); 
  } else { 
    fUpwindQuad[65] = ser_5x_p2_surfx3_quad_65_l(fEdge); 
  } 
  if (0.5031152949374527*alpha[36]-0.375*alpha[26]+0.2999999999999999*alpha[25]+0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]-0.45*alpha[8]+0.3354101966249685*alpha[4]-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[66] = ser_5x_p2_surfx3_quad_66_r(fSkin); 
  } else { 
    fUpwindQuad[66] = ser_5x_p2_surfx3_quad_66_l(fEdge); 
  } 
  if ((-0.375*(alpha[26]+alpha[25]))-0.2795084971874737*(alpha[12]+alpha[11])+0.3354101966249685*alpha[4]+0.25*alpha[0] > 0) { 
    fUpwindQuad[67] = ser_5x_p2_surfx3_quad_67_r(fSkin); 
  } else { 
    fUpwindQuad[67] = ser_5x_p2_surfx3_quad_67_l(fEdge); 
  } 
  if ((-0.5031152949374527*alpha[36])-0.375*alpha[26]+0.2999999999999999*alpha[25]-0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]+0.45*alpha[8]+0.3354101966249685*(alpha[4]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[68] = ser_5x_p2_surfx3_quad_68_r(fSkin); 
  } else { 
    fUpwindQuad[68] = ser_5x_p2_surfx3_quad_68_l(fEdge); 
  } 
  if ((-0.4024922359499621*alpha[36])+0.4024922359499621*alpha[35]+0.2999999999999999*(alpha[26]+alpha[25])-0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]-0.6037383539249431*alpha[16]+0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[9]-0.45*(alpha[8]+alpha[5])+0.3354101966249685*(alpha[4]+alpha[2])-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[69] = ser_5x_p2_surfx3_quad_69_r(fSkin); 
  } else { 
    fUpwindQuad[69] = ser_5x_p2_surfx3_quad_69_l(fEdge); 
  } 
  if ((-0.5031152949374527*alpha[35])+0.2999999999999999*alpha[26]-0.375*(alpha[25]+alpha[19])+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]+0.45*alpha[9]+0.3354101966249685*(alpha[4]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad[70] = ser_5x_p2_surfx3_quad_70_r(fSkin); 
  } else { 
    fUpwindQuad[70] = ser_5x_p2_surfx3_quad_70_l(fEdge); 
  } 
  if (0.4024922359499621*(alpha[36]+alpha[35])+0.2999999999999999*(alpha[26]+alpha[25])+0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]+0.6037383539249431*alpha[16]+0.223606797749979*(alpha[12]+alpha[11])+0.45*(alpha[9]+alpha[8]+alpha[5])+0.3354101966249685*(alpha[4]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[71] = ser_5x_p2_surfx3_quad_71_r(fSkin); 
  } else { 
    fUpwindQuad[71] = ser_5x_p2_surfx3_quad_71_l(fEdge); 
  } 
  if ((-0.4024922359499621*(alpha[36]+alpha[35]+alpha[33]+alpha[32]))+0.2999999999999999*(alpha[26]+alpha[25]+alpha[22]+alpha[21])-0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]+0.6037383539249431*(alpha[16]+alpha[15])+0.223606797749979*(alpha[12]+alpha[11])-0.45*(alpha[9]+alpha[8]+alpha[7]+alpha[6])+0.45*alpha[5]+0.3354101966249685*(alpha[4]+alpha[3])-0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[72] = ser_5x_p2_surfx3_quad_72_r(fSkin); 
  } else { 
    fUpwindQuad[72] = ser_5x_p2_surfx3_quad_72_l(fEdge); 
  } 
  if (0.5031152949374527*(alpha[35]+alpha[32])+0.2999999999999999*alpha[26]-0.375*alpha[25]+0.2999999999999999*alpha[22]-0.375*alpha[21]+0.375*alpha[19]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]-0.45*(alpha[9]+alpha[7])+0.3354101966249685*(alpha[4]+alpha[3])-0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[73] = ser_5x_p2_surfx3_quad_73_r(fSkin); 
  } else { 
    fUpwindQuad[73] = ser_5x_p2_surfx3_quad_73_l(fEdge); 
  } 
  if (0.4024922359499621*alpha[36]-0.4024922359499621*alpha[35]+0.4024922359499621*alpha[33]-0.4024922359499621*alpha[32]+0.2999999999999999*(alpha[26]+alpha[25]+alpha[22]+alpha[21])+0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]-0.6037383539249431*(alpha[16]+alpha[15])+0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[9]+0.45*alpha[8]-0.45*alpha[7]+0.45*alpha[6]-0.45*alpha[5]+0.3354101966249685*(alpha[4]+alpha[3])-0.3354101966249685*alpha[2]+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[74] = ser_5x_p2_surfx3_quad_74_r(fSkin); 
  } else { 
    fUpwindQuad[74] = ser_5x_p2_surfx3_quad_74_l(fEdge); 
  } 
  if (0.5031152949374527*(alpha[36]+alpha[33])-0.375*alpha[26]+0.2999999999999999*alpha[25]-0.375*alpha[22]+0.2999999999999999*alpha[21]+0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]-0.45*(alpha[8]+alpha[6])+0.3354101966249685*(alpha[4]+alpha[3])-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[75] = ser_5x_p2_surfx3_quad_75_r(fSkin); 
  } else { 
    fUpwindQuad[75] = ser_5x_p2_surfx3_quad_75_l(fEdge); 
  } 
  if ((-0.375*(alpha[26]+alpha[25]+alpha[22]+alpha[21]))-0.2795084971874737*(alpha[12]+alpha[11])+0.3354101966249685*(alpha[4]+alpha[3])+0.25*alpha[0] > 0) { 
    fUpwindQuad[76] = ser_5x_p2_surfx3_quad_76_r(fSkin); 
  } else { 
    fUpwindQuad[76] = ser_5x_p2_surfx3_quad_76_l(fEdge); 
  } 
  if ((-0.5031152949374527*(alpha[36]+alpha[33]))-0.375*alpha[26]+0.2999999999999999*alpha[25]-0.375*alpha[22]+0.2999999999999999*alpha[21]-0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]+0.45*(alpha[8]+alpha[6])+0.3354101966249685*(alpha[4]+alpha[3]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[77] = ser_5x_p2_surfx3_quad_77_r(fSkin); 
  } else { 
    fUpwindQuad[77] = ser_5x_p2_surfx3_quad_77_l(fEdge); 
  } 
  if ((-0.4024922359499621*alpha[36])+0.4024922359499621*alpha[35]-0.4024922359499621*alpha[33]+0.4024922359499621*alpha[32]+0.2999999999999999*(alpha[26]+alpha[25]+alpha[22]+alpha[21])-0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]-0.6037383539249431*(alpha[16]+alpha[15])+0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[9]-0.45*alpha[8]+0.45*alpha[7]-0.45*(alpha[6]+alpha[5])+0.3354101966249685*(alpha[4]+alpha[3]+alpha[2])-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[78] = ser_5x_p2_surfx3_quad_78_r(fSkin); 
  } else { 
    fUpwindQuad[78] = ser_5x_p2_surfx3_quad_78_l(fEdge); 
  } 
  if ((-0.5031152949374527*(alpha[35]+alpha[32]))+0.2999999999999999*alpha[26]-0.375*alpha[25]+0.2999999999999999*alpha[22]-0.375*(alpha[21]+alpha[19])+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]+0.45*(alpha[9]+alpha[7])+0.3354101966249685*(alpha[4]+alpha[3]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad[79] = ser_5x_p2_surfx3_quad_79_r(fSkin); 
  } else { 
    fUpwindQuad[79] = ser_5x_p2_surfx3_quad_79_l(fEdge); 
  } 
  if (0.4024922359499621*(alpha[36]+alpha[35]+alpha[33]+alpha[32])+0.2999999999999999*(alpha[26]+alpha[25]+alpha[22]+alpha[21])+0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]+0.6037383539249431*(alpha[16]+alpha[15])+0.223606797749979*(alpha[12]+alpha[11])+0.45*(alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5])+0.3354101966249685*(alpha[4]+alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[80] = ser_5x_p2_surfx3_quad_80_r(fSkin); 
  } else { 
    fUpwindQuad[80] = ser_5x_p2_surfx3_quad_80_l(fEdge); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_5x_p2_upwind(fUpwindQuad, fUpwind); 

  Ghat[0] += 0.25*(alpha[36]*fUpwind[36]+alpha[35]*fUpwind[35]+alpha[33]*fUpwind[33]+alpha[32]*fUpwind[32]+alpha[26]*fUpwind[26]+alpha[25]*fUpwind[25]+alpha[22]*fUpwind[22]+alpha[21]*fUpwind[21]+alpha[20]*fUpwind[20]+alpha[19]*fUpwind[19]+alpha[16]*fUpwind[16]+alpha[15]*fUpwind[15]+alpha[12]*fUpwind[12]+alpha[11]*fUpwind[11]+alpha[9]*fUpwind[9]+alpha[8]*fUpwind[8]+alpha[7]*fUpwind[7]+alpha[6]*fUpwind[6]+alpha[5]*fUpwind[5]+alpha[4]*fUpwind[4]+alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] += 0.2500000000000001*(alpha[26]*fUpwind[36]+fUpwind[26]*alpha[36])+0.223606797749979*(alpha[16]*fUpwind[35]+fUpwind[16]*alpha[35])+0.2500000000000001*(alpha[22]*fUpwind[33]+fUpwind[22]*alpha[33])+0.223606797749979*(alpha[15]*fUpwind[32]+fUpwind[15]*alpha[32])+0.223606797749979*(alpha[8]*fUpwind[25]+fUpwind[8]*alpha[25]+alpha[6]*fUpwind[21]+fUpwind[6]*alpha[21])+0.2500000000000001*(alpha[12]*fUpwind[20]+fUpwind[12]*alpha[20])+0.223606797749979*(alpha[5]*fUpwind[19]+fUpwind[5]*alpha[19])+0.25*(alpha[9]*fUpwind[16]+fUpwind[9]*alpha[16]+alpha[7]*fUpwind[15]+fUpwind[7]*alpha[15])+0.223606797749979*(alpha[1]*fUpwind[11]+fUpwind[1]*alpha[11])+0.25*(alpha[4]*fUpwind[8]+fUpwind[4]*alpha[8]+alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6]+alpha[2]*fUpwind[5]+fUpwind[2]*alpha[5]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] += 0.223606797749979*(alpha[16]*fUpwind[36]+fUpwind[16]*alpha[36])+0.2500000000000001*(alpha[25]*fUpwind[35]+fUpwind[25]*alpha[35])+0.223606797749979*(alpha[15]*fUpwind[33]+fUpwind[15]*alpha[33])+0.2500000000000001*(alpha[21]*fUpwind[32]+fUpwind[21]*alpha[32])+0.223606797749979*(alpha[9]*fUpwind[26]+fUpwind[9]*alpha[26]+alpha[7]*fUpwind[22]+fUpwind[7]*alpha[22]+alpha[5]*fUpwind[20]+fUpwind[5]*alpha[20])+0.2500000000000001*(alpha[11]*fUpwind[19]+fUpwind[11]*alpha[19])+0.25*(alpha[8]*fUpwind[16]+fUpwind[8]*alpha[16]+alpha[6]*fUpwind[15]+fUpwind[6]*alpha[15])+0.223606797749979*(alpha[2]*fUpwind[12]+fUpwind[2]*alpha[12])+0.25*(alpha[4]*fUpwind[9]+fUpwind[4]*alpha[9]+alpha[3]*fUpwind[7]+fUpwind[3]*alpha[7]+alpha[1]*fUpwind[5]+fUpwind[1]*alpha[5]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] += 0.2500000000000001*(alpha[36]*fUpwind[45]+alpha[35]*fUpwind[44]+alpha[26]*fUpwind[38]+alpha[25]*fUpwind[37])+0.223606797749979*alpha[15]*fUpwind[34]+0.2500000000000001*(alpha[20]*fUpwind[33]+fUpwind[20]*alpha[33]+alpha[19]*fUpwind[32]+fUpwind[19]*alpha[32])+0.25*alpha[16]*fUpwind[31]+0.223606797749979*(alpha[7]*fUpwind[24]+alpha[6]*fUpwind[23])+0.2500000000000001*(alpha[12]*fUpwind[22]+fUpwind[12]*alpha[22]+alpha[11]*fUpwind[21]+fUpwind[11]*alpha[21])+0.25*(alpha[9]*fUpwind[18]+alpha[8]*fUpwind[17]+alpha[5]*fUpwind[15]+fUpwind[5]*alpha[15])+0.223606797749979*alpha[3]*fUpwind[13]+0.25*(alpha[4]*fUpwind[10]+alpha[2]*fUpwind[7]+fUpwind[2]*alpha[7]+alpha[1]*fUpwind[6]+fUpwind[1]*alpha[6]+alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]); 
  Ghat[4] += 0.2500000000000001*(alpha[33]*fUpwind[45]+alpha[32]*fUpwind[44])+0.223606797749979*alpha[16]*fUpwind[41]+0.2500000000000001*(alpha[22]*fUpwind[38]+alpha[21]*fUpwind[37]+alpha[20]*fUpwind[36]+fUpwind[20]*alpha[36]+alpha[19]*fUpwind[35]+fUpwind[19]*alpha[35])+0.25*alpha[15]*fUpwind[31]+0.223606797749979*(alpha[9]*fUpwind[29]+alpha[8]*fUpwind[28])+0.2500000000000001*(alpha[12]*fUpwind[26]+fUpwind[12]*alpha[26]+alpha[11]*fUpwind[25]+fUpwind[11]*alpha[25])+0.25*(alpha[7]*fUpwind[18]+alpha[6]*fUpwind[17]+alpha[5]*fUpwind[16]+fUpwind[5]*alpha[16])+0.223606797749979*alpha[4]*fUpwind[14]+0.25*(alpha[3]*fUpwind[10]+alpha[2]*fUpwind[9]+fUpwind[2]*alpha[9]+alpha[1]*fUpwind[8]+fUpwind[1]*alpha[8]+alpha[0]*fUpwind[4]+fUpwind[0]*alpha[4]); 
  Ghat[5] += (0.2*alpha[35]+0.223606797749979*alpha[9])*fUpwind[36]+0.2*fUpwind[35]*alpha[36]+0.223606797749979*(fUpwind[9]*alpha[36]+alpha[8]*fUpwind[35]+fUpwind[8]*alpha[35])+(0.2*alpha[32]+0.223606797749979*alpha[7])*fUpwind[33]+0.2*fUpwind[32]*alpha[33]+0.223606797749979*(fUpwind[7]*alpha[33]+alpha[6]*fUpwind[32]+fUpwind[6]*alpha[32])+0.223606797749979*(alpha[16]*fUpwind[26]+fUpwind[16]*alpha[26]+alpha[16]*fUpwind[25]+fUpwind[16]*alpha[25]+alpha[15]*fUpwind[22]+fUpwind[15]*alpha[22]+alpha[15]*fUpwind[21]+fUpwind[15]*alpha[21])+(0.2*alpha[19]+0.223606797749979*alpha[2])*fUpwind[20]+0.2*fUpwind[19]*alpha[20]+0.223606797749979*(fUpwind[2]*alpha[20]+alpha[1]*fUpwind[19]+fUpwind[1]*alpha[19])+0.25*(alpha[4]*fUpwind[16]+fUpwind[4]*alpha[16]+alpha[3]*fUpwind[15]+fUpwind[3]*alpha[15])+0.223606797749979*(alpha[5]*fUpwind[12]+fUpwind[5]*alpha[12]+alpha[5]*fUpwind[11]+fUpwind[5]*alpha[11])+0.25*(alpha[8]*fUpwind[9]+fUpwind[8]*alpha[9]+alpha[6]*fUpwind[7]+fUpwind[6]*alpha[7]+alpha[0]*fUpwind[5]+fUpwind[0]*alpha[5]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[6] += 0.25*alpha[26]*fUpwind[45]+0.223606797749979*alpha[16]*fUpwind[44]+0.25*alpha[36]*fUpwind[38]+0.223606797749979*(alpha[8]*fUpwind[37]+fUpwind[31]*alpha[35])+(0.2*alpha[32]+0.223606797749979*alpha[7])*fUpwind[34]+0.25*(alpha[12]*fUpwind[33]+fUpwind[12]*alpha[33])+0.223606797749979*(alpha[5]*fUpwind[32]+fUpwind[5]*alpha[32])+0.25*alpha[9]*fUpwind[31]+0.223606797749979*(fUpwind[17]*alpha[25]+alpha[15]*fUpwind[24])+(0.2*alpha[21]+0.223606797749979*alpha[3])*fUpwind[23]+0.25*(alpha[20]*fUpwind[22]+fUpwind[20]*alpha[22])+0.223606797749979*(alpha[1]*fUpwind[21]+fUpwind[1]*alpha[21]+alpha[15]*fUpwind[19]+fUpwind[15]*alpha[19])+0.25*(alpha[16]*fUpwind[18]+alpha[4]*fUpwind[17]+alpha[2]*fUpwind[15]+fUpwind[2]*alpha[15])+0.223606797749979*(alpha[6]*(fUpwind[13]+fUpwind[11])+fUpwind[6]*alpha[11])+0.25*(alpha[8]*fUpwind[10]+alpha[5]*fUpwind[7]+fUpwind[5]*alpha[7]+alpha[0]*fUpwind[6]+fUpwind[0]*alpha[6]+alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]); 
  Ghat[7] += 0.223606797749979*alpha[16]*fUpwind[45]+0.25*alpha[25]*fUpwind[44]+0.223606797749979*alpha[9]*fUpwind[38]+0.25*alpha[35]*fUpwind[37]+0.223606797749979*fUpwind[31]*alpha[36]+0.2*alpha[33]*fUpwind[34]+0.223606797749979*(alpha[6]*fUpwind[34]+alpha[5]*fUpwind[33]+fUpwind[5]*alpha[33])+0.25*(alpha[11]*fUpwind[32]+fUpwind[11]*alpha[32]+alpha[8]*fUpwind[31])+0.223606797749979*fUpwind[18]*alpha[26]+0.2*alpha[22]*fUpwind[24]+0.223606797749979*(alpha[3]*fUpwind[24]+alpha[15]*fUpwind[23]+alpha[2]*fUpwind[22]+fUpwind[2]*alpha[22])+0.25*(alpha[19]*fUpwind[21]+fUpwind[19]*alpha[21])+0.223606797749979*(alpha[15]*fUpwind[20]+fUpwind[15]*alpha[20])+0.25*(alpha[4]*fUpwind[18]+alpha[16]*fUpwind[17]+alpha[1]*fUpwind[15]+fUpwind[1]*alpha[15])+0.223606797749979*(alpha[7]*(fUpwind[13]+fUpwind[12])+fUpwind[7]*alpha[12])+0.25*(alpha[9]*fUpwind[10]+alpha[0]*fUpwind[7]+fUpwind[0]*alpha[7]+alpha[5]*fUpwind[6]+fUpwind[5]*alpha[6]+alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]); 
  Ghat[8] += 0.25*alpha[22]*fUpwind[45]+0.223606797749979*alpha[15]*fUpwind[44]+(0.2*alpha[35]+0.223606797749979*alpha[9])*fUpwind[41]+0.25*alpha[33]*fUpwind[38]+0.223606797749979*alpha[6]*fUpwind[37]+0.25*(alpha[12]*fUpwind[36]+fUpwind[12]*alpha[36])+0.223606797749979*(alpha[5]*fUpwind[35]+fUpwind[5]*alpha[35])+fUpwind[31]*(0.223606797749979*alpha[32]+0.25*alpha[7])+0.223606797749979*alpha[16]*fUpwind[29]+(0.2*alpha[25]+0.223606797749979*alpha[4])*fUpwind[28]+0.25*(alpha[20]*fUpwind[26]+fUpwind[20]*alpha[26])+0.223606797749979*(alpha[1]*fUpwind[25]+fUpwind[1]*alpha[25]+fUpwind[17]*alpha[21]+alpha[16]*fUpwind[19]+fUpwind[16]*alpha[19])+0.25*(alpha[15]*fUpwind[18]+alpha[3]*fUpwind[17]+alpha[2]*fUpwind[16]+fUpwind[2]*alpha[16])+0.223606797749979*(alpha[8]*(fUpwind[14]+fUpwind[11])+fUpwind[8]*alpha[11])+0.25*(alpha[6]*fUpwind[10]+alpha[5]*fUpwind[9]+fUpwind[5]*alpha[9]+alpha[0]*fUpwind[8]+fUpwind[0]*alpha[8]+alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]); 
  Ghat[9] += 0.223606797749979*alpha[15]*fUpwind[45]+0.25*alpha[21]*fUpwind[44]+0.2*alpha[36]*fUpwind[41]+0.223606797749979*(alpha[8]*fUpwind[41]+alpha[7]*fUpwind[38])+0.25*alpha[32]*fUpwind[37]+0.223606797749979*(alpha[5]*fUpwind[36]+fUpwind[5]*alpha[36])+0.25*(alpha[11]*fUpwind[35]+fUpwind[11]*alpha[35])+fUpwind[31]*(0.223606797749979*alpha[33]+0.25*alpha[6])+0.2*alpha[26]*fUpwind[29]+0.223606797749979*(alpha[4]*fUpwind[29]+alpha[16]*fUpwind[28]+alpha[2]*fUpwind[26]+fUpwind[2]*alpha[26])+0.25*(alpha[19]*fUpwind[25]+fUpwind[19]*alpha[25])+0.223606797749979*(fUpwind[18]*alpha[22]+alpha[16]*fUpwind[20]+fUpwind[16]*alpha[20])+0.25*(alpha[3]*fUpwind[18]+alpha[15]*fUpwind[17]+alpha[1]*fUpwind[16]+fUpwind[1]*alpha[16])+0.223606797749979*(alpha[9]*(fUpwind[14]+fUpwind[12])+fUpwind[9]*alpha[12])+0.25*(alpha[7]*fUpwind[10]+alpha[0]*fUpwind[9]+fUpwind[0]*alpha[9]+alpha[5]*fUpwind[8]+fUpwind[5]*alpha[8]+alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]); 
  Ghat[10] += 0.223606797749979*(alpha[16]*fUpwind[47]+alpha[15]*fUpwind[46])+0.25*(alpha[20]*fUpwind[45]+alpha[19]*fUpwind[44])+0.223606797749979*(alpha[9]*fUpwind[43]+alpha[8]*fUpwind[42]+alpha[7]*fUpwind[40]+alpha[6]*fUpwind[39])+0.25*(alpha[12]*fUpwind[38]+alpha[11]*fUpwind[37]+alpha[33]*fUpwind[36]+fUpwind[33]*alpha[36]+alpha[32]*fUpwind[35]+fUpwind[32]*alpha[35]+alpha[5]*fUpwind[31])+0.223606797749979*(alpha[4]*fUpwind[30]+alpha[3]*fUpwind[27])+0.25*(alpha[22]*fUpwind[26]+fUpwind[22]*alpha[26]+alpha[21]*fUpwind[25]+fUpwind[21]*alpha[25]+alpha[2]*fUpwind[18]+alpha[1]*fUpwind[17]+alpha[15]*fUpwind[16]+fUpwind[15]*alpha[16]+alpha[0]*fUpwind[10]+alpha[7]*fUpwind[9]+fUpwind[7]*alpha[9]+alpha[6]*fUpwind[8]+fUpwind[6]*alpha[8]+alpha[3]*fUpwind[4]+fUpwind[3]*alpha[4]); 
  Ghat[11] += 0.223606797749979*alpha[36]*fUpwind[36]+0.159719141249985*alpha[35]*fUpwind[35]+0.25*(alpha[9]*fUpwind[35]+fUpwind[9]*alpha[35])+0.223606797749979*alpha[33]*fUpwind[33]+0.159719141249985*alpha[32]*fUpwind[32]+0.25*(alpha[7]*fUpwind[32]+fUpwind[7]*alpha[32])+0.159719141249985*alpha[25]*fUpwind[25]+0.2500000000000001*(alpha[4]*fUpwind[25]+fUpwind[4]*alpha[25])+0.159719141249985*alpha[21]*fUpwind[21]+0.2500000000000001*(alpha[3]*fUpwind[21]+fUpwind[3]*alpha[21])+0.223606797749979*alpha[20]*fUpwind[20]+0.159719141249985*alpha[19]*fUpwind[19]+0.2500000000000001*(alpha[2]*fUpwind[19]+fUpwind[2]*alpha[19])+0.223606797749979*(alpha[16]*fUpwind[16]+alpha[15]*fUpwind[15])+0.159719141249985*alpha[11]*fUpwind[11]+0.25*(alpha[0]*fUpwind[11]+fUpwind[0]*alpha[11])+0.223606797749979*(alpha[8]*fUpwind[8]+alpha[6]*fUpwind[6]+alpha[5]*fUpwind[5]+alpha[1]*fUpwind[1]); 
  Ghat[12] += 0.159719141249985*alpha[36]*fUpwind[36]+0.25*(alpha[8]*fUpwind[36]+fUpwind[8]*alpha[36])+0.223606797749979*alpha[35]*fUpwind[35]+0.159719141249985*alpha[33]*fUpwind[33]+0.25*(alpha[6]*fUpwind[33]+fUpwind[6]*alpha[33])+0.223606797749979*alpha[32]*fUpwind[32]+0.159719141249985*alpha[26]*fUpwind[26]+0.2500000000000001*(alpha[4]*fUpwind[26]+fUpwind[4]*alpha[26])+0.159719141249985*alpha[22]*fUpwind[22]+0.2500000000000001*(alpha[3]*fUpwind[22]+fUpwind[3]*alpha[22])+0.159719141249985*alpha[20]*fUpwind[20]+0.2500000000000001*(alpha[1]*fUpwind[20]+fUpwind[1]*alpha[20])+0.223606797749979*(alpha[19]*fUpwind[19]+alpha[16]*fUpwind[16]+alpha[15]*fUpwind[15])+0.159719141249985*alpha[12]*fUpwind[12]+0.25*(alpha[0]*fUpwind[12]+fUpwind[0]*alpha[12])+0.223606797749979*(alpha[9]*fUpwind[9]+alpha[7]*fUpwind[7]+alpha[5]*fUpwind[5]+alpha[2]*fUpwind[2]); 
  Ghat[13] += 0.2500000000000001*alpha[16]*fUpwind[46]+0.25*(alpha[9]*fUpwind[40]+alpha[8]*fUpwind[39]+alpha[5]*fUpwind[34])+0.223606797749979*(alpha[33]*fUpwind[33]+alpha[32]*fUpwind[32])+0.2500000000000001*(alpha[4]*fUpwind[27]+alpha[2]*fUpwind[24]+alpha[1]*fUpwind[23])+0.223606797749979*(alpha[22]*fUpwind[22]+alpha[21]*fUpwind[21]+alpha[15]*fUpwind[15])+0.25*alpha[0]*fUpwind[13]+0.223606797749979*(alpha[7]*fUpwind[7]+alpha[6]*fUpwind[6]+alpha[3]*fUpwind[3]); 
  Ghat[14] += 0.2500000000000001*alpha[15]*fUpwind[47]+0.25*(alpha[7]*fUpwind[43]+alpha[6]*fUpwind[42]+alpha[5]*fUpwind[41])+0.223606797749979*(alpha[36]*fUpwind[36]+alpha[35]*fUpwind[35])+0.2500000000000001*(alpha[3]*fUpwind[30]+alpha[2]*fUpwind[29]+alpha[1]*fUpwind[28])+0.223606797749979*(alpha[26]*fUpwind[26]+alpha[25]*fUpwind[25]+alpha[16]*fUpwind[16])+0.25*alpha[0]*fUpwind[14]+0.223606797749979*(alpha[9]*fUpwind[9]+alpha[8]*fUpwind[8]+alpha[4]*fUpwind[4]); 
  Ghat[15] += (0.2*alpha[35]+0.223606797749979*alpha[9])*fUpwind[45]+(0.2*alpha[36]+0.223606797749979*alpha[8])*fUpwind[44]+0.223606797749979*(alpha[16]*(fUpwind[38]+fUpwind[37])+fUpwind[18]*alpha[36]+fUpwind[17]*alpha[35])+(0.2*(alpha[22]+alpha[21])+0.223606797749979*alpha[3])*fUpwind[34]+(0.2*alpha[19]+0.223606797749979*alpha[2])*fUpwind[33]+(0.2*(fUpwind[24]+fUpwind[19])+0.223606797749979*fUpwind[2])*alpha[33]+(0.2*alpha[20]+0.223606797749979*alpha[1])*fUpwind[32]+(0.2*(fUpwind[23]+fUpwind[20])+0.223606797749979*fUpwind[1])*alpha[32]+(0.223606797749979*(alpha[26]+alpha[25])+0.25*alpha[4])*fUpwind[31]+0.223606797749979*(alpha[6]*fUpwind[24]+alpha[7]*fUpwind[23]+alpha[5]*fUpwind[22]+fUpwind[5]*alpha[22]+alpha[5]*fUpwind[21]+fUpwind[5]*alpha[21]+alpha[7]*fUpwind[20]+fUpwind[7]*alpha[20]+alpha[6]*fUpwind[19]+fUpwind[6]*alpha[19])+0.25*(alpha[8]*fUpwind[18]+alpha[9]*fUpwind[17]+fUpwind[10]*alpha[16])+(0.223606797749979*(alpha[12]+alpha[11])+0.25*alpha[0])*fUpwind[15]+0.223606797749979*(fUpwind[13]+fUpwind[12]+fUpwind[11])*alpha[15]+0.25*(fUpwind[0]*alpha[15]+alpha[1]*fUpwind[7]+fUpwind[1]*alpha[7]+alpha[2]*fUpwind[6]+fUpwind[2]*alpha[6]+alpha[3]*fUpwind[5]+fUpwind[3]*alpha[5]); 
  Ghat[16] += (0.2*alpha[32]+0.223606797749979*alpha[7])*fUpwind[45]+(0.2*alpha[33]+0.223606797749979*alpha[6])*fUpwind[44]+0.2*(alpha[26]+alpha[25])*fUpwind[41]+0.223606797749979*(alpha[4]*fUpwind[41]+alpha[15]*(fUpwind[38]+fUpwind[37]))+(0.2*alpha[19]+0.223606797749979*alpha[2])*fUpwind[36]+(0.2*(fUpwind[29]+fUpwind[19])+0.223606797749979*fUpwind[2])*alpha[36]+(0.2*alpha[20]+0.223606797749979*alpha[1])*fUpwind[35]+0.2*(fUpwind[28]+fUpwind[20])*alpha[35]+0.223606797749979*(fUpwind[1]*alpha[35]+fUpwind[18]*alpha[33]+fUpwind[17]*alpha[32])+(0.223606797749979*(alpha[22]+alpha[21])+0.25*alpha[3])*fUpwind[31]+0.223606797749979*(alpha[8]*fUpwind[29]+alpha[9]*fUpwind[28]+alpha[5]*fUpwind[26]+fUpwind[5]*alpha[26]+alpha[5]*fUpwind[25]+fUpwind[5]*alpha[25]+alpha[9]*fUpwind[20]+fUpwind[9]*alpha[20]+alpha[8]*fUpwind[19]+fUpwind[8]*alpha[19])+0.25*(alpha[6]*fUpwind[18]+alpha[7]*fUpwind[17])+(0.223606797749979*(alpha[12]+alpha[11])+0.25*alpha[0])*fUpwind[16]+0.223606797749979*(fUpwind[14]+fUpwind[12]+fUpwind[11])*alpha[16]+0.25*(fUpwind[0]*alpha[16]+fUpwind[10]*alpha[15]+alpha[1]*fUpwind[9]+fUpwind[1]*alpha[9]+alpha[2]*fUpwind[8]+fUpwind[2]*alpha[8]+alpha[4]*fUpwind[5]+fUpwind[4]*alpha[5]); 
  Ghat[17] += (0.2*alpha[35]+0.223606797749979*alpha[9])*fUpwind[47]+(0.2*alpha[32]+0.223606797749979*alpha[7])*fUpwind[46]+0.2500000000000001*alpha[12]*fUpwind[45]+0.223606797749979*alpha[5]*fUpwind[44]+0.223606797749979*alpha[16]*fUpwind[43]+0.2*alpha[25]*fUpwind[42]+0.223606797749979*(alpha[4]*fUpwind[42]+alpha[15]*fUpwind[40])+(0.2*alpha[21]+0.223606797749979*alpha[3])*fUpwind[39]+0.2500000000000001*alpha[20]*fUpwind[38]+0.223606797749979*alpha[1]*fUpwind[37]+0.2500000000000001*(alpha[22]*fUpwind[36]+fUpwind[22]*alpha[36])+0.223606797749979*(alpha[15]*fUpwind[35]+fUpwind[15]*alpha[35])+0.2500000000000001*(alpha[26]*fUpwind[33]+fUpwind[26]*alpha[33])+0.223606797749979*(alpha[16]*fUpwind[32]+fUpwind[16]*alpha[32])+(0.223606797749979*alpha[19]+0.25*alpha[2])*fUpwind[31]+0.223606797749979*(alpha[8]*fUpwind[30]+alpha[6]*(fUpwind[27]+fUpwind[25])+fUpwind[6]*alpha[25]+alpha[8]*fUpwind[21]+fUpwind[8]*alpha[21])+0.25*alpha[5]*fUpwind[18]+0.223606797749979*alpha[11]*fUpwind[17]+0.25*(alpha[0]*fUpwind[17]+alpha[7]*fUpwind[16]+fUpwind[7]*alpha[16]+alpha[9]*fUpwind[15]+fUpwind[9]*alpha[15]+alpha[1]*fUpwind[10]+alpha[3]*fUpwind[8]+fUpwind[3]*alpha[8]+alpha[4]*fUpwind[6]+fUpwind[4]*alpha[6]); 
  Ghat[18] += (0.2*alpha[36]+0.223606797749979*alpha[8])*fUpwind[47]+0.2*alpha[33]*fUpwind[46]+0.223606797749979*(alpha[6]*fUpwind[46]+alpha[5]*fUpwind[45])+0.2500000000000001*alpha[11]*fUpwind[44]+0.2*alpha[26]*fUpwind[43]+0.223606797749979*(alpha[4]*fUpwind[43]+alpha[16]*fUpwind[42])+0.2*alpha[22]*fUpwind[40]+0.223606797749979*(alpha[3]*fUpwind[40]+alpha[15]*fUpwind[39]+alpha[2]*fUpwind[38])+0.2500000000000001*alpha[19]*fUpwind[37]+0.223606797749979*(alpha[15]*fUpwind[36]+fUpwind[15]*alpha[36])+0.2500000000000001*(alpha[21]*fUpwind[35]+fUpwind[21]*alpha[35])+0.223606797749979*(alpha[16]*fUpwind[33]+fUpwind[16]*alpha[33])+0.2500000000000001*(alpha[25]*fUpwind[32]+fUpwind[25]*alpha[32])+(0.223606797749979*alpha[20]+0.25*alpha[1])*fUpwind[31]+0.223606797749979*(alpha[9]*fUpwind[30]+alpha[7]*(fUpwind[27]+fUpwind[26])+fUpwind[7]*alpha[26]+alpha[9]*fUpwind[22]+fUpwind[9]*alpha[22])+0.223606797749979*alpha[12]*fUpwind[18]+0.25*(alpha[0]*fUpwind[18]+alpha[5]*fUpwind[17]+alpha[6]*fUpwind[16]+fUpwind[6]*alpha[16]+alpha[8]*fUpwind[15]+fUpwind[8]*alpha[15]+alpha[2]*fUpwind[10]+alpha[3]*fUpwind[9]+fUpwind[3]*alpha[9]+alpha[4]*fUpwind[7]+fUpwind[4]*alpha[7]); 
  Ghat[19] += 0.2*(alpha[16]*fUpwind[36]+fUpwind[16]*alpha[36])+(0.223606797749979*alpha[26]+0.159719141249985*alpha[25]+0.2500000000000001*alpha[4])*fUpwind[35]+(0.223606797749979*fUpwind[26]+0.159719141249985*fUpwind[25]+0.2500000000000001*fUpwind[4])*alpha[35]+0.2*(alpha[15]*fUpwind[33]+fUpwind[15]*alpha[33])+(0.223606797749979*alpha[22]+0.159719141249985*alpha[21]+0.2500000000000001*alpha[3])*fUpwind[32]+(0.223606797749979*fUpwind[22]+0.159719141249985*fUpwind[21]+0.2500000000000001*fUpwind[3])*alpha[32]+0.25*(alpha[9]*fUpwind[25]+fUpwind[9]*alpha[25]+alpha[7]*fUpwind[21]+fUpwind[7]*alpha[21])+0.2*(alpha[5]*fUpwind[20]+fUpwind[5]*alpha[20])+(0.223606797749979*alpha[12]+0.159719141249985*alpha[11]+0.25*alpha[0])*fUpwind[19]+(0.223606797749979*fUpwind[12]+0.159719141249985*fUpwind[11]+0.25*fUpwind[0])*alpha[19]+0.223606797749979*(alpha[8]*fUpwind[16]+fUpwind[8]*alpha[16]+alpha[6]*fUpwind[15]+fUpwind[6]*alpha[15])+0.2500000000000001*(alpha[2]*fUpwind[11]+fUpwind[2]*alpha[11])+0.223606797749979*(alpha[1]*fUpwind[5]+fUpwind[1]*alpha[5]); 
  Ghat[20] += (0.159719141249985*alpha[26]+0.223606797749979*alpha[25]+0.2500000000000001*alpha[4])*fUpwind[36]+(0.159719141249985*fUpwind[26]+0.223606797749979*fUpwind[25]+0.2500000000000001*fUpwind[4])*alpha[36]+0.2*(alpha[16]*fUpwind[35]+fUpwind[16]*alpha[35])+(0.159719141249985*alpha[22]+0.223606797749979*alpha[21]+0.2500000000000001*alpha[3])*fUpwind[33]+(0.159719141249985*fUpwind[22]+0.223606797749979*fUpwind[21]+0.2500000000000001*fUpwind[3])*alpha[33]+0.2*(alpha[15]*fUpwind[32]+fUpwind[15]*alpha[32])+0.25*(alpha[8]*fUpwind[26]+fUpwind[8]*alpha[26]+alpha[6]*fUpwind[22]+fUpwind[6]*alpha[22])+(0.159719141249985*alpha[12]+0.223606797749979*alpha[11]+0.25*alpha[0])*fUpwind[20]+(0.159719141249985*fUpwind[12]+0.223606797749979*fUpwind[11]+0.25*fUpwind[0])*alpha[20]+0.2*(alpha[5]*fUpwind[19]+fUpwind[5]*alpha[19])+0.223606797749979*(alpha[9]*fUpwind[16]+fUpwind[9]*alpha[16]+alpha[7]*fUpwind[15]+fUpwind[7]*alpha[15])+0.2500000000000001*(alpha[1]*fUpwind[12]+fUpwind[1]*alpha[12])+0.223606797749979*(alpha[2]*fUpwind[5]+fUpwind[2]*alpha[5]); 
  Ghat[21] += 0.223606797749979*alpha[36]*fUpwind[45]+(0.159719141249985*alpha[35]+0.25*alpha[9])*fUpwind[44]+0.159719141249985*alpha[25]*fUpwind[37]+0.2500000000000001*(alpha[4]*fUpwind[37]+fUpwind[18]*alpha[35])+0.2*alpha[15]*fUpwind[34]+0.223606797749979*(alpha[20]*fUpwind[33]+fUpwind[20]*alpha[33])+(0.159719141249985*alpha[19]+0.2500000000000001*alpha[2])*fUpwind[32]+(0.223606797749979*fUpwind[24]+0.159719141249985*fUpwind[19]+0.2500000000000001*fUpwind[2])*alpha[32]+0.223606797749979*alpha[16]*fUpwind[31]+0.25*fUpwind[10]*alpha[25]+0.2*alpha[6]*fUpwind[23]+(0.159719141249985*alpha[11]+0.25*alpha[0])*fUpwind[21]+(0.223606797749979*fUpwind[13]+0.159719141249985*fUpwind[11])*alpha[21]+0.25*(fUpwind[0]*alpha[21]+alpha[7]*fUpwind[19]+fUpwind[7]*alpha[19])+0.223606797749979*(alpha[8]*fUpwind[17]+alpha[5]*fUpwind[15]+fUpwind[5]*alpha[15])+0.2500000000000001*(alpha[3]*fUpwind[11]+fUpwind[3]*alpha[11])+0.223606797749979*(alpha[1]*fUpwind[6]+fUpwind[1]*alpha[6]); 
  Ghat[22] += (0.159719141249985*alpha[36]+0.25*alpha[8])*fUpwind[45]+0.223606797749979*alpha[35]*fUpwind[44]+0.159719141249985*alpha[26]*fUpwind[38]+0.2500000000000001*(alpha[4]*fUpwind[38]+fUpwind[17]*alpha[36])+0.2*alpha[15]*fUpwind[34]+(0.159719141249985*alpha[20]+0.2500000000000001*alpha[1])*fUpwind[33]+(0.223606797749979*fUpwind[23]+0.159719141249985*fUpwind[20]+0.2500000000000001*fUpwind[1])*alpha[33]+0.223606797749979*(alpha[19]*fUpwind[32]+fUpwind[19]*alpha[32])+0.223606797749979*alpha[16]*fUpwind[31]+0.25*fUpwind[10]*alpha[26]+0.2*alpha[7]*fUpwind[24]+(0.159719141249985*alpha[12]+0.25*alpha[0])*fUpwind[22]+(0.223606797749979*fUpwind[13]+0.159719141249985*fUpwind[12])*alpha[22]+0.25*(fUpwind[0]*alpha[22]+alpha[6]*fUpwind[20]+fUpwind[6]*alpha[20])+0.223606797749979*(alpha[9]*fUpwind[18]+alpha[5]*fUpwind[15]+fUpwind[5]*alpha[15])+0.2500000000000001*(alpha[3]*fUpwind[12]+fUpwind[3]*alpha[12])+0.223606797749979*(alpha[2]*fUpwind[7]+fUpwind[2]*alpha[7]); 
  Ghat[23] += (0.223606797749979*alpha[35]+0.25*alpha[9])*fUpwind[46]+0.2500000000000001*alpha[16]*fUpwind[40]+(0.223606797749979*alpha[25]+0.2500000000000001*alpha[4])*fUpwind[39]+(0.223606797749979*alpha[19]+0.2500000000000001*alpha[2])*fUpwind[34]+0.223606797749979*(alpha[22]*fUpwind[33]+fUpwind[22]*alpha[33])+0.2*(alpha[15]*fUpwind[32]+fUpwind[15]*alpha[32])+0.25*(alpha[8]*fUpwind[27]+alpha[5]*fUpwind[24])+(0.223606797749979*alpha[11]+0.25*alpha[0])*fUpwind[23]+0.2*(alpha[6]*fUpwind[21]+fUpwind[6]*alpha[21])+0.223606797749979*(alpha[7]*fUpwind[15]+fUpwind[7]*alpha[15])+0.2500000000000001*alpha[1]*fUpwind[13]+0.223606797749979*(alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6]); 
  Ghat[24] += (0.223606797749979*alpha[36]+0.25*alpha[8])*fUpwind[46]+0.223606797749979*alpha[26]*fUpwind[40]+0.2500000000000001*(alpha[4]*fUpwind[40]+alpha[16]*fUpwind[39])+(0.223606797749979*alpha[20]+0.2500000000000001*alpha[1])*fUpwind[34]+0.2*(alpha[15]*fUpwind[33]+fUpwind[15]*alpha[33])+0.223606797749979*(alpha[21]*fUpwind[32]+fUpwind[21]*alpha[32])+0.25*alpha[9]*fUpwind[27]+0.223606797749979*alpha[12]*fUpwind[24]+0.25*(alpha[0]*fUpwind[24]+alpha[5]*fUpwind[23])+0.2*(alpha[7]*fUpwind[22]+fUpwind[7]*alpha[22])+0.223606797749979*(alpha[6]*fUpwind[15]+fUpwind[6]*alpha[15])+0.2500000000000001*alpha[2]*fUpwind[13]+0.223606797749979*(alpha[3]*fUpwind[7]+fUpwind[3]*alpha[7]); 
  Ghat[25] += 0.223606797749979*alpha[33]*fUpwind[45]+(0.159719141249985*alpha[32]+0.25*alpha[7])*fUpwind[44]+0.2*alpha[16]*fUpwind[41]+(0.159719141249985*alpha[21]+0.2500000000000001*alpha[3])*fUpwind[37]+0.223606797749979*(alpha[20]*fUpwind[36]+fUpwind[20]*alpha[36])+(0.159719141249985*alpha[19]+0.2500000000000001*alpha[2])*fUpwind[35]+(0.223606797749979*fUpwind[29]+0.159719141249985*fUpwind[19])*alpha[35]+0.2500000000000001*(fUpwind[2]*alpha[35]+fUpwind[18]*alpha[32])+0.223606797749979*alpha[15]*fUpwind[31]+0.2*alpha[8]*fUpwind[28]+(0.159719141249985*alpha[11]+0.25*alpha[0])*fUpwind[25]+(0.223606797749979*fUpwind[14]+0.159719141249985*fUpwind[11])*alpha[25]+0.25*(fUpwind[0]*alpha[25]+fUpwind[10]*alpha[21]+alpha[9]*fUpwind[19]+fUpwind[9]*alpha[19])+0.223606797749979*(alpha[6]*fUpwind[17]+alpha[5]*fUpwind[16]+fUpwind[5]*alpha[16])+0.2500000000000001*(alpha[4]*fUpwind[11]+fUpwind[4]*alpha[11])+0.223606797749979*(alpha[1]*fUpwind[8]+fUpwind[1]*alpha[8]); 
  Ghat[26] += (0.159719141249985*alpha[33]+0.25*alpha[6])*fUpwind[45]+0.223606797749979*alpha[32]*fUpwind[44]+0.2*alpha[16]*fUpwind[41]+(0.159719141249985*alpha[22]+0.2500000000000001*alpha[3])*fUpwind[38]+(0.159719141249985*alpha[20]+0.2500000000000001*alpha[1])*fUpwind[36]+(0.223606797749979*fUpwind[28]+0.159719141249985*fUpwind[20]+0.2500000000000001*fUpwind[1])*alpha[36]+0.223606797749979*(alpha[19]*fUpwind[35]+fUpwind[19]*alpha[35])+0.2500000000000001*fUpwind[17]*alpha[33]+0.223606797749979*alpha[15]*fUpwind[31]+0.2*alpha[9]*fUpwind[29]+(0.159719141249985*alpha[12]+0.25*alpha[0])*fUpwind[26]+(0.223606797749979*fUpwind[14]+0.159719141249985*fUpwind[12])*alpha[26]+0.25*(fUpwind[0]*alpha[26]+fUpwind[10]*alpha[22]+alpha[8]*fUpwind[20]+fUpwind[8]*alpha[20])+0.223606797749979*(alpha[7]*fUpwind[18]+alpha[5]*fUpwind[16]+fUpwind[5]*alpha[16])+0.2500000000000001*(alpha[4]*fUpwind[12]+fUpwind[4]*alpha[12])+0.223606797749979*(alpha[2]*fUpwind[9]+fUpwind[2]*alpha[9]); 
  Ghat[27] += 0.25*alpha[5]*fUpwind[46]+0.223606797749979*(alpha[33]*fUpwind[45]+alpha[32]*fUpwind[44])+0.2500000000000001*(alpha[2]*fUpwind[40]+alpha[1]*fUpwind[39])+0.223606797749979*(alpha[22]*fUpwind[38]+alpha[21]*fUpwind[37])+0.2500000000000001*alpha[16]*fUpwind[34]+0.223606797749979*alpha[15]*fUpwind[31]+0.25*(alpha[0]*fUpwind[27]+alpha[9]*fUpwind[24]+alpha[8]*fUpwind[23])+0.223606797749979*(alpha[7]*fUpwind[18]+alpha[6]*fUpwind[17])+0.2500000000000001*alpha[4]*fUpwind[13]+0.223606797749979*alpha[3]*fUpwind[10]; 
  Ghat[28] += (0.223606797749979*alpha[32]+0.25*alpha[7])*fUpwind[47]+0.2500000000000001*alpha[15]*fUpwind[43]+(0.223606797749979*alpha[21]+0.2500000000000001*alpha[3])*fUpwind[42]+(0.223606797749979*alpha[19]+0.2500000000000001*alpha[2])*fUpwind[41]+0.223606797749979*(alpha[26]*fUpwind[36]+fUpwind[26]*alpha[36])+0.2*(alpha[16]*fUpwind[35]+fUpwind[16]*alpha[35])+0.25*(alpha[6]*fUpwind[30]+alpha[5]*fUpwind[29])+(0.223606797749979*alpha[11]+0.25*alpha[0])*fUpwind[28]+0.2*(alpha[8]*fUpwind[25]+fUpwind[8]*alpha[25])+0.223606797749979*(alpha[9]*fUpwind[16]+fUpwind[9]*alpha[16])+0.2500000000000001*alpha[1]*fUpwind[14]+0.223606797749979*(alpha[4]*fUpwind[8]+fUpwind[4]*alpha[8]); 
  Ghat[29] += (0.223606797749979*alpha[33]+0.25*alpha[6])*fUpwind[47]+0.223606797749979*alpha[22]*fUpwind[43]+0.2500000000000001*(alpha[3]*fUpwind[43]+alpha[15]*fUpwind[42])+(0.223606797749979*alpha[20]+0.2500000000000001*alpha[1])*fUpwind[41]+0.2*(alpha[16]*fUpwind[36]+fUpwind[16]*alpha[36])+0.223606797749979*(alpha[25]*fUpwind[35]+fUpwind[25]*alpha[35])+0.25*alpha[7]*fUpwind[30]+0.223606797749979*alpha[12]*fUpwind[29]+0.25*(alpha[0]*fUpwind[29]+alpha[5]*fUpwind[28])+0.2*(alpha[9]*fUpwind[26]+fUpwind[9]*alpha[26])+0.223606797749979*(alpha[8]*fUpwind[16]+fUpwind[8]*alpha[16])+0.2500000000000001*alpha[2]*fUpwind[14]+0.223606797749979*(alpha[4]*fUpwind[9]+fUpwind[4]*alpha[9]); 
  Ghat[30] += 0.25*alpha[5]*fUpwind[47]+0.223606797749979*(alpha[36]*fUpwind[45]+alpha[35]*fUpwind[44])+0.2500000000000001*(alpha[2]*fUpwind[43]+alpha[1]*fUpwind[42]+alpha[15]*fUpwind[41])+0.223606797749979*(alpha[26]*fUpwind[38]+alpha[25]*fUpwind[37])+0.223606797749979*alpha[16]*fUpwind[31]+0.25*(alpha[0]*fUpwind[30]+alpha[7]*fUpwind[29]+alpha[6]*fUpwind[28])+0.223606797749979*(alpha[9]*fUpwind[18]+alpha[8]*fUpwind[17])+0.2500000000000001*alpha[3]*fUpwind[14]+0.223606797749979*alpha[4]*fUpwind[10]; 
  Ghat[31] += (0.2*(alpha[26]+alpha[25])+0.223606797749979*alpha[4])*fUpwind[47]+(0.2*(alpha[22]+alpha[21])+0.223606797749979*alpha[3])*fUpwind[46]+(0.2*alpha[19]+0.223606797749979*alpha[2])*fUpwind[45]+(0.2*alpha[20]+0.223606797749979*alpha[1])*fUpwind[44]+(0.2*alpha[36]+0.223606797749979*alpha[8])*fUpwind[43]+(0.2*alpha[35]+0.223606797749979*alpha[9])*fUpwind[42]+(0.2*alpha[33]+0.223606797749979*alpha[6])*fUpwind[40]+0.2*alpha[32]*fUpwind[39]+0.223606797749979*(alpha[7]*fUpwind[39]+alpha[5]*(fUpwind[38]+fUpwind[37]))+(0.2*alpha[32]+0.223606797749979*alpha[7])*fUpwind[36]+(0.2*fUpwind[32]+0.223606797749979*fUpwind[7])*alpha[36]+(0.2*alpha[33]+0.223606797749979*alpha[6])*fUpwind[35]+0.2*fUpwind[33]*alpha[35]+0.223606797749979*(fUpwind[6]*alpha[35]+alpha[9]*fUpwind[33]+fUpwind[9]*alpha[33]+alpha[8]*fUpwind[32]+fUpwind[8]*alpha[32])+(0.223606797749979*(alpha[12]+alpha[11])+0.25*alpha[0])*fUpwind[31]+0.223606797749979*(alpha[16]*fUpwind[30]+alpha[15]*(fUpwind[27]+fUpwind[26])+fUpwind[15]*alpha[26]+alpha[15]*fUpwind[25]+fUpwind[15]*alpha[25]+alpha[16]*fUpwind[22]+fUpwind[16]*alpha[22]+alpha[16]*fUpwind[21]+fUpwind[16]*alpha[21]+fUpwind[18]*alpha[20]+fUpwind[17]*alpha[19])+0.25*(alpha[1]*fUpwind[18]+alpha[2]*fUpwind[17]+alpha[3]*fUpwind[16]+fUpwind[3]*alpha[16]+alpha[4]*fUpwind[15]+fUpwind[4]*alpha[15]+alpha[5]*fUpwind[10]+alpha[6]*fUpwind[9]+fUpwind[6]*alpha[9]+alpha[7]*fUpwind[8]+fUpwind[7]*alpha[8]); 
  Ghat[32] += 0.2*alpha[16]*fUpwind[45]+(0.223606797749979*alpha[26]+0.159719141249985*alpha[25]+0.2500000000000001*alpha[4])*fUpwind[44]+0.223606797749979*alpha[35]*fUpwind[38]+(0.159719141249985*alpha[35]+0.25*alpha[9])*fUpwind[37]+0.2*fUpwind[31]*alpha[36]+0.25*fUpwind[10]*alpha[35]+0.1788854381999831*alpha[33]*fUpwind[34]+0.2*(alpha[6]*fUpwind[34]+alpha[5]*fUpwind[33]+fUpwind[5]*alpha[33])+(0.223606797749979*alpha[12]+0.159719141249985*alpha[11]+0.25*alpha[0])*fUpwind[32]+(0.223606797749979*(fUpwind[13]+fUpwind[12])+0.159719141249985*fUpwind[11]+0.25*fUpwind[0])*alpha[32]+0.223606797749979*alpha[8]*fUpwind[31]+0.2500000000000001*fUpwind[18]*alpha[25]+0.223606797749979*alpha[21]*fUpwind[24]+0.2*alpha[15]*fUpwind[23]+0.223606797749979*(alpha[19]*fUpwind[22]+fUpwind[19]*alpha[22])+(0.159719141249985*alpha[19]+0.2500000000000001*alpha[2])*fUpwind[21]+(0.159719141249985*fUpwind[19]+0.2500000000000001*fUpwind[2])*alpha[21]+0.2*(alpha[15]*fUpwind[20]+fUpwind[15]*alpha[20])+0.2500000000000001*(alpha[3]*fUpwind[19]+fUpwind[3]*alpha[19])+0.223606797749979*(alpha[16]*fUpwind[17]+alpha[1]*fUpwind[15]+fUpwind[1]*alpha[15])+0.25*(alpha[7]*fUpwind[11]+fUpwind[7]*alpha[11])+0.223606797749979*(alpha[5]*fUpwind[6]+fUpwind[5]*alpha[6]); 
  Ghat[33] += (0.159719141249985*alpha[26]+0.223606797749979*alpha[25]+0.2500000000000001*alpha[4])*fUpwind[45]+0.2*alpha[16]*fUpwind[44]+(0.159719141249985*alpha[36]+0.25*alpha[8])*fUpwind[38]+alpha[36]*(0.223606797749979*fUpwind[37]+0.25*fUpwind[10])+0.2*fUpwind[31]*alpha[35]+(0.1788854381999831*alpha[32]+0.2*alpha[7])*fUpwind[34]+(0.159719141249985*alpha[12]+0.223606797749979*alpha[11]+0.25*alpha[0])*fUpwind[33]+(0.223606797749979*fUpwind[13]+0.159719141249985*fUpwind[12]+0.223606797749979*fUpwind[11]+0.25*fUpwind[0])*alpha[33]+0.2*(alpha[5]*fUpwind[32]+fUpwind[5]*alpha[32])+0.223606797749979*alpha[9]*fUpwind[31]+0.2500000000000001*fUpwind[17]*alpha[26]+0.2*alpha[15]*fUpwind[24]+0.223606797749979*alpha[22]*fUpwind[23]+(0.159719141249985*alpha[20]+0.2500000000000001*alpha[1])*fUpwind[22]+(0.159719141249985*fUpwind[20]+0.2500000000000001*fUpwind[1])*alpha[22]+0.223606797749979*(alpha[20]*fUpwind[21]+fUpwind[20]*alpha[21])+0.2500000000000001*(alpha[3]*fUpwind[20]+fUpwind[3]*alpha[20])+0.2*(alpha[15]*fUpwind[19]+fUpwind[15]*alpha[19])+0.223606797749979*(alpha[16]*fUpwind[18]+alpha[2]*fUpwind[15]+fUpwind[2]*alpha[15])+0.25*(alpha[6]*fUpwind[12]+fUpwind[6]*alpha[12])+0.223606797749979*(alpha[5]*fUpwind[7]+fUpwind[5]*alpha[7]); 
  Ghat[34] += (0.223606797749979*(alpha[26]+alpha[25])+0.2500000000000001*alpha[4])*fUpwind[46]+(0.223606797749979*alpha[36]+0.25*alpha[8])*fUpwind[40]+(0.223606797749979*alpha[35]+0.25*alpha[9])*fUpwind[39]+(0.223606797749979*(alpha[12]+alpha[11])+0.25*alpha[0])*fUpwind[34]+(0.1788854381999831*alpha[32]+0.2*alpha[7])*fUpwind[33]+0.1788854381999831*fUpwind[32]*alpha[33]+0.2*(fUpwind[7]*alpha[33]+alpha[6]*fUpwind[32]+fUpwind[6]*alpha[32])+0.2500000000000001*alpha[16]*fUpwind[27]+(0.223606797749979*alpha[20]+0.2500000000000001*alpha[1])*fUpwind[24]+(0.223606797749979*alpha[19]+0.2500000000000001*alpha[2])*fUpwind[23]+0.2*(alpha[15]*fUpwind[22]+fUpwind[15]*alpha[22]+alpha[15]*fUpwind[21]+fUpwind[15]*alpha[21])+0.223606797749979*(alpha[3]*fUpwind[15]+fUpwind[3]*alpha[15])+0.25*alpha[5]*fUpwind[13]+0.223606797749979*(alpha[6]*fUpwind[7]+fUpwind[6]*alpha[7]); 
  Ghat[35] += 0.2*alpha[15]*fUpwind[45]+(0.223606797749979*alpha[22]+0.159719141249985*alpha[21]+0.2500000000000001*alpha[3])*fUpwind[44]+(0.1788854381999831*alpha[36]+0.2*alpha[8])*fUpwind[41]+0.223606797749979*alpha[32]*fUpwind[38]+(0.159719141249985*alpha[32]+0.25*alpha[7])*fUpwind[37]+0.2*(alpha[5]*fUpwind[36]+fUpwind[5]*alpha[36])+(0.223606797749979*alpha[12]+0.159719141249985*alpha[11]+0.25*alpha[0])*fUpwind[35]+(0.223606797749979*(fUpwind[14]+fUpwind[12])+0.159719141249985*fUpwind[11]+0.25*fUpwind[0])*alpha[35]+0.2*fUpwind[31]*alpha[33]+0.25*fUpwind[10]*alpha[32]+0.223606797749979*(alpha[6]*fUpwind[31]+alpha[25]*fUpwind[29])+0.2*alpha[16]*fUpwind[28]+0.223606797749979*(alpha[19]*fUpwind[26]+fUpwind[19]*alpha[26])+(0.159719141249985*alpha[19]+0.2500000000000001*alpha[2])*fUpwind[25]+0.159719141249985*fUpwind[19]*alpha[25]+0.2500000000000001*(fUpwind[2]*alpha[25]+fUpwind[18]*alpha[21])+0.2*(alpha[16]*fUpwind[20]+fUpwind[16]*alpha[20])+0.2500000000000001*(alpha[4]*fUpwind[19]+fUpwind[4]*alpha[19])+0.223606797749979*(alpha[15]*fUpwind[17]+alpha[1]*fUpwind[16]+fUpwind[1]*alpha[16])+0.25*(alpha[9]*fUpwind[11]+fUpwind[9]*alpha[11])+0.223606797749979*(alpha[5]*fUpwind[8]+fUpwind[5]*alpha[8]); 
  Ghat[36] += (0.159719141249985*alpha[22]+0.223606797749979*alpha[21]+0.2500000000000001*alpha[3])*fUpwind[45]+0.2*alpha[15]*fUpwind[44]+(0.1788854381999831*alpha[35]+0.2*alpha[9])*fUpwind[41]+(0.159719141249985*alpha[33]+0.25*alpha[6])*fUpwind[38]+0.223606797749979*alpha[33]*fUpwind[37]+(0.159719141249985*alpha[12]+0.223606797749979*alpha[11]+0.25*alpha[0])*fUpwind[36]+(0.223606797749979*fUpwind[14]+0.159719141249985*fUpwind[12]+0.223606797749979*fUpwind[11]+0.25*fUpwind[0])*alpha[36]+0.2*(alpha[5]*fUpwind[35]+fUpwind[5]*alpha[35])+0.25*fUpwind[10]*alpha[33]+fUpwind[31]*(0.2*alpha[32]+0.223606797749979*alpha[7])+0.2*alpha[16]*fUpwind[29]+0.223606797749979*alpha[26]*fUpwind[28]+(0.159719141249985*alpha[20]+0.2500000000000001*alpha[1])*fUpwind[26]+(0.159719141249985*fUpwind[20]+0.2500000000000001*fUpwind[1])*alpha[26]+0.223606797749979*(alpha[20]*fUpwind[25]+fUpwind[20]*alpha[25])+0.2500000000000001*(fUpwind[17]*alpha[22]+alpha[4]*fUpwind[20]+fUpwind[4]*alpha[20])+0.2*(alpha[16]*fUpwind[19]+fUpwind[16]*alpha[19])+0.223606797749979*(alpha[15]*fUpwind[18]+alpha[2]*fUpwind[16]+fUpwind[2]*alpha[16])+0.25*(alpha[8]*fUpwind[12]+fUpwind[8]*alpha[12])+0.223606797749979*(alpha[5]*fUpwind[9]+fUpwind[5]*alpha[9]); 
  Ghat[37] += 0.2*(alpha[16]*fUpwind[47]+alpha[15]*fUpwind[46])+0.223606797749979*alpha[20]*fUpwind[45]+(0.159719141249985*alpha[19]+0.2500000000000001*alpha[2])*fUpwind[44]+0.223606797749979*alpha[35]*fUpwind[43]+0.2*alpha[8]*fUpwind[42]+0.223606797749979*alpha[32]*fUpwind[40]+0.2*alpha[6]*fUpwind[39]+(0.159719141249985*alpha[11]+0.25*alpha[0])*fUpwind[37]+0.223606797749979*(alpha[33]*fUpwind[36]+fUpwind[33]*alpha[36])+(0.159719141249985*alpha[32]+0.25*alpha[7])*fUpwind[35]+0.159719141249985*fUpwind[32]*alpha[35]+0.25*(fUpwind[7]*alpha[35]+alpha[9]*fUpwind[32]+fUpwind[9]*alpha[32])+0.223606797749979*(alpha[5]*fUpwind[31]+alpha[25]*fUpwind[30]+alpha[21]*fUpwind[27])+(0.159719141249985*alpha[21]+0.2500000000000001*alpha[3])*fUpwind[25]+0.159719141249985*fUpwind[21]*alpha[25]+0.2500000000000001*(fUpwind[3]*alpha[25]+alpha[4]*fUpwind[21]+fUpwind[4]*alpha[21]+fUpwind[18]*alpha[19])+0.223606797749979*(alpha[1]*fUpwind[17]+alpha[15]*fUpwind[16]+fUpwind[15]*alpha[16])+0.25*fUpwind[10]*alpha[11]+0.223606797749979*(alpha[6]*fUpwind[8]+fUpwind[6]*alpha[8]); 
  Ghat[38] += 0.2*(alpha[16]*fUpwind[47]+alpha[15]*fUpwind[46])+(0.159719141249985*alpha[20]+0.2500000000000001*alpha[1])*fUpwind[45]+0.223606797749979*alpha[19]*fUpwind[44]+0.2*alpha[9]*fUpwind[43]+0.223606797749979*alpha[36]*fUpwind[42]+0.2*alpha[7]*fUpwind[40]+0.223606797749979*alpha[33]*fUpwind[39]+(0.159719141249985*alpha[12]+0.25*alpha[0])*fUpwind[38]+(0.159719141249985*alpha[33]+0.25*alpha[6])*fUpwind[36]+(0.159719141249985*fUpwind[33]+0.25*fUpwind[6])*alpha[36]+0.223606797749979*(alpha[32]*fUpwind[35]+fUpwind[32]*alpha[35])+0.25*(alpha[8]*fUpwind[33]+fUpwind[8]*alpha[33])+0.223606797749979*(alpha[5]*fUpwind[31]+alpha[26]*fUpwind[30]+alpha[22]*fUpwind[27])+(0.159719141249985*alpha[22]+0.2500000000000001*alpha[3])*fUpwind[26]+0.159719141249985*fUpwind[22]*alpha[26]+0.2500000000000001*(fUpwind[3]*alpha[26]+alpha[4]*fUpwind[22]+fUpwind[4]*alpha[22]+fUpwind[17]*alpha[20])+0.223606797749979*(alpha[2]*fUpwind[18]+alpha[15]*fUpwind[16]+fUpwind[15]*alpha[16])+0.25*fUpwind[10]*alpha[12]+0.223606797749979*(alpha[7]*fUpwind[9]+fUpwind[7]*alpha[9]); 
  Ghat[39] += (0.223606797749979*alpha[19]+0.2500000000000001*alpha[2])*fUpwind[46]+0.223606797749979*alpha[22]*fUpwind[45]+0.2*alpha[15]*fUpwind[44]+0.25*alpha[5]*fUpwind[40]+(0.223606797749979*alpha[11]+0.25*alpha[0])*fUpwind[39]+0.223606797749979*alpha[33]*fUpwind[38]+0.2*alpha[6]*fUpwind[37]+fUpwind[34]*(0.223606797749979*alpha[35]+0.25*alpha[9])+fUpwind[31]*(0.2*alpha[32]+0.223606797749979*alpha[7])+0.2500000000000001*alpha[1]*fUpwind[27]+0.223606797749979*fUpwind[23]*alpha[25]+0.2500000000000001*(alpha[16]*fUpwind[24]+alpha[4]*fUpwind[23])+0.2*fUpwind[17]*alpha[21]+0.223606797749979*(alpha[15]*fUpwind[18]+alpha[3]*fUpwind[17])+0.25*alpha[8]*fUpwind[13]+0.223606797749979*alpha[6]*fUpwind[10]; 
  Ghat[40] += (0.223606797749979*alpha[20]+0.2500000000000001*alpha[1])*fUpwind[46]+0.2*alpha[15]*fUpwind[45]+0.223606797749979*(alpha[21]*fUpwind[44]+alpha[12]*fUpwind[40])+0.25*(alpha[0]*fUpwind[40]+alpha[5]*fUpwind[39])+0.2*alpha[7]*fUpwind[38]+0.223606797749979*alpha[32]*fUpwind[37]+fUpwind[34]*(0.223606797749979*alpha[36]+0.25*alpha[8])+fUpwind[31]*(0.2*alpha[33]+0.223606797749979*alpha[6])+0.2500000000000001*alpha[2]*fUpwind[27]+0.223606797749979*fUpwind[24]*alpha[26]+0.2500000000000001*(alpha[4]*fUpwind[24]+alpha[16]*fUpwind[23])+0.2*fUpwind[18]*alpha[22]+0.223606797749979*(alpha[3]*fUpwind[18]+alpha[15]*fUpwind[17])+0.25*alpha[9]*fUpwind[13]+0.223606797749979*alpha[7]*fUpwind[10]; 
  Ghat[41] += (0.223606797749979*(alpha[22]+alpha[21])+0.2500000000000001*alpha[3])*fUpwind[47]+(0.223606797749979*alpha[33]+0.25*alpha[6])*fUpwind[43]+(0.223606797749979*alpha[32]+0.25*alpha[7])*fUpwind[42]+(0.223606797749979*(alpha[12]+alpha[11])+0.25*alpha[0])*fUpwind[41]+(0.1788854381999831*alpha[35]+0.2*alpha[9])*fUpwind[36]+0.1788854381999831*fUpwind[35]*alpha[36]+0.2*(fUpwind[9]*alpha[36]+alpha[8]*fUpwind[35]+fUpwind[8]*alpha[35])+0.2500000000000001*alpha[15]*fUpwind[30]+(0.223606797749979*alpha[20]+0.2500000000000001*alpha[1])*fUpwind[29]+(0.223606797749979*alpha[19]+0.2500000000000001*alpha[2])*fUpwind[28]+0.2*(alpha[16]*fUpwind[26]+fUpwind[16]*alpha[26]+alpha[16]*fUpwind[25]+fUpwind[16]*alpha[25])+0.223606797749979*(alpha[4]*fUpwind[16]+fUpwind[4]*alpha[16])+0.25*alpha[5]*fUpwind[14]+0.223606797749979*(alpha[8]*fUpwind[9]+fUpwind[8]*alpha[9]); 
  Ghat[42] += (0.223606797749979*alpha[19]+0.2500000000000001*alpha[2])*fUpwind[47]+0.223606797749979*alpha[26]*fUpwind[45]+0.2*alpha[16]*fUpwind[44]+0.25*alpha[5]*fUpwind[43]+(0.223606797749979*alpha[11]+0.25*alpha[0])*fUpwind[42]+(0.223606797749979*alpha[32]+0.25*alpha[7])*fUpwind[41]+0.223606797749979*alpha[36]*fUpwind[38]+0.2*alpha[8]*fUpwind[37]+fUpwind[31]*(0.2*alpha[35]+0.223606797749979*alpha[9])+0.2500000000000001*(alpha[1]*fUpwind[30]+alpha[15]*fUpwind[29])+(0.223606797749979*alpha[21]+0.2500000000000001*alpha[3])*fUpwind[28]+0.2*fUpwind[17]*alpha[25]+0.223606797749979*(alpha[16]*fUpwind[18]+alpha[4]*fUpwind[17])+0.25*alpha[6]*fUpwind[14]+0.223606797749979*alpha[8]*fUpwind[10]; 
  Ghat[43] += (0.223606797749979*alpha[20]+0.2500000000000001*alpha[1])*fUpwind[47]+0.2*alpha[16]*fUpwind[45]+0.223606797749979*(alpha[25]*fUpwind[44]+alpha[12]*fUpwind[43])+0.25*(alpha[0]*fUpwind[43]+alpha[5]*fUpwind[42])+(0.223606797749979*alpha[33]+0.25*alpha[6])*fUpwind[41]+0.2*alpha[9]*fUpwind[38]+0.223606797749979*alpha[35]*fUpwind[37]+fUpwind[31]*(0.2*alpha[36]+0.223606797749979*alpha[8])+0.2500000000000001*alpha[2]*fUpwind[30]+0.223606797749979*alpha[22]*fUpwind[29]+0.2500000000000001*(alpha[3]*fUpwind[29]+alpha[15]*fUpwind[28])+0.2*fUpwind[18]*alpha[26]+0.223606797749979*(alpha[4]*fUpwind[18]+alpha[16]*fUpwind[17])+0.25*alpha[7]*fUpwind[14]+0.223606797749979*alpha[9]*fUpwind[10]; 
  Ghat[44] += (0.1788854381999831*alpha[36]+0.2*alpha[8])*fUpwind[47]+0.1788854381999831*alpha[33]*fUpwind[46]+0.2*(alpha[6]*fUpwind[46]+alpha[5]*fUpwind[45])+(0.223606797749979*alpha[12]+0.159719141249985*alpha[11]+0.25*alpha[0])*fUpwind[44]+0.223606797749979*alpha[25]*fUpwind[43]+0.2*alpha[16]*fUpwind[42]+0.223606797749979*alpha[21]*fUpwind[40]+0.2*alpha[15]*fUpwind[39]+0.223606797749979*alpha[19]*fUpwind[38]+(0.159719141249985*alpha[19]+0.2500000000000001*alpha[2])*fUpwind[37]+0.2*(alpha[15]*fUpwind[36]+fUpwind[15]*alpha[36])+(0.223606797749979*alpha[22]+0.159719141249985*alpha[21]+0.2500000000000001*alpha[3])*fUpwind[35]+(0.223606797749979*(fUpwind[30]+fUpwind[22])+0.159719141249985*fUpwind[21]+0.2500000000000001*fUpwind[3])*alpha[35]+0.2*(alpha[16]*fUpwind[33]+fUpwind[16]*alpha[33])+(0.223606797749979*alpha[26]+0.159719141249985*alpha[25]+0.2500000000000001*alpha[4])*fUpwind[32]+(0.223606797749979*(fUpwind[27]+fUpwind[26])+0.159719141249985*fUpwind[25]+0.2500000000000001*fUpwind[4])*alpha[32]+(0.2*alpha[20]+0.223606797749979*alpha[1])*fUpwind[31]+0.25*(alpha[7]*fUpwind[25]+fUpwind[7]*alpha[25]+alpha[9]*fUpwind[21]+fUpwind[9]*alpha[21]+fUpwind[10]*alpha[19])+0.2500000000000001*alpha[11]*fUpwind[18]+0.223606797749979*(alpha[5]*fUpwind[17]+alpha[6]*fUpwind[16]+fUpwind[6]*alpha[16]+alpha[8]*fUpwind[15]+fUpwind[8]*alpha[15]); 
  Ghat[45] += (0.1788854381999831*alpha[35]+0.2*alpha[9])*fUpwind[47]+(0.1788854381999831*alpha[32]+0.2*alpha[7])*fUpwind[46]+(0.159719141249985*alpha[12]+0.223606797749979*alpha[11]+0.25*alpha[0])*fUpwind[45]+0.2*(alpha[5]*fUpwind[44]+alpha[16]*fUpwind[43])+0.223606797749979*alpha[26]*fUpwind[42]+0.2*alpha[15]*fUpwind[40]+0.223606797749979*alpha[22]*fUpwind[39]+(0.159719141249985*alpha[20]+0.2500000000000001*alpha[1])*fUpwind[38]+0.223606797749979*alpha[20]*fUpwind[37]+(0.159719141249985*alpha[22]+0.223606797749979*alpha[21]+0.2500000000000001*alpha[3])*fUpwind[36]+(0.223606797749979*fUpwind[30]+0.159719141249985*fUpwind[22]+0.223606797749979*fUpwind[21]+0.2500000000000001*fUpwind[3])*alpha[36]+0.2*(alpha[15]*fUpwind[35]+fUpwind[15]*alpha[35])+(0.159719141249985*alpha[26]+0.223606797749979*alpha[25]+0.2500000000000001*alpha[4])*fUpwind[33]+(0.223606797749979*fUpwind[27]+0.159719141249985*fUpwind[26]+0.223606797749979*fUpwind[25]+0.2500000000000001*fUpwind[4])*alpha[33]+0.2*(alpha[16]*fUpwind[32]+fUpwind[16]*alpha[32])+(0.2*alpha[19]+0.223606797749979*alpha[2])*fUpwind[31]+0.25*(alpha[6]*fUpwind[26]+fUpwind[6]*alpha[26]+alpha[8]*fUpwind[22]+fUpwind[8]*alpha[22]+fUpwind[10]*alpha[20])+0.223606797749979*alpha[5]*fUpwind[18]+0.2500000000000001*alpha[12]*fUpwind[17]+0.223606797749979*(alpha[7]*fUpwind[16]+fUpwind[7]*alpha[16]+alpha[9]*fUpwind[15]+fUpwind[9]*alpha[15]); 
  Ghat[46] += (0.223606797749979*(alpha[12]+alpha[11])+0.25*alpha[0])*fUpwind[46]+(0.1788854381999831*alpha[32]+0.2*alpha[7])*fUpwind[45]+(0.1788854381999831*alpha[33]+0.2*alpha[6])*fUpwind[44]+(0.223606797749979*alpha[20]+0.2500000000000001*alpha[1])*fUpwind[40]+(0.223606797749979*alpha[19]+0.2500000000000001*alpha[2])*fUpwind[39]+0.2*alpha[15]*(fUpwind[38]+fUpwind[37])+0.223606797749979*(fUpwind[24]*alpha[36]+fUpwind[23]*alpha[35])+(0.223606797749979*(alpha[26]+alpha[25])+0.2500000000000001*alpha[4])*fUpwind[34]+0.2*(fUpwind[18]*alpha[33]+fUpwind[17]*alpha[32])+(0.2*(alpha[22]+alpha[21])+0.223606797749979*alpha[3])*fUpwind[31]+0.25*(alpha[5]*fUpwind[27]+alpha[8]*fUpwind[24]+alpha[9]*fUpwind[23])+0.223606797749979*(alpha[6]*fUpwind[18]+alpha[7]*fUpwind[17])+0.2500000000000001*fUpwind[13]*alpha[16]+0.223606797749979*fUpwind[10]*alpha[15]; 
  Ghat[47] += (0.223606797749979*(alpha[12]+alpha[11])+0.25*alpha[0])*fUpwind[47]+(0.1788854381999831*alpha[35]+0.2*alpha[9])*fUpwind[45]+(0.1788854381999831*alpha[36]+0.2*alpha[8])*fUpwind[44]+(0.223606797749979*alpha[20]+0.2500000000000001*alpha[1])*fUpwind[43]+(0.223606797749979*alpha[19]+0.2500000000000001*alpha[2])*fUpwind[42]+(0.223606797749979*(alpha[22]+alpha[21])+0.2500000000000001*alpha[3])*fUpwind[41]+0.2*(alpha[16]*(fUpwind[38]+fUpwind[37])+fUpwind[18]*alpha[36]+fUpwind[17]*alpha[35])+0.223606797749979*(fUpwind[29]*alpha[33]+fUpwind[28]*alpha[32])+(0.2*(alpha[26]+alpha[25])+0.223606797749979*alpha[4])*fUpwind[31]+0.25*(alpha[5]*fUpwind[30]+alpha[6]*fUpwind[29]+alpha[7]*fUpwind[28])+0.223606797749979*(alpha[8]*fUpwind[18]+alpha[9]*fUpwind[17]+fUpwind[10]*alpha[16])+0.2500000000000001*fUpwind[14]*alpha[15]; 

  out[0] += -0.7071067811865475*Ghat[0]*dv10; 
  out[1] += -0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -0.7071067811865475*Ghat[2]*dv10; 
  out[3] += -1.224744871391589*Ghat[0]*dv10; 
  out[4] += -0.7071067811865475*Ghat[3]*dv10; 
  out[5] += -0.7071067811865475*Ghat[4]*dv10; 
  out[6] += -0.7071067811865475*Ghat[5]*dv10; 
  out[7] += -1.224744871391589*Ghat[1]*dv10; 
  out[8] += -1.224744871391589*Ghat[2]*dv10; 
  out[9] += -0.7071067811865475*Ghat[6]*dv10; 
  out[10] += -0.7071067811865475*Ghat[7]*dv10; 
  out[11] += -1.224744871391589*Ghat[3]*dv10; 
  out[12] += -0.7071067811865475*Ghat[8]*dv10; 
  out[13] += -0.7071067811865475*Ghat[9]*dv10; 
  out[14] += -1.224744871391589*Ghat[4]*dv10; 
  out[15] += -0.7071067811865475*Ghat[10]*dv10; 
  out[16] += -0.7071067811865475*Ghat[11]*dv10; 
  out[17] += -0.7071067811865475*Ghat[12]*dv10; 
  out[18] += -1.58113883008419*Ghat[0]*dv10; 
  out[19] += -0.7071067811865475*Ghat[13]*dv10; 
  out[20] += -0.7071067811865475*Ghat[14]*dv10; 
  out[21] += -1.224744871391589*Ghat[5]*dv10; 
  out[22] += -0.7071067811865475*Ghat[15]*dv10; 
  out[23] += -1.224744871391589*Ghat[6]*dv10; 
  out[24] += -1.224744871391589*Ghat[7]*dv10; 
  out[25] += -0.7071067811865475*Ghat[16]*dv10; 
  out[26] += -1.224744871391589*Ghat[8]*dv10; 
  out[27] += -1.224744871391589*Ghat[9]*dv10; 
  out[28] += -0.7071067811865475*Ghat[17]*dv10; 
  out[29] += -0.7071067811865475*Ghat[18]*dv10; 
  out[30] += -1.224744871391589*Ghat[10]*dv10; 
  out[31] += -0.7071067811865475*Ghat[19]*dv10; 
  out[32] += -0.7071067811865475*Ghat[20]*dv10; 
  out[33] += -1.224744871391589*Ghat[11]*dv10; 
  out[34] += -1.224744871391589*Ghat[12]*dv10; 
  out[35] += -1.58113883008419*Ghat[1]*dv10; 
  out[36] += -1.58113883008419*Ghat[2]*dv10; 
  out[37] += -0.7071067811865475*Ghat[21]*dv10; 
  out[38] += -0.7071067811865475*Ghat[22]*dv10; 
  out[39] += -1.58113883008419*Ghat[3]*dv10; 
  out[40] += -0.7071067811865475*Ghat[23]*dv10; 
  out[41] += -0.7071067811865475*Ghat[24]*dv10; 
  out[42] += -1.224744871391589*Ghat[13]*dv10; 
  out[43] += -0.7071067811865475*Ghat[25]*dv10; 
  out[44] += -0.7071067811865475*Ghat[26]*dv10; 
  out[45] += -1.58113883008419*Ghat[4]*dv10; 
  out[46] += -0.7071067811865475*Ghat[27]*dv10; 
  out[47] += -0.7071067811865475*Ghat[28]*dv10; 
  out[48] += -0.7071067811865475*Ghat[29]*dv10; 
  out[49] += -1.224744871391589*Ghat[14]*dv10; 
  out[50] += -0.7071067811865475*Ghat[30]*dv10; 
  out[51] += -1.224744871391589*Ghat[15]*dv10; 
  out[52] += -1.224744871391589*Ghat[16]*dv10; 
  out[53] += -0.7071067811865475*Ghat[31]*dv10; 
  out[54] += -1.224744871391589*Ghat[17]*dv10; 
  out[55] += -1.224744871391589*Ghat[18]*dv10; 
  out[56] += -1.224744871391589*Ghat[19]*dv10; 
  out[57] += -1.224744871391589*Ghat[20]*dv10; 
  out[58] += -1.58113883008419*Ghat[5]*dv10; 
  out[59] += -0.7071067811865475*Ghat[32]*dv10; 
  out[60] += -0.7071067811865475*Ghat[33]*dv10; 
  out[61] += -1.224744871391589*Ghat[21]*dv10; 
  out[62] += -1.224744871391589*Ghat[22]*dv10; 
  out[63] += -1.58113883008419*Ghat[6]*dv10; 
  out[64] += -1.58113883008419*Ghat[7]*dv10; 
  out[65] += -0.7071067811865475*Ghat[34]*dv10; 
  out[66] += -1.224744871391589*Ghat[23]*dv10; 
  out[67] += -1.224744871391589*Ghat[24]*dv10; 
  out[68] += -0.7071067811865475*Ghat[35]*dv10; 
  out[69] += -0.7071067811865475*Ghat[36]*dv10; 
  out[70] += -1.224744871391589*Ghat[25]*dv10; 
  out[71] += -1.224744871391589*Ghat[26]*dv10; 
  out[72] += -1.58113883008419*Ghat[8]*dv10; 
  out[73] += -1.58113883008419*Ghat[9]*dv10; 
  out[74] += -0.7071067811865475*Ghat[37]*dv10; 
  out[75] += -0.7071067811865475*Ghat[38]*dv10; 
  out[76] += -1.58113883008419*Ghat[10]*dv10; 
  out[77] += -0.7071067811865475*Ghat[39]*dv10; 
  out[78] += -0.7071067811865475*Ghat[40]*dv10; 
  out[79] += -1.224744871391589*Ghat[27]*dv10; 
  out[80] += -0.7071067811865475*Ghat[41]*dv10; 
  out[81] += -1.224744871391589*Ghat[28]*dv10; 
  out[82] += -1.224744871391589*Ghat[29]*dv10; 
  out[83] += -0.7071067811865475*Ghat[42]*dv10; 
  out[84] += -0.7071067811865475*Ghat[43]*dv10; 
  out[85] += -1.224744871391589*Ghat[30]*dv10; 
  out[86] += -1.224744871391589*Ghat[31]*dv10; 
  out[87] += -1.224744871391589*Ghat[32]*dv10; 
  out[88] += -1.224744871391589*Ghat[33]*dv10; 
  out[89] += -1.58113883008419*Ghat[15]*dv10; 
  out[90] += -1.224744871391589*Ghat[34]*dv10; 
  out[91] += -1.224744871391589*Ghat[35]*dv10; 
  out[92] += -1.224744871391589*Ghat[36]*dv10; 
  out[93] += -1.58113883008419*Ghat[16]*dv10; 
  out[94] += -0.7071067811865475*Ghat[44]*dv10; 
  out[95] += -0.7071067811865475*Ghat[45]*dv10; 
  out[96] += -1.224744871391589*Ghat[37]*dv10; 
  out[97] += -1.224744871391589*Ghat[38]*dv10; 
  out[98] += -1.58113883008419*Ghat[17]*dv10; 
  out[99] += -1.58113883008419*Ghat[18]*dv10; 
  out[100] += -0.7071067811865475*Ghat[46]*dv10; 
  out[101] += -1.224744871391589*Ghat[39]*dv10; 
  out[102] += -1.224744871391589*Ghat[40]*dv10; 
  out[103] += -1.224744871391589*Ghat[41]*dv10; 
  out[104] += -0.7071067811865475*Ghat[47]*dv10; 
  out[105] += -1.224744871391589*Ghat[42]*dv10; 
  out[106] += -1.224744871391589*Ghat[43]*dv10; 
  out[107] += -1.224744871391589*Ghat[44]*dv10; 
  out[108] += -1.224744871391589*Ghat[45]*dv10; 
  out[109] += -1.58113883008419*Ghat[31]*dv10; 
  out[110] += -1.224744871391589*Ghat[46]*dv10; 
  out[111] += -1.224744871391589*Ghat[47]*dv10; 

  } else { 

  if (0.4024922359499621*(alpha[36]+alpha[35]+alpha[33]+alpha[32])-0.2999999999999999*(alpha[26]+alpha[25]+alpha[22]+alpha[21])-0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]-0.6037383539249431*(alpha[16]+alpha[15])+0.223606797749979*(alpha[12]+alpha[11])+0.45*(alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5])-0.3354101966249685*(alpha[4]+alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[0] = ser_5x_p2_surfx3_quad_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_5x_p2_surfx3_quad_0_l(fSkin); 
  } 
  if ((-0.5031152949374527*(alpha[35]+alpha[32]))-0.2999999999999999*alpha[26]+0.375*alpha[25]-0.2999999999999999*alpha[22]+0.375*(alpha[21]+alpha[19])+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]+0.45*(alpha[9]+alpha[7])-0.3354101966249685*(alpha[4]+alpha[3]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad[1] = ser_5x_p2_surfx3_quad_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = ser_5x_p2_surfx3_quad_1_l(fSkin); 
  } 
  if ((-0.4024922359499621*alpha[36])+0.4024922359499621*alpha[35]-0.4024922359499621*alpha[33]+0.4024922359499621*alpha[32]-0.2999999999999999*(alpha[26]+alpha[25]+alpha[22]+alpha[21])+0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]+0.6037383539249431*(alpha[16]+alpha[15])+0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[9]-0.45*alpha[8]+0.45*alpha[7]-0.45*(alpha[6]+alpha[5])-0.3354101966249685*(alpha[4]+alpha[3]+alpha[2])+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[2] = ser_5x_p2_surfx3_quad_2_r(fEdge); 
  } else { 
    fUpwindQuad[2] = ser_5x_p2_surfx3_quad_2_l(fSkin); 
  } 
  if ((-0.5031152949374527*(alpha[36]+alpha[33]))+0.375*alpha[26]-0.2999999999999999*alpha[25]+0.375*alpha[22]-0.2999999999999999*alpha[21]+0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]+0.45*(alpha[8]+alpha[6])-0.3354101966249685*(alpha[4]+alpha[3]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[3] = ser_5x_p2_surfx3_quad_3_r(fEdge); 
  } else { 
    fUpwindQuad[3] = ser_5x_p2_surfx3_quad_3_l(fSkin); 
  } 
  if (0.375*(alpha[26]+alpha[25]+alpha[22]+alpha[21])-0.2795084971874737*(alpha[12]+alpha[11])-0.3354101966249685*(alpha[4]+alpha[3])+0.25*alpha[0] > 0) { 
    fUpwindQuad[4] = ser_5x_p2_surfx3_quad_4_r(fEdge); 
  } else { 
    fUpwindQuad[4] = ser_5x_p2_surfx3_quad_4_l(fSkin); 
  } 
  if (0.5031152949374527*(alpha[36]+alpha[33])+0.375*alpha[26]-0.2999999999999999*alpha[25]+0.375*alpha[22]-0.2999999999999999*alpha[21]-0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]-0.45*(alpha[8]+alpha[6])-0.3354101966249685*(alpha[4]+alpha[3])+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[5] = ser_5x_p2_surfx3_quad_5_r(fEdge); 
  } else { 
    fUpwindQuad[5] = ser_5x_p2_surfx3_quad_5_l(fSkin); 
  } 
  if (0.4024922359499621*alpha[36]-0.4024922359499621*alpha[35]+0.4024922359499621*alpha[33]-0.4024922359499621*alpha[32]-0.2999999999999999*(alpha[26]+alpha[25]+alpha[22]+alpha[21])-0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]+0.6037383539249431*(alpha[16]+alpha[15])+0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[9]+0.45*alpha[8]-0.45*alpha[7]+0.45*alpha[6]-0.45*alpha[5]-0.3354101966249685*(alpha[4]+alpha[3])+0.3354101966249685*alpha[2]-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[6] = ser_5x_p2_surfx3_quad_6_r(fEdge); 
  } else { 
    fUpwindQuad[6] = ser_5x_p2_surfx3_quad_6_l(fSkin); 
  } 
  if (0.5031152949374527*(alpha[35]+alpha[32])-0.2999999999999999*alpha[26]+0.375*alpha[25]-0.2999999999999999*alpha[22]+0.375*alpha[21]-0.375*alpha[19]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]-0.45*(alpha[9]+alpha[7])-0.3354101966249685*(alpha[4]+alpha[3])+0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[7] = ser_5x_p2_surfx3_quad_7_r(fEdge); 
  } else { 
    fUpwindQuad[7] = ser_5x_p2_surfx3_quad_7_l(fSkin); 
  } 
  if ((-0.4024922359499621*(alpha[36]+alpha[35]+alpha[33]+alpha[32]))-0.2999999999999999*(alpha[26]+alpha[25]+alpha[22]+alpha[21])+0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]-0.6037383539249431*(alpha[16]+alpha[15])+0.223606797749979*(alpha[12]+alpha[11])-0.45*(alpha[9]+alpha[8]+alpha[7]+alpha[6])+0.45*alpha[5]-0.3354101966249685*(alpha[4]+alpha[3])+0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[8] = ser_5x_p2_surfx3_quad_8_r(fEdge); 
  } else { 
    fUpwindQuad[8] = ser_5x_p2_surfx3_quad_8_l(fSkin); 
  } 
  if (0.4024922359499621*(alpha[36]+alpha[35])-0.2999999999999999*(alpha[26]+alpha[25])-0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]-0.6037383539249431*alpha[16]+0.223606797749979*(alpha[12]+alpha[11])+0.45*(alpha[9]+alpha[8]+alpha[5])-0.3354101966249685*(alpha[4]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[9] = ser_5x_p2_surfx3_quad_9_r(fEdge); 
  } else { 
    fUpwindQuad[9] = ser_5x_p2_surfx3_quad_9_l(fSkin); 
  } 
  if ((-0.5031152949374527*alpha[35])-0.2999999999999999*alpha[26]+0.375*(alpha[25]+alpha[19])+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]+0.45*alpha[9]-0.3354101966249685*(alpha[4]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad[10] = ser_5x_p2_surfx3_quad_10_r(fEdge); 
  } else { 
    fUpwindQuad[10] = ser_5x_p2_surfx3_quad_10_l(fSkin); 
  } 
  if ((-0.4024922359499621*alpha[36])+0.4024922359499621*alpha[35]-0.2999999999999999*(alpha[26]+alpha[25])+0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]+0.6037383539249431*alpha[16]+0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[9]-0.45*(alpha[8]+alpha[5])-0.3354101966249685*(alpha[4]+alpha[2])+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[11] = ser_5x_p2_surfx3_quad_11_r(fEdge); 
  } else { 
    fUpwindQuad[11] = ser_5x_p2_surfx3_quad_11_l(fSkin); 
  } 
  if ((-0.5031152949374527*alpha[36])+0.375*alpha[26]-0.2999999999999999*alpha[25]+0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]+0.45*alpha[8]-0.3354101966249685*(alpha[4]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[12] = ser_5x_p2_surfx3_quad_12_r(fEdge); 
  } else { 
    fUpwindQuad[12] = ser_5x_p2_surfx3_quad_12_l(fSkin); 
  } 
  if (0.375*(alpha[26]+alpha[25])-0.2795084971874737*(alpha[12]+alpha[11])-0.3354101966249685*alpha[4]+0.25*alpha[0] > 0) { 
    fUpwindQuad[13] = ser_5x_p2_surfx3_quad_13_r(fEdge); 
  } else { 
    fUpwindQuad[13] = ser_5x_p2_surfx3_quad_13_l(fSkin); 
  } 
  if (0.5031152949374527*alpha[36]+0.375*alpha[26]-0.2999999999999999*alpha[25]-0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]-0.45*alpha[8]-0.3354101966249685*alpha[4]+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[14] = ser_5x_p2_surfx3_quad_14_r(fEdge); 
  } else { 
    fUpwindQuad[14] = ser_5x_p2_surfx3_quad_14_l(fSkin); 
  } 
  if (0.4024922359499621*alpha[36]-0.4024922359499621*alpha[35]-0.2999999999999999*(alpha[26]+alpha[25])-0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]+0.6037383539249431*alpha[16]+0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[9]+0.45*alpha[8]-0.45*alpha[5]-0.3354101966249685*alpha[4]+0.3354101966249685*alpha[2]-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[15] = ser_5x_p2_surfx3_quad_15_r(fEdge); 
  } else { 
    fUpwindQuad[15] = ser_5x_p2_surfx3_quad_15_l(fSkin); 
  } 
  if (0.5031152949374527*alpha[35]-0.2999999999999999*alpha[26]+0.375*alpha[25]-0.375*alpha[19]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]-0.45*alpha[9]-0.3354101966249685*alpha[4]+0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[16] = ser_5x_p2_surfx3_quad_16_r(fEdge); 
  } else { 
    fUpwindQuad[16] = ser_5x_p2_surfx3_quad_16_l(fSkin); 
  } 
  if ((-0.4024922359499621*(alpha[36]+alpha[35]))-0.2999999999999999*(alpha[26]+alpha[25])+0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]-0.6037383539249431*alpha[16]+0.223606797749979*(alpha[12]+alpha[11])-0.45*(alpha[9]+alpha[8])+0.45*alpha[5]-0.3354101966249685*alpha[4]+0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[17] = ser_5x_p2_surfx3_quad_17_r(fEdge); 
  } else { 
    fUpwindQuad[17] = ser_5x_p2_surfx3_quad_17_l(fSkin); 
  } 
  if (0.4024922359499621*(alpha[36]+alpha[35])-0.4024922359499621*(alpha[33]+alpha[32])-0.2999999999999999*(alpha[26]+alpha[25])+0.2999999999999999*(alpha[22]+alpha[21])-0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]-0.6037383539249431*alpha[16]+0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])+0.45*(alpha[9]+alpha[8])-0.45*(alpha[7]+alpha[6])+0.45*alpha[5]-0.3354101966249685*alpha[4]+0.3354101966249685*alpha[3]-0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[18] = ser_5x_p2_surfx3_quad_18_r(fEdge); 
  } else { 
    fUpwindQuad[18] = ser_5x_p2_surfx3_quad_18_l(fSkin); 
  } 
  if ((-0.5031152949374527*alpha[35])+0.5031152949374527*alpha[32]-0.2999999999999999*alpha[26]+0.375*alpha[25]+0.2999999999999999*alpha[22]-0.375*alpha[21]+0.375*alpha[19]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]+0.45*alpha[9]-0.45*alpha[7]-0.3354101966249685*alpha[4]+0.3354101966249685*alpha[3]-0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[19] = ser_5x_p2_surfx3_quad_19_r(fEdge); 
  } else { 
    fUpwindQuad[19] = ser_5x_p2_surfx3_quad_19_l(fSkin); 
  } 
  if ((-0.4024922359499621*alpha[36])+0.4024922359499621*(alpha[35]+alpha[33])-0.4024922359499621*alpha[32]-0.2999999999999999*(alpha[26]+alpha[25])+0.2999999999999999*(alpha[22]+alpha[21])+0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]+0.6037383539249431*alpha[16]-0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[9]-0.45*(alpha[8]+alpha[7])+0.45*alpha[6]-0.45*alpha[5]-0.3354101966249685*alpha[4]+0.3354101966249685*alpha[3]-0.3354101966249685*alpha[2]+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[20] = ser_5x_p2_surfx3_quad_20_r(fEdge); 
  } else { 
    fUpwindQuad[20] = ser_5x_p2_surfx3_quad_20_l(fSkin); 
  } 
  if ((-0.5031152949374527*alpha[36])+0.5031152949374527*alpha[33]+0.375*alpha[26]-0.2999999999999999*alpha[25]-0.375*alpha[22]+0.2999999999999999*alpha[21]+0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]+0.45*alpha[8]-0.45*alpha[6]-0.3354101966249685*alpha[4]+0.3354101966249685*alpha[3]-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[21] = ser_5x_p2_surfx3_quad_21_r(fEdge); 
  } else { 
    fUpwindQuad[21] = ser_5x_p2_surfx3_quad_21_l(fSkin); 
  } 
  if (0.375*(alpha[26]+alpha[25])-0.375*(alpha[22]+alpha[21])-0.2795084971874737*(alpha[12]+alpha[11])-0.3354101966249685*alpha[4]+0.3354101966249685*alpha[3]+0.25*alpha[0] > 0) { 
    fUpwindQuad[22] = ser_5x_p2_surfx3_quad_22_r(fEdge); 
  } else { 
    fUpwindQuad[22] = ser_5x_p2_surfx3_quad_22_l(fSkin); 
  } 
  if (0.5031152949374527*alpha[36]-0.5031152949374527*alpha[33]+0.375*alpha[26]-0.2999999999999999*alpha[25]-0.375*alpha[22]+0.2999999999999999*alpha[21]-0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]-0.45*alpha[8]+0.45*alpha[6]-0.3354101966249685*alpha[4]+0.3354101966249685*(alpha[3]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[23] = ser_5x_p2_surfx3_quad_23_r(fEdge); 
  } else { 
    fUpwindQuad[23] = ser_5x_p2_surfx3_quad_23_l(fSkin); 
  } 
  if (0.4024922359499621*alpha[36]-0.4024922359499621*(alpha[35]+alpha[33])+0.4024922359499621*alpha[32]-0.2999999999999999*(alpha[26]+alpha[25])+0.2999999999999999*(alpha[22]+alpha[21])-0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]+0.6037383539249431*alpha[16]-0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[9]+0.45*(alpha[8]+alpha[7])-0.45*(alpha[6]+alpha[5])-0.3354101966249685*alpha[4]+0.3354101966249685*(alpha[3]+alpha[2])-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[24] = ser_5x_p2_surfx3_quad_24_r(fEdge); 
  } else { 
    fUpwindQuad[24] = ser_5x_p2_surfx3_quad_24_l(fSkin); 
  } 
  if (0.5031152949374527*alpha[35]-0.5031152949374527*alpha[32]-0.2999999999999999*alpha[26]+0.375*alpha[25]+0.2999999999999999*alpha[22]-0.375*(alpha[21]+alpha[19])+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]-0.45*alpha[9]+0.45*alpha[7]-0.3354101966249685*alpha[4]+0.3354101966249685*(alpha[3]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad[25] = ser_5x_p2_surfx3_quad_25_r(fEdge); 
  } else { 
    fUpwindQuad[25] = ser_5x_p2_surfx3_quad_25_l(fSkin); 
  } 
  if ((-0.4024922359499621*(alpha[36]+alpha[35]))+0.4024922359499621*(alpha[33]+alpha[32])-0.2999999999999999*(alpha[26]+alpha[25])+0.2999999999999999*(alpha[22]+alpha[21])+0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]-0.6037383539249431*alpha[16]+0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])-0.45*(alpha[9]+alpha[8])+0.45*(alpha[7]+alpha[6]+alpha[5])-0.3354101966249685*alpha[4]+0.3354101966249685*(alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[26] = ser_5x_p2_surfx3_quad_26_r(fEdge); 
  } else { 
    fUpwindQuad[26] = ser_5x_p2_surfx3_quad_26_l(fSkin); 
  } 
  if (0.4024922359499621*(alpha[33]+alpha[32])-0.2999999999999999*(alpha[22]+alpha[21])-0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]-0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])+0.45*(alpha[7]+alpha[6]+alpha[5])-0.3354101966249685*(alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[27] = ser_5x_p2_surfx3_quad_27_r(fEdge); 
  } else { 
    fUpwindQuad[27] = ser_5x_p2_surfx3_quad_27_l(fSkin); 
  } 
  if ((-0.5031152949374527*alpha[32])-0.2999999999999999*alpha[22]+0.375*(alpha[21]+alpha[19])+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]+0.45*alpha[7]-0.3354101966249685*(alpha[3]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad[28] = ser_5x_p2_surfx3_quad_28_r(fEdge); 
  } else { 
    fUpwindQuad[28] = ser_5x_p2_surfx3_quad_28_l(fSkin); 
  } 
  if ((-0.4024922359499621*alpha[33])+0.4024922359499621*alpha[32]-0.2999999999999999*(alpha[22]+alpha[21])+0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]+0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[7]-0.45*(alpha[6]+alpha[5])-0.3354101966249685*(alpha[3]+alpha[2])+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[29] = ser_5x_p2_surfx3_quad_29_r(fEdge); 
  } else { 
    fUpwindQuad[29] = ser_5x_p2_surfx3_quad_29_l(fSkin); 
  } 
  if ((-0.5031152949374527*alpha[33])+0.375*alpha[22]-0.2999999999999999*alpha[21]+0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]+0.45*alpha[6]-0.3354101966249685*(alpha[3]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[30] = ser_5x_p2_surfx3_quad_30_r(fEdge); 
  } else { 
    fUpwindQuad[30] = ser_5x_p2_surfx3_quad_30_l(fSkin); 
  } 
  if (0.375*(alpha[22]+alpha[21])-0.2795084971874737*(alpha[12]+alpha[11])-0.3354101966249685*alpha[3]+0.25*alpha[0] > 0) { 
    fUpwindQuad[31] = ser_5x_p2_surfx3_quad_31_r(fEdge); 
  } else { 
    fUpwindQuad[31] = ser_5x_p2_surfx3_quad_31_l(fSkin); 
  } 
  if (0.5031152949374527*alpha[33]+0.375*alpha[22]-0.2999999999999999*alpha[21]-0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]-0.45*alpha[6]-0.3354101966249685*alpha[3]+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[32] = ser_5x_p2_surfx3_quad_32_r(fEdge); 
  } else { 
    fUpwindQuad[32] = ser_5x_p2_surfx3_quad_32_l(fSkin); 
  } 
  if (0.4024922359499621*alpha[33]-0.4024922359499621*alpha[32]-0.2999999999999999*(alpha[22]+alpha[21])-0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]+0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[7]+0.45*alpha[6]-0.45*alpha[5]-0.3354101966249685*alpha[3]+0.3354101966249685*alpha[2]-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[33] = ser_5x_p2_surfx3_quad_33_r(fEdge); 
  } else { 
    fUpwindQuad[33] = ser_5x_p2_surfx3_quad_33_l(fSkin); 
  } 
  if (0.5031152949374527*alpha[32]-0.2999999999999999*alpha[22]+0.375*alpha[21]-0.375*alpha[19]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]-0.45*alpha[7]-0.3354101966249685*alpha[3]+0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[34] = ser_5x_p2_surfx3_quad_34_r(fEdge); 
  } else { 
    fUpwindQuad[34] = ser_5x_p2_surfx3_quad_34_l(fSkin); 
  } 
  if ((-0.4024922359499621*(alpha[33]+alpha[32]))-0.2999999999999999*(alpha[22]+alpha[21])+0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]-0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])-0.45*(alpha[7]+alpha[6])+0.45*alpha[5]-0.3354101966249685*alpha[3]+0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[35] = ser_5x_p2_surfx3_quad_35_r(fEdge); 
  } else { 
    fUpwindQuad[35] = ser_5x_p2_surfx3_quad_35_l(fSkin); 
  } 
  if ((-0.2999999999999998*alpha[20])-0.2999999999999999*alpha[19]+0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[5]-0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[36] = ser_5x_p2_surfx3_quad_36_r(fEdge); 
  } else { 
    fUpwindQuad[36] = ser_5x_p2_surfx3_quad_36_l(fSkin); 
  } 
  if (0.375*alpha[19]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]-0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[37] = ser_5x_p2_surfx3_quad_37_r(fEdge); 
  } else { 
    fUpwindQuad[37] = ser_5x_p2_surfx3_quad_37_l(fSkin); 
  } 
  if (0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]+0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[5]-0.3354101966249685*alpha[2]+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[38] = ser_5x_p2_surfx3_quad_38_r(fEdge); 
  } else { 
    fUpwindQuad[38] = ser_5x_p2_surfx3_quad_38_l(fSkin); 
  } 
  if (0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[39] = ser_5x_p2_surfx3_quad_39_r(fEdge); 
  } else { 
    fUpwindQuad[39] = ser_5x_p2_surfx3_quad_39_l(fSkin); 
  } 
  if (0.25*alpha[0]-0.2795084971874737*(alpha[12]+alpha[11]) > 0) { 
    fUpwindQuad[40] = ser_5x_p2_surfx3_quad_40_r(fEdge); 
  } else { 
    fUpwindQuad[40] = ser_5x_p2_surfx3_quad_40_l(fSkin); 
  } 
  if ((-0.375*alpha[20])-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[41] = ser_5x_p2_surfx3_quad_41_r(fEdge); 
  } else { 
    fUpwindQuad[41] = ser_5x_p2_surfx3_quad_41_l(fSkin); 
  } 
  if ((-0.2999999999999998*alpha[20])+0.2999999999999999*alpha[19]+0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[5]+0.3354101966249685*alpha[2]-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[42] = ser_5x_p2_surfx3_quad_42_r(fEdge); 
  } else { 
    fUpwindQuad[42] = ser_5x_p2_surfx3_quad_42_l(fSkin); 
  } 
  if ((-0.375*alpha[19])+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]+0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[43] = ser_5x_p2_surfx3_quad_43_r(fEdge); 
  } else { 
    fUpwindQuad[43] = ser_5x_p2_surfx3_quad_43_l(fSkin); 
  } 
  if (0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]+0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[5]+0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[44] = ser_5x_p2_surfx3_quad_44_r(fEdge); 
  } else { 
    fUpwindQuad[44] = ser_5x_p2_surfx3_quad_44_l(fSkin); 
  } 
  if ((-0.4024922359499621*(alpha[33]+alpha[32]))+0.2999999999999999*(alpha[22]+alpha[21])-0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]+0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])-0.45*(alpha[7]+alpha[6])+0.45*alpha[5]+0.3354101966249685*alpha[3]-0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[45] = ser_5x_p2_surfx3_quad_45_r(fEdge); 
  } else { 
    fUpwindQuad[45] = ser_5x_p2_surfx3_quad_45_l(fSkin); 
  } 
  if (0.5031152949374527*alpha[32]+0.2999999999999999*alpha[22]-0.375*alpha[21]+0.375*alpha[19]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]-0.45*alpha[7]+0.3354101966249685*alpha[3]-0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[46] = ser_5x_p2_surfx3_quad_46_r(fEdge); 
  } else { 
    fUpwindQuad[46] = ser_5x_p2_surfx3_quad_46_l(fSkin); 
  } 
  if (0.4024922359499621*alpha[33]-0.4024922359499621*alpha[32]+0.2999999999999999*(alpha[22]+alpha[21])+0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]-0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[7]+0.45*alpha[6]-0.45*alpha[5]+0.3354101966249685*alpha[3]-0.3354101966249685*alpha[2]+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[47] = ser_5x_p2_surfx3_quad_47_r(fEdge); 
  } else { 
    fUpwindQuad[47] = ser_5x_p2_surfx3_quad_47_l(fSkin); 
  } 
  if (0.5031152949374527*alpha[33]-0.375*alpha[22]+0.2999999999999999*alpha[21]+0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]-0.45*alpha[6]+0.3354101966249685*alpha[3]-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[48] = ser_5x_p2_surfx3_quad_48_r(fEdge); 
  } else { 
    fUpwindQuad[48] = ser_5x_p2_surfx3_quad_48_l(fSkin); 
  } 
  if ((-0.375*(alpha[22]+alpha[21]))-0.2795084971874737*(alpha[12]+alpha[11])+0.3354101966249685*alpha[3]+0.25*alpha[0] > 0) { 
    fUpwindQuad[49] = ser_5x_p2_surfx3_quad_49_r(fEdge); 
  } else { 
    fUpwindQuad[49] = ser_5x_p2_surfx3_quad_49_l(fSkin); 
  } 
  if ((-0.5031152949374527*alpha[33])-0.375*alpha[22]+0.2999999999999999*alpha[21]-0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]+0.45*alpha[6]+0.3354101966249685*(alpha[3]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[50] = ser_5x_p2_surfx3_quad_50_r(fEdge); 
  } else { 
    fUpwindQuad[50] = ser_5x_p2_surfx3_quad_50_l(fSkin); 
  } 
  if ((-0.4024922359499621*alpha[33])+0.4024922359499621*alpha[32]+0.2999999999999999*(alpha[22]+alpha[21])-0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]-0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[7]-0.45*(alpha[6]+alpha[5])+0.3354101966249685*(alpha[3]+alpha[2])-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[51] = ser_5x_p2_surfx3_quad_51_r(fEdge); 
  } else { 
    fUpwindQuad[51] = ser_5x_p2_surfx3_quad_51_l(fSkin); 
  } 
  if ((-0.5031152949374527*alpha[32])+0.2999999999999999*alpha[22]-0.375*(alpha[21]+alpha[19])+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]+0.45*alpha[7]+0.3354101966249685*(alpha[3]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad[52] = ser_5x_p2_surfx3_quad_52_r(fEdge); 
  } else { 
    fUpwindQuad[52] = ser_5x_p2_surfx3_quad_52_l(fSkin); 
  } 
  if (0.4024922359499621*(alpha[33]+alpha[32])+0.2999999999999999*(alpha[22]+alpha[21])+0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]+0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])+0.45*(alpha[7]+alpha[6]+alpha[5])+0.3354101966249685*(alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[53] = ser_5x_p2_surfx3_quad_53_r(fEdge); 
  } else { 
    fUpwindQuad[53] = ser_5x_p2_surfx3_quad_53_l(fSkin); 
  } 
  if ((-0.4024922359499621*(alpha[36]+alpha[35]))+0.4024922359499621*(alpha[33]+alpha[32])+0.2999999999999999*(alpha[26]+alpha[25])-0.2999999999999999*(alpha[22]+alpha[21])-0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]+0.6037383539249431*alpha[16]-0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])-0.45*(alpha[9]+alpha[8])+0.45*(alpha[7]+alpha[6]+alpha[5])+0.3354101966249685*alpha[4]-0.3354101966249685*(alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[54] = ser_5x_p2_surfx3_quad_54_r(fEdge); 
  } else { 
    fUpwindQuad[54] = ser_5x_p2_surfx3_quad_54_l(fSkin); 
  } 
  if (0.5031152949374527*alpha[35]-0.5031152949374527*alpha[32]+0.2999999999999999*alpha[26]-0.375*alpha[25]-0.2999999999999999*alpha[22]+0.375*(alpha[21]+alpha[19])+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]-0.45*alpha[9]+0.45*alpha[7]+0.3354101966249685*alpha[4]-0.3354101966249685*(alpha[3]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad[55] = ser_5x_p2_surfx3_quad_55_r(fEdge); 
  } else { 
    fUpwindQuad[55] = ser_5x_p2_surfx3_quad_55_l(fSkin); 
  } 
  if (0.4024922359499621*alpha[36]-0.4024922359499621*(alpha[35]+alpha[33])+0.4024922359499621*alpha[32]+0.2999999999999999*(alpha[26]+alpha[25])-0.2999999999999999*(alpha[22]+alpha[21])+0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]-0.6037383539249431*alpha[16]+0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[9]+0.45*(alpha[8]+alpha[7])-0.45*(alpha[6]+alpha[5])+0.3354101966249685*alpha[4]-0.3354101966249685*(alpha[3]+alpha[2])+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[56] = ser_5x_p2_surfx3_quad_56_r(fEdge); 
  } else { 
    fUpwindQuad[56] = ser_5x_p2_surfx3_quad_56_l(fSkin); 
  } 
  if (0.5031152949374527*alpha[36]-0.5031152949374527*alpha[33]-0.375*alpha[26]+0.2999999999999999*alpha[25]+0.375*alpha[22]-0.2999999999999999*alpha[21]+0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]-0.45*alpha[8]+0.45*alpha[6]+0.3354101966249685*alpha[4]-0.3354101966249685*(alpha[3]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[57] = ser_5x_p2_surfx3_quad_57_r(fEdge); 
  } else { 
    fUpwindQuad[57] = ser_5x_p2_surfx3_quad_57_l(fSkin); 
  } 
  if ((-0.375*(alpha[26]+alpha[25]))+0.375*(alpha[22]+alpha[21])-0.2795084971874737*(alpha[12]+alpha[11])+0.3354101966249685*alpha[4]-0.3354101966249685*alpha[3]+0.25*alpha[0] > 0) { 
    fUpwindQuad[58] = ser_5x_p2_surfx3_quad_58_r(fEdge); 
  } else { 
    fUpwindQuad[58] = ser_5x_p2_surfx3_quad_58_l(fSkin); 
  } 
  if ((-0.5031152949374527*alpha[36])+0.5031152949374527*alpha[33]-0.375*alpha[26]+0.2999999999999999*alpha[25]+0.375*alpha[22]-0.2999999999999999*alpha[21]-0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]+0.45*alpha[8]-0.45*alpha[6]+0.3354101966249685*alpha[4]-0.3354101966249685*alpha[3]+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[59] = ser_5x_p2_surfx3_quad_59_r(fEdge); 
  } else { 
    fUpwindQuad[59] = ser_5x_p2_surfx3_quad_59_l(fSkin); 
  } 
  if ((-0.4024922359499621*alpha[36])+0.4024922359499621*(alpha[35]+alpha[33])-0.4024922359499621*alpha[32]+0.2999999999999999*(alpha[26]+alpha[25])-0.2999999999999999*(alpha[22]+alpha[21])-0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]-0.6037383539249431*alpha[16]+0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[9]-0.45*(alpha[8]+alpha[7])+0.45*alpha[6]-0.45*alpha[5]+0.3354101966249685*alpha[4]-0.3354101966249685*alpha[3]+0.3354101966249685*alpha[2]-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[60] = ser_5x_p2_surfx3_quad_60_r(fEdge); 
  } else { 
    fUpwindQuad[60] = ser_5x_p2_surfx3_quad_60_l(fSkin); 
  } 
  if ((-0.5031152949374527*alpha[35])+0.5031152949374527*alpha[32]+0.2999999999999999*alpha[26]-0.375*alpha[25]-0.2999999999999999*alpha[22]+0.375*alpha[21]-0.375*alpha[19]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]+0.45*alpha[9]-0.45*alpha[7]+0.3354101966249685*alpha[4]-0.3354101966249685*alpha[3]+0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[61] = ser_5x_p2_surfx3_quad_61_r(fEdge); 
  } else { 
    fUpwindQuad[61] = ser_5x_p2_surfx3_quad_61_l(fSkin); 
  } 
  if (0.4024922359499621*(alpha[36]+alpha[35])-0.4024922359499621*(alpha[33]+alpha[32])+0.2999999999999999*(alpha[26]+alpha[25])-0.2999999999999999*(alpha[22]+alpha[21])+0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]+0.6037383539249431*alpha[16]-0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])+0.45*(alpha[9]+alpha[8])-0.45*(alpha[7]+alpha[6])+0.45*alpha[5]+0.3354101966249685*alpha[4]-0.3354101966249685*alpha[3]+0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[62] = ser_5x_p2_surfx3_quad_62_r(fEdge); 
  } else { 
    fUpwindQuad[62] = ser_5x_p2_surfx3_quad_62_l(fSkin); 
  } 
  if ((-0.4024922359499621*(alpha[36]+alpha[35]))+0.2999999999999999*(alpha[26]+alpha[25])-0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]+0.6037383539249431*alpha[16]+0.223606797749979*(alpha[12]+alpha[11])-0.45*(alpha[9]+alpha[8])+0.45*alpha[5]+0.3354101966249685*alpha[4]-0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[63] = ser_5x_p2_surfx3_quad_63_r(fEdge); 
  } else { 
    fUpwindQuad[63] = ser_5x_p2_surfx3_quad_63_l(fSkin); 
  } 
  if (0.5031152949374527*alpha[35]+0.2999999999999999*alpha[26]-0.375*alpha[25]+0.375*alpha[19]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]-0.45*alpha[9]+0.3354101966249685*alpha[4]-0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[64] = ser_5x_p2_surfx3_quad_64_r(fEdge); 
  } else { 
    fUpwindQuad[64] = ser_5x_p2_surfx3_quad_64_l(fSkin); 
  } 
  if (0.4024922359499621*alpha[36]-0.4024922359499621*alpha[35]+0.2999999999999999*(alpha[26]+alpha[25])+0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]-0.6037383539249431*alpha[16]+0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[9]+0.45*alpha[8]-0.45*alpha[5]+0.3354101966249685*alpha[4]-0.3354101966249685*alpha[2]+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[65] = ser_5x_p2_surfx3_quad_65_r(fEdge); 
  } else { 
    fUpwindQuad[65] = ser_5x_p2_surfx3_quad_65_l(fSkin); 
  } 
  if (0.5031152949374527*alpha[36]-0.375*alpha[26]+0.2999999999999999*alpha[25]+0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]-0.45*alpha[8]+0.3354101966249685*alpha[4]-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[66] = ser_5x_p2_surfx3_quad_66_r(fEdge); 
  } else { 
    fUpwindQuad[66] = ser_5x_p2_surfx3_quad_66_l(fSkin); 
  } 
  if ((-0.375*(alpha[26]+alpha[25]))-0.2795084971874737*(alpha[12]+alpha[11])+0.3354101966249685*alpha[4]+0.25*alpha[0] > 0) { 
    fUpwindQuad[67] = ser_5x_p2_surfx3_quad_67_r(fEdge); 
  } else { 
    fUpwindQuad[67] = ser_5x_p2_surfx3_quad_67_l(fSkin); 
  } 
  if ((-0.5031152949374527*alpha[36])-0.375*alpha[26]+0.2999999999999999*alpha[25]-0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]+0.45*alpha[8]+0.3354101966249685*(alpha[4]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[68] = ser_5x_p2_surfx3_quad_68_r(fEdge); 
  } else { 
    fUpwindQuad[68] = ser_5x_p2_surfx3_quad_68_l(fSkin); 
  } 
  if ((-0.4024922359499621*alpha[36])+0.4024922359499621*alpha[35]+0.2999999999999999*(alpha[26]+alpha[25])-0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]-0.6037383539249431*alpha[16]+0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[9]-0.45*(alpha[8]+alpha[5])+0.3354101966249685*(alpha[4]+alpha[2])-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[69] = ser_5x_p2_surfx3_quad_69_r(fEdge); 
  } else { 
    fUpwindQuad[69] = ser_5x_p2_surfx3_quad_69_l(fSkin); 
  } 
  if ((-0.5031152949374527*alpha[35])+0.2999999999999999*alpha[26]-0.375*(alpha[25]+alpha[19])+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]+0.45*alpha[9]+0.3354101966249685*(alpha[4]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad[70] = ser_5x_p2_surfx3_quad_70_r(fEdge); 
  } else { 
    fUpwindQuad[70] = ser_5x_p2_surfx3_quad_70_l(fSkin); 
  } 
  if (0.4024922359499621*(alpha[36]+alpha[35])+0.2999999999999999*(alpha[26]+alpha[25])+0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]+0.6037383539249431*alpha[16]+0.223606797749979*(alpha[12]+alpha[11])+0.45*(alpha[9]+alpha[8]+alpha[5])+0.3354101966249685*(alpha[4]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[71] = ser_5x_p2_surfx3_quad_71_r(fEdge); 
  } else { 
    fUpwindQuad[71] = ser_5x_p2_surfx3_quad_71_l(fSkin); 
  } 
  if ((-0.4024922359499621*(alpha[36]+alpha[35]+alpha[33]+alpha[32]))+0.2999999999999999*(alpha[26]+alpha[25]+alpha[22]+alpha[21])-0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]+0.6037383539249431*(alpha[16]+alpha[15])+0.223606797749979*(alpha[12]+alpha[11])-0.45*(alpha[9]+alpha[8]+alpha[7]+alpha[6])+0.45*alpha[5]+0.3354101966249685*(alpha[4]+alpha[3])-0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[72] = ser_5x_p2_surfx3_quad_72_r(fEdge); 
  } else { 
    fUpwindQuad[72] = ser_5x_p2_surfx3_quad_72_l(fSkin); 
  } 
  if (0.5031152949374527*(alpha[35]+alpha[32])+0.2999999999999999*alpha[26]-0.375*alpha[25]+0.2999999999999999*alpha[22]-0.375*alpha[21]+0.375*alpha[19]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]-0.45*(alpha[9]+alpha[7])+0.3354101966249685*(alpha[4]+alpha[3])-0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[73] = ser_5x_p2_surfx3_quad_73_r(fEdge); 
  } else { 
    fUpwindQuad[73] = ser_5x_p2_surfx3_quad_73_l(fSkin); 
  } 
  if (0.4024922359499621*alpha[36]-0.4024922359499621*alpha[35]+0.4024922359499621*alpha[33]-0.4024922359499621*alpha[32]+0.2999999999999999*(alpha[26]+alpha[25]+alpha[22]+alpha[21])+0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]-0.6037383539249431*(alpha[16]+alpha[15])+0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[9]+0.45*alpha[8]-0.45*alpha[7]+0.45*alpha[6]-0.45*alpha[5]+0.3354101966249685*(alpha[4]+alpha[3])-0.3354101966249685*alpha[2]+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[74] = ser_5x_p2_surfx3_quad_74_r(fEdge); 
  } else { 
    fUpwindQuad[74] = ser_5x_p2_surfx3_quad_74_l(fSkin); 
  } 
  if (0.5031152949374527*(alpha[36]+alpha[33])-0.375*alpha[26]+0.2999999999999999*alpha[25]-0.375*alpha[22]+0.2999999999999999*alpha[21]+0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]-0.45*(alpha[8]+alpha[6])+0.3354101966249685*(alpha[4]+alpha[3])-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[75] = ser_5x_p2_surfx3_quad_75_r(fEdge); 
  } else { 
    fUpwindQuad[75] = ser_5x_p2_surfx3_quad_75_l(fSkin); 
  } 
  if ((-0.375*(alpha[26]+alpha[25]+alpha[22]+alpha[21]))-0.2795084971874737*(alpha[12]+alpha[11])+0.3354101966249685*(alpha[4]+alpha[3])+0.25*alpha[0] > 0) { 
    fUpwindQuad[76] = ser_5x_p2_surfx3_quad_76_r(fEdge); 
  } else { 
    fUpwindQuad[76] = ser_5x_p2_surfx3_quad_76_l(fSkin); 
  } 
  if ((-0.5031152949374527*(alpha[36]+alpha[33]))-0.375*alpha[26]+0.2999999999999999*alpha[25]-0.375*alpha[22]+0.2999999999999999*alpha[21]-0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]+0.45*(alpha[8]+alpha[6])+0.3354101966249685*(alpha[4]+alpha[3]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[77] = ser_5x_p2_surfx3_quad_77_r(fEdge); 
  } else { 
    fUpwindQuad[77] = ser_5x_p2_surfx3_quad_77_l(fSkin); 
  } 
  if ((-0.4024922359499621*alpha[36])+0.4024922359499621*alpha[35]-0.4024922359499621*alpha[33]+0.4024922359499621*alpha[32]+0.2999999999999999*(alpha[26]+alpha[25]+alpha[22]+alpha[21])-0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]-0.6037383539249431*(alpha[16]+alpha[15])+0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[9]-0.45*alpha[8]+0.45*alpha[7]-0.45*(alpha[6]+alpha[5])+0.3354101966249685*(alpha[4]+alpha[3]+alpha[2])-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[78] = ser_5x_p2_surfx3_quad_78_r(fEdge); 
  } else { 
    fUpwindQuad[78] = ser_5x_p2_surfx3_quad_78_l(fSkin); 
  } 
  if ((-0.5031152949374527*(alpha[35]+alpha[32]))+0.2999999999999999*alpha[26]-0.375*alpha[25]+0.2999999999999999*alpha[22]-0.375*(alpha[21]+alpha[19])+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]+0.45*(alpha[9]+alpha[7])+0.3354101966249685*(alpha[4]+alpha[3]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad[79] = ser_5x_p2_surfx3_quad_79_r(fEdge); 
  } else { 
    fUpwindQuad[79] = ser_5x_p2_surfx3_quad_79_l(fSkin); 
  } 
  if (0.4024922359499621*(alpha[36]+alpha[35]+alpha[33]+alpha[32])+0.2999999999999999*(alpha[26]+alpha[25]+alpha[22]+alpha[21])+0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]+0.6037383539249431*(alpha[16]+alpha[15])+0.223606797749979*(alpha[12]+alpha[11])+0.45*(alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5])+0.3354101966249685*(alpha[4]+alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[80] = ser_5x_p2_surfx3_quad_80_r(fEdge); 
  } else { 
    fUpwindQuad[80] = ser_5x_p2_surfx3_quad_80_l(fSkin); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_5x_p2_upwind(fUpwindQuad, fUpwind); 

  Ghat[0] += 0.25*(alpha[36]*fUpwind[36]+alpha[35]*fUpwind[35]+alpha[33]*fUpwind[33]+alpha[32]*fUpwind[32]+alpha[26]*fUpwind[26]+alpha[25]*fUpwind[25]+alpha[22]*fUpwind[22]+alpha[21]*fUpwind[21]+alpha[20]*fUpwind[20]+alpha[19]*fUpwind[19]+alpha[16]*fUpwind[16]+alpha[15]*fUpwind[15]+alpha[12]*fUpwind[12]+alpha[11]*fUpwind[11]+alpha[9]*fUpwind[9]+alpha[8]*fUpwind[8]+alpha[7]*fUpwind[7]+alpha[6]*fUpwind[6]+alpha[5]*fUpwind[5]+alpha[4]*fUpwind[4]+alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] += 0.2500000000000001*(alpha[26]*fUpwind[36]+fUpwind[26]*alpha[36])+0.223606797749979*(alpha[16]*fUpwind[35]+fUpwind[16]*alpha[35])+0.2500000000000001*(alpha[22]*fUpwind[33]+fUpwind[22]*alpha[33])+0.223606797749979*(alpha[15]*fUpwind[32]+fUpwind[15]*alpha[32])+0.223606797749979*(alpha[8]*fUpwind[25]+fUpwind[8]*alpha[25]+alpha[6]*fUpwind[21]+fUpwind[6]*alpha[21])+0.2500000000000001*(alpha[12]*fUpwind[20]+fUpwind[12]*alpha[20])+0.223606797749979*(alpha[5]*fUpwind[19]+fUpwind[5]*alpha[19])+0.25*(alpha[9]*fUpwind[16]+fUpwind[9]*alpha[16]+alpha[7]*fUpwind[15]+fUpwind[7]*alpha[15])+0.223606797749979*(alpha[1]*fUpwind[11]+fUpwind[1]*alpha[11])+0.25*(alpha[4]*fUpwind[8]+fUpwind[4]*alpha[8]+alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6]+alpha[2]*fUpwind[5]+fUpwind[2]*alpha[5]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] += 0.223606797749979*(alpha[16]*fUpwind[36]+fUpwind[16]*alpha[36])+0.2500000000000001*(alpha[25]*fUpwind[35]+fUpwind[25]*alpha[35])+0.223606797749979*(alpha[15]*fUpwind[33]+fUpwind[15]*alpha[33])+0.2500000000000001*(alpha[21]*fUpwind[32]+fUpwind[21]*alpha[32])+0.223606797749979*(alpha[9]*fUpwind[26]+fUpwind[9]*alpha[26]+alpha[7]*fUpwind[22]+fUpwind[7]*alpha[22]+alpha[5]*fUpwind[20]+fUpwind[5]*alpha[20])+0.2500000000000001*(alpha[11]*fUpwind[19]+fUpwind[11]*alpha[19])+0.25*(alpha[8]*fUpwind[16]+fUpwind[8]*alpha[16]+alpha[6]*fUpwind[15]+fUpwind[6]*alpha[15])+0.223606797749979*(alpha[2]*fUpwind[12]+fUpwind[2]*alpha[12])+0.25*(alpha[4]*fUpwind[9]+fUpwind[4]*alpha[9]+alpha[3]*fUpwind[7]+fUpwind[3]*alpha[7]+alpha[1]*fUpwind[5]+fUpwind[1]*alpha[5]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] += 0.2500000000000001*(alpha[36]*fUpwind[45]+alpha[35]*fUpwind[44]+alpha[26]*fUpwind[38]+alpha[25]*fUpwind[37])+0.223606797749979*alpha[15]*fUpwind[34]+0.2500000000000001*(alpha[20]*fUpwind[33]+fUpwind[20]*alpha[33]+alpha[19]*fUpwind[32]+fUpwind[19]*alpha[32])+0.25*alpha[16]*fUpwind[31]+0.223606797749979*(alpha[7]*fUpwind[24]+alpha[6]*fUpwind[23])+0.2500000000000001*(alpha[12]*fUpwind[22]+fUpwind[12]*alpha[22]+alpha[11]*fUpwind[21]+fUpwind[11]*alpha[21])+0.25*(alpha[9]*fUpwind[18]+alpha[8]*fUpwind[17]+alpha[5]*fUpwind[15]+fUpwind[5]*alpha[15])+0.223606797749979*alpha[3]*fUpwind[13]+0.25*(alpha[4]*fUpwind[10]+alpha[2]*fUpwind[7]+fUpwind[2]*alpha[7]+alpha[1]*fUpwind[6]+fUpwind[1]*alpha[6]+alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]); 
  Ghat[4] += 0.2500000000000001*(alpha[33]*fUpwind[45]+alpha[32]*fUpwind[44])+0.223606797749979*alpha[16]*fUpwind[41]+0.2500000000000001*(alpha[22]*fUpwind[38]+alpha[21]*fUpwind[37]+alpha[20]*fUpwind[36]+fUpwind[20]*alpha[36]+alpha[19]*fUpwind[35]+fUpwind[19]*alpha[35])+0.25*alpha[15]*fUpwind[31]+0.223606797749979*(alpha[9]*fUpwind[29]+alpha[8]*fUpwind[28])+0.2500000000000001*(alpha[12]*fUpwind[26]+fUpwind[12]*alpha[26]+alpha[11]*fUpwind[25]+fUpwind[11]*alpha[25])+0.25*(alpha[7]*fUpwind[18]+alpha[6]*fUpwind[17]+alpha[5]*fUpwind[16]+fUpwind[5]*alpha[16])+0.223606797749979*alpha[4]*fUpwind[14]+0.25*(alpha[3]*fUpwind[10]+alpha[2]*fUpwind[9]+fUpwind[2]*alpha[9]+alpha[1]*fUpwind[8]+fUpwind[1]*alpha[8]+alpha[0]*fUpwind[4]+fUpwind[0]*alpha[4]); 
  Ghat[5] += (0.2*alpha[35]+0.223606797749979*alpha[9])*fUpwind[36]+0.2*fUpwind[35]*alpha[36]+0.223606797749979*(fUpwind[9]*alpha[36]+alpha[8]*fUpwind[35]+fUpwind[8]*alpha[35])+(0.2*alpha[32]+0.223606797749979*alpha[7])*fUpwind[33]+0.2*fUpwind[32]*alpha[33]+0.223606797749979*(fUpwind[7]*alpha[33]+alpha[6]*fUpwind[32]+fUpwind[6]*alpha[32])+0.223606797749979*(alpha[16]*fUpwind[26]+fUpwind[16]*alpha[26]+alpha[16]*fUpwind[25]+fUpwind[16]*alpha[25]+alpha[15]*fUpwind[22]+fUpwind[15]*alpha[22]+alpha[15]*fUpwind[21]+fUpwind[15]*alpha[21])+(0.2*alpha[19]+0.223606797749979*alpha[2])*fUpwind[20]+0.2*fUpwind[19]*alpha[20]+0.223606797749979*(fUpwind[2]*alpha[20]+alpha[1]*fUpwind[19]+fUpwind[1]*alpha[19])+0.25*(alpha[4]*fUpwind[16]+fUpwind[4]*alpha[16]+alpha[3]*fUpwind[15]+fUpwind[3]*alpha[15])+0.223606797749979*(alpha[5]*fUpwind[12]+fUpwind[5]*alpha[12]+alpha[5]*fUpwind[11]+fUpwind[5]*alpha[11])+0.25*(alpha[8]*fUpwind[9]+fUpwind[8]*alpha[9]+alpha[6]*fUpwind[7]+fUpwind[6]*alpha[7]+alpha[0]*fUpwind[5]+fUpwind[0]*alpha[5]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[6] += 0.25*alpha[26]*fUpwind[45]+0.223606797749979*alpha[16]*fUpwind[44]+0.25*alpha[36]*fUpwind[38]+0.223606797749979*(alpha[8]*fUpwind[37]+fUpwind[31]*alpha[35])+(0.2*alpha[32]+0.223606797749979*alpha[7])*fUpwind[34]+0.25*(alpha[12]*fUpwind[33]+fUpwind[12]*alpha[33])+0.223606797749979*(alpha[5]*fUpwind[32]+fUpwind[5]*alpha[32])+0.25*alpha[9]*fUpwind[31]+0.223606797749979*(fUpwind[17]*alpha[25]+alpha[15]*fUpwind[24])+(0.2*alpha[21]+0.223606797749979*alpha[3])*fUpwind[23]+0.25*(alpha[20]*fUpwind[22]+fUpwind[20]*alpha[22])+0.223606797749979*(alpha[1]*fUpwind[21]+fUpwind[1]*alpha[21]+alpha[15]*fUpwind[19]+fUpwind[15]*alpha[19])+0.25*(alpha[16]*fUpwind[18]+alpha[4]*fUpwind[17]+alpha[2]*fUpwind[15]+fUpwind[2]*alpha[15])+0.223606797749979*(alpha[6]*(fUpwind[13]+fUpwind[11])+fUpwind[6]*alpha[11])+0.25*(alpha[8]*fUpwind[10]+alpha[5]*fUpwind[7]+fUpwind[5]*alpha[7]+alpha[0]*fUpwind[6]+fUpwind[0]*alpha[6]+alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]); 
  Ghat[7] += 0.223606797749979*alpha[16]*fUpwind[45]+0.25*alpha[25]*fUpwind[44]+0.223606797749979*alpha[9]*fUpwind[38]+0.25*alpha[35]*fUpwind[37]+0.223606797749979*fUpwind[31]*alpha[36]+0.2*alpha[33]*fUpwind[34]+0.223606797749979*(alpha[6]*fUpwind[34]+alpha[5]*fUpwind[33]+fUpwind[5]*alpha[33])+0.25*(alpha[11]*fUpwind[32]+fUpwind[11]*alpha[32]+alpha[8]*fUpwind[31])+0.223606797749979*fUpwind[18]*alpha[26]+0.2*alpha[22]*fUpwind[24]+0.223606797749979*(alpha[3]*fUpwind[24]+alpha[15]*fUpwind[23]+alpha[2]*fUpwind[22]+fUpwind[2]*alpha[22])+0.25*(alpha[19]*fUpwind[21]+fUpwind[19]*alpha[21])+0.223606797749979*(alpha[15]*fUpwind[20]+fUpwind[15]*alpha[20])+0.25*(alpha[4]*fUpwind[18]+alpha[16]*fUpwind[17]+alpha[1]*fUpwind[15]+fUpwind[1]*alpha[15])+0.223606797749979*(alpha[7]*(fUpwind[13]+fUpwind[12])+fUpwind[7]*alpha[12])+0.25*(alpha[9]*fUpwind[10]+alpha[0]*fUpwind[7]+fUpwind[0]*alpha[7]+alpha[5]*fUpwind[6]+fUpwind[5]*alpha[6]+alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]); 
  Ghat[8] += 0.25*alpha[22]*fUpwind[45]+0.223606797749979*alpha[15]*fUpwind[44]+(0.2*alpha[35]+0.223606797749979*alpha[9])*fUpwind[41]+0.25*alpha[33]*fUpwind[38]+0.223606797749979*alpha[6]*fUpwind[37]+0.25*(alpha[12]*fUpwind[36]+fUpwind[12]*alpha[36])+0.223606797749979*(alpha[5]*fUpwind[35]+fUpwind[5]*alpha[35])+fUpwind[31]*(0.223606797749979*alpha[32]+0.25*alpha[7])+0.223606797749979*alpha[16]*fUpwind[29]+(0.2*alpha[25]+0.223606797749979*alpha[4])*fUpwind[28]+0.25*(alpha[20]*fUpwind[26]+fUpwind[20]*alpha[26])+0.223606797749979*(alpha[1]*fUpwind[25]+fUpwind[1]*alpha[25]+fUpwind[17]*alpha[21]+alpha[16]*fUpwind[19]+fUpwind[16]*alpha[19])+0.25*(alpha[15]*fUpwind[18]+alpha[3]*fUpwind[17]+alpha[2]*fUpwind[16]+fUpwind[2]*alpha[16])+0.223606797749979*(alpha[8]*(fUpwind[14]+fUpwind[11])+fUpwind[8]*alpha[11])+0.25*(alpha[6]*fUpwind[10]+alpha[5]*fUpwind[9]+fUpwind[5]*alpha[9]+alpha[0]*fUpwind[8]+fUpwind[0]*alpha[8]+alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]); 
  Ghat[9] += 0.223606797749979*alpha[15]*fUpwind[45]+0.25*alpha[21]*fUpwind[44]+0.2*alpha[36]*fUpwind[41]+0.223606797749979*(alpha[8]*fUpwind[41]+alpha[7]*fUpwind[38])+0.25*alpha[32]*fUpwind[37]+0.223606797749979*(alpha[5]*fUpwind[36]+fUpwind[5]*alpha[36])+0.25*(alpha[11]*fUpwind[35]+fUpwind[11]*alpha[35])+fUpwind[31]*(0.223606797749979*alpha[33]+0.25*alpha[6])+0.2*alpha[26]*fUpwind[29]+0.223606797749979*(alpha[4]*fUpwind[29]+alpha[16]*fUpwind[28]+alpha[2]*fUpwind[26]+fUpwind[2]*alpha[26])+0.25*(alpha[19]*fUpwind[25]+fUpwind[19]*alpha[25])+0.223606797749979*(fUpwind[18]*alpha[22]+alpha[16]*fUpwind[20]+fUpwind[16]*alpha[20])+0.25*(alpha[3]*fUpwind[18]+alpha[15]*fUpwind[17]+alpha[1]*fUpwind[16]+fUpwind[1]*alpha[16])+0.223606797749979*(alpha[9]*(fUpwind[14]+fUpwind[12])+fUpwind[9]*alpha[12])+0.25*(alpha[7]*fUpwind[10]+alpha[0]*fUpwind[9]+fUpwind[0]*alpha[9]+alpha[5]*fUpwind[8]+fUpwind[5]*alpha[8]+alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]); 
  Ghat[10] += 0.223606797749979*(alpha[16]*fUpwind[47]+alpha[15]*fUpwind[46])+0.25*(alpha[20]*fUpwind[45]+alpha[19]*fUpwind[44])+0.223606797749979*(alpha[9]*fUpwind[43]+alpha[8]*fUpwind[42]+alpha[7]*fUpwind[40]+alpha[6]*fUpwind[39])+0.25*(alpha[12]*fUpwind[38]+alpha[11]*fUpwind[37]+alpha[33]*fUpwind[36]+fUpwind[33]*alpha[36]+alpha[32]*fUpwind[35]+fUpwind[32]*alpha[35]+alpha[5]*fUpwind[31])+0.223606797749979*(alpha[4]*fUpwind[30]+alpha[3]*fUpwind[27])+0.25*(alpha[22]*fUpwind[26]+fUpwind[22]*alpha[26]+alpha[21]*fUpwind[25]+fUpwind[21]*alpha[25]+alpha[2]*fUpwind[18]+alpha[1]*fUpwind[17]+alpha[15]*fUpwind[16]+fUpwind[15]*alpha[16]+alpha[0]*fUpwind[10]+alpha[7]*fUpwind[9]+fUpwind[7]*alpha[9]+alpha[6]*fUpwind[8]+fUpwind[6]*alpha[8]+alpha[3]*fUpwind[4]+fUpwind[3]*alpha[4]); 
  Ghat[11] += 0.223606797749979*alpha[36]*fUpwind[36]+0.159719141249985*alpha[35]*fUpwind[35]+0.25*(alpha[9]*fUpwind[35]+fUpwind[9]*alpha[35])+0.223606797749979*alpha[33]*fUpwind[33]+0.159719141249985*alpha[32]*fUpwind[32]+0.25*(alpha[7]*fUpwind[32]+fUpwind[7]*alpha[32])+0.159719141249985*alpha[25]*fUpwind[25]+0.2500000000000001*(alpha[4]*fUpwind[25]+fUpwind[4]*alpha[25])+0.159719141249985*alpha[21]*fUpwind[21]+0.2500000000000001*(alpha[3]*fUpwind[21]+fUpwind[3]*alpha[21])+0.223606797749979*alpha[20]*fUpwind[20]+0.159719141249985*alpha[19]*fUpwind[19]+0.2500000000000001*(alpha[2]*fUpwind[19]+fUpwind[2]*alpha[19])+0.223606797749979*(alpha[16]*fUpwind[16]+alpha[15]*fUpwind[15])+0.159719141249985*alpha[11]*fUpwind[11]+0.25*(alpha[0]*fUpwind[11]+fUpwind[0]*alpha[11])+0.223606797749979*(alpha[8]*fUpwind[8]+alpha[6]*fUpwind[6]+alpha[5]*fUpwind[5]+alpha[1]*fUpwind[1]); 
  Ghat[12] += 0.159719141249985*alpha[36]*fUpwind[36]+0.25*(alpha[8]*fUpwind[36]+fUpwind[8]*alpha[36])+0.223606797749979*alpha[35]*fUpwind[35]+0.159719141249985*alpha[33]*fUpwind[33]+0.25*(alpha[6]*fUpwind[33]+fUpwind[6]*alpha[33])+0.223606797749979*alpha[32]*fUpwind[32]+0.159719141249985*alpha[26]*fUpwind[26]+0.2500000000000001*(alpha[4]*fUpwind[26]+fUpwind[4]*alpha[26])+0.159719141249985*alpha[22]*fUpwind[22]+0.2500000000000001*(alpha[3]*fUpwind[22]+fUpwind[3]*alpha[22])+0.159719141249985*alpha[20]*fUpwind[20]+0.2500000000000001*(alpha[1]*fUpwind[20]+fUpwind[1]*alpha[20])+0.223606797749979*(alpha[19]*fUpwind[19]+alpha[16]*fUpwind[16]+alpha[15]*fUpwind[15])+0.159719141249985*alpha[12]*fUpwind[12]+0.25*(alpha[0]*fUpwind[12]+fUpwind[0]*alpha[12])+0.223606797749979*(alpha[9]*fUpwind[9]+alpha[7]*fUpwind[7]+alpha[5]*fUpwind[5]+alpha[2]*fUpwind[2]); 
  Ghat[13] += 0.2500000000000001*alpha[16]*fUpwind[46]+0.25*(alpha[9]*fUpwind[40]+alpha[8]*fUpwind[39]+alpha[5]*fUpwind[34])+0.223606797749979*(alpha[33]*fUpwind[33]+alpha[32]*fUpwind[32])+0.2500000000000001*(alpha[4]*fUpwind[27]+alpha[2]*fUpwind[24]+alpha[1]*fUpwind[23])+0.223606797749979*(alpha[22]*fUpwind[22]+alpha[21]*fUpwind[21]+alpha[15]*fUpwind[15])+0.25*alpha[0]*fUpwind[13]+0.223606797749979*(alpha[7]*fUpwind[7]+alpha[6]*fUpwind[6]+alpha[3]*fUpwind[3]); 
  Ghat[14] += 0.2500000000000001*alpha[15]*fUpwind[47]+0.25*(alpha[7]*fUpwind[43]+alpha[6]*fUpwind[42]+alpha[5]*fUpwind[41])+0.223606797749979*(alpha[36]*fUpwind[36]+alpha[35]*fUpwind[35])+0.2500000000000001*(alpha[3]*fUpwind[30]+alpha[2]*fUpwind[29]+alpha[1]*fUpwind[28])+0.223606797749979*(alpha[26]*fUpwind[26]+alpha[25]*fUpwind[25]+alpha[16]*fUpwind[16])+0.25*alpha[0]*fUpwind[14]+0.223606797749979*(alpha[9]*fUpwind[9]+alpha[8]*fUpwind[8]+alpha[4]*fUpwind[4]); 
  Ghat[15] += (0.2*alpha[35]+0.223606797749979*alpha[9])*fUpwind[45]+(0.2*alpha[36]+0.223606797749979*alpha[8])*fUpwind[44]+0.223606797749979*(alpha[16]*(fUpwind[38]+fUpwind[37])+fUpwind[18]*alpha[36]+fUpwind[17]*alpha[35])+(0.2*(alpha[22]+alpha[21])+0.223606797749979*alpha[3])*fUpwind[34]+(0.2*alpha[19]+0.223606797749979*alpha[2])*fUpwind[33]+(0.2*(fUpwind[24]+fUpwind[19])+0.223606797749979*fUpwind[2])*alpha[33]+(0.2*alpha[20]+0.223606797749979*alpha[1])*fUpwind[32]+(0.2*(fUpwind[23]+fUpwind[20])+0.223606797749979*fUpwind[1])*alpha[32]+(0.223606797749979*(alpha[26]+alpha[25])+0.25*alpha[4])*fUpwind[31]+0.223606797749979*(alpha[6]*fUpwind[24]+alpha[7]*fUpwind[23]+alpha[5]*fUpwind[22]+fUpwind[5]*alpha[22]+alpha[5]*fUpwind[21]+fUpwind[5]*alpha[21]+alpha[7]*fUpwind[20]+fUpwind[7]*alpha[20]+alpha[6]*fUpwind[19]+fUpwind[6]*alpha[19])+0.25*(alpha[8]*fUpwind[18]+alpha[9]*fUpwind[17]+fUpwind[10]*alpha[16])+(0.223606797749979*(alpha[12]+alpha[11])+0.25*alpha[0])*fUpwind[15]+0.223606797749979*(fUpwind[13]+fUpwind[12]+fUpwind[11])*alpha[15]+0.25*(fUpwind[0]*alpha[15]+alpha[1]*fUpwind[7]+fUpwind[1]*alpha[7]+alpha[2]*fUpwind[6]+fUpwind[2]*alpha[6]+alpha[3]*fUpwind[5]+fUpwind[3]*alpha[5]); 
  Ghat[16] += (0.2*alpha[32]+0.223606797749979*alpha[7])*fUpwind[45]+(0.2*alpha[33]+0.223606797749979*alpha[6])*fUpwind[44]+0.2*(alpha[26]+alpha[25])*fUpwind[41]+0.223606797749979*(alpha[4]*fUpwind[41]+alpha[15]*(fUpwind[38]+fUpwind[37]))+(0.2*alpha[19]+0.223606797749979*alpha[2])*fUpwind[36]+(0.2*(fUpwind[29]+fUpwind[19])+0.223606797749979*fUpwind[2])*alpha[36]+(0.2*alpha[20]+0.223606797749979*alpha[1])*fUpwind[35]+0.2*(fUpwind[28]+fUpwind[20])*alpha[35]+0.223606797749979*(fUpwind[1]*alpha[35]+fUpwind[18]*alpha[33]+fUpwind[17]*alpha[32])+(0.223606797749979*(alpha[22]+alpha[21])+0.25*alpha[3])*fUpwind[31]+0.223606797749979*(alpha[8]*fUpwind[29]+alpha[9]*fUpwind[28]+alpha[5]*fUpwind[26]+fUpwind[5]*alpha[26]+alpha[5]*fUpwind[25]+fUpwind[5]*alpha[25]+alpha[9]*fUpwind[20]+fUpwind[9]*alpha[20]+alpha[8]*fUpwind[19]+fUpwind[8]*alpha[19])+0.25*(alpha[6]*fUpwind[18]+alpha[7]*fUpwind[17])+(0.223606797749979*(alpha[12]+alpha[11])+0.25*alpha[0])*fUpwind[16]+0.223606797749979*(fUpwind[14]+fUpwind[12]+fUpwind[11])*alpha[16]+0.25*(fUpwind[0]*alpha[16]+fUpwind[10]*alpha[15]+alpha[1]*fUpwind[9]+fUpwind[1]*alpha[9]+alpha[2]*fUpwind[8]+fUpwind[2]*alpha[8]+alpha[4]*fUpwind[5]+fUpwind[4]*alpha[5]); 
  Ghat[17] += (0.2*alpha[35]+0.223606797749979*alpha[9])*fUpwind[47]+(0.2*alpha[32]+0.223606797749979*alpha[7])*fUpwind[46]+0.2500000000000001*alpha[12]*fUpwind[45]+0.223606797749979*alpha[5]*fUpwind[44]+0.223606797749979*alpha[16]*fUpwind[43]+0.2*alpha[25]*fUpwind[42]+0.223606797749979*(alpha[4]*fUpwind[42]+alpha[15]*fUpwind[40])+(0.2*alpha[21]+0.223606797749979*alpha[3])*fUpwind[39]+0.2500000000000001*alpha[20]*fUpwind[38]+0.223606797749979*alpha[1]*fUpwind[37]+0.2500000000000001*(alpha[22]*fUpwind[36]+fUpwind[22]*alpha[36])+0.223606797749979*(alpha[15]*fUpwind[35]+fUpwind[15]*alpha[35])+0.2500000000000001*(alpha[26]*fUpwind[33]+fUpwind[26]*alpha[33])+0.223606797749979*(alpha[16]*fUpwind[32]+fUpwind[16]*alpha[32])+(0.223606797749979*alpha[19]+0.25*alpha[2])*fUpwind[31]+0.223606797749979*(alpha[8]*fUpwind[30]+alpha[6]*(fUpwind[27]+fUpwind[25])+fUpwind[6]*alpha[25]+alpha[8]*fUpwind[21]+fUpwind[8]*alpha[21])+0.25*alpha[5]*fUpwind[18]+0.223606797749979*alpha[11]*fUpwind[17]+0.25*(alpha[0]*fUpwind[17]+alpha[7]*fUpwind[16]+fUpwind[7]*alpha[16]+alpha[9]*fUpwind[15]+fUpwind[9]*alpha[15]+alpha[1]*fUpwind[10]+alpha[3]*fUpwind[8]+fUpwind[3]*alpha[8]+alpha[4]*fUpwind[6]+fUpwind[4]*alpha[6]); 
  Ghat[18] += (0.2*alpha[36]+0.223606797749979*alpha[8])*fUpwind[47]+0.2*alpha[33]*fUpwind[46]+0.223606797749979*(alpha[6]*fUpwind[46]+alpha[5]*fUpwind[45])+0.2500000000000001*alpha[11]*fUpwind[44]+0.2*alpha[26]*fUpwind[43]+0.223606797749979*(alpha[4]*fUpwind[43]+alpha[16]*fUpwind[42])+0.2*alpha[22]*fUpwind[40]+0.223606797749979*(alpha[3]*fUpwind[40]+alpha[15]*fUpwind[39]+alpha[2]*fUpwind[38])+0.2500000000000001*alpha[19]*fUpwind[37]+0.223606797749979*(alpha[15]*fUpwind[36]+fUpwind[15]*alpha[36])+0.2500000000000001*(alpha[21]*fUpwind[35]+fUpwind[21]*alpha[35])+0.223606797749979*(alpha[16]*fUpwind[33]+fUpwind[16]*alpha[33])+0.2500000000000001*(alpha[25]*fUpwind[32]+fUpwind[25]*alpha[32])+(0.223606797749979*alpha[20]+0.25*alpha[1])*fUpwind[31]+0.223606797749979*(alpha[9]*fUpwind[30]+alpha[7]*(fUpwind[27]+fUpwind[26])+fUpwind[7]*alpha[26]+alpha[9]*fUpwind[22]+fUpwind[9]*alpha[22])+0.223606797749979*alpha[12]*fUpwind[18]+0.25*(alpha[0]*fUpwind[18]+alpha[5]*fUpwind[17]+alpha[6]*fUpwind[16]+fUpwind[6]*alpha[16]+alpha[8]*fUpwind[15]+fUpwind[8]*alpha[15]+alpha[2]*fUpwind[10]+alpha[3]*fUpwind[9]+fUpwind[3]*alpha[9]+alpha[4]*fUpwind[7]+fUpwind[4]*alpha[7]); 
  Ghat[19] += 0.2*(alpha[16]*fUpwind[36]+fUpwind[16]*alpha[36])+(0.223606797749979*alpha[26]+0.159719141249985*alpha[25]+0.2500000000000001*alpha[4])*fUpwind[35]+(0.223606797749979*fUpwind[26]+0.159719141249985*fUpwind[25]+0.2500000000000001*fUpwind[4])*alpha[35]+0.2*(alpha[15]*fUpwind[33]+fUpwind[15]*alpha[33])+(0.223606797749979*alpha[22]+0.159719141249985*alpha[21]+0.2500000000000001*alpha[3])*fUpwind[32]+(0.223606797749979*fUpwind[22]+0.159719141249985*fUpwind[21]+0.2500000000000001*fUpwind[3])*alpha[32]+0.25*(alpha[9]*fUpwind[25]+fUpwind[9]*alpha[25]+alpha[7]*fUpwind[21]+fUpwind[7]*alpha[21])+0.2*(alpha[5]*fUpwind[20]+fUpwind[5]*alpha[20])+(0.223606797749979*alpha[12]+0.159719141249985*alpha[11]+0.25*alpha[0])*fUpwind[19]+(0.223606797749979*fUpwind[12]+0.159719141249985*fUpwind[11]+0.25*fUpwind[0])*alpha[19]+0.223606797749979*(alpha[8]*fUpwind[16]+fUpwind[8]*alpha[16]+alpha[6]*fUpwind[15]+fUpwind[6]*alpha[15])+0.2500000000000001*(alpha[2]*fUpwind[11]+fUpwind[2]*alpha[11])+0.223606797749979*(alpha[1]*fUpwind[5]+fUpwind[1]*alpha[5]); 
  Ghat[20] += (0.159719141249985*alpha[26]+0.223606797749979*alpha[25]+0.2500000000000001*alpha[4])*fUpwind[36]+(0.159719141249985*fUpwind[26]+0.223606797749979*fUpwind[25]+0.2500000000000001*fUpwind[4])*alpha[36]+0.2*(alpha[16]*fUpwind[35]+fUpwind[16]*alpha[35])+(0.159719141249985*alpha[22]+0.223606797749979*alpha[21]+0.2500000000000001*alpha[3])*fUpwind[33]+(0.159719141249985*fUpwind[22]+0.223606797749979*fUpwind[21]+0.2500000000000001*fUpwind[3])*alpha[33]+0.2*(alpha[15]*fUpwind[32]+fUpwind[15]*alpha[32])+0.25*(alpha[8]*fUpwind[26]+fUpwind[8]*alpha[26]+alpha[6]*fUpwind[22]+fUpwind[6]*alpha[22])+(0.159719141249985*alpha[12]+0.223606797749979*alpha[11]+0.25*alpha[0])*fUpwind[20]+(0.159719141249985*fUpwind[12]+0.223606797749979*fUpwind[11]+0.25*fUpwind[0])*alpha[20]+0.2*(alpha[5]*fUpwind[19]+fUpwind[5]*alpha[19])+0.223606797749979*(alpha[9]*fUpwind[16]+fUpwind[9]*alpha[16]+alpha[7]*fUpwind[15]+fUpwind[7]*alpha[15])+0.2500000000000001*(alpha[1]*fUpwind[12]+fUpwind[1]*alpha[12])+0.223606797749979*(alpha[2]*fUpwind[5]+fUpwind[2]*alpha[5]); 
  Ghat[21] += 0.223606797749979*alpha[36]*fUpwind[45]+(0.159719141249985*alpha[35]+0.25*alpha[9])*fUpwind[44]+0.159719141249985*alpha[25]*fUpwind[37]+0.2500000000000001*(alpha[4]*fUpwind[37]+fUpwind[18]*alpha[35])+0.2*alpha[15]*fUpwind[34]+0.223606797749979*(alpha[20]*fUpwind[33]+fUpwind[20]*alpha[33])+(0.159719141249985*alpha[19]+0.2500000000000001*alpha[2])*fUpwind[32]+(0.223606797749979*fUpwind[24]+0.159719141249985*fUpwind[19]+0.2500000000000001*fUpwind[2])*alpha[32]+0.223606797749979*alpha[16]*fUpwind[31]+0.25*fUpwind[10]*alpha[25]+0.2*alpha[6]*fUpwind[23]+(0.159719141249985*alpha[11]+0.25*alpha[0])*fUpwind[21]+(0.223606797749979*fUpwind[13]+0.159719141249985*fUpwind[11])*alpha[21]+0.25*(fUpwind[0]*alpha[21]+alpha[7]*fUpwind[19]+fUpwind[7]*alpha[19])+0.223606797749979*(alpha[8]*fUpwind[17]+alpha[5]*fUpwind[15]+fUpwind[5]*alpha[15])+0.2500000000000001*(alpha[3]*fUpwind[11]+fUpwind[3]*alpha[11])+0.223606797749979*(alpha[1]*fUpwind[6]+fUpwind[1]*alpha[6]); 
  Ghat[22] += (0.159719141249985*alpha[36]+0.25*alpha[8])*fUpwind[45]+0.223606797749979*alpha[35]*fUpwind[44]+0.159719141249985*alpha[26]*fUpwind[38]+0.2500000000000001*(alpha[4]*fUpwind[38]+fUpwind[17]*alpha[36])+0.2*alpha[15]*fUpwind[34]+(0.159719141249985*alpha[20]+0.2500000000000001*alpha[1])*fUpwind[33]+(0.223606797749979*fUpwind[23]+0.159719141249985*fUpwind[20]+0.2500000000000001*fUpwind[1])*alpha[33]+0.223606797749979*(alpha[19]*fUpwind[32]+fUpwind[19]*alpha[32])+0.223606797749979*alpha[16]*fUpwind[31]+0.25*fUpwind[10]*alpha[26]+0.2*alpha[7]*fUpwind[24]+(0.159719141249985*alpha[12]+0.25*alpha[0])*fUpwind[22]+(0.223606797749979*fUpwind[13]+0.159719141249985*fUpwind[12])*alpha[22]+0.25*(fUpwind[0]*alpha[22]+alpha[6]*fUpwind[20]+fUpwind[6]*alpha[20])+0.223606797749979*(alpha[9]*fUpwind[18]+alpha[5]*fUpwind[15]+fUpwind[5]*alpha[15])+0.2500000000000001*(alpha[3]*fUpwind[12]+fUpwind[3]*alpha[12])+0.223606797749979*(alpha[2]*fUpwind[7]+fUpwind[2]*alpha[7]); 
  Ghat[23] += (0.223606797749979*alpha[35]+0.25*alpha[9])*fUpwind[46]+0.2500000000000001*alpha[16]*fUpwind[40]+(0.223606797749979*alpha[25]+0.2500000000000001*alpha[4])*fUpwind[39]+(0.223606797749979*alpha[19]+0.2500000000000001*alpha[2])*fUpwind[34]+0.223606797749979*(alpha[22]*fUpwind[33]+fUpwind[22]*alpha[33])+0.2*(alpha[15]*fUpwind[32]+fUpwind[15]*alpha[32])+0.25*(alpha[8]*fUpwind[27]+alpha[5]*fUpwind[24])+(0.223606797749979*alpha[11]+0.25*alpha[0])*fUpwind[23]+0.2*(alpha[6]*fUpwind[21]+fUpwind[6]*alpha[21])+0.223606797749979*(alpha[7]*fUpwind[15]+fUpwind[7]*alpha[15])+0.2500000000000001*alpha[1]*fUpwind[13]+0.223606797749979*(alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6]); 
  Ghat[24] += (0.223606797749979*alpha[36]+0.25*alpha[8])*fUpwind[46]+0.223606797749979*alpha[26]*fUpwind[40]+0.2500000000000001*(alpha[4]*fUpwind[40]+alpha[16]*fUpwind[39])+(0.223606797749979*alpha[20]+0.2500000000000001*alpha[1])*fUpwind[34]+0.2*(alpha[15]*fUpwind[33]+fUpwind[15]*alpha[33])+0.223606797749979*(alpha[21]*fUpwind[32]+fUpwind[21]*alpha[32])+0.25*alpha[9]*fUpwind[27]+0.223606797749979*alpha[12]*fUpwind[24]+0.25*(alpha[0]*fUpwind[24]+alpha[5]*fUpwind[23])+0.2*(alpha[7]*fUpwind[22]+fUpwind[7]*alpha[22])+0.223606797749979*(alpha[6]*fUpwind[15]+fUpwind[6]*alpha[15])+0.2500000000000001*alpha[2]*fUpwind[13]+0.223606797749979*(alpha[3]*fUpwind[7]+fUpwind[3]*alpha[7]); 
  Ghat[25] += 0.223606797749979*alpha[33]*fUpwind[45]+(0.159719141249985*alpha[32]+0.25*alpha[7])*fUpwind[44]+0.2*alpha[16]*fUpwind[41]+(0.159719141249985*alpha[21]+0.2500000000000001*alpha[3])*fUpwind[37]+0.223606797749979*(alpha[20]*fUpwind[36]+fUpwind[20]*alpha[36])+(0.159719141249985*alpha[19]+0.2500000000000001*alpha[2])*fUpwind[35]+(0.223606797749979*fUpwind[29]+0.159719141249985*fUpwind[19])*alpha[35]+0.2500000000000001*(fUpwind[2]*alpha[35]+fUpwind[18]*alpha[32])+0.223606797749979*alpha[15]*fUpwind[31]+0.2*alpha[8]*fUpwind[28]+(0.159719141249985*alpha[11]+0.25*alpha[0])*fUpwind[25]+(0.223606797749979*fUpwind[14]+0.159719141249985*fUpwind[11])*alpha[25]+0.25*(fUpwind[0]*alpha[25]+fUpwind[10]*alpha[21]+alpha[9]*fUpwind[19]+fUpwind[9]*alpha[19])+0.223606797749979*(alpha[6]*fUpwind[17]+alpha[5]*fUpwind[16]+fUpwind[5]*alpha[16])+0.2500000000000001*(alpha[4]*fUpwind[11]+fUpwind[4]*alpha[11])+0.223606797749979*(alpha[1]*fUpwind[8]+fUpwind[1]*alpha[8]); 
  Ghat[26] += (0.159719141249985*alpha[33]+0.25*alpha[6])*fUpwind[45]+0.223606797749979*alpha[32]*fUpwind[44]+0.2*alpha[16]*fUpwind[41]+(0.159719141249985*alpha[22]+0.2500000000000001*alpha[3])*fUpwind[38]+(0.159719141249985*alpha[20]+0.2500000000000001*alpha[1])*fUpwind[36]+(0.223606797749979*fUpwind[28]+0.159719141249985*fUpwind[20]+0.2500000000000001*fUpwind[1])*alpha[36]+0.223606797749979*(alpha[19]*fUpwind[35]+fUpwind[19]*alpha[35])+0.2500000000000001*fUpwind[17]*alpha[33]+0.223606797749979*alpha[15]*fUpwind[31]+0.2*alpha[9]*fUpwind[29]+(0.159719141249985*alpha[12]+0.25*alpha[0])*fUpwind[26]+(0.223606797749979*fUpwind[14]+0.159719141249985*fUpwind[12])*alpha[26]+0.25*(fUpwind[0]*alpha[26]+fUpwind[10]*alpha[22]+alpha[8]*fUpwind[20]+fUpwind[8]*alpha[20])+0.223606797749979*(alpha[7]*fUpwind[18]+alpha[5]*fUpwind[16]+fUpwind[5]*alpha[16])+0.2500000000000001*(alpha[4]*fUpwind[12]+fUpwind[4]*alpha[12])+0.223606797749979*(alpha[2]*fUpwind[9]+fUpwind[2]*alpha[9]); 
  Ghat[27] += 0.25*alpha[5]*fUpwind[46]+0.223606797749979*(alpha[33]*fUpwind[45]+alpha[32]*fUpwind[44])+0.2500000000000001*(alpha[2]*fUpwind[40]+alpha[1]*fUpwind[39])+0.223606797749979*(alpha[22]*fUpwind[38]+alpha[21]*fUpwind[37])+0.2500000000000001*alpha[16]*fUpwind[34]+0.223606797749979*alpha[15]*fUpwind[31]+0.25*(alpha[0]*fUpwind[27]+alpha[9]*fUpwind[24]+alpha[8]*fUpwind[23])+0.223606797749979*(alpha[7]*fUpwind[18]+alpha[6]*fUpwind[17])+0.2500000000000001*alpha[4]*fUpwind[13]+0.223606797749979*alpha[3]*fUpwind[10]; 
  Ghat[28] += (0.223606797749979*alpha[32]+0.25*alpha[7])*fUpwind[47]+0.2500000000000001*alpha[15]*fUpwind[43]+(0.223606797749979*alpha[21]+0.2500000000000001*alpha[3])*fUpwind[42]+(0.223606797749979*alpha[19]+0.2500000000000001*alpha[2])*fUpwind[41]+0.223606797749979*(alpha[26]*fUpwind[36]+fUpwind[26]*alpha[36])+0.2*(alpha[16]*fUpwind[35]+fUpwind[16]*alpha[35])+0.25*(alpha[6]*fUpwind[30]+alpha[5]*fUpwind[29])+(0.223606797749979*alpha[11]+0.25*alpha[0])*fUpwind[28]+0.2*(alpha[8]*fUpwind[25]+fUpwind[8]*alpha[25])+0.223606797749979*(alpha[9]*fUpwind[16]+fUpwind[9]*alpha[16])+0.2500000000000001*alpha[1]*fUpwind[14]+0.223606797749979*(alpha[4]*fUpwind[8]+fUpwind[4]*alpha[8]); 
  Ghat[29] += (0.223606797749979*alpha[33]+0.25*alpha[6])*fUpwind[47]+0.223606797749979*alpha[22]*fUpwind[43]+0.2500000000000001*(alpha[3]*fUpwind[43]+alpha[15]*fUpwind[42])+(0.223606797749979*alpha[20]+0.2500000000000001*alpha[1])*fUpwind[41]+0.2*(alpha[16]*fUpwind[36]+fUpwind[16]*alpha[36])+0.223606797749979*(alpha[25]*fUpwind[35]+fUpwind[25]*alpha[35])+0.25*alpha[7]*fUpwind[30]+0.223606797749979*alpha[12]*fUpwind[29]+0.25*(alpha[0]*fUpwind[29]+alpha[5]*fUpwind[28])+0.2*(alpha[9]*fUpwind[26]+fUpwind[9]*alpha[26])+0.223606797749979*(alpha[8]*fUpwind[16]+fUpwind[8]*alpha[16])+0.2500000000000001*alpha[2]*fUpwind[14]+0.223606797749979*(alpha[4]*fUpwind[9]+fUpwind[4]*alpha[9]); 
  Ghat[30] += 0.25*alpha[5]*fUpwind[47]+0.223606797749979*(alpha[36]*fUpwind[45]+alpha[35]*fUpwind[44])+0.2500000000000001*(alpha[2]*fUpwind[43]+alpha[1]*fUpwind[42]+alpha[15]*fUpwind[41])+0.223606797749979*(alpha[26]*fUpwind[38]+alpha[25]*fUpwind[37])+0.223606797749979*alpha[16]*fUpwind[31]+0.25*(alpha[0]*fUpwind[30]+alpha[7]*fUpwind[29]+alpha[6]*fUpwind[28])+0.223606797749979*(alpha[9]*fUpwind[18]+alpha[8]*fUpwind[17])+0.2500000000000001*alpha[3]*fUpwind[14]+0.223606797749979*alpha[4]*fUpwind[10]; 
  Ghat[31] += (0.2*(alpha[26]+alpha[25])+0.223606797749979*alpha[4])*fUpwind[47]+(0.2*(alpha[22]+alpha[21])+0.223606797749979*alpha[3])*fUpwind[46]+(0.2*alpha[19]+0.223606797749979*alpha[2])*fUpwind[45]+(0.2*alpha[20]+0.223606797749979*alpha[1])*fUpwind[44]+(0.2*alpha[36]+0.223606797749979*alpha[8])*fUpwind[43]+(0.2*alpha[35]+0.223606797749979*alpha[9])*fUpwind[42]+(0.2*alpha[33]+0.223606797749979*alpha[6])*fUpwind[40]+0.2*alpha[32]*fUpwind[39]+0.223606797749979*(alpha[7]*fUpwind[39]+alpha[5]*(fUpwind[38]+fUpwind[37]))+(0.2*alpha[32]+0.223606797749979*alpha[7])*fUpwind[36]+(0.2*fUpwind[32]+0.223606797749979*fUpwind[7])*alpha[36]+(0.2*alpha[33]+0.223606797749979*alpha[6])*fUpwind[35]+0.2*fUpwind[33]*alpha[35]+0.223606797749979*(fUpwind[6]*alpha[35]+alpha[9]*fUpwind[33]+fUpwind[9]*alpha[33]+alpha[8]*fUpwind[32]+fUpwind[8]*alpha[32])+(0.223606797749979*(alpha[12]+alpha[11])+0.25*alpha[0])*fUpwind[31]+0.223606797749979*(alpha[16]*fUpwind[30]+alpha[15]*(fUpwind[27]+fUpwind[26])+fUpwind[15]*alpha[26]+alpha[15]*fUpwind[25]+fUpwind[15]*alpha[25]+alpha[16]*fUpwind[22]+fUpwind[16]*alpha[22]+alpha[16]*fUpwind[21]+fUpwind[16]*alpha[21]+fUpwind[18]*alpha[20]+fUpwind[17]*alpha[19])+0.25*(alpha[1]*fUpwind[18]+alpha[2]*fUpwind[17]+alpha[3]*fUpwind[16]+fUpwind[3]*alpha[16]+alpha[4]*fUpwind[15]+fUpwind[4]*alpha[15]+alpha[5]*fUpwind[10]+alpha[6]*fUpwind[9]+fUpwind[6]*alpha[9]+alpha[7]*fUpwind[8]+fUpwind[7]*alpha[8]); 
  Ghat[32] += 0.2*alpha[16]*fUpwind[45]+(0.223606797749979*alpha[26]+0.159719141249985*alpha[25]+0.2500000000000001*alpha[4])*fUpwind[44]+0.223606797749979*alpha[35]*fUpwind[38]+(0.159719141249985*alpha[35]+0.25*alpha[9])*fUpwind[37]+0.2*fUpwind[31]*alpha[36]+0.25*fUpwind[10]*alpha[35]+0.1788854381999831*alpha[33]*fUpwind[34]+0.2*(alpha[6]*fUpwind[34]+alpha[5]*fUpwind[33]+fUpwind[5]*alpha[33])+(0.223606797749979*alpha[12]+0.159719141249985*alpha[11]+0.25*alpha[0])*fUpwind[32]+(0.223606797749979*(fUpwind[13]+fUpwind[12])+0.159719141249985*fUpwind[11]+0.25*fUpwind[0])*alpha[32]+0.223606797749979*alpha[8]*fUpwind[31]+0.2500000000000001*fUpwind[18]*alpha[25]+0.223606797749979*alpha[21]*fUpwind[24]+0.2*alpha[15]*fUpwind[23]+0.223606797749979*(alpha[19]*fUpwind[22]+fUpwind[19]*alpha[22])+(0.159719141249985*alpha[19]+0.2500000000000001*alpha[2])*fUpwind[21]+(0.159719141249985*fUpwind[19]+0.2500000000000001*fUpwind[2])*alpha[21]+0.2*(alpha[15]*fUpwind[20]+fUpwind[15]*alpha[20])+0.2500000000000001*(alpha[3]*fUpwind[19]+fUpwind[3]*alpha[19])+0.223606797749979*(alpha[16]*fUpwind[17]+alpha[1]*fUpwind[15]+fUpwind[1]*alpha[15])+0.25*(alpha[7]*fUpwind[11]+fUpwind[7]*alpha[11])+0.223606797749979*(alpha[5]*fUpwind[6]+fUpwind[5]*alpha[6]); 
  Ghat[33] += (0.159719141249985*alpha[26]+0.223606797749979*alpha[25]+0.2500000000000001*alpha[4])*fUpwind[45]+0.2*alpha[16]*fUpwind[44]+(0.159719141249985*alpha[36]+0.25*alpha[8])*fUpwind[38]+alpha[36]*(0.223606797749979*fUpwind[37]+0.25*fUpwind[10])+0.2*fUpwind[31]*alpha[35]+(0.1788854381999831*alpha[32]+0.2*alpha[7])*fUpwind[34]+(0.159719141249985*alpha[12]+0.223606797749979*alpha[11]+0.25*alpha[0])*fUpwind[33]+(0.223606797749979*fUpwind[13]+0.159719141249985*fUpwind[12]+0.223606797749979*fUpwind[11]+0.25*fUpwind[0])*alpha[33]+0.2*(alpha[5]*fUpwind[32]+fUpwind[5]*alpha[32])+0.223606797749979*alpha[9]*fUpwind[31]+0.2500000000000001*fUpwind[17]*alpha[26]+0.2*alpha[15]*fUpwind[24]+0.223606797749979*alpha[22]*fUpwind[23]+(0.159719141249985*alpha[20]+0.2500000000000001*alpha[1])*fUpwind[22]+(0.159719141249985*fUpwind[20]+0.2500000000000001*fUpwind[1])*alpha[22]+0.223606797749979*(alpha[20]*fUpwind[21]+fUpwind[20]*alpha[21])+0.2500000000000001*(alpha[3]*fUpwind[20]+fUpwind[3]*alpha[20])+0.2*(alpha[15]*fUpwind[19]+fUpwind[15]*alpha[19])+0.223606797749979*(alpha[16]*fUpwind[18]+alpha[2]*fUpwind[15]+fUpwind[2]*alpha[15])+0.25*(alpha[6]*fUpwind[12]+fUpwind[6]*alpha[12])+0.223606797749979*(alpha[5]*fUpwind[7]+fUpwind[5]*alpha[7]); 
  Ghat[34] += (0.223606797749979*(alpha[26]+alpha[25])+0.2500000000000001*alpha[4])*fUpwind[46]+(0.223606797749979*alpha[36]+0.25*alpha[8])*fUpwind[40]+(0.223606797749979*alpha[35]+0.25*alpha[9])*fUpwind[39]+(0.223606797749979*(alpha[12]+alpha[11])+0.25*alpha[0])*fUpwind[34]+(0.1788854381999831*alpha[32]+0.2*alpha[7])*fUpwind[33]+0.1788854381999831*fUpwind[32]*alpha[33]+0.2*(fUpwind[7]*alpha[33]+alpha[6]*fUpwind[32]+fUpwind[6]*alpha[32])+0.2500000000000001*alpha[16]*fUpwind[27]+(0.223606797749979*alpha[20]+0.2500000000000001*alpha[1])*fUpwind[24]+(0.223606797749979*alpha[19]+0.2500000000000001*alpha[2])*fUpwind[23]+0.2*(alpha[15]*fUpwind[22]+fUpwind[15]*alpha[22]+alpha[15]*fUpwind[21]+fUpwind[15]*alpha[21])+0.223606797749979*(alpha[3]*fUpwind[15]+fUpwind[3]*alpha[15])+0.25*alpha[5]*fUpwind[13]+0.223606797749979*(alpha[6]*fUpwind[7]+fUpwind[6]*alpha[7]); 
  Ghat[35] += 0.2*alpha[15]*fUpwind[45]+(0.223606797749979*alpha[22]+0.159719141249985*alpha[21]+0.2500000000000001*alpha[3])*fUpwind[44]+(0.1788854381999831*alpha[36]+0.2*alpha[8])*fUpwind[41]+0.223606797749979*alpha[32]*fUpwind[38]+(0.159719141249985*alpha[32]+0.25*alpha[7])*fUpwind[37]+0.2*(alpha[5]*fUpwind[36]+fUpwind[5]*alpha[36])+(0.223606797749979*alpha[12]+0.159719141249985*alpha[11]+0.25*alpha[0])*fUpwind[35]+(0.223606797749979*(fUpwind[14]+fUpwind[12])+0.159719141249985*fUpwind[11]+0.25*fUpwind[0])*alpha[35]+0.2*fUpwind[31]*alpha[33]+0.25*fUpwind[10]*alpha[32]+0.223606797749979*(alpha[6]*fUpwind[31]+alpha[25]*fUpwind[29])+0.2*alpha[16]*fUpwind[28]+0.223606797749979*(alpha[19]*fUpwind[26]+fUpwind[19]*alpha[26])+(0.159719141249985*alpha[19]+0.2500000000000001*alpha[2])*fUpwind[25]+0.159719141249985*fUpwind[19]*alpha[25]+0.2500000000000001*(fUpwind[2]*alpha[25]+fUpwind[18]*alpha[21])+0.2*(alpha[16]*fUpwind[20]+fUpwind[16]*alpha[20])+0.2500000000000001*(alpha[4]*fUpwind[19]+fUpwind[4]*alpha[19])+0.223606797749979*(alpha[15]*fUpwind[17]+alpha[1]*fUpwind[16]+fUpwind[1]*alpha[16])+0.25*(alpha[9]*fUpwind[11]+fUpwind[9]*alpha[11])+0.223606797749979*(alpha[5]*fUpwind[8]+fUpwind[5]*alpha[8]); 
  Ghat[36] += (0.159719141249985*alpha[22]+0.223606797749979*alpha[21]+0.2500000000000001*alpha[3])*fUpwind[45]+0.2*alpha[15]*fUpwind[44]+(0.1788854381999831*alpha[35]+0.2*alpha[9])*fUpwind[41]+(0.159719141249985*alpha[33]+0.25*alpha[6])*fUpwind[38]+0.223606797749979*alpha[33]*fUpwind[37]+(0.159719141249985*alpha[12]+0.223606797749979*alpha[11]+0.25*alpha[0])*fUpwind[36]+(0.223606797749979*fUpwind[14]+0.159719141249985*fUpwind[12]+0.223606797749979*fUpwind[11]+0.25*fUpwind[0])*alpha[36]+0.2*(alpha[5]*fUpwind[35]+fUpwind[5]*alpha[35])+0.25*fUpwind[10]*alpha[33]+fUpwind[31]*(0.2*alpha[32]+0.223606797749979*alpha[7])+0.2*alpha[16]*fUpwind[29]+0.223606797749979*alpha[26]*fUpwind[28]+(0.159719141249985*alpha[20]+0.2500000000000001*alpha[1])*fUpwind[26]+(0.159719141249985*fUpwind[20]+0.2500000000000001*fUpwind[1])*alpha[26]+0.223606797749979*(alpha[20]*fUpwind[25]+fUpwind[20]*alpha[25])+0.2500000000000001*(fUpwind[17]*alpha[22]+alpha[4]*fUpwind[20]+fUpwind[4]*alpha[20])+0.2*(alpha[16]*fUpwind[19]+fUpwind[16]*alpha[19])+0.223606797749979*(alpha[15]*fUpwind[18]+alpha[2]*fUpwind[16]+fUpwind[2]*alpha[16])+0.25*(alpha[8]*fUpwind[12]+fUpwind[8]*alpha[12])+0.223606797749979*(alpha[5]*fUpwind[9]+fUpwind[5]*alpha[9]); 
  Ghat[37] += 0.2*(alpha[16]*fUpwind[47]+alpha[15]*fUpwind[46])+0.223606797749979*alpha[20]*fUpwind[45]+(0.159719141249985*alpha[19]+0.2500000000000001*alpha[2])*fUpwind[44]+0.223606797749979*alpha[35]*fUpwind[43]+0.2*alpha[8]*fUpwind[42]+0.223606797749979*alpha[32]*fUpwind[40]+0.2*alpha[6]*fUpwind[39]+(0.159719141249985*alpha[11]+0.25*alpha[0])*fUpwind[37]+0.223606797749979*(alpha[33]*fUpwind[36]+fUpwind[33]*alpha[36])+(0.159719141249985*alpha[32]+0.25*alpha[7])*fUpwind[35]+0.159719141249985*fUpwind[32]*alpha[35]+0.25*(fUpwind[7]*alpha[35]+alpha[9]*fUpwind[32]+fUpwind[9]*alpha[32])+0.223606797749979*(alpha[5]*fUpwind[31]+alpha[25]*fUpwind[30]+alpha[21]*fUpwind[27])+(0.159719141249985*alpha[21]+0.2500000000000001*alpha[3])*fUpwind[25]+0.159719141249985*fUpwind[21]*alpha[25]+0.2500000000000001*(fUpwind[3]*alpha[25]+alpha[4]*fUpwind[21]+fUpwind[4]*alpha[21]+fUpwind[18]*alpha[19])+0.223606797749979*(alpha[1]*fUpwind[17]+alpha[15]*fUpwind[16]+fUpwind[15]*alpha[16])+0.25*fUpwind[10]*alpha[11]+0.223606797749979*(alpha[6]*fUpwind[8]+fUpwind[6]*alpha[8]); 
  Ghat[38] += 0.2*(alpha[16]*fUpwind[47]+alpha[15]*fUpwind[46])+(0.159719141249985*alpha[20]+0.2500000000000001*alpha[1])*fUpwind[45]+0.223606797749979*alpha[19]*fUpwind[44]+0.2*alpha[9]*fUpwind[43]+0.223606797749979*alpha[36]*fUpwind[42]+0.2*alpha[7]*fUpwind[40]+0.223606797749979*alpha[33]*fUpwind[39]+(0.159719141249985*alpha[12]+0.25*alpha[0])*fUpwind[38]+(0.159719141249985*alpha[33]+0.25*alpha[6])*fUpwind[36]+(0.159719141249985*fUpwind[33]+0.25*fUpwind[6])*alpha[36]+0.223606797749979*(alpha[32]*fUpwind[35]+fUpwind[32]*alpha[35])+0.25*(alpha[8]*fUpwind[33]+fUpwind[8]*alpha[33])+0.223606797749979*(alpha[5]*fUpwind[31]+alpha[26]*fUpwind[30]+alpha[22]*fUpwind[27])+(0.159719141249985*alpha[22]+0.2500000000000001*alpha[3])*fUpwind[26]+0.159719141249985*fUpwind[22]*alpha[26]+0.2500000000000001*(fUpwind[3]*alpha[26]+alpha[4]*fUpwind[22]+fUpwind[4]*alpha[22]+fUpwind[17]*alpha[20])+0.223606797749979*(alpha[2]*fUpwind[18]+alpha[15]*fUpwind[16]+fUpwind[15]*alpha[16])+0.25*fUpwind[10]*alpha[12]+0.223606797749979*(alpha[7]*fUpwind[9]+fUpwind[7]*alpha[9]); 
  Ghat[39] += (0.223606797749979*alpha[19]+0.2500000000000001*alpha[2])*fUpwind[46]+0.223606797749979*alpha[22]*fUpwind[45]+0.2*alpha[15]*fUpwind[44]+0.25*alpha[5]*fUpwind[40]+(0.223606797749979*alpha[11]+0.25*alpha[0])*fUpwind[39]+0.223606797749979*alpha[33]*fUpwind[38]+0.2*alpha[6]*fUpwind[37]+fUpwind[34]*(0.223606797749979*alpha[35]+0.25*alpha[9])+fUpwind[31]*(0.2*alpha[32]+0.223606797749979*alpha[7])+0.2500000000000001*alpha[1]*fUpwind[27]+0.223606797749979*fUpwind[23]*alpha[25]+0.2500000000000001*(alpha[16]*fUpwind[24]+alpha[4]*fUpwind[23])+0.2*fUpwind[17]*alpha[21]+0.223606797749979*(alpha[15]*fUpwind[18]+alpha[3]*fUpwind[17])+0.25*alpha[8]*fUpwind[13]+0.223606797749979*alpha[6]*fUpwind[10]; 
  Ghat[40] += (0.223606797749979*alpha[20]+0.2500000000000001*alpha[1])*fUpwind[46]+0.2*alpha[15]*fUpwind[45]+0.223606797749979*(alpha[21]*fUpwind[44]+alpha[12]*fUpwind[40])+0.25*(alpha[0]*fUpwind[40]+alpha[5]*fUpwind[39])+0.2*alpha[7]*fUpwind[38]+0.223606797749979*alpha[32]*fUpwind[37]+fUpwind[34]*(0.223606797749979*alpha[36]+0.25*alpha[8])+fUpwind[31]*(0.2*alpha[33]+0.223606797749979*alpha[6])+0.2500000000000001*alpha[2]*fUpwind[27]+0.223606797749979*fUpwind[24]*alpha[26]+0.2500000000000001*(alpha[4]*fUpwind[24]+alpha[16]*fUpwind[23])+0.2*fUpwind[18]*alpha[22]+0.223606797749979*(alpha[3]*fUpwind[18]+alpha[15]*fUpwind[17])+0.25*alpha[9]*fUpwind[13]+0.223606797749979*alpha[7]*fUpwind[10]; 
  Ghat[41] += (0.223606797749979*(alpha[22]+alpha[21])+0.2500000000000001*alpha[3])*fUpwind[47]+(0.223606797749979*alpha[33]+0.25*alpha[6])*fUpwind[43]+(0.223606797749979*alpha[32]+0.25*alpha[7])*fUpwind[42]+(0.223606797749979*(alpha[12]+alpha[11])+0.25*alpha[0])*fUpwind[41]+(0.1788854381999831*alpha[35]+0.2*alpha[9])*fUpwind[36]+0.1788854381999831*fUpwind[35]*alpha[36]+0.2*(fUpwind[9]*alpha[36]+alpha[8]*fUpwind[35]+fUpwind[8]*alpha[35])+0.2500000000000001*alpha[15]*fUpwind[30]+(0.223606797749979*alpha[20]+0.2500000000000001*alpha[1])*fUpwind[29]+(0.223606797749979*alpha[19]+0.2500000000000001*alpha[2])*fUpwind[28]+0.2*(alpha[16]*fUpwind[26]+fUpwind[16]*alpha[26]+alpha[16]*fUpwind[25]+fUpwind[16]*alpha[25])+0.223606797749979*(alpha[4]*fUpwind[16]+fUpwind[4]*alpha[16])+0.25*alpha[5]*fUpwind[14]+0.223606797749979*(alpha[8]*fUpwind[9]+fUpwind[8]*alpha[9]); 
  Ghat[42] += (0.223606797749979*alpha[19]+0.2500000000000001*alpha[2])*fUpwind[47]+0.223606797749979*alpha[26]*fUpwind[45]+0.2*alpha[16]*fUpwind[44]+0.25*alpha[5]*fUpwind[43]+(0.223606797749979*alpha[11]+0.25*alpha[0])*fUpwind[42]+(0.223606797749979*alpha[32]+0.25*alpha[7])*fUpwind[41]+0.223606797749979*alpha[36]*fUpwind[38]+0.2*alpha[8]*fUpwind[37]+fUpwind[31]*(0.2*alpha[35]+0.223606797749979*alpha[9])+0.2500000000000001*(alpha[1]*fUpwind[30]+alpha[15]*fUpwind[29])+(0.223606797749979*alpha[21]+0.2500000000000001*alpha[3])*fUpwind[28]+0.2*fUpwind[17]*alpha[25]+0.223606797749979*(alpha[16]*fUpwind[18]+alpha[4]*fUpwind[17])+0.25*alpha[6]*fUpwind[14]+0.223606797749979*alpha[8]*fUpwind[10]; 
  Ghat[43] += (0.223606797749979*alpha[20]+0.2500000000000001*alpha[1])*fUpwind[47]+0.2*alpha[16]*fUpwind[45]+0.223606797749979*(alpha[25]*fUpwind[44]+alpha[12]*fUpwind[43])+0.25*(alpha[0]*fUpwind[43]+alpha[5]*fUpwind[42])+(0.223606797749979*alpha[33]+0.25*alpha[6])*fUpwind[41]+0.2*alpha[9]*fUpwind[38]+0.223606797749979*alpha[35]*fUpwind[37]+fUpwind[31]*(0.2*alpha[36]+0.223606797749979*alpha[8])+0.2500000000000001*alpha[2]*fUpwind[30]+0.223606797749979*alpha[22]*fUpwind[29]+0.2500000000000001*(alpha[3]*fUpwind[29]+alpha[15]*fUpwind[28])+0.2*fUpwind[18]*alpha[26]+0.223606797749979*(alpha[4]*fUpwind[18]+alpha[16]*fUpwind[17])+0.25*alpha[7]*fUpwind[14]+0.223606797749979*alpha[9]*fUpwind[10]; 
  Ghat[44] += (0.1788854381999831*alpha[36]+0.2*alpha[8])*fUpwind[47]+0.1788854381999831*alpha[33]*fUpwind[46]+0.2*(alpha[6]*fUpwind[46]+alpha[5]*fUpwind[45])+(0.223606797749979*alpha[12]+0.159719141249985*alpha[11]+0.25*alpha[0])*fUpwind[44]+0.223606797749979*alpha[25]*fUpwind[43]+0.2*alpha[16]*fUpwind[42]+0.223606797749979*alpha[21]*fUpwind[40]+0.2*alpha[15]*fUpwind[39]+0.223606797749979*alpha[19]*fUpwind[38]+(0.159719141249985*alpha[19]+0.2500000000000001*alpha[2])*fUpwind[37]+0.2*(alpha[15]*fUpwind[36]+fUpwind[15]*alpha[36])+(0.223606797749979*alpha[22]+0.159719141249985*alpha[21]+0.2500000000000001*alpha[3])*fUpwind[35]+(0.223606797749979*(fUpwind[30]+fUpwind[22])+0.159719141249985*fUpwind[21]+0.2500000000000001*fUpwind[3])*alpha[35]+0.2*(alpha[16]*fUpwind[33]+fUpwind[16]*alpha[33])+(0.223606797749979*alpha[26]+0.159719141249985*alpha[25]+0.2500000000000001*alpha[4])*fUpwind[32]+(0.223606797749979*(fUpwind[27]+fUpwind[26])+0.159719141249985*fUpwind[25]+0.2500000000000001*fUpwind[4])*alpha[32]+(0.2*alpha[20]+0.223606797749979*alpha[1])*fUpwind[31]+0.25*(alpha[7]*fUpwind[25]+fUpwind[7]*alpha[25]+alpha[9]*fUpwind[21]+fUpwind[9]*alpha[21]+fUpwind[10]*alpha[19])+0.2500000000000001*alpha[11]*fUpwind[18]+0.223606797749979*(alpha[5]*fUpwind[17]+alpha[6]*fUpwind[16]+fUpwind[6]*alpha[16]+alpha[8]*fUpwind[15]+fUpwind[8]*alpha[15]); 
  Ghat[45] += (0.1788854381999831*alpha[35]+0.2*alpha[9])*fUpwind[47]+(0.1788854381999831*alpha[32]+0.2*alpha[7])*fUpwind[46]+(0.159719141249985*alpha[12]+0.223606797749979*alpha[11]+0.25*alpha[0])*fUpwind[45]+0.2*(alpha[5]*fUpwind[44]+alpha[16]*fUpwind[43])+0.223606797749979*alpha[26]*fUpwind[42]+0.2*alpha[15]*fUpwind[40]+0.223606797749979*alpha[22]*fUpwind[39]+(0.159719141249985*alpha[20]+0.2500000000000001*alpha[1])*fUpwind[38]+0.223606797749979*alpha[20]*fUpwind[37]+(0.159719141249985*alpha[22]+0.223606797749979*alpha[21]+0.2500000000000001*alpha[3])*fUpwind[36]+(0.223606797749979*fUpwind[30]+0.159719141249985*fUpwind[22]+0.223606797749979*fUpwind[21]+0.2500000000000001*fUpwind[3])*alpha[36]+0.2*(alpha[15]*fUpwind[35]+fUpwind[15]*alpha[35])+(0.159719141249985*alpha[26]+0.223606797749979*alpha[25]+0.2500000000000001*alpha[4])*fUpwind[33]+(0.223606797749979*fUpwind[27]+0.159719141249985*fUpwind[26]+0.223606797749979*fUpwind[25]+0.2500000000000001*fUpwind[4])*alpha[33]+0.2*(alpha[16]*fUpwind[32]+fUpwind[16]*alpha[32])+(0.2*alpha[19]+0.223606797749979*alpha[2])*fUpwind[31]+0.25*(alpha[6]*fUpwind[26]+fUpwind[6]*alpha[26]+alpha[8]*fUpwind[22]+fUpwind[8]*alpha[22]+fUpwind[10]*alpha[20])+0.223606797749979*alpha[5]*fUpwind[18]+0.2500000000000001*alpha[12]*fUpwind[17]+0.223606797749979*(alpha[7]*fUpwind[16]+fUpwind[7]*alpha[16]+alpha[9]*fUpwind[15]+fUpwind[9]*alpha[15]); 
  Ghat[46] += (0.223606797749979*(alpha[12]+alpha[11])+0.25*alpha[0])*fUpwind[46]+(0.1788854381999831*alpha[32]+0.2*alpha[7])*fUpwind[45]+(0.1788854381999831*alpha[33]+0.2*alpha[6])*fUpwind[44]+(0.223606797749979*alpha[20]+0.2500000000000001*alpha[1])*fUpwind[40]+(0.223606797749979*alpha[19]+0.2500000000000001*alpha[2])*fUpwind[39]+0.2*alpha[15]*(fUpwind[38]+fUpwind[37])+0.223606797749979*(fUpwind[24]*alpha[36]+fUpwind[23]*alpha[35])+(0.223606797749979*(alpha[26]+alpha[25])+0.2500000000000001*alpha[4])*fUpwind[34]+0.2*(fUpwind[18]*alpha[33]+fUpwind[17]*alpha[32])+(0.2*(alpha[22]+alpha[21])+0.223606797749979*alpha[3])*fUpwind[31]+0.25*(alpha[5]*fUpwind[27]+alpha[8]*fUpwind[24]+alpha[9]*fUpwind[23])+0.223606797749979*(alpha[6]*fUpwind[18]+alpha[7]*fUpwind[17])+0.2500000000000001*fUpwind[13]*alpha[16]+0.223606797749979*fUpwind[10]*alpha[15]; 
  Ghat[47] += (0.223606797749979*(alpha[12]+alpha[11])+0.25*alpha[0])*fUpwind[47]+(0.1788854381999831*alpha[35]+0.2*alpha[9])*fUpwind[45]+(0.1788854381999831*alpha[36]+0.2*alpha[8])*fUpwind[44]+(0.223606797749979*alpha[20]+0.2500000000000001*alpha[1])*fUpwind[43]+(0.223606797749979*alpha[19]+0.2500000000000001*alpha[2])*fUpwind[42]+(0.223606797749979*(alpha[22]+alpha[21])+0.2500000000000001*alpha[3])*fUpwind[41]+0.2*(alpha[16]*(fUpwind[38]+fUpwind[37])+fUpwind[18]*alpha[36]+fUpwind[17]*alpha[35])+0.223606797749979*(fUpwind[29]*alpha[33]+fUpwind[28]*alpha[32])+(0.2*(alpha[26]+alpha[25])+0.223606797749979*alpha[4])*fUpwind[31]+0.25*(alpha[5]*fUpwind[30]+alpha[6]*fUpwind[29]+alpha[7]*fUpwind[28])+0.223606797749979*(alpha[8]*fUpwind[18]+alpha[9]*fUpwind[17]+fUpwind[10]*alpha[16])+0.2500000000000001*fUpwind[14]*alpha[15]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv10; 
  out[1] += 0.7071067811865475*Ghat[1]*dv10; 
  out[2] += 0.7071067811865475*Ghat[2]*dv10; 
  out[3] += -1.224744871391589*Ghat[0]*dv10; 
  out[4] += 0.7071067811865475*Ghat[3]*dv10; 
  out[5] += 0.7071067811865475*Ghat[4]*dv10; 
  out[6] += 0.7071067811865475*Ghat[5]*dv10; 
  out[7] += -1.224744871391589*Ghat[1]*dv10; 
  out[8] += -1.224744871391589*Ghat[2]*dv10; 
  out[9] += 0.7071067811865475*Ghat[6]*dv10; 
  out[10] += 0.7071067811865475*Ghat[7]*dv10; 
  out[11] += -1.224744871391589*Ghat[3]*dv10; 
  out[12] += 0.7071067811865475*Ghat[8]*dv10; 
  out[13] += 0.7071067811865475*Ghat[9]*dv10; 
  out[14] += -1.224744871391589*Ghat[4]*dv10; 
  out[15] += 0.7071067811865475*Ghat[10]*dv10; 
  out[16] += 0.7071067811865475*Ghat[11]*dv10; 
  out[17] += 0.7071067811865475*Ghat[12]*dv10; 
  out[18] += 1.58113883008419*Ghat[0]*dv10; 
  out[19] += 0.7071067811865475*Ghat[13]*dv10; 
  out[20] += 0.7071067811865475*Ghat[14]*dv10; 
  out[21] += -1.224744871391589*Ghat[5]*dv10; 
  out[22] += 0.7071067811865475*Ghat[15]*dv10; 
  out[23] += -1.224744871391589*Ghat[6]*dv10; 
  out[24] += -1.224744871391589*Ghat[7]*dv10; 
  out[25] += 0.7071067811865475*Ghat[16]*dv10; 
  out[26] += -1.224744871391589*Ghat[8]*dv10; 
  out[27] += -1.224744871391589*Ghat[9]*dv10; 
  out[28] += 0.7071067811865475*Ghat[17]*dv10; 
  out[29] += 0.7071067811865475*Ghat[18]*dv10; 
  out[30] += -1.224744871391589*Ghat[10]*dv10; 
  out[31] += 0.7071067811865475*Ghat[19]*dv10; 
  out[32] += 0.7071067811865475*Ghat[20]*dv10; 
  out[33] += -1.224744871391589*Ghat[11]*dv10; 
  out[34] += -1.224744871391589*Ghat[12]*dv10; 
  out[35] += 1.58113883008419*Ghat[1]*dv10; 
  out[36] += 1.58113883008419*Ghat[2]*dv10; 
  out[37] += 0.7071067811865475*Ghat[21]*dv10; 
  out[38] += 0.7071067811865475*Ghat[22]*dv10; 
  out[39] += 1.58113883008419*Ghat[3]*dv10; 
  out[40] += 0.7071067811865475*Ghat[23]*dv10; 
  out[41] += 0.7071067811865475*Ghat[24]*dv10; 
  out[42] += -1.224744871391589*Ghat[13]*dv10; 
  out[43] += 0.7071067811865475*Ghat[25]*dv10; 
  out[44] += 0.7071067811865475*Ghat[26]*dv10; 
  out[45] += 1.58113883008419*Ghat[4]*dv10; 
  out[46] += 0.7071067811865475*Ghat[27]*dv10; 
  out[47] += 0.7071067811865475*Ghat[28]*dv10; 
  out[48] += 0.7071067811865475*Ghat[29]*dv10; 
  out[49] += -1.224744871391589*Ghat[14]*dv10; 
  out[50] += 0.7071067811865475*Ghat[30]*dv10; 
  out[51] += -1.224744871391589*Ghat[15]*dv10; 
  out[52] += -1.224744871391589*Ghat[16]*dv10; 
  out[53] += 0.7071067811865475*Ghat[31]*dv10; 
  out[54] += -1.224744871391589*Ghat[17]*dv10; 
  out[55] += -1.224744871391589*Ghat[18]*dv10; 
  out[56] += -1.224744871391589*Ghat[19]*dv10; 
  out[57] += -1.224744871391589*Ghat[20]*dv10; 
  out[58] += 1.58113883008419*Ghat[5]*dv10; 
  out[59] += 0.7071067811865475*Ghat[32]*dv10; 
  out[60] += 0.7071067811865475*Ghat[33]*dv10; 
  out[61] += -1.224744871391589*Ghat[21]*dv10; 
  out[62] += -1.224744871391589*Ghat[22]*dv10; 
  out[63] += 1.58113883008419*Ghat[6]*dv10; 
  out[64] += 1.58113883008419*Ghat[7]*dv10; 
  out[65] += 0.7071067811865475*Ghat[34]*dv10; 
  out[66] += -1.224744871391589*Ghat[23]*dv10; 
  out[67] += -1.224744871391589*Ghat[24]*dv10; 
  out[68] += 0.7071067811865475*Ghat[35]*dv10; 
  out[69] += 0.7071067811865475*Ghat[36]*dv10; 
  out[70] += -1.224744871391589*Ghat[25]*dv10; 
  out[71] += -1.224744871391589*Ghat[26]*dv10; 
  out[72] += 1.58113883008419*Ghat[8]*dv10; 
  out[73] += 1.58113883008419*Ghat[9]*dv10; 
  out[74] += 0.7071067811865475*Ghat[37]*dv10; 
  out[75] += 0.7071067811865475*Ghat[38]*dv10; 
  out[76] += 1.58113883008419*Ghat[10]*dv10; 
  out[77] += 0.7071067811865475*Ghat[39]*dv10; 
  out[78] += 0.7071067811865475*Ghat[40]*dv10; 
  out[79] += -1.224744871391589*Ghat[27]*dv10; 
  out[80] += 0.7071067811865475*Ghat[41]*dv10; 
  out[81] += -1.224744871391589*Ghat[28]*dv10; 
  out[82] += -1.224744871391589*Ghat[29]*dv10; 
  out[83] += 0.7071067811865475*Ghat[42]*dv10; 
  out[84] += 0.7071067811865475*Ghat[43]*dv10; 
  out[85] += -1.224744871391589*Ghat[30]*dv10; 
  out[86] += -1.224744871391589*Ghat[31]*dv10; 
  out[87] += -1.224744871391589*Ghat[32]*dv10; 
  out[88] += -1.224744871391589*Ghat[33]*dv10; 
  out[89] += 1.58113883008419*Ghat[15]*dv10; 
  out[90] += -1.224744871391589*Ghat[34]*dv10; 
  out[91] += -1.224744871391589*Ghat[35]*dv10; 
  out[92] += -1.224744871391589*Ghat[36]*dv10; 
  out[93] += 1.58113883008419*Ghat[16]*dv10; 
  out[94] += 0.7071067811865475*Ghat[44]*dv10; 
  out[95] += 0.7071067811865475*Ghat[45]*dv10; 
  out[96] += -1.224744871391589*Ghat[37]*dv10; 
  out[97] += -1.224744871391589*Ghat[38]*dv10; 
  out[98] += 1.58113883008419*Ghat[17]*dv10; 
  out[99] += 1.58113883008419*Ghat[18]*dv10; 
  out[100] += 0.7071067811865475*Ghat[46]*dv10; 
  out[101] += -1.224744871391589*Ghat[39]*dv10; 
  out[102] += -1.224744871391589*Ghat[40]*dv10; 
  out[103] += -1.224744871391589*Ghat[41]*dv10; 
  out[104] += 0.7071067811865475*Ghat[47]*dv10; 
  out[105] += -1.224744871391589*Ghat[42]*dv10; 
  out[106] += -1.224744871391589*Ghat[43]*dv10; 
  out[107] += -1.224744871391589*Ghat[44]*dv10; 
  out[108] += -1.224744871391589*Ghat[45]*dv10; 
  out[109] += 1.58113883008419*Ghat[31]*dv10; 
  out[110] += -1.224744871391589*Ghat[46]*dv10; 
  out[111] += -1.224744871391589*Ghat[47]*dv10; 

  } 
} 
