#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_5x_p2_surfx3_quad.h> 
#include <gkyl_basis_ser_5x_p2_upwind.h> 
GKYL_CU_DH void vlasov_surfvx_2x3v_ser_p2(const double *w, const double *dxv, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // qmem:      q/m*EM fields.
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 
  const double dv10 = 2/dxv[2]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double dv3 = dxv[4], wv3 = w[4]; 
  const double *E0 = &qmem[0]; 
  const double *B0 = &qmem[24]; 
  const double *B1 = &qmem[32]; 
  const double *B2 = &qmem[40]; 

  double alpha[48] = {0.0}; 

  alpha[0] = (-2.0*B1[0]*wv3)+2.0*B2[0]*wv2+2.0*E0[0]; 
  alpha[1] = (-2.0*B1[1]*wv3)+2.0*B2[1]*wv2+2.0*E0[1]; 
  alpha[2] = (-2.0*B1[2]*wv3)+2.0*B2[2]*wv2+2.0*E0[2]; 
  alpha[3] = 0.5773502691896258*B2[0]*dv2; 
  alpha[4] = -0.5773502691896258*B1[0]*dv3; 
  alpha[5] = (-2.0*B1[3]*wv3)+2.0*B2[3]*wv2+2.0*E0[3]; 
  alpha[6] = 0.5773502691896258*B2[1]*dv2; 
  alpha[7] = 0.5773502691896258*B2[2]*dv2; 
  alpha[8] = -0.5773502691896258*B1[1]*dv3; 
  alpha[9] = -0.5773502691896258*B1[2]*dv3; 
  alpha[11] = (-2.0*B1[4]*wv3)+2.0*B2[4]*wv2+2.0*E0[4]; 
  alpha[12] = (-2.0*B1[5]*wv3)+2.0*B2[5]*wv2+2.0*E0[5]; 
  alpha[15] = 0.5773502691896258*B2[3]*dv2; 
  alpha[16] = -0.5773502691896258*B1[3]*dv3; 
  alpha[19] = (-2.0*B1[6]*wv3)+2.0*B2[6]*wv2+2.0*E0[6]; 
  alpha[20] = (-2.0*B1[7]*wv3)+2.0*B2[7]*wv2+2.0*E0[7]; 
  alpha[21] = 0.5773502691896257*B2[4]*dv2; 
  alpha[22] = 0.5773502691896257*B2[5]*dv2; 
  alpha[25] = -0.5773502691896257*B1[4]*dv3; 
  alpha[26] = -0.5773502691896257*B1[5]*dv3; 
  alpha[32] = 0.5773502691896257*B2[6]*dv2; 
  alpha[33] = 0.5773502691896257*B2[7]*dv2; 
  alpha[35] = -0.5773502691896257*B1[6]*dv3; 
  alpha[36] = -0.5773502691896257*B1[7]*dv3; 

  double fUpwindQuad_l[81] = {0.0};
  double fUpwindQuad_r[81] = {0.0};
  double fUpwind_l[48] = {0.0};;
  double fUpwind_r[48] = {0.0};
  double Ghat_l[48] = {0.0}; 
  double Ghat_r[48] = {0.0}; 

  if (0.4024922359499621*(alpha[36]+alpha[35]+alpha[33]+alpha[32])-0.2999999999999999*(alpha[26]+alpha[25]+alpha[22]+alpha[21])-0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]-0.6037383539249431*(alpha[16]+alpha[15])+0.223606797749979*(alpha[12]+alpha[11])+0.45*(alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5])-0.3354101966249685*(alpha[4]+alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[0] = ser_5x_p2_surfx3_quad_0_r(fl); 
    fUpwindQuad_r[0] = ser_5x_p2_surfx3_quad_0_r(fc); 
  } else { 
    fUpwindQuad_l[0] = ser_5x_p2_surfx3_quad_0_l(fc); 
    fUpwindQuad_r[0] = ser_5x_p2_surfx3_quad_0_l(fr); 
  } 
  if ((-0.5031152949374527*(alpha[35]+alpha[32]))-0.2999999999999999*alpha[26]+0.375*alpha[25]-0.2999999999999999*alpha[22]+0.375*(alpha[21]+alpha[19])+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]+0.45*(alpha[9]+alpha[7])-0.3354101966249685*(alpha[4]+alpha[3]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[1] = ser_5x_p2_surfx3_quad_1_r(fl); 
    fUpwindQuad_r[1] = ser_5x_p2_surfx3_quad_1_r(fc); 
  } else { 
    fUpwindQuad_l[1] = ser_5x_p2_surfx3_quad_1_l(fc); 
    fUpwindQuad_r[1] = ser_5x_p2_surfx3_quad_1_l(fr); 
  } 
  if ((-0.4024922359499621*alpha[36])+0.4024922359499621*alpha[35]-0.4024922359499621*alpha[33]+0.4024922359499621*alpha[32]-0.2999999999999999*(alpha[26]+alpha[25]+alpha[22]+alpha[21])+0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]+0.6037383539249431*(alpha[16]+alpha[15])+0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[9]-0.45*alpha[8]+0.45*alpha[7]-0.45*(alpha[6]+alpha[5])-0.3354101966249685*(alpha[4]+alpha[3]+alpha[2])+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[2] = ser_5x_p2_surfx3_quad_2_r(fl); 
    fUpwindQuad_r[2] = ser_5x_p2_surfx3_quad_2_r(fc); 
  } else { 
    fUpwindQuad_l[2] = ser_5x_p2_surfx3_quad_2_l(fc); 
    fUpwindQuad_r[2] = ser_5x_p2_surfx3_quad_2_l(fr); 
  } 
  if ((-0.5031152949374527*(alpha[36]+alpha[33]))+0.375*alpha[26]-0.2999999999999999*alpha[25]+0.375*alpha[22]-0.2999999999999999*alpha[21]+0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]+0.45*(alpha[8]+alpha[6])-0.3354101966249685*(alpha[4]+alpha[3]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[3] = ser_5x_p2_surfx3_quad_3_r(fl); 
    fUpwindQuad_r[3] = ser_5x_p2_surfx3_quad_3_r(fc); 
  } else { 
    fUpwindQuad_l[3] = ser_5x_p2_surfx3_quad_3_l(fc); 
    fUpwindQuad_r[3] = ser_5x_p2_surfx3_quad_3_l(fr); 
  } 
  if (0.375*(alpha[26]+alpha[25]+alpha[22]+alpha[21])-0.2795084971874737*(alpha[12]+alpha[11])-0.3354101966249685*(alpha[4]+alpha[3])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[4] = ser_5x_p2_surfx3_quad_4_r(fl); 
    fUpwindQuad_r[4] = ser_5x_p2_surfx3_quad_4_r(fc); 
  } else { 
    fUpwindQuad_l[4] = ser_5x_p2_surfx3_quad_4_l(fc); 
    fUpwindQuad_r[4] = ser_5x_p2_surfx3_quad_4_l(fr); 
  } 
  if (0.5031152949374527*(alpha[36]+alpha[33])+0.375*alpha[26]-0.2999999999999999*alpha[25]+0.375*alpha[22]-0.2999999999999999*alpha[21]-0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]-0.45*(alpha[8]+alpha[6])-0.3354101966249685*(alpha[4]+alpha[3])+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[5] = ser_5x_p2_surfx3_quad_5_r(fl); 
    fUpwindQuad_r[5] = ser_5x_p2_surfx3_quad_5_r(fc); 
  } else { 
    fUpwindQuad_l[5] = ser_5x_p2_surfx3_quad_5_l(fc); 
    fUpwindQuad_r[5] = ser_5x_p2_surfx3_quad_5_l(fr); 
  } 
  if (0.4024922359499621*alpha[36]-0.4024922359499621*alpha[35]+0.4024922359499621*alpha[33]-0.4024922359499621*alpha[32]-0.2999999999999999*(alpha[26]+alpha[25]+alpha[22]+alpha[21])-0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]+0.6037383539249431*(alpha[16]+alpha[15])+0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[9]+0.45*alpha[8]-0.45*alpha[7]+0.45*alpha[6]-0.45*alpha[5]-0.3354101966249685*(alpha[4]+alpha[3])+0.3354101966249685*alpha[2]-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[6] = ser_5x_p2_surfx3_quad_6_r(fl); 
    fUpwindQuad_r[6] = ser_5x_p2_surfx3_quad_6_r(fc); 
  } else { 
    fUpwindQuad_l[6] = ser_5x_p2_surfx3_quad_6_l(fc); 
    fUpwindQuad_r[6] = ser_5x_p2_surfx3_quad_6_l(fr); 
  } 
  if (0.5031152949374527*(alpha[35]+alpha[32])-0.2999999999999999*alpha[26]+0.375*alpha[25]-0.2999999999999999*alpha[22]+0.375*alpha[21]-0.375*alpha[19]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]-0.45*(alpha[9]+alpha[7])-0.3354101966249685*(alpha[4]+alpha[3])+0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[7] = ser_5x_p2_surfx3_quad_7_r(fl); 
    fUpwindQuad_r[7] = ser_5x_p2_surfx3_quad_7_r(fc); 
  } else { 
    fUpwindQuad_l[7] = ser_5x_p2_surfx3_quad_7_l(fc); 
    fUpwindQuad_r[7] = ser_5x_p2_surfx3_quad_7_l(fr); 
  } 
  if ((-0.4024922359499621*(alpha[36]+alpha[35]+alpha[33]+alpha[32]))-0.2999999999999999*(alpha[26]+alpha[25]+alpha[22]+alpha[21])+0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]-0.6037383539249431*(alpha[16]+alpha[15])+0.223606797749979*(alpha[12]+alpha[11])-0.45*(alpha[9]+alpha[8]+alpha[7]+alpha[6])+0.45*alpha[5]-0.3354101966249685*(alpha[4]+alpha[3])+0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[8] = ser_5x_p2_surfx3_quad_8_r(fl); 
    fUpwindQuad_r[8] = ser_5x_p2_surfx3_quad_8_r(fc); 
  } else { 
    fUpwindQuad_l[8] = ser_5x_p2_surfx3_quad_8_l(fc); 
    fUpwindQuad_r[8] = ser_5x_p2_surfx3_quad_8_l(fr); 
  } 
  if (0.4024922359499621*(alpha[36]+alpha[35])-0.2999999999999999*(alpha[26]+alpha[25])-0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]-0.6037383539249431*alpha[16]+0.223606797749979*(alpha[12]+alpha[11])+0.45*(alpha[9]+alpha[8]+alpha[5])-0.3354101966249685*(alpha[4]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[9] = ser_5x_p2_surfx3_quad_9_r(fl); 
    fUpwindQuad_r[9] = ser_5x_p2_surfx3_quad_9_r(fc); 
  } else { 
    fUpwindQuad_l[9] = ser_5x_p2_surfx3_quad_9_l(fc); 
    fUpwindQuad_r[9] = ser_5x_p2_surfx3_quad_9_l(fr); 
  } 
  if ((-0.5031152949374527*alpha[35])-0.2999999999999999*alpha[26]+0.375*(alpha[25]+alpha[19])+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]+0.45*alpha[9]-0.3354101966249685*(alpha[4]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[10] = ser_5x_p2_surfx3_quad_10_r(fl); 
    fUpwindQuad_r[10] = ser_5x_p2_surfx3_quad_10_r(fc); 
  } else { 
    fUpwindQuad_l[10] = ser_5x_p2_surfx3_quad_10_l(fc); 
    fUpwindQuad_r[10] = ser_5x_p2_surfx3_quad_10_l(fr); 
  } 
  if ((-0.4024922359499621*alpha[36])+0.4024922359499621*alpha[35]-0.2999999999999999*(alpha[26]+alpha[25])+0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]+0.6037383539249431*alpha[16]+0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[9]-0.45*(alpha[8]+alpha[5])-0.3354101966249685*(alpha[4]+alpha[2])+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[11] = ser_5x_p2_surfx3_quad_11_r(fl); 
    fUpwindQuad_r[11] = ser_5x_p2_surfx3_quad_11_r(fc); 
  } else { 
    fUpwindQuad_l[11] = ser_5x_p2_surfx3_quad_11_l(fc); 
    fUpwindQuad_r[11] = ser_5x_p2_surfx3_quad_11_l(fr); 
  } 
  if ((-0.5031152949374527*alpha[36])+0.375*alpha[26]-0.2999999999999999*alpha[25]+0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]+0.45*alpha[8]-0.3354101966249685*(alpha[4]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[12] = ser_5x_p2_surfx3_quad_12_r(fl); 
    fUpwindQuad_r[12] = ser_5x_p2_surfx3_quad_12_r(fc); 
  } else { 
    fUpwindQuad_l[12] = ser_5x_p2_surfx3_quad_12_l(fc); 
    fUpwindQuad_r[12] = ser_5x_p2_surfx3_quad_12_l(fr); 
  } 
  if (0.375*(alpha[26]+alpha[25])-0.2795084971874737*(alpha[12]+alpha[11])-0.3354101966249685*alpha[4]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[13] = ser_5x_p2_surfx3_quad_13_r(fl); 
    fUpwindQuad_r[13] = ser_5x_p2_surfx3_quad_13_r(fc); 
  } else { 
    fUpwindQuad_l[13] = ser_5x_p2_surfx3_quad_13_l(fc); 
    fUpwindQuad_r[13] = ser_5x_p2_surfx3_quad_13_l(fr); 
  } 
  if (0.5031152949374527*alpha[36]+0.375*alpha[26]-0.2999999999999999*alpha[25]-0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]-0.45*alpha[8]-0.3354101966249685*alpha[4]+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[14] = ser_5x_p2_surfx3_quad_14_r(fl); 
    fUpwindQuad_r[14] = ser_5x_p2_surfx3_quad_14_r(fc); 
  } else { 
    fUpwindQuad_l[14] = ser_5x_p2_surfx3_quad_14_l(fc); 
    fUpwindQuad_r[14] = ser_5x_p2_surfx3_quad_14_l(fr); 
  } 
  if (0.4024922359499621*alpha[36]-0.4024922359499621*alpha[35]-0.2999999999999999*(alpha[26]+alpha[25])-0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]+0.6037383539249431*alpha[16]+0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[9]+0.45*alpha[8]-0.45*alpha[5]-0.3354101966249685*alpha[4]+0.3354101966249685*alpha[2]-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[15] = ser_5x_p2_surfx3_quad_15_r(fl); 
    fUpwindQuad_r[15] = ser_5x_p2_surfx3_quad_15_r(fc); 
  } else { 
    fUpwindQuad_l[15] = ser_5x_p2_surfx3_quad_15_l(fc); 
    fUpwindQuad_r[15] = ser_5x_p2_surfx3_quad_15_l(fr); 
  } 
  if (0.5031152949374527*alpha[35]-0.2999999999999999*alpha[26]+0.375*alpha[25]-0.375*alpha[19]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]-0.45*alpha[9]-0.3354101966249685*alpha[4]+0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[16] = ser_5x_p2_surfx3_quad_16_r(fl); 
    fUpwindQuad_r[16] = ser_5x_p2_surfx3_quad_16_r(fc); 
  } else { 
    fUpwindQuad_l[16] = ser_5x_p2_surfx3_quad_16_l(fc); 
    fUpwindQuad_r[16] = ser_5x_p2_surfx3_quad_16_l(fr); 
  } 
  if ((-0.4024922359499621*(alpha[36]+alpha[35]))-0.2999999999999999*(alpha[26]+alpha[25])+0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]-0.6037383539249431*alpha[16]+0.223606797749979*(alpha[12]+alpha[11])-0.45*(alpha[9]+alpha[8])+0.45*alpha[5]-0.3354101966249685*alpha[4]+0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[17] = ser_5x_p2_surfx3_quad_17_r(fl); 
    fUpwindQuad_r[17] = ser_5x_p2_surfx3_quad_17_r(fc); 
  } else { 
    fUpwindQuad_l[17] = ser_5x_p2_surfx3_quad_17_l(fc); 
    fUpwindQuad_r[17] = ser_5x_p2_surfx3_quad_17_l(fr); 
  } 
  if (0.4024922359499621*(alpha[36]+alpha[35])-0.4024922359499621*(alpha[33]+alpha[32])-0.2999999999999999*(alpha[26]+alpha[25])+0.2999999999999999*(alpha[22]+alpha[21])-0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]-0.6037383539249431*alpha[16]+0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])+0.45*(alpha[9]+alpha[8])-0.45*(alpha[7]+alpha[6])+0.45*alpha[5]-0.3354101966249685*alpha[4]+0.3354101966249685*alpha[3]-0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[18] = ser_5x_p2_surfx3_quad_18_r(fl); 
    fUpwindQuad_r[18] = ser_5x_p2_surfx3_quad_18_r(fc); 
  } else { 
    fUpwindQuad_l[18] = ser_5x_p2_surfx3_quad_18_l(fc); 
    fUpwindQuad_r[18] = ser_5x_p2_surfx3_quad_18_l(fr); 
  } 
  if ((-0.5031152949374527*alpha[35])+0.5031152949374527*alpha[32]-0.2999999999999999*alpha[26]+0.375*alpha[25]+0.2999999999999999*alpha[22]-0.375*alpha[21]+0.375*alpha[19]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]+0.45*alpha[9]-0.45*alpha[7]-0.3354101966249685*alpha[4]+0.3354101966249685*alpha[3]-0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[19] = ser_5x_p2_surfx3_quad_19_r(fl); 
    fUpwindQuad_r[19] = ser_5x_p2_surfx3_quad_19_r(fc); 
  } else { 
    fUpwindQuad_l[19] = ser_5x_p2_surfx3_quad_19_l(fc); 
    fUpwindQuad_r[19] = ser_5x_p2_surfx3_quad_19_l(fr); 
  } 
  if ((-0.4024922359499621*alpha[36])+0.4024922359499621*(alpha[35]+alpha[33])-0.4024922359499621*alpha[32]-0.2999999999999999*(alpha[26]+alpha[25])+0.2999999999999999*(alpha[22]+alpha[21])+0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]+0.6037383539249431*alpha[16]-0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[9]-0.45*(alpha[8]+alpha[7])+0.45*alpha[6]-0.45*alpha[5]-0.3354101966249685*alpha[4]+0.3354101966249685*alpha[3]-0.3354101966249685*alpha[2]+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[20] = ser_5x_p2_surfx3_quad_20_r(fl); 
    fUpwindQuad_r[20] = ser_5x_p2_surfx3_quad_20_r(fc); 
  } else { 
    fUpwindQuad_l[20] = ser_5x_p2_surfx3_quad_20_l(fc); 
    fUpwindQuad_r[20] = ser_5x_p2_surfx3_quad_20_l(fr); 
  } 
  if ((-0.5031152949374527*alpha[36])+0.5031152949374527*alpha[33]+0.375*alpha[26]-0.2999999999999999*alpha[25]-0.375*alpha[22]+0.2999999999999999*alpha[21]+0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]+0.45*alpha[8]-0.45*alpha[6]-0.3354101966249685*alpha[4]+0.3354101966249685*alpha[3]-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[21] = ser_5x_p2_surfx3_quad_21_r(fl); 
    fUpwindQuad_r[21] = ser_5x_p2_surfx3_quad_21_r(fc); 
  } else { 
    fUpwindQuad_l[21] = ser_5x_p2_surfx3_quad_21_l(fc); 
    fUpwindQuad_r[21] = ser_5x_p2_surfx3_quad_21_l(fr); 
  } 
  if (0.375*(alpha[26]+alpha[25])-0.375*(alpha[22]+alpha[21])-0.2795084971874737*(alpha[12]+alpha[11])-0.3354101966249685*alpha[4]+0.3354101966249685*alpha[3]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[22] = ser_5x_p2_surfx3_quad_22_r(fl); 
    fUpwindQuad_r[22] = ser_5x_p2_surfx3_quad_22_r(fc); 
  } else { 
    fUpwindQuad_l[22] = ser_5x_p2_surfx3_quad_22_l(fc); 
    fUpwindQuad_r[22] = ser_5x_p2_surfx3_quad_22_l(fr); 
  } 
  if (0.5031152949374527*alpha[36]-0.5031152949374527*alpha[33]+0.375*alpha[26]-0.2999999999999999*alpha[25]-0.375*alpha[22]+0.2999999999999999*alpha[21]-0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]-0.45*alpha[8]+0.45*alpha[6]-0.3354101966249685*alpha[4]+0.3354101966249685*(alpha[3]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[23] = ser_5x_p2_surfx3_quad_23_r(fl); 
    fUpwindQuad_r[23] = ser_5x_p2_surfx3_quad_23_r(fc); 
  } else { 
    fUpwindQuad_l[23] = ser_5x_p2_surfx3_quad_23_l(fc); 
    fUpwindQuad_r[23] = ser_5x_p2_surfx3_quad_23_l(fr); 
  } 
  if (0.4024922359499621*alpha[36]-0.4024922359499621*(alpha[35]+alpha[33])+0.4024922359499621*alpha[32]-0.2999999999999999*(alpha[26]+alpha[25])+0.2999999999999999*(alpha[22]+alpha[21])-0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]+0.6037383539249431*alpha[16]-0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[9]+0.45*(alpha[8]+alpha[7])-0.45*(alpha[6]+alpha[5])-0.3354101966249685*alpha[4]+0.3354101966249685*(alpha[3]+alpha[2])-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[24] = ser_5x_p2_surfx3_quad_24_r(fl); 
    fUpwindQuad_r[24] = ser_5x_p2_surfx3_quad_24_r(fc); 
  } else { 
    fUpwindQuad_l[24] = ser_5x_p2_surfx3_quad_24_l(fc); 
    fUpwindQuad_r[24] = ser_5x_p2_surfx3_quad_24_l(fr); 
  } 
  if (0.5031152949374527*alpha[35]-0.5031152949374527*alpha[32]-0.2999999999999999*alpha[26]+0.375*alpha[25]+0.2999999999999999*alpha[22]-0.375*(alpha[21]+alpha[19])+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]-0.45*alpha[9]+0.45*alpha[7]-0.3354101966249685*alpha[4]+0.3354101966249685*(alpha[3]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[25] = ser_5x_p2_surfx3_quad_25_r(fl); 
    fUpwindQuad_r[25] = ser_5x_p2_surfx3_quad_25_r(fc); 
  } else { 
    fUpwindQuad_l[25] = ser_5x_p2_surfx3_quad_25_l(fc); 
    fUpwindQuad_r[25] = ser_5x_p2_surfx3_quad_25_l(fr); 
  } 
  if ((-0.4024922359499621*(alpha[36]+alpha[35]))+0.4024922359499621*(alpha[33]+alpha[32])-0.2999999999999999*(alpha[26]+alpha[25])+0.2999999999999999*(alpha[22]+alpha[21])+0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]-0.6037383539249431*alpha[16]+0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])-0.45*(alpha[9]+alpha[8])+0.45*(alpha[7]+alpha[6]+alpha[5])-0.3354101966249685*alpha[4]+0.3354101966249685*(alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[26] = ser_5x_p2_surfx3_quad_26_r(fl); 
    fUpwindQuad_r[26] = ser_5x_p2_surfx3_quad_26_r(fc); 
  } else { 
    fUpwindQuad_l[26] = ser_5x_p2_surfx3_quad_26_l(fc); 
    fUpwindQuad_r[26] = ser_5x_p2_surfx3_quad_26_l(fr); 
  } 
  if (0.4024922359499621*(alpha[33]+alpha[32])-0.2999999999999999*(alpha[22]+alpha[21])-0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]-0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])+0.45*(alpha[7]+alpha[6]+alpha[5])-0.3354101966249685*(alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[27] = ser_5x_p2_surfx3_quad_27_r(fl); 
    fUpwindQuad_r[27] = ser_5x_p2_surfx3_quad_27_r(fc); 
  } else { 
    fUpwindQuad_l[27] = ser_5x_p2_surfx3_quad_27_l(fc); 
    fUpwindQuad_r[27] = ser_5x_p2_surfx3_quad_27_l(fr); 
  } 
  if ((-0.5031152949374527*alpha[32])-0.2999999999999999*alpha[22]+0.375*(alpha[21]+alpha[19])+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]+0.45*alpha[7]-0.3354101966249685*(alpha[3]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[28] = ser_5x_p2_surfx3_quad_28_r(fl); 
    fUpwindQuad_r[28] = ser_5x_p2_surfx3_quad_28_r(fc); 
  } else { 
    fUpwindQuad_l[28] = ser_5x_p2_surfx3_quad_28_l(fc); 
    fUpwindQuad_r[28] = ser_5x_p2_surfx3_quad_28_l(fr); 
  } 
  if ((-0.4024922359499621*alpha[33])+0.4024922359499621*alpha[32]-0.2999999999999999*(alpha[22]+alpha[21])+0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]+0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[7]-0.45*(alpha[6]+alpha[5])-0.3354101966249685*(alpha[3]+alpha[2])+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[29] = ser_5x_p2_surfx3_quad_29_r(fl); 
    fUpwindQuad_r[29] = ser_5x_p2_surfx3_quad_29_r(fc); 
  } else { 
    fUpwindQuad_l[29] = ser_5x_p2_surfx3_quad_29_l(fc); 
    fUpwindQuad_r[29] = ser_5x_p2_surfx3_quad_29_l(fr); 
  } 
  if ((-0.5031152949374527*alpha[33])+0.375*alpha[22]-0.2999999999999999*alpha[21]+0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]+0.45*alpha[6]-0.3354101966249685*(alpha[3]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[30] = ser_5x_p2_surfx3_quad_30_r(fl); 
    fUpwindQuad_r[30] = ser_5x_p2_surfx3_quad_30_r(fc); 
  } else { 
    fUpwindQuad_l[30] = ser_5x_p2_surfx3_quad_30_l(fc); 
    fUpwindQuad_r[30] = ser_5x_p2_surfx3_quad_30_l(fr); 
  } 
  if (0.375*(alpha[22]+alpha[21])-0.2795084971874737*(alpha[12]+alpha[11])-0.3354101966249685*alpha[3]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[31] = ser_5x_p2_surfx3_quad_31_r(fl); 
    fUpwindQuad_r[31] = ser_5x_p2_surfx3_quad_31_r(fc); 
  } else { 
    fUpwindQuad_l[31] = ser_5x_p2_surfx3_quad_31_l(fc); 
    fUpwindQuad_r[31] = ser_5x_p2_surfx3_quad_31_l(fr); 
  } 
  if (0.5031152949374527*alpha[33]+0.375*alpha[22]-0.2999999999999999*alpha[21]-0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]-0.45*alpha[6]-0.3354101966249685*alpha[3]+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[32] = ser_5x_p2_surfx3_quad_32_r(fl); 
    fUpwindQuad_r[32] = ser_5x_p2_surfx3_quad_32_r(fc); 
  } else { 
    fUpwindQuad_l[32] = ser_5x_p2_surfx3_quad_32_l(fc); 
    fUpwindQuad_r[32] = ser_5x_p2_surfx3_quad_32_l(fr); 
  } 
  if (0.4024922359499621*alpha[33]-0.4024922359499621*alpha[32]-0.2999999999999999*(alpha[22]+alpha[21])-0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]+0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[7]+0.45*alpha[6]-0.45*alpha[5]-0.3354101966249685*alpha[3]+0.3354101966249685*alpha[2]-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[33] = ser_5x_p2_surfx3_quad_33_r(fl); 
    fUpwindQuad_r[33] = ser_5x_p2_surfx3_quad_33_r(fc); 
  } else { 
    fUpwindQuad_l[33] = ser_5x_p2_surfx3_quad_33_l(fc); 
    fUpwindQuad_r[33] = ser_5x_p2_surfx3_quad_33_l(fr); 
  } 
  if (0.5031152949374527*alpha[32]-0.2999999999999999*alpha[22]+0.375*alpha[21]-0.375*alpha[19]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]-0.45*alpha[7]-0.3354101966249685*alpha[3]+0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[34] = ser_5x_p2_surfx3_quad_34_r(fl); 
    fUpwindQuad_r[34] = ser_5x_p2_surfx3_quad_34_r(fc); 
  } else { 
    fUpwindQuad_l[34] = ser_5x_p2_surfx3_quad_34_l(fc); 
    fUpwindQuad_r[34] = ser_5x_p2_surfx3_quad_34_l(fr); 
  } 
  if ((-0.4024922359499621*(alpha[33]+alpha[32]))-0.2999999999999999*(alpha[22]+alpha[21])+0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]-0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])-0.45*(alpha[7]+alpha[6])+0.45*alpha[5]-0.3354101966249685*alpha[3]+0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[35] = ser_5x_p2_surfx3_quad_35_r(fl); 
    fUpwindQuad_r[35] = ser_5x_p2_surfx3_quad_35_r(fc); 
  } else { 
    fUpwindQuad_l[35] = ser_5x_p2_surfx3_quad_35_l(fc); 
    fUpwindQuad_r[35] = ser_5x_p2_surfx3_quad_35_l(fr); 
  } 
  if ((-0.2999999999999998*alpha[20])-0.2999999999999999*alpha[19]+0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[5]-0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[36] = ser_5x_p2_surfx3_quad_36_r(fl); 
    fUpwindQuad_r[36] = ser_5x_p2_surfx3_quad_36_r(fc); 
  } else { 
    fUpwindQuad_l[36] = ser_5x_p2_surfx3_quad_36_l(fc); 
    fUpwindQuad_r[36] = ser_5x_p2_surfx3_quad_36_l(fr); 
  } 
  if (0.375*alpha[19]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]-0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[37] = ser_5x_p2_surfx3_quad_37_r(fl); 
    fUpwindQuad_r[37] = ser_5x_p2_surfx3_quad_37_r(fc); 
  } else { 
    fUpwindQuad_l[37] = ser_5x_p2_surfx3_quad_37_l(fc); 
    fUpwindQuad_r[37] = ser_5x_p2_surfx3_quad_37_l(fr); 
  } 
  if (0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]+0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[5]-0.3354101966249685*alpha[2]+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[38] = ser_5x_p2_surfx3_quad_38_r(fl); 
    fUpwindQuad_r[38] = ser_5x_p2_surfx3_quad_38_r(fc); 
  } else { 
    fUpwindQuad_l[38] = ser_5x_p2_surfx3_quad_38_l(fc); 
    fUpwindQuad_r[38] = ser_5x_p2_surfx3_quad_38_l(fr); 
  } 
  if (0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[39] = ser_5x_p2_surfx3_quad_39_r(fl); 
    fUpwindQuad_r[39] = ser_5x_p2_surfx3_quad_39_r(fc); 
  } else { 
    fUpwindQuad_l[39] = ser_5x_p2_surfx3_quad_39_l(fc); 
    fUpwindQuad_r[39] = ser_5x_p2_surfx3_quad_39_l(fr); 
  } 
  if (0.25*alpha[0]-0.2795084971874737*(alpha[12]+alpha[11]) > 0) { 
    fUpwindQuad_l[40] = ser_5x_p2_surfx3_quad_40_r(fl); 
    fUpwindQuad_r[40] = ser_5x_p2_surfx3_quad_40_r(fc); 
  } else { 
    fUpwindQuad_l[40] = ser_5x_p2_surfx3_quad_40_l(fc); 
    fUpwindQuad_r[40] = ser_5x_p2_surfx3_quad_40_l(fr); 
  } 
  if ((-0.375*alpha[20])-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[41] = ser_5x_p2_surfx3_quad_41_r(fl); 
    fUpwindQuad_r[41] = ser_5x_p2_surfx3_quad_41_r(fc); 
  } else { 
    fUpwindQuad_l[41] = ser_5x_p2_surfx3_quad_41_l(fc); 
    fUpwindQuad_r[41] = ser_5x_p2_surfx3_quad_41_l(fr); 
  } 
  if ((-0.2999999999999998*alpha[20])+0.2999999999999999*alpha[19]+0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[5]+0.3354101966249685*alpha[2]-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[42] = ser_5x_p2_surfx3_quad_42_r(fl); 
    fUpwindQuad_r[42] = ser_5x_p2_surfx3_quad_42_r(fc); 
  } else { 
    fUpwindQuad_l[42] = ser_5x_p2_surfx3_quad_42_l(fc); 
    fUpwindQuad_r[42] = ser_5x_p2_surfx3_quad_42_l(fr); 
  } 
  if ((-0.375*alpha[19])+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]+0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[43] = ser_5x_p2_surfx3_quad_43_r(fl); 
    fUpwindQuad_r[43] = ser_5x_p2_surfx3_quad_43_r(fc); 
  } else { 
    fUpwindQuad_l[43] = ser_5x_p2_surfx3_quad_43_l(fc); 
    fUpwindQuad_r[43] = ser_5x_p2_surfx3_quad_43_l(fr); 
  } 
  if (0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]+0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[5]+0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[44] = ser_5x_p2_surfx3_quad_44_r(fl); 
    fUpwindQuad_r[44] = ser_5x_p2_surfx3_quad_44_r(fc); 
  } else { 
    fUpwindQuad_l[44] = ser_5x_p2_surfx3_quad_44_l(fc); 
    fUpwindQuad_r[44] = ser_5x_p2_surfx3_quad_44_l(fr); 
  } 
  if ((-0.4024922359499621*(alpha[33]+alpha[32]))+0.2999999999999999*(alpha[22]+alpha[21])-0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]+0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])-0.45*(alpha[7]+alpha[6])+0.45*alpha[5]+0.3354101966249685*alpha[3]-0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[45] = ser_5x_p2_surfx3_quad_45_r(fl); 
    fUpwindQuad_r[45] = ser_5x_p2_surfx3_quad_45_r(fc); 
  } else { 
    fUpwindQuad_l[45] = ser_5x_p2_surfx3_quad_45_l(fc); 
    fUpwindQuad_r[45] = ser_5x_p2_surfx3_quad_45_l(fr); 
  } 
  if (0.5031152949374527*alpha[32]+0.2999999999999999*alpha[22]-0.375*alpha[21]+0.375*alpha[19]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]-0.45*alpha[7]+0.3354101966249685*alpha[3]-0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[46] = ser_5x_p2_surfx3_quad_46_r(fl); 
    fUpwindQuad_r[46] = ser_5x_p2_surfx3_quad_46_r(fc); 
  } else { 
    fUpwindQuad_l[46] = ser_5x_p2_surfx3_quad_46_l(fc); 
    fUpwindQuad_r[46] = ser_5x_p2_surfx3_quad_46_l(fr); 
  } 
  if (0.4024922359499621*alpha[33]-0.4024922359499621*alpha[32]+0.2999999999999999*(alpha[22]+alpha[21])+0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]-0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[7]+0.45*alpha[6]-0.45*alpha[5]+0.3354101966249685*alpha[3]-0.3354101966249685*alpha[2]+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[47] = ser_5x_p2_surfx3_quad_47_r(fl); 
    fUpwindQuad_r[47] = ser_5x_p2_surfx3_quad_47_r(fc); 
  } else { 
    fUpwindQuad_l[47] = ser_5x_p2_surfx3_quad_47_l(fc); 
    fUpwindQuad_r[47] = ser_5x_p2_surfx3_quad_47_l(fr); 
  } 
  if (0.5031152949374527*alpha[33]-0.375*alpha[22]+0.2999999999999999*alpha[21]+0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]-0.45*alpha[6]+0.3354101966249685*alpha[3]-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[48] = ser_5x_p2_surfx3_quad_48_r(fl); 
    fUpwindQuad_r[48] = ser_5x_p2_surfx3_quad_48_r(fc); 
  } else { 
    fUpwindQuad_l[48] = ser_5x_p2_surfx3_quad_48_l(fc); 
    fUpwindQuad_r[48] = ser_5x_p2_surfx3_quad_48_l(fr); 
  } 
  if ((-0.375*(alpha[22]+alpha[21]))-0.2795084971874737*(alpha[12]+alpha[11])+0.3354101966249685*alpha[3]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[49] = ser_5x_p2_surfx3_quad_49_r(fl); 
    fUpwindQuad_r[49] = ser_5x_p2_surfx3_quad_49_r(fc); 
  } else { 
    fUpwindQuad_l[49] = ser_5x_p2_surfx3_quad_49_l(fc); 
    fUpwindQuad_r[49] = ser_5x_p2_surfx3_quad_49_l(fr); 
  } 
  if ((-0.5031152949374527*alpha[33])-0.375*alpha[22]+0.2999999999999999*alpha[21]-0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]+0.45*alpha[6]+0.3354101966249685*(alpha[3]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[50] = ser_5x_p2_surfx3_quad_50_r(fl); 
    fUpwindQuad_r[50] = ser_5x_p2_surfx3_quad_50_r(fc); 
  } else { 
    fUpwindQuad_l[50] = ser_5x_p2_surfx3_quad_50_l(fc); 
    fUpwindQuad_r[50] = ser_5x_p2_surfx3_quad_50_l(fr); 
  } 
  if ((-0.4024922359499621*alpha[33])+0.4024922359499621*alpha[32]+0.2999999999999999*(alpha[22]+alpha[21])-0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]-0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[7]-0.45*(alpha[6]+alpha[5])+0.3354101966249685*(alpha[3]+alpha[2])-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[51] = ser_5x_p2_surfx3_quad_51_r(fl); 
    fUpwindQuad_r[51] = ser_5x_p2_surfx3_quad_51_r(fc); 
  } else { 
    fUpwindQuad_l[51] = ser_5x_p2_surfx3_quad_51_l(fc); 
    fUpwindQuad_r[51] = ser_5x_p2_surfx3_quad_51_l(fr); 
  } 
  if ((-0.5031152949374527*alpha[32])+0.2999999999999999*alpha[22]-0.375*(alpha[21]+alpha[19])+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]+0.45*alpha[7]+0.3354101966249685*(alpha[3]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[52] = ser_5x_p2_surfx3_quad_52_r(fl); 
    fUpwindQuad_r[52] = ser_5x_p2_surfx3_quad_52_r(fc); 
  } else { 
    fUpwindQuad_l[52] = ser_5x_p2_surfx3_quad_52_l(fc); 
    fUpwindQuad_r[52] = ser_5x_p2_surfx3_quad_52_l(fr); 
  } 
  if (0.4024922359499621*(alpha[33]+alpha[32])+0.2999999999999999*(alpha[22]+alpha[21])+0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]+0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])+0.45*(alpha[7]+alpha[6]+alpha[5])+0.3354101966249685*(alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[53] = ser_5x_p2_surfx3_quad_53_r(fl); 
    fUpwindQuad_r[53] = ser_5x_p2_surfx3_quad_53_r(fc); 
  } else { 
    fUpwindQuad_l[53] = ser_5x_p2_surfx3_quad_53_l(fc); 
    fUpwindQuad_r[53] = ser_5x_p2_surfx3_quad_53_l(fr); 
  } 
  if ((-0.4024922359499621*(alpha[36]+alpha[35]))+0.4024922359499621*(alpha[33]+alpha[32])+0.2999999999999999*(alpha[26]+alpha[25])-0.2999999999999999*(alpha[22]+alpha[21])-0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]+0.6037383539249431*alpha[16]-0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])-0.45*(alpha[9]+alpha[8])+0.45*(alpha[7]+alpha[6]+alpha[5])+0.3354101966249685*alpha[4]-0.3354101966249685*(alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[54] = ser_5x_p2_surfx3_quad_54_r(fl); 
    fUpwindQuad_r[54] = ser_5x_p2_surfx3_quad_54_r(fc); 
  } else { 
    fUpwindQuad_l[54] = ser_5x_p2_surfx3_quad_54_l(fc); 
    fUpwindQuad_r[54] = ser_5x_p2_surfx3_quad_54_l(fr); 
  } 
  if (0.5031152949374527*alpha[35]-0.5031152949374527*alpha[32]+0.2999999999999999*alpha[26]-0.375*alpha[25]-0.2999999999999999*alpha[22]+0.375*(alpha[21]+alpha[19])+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]-0.45*alpha[9]+0.45*alpha[7]+0.3354101966249685*alpha[4]-0.3354101966249685*(alpha[3]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[55] = ser_5x_p2_surfx3_quad_55_r(fl); 
    fUpwindQuad_r[55] = ser_5x_p2_surfx3_quad_55_r(fc); 
  } else { 
    fUpwindQuad_l[55] = ser_5x_p2_surfx3_quad_55_l(fc); 
    fUpwindQuad_r[55] = ser_5x_p2_surfx3_quad_55_l(fr); 
  } 
  if (0.4024922359499621*alpha[36]-0.4024922359499621*(alpha[35]+alpha[33])+0.4024922359499621*alpha[32]+0.2999999999999999*(alpha[26]+alpha[25])-0.2999999999999999*(alpha[22]+alpha[21])+0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]-0.6037383539249431*alpha[16]+0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[9]+0.45*(alpha[8]+alpha[7])-0.45*(alpha[6]+alpha[5])+0.3354101966249685*alpha[4]-0.3354101966249685*(alpha[3]+alpha[2])+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[56] = ser_5x_p2_surfx3_quad_56_r(fl); 
    fUpwindQuad_r[56] = ser_5x_p2_surfx3_quad_56_r(fc); 
  } else { 
    fUpwindQuad_l[56] = ser_5x_p2_surfx3_quad_56_l(fc); 
    fUpwindQuad_r[56] = ser_5x_p2_surfx3_quad_56_l(fr); 
  } 
  if (0.5031152949374527*alpha[36]-0.5031152949374527*alpha[33]-0.375*alpha[26]+0.2999999999999999*alpha[25]+0.375*alpha[22]-0.2999999999999999*alpha[21]+0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]-0.45*alpha[8]+0.45*alpha[6]+0.3354101966249685*alpha[4]-0.3354101966249685*(alpha[3]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[57] = ser_5x_p2_surfx3_quad_57_r(fl); 
    fUpwindQuad_r[57] = ser_5x_p2_surfx3_quad_57_r(fc); 
  } else { 
    fUpwindQuad_l[57] = ser_5x_p2_surfx3_quad_57_l(fc); 
    fUpwindQuad_r[57] = ser_5x_p2_surfx3_quad_57_l(fr); 
  } 
  if ((-0.375*(alpha[26]+alpha[25]))+0.375*(alpha[22]+alpha[21])-0.2795084971874737*(alpha[12]+alpha[11])+0.3354101966249685*alpha[4]-0.3354101966249685*alpha[3]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[58] = ser_5x_p2_surfx3_quad_58_r(fl); 
    fUpwindQuad_r[58] = ser_5x_p2_surfx3_quad_58_r(fc); 
  } else { 
    fUpwindQuad_l[58] = ser_5x_p2_surfx3_quad_58_l(fc); 
    fUpwindQuad_r[58] = ser_5x_p2_surfx3_quad_58_l(fr); 
  } 
  if ((-0.5031152949374527*alpha[36])+0.5031152949374527*alpha[33]-0.375*alpha[26]+0.2999999999999999*alpha[25]+0.375*alpha[22]-0.2999999999999999*alpha[21]-0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]+0.45*alpha[8]-0.45*alpha[6]+0.3354101966249685*alpha[4]-0.3354101966249685*alpha[3]+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[59] = ser_5x_p2_surfx3_quad_59_r(fl); 
    fUpwindQuad_r[59] = ser_5x_p2_surfx3_quad_59_r(fc); 
  } else { 
    fUpwindQuad_l[59] = ser_5x_p2_surfx3_quad_59_l(fc); 
    fUpwindQuad_r[59] = ser_5x_p2_surfx3_quad_59_l(fr); 
  } 
  if ((-0.4024922359499621*alpha[36])+0.4024922359499621*(alpha[35]+alpha[33])-0.4024922359499621*alpha[32]+0.2999999999999999*(alpha[26]+alpha[25])-0.2999999999999999*(alpha[22]+alpha[21])-0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]-0.6037383539249431*alpha[16]+0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[9]-0.45*(alpha[8]+alpha[7])+0.45*alpha[6]-0.45*alpha[5]+0.3354101966249685*alpha[4]-0.3354101966249685*alpha[3]+0.3354101966249685*alpha[2]-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[60] = ser_5x_p2_surfx3_quad_60_r(fl); 
    fUpwindQuad_r[60] = ser_5x_p2_surfx3_quad_60_r(fc); 
  } else { 
    fUpwindQuad_l[60] = ser_5x_p2_surfx3_quad_60_l(fc); 
    fUpwindQuad_r[60] = ser_5x_p2_surfx3_quad_60_l(fr); 
  } 
  if ((-0.5031152949374527*alpha[35])+0.5031152949374527*alpha[32]+0.2999999999999999*alpha[26]-0.375*alpha[25]-0.2999999999999999*alpha[22]+0.375*alpha[21]-0.375*alpha[19]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]+0.45*alpha[9]-0.45*alpha[7]+0.3354101966249685*alpha[4]-0.3354101966249685*alpha[3]+0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[61] = ser_5x_p2_surfx3_quad_61_r(fl); 
    fUpwindQuad_r[61] = ser_5x_p2_surfx3_quad_61_r(fc); 
  } else { 
    fUpwindQuad_l[61] = ser_5x_p2_surfx3_quad_61_l(fc); 
    fUpwindQuad_r[61] = ser_5x_p2_surfx3_quad_61_l(fr); 
  } 
  if (0.4024922359499621*(alpha[36]+alpha[35])-0.4024922359499621*(alpha[33]+alpha[32])+0.2999999999999999*(alpha[26]+alpha[25])-0.2999999999999999*(alpha[22]+alpha[21])+0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]+0.6037383539249431*alpha[16]-0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])+0.45*(alpha[9]+alpha[8])-0.45*(alpha[7]+alpha[6])+0.45*alpha[5]+0.3354101966249685*alpha[4]-0.3354101966249685*alpha[3]+0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[62] = ser_5x_p2_surfx3_quad_62_r(fl); 
    fUpwindQuad_r[62] = ser_5x_p2_surfx3_quad_62_r(fc); 
  } else { 
    fUpwindQuad_l[62] = ser_5x_p2_surfx3_quad_62_l(fc); 
    fUpwindQuad_r[62] = ser_5x_p2_surfx3_quad_62_l(fr); 
  } 
  if ((-0.4024922359499621*(alpha[36]+alpha[35]))+0.2999999999999999*(alpha[26]+alpha[25])-0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]+0.6037383539249431*alpha[16]+0.223606797749979*(alpha[12]+alpha[11])-0.45*(alpha[9]+alpha[8])+0.45*alpha[5]+0.3354101966249685*alpha[4]-0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[63] = ser_5x_p2_surfx3_quad_63_r(fl); 
    fUpwindQuad_r[63] = ser_5x_p2_surfx3_quad_63_r(fc); 
  } else { 
    fUpwindQuad_l[63] = ser_5x_p2_surfx3_quad_63_l(fc); 
    fUpwindQuad_r[63] = ser_5x_p2_surfx3_quad_63_l(fr); 
  } 
  if (0.5031152949374527*alpha[35]+0.2999999999999999*alpha[26]-0.375*alpha[25]+0.375*alpha[19]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]-0.45*alpha[9]+0.3354101966249685*alpha[4]-0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[64] = ser_5x_p2_surfx3_quad_64_r(fl); 
    fUpwindQuad_r[64] = ser_5x_p2_surfx3_quad_64_r(fc); 
  } else { 
    fUpwindQuad_l[64] = ser_5x_p2_surfx3_quad_64_l(fc); 
    fUpwindQuad_r[64] = ser_5x_p2_surfx3_quad_64_l(fr); 
  } 
  if (0.4024922359499621*alpha[36]-0.4024922359499621*alpha[35]+0.2999999999999999*(alpha[26]+alpha[25])+0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]-0.6037383539249431*alpha[16]+0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[9]+0.45*alpha[8]-0.45*alpha[5]+0.3354101966249685*alpha[4]-0.3354101966249685*alpha[2]+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[65] = ser_5x_p2_surfx3_quad_65_r(fl); 
    fUpwindQuad_r[65] = ser_5x_p2_surfx3_quad_65_r(fc); 
  } else { 
    fUpwindQuad_l[65] = ser_5x_p2_surfx3_quad_65_l(fc); 
    fUpwindQuad_r[65] = ser_5x_p2_surfx3_quad_65_l(fr); 
  } 
  if (0.5031152949374527*alpha[36]-0.375*alpha[26]+0.2999999999999999*alpha[25]+0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]-0.45*alpha[8]+0.3354101966249685*alpha[4]-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[66] = ser_5x_p2_surfx3_quad_66_r(fl); 
    fUpwindQuad_r[66] = ser_5x_p2_surfx3_quad_66_r(fc); 
  } else { 
    fUpwindQuad_l[66] = ser_5x_p2_surfx3_quad_66_l(fc); 
    fUpwindQuad_r[66] = ser_5x_p2_surfx3_quad_66_l(fr); 
  } 
  if ((-0.375*(alpha[26]+alpha[25]))-0.2795084971874737*(alpha[12]+alpha[11])+0.3354101966249685*alpha[4]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[67] = ser_5x_p2_surfx3_quad_67_r(fl); 
    fUpwindQuad_r[67] = ser_5x_p2_surfx3_quad_67_r(fc); 
  } else { 
    fUpwindQuad_l[67] = ser_5x_p2_surfx3_quad_67_l(fc); 
    fUpwindQuad_r[67] = ser_5x_p2_surfx3_quad_67_l(fr); 
  } 
  if ((-0.5031152949374527*alpha[36])-0.375*alpha[26]+0.2999999999999999*alpha[25]-0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]+0.45*alpha[8]+0.3354101966249685*(alpha[4]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[68] = ser_5x_p2_surfx3_quad_68_r(fl); 
    fUpwindQuad_r[68] = ser_5x_p2_surfx3_quad_68_r(fc); 
  } else { 
    fUpwindQuad_l[68] = ser_5x_p2_surfx3_quad_68_l(fc); 
    fUpwindQuad_r[68] = ser_5x_p2_surfx3_quad_68_l(fr); 
  } 
  if ((-0.4024922359499621*alpha[36])+0.4024922359499621*alpha[35]+0.2999999999999999*(alpha[26]+alpha[25])-0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]-0.6037383539249431*alpha[16]+0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[9]-0.45*(alpha[8]+alpha[5])+0.3354101966249685*(alpha[4]+alpha[2])-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[69] = ser_5x_p2_surfx3_quad_69_r(fl); 
    fUpwindQuad_r[69] = ser_5x_p2_surfx3_quad_69_r(fc); 
  } else { 
    fUpwindQuad_l[69] = ser_5x_p2_surfx3_quad_69_l(fc); 
    fUpwindQuad_r[69] = ser_5x_p2_surfx3_quad_69_l(fr); 
  } 
  if ((-0.5031152949374527*alpha[35])+0.2999999999999999*alpha[26]-0.375*(alpha[25]+alpha[19])+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]+0.45*alpha[9]+0.3354101966249685*(alpha[4]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[70] = ser_5x_p2_surfx3_quad_70_r(fl); 
    fUpwindQuad_r[70] = ser_5x_p2_surfx3_quad_70_r(fc); 
  } else { 
    fUpwindQuad_l[70] = ser_5x_p2_surfx3_quad_70_l(fc); 
    fUpwindQuad_r[70] = ser_5x_p2_surfx3_quad_70_l(fr); 
  } 
  if (0.4024922359499621*(alpha[36]+alpha[35])+0.2999999999999999*(alpha[26]+alpha[25])+0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]+0.6037383539249431*alpha[16]+0.223606797749979*(alpha[12]+alpha[11])+0.45*(alpha[9]+alpha[8]+alpha[5])+0.3354101966249685*(alpha[4]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[71] = ser_5x_p2_surfx3_quad_71_r(fl); 
    fUpwindQuad_r[71] = ser_5x_p2_surfx3_quad_71_r(fc); 
  } else { 
    fUpwindQuad_l[71] = ser_5x_p2_surfx3_quad_71_l(fc); 
    fUpwindQuad_r[71] = ser_5x_p2_surfx3_quad_71_l(fr); 
  } 
  if ((-0.4024922359499621*(alpha[36]+alpha[35]+alpha[33]+alpha[32]))+0.2999999999999999*(alpha[26]+alpha[25]+alpha[22]+alpha[21])-0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]+0.6037383539249431*(alpha[16]+alpha[15])+0.223606797749979*(alpha[12]+alpha[11])-0.45*(alpha[9]+alpha[8]+alpha[7]+alpha[6])+0.45*alpha[5]+0.3354101966249685*(alpha[4]+alpha[3])-0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[72] = ser_5x_p2_surfx3_quad_72_r(fl); 
    fUpwindQuad_r[72] = ser_5x_p2_surfx3_quad_72_r(fc); 
  } else { 
    fUpwindQuad_l[72] = ser_5x_p2_surfx3_quad_72_l(fc); 
    fUpwindQuad_r[72] = ser_5x_p2_surfx3_quad_72_l(fr); 
  } 
  if (0.5031152949374527*(alpha[35]+alpha[32])+0.2999999999999999*alpha[26]-0.375*alpha[25]+0.2999999999999999*alpha[22]-0.375*alpha[21]+0.375*alpha[19]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]-0.45*(alpha[9]+alpha[7])+0.3354101966249685*(alpha[4]+alpha[3])-0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[73] = ser_5x_p2_surfx3_quad_73_r(fl); 
    fUpwindQuad_r[73] = ser_5x_p2_surfx3_quad_73_r(fc); 
  } else { 
    fUpwindQuad_l[73] = ser_5x_p2_surfx3_quad_73_l(fc); 
    fUpwindQuad_r[73] = ser_5x_p2_surfx3_quad_73_l(fr); 
  } 
  if (0.4024922359499621*alpha[36]-0.4024922359499621*alpha[35]+0.4024922359499621*alpha[33]-0.4024922359499621*alpha[32]+0.2999999999999999*(alpha[26]+alpha[25]+alpha[22]+alpha[21])+0.2999999999999998*alpha[20]-0.2999999999999999*alpha[19]-0.6037383539249431*(alpha[16]+alpha[15])+0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[9]+0.45*alpha[8]-0.45*alpha[7]+0.45*alpha[6]-0.45*alpha[5]+0.3354101966249685*(alpha[4]+alpha[3])-0.3354101966249685*alpha[2]+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[74] = ser_5x_p2_surfx3_quad_74_r(fl); 
    fUpwindQuad_r[74] = ser_5x_p2_surfx3_quad_74_r(fc); 
  } else { 
    fUpwindQuad_l[74] = ser_5x_p2_surfx3_quad_74_l(fc); 
    fUpwindQuad_r[74] = ser_5x_p2_surfx3_quad_74_l(fr); 
  } 
  if (0.5031152949374527*(alpha[36]+alpha[33])-0.375*alpha[26]+0.2999999999999999*alpha[25]-0.375*alpha[22]+0.2999999999999999*alpha[21]+0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]-0.45*(alpha[8]+alpha[6])+0.3354101966249685*(alpha[4]+alpha[3])-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[75] = ser_5x_p2_surfx3_quad_75_r(fl); 
    fUpwindQuad_r[75] = ser_5x_p2_surfx3_quad_75_r(fc); 
  } else { 
    fUpwindQuad_l[75] = ser_5x_p2_surfx3_quad_75_l(fc); 
    fUpwindQuad_r[75] = ser_5x_p2_surfx3_quad_75_l(fr); 
  } 
  if ((-0.375*(alpha[26]+alpha[25]+alpha[22]+alpha[21]))-0.2795084971874737*(alpha[12]+alpha[11])+0.3354101966249685*(alpha[4]+alpha[3])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[76] = ser_5x_p2_surfx3_quad_76_r(fl); 
    fUpwindQuad_r[76] = ser_5x_p2_surfx3_quad_76_r(fc); 
  } else { 
    fUpwindQuad_l[76] = ser_5x_p2_surfx3_quad_76_l(fc); 
    fUpwindQuad_r[76] = ser_5x_p2_surfx3_quad_76_l(fr); 
  } 
  if ((-0.5031152949374527*(alpha[36]+alpha[33]))-0.375*alpha[26]+0.2999999999999999*alpha[25]-0.375*alpha[22]+0.2999999999999999*alpha[21]-0.375*alpha[20]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]+0.45*(alpha[8]+alpha[6])+0.3354101966249685*(alpha[4]+alpha[3]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[77] = ser_5x_p2_surfx3_quad_77_r(fl); 
    fUpwindQuad_r[77] = ser_5x_p2_surfx3_quad_77_r(fc); 
  } else { 
    fUpwindQuad_l[77] = ser_5x_p2_surfx3_quad_77_l(fc); 
    fUpwindQuad_r[77] = ser_5x_p2_surfx3_quad_77_l(fr); 
  } 
  if ((-0.4024922359499621*alpha[36])+0.4024922359499621*alpha[35]-0.4024922359499621*alpha[33]+0.4024922359499621*alpha[32]+0.2999999999999999*(alpha[26]+alpha[25]+alpha[22]+alpha[21])-0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]-0.6037383539249431*(alpha[16]+alpha[15])+0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[9]-0.45*alpha[8]+0.45*alpha[7]-0.45*(alpha[6]+alpha[5])+0.3354101966249685*(alpha[4]+alpha[3]+alpha[2])-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[78] = ser_5x_p2_surfx3_quad_78_r(fl); 
    fUpwindQuad_r[78] = ser_5x_p2_surfx3_quad_78_r(fc); 
  } else { 
    fUpwindQuad_l[78] = ser_5x_p2_surfx3_quad_78_l(fc); 
    fUpwindQuad_r[78] = ser_5x_p2_surfx3_quad_78_l(fr); 
  } 
  if ((-0.5031152949374527*(alpha[35]+alpha[32]))+0.2999999999999999*alpha[26]-0.375*alpha[25]+0.2999999999999999*alpha[22]-0.375*(alpha[21]+alpha[19])+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]+0.45*(alpha[9]+alpha[7])+0.3354101966249685*(alpha[4]+alpha[3]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[79] = ser_5x_p2_surfx3_quad_79_r(fl); 
    fUpwindQuad_r[79] = ser_5x_p2_surfx3_quad_79_r(fc); 
  } else { 
    fUpwindQuad_l[79] = ser_5x_p2_surfx3_quad_79_l(fc); 
    fUpwindQuad_r[79] = ser_5x_p2_surfx3_quad_79_l(fr); 
  } 
  if (0.4024922359499621*(alpha[36]+alpha[35]+alpha[33]+alpha[32])+0.2999999999999999*(alpha[26]+alpha[25]+alpha[22]+alpha[21])+0.2999999999999998*alpha[20]+0.2999999999999999*alpha[19]+0.6037383539249431*(alpha[16]+alpha[15])+0.223606797749979*(alpha[12]+alpha[11])+0.45*(alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5])+0.3354101966249685*(alpha[4]+alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[80] = ser_5x_p2_surfx3_quad_80_r(fl); 
    fUpwindQuad_r[80] = ser_5x_p2_surfx3_quad_80_r(fc); 
  } else { 
    fUpwindQuad_l[80] = ser_5x_p2_surfx3_quad_80_l(fc); 
    fUpwindQuad_r[80] = ser_5x_p2_surfx3_quad_80_l(fr); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_5x_p2_upwind(fUpwindQuad_l, fUpwind_l); 
  ser_5x_p2_upwind(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] += 0.25*(alpha[36]*fUpwind_l[36]+alpha[35]*fUpwind_l[35]+alpha[33]*fUpwind_l[33]+alpha[32]*fUpwind_l[32]+alpha[26]*fUpwind_l[26]+alpha[25]*fUpwind_l[25]+alpha[22]*fUpwind_l[22]+alpha[21]*fUpwind_l[21]+alpha[20]*fUpwind_l[20]+alpha[19]*fUpwind_l[19]+alpha[16]*fUpwind_l[16]+alpha[15]*fUpwind_l[15]+alpha[12]*fUpwind_l[12]+alpha[11]*fUpwind_l[11]+alpha[9]*fUpwind_l[9]+alpha[8]*fUpwind_l[8]+alpha[7]*fUpwind_l[7]+alpha[6]*fUpwind_l[6]+alpha[5]*fUpwind_l[5]+alpha[4]*fUpwind_l[4]+alpha[3]*fUpwind_l[3]+alpha[2]*fUpwind_l[2]+alpha[1]*fUpwind_l[1]+alpha[0]*fUpwind_l[0]); 
  Ghat_l[1] += 0.2500000000000001*(alpha[26]*fUpwind_l[36]+fUpwind_l[26]*alpha[36])+0.223606797749979*(alpha[16]*fUpwind_l[35]+fUpwind_l[16]*alpha[35])+0.2500000000000001*(alpha[22]*fUpwind_l[33]+fUpwind_l[22]*alpha[33])+0.223606797749979*(alpha[15]*fUpwind_l[32]+fUpwind_l[15]*alpha[32])+0.223606797749979*(alpha[8]*fUpwind_l[25]+fUpwind_l[8]*alpha[25]+alpha[6]*fUpwind_l[21]+fUpwind_l[6]*alpha[21])+0.2500000000000001*(alpha[12]*fUpwind_l[20]+fUpwind_l[12]*alpha[20])+0.223606797749979*(alpha[5]*fUpwind_l[19]+fUpwind_l[5]*alpha[19])+0.25*(alpha[9]*fUpwind_l[16]+fUpwind_l[9]*alpha[16]+alpha[7]*fUpwind_l[15]+fUpwind_l[7]*alpha[15])+0.223606797749979*(alpha[1]*fUpwind_l[11]+fUpwind_l[1]*alpha[11])+0.25*(alpha[4]*fUpwind_l[8]+fUpwind_l[4]*alpha[8]+alpha[3]*fUpwind_l[6]+fUpwind_l[3]*alpha[6]+alpha[2]*fUpwind_l[5]+fUpwind_l[2]*alpha[5]+alpha[0]*fUpwind_l[1]+fUpwind_l[0]*alpha[1]); 
  Ghat_l[2] += 0.223606797749979*(alpha[16]*fUpwind_l[36]+fUpwind_l[16]*alpha[36])+0.2500000000000001*(alpha[25]*fUpwind_l[35]+fUpwind_l[25]*alpha[35])+0.223606797749979*(alpha[15]*fUpwind_l[33]+fUpwind_l[15]*alpha[33])+0.2500000000000001*(alpha[21]*fUpwind_l[32]+fUpwind_l[21]*alpha[32])+0.223606797749979*(alpha[9]*fUpwind_l[26]+fUpwind_l[9]*alpha[26]+alpha[7]*fUpwind_l[22]+fUpwind_l[7]*alpha[22]+alpha[5]*fUpwind_l[20]+fUpwind_l[5]*alpha[20])+0.2500000000000001*(alpha[11]*fUpwind_l[19]+fUpwind_l[11]*alpha[19])+0.25*(alpha[8]*fUpwind_l[16]+fUpwind_l[8]*alpha[16]+alpha[6]*fUpwind_l[15]+fUpwind_l[6]*alpha[15])+0.223606797749979*(alpha[2]*fUpwind_l[12]+fUpwind_l[2]*alpha[12])+0.25*(alpha[4]*fUpwind_l[9]+fUpwind_l[4]*alpha[9]+alpha[3]*fUpwind_l[7]+fUpwind_l[3]*alpha[7]+alpha[1]*fUpwind_l[5]+fUpwind_l[1]*alpha[5]+alpha[0]*fUpwind_l[2]+fUpwind_l[0]*alpha[2]); 
  Ghat_l[3] += 0.2500000000000001*(alpha[36]*fUpwind_l[45]+alpha[35]*fUpwind_l[44]+alpha[26]*fUpwind_l[38]+alpha[25]*fUpwind_l[37])+0.223606797749979*alpha[15]*fUpwind_l[34]+0.2500000000000001*(alpha[20]*fUpwind_l[33]+fUpwind_l[20]*alpha[33]+alpha[19]*fUpwind_l[32]+fUpwind_l[19]*alpha[32])+0.25*alpha[16]*fUpwind_l[31]+0.223606797749979*(alpha[7]*fUpwind_l[24]+alpha[6]*fUpwind_l[23])+0.2500000000000001*(alpha[12]*fUpwind_l[22]+fUpwind_l[12]*alpha[22]+alpha[11]*fUpwind_l[21]+fUpwind_l[11]*alpha[21])+0.25*(alpha[9]*fUpwind_l[18]+alpha[8]*fUpwind_l[17]+alpha[5]*fUpwind_l[15]+fUpwind_l[5]*alpha[15])+0.223606797749979*alpha[3]*fUpwind_l[13]+0.25*(alpha[4]*fUpwind_l[10]+alpha[2]*fUpwind_l[7]+fUpwind_l[2]*alpha[7]+alpha[1]*fUpwind_l[6]+fUpwind_l[1]*alpha[6]+alpha[0]*fUpwind_l[3]+fUpwind_l[0]*alpha[3]); 
  Ghat_l[4] += 0.2500000000000001*(alpha[33]*fUpwind_l[45]+alpha[32]*fUpwind_l[44])+0.223606797749979*alpha[16]*fUpwind_l[41]+0.2500000000000001*(alpha[22]*fUpwind_l[38]+alpha[21]*fUpwind_l[37]+alpha[20]*fUpwind_l[36]+fUpwind_l[20]*alpha[36]+alpha[19]*fUpwind_l[35]+fUpwind_l[19]*alpha[35])+0.25*alpha[15]*fUpwind_l[31]+0.223606797749979*(alpha[9]*fUpwind_l[29]+alpha[8]*fUpwind_l[28])+0.2500000000000001*(alpha[12]*fUpwind_l[26]+fUpwind_l[12]*alpha[26]+alpha[11]*fUpwind_l[25]+fUpwind_l[11]*alpha[25])+0.25*(alpha[7]*fUpwind_l[18]+alpha[6]*fUpwind_l[17]+alpha[5]*fUpwind_l[16]+fUpwind_l[5]*alpha[16])+0.223606797749979*alpha[4]*fUpwind_l[14]+0.25*(alpha[3]*fUpwind_l[10]+alpha[2]*fUpwind_l[9]+fUpwind_l[2]*alpha[9]+alpha[1]*fUpwind_l[8]+fUpwind_l[1]*alpha[8]+alpha[0]*fUpwind_l[4]+fUpwind_l[0]*alpha[4]); 
  Ghat_l[5] += (0.2*alpha[35]+0.223606797749979*alpha[9])*fUpwind_l[36]+0.2*fUpwind_l[35]*alpha[36]+0.223606797749979*(fUpwind_l[9]*alpha[36]+alpha[8]*fUpwind_l[35]+fUpwind_l[8]*alpha[35])+(0.2*alpha[32]+0.223606797749979*alpha[7])*fUpwind_l[33]+0.2*fUpwind_l[32]*alpha[33]+0.223606797749979*(fUpwind_l[7]*alpha[33]+alpha[6]*fUpwind_l[32]+fUpwind_l[6]*alpha[32])+0.223606797749979*(alpha[16]*fUpwind_l[26]+fUpwind_l[16]*alpha[26]+alpha[16]*fUpwind_l[25]+fUpwind_l[16]*alpha[25]+alpha[15]*fUpwind_l[22]+fUpwind_l[15]*alpha[22]+alpha[15]*fUpwind_l[21]+fUpwind_l[15]*alpha[21])+(0.2*alpha[19]+0.223606797749979*alpha[2])*fUpwind_l[20]+0.2*fUpwind_l[19]*alpha[20]+0.223606797749979*(fUpwind_l[2]*alpha[20]+alpha[1]*fUpwind_l[19]+fUpwind_l[1]*alpha[19])+0.25*(alpha[4]*fUpwind_l[16]+fUpwind_l[4]*alpha[16]+alpha[3]*fUpwind_l[15]+fUpwind_l[3]*alpha[15])+0.223606797749979*(alpha[5]*fUpwind_l[12]+fUpwind_l[5]*alpha[12]+alpha[5]*fUpwind_l[11]+fUpwind_l[5]*alpha[11])+0.25*(alpha[8]*fUpwind_l[9]+fUpwind_l[8]*alpha[9]+alpha[6]*fUpwind_l[7]+fUpwind_l[6]*alpha[7]+alpha[0]*fUpwind_l[5]+fUpwind_l[0]*alpha[5]+alpha[1]*fUpwind_l[2]+fUpwind_l[1]*alpha[2]); 
  Ghat_l[6] += 0.25*alpha[26]*fUpwind_l[45]+0.223606797749979*alpha[16]*fUpwind_l[44]+0.25*alpha[36]*fUpwind_l[38]+0.223606797749979*(alpha[8]*fUpwind_l[37]+fUpwind_l[31]*alpha[35])+(0.2*alpha[32]+0.223606797749979*alpha[7])*fUpwind_l[34]+0.25*(alpha[12]*fUpwind_l[33]+fUpwind_l[12]*alpha[33])+0.223606797749979*(alpha[5]*fUpwind_l[32]+fUpwind_l[5]*alpha[32])+0.25*alpha[9]*fUpwind_l[31]+0.223606797749979*(fUpwind_l[17]*alpha[25]+alpha[15]*fUpwind_l[24])+(0.2*alpha[21]+0.223606797749979*alpha[3])*fUpwind_l[23]+0.25*(alpha[20]*fUpwind_l[22]+fUpwind_l[20]*alpha[22])+0.223606797749979*(alpha[1]*fUpwind_l[21]+fUpwind_l[1]*alpha[21]+alpha[15]*fUpwind_l[19]+fUpwind_l[15]*alpha[19])+0.25*(alpha[16]*fUpwind_l[18]+alpha[4]*fUpwind_l[17]+alpha[2]*fUpwind_l[15]+fUpwind_l[2]*alpha[15])+0.223606797749979*(alpha[6]*(fUpwind_l[13]+fUpwind_l[11])+fUpwind_l[6]*alpha[11])+0.25*(alpha[8]*fUpwind_l[10]+alpha[5]*fUpwind_l[7]+fUpwind_l[5]*alpha[7]+alpha[0]*fUpwind_l[6]+fUpwind_l[0]*alpha[6]+alpha[1]*fUpwind_l[3]+fUpwind_l[1]*alpha[3]); 
  Ghat_l[7] += 0.223606797749979*alpha[16]*fUpwind_l[45]+0.25*alpha[25]*fUpwind_l[44]+0.223606797749979*alpha[9]*fUpwind_l[38]+0.25*alpha[35]*fUpwind_l[37]+0.223606797749979*fUpwind_l[31]*alpha[36]+0.2*alpha[33]*fUpwind_l[34]+0.223606797749979*(alpha[6]*fUpwind_l[34]+alpha[5]*fUpwind_l[33]+fUpwind_l[5]*alpha[33])+0.25*(alpha[11]*fUpwind_l[32]+fUpwind_l[11]*alpha[32]+alpha[8]*fUpwind_l[31])+0.223606797749979*fUpwind_l[18]*alpha[26]+0.2*alpha[22]*fUpwind_l[24]+0.223606797749979*(alpha[3]*fUpwind_l[24]+alpha[15]*fUpwind_l[23]+alpha[2]*fUpwind_l[22]+fUpwind_l[2]*alpha[22])+0.25*(alpha[19]*fUpwind_l[21]+fUpwind_l[19]*alpha[21])+0.223606797749979*(alpha[15]*fUpwind_l[20]+fUpwind_l[15]*alpha[20])+0.25*(alpha[4]*fUpwind_l[18]+alpha[16]*fUpwind_l[17]+alpha[1]*fUpwind_l[15]+fUpwind_l[1]*alpha[15])+0.223606797749979*(alpha[7]*(fUpwind_l[13]+fUpwind_l[12])+fUpwind_l[7]*alpha[12])+0.25*(alpha[9]*fUpwind_l[10]+alpha[0]*fUpwind_l[7]+fUpwind_l[0]*alpha[7]+alpha[5]*fUpwind_l[6]+fUpwind_l[5]*alpha[6]+alpha[2]*fUpwind_l[3]+fUpwind_l[2]*alpha[3]); 
  Ghat_l[8] += 0.25*alpha[22]*fUpwind_l[45]+0.223606797749979*alpha[15]*fUpwind_l[44]+(0.2*alpha[35]+0.223606797749979*alpha[9])*fUpwind_l[41]+0.25*alpha[33]*fUpwind_l[38]+0.223606797749979*alpha[6]*fUpwind_l[37]+0.25*(alpha[12]*fUpwind_l[36]+fUpwind_l[12]*alpha[36])+0.223606797749979*(alpha[5]*fUpwind_l[35]+fUpwind_l[5]*alpha[35])+fUpwind_l[31]*(0.223606797749979*alpha[32]+0.25*alpha[7])+0.223606797749979*alpha[16]*fUpwind_l[29]+(0.2*alpha[25]+0.223606797749979*alpha[4])*fUpwind_l[28]+0.25*(alpha[20]*fUpwind_l[26]+fUpwind_l[20]*alpha[26])+0.223606797749979*(alpha[1]*fUpwind_l[25]+fUpwind_l[1]*alpha[25]+fUpwind_l[17]*alpha[21]+alpha[16]*fUpwind_l[19]+fUpwind_l[16]*alpha[19])+0.25*(alpha[15]*fUpwind_l[18]+alpha[3]*fUpwind_l[17]+alpha[2]*fUpwind_l[16]+fUpwind_l[2]*alpha[16])+0.223606797749979*(alpha[8]*(fUpwind_l[14]+fUpwind_l[11])+fUpwind_l[8]*alpha[11])+0.25*(alpha[6]*fUpwind_l[10]+alpha[5]*fUpwind_l[9]+fUpwind_l[5]*alpha[9]+alpha[0]*fUpwind_l[8]+fUpwind_l[0]*alpha[8]+alpha[1]*fUpwind_l[4]+fUpwind_l[1]*alpha[4]); 
  Ghat_l[9] += 0.223606797749979*alpha[15]*fUpwind_l[45]+0.25*alpha[21]*fUpwind_l[44]+0.2*alpha[36]*fUpwind_l[41]+0.223606797749979*(alpha[8]*fUpwind_l[41]+alpha[7]*fUpwind_l[38])+0.25*alpha[32]*fUpwind_l[37]+0.223606797749979*(alpha[5]*fUpwind_l[36]+fUpwind_l[5]*alpha[36])+0.25*(alpha[11]*fUpwind_l[35]+fUpwind_l[11]*alpha[35])+fUpwind_l[31]*(0.223606797749979*alpha[33]+0.25*alpha[6])+0.2*alpha[26]*fUpwind_l[29]+0.223606797749979*(alpha[4]*fUpwind_l[29]+alpha[16]*fUpwind_l[28]+alpha[2]*fUpwind_l[26]+fUpwind_l[2]*alpha[26])+0.25*(alpha[19]*fUpwind_l[25]+fUpwind_l[19]*alpha[25])+0.223606797749979*(fUpwind_l[18]*alpha[22]+alpha[16]*fUpwind_l[20]+fUpwind_l[16]*alpha[20])+0.25*(alpha[3]*fUpwind_l[18]+alpha[15]*fUpwind_l[17]+alpha[1]*fUpwind_l[16]+fUpwind_l[1]*alpha[16])+0.223606797749979*(alpha[9]*(fUpwind_l[14]+fUpwind_l[12])+fUpwind_l[9]*alpha[12])+0.25*(alpha[7]*fUpwind_l[10]+alpha[0]*fUpwind_l[9]+fUpwind_l[0]*alpha[9]+alpha[5]*fUpwind_l[8]+fUpwind_l[5]*alpha[8]+alpha[2]*fUpwind_l[4]+fUpwind_l[2]*alpha[4]); 
  Ghat_l[10] += 0.223606797749979*(alpha[16]*fUpwind_l[47]+alpha[15]*fUpwind_l[46])+0.25*(alpha[20]*fUpwind_l[45]+alpha[19]*fUpwind_l[44])+0.223606797749979*(alpha[9]*fUpwind_l[43]+alpha[8]*fUpwind_l[42]+alpha[7]*fUpwind_l[40]+alpha[6]*fUpwind_l[39])+0.25*(alpha[12]*fUpwind_l[38]+alpha[11]*fUpwind_l[37]+alpha[33]*fUpwind_l[36]+fUpwind_l[33]*alpha[36]+alpha[32]*fUpwind_l[35]+fUpwind_l[32]*alpha[35]+alpha[5]*fUpwind_l[31])+0.223606797749979*(alpha[4]*fUpwind_l[30]+alpha[3]*fUpwind_l[27])+0.25*(alpha[22]*fUpwind_l[26]+fUpwind_l[22]*alpha[26]+alpha[21]*fUpwind_l[25]+fUpwind_l[21]*alpha[25]+alpha[2]*fUpwind_l[18]+alpha[1]*fUpwind_l[17]+alpha[15]*fUpwind_l[16]+fUpwind_l[15]*alpha[16]+alpha[0]*fUpwind_l[10]+alpha[7]*fUpwind_l[9]+fUpwind_l[7]*alpha[9]+alpha[6]*fUpwind_l[8]+fUpwind_l[6]*alpha[8]+alpha[3]*fUpwind_l[4]+fUpwind_l[3]*alpha[4]); 
  Ghat_l[11] += 0.223606797749979*alpha[36]*fUpwind_l[36]+0.159719141249985*alpha[35]*fUpwind_l[35]+0.25*(alpha[9]*fUpwind_l[35]+fUpwind_l[9]*alpha[35])+0.223606797749979*alpha[33]*fUpwind_l[33]+0.159719141249985*alpha[32]*fUpwind_l[32]+0.25*(alpha[7]*fUpwind_l[32]+fUpwind_l[7]*alpha[32])+0.159719141249985*alpha[25]*fUpwind_l[25]+0.2500000000000001*(alpha[4]*fUpwind_l[25]+fUpwind_l[4]*alpha[25])+0.159719141249985*alpha[21]*fUpwind_l[21]+0.2500000000000001*(alpha[3]*fUpwind_l[21]+fUpwind_l[3]*alpha[21])+0.223606797749979*alpha[20]*fUpwind_l[20]+0.159719141249985*alpha[19]*fUpwind_l[19]+0.2500000000000001*(alpha[2]*fUpwind_l[19]+fUpwind_l[2]*alpha[19])+0.223606797749979*(alpha[16]*fUpwind_l[16]+alpha[15]*fUpwind_l[15])+0.159719141249985*alpha[11]*fUpwind_l[11]+0.25*(alpha[0]*fUpwind_l[11]+fUpwind_l[0]*alpha[11])+0.223606797749979*(alpha[8]*fUpwind_l[8]+alpha[6]*fUpwind_l[6]+alpha[5]*fUpwind_l[5]+alpha[1]*fUpwind_l[1]); 
  Ghat_l[12] += 0.159719141249985*alpha[36]*fUpwind_l[36]+0.25*(alpha[8]*fUpwind_l[36]+fUpwind_l[8]*alpha[36])+0.223606797749979*alpha[35]*fUpwind_l[35]+0.159719141249985*alpha[33]*fUpwind_l[33]+0.25*(alpha[6]*fUpwind_l[33]+fUpwind_l[6]*alpha[33])+0.223606797749979*alpha[32]*fUpwind_l[32]+0.159719141249985*alpha[26]*fUpwind_l[26]+0.2500000000000001*(alpha[4]*fUpwind_l[26]+fUpwind_l[4]*alpha[26])+0.159719141249985*alpha[22]*fUpwind_l[22]+0.2500000000000001*(alpha[3]*fUpwind_l[22]+fUpwind_l[3]*alpha[22])+0.159719141249985*alpha[20]*fUpwind_l[20]+0.2500000000000001*(alpha[1]*fUpwind_l[20]+fUpwind_l[1]*alpha[20])+0.223606797749979*(alpha[19]*fUpwind_l[19]+alpha[16]*fUpwind_l[16]+alpha[15]*fUpwind_l[15])+0.159719141249985*alpha[12]*fUpwind_l[12]+0.25*(alpha[0]*fUpwind_l[12]+fUpwind_l[0]*alpha[12])+0.223606797749979*(alpha[9]*fUpwind_l[9]+alpha[7]*fUpwind_l[7]+alpha[5]*fUpwind_l[5]+alpha[2]*fUpwind_l[2]); 
  Ghat_l[13] += 0.2500000000000001*alpha[16]*fUpwind_l[46]+0.25*(alpha[9]*fUpwind_l[40]+alpha[8]*fUpwind_l[39]+alpha[5]*fUpwind_l[34])+0.223606797749979*(alpha[33]*fUpwind_l[33]+alpha[32]*fUpwind_l[32])+0.2500000000000001*(alpha[4]*fUpwind_l[27]+alpha[2]*fUpwind_l[24]+alpha[1]*fUpwind_l[23])+0.223606797749979*(alpha[22]*fUpwind_l[22]+alpha[21]*fUpwind_l[21]+alpha[15]*fUpwind_l[15])+0.25*alpha[0]*fUpwind_l[13]+0.223606797749979*(alpha[7]*fUpwind_l[7]+alpha[6]*fUpwind_l[6]+alpha[3]*fUpwind_l[3]); 
  Ghat_l[14] += 0.2500000000000001*alpha[15]*fUpwind_l[47]+0.25*(alpha[7]*fUpwind_l[43]+alpha[6]*fUpwind_l[42]+alpha[5]*fUpwind_l[41])+0.223606797749979*(alpha[36]*fUpwind_l[36]+alpha[35]*fUpwind_l[35])+0.2500000000000001*(alpha[3]*fUpwind_l[30]+alpha[2]*fUpwind_l[29]+alpha[1]*fUpwind_l[28])+0.223606797749979*(alpha[26]*fUpwind_l[26]+alpha[25]*fUpwind_l[25]+alpha[16]*fUpwind_l[16])+0.25*alpha[0]*fUpwind_l[14]+0.223606797749979*(alpha[9]*fUpwind_l[9]+alpha[8]*fUpwind_l[8]+alpha[4]*fUpwind_l[4]); 
  Ghat_l[15] += (0.2*alpha[35]+0.223606797749979*alpha[9])*fUpwind_l[45]+(0.2*alpha[36]+0.223606797749979*alpha[8])*fUpwind_l[44]+0.223606797749979*(alpha[16]*(fUpwind_l[38]+fUpwind_l[37])+fUpwind_l[18]*alpha[36]+fUpwind_l[17]*alpha[35])+(0.2*(alpha[22]+alpha[21])+0.223606797749979*alpha[3])*fUpwind_l[34]+(0.2*alpha[19]+0.223606797749979*alpha[2])*fUpwind_l[33]+(0.2*(fUpwind_l[24]+fUpwind_l[19])+0.223606797749979*fUpwind_l[2])*alpha[33]+(0.2*alpha[20]+0.223606797749979*alpha[1])*fUpwind_l[32]+(0.2*(fUpwind_l[23]+fUpwind_l[20])+0.223606797749979*fUpwind_l[1])*alpha[32]+(0.223606797749979*(alpha[26]+alpha[25])+0.25*alpha[4])*fUpwind_l[31]+0.223606797749979*(alpha[6]*fUpwind_l[24]+alpha[7]*fUpwind_l[23]+alpha[5]*fUpwind_l[22]+fUpwind_l[5]*alpha[22]+alpha[5]*fUpwind_l[21]+fUpwind_l[5]*alpha[21]+alpha[7]*fUpwind_l[20]+fUpwind_l[7]*alpha[20]+alpha[6]*fUpwind_l[19]+fUpwind_l[6]*alpha[19])+0.25*(alpha[8]*fUpwind_l[18]+alpha[9]*fUpwind_l[17]+fUpwind_l[10]*alpha[16])+(0.223606797749979*(alpha[12]+alpha[11])+0.25*alpha[0])*fUpwind_l[15]+0.223606797749979*(fUpwind_l[13]+fUpwind_l[12]+fUpwind_l[11])*alpha[15]+0.25*(fUpwind_l[0]*alpha[15]+alpha[1]*fUpwind_l[7]+fUpwind_l[1]*alpha[7]+alpha[2]*fUpwind_l[6]+fUpwind_l[2]*alpha[6]+alpha[3]*fUpwind_l[5]+fUpwind_l[3]*alpha[5]); 
  Ghat_l[16] += (0.2*alpha[32]+0.223606797749979*alpha[7])*fUpwind_l[45]+(0.2*alpha[33]+0.223606797749979*alpha[6])*fUpwind_l[44]+0.2*(alpha[26]+alpha[25])*fUpwind_l[41]+0.223606797749979*(alpha[4]*fUpwind_l[41]+alpha[15]*(fUpwind_l[38]+fUpwind_l[37]))+(0.2*alpha[19]+0.223606797749979*alpha[2])*fUpwind_l[36]+(0.2*(fUpwind_l[29]+fUpwind_l[19])+0.223606797749979*fUpwind_l[2])*alpha[36]+(0.2*alpha[20]+0.223606797749979*alpha[1])*fUpwind_l[35]+0.2*(fUpwind_l[28]+fUpwind_l[20])*alpha[35]+0.223606797749979*(fUpwind_l[1]*alpha[35]+fUpwind_l[18]*alpha[33]+fUpwind_l[17]*alpha[32])+(0.223606797749979*(alpha[22]+alpha[21])+0.25*alpha[3])*fUpwind_l[31]+0.223606797749979*(alpha[8]*fUpwind_l[29]+alpha[9]*fUpwind_l[28]+alpha[5]*fUpwind_l[26]+fUpwind_l[5]*alpha[26]+alpha[5]*fUpwind_l[25]+fUpwind_l[5]*alpha[25]+alpha[9]*fUpwind_l[20]+fUpwind_l[9]*alpha[20]+alpha[8]*fUpwind_l[19]+fUpwind_l[8]*alpha[19])+0.25*(alpha[6]*fUpwind_l[18]+alpha[7]*fUpwind_l[17])+(0.223606797749979*(alpha[12]+alpha[11])+0.25*alpha[0])*fUpwind_l[16]+0.223606797749979*(fUpwind_l[14]+fUpwind_l[12]+fUpwind_l[11])*alpha[16]+0.25*(fUpwind_l[0]*alpha[16]+fUpwind_l[10]*alpha[15]+alpha[1]*fUpwind_l[9]+fUpwind_l[1]*alpha[9]+alpha[2]*fUpwind_l[8]+fUpwind_l[2]*alpha[8]+alpha[4]*fUpwind_l[5]+fUpwind_l[4]*alpha[5]); 
  Ghat_l[17] += (0.2*alpha[35]+0.223606797749979*alpha[9])*fUpwind_l[47]+(0.2*alpha[32]+0.223606797749979*alpha[7])*fUpwind_l[46]+0.2500000000000001*alpha[12]*fUpwind_l[45]+0.223606797749979*alpha[5]*fUpwind_l[44]+0.223606797749979*alpha[16]*fUpwind_l[43]+0.2*alpha[25]*fUpwind_l[42]+0.223606797749979*(alpha[4]*fUpwind_l[42]+alpha[15]*fUpwind_l[40])+(0.2*alpha[21]+0.223606797749979*alpha[3])*fUpwind_l[39]+0.2500000000000001*alpha[20]*fUpwind_l[38]+0.223606797749979*alpha[1]*fUpwind_l[37]+0.2500000000000001*(alpha[22]*fUpwind_l[36]+fUpwind_l[22]*alpha[36])+0.223606797749979*(alpha[15]*fUpwind_l[35]+fUpwind_l[15]*alpha[35])+0.2500000000000001*(alpha[26]*fUpwind_l[33]+fUpwind_l[26]*alpha[33])+0.223606797749979*(alpha[16]*fUpwind_l[32]+fUpwind_l[16]*alpha[32])+(0.223606797749979*alpha[19]+0.25*alpha[2])*fUpwind_l[31]+0.223606797749979*(alpha[8]*fUpwind_l[30]+alpha[6]*(fUpwind_l[27]+fUpwind_l[25])+fUpwind_l[6]*alpha[25]+alpha[8]*fUpwind_l[21]+fUpwind_l[8]*alpha[21])+0.25*alpha[5]*fUpwind_l[18]+0.223606797749979*alpha[11]*fUpwind_l[17]+0.25*(alpha[0]*fUpwind_l[17]+alpha[7]*fUpwind_l[16]+fUpwind_l[7]*alpha[16]+alpha[9]*fUpwind_l[15]+fUpwind_l[9]*alpha[15]+alpha[1]*fUpwind_l[10]+alpha[3]*fUpwind_l[8]+fUpwind_l[3]*alpha[8]+alpha[4]*fUpwind_l[6]+fUpwind_l[4]*alpha[6]); 
  Ghat_l[18] += (0.2*alpha[36]+0.223606797749979*alpha[8])*fUpwind_l[47]+0.2*alpha[33]*fUpwind_l[46]+0.223606797749979*(alpha[6]*fUpwind_l[46]+alpha[5]*fUpwind_l[45])+0.2500000000000001*alpha[11]*fUpwind_l[44]+0.2*alpha[26]*fUpwind_l[43]+0.223606797749979*(alpha[4]*fUpwind_l[43]+alpha[16]*fUpwind_l[42])+0.2*alpha[22]*fUpwind_l[40]+0.223606797749979*(alpha[3]*fUpwind_l[40]+alpha[15]*fUpwind_l[39]+alpha[2]*fUpwind_l[38])+0.2500000000000001*alpha[19]*fUpwind_l[37]+0.223606797749979*(alpha[15]*fUpwind_l[36]+fUpwind_l[15]*alpha[36])+0.2500000000000001*(alpha[21]*fUpwind_l[35]+fUpwind_l[21]*alpha[35])+0.223606797749979*(alpha[16]*fUpwind_l[33]+fUpwind_l[16]*alpha[33])+0.2500000000000001*(alpha[25]*fUpwind_l[32]+fUpwind_l[25]*alpha[32])+(0.223606797749979*alpha[20]+0.25*alpha[1])*fUpwind_l[31]+0.223606797749979*(alpha[9]*fUpwind_l[30]+alpha[7]*(fUpwind_l[27]+fUpwind_l[26])+fUpwind_l[7]*alpha[26]+alpha[9]*fUpwind_l[22]+fUpwind_l[9]*alpha[22])+0.223606797749979*alpha[12]*fUpwind_l[18]+0.25*(alpha[0]*fUpwind_l[18]+alpha[5]*fUpwind_l[17]+alpha[6]*fUpwind_l[16]+fUpwind_l[6]*alpha[16]+alpha[8]*fUpwind_l[15]+fUpwind_l[8]*alpha[15]+alpha[2]*fUpwind_l[10]+alpha[3]*fUpwind_l[9]+fUpwind_l[3]*alpha[9]+alpha[4]*fUpwind_l[7]+fUpwind_l[4]*alpha[7]); 
  Ghat_l[19] += 0.2*(alpha[16]*fUpwind_l[36]+fUpwind_l[16]*alpha[36])+(0.223606797749979*alpha[26]+0.159719141249985*alpha[25]+0.2500000000000001*alpha[4])*fUpwind_l[35]+(0.223606797749979*fUpwind_l[26]+0.159719141249985*fUpwind_l[25]+0.2500000000000001*fUpwind_l[4])*alpha[35]+0.2*(alpha[15]*fUpwind_l[33]+fUpwind_l[15]*alpha[33])+(0.223606797749979*alpha[22]+0.159719141249985*alpha[21]+0.2500000000000001*alpha[3])*fUpwind_l[32]+(0.223606797749979*fUpwind_l[22]+0.159719141249985*fUpwind_l[21]+0.2500000000000001*fUpwind_l[3])*alpha[32]+0.25*(alpha[9]*fUpwind_l[25]+fUpwind_l[9]*alpha[25]+alpha[7]*fUpwind_l[21]+fUpwind_l[7]*alpha[21])+0.2*(alpha[5]*fUpwind_l[20]+fUpwind_l[5]*alpha[20])+(0.223606797749979*alpha[12]+0.159719141249985*alpha[11]+0.25*alpha[0])*fUpwind_l[19]+(0.223606797749979*fUpwind_l[12]+0.159719141249985*fUpwind_l[11]+0.25*fUpwind_l[0])*alpha[19]+0.223606797749979*(alpha[8]*fUpwind_l[16]+fUpwind_l[8]*alpha[16]+alpha[6]*fUpwind_l[15]+fUpwind_l[6]*alpha[15])+0.2500000000000001*(alpha[2]*fUpwind_l[11]+fUpwind_l[2]*alpha[11])+0.223606797749979*(alpha[1]*fUpwind_l[5]+fUpwind_l[1]*alpha[5]); 
  Ghat_l[20] += (0.159719141249985*alpha[26]+0.223606797749979*alpha[25]+0.2500000000000001*alpha[4])*fUpwind_l[36]+(0.159719141249985*fUpwind_l[26]+0.223606797749979*fUpwind_l[25]+0.2500000000000001*fUpwind_l[4])*alpha[36]+0.2*(alpha[16]*fUpwind_l[35]+fUpwind_l[16]*alpha[35])+(0.159719141249985*alpha[22]+0.223606797749979*alpha[21]+0.2500000000000001*alpha[3])*fUpwind_l[33]+(0.159719141249985*fUpwind_l[22]+0.223606797749979*fUpwind_l[21]+0.2500000000000001*fUpwind_l[3])*alpha[33]+0.2*(alpha[15]*fUpwind_l[32]+fUpwind_l[15]*alpha[32])+0.25*(alpha[8]*fUpwind_l[26]+fUpwind_l[8]*alpha[26]+alpha[6]*fUpwind_l[22]+fUpwind_l[6]*alpha[22])+(0.159719141249985*alpha[12]+0.223606797749979*alpha[11]+0.25*alpha[0])*fUpwind_l[20]+(0.159719141249985*fUpwind_l[12]+0.223606797749979*fUpwind_l[11]+0.25*fUpwind_l[0])*alpha[20]+0.2*(alpha[5]*fUpwind_l[19]+fUpwind_l[5]*alpha[19])+0.223606797749979*(alpha[9]*fUpwind_l[16]+fUpwind_l[9]*alpha[16]+alpha[7]*fUpwind_l[15]+fUpwind_l[7]*alpha[15])+0.2500000000000001*(alpha[1]*fUpwind_l[12]+fUpwind_l[1]*alpha[12])+0.223606797749979*(alpha[2]*fUpwind_l[5]+fUpwind_l[2]*alpha[5]); 
  Ghat_l[21] += 0.223606797749979*alpha[36]*fUpwind_l[45]+(0.159719141249985*alpha[35]+0.25*alpha[9])*fUpwind_l[44]+0.159719141249985*alpha[25]*fUpwind_l[37]+0.2500000000000001*(alpha[4]*fUpwind_l[37]+fUpwind_l[18]*alpha[35])+0.2*alpha[15]*fUpwind_l[34]+0.223606797749979*(alpha[20]*fUpwind_l[33]+fUpwind_l[20]*alpha[33])+(0.159719141249985*alpha[19]+0.2500000000000001*alpha[2])*fUpwind_l[32]+(0.223606797749979*fUpwind_l[24]+0.159719141249985*fUpwind_l[19]+0.2500000000000001*fUpwind_l[2])*alpha[32]+0.223606797749979*alpha[16]*fUpwind_l[31]+0.25*fUpwind_l[10]*alpha[25]+0.2*alpha[6]*fUpwind_l[23]+(0.159719141249985*alpha[11]+0.25*alpha[0])*fUpwind_l[21]+(0.223606797749979*fUpwind_l[13]+0.159719141249985*fUpwind_l[11])*alpha[21]+0.25*(fUpwind_l[0]*alpha[21]+alpha[7]*fUpwind_l[19]+fUpwind_l[7]*alpha[19])+0.223606797749979*(alpha[8]*fUpwind_l[17]+alpha[5]*fUpwind_l[15]+fUpwind_l[5]*alpha[15])+0.2500000000000001*(alpha[3]*fUpwind_l[11]+fUpwind_l[3]*alpha[11])+0.223606797749979*(alpha[1]*fUpwind_l[6]+fUpwind_l[1]*alpha[6]); 
  Ghat_l[22] += (0.159719141249985*alpha[36]+0.25*alpha[8])*fUpwind_l[45]+0.223606797749979*alpha[35]*fUpwind_l[44]+0.159719141249985*alpha[26]*fUpwind_l[38]+0.2500000000000001*(alpha[4]*fUpwind_l[38]+fUpwind_l[17]*alpha[36])+0.2*alpha[15]*fUpwind_l[34]+(0.159719141249985*alpha[20]+0.2500000000000001*alpha[1])*fUpwind_l[33]+(0.223606797749979*fUpwind_l[23]+0.159719141249985*fUpwind_l[20]+0.2500000000000001*fUpwind_l[1])*alpha[33]+0.223606797749979*(alpha[19]*fUpwind_l[32]+fUpwind_l[19]*alpha[32])+0.223606797749979*alpha[16]*fUpwind_l[31]+0.25*fUpwind_l[10]*alpha[26]+0.2*alpha[7]*fUpwind_l[24]+(0.159719141249985*alpha[12]+0.25*alpha[0])*fUpwind_l[22]+(0.223606797749979*fUpwind_l[13]+0.159719141249985*fUpwind_l[12])*alpha[22]+0.25*(fUpwind_l[0]*alpha[22]+alpha[6]*fUpwind_l[20]+fUpwind_l[6]*alpha[20])+0.223606797749979*(alpha[9]*fUpwind_l[18]+alpha[5]*fUpwind_l[15]+fUpwind_l[5]*alpha[15])+0.2500000000000001*(alpha[3]*fUpwind_l[12]+fUpwind_l[3]*alpha[12])+0.223606797749979*(alpha[2]*fUpwind_l[7]+fUpwind_l[2]*alpha[7]); 
  Ghat_l[23] += (0.223606797749979*alpha[35]+0.25*alpha[9])*fUpwind_l[46]+0.2500000000000001*alpha[16]*fUpwind_l[40]+(0.223606797749979*alpha[25]+0.2500000000000001*alpha[4])*fUpwind_l[39]+(0.223606797749979*alpha[19]+0.2500000000000001*alpha[2])*fUpwind_l[34]+0.223606797749979*(alpha[22]*fUpwind_l[33]+fUpwind_l[22]*alpha[33])+0.2*(alpha[15]*fUpwind_l[32]+fUpwind_l[15]*alpha[32])+0.25*(alpha[8]*fUpwind_l[27]+alpha[5]*fUpwind_l[24])+(0.223606797749979*alpha[11]+0.25*alpha[0])*fUpwind_l[23]+0.2*(alpha[6]*fUpwind_l[21]+fUpwind_l[6]*alpha[21])+0.223606797749979*(alpha[7]*fUpwind_l[15]+fUpwind_l[7]*alpha[15])+0.2500000000000001*alpha[1]*fUpwind_l[13]+0.223606797749979*(alpha[3]*fUpwind_l[6]+fUpwind_l[3]*alpha[6]); 
  Ghat_l[24] += (0.223606797749979*alpha[36]+0.25*alpha[8])*fUpwind_l[46]+0.223606797749979*alpha[26]*fUpwind_l[40]+0.2500000000000001*(alpha[4]*fUpwind_l[40]+alpha[16]*fUpwind_l[39])+(0.223606797749979*alpha[20]+0.2500000000000001*alpha[1])*fUpwind_l[34]+0.2*(alpha[15]*fUpwind_l[33]+fUpwind_l[15]*alpha[33])+0.223606797749979*(alpha[21]*fUpwind_l[32]+fUpwind_l[21]*alpha[32])+0.25*alpha[9]*fUpwind_l[27]+0.223606797749979*alpha[12]*fUpwind_l[24]+0.25*(alpha[0]*fUpwind_l[24]+alpha[5]*fUpwind_l[23])+0.2*(alpha[7]*fUpwind_l[22]+fUpwind_l[7]*alpha[22])+0.223606797749979*(alpha[6]*fUpwind_l[15]+fUpwind_l[6]*alpha[15])+0.2500000000000001*alpha[2]*fUpwind_l[13]+0.223606797749979*(alpha[3]*fUpwind_l[7]+fUpwind_l[3]*alpha[7]); 
  Ghat_l[25] += 0.223606797749979*alpha[33]*fUpwind_l[45]+(0.159719141249985*alpha[32]+0.25*alpha[7])*fUpwind_l[44]+0.2*alpha[16]*fUpwind_l[41]+(0.159719141249985*alpha[21]+0.2500000000000001*alpha[3])*fUpwind_l[37]+0.223606797749979*(alpha[20]*fUpwind_l[36]+fUpwind_l[20]*alpha[36])+(0.159719141249985*alpha[19]+0.2500000000000001*alpha[2])*fUpwind_l[35]+(0.223606797749979*fUpwind_l[29]+0.159719141249985*fUpwind_l[19])*alpha[35]+0.2500000000000001*(fUpwind_l[2]*alpha[35]+fUpwind_l[18]*alpha[32])+0.223606797749979*alpha[15]*fUpwind_l[31]+0.2*alpha[8]*fUpwind_l[28]+(0.159719141249985*alpha[11]+0.25*alpha[0])*fUpwind_l[25]+(0.223606797749979*fUpwind_l[14]+0.159719141249985*fUpwind_l[11])*alpha[25]+0.25*(fUpwind_l[0]*alpha[25]+fUpwind_l[10]*alpha[21]+alpha[9]*fUpwind_l[19]+fUpwind_l[9]*alpha[19])+0.223606797749979*(alpha[6]*fUpwind_l[17]+alpha[5]*fUpwind_l[16]+fUpwind_l[5]*alpha[16])+0.2500000000000001*(alpha[4]*fUpwind_l[11]+fUpwind_l[4]*alpha[11])+0.223606797749979*(alpha[1]*fUpwind_l[8]+fUpwind_l[1]*alpha[8]); 
  Ghat_l[26] += (0.159719141249985*alpha[33]+0.25*alpha[6])*fUpwind_l[45]+0.223606797749979*alpha[32]*fUpwind_l[44]+0.2*alpha[16]*fUpwind_l[41]+(0.159719141249985*alpha[22]+0.2500000000000001*alpha[3])*fUpwind_l[38]+(0.159719141249985*alpha[20]+0.2500000000000001*alpha[1])*fUpwind_l[36]+(0.223606797749979*fUpwind_l[28]+0.159719141249985*fUpwind_l[20]+0.2500000000000001*fUpwind_l[1])*alpha[36]+0.223606797749979*(alpha[19]*fUpwind_l[35]+fUpwind_l[19]*alpha[35])+0.2500000000000001*fUpwind_l[17]*alpha[33]+0.223606797749979*alpha[15]*fUpwind_l[31]+0.2*alpha[9]*fUpwind_l[29]+(0.159719141249985*alpha[12]+0.25*alpha[0])*fUpwind_l[26]+(0.223606797749979*fUpwind_l[14]+0.159719141249985*fUpwind_l[12])*alpha[26]+0.25*(fUpwind_l[0]*alpha[26]+fUpwind_l[10]*alpha[22]+alpha[8]*fUpwind_l[20]+fUpwind_l[8]*alpha[20])+0.223606797749979*(alpha[7]*fUpwind_l[18]+alpha[5]*fUpwind_l[16]+fUpwind_l[5]*alpha[16])+0.2500000000000001*(alpha[4]*fUpwind_l[12]+fUpwind_l[4]*alpha[12])+0.223606797749979*(alpha[2]*fUpwind_l[9]+fUpwind_l[2]*alpha[9]); 
  Ghat_l[27] += 0.25*alpha[5]*fUpwind_l[46]+0.223606797749979*(alpha[33]*fUpwind_l[45]+alpha[32]*fUpwind_l[44])+0.2500000000000001*(alpha[2]*fUpwind_l[40]+alpha[1]*fUpwind_l[39])+0.223606797749979*(alpha[22]*fUpwind_l[38]+alpha[21]*fUpwind_l[37])+0.2500000000000001*alpha[16]*fUpwind_l[34]+0.223606797749979*alpha[15]*fUpwind_l[31]+0.25*(alpha[0]*fUpwind_l[27]+alpha[9]*fUpwind_l[24]+alpha[8]*fUpwind_l[23])+0.223606797749979*(alpha[7]*fUpwind_l[18]+alpha[6]*fUpwind_l[17])+0.2500000000000001*alpha[4]*fUpwind_l[13]+0.223606797749979*alpha[3]*fUpwind_l[10]; 
  Ghat_l[28] += (0.223606797749979*alpha[32]+0.25*alpha[7])*fUpwind_l[47]+0.2500000000000001*alpha[15]*fUpwind_l[43]+(0.223606797749979*alpha[21]+0.2500000000000001*alpha[3])*fUpwind_l[42]+(0.223606797749979*alpha[19]+0.2500000000000001*alpha[2])*fUpwind_l[41]+0.223606797749979*(alpha[26]*fUpwind_l[36]+fUpwind_l[26]*alpha[36])+0.2*(alpha[16]*fUpwind_l[35]+fUpwind_l[16]*alpha[35])+0.25*(alpha[6]*fUpwind_l[30]+alpha[5]*fUpwind_l[29])+(0.223606797749979*alpha[11]+0.25*alpha[0])*fUpwind_l[28]+0.2*(alpha[8]*fUpwind_l[25]+fUpwind_l[8]*alpha[25])+0.223606797749979*(alpha[9]*fUpwind_l[16]+fUpwind_l[9]*alpha[16])+0.2500000000000001*alpha[1]*fUpwind_l[14]+0.223606797749979*(alpha[4]*fUpwind_l[8]+fUpwind_l[4]*alpha[8]); 
  Ghat_l[29] += (0.223606797749979*alpha[33]+0.25*alpha[6])*fUpwind_l[47]+0.223606797749979*alpha[22]*fUpwind_l[43]+0.2500000000000001*(alpha[3]*fUpwind_l[43]+alpha[15]*fUpwind_l[42])+(0.223606797749979*alpha[20]+0.2500000000000001*alpha[1])*fUpwind_l[41]+0.2*(alpha[16]*fUpwind_l[36]+fUpwind_l[16]*alpha[36])+0.223606797749979*(alpha[25]*fUpwind_l[35]+fUpwind_l[25]*alpha[35])+0.25*alpha[7]*fUpwind_l[30]+0.223606797749979*alpha[12]*fUpwind_l[29]+0.25*(alpha[0]*fUpwind_l[29]+alpha[5]*fUpwind_l[28])+0.2*(alpha[9]*fUpwind_l[26]+fUpwind_l[9]*alpha[26])+0.223606797749979*(alpha[8]*fUpwind_l[16]+fUpwind_l[8]*alpha[16])+0.2500000000000001*alpha[2]*fUpwind_l[14]+0.223606797749979*(alpha[4]*fUpwind_l[9]+fUpwind_l[4]*alpha[9]); 
  Ghat_l[30] += 0.25*alpha[5]*fUpwind_l[47]+0.223606797749979*(alpha[36]*fUpwind_l[45]+alpha[35]*fUpwind_l[44])+0.2500000000000001*(alpha[2]*fUpwind_l[43]+alpha[1]*fUpwind_l[42]+alpha[15]*fUpwind_l[41])+0.223606797749979*(alpha[26]*fUpwind_l[38]+alpha[25]*fUpwind_l[37])+0.223606797749979*alpha[16]*fUpwind_l[31]+0.25*(alpha[0]*fUpwind_l[30]+alpha[7]*fUpwind_l[29]+alpha[6]*fUpwind_l[28])+0.223606797749979*(alpha[9]*fUpwind_l[18]+alpha[8]*fUpwind_l[17])+0.2500000000000001*alpha[3]*fUpwind_l[14]+0.223606797749979*alpha[4]*fUpwind_l[10]; 
  Ghat_l[31] += (0.2*(alpha[26]+alpha[25])+0.223606797749979*alpha[4])*fUpwind_l[47]+(0.2*(alpha[22]+alpha[21])+0.223606797749979*alpha[3])*fUpwind_l[46]+(0.2*alpha[19]+0.223606797749979*alpha[2])*fUpwind_l[45]+(0.2*alpha[20]+0.223606797749979*alpha[1])*fUpwind_l[44]+(0.2*alpha[36]+0.223606797749979*alpha[8])*fUpwind_l[43]+(0.2*alpha[35]+0.223606797749979*alpha[9])*fUpwind_l[42]+(0.2*alpha[33]+0.223606797749979*alpha[6])*fUpwind_l[40]+0.2*alpha[32]*fUpwind_l[39]+0.223606797749979*(alpha[7]*fUpwind_l[39]+alpha[5]*(fUpwind_l[38]+fUpwind_l[37]))+(0.2*alpha[32]+0.223606797749979*alpha[7])*fUpwind_l[36]+(0.2*fUpwind_l[32]+0.223606797749979*fUpwind_l[7])*alpha[36]+(0.2*alpha[33]+0.223606797749979*alpha[6])*fUpwind_l[35]+0.2*fUpwind_l[33]*alpha[35]+0.223606797749979*(fUpwind_l[6]*alpha[35]+alpha[9]*fUpwind_l[33]+fUpwind_l[9]*alpha[33]+alpha[8]*fUpwind_l[32]+fUpwind_l[8]*alpha[32])+(0.223606797749979*(alpha[12]+alpha[11])+0.25*alpha[0])*fUpwind_l[31]+0.223606797749979*(alpha[16]*fUpwind_l[30]+alpha[15]*(fUpwind_l[27]+fUpwind_l[26])+fUpwind_l[15]*alpha[26]+alpha[15]*fUpwind_l[25]+fUpwind_l[15]*alpha[25]+alpha[16]*fUpwind_l[22]+fUpwind_l[16]*alpha[22]+alpha[16]*fUpwind_l[21]+fUpwind_l[16]*alpha[21]+fUpwind_l[18]*alpha[20]+fUpwind_l[17]*alpha[19])+0.25*(alpha[1]*fUpwind_l[18]+alpha[2]*fUpwind_l[17]+alpha[3]*fUpwind_l[16]+fUpwind_l[3]*alpha[16]+alpha[4]*fUpwind_l[15]+fUpwind_l[4]*alpha[15]+alpha[5]*fUpwind_l[10]+alpha[6]*fUpwind_l[9]+fUpwind_l[6]*alpha[9]+alpha[7]*fUpwind_l[8]+fUpwind_l[7]*alpha[8]); 
  Ghat_l[32] += 0.2*alpha[16]*fUpwind_l[45]+(0.223606797749979*alpha[26]+0.159719141249985*alpha[25]+0.2500000000000001*alpha[4])*fUpwind_l[44]+0.223606797749979*alpha[35]*fUpwind_l[38]+(0.159719141249985*alpha[35]+0.25*alpha[9])*fUpwind_l[37]+0.2*fUpwind_l[31]*alpha[36]+0.25*fUpwind_l[10]*alpha[35]+0.1788854381999831*alpha[33]*fUpwind_l[34]+0.2*(alpha[6]*fUpwind_l[34]+alpha[5]*fUpwind_l[33]+fUpwind_l[5]*alpha[33])+(0.223606797749979*alpha[12]+0.159719141249985*alpha[11]+0.25*alpha[0])*fUpwind_l[32]+(0.223606797749979*(fUpwind_l[13]+fUpwind_l[12])+0.159719141249985*fUpwind_l[11]+0.25*fUpwind_l[0])*alpha[32]+0.223606797749979*alpha[8]*fUpwind_l[31]+0.2500000000000001*fUpwind_l[18]*alpha[25]+0.223606797749979*alpha[21]*fUpwind_l[24]+0.2*alpha[15]*fUpwind_l[23]+0.223606797749979*(alpha[19]*fUpwind_l[22]+fUpwind_l[19]*alpha[22])+(0.159719141249985*alpha[19]+0.2500000000000001*alpha[2])*fUpwind_l[21]+(0.159719141249985*fUpwind_l[19]+0.2500000000000001*fUpwind_l[2])*alpha[21]+0.2*(alpha[15]*fUpwind_l[20]+fUpwind_l[15]*alpha[20])+0.2500000000000001*(alpha[3]*fUpwind_l[19]+fUpwind_l[3]*alpha[19])+0.223606797749979*(alpha[16]*fUpwind_l[17]+alpha[1]*fUpwind_l[15]+fUpwind_l[1]*alpha[15])+0.25*(alpha[7]*fUpwind_l[11]+fUpwind_l[7]*alpha[11])+0.223606797749979*(alpha[5]*fUpwind_l[6]+fUpwind_l[5]*alpha[6]); 
  Ghat_l[33] += (0.159719141249985*alpha[26]+0.223606797749979*alpha[25]+0.2500000000000001*alpha[4])*fUpwind_l[45]+0.2*alpha[16]*fUpwind_l[44]+(0.159719141249985*alpha[36]+0.25*alpha[8])*fUpwind_l[38]+alpha[36]*(0.223606797749979*fUpwind_l[37]+0.25*fUpwind_l[10])+0.2*fUpwind_l[31]*alpha[35]+(0.1788854381999831*alpha[32]+0.2*alpha[7])*fUpwind_l[34]+(0.159719141249985*alpha[12]+0.223606797749979*alpha[11]+0.25*alpha[0])*fUpwind_l[33]+(0.223606797749979*fUpwind_l[13]+0.159719141249985*fUpwind_l[12]+0.223606797749979*fUpwind_l[11]+0.25*fUpwind_l[0])*alpha[33]+0.2*(alpha[5]*fUpwind_l[32]+fUpwind_l[5]*alpha[32])+0.223606797749979*alpha[9]*fUpwind_l[31]+0.2500000000000001*fUpwind_l[17]*alpha[26]+0.2*alpha[15]*fUpwind_l[24]+0.223606797749979*alpha[22]*fUpwind_l[23]+(0.159719141249985*alpha[20]+0.2500000000000001*alpha[1])*fUpwind_l[22]+(0.159719141249985*fUpwind_l[20]+0.2500000000000001*fUpwind_l[1])*alpha[22]+0.223606797749979*(alpha[20]*fUpwind_l[21]+fUpwind_l[20]*alpha[21])+0.2500000000000001*(alpha[3]*fUpwind_l[20]+fUpwind_l[3]*alpha[20])+0.2*(alpha[15]*fUpwind_l[19]+fUpwind_l[15]*alpha[19])+0.223606797749979*(alpha[16]*fUpwind_l[18]+alpha[2]*fUpwind_l[15]+fUpwind_l[2]*alpha[15])+0.25*(alpha[6]*fUpwind_l[12]+fUpwind_l[6]*alpha[12])+0.223606797749979*(alpha[5]*fUpwind_l[7]+fUpwind_l[5]*alpha[7]); 
  Ghat_l[34] += (0.223606797749979*(alpha[26]+alpha[25])+0.2500000000000001*alpha[4])*fUpwind_l[46]+(0.223606797749979*alpha[36]+0.25*alpha[8])*fUpwind_l[40]+(0.223606797749979*alpha[35]+0.25*alpha[9])*fUpwind_l[39]+(0.223606797749979*(alpha[12]+alpha[11])+0.25*alpha[0])*fUpwind_l[34]+(0.1788854381999831*alpha[32]+0.2*alpha[7])*fUpwind_l[33]+0.1788854381999831*fUpwind_l[32]*alpha[33]+0.2*(fUpwind_l[7]*alpha[33]+alpha[6]*fUpwind_l[32]+fUpwind_l[6]*alpha[32])+0.2500000000000001*alpha[16]*fUpwind_l[27]+(0.223606797749979*alpha[20]+0.2500000000000001*alpha[1])*fUpwind_l[24]+(0.223606797749979*alpha[19]+0.2500000000000001*alpha[2])*fUpwind_l[23]+0.2*(alpha[15]*fUpwind_l[22]+fUpwind_l[15]*alpha[22]+alpha[15]*fUpwind_l[21]+fUpwind_l[15]*alpha[21])+0.223606797749979*(alpha[3]*fUpwind_l[15]+fUpwind_l[3]*alpha[15])+0.25*alpha[5]*fUpwind_l[13]+0.223606797749979*(alpha[6]*fUpwind_l[7]+fUpwind_l[6]*alpha[7]); 
  Ghat_l[35] += 0.2*alpha[15]*fUpwind_l[45]+(0.223606797749979*alpha[22]+0.159719141249985*alpha[21]+0.2500000000000001*alpha[3])*fUpwind_l[44]+(0.1788854381999831*alpha[36]+0.2*alpha[8])*fUpwind_l[41]+0.223606797749979*alpha[32]*fUpwind_l[38]+(0.159719141249985*alpha[32]+0.25*alpha[7])*fUpwind_l[37]+0.2*(alpha[5]*fUpwind_l[36]+fUpwind_l[5]*alpha[36])+(0.223606797749979*alpha[12]+0.159719141249985*alpha[11]+0.25*alpha[0])*fUpwind_l[35]+(0.223606797749979*(fUpwind_l[14]+fUpwind_l[12])+0.159719141249985*fUpwind_l[11]+0.25*fUpwind_l[0])*alpha[35]+0.2*fUpwind_l[31]*alpha[33]+0.25*fUpwind_l[10]*alpha[32]+0.223606797749979*(alpha[6]*fUpwind_l[31]+alpha[25]*fUpwind_l[29])+0.2*alpha[16]*fUpwind_l[28]+0.223606797749979*(alpha[19]*fUpwind_l[26]+fUpwind_l[19]*alpha[26])+(0.159719141249985*alpha[19]+0.2500000000000001*alpha[2])*fUpwind_l[25]+0.159719141249985*fUpwind_l[19]*alpha[25]+0.2500000000000001*(fUpwind_l[2]*alpha[25]+fUpwind_l[18]*alpha[21])+0.2*(alpha[16]*fUpwind_l[20]+fUpwind_l[16]*alpha[20])+0.2500000000000001*(alpha[4]*fUpwind_l[19]+fUpwind_l[4]*alpha[19])+0.223606797749979*(alpha[15]*fUpwind_l[17]+alpha[1]*fUpwind_l[16]+fUpwind_l[1]*alpha[16])+0.25*(alpha[9]*fUpwind_l[11]+fUpwind_l[9]*alpha[11])+0.223606797749979*(alpha[5]*fUpwind_l[8]+fUpwind_l[5]*alpha[8]); 
  Ghat_l[36] += (0.159719141249985*alpha[22]+0.223606797749979*alpha[21]+0.2500000000000001*alpha[3])*fUpwind_l[45]+0.2*alpha[15]*fUpwind_l[44]+(0.1788854381999831*alpha[35]+0.2*alpha[9])*fUpwind_l[41]+(0.159719141249985*alpha[33]+0.25*alpha[6])*fUpwind_l[38]+0.223606797749979*alpha[33]*fUpwind_l[37]+(0.159719141249985*alpha[12]+0.223606797749979*alpha[11]+0.25*alpha[0])*fUpwind_l[36]+(0.223606797749979*fUpwind_l[14]+0.159719141249985*fUpwind_l[12]+0.223606797749979*fUpwind_l[11]+0.25*fUpwind_l[0])*alpha[36]+0.2*(alpha[5]*fUpwind_l[35]+fUpwind_l[5]*alpha[35])+0.25*fUpwind_l[10]*alpha[33]+fUpwind_l[31]*(0.2*alpha[32]+0.223606797749979*alpha[7])+0.2*alpha[16]*fUpwind_l[29]+0.223606797749979*alpha[26]*fUpwind_l[28]+(0.159719141249985*alpha[20]+0.2500000000000001*alpha[1])*fUpwind_l[26]+(0.159719141249985*fUpwind_l[20]+0.2500000000000001*fUpwind_l[1])*alpha[26]+0.223606797749979*(alpha[20]*fUpwind_l[25]+fUpwind_l[20]*alpha[25])+0.2500000000000001*(fUpwind_l[17]*alpha[22]+alpha[4]*fUpwind_l[20]+fUpwind_l[4]*alpha[20])+0.2*(alpha[16]*fUpwind_l[19]+fUpwind_l[16]*alpha[19])+0.223606797749979*(alpha[15]*fUpwind_l[18]+alpha[2]*fUpwind_l[16]+fUpwind_l[2]*alpha[16])+0.25*(alpha[8]*fUpwind_l[12]+fUpwind_l[8]*alpha[12])+0.223606797749979*(alpha[5]*fUpwind_l[9]+fUpwind_l[5]*alpha[9]); 
  Ghat_l[37] += 0.2*(alpha[16]*fUpwind_l[47]+alpha[15]*fUpwind_l[46])+0.223606797749979*alpha[20]*fUpwind_l[45]+(0.159719141249985*alpha[19]+0.2500000000000001*alpha[2])*fUpwind_l[44]+0.223606797749979*alpha[35]*fUpwind_l[43]+0.2*alpha[8]*fUpwind_l[42]+0.223606797749979*alpha[32]*fUpwind_l[40]+0.2*alpha[6]*fUpwind_l[39]+(0.159719141249985*alpha[11]+0.25*alpha[0])*fUpwind_l[37]+0.223606797749979*(alpha[33]*fUpwind_l[36]+fUpwind_l[33]*alpha[36])+(0.159719141249985*alpha[32]+0.25*alpha[7])*fUpwind_l[35]+0.159719141249985*fUpwind_l[32]*alpha[35]+0.25*(fUpwind_l[7]*alpha[35]+alpha[9]*fUpwind_l[32]+fUpwind_l[9]*alpha[32])+0.223606797749979*(alpha[5]*fUpwind_l[31]+alpha[25]*fUpwind_l[30]+alpha[21]*fUpwind_l[27])+(0.159719141249985*alpha[21]+0.2500000000000001*alpha[3])*fUpwind_l[25]+0.159719141249985*fUpwind_l[21]*alpha[25]+0.2500000000000001*(fUpwind_l[3]*alpha[25]+alpha[4]*fUpwind_l[21]+fUpwind_l[4]*alpha[21]+fUpwind_l[18]*alpha[19])+0.223606797749979*(alpha[1]*fUpwind_l[17]+alpha[15]*fUpwind_l[16]+fUpwind_l[15]*alpha[16])+0.25*fUpwind_l[10]*alpha[11]+0.223606797749979*(alpha[6]*fUpwind_l[8]+fUpwind_l[6]*alpha[8]); 
  Ghat_l[38] += 0.2*(alpha[16]*fUpwind_l[47]+alpha[15]*fUpwind_l[46])+(0.159719141249985*alpha[20]+0.2500000000000001*alpha[1])*fUpwind_l[45]+0.223606797749979*alpha[19]*fUpwind_l[44]+0.2*alpha[9]*fUpwind_l[43]+0.223606797749979*alpha[36]*fUpwind_l[42]+0.2*alpha[7]*fUpwind_l[40]+0.223606797749979*alpha[33]*fUpwind_l[39]+(0.159719141249985*alpha[12]+0.25*alpha[0])*fUpwind_l[38]+(0.159719141249985*alpha[33]+0.25*alpha[6])*fUpwind_l[36]+(0.159719141249985*fUpwind_l[33]+0.25*fUpwind_l[6])*alpha[36]+0.223606797749979*(alpha[32]*fUpwind_l[35]+fUpwind_l[32]*alpha[35])+0.25*(alpha[8]*fUpwind_l[33]+fUpwind_l[8]*alpha[33])+0.223606797749979*(alpha[5]*fUpwind_l[31]+alpha[26]*fUpwind_l[30]+alpha[22]*fUpwind_l[27])+(0.159719141249985*alpha[22]+0.2500000000000001*alpha[3])*fUpwind_l[26]+0.159719141249985*fUpwind_l[22]*alpha[26]+0.2500000000000001*(fUpwind_l[3]*alpha[26]+alpha[4]*fUpwind_l[22]+fUpwind_l[4]*alpha[22]+fUpwind_l[17]*alpha[20])+0.223606797749979*(alpha[2]*fUpwind_l[18]+alpha[15]*fUpwind_l[16]+fUpwind_l[15]*alpha[16])+0.25*fUpwind_l[10]*alpha[12]+0.223606797749979*(alpha[7]*fUpwind_l[9]+fUpwind_l[7]*alpha[9]); 
  Ghat_l[39] += (0.223606797749979*alpha[19]+0.2500000000000001*alpha[2])*fUpwind_l[46]+0.223606797749979*alpha[22]*fUpwind_l[45]+0.2*alpha[15]*fUpwind_l[44]+0.25*alpha[5]*fUpwind_l[40]+(0.223606797749979*alpha[11]+0.25*alpha[0])*fUpwind_l[39]+0.223606797749979*alpha[33]*fUpwind_l[38]+0.2*alpha[6]*fUpwind_l[37]+fUpwind_l[34]*(0.223606797749979*alpha[35]+0.25*alpha[9])+fUpwind_l[31]*(0.2*alpha[32]+0.223606797749979*alpha[7])+0.2500000000000001*alpha[1]*fUpwind_l[27]+0.223606797749979*fUpwind_l[23]*alpha[25]+0.2500000000000001*(alpha[16]*fUpwind_l[24]+alpha[4]*fUpwind_l[23])+0.2*fUpwind_l[17]*alpha[21]+0.223606797749979*(alpha[15]*fUpwind_l[18]+alpha[3]*fUpwind_l[17])+0.25*alpha[8]*fUpwind_l[13]+0.223606797749979*alpha[6]*fUpwind_l[10]; 
  Ghat_l[40] += (0.223606797749979*alpha[20]+0.2500000000000001*alpha[1])*fUpwind_l[46]+0.2*alpha[15]*fUpwind_l[45]+0.223606797749979*(alpha[21]*fUpwind_l[44]+alpha[12]*fUpwind_l[40])+0.25*(alpha[0]*fUpwind_l[40]+alpha[5]*fUpwind_l[39])+0.2*alpha[7]*fUpwind_l[38]+0.223606797749979*alpha[32]*fUpwind_l[37]+fUpwind_l[34]*(0.223606797749979*alpha[36]+0.25*alpha[8])+fUpwind_l[31]*(0.2*alpha[33]+0.223606797749979*alpha[6])+0.2500000000000001*alpha[2]*fUpwind_l[27]+0.223606797749979*fUpwind_l[24]*alpha[26]+0.2500000000000001*(alpha[4]*fUpwind_l[24]+alpha[16]*fUpwind_l[23])+0.2*fUpwind_l[18]*alpha[22]+0.223606797749979*(alpha[3]*fUpwind_l[18]+alpha[15]*fUpwind_l[17])+0.25*alpha[9]*fUpwind_l[13]+0.223606797749979*alpha[7]*fUpwind_l[10]; 
  Ghat_l[41] += (0.223606797749979*(alpha[22]+alpha[21])+0.2500000000000001*alpha[3])*fUpwind_l[47]+(0.223606797749979*alpha[33]+0.25*alpha[6])*fUpwind_l[43]+(0.223606797749979*alpha[32]+0.25*alpha[7])*fUpwind_l[42]+(0.223606797749979*(alpha[12]+alpha[11])+0.25*alpha[0])*fUpwind_l[41]+(0.1788854381999831*alpha[35]+0.2*alpha[9])*fUpwind_l[36]+0.1788854381999831*fUpwind_l[35]*alpha[36]+0.2*(fUpwind_l[9]*alpha[36]+alpha[8]*fUpwind_l[35]+fUpwind_l[8]*alpha[35])+0.2500000000000001*alpha[15]*fUpwind_l[30]+(0.223606797749979*alpha[20]+0.2500000000000001*alpha[1])*fUpwind_l[29]+(0.223606797749979*alpha[19]+0.2500000000000001*alpha[2])*fUpwind_l[28]+0.2*(alpha[16]*fUpwind_l[26]+fUpwind_l[16]*alpha[26]+alpha[16]*fUpwind_l[25]+fUpwind_l[16]*alpha[25])+0.223606797749979*(alpha[4]*fUpwind_l[16]+fUpwind_l[4]*alpha[16])+0.25*alpha[5]*fUpwind_l[14]+0.223606797749979*(alpha[8]*fUpwind_l[9]+fUpwind_l[8]*alpha[9]); 
  Ghat_l[42] += (0.223606797749979*alpha[19]+0.2500000000000001*alpha[2])*fUpwind_l[47]+0.223606797749979*alpha[26]*fUpwind_l[45]+0.2*alpha[16]*fUpwind_l[44]+0.25*alpha[5]*fUpwind_l[43]+(0.223606797749979*alpha[11]+0.25*alpha[0])*fUpwind_l[42]+(0.223606797749979*alpha[32]+0.25*alpha[7])*fUpwind_l[41]+0.223606797749979*alpha[36]*fUpwind_l[38]+0.2*alpha[8]*fUpwind_l[37]+fUpwind_l[31]*(0.2*alpha[35]+0.223606797749979*alpha[9])+0.2500000000000001*(alpha[1]*fUpwind_l[30]+alpha[15]*fUpwind_l[29])+(0.223606797749979*alpha[21]+0.2500000000000001*alpha[3])*fUpwind_l[28]+0.2*fUpwind_l[17]*alpha[25]+0.223606797749979*(alpha[16]*fUpwind_l[18]+alpha[4]*fUpwind_l[17])+0.25*alpha[6]*fUpwind_l[14]+0.223606797749979*alpha[8]*fUpwind_l[10]; 
  Ghat_l[43] += (0.223606797749979*alpha[20]+0.2500000000000001*alpha[1])*fUpwind_l[47]+0.2*alpha[16]*fUpwind_l[45]+0.223606797749979*(alpha[25]*fUpwind_l[44]+alpha[12]*fUpwind_l[43])+0.25*(alpha[0]*fUpwind_l[43]+alpha[5]*fUpwind_l[42])+(0.223606797749979*alpha[33]+0.25*alpha[6])*fUpwind_l[41]+0.2*alpha[9]*fUpwind_l[38]+0.223606797749979*alpha[35]*fUpwind_l[37]+fUpwind_l[31]*(0.2*alpha[36]+0.223606797749979*alpha[8])+0.2500000000000001*alpha[2]*fUpwind_l[30]+0.223606797749979*alpha[22]*fUpwind_l[29]+0.2500000000000001*(alpha[3]*fUpwind_l[29]+alpha[15]*fUpwind_l[28])+0.2*fUpwind_l[18]*alpha[26]+0.223606797749979*(alpha[4]*fUpwind_l[18]+alpha[16]*fUpwind_l[17])+0.25*alpha[7]*fUpwind_l[14]+0.223606797749979*alpha[9]*fUpwind_l[10]; 
  Ghat_l[44] += (0.1788854381999831*alpha[36]+0.2*alpha[8])*fUpwind_l[47]+0.1788854381999831*alpha[33]*fUpwind_l[46]+0.2*(alpha[6]*fUpwind_l[46]+alpha[5]*fUpwind_l[45])+(0.223606797749979*alpha[12]+0.159719141249985*alpha[11]+0.25*alpha[0])*fUpwind_l[44]+0.223606797749979*alpha[25]*fUpwind_l[43]+0.2*alpha[16]*fUpwind_l[42]+0.223606797749979*alpha[21]*fUpwind_l[40]+0.2*alpha[15]*fUpwind_l[39]+0.223606797749979*alpha[19]*fUpwind_l[38]+(0.159719141249985*alpha[19]+0.2500000000000001*alpha[2])*fUpwind_l[37]+0.2*(alpha[15]*fUpwind_l[36]+fUpwind_l[15]*alpha[36])+(0.223606797749979*alpha[22]+0.159719141249985*alpha[21]+0.2500000000000001*alpha[3])*fUpwind_l[35]+(0.223606797749979*(fUpwind_l[30]+fUpwind_l[22])+0.159719141249985*fUpwind_l[21]+0.2500000000000001*fUpwind_l[3])*alpha[35]+0.2*(alpha[16]*fUpwind_l[33]+fUpwind_l[16]*alpha[33])+(0.223606797749979*alpha[26]+0.159719141249985*alpha[25]+0.2500000000000001*alpha[4])*fUpwind_l[32]+(0.223606797749979*(fUpwind_l[27]+fUpwind_l[26])+0.159719141249985*fUpwind_l[25]+0.2500000000000001*fUpwind_l[4])*alpha[32]+(0.2*alpha[20]+0.223606797749979*alpha[1])*fUpwind_l[31]+0.25*(alpha[7]*fUpwind_l[25]+fUpwind_l[7]*alpha[25]+alpha[9]*fUpwind_l[21]+fUpwind_l[9]*alpha[21]+fUpwind_l[10]*alpha[19])+0.2500000000000001*alpha[11]*fUpwind_l[18]+0.223606797749979*(alpha[5]*fUpwind_l[17]+alpha[6]*fUpwind_l[16]+fUpwind_l[6]*alpha[16]+alpha[8]*fUpwind_l[15]+fUpwind_l[8]*alpha[15]); 
  Ghat_l[45] += (0.1788854381999831*alpha[35]+0.2*alpha[9])*fUpwind_l[47]+(0.1788854381999831*alpha[32]+0.2*alpha[7])*fUpwind_l[46]+(0.159719141249985*alpha[12]+0.223606797749979*alpha[11]+0.25*alpha[0])*fUpwind_l[45]+0.2*(alpha[5]*fUpwind_l[44]+alpha[16]*fUpwind_l[43])+0.223606797749979*alpha[26]*fUpwind_l[42]+0.2*alpha[15]*fUpwind_l[40]+0.223606797749979*alpha[22]*fUpwind_l[39]+(0.159719141249985*alpha[20]+0.2500000000000001*alpha[1])*fUpwind_l[38]+0.223606797749979*alpha[20]*fUpwind_l[37]+(0.159719141249985*alpha[22]+0.223606797749979*alpha[21]+0.2500000000000001*alpha[3])*fUpwind_l[36]+(0.223606797749979*fUpwind_l[30]+0.159719141249985*fUpwind_l[22]+0.223606797749979*fUpwind_l[21]+0.2500000000000001*fUpwind_l[3])*alpha[36]+0.2*(alpha[15]*fUpwind_l[35]+fUpwind_l[15]*alpha[35])+(0.159719141249985*alpha[26]+0.223606797749979*alpha[25]+0.2500000000000001*alpha[4])*fUpwind_l[33]+(0.223606797749979*fUpwind_l[27]+0.159719141249985*fUpwind_l[26]+0.223606797749979*fUpwind_l[25]+0.2500000000000001*fUpwind_l[4])*alpha[33]+0.2*(alpha[16]*fUpwind_l[32]+fUpwind_l[16]*alpha[32])+(0.2*alpha[19]+0.223606797749979*alpha[2])*fUpwind_l[31]+0.25*(alpha[6]*fUpwind_l[26]+fUpwind_l[6]*alpha[26]+alpha[8]*fUpwind_l[22]+fUpwind_l[8]*alpha[22]+fUpwind_l[10]*alpha[20])+0.223606797749979*alpha[5]*fUpwind_l[18]+0.2500000000000001*alpha[12]*fUpwind_l[17]+0.223606797749979*(alpha[7]*fUpwind_l[16]+fUpwind_l[7]*alpha[16]+alpha[9]*fUpwind_l[15]+fUpwind_l[9]*alpha[15]); 
  Ghat_l[46] += (0.223606797749979*(alpha[12]+alpha[11])+0.25*alpha[0])*fUpwind_l[46]+(0.1788854381999831*alpha[32]+0.2*alpha[7])*fUpwind_l[45]+(0.1788854381999831*alpha[33]+0.2*alpha[6])*fUpwind_l[44]+(0.223606797749979*alpha[20]+0.2500000000000001*alpha[1])*fUpwind_l[40]+(0.223606797749979*alpha[19]+0.2500000000000001*alpha[2])*fUpwind_l[39]+0.2*alpha[15]*(fUpwind_l[38]+fUpwind_l[37])+0.223606797749979*(fUpwind_l[24]*alpha[36]+fUpwind_l[23]*alpha[35])+(0.223606797749979*(alpha[26]+alpha[25])+0.2500000000000001*alpha[4])*fUpwind_l[34]+0.2*(fUpwind_l[18]*alpha[33]+fUpwind_l[17]*alpha[32])+(0.2*(alpha[22]+alpha[21])+0.223606797749979*alpha[3])*fUpwind_l[31]+0.25*(alpha[5]*fUpwind_l[27]+alpha[8]*fUpwind_l[24]+alpha[9]*fUpwind_l[23])+0.223606797749979*(alpha[6]*fUpwind_l[18]+alpha[7]*fUpwind_l[17])+0.2500000000000001*fUpwind_l[13]*alpha[16]+0.223606797749979*fUpwind_l[10]*alpha[15]; 
  Ghat_l[47] += (0.223606797749979*(alpha[12]+alpha[11])+0.25*alpha[0])*fUpwind_l[47]+(0.1788854381999831*alpha[35]+0.2*alpha[9])*fUpwind_l[45]+(0.1788854381999831*alpha[36]+0.2*alpha[8])*fUpwind_l[44]+(0.223606797749979*alpha[20]+0.2500000000000001*alpha[1])*fUpwind_l[43]+(0.223606797749979*alpha[19]+0.2500000000000001*alpha[2])*fUpwind_l[42]+(0.223606797749979*(alpha[22]+alpha[21])+0.2500000000000001*alpha[3])*fUpwind_l[41]+0.2*(alpha[16]*(fUpwind_l[38]+fUpwind_l[37])+fUpwind_l[18]*alpha[36]+fUpwind_l[17]*alpha[35])+0.223606797749979*(fUpwind_l[29]*alpha[33]+fUpwind_l[28]*alpha[32])+(0.2*(alpha[26]+alpha[25])+0.223606797749979*alpha[4])*fUpwind_l[31]+0.25*(alpha[5]*fUpwind_l[30]+alpha[6]*fUpwind_l[29]+alpha[7]*fUpwind_l[28])+0.223606797749979*(alpha[8]*fUpwind_l[18]+alpha[9]*fUpwind_l[17]+fUpwind_l[10]*alpha[16])+0.2500000000000001*fUpwind_l[14]*alpha[15]; 

  Ghat_r[0] += 0.25*(alpha[36]*fUpwind_r[36]+alpha[35]*fUpwind_r[35]+alpha[33]*fUpwind_r[33]+alpha[32]*fUpwind_r[32]+alpha[26]*fUpwind_r[26]+alpha[25]*fUpwind_r[25]+alpha[22]*fUpwind_r[22]+alpha[21]*fUpwind_r[21]+alpha[20]*fUpwind_r[20]+alpha[19]*fUpwind_r[19]+alpha[16]*fUpwind_r[16]+alpha[15]*fUpwind_r[15]+alpha[12]*fUpwind_r[12]+alpha[11]*fUpwind_r[11]+alpha[9]*fUpwind_r[9]+alpha[8]*fUpwind_r[8]+alpha[7]*fUpwind_r[7]+alpha[6]*fUpwind_r[6]+alpha[5]*fUpwind_r[5]+alpha[4]*fUpwind_r[4]+alpha[3]*fUpwind_r[3]+alpha[2]*fUpwind_r[2]+alpha[1]*fUpwind_r[1]+alpha[0]*fUpwind_r[0]); 
  Ghat_r[1] += 0.2500000000000001*(alpha[26]*fUpwind_r[36]+fUpwind_r[26]*alpha[36])+0.223606797749979*(alpha[16]*fUpwind_r[35]+fUpwind_r[16]*alpha[35])+0.2500000000000001*(alpha[22]*fUpwind_r[33]+fUpwind_r[22]*alpha[33])+0.223606797749979*(alpha[15]*fUpwind_r[32]+fUpwind_r[15]*alpha[32])+0.223606797749979*(alpha[8]*fUpwind_r[25]+fUpwind_r[8]*alpha[25]+alpha[6]*fUpwind_r[21]+fUpwind_r[6]*alpha[21])+0.2500000000000001*(alpha[12]*fUpwind_r[20]+fUpwind_r[12]*alpha[20])+0.223606797749979*(alpha[5]*fUpwind_r[19]+fUpwind_r[5]*alpha[19])+0.25*(alpha[9]*fUpwind_r[16]+fUpwind_r[9]*alpha[16]+alpha[7]*fUpwind_r[15]+fUpwind_r[7]*alpha[15])+0.223606797749979*(alpha[1]*fUpwind_r[11]+fUpwind_r[1]*alpha[11])+0.25*(alpha[4]*fUpwind_r[8]+fUpwind_r[4]*alpha[8]+alpha[3]*fUpwind_r[6]+fUpwind_r[3]*alpha[6]+alpha[2]*fUpwind_r[5]+fUpwind_r[2]*alpha[5]+alpha[0]*fUpwind_r[1]+fUpwind_r[0]*alpha[1]); 
  Ghat_r[2] += 0.223606797749979*(alpha[16]*fUpwind_r[36]+fUpwind_r[16]*alpha[36])+0.2500000000000001*(alpha[25]*fUpwind_r[35]+fUpwind_r[25]*alpha[35])+0.223606797749979*(alpha[15]*fUpwind_r[33]+fUpwind_r[15]*alpha[33])+0.2500000000000001*(alpha[21]*fUpwind_r[32]+fUpwind_r[21]*alpha[32])+0.223606797749979*(alpha[9]*fUpwind_r[26]+fUpwind_r[9]*alpha[26]+alpha[7]*fUpwind_r[22]+fUpwind_r[7]*alpha[22]+alpha[5]*fUpwind_r[20]+fUpwind_r[5]*alpha[20])+0.2500000000000001*(alpha[11]*fUpwind_r[19]+fUpwind_r[11]*alpha[19])+0.25*(alpha[8]*fUpwind_r[16]+fUpwind_r[8]*alpha[16]+alpha[6]*fUpwind_r[15]+fUpwind_r[6]*alpha[15])+0.223606797749979*(alpha[2]*fUpwind_r[12]+fUpwind_r[2]*alpha[12])+0.25*(alpha[4]*fUpwind_r[9]+fUpwind_r[4]*alpha[9]+alpha[3]*fUpwind_r[7]+fUpwind_r[3]*alpha[7]+alpha[1]*fUpwind_r[5]+fUpwind_r[1]*alpha[5]+alpha[0]*fUpwind_r[2]+fUpwind_r[0]*alpha[2]); 
  Ghat_r[3] += 0.2500000000000001*(alpha[36]*fUpwind_r[45]+alpha[35]*fUpwind_r[44]+alpha[26]*fUpwind_r[38]+alpha[25]*fUpwind_r[37])+0.223606797749979*alpha[15]*fUpwind_r[34]+0.2500000000000001*(alpha[20]*fUpwind_r[33]+fUpwind_r[20]*alpha[33]+alpha[19]*fUpwind_r[32]+fUpwind_r[19]*alpha[32])+0.25*alpha[16]*fUpwind_r[31]+0.223606797749979*(alpha[7]*fUpwind_r[24]+alpha[6]*fUpwind_r[23])+0.2500000000000001*(alpha[12]*fUpwind_r[22]+fUpwind_r[12]*alpha[22]+alpha[11]*fUpwind_r[21]+fUpwind_r[11]*alpha[21])+0.25*(alpha[9]*fUpwind_r[18]+alpha[8]*fUpwind_r[17]+alpha[5]*fUpwind_r[15]+fUpwind_r[5]*alpha[15])+0.223606797749979*alpha[3]*fUpwind_r[13]+0.25*(alpha[4]*fUpwind_r[10]+alpha[2]*fUpwind_r[7]+fUpwind_r[2]*alpha[7]+alpha[1]*fUpwind_r[6]+fUpwind_r[1]*alpha[6]+alpha[0]*fUpwind_r[3]+fUpwind_r[0]*alpha[3]); 
  Ghat_r[4] += 0.2500000000000001*(alpha[33]*fUpwind_r[45]+alpha[32]*fUpwind_r[44])+0.223606797749979*alpha[16]*fUpwind_r[41]+0.2500000000000001*(alpha[22]*fUpwind_r[38]+alpha[21]*fUpwind_r[37]+alpha[20]*fUpwind_r[36]+fUpwind_r[20]*alpha[36]+alpha[19]*fUpwind_r[35]+fUpwind_r[19]*alpha[35])+0.25*alpha[15]*fUpwind_r[31]+0.223606797749979*(alpha[9]*fUpwind_r[29]+alpha[8]*fUpwind_r[28])+0.2500000000000001*(alpha[12]*fUpwind_r[26]+fUpwind_r[12]*alpha[26]+alpha[11]*fUpwind_r[25]+fUpwind_r[11]*alpha[25])+0.25*(alpha[7]*fUpwind_r[18]+alpha[6]*fUpwind_r[17]+alpha[5]*fUpwind_r[16]+fUpwind_r[5]*alpha[16])+0.223606797749979*alpha[4]*fUpwind_r[14]+0.25*(alpha[3]*fUpwind_r[10]+alpha[2]*fUpwind_r[9]+fUpwind_r[2]*alpha[9]+alpha[1]*fUpwind_r[8]+fUpwind_r[1]*alpha[8]+alpha[0]*fUpwind_r[4]+fUpwind_r[0]*alpha[4]); 
  Ghat_r[5] += (0.2*alpha[35]+0.223606797749979*alpha[9])*fUpwind_r[36]+0.2*fUpwind_r[35]*alpha[36]+0.223606797749979*(fUpwind_r[9]*alpha[36]+alpha[8]*fUpwind_r[35]+fUpwind_r[8]*alpha[35])+(0.2*alpha[32]+0.223606797749979*alpha[7])*fUpwind_r[33]+0.2*fUpwind_r[32]*alpha[33]+0.223606797749979*(fUpwind_r[7]*alpha[33]+alpha[6]*fUpwind_r[32]+fUpwind_r[6]*alpha[32])+0.223606797749979*(alpha[16]*fUpwind_r[26]+fUpwind_r[16]*alpha[26]+alpha[16]*fUpwind_r[25]+fUpwind_r[16]*alpha[25]+alpha[15]*fUpwind_r[22]+fUpwind_r[15]*alpha[22]+alpha[15]*fUpwind_r[21]+fUpwind_r[15]*alpha[21])+(0.2*alpha[19]+0.223606797749979*alpha[2])*fUpwind_r[20]+0.2*fUpwind_r[19]*alpha[20]+0.223606797749979*(fUpwind_r[2]*alpha[20]+alpha[1]*fUpwind_r[19]+fUpwind_r[1]*alpha[19])+0.25*(alpha[4]*fUpwind_r[16]+fUpwind_r[4]*alpha[16]+alpha[3]*fUpwind_r[15]+fUpwind_r[3]*alpha[15])+0.223606797749979*(alpha[5]*fUpwind_r[12]+fUpwind_r[5]*alpha[12]+alpha[5]*fUpwind_r[11]+fUpwind_r[5]*alpha[11])+0.25*(alpha[8]*fUpwind_r[9]+fUpwind_r[8]*alpha[9]+alpha[6]*fUpwind_r[7]+fUpwind_r[6]*alpha[7]+alpha[0]*fUpwind_r[5]+fUpwind_r[0]*alpha[5]+alpha[1]*fUpwind_r[2]+fUpwind_r[1]*alpha[2]); 
  Ghat_r[6] += 0.25*alpha[26]*fUpwind_r[45]+0.223606797749979*alpha[16]*fUpwind_r[44]+0.25*alpha[36]*fUpwind_r[38]+0.223606797749979*(alpha[8]*fUpwind_r[37]+fUpwind_r[31]*alpha[35])+(0.2*alpha[32]+0.223606797749979*alpha[7])*fUpwind_r[34]+0.25*(alpha[12]*fUpwind_r[33]+fUpwind_r[12]*alpha[33])+0.223606797749979*(alpha[5]*fUpwind_r[32]+fUpwind_r[5]*alpha[32])+0.25*alpha[9]*fUpwind_r[31]+0.223606797749979*(fUpwind_r[17]*alpha[25]+alpha[15]*fUpwind_r[24])+(0.2*alpha[21]+0.223606797749979*alpha[3])*fUpwind_r[23]+0.25*(alpha[20]*fUpwind_r[22]+fUpwind_r[20]*alpha[22])+0.223606797749979*(alpha[1]*fUpwind_r[21]+fUpwind_r[1]*alpha[21]+alpha[15]*fUpwind_r[19]+fUpwind_r[15]*alpha[19])+0.25*(alpha[16]*fUpwind_r[18]+alpha[4]*fUpwind_r[17]+alpha[2]*fUpwind_r[15]+fUpwind_r[2]*alpha[15])+0.223606797749979*(alpha[6]*(fUpwind_r[13]+fUpwind_r[11])+fUpwind_r[6]*alpha[11])+0.25*(alpha[8]*fUpwind_r[10]+alpha[5]*fUpwind_r[7]+fUpwind_r[5]*alpha[7]+alpha[0]*fUpwind_r[6]+fUpwind_r[0]*alpha[6]+alpha[1]*fUpwind_r[3]+fUpwind_r[1]*alpha[3]); 
  Ghat_r[7] += 0.223606797749979*alpha[16]*fUpwind_r[45]+0.25*alpha[25]*fUpwind_r[44]+0.223606797749979*alpha[9]*fUpwind_r[38]+0.25*alpha[35]*fUpwind_r[37]+0.223606797749979*fUpwind_r[31]*alpha[36]+0.2*alpha[33]*fUpwind_r[34]+0.223606797749979*(alpha[6]*fUpwind_r[34]+alpha[5]*fUpwind_r[33]+fUpwind_r[5]*alpha[33])+0.25*(alpha[11]*fUpwind_r[32]+fUpwind_r[11]*alpha[32]+alpha[8]*fUpwind_r[31])+0.223606797749979*fUpwind_r[18]*alpha[26]+0.2*alpha[22]*fUpwind_r[24]+0.223606797749979*(alpha[3]*fUpwind_r[24]+alpha[15]*fUpwind_r[23]+alpha[2]*fUpwind_r[22]+fUpwind_r[2]*alpha[22])+0.25*(alpha[19]*fUpwind_r[21]+fUpwind_r[19]*alpha[21])+0.223606797749979*(alpha[15]*fUpwind_r[20]+fUpwind_r[15]*alpha[20])+0.25*(alpha[4]*fUpwind_r[18]+alpha[16]*fUpwind_r[17]+alpha[1]*fUpwind_r[15]+fUpwind_r[1]*alpha[15])+0.223606797749979*(alpha[7]*(fUpwind_r[13]+fUpwind_r[12])+fUpwind_r[7]*alpha[12])+0.25*(alpha[9]*fUpwind_r[10]+alpha[0]*fUpwind_r[7]+fUpwind_r[0]*alpha[7]+alpha[5]*fUpwind_r[6]+fUpwind_r[5]*alpha[6]+alpha[2]*fUpwind_r[3]+fUpwind_r[2]*alpha[3]); 
  Ghat_r[8] += 0.25*alpha[22]*fUpwind_r[45]+0.223606797749979*alpha[15]*fUpwind_r[44]+(0.2*alpha[35]+0.223606797749979*alpha[9])*fUpwind_r[41]+0.25*alpha[33]*fUpwind_r[38]+0.223606797749979*alpha[6]*fUpwind_r[37]+0.25*(alpha[12]*fUpwind_r[36]+fUpwind_r[12]*alpha[36])+0.223606797749979*(alpha[5]*fUpwind_r[35]+fUpwind_r[5]*alpha[35])+fUpwind_r[31]*(0.223606797749979*alpha[32]+0.25*alpha[7])+0.223606797749979*alpha[16]*fUpwind_r[29]+(0.2*alpha[25]+0.223606797749979*alpha[4])*fUpwind_r[28]+0.25*(alpha[20]*fUpwind_r[26]+fUpwind_r[20]*alpha[26])+0.223606797749979*(alpha[1]*fUpwind_r[25]+fUpwind_r[1]*alpha[25]+fUpwind_r[17]*alpha[21]+alpha[16]*fUpwind_r[19]+fUpwind_r[16]*alpha[19])+0.25*(alpha[15]*fUpwind_r[18]+alpha[3]*fUpwind_r[17]+alpha[2]*fUpwind_r[16]+fUpwind_r[2]*alpha[16])+0.223606797749979*(alpha[8]*(fUpwind_r[14]+fUpwind_r[11])+fUpwind_r[8]*alpha[11])+0.25*(alpha[6]*fUpwind_r[10]+alpha[5]*fUpwind_r[9]+fUpwind_r[5]*alpha[9]+alpha[0]*fUpwind_r[8]+fUpwind_r[0]*alpha[8]+alpha[1]*fUpwind_r[4]+fUpwind_r[1]*alpha[4]); 
  Ghat_r[9] += 0.223606797749979*alpha[15]*fUpwind_r[45]+0.25*alpha[21]*fUpwind_r[44]+0.2*alpha[36]*fUpwind_r[41]+0.223606797749979*(alpha[8]*fUpwind_r[41]+alpha[7]*fUpwind_r[38])+0.25*alpha[32]*fUpwind_r[37]+0.223606797749979*(alpha[5]*fUpwind_r[36]+fUpwind_r[5]*alpha[36])+0.25*(alpha[11]*fUpwind_r[35]+fUpwind_r[11]*alpha[35])+fUpwind_r[31]*(0.223606797749979*alpha[33]+0.25*alpha[6])+0.2*alpha[26]*fUpwind_r[29]+0.223606797749979*(alpha[4]*fUpwind_r[29]+alpha[16]*fUpwind_r[28]+alpha[2]*fUpwind_r[26]+fUpwind_r[2]*alpha[26])+0.25*(alpha[19]*fUpwind_r[25]+fUpwind_r[19]*alpha[25])+0.223606797749979*(fUpwind_r[18]*alpha[22]+alpha[16]*fUpwind_r[20]+fUpwind_r[16]*alpha[20])+0.25*(alpha[3]*fUpwind_r[18]+alpha[15]*fUpwind_r[17]+alpha[1]*fUpwind_r[16]+fUpwind_r[1]*alpha[16])+0.223606797749979*(alpha[9]*(fUpwind_r[14]+fUpwind_r[12])+fUpwind_r[9]*alpha[12])+0.25*(alpha[7]*fUpwind_r[10]+alpha[0]*fUpwind_r[9]+fUpwind_r[0]*alpha[9]+alpha[5]*fUpwind_r[8]+fUpwind_r[5]*alpha[8]+alpha[2]*fUpwind_r[4]+fUpwind_r[2]*alpha[4]); 
  Ghat_r[10] += 0.223606797749979*(alpha[16]*fUpwind_r[47]+alpha[15]*fUpwind_r[46])+0.25*(alpha[20]*fUpwind_r[45]+alpha[19]*fUpwind_r[44])+0.223606797749979*(alpha[9]*fUpwind_r[43]+alpha[8]*fUpwind_r[42]+alpha[7]*fUpwind_r[40]+alpha[6]*fUpwind_r[39])+0.25*(alpha[12]*fUpwind_r[38]+alpha[11]*fUpwind_r[37]+alpha[33]*fUpwind_r[36]+fUpwind_r[33]*alpha[36]+alpha[32]*fUpwind_r[35]+fUpwind_r[32]*alpha[35]+alpha[5]*fUpwind_r[31])+0.223606797749979*(alpha[4]*fUpwind_r[30]+alpha[3]*fUpwind_r[27])+0.25*(alpha[22]*fUpwind_r[26]+fUpwind_r[22]*alpha[26]+alpha[21]*fUpwind_r[25]+fUpwind_r[21]*alpha[25]+alpha[2]*fUpwind_r[18]+alpha[1]*fUpwind_r[17]+alpha[15]*fUpwind_r[16]+fUpwind_r[15]*alpha[16]+alpha[0]*fUpwind_r[10]+alpha[7]*fUpwind_r[9]+fUpwind_r[7]*alpha[9]+alpha[6]*fUpwind_r[8]+fUpwind_r[6]*alpha[8]+alpha[3]*fUpwind_r[4]+fUpwind_r[3]*alpha[4]); 
  Ghat_r[11] += 0.223606797749979*alpha[36]*fUpwind_r[36]+0.159719141249985*alpha[35]*fUpwind_r[35]+0.25*(alpha[9]*fUpwind_r[35]+fUpwind_r[9]*alpha[35])+0.223606797749979*alpha[33]*fUpwind_r[33]+0.159719141249985*alpha[32]*fUpwind_r[32]+0.25*(alpha[7]*fUpwind_r[32]+fUpwind_r[7]*alpha[32])+0.159719141249985*alpha[25]*fUpwind_r[25]+0.2500000000000001*(alpha[4]*fUpwind_r[25]+fUpwind_r[4]*alpha[25])+0.159719141249985*alpha[21]*fUpwind_r[21]+0.2500000000000001*(alpha[3]*fUpwind_r[21]+fUpwind_r[3]*alpha[21])+0.223606797749979*alpha[20]*fUpwind_r[20]+0.159719141249985*alpha[19]*fUpwind_r[19]+0.2500000000000001*(alpha[2]*fUpwind_r[19]+fUpwind_r[2]*alpha[19])+0.223606797749979*(alpha[16]*fUpwind_r[16]+alpha[15]*fUpwind_r[15])+0.159719141249985*alpha[11]*fUpwind_r[11]+0.25*(alpha[0]*fUpwind_r[11]+fUpwind_r[0]*alpha[11])+0.223606797749979*(alpha[8]*fUpwind_r[8]+alpha[6]*fUpwind_r[6]+alpha[5]*fUpwind_r[5]+alpha[1]*fUpwind_r[1]); 
  Ghat_r[12] += 0.159719141249985*alpha[36]*fUpwind_r[36]+0.25*(alpha[8]*fUpwind_r[36]+fUpwind_r[8]*alpha[36])+0.223606797749979*alpha[35]*fUpwind_r[35]+0.159719141249985*alpha[33]*fUpwind_r[33]+0.25*(alpha[6]*fUpwind_r[33]+fUpwind_r[6]*alpha[33])+0.223606797749979*alpha[32]*fUpwind_r[32]+0.159719141249985*alpha[26]*fUpwind_r[26]+0.2500000000000001*(alpha[4]*fUpwind_r[26]+fUpwind_r[4]*alpha[26])+0.159719141249985*alpha[22]*fUpwind_r[22]+0.2500000000000001*(alpha[3]*fUpwind_r[22]+fUpwind_r[3]*alpha[22])+0.159719141249985*alpha[20]*fUpwind_r[20]+0.2500000000000001*(alpha[1]*fUpwind_r[20]+fUpwind_r[1]*alpha[20])+0.223606797749979*(alpha[19]*fUpwind_r[19]+alpha[16]*fUpwind_r[16]+alpha[15]*fUpwind_r[15])+0.159719141249985*alpha[12]*fUpwind_r[12]+0.25*(alpha[0]*fUpwind_r[12]+fUpwind_r[0]*alpha[12])+0.223606797749979*(alpha[9]*fUpwind_r[9]+alpha[7]*fUpwind_r[7]+alpha[5]*fUpwind_r[5]+alpha[2]*fUpwind_r[2]); 
  Ghat_r[13] += 0.2500000000000001*alpha[16]*fUpwind_r[46]+0.25*(alpha[9]*fUpwind_r[40]+alpha[8]*fUpwind_r[39]+alpha[5]*fUpwind_r[34])+0.223606797749979*(alpha[33]*fUpwind_r[33]+alpha[32]*fUpwind_r[32])+0.2500000000000001*(alpha[4]*fUpwind_r[27]+alpha[2]*fUpwind_r[24]+alpha[1]*fUpwind_r[23])+0.223606797749979*(alpha[22]*fUpwind_r[22]+alpha[21]*fUpwind_r[21]+alpha[15]*fUpwind_r[15])+0.25*alpha[0]*fUpwind_r[13]+0.223606797749979*(alpha[7]*fUpwind_r[7]+alpha[6]*fUpwind_r[6]+alpha[3]*fUpwind_r[3]); 
  Ghat_r[14] += 0.2500000000000001*alpha[15]*fUpwind_r[47]+0.25*(alpha[7]*fUpwind_r[43]+alpha[6]*fUpwind_r[42]+alpha[5]*fUpwind_r[41])+0.223606797749979*(alpha[36]*fUpwind_r[36]+alpha[35]*fUpwind_r[35])+0.2500000000000001*(alpha[3]*fUpwind_r[30]+alpha[2]*fUpwind_r[29]+alpha[1]*fUpwind_r[28])+0.223606797749979*(alpha[26]*fUpwind_r[26]+alpha[25]*fUpwind_r[25]+alpha[16]*fUpwind_r[16])+0.25*alpha[0]*fUpwind_r[14]+0.223606797749979*(alpha[9]*fUpwind_r[9]+alpha[8]*fUpwind_r[8]+alpha[4]*fUpwind_r[4]); 
  Ghat_r[15] += (0.2*alpha[35]+0.223606797749979*alpha[9])*fUpwind_r[45]+(0.2*alpha[36]+0.223606797749979*alpha[8])*fUpwind_r[44]+0.223606797749979*(alpha[16]*(fUpwind_r[38]+fUpwind_r[37])+fUpwind_r[18]*alpha[36]+fUpwind_r[17]*alpha[35])+(0.2*(alpha[22]+alpha[21])+0.223606797749979*alpha[3])*fUpwind_r[34]+(0.2*alpha[19]+0.223606797749979*alpha[2])*fUpwind_r[33]+(0.2*(fUpwind_r[24]+fUpwind_r[19])+0.223606797749979*fUpwind_r[2])*alpha[33]+(0.2*alpha[20]+0.223606797749979*alpha[1])*fUpwind_r[32]+(0.2*(fUpwind_r[23]+fUpwind_r[20])+0.223606797749979*fUpwind_r[1])*alpha[32]+(0.223606797749979*(alpha[26]+alpha[25])+0.25*alpha[4])*fUpwind_r[31]+0.223606797749979*(alpha[6]*fUpwind_r[24]+alpha[7]*fUpwind_r[23]+alpha[5]*fUpwind_r[22]+fUpwind_r[5]*alpha[22]+alpha[5]*fUpwind_r[21]+fUpwind_r[5]*alpha[21]+alpha[7]*fUpwind_r[20]+fUpwind_r[7]*alpha[20]+alpha[6]*fUpwind_r[19]+fUpwind_r[6]*alpha[19])+0.25*(alpha[8]*fUpwind_r[18]+alpha[9]*fUpwind_r[17]+fUpwind_r[10]*alpha[16])+(0.223606797749979*(alpha[12]+alpha[11])+0.25*alpha[0])*fUpwind_r[15]+0.223606797749979*(fUpwind_r[13]+fUpwind_r[12]+fUpwind_r[11])*alpha[15]+0.25*(fUpwind_r[0]*alpha[15]+alpha[1]*fUpwind_r[7]+fUpwind_r[1]*alpha[7]+alpha[2]*fUpwind_r[6]+fUpwind_r[2]*alpha[6]+alpha[3]*fUpwind_r[5]+fUpwind_r[3]*alpha[5]); 
  Ghat_r[16] += (0.2*alpha[32]+0.223606797749979*alpha[7])*fUpwind_r[45]+(0.2*alpha[33]+0.223606797749979*alpha[6])*fUpwind_r[44]+0.2*(alpha[26]+alpha[25])*fUpwind_r[41]+0.223606797749979*(alpha[4]*fUpwind_r[41]+alpha[15]*(fUpwind_r[38]+fUpwind_r[37]))+(0.2*alpha[19]+0.223606797749979*alpha[2])*fUpwind_r[36]+(0.2*(fUpwind_r[29]+fUpwind_r[19])+0.223606797749979*fUpwind_r[2])*alpha[36]+(0.2*alpha[20]+0.223606797749979*alpha[1])*fUpwind_r[35]+0.2*(fUpwind_r[28]+fUpwind_r[20])*alpha[35]+0.223606797749979*(fUpwind_r[1]*alpha[35]+fUpwind_r[18]*alpha[33]+fUpwind_r[17]*alpha[32])+(0.223606797749979*(alpha[22]+alpha[21])+0.25*alpha[3])*fUpwind_r[31]+0.223606797749979*(alpha[8]*fUpwind_r[29]+alpha[9]*fUpwind_r[28]+alpha[5]*fUpwind_r[26]+fUpwind_r[5]*alpha[26]+alpha[5]*fUpwind_r[25]+fUpwind_r[5]*alpha[25]+alpha[9]*fUpwind_r[20]+fUpwind_r[9]*alpha[20]+alpha[8]*fUpwind_r[19]+fUpwind_r[8]*alpha[19])+0.25*(alpha[6]*fUpwind_r[18]+alpha[7]*fUpwind_r[17])+(0.223606797749979*(alpha[12]+alpha[11])+0.25*alpha[0])*fUpwind_r[16]+0.223606797749979*(fUpwind_r[14]+fUpwind_r[12]+fUpwind_r[11])*alpha[16]+0.25*(fUpwind_r[0]*alpha[16]+fUpwind_r[10]*alpha[15]+alpha[1]*fUpwind_r[9]+fUpwind_r[1]*alpha[9]+alpha[2]*fUpwind_r[8]+fUpwind_r[2]*alpha[8]+alpha[4]*fUpwind_r[5]+fUpwind_r[4]*alpha[5]); 
  Ghat_r[17] += (0.2*alpha[35]+0.223606797749979*alpha[9])*fUpwind_r[47]+(0.2*alpha[32]+0.223606797749979*alpha[7])*fUpwind_r[46]+0.2500000000000001*alpha[12]*fUpwind_r[45]+0.223606797749979*alpha[5]*fUpwind_r[44]+0.223606797749979*alpha[16]*fUpwind_r[43]+0.2*alpha[25]*fUpwind_r[42]+0.223606797749979*(alpha[4]*fUpwind_r[42]+alpha[15]*fUpwind_r[40])+(0.2*alpha[21]+0.223606797749979*alpha[3])*fUpwind_r[39]+0.2500000000000001*alpha[20]*fUpwind_r[38]+0.223606797749979*alpha[1]*fUpwind_r[37]+0.2500000000000001*(alpha[22]*fUpwind_r[36]+fUpwind_r[22]*alpha[36])+0.223606797749979*(alpha[15]*fUpwind_r[35]+fUpwind_r[15]*alpha[35])+0.2500000000000001*(alpha[26]*fUpwind_r[33]+fUpwind_r[26]*alpha[33])+0.223606797749979*(alpha[16]*fUpwind_r[32]+fUpwind_r[16]*alpha[32])+(0.223606797749979*alpha[19]+0.25*alpha[2])*fUpwind_r[31]+0.223606797749979*(alpha[8]*fUpwind_r[30]+alpha[6]*(fUpwind_r[27]+fUpwind_r[25])+fUpwind_r[6]*alpha[25]+alpha[8]*fUpwind_r[21]+fUpwind_r[8]*alpha[21])+0.25*alpha[5]*fUpwind_r[18]+0.223606797749979*alpha[11]*fUpwind_r[17]+0.25*(alpha[0]*fUpwind_r[17]+alpha[7]*fUpwind_r[16]+fUpwind_r[7]*alpha[16]+alpha[9]*fUpwind_r[15]+fUpwind_r[9]*alpha[15]+alpha[1]*fUpwind_r[10]+alpha[3]*fUpwind_r[8]+fUpwind_r[3]*alpha[8]+alpha[4]*fUpwind_r[6]+fUpwind_r[4]*alpha[6]); 
  Ghat_r[18] += (0.2*alpha[36]+0.223606797749979*alpha[8])*fUpwind_r[47]+0.2*alpha[33]*fUpwind_r[46]+0.223606797749979*(alpha[6]*fUpwind_r[46]+alpha[5]*fUpwind_r[45])+0.2500000000000001*alpha[11]*fUpwind_r[44]+0.2*alpha[26]*fUpwind_r[43]+0.223606797749979*(alpha[4]*fUpwind_r[43]+alpha[16]*fUpwind_r[42])+0.2*alpha[22]*fUpwind_r[40]+0.223606797749979*(alpha[3]*fUpwind_r[40]+alpha[15]*fUpwind_r[39]+alpha[2]*fUpwind_r[38])+0.2500000000000001*alpha[19]*fUpwind_r[37]+0.223606797749979*(alpha[15]*fUpwind_r[36]+fUpwind_r[15]*alpha[36])+0.2500000000000001*(alpha[21]*fUpwind_r[35]+fUpwind_r[21]*alpha[35])+0.223606797749979*(alpha[16]*fUpwind_r[33]+fUpwind_r[16]*alpha[33])+0.2500000000000001*(alpha[25]*fUpwind_r[32]+fUpwind_r[25]*alpha[32])+(0.223606797749979*alpha[20]+0.25*alpha[1])*fUpwind_r[31]+0.223606797749979*(alpha[9]*fUpwind_r[30]+alpha[7]*(fUpwind_r[27]+fUpwind_r[26])+fUpwind_r[7]*alpha[26]+alpha[9]*fUpwind_r[22]+fUpwind_r[9]*alpha[22])+0.223606797749979*alpha[12]*fUpwind_r[18]+0.25*(alpha[0]*fUpwind_r[18]+alpha[5]*fUpwind_r[17]+alpha[6]*fUpwind_r[16]+fUpwind_r[6]*alpha[16]+alpha[8]*fUpwind_r[15]+fUpwind_r[8]*alpha[15]+alpha[2]*fUpwind_r[10]+alpha[3]*fUpwind_r[9]+fUpwind_r[3]*alpha[9]+alpha[4]*fUpwind_r[7]+fUpwind_r[4]*alpha[7]); 
  Ghat_r[19] += 0.2*(alpha[16]*fUpwind_r[36]+fUpwind_r[16]*alpha[36])+(0.223606797749979*alpha[26]+0.159719141249985*alpha[25]+0.2500000000000001*alpha[4])*fUpwind_r[35]+(0.223606797749979*fUpwind_r[26]+0.159719141249985*fUpwind_r[25]+0.2500000000000001*fUpwind_r[4])*alpha[35]+0.2*(alpha[15]*fUpwind_r[33]+fUpwind_r[15]*alpha[33])+(0.223606797749979*alpha[22]+0.159719141249985*alpha[21]+0.2500000000000001*alpha[3])*fUpwind_r[32]+(0.223606797749979*fUpwind_r[22]+0.159719141249985*fUpwind_r[21]+0.2500000000000001*fUpwind_r[3])*alpha[32]+0.25*(alpha[9]*fUpwind_r[25]+fUpwind_r[9]*alpha[25]+alpha[7]*fUpwind_r[21]+fUpwind_r[7]*alpha[21])+0.2*(alpha[5]*fUpwind_r[20]+fUpwind_r[5]*alpha[20])+(0.223606797749979*alpha[12]+0.159719141249985*alpha[11]+0.25*alpha[0])*fUpwind_r[19]+(0.223606797749979*fUpwind_r[12]+0.159719141249985*fUpwind_r[11]+0.25*fUpwind_r[0])*alpha[19]+0.223606797749979*(alpha[8]*fUpwind_r[16]+fUpwind_r[8]*alpha[16]+alpha[6]*fUpwind_r[15]+fUpwind_r[6]*alpha[15])+0.2500000000000001*(alpha[2]*fUpwind_r[11]+fUpwind_r[2]*alpha[11])+0.223606797749979*(alpha[1]*fUpwind_r[5]+fUpwind_r[1]*alpha[5]); 
  Ghat_r[20] += (0.159719141249985*alpha[26]+0.223606797749979*alpha[25]+0.2500000000000001*alpha[4])*fUpwind_r[36]+(0.159719141249985*fUpwind_r[26]+0.223606797749979*fUpwind_r[25]+0.2500000000000001*fUpwind_r[4])*alpha[36]+0.2*(alpha[16]*fUpwind_r[35]+fUpwind_r[16]*alpha[35])+(0.159719141249985*alpha[22]+0.223606797749979*alpha[21]+0.2500000000000001*alpha[3])*fUpwind_r[33]+(0.159719141249985*fUpwind_r[22]+0.223606797749979*fUpwind_r[21]+0.2500000000000001*fUpwind_r[3])*alpha[33]+0.2*(alpha[15]*fUpwind_r[32]+fUpwind_r[15]*alpha[32])+0.25*(alpha[8]*fUpwind_r[26]+fUpwind_r[8]*alpha[26]+alpha[6]*fUpwind_r[22]+fUpwind_r[6]*alpha[22])+(0.159719141249985*alpha[12]+0.223606797749979*alpha[11]+0.25*alpha[0])*fUpwind_r[20]+(0.159719141249985*fUpwind_r[12]+0.223606797749979*fUpwind_r[11]+0.25*fUpwind_r[0])*alpha[20]+0.2*(alpha[5]*fUpwind_r[19]+fUpwind_r[5]*alpha[19])+0.223606797749979*(alpha[9]*fUpwind_r[16]+fUpwind_r[9]*alpha[16]+alpha[7]*fUpwind_r[15]+fUpwind_r[7]*alpha[15])+0.2500000000000001*(alpha[1]*fUpwind_r[12]+fUpwind_r[1]*alpha[12])+0.223606797749979*(alpha[2]*fUpwind_r[5]+fUpwind_r[2]*alpha[5]); 
  Ghat_r[21] += 0.223606797749979*alpha[36]*fUpwind_r[45]+(0.159719141249985*alpha[35]+0.25*alpha[9])*fUpwind_r[44]+0.159719141249985*alpha[25]*fUpwind_r[37]+0.2500000000000001*(alpha[4]*fUpwind_r[37]+fUpwind_r[18]*alpha[35])+0.2*alpha[15]*fUpwind_r[34]+0.223606797749979*(alpha[20]*fUpwind_r[33]+fUpwind_r[20]*alpha[33])+(0.159719141249985*alpha[19]+0.2500000000000001*alpha[2])*fUpwind_r[32]+(0.223606797749979*fUpwind_r[24]+0.159719141249985*fUpwind_r[19]+0.2500000000000001*fUpwind_r[2])*alpha[32]+0.223606797749979*alpha[16]*fUpwind_r[31]+0.25*fUpwind_r[10]*alpha[25]+0.2*alpha[6]*fUpwind_r[23]+(0.159719141249985*alpha[11]+0.25*alpha[0])*fUpwind_r[21]+(0.223606797749979*fUpwind_r[13]+0.159719141249985*fUpwind_r[11])*alpha[21]+0.25*(fUpwind_r[0]*alpha[21]+alpha[7]*fUpwind_r[19]+fUpwind_r[7]*alpha[19])+0.223606797749979*(alpha[8]*fUpwind_r[17]+alpha[5]*fUpwind_r[15]+fUpwind_r[5]*alpha[15])+0.2500000000000001*(alpha[3]*fUpwind_r[11]+fUpwind_r[3]*alpha[11])+0.223606797749979*(alpha[1]*fUpwind_r[6]+fUpwind_r[1]*alpha[6]); 
  Ghat_r[22] += (0.159719141249985*alpha[36]+0.25*alpha[8])*fUpwind_r[45]+0.223606797749979*alpha[35]*fUpwind_r[44]+0.159719141249985*alpha[26]*fUpwind_r[38]+0.2500000000000001*(alpha[4]*fUpwind_r[38]+fUpwind_r[17]*alpha[36])+0.2*alpha[15]*fUpwind_r[34]+(0.159719141249985*alpha[20]+0.2500000000000001*alpha[1])*fUpwind_r[33]+(0.223606797749979*fUpwind_r[23]+0.159719141249985*fUpwind_r[20]+0.2500000000000001*fUpwind_r[1])*alpha[33]+0.223606797749979*(alpha[19]*fUpwind_r[32]+fUpwind_r[19]*alpha[32])+0.223606797749979*alpha[16]*fUpwind_r[31]+0.25*fUpwind_r[10]*alpha[26]+0.2*alpha[7]*fUpwind_r[24]+(0.159719141249985*alpha[12]+0.25*alpha[0])*fUpwind_r[22]+(0.223606797749979*fUpwind_r[13]+0.159719141249985*fUpwind_r[12])*alpha[22]+0.25*(fUpwind_r[0]*alpha[22]+alpha[6]*fUpwind_r[20]+fUpwind_r[6]*alpha[20])+0.223606797749979*(alpha[9]*fUpwind_r[18]+alpha[5]*fUpwind_r[15]+fUpwind_r[5]*alpha[15])+0.2500000000000001*(alpha[3]*fUpwind_r[12]+fUpwind_r[3]*alpha[12])+0.223606797749979*(alpha[2]*fUpwind_r[7]+fUpwind_r[2]*alpha[7]); 
  Ghat_r[23] += (0.223606797749979*alpha[35]+0.25*alpha[9])*fUpwind_r[46]+0.2500000000000001*alpha[16]*fUpwind_r[40]+(0.223606797749979*alpha[25]+0.2500000000000001*alpha[4])*fUpwind_r[39]+(0.223606797749979*alpha[19]+0.2500000000000001*alpha[2])*fUpwind_r[34]+0.223606797749979*(alpha[22]*fUpwind_r[33]+fUpwind_r[22]*alpha[33])+0.2*(alpha[15]*fUpwind_r[32]+fUpwind_r[15]*alpha[32])+0.25*(alpha[8]*fUpwind_r[27]+alpha[5]*fUpwind_r[24])+(0.223606797749979*alpha[11]+0.25*alpha[0])*fUpwind_r[23]+0.2*(alpha[6]*fUpwind_r[21]+fUpwind_r[6]*alpha[21])+0.223606797749979*(alpha[7]*fUpwind_r[15]+fUpwind_r[7]*alpha[15])+0.2500000000000001*alpha[1]*fUpwind_r[13]+0.223606797749979*(alpha[3]*fUpwind_r[6]+fUpwind_r[3]*alpha[6]); 
  Ghat_r[24] += (0.223606797749979*alpha[36]+0.25*alpha[8])*fUpwind_r[46]+0.223606797749979*alpha[26]*fUpwind_r[40]+0.2500000000000001*(alpha[4]*fUpwind_r[40]+alpha[16]*fUpwind_r[39])+(0.223606797749979*alpha[20]+0.2500000000000001*alpha[1])*fUpwind_r[34]+0.2*(alpha[15]*fUpwind_r[33]+fUpwind_r[15]*alpha[33])+0.223606797749979*(alpha[21]*fUpwind_r[32]+fUpwind_r[21]*alpha[32])+0.25*alpha[9]*fUpwind_r[27]+0.223606797749979*alpha[12]*fUpwind_r[24]+0.25*(alpha[0]*fUpwind_r[24]+alpha[5]*fUpwind_r[23])+0.2*(alpha[7]*fUpwind_r[22]+fUpwind_r[7]*alpha[22])+0.223606797749979*(alpha[6]*fUpwind_r[15]+fUpwind_r[6]*alpha[15])+0.2500000000000001*alpha[2]*fUpwind_r[13]+0.223606797749979*(alpha[3]*fUpwind_r[7]+fUpwind_r[3]*alpha[7]); 
  Ghat_r[25] += 0.223606797749979*alpha[33]*fUpwind_r[45]+(0.159719141249985*alpha[32]+0.25*alpha[7])*fUpwind_r[44]+0.2*alpha[16]*fUpwind_r[41]+(0.159719141249985*alpha[21]+0.2500000000000001*alpha[3])*fUpwind_r[37]+0.223606797749979*(alpha[20]*fUpwind_r[36]+fUpwind_r[20]*alpha[36])+(0.159719141249985*alpha[19]+0.2500000000000001*alpha[2])*fUpwind_r[35]+(0.223606797749979*fUpwind_r[29]+0.159719141249985*fUpwind_r[19])*alpha[35]+0.2500000000000001*(fUpwind_r[2]*alpha[35]+fUpwind_r[18]*alpha[32])+0.223606797749979*alpha[15]*fUpwind_r[31]+0.2*alpha[8]*fUpwind_r[28]+(0.159719141249985*alpha[11]+0.25*alpha[0])*fUpwind_r[25]+(0.223606797749979*fUpwind_r[14]+0.159719141249985*fUpwind_r[11])*alpha[25]+0.25*(fUpwind_r[0]*alpha[25]+fUpwind_r[10]*alpha[21]+alpha[9]*fUpwind_r[19]+fUpwind_r[9]*alpha[19])+0.223606797749979*(alpha[6]*fUpwind_r[17]+alpha[5]*fUpwind_r[16]+fUpwind_r[5]*alpha[16])+0.2500000000000001*(alpha[4]*fUpwind_r[11]+fUpwind_r[4]*alpha[11])+0.223606797749979*(alpha[1]*fUpwind_r[8]+fUpwind_r[1]*alpha[8]); 
  Ghat_r[26] += (0.159719141249985*alpha[33]+0.25*alpha[6])*fUpwind_r[45]+0.223606797749979*alpha[32]*fUpwind_r[44]+0.2*alpha[16]*fUpwind_r[41]+(0.159719141249985*alpha[22]+0.2500000000000001*alpha[3])*fUpwind_r[38]+(0.159719141249985*alpha[20]+0.2500000000000001*alpha[1])*fUpwind_r[36]+(0.223606797749979*fUpwind_r[28]+0.159719141249985*fUpwind_r[20]+0.2500000000000001*fUpwind_r[1])*alpha[36]+0.223606797749979*(alpha[19]*fUpwind_r[35]+fUpwind_r[19]*alpha[35])+0.2500000000000001*fUpwind_r[17]*alpha[33]+0.223606797749979*alpha[15]*fUpwind_r[31]+0.2*alpha[9]*fUpwind_r[29]+(0.159719141249985*alpha[12]+0.25*alpha[0])*fUpwind_r[26]+(0.223606797749979*fUpwind_r[14]+0.159719141249985*fUpwind_r[12])*alpha[26]+0.25*(fUpwind_r[0]*alpha[26]+fUpwind_r[10]*alpha[22]+alpha[8]*fUpwind_r[20]+fUpwind_r[8]*alpha[20])+0.223606797749979*(alpha[7]*fUpwind_r[18]+alpha[5]*fUpwind_r[16]+fUpwind_r[5]*alpha[16])+0.2500000000000001*(alpha[4]*fUpwind_r[12]+fUpwind_r[4]*alpha[12])+0.223606797749979*(alpha[2]*fUpwind_r[9]+fUpwind_r[2]*alpha[9]); 
  Ghat_r[27] += 0.25*alpha[5]*fUpwind_r[46]+0.223606797749979*(alpha[33]*fUpwind_r[45]+alpha[32]*fUpwind_r[44])+0.2500000000000001*(alpha[2]*fUpwind_r[40]+alpha[1]*fUpwind_r[39])+0.223606797749979*(alpha[22]*fUpwind_r[38]+alpha[21]*fUpwind_r[37])+0.2500000000000001*alpha[16]*fUpwind_r[34]+0.223606797749979*alpha[15]*fUpwind_r[31]+0.25*(alpha[0]*fUpwind_r[27]+alpha[9]*fUpwind_r[24]+alpha[8]*fUpwind_r[23])+0.223606797749979*(alpha[7]*fUpwind_r[18]+alpha[6]*fUpwind_r[17])+0.2500000000000001*alpha[4]*fUpwind_r[13]+0.223606797749979*alpha[3]*fUpwind_r[10]; 
  Ghat_r[28] += (0.223606797749979*alpha[32]+0.25*alpha[7])*fUpwind_r[47]+0.2500000000000001*alpha[15]*fUpwind_r[43]+(0.223606797749979*alpha[21]+0.2500000000000001*alpha[3])*fUpwind_r[42]+(0.223606797749979*alpha[19]+0.2500000000000001*alpha[2])*fUpwind_r[41]+0.223606797749979*(alpha[26]*fUpwind_r[36]+fUpwind_r[26]*alpha[36])+0.2*(alpha[16]*fUpwind_r[35]+fUpwind_r[16]*alpha[35])+0.25*(alpha[6]*fUpwind_r[30]+alpha[5]*fUpwind_r[29])+(0.223606797749979*alpha[11]+0.25*alpha[0])*fUpwind_r[28]+0.2*(alpha[8]*fUpwind_r[25]+fUpwind_r[8]*alpha[25])+0.223606797749979*(alpha[9]*fUpwind_r[16]+fUpwind_r[9]*alpha[16])+0.2500000000000001*alpha[1]*fUpwind_r[14]+0.223606797749979*(alpha[4]*fUpwind_r[8]+fUpwind_r[4]*alpha[8]); 
  Ghat_r[29] += (0.223606797749979*alpha[33]+0.25*alpha[6])*fUpwind_r[47]+0.223606797749979*alpha[22]*fUpwind_r[43]+0.2500000000000001*(alpha[3]*fUpwind_r[43]+alpha[15]*fUpwind_r[42])+(0.223606797749979*alpha[20]+0.2500000000000001*alpha[1])*fUpwind_r[41]+0.2*(alpha[16]*fUpwind_r[36]+fUpwind_r[16]*alpha[36])+0.223606797749979*(alpha[25]*fUpwind_r[35]+fUpwind_r[25]*alpha[35])+0.25*alpha[7]*fUpwind_r[30]+0.223606797749979*alpha[12]*fUpwind_r[29]+0.25*(alpha[0]*fUpwind_r[29]+alpha[5]*fUpwind_r[28])+0.2*(alpha[9]*fUpwind_r[26]+fUpwind_r[9]*alpha[26])+0.223606797749979*(alpha[8]*fUpwind_r[16]+fUpwind_r[8]*alpha[16])+0.2500000000000001*alpha[2]*fUpwind_r[14]+0.223606797749979*(alpha[4]*fUpwind_r[9]+fUpwind_r[4]*alpha[9]); 
  Ghat_r[30] += 0.25*alpha[5]*fUpwind_r[47]+0.223606797749979*(alpha[36]*fUpwind_r[45]+alpha[35]*fUpwind_r[44])+0.2500000000000001*(alpha[2]*fUpwind_r[43]+alpha[1]*fUpwind_r[42]+alpha[15]*fUpwind_r[41])+0.223606797749979*(alpha[26]*fUpwind_r[38]+alpha[25]*fUpwind_r[37])+0.223606797749979*alpha[16]*fUpwind_r[31]+0.25*(alpha[0]*fUpwind_r[30]+alpha[7]*fUpwind_r[29]+alpha[6]*fUpwind_r[28])+0.223606797749979*(alpha[9]*fUpwind_r[18]+alpha[8]*fUpwind_r[17])+0.2500000000000001*alpha[3]*fUpwind_r[14]+0.223606797749979*alpha[4]*fUpwind_r[10]; 
  Ghat_r[31] += (0.2*(alpha[26]+alpha[25])+0.223606797749979*alpha[4])*fUpwind_r[47]+(0.2*(alpha[22]+alpha[21])+0.223606797749979*alpha[3])*fUpwind_r[46]+(0.2*alpha[19]+0.223606797749979*alpha[2])*fUpwind_r[45]+(0.2*alpha[20]+0.223606797749979*alpha[1])*fUpwind_r[44]+(0.2*alpha[36]+0.223606797749979*alpha[8])*fUpwind_r[43]+(0.2*alpha[35]+0.223606797749979*alpha[9])*fUpwind_r[42]+(0.2*alpha[33]+0.223606797749979*alpha[6])*fUpwind_r[40]+0.2*alpha[32]*fUpwind_r[39]+0.223606797749979*(alpha[7]*fUpwind_r[39]+alpha[5]*(fUpwind_r[38]+fUpwind_r[37]))+(0.2*alpha[32]+0.223606797749979*alpha[7])*fUpwind_r[36]+(0.2*fUpwind_r[32]+0.223606797749979*fUpwind_r[7])*alpha[36]+(0.2*alpha[33]+0.223606797749979*alpha[6])*fUpwind_r[35]+0.2*fUpwind_r[33]*alpha[35]+0.223606797749979*(fUpwind_r[6]*alpha[35]+alpha[9]*fUpwind_r[33]+fUpwind_r[9]*alpha[33]+alpha[8]*fUpwind_r[32]+fUpwind_r[8]*alpha[32])+(0.223606797749979*(alpha[12]+alpha[11])+0.25*alpha[0])*fUpwind_r[31]+0.223606797749979*(alpha[16]*fUpwind_r[30]+alpha[15]*(fUpwind_r[27]+fUpwind_r[26])+fUpwind_r[15]*alpha[26]+alpha[15]*fUpwind_r[25]+fUpwind_r[15]*alpha[25]+alpha[16]*fUpwind_r[22]+fUpwind_r[16]*alpha[22]+alpha[16]*fUpwind_r[21]+fUpwind_r[16]*alpha[21]+fUpwind_r[18]*alpha[20]+fUpwind_r[17]*alpha[19])+0.25*(alpha[1]*fUpwind_r[18]+alpha[2]*fUpwind_r[17]+alpha[3]*fUpwind_r[16]+fUpwind_r[3]*alpha[16]+alpha[4]*fUpwind_r[15]+fUpwind_r[4]*alpha[15]+alpha[5]*fUpwind_r[10]+alpha[6]*fUpwind_r[9]+fUpwind_r[6]*alpha[9]+alpha[7]*fUpwind_r[8]+fUpwind_r[7]*alpha[8]); 
  Ghat_r[32] += 0.2*alpha[16]*fUpwind_r[45]+(0.223606797749979*alpha[26]+0.159719141249985*alpha[25]+0.2500000000000001*alpha[4])*fUpwind_r[44]+0.223606797749979*alpha[35]*fUpwind_r[38]+(0.159719141249985*alpha[35]+0.25*alpha[9])*fUpwind_r[37]+0.2*fUpwind_r[31]*alpha[36]+0.25*fUpwind_r[10]*alpha[35]+0.1788854381999831*alpha[33]*fUpwind_r[34]+0.2*(alpha[6]*fUpwind_r[34]+alpha[5]*fUpwind_r[33]+fUpwind_r[5]*alpha[33])+(0.223606797749979*alpha[12]+0.159719141249985*alpha[11]+0.25*alpha[0])*fUpwind_r[32]+(0.223606797749979*(fUpwind_r[13]+fUpwind_r[12])+0.159719141249985*fUpwind_r[11]+0.25*fUpwind_r[0])*alpha[32]+0.223606797749979*alpha[8]*fUpwind_r[31]+0.2500000000000001*fUpwind_r[18]*alpha[25]+0.223606797749979*alpha[21]*fUpwind_r[24]+0.2*alpha[15]*fUpwind_r[23]+0.223606797749979*(alpha[19]*fUpwind_r[22]+fUpwind_r[19]*alpha[22])+(0.159719141249985*alpha[19]+0.2500000000000001*alpha[2])*fUpwind_r[21]+(0.159719141249985*fUpwind_r[19]+0.2500000000000001*fUpwind_r[2])*alpha[21]+0.2*(alpha[15]*fUpwind_r[20]+fUpwind_r[15]*alpha[20])+0.2500000000000001*(alpha[3]*fUpwind_r[19]+fUpwind_r[3]*alpha[19])+0.223606797749979*(alpha[16]*fUpwind_r[17]+alpha[1]*fUpwind_r[15]+fUpwind_r[1]*alpha[15])+0.25*(alpha[7]*fUpwind_r[11]+fUpwind_r[7]*alpha[11])+0.223606797749979*(alpha[5]*fUpwind_r[6]+fUpwind_r[5]*alpha[6]); 
  Ghat_r[33] += (0.159719141249985*alpha[26]+0.223606797749979*alpha[25]+0.2500000000000001*alpha[4])*fUpwind_r[45]+0.2*alpha[16]*fUpwind_r[44]+(0.159719141249985*alpha[36]+0.25*alpha[8])*fUpwind_r[38]+alpha[36]*(0.223606797749979*fUpwind_r[37]+0.25*fUpwind_r[10])+0.2*fUpwind_r[31]*alpha[35]+(0.1788854381999831*alpha[32]+0.2*alpha[7])*fUpwind_r[34]+(0.159719141249985*alpha[12]+0.223606797749979*alpha[11]+0.25*alpha[0])*fUpwind_r[33]+(0.223606797749979*fUpwind_r[13]+0.159719141249985*fUpwind_r[12]+0.223606797749979*fUpwind_r[11]+0.25*fUpwind_r[0])*alpha[33]+0.2*(alpha[5]*fUpwind_r[32]+fUpwind_r[5]*alpha[32])+0.223606797749979*alpha[9]*fUpwind_r[31]+0.2500000000000001*fUpwind_r[17]*alpha[26]+0.2*alpha[15]*fUpwind_r[24]+0.223606797749979*alpha[22]*fUpwind_r[23]+(0.159719141249985*alpha[20]+0.2500000000000001*alpha[1])*fUpwind_r[22]+(0.159719141249985*fUpwind_r[20]+0.2500000000000001*fUpwind_r[1])*alpha[22]+0.223606797749979*(alpha[20]*fUpwind_r[21]+fUpwind_r[20]*alpha[21])+0.2500000000000001*(alpha[3]*fUpwind_r[20]+fUpwind_r[3]*alpha[20])+0.2*(alpha[15]*fUpwind_r[19]+fUpwind_r[15]*alpha[19])+0.223606797749979*(alpha[16]*fUpwind_r[18]+alpha[2]*fUpwind_r[15]+fUpwind_r[2]*alpha[15])+0.25*(alpha[6]*fUpwind_r[12]+fUpwind_r[6]*alpha[12])+0.223606797749979*(alpha[5]*fUpwind_r[7]+fUpwind_r[5]*alpha[7]); 
  Ghat_r[34] += (0.223606797749979*(alpha[26]+alpha[25])+0.2500000000000001*alpha[4])*fUpwind_r[46]+(0.223606797749979*alpha[36]+0.25*alpha[8])*fUpwind_r[40]+(0.223606797749979*alpha[35]+0.25*alpha[9])*fUpwind_r[39]+(0.223606797749979*(alpha[12]+alpha[11])+0.25*alpha[0])*fUpwind_r[34]+(0.1788854381999831*alpha[32]+0.2*alpha[7])*fUpwind_r[33]+0.1788854381999831*fUpwind_r[32]*alpha[33]+0.2*(fUpwind_r[7]*alpha[33]+alpha[6]*fUpwind_r[32]+fUpwind_r[6]*alpha[32])+0.2500000000000001*alpha[16]*fUpwind_r[27]+(0.223606797749979*alpha[20]+0.2500000000000001*alpha[1])*fUpwind_r[24]+(0.223606797749979*alpha[19]+0.2500000000000001*alpha[2])*fUpwind_r[23]+0.2*(alpha[15]*fUpwind_r[22]+fUpwind_r[15]*alpha[22]+alpha[15]*fUpwind_r[21]+fUpwind_r[15]*alpha[21])+0.223606797749979*(alpha[3]*fUpwind_r[15]+fUpwind_r[3]*alpha[15])+0.25*alpha[5]*fUpwind_r[13]+0.223606797749979*(alpha[6]*fUpwind_r[7]+fUpwind_r[6]*alpha[7]); 
  Ghat_r[35] += 0.2*alpha[15]*fUpwind_r[45]+(0.223606797749979*alpha[22]+0.159719141249985*alpha[21]+0.2500000000000001*alpha[3])*fUpwind_r[44]+(0.1788854381999831*alpha[36]+0.2*alpha[8])*fUpwind_r[41]+0.223606797749979*alpha[32]*fUpwind_r[38]+(0.159719141249985*alpha[32]+0.25*alpha[7])*fUpwind_r[37]+0.2*(alpha[5]*fUpwind_r[36]+fUpwind_r[5]*alpha[36])+(0.223606797749979*alpha[12]+0.159719141249985*alpha[11]+0.25*alpha[0])*fUpwind_r[35]+(0.223606797749979*(fUpwind_r[14]+fUpwind_r[12])+0.159719141249985*fUpwind_r[11]+0.25*fUpwind_r[0])*alpha[35]+0.2*fUpwind_r[31]*alpha[33]+0.25*fUpwind_r[10]*alpha[32]+0.223606797749979*(alpha[6]*fUpwind_r[31]+alpha[25]*fUpwind_r[29])+0.2*alpha[16]*fUpwind_r[28]+0.223606797749979*(alpha[19]*fUpwind_r[26]+fUpwind_r[19]*alpha[26])+(0.159719141249985*alpha[19]+0.2500000000000001*alpha[2])*fUpwind_r[25]+0.159719141249985*fUpwind_r[19]*alpha[25]+0.2500000000000001*(fUpwind_r[2]*alpha[25]+fUpwind_r[18]*alpha[21])+0.2*(alpha[16]*fUpwind_r[20]+fUpwind_r[16]*alpha[20])+0.2500000000000001*(alpha[4]*fUpwind_r[19]+fUpwind_r[4]*alpha[19])+0.223606797749979*(alpha[15]*fUpwind_r[17]+alpha[1]*fUpwind_r[16]+fUpwind_r[1]*alpha[16])+0.25*(alpha[9]*fUpwind_r[11]+fUpwind_r[9]*alpha[11])+0.223606797749979*(alpha[5]*fUpwind_r[8]+fUpwind_r[5]*alpha[8]); 
  Ghat_r[36] += (0.159719141249985*alpha[22]+0.223606797749979*alpha[21]+0.2500000000000001*alpha[3])*fUpwind_r[45]+0.2*alpha[15]*fUpwind_r[44]+(0.1788854381999831*alpha[35]+0.2*alpha[9])*fUpwind_r[41]+(0.159719141249985*alpha[33]+0.25*alpha[6])*fUpwind_r[38]+0.223606797749979*alpha[33]*fUpwind_r[37]+(0.159719141249985*alpha[12]+0.223606797749979*alpha[11]+0.25*alpha[0])*fUpwind_r[36]+(0.223606797749979*fUpwind_r[14]+0.159719141249985*fUpwind_r[12]+0.223606797749979*fUpwind_r[11]+0.25*fUpwind_r[0])*alpha[36]+0.2*(alpha[5]*fUpwind_r[35]+fUpwind_r[5]*alpha[35])+0.25*fUpwind_r[10]*alpha[33]+fUpwind_r[31]*(0.2*alpha[32]+0.223606797749979*alpha[7])+0.2*alpha[16]*fUpwind_r[29]+0.223606797749979*alpha[26]*fUpwind_r[28]+(0.159719141249985*alpha[20]+0.2500000000000001*alpha[1])*fUpwind_r[26]+(0.159719141249985*fUpwind_r[20]+0.2500000000000001*fUpwind_r[1])*alpha[26]+0.223606797749979*(alpha[20]*fUpwind_r[25]+fUpwind_r[20]*alpha[25])+0.2500000000000001*(fUpwind_r[17]*alpha[22]+alpha[4]*fUpwind_r[20]+fUpwind_r[4]*alpha[20])+0.2*(alpha[16]*fUpwind_r[19]+fUpwind_r[16]*alpha[19])+0.223606797749979*(alpha[15]*fUpwind_r[18]+alpha[2]*fUpwind_r[16]+fUpwind_r[2]*alpha[16])+0.25*(alpha[8]*fUpwind_r[12]+fUpwind_r[8]*alpha[12])+0.223606797749979*(alpha[5]*fUpwind_r[9]+fUpwind_r[5]*alpha[9]); 
  Ghat_r[37] += 0.2*(alpha[16]*fUpwind_r[47]+alpha[15]*fUpwind_r[46])+0.223606797749979*alpha[20]*fUpwind_r[45]+(0.159719141249985*alpha[19]+0.2500000000000001*alpha[2])*fUpwind_r[44]+0.223606797749979*alpha[35]*fUpwind_r[43]+0.2*alpha[8]*fUpwind_r[42]+0.223606797749979*alpha[32]*fUpwind_r[40]+0.2*alpha[6]*fUpwind_r[39]+(0.159719141249985*alpha[11]+0.25*alpha[0])*fUpwind_r[37]+0.223606797749979*(alpha[33]*fUpwind_r[36]+fUpwind_r[33]*alpha[36])+(0.159719141249985*alpha[32]+0.25*alpha[7])*fUpwind_r[35]+0.159719141249985*fUpwind_r[32]*alpha[35]+0.25*(fUpwind_r[7]*alpha[35]+alpha[9]*fUpwind_r[32]+fUpwind_r[9]*alpha[32])+0.223606797749979*(alpha[5]*fUpwind_r[31]+alpha[25]*fUpwind_r[30]+alpha[21]*fUpwind_r[27])+(0.159719141249985*alpha[21]+0.2500000000000001*alpha[3])*fUpwind_r[25]+0.159719141249985*fUpwind_r[21]*alpha[25]+0.2500000000000001*(fUpwind_r[3]*alpha[25]+alpha[4]*fUpwind_r[21]+fUpwind_r[4]*alpha[21]+fUpwind_r[18]*alpha[19])+0.223606797749979*(alpha[1]*fUpwind_r[17]+alpha[15]*fUpwind_r[16]+fUpwind_r[15]*alpha[16])+0.25*fUpwind_r[10]*alpha[11]+0.223606797749979*(alpha[6]*fUpwind_r[8]+fUpwind_r[6]*alpha[8]); 
  Ghat_r[38] += 0.2*(alpha[16]*fUpwind_r[47]+alpha[15]*fUpwind_r[46])+(0.159719141249985*alpha[20]+0.2500000000000001*alpha[1])*fUpwind_r[45]+0.223606797749979*alpha[19]*fUpwind_r[44]+0.2*alpha[9]*fUpwind_r[43]+0.223606797749979*alpha[36]*fUpwind_r[42]+0.2*alpha[7]*fUpwind_r[40]+0.223606797749979*alpha[33]*fUpwind_r[39]+(0.159719141249985*alpha[12]+0.25*alpha[0])*fUpwind_r[38]+(0.159719141249985*alpha[33]+0.25*alpha[6])*fUpwind_r[36]+(0.159719141249985*fUpwind_r[33]+0.25*fUpwind_r[6])*alpha[36]+0.223606797749979*(alpha[32]*fUpwind_r[35]+fUpwind_r[32]*alpha[35])+0.25*(alpha[8]*fUpwind_r[33]+fUpwind_r[8]*alpha[33])+0.223606797749979*(alpha[5]*fUpwind_r[31]+alpha[26]*fUpwind_r[30]+alpha[22]*fUpwind_r[27])+(0.159719141249985*alpha[22]+0.2500000000000001*alpha[3])*fUpwind_r[26]+0.159719141249985*fUpwind_r[22]*alpha[26]+0.2500000000000001*(fUpwind_r[3]*alpha[26]+alpha[4]*fUpwind_r[22]+fUpwind_r[4]*alpha[22]+fUpwind_r[17]*alpha[20])+0.223606797749979*(alpha[2]*fUpwind_r[18]+alpha[15]*fUpwind_r[16]+fUpwind_r[15]*alpha[16])+0.25*fUpwind_r[10]*alpha[12]+0.223606797749979*(alpha[7]*fUpwind_r[9]+fUpwind_r[7]*alpha[9]); 
  Ghat_r[39] += (0.223606797749979*alpha[19]+0.2500000000000001*alpha[2])*fUpwind_r[46]+0.223606797749979*alpha[22]*fUpwind_r[45]+0.2*alpha[15]*fUpwind_r[44]+0.25*alpha[5]*fUpwind_r[40]+(0.223606797749979*alpha[11]+0.25*alpha[0])*fUpwind_r[39]+0.223606797749979*alpha[33]*fUpwind_r[38]+0.2*alpha[6]*fUpwind_r[37]+fUpwind_r[34]*(0.223606797749979*alpha[35]+0.25*alpha[9])+fUpwind_r[31]*(0.2*alpha[32]+0.223606797749979*alpha[7])+0.2500000000000001*alpha[1]*fUpwind_r[27]+0.223606797749979*fUpwind_r[23]*alpha[25]+0.2500000000000001*(alpha[16]*fUpwind_r[24]+alpha[4]*fUpwind_r[23])+0.2*fUpwind_r[17]*alpha[21]+0.223606797749979*(alpha[15]*fUpwind_r[18]+alpha[3]*fUpwind_r[17])+0.25*alpha[8]*fUpwind_r[13]+0.223606797749979*alpha[6]*fUpwind_r[10]; 
  Ghat_r[40] += (0.223606797749979*alpha[20]+0.2500000000000001*alpha[1])*fUpwind_r[46]+0.2*alpha[15]*fUpwind_r[45]+0.223606797749979*(alpha[21]*fUpwind_r[44]+alpha[12]*fUpwind_r[40])+0.25*(alpha[0]*fUpwind_r[40]+alpha[5]*fUpwind_r[39])+0.2*alpha[7]*fUpwind_r[38]+0.223606797749979*alpha[32]*fUpwind_r[37]+fUpwind_r[34]*(0.223606797749979*alpha[36]+0.25*alpha[8])+fUpwind_r[31]*(0.2*alpha[33]+0.223606797749979*alpha[6])+0.2500000000000001*alpha[2]*fUpwind_r[27]+0.223606797749979*fUpwind_r[24]*alpha[26]+0.2500000000000001*(alpha[4]*fUpwind_r[24]+alpha[16]*fUpwind_r[23])+0.2*fUpwind_r[18]*alpha[22]+0.223606797749979*(alpha[3]*fUpwind_r[18]+alpha[15]*fUpwind_r[17])+0.25*alpha[9]*fUpwind_r[13]+0.223606797749979*alpha[7]*fUpwind_r[10]; 
  Ghat_r[41] += (0.223606797749979*(alpha[22]+alpha[21])+0.2500000000000001*alpha[3])*fUpwind_r[47]+(0.223606797749979*alpha[33]+0.25*alpha[6])*fUpwind_r[43]+(0.223606797749979*alpha[32]+0.25*alpha[7])*fUpwind_r[42]+(0.223606797749979*(alpha[12]+alpha[11])+0.25*alpha[0])*fUpwind_r[41]+(0.1788854381999831*alpha[35]+0.2*alpha[9])*fUpwind_r[36]+0.1788854381999831*fUpwind_r[35]*alpha[36]+0.2*(fUpwind_r[9]*alpha[36]+alpha[8]*fUpwind_r[35]+fUpwind_r[8]*alpha[35])+0.2500000000000001*alpha[15]*fUpwind_r[30]+(0.223606797749979*alpha[20]+0.2500000000000001*alpha[1])*fUpwind_r[29]+(0.223606797749979*alpha[19]+0.2500000000000001*alpha[2])*fUpwind_r[28]+0.2*(alpha[16]*fUpwind_r[26]+fUpwind_r[16]*alpha[26]+alpha[16]*fUpwind_r[25]+fUpwind_r[16]*alpha[25])+0.223606797749979*(alpha[4]*fUpwind_r[16]+fUpwind_r[4]*alpha[16])+0.25*alpha[5]*fUpwind_r[14]+0.223606797749979*(alpha[8]*fUpwind_r[9]+fUpwind_r[8]*alpha[9]); 
  Ghat_r[42] += (0.223606797749979*alpha[19]+0.2500000000000001*alpha[2])*fUpwind_r[47]+0.223606797749979*alpha[26]*fUpwind_r[45]+0.2*alpha[16]*fUpwind_r[44]+0.25*alpha[5]*fUpwind_r[43]+(0.223606797749979*alpha[11]+0.25*alpha[0])*fUpwind_r[42]+(0.223606797749979*alpha[32]+0.25*alpha[7])*fUpwind_r[41]+0.223606797749979*alpha[36]*fUpwind_r[38]+0.2*alpha[8]*fUpwind_r[37]+fUpwind_r[31]*(0.2*alpha[35]+0.223606797749979*alpha[9])+0.2500000000000001*(alpha[1]*fUpwind_r[30]+alpha[15]*fUpwind_r[29])+(0.223606797749979*alpha[21]+0.2500000000000001*alpha[3])*fUpwind_r[28]+0.2*fUpwind_r[17]*alpha[25]+0.223606797749979*(alpha[16]*fUpwind_r[18]+alpha[4]*fUpwind_r[17])+0.25*alpha[6]*fUpwind_r[14]+0.223606797749979*alpha[8]*fUpwind_r[10]; 
  Ghat_r[43] += (0.223606797749979*alpha[20]+0.2500000000000001*alpha[1])*fUpwind_r[47]+0.2*alpha[16]*fUpwind_r[45]+0.223606797749979*(alpha[25]*fUpwind_r[44]+alpha[12]*fUpwind_r[43])+0.25*(alpha[0]*fUpwind_r[43]+alpha[5]*fUpwind_r[42])+(0.223606797749979*alpha[33]+0.25*alpha[6])*fUpwind_r[41]+0.2*alpha[9]*fUpwind_r[38]+0.223606797749979*alpha[35]*fUpwind_r[37]+fUpwind_r[31]*(0.2*alpha[36]+0.223606797749979*alpha[8])+0.2500000000000001*alpha[2]*fUpwind_r[30]+0.223606797749979*alpha[22]*fUpwind_r[29]+0.2500000000000001*(alpha[3]*fUpwind_r[29]+alpha[15]*fUpwind_r[28])+0.2*fUpwind_r[18]*alpha[26]+0.223606797749979*(alpha[4]*fUpwind_r[18]+alpha[16]*fUpwind_r[17])+0.25*alpha[7]*fUpwind_r[14]+0.223606797749979*alpha[9]*fUpwind_r[10]; 
  Ghat_r[44] += (0.1788854381999831*alpha[36]+0.2*alpha[8])*fUpwind_r[47]+0.1788854381999831*alpha[33]*fUpwind_r[46]+0.2*(alpha[6]*fUpwind_r[46]+alpha[5]*fUpwind_r[45])+(0.223606797749979*alpha[12]+0.159719141249985*alpha[11]+0.25*alpha[0])*fUpwind_r[44]+0.223606797749979*alpha[25]*fUpwind_r[43]+0.2*alpha[16]*fUpwind_r[42]+0.223606797749979*alpha[21]*fUpwind_r[40]+0.2*alpha[15]*fUpwind_r[39]+0.223606797749979*alpha[19]*fUpwind_r[38]+(0.159719141249985*alpha[19]+0.2500000000000001*alpha[2])*fUpwind_r[37]+0.2*(alpha[15]*fUpwind_r[36]+fUpwind_r[15]*alpha[36])+(0.223606797749979*alpha[22]+0.159719141249985*alpha[21]+0.2500000000000001*alpha[3])*fUpwind_r[35]+(0.223606797749979*(fUpwind_r[30]+fUpwind_r[22])+0.159719141249985*fUpwind_r[21]+0.2500000000000001*fUpwind_r[3])*alpha[35]+0.2*(alpha[16]*fUpwind_r[33]+fUpwind_r[16]*alpha[33])+(0.223606797749979*alpha[26]+0.159719141249985*alpha[25]+0.2500000000000001*alpha[4])*fUpwind_r[32]+(0.223606797749979*(fUpwind_r[27]+fUpwind_r[26])+0.159719141249985*fUpwind_r[25]+0.2500000000000001*fUpwind_r[4])*alpha[32]+(0.2*alpha[20]+0.223606797749979*alpha[1])*fUpwind_r[31]+0.25*(alpha[7]*fUpwind_r[25]+fUpwind_r[7]*alpha[25]+alpha[9]*fUpwind_r[21]+fUpwind_r[9]*alpha[21]+fUpwind_r[10]*alpha[19])+0.2500000000000001*alpha[11]*fUpwind_r[18]+0.223606797749979*(alpha[5]*fUpwind_r[17]+alpha[6]*fUpwind_r[16]+fUpwind_r[6]*alpha[16]+alpha[8]*fUpwind_r[15]+fUpwind_r[8]*alpha[15]); 
  Ghat_r[45] += (0.1788854381999831*alpha[35]+0.2*alpha[9])*fUpwind_r[47]+(0.1788854381999831*alpha[32]+0.2*alpha[7])*fUpwind_r[46]+(0.159719141249985*alpha[12]+0.223606797749979*alpha[11]+0.25*alpha[0])*fUpwind_r[45]+0.2*(alpha[5]*fUpwind_r[44]+alpha[16]*fUpwind_r[43])+0.223606797749979*alpha[26]*fUpwind_r[42]+0.2*alpha[15]*fUpwind_r[40]+0.223606797749979*alpha[22]*fUpwind_r[39]+(0.159719141249985*alpha[20]+0.2500000000000001*alpha[1])*fUpwind_r[38]+0.223606797749979*alpha[20]*fUpwind_r[37]+(0.159719141249985*alpha[22]+0.223606797749979*alpha[21]+0.2500000000000001*alpha[3])*fUpwind_r[36]+(0.223606797749979*fUpwind_r[30]+0.159719141249985*fUpwind_r[22]+0.223606797749979*fUpwind_r[21]+0.2500000000000001*fUpwind_r[3])*alpha[36]+0.2*(alpha[15]*fUpwind_r[35]+fUpwind_r[15]*alpha[35])+(0.159719141249985*alpha[26]+0.223606797749979*alpha[25]+0.2500000000000001*alpha[4])*fUpwind_r[33]+(0.223606797749979*fUpwind_r[27]+0.159719141249985*fUpwind_r[26]+0.223606797749979*fUpwind_r[25]+0.2500000000000001*fUpwind_r[4])*alpha[33]+0.2*(alpha[16]*fUpwind_r[32]+fUpwind_r[16]*alpha[32])+(0.2*alpha[19]+0.223606797749979*alpha[2])*fUpwind_r[31]+0.25*(alpha[6]*fUpwind_r[26]+fUpwind_r[6]*alpha[26]+alpha[8]*fUpwind_r[22]+fUpwind_r[8]*alpha[22]+fUpwind_r[10]*alpha[20])+0.223606797749979*alpha[5]*fUpwind_r[18]+0.2500000000000001*alpha[12]*fUpwind_r[17]+0.223606797749979*(alpha[7]*fUpwind_r[16]+fUpwind_r[7]*alpha[16]+alpha[9]*fUpwind_r[15]+fUpwind_r[9]*alpha[15]); 
  Ghat_r[46] += (0.223606797749979*(alpha[12]+alpha[11])+0.25*alpha[0])*fUpwind_r[46]+(0.1788854381999831*alpha[32]+0.2*alpha[7])*fUpwind_r[45]+(0.1788854381999831*alpha[33]+0.2*alpha[6])*fUpwind_r[44]+(0.223606797749979*alpha[20]+0.2500000000000001*alpha[1])*fUpwind_r[40]+(0.223606797749979*alpha[19]+0.2500000000000001*alpha[2])*fUpwind_r[39]+0.2*alpha[15]*(fUpwind_r[38]+fUpwind_r[37])+0.223606797749979*(fUpwind_r[24]*alpha[36]+fUpwind_r[23]*alpha[35])+(0.223606797749979*(alpha[26]+alpha[25])+0.2500000000000001*alpha[4])*fUpwind_r[34]+0.2*(fUpwind_r[18]*alpha[33]+fUpwind_r[17]*alpha[32])+(0.2*(alpha[22]+alpha[21])+0.223606797749979*alpha[3])*fUpwind_r[31]+0.25*(alpha[5]*fUpwind_r[27]+alpha[8]*fUpwind_r[24]+alpha[9]*fUpwind_r[23])+0.223606797749979*(alpha[6]*fUpwind_r[18]+alpha[7]*fUpwind_r[17])+0.2500000000000001*fUpwind_r[13]*alpha[16]+0.223606797749979*fUpwind_r[10]*alpha[15]; 
  Ghat_r[47] += (0.223606797749979*(alpha[12]+alpha[11])+0.25*alpha[0])*fUpwind_r[47]+(0.1788854381999831*alpha[35]+0.2*alpha[9])*fUpwind_r[45]+(0.1788854381999831*alpha[36]+0.2*alpha[8])*fUpwind_r[44]+(0.223606797749979*alpha[20]+0.2500000000000001*alpha[1])*fUpwind_r[43]+(0.223606797749979*alpha[19]+0.2500000000000001*alpha[2])*fUpwind_r[42]+(0.223606797749979*(alpha[22]+alpha[21])+0.2500000000000001*alpha[3])*fUpwind_r[41]+0.2*(alpha[16]*(fUpwind_r[38]+fUpwind_r[37])+fUpwind_r[18]*alpha[36]+fUpwind_r[17]*alpha[35])+0.223606797749979*(fUpwind_r[29]*alpha[33]+fUpwind_r[28]*alpha[32])+(0.2*(alpha[26]+alpha[25])+0.223606797749979*alpha[4])*fUpwind_r[31]+0.25*(alpha[5]*fUpwind_r[30]+alpha[6]*fUpwind_r[29]+alpha[7]*fUpwind_r[28])+0.223606797749979*(alpha[8]*fUpwind_r[18]+alpha[9]*fUpwind_r[17]+fUpwind_r[10]*alpha[16])+0.2500000000000001*fUpwind_r[14]*alpha[15]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv10; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv10; 
  out[2] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv10; 
  out[3] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv10; 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv10; 
  out[5] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dv10; 
  out[6] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dv10; 
  out[7] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv10; 
  out[8] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv10; 
  out[9] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dv10; 
  out[10] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dv10; 
  out[11] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv10; 
  out[12] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*dv10; 
  out[13] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*dv10; 
  out[14] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dv10; 
  out[15] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*dv10; 
  out[16] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*dv10; 
  out[17] += (0.7071067811865475*Ghat_l[12]-0.7071067811865475*Ghat_r[12])*dv10; 
  out[18] += (1.58113883008419*Ghat_l[0]-1.58113883008419*Ghat_r[0])*dv10; 
  out[19] += (0.7071067811865475*Ghat_l[13]-0.7071067811865475*Ghat_r[13])*dv10; 
  out[20] += (0.7071067811865475*Ghat_l[14]-0.7071067811865475*Ghat_r[14])*dv10; 
  out[21] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dv10; 
  out[22] += (0.7071067811865475*Ghat_l[15]-0.7071067811865475*Ghat_r[15])*dv10; 
  out[23] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dv10; 
  out[24] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dv10; 
  out[25] += (0.7071067811865475*Ghat_l[16]-0.7071067811865475*Ghat_r[16])*dv10; 
  out[26] += -1.224744871391589*(Ghat_r[8]+Ghat_l[8])*dv10; 
  out[27] += -1.224744871391589*(Ghat_r[9]+Ghat_l[9])*dv10; 
  out[28] += (0.7071067811865475*Ghat_l[17]-0.7071067811865475*Ghat_r[17])*dv10; 
  out[29] += (0.7071067811865475*Ghat_l[18]-0.7071067811865475*Ghat_r[18])*dv10; 
  out[30] += -1.224744871391589*(Ghat_r[10]+Ghat_l[10])*dv10; 
  out[31] += (0.7071067811865475*Ghat_l[19]-0.7071067811865475*Ghat_r[19])*dv10; 
  out[32] += (0.7071067811865475*Ghat_l[20]-0.7071067811865475*Ghat_r[20])*dv10; 
  out[33] += -1.224744871391589*(Ghat_r[11]+Ghat_l[11])*dv10; 
  out[34] += -1.224744871391589*(Ghat_r[12]+Ghat_l[12])*dv10; 
  out[35] += (1.58113883008419*Ghat_l[1]-1.58113883008419*Ghat_r[1])*dv10; 
  out[36] += (1.58113883008419*Ghat_l[2]-1.58113883008419*Ghat_r[2])*dv10; 
  out[37] += (0.7071067811865475*Ghat_l[21]-0.7071067811865475*Ghat_r[21])*dv10; 
  out[38] += (0.7071067811865475*Ghat_l[22]-0.7071067811865475*Ghat_r[22])*dv10; 
  out[39] += (1.58113883008419*Ghat_l[3]-1.58113883008419*Ghat_r[3])*dv10; 
  out[40] += (0.7071067811865475*Ghat_l[23]-0.7071067811865475*Ghat_r[23])*dv10; 
  out[41] += (0.7071067811865475*Ghat_l[24]-0.7071067811865475*Ghat_r[24])*dv10; 
  out[42] += -1.224744871391589*(Ghat_r[13]+Ghat_l[13])*dv10; 
  out[43] += (0.7071067811865475*Ghat_l[25]-0.7071067811865475*Ghat_r[25])*dv10; 
  out[44] += (0.7071067811865475*Ghat_l[26]-0.7071067811865475*Ghat_r[26])*dv10; 
  out[45] += (1.58113883008419*Ghat_l[4]-1.58113883008419*Ghat_r[4])*dv10; 
  out[46] += (0.7071067811865475*Ghat_l[27]-0.7071067811865475*Ghat_r[27])*dv10; 
  out[47] += (0.7071067811865475*Ghat_l[28]-0.7071067811865475*Ghat_r[28])*dv10; 
  out[48] += (0.7071067811865475*Ghat_l[29]-0.7071067811865475*Ghat_r[29])*dv10; 
  out[49] += -1.224744871391589*(Ghat_r[14]+Ghat_l[14])*dv10; 
  out[50] += (0.7071067811865475*Ghat_l[30]-0.7071067811865475*Ghat_r[30])*dv10; 
  out[51] += -1.224744871391589*(Ghat_r[15]+Ghat_l[15])*dv10; 
  out[52] += -1.224744871391589*(Ghat_r[16]+Ghat_l[16])*dv10; 
  out[53] += (0.7071067811865475*Ghat_l[31]-0.7071067811865475*Ghat_r[31])*dv10; 
  out[54] += -1.224744871391589*(Ghat_r[17]+Ghat_l[17])*dv10; 
  out[55] += -1.224744871391589*(Ghat_r[18]+Ghat_l[18])*dv10; 
  out[56] += -1.224744871391589*(Ghat_r[19]+Ghat_l[19])*dv10; 
  out[57] += -1.224744871391589*(Ghat_r[20]+Ghat_l[20])*dv10; 
  out[58] += (1.58113883008419*Ghat_l[5]-1.58113883008419*Ghat_r[5])*dv10; 
  out[59] += (0.7071067811865475*Ghat_l[32]-0.7071067811865475*Ghat_r[32])*dv10; 
  out[60] += (0.7071067811865475*Ghat_l[33]-0.7071067811865475*Ghat_r[33])*dv10; 
  out[61] += -1.224744871391589*(Ghat_r[21]+Ghat_l[21])*dv10; 
  out[62] += -1.224744871391589*(Ghat_r[22]+Ghat_l[22])*dv10; 
  out[63] += (1.58113883008419*Ghat_l[6]-1.58113883008419*Ghat_r[6])*dv10; 
  out[64] += (1.58113883008419*Ghat_l[7]-1.58113883008419*Ghat_r[7])*dv10; 
  out[65] += (0.7071067811865475*Ghat_l[34]-0.7071067811865475*Ghat_r[34])*dv10; 
  out[66] += -1.224744871391589*(Ghat_r[23]+Ghat_l[23])*dv10; 
  out[67] += -1.224744871391589*(Ghat_r[24]+Ghat_l[24])*dv10; 
  out[68] += (0.7071067811865475*Ghat_l[35]-0.7071067811865475*Ghat_r[35])*dv10; 
  out[69] += (0.7071067811865475*Ghat_l[36]-0.7071067811865475*Ghat_r[36])*dv10; 
  out[70] += -1.224744871391589*(Ghat_r[25]+Ghat_l[25])*dv10; 
  out[71] += -1.224744871391589*(Ghat_r[26]+Ghat_l[26])*dv10; 
  out[72] += (1.58113883008419*Ghat_l[8]-1.58113883008419*Ghat_r[8])*dv10; 
  out[73] += (1.58113883008419*Ghat_l[9]-1.58113883008419*Ghat_r[9])*dv10; 
  out[74] += (0.7071067811865475*Ghat_l[37]-0.7071067811865475*Ghat_r[37])*dv10; 
  out[75] += (0.7071067811865475*Ghat_l[38]-0.7071067811865475*Ghat_r[38])*dv10; 
  out[76] += (1.58113883008419*Ghat_l[10]-1.58113883008419*Ghat_r[10])*dv10; 
  out[77] += (0.7071067811865475*Ghat_l[39]-0.7071067811865475*Ghat_r[39])*dv10; 
  out[78] += (0.7071067811865475*Ghat_l[40]-0.7071067811865475*Ghat_r[40])*dv10; 
  out[79] += -1.224744871391589*(Ghat_r[27]+Ghat_l[27])*dv10; 
  out[80] += (0.7071067811865475*Ghat_l[41]-0.7071067811865475*Ghat_r[41])*dv10; 
  out[81] += -1.224744871391589*(Ghat_r[28]+Ghat_l[28])*dv10; 
  out[82] += -1.224744871391589*(Ghat_r[29]+Ghat_l[29])*dv10; 
  out[83] += (0.7071067811865475*Ghat_l[42]-0.7071067811865475*Ghat_r[42])*dv10; 
  out[84] += (0.7071067811865475*Ghat_l[43]-0.7071067811865475*Ghat_r[43])*dv10; 
  out[85] += -1.224744871391589*(Ghat_r[30]+Ghat_l[30])*dv10; 
  out[86] += -1.224744871391589*(Ghat_r[31]+Ghat_l[31])*dv10; 
  out[87] += -1.224744871391589*(Ghat_r[32]+Ghat_l[32])*dv10; 
  out[88] += -1.224744871391589*(Ghat_r[33]+Ghat_l[33])*dv10; 
  out[89] += (1.58113883008419*Ghat_l[15]-1.58113883008419*Ghat_r[15])*dv10; 
  out[90] += -1.224744871391589*(Ghat_r[34]+Ghat_l[34])*dv10; 
  out[91] += -1.224744871391589*(Ghat_r[35]+Ghat_l[35])*dv10; 
  out[92] += -1.224744871391589*(Ghat_r[36]+Ghat_l[36])*dv10; 
  out[93] += (1.58113883008419*Ghat_l[16]-1.58113883008419*Ghat_r[16])*dv10; 
  out[94] += (0.7071067811865475*Ghat_l[44]-0.7071067811865475*Ghat_r[44])*dv10; 
  out[95] += (0.7071067811865475*Ghat_l[45]-0.7071067811865475*Ghat_r[45])*dv10; 
  out[96] += -1.224744871391589*(Ghat_r[37]+Ghat_l[37])*dv10; 
  out[97] += -1.224744871391589*(Ghat_r[38]+Ghat_l[38])*dv10; 
  out[98] += (1.58113883008419*Ghat_l[17]-1.58113883008419*Ghat_r[17])*dv10; 
  out[99] += (1.58113883008419*Ghat_l[18]-1.58113883008419*Ghat_r[18])*dv10; 
  out[100] += (0.7071067811865475*Ghat_l[46]-0.7071067811865475*Ghat_r[46])*dv10; 
  out[101] += -1.224744871391589*(Ghat_r[39]+Ghat_l[39])*dv10; 
  out[102] += -1.224744871391589*(Ghat_r[40]+Ghat_l[40])*dv10; 
  out[103] += -1.224744871391589*(Ghat_r[41]+Ghat_l[41])*dv10; 
  out[104] += (0.7071067811865475*Ghat_l[47]-0.7071067811865475*Ghat_r[47])*dv10; 
  out[105] += -1.224744871391589*(Ghat_r[42]+Ghat_l[42])*dv10; 
  out[106] += -1.224744871391589*(Ghat_r[43]+Ghat_l[43])*dv10; 
  out[107] += -1.224744871391589*(Ghat_r[44]+Ghat_l[44])*dv10; 
  out[108] += -1.224744871391589*(Ghat_r[45]+Ghat_l[45])*dv10; 
  out[109] += (1.58113883008419*Ghat_l[31]-1.58113883008419*Ghat_r[31])*dv10; 
  out[110] += -1.224744871391589*(Ghat_r[46]+Ghat_l[46])*dv10; 
  out[111] += -1.224744871391589*(Ghat_r[47]+Ghat_l[47])*dv10; 

} 
