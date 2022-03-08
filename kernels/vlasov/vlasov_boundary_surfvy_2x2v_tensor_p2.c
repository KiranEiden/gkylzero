#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_tensor_4x_p2_surfx4_quad.h> 
#include <gkyl_basis_tensor_4x_p2_upwind.h> 
GKYL_CU_DH void vlasov_boundary_surfvy_2x2v_tensor_p2(const double *w, const double *dxv, const double *qmem, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // qmem:        q/m*EM fields.
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv11 = 2/dxv[3]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double *E1 = &qmem[9]; 
  const double *B2 = &qmem[45]; 

  double alpha[27] = {0.0}; 

  alpha[0] = 1.414213562373095*E1[0]-1.414213562373095*B2[0]*wv1; 
  alpha[1] = 1.414213562373095*E1[1]-1.414213562373095*B2[1]*wv1; 
  alpha[2] = 1.414213562373095*E1[2]-1.414213562373095*B2[2]*wv1; 
  alpha[3] = -0.408248290463863*B2[0]*dv1; 
  alpha[4] = 1.414213562373095*E1[3]-1.414213562373095*B2[3]*wv1; 
  alpha[5] = -0.408248290463863*B2[1]*dv1; 
  alpha[6] = -0.408248290463863*B2[2]*dv1; 
  alpha[7] = 1.414213562373095*E1[4]-1.414213562373095*B2[4]*wv1; 
  alpha[8] = 1.414213562373095*E1[5]-1.414213562373095*B2[5]*wv1; 
  alpha[10] = -0.408248290463863*B2[3]*dv1; 
  alpha[11] = 1.414213562373095*E1[6]-1.414213562373095*B2[6]*wv1; 
  alpha[12] = 1.414213562373095*E1[7]-1.414213562373095*B2[7]*wv1; 
  alpha[13] = -0.408248290463863*B2[4]*dv1; 
  alpha[14] = -0.408248290463863*B2[5]*dv1; 
  alpha[17] = -0.408248290463863*B2[6]*dv1; 
  alpha[18] = -0.408248290463863*B2[7]*dv1; 
  alpha[20] = 1.414213562373095*E1[8]-1.414213562373095*B2[8]*wv1; 
  alpha[23] = -0.408248290463863*B2[8]*dv1; 

  double fUpwindQuad[27] = {0.0};
  double fUpwind[27] = {0.0};
  double Ghat[27] = {0.0}; 

  if (edge == -1) { 

  if ((-0.3794733192202055*alpha[23])+0.2828427124746191*alpha[20]+0.5692099788303082*(alpha[18]+alpha[17])-0.4242640687119285*(alpha[14]+alpha[13])-0.4242640687119281*alpha[12]-0.4242640687119285*alpha[11]-0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*(alpha[6]+alpha[5]+alpha[4])-0.4743416490252568*(alpha[3]+alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[0] = tensor_4x_p2_surfx4_quad_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = tensor_4x_p2_surfx4_quad_0_l(fEdge); 
  } 
  if (0.4743416490252568*alpha[23]-0.3535533905932737*alpha[20]-0.711512473537885*alpha[17]-0.4242640687119285*alpha[14]+0.5303300858899104*(alpha[13]+alpha[11])+0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]+0.6363961030678926*alpha[6]-0.4743416490252568*(alpha[3]+alpha[2])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[1] = tensor_4x_p2_surfx4_quad_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = tensor_4x_p2_surfx4_quad_1_l(fEdge); 
  } 
  if ((-0.3794733192202055*alpha[23])+0.2828427124746191*alpha[20]-0.5692099788303082*alpha[18]+0.5692099788303082*alpha[17]-0.4242640687119285*(alpha[14]+alpha[13])+0.4242640687119281*alpha[12]-0.4242640687119285*alpha[11]+0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*alpha[6]-0.6363961030678926*(alpha[5]+alpha[4])-0.4743416490252568*(alpha[3]+alpha[2])+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[2] = tensor_4x_p2_surfx4_quad_2_r(fSkin); 
  } else { 
    fUpwindQuad[2] = tensor_4x_p2_surfx4_quad_2_l(fEdge); 
  } 
  if (0.4743416490252568*alpha[23]-0.3535533905932737*alpha[20]-0.711512473537885*alpha[18]+0.5303300858899104*alpha[14]-0.4242640687119285*alpha[13]+0.5303300858899104*alpha[12]-0.3952847075210473*alpha[8]+0.3162277660168379*alpha[7]+0.6363961030678926*alpha[5]-0.4743416490252568*(alpha[3]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[3] = tensor_4x_p2_surfx4_quad_3_r(fSkin); 
  } else { 
    fUpwindQuad[3] = tensor_4x_p2_surfx4_quad_3_l(fEdge); 
  } 
  if ((-0.592927061281571*alpha[23])+0.441941738241592*alpha[20]+0.5303300858899104*(alpha[14]+alpha[13])-0.3952847075210473*(alpha[8]+alpha[7])-0.4743416490252568*alpha[3]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[4] = tensor_4x_p2_surfx4_quad_4_r(fSkin); 
  } else { 
    fUpwindQuad[4] = tensor_4x_p2_surfx4_quad_4_l(fEdge); 
  } 
  if (0.4743416490252568*alpha[23]-0.3535533905932737*alpha[20]+0.711512473537885*alpha[18]+0.5303300858899104*alpha[14]-0.4242640687119285*alpha[13]-0.5303300858899104*alpha[12]-0.3952847075210473*alpha[8]+0.3162277660168379*alpha[7]-0.6363961030678926*alpha[5]-0.4743416490252568*alpha[3]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[5] = tensor_4x_p2_surfx4_quad_5_r(fSkin); 
  } else { 
    fUpwindQuad[5] = tensor_4x_p2_surfx4_quad_5_l(fEdge); 
  } 
  if ((-0.3794733192202055*alpha[23])+0.2828427124746191*alpha[20]+0.5692099788303082*alpha[18]-0.5692099788303082*alpha[17]-0.4242640687119285*(alpha[14]+alpha[13])-0.4242640687119281*alpha[12]+0.4242640687119285*alpha[11]+0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*alpha[6]+0.6363961030678926*alpha[5]-0.6363961030678926*alpha[4]-0.4743416490252568*alpha[3]+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[6] = tensor_4x_p2_surfx4_quad_6_r(fSkin); 
  } else { 
    fUpwindQuad[6] = tensor_4x_p2_surfx4_quad_6_l(fEdge); 
  } 
  if (0.4743416490252568*alpha[23]-0.3535533905932737*alpha[20]+0.711512473537885*alpha[17]-0.4242640687119285*alpha[14]+0.5303300858899104*alpha[13]-0.5303300858899104*alpha[11]+0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]-0.6363961030678926*alpha[6]-0.4743416490252568*alpha[3]+0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[7] = tensor_4x_p2_surfx4_quad_7_r(fSkin); 
  } else { 
    fUpwindQuad[7] = tensor_4x_p2_surfx4_quad_7_l(fEdge); 
  } 
  if ((-0.3794733192202055*alpha[23])+0.2828427124746191*alpha[20]-0.5692099788303082*(alpha[18]+alpha[17])-0.4242640687119285*(alpha[14]+alpha[13])+0.4242640687119281*alpha[12]+0.4242640687119285*alpha[11]-0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*(alpha[6]+alpha[5])+0.6363961030678926*alpha[4]-0.4743416490252568*alpha[3]+0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[8] = tensor_4x_p2_surfx4_quad_8_r(fSkin); 
  } else { 
    fUpwindQuad[8] = tensor_4x_p2_surfx4_quad_8_l(fEdge); 
  } 
  if (0.2828427124746191*alpha[20]-0.4242640687119281*alpha[12]-0.4242640687119285*alpha[11]+0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*alpha[4]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[9] = tensor_4x_p2_surfx4_quad_9_r(fSkin); 
  } else { 
    fUpwindQuad[9] = tensor_4x_p2_surfx4_quad_9_l(fEdge); 
  } 
  if ((-0.3535533905932737*alpha[20])+0.5303300858899104*alpha[11]+0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]-0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[10] = tensor_4x_p2_surfx4_quad_10_r(fSkin); 
  } else { 
    fUpwindQuad[10] = tensor_4x_p2_surfx4_quad_10_l(fEdge); 
  } 
  if (0.2828427124746191*alpha[20]+0.4242640687119281*alpha[12]-0.4242640687119285*alpha[11]+0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*alpha[4]-0.4743416490252568*alpha[2]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[11] = tensor_4x_p2_surfx4_quad_11_r(fSkin); 
  } else { 
    fUpwindQuad[11] = tensor_4x_p2_surfx4_quad_11_l(fEdge); 
  } 
  if ((-0.3535533905932737*alpha[20])+0.5303300858899104*alpha[12]-0.3952847075210473*alpha[8]+0.3162277660168379*alpha[7]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[12] = tensor_4x_p2_surfx4_quad_12_r(fSkin); 
  } else { 
    fUpwindQuad[12] = tensor_4x_p2_surfx4_quad_12_l(fEdge); 
  } 
  if (0.441941738241592*alpha[20]-0.3952847075210473*(alpha[8]+alpha[7])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[13] = tensor_4x_p2_surfx4_quad_13_r(fSkin); 
  } else { 
    fUpwindQuad[13] = tensor_4x_p2_surfx4_quad_13_l(fEdge); 
  } 
  if ((-0.3535533905932737*alpha[20])-0.5303300858899104*alpha[12]-0.3952847075210473*alpha[8]+0.3162277660168379*alpha[7]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[14] = tensor_4x_p2_surfx4_quad_14_r(fSkin); 
  } else { 
    fUpwindQuad[14] = tensor_4x_p2_surfx4_quad_14_l(fEdge); 
  } 
  if (0.2828427124746191*alpha[20]-0.4242640687119281*alpha[12]+0.4242640687119285*alpha[11]+0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*alpha[4]+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[15] = tensor_4x_p2_surfx4_quad_15_r(fSkin); 
  } else { 
    fUpwindQuad[15] = tensor_4x_p2_surfx4_quad_15_l(fEdge); 
  } 
  if ((-0.3535533905932737*alpha[20])-0.5303300858899104*alpha[11]+0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]+0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[16] = tensor_4x_p2_surfx4_quad_16_r(fSkin); 
  } else { 
    fUpwindQuad[16] = tensor_4x_p2_surfx4_quad_16_l(fEdge); 
  } 
  if (0.2828427124746191*alpha[20]+0.4242640687119281*alpha[12]+0.4242640687119285*alpha[11]+0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*alpha[4]+0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[17] = tensor_4x_p2_surfx4_quad_17_r(fSkin); 
  } else { 
    fUpwindQuad[17] = tensor_4x_p2_surfx4_quad_17_l(fEdge); 
  } 
  if (0.3794733192202055*alpha[23]+0.2828427124746191*alpha[20]-0.5692099788303082*(alpha[18]+alpha[17])+0.4242640687119285*(alpha[14]+alpha[13])-0.4242640687119281*alpha[12]-0.4242640687119285*alpha[11]+0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*(alpha[6]+alpha[5])+0.6363961030678926*alpha[4]+0.4743416490252568*alpha[3]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[18] = tensor_4x_p2_surfx4_quad_18_r(fSkin); 
  } else { 
    fUpwindQuad[18] = tensor_4x_p2_surfx4_quad_18_l(fEdge); 
  } 
  if ((-0.4743416490252568*alpha[23])-0.3535533905932737*alpha[20]+0.711512473537885*alpha[17]+0.4242640687119285*alpha[14]-0.5303300858899104*alpha[13]+0.5303300858899104*alpha[11]+0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]-0.6363961030678926*alpha[6]+0.4743416490252568*alpha[3]-0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[19] = tensor_4x_p2_surfx4_quad_19_r(fSkin); 
  } else { 
    fUpwindQuad[19] = tensor_4x_p2_surfx4_quad_19_l(fEdge); 
  } 
  if (0.3794733192202055*alpha[23]+0.2828427124746191*alpha[20]+0.5692099788303082*alpha[18]-0.5692099788303082*alpha[17]+0.4242640687119285*(alpha[14]+alpha[13])+0.4242640687119281*alpha[12]-0.4242640687119285*alpha[11]-0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*alpha[6]+0.6363961030678926*alpha[5]-0.6363961030678926*alpha[4]+0.4743416490252568*alpha[3]-0.4743416490252568*alpha[2]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[20] = tensor_4x_p2_surfx4_quad_20_r(fSkin); 
  } else { 
    fUpwindQuad[20] = tensor_4x_p2_surfx4_quad_20_l(fEdge); 
  } 
  if ((-0.4743416490252568*alpha[23])-0.3535533905932737*alpha[20]+0.711512473537885*alpha[18]-0.5303300858899104*alpha[14]+0.4242640687119285*alpha[13]+0.5303300858899104*alpha[12]-0.3952847075210473*alpha[8]+0.3162277660168379*alpha[7]-0.6363961030678926*alpha[5]+0.4743416490252568*alpha[3]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[21] = tensor_4x_p2_surfx4_quad_21_r(fSkin); 
  } else { 
    fUpwindQuad[21] = tensor_4x_p2_surfx4_quad_21_l(fEdge); 
  } 
  if (0.592927061281571*alpha[23]+0.441941738241592*alpha[20]-0.5303300858899104*(alpha[14]+alpha[13])-0.3952847075210473*(alpha[8]+alpha[7])+0.4743416490252568*alpha[3]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[22] = tensor_4x_p2_surfx4_quad_22_r(fSkin); 
  } else { 
    fUpwindQuad[22] = tensor_4x_p2_surfx4_quad_22_l(fEdge); 
  } 
  if ((-0.4743416490252568*alpha[23])-0.3535533905932737*alpha[20]-0.711512473537885*alpha[18]-0.5303300858899104*alpha[14]+0.4242640687119285*alpha[13]-0.5303300858899104*alpha[12]-0.3952847075210473*alpha[8]+0.3162277660168379*alpha[7]+0.6363961030678926*alpha[5]+0.4743416490252568*(alpha[3]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[23] = tensor_4x_p2_surfx4_quad_23_r(fSkin); 
  } else { 
    fUpwindQuad[23] = tensor_4x_p2_surfx4_quad_23_l(fEdge); 
  } 
  if (0.3794733192202055*alpha[23]+0.2828427124746191*alpha[20]-0.5692099788303082*alpha[18]+0.5692099788303082*alpha[17]+0.4242640687119285*(alpha[14]+alpha[13])-0.4242640687119281*alpha[12]+0.4242640687119285*alpha[11]-0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*alpha[6]-0.6363961030678926*(alpha[5]+alpha[4])+0.4743416490252568*(alpha[3]+alpha[2])-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[24] = tensor_4x_p2_surfx4_quad_24_r(fSkin); 
  } else { 
    fUpwindQuad[24] = tensor_4x_p2_surfx4_quad_24_l(fEdge); 
  } 
  if ((-0.4743416490252568*alpha[23])-0.3535533905932737*alpha[20]-0.711512473537885*alpha[17]+0.4242640687119285*alpha[14]-0.5303300858899104*(alpha[13]+alpha[11])+0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]+0.6363961030678926*alpha[6]+0.4743416490252568*(alpha[3]+alpha[2])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[25] = tensor_4x_p2_surfx4_quad_25_r(fSkin); 
  } else { 
    fUpwindQuad[25] = tensor_4x_p2_surfx4_quad_25_l(fEdge); 
  } 
  if (0.3794733192202055*alpha[23]+0.2828427124746191*alpha[20]+0.5692099788303082*(alpha[18]+alpha[17])+0.4242640687119285*(alpha[14]+alpha[13])+0.4242640687119281*alpha[12]+0.4242640687119285*alpha[11]+0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*(alpha[6]+alpha[5]+alpha[4])+0.4743416490252568*(alpha[3]+alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[26] = tensor_4x_p2_surfx4_quad_26_r(fSkin); 
  } else { 
    fUpwindQuad[26] = tensor_4x_p2_surfx4_quad_26_l(fEdge); 
  } 

  // Project nodal basis back onto modal basis. 
  tensor_4x_p2_upwind(fUpwindQuad, fUpwind); 

  Ghat[0] += 0.3535533905932737*(alpha[23]*fUpwind[23]+alpha[20]*fUpwind[20]+alpha[18]*fUpwind[18]+alpha[17]*fUpwind[17]+alpha[14]*fUpwind[14]+alpha[13]*fUpwind[13]+alpha[12]*fUpwind[12]+alpha[11]*fUpwind[11]+alpha[10]*fUpwind[10]+alpha[8]*fUpwind[8]+alpha[7]*fUpwind[7]+alpha[6]*fUpwind[6]+alpha[5]*fUpwind[5]+alpha[4]*fUpwind[4]+alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] += 0.3162277660168379*(alpha[18]*fUpwind[23]+fUpwind[18]*alpha[23])+0.3162277660168379*(alpha[12]*fUpwind[20]+fUpwind[12]*alpha[20])+0.3535533905932737*(alpha[14]*fUpwind[18]+fUpwind[14]*alpha[18])+0.3162277660168379*(alpha[10]*fUpwind[17]+fUpwind[10]*alpha[17])+0.3162277660168379*(alpha[5]*fUpwind[13]+fUpwind[5]*alpha[13])+0.3535533905932737*(alpha[8]*fUpwind[12]+fUpwind[8]*alpha[12])+0.3162277660168379*(alpha[4]*fUpwind[11]+fUpwind[4]*alpha[11])+0.3535533905932737*(alpha[6]*fUpwind[10]+fUpwind[6]*alpha[10])+0.3162277660168379*(alpha[1]*fUpwind[7]+fUpwind[1]*alpha[7])+0.3535533905932737*(alpha[3]*fUpwind[5]+fUpwind[3]*alpha[5]+alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] += 0.3162277660168379*(alpha[17]*fUpwind[23]+fUpwind[17]*alpha[23])+0.3162277660168379*(alpha[11]*fUpwind[20]+fUpwind[11]*alpha[20])+0.3162277660168379*(alpha[10]*fUpwind[18]+fUpwind[10]*alpha[18])+0.3535533905932737*(alpha[13]*fUpwind[17]+fUpwind[13]*alpha[17])+0.3162277660168379*(alpha[6]*fUpwind[14]+fUpwind[6]*alpha[14]+alpha[4]*fUpwind[12]+fUpwind[4]*alpha[12])+0.3535533905932737*(alpha[7]*fUpwind[11]+fUpwind[7]*alpha[11]+alpha[5]*fUpwind[10]+fUpwind[5]*alpha[10])+0.3162277660168379*(alpha[2]*fUpwind[8]+fUpwind[2]*alpha[8])+0.3535533905932737*(alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6]+alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] += 0.3162277660168379*(alpha[23]*fUpwind[26]+alpha[18]*fUpwind[25]+alpha[17]*fUpwind[24])+0.3535533905932737*(alpha[20]*fUpwind[23]+fUpwind[20]*alpha[23])+0.3162277660168379*(alpha[14]*fUpwind[22]+alpha[13]*fUpwind[21])+0.3162277660168379*alpha[10]*fUpwind[19]+0.3535533905932737*(alpha[12]*fUpwind[18]+fUpwind[12]*alpha[18]+alpha[11]*fUpwind[17]+fUpwind[11]*alpha[17])+0.3162277660168379*(alpha[6]*fUpwind[16]+alpha[5]*fUpwind[15])+0.3535533905932737*(alpha[8]*fUpwind[14]+fUpwind[8]*alpha[14]+alpha[7]*fUpwind[13]+fUpwind[7]*alpha[13]+alpha[4]*fUpwind[10]+fUpwind[4]*alpha[10])+0.3162277660168379*alpha[3]*fUpwind[9]+0.3535533905932737*(alpha[2]*fUpwind[6]+fUpwind[2]*alpha[6]+alpha[1]*fUpwind[5]+fUpwind[1]*alpha[5]+alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]); 
  Ghat[4] += 0.2828427124746191*(alpha[10]*fUpwind[23]+fUpwind[10]*alpha[23]+alpha[4]*fUpwind[20]+fUpwind[4]*alpha[20])+(0.2828427124746191*alpha[17]+0.3162277660168379*alpha[6])*fUpwind[18]+0.2828427124746191*fUpwind[17]*alpha[18]+0.3162277660168379*(fUpwind[6]*alpha[18]+alpha[5]*fUpwind[17]+fUpwind[5]*alpha[17])+0.3162277660168379*(alpha[10]*fUpwind[14]+fUpwind[10]*alpha[14]+alpha[10]*fUpwind[13]+fUpwind[10]*alpha[13])+(0.2828427124746191*alpha[11]+0.3162277660168379*alpha[2])*fUpwind[12]+0.2828427124746191*fUpwind[11]*alpha[12]+0.3162277660168379*(fUpwind[2]*alpha[12]+alpha[1]*fUpwind[11]+fUpwind[1]*alpha[11])+0.3535533905932737*(alpha[3]*fUpwind[10]+fUpwind[3]*alpha[10])+0.3162277660168379*(alpha[4]*fUpwind[8]+fUpwind[4]*alpha[8]+alpha[4]*fUpwind[7]+fUpwind[4]*alpha[7])+0.3535533905932737*(alpha[5]*fUpwind[6]+fUpwind[5]*alpha[6]+alpha[0]*fUpwind[4]+fUpwind[0]*alpha[4]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[5] += 0.2828427124746191*alpha[18]*fUpwind[26]+(0.2828427124746191*alpha[23]+0.3162277660168379*alpha[14])*fUpwind[25]+0.2828427124746191*alpha[10]*fUpwind[24]+0.3162277660168379*(alpha[12]*fUpwind[23]+fUpwind[12]*alpha[23])+0.3162277660168379*alpha[18]*fUpwind[22]+0.2828427124746191*alpha[5]*fUpwind[21]+0.3162277660168379*(alpha[18]*fUpwind[20]+fUpwind[18]*alpha[20])+(0.2828427124746191*alpha[17]+0.3162277660168379*alpha[6])*fUpwind[19]+0.3535533905932737*(alpha[8]*fUpwind[18]+fUpwind[8]*alpha[18])+0.3162277660168379*(alpha[4]*fUpwind[17]+fUpwind[4]*alpha[17])+0.3162277660168379*alpha[10]*fUpwind[16]+(0.2828427124746191*alpha[13]+0.3162277660168379*alpha[3])*fUpwind[15]+0.3535533905932737*(alpha[12]*fUpwind[14]+fUpwind[12]*alpha[14])+0.3162277660168379*(alpha[1]*fUpwind[13]+fUpwind[1]*alpha[13]+alpha[10]*fUpwind[11]+fUpwind[10]*alpha[11])+0.3535533905932737*(alpha[2]*fUpwind[10]+fUpwind[2]*alpha[10])+0.3162277660168379*(alpha[5]*(fUpwind[9]+fUpwind[7])+fUpwind[5]*alpha[7])+0.3535533905932737*(alpha[4]*fUpwind[6]+fUpwind[4]*alpha[6]+alpha[0]*fUpwind[5]+fUpwind[0]*alpha[5]+alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]); 
  Ghat[6] += 0.2828427124746191*(alpha[17]*fUpwind[26]+alpha[10]*fUpwind[25]+alpha[23]*fUpwind[24])+0.3162277660168379*(alpha[13]*fUpwind[24]+alpha[11]*fUpwind[23]+fUpwind[11]*alpha[23])+0.2828427124746191*alpha[6]*fUpwind[22]+0.3162277660168379*(alpha[17]*(fUpwind[21]+fUpwind[20])+fUpwind[17]*alpha[20])+0.2828427124746191*alpha[18]*fUpwind[19]+0.3162277660168379*(alpha[5]*fUpwind[19]+alpha[4]*fUpwind[18]+fUpwind[4]*alpha[18])+0.3535533905932737*(alpha[7]*fUpwind[17]+fUpwind[7]*alpha[17])+0.2828427124746191*alpha[14]*fUpwind[16]+0.3162277660168379*(alpha[3]*fUpwind[16]+alpha[10]*fUpwind[15]+alpha[2]*fUpwind[14]+fUpwind[2]*alpha[14])+0.3535533905932737*(alpha[11]*fUpwind[13]+fUpwind[11]*alpha[13])+0.3162277660168379*(alpha[10]*fUpwind[12]+fUpwind[10]*alpha[12])+0.3535533905932737*(alpha[1]*fUpwind[10]+fUpwind[1]*alpha[10])+0.3162277660168379*(alpha[6]*(fUpwind[9]+fUpwind[8])+fUpwind[6]*alpha[8])+0.3535533905932737*(alpha[0]*fUpwind[6]+fUpwind[0]*alpha[6]+alpha[4]*fUpwind[5]+fUpwind[4]*alpha[5]+alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]); 
  Ghat[7] += 0.2258769757263128*alpha[23]*fUpwind[23]+0.3535533905932737*(alpha[14]*fUpwind[23]+fUpwind[14]*alpha[23])+0.2258769757263128*alpha[20]*fUpwind[20]+0.3535533905932737*(alpha[8]*fUpwind[20]+fUpwind[8]*alpha[20])+0.3162277660168379*alpha[18]*fUpwind[18]+0.2258769757263128*alpha[17]*fUpwind[17]+0.3535533905932737*(alpha[6]*fUpwind[17]+fUpwind[6]*alpha[17])+0.2258769757263128*alpha[13]*fUpwind[13]+0.3535533905932737*(alpha[3]*fUpwind[13]+fUpwind[3]*alpha[13])+0.3162277660168379*alpha[12]*fUpwind[12]+0.2258769757263128*alpha[11]*fUpwind[11]+0.3535533905932737*(alpha[2]*fUpwind[11]+fUpwind[2]*alpha[11])+0.3162277660168379*alpha[10]*fUpwind[10]+0.2258769757263128*alpha[7]*fUpwind[7]+0.3535533905932737*(alpha[0]*fUpwind[7]+fUpwind[0]*alpha[7])+0.3162277660168379*(alpha[5]*fUpwind[5]+alpha[4]*fUpwind[4]+alpha[1]*fUpwind[1]); 
  Ghat[8] += 0.2258769757263128*alpha[23]*fUpwind[23]+0.3535533905932737*(alpha[13]*fUpwind[23]+fUpwind[13]*alpha[23])+0.2258769757263128*alpha[20]*fUpwind[20]+0.3535533905932737*(alpha[7]*fUpwind[20]+fUpwind[7]*alpha[20])+0.2258769757263128*alpha[18]*fUpwind[18]+0.3535533905932737*(alpha[5]*fUpwind[18]+fUpwind[5]*alpha[18])+0.3162277660168379*alpha[17]*fUpwind[17]+0.2258769757263128*alpha[14]*fUpwind[14]+0.3535533905932737*(alpha[3]*fUpwind[14]+fUpwind[3]*alpha[14])+0.2258769757263128*alpha[12]*fUpwind[12]+0.3535533905932737*(alpha[1]*fUpwind[12]+fUpwind[1]*alpha[12])+0.3162277660168379*(alpha[11]*fUpwind[11]+alpha[10]*fUpwind[10])+0.2258769757263128*alpha[8]*fUpwind[8]+0.3535533905932737*(alpha[0]*fUpwind[8]+fUpwind[0]*alpha[8])+0.3162277660168379*(alpha[6]*fUpwind[6]+alpha[4]*fUpwind[4]+alpha[2]*fUpwind[2]); 
  Ghat[9] += 0.3535533905932737*(alpha[20]*fUpwind[26]+alpha[12]*fUpwind[25]+alpha[11]*fUpwind[24])+0.3162277660168379*alpha[23]*fUpwind[23]+0.3535533905932737*(alpha[8]*fUpwind[22]+alpha[7]*fUpwind[21]+alpha[4]*fUpwind[19])+0.3162277660168379*(alpha[18]*fUpwind[18]+alpha[17]*fUpwind[17])+0.3535533905932737*(alpha[2]*fUpwind[16]+alpha[1]*fUpwind[15])+0.3162277660168379*(alpha[14]*fUpwind[14]+alpha[13]*fUpwind[13]+alpha[10]*fUpwind[10])+0.3535533905932737*alpha[0]*fUpwind[9]+0.3162277660168379*(alpha[6]*fUpwind[6]+alpha[5]*fUpwind[5]+alpha[3]*fUpwind[3]); 
  Ghat[10] += 0.2529822128134704*alpha[10]*fUpwind[26]+(0.2529822128134704*alpha[17]+0.2828427124746191*alpha[6])*fUpwind[25]+0.2529822128134704*alpha[18]*fUpwind[24]+0.2828427124746191*(alpha[5]*fUpwind[24]+alpha[4]*fUpwind[23])+0.2529822128134704*fUpwind[19]*alpha[23]+0.2828427124746191*(fUpwind[4]*alpha[23]+alpha[10]*(fUpwind[22]+fUpwind[21]+fUpwind[20])+fUpwind[10]*alpha[20])+(0.282842712474619*(alpha[14]+alpha[13])+0.3162277660168379*alpha[3])*fUpwind[19]+(0.282842712474619*alpha[11]+0.3162277660168379*alpha[2])*fUpwind[18]+(0.282842712474619*(fUpwind[16]+fUpwind[11])+0.3162277660168379*fUpwind[2])*alpha[18]+(0.282842712474619*alpha[12]+0.3162277660168379*alpha[1])*fUpwind[17]+(0.282842712474619*(fUpwind[15]+fUpwind[12])+0.3162277660168379*fUpwind[1])*alpha[17]+0.3162277660168379*(alpha[5]*fUpwind[16]+alpha[6]*fUpwind[15]+alpha[4]*fUpwind[14]+fUpwind[4]*alpha[14]+alpha[4]*fUpwind[13]+fUpwind[4]*alpha[13]+alpha[6]*fUpwind[12]+fUpwind[6]*alpha[12]+alpha[5]*fUpwind[11]+fUpwind[5]*alpha[11])+(0.3162277660168379*(alpha[8]+alpha[7])+0.3535533905932737*alpha[0])*fUpwind[10]+0.3162277660168379*(fUpwind[9]+fUpwind[8]+fUpwind[7])*alpha[10]+0.3535533905932737*(fUpwind[0]*alpha[10]+alpha[1]*fUpwind[6]+fUpwind[1]*alpha[6]+alpha[2]*fUpwind[5]+fUpwind[2]*alpha[5]+alpha[3]*fUpwind[4]+fUpwind[3]*alpha[4]); 
  Ghat[11] += (0.2020305089104422*alpha[17]+0.3162277660168379*alpha[6])*fUpwind[23]+(0.2020305089104422*fUpwind[17]+0.3162277660168379*fUpwind[6])*alpha[23]+(0.2020305089104422*alpha[11]+0.3162277660168379*alpha[2])*fUpwind[20]+(0.2020305089104422*fUpwind[11]+0.3162277660168379*fUpwind[2])*alpha[20]+0.282842712474619*(alpha[10]*fUpwind[18]+fUpwind[10]*alpha[18])+(0.3162277660168379*alpha[14]+0.2258769757263128*alpha[13]+0.3535533905932737*alpha[3])*fUpwind[17]+(0.3162277660168379*fUpwind[14]+0.2258769757263128*fUpwind[13])*alpha[17]+0.3535533905932737*(fUpwind[3]*alpha[17]+alpha[6]*fUpwind[13]+fUpwind[6]*alpha[13])+0.2828427124746191*(alpha[4]*fUpwind[12]+fUpwind[4]*alpha[12])+(0.3162277660168379*alpha[8]+0.2258769757263128*alpha[7]+0.3535533905932737*alpha[0])*fUpwind[11]+(0.3162277660168379*fUpwind[8]+0.2258769757263128*fUpwind[7]+0.3535533905932737*fUpwind[0])*alpha[11]+0.3162277660168379*(alpha[5]*fUpwind[10]+fUpwind[5]*alpha[10])+0.3535533905932737*(alpha[2]*fUpwind[7]+fUpwind[2]*alpha[7])+0.3162277660168379*(alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]); 
  Ghat[12] += (0.2020305089104422*alpha[18]+0.3162277660168379*alpha[5])*fUpwind[23]+(0.2020305089104422*fUpwind[18]+0.3162277660168379*fUpwind[5])*alpha[23]+(0.2020305089104422*alpha[12]+0.3162277660168379*alpha[1])*fUpwind[20]+(0.2020305089104422*fUpwind[12]+0.3162277660168379*fUpwind[1])*alpha[20]+(0.2258769757263128*alpha[14]+0.3162277660168379*alpha[13]+0.3535533905932737*alpha[3])*fUpwind[18]+(0.2258769757263128*fUpwind[14]+0.3162277660168379*fUpwind[13]+0.3535533905932737*fUpwind[3])*alpha[18]+0.282842712474619*(alpha[10]*fUpwind[17]+fUpwind[10]*alpha[17])+0.3535533905932737*(alpha[5]*fUpwind[14]+fUpwind[5]*alpha[14])+(0.2258769757263128*alpha[8]+0.3162277660168379*alpha[7]+0.3535533905932737*alpha[0])*fUpwind[12]+(0.2258769757263128*fUpwind[8]+0.3162277660168379*fUpwind[7]+0.3535533905932737*fUpwind[0])*alpha[12]+0.2828427124746191*(alpha[4]*fUpwind[11]+fUpwind[4]*alpha[11])+0.3162277660168379*(alpha[6]*fUpwind[10]+fUpwind[6]*alpha[10])+0.3535533905932737*(alpha[1]*fUpwind[8]+fUpwind[1]*alpha[8])+0.3162277660168379*(alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]); 
  Ghat[13] += (0.2020305089104422*alpha[23]+0.3162277660168379*alpha[14])*fUpwind[26]+0.282842712474619*alpha[18]*fUpwind[25]+(0.2020305089104422*alpha[17]+0.3162277660168379*alpha[6])*fUpwind[24]+(0.2258769757263128*alpha[20]+0.3535533905932737*alpha[8])*fUpwind[23]+(0.3162277660168379*fUpwind[22]+0.2258769757263128*fUpwind[20]+0.3535533905932737*fUpwind[8])*alpha[23]+(0.2020305089104422*alpha[13]+0.3162277660168379*alpha[3])*fUpwind[21]+0.3535533905932737*(alpha[14]*fUpwind[20]+fUpwind[14]*alpha[20])+0.282842712474619*alpha[10]*fUpwind[19]+0.3162277660168379*(alpha[12]*fUpwind[18]+fUpwind[12]*alpha[18])+(0.2258769757263128*alpha[11]+0.3535533905932737*alpha[2])*fUpwind[17]+(0.3162277660168379*fUpwind[16]+0.2258769757263128*fUpwind[11]+0.3535533905932737*fUpwind[2])*alpha[17]+0.2828427124746191*alpha[5]*fUpwind[15]+(0.2258769757263128*alpha[7]+0.3535533905932737*alpha[0])*fUpwind[13]+(0.3162277660168379*fUpwind[9]+0.2258769757263128*fUpwind[7])*alpha[13]+0.3535533905932737*(fUpwind[0]*alpha[13]+alpha[6]*fUpwind[11]+fUpwind[6]*alpha[11])+0.3162277660168379*(alpha[4]*fUpwind[10]+fUpwind[4]*alpha[10])+0.3535533905932737*(alpha[3]*fUpwind[7]+fUpwind[3]*alpha[7])+0.3162277660168379*(alpha[1]*fUpwind[5]+fUpwind[1]*alpha[5]); 
  Ghat[14] += (0.2020305089104422*alpha[23]+0.3162277660168379*alpha[13])*fUpwind[26]+(0.2020305089104422*alpha[18]+0.3162277660168379*alpha[5])*fUpwind[25]+0.282842712474619*alpha[17]*fUpwind[24]+(0.2258769757263128*alpha[20]+0.3535533905932737*alpha[7])*fUpwind[23]+(0.3162277660168379*fUpwind[21]+0.2258769757263128*fUpwind[20]+0.3535533905932737*fUpwind[7])*alpha[23]+(0.2020305089104422*alpha[14]+0.3162277660168379*alpha[3])*fUpwind[22]+0.3535533905932737*(alpha[13]*fUpwind[20]+fUpwind[13]*alpha[20])+0.282842712474619*alpha[10]*fUpwind[19]+(0.2258769757263128*alpha[12]+0.3535533905932737*alpha[1])*fUpwind[18]+(0.3162277660168379*fUpwind[15]+0.2258769757263128*fUpwind[12]+0.3535533905932737*fUpwind[1])*alpha[18]+0.3162277660168379*(alpha[11]*fUpwind[17]+fUpwind[11]*alpha[17])+0.2828427124746191*alpha[6]*fUpwind[16]+(0.2258769757263128*alpha[8]+0.3535533905932737*alpha[0])*fUpwind[14]+(0.3162277660168379*fUpwind[9]+0.2258769757263128*fUpwind[8])*alpha[14]+0.3535533905932737*(fUpwind[0]*alpha[14]+alpha[5]*fUpwind[12]+fUpwind[5]*alpha[12])+0.3162277660168379*(alpha[4]*fUpwind[10]+fUpwind[4]*alpha[10])+0.3535533905932737*(alpha[3]*fUpwind[8]+fUpwind[3]*alpha[8])+0.3162277660168379*(alpha[2]*fUpwind[6]+fUpwind[2]*alpha[6]); 
  Ghat[15] += 0.3162277660168379*alpha[12]*fUpwind[26]+(0.3162277660168379*alpha[20]+0.3535533905932737*alpha[8])*fUpwind[25]+0.3162277660168379*alpha[4]*fUpwind[24]+0.282842712474619*(alpha[18]*fUpwind[23]+fUpwind[18]*alpha[23])+0.3535533905932737*alpha[12]*fUpwind[22]+0.3162277660168379*alpha[1]*fUpwind[21]+(0.3162277660168379*alpha[11]+0.3535533905932737*alpha[2])*fUpwind[19]+0.3162277660168379*(alpha[14]*fUpwind[18]+fUpwind[14]*alpha[18])+0.282842712474619*(alpha[10]*fUpwind[17]+fUpwind[10]*alpha[17])+0.3535533905932737*alpha[4]*fUpwind[16]+(0.3162277660168379*alpha[7]+0.3535533905932737*alpha[0])*fUpwind[15]+0.2828427124746191*(alpha[5]*fUpwind[13]+fUpwind[5]*alpha[13])+0.3162277660168379*(alpha[6]*fUpwind[10]+fUpwind[6]*alpha[10])+0.3535533905932737*alpha[1]*fUpwind[9]+0.3162277660168379*(alpha[3]*fUpwind[5]+fUpwind[3]*alpha[5]); 
  Ghat[16] += 0.3162277660168379*alpha[11]*fUpwind[26]+0.3162277660168379*alpha[4]*fUpwind[25]+(0.3162277660168379*alpha[20]+0.3535533905932737*alpha[7])*fUpwind[24]+0.282842712474619*(alpha[17]*fUpwind[23]+fUpwind[17]*alpha[23])+0.3162277660168379*alpha[2]*fUpwind[22]+0.3535533905932737*alpha[11]*fUpwind[21]+(0.3162277660168379*alpha[12]+0.3535533905932737*alpha[1])*fUpwind[19]+0.282842712474619*(alpha[10]*fUpwind[18]+fUpwind[10]*alpha[18])+0.3162277660168379*(alpha[13]*fUpwind[17]+fUpwind[13]*alpha[17]+alpha[8]*fUpwind[16])+0.3535533905932737*(alpha[0]*fUpwind[16]+alpha[4]*fUpwind[15])+0.2828427124746191*(alpha[6]*fUpwind[14]+fUpwind[6]*alpha[14])+0.3162277660168379*(alpha[5]*fUpwind[10]+fUpwind[5]*alpha[10])+0.3535533905932737*alpha[2]*fUpwind[9]+0.3162277660168379*(alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6]); 
  Ghat[17] += (0.1807015805810503*alpha[17]+0.2828427124746191*alpha[6])*fUpwind[26]+0.2529822128134704*alpha[10]*fUpwind[25]+(0.1807015805810503*alpha[23]+0.282842712474619*alpha[14]+0.2020305089104422*alpha[13]+0.3162277660168379*alpha[3])*fUpwind[24]+(0.2020305089104422*alpha[11]+0.3162277660168379*alpha[2])*fUpwind[23]+(0.282842712474619*fUpwind[16]+0.2020305089104422*fUpwind[11]+0.3162277660168379*fUpwind[2])*alpha[23]+0.2828427124746191*alpha[17]*fUpwind[22]+(0.2020305089104422*alpha[17]+0.3162277660168379*alpha[6])*fUpwind[21]+(0.2020305089104422*alpha[17]+0.3162277660168379*alpha[6])*fUpwind[20]+(0.2020305089104422*fUpwind[17]+0.3162277660168379*fUpwind[6])*alpha[20]+0.2529822128134704*alpha[18]*fUpwind[19]+0.2828427124746191*(alpha[5]*fUpwind[19]+alpha[4]*fUpwind[18]+fUpwind[4]*alpha[18])+(0.3162277660168379*alpha[8]+0.2258769757263128*alpha[7]+0.3535533905932737*alpha[0])*fUpwind[17]+(0.3162277660168379*(fUpwind[9]+fUpwind[8])+0.2258769757263128*fUpwind[7]+0.3535533905932737*fUpwind[0])*alpha[17]+0.3162277660168379*alpha[13]*fUpwind[16]+0.282842712474619*alpha[10]*fUpwind[15]+0.3162277660168379*(alpha[11]*fUpwind[14]+fUpwind[11]*alpha[14])+(0.2258769757263128*alpha[11]+0.3535533905932737*alpha[2])*fUpwind[13]+(0.2258769757263128*fUpwind[11]+0.3535533905932737*fUpwind[2])*alpha[13]+0.282842712474619*(alpha[10]*fUpwind[12]+fUpwind[10]*alpha[12])+0.3535533905932737*(alpha[3]*fUpwind[11]+fUpwind[3]*alpha[11])+0.3162277660168379*(alpha[1]*fUpwind[10]+fUpwind[1]*alpha[10])+0.3535533905932737*(alpha[6]*fUpwind[7]+fUpwind[6]*alpha[7])+0.3162277660168379*(alpha[4]*fUpwind[5]+fUpwind[4]*alpha[5]); 
  Ghat[18] += (0.1807015805810503*alpha[18]+0.2828427124746191*alpha[5])*fUpwind[26]+(0.1807015805810503*alpha[23]+0.2020305089104422*alpha[14]+0.282842712474619*alpha[13]+0.3162277660168379*alpha[3])*fUpwind[25]+0.2529822128134704*alpha[10]*fUpwind[24]+(0.2020305089104422*alpha[12]+0.3162277660168379*alpha[1])*fUpwind[23]+(0.282842712474619*fUpwind[15]+0.2020305089104422*fUpwind[12]+0.3162277660168379*fUpwind[1])*alpha[23]+(0.2020305089104422*alpha[18]+0.3162277660168379*alpha[5])*fUpwind[22]+0.2828427124746191*alpha[18]*fUpwind[21]+(0.2020305089104422*alpha[18]+0.3162277660168379*alpha[5])*fUpwind[20]+(0.2020305089104422*fUpwind[18]+0.3162277660168379*fUpwind[5])*alpha[20]+(0.2529822128134704*alpha[17]+0.2828427124746191*alpha[6])*fUpwind[19]+(0.2258769757263128*alpha[8]+0.3162277660168379*alpha[7]+0.3535533905932737*alpha[0])*fUpwind[18]+(0.3162277660168379*fUpwind[9]+0.2258769757263128*fUpwind[8]+0.3162277660168379*fUpwind[7]+0.3535533905932737*fUpwind[0])*alpha[18]+0.2828427124746191*(alpha[4]*fUpwind[17]+fUpwind[4]*alpha[17])+0.282842712474619*alpha[10]*fUpwind[16]+0.3162277660168379*alpha[14]*fUpwind[15]+(0.2258769757263128*alpha[12]+0.3535533905932737*alpha[1])*fUpwind[14]+(0.2258769757263128*fUpwind[12]+0.3535533905932737*fUpwind[1])*alpha[14]+0.3162277660168379*(alpha[12]*fUpwind[13]+fUpwind[12]*alpha[13])+0.3535533905932737*(alpha[3]*fUpwind[12]+fUpwind[3]*alpha[12])+0.282842712474619*(alpha[10]*fUpwind[11]+fUpwind[10]*alpha[11])+0.3162277660168379*(alpha[2]*fUpwind[10]+fUpwind[2]*alpha[10])+0.3535533905932737*(alpha[5]*fUpwind[8]+fUpwind[5]*alpha[8])+0.3162277660168379*(alpha[4]*fUpwind[6]+fUpwind[4]*alpha[6]); 
  Ghat[19] += 0.2828427124746191*alpha[4]*fUpwind[26]+(0.282842712474619*alpha[11]+0.3162277660168379*alpha[2])*fUpwind[25]+(0.282842712474619*alpha[12]+0.3162277660168379*alpha[1])*fUpwind[24]+0.2529822128134704*(alpha[10]*fUpwind[23]+fUpwind[10]*alpha[23])+0.3162277660168379*alpha[4]*(fUpwind[22]+fUpwind[21])+fUpwind[19]*(0.2828427124746191*alpha[20]+0.3162277660168379*(alpha[8]+alpha[7])+0.3535533905932737*alpha[0])+(0.2529822128134704*alpha[17]+0.2828427124746191*alpha[6])*fUpwind[18]+0.2529822128134704*fUpwind[17]*alpha[18]+0.2828427124746191*(fUpwind[6]*alpha[18]+alpha[5]*fUpwind[17]+fUpwind[5]*alpha[17])+(0.3162277660168379*alpha[12]+0.3535533905932737*alpha[1])*fUpwind[16]+(0.3162277660168379*alpha[11]+0.3535533905932737*alpha[2])*fUpwind[15]+0.282842712474619*(alpha[10]*fUpwind[14]+fUpwind[10]*alpha[14]+alpha[10]*fUpwind[13]+fUpwind[10]*alpha[13])+0.3162277660168379*(alpha[3]*fUpwind[10]+fUpwind[3]*alpha[10])+0.3535533905932737*alpha[4]*fUpwind[9]+0.3162277660168379*(alpha[5]*fUpwind[6]+fUpwind[5]*alpha[6]); 
  Ghat[20] += (0.1443075063646015*alpha[23]+0.2258769757263128*(alpha[14]+alpha[13])+0.3535533905932737*alpha[3])*fUpwind[23]+(0.2258769757263128*(fUpwind[14]+fUpwind[13])+0.3535533905932737*fUpwind[3])*alpha[23]+(0.1443075063646015*alpha[20]+0.2258769757263128*(alpha[8]+alpha[7])+0.3535533905932737*alpha[0])*fUpwind[20]+(0.2258769757263128*(fUpwind[8]+fUpwind[7])+0.3535533905932737*fUpwind[0])*alpha[20]+0.2020305089104422*alpha[18]*fUpwind[18]+0.3162277660168379*(alpha[5]*fUpwind[18]+fUpwind[5]*alpha[18])+0.2020305089104422*alpha[17]*fUpwind[17]+0.3162277660168379*(alpha[6]*fUpwind[17]+fUpwind[6]*alpha[17])+0.3535533905932737*(alpha[13]*fUpwind[14]+fUpwind[13]*alpha[14])+0.2020305089104422*alpha[12]*fUpwind[12]+0.3162277660168379*(alpha[1]*fUpwind[12]+fUpwind[1]*alpha[12])+0.2020305089104422*alpha[11]*fUpwind[11]+0.3162277660168379*(alpha[2]*fUpwind[11]+fUpwind[2]*alpha[11])+0.2828427124746191*alpha[10]*fUpwind[10]+0.3535533905932737*(alpha[7]*fUpwind[8]+fUpwind[7]*alpha[8])+0.2828427124746191*alpha[4]*fUpwind[4]; 
  Ghat[21] += (0.2258769757263128*alpha[20]+0.3535533905932737*alpha[8])*fUpwind[26]+0.3162277660168379*alpha[12]*fUpwind[25]+(0.2258769757263128*alpha[11]+0.3535533905932737*alpha[2])*fUpwind[24]+0.2020305089104422*alpha[23]*fUpwind[23]+0.3162277660168379*(alpha[14]*fUpwind[23]+fUpwind[14]*alpha[23])+0.3535533905932737*alpha[20]*fUpwind[22]+(0.2258769757263128*alpha[7]+0.3535533905932737*alpha[0])*fUpwind[21]+0.3162277660168379*alpha[4]*fUpwind[19]+0.2828427124746191*alpha[18]*fUpwind[18]+0.2020305089104422*alpha[17]*fUpwind[17]+0.3162277660168379*(alpha[6]*fUpwind[17]+fUpwind[6]*alpha[17])+0.3535533905932737*alpha[11]*fUpwind[16]+0.3162277660168379*alpha[1]*fUpwind[15]+0.2020305089104422*alpha[13]*fUpwind[13]+0.3162277660168379*(alpha[3]*fUpwind[13]+fUpwind[3]*alpha[13])+0.2828427124746191*alpha[10]*fUpwind[10]+0.3535533905932737*alpha[7]*fUpwind[9]+0.2828427124746191*alpha[5]*fUpwind[5]; 
  Ghat[22] += (0.2258769757263128*alpha[20]+0.3535533905932737*alpha[7])*fUpwind[26]+(0.2258769757263128*alpha[12]+0.3535533905932737*alpha[1])*fUpwind[25]+0.3162277660168379*alpha[11]*fUpwind[24]+0.2020305089104422*alpha[23]*fUpwind[23]+0.3162277660168379*(alpha[13]*fUpwind[23]+fUpwind[13]*alpha[23])+0.2258769757263128*alpha[8]*fUpwind[22]+0.3535533905932737*(alpha[0]*fUpwind[22]+alpha[20]*fUpwind[21])+0.3162277660168379*alpha[4]*fUpwind[19]+0.2020305089104422*alpha[18]*fUpwind[18]+0.3162277660168379*(alpha[5]*fUpwind[18]+fUpwind[5]*alpha[18])+0.2828427124746191*alpha[17]*fUpwind[17]+0.3162277660168379*alpha[2]*fUpwind[16]+0.3535533905932737*alpha[12]*fUpwind[15]+0.2020305089104422*alpha[14]*fUpwind[14]+0.3162277660168379*(alpha[3]*fUpwind[14]+fUpwind[3]*alpha[14])+0.2828427124746191*alpha[10]*fUpwind[10]+0.3535533905932737*alpha[8]*fUpwind[9]+0.2828427124746191*alpha[6]*fUpwind[6]; 
  Ghat[23] += (0.1290725575578931*alpha[23]+0.2020305089104422*(alpha[14]+alpha[13])+0.3162277660168379*alpha[3])*fUpwind[26]+(0.1807015805810503*alpha[18]+0.2828427124746191*alpha[5])*fUpwind[25]+(0.1807015805810503*alpha[17]+0.2828427124746191*alpha[6])*fUpwind[24]+(0.1443075063646015*alpha[20]+0.2258769757263128*(alpha[8]+alpha[7])+0.3535533905932737*alpha[0])*fUpwind[23]+(0.2020305089104422*(fUpwind[22]+fUpwind[21])+0.1443075063646015*fUpwind[20]+0.3162277660168379*fUpwind[9]+0.2258769757263128*(fUpwind[8]+fUpwind[7])+0.3535533905932737*fUpwind[0])*alpha[23]+0.3162277660168379*(alpha[13]*fUpwind[22]+alpha[14]*fUpwind[21])+(0.2258769757263128*(alpha[14]+alpha[13])+0.3535533905932737*alpha[3])*fUpwind[20]+(0.2258769757263128*(fUpwind[14]+fUpwind[13])+0.3535533905932737*fUpwind[3])*alpha[20]+0.2529822128134704*alpha[10]*fUpwind[19]+(0.2020305089104422*alpha[12]+0.3162277660168379*alpha[1])*fUpwind[18]+(0.282842712474619*fUpwind[15]+0.2020305089104422*fUpwind[12]+0.3162277660168379*fUpwind[1])*alpha[18]+(0.2020305089104422*alpha[11]+0.3162277660168379*alpha[2])*fUpwind[17]+(0.282842712474619*fUpwind[16]+0.2020305089104422*fUpwind[11]+0.3162277660168379*fUpwind[2])*alpha[17]+0.3535533905932737*(alpha[7]*fUpwind[14]+fUpwind[7]*alpha[14]+alpha[8]*fUpwind[13]+fUpwind[8]*alpha[13])+0.3162277660168379*(alpha[5]*fUpwind[12]+fUpwind[5]*alpha[12]+alpha[6]*fUpwind[11]+fUpwind[6]*alpha[11])+0.2828427124746191*(alpha[4]*fUpwind[10]+fUpwind[4]*alpha[10]); 
  Ghat[24] += (0.2020305089104422*alpha[11]+0.3162277660168379*alpha[2])*fUpwind[26]+0.2828427124746191*alpha[4]*fUpwind[25]+(0.2020305089104422*alpha[20]+0.3162277660168379*alpha[8]+0.2258769757263128*alpha[7]+0.3535533905932737*alpha[0])*fUpwind[24]+(0.1807015805810503*alpha[17]+0.2828427124746191*alpha[6])*fUpwind[23]+(0.1807015805810503*fUpwind[17]+0.2828427124746191*fUpwind[6])*alpha[23]+0.3162277660168379*alpha[11]*fUpwind[22]+(0.2258769757263128*alpha[11]+0.3535533905932737*alpha[2])*fUpwind[21]+0.3162277660168379*fUpwind[16]*alpha[20]+(0.282842712474619*alpha[12]+0.3162277660168379*alpha[1])*fUpwind[19]+0.2529822128134704*(alpha[10]*fUpwind[18]+fUpwind[10]*alpha[18])+(0.282842712474619*alpha[14]+0.2020305089104422*alpha[13]+0.3162277660168379*alpha[3])*fUpwind[17]+(0.282842712474619*fUpwind[14]+0.2020305089104422*fUpwind[13]+0.3162277660168379*fUpwind[3])*alpha[17]+0.3535533905932737*alpha[7]*fUpwind[16]+0.3162277660168379*(alpha[4]*fUpwind[15]+alpha[6]*fUpwind[13]+fUpwind[6]*alpha[13])+0.3535533905932737*fUpwind[9]*alpha[11]+0.2828427124746191*(alpha[5]*fUpwind[10]+fUpwind[5]*alpha[10]); 
  Ghat[25] += (0.2020305089104422*alpha[12]+0.3162277660168379*alpha[1])*fUpwind[26]+(0.2020305089104422*alpha[20]+0.2258769757263128*alpha[8]+0.3162277660168379*alpha[7]+0.3535533905932737*alpha[0])*fUpwind[25]+0.2828427124746191*alpha[4]*fUpwind[24]+(0.1807015805810503*alpha[18]+0.2828427124746191*alpha[5])*fUpwind[23]+(0.1807015805810503*fUpwind[18]+0.2828427124746191*fUpwind[5])*alpha[23]+(0.2258769757263128*alpha[12]+0.3535533905932737*alpha[1])*fUpwind[22]+0.3162277660168379*(alpha[12]*fUpwind[21]+fUpwind[15]*alpha[20])+(0.282842712474619*alpha[11]+0.3162277660168379*alpha[2])*fUpwind[19]+(0.2020305089104422*alpha[14]+0.282842712474619*alpha[13]+0.3162277660168379*alpha[3])*fUpwind[18]+(0.2020305089104422*fUpwind[14]+0.282842712474619*fUpwind[13]+0.3162277660168379*fUpwind[3])*alpha[18]+0.2529822128134704*(alpha[10]*fUpwind[17]+fUpwind[10]*alpha[17])+0.3162277660168379*alpha[4]*fUpwind[16]+0.3535533905932737*alpha[8]*fUpwind[15]+0.3162277660168379*(alpha[5]*fUpwind[14]+fUpwind[5]*alpha[14])+0.3535533905932737*fUpwind[9]*alpha[12]+0.2828427124746191*(alpha[6]*fUpwind[10]+fUpwind[6]*alpha[10]); 
  Ghat[26] += (0.1443075063646015*alpha[20]+0.2258769757263128*(alpha[8]+alpha[7])+0.3535533905932737*alpha[0])*fUpwind[26]+(0.2020305089104422*alpha[12]+0.3162277660168379*alpha[1])*fUpwind[25]+(0.2020305089104422*alpha[11]+0.3162277660168379*alpha[2])*fUpwind[24]+(0.1290725575578931*alpha[23]+0.2020305089104422*(alpha[14]+alpha[13])+0.3162277660168379*alpha[3])*fUpwind[23]+(0.2020305089104422*(fUpwind[14]+fUpwind[13])+0.3162277660168379*fUpwind[3])*alpha[23]+(0.2258769757263128*alpha[20]+0.3535533905932737*alpha[7])*fUpwind[22]+0.2258769757263128*alpha[20]*fUpwind[21]+0.3535533905932737*(alpha[8]*fUpwind[21]+fUpwind[9]*alpha[20])+0.2828427124746191*alpha[4]*fUpwind[19]+0.1807015805810503*alpha[18]*fUpwind[18]+0.2828427124746191*(alpha[5]*fUpwind[18]+fUpwind[5]*alpha[18])+0.1807015805810503*alpha[17]*fUpwind[17]+0.2828427124746191*(alpha[6]*fUpwind[17]+fUpwind[6]*alpha[17])+0.3162277660168379*(alpha[11]*fUpwind[16]+alpha[12]*fUpwind[15]+alpha[13]*fUpwind[14]+fUpwind[13]*alpha[14])+0.2529822128134704*alpha[10]*fUpwind[10]; 

  out[0] += -0.7071067811865475*Ghat[0]*dv11; 
  out[1] += -0.7071067811865475*Ghat[1]*dv11; 
  out[2] += -0.7071067811865475*Ghat[2]*dv11; 
  out[3] += -0.7071067811865475*Ghat[3]*dv11; 
  out[4] += -1.224744871391589*Ghat[0]*dv11; 
  out[5] += -0.7071067811865475*Ghat[4]*dv11; 
  out[6] += -0.7071067811865475*Ghat[5]*dv11; 
  out[7] += -0.7071067811865475*Ghat[6]*dv11; 
  out[8] += -1.224744871391589*Ghat[1]*dv11; 
  out[9] += -1.224744871391589*Ghat[2]*dv11; 
  out[10] += -1.224744871391589*Ghat[3]*dv11; 
  out[11] += -0.7071067811865475*Ghat[7]*dv11; 
  out[12] += -0.7071067811865475*Ghat[8]*dv11; 
  out[13] += -0.7071067811865475*Ghat[9]*dv11; 
  out[14] += -1.58113883008419*Ghat[0]*dv11; 
  out[15] += -0.7071067811865475*Ghat[10]*dv11; 
  out[16] += -1.224744871391589*Ghat[4]*dv11; 
  out[17] += -1.224744871391589*Ghat[5]*dv11; 
  out[18] += -1.224744871391589*Ghat[6]*dv11; 
  out[19] += -0.7071067811865475*Ghat[11]*dv11; 
  out[20] += -0.7071067811865475*Ghat[12]*dv11; 
  out[21] += -0.7071067811865475*Ghat[13]*dv11; 
  out[22] += -0.7071067811865475*Ghat[14]*dv11; 
  out[23] += -0.7071067811865475*Ghat[15]*dv11; 
  out[24] += -0.7071067811865475*Ghat[16]*dv11; 
  out[25] += -1.224744871391589*Ghat[7]*dv11; 
  out[26] += -1.224744871391589*Ghat[8]*dv11; 
  out[27] += -1.224744871391589*Ghat[9]*dv11; 
  out[28] += -1.58113883008419*Ghat[1]*dv11; 
  out[29] += -1.58113883008419*Ghat[2]*dv11; 
  out[30] += -1.58113883008419*Ghat[3]*dv11; 
  out[31] += -1.224744871391589*Ghat[10]*dv11; 
  out[32] += -0.7071067811865475*Ghat[17]*dv11; 
  out[33] += -0.7071067811865475*Ghat[18]*dv11; 
  out[34] += -0.7071067811865475*Ghat[19]*dv11; 
  out[35] += -1.224744871391589*Ghat[11]*dv11; 
  out[36] += -1.224744871391589*Ghat[12]*dv11; 
  out[37] += -1.224744871391589*Ghat[13]*dv11; 
  out[38] += -1.224744871391589*Ghat[14]*dv11; 
  out[39] += -1.224744871391589*Ghat[15]*dv11; 
  out[40] += -1.224744871391589*Ghat[16]*dv11; 
  out[41] += -1.58113883008419*Ghat[4]*dv11; 
  out[42] += -1.58113883008419*Ghat[5]*dv11; 
  out[43] += -1.58113883008419*Ghat[6]*dv11; 
  out[44] += -0.7071067811865475*Ghat[20]*dv11; 
  out[45] += -0.7071067811865475*Ghat[21]*dv11; 
  out[46] += -0.7071067811865475*Ghat[22]*dv11; 
  out[47] += -1.58113883008419*Ghat[7]*dv11; 
  out[48] += -1.58113883008419*Ghat[8]*dv11; 
  out[49] += -1.58113883008419*Ghat[9]*dv11; 
  out[50] += -1.224744871391589*Ghat[17]*dv11; 
  out[51] += -1.224744871391589*Ghat[18]*dv11; 
  out[52] += -1.224744871391589*Ghat[19]*dv11; 
  out[53] += -1.58113883008419*Ghat[10]*dv11; 
  out[54] += -0.7071067811865475*Ghat[23]*dv11; 
  out[55] += -0.7071067811865475*Ghat[24]*dv11; 
  out[56] += -0.7071067811865475*Ghat[25]*dv11; 
  out[57] += -1.224744871391589*Ghat[20]*dv11; 
  out[58] += -1.224744871391589*Ghat[21]*dv11; 
  out[59] += -1.224744871391589*Ghat[22]*dv11; 
  out[60] += -1.58113883008419*Ghat[11]*dv11; 
  out[61] += -1.58113883008419*Ghat[12]*dv11; 
  out[62] += -1.58113883008419*Ghat[13]*dv11; 
  out[63] += -1.58113883008419*Ghat[14]*dv11; 
  out[64] += -1.58113883008419*Ghat[15]*dv11; 
  out[65] += -1.58113883008419*Ghat[16]*dv11; 
  out[66] += -1.224744871391589*Ghat[23]*dv11; 
  out[67] += -1.224744871391589*Ghat[24]*dv11; 
  out[68] += -1.224744871391589*Ghat[25]*dv11; 
  out[69] += -1.58113883008419*Ghat[17]*dv11; 
  out[70] += -1.58113883008419*Ghat[18]*dv11; 
  out[71] += -1.58113883008419*Ghat[19]*dv11; 
  out[72] += -0.7071067811865475*Ghat[26]*dv11; 
  out[73] += -1.58113883008419*Ghat[20]*dv11; 
  out[74] += -1.58113883008419*Ghat[21]*dv11; 
  out[75] += -1.58113883008419*Ghat[22]*dv11; 
  out[76] += -1.224744871391589*Ghat[26]*dv11; 
  out[77] += -1.58113883008419*Ghat[23]*dv11; 
  out[78] += -1.58113883008419*Ghat[24]*dv11; 
  out[79] += -1.58113883008419*Ghat[25]*dv11; 
  out[80] += -1.58113883008419*Ghat[26]*dv11; 

  } else { 

  if ((-0.3794733192202055*alpha[23])+0.2828427124746191*alpha[20]+0.5692099788303082*(alpha[18]+alpha[17])-0.4242640687119285*(alpha[14]+alpha[13])-0.4242640687119281*alpha[12]-0.4242640687119285*alpha[11]-0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*(alpha[6]+alpha[5]+alpha[4])-0.4743416490252568*(alpha[3]+alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[0] = tensor_4x_p2_surfx4_quad_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = tensor_4x_p2_surfx4_quad_0_l(fSkin); 
  } 
  if (0.4743416490252568*alpha[23]-0.3535533905932737*alpha[20]-0.711512473537885*alpha[17]-0.4242640687119285*alpha[14]+0.5303300858899104*(alpha[13]+alpha[11])+0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]+0.6363961030678926*alpha[6]-0.4743416490252568*(alpha[3]+alpha[2])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[1] = tensor_4x_p2_surfx4_quad_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = tensor_4x_p2_surfx4_quad_1_l(fSkin); 
  } 
  if ((-0.3794733192202055*alpha[23])+0.2828427124746191*alpha[20]-0.5692099788303082*alpha[18]+0.5692099788303082*alpha[17]-0.4242640687119285*(alpha[14]+alpha[13])+0.4242640687119281*alpha[12]-0.4242640687119285*alpha[11]+0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*alpha[6]-0.6363961030678926*(alpha[5]+alpha[4])-0.4743416490252568*(alpha[3]+alpha[2])+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[2] = tensor_4x_p2_surfx4_quad_2_r(fEdge); 
  } else { 
    fUpwindQuad[2] = tensor_4x_p2_surfx4_quad_2_l(fSkin); 
  } 
  if (0.4743416490252568*alpha[23]-0.3535533905932737*alpha[20]-0.711512473537885*alpha[18]+0.5303300858899104*alpha[14]-0.4242640687119285*alpha[13]+0.5303300858899104*alpha[12]-0.3952847075210473*alpha[8]+0.3162277660168379*alpha[7]+0.6363961030678926*alpha[5]-0.4743416490252568*(alpha[3]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[3] = tensor_4x_p2_surfx4_quad_3_r(fEdge); 
  } else { 
    fUpwindQuad[3] = tensor_4x_p2_surfx4_quad_3_l(fSkin); 
  } 
  if ((-0.592927061281571*alpha[23])+0.441941738241592*alpha[20]+0.5303300858899104*(alpha[14]+alpha[13])-0.3952847075210473*(alpha[8]+alpha[7])-0.4743416490252568*alpha[3]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[4] = tensor_4x_p2_surfx4_quad_4_r(fEdge); 
  } else { 
    fUpwindQuad[4] = tensor_4x_p2_surfx4_quad_4_l(fSkin); 
  } 
  if (0.4743416490252568*alpha[23]-0.3535533905932737*alpha[20]+0.711512473537885*alpha[18]+0.5303300858899104*alpha[14]-0.4242640687119285*alpha[13]-0.5303300858899104*alpha[12]-0.3952847075210473*alpha[8]+0.3162277660168379*alpha[7]-0.6363961030678926*alpha[5]-0.4743416490252568*alpha[3]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[5] = tensor_4x_p2_surfx4_quad_5_r(fEdge); 
  } else { 
    fUpwindQuad[5] = tensor_4x_p2_surfx4_quad_5_l(fSkin); 
  } 
  if ((-0.3794733192202055*alpha[23])+0.2828427124746191*alpha[20]+0.5692099788303082*alpha[18]-0.5692099788303082*alpha[17]-0.4242640687119285*(alpha[14]+alpha[13])-0.4242640687119281*alpha[12]+0.4242640687119285*alpha[11]+0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*alpha[6]+0.6363961030678926*alpha[5]-0.6363961030678926*alpha[4]-0.4743416490252568*alpha[3]+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[6] = tensor_4x_p2_surfx4_quad_6_r(fEdge); 
  } else { 
    fUpwindQuad[6] = tensor_4x_p2_surfx4_quad_6_l(fSkin); 
  } 
  if (0.4743416490252568*alpha[23]-0.3535533905932737*alpha[20]+0.711512473537885*alpha[17]-0.4242640687119285*alpha[14]+0.5303300858899104*alpha[13]-0.5303300858899104*alpha[11]+0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]-0.6363961030678926*alpha[6]-0.4743416490252568*alpha[3]+0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[7] = tensor_4x_p2_surfx4_quad_7_r(fEdge); 
  } else { 
    fUpwindQuad[7] = tensor_4x_p2_surfx4_quad_7_l(fSkin); 
  } 
  if ((-0.3794733192202055*alpha[23])+0.2828427124746191*alpha[20]-0.5692099788303082*(alpha[18]+alpha[17])-0.4242640687119285*(alpha[14]+alpha[13])+0.4242640687119281*alpha[12]+0.4242640687119285*alpha[11]-0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*(alpha[6]+alpha[5])+0.6363961030678926*alpha[4]-0.4743416490252568*alpha[3]+0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[8] = tensor_4x_p2_surfx4_quad_8_r(fEdge); 
  } else { 
    fUpwindQuad[8] = tensor_4x_p2_surfx4_quad_8_l(fSkin); 
  } 
  if (0.2828427124746191*alpha[20]-0.4242640687119281*alpha[12]-0.4242640687119285*alpha[11]+0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*alpha[4]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[9] = tensor_4x_p2_surfx4_quad_9_r(fEdge); 
  } else { 
    fUpwindQuad[9] = tensor_4x_p2_surfx4_quad_9_l(fSkin); 
  } 
  if ((-0.3535533905932737*alpha[20])+0.5303300858899104*alpha[11]+0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]-0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[10] = tensor_4x_p2_surfx4_quad_10_r(fEdge); 
  } else { 
    fUpwindQuad[10] = tensor_4x_p2_surfx4_quad_10_l(fSkin); 
  } 
  if (0.2828427124746191*alpha[20]+0.4242640687119281*alpha[12]-0.4242640687119285*alpha[11]+0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*alpha[4]-0.4743416490252568*alpha[2]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[11] = tensor_4x_p2_surfx4_quad_11_r(fEdge); 
  } else { 
    fUpwindQuad[11] = tensor_4x_p2_surfx4_quad_11_l(fSkin); 
  } 
  if ((-0.3535533905932737*alpha[20])+0.5303300858899104*alpha[12]-0.3952847075210473*alpha[8]+0.3162277660168379*alpha[7]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[12] = tensor_4x_p2_surfx4_quad_12_r(fEdge); 
  } else { 
    fUpwindQuad[12] = tensor_4x_p2_surfx4_quad_12_l(fSkin); 
  } 
  if (0.441941738241592*alpha[20]-0.3952847075210473*(alpha[8]+alpha[7])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[13] = tensor_4x_p2_surfx4_quad_13_r(fEdge); 
  } else { 
    fUpwindQuad[13] = tensor_4x_p2_surfx4_quad_13_l(fSkin); 
  } 
  if ((-0.3535533905932737*alpha[20])-0.5303300858899104*alpha[12]-0.3952847075210473*alpha[8]+0.3162277660168379*alpha[7]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[14] = tensor_4x_p2_surfx4_quad_14_r(fEdge); 
  } else { 
    fUpwindQuad[14] = tensor_4x_p2_surfx4_quad_14_l(fSkin); 
  } 
  if (0.2828427124746191*alpha[20]-0.4242640687119281*alpha[12]+0.4242640687119285*alpha[11]+0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*alpha[4]+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[15] = tensor_4x_p2_surfx4_quad_15_r(fEdge); 
  } else { 
    fUpwindQuad[15] = tensor_4x_p2_surfx4_quad_15_l(fSkin); 
  } 
  if ((-0.3535533905932737*alpha[20])-0.5303300858899104*alpha[11]+0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]+0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[16] = tensor_4x_p2_surfx4_quad_16_r(fEdge); 
  } else { 
    fUpwindQuad[16] = tensor_4x_p2_surfx4_quad_16_l(fSkin); 
  } 
  if (0.2828427124746191*alpha[20]+0.4242640687119281*alpha[12]+0.4242640687119285*alpha[11]+0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*alpha[4]+0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[17] = tensor_4x_p2_surfx4_quad_17_r(fEdge); 
  } else { 
    fUpwindQuad[17] = tensor_4x_p2_surfx4_quad_17_l(fSkin); 
  } 
  if (0.3794733192202055*alpha[23]+0.2828427124746191*alpha[20]-0.5692099788303082*(alpha[18]+alpha[17])+0.4242640687119285*(alpha[14]+alpha[13])-0.4242640687119281*alpha[12]-0.4242640687119285*alpha[11]+0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*(alpha[6]+alpha[5])+0.6363961030678926*alpha[4]+0.4743416490252568*alpha[3]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[18] = tensor_4x_p2_surfx4_quad_18_r(fEdge); 
  } else { 
    fUpwindQuad[18] = tensor_4x_p2_surfx4_quad_18_l(fSkin); 
  } 
  if ((-0.4743416490252568*alpha[23])-0.3535533905932737*alpha[20]+0.711512473537885*alpha[17]+0.4242640687119285*alpha[14]-0.5303300858899104*alpha[13]+0.5303300858899104*alpha[11]+0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]-0.6363961030678926*alpha[6]+0.4743416490252568*alpha[3]-0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[19] = tensor_4x_p2_surfx4_quad_19_r(fEdge); 
  } else { 
    fUpwindQuad[19] = tensor_4x_p2_surfx4_quad_19_l(fSkin); 
  } 
  if (0.3794733192202055*alpha[23]+0.2828427124746191*alpha[20]+0.5692099788303082*alpha[18]-0.5692099788303082*alpha[17]+0.4242640687119285*(alpha[14]+alpha[13])+0.4242640687119281*alpha[12]-0.4242640687119285*alpha[11]-0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*alpha[6]+0.6363961030678926*alpha[5]-0.6363961030678926*alpha[4]+0.4743416490252568*alpha[3]-0.4743416490252568*alpha[2]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[20] = tensor_4x_p2_surfx4_quad_20_r(fEdge); 
  } else { 
    fUpwindQuad[20] = tensor_4x_p2_surfx4_quad_20_l(fSkin); 
  } 
  if ((-0.4743416490252568*alpha[23])-0.3535533905932737*alpha[20]+0.711512473537885*alpha[18]-0.5303300858899104*alpha[14]+0.4242640687119285*alpha[13]+0.5303300858899104*alpha[12]-0.3952847075210473*alpha[8]+0.3162277660168379*alpha[7]-0.6363961030678926*alpha[5]+0.4743416490252568*alpha[3]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[21] = tensor_4x_p2_surfx4_quad_21_r(fEdge); 
  } else { 
    fUpwindQuad[21] = tensor_4x_p2_surfx4_quad_21_l(fSkin); 
  } 
  if (0.592927061281571*alpha[23]+0.441941738241592*alpha[20]-0.5303300858899104*(alpha[14]+alpha[13])-0.3952847075210473*(alpha[8]+alpha[7])+0.4743416490252568*alpha[3]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[22] = tensor_4x_p2_surfx4_quad_22_r(fEdge); 
  } else { 
    fUpwindQuad[22] = tensor_4x_p2_surfx4_quad_22_l(fSkin); 
  } 
  if ((-0.4743416490252568*alpha[23])-0.3535533905932737*alpha[20]-0.711512473537885*alpha[18]-0.5303300858899104*alpha[14]+0.4242640687119285*alpha[13]-0.5303300858899104*alpha[12]-0.3952847075210473*alpha[8]+0.3162277660168379*alpha[7]+0.6363961030678926*alpha[5]+0.4743416490252568*(alpha[3]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[23] = tensor_4x_p2_surfx4_quad_23_r(fEdge); 
  } else { 
    fUpwindQuad[23] = tensor_4x_p2_surfx4_quad_23_l(fSkin); 
  } 
  if (0.3794733192202055*alpha[23]+0.2828427124746191*alpha[20]-0.5692099788303082*alpha[18]+0.5692099788303082*alpha[17]+0.4242640687119285*(alpha[14]+alpha[13])-0.4242640687119281*alpha[12]+0.4242640687119285*alpha[11]-0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*alpha[6]-0.6363961030678926*(alpha[5]+alpha[4])+0.4743416490252568*(alpha[3]+alpha[2])-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[24] = tensor_4x_p2_surfx4_quad_24_r(fEdge); 
  } else { 
    fUpwindQuad[24] = tensor_4x_p2_surfx4_quad_24_l(fSkin); 
  } 
  if ((-0.4743416490252568*alpha[23])-0.3535533905932737*alpha[20]-0.711512473537885*alpha[17]+0.4242640687119285*alpha[14]-0.5303300858899104*(alpha[13]+alpha[11])+0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]+0.6363961030678926*alpha[6]+0.4743416490252568*(alpha[3]+alpha[2])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[25] = tensor_4x_p2_surfx4_quad_25_r(fEdge); 
  } else { 
    fUpwindQuad[25] = tensor_4x_p2_surfx4_quad_25_l(fSkin); 
  } 
  if (0.3794733192202055*alpha[23]+0.2828427124746191*alpha[20]+0.5692099788303082*(alpha[18]+alpha[17])+0.4242640687119285*(alpha[14]+alpha[13])+0.4242640687119281*alpha[12]+0.4242640687119285*alpha[11]+0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*(alpha[6]+alpha[5]+alpha[4])+0.4743416490252568*(alpha[3]+alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[26] = tensor_4x_p2_surfx4_quad_26_r(fEdge); 
  } else { 
    fUpwindQuad[26] = tensor_4x_p2_surfx4_quad_26_l(fSkin); 
  } 

  // Project nodal basis back onto modal basis. 
  tensor_4x_p2_upwind(fUpwindQuad, fUpwind); 

  Ghat[0] += 0.3535533905932737*(alpha[23]*fUpwind[23]+alpha[20]*fUpwind[20]+alpha[18]*fUpwind[18]+alpha[17]*fUpwind[17]+alpha[14]*fUpwind[14]+alpha[13]*fUpwind[13]+alpha[12]*fUpwind[12]+alpha[11]*fUpwind[11]+alpha[10]*fUpwind[10]+alpha[8]*fUpwind[8]+alpha[7]*fUpwind[7]+alpha[6]*fUpwind[6]+alpha[5]*fUpwind[5]+alpha[4]*fUpwind[4]+alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] += 0.3162277660168379*(alpha[18]*fUpwind[23]+fUpwind[18]*alpha[23])+0.3162277660168379*(alpha[12]*fUpwind[20]+fUpwind[12]*alpha[20])+0.3535533905932737*(alpha[14]*fUpwind[18]+fUpwind[14]*alpha[18])+0.3162277660168379*(alpha[10]*fUpwind[17]+fUpwind[10]*alpha[17])+0.3162277660168379*(alpha[5]*fUpwind[13]+fUpwind[5]*alpha[13])+0.3535533905932737*(alpha[8]*fUpwind[12]+fUpwind[8]*alpha[12])+0.3162277660168379*(alpha[4]*fUpwind[11]+fUpwind[4]*alpha[11])+0.3535533905932737*(alpha[6]*fUpwind[10]+fUpwind[6]*alpha[10])+0.3162277660168379*(alpha[1]*fUpwind[7]+fUpwind[1]*alpha[7])+0.3535533905932737*(alpha[3]*fUpwind[5]+fUpwind[3]*alpha[5]+alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] += 0.3162277660168379*(alpha[17]*fUpwind[23]+fUpwind[17]*alpha[23])+0.3162277660168379*(alpha[11]*fUpwind[20]+fUpwind[11]*alpha[20])+0.3162277660168379*(alpha[10]*fUpwind[18]+fUpwind[10]*alpha[18])+0.3535533905932737*(alpha[13]*fUpwind[17]+fUpwind[13]*alpha[17])+0.3162277660168379*(alpha[6]*fUpwind[14]+fUpwind[6]*alpha[14]+alpha[4]*fUpwind[12]+fUpwind[4]*alpha[12])+0.3535533905932737*(alpha[7]*fUpwind[11]+fUpwind[7]*alpha[11]+alpha[5]*fUpwind[10]+fUpwind[5]*alpha[10])+0.3162277660168379*(alpha[2]*fUpwind[8]+fUpwind[2]*alpha[8])+0.3535533905932737*(alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6]+alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] += 0.3162277660168379*(alpha[23]*fUpwind[26]+alpha[18]*fUpwind[25]+alpha[17]*fUpwind[24])+0.3535533905932737*(alpha[20]*fUpwind[23]+fUpwind[20]*alpha[23])+0.3162277660168379*(alpha[14]*fUpwind[22]+alpha[13]*fUpwind[21])+0.3162277660168379*alpha[10]*fUpwind[19]+0.3535533905932737*(alpha[12]*fUpwind[18]+fUpwind[12]*alpha[18]+alpha[11]*fUpwind[17]+fUpwind[11]*alpha[17])+0.3162277660168379*(alpha[6]*fUpwind[16]+alpha[5]*fUpwind[15])+0.3535533905932737*(alpha[8]*fUpwind[14]+fUpwind[8]*alpha[14]+alpha[7]*fUpwind[13]+fUpwind[7]*alpha[13]+alpha[4]*fUpwind[10]+fUpwind[4]*alpha[10])+0.3162277660168379*alpha[3]*fUpwind[9]+0.3535533905932737*(alpha[2]*fUpwind[6]+fUpwind[2]*alpha[6]+alpha[1]*fUpwind[5]+fUpwind[1]*alpha[5]+alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]); 
  Ghat[4] += 0.2828427124746191*(alpha[10]*fUpwind[23]+fUpwind[10]*alpha[23]+alpha[4]*fUpwind[20]+fUpwind[4]*alpha[20])+(0.2828427124746191*alpha[17]+0.3162277660168379*alpha[6])*fUpwind[18]+0.2828427124746191*fUpwind[17]*alpha[18]+0.3162277660168379*(fUpwind[6]*alpha[18]+alpha[5]*fUpwind[17]+fUpwind[5]*alpha[17])+0.3162277660168379*(alpha[10]*fUpwind[14]+fUpwind[10]*alpha[14]+alpha[10]*fUpwind[13]+fUpwind[10]*alpha[13])+(0.2828427124746191*alpha[11]+0.3162277660168379*alpha[2])*fUpwind[12]+0.2828427124746191*fUpwind[11]*alpha[12]+0.3162277660168379*(fUpwind[2]*alpha[12]+alpha[1]*fUpwind[11]+fUpwind[1]*alpha[11])+0.3535533905932737*(alpha[3]*fUpwind[10]+fUpwind[3]*alpha[10])+0.3162277660168379*(alpha[4]*fUpwind[8]+fUpwind[4]*alpha[8]+alpha[4]*fUpwind[7]+fUpwind[4]*alpha[7])+0.3535533905932737*(alpha[5]*fUpwind[6]+fUpwind[5]*alpha[6]+alpha[0]*fUpwind[4]+fUpwind[0]*alpha[4]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[5] += 0.2828427124746191*alpha[18]*fUpwind[26]+(0.2828427124746191*alpha[23]+0.3162277660168379*alpha[14])*fUpwind[25]+0.2828427124746191*alpha[10]*fUpwind[24]+0.3162277660168379*(alpha[12]*fUpwind[23]+fUpwind[12]*alpha[23])+0.3162277660168379*alpha[18]*fUpwind[22]+0.2828427124746191*alpha[5]*fUpwind[21]+0.3162277660168379*(alpha[18]*fUpwind[20]+fUpwind[18]*alpha[20])+(0.2828427124746191*alpha[17]+0.3162277660168379*alpha[6])*fUpwind[19]+0.3535533905932737*(alpha[8]*fUpwind[18]+fUpwind[8]*alpha[18])+0.3162277660168379*(alpha[4]*fUpwind[17]+fUpwind[4]*alpha[17])+0.3162277660168379*alpha[10]*fUpwind[16]+(0.2828427124746191*alpha[13]+0.3162277660168379*alpha[3])*fUpwind[15]+0.3535533905932737*(alpha[12]*fUpwind[14]+fUpwind[12]*alpha[14])+0.3162277660168379*(alpha[1]*fUpwind[13]+fUpwind[1]*alpha[13]+alpha[10]*fUpwind[11]+fUpwind[10]*alpha[11])+0.3535533905932737*(alpha[2]*fUpwind[10]+fUpwind[2]*alpha[10])+0.3162277660168379*(alpha[5]*(fUpwind[9]+fUpwind[7])+fUpwind[5]*alpha[7])+0.3535533905932737*(alpha[4]*fUpwind[6]+fUpwind[4]*alpha[6]+alpha[0]*fUpwind[5]+fUpwind[0]*alpha[5]+alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]); 
  Ghat[6] += 0.2828427124746191*(alpha[17]*fUpwind[26]+alpha[10]*fUpwind[25]+alpha[23]*fUpwind[24])+0.3162277660168379*(alpha[13]*fUpwind[24]+alpha[11]*fUpwind[23]+fUpwind[11]*alpha[23])+0.2828427124746191*alpha[6]*fUpwind[22]+0.3162277660168379*(alpha[17]*(fUpwind[21]+fUpwind[20])+fUpwind[17]*alpha[20])+0.2828427124746191*alpha[18]*fUpwind[19]+0.3162277660168379*(alpha[5]*fUpwind[19]+alpha[4]*fUpwind[18]+fUpwind[4]*alpha[18])+0.3535533905932737*(alpha[7]*fUpwind[17]+fUpwind[7]*alpha[17])+0.2828427124746191*alpha[14]*fUpwind[16]+0.3162277660168379*(alpha[3]*fUpwind[16]+alpha[10]*fUpwind[15]+alpha[2]*fUpwind[14]+fUpwind[2]*alpha[14])+0.3535533905932737*(alpha[11]*fUpwind[13]+fUpwind[11]*alpha[13])+0.3162277660168379*(alpha[10]*fUpwind[12]+fUpwind[10]*alpha[12])+0.3535533905932737*(alpha[1]*fUpwind[10]+fUpwind[1]*alpha[10])+0.3162277660168379*(alpha[6]*(fUpwind[9]+fUpwind[8])+fUpwind[6]*alpha[8])+0.3535533905932737*(alpha[0]*fUpwind[6]+fUpwind[0]*alpha[6]+alpha[4]*fUpwind[5]+fUpwind[4]*alpha[5]+alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]); 
  Ghat[7] += 0.2258769757263128*alpha[23]*fUpwind[23]+0.3535533905932737*(alpha[14]*fUpwind[23]+fUpwind[14]*alpha[23])+0.2258769757263128*alpha[20]*fUpwind[20]+0.3535533905932737*(alpha[8]*fUpwind[20]+fUpwind[8]*alpha[20])+0.3162277660168379*alpha[18]*fUpwind[18]+0.2258769757263128*alpha[17]*fUpwind[17]+0.3535533905932737*(alpha[6]*fUpwind[17]+fUpwind[6]*alpha[17])+0.2258769757263128*alpha[13]*fUpwind[13]+0.3535533905932737*(alpha[3]*fUpwind[13]+fUpwind[3]*alpha[13])+0.3162277660168379*alpha[12]*fUpwind[12]+0.2258769757263128*alpha[11]*fUpwind[11]+0.3535533905932737*(alpha[2]*fUpwind[11]+fUpwind[2]*alpha[11])+0.3162277660168379*alpha[10]*fUpwind[10]+0.2258769757263128*alpha[7]*fUpwind[7]+0.3535533905932737*(alpha[0]*fUpwind[7]+fUpwind[0]*alpha[7])+0.3162277660168379*(alpha[5]*fUpwind[5]+alpha[4]*fUpwind[4]+alpha[1]*fUpwind[1]); 
  Ghat[8] += 0.2258769757263128*alpha[23]*fUpwind[23]+0.3535533905932737*(alpha[13]*fUpwind[23]+fUpwind[13]*alpha[23])+0.2258769757263128*alpha[20]*fUpwind[20]+0.3535533905932737*(alpha[7]*fUpwind[20]+fUpwind[7]*alpha[20])+0.2258769757263128*alpha[18]*fUpwind[18]+0.3535533905932737*(alpha[5]*fUpwind[18]+fUpwind[5]*alpha[18])+0.3162277660168379*alpha[17]*fUpwind[17]+0.2258769757263128*alpha[14]*fUpwind[14]+0.3535533905932737*(alpha[3]*fUpwind[14]+fUpwind[3]*alpha[14])+0.2258769757263128*alpha[12]*fUpwind[12]+0.3535533905932737*(alpha[1]*fUpwind[12]+fUpwind[1]*alpha[12])+0.3162277660168379*(alpha[11]*fUpwind[11]+alpha[10]*fUpwind[10])+0.2258769757263128*alpha[8]*fUpwind[8]+0.3535533905932737*(alpha[0]*fUpwind[8]+fUpwind[0]*alpha[8])+0.3162277660168379*(alpha[6]*fUpwind[6]+alpha[4]*fUpwind[4]+alpha[2]*fUpwind[2]); 
  Ghat[9] += 0.3535533905932737*(alpha[20]*fUpwind[26]+alpha[12]*fUpwind[25]+alpha[11]*fUpwind[24])+0.3162277660168379*alpha[23]*fUpwind[23]+0.3535533905932737*(alpha[8]*fUpwind[22]+alpha[7]*fUpwind[21]+alpha[4]*fUpwind[19])+0.3162277660168379*(alpha[18]*fUpwind[18]+alpha[17]*fUpwind[17])+0.3535533905932737*(alpha[2]*fUpwind[16]+alpha[1]*fUpwind[15])+0.3162277660168379*(alpha[14]*fUpwind[14]+alpha[13]*fUpwind[13]+alpha[10]*fUpwind[10])+0.3535533905932737*alpha[0]*fUpwind[9]+0.3162277660168379*(alpha[6]*fUpwind[6]+alpha[5]*fUpwind[5]+alpha[3]*fUpwind[3]); 
  Ghat[10] += 0.2529822128134704*alpha[10]*fUpwind[26]+(0.2529822128134704*alpha[17]+0.2828427124746191*alpha[6])*fUpwind[25]+0.2529822128134704*alpha[18]*fUpwind[24]+0.2828427124746191*(alpha[5]*fUpwind[24]+alpha[4]*fUpwind[23])+0.2529822128134704*fUpwind[19]*alpha[23]+0.2828427124746191*(fUpwind[4]*alpha[23]+alpha[10]*(fUpwind[22]+fUpwind[21]+fUpwind[20])+fUpwind[10]*alpha[20])+(0.282842712474619*(alpha[14]+alpha[13])+0.3162277660168379*alpha[3])*fUpwind[19]+(0.282842712474619*alpha[11]+0.3162277660168379*alpha[2])*fUpwind[18]+(0.282842712474619*(fUpwind[16]+fUpwind[11])+0.3162277660168379*fUpwind[2])*alpha[18]+(0.282842712474619*alpha[12]+0.3162277660168379*alpha[1])*fUpwind[17]+(0.282842712474619*(fUpwind[15]+fUpwind[12])+0.3162277660168379*fUpwind[1])*alpha[17]+0.3162277660168379*(alpha[5]*fUpwind[16]+alpha[6]*fUpwind[15]+alpha[4]*fUpwind[14]+fUpwind[4]*alpha[14]+alpha[4]*fUpwind[13]+fUpwind[4]*alpha[13]+alpha[6]*fUpwind[12]+fUpwind[6]*alpha[12]+alpha[5]*fUpwind[11]+fUpwind[5]*alpha[11])+(0.3162277660168379*(alpha[8]+alpha[7])+0.3535533905932737*alpha[0])*fUpwind[10]+0.3162277660168379*(fUpwind[9]+fUpwind[8]+fUpwind[7])*alpha[10]+0.3535533905932737*(fUpwind[0]*alpha[10]+alpha[1]*fUpwind[6]+fUpwind[1]*alpha[6]+alpha[2]*fUpwind[5]+fUpwind[2]*alpha[5]+alpha[3]*fUpwind[4]+fUpwind[3]*alpha[4]); 
  Ghat[11] += (0.2020305089104422*alpha[17]+0.3162277660168379*alpha[6])*fUpwind[23]+(0.2020305089104422*fUpwind[17]+0.3162277660168379*fUpwind[6])*alpha[23]+(0.2020305089104422*alpha[11]+0.3162277660168379*alpha[2])*fUpwind[20]+(0.2020305089104422*fUpwind[11]+0.3162277660168379*fUpwind[2])*alpha[20]+0.282842712474619*(alpha[10]*fUpwind[18]+fUpwind[10]*alpha[18])+(0.3162277660168379*alpha[14]+0.2258769757263128*alpha[13]+0.3535533905932737*alpha[3])*fUpwind[17]+(0.3162277660168379*fUpwind[14]+0.2258769757263128*fUpwind[13])*alpha[17]+0.3535533905932737*(fUpwind[3]*alpha[17]+alpha[6]*fUpwind[13]+fUpwind[6]*alpha[13])+0.2828427124746191*(alpha[4]*fUpwind[12]+fUpwind[4]*alpha[12])+(0.3162277660168379*alpha[8]+0.2258769757263128*alpha[7]+0.3535533905932737*alpha[0])*fUpwind[11]+(0.3162277660168379*fUpwind[8]+0.2258769757263128*fUpwind[7]+0.3535533905932737*fUpwind[0])*alpha[11]+0.3162277660168379*(alpha[5]*fUpwind[10]+fUpwind[5]*alpha[10])+0.3535533905932737*(alpha[2]*fUpwind[7]+fUpwind[2]*alpha[7])+0.3162277660168379*(alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]); 
  Ghat[12] += (0.2020305089104422*alpha[18]+0.3162277660168379*alpha[5])*fUpwind[23]+(0.2020305089104422*fUpwind[18]+0.3162277660168379*fUpwind[5])*alpha[23]+(0.2020305089104422*alpha[12]+0.3162277660168379*alpha[1])*fUpwind[20]+(0.2020305089104422*fUpwind[12]+0.3162277660168379*fUpwind[1])*alpha[20]+(0.2258769757263128*alpha[14]+0.3162277660168379*alpha[13]+0.3535533905932737*alpha[3])*fUpwind[18]+(0.2258769757263128*fUpwind[14]+0.3162277660168379*fUpwind[13]+0.3535533905932737*fUpwind[3])*alpha[18]+0.282842712474619*(alpha[10]*fUpwind[17]+fUpwind[10]*alpha[17])+0.3535533905932737*(alpha[5]*fUpwind[14]+fUpwind[5]*alpha[14])+(0.2258769757263128*alpha[8]+0.3162277660168379*alpha[7]+0.3535533905932737*alpha[0])*fUpwind[12]+(0.2258769757263128*fUpwind[8]+0.3162277660168379*fUpwind[7]+0.3535533905932737*fUpwind[0])*alpha[12]+0.2828427124746191*(alpha[4]*fUpwind[11]+fUpwind[4]*alpha[11])+0.3162277660168379*(alpha[6]*fUpwind[10]+fUpwind[6]*alpha[10])+0.3535533905932737*(alpha[1]*fUpwind[8]+fUpwind[1]*alpha[8])+0.3162277660168379*(alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]); 
  Ghat[13] += (0.2020305089104422*alpha[23]+0.3162277660168379*alpha[14])*fUpwind[26]+0.282842712474619*alpha[18]*fUpwind[25]+(0.2020305089104422*alpha[17]+0.3162277660168379*alpha[6])*fUpwind[24]+(0.2258769757263128*alpha[20]+0.3535533905932737*alpha[8])*fUpwind[23]+(0.3162277660168379*fUpwind[22]+0.2258769757263128*fUpwind[20]+0.3535533905932737*fUpwind[8])*alpha[23]+(0.2020305089104422*alpha[13]+0.3162277660168379*alpha[3])*fUpwind[21]+0.3535533905932737*(alpha[14]*fUpwind[20]+fUpwind[14]*alpha[20])+0.282842712474619*alpha[10]*fUpwind[19]+0.3162277660168379*(alpha[12]*fUpwind[18]+fUpwind[12]*alpha[18])+(0.2258769757263128*alpha[11]+0.3535533905932737*alpha[2])*fUpwind[17]+(0.3162277660168379*fUpwind[16]+0.2258769757263128*fUpwind[11]+0.3535533905932737*fUpwind[2])*alpha[17]+0.2828427124746191*alpha[5]*fUpwind[15]+(0.2258769757263128*alpha[7]+0.3535533905932737*alpha[0])*fUpwind[13]+(0.3162277660168379*fUpwind[9]+0.2258769757263128*fUpwind[7])*alpha[13]+0.3535533905932737*(fUpwind[0]*alpha[13]+alpha[6]*fUpwind[11]+fUpwind[6]*alpha[11])+0.3162277660168379*(alpha[4]*fUpwind[10]+fUpwind[4]*alpha[10])+0.3535533905932737*(alpha[3]*fUpwind[7]+fUpwind[3]*alpha[7])+0.3162277660168379*(alpha[1]*fUpwind[5]+fUpwind[1]*alpha[5]); 
  Ghat[14] += (0.2020305089104422*alpha[23]+0.3162277660168379*alpha[13])*fUpwind[26]+(0.2020305089104422*alpha[18]+0.3162277660168379*alpha[5])*fUpwind[25]+0.282842712474619*alpha[17]*fUpwind[24]+(0.2258769757263128*alpha[20]+0.3535533905932737*alpha[7])*fUpwind[23]+(0.3162277660168379*fUpwind[21]+0.2258769757263128*fUpwind[20]+0.3535533905932737*fUpwind[7])*alpha[23]+(0.2020305089104422*alpha[14]+0.3162277660168379*alpha[3])*fUpwind[22]+0.3535533905932737*(alpha[13]*fUpwind[20]+fUpwind[13]*alpha[20])+0.282842712474619*alpha[10]*fUpwind[19]+(0.2258769757263128*alpha[12]+0.3535533905932737*alpha[1])*fUpwind[18]+(0.3162277660168379*fUpwind[15]+0.2258769757263128*fUpwind[12]+0.3535533905932737*fUpwind[1])*alpha[18]+0.3162277660168379*(alpha[11]*fUpwind[17]+fUpwind[11]*alpha[17])+0.2828427124746191*alpha[6]*fUpwind[16]+(0.2258769757263128*alpha[8]+0.3535533905932737*alpha[0])*fUpwind[14]+(0.3162277660168379*fUpwind[9]+0.2258769757263128*fUpwind[8])*alpha[14]+0.3535533905932737*(fUpwind[0]*alpha[14]+alpha[5]*fUpwind[12]+fUpwind[5]*alpha[12])+0.3162277660168379*(alpha[4]*fUpwind[10]+fUpwind[4]*alpha[10])+0.3535533905932737*(alpha[3]*fUpwind[8]+fUpwind[3]*alpha[8])+0.3162277660168379*(alpha[2]*fUpwind[6]+fUpwind[2]*alpha[6]); 
  Ghat[15] += 0.3162277660168379*alpha[12]*fUpwind[26]+(0.3162277660168379*alpha[20]+0.3535533905932737*alpha[8])*fUpwind[25]+0.3162277660168379*alpha[4]*fUpwind[24]+0.282842712474619*(alpha[18]*fUpwind[23]+fUpwind[18]*alpha[23])+0.3535533905932737*alpha[12]*fUpwind[22]+0.3162277660168379*alpha[1]*fUpwind[21]+(0.3162277660168379*alpha[11]+0.3535533905932737*alpha[2])*fUpwind[19]+0.3162277660168379*(alpha[14]*fUpwind[18]+fUpwind[14]*alpha[18])+0.282842712474619*(alpha[10]*fUpwind[17]+fUpwind[10]*alpha[17])+0.3535533905932737*alpha[4]*fUpwind[16]+(0.3162277660168379*alpha[7]+0.3535533905932737*alpha[0])*fUpwind[15]+0.2828427124746191*(alpha[5]*fUpwind[13]+fUpwind[5]*alpha[13])+0.3162277660168379*(alpha[6]*fUpwind[10]+fUpwind[6]*alpha[10])+0.3535533905932737*alpha[1]*fUpwind[9]+0.3162277660168379*(alpha[3]*fUpwind[5]+fUpwind[3]*alpha[5]); 
  Ghat[16] += 0.3162277660168379*alpha[11]*fUpwind[26]+0.3162277660168379*alpha[4]*fUpwind[25]+(0.3162277660168379*alpha[20]+0.3535533905932737*alpha[7])*fUpwind[24]+0.282842712474619*(alpha[17]*fUpwind[23]+fUpwind[17]*alpha[23])+0.3162277660168379*alpha[2]*fUpwind[22]+0.3535533905932737*alpha[11]*fUpwind[21]+(0.3162277660168379*alpha[12]+0.3535533905932737*alpha[1])*fUpwind[19]+0.282842712474619*(alpha[10]*fUpwind[18]+fUpwind[10]*alpha[18])+0.3162277660168379*(alpha[13]*fUpwind[17]+fUpwind[13]*alpha[17]+alpha[8]*fUpwind[16])+0.3535533905932737*(alpha[0]*fUpwind[16]+alpha[4]*fUpwind[15])+0.2828427124746191*(alpha[6]*fUpwind[14]+fUpwind[6]*alpha[14])+0.3162277660168379*(alpha[5]*fUpwind[10]+fUpwind[5]*alpha[10])+0.3535533905932737*alpha[2]*fUpwind[9]+0.3162277660168379*(alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6]); 
  Ghat[17] += (0.1807015805810503*alpha[17]+0.2828427124746191*alpha[6])*fUpwind[26]+0.2529822128134704*alpha[10]*fUpwind[25]+(0.1807015805810503*alpha[23]+0.282842712474619*alpha[14]+0.2020305089104422*alpha[13]+0.3162277660168379*alpha[3])*fUpwind[24]+(0.2020305089104422*alpha[11]+0.3162277660168379*alpha[2])*fUpwind[23]+(0.282842712474619*fUpwind[16]+0.2020305089104422*fUpwind[11]+0.3162277660168379*fUpwind[2])*alpha[23]+0.2828427124746191*alpha[17]*fUpwind[22]+(0.2020305089104422*alpha[17]+0.3162277660168379*alpha[6])*fUpwind[21]+(0.2020305089104422*alpha[17]+0.3162277660168379*alpha[6])*fUpwind[20]+(0.2020305089104422*fUpwind[17]+0.3162277660168379*fUpwind[6])*alpha[20]+0.2529822128134704*alpha[18]*fUpwind[19]+0.2828427124746191*(alpha[5]*fUpwind[19]+alpha[4]*fUpwind[18]+fUpwind[4]*alpha[18])+(0.3162277660168379*alpha[8]+0.2258769757263128*alpha[7]+0.3535533905932737*alpha[0])*fUpwind[17]+(0.3162277660168379*(fUpwind[9]+fUpwind[8])+0.2258769757263128*fUpwind[7]+0.3535533905932737*fUpwind[0])*alpha[17]+0.3162277660168379*alpha[13]*fUpwind[16]+0.282842712474619*alpha[10]*fUpwind[15]+0.3162277660168379*(alpha[11]*fUpwind[14]+fUpwind[11]*alpha[14])+(0.2258769757263128*alpha[11]+0.3535533905932737*alpha[2])*fUpwind[13]+(0.2258769757263128*fUpwind[11]+0.3535533905932737*fUpwind[2])*alpha[13]+0.282842712474619*(alpha[10]*fUpwind[12]+fUpwind[10]*alpha[12])+0.3535533905932737*(alpha[3]*fUpwind[11]+fUpwind[3]*alpha[11])+0.3162277660168379*(alpha[1]*fUpwind[10]+fUpwind[1]*alpha[10])+0.3535533905932737*(alpha[6]*fUpwind[7]+fUpwind[6]*alpha[7])+0.3162277660168379*(alpha[4]*fUpwind[5]+fUpwind[4]*alpha[5]); 
  Ghat[18] += (0.1807015805810503*alpha[18]+0.2828427124746191*alpha[5])*fUpwind[26]+(0.1807015805810503*alpha[23]+0.2020305089104422*alpha[14]+0.282842712474619*alpha[13]+0.3162277660168379*alpha[3])*fUpwind[25]+0.2529822128134704*alpha[10]*fUpwind[24]+(0.2020305089104422*alpha[12]+0.3162277660168379*alpha[1])*fUpwind[23]+(0.282842712474619*fUpwind[15]+0.2020305089104422*fUpwind[12]+0.3162277660168379*fUpwind[1])*alpha[23]+(0.2020305089104422*alpha[18]+0.3162277660168379*alpha[5])*fUpwind[22]+0.2828427124746191*alpha[18]*fUpwind[21]+(0.2020305089104422*alpha[18]+0.3162277660168379*alpha[5])*fUpwind[20]+(0.2020305089104422*fUpwind[18]+0.3162277660168379*fUpwind[5])*alpha[20]+(0.2529822128134704*alpha[17]+0.2828427124746191*alpha[6])*fUpwind[19]+(0.2258769757263128*alpha[8]+0.3162277660168379*alpha[7]+0.3535533905932737*alpha[0])*fUpwind[18]+(0.3162277660168379*fUpwind[9]+0.2258769757263128*fUpwind[8]+0.3162277660168379*fUpwind[7]+0.3535533905932737*fUpwind[0])*alpha[18]+0.2828427124746191*(alpha[4]*fUpwind[17]+fUpwind[4]*alpha[17])+0.282842712474619*alpha[10]*fUpwind[16]+0.3162277660168379*alpha[14]*fUpwind[15]+(0.2258769757263128*alpha[12]+0.3535533905932737*alpha[1])*fUpwind[14]+(0.2258769757263128*fUpwind[12]+0.3535533905932737*fUpwind[1])*alpha[14]+0.3162277660168379*(alpha[12]*fUpwind[13]+fUpwind[12]*alpha[13])+0.3535533905932737*(alpha[3]*fUpwind[12]+fUpwind[3]*alpha[12])+0.282842712474619*(alpha[10]*fUpwind[11]+fUpwind[10]*alpha[11])+0.3162277660168379*(alpha[2]*fUpwind[10]+fUpwind[2]*alpha[10])+0.3535533905932737*(alpha[5]*fUpwind[8]+fUpwind[5]*alpha[8])+0.3162277660168379*(alpha[4]*fUpwind[6]+fUpwind[4]*alpha[6]); 
  Ghat[19] += 0.2828427124746191*alpha[4]*fUpwind[26]+(0.282842712474619*alpha[11]+0.3162277660168379*alpha[2])*fUpwind[25]+(0.282842712474619*alpha[12]+0.3162277660168379*alpha[1])*fUpwind[24]+0.2529822128134704*(alpha[10]*fUpwind[23]+fUpwind[10]*alpha[23])+0.3162277660168379*alpha[4]*(fUpwind[22]+fUpwind[21])+fUpwind[19]*(0.2828427124746191*alpha[20]+0.3162277660168379*(alpha[8]+alpha[7])+0.3535533905932737*alpha[0])+(0.2529822128134704*alpha[17]+0.2828427124746191*alpha[6])*fUpwind[18]+0.2529822128134704*fUpwind[17]*alpha[18]+0.2828427124746191*(fUpwind[6]*alpha[18]+alpha[5]*fUpwind[17]+fUpwind[5]*alpha[17])+(0.3162277660168379*alpha[12]+0.3535533905932737*alpha[1])*fUpwind[16]+(0.3162277660168379*alpha[11]+0.3535533905932737*alpha[2])*fUpwind[15]+0.282842712474619*(alpha[10]*fUpwind[14]+fUpwind[10]*alpha[14]+alpha[10]*fUpwind[13]+fUpwind[10]*alpha[13])+0.3162277660168379*(alpha[3]*fUpwind[10]+fUpwind[3]*alpha[10])+0.3535533905932737*alpha[4]*fUpwind[9]+0.3162277660168379*(alpha[5]*fUpwind[6]+fUpwind[5]*alpha[6]); 
  Ghat[20] += (0.1443075063646015*alpha[23]+0.2258769757263128*(alpha[14]+alpha[13])+0.3535533905932737*alpha[3])*fUpwind[23]+(0.2258769757263128*(fUpwind[14]+fUpwind[13])+0.3535533905932737*fUpwind[3])*alpha[23]+(0.1443075063646015*alpha[20]+0.2258769757263128*(alpha[8]+alpha[7])+0.3535533905932737*alpha[0])*fUpwind[20]+(0.2258769757263128*(fUpwind[8]+fUpwind[7])+0.3535533905932737*fUpwind[0])*alpha[20]+0.2020305089104422*alpha[18]*fUpwind[18]+0.3162277660168379*(alpha[5]*fUpwind[18]+fUpwind[5]*alpha[18])+0.2020305089104422*alpha[17]*fUpwind[17]+0.3162277660168379*(alpha[6]*fUpwind[17]+fUpwind[6]*alpha[17])+0.3535533905932737*(alpha[13]*fUpwind[14]+fUpwind[13]*alpha[14])+0.2020305089104422*alpha[12]*fUpwind[12]+0.3162277660168379*(alpha[1]*fUpwind[12]+fUpwind[1]*alpha[12])+0.2020305089104422*alpha[11]*fUpwind[11]+0.3162277660168379*(alpha[2]*fUpwind[11]+fUpwind[2]*alpha[11])+0.2828427124746191*alpha[10]*fUpwind[10]+0.3535533905932737*(alpha[7]*fUpwind[8]+fUpwind[7]*alpha[8])+0.2828427124746191*alpha[4]*fUpwind[4]; 
  Ghat[21] += (0.2258769757263128*alpha[20]+0.3535533905932737*alpha[8])*fUpwind[26]+0.3162277660168379*alpha[12]*fUpwind[25]+(0.2258769757263128*alpha[11]+0.3535533905932737*alpha[2])*fUpwind[24]+0.2020305089104422*alpha[23]*fUpwind[23]+0.3162277660168379*(alpha[14]*fUpwind[23]+fUpwind[14]*alpha[23])+0.3535533905932737*alpha[20]*fUpwind[22]+(0.2258769757263128*alpha[7]+0.3535533905932737*alpha[0])*fUpwind[21]+0.3162277660168379*alpha[4]*fUpwind[19]+0.2828427124746191*alpha[18]*fUpwind[18]+0.2020305089104422*alpha[17]*fUpwind[17]+0.3162277660168379*(alpha[6]*fUpwind[17]+fUpwind[6]*alpha[17])+0.3535533905932737*alpha[11]*fUpwind[16]+0.3162277660168379*alpha[1]*fUpwind[15]+0.2020305089104422*alpha[13]*fUpwind[13]+0.3162277660168379*(alpha[3]*fUpwind[13]+fUpwind[3]*alpha[13])+0.2828427124746191*alpha[10]*fUpwind[10]+0.3535533905932737*alpha[7]*fUpwind[9]+0.2828427124746191*alpha[5]*fUpwind[5]; 
  Ghat[22] += (0.2258769757263128*alpha[20]+0.3535533905932737*alpha[7])*fUpwind[26]+(0.2258769757263128*alpha[12]+0.3535533905932737*alpha[1])*fUpwind[25]+0.3162277660168379*alpha[11]*fUpwind[24]+0.2020305089104422*alpha[23]*fUpwind[23]+0.3162277660168379*(alpha[13]*fUpwind[23]+fUpwind[13]*alpha[23])+0.2258769757263128*alpha[8]*fUpwind[22]+0.3535533905932737*(alpha[0]*fUpwind[22]+alpha[20]*fUpwind[21])+0.3162277660168379*alpha[4]*fUpwind[19]+0.2020305089104422*alpha[18]*fUpwind[18]+0.3162277660168379*(alpha[5]*fUpwind[18]+fUpwind[5]*alpha[18])+0.2828427124746191*alpha[17]*fUpwind[17]+0.3162277660168379*alpha[2]*fUpwind[16]+0.3535533905932737*alpha[12]*fUpwind[15]+0.2020305089104422*alpha[14]*fUpwind[14]+0.3162277660168379*(alpha[3]*fUpwind[14]+fUpwind[3]*alpha[14])+0.2828427124746191*alpha[10]*fUpwind[10]+0.3535533905932737*alpha[8]*fUpwind[9]+0.2828427124746191*alpha[6]*fUpwind[6]; 
  Ghat[23] += (0.1290725575578931*alpha[23]+0.2020305089104422*(alpha[14]+alpha[13])+0.3162277660168379*alpha[3])*fUpwind[26]+(0.1807015805810503*alpha[18]+0.2828427124746191*alpha[5])*fUpwind[25]+(0.1807015805810503*alpha[17]+0.2828427124746191*alpha[6])*fUpwind[24]+(0.1443075063646015*alpha[20]+0.2258769757263128*(alpha[8]+alpha[7])+0.3535533905932737*alpha[0])*fUpwind[23]+(0.2020305089104422*(fUpwind[22]+fUpwind[21])+0.1443075063646015*fUpwind[20]+0.3162277660168379*fUpwind[9]+0.2258769757263128*(fUpwind[8]+fUpwind[7])+0.3535533905932737*fUpwind[0])*alpha[23]+0.3162277660168379*(alpha[13]*fUpwind[22]+alpha[14]*fUpwind[21])+(0.2258769757263128*(alpha[14]+alpha[13])+0.3535533905932737*alpha[3])*fUpwind[20]+(0.2258769757263128*(fUpwind[14]+fUpwind[13])+0.3535533905932737*fUpwind[3])*alpha[20]+0.2529822128134704*alpha[10]*fUpwind[19]+(0.2020305089104422*alpha[12]+0.3162277660168379*alpha[1])*fUpwind[18]+(0.282842712474619*fUpwind[15]+0.2020305089104422*fUpwind[12]+0.3162277660168379*fUpwind[1])*alpha[18]+(0.2020305089104422*alpha[11]+0.3162277660168379*alpha[2])*fUpwind[17]+(0.282842712474619*fUpwind[16]+0.2020305089104422*fUpwind[11]+0.3162277660168379*fUpwind[2])*alpha[17]+0.3535533905932737*(alpha[7]*fUpwind[14]+fUpwind[7]*alpha[14]+alpha[8]*fUpwind[13]+fUpwind[8]*alpha[13])+0.3162277660168379*(alpha[5]*fUpwind[12]+fUpwind[5]*alpha[12]+alpha[6]*fUpwind[11]+fUpwind[6]*alpha[11])+0.2828427124746191*(alpha[4]*fUpwind[10]+fUpwind[4]*alpha[10]); 
  Ghat[24] += (0.2020305089104422*alpha[11]+0.3162277660168379*alpha[2])*fUpwind[26]+0.2828427124746191*alpha[4]*fUpwind[25]+(0.2020305089104422*alpha[20]+0.3162277660168379*alpha[8]+0.2258769757263128*alpha[7]+0.3535533905932737*alpha[0])*fUpwind[24]+(0.1807015805810503*alpha[17]+0.2828427124746191*alpha[6])*fUpwind[23]+(0.1807015805810503*fUpwind[17]+0.2828427124746191*fUpwind[6])*alpha[23]+0.3162277660168379*alpha[11]*fUpwind[22]+(0.2258769757263128*alpha[11]+0.3535533905932737*alpha[2])*fUpwind[21]+0.3162277660168379*fUpwind[16]*alpha[20]+(0.282842712474619*alpha[12]+0.3162277660168379*alpha[1])*fUpwind[19]+0.2529822128134704*(alpha[10]*fUpwind[18]+fUpwind[10]*alpha[18])+(0.282842712474619*alpha[14]+0.2020305089104422*alpha[13]+0.3162277660168379*alpha[3])*fUpwind[17]+(0.282842712474619*fUpwind[14]+0.2020305089104422*fUpwind[13]+0.3162277660168379*fUpwind[3])*alpha[17]+0.3535533905932737*alpha[7]*fUpwind[16]+0.3162277660168379*(alpha[4]*fUpwind[15]+alpha[6]*fUpwind[13]+fUpwind[6]*alpha[13])+0.3535533905932737*fUpwind[9]*alpha[11]+0.2828427124746191*(alpha[5]*fUpwind[10]+fUpwind[5]*alpha[10]); 
  Ghat[25] += (0.2020305089104422*alpha[12]+0.3162277660168379*alpha[1])*fUpwind[26]+(0.2020305089104422*alpha[20]+0.2258769757263128*alpha[8]+0.3162277660168379*alpha[7]+0.3535533905932737*alpha[0])*fUpwind[25]+0.2828427124746191*alpha[4]*fUpwind[24]+(0.1807015805810503*alpha[18]+0.2828427124746191*alpha[5])*fUpwind[23]+(0.1807015805810503*fUpwind[18]+0.2828427124746191*fUpwind[5])*alpha[23]+(0.2258769757263128*alpha[12]+0.3535533905932737*alpha[1])*fUpwind[22]+0.3162277660168379*(alpha[12]*fUpwind[21]+fUpwind[15]*alpha[20])+(0.282842712474619*alpha[11]+0.3162277660168379*alpha[2])*fUpwind[19]+(0.2020305089104422*alpha[14]+0.282842712474619*alpha[13]+0.3162277660168379*alpha[3])*fUpwind[18]+(0.2020305089104422*fUpwind[14]+0.282842712474619*fUpwind[13]+0.3162277660168379*fUpwind[3])*alpha[18]+0.2529822128134704*(alpha[10]*fUpwind[17]+fUpwind[10]*alpha[17])+0.3162277660168379*alpha[4]*fUpwind[16]+0.3535533905932737*alpha[8]*fUpwind[15]+0.3162277660168379*(alpha[5]*fUpwind[14]+fUpwind[5]*alpha[14])+0.3535533905932737*fUpwind[9]*alpha[12]+0.2828427124746191*(alpha[6]*fUpwind[10]+fUpwind[6]*alpha[10]); 
  Ghat[26] += (0.1443075063646015*alpha[20]+0.2258769757263128*(alpha[8]+alpha[7])+0.3535533905932737*alpha[0])*fUpwind[26]+(0.2020305089104422*alpha[12]+0.3162277660168379*alpha[1])*fUpwind[25]+(0.2020305089104422*alpha[11]+0.3162277660168379*alpha[2])*fUpwind[24]+(0.1290725575578931*alpha[23]+0.2020305089104422*(alpha[14]+alpha[13])+0.3162277660168379*alpha[3])*fUpwind[23]+(0.2020305089104422*(fUpwind[14]+fUpwind[13])+0.3162277660168379*fUpwind[3])*alpha[23]+(0.2258769757263128*alpha[20]+0.3535533905932737*alpha[7])*fUpwind[22]+0.2258769757263128*alpha[20]*fUpwind[21]+0.3535533905932737*(alpha[8]*fUpwind[21]+fUpwind[9]*alpha[20])+0.2828427124746191*alpha[4]*fUpwind[19]+0.1807015805810503*alpha[18]*fUpwind[18]+0.2828427124746191*(alpha[5]*fUpwind[18]+fUpwind[5]*alpha[18])+0.1807015805810503*alpha[17]*fUpwind[17]+0.2828427124746191*(alpha[6]*fUpwind[17]+fUpwind[6]*alpha[17])+0.3162277660168379*(alpha[11]*fUpwind[16]+alpha[12]*fUpwind[15]+alpha[13]*fUpwind[14]+fUpwind[13]*alpha[14])+0.2529822128134704*alpha[10]*fUpwind[10]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv11; 
  out[1] += 0.7071067811865475*Ghat[1]*dv11; 
  out[2] += 0.7071067811865475*Ghat[2]*dv11; 
  out[3] += 0.7071067811865475*Ghat[3]*dv11; 
  out[4] += -1.224744871391589*Ghat[0]*dv11; 
  out[5] += 0.7071067811865475*Ghat[4]*dv11; 
  out[6] += 0.7071067811865475*Ghat[5]*dv11; 
  out[7] += 0.7071067811865475*Ghat[6]*dv11; 
  out[8] += -1.224744871391589*Ghat[1]*dv11; 
  out[9] += -1.224744871391589*Ghat[2]*dv11; 
  out[10] += -1.224744871391589*Ghat[3]*dv11; 
  out[11] += 0.7071067811865475*Ghat[7]*dv11; 
  out[12] += 0.7071067811865475*Ghat[8]*dv11; 
  out[13] += 0.7071067811865475*Ghat[9]*dv11; 
  out[14] += 1.58113883008419*Ghat[0]*dv11; 
  out[15] += 0.7071067811865475*Ghat[10]*dv11; 
  out[16] += -1.224744871391589*Ghat[4]*dv11; 
  out[17] += -1.224744871391589*Ghat[5]*dv11; 
  out[18] += -1.224744871391589*Ghat[6]*dv11; 
  out[19] += 0.7071067811865475*Ghat[11]*dv11; 
  out[20] += 0.7071067811865475*Ghat[12]*dv11; 
  out[21] += 0.7071067811865475*Ghat[13]*dv11; 
  out[22] += 0.7071067811865475*Ghat[14]*dv11; 
  out[23] += 0.7071067811865475*Ghat[15]*dv11; 
  out[24] += 0.7071067811865475*Ghat[16]*dv11; 
  out[25] += -1.224744871391589*Ghat[7]*dv11; 
  out[26] += -1.224744871391589*Ghat[8]*dv11; 
  out[27] += -1.224744871391589*Ghat[9]*dv11; 
  out[28] += 1.58113883008419*Ghat[1]*dv11; 
  out[29] += 1.58113883008419*Ghat[2]*dv11; 
  out[30] += 1.58113883008419*Ghat[3]*dv11; 
  out[31] += -1.224744871391589*Ghat[10]*dv11; 
  out[32] += 0.7071067811865475*Ghat[17]*dv11; 
  out[33] += 0.7071067811865475*Ghat[18]*dv11; 
  out[34] += 0.7071067811865475*Ghat[19]*dv11; 
  out[35] += -1.224744871391589*Ghat[11]*dv11; 
  out[36] += -1.224744871391589*Ghat[12]*dv11; 
  out[37] += -1.224744871391589*Ghat[13]*dv11; 
  out[38] += -1.224744871391589*Ghat[14]*dv11; 
  out[39] += -1.224744871391589*Ghat[15]*dv11; 
  out[40] += -1.224744871391589*Ghat[16]*dv11; 
  out[41] += 1.58113883008419*Ghat[4]*dv11; 
  out[42] += 1.58113883008419*Ghat[5]*dv11; 
  out[43] += 1.58113883008419*Ghat[6]*dv11; 
  out[44] += 0.7071067811865475*Ghat[20]*dv11; 
  out[45] += 0.7071067811865475*Ghat[21]*dv11; 
  out[46] += 0.7071067811865475*Ghat[22]*dv11; 
  out[47] += 1.58113883008419*Ghat[7]*dv11; 
  out[48] += 1.58113883008419*Ghat[8]*dv11; 
  out[49] += 1.58113883008419*Ghat[9]*dv11; 
  out[50] += -1.224744871391589*Ghat[17]*dv11; 
  out[51] += -1.224744871391589*Ghat[18]*dv11; 
  out[52] += -1.224744871391589*Ghat[19]*dv11; 
  out[53] += 1.58113883008419*Ghat[10]*dv11; 
  out[54] += 0.7071067811865475*Ghat[23]*dv11; 
  out[55] += 0.7071067811865475*Ghat[24]*dv11; 
  out[56] += 0.7071067811865475*Ghat[25]*dv11; 
  out[57] += -1.224744871391589*Ghat[20]*dv11; 
  out[58] += -1.224744871391589*Ghat[21]*dv11; 
  out[59] += -1.224744871391589*Ghat[22]*dv11; 
  out[60] += 1.58113883008419*Ghat[11]*dv11; 
  out[61] += 1.58113883008419*Ghat[12]*dv11; 
  out[62] += 1.58113883008419*Ghat[13]*dv11; 
  out[63] += 1.58113883008419*Ghat[14]*dv11; 
  out[64] += 1.58113883008419*Ghat[15]*dv11; 
  out[65] += 1.58113883008419*Ghat[16]*dv11; 
  out[66] += -1.224744871391589*Ghat[23]*dv11; 
  out[67] += -1.224744871391589*Ghat[24]*dv11; 
  out[68] += -1.224744871391589*Ghat[25]*dv11; 
  out[69] += 1.58113883008419*Ghat[17]*dv11; 
  out[70] += 1.58113883008419*Ghat[18]*dv11; 
  out[71] += 1.58113883008419*Ghat[19]*dv11; 
  out[72] += 0.7071067811865475*Ghat[26]*dv11; 
  out[73] += 1.58113883008419*Ghat[20]*dv11; 
  out[74] += 1.58113883008419*Ghat[21]*dv11; 
  out[75] += 1.58113883008419*Ghat[22]*dv11; 
  out[76] += -1.224744871391589*Ghat[26]*dv11; 
  out[77] += 1.58113883008419*Ghat[23]*dv11; 
  out[78] += 1.58113883008419*Ghat[24]*dv11; 
  out[79] += 1.58113883008419*Ghat[25]*dv11; 
  out[80] += 1.58113883008419*Ghat[26]*dv11; 

  } 
} 
