#include <gkyl_vlasov_sr_kernels.h> 
#include <gkyl_basis_ser_6x_p1_surfx5_eval_quad.h> 
#include <gkyl_basis_ser_6x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_sr_boundary_surfvy_3x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // p_over_gamma: p/gamma (velocity).
  // qmem:        q/m*EM fields.
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv11 = 2/dxv[4]; 
  const double dv1 = dxv[3], wv1 = w[3]; 
  const double dv2 = dxv[4], wv2 = w[4]; 
  const double dv3 = dxv[5], wv3 = w[5]; 
  const double *E1 = &qmem[8]; 
  const double *p0_over_gamma = &p_over_gamma[0]; 
  const double *B0 = &qmem[24]; 
  const double *p1_over_gamma = &p_over_gamma[8]; 
  const double *B1 = &qmem[32]; 
  const double *p2_over_gamma = &p_over_gamma[16]; 
  const double *B2 = &qmem[40]; 

  double alpha[32] = {0.0}; 

  double fUpwindQuad[32] = {0.0};
  double fUpwind[32] = {0.0};
  double Ghat[32] = {0.0}; 

  if (edge == -1) { 

  alpha[0] = 1.224744871391589*B0[0]*p2_over_gamma[2]-1.224744871391589*B2[0]*p0_over_gamma[2]+0.7071067811865475*B0[0]*p2_over_gamma[0]-0.7071067811865475*B2[0]*p0_over_gamma[0]+2.0*E1[0]; 
  alpha[1] = 1.224744871391589*B0[1]*p2_over_gamma[2]-1.224744871391589*B2[1]*p0_over_gamma[2]+2.0*E1[1]-0.7071067811865475*p0_over_gamma[0]*B2[1]+0.7071067811865475*p2_over_gamma[0]*B0[1]; 
  alpha[2] = 1.224744871391589*B0[2]*p2_over_gamma[2]-1.224744871391589*B2[2]*p0_over_gamma[2]+2.0*E1[2]-0.7071067811865475*p0_over_gamma[0]*B2[2]+0.7071067811865475*p2_over_gamma[0]*B0[2]; 
  alpha[3] = 2.0*E1[3]-1.224744871391589*p0_over_gamma[2]*B2[3]-0.7071067811865475*p0_over_gamma[0]*B2[3]+1.224744871391589*p2_over_gamma[2]*B0[3]+0.7071067811865475*p2_over_gamma[0]*B0[3]; 
  alpha[4] = 1.224744871391589*B0[0]*p2_over_gamma[4]-1.224744871391589*B2[0]*p0_over_gamma[4]+0.7071067811865475*B0[0]*p2_over_gamma[1]-0.7071067811865475*B2[0]*p0_over_gamma[1]; 
  alpha[5] = 1.224744871391589*B0[0]*p2_over_gamma[6]-1.224744871391589*B2[0]*p0_over_gamma[6]+0.7071067811865475*B0[0]*p2_over_gamma[3]-0.7071067811865475*B2[0]*p0_over_gamma[3]; 
  alpha[6] = 2.0*E1[4]-1.224744871391589*p0_over_gamma[2]*B2[4]-0.7071067811865475*p0_over_gamma[0]*B2[4]+1.224744871391589*p2_over_gamma[2]*B0[4]+0.7071067811865475*p2_over_gamma[0]*B0[4]; 
  alpha[7] = 2.0*E1[5]-1.224744871391589*p0_over_gamma[2]*B2[5]-0.7071067811865475*p0_over_gamma[0]*B2[5]+1.224744871391589*p2_over_gamma[2]*B0[5]+0.7071067811865475*p2_over_gamma[0]*B0[5]; 
  alpha[8] = 2.0*E1[6]-1.224744871391589*p0_over_gamma[2]*B2[6]-0.7071067811865475*p0_over_gamma[0]*B2[6]+1.224744871391589*p2_over_gamma[2]*B0[6]+0.7071067811865475*p2_over_gamma[0]*B0[6]; 
  alpha[9] = 1.224744871391589*B0[1]*p2_over_gamma[4]-1.224744871391589*B2[1]*p0_over_gamma[4]+0.7071067811865475*B0[1]*p2_over_gamma[1]-0.7071067811865475*B2[1]*p0_over_gamma[1]; 
  alpha[10] = 1.224744871391589*B0[2]*p2_over_gamma[4]-1.224744871391589*B2[2]*p0_over_gamma[4]-0.7071067811865475*p0_over_gamma[1]*B2[2]+0.7071067811865475*p2_over_gamma[1]*B0[2]; 
  alpha[11] = 1.224744871391589*B0[3]*p2_over_gamma[4]-1.224744871391589*B2[3]*p0_over_gamma[4]-0.7071067811865475*p0_over_gamma[1]*B2[3]+0.7071067811865475*p2_over_gamma[1]*B0[3]; 
  alpha[12] = 1.224744871391589*B0[1]*p2_over_gamma[6]-1.224744871391589*B2[1]*p0_over_gamma[6]+0.7071067811865475*B0[1]*p2_over_gamma[3]-0.7071067811865475*B2[1]*p0_over_gamma[3]; 
  alpha[13] = 1.224744871391589*B0[2]*p2_over_gamma[6]-1.224744871391589*B2[2]*p0_over_gamma[6]+0.7071067811865475*B0[2]*p2_over_gamma[3]-0.7071067811865475*B2[2]*p0_over_gamma[3]; 
  alpha[14] = 1.224744871391589*B0[3]*p2_over_gamma[6]-1.224744871391589*B2[3]*p0_over_gamma[6]+0.7071067811865475*B0[3]*p2_over_gamma[3]-0.7071067811865475*B2[3]*p0_over_gamma[3]; 
  alpha[15] = 1.224744871391589*B0[0]*p2_over_gamma[7]-1.224744871391589*B2[0]*p0_over_gamma[7]+0.7071067811865475*B0[0]*p2_over_gamma[5]-0.7071067811865475*B2[0]*p0_over_gamma[5]; 
  alpha[16] = 2.0*E1[7]-1.224744871391589*p0_over_gamma[2]*B2[7]-0.7071067811865475*p0_over_gamma[0]*B2[7]+1.224744871391589*p2_over_gamma[2]*B0[7]+0.7071067811865475*p2_over_gamma[0]*B0[7]; 
  alpha[17] = 1.224744871391589*B0[4]*p2_over_gamma[4]-1.224744871391589*B2[4]*p0_over_gamma[4]-0.7071067811865475*p0_over_gamma[1]*B2[4]+0.7071067811865475*p2_over_gamma[1]*B0[4]; 
  alpha[18] = (-1.224744871391589*p0_over_gamma[4]*B2[5])-0.7071067811865475*p0_over_gamma[1]*B2[5]+1.224744871391589*p2_over_gamma[4]*B0[5]+0.7071067811865475*p2_over_gamma[1]*B0[5]; 
  alpha[19] = (-1.224744871391589*p0_over_gamma[4]*B2[6])-0.7071067811865475*p0_over_gamma[1]*B2[6]+1.224744871391589*p2_over_gamma[4]*B0[6]+0.7071067811865475*p2_over_gamma[1]*B0[6]; 
  alpha[20] = 1.224744871391589*B0[4]*p2_over_gamma[6]-1.224744871391589*B2[4]*p0_over_gamma[6]-0.7071067811865475*p0_over_gamma[3]*B2[4]+0.7071067811865475*p2_over_gamma[3]*B0[4]; 
  alpha[21] = 1.224744871391589*B0[5]*p2_over_gamma[6]-1.224744871391589*B2[5]*p0_over_gamma[6]-0.7071067811865475*p0_over_gamma[3]*B2[5]+0.7071067811865475*p2_over_gamma[3]*B0[5]; 
  alpha[22] = 1.224744871391589*B0[6]*p2_over_gamma[6]-1.224744871391589*B2[6]*p0_over_gamma[6]-0.7071067811865475*p0_over_gamma[3]*B2[6]+0.7071067811865475*p2_over_gamma[3]*B0[6]; 
  alpha[23] = 1.224744871391589*B0[1]*p2_over_gamma[7]-1.224744871391589*B2[1]*p0_over_gamma[7]+0.7071067811865475*B0[1]*p2_over_gamma[5]-0.7071067811865475*B2[1]*p0_over_gamma[5]; 
  alpha[24] = 1.224744871391589*B0[2]*p2_over_gamma[7]-1.224744871391589*B2[2]*p0_over_gamma[7]+0.7071067811865475*B0[2]*p2_over_gamma[5]-0.7071067811865475*B2[2]*p0_over_gamma[5]; 
  alpha[25] = 1.224744871391589*B0[3]*p2_over_gamma[7]-1.224744871391589*B2[3]*p0_over_gamma[7]+0.7071067811865475*B0[3]*p2_over_gamma[5]-0.7071067811865475*B2[3]*p0_over_gamma[5]; 
  alpha[26] = (-1.224744871391589*p0_over_gamma[4]*B2[7])-0.7071067811865475*p0_over_gamma[1]*B2[7]+1.224744871391589*p2_over_gamma[4]*B0[7]+0.7071067811865475*p2_over_gamma[1]*B0[7]; 
  alpha[27] = (-1.224744871391589*p0_over_gamma[6]*B2[7])-0.7071067811865475*p0_over_gamma[3]*B2[7]+1.224744871391589*p2_over_gamma[6]*B0[7]+0.7071067811865475*p2_over_gamma[3]*B0[7]; 
  alpha[28] = 1.224744871391589*B0[4]*p2_over_gamma[7]-1.224744871391589*B2[4]*p0_over_gamma[7]+0.7071067811865475*B0[4]*p2_over_gamma[5]-0.7071067811865475*B2[4]*p0_over_gamma[5]; 
  alpha[29] = 1.224744871391589*B0[5]*p2_over_gamma[7]-1.224744871391589*B2[5]*p0_over_gamma[7]+0.7071067811865475*B0[5]*p2_over_gamma[5]-0.7071067811865475*B2[5]*p0_over_gamma[5]; 
  alpha[30] = 1.224744871391589*B0[6]*p2_over_gamma[7]-1.224744871391589*B2[6]*p0_over_gamma[7]-0.7071067811865475*p0_over_gamma[5]*B2[6]+0.7071067811865475*p2_over_gamma[5]*B0[6]; 
  alpha[31] = 1.224744871391589*B0[7]*p2_over_gamma[7]-1.224744871391589*B2[7]*p0_over_gamma[7]-0.7071067811865475*p0_over_gamma[5]*B2[7]+0.7071067811865475*p2_over_gamma[5]*B0[7]; 

  if ((-alpha[31])+alpha[30]+alpha[29]+alpha[28]+alpha[27]+alpha[26]-alpha[25]-alpha[24]-alpha[23]-alpha[22]-alpha[21]-alpha[20]-alpha[19]-alpha[18]-alpha[17]-alpha[16]+alpha[15]+alpha[14]+alpha[13]+alpha[12]+alpha[11]+alpha[10]+alpha[9]+alpha[8]+alpha[7]+alpha[6]-alpha[5]-alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[0] = ser_6x_p1_surfx5_eval_quad_node_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_6x_p1_surfx5_eval_quad_node_0_l(fEdge); 
  } 
  if (alpha[31]-alpha[30]-alpha[29]-alpha[28]-alpha[27]+alpha[26]+alpha[25]+alpha[24]+alpha[23]+alpha[22]+alpha[21]+alpha[20]-alpha[19]-alpha[18]-alpha[17]-alpha[16]-alpha[15]-alpha[14]-alpha[13]-alpha[12]+alpha[11]+alpha[10]+alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]-alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[1] = ser_6x_p1_surfx5_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = ser_6x_p1_surfx5_eval_quad_node_1_l(fEdge); 
  } 
  if (alpha[31]-alpha[30]-alpha[29]-alpha[28]+alpha[27]-alpha[26]+alpha[25]+alpha[24]+alpha[23]-alpha[22]-alpha[21]-alpha[20]+alpha[19]+alpha[18]+alpha[17]-alpha[16]-alpha[15]+alpha[14]+alpha[13]+alpha[12]-alpha[11]-alpha[10]-alpha[9]+alpha[8]+alpha[7]+alpha[6]-alpha[5]+alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[2] = ser_6x_p1_surfx5_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[2] = ser_6x_p1_surfx5_eval_quad_node_2_l(fEdge); 
  } 
  if ((-alpha[31])+alpha[30]+alpha[29]+alpha[28]-alpha[27]-alpha[26]-alpha[25]-alpha[24]-alpha[23]+alpha[22]+alpha[21]+alpha[20]+alpha[19]+alpha[18]+alpha[17]-alpha[16]+alpha[15]-alpha[14]-alpha[13]-alpha[12]-alpha[11]-alpha[10]-alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]+alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[3] = ser_6x_p1_surfx5_eval_quad_node_3_r(fSkin); 
  } else { 
    fUpwindQuad[3] = ser_6x_p1_surfx5_eval_quad_node_3_l(fEdge); 
  } 
  if (alpha[31]-alpha[30]-alpha[29]+alpha[28]-alpha[27]-alpha[26]+alpha[25]-alpha[24]-alpha[23]+alpha[22]+alpha[21]-alpha[20]+alpha[19]+alpha[18]-alpha[17]+alpha[16]+alpha[15]-alpha[14]+alpha[13]+alpha[12]-alpha[11]+alpha[10]+alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]-alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[4] = ser_6x_p1_surfx5_eval_quad_node_4_r(fSkin); 
  } else { 
    fUpwindQuad[4] = ser_6x_p1_surfx5_eval_quad_node_4_l(fEdge); 
  } 
  if ((-alpha[31])+alpha[30]+alpha[29]-alpha[28]+alpha[27]-alpha[26]-alpha[25]+alpha[24]+alpha[23]-alpha[22]-alpha[21]+alpha[20]+alpha[19]+alpha[18]-alpha[17]+alpha[16]-alpha[15]+alpha[14]-alpha[13]-alpha[12]-alpha[11]+alpha[10]+alpha[9]-alpha[8]-alpha[7]+alpha[6]+alpha[5]-alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[5] = ser_6x_p1_surfx5_eval_quad_node_5_r(fSkin); 
  } else { 
    fUpwindQuad[5] = ser_6x_p1_surfx5_eval_quad_node_5_l(fEdge); 
  } 
  if ((-alpha[31])+alpha[30]+alpha[29]-alpha[28]-alpha[27]+alpha[26]-alpha[25]+alpha[24]+alpha[23]+alpha[22]+alpha[21]-alpha[20]-alpha[19]-alpha[18]+alpha[17]+alpha[16]-alpha[15]-alpha[14]+alpha[13]+alpha[12]+alpha[11]-alpha[10]-alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]+alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[6] = ser_6x_p1_surfx5_eval_quad_node_6_r(fSkin); 
  } else { 
    fUpwindQuad[6] = ser_6x_p1_surfx5_eval_quad_node_6_l(fEdge); 
  } 
  if (alpha[31]-alpha[30]-alpha[29]+alpha[28]+alpha[27]+alpha[26]+alpha[25]-alpha[24]-alpha[23]-alpha[22]-alpha[21]+alpha[20]-alpha[19]-alpha[18]+alpha[17]+alpha[16]+alpha[15]+alpha[14]-alpha[13]-alpha[12]+alpha[11]-alpha[10]-alpha[9]-alpha[8]-alpha[7]+alpha[6]+alpha[5]+alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[7] = ser_6x_p1_surfx5_eval_quad_node_7_r(fSkin); 
  } else { 
    fUpwindQuad[7] = ser_6x_p1_surfx5_eval_quad_node_7_l(fEdge); 
  } 
  if (alpha[31]-alpha[30]+alpha[29]-alpha[28]-alpha[27]-alpha[26]-alpha[25]+alpha[24]-alpha[23]+alpha[22]-alpha[21]+alpha[20]+alpha[19]-alpha[18]+alpha[17]+alpha[16]+alpha[15]+alpha[14]-alpha[13]+alpha[12]+alpha[11]-alpha[10]+alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]-alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[8] = ser_6x_p1_surfx5_eval_quad_node_8_r(fSkin); 
  } else { 
    fUpwindQuad[8] = ser_6x_p1_surfx5_eval_quad_node_8_l(fEdge); 
  } 
  if ((-alpha[31])+alpha[30]-alpha[29]+alpha[28]+alpha[27]-alpha[26]+alpha[25]-alpha[24]+alpha[23]-alpha[22]+alpha[21]-alpha[20]+alpha[19]-alpha[18]+alpha[17]+alpha[16]-alpha[15]-alpha[14]+alpha[13]-alpha[12]+alpha[11]-alpha[10]+alpha[9]-alpha[8]+alpha[7]-alpha[6]+alpha[5]-alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[9] = ser_6x_p1_surfx5_eval_quad_node_9_r(fSkin); 
  } else { 
    fUpwindQuad[9] = ser_6x_p1_surfx5_eval_quad_node_9_l(fEdge); 
  } 
  if ((-alpha[31])+alpha[30]-alpha[29]+alpha[28]-alpha[27]+alpha[26]+alpha[25]-alpha[24]+alpha[23]+alpha[22]-alpha[21]+alpha[20]-alpha[19]+alpha[18]-alpha[17]+alpha[16]-alpha[15]+alpha[14]-alpha[13]+alpha[12]-alpha[11]+alpha[10]-alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]+alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[10] = ser_6x_p1_surfx5_eval_quad_node_10_r(fSkin); 
  } else { 
    fUpwindQuad[10] = ser_6x_p1_surfx5_eval_quad_node_10_l(fEdge); 
  } 
  if (alpha[31]-alpha[30]+alpha[29]-alpha[28]+alpha[27]+alpha[26]-alpha[25]+alpha[24]-alpha[23]-alpha[22]+alpha[21]-alpha[20]-alpha[19]+alpha[18]-alpha[17]+alpha[16]+alpha[15]-alpha[14]+alpha[13]-alpha[12]-alpha[11]+alpha[10]-alpha[9]-alpha[8]+alpha[7]-alpha[6]+alpha[5]+alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[11] = ser_6x_p1_surfx5_eval_quad_node_11_r(fSkin); 
  } else { 
    fUpwindQuad[11] = ser_6x_p1_surfx5_eval_quad_node_11_l(fEdge); 
  } 
  if ((-alpha[31])+alpha[30]-alpha[29]-alpha[28]+alpha[27]+alpha[26]+alpha[25]+alpha[24]-alpha[23]-alpha[22]+alpha[21]+alpha[20]-alpha[19]+alpha[18]+alpha[17]-alpha[16]+alpha[15]-alpha[14]-alpha[13]+alpha[12]-alpha[11]-alpha[10]+alpha[9]+alpha[8]-alpha[7]-alpha[6]-alpha[5]-alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[12] = ser_6x_p1_surfx5_eval_quad_node_12_r(fSkin); 
  } else { 
    fUpwindQuad[12] = ser_6x_p1_surfx5_eval_quad_node_12_l(fEdge); 
  } 
  if (alpha[31]-alpha[30]+alpha[29]+alpha[28]-alpha[27]+alpha[26]-alpha[25]-alpha[24]+alpha[23]+alpha[22]-alpha[21]-alpha[20]-alpha[19]+alpha[18]+alpha[17]-alpha[16]-alpha[15]+alpha[14]+alpha[13]-alpha[12]-alpha[11]-alpha[10]+alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]-alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[13] = ser_6x_p1_surfx5_eval_quad_node_13_r(fSkin); 
  } else { 
    fUpwindQuad[13] = ser_6x_p1_surfx5_eval_quad_node_13_l(fEdge); 
  } 
  if (alpha[31]-alpha[30]+alpha[29]+alpha[28]+alpha[27]-alpha[26]-alpha[25]-alpha[24]+alpha[23]-alpha[22]+alpha[21]+alpha[20]+alpha[19]-alpha[18]-alpha[17]-alpha[16]-alpha[15]-alpha[14]-alpha[13]+alpha[12]+alpha[11]+alpha[10]-alpha[9]+alpha[8]-alpha[7]-alpha[6]-alpha[5]+alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[14] = ser_6x_p1_surfx5_eval_quad_node_14_r(fSkin); 
  } else { 
    fUpwindQuad[14] = ser_6x_p1_surfx5_eval_quad_node_14_l(fEdge); 
  } 
  if ((-alpha[31])+alpha[30]-alpha[29]-alpha[28]-alpha[27]-alpha[26]+alpha[25]+alpha[24]-alpha[23]+alpha[22]-alpha[21]-alpha[20]+alpha[19]-alpha[18]-alpha[17]-alpha[16]+alpha[15]+alpha[14]+alpha[13]-alpha[12]+alpha[11]+alpha[10]-alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[15] = ser_6x_p1_surfx5_eval_quad_node_15_r(fSkin); 
  } else { 
    fUpwindQuad[15] = ser_6x_p1_surfx5_eval_quad_node_15_l(fEdge); 
  } 
  if (alpha[31]+alpha[30]-alpha[29]-alpha[28]-alpha[27]-alpha[26]-alpha[25]-alpha[24]+alpha[23]-alpha[22]+alpha[21]+alpha[20]-alpha[19]+alpha[18]+alpha[17]+alpha[16]+alpha[15]+alpha[14]+alpha[13]-alpha[12]+alpha[11]+alpha[10]-alpha[9]+alpha[8]-alpha[7]-alpha[6]-alpha[5]-alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[16] = ser_6x_p1_surfx5_eval_quad_node_16_r(fSkin); 
  } else { 
    fUpwindQuad[16] = ser_6x_p1_surfx5_eval_quad_node_16_l(fEdge); 
  } 
  if ((-alpha[31])-alpha[30]+alpha[29]+alpha[28]+alpha[27]-alpha[26]+alpha[25]+alpha[24]-alpha[23]+alpha[22]-alpha[21]-alpha[20]-alpha[19]+alpha[18]+alpha[17]+alpha[16]-alpha[15]-alpha[14]-alpha[13]+alpha[12]+alpha[11]+alpha[10]-alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]-alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[17] = ser_6x_p1_surfx5_eval_quad_node_17_r(fSkin); 
  } else { 
    fUpwindQuad[17] = ser_6x_p1_surfx5_eval_quad_node_17_l(fEdge); 
  } 
  if ((-alpha[31])-alpha[30]+alpha[29]+alpha[28]-alpha[27]+alpha[26]+alpha[25]+alpha[24]-alpha[23]-alpha[22]+alpha[21]+alpha[20]+alpha[19]-alpha[18]-alpha[17]+alpha[16]-alpha[15]+alpha[14]+alpha[13]-alpha[12]-alpha[11]-alpha[10]+alpha[9]+alpha[8]-alpha[7]-alpha[6]-alpha[5]+alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[18] = ser_6x_p1_surfx5_eval_quad_node_18_r(fSkin); 
  } else { 
    fUpwindQuad[18] = ser_6x_p1_surfx5_eval_quad_node_18_l(fEdge); 
  } 
  if (alpha[31]+alpha[30]-alpha[29]-alpha[28]+alpha[27]+alpha[26]-alpha[25]-alpha[24]+alpha[23]+alpha[22]-alpha[21]-alpha[20]+alpha[19]-alpha[18]-alpha[17]+alpha[16]+alpha[15]-alpha[14]-alpha[13]+alpha[12]-alpha[11]-alpha[10]+alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]+alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[19] = ser_6x_p1_surfx5_eval_quad_node_19_r(fSkin); 
  } else { 
    fUpwindQuad[19] = ser_6x_p1_surfx5_eval_quad_node_19_l(fEdge); 
  } 
  if ((-alpha[31])-alpha[30]+alpha[29]-alpha[28]+alpha[27]+alpha[26]+alpha[25]-alpha[24]+alpha[23]+alpha[22]-alpha[21]+alpha[20]+alpha[19]-alpha[18]+alpha[17]-alpha[16]+alpha[15]-alpha[14]+alpha[13]-alpha[12]-alpha[11]+alpha[10]-alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]-alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[20] = ser_6x_p1_surfx5_eval_quad_node_20_r(fSkin); 
  } else { 
    fUpwindQuad[20] = ser_6x_p1_surfx5_eval_quad_node_20_l(fEdge); 
  } 
  if (alpha[31]+alpha[30]-alpha[29]+alpha[28]-alpha[27]+alpha[26]-alpha[25]+alpha[24]-alpha[23]-alpha[22]+alpha[21]-alpha[20]+alpha[19]-alpha[18]+alpha[17]-alpha[16]-alpha[15]+alpha[14]-alpha[13]+alpha[12]-alpha[11]+alpha[10]-alpha[9]-alpha[8]+alpha[7]-alpha[6]+alpha[5]-alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[21] = ser_6x_p1_surfx5_eval_quad_node_21_r(fSkin); 
  } else { 
    fUpwindQuad[21] = ser_6x_p1_surfx5_eval_quad_node_21_l(fEdge); 
  } 
  if (alpha[31]+alpha[30]-alpha[29]+alpha[28]+alpha[27]-alpha[26]-alpha[25]+alpha[24]-alpha[23]+alpha[22]-alpha[21]+alpha[20]-alpha[19]+alpha[18]-alpha[17]-alpha[16]-alpha[15]-alpha[14]+alpha[13]-alpha[12]+alpha[11]-alpha[10]+alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]+alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[22] = ser_6x_p1_surfx5_eval_quad_node_22_r(fSkin); 
  } else { 
    fUpwindQuad[22] = ser_6x_p1_surfx5_eval_quad_node_22_l(fEdge); 
  } 
  if ((-alpha[31])-alpha[30]+alpha[29]-alpha[28]-alpha[27]-alpha[26]+alpha[25]-alpha[24]+alpha[23]-alpha[22]+alpha[21]-alpha[20]-alpha[19]+alpha[18]-alpha[17]-alpha[16]+alpha[15]+alpha[14]-alpha[13]+alpha[12]+alpha[11]-alpha[10]+alpha[9]-alpha[8]+alpha[7]-alpha[6]+alpha[5]+alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[23] = ser_6x_p1_surfx5_eval_quad_node_23_r(fSkin); 
  } else { 
    fUpwindQuad[23] = ser_6x_p1_surfx5_eval_quad_node_23_l(fEdge); 
  } 
  if ((-alpha[31])-alpha[30]-alpha[29]+alpha[28]+alpha[27]+alpha[26]-alpha[25]+alpha[24]+alpha[23]+alpha[22]+alpha[21]-alpha[20]+alpha[19]+alpha[18]-alpha[17]-alpha[16]+alpha[15]+alpha[14]-alpha[13]-alpha[12]+alpha[11]-alpha[10]-alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]-alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[24] = ser_6x_p1_surfx5_eval_quad_node_24_r(fSkin); 
  } else { 
    fUpwindQuad[24] = ser_6x_p1_surfx5_eval_quad_node_24_l(fEdge); 
  } 
  if (alpha[31]+alpha[30]+alpha[29]-alpha[28]-alpha[27]+alpha[26]+alpha[25]-alpha[24]-alpha[23]-alpha[22]-alpha[21]+alpha[20]+alpha[19]+alpha[18]-alpha[17]-alpha[16]-alpha[15]-alpha[14]+alpha[13]+alpha[12]+alpha[11]-alpha[10]-alpha[9]-alpha[8]-alpha[7]+alpha[6]+alpha[5]-alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[25] = ser_6x_p1_surfx5_eval_quad_node_25_r(fSkin); 
  } else { 
    fUpwindQuad[25] = ser_6x_p1_surfx5_eval_quad_node_25_l(fEdge); 
  } 
  if (alpha[31]+alpha[30]+alpha[29]-alpha[28]+alpha[27]-alpha[26]+alpha[25]-alpha[24]-alpha[23]+alpha[22]+alpha[21]-alpha[20]-alpha[19]-alpha[18]+alpha[17]-alpha[16]-alpha[15]+alpha[14]-alpha[13]-alpha[12]-alpha[11]+alpha[10]+alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]+alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[26] = ser_6x_p1_surfx5_eval_quad_node_26_r(fSkin); 
  } else { 
    fUpwindQuad[26] = ser_6x_p1_surfx5_eval_quad_node_26_l(fEdge); 
  } 
  if ((-alpha[31])-alpha[30]-alpha[29]+alpha[28]-alpha[27]-alpha[26]-alpha[25]+alpha[24]+alpha[23]-alpha[22]-alpha[21]+alpha[20]-alpha[19]-alpha[18]+alpha[17]-alpha[16]+alpha[15]-alpha[14]+alpha[13]+alpha[12]-alpha[11]+alpha[10]+alpha[9]-alpha[8]-alpha[7]+alpha[6]+alpha[5]+alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[27] = ser_6x_p1_surfx5_eval_quad_node_27_r(fSkin); 
  } else { 
    fUpwindQuad[27] = ser_6x_p1_surfx5_eval_quad_node_27_l(fEdge); 
  } 
  if (alpha[31]+alpha[30]+alpha[29]+alpha[28]-alpha[27]-alpha[26]+alpha[25]+alpha[24]+alpha[23]-alpha[22]-alpha[21]-alpha[20]-alpha[19]-alpha[18]-alpha[17]+alpha[16]+alpha[15]-alpha[14]-alpha[13]-alpha[12]-alpha[11]-alpha[10]-alpha[9]+alpha[8]+alpha[7]+alpha[6]-alpha[5]-alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[28] = ser_6x_p1_surfx5_eval_quad_node_28_r(fSkin); 
  } else { 
    fUpwindQuad[28] = ser_6x_p1_surfx5_eval_quad_node_28_l(fEdge); 
  } 
  if ((-alpha[31])-alpha[30]-alpha[29]-alpha[28]+alpha[27]-alpha[26]-alpha[25]-alpha[24]-alpha[23]+alpha[22]+alpha[21]+alpha[20]-alpha[19]-alpha[18]-alpha[17]+alpha[16]-alpha[15]+alpha[14]+alpha[13]+alpha[12]-alpha[11]-alpha[10]-alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]-alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[29] = ser_6x_p1_surfx5_eval_quad_node_29_r(fSkin); 
  } else { 
    fUpwindQuad[29] = ser_6x_p1_surfx5_eval_quad_node_29_l(fEdge); 
  } 
  if ((-alpha[31])-alpha[30]-alpha[29]-alpha[28]-alpha[27]+alpha[26]-alpha[25]-alpha[24]-alpha[23]-alpha[22]-alpha[21]-alpha[20]+alpha[19]+alpha[18]+alpha[17]+alpha[16]-alpha[15]-alpha[14]-alpha[13]-alpha[12]+alpha[11]+alpha[10]+alpha[9]+alpha[8]+alpha[7]+alpha[6]-alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[30] = ser_6x_p1_surfx5_eval_quad_node_30_r(fSkin); 
  } else { 
    fUpwindQuad[30] = ser_6x_p1_surfx5_eval_quad_node_30_l(fEdge); 
  } 
  if (alpha[31]+alpha[30]+alpha[29]+alpha[28]+alpha[27]+alpha[26]+alpha[25]+alpha[24]+alpha[23]+alpha[22]+alpha[21]+alpha[20]+alpha[19]+alpha[18]+alpha[17]+alpha[16]+alpha[15]+alpha[14]+alpha[13]+alpha[12]+alpha[11]+alpha[10]+alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[31] = ser_6x_p1_surfx5_eval_quad_node_31_r(fSkin); 
  } else { 
    fUpwindQuad[31] = ser_6x_p1_surfx5_eval_quad_node_31_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_6x_p1_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.1767766952966368*(alpha[31]*fUpwind[31]+alpha[30]*fUpwind[30]+alpha[29]*fUpwind[29]+alpha[28]*fUpwind[28]+alpha[27]*fUpwind[27]+alpha[26]*fUpwind[26]+alpha[25]*fUpwind[25]+alpha[24]*fUpwind[24]+alpha[23]*fUpwind[23]+alpha[22]*fUpwind[22]+alpha[21]*fUpwind[21]+alpha[20]*fUpwind[20]+alpha[19]*fUpwind[19]+alpha[18]*fUpwind[18]+alpha[17]*fUpwind[17]+alpha[16]*fUpwind[16]+alpha[15]*fUpwind[15]+alpha[14]*fUpwind[14]+alpha[13]*fUpwind[13]+alpha[12]*fUpwind[12]+alpha[11]*fUpwind[11]+alpha[10]*fUpwind[10]+alpha[9]*fUpwind[9]+alpha[8]*fUpwind[8]+alpha[7]*fUpwind[7]+alpha[6]*fUpwind[6]+alpha[5]*fUpwind[5]+alpha[4]*fUpwind[4]+alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] = 0.1767766952966368*(alpha[30]*fUpwind[31]+fUpwind[30]*alpha[31]+alpha[25]*fUpwind[29]+fUpwind[25]*alpha[29]+alpha[24]*fUpwind[28]+fUpwind[24]*alpha[28]+alpha[22]*fUpwind[27]+fUpwind[22]*alpha[27]+alpha[19]*fUpwind[26]+fUpwind[19]*alpha[26]+alpha[15]*fUpwind[23]+fUpwind[15]*alpha[23]+alpha[14]*fUpwind[21]+fUpwind[14]*alpha[21]+alpha[13]*fUpwind[20]+fUpwind[13]*alpha[20]+alpha[11]*fUpwind[18]+fUpwind[11]*alpha[18]+alpha[10]*fUpwind[17]+fUpwind[10]*alpha[17]+alpha[8]*fUpwind[16]+fUpwind[8]*alpha[16]+alpha[5]*fUpwind[12]+fUpwind[5]*alpha[12]+alpha[4]*fUpwind[9]+fUpwind[4]*alpha[9]+alpha[3]*fUpwind[7]+fUpwind[3]*alpha[7]+alpha[2]*fUpwind[6]+fUpwind[2]*alpha[6]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] = 0.1767766952966368*(alpha[29]*fUpwind[31]+fUpwind[29]*alpha[31]+alpha[25]*fUpwind[30]+fUpwind[25]*alpha[30]+alpha[23]*fUpwind[28]+fUpwind[23]*alpha[28]+alpha[21]*fUpwind[27]+fUpwind[21]*alpha[27]+alpha[18]*fUpwind[26]+fUpwind[18]*alpha[26]+alpha[15]*fUpwind[24]+fUpwind[15]*alpha[24]+alpha[14]*fUpwind[22]+fUpwind[14]*alpha[22]+alpha[12]*fUpwind[20]+fUpwind[12]*alpha[20]+alpha[11]*fUpwind[19]+fUpwind[11]*alpha[19]+alpha[9]*fUpwind[17]+fUpwind[9]*alpha[17]+alpha[7]*fUpwind[16]+fUpwind[7]*alpha[16]+alpha[5]*fUpwind[13]+fUpwind[5]*alpha[13]+alpha[4]*fUpwind[10]+fUpwind[4]*alpha[10]+alpha[3]*fUpwind[8]+fUpwind[3]*alpha[8]+alpha[1]*fUpwind[6]+fUpwind[1]*alpha[6]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] = 0.1767766952966368*(alpha[28]*fUpwind[31]+fUpwind[28]*alpha[31]+alpha[24]*fUpwind[30]+fUpwind[24]*alpha[30]+alpha[23]*fUpwind[29]+fUpwind[23]*alpha[29]+alpha[20]*fUpwind[27]+fUpwind[20]*alpha[27]+alpha[17]*fUpwind[26]+fUpwind[17]*alpha[26]+alpha[15]*fUpwind[25]+fUpwind[15]*alpha[25]+alpha[13]*fUpwind[22]+fUpwind[13]*alpha[22]+alpha[12]*fUpwind[21]+fUpwind[12]*alpha[21]+alpha[10]*fUpwind[19]+fUpwind[10]*alpha[19]+alpha[9]*fUpwind[18]+fUpwind[9]*alpha[18]+alpha[6]*fUpwind[16]+fUpwind[6]*alpha[16]+alpha[5]*fUpwind[14]+fUpwind[5]*alpha[14]+alpha[4]*fUpwind[11]+fUpwind[4]*alpha[11]+alpha[2]*fUpwind[8]+fUpwind[2]*alpha[8]+alpha[1]*fUpwind[7]+fUpwind[1]*alpha[7]+alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]); 
  Ghat[4] = 0.1767766952966368*(alpha[27]*fUpwind[31]+fUpwind[27]*alpha[31]+alpha[22]*fUpwind[30]+fUpwind[22]*alpha[30]+alpha[21]*fUpwind[29]+fUpwind[21]*alpha[29]+alpha[20]*fUpwind[28]+fUpwind[20]*alpha[28]+alpha[16]*fUpwind[26]+fUpwind[16]*alpha[26]+alpha[14]*fUpwind[25]+fUpwind[14]*alpha[25]+alpha[13]*fUpwind[24]+fUpwind[13]*alpha[24]+alpha[12]*fUpwind[23]+fUpwind[12]*alpha[23]+alpha[8]*fUpwind[19]+fUpwind[8]*alpha[19]+alpha[7]*fUpwind[18]+fUpwind[7]*alpha[18]+alpha[6]*fUpwind[17]+fUpwind[6]*alpha[17]+alpha[5]*fUpwind[15]+fUpwind[5]*alpha[15]+alpha[3]*fUpwind[11]+fUpwind[3]*alpha[11]+alpha[2]*fUpwind[10]+fUpwind[2]*alpha[10]+alpha[1]*fUpwind[9]+fUpwind[1]*alpha[9]+alpha[0]*fUpwind[4]+fUpwind[0]*alpha[4]); 
  Ghat[5] = 0.1767766952966368*(alpha[26]*fUpwind[31]+fUpwind[26]*alpha[31]+alpha[19]*fUpwind[30]+fUpwind[19]*alpha[30]+alpha[18]*fUpwind[29]+fUpwind[18]*alpha[29]+alpha[17]*fUpwind[28]+fUpwind[17]*alpha[28]+alpha[16]*fUpwind[27]+fUpwind[16]*alpha[27]+alpha[11]*fUpwind[25]+fUpwind[11]*alpha[25]+alpha[10]*fUpwind[24]+fUpwind[10]*alpha[24]+alpha[9]*fUpwind[23]+fUpwind[9]*alpha[23]+alpha[8]*fUpwind[22]+fUpwind[8]*alpha[22]+alpha[7]*fUpwind[21]+fUpwind[7]*alpha[21]+alpha[6]*fUpwind[20]+fUpwind[6]*alpha[20]+alpha[4]*fUpwind[15]+fUpwind[4]*alpha[15]+alpha[3]*fUpwind[14]+fUpwind[3]*alpha[14]+alpha[2]*fUpwind[13]+fUpwind[2]*alpha[13]+alpha[1]*fUpwind[12]+fUpwind[1]*alpha[12]+alpha[0]*fUpwind[5]+fUpwind[0]*alpha[5]); 
  Ghat[6] = 0.1767766952966368*(alpha[25]*fUpwind[31]+fUpwind[25]*alpha[31]+alpha[29]*fUpwind[30]+fUpwind[29]*alpha[30]+alpha[15]*fUpwind[28]+fUpwind[15]*alpha[28]+alpha[14]*fUpwind[27]+fUpwind[14]*alpha[27]+alpha[11]*fUpwind[26]+fUpwind[11]*alpha[26]+alpha[23]*fUpwind[24]+fUpwind[23]*alpha[24]+alpha[21]*fUpwind[22]+fUpwind[21]*alpha[22]+alpha[5]*fUpwind[20]+fUpwind[5]*alpha[20]+alpha[18]*fUpwind[19]+fUpwind[18]*alpha[19]+alpha[4]*fUpwind[17]+fUpwind[4]*alpha[17]+alpha[3]*fUpwind[16]+fUpwind[3]*alpha[16]+alpha[12]*fUpwind[13]+fUpwind[12]*alpha[13]+alpha[9]*fUpwind[10]+fUpwind[9]*alpha[10]+alpha[7]*fUpwind[8]+fUpwind[7]*alpha[8]+alpha[0]*fUpwind[6]+fUpwind[0]*alpha[6]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[7] = 0.1767766952966368*(alpha[24]*fUpwind[31]+fUpwind[24]*alpha[31]+alpha[28]*fUpwind[30]+fUpwind[28]*alpha[30]+alpha[15]*fUpwind[29]+fUpwind[15]*alpha[29]+alpha[13]*fUpwind[27]+fUpwind[13]*alpha[27]+alpha[10]*fUpwind[26]+fUpwind[10]*alpha[26]+alpha[23]*fUpwind[25]+fUpwind[23]*alpha[25]+alpha[20]*fUpwind[22]+fUpwind[20]*alpha[22]+alpha[5]*fUpwind[21]+fUpwind[5]*alpha[21]+alpha[17]*fUpwind[19]+fUpwind[17]*alpha[19]+alpha[4]*fUpwind[18]+fUpwind[4]*alpha[18]+alpha[2]*fUpwind[16]+fUpwind[2]*alpha[16]+alpha[12]*fUpwind[14]+fUpwind[12]*alpha[14]+alpha[9]*fUpwind[11]+fUpwind[9]*alpha[11]+alpha[6]*fUpwind[8]+fUpwind[6]*alpha[8]+alpha[0]*fUpwind[7]+fUpwind[0]*alpha[7]+alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]); 
  Ghat[8] = 0.1767766952966368*(alpha[23]*fUpwind[31]+fUpwind[23]*alpha[31]+alpha[15]*fUpwind[30]+fUpwind[15]*alpha[30]+alpha[28]*fUpwind[29]+fUpwind[28]*alpha[29]+alpha[12]*fUpwind[27]+fUpwind[12]*alpha[27]+alpha[9]*fUpwind[26]+fUpwind[9]*alpha[26]+alpha[24]*fUpwind[25]+fUpwind[24]*alpha[25]+alpha[5]*fUpwind[22]+fUpwind[5]*alpha[22]+alpha[20]*fUpwind[21]+fUpwind[20]*alpha[21]+alpha[4]*fUpwind[19]+fUpwind[4]*alpha[19]+alpha[17]*fUpwind[18]+fUpwind[17]*alpha[18]+alpha[1]*fUpwind[16]+fUpwind[1]*alpha[16]+alpha[13]*fUpwind[14]+fUpwind[13]*alpha[14]+alpha[10]*fUpwind[11]+fUpwind[10]*alpha[11]+alpha[0]*fUpwind[8]+fUpwind[0]*alpha[8]+alpha[6]*fUpwind[7]+fUpwind[6]*alpha[7]+alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]); 
  Ghat[9] = 0.1767766952966368*(alpha[22]*fUpwind[31]+fUpwind[22]*alpha[31]+alpha[27]*fUpwind[30]+fUpwind[27]*alpha[30]+alpha[14]*fUpwind[29]+fUpwind[14]*alpha[29]+alpha[13]*fUpwind[28]+fUpwind[13]*alpha[28]+alpha[8]*fUpwind[26]+fUpwind[8]*alpha[26]+alpha[21]*fUpwind[25]+fUpwind[21]*alpha[25]+alpha[20]*fUpwind[24]+fUpwind[20]*alpha[24]+alpha[5]*fUpwind[23]+fUpwind[5]*alpha[23]+alpha[16]*fUpwind[19]+fUpwind[16]*alpha[19]+alpha[3]*fUpwind[18]+fUpwind[3]*alpha[18]+alpha[2]*fUpwind[17]+fUpwind[2]*alpha[17]+alpha[12]*fUpwind[15]+fUpwind[12]*alpha[15]+alpha[7]*fUpwind[11]+fUpwind[7]*alpha[11]+alpha[6]*fUpwind[10]+fUpwind[6]*alpha[10]+alpha[0]*fUpwind[9]+fUpwind[0]*alpha[9]+alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]); 
  Ghat[10] = 0.1767766952966368*(alpha[21]*fUpwind[31]+fUpwind[21]*alpha[31]+alpha[14]*fUpwind[30]+fUpwind[14]*alpha[30]+alpha[27]*fUpwind[29]+fUpwind[27]*alpha[29]+alpha[12]*fUpwind[28]+fUpwind[12]*alpha[28]+alpha[7]*fUpwind[26]+fUpwind[7]*alpha[26]+alpha[22]*fUpwind[25]+fUpwind[22]*alpha[25]+alpha[5]*fUpwind[24]+fUpwind[5]*alpha[24]+alpha[20]*fUpwind[23]+fUpwind[20]*alpha[23]+alpha[3]*fUpwind[19]+fUpwind[3]*alpha[19]+alpha[16]*fUpwind[18]+fUpwind[16]*alpha[18]+alpha[1]*fUpwind[17]+fUpwind[1]*alpha[17]+alpha[13]*fUpwind[15]+fUpwind[13]*alpha[15]+alpha[8]*fUpwind[11]+fUpwind[8]*alpha[11]+alpha[0]*fUpwind[10]+fUpwind[0]*alpha[10]+alpha[6]*fUpwind[9]+fUpwind[6]*alpha[9]+alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]); 
  Ghat[11] = 0.1767766952966368*(alpha[20]*fUpwind[31]+fUpwind[20]*alpha[31]+alpha[13]*fUpwind[30]+fUpwind[13]*alpha[30]+alpha[12]*fUpwind[29]+fUpwind[12]*alpha[29]+alpha[27]*fUpwind[28]+fUpwind[27]*alpha[28]+alpha[6]*fUpwind[26]+fUpwind[6]*alpha[26]+alpha[5]*fUpwind[25]+fUpwind[5]*alpha[25]+alpha[22]*fUpwind[24]+fUpwind[22]*alpha[24]+alpha[21]*fUpwind[23]+fUpwind[21]*alpha[23]+alpha[2]*fUpwind[19]+fUpwind[2]*alpha[19]+alpha[1]*fUpwind[18]+fUpwind[1]*alpha[18]+alpha[16]*fUpwind[17]+fUpwind[16]*alpha[17]+alpha[14]*fUpwind[15]+fUpwind[14]*alpha[15]+alpha[0]*fUpwind[11]+fUpwind[0]*alpha[11]+alpha[8]*fUpwind[10]+fUpwind[8]*alpha[10]+alpha[7]*fUpwind[9]+fUpwind[7]*alpha[9]+alpha[3]*fUpwind[4]+fUpwind[3]*alpha[4]); 
  Ghat[12] = 0.1767766952966368*(alpha[19]*fUpwind[31]+fUpwind[19]*alpha[31]+alpha[26]*fUpwind[30]+fUpwind[26]*alpha[30]+alpha[11]*fUpwind[29]+fUpwind[11]*alpha[29]+alpha[10]*fUpwind[28]+fUpwind[10]*alpha[28]+alpha[8]*fUpwind[27]+fUpwind[8]*alpha[27]+alpha[18]*fUpwind[25]+fUpwind[18]*alpha[25]+alpha[17]*fUpwind[24]+fUpwind[17]*alpha[24]+alpha[4]*fUpwind[23]+fUpwind[4]*alpha[23]+alpha[16]*fUpwind[22]+fUpwind[16]*alpha[22]+alpha[3]*fUpwind[21]+fUpwind[3]*alpha[21]+alpha[2]*fUpwind[20]+fUpwind[2]*alpha[20]+alpha[9]*fUpwind[15]+fUpwind[9]*alpha[15]+alpha[7]*fUpwind[14]+fUpwind[7]*alpha[14]+alpha[6]*fUpwind[13]+fUpwind[6]*alpha[13]+alpha[0]*fUpwind[12]+fUpwind[0]*alpha[12]+alpha[1]*fUpwind[5]+fUpwind[1]*alpha[5]); 
  Ghat[13] = 0.1767766952966368*(alpha[18]*fUpwind[31]+fUpwind[18]*alpha[31]+alpha[11]*fUpwind[30]+fUpwind[11]*alpha[30]+alpha[26]*fUpwind[29]+fUpwind[26]*alpha[29]+alpha[9]*fUpwind[28]+fUpwind[9]*alpha[28]+alpha[7]*fUpwind[27]+fUpwind[7]*alpha[27]+alpha[19]*fUpwind[25]+fUpwind[19]*alpha[25]+alpha[4]*fUpwind[24]+fUpwind[4]*alpha[24]+alpha[17]*fUpwind[23]+fUpwind[17]*alpha[23]+alpha[3]*fUpwind[22]+fUpwind[3]*alpha[22]+alpha[16]*fUpwind[21]+fUpwind[16]*alpha[21]+alpha[1]*fUpwind[20]+fUpwind[1]*alpha[20]+alpha[10]*fUpwind[15]+fUpwind[10]*alpha[15]+alpha[8]*fUpwind[14]+fUpwind[8]*alpha[14]+alpha[0]*fUpwind[13]+fUpwind[0]*alpha[13]+alpha[6]*fUpwind[12]+fUpwind[6]*alpha[12]+alpha[2]*fUpwind[5]+fUpwind[2]*alpha[5]); 
  Ghat[14] = 0.1767766952966368*(alpha[17]*fUpwind[31]+fUpwind[17]*alpha[31]+alpha[10]*fUpwind[30]+fUpwind[10]*alpha[30]+alpha[9]*fUpwind[29]+fUpwind[9]*alpha[29]+alpha[26]*fUpwind[28]+fUpwind[26]*alpha[28]+alpha[6]*fUpwind[27]+fUpwind[6]*alpha[27]+alpha[4]*fUpwind[25]+fUpwind[4]*alpha[25]+alpha[19]*fUpwind[24]+fUpwind[19]*alpha[24]+alpha[18]*fUpwind[23]+fUpwind[18]*alpha[23]+alpha[2]*fUpwind[22]+fUpwind[2]*alpha[22]+alpha[1]*fUpwind[21]+fUpwind[1]*alpha[21]+alpha[16]*fUpwind[20]+fUpwind[16]*alpha[20]+alpha[11]*fUpwind[15]+fUpwind[11]*alpha[15]+alpha[0]*fUpwind[14]+fUpwind[0]*alpha[14]+alpha[8]*fUpwind[13]+fUpwind[8]*alpha[13]+alpha[7]*fUpwind[12]+fUpwind[7]*alpha[12]+alpha[3]*fUpwind[5]+fUpwind[3]*alpha[5]); 
  Ghat[15] = 0.1767766952966368*(alpha[16]*fUpwind[31]+fUpwind[16]*alpha[31]+alpha[8]*fUpwind[30]+fUpwind[8]*alpha[30]+alpha[7]*fUpwind[29]+fUpwind[7]*alpha[29]+alpha[6]*fUpwind[28]+fUpwind[6]*alpha[28]+alpha[26]*fUpwind[27]+fUpwind[26]*alpha[27]+alpha[3]*fUpwind[25]+fUpwind[3]*alpha[25]+alpha[2]*fUpwind[24]+fUpwind[2]*alpha[24]+alpha[1]*fUpwind[23]+fUpwind[1]*alpha[23]+alpha[19]*fUpwind[22]+fUpwind[19]*alpha[22]+alpha[18]*fUpwind[21]+fUpwind[18]*alpha[21]+alpha[17]*fUpwind[20]+fUpwind[17]*alpha[20]+alpha[0]*fUpwind[15]+fUpwind[0]*alpha[15]+alpha[11]*fUpwind[14]+fUpwind[11]*alpha[14]+alpha[10]*fUpwind[13]+fUpwind[10]*alpha[13]+alpha[9]*fUpwind[12]+fUpwind[9]*alpha[12]+alpha[4]*fUpwind[5]+fUpwind[4]*alpha[5]); 
  Ghat[16] = 0.1767766952966368*(alpha[15]*fUpwind[31]+fUpwind[15]*alpha[31]+alpha[23]*fUpwind[30]+fUpwind[23]*alpha[30]+alpha[24]*fUpwind[29]+fUpwind[24]*alpha[29]+alpha[25]*fUpwind[28]+fUpwind[25]*alpha[28]+alpha[5]*fUpwind[27]+fUpwind[5]*alpha[27]+alpha[4]*fUpwind[26]+fUpwind[4]*alpha[26]+alpha[12]*fUpwind[22]+fUpwind[12]*alpha[22]+alpha[13]*fUpwind[21]+fUpwind[13]*alpha[21]+alpha[14]*fUpwind[20]+fUpwind[14]*alpha[20]+alpha[9]*fUpwind[19]+fUpwind[9]*alpha[19]+alpha[10]*fUpwind[18]+fUpwind[10]*alpha[18]+alpha[11]*fUpwind[17]+fUpwind[11]*alpha[17]+alpha[0]*fUpwind[16]+fUpwind[0]*alpha[16]+alpha[1]*fUpwind[8]+fUpwind[1]*alpha[8]+alpha[2]*fUpwind[7]+fUpwind[2]*alpha[7]+alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6]); 
  Ghat[17] = 0.1767766952966368*(alpha[14]*fUpwind[31]+fUpwind[14]*alpha[31]+alpha[21]*fUpwind[30]+fUpwind[21]*alpha[30]+alpha[22]*fUpwind[29]+fUpwind[22]*alpha[29]+alpha[5]*fUpwind[28]+fUpwind[5]*alpha[28]+alpha[25]*fUpwind[27]+fUpwind[25]*alpha[27]+alpha[3]*fUpwind[26]+fUpwind[3]*alpha[26]+alpha[12]*fUpwind[24]+fUpwind[12]*alpha[24]+alpha[13]*fUpwind[23]+fUpwind[13]*alpha[23]+alpha[15]*fUpwind[20]+fUpwind[15]*alpha[20]+alpha[7]*fUpwind[19]+fUpwind[7]*alpha[19]+alpha[8]*fUpwind[18]+fUpwind[8]*alpha[18]+alpha[0]*fUpwind[17]+fUpwind[0]*alpha[17]+alpha[11]*fUpwind[16]+fUpwind[11]*alpha[16]+alpha[1]*fUpwind[10]+fUpwind[1]*alpha[10]+alpha[2]*fUpwind[9]+fUpwind[2]*alpha[9]+alpha[4]*fUpwind[6]+fUpwind[4]*alpha[6]); 
  Ghat[18] = 0.1767766952966368*(alpha[13]*fUpwind[31]+fUpwind[13]*alpha[31]+alpha[20]*fUpwind[30]+fUpwind[20]*alpha[30]+alpha[5]*fUpwind[29]+fUpwind[5]*alpha[29]+alpha[22]*fUpwind[28]+fUpwind[22]*alpha[28]+alpha[24]*fUpwind[27]+fUpwind[24]*alpha[27]+alpha[2]*fUpwind[26]+fUpwind[2]*alpha[26]+alpha[12]*fUpwind[25]+fUpwind[12]*alpha[25]+alpha[14]*fUpwind[23]+fUpwind[14]*alpha[23]+alpha[15]*fUpwind[21]+fUpwind[15]*alpha[21]+alpha[6]*fUpwind[19]+fUpwind[6]*alpha[19]+alpha[0]*fUpwind[18]+fUpwind[0]*alpha[18]+alpha[8]*fUpwind[17]+fUpwind[8]*alpha[17]+alpha[10]*fUpwind[16]+fUpwind[10]*alpha[16]+alpha[1]*fUpwind[11]+fUpwind[1]*alpha[11]+alpha[3]*fUpwind[9]+fUpwind[3]*alpha[9]+alpha[4]*fUpwind[7]+fUpwind[4]*alpha[7]); 
  Ghat[19] = 0.1767766952966368*(alpha[12]*fUpwind[31]+fUpwind[12]*alpha[31]+alpha[5]*fUpwind[30]+fUpwind[5]*alpha[30]+alpha[20]*fUpwind[29]+fUpwind[20]*alpha[29]+alpha[21]*fUpwind[28]+fUpwind[21]*alpha[28]+alpha[23]*fUpwind[27]+fUpwind[23]*alpha[27]+alpha[1]*fUpwind[26]+fUpwind[1]*alpha[26]+alpha[13]*fUpwind[25]+fUpwind[13]*alpha[25]+alpha[14]*fUpwind[24]+fUpwind[14]*alpha[24]+alpha[15]*fUpwind[22]+fUpwind[15]*alpha[22]+alpha[0]*fUpwind[19]+fUpwind[0]*alpha[19]+alpha[6]*fUpwind[18]+fUpwind[6]*alpha[18]+alpha[7]*fUpwind[17]+fUpwind[7]*alpha[17]+alpha[9]*fUpwind[16]+fUpwind[9]*alpha[16]+alpha[2]*fUpwind[11]+fUpwind[2]*alpha[11]+alpha[3]*fUpwind[10]+fUpwind[3]*alpha[10]+alpha[4]*fUpwind[8]+fUpwind[4]*alpha[8]); 
  Ghat[20] = 0.1767766952966368*(alpha[11]*fUpwind[31]+fUpwind[11]*alpha[31]+alpha[18]*fUpwind[30]+fUpwind[18]*alpha[30]+alpha[19]*fUpwind[29]+fUpwind[19]*alpha[29]+alpha[4]*fUpwind[28]+fUpwind[4]*alpha[28]+alpha[3]*fUpwind[27]+fUpwind[3]*alpha[27]+alpha[25]*fUpwind[26]+fUpwind[25]*alpha[26]+alpha[9]*fUpwind[24]+fUpwind[9]*alpha[24]+alpha[10]*fUpwind[23]+fUpwind[10]*alpha[23]+alpha[7]*fUpwind[22]+fUpwind[7]*alpha[22]+alpha[8]*fUpwind[21]+fUpwind[8]*alpha[21]+alpha[0]*fUpwind[20]+fUpwind[0]*alpha[20]+alpha[15]*fUpwind[17]+fUpwind[15]*alpha[17]+alpha[14]*fUpwind[16]+fUpwind[14]*alpha[16]+alpha[1]*fUpwind[13]+fUpwind[1]*alpha[13]+alpha[2]*fUpwind[12]+fUpwind[2]*alpha[12]+alpha[5]*fUpwind[6]+fUpwind[5]*alpha[6]); 
  Ghat[21] = 0.1767766952966368*(alpha[10]*fUpwind[31]+fUpwind[10]*alpha[31]+alpha[17]*fUpwind[30]+fUpwind[17]*alpha[30]+alpha[4]*fUpwind[29]+fUpwind[4]*alpha[29]+alpha[19]*fUpwind[28]+fUpwind[19]*alpha[28]+alpha[2]*fUpwind[27]+fUpwind[2]*alpha[27]+alpha[24]*fUpwind[26]+fUpwind[24]*alpha[26]+alpha[9]*fUpwind[25]+fUpwind[9]*alpha[25]+alpha[11]*fUpwind[23]+fUpwind[11]*alpha[23]+alpha[6]*fUpwind[22]+fUpwind[6]*alpha[22]+alpha[0]*fUpwind[21]+fUpwind[0]*alpha[21]+alpha[8]*fUpwind[20]+fUpwind[8]*alpha[20]+alpha[15]*fUpwind[18]+fUpwind[15]*alpha[18]+alpha[13]*fUpwind[16]+fUpwind[13]*alpha[16]+alpha[1]*fUpwind[14]+fUpwind[1]*alpha[14]+alpha[3]*fUpwind[12]+fUpwind[3]*alpha[12]+alpha[5]*fUpwind[7]+fUpwind[5]*alpha[7]); 
  Ghat[22] = 0.1767766952966368*(alpha[9]*fUpwind[31]+fUpwind[9]*alpha[31]+alpha[4]*fUpwind[30]+fUpwind[4]*alpha[30]+alpha[17]*fUpwind[29]+fUpwind[17]*alpha[29]+alpha[18]*fUpwind[28]+fUpwind[18]*alpha[28]+alpha[1]*fUpwind[27]+fUpwind[1]*alpha[27]+alpha[23]*fUpwind[26]+fUpwind[23]*alpha[26]+alpha[10]*fUpwind[25]+fUpwind[10]*alpha[25]+alpha[11]*fUpwind[24]+fUpwind[11]*alpha[24]+alpha[0]*fUpwind[22]+fUpwind[0]*alpha[22]+alpha[6]*fUpwind[21]+fUpwind[6]*alpha[21]+alpha[7]*fUpwind[20]+fUpwind[7]*alpha[20]+alpha[15]*fUpwind[19]+fUpwind[15]*alpha[19]+alpha[12]*fUpwind[16]+fUpwind[12]*alpha[16]+alpha[2]*fUpwind[14]+fUpwind[2]*alpha[14]+alpha[3]*fUpwind[13]+fUpwind[3]*alpha[13]+alpha[5]*fUpwind[8]+fUpwind[5]*alpha[8]); 
  Ghat[23] = 0.1767766952966368*(alpha[8]*fUpwind[31]+fUpwind[8]*alpha[31]+alpha[16]*fUpwind[30]+fUpwind[16]*alpha[30]+alpha[3]*fUpwind[29]+fUpwind[3]*alpha[29]+alpha[2]*fUpwind[28]+fUpwind[2]*alpha[28]+alpha[19]*fUpwind[27]+fUpwind[19]*alpha[27]+alpha[22]*fUpwind[26]+fUpwind[22]*alpha[26]+alpha[7]*fUpwind[25]+fUpwind[7]*alpha[25]+alpha[6]*fUpwind[24]+fUpwind[6]*alpha[24]+alpha[0]*fUpwind[23]+fUpwind[0]*alpha[23]+alpha[11]*fUpwind[21]+fUpwind[11]*alpha[21]+alpha[10]*fUpwind[20]+fUpwind[10]*alpha[20]+alpha[14]*fUpwind[18]+fUpwind[14]*alpha[18]+alpha[13]*fUpwind[17]+fUpwind[13]*alpha[17]+alpha[1]*fUpwind[15]+fUpwind[1]*alpha[15]+alpha[4]*fUpwind[12]+fUpwind[4]*alpha[12]+alpha[5]*fUpwind[9]+fUpwind[5]*alpha[9]); 
  Ghat[24] = 0.1767766952966368*(alpha[7]*fUpwind[31]+fUpwind[7]*alpha[31]+alpha[3]*fUpwind[30]+fUpwind[3]*alpha[30]+alpha[16]*fUpwind[29]+fUpwind[16]*alpha[29]+alpha[1]*fUpwind[28]+fUpwind[1]*alpha[28]+alpha[18]*fUpwind[27]+fUpwind[18]*alpha[27]+alpha[21]*fUpwind[26]+fUpwind[21]*alpha[26]+alpha[8]*fUpwind[25]+fUpwind[8]*alpha[25]+alpha[0]*fUpwind[24]+fUpwind[0]*alpha[24]+alpha[6]*fUpwind[23]+fUpwind[6]*alpha[23]+alpha[11]*fUpwind[22]+fUpwind[11]*alpha[22]+alpha[9]*fUpwind[20]+fUpwind[9]*alpha[20]+alpha[14]*fUpwind[19]+fUpwind[14]*alpha[19]+alpha[12]*fUpwind[17]+fUpwind[12]*alpha[17]+alpha[2]*fUpwind[15]+fUpwind[2]*alpha[15]+alpha[4]*fUpwind[13]+fUpwind[4]*alpha[13]+alpha[5]*fUpwind[10]+fUpwind[5]*alpha[10]); 
  Ghat[25] = 0.1767766952966368*(alpha[6]*fUpwind[31]+fUpwind[6]*alpha[31]+alpha[2]*fUpwind[30]+fUpwind[2]*alpha[30]+alpha[1]*fUpwind[29]+fUpwind[1]*alpha[29]+alpha[16]*fUpwind[28]+fUpwind[16]*alpha[28]+alpha[17]*fUpwind[27]+fUpwind[17]*alpha[27]+alpha[20]*fUpwind[26]+fUpwind[20]*alpha[26]+alpha[0]*fUpwind[25]+fUpwind[0]*alpha[25]+alpha[8]*fUpwind[24]+fUpwind[8]*alpha[24]+alpha[7]*fUpwind[23]+fUpwind[7]*alpha[23]+alpha[10]*fUpwind[22]+fUpwind[10]*alpha[22]+alpha[9]*fUpwind[21]+fUpwind[9]*alpha[21]+alpha[13]*fUpwind[19]+fUpwind[13]*alpha[19]+alpha[12]*fUpwind[18]+fUpwind[12]*alpha[18]+alpha[3]*fUpwind[15]+fUpwind[3]*alpha[15]+alpha[4]*fUpwind[14]+fUpwind[4]*alpha[14]+alpha[5]*fUpwind[11]+fUpwind[5]*alpha[11]); 
  Ghat[26] = 0.1767766952966368*(alpha[5]*fUpwind[31]+fUpwind[5]*alpha[31]+alpha[12]*fUpwind[30]+fUpwind[12]*alpha[30]+alpha[13]*fUpwind[29]+fUpwind[13]*alpha[29]+alpha[14]*fUpwind[28]+fUpwind[14]*alpha[28]+alpha[15]*fUpwind[27]+fUpwind[15]*alpha[27]+alpha[0]*fUpwind[26]+fUpwind[0]*alpha[26]+alpha[20]*fUpwind[25]+fUpwind[20]*alpha[25]+alpha[21]*fUpwind[24]+fUpwind[21]*alpha[24]+alpha[22]*fUpwind[23]+fUpwind[22]*alpha[23]+alpha[1]*fUpwind[19]+fUpwind[1]*alpha[19]+alpha[2]*fUpwind[18]+fUpwind[2]*alpha[18]+alpha[3]*fUpwind[17]+fUpwind[3]*alpha[17]+alpha[4]*fUpwind[16]+fUpwind[4]*alpha[16]+alpha[6]*fUpwind[11]+fUpwind[6]*alpha[11]+alpha[7]*fUpwind[10]+fUpwind[7]*alpha[10]+alpha[8]*fUpwind[9]+fUpwind[8]*alpha[9]); 
  Ghat[27] = 0.1767766952966368*(alpha[4]*fUpwind[31]+fUpwind[4]*alpha[31]+alpha[9]*fUpwind[30]+fUpwind[9]*alpha[30]+alpha[10]*fUpwind[29]+fUpwind[10]*alpha[29]+alpha[11]*fUpwind[28]+fUpwind[11]*alpha[28]+alpha[0]*fUpwind[27]+fUpwind[0]*alpha[27]+alpha[15]*fUpwind[26]+fUpwind[15]*alpha[26]+alpha[17]*fUpwind[25]+fUpwind[17]*alpha[25]+alpha[18]*fUpwind[24]+fUpwind[18]*alpha[24]+alpha[19]*fUpwind[23]+fUpwind[19]*alpha[23]+alpha[1]*fUpwind[22]+fUpwind[1]*alpha[22]+alpha[2]*fUpwind[21]+fUpwind[2]*alpha[21]+alpha[3]*fUpwind[20]+fUpwind[3]*alpha[20]+alpha[5]*fUpwind[16]+fUpwind[5]*alpha[16]+alpha[6]*fUpwind[14]+fUpwind[6]*alpha[14]+alpha[7]*fUpwind[13]+fUpwind[7]*alpha[13]+alpha[8]*fUpwind[12]+fUpwind[8]*alpha[12]); 
  Ghat[28] = 0.1767766952966368*(alpha[3]*fUpwind[31]+fUpwind[3]*alpha[31]+alpha[7]*fUpwind[30]+fUpwind[7]*alpha[30]+alpha[8]*fUpwind[29]+fUpwind[8]*alpha[29]+alpha[0]*fUpwind[28]+fUpwind[0]*alpha[28]+alpha[11]*fUpwind[27]+fUpwind[11]*alpha[27]+alpha[14]*fUpwind[26]+fUpwind[14]*alpha[26]+alpha[16]*fUpwind[25]+fUpwind[16]*alpha[25]+alpha[1]*fUpwind[24]+fUpwind[1]*alpha[24]+alpha[2]*fUpwind[23]+fUpwind[2]*alpha[23]+alpha[18]*fUpwind[22]+fUpwind[18]*alpha[22]+alpha[19]*fUpwind[21]+fUpwind[19]*alpha[21]+alpha[4]*fUpwind[20]+fUpwind[4]*alpha[20]+alpha[5]*fUpwind[17]+fUpwind[5]*alpha[17]+alpha[6]*fUpwind[15]+fUpwind[6]*alpha[15]+alpha[9]*fUpwind[13]+fUpwind[9]*alpha[13]+alpha[10]*fUpwind[12]+fUpwind[10]*alpha[12]); 
  Ghat[29] = 0.1767766952966368*(alpha[2]*fUpwind[31]+fUpwind[2]*alpha[31]+alpha[6]*fUpwind[30]+fUpwind[6]*alpha[30]+alpha[0]*fUpwind[29]+fUpwind[0]*alpha[29]+alpha[8]*fUpwind[28]+fUpwind[8]*alpha[28]+alpha[10]*fUpwind[27]+fUpwind[10]*alpha[27]+alpha[13]*fUpwind[26]+fUpwind[13]*alpha[26]+alpha[1]*fUpwind[25]+fUpwind[1]*alpha[25]+alpha[16]*fUpwind[24]+fUpwind[16]*alpha[24]+alpha[3]*fUpwind[23]+fUpwind[3]*alpha[23]+alpha[17]*fUpwind[22]+fUpwind[17]*alpha[22]+alpha[4]*fUpwind[21]+fUpwind[4]*alpha[21]+alpha[19]*fUpwind[20]+fUpwind[19]*alpha[20]+alpha[5]*fUpwind[18]+fUpwind[5]*alpha[18]+alpha[7]*fUpwind[15]+fUpwind[7]*alpha[15]+alpha[9]*fUpwind[14]+fUpwind[9]*alpha[14]+alpha[11]*fUpwind[12]+fUpwind[11]*alpha[12]); 
  Ghat[30] = 0.1767766952966368*(alpha[1]*fUpwind[31]+fUpwind[1]*alpha[31]+alpha[0]*fUpwind[30]+fUpwind[0]*alpha[30]+alpha[6]*fUpwind[29]+fUpwind[6]*alpha[29]+alpha[7]*fUpwind[28]+fUpwind[7]*alpha[28]+alpha[9]*fUpwind[27]+fUpwind[9]*alpha[27]+alpha[12]*fUpwind[26]+fUpwind[12]*alpha[26]+alpha[2]*fUpwind[25]+fUpwind[2]*alpha[25]+alpha[3]*fUpwind[24]+fUpwind[3]*alpha[24]+alpha[16]*fUpwind[23]+fUpwind[16]*alpha[23]+alpha[4]*fUpwind[22]+fUpwind[4]*alpha[22]+alpha[17]*fUpwind[21]+fUpwind[17]*alpha[21]+alpha[18]*fUpwind[20]+fUpwind[18]*alpha[20]+alpha[5]*fUpwind[19]+fUpwind[5]*alpha[19]+alpha[8]*fUpwind[15]+fUpwind[8]*alpha[15]+alpha[10]*fUpwind[14]+fUpwind[10]*alpha[14]+alpha[11]*fUpwind[13]+fUpwind[11]*alpha[13]); 
  Ghat[31] = 0.1767766952966368*(alpha[0]*fUpwind[31]+fUpwind[0]*alpha[31]+alpha[1]*fUpwind[30]+fUpwind[1]*alpha[30]+alpha[2]*fUpwind[29]+fUpwind[2]*alpha[29]+alpha[3]*fUpwind[28]+fUpwind[3]*alpha[28]+alpha[4]*fUpwind[27]+fUpwind[4]*alpha[27]+alpha[5]*fUpwind[26]+fUpwind[5]*alpha[26]+alpha[6]*fUpwind[25]+fUpwind[6]*alpha[25]+alpha[7]*fUpwind[24]+fUpwind[7]*alpha[24]+alpha[8]*fUpwind[23]+fUpwind[8]*alpha[23]+alpha[9]*fUpwind[22]+fUpwind[9]*alpha[22]+alpha[10]*fUpwind[21]+fUpwind[10]*alpha[21]+alpha[11]*fUpwind[20]+fUpwind[11]*alpha[20]+alpha[12]*fUpwind[19]+fUpwind[12]*alpha[19]+alpha[13]*fUpwind[18]+fUpwind[13]*alpha[18]+alpha[14]*fUpwind[17]+fUpwind[14]*alpha[17]+alpha[15]*fUpwind[16]+fUpwind[15]*alpha[16]); 

  out[0] += -0.7071067811865475*Ghat[0]*dv11; 
  out[1] += -0.7071067811865475*Ghat[1]*dv11; 
  out[2] += -0.7071067811865475*Ghat[2]*dv11; 
  out[3] += -0.7071067811865475*Ghat[3]*dv11; 
  out[4] += -0.7071067811865475*Ghat[4]*dv11; 
  out[5] += -1.224744871391589*Ghat[0]*dv11; 
  out[6] += -0.7071067811865475*Ghat[5]*dv11; 
  out[7] += -0.7071067811865475*Ghat[6]*dv11; 
  out[8] += -0.7071067811865475*Ghat[7]*dv11; 
  out[9] += -0.7071067811865475*Ghat[8]*dv11; 
  out[10] += -0.7071067811865475*Ghat[9]*dv11; 
  out[11] += -0.7071067811865475*Ghat[10]*dv11; 
  out[12] += -0.7071067811865475*Ghat[11]*dv11; 
  out[13] += -1.224744871391589*Ghat[1]*dv11; 
  out[14] += -1.224744871391589*Ghat[2]*dv11; 
  out[15] += -1.224744871391589*Ghat[3]*dv11; 
  out[16] += -1.224744871391589*Ghat[4]*dv11; 
  out[17] += -0.7071067811865475*Ghat[12]*dv11; 
  out[18] += -0.7071067811865475*Ghat[13]*dv11; 
  out[19] += -0.7071067811865475*Ghat[14]*dv11; 
  out[20] += -0.7071067811865475*Ghat[15]*dv11; 
  out[21] += -1.224744871391589*Ghat[5]*dv11; 
  out[22] += -0.7071067811865475*Ghat[16]*dv11; 
  out[23] += -0.7071067811865475*Ghat[17]*dv11; 
  out[24] += -0.7071067811865475*Ghat[18]*dv11; 
  out[25] += -0.7071067811865475*Ghat[19]*dv11; 
  out[26] += -1.224744871391589*Ghat[6]*dv11; 
  out[27] += -1.224744871391589*Ghat[7]*dv11; 
  out[28] += -1.224744871391589*Ghat[8]*dv11; 
  out[29] += -1.224744871391589*Ghat[9]*dv11; 
  out[30] += -1.224744871391589*Ghat[10]*dv11; 
  out[31] += -1.224744871391589*Ghat[11]*dv11; 
  out[32] += -0.7071067811865475*Ghat[20]*dv11; 
  out[33] += -0.7071067811865475*Ghat[21]*dv11; 
  out[34] += -0.7071067811865475*Ghat[22]*dv11; 
  out[35] += -0.7071067811865475*Ghat[23]*dv11; 
  out[36] += -0.7071067811865475*Ghat[24]*dv11; 
  out[37] += -0.7071067811865475*Ghat[25]*dv11; 
  out[38] += -1.224744871391589*Ghat[12]*dv11; 
  out[39] += -1.224744871391589*Ghat[13]*dv11; 
  out[40] += -1.224744871391589*Ghat[14]*dv11; 
  out[41] += -1.224744871391589*Ghat[15]*dv11; 
  out[42] += -0.7071067811865475*Ghat[26]*dv11; 
  out[43] += -1.224744871391589*Ghat[16]*dv11; 
  out[44] += -1.224744871391589*Ghat[17]*dv11; 
  out[45] += -1.224744871391589*Ghat[18]*dv11; 
  out[46] += -1.224744871391589*Ghat[19]*dv11; 
  out[47] += -0.7071067811865475*Ghat[27]*dv11; 
  out[48] += -0.7071067811865475*Ghat[28]*dv11; 
  out[49] += -0.7071067811865475*Ghat[29]*dv11; 
  out[50] += -0.7071067811865475*Ghat[30]*dv11; 
  out[51] += -1.224744871391589*Ghat[20]*dv11; 
  out[52] += -1.224744871391589*Ghat[21]*dv11; 
  out[53] += -1.224744871391589*Ghat[22]*dv11; 
  out[54] += -1.224744871391589*Ghat[23]*dv11; 
  out[55] += -1.224744871391589*Ghat[24]*dv11; 
  out[56] += -1.224744871391589*Ghat[25]*dv11; 
  out[57] += -1.224744871391589*Ghat[26]*dv11; 
  out[58] += -0.7071067811865475*Ghat[31]*dv11; 
  out[59] += -1.224744871391589*Ghat[27]*dv11; 
  out[60] += -1.224744871391589*Ghat[28]*dv11; 
  out[61] += -1.224744871391589*Ghat[29]*dv11; 
  out[62] += -1.224744871391589*Ghat[30]*dv11; 
  out[63] += -1.224744871391589*Ghat[31]*dv11; 

  } else { 

  alpha[0] = (-1.224744871391589*B0[0]*p2_over_gamma[2])+1.224744871391589*B2[0]*p0_over_gamma[2]+0.7071067811865475*B0[0]*p2_over_gamma[0]-0.7071067811865475*B2[0]*p0_over_gamma[0]+2.0*E1[0]; 
  alpha[1] = (-1.224744871391589*B0[1]*p2_over_gamma[2])+1.224744871391589*B2[1]*p0_over_gamma[2]+2.0*E1[1]-0.7071067811865475*p0_over_gamma[0]*B2[1]+0.7071067811865475*p2_over_gamma[0]*B0[1]; 
  alpha[2] = (-1.224744871391589*B0[2]*p2_over_gamma[2])+1.224744871391589*B2[2]*p0_over_gamma[2]+2.0*E1[2]-0.7071067811865475*p0_over_gamma[0]*B2[2]+0.7071067811865475*p2_over_gamma[0]*B0[2]; 
  alpha[3] = 2.0*E1[3]+1.224744871391589*p0_over_gamma[2]*B2[3]-0.7071067811865475*p0_over_gamma[0]*B2[3]-1.224744871391589*p2_over_gamma[2]*B0[3]+0.7071067811865475*p2_over_gamma[0]*B0[3]; 
  alpha[4] = (-1.224744871391589*B0[0]*p2_over_gamma[4])+1.224744871391589*B2[0]*p0_over_gamma[4]+0.7071067811865475*B0[0]*p2_over_gamma[1]-0.7071067811865475*B2[0]*p0_over_gamma[1]; 
  alpha[5] = (-1.224744871391589*B0[0]*p2_over_gamma[6])+1.224744871391589*B2[0]*p0_over_gamma[6]+0.7071067811865475*B0[0]*p2_over_gamma[3]-0.7071067811865475*B2[0]*p0_over_gamma[3]; 
  alpha[6] = 2.0*E1[4]+1.224744871391589*p0_over_gamma[2]*B2[4]-0.7071067811865475*p0_over_gamma[0]*B2[4]-1.224744871391589*p2_over_gamma[2]*B0[4]+0.7071067811865475*p2_over_gamma[0]*B0[4]; 
  alpha[7] = 2.0*E1[5]+1.224744871391589*p0_over_gamma[2]*B2[5]-0.7071067811865475*p0_over_gamma[0]*B2[5]-1.224744871391589*p2_over_gamma[2]*B0[5]+0.7071067811865475*p2_over_gamma[0]*B0[5]; 
  alpha[8] = 2.0*E1[6]+1.224744871391589*p0_over_gamma[2]*B2[6]-0.7071067811865475*p0_over_gamma[0]*B2[6]-1.224744871391589*p2_over_gamma[2]*B0[6]+0.7071067811865475*p2_over_gamma[0]*B0[6]; 
  alpha[9] = (-1.224744871391589*B0[1]*p2_over_gamma[4])+1.224744871391589*B2[1]*p0_over_gamma[4]+0.7071067811865475*B0[1]*p2_over_gamma[1]-0.7071067811865475*B2[1]*p0_over_gamma[1]; 
  alpha[10] = (-1.224744871391589*B0[2]*p2_over_gamma[4])+1.224744871391589*B2[2]*p0_over_gamma[4]-0.7071067811865475*p0_over_gamma[1]*B2[2]+0.7071067811865475*p2_over_gamma[1]*B0[2]; 
  alpha[11] = (-1.224744871391589*B0[3]*p2_over_gamma[4])+1.224744871391589*B2[3]*p0_over_gamma[4]-0.7071067811865475*p0_over_gamma[1]*B2[3]+0.7071067811865475*p2_over_gamma[1]*B0[3]; 
  alpha[12] = (-1.224744871391589*B0[1]*p2_over_gamma[6])+1.224744871391589*B2[1]*p0_over_gamma[6]+0.7071067811865475*B0[1]*p2_over_gamma[3]-0.7071067811865475*B2[1]*p0_over_gamma[3]; 
  alpha[13] = (-1.224744871391589*B0[2]*p2_over_gamma[6])+1.224744871391589*B2[2]*p0_over_gamma[6]+0.7071067811865475*B0[2]*p2_over_gamma[3]-0.7071067811865475*B2[2]*p0_over_gamma[3]; 
  alpha[14] = (-1.224744871391589*B0[3]*p2_over_gamma[6])+1.224744871391589*B2[3]*p0_over_gamma[6]+0.7071067811865475*B0[3]*p2_over_gamma[3]-0.7071067811865475*B2[3]*p0_over_gamma[3]; 
  alpha[15] = (-1.224744871391589*B0[0]*p2_over_gamma[7])+1.224744871391589*B2[0]*p0_over_gamma[7]+0.7071067811865475*B0[0]*p2_over_gamma[5]-0.7071067811865475*B2[0]*p0_over_gamma[5]; 
  alpha[16] = 2.0*E1[7]+1.224744871391589*p0_over_gamma[2]*B2[7]-0.7071067811865475*p0_over_gamma[0]*B2[7]-1.224744871391589*p2_over_gamma[2]*B0[7]+0.7071067811865475*p2_over_gamma[0]*B0[7]; 
  alpha[17] = (-1.224744871391589*B0[4]*p2_over_gamma[4])+1.224744871391589*B2[4]*p0_over_gamma[4]-0.7071067811865475*p0_over_gamma[1]*B2[4]+0.7071067811865475*p2_over_gamma[1]*B0[4]; 
  alpha[18] = 1.224744871391589*p0_over_gamma[4]*B2[5]-0.7071067811865475*p0_over_gamma[1]*B2[5]-1.224744871391589*p2_over_gamma[4]*B0[5]+0.7071067811865475*p2_over_gamma[1]*B0[5]; 
  alpha[19] = 1.224744871391589*p0_over_gamma[4]*B2[6]-0.7071067811865475*p0_over_gamma[1]*B2[6]-1.224744871391589*p2_over_gamma[4]*B0[6]+0.7071067811865475*p2_over_gamma[1]*B0[6]; 
  alpha[20] = (-1.224744871391589*B0[4]*p2_over_gamma[6])+1.224744871391589*B2[4]*p0_over_gamma[6]-0.7071067811865475*p0_over_gamma[3]*B2[4]+0.7071067811865475*p2_over_gamma[3]*B0[4]; 
  alpha[21] = (-1.224744871391589*B0[5]*p2_over_gamma[6])+1.224744871391589*B2[5]*p0_over_gamma[6]-0.7071067811865475*p0_over_gamma[3]*B2[5]+0.7071067811865475*p2_over_gamma[3]*B0[5]; 
  alpha[22] = (-1.224744871391589*B0[6]*p2_over_gamma[6])+1.224744871391589*B2[6]*p0_over_gamma[6]-0.7071067811865475*p0_over_gamma[3]*B2[6]+0.7071067811865475*p2_over_gamma[3]*B0[6]; 
  alpha[23] = (-1.224744871391589*B0[1]*p2_over_gamma[7])+1.224744871391589*B2[1]*p0_over_gamma[7]+0.7071067811865475*B0[1]*p2_over_gamma[5]-0.7071067811865475*B2[1]*p0_over_gamma[5]; 
  alpha[24] = (-1.224744871391589*B0[2]*p2_over_gamma[7])+1.224744871391589*B2[2]*p0_over_gamma[7]+0.7071067811865475*B0[2]*p2_over_gamma[5]-0.7071067811865475*B2[2]*p0_over_gamma[5]; 
  alpha[25] = (-1.224744871391589*B0[3]*p2_over_gamma[7])+1.224744871391589*B2[3]*p0_over_gamma[7]+0.7071067811865475*B0[3]*p2_over_gamma[5]-0.7071067811865475*B2[3]*p0_over_gamma[5]; 
  alpha[26] = 1.224744871391589*p0_over_gamma[4]*B2[7]-0.7071067811865475*p0_over_gamma[1]*B2[7]-1.224744871391589*p2_over_gamma[4]*B0[7]+0.7071067811865475*p2_over_gamma[1]*B0[7]; 
  alpha[27] = 1.224744871391589*p0_over_gamma[6]*B2[7]-0.7071067811865475*p0_over_gamma[3]*B2[7]-1.224744871391589*p2_over_gamma[6]*B0[7]+0.7071067811865475*p2_over_gamma[3]*B0[7]; 
  alpha[28] = (-1.224744871391589*B0[4]*p2_over_gamma[7])+1.224744871391589*B2[4]*p0_over_gamma[7]+0.7071067811865475*B0[4]*p2_over_gamma[5]-0.7071067811865475*B2[4]*p0_over_gamma[5]; 
  alpha[29] = (-1.224744871391589*B0[5]*p2_over_gamma[7])+1.224744871391589*B2[5]*p0_over_gamma[7]+0.7071067811865475*B0[5]*p2_over_gamma[5]-0.7071067811865475*B2[5]*p0_over_gamma[5]; 
  alpha[30] = (-1.224744871391589*B0[6]*p2_over_gamma[7])+1.224744871391589*B2[6]*p0_over_gamma[7]-0.7071067811865475*p0_over_gamma[5]*B2[6]+0.7071067811865475*p2_over_gamma[5]*B0[6]; 
  alpha[31] = (-1.224744871391589*B0[7]*p2_over_gamma[7])+1.224744871391589*B2[7]*p0_over_gamma[7]-0.7071067811865475*p0_over_gamma[5]*B2[7]+0.7071067811865475*p2_over_gamma[5]*B0[7]; 

  if ((-alpha[31])+alpha[30]+alpha[29]+alpha[28]+alpha[27]+alpha[26]-alpha[25]-alpha[24]-alpha[23]-alpha[22]-alpha[21]-alpha[20]-alpha[19]-alpha[18]-alpha[17]-alpha[16]+alpha[15]+alpha[14]+alpha[13]+alpha[12]+alpha[11]+alpha[10]+alpha[9]+alpha[8]+alpha[7]+alpha[6]-alpha[5]-alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[0] = ser_6x_p1_surfx5_eval_quad_node_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_6x_p1_surfx5_eval_quad_node_0_l(fSkin); 
  } 
  if (alpha[31]-alpha[30]-alpha[29]-alpha[28]-alpha[27]+alpha[26]+alpha[25]+alpha[24]+alpha[23]+alpha[22]+alpha[21]+alpha[20]-alpha[19]-alpha[18]-alpha[17]-alpha[16]-alpha[15]-alpha[14]-alpha[13]-alpha[12]+alpha[11]+alpha[10]+alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]-alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[1] = ser_6x_p1_surfx5_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = ser_6x_p1_surfx5_eval_quad_node_1_l(fSkin); 
  } 
  if (alpha[31]-alpha[30]-alpha[29]-alpha[28]+alpha[27]-alpha[26]+alpha[25]+alpha[24]+alpha[23]-alpha[22]-alpha[21]-alpha[20]+alpha[19]+alpha[18]+alpha[17]-alpha[16]-alpha[15]+alpha[14]+alpha[13]+alpha[12]-alpha[11]-alpha[10]-alpha[9]+alpha[8]+alpha[7]+alpha[6]-alpha[5]+alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[2] = ser_6x_p1_surfx5_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[2] = ser_6x_p1_surfx5_eval_quad_node_2_l(fSkin); 
  } 
  if ((-alpha[31])+alpha[30]+alpha[29]+alpha[28]-alpha[27]-alpha[26]-alpha[25]-alpha[24]-alpha[23]+alpha[22]+alpha[21]+alpha[20]+alpha[19]+alpha[18]+alpha[17]-alpha[16]+alpha[15]-alpha[14]-alpha[13]-alpha[12]-alpha[11]-alpha[10]-alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]+alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[3] = ser_6x_p1_surfx5_eval_quad_node_3_r(fEdge); 
  } else { 
    fUpwindQuad[3] = ser_6x_p1_surfx5_eval_quad_node_3_l(fSkin); 
  } 
  if (alpha[31]-alpha[30]-alpha[29]+alpha[28]-alpha[27]-alpha[26]+alpha[25]-alpha[24]-alpha[23]+alpha[22]+alpha[21]-alpha[20]+alpha[19]+alpha[18]-alpha[17]+alpha[16]+alpha[15]-alpha[14]+alpha[13]+alpha[12]-alpha[11]+alpha[10]+alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]-alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[4] = ser_6x_p1_surfx5_eval_quad_node_4_r(fEdge); 
  } else { 
    fUpwindQuad[4] = ser_6x_p1_surfx5_eval_quad_node_4_l(fSkin); 
  } 
  if ((-alpha[31])+alpha[30]+alpha[29]-alpha[28]+alpha[27]-alpha[26]-alpha[25]+alpha[24]+alpha[23]-alpha[22]-alpha[21]+alpha[20]+alpha[19]+alpha[18]-alpha[17]+alpha[16]-alpha[15]+alpha[14]-alpha[13]-alpha[12]-alpha[11]+alpha[10]+alpha[9]-alpha[8]-alpha[7]+alpha[6]+alpha[5]-alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[5] = ser_6x_p1_surfx5_eval_quad_node_5_r(fEdge); 
  } else { 
    fUpwindQuad[5] = ser_6x_p1_surfx5_eval_quad_node_5_l(fSkin); 
  } 
  if ((-alpha[31])+alpha[30]+alpha[29]-alpha[28]-alpha[27]+alpha[26]-alpha[25]+alpha[24]+alpha[23]+alpha[22]+alpha[21]-alpha[20]-alpha[19]-alpha[18]+alpha[17]+alpha[16]-alpha[15]-alpha[14]+alpha[13]+alpha[12]+alpha[11]-alpha[10]-alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]+alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[6] = ser_6x_p1_surfx5_eval_quad_node_6_r(fEdge); 
  } else { 
    fUpwindQuad[6] = ser_6x_p1_surfx5_eval_quad_node_6_l(fSkin); 
  } 
  if (alpha[31]-alpha[30]-alpha[29]+alpha[28]+alpha[27]+alpha[26]+alpha[25]-alpha[24]-alpha[23]-alpha[22]-alpha[21]+alpha[20]-alpha[19]-alpha[18]+alpha[17]+alpha[16]+alpha[15]+alpha[14]-alpha[13]-alpha[12]+alpha[11]-alpha[10]-alpha[9]-alpha[8]-alpha[7]+alpha[6]+alpha[5]+alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[7] = ser_6x_p1_surfx5_eval_quad_node_7_r(fEdge); 
  } else { 
    fUpwindQuad[7] = ser_6x_p1_surfx5_eval_quad_node_7_l(fSkin); 
  } 
  if (alpha[31]-alpha[30]+alpha[29]-alpha[28]-alpha[27]-alpha[26]-alpha[25]+alpha[24]-alpha[23]+alpha[22]-alpha[21]+alpha[20]+alpha[19]-alpha[18]+alpha[17]+alpha[16]+alpha[15]+alpha[14]-alpha[13]+alpha[12]+alpha[11]-alpha[10]+alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]-alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[8] = ser_6x_p1_surfx5_eval_quad_node_8_r(fEdge); 
  } else { 
    fUpwindQuad[8] = ser_6x_p1_surfx5_eval_quad_node_8_l(fSkin); 
  } 
  if ((-alpha[31])+alpha[30]-alpha[29]+alpha[28]+alpha[27]-alpha[26]+alpha[25]-alpha[24]+alpha[23]-alpha[22]+alpha[21]-alpha[20]+alpha[19]-alpha[18]+alpha[17]+alpha[16]-alpha[15]-alpha[14]+alpha[13]-alpha[12]+alpha[11]-alpha[10]+alpha[9]-alpha[8]+alpha[7]-alpha[6]+alpha[5]-alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[9] = ser_6x_p1_surfx5_eval_quad_node_9_r(fEdge); 
  } else { 
    fUpwindQuad[9] = ser_6x_p1_surfx5_eval_quad_node_9_l(fSkin); 
  } 
  if ((-alpha[31])+alpha[30]-alpha[29]+alpha[28]-alpha[27]+alpha[26]+alpha[25]-alpha[24]+alpha[23]+alpha[22]-alpha[21]+alpha[20]-alpha[19]+alpha[18]-alpha[17]+alpha[16]-alpha[15]+alpha[14]-alpha[13]+alpha[12]-alpha[11]+alpha[10]-alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]+alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[10] = ser_6x_p1_surfx5_eval_quad_node_10_r(fEdge); 
  } else { 
    fUpwindQuad[10] = ser_6x_p1_surfx5_eval_quad_node_10_l(fSkin); 
  } 
  if (alpha[31]-alpha[30]+alpha[29]-alpha[28]+alpha[27]+alpha[26]-alpha[25]+alpha[24]-alpha[23]-alpha[22]+alpha[21]-alpha[20]-alpha[19]+alpha[18]-alpha[17]+alpha[16]+alpha[15]-alpha[14]+alpha[13]-alpha[12]-alpha[11]+alpha[10]-alpha[9]-alpha[8]+alpha[7]-alpha[6]+alpha[5]+alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[11] = ser_6x_p1_surfx5_eval_quad_node_11_r(fEdge); 
  } else { 
    fUpwindQuad[11] = ser_6x_p1_surfx5_eval_quad_node_11_l(fSkin); 
  } 
  if ((-alpha[31])+alpha[30]-alpha[29]-alpha[28]+alpha[27]+alpha[26]+alpha[25]+alpha[24]-alpha[23]-alpha[22]+alpha[21]+alpha[20]-alpha[19]+alpha[18]+alpha[17]-alpha[16]+alpha[15]-alpha[14]-alpha[13]+alpha[12]-alpha[11]-alpha[10]+alpha[9]+alpha[8]-alpha[7]-alpha[6]-alpha[5]-alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[12] = ser_6x_p1_surfx5_eval_quad_node_12_r(fEdge); 
  } else { 
    fUpwindQuad[12] = ser_6x_p1_surfx5_eval_quad_node_12_l(fSkin); 
  } 
  if (alpha[31]-alpha[30]+alpha[29]+alpha[28]-alpha[27]+alpha[26]-alpha[25]-alpha[24]+alpha[23]+alpha[22]-alpha[21]-alpha[20]-alpha[19]+alpha[18]+alpha[17]-alpha[16]-alpha[15]+alpha[14]+alpha[13]-alpha[12]-alpha[11]-alpha[10]+alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]-alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[13] = ser_6x_p1_surfx5_eval_quad_node_13_r(fEdge); 
  } else { 
    fUpwindQuad[13] = ser_6x_p1_surfx5_eval_quad_node_13_l(fSkin); 
  } 
  if (alpha[31]-alpha[30]+alpha[29]+alpha[28]+alpha[27]-alpha[26]-alpha[25]-alpha[24]+alpha[23]-alpha[22]+alpha[21]+alpha[20]+alpha[19]-alpha[18]-alpha[17]-alpha[16]-alpha[15]-alpha[14]-alpha[13]+alpha[12]+alpha[11]+alpha[10]-alpha[9]+alpha[8]-alpha[7]-alpha[6]-alpha[5]+alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[14] = ser_6x_p1_surfx5_eval_quad_node_14_r(fEdge); 
  } else { 
    fUpwindQuad[14] = ser_6x_p1_surfx5_eval_quad_node_14_l(fSkin); 
  } 
  if ((-alpha[31])+alpha[30]-alpha[29]-alpha[28]-alpha[27]-alpha[26]+alpha[25]+alpha[24]-alpha[23]+alpha[22]-alpha[21]-alpha[20]+alpha[19]-alpha[18]-alpha[17]-alpha[16]+alpha[15]+alpha[14]+alpha[13]-alpha[12]+alpha[11]+alpha[10]-alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[15] = ser_6x_p1_surfx5_eval_quad_node_15_r(fEdge); 
  } else { 
    fUpwindQuad[15] = ser_6x_p1_surfx5_eval_quad_node_15_l(fSkin); 
  } 
  if (alpha[31]+alpha[30]-alpha[29]-alpha[28]-alpha[27]-alpha[26]-alpha[25]-alpha[24]+alpha[23]-alpha[22]+alpha[21]+alpha[20]-alpha[19]+alpha[18]+alpha[17]+alpha[16]+alpha[15]+alpha[14]+alpha[13]-alpha[12]+alpha[11]+alpha[10]-alpha[9]+alpha[8]-alpha[7]-alpha[6]-alpha[5]-alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[16] = ser_6x_p1_surfx5_eval_quad_node_16_r(fEdge); 
  } else { 
    fUpwindQuad[16] = ser_6x_p1_surfx5_eval_quad_node_16_l(fSkin); 
  } 
  if ((-alpha[31])-alpha[30]+alpha[29]+alpha[28]+alpha[27]-alpha[26]+alpha[25]+alpha[24]-alpha[23]+alpha[22]-alpha[21]-alpha[20]-alpha[19]+alpha[18]+alpha[17]+alpha[16]-alpha[15]-alpha[14]-alpha[13]+alpha[12]+alpha[11]+alpha[10]-alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]-alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[17] = ser_6x_p1_surfx5_eval_quad_node_17_r(fEdge); 
  } else { 
    fUpwindQuad[17] = ser_6x_p1_surfx5_eval_quad_node_17_l(fSkin); 
  } 
  if ((-alpha[31])-alpha[30]+alpha[29]+alpha[28]-alpha[27]+alpha[26]+alpha[25]+alpha[24]-alpha[23]-alpha[22]+alpha[21]+alpha[20]+alpha[19]-alpha[18]-alpha[17]+alpha[16]-alpha[15]+alpha[14]+alpha[13]-alpha[12]-alpha[11]-alpha[10]+alpha[9]+alpha[8]-alpha[7]-alpha[6]-alpha[5]+alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[18] = ser_6x_p1_surfx5_eval_quad_node_18_r(fEdge); 
  } else { 
    fUpwindQuad[18] = ser_6x_p1_surfx5_eval_quad_node_18_l(fSkin); 
  } 
  if (alpha[31]+alpha[30]-alpha[29]-alpha[28]+alpha[27]+alpha[26]-alpha[25]-alpha[24]+alpha[23]+alpha[22]-alpha[21]-alpha[20]+alpha[19]-alpha[18]-alpha[17]+alpha[16]+alpha[15]-alpha[14]-alpha[13]+alpha[12]-alpha[11]-alpha[10]+alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]+alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[19] = ser_6x_p1_surfx5_eval_quad_node_19_r(fEdge); 
  } else { 
    fUpwindQuad[19] = ser_6x_p1_surfx5_eval_quad_node_19_l(fSkin); 
  } 
  if ((-alpha[31])-alpha[30]+alpha[29]-alpha[28]+alpha[27]+alpha[26]+alpha[25]-alpha[24]+alpha[23]+alpha[22]-alpha[21]+alpha[20]+alpha[19]-alpha[18]+alpha[17]-alpha[16]+alpha[15]-alpha[14]+alpha[13]-alpha[12]-alpha[11]+alpha[10]-alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]-alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[20] = ser_6x_p1_surfx5_eval_quad_node_20_r(fEdge); 
  } else { 
    fUpwindQuad[20] = ser_6x_p1_surfx5_eval_quad_node_20_l(fSkin); 
  } 
  if (alpha[31]+alpha[30]-alpha[29]+alpha[28]-alpha[27]+alpha[26]-alpha[25]+alpha[24]-alpha[23]-alpha[22]+alpha[21]-alpha[20]+alpha[19]-alpha[18]+alpha[17]-alpha[16]-alpha[15]+alpha[14]-alpha[13]+alpha[12]-alpha[11]+alpha[10]-alpha[9]-alpha[8]+alpha[7]-alpha[6]+alpha[5]-alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[21] = ser_6x_p1_surfx5_eval_quad_node_21_r(fEdge); 
  } else { 
    fUpwindQuad[21] = ser_6x_p1_surfx5_eval_quad_node_21_l(fSkin); 
  } 
  if (alpha[31]+alpha[30]-alpha[29]+alpha[28]+alpha[27]-alpha[26]-alpha[25]+alpha[24]-alpha[23]+alpha[22]-alpha[21]+alpha[20]-alpha[19]+alpha[18]-alpha[17]-alpha[16]-alpha[15]-alpha[14]+alpha[13]-alpha[12]+alpha[11]-alpha[10]+alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]+alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[22] = ser_6x_p1_surfx5_eval_quad_node_22_r(fEdge); 
  } else { 
    fUpwindQuad[22] = ser_6x_p1_surfx5_eval_quad_node_22_l(fSkin); 
  } 
  if ((-alpha[31])-alpha[30]+alpha[29]-alpha[28]-alpha[27]-alpha[26]+alpha[25]-alpha[24]+alpha[23]-alpha[22]+alpha[21]-alpha[20]-alpha[19]+alpha[18]-alpha[17]-alpha[16]+alpha[15]+alpha[14]-alpha[13]+alpha[12]+alpha[11]-alpha[10]+alpha[9]-alpha[8]+alpha[7]-alpha[6]+alpha[5]+alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[23] = ser_6x_p1_surfx5_eval_quad_node_23_r(fEdge); 
  } else { 
    fUpwindQuad[23] = ser_6x_p1_surfx5_eval_quad_node_23_l(fSkin); 
  } 
  if ((-alpha[31])-alpha[30]-alpha[29]+alpha[28]+alpha[27]+alpha[26]-alpha[25]+alpha[24]+alpha[23]+alpha[22]+alpha[21]-alpha[20]+alpha[19]+alpha[18]-alpha[17]-alpha[16]+alpha[15]+alpha[14]-alpha[13]-alpha[12]+alpha[11]-alpha[10]-alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]-alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[24] = ser_6x_p1_surfx5_eval_quad_node_24_r(fEdge); 
  } else { 
    fUpwindQuad[24] = ser_6x_p1_surfx5_eval_quad_node_24_l(fSkin); 
  } 
  if (alpha[31]+alpha[30]+alpha[29]-alpha[28]-alpha[27]+alpha[26]+alpha[25]-alpha[24]-alpha[23]-alpha[22]-alpha[21]+alpha[20]+alpha[19]+alpha[18]-alpha[17]-alpha[16]-alpha[15]-alpha[14]+alpha[13]+alpha[12]+alpha[11]-alpha[10]-alpha[9]-alpha[8]-alpha[7]+alpha[6]+alpha[5]-alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[25] = ser_6x_p1_surfx5_eval_quad_node_25_r(fEdge); 
  } else { 
    fUpwindQuad[25] = ser_6x_p1_surfx5_eval_quad_node_25_l(fSkin); 
  } 
  if (alpha[31]+alpha[30]+alpha[29]-alpha[28]+alpha[27]-alpha[26]+alpha[25]-alpha[24]-alpha[23]+alpha[22]+alpha[21]-alpha[20]-alpha[19]-alpha[18]+alpha[17]-alpha[16]-alpha[15]+alpha[14]-alpha[13]-alpha[12]-alpha[11]+alpha[10]+alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]+alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[26] = ser_6x_p1_surfx5_eval_quad_node_26_r(fEdge); 
  } else { 
    fUpwindQuad[26] = ser_6x_p1_surfx5_eval_quad_node_26_l(fSkin); 
  } 
  if ((-alpha[31])-alpha[30]-alpha[29]+alpha[28]-alpha[27]-alpha[26]-alpha[25]+alpha[24]+alpha[23]-alpha[22]-alpha[21]+alpha[20]-alpha[19]-alpha[18]+alpha[17]-alpha[16]+alpha[15]-alpha[14]+alpha[13]+alpha[12]-alpha[11]+alpha[10]+alpha[9]-alpha[8]-alpha[7]+alpha[6]+alpha[5]+alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[27] = ser_6x_p1_surfx5_eval_quad_node_27_r(fEdge); 
  } else { 
    fUpwindQuad[27] = ser_6x_p1_surfx5_eval_quad_node_27_l(fSkin); 
  } 
  if (alpha[31]+alpha[30]+alpha[29]+alpha[28]-alpha[27]-alpha[26]+alpha[25]+alpha[24]+alpha[23]-alpha[22]-alpha[21]-alpha[20]-alpha[19]-alpha[18]-alpha[17]+alpha[16]+alpha[15]-alpha[14]-alpha[13]-alpha[12]-alpha[11]-alpha[10]-alpha[9]+alpha[8]+alpha[7]+alpha[6]-alpha[5]-alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[28] = ser_6x_p1_surfx5_eval_quad_node_28_r(fEdge); 
  } else { 
    fUpwindQuad[28] = ser_6x_p1_surfx5_eval_quad_node_28_l(fSkin); 
  } 
  if ((-alpha[31])-alpha[30]-alpha[29]-alpha[28]+alpha[27]-alpha[26]-alpha[25]-alpha[24]-alpha[23]+alpha[22]+alpha[21]+alpha[20]-alpha[19]-alpha[18]-alpha[17]+alpha[16]-alpha[15]+alpha[14]+alpha[13]+alpha[12]-alpha[11]-alpha[10]-alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]-alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[29] = ser_6x_p1_surfx5_eval_quad_node_29_r(fEdge); 
  } else { 
    fUpwindQuad[29] = ser_6x_p1_surfx5_eval_quad_node_29_l(fSkin); 
  } 
  if ((-alpha[31])-alpha[30]-alpha[29]-alpha[28]-alpha[27]+alpha[26]-alpha[25]-alpha[24]-alpha[23]-alpha[22]-alpha[21]-alpha[20]+alpha[19]+alpha[18]+alpha[17]+alpha[16]-alpha[15]-alpha[14]-alpha[13]-alpha[12]+alpha[11]+alpha[10]+alpha[9]+alpha[8]+alpha[7]+alpha[6]-alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[30] = ser_6x_p1_surfx5_eval_quad_node_30_r(fEdge); 
  } else { 
    fUpwindQuad[30] = ser_6x_p1_surfx5_eval_quad_node_30_l(fSkin); 
  } 
  if (alpha[31]+alpha[30]+alpha[29]+alpha[28]+alpha[27]+alpha[26]+alpha[25]+alpha[24]+alpha[23]+alpha[22]+alpha[21]+alpha[20]+alpha[19]+alpha[18]+alpha[17]+alpha[16]+alpha[15]+alpha[14]+alpha[13]+alpha[12]+alpha[11]+alpha[10]+alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[31] = ser_6x_p1_surfx5_eval_quad_node_31_r(fEdge); 
  } else { 
    fUpwindQuad[31] = ser_6x_p1_surfx5_eval_quad_node_31_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_6x_p1_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.1767766952966368*(alpha[31]*fUpwind[31]+alpha[30]*fUpwind[30]+alpha[29]*fUpwind[29]+alpha[28]*fUpwind[28]+alpha[27]*fUpwind[27]+alpha[26]*fUpwind[26]+alpha[25]*fUpwind[25]+alpha[24]*fUpwind[24]+alpha[23]*fUpwind[23]+alpha[22]*fUpwind[22]+alpha[21]*fUpwind[21]+alpha[20]*fUpwind[20]+alpha[19]*fUpwind[19]+alpha[18]*fUpwind[18]+alpha[17]*fUpwind[17]+alpha[16]*fUpwind[16]+alpha[15]*fUpwind[15]+alpha[14]*fUpwind[14]+alpha[13]*fUpwind[13]+alpha[12]*fUpwind[12]+alpha[11]*fUpwind[11]+alpha[10]*fUpwind[10]+alpha[9]*fUpwind[9]+alpha[8]*fUpwind[8]+alpha[7]*fUpwind[7]+alpha[6]*fUpwind[6]+alpha[5]*fUpwind[5]+alpha[4]*fUpwind[4]+alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] = 0.1767766952966368*(alpha[30]*fUpwind[31]+fUpwind[30]*alpha[31]+alpha[25]*fUpwind[29]+fUpwind[25]*alpha[29]+alpha[24]*fUpwind[28]+fUpwind[24]*alpha[28]+alpha[22]*fUpwind[27]+fUpwind[22]*alpha[27]+alpha[19]*fUpwind[26]+fUpwind[19]*alpha[26]+alpha[15]*fUpwind[23]+fUpwind[15]*alpha[23]+alpha[14]*fUpwind[21]+fUpwind[14]*alpha[21]+alpha[13]*fUpwind[20]+fUpwind[13]*alpha[20]+alpha[11]*fUpwind[18]+fUpwind[11]*alpha[18]+alpha[10]*fUpwind[17]+fUpwind[10]*alpha[17]+alpha[8]*fUpwind[16]+fUpwind[8]*alpha[16]+alpha[5]*fUpwind[12]+fUpwind[5]*alpha[12]+alpha[4]*fUpwind[9]+fUpwind[4]*alpha[9]+alpha[3]*fUpwind[7]+fUpwind[3]*alpha[7]+alpha[2]*fUpwind[6]+fUpwind[2]*alpha[6]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] = 0.1767766952966368*(alpha[29]*fUpwind[31]+fUpwind[29]*alpha[31]+alpha[25]*fUpwind[30]+fUpwind[25]*alpha[30]+alpha[23]*fUpwind[28]+fUpwind[23]*alpha[28]+alpha[21]*fUpwind[27]+fUpwind[21]*alpha[27]+alpha[18]*fUpwind[26]+fUpwind[18]*alpha[26]+alpha[15]*fUpwind[24]+fUpwind[15]*alpha[24]+alpha[14]*fUpwind[22]+fUpwind[14]*alpha[22]+alpha[12]*fUpwind[20]+fUpwind[12]*alpha[20]+alpha[11]*fUpwind[19]+fUpwind[11]*alpha[19]+alpha[9]*fUpwind[17]+fUpwind[9]*alpha[17]+alpha[7]*fUpwind[16]+fUpwind[7]*alpha[16]+alpha[5]*fUpwind[13]+fUpwind[5]*alpha[13]+alpha[4]*fUpwind[10]+fUpwind[4]*alpha[10]+alpha[3]*fUpwind[8]+fUpwind[3]*alpha[8]+alpha[1]*fUpwind[6]+fUpwind[1]*alpha[6]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] = 0.1767766952966368*(alpha[28]*fUpwind[31]+fUpwind[28]*alpha[31]+alpha[24]*fUpwind[30]+fUpwind[24]*alpha[30]+alpha[23]*fUpwind[29]+fUpwind[23]*alpha[29]+alpha[20]*fUpwind[27]+fUpwind[20]*alpha[27]+alpha[17]*fUpwind[26]+fUpwind[17]*alpha[26]+alpha[15]*fUpwind[25]+fUpwind[15]*alpha[25]+alpha[13]*fUpwind[22]+fUpwind[13]*alpha[22]+alpha[12]*fUpwind[21]+fUpwind[12]*alpha[21]+alpha[10]*fUpwind[19]+fUpwind[10]*alpha[19]+alpha[9]*fUpwind[18]+fUpwind[9]*alpha[18]+alpha[6]*fUpwind[16]+fUpwind[6]*alpha[16]+alpha[5]*fUpwind[14]+fUpwind[5]*alpha[14]+alpha[4]*fUpwind[11]+fUpwind[4]*alpha[11]+alpha[2]*fUpwind[8]+fUpwind[2]*alpha[8]+alpha[1]*fUpwind[7]+fUpwind[1]*alpha[7]+alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]); 
  Ghat[4] = 0.1767766952966368*(alpha[27]*fUpwind[31]+fUpwind[27]*alpha[31]+alpha[22]*fUpwind[30]+fUpwind[22]*alpha[30]+alpha[21]*fUpwind[29]+fUpwind[21]*alpha[29]+alpha[20]*fUpwind[28]+fUpwind[20]*alpha[28]+alpha[16]*fUpwind[26]+fUpwind[16]*alpha[26]+alpha[14]*fUpwind[25]+fUpwind[14]*alpha[25]+alpha[13]*fUpwind[24]+fUpwind[13]*alpha[24]+alpha[12]*fUpwind[23]+fUpwind[12]*alpha[23]+alpha[8]*fUpwind[19]+fUpwind[8]*alpha[19]+alpha[7]*fUpwind[18]+fUpwind[7]*alpha[18]+alpha[6]*fUpwind[17]+fUpwind[6]*alpha[17]+alpha[5]*fUpwind[15]+fUpwind[5]*alpha[15]+alpha[3]*fUpwind[11]+fUpwind[3]*alpha[11]+alpha[2]*fUpwind[10]+fUpwind[2]*alpha[10]+alpha[1]*fUpwind[9]+fUpwind[1]*alpha[9]+alpha[0]*fUpwind[4]+fUpwind[0]*alpha[4]); 
  Ghat[5] = 0.1767766952966368*(alpha[26]*fUpwind[31]+fUpwind[26]*alpha[31]+alpha[19]*fUpwind[30]+fUpwind[19]*alpha[30]+alpha[18]*fUpwind[29]+fUpwind[18]*alpha[29]+alpha[17]*fUpwind[28]+fUpwind[17]*alpha[28]+alpha[16]*fUpwind[27]+fUpwind[16]*alpha[27]+alpha[11]*fUpwind[25]+fUpwind[11]*alpha[25]+alpha[10]*fUpwind[24]+fUpwind[10]*alpha[24]+alpha[9]*fUpwind[23]+fUpwind[9]*alpha[23]+alpha[8]*fUpwind[22]+fUpwind[8]*alpha[22]+alpha[7]*fUpwind[21]+fUpwind[7]*alpha[21]+alpha[6]*fUpwind[20]+fUpwind[6]*alpha[20]+alpha[4]*fUpwind[15]+fUpwind[4]*alpha[15]+alpha[3]*fUpwind[14]+fUpwind[3]*alpha[14]+alpha[2]*fUpwind[13]+fUpwind[2]*alpha[13]+alpha[1]*fUpwind[12]+fUpwind[1]*alpha[12]+alpha[0]*fUpwind[5]+fUpwind[0]*alpha[5]); 
  Ghat[6] = 0.1767766952966368*(alpha[25]*fUpwind[31]+fUpwind[25]*alpha[31]+alpha[29]*fUpwind[30]+fUpwind[29]*alpha[30]+alpha[15]*fUpwind[28]+fUpwind[15]*alpha[28]+alpha[14]*fUpwind[27]+fUpwind[14]*alpha[27]+alpha[11]*fUpwind[26]+fUpwind[11]*alpha[26]+alpha[23]*fUpwind[24]+fUpwind[23]*alpha[24]+alpha[21]*fUpwind[22]+fUpwind[21]*alpha[22]+alpha[5]*fUpwind[20]+fUpwind[5]*alpha[20]+alpha[18]*fUpwind[19]+fUpwind[18]*alpha[19]+alpha[4]*fUpwind[17]+fUpwind[4]*alpha[17]+alpha[3]*fUpwind[16]+fUpwind[3]*alpha[16]+alpha[12]*fUpwind[13]+fUpwind[12]*alpha[13]+alpha[9]*fUpwind[10]+fUpwind[9]*alpha[10]+alpha[7]*fUpwind[8]+fUpwind[7]*alpha[8]+alpha[0]*fUpwind[6]+fUpwind[0]*alpha[6]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[7] = 0.1767766952966368*(alpha[24]*fUpwind[31]+fUpwind[24]*alpha[31]+alpha[28]*fUpwind[30]+fUpwind[28]*alpha[30]+alpha[15]*fUpwind[29]+fUpwind[15]*alpha[29]+alpha[13]*fUpwind[27]+fUpwind[13]*alpha[27]+alpha[10]*fUpwind[26]+fUpwind[10]*alpha[26]+alpha[23]*fUpwind[25]+fUpwind[23]*alpha[25]+alpha[20]*fUpwind[22]+fUpwind[20]*alpha[22]+alpha[5]*fUpwind[21]+fUpwind[5]*alpha[21]+alpha[17]*fUpwind[19]+fUpwind[17]*alpha[19]+alpha[4]*fUpwind[18]+fUpwind[4]*alpha[18]+alpha[2]*fUpwind[16]+fUpwind[2]*alpha[16]+alpha[12]*fUpwind[14]+fUpwind[12]*alpha[14]+alpha[9]*fUpwind[11]+fUpwind[9]*alpha[11]+alpha[6]*fUpwind[8]+fUpwind[6]*alpha[8]+alpha[0]*fUpwind[7]+fUpwind[0]*alpha[7]+alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]); 
  Ghat[8] = 0.1767766952966368*(alpha[23]*fUpwind[31]+fUpwind[23]*alpha[31]+alpha[15]*fUpwind[30]+fUpwind[15]*alpha[30]+alpha[28]*fUpwind[29]+fUpwind[28]*alpha[29]+alpha[12]*fUpwind[27]+fUpwind[12]*alpha[27]+alpha[9]*fUpwind[26]+fUpwind[9]*alpha[26]+alpha[24]*fUpwind[25]+fUpwind[24]*alpha[25]+alpha[5]*fUpwind[22]+fUpwind[5]*alpha[22]+alpha[20]*fUpwind[21]+fUpwind[20]*alpha[21]+alpha[4]*fUpwind[19]+fUpwind[4]*alpha[19]+alpha[17]*fUpwind[18]+fUpwind[17]*alpha[18]+alpha[1]*fUpwind[16]+fUpwind[1]*alpha[16]+alpha[13]*fUpwind[14]+fUpwind[13]*alpha[14]+alpha[10]*fUpwind[11]+fUpwind[10]*alpha[11]+alpha[0]*fUpwind[8]+fUpwind[0]*alpha[8]+alpha[6]*fUpwind[7]+fUpwind[6]*alpha[7]+alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]); 
  Ghat[9] = 0.1767766952966368*(alpha[22]*fUpwind[31]+fUpwind[22]*alpha[31]+alpha[27]*fUpwind[30]+fUpwind[27]*alpha[30]+alpha[14]*fUpwind[29]+fUpwind[14]*alpha[29]+alpha[13]*fUpwind[28]+fUpwind[13]*alpha[28]+alpha[8]*fUpwind[26]+fUpwind[8]*alpha[26]+alpha[21]*fUpwind[25]+fUpwind[21]*alpha[25]+alpha[20]*fUpwind[24]+fUpwind[20]*alpha[24]+alpha[5]*fUpwind[23]+fUpwind[5]*alpha[23]+alpha[16]*fUpwind[19]+fUpwind[16]*alpha[19]+alpha[3]*fUpwind[18]+fUpwind[3]*alpha[18]+alpha[2]*fUpwind[17]+fUpwind[2]*alpha[17]+alpha[12]*fUpwind[15]+fUpwind[12]*alpha[15]+alpha[7]*fUpwind[11]+fUpwind[7]*alpha[11]+alpha[6]*fUpwind[10]+fUpwind[6]*alpha[10]+alpha[0]*fUpwind[9]+fUpwind[0]*alpha[9]+alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]); 
  Ghat[10] = 0.1767766952966368*(alpha[21]*fUpwind[31]+fUpwind[21]*alpha[31]+alpha[14]*fUpwind[30]+fUpwind[14]*alpha[30]+alpha[27]*fUpwind[29]+fUpwind[27]*alpha[29]+alpha[12]*fUpwind[28]+fUpwind[12]*alpha[28]+alpha[7]*fUpwind[26]+fUpwind[7]*alpha[26]+alpha[22]*fUpwind[25]+fUpwind[22]*alpha[25]+alpha[5]*fUpwind[24]+fUpwind[5]*alpha[24]+alpha[20]*fUpwind[23]+fUpwind[20]*alpha[23]+alpha[3]*fUpwind[19]+fUpwind[3]*alpha[19]+alpha[16]*fUpwind[18]+fUpwind[16]*alpha[18]+alpha[1]*fUpwind[17]+fUpwind[1]*alpha[17]+alpha[13]*fUpwind[15]+fUpwind[13]*alpha[15]+alpha[8]*fUpwind[11]+fUpwind[8]*alpha[11]+alpha[0]*fUpwind[10]+fUpwind[0]*alpha[10]+alpha[6]*fUpwind[9]+fUpwind[6]*alpha[9]+alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]); 
  Ghat[11] = 0.1767766952966368*(alpha[20]*fUpwind[31]+fUpwind[20]*alpha[31]+alpha[13]*fUpwind[30]+fUpwind[13]*alpha[30]+alpha[12]*fUpwind[29]+fUpwind[12]*alpha[29]+alpha[27]*fUpwind[28]+fUpwind[27]*alpha[28]+alpha[6]*fUpwind[26]+fUpwind[6]*alpha[26]+alpha[5]*fUpwind[25]+fUpwind[5]*alpha[25]+alpha[22]*fUpwind[24]+fUpwind[22]*alpha[24]+alpha[21]*fUpwind[23]+fUpwind[21]*alpha[23]+alpha[2]*fUpwind[19]+fUpwind[2]*alpha[19]+alpha[1]*fUpwind[18]+fUpwind[1]*alpha[18]+alpha[16]*fUpwind[17]+fUpwind[16]*alpha[17]+alpha[14]*fUpwind[15]+fUpwind[14]*alpha[15]+alpha[0]*fUpwind[11]+fUpwind[0]*alpha[11]+alpha[8]*fUpwind[10]+fUpwind[8]*alpha[10]+alpha[7]*fUpwind[9]+fUpwind[7]*alpha[9]+alpha[3]*fUpwind[4]+fUpwind[3]*alpha[4]); 
  Ghat[12] = 0.1767766952966368*(alpha[19]*fUpwind[31]+fUpwind[19]*alpha[31]+alpha[26]*fUpwind[30]+fUpwind[26]*alpha[30]+alpha[11]*fUpwind[29]+fUpwind[11]*alpha[29]+alpha[10]*fUpwind[28]+fUpwind[10]*alpha[28]+alpha[8]*fUpwind[27]+fUpwind[8]*alpha[27]+alpha[18]*fUpwind[25]+fUpwind[18]*alpha[25]+alpha[17]*fUpwind[24]+fUpwind[17]*alpha[24]+alpha[4]*fUpwind[23]+fUpwind[4]*alpha[23]+alpha[16]*fUpwind[22]+fUpwind[16]*alpha[22]+alpha[3]*fUpwind[21]+fUpwind[3]*alpha[21]+alpha[2]*fUpwind[20]+fUpwind[2]*alpha[20]+alpha[9]*fUpwind[15]+fUpwind[9]*alpha[15]+alpha[7]*fUpwind[14]+fUpwind[7]*alpha[14]+alpha[6]*fUpwind[13]+fUpwind[6]*alpha[13]+alpha[0]*fUpwind[12]+fUpwind[0]*alpha[12]+alpha[1]*fUpwind[5]+fUpwind[1]*alpha[5]); 
  Ghat[13] = 0.1767766952966368*(alpha[18]*fUpwind[31]+fUpwind[18]*alpha[31]+alpha[11]*fUpwind[30]+fUpwind[11]*alpha[30]+alpha[26]*fUpwind[29]+fUpwind[26]*alpha[29]+alpha[9]*fUpwind[28]+fUpwind[9]*alpha[28]+alpha[7]*fUpwind[27]+fUpwind[7]*alpha[27]+alpha[19]*fUpwind[25]+fUpwind[19]*alpha[25]+alpha[4]*fUpwind[24]+fUpwind[4]*alpha[24]+alpha[17]*fUpwind[23]+fUpwind[17]*alpha[23]+alpha[3]*fUpwind[22]+fUpwind[3]*alpha[22]+alpha[16]*fUpwind[21]+fUpwind[16]*alpha[21]+alpha[1]*fUpwind[20]+fUpwind[1]*alpha[20]+alpha[10]*fUpwind[15]+fUpwind[10]*alpha[15]+alpha[8]*fUpwind[14]+fUpwind[8]*alpha[14]+alpha[0]*fUpwind[13]+fUpwind[0]*alpha[13]+alpha[6]*fUpwind[12]+fUpwind[6]*alpha[12]+alpha[2]*fUpwind[5]+fUpwind[2]*alpha[5]); 
  Ghat[14] = 0.1767766952966368*(alpha[17]*fUpwind[31]+fUpwind[17]*alpha[31]+alpha[10]*fUpwind[30]+fUpwind[10]*alpha[30]+alpha[9]*fUpwind[29]+fUpwind[9]*alpha[29]+alpha[26]*fUpwind[28]+fUpwind[26]*alpha[28]+alpha[6]*fUpwind[27]+fUpwind[6]*alpha[27]+alpha[4]*fUpwind[25]+fUpwind[4]*alpha[25]+alpha[19]*fUpwind[24]+fUpwind[19]*alpha[24]+alpha[18]*fUpwind[23]+fUpwind[18]*alpha[23]+alpha[2]*fUpwind[22]+fUpwind[2]*alpha[22]+alpha[1]*fUpwind[21]+fUpwind[1]*alpha[21]+alpha[16]*fUpwind[20]+fUpwind[16]*alpha[20]+alpha[11]*fUpwind[15]+fUpwind[11]*alpha[15]+alpha[0]*fUpwind[14]+fUpwind[0]*alpha[14]+alpha[8]*fUpwind[13]+fUpwind[8]*alpha[13]+alpha[7]*fUpwind[12]+fUpwind[7]*alpha[12]+alpha[3]*fUpwind[5]+fUpwind[3]*alpha[5]); 
  Ghat[15] = 0.1767766952966368*(alpha[16]*fUpwind[31]+fUpwind[16]*alpha[31]+alpha[8]*fUpwind[30]+fUpwind[8]*alpha[30]+alpha[7]*fUpwind[29]+fUpwind[7]*alpha[29]+alpha[6]*fUpwind[28]+fUpwind[6]*alpha[28]+alpha[26]*fUpwind[27]+fUpwind[26]*alpha[27]+alpha[3]*fUpwind[25]+fUpwind[3]*alpha[25]+alpha[2]*fUpwind[24]+fUpwind[2]*alpha[24]+alpha[1]*fUpwind[23]+fUpwind[1]*alpha[23]+alpha[19]*fUpwind[22]+fUpwind[19]*alpha[22]+alpha[18]*fUpwind[21]+fUpwind[18]*alpha[21]+alpha[17]*fUpwind[20]+fUpwind[17]*alpha[20]+alpha[0]*fUpwind[15]+fUpwind[0]*alpha[15]+alpha[11]*fUpwind[14]+fUpwind[11]*alpha[14]+alpha[10]*fUpwind[13]+fUpwind[10]*alpha[13]+alpha[9]*fUpwind[12]+fUpwind[9]*alpha[12]+alpha[4]*fUpwind[5]+fUpwind[4]*alpha[5]); 
  Ghat[16] = 0.1767766952966368*(alpha[15]*fUpwind[31]+fUpwind[15]*alpha[31]+alpha[23]*fUpwind[30]+fUpwind[23]*alpha[30]+alpha[24]*fUpwind[29]+fUpwind[24]*alpha[29]+alpha[25]*fUpwind[28]+fUpwind[25]*alpha[28]+alpha[5]*fUpwind[27]+fUpwind[5]*alpha[27]+alpha[4]*fUpwind[26]+fUpwind[4]*alpha[26]+alpha[12]*fUpwind[22]+fUpwind[12]*alpha[22]+alpha[13]*fUpwind[21]+fUpwind[13]*alpha[21]+alpha[14]*fUpwind[20]+fUpwind[14]*alpha[20]+alpha[9]*fUpwind[19]+fUpwind[9]*alpha[19]+alpha[10]*fUpwind[18]+fUpwind[10]*alpha[18]+alpha[11]*fUpwind[17]+fUpwind[11]*alpha[17]+alpha[0]*fUpwind[16]+fUpwind[0]*alpha[16]+alpha[1]*fUpwind[8]+fUpwind[1]*alpha[8]+alpha[2]*fUpwind[7]+fUpwind[2]*alpha[7]+alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6]); 
  Ghat[17] = 0.1767766952966368*(alpha[14]*fUpwind[31]+fUpwind[14]*alpha[31]+alpha[21]*fUpwind[30]+fUpwind[21]*alpha[30]+alpha[22]*fUpwind[29]+fUpwind[22]*alpha[29]+alpha[5]*fUpwind[28]+fUpwind[5]*alpha[28]+alpha[25]*fUpwind[27]+fUpwind[25]*alpha[27]+alpha[3]*fUpwind[26]+fUpwind[3]*alpha[26]+alpha[12]*fUpwind[24]+fUpwind[12]*alpha[24]+alpha[13]*fUpwind[23]+fUpwind[13]*alpha[23]+alpha[15]*fUpwind[20]+fUpwind[15]*alpha[20]+alpha[7]*fUpwind[19]+fUpwind[7]*alpha[19]+alpha[8]*fUpwind[18]+fUpwind[8]*alpha[18]+alpha[0]*fUpwind[17]+fUpwind[0]*alpha[17]+alpha[11]*fUpwind[16]+fUpwind[11]*alpha[16]+alpha[1]*fUpwind[10]+fUpwind[1]*alpha[10]+alpha[2]*fUpwind[9]+fUpwind[2]*alpha[9]+alpha[4]*fUpwind[6]+fUpwind[4]*alpha[6]); 
  Ghat[18] = 0.1767766952966368*(alpha[13]*fUpwind[31]+fUpwind[13]*alpha[31]+alpha[20]*fUpwind[30]+fUpwind[20]*alpha[30]+alpha[5]*fUpwind[29]+fUpwind[5]*alpha[29]+alpha[22]*fUpwind[28]+fUpwind[22]*alpha[28]+alpha[24]*fUpwind[27]+fUpwind[24]*alpha[27]+alpha[2]*fUpwind[26]+fUpwind[2]*alpha[26]+alpha[12]*fUpwind[25]+fUpwind[12]*alpha[25]+alpha[14]*fUpwind[23]+fUpwind[14]*alpha[23]+alpha[15]*fUpwind[21]+fUpwind[15]*alpha[21]+alpha[6]*fUpwind[19]+fUpwind[6]*alpha[19]+alpha[0]*fUpwind[18]+fUpwind[0]*alpha[18]+alpha[8]*fUpwind[17]+fUpwind[8]*alpha[17]+alpha[10]*fUpwind[16]+fUpwind[10]*alpha[16]+alpha[1]*fUpwind[11]+fUpwind[1]*alpha[11]+alpha[3]*fUpwind[9]+fUpwind[3]*alpha[9]+alpha[4]*fUpwind[7]+fUpwind[4]*alpha[7]); 
  Ghat[19] = 0.1767766952966368*(alpha[12]*fUpwind[31]+fUpwind[12]*alpha[31]+alpha[5]*fUpwind[30]+fUpwind[5]*alpha[30]+alpha[20]*fUpwind[29]+fUpwind[20]*alpha[29]+alpha[21]*fUpwind[28]+fUpwind[21]*alpha[28]+alpha[23]*fUpwind[27]+fUpwind[23]*alpha[27]+alpha[1]*fUpwind[26]+fUpwind[1]*alpha[26]+alpha[13]*fUpwind[25]+fUpwind[13]*alpha[25]+alpha[14]*fUpwind[24]+fUpwind[14]*alpha[24]+alpha[15]*fUpwind[22]+fUpwind[15]*alpha[22]+alpha[0]*fUpwind[19]+fUpwind[0]*alpha[19]+alpha[6]*fUpwind[18]+fUpwind[6]*alpha[18]+alpha[7]*fUpwind[17]+fUpwind[7]*alpha[17]+alpha[9]*fUpwind[16]+fUpwind[9]*alpha[16]+alpha[2]*fUpwind[11]+fUpwind[2]*alpha[11]+alpha[3]*fUpwind[10]+fUpwind[3]*alpha[10]+alpha[4]*fUpwind[8]+fUpwind[4]*alpha[8]); 
  Ghat[20] = 0.1767766952966368*(alpha[11]*fUpwind[31]+fUpwind[11]*alpha[31]+alpha[18]*fUpwind[30]+fUpwind[18]*alpha[30]+alpha[19]*fUpwind[29]+fUpwind[19]*alpha[29]+alpha[4]*fUpwind[28]+fUpwind[4]*alpha[28]+alpha[3]*fUpwind[27]+fUpwind[3]*alpha[27]+alpha[25]*fUpwind[26]+fUpwind[25]*alpha[26]+alpha[9]*fUpwind[24]+fUpwind[9]*alpha[24]+alpha[10]*fUpwind[23]+fUpwind[10]*alpha[23]+alpha[7]*fUpwind[22]+fUpwind[7]*alpha[22]+alpha[8]*fUpwind[21]+fUpwind[8]*alpha[21]+alpha[0]*fUpwind[20]+fUpwind[0]*alpha[20]+alpha[15]*fUpwind[17]+fUpwind[15]*alpha[17]+alpha[14]*fUpwind[16]+fUpwind[14]*alpha[16]+alpha[1]*fUpwind[13]+fUpwind[1]*alpha[13]+alpha[2]*fUpwind[12]+fUpwind[2]*alpha[12]+alpha[5]*fUpwind[6]+fUpwind[5]*alpha[6]); 
  Ghat[21] = 0.1767766952966368*(alpha[10]*fUpwind[31]+fUpwind[10]*alpha[31]+alpha[17]*fUpwind[30]+fUpwind[17]*alpha[30]+alpha[4]*fUpwind[29]+fUpwind[4]*alpha[29]+alpha[19]*fUpwind[28]+fUpwind[19]*alpha[28]+alpha[2]*fUpwind[27]+fUpwind[2]*alpha[27]+alpha[24]*fUpwind[26]+fUpwind[24]*alpha[26]+alpha[9]*fUpwind[25]+fUpwind[9]*alpha[25]+alpha[11]*fUpwind[23]+fUpwind[11]*alpha[23]+alpha[6]*fUpwind[22]+fUpwind[6]*alpha[22]+alpha[0]*fUpwind[21]+fUpwind[0]*alpha[21]+alpha[8]*fUpwind[20]+fUpwind[8]*alpha[20]+alpha[15]*fUpwind[18]+fUpwind[15]*alpha[18]+alpha[13]*fUpwind[16]+fUpwind[13]*alpha[16]+alpha[1]*fUpwind[14]+fUpwind[1]*alpha[14]+alpha[3]*fUpwind[12]+fUpwind[3]*alpha[12]+alpha[5]*fUpwind[7]+fUpwind[5]*alpha[7]); 
  Ghat[22] = 0.1767766952966368*(alpha[9]*fUpwind[31]+fUpwind[9]*alpha[31]+alpha[4]*fUpwind[30]+fUpwind[4]*alpha[30]+alpha[17]*fUpwind[29]+fUpwind[17]*alpha[29]+alpha[18]*fUpwind[28]+fUpwind[18]*alpha[28]+alpha[1]*fUpwind[27]+fUpwind[1]*alpha[27]+alpha[23]*fUpwind[26]+fUpwind[23]*alpha[26]+alpha[10]*fUpwind[25]+fUpwind[10]*alpha[25]+alpha[11]*fUpwind[24]+fUpwind[11]*alpha[24]+alpha[0]*fUpwind[22]+fUpwind[0]*alpha[22]+alpha[6]*fUpwind[21]+fUpwind[6]*alpha[21]+alpha[7]*fUpwind[20]+fUpwind[7]*alpha[20]+alpha[15]*fUpwind[19]+fUpwind[15]*alpha[19]+alpha[12]*fUpwind[16]+fUpwind[12]*alpha[16]+alpha[2]*fUpwind[14]+fUpwind[2]*alpha[14]+alpha[3]*fUpwind[13]+fUpwind[3]*alpha[13]+alpha[5]*fUpwind[8]+fUpwind[5]*alpha[8]); 
  Ghat[23] = 0.1767766952966368*(alpha[8]*fUpwind[31]+fUpwind[8]*alpha[31]+alpha[16]*fUpwind[30]+fUpwind[16]*alpha[30]+alpha[3]*fUpwind[29]+fUpwind[3]*alpha[29]+alpha[2]*fUpwind[28]+fUpwind[2]*alpha[28]+alpha[19]*fUpwind[27]+fUpwind[19]*alpha[27]+alpha[22]*fUpwind[26]+fUpwind[22]*alpha[26]+alpha[7]*fUpwind[25]+fUpwind[7]*alpha[25]+alpha[6]*fUpwind[24]+fUpwind[6]*alpha[24]+alpha[0]*fUpwind[23]+fUpwind[0]*alpha[23]+alpha[11]*fUpwind[21]+fUpwind[11]*alpha[21]+alpha[10]*fUpwind[20]+fUpwind[10]*alpha[20]+alpha[14]*fUpwind[18]+fUpwind[14]*alpha[18]+alpha[13]*fUpwind[17]+fUpwind[13]*alpha[17]+alpha[1]*fUpwind[15]+fUpwind[1]*alpha[15]+alpha[4]*fUpwind[12]+fUpwind[4]*alpha[12]+alpha[5]*fUpwind[9]+fUpwind[5]*alpha[9]); 
  Ghat[24] = 0.1767766952966368*(alpha[7]*fUpwind[31]+fUpwind[7]*alpha[31]+alpha[3]*fUpwind[30]+fUpwind[3]*alpha[30]+alpha[16]*fUpwind[29]+fUpwind[16]*alpha[29]+alpha[1]*fUpwind[28]+fUpwind[1]*alpha[28]+alpha[18]*fUpwind[27]+fUpwind[18]*alpha[27]+alpha[21]*fUpwind[26]+fUpwind[21]*alpha[26]+alpha[8]*fUpwind[25]+fUpwind[8]*alpha[25]+alpha[0]*fUpwind[24]+fUpwind[0]*alpha[24]+alpha[6]*fUpwind[23]+fUpwind[6]*alpha[23]+alpha[11]*fUpwind[22]+fUpwind[11]*alpha[22]+alpha[9]*fUpwind[20]+fUpwind[9]*alpha[20]+alpha[14]*fUpwind[19]+fUpwind[14]*alpha[19]+alpha[12]*fUpwind[17]+fUpwind[12]*alpha[17]+alpha[2]*fUpwind[15]+fUpwind[2]*alpha[15]+alpha[4]*fUpwind[13]+fUpwind[4]*alpha[13]+alpha[5]*fUpwind[10]+fUpwind[5]*alpha[10]); 
  Ghat[25] = 0.1767766952966368*(alpha[6]*fUpwind[31]+fUpwind[6]*alpha[31]+alpha[2]*fUpwind[30]+fUpwind[2]*alpha[30]+alpha[1]*fUpwind[29]+fUpwind[1]*alpha[29]+alpha[16]*fUpwind[28]+fUpwind[16]*alpha[28]+alpha[17]*fUpwind[27]+fUpwind[17]*alpha[27]+alpha[20]*fUpwind[26]+fUpwind[20]*alpha[26]+alpha[0]*fUpwind[25]+fUpwind[0]*alpha[25]+alpha[8]*fUpwind[24]+fUpwind[8]*alpha[24]+alpha[7]*fUpwind[23]+fUpwind[7]*alpha[23]+alpha[10]*fUpwind[22]+fUpwind[10]*alpha[22]+alpha[9]*fUpwind[21]+fUpwind[9]*alpha[21]+alpha[13]*fUpwind[19]+fUpwind[13]*alpha[19]+alpha[12]*fUpwind[18]+fUpwind[12]*alpha[18]+alpha[3]*fUpwind[15]+fUpwind[3]*alpha[15]+alpha[4]*fUpwind[14]+fUpwind[4]*alpha[14]+alpha[5]*fUpwind[11]+fUpwind[5]*alpha[11]); 
  Ghat[26] = 0.1767766952966368*(alpha[5]*fUpwind[31]+fUpwind[5]*alpha[31]+alpha[12]*fUpwind[30]+fUpwind[12]*alpha[30]+alpha[13]*fUpwind[29]+fUpwind[13]*alpha[29]+alpha[14]*fUpwind[28]+fUpwind[14]*alpha[28]+alpha[15]*fUpwind[27]+fUpwind[15]*alpha[27]+alpha[0]*fUpwind[26]+fUpwind[0]*alpha[26]+alpha[20]*fUpwind[25]+fUpwind[20]*alpha[25]+alpha[21]*fUpwind[24]+fUpwind[21]*alpha[24]+alpha[22]*fUpwind[23]+fUpwind[22]*alpha[23]+alpha[1]*fUpwind[19]+fUpwind[1]*alpha[19]+alpha[2]*fUpwind[18]+fUpwind[2]*alpha[18]+alpha[3]*fUpwind[17]+fUpwind[3]*alpha[17]+alpha[4]*fUpwind[16]+fUpwind[4]*alpha[16]+alpha[6]*fUpwind[11]+fUpwind[6]*alpha[11]+alpha[7]*fUpwind[10]+fUpwind[7]*alpha[10]+alpha[8]*fUpwind[9]+fUpwind[8]*alpha[9]); 
  Ghat[27] = 0.1767766952966368*(alpha[4]*fUpwind[31]+fUpwind[4]*alpha[31]+alpha[9]*fUpwind[30]+fUpwind[9]*alpha[30]+alpha[10]*fUpwind[29]+fUpwind[10]*alpha[29]+alpha[11]*fUpwind[28]+fUpwind[11]*alpha[28]+alpha[0]*fUpwind[27]+fUpwind[0]*alpha[27]+alpha[15]*fUpwind[26]+fUpwind[15]*alpha[26]+alpha[17]*fUpwind[25]+fUpwind[17]*alpha[25]+alpha[18]*fUpwind[24]+fUpwind[18]*alpha[24]+alpha[19]*fUpwind[23]+fUpwind[19]*alpha[23]+alpha[1]*fUpwind[22]+fUpwind[1]*alpha[22]+alpha[2]*fUpwind[21]+fUpwind[2]*alpha[21]+alpha[3]*fUpwind[20]+fUpwind[3]*alpha[20]+alpha[5]*fUpwind[16]+fUpwind[5]*alpha[16]+alpha[6]*fUpwind[14]+fUpwind[6]*alpha[14]+alpha[7]*fUpwind[13]+fUpwind[7]*alpha[13]+alpha[8]*fUpwind[12]+fUpwind[8]*alpha[12]); 
  Ghat[28] = 0.1767766952966368*(alpha[3]*fUpwind[31]+fUpwind[3]*alpha[31]+alpha[7]*fUpwind[30]+fUpwind[7]*alpha[30]+alpha[8]*fUpwind[29]+fUpwind[8]*alpha[29]+alpha[0]*fUpwind[28]+fUpwind[0]*alpha[28]+alpha[11]*fUpwind[27]+fUpwind[11]*alpha[27]+alpha[14]*fUpwind[26]+fUpwind[14]*alpha[26]+alpha[16]*fUpwind[25]+fUpwind[16]*alpha[25]+alpha[1]*fUpwind[24]+fUpwind[1]*alpha[24]+alpha[2]*fUpwind[23]+fUpwind[2]*alpha[23]+alpha[18]*fUpwind[22]+fUpwind[18]*alpha[22]+alpha[19]*fUpwind[21]+fUpwind[19]*alpha[21]+alpha[4]*fUpwind[20]+fUpwind[4]*alpha[20]+alpha[5]*fUpwind[17]+fUpwind[5]*alpha[17]+alpha[6]*fUpwind[15]+fUpwind[6]*alpha[15]+alpha[9]*fUpwind[13]+fUpwind[9]*alpha[13]+alpha[10]*fUpwind[12]+fUpwind[10]*alpha[12]); 
  Ghat[29] = 0.1767766952966368*(alpha[2]*fUpwind[31]+fUpwind[2]*alpha[31]+alpha[6]*fUpwind[30]+fUpwind[6]*alpha[30]+alpha[0]*fUpwind[29]+fUpwind[0]*alpha[29]+alpha[8]*fUpwind[28]+fUpwind[8]*alpha[28]+alpha[10]*fUpwind[27]+fUpwind[10]*alpha[27]+alpha[13]*fUpwind[26]+fUpwind[13]*alpha[26]+alpha[1]*fUpwind[25]+fUpwind[1]*alpha[25]+alpha[16]*fUpwind[24]+fUpwind[16]*alpha[24]+alpha[3]*fUpwind[23]+fUpwind[3]*alpha[23]+alpha[17]*fUpwind[22]+fUpwind[17]*alpha[22]+alpha[4]*fUpwind[21]+fUpwind[4]*alpha[21]+alpha[19]*fUpwind[20]+fUpwind[19]*alpha[20]+alpha[5]*fUpwind[18]+fUpwind[5]*alpha[18]+alpha[7]*fUpwind[15]+fUpwind[7]*alpha[15]+alpha[9]*fUpwind[14]+fUpwind[9]*alpha[14]+alpha[11]*fUpwind[12]+fUpwind[11]*alpha[12]); 
  Ghat[30] = 0.1767766952966368*(alpha[1]*fUpwind[31]+fUpwind[1]*alpha[31]+alpha[0]*fUpwind[30]+fUpwind[0]*alpha[30]+alpha[6]*fUpwind[29]+fUpwind[6]*alpha[29]+alpha[7]*fUpwind[28]+fUpwind[7]*alpha[28]+alpha[9]*fUpwind[27]+fUpwind[9]*alpha[27]+alpha[12]*fUpwind[26]+fUpwind[12]*alpha[26]+alpha[2]*fUpwind[25]+fUpwind[2]*alpha[25]+alpha[3]*fUpwind[24]+fUpwind[3]*alpha[24]+alpha[16]*fUpwind[23]+fUpwind[16]*alpha[23]+alpha[4]*fUpwind[22]+fUpwind[4]*alpha[22]+alpha[17]*fUpwind[21]+fUpwind[17]*alpha[21]+alpha[18]*fUpwind[20]+fUpwind[18]*alpha[20]+alpha[5]*fUpwind[19]+fUpwind[5]*alpha[19]+alpha[8]*fUpwind[15]+fUpwind[8]*alpha[15]+alpha[10]*fUpwind[14]+fUpwind[10]*alpha[14]+alpha[11]*fUpwind[13]+fUpwind[11]*alpha[13]); 
  Ghat[31] = 0.1767766952966368*(alpha[0]*fUpwind[31]+fUpwind[0]*alpha[31]+alpha[1]*fUpwind[30]+fUpwind[1]*alpha[30]+alpha[2]*fUpwind[29]+fUpwind[2]*alpha[29]+alpha[3]*fUpwind[28]+fUpwind[3]*alpha[28]+alpha[4]*fUpwind[27]+fUpwind[4]*alpha[27]+alpha[5]*fUpwind[26]+fUpwind[5]*alpha[26]+alpha[6]*fUpwind[25]+fUpwind[6]*alpha[25]+alpha[7]*fUpwind[24]+fUpwind[7]*alpha[24]+alpha[8]*fUpwind[23]+fUpwind[8]*alpha[23]+alpha[9]*fUpwind[22]+fUpwind[9]*alpha[22]+alpha[10]*fUpwind[21]+fUpwind[10]*alpha[21]+alpha[11]*fUpwind[20]+fUpwind[11]*alpha[20]+alpha[12]*fUpwind[19]+fUpwind[12]*alpha[19]+alpha[13]*fUpwind[18]+fUpwind[13]*alpha[18]+alpha[14]*fUpwind[17]+fUpwind[14]*alpha[17]+alpha[15]*fUpwind[16]+fUpwind[15]*alpha[16]); 

  out[0] += 0.7071067811865475*Ghat[0]*dv11; 
  out[1] += 0.7071067811865475*Ghat[1]*dv11; 
  out[2] += 0.7071067811865475*Ghat[2]*dv11; 
  out[3] += 0.7071067811865475*Ghat[3]*dv11; 
  out[4] += 0.7071067811865475*Ghat[4]*dv11; 
  out[5] += -1.224744871391589*Ghat[0]*dv11; 
  out[6] += 0.7071067811865475*Ghat[5]*dv11; 
  out[7] += 0.7071067811865475*Ghat[6]*dv11; 
  out[8] += 0.7071067811865475*Ghat[7]*dv11; 
  out[9] += 0.7071067811865475*Ghat[8]*dv11; 
  out[10] += 0.7071067811865475*Ghat[9]*dv11; 
  out[11] += 0.7071067811865475*Ghat[10]*dv11; 
  out[12] += 0.7071067811865475*Ghat[11]*dv11; 
  out[13] += -1.224744871391589*Ghat[1]*dv11; 
  out[14] += -1.224744871391589*Ghat[2]*dv11; 
  out[15] += -1.224744871391589*Ghat[3]*dv11; 
  out[16] += -1.224744871391589*Ghat[4]*dv11; 
  out[17] += 0.7071067811865475*Ghat[12]*dv11; 
  out[18] += 0.7071067811865475*Ghat[13]*dv11; 
  out[19] += 0.7071067811865475*Ghat[14]*dv11; 
  out[20] += 0.7071067811865475*Ghat[15]*dv11; 
  out[21] += -1.224744871391589*Ghat[5]*dv11; 
  out[22] += 0.7071067811865475*Ghat[16]*dv11; 
  out[23] += 0.7071067811865475*Ghat[17]*dv11; 
  out[24] += 0.7071067811865475*Ghat[18]*dv11; 
  out[25] += 0.7071067811865475*Ghat[19]*dv11; 
  out[26] += -1.224744871391589*Ghat[6]*dv11; 
  out[27] += -1.224744871391589*Ghat[7]*dv11; 
  out[28] += -1.224744871391589*Ghat[8]*dv11; 
  out[29] += -1.224744871391589*Ghat[9]*dv11; 
  out[30] += -1.224744871391589*Ghat[10]*dv11; 
  out[31] += -1.224744871391589*Ghat[11]*dv11; 
  out[32] += 0.7071067811865475*Ghat[20]*dv11; 
  out[33] += 0.7071067811865475*Ghat[21]*dv11; 
  out[34] += 0.7071067811865475*Ghat[22]*dv11; 
  out[35] += 0.7071067811865475*Ghat[23]*dv11; 
  out[36] += 0.7071067811865475*Ghat[24]*dv11; 
  out[37] += 0.7071067811865475*Ghat[25]*dv11; 
  out[38] += -1.224744871391589*Ghat[12]*dv11; 
  out[39] += -1.224744871391589*Ghat[13]*dv11; 
  out[40] += -1.224744871391589*Ghat[14]*dv11; 
  out[41] += -1.224744871391589*Ghat[15]*dv11; 
  out[42] += 0.7071067811865475*Ghat[26]*dv11; 
  out[43] += -1.224744871391589*Ghat[16]*dv11; 
  out[44] += -1.224744871391589*Ghat[17]*dv11; 
  out[45] += -1.224744871391589*Ghat[18]*dv11; 
  out[46] += -1.224744871391589*Ghat[19]*dv11; 
  out[47] += 0.7071067811865475*Ghat[27]*dv11; 
  out[48] += 0.7071067811865475*Ghat[28]*dv11; 
  out[49] += 0.7071067811865475*Ghat[29]*dv11; 
  out[50] += 0.7071067811865475*Ghat[30]*dv11; 
  out[51] += -1.224744871391589*Ghat[20]*dv11; 
  out[52] += -1.224744871391589*Ghat[21]*dv11; 
  out[53] += -1.224744871391589*Ghat[22]*dv11; 
  out[54] += -1.224744871391589*Ghat[23]*dv11; 
  out[55] += -1.224744871391589*Ghat[24]*dv11; 
  out[56] += -1.224744871391589*Ghat[25]*dv11; 
  out[57] += -1.224744871391589*Ghat[26]*dv11; 
  out[58] += 0.7071067811865475*Ghat[31]*dv11; 
  out[59] += -1.224744871391589*Ghat[27]*dv11; 
  out[60] += -1.224744871391589*Ghat[28]*dv11; 
  out[61] += -1.224744871391589*Ghat[29]*dv11; 
  out[62] += -1.224744871391589*Ghat[30]*dv11; 
  out[63] += -1.224744871391589*Ghat[31]*dv11; 

  } 
} 
