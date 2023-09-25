#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_hyb_3x3v_p1_surfx4_eval_quad.h> 
#include <gkyl_basis_hyb_3x3v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_poisson_extem_boundary_surfvx_3x3v_ser_p1(const double *w, const double *dxv, const double *field, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // field:       potentials, including external (scaled by appropriate factors).
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv10 = 2/dxv[3]; 
  const double dv1 = dxv[3], wv1 = w[3]; 
  const double dv2 = dxv[4], wv2 = w[4]; 
  const double dv3 = dxv[5], wv3 = w[5]; 
  const double *phi = &field[0]; 
  const double dx10 = 2/dxv[0]; 
  const double dx11 = 2/dxv[1]; 
  const double dx12 = 2/dxv[2]; 
  const double *A0 = &field[8]; 
  const double *A1 = &field[16]; 
  const double *A2 = &field[24]; 
  double alpha[64] = {0.0}; 

  alpha[0] = (-3.464101615137754*A0[3]*dx12*wv3)+3.464101615137754*A2[1]*dx10*wv3-3.464101615137754*A0[2]*dx11*wv2+3.464101615137754*A1[1]*dx10*wv2-3.464101615137754*phi[1]*dx10; 
  alpha[1] = (-3.464101615137754*A0[5]*dx12*wv3)-3.464101615137754*A0[4]*dx11*wv2; 
  alpha[2] = (-3.464101615137754*A0[6]*dx12*wv3)+3.464101615137754*A2[4]*dx10*wv3+3.464101615137754*A1[4]*dx10*wv2-3.464101615137754*phi[4]*dx10; 
  alpha[3] = 3.464101615137754*A2[5]*dx10*wv3-3.464101615137754*A0[6]*dx11*wv2+3.464101615137754*A1[5]*dx10*wv2-3.464101615137754*phi[5]*dx10; 
  alpha[4] = A1[1]*dv2*dx10-1.0*A0[2]*dv2*dx11; 
  alpha[5] = A2[1]*dv3*dx10-1.0*A0[3]*dv3*dx12; 
  alpha[6] = -3.464101615137754*A0[7]*dx12*wv3; 
  alpha[7] = -3.464101615137754*A0[7]*dx11*wv2; 
  alpha[8] = 3.464101615137754*A2[7]*dx10*wv3+3.464101615137754*A1[7]*dx10*wv2-3.464101615137754*phi[7]*dx10; 
  alpha[9] = -1.0*A0[4]*dv2*dx11; 
  alpha[10] = A1[4]*dv2*dx10; 
  alpha[11] = A1[5]*dv2*dx10-1.0*A0[6]*dv2*dx11; 
  alpha[12] = -1.0*A0[5]*dv3*dx12; 
  alpha[13] = A2[4]*dv3*dx10-1.0*A0[6]*dv3*dx12; 
  alpha[14] = A2[5]*dv3*dx10; 
  alpha[18] = -1.0*A0[7]*dv2*dx11; 
  alpha[19] = A1[7]*dv2*dx10; 
  alpha[20] = -1.0*A0[7]*dv3*dx12; 
  alpha[22] = A2[7]*dv3*dx10; 

  double fUpwindQuad[72] = {0.0};
  double fUpwind[64] = {0.0};
  double Ghat[64] = {0.0}; 

  if (edge == -1) { 

  if ((-0.2371708245126284*(alpha[22]+alpha[20]+alpha[19]+alpha[18]))+0.2371708245126284*(alpha[14]+alpha[13]+alpha[12]+alpha[11]+alpha[10]+alpha[9])+0.1767766952966368*(alpha[8]+alpha[7]+alpha[6])-0.2371708245126284*(alpha[5]+alpha[4])-0.1767766952966368*(alpha[3]+alpha[2]+alpha[1])+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[0] = hyb_3x3v_p1_surfx4_eval_quad_node_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = hyb_3x3v_p1_surfx4_eval_quad_node_0_l(fEdge); 
  } 
  if ((-0.2371708245126284*(alpha[19]+alpha[18]))+0.2371708245126284*(alpha[11]+alpha[10]+alpha[9])+0.1767766952966368*(alpha[8]+alpha[7]+alpha[6])-0.2371708245126284*alpha[4]-0.1767766952966368*(alpha[3]+alpha[2]+alpha[1])+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[1] = hyb_3x3v_p1_surfx4_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = hyb_3x3v_p1_surfx4_eval_quad_node_1_l(fEdge); 
  } 
  if (0.2371708245126284*(alpha[22]+alpha[20])-0.2371708245126284*(alpha[19]+alpha[18]+alpha[14]+alpha[13]+alpha[12])+0.2371708245126284*(alpha[11]+alpha[10]+alpha[9])+0.1767766952966368*(alpha[8]+alpha[7]+alpha[6])+0.2371708245126284*alpha[5]-0.2371708245126284*alpha[4]-0.1767766952966368*(alpha[3]+alpha[2]+alpha[1])+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[2] = hyb_3x3v_p1_surfx4_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[2] = hyb_3x3v_p1_surfx4_eval_quad_node_2_l(fEdge); 
  } 
  if ((-0.2371708245126284*(alpha[22]+alpha[20]))+0.2371708245126284*(alpha[14]+alpha[13]+alpha[12])+0.1767766952966368*(alpha[8]+alpha[7]+alpha[6])-0.2371708245126284*alpha[5]-0.1767766952966368*(alpha[3]+alpha[2]+alpha[1])+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[3] = hyb_3x3v_p1_surfx4_eval_quad_node_3_r(fSkin); 
  } else { 
    fUpwindQuad[3] = hyb_3x3v_p1_surfx4_eval_quad_node_3_l(fEdge); 
  } 
  if (0.1767766952966368*(alpha[8]+alpha[7]+alpha[6])-0.1767766952966368*(alpha[3]+alpha[2]+alpha[1])+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[4] = hyb_3x3v_p1_surfx4_eval_quad_node_4_r(fSkin); 
  } else { 
    fUpwindQuad[4] = hyb_3x3v_p1_surfx4_eval_quad_node_4_l(fEdge); 
  } 
  if (0.2371708245126284*(alpha[22]+alpha[20])-0.2371708245126284*(alpha[14]+alpha[13]+alpha[12])+0.1767766952966368*(alpha[8]+alpha[7]+alpha[6])+0.2371708245126284*alpha[5]-0.1767766952966368*(alpha[3]+alpha[2]+alpha[1])+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[5] = hyb_3x3v_p1_surfx4_eval_quad_node_5_r(fSkin); 
  } else { 
    fUpwindQuad[5] = hyb_3x3v_p1_surfx4_eval_quad_node_5_l(fEdge); 
  } 
  if ((-0.2371708245126284*(alpha[22]+alpha[20]))+0.2371708245126284*(alpha[19]+alpha[18]+alpha[14]+alpha[13]+alpha[12])-0.2371708245126284*(alpha[11]+alpha[10]+alpha[9])+0.1767766952966368*(alpha[8]+alpha[7]+alpha[6])-0.2371708245126284*alpha[5]+0.2371708245126284*alpha[4]-0.1767766952966368*(alpha[3]+alpha[2]+alpha[1])+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[6] = hyb_3x3v_p1_surfx4_eval_quad_node_6_r(fSkin); 
  } else { 
    fUpwindQuad[6] = hyb_3x3v_p1_surfx4_eval_quad_node_6_l(fEdge); 
  } 
  if (0.2371708245126284*(alpha[19]+alpha[18])-0.2371708245126284*(alpha[11]+alpha[10]+alpha[9])+0.1767766952966368*(alpha[8]+alpha[7]+alpha[6])+0.2371708245126284*alpha[4]-0.1767766952966368*(alpha[3]+alpha[2]+alpha[1])+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[7] = hyb_3x3v_p1_surfx4_eval_quad_node_7_r(fSkin); 
  } else { 
    fUpwindQuad[7] = hyb_3x3v_p1_surfx4_eval_quad_node_7_l(fEdge); 
  } 
  if (0.2371708245126284*(alpha[22]+alpha[20]+alpha[19]+alpha[18])-0.2371708245126284*(alpha[14]+alpha[13]+alpha[12]+alpha[11]+alpha[10]+alpha[9])+0.1767766952966368*(alpha[8]+alpha[7]+alpha[6])+0.2371708245126284*(alpha[5]+alpha[4])-0.1767766952966368*(alpha[3]+alpha[2]+alpha[1])+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[8] = hyb_3x3v_p1_surfx4_eval_quad_node_8_r(fSkin); 
  } else { 
    fUpwindQuad[8] = hyb_3x3v_p1_surfx4_eval_quad_node_8_l(fEdge); 
  } 
  if (0.2371708245126284*alpha[22]-0.2371708245126284*alpha[20]+0.2371708245126284*(alpha[19]+alpha[18])-0.2371708245126284*alpha[14]+0.2371708245126284*(alpha[13]+alpha[12])-0.2371708245126284*alpha[11]+0.2371708245126284*(alpha[10]+alpha[9])-0.1767766952966368*(alpha[8]+alpha[7])+0.1767766952966368*alpha[6]-0.2371708245126284*(alpha[5]+alpha[4])+0.1767766952966368*alpha[3]-0.1767766952966368*(alpha[2]+alpha[1])+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[9] = hyb_3x3v_p1_surfx4_eval_quad_node_9_r(fSkin); 
  } else { 
    fUpwindQuad[9] = hyb_3x3v_p1_surfx4_eval_quad_node_9_l(fEdge); 
  } 
  if (0.2371708245126284*(alpha[19]+alpha[18])-0.2371708245126284*alpha[11]+0.2371708245126284*(alpha[10]+alpha[9])-0.1767766952966368*(alpha[8]+alpha[7])+0.1767766952966368*alpha[6]-0.2371708245126284*alpha[4]+0.1767766952966368*alpha[3]-0.1767766952966368*(alpha[2]+alpha[1])+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[10] = hyb_3x3v_p1_surfx4_eval_quad_node_10_r(fSkin); 
  } else { 
    fUpwindQuad[10] = hyb_3x3v_p1_surfx4_eval_quad_node_10_l(fEdge); 
  } 
  if ((-0.2371708245126284*alpha[22])+0.2371708245126284*(alpha[20]+alpha[19]+alpha[18]+alpha[14])-0.2371708245126284*(alpha[13]+alpha[12]+alpha[11])+0.2371708245126284*(alpha[10]+alpha[9])-0.1767766952966368*(alpha[8]+alpha[7])+0.1767766952966368*alpha[6]+0.2371708245126284*alpha[5]-0.2371708245126284*alpha[4]+0.1767766952966368*alpha[3]-0.1767766952966368*(alpha[2]+alpha[1])+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[11] = hyb_3x3v_p1_surfx4_eval_quad_node_11_r(fSkin); 
  } else { 
    fUpwindQuad[11] = hyb_3x3v_p1_surfx4_eval_quad_node_11_l(fEdge); 
  } 
  if (0.2371708245126284*alpha[22]-0.2371708245126284*(alpha[20]+alpha[14])+0.2371708245126284*(alpha[13]+alpha[12])-0.1767766952966368*(alpha[8]+alpha[7])+0.1767766952966368*alpha[6]-0.2371708245126284*alpha[5]+0.1767766952966368*alpha[3]-0.1767766952966368*(alpha[2]+alpha[1])+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[12] = hyb_3x3v_p1_surfx4_eval_quad_node_12_r(fSkin); 
  } else { 
    fUpwindQuad[12] = hyb_3x3v_p1_surfx4_eval_quad_node_12_l(fEdge); 
  } 
  if ((-0.1767766952966368*(alpha[8]+alpha[7]))+0.1767766952966368*(alpha[6]+alpha[3])-0.1767766952966368*(alpha[2]+alpha[1])+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[13] = hyb_3x3v_p1_surfx4_eval_quad_node_13_r(fSkin); 
  } else { 
    fUpwindQuad[13] = hyb_3x3v_p1_surfx4_eval_quad_node_13_l(fEdge); 
  } 
  if ((-0.2371708245126284*alpha[22])+0.2371708245126284*(alpha[20]+alpha[14])-0.2371708245126284*(alpha[13]+alpha[12])-0.1767766952966368*(alpha[8]+alpha[7])+0.1767766952966368*alpha[6]+0.2371708245126284*alpha[5]+0.1767766952966368*alpha[3]-0.1767766952966368*(alpha[2]+alpha[1])+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[14] = hyb_3x3v_p1_surfx4_eval_quad_node_14_r(fSkin); 
  } else { 
    fUpwindQuad[14] = hyb_3x3v_p1_surfx4_eval_quad_node_14_l(fEdge); 
  } 
  if (0.2371708245126284*alpha[22]-0.2371708245126284*(alpha[20]+alpha[19]+alpha[18]+alpha[14])+0.2371708245126284*(alpha[13]+alpha[12]+alpha[11])-0.2371708245126284*(alpha[10]+alpha[9])-0.1767766952966368*(alpha[8]+alpha[7])+0.1767766952966368*alpha[6]-0.2371708245126284*alpha[5]+0.2371708245126284*alpha[4]+0.1767766952966368*alpha[3]-0.1767766952966368*(alpha[2]+alpha[1])+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[15] = hyb_3x3v_p1_surfx4_eval_quad_node_15_r(fSkin); 
  } else { 
    fUpwindQuad[15] = hyb_3x3v_p1_surfx4_eval_quad_node_15_l(fEdge); 
  } 
  if ((-0.2371708245126284*(alpha[19]+alpha[18]))+0.2371708245126284*alpha[11]-0.2371708245126284*(alpha[10]+alpha[9])-0.1767766952966368*(alpha[8]+alpha[7])+0.1767766952966368*alpha[6]+0.2371708245126284*alpha[4]+0.1767766952966368*alpha[3]-0.1767766952966368*(alpha[2]+alpha[1])+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[16] = hyb_3x3v_p1_surfx4_eval_quad_node_16_r(fSkin); 
  } else { 
    fUpwindQuad[16] = hyb_3x3v_p1_surfx4_eval_quad_node_16_l(fEdge); 
  } 
  if ((-0.2371708245126284*alpha[22])+0.2371708245126284*alpha[20]-0.2371708245126284*(alpha[19]+alpha[18])+0.2371708245126284*alpha[14]-0.2371708245126284*(alpha[13]+alpha[12])+0.2371708245126284*alpha[11]-0.2371708245126284*(alpha[10]+alpha[9])-0.1767766952966368*(alpha[8]+alpha[7])+0.1767766952966368*alpha[6]+0.2371708245126284*(alpha[5]+alpha[4])+0.1767766952966368*alpha[3]-0.1767766952966368*(alpha[2]+alpha[1])+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[17] = hyb_3x3v_p1_surfx4_eval_quad_node_17_r(fSkin); 
  } else { 
    fUpwindQuad[17] = hyb_3x3v_p1_surfx4_eval_quad_node_17_l(fEdge); 
  } 
  if (0.2371708245126284*(alpha[22]+alpha[20]+alpha[19])-0.2371708245126284*alpha[18]+0.2371708245126284*alpha[14]-0.2371708245126284*alpha[13]+0.2371708245126284*(alpha[12]+alpha[11])-0.2371708245126284*alpha[10]+0.2371708245126284*alpha[9]-0.1767766952966368*alpha[8]+0.1767766952966368*alpha[7]-0.1767766952966368*alpha[6]-0.2371708245126284*(alpha[5]+alpha[4])-0.1767766952966368*alpha[3]+0.1767766952966368*alpha[2]-0.1767766952966368*alpha[1]+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[18] = hyb_3x3v_p1_surfx4_eval_quad_node_18_r(fSkin); 
  } else { 
    fUpwindQuad[18] = hyb_3x3v_p1_surfx4_eval_quad_node_18_l(fEdge); 
  } 
  if (0.2371708245126284*alpha[19]-0.2371708245126284*alpha[18]+0.2371708245126284*alpha[11]-0.2371708245126284*alpha[10]+0.2371708245126284*alpha[9]-0.1767766952966368*alpha[8]+0.1767766952966368*alpha[7]-0.1767766952966368*alpha[6]-0.2371708245126284*alpha[4]-0.1767766952966368*alpha[3]+0.1767766952966368*alpha[2]-0.1767766952966368*alpha[1]+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[19] = hyb_3x3v_p1_surfx4_eval_quad_node_19_r(fSkin); 
  } else { 
    fUpwindQuad[19] = hyb_3x3v_p1_surfx4_eval_quad_node_19_l(fEdge); 
  } 
  if ((-0.2371708245126284*(alpha[22]+alpha[20]))+0.2371708245126284*alpha[19]-0.2371708245126284*(alpha[18]+alpha[14])+0.2371708245126284*alpha[13]-0.2371708245126284*alpha[12]+0.2371708245126284*alpha[11]-0.2371708245126284*alpha[10]+0.2371708245126284*alpha[9]-0.1767766952966368*alpha[8]+0.1767766952966368*alpha[7]-0.1767766952966368*alpha[6]+0.2371708245126284*alpha[5]-0.2371708245126284*alpha[4]-0.1767766952966368*alpha[3]+0.1767766952966368*alpha[2]-0.1767766952966368*alpha[1]+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[20] = hyb_3x3v_p1_surfx4_eval_quad_node_20_r(fSkin); 
  } else { 
    fUpwindQuad[20] = hyb_3x3v_p1_surfx4_eval_quad_node_20_l(fEdge); 
  } 
  if (0.2371708245126284*(alpha[22]+alpha[20]+alpha[14])-0.2371708245126284*alpha[13]+0.2371708245126284*alpha[12]-0.1767766952966368*alpha[8]+0.1767766952966368*alpha[7]-0.1767766952966368*alpha[6]-0.2371708245126284*alpha[5]-0.1767766952966368*alpha[3]+0.1767766952966368*alpha[2]-0.1767766952966368*alpha[1]+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[21] = hyb_3x3v_p1_surfx4_eval_quad_node_21_r(fSkin); 
  } else { 
    fUpwindQuad[21] = hyb_3x3v_p1_surfx4_eval_quad_node_21_l(fEdge); 
  } 
  if ((-0.1767766952966368*alpha[8])+0.1767766952966368*alpha[7]-0.1767766952966368*(alpha[6]+alpha[3])+0.1767766952966368*alpha[2]-0.1767766952966368*alpha[1]+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[22] = hyb_3x3v_p1_surfx4_eval_quad_node_22_r(fSkin); 
  } else { 
    fUpwindQuad[22] = hyb_3x3v_p1_surfx4_eval_quad_node_22_l(fEdge); 
  } 
  if ((-0.2371708245126284*(alpha[22]+alpha[20]+alpha[14]))+0.2371708245126284*alpha[13]-0.2371708245126284*alpha[12]-0.1767766952966368*alpha[8]+0.1767766952966368*alpha[7]-0.1767766952966368*alpha[6]+0.2371708245126284*alpha[5]-0.1767766952966368*alpha[3]+0.1767766952966368*alpha[2]-0.1767766952966368*alpha[1]+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[23] = hyb_3x3v_p1_surfx4_eval_quad_node_23_r(fSkin); 
  } else { 
    fUpwindQuad[23] = hyb_3x3v_p1_surfx4_eval_quad_node_23_l(fEdge); 
  } 
  if (0.2371708245126284*(alpha[22]+alpha[20])-0.2371708245126284*alpha[19]+0.2371708245126284*(alpha[18]+alpha[14])-0.2371708245126284*alpha[13]+0.2371708245126284*alpha[12]-0.2371708245126284*alpha[11]+0.2371708245126284*alpha[10]-0.2371708245126284*alpha[9]-0.1767766952966368*alpha[8]+0.1767766952966368*alpha[7]-0.1767766952966368*alpha[6]-0.2371708245126284*alpha[5]+0.2371708245126284*alpha[4]-0.1767766952966368*alpha[3]+0.1767766952966368*alpha[2]-0.1767766952966368*alpha[1]+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[24] = hyb_3x3v_p1_surfx4_eval_quad_node_24_r(fSkin); 
  } else { 
    fUpwindQuad[24] = hyb_3x3v_p1_surfx4_eval_quad_node_24_l(fEdge); 
  } 
  if ((-0.2371708245126284*alpha[19])+0.2371708245126284*alpha[18]-0.2371708245126284*alpha[11]+0.2371708245126284*alpha[10]-0.2371708245126284*alpha[9]-0.1767766952966368*alpha[8]+0.1767766952966368*alpha[7]-0.1767766952966368*alpha[6]+0.2371708245126284*alpha[4]-0.1767766952966368*alpha[3]+0.1767766952966368*alpha[2]-0.1767766952966368*alpha[1]+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[25] = hyb_3x3v_p1_surfx4_eval_quad_node_25_r(fSkin); 
  } else { 
    fUpwindQuad[25] = hyb_3x3v_p1_surfx4_eval_quad_node_25_l(fEdge); 
  } 
  if ((-0.2371708245126284*(alpha[22]+alpha[20]+alpha[19]))+0.2371708245126284*alpha[18]-0.2371708245126284*alpha[14]+0.2371708245126284*alpha[13]-0.2371708245126284*(alpha[12]+alpha[11])+0.2371708245126284*alpha[10]-0.2371708245126284*alpha[9]-0.1767766952966368*alpha[8]+0.1767766952966368*alpha[7]-0.1767766952966368*alpha[6]+0.2371708245126284*(alpha[5]+alpha[4])-0.1767766952966368*alpha[3]+0.1767766952966368*alpha[2]-0.1767766952966368*alpha[1]+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[26] = hyb_3x3v_p1_surfx4_eval_quad_node_26_r(fSkin); 
  } else { 
    fUpwindQuad[26] = hyb_3x3v_p1_surfx4_eval_quad_node_26_l(fEdge); 
  } 
  if ((-0.2371708245126284*alpha[22])+0.2371708245126284*alpha[20]-0.2371708245126284*alpha[19]+0.2371708245126284*alpha[18]-0.2371708245126284*(alpha[14]+alpha[13])+0.2371708245126284*alpha[12]-0.2371708245126284*(alpha[11]+alpha[10])+0.2371708245126284*alpha[9]+0.1767766952966368*alpha[8]-0.1767766952966368*(alpha[7]+alpha[6])-0.2371708245126284*(alpha[5]+alpha[4])+0.1767766952966368*(alpha[3]+alpha[2])-0.1767766952966368*alpha[1]+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[27] = hyb_3x3v_p1_surfx4_eval_quad_node_27_r(fSkin); 
  } else { 
    fUpwindQuad[27] = hyb_3x3v_p1_surfx4_eval_quad_node_27_l(fEdge); 
  } 
  if ((-0.2371708245126284*alpha[19])+0.2371708245126284*alpha[18]-0.2371708245126284*(alpha[11]+alpha[10])+0.2371708245126284*alpha[9]+0.1767766952966368*alpha[8]-0.1767766952966368*(alpha[7]+alpha[6])-0.2371708245126284*alpha[4]+0.1767766952966368*(alpha[3]+alpha[2])-0.1767766952966368*alpha[1]+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[28] = hyb_3x3v_p1_surfx4_eval_quad_node_28_r(fSkin); 
  } else { 
    fUpwindQuad[28] = hyb_3x3v_p1_surfx4_eval_quad_node_28_l(fEdge); 
  } 
  if (0.2371708245126284*alpha[22]-0.2371708245126284*(alpha[20]+alpha[19])+0.2371708245126284*(alpha[18]+alpha[14]+alpha[13])-0.2371708245126284*(alpha[12]+alpha[11]+alpha[10])+0.2371708245126284*alpha[9]+0.1767766952966368*alpha[8]-0.1767766952966368*(alpha[7]+alpha[6])+0.2371708245126284*alpha[5]-0.2371708245126284*alpha[4]+0.1767766952966368*(alpha[3]+alpha[2])-0.1767766952966368*alpha[1]+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[29] = hyb_3x3v_p1_surfx4_eval_quad_node_29_r(fSkin); 
  } else { 
    fUpwindQuad[29] = hyb_3x3v_p1_surfx4_eval_quad_node_29_l(fEdge); 
  } 
  if ((-0.2371708245126284*alpha[22])+0.2371708245126284*alpha[20]-0.2371708245126284*(alpha[14]+alpha[13])+0.2371708245126284*alpha[12]+0.1767766952966368*alpha[8]-0.1767766952966368*(alpha[7]+alpha[6])-0.2371708245126284*alpha[5]+0.1767766952966368*(alpha[3]+alpha[2])-0.1767766952966368*alpha[1]+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[30] = hyb_3x3v_p1_surfx4_eval_quad_node_30_r(fSkin); 
  } else { 
    fUpwindQuad[30] = hyb_3x3v_p1_surfx4_eval_quad_node_30_l(fEdge); 
  } 
  if (0.1767766952966368*alpha[8]-0.1767766952966368*(alpha[7]+alpha[6])+0.1767766952966368*(alpha[3]+alpha[2])-0.1767766952966368*alpha[1]+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[31] = hyb_3x3v_p1_surfx4_eval_quad_node_31_r(fSkin); 
  } else { 
    fUpwindQuad[31] = hyb_3x3v_p1_surfx4_eval_quad_node_31_l(fEdge); 
  } 
  if (0.2371708245126284*alpha[22]-0.2371708245126284*alpha[20]+0.2371708245126284*(alpha[14]+alpha[13])-0.2371708245126284*alpha[12]+0.1767766952966368*alpha[8]-0.1767766952966368*(alpha[7]+alpha[6])+0.2371708245126284*alpha[5]+0.1767766952966368*(alpha[3]+alpha[2])-0.1767766952966368*alpha[1]+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[32] = hyb_3x3v_p1_surfx4_eval_quad_node_32_r(fSkin); 
  } else { 
    fUpwindQuad[32] = hyb_3x3v_p1_surfx4_eval_quad_node_32_l(fEdge); 
  } 
  if ((-0.2371708245126284*alpha[22])+0.2371708245126284*(alpha[20]+alpha[19])-0.2371708245126284*(alpha[18]+alpha[14]+alpha[13])+0.2371708245126284*(alpha[12]+alpha[11]+alpha[10])-0.2371708245126284*alpha[9]+0.1767766952966368*alpha[8]-0.1767766952966368*(alpha[7]+alpha[6])-0.2371708245126284*alpha[5]+0.2371708245126284*alpha[4]+0.1767766952966368*(alpha[3]+alpha[2])-0.1767766952966368*alpha[1]+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[33] = hyb_3x3v_p1_surfx4_eval_quad_node_33_r(fSkin); 
  } else { 
    fUpwindQuad[33] = hyb_3x3v_p1_surfx4_eval_quad_node_33_l(fEdge); 
  } 
  if (0.2371708245126284*alpha[19]-0.2371708245126284*alpha[18]+0.2371708245126284*(alpha[11]+alpha[10])-0.2371708245126284*alpha[9]+0.1767766952966368*alpha[8]-0.1767766952966368*(alpha[7]+alpha[6])+0.2371708245126284*alpha[4]+0.1767766952966368*(alpha[3]+alpha[2])-0.1767766952966368*alpha[1]+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[34] = hyb_3x3v_p1_surfx4_eval_quad_node_34_r(fSkin); 
  } else { 
    fUpwindQuad[34] = hyb_3x3v_p1_surfx4_eval_quad_node_34_l(fEdge); 
  } 
  if (0.2371708245126284*alpha[22]-0.2371708245126284*alpha[20]+0.2371708245126284*alpha[19]-0.2371708245126284*alpha[18]+0.2371708245126284*(alpha[14]+alpha[13])-0.2371708245126284*alpha[12]+0.2371708245126284*(alpha[11]+alpha[10])-0.2371708245126284*alpha[9]+0.1767766952966368*alpha[8]-0.1767766952966368*(alpha[7]+alpha[6])+0.2371708245126284*(alpha[5]+alpha[4])+0.1767766952966368*(alpha[3]+alpha[2])-0.1767766952966368*alpha[1]+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[35] = hyb_3x3v_p1_surfx4_eval_quad_node_35_r(fSkin); 
  } else { 
    fUpwindQuad[35] = hyb_3x3v_p1_surfx4_eval_quad_node_35_l(fEdge); 
  } 
  if ((-0.2371708245126284*alpha[22])+0.2371708245126284*alpha[20]-0.2371708245126284*alpha[19]+0.2371708245126284*(alpha[18]+alpha[14]+alpha[13])-0.2371708245126284*alpha[12]+0.2371708245126284*(alpha[11]+alpha[10])-0.2371708245126284*alpha[9]+0.1767766952966368*alpha[8]-0.1767766952966368*(alpha[7]+alpha[6])-0.2371708245126284*(alpha[5]+alpha[4])-0.1767766952966368*(alpha[3]+alpha[2])+0.1767766952966368*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[36] = hyb_3x3v_p1_surfx4_eval_quad_node_36_r(fSkin); 
  } else { 
    fUpwindQuad[36] = hyb_3x3v_p1_surfx4_eval_quad_node_36_l(fEdge); 
  } 
  if ((-0.2371708245126284*alpha[19])+0.2371708245126284*(alpha[18]+alpha[11]+alpha[10])-0.2371708245126284*alpha[9]+0.1767766952966368*alpha[8]-0.1767766952966368*(alpha[7]+alpha[6])-0.2371708245126284*alpha[4]-0.1767766952966368*(alpha[3]+alpha[2])+0.1767766952966368*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[37] = hyb_3x3v_p1_surfx4_eval_quad_node_37_r(fSkin); 
  } else { 
    fUpwindQuad[37] = hyb_3x3v_p1_surfx4_eval_quad_node_37_l(fEdge); 
  } 
  if (0.2371708245126284*alpha[22]-0.2371708245126284*(alpha[20]+alpha[19])+0.2371708245126284*alpha[18]-0.2371708245126284*(alpha[14]+alpha[13])+0.2371708245126284*(alpha[12]+alpha[11]+alpha[10])-0.2371708245126284*alpha[9]+0.1767766952966368*alpha[8]-0.1767766952966368*(alpha[7]+alpha[6])+0.2371708245126284*alpha[5]-0.2371708245126284*alpha[4]-0.1767766952966368*(alpha[3]+alpha[2])+0.1767766952966368*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[38] = hyb_3x3v_p1_surfx4_eval_quad_node_38_r(fSkin); 
  } else { 
    fUpwindQuad[38] = hyb_3x3v_p1_surfx4_eval_quad_node_38_l(fEdge); 
  } 
  if ((-0.2371708245126284*alpha[22])+0.2371708245126284*(alpha[20]+alpha[14]+alpha[13])-0.2371708245126284*alpha[12]+0.1767766952966368*alpha[8]-0.1767766952966368*(alpha[7]+alpha[6])-0.2371708245126284*alpha[5]-0.1767766952966368*(alpha[3]+alpha[2])+0.1767766952966368*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[39] = hyb_3x3v_p1_surfx4_eval_quad_node_39_r(fSkin); 
  } else { 
    fUpwindQuad[39] = hyb_3x3v_p1_surfx4_eval_quad_node_39_l(fEdge); 
  } 
  if (0.1767766952966368*alpha[8]-0.1767766952966368*(alpha[7]+alpha[6]+alpha[3]+alpha[2])+0.1767766952966368*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[40] = hyb_3x3v_p1_surfx4_eval_quad_node_40_r(fSkin); 
  } else { 
    fUpwindQuad[40] = hyb_3x3v_p1_surfx4_eval_quad_node_40_l(fEdge); 
  } 
  if (0.2371708245126284*alpha[22]-0.2371708245126284*(alpha[20]+alpha[14]+alpha[13])+0.2371708245126284*alpha[12]+0.1767766952966368*alpha[8]-0.1767766952966368*(alpha[7]+alpha[6])+0.2371708245126284*alpha[5]-0.1767766952966368*(alpha[3]+alpha[2])+0.1767766952966368*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[41] = hyb_3x3v_p1_surfx4_eval_quad_node_41_r(fSkin); 
  } else { 
    fUpwindQuad[41] = hyb_3x3v_p1_surfx4_eval_quad_node_41_l(fEdge); 
  } 
  if ((-0.2371708245126284*alpha[22])+0.2371708245126284*(alpha[20]+alpha[19])-0.2371708245126284*alpha[18]+0.2371708245126284*(alpha[14]+alpha[13])-0.2371708245126284*(alpha[12]+alpha[11]+alpha[10])+0.2371708245126284*alpha[9]+0.1767766952966368*alpha[8]-0.1767766952966368*(alpha[7]+alpha[6])-0.2371708245126284*alpha[5]+0.2371708245126284*alpha[4]-0.1767766952966368*(alpha[3]+alpha[2])+0.1767766952966368*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[42] = hyb_3x3v_p1_surfx4_eval_quad_node_42_r(fSkin); 
  } else { 
    fUpwindQuad[42] = hyb_3x3v_p1_surfx4_eval_quad_node_42_l(fEdge); 
  } 
  if (0.2371708245126284*alpha[19]-0.2371708245126284*(alpha[18]+alpha[11]+alpha[10])+0.2371708245126284*alpha[9]+0.1767766952966368*alpha[8]-0.1767766952966368*(alpha[7]+alpha[6])+0.2371708245126284*alpha[4]-0.1767766952966368*(alpha[3]+alpha[2])+0.1767766952966368*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[43] = hyb_3x3v_p1_surfx4_eval_quad_node_43_r(fSkin); 
  } else { 
    fUpwindQuad[43] = hyb_3x3v_p1_surfx4_eval_quad_node_43_l(fEdge); 
  } 
  if (0.2371708245126284*alpha[22]-0.2371708245126284*alpha[20]+0.2371708245126284*alpha[19]-0.2371708245126284*(alpha[18]+alpha[14]+alpha[13])+0.2371708245126284*alpha[12]-0.2371708245126284*(alpha[11]+alpha[10])+0.2371708245126284*alpha[9]+0.1767766952966368*alpha[8]-0.1767766952966368*(alpha[7]+alpha[6])+0.2371708245126284*(alpha[5]+alpha[4])-0.1767766952966368*(alpha[3]+alpha[2])+0.1767766952966368*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[44] = hyb_3x3v_p1_surfx4_eval_quad_node_44_r(fSkin); 
  } else { 
    fUpwindQuad[44] = hyb_3x3v_p1_surfx4_eval_quad_node_44_l(fEdge); 
  } 
  if (0.2371708245126284*(alpha[22]+alpha[20]+alpha[19])-0.2371708245126284*(alpha[18]+alpha[14])+0.2371708245126284*alpha[13]-0.2371708245126284*(alpha[12]+alpha[11])+0.2371708245126284*alpha[10]-0.2371708245126284*alpha[9]-0.1767766952966368*alpha[8]+0.1767766952966368*alpha[7]-0.1767766952966368*alpha[6]-0.2371708245126284*(alpha[5]+alpha[4])+0.1767766952966368*alpha[3]-0.1767766952966368*alpha[2]+0.1767766952966368*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[45] = hyb_3x3v_p1_surfx4_eval_quad_node_45_r(fSkin); 
  } else { 
    fUpwindQuad[45] = hyb_3x3v_p1_surfx4_eval_quad_node_45_l(fEdge); 
  } 
  if (0.2371708245126284*alpha[19]-0.2371708245126284*(alpha[18]+alpha[11])+0.2371708245126284*alpha[10]-0.2371708245126284*alpha[9]-0.1767766952966368*alpha[8]+0.1767766952966368*alpha[7]-0.1767766952966368*alpha[6]-0.2371708245126284*alpha[4]+0.1767766952966368*alpha[3]-0.1767766952966368*alpha[2]+0.1767766952966368*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[46] = hyb_3x3v_p1_surfx4_eval_quad_node_46_r(fSkin); 
  } else { 
    fUpwindQuad[46] = hyb_3x3v_p1_surfx4_eval_quad_node_46_l(fEdge); 
  } 
  if ((-0.2371708245126284*(alpha[22]+alpha[20]))+0.2371708245126284*alpha[19]-0.2371708245126284*alpha[18]+0.2371708245126284*alpha[14]-0.2371708245126284*alpha[13]+0.2371708245126284*alpha[12]-0.2371708245126284*alpha[11]+0.2371708245126284*alpha[10]-0.2371708245126284*alpha[9]-0.1767766952966368*alpha[8]+0.1767766952966368*alpha[7]-0.1767766952966368*alpha[6]+0.2371708245126284*alpha[5]-0.2371708245126284*alpha[4]+0.1767766952966368*alpha[3]-0.1767766952966368*alpha[2]+0.1767766952966368*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[47] = hyb_3x3v_p1_surfx4_eval_quad_node_47_r(fSkin); 
  } else { 
    fUpwindQuad[47] = hyb_3x3v_p1_surfx4_eval_quad_node_47_l(fEdge); 
  } 
  if (0.2371708245126284*(alpha[22]+alpha[20])-0.2371708245126284*alpha[14]+0.2371708245126284*alpha[13]-0.2371708245126284*alpha[12]-0.1767766952966368*alpha[8]+0.1767766952966368*alpha[7]-0.1767766952966368*alpha[6]-0.2371708245126284*alpha[5]+0.1767766952966368*alpha[3]-0.1767766952966368*alpha[2]+0.1767766952966368*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[48] = hyb_3x3v_p1_surfx4_eval_quad_node_48_r(fSkin); 
  } else { 
    fUpwindQuad[48] = hyb_3x3v_p1_surfx4_eval_quad_node_48_l(fEdge); 
  } 
  if ((-0.1767766952966368*alpha[8])+0.1767766952966368*alpha[7]-0.1767766952966368*alpha[6]+0.1767766952966368*alpha[3]-0.1767766952966368*alpha[2]+0.1767766952966368*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[49] = hyb_3x3v_p1_surfx4_eval_quad_node_49_r(fSkin); 
  } else { 
    fUpwindQuad[49] = hyb_3x3v_p1_surfx4_eval_quad_node_49_l(fEdge); 
  } 
  if ((-0.2371708245126284*(alpha[22]+alpha[20]))+0.2371708245126284*alpha[14]-0.2371708245126284*alpha[13]+0.2371708245126284*alpha[12]-0.1767766952966368*alpha[8]+0.1767766952966368*alpha[7]-0.1767766952966368*alpha[6]+0.2371708245126284*alpha[5]+0.1767766952966368*alpha[3]-0.1767766952966368*alpha[2]+0.1767766952966368*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[50] = hyb_3x3v_p1_surfx4_eval_quad_node_50_r(fSkin); 
  } else { 
    fUpwindQuad[50] = hyb_3x3v_p1_surfx4_eval_quad_node_50_l(fEdge); 
  } 
  if (0.2371708245126284*(alpha[22]+alpha[20])-0.2371708245126284*alpha[19]+0.2371708245126284*alpha[18]-0.2371708245126284*alpha[14]+0.2371708245126284*alpha[13]-0.2371708245126284*alpha[12]+0.2371708245126284*alpha[11]-0.2371708245126284*alpha[10]+0.2371708245126284*alpha[9]-0.1767766952966368*alpha[8]+0.1767766952966368*alpha[7]-0.1767766952966368*alpha[6]-0.2371708245126284*alpha[5]+0.2371708245126284*alpha[4]+0.1767766952966368*alpha[3]-0.1767766952966368*alpha[2]+0.1767766952966368*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[51] = hyb_3x3v_p1_surfx4_eval_quad_node_51_r(fSkin); 
  } else { 
    fUpwindQuad[51] = hyb_3x3v_p1_surfx4_eval_quad_node_51_l(fEdge); 
  } 
  if ((-0.2371708245126284*alpha[19])+0.2371708245126284*(alpha[18]+alpha[11])-0.2371708245126284*alpha[10]+0.2371708245126284*alpha[9]-0.1767766952966368*alpha[8]+0.1767766952966368*alpha[7]-0.1767766952966368*alpha[6]+0.2371708245126284*alpha[4]+0.1767766952966368*alpha[3]-0.1767766952966368*alpha[2]+0.1767766952966368*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[52] = hyb_3x3v_p1_surfx4_eval_quad_node_52_r(fSkin); 
  } else { 
    fUpwindQuad[52] = hyb_3x3v_p1_surfx4_eval_quad_node_52_l(fEdge); 
  } 
  if ((-0.2371708245126284*(alpha[22]+alpha[20]+alpha[19]))+0.2371708245126284*(alpha[18]+alpha[14])-0.2371708245126284*alpha[13]+0.2371708245126284*(alpha[12]+alpha[11])-0.2371708245126284*alpha[10]+0.2371708245126284*alpha[9]-0.1767766952966368*alpha[8]+0.1767766952966368*alpha[7]-0.1767766952966368*alpha[6]+0.2371708245126284*(alpha[5]+alpha[4])+0.1767766952966368*alpha[3]-0.1767766952966368*alpha[2]+0.1767766952966368*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[53] = hyb_3x3v_p1_surfx4_eval_quad_node_53_r(fSkin); 
  } else { 
    fUpwindQuad[53] = hyb_3x3v_p1_surfx4_eval_quad_node_53_l(fEdge); 
  } 
  if (0.2371708245126284*alpha[22]-0.2371708245126284*alpha[20]+0.2371708245126284*(alpha[19]+alpha[18]+alpha[14])-0.2371708245126284*(alpha[13]+alpha[12])+0.2371708245126284*alpha[11]-0.2371708245126284*(alpha[10]+alpha[9])-0.1767766952966368*(alpha[8]+alpha[7])+0.1767766952966368*alpha[6]-0.2371708245126284*(alpha[5]+alpha[4])-0.1767766952966368*alpha[3]+0.1767766952966368*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[54] = hyb_3x3v_p1_surfx4_eval_quad_node_54_r(fSkin); 
  } else { 
    fUpwindQuad[54] = hyb_3x3v_p1_surfx4_eval_quad_node_54_l(fEdge); 
  } 
  if (0.2371708245126284*(alpha[19]+alpha[18]+alpha[11])-0.2371708245126284*(alpha[10]+alpha[9])-0.1767766952966368*(alpha[8]+alpha[7])+0.1767766952966368*alpha[6]-0.2371708245126284*alpha[4]-0.1767766952966368*alpha[3]+0.1767766952966368*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[55] = hyb_3x3v_p1_surfx4_eval_quad_node_55_r(fSkin); 
  } else { 
    fUpwindQuad[55] = hyb_3x3v_p1_surfx4_eval_quad_node_55_l(fEdge); 
  } 
  if ((-0.2371708245126284*alpha[22])+0.2371708245126284*(alpha[20]+alpha[19]+alpha[18])-0.2371708245126284*alpha[14]+0.2371708245126284*(alpha[13]+alpha[12]+alpha[11])-0.2371708245126284*(alpha[10]+alpha[9])-0.1767766952966368*(alpha[8]+alpha[7])+0.1767766952966368*alpha[6]+0.2371708245126284*alpha[5]-0.2371708245126284*alpha[4]-0.1767766952966368*alpha[3]+0.1767766952966368*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[56] = hyb_3x3v_p1_surfx4_eval_quad_node_56_r(fSkin); 
  } else { 
    fUpwindQuad[56] = hyb_3x3v_p1_surfx4_eval_quad_node_56_l(fEdge); 
  } 
  if (0.2371708245126284*alpha[22]-0.2371708245126284*alpha[20]+0.2371708245126284*alpha[14]-0.2371708245126284*(alpha[13]+alpha[12])-0.1767766952966368*(alpha[8]+alpha[7])+0.1767766952966368*alpha[6]-0.2371708245126284*alpha[5]-0.1767766952966368*alpha[3]+0.1767766952966368*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[57] = hyb_3x3v_p1_surfx4_eval_quad_node_57_r(fSkin); 
  } else { 
    fUpwindQuad[57] = hyb_3x3v_p1_surfx4_eval_quad_node_57_l(fEdge); 
  } 
  if ((-0.1767766952966368*(alpha[8]+alpha[7]))+0.1767766952966368*alpha[6]-0.1767766952966368*alpha[3]+0.1767766952966368*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[58] = hyb_3x3v_p1_surfx4_eval_quad_node_58_r(fSkin); 
  } else { 
    fUpwindQuad[58] = hyb_3x3v_p1_surfx4_eval_quad_node_58_l(fEdge); 
  } 
  if ((-0.2371708245126284*alpha[22])+0.2371708245126284*alpha[20]-0.2371708245126284*alpha[14]+0.2371708245126284*(alpha[13]+alpha[12])-0.1767766952966368*(alpha[8]+alpha[7])+0.1767766952966368*alpha[6]+0.2371708245126284*alpha[5]-0.1767766952966368*alpha[3]+0.1767766952966368*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[59] = hyb_3x3v_p1_surfx4_eval_quad_node_59_r(fSkin); 
  } else { 
    fUpwindQuad[59] = hyb_3x3v_p1_surfx4_eval_quad_node_59_l(fEdge); 
  } 
  if (0.2371708245126284*alpha[22]-0.2371708245126284*(alpha[20]+alpha[19]+alpha[18])+0.2371708245126284*alpha[14]-0.2371708245126284*(alpha[13]+alpha[12]+alpha[11])+0.2371708245126284*(alpha[10]+alpha[9])-0.1767766952966368*(alpha[8]+alpha[7])+0.1767766952966368*alpha[6]-0.2371708245126284*alpha[5]+0.2371708245126284*alpha[4]-0.1767766952966368*alpha[3]+0.1767766952966368*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[60] = hyb_3x3v_p1_surfx4_eval_quad_node_60_r(fSkin); 
  } else { 
    fUpwindQuad[60] = hyb_3x3v_p1_surfx4_eval_quad_node_60_l(fEdge); 
  } 
  if ((-0.2371708245126284*(alpha[19]+alpha[18]+alpha[11]))+0.2371708245126284*(alpha[10]+alpha[9])-0.1767766952966368*(alpha[8]+alpha[7])+0.1767766952966368*alpha[6]+0.2371708245126284*alpha[4]-0.1767766952966368*alpha[3]+0.1767766952966368*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[61] = hyb_3x3v_p1_surfx4_eval_quad_node_61_r(fSkin); 
  } else { 
    fUpwindQuad[61] = hyb_3x3v_p1_surfx4_eval_quad_node_61_l(fEdge); 
  } 
  if ((-0.2371708245126284*alpha[22])+0.2371708245126284*alpha[20]-0.2371708245126284*(alpha[19]+alpha[18]+alpha[14])+0.2371708245126284*(alpha[13]+alpha[12])-0.2371708245126284*alpha[11]+0.2371708245126284*(alpha[10]+alpha[9])-0.1767766952966368*(alpha[8]+alpha[7])+0.1767766952966368*alpha[6]+0.2371708245126284*(alpha[5]+alpha[4])-0.1767766952966368*alpha[3]+0.1767766952966368*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[62] = hyb_3x3v_p1_surfx4_eval_quad_node_62_r(fSkin); 
  } else { 
    fUpwindQuad[62] = hyb_3x3v_p1_surfx4_eval_quad_node_62_l(fEdge); 
  } 
  if ((-0.2371708245126284*(alpha[22]+alpha[20]+alpha[19]+alpha[18]+alpha[14]+alpha[13]+alpha[12]+alpha[11]+alpha[10]+alpha[9]))+0.1767766952966368*(alpha[8]+alpha[7]+alpha[6])-0.2371708245126284*(alpha[5]+alpha[4])+0.1767766952966368*(alpha[3]+alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[63] = hyb_3x3v_p1_surfx4_eval_quad_node_63_r(fSkin); 
  } else { 
    fUpwindQuad[63] = hyb_3x3v_p1_surfx4_eval_quad_node_63_l(fEdge); 
  } 
  if ((-0.2371708245126284*(alpha[19]+alpha[18]+alpha[11]+alpha[10]+alpha[9]))+0.1767766952966368*(alpha[8]+alpha[7]+alpha[6])-0.2371708245126284*alpha[4]+0.1767766952966368*(alpha[3]+alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[64] = hyb_3x3v_p1_surfx4_eval_quad_node_64_r(fSkin); 
  } else { 
    fUpwindQuad[64] = hyb_3x3v_p1_surfx4_eval_quad_node_64_l(fEdge); 
  } 
  if (0.2371708245126284*(alpha[22]+alpha[20])-0.2371708245126284*(alpha[19]+alpha[18])+0.2371708245126284*(alpha[14]+alpha[13]+alpha[12])-0.2371708245126284*(alpha[11]+alpha[10]+alpha[9])+0.1767766952966368*(alpha[8]+alpha[7]+alpha[6])+0.2371708245126284*alpha[5]-0.2371708245126284*alpha[4]+0.1767766952966368*(alpha[3]+alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[65] = hyb_3x3v_p1_surfx4_eval_quad_node_65_r(fSkin); 
  } else { 
    fUpwindQuad[65] = hyb_3x3v_p1_surfx4_eval_quad_node_65_l(fEdge); 
  } 
  if ((-0.2371708245126284*(alpha[22]+alpha[20]+alpha[14]+alpha[13]+alpha[12]))+0.1767766952966368*(alpha[8]+alpha[7]+alpha[6])-0.2371708245126284*alpha[5]+0.1767766952966368*(alpha[3]+alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[66] = hyb_3x3v_p1_surfx4_eval_quad_node_66_r(fSkin); 
  } else { 
    fUpwindQuad[66] = hyb_3x3v_p1_surfx4_eval_quad_node_66_l(fEdge); 
  } 
  if (0.1767766952966368*(alpha[8]+alpha[7]+alpha[6]+alpha[3]+alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[67] = hyb_3x3v_p1_surfx4_eval_quad_node_67_r(fSkin); 
  } else { 
    fUpwindQuad[67] = hyb_3x3v_p1_surfx4_eval_quad_node_67_l(fEdge); 
  } 
  if (0.2371708245126284*(alpha[22]+alpha[20]+alpha[14]+alpha[13]+alpha[12])+0.1767766952966368*(alpha[8]+alpha[7]+alpha[6])+0.2371708245126284*alpha[5]+0.1767766952966368*(alpha[3]+alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[68] = hyb_3x3v_p1_surfx4_eval_quad_node_68_r(fSkin); 
  } else { 
    fUpwindQuad[68] = hyb_3x3v_p1_surfx4_eval_quad_node_68_l(fEdge); 
  } 
  if ((-0.2371708245126284*(alpha[22]+alpha[20]))+0.2371708245126284*(alpha[19]+alpha[18])-0.2371708245126284*(alpha[14]+alpha[13]+alpha[12])+0.2371708245126284*(alpha[11]+alpha[10]+alpha[9])+0.1767766952966368*(alpha[8]+alpha[7]+alpha[6])-0.2371708245126284*alpha[5]+0.2371708245126284*alpha[4]+0.1767766952966368*(alpha[3]+alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[69] = hyb_3x3v_p1_surfx4_eval_quad_node_69_r(fSkin); 
  } else { 
    fUpwindQuad[69] = hyb_3x3v_p1_surfx4_eval_quad_node_69_l(fEdge); 
  } 
  if (0.2371708245126284*(alpha[19]+alpha[18]+alpha[11]+alpha[10]+alpha[9])+0.1767766952966368*(alpha[8]+alpha[7]+alpha[6])+0.2371708245126284*alpha[4]+0.1767766952966368*(alpha[3]+alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[70] = hyb_3x3v_p1_surfx4_eval_quad_node_70_r(fSkin); 
  } else { 
    fUpwindQuad[70] = hyb_3x3v_p1_surfx4_eval_quad_node_70_l(fEdge); 
  } 
  if (0.2371708245126284*(alpha[22]+alpha[20]+alpha[19]+alpha[18]+alpha[14]+alpha[13]+alpha[12]+alpha[11]+alpha[10]+alpha[9])+0.1767766952966368*(alpha[8]+alpha[7]+alpha[6])+0.2371708245126284*(alpha[5]+alpha[4])+0.1767766952966368*(alpha[3]+alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[71] = hyb_3x3v_p1_surfx4_eval_quad_node_71_r(fSkin); 
  } else { 
    fUpwindQuad[71] = hyb_3x3v_p1_surfx4_eval_quad_node_71_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_3x3v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.1767766952966368*alpha[22]*fUpwind[22]+0.1767766952966368*alpha[20]*fUpwind[20]+0.1767766952966368*alpha[19]*fUpwind[19]+0.1767766952966368*alpha[18]*fUpwind[18]+0.1767766952966368*alpha[14]*fUpwind[14]+0.1767766952966368*alpha[13]*fUpwind[13]+0.1767766952966368*alpha[12]*fUpwind[12]+0.1767766952966368*alpha[11]*fUpwind[11]+0.1767766952966368*alpha[10]*fUpwind[10]+0.1767766952966368*alpha[9]*fUpwind[9]+0.1767766952966368*alpha[8]*fUpwind[8]+0.1767766952966368*alpha[7]*fUpwind[7]+0.1767766952966368*alpha[6]*fUpwind[6]+0.1767766952966368*alpha[5]*fUpwind[5]+0.1767766952966368*alpha[4]*fUpwind[4]+0.1767766952966368*alpha[3]*fUpwind[3]+0.1767766952966368*alpha[2]*fUpwind[2]+0.1767766952966368*alpha[1]*fUpwind[1]+0.1767766952966368*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.1767766952966368*alpha[22]*fUpwind[27]+0.1767766952966368*alpha[19]*fUpwind[26]+0.1767766952966368*alpha[14]*fUpwind[21]+0.1767766952966368*alpha[13]*fUpwind[20]+0.1767766952966368*fUpwind[13]*alpha[20]+0.1767766952966368*alpha[11]*fUpwind[18]+0.1767766952966368*fUpwind[11]*alpha[18]+0.1767766952966368*alpha[10]*fUpwind[17]+0.1767766952966368*alpha[8]*fUpwind[16]+0.1767766952966368*alpha[5]*fUpwind[12]+0.1767766952966368*fUpwind[5]*alpha[12]+0.1767766952966368*alpha[4]*fUpwind[9]+0.1767766952966368*fUpwind[4]*alpha[9]+0.1767766952966368*alpha[3]*fUpwind[7]+0.1767766952966368*fUpwind[3]*alpha[7]+0.1767766952966368*alpha[2]*fUpwind[6]+0.1767766952966368*fUpwind[2]*alpha[6]+0.1767766952966368*alpha[0]*fUpwind[1]+0.1767766952966368*fUpwind[0]*alpha[1]; 
  Ghat[2] = 0.1767766952966368*alpha[18]*fUpwind[26]+0.1767766952966368*alpha[14]*fUpwind[22]+0.1767766952966368*fUpwind[14]*alpha[22]+0.1767766952966368*alpha[12]*fUpwind[20]+0.1767766952966368*fUpwind[12]*alpha[20]+0.1767766952966368*alpha[11]*fUpwind[19]+0.1767766952966368*fUpwind[11]*alpha[19]+0.1767766952966368*alpha[9]*fUpwind[17]+0.1767766952966368*alpha[7]*fUpwind[16]+0.1767766952966368*alpha[5]*fUpwind[13]+0.1767766952966368*fUpwind[5]*alpha[13]+0.1767766952966368*alpha[4]*fUpwind[10]+0.1767766952966368*fUpwind[4]*alpha[10]+0.1767766952966368*alpha[3]*fUpwind[8]+0.1767766952966368*fUpwind[3]*alpha[8]+0.1767766952966368*alpha[1]*fUpwind[6]+0.1767766952966368*fUpwind[1]*alpha[6]+0.1767766952966368*alpha[0]*fUpwind[2]+0.1767766952966368*fUpwind[0]*alpha[2]; 
  Ghat[3] = 0.1767766952966368*alpha[20]*fUpwind[27]+0.1767766952966368*alpha[13]*fUpwind[22]+0.1767766952966368*fUpwind[13]*alpha[22]+0.1767766952966368*alpha[12]*fUpwind[21]+0.1767766952966368*alpha[10]*fUpwind[19]+0.1767766952966368*fUpwind[10]*alpha[19]+0.1767766952966368*alpha[9]*fUpwind[18]+0.1767766952966368*fUpwind[9]*alpha[18]+0.1767766952966368*alpha[6]*fUpwind[16]+0.1767766952966368*alpha[5]*fUpwind[14]+0.1767766952966368*fUpwind[5]*alpha[14]+0.1767766952966368*alpha[4]*fUpwind[11]+0.1767766952966368*fUpwind[4]*alpha[11]+0.1767766952966368*alpha[2]*fUpwind[8]+0.1767766952966368*fUpwind[2]*alpha[8]+0.1767766952966368*alpha[1]*fUpwind[7]+0.1767766952966368*fUpwind[1]*alpha[7]+0.1767766952966368*alpha[0]*fUpwind[3]+0.1767766952966368*fUpwind[0]*alpha[3]; 
  Ghat[4] = 0.1581138830084189*alpha[19]*fUpwind[39]+0.1581138830084189*alpha[18]*fUpwind[38]+0.1581138830084189*alpha[11]*fUpwind[35]+0.1581138830084189*alpha[10]*fUpwind[34]+0.1581138830084189*alpha[9]*fUpwind[33]+0.1581138830084189*alpha[4]*fUpwind[32]+0.1767766952966368*alpha[22]*fUpwind[30]+0.1767766952966368*alpha[20]*fUpwind[28]+0.1767766952966368*alpha[14]*fUpwind[25]+0.1767766952966368*alpha[13]*fUpwind[24]+0.1767766952966368*alpha[12]*fUpwind[23]+0.1767766952966368*alpha[8]*fUpwind[19]+0.1767766952966368*fUpwind[8]*alpha[19]+0.1767766952966368*alpha[7]*fUpwind[18]+0.1767766952966368*fUpwind[7]*alpha[18]+0.1767766952966368*alpha[6]*fUpwind[17]+0.1767766952966368*alpha[5]*fUpwind[15]+0.1767766952966368*alpha[3]*fUpwind[11]+0.1767766952966368*fUpwind[3]*alpha[11]+0.1767766952966368*alpha[2]*fUpwind[10]+0.1767766952966368*fUpwind[2]*alpha[10]+0.1767766952966368*alpha[1]*fUpwind[9]+0.1767766952966368*fUpwind[1]*alpha[9]+0.1767766952966368*alpha[0]*fUpwind[4]+0.1767766952966368*fUpwind[0]*alpha[4]; 
  Ghat[5] = 0.1581138830084189*alpha[22]*fUpwind[55]+0.1581138830084189*alpha[20]*fUpwind[53]+0.1581138830084189*alpha[14]*fUpwind[51]+0.1581138830084189*alpha[13]*fUpwind[50]+0.1581138830084189*alpha[12]*fUpwind[49]+0.1581138830084189*alpha[5]*fUpwind[48]+0.1767766952966368*alpha[19]*fUpwind[30]+0.1767766952966368*alpha[18]*fUpwind[29]+0.1767766952966368*alpha[11]*fUpwind[25]+0.1767766952966368*alpha[10]*fUpwind[24]+0.1767766952966368*alpha[9]*fUpwind[23]+0.1767766952966368*alpha[8]*fUpwind[22]+0.1767766952966368*fUpwind[8]*alpha[22]+0.1767766952966368*alpha[7]*fUpwind[21]+0.1767766952966368*alpha[6]*fUpwind[20]+0.1767766952966368*fUpwind[6]*alpha[20]+0.1767766952966368*alpha[4]*fUpwind[15]+0.1767766952966368*alpha[3]*fUpwind[14]+0.1767766952966368*fUpwind[3]*alpha[14]+0.1767766952966368*alpha[2]*fUpwind[13]+0.1767766952966368*fUpwind[2]*alpha[13]+0.1767766952966368*alpha[1]*fUpwind[12]+0.1767766952966368*fUpwind[1]*alpha[12]+0.1767766952966368*alpha[0]*fUpwind[5]+0.1767766952966368*fUpwind[0]*alpha[5]; 
  Ghat[6] = 0.1767766952966368*alpha[14]*fUpwind[27]+0.1767766952966368*alpha[11]*fUpwind[26]+0.1767766952966368*fUpwind[21]*alpha[22]+0.1767766952966368*alpha[5]*fUpwind[20]+0.1767766952966368*fUpwind[5]*alpha[20]+0.1767766952966368*alpha[18]*fUpwind[19]+0.1767766952966368*fUpwind[18]*alpha[19]+0.1767766952966368*alpha[4]*fUpwind[17]+0.1767766952966368*alpha[3]*fUpwind[16]+0.1767766952966368*alpha[12]*fUpwind[13]+0.1767766952966368*fUpwind[12]*alpha[13]+0.1767766952966368*alpha[9]*fUpwind[10]+0.1767766952966368*fUpwind[9]*alpha[10]+0.1767766952966368*alpha[7]*fUpwind[8]+0.1767766952966368*fUpwind[7]*alpha[8]+0.1767766952966368*alpha[0]*fUpwind[6]+0.1767766952966368*fUpwind[0]*alpha[6]+0.1767766952966368*alpha[1]*fUpwind[2]+0.1767766952966368*fUpwind[1]*alpha[2]; 
  Ghat[7] = 0.1767766952966368*alpha[13]*fUpwind[27]+0.1767766952966368*alpha[10]*fUpwind[26]+0.1767766952966368*alpha[20]*fUpwind[22]+0.1767766952966368*fUpwind[20]*alpha[22]+0.1767766952966368*alpha[5]*fUpwind[21]+0.1767766952966368*fUpwind[17]*alpha[19]+0.1767766952966368*alpha[4]*fUpwind[18]+0.1767766952966368*fUpwind[4]*alpha[18]+0.1767766952966368*alpha[2]*fUpwind[16]+0.1767766952966368*alpha[12]*fUpwind[14]+0.1767766952966368*fUpwind[12]*alpha[14]+0.1767766952966368*alpha[9]*fUpwind[11]+0.1767766952966368*fUpwind[9]*alpha[11]+0.1767766952966368*alpha[6]*fUpwind[8]+0.1767766952966368*fUpwind[6]*alpha[8]+0.1767766952966368*alpha[0]*fUpwind[7]+0.1767766952966368*fUpwind[0]*alpha[7]+0.1767766952966368*alpha[1]*fUpwind[3]+0.1767766952966368*fUpwind[1]*alpha[3]; 
  Ghat[8] = 0.1767766952966368*alpha[12]*fUpwind[27]+0.1767766952966368*alpha[9]*fUpwind[26]+0.1767766952966368*alpha[5]*fUpwind[22]+0.1767766952966368*fUpwind[5]*alpha[22]+0.1767766952966368*alpha[20]*fUpwind[21]+0.1767766952966368*alpha[4]*fUpwind[19]+0.1767766952966368*fUpwind[4]*alpha[19]+0.1767766952966368*fUpwind[17]*alpha[18]+0.1767766952966368*alpha[1]*fUpwind[16]+0.1767766952966368*alpha[13]*fUpwind[14]+0.1767766952966368*fUpwind[13]*alpha[14]+0.1767766952966368*alpha[10]*fUpwind[11]+0.1767766952966368*fUpwind[10]*alpha[11]+0.1767766952966368*alpha[0]*fUpwind[8]+0.1767766952966368*fUpwind[0]*alpha[8]+0.1767766952966368*alpha[6]*fUpwind[7]+0.1767766952966368*fUpwind[6]*alpha[7]+0.1767766952966368*alpha[2]*fUpwind[3]+0.1767766952966368*fUpwind[2]*alpha[3]; 
  Ghat[9] = 0.1581138830084189*alpha[19]*fUpwind[43]+0.1581138830084189*alpha[11]*fUpwind[38]+0.1581138830084189*alpha[10]*fUpwind[37]+0.1581138830084189*alpha[18]*fUpwind[35]+0.1581138830084189*alpha[4]*fUpwind[33]+0.1581138830084189*alpha[9]*fUpwind[32]+0.1767766952966368*alpha[22]*fUpwind[31]+0.1767766952966368*alpha[14]*fUpwind[29]+0.1767766952966368*alpha[13]*fUpwind[28]+0.1767766952966368*alpha[8]*fUpwind[26]+0.1767766952966368*alpha[20]*fUpwind[24]+0.1767766952966368*alpha[5]*fUpwind[23]+0.1767766952966368*fUpwind[16]*alpha[19]+0.1767766952966368*alpha[3]*fUpwind[18]+0.1767766952966368*fUpwind[3]*alpha[18]+0.1767766952966368*alpha[2]*fUpwind[17]+0.1767766952966368*alpha[12]*fUpwind[15]+0.1767766952966368*alpha[7]*fUpwind[11]+0.1767766952966368*fUpwind[7]*alpha[11]+0.1767766952966368*alpha[6]*fUpwind[10]+0.1767766952966368*fUpwind[6]*alpha[10]+0.1767766952966368*alpha[0]*fUpwind[9]+0.1767766952966368*fUpwind[0]*alpha[9]+0.1767766952966368*alpha[1]*fUpwind[4]+0.1767766952966368*fUpwind[1]*alpha[4]; 
  Ghat[10] = 0.1581138830084189*alpha[18]*fUpwind[43]+0.1581138830084189*alpha[11]*fUpwind[39]+0.1581138830084189*alpha[9]*fUpwind[37]+0.1581138830084189*alpha[19]*fUpwind[35]+0.1581138830084189*alpha[4]*fUpwind[34]+0.1581138830084189*alpha[10]*fUpwind[32]+0.1767766952966368*alpha[14]*fUpwind[30]+0.1767766952966368*alpha[12]*fUpwind[28]+0.1767766952966368*alpha[7]*fUpwind[26]+0.1767766952966368*alpha[22]*fUpwind[25]+0.1767766952966368*alpha[5]*fUpwind[24]+0.1767766952966368*alpha[20]*fUpwind[23]+0.1767766952966368*alpha[3]*fUpwind[19]+0.1767766952966368*fUpwind[3]*alpha[19]+0.1767766952966368*fUpwind[16]*alpha[18]+0.1767766952966368*alpha[1]*fUpwind[17]+0.1767766952966368*alpha[13]*fUpwind[15]+0.1767766952966368*alpha[8]*fUpwind[11]+0.1767766952966368*fUpwind[8]*alpha[11]+0.1767766952966368*alpha[0]*fUpwind[10]+0.1767766952966368*fUpwind[0]*alpha[10]+0.1767766952966368*alpha[6]*fUpwind[9]+0.1767766952966368*fUpwind[6]*alpha[9]+0.1767766952966368*alpha[2]*fUpwind[4]+0.1767766952966368*fUpwind[2]*alpha[4]; 
  Ghat[11] = 0.1581138830084189*alpha[10]*fUpwind[39]+0.1581138830084189*alpha[9]*fUpwind[38]+0.1581138830084189*alpha[4]*fUpwind[35]+0.1581138830084189*alpha[19]*fUpwind[34]+0.1581138830084189*alpha[18]*fUpwind[33]+0.1581138830084189*alpha[11]*fUpwind[32]+0.1767766952966368*alpha[20]*fUpwind[31]+0.1767766952966368*alpha[13]*fUpwind[30]+0.1767766952966368*alpha[12]*fUpwind[29]+0.1767766952966368*alpha[6]*fUpwind[26]+0.1767766952966368*alpha[5]*fUpwind[25]+0.1767766952966368*alpha[22]*fUpwind[24]+0.1767766952966368*alpha[2]*fUpwind[19]+0.1767766952966368*fUpwind[2]*alpha[19]+0.1767766952966368*alpha[1]*fUpwind[18]+0.1767766952966368*fUpwind[1]*alpha[18]+0.1767766952966368*alpha[14]*fUpwind[15]+0.1767766952966368*alpha[0]*fUpwind[11]+0.1767766952966368*fUpwind[0]*alpha[11]+0.1767766952966368*alpha[8]*fUpwind[10]+0.1767766952966368*fUpwind[8]*alpha[10]+0.1767766952966368*alpha[7]*fUpwind[9]+0.1767766952966368*fUpwind[7]*alpha[9]+0.1767766952966368*alpha[3]*fUpwind[4]+0.1767766952966368*fUpwind[3]*alpha[4]; 
  Ghat[12] = 0.1581138830084189*alpha[22]*fUpwind[59]+0.1581138830084189*alpha[14]*fUpwind[54]+0.1581138830084189*alpha[13]*fUpwind[53]+0.1581138830084189*alpha[20]*fUpwind[50]+0.1581138830084189*alpha[5]*fUpwind[49]+0.1581138830084189*alpha[12]*fUpwind[48]+0.1767766952966368*alpha[19]*fUpwind[31]+0.1767766952966368*alpha[11]*fUpwind[29]+0.1767766952966368*alpha[10]*fUpwind[28]+0.1767766952966368*alpha[8]*fUpwind[27]+0.1767766952966368*alpha[18]*fUpwind[25]+0.1767766952966368*alpha[4]*fUpwind[23]+0.1767766952966368*fUpwind[16]*alpha[22]+0.1767766952966368*alpha[3]*fUpwind[21]+0.1767766952966368*alpha[2]*fUpwind[20]+0.1767766952966368*fUpwind[2]*alpha[20]+0.1767766952966368*alpha[9]*fUpwind[15]+0.1767766952966368*alpha[7]*fUpwind[14]+0.1767766952966368*fUpwind[7]*alpha[14]+0.1767766952966368*alpha[6]*fUpwind[13]+0.1767766952966368*fUpwind[6]*alpha[13]+0.1767766952966368*alpha[0]*fUpwind[12]+0.1767766952966368*fUpwind[0]*alpha[12]+0.1767766952966368*alpha[1]*fUpwind[5]+0.1767766952966368*fUpwind[1]*alpha[5]; 
  Ghat[13] = 0.1581138830084189*alpha[14]*fUpwind[55]+0.1581138830084189*alpha[12]*fUpwind[53]+0.1581138830084189*alpha[22]*fUpwind[51]+0.1581138830084189*alpha[5]*fUpwind[50]+0.1581138830084189*alpha[20]*fUpwind[49]+0.1581138830084189*alpha[13]*fUpwind[48]+0.1767766952966368*alpha[18]*fUpwind[31]+0.1767766952966368*alpha[11]*fUpwind[30]+0.1767766952966368*alpha[9]*fUpwind[28]+0.1767766952966368*alpha[7]*fUpwind[27]+0.1767766952966368*alpha[19]*fUpwind[25]+0.1767766952966368*alpha[4]*fUpwind[24]+0.1767766952966368*alpha[3]*fUpwind[22]+0.1767766952966368*fUpwind[3]*alpha[22]+0.1767766952966368*alpha[1]*fUpwind[20]+0.1767766952966368*fUpwind[1]*alpha[20]+0.1767766952966368*alpha[10]*fUpwind[15]+0.1767766952966368*alpha[8]*fUpwind[14]+0.1767766952966368*fUpwind[8]*alpha[14]+0.1767766952966368*alpha[0]*fUpwind[13]+0.1767766952966368*fUpwind[0]*alpha[13]+0.1767766952966368*alpha[6]*fUpwind[12]+0.1767766952966368*fUpwind[6]*alpha[12]+0.1767766952966368*alpha[2]*fUpwind[5]+0.1767766952966368*fUpwind[2]*alpha[5]; 
  Ghat[14] = 0.1581138830084189*alpha[20]*fUpwind[59]+0.1581138830084189*alpha[13]*fUpwind[55]+0.1581138830084189*alpha[12]*fUpwind[54]+0.1581138830084189*alpha[5]*fUpwind[51]+0.1581138830084189*alpha[22]*fUpwind[50]+0.1581138830084189*alpha[14]*fUpwind[48]+0.1767766952966368*alpha[10]*fUpwind[30]+0.1767766952966368*alpha[9]*fUpwind[29]+0.1767766952966368*alpha[6]*fUpwind[27]+0.1767766952966368*alpha[4]*fUpwind[25]+0.1767766952966368*alpha[19]*fUpwind[24]+0.1767766952966368*alpha[18]*fUpwind[23]+0.1767766952966368*alpha[2]*fUpwind[22]+0.1767766952966368*fUpwind[2]*alpha[22]+0.1767766952966368*alpha[1]*fUpwind[21]+0.1767766952966368*fUpwind[16]*alpha[20]+0.1767766952966368*alpha[11]*fUpwind[15]+0.1767766952966368*alpha[0]*fUpwind[14]+0.1767766952966368*fUpwind[0]*alpha[14]+0.1767766952966368*alpha[8]*fUpwind[13]+0.1767766952966368*fUpwind[8]*alpha[13]+0.1767766952966368*alpha[7]*fUpwind[12]+0.1767766952966368*fUpwind[7]*alpha[12]+0.1767766952966368*alpha[3]*fUpwind[5]+0.1767766952966368*fUpwind[3]*alpha[5]; 
  Ghat[15] = 0.1581138830084189*alpha[22]*fUpwind[62]+0.1581138830084189*alpha[20]*fUpwind[60]+0.1581138830084189*alpha[14]*fUpwind[58]+0.1581138830084189*alpha[13]*fUpwind[57]+0.1581138830084189*alpha[12]*fUpwind[56]+0.1581138830084189*alpha[5]*fUpwind[52]+0.1581138830084189*alpha[19]*fUpwind[46]+0.1581138830084189*alpha[18]*fUpwind[45]+0.1581138830084189*alpha[11]*fUpwind[42]+0.1581138830084189*alpha[10]*fUpwind[41]+0.1581138830084189*alpha[9]*fUpwind[40]+0.1581138830084189*alpha[4]*fUpwind[36]+0.1767766952966368*alpha[8]*fUpwind[30]+0.1767766952966368*alpha[7]*fUpwind[29]+0.1767766952966368*alpha[6]*fUpwind[28]+0.1767766952966368*alpha[3]*fUpwind[25]+0.1767766952966368*alpha[2]*fUpwind[24]+0.1767766952966368*alpha[1]*fUpwind[23]+0.1767766952966368*alpha[19]*fUpwind[22]+0.1767766952966368*fUpwind[19]*alpha[22]+0.1767766952966368*alpha[18]*fUpwind[21]+0.1767766952966368*fUpwind[17]*alpha[20]+0.1767766952966368*alpha[0]*fUpwind[15]+0.1767766952966368*alpha[11]*fUpwind[14]+0.1767766952966368*fUpwind[11]*alpha[14]+0.1767766952966368*alpha[10]*fUpwind[13]+0.1767766952966368*fUpwind[10]*alpha[13]+0.1767766952966368*alpha[9]*fUpwind[12]+0.1767766952966368*fUpwind[9]*alpha[12]+0.1767766952966368*alpha[4]*fUpwind[5]+0.1767766952966368*fUpwind[4]*alpha[5]; 
  Ghat[16] = 0.1767766952966368*alpha[5]*fUpwind[27]+0.1767766952966368*alpha[4]*fUpwind[26]+0.1767766952966368*alpha[12]*fUpwind[22]+0.1767766952966368*fUpwind[12]*alpha[22]+0.1767766952966368*alpha[13]*fUpwind[21]+0.1767766952966368*alpha[14]*fUpwind[20]+0.1767766952966368*fUpwind[14]*alpha[20]+0.1767766952966368*alpha[9]*fUpwind[19]+0.1767766952966368*fUpwind[9]*alpha[19]+0.1767766952966368*alpha[10]*fUpwind[18]+0.1767766952966368*fUpwind[10]*alpha[18]+0.1767766952966368*alpha[11]*fUpwind[17]+0.1767766952966368*alpha[0]*fUpwind[16]+0.1767766952966368*alpha[1]*fUpwind[8]+0.1767766952966368*fUpwind[1]*alpha[8]+0.1767766952966368*alpha[2]*fUpwind[7]+0.1767766952966368*fUpwind[2]*alpha[7]+0.1767766952966368*alpha[3]*fUpwind[6]+0.1767766952966368*fUpwind[3]*alpha[6]; 
  Ghat[17] = 0.1581138830084189*alpha[11]*fUpwind[43]+0.1581138830084189*alpha[18]*fUpwind[39]+0.1581138830084189*alpha[19]*fUpwind[38]+0.1581138830084189*alpha[4]*fUpwind[37]+0.1581138830084189*alpha[9]*fUpwind[34]+0.1581138830084189*alpha[10]*fUpwind[33]+0.1767766952966368*alpha[14]*fUpwind[31]+0.1767766952966368*alpha[22]*fUpwind[29]+0.1767766952966368*alpha[5]*fUpwind[28]+0.1767766952966368*alpha[3]*fUpwind[26]+0.1767766952966368*alpha[12]*fUpwind[24]+0.1767766952966368*alpha[13]*fUpwind[23]+0.1767766952966368*fUpwind[15]*alpha[20]+0.1767766952966368*alpha[7]*fUpwind[19]+0.1767766952966368*fUpwind[7]*alpha[19]+0.1767766952966368*alpha[8]*fUpwind[18]+0.1767766952966368*fUpwind[8]*alpha[18]+0.1767766952966368*alpha[0]*fUpwind[17]+0.1767766952966368*alpha[11]*fUpwind[16]+0.1767766952966368*alpha[1]*fUpwind[10]+0.1767766952966368*fUpwind[1]*alpha[10]+0.1767766952966368*alpha[2]*fUpwind[9]+0.1767766952966368*fUpwind[2]*alpha[9]+0.1767766952966368*alpha[4]*fUpwind[6]+0.1767766952966368*fUpwind[4]*alpha[6]; 
  Ghat[18] = 0.1581138830084189*alpha[10]*fUpwind[43]+0.1581138830084189*alpha[4]*fUpwind[38]+0.1581138830084189*alpha[19]*fUpwind[37]+0.1581138830084189*alpha[9]*fUpwind[35]+0.1581138830084189*alpha[11]*fUpwind[33]+0.1581138830084189*alpha[18]*fUpwind[32]+0.1767766952966368*alpha[13]*fUpwind[31]+0.1767766952966368*alpha[20]*fUpwind[30]+0.1767766952966368*alpha[5]*fUpwind[29]+0.1767766952966368*alpha[22]*fUpwind[28]+0.1767766952966368*alpha[2]*fUpwind[26]+0.1767766952966368*alpha[12]*fUpwind[25]+0.1767766952966368*alpha[14]*fUpwind[23]+0.1767766952966368*alpha[6]*fUpwind[19]+0.1767766952966368*fUpwind[6]*alpha[19]+0.1767766952966368*alpha[0]*fUpwind[18]+0.1767766952966368*fUpwind[0]*alpha[18]+0.1767766952966368*alpha[8]*fUpwind[17]+0.1767766952966368*alpha[10]*fUpwind[16]+0.1767766952966368*alpha[1]*fUpwind[11]+0.1767766952966368*fUpwind[1]*alpha[11]+0.1767766952966368*alpha[3]*fUpwind[9]+0.1767766952966368*fUpwind[3]*alpha[9]+0.1767766952966368*alpha[4]*fUpwind[7]+0.1767766952966368*fUpwind[4]*alpha[7]; 
  Ghat[19] = 0.1581138830084189*alpha[9]*fUpwind[43]+0.1581138830084189*alpha[4]*fUpwind[39]+0.1581138830084189*alpha[18]*fUpwind[37]+0.1581138830084189*alpha[10]*fUpwind[35]+0.1581138830084189*alpha[11]*fUpwind[34]+0.1581138830084189*alpha[19]*fUpwind[32]+0.1767766952966368*alpha[12]*fUpwind[31]+0.1767766952966368*alpha[5]*fUpwind[30]+0.1767766952966368*alpha[20]*fUpwind[29]+0.1767766952966368*alpha[1]*fUpwind[26]+0.1767766952966368*alpha[13]*fUpwind[25]+0.1767766952966368*alpha[14]*fUpwind[24]+0.1767766952966368*fUpwind[15]*alpha[22]+0.1767766952966368*alpha[0]*fUpwind[19]+0.1767766952966368*fUpwind[0]*alpha[19]+0.1767766952966368*alpha[6]*fUpwind[18]+0.1767766952966368*fUpwind[6]*alpha[18]+0.1767766952966368*alpha[7]*fUpwind[17]+0.1767766952966368*alpha[9]*fUpwind[16]+0.1767766952966368*alpha[2]*fUpwind[11]+0.1767766952966368*fUpwind[2]*alpha[11]+0.1767766952966368*alpha[3]*fUpwind[10]+0.1767766952966368*fUpwind[3]*alpha[10]+0.1767766952966368*alpha[4]*fUpwind[8]+0.1767766952966368*fUpwind[4]*alpha[8]; 
  Ghat[20] = 0.1581138830084189*alpha[14]*fUpwind[59]+0.1581138830084189*alpha[22]*fUpwind[54]+0.1581138830084189*alpha[5]*fUpwind[53]+0.1581138830084189*alpha[12]*fUpwind[50]+0.1581138830084189*alpha[13]*fUpwind[49]+0.1581138830084189*alpha[20]*fUpwind[48]+0.1767766952966368*alpha[11]*fUpwind[31]+0.1767766952966368*alpha[18]*fUpwind[30]+0.1767766952966368*alpha[19]*fUpwind[29]+0.1767766952966368*alpha[4]*fUpwind[28]+0.1767766952966368*alpha[3]*fUpwind[27]+0.1767766952966368*alpha[9]*fUpwind[24]+0.1767766952966368*alpha[10]*fUpwind[23]+0.1767766952966368*alpha[7]*fUpwind[22]+0.1767766952966368*fUpwind[7]*alpha[22]+0.1767766952966368*alpha[8]*fUpwind[21]+0.1767766952966368*alpha[0]*fUpwind[20]+0.1767766952966368*fUpwind[0]*alpha[20]+0.1767766952966368*alpha[14]*fUpwind[16]+0.1767766952966368*alpha[1]*fUpwind[13]+0.1767766952966368*fUpwind[1]*alpha[13]+0.1767766952966368*alpha[2]*fUpwind[12]+0.1767766952966368*fUpwind[2]*alpha[12]+0.1767766952966368*alpha[5]*fUpwind[6]+0.1767766952966368*fUpwind[5]*alpha[6]; 
  Ghat[21] = 0.1581138830084189*alpha[13]*fUpwind[59]+0.1581138830084189*alpha[20]*fUpwind[55]+0.1581138830084189*alpha[5]*fUpwind[54]+0.1581138830084189*alpha[22]*fUpwind[53]+0.1581138830084189*alpha[12]*fUpwind[51]+0.1581138830084189*alpha[14]*fUpwind[49]+0.1767766952966368*alpha[10]*fUpwind[31]+0.1767766952966368*alpha[4]*fUpwind[29]+0.1767766952966368*alpha[19]*fUpwind[28]+0.1767766952966368*alpha[2]*fUpwind[27]+0.1767766952966368*alpha[9]*fUpwind[25]+0.1767766952966368*alpha[11]*fUpwind[23]+0.1767766952966368*alpha[6]*fUpwind[22]+0.1767766952966368*fUpwind[6]*alpha[22]+0.1767766952966368*alpha[0]*fUpwind[21]+0.1767766952966368*alpha[8]*fUpwind[20]+0.1767766952966368*fUpwind[8]*alpha[20]+0.1767766952966368*fUpwind[15]*alpha[18]+0.1767766952966368*alpha[13]*fUpwind[16]+0.1767766952966368*alpha[1]*fUpwind[14]+0.1767766952966368*fUpwind[1]*alpha[14]+0.1767766952966368*alpha[3]*fUpwind[12]+0.1767766952966368*fUpwind[3]*alpha[12]+0.1767766952966368*alpha[5]*fUpwind[7]+0.1767766952966368*fUpwind[5]*alpha[7]; 
  Ghat[22] = 0.1581138830084189*alpha[12]*fUpwind[59]+0.1581138830084189*alpha[5]*fUpwind[55]+0.1581138830084189*alpha[20]*fUpwind[54]+0.1581138830084189*alpha[13]*fUpwind[51]+0.1581138830084189*alpha[14]*fUpwind[50]+0.1581138830084189*alpha[22]*fUpwind[48]+0.1767766952966368*alpha[9]*fUpwind[31]+0.1767766952966368*alpha[4]*fUpwind[30]+0.1767766952966368*alpha[18]*fUpwind[28]+0.1767766952966368*alpha[1]*fUpwind[27]+0.1767766952966368*alpha[10]*fUpwind[25]+0.1767766952966368*alpha[11]*fUpwind[24]+0.1767766952966368*alpha[0]*fUpwind[22]+0.1767766952966368*fUpwind[0]*alpha[22]+0.1767766952966368*alpha[6]*fUpwind[21]+0.1767766952966368*alpha[7]*fUpwind[20]+0.1767766952966368*fUpwind[7]*alpha[20]+0.1767766952966368*fUpwind[15]*alpha[19]+0.1767766952966368*alpha[12]*fUpwind[16]+0.1767766952966368*alpha[2]*fUpwind[14]+0.1767766952966368*fUpwind[2]*alpha[14]+0.1767766952966368*alpha[3]*fUpwind[13]+0.1767766952966368*fUpwind[3]*alpha[13]+0.1767766952966368*alpha[5]*fUpwind[8]+0.1767766952966368*fUpwind[5]*alpha[8]; 
  Ghat[23] = 0.1581138830084189*alpha[22]*fUpwind[63]+0.1581138830084189*alpha[14]*fUpwind[61]+0.1581138830084189*alpha[13]*fUpwind[60]+0.1581138830084189*alpha[20]*fUpwind[57]+0.1581138830084189*alpha[5]*fUpwind[56]+0.1581138830084189*alpha[12]*fUpwind[52]+0.1581138830084189*alpha[19]*fUpwind[47]+0.1581138830084189*alpha[11]*fUpwind[45]+0.1581138830084189*alpha[10]*fUpwind[44]+0.1581138830084189*alpha[18]*fUpwind[42]+0.1581138830084189*alpha[4]*fUpwind[40]+0.1581138830084189*alpha[9]*fUpwind[36]+0.1767766952966368*alpha[8]*fUpwind[31]+0.1767766952966368*alpha[3]*fUpwind[29]+0.1767766952966368*alpha[2]*fUpwind[28]+0.1767766952966368*alpha[19]*fUpwind[27]+0.1767766952966368*alpha[22]*fUpwind[26]+0.1767766952966368*alpha[7]*fUpwind[25]+0.1767766952966368*alpha[6]*fUpwind[24]+0.1767766952966368*alpha[0]*fUpwind[23]+0.1767766952966368*alpha[11]*fUpwind[21]+0.1767766952966368*alpha[10]*fUpwind[20]+0.1767766952966368*fUpwind[10]*alpha[20]+0.1767766952966368*alpha[14]*fUpwind[18]+0.1767766952966368*fUpwind[14]*alpha[18]+0.1767766952966368*alpha[13]*fUpwind[17]+0.1767766952966368*alpha[1]*fUpwind[15]+0.1767766952966368*alpha[4]*fUpwind[12]+0.1767766952966368*fUpwind[4]*alpha[12]+0.1767766952966368*alpha[5]*fUpwind[9]+0.1767766952966368*fUpwind[5]*alpha[9]; 
  Ghat[24] = 0.1581138830084189*alpha[14]*fUpwind[62]+0.1581138830084189*alpha[12]*fUpwind[60]+0.1581138830084189*alpha[22]*fUpwind[58]+0.1581138830084189*alpha[5]*fUpwind[57]+0.1581138830084189*alpha[20]*fUpwind[56]+0.1581138830084189*alpha[13]*fUpwind[52]+0.1581138830084189*alpha[18]*fUpwind[47]+0.1581138830084189*alpha[11]*fUpwind[46]+0.1581138830084189*alpha[9]*fUpwind[44]+0.1581138830084189*alpha[19]*fUpwind[42]+0.1581138830084189*alpha[4]*fUpwind[41]+0.1581138830084189*alpha[10]*fUpwind[36]+0.1767766952966368*alpha[7]*fUpwind[31]+0.1767766952966368*alpha[3]*fUpwind[30]+0.1767766952966368*alpha[1]*fUpwind[28]+0.1767766952966368*alpha[18]*fUpwind[27]+0.1767766952966368*alpha[8]*fUpwind[25]+0.1767766952966368*alpha[0]*fUpwind[24]+0.1767766952966368*alpha[6]*fUpwind[23]+0.1767766952966368*alpha[11]*fUpwind[22]+0.1767766952966368*fUpwind[11]*alpha[22]+0.1767766952966368*alpha[9]*fUpwind[20]+0.1767766952966368*fUpwind[9]*alpha[20]+0.1767766952966368*alpha[14]*fUpwind[19]+0.1767766952966368*fUpwind[14]*alpha[19]+0.1767766952966368*alpha[12]*fUpwind[17]+0.1767766952966368*alpha[2]*fUpwind[15]+0.1767766952966368*alpha[4]*fUpwind[13]+0.1767766952966368*fUpwind[4]*alpha[13]+0.1767766952966368*alpha[5]*fUpwind[10]+0.1767766952966368*fUpwind[5]*alpha[10]; 
  Ghat[25] = 0.1581138830084189*alpha[20]*fUpwind[63]+0.1581138830084189*alpha[13]*fUpwind[62]+0.1581138830084189*alpha[12]*fUpwind[61]+0.1581138830084189*alpha[5]*fUpwind[58]+0.1581138830084189*alpha[22]*fUpwind[57]+0.1581138830084189*alpha[14]*fUpwind[52]+0.1581138830084189*alpha[10]*fUpwind[46]+0.1581138830084189*alpha[9]*fUpwind[45]+0.1581138830084189*alpha[4]*fUpwind[42]+0.1581138830084189*alpha[19]*fUpwind[41]+0.1581138830084189*alpha[18]*fUpwind[40]+0.1581138830084189*alpha[11]*fUpwind[36]+0.1767766952966368*alpha[6]*fUpwind[31]+0.1767766952966368*alpha[2]*fUpwind[30]+0.1767766952966368*alpha[1]*fUpwind[29]+0.1767766952966368*alpha[20]*fUpwind[26]+0.1767766952966368*alpha[0]*fUpwind[25]+0.1767766952966368*alpha[8]*fUpwind[24]+0.1767766952966368*alpha[7]*fUpwind[23]+0.1767766952966368*alpha[10]*fUpwind[22]+0.1767766952966368*fUpwind[10]*alpha[22]+0.1767766952966368*alpha[9]*fUpwind[21]+0.1767766952966368*alpha[13]*fUpwind[19]+0.1767766952966368*fUpwind[13]*alpha[19]+0.1767766952966368*alpha[12]*fUpwind[18]+0.1767766952966368*fUpwind[12]*alpha[18]+0.1767766952966368*alpha[3]*fUpwind[15]+0.1767766952966368*alpha[4]*fUpwind[14]+0.1767766952966368*fUpwind[4]*alpha[14]+0.1767766952966368*alpha[5]*fUpwind[11]+0.1767766952966368*fUpwind[5]*alpha[11]; 
  Ghat[26] = 0.1581138830084189*alpha[4]*fUpwind[43]+0.1581138830084189*alpha[9]*fUpwind[39]+0.1581138830084189*alpha[10]*fUpwind[38]+0.1581138830084189*alpha[11]*fUpwind[37]+0.1581138830084189*alpha[18]*fUpwind[34]+0.1581138830084189*alpha[19]*fUpwind[33]+0.1767766952966368*alpha[5]*fUpwind[31]+0.1767766952966368*alpha[12]*fUpwind[30]+0.1767766952966368*alpha[13]*fUpwind[29]+0.1767766952966368*alpha[14]*fUpwind[28]+0.1767766952966368*alpha[0]*fUpwind[26]+0.1767766952966368*alpha[20]*fUpwind[25]+0.1767766952966368*alpha[22]*fUpwind[23]+0.1767766952966368*alpha[1]*fUpwind[19]+0.1767766952966368*fUpwind[1]*alpha[19]+0.1767766952966368*alpha[2]*fUpwind[18]+0.1767766952966368*fUpwind[2]*alpha[18]+0.1767766952966368*alpha[3]*fUpwind[17]+0.1767766952966368*alpha[4]*fUpwind[16]+0.1767766952966368*alpha[6]*fUpwind[11]+0.1767766952966368*fUpwind[6]*alpha[11]+0.1767766952966368*alpha[7]*fUpwind[10]+0.1767766952966368*fUpwind[7]*alpha[10]+0.1767766952966368*alpha[8]*fUpwind[9]+0.1767766952966368*fUpwind[8]*alpha[9]; 
  Ghat[27] = 0.1581138830084189*alpha[5]*fUpwind[59]+0.1581138830084189*alpha[12]*fUpwind[55]+0.1581138830084189*alpha[13]*fUpwind[54]+0.1581138830084189*alpha[14]*fUpwind[53]+0.1581138830084189*alpha[20]*fUpwind[51]+0.1581138830084189*alpha[22]*fUpwind[49]+0.1767766952966368*alpha[4]*fUpwind[31]+0.1767766952966368*alpha[9]*fUpwind[30]+0.1767766952966368*alpha[10]*fUpwind[29]+0.1767766952966368*alpha[11]*fUpwind[28]+0.1767766952966368*alpha[0]*fUpwind[27]+0.1767766952966368*alpha[18]*fUpwind[24]+0.1767766952966368*alpha[19]*fUpwind[23]+0.1767766952966368*alpha[1]*fUpwind[22]+0.1767766952966368*fUpwind[1]*alpha[22]+0.1767766952966368*alpha[2]*fUpwind[21]+0.1767766952966368*alpha[3]*fUpwind[20]+0.1767766952966368*fUpwind[3]*alpha[20]+0.1767766952966368*alpha[5]*fUpwind[16]+0.1767766952966368*alpha[6]*fUpwind[14]+0.1767766952966368*fUpwind[6]*alpha[14]+0.1767766952966368*alpha[7]*fUpwind[13]+0.1767766952966368*fUpwind[7]*alpha[13]+0.1767766952966368*alpha[8]*fUpwind[12]+0.1767766952966368*fUpwind[8]*alpha[12]; 
  Ghat[28] = 0.1581138830084189*alpha[14]*fUpwind[63]+0.1581138830084189*alpha[22]*fUpwind[61]+0.1581138830084189*alpha[5]*fUpwind[60]+0.1581138830084189*alpha[12]*fUpwind[57]+0.1581138830084189*alpha[13]*fUpwind[56]+0.1581138830084189*alpha[20]*fUpwind[52]+0.1581138830084189*alpha[11]*fUpwind[47]+0.1581138830084189*alpha[18]*fUpwind[46]+0.1581138830084189*alpha[19]*fUpwind[45]+0.1581138830084189*alpha[4]*fUpwind[44]+0.1581138830084189*alpha[9]*fUpwind[41]+0.1581138830084189*alpha[10]*fUpwind[40]+0.1767766952966368*alpha[3]*fUpwind[31]+0.1767766952966368*alpha[7]*fUpwind[30]+0.1767766952966368*alpha[8]*fUpwind[29]+0.1767766952966368*alpha[0]*fUpwind[28]+0.1767766952966368*alpha[11]*fUpwind[27]+0.1767766952966368*alpha[14]*fUpwind[26]+0.1767766952966368*alpha[1]*fUpwind[24]+0.1767766952966368*alpha[2]*fUpwind[23]+0.1767766952966368*alpha[18]*fUpwind[22]+0.1767766952966368*fUpwind[18]*alpha[22]+0.1767766952966368*alpha[19]*fUpwind[21]+0.1767766952966368*alpha[4]*fUpwind[20]+0.1767766952966368*fUpwind[4]*alpha[20]+0.1767766952966368*alpha[5]*fUpwind[17]+0.1767766952966368*alpha[6]*fUpwind[15]+0.1767766952966368*alpha[9]*fUpwind[13]+0.1767766952966368*fUpwind[9]*alpha[13]+0.1767766952966368*alpha[10]*fUpwind[12]+0.1767766952966368*fUpwind[10]*alpha[12]; 
  Ghat[29] = 0.1581138830084189*alpha[13]*fUpwind[63]+0.1581138830084189*alpha[20]*fUpwind[62]+0.1581138830084189*alpha[5]*fUpwind[61]+0.1581138830084189*alpha[22]*fUpwind[60]+0.1581138830084189*alpha[12]*fUpwind[58]+0.1581138830084189*alpha[14]*fUpwind[56]+0.1581138830084189*alpha[10]*fUpwind[47]+0.1581138830084189*alpha[4]*fUpwind[45]+0.1581138830084189*alpha[19]*fUpwind[44]+0.1581138830084189*alpha[9]*fUpwind[42]+0.1581138830084189*alpha[11]*fUpwind[40]+0.1581138830084189*alpha[18]*fUpwind[36]+0.1767766952966368*alpha[2]*fUpwind[31]+0.1767766952966368*alpha[6]*fUpwind[30]+0.1767766952966368*alpha[0]*fUpwind[29]+0.1767766952966368*alpha[8]*fUpwind[28]+0.1767766952966368*alpha[10]*fUpwind[27]+0.1767766952966368*alpha[13]*fUpwind[26]+0.1767766952966368*alpha[1]*fUpwind[25]+0.1767766952966368*alpha[3]*fUpwind[23]+0.1767766952966368*fUpwind[17]*alpha[22]+0.1767766952966368*alpha[4]*fUpwind[21]+0.1767766952966368*alpha[19]*fUpwind[20]+0.1767766952966368*fUpwind[19]*alpha[20]+0.1767766952966368*alpha[5]*fUpwind[18]+0.1767766952966368*fUpwind[5]*alpha[18]+0.1767766952966368*alpha[7]*fUpwind[15]+0.1767766952966368*alpha[9]*fUpwind[14]+0.1767766952966368*fUpwind[9]*alpha[14]+0.1767766952966368*alpha[11]*fUpwind[12]+0.1767766952966368*fUpwind[11]*alpha[12]; 
  Ghat[30] = 0.1581138830084189*alpha[12]*fUpwind[63]+0.1581138830084189*alpha[5]*fUpwind[62]+0.1581138830084189*alpha[20]*fUpwind[61]+0.1581138830084189*alpha[13]*fUpwind[58]+0.1581138830084189*alpha[14]*fUpwind[57]+0.1581138830084189*alpha[22]*fUpwind[52]+0.1581138830084189*alpha[9]*fUpwind[47]+0.1581138830084189*alpha[4]*fUpwind[46]+0.1581138830084189*alpha[18]*fUpwind[44]+0.1581138830084189*alpha[10]*fUpwind[42]+0.1581138830084189*alpha[11]*fUpwind[41]+0.1581138830084189*alpha[19]*fUpwind[36]+0.1767766952966368*alpha[1]*fUpwind[31]+0.1767766952966368*alpha[0]*fUpwind[30]+0.1767766952966368*alpha[6]*fUpwind[29]+0.1767766952966368*alpha[7]*fUpwind[28]+0.1767766952966368*alpha[9]*fUpwind[27]+0.1767766952966368*alpha[12]*fUpwind[26]+0.1767766952966368*alpha[2]*fUpwind[25]+0.1767766952966368*alpha[3]*fUpwind[24]+0.1767766952966368*alpha[4]*fUpwind[22]+0.1767766952966368*fUpwind[4]*alpha[22]+0.1767766952966368*alpha[18]*fUpwind[20]+0.1767766952966368*fUpwind[18]*alpha[20]+0.1767766952966368*alpha[5]*fUpwind[19]+0.1767766952966368*fUpwind[5]*alpha[19]+0.1767766952966368*alpha[8]*fUpwind[15]+0.1767766952966368*alpha[10]*fUpwind[14]+0.1767766952966368*fUpwind[10]*alpha[14]+0.1767766952966368*alpha[11]*fUpwind[13]+0.1767766952966368*fUpwind[11]*alpha[13]; 
  Ghat[31] = 0.1581138830084189*alpha[5]*fUpwind[63]+0.1581138830084189*alpha[12]*fUpwind[62]+0.1581138830084189*alpha[13]*fUpwind[61]+0.1581138830084189*alpha[14]*fUpwind[60]+0.1581138830084189*alpha[20]*fUpwind[58]+0.1581138830084189*alpha[22]*fUpwind[56]+0.1581138830084189*alpha[4]*fUpwind[47]+0.1581138830084189*alpha[9]*fUpwind[46]+0.1581138830084189*alpha[10]*fUpwind[45]+0.1581138830084189*alpha[11]*fUpwind[44]+0.1581138830084189*alpha[18]*fUpwind[41]+0.1581138830084189*alpha[19]*fUpwind[40]+0.1767766952966368*alpha[0]*fUpwind[31]+0.1767766952966368*alpha[1]*fUpwind[30]+0.1767766952966368*alpha[2]*fUpwind[29]+0.1767766952966368*alpha[3]*fUpwind[28]+0.1767766952966368*alpha[4]*fUpwind[27]+0.1767766952966368*alpha[5]*fUpwind[26]+0.1767766952966368*alpha[6]*fUpwind[25]+0.1767766952966368*alpha[7]*fUpwind[24]+0.1767766952966368*alpha[8]*fUpwind[23]+0.1767766952966368*alpha[9]*fUpwind[22]+0.1767766952966368*fUpwind[9]*alpha[22]+0.1767766952966368*alpha[10]*fUpwind[21]+0.1767766952966368*alpha[11]*fUpwind[20]+0.1767766952966368*fUpwind[11]*alpha[20]+0.1767766952966368*alpha[12]*fUpwind[19]+0.1767766952966368*fUpwind[12]*alpha[19]+0.1767766952966368*alpha[13]*fUpwind[18]+0.1767766952966368*fUpwind[13]*alpha[18]+0.1767766952966368*alpha[14]*fUpwind[17]; 
  Ghat[32] = 0.1767766952966368*alpha[22]*fUpwind[46]+0.1767766952966368*alpha[20]*fUpwind[44]+0.1767766952966368*alpha[14]*fUpwind[42]+0.1767766952966368*alpha[13]*fUpwind[41]+0.1767766952966368*alpha[12]*fUpwind[40]+0.1767766952966368*alpha[8]*fUpwind[39]+0.1767766952966368*alpha[7]*fUpwind[38]+0.1767766952966368*alpha[6]*fUpwind[37]+0.1767766952966368*alpha[5]*fUpwind[36]+0.1767766952966368*alpha[3]*fUpwind[35]+0.1767766952966368*alpha[2]*fUpwind[34]+0.1767766952966368*alpha[1]*fUpwind[33]+0.1767766952966368*alpha[0]*fUpwind[32]+0.1581138830084189*alpha[19]*fUpwind[19]+0.1581138830084189*alpha[18]*fUpwind[18]+0.1581138830084189*alpha[11]*fUpwind[11]+0.1581138830084189*alpha[10]*fUpwind[10]+0.1581138830084189*alpha[9]*fUpwind[9]+0.1581138830084189*alpha[4]*fUpwind[4]; 
  Ghat[33] = 0.1767766952966368*alpha[22]*fUpwind[47]+0.1767766952966368*alpha[14]*fUpwind[45]+0.1767766952966368*alpha[13]*fUpwind[44]+0.1767766952966368*alpha[8]*fUpwind[43]+0.1767766952966368*alpha[20]*fUpwind[41]+0.1767766952966368*alpha[5]*fUpwind[40]+0.1767766952966368*alpha[3]*fUpwind[38]+0.1767766952966368*alpha[2]*fUpwind[37]+0.1767766952966368*alpha[12]*fUpwind[36]+0.1767766952966368*alpha[7]*fUpwind[35]+0.1767766952966368*alpha[6]*fUpwind[34]+0.1767766952966368*alpha[0]*fUpwind[33]+0.1767766952966368*alpha[1]*fUpwind[32]+0.1581138830084189*alpha[19]*fUpwind[26]+0.1581138830084189*alpha[11]*fUpwind[18]+0.1581138830084189*fUpwind[11]*alpha[18]+0.1581138830084189*alpha[10]*fUpwind[17]+0.1581138830084189*alpha[4]*fUpwind[9]+0.1581138830084189*fUpwind[4]*alpha[9]; 
  Ghat[34] = 0.1767766952966368*alpha[14]*fUpwind[46]+0.1767766952966368*alpha[12]*fUpwind[44]+0.1767766952966368*alpha[7]*fUpwind[43]+0.1767766952966368*alpha[22]*fUpwind[42]+0.1767766952966368*alpha[5]*fUpwind[41]+0.1767766952966368*alpha[20]*fUpwind[40]+0.1767766952966368*alpha[3]*fUpwind[39]+0.1767766952966368*alpha[1]*fUpwind[37]+0.1767766952966368*alpha[13]*fUpwind[36]+0.1767766952966368*alpha[8]*fUpwind[35]+0.1767766952966368*alpha[0]*fUpwind[34]+0.1767766952966368*alpha[6]*fUpwind[33]+0.1767766952966368*alpha[2]*fUpwind[32]+0.1581138830084189*alpha[18]*fUpwind[26]+0.1581138830084189*alpha[11]*fUpwind[19]+0.1581138830084189*fUpwind[11]*alpha[19]+0.1581138830084189*alpha[9]*fUpwind[17]+0.1581138830084189*alpha[4]*fUpwind[10]+0.1581138830084189*fUpwind[4]*alpha[10]; 
  Ghat[35] = 0.1767766952966368*alpha[20]*fUpwind[47]+0.1767766952966368*alpha[13]*fUpwind[46]+0.1767766952966368*alpha[12]*fUpwind[45]+0.1767766952966368*alpha[6]*fUpwind[43]+0.1767766952966368*alpha[5]*fUpwind[42]+0.1767766952966368*alpha[22]*fUpwind[41]+0.1767766952966368*alpha[2]*fUpwind[39]+0.1767766952966368*alpha[1]*fUpwind[38]+0.1767766952966368*alpha[14]*fUpwind[36]+0.1767766952966368*alpha[0]*fUpwind[35]+0.1767766952966368*alpha[8]*fUpwind[34]+0.1767766952966368*alpha[7]*fUpwind[33]+0.1767766952966368*alpha[3]*fUpwind[32]+0.1581138830084189*alpha[10]*fUpwind[19]+0.1581138830084189*fUpwind[10]*alpha[19]+0.1581138830084189*alpha[9]*fUpwind[18]+0.1581138830084189*fUpwind[9]*alpha[18]+0.1581138830084189*alpha[4]*fUpwind[11]+0.1581138830084189*fUpwind[4]*alpha[11]; 
  Ghat[36] = 0.1767766952966368*alpha[8]*fUpwind[46]+0.1767766952966368*alpha[7]*fUpwind[45]+0.1767766952966368*alpha[6]*fUpwind[44]+0.1767766952966368*alpha[3]*fUpwind[42]+0.1767766952966368*alpha[2]*fUpwind[41]+0.1767766952966368*alpha[1]*fUpwind[40]+0.1767766952966368*alpha[22]*fUpwind[39]+0.1767766952966368*alpha[20]*fUpwind[37]+0.1767766952966368*alpha[0]*fUpwind[36]+0.1767766952966368*alpha[14]*fUpwind[35]+0.1767766952966368*alpha[13]*fUpwind[34]+0.1767766952966368*alpha[12]*fUpwind[33]+0.1767766952966368*alpha[5]*fUpwind[32]+0.1581138830084189*alpha[19]*fUpwind[30]+0.1581138830084189*alpha[18]*fUpwind[29]+0.1581138830084189*alpha[11]*fUpwind[25]+0.1581138830084189*alpha[10]*fUpwind[24]+0.1581138830084189*alpha[9]*fUpwind[23]+0.1581138830084189*alpha[4]*fUpwind[15]; 
  Ghat[37] = 0.1767766952966368*alpha[14]*fUpwind[47]+0.1767766952966368*alpha[22]*fUpwind[45]+0.1767766952966368*alpha[5]*fUpwind[44]+0.1767766952966368*alpha[3]*fUpwind[43]+0.1767766952966368*alpha[12]*fUpwind[41]+0.1767766952966368*alpha[13]*fUpwind[40]+0.1767766952966368*alpha[7]*fUpwind[39]+0.1767766952966368*alpha[8]*fUpwind[38]+0.1767766952966368*alpha[0]*fUpwind[37]+0.1767766952966368*alpha[20]*fUpwind[36]+0.1767766952966368*alpha[1]*fUpwind[34]+0.1767766952966368*alpha[2]*fUpwind[33]+0.1767766952966368*alpha[6]*fUpwind[32]+0.1581138830084189*alpha[11]*fUpwind[26]+0.1581138830084189*alpha[18]*fUpwind[19]+0.1581138830084189*fUpwind[18]*alpha[19]+0.1581138830084189*alpha[4]*fUpwind[17]+0.1581138830084189*alpha[9]*fUpwind[10]+0.1581138830084189*fUpwind[9]*alpha[10]; 
  Ghat[38] = 0.1767766952966368*alpha[13]*fUpwind[47]+0.1767766952966368*alpha[20]*fUpwind[46]+0.1767766952966368*alpha[5]*fUpwind[45]+0.1767766952966368*alpha[22]*fUpwind[44]+0.1767766952966368*alpha[2]*fUpwind[43]+0.1767766952966368*alpha[12]*fUpwind[42]+0.1767766952966368*alpha[14]*fUpwind[40]+0.1767766952966368*alpha[6]*fUpwind[39]+0.1767766952966368*alpha[0]*fUpwind[38]+0.1767766952966368*alpha[8]*fUpwind[37]+0.1767766952966368*alpha[1]*fUpwind[35]+0.1767766952966368*alpha[3]*fUpwind[33]+0.1767766952966368*alpha[7]*fUpwind[32]+0.1581138830084189*alpha[10]*fUpwind[26]+0.1581138830084189*fUpwind[17]*alpha[19]+0.1581138830084189*alpha[4]*fUpwind[18]+0.1581138830084189*fUpwind[4]*alpha[18]+0.1581138830084189*alpha[9]*fUpwind[11]+0.1581138830084189*fUpwind[9]*alpha[11]; 
  Ghat[39] = 0.1767766952966368*alpha[12]*fUpwind[47]+0.1767766952966368*alpha[5]*fUpwind[46]+0.1767766952966368*alpha[20]*fUpwind[45]+0.1767766952966368*alpha[1]*fUpwind[43]+0.1767766952966368*alpha[13]*fUpwind[42]+0.1767766952966368*alpha[14]*fUpwind[41]+0.1767766952966368*alpha[0]*fUpwind[39]+0.1767766952966368*alpha[6]*fUpwind[38]+0.1767766952966368*alpha[7]*fUpwind[37]+0.1767766952966368*alpha[22]*fUpwind[36]+0.1767766952966368*alpha[2]*fUpwind[35]+0.1767766952966368*alpha[3]*fUpwind[34]+0.1767766952966368*alpha[8]*fUpwind[32]+0.1581138830084189*alpha[9]*fUpwind[26]+0.1581138830084189*alpha[4]*fUpwind[19]+0.1581138830084189*fUpwind[4]*alpha[19]+0.1581138830084189*fUpwind[17]*alpha[18]+0.1581138830084189*alpha[10]*fUpwind[11]+0.1581138830084189*fUpwind[10]*alpha[11]; 
  Ghat[40] = 0.1767766952966368*alpha[8]*fUpwind[47]+0.1767766952966368*alpha[3]*fUpwind[45]+0.1767766952966368*alpha[2]*fUpwind[44]+0.1767766952966368*alpha[22]*fUpwind[43]+0.1767766952966368*alpha[7]*fUpwind[42]+0.1767766952966368*alpha[6]*fUpwind[41]+0.1767766952966368*alpha[0]*fUpwind[40]+0.1767766952966368*alpha[14]*fUpwind[38]+0.1767766952966368*alpha[13]*fUpwind[37]+0.1767766952966368*alpha[1]*fUpwind[36]+0.1767766952966368*alpha[20]*fUpwind[34]+0.1767766952966368*alpha[5]*fUpwind[33]+0.1767766952966368*alpha[12]*fUpwind[32]+0.1581138830084189*alpha[19]*fUpwind[31]+0.1581138830084189*alpha[11]*fUpwind[29]+0.1581138830084189*alpha[10]*fUpwind[28]+0.1581138830084189*alpha[18]*fUpwind[25]+0.1581138830084189*alpha[4]*fUpwind[23]+0.1581138830084189*alpha[9]*fUpwind[15]; 
  Ghat[41] = 0.1767766952966368*alpha[7]*fUpwind[47]+0.1767766952966368*alpha[3]*fUpwind[46]+0.1767766952966368*alpha[1]*fUpwind[44]+0.1767766952966368*alpha[8]*fUpwind[42]+0.1767766952966368*alpha[0]*fUpwind[41]+0.1767766952966368*alpha[6]*fUpwind[40]+0.1767766952966368*alpha[14]*fUpwind[39]+0.1767766952966368*alpha[12]*fUpwind[37]+0.1767766952966368*alpha[2]*fUpwind[36]+0.1767766952966368*alpha[22]*fUpwind[35]+0.1767766952966368*alpha[5]*fUpwind[34]+0.1767766952966368*alpha[20]*fUpwind[33]+0.1767766952966368*alpha[13]*fUpwind[32]+0.1581138830084189*alpha[18]*fUpwind[31]+0.1581138830084189*alpha[11]*fUpwind[30]+0.1581138830084189*alpha[9]*fUpwind[28]+0.1581138830084189*alpha[19]*fUpwind[25]+0.1581138830084189*alpha[4]*fUpwind[24]+0.1581138830084189*alpha[10]*fUpwind[15]; 
  Ghat[42] = 0.1767766952966368*alpha[6]*fUpwind[47]+0.1767766952966368*alpha[2]*fUpwind[46]+0.1767766952966368*alpha[1]*fUpwind[45]+0.1767766952966368*alpha[20]*fUpwind[43]+0.1767766952966368*alpha[0]*fUpwind[42]+0.1767766952966368*alpha[8]*fUpwind[41]+0.1767766952966368*alpha[7]*fUpwind[40]+0.1767766952966368*alpha[13]*fUpwind[39]+0.1767766952966368*alpha[12]*fUpwind[38]+0.1767766952966368*alpha[3]*fUpwind[36]+0.1767766952966368*alpha[5]*fUpwind[35]+0.1767766952966368*alpha[22]*fUpwind[34]+0.1767766952966368*alpha[14]*fUpwind[32]+0.1581138830084189*alpha[10]*fUpwind[30]+0.1581138830084189*alpha[9]*fUpwind[29]+0.1581138830084189*alpha[4]*fUpwind[25]+0.1581138830084189*alpha[19]*fUpwind[24]+0.1581138830084189*alpha[18]*fUpwind[23]+0.1581138830084189*alpha[11]*fUpwind[15]; 
  Ghat[43] = 0.1767766952966368*alpha[5]*fUpwind[47]+0.1767766952966368*alpha[12]*fUpwind[46]+0.1767766952966368*alpha[13]*fUpwind[45]+0.1767766952966368*alpha[14]*fUpwind[44]+0.1767766952966368*alpha[0]*fUpwind[43]+0.1767766952966368*alpha[20]*fUpwind[42]+0.1767766952966368*alpha[22]*fUpwind[40]+0.1767766952966368*alpha[1]*fUpwind[39]+0.1767766952966368*alpha[2]*fUpwind[38]+0.1767766952966368*alpha[3]*fUpwind[37]+0.1767766952966368*alpha[6]*fUpwind[35]+0.1767766952966368*alpha[7]*fUpwind[34]+0.1767766952966368*alpha[8]*fUpwind[33]+0.1581138830084189*alpha[4]*fUpwind[26]+0.1581138830084189*alpha[9]*fUpwind[19]+0.1581138830084189*fUpwind[9]*alpha[19]+0.1581138830084189*alpha[10]*fUpwind[18]+0.1581138830084189*fUpwind[10]*alpha[18]+0.1581138830084189*alpha[11]*fUpwind[17]; 
  Ghat[44] = 0.1767766952966368*alpha[3]*fUpwind[47]+0.1767766952966368*alpha[7]*fUpwind[46]+0.1767766952966368*alpha[8]*fUpwind[45]+0.1767766952966368*alpha[0]*fUpwind[44]+0.1767766952966368*alpha[14]*fUpwind[43]+0.1767766952966368*alpha[1]*fUpwind[41]+0.1767766952966368*alpha[2]*fUpwind[40]+0.1767766952966368*alpha[22]*fUpwind[38]+0.1767766952966368*alpha[5]*fUpwind[37]+0.1767766952966368*alpha[6]*fUpwind[36]+0.1767766952966368*alpha[12]*fUpwind[34]+0.1767766952966368*alpha[13]*fUpwind[33]+0.1767766952966368*alpha[20]*fUpwind[32]+0.1581138830084189*alpha[11]*fUpwind[31]+0.1581138830084189*alpha[18]*fUpwind[30]+0.1581138830084189*alpha[19]*fUpwind[29]+0.1581138830084189*alpha[4]*fUpwind[28]+0.1581138830084189*alpha[9]*fUpwind[24]+0.1581138830084189*alpha[10]*fUpwind[23]; 
  Ghat[45] = 0.1767766952966368*alpha[2]*fUpwind[47]+0.1767766952966368*alpha[6]*fUpwind[46]+0.1767766952966368*alpha[0]*fUpwind[45]+0.1767766952966368*alpha[8]*fUpwind[44]+0.1767766952966368*alpha[13]*fUpwind[43]+0.1767766952966368*alpha[1]*fUpwind[42]+0.1767766952966368*alpha[3]*fUpwind[40]+0.1767766952966368*alpha[20]*fUpwind[39]+0.1767766952966368*alpha[5]*fUpwind[38]+0.1767766952966368*alpha[22]*fUpwind[37]+0.1767766952966368*alpha[7]*fUpwind[36]+0.1767766952966368*alpha[12]*fUpwind[35]+0.1767766952966368*alpha[14]*fUpwind[33]+0.1581138830084189*alpha[10]*fUpwind[31]+0.1581138830084189*alpha[4]*fUpwind[29]+0.1581138830084189*alpha[19]*fUpwind[28]+0.1581138830084189*alpha[9]*fUpwind[25]+0.1581138830084189*alpha[11]*fUpwind[23]+0.1581138830084189*fUpwind[15]*alpha[18]; 
  Ghat[46] = 0.1767766952966368*alpha[1]*fUpwind[47]+0.1767766952966368*alpha[0]*fUpwind[46]+0.1767766952966368*alpha[6]*fUpwind[45]+0.1767766952966368*alpha[7]*fUpwind[44]+0.1767766952966368*alpha[12]*fUpwind[43]+0.1767766952966368*alpha[2]*fUpwind[42]+0.1767766952966368*alpha[3]*fUpwind[41]+0.1767766952966368*alpha[5]*fUpwind[39]+0.1767766952966368*alpha[20]*fUpwind[38]+0.1767766952966368*alpha[8]*fUpwind[36]+0.1767766952966368*alpha[13]*fUpwind[35]+0.1767766952966368*alpha[14]*fUpwind[34]+0.1767766952966368*alpha[22]*fUpwind[32]+0.1581138830084189*alpha[9]*fUpwind[31]+0.1581138830084189*alpha[4]*fUpwind[30]+0.1581138830084189*alpha[18]*fUpwind[28]+0.1581138830084189*alpha[10]*fUpwind[25]+0.1581138830084189*alpha[11]*fUpwind[24]+0.1581138830084189*fUpwind[15]*alpha[19]; 
  Ghat[47] = 0.1767766952966368*alpha[0]*fUpwind[47]+0.1767766952966368*alpha[1]*fUpwind[46]+0.1767766952966368*alpha[2]*fUpwind[45]+0.1767766952966368*alpha[3]*fUpwind[44]+0.1767766952966368*alpha[5]*fUpwind[43]+0.1767766952966368*alpha[6]*fUpwind[42]+0.1767766952966368*alpha[7]*fUpwind[41]+0.1767766952966368*alpha[8]*fUpwind[40]+0.1767766952966368*alpha[12]*fUpwind[39]+0.1767766952966368*alpha[13]*fUpwind[38]+0.1767766952966368*alpha[14]*fUpwind[37]+0.1767766952966368*alpha[20]*fUpwind[35]+0.1767766952966368*alpha[22]*fUpwind[33]+0.1581138830084189*alpha[4]*fUpwind[31]+0.1581138830084189*alpha[9]*fUpwind[30]+0.1581138830084189*alpha[10]*fUpwind[29]+0.1581138830084189*alpha[11]*fUpwind[28]+0.1581138830084189*alpha[18]*fUpwind[24]+0.1581138830084189*alpha[19]*fUpwind[23]; 
  Ghat[48] = 0.1767766952966368*alpha[19]*fUpwind[62]+0.1767766952966368*alpha[18]*fUpwind[61]+0.1767766952966368*alpha[11]*fUpwind[58]+0.1767766952966368*alpha[10]*fUpwind[57]+0.1767766952966368*alpha[9]*fUpwind[56]+0.1767766952966368*alpha[8]*fUpwind[55]+0.1767766952966368*alpha[7]*fUpwind[54]+0.1767766952966368*alpha[6]*fUpwind[53]+0.1767766952966368*alpha[4]*fUpwind[52]+0.1767766952966368*alpha[3]*fUpwind[51]+0.1767766952966368*alpha[2]*fUpwind[50]+0.1767766952966368*alpha[1]*fUpwind[49]+0.1767766952966368*alpha[0]*fUpwind[48]+0.1581138830084189*alpha[22]*fUpwind[22]+0.1581138830084189*alpha[20]*fUpwind[20]+0.1581138830084189*alpha[14]*fUpwind[14]+0.1581138830084189*alpha[13]*fUpwind[13]+0.1581138830084189*alpha[12]*fUpwind[12]+0.1581138830084189*alpha[5]*fUpwind[5]; 
  Ghat[49] = 0.1767766952966368*alpha[19]*fUpwind[63]+0.1767766952966368*alpha[11]*fUpwind[61]+0.1767766952966368*alpha[10]*fUpwind[60]+0.1767766952966368*alpha[8]*fUpwind[59]+0.1767766952966368*alpha[18]*fUpwind[58]+0.1767766952966368*alpha[4]*fUpwind[56]+0.1767766952966368*alpha[3]*fUpwind[54]+0.1767766952966368*alpha[2]*fUpwind[53]+0.1767766952966368*alpha[9]*fUpwind[52]+0.1767766952966368*alpha[7]*fUpwind[51]+0.1767766952966368*alpha[6]*fUpwind[50]+0.1767766952966368*alpha[0]*fUpwind[49]+0.1767766952966368*alpha[1]*fUpwind[48]+0.1581138830084189*alpha[22]*fUpwind[27]+0.1581138830084189*alpha[14]*fUpwind[21]+0.1581138830084189*alpha[13]*fUpwind[20]+0.1581138830084189*fUpwind[13]*alpha[20]+0.1581138830084189*alpha[5]*fUpwind[12]+0.1581138830084189*fUpwind[5]*alpha[12]; 
  Ghat[50] = 0.1767766952966368*alpha[18]*fUpwind[63]+0.1767766952966368*alpha[11]*fUpwind[62]+0.1767766952966368*alpha[9]*fUpwind[60]+0.1767766952966368*alpha[7]*fUpwind[59]+0.1767766952966368*alpha[19]*fUpwind[58]+0.1767766952966368*alpha[4]*fUpwind[57]+0.1767766952966368*alpha[3]*fUpwind[55]+0.1767766952966368*alpha[1]*fUpwind[53]+0.1767766952966368*alpha[10]*fUpwind[52]+0.1767766952966368*alpha[8]*fUpwind[51]+0.1767766952966368*alpha[0]*fUpwind[50]+0.1767766952966368*alpha[6]*fUpwind[49]+0.1767766952966368*alpha[2]*fUpwind[48]+0.1581138830084189*alpha[14]*fUpwind[22]+0.1581138830084189*fUpwind[14]*alpha[22]+0.1581138830084189*alpha[12]*fUpwind[20]+0.1581138830084189*fUpwind[12]*alpha[20]+0.1581138830084189*alpha[5]*fUpwind[13]+0.1581138830084189*fUpwind[5]*alpha[13]; 
  Ghat[51] = 0.1767766952966368*alpha[10]*fUpwind[62]+0.1767766952966368*alpha[9]*fUpwind[61]+0.1767766952966368*alpha[6]*fUpwind[59]+0.1767766952966368*alpha[4]*fUpwind[58]+0.1767766952966368*alpha[19]*fUpwind[57]+0.1767766952966368*alpha[18]*fUpwind[56]+0.1767766952966368*alpha[2]*fUpwind[55]+0.1767766952966368*alpha[1]*fUpwind[54]+0.1767766952966368*alpha[11]*fUpwind[52]+0.1767766952966368*alpha[0]*fUpwind[51]+0.1767766952966368*alpha[8]*fUpwind[50]+0.1767766952966368*alpha[7]*fUpwind[49]+0.1767766952966368*alpha[3]*fUpwind[48]+0.1581138830084189*alpha[20]*fUpwind[27]+0.1581138830084189*alpha[13]*fUpwind[22]+0.1581138830084189*fUpwind[13]*alpha[22]+0.1581138830084189*alpha[12]*fUpwind[21]+0.1581138830084189*alpha[5]*fUpwind[14]+0.1581138830084189*fUpwind[5]*alpha[14]; 
  Ghat[52] = 0.1767766952966368*alpha[8]*fUpwind[62]+0.1767766952966368*alpha[7]*fUpwind[61]+0.1767766952966368*alpha[6]*fUpwind[60]+0.1767766952966368*alpha[3]*fUpwind[58]+0.1767766952966368*alpha[2]*fUpwind[57]+0.1767766952966368*alpha[1]*fUpwind[56]+0.1767766952966368*alpha[19]*fUpwind[55]+0.1767766952966368*alpha[18]*fUpwind[54]+0.1767766952966368*alpha[0]*fUpwind[52]+0.1767766952966368*alpha[11]*fUpwind[51]+0.1767766952966368*alpha[10]*fUpwind[50]+0.1767766952966368*alpha[9]*fUpwind[49]+0.1767766952966368*alpha[4]*fUpwind[48]+0.1581138830084189*alpha[22]*fUpwind[30]+0.1581138830084189*alpha[20]*fUpwind[28]+0.1581138830084189*alpha[14]*fUpwind[25]+0.1581138830084189*alpha[13]*fUpwind[24]+0.1581138830084189*alpha[12]*fUpwind[23]+0.1581138830084189*alpha[5]*fUpwind[15]; 
  Ghat[53] = 0.1767766952966368*alpha[11]*fUpwind[63]+0.1767766952966368*alpha[18]*fUpwind[62]+0.1767766952966368*alpha[19]*fUpwind[61]+0.1767766952966368*alpha[4]*fUpwind[60]+0.1767766952966368*alpha[3]*fUpwind[59]+0.1767766952966368*alpha[9]*fUpwind[57]+0.1767766952966368*alpha[10]*fUpwind[56]+0.1767766952966368*alpha[7]*fUpwind[55]+0.1767766952966368*alpha[8]*fUpwind[54]+0.1767766952966368*alpha[0]*fUpwind[53]+0.1767766952966368*alpha[1]*fUpwind[50]+0.1767766952966368*alpha[2]*fUpwind[49]+0.1767766952966368*alpha[6]*fUpwind[48]+0.1581138830084189*alpha[14]*fUpwind[27]+0.1581138830084189*fUpwind[21]*alpha[22]+0.1581138830084189*alpha[5]*fUpwind[20]+0.1581138830084189*fUpwind[5]*alpha[20]+0.1581138830084189*alpha[12]*fUpwind[13]+0.1581138830084189*fUpwind[12]*alpha[13]; 
  Ghat[54] = 0.1767766952966368*alpha[10]*fUpwind[63]+0.1767766952966368*alpha[4]*fUpwind[61]+0.1767766952966368*alpha[19]*fUpwind[60]+0.1767766952966368*alpha[2]*fUpwind[59]+0.1767766952966368*alpha[9]*fUpwind[58]+0.1767766952966368*alpha[11]*fUpwind[56]+0.1767766952966368*alpha[6]*fUpwind[55]+0.1767766952966368*alpha[0]*fUpwind[54]+0.1767766952966368*alpha[8]*fUpwind[53]+0.1767766952966368*alpha[18]*fUpwind[52]+0.1767766952966368*alpha[1]*fUpwind[51]+0.1767766952966368*alpha[3]*fUpwind[49]+0.1767766952966368*alpha[7]*fUpwind[48]+0.1581138830084189*alpha[13]*fUpwind[27]+0.1581138830084189*alpha[20]*fUpwind[22]+0.1581138830084189*fUpwind[20]*alpha[22]+0.1581138830084189*alpha[5]*fUpwind[21]+0.1581138830084189*alpha[12]*fUpwind[14]+0.1581138830084189*fUpwind[12]*alpha[14]; 
  Ghat[55] = 0.1767766952966368*alpha[9]*fUpwind[63]+0.1767766952966368*alpha[4]*fUpwind[62]+0.1767766952966368*alpha[18]*fUpwind[60]+0.1767766952966368*alpha[1]*fUpwind[59]+0.1767766952966368*alpha[10]*fUpwind[58]+0.1767766952966368*alpha[11]*fUpwind[57]+0.1767766952966368*alpha[0]*fUpwind[55]+0.1767766952966368*alpha[6]*fUpwind[54]+0.1767766952966368*alpha[7]*fUpwind[53]+0.1767766952966368*alpha[19]*fUpwind[52]+0.1767766952966368*alpha[2]*fUpwind[51]+0.1767766952966368*alpha[3]*fUpwind[50]+0.1767766952966368*alpha[8]*fUpwind[48]+0.1581138830084189*alpha[12]*fUpwind[27]+0.1581138830084189*alpha[5]*fUpwind[22]+0.1581138830084189*fUpwind[5]*alpha[22]+0.1581138830084189*alpha[20]*fUpwind[21]+0.1581138830084189*alpha[13]*fUpwind[14]+0.1581138830084189*fUpwind[13]*alpha[14]; 
  Ghat[56] = 0.1767766952966368*alpha[8]*fUpwind[63]+0.1767766952966368*alpha[3]*fUpwind[61]+0.1767766952966368*alpha[2]*fUpwind[60]+0.1767766952966368*alpha[19]*fUpwind[59]+0.1767766952966368*alpha[7]*fUpwind[58]+0.1767766952966368*alpha[6]*fUpwind[57]+0.1767766952966368*alpha[0]*fUpwind[56]+0.1767766952966368*alpha[11]*fUpwind[54]+0.1767766952966368*alpha[10]*fUpwind[53]+0.1767766952966368*alpha[1]*fUpwind[52]+0.1767766952966368*alpha[18]*fUpwind[51]+0.1767766952966368*alpha[4]*fUpwind[49]+0.1767766952966368*alpha[9]*fUpwind[48]+0.1581138830084189*alpha[22]*fUpwind[31]+0.1581138830084189*alpha[14]*fUpwind[29]+0.1581138830084189*alpha[13]*fUpwind[28]+0.1581138830084189*alpha[20]*fUpwind[24]+0.1581138830084189*alpha[5]*fUpwind[23]+0.1581138830084189*alpha[12]*fUpwind[15]; 
  Ghat[57] = 0.1767766952966368*alpha[7]*fUpwind[63]+0.1767766952966368*alpha[3]*fUpwind[62]+0.1767766952966368*alpha[1]*fUpwind[60]+0.1767766952966368*alpha[18]*fUpwind[59]+0.1767766952966368*alpha[8]*fUpwind[58]+0.1767766952966368*alpha[0]*fUpwind[57]+0.1767766952966368*alpha[6]*fUpwind[56]+0.1767766952966368*alpha[11]*fUpwind[55]+0.1767766952966368*alpha[9]*fUpwind[53]+0.1767766952966368*alpha[2]*fUpwind[52]+0.1767766952966368*alpha[19]*fUpwind[51]+0.1767766952966368*alpha[4]*fUpwind[50]+0.1767766952966368*alpha[10]*fUpwind[48]+0.1581138830084189*alpha[14]*fUpwind[30]+0.1581138830084189*alpha[12]*fUpwind[28]+0.1581138830084189*alpha[22]*fUpwind[25]+0.1581138830084189*alpha[5]*fUpwind[24]+0.1581138830084189*alpha[20]*fUpwind[23]+0.1581138830084189*alpha[13]*fUpwind[15]; 
  Ghat[58] = 0.1767766952966368*alpha[6]*fUpwind[63]+0.1767766952966368*alpha[2]*fUpwind[62]+0.1767766952966368*alpha[1]*fUpwind[61]+0.1767766952966368*alpha[0]*fUpwind[58]+0.1767766952966368*alpha[8]*fUpwind[57]+0.1767766952966368*alpha[7]*fUpwind[56]+0.1767766952966368*alpha[10]*fUpwind[55]+0.1767766952966368*alpha[9]*fUpwind[54]+0.1767766952966368*alpha[3]*fUpwind[52]+0.1767766952966368*alpha[4]*fUpwind[51]+0.1767766952966368*alpha[19]*fUpwind[50]+0.1767766952966368*alpha[18]*fUpwind[49]+0.1767766952966368*alpha[11]*fUpwind[48]+0.1581138830084189*alpha[20]*fUpwind[31]+0.1581138830084189*alpha[13]*fUpwind[30]+0.1581138830084189*alpha[12]*fUpwind[29]+0.1581138830084189*alpha[5]*fUpwind[25]+0.1581138830084189*alpha[22]*fUpwind[24]+0.1581138830084189*alpha[14]*fUpwind[15]; 
  Ghat[59] = 0.1767766952966368*alpha[4]*fUpwind[63]+0.1767766952966368*alpha[9]*fUpwind[62]+0.1767766952966368*alpha[10]*fUpwind[61]+0.1767766952966368*alpha[11]*fUpwind[60]+0.1767766952966368*alpha[0]*fUpwind[59]+0.1767766952966368*alpha[18]*fUpwind[57]+0.1767766952966368*alpha[19]*fUpwind[56]+0.1767766952966368*alpha[1]*fUpwind[55]+0.1767766952966368*alpha[2]*fUpwind[54]+0.1767766952966368*alpha[3]*fUpwind[53]+0.1767766952966368*alpha[6]*fUpwind[51]+0.1767766952966368*alpha[7]*fUpwind[50]+0.1767766952966368*alpha[8]*fUpwind[49]+0.1581138830084189*alpha[5]*fUpwind[27]+0.1581138830084189*alpha[12]*fUpwind[22]+0.1581138830084189*fUpwind[12]*alpha[22]+0.1581138830084189*alpha[13]*fUpwind[21]+0.1581138830084189*alpha[14]*fUpwind[20]+0.1581138830084189*fUpwind[14]*alpha[20]; 
  Ghat[60] = 0.1767766952966368*alpha[3]*fUpwind[63]+0.1767766952966368*alpha[7]*fUpwind[62]+0.1767766952966368*alpha[8]*fUpwind[61]+0.1767766952966368*alpha[0]*fUpwind[60]+0.1767766952966368*alpha[11]*fUpwind[59]+0.1767766952966368*alpha[1]*fUpwind[57]+0.1767766952966368*alpha[2]*fUpwind[56]+0.1767766952966368*alpha[18]*fUpwind[55]+0.1767766952966368*alpha[19]*fUpwind[54]+0.1767766952966368*alpha[4]*fUpwind[53]+0.1767766952966368*alpha[6]*fUpwind[52]+0.1767766952966368*alpha[9]*fUpwind[50]+0.1767766952966368*alpha[10]*fUpwind[49]+0.1581138830084189*alpha[14]*fUpwind[31]+0.1581138830084189*alpha[22]*fUpwind[29]+0.1581138830084189*alpha[5]*fUpwind[28]+0.1581138830084189*alpha[12]*fUpwind[24]+0.1581138830084189*alpha[13]*fUpwind[23]+0.1581138830084189*fUpwind[15]*alpha[20]; 
  Ghat[61] = 0.1767766952966368*alpha[2]*fUpwind[63]+0.1767766952966368*alpha[6]*fUpwind[62]+0.1767766952966368*alpha[0]*fUpwind[61]+0.1767766952966368*alpha[8]*fUpwind[60]+0.1767766952966368*alpha[10]*fUpwind[59]+0.1767766952966368*alpha[1]*fUpwind[58]+0.1767766952966368*alpha[3]*fUpwind[56]+0.1767766952966368*alpha[4]*fUpwind[54]+0.1767766952966368*alpha[19]*fUpwind[53]+0.1767766952966368*alpha[7]*fUpwind[52]+0.1767766952966368*alpha[9]*fUpwind[51]+0.1767766952966368*alpha[11]*fUpwind[49]+0.1767766952966368*alpha[18]*fUpwind[48]+0.1581138830084189*alpha[13]*fUpwind[31]+0.1581138830084189*alpha[20]*fUpwind[30]+0.1581138830084189*alpha[5]*fUpwind[29]+0.1581138830084189*alpha[22]*fUpwind[28]+0.1581138830084189*alpha[12]*fUpwind[25]+0.1581138830084189*alpha[14]*fUpwind[23]; 
  Ghat[62] = 0.1767766952966368*alpha[1]*fUpwind[63]+0.1767766952966368*alpha[0]*fUpwind[62]+0.1767766952966368*alpha[6]*fUpwind[61]+0.1767766952966368*alpha[7]*fUpwind[60]+0.1767766952966368*alpha[9]*fUpwind[59]+0.1767766952966368*alpha[2]*fUpwind[58]+0.1767766952966368*alpha[3]*fUpwind[57]+0.1767766952966368*alpha[4]*fUpwind[55]+0.1767766952966368*alpha[18]*fUpwind[53]+0.1767766952966368*alpha[8]*fUpwind[52]+0.1767766952966368*alpha[10]*fUpwind[51]+0.1767766952966368*alpha[11]*fUpwind[50]+0.1767766952966368*alpha[19]*fUpwind[48]+0.1581138830084189*alpha[12]*fUpwind[31]+0.1581138830084189*alpha[5]*fUpwind[30]+0.1581138830084189*alpha[20]*fUpwind[29]+0.1581138830084189*alpha[13]*fUpwind[25]+0.1581138830084189*alpha[14]*fUpwind[24]+0.1581138830084189*fUpwind[15]*alpha[22]; 
  Ghat[63] = 0.1767766952966368*alpha[0]*fUpwind[63]+0.1767766952966368*alpha[1]*fUpwind[62]+0.1767766952966368*alpha[2]*fUpwind[61]+0.1767766952966368*alpha[3]*fUpwind[60]+0.1767766952966368*alpha[4]*fUpwind[59]+0.1767766952966368*alpha[6]*fUpwind[58]+0.1767766952966368*alpha[7]*fUpwind[57]+0.1767766952966368*alpha[8]*fUpwind[56]+0.1767766952966368*alpha[9]*fUpwind[55]+0.1767766952966368*alpha[10]*fUpwind[54]+0.1767766952966368*alpha[11]*fUpwind[53]+0.1767766952966368*alpha[18]*fUpwind[50]+0.1767766952966368*alpha[19]*fUpwind[49]+0.1581138830084189*alpha[5]*fUpwind[31]+0.1581138830084189*alpha[12]*fUpwind[30]+0.1581138830084189*alpha[13]*fUpwind[29]+0.1581138830084189*alpha[14]*fUpwind[28]+0.1581138830084189*alpha[20]*fUpwind[25]+0.1581138830084189*alpha[22]*fUpwind[23]; 

  out[0] += -0.7071067811865475*Ghat[0]*dv10; 
  out[1] += -0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -0.7071067811865475*Ghat[2]*dv10; 
  out[3] += -0.7071067811865475*Ghat[3]*dv10; 
  out[4] += -1.224744871391589*Ghat[0]*dv10; 
  out[5] += -0.7071067811865475*Ghat[4]*dv10; 
  out[6] += -0.7071067811865475*Ghat[5]*dv10; 
  out[7] += -0.7071067811865475*Ghat[6]*dv10; 
  out[8] += -0.7071067811865475*Ghat[7]*dv10; 
  out[9] += -0.7071067811865475*Ghat[8]*dv10; 
  out[10] += -1.224744871391589*Ghat[1]*dv10; 
  out[11] += -1.224744871391589*Ghat[2]*dv10; 
  out[12] += -1.224744871391589*Ghat[3]*dv10; 
  out[13] += -0.7071067811865475*Ghat[9]*dv10; 
  out[14] += -0.7071067811865475*Ghat[10]*dv10; 
  out[15] += -0.7071067811865475*Ghat[11]*dv10; 
  out[16] += -1.224744871391589*Ghat[4]*dv10; 
  out[17] += -0.7071067811865475*Ghat[12]*dv10; 
  out[18] += -0.7071067811865475*Ghat[13]*dv10; 
  out[19] += -0.7071067811865475*Ghat[14]*dv10; 
  out[20] += -1.224744871391589*Ghat[5]*dv10; 
  out[21] += -0.7071067811865475*Ghat[15]*dv10; 
  out[22] += -0.7071067811865475*Ghat[16]*dv10; 
  out[23] += -1.224744871391589*Ghat[6]*dv10; 
  out[24] += -1.224744871391589*Ghat[7]*dv10; 
  out[25] += -1.224744871391589*Ghat[8]*dv10; 
  out[26] += -0.7071067811865475*Ghat[17]*dv10; 
  out[27] += -0.7071067811865475*Ghat[18]*dv10; 
  out[28] += -0.7071067811865475*Ghat[19]*dv10; 
  out[29] += -1.224744871391589*Ghat[9]*dv10; 
  out[30] += -1.224744871391589*Ghat[10]*dv10; 
  out[31] += -1.224744871391589*Ghat[11]*dv10; 
  out[32] += -0.7071067811865475*Ghat[20]*dv10; 
  out[33] += -0.7071067811865475*Ghat[21]*dv10; 
  out[34] += -0.7071067811865475*Ghat[22]*dv10; 
  out[35] += -1.224744871391589*Ghat[12]*dv10; 
  out[36] += -1.224744871391589*Ghat[13]*dv10; 
  out[37] += -1.224744871391589*Ghat[14]*dv10; 
  out[38] += -0.7071067811865475*Ghat[23]*dv10; 
  out[39] += -0.7071067811865475*Ghat[24]*dv10; 
  out[40] += -0.7071067811865475*Ghat[25]*dv10; 
  out[41] += -1.224744871391589*Ghat[15]*dv10; 
  out[42] += -1.224744871391589*Ghat[16]*dv10; 
  out[43] += -0.7071067811865475*Ghat[26]*dv10; 
  out[44] += -1.224744871391589*Ghat[17]*dv10; 
  out[45] += -1.224744871391589*Ghat[18]*dv10; 
  out[46] += -1.224744871391589*Ghat[19]*dv10; 
  out[47] += -0.7071067811865475*Ghat[27]*dv10; 
  out[48] += -1.224744871391589*Ghat[20]*dv10; 
  out[49] += -1.224744871391589*Ghat[21]*dv10; 
  out[50] += -1.224744871391589*Ghat[22]*dv10; 
  out[51] += -0.7071067811865475*Ghat[28]*dv10; 
  out[52] += -0.7071067811865475*Ghat[29]*dv10; 
  out[53] += -0.7071067811865475*Ghat[30]*dv10; 
  out[54] += -1.224744871391589*Ghat[23]*dv10; 
  out[55] += -1.224744871391589*Ghat[24]*dv10; 
  out[56] += -1.224744871391589*Ghat[25]*dv10; 
  out[57] += -1.224744871391589*Ghat[26]*dv10; 
  out[58] += -1.224744871391589*Ghat[27]*dv10; 
  out[59] += -0.7071067811865475*Ghat[31]*dv10; 
  out[60] += -1.224744871391589*Ghat[28]*dv10; 
  out[61] += -1.224744871391589*Ghat[29]*dv10; 
  out[62] += -1.224744871391589*Ghat[30]*dv10; 
  out[63] += -1.224744871391589*Ghat[31]*dv10; 
  out[64] += -1.58113883008419*Ghat[0]*dv10; 
  out[65] += -1.58113883008419*Ghat[1]*dv10; 
  out[66] += -1.58113883008419*Ghat[2]*dv10; 
  out[67] += -1.58113883008419*Ghat[3]*dv10; 
  out[68] += -1.58113883008419*Ghat[4]*dv10; 
  out[69] += -1.58113883008419*Ghat[5]*dv10; 
  out[70] += -1.58113883008419*Ghat[6]*dv10; 
  out[71] += -1.58113883008419*Ghat[7]*dv10; 
  out[72] += -1.58113883008419*Ghat[8]*dv10; 
  out[73] += -1.58113883008419*Ghat[9]*dv10; 
  out[74] += -1.58113883008419*Ghat[10]*dv10; 
  out[75] += -1.58113883008419*Ghat[11]*dv10; 
  out[76] += -1.58113883008419*Ghat[12]*dv10; 
  out[77] += -1.58113883008419*Ghat[13]*dv10; 
  out[78] += -1.58113883008419*Ghat[14]*dv10; 
  out[79] += -1.58113883008419*Ghat[15]*dv10; 
  out[80] += -1.58113883008419*Ghat[16]*dv10; 
  out[81] += -1.58113883008419*Ghat[17]*dv10; 
  out[82] += -1.58113883008419*Ghat[18]*dv10; 
  out[83] += -1.58113883008419*Ghat[19]*dv10; 
  out[84] += -1.58113883008419*Ghat[20]*dv10; 
  out[85] += -1.58113883008419*Ghat[21]*dv10; 
  out[86] += -1.58113883008419*Ghat[22]*dv10; 
  out[87] += -1.58113883008419*Ghat[23]*dv10; 
  out[88] += -1.58113883008419*Ghat[24]*dv10; 
  out[89] += -1.58113883008419*Ghat[25]*dv10; 
  out[90] += -1.58113883008419*Ghat[26]*dv10; 
  out[91] += -1.58113883008419*Ghat[27]*dv10; 
  out[92] += -1.58113883008419*Ghat[28]*dv10; 
  out[93] += -1.58113883008419*Ghat[29]*dv10; 
  out[94] += -1.58113883008419*Ghat[30]*dv10; 
  out[95] += -1.58113883008419*Ghat[31]*dv10; 
  out[96] += -0.7071067811865475*Ghat[32]*dv10; 
  out[97] += -0.7071067811865475*Ghat[33]*dv10; 
  out[98] += -0.7071067811865475*Ghat[34]*dv10; 
  out[99] += -0.7071067811865475*Ghat[35]*dv10; 
  out[100] += -1.224744871391589*Ghat[32]*dv10; 
  out[101] += -0.7071067811865475*Ghat[36]*dv10; 
  out[102] += -0.7071067811865475*Ghat[37]*dv10; 
  out[103] += -0.7071067811865475*Ghat[38]*dv10; 
  out[104] += -0.7071067811865475*Ghat[39]*dv10; 
  out[105] += -1.224744871391589*Ghat[33]*dv10; 
  out[106] += -1.224744871391589*Ghat[34]*dv10; 
  out[107] += -1.224744871391589*Ghat[35]*dv10; 
  out[108] += -0.7071067811865475*Ghat[40]*dv10; 
  out[109] += -0.7071067811865475*Ghat[41]*dv10; 
  out[110] += -0.7071067811865475*Ghat[42]*dv10; 
  out[111] += -1.224744871391589*Ghat[36]*dv10; 
  out[112] += -0.7071067811865475*Ghat[43]*dv10; 
  out[113] += -1.224744871391589*Ghat[37]*dv10; 
  out[114] += -1.224744871391589*Ghat[38]*dv10; 
  out[115] += -1.224744871391589*Ghat[39]*dv10; 
  out[116] += -0.7071067811865475*Ghat[44]*dv10; 
  out[117] += -0.7071067811865475*Ghat[45]*dv10; 
  out[118] += -0.7071067811865475*Ghat[46]*dv10; 
  out[119] += -1.224744871391589*Ghat[40]*dv10; 
  out[120] += -1.224744871391589*Ghat[41]*dv10; 
  out[121] += -1.224744871391589*Ghat[42]*dv10; 
  out[122] += -1.224744871391589*Ghat[43]*dv10; 
  out[123] += -0.7071067811865475*Ghat[47]*dv10; 
  out[124] += -1.224744871391589*Ghat[44]*dv10; 
  out[125] += -1.224744871391589*Ghat[45]*dv10; 
  out[126] += -1.224744871391589*Ghat[46]*dv10; 
  out[127] += -1.224744871391589*Ghat[47]*dv10; 
  out[128] += -0.7071067811865475*Ghat[48]*dv10; 
  out[129] += -0.7071067811865475*Ghat[49]*dv10; 
  out[130] += -0.7071067811865475*Ghat[50]*dv10; 
  out[131] += -0.7071067811865475*Ghat[51]*dv10; 
  out[132] += -1.224744871391589*Ghat[48]*dv10; 
  out[133] += -0.7071067811865475*Ghat[52]*dv10; 
  out[134] += -0.7071067811865475*Ghat[53]*dv10; 
  out[135] += -0.7071067811865475*Ghat[54]*dv10; 
  out[136] += -0.7071067811865475*Ghat[55]*dv10; 
  out[137] += -1.224744871391589*Ghat[49]*dv10; 
  out[138] += -1.224744871391589*Ghat[50]*dv10; 
  out[139] += -1.224744871391589*Ghat[51]*dv10; 
  out[140] += -0.7071067811865475*Ghat[56]*dv10; 
  out[141] += -0.7071067811865475*Ghat[57]*dv10; 
  out[142] += -0.7071067811865475*Ghat[58]*dv10; 
  out[143] += -1.224744871391589*Ghat[52]*dv10; 
  out[144] += -0.7071067811865475*Ghat[59]*dv10; 
  out[145] += -1.224744871391589*Ghat[53]*dv10; 
  out[146] += -1.224744871391589*Ghat[54]*dv10; 
  out[147] += -1.224744871391589*Ghat[55]*dv10; 
  out[148] += -0.7071067811865475*Ghat[60]*dv10; 
  out[149] += -0.7071067811865475*Ghat[61]*dv10; 
  out[150] += -0.7071067811865475*Ghat[62]*dv10; 
  out[151] += -1.224744871391589*Ghat[56]*dv10; 
  out[152] += -1.224744871391589*Ghat[57]*dv10; 
  out[153] += -1.224744871391589*Ghat[58]*dv10; 
  out[154] += -1.224744871391589*Ghat[59]*dv10; 
  out[155] += -0.7071067811865475*Ghat[63]*dv10; 
  out[156] += -1.224744871391589*Ghat[60]*dv10; 
  out[157] += -1.224744871391589*Ghat[61]*dv10; 
  out[158] += -1.224744871391589*Ghat[62]*dv10; 
  out[159] += -1.224744871391589*Ghat[63]*dv10; 

  } else { 

  if ((-0.2371708245126284*(alpha[22]+alpha[20]+alpha[19]+alpha[18]))+0.2371708245126284*(alpha[14]+alpha[13]+alpha[12]+alpha[11]+alpha[10]+alpha[9])+0.1767766952966368*(alpha[8]+alpha[7]+alpha[6])-0.2371708245126284*(alpha[5]+alpha[4])-0.1767766952966368*(alpha[3]+alpha[2]+alpha[1])+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[0] = hyb_3x3v_p1_surfx4_eval_quad_node_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = hyb_3x3v_p1_surfx4_eval_quad_node_0_l(fSkin); 
  } 
  if ((-0.2371708245126284*(alpha[19]+alpha[18]))+0.2371708245126284*(alpha[11]+alpha[10]+alpha[9])+0.1767766952966368*(alpha[8]+alpha[7]+alpha[6])-0.2371708245126284*alpha[4]-0.1767766952966368*(alpha[3]+alpha[2]+alpha[1])+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[1] = hyb_3x3v_p1_surfx4_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = hyb_3x3v_p1_surfx4_eval_quad_node_1_l(fSkin); 
  } 
  if (0.2371708245126284*(alpha[22]+alpha[20])-0.2371708245126284*(alpha[19]+alpha[18]+alpha[14]+alpha[13]+alpha[12])+0.2371708245126284*(alpha[11]+alpha[10]+alpha[9])+0.1767766952966368*(alpha[8]+alpha[7]+alpha[6])+0.2371708245126284*alpha[5]-0.2371708245126284*alpha[4]-0.1767766952966368*(alpha[3]+alpha[2]+alpha[1])+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[2] = hyb_3x3v_p1_surfx4_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[2] = hyb_3x3v_p1_surfx4_eval_quad_node_2_l(fSkin); 
  } 
  if ((-0.2371708245126284*(alpha[22]+alpha[20]))+0.2371708245126284*(alpha[14]+alpha[13]+alpha[12])+0.1767766952966368*(alpha[8]+alpha[7]+alpha[6])-0.2371708245126284*alpha[5]-0.1767766952966368*(alpha[3]+alpha[2]+alpha[1])+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[3] = hyb_3x3v_p1_surfx4_eval_quad_node_3_r(fEdge); 
  } else { 
    fUpwindQuad[3] = hyb_3x3v_p1_surfx4_eval_quad_node_3_l(fSkin); 
  } 
  if (0.1767766952966368*(alpha[8]+alpha[7]+alpha[6])-0.1767766952966368*(alpha[3]+alpha[2]+alpha[1])+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[4] = hyb_3x3v_p1_surfx4_eval_quad_node_4_r(fEdge); 
  } else { 
    fUpwindQuad[4] = hyb_3x3v_p1_surfx4_eval_quad_node_4_l(fSkin); 
  } 
  if (0.2371708245126284*(alpha[22]+alpha[20])-0.2371708245126284*(alpha[14]+alpha[13]+alpha[12])+0.1767766952966368*(alpha[8]+alpha[7]+alpha[6])+0.2371708245126284*alpha[5]-0.1767766952966368*(alpha[3]+alpha[2]+alpha[1])+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[5] = hyb_3x3v_p1_surfx4_eval_quad_node_5_r(fEdge); 
  } else { 
    fUpwindQuad[5] = hyb_3x3v_p1_surfx4_eval_quad_node_5_l(fSkin); 
  } 
  if ((-0.2371708245126284*(alpha[22]+alpha[20]))+0.2371708245126284*(alpha[19]+alpha[18]+alpha[14]+alpha[13]+alpha[12])-0.2371708245126284*(alpha[11]+alpha[10]+alpha[9])+0.1767766952966368*(alpha[8]+alpha[7]+alpha[6])-0.2371708245126284*alpha[5]+0.2371708245126284*alpha[4]-0.1767766952966368*(alpha[3]+alpha[2]+alpha[1])+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[6] = hyb_3x3v_p1_surfx4_eval_quad_node_6_r(fEdge); 
  } else { 
    fUpwindQuad[6] = hyb_3x3v_p1_surfx4_eval_quad_node_6_l(fSkin); 
  } 
  if (0.2371708245126284*(alpha[19]+alpha[18])-0.2371708245126284*(alpha[11]+alpha[10]+alpha[9])+0.1767766952966368*(alpha[8]+alpha[7]+alpha[6])+0.2371708245126284*alpha[4]-0.1767766952966368*(alpha[3]+alpha[2]+alpha[1])+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[7] = hyb_3x3v_p1_surfx4_eval_quad_node_7_r(fEdge); 
  } else { 
    fUpwindQuad[7] = hyb_3x3v_p1_surfx4_eval_quad_node_7_l(fSkin); 
  } 
  if (0.2371708245126284*(alpha[22]+alpha[20]+alpha[19]+alpha[18])-0.2371708245126284*(alpha[14]+alpha[13]+alpha[12]+alpha[11]+alpha[10]+alpha[9])+0.1767766952966368*(alpha[8]+alpha[7]+alpha[6])+0.2371708245126284*(alpha[5]+alpha[4])-0.1767766952966368*(alpha[3]+alpha[2]+alpha[1])+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[8] = hyb_3x3v_p1_surfx4_eval_quad_node_8_r(fEdge); 
  } else { 
    fUpwindQuad[8] = hyb_3x3v_p1_surfx4_eval_quad_node_8_l(fSkin); 
  } 
  if (0.2371708245126284*alpha[22]-0.2371708245126284*alpha[20]+0.2371708245126284*(alpha[19]+alpha[18])-0.2371708245126284*alpha[14]+0.2371708245126284*(alpha[13]+alpha[12])-0.2371708245126284*alpha[11]+0.2371708245126284*(alpha[10]+alpha[9])-0.1767766952966368*(alpha[8]+alpha[7])+0.1767766952966368*alpha[6]-0.2371708245126284*(alpha[5]+alpha[4])+0.1767766952966368*alpha[3]-0.1767766952966368*(alpha[2]+alpha[1])+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[9] = hyb_3x3v_p1_surfx4_eval_quad_node_9_r(fEdge); 
  } else { 
    fUpwindQuad[9] = hyb_3x3v_p1_surfx4_eval_quad_node_9_l(fSkin); 
  } 
  if (0.2371708245126284*(alpha[19]+alpha[18])-0.2371708245126284*alpha[11]+0.2371708245126284*(alpha[10]+alpha[9])-0.1767766952966368*(alpha[8]+alpha[7])+0.1767766952966368*alpha[6]-0.2371708245126284*alpha[4]+0.1767766952966368*alpha[3]-0.1767766952966368*(alpha[2]+alpha[1])+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[10] = hyb_3x3v_p1_surfx4_eval_quad_node_10_r(fEdge); 
  } else { 
    fUpwindQuad[10] = hyb_3x3v_p1_surfx4_eval_quad_node_10_l(fSkin); 
  } 
  if ((-0.2371708245126284*alpha[22])+0.2371708245126284*(alpha[20]+alpha[19]+alpha[18]+alpha[14])-0.2371708245126284*(alpha[13]+alpha[12]+alpha[11])+0.2371708245126284*(alpha[10]+alpha[9])-0.1767766952966368*(alpha[8]+alpha[7])+0.1767766952966368*alpha[6]+0.2371708245126284*alpha[5]-0.2371708245126284*alpha[4]+0.1767766952966368*alpha[3]-0.1767766952966368*(alpha[2]+alpha[1])+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[11] = hyb_3x3v_p1_surfx4_eval_quad_node_11_r(fEdge); 
  } else { 
    fUpwindQuad[11] = hyb_3x3v_p1_surfx4_eval_quad_node_11_l(fSkin); 
  } 
  if (0.2371708245126284*alpha[22]-0.2371708245126284*(alpha[20]+alpha[14])+0.2371708245126284*(alpha[13]+alpha[12])-0.1767766952966368*(alpha[8]+alpha[7])+0.1767766952966368*alpha[6]-0.2371708245126284*alpha[5]+0.1767766952966368*alpha[3]-0.1767766952966368*(alpha[2]+alpha[1])+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[12] = hyb_3x3v_p1_surfx4_eval_quad_node_12_r(fEdge); 
  } else { 
    fUpwindQuad[12] = hyb_3x3v_p1_surfx4_eval_quad_node_12_l(fSkin); 
  } 
  if ((-0.1767766952966368*(alpha[8]+alpha[7]))+0.1767766952966368*(alpha[6]+alpha[3])-0.1767766952966368*(alpha[2]+alpha[1])+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[13] = hyb_3x3v_p1_surfx4_eval_quad_node_13_r(fEdge); 
  } else { 
    fUpwindQuad[13] = hyb_3x3v_p1_surfx4_eval_quad_node_13_l(fSkin); 
  } 
  if ((-0.2371708245126284*alpha[22])+0.2371708245126284*(alpha[20]+alpha[14])-0.2371708245126284*(alpha[13]+alpha[12])-0.1767766952966368*(alpha[8]+alpha[7])+0.1767766952966368*alpha[6]+0.2371708245126284*alpha[5]+0.1767766952966368*alpha[3]-0.1767766952966368*(alpha[2]+alpha[1])+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[14] = hyb_3x3v_p1_surfx4_eval_quad_node_14_r(fEdge); 
  } else { 
    fUpwindQuad[14] = hyb_3x3v_p1_surfx4_eval_quad_node_14_l(fSkin); 
  } 
  if (0.2371708245126284*alpha[22]-0.2371708245126284*(alpha[20]+alpha[19]+alpha[18]+alpha[14])+0.2371708245126284*(alpha[13]+alpha[12]+alpha[11])-0.2371708245126284*(alpha[10]+alpha[9])-0.1767766952966368*(alpha[8]+alpha[7])+0.1767766952966368*alpha[6]-0.2371708245126284*alpha[5]+0.2371708245126284*alpha[4]+0.1767766952966368*alpha[3]-0.1767766952966368*(alpha[2]+alpha[1])+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[15] = hyb_3x3v_p1_surfx4_eval_quad_node_15_r(fEdge); 
  } else { 
    fUpwindQuad[15] = hyb_3x3v_p1_surfx4_eval_quad_node_15_l(fSkin); 
  } 
  if ((-0.2371708245126284*(alpha[19]+alpha[18]))+0.2371708245126284*alpha[11]-0.2371708245126284*(alpha[10]+alpha[9])-0.1767766952966368*(alpha[8]+alpha[7])+0.1767766952966368*alpha[6]+0.2371708245126284*alpha[4]+0.1767766952966368*alpha[3]-0.1767766952966368*(alpha[2]+alpha[1])+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[16] = hyb_3x3v_p1_surfx4_eval_quad_node_16_r(fEdge); 
  } else { 
    fUpwindQuad[16] = hyb_3x3v_p1_surfx4_eval_quad_node_16_l(fSkin); 
  } 
  if ((-0.2371708245126284*alpha[22])+0.2371708245126284*alpha[20]-0.2371708245126284*(alpha[19]+alpha[18])+0.2371708245126284*alpha[14]-0.2371708245126284*(alpha[13]+alpha[12])+0.2371708245126284*alpha[11]-0.2371708245126284*(alpha[10]+alpha[9])-0.1767766952966368*(alpha[8]+alpha[7])+0.1767766952966368*alpha[6]+0.2371708245126284*(alpha[5]+alpha[4])+0.1767766952966368*alpha[3]-0.1767766952966368*(alpha[2]+alpha[1])+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[17] = hyb_3x3v_p1_surfx4_eval_quad_node_17_r(fEdge); 
  } else { 
    fUpwindQuad[17] = hyb_3x3v_p1_surfx4_eval_quad_node_17_l(fSkin); 
  } 
  if (0.2371708245126284*(alpha[22]+alpha[20]+alpha[19])-0.2371708245126284*alpha[18]+0.2371708245126284*alpha[14]-0.2371708245126284*alpha[13]+0.2371708245126284*(alpha[12]+alpha[11])-0.2371708245126284*alpha[10]+0.2371708245126284*alpha[9]-0.1767766952966368*alpha[8]+0.1767766952966368*alpha[7]-0.1767766952966368*alpha[6]-0.2371708245126284*(alpha[5]+alpha[4])-0.1767766952966368*alpha[3]+0.1767766952966368*alpha[2]-0.1767766952966368*alpha[1]+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[18] = hyb_3x3v_p1_surfx4_eval_quad_node_18_r(fEdge); 
  } else { 
    fUpwindQuad[18] = hyb_3x3v_p1_surfx4_eval_quad_node_18_l(fSkin); 
  } 
  if (0.2371708245126284*alpha[19]-0.2371708245126284*alpha[18]+0.2371708245126284*alpha[11]-0.2371708245126284*alpha[10]+0.2371708245126284*alpha[9]-0.1767766952966368*alpha[8]+0.1767766952966368*alpha[7]-0.1767766952966368*alpha[6]-0.2371708245126284*alpha[4]-0.1767766952966368*alpha[3]+0.1767766952966368*alpha[2]-0.1767766952966368*alpha[1]+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[19] = hyb_3x3v_p1_surfx4_eval_quad_node_19_r(fEdge); 
  } else { 
    fUpwindQuad[19] = hyb_3x3v_p1_surfx4_eval_quad_node_19_l(fSkin); 
  } 
  if ((-0.2371708245126284*(alpha[22]+alpha[20]))+0.2371708245126284*alpha[19]-0.2371708245126284*(alpha[18]+alpha[14])+0.2371708245126284*alpha[13]-0.2371708245126284*alpha[12]+0.2371708245126284*alpha[11]-0.2371708245126284*alpha[10]+0.2371708245126284*alpha[9]-0.1767766952966368*alpha[8]+0.1767766952966368*alpha[7]-0.1767766952966368*alpha[6]+0.2371708245126284*alpha[5]-0.2371708245126284*alpha[4]-0.1767766952966368*alpha[3]+0.1767766952966368*alpha[2]-0.1767766952966368*alpha[1]+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[20] = hyb_3x3v_p1_surfx4_eval_quad_node_20_r(fEdge); 
  } else { 
    fUpwindQuad[20] = hyb_3x3v_p1_surfx4_eval_quad_node_20_l(fSkin); 
  } 
  if (0.2371708245126284*(alpha[22]+alpha[20]+alpha[14])-0.2371708245126284*alpha[13]+0.2371708245126284*alpha[12]-0.1767766952966368*alpha[8]+0.1767766952966368*alpha[7]-0.1767766952966368*alpha[6]-0.2371708245126284*alpha[5]-0.1767766952966368*alpha[3]+0.1767766952966368*alpha[2]-0.1767766952966368*alpha[1]+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[21] = hyb_3x3v_p1_surfx4_eval_quad_node_21_r(fEdge); 
  } else { 
    fUpwindQuad[21] = hyb_3x3v_p1_surfx4_eval_quad_node_21_l(fSkin); 
  } 
  if ((-0.1767766952966368*alpha[8])+0.1767766952966368*alpha[7]-0.1767766952966368*(alpha[6]+alpha[3])+0.1767766952966368*alpha[2]-0.1767766952966368*alpha[1]+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[22] = hyb_3x3v_p1_surfx4_eval_quad_node_22_r(fEdge); 
  } else { 
    fUpwindQuad[22] = hyb_3x3v_p1_surfx4_eval_quad_node_22_l(fSkin); 
  } 
  if ((-0.2371708245126284*(alpha[22]+alpha[20]+alpha[14]))+0.2371708245126284*alpha[13]-0.2371708245126284*alpha[12]-0.1767766952966368*alpha[8]+0.1767766952966368*alpha[7]-0.1767766952966368*alpha[6]+0.2371708245126284*alpha[5]-0.1767766952966368*alpha[3]+0.1767766952966368*alpha[2]-0.1767766952966368*alpha[1]+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[23] = hyb_3x3v_p1_surfx4_eval_quad_node_23_r(fEdge); 
  } else { 
    fUpwindQuad[23] = hyb_3x3v_p1_surfx4_eval_quad_node_23_l(fSkin); 
  } 
  if (0.2371708245126284*(alpha[22]+alpha[20])-0.2371708245126284*alpha[19]+0.2371708245126284*(alpha[18]+alpha[14])-0.2371708245126284*alpha[13]+0.2371708245126284*alpha[12]-0.2371708245126284*alpha[11]+0.2371708245126284*alpha[10]-0.2371708245126284*alpha[9]-0.1767766952966368*alpha[8]+0.1767766952966368*alpha[7]-0.1767766952966368*alpha[6]-0.2371708245126284*alpha[5]+0.2371708245126284*alpha[4]-0.1767766952966368*alpha[3]+0.1767766952966368*alpha[2]-0.1767766952966368*alpha[1]+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[24] = hyb_3x3v_p1_surfx4_eval_quad_node_24_r(fEdge); 
  } else { 
    fUpwindQuad[24] = hyb_3x3v_p1_surfx4_eval_quad_node_24_l(fSkin); 
  } 
  if ((-0.2371708245126284*alpha[19])+0.2371708245126284*alpha[18]-0.2371708245126284*alpha[11]+0.2371708245126284*alpha[10]-0.2371708245126284*alpha[9]-0.1767766952966368*alpha[8]+0.1767766952966368*alpha[7]-0.1767766952966368*alpha[6]+0.2371708245126284*alpha[4]-0.1767766952966368*alpha[3]+0.1767766952966368*alpha[2]-0.1767766952966368*alpha[1]+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[25] = hyb_3x3v_p1_surfx4_eval_quad_node_25_r(fEdge); 
  } else { 
    fUpwindQuad[25] = hyb_3x3v_p1_surfx4_eval_quad_node_25_l(fSkin); 
  } 
  if ((-0.2371708245126284*(alpha[22]+alpha[20]+alpha[19]))+0.2371708245126284*alpha[18]-0.2371708245126284*alpha[14]+0.2371708245126284*alpha[13]-0.2371708245126284*(alpha[12]+alpha[11])+0.2371708245126284*alpha[10]-0.2371708245126284*alpha[9]-0.1767766952966368*alpha[8]+0.1767766952966368*alpha[7]-0.1767766952966368*alpha[6]+0.2371708245126284*(alpha[5]+alpha[4])-0.1767766952966368*alpha[3]+0.1767766952966368*alpha[2]-0.1767766952966368*alpha[1]+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[26] = hyb_3x3v_p1_surfx4_eval_quad_node_26_r(fEdge); 
  } else { 
    fUpwindQuad[26] = hyb_3x3v_p1_surfx4_eval_quad_node_26_l(fSkin); 
  } 
  if ((-0.2371708245126284*alpha[22])+0.2371708245126284*alpha[20]-0.2371708245126284*alpha[19]+0.2371708245126284*alpha[18]-0.2371708245126284*(alpha[14]+alpha[13])+0.2371708245126284*alpha[12]-0.2371708245126284*(alpha[11]+alpha[10])+0.2371708245126284*alpha[9]+0.1767766952966368*alpha[8]-0.1767766952966368*(alpha[7]+alpha[6])-0.2371708245126284*(alpha[5]+alpha[4])+0.1767766952966368*(alpha[3]+alpha[2])-0.1767766952966368*alpha[1]+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[27] = hyb_3x3v_p1_surfx4_eval_quad_node_27_r(fEdge); 
  } else { 
    fUpwindQuad[27] = hyb_3x3v_p1_surfx4_eval_quad_node_27_l(fSkin); 
  } 
  if ((-0.2371708245126284*alpha[19])+0.2371708245126284*alpha[18]-0.2371708245126284*(alpha[11]+alpha[10])+0.2371708245126284*alpha[9]+0.1767766952966368*alpha[8]-0.1767766952966368*(alpha[7]+alpha[6])-0.2371708245126284*alpha[4]+0.1767766952966368*(alpha[3]+alpha[2])-0.1767766952966368*alpha[1]+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[28] = hyb_3x3v_p1_surfx4_eval_quad_node_28_r(fEdge); 
  } else { 
    fUpwindQuad[28] = hyb_3x3v_p1_surfx4_eval_quad_node_28_l(fSkin); 
  } 
  if (0.2371708245126284*alpha[22]-0.2371708245126284*(alpha[20]+alpha[19])+0.2371708245126284*(alpha[18]+alpha[14]+alpha[13])-0.2371708245126284*(alpha[12]+alpha[11]+alpha[10])+0.2371708245126284*alpha[9]+0.1767766952966368*alpha[8]-0.1767766952966368*(alpha[7]+alpha[6])+0.2371708245126284*alpha[5]-0.2371708245126284*alpha[4]+0.1767766952966368*(alpha[3]+alpha[2])-0.1767766952966368*alpha[1]+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[29] = hyb_3x3v_p1_surfx4_eval_quad_node_29_r(fEdge); 
  } else { 
    fUpwindQuad[29] = hyb_3x3v_p1_surfx4_eval_quad_node_29_l(fSkin); 
  } 
  if ((-0.2371708245126284*alpha[22])+0.2371708245126284*alpha[20]-0.2371708245126284*(alpha[14]+alpha[13])+0.2371708245126284*alpha[12]+0.1767766952966368*alpha[8]-0.1767766952966368*(alpha[7]+alpha[6])-0.2371708245126284*alpha[5]+0.1767766952966368*(alpha[3]+alpha[2])-0.1767766952966368*alpha[1]+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[30] = hyb_3x3v_p1_surfx4_eval_quad_node_30_r(fEdge); 
  } else { 
    fUpwindQuad[30] = hyb_3x3v_p1_surfx4_eval_quad_node_30_l(fSkin); 
  } 
  if (0.1767766952966368*alpha[8]-0.1767766952966368*(alpha[7]+alpha[6])+0.1767766952966368*(alpha[3]+alpha[2])-0.1767766952966368*alpha[1]+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[31] = hyb_3x3v_p1_surfx4_eval_quad_node_31_r(fEdge); 
  } else { 
    fUpwindQuad[31] = hyb_3x3v_p1_surfx4_eval_quad_node_31_l(fSkin); 
  } 
  if (0.2371708245126284*alpha[22]-0.2371708245126284*alpha[20]+0.2371708245126284*(alpha[14]+alpha[13])-0.2371708245126284*alpha[12]+0.1767766952966368*alpha[8]-0.1767766952966368*(alpha[7]+alpha[6])+0.2371708245126284*alpha[5]+0.1767766952966368*(alpha[3]+alpha[2])-0.1767766952966368*alpha[1]+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[32] = hyb_3x3v_p1_surfx4_eval_quad_node_32_r(fEdge); 
  } else { 
    fUpwindQuad[32] = hyb_3x3v_p1_surfx4_eval_quad_node_32_l(fSkin); 
  } 
  if ((-0.2371708245126284*alpha[22])+0.2371708245126284*(alpha[20]+alpha[19])-0.2371708245126284*(alpha[18]+alpha[14]+alpha[13])+0.2371708245126284*(alpha[12]+alpha[11]+alpha[10])-0.2371708245126284*alpha[9]+0.1767766952966368*alpha[8]-0.1767766952966368*(alpha[7]+alpha[6])-0.2371708245126284*alpha[5]+0.2371708245126284*alpha[4]+0.1767766952966368*(alpha[3]+alpha[2])-0.1767766952966368*alpha[1]+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[33] = hyb_3x3v_p1_surfx4_eval_quad_node_33_r(fEdge); 
  } else { 
    fUpwindQuad[33] = hyb_3x3v_p1_surfx4_eval_quad_node_33_l(fSkin); 
  } 
  if (0.2371708245126284*alpha[19]-0.2371708245126284*alpha[18]+0.2371708245126284*(alpha[11]+alpha[10])-0.2371708245126284*alpha[9]+0.1767766952966368*alpha[8]-0.1767766952966368*(alpha[7]+alpha[6])+0.2371708245126284*alpha[4]+0.1767766952966368*(alpha[3]+alpha[2])-0.1767766952966368*alpha[1]+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[34] = hyb_3x3v_p1_surfx4_eval_quad_node_34_r(fEdge); 
  } else { 
    fUpwindQuad[34] = hyb_3x3v_p1_surfx4_eval_quad_node_34_l(fSkin); 
  } 
  if (0.2371708245126284*alpha[22]-0.2371708245126284*alpha[20]+0.2371708245126284*alpha[19]-0.2371708245126284*alpha[18]+0.2371708245126284*(alpha[14]+alpha[13])-0.2371708245126284*alpha[12]+0.2371708245126284*(alpha[11]+alpha[10])-0.2371708245126284*alpha[9]+0.1767766952966368*alpha[8]-0.1767766952966368*(alpha[7]+alpha[6])+0.2371708245126284*(alpha[5]+alpha[4])+0.1767766952966368*(alpha[3]+alpha[2])-0.1767766952966368*alpha[1]+0.1767766952966368*alpha[0] > 0) { 
    fUpwindQuad[35] = hyb_3x3v_p1_surfx4_eval_quad_node_35_r(fEdge); 
  } else { 
    fUpwindQuad[35] = hyb_3x3v_p1_surfx4_eval_quad_node_35_l(fSkin); 
  } 
  if ((-0.2371708245126284*alpha[22])+0.2371708245126284*alpha[20]-0.2371708245126284*alpha[19]+0.2371708245126284*(alpha[18]+alpha[14]+alpha[13])-0.2371708245126284*alpha[12]+0.2371708245126284*(alpha[11]+alpha[10])-0.2371708245126284*alpha[9]+0.1767766952966368*alpha[8]-0.1767766952966368*(alpha[7]+alpha[6])-0.2371708245126284*(alpha[5]+alpha[4])-0.1767766952966368*(alpha[3]+alpha[2])+0.1767766952966368*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[36] = hyb_3x3v_p1_surfx4_eval_quad_node_36_r(fEdge); 
  } else { 
    fUpwindQuad[36] = hyb_3x3v_p1_surfx4_eval_quad_node_36_l(fSkin); 
  } 
  if ((-0.2371708245126284*alpha[19])+0.2371708245126284*(alpha[18]+alpha[11]+alpha[10])-0.2371708245126284*alpha[9]+0.1767766952966368*alpha[8]-0.1767766952966368*(alpha[7]+alpha[6])-0.2371708245126284*alpha[4]-0.1767766952966368*(alpha[3]+alpha[2])+0.1767766952966368*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[37] = hyb_3x3v_p1_surfx4_eval_quad_node_37_r(fEdge); 
  } else { 
    fUpwindQuad[37] = hyb_3x3v_p1_surfx4_eval_quad_node_37_l(fSkin); 
  } 
  if (0.2371708245126284*alpha[22]-0.2371708245126284*(alpha[20]+alpha[19])+0.2371708245126284*alpha[18]-0.2371708245126284*(alpha[14]+alpha[13])+0.2371708245126284*(alpha[12]+alpha[11]+alpha[10])-0.2371708245126284*alpha[9]+0.1767766952966368*alpha[8]-0.1767766952966368*(alpha[7]+alpha[6])+0.2371708245126284*alpha[5]-0.2371708245126284*alpha[4]-0.1767766952966368*(alpha[3]+alpha[2])+0.1767766952966368*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[38] = hyb_3x3v_p1_surfx4_eval_quad_node_38_r(fEdge); 
  } else { 
    fUpwindQuad[38] = hyb_3x3v_p1_surfx4_eval_quad_node_38_l(fSkin); 
  } 
  if ((-0.2371708245126284*alpha[22])+0.2371708245126284*(alpha[20]+alpha[14]+alpha[13])-0.2371708245126284*alpha[12]+0.1767766952966368*alpha[8]-0.1767766952966368*(alpha[7]+alpha[6])-0.2371708245126284*alpha[5]-0.1767766952966368*(alpha[3]+alpha[2])+0.1767766952966368*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[39] = hyb_3x3v_p1_surfx4_eval_quad_node_39_r(fEdge); 
  } else { 
    fUpwindQuad[39] = hyb_3x3v_p1_surfx4_eval_quad_node_39_l(fSkin); 
  } 
  if (0.1767766952966368*alpha[8]-0.1767766952966368*(alpha[7]+alpha[6]+alpha[3]+alpha[2])+0.1767766952966368*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[40] = hyb_3x3v_p1_surfx4_eval_quad_node_40_r(fEdge); 
  } else { 
    fUpwindQuad[40] = hyb_3x3v_p1_surfx4_eval_quad_node_40_l(fSkin); 
  } 
  if (0.2371708245126284*alpha[22]-0.2371708245126284*(alpha[20]+alpha[14]+alpha[13])+0.2371708245126284*alpha[12]+0.1767766952966368*alpha[8]-0.1767766952966368*(alpha[7]+alpha[6])+0.2371708245126284*alpha[5]-0.1767766952966368*(alpha[3]+alpha[2])+0.1767766952966368*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[41] = hyb_3x3v_p1_surfx4_eval_quad_node_41_r(fEdge); 
  } else { 
    fUpwindQuad[41] = hyb_3x3v_p1_surfx4_eval_quad_node_41_l(fSkin); 
  } 
  if ((-0.2371708245126284*alpha[22])+0.2371708245126284*(alpha[20]+alpha[19])-0.2371708245126284*alpha[18]+0.2371708245126284*(alpha[14]+alpha[13])-0.2371708245126284*(alpha[12]+alpha[11]+alpha[10])+0.2371708245126284*alpha[9]+0.1767766952966368*alpha[8]-0.1767766952966368*(alpha[7]+alpha[6])-0.2371708245126284*alpha[5]+0.2371708245126284*alpha[4]-0.1767766952966368*(alpha[3]+alpha[2])+0.1767766952966368*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[42] = hyb_3x3v_p1_surfx4_eval_quad_node_42_r(fEdge); 
  } else { 
    fUpwindQuad[42] = hyb_3x3v_p1_surfx4_eval_quad_node_42_l(fSkin); 
  } 
  if (0.2371708245126284*alpha[19]-0.2371708245126284*(alpha[18]+alpha[11]+alpha[10])+0.2371708245126284*alpha[9]+0.1767766952966368*alpha[8]-0.1767766952966368*(alpha[7]+alpha[6])+0.2371708245126284*alpha[4]-0.1767766952966368*(alpha[3]+alpha[2])+0.1767766952966368*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[43] = hyb_3x3v_p1_surfx4_eval_quad_node_43_r(fEdge); 
  } else { 
    fUpwindQuad[43] = hyb_3x3v_p1_surfx4_eval_quad_node_43_l(fSkin); 
  } 
  if (0.2371708245126284*alpha[22]-0.2371708245126284*alpha[20]+0.2371708245126284*alpha[19]-0.2371708245126284*(alpha[18]+alpha[14]+alpha[13])+0.2371708245126284*alpha[12]-0.2371708245126284*(alpha[11]+alpha[10])+0.2371708245126284*alpha[9]+0.1767766952966368*alpha[8]-0.1767766952966368*(alpha[7]+alpha[6])+0.2371708245126284*(alpha[5]+alpha[4])-0.1767766952966368*(alpha[3]+alpha[2])+0.1767766952966368*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[44] = hyb_3x3v_p1_surfx4_eval_quad_node_44_r(fEdge); 
  } else { 
    fUpwindQuad[44] = hyb_3x3v_p1_surfx4_eval_quad_node_44_l(fSkin); 
  } 
  if (0.2371708245126284*(alpha[22]+alpha[20]+alpha[19])-0.2371708245126284*(alpha[18]+alpha[14])+0.2371708245126284*alpha[13]-0.2371708245126284*(alpha[12]+alpha[11])+0.2371708245126284*alpha[10]-0.2371708245126284*alpha[9]-0.1767766952966368*alpha[8]+0.1767766952966368*alpha[7]-0.1767766952966368*alpha[6]-0.2371708245126284*(alpha[5]+alpha[4])+0.1767766952966368*alpha[3]-0.1767766952966368*alpha[2]+0.1767766952966368*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[45] = hyb_3x3v_p1_surfx4_eval_quad_node_45_r(fEdge); 
  } else { 
    fUpwindQuad[45] = hyb_3x3v_p1_surfx4_eval_quad_node_45_l(fSkin); 
  } 
  if (0.2371708245126284*alpha[19]-0.2371708245126284*(alpha[18]+alpha[11])+0.2371708245126284*alpha[10]-0.2371708245126284*alpha[9]-0.1767766952966368*alpha[8]+0.1767766952966368*alpha[7]-0.1767766952966368*alpha[6]-0.2371708245126284*alpha[4]+0.1767766952966368*alpha[3]-0.1767766952966368*alpha[2]+0.1767766952966368*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[46] = hyb_3x3v_p1_surfx4_eval_quad_node_46_r(fEdge); 
  } else { 
    fUpwindQuad[46] = hyb_3x3v_p1_surfx4_eval_quad_node_46_l(fSkin); 
  } 
  if ((-0.2371708245126284*(alpha[22]+alpha[20]))+0.2371708245126284*alpha[19]-0.2371708245126284*alpha[18]+0.2371708245126284*alpha[14]-0.2371708245126284*alpha[13]+0.2371708245126284*alpha[12]-0.2371708245126284*alpha[11]+0.2371708245126284*alpha[10]-0.2371708245126284*alpha[9]-0.1767766952966368*alpha[8]+0.1767766952966368*alpha[7]-0.1767766952966368*alpha[6]+0.2371708245126284*alpha[5]-0.2371708245126284*alpha[4]+0.1767766952966368*alpha[3]-0.1767766952966368*alpha[2]+0.1767766952966368*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[47] = hyb_3x3v_p1_surfx4_eval_quad_node_47_r(fEdge); 
  } else { 
    fUpwindQuad[47] = hyb_3x3v_p1_surfx4_eval_quad_node_47_l(fSkin); 
  } 
  if (0.2371708245126284*(alpha[22]+alpha[20])-0.2371708245126284*alpha[14]+0.2371708245126284*alpha[13]-0.2371708245126284*alpha[12]-0.1767766952966368*alpha[8]+0.1767766952966368*alpha[7]-0.1767766952966368*alpha[6]-0.2371708245126284*alpha[5]+0.1767766952966368*alpha[3]-0.1767766952966368*alpha[2]+0.1767766952966368*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[48] = hyb_3x3v_p1_surfx4_eval_quad_node_48_r(fEdge); 
  } else { 
    fUpwindQuad[48] = hyb_3x3v_p1_surfx4_eval_quad_node_48_l(fSkin); 
  } 
  if ((-0.1767766952966368*alpha[8])+0.1767766952966368*alpha[7]-0.1767766952966368*alpha[6]+0.1767766952966368*alpha[3]-0.1767766952966368*alpha[2]+0.1767766952966368*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[49] = hyb_3x3v_p1_surfx4_eval_quad_node_49_r(fEdge); 
  } else { 
    fUpwindQuad[49] = hyb_3x3v_p1_surfx4_eval_quad_node_49_l(fSkin); 
  } 
  if ((-0.2371708245126284*(alpha[22]+alpha[20]))+0.2371708245126284*alpha[14]-0.2371708245126284*alpha[13]+0.2371708245126284*alpha[12]-0.1767766952966368*alpha[8]+0.1767766952966368*alpha[7]-0.1767766952966368*alpha[6]+0.2371708245126284*alpha[5]+0.1767766952966368*alpha[3]-0.1767766952966368*alpha[2]+0.1767766952966368*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[50] = hyb_3x3v_p1_surfx4_eval_quad_node_50_r(fEdge); 
  } else { 
    fUpwindQuad[50] = hyb_3x3v_p1_surfx4_eval_quad_node_50_l(fSkin); 
  } 
  if (0.2371708245126284*(alpha[22]+alpha[20])-0.2371708245126284*alpha[19]+0.2371708245126284*alpha[18]-0.2371708245126284*alpha[14]+0.2371708245126284*alpha[13]-0.2371708245126284*alpha[12]+0.2371708245126284*alpha[11]-0.2371708245126284*alpha[10]+0.2371708245126284*alpha[9]-0.1767766952966368*alpha[8]+0.1767766952966368*alpha[7]-0.1767766952966368*alpha[6]-0.2371708245126284*alpha[5]+0.2371708245126284*alpha[4]+0.1767766952966368*alpha[3]-0.1767766952966368*alpha[2]+0.1767766952966368*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[51] = hyb_3x3v_p1_surfx4_eval_quad_node_51_r(fEdge); 
  } else { 
    fUpwindQuad[51] = hyb_3x3v_p1_surfx4_eval_quad_node_51_l(fSkin); 
  } 
  if ((-0.2371708245126284*alpha[19])+0.2371708245126284*(alpha[18]+alpha[11])-0.2371708245126284*alpha[10]+0.2371708245126284*alpha[9]-0.1767766952966368*alpha[8]+0.1767766952966368*alpha[7]-0.1767766952966368*alpha[6]+0.2371708245126284*alpha[4]+0.1767766952966368*alpha[3]-0.1767766952966368*alpha[2]+0.1767766952966368*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[52] = hyb_3x3v_p1_surfx4_eval_quad_node_52_r(fEdge); 
  } else { 
    fUpwindQuad[52] = hyb_3x3v_p1_surfx4_eval_quad_node_52_l(fSkin); 
  } 
  if ((-0.2371708245126284*(alpha[22]+alpha[20]+alpha[19]))+0.2371708245126284*(alpha[18]+alpha[14])-0.2371708245126284*alpha[13]+0.2371708245126284*(alpha[12]+alpha[11])-0.2371708245126284*alpha[10]+0.2371708245126284*alpha[9]-0.1767766952966368*alpha[8]+0.1767766952966368*alpha[7]-0.1767766952966368*alpha[6]+0.2371708245126284*(alpha[5]+alpha[4])+0.1767766952966368*alpha[3]-0.1767766952966368*alpha[2]+0.1767766952966368*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[53] = hyb_3x3v_p1_surfx4_eval_quad_node_53_r(fEdge); 
  } else { 
    fUpwindQuad[53] = hyb_3x3v_p1_surfx4_eval_quad_node_53_l(fSkin); 
  } 
  if (0.2371708245126284*alpha[22]-0.2371708245126284*alpha[20]+0.2371708245126284*(alpha[19]+alpha[18]+alpha[14])-0.2371708245126284*(alpha[13]+alpha[12])+0.2371708245126284*alpha[11]-0.2371708245126284*(alpha[10]+alpha[9])-0.1767766952966368*(alpha[8]+alpha[7])+0.1767766952966368*alpha[6]-0.2371708245126284*(alpha[5]+alpha[4])-0.1767766952966368*alpha[3]+0.1767766952966368*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[54] = hyb_3x3v_p1_surfx4_eval_quad_node_54_r(fEdge); 
  } else { 
    fUpwindQuad[54] = hyb_3x3v_p1_surfx4_eval_quad_node_54_l(fSkin); 
  } 
  if (0.2371708245126284*(alpha[19]+alpha[18]+alpha[11])-0.2371708245126284*(alpha[10]+alpha[9])-0.1767766952966368*(alpha[8]+alpha[7])+0.1767766952966368*alpha[6]-0.2371708245126284*alpha[4]-0.1767766952966368*alpha[3]+0.1767766952966368*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[55] = hyb_3x3v_p1_surfx4_eval_quad_node_55_r(fEdge); 
  } else { 
    fUpwindQuad[55] = hyb_3x3v_p1_surfx4_eval_quad_node_55_l(fSkin); 
  } 
  if ((-0.2371708245126284*alpha[22])+0.2371708245126284*(alpha[20]+alpha[19]+alpha[18])-0.2371708245126284*alpha[14]+0.2371708245126284*(alpha[13]+alpha[12]+alpha[11])-0.2371708245126284*(alpha[10]+alpha[9])-0.1767766952966368*(alpha[8]+alpha[7])+0.1767766952966368*alpha[6]+0.2371708245126284*alpha[5]-0.2371708245126284*alpha[4]-0.1767766952966368*alpha[3]+0.1767766952966368*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[56] = hyb_3x3v_p1_surfx4_eval_quad_node_56_r(fEdge); 
  } else { 
    fUpwindQuad[56] = hyb_3x3v_p1_surfx4_eval_quad_node_56_l(fSkin); 
  } 
  if (0.2371708245126284*alpha[22]-0.2371708245126284*alpha[20]+0.2371708245126284*alpha[14]-0.2371708245126284*(alpha[13]+alpha[12])-0.1767766952966368*(alpha[8]+alpha[7])+0.1767766952966368*alpha[6]-0.2371708245126284*alpha[5]-0.1767766952966368*alpha[3]+0.1767766952966368*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[57] = hyb_3x3v_p1_surfx4_eval_quad_node_57_r(fEdge); 
  } else { 
    fUpwindQuad[57] = hyb_3x3v_p1_surfx4_eval_quad_node_57_l(fSkin); 
  } 
  if ((-0.1767766952966368*(alpha[8]+alpha[7]))+0.1767766952966368*alpha[6]-0.1767766952966368*alpha[3]+0.1767766952966368*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[58] = hyb_3x3v_p1_surfx4_eval_quad_node_58_r(fEdge); 
  } else { 
    fUpwindQuad[58] = hyb_3x3v_p1_surfx4_eval_quad_node_58_l(fSkin); 
  } 
  if ((-0.2371708245126284*alpha[22])+0.2371708245126284*alpha[20]-0.2371708245126284*alpha[14]+0.2371708245126284*(alpha[13]+alpha[12])-0.1767766952966368*(alpha[8]+alpha[7])+0.1767766952966368*alpha[6]+0.2371708245126284*alpha[5]-0.1767766952966368*alpha[3]+0.1767766952966368*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[59] = hyb_3x3v_p1_surfx4_eval_quad_node_59_r(fEdge); 
  } else { 
    fUpwindQuad[59] = hyb_3x3v_p1_surfx4_eval_quad_node_59_l(fSkin); 
  } 
  if (0.2371708245126284*alpha[22]-0.2371708245126284*(alpha[20]+alpha[19]+alpha[18])+0.2371708245126284*alpha[14]-0.2371708245126284*(alpha[13]+alpha[12]+alpha[11])+0.2371708245126284*(alpha[10]+alpha[9])-0.1767766952966368*(alpha[8]+alpha[7])+0.1767766952966368*alpha[6]-0.2371708245126284*alpha[5]+0.2371708245126284*alpha[4]-0.1767766952966368*alpha[3]+0.1767766952966368*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[60] = hyb_3x3v_p1_surfx4_eval_quad_node_60_r(fEdge); 
  } else { 
    fUpwindQuad[60] = hyb_3x3v_p1_surfx4_eval_quad_node_60_l(fSkin); 
  } 
  if ((-0.2371708245126284*(alpha[19]+alpha[18]+alpha[11]))+0.2371708245126284*(alpha[10]+alpha[9])-0.1767766952966368*(alpha[8]+alpha[7])+0.1767766952966368*alpha[6]+0.2371708245126284*alpha[4]-0.1767766952966368*alpha[3]+0.1767766952966368*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[61] = hyb_3x3v_p1_surfx4_eval_quad_node_61_r(fEdge); 
  } else { 
    fUpwindQuad[61] = hyb_3x3v_p1_surfx4_eval_quad_node_61_l(fSkin); 
  } 
  if ((-0.2371708245126284*alpha[22])+0.2371708245126284*alpha[20]-0.2371708245126284*(alpha[19]+alpha[18]+alpha[14])+0.2371708245126284*(alpha[13]+alpha[12])-0.2371708245126284*alpha[11]+0.2371708245126284*(alpha[10]+alpha[9])-0.1767766952966368*(alpha[8]+alpha[7])+0.1767766952966368*alpha[6]+0.2371708245126284*(alpha[5]+alpha[4])-0.1767766952966368*alpha[3]+0.1767766952966368*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[62] = hyb_3x3v_p1_surfx4_eval_quad_node_62_r(fEdge); 
  } else { 
    fUpwindQuad[62] = hyb_3x3v_p1_surfx4_eval_quad_node_62_l(fSkin); 
  } 
  if ((-0.2371708245126284*(alpha[22]+alpha[20]+alpha[19]+alpha[18]+alpha[14]+alpha[13]+alpha[12]+alpha[11]+alpha[10]+alpha[9]))+0.1767766952966368*(alpha[8]+alpha[7]+alpha[6])-0.2371708245126284*(alpha[5]+alpha[4])+0.1767766952966368*(alpha[3]+alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[63] = hyb_3x3v_p1_surfx4_eval_quad_node_63_r(fEdge); 
  } else { 
    fUpwindQuad[63] = hyb_3x3v_p1_surfx4_eval_quad_node_63_l(fSkin); 
  } 
  if ((-0.2371708245126284*(alpha[19]+alpha[18]+alpha[11]+alpha[10]+alpha[9]))+0.1767766952966368*(alpha[8]+alpha[7]+alpha[6])-0.2371708245126284*alpha[4]+0.1767766952966368*(alpha[3]+alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[64] = hyb_3x3v_p1_surfx4_eval_quad_node_64_r(fEdge); 
  } else { 
    fUpwindQuad[64] = hyb_3x3v_p1_surfx4_eval_quad_node_64_l(fSkin); 
  } 
  if (0.2371708245126284*(alpha[22]+alpha[20])-0.2371708245126284*(alpha[19]+alpha[18])+0.2371708245126284*(alpha[14]+alpha[13]+alpha[12])-0.2371708245126284*(alpha[11]+alpha[10]+alpha[9])+0.1767766952966368*(alpha[8]+alpha[7]+alpha[6])+0.2371708245126284*alpha[5]-0.2371708245126284*alpha[4]+0.1767766952966368*(alpha[3]+alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[65] = hyb_3x3v_p1_surfx4_eval_quad_node_65_r(fEdge); 
  } else { 
    fUpwindQuad[65] = hyb_3x3v_p1_surfx4_eval_quad_node_65_l(fSkin); 
  } 
  if ((-0.2371708245126284*(alpha[22]+alpha[20]+alpha[14]+alpha[13]+alpha[12]))+0.1767766952966368*(alpha[8]+alpha[7]+alpha[6])-0.2371708245126284*alpha[5]+0.1767766952966368*(alpha[3]+alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[66] = hyb_3x3v_p1_surfx4_eval_quad_node_66_r(fEdge); 
  } else { 
    fUpwindQuad[66] = hyb_3x3v_p1_surfx4_eval_quad_node_66_l(fSkin); 
  } 
  if (0.1767766952966368*(alpha[8]+alpha[7]+alpha[6]+alpha[3]+alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[67] = hyb_3x3v_p1_surfx4_eval_quad_node_67_r(fEdge); 
  } else { 
    fUpwindQuad[67] = hyb_3x3v_p1_surfx4_eval_quad_node_67_l(fSkin); 
  } 
  if (0.2371708245126284*(alpha[22]+alpha[20]+alpha[14]+alpha[13]+alpha[12])+0.1767766952966368*(alpha[8]+alpha[7]+alpha[6])+0.2371708245126284*alpha[5]+0.1767766952966368*(alpha[3]+alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[68] = hyb_3x3v_p1_surfx4_eval_quad_node_68_r(fEdge); 
  } else { 
    fUpwindQuad[68] = hyb_3x3v_p1_surfx4_eval_quad_node_68_l(fSkin); 
  } 
  if ((-0.2371708245126284*(alpha[22]+alpha[20]))+0.2371708245126284*(alpha[19]+alpha[18])-0.2371708245126284*(alpha[14]+alpha[13]+alpha[12])+0.2371708245126284*(alpha[11]+alpha[10]+alpha[9])+0.1767766952966368*(alpha[8]+alpha[7]+alpha[6])-0.2371708245126284*alpha[5]+0.2371708245126284*alpha[4]+0.1767766952966368*(alpha[3]+alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[69] = hyb_3x3v_p1_surfx4_eval_quad_node_69_r(fEdge); 
  } else { 
    fUpwindQuad[69] = hyb_3x3v_p1_surfx4_eval_quad_node_69_l(fSkin); 
  } 
  if (0.2371708245126284*(alpha[19]+alpha[18]+alpha[11]+alpha[10]+alpha[9])+0.1767766952966368*(alpha[8]+alpha[7]+alpha[6])+0.2371708245126284*alpha[4]+0.1767766952966368*(alpha[3]+alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[70] = hyb_3x3v_p1_surfx4_eval_quad_node_70_r(fEdge); 
  } else { 
    fUpwindQuad[70] = hyb_3x3v_p1_surfx4_eval_quad_node_70_l(fSkin); 
  } 
  if (0.2371708245126284*(alpha[22]+alpha[20]+alpha[19]+alpha[18]+alpha[14]+alpha[13]+alpha[12]+alpha[11]+alpha[10]+alpha[9])+0.1767766952966368*(alpha[8]+alpha[7]+alpha[6])+0.2371708245126284*(alpha[5]+alpha[4])+0.1767766952966368*(alpha[3]+alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[71] = hyb_3x3v_p1_surfx4_eval_quad_node_71_r(fEdge); 
  } else { 
    fUpwindQuad[71] = hyb_3x3v_p1_surfx4_eval_quad_node_71_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_3x3v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.1767766952966368*alpha[22]*fUpwind[22]+0.1767766952966368*alpha[20]*fUpwind[20]+0.1767766952966368*alpha[19]*fUpwind[19]+0.1767766952966368*alpha[18]*fUpwind[18]+0.1767766952966368*alpha[14]*fUpwind[14]+0.1767766952966368*alpha[13]*fUpwind[13]+0.1767766952966368*alpha[12]*fUpwind[12]+0.1767766952966368*alpha[11]*fUpwind[11]+0.1767766952966368*alpha[10]*fUpwind[10]+0.1767766952966368*alpha[9]*fUpwind[9]+0.1767766952966368*alpha[8]*fUpwind[8]+0.1767766952966368*alpha[7]*fUpwind[7]+0.1767766952966368*alpha[6]*fUpwind[6]+0.1767766952966368*alpha[5]*fUpwind[5]+0.1767766952966368*alpha[4]*fUpwind[4]+0.1767766952966368*alpha[3]*fUpwind[3]+0.1767766952966368*alpha[2]*fUpwind[2]+0.1767766952966368*alpha[1]*fUpwind[1]+0.1767766952966368*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.1767766952966368*alpha[22]*fUpwind[27]+0.1767766952966368*alpha[19]*fUpwind[26]+0.1767766952966368*alpha[14]*fUpwind[21]+0.1767766952966368*alpha[13]*fUpwind[20]+0.1767766952966368*fUpwind[13]*alpha[20]+0.1767766952966368*alpha[11]*fUpwind[18]+0.1767766952966368*fUpwind[11]*alpha[18]+0.1767766952966368*alpha[10]*fUpwind[17]+0.1767766952966368*alpha[8]*fUpwind[16]+0.1767766952966368*alpha[5]*fUpwind[12]+0.1767766952966368*fUpwind[5]*alpha[12]+0.1767766952966368*alpha[4]*fUpwind[9]+0.1767766952966368*fUpwind[4]*alpha[9]+0.1767766952966368*alpha[3]*fUpwind[7]+0.1767766952966368*fUpwind[3]*alpha[7]+0.1767766952966368*alpha[2]*fUpwind[6]+0.1767766952966368*fUpwind[2]*alpha[6]+0.1767766952966368*alpha[0]*fUpwind[1]+0.1767766952966368*fUpwind[0]*alpha[1]; 
  Ghat[2] = 0.1767766952966368*alpha[18]*fUpwind[26]+0.1767766952966368*alpha[14]*fUpwind[22]+0.1767766952966368*fUpwind[14]*alpha[22]+0.1767766952966368*alpha[12]*fUpwind[20]+0.1767766952966368*fUpwind[12]*alpha[20]+0.1767766952966368*alpha[11]*fUpwind[19]+0.1767766952966368*fUpwind[11]*alpha[19]+0.1767766952966368*alpha[9]*fUpwind[17]+0.1767766952966368*alpha[7]*fUpwind[16]+0.1767766952966368*alpha[5]*fUpwind[13]+0.1767766952966368*fUpwind[5]*alpha[13]+0.1767766952966368*alpha[4]*fUpwind[10]+0.1767766952966368*fUpwind[4]*alpha[10]+0.1767766952966368*alpha[3]*fUpwind[8]+0.1767766952966368*fUpwind[3]*alpha[8]+0.1767766952966368*alpha[1]*fUpwind[6]+0.1767766952966368*fUpwind[1]*alpha[6]+0.1767766952966368*alpha[0]*fUpwind[2]+0.1767766952966368*fUpwind[0]*alpha[2]; 
  Ghat[3] = 0.1767766952966368*alpha[20]*fUpwind[27]+0.1767766952966368*alpha[13]*fUpwind[22]+0.1767766952966368*fUpwind[13]*alpha[22]+0.1767766952966368*alpha[12]*fUpwind[21]+0.1767766952966368*alpha[10]*fUpwind[19]+0.1767766952966368*fUpwind[10]*alpha[19]+0.1767766952966368*alpha[9]*fUpwind[18]+0.1767766952966368*fUpwind[9]*alpha[18]+0.1767766952966368*alpha[6]*fUpwind[16]+0.1767766952966368*alpha[5]*fUpwind[14]+0.1767766952966368*fUpwind[5]*alpha[14]+0.1767766952966368*alpha[4]*fUpwind[11]+0.1767766952966368*fUpwind[4]*alpha[11]+0.1767766952966368*alpha[2]*fUpwind[8]+0.1767766952966368*fUpwind[2]*alpha[8]+0.1767766952966368*alpha[1]*fUpwind[7]+0.1767766952966368*fUpwind[1]*alpha[7]+0.1767766952966368*alpha[0]*fUpwind[3]+0.1767766952966368*fUpwind[0]*alpha[3]; 
  Ghat[4] = 0.1581138830084189*alpha[19]*fUpwind[39]+0.1581138830084189*alpha[18]*fUpwind[38]+0.1581138830084189*alpha[11]*fUpwind[35]+0.1581138830084189*alpha[10]*fUpwind[34]+0.1581138830084189*alpha[9]*fUpwind[33]+0.1581138830084189*alpha[4]*fUpwind[32]+0.1767766952966368*alpha[22]*fUpwind[30]+0.1767766952966368*alpha[20]*fUpwind[28]+0.1767766952966368*alpha[14]*fUpwind[25]+0.1767766952966368*alpha[13]*fUpwind[24]+0.1767766952966368*alpha[12]*fUpwind[23]+0.1767766952966368*alpha[8]*fUpwind[19]+0.1767766952966368*fUpwind[8]*alpha[19]+0.1767766952966368*alpha[7]*fUpwind[18]+0.1767766952966368*fUpwind[7]*alpha[18]+0.1767766952966368*alpha[6]*fUpwind[17]+0.1767766952966368*alpha[5]*fUpwind[15]+0.1767766952966368*alpha[3]*fUpwind[11]+0.1767766952966368*fUpwind[3]*alpha[11]+0.1767766952966368*alpha[2]*fUpwind[10]+0.1767766952966368*fUpwind[2]*alpha[10]+0.1767766952966368*alpha[1]*fUpwind[9]+0.1767766952966368*fUpwind[1]*alpha[9]+0.1767766952966368*alpha[0]*fUpwind[4]+0.1767766952966368*fUpwind[0]*alpha[4]; 
  Ghat[5] = 0.1581138830084189*alpha[22]*fUpwind[55]+0.1581138830084189*alpha[20]*fUpwind[53]+0.1581138830084189*alpha[14]*fUpwind[51]+0.1581138830084189*alpha[13]*fUpwind[50]+0.1581138830084189*alpha[12]*fUpwind[49]+0.1581138830084189*alpha[5]*fUpwind[48]+0.1767766952966368*alpha[19]*fUpwind[30]+0.1767766952966368*alpha[18]*fUpwind[29]+0.1767766952966368*alpha[11]*fUpwind[25]+0.1767766952966368*alpha[10]*fUpwind[24]+0.1767766952966368*alpha[9]*fUpwind[23]+0.1767766952966368*alpha[8]*fUpwind[22]+0.1767766952966368*fUpwind[8]*alpha[22]+0.1767766952966368*alpha[7]*fUpwind[21]+0.1767766952966368*alpha[6]*fUpwind[20]+0.1767766952966368*fUpwind[6]*alpha[20]+0.1767766952966368*alpha[4]*fUpwind[15]+0.1767766952966368*alpha[3]*fUpwind[14]+0.1767766952966368*fUpwind[3]*alpha[14]+0.1767766952966368*alpha[2]*fUpwind[13]+0.1767766952966368*fUpwind[2]*alpha[13]+0.1767766952966368*alpha[1]*fUpwind[12]+0.1767766952966368*fUpwind[1]*alpha[12]+0.1767766952966368*alpha[0]*fUpwind[5]+0.1767766952966368*fUpwind[0]*alpha[5]; 
  Ghat[6] = 0.1767766952966368*alpha[14]*fUpwind[27]+0.1767766952966368*alpha[11]*fUpwind[26]+0.1767766952966368*fUpwind[21]*alpha[22]+0.1767766952966368*alpha[5]*fUpwind[20]+0.1767766952966368*fUpwind[5]*alpha[20]+0.1767766952966368*alpha[18]*fUpwind[19]+0.1767766952966368*fUpwind[18]*alpha[19]+0.1767766952966368*alpha[4]*fUpwind[17]+0.1767766952966368*alpha[3]*fUpwind[16]+0.1767766952966368*alpha[12]*fUpwind[13]+0.1767766952966368*fUpwind[12]*alpha[13]+0.1767766952966368*alpha[9]*fUpwind[10]+0.1767766952966368*fUpwind[9]*alpha[10]+0.1767766952966368*alpha[7]*fUpwind[8]+0.1767766952966368*fUpwind[7]*alpha[8]+0.1767766952966368*alpha[0]*fUpwind[6]+0.1767766952966368*fUpwind[0]*alpha[6]+0.1767766952966368*alpha[1]*fUpwind[2]+0.1767766952966368*fUpwind[1]*alpha[2]; 
  Ghat[7] = 0.1767766952966368*alpha[13]*fUpwind[27]+0.1767766952966368*alpha[10]*fUpwind[26]+0.1767766952966368*alpha[20]*fUpwind[22]+0.1767766952966368*fUpwind[20]*alpha[22]+0.1767766952966368*alpha[5]*fUpwind[21]+0.1767766952966368*fUpwind[17]*alpha[19]+0.1767766952966368*alpha[4]*fUpwind[18]+0.1767766952966368*fUpwind[4]*alpha[18]+0.1767766952966368*alpha[2]*fUpwind[16]+0.1767766952966368*alpha[12]*fUpwind[14]+0.1767766952966368*fUpwind[12]*alpha[14]+0.1767766952966368*alpha[9]*fUpwind[11]+0.1767766952966368*fUpwind[9]*alpha[11]+0.1767766952966368*alpha[6]*fUpwind[8]+0.1767766952966368*fUpwind[6]*alpha[8]+0.1767766952966368*alpha[0]*fUpwind[7]+0.1767766952966368*fUpwind[0]*alpha[7]+0.1767766952966368*alpha[1]*fUpwind[3]+0.1767766952966368*fUpwind[1]*alpha[3]; 
  Ghat[8] = 0.1767766952966368*alpha[12]*fUpwind[27]+0.1767766952966368*alpha[9]*fUpwind[26]+0.1767766952966368*alpha[5]*fUpwind[22]+0.1767766952966368*fUpwind[5]*alpha[22]+0.1767766952966368*alpha[20]*fUpwind[21]+0.1767766952966368*alpha[4]*fUpwind[19]+0.1767766952966368*fUpwind[4]*alpha[19]+0.1767766952966368*fUpwind[17]*alpha[18]+0.1767766952966368*alpha[1]*fUpwind[16]+0.1767766952966368*alpha[13]*fUpwind[14]+0.1767766952966368*fUpwind[13]*alpha[14]+0.1767766952966368*alpha[10]*fUpwind[11]+0.1767766952966368*fUpwind[10]*alpha[11]+0.1767766952966368*alpha[0]*fUpwind[8]+0.1767766952966368*fUpwind[0]*alpha[8]+0.1767766952966368*alpha[6]*fUpwind[7]+0.1767766952966368*fUpwind[6]*alpha[7]+0.1767766952966368*alpha[2]*fUpwind[3]+0.1767766952966368*fUpwind[2]*alpha[3]; 
  Ghat[9] = 0.1581138830084189*alpha[19]*fUpwind[43]+0.1581138830084189*alpha[11]*fUpwind[38]+0.1581138830084189*alpha[10]*fUpwind[37]+0.1581138830084189*alpha[18]*fUpwind[35]+0.1581138830084189*alpha[4]*fUpwind[33]+0.1581138830084189*alpha[9]*fUpwind[32]+0.1767766952966368*alpha[22]*fUpwind[31]+0.1767766952966368*alpha[14]*fUpwind[29]+0.1767766952966368*alpha[13]*fUpwind[28]+0.1767766952966368*alpha[8]*fUpwind[26]+0.1767766952966368*alpha[20]*fUpwind[24]+0.1767766952966368*alpha[5]*fUpwind[23]+0.1767766952966368*fUpwind[16]*alpha[19]+0.1767766952966368*alpha[3]*fUpwind[18]+0.1767766952966368*fUpwind[3]*alpha[18]+0.1767766952966368*alpha[2]*fUpwind[17]+0.1767766952966368*alpha[12]*fUpwind[15]+0.1767766952966368*alpha[7]*fUpwind[11]+0.1767766952966368*fUpwind[7]*alpha[11]+0.1767766952966368*alpha[6]*fUpwind[10]+0.1767766952966368*fUpwind[6]*alpha[10]+0.1767766952966368*alpha[0]*fUpwind[9]+0.1767766952966368*fUpwind[0]*alpha[9]+0.1767766952966368*alpha[1]*fUpwind[4]+0.1767766952966368*fUpwind[1]*alpha[4]; 
  Ghat[10] = 0.1581138830084189*alpha[18]*fUpwind[43]+0.1581138830084189*alpha[11]*fUpwind[39]+0.1581138830084189*alpha[9]*fUpwind[37]+0.1581138830084189*alpha[19]*fUpwind[35]+0.1581138830084189*alpha[4]*fUpwind[34]+0.1581138830084189*alpha[10]*fUpwind[32]+0.1767766952966368*alpha[14]*fUpwind[30]+0.1767766952966368*alpha[12]*fUpwind[28]+0.1767766952966368*alpha[7]*fUpwind[26]+0.1767766952966368*alpha[22]*fUpwind[25]+0.1767766952966368*alpha[5]*fUpwind[24]+0.1767766952966368*alpha[20]*fUpwind[23]+0.1767766952966368*alpha[3]*fUpwind[19]+0.1767766952966368*fUpwind[3]*alpha[19]+0.1767766952966368*fUpwind[16]*alpha[18]+0.1767766952966368*alpha[1]*fUpwind[17]+0.1767766952966368*alpha[13]*fUpwind[15]+0.1767766952966368*alpha[8]*fUpwind[11]+0.1767766952966368*fUpwind[8]*alpha[11]+0.1767766952966368*alpha[0]*fUpwind[10]+0.1767766952966368*fUpwind[0]*alpha[10]+0.1767766952966368*alpha[6]*fUpwind[9]+0.1767766952966368*fUpwind[6]*alpha[9]+0.1767766952966368*alpha[2]*fUpwind[4]+0.1767766952966368*fUpwind[2]*alpha[4]; 
  Ghat[11] = 0.1581138830084189*alpha[10]*fUpwind[39]+0.1581138830084189*alpha[9]*fUpwind[38]+0.1581138830084189*alpha[4]*fUpwind[35]+0.1581138830084189*alpha[19]*fUpwind[34]+0.1581138830084189*alpha[18]*fUpwind[33]+0.1581138830084189*alpha[11]*fUpwind[32]+0.1767766952966368*alpha[20]*fUpwind[31]+0.1767766952966368*alpha[13]*fUpwind[30]+0.1767766952966368*alpha[12]*fUpwind[29]+0.1767766952966368*alpha[6]*fUpwind[26]+0.1767766952966368*alpha[5]*fUpwind[25]+0.1767766952966368*alpha[22]*fUpwind[24]+0.1767766952966368*alpha[2]*fUpwind[19]+0.1767766952966368*fUpwind[2]*alpha[19]+0.1767766952966368*alpha[1]*fUpwind[18]+0.1767766952966368*fUpwind[1]*alpha[18]+0.1767766952966368*alpha[14]*fUpwind[15]+0.1767766952966368*alpha[0]*fUpwind[11]+0.1767766952966368*fUpwind[0]*alpha[11]+0.1767766952966368*alpha[8]*fUpwind[10]+0.1767766952966368*fUpwind[8]*alpha[10]+0.1767766952966368*alpha[7]*fUpwind[9]+0.1767766952966368*fUpwind[7]*alpha[9]+0.1767766952966368*alpha[3]*fUpwind[4]+0.1767766952966368*fUpwind[3]*alpha[4]; 
  Ghat[12] = 0.1581138830084189*alpha[22]*fUpwind[59]+0.1581138830084189*alpha[14]*fUpwind[54]+0.1581138830084189*alpha[13]*fUpwind[53]+0.1581138830084189*alpha[20]*fUpwind[50]+0.1581138830084189*alpha[5]*fUpwind[49]+0.1581138830084189*alpha[12]*fUpwind[48]+0.1767766952966368*alpha[19]*fUpwind[31]+0.1767766952966368*alpha[11]*fUpwind[29]+0.1767766952966368*alpha[10]*fUpwind[28]+0.1767766952966368*alpha[8]*fUpwind[27]+0.1767766952966368*alpha[18]*fUpwind[25]+0.1767766952966368*alpha[4]*fUpwind[23]+0.1767766952966368*fUpwind[16]*alpha[22]+0.1767766952966368*alpha[3]*fUpwind[21]+0.1767766952966368*alpha[2]*fUpwind[20]+0.1767766952966368*fUpwind[2]*alpha[20]+0.1767766952966368*alpha[9]*fUpwind[15]+0.1767766952966368*alpha[7]*fUpwind[14]+0.1767766952966368*fUpwind[7]*alpha[14]+0.1767766952966368*alpha[6]*fUpwind[13]+0.1767766952966368*fUpwind[6]*alpha[13]+0.1767766952966368*alpha[0]*fUpwind[12]+0.1767766952966368*fUpwind[0]*alpha[12]+0.1767766952966368*alpha[1]*fUpwind[5]+0.1767766952966368*fUpwind[1]*alpha[5]; 
  Ghat[13] = 0.1581138830084189*alpha[14]*fUpwind[55]+0.1581138830084189*alpha[12]*fUpwind[53]+0.1581138830084189*alpha[22]*fUpwind[51]+0.1581138830084189*alpha[5]*fUpwind[50]+0.1581138830084189*alpha[20]*fUpwind[49]+0.1581138830084189*alpha[13]*fUpwind[48]+0.1767766952966368*alpha[18]*fUpwind[31]+0.1767766952966368*alpha[11]*fUpwind[30]+0.1767766952966368*alpha[9]*fUpwind[28]+0.1767766952966368*alpha[7]*fUpwind[27]+0.1767766952966368*alpha[19]*fUpwind[25]+0.1767766952966368*alpha[4]*fUpwind[24]+0.1767766952966368*alpha[3]*fUpwind[22]+0.1767766952966368*fUpwind[3]*alpha[22]+0.1767766952966368*alpha[1]*fUpwind[20]+0.1767766952966368*fUpwind[1]*alpha[20]+0.1767766952966368*alpha[10]*fUpwind[15]+0.1767766952966368*alpha[8]*fUpwind[14]+0.1767766952966368*fUpwind[8]*alpha[14]+0.1767766952966368*alpha[0]*fUpwind[13]+0.1767766952966368*fUpwind[0]*alpha[13]+0.1767766952966368*alpha[6]*fUpwind[12]+0.1767766952966368*fUpwind[6]*alpha[12]+0.1767766952966368*alpha[2]*fUpwind[5]+0.1767766952966368*fUpwind[2]*alpha[5]; 
  Ghat[14] = 0.1581138830084189*alpha[20]*fUpwind[59]+0.1581138830084189*alpha[13]*fUpwind[55]+0.1581138830084189*alpha[12]*fUpwind[54]+0.1581138830084189*alpha[5]*fUpwind[51]+0.1581138830084189*alpha[22]*fUpwind[50]+0.1581138830084189*alpha[14]*fUpwind[48]+0.1767766952966368*alpha[10]*fUpwind[30]+0.1767766952966368*alpha[9]*fUpwind[29]+0.1767766952966368*alpha[6]*fUpwind[27]+0.1767766952966368*alpha[4]*fUpwind[25]+0.1767766952966368*alpha[19]*fUpwind[24]+0.1767766952966368*alpha[18]*fUpwind[23]+0.1767766952966368*alpha[2]*fUpwind[22]+0.1767766952966368*fUpwind[2]*alpha[22]+0.1767766952966368*alpha[1]*fUpwind[21]+0.1767766952966368*fUpwind[16]*alpha[20]+0.1767766952966368*alpha[11]*fUpwind[15]+0.1767766952966368*alpha[0]*fUpwind[14]+0.1767766952966368*fUpwind[0]*alpha[14]+0.1767766952966368*alpha[8]*fUpwind[13]+0.1767766952966368*fUpwind[8]*alpha[13]+0.1767766952966368*alpha[7]*fUpwind[12]+0.1767766952966368*fUpwind[7]*alpha[12]+0.1767766952966368*alpha[3]*fUpwind[5]+0.1767766952966368*fUpwind[3]*alpha[5]; 
  Ghat[15] = 0.1581138830084189*alpha[22]*fUpwind[62]+0.1581138830084189*alpha[20]*fUpwind[60]+0.1581138830084189*alpha[14]*fUpwind[58]+0.1581138830084189*alpha[13]*fUpwind[57]+0.1581138830084189*alpha[12]*fUpwind[56]+0.1581138830084189*alpha[5]*fUpwind[52]+0.1581138830084189*alpha[19]*fUpwind[46]+0.1581138830084189*alpha[18]*fUpwind[45]+0.1581138830084189*alpha[11]*fUpwind[42]+0.1581138830084189*alpha[10]*fUpwind[41]+0.1581138830084189*alpha[9]*fUpwind[40]+0.1581138830084189*alpha[4]*fUpwind[36]+0.1767766952966368*alpha[8]*fUpwind[30]+0.1767766952966368*alpha[7]*fUpwind[29]+0.1767766952966368*alpha[6]*fUpwind[28]+0.1767766952966368*alpha[3]*fUpwind[25]+0.1767766952966368*alpha[2]*fUpwind[24]+0.1767766952966368*alpha[1]*fUpwind[23]+0.1767766952966368*alpha[19]*fUpwind[22]+0.1767766952966368*fUpwind[19]*alpha[22]+0.1767766952966368*alpha[18]*fUpwind[21]+0.1767766952966368*fUpwind[17]*alpha[20]+0.1767766952966368*alpha[0]*fUpwind[15]+0.1767766952966368*alpha[11]*fUpwind[14]+0.1767766952966368*fUpwind[11]*alpha[14]+0.1767766952966368*alpha[10]*fUpwind[13]+0.1767766952966368*fUpwind[10]*alpha[13]+0.1767766952966368*alpha[9]*fUpwind[12]+0.1767766952966368*fUpwind[9]*alpha[12]+0.1767766952966368*alpha[4]*fUpwind[5]+0.1767766952966368*fUpwind[4]*alpha[5]; 
  Ghat[16] = 0.1767766952966368*alpha[5]*fUpwind[27]+0.1767766952966368*alpha[4]*fUpwind[26]+0.1767766952966368*alpha[12]*fUpwind[22]+0.1767766952966368*fUpwind[12]*alpha[22]+0.1767766952966368*alpha[13]*fUpwind[21]+0.1767766952966368*alpha[14]*fUpwind[20]+0.1767766952966368*fUpwind[14]*alpha[20]+0.1767766952966368*alpha[9]*fUpwind[19]+0.1767766952966368*fUpwind[9]*alpha[19]+0.1767766952966368*alpha[10]*fUpwind[18]+0.1767766952966368*fUpwind[10]*alpha[18]+0.1767766952966368*alpha[11]*fUpwind[17]+0.1767766952966368*alpha[0]*fUpwind[16]+0.1767766952966368*alpha[1]*fUpwind[8]+0.1767766952966368*fUpwind[1]*alpha[8]+0.1767766952966368*alpha[2]*fUpwind[7]+0.1767766952966368*fUpwind[2]*alpha[7]+0.1767766952966368*alpha[3]*fUpwind[6]+0.1767766952966368*fUpwind[3]*alpha[6]; 
  Ghat[17] = 0.1581138830084189*alpha[11]*fUpwind[43]+0.1581138830084189*alpha[18]*fUpwind[39]+0.1581138830084189*alpha[19]*fUpwind[38]+0.1581138830084189*alpha[4]*fUpwind[37]+0.1581138830084189*alpha[9]*fUpwind[34]+0.1581138830084189*alpha[10]*fUpwind[33]+0.1767766952966368*alpha[14]*fUpwind[31]+0.1767766952966368*alpha[22]*fUpwind[29]+0.1767766952966368*alpha[5]*fUpwind[28]+0.1767766952966368*alpha[3]*fUpwind[26]+0.1767766952966368*alpha[12]*fUpwind[24]+0.1767766952966368*alpha[13]*fUpwind[23]+0.1767766952966368*fUpwind[15]*alpha[20]+0.1767766952966368*alpha[7]*fUpwind[19]+0.1767766952966368*fUpwind[7]*alpha[19]+0.1767766952966368*alpha[8]*fUpwind[18]+0.1767766952966368*fUpwind[8]*alpha[18]+0.1767766952966368*alpha[0]*fUpwind[17]+0.1767766952966368*alpha[11]*fUpwind[16]+0.1767766952966368*alpha[1]*fUpwind[10]+0.1767766952966368*fUpwind[1]*alpha[10]+0.1767766952966368*alpha[2]*fUpwind[9]+0.1767766952966368*fUpwind[2]*alpha[9]+0.1767766952966368*alpha[4]*fUpwind[6]+0.1767766952966368*fUpwind[4]*alpha[6]; 
  Ghat[18] = 0.1581138830084189*alpha[10]*fUpwind[43]+0.1581138830084189*alpha[4]*fUpwind[38]+0.1581138830084189*alpha[19]*fUpwind[37]+0.1581138830084189*alpha[9]*fUpwind[35]+0.1581138830084189*alpha[11]*fUpwind[33]+0.1581138830084189*alpha[18]*fUpwind[32]+0.1767766952966368*alpha[13]*fUpwind[31]+0.1767766952966368*alpha[20]*fUpwind[30]+0.1767766952966368*alpha[5]*fUpwind[29]+0.1767766952966368*alpha[22]*fUpwind[28]+0.1767766952966368*alpha[2]*fUpwind[26]+0.1767766952966368*alpha[12]*fUpwind[25]+0.1767766952966368*alpha[14]*fUpwind[23]+0.1767766952966368*alpha[6]*fUpwind[19]+0.1767766952966368*fUpwind[6]*alpha[19]+0.1767766952966368*alpha[0]*fUpwind[18]+0.1767766952966368*fUpwind[0]*alpha[18]+0.1767766952966368*alpha[8]*fUpwind[17]+0.1767766952966368*alpha[10]*fUpwind[16]+0.1767766952966368*alpha[1]*fUpwind[11]+0.1767766952966368*fUpwind[1]*alpha[11]+0.1767766952966368*alpha[3]*fUpwind[9]+0.1767766952966368*fUpwind[3]*alpha[9]+0.1767766952966368*alpha[4]*fUpwind[7]+0.1767766952966368*fUpwind[4]*alpha[7]; 
  Ghat[19] = 0.1581138830084189*alpha[9]*fUpwind[43]+0.1581138830084189*alpha[4]*fUpwind[39]+0.1581138830084189*alpha[18]*fUpwind[37]+0.1581138830084189*alpha[10]*fUpwind[35]+0.1581138830084189*alpha[11]*fUpwind[34]+0.1581138830084189*alpha[19]*fUpwind[32]+0.1767766952966368*alpha[12]*fUpwind[31]+0.1767766952966368*alpha[5]*fUpwind[30]+0.1767766952966368*alpha[20]*fUpwind[29]+0.1767766952966368*alpha[1]*fUpwind[26]+0.1767766952966368*alpha[13]*fUpwind[25]+0.1767766952966368*alpha[14]*fUpwind[24]+0.1767766952966368*fUpwind[15]*alpha[22]+0.1767766952966368*alpha[0]*fUpwind[19]+0.1767766952966368*fUpwind[0]*alpha[19]+0.1767766952966368*alpha[6]*fUpwind[18]+0.1767766952966368*fUpwind[6]*alpha[18]+0.1767766952966368*alpha[7]*fUpwind[17]+0.1767766952966368*alpha[9]*fUpwind[16]+0.1767766952966368*alpha[2]*fUpwind[11]+0.1767766952966368*fUpwind[2]*alpha[11]+0.1767766952966368*alpha[3]*fUpwind[10]+0.1767766952966368*fUpwind[3]*alpha[10]+0.1767766952966368*alpha[4]*fUpwind[8]+0.1767766952966368*fUpwind[4]*alpha[8]; 
  Ghat[20] = 0.1581138830084189*alpha[14]*fUpwind[59]+0.1581138830084189*alpha[22]*fUpwind[54]+0.1581138830084189*alpha[5]*fUpwind[53]+0.1581138830084189*alpha[12]*fUpwind[50]+0.1581138830084189*alpha[13]*fUpwind[49]+0.1581138830084189*alpha[20]*fUpwind[48]+0.1767766952966368*alpha[11]*fUpwind[31]+0.1767766952966368*alpha[18]*fUpwind[30]+0.1767766952966368*alpha[19]*fUpwind[29]+0.1767766952966368*alpha[4]*fUpwind[28]+0.1767766952966368*alpha[3]*fUpwind[27]+0.1767766952966368*alpha[9]*fUpwind[24]+0.1767766952966368*alpha[10]*fUpwind[23]+0.1767766952966368*alpha[7]*fUpwind[22]+0.1767766952966368*fUpwind[7]*alpha[22]+0.1767766952966368*alpha[8]*fUpwind[21]+0.1767766952966368*alpha[0]*fUpwind[20]+0.1767766952966368*fUpwind[0]*alpha[20]+0.1767766952966368*alpha[14]*fUpwind[16]+0.1767766952966368*alpha[1]*fUpwind[13]+0.1767766952966368*fUpwind[1]*alpha[13]+0.1767766952966368*alpha[2]*fUpwind[12]+0.1767766952966368*fUpwind[2]*alpha[12]+0.1767766952966368*alpha[5]*fUpwind[6]+0.1767766952966368*fUpwind[5]*alpha[6]; 
  Ghat[21] = 0.1581138830084189*alpha[13]*fUpwind[59]+0.1581138830084189*alpha[20]*fUpwind[55]+0.1581138830084189*alpha[5]*fUpwind[54]+0.1581138830084189*alpha[22]*fUpwind[53]+0.1581138830084189*alpha[12]*fUpwind[51]+0.1581138830084189*alpha[14]*fUpwind[49]+0.1767766952966368*alpha[10]*fUpwind[31]+0.1767766952966368*alpha[4]*fUpwind[29]+0.1767766952966368*alpha[19]*fUpwind[28]+0.1767766952966368*alpha[2]*fUpwind[27]+0.1767766952966368*alpha[9]*fUpwind[25]+0.1767766952966368*alpha[11]*fUpwind[23]+0.1767766952966368*alpha[6]*fUpwind[22]+0.1767766952966368*fUpwind[6]*alpha[22]+0.1767766952966368*alpha[0]*fUpwind[21]+0.1767766952966368*alpha[8]*fUpwind[20]+0.1767766952966368*fUpwind[8]*alpha[20]+0.1767766952966368*fUpwind[15]*alpha[18]+0.1767766952966368*alpha[13]*fUpwind[16]+0.1767766952966368*alpha[1]*fUpwind[14]+0.1767766952966368*fUpwind[1]*alpha[14]+0.1767766952966368*alpha[3]*fUpwind[12]+0.1767766952966368*fUpwind[3]*alpha[12]+0.1767766952966368*alpha[5]*fUpwind[7]+0.1767766952966368*fUpwind[5]*alpha[7]; 
  Ghat[22] = 0.1581138830084189*alpha[12]*fUpwind[59]+0.1581138830084189*alpha[5]*fUpwind[55]+0.1581138830084189*alpha[20]*fUpwind[54]+0.1581138830084189*alpha[13]*fUpwind[51]+0.1581138830084189*alpha[14]*fUpwind[50]+0.1581138830084189*alpha[22]*fUpwind[48]+0.1767766952966368*alpha[9]*fUpwind[31]+0.1767766952966368*alpha[4]*fUpwind[30]+0.1767766952966368*alpha[18]*fUpwind[28]+0.1767766952966368*alpha[1]*fUpwind[27]+0.1767766952966368*alpha[10]*fUpwind[25]+0.1767766952966368*alpha[11]*fUpwind[24]+0.1767766952966368*alpha[0]*fUpwind[22]+0.1767766952966368*fUpwind[0]*alpha[22]+0.1767766952966368*alpha[6]*fUpwind[21]+0.1767766952966368*alpha[7]*fUpwind[20]+0.1767766952966368*fUpwind[7]*alpha[20]+0.1767766952966368*fUpwind[15]*alpha[19]+0.1767766952966368*alpha[12]*fUpwind[16]+0.1767766952966368*alpha[2]*fUpwind[14]+0.1767766952966368*fUpwind[2]*alpha[14]+0.1767766952966368*alpha[3]*fUpwind[13]+0.1767766952966368*fUpwind[3]*alpha[13]+0.1767766952966368*alpha[5]*fUpwind[8]+0.1767766952966368*fUpwind[5]*alpha[8]; 
  Ghat[23] = 0.1581138830084189*alpha[22]*fUpwind[63]+0.1581138830084189*alpha[14]*fUpwind[61]+0.1581138830084189*alpha[13]*fUpwind[60]+0.1581138830084189*alpha[20]*fUpwind[57]+0.1581138830084189*alpha[5]*fUpwind[56]+0.1581138830084189*alpha[12]*fUpwind[52]+0.1581138830084189*alpha[19]*fUpwind[47]+0.1581138830084189*alpha[11]*fUpwind[45]+0.1581138830084189*alpha[10]*fUpwind[44]+0.1581138830084189*alpha[18]*fUpwind[42]+0.1581138830084189*alpha[4]*fUpwind[40]+0.1581138830084189*alpha[9]*fUpwind[36]+0.1767766952966368*alpha[8]*fUpwind[31]+0.1767766952966368*alpha[3]*fUpwind[29]+0.1767766952966368*alpha[2]*fUpwind[28]+0.1767766952966368*alpha[19]*fUpwind[27]+0.1767766952966368*alpha[22]*fUpwind[26]+0.1767766952966368*alpha[7]*fUpwind[25]+0.1767766952966368*alpha[6]*fUpwind[24]+0.1767766952966368*alpha[0]*fUpwind[23]+0.1767766952966368*alpha[11]*fUpwind[21]+0.1767766952966368*alpha[10]*fUpwind[20]+0.1767766952966368*fUpwind[10]*alpha[20]+0.1767766952966368*alpha[14]*fUpwind[18]+0.1767766952966368*fUpwind[14]*alpha[18]+0.1767766952966368*alpha[13]*fUpwind[17]+0.1767766952966368*alpha[1]*fUpwind[15]+0.1767766952966368*alpha[4]*fUpwind[12]+0.1767766952966368*fUpwind[4]*alpha[12]+0.1767766952966368*alpha[5]*fUpwind[9]+0.1767766952966368*fUpwind[5]*alpha[9]; 
  Ghat[24] = 0.1581138830084189*alpha[14]*fUpwind[62]+0.1581138830084189*alpha[12]*fUpwind[60]+0.1581138830084189*alpha[22]*fUpwind[58]+0.1581138830084189*alpha[5]*fUpwind[57]+0.1581138830084189*alpha[20]*fUpwind[56]+0.1581138830084189*alpha[13]*fUpwind[52]+0.1581138830084189*alpha[18]*fUpwind[47]+0.1581138830084189*alpha[11]*fUpwind[46]+0.1581138830084189*alpha[9]*fUpwind[44]+0.1581138830084189*alpha[19]*fUpwind[42]+0.1581138830084189*alpha[4]*fUpwind[41]+0.1581138830084189*alpha[10]*fUpwind[36]+0.1767766952966368*alpha[7]*fUpwind[31]+0.1767766952966368*alpha[3]*fUpwind[30]+0.1767766952966368*alpha[1]*fUpwind[28]+0.1767766952966368*alpha[18]*fUpwind[27]+0.1767766952966368*alpha[8]*fUpwind[25]+0.1767766952966368*alpha[0]*fUpwind[24]+0.1767766952966368*alpha[6]*fUpwind[23]+0.1767766952966368*alpha[11]*fUpwind[22]+0.1767766952966368*fUpwind[11]*alpha[22]+0.1767766952966368*alpha[9]*fUpwind[20]+0.1767766952966368*fUpwind[9]*alpha[20]+0.1767766952966368*alpha[14]*fUpwind[19]+0.1767766952966368*fUpwind[14]*alpha[19]+0.1767766952966368*alpha[12]*fUpwind[17]+0.1767766952966368*alpha[2]*fUpwind[15]+0.1767766952966368*alpha[4]*fUpwind[13]+0.1767766952966368*fUpwind[4]*alpha[13]+0.1767766952966368*alpha[5]*fUpwind[10]+0.1767766952966368*fUpwind[5]*alpha[10]; 
  Ghat[25] = 0.1581138830084189*alpha[20]*fUpwind[63]+0.1581138830084189*alpha[13]*fUpwind[62]+0.1581138830084189*alpha[12]*fUpwind[61]+0.1581138830084189*alpha[5]*fUpwind[58]+0.1581138830084189*alpha[22]*fUpwind[57]+0.1581138830084189*alpha[14]*fUpwind[52]+0.1581138830084189*alpha[10]*fUpwind[46]+0.1581138830084189*alpha[9]*fUpwind[45]+0.1581138830084189*alpha[4]*fUpwind[42]+0.1581138830084189*alpha[19]*fUpwind[41]+0.1581138830084189*alpha[18]*fUpwind[40]+0.1581138830084189*alpha[11]*fUpwind[36]+0.1767766952966368*alpha[6]*fUpwind[31]+0.1767766952966368*alpha[2]*fUpwind[30]+0.1767766952966368*alpha[1]*fUpwind[29]+0.1767766952966368*alpha[20]*fUpwind[26]+0.1767766952966368*alpha[0]*fUpwind[25]+0.1767766952966368*alpha[8]*fUpwind[24]+0.1767766952966368*alpha[7]*fUpwind[23]+0.1767766952966368*alpha[10]*fUpwind[22]+0.1767766952966368*fUpwind[10]*alpha[22]+0.1767766952966368*alpha[9]*fUpwind[21]+0.1767766952966368*alpha[13]*fUpwind[19]+0.1767766952966368*fUpwind[13]*alpha[19]+0.1767766952966368*alpha[12]*fUpwind[18]+0.1767766952966368*fUpwind[12]*alpha[18]+0.1767766952966368*alpha[3]*fUpwind[15]+0.1767766952966368*alpha[4]*fUpwind[14]+0.1767766952966368*fUpwind[4]*alpha[14]+0.1767766952966368*alpha[5]*fUpwind[11]+0.1767766952966368*fUpwind[5]*alpha[11]; 
  Ghat[26] = 0.1581138830084189*alpha[4]*fUpwind[43]+0.1581138830084189*alpha[9]*fUpwind[39]+0.1581138830084189*alpha[10]*fUpwind[38]+0.1581138830084189*alpha[11]*fUpwind[37]+0.1581138830084189*alpha[18]*fUpwind[34]+0.1581138830084189*alpha[19]*fUpwind[33]+0.1767766952966368*alpha[5]*fUpwind[31]+0.1767766952966368*alpha[12]*fUpwind[30]+0.1767766952966368*alpha[13]*fUpwind[29]+0.1767766952966368*alpha[14]*fUpwind[28]+0.1767766952966368*alpha[0]*fUpwind[26]+0.1767766952966368*alpha[20]*fUpwind[25]+0.1767766952966368*alpha[22]*fUpwind[23]+0.1767766952966368*alpha[1]*fUpwind[19]+0.1767766952966368*fUpwind[1]*alpha[19]+0.1767766952966368*alpha[2]*fUpwind[18]+0.1767766952966368*fUpwind[2]*alpha[18]+0.1767766952966368*alpha[3]*fUpwind[17]+0.1767766952966368*alpha[4]*fUpwind[16]+0.1767766952966368*alpha[6]*fUpwind[11]+0.1767766952966368*fUpwind[6]*alpha[11]+0.1767766952966368*alpha[7]*fUpwind[10]+0.1767766952966368*fUpwind[7]*alpha[10]+0.1767766952966368*alpha[8]*fUpwind[9]+0.1767766952966368*fUpwind[8]*alpha[9]; 
  Ghat[27] = 0.1581138830084189*alpha[5]*fUpwind[59]+0.1581138830084189*alpha[12]*fUpwind[55]+0.1581138830084189*alpha[13]*fUpwind[54]+0.1581138830084189*alpha[14]*fUpwind[53]+0.1581138830084189*alpha[20]*fUpwind[51]+0.1581138830084189*alpha[22]*fUpwind[49]+0.1767766952966368*alpha[4]*fUpwind[31]+0.1767766952966368*alpha[9]*fUpwind[30]+0.1767766952966368*alpha[10]*fUpwind[29]+0.1767766952966368*alpha[11]*fUpwind[28]+0.1767766952966368*alpha[0]*fUpwind[27]+0.1767766952966368*alpha[18]*fUpwind[24]+0.1767766952966368*alpha[19]*fUpwind[23]+0.1767766952966368*alpha[1]*fUpwind[22]+0.1767766952966368*fUpwind[1]*alpha[22]+0.1767766952966368*alpha[2]*fUpwind[21]+0.1767766952966368*alpha[3]*fUpwind[20]+0.1767766952966368*fUpwind[3]*alpha[20]+0.1767766952966368*alpha[5]*fUpwind[16]+0.1767766952966368*alpha[6]*fUpwind[14]+0.1767766952966368*fUpwind[6]*alpha[14]+0.1767766952966368*alpha[7]*fUpwind[13]+0.1767766952966368*fUpwind[7]*alpha[13]+0.1767766952966368*alpha[8]*fUpwind[12]+0.1767766952966368*fUpwind[8]*alpha[12]; 
  Ghat[28] = 0.1581138830084189*alpha[14]*fUpwind[63]+0.1581138830084189*alpha[22]*fUpwind[61]+0.1581138830084189*alpha[5]*fUpwind[60]+0.1581138830084189*alpha[12]*fUpwind[57]+0.1581138830084189*alpha[13]*fUpwind[56]+0.1581138830084189*alpha[20]*fUpwind[52]+0.1581138830084189*alpha[11]*fUpwind[47]+0.1581138830084189*alpha[18]*fUpwind[46]+0.1581138830084189*alpha[19]*fUpwind[45]+0.1581138830084189*alpha[4]*fUpwind[44]+0.1581138830084189*alpha[9]*fUpwind[41]+0.1581138830084189*alpha[10]*fUpwind[40]+0.1767766952966368*alpha[3]*fUpwind[31]+0.1767766952966368*alpha[7]*fUpwind[30]+0.1767766952966368*alpha[8]*fUpwind[29]+0.1767766952966368*alpha[0]*fUpwind[28]+0.1767766952966368*alpha[11]*fUpwind[27]+0.1767766952966368*alpha[14]*fUpwind[26]+0.1767766952966368*alpha[1]*fUpwind[24]+0.1767766952966368*alpha[2]*fUpwind[23]+0.1767766952966368*alpha[18]*fUpwind[22]+0.1767766952966368*fUpwind[18]*alpha[22]+0.1767766952966368*alpha[19]*fUpwind[21]+0.1767766952966368*alpha[4]*fUpwind[20]+0.1767766952966368*fUpwind[4]*alpha[20]+0.1767766952966368*alpha[5]*fUpwind[17]+0.1767766952966368*alpha[6]*fUpwind[15]+0.1767766952966368*alpha[9]*fUpwind[13]+0.1767766952966368*fUpwind[9]*alpha[13]+0.1767766952966368*alpha[10]*fUpwind[12]+0.1767766952966368*fUpwind[10]*alpha[12]; 
  Ghat[29] = 0.1581138830084189*alpha[13]*fUpwind[63]+0.1581138830084189*alpha[20]*fUpwind[62]+0.1581138830084189*alpha[5]*fUpwind[61]+0.1581138830084189*alpha[22]*fUpwind[60]+0.1581138830084189*alpha[12]*fUpwind[58]+0.1581138830084189*alpha[14]*fUpwind[56]+0.1581138830084189*alpha[10]*fUpwind[47]+0.1581138830084189*alpha[4]*fUpwind[45]+0.1581138830084189*alpha[19]*fUpwind[44]+0.1581138830084189*alpha[9]*fUpwind[42]+0.1581138830084189*alpha[11]*fUpwind[40]+0.1581138830084189*alpha[18]*fUpwind[36]+0.1767766952966368*alpha[2]*fUpwind[31]+0.1767766952966368*alpha[6]*fUpwind[30]+0.1767766952966368*alpha[0]*fUpwind[29]+0.1767766952966368*alpha[8]*fUpwind[28]+0.1767766952966368*alpha[10]*fUpwind[27]+0.1767766952966368*alpha[13]*fUpwind[26]+0.1767766952966368*alpha[1]*fUpwind[25]+0.1767766952966368*alpha[3]*fUpwind[23]+0.1767766952966368*fUpwind[17]*alpha[22]+0.1767766952966368*alpha[4]*fUpwind[21]+0.1767766952966368*alpha[19]*fUpwind[20]+0.1767766952966368*fUpwind[19]*alpha[20]+0.1767766952966368*alpha[5]*fUpwind[18]+0.1767766952966368*fUpwind[5]*alpha[18]+0.1767766952966368*alpha[7]*fUpwind[15]+0.1767766952966368*alpha[9]*fUpwind[14]+0.1767766952966368*fUpwind[9]*alpha[14]+0.1767766952966368*alpha[11]*fUpwind[12]+0.1767766952966368*fUpwind[11]*alpha[12]; 
  Ghat[30] = 0.1581138830084189*alpha[12]*fUpwind[63]+0.1581138830084189*alpha[5]*fUpwind[62]+0.1581138830084189*alpha[20]*fUpwind[61]+0.1581138830084189*alpha[13]*fUpwind[58]+0.1581138830084189*alpha[14]*fUpwind[57]+0.1581138830084189*alpha[22]*fUpwind[52]+0.1581138830084189*alpha[9]*fUpwind[47]+0.1581138830084189*alpha[4]*fUpwind[46]+0.1581138830084189*alpha[18]*fUpwind[44]+0.1581138830084189*alpha[10]*fUpwind[42]+0.1581138830084189*alpha[11]*fUpwind[41]+0.1581138830084189*alpha[19]*fUpwind[36]+0.1767766952966368*alpha[1]*fUpwind[31]+0.1767766952966368*alpha[0]*fUpwind[30]+0.1767766952966368*alpha[6]*fUpwind[29]+0.1767766952966368*alpha[7]*fUpwind[28]+0.1767766952966368*alpha[9]*fUpwind[27]+0.1767766952966368*alpha[12]*fUpwind[26]+0.1767766952966368*alpha[2]*fUpwind[25]+0.1767766952966368*alpha[3]*fUpwind[24]+0.1767766952966368*alpha[4]*fUpwind[22]+0.1767766952966368*fUpwind[4]*alpha[22]+0.1767766952966368*alpha[18]*fUpwind[20]+0.1767766952966368*fUpwind[18]*alpha[20]+0.1767766952966368*alpha[5]*fUpwind[19]+0.1767766952966368*fUpwind[5]*alpha[19]+0.1767766952966368*alpha[8]*fUpwind[15]+0.1767766952966368*alpha[10]*fUpwind[14]+0.1767766952966368*fUpwind[10]*alpha[14]+0.1767766952966368*alpha[11]*fUpwind[13]+0.1767766952966368*fUpwind[11]*alpha[13]; 
  Ghat[31] = 0.1581138830084189*alpha[5]*fUpwind[63]+0.1581138830084189*alpha[12]*fUpwind[62]+0.1581138830084189*alpha[13]*fUpwind[61]+0.1581138830084189*alpha[14]*fUpwind[60]+0.1581138830084189*alpha[20]*fUpwind[58]+0.1581138830084189*alpha[22]*fUpwind[56]+0.1581138830084189*alpha[4]*fUpwind[47]+0.1581138830084189*alpha[9]*fUpwind[46]+0.1581138830084189*alpha[10]*fUpwind[45]+0.1581138830084189*alpha[11]*fUpwind[44]+0.1581138830084189*alpha[18]*fUpwind[41]+0.1581138830084189*alpha[19]*fUpwind[40]+0.1767766952966368*alpha[0]*fUpwind[31]+0.1767766952966368*alpha[1]*fUpwind[30]+0.1767766952966368*alpha[2]*fUpwind[29]+0.1767766952966368*alpha[3]*fUpwind[28]+0.1767766952966368*alpha[4]*fUpwind[27]+0.1767766952966368*alpha[5]*fUpwind[26]+0.1767766952966368*alpha[6]*fUpwind[25]+0.1767766952966368*alpha[7]*fUpwind[24]+0.1767766952966368*alpha[8]*fUpwind[23]+0.1767766952966368*alpha[9]*fUpwind[22]+0.1767766952966368*fUpwind[9]*alpha[22]+0.1767766952966368*alpha[10]*fUpwind[21]+0.1767766952966368*alpha[11]*fUpwind[20]+0.1767766952966368*fUpwind[11]*alpha[20]+0.1767766952966368*alpha[12]*fUpwind[19]+0.1767766952966368*fUpwind[12]*alpha[19]+0.1767766952966368*alpha[13]*fUpwind[18]+0.1767766952966368*fUpwind[13]*alpha[18]+0.1767766952966368*alpha[14]*fUpwind[17]; 
  Ghat[32] = 0.1767766952966368*alpha[22]*fUpwind[46]+0.1767766952966368*alpha[20]*fUpwind[44]+0.1767766952966368*alpha[14]*fUpwind[42]+0.1767766952966368*alpha[13]*fUpwind[41]+0.1767766952966368*alpha[12]*fUpwind[40]+0.1767766952966368*alpha[8]*fUpwind[39]+0.1767766952966368*alpha[7]*fUpwind[38]+0.1767766952966368*alpha[6]*fUpwind[37]+0.1767766952966368*alpha[5]*fUpwind[36]+0.1767766952966368*alpha[3]*fUpwind[35]+0.1767766952966368*alpha[2]*fUpwind[34]+0.1767766952966368*alpha[1]*fUpwind[33]+0.1767766952966368*alpha[0]*fUpwind[32]+0.1581138830084189*alpha[19]*fUpwind[19]+0.1581138830084189*alpha[18]*fUpwind[18]+0.1581138830084189*alpha[11]*fUpwind[11]+0.1581138830084189*alpha[10]*fUpwind[10]+0.1581138830084189*alpha[9]*fUpwind[9]+0.1581138830084189*alpha[4]*fUpwind[4]; 
  Ghat[33] = 0.1767766952966368*alpha[22]*fUpwind[47]+0.1767766952966368*alpha[14]*fUpwind[45]+0.1767766952966368*alpha[13]*fUpwind[44]+0.1767766952966368*alpha[8]*fUpwind[43]+0.1767766952966368*alpha[20]*fUpwind[41]+0.1767766952966368*alpha[5]*fUpwind[40]+0.1767766952966368*alpha[3]*fUpwind[38]+0.1767766952966368*alpha[2]*fUpwind[37]+0.1767766952966368*alpha[12]*fUpwind[36]+0.1767766952966368*alpha[7]*fUpwind[35]+0.1767766952966368*alpha[6]*fUpwind[34]+0.1767766952966368*alpha[0]*fUpwind[33]+0.1767766952966368*alpha[1]*fUpwind[32]+0.1581138830084189*alpha[19]*fUpwind[26]+0.1581138830084189*alpha[11]*fUpwind[18]+0.1581138830084189*fUpwind[11]*alpha[18]+0.1581138830084189*alpha[10]*fUpwind[17]+0.1581138830084189*alpha[4]*fUpwind[9]+0.1581138830084189*fUpwind[4]*alpha[9]; 
  Ghat[34] = 0.1767766952966368*alpha[14]*fUpwind[46]+0.1767766952966368*alpha[12]*fUpwind[44]+0.1767766952966368*alpha[7]*fUpwind[43]+0.1767766952966368*alpha[22]*fUpwind[42]+0.1767766952966368*alpha[5]*fUpwind[41]+0.1767766952966368*alpha[20]*fUpwind[40]+0.1767766952966368*alpha[3]*fUpwind[39]+0.1767766952966368*alpha[1]*fUpwind[37]+0.1767766952966368*alpha[13]*fUpwind[36]+0.1767766952966368*alpha[8]*fUpwind[35]+0.1767766952966368*alpha[0]*fUpwind[34]+0.1767766952966368*alpha[6]*fUpwind[33]+0.1767766952966368*alpha[2]*fUpwind[32]+0.1581138830084189*alpha[18]*fUpwind[26]+0.1581138830084189*alpha[11]*fUpwind[19]+0.1581138830084189*fUpwind[11]*alpha[19]+0.1581138830084189*alpha[9]*fUpwind[17]+0.1581138830084189*alpha[4]*fUpwind[10]+0.1581138830084189*fUpwind[4]*alpha[10]; 
  Ghat[35] = 0.1767766952966368*alpha[20]*fUpwind[47]+0.1767766952966368*alpha[13]*fUpwind[46]+0.1767766952966368*alpha[12]*fUpwind[45]+0.1767766952966368*alpha[6]*fUpwind[43]+0.1767766952966368*alpha[5]*fUpwind[42]+0.1767766952966368*alpha[22]*fUpwind[41]+0.1767766952966368*alpha[2]*fUpwind[39]+0.1767766952966368*alpha[1]*fUpwind[38]+0.1767766952966368*alpha[14]*fUpwind[36]+0.1767766952966368*alpha[0]*fUpwind[35]+0.1767766952966368*alpha[8]*fUpwind[34]+0.1767766952966368*alpha[7]*fUpwind[33]+0.1767766952966368*alpha[3]*fUpwind[32]+0.1581138830084189*alpha[10]*fUpwind[19]+0.1581138830084189*fUpwind[10]*alpha[19]+0.1581138830084189*alpha[9]*fUpwind[18]+0.1581138830084189*fUpwind[9]*alpha[18]+0.1581138830084189*alpha[4]*fUpwind[11]+0.1581138830084189*fUpwind[4]*alpha[11]; 
  Ghat[36] = 0.1767766952966368*alpha[8]*fUpwind[46]+0.1767766952966368*alpha[7]*fUpwind[45]+0.1767766952966368*alpha[6]*fUpwind[44]+0.1767766952966368*alpha[3]*fUpwind[42]+0.1767766952966368*alpha[2]*fUpwind[41]+0.1767766952966368*alpha[1]*fUpwind[40]+0.1767766952966368*alpha[22]*fUpwind[39]+0.1767766952966368*alpha[20]*fUpwind[37]+0.1767766952966368*alpha[0]*fUpwind[36]+0.1767766952966368*alpha[14]*fUpwind[35]+0.1767766952966368*alpha[13]*fUpwind[34]+0.1767766952966368*alpha[12]*fUpwind[33]+0.1767766952966368*alpha[5]*fUpwind[32]+0.1581138830084189*alpha[19]*fUpwind[30]+0.1581138830084189*alpha[18]*fUpwind[29]+0.1581138830084189*alpha[11]*fUpwind[25]+0.1581138830084189*alpha[10]*fUpwind[24]+0.1581138830084189*alpha[9]*fUpwind[23]+0.1581138830084189*alpha[4]*fUpwind[15]; 
  Ghat[37] = 0.1767766952966368*alpha[14]*fUpwind[47]+0.1767766952966368*alpha[22]*fUpwind[45]+0.1767766952966368*alpha[5]*fUpwind[44]+0.1767766952966368*alpha[3]*fUpwind[43]+0.1767766952966368*alpha[12]*fUpwind[41]+0.1767766952966368*alpha[13]*fUpwind[40]+0.1767766952966368*alpha[7]*fUpwind[39]+0.1767766952966368*alpha[8]*fUpwind[38]+0.1767766952966368*alpha[0]*fUpwind[37]+0.1767766952966368*alpha[20]*fUpwind[36]+0.1767766952966368*alpha[1]*fUpwind[34]+0.1767766952966368*alpha[2]*fUpwind[33]+0.1767766952966368*alpha[6]*fUpwind[32]+0.1581138830084189*alpha[11]*fUpwind[26]+0.1581138830084189*alpha[18]*fUpwind[19]+0.1581138830084189*fUpwind[18]*alpha[19]+0.1581138830084189*alpha[4]*fUpwind[17]+0.1581138830084189*alpha[9]*fUpwind[10]+0.1581138830084189*fUpwind[9]*alpha[10]; 
  Ghat[38] = 0.1767766952966368*alpha[13]*fUpwind[47]+0.1767766952966368*alpha[20]*fUpwind[46]+0.1767766952966368*alpha[5]*fUpwind[45]+0.1767766952966368*alpha[22]*fUpwind[44]+0.1767766952966368*alpha[2]*fUpwind[43]+0.1767766952966368*alpha[12]*fUpwind[42]+0.1767766952966368*alpha[14]*fUpwind[40]+0.1767766952966368*alpha[6]*fUpwind[39]+0.1767766952966368*alpha[0]*fUpwind[38]+0.1767766952966368*alpha[8]*fUpwind[37]+0.1767766952966368*alpha[1]*fUpwind[35]+0.1767766952966368*alpha[3]*fUpwind[33]+0.1767766952966368*alpha[7]*fUpwind[32]+0.1581138830084189*alpha[10]*fUpwind[26]+0.1581138830084189*fUpwind[17]*alpha[19]+0.1581138830084189*alpha[4]*fUpwind[18]+0.1581138830084189*fUpwind[4]*alpha[18]+0.1581138830084189*alpha[9]*fUpwind[11]+0.1581138830084189*fUpwind[9]*alpha[11]; 
  Ghat[39] = 0.1767766952966368*alpha[12]*fUpwind[47]+0.1767766952966368*alpha[5]*fUpwind[46]+0.1767766952966368*alpha[20]*fUpwind[45]+0.1767766952966368*alpha[1]*fUpwind[43]+0.1767766952966368*alpha[13]*fUpwind[42]+0.1767766952966368*alpha[14]*fUpwind[41]+0.1767766952966368*alpha[0]*fUpwind[39]+0.1767766952966368*alpha[6]*fUpwind[38]+0.1767766952966368*alpha[7]*fUpwind[37]+0.1767766952966368*alpha[22]*fUpwind[36]+0.1767766952966368*alpha[2]*fUpwind[35]+0.1767766952966368*alpha[3]*fUpwind[34]+0.1767766952966368*alpha[8]*fUpwind[32]+0.1581138830084189*alpha[9]*fUpwind[26]+0.1581138830084189*alpha[4]*fUpwind[19]+0.1581138830084189*fUpwind[4]*alpha[19]+0.1581138830084189*fUpwind[17]*alpha[18]+0.1581138830084189*alpha[10]*fUpwind[11]+0.1581138830084189*fUpwind[10]*alpha[11]; 
  Ghat[40] = 0.1767766952966368*alpha[8]*fUpwind[47]+0.1767766952966368*alpha[3]*fUpwind[45]+0.1767766952966368*alpha[2]*fUpwind[44]+0.1767766952966368*alpha[22]*fUpwind[43]+0.1767766952966368*alpha[7]*fUpwind[42]+0.1767766952966368*alpha[6]*fUpwind[41]+0.1767766952966368*alpha[0]*fUpwind[40]+0.1767766952966368*alpha[14]*fUpwind[38]+0.1767766952966368*alpha[13]*fUpwind[37]+0.1767766952966368*alpha[1]*fUpwind[36]+0.1767766952966368*alpha[20]*fUpwind[34]+0.1767766952966368*alpha[5]*fUpwind[33]+0.1767766952966368*alpha[12]*fUpwind[32]+0.1581138830084189*alpha[19]*fUpwind[31]+0.1581138830084189*alpha[11]*fUpwind[29]+0.1581138830084189*alpha[10]*fUpwind[28]+0.1581138830084189*alpha[18]*fUpwind[25]+0.1581138830084189*alpha[4]*fUpwind[23]+0.1581138830084189*alpha[9]*fUpwind[15]; 
  Ghat[41] = 0.1767766952966368*alpha[7]*fUpwind[47]+0.1767766952966368*alpha[3]*fUpwind[46]+0.1767766952966368*alpha[1]*fUpwind[44]+0.1767766952966368*alpha[8]*fUpwind[42]+0.1767766952966368*alpha[0]*fUpwind[41]+0.1767766952966368*alpha[6]*fUpwind[40]+0.1767766952966368*alpha[14]*fUpwind[39]+0.1767766952966368*alpha[12]*fUpwind[37]+0.1767766952966368*alpha[2]*fUpwind[36]+0.1767766952966368*alpha[22]*fUpwind[35]+0.1767766952966368*alpha[5]*fUpwind[34]+0.1767766952966368*alpha[20]*fUpwind[33]+0.1767766952966368*alpha[13]*fUpwind[32]+0.1581138830084189*alpha[18]*fUpwind[31]+0.1581138830084189*alpha[11]*fUpwind[30]+0.1581138830084189*alpha[9]*fUpwind[28]+0.1581138830084189*alpha[19]*fUpwind[25]+0.1581138830084189*alpha[4]*fUpwind[24]+0.1581138830084189*alpha[10]*fUpwind[15]; 
  Ghat[42] = 0.1767766952966368*alpha[6]*fUpwind[47]+0.1767766952966368*alpha[2]*fUpwind[46]+0.1767766952966368*alpha[1]*fUpwind[45]+0.1767766952966368*alpha[20]*fUpwind[43]+0.1767766952966368*alpha[0]*fUpwind[42]+0.1767766952966368*alpha[8]*fUpwind[41]+0.1767766952966368*alpha[7]*fUpwind[40]+0.1767766952966368*alpha[13]*fUpwind[39]+0.1767766952966368*alpha[12]*fUpwind[38]+0.1767766952966368*alpha[3]*fUpwind[36]+0.1767766952966368*alpha[5]*fUpwind[35]+0.1767766952966368*alpha[22]*fUpwind[34]+0.1767766952966368*alpha[14]*fUpwind[32]+0.1581138830084189*alpha[10]*fUpwind[30]+0.1581138830084189*alpha[9]*fUpwind[29]+0.1581138830084189*alpha[4]*fUpwind[25]+0.1581138830084189*alpha[19]*fUpwind[24]+0.1581138830084189*alpha[18]*fUpwind[23]+0.1581138830084189*alpha[11]*fUpwind[15]; 
  Ghat[43] = 0.1767766952966368*alpha[5]*fUpwind[47]+0.1767766952966368*alpha[12]*fUpwind[46]+0.1767766952966368*alpha[13]*fUpwind[45]+0.1767766952966368*alpha[14]*fUpwind[44]+0.1767766952966368*alpha[0]*fUpwind[43]+0.1767766952966368*alpha[20]*fUpwind[42]+0.1767766952966368*alpha[22]*fUpwind[40]+0.1767766952966368*alpha[1]*fUpwind[39]+0.1767766952966368*alpha[2]*fUpwind[38]+0.1767766952966368*alpha[3]*fUpwind[37]+0.1767766952966368*alpha[6]*fUpwind[35]+0.1767766952966368*alpha[7]*fUpwind[34]+0.1767766952966368*alpha[8]*fUpwind[33]+0.1581138830084189*alpha[4]*fUpwind[26]+0.1581138830084189*alpha[9]*fUpwind[19]+0.1581138830084189*fUpwind[9]*alpha[19]+0.1581138830084189*alpha[10]*fUpwind[18]+0.1581138830084189*fUpwind[10]*alpha[18]+0.1581138830084189*alpha[11]*fUpwind[17]; 
  Ghat[44] = 0.1767766952966368*alpha[3]*fUpwind[47]+0.1767766952966368*alpha[7]*fUpwind[46]+0.1767766952966368*alpha[8]*fUpwind[45]+0.1767766952966368*alpha[0]*fUpwind[44]+0.1767766952966368*alpha[14]*fUpwind[43]+0.1767766952966368*alpha[1]*fUpwind[41]+0.1767766952966368*alpha[2]*fUpwind[40]+0.1767766952966368*alpha[22]*fUpwind[38]+0.1767766952966368*alpha[5]*fUpwind[37]+0.1767766952966368*alpha[6]*fUpwind[36]+0.1767766952966368*alpha[12]*fUpwind[34]+0.1767766952966368*alpha[13]*fUpwind[33]+0.1767766952966368*alpha[20]*fUpwind[32]+0.1581138830084189*alpha[11]*fUpwind[31]+0.1581138830084189*alpha[18]*fUpwind[30]+0.1581138830084189*alpha[19]*fUpwind[29]+0.1581138830084189*alpha[4]*fUpwind[28]+0.1581138830084189*alpha[9]*fUpwind[24]+0.1581138830084189*alpha[10]*fUpwind[23]; 
  Ghat[45] = 0.1767766952966368*alpha[2]*fUpwind[47]+0.1767766952966368*alpha[6]*fUpwind[46]+0.1767766952966368*alpha[0]*fUpwind[45]+0.1767766952966368*alpha[8]*fUpwind[44]+0.1767766952966368*alpha[13]*fUpwind[43]+0.1767766952966368*alpha[1]*fUpwind[42]+0.1767766952966368*alpha[3]*fUpwind[40]+0.1767766952966368*alpha[20]*fUpwind[39]+0.1767766952966368*alpha[5]*fUpwind[38]+0.1767766952966368*alpha[22]*fUpwind[37]+0.1767766952966368*alpha[7]*fUpwind[36]+0.1767766952966368*alpha[12]*fUpwind[35]+0.1767766952966368*alpha[14]*fUpwind[33]+0.1581138830084189*alpha[10]*fUpwind[31]+0.1581138830084189*alpha[4]*fUpwind[29]+0.1581138830084189*alpha[19]*fUpwind[28]+0.1581138830084189*alpha[9]*fUpwind[25]+0.1581138830084189*alpha[11]*fUpwind[23]+0.1581138830084189*fUpwind[15]*alpha[18]; 
  Ghat[46] = 0.1767766952966368*alpha[1]*fUpwind[47]+0.1767766952966368*alpha[0]*fUpwind[46]+0.1767766952966368*alpha[6]*fUpwind[45]+0.1767766952966368*alpha[7]*fUpwind[44]+0.1767766952966368*alpha[12]*fUpwind[43]+0.1767766952966368*alpha[2]*fUpwind[42]+0.1767766952966368*alpha[3]*fUpwind[41]+0.1767766952966368*alpha[5]*fUpwind[39]+0.1767766952966368*alpha[20]*fUpwind[38]+0.1767766952966368*alpha[8]*fUpwind[36]+0.1767766952966368*alpha[13]*fUpwind[35]+0.1767766952966368*alpha[14]*fUpwind[34]+0.1767766952966368*alpha[22]*fUpwind[32]+0.1581138830084189*alpha[9]*fUpwind[31]+0.1581138830084189*alpha[4]*fUpwind[30]+0.1581138830084189*alpha[18]*fUpwind[28]+0.1581138830084189*alpha[10]*fUpwind[25]+0.1581138830084189*alpha[11]*fUpwind[24]+0.1581138830084189*fUpwind[15]*alpha[19]; 
  Ghat[47] = 0.1767766952966368*alpha[0]*fUpwind[47]+0.1767766952966368*alpha[1]*fUpwind[46]+0.1767766952966368*alpha[2]*fUpwind[45]+0.1767766952966368*alpha[3]*fUpwind[44]+0.1767766952966368*alpha[5]*fUpwind[43]+0.1767766952966368*alpha[6]*fUpwind[42]+0.1767766952966368*alpha[7]*fUpwind[41]+0.1767766952966368*alpha[8]*fUpwind[40]+0.1767766952966368*alpha[12]*fUpwind[39]+0.1767766952966368*alpha[13]*fUpwind[38]+0.1767766952966368*alpha[14]*fUpwind[37]+0.1767766952966368*alpha[20]*fUpwind[35]+0.1767766952966368*alpha[22]*fUpwind[33]+0.1581138830084189*alpha[4]*fUpwind[31]+0.1581138830084189*alpha[9]*fUpwind[30]+0.1581138830084189*alpha[10]*fUpwind[29]+0.1581138830084189*alpha[11]*fUpwind[28]+0.1581138830084189*alpha[18]*fUpwind[24]+0.1581138830084189*alpha[19]*fUpwind[23]; 
  Ghat[48] = 0.1767766952966368*alpha[19]*fUpwind[62]+0.1767766952966368*alpha[18]*fUpwind[61]+0.1767766952966368*alpha[11]*fUpwind[58]+0.1767766952966368*alpha[10]*fUpwind[57]+0.1767766952966368*alpha[9]*fUpwind[56]+0.1767766952966368*alpha[8]*fUpwind[55]+0.1767766952966368*alpha[7]*fUpwind[54]+0.1767766952966368*alpha[6]*fUpwind[53]+0.1767766952966368*alpha[4]*fUpwind[52]+0.1767766952966368*alpha[3]*fUpwind[51]+0.1767766952966368*alpha[2]*fUpwind[50]+0.1767766952966368*alpha[1]*fUpwind[49]+0.1767766952966368*alpha[0]*fUpwind[48]+0.1581138830084189*alpha[22]*fUpwind[22]+0.1581138830084189*alpha[20]*fUpwind[20]+0.1581138830084189*alpha[14]*fUpwind[14]+0.1581138830084189*alpha[13]*fUpwind[13]+0.1581138830084189*alpha[12]*fUpwind[12]+0.1581138830084189*alpha[5]*fUpwind[5]; 
  Ghat[49] = 0.1767766952966368*alpha[19]*fUpwind[63]+0.1767766952966368*alpha[11]*fUpwind[61]+0.1767766952966368*alpha[10]*fUpwind[60]+0.1767766952966368*alpha[8]*fUpwind[59]+0.1767766952966368*alpha[18]*fUpwind[58]+0.1767766952966368*alpha[4]*fUpwind[56]+0.1767766952966368*alpha[3]*fUpwind[54]+0.1767766952966368*alpha[2]*fUpwind[53]+0.1767766952966368*alpha[9]*fUpwind[52]+0.1767766952966368*alpha[7]*fUpwind[51]+0.1767766952966368*alpha[6]*fUpwind[50]+0.1767766952966368*alpha[0]*fUpwind[49]+0.1767766952966368*alpha[1]*fUpwind[48]+0.1581138830084189*alpha[22]*fUpwind[27]+0.1581138830084189*alpha[14]*fUpwind[21]+0.1581138830084189*alpha[13]*fUpwind[20]+0.1581138830084189*fUpwind[13]*alpha[20]+0.1581138830084189*alpha[5]*fUpwind[12]+0.1581138830084189*fUpwind[5]*alpha[12]; 
  Ghat[50] = 0.1767766952966368*alpha[18]*fUpwind[63]+0.1767766952966368*alpha[11]*fUpwind[62]+0.1767766952966368*alpha[9]*fUpwind[60]+0.1767766952966368*alpha[7]*fUpwind[59]+0.1767766952966368*alpha[19]*fUpwind[58]+0.1767766952966368*alpha[4]*fUpwind[57]+0.1767766952966368*alpha[3]*fUpwind[55]+0.1767766952966368*alpha[1]*fUpwind[53]+0.1767766952966368*alpha[10]*fUpwind[52]+0.1767766952966368*alpha[8]*fUpwind[51]+0.1767766952966368*alpha[0]*fUpwind[50]+0.1767766952966368*alpha[6]*fUpwind[49]+0.1767766952966368*alpha[2]*fUpwind[48]+0.1581138830084189*alpha[14]*fUpwind[22]+0.1581138830084189*fUpwind[14]*alpha[22]+0.1581138830084189*alpha[12]*fUpwind[20]+0.1581138830084189*fUpwind[12]*alpha[20]+0.1581138830084189*alpha[5]*fUpwind[13]+0.1581138830084189*fUpwind[5]*alpha[13]; 
  Ghat[51] = 0.1767766952966368*alpha[10]*fUpwind[62]+0.1767766952966368*alpha[9]*fUpwind[61]+0.1767766952966368*alpha[6]*fUpwind[59]+0.1767766952966368*alpha[4]*fUpwind[58]+0.1767766952966368*alpha[19]*fUpwind[57]+0.1767766952966368*alpha[18]*fUpwind[56]+0.1767766952966368*alpha[2]*fUpwind[55]+0.1767766952966368*alpha[1]*fUpwind[54]+0.1767766952966368*alpha[11]*fUpwind[52]+0.1767766952966368*alpha[0]*fUpwind[51]+0.1767766952966368*alpha[8]*fUpwind[50]+0.1767766952966368*alpha[7]*fUpwind[49]+0.1767766952966368*alpha[3]*fUpwind[48]+0.1581138830084189*alpha[20]*fUpwind[27]+0.1581138830084189*alpha[13]*fUpwind[22]+0.1581138830084189*fUpwind[13]*alpha[22]+0.1581138830084189*alpha[12]*fUpwind[21]+0.1581138830084189*alpha[5]*fUpwind[14]+0.1581138830084189*fUpwind[5]*alpha[14]; 
  Ghat[52] = 0.1767766952966368*alpha[8]*fUpwind[62]+0.1767766952966368*alpha[7]*fUpwind[61]+0.1767766952966368*alpha[6]*fUpwind[60]+0.1767766952966368*alpha[3]*fUpwind[58]+0.1767766952966368*alpha[2]*fUpwind[57]+0.1767766952966368*alpha[1]*fUpwind[56]+0.1767766952966368*alpha[19]*fUpwind[55]+0.1767766952966368*alpha[18]*fUpwind[54]+0.1767766952966368*alpha[0]*fUpwind[52]+0.1767766952966368*alpha[11]*fUpwind[51]+0.1767766952966368*alpha[10]*fUpwind[50]+0.1767766952966368*alpha[9]*fUpwind[49]+0.1767766952966368*alpha[4]*fUpwind[48]+0.1581138830084189*alpha[22]*fUpwind[30]+0.1581138830084189*alpha[20]*fUpwind[28]+0.1581138830084189*alpha[14]*fUpwind[25]+0.1581138830084189*alpha[13]*fUpwind[24]+0.1581138830084189*alpha[12]*fUpwind[23]+0.1581138830084189*alpha[5]*fUpwind[15]; 
  Ghat[53] = 0.1767766952966368*alpha[11]*fUpwind[63]+0.1767766952966368*alpha[18]*fUpwind[62]+0.1767766952966368*alpha[19]*fUpwind[61]+0.1767766952966368*alpha[4]*fUpwind[60]+0.1767766952966368*alpha[3]*fUpwind[59]+0.1767766952966368*alpha[9]*fUpwind[57]+0.1767766952966368*alpha[10]*fUpwind[56]+0.1767766952966368*alpha[7]*fUpwind[55]+0.1767766952966368*alpha[8]*fUpwind[54]+0.1767766952966368*alpha[0]*fUpwind[53]+0.1767766952966368*alpha[1]*fUpwind[50]+0.1767766952966368*alpha[2]*fUpwind[49]+0.1767766952966368*alpha[6]*fUpwind[48]+0.1581138830084189*alpha[14]*fUpwind[27]+0.1581138830084189*fUpwind[21]*alpha[22]+0.1581138830084189*alpha[5]*fUpwind[20]+0.1581138830084189*fUpwind[5]*alpha[20]+0.1581138830084189*alpha[12]*fUpwind[13]+0.1581138830084189*fUpwind[12]*alpha[13]; 
  Ghat[54] = 0.1767766952966368*alpha[10]*fUpwind[63]+0.1767766952966368*alpha[4]*fUpwind[61]+0.1767766952966368*alpha[19]*fUpwind[60]+0.1767766952966368*alpha[2]*fUpwind[59]+0.1767766952966368*alpha[9]*fUpwind[58]+0.1767766952966368*alpha[11]*fUpwind[56]+0.1767766952966368*alpha[6]*fUpwind[55]+0.1767766952966368*alpha[0]*fUpwind[54]+0.1767766952966368*alpha[8]*fUpwind[53]+0.1767766952966368*alpha[18]*fUpwind[52]+0.1767766952966368*alpha[1]*fUpwind[51]+0.1767766952966368*alpha[3]*fUpwind[49]+0.1767766952966368*alpha[7]*fUpwind[48]+0.1581138830084189*alpha[13]*fUpwind[27]+0.1581138830084189*alpha[20]*fUpwind[22]+0.1581138830084189*fUpwind[20]*alpha[22]+0.1581138830084189*alpha[5]*fUpwind[21]+0.1581138830084189*alpha[12]*fUpwind[14]+0.1581138830084189*fUpwind[12]*alpha[14]; 
  Ghat[55] = 0.1767766952966368*alpha[9]*fUpwind[63]+0.1767766952966368*alpha[4]*fUpwind[62]+0.1767766952966368*alpha[18]*fUpwind[60]+0.1767766952966368*alpha[1]*fUpwind[59]+0.1767766952966368*alpha[10]*fUpwind[58]+0.1767766952966368*alpha[11]*fUpwind[57]+0.1767766952966368*alpha[0]*fUpwind[55]+0.1767766952966368*alpha[6]*fUpwind[54]+0.1767766952966368*alpha[7]*fUpwind[53]+0.1767766952966368*alpha[19]*fUpwind[52]+0.1767766952966368*alpha[2]*fUpwind[51]+0.1767766952966368*alpha[3]*fUpwind[50]+0.1767766952966368*alpha[8]*fUpwind[48]+0.1581138830084189*alpha[12]*fUpwind[27]+0.1581138830084189*alpha[5]*fUpwind[22]+0.1581138830084189*fUpwind[5]*alpha[22]+0.1581138830084189*alpha[20]*fUpwind[21]+0.1581138830084189*alpha[13]*fUpwind[14]+0.1581138830084189*fUpwind[13]*alpha[14]; 
  Ghat[56] = 0.1767766952966368*alpha[8]*fUpwind[63]+0.1767766952966368*alpha[3]*fUpwind[61]+0.1767766952966368*alpha[2]*fUpwind[60]+0.1767766952966368*alpha[19]*fUpwind[59]+0.1767766952966368*alpha[7]*fUpwind[58]+0.1767766952966368*alpha[6]*fUpwind[57]+0.1767766952966368*alpha[0]*fUpwind[56]+0.1767766952966368*alpha[11]*fUpwind[54]+0.1767766952966368*alpha[10]*fUpwind[53]+0.1767766952966368*alpha[1]*fUpwind[52]+0.1767766952966368*alpha[18]*fUpwind[51]+0.1767766952966368*alpha[4]*fUpwind[49]+0.1767766952966368*alpha[9]*fUpwind[48]+0.1581138830084189*alpha[22]*fUpwind[31]+0.1581138830084189*alpha[14]*fUpwind[29]+0.1581138830084189*alpha[13]*fUpwind[28]+0.1581138830084189*alpha[20]*fUpwind[24]+0.1581138830084189*alpha[5]*fUpwind[23]+0.1581138830084189*alpha[12]*fUpwind[15]; 
  Ghat[57] = 0.1767766952966368*alpha[7]*fUpwind[63]+0.1767766952966368*alpha[3]*fUpwind[62]+0.1767766952966368*alpha[1]*fUpwind[60]+0.1767766952966368*alpha[18]*fUpwind[59]+0.1767766952966368*alpha[8]*fUpwind[58]+0.1767766952966368*alpha[0]*fUpwind[57]+0.1767766952966368*alpha[6]*fUpwind[56]+0.1767766952966368*alpha[11]*fUpwind[55]+0.1767766952966368*alpha[9]*fUpwind[53]+0.1767766952966368*alpha[2]*fUpwind[52]+0.1767766952966368*alpha[19]*fUpwind[51]+0.1767766952966368*alpha[4]*fUpwind[50]+0.1767766952966368*alpha[10]*fUpwind[48]+0.1581138830084189*alpha[14]*fUpwind[30]+0.1581138830084189*alpha[12]*fUpwind[28]+0.1581138830084189*alpha[22]*fUpwind[25]+0.1581138830084189*alpha[5]*fUpwind[24]+0.1581138830084189*alpha[20]*fUpwind[23]+0.1581138830084189*alpha[13]*fUpwind[15]; 
  Ghat[58] = 0.1767766952966368*alpha[6]*fUpwind[63]+0.1767766952966368*alpha[2]*fUpwind[62]+0.1767766952966368*alpha[1]*fUpwind[61]+0.1767766952966368*alpha[0]*fUpwind[58]+0.1767766952966368*alpha[8]*fUpwind[57]+0.1767766952966368*alpha[7]*fUpwind[56]+0.1767766952966368*alpha[10]*fUpwind[55]+0.1767766952966368*alpha[9]*fUpwind[54]+0.1767766952966368*alpha[3]*fUpwind[52]+0.1767766952966368*alpha[4]*fUpwind[51]+0.1767766952966368*alpha[19]*fUpwind[50]+0.1767766952966368*alpha[18]*fUpwind[49]+0.1767766952966368*alpha[11]*fUpwind[48]+0.1581138830084189*alpha[20]*fUpwind[31]+0.1581138830084189*alpha[13]*fUpwind[30]+0.1581138830084189*alpha[12]*fUpwind[29]+0.1581138830084189*alpha[5]*fUpwind[25]+0.1581138830084189*alpha[22]*fUpwind[24]+0.1581138830084189*alpha[14]*fUpwind[15]; 
  Ghat[59] = 0.1767766952966368*alpha[4]*fUpwind[63]+0.1767766952966368*alpha[9]*fUpwind[62]+0.1767766952966368*alpha[10]*fUpwind[61]+0.1767766952966368*alpha[11]*fUpwind[60]+0.1767766952966368*alpha[0]*fUpwind[59]+0.1767766952966368*alpha[18]*fUpwind[57]+0.1767766952966368*alpha[19]*fUpwind[56]+0.1767766952966368*alpha[1]*fUpwind[55]+0.1767766952966368*alpha[2]*fUpwind[54]+0.1767766952966368*alpha[3]*fUpwind[53]+0.1767766952966368*alpha[6]*fUpwind[51]+0.1767766952966368*alpha[7]*fUpwind[50]+0.1767766952966368*alpha[8]*fUpwind[49]+0.1581138830084189*alpha[5]*fUpwind[27]+0.1581138830084189*alpha[12]*fUpwind[22]+0.1581138830084189*fUpwind[12]*alpha[22]+0.1581138830084189*alpha[13]*fUpwind[21]+0.1581138830084189*alpha[14]*fUpwind[20]+0.1581138830084189*fUpwind[14]*alpha[20]; 
  Ghat[60] = 0.1767766952966368*alpha[3]*fUpwind[63]+0.1767766952966368*alpha[7]*fUpwind[62]+0.1767766952966368*alpha[8]*fUpwind[61]+0.1767766952966368*alpha[0]*fUpwind[60]+0.1767766952966368*alpha[11]*fUpwind[59]+0.1767766952966368*alpha[1]*fUpwind[57]+0.1767766952966368*alpha[2]*fUpwind[56]+0.1767766952966368*alpha[18]*fUpwind[55]+0.1767766952966368*alpha[19]*fUpwind[54]+0.1767766952966368*alpha[4]*fUpwind[53]+0.1767766952966368*alpha[6]*fUpwind[52]+0.1767766952966368*alpha[9]*fUpwind[50]+0.1767766952966368*alpha[10]*fUpwind[49]+0.1581138830084189*alpha[14]*fUpwind[31]+0.1581138830084189*alpha[22]*fUpwind[29]+0.1581138830084189*alpha[5]*fUpwind[28]+0.1581138830084189*alpha[12]*fUpwind[24]+0.1581138830084189*alpha[13]*fUpwind[23]+0.1581138830084189*fUpwind[15]*alpha[20]; 
  Ghat[61] = 0.1767766952966368*alpha[2]*fUpwind[63]+0.1767766952966368*alpha[6]*fUpwind[62]+0.1767766952966368*alpha[0]*fUpwind[61]+0.1767766952966368*alpha[8]*fUpwind[60]+0.1767766952966368*alpha[10]*fUpwind[59]+0.1767766952966368*alpha[1]*fUpwind[58]+0.1767766952966368*alpha[3]*fUpwind[56]+0.1767766952966368*alpha[4]*fUpwind[54]+0.1767766952966368*alpha[19]*fUpwind[53]+0.1767766952966368*alpha[7]*fUpwind[52]+0.1767766952966368*alpha[9]*fUpwind[51]+0.1767766952966368*alpha[11]*fUpwind[49]+0.1767766952966368*alpha[18]*fUpwind[48]+0.1581138830084189*alpha[13]*fUpwind[31]+0.1581138830084189*alpha[20]*fUpwind[30]+0.1581138830084189*alpha[5]*fUpwind[29]+0.1581138830084189*alpha[22]*fUpwind[28]+0.1581138830084189*alpha[12]*fUpwind[25]+0.1581138830084189*alpha[14]*fUpwind[23]; 
  Ghat[62] = 0.1767766952966368*alpha[1]*fUpwind[63]+0.1767766952966368*alpha[0]*fUpwind[62]+0.1767766952966368*alpha[6]*fUpwind[61]+0.1767766952966368*alpha[7]*fUpwind[60]+0.1767766952966368*alpha[9]*fUpwind[59]+0.1767766952966368*alpha[2]*fUpwind[58]+0.1767766952966368*alpha[3]*fUpwind[57]+0.1767766952966368*alpha[4]*fUpwind[55]+0.1767766952966368*alpha[18]*fUpwind[53]+0.1767766952966368*alpha[8]*fUpwind[52]+0.1767766952966368*alpha[10]*fUpwind[51]+0.1767766952966368*alpha[11]*fUpwind[50]+0.1767766952966368*alpha[19]*fUpwind[48]+0.1581138830084189*alpha[12]*fUpwind[31]+0.1581138830084189*alpha[5]*fUpwind[30]+0.1581138830084189*alpha[20]*fUpwind[29]+0.1581138830084189*alpha[13]*fUpwind[25]+0.1581138830084189*alpha[14]*fUpwind[24]+0.1581138830084189*fUpwind[15]*alpha[22]; 
  Ghat[63] = 0.1767766952966368*alpha[0]*fUpwind[63]+0.1767766952966368*alpha[1]*fUpwind[62]+0.1767766952966368*alpha[2]*fUpwind[61]+0.1767766952966368*alpha[3]*fUpwind[60]+0.1767766952966368*alpha[4]*fUpwind[59]+0.1767766952966368*alpha[6]*fUpwind[58]+0.1767766952966368*alpha[7]*fUpwind[57]+0.1767766952966368*alpha[8]*fUpwind[56]+0.1767766952966368*alpha[9]*fUpwind[55]+0.1767766952966368*alpha[10]*fUpwind[54]+0.1767766952966368*alpha[11]*fUpwind[53]+0.1767766952966368*alpha[18]*fUpwind[50]+0.1767766952966368*alpha[19]*fUpwind[49]+0.1581138830084189*alpha[5]*fUpwind[31]+0.1581138830084189*alpha[12]*fUpwind[30]+0.1581138830084189*alpha[13]*fUpwind[29]+0.1581138830084189*alpha[14]*fUpwind[28]+0.1581138830084189*alpha[20]*fUpwind[25]+0.1581138830084189*alpha[22]*fUpwind[23]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv10; 
  out[1] += 0.7071067811865475*Ghat[1]*dv10; 
  out[2] += 0.7071067811865475*Ghat[2]*dv10; 
  out[3] += 0.7071067811865475*Ghat[3]*dv10; 
  out[4] += -1.224744871391589*Ghat[0]*dv10; 
  out[5] += 0.7071067811865475*Ghat[4]*dv10; 
  out[6] += 0.7071067811865475*Ghat[5]*dv10; 
  out[7] += 0.7071067811865475*Ghat[6]*dv10; 
  out[8] += 0.7071067811865475*Ghat[7]*dv10; 
  out[9] += 0.7071067811865475*Ghat[8]*dv10; 
  out[10] += -1.224744871391589*Ghat[1]*dv10; 
  out[11] += -1.224744871391589*Ghat[2]*dv10; 
  out[12] += -1.224744871391589*Ghat[3]*dv10; 
  out[13] += 0.7071067811865475*Ghat[9]*dv10; 
  out[14] += 0.7071067811865475*Ghat[10]*dv10; 
  out[15] += 0.7071067811865475*Ghat[11]*dv10; 
  out[16] += -1.224744871391589*Ghat[4]*dv10; 
  out[17] += 0.7071067811865475*Ghat[12]*dv10; 
  out[18] += 0.7071067811865475*Ghat[13]*dv10; 
  out[19] += 0.7071067811865475*Ghat[14]*dv10; 
  out[20] += -1.224744871391589*Ghat[5]*dv10; 
  out[21] += 0.7071067811865475*Ghat[15]*dv10; 
  out[22] += 0.7071067811865475*Ghat[16]*dv10; 
  out[23] += -1.224744871391589*Ghat[6]*dv10; 
  out[24] += -1.224744871391589*Ghat[7]*dv10; 
  out[25] += -1.224744871391589*Ghat[8]*dv10; 
  out[26] += 0.7071067811865475*Ghat[17]*dv10; 
  out[27] += 0.7071067811865475*Ghat[18]*dv10; 
  out[28] += 0.7071067811865475*Ghat[19]*dv10; 
  out[29] += -1.224744871391589*Ghat[9]*dv10; 
  out[30] += -1.224744871391589*Ghat[10]*dv10; 
  out[31] += -1.224744871391589*Ghat[11]*dv10; 
  out[32] += 0.7071067811865475*Ghat[20]*dv10; 
  out[33] += 0.7071067811865475*Ghat[21]*dv10; 
  out[34] += 0.7071067811865475*Ghat[22]*dv10; 
  out[35] += -1.224744871391589*Ghat[12]*dv10; 
  out[36] += -1.224744871391589*Ghat[13]*dv10; 
  out[37] += -1.224744871391589*Ghat[14]*dv10; 
  out[38] += 0.7071067811865475*Ghat[23]*dv10; 
  out[39] += 0.7071067811865475*Ghat[24]*dv10; 
  out[40] += 0.7071067811865475*Ghat[25]*dv10; 
  out[41] += -1.224744871391589*Ghat[15]*dv10; 
  out[42] += -1.224744871391589*Ghat[16]*dv10; 
  out[43] += 0.7071067811865475*Ghat[26]*dv10; 
  out[44] += -1.224744871391589*Ghat[17]*dv10; 
  out[45] += -1.224744871391589*Ghat[18]*dv10; 
  out[46] += -1.224744871391589*Ghat[19]*dv10; 
  out[47] += 0.7071067811865475*Ghat[27]*dv10; 
  out[48] += -1.224744871391589*Ghat[20]*dv10; 
  out[49] += -1.224744871391589*Ghat[21]*dv10; 
  out[50] += -1.224744871391589*Ghat[22]*dv10; 
  out[51] += 0.7071067811865475*Ghat[28]*dv10; 
  out[52] += 0.7071067811865475*Ghat[29]*dv10; 
  out[53] += 0.7071067811865475*Ghat[30]*dv10; 
  out[54] += -1.224744871391589*Ghat[23]*dv10; 
  out[55] += -1.224744871391589*Ghat[24]*dv10; 
  out[56] += -1.224744871391589*Ghat[25]*dv10; 
  out[57] += -1.224744871391589*Ghat[26]*dv10; 
  out[58] += -1.224744871391589*Ghat[27]*dv10; 
  out[59] += 0.7071067811865475*Ghat[31]*dv10; 
  out[60] += -1.224744871391589*Ghat[28]*dv10; 
  out[61] += -1.224744871391589*Ghat[29]*dv10; 
  out[62] += -1.224744871391589*Ghat[30]*dv10; 
  out[63] += -1.224744871391589*Ghat[31]*dv10; 
  out[64] += 1.58113883008419*Ghat[0]*dv10; 
  out[65] += 1.58113883008419*Ghat[1]*dv10; 
  out[66] += 1.58113883008419*Ghat[2]*dv10; 
  out[67] += 1.58113883008419*Ghat[3]*dv10; 
  out[68] += 1.58113883008419*Ghat[4]*dv10; 
  out[69] += 1.58113883008419*Ghat[5]*dv10; 
  out[70] += 1.58113883008419*Ghat[6]*dv10; 
  out[71] += 1.58113883008419*Ghat[7]*dv10; 
  out[72] += 1.58113883008419*Ghat[8]*dv10; 
  out[73] += 1.58113883008419*Ghat[9]*dv10; 
  out[74] += 1.58113883008419*Ghat[10]*dv10; 
  out[75] += 1.58113883008419*Ghat[11]*dv10; 
  out[76] += 1.58113883008419*Ghat[12]*dv10; 
  out[77] += 1.58113883008419*Ghat[13]*dv10; 
  out[78] += 1.58113883008419*Ghat[14]*dv10; 
  out[79] += 1.58113883008419*Ghat[15]*dv10; 
  out[80] += 1.58113883008419*Ghat[16]*dv10; 
  out[81] += 1.58113883008419*Ghat[17]*dv10; 
  out[82] += 1.58113883008419*Ghat[18]*dv10; 
  out[83] += 1.58113883008419*Ghat[19]*dv10; 
  out[84] += 1.58113883008419*Ghat[20]*dv10; 
  out[85] += 1.58113883008419*Ghat[21]*dv10; 
  out[86] += 1.58113883008419*Ghat[22]*dv10; 
  out[87] += 1.58113883008419*Ghat[23]*dv10; 
  out[88] += 1.58113883008419*Ghat[24]*dv10; 
  out[89] += 1.58113883008419*Ghat[25]*dv10; 
  out[90] += 1.58113883008419*Ghat[26]*dv10; 
  out[91] += 1.58113883008419*Ghat[27]*dv10; 
  out[92] += 1.58113883008419*Ghat[28]*dv10; 
  out[93] += 1.58113883008419*Ghat[29]*dv10; 
  out[94] += 1.58113883008419*Ghat[30]*dv10; 
  out[95] += 1.58113883008419*Ghat[31]*dv10; 
  out[96] += 0.7071067811865475*Ghat[32]*dv10; 
  out[97] += 0.7071067811865475*Ghat[33]*dv10; 
  out[98] += 0.7071067811865475*Ghat[34]*dv10; 
  out[99] += 0.7071067811865475*Ghat[35]*dv10; 
  out[100] += -1.224744871391589*Ghat[32]*dv10; 
  out[101] += 0.7071067811865475*Ghat[36]*dv10; 
  out[102] += 0.7071067811865475*Ghat[37]*dv10; 
  out[103] += 0.7071067811865475*Ghat[38]*dv10; 
  out[104] += 0.7071067811865475*Ghat[39]*dv10; 
  out[105] += -1.224744871391589*Ghat[33]*dv10; 
  out[106] += -1.224744871391589*Ghat[34]*dv10; 
  out[107] += -1.224744871391589*Ghat[35]*dv10; 
  out[108] += 0.7071067811865475*Ghat[40]*dv10; 
  out[109] += 0.7071067811865475*Ghat[41]*dv10; 
  out[110] += 0.7071067811865475*Ghat[42]*dv10; 
  out[111] += -1.224744871391589*Ghat[36]*dv10; 
  out[112] += 0.7071067811865475*Ghat[43]*dv10; 
  out[113] += -1.224744871391589*Ghat[37]*dv10; 
  out[114] += -1.224744871391589*Ghat[38]*dv10; 
  out[115] += -1.224744871391589*Ghat[39]*dv10; 
  out[116] += 0.7071067811865475*Ghat[44]*dv10; 
  out[117] += 0.7071067811865475*Ghat[45]*dv10; 
  out[118] += 0.7071067811865475*Ghat[46]*dv10; 
  out[119] += -1.224744871391589*Ghat[40]*dv10; 
  out[120] += -1.224744871391589*Ghat[41]*dv10; 
  out[121] += -1.224744871391589*Ghat[42]*dv10; 
  out[122] += -1.224744871391589*Ghat[43]*dv10; 
  out[123] += 0.7071067811865475*Ghat[47]*dv10; 
  out[124] += -1.224744871391589*Ghat[44]*dv10; 
  out[125] += -1.224744871391589*Ghat[45]*dv10; 
  out[126] += -1.224744871391589*Ghat[46]*dv10; 
  out[127] += -1.224744871391589*Ghat[47]*dv10; 
  out[128] += 0.7071067811865475*Ghat[48]*dv10; 
  out[129] += 0.7071067811865475*Ghat[49]*dv10; 
  out[130] += 0.7071067811865475*Ghat[50]*dv10; 
  out[131] += 0.7071067811865475*Ghat[51]*dv10; 
  out[132] += -1.224744871391589*Ghat[48]*dv10; 
  out[133] += 0.7071067811865475*Ghat[52]*dv10; 
  out[134] += 0.7071067811865475*Ghat[53]*dv10; 
  out[135] += 0.7071067811865475*Ghat[54]*dv10; 
  out[136] += 0.7071067811865475*Ghat[55]*dv10; 
  out[137] += -1.224744871391589*Ghat[49]*dv10; 
  out[138] += -1.224744871391589*Ghat[50]*dv10; 
  out[139] += -1.224744871391589*Ghat[51]*dv10; 
  out[140] += 0.7071067811865475*Ghat[56]*dv10; 
  out[141] += 0.7071067811865475*Ghat[57]*dv10; 
  out[142] += 0.7071067811865475*Ghat[58]*dv10; 
  out[143] += -1.224744871391589*Ghat[52]*dv10; 
  out[144] += 0.7071067811865475*Ghat[59]*dv10; 
  out[145] += -1.224744871391589*Ghat[53]*dv10; 
  out[146] += -1.224744871391589*Ghat[54]*dv10; 
  out[147] += -1.224744871391589*Ghat[55]*dv10; 
  out[148] += 0.7071067811865475*Ghat[60]*dv10; 
  out[149] += 0.7071067811865475*Ghat[61]*dv10; 
  out[150] += 0.7071067811865475*Ghat[62]*dv10; 
  out[151] += -1.224744871391589*Ghat[56]*dv10; 
  out[152] += -1.224744871391589*Ghat[57]*dv10; 
  out[153] += -1.224744871391589*Ghat[58]*dv10; 
  out[154] += -1.224744871391589*Ghat[59]*dv10; 
  out[155] += 0.7071067811865475*Ghat[63]*dv10; 
  out[156] += -1.224744871391589*Ghat[60]*dv10; 
  out[157] += -1.224744871391589*Ghat[61]*dv10; 
  out[158] += -1.224744871391589*Ghat[62]*dv10; 
  out[159] += -1.224744871391589*Ghat[63]*dv10; 

  } 
  return 0.;

} 
