#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_hyb_2x2v_p1_surfx4_eval_quad.h> 
#include <gkyl_basis_hyb_2x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_poisson_boundary_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // field:       potential (scaled by appropriate factors).
  // ext_field:   vector potential (scaled by appropriate factors). Unused in pure Vlasov-Poisson. 
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv11 = 2/dxv[3]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double *phi = &field[0]; 
  const double dx10 = 2/dxv[0]; 
  const double dx11 = 2/dxv[1]; 
  double alpha[12] = {0.0}; 

  alpha[0] = -2.449489742783178*phi[2]*dx11; 
  alpha[1] = -2.449489742783178*phi[3]*dx11; 

  double fUpwindQuad[12] = {0.0};
  double fUpwind[12] = {0.0};
  double Ghat[12] = {0.0}; 

  if (edge == -1) { 

  if (0.3535533905932737*alpha[0]-0.3535533905932737*alpha[1] > 0) { 
    fUpwindQuad[0] = hyb_2x2v_p1_surfx4_eval_quad_node_0_r(fSkin); 
    fUpwindQuad[1] = hyb_2x2v_p1_surfx4_eval_quad_node_1_r(fSkin); 
    fUpwindQuad[2] = hyb_2x2v_p1_surfx4_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[0] = hyb_2x2v_p1_surfx4_eval_quad_node_0_l(fEdge); 
    fUpwindQuad[1] = hyb_2x2v_p1_surfx4_eval_quad_node_1_l(fEdge); 
    fUpwindQuad[2] = hyb_2x2v_p1_surfx4_eval_quad_node_2_l(fEdge); 
  } 
  if (0.3535533905932737*alpha[0]-0.3535533905932737*alpha[1] > 0) { 
    fUpwindQuad[3] = hyb_2x2v_p1_surfx4_eval_quad_node_3_r(fSkin); 
    fUpwindQuad[4] = hyb_2x2v_p1_surfx4_eval_quad_node_4_r(fSkin); 
    fUpwindQuad[5] = hyb_2x2v_p1_surfx4_eval_quad_node_5_r(fSkin); 
  } else { 
    fUpwindQuad[3] = hyb_2x2v_p1_surfx4_eval_quad_node_3_l(fEdge); 
    fUpwindQuad[4] = hyb_2x2v_p1_surfx4_eval_quad_node_4_l(fEdge); 
    fUpwindQuad[5] = hyb_2x2v_p1_surfx4_eval_quad_node_5_l(fEdge); 
  } 
  if (0.3535533905932737*alpha[0]-0.3535533905932737*alpha[1] > 0) { 
    fUpwindQuad[6] = hyb_2x2v_p1_surfx4_eval_quad_node_6_r(fSkin); 
    fUpwindQuad[7] = hyb_2x2v_p1_surfx4_eval_quad_node_7_r(fSkin); 
    fUpwindQuad[8] = hyb_2x2v_p1_surfx4_eval_quad_node_8_r(fSkin); 
  } else { 
    fUpwindQuad[6] = hyb_2x2v_p1_surfx4_eval_quad_node_6_l(fEdge); 
    fUpwindQuad[7] = hyb_2x2v_p1_surfx4_eval_quad_node_7_l(fEdge); 
    fUpwindQuad[8] = hyb_2x2v_p1_surfx4_eval_quad_node_8_l(fEdge); 
  } 
  if (0.3535533905932737*alpha[0]-0.3535533905932737*alpha[1] > 0) { 
    fUpwindQuad[9] = hyb_2x2v_p1_surfx4_eval_quad_node_9_r(fSkin); 
    fUpwindQuad[10] = hyb_2x2v_p1_surfx4_eval_quad_node_10_r(fSkin); 
    fUpwindQuad[11] = hyb_2x2v_p1_surfx4_eval_quad_node_11_r(fSkin); 
  } else { 
    fUpwindQuad[9] = hyb_2x2v_p1_surfx4_eval_quad_node_9_l(fEdge); 
    fUpwindQuad[10] = hyb_2x2v_p1_surfx4_eval_quad_node_10_l(fEdge); 
    fUpwindQuad[11] = hyb_2x2v_p1_surfx4_eval_quad_node_11_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_2x2v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.3535533905932737*alpha[1]*fUpwind[1]+0.3535533905932737*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.3535533905932737*alpha[0]*fUpwind[1]+0.3535533905932737*fUpwind[0]*alpha[1]; 
  Ghat[2] = 0.3535533905932737*alpha[1]*fUpwind[4]+0.3535533905932737*alpha[0]*fUpwind[2]; 
  Ghat[3] = 0.3535533905932737*alpha[1]*fUpwind[5]+0.3535533905932737*alpha[0]*fUpwind[3]; 
  Ghat[4] = 0.3535533905932737*alpha[0]*fUpwind[4]+0.3535533905932737*alpha[1]*fUpwind[2]; 
  Ghat[5] = 0.3535533905932737*alpha[0]*fUpwind[5]+0.3535533905932737*alpha[1]*fUpwind[3]; 
  Ghat[6] = 0.3535533905932737*alpha[1]*fUpwind[7]+0.3535533905932737*alpha[0]*fUpwind[6]; 
  Ghat[7] = 0.3535533905932737*alpha[0]*fUpwind[7]+0.3535533905932737*alpha[1]*fUpwind[6]; 
  Ghat[8] = 0.3535533905932737*alpha[1]*fUpwind[9]+0.3535533905932737*alpha[0]*fUpwind[8]; 
  Ghat[9] = 0.3535533905932737*alpha[0]*fUpwind[9]+0.3535533905932737*alpha[1]*fUpwind[8]; 
  Ghat[10] = 0.3535533905932737*alpha[1]*fUpwind[11]+0.3535533905932737*alpha[0]*fUpwind[10]; 
  Ghat[11] = 0.3535533905932737*alpha[0]*fUpwind[11]+0.3535533905932737*alpha[1]*fUpwind[10]; 

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
  out[12] += -1.224744871391589*Ghat[4]*dv11; 
  out[13] += -1.224744871391589*Ghat[5]*dv11; 
  out[14] += -1.224744871391589*Ghat[6]*dv11; 
  out[15] += -1.224744871391589*Ghat[7]*dv11; 
  out[16] += -0.7071067811865475*Ghat[8]*dv11; 
  out[17] += -0.7071067811865475*Ghat[9]*dv11; 
  out[18] += -0.7071067811865475*Ghat[10]*dv11; 
  out[19] += -1.224744871391589*Ghat[8]*dv11; 
  out[20] += -0.7071067811865475*Ghat[11]*dv11; 
  out[21] += -1.224744871391589*Ghat[9]*dv11; 
  out[22] += -1.224744871391589*Ghat[10]*dv11; 
  out[23] += -1.224744871391589*Ghat[11]*dv11; 
  out[24] += -1.58113883008419*Ghat[0]*dv11; 
  out[25] += -1.58113883008419*Ghat[1]*dv11; 
  out[26] += -1.58113883008419*Ghat[2]*dv11; 
  out[27] += -1.58113883008419*Ghat[3]*dv11; 
  out[28] += -1.58113883008419*Ghat[4]*dv11; 
  out[29] += -1.58113883008419*Ghat[5]*dv11; 
  out[30] += -1.58113883008419*Ghat[6]*dv11; 
  out[31] += -1.58113883008419*Ghat[7]*dv11; 

  } else { 

  if (0.3535533905932737*alpha[0]-0.3535533905932737*alpha[1] > 0) { 
    fUpwindQuad[0] = hyb_2x2v_p1_surfx4_eval_quad_node_0_r(fEdge); 
    fUpwindQuad[1] = hyb_2x2v_p1_surfx4_eval_quad_node_1_r(fEdge); 
    fUpwindQuad[2] = hyb_2x2v_p1_surfx4_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[0] = hyb_2x2v_p1_surfx4_eval_quad_node_0_l(fSkin); 
    fUpwindQuad[1] = hyb_2x2v_p1_surfx4_eval_quad_node_1_l(fSkin); 
    fUpwindQuad[2] = hyb_2x2v_p1_surfx4_eval_quad_node_2_l(fSkin); 
  } 
  if (0.3535533905932737*alpha[0]-0.3535533905932737*alpha[1] > 0) { 
    fUpwindQuad[3] = hyb_2x2v_p1_surfx4_eval_quad_node_3_r(fEdge); 
    fUpwindQuad[4] = hyb_2x2v_p1_surfx4_eval_quad_node_4_r(fEdge); 
    fUpwindQuad[5] = hyb_2x2v_p1_surfx4_eval_quad_node_5_r(fEdge); 
  } else { 
    fUpwindQuad[3] = hyb_2x2v_p1_surfx4_eval_quad_node_3_l(fSkin); 
    fUpwindQuad[4] = hyb_2x2v_p1_surfx4_eval_quad_node_4_l(fSkin); 
    fUpwindQuad[5] = hyb_2x2v_p1_surfx4_eval_quad_node_5_l(fSkin); 
  } 
  if (0.3535533905932737*alpha[0]-0.3535533905932737*alpha[1] > 0) { 
    fUpwindQuad[6] = hyb_2x2v_p1_surfx4_eval_quad_node_6_r(fEdge); 
    fUpwindQuad[7] = hyb_2x2v_p1_surfx4_eval_quad_node_7_r(fEdge); 
    fUpwindQuad[8] = hyb_2x2v_p1_surfx4_eval_quad_node_8_r(fEdge); 
  } else { 
    fUpwindQuad[6] = hyb_2x2v_p1_surfx4_eval_quad_node_6_l(fSkin); 
    fUpwindQuad[7] = hyb_2x2v_p1_surfx4_eval_quad_node_7_l(fSkin); 
    fUpwindQuad[8] = hyb_2x2v_p1_surfx4_eval_quad_node_8_l(fSkin); 
  } 
  if (0.3535533905932737*alpha[0]-0.3535533905932737*alpha[1] > 0) { 
    fUpwindQuad[9] = hyb_2x2v_p1_surfx4_eval_quad_node_9_r(fEdge); 
    fUpwindQuad[10] = hyb_2x2v_p1_surfx4_eval_quad_node_10_r(fEdge); 
    fUpwindQuad[11] = hyb_2x2v_p1_surfx4_eval_quad_node_11_r(fEdge); 
  } else { 
    fUpwindQuad[9] = hyb_2x2v_p1_surfx4_eval_quad_node_9_l(fSkin); 
    fUpwindQuad[10] = hyb_2x2v_p1_surfx4_eval_quad_node_10_l(fSkin); 
    fUpwindQuad[11] = hyb_2x2v_p1_surfx4_eval_quad_node_11_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_2x2v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.3535533905932737*alpha[1]*fUpwind[1]+0.3535533905932737*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.3535533905932737*alpha[0]*fUpwind[1]+0.3535533905932737*fUpwind[0]*alpha[1]; 
  Ghat[2] = 0.3535533905932737*alpha[1]*fUpwind[4]+0.3535533905932737*alpha[0]*fUpwind[2]; 
  Ghat[3] = 0.3535533905932737*alpha[1]*fUpwind[5]+0.3535533905932737*alpha[0]*fUpwind[3]; 
  Ghat[4] = 0.3535533905932737*alpha[0]*fUpwind[4]+0.3535533905932737*alpha[1]*fUpwind[2]; 
  Ghat[5] = 0.3535533905932737*alpha[0]*fUpwind[5]+0.3535533905932737*alpha[1]*fUpwind[3]; 
  Ghat[6] = 0.3535533905932737*alpha[1]*fUpwind[7]+0.3535533905932737*alpha[0]*fUpwind[6]; 
  Ghat[7] = 0.3535533905932737*alpha[0]*fUpwind[7]+0.3535533905932737*alpha[1]*fUpwind[6]; 
  Ghat[8] = 0.3535533905932737*alpha[1]*fUpwind[9]+0.3535533905932737*alpha[0]*fUpwind[8]; 
  Ghat[9] = 0.3535533905932737*alpha[0]*fUpwind[9]+0.3535533905932737*alpha[1]*fUpwind[8]; 
  Ghat[10] = 0.3535533905932737*alpha[1]*fUpwind[11]+0.3535533905932737*alpha[0]*fUpwind[10]; 
  Ghat[11] = 0.3535533905932737*alpha[0]*fUpwind[11]+0.3535533905932737*alpha[1]*fUpwind[10]; 

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
  out[12] += -1.224744871391589*Ghat[4]*dv11; 
  out[13] += -1.224744871391589*Ghat[5]*dv11; 
  out[14] += -1.224744871391589*Ghat[6]*dv11; 
  out[15] += -1.224744871391589*Ghat[7]*dv11; 
  out[16] += 0.7071067811865475*Ghat[8]*dv11; 
  out[17] += 0.7071067811865475*Ghat[9]*dv11; 
  out[18] += 0.7071067811865475*Ghat[10]*dv11; 
  out[19] += -1.224744871391589*Ghat[8]*dv11; 
  out[20] += 0.7071067811865475*Ghat[11]*dv11; 
  out[21] += -1.224744871391589*Ghat[9]*dv11; 
  out[22] += -1.224744871391589*Ghat[10]*dv11; 
  out[23] += -1.224744871391589*Ghat[11]*dv11; 
  out[24] += 1.58113883008419*Ghat[0]*dv11; 
  out[25] += 1.58113883008419*Ghat[1]*dv11; 
  out[26] += 1.58113883008419*Ghat[2]*dv11; 
  out[27] += 1.58113883008419*Ghat[3]*dv11; 
  out[28] += 1.58113883008419*Ghat[4]*dv11; 
  out[29] += 1.58113883008419*Ghat[5]*dv11; 
  out[30] += 1.58113883008419*Ghat[6]*dv11; 
  out[31] += 1.58113883008419*Ghat[7]*dv11; 

  } 
  return 0.;

} 
