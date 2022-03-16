#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_4x_p1_surfx4_eval_quad.h> 
#include <gkyl_basis_ser_4x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_poisson_boundary_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double *fac_phi, const double *vecA, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // fac_phi:     potential (scaled by appropriate factors).
  // vecA:        vector potential (scaled by appropriate factors). Unused in pure Vlasov-Poisson. 
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv11 = 2/dxv[3]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double *phi = &fac_phi[0]; 
  const double dx10 = 2/dxv[0]; 
  const double dx11 = 2/dxv[1]; 
  double alpha[8] = {0.0}; 

  alpha[0] = -2.449489742783178*phi[2]*dx11; 
  alpha[1] = -2.449489742783178*phi[3]*dx11; 

  double fUpwindQuad[8] = {0.0};
  double fUpwind[8] = {0.0};
  double Ghat[8] = {0.0}; 

  if (edge == -1) { 

  if (alpha[0]-alpha[1] > 0) { 
    fUpwindQuad[0] = ser_4x_p1_surfx4_eval_quad_node_0_r(fSkin); 
    fUpwindQuad[1] = ser_4x_p1_surfx4_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_4x_p1_surfx4_eval_quad_node_0_l(fEdge); 
    fUpwindQuad[1] = ser_4x_p1_surfx4_eval_quad_node_1_l(fEdge); 
  } 
  if (alpha[0]-alpha[1] > 0) { 
    fUpwindQuad[2] = ser_4x_p1_surfx4_eval_quad_node_2_r(fSkin); 
    fUpwindQuad[3] = ser_4x_p1_surfx4_eval_quad_node_3_r(fSkin); 
  } else { 
    fUpwindQuad[2] = ser_4x_p1_surfx4_eval_quad_node_2_l(fEdge); 
    fUpwindQuad[3] = ser_4x_p1_surfx4_eval_quad_node_3_l(fEdge); 
  } 
  if (alpha[0]-alpha[1] > 0) { 
    fUpwindQuad[4] = ser_4x_p1_surfx4_eval_quad_node_4_r(fSkin); 
    fUpwindQuad[5] = ser_4x_p1_surfx4_eval_quad_node_5_r(fSkin); 
  } else { 
    fUpwindQuad[4] = ser_4x_p1_surfx4_eval_quad_node_4_l(fEdge); 
    fUpwindQuad[5] = ser_4x_p1_surfx4_eval_quad_node_5_l(fEdge); 
  } 
  if (alpha[0]-alpha[1] > 0) { 
    fUpwindQuad[6] = ser_4x_p1_surfx4_eval_quad_node_6_r(fSkin); 
    fUpwindQuad[7] = ser_4x_p1_surfx4_eval_quad_node_7_r(fSkin); 
  } else { 
    fUpwindQuad[6] = ser_4x_p1_surfx4_eval_quad_node_6_l(fEdge); 
    fUpwindQuad[7] = ser_4x_p1_surfx4_eval_quad_node_7_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_4x_p1_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.3535533905932737*(alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] = 0.3535533905932737*(alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] = 0.3535533905932737*(alpha[1]*fUpwind[4]+alpha[0]*fUpwind[2]); 
  Ghat[3] = 0.3535533905932737*(alpha[1]*fUpwind[5]+alpha[0]*fUpwind[3]); 
  Ghat[4] = 0.3535533905932737*(alpha[0]*fUpwind[4]+alpha[1]*fUpwind[2]); 
  Ghat[5] = 0.3535533905932737*(alpha[0]*fUpwind[5]+alpha[1]*fUpwind[3]); 
  Ghat[6] = 0.3535533905932737*(alpha[1]*fUpwind[7]+alpha[0]*fUpwind[6]); 
  Ghat[7] = 0.3535533905932737*(alpha[0]*fUpwind[7]+alpha[1]*fUpwind[6]); 

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

  } else { 

  if (alpha[0]-alpha[1] > 0) { 
    fUpwindQuad[0] = ser_4x_p1_surfx4_eval_quad_node_0_r(fEdge); 
    fUpwindQuad[1] = ser_4x_p1_surfx4_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_4x_p1_surfx4_eval_quad_node_0_l(fSkin); 
    fUpwindQuad[1] = ser_4x_p1_surfx4_eval_quad_node_1_l(fSkin); 
  } 
  if (alpha[0]-alpha[1] > 0) { 
    fUpwindQuad[2] = ser_4x_p1_surfx4_eval_quad_node_2_r(fEdge); 
    fUpwindQuad[3] = ser_4x_p1_surfx4_eval_quad_node_3_r(fEdge); 
  } else { 
    fUpwindQuad[2] = ser_4x_p1_surfx4_eval_quad_node_2_l(fSkin); 
    fUpwindQuad[3] = ser_4x_p1_surfx4_eval_quad_node_3_l(fSkin); 
  } 
  if (alpha[0]-alpha[1] > 0) { 
    fUpwindQuad[4] = ser_4x_p1_surfx4_eval_quad_node_4_r(fEdge); 
    fUpwindQuad[5] = ser_4x_p1_surfx4_eval_quad_node_5_r(fEdge); 
  } else { 
    fUpwindQuad[4] = ser_4x_p1_surfx4_eval_quad_node_4_l(fSkin); 
    fUpwindQuad[5] = ser_4x_p1_surfx4_eval_quad_node_5_l(fSkin); 
  } 
  if (alpha[0]-alpha[1] > 0) { 
    fUpwindQuad[6] = ser_4x_p1_surfx4_eval_quad_node_6_r(fEdge); 
    fUpwindQuad[7] = ser_4x_p1_surfx4_eval_quad_node_7_r(fEdge); 
  } else { 
    fUpwindQuad[6] = ser_4x_p1_surfx4_eval_quad_node_6_l(fSkin); 
    fUpwindQuad[7] = ser_4x_p1_surfx4_eval_quad_node_7_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_4x_p1_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.3535533905932737*(alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] = 0.3535533905932737*(alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] = 0.3535533905932737*(alpha[1]*fUpwind[4]+alpha[0]*fUpwind[2]); 
  Ghat[3] = 0.3535533905932737*(alpha[1]*fUpwind[5]+alpha[0]*fUpwind[3]); 
  Ghat[4] = 0.3535533905932737*(alpha[0]*fUpwind[4]+alpha[1]*fUpwind[2]); 
  Ghat[5] = 0.3535533905932737*(alpha[0]*fUpwind[5]+alpha[1]*fUpwind[3]); 
  Ghat[6] = 0.3535533905932737*(alpha[1]*fUpwind[7]+alpha[0]*fUpwind[6]); 
  Ghat[7] = 0.3535533905932737*(alpha[0]*fUpwind[7]+alpha[1]*fUpwind[6]); 

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

  } 
} 
