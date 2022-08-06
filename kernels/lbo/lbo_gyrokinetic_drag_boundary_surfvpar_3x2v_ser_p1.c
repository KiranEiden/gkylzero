#include <gkyl_lbo_gyrokinetic_kernels.h> 
#include <gkyl_basis_gkhyb_3x2v_p1_surfx4_eval_quad.h> 
#include <gkyl_basis_gkhyb_3x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void lbo_gyrokinetic_drag_boundary_surfvpar_3x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[5]:     cell-center coordinates. 
  // dxv[5]:   cell spacing. 
  // m_:        species mass.
  // bmag_inv:  1/(magnetic field magnitude). 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum[16]:sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[3]; 


  double alphaDrSurf[16] = {0.0}; 
  double fUpwindQuad[16] = {0.0};
  double fUpwind[16] = {0.0};;
  double Ghat[16] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*(2.828427124746191*w[3]+1.414213562373095*dxv[3])-2.828427124746191*nuUSum[0]); 
  alphaDrSurf[1] = 0.5*(nuSum[1]*(2.828427124746191*w[3]+1.414213562373095*dxv[3])-2.828427124746191*nuUSum[1]); 
  alphaDrSurf[2] = 0.5*(nuSum[2]*(2.828427124746191*w[3]+1.414213562373095*dxv[3])-2.828427124746191*nuUSum[2]); 
  alphaDrSurf[3] = 0.5*(2.828427124746191*nuSum[3]*w[3]-2.828427124746191*nuUSum[3]+1.414213562373095*dxv[3]*nuSum[3]); 
  alphaDrSurf[5] = -0.5*(2.828427124746191*nuUSum[4]+((-2.828427124746191*w[3])-1.414213562373095*dxv[3])*nuSum[4]); 
  alphaDrSurf[6] = -0.5*(2.828427124746191*nuUSum[5]+((-2.828427124746191*w[3])-1.414213562373095*dxv[3])*nuSum[5]); 
  alphaDrSurf[7] = -0.5*(2.828427124746191*nuUSum[6]+((-2.828427124746191*w[3])-1.414213562373095*dxv[3])*nuSum[6]); 
  alphaDrSurf[11] = -0.5*(2.828427124746191*nuUSum[7]+((-2.828427124746191*w[3])-1.414213562373095*dxv[3])*nuSum[7]); 

  if ((-alphaDrSurf[11])+alphaDrSurf[7]+alphaDrSurf[6]+alphaDrSurf[5]-alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = gkhyb_3x2v_p1_surfx4_eval_quad_node_0_r(fSkin); 
    fUpwindQuad[1] = gkhyb_3x2v_p1_surfx4_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[0] = gkhyb_3x2v_p1_surfx4_eval_quad_node_0_l(fEdge); 
    fUpwindQuad[1] = gkhyb_3x2v_p1_surfx4_eval_quad_node_1_l(fEdge); 
  } 
  if ((-alphaDrSurf[11])+alphaDrSurf[7]+alphaDrSurf[6]+alphaDrSurf[5]-alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[2] = gkhyb_3x2v_p1_surfx4_eval_quad_node_2_r(fSkin); 
    fUpwindQuad[3] = gkhyb_3x2v_p1_surfx4_eval_quad_node_3_r(fSkin); 
  } else { 
    fUpwindQuad[2] = gkhyb_3x2v_p1_surfx4_eval_quad_node_2_l(fEdge); 
    fUpwindQuad[3] = gkhyb_3x2v_p1_surfx4_eval_quad_node_3_l(fEdge); 
  } 
  if (alphaDrSurf[11]-alphaDrSurf[7]-alphaDrSurf[6]+alphaDrSurf[5]+alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[4] = gkhyb_3x2v_p1_surfx4_eval_quad_node_4_r(fSkin); 
    fUpwindQuad[5] = gkhyb_3x2v_p1_surfx4_eval_quad_node_5_r(fSkin); 
  } else { 
    fUpwindQuad[4] = gkhyb_3x2v_p1_surfx4_eval_quad_node_4_l(fEdge); 
    fUpwindQuad[5] = gkhyb_3x2v_p1_surfx4_eval_quad_node_5_l(fEdge); 
  } 
  if (alphaDrSurf[11]-alphaDrSurf[7]-alphaDrSurf[6]+alphaDrSurf[5]+alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[6] = gkhyb_3x2v_p1_surfx4_eval_quad_node_6_r(fSkin); 
    fUpwindQuad[7] = gkhyb_3x2v_p1_surfx4_eval_quad_node_7_r(fSkin); 
  } else { 
    fUpwindQuad[6] = gkhyb_3x2v_p1_surfx4_eval_quad_node_6_l(fEdge); 
    fUpwindQuad[7] = gkhyb_3x2v_p1_surfx4_eval_quad_node_7_l(fEdge); 
  } 
  if (alphaDrSurf[11]-alphaDrSurf[7]+alphaDrSurf[6]-alphaDrSurf[5]-alphaDrSurf[3]+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[8] = gkhyb_3x2v_p1_surfx4_eval_quad_node_8_r(fSkin); 
    fUpwindQuad[9] = gkhyb_3x2v_p1_surfx4_eval_quad_node_9_r(fSkin); 
  } else { 
    fUpwindQuad[8] = gkhyb_3x2v_p1_surfx4_eval_quad_node_8_l(fEdge); 
    fUpwindQuad[9] = gkhyb_3x2v_p1_surfx4_eval_quad_node_9_l(fEdge); 
  } 
  if (alphaDrSurf[11]-alphaDrSurf[7]+alphaDrSurf[6]-alphaDrSurf[5]-alphaDrSurf[3]+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[10] = gkhyb_3x2v_p1_surfx4_eval_quad_node_10_r(fSkin); 
    fUpwindQuad[11] = gkhyb_3x2v_p1_surfx4_eval_quad_node_11_r(fSkin); 
  } else { 
    fUpwindQuad[10] = gkhyb_3x2v_p1_surfx4_eval_quad_node_10_l(fEdge); 
    fUpwindQuad[11] = gkhyb_3x2v_p1_surfx4_eval_quad_node_11_l(fEdge); 
  } 
  if ((-alphaDrSurf[11])+alphaDrSurf[7]-alphaDrSurf[6]-alphaDrSurf[5]+alphaDrSurf[3]+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[12] = gkhyb_3x2v_p1_surfx4_eval_quad_node_12_r(fSkin); 
    fUpwindQuad[13] = gkhyb_3x2v_p1_surfx4_eval_quad_node_13_r(fSkin); 
  } else { 
    fUpwindQuad[12] = gkhyb_3x2v_p1_surfx4_eval_quad_node_12_l(fEdge); 
    fUpwindQuad[13] = gkhyb_3x2v_p1_surfx4_eval_quad_node_13_l(fEdge); 
  } 
  if ((-alphaDrSurf[11])+alphaDrSurf[7]-alphaDrSurf[6]-alphaDrSurf[5]+alphaDrSurf[3]+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[14] = gkhyb_3x2v_p1_surfx4_eval_quad_node_14_r(fSkin); 
    fUpwindQuad[15] = gkhyb_3x2v_p1_surfx4_eval_quad_node_15_r(fSkin); 
  } else { 
    fUpwindQuad[14] = gkhyb_3x2v_p1_surfx4_eval_quad_node_14_l(fEdge); 
    fUpwindQuad[15] = gkhyb_3x2v_p1_surfx4_eval_quad_node_15_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  gkhyb_3x2v_p1_vpardir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.25*(alphaDrSurf[11]*fUpwind[11]+alphaDrSurf[7]*fUpwind[7]+alphaDrSurf[6]*fUpwind[6]+alphaDrSurf[5]*fUpwind[5]+alphaDrSurf[3]*fUpwind[3]+alphaDrSurf[2]*fUpwind[2]+alphaDrSurf[1]*fUpwind[1]+alphaDrSurf[0]*fUpwind[0]); 
  Ghat[1] = 0.25*(alphaDrSurf[7]*fUpwind[11]+fUpwind[7]*alphaDrSurf[11]+alphaDrSurf[3]*fUpwind[6]+fUpwind[3]*alphaDrSurf[6]+alphaDrSurf[2]*fUpwind[5]+fUpwind[2]*alphaDrSurf[5]+alphaDrSurf[0]*fUpwind[1]+fUpwind[0]*alphaDrSurf[1]); 
  Ghat[2] = 0.25*(alphaDrSurf[6]*fUpwind[11]+fUpwind[6]*alphaDrSurf[11]+alphaDrSurf[3]*fUpwind[7]+fUpwind[3]*alphaDrSurf[7]+alphaDrSurf[1]*fUpwind[5]+fUpwind[1]*alphaDrSurf[5]+alphaDrSurf[0]*fUpwind[2]+fUpwind[0]*alphaDrSurf[2]); 
  Ghat[3] = 0.25*(alphaDrSurf[5]*fUpwind[11]+fUpwind[5]*alphaDrSurf[11]+alphaDrSurf[2]*fUpwind[7]+fUpwind[2]*alphaDrSurf[7]+alphaDrSurf[1]*fUpwind[6]+fUpwind[1]*alphaDrSurf[6]+alphaDrSurf[0]*fUpwind[3]+fUpwind[0]*alphaDrSurf[3]); 
  Ghat[4] = 0.25*(alphaDrSurf[11]*fUpwind[15]+alphaDrSurf[7]*fUpwind[14]+alphaDrSurf[6]*fUpwind[13]+alphaDrSurf[5]*fUpwind[12]+alphaDrSurf[3]*fUpwind[10]+alphaDrSurf[2]*fUpwind[9]+alphaDrSurf[1]*fUpwind[8]+alphaDrSurf[0]*fUpwind[4]); 
  Ghat[5] = 0.25*(alphaDrSurf[3]*fUpwind[11]+fUpwind[3]*alphaDrSurf[11]+alphaDrSurf[6]*fUpwind[7]+fUpwind[6]*alphaDrSurf[7]+alphaDrSurf[0]*fUpwind[5]+fUpwind[0]*alphaDrSurf[5]+alphaDrSurf[1]*fUpwind[2]+fUpwind[1]*alphaDrSurf[2]); 
  Ghat[6] = 0.25*(alphaDrSurf[2]*fUpwind[11]+fUpwind[2]*alphaDrSurf[11]+alphaDrSurf[5]*fUpwind[7]+fUpwind[5]*alphaDrSurf[7]+alphaDrSurf[0]*fUpwind[6]+fUpwind[0]*alphaDrSurf[6]+alphaDrSurf[1]*fUpwind[3]+fUpwind[1]*alphaDrSurf[3]); 
  Ghat[7] = 0.25*(alphaDrSurf[1]*fUpwind[11]+fUpwind[1]*alphaDrSurf[11]+alphaDrSurf[0]*fUpwind[7]+fUpwind[0]*alphaDrSurf[7]+alphaDrSurf[5]*fUpwind[6]+fUpwind[5]*alphaDrSurf[6]+alphaDrSurf[2]*fUpwind[3]+fUpwind[2]*alphaDrSurf[3]); 
  Ghat[8] = 0.25*(alphaDrSurf[7]*fUpwind[15]+alphaDrSurf[11]*fUpwind[14]+alphaDrSurf[3]*fUpwind[13]+alphaDrSurf[2]*fUpwind[12]+alphaDrSurf[6]*fUpwind[10]+alphaDrSurf[5]*fUpwind[9]+alphaDrSurf[0]*fUpwind[8]+alphaDrSurf[1]*fUpwind[4]); 
  Ghat[9] = 0.25*(alphaDrSurf[6]*fUpwind[15]+alphaDrSurf[3]*fUpwind[14]+alphaDrSurf[11]*fUpwind[13]+alphaDrSurf[1]*fUpwind[12]+alphaDrSurf[7]*fUpwind[10]+alphaDrSurf[0]*fUpwind[9]+alphaDrSurf[5]*fUpwind[8]+alphaDrSurf[2]*fUpwind[4]); 
  Ghat[10] = 0.25*(alphaDrSurf[5]*fUpwind[15]+alphaDrSurf[2]*fUpwind[14]+alphaDrSurf[1]*fUpwind[13]+alphaDrSurf[11]*fUpwind[12]+alphaDrSurf[0]*fUpwind[10]+alphaDrSurf[7]*fUpwind[9]+alphaDrSurf[6]*fUpwind[8]+alphaDrSurf[3]*fUpwind[4]); 
  Ghat[11] = 0.25*(alphaDrSurf[0]*fUpwind[11]+fUpwind[0]*alphaDrSurf[11]+alphaDrSurf[1]*fUpwind[7]+fUpwind[1]*alphaDrSurf[7]+alphaDrSurf[2]*fUpwind[6]+fUpwind[2]*alphaDrSurf[6]+alphaDrSurf[3]*fUpwind[5]+fUpwind[3]*alphaDrSurf[5]); 
  Ghat[12] = 0.25*(alphaDrSurf[3]*fUpwind[15]+alphaDrSurf[6]*fUpwind[14]+alphaDrSurf[7]*fUpwind[13]+alphaDrSurf[0]*fUpwind[12]+fUpwind[10]*alphaDrSurf[11]+alphaDrSurf[1]*fUpwind[9]+alphaDrSurf[2]*fUpwind[8]+fUpwind[4]*alphaDrSurf[5]); 
  Ghat[13] = 0.25*(alphaDrSurf[2]*fUpwind[15]+alphaDrSurf[5]*fUpwind[14]+alphaDrSurf[0]*fUpwind[13]+alphaDrSurf[7]*fUpwind[12]+fUpwind[9]*alphaDrSurf[11]+alphaDrSurf[1]*fUpwind[10]+alphaDrSurf[3]*fUpwind[8]+fUpwind[4]*alphaDrSurf[6]); 
  Ghat[14] = 0.25*(alphaDrSurf[1]*fUpwind[15]+alphaDrSurf[0]*fUpwind[14]+alphaDrSurf[5]*fUpwind[13]+alphaDrSurf[6]*fUpwind[12]+fUpwind[8]*alphaDrSurf[11]+alphaDrSurf[2]*fUpwind[10]+alphaDrSurf[3]*fUpwind[9]+fUpwind[4]*alphaDrSurf[7]); 
  Ghat[15] = 0.25*(alphaDrSurf[0]*fUpwind[15]+alphaDrSurf[1]*fUpwind[14]+alphaDrSurf[2]*fUpwind[13]+alphaDrSurf[3]*fUpwind[12]+fUpwind[4]*alphaDrSurf[11]+alphaDrSurf[5]*fUpwind[10]+alphaDrSurf[6]*fUpwind[9]+alphaDrSurf[7]*fUpwind[8]); 

  out[0] += 0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += 0.7071067811865475*Ghat[3]*rdv2; 
  out[4] += 1.224744871391589*Ghat[0]*rdv2; 
  out[5] += 0.7071067811865475*Ghat[4]*rdv2; 
  out[6] += 0.7071067811865475*Ghat[5]*rdv2; 
  out[7] += 0.7071067811865475*Ghat[6]*rdv2; 
  out[8] += 0.7071067811865475*Ghat[7]*rdv2; 
  out[9] += 1.224744871391589*Ghat[1]*rdv2; 
  out[10] += 1.224744871391589*Ghat[2]*rdv2; 
  out[11] += 1.224744871391589*Ghat[3]*rdv2; 
  out[12] += 0.7071067811865475*Ghat[8]*rdv2; 
  out[13] += 0.7071067811865475*Ghat[9]*rdv2; 
  out[14] += 0.7071067811865475*Ghat[10]*rdv2; 
  out[15] += 1.224744871391589*Ghat[4]*rdv2; 
  out[16] += 0.7071067811865475*Ghat[11]*rdv2; 
  out[17] += 1.224744871391589*Ghat[5]*rdv2; 
  out[18] += 1.224744871391589*Ghat[6]*rdv2; 
  out[19] += 1.224744871391589*Ghat[7]*rdv2; 
  out[20] += 0.7071067811865475*Ghat[12]*rdv2; 
  out[21] += 0.7071067811865475*Ghat[13]*rdv2; 
  out[22] += 0.7071067811865475*Ghat[14]*rdv2; 
  out[23] += 1.224744871391589*Ghat[8]*rdv2; 
  out[24] += 1.224744871391589*Ghat[9]*rdv2; 
  out[25] += 1.224744871391589*Ghat[10]*rdv2; 
  out[26] += 1.224744871391589*Ghat[11]*rdv2; 
  out[27] += 0.7071067811865475*Ghat[15]*rdv2; 
  out[28] += 1.224744871391589*Ghat[12]*rdv2; 
  out[29] += 1.224744871391589*Ghat[13]*rdv2; 
  out[30] += 1.224744871391589*Ghat[14]*rdv2; 
  out[31] += 1.224744871391589*Ghat[15]*rdv2; 
  out[32] += 1.58113883008419*Ghat[0]*rdv2; 
  out[33] += 1.58113883008419*Ghat[1]*rdv2; 
  out[34] += 1.58113883008419*Ghat[2]*rdv2; 
  out[35] += 1.58113883008419*Ghat[3]*rdv2; 
  out[36] += 1.58113883008419*Ghat[4]*rdv2; 
  out[37] += 1.58113883008419*Ghat[5]*rdv2; 
  out[38] += 1.58113883008419*Ghat[6]*rdv2; 
  out[39] += 1.58113883008419*Ghat[7]*rdv2; 
  out[40] += 1.58113883008419*Ghat[8]*rdv2; 
  out[41] += 1.58113883008419*Ghat[9]*rdv2; 
  out[42] += 1.58113883008419*Ghat[10]*rdv2; 
  out[43] += 1.58113883008419*Ghat[11]*rdv2; 
  out[44] += 1.58113883008419*Ghat[12]*rdv2; 
  out[45] += 1.58113883008419*Ghat[13]*rdv2; 
  out[46] += 1.58113883008419*Ghat[14]*rdv2; 
  out[47] += 1.58113883008419*Ghat[15]*rdv2; 

  } else { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*(2.828427124746191*w[3]-1.414213562373095*dxv[3])-2.828427124746191*nuUSum[0]); 
  alphaDrSurf[1] = 0.5*(nuSum[1]*(2.828427124746191*w[3]-1.414213562373095*dxv[3])-2.828427124746191*nuUSum[1]); 
  alphaDrSurf[2] = 0.5*(nuSum[2]*(2.828427124746191*w[3]-1.414213562373095*dxv[3])-2.828427124746191*nuUSum[2]); 
  alphaDrSurf[3] = 0.5*(2.828427124746191*nuSum[3]*w[3]-2.828427124746191*nuUSum[3]-1.414213562373095*dxv[3]*nuSum[3]); 
  alphaDrSurf[5] = -0.5*(2.828427124746191*nuUSum[4]+(1.414213562373095*dxv[3]-2.828427124746191*w[3])*nuSum[4]); 
  alphaDrSurf[6] = -0.5*(2.828427124746191*nuUSum[5]+(1.414213562373095*dxv[3]-2.828427124746191*w[3])*nuSum[5]); 
  alphaDrSurf[7] = -0.5*(2.828427124746191*nuUSum[6]+(1.414213562373095*dxv[3]-2.828427124746191*w[3])*nuSum[6]); 
  alphaDrSurf[11] = -0.5*(2.828427124746191*nuUSum[7]+(1.414213562373095*dxv[3]-2.828427124746191*w[3])*nuSum[7]); 

  if ((-alphaDrSurf[11])+alphaDrSurf[7]+alphaDrSurf[6]+alphaDrSurf[5]-alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = gkhyb_3x2v_p1_surfx4_eval_quad_node_0_r(fEdge); 
    fUpwindQuad[1] = gkhyb_3x2v_p1_surfx4_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[0] = gkhyb_3x2v_p1_surfx4_eval_quad_node_0_l(fSkin); 
    fUpwindQuad[1] = gkhyb_3x2v_p1_surfx4_eval_quad_node_1_l(fSkin); 
  } 
  if ((-alphaDrSurf[11])+alphaDrSurf[7]+alphaDrSurf[6]+alphaDrSurf[5]-alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[2] = gkhyb_3x2v_p1_surfx4_eval_quad_node_2_r(fEdge); 
    fUpwindQuad[3] = gkhyb_3x2v_p1_surfx4_eval_quad_node_3_r(fEdge); 
  } else { 
    fUpwindQuad[2] = gkhyb_3x2v_p1_surfx4_eval_quad_node_2_l(fSkin); 
    fUpwindQuad[3] = gkhyb_3x2v_p1_surfx4_eval_quad_node_3_l(fSkin); 
  } 
  if (alphaDrSurf[11]-alphaDrSurf[7]-alphaDrSurf[6]+alphaDrSurf[5]+alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[4] = gkhyb_3x2v_p1_surfx4_eval_quad_node_4_r(fEdge); 
    fUpwindQuad[5] = gkhyb_3x2v_p1_surfx4_eval_quad_node_5_r(fEdge); 
  } else { 
    fUpwindQuad[4] = gkhyb_3x2v_p1_surfx4_eval_quad_node_4_l(fSkin); 
    fUpwindQuad[5] = gkhyb_3x2v_p1_surfx4_eval_quad_node_5_l(fSkin); 
  } 
  if (alphaDrSurf[11]-alphaDrSurf[7]-alphaDrSurf[6]+alphaDrSurf[5]+alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[6] = gkhyb_3x2v_p1_surfx4_eval_quad_node_6_r(fEdge); 
    fUpwindQuad[7] = gkhyb_3x2v_p1_surfx4_eval_quad_node_7_r(fEdge); 
  } else { 
    fUpwindQuad[6] = gkhyb_3x2v_p1_surfx4_eval_quad_node_6_l(fSkin); 
    fUpwindQuad[7] = gkhyb_3x2v_p1_surfx4_eval_quad_node_7_l(fSkin); 
  } 
  if (alphaDrSurf[11]-alphaDrSurf[7]+alphaDrSurf[6]-alphaDrSurf[5]-alphaDrSurf[3]+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[8] = gkhyb_3x2v_p1_surfx4_eval_quad_node_8_r(fEdge); 
    fUpwindQuad[9] = gkhyb_3x2v_p1_surfx4_eval_quad_node_9_r(fEdge); 
  } else { 
    fUpwindQuad[8] = gkhyb_3x2v_p1_surfx4_eval_quad_node_8_l(fSkin); 
    fUpwindQuad[9] = gkhyb_3x2v_p1_surfx4_eval_quad_node_9_l(fSkin); 
  } 
  if (alphaDrSurf[11]-alphaDrSurf[7]+alphaDrSurf[6]-alphaDrSurf[5]-alphaDrSurf[3]+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[10] = gkhyb_3x2v_p1_surfx4_eval_quad_node_10_r(fEdge); 
    fUpwindQuad[11] = gkhyb_3x2v_p1_surfx4_eval_quad_node_11_r(fEdge); 
  } else { 
    fUpwindQuad[10] = gkhyb_3x2v_p1_surfx4_eval_quad_node_10_l(fSkin); 
    fUpwindQuad[11] = gkhyb_3x2v_p1_surfx4_eval_quad_node_11_l(fSkin); 
  } 
  if ((-alphaDrSurf[11])+alphaDrSurf[7]-alphaDrSurf[6]-alphaDrSurf[5]+alphaDrSurf[3]+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[12] = gkhyb_3x2v_p1_surfx4_eval_quad_node_12_r(fEdge); 
    fUpwindQuad[13] = gkhyb_3x2v_p1_surfx4_eval_quad_node_13_r(fEdge); 
  } else { 
    fUpwindQuad[12] = gkhyb_3x2v_p1_surfx4_eval_quad_node_12_l(fSkin); 
    fUpwindQuad[13] = gkhyb_3x2v_p1_surfx4_eval_quad_node_13_l(fSkin); 
  } 
  if ((-alphaDrSurf[11])+alphaDrSurf[7]-alphaDrSurf[6]-alphaDrSurf[5]+alphaDrSurf[3]+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[14] = gkhyb_3x2v_p1_surfx4_eval_quad_node_14_r(fEdge); 
    fUpwindQuad[15] = gkhyb_3x2v_p1_surfx4_eval_quad_node_15_r(fEdge); 
  } else { 
    fUpwindQuad[14] = gkhyb_3x2v_p1_surfx4_eval_quad_node_14_l(fSkin); 
    fUpwindQuad[15] = gkhyb_3x2v_p1_surfx4_eval_quad_node_15_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  gkhyb_3x2v_p1_vpardir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.25*(alphaDrSurf[11]*fUpwind[11]+alphaDrSurf[7]*fUpwind[7]+alphaDrSurf[6]*fUpwind[6]+alphaDrSurf[5]*fUpwind[5]+alphaDrSurf[3]*fUpwind[3]+alphaDrSurf[2]*fUpwind[2]+alphaDrSurf[1]*fUpwind[1]+alphaDrSurf[0]*fUpwind[0]); 
  Ghat[1] = 0.25*(alphaDrSurf[7]*fUpwind[11]+fUpwind[7]*alphaDrSurf[11]+alphaDrSurf[3]*fUpwind[6]+fUpwind[3]*alphaDrSurf[6]+alphaDrSurf[2]*fUpwind[5]+fUpwind[2]*alphaDrSurf[5]+alphaDrSurf[0]*fUpwind[1]+fUpwind[0]*alphaDrSurf[1]); 
  Ghat[2] = 0.25*(alphaDrSurf[6]*fUpwind[11]+fUpwind[6]*alphaDrSurf[11]+alphaDrSurf[3]*fUpwind[7]+fUpwind[3]*alphaDrSurf[7]+alphaDrSurf[1]*fUpwind[5]+fUpwind[1]*alphaDrSurf[5]+alphaDrSurf[0]*fUpwind[2]+fUpwind[0]*alphaDrSurf[2]); 
  Ghat[3] = 0.25*(alphaDrSurf[5]*fUpwind[11]+fUpwind[5]*alphaDrSurf[11]+alphaDrSurf[2]*fUpwind[7]+fUpwind[2]*alphaDrSurf[7]+alphaDrSurf[1]*fUpwind[6]+fUpwind[1]*alphaDrSurf[6]+alphaDrSurf[0]*fUpwind[3]+fUpwind[0]*alphaDrSurf[3]); 
  Ghat[4] = 0.25*(alphaDrSurf[11]*fUpwind[15]+alphaDrSurf[7]*fUpwind[14]+alphaDrSurf[6]*fUpwind[13]+alphaDrSurf[5]*fUpwind[12]+alphaDrSurf[3]*fUpwind[10]+alphaDrSurf[2]*fUpwind[9]+alphaDrSurf[1]*fUpwind[8]+alphaDrSurf[0]*fUpwind[4]); 
  Ghat[5] = 0.25*(alphaDrSurf[3]*fUpwind[11]+fUpwind[3]*alphaDrSurf[11]+alphaDrSurf[6]*fUpwind[7]+fUpwind[6]*alphaDrSurf[7]+alphaDrSurf[0]*fUpwind[5]+fUpwind[0]*alphaDrSurf[5]+alphaDrSurf[1]*fUpwind[2]+fUpwind[1]*alphaDrSurf[2]); 
  Ghat[6] = 0.25*(alphaDrSurf[2]*fUpwind[11]+fUpwind[2]*alphaDrSurf[11]+alphaDrSurf[5]*fUpwind[7]+fUpwind[5]*alphaDrSurf[7]+alphaDrSurf[0]*fUpwind[6]+fUpwind[0]*alphaDrSurf[6]+alphaDrSurf[1]*fUpwind[3]+fUpwind[1]*alphaDrSurf[3]); 
  Ghat[7] = 0.25*(alphaDrSurf[1]*fUpwind[11]+fUpwind[1]*alphaDrSurf[11]+alphaDrSurf[0]*fUpwind[7]+fUpwind[0]*alphaDrSurf[7]+alphaDrSurf[5]*fUpwind[6]+fUpwind[5]*alphaDrSurf[6]+alphaDrSurf[2]*fUpwind[3]+fUpwind[2]*alphaDrSurf[3]); 
  Ghat[8] = 0.25*(alphaDrSurf[7]*fUpwind[15]+alphaDrSurf[11]*fUpwind[14]+alphaDrSurf[3]*fUpwind[13]+alphaDrSurf[2]*fUpwind[12]+alphaDrSurf[6]*fUpwind[10]+alphaDrSurf[5]*fUpwind[9]+alphaDrSurf[0]*fUpwind[8]+alphaDrSurf[1]*fUpwind[4]); 
  Ghat[9] = 0.25*(alphaDrSurf[6]*fUpwind[15]+alphaDrSurf[3]*fUpwind[14]+alphaDrSurf[11]*fUpwind[13]+alphaDrSurf[1]*fUpwind[12]+alphaDrSurf[7]*fUpwind[10]+alphaDrSurf[0]*fUpwind[9]+alphaDrSurf[5]*fUpwind[8]+alphaDrSurf[2]*fUpwind[4]); 
  Ghat[10] = 0.25*(alphaDrSurf[5]*fUpwind[15]+alphaDrSurf[2]*fUpwind[14]+alphaDrSurf[1]*fUpwind[13]+alphaDrSurf[11]*fUpwind[12]+alphaDrSurf[0]*fUpwind[10]+alphaDrSurf[7]*fUpwind[9]+alphaDrSurf[6]*fUpwind[8]+alphaDrSurf[3]*fUpwind[4]); 
  Ghat[11] = 0.25*(alphaDrSurf[0]*fUpwind[11]+fUpwind[0]*alphaDrSurf[11]+alphaDrSurf[1]*fUpwind[7]+fUpwind[1]*alphaDrSurf[7]+alphaDrSurf[2]*fUpwind[6]+fUpwind[2]*alphaDrSurf[6]+alphaDrSurf[3]*fUpwind[5]+fUpwind[3]*alphaDrSurf[5]); 
  Ghat[12] = 0.25*(alphaDrSurf[3]*fUpwind[15]+alphaDrSurf[6]*fUpwind[14]+alphaDrSurf[7]*fUpwind[13]+alphaDrSurf[0]*fUpwind[12]+fUpwind[10]*alphaDrSurf[11]+alphaDrSurf[1]*fUpwind[9]+alphaDrSurf[2]*fUpwind[8]+fUpwind[4]*alphaDrSurf[5]); 
  Ghat[13] = 0.25*(alphaDrSurf[2]*fUpwind[15]+alphaDrSurf[5]*fUpwind[14]+alphaDrSurf[0]*fUpwind[13]+alphaDrSurf[7]*fUpwind[12]+fUpwind[9]*alphaDrSurf[11]+alphaDrSurf[1]*fUpwind[10]+alphaDrSurf[3]*fUpwind[8]+fUpwind[4]*alphaDrSurf[6]); 
  Ghat[14] = 0.25*(alphaDrSurf[1]*fUpwind[15]+alphaDrSurf[0]*fUpwind[14]+alphaDrSurf[5]*fUpwind[13]+alphaDrSurf[6]*fUpwind[12]+fUpwind[8]*alphaDrSurf[11]+alphaDrSurf[2]*fUpwind[10]+alphaDrSurf[3]*fUpwind[9]+fUpwind[4]*alphaDrSurf[7]); 
  Ghat[15] = 0.25*(alphaDrSurf[0]*fUpwind[15]+alphaDrSurf[1]*fUpwind[14]+alphaDrSurf[2]*fUpwind[13]+alphaDrSurf[3]*fUpwind[12]+fUpwind[4]*alphaDrSurf[11]+alphaDrSurf[5]*fUpwind[10]+alphaDrSurf[6]*fUpwind[9]+alphaDrSurf[7]*fUpwind[8]); 

  out[0] += -0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += -0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += -0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += -0.7071067811865475*Ghat[3]*rdv2; 
  out[4] += 1.224744871391589*Ghat[0]*rdv2; 
  out[5] += -0.7071067811865475*Ghat[4]*rdv2; 
  out[6] += -0.7071067811865475*Ghat[5]*rdv2; 
  out[7] += -0.7071067811865475*Ghat[6]*rdv2; 
  out[8] += -0.7071067811865475*Ghat[7]*rdv2; 
  out[9] += 1.224744871391589*Ghat[1]*rdv2; 
  out[10] += 1.224744871391589*Ghat[2]*rdv2; 
  out[11] += 1.224744871391589*Ghat[3]*rdv2; 
  out[12] += -0.7071067811865475*Ghat[8]*rdv2; 
  out[13] += -0.7071067811865475*Ghat[9]*rdv2; 
  out[14] += -0.7071067811865475*Ghat[10]*rdv2; 
  out[15] += 1.224744871391589*Ghat[4]*rdv2; 
  out[16] += -0.7071067811865475*Ghat[11]*rdv2; 
  out[17] += 1.224744871391589*Ghat[5]*rdv2; 
  out[18] += 1.224744871391589*Ghat[6]*rdv2; 
  out[19] += 1.224744871391589*Ghat[7]*rdv2; 
  out[20] += -0.7071067811865475*Ghat[12]*rdv2; 
  out[21] += -0.7071067811865475*Ghat[13]*rdv2; 
  out[22] += -0.7071067811865475*Ghat[14]*rdv2; 
  out[23] += 1.224744871391589*Ghat[8]*rdv2; 
  out[24] += 1.224744871391589*Ghat[9]*rdv2; 
  out[25] += 1.224744871391589*Ghat[10]*rdv2; 
  out[26] += 1.224744871391589*Ghat[11]*rdv2; 
  out[27] += -0.7071067811865475*Ghat[15]*rdv2; 
  out[28] += 1.224744871391589*Ghat[12]*rdv2; 
  out[29] += 1.224744871391589*Ghat[13]*rdv2; 
  out[30] += 1.224744871391589*Ghat[14]*rdv2; 
  out[31] += 1.224744871391589*Ghat[15]*rdv2; 
  out[32] += -1.58113883008419*Ghat[0]*rdv2; 
  out[33] += -1.58113883008419*Ghat[1]*rdv2; 
  out[34] += -1.58113883008419*Ghat[2]*rdv2; 
  out[35] += -1.58113883008419*Ghat[3]*rdv2; 
  out[36] += -1.58113883008419*Ghat[4]*rdv2; 
  out[37] += -1.58113883008419*Ghat[5]*rdv2; 
  out[38] += -1.58113883008419*Ghat[6]*rdv2; 
  out[39] += -1.58113883008419*Ghat[7]*rdv2; 
  out[40] += -1.58113883008419*Ghat[8]*rdv2; 
  out[41] += -1.58113883008419*Ghat[9]*rdv2; 
  out[42] += -1.58113883008419*Ghat[10]*rdv2; 
  out[43] += -1.58113883008419*Ghat[11]*rdv2; 
  out[44] += -1.58113883008419*Ghat[12]*rdv2; 
  out[45] += -1.58113883008419*Ghat[13]*rdv2; 
  out[46] += -1.58113883008419*Ghat[14]*rdv2; 
  out[47] += -1.58113883008419*Ghat[15]*rdv2; 

  } 
} 
