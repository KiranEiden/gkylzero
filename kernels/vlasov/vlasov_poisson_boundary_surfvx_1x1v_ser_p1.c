#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_hyb_1x1v_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_hyb_1x1v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_poisson_boundary_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *fac_phi, const double *vecA, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // fac_phi:     potential (scaled by appropriate factors).
  // vecA:        vector potential (scaled by appropriate factors). Unused in pure Vlasov-Poisson. 
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv10 = 2/dxv[1]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double *phi = &fac_phi[0]; 
  const double dx10 = 2/dxv[0]; 
  double alpha[2] = {0.0}; 

  alpha[0] = -1.732050807568877*phi[1]*dx10; 

  double fUpwindQuad[2] = {0.0};
  double fUpwind[2] = {0.0};
  double Ghat[2] = {0.0}; 

  if (edge == -1) { 

  if (alpha[0] > 0) { 
    fUpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(fEdge); 
  } 
  if (alpha[0] > 0) { 
    fUpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.7071067811865475*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.7071067811865475*alpha[0]*fUpwind[1]; 

  out[0] += -0.7071067811865475*Ghat[0]*dv10; 
  out[1] += -0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -1.224744871391589*Ghat[0]*dv10; 
  out[3] += -1.224744871391589*Ghat[1]*dv10; 
  out[4] += -1.58113883008419*Ghat[0]*dv10; 
  out[5] += -1.58113883008419*Ghat[1]*dv10; 

  } else { 

  if (alpha[0] > 0) { 
    fUpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(fSkin); 
  } 
  if (alpha[0] > 0) { 
    fUpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.7071067811865475*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.7071067811865475*alpha[0]*fUpwind[1]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv10; 
  out[1] += 0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -1.224744871391589*Ghat[0]*dv10; 
  out[3] += -1.224744871391589*Ghat[1]*dv10; 
  out[4] += 1.58113883008419*Ghat[0]*dv10; 
  out[5] += 1.58113883008419*Ghat[1]*dv10; 

  } 
} 
