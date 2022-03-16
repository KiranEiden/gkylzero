#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_ser_3x_p1_surfx3_eval_quad.h> 
#include <gkyl_basis_ser_3x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void lbo_vlasov_drag_surfvy_1x2v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[3]:         cell-center coordinates. 
  // dxv[3]:       cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[4]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      distribution function in cells 
  // out:           incremented distribution function in cell 
  double rdv2 = 2.0/dxv[2]; 

  const double *sumNuUy = &nuUSum[2]; 

  double alphaDrSurf_l[4] = {0.0}; 
  alphaDrSurf_l[0] = 1.414213562373095*nuSum[0]*w[2]-0.7071067811865475*nuSum[0]*dxv[2]-1.414213562373095*sumNuUy[0]; 
  alphaDrSurf_l[1] = 1.414213562373095*nuSum[1]*w[2]-0.7071067811865475*nuSum[1]*dxv[2]-1.414213562373095*sumNuUy[1]; 

  double alphaDrSurf_r[4] = {0.0}; 
  alphaDrSurf_r[0] = 1.414213562373095*nuSum[0]*w[2]+0.7071067811865475*nuSum[0]*dxv[2]-1.414213562373095*sumNuUy[0]; 
  alphaDrSurf_r[1] = 1.414213562373095*nuSum[1]*w[2]+0.7071067811865475*nuSum[1]*dxv[2]-1.414213562373095*sumNuUy[1]; 

  double fUpwindQuad_l[4] = {0.0};
  double fUpwindQuad_r[4] = {0.0};
  double fUpwind_l[4] = {0.0};
  double fUpwind_r[4] = {0.0};
  double Ghat_l[4] = {0.0}; 
  double Ghat_r[4] = {0.0}; 

  if (alphaDrSurf_l[0]-alphaDrSurf_l[1] < 0) { 
    fUpwindQuad_l[0] = ser_3x_p1_surfx3_eval_quad_node_0_r(fl); 
    fUpwindQuad_l[1] = ser_3x_p1_surfx3_eval_quad_node_1_r(fl); 
  } else { 
    fUpwindQuad_l[0] = ser_3x_p1_surfx3_eval_quad_node_0_l(fc); 
    fUpwindQuad_l[1] = ser_3x_p1_surfx3_eval_quad_node_1_l(fc); 
  } 
  if (alphaDrSurf_r[0]-alphaDrSurf_r[1] < 0) { 
    fUpwindQuad_r[0] = ser_3x_p1_surfx3_eval_quad_node_0_r(fc); 
    fUpwindQuad_r[1] = ser_3x_p1_surfx3_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_r[0] = ser_3x_p1_surfx3_eval_quad_node_0_l(fr); 
    fUpwindQuad_r[1] = ser_3x_p1_surfx3_eval_quad_node_1_l(fr); 
  } 
  if (alphaDrSurf_l[0]-alphaDrSurf_l[1] < 0) { 
    fUpwindQuad_l[2] = ser_3x_p1_surfx3_eval_quad_node_2_r(fl); 
    fUpwindQuad_l[3] = ser_3x_p1_surfx3_eval_quad_node_3_r(fl); 
  } else { 
    fUpwindQuad_l[2] = ser_3x_p1_surfx3_eval_quad_node_2_l(fc); 
    fUpwindQuad_l[3] = ser_3x_p1_surfx3_eval_quad_node_3_l(fc); 
  } 
  if (alphaDrSurf_r[0]-alphaDrSurf_r[1] < 0) { 
    fUpwindQuad_r[2] = ser_3x_p1_surfx3_eval_quad_node_2_r(fc); 
    fUpwindQuad_r[3] = ser_3x_p1_surfx3_eval_quad_node_3_r(fc); 
  } else { 
    fUpwindQuad_r[2] = ser_3x_p1_surfx3_eval_quad_node_2_l(fr); 
    fUpwindQuad_r[3] = ser_3x_p1_surfx3_eval_quad_node_3_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_3x_p1_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  ser_3x_p1_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.5*(alphaDrSurf_l[1]*fUpwind_l[1]+alphaDrSurf_l[0]*fUpwind_l[0]); 
  Ghat_l[1] = 0.5*(alphaDrSurf_l[0]*fUpwind_l[1]+fUpwind_l[0]*alphaDrSurf_l[1]); 
  Ghat_l[2] = 0.5*(alphaDrSurf_l[1]*fUpwind_l[3]+alphaDrSurf_l[0]*fUpwind_l[2]); 
  Ghat_l[3] = 0.5*(alphaDrSurf_l[0]*fUpwind_l[3]+alphaDrSurf_l[1]*fUpwind_l[2]); 

  Ghat_r[0] = 0.5*(alphaDrSurf_r[1]*fUpwind_r[1]+alphaDrSurf_r[0]*fUpwind_r[0]); 
  Ghat_r[1] = 0.5*(alphaDrSurf_r[0]*fUpwind_r[1]+fUpwind_r[0]*alphaDrSurf_r[1]); 
  Ghat_r[2] = 0.5*(alphaDrSurf_r[1]*fUpwind_r[3]+alphaDrSurf_r[0]*fUpwind_r[2]); 
  Ghat_r[3] = 0.5*(alphaDrSurf_r[0]*fUpwind_r[3]+alphaDrSurf_r[1]*fUpwind_r[2]); 

  out[0] += 0.7071067811865475*Ghat_r[0]*rdv2-0.7071067811865475*Ghat_l[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat_r[1]*rdv2-0.7071067811865475*Ghat_l[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat_r[2]*rdv2-0.7071067811865475*Ghat_l[2]*rdv2; 
  out[3] += 1.224744871391589*Ghat_r[0]*rdv2+1.224744871391589*Ghat_l[0]*rdv2; 
  out[4] += 0.7071067811865475*Ghat_r[3]*rdv2-0.7071067811865475*Ghat_l[3]*rdv2; 
  out[5] += 1.224744871391589*Ghat_r[1]*rdv2+1.224744871391589*Ghat_l[1]*rdv2; 
  out[6] += 1.224744871391589*Ghat_r[2]*rdv2+1.224744871391589*Ghat_l[2]*rdv2; 
  out[7] += 1.224744871391589*Ghat_r[3]*rdv2+1.224744871391589*Ghat_l[3]*rdv2; 
} 