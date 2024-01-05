#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_hyb_1x2v_p1_surfx3_eval_quad.h> 
#include <gkyl_basis_hyb_1x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_surfvy_1x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // field:     q/m*EM fields.
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 
  const double dv11 = 2/dxv[2]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double *E1 = &field[2]; 
  const double *B2 = &field[10]; 

  double alpha[6] = {0.0}; 

  alpha[0] = 1.414213562373095*E1[0]-1.414213562373095*B2[0]*wv1; 
  alpha[1] = 1.414213562373095*E1[1]-1.414213562373095*B2[1]*wv1; 
  alpha[2] = -0.408248290463863*B2[0]*dv1; 
  alpha[3] = -0.408248290463863*B2[1]*dv1; 

  double fUpwindQuad_l[6] = {0.0};
  double fUpwindQuad_r[6] = {0.0};
  double fUpwind_l[6] = {0.0};;
  double fUpwind_r[6] = {0.0};
  double Ghat_l[6] = {0.0}; 
  double Ghat_r[6] = {0.0}; 

  if (0.6708203932499369*alpha[3]-0.6708203932499369*alpha[2]-0.5*alpha[1]+0.5*alpha[0] > 0) { 
    fUpwindQuad_l[0] = hyb_1x2v_p1_surfx3_eval_quad_node_0_r(fl); 
    fUpwindQuad_r[0] = hyb_1x2v_p1_surfx3_eval_quad_node_0_r(fc); 
  } else { 
    fUpwindQuad_l[0] = hyb_1x2v_p1_surfx3_eval_quad_node_0_l(fc); 
    fUpwindQuad_r[0] = hyb_1x2v_p1_surfx3_eval_quad_node_0_l(fr); 
  } 
  if (0.5*alpha[0]-0.5*alpha[1] > 0) { 
    fUpwindQuad_l[1] = hyb_1x2v_p1_surfx3_eval_quad_node_1_r(fl); 
    fUpwindQuad_r[1] = hyb_1x2v_p1_surfx3_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_l[1] = hyb_1x2v_p1_surfx3_eval_quad_node_1_l(fc); 
    fUpwindQuad_r[1] = hyb_1x2v_p1_surfx3_eval_quad_node_1_l(fr); 
  } 
  if ((-0.6708203932499369*alpha[3])+0.6708203932499369*alpha[2]-0.5*alpha[1]+0.5*alpha[0] > 0) { 
    fUpwindQuad_l[2] = hyb_1x2v_p1_surfx3_eval_quad_node_2_r(fl); 
    fUpwindQuad_r[2] = hyb_1x2v_p1_surfx3_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_l[2] = hyb_1x2v_p1_surfx3_eval_quad_node_2_l(fc); 
    fUpwindQuad_r[2] = hyb_1x2v_p1_surfx3_eval_quad_node_2_l(fr); 
  } 
  if (0.5*(alpha[1]+alpha[0])-0.6708203932499369*(alpha[3]+alpha[2]) > 0) { 
    fUpwindQuad_l[3] = hyb_1x2v_p1_surfx3_eval_quad_node_3_r(fl); 
    fUpwindQuad_r[3] = hyb_1x2v_p1_surfx3_eval_quad_node_3_r(fc); 
  } else { 
    fUpwindQuad_l[3] = hyb_1x2v_p1_surfx3_eval_quad_node_3_l(fc); 
    fUpwindQuad_r[3] = hyb_1x2v_p1_surfx3_eval_quad_node_3_l(fr); 
  } 
  if (0.5*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[4] = hyb_1x2v_p1_surfx3_eval_quad_node_4_r(fl); 
    fUpwindQuad_r[4] = hyb_1x2v_p1_surfx3_eval_quad_node_4_r(fc); 
  } else { 
    fUpwindQuad_l[4] = hyb_1x2v_p1_surfx3_eval_quad_node_4_l(fc); 
    fUpwindQuad_r[4] = hyb_1x2v_p1_surfx3_eval_quad_node_4_l(fr); 
  } 
  if (0.6708203932499369*(alpha[3]+alpha[2])+0.5*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[5] = hyb_1x2v_p1_surfx3_eval_quad_node_5_r(fl); 
    fUpwindQuad_r[5] = hyb_1x2v_p1_surfx3_eval_quad_node_5_r(fc); 
  } else { 
    fUpwindQuad_l[5] = hyb_1x2v_p1_surfx3_eval_quad_node_5_l(fc); 
    fUpwindQuad_r[5] = hyb_1x2v_p1_surfx3_eval_quad_node_5_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_1x2v_p1_vdir_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  hyb_1x2v_p1_vdir_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.5*alpha[3]*fUpwind_l[3]+0.5*alpha[2]*fUpwind_l[2]+0.5*alpha[1]*fUpwind_l[1]+0.5*alpha[0]*fUpwind_l[0]; 
  Ghat_l[1] = 0.5*alpha[2]*fUpwind_l[3]+0.5*fUpwind_l[2]*alpha[3]+0.5*alpha[0]*fUpwind_l[1]+0.5*fUpwind_l[0]*alpha[1]; 
  Ghat_l[2] = 0.447213595499958*alpha[3]*fUpwind_l[5]+0.4472135954999579*alpha[2]*fUpwind_l[4]+0.5*alpha[1]*fUpwind_l[3]+0.5*fUpwind_l[1]*alpha[3]+0.5*alpha[0]*fUpwind_l[2]+0.5*fUpwind_l[0]*alpha[2]; 
  Ghat_l[3] = 0.447213595499958*alpha[2]*fUpwind_l[5]+0.4472135954999579*alpha[3]*fUpwind_l[4]+0.5*alpha[0]*fUpwind_l[3]+0.5*fUpwind_l[0]*alpha[3]+0.5*alpha[1]*fUpwind_l[2]+0.5*fUpwind_l[1]*alpha[2]; 
  Ghat_l[4] = 0.5000000000000001*alpha[1]*fUpwind_l[5]+0.5*alpha[0]*fUpwind_l[4]+0.4472135954999579*alpha[3]*fUpwind_l[3]+0.4472135954999579*alpha[2]*fUpwind_l[2]; 
  Ghat_l[5] = 0.5*alpha[0]*fUpwind_l[5]+0.5000000000000001*alpha[1]*fUpwind_l[4]+0.447213595499958*alpha[2]*fUpwind_l[3]+0.447213595499958*fUpwind_l[2]*alpha[3]; 

  Ghat_r[0] = 0.5*alpha[3]*fUpwind_r[3]+0.5*alpha[2]*fUpwind_r[2]+0.5*alpha[1]*fUpwind_r[1]+0.5*alpha[0]*fUpwind_r[0]; 
  Ghat_r[1] = 0.5*alpha[2]*fUpwind_r[3]+0.5*fUpwind_r[2]*alpha[3]+0.5*alpha[0]*fUpwind_r[1]+0.5*fUpwind_r[0]*alpha[1]; 
  Ghat_r[2] = 0.447213595499958*alpha[3]*fUpwind_r[5]+0.4472135954999579*alpha[2]*fUpwind_r[4]+0.5*alpha[1]*fUpwind_r[3]+0.5*fUpwind_r[1]*alpha[3]+0.5*alpha[0]*fUpwind_r[2]+0.5*fUpwind_r[0]*alpha[2]; 
  Ghat_r[3] = 0.447213595499958*alpha[2]*fUpwind_r[5]+0.4472135954999579*alpha[3]*fUpwind_r[4]+0.5*alpha[0]*fUpwind_r[3]+0.5*fUpwind_r[0]*alpha[3]+0.5*alpha[1]*fUpwind_r[2]+0.5*fUpwind_r[1]*alpha[2]; 
  Ghat_r[4] = 0.5000000000000001*alpha[1]*fUpwind_r[5]+0.5*alpha[0]*fUpwind_r[4]+0.4472135954999579*alpha[3]*fUpwind_r[3]+0.4472135954999579*alpha[2]*fUpwind_r[2]; 
  Ghat_r[5] = 0.5*alpha[0]*fUpwind_r[5]+0.5000000000000001*alpha[1]*fUpwind_r[4]+0.447213595499958*alpha[2]*fUpwind_r[3]+0.447213595499958*fUpwind_r[2]*alpha[3]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv11; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv11; 
  out[2] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv11; 
  out[3] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv11; 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv11; 
  out[5] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv11; 
  out[6] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv11; 
  out[7] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv11; 
  out[8] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dv11; 
  out[9] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dv11; 
  out[10] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dv11; 
  out[11] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dv11; 
  out[12] += (1.58113883008419*Ghat_l[0]-1.58113883008419*Ghat_r[0])*dv11; 
  out[13] += (1.58113883008419*Ghat_l[1]-1.58113883008419*Ghat_r[1])*dv11; 
  out[14] += (1.58113883008419*Ghat_l[2]-1.58113883008419*Ghat_r[2])*dv11; 
  out[15] += (1.58113883008419*Ghat_l[3]-1.58113883008419*Ghat_r[3])*dv11; 

  return 0.;

} 
