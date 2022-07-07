#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_4x_p1_surfx3_eval_quad.h> 
#include <gkyl_basis_ser_4x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_surfvy_1x3v_ser_p1(const double *w, const double *dxv, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // qmem:      q/m*EM fields.
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 
  const double dv11 = 2/dxv[2]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double dv3 = dxv[3], wv3 = w[3]; 
  const double *E1 = &qmem[2]; 
  const double *B0 = &qmem[6]; 
  const double *B1 = &qmem[8]; 
  const double *B2 = &qmem[10]; 

  double alpha[16] = {0.0}; 

  alpha[0] = 2.0*B0[0]*wv3-2.0*B2[0]*wv1+2.0*E1[0]; 
  alpha[1] = 2.0*B0[1]*wv3-2.0*B2[1]*wv1+2.0*E1[1]; 
  alpha[2] = -0.5773502691896258*B2[0]*dv1; 
  alpha[3] = 0.5773502691896258*B0[0]*dv3; 
  alpha[4] = -0.5773502691896258*B2[1]*dv1; 
  alpha[5] = 0.5773502691896258*B0[1]*dv3; 

  double fUpwindQuad_l[8] = {0.0};
  double fUpwindQuad_r[8] = {0.0};
  double fUpwind_l[16] = {0.0};;
  double fUpwind_r[16] = {0.0};
  double Ghat_l[16] = {0.0}; 
  double Ghat_r[16] = {0.0}; 

  if (alpha[5]+alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad_l[0] = ser_4x_p1_surfx3_eval_quad_node_0_r(fl); 
    fUpwindQuad_r[0] = ser_4x_p1_surfx3_eval_quad_node_0_r(fc); 
  } else { 
    fUpwindQuad_l[0] = ser_4x_p1_surfx3_eval_quad_node_0_l(fc); 
    fUpwindQuad_r[0] = ser_4x_p1_surfx3_eval_quad_node_0_l(fr); 
  } 
  if ((-alpha[5])+alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad_l[1] = ser_4x_p1_surfx3_eval_quad_node_1_r(fl); 
    fUpwindQuad_r[1] = ser_4x_p1_surfx3_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_l[1] = ser_4x_p1_surfx3_eval_quad_node_1_l(fc); 
    fUpwindQuad_r[1] = ser_4x_p1_surfx3_eval_quad_node_1_l(fr); 
  } 
  if (alpha[5]-alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad_l[2] = ser_4x_p1_surfx3_eval_quad_node_2_r(fl); 
    fUpwindQuad_r[2] = ser_4x_p1_surfx3_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_l[2] = ser_4x_p1_surfx3_eval_quad_node_2_l(fc); 
    fUpwindQuad_r[2] = ser_4x_p1_surfx3_eval_quad_node_2_l(fr); 
  } 
  if ((-alpha[5])-alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad_l[3] = ser_4x_p1_surfx3_eval_quad_node_3_r(fl); 
    fUpwindQuad_r[3] = ser_4x_p1_surfx3_eval_quad_node_3_r(fc); 
  } else { 
    fUpwindQuad_l[3] = ser_4x_p1_surfx3_eval_quad_node_3_l(fc); 
    fUpwindQuad_r[3] = ser_4x_p1_surfx3_eval_quad_node_3_l(fr); 
  } 
  if ((-alpha[5])-alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad_l[4] = ser_4x_p1_surfx3_eval_quad_node_4_r(fl); 
    fUpwindQuad_r[4] = ser_4x_p1_surfx3_eval_quad_node_4_r(fc); 
  } else { 
    fUpwindQuad_l[4] = ser_4x_p1_surfx3_eval_quad_node_4_l(fc); 
    fUpwindQuad_r[4] = ser_4x_p1_surfx3_eval_quad_node_4_l(fr); 
  } 
  if (alpha[5]-alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad_l[5] = ser_4x_p1_surfx3_eval_quad_node_5_r(fl); 
    fUpwindQuad_r[5] = ser_4x_p1_surfx3_eval_quad_node_5_r(fc); 
  } else { 
    fUpwindQuad_l[5] = ser_4x_p1_surfx3_eval_quad_node_5_l(fc); 
    fUpwindQuad_r[5] = ser_4x_p1_surfx3_eval_quad_node_5_l(fr); 
  } 
  if ((-alpha[5])+alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad_l[6] = ser_4x_p1_surfx3_eval_quad_node_6_r(fl); 
    fUpwindQuad_r[6] = ser_4x_p1_surfx3_eval_quad_node_6_r(fc); 
  } else { 
    fUpwindQuad_l[6] = ser_4x_p1_surfx3_eval_quad_node_6_l(fc); 
    fUpwindQuad_r[6] = ser_4x_p1_surfx3_eval_quad_node_6_l(fr); 
  } 
  if (alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad_l[7] = ser_4x_p1_surfx3_eval_quad_node_7_r(fl); 
    fUpwindQuad_r[7] = ser_4x_p1_surfx3_eval_quad_node_7_r(fc); 
  } else { 
    fUpwindQuad_l[7] = ser_4x_p1_surfx3_eval_quad_node_7_l(fc); 
    fUpwindQuad_r[7] = ser_4x_p1_surfx3_eval_quad_node_7_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_4x_p1_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  ser_4x_p1_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.3535533905932737*(alpha[5]*fUpwind_l[5]+alpha[4]*fUpwind_l[4]+alpha[3]*fUpwind_l[3]+alpha[2]*fUpwind_l[2]+alpha[1]*fUpwind_l[1]+alpha[0]*fUpwind_l[0]); 
  Ghat_l[1] = 0.3535533905932737*(alpha[3]*fUpwind_l[5]+fUpwind_l[3]*alpha[5]+alpha[2]*fUpwind_l[4]+fUpwind_l[2]*alpha[4]+alpha[0]*fUpwind_l[1]+fUpwind_l[0]*alpha[1]); 
  Ghat_l[2] = 0.02357022603955158*(13.41640786499874*alpha[4]*fUpwind_l[9]+13.41640786499874*alpha[2]*fUpwind_l[8]+15.0*(alpha[5]*fUpwind_l[7]+alpha[3]*fUpwind_l[6]+alpha[1]*fUpwind_l[4]+fUpwind_l[1]*alpha[4]+alpha[0]*fUpwind_l[2]+fUpwind_l[0]*alpha[2])); 
  Ghat_l[3] = 0.02357022603955158*(13.41640786499874*alpha[5]*fUpwind_l[13]+13.41640786499874*alpha[3]*fUpwind_l[12]+15.0*(alpha[4]*fUpwind_l[7]+alpha[2]*fUpwind_l[6]+alpha[1]*fUpwind_l[5]+fUpwind_l[1]*alpha[5]+alpha[0]*fUpwind_l[3]+fUpwind_l[0]*alpha[3])); 
  Ghat_l[4] = 0.02357022603955158*(13.41640786499874*alpha[2]*fUpwind_l[9]+13.41640786499874*alpha[4]*fUpwind_l[8]+15.0*(alpha[3]*fUpwind_l[7]+alpha[5]*fUpwind_l[6]+alpha[0]*fUpwind_l[4]+fUpwind_l[0]*alpha[4]+alpha[1]*fUpwind_l[2]+fUpwind_l[1]*alpha[2])); 
  Ghat_l[5] = 0.02357022603955158*(13.41640786499874*alpha[3]*fUpwind_l[13]+13.41640786499874*alpha[5]*fUpwind_l[12]+15.0*(alpha[2]*fUpwind_l[7]+alpha[4]*fUpwind_l[6]+alpha[0]*fUpwind_l[5]+fUpwind_l[0]*alpha[5]+alpha[1]*fUpwind_l[3]+fUpwind_l[1]*alpha[3])); 
  Ghat_l[6] = 0.02357022603955158*(13.41640786499874*alpha[5]*fUpwind_l[15]+13.41640786499874*alpha[3]*fUpwind_l[14]+13.41640786499874*alpha[4]*fUpwind_l[11]+13.41640786499874*alpha[2]*fUpwind_l[10]+15.0*(alpha[1]*fUpwind_l[7]+alpha[0]*fUpwind_l[6]+alpha[4]*fUpwind_l[5]+fUpwind_l[4]*alpha[5]+alpha[2]*fUpwind_l[3]+fUpwind_l[2]*alpha[3])); 
  Ghat_l[7] = 0.02357022603955158*(13.41640786499874*alpha[3]*fUpwind_l[15]+13.41640786499874*alpha[5]*fUpwind_l[14]+13.41640786499874*alpha[2]*fUpwind_l[11]+13.41640786499874*alpha[4]*fUpwind_l[10]+15.0*(alpha[0]*fUpwind_l[7]+alpha[1]*fUpwind_l[6]+alpha[2]*fUpwind_l[5]+fUpwind_l[2]*alpha[5]+alpha[3]*fUpwind_l[4]+fUpwind_l[3]*alpha[4])); 
  Ghat_l[8] = 0.02357022603955158*(15.0*alpha[5]*fUpwind_l[11]+15.0*(alpha[3]*fUpwind_l[10]+alpha[1]*fUpwind_l[9])+15.0*alpha[0]*fUpwind_l[8]+13.41640786499874*(alpha[4]*fUpwind_l[4]+alpha[2]*fUpwind_l[2])); 
  Ghat_l[9] = 0.02357022603955158*(15.0*alpha[3]*fUpwind_l[11]+15.0*(alpha[5]*fUpwind_l[10]+alpha[0]*fUpwind_l[9])+15.0*alpha[1]*fUpwind_l[8]+13.41640786499874*(alpha[2]*fUpwind_l[4]+fUpwind_l[2]*alpha[4])); 
  Ghat_l[10] = 0.02357022603955158*(15.0*alpha[1]*fUpwind_l[11]+15.0*(alpha[0]*fUpwind_l[10]+alpha[5]*fUpwind_l[9])+15.0*alpha[3]*fUpwind_l[8]+13.41640786499874*(alpha[4]*fUpwind_l[7]+alpha[2]*fUpwind_l[6])); 
  Ghat_l[11] = 0.02357022603955158*(15.0*alpha[0]*fUpwind_l[11]+15.0*(alpha[1]*fUpwind_l[10]+alpha[3]*fUpwind_l[9])+15.0*alpha[5]*fUpwind_l[8]+13.41640786499874*(alpha[2]*fUpwind_l[7]+alpha[4]*fUpwind_l[6])); 
  Ghat_l[12] = 0.02357022603955158*(15.0*alpha[4]*fUpwind_l[15]+15.0*(alpha[2]*fUpwind_l[14]+alpha[1]*fUpwind_l[13])+15.0*alpha[0]*fUpwind_l[12]+13.41640786499874*(alpha[5]*fUpwind_l[5]+alpha[3]*fUpwind_l[3])); 
  Ghat_l[13] = 0.02357022603955158*(15.0*alpha[2]*fUpwind_l[15]+15.0*(alpha[4]*fUpwind_l[14]+alpha[0]*fUpwind_l[13])+15.0*alpha[1]*fUpwind_l[12]+13.41640786499874*(alpha[3]*fUpwind_l[5]+fUpwind_l[3]*alpha[5])); 
  Ghat_l[14] = 0.02357022603955158*(15.0*alpha[1]*fUpwind_l[15]+15.0*(alpha[0]*fUpwind_l[14]+alpha[4]*fUpwind_l[13])+15.0*alpha[2]*fUpwind_l[12]+13.41640786499874*(alpha[5]*fUpwind_l[7]+alpha[3]*fUpwind_l[6])); 
  Ghat_l[15] = 0.02357022603955158*(15.0*alpha[0]*fUpwind_l[15]+15.0*(alpha[1]*fUpwind_l[14]+alpha[2]*fUpwind_l[13])+15.0*alpha[4]*fUpwind_l[12]+13.41640786499874*(alpha[3]*fUpwind_l[7]+alpha[5]*fUpwind_l[6])); 

  Ghat_r[0] = 0.3535533905932737*(alpha[5]*fUpwind_r[5]+alpha[4]*fUpwind_r[4]+alpha[3]*fUpwind_r[3]+alpha[2]*fUpwind_r[2]+alpha[1]*fUpwind_r[1]+alpha[0]*fUpwind_r[0]); 
  Ghat_r[1] = 0.3535533905932737*(alpha[3]*fUpwind_r[5]+fUpwind_r[3]*alpha[5]+alpha[2]*fUpwind_r[4]+fUpwind_r[2]*alpha[4]+alpha[0]*fUpwind_r[1]+fUpwind_r[0]*alpha[1]); 
  Ghat_r[2] = 0.02357022603955158*(13.41640786499874*alpha[4]*fUpwind_r[9]+13.41640786499874*alpha[2]*fUpwind_r[8]+15.0*(alpha[5]*fUpwind_r[7]+alpha[3]*fUpwind_r[6]+alpha[1]*fUpwind_r[4]+fUpwind_r[1]*alpha[4]+alpha[0]*fUpwind_r[2]+fUpwind_r[0]*alpha[2])); 
  Ghat_r[3] = 0.02357022603955158*(13.41640786499874*alpha[5]*fUpwind_r[13]+13.41640786499874*alpha[3]*fUpwind_r[12]+15.0*(alpha[4]*fUpwind_r[7]+alpha[2]*fUpwind_r[6]+alpha[1]*fUpwind_r[5]+fUpwind_r[1]*alpha[5]+alpha[0]*fUpwind_r[3]+fUpwind_r[0]*alpha[3])); 
  Ghat_r[4] = 0.02357022603955158*(13.41640786499874*alpha[2]*fUpwind_r[9]+13.41640786499874*alpha[4]*fUpwind_r[8]+15.0*(alpha[3]*fUpwind_r[7]+alpha[5]*fUpwind_r[6]+alpha[0]*fUpwind_r[4]+fUpwind_r[0]*alpha[4]+alpha[1]*fUpwind_r[2]+fUpwind_r[1]*alpha[2])); 
  Ghat_r[5] = 0.02357022603955158*(13.41640786499874*alpha[3]*fUpwind_r[13]+13.41640786499874*alpha[5]*fUpwind_r[12]+15.0*(alpha[2]*fUpwind_r[7]+alpha[4]*fUpwind_r[6]+alpha[0]*fUpwind_r[5]+fUpwind_r[0]*alpha[5]+alpha[1]*fUpwind_r[3]+fUpwind_r[1]*alpha[3])); 
  Ghat_r[6] = 0.02357022603955158*(13.41640786499874*alpha[5]*fUpwind_r[15]+13.41640786499874*alpha[3]*fUpwind_r[14]+13.41640786499874*alpha[4]*fUpwind_r[11]+13.41640786499874*alpha[2]*fUpwind_r[10]+15.0*(alpha[1]*fUpwind_r[7]+alpha[0]*fUpwind_r[6]+alpha[4]*fUpwind_r[5]+fUpwind_r[4]*alpha[5]+alpha[2]*fUpwind_r[3]+fUpwind_r[2]*alpha[3])); 
  Ghat_r[7] = 0.02357022603955158*(13.41640786499874*alpha[3]*fUpwind_r[15]+13.41640786499874*alpha[5]*fUpwind_r[14]+13.41640786499874*alpha[2]*fUpwind_r[11]+13.41640786499874*alpha[4]*fUpwind_r[10]+15.0*(alpha[0]*fUpwind_r[7]+alpha[1]*fUpwind_r[6]+alpha[2]*fUpwind_r[5]+fUpwind_r[2]*alpha[5]+alpha[3]*fUpwind_r[4]+fUpwind_r[3]*alpha[4])); 
  Ghat_r[8] = 0.02357022603955158*(15.0*alpha[5]*fUpwind_r[11]+15.0*(alpha[3]*fUpwind_r[10]+alpha[1]*fUpwind_r[9])+15.0*alpha[0]*fUpwind_r[8]+13.41640786499874*(alpha[4]*fUpwind_r[4]+alpha[2]*fUpwind_r[2])); 
  Ghat_r[9] = 0.02357022603955158*(15.0*alpha[3]*fUpwind_r[11]+15.0*(alpha[5]*fUpwind_r[10]+alpha[0]*fUpwind_r[9])+15.0*alpha[1]*fUpwind_r[8]+13.41640786499874*(alpha[2]*fUpwind_r[4]+fUpwind_r[2]*alpha[4])); 
  Ghat_r[10] = 0.02357022603955158*(15.0*alpha[1]*fUpwind_r[11]+15.0*(alpha[0]*fUpwind_r[10]+alpha[5]*fUpwind_r[9])+15.0*alpha[3]*fUpwind_r[8]+13.41640786499874*(alpha[4]*fUpwind_r[7]+alpha[2]*fUpwind_r[6])); 
  Ghat_r[11] = 0.02357022603955158*(15.0*alpha[0]*fUpwind_r[11]+15.0*(alpha[1]*fUpwind_r[10]+alpha[3]*fUpwind_r[9])+15.0*alpha[5]*fUpwind_r[8]+13.41640786499874*(alpha[2]*fUpwind_r[7]+alpha[4]*fUpwind_r[6])); 
  Ghat_r[12] = 0.02357022603955158*(15.0*alpha[4]*fUpwind_r[15]+15.0*(alpha[2]*fUpwind_r[14]+alpha[1]*fUpwind_r[13])+15.0*alpha[0]*fUpwind_r[12]+13.41640786499874*(alpha[5]*fUpwind_r[5]+alpha[3]*fUpwind_r[3])); 
  Ghat_r[13] = 0.02357022603955158*(15.0*alpha[2]*fUpwind_r[15]+15.0*(alpha[4]*fUpwind_r[14]+alpha[0]*fUpwind_r[13])+15.0*alpha[1]*fUpwind_r[12]+13.41640786499874*(alpha[3]*fUpwind_r[5]+fUpwind_r[3]*alpha[5])); 
  Ghat_r[14] = 0.02357022603955158*(15.0*alpha[1]*fUpwind_r[15]+15.0*(alpha[0]*fUpwind_r[14]+alpha[4]*fUpwind_r[13])+15.0*alpha[2]*fUpwind_r[12]+13.41640786499874*(alpha[5]*fUpwind_r[7]+alpha[3]*fUpwind_r[6])); 
  Ghat_r[15] = 0.02357022603955158*(15.0*alpha[0]*fUpwind_r[15]+15.0*(alpha[1]*fUpwind_r[14]+alpha[2]*fUpwind_r[13])+15.0*alpha[4]*fUpwind_r[12]+13.41640786499874*(alpha[3]*fUpwind_r[7]+alpha[5]*fUpwind_r[6])); 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv11; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv11; 
  out[2] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv11; 
  out[3] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv11; 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv11; 
  out[5] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dv11; 
  out[6] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv11; 
  out[7] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv11; 
  out[8] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dv11; 
  out[9] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dv11; 
  out[10] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv11; 
  out[11] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dv11; 
  out[12] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dv11; 
  out[13] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dv11; 
  out[14] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dv11; 
  out[15] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dv11; 
  out[16] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*dv11; 
  out[17] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*dv11; 
  out[18] += -1.224744871391589*(Ghat_r[8]+Ghat_l[8])*dv11; 
  out[19] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*dv11; 
  out[20] += -1.224744871391589*(Ghat_r[9]+Ghat_l[9])*dv11; 
  out[21] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*dv11; 
  out[22] += -1.224744871391589*(Ghat_r[10]+Ghat_l[10])*dv11; 
  out[23] += -1.224744871391589*(Ghat_r[11]+Ghat_l[11])*dv11; 
  out[24] += (1.58113883008419*Ghat_l[0]-1.58113883008419*Ghat_r[0])*dv11; 
  out[25] += (1.58113883008419*Ghat_l[1]-1.58113883008419*Ghat_r[1])*dv11; 
  out[26] += (1.58113883008419*Ghat_l[2]-1.58113883008419*Ghat_r[2])*dv11; 
  out[27] += (1.58113883008419*Ghat_l[3]-1.58113883008419*Ghat_r[3])*dv11; 
  out[28] += (1.58113883008419*Ghat_l[4]-1.58113883008419*Ghat_r[4])*dv11; 
  out[29] += (1.58113883008419*Ghat_l[5]-1.58113883008419*Ghat_r[5])*dv11; 
  out[30] += (1.58113883008419*Ghat_l[6]-1.58113883008419*Ghat_r[6])*dv11; 
  out[31] += (1.58113883008419*Ghat_l[7]-1.58113883008419*Ghat_r[7])*dv11; 
  out[32] += (0.7071067811865475*Ghat_l[12]-0.7071067811865475*Ghat_r[12])*dv11; 
  out[33] += (0.7071067811865475*Ghat_l[13]-0.7071067811865475*Ghat_r[13])*dv11; 
  out[34] += (0.7071067811865475*Ghat_l[14]-0.7071067811865475*Ghat_r[14])*dv11; 
  out[35] += -1.224744871391589*(Ghat_r[12]+Ghat_l[12])*dv11; 
  out[36] += (0.7071067811865475*Ghat_l[15]-0.7071067811865475*Ghat_r[15])*dv11; 
  out[37] += -1.224744871391589*(Ghat_r[13]+Ghat_l[13])*dv11; 
  out[38] += -1.224744871391589*(Ghat_r[14]+Ghat_l[14])*dv11; 
  out[39] += -1.224744871391589*(Ghat_r[15]+Ghat_l[15])*dv11; 

} 
