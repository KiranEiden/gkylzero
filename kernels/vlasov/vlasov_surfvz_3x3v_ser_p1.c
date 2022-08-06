#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_hyb_3x3v_p1_surfx6_eval_quad.h> 
#include <gkyl_basis_hyb_3x3v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_surfvz_3x3v_ser_p1(const double *w, const double *dxv, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // qmem:      q/m*EM fields.
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 
  const double dv12 = 2/dxv[5]; 
  const double dv1 = dxv[3], wv1 = w[3]; 
  const double dv2 = dxv[4], wv2 = w[4]; 
  const double dv3 = dxv[5], wv3 = w[5]; 
  const double *E2 = &qmem[16]; 
  const double *B0 = &qmem[24]; 
  const double *B1 = &qmem[32]; 
  const double *B2 = &qmem[40]; 

  double alpha[64] = {0.0}; 

  alpha[0] = (-2.0*B0[0]*wv2)+2.0*B1[0]*wv1+2.0*E2[0]; 
  alpha[1] = (-2.0*B0[1]*wv2)+2.0*B1[1]*wv1+2.0*E2[1]; 
  alpha[2] = (-2.0*B0[2]*wv2)+2.0*B1[2]*wv1+2.0*E2[2]; 
  alpha[3] = (-2.0*B0[3]*wv2)+2.0*B1[3]*wv1+2.0*E2[3]; 
  alpha[4] = 0.5773502691896258*B1[0]*dv1; 
  alpha[5] = -0.5773502691896258*B0[0]*dv2; 
  alpha[6] = (-2.0*B0[4]*wv2)+2.0*B1[4]*wv1+2.0*E2[4]; 
  alpha[7] = (-2.0*B0[5]*wv2)+2.0*B1[5]*wv1+2.0*E2[5]; 
  alpha[8] = (-2.0*B0[6]*wv2)+2.0*B1[6]*wv1+2.0*E2[6]; 
  alpha[9] = 0.5773502691896258*B1[1]*dv1; 
  alpha[10] = 0.5773502691896258*B1[2]*dv1; 
  alpha[11] = 0.5773502691896258*B1[3]*dv1; 
  alpha[12] = -0.5773502691896258*B0[1]*dv2; 
  alpha[13] = -0.5773502691896258*B0[2]*dv2; 
  alpha[14] = -0.5773502691896258*B0[3]*dv2; 
  alpha[16] = (-2.0*B0[7]*wv2)+2.0*B1[7]*wv1+2.0*E2[7]; 
  alpha[17] = 0.5773502691896258*B1[4]*dv1; 
  alpha[18] = 0.5773502691896258*B1[5]*dv1; 
  alpha[19] = 0.5773502691896258*B1[6]*dv1; 
  alpha[20] = -0.5773502691896258*B0[4]*dv2; 
  alpha[21] = -0.5773502691896258*B0[5]*dv2; 
  alpha[22] = -0.5773502691896258*B0[6]*dv2; 
  alpha[26] = 0.5773502691896258*B1[7]*dv1; 
  alpha[27] = -0.5773502691896258*B0[7]*dv2; 

  double fUpwindQuad_l[72] = {0.0};
  double fUpwindQuad_r[72] = {0.0};
  double fUpwind_l[64] = {0.0};;
  double fUpwind_r[64] = {0.0};
  double Ghat_l[64] = {0.0}; 
  double Ghat_r[64] = {0.0}; 

  if (188887112823866*alpha[27]+188887112823866*alpha[26]-188887112823866*alpha[22]-188887112823866*alpha[21]-188887112823866*alpha[20]-188887112823866*alpha[19]-188887112823866*alpha[18]-188887112823866*alpha[17]-140788141449279*alpha[16]+188887112823866*alpha[14]+188887112823866*alpha[13]+188887112823866*alpha[12]+188887112823866*alpha[11]+188887112823866*alpha[10]+188887112823866*alpha[9]+140788141449279*alpha[8]+140788141449279*alpha[7]+140788141449279*alpha[6]-188887112823866*alpha[5]-188887112823866*alpha[4]-140788141449279*alpha[3]-140788141449279*alpha[2]-140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[0] = hyb_3x3v_p1_surfx6_eval_quad_node_0_r(fl); 
    fUpwindQuad_r[0] = hyb_3x3v_p1_surfx6_eval_quad_node_0_r(fc); 
  } else { 
    fUpwindQuad_l[0] = hyb_3x3v_p1_surfx6_eval_quad_node_0_l(fc); 
    fUpwindQuad_r[0] = hyb_3x3v_p1_surfx6_eval_quad_node_0_l(fr); 
  } 
  if (188887112823866*alpha[26]-188887112823866*alpha[19]-188887112823866*alpha[18]-188887112823866*alpha[17]-140788141449279*alpha[16]+188887112823866*alpha[11]+188887112823866*alpha[10]+188887112823866*alpha[9]+140788141449279*alpha[8]+140788141449279*alpha[7]+140788141449279*alpha[6]-188887112823866*alpha[4]-140788141449279*alpha[3]-140788141449279*alpha[2]-140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[1] = hyb_3x3v_p1_surfx6_eval_quad_node_1_r(fl); 
    fUpwindQuad_r[1] = hyb_3x3v_p1_surfx6_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_l[1] = hyb_3x3v_p1_surfx6_eval_quad_node_1_l(fc); 
    fUpwindQuad_r[1] = hyb_3x3v_p1_surfx6_eval_quad_node_1_l(fr); 
  } 
  if ((-188887112823866*alpha[27])+188887112823866*alpha[26]+188887112823866*alpha[22]+188887112823866*alpha[21]+188887112823866*alpha[20]-188887112823866*alpha[19]-188887112823866*alpha[18]-188887112823866*alpha[17]-140788141449279*alpha[16]-188887112823866*alpha[14]-188887112823866*alpha[13]-188887112823866*alpha[12]+188887112823866*alpha[11]+188887112823866*alpha[10]+188887112823866*alpha[9]+140788141449279*alpha[8]+140788141449279*alpha[7]+140788141449279*alpha[6]+188887112823866*alpha[5]-188887112823866*alpha[4]-140788141449279*alpha[3]-140788141449279*alpha[2]-140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[2] = hyb_3x3v_p1_surfx6_eval_quad_node_2_r(fl); 
    fUpwindQuad_r[2] = hyb_3x3v_p1_surfx6_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_l[2] = hyb_3x3v_p1_surfx6_eval_quad_node_2_l(fc); 
    fUpwindQuad_r[2] = hyb_3x3v_p1_surfx6_eval_quad_node_2_l(fr); 
  } 
  if (188887112823866*alpha[27]-188887112823866*alpha[22]-188887112823866*alpha[21]-188887112823866*alpha[20]-140788141449279*alpha[16]+188887112823866*alpha[14]+188887112823866*alpha[13]+188887112823866*alpha[12]+140788141449279*alpha[8]+140788141449279*alpha[7]+140788141449279*alpha[6]-188887112823866*alpha[5]-140788141449279*alpha[3]-140788141449279*alpha[2]-140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[3] = hyb_3x3v_p1_surfx6_eval_quad_node_3_r(fl); 
    fUpwindQuad_r[3] = hyb_3x3v_p1_surfx6_eval_quad_node_3_r(fc); 
  } else { 
    fUpwindQuad_l[3] = hyb_3x3v_p1_surfx6_eval_quad_node_3_l(fc); 
    fUpwindQuad_r[3] = hyb_3x3v_p1_surfx6_eval_quad_node_3_l(fr); 
  } 
  if ((-140788141449279*alpha[16])+140788141449279*alpha[8]+140788141449279*alpha[7]+140788141449279*alpha[6]-140788141449279*alpha[3]-140788141449279*alpha[2]-140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[4] = hyb_3x3v_p1_surfx6_eval_quad_node_4_r(fl); 
    fUpwindQuad_r[4] = hyb_3x3v_p1_surfx6_eval_quad_node_4_r(fc); 
  } else { 
    fUpwindQuad_l[4] = hyb_3x3v_p1_surfx6_eval_quad_node_4_l(fc); 
    fUpwindQuad_r[4] = hyb_3x3v_p1_surfx6_eval_quad_node_4_l(fr); 
  } 
  if ((-188887112823866*alpha[27])+188887112823866*alpha[22]+188887112823866*alpha[21]+188887112823866*alpha[20]-140788141449279*alpha[16]-188887112823866*alpha[14]-188887112823866*alpha[13]-188887112823866*alpha[12]+140788141449279*alpha[8]+140788141449279*alpha[7]+140788141449279*alpha[6]+188887112823866*alpha[5]-140788141449279*alpha[3]-140788141449279*alpha[2]-140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[5] = hyb_3x3v_p1_surfx6_eval_quad_node_5_r(fl); 
    fUpwindQuad_r[5] = hyb_3x3v_p1_surfx6_eval_quad_node_5_r(fc); 
  } else { 
    fUpwindQuad_l[5] = hyb_3x3v_p1_surfx6_eval_quad_node_5_l(fc); 
    fUpwindQuad_r[5] = hyb_3x3v_p1_surfx6_eval_quad_node_5_l(fr); 
  } 
  if (188887112823866*alpha[27]-188887112823866*alpha[26]-188887112823866*alpha[22]-188887112823866*alpha[21]-188887112823866*alpha[20]+188887112823866*alpha[19]+188887112823866*alpha[18]+188887112823866*alpha[17]-140788141449279*alpha[16]+188887112823866*alpha[14]+188887112823866*alpha[13]+188887112823866*alpha[12]-188887112823866*alpha[11]-188887112823866*alpha[10]-188887112823866*alpha[9]+140788141449279*alpha[8]+140788141449279*alpha[7]+140788141449279*alpha[6]-188887112823866*alpha[5]+188887112823866*alpha[4]-140788141449279*alpha[3]-140788141449279*alpha[2]-140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[6] = hyb_3x3v_p1_surfx6_eval_quad_node_6_r(fl); 
    fUpwindQuad_r[6] = hyb_3x3v_p1_surfx6_eval_quad_node_6_r(fc); 
  } else { 
    fUpwindQuad_l[6] = hyb_3x3v_p1_surfx6_eval_quad_node_6_l(fc); 
    fUpwindQuad_r[6] = hyb_3x3v_p1_surfx6_eval_quad_node_6_l(fr); 
  } 
  if ((-188887112823866*alpha[26])+188887112823866*alpha[19]+188887112823866*alpha[18]+188887112823866*alpha[17]-140788141449279*alpha[16]-188887112823866*alpha[11]-188887112823866*alpha[10]-188887112823866*alpha[9]+140788141449279*alpha[8]+140788141449279*alpha[7]+140788141449279*alpha[6]+188887112823866*alpha[4]-140788141449279*alpha[3]-140788141449279*alpha[2]-140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[7] = hyb_3x3v_p1_surfx6_eval_quad_node_7_r(fl); 
    fUpwindQuad_r[7] = hyb_3x3v_p1_surfx6_eval_quad_node_7_r(fc); 
  } else { 
    fUpwindQuad_l[7] = hyb_3x3v_p1_surfx6_eval_quad_node_7_l(fc); 
    fUpwindQuad_r[7] = hyb_3x3v_p1_surfx6_eval_quad_node_7_l(fr); 
  } 
  if ((-188887112823866*alpha[27])-188887112823866*alpha[26]+188887112823866*alpha[22]+188887112823866*alpha[21]+188887112823866*alpha[20]+188887112823866*alpha[19]+188887112823866*alpha[18]+188887112823866*alpha[17]-140788141449279*alpha[16]-188887112823866*alpha[14]-188887112823866*alpha[13]-188887112823866*alpha[12]-188887112823866*alpha[11]-188887112823866*alpha[10]-188887112823866*alpha[9]+140788141449279*alpha[8]+140788141449279*alpha[7]+140788141449279*alpha[6]+188887112823866*alpha[5]+188887112823866*alpha[4]-140788141449279*alpha[3]-140788141449279*alpha[2]-140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[8] = hyb_3x3v_p1_surfx6_eval_quad_node_8_r(fl); 
    fUpwindQuad_r[8] = hyb_3x3v_p1_surfx6_eval_quad_node_8_r(fc); 
  } else { 
    fUpwindQuad_l[8] = hyb_3x3v_p1_surfx6_eval_quad_node_8_l(fc); 
    fUpwindQuad_r[8] = hyb_3x3v_p1_surfx6_eval_quad_node_8_l(fr); 
  } 
  if ((-188887112823866*alpha[27])-188887112823866*alpha[26]+188887112823866*alpha[22]+188887112823866*alpha[21]-188887112823866*alpha[20]+188887112823866*alpha[19]+188887112823866*alpha[18]-188887112823866*alpha[17]+140788141449279*alpha[16]-188887112823866*alpha[14]+188887112823866*alpha[13]+188887112823866*alpha[12]-188887112823866*alpha[11]+188887112823866*alpha[10]+188887112823866*alpha[9]-140788141449279*alpha[8]-140788141449279*alpha[7]+140788141449279*alpha[6]-188887112823866*alpha[5]-188887112823866*alpha[4]+140788141449279*alpha[3]-140788141449279*alpha[2]-140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[9] = hyb_3x3v_p1_surfx6_eval_quad_node_9_r(fl); 
    fUpwindQuad_r[9] = hyb_3x3v_p1_surfx6_eval_quad_node_9_r(fc); 
  } else { 
    fUpwindQuad_l[9] = hyb_3x3v_p1_surfx6_eval_quad_node_9_l(fc); 
    fUpwindQuad_r[9] = hyb_3x3v_p1_surfx6_eval_quad_node_9_l(fr); 
  } 
  if ((-188887112823866*alpha[26])+188887112823866*alpha[19]+188887112823866*alpha[18]-188887112823866*alpha[17]+140788141449279*alpha[16]-188887112823866*alpha[11]+188887112823866*alpha[10]+188887112823866*alpha[9]-140788141449279*alpha[8]-140788141449279*alpha[7]+140788141449279*alpha[6]-188887112823866*alpha[4]+140788141449279*alpha[3]-140788141449279*alpha[2]-140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[10] = hyb_3x3v_p1_surfx6_eval_quad_node_10_r(fl); 
    fUpwindQuad_r[10] = hyb_3x3v_p1_surfx6_eval_quad_node_10_r(fc); 
  } else { 
    fUpwindQuad_l[10] = hyb_3x3v_p1_surfx6_eval_quad_node_10_l(fc); 
    fUpwindQuad_r[10] = hyb_3x3v_p1_surfx6_eval_quad_node_10_l(fr); 
  } 
  if (188887112823866*alpha[27]-188887112823866*alpha[26]-188887112823866*alpha[22]-188887112823866*alpha[21]+188887112823866*alpha[20]+188887112823866*alpha[19]+188887112823866*alpha[18]-188887112823866*alpha[17]+140788141449279*alpha[16]+188887112823866*alpha[14]-188887112823866*alpha[13]-188887112823866*alpha[12]-188887112823866*alpha[11]+188887112823866*alpha[10]+188887112823866*alpha[9]-140788141449279*alpha[8]-140788141449279*alpha[7]+140788141449279*alpha[6]+188887112823866*alpha[5]-188887112823866*alpha[4]+140788141449279*alpha[3]-140788141449279*alpha[2]-140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[11] = hyb_3x3v_p1_surfx6_eval_quad_node_11_r(fl); 
    fUpwindQuad_r[11] = hyb_3x3v_p1_surfx6_eval_quad_node_11_r(fc); 
  } else { 
    fUpwindQuad_l[11] = hyb_3x3v_p1_surfx6_eval_quad_node_11_l(fc); 
    fUpwindQuad_r[11] = hyb_3x3v_p1_surfx6_eval_quad_node_11_l(fr); 
  } 
  if ((-188887112823866*alpha[27])+188887112823866*alpha[22]+188887112823866*alpha[21]-188887112823866*alpha[20]+140788141449279*alpha[16]-188887112823866*alpha[14]+188887112823866*alpha[13]+188887112823866*alpha[12]-140788141449279*alpha[8]-140788141449279*alpha[7]+140788141449279*alpha[6]-188887112823866*alpha[5]+140788141449279*alpha[3]-140788141449279*alpha[2]-140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[12] = hyb_3x3v_p1_surfx6_eval_quad_node_12_r(fl); 
    fUpwindQuad_r[12] = hyb_3x3v_p1_surfx6_eval_quad_node_12_r(fc); 
  } else { 
    fUpwindQuad_l[12] = hyb_3x3v_p1_surfx6_eval_quad_node_12_l(fc); 
    fUpwindQuad_r[12] = hyb_3x3v_p1_surfx6_eval_quad_node_12_l(fr); 
  } 
  if (140788141449279*alpha[16]-140788141449279*alpha[8]-140788141449279*alpha[7]+140788141449279*alpha[6]+140788141449279*alpha[3]-140788141449279*alpha[2]-140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[13] = hyb_3x3v_p1_surfx6_eval_quad_node_13_r(fl); 
    fUpwindQuad_r[13] = hyb_3x3v_p1_surfx6_eval_quad_node_13_r(fc); 
  } else { 
    fUpwindQuad_l[13] = hyb_3x3v_p1_surfx6_eval_quad_node_13_l(fc); 
    fUpwindQuad_r[13] = hyb_3x3v_p1_surfx6_eval_quad_node_13_l(fr); 
  } 
  if (188887112823866*alpha[27]-188887112823866*alpha[22]-188887112823866*alpha[21]+188887112823866*alpha[20]+140788141449279*alpha[16]+188887112823866*alpha[14]-188887112823866*alpha[13]-188887112823866*alpha[12]-140788141449279*alpha[8]-140788141449279*alpha[7]+140788141449279*alpha[6]+188887112823866*alpha[5]+140788141449279*alpha[3]-140788141449279*alpha[2]-140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[14] = hyb_3x3v_p1_surfx6_eval_quad_node_14_r(fl); 
    fUpwindQuad_r[14] = hyb_3x3v_p1_surfx6_eval_quad_node_14_r(fc); 
  } else { 
    fUpwindQuad_l[14] = hyb_3x3v_p1_surfx6_eval_quad_node_14_l(fc); 
    fUpwindQuad_r[14] = hyb_3x3v_p1_surfx6_eval_quad_node_14_l(fr); 
  } 
  if ((-188887112823866*alpha[27])+188887112823866*alpha[26]+188887112823866*alpha[22]+188887112823866*alpha[21]-188887112823866*alpha[20]-188887112823866*alpha[19]-188887112823866*alpha[18]+188887112823866*alpha[17]+140788141449279*alpha[16]-188887112823866*alpha[14]+188887112823866*alpha[13]+188887112823866*alpha[12]+188887112823866*alpha[11]-188887112823866*alpha[10]-188887112823866*alpha[9]-140788141449279*alpha[8]-140788141449279*alpha[7]+140788141449279*alpha[6]-188887112823866*alpha[5]+188887112823866*alpha[4]+140788141449279*alpha[3]-140788141449279*alpha[2]-140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[15] = hyb_3x3v_p1_surfx6_eval_quad_node_15_r(fl); 
    fUpwindQuad_r[15] = hyb_3x3v_p1_surfx6_eval_quad_node_15_r(fc); 
  } else { 
    fUpwindQuad_l[15] = hyb_3x3v_p1_surfx6_eval_quad_node_15_l(fc); 
    fUpwindQuad_r[15] = hyb_3x3v_p1_surfx6_eval_quad_node_15_l(fr); 
  } 
  if (188887112823866*alpha[26]-188887112823866*alpha[19]-188887112823866*alpha[18]+188887112823866*alpha[17]+140788141449279*alpha[16]+188887112823866*alpha[11]-188887112823866*alpha[10]-188887112823866*alpha[9]-140788141449279*alpha[8]-140788141449279*alpha[7]+140788141449279*alpha[6]+188887112823866*alpha[4]+140788141449279*alpha[3]-140788141449279*alpha[2]-140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[16] = hyb_3x3v_p1_surfx6_eval_quad_node_16_r(fl); 
    fUpwindQuad_r[16] = hyb_3x3v_p1_surfx6_eval_quad_node_16_r(fc); 
  } else { 
    fUpwindQuad_l[16] = hyb_3x3v_p1_surfx6_eval_quad_node_16_l(fc); 
    fUpwindQuad_r[16] = hyb_3x3v_p1_surfx6_eval_quad_node_16_l(fr); 
  } 
  if (188887112823866*alpha[27]+188887112823866*alpha[26]-188887112823866*alpha[22]-188887112823866*alpha[21]+188887112823866*alpha[20]-188887112823866*alpha[19]-188887112823866*alpha[18]+188887112823866*alpha[17]+140788141449279*alpha[16]+188887112823866*alpha[14]-188887112823866*alpha[13]-188887112823866*alpha[12]+188887112823866*alpha[11]-188887112823866*alpha[10]-188887112823866*alpha[9]-140788141449279*alpha[8]-140788141449279*alpha[7]+140788141449279*alpha[6]+188887112823866*alpha[5]+188887112823866*alpha[4]+140788141449279*alpha[3]-140788141449279*alpha[2]-140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[17] = hyb_3x3v_p1_surfx6_eval_quad_node_17_r(fl); 
    fUpwindQuad_r[17] = hyb_3x3v_p1_surfx6_eval_quad_node_17_r(fc); 
  } else { 
    fUpwindQuad_l[17] = hyb_3x3v_p1_surfx6_eval_quad_node_17_l(fc); 
    fUpwindQuad_r[17] = hyb_3x3v_p1_surfx6_eval_quad_node_17_l(fr); 
  } 
  if ((-188887112823866*alpha[27])-188887112823866*alpha[26]+188887112823866*alpha[22]-188887112823866*alpha[21]+188887112823866*alpha[20]+188887112823866*alpha[19]-188887112823866*alpha[18]+188887112823866*alpha[17]+140788141449279*alpha[16]+188887112823866*alpha[14]-188887112823866*alpha[13]+188887112823866*alpha[12]+188887112823866*alpha[11]-188887112823866*alpha[10]+188887112823866*alpha[9]-140788141449279*alpha[8]+140788141449279*alpha[7]-140788141449279*alpha[6]-188887112823866*alpha[5]-188887112823866*alpha[4]-140788141449279*alpha[3]+140788141449279*alpha[2]-140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[18] = hyb_3x3v_p1_surfx6_eval_quad_node_18_r(fl); 
    fUpwindQuad_r[18] = hyb_3x3v_p1_surfx6_eval_quad_node_18_r(fc); 
  } else { 
    fUpwindQuad_l[18] = hyb_3x3v_p1_surfx6_eval_quad_node_18_l(fc); 
    fUpwindQuad_r[18] = hyb_3x3v_p1_surfx6_eval_quad_node_18_l(fr); 
  } 
  if ((-188887112823866*alpha[26])+188887112823866*alpha[19]-188887112823866*alpha[18]+188887112823866*alpha[17]+140788141449279*alpha[16]+188887112823866*alpha[11]-188887112823866*alpha[10]+188887112823866*alpha[9]-140788141449279*alpha[8]+140788141449279*alpha[7]-140788141449279*alpha[6]-188887112823866*alpha[4]-140788141449279*alpha[3]+140788141449279*alpha[2]-140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[19] = hyb_3x3v_p1_surfx6_eval_quad_node_19_r(fl); 
    fUpwindQuad_r[19] = hyb_3x3v_p1_surfx6_eval_quad_node_19_r(fc); 
  } else { 
    fUpwindQuad_l[19] = hyb_3x3v_p1_surfx6_eval_quad_node_19_l(fc); 
    fUpwindQuad_r[19] = hyb_3x3v_p1_surfx6_eval_quad_node_19_l(fr); 
  } 
  if (188887112823866*alpha[27]-188887112823866*alpha[26]-188887112823866*alpha[22]+188887112823866*alpha[21]-188887112823866*alpha[20]+188887112823866*alpha[19]-188887112823866*alpha[18]+188887112823866*alpha[17]+140788141449279*alpha[16]-188887112823866*alpha[14]+188887112823866*alpha[13]-188887112823866*alpha[12]+188887112823866*alpha[11]-188887112823866*alpha[10]+188887112823866*alpha[9]-140788141449279*alpha[8]+140788141449279*alpha[7]-140788141449279*alpha[6]+188887112823866*alpha[5]-188887112823866*alpha[4]-140788141449279*alpha[3]+140788141449279*alpha[2]-140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[20] = hyb_3x3v_p1_surfx6_eval_quad_node_20_r(fl); 
    fUpwindQuad_r[20] = hyb_3x3v_p1_surfx6_eval_quad_node_20_r(fc); 
  } else { 
    fUpwindQuad_l[20] = hyb_3x3v_p1_surfx6_eval_quad_node_20_l(fc); 
    fUpwindQuad_r[20] = hyb_3x3v_p1_surfx6_eval_quad_node_20_l(fr); 
  } 
  if ((-188887112823866*alpha[27])+188887112823866*alpha[22]-188887112823866*alpha[21]+188887112823866*alpha[20]+140788141449279*alpha[16]+188887112823866*alpha[14]-188887112823866*alpha[13]+188887112823866*alpha[12]-140788141449279*alpha[8]+140788141449279*alpha[7]-140788141449279*alpha[6]-188887112823866*alpha[5]-140788141449279*alpha[3]+140788141449279*alpha[2]-140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[21] = hyb_3x3v_p1_surfx6_eval_quad_node_21_r(fl); 
    fUpwindQuad_r[21] = hyb_3x3v_p1_surfx6_eval_quad_node_21_r(fc); 
  } else { 
    fUpwindQuad_l[21] = hyb_3x3v_p1_surfx6_eval_quad_node_21_l(fc); 
    fUpwindQuad_r[21] = hyb_3x3v_p1_surfx6_eval_quad_node_21_l(fr); 
  } 
  if (140788141449279*alpha[16]-140788141449279*alpha[8]+140788141449279*alpha[7]-140788141449279*alpha[6]-140788141449279*alpha[3]+140788141449279*alpha[2]-140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[22] = hyb_3x3v_p1_surfx6_eval_quad_node_22_r(fl); 
    fUpwindQuad_r[22] = hyb_3x3v_p1_surfx6_eval_quad_node_22_r(fc); 
  } else { 
    fUpwindQuad_l[22] = hyb_3x3v_p1_surfx6_eval_quad_node_22_l(fc); 
    fUpwindQuad_r[22] = hyb_3x3v_p1_surfx6_eval_quad_node_22_l(fr); 
  } 
  if (188887112823866*alpha[27]-188887112823866*alpha[22]+188887112823866*alpha[21]-188887112823866*alpha[20]+140788141449279*alpha[16]-188887112823866*alpha[14]+188887112823866*alpha[13]-188887112823866*alpha[12]-140788141449279*alpha[8]+140788141449279*alpha[7]-140788141449279*alpha[6]+188887112823866*alpha[5]-140788141449279*alpha[3]+140788141449279*alpha[2]-140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[23] = hyb_3x3v_p1_surfx6_eval_quad_node_23_r(fl); 
    fUpwindQuad_r[23] = hyb_3x3v_p1_surfx6_eval_quad_node_23_r(fc); 
  } else { 
    fUpwindQuad_l[23] = hyb_3x3v_p1_surfx6_eval_quad_node_23_l(fc); 
    fUpwindQuad_r[23] = hyb_3x3v_p1_surfx6_eval_quad_node_23_l(fr); 
  } 
  if ((-188887112823866*alpha[27])+188887112823866*alpha[26]+188887112823866*alpha[22]-188887112823866*alpha[21]+188887112823866*alpha[20]-188887112823866*alpha[19]+188887112823866*alpha[18]-188887112823866*alpha[17]+140788141449279*alpha[16]+188887112823866*alpha[14]-188887112823866*alpha[13]+188887112823866*alpha[12]-188887112823866*alpha[11]+188887112823866*alpha[10]-188887112823866*alpha[9]-140788141449279*alpha[8]+140788141449279*alpha[7]-140788141449279*alpha[6]-188887112823866*alpha[5]+188887112823866*alpha[4]-140788141449279*alpha[3]+140788141449279*alpha[2]-140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[24] = hyb_3x3v_p1_surfx6_eval_quad_node_24_r(fl); 
    fUpwindQuad_r[24] = hyb_3x3v_p1_surfx6_eval_quad_node_24_r(fc); 
  } else { 
    fUpwindQuad_l[24] = hyb_3x3v_p1_surfx6_eval_quad_node_24_l(fc); 
    fUpwindQuad_r[24] = hyb_3x3v_p1_surfx6_eval_quad_node_24_l(fr); 
  } 
  if (188887112823866*alpha[26]-188887112823866*alpha[19]+188887112823866*alpha[18]-188887112823866*alpha[17]+140788141449279*alpha[16]-188887112823866*alpha[11]+188887112823866*alpha[10]-188887112823866*alpha[9]-140788141449279*alpha[8]+140788141449279*alpha[7]-140788141449279*alpha[6]+188887112823866*alpha[4]-140788141449279*alpha[3]+140788141449279*alpha[2]-140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[25] = hyb_3x3v_p1_surfx6_eval_quad_node_25_r(fl); 
    fUpwindQuad_r[25] = hyb_3x3v_p1_surfx6_eval_quad_node_25_r(fc); 
  } else { 
    fUpwindQuad_l[25] = hyb_3x3v_p1_surfx6_eval_quad_node_25_l(fc); 
    fUpwindQuad_r[25] = hyb_3x3v_p1_surfx6_eval_quad_node_25_l(fr); 
  } 
  if (188887112823866*alpha[27]+188887112823866*alpha[26]-188887112823866*alpha[22]+188887112823866*alpha[21]-188887112823866*alpha[20]-188887112823866*alpha[19]+188887112823866*alpha[18]-188887112823866*alpha[17]+140788141449279*alpha[16]-188887112823866*alpha[14]+188887112823866*alpha[13]-188887112823866*alpha[12]-188887112823866*alpha[11]+188887112823866*alpha[10]-188887112823866*alpha[9]-140788141449279*alpha[8]+140788141449279*alpha[7]-140788141449279*alpha[6]+188887112823866*alpha[5]+188887112823866*alpha[4]-140788141449279*alpha[3]+140788141449279*alpha[2]-140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[26] = hyb_3x3v_p1_surfx6_eval_quad_node_26_r(fl); 
    fUpwindQuad_r[26] = hyb_3x3v_p1_surfx6_eval_quad_node_26_r(fc); 
  } else { 
    fUpwindQuad_l[26] = hyb_3x3v_p1_surfx6_eval_quad_node_26_l(fc); 
    fUpwindQuad_r[26] = hyb_3x3v_p1_surfx6_eval_quad_node_26_l(fr); 
  } 
  if (188887112823866*alpha[27]+188887112823866*alpha[26]-188887112823866*alpha[22]+188887112823866*alpha[21]+188887112823866*alpha[20]-188887112823866*alpha[19]+188887112823866*alpha[18]+188887112823866*alpha[17]-140788141449279*alpha[16]-188887112823866*alpha[14]-188887112823866*alpha[13]+188887112823866*alpha[12]-188887112823866*alpha[11]-188887112823866*alpha[10]+188887112823866*alpha[9]+140788141449279*alpha[8]-140788141449279*alpha[7]-140788141449279*alpha[6]-188887112823866*alpha[5]-188887112823866*alpha[4]+140788141449279*alpha[3]+140788141449279*alpha[2]-140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[27] = hyb_3x3v_p1_surfx6_eval_quad_node_27_r(fl); 
    fUpwindQuad_r[27] = hyb_3x3v_p1_surfx6_eval_quad_node_27_r(fc); 
  } else { 
    fUpwindQuad_l[27] = hyb_3x3v_p1_surfx6_eval_quad_node_27_l(fc); 
    fUpwindQuad_r[27] = hyb_3x3v_p1_surfx6_eval_quad_node_27_l(fr); 
  } 
  if (188887112823866*alpha[26]-188887112823866*alpha[19]+188887112823866*alpha[18]+188887112823866*alpha[17]-140788141449279*alpha[16]-188887112823866*alpha[11]-188887112823866*alpha[10]+188887112823866*alpha[9]+140788141449279*alpha[8]-140788141449279*alpha[7]-140788141449279*alpha[6]-188887112823866*alpha[4]+140788141449279*alpha[3]+140788141449279*alpha[2]-140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[28] = hyb_3x3v_p1_surfx6_eval_quad_node_28_r(fl); 
    fUpwindQuad_r[28] = hyb_3x3v_p1_surfx6_eval_quad_node_28_r(fc); 
  } else { 
    fUpwindQuad_l[28] = hyb_3x3v_p1_surfx6_eval_quad_node_28_l(fc); 
    fUpwindQuad_r[28] = hyb_3x3v_p1_surfx6_eval_quad_node_28_l(fr); 
  } 
  if ((-188887112823866*alpha[27])+188887112823866*alpha[26]+188887112823866*alpha[22]-188887112823866*alpha[21]-188887112823866*alpha[20]-188887112823866*alpha[19]+188887112823866*alpha[18]+188887112823866*alpha[17]-140788141449279*alpha[16]+188887112823866*alpha[14]+188887112823866*alpha[13]-188887112823866*alpha[12]-188887112823866*alpha[11]-188887112823866*alpha[10]+188887112823866*alpha[9]+140788141449279*alpha[8]-140788141449279*alpha[7]-140788141449279*alpha[6]+188887112823866*alpha[5]-188887112823866*alpha[4]+140788141449279*alpha[3]+140788141449279*alpha[2]-140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[29] = hyb_3x3v_p1_surfx6_eval_quad_node_29_r(fl); 
    fUpwindQuad_r[29] = hyb_3x3v_p1_surfx6_eval_quad_node_29_r(fc); 
  } else { 
    fUpwindQuad_l[29] = hyb_3x3v_p1_surfx6_eval_quad_node_29_l(fc); 
    fUpwindQuad_r[29] = hyb_3x3v_p1_surfx6_eval_quad_node_29_l(fr); 
  } 
  if (188887112823866*alpha[27]-188887112823866*alpha[22]+188887112823866*alpha[21]+188887112823866*alpha[20]-140788141449279*alpha[16]-188887112823866*alpha[14]-188887112823866*alpha[13]+188887112823866*alpha[12]+140788141449279*alpha[8]-140788141449279*alpha[7]-140788141449279*alpha[6]-188887112823866*alpha[5]+140788141449279*alpha[3]+140788141449279*alpha[2]-140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[30] = hyb_3x3v_p1_surfx6_eval_quad_node_30_r(fl); 
    fUpwindQuad_r[30] = hyb_3x3v_p1_surfx6_eval_quad_node_30_r(fc); 
  } else { 
    fUpwindQuad_l[30] = hyb_3x3v_p1_surfx6_eval_quad_node_30_l(fc); 
    fUpwindQuad_r[30] = hyb_3x3v_p1_surfx6_eval_quad_node_30_l(fr); 
  } 
  if ((-140788141449279*alpha[16])+140788141449279*alpha[8]-140788141449279*alpha[7]-140788141449279*alpha[6]+140788141449279*alpha[3]+140788141449279*alpha[2]-140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[31] = hyb_3x3v_p1_surfx6_eval_quad_node_31_r(fl); 
    fUpwindQuad_r[31] = hyb_3x3v_p1_surfx6_eval_quad_node_31_r(fc); 
  } else { 
    fUpwindQuad_l[31] = hyb_3x3v_p1_surfx6_eval_quad_node_31_l(fc); 
    fUpwindQuad_r[31] = hyb_3x3v_p1_surfx6_eval_quad_node_31_l(fr); 
  } 
  if ((-188887112823866*alpha[27])+188887112823866*alpha[22]-188887112823866*alpha[21]-188887112823866*alpha[20]-140788141449279*alpha[16]+188887112823866*alpha[14]+188887112823866*alpha[13]-188887112823866*alpha[12]+140788141449279*alpha[8]-140788141449279*alpha[7]-140788141449279*alpha[6]+188887112823866*alpha[5]+140788141449279*alpha[3]+140788141449279*alpha[2]-140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[32] = hyb_3x3v_p1_surfx6_eval_quad_node_32_r(fl); 
    fUpwindQuad_r[32] = hyb_3x3v_p1_surfx6_eval_quad_node_32_r(fc); 
  } else { 
    fUpwindQuad_l[32] = hyb_3x3v_p1_surfx6_eval_quad_node_32_l(fc); 
    fUpwindQuad_r[32] = hyb_3x3v_p1_surfx6_eval_quad_node_32_l(fr); 
  } 
  if (188887112823866*alpha[27]-188887112823866*alpha[26]-188887112823866*alpha[22]+188887112823866*alpha[21]+188887112823866*alpha[20]+188887112823866*alpha[19]-188887112823866*alpha[18]-188887112823866*alpha[17]-140788141449279*alpha[16]-188887112823866*alpha[14]-188887112823866*alpha[13]+188887112823866*alpha[12]+188887112823866*alpha[11]+188887112823866*alpha[10]-188887112823866*alpha[9]+140788141449279*alpha[8]-140788141449279*alpha[7]-140788141449279*alpha[6]-188887112823866*alpha[5]+188887112823866*alpha[4]+140788141449279*alpha[3]+140788141449279*alpha[2]-140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[33] = hyb_3x3v_p1_surfx6_eval_quad_node_33_r(fl); 
    fUpwindQuad_r[33] = hyb_3x3v_p1_surfx6_eval_quad_node_33_r(fc); 
  } else { 
    fUpwindQuad_l[33] = hyb_3x3v_p1_surfx6_eval_quad_node_33_l(fc); 
    fUpwindQuad_r[33] = hyb_3x3v_p1_surfx6_eval_quad_node_33_l(fr); 
  } 
  if ((-188887112823866*alpha[26])+188887112823866*alpha[19]-188887112823866*alpha[18]-188887112823866*alpha[17]-140788141449279*alpha[16]+188887112823866*alpha[11]+188887112823866*alpha[10]-188887112823866*alpha[9]+140788141449279*alpha[8]-140788141449279*alpha[7]-140788141449279*alpha[6]+188887112823866*alpha[4]+140788141449279*alpha[3]+140788141449279*alpha[2]-140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[34] = hyb_3x3v_p1_surfx6_eval_quad_node_34_r(fl); 
    fUpwindQuad_r[34] = hyb_3x3v_p1_surfx6_eval_quad_node_34_r(fc); 
  } else { 
    fUpwindQuad_l[34] = hyb_3x3v_p1_surfx6_eval_quad_node_34_l(fc); 
    fUpwindQuad_r[34] = hyb_3x3v_p1_surfx6_eval_quad_node_34_l(fr); 
  } 
  if ((-188887112823866*alpha[27])-188887112823866*alpha[26]+188887112823866*alpha[22]-188887112823866*alpha[21]-188887112823866*alpha[20]+188887112823866*alpha[19]-188887112823866*alpha[18]-188887112823866*alpha[17]-140788141449279*alpha[16]+188887112823866*alpha[14]+188887112823866*alpha[13]-188887112823866*alpha[12]+188887112823866*alpha[11]+188887112823866*alpha[10]-188887112823866*alpha[9]+140788141449279*alpha[8]-140788141449279*alpha[7]-140788141449279*alpha[6]+188887112823866*alpha[5]+188887112823866*alpha[4]+140788141449279*alpha[3]+140788141449279*alpha[2]-140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[35] = hyb_3x3v_p1_surfx6_eval_quad_node_35_r(fl); 
    fUpwindQuad_r[35] = hyb_3x3v_p1_surfx6_eval_quad_node_35_r(fc); 
  } else { 
    fUpwindQuad_l[35] = hyb_3x3v_p1_surfx6_eval_quad_node_35_l(fc); 
    fUpwindQuad_r[35] = hyb_3x3v_p1_surfx6_eval_quad_node_35_l(fr); 
  } 
  if ((-188887112823866*alpha[27])-188887112823866*alpha[26]-188887112823866*alpha[22]+188887112823866*alpha[21]+188887112823866*alpha[20]-188887112823866*alpha[19]+188887112823866*alpha[18]+188887112823866*alpha[17]+140788141449279*alpha[16]+188887112823866*alpha[14]+188887112823866*alpha[13]-188887112823866*alpha[12]+188887112823866*alpha[11]+188887112823866*alpha[10]-188887112823866*alpha[9]+140788141449279*alpha[8]-140788141449279*alpha[7]-140788141449279*alpha[6]-188887112823866*alpha[5]-188887112823866*alpha[4]-140788141449279*alpha[3]-140788141449279*alpha[2]+140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[36] = hyb_3x3v_p1_surfx6_eval_quad_node_36_r(fl); 
    fUpwindQuad_r[36] = hyb_3x3v_p1_surfx6_eval_quad_node_36_r(fc); 
  } else { 
    fUpwindQuad_l[36] = hyb_3x3v_p1_surfx6_eval_quad_node_36_l(fc); 
    fUpwindQuad_r[36] = hyb_3x3v_p1_surfx6_eval_quad_node_36_l(fr); 
  } 
  if ((-188887112823866*alpha[26])-188887112823866*alpha[19]+188887112823866*alpha[18]+188887112823866*alpha[17]+140788141449279*alpha[16]+188887112823866*alpha[11]+188887112823866*alpha[10]-188887112823866*alpha[9]+140788141449279*alpha[8]-140788141449279*alpha[7]-140788141449279*alpha[6]-188887112823866*alpha[4]-140788141449279*alpha[3]-140788141449279*alpha[2]+140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[37] = hyb_3x3v_p1_surfx6_eval_quad_node_37_r(fl); 
    fUpwindQuad_r[37] = hyb_3x3v_p1_surfx6_eval_quad_node_37_r(fc); 
  } else { 
    fUpwindQuad_l[37] = hyb_3x3v_p1_surfx6_eval_quad_node_37_l(fc); 
    fUpwindQuad_r[37] = hyb_3x3v_p1_surfx6_eval_quad_node_37_l(fr); 
  } 
  if (188887112823866*alpha[27]-188887112823866*alpha[26]+188887112823866*alpha[22]-188887112823866*alpha[21]-188887112823866*alpha[20]-188887112823866*alpha[19]+188887112823866*alpha[18]+188887112823866*alpha[17]+140788141449279*alpha[16]-188887112823866*alpha[14]-188887112823866*alpha[13]+188887112823866*alpha[12]+188887112823866*alpha[11]+188887112823866*alpha[10]-188887112823866*alpha[9]+140788141449279*alpha[8]-140788141449279*alpha[7]-140788141449279*alpha[6]+188887112823866*alpha[5]-188887112823866*alpha[4]-140788141449279*alpha[3]-140788141449279*alpha[2]+140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[38] = hyb_3x3v_p1_surfx6_eval_quad_node_38_r(fl); 
    fUpwindQuad_r[38] = hyb_3x3v_p1_surfx6_eval_quad_node_38_r(fc); 
  } else { 
    fUpwindQuad_l[38] = hyb_3x3v_p1_surfx6_eval_quad_node_38_l(fc); 
    fUpwindQuad_r[38] = hyb_3x3v_p1_surfx6_eval_quad_node_38_l(fr); 
  } 
  if ((-188887112823866*alpha[27])-188887112823866*alpha[22]+188887112823866*alpha[21]+188887112823866*alpha[20]+140788141449279*alpha[16]+188887112823866*alpha[14]+188887112823866*alpha[13]-188887112823866*alpha[12]+140788141449279*alpha[8]-140788141449279*alpha[7]-140788141449279*alpha[6]-188887112823866*alpha[5]-140788141449279*alpha[3]-140788141449279*alpha[2]+140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[39] = hyb_3x3v_p1_surfx6_eval_quad_node_39_r(fl); 
    fUpwindQuad_r[39] = hyb_3x3v_p1_surfx6_eval_quad_node_39_r(fc); 
  } else { 
    fUpwindQuad_l[39] = hyb_3x3v_p1_surfx6_eval_quad_node_39_l(fc); 
    fUpwindQuad_r[39] = hyb_3x3v_p1_surfx6_eval_quad_node_39_l(fr); 
  } 
  if (140788141449279*alpha[16]+140788141449279*alpha[8]-140788141449279*alpha[7]-140788141449279*alpha[6]-140788141449279*alpha[3]-140788141449279*alpha[2]+140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[40] = hyb_3x3v_p1_surfx6_eval_quad_node_40_r(fl); 
    fUpwindQuad_r[40] = hyb_3x3v_p1_surfx6_eval_quad_node_40_r(fc); 
  } else { 
    fUpwindQuad_l[40] = hyb_3x3v_p1_surfx6_eval_quad_node_40_l(fc); 
    fUpwindQuad_r[40] = hyb_3x3v_p1_surfx6_eval_quad_node_40_l(fr); 
  } 
  if (188887112823866*alpha[27]+188887112823866*alpha[22]-188887112823866*alpha[21]-188887112823866*alpha[20]+140788141449279*alpha[16]-188887112823866*alpha[14]-188887112823866*alpha[13]+188887112823866*alpha[12]+140788141449279*alpha[8]-140788141449279*alpha[7]-140788141449279*alpha[6]+188887112823866*alpha[5]-140788141449279*alpha[3]-140788141449279*alpha[2]+140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[41] = hyb_3x3v_p1_surfx6_eval_quad_node_41_r(fl); 
    fUpwindQuad_r[41] = hyb_3x3v_p1_surfx6_eval_quad_node_41_r(fc); 
  } else { 
    fUpwindQuad_l[41] = hyb_3x3v_p1_surfx6_eval_quad_node_41_l(fc); 
    fUpwindQuad_r[41] = hyb_3x3v_p1_surfx6_eval_quad_node_41_l(fr); 
  } 
  if ((-188887112823866*alpha[27])+188887112823866*alpha[26]-188887112823866*alpha[22]+188887112823866*alpha[21]+188887112823866*alpha[20]+188887112823866*alpha[19]-188887112823866*alpha[18]-188887112823866*alpha[17]+140788141449279*alpha[16]+188887112823866*alpha[14]+188887112823866*alpha[13]-188887112823866*alpha[12]-188887112823866*alpha[11]-188887112823866*alpha[10]+188887112823866*alpha[9]+140788141449279*alpha[8]-140788141449279*alpha[7]-140788141449279*alpha[6]-188887112823866*alpha[5]+188887112823866*alpha[4]-140788141449279*alpha[3]-140788141449279*alpha[2]+140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[42] = hyb_3x3v_p1_surfx6_eval_quad_node_42_r(fl); 
    fUpwindQuad_r[42] = hyb_3x3v_p1_surfx6_eval_quad_node_42_r(fc); 
  } else { 
    fUpwindQuad_l[42] = hyb_3x3v_p1_surfx6_eval_quad_node_42_l(fc); 
    fUpwindQuad_r[42] = hyb_3x3v_p1_surfx6_eval_quad_node_42_l(fr); 
  } 
  if (188887112823866*alpha[26]+188887112823866*alpha[19]-188887112823866*alpha[18]-188887112823866*alpha[17]+140788141449279*alpha[16]-188887112823866*alpha[11]-188887112823866*alpha[10]+188887112823866*alpha[9]+140788141449279*alpha[8]-140788141449279*alpha[7]-140788141449279*alpha[6]+188887112823866*alpha[4]-140788141449279*alpha[3]-140788141449279*alpha[2]+140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[43] = hyb_3x3v_p1_surfx6_eval_quad_node_43_r(fl); 
    fUpwindQuad_r[43] = hyb_3x3v_p1_surfx6_eval_quad_node_43_r(fc); 
  } else { 
    fUpwindQuad_l[43] = hyb_3x3v_p1_surfx6_eval_quad_node_43_l(fc); 
    fUpwindQuad_r[43] = hyb_3x3v_p1_surfx6_eval_quad_node_43_l(fr); 
  } 
  if (188887112823866*alpha[27]+188887112823866*alpha[26]+188887112823866*alpha[22]-188887112823866*alpha[21]-188887112823866*alpha[20]+188887112823866*alpha[19]-188887112823866*alpha[18]-188887112823866*alpha[17]+140788141449279*alpha[16]-188887112823866*alpha[14]-188887112823866*alpha[13]+188887112823866*alpha[12]-188887112823866*alpha[11]-188887112823866*alpha[10]+188887112823866*alpha[9]+140788141449279*alpha[8]-140788141449279*alpha[7]-140788141449279*alpha[6]+188887112823866*alpha[5]+188887112823866*alpha[4]-140788141449279*alpha[3]-140788141449279*alpha[2]+140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[44] = hyb_3x3v_p1_surfx6_eval_quad_node_44_r(fl); 
    fUpwindQuad_r[44] = hyb_3x3v_p1_surfx6_eval_quad_node_44_r(fc); 
  } else { 
    fUpwindQuad_l[44] = hyb_3x3v_p1_surfx6_eval_quad_node_44_l(fc); 
    fUpwindQuad_r[44] = hyb_3x3v_p1_surfx6_eval_quad_node_44_l(fr); 
  } 
  if (188887112823866*alpha[27]+188887112823866*alpha[26]+188887112823866*alpha[22]-188887112823866*alpha[21]+188887112823866*alpha[20]+188887112823866*alpha[19]-188887112823866*alpha[18]+188887112823866*alpha[17]-140788141449279*alpha[16]-188887112823866*alpha[14]+188887112823866*alpha[13]-188887112823866*alpha[12]-188887112823866*alpha[11]+188887112823866*alpha[10]-188887112823866*alpha[9]-140788141449279*alpha[8]+140788141449279*alpha[7]-140788141449279*alpha[6]-188887112823866*alpha[5]-188887112823866*alpha[4]+140788141449279*alpha[3]-140788141449279*alpha[2]+140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[45] = hyb_3x3v_p1_surfx6_eval_quad_node_45_r(fl); 
    fUpwindQuad_r[45] = hyb_3x3v_p1_surfx6_eval_quad_node_45_r(fc); 
  } else { 
    fUpwindQuad_l[45] = hyb_3x3v_p1_surfx6_eval_quad_node_45_l(fc); 
    fUpwindQuad_r[45] = hyb_3x3v_p1_surfx6_eval_quad_node_45_l(fr); 
  } 
  if (188887112823866*alpha[26]+188887112823866*alpha[19]-188887112823866*alpha[18]+188887112823866*alpha[17]-140788141449279*alpha[16]-188887112823866*alpha[11]+188887112823866*alpha[10]-188887112823866*alpha[9]-140788141449279*alpha[8]+140788141449279*alpha[7]-140788141449279*alpha[6]-188887112823866*alpha[4]+140788141449279*alpha[3]-140788141449279*alpha[2]+140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[46] = hyb_3x3v_p1_surfx6_eval_quad_node_46_r(fl); 
    fUpwindQuad_r[46] = hyb_3x3v_p1_surfx6_eval_quad_node_46_r(fc); 
  } else { 
    fUpwindQuad_l[46] = hyb_3x3v_p1_surfx6_eval_quad_node_46_l(fc); 
    fUpwindQuad_r[46] = hyb_3x3v_p1_surfx6_eval_quad_node_46_l(fr); 
  } 
  if ((-188887112823866*alpha[27])+188887112823866*alpha[26]-188887112823866*alpha[22]+188887112823866*alpha[21]-188887112823866*alpha[20]+188887112823866*alpha[19]-188887112823866*alpha[18]+188887112823866*alpha[17]-140788141449279*alpha[16]+188887112823866*alpha[14]-188887112823866*alpha[13]+188887112823866*alpha[12]-188887112823866*alpha[11]+188887112823866*alpha[10]-188887112823866*alpha[9]-140788141449279*alpha[8]+140788141449279*alpha[7]-140788141449279*alpha[6]+188887112823866*alpha[5]-188887112823866*alpha[4]+140788141449279*alpha[3]-140788141449279*alpha[2]+140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[47] = hyb_3x3v_p1_surfx6_eval_quad_node_47_r(fl); 
    fUpwindQuad_r[47] = hyb_3x3v_p1_surfx6_eval_quad_node_47_r(fc); 
  } else { 
    fUpwindQuad_l[47] = hyb_3x3v_p1_surfx6_eval_quad_node_47_l(fc); 
    fUpwindQuad_r[47] = hyb_3x3v_p1_surfx6_eval_quad_node_47_l(fr); 
  } 
  if (188887112823866*alpha[27]+188887112823866*alpha[22]-188887112823866*alpha[21]+188887112823866*alpha[20]-140788141449279*alpha[16]-188887112823866*alpha[14]+188887112823866*alpha[13]-188887112823866*alpha[12]-140788141449279*alpha[8]+140788141449279*alpha[7]-140788141449279*alpha[6]-188887112823866*alpha[5]+140788141449279*alpha[3]-140788141449279*alpha[2]+140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[48] = hyb_3x3v_p1_surfx6_eval_quad_node_48_r(fl); 
    fUpwindQuad_r[48] = hyb_3x3v_p1_surfx6_eval_quad_node_48_r(fc); 
  } else { 
    fUpwindQuad_l[48] = hyb_3x3v_p1_surfx6_eval_quad_node_48_l(fc); 
    fUpwindQuad_r[48] = hyb_3x3v_p1_surfx6_eval_quad_node_48_l(fr); 
  } 
  if ((-140788141449279*alpha[16])-140788141449279*alpha[8]+140788141449279*alpha[7]-140788141449279*alpha[6]+140788141449279*alpha[3]-140788141449279*alpha[2]+140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[49] = hyb_3x3v_p1_surfx6_eval_quad_node_49_r(fl); 
    fUpwindQuad_r[49] = hyb_3x3v_p1_surfx6_eval_quad_node_49_r(fc); 
  } else { 
    fUpwindQuad_l[49] = hyb_3x3v_p1_surfx6_eval_quad_node_49_l(fc); 
    fUpwindQuad_r[49] = hyb_3x3v_p1_surfx6_eval_quad_node_49_l(fr); 
  } 
  if ((-188887112823866*alpha[27])-188887112823866*alpha[22]+188887112823866*alpha[21]-188887112823866*alpha[20]-140788141449279*alpha[16]+188887112823866*alpha[14]-188887112823866*alpha[13]+188887112823866*alpha[12]-140788141449279*alpha[8]+140788141449279*alpha[7]-140788141449279*alpha[6]+188887112823866*alpha[5]+140788141449279*alpha[3]-140788141449279*alpha[2]+140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[50] = hyb_3x3v_p1_surfx6_eval_quad_node_50_r(fl); 
    fUpwindQuad_r[50] = hyb_3x3v_p1_surfx6_eval_quad_node_50_r(fc); 
  } else { 
    fUpwindQuad_l[50] = hyb_3x3v_p1_surfx6_eval_quad_node_50_l(fc); 
    fUpwindQuad_r[50] = hyb_3x3v_p1_surfx6_eval_quad_node_50_l(fr); 
  } 
  if (188887112823866*alpha[27]-188887112823866*alpha[26]+188887112823866*alpha[22]-188887112823866*alpha[21]+188887112823866*alpha[20]-188887112823866*alpha[19]+188887112823866*alpha[18]-188887112823866*alpha[17]-140788141449279*alpha[16]-188887112823866*alpha[14]+188887112823866*alpha[13]-188887112823866*alpha[12]+188887112823866*alpha[11]-188887112823866*alpha[10]+188887112823866*alpha[9]-140788141449279*alpha[8]+140788141449279*alpha[7]-140788141449279*alpha[6]-188887112823866*alpha[5]+188887112823866*alpha[4]+140788141449279*alpha[3]-140788141449279*alpha[2]+140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[51] = hyb_3x3v_p1_surfx6_eval_quad_node_51_r(fl); 
    fUpwindQuad_r[51] = hyb_3x3v_p1_surfx6_eval_quad_node_51_r(fc); 
  } else { 
    fUpwindQuad_l[51] = hyb_3x3v_p1_surfx6_eval_quad_node_51_l(fc); 
    fUpwindQuad_r[51] = hyb_3x3v_p1_surfx6_eval_quad_node_51_l(fr); 
  } 
  if ((-188887112823866*alpha[26])-188887112823866*alpha[19]+188887112823866*alpha[18]-188887112823866*alpha[17]-140788141449279*alpha[16]+188887112823866*alpha[11]-188887112823866*alpha[10]+188887112823866*alpha[9]-140788141449279*alpha[8]+140788141449279*alpha[7]-140788141449279*alpha[6]+188887112823866*alpha[4]+140788141449279*alpha[3]-140788141449279*alpha[2]+140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[52] = hyb_3x3v_p1_surfx6_eval_quad_node_52_r(fl); 
    fUpwindQuad_r[52] = hyb_3x3v_p1_surfx6_eval_quad_node_52_r(fc); 
  } else { 
    fUpwindQuad_l[52] = hyb_3x3v_p1_surfx6_eval_quad_node_52_l(fc); 
    fUpwindQuad_r[52] = hyb_3x3v_p1_surfx6_eval_quad_node_52_l(fr); 
  } 
  if ((-188887112823866*alpha[27])-188887112823866*alpha[26]-188887112823866*alpha[22]+188887112823866*alpha[21]-188887112823866*alpha[20]-188887112823866*alpha[19]+188887112823866*alpha[18]-188887112823866*alpha[17]-140788141449279*alpha[16]+188887112823866*alpha[14]-188887112823866*alpha[13]+188887112823866*alpha[12]+188887112823866*alpha[11]-188887112823866*alpha[10]+188887112823866*alpha[9]-140788141449279*alpha[8]+140788141449279*alpha[7]-140788141449279*alpha[6]+188887112823866*alpha[5]+188887112823866*alpha[4]+140788141449279*alpha[3]-140788141449279*alpha[2]+140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[53] = hyb_3x3v_p1_surfx6_eval_quad_node_53_r(fl); 
    fUpwindQuad_r[53] = hyb_3x3v_p1_surfx6_eval_quad_node_53_r(fc); 
  } else { 
    fUpwindQuad_l[53] = hyb_3x3v_p1_surfx6_eval_quad_node_53_l(fc); 
    fUpwindQuad_r[53] = hyb_3x3v_p1_surfx6_eval_quad_node_53_l(fr); 
  } 
  if (188887112823866*alpha[27]+188887112823866*alpha[26]+188887112823866*alpha[22]+188887112823866*alpha[21]-188887112823866*alpha[20]+188887112823866*alpha[19]+188887112823866*alpha[18]-188887112823866*alpha[17]-140788141449279*alpha[16]+188887112823866*alpha[14]-188887112823866*alpha[13]-188887112823866*alpha[12]+188887112823866*alpha[11]-188887112823866*alpha[10]-188887112823866*alpha[9]-140788141449279*alpha[8]-140788141449279*alpha[7]+140788141449279*alpha[6]-188887112823866*alpha[5]-188887112823866*alpha[4]-140788141449279*alpha[3]+140788141449279*alpha[2]+140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[54] = hyb_3x3v_p1_surfx6_eval_quad_node_54_r(fl); 
    fUpwindQuad_r[54] = hyb_3x3v_p1_surfx6_eval_quad_node_54_r(fc); 
  } else { 
    fUpwindQuad_l[54] = hyb_3x3v_p1_surfx6_eval_quad_node_54_l(fc); 
    fUpwindQuad_r[54] = hyb_3x3v_p1_surfx6_eval_quad_node_54_l(fr); 
  } 
  if (188887112823866*alpha[26]+188887112823866*alpha[19]+188887112823866*alpha[18]-188887112823866*alpha[17]-140788141449279*alpha[16]+188887112823866*alpha[11]-188887112823866*alpha[10]-188887112823866*alpha[9]-140788141449279*alpha[8]-140788141449279*alpha[7]+140788141449279*alpha[6]-188887112823866*alpha[4]-140788141449279*alpha[3]+140788141449279*alpha[2]+140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[55] = hyb_3x3v_p1_surfx6_eval_quad_node_55_r(fl); 
    fUpwindQuad_r[55] = hyb_3x3v_p1_surfx6_eval_quad_node_55_r(fc); 
  } else { 
    fUpwindQuad_l[55] = hyb_3x3v_p1_surfx6_eval_quad_node_55_l(fc); 
    fUpwindQuad_r[55] = hyb_3x3v_p1_surfx6_eval_quad_node_55_l(fr); 
  } 
  if ((-188887112823866*alpha[27])+188887112823866*alpha[26]-188887112823866*alpha[22]-188887112823866*alpha[21]+188887112823866*alpha[20]+188887112823866*alpha[19]+188887112823866*alpha[18]-188887112823866*alpha[17]-140788141449279*alpha[16]-188887112823866*alpha[14]+188887112823866*alpha[13]+188887112823866*alpha[12]+188887112823866*alpha[11]-188887112823866*alpha[10]-188887112823866*alpha[9]-140788141449279*alpha[8]-140788141449279*alpha[7]+140788141449279*alpha[6]+188887112823866*alpha[5]-188887112823866*alpha[4]-140788141449279*alpha[3]+140788141449279*alpha[2]+140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[56] = hyb_3x3v_p1_surfx6_eval_quad_node_56_r(fl); 
    fUpwindQuad_r[56] = hyb_3x3v_p1_surfx6_eval_quad_node_56_r(fc); 
  } else { 
    fUpwindQuad_l[56] = hyb_3x3v_p1_surfx6_eval_quad_node_56_l(fc); 
    fUpwindQuad_r[56] = hyb_3x3v_p1_surfx6_eval_quad_node_56_l(fr); 
  } 
  if (188887112823866*alpha[27]+188887112823866*alpha[22]+188887112823866*alpha[21]-188887112823866*alpha[20]-140788141449279*alpha[16]+188887112823866*alpha[14]-188887112823866*alpha[13]-188887112823866*alpha[12]-140788141449279*alpha[8]-140788141449279*alpha[7]+140788141449279*alpha[6]-188887112823866*alpha[5]-140788141449279*alpha[3]+140788141449279*alpha[2]+140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[57] = hyb_3x3v_p1_surfx6_eval_quad_node_57_r(fl); 
    fUpwindQuad_r[57] = hyb_3x3v_p1_surfx6_eval_quad_node_57_r(fc); 
  } else { 
    fUpwindQuad_l[57] = hyb_3x3v_p1_surfx6_eval_quad_node_57_l(fc); 
    fUpwindQuad_r[57] = hyb_3x3v_p1_surfx6_eval_quad_node_57_l(fr); 
  } 
  if ((-140788141449279*alpha[16])-140788141449279*alpha[8]-140788141449279*alpha[7]+140788141449279*alpha[6]-140788141449279*alpha[3]+140788141449279*alpha[2]+140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[58] = hyb_3x3v_p1_surfx6_eval_quad_node_58_r(fl); 
    fUpwindQuad_r[58] = hyb_3x3v_p1_surfx6_eval_quad_node_58_r(fc); 
  } else { 
    fUpwindQuad_l[58] = hyb_3x3v_p1_surfx6_eval_quad_node_58_l(fc); 
    fUpwindQuad_r[58] = hyb_3x3v_p1_surfx6_eval_quad_node_58_l(fr); 
  } 
  if ((-188887112823866*alpha[27])-188887112823866*alpha[22]-188887112823866*alpha[21]+188887112823866*alpha[20]-140788141449279*alpha[16]-188887112823866*alpha[14]+188887112823866*alpha[13]+188887112823866*alpha[12]-140788141449279*alpha[8]-140788141449279*alpha[7]+140788141449279*alpha[6]+188887112823866*alpha[5]-140788141449279*alpha[3]+140788141449279*alpha[2]+140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[59] = hyb_3x3v_p1_surfx6_eval_quad_node_59_r(fl); 
    fUpwindQuad_r[59] = hyb_3x3v_p1_surfx6_eval_quad_node_59_r(fc); 
  } else { 
    fUpwindQuad_l[59] = hyb_3x3v_p1_surfx6_eval_quad_node_59_l(fc); 
    fUpwindQuad_r[59] = hyb_3x3v_p1_surfx6_eval_quad_node_59_l(fr); 
  } 
  if (188887112823866*alpha[27]-188887112823866*alpha[26]+188887112823866*alpha[22]+188887112823866*alpha[21]-188887112823866*alpha[20]-188887112823866*alpha[19]-188887112823866*alpha[18]+188887112823866*alpha[17]-140788141449279*alpha[16]+188887112823866*alpha[14]-188887112823866*alpha[13]-188887112823866*alpha[12]-188887112823866*alpha[11]+188887112823866*alpha[10]+188887112823866*alpha[9]-140788141449279*alpha[8]-140788141449279*alpha[7]+140788141449279*alpha[6]-188887112823866*alpha[5]+188887112823866*alpha[4]-140788141449279*alpha[3]+140788141449279*alpha[2]+140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[60] = hyb_3x3v_p1_surfx6_eval_quad_node_60_r(fl); 
    fUpwindQuad_r[60] = hyb_3x3v_p1_surfx6_eval_quad_node_60_r(fc); 
  } else { 
    fUpwindQuad_l[60] = hyb_3x3v_p1_surfx6_eval_quad_node_60_l(fc); 
    fUpwindQuad_r[60] = hyb_3x3v_p1_surfx6_eval_quad_node_60_l(fr); 
  } 
  if ((-188887112823866*alpha[26])-188887112823866*alpha[19]-188887112823866*alpha[18]+188887112823866*alpha[17]-140788141449279*alpha[16]-188887112823866*alpha[11]+188887112823866*alpha[10]+188887112823866*alpha[9]-140788141449279*alpha[8]-140788141449279*alpha[7]+140788141449279*alpha[6]+188887112823866*alpha[4]-140788141449279*alpha[3]+140788141449279*alpha[2]+140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[61] = hyb_3x3v_p1_surfx6_eval_quad_node_61_r(fl); 
    fUpwindQuad_r[61] = hyb_3x3v_p1_surfx6_eval_quad_node_61_r(fc); 
  } else { 
    fUpwindQuad_l[61] = hyb_3x3v_p1_surfx6_eval_quad_node_61_l(fc); 
    fUpwindQuad_r[61] = hyb_3x3v_p1_surfx6_eval_quad_node_61_l(fr); 
  } 
  if ((-188887112823866*alpha[27])-188887112823866*alpha[26]-188887112823866*alpha[22]-188887112823866*alpha[21]+188887112823866*alpha[20]-188887112823866*alpha[19]-188887112823866*alpha[18]+188887112823866*alpha[17]-140788141449279*alpha[16]-188887112823866*alpha[14]+188887112823866*alpha[13]+188887112823866*alpha[12]-188887112823866*alpha[11]+188887112823866*alpha[10]+188887112823866*alpha[9]-140788141449279*alpha[8]-140788141449279*alpha[7]+140788141449279*alpha[6]+188887112823866*alpha[5]+188887112823866*alpha[4]-140788141449279*alpha[3]+140788141449279*alpha[2]+140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[62] = hyb_3x3v_p1_surfx6_eval_quad_node_62_r(fl); 
    fUpwindQuad_r[62] = hyb_3x3v_p1_surfx6_eval_quad_node_62_r(fc); 
  } else { 
    fUpwindQuad_l[62] = hyb_3x3v_p1_surfx6_eval_quad_node_62_l(fc); 
    fUpwindQuad_r[62] = hyb_3x3v_p1_surfx6_eval_quad_node_62_l(fr); 
  } 
  if ((-188887112823866*alpha[27])-188887112823866*alpha[26]-188887112823866*alpha[22]-188887112823866*alpha[21]-188887112823866*alpha[20]-188887112823866*alpha[19]-188887112823866*alpha[18]-188887112823866*alpha[17]+140788141449279*alpha[16]-188887112823866*alpha[14]-188887112823866*alpha[13]-188887112823866*alpha[12]-188887112823866*alpha[11]-188887112823866*alpha[10]-188887112823866*alpha[9]+140788141449279*alpha[8]+140788141449279*alpha[7]+140788141449279*alpha[6]-188887112823866*alpha[5]-188887112823866*alpha[4]+140788141449279*alpha[3]+140788141449279*alpha[2]+140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[63] = hyb_3x3v_p1_surfx6_eval_quad_node_63_r(fl); 
    fUpwindQuad_r[63] = hyb_3x3v_p1_surfx6_eval_quad_node_63_r(fc); 
  } else { 
    fUpwindQuad_l[63] = hyb_3x3v_p1_surfx6_eval_quad_node_63_l(fc); 
    fUpwindQuad_r[63] = hyb_3x3v_p1_surfx6_eval_quad_node_63_l(fr); 
  } 
  if ((-188887112823866*alpha[26])-188887112823866*alpha[19]-188887112823866*alpha[18]-188887112823866*alpha[17]+140788141449279*alpha[16]-188887112823866*alpha[11]-188887112823866*alpha[10]-188887112823866*alpha[9]+140788141449279*alpha[8]+140788141449279*alpha[7]+140788141449279*alpha[6]-188887112823866*alpha[4]+140788141449279*alpha[3]+140788141449279*alpha[2]+140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[64] = hyb_3x3v_p1_surfx6_eval_quad_node_64_r(fl); 
    fUpwindQuad_r[64] = hyb_3x3v_p1_surfx6_eval_quad_node_64_r(fc); 
  } else { 
    fUpwindQuad_l[64] = hyb_3x3v_p1_surfx6_eval_quad_node_64_l(fc); 
    fUpwindQuad_r[64] = hyb_3x3v_p1_surfx6_eval_quad_node_64_l(fr); 
  } 
  if (188887112823866*alpha[27]-188887112823866*alpha[26]+188887112823866*alpha[22]+188887112823866*alpha[21]+188887112823866*alpha[20]-188887112823866*alpha[19]-188887112823866*alpha[18]-188887112823866*alpha[17]+140788141449279*alpha[16]+188887112823866*alpha[14]+188887112823866*alpha[13]+188887112823866*alpha[12]-188887112823866*alpha[11]-188887112823866*alpha[10]-188887112823866*alpha[9]+140788141449279*alpha[8]+140788141449279*alpha[7]+140788141449279*alpha[6]+188887112823866*alpha[5]-188887112823866*alpha[4]+140788141449279*alpha[3]+140788141449279*alpha[2]+140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[65] = hyb_3x3v_p1_surfx6_eval_quad_node_65_r(fl); 
    fUpwindQuad_r[65] = hyb_3x3v_p1_surfx6_eval_quad_node_65_r(fc); 
  } else { 
    fUpwindQuad_l[65] = hyb_3x3v_p1_surfx6_eval_quad_node_65_l(fc); 
    fUpwindQuad_r[65] = hyb_3x3v_p1_surfx6_eval_quad_node_65_l(fr); 
  } 
  if ((-188887112823866*alpha[27])-188887112823866*alpha[22]-188887112823866*alpha[21]-188887112823866*alpha[20]+140788141449279*alpha[16]-188887112823866*alpha[14]-188887112823866*alpha[13]-188887112823866*alpha[12]+140788141449279*alpha[8]+140788141449279*alpha[7]+140788141449279*alpha[6]-188887112823866*alpha[5]+140788141449279*alpha[3]+140788141449279*alpha[2]+140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[66] = hyb_3x3v_p1_surfx6_eval_quad_node_66_r(fl); 
    fUpwindQuad_r[66] = hyb_3x3v_p1_surfx6_eval_quad_node_66_r(fc); 
  } else { 
    fUpwindQuad_l[66] = hyb_3x3v_p1_surfx6_eval_quad_node_66_l(fc); 
    fUpwindQuad_r[66] = hyb_3x3v_p1_surfx6_eval_quad_node_66_l(fr); 
  } 
  if (140788141449279*alpha[16]+140788141449279*alpha[8]+140788141449279*alpha[7]+140788141449279*alpha[6]+140788141449279*alpha[3]+140788141449279*alpha[2]+140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[67] = hyb_3x3v_p1_surfx6_eval_quad_node_67_r(fl); 
    fUpwindQuad_r[67] = hyb_3x3v_p1_surfx6_eval_quad_node_67_r(fc); 
  } else { 
    fUpwindQuad_l[67] = hyb_3x3v_p1_surfx6_eval_quad_node_67_l(fc); 
    fUpwindQuad_r[67] = hyb_3x3v_p1_surfx6_eval_quad_node_67_l(fr); 
  } 
  if (188887112823866*alpha[27]+188887112823866*alpha[22]+188887112823866*alpha[21]+188887112823866*alpha[20]+140788141449279*alpha[16]+188887112823866*alpha[14]+188887112823866*alpha[13]+188887112823866*alpha[12]+140788141449279*alpha[8]+140788141449279*alpha[7]+140788141449279*alpha[6]+188887112823866*alpha[5]+140788141449279*alpha[3]+140788141449279*alpha[2]+140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[68] = hyb_3x3v_p1_surfx6_eval_quad_node_68_r(fl); 
    fUpwindQuad_r[68] = hyb_3x3v_p1_surfx6_eval_quad_node_68_r(fc); 
  } else { 
    fUpwindQuad_l[68] = hyb_3x3v_p1_surfx6_eval_quad_node_68_l(fc); 
    fUpwindQuad_r[68] = hyb_3x3v_p1_surfx6_eval_quad_node_68_l(fr); 
  } 
  if ((-188887112823866*alpha[27])+188887112823866*alpha[26]-188887112823866*alpha[22]-188887112823866*alpha[21]-188887112823866*alpha[20]+188887112823866*alpha[19]+188887112823866*alpha[18]+188887112823866*alpha[17]+140788141449279*alpha[16]-188887112823866*alpha[14]-188887112823866*alpha[13]-188887112823866*alpha[12]+188887112823866*alpha[11]+188887112823866*alpha[10]+188887112823866*alpha[9]+140788141449279*alpha[8]+140788141449279*alpha[7]+140788141449279*alpha[6]-188887112823866*alpha[5]+188887112823866*alpha[4]+140788141449279*alpha[3]+140788141449279*alpha[2]+140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[69] = hyb_3x3v_p1_surfx6_eval_quad_node_69_r(fl); 
    fUpwindQuad_r[69] = hyb_3x3v_p1_surfx6_eval_quad_node_69_r(fc); 
  } else { 
    fUpwindQuad_l[69] = hyb_3x3v_p1_surfx6_eval_quad_node_69_l(fc); 
    fUpwindQuad_r[69] = hyb_3x3v_p1_surfx6_eval_quad_node_69_l(fr); 
  } 
  if (188887112823866*alpha[26]+188887112823866*alpha[19]+188887112823866*alpha[18]+188887112823866*alpha[17]+140788141449279*alpha[16]+188887112823866*alpha[11]+188887112823866*alpha[10]+188887112823866*alpha[9]+140788141449279*alpha[8]+140788141449279*alpha[7]+140788141449279*alpha[6]+188887112823866*alpha[4]+140788141449279*alpha[3]+140788141449279*alpha[2]+140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[70] = hyb_3x3v_p1_surfx6_eval_quad_node_70_r(fl); 
    fUpwindQuad_r[70] = hyb_3x3v_p1_surfx6_eval_quad_node_70_r(fc); 
  } else { 
    fUpwindQuad_l[70] = hyb_3x3v_p1_surfx6_eval_quad_node_70_l(fc); 
    fUpwindQuad_r[70] = hyb_3x3v_p1_surfx6_eval_quad_node_70_l(fr); 
  } 
  if (188887112823866*alpha[27]+188887112823866*alpha[26]+188887112823866*alpha[22]+188887112823866*alpha[21]+188887112823866*alpha[20]+188887112823866*alpha[19]+188887112823866*alpha[18]+188887112823866*alpha[17]+140788141449279*alpha[16]+188887112823866*alpha[14]+188887112823866*alpha[13]+188887112823866*alpha[12]+188887112823866*alpha[11]+188887112823866*alpha[10]+188887112823866*alpha[9]+140788141449279*alpha[8]+140788141449279*alpha[7]+140788141449279*alpha[6]+188887112823866*alpha[5]+188887112823866*alpha[4]+140788141449279*alpha[3]+140788141449279*alpha[2]+140788141449279*alpha[1]+140788141449279*alpha[0] > 0) { 
    fUpwindQuad_l[71] = hyb_3x3v_p1_surfx6_eval_quad_node_71_r(fl); 
    fUpwindQuad_r[71] = hyb_3x3v_p1_surfx6_eval_quad_node_71_r(fc); 
  } else { 
    fUpwindQuad_l[71] = hyb_3x3v_p1_surfx6_eval_quad_node_71_l(fc); 
    fUpwindQuad_r[71] = hyb_3x3v_p1_surfx6_eval_quad_node_71_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_3x3v_p1_vdir_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  hyb_3x3v_p1_vdir_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.1767766952966368*(alpha[27]*fUpwind_l[27]+alpha[26]*fUpwind_l[26]+alpha[22]*fUpwind_l[22]+alpha[21]*fUpwind_l[21]+alpha[20]*fUpwind_l[20]+alpha[19]*fUpwind_l[19]+alpha[18]*fUpwind_l[18]+alpha[17]*fUpwind_l[17]+alpha[16]*fUpwind_l[16]+alpha[14]*fUpwind_l[14]+alpha[13]*fUpwind_l[13]+alpha[12]*fUpwind_l[12]+alpha[11]*fUpwind_l[11]+alpha[10]*fUpwind_l[10]+alpha[9]*fUpwind_l[9]+alpha[8]*fUpwind_l[8]+alpha[7]*fUpwind_l[7]+alpha[6]*fUpwind_l[6]+alpha[5]*fUpwind_l[5]+alpha[4]*fUpwind_l[4]+alpha[3]*fUpwind_l[3]+alpha[2]*fUpwind_l[2]+alpha[1]*fUpwind_l[1]+alpha[0]*fUpwind_l[0]); 
  Ghat_l[1] = 0.1767766952966368*(alpha[22]*fUpwind_l[27]+fUpwind_l[22]*alpha[27]+alpha[19]*fUpwind_l[26]+fUpwind_l[19]*alpha[26]+alpha[14]*fUpwind_l[21]+fUpwind_l[14]*alpha[21]+alpha[13]*fUpwind_l[20]+fUpwind_l[13]*alpha[20]+alpha[11]*fUpwind_l[18]+fUpwind_l[11]*alpha[18]+alpha[10]*fUpwind_l[17]+fUpwind_l[10]*alpha[17]+alpha[8]*fUpwind_l[16]+fUpwind_l[8]*alpha[16]+alpha[5]*fUpwind_l[12]+fUpwind_l[5]*alpha[12]+alpha[4]*fUpwind_l[9]+fUpwind_l[4]*alpha[9]+alpha[3]*fUpwind_l[7]+fUpwind_l[3]*alpha[7]+alpha[2]*fUpwind_l[6]+fUpwind_l[2]*alpha[6]+alpha[0]*fUpwind_l[1]+fUpwind_l[0]*alpha[1]); 
  Ghat_l[2] = 0.1767766952966368*(alpha[21]*fUpwind_l[27]+fUpwind_l[21]*alpha[27]+alpha[18]*fUpwind_l[26]+fUpwind_l[18]*alpha[26]+alpha[14]*fUpwind_l[22]+fUpwind_l[14]*alpha[22]+alpha[12]*fUpwind_l[20]+fUpwind_l[12]*alpha[20]+alpha[11]*fUpwind_l[19]+fUpwind_l[11]*alpha[19]+alpha[9]*fUpwind_l[17]+fUpwind_l[9]*alpha[17]+alpha[7]*fUpwind_l[16]+fUpwind_l[7]*alpha[16]+alpha[5]*fUpwind_l[13]+fUpwind_l[5]*alpha[13]+alpha[4]*fUpwind_l[10]+fUpwind_l[4]*alpha[10]+alpha[3]*fUpwind_l[8]+fUpwind_l[3]*alpha[8]+alpha[1]*fUpwind_l[6]+fUpwind_l[1]*alpha[6]+alpha[0]*fUpwind_l[2]+fUpwind_l[0]*alpha[2]); 
  Ghat_l[3] = 0.1767766952966368*(alpha[20]*fUpwind_l[27]+fUpwind_l[20]*alpha[27]+alpha[17]*fUpwind_l[26]+fUpwind_l[17]*alpha[26]+alpha[13]*fUpwind_l[22]+fUpwind_l[13]*alpha[22]+alpha[12]*fUpwind_l[21]+fUpwind_l[12]*alpha[21]+alpha[10]*fUpwind_l[19]+fUpwind_l[10]*alpha[19]+alpha[9]*fUpwind_l[18]+fUpwind_l[9]*alpha[18]+alpha[6]*fUpwind_l[16]+fUpwind_l[6]*alpha[16]+alpha[5]*fUpwind_l[14]+fUpwind_l[5]*alpha[14]+alpha[4]*fUpwind_l[11]+fUpwind_l[4]*alpha[11]+alpha[2]*fUpwind_l[8]+fUpwind_l[2]*alpha[8]+alpha[1]*fUpwind_l[7]+fUpwind_l[1]*alpha[7]+alpha[0]*fUpwind_l[3]+fUpwind_l[0]*alpha[3]); 
  Ghat_l[4] = 0.01178511301977579*(13.41640786499874*alpha[26]*fUpwind_l[43]+13.41640786499874*(alpha[19]*fUpwind_l[39]+alpha[18]*fUpwind_l[38]+alpha[17]*fUpwind_l[37])+13.41640786499874*(alpha[11]*fUpwind_l[35]+alpha[10]*fUpwind_l[34]+alpha[9]*fUpwind_l[33])+13.41640786499874*alpha[4]*fUpwind_l[32]+15.0*(alpha[27]*fUpwind_l[31]+alpha[22]*fUpwind_l[30]+alpha[21]*fUpwind_l[29]+alpha[20]*fUpwind_l[28]+alpha[16]*fUpwind_l[26]+fUpwind_l[16]*alpha[26]+alpha[14]*fUpwind_l[25]+alpha[13]*fUpwind_l[24]+alpha[12]*fUpwind_l[23]+alpha[8]*fUpwind_l[19]+fUpwind_l[8]*alpha[19]+alpha[7]*fUpwind_l[18]+fUpwind_l[7]*alpha[18]+alpha[6]*fUpwind_l[17]+fUpwind_l[6]*alpha[17]+alpha[5]*fUpwind_l[15]+alpha[3]*fUpwind_l[11]+fUpwind_l[3]*alpha[11]+alpha[2]*fUpwind_l[10]+fUpwind_l[2]*alpha[10]+alpha[1]*fUpwind_l[9]+fUpwind_l[1]*alpha[9]+alpha[0]*fUpwind_l[4]+fUpwind_l[0]*alpha[4])); 
  Ghat_l[5] = 0.01178511301977579*(13.41640786499874*alpha[27]*fUpwind_l[59]+13.41640786499874*(alpha[22]*fUpwind_l[55]+alpha[21]*fUpwind_l[54]+alpha[20]*fUpwind_l[53])+13.41640786499874*(alpha[14]*fUpwind_l[51]+alpha[13]*fUpwind_l[50]+alpha[12]*fUpwind_l[49])+13.41640786499874*alpha[5]*fUpwind_l[48]+15.0*(alpha[26]*fUpwind_l[31]+alpha[19]*fUpwind_l[30]+alpha[18]*fUpwind_l[29]+alpha[17]*fUpwind_l[28]+alpha[16]*fUpwind_l[27]+fUpwind_l[16]*alpha[27]+alpha[11]*fUpwind_l[25]+alpha[10]*fUpwind_l[24]+alpha[9]*fUpwind_l[23]+alpha[8]*fUpwind_l[22]+fUpwind_l[8]*alpha[22]+alpha[7]*fUpwind_l[21]+fUpwind_l[7]*alpha[21]+alpha[6]*fUpwind_l[20]+fUpwind_l[6]*alpha[20]+alpha[4]*fUpwind_l[15]+alpha[3]*fUpwind_l[14]+fUpwind_l[3]*alpha[14]+alpha[2]*fUpwind_l[13]+fUpwind_l[2]*alpha[13]+alpha[1]*fUpwind_l[12]+fUpwind_l[1]*alpha[12]+alpha[0]*fUpwind_l[5]+fUpwind_l[0]*alpha[5])); 
  Ghat_l[6] = 0.1767766952966368*(alpha[14]*fUpwind_l[27]+fUpwind_l[14]*alpha[27]+alpha[11]*fUpwind_l[26]+fUpwind_l[11]*alpha[26]+alpha[21]*fUpwind_l[22]+fUpwind_l[21]*alpha[22]+alpha[5]*fUpwind_l[20]+fUpwind_l[5]*alpha[20]+alpha[18]*fUpwind_l[19]+fUpwind_l[18]*alpha[19]+alpha[4]*fUpwind_l[17]+fUpwind_l[4]*alpha[17]+alpha[3]*fUpwind_l[16]+fUpwind_l[3]*alpha[16]+alpha[12]*fUpwind_l[13]+fUpwind_l[12]*alpha[13]+alpha[9]*fUpwind_l[10]+fUpwind_l[9]*alpha[10]+alpha[7]*fUpwind_l[8]+fUpwind_l[7]*alpha[8]+alpha[0]*fUpwind_l[6]+fUpwind_l[0]*alpha[6]+alpha[1]*fUpwind_l[2]+fUpwind_l[1]*alpha[2]); 
  Ghat_l[7] = 0.1767766952966368*(alpha[13]*fUpwind_l[27]+fUpwind_l[13]*alpha[27]+alpha[10]*fUpwind_l[26]+fUpwind_l[10]*alpha[26]+alpha[20]*fUpwind_l[22]+fUpwind_l[20]*alpha[22]+alpha[5]*fUpwind_l[21]+fUpwind_l[5]*alpha[21]+alpha[17]*fUpwind_l[19]+fUpwind_l[17]*alpha[19]+alpha[4]*fUpwind_l[18]+fUpwind_l[4]*alpha[18]+alpha[2]*fUpwind_l[16]+fUpwind_l[2]*alpha[16]+alpha[12]*fUpwind_l[14]+fUpwind_l[12]*alpha[14]+alpha[9]*fUpwind_l[11]+fUpwind_l[9]*alpha[11]+alpha[6]*fUpwind_l[8]+fUpwind_l[6]*alpha[8]+alpha[0]*fUpwind_l[7]+fUpwind_l[0]*alpha[7]+alpha[1]*fUpwind_l[3]+fUpwind_l[1]*alpha[3]); 
  Ghat_l[8] = 0.1767766952966368*(alpha[12]*fUpwind_l[27]+fUpwind_l[12]*alpha[27]+alpha[9]*fUpwind_l[26]+fUpwind_l[9]*alpha[26]+alpha[5]*fUpwind_l[22]+fUpwind_l[5]*alpha[22]+alpha[20]*fUpwind_l[21]+fUpwind_l[20]*alpha[21]+alpha[4]*fUpwind_l[19]+fUpwind_l[4]*alpha[19]+alpha[17]*fUpwind_l[18]+fUpwind_l[17]*alpha[18]+alpha[1]*fUpwind_l[16]+fUpwind_l[1]*alpha[16]+alpha[13]*fUpwind_l[14]+fUpwind_l[13]*alpha[14]+alpha[10]*fUpwind_l[11]+fUpwind_l[10]*alpha[11]+alpha[0]*fUpwind_l[8]+fUpwind_l[0]*alpha[8]+alpha[6]*fUpwind_l[7]+fUpwind_l[6]*alpha[7]+alpha[2]*fUpwind_l[3]+fUpwind_l[2]*alpha[3]); 
  Ghat_l[9] = 0.01178511301977579*(13.41640786499874*alpha[19]*fUpwind_l[43]+13.41640786499874*(alpha[26]*fUpwind_l[39]+alpha[11]*fUpwind_l[38]+alpha[10]*fUpwind_l[37])+13.41640786499874*(alpha[18]*fUpwind_l[35]+alpha[17]*fUpwind_l[34]+alpha[4]*fUpwind_l[33])+13.41640786499874*alpha[9]*fUpwind_l[32]+15.0*(alpha[22]*fUpwind_l[31]+alpha[27]*fUpwind_l[30]+alpha[14]*fUpwind_l[29]+alpha[13]*fUpwind_l[28]+alpha[8]*fUpwind_l[26]+fUpwind_l[8]*alpha[26]+alpha[21]*fUpwind_l[25]+alpha[20]*fUpwind_l[24]+alpha[5]*fUpwind_l[23]+alpha[16]*fUpwind_l[19]+fUpwind_l[16]*alpha[19]+alpha[3]*fUpwind_l[18]+fUpwind_l[3]*alpha[18]+alpha[2]*fUpwind_l[17]+fUpwind_l[2]*alpha[17]+alpha[12]*fUpwind_l[15]+alpha[7]*fUpwind_l[11]+fUpwind_l[7]*alpha[11]+alpha[6]*fUpwind_l[10]+fUpwind_l[6]*alpha[10]+alpha[0]*fUpwind_l[9]+fUpwind_l[0]*alpha[9]+alpha[1]*fUpwind_l[4]+fUpwind_l[1]*alpha[4])); 
  Ghat_l[10] = 0.01178511301977579*(13.41640786499874*alpha[18]*fUpwind_l[43]+13.41640786499874*(alpha[11]*fUpwind_l[39]+alpha[26]*fUpwind_l[38]+alpha[9]*fUpwind_l[37])+13.41640786499874*(alpha[19]*fUpwind_l[35]+alpha[4]*fUpwind_l[34]+alpha[17]*fUpwind_l[33])+13.41640786499874*alpha[10]*fUpwind_l[32]+15.0*(alpha[21]*fUpwind_l[31]+alpha[14]*fUpwind_l[30]+alpha[27]*fUpwind_l[29]+alpha[12]*fUpwind_l[28]+alpha[7]*fUpwind_l[26]+fUpwind_l[7]*alpha[26]+alpha[22]*fUpwind_l[25]+alpha[5]*fUpwind_l[24]+alpha[20]*fUpwind_l[23]+alpha[3]*fUpwind_l[19]+fUpwind_l[3]*alpha[19]+alpha[16]*fUpwind_l[18]+fUpwind_l[16]*alpha[18]+alpha[1]*fUpwind_l[17]+fUpwind_l[1]*alpha[17]+alpha[13]*fUpwind_l[15]+alpha[8]*fUpwind_l[11]+fUpwind_l[8]*alpha[11]+alpha[0]*fUpwind_l[10]+fUpwind_l[0]*alpha[10]+alpha[6]*fUpwind_l[9]+fUpwind_l[6]*alpha[9]+alpha[2]*fUpwind_l[4]+fUpwind_l[2]*alpha[4])); 
  Ghat_l[11] = 0.01178511301977579*(13.41640786499874*alpha[17]*fUpwind_l[43]+13.41640786499874*(alpha[10]*fUpwind_l[39]+alpha[9]*fUpwind_l[38]+alpha[26]*fUpwind_l[37])+13.41640786499874*(alpha[4]*fUpwind_l[35]+alpha[19]*fUpwind_l[34]+alpha[18]*fUpwind_l[33])+13.41640786499874*alpha[11]*fUpwind_l[32]+15.0*(alpha[20]*fUpwind_l[31]+alpha[13]*fUpwind_l[30]+alpha[12]*fUpwind_l[29]+alpha[27]*fUpwind_l[28]+alpha[6]*fUpwind_l[26]+fUpwind_l[6]*alpha[26]+alpha[5]*fUpwind_l[25]+alpha[22]*fUpwind_l[24]+alpha[21]*fUpwind_l[23]+alpha[2]*fUpwind_l[19]+fUpwind_l[2]*alpha[19]+alpha[1]*fUpwind_l[18]+fUpwind_l[1]*alpha[18]+alpha[16]*fUpwind_l[17]+fUpwind_l[16]*alpha[17]+alpha[14]*fUpwind_l[15]+alpha[0]*fUpwind_l[11]+fUpwind_l[0]*alpha[11]+alpha[8]*fUpwind_l[10]+fUpwind_l[8]*alpha[10]+alpha[7]*fUpwind_l[9]+fUpwind_l[7]*alpha[9]+alpha[3]*fUpwind_l[4]+fUpwind_l[3]*alpha[4])); 
  Ghat_l[12] = 0.01178511301977579*(13.41640786499874*alpha[22]*fUpwind_l[59]+13.41640786499874*(alpha[27]*fUpwind_l[55]+alpha[14]*fUpwind_l[54]+alpha[13]*fUpwind_l[53])+13.41640786499874*(alpha[21]*fUpwind_l[51]+alpha[20]*fUpwind_l[50]+alpha[5]*fUpwind_l[49])+13.41640786499874*alpha[12]*fUpwind_l[48]+15.0*(alpha[19]*fUpwind_l[31]+alpha[26]*fUpwind_l[30]+alpha[11]*fUpwind_l[29]+alpha[10]*fUpwind_l[28]+alpha[8]*fUpwind_l[27]+fUpwind_l[8]*alpha[27]+alpha[18]*fUpwind_l[25]+alpha[17]*fUpwind_l[24]+alpha[4]*fUpwind_l[23]+alpha[16]*fUpwind_l[22]+fUpwind_l[16]*alpha[22]+alpha[3]*fUpwind_l[21]+fUpwind_l[3]*alpha[21]+alpha[2]*fUpwind_l[20]+fUpwind_l[2]*alpha[20]+alpha[9]*fUpwind_l[15]+alpha[7]*fUpwind_l[14]+fUpwind_l[7]*alpha[14]+alpha[6]*fUpwind_l[13]+fUpwind_l[6]*alpha[13]+alpha[0]*fUpwind_l[12]+fUpwind_l[0]*alpha[12]+alpha[1]*fUpwind_l[5]+fUpwind_l[1]*alpha[5])); 
  Ghat_l[13] = 0.01178511301977579*(13.41640786499874*alpha[21]*fUpwind_l[59]+13.41640786499874*(alpha[14]*fUpwind_l[55]+alpha[27]*fUpwind_l[54]+alpha[12]*fUpwind_l[53])+13.41640786499874*(alpha[22]*fUpwind_l[51]+alpha[5]*fUpwind_l[50]+alpha[20]*fUpwind_l[49])+13.41640786499874*alpha[13]*fUpwind_l[48]+15.0*(alpha[18]*fUpwind_l[31]+alpha[11]*fUpwind_l[30]+alpha[26]*fUpwind_l[29]+alpha[9]*fUpwind_l[28]+alpha[7]*fUpwind_l[27]+fUpwind_l[7]*alpha[27]+alpha[19]*fUpwind_l[25]+alpha[4]*fUpwind_l[24]+alpha[17]*fUpwind_l[23]+alpha[3]*fUpwind_l[22]+fUpwind_l[3]*alpha[22]+alpha[16]*fUpwind_l[21]+fUpwind_l[16]*alpha[21]+alpha[1]*fUpwind_l[20]+fUpwind_l[1]*alpha[20]+alpha[10]*fUpwind_l[15]+alpha[8]*fUpwind_l[14]+fUpwind_l[8]*alpha[14]+alpha[0]*fUpwind_l[13]+fUpwind_l[0]*alpha[13]+alpha[6]*fUpwind_l[12]+fUpwind_l[6]*alpha[12]+alpha[2]*fUpwind_l[5]+fUpwind_l[2]*alpha[5])); 
  Ghat_l[14] = 0.01178511301977579*(13.41640786499874*alpha[20]*fUpwind_l[59]+13.41640786499874*(alpha[13]*fUpwind_l[55]+alpha[12]*fUpwind_l[54]+alpha[27]*fUpwind_l[53])+13.41640786499874*(alpha[5]*fUpwind_l[51]+alpha[22]*fUpwind_l[50]+alpha[21]*fUpwind_l[49])+13.41640786499874*alpha[14]*fUpwind_l[48]+15.0*(alpha[17]*fUpwind_l[31]+alpha[10]*fUpwind_l[30]+alpha[9]*fUpwind_l[29]+alpha[26]*fUpwind_l[28]+alpha[6]*fUpwind_l[27]+fUpwind_l[6]*alpha[27]+alpha[4]*fUpwind_l[25]+alpha[19]*fUpwind_l[24]+alpha[18]*fUpwind_l[23]+alpha[2]*fUpwind_l[22]+fUpwind_l[2]*alpha[22]+alpha[1]*fUpwind_l[21]+fUpwind_l[1]*alpha[21]+alpha[16]*fUpwind_l[20]+fUpwind_l[16]*alpha[20]+alpha[11]*fUpwind_l[15]+alpha[0]*fUpwind_l[14]+fUpwind_l[0]*alpha[14]+alpha[8]*fUpwind_l[13]+fUpwind_l[8]*alpha[13]+alpha[7]*fUpwind_l[12]+fUpwind_l[7]*alpha[12]+alpha[3]*fUpwind_l[5]+fUpwind_l[3]*alpha[5])); 
  Ghat_l[15] = 0.01178511301977579*(13.41640786499874*alpha[27]*fUpwind_l[63]+13.41640786499874*(alpha[22]*fUpwind_l[62]+alpha[21]*fUpwind_l[61]+alpha[20]*fUpwind_l[60])+13.41640786499874*(alpha[14]*fUpwind_l[58]+alpha[13]*fUpwind_l[57]+alpha[12]*fUpwind_l[56])+13.41640786499874*alpha[5]*fUpwind_l[52]+13.41640786499874*alpha[26]*fUpwind_l[47]+13.41640786499874*(alpha[19]*fUpwind_l[46]+alpha[18]*fUpwind_l[45]+alpha[17]*fUpwind_l[44])+13.41640786499874*(alpha[11]*fUpwind_l[42]+alpha[10]*fUpwind_l[41]+alpha[9]*fUpwind_l[40])+13.41640786499874*alpha[4]*fUpwind_l[36]+15.0*(alpha[16]*fUpwind_l[31]+alpha[8]*fUpwind_l[30]+alpha[7]*fUpwind_l[29]+alpha[6]*fUpwind_l[28]+alpha[26]*fUpwind_l[27]+fUpwind_l[26]*alpha[27]+alpha[3]*fUpwind_l[25]+alpha[2]*fUpwind_l[24]+alpha[1]*fUpwind_l[23]+alpha[19]*fUpwind_l[22]+fUpwind_l[19]*alpha[22]+alpha[18]*fUpwind_l[21]+fUpwind_l[18]*alpha[21]+alpha[17]*fUpwind_l[20]+fUpwind_l[17]*alpha[20]+alpha[0]*fUpwind_l[15]+alpha[11]*fUpwind_l[14]+fUpwind_l[11]*alpha[14]+alpha[10]*fUpwind_l[13]+fUpwind_l[10]*alpha[13]+alpha[9]*fUpwind_l[12]+fUpwind_l[9]*alpha[12]+alpha[4]*fUpwind_l[5]+fUpwind_l[4]*alpha[5])); 
  Ghat_l[16] = 0.1767766952966368*(alpha[5]*fUpwind_l[27]+fUpwind_l[5]*alpha[27]+alpha[4]*fUpwind_l[26]+fUpwind_l[4]*alpha[26]+alpha[12]*fUpwind_l[22]+fUpwind_l[12]*alpha[22]+alpha[13]*fUpwind_l[21]+fUpwind_l[13]*alpha[21]+alpha[14]*fUpwind_l[20]+fUpwind_l[14]*alpha[20]+alpha[9]*fUpwind_l[19]+fUpwind_l[9]*alpha[19]+alpha[10]*fUpwind_l[18]+fUpwind_l[10]*alpha[18]+alpha[11]*fUpwind_l[17]+fUpwind_l[11]*alpha[17]+alpha[0]*fUpwind_l[16]+fUpwind_l[0]*alpha[16]+alpha[1]*fUpwind_l[8]+fUpwind_l[1]*alpha[8]+alpha[2]*fUpwind_l[7]+fUpwind_l[2]*alpha[7]+alpha[3]*fUpwind_l[6]+fUpwind_l[3]*alpha[6]); 
  Ghat_l[17] = 0.01178511301977579*(13.41640786499874*alpha[11]*fUpwind_l[43]+13.41640786499874*(alpha[18]*fUpwind_l[39]+alpha[19]*fUpwind_l[38]+alpha[4]*fUpwind_l[37])+13.41640786499874*(alpha[26]*fUpwind_l[35]+alpha[9]*fUpwind_l[34]+alpha[10]*fUpwind_l[33])+13.41640786499874*alpha[17]*fUpwind_l[32]+15.0*(alpha[14]*fUpwind_l[31]+alpha[21]*fUpwind_l[30]+alpha[22]*fUpwind_l[29]+alpha[5]*fUpwind_l[28]+fUpwind_l[25]*alpha[27]+alpha[3]*fUpwind_l[26]+fUpwind_l[3]*alpha[26]+alpha[12]*fUpwind_l[24]+alpha[13]*fUpwind_l[23]+fUpwind_l[15]*alpha[20]+alpha[7]*fUpwind_l[19]+fUpwind_l[7]*alpha[19]+alpha[8]*fUpwind_l[18]+fUpwind_l[8]*alpha[18]+alpha[0]*fUpwind_l[17]+fUpwind_l[0]*alpha[17]+alpha[11]*fUpwind_l[16]+fUpwind_l[11]*alpha[16]+alpha[1]*fUpwind_l[10]+fUpwind_l[1]*alpha[10]+alpha[2]*fUpwind_l[9]+fUpwind_l[2]*alpha[9]+alpha[4]*fUpwind_l[6]+fUpwind_l[4]*alpha[6])); 
  Ghat_l[18] = 0.01178511301977579*(13.41640786499874*alpha[10]*fUpwind_l[43]+13.41640786499874*(alpha[17]*fUpwind_l[39]+alpha[4]*fUpwind_l[38]+alpha[19]*fUpwind_l[37])+13.41640786499874*(alpha[9]*fUpwind_l[35]+alpha[26]*fUpwind_l[34]+alpha[11]*fUpwind_l[33])+13.41640786499874*alpha[18]*fUpwind_l[32]+15.0*(alpha[13]*fUpwind_l[31]+alpha[20]*fUpwind_l[30]+alpha[5]*fUpwind_l[29]+alpha[22]*fUpwind_l[28]+fUpwind_l[24]*alpha[27]+alpha[2]*fUpwind_l[26]+fUpwind_l[2]*alpha[26]+alpha[12]*fUpwind_l[25]+alpha[14]*fUpwind_l[23]+fUpwind_l[15]*alpha[21]+alpha[6]*fUpwind_l[19]+fUpwind_l[6]*alpha[19]+alpha[0]*fUpwind_l[18]+fUpwind_l[0]*alpha[18]+alpha[8]*fUpwind_l[17]+fUpwind_l[8]*alpha[17]+alpha[10]*fUpwind_l[16]+fUpwind_l[10]*alpha[16]+alpha[1]*fUpwind_l[11]+fUpwind_l[1]*alpha[11]+alpha[3]*fUpwind_l[9]+fUpwind_l[3]*alpha[9]+alpha[4]*fUpwind_l[7]+fUpwind_l[4]*alpha[7])); 
  Ghat_l[19] = 0.01178511301977579*(13.41640786499874*alpha[9]*fUpwind_l[43]+13.41640786499874*(alpha[4]*fUpwind_l[39]+alpha[17]*fUpwind_l[38]+alpha[18]*fUpwind_l[37])+13.41640786499874*(alpha[10]*fUpwind_l[35]+alpha[11]*fUpwind_l[34]+alpha[26]*fUpwind_l[33])+13.41640786499874*alpha[19]*fUpwind_l[32]+15.0*(alpha[12]*fUpwind_l[31]+alpha[5]*fUpwind_l[30]+alpha[20]*fUpwind_l[29]+alpha[21]*fUpwind_l[28]+fUpwind_l[23]*alpha[27]+alpha[1]*fUpwind_l[26]+fUpwind_l[1]*alpha[26]+alpha[13]*fUpwind_l[25]+alpha[14]*fUpwind_l[24]+fUpwind_l[15]*alpha[22]+alpha[0]*fUpwind_l[19]+fUpwind_l[0]*alpha[19]+alpha[6]*fUpwind_l[18]+fUpwind_l[6]*alpha[18]+alpha[7]*fUpwind_l[17]+fUpwind_l[7]*alpha[17]+alpha[9]*fUpwind_l[16]+fUpwind_l[9]*alpha[16]+alpha[2]*fUpwind_l[11]+fUpwind_l[2]*alpha[11]+alpha[3]*fUpwind_l[10]+fUpwind_l[3]*alpha[10]+alpha[4]*fUpwind_l[8]+fUpwind_l[4]*alpha[8])); 
  Ghat_l[20] = 0.01178511301977579*(13.41640786499874*alpha[14]*fUpwind_l[59]+13.41640786499874*(alpha[21]*fUpwind_l[55]+alpha[22]*fUpwind_l[54]+alpha[5]*fUpwind_l[53])+13.41640786499874*(alpha[27]*fUpwind_l[51]+alpha[12]*fUpwind_l[50]+alpha[13]*fUpwind_l[49])+13.41640786499874*alpha[20]*fUpwind_l[48]+15.0*(alpha[11]*fUpwind_l[31]+alpha[18]*fUpwind_l[30]+alpha[19]*fUpwind_l[29]+alpha[4]*fUpwind_l[28]+alpha[3]*fUpwind_l[27]+fUpwind_l[3]*alpha[27]+fUpwind_l[25]*alpha[26]+alpha[9]*fUpwind_l[24]+alpha[10]*fUpwind_l[23]+alpha[7]*fUpwind_l[22]+fUpwind_l[7]*alpha[22]+alpha[8]*fUpwind_l[21]+fUpwind_l[8]*alpha[21]+alpha[0]*fUpwind_l[20]+fUpwind_l[0]*alpha[20]+fUpwind_l[15]*alpha[17]+alpha[14]*fUpwind_l[16]+fUpwind_l[14]*alpha[16]+alpha[1]*fUpwind_l[13]+fUpwind_l[1]*alpha[13]+alpha[2]*fUpwind_l[12]+fUpwind_l[2]*alpha[12]+alpha[5]*fUpwind_l[6]+fUpwind_l[5]*alpha[6])); 
  Ghat_l[21] = 0.01178511301977579*(13.41640786499874*alpha[13]*fUpwind_l[59]+13.41640786499874*(alpha[20]*fUpwind_l[55]+alpha[5]*fUpwind_l[54]+alpha[22]*fUpwind_l[53])+13.41640786499874*(alpha[12]*fUpwind_l[51]+alpha[27]*fUpwind_l[50]+alpha[14]*fUpwind_l[49])+13.41640786499874*alpha[21]*fUpwind_l[48]+15.0*(alpha[10]*fUpwind_l[31]+alpha[17]*fUpwind_l[30]+alpha[4]*fUpwind_l[29]+alpha[19]*fUpwind_l[28]+alpha[2]*fUpwind_l[27]+fUpwind_l[2]*alpha[27]+fUpwind_l[24]*alpha[26]+alpha[9]*fUpwind_l[25]+alpha[11]*fUpwind_l[23]+alpha[6]*fUpwind_l[22]+fUpwind_l[6]*alpha[22]+alpha[0]*fUpwind_l[21]+fUpwind_l[0]*alpha[21]+alpha[8]*fUpwind_l[20]+fUpwind_l[8]*alpha[20]+fUpwind_l[15]*alpha[18]+alpha[13]*fUpwind_l[16]+fUpwind_l[13]*alpha[16]+alpha[1]*fUpwind_l[14]+fUpwind_l[1]*alpha[14]+alpha[3]*fUpwind_l[12]+fUpwind_l[3]*alpha[12]+alpha[5]*fUpwind_l[7]+fUpwind_l[5]*alpha[7])); 
  Ghat_l[22] = 0.01178511301977579*(13.41640786499874*alpha[12]*fUpwind_l[59]+13.41640786499874*(alpha[5]*fUpwind_l[55]+alpha[20]*fUpwind_l[54]+alpha[21]*fUpwind_l[53])+13.41640786499874*(alpha[13]*fUpwind_l[51]+alpha[14]*fUpwind_l[50]+alpha[27]*fUpwind_l[49])+13.41640786499874*alpha[22]*fUpwind_l[48]+15.0*(alpha[9]*fUpwind_l[31]+alpha[4]*fUpwind_l[30]+alpha[17]*fUpwind_l[29]+alpha[18]*fUpwind_l[28]+alpha[1]*fUpwind_l[27]+fUpwind_l[1]*alpha[27]+fUpwind_l[23]*alpha[26]+alpha[10]*fUpwind_l[25]+alpha[11]*fUpwind_l[24]+alpha[0]*fUpwind_l[22]+fUpwind_l[0]*alpha[22]+alpha[6]*fUpwind_l[21]+fUpwind_l[6]*alpha[21]+alpha[7]*fUpwind_l[20]+fUpwind_l[7]*alpha[20]+fUpwind_l[15]*alpha[19]+alpha[12]*fUpwind_l[16]+fUpwind_l[12]*alpha[16]+alpha[2]*fUpwind_l[14]+fUpwind_l[2]*alpha[14]+alpha[3]*fUpwind_l[13]+fUpwind_l[3]*alpha[13]+alpha[5]*fUpwind_l[8]+fUpwind_l[5]*alpha[8])); 
  Ghat_l[23] = 0.01178511301977579*(13.41640786499874*alpha[22]*fUpwind_l[63]+13.41640786499874*(alpha[27]*fUpwind_l[62]+alpha[14]*fUpwind_l[61]+alpha[13]*fUpwind_l[60])+13.41640786499874*(alpha[21]*fUpwind_l[58]+alpha[20]*fUpwind_l[57]+alpha[5]*fUpwind_l[56])+13.41640786499874*alpha[12]*fUpwind_l[52]+13.41640786499874*alpha[19]*fUpwind_l[47]+13.41640786499874*(alpha[26]*fUpwind_l[46]+alpha[11]*fUpwind_l[45]+alpha[10]*fUpwind_l[44])+13.41640786499874*(alpha[18]*fUpwind_l[42]+alpha[17]*fUpwind_l[41]+alpha[4]*fUpwind_l[40])+13.41640786499874*alpha[9]*fUpwind_l[36]+15.0*(alpha[8]*fUpwind_l[31]+alpha[16]*fUpwind_l[30]+alpha[3]*fUpwind_l[29]+alpha[2]*fUpwind_l[28]+alpha[19]*fUpwind_l[27]+fUpwind_l[19]*alpha[27]+alpha[22]*fUpwind_l[26]+fUpwind_l[22]*alpha[26]+alpha[7]*fUpwind_l[25]+alpha[6]*fUpwind_l[24]+alpha[0]*fUpwind_l[23]+alpha[11]*fUpwind_l[21]+fUpwind_l[11]*alpha[21]+alpha[10]*fUpwind_l[20]+fUpwind_l[10]*alpha[20]+alpha[14]*fUpwind_l[18]+fUpwind_l[14]*alpha[18]+alpha[13]*fUpwind_l[17]+fUpwind_l[13]*alpha[17]+alpha[1]*fUpwind_l[15]+alpha[4]*fUpwind_l[12]+fUpwind_l[4]*alpha[12]+alpha[5]*fUpwind_l[9]+fUpwind_l[5]*alpha[9])); 
  Ghat_l[24] = 0.01178511301977579*(13.41640786499874*alpha[21]*fUpwind_l[63]+13.41640786499874*(alpha[14]*fUpwind_l[62]+alpha[27]*fUpwind_l[61]+alpha[12]*fUpwind_l[60])+13.41640786499874*(alpha[22]*fUpwind_l[58]+alpha[5]*fUpwind_l[57]+alpha[20]*fUpwind_l[56])+13.41640786499874*alpha[13]*fUpwind_l[52]+13.41640786499874*alpha[18]*fUpwind_l[47]+13.41640786499874*(alpha[11]*fUpwind_l[46]+alpha[26]*fUpwind_l[45]+alpha[9]*fUpwind_l[44])+13.41640786499874*(alpha[19]*fUpwind_l[42]+alpha[4]*fUpwind_l[41]+alpha[17]*fUpwind_l[40])+13.41640786499874*alpha[10]*fUpwind_l[36]+15.0*(alpha[7]*fUpwind_l[31]+alpha[3]*fUpwind_l[30]+alpha[16]*fUpwind_l[29]+alpha[1]*fUpwind_l[28]+alpha[18]*fUpwind_l[27]+fUpwind_l[18]*alpha[27]+alpha[21]*fUpwind_l[26]+fUpwind_l[21]*alpha[26]+alpha[8]*fUpwind_l[25]+alpha[0]*fUpwind_l[24]+alpha[6]*fUpwind_l[23]+alpha[11]*fUpwind_l[22]+fUpwind_l[11]*alpha[22]+alpha[9]*fUpwind_l[20]+fUpwind_l[9]*alpha[20]+alpha[14]*fUpwind_l[19]+fUpwind_l[14]*alpha[19]+alpha[12]*fUpwind_l[17]+fUpwind_l[12]*alpha[17]+alpha[2]*fUpwind_l[15]+alpha[4]*fUpwind_l[13]+fUpwind_l[4]*alpha[13]+alpha[5]*fUpwind_l[10]+fUpwind_l[5]*alpha[10])); 
  Ghat_l[25] = 0.01178511301977579*(13.41640786499874*alpha[20]*fUpwind_l[63]+13.41640786499874*(alpha[13]*fUpwind_l[62]+alpha[12]*fUpwind_l[61]+alpha[27]*fUpwind_l[60])+13.41640786499874*(alpha[5]*fUpwind_l[58]+alpha[22]*fUpwind_l[57]+alpha[21]*fUpwind_l[56])+13.41640786499874*alpha[14]*fUpwind_l[52]+13.41640786499874*alpha[17]*fUpwind_l[47]+13.41640786499874*(alpha[10]*fUpwind_l[46]+alpha[9]*fUpwind_l[45]+alpha[26]*fUpwind_l[44])+13.41640786499874*(alpha[4]*fUpwind_l[42]+alpha[19]*fUpwind_l[41]+alpha[18]*fUpwind_l[40])+13.41640786499874*alpha[11]*fUpwind_l[36]+15.0*(alpha[6]*fUpwind_l[31]+alpha[2]*fUpwind_l[30]+alpha[1]*fUpwind_l[29]+alpha[16]*fUpwind_l[28]+alpha[17]*fUpwind_l[27]+fUpwind_l[17]*alpha[27]+alpha[20]*fUpwind_l[26]+fUpwind_l[20]*alpha[26]+alpha[0]*fUpwind_l[25]+alpha[8]*fUpwind_l[24]+alpha[7]*fUpwind_l[23]+alpha[10]*fUpwind_l[22]+fUpwind_l[10]*alpha[22]+alpha[9]*fUpwind_l[21]+fUpwind_l[9]*alpha[21]+alpha[13]*fUpwind_l[19]+fUpwind_l[13]*alpha[19]+alpha[12]*fUpwind_l[18]+fUpwind_l[12]*alpha[18]+alpha[3]*fUpwind_l[15]+alpha[4]*fUpwind_l[14]+fUpwind_l[4]*alpha[14]+alpha[5]*fUpwind_l[11]+fUpwind_l[5]*alpha[11])); 
  Ghat_l[26] = 0.01178511301977579*(13.41640786499874*alpha[4]*fUpwind_l[43]+13.41640786499874*(alpha[9]*fUpwind_l[39]+alpha[10]*fUpwind_l[38]+alpha[11]*fUpwind_l[37])+13.41640786499874*(alpha[17]*fUpwind_l[35]+alpha[18]*fUpwind_l[34]+alpha[19]*fUpwind_l[33])+13.41640786499874*alpha[26]*fUpwind_l[32]+15.0*(alpha[5]*fUpwind_l[31]+alpha[12]*fUpwind_l[30]+alpha[13]*fUpwind_l[29]+alpha[14]*fUpwind_l[28]+fUpwind_l[15]*alpha[27]+alpha[0]*fUpwind_l[26]+fUpwind_l[0]*alpha[26]+alpha[20]*fUpwind_l[25]+alpha[21]*fUpwind_l[24]+alpha[22]*fUpwind_l[23]+alpha[1]*fUpwind_l[19]+fUpwind_l[1]*alpha[19]+alpha[2]*fUpwind_l[18]+fUpwind_l[2]*alpha[18]+alpha[3]*fUpwind_l[17]+fUpwind_l[3]*alpha[17]+alpha[4]*fUpwind_l[16]+fUpwind_l[4]*alpha[16]+alpha[6]*fUpwind_l[11]+fUpwind_l[6]*alpha[11]+alpha[7]*fUpwind_l[10]+fUpwind_l[7]*alpha[10]+alpha[8]*fUpwind_l[9]+fUpwind_l[8]*alpha[9])); 
  Ghat_l[27] = 0.01178511301977579*(13.41640786499874*alpha[5]*fUpwind_l[59]+13.41640786499874*(alpha[12]*fUpwind_l[55]+alpha[13]*fUpwind_l[54]+alpha[14]*fUpwind_l[53])+13.41640786499874*(alpha[20]*fUpwind_l[51]+alpha[21]*fUpwind_l[50]+alpha[22]*fUpwind_l[49])+13.41640786499874*alpha[27]*fUpwind_l[48]+15.0*(alpha[4]*fUpwind_l[31]+alpha[9]*fUpwind_l[30]+alpha[10]*fUpwind_l[29]+alpha[11]*fUpwind_l[28]+alpha[0]*fUpwind_l[27]+fUpwind_l[0]*alpha[27]+fUpwind_l[15]*alpha[26]+alpha[17]*fUpwind_l[25]+alpha[18]*fUpwind_l[24]+alpha[19]*fUpwind_l[23]+alpha[1]*fUpwind_l[22]+fUpwind_l[1]*alpha[22]+alpha[2]*fUpwind_l[21]+fUpwind_l[2]*alpha[21]+alpha[3]*fUpwind_l[20]+fUpwind_l[3]*alpha[20]+alpha[5]*fUpwind_l[16]+fUpwind_l[5]*alpha[16]+alpha[6]*fUpwind_l[14]+fUpwind_l[6]*alpha[14]+alpha[7]*fUpwind_l[13]+fUpwind_l[7]*alpha[13]+alpha[8]*fUpwind_l[12]+fUpwind_l[8]*alpha[12])); 
  Ghat_l[28] = 0.01178511301977579*(13.41640786499874*alpha[14]*fUpwind_l[63]+13.41640786499874*(alpha[21]*fUpwind_l[62]+alpha[22]*fUpwind_l[61]+alpha[5]*fUpwind_l[60])+13.41640786499874*(alpha[27]*fUpwind_l[58]+alpha[12]*fUpwind_l[57]+alpha[13]*fUpwind_l[56])+13.41640786499874*alpha[20]*fUpwind_l[52]+13.41640786499874*alpha[11]*fUpwind_l[47]+13.41640786499874*(alpha[18]*fUpwind_l[46]+alpha[19]*fUpwind_l[45]+alpha[4]*fUpwind_l[44])+13.41640786499874*(alpha[26]*fUpwind_l[42]+alpha[9]*fUpwind_l[41]+alpha[10]*fUpwind_l[40])+13.41640786499874*alpha[17]*fUpwind_l[36]+15.0*(alpha[3]*fUpwind_l[31]+alpha[7]*fUpwind_l[30]+alpha[8]*fUpwind_l[29]+alpha[0]*fUpwind_l[28]+alpha[11]*fUpwind_l[27]+fUpwind_l[11]*alpha[27]+alpha[14]*fUpwind_l[26]+fUpwind_l[14]*alpha[26]+alpha[16]*fUpwind_l[25]+alpha[1]*fUpwind_l[24]+alpha[2]*fUpwind_l[23]+alpha[18]*fUpwind_l[22]+fUpwind_l[18]*alpha[22]+alpha[19]*fUpwind_l[21]+fUpwind_l[19]*alpha[21]+alpha[4]*fUpwind_l[20]+fUpwind_l[4]*alpha[20]+alpha[5]*fUpwind_l[17]+fUpwind_l[5]*alpha[17]+alpha[6]*fUpwind_l[15]+alpha[9]*fUpwind_l[13]+fUpwind_l[9]*alpha[13]+alpha[10]*fUpwind_l[12]+fUpwind_l[10]*alpha[12])); 
  Ghat_l[29] = 0.01178511301977579*(13.41640786499874*alpha[13]*fUpwind_l[63]+13.41640786499874*(alpha[20]*fUpwind_l[62]+alpha[5]*fUpwind_l[61]+alpha[22]*fUpwind_l[60])+13.41640786499874*(alpha[12]*fUpwind_l[58]+alpha[27]*fUpwind_l[57]+alpha[14]*fUpwind_l[56])+13.41640786499874*alpha[21]*fUpwind_l[52]+13.41640786499874*alpha[10]*fUpwind_l[47]+13.41640786499874*(alpha[17]*fUpwind_l[46]+alpha[4]*fUpwind_l[45]+alpha[19]*fUpwind_l[44])+13.41640786499874*(alpha[9]*fUpwind_l[42]+alpha[26]*fUpwind_l[41]+alpha[11]*fUpwind_l[40])+13.41640786499874*alpha[18]*fUpwind_l[36]+15.0*(alpha[2]*fUpwind_l[31]+alpha[6]*fUpwind_l[30]+alpha[0]*fUpwind_l[29]+alpha[8]*fUpwind_l[28]+alpha[10]*fUpwind_l[27]+fUpwind_l[10]*alpha[27]+alpha[13]*fUpwind_l[26]+fUpwind_l[13]*alpha[26]+alpha[1]*fUpwind_l[25]+alpha[16]*fUpwind_l[24]+alpha[3]*fUpwind_l[23]+alpha[17]*fUpwind_l[22]+fUpwind_l[17]*alpha[22]+alpha[4]*fUpwind_l[21]+fUpwind_l[4]*alpha[21]+alpha[19]*fUpwind_l[20]+fUpwind_l[19]*alpha[20]+alpha[5]*fUpwind_l[18]+fUpwind_l[5]*alpha[18]+alpha[7]*fUpwind_l[15]+alpha[9]*fUpwind_l[14]+fUpwind_l[9]*alpha[14]+alpha[11]*fUpwind_l[12]+fUpwind_l[11]*alpha[12])); 
  Ghat_l[30] = 0.01178511301977579*(13.41640786499874*alpha[12]*fUpwind_l[63]+13.41640786499874*(alpha[5]*fUpwind_l[62]+alpha[20]*fUpwind_l[61]+alpha[21]*fUpwind_l[60])+13.41640786499874*(alpha[13]*fUpwind_l[58]+alpha[14]*fUpwind_l[57]+alpha[27]*fUpwind_l[56])+13.41640786499874*alpha[22]*fUpwind_l[52]+13.41640786499874*alpha[9]*fUpwind_l[47]+13.41640786499874*(alpha[4]*fUpwind_l[46]+alpha[17]*fUpwind_l[45]+alpha[18]*fUpwind_l[44])+13.41640786499874*(alpha[10]*fUpwind_l[42]+alpha[11]*fUpwind_l[41]+alpha[26]*fUpwind_l[40])+13.41640786499874*alpha[19]*fUpwind_l[36]+15.0*(alpha[1]*fUpwind_l[31]+alpha[0]*fUpwind_l[30]+alpha[6]*fUpwind_l[29]+alpha[7]*fUpwind_l[28]+alpha[9]*fUpwind_l[27]+fUpwind_l[9]*alpha[27]+alpha[12]*fUpwind_l[26]+fUpwind_l[12]*alpha[26]+alpha[2]*fUpwind_l[25]+alpha[3]*fUpwind_l[24]+alpha[16]*fUpwind_l[23]+alpha[4]*fUpwind_l[22]+fUpwind_l[4]*alpha[22]+alpha[17]*fUpwind_l[21]+fUpwind_l[17]*alpha[21]+alpha[18]*fUpwind_l[20]+fUpwind_l[18]*alpha[20]+alpha[5]*fUpwind_l[19]+fUpwind_l[5]*alpha[19]+alpha[8]*fUpwind_l[15]+alpha[10]*fUpwind_l[14]+fUpwind_l[10]*alpha[14]+alpha[11]*fUpwind_l[13]+fUpwind_l[11]*alpha[13])); 
  Ghat_l[31] = 0.01178511301977579*(13.41640786499874*alpha[5]*fUpwind_l[63]+13.41640786499874*(alpha[12]*fUpwind_l[62]+alpha[13]*fUpwind_l[61]+alpha[14]*fUpwind_l[60])+13.41640786499874*(alpha[20]*fUpwind_l[58]+alpha[21]*fUpwind_l[57]+alpha[22]*fUpwind_l[56])+13.41640786499874*alpha[27]*fUpwind_l[52]+13.41640786499874*alpha[4]*fUpwind_l[47]+13.41640786499874*(alpha[9]*fUpwind_l[46]+alpha[10]*fUpwind_l[45]+alpha[11]*fUpwind_l[44])+13.41640786499874*(alpha[17]*fUpwind_l[42]+alpha[18]*fUpwind_l[41]+alpha[19]*fUpwind_l[40])+13.41640786499874*alpha[26]*fUpwind_l[36]+15.0*(alpha[0]*fUpwind_l[31]+alpha[1]*fUpwind_l[30]+alpha[2]*fUpwind_l[29]+alpha[3]*fUpwind_l[28]+alpha[4]*fUpwind_l[27]+fUpwind_l[4]*alpha[27]+alpha[5]*fUpwind_l[26]+fUpwind_l[5]*alpha[26]+alpha[6]*fUpwind_l[25]+alpha[7]*fUpwind_l[24]+alpha[8]*fUpwind_l[23]+alpha[9]*fUpwind_l[22]+fUpwind_l[9]*alpha[22]+alpha[10]*fUpwind_l[21]+fUpwind_l[10]*alpha[21]+alpha[11]*fUpwind_l[20]+fUpwind_l[11]*alpha[20]+alpha[12]*fUpwind_l[19]+fUpwind_l[12]*alpha[19]+alpha[13]*fUpwind_l[18]+fUpwind_l[13]*alpha[18]+alpha[14]*fUpwind_l[17]+fUpwind_l[14]*alpha[17]+fUpwind_l[15]*alpha[16])); 
  Ghat_l[32] = 0.01178511301977579*(15.0*alpha[27]*fUpwind_l[47]+15.0*(alpha[22]*fUpwind_l[46]+alpha[21]*fUpwind_l[45]+alpha[20]*fUpwind_l[44]+alpha[16]*fUpwind_l[43])+15.0*(alpha[14]*fUpwind_l[42]+alpha[13]*fUpwind_l[41]+alpha[12]*fUpwind_l[40]+alpha[8]*fUpwind_l[39]+alpha[7]*fUpwind_l[38]+alpha[6]*fUpwind_l[37])+15.0*(alpha[5]*fUpwind_l[36]+alpha[3]*fUpwind_l[35]+alpha[2]*fUpwind_l[34]+alpha[1]*fUpwind_l[33])+15.0*alpha[0]*fUpwind_l[32]+13.41640786499874*(alpha[26]*fUpwind_l[26]+alpha[19]*fUpwind_l[19]+alpha[18]*fUpwind_l[18]+alpha[17]*fUpwind_l[17]+alpha[11]*fUpwind_l[11]+alpha[10]*fUpwind_l[10]+alpha[9]*fUpwind_l[9]+alpha[4]*fUpwind_l[4])); 
  Ghat_l[33] = 0.01178511301977579*(15.0*alpha[22]*fUpwind_l[47]+15.0*(alpha[27]*fUpwind_l[46]+alpha[14]*fUpwind_l[45]+alpha[13]*fUpwind_l[44]+alpha[8]*fUpwind_l[43])+15.0*(alpha[21]*fUpwind_l[42]+alpha[20]*fUpwind_l[41]+alpha[5]*fUpwind_l[40]+alpha[16]*fUpwind_l[39]+alpha[3]*fUpwind_l[38]+alpha[2]*fUpwind_l[37])+15.0*(alpha[12]*fUpwind_l[36]+alpha[7]*fUpwind_l[35]+alpha[6]*fUpwind_l[34]+alpha[0]*fUpwind_l[33])+15.0*alpha[1]*fUpwind_l[32]+13.41640786499874*(alpha[19]*fUpwind_l[26]+fUpwind_l[19]*alpha[26]+alpha[11]*fUpwind_l[18]+fUpwind_l[11]*alpha[18]+alpha[10]*fUpwind_l[17]+fUpwind_l[10]*alpha[17]+alpha[4]*fUpwind_l[9]+fUpwind_l[4]*alpha[9])); 
  Ghat_l[34] = 0.01178511301977579*(15.0*alpha[21]*fUpwind_l[47]+15.0*(alpha[14]*fUpwind_l[46]+alpha[27]*fUpwind_l[45]+alpha[12]*fUpwind_l[44]+alpha[7]*fUpwind_l[43])+15.0*(alpha[22]*fUpwind_l[42]+alpha[5]*fUpwind_l[41]+alpha[20]*fUpwind_l[40]+alpha[3]*fUpwind_l[39]+alpha[16]*fUpwind_l[38]+alpha[1]*fUpwind_l[37])+15.0*(alpha[13]*fUpwind_l[36]+alpha[8]*fUpwind_l[35]+alpha[0]*fUpwind_l[34]+alpha[6]*fUpwind_l[33])+15.0*alpha[2]*fUpwind_l[32]+13.41640786499874*(alpha[18]*fUpwind_l[26]+fUpwind_l[18]*alpha[26]+alpha[11]*fUpwind_l[19]+fUpwind_l[11]*alpha[19]+alpha[9]*fUpwind_l[17]+fUpwind_l[9]*alpha[17]+alpha[4]*fUpwind_l[10]+fUpwind_l[4]*alpha[10])); 
  Ghat_l[35] = 0.01178511301977579*(15.0*alpha[20]*fUpwind_l[47]+15.0*(alpha[13]*fUpwind_l[46]+alpha[12]*fUpwind_l[45]+alpha[27]*fUpwind_l[44]+alpha[6]*fUpwind_l[43])+15.0*(alpha[5]*fUpwind_l[42]+alpha[22]*fUpwind_l[41]+alpha[21]*fUpwind_l[40]+alpha[2]*fUpwind_l[39]+alpha[1]*fUpwind_l[38]+alpha[16]*fUpwind_l[37])+15.0*(alpha[14]*fUpwind_l[36]+alpha[0]*fUpwind_l[35]+alpha[8]*fUpwind_l[34]+alpha[7]*fUpwind_l[33])+15.0*alpha[3]*fUpwind_l[32]+13.41640786499874*(alpha[17]*fUpwind_l[26]+fUpwind_l[17]*alpha[26]+alpha[10]*fUpwind_l[19]+fUpwind_l[10]*alpha[19]+alpha[9]*fUpwind_l[18]+fUpwind_l[9]*alpha[18]+alpha[4]*fUpwind_l[11]+fUpwind_l[4]*alpha[11])); 
  Ghat_l[36] = 0.01178511301977579*(15.0*alpha[16]*fUpwind_l[47]+15.0*(alpha[8]*fUpwind_l[46]+alpha[7]*fUpwind_l[45]+alpha[6]*fUpwind_l[44]+alpha[27]*fUpwind_l[43])+15.0*(alpha[3]*fUpwind_l[42]+alpha[2]*fUpwind_l[41]+alpha[1]*fUpwind_l[40]+alpha[22]*fUpwind_l[39]+alpha[21]*fUpwind_l[38]+alpha[20]*fUpwind_l[37])+15.0*(alpha[0]*fUpwind_l[36]+alpha[14]*fUpwind_l[35]+alpha[13]*fUpwind_l[34]+alpha[12]*fUpwind_l[33])+15.0*alpha[5]*fUpwind_l[32]+13.41640786499874*(alpha[26]*fUpwind_l[31]+alpha[19]*fUpwind_l[30]+alpha[18]*fUpwind_l[29]+alpha[17]*fUpwind_l[28]+alpha[11]*fUpwind_l[25]+alpha[10]*fUpwind_l[24]+alpha[9]*fUpwind_l[23]+alpha[4]*fUpwind_l[15])); 
  Ghat_l[37] = 0.01178511301977579*(15.0*alpha[14]*fUpwind_l[47]+15.0*(alpha[21]*fUpwind_l[46]+alpha[22]*fUpwind_l[45]+alpha[5]*fUpwind_l[44]+alpha[3]*fUpwind_l[43])+15.0*(alpha[27]*fUpwind_l[42]+alpha[12]*fUpwind_l[41]+alpha[13]*fUpwind_l[40]+alpha[7]*fUpwind_l[39]+alpha[8]*fUpwind_l[38]+alpha[0]*fUpwind_l[37])+15.0*(alpha[20]*fUpwind_l[36]+alpha[16]*fUpwind_l[35]+alpha[1]*fUpwind_l[34]+alpha[2]*fUpwind_l[33])+15.0*alpha[6]*fUpwind_l[32]+13.41640786499874*(alpha[11]*fUpwind_l[26]+fUpwind_l[11]*alpha[26]+alpha[18]*fUpwind_l[19]+fUpwind_l[18]*alpha[19]+alpha[4]*fUpwind_l[17]+fUpwind_l[4]*alpha[17]+alpha[9]*fUpwind_l[10]+fUpwind_l[9]*alpha[10])); 
  Ghat_l[38] = 0.01178511301977579*(15.0*alpha[13]*fUpwind_l[47]+15.0*(alpha[20]*fUpwind_l[46]+alpha[5]*fUpwind_l[45]+alpha[22]*fUpwind_l[44]+alpha[2]*fUpwind_l[43])+15.0*(alpha[12]*fUpwind_l[42]+alpha[27]*fUpwind_l[41]+alpha[14]*fUpwind_l[40]+alpha[6]*fUpwind_l[39]+alpha[0]*fUpwind_l[38]+alpha[8]*fUpwind_l[37])+15.0*(alpha[21]*fUpwind_l[36]+alpha[1]*fUpwind_l[35]+alpha[16]*fUpwind_l[34]+alpha[3]*fUpwind_l[33])+15.0*alpha[7]*fUpwind_l[32]+13.41640786499874*(alpha[10]*fUpwind_l[26]+fUpwind_l[10]*alpha[26]+alpha[17]*fUpwind_l[19]+fUpwind_l[17]*alpha[19]+alpha[4]*fUpwind_l[18]+fUpwind_l[4]*alpha[18]+alpha[9]*fUpwind_l[11]+fUpwind_l[9]*alpha[11])); 
  Ghat_l[39] = 0.01178511301977579*(15.0*alpha[12]*fUpwind_l[47]+15.0*(alpha[5]*fUpwind_l[46]+alpha[20]*fUpwind_l[45]+alpha[21]*fUpwind_l[44]+alpha[1]*fUpwind_l[43])+15.0*(alpha[13]*fUpwind_l[42]+alpha[14]*fUpwind_l[41]+alpha[27]*fUpwind_l[40]+alpha[0]*fUpwind_l[39]+alpha[6]*fUpwind_l[38]+alpha[7]*fUpwind_l[37])+15.0*(alpha[22]*fUpwind_l[36]+alpha[2]*fUpwind_l[35]+alpha[3]*fUpwind_l[34]+alpha[16]*fUpwind_l[33])+15.0*alpha[8]*fUpwind_l[32]+13.41640786499874*(alpha[9]*fUpwind_l[26]+fUpwind_l[9]*alpha[26]+alpha[4]*fUpwind_l[19]+fUpwind_l[4]*alpha[19]+alpha[17]*fUpwind_l[18]+fUpwind_l[17]*alpha[18]+alpha[10]*fUpwind_l[11]+fUpwind_l[10]*alpha[11])); 
  Ghat_l[40] = 0.01178511301977579*(15.0*alpha[8]*fUpwind_l[47]+15.0*(alpha[16]*fUpwind_l[46]+alpha[3]*fUpwind_l[45]+alpha[2]*fUpwind_l[44]+alpha[22]*fUpwind_l[43])+15.0*(alpha[7]*fUpwind_l[42]+alpha[6]*fUpwind_l[41]+alpha[0]*fUpwind_l[40]+alpha[27]*fUpwind_l[39]+alpha[14]*fUpwind_l[38]+alpha[13]*fUpwind_l[37])+15.0*(alpha[1]*fUpwind_l[36]+alpha[21]*fUpwind_l[35]+alpha[20]*fUpwind_l[34]+alpha[5]*fUpwind_l[33])+15.0*alpha[12]*fUpwind_l[32]+13.41640786499874*(alpha[19]*fUpwind_l[31]+alpha[26]*fUpwind_l[30]+alpha[11]*fUpwind_l[29]+alpha[10]*fUpwind_l[28]+alpha[18]*fUpwind_l[25]+alpha[17]*fUpwind_l[24]+alpha[4]*fUpwind_l[23]+alpha[9]*fUpwind_l[15])); 
  Ghat_l[41] = 0.01178511301977579*(15.0*alpha[7]*fUpwind_l[47]+15.0*(alpha[3]*fUpwind_l[46]+alpha[16]*fUpwind_l[45]+alpha[1]*fUpwind_l[44]+alpha[21]*fUpwind_l[43])+15.0*(alpha[8]*fUpwind_l[42]+alpha[0]*fUpwind_l[41]+alpha[6]*fUpwind_l[40]+alpha[14]*fUpwind_l[39]+alpha[27]*fUpwind_l[38]+alpha[12]*fUpwind_l[37])+15.0*(alpha[2]*fUpwind_l[36]+alpha[22]*fUpwind_l[35]+alpha[5]*fUpwind_l[34]+alpha[20]*fUpwind_l[33])+15.0*alpha[13]*fUpwind_l[32]+13.41640786499874*(alpha[18]*fUpwind_l[31]+alpha[11]*fUpwind_l[30]+alpha[26]*fUpwind_l[29]+alpha[9]*fUpwind_l[28]+alpha[19]*fUpwind_l[25]+alpha[4]*fUpwind_l[24]+alpha[17]*fUpwind_l[23]+alpha[10]*fUpwind_l[15])); 
  Ghat_l[42] = 0.01178511301977579*(15.0*alpha[6]*fUpwind_l[47]+15.0*(alpha[2]*fUpwind_l[46]+alpha[1]*fUpwind_l[45]+alpha[16]*fUpwind_l[44]+alpha[20]*fUpwind_l[43])+15.0*(alpha[0]*fUpwind_l[42]+alpha[8]*fUpwind_l[41]+alpha[7]*fUpwind_l[40]+alpha[13]*fUpwind_l[39]+alpha[12]*fUpwind_l[38]+alpha[27]*fUpwind_l[37])+15.0*(alpha[3]*fUpwind_l[36]+alpha[5]*fUpwind_l[35]+alpha[22]*fUpwind_l[34]+alpha[21]*fUpwind_l[33])+15.0*alpha[14]*fUpwind_l[32]+13.41640786499874*(alpha[17]*fUpwind_l[31]+alpha[10]*fUpwind_l[30]+alpha[9]*fUpwind_l[29]+alpha[26]*fUpwind_l[28]+alpha[4]*fUpwind_l[25]+alpha[19]*fUpwind_l[24]+alpha[18]*fUpwind_l[23]+alpha[11]*fUpwind_l[15])); 
  Ghat_l[43] = 0.01178511301977579*(15.0*alpha[5]*fUpwind_l[47]+15.0*(alpha[12]*fUpwind_l[46]+alpha[13]*fUpwind_l[45]+alpha[14]*fUpwind_l[44]+alpha[0]*fUpwind_l[43])+15.0*(alpha[20]*fUpwind_l[42]+alpha[21]*fUpwind_l[41]+alpha[22]*fUpwind_l[40]+alpha[1]*fUpwind_l[39]+alpha[2]*fUpwind_l[38]+alpha[3]*fUpwind_l[37])+15.0*(alpha[27]*fUpwind_l[36]+alpha[6]*fUpwind_l[35]+alpha[7]*fUpwind_l[34]+alpha[8]*fUpwind_l[33])+15.0*alpha[16]*fUpwind_l[32]+13.41640786499874*(alpha[4]*fUpwind_l[26]+fUpwind_l[4]*alpha[26]+alpha[9]*fUpwind_l[19]+fUpwind_l[9]*alpha[19]+alpha[10]*fUpwind_l[18]+fUpwind_l[10]*alpha[18]+alpha[11]*fUpwind_l[17]+fUpwind_l[11]*alpha[17])); 
  Ghat_l[44] = 0.01178511301977579*(15.0*alpha[3]*fUpwind_l[47]+15.0*(alpha[7]*fUpwind_l[46]+alpha[8]*fUpwind_l[45]+alpha[0]*fUpwind_l[44]+alpha[14]*fUpwind_l[43])+15.0*(alpha[16]*fUpwind_l[42]+alpha[1]*fUpwind_l[41]+alpha[2]*fUpwind_l[40]+alpha[21]*fUpwind_l[39]+alpha[22]*fUpwind_l[38]+alpha[5]*fUpwind_l[37])+15.0*(alpha[6]*fUpwind_l[36]+alpha[27]*fUpwind_l[35]+alpha[12]*fUpwind_l[34]+alpha[13]*fUpwind_l[33])+15.0*alpha[20]*fUpwind_l[32]+13.41640786499874*(alpha[11]*fUpwind_l[31]+alpha[18]*fUpwind_l[30]+alpha[19]*fUpwind_l[29]+alpha[4]*fUpwind_l[28]+fUpwind_l[25]*alpha[26]+alpha[9]*fUpwind_l[24]+alpha[10]*fUpwind_l[23]+fUpwind_l[15]*alpha[17])); 
  Ghat_l[45] = 0.01178511301977579*(15.0*alpha[2]*fUpwind_l[47]+15.0*(alpha[6]*fUpwind_l[46]+alpha[0]*fUpwind_l[45]+alpha[8]*fUpwind_l[44]+alpha[13]*fUpwind_l[43])+15.0*(alpha[1]*fUpwind_l[42]+alpha[16]*fUpwind_l[41]+alpha[3]*fUpwind_l[40]+alpha[20]*fUpwind_l[39]+alpha[5]*fUpwind_l[38]+alpha[22]*fUpwind_l[37])+15.0*(alpha[7]*fUpwind_l[36]+alpha[12]*fUpwind_l[35]+alpha[27]*fUpwind_l[34]+alpha[14]*fUpwind_l[33])+15.0*alpha[21]*fUpwind_l[32]+13.41640786499874*(alpha[10]*fUpwind_l[31]+alpha[17]*fUpwind_l[30]+alpha[4]*fUpwind_l[29]+alpha[19]*fUpwind_l[28]+fUpwind_l[24]*alpha[26]+alpha[9]*fUpwind_l[25]+alpha[11]*fUpwind_l[23]+fUpwind_l[15]*alpha[18])); 
  Ghat_l[46] = 0.01178511301977579*(15.0*alpha[1]*fUpwind_l[47]+15.0*(alpha[0]*fUpwind_l[46]+alpha[6]*fUpwind_l[45]+alpha[7]*fUpwind_l[44]+alpha[12]*fUpwind_l[43])+15.0*(alpha[2]*fUpwind_l[42]+alpha[3]*fUpwind_l[41]+alpha[16]*fUpwind_l[40]+alpha[5]*fUpwind_l[39]+alpha[20]*fUpwind_l[38]+alpha[21]*fUpwind_l[37])+15.0*(alpha[8]*fUpwind_l[36]+alpha[13]*fUpwind_l[35]+alpha[14]*fUpwind_l[34]+alpha[27]*fUpwind_l[33])+15.0*alpha[22]*fUpwind_l[32]+13.41640786499874*(alpha[9]*fUpwind_l[31]+alpha[4]*fUpwind_l[30]+alpha[17]*fUpwind_l[29]+alpha[18]*fUpwind_l[28]+fUpwind_l[23]*alpha[26]+alpha[10]*fUpwind_l[25]+alpha[11]*fUpwind_l[24]+fUpwind_l[15]*alpha[19])); 
  Ghat_l[47] = 0.01178511301977579*(15.0*alpha[0]*fUpwind_l[47]+15.0*(alpha[1]*fUpwind_l[46]+alpha[2]*fUpwind_l[45]+alpha[3]*fUpwind_l[44]+alpha[5]*fUpwind_l[43])+15.0*(alpha[6]*fUpwind_l[42]+alpha[7]*fUpwind_l[41]+alpha[8]*fUpwind_l[40]+alpha[12]*fUpwind_l[39]+alpha[13]*fUpwind_l[38]+alpha[14]*fUpwind_l[37])+15.0*(alpha[16]*fUpwind_l[36]+alpha[20]*fUpwind_l[35]+alpha[21]*fUpwind_l[34]+alpha[22]*fUpwind_l[33])+15.0*alpha[27]*fUpwind_l[32]+13.41640786499874*(alpha[4]*fUpwind_l[31]+alpha[9]*fUpwind_l[30]+alpha[10]*fUpwind_l[29]+alpha[11]*fUpwind_l[28]+fUpwind_l[15]*alpha[26]+alpha[17]*fUpwind_l[25]+alpha[18]*fUpwind_l[24]+alpha[19]*fUpwind_l[23])); 
  Ghat_l[48] = 0.01178511301977579*(15.0*alpha[26]*fUpwind_l[63]+15.0*(alpha[19]*fUpwind_l[62]+alpha[18]*fUpwind_l[61]+alpha[17]*fUpwind_l[60]+alpha[16]*fUpwind_l[59])+15.0*(alpha[11]*fUpwind_l[58]+alpha[10]*fUpwind_l[57]+alpha[9]*fUpwind_l[56]+alpha[8]*fUpwind_l[55]+alpha[7]*fUpwind_l[54]+alpha[6]*fUpwind_l[53])+15.0*(alpha[4]*fUpwind_l[52]+alpha[3]*fUpwind_l[51]+alpha[2]*fUpwind_l[50]+alpha[1]*fUpwind_l[49])+15.0*alpha[0]*fUpwind_l[48]+13.41640786499874*(alpha[27]*fUpwind_l[27]+alpha[22]*fUpwind_l[22]+alpha[21]*fUpwind_l[21]+alpha[20]*fUpwind_l[20]+alpha[14]*fUpwind_l[14]+alpha[13]*fUpwind_l[13]+alpha[12]*fUpwind_l[12]+alpha[5]*fUpwind_l[5])); 
  Ghat_l[49] = 0.01178511301977579*(15.0*alpha[19]*fUpwind_l[63]+15.0*(alpha[26]*fUpwind_l[62]+alpha[11]*fUpwind_l[61]+alpha[10]*fUpwind_l[60]+alpha[8]*fUpwind_l[59])+15.0*(alpha[18]*fUpwind_l[58]+alpha[17]*fUpwind_l[57]+alpha[4]*fUpwind_l[56]+alpha[16]*fUpwind_l[55]+alpha[3]*fUpwind_l[54]+alpha[2]*fUpwind_l[53])+15.0*(alpha[9]*fUpwind_l[52]+alpha[7]*fUpwind_l[51]+alpha[6]*fUpwind_l[50]+alpha[0]*fUpwind_l[49])+15.0*alpha[1]*fUpwind_l[48]+13.41640786499874*(alpha[22]*fUpwind_l[27]+fUpwind_l[22]*alpha[27]+alpha[14]*fUpwind_l[21]+fUpwind_l[14]*alpha[21]+alpha[13]*fUpwind_l[20]+fUpwind_l[13]*alpha[20]+alpha[5]*fUpwind_l[12]+fUpwind_l[5]*alpha[12])); 
  Ghat_l[50] = 0.01178511301977579*(15.0*alpha[18]*fUpwind_l[63]+15.0*(alpha[11]*fUpwind_l[62]+alpha[26]*fUpwind_l[61]+alpha[9]*fUpwind_l[60]+alpha[7]*fUpwind_l[59])+15.0*(alpha[19]*fUpwind_l[58]+alpha[4]*fUpwind_l[57]+alpha[17]*fUpwind_l[56]+alpha[3]*fUpwind_l[55]+alpha[16]*fUpwind_l[54]+alpha[1]*fUpwind_l[53])+15.0*(alpha[10]*fUpwind_l[52]+alpha[8]*fUpwind_l[51]+alpha[0]*fUpwind_l[50]+alpha[6]*fUpwind_l[49])+15.0*alpha[2]*fUpwind_l[48]+13.41640786499874*(alpha[21]*fUpwind_l[27]+fUpwind_l[21]*alpha[27]+alpha[14]*fUpwind_l[22]+fUpwind_l[14]*alpha[22]+alpha[12]*fUpwind_l[20]+fUpwind_l[12]*alpha[20]+alpha[5]*fUpwind_l[13]+fUpwind_l[5]*alpha[13])); 
  Ghat_l[51] = 0.01178511301977579*(15.0*alpha[17]*fUpwind_l[63]+15.0*(alpha[10]*fUpwind_l[62]+alpha[9]*fUpwind_l[61]+alpha[26]*fUpwind_l[60]+alpha[6]*fUpwind_l[59])+15.0*(alpha[4]*fUpwind_l[58]+alpha[19]*fUpwind_l[57]+alpha[18]*fUpwind_l[56]+alpha[2]*fUpwind_l[55]+alpha[1]*fUpwind_l[54]+alpha[16]*fUpwind_l[53])+15.0*(alpha[11]*fUpwind_l[52]+alpha[0]*fUpwind_l[51]+alpha[8]*fUpwind_l[50]+alpha[7]*fUpwind_l[49])+15.0*alpha[3]*fUpwind_l[48]+13.41640786499874*(alpha[20]*fUpwind_l[27]+fUpwind_l[20]*alpha[27]+alpha[13]*fUpwind_l[22]+fUpwind_l[13]*alpha[22]+alpha[12]*fUpwind_l[21]+fUpwind_l[12]*alpha[21]+alpha[5]*fUpwind_l[14]+fUpwind_l[5]*alpha[14])); 
  Ghat_l[52] = 0.01178511301977579*(15.0*alpha[16]*fUpwind_l[63]+15.0*(alpha[8]*fUpwind_l[62]+alpha[7]*fUpwind_l[61]+alpha[6]*fUpwind_l[60]+alpha[26]*fUpwind_l[59])+15.0*(alpha[3]*fUpwind_l[58]+alpha[2]*fUpwind_l[57]+alpha[1]*fUpwind_l[56]+alpha[19]*fUpwind_l[55]+alpha[18]*fUpwind_l[54]+alpha[17]*fUpwind_l[53])+15.0*(alpha[0]*fUpwind_l[52]+alpha[11]*fUpwind_l[51]+alpha[10]*fUpwind_l[50]+alpha[9]*fUpwind_l[49])+15.0*alpha[4]*fUpwind_l[48]+13.41640786499874*(alpha[27]*fUpwind_l[31]+alpha[22]*fUpwind_l[30]+alpha[21]*fUpwind_l[29]+alpha[20]*fUpwind_l[28]+alpha[14]*fUpwind_l[25]+alpha[13]*fUpwind_l[24]+alpha[12]*fUpwind_l[23]+alpha[5]*fUpwind_l[15])); 
  Ghat_l[53] = 0.01178511301977579*(15.0*alpha[11]*fUpwind_l[63]+15.0*(alpha[18]*fUpwind_l[62]+alpha[19]*fUpwind_l[61]+alpha[4]*fUpwind_l[60]+alpha[3]*fUpwind_l[59])+15.0*(alpha[26]*fUpwind_l[58]+alpha[9]*fUpwind_l[57]+alpha[10]*fUpwind_l[56]+alpha[7]*fUpwind_l[55]+alpha[8]*fUpwind_l[54]+alpha[0]*fUpwind_l[53])+15.0*(alpha[17]*fUpwind_l[52]+alpha[16]*fUpwind_l[51]+alpha[1]*fUpwind_l[50]+alpha[2]*fUpwind_l[49])+15.0*alpha[6]*fUpwind_l[48]+13.41640786499874*(alpha[14]*fUpwind_l[27]+fUpwind_l[14]*alpha[27]+alpha[21]*fUpwind_l[22]+fUpwind_l[21]*alpha[22]+alpha[5]*fUpwind_l[20]+fUpwind_l[5]*alpha[20]+alpha[12]*fUpwind_l[13]+fUpwind_l[12]*alpha[13])); 
  Ghat_l[54] = 0.01178511301977579*(15.0*alpha[10]*fUpwind_l[63]+15.0*(alpha[17]*fUpwind_l[62]+alpha[4]*fUpwind_l[61]+alpha[19]*fUpwind_l[60]+alpha[2]*fUpwind_l[59])+15.0*(alpha[9]*fUpwind_l[58]+alpha[26]*fUpwind_l[57]+alpha[11]*fUpwind_l[56]+alpha[6]*fUpwind_l[55]+alpha[0]*fUpwind_l[54]+alpha[8]*fUpwind_l[53])+15.0*(alpha[18]*fUpwind_l[52]+alpha[1]*fUpwind_l[51]+alpha[16]*fUpwind_l[50]+alpha[3]*fUpwind_l[49])+15.0*alpha[7]*fUpwind_l[48]+13.41640786499874*(alpha[13]*fUpwind_l[27]+fUpwind_l[13]*alpha[27]+alpha[20]*fUpwind_l[22]+fUpwind_l[20]*alpha[22]+alpha[5]*fUpwind_l[21]+fUpwind_l[5]*alpha[21]+alpha[12]*fUpwind_l[14]+fUpwind_l[12]*alpha[14])); 
  Ghat_l[55] = 0.01178511301977579*(15.0*alpha[9]*fUpwind_l[63]+15.0*(alpha[4]*fUpwind_l[62]+alpha[17]*fUpwind_l[61]+alpha[18]*fUpwind_l[60]+alpha[1]*fUpwind_l[59])+15.0*(alpha[10]*fUpwind_l[58]+alpha[11]*fUpwind_l[57]+alpha[26]*fUpwind_l[56]+alpha[0]*fUpwind_l[55]+alpha[6]*fUpwind_l[54]+alpha[7]*fUpwind_l[53])+15.0*(alpha[19]*fUpwind_l[52]+alpha[2]*fUpwind_l[51]+alpha[3]*fUpwind_l[50]+alpha[16]*fUpwind_l[49])+15.0*alpha[8]*fUpwind_l[48]+13.41640786499874*(alpha[12]*fUpwind_l[27]+fUpwind_l[12]*alpha[27]+alpha[5]*fUpwind_l[22]+fUpwind_l[5]*alpha[22]+alpha[20]*fUpwind_l[21]+fUpwind_l[20]*alpha[21]+alpha[13]*fUpwind_l[14]+fUpwind_l[13]*alpha[14])); 
  Ghat_l[56] = 0.01178511301977579*(15.0*alpha[8]*fUpwind_l[63]+15.0*(alpha[16]*fUpwind_l[62]+alpha[3]*fUpwind_l[61]+alpha[2]*fUpwind_l[60]+alpha[19]*fUpwind_l[59])+15.0*(alpha[7]*fUpwind_l[58]+alpha[6]*fUpwind_l[57]+alpha[0]*fUpwind_l[56]+alpha[26]*fUpwind_l[55]+alpha[11]*fUpwind_l[54]+alpha[10]*fUpwind_l[53])+15.0*(alpha[1]*fUpwind_l[52]+alpha[18]*fUpwind_l[51]+alpha[17]*fUpwind_l[50]+alpha[4]*fUpwind_l[49])+15.0*alpha[9]*fUpwind_l[48]+13.41640786499874*(alpha[22]*fUpwind_l[31]+alpha[27]*fUpwind_l[30]+alpha[14]*fUpwind_l[29]+alpha[13]*fUpwind_l[28]+alpha[21]*fUpwind_l[25]+alpha[20]*fUpwind_l[24]+alpha[5]*fUpwind_l[23]+alpha[12]*fUpwind_l[15])); 
  Ghat_l[57] = 0.01178511301977579*(15.0*alpha[7]*fUpwind_l[63]+15.0*(alpha[3]*fUpwind_l[62]+alpha[16]*fUpwind_l[61]+alpha[1]*fUpwind_l[60]+alpha[18]*fUpwind_l[59])+15.0*(alpha[8]*fUpwind_l[58]+alpha[0]*fUpwind_l[57]+alpha[6]*fUpwind_l[56]+alpha[11]*fUpwind_l[55]+alpha[26]*fUpwind_l[54]+alpha[9]*fUpwind_l[53])+15.0*(alpha[2]*fUpwind_l[52]+alpha[19]*fUpwind_l[51]+alpha[4]*fUpwind_l[50]+alpha[17]*fUpwind_l[49])+15.0*alpha[10]*fUpwind_l[48]+13.41640786499874*(alpha[21]*fUpwind_l[31]+alpha[14]*fUpwind_l[30]+alpha[27]*fUpwind_l[29]+alpha[12]*fUpwind_l[28]+alpha[22]*fUpwind_l[25]+alpha[5]*fUpwind_l[24]+alpha[20]*fUpwind_l[23]+alpha[13]*fUpwind_l[15])); 
  Ghat_l[58] = 0.01178511301977579*(15.0*alpha[6]*fUpwind_l[63]+15.0*(alpha[2]*fUpwind_l[62]+alpha[1]*fUpwind_l[61]+alpha[16]*fUpwind_l[60]+alpha[17]*fUpwind_l[59])+15.0*(alpha[0]*fUpwind_l[58]+alpha[8]*fUpwind_l[57]+alpha[7]*fUpwind_l[56]+alpha[10]*fUpwind_l[55]+alpha[9]*fUpwind_l[54]+alpha[26]*fUpwind_l[53])+15.0*(alpha[3]*fUpwind_l[52]+alpha[4]*fUpwind_l[51]+alpha[19]*fUpwind_l[50]+alpha[18]*fUpwind_l[49])+15.0*alpha[11]*fUpwind_l[48]+13.41640786499874*(alpha[20]*fUpwind_l[31]+alpha[13]*fUpwind_l[30]+alpha[12]*fUpwind_l[29]+alpha[27]*fUpwind_l[28]+alpha[5]*fUpwind_l[25]+alpha[22]*fUpwind_l[24]+alpha[21]*fUpwind_l[23]+alpha[14]*fUpwind_l[15])); 
  Ghat_l[59] = 0.01178511301977579*(15.0*alpha[4]*fUpwind_l[63]+15.0*(alpha[9]*fUpwind_l[62]+alpha[10]*fUpwind_l[61]+alpha[11]*fUpwind_l[60]+alpha[0]*fUpwind_l[59])+15.0*(alpha[17]*fUpwind_l[58]+alpha[18]*fUpwind_l[57]+alpha[19]*fUpwind_l[56]+alpha[1]*fUpwind_l[55]+alpha[2]*fUpwind_l[54]+alpha[3]*fUpwind_l[53])+15.0*(alpha[26]*fUpwind_l[52]+alpha[6]*fUpwind_l[51]+alpha[7]*fUpwind_l[50]+alpha[8]*fUpwind_l[49])+15.0*alpha[16]*fUpwind_l[48]+13.41640786499874*(alpha[5]*fUpwind_l[27]+fUpwind_l[5]*alpha[27]+alpha[12]*fUpwind_l[22]+fUpwind_l[12]*alpha[22]+alpha[13]*fUpwind_l[21]+fUpwind_l[13]*alpha[21]+alpha[14]*fUpwind_l[20]+fUpwind_l[14]*alpha[20])); 
  Ghat_l[60] = 0.01178511301977579*(15.0*alpha[3]*fUpwind_l[63]+15.0*(alpha[7]*fUpwind_l[62]+alpha[8]*fUpwind_l[61]+alpha[0]*fUpwind_l[60]+alpha[11]*fUpwind_l[59])+15.0*(alpha[16]*fUpwind_l[58]+alpha[1]*fUpwind_l[57]+alpha[2]*fUpwind_l[56]+alpha[18]*fUpwind_l[55]+alpha[19]*fUpwind_l[54]+alpha[4]*fUpwind_l[53])+15.0*(alpha[6]*fUpwind_l[52]+alpha[26]*fUpwind_l[51]+alpha[9]*fUpwind_l[50]+alpha[10]*fUpwind_l[49])+15.0*alpha[17]*fUpwind_l[48]+13.41640786499874*(alpha[14]*fUpwind_l[31]+alpha[21]*fUpwind_l[30]+alpha[22]*fUpwind_l[29]+alpha[5]*fUpwind_l[28]+fUpwind_l[25]*alpha[27]+alpha[12]*fUpwind_l[24]+alpha[13]*fUpwind_l[23]+fUpwind_l[15]*alpha[20])); 
  Ghat_l[61] = 0.01178511301977579*(15.0*alpha[2]*fUpwind_l[63]+15.0*(alpha[6]*fUpwind_l[62]+alpha[0]*fUpwind_l[61]+alpha[8]*fUpwind_l[60]+alpha[10]*fUpwind_l[59])+15.0*(alpha[1]*fUpwind_l[58]+alpha[16]*fUpwind_l[57]+alpha[3]*fUpwind_l[56]+alpha[17]*fUpwind_l[55]+alpha[4]*fUpwind_l[54]+alpha[19]*fUpwind_l[53])+15.0*(alpha[7]*fUpwind_l[52]+alpha[9]*fUpwind_l[51]+alpha[26]*fUpwind_l[50]+alpha[11]*fUpwind_l[49])+15.0*alpha[18]*fUpwind_l[48]+13.41640786499874*(alpha[13]*fUpwind_l[31]+alpha[20]*fUpwind_l[30]+alpha[5]*fUpwind_l[29]+alpha[22]*fUpwind_l[28]+fUpwind_l[24]*alpha[27]+alpha[12]*fUpwind_l[25]+alpha[14]*fUpwind_l[23]+fUpwind_l[15]*alpha[21])); 
  Ghat_l[62] = 0.01178511301977579*(15.0*alpha[1]*fUpwind_l[63]+15.0*(alpha[0]*fUpwind_l[62]+alpha[6]*fUpwind_l[61]+alpha[7]*fUpwind_l[60]+alpha[9]*fUpwind_l[59])+15.0*(alpha[2]*fUpwind_l[58]+alpha[3]*fUpwind_l[57]+alpha[16]*fUpwind_l[56]+alpha[4]*fUpwind_l[55]+alpha[17]*fUpwind_l[54]+alpha[18]*fUpwind_l[53])+15.0*(alpha[8]*fUpwind_l[52]+alpha[10]*fUpwind_l[51]+alpha[11]*fUpwind_l[50]+alpha[26]*fUpwind_l[49])+15.0*alpha[19]*fUpwind_l[48]+13.41640786499874*(alpha[12]*fUpwind_l[31]+alpha[5]*fUpwind_l[30]+alpha[20]*fUpwind_l[29]+alpha[21]*fUpwind_l[28]+fUpwind_l[23]*alpha[27]+alpha[13]*fUpwind_l[25]+alpha[14]*fUpwind_l[24]+fUpwind_l[15]*alpha[22])); 
  Ghat_l[63] = 0.01178511301977579*(15.0*alpha[0]*fUpwind_l[63]+15.0*(alpha[1]*fUpwind_l[62]+alpha[2]*fUpwind_l[61]+alpha[3]*fUpwind_l[60]+alpha[4]*fUpwind_l[59])+15.0*(alpha[6]*fUpwind_l[58]+alpha[7]*fUpwind_l[57]+alpha[8]*fUpwind_l[56]+alpha[9]*fUpwind_l[55]+alpha[10]*fUpwind_l[54]+alpha[11]*fUpwind_l[53])+15.0*(alpha[16]*fUpwind_l[52]+alpha[17]*fUpwind_l[51]+alpha[18]*fUpwind_l[50]+alpha[19]*fUpwind_l[49])+15.0*alpha[26]*fUpwind_l[48]+13.41640786499874*(alpha[5]*fUpwind_l[31]+alpha[12]*fUpwind_l[30]+alpha[13]*fUpwind_l[29]+alpha[14]*fUpwind_l[28]+fUpwind_l[15]*alpha[27]+alpha[20]*fUpwind_l[25]+alpha[21]*fUpwind_l[24]+alpha[22]*fUpwind_l[23])); 

  Ghat_r[0] = 0.1767766952966368*(alpha[27]*fUpwind_r[27]+alpha[26]*fUpwind_r[26]+alpha[22]*fUpwind_r[22]+alpha[21]*fUpwind_r[21]+alpha[20]*fUpwind_r[20]+alpha[19]*fUpwind_r[19]+alpha[18]*fUpwind_r[18]+alpha[17]*fUpwind_r[17]+alpha[16]*fUpwind_r[16]+alpha[14]*fUpwind_r[14]+alpha[13]*fUpwind_r[13]+alpha[12]*fUpwind_r[12]+alpha[11]*fUpwind_r[11]+alpha[10]*fUpwind_r[10]+alpha[9]*fUpwind_r[9]+alpha[8]*fUpwind_r[8]+alpha[7]*fUpwind_r[7]+alpha[6]*fUpwind_r[6]+alpha[5]*fUpwind_r[5]+alpha[4]*fUpwind_r[4]+alpha[3]*fUpwind_r[3]+alpha[2]*fUpwind_r[2]+alpha[1]*fUpwind_r[1]+alpha[0]*fUpwind_r[0]); 
  Ghat_r[1] = 0.1767766952966368*(alpha[22]*fUpwind_r[27]+fUpwind_r[22]*alpha[27]+alpha[19]*fUpwind_r[26]+fUpwind_r[19]*alpha[26]+alpha[14]*fUpwind_r[21]+fUpwind_r[14]*alpha[21]+alpha[13]*fUpwind_r[20]+fUpwind_r[13]*alpha[20]+alpha[11]*fUpwind_r[18]+fUpwind_r[11]*alpha[18]+alpha[10]*fUpwind_r[17]+fUpwind_r[10]*alpha[17]+alpha[8]*fUpwind_r[16]+fUpwind_r[8]*alpha[16]+alpha[5]*fUpwind_r[12]+fUpwind_r[5]*alpha[12]+alpha[4]*fUpwind_r[9]+fUpwind_r[4]*alpha[9]+alpha[3]*fUpwind_r[7]+fUpwind_r[3]*alpha[7]+alpha[2]*fUpwind_r[6]+fUpwind_r[2]*alpha[6]+alpha[0]*fUpwind_r[1]+fUpwind_r[0]*alpha[1]); 
  Ghat_r[2] = 0.1767766952966368*(alpha[21]*fUpwind_r[27]+fUpwind_r[21]*alpha[27]+alpha[18]*fUpwind_r[26]+fUpwind_r[18]*alpha[26]+alpha[14]*fUpwind_r[22]+fUpwind_r[14]*alpha[22]+alpha[12]*fUpwind_r[20]+fUpwind_r[12]*alpha[20]+alpha[11]*fUpwind_r[19]+fUpwind_r[11]*alpha[19]+alpha[9]*fUpwind_r[17]+fUpwind_r[9]*alpha[17]+alpha[7]*fUpwind_r[16]+fUpwind_r[7]*alpha[16]+alpha[5]*fUpwind_r[13]+fUpwind_r[5]*alpha[13]+alpha[4]*fUpwind_r[10]+fUpwind_r[4]*alpha[10]+alpha[3]*fUpwind_r[8]+fUpwind_r[3]*alpha[8]+alpha[1]*fUpwind_r[6]+fUpwind_r[1]*alpha[6]+alpha[0]*fUpwind_r[2]+fUpwind_r[0]*alpha[2]); 
  Ghat_r[3] = 0.1767766952966368*(alpha[20]*fUpwind_r[27]+fUpwind_r[20]*alpha[27]+alpha[17]*fUpwind_r[26]+fUpwind_r[17]*alpha[26]+alpha[13]*fUpwind_r[22]+fUpwind_r[13]*alpha[22]+alpha[12]*fUpwind_r[21]+fUpwind_r[12]*alpha[21]+alpha[10]*fUpwind_r[19]+fUpwind_r[10]*alpha[19]+alpha[9]*fUpwind_r[18]+fUpwind_r[9]*alpha[18]+alpha[6]*fUpwind_r[16]+fUpwind_r[6]*alpha[16]+alpha[5]*fUpwind_r[14]+fUpwind_r[5]*alpha[14]+alpha[4]*fUpwind_r[11]+fUpwind_r[4]*alpha[11]+alpha[2]*fUpwind_r[8]+fUpwind_r[2]*alpha[8]+alpha[1]*fUpwind_r[7]+fUpwind_r[1]*alpha[7]+alpha[0]*fUpwind_r[3]+fUpwind_r[0]*alpha[3]); 
  Ghat_r[4] = 0.01178511301977579*(13.41640786499874*alpha[26]*fUpwind_r[43]+13.41640786499874*(alpha[19]*fUpwind_r[39]+alpha[18]*fUpwind_r[38]+alpha[17]*fUpwind_r[37])+13.41640786499874*(alpha[11]*fUpwind_r[35]+alpha[10]*fUpwind_r[34]+alpha[9]*fUpwind_r[33])+13.41640786499874*alpha[4]*fUpwind_r[32]+15.0*(alpha[27]*fUpwind_r[31]+alpha[22]*fUpwind_r[30]+alpha[21]*fUpwind_r[29]+alpha[20]*fUpwind_r[28]+alpha[16]*fUpwind_r[26]+fUpwind_r[16]*alpha[26]+alpha[14]*fUpwind_r[25]+alpha[13]*fUpwind_r[24]+alpha[12]*fUpwind_r[23]+alpha[8]*fUpwind_r[19]+fUpwind_r[8]*alpha[19]+alpha[7]*fUpwind_r[18]+fUpwind_r[7]*alpha[18]+alpha[6]*fUpwind_r[17]+fUpwind_r[6]*alpha[17]+alpha[5]*fUpwind_r[15]+alpha[3]*fUpwind_r[11]+fUpwind_r[3]*alpha[11]+alpha[2]*fUpwind_r[10]+fUpwind_r[2]*alpha[10]+alpha[1]*fUpwind_r[9]+fUpwind_r[1]*alpha[9]+alpha[0]*fUpwind_r[4]+fUpwind_r[0]*alpha[4])); 
  Ghat_r[5] = 0.01178511301977579*(13.41640786499874*alpha[27]*fUpwind_r[59]+13.41640786499874*(alpha[22]*fUpwind_r[55]+alpha[21]*fUpwind_r[54]+alpha[20]*fUpwind_r[53])+13.41640786499874*(alpha[14]*fUpwind_r[51]+alpha[13]*fUpwind_r[50]+alpha[12]*fUpwind_r[49])+13.41640786499874*alpha[5]*fUpwind_r[48]+15.0*(alpha[26]*fUpwind_r[31]+alpha[19]*fUpwind_r[30]+alpha[18]*fUpwind_r[29]+alpha[17]*fUpwind_r[28]+alpha[16]*fUpwind_r[27]+fUpwind_r[16]*alpha[27]+alpha[11]*fUpwind_r[25]+alpha[10]*fUpwind_r[24]+alpha[9]*fUpwind_r[23]+alpha[8]*fUpwind_r[22]+fUpwind_r[8]*alpha[22]+alpha[7]*fUpwind_r[21]+fUpwind_r[7]*alpha[21]+alpha[6]*fUpwind_r[20]+fUpwind_r[6]*alpha[20]+alpha[4]*fUpwind_r[15]+alpha[3]*fUpwind_r[14]+fUpwind_r[3]*alpha[14]+alpha[2]*fUpwind_r[13]+fUpwind_r[2]*alpha[13]+alpha[1]*fUpwind_r[12]+fUpwind_r[1]*alpha[12]+alpha[0]*fUpwind_r[5]+fUpwind_r[0]*alpha[5])); 
  Ghat_r[6] = 0.1767766952966368*(alpha[14]*fUpwind_r[27]+fUpwind_r[14]*alpha[27]+alpha[11]*fUpwind_r[26]+fUpwind_r[11]*alpha[26]+alpha[21]*fUpwind_r[22]+fUpwind_r[21]*alpha[22]+alpha[5]*fUpwind_r[20]+fUpwind_r[5]*alpha[20]+alpha[18]*fUpwind_r[19]+fUpwind_r[18]*alpha[19]+alpha[4]*fUpwind_r[17]+fUpwind_r[4]*alpha[17]+alpha[3]*fUpwind_r[16]+fUpwind_r[3]*alpha[16]+alpha[12]*fUpwind_r[13]+fUpwind_r[12]*alpha[13]+alpha[9]*fUpwind_r[10]+fUpwind_r[9]*alpha[10]+alpha[7]*fUpwind_r[8]+fUpwind_r[7]*alpha[8]+alpha[0]*fUpwind_r[6]+fUpwind_r[0]*alpha[6]+alpha[1]*fUpwind_r[2]+fUpwind_r[1]*alpha[2]); 
  Ghat_r[7] = 0.1767766952966368*(alpha[13]*fUpwind_r[27]+fUpwind_r[13]*alpha[27]+alpha[10]*fUpwind_r[26]+fUpwind_r[10]*alpha[26]+alpha[20]*fUpwind_r[22]+fUpwind_r[20]*alpha[22]+alpha[5]*fUpwind_r[21]+fUpwind_r[5]*alpha[21]+alpha[17]*fUpwind_r[19]+fUpwind_r[17]*alpha[19]+alpha[4]*fUpwind_r[18]+fUpwind_r[4]*alpha[18]+alpha[2]*fUpwind_r[16]+fUpwind_r[2]*alpha[16]+alpha[12]*fUpwind_r[14]+fUpwind_r[12]*alpha[14]+alpha[9]*fUpwind_r[11]+fUpwind_r[9]*alpha[11]+alpha[6]*fUpwind_r[8]+fUpwind_r[6]*alpha[8]+alpha[0]*fUpwind_r[7]+fUpwind_r[0]*alpha[7]+alpha[1]*fUpwind_r[3]+fUpwind_r[1]*alpha[3]); 
  Ghat_r[8] = 0.1767766952966368*(alpha[12]*fUpwind_r[27]+fUpwind_r[12]*alpha[27]+alpha[9]*fUpwind_r[26]+fUpwind_r[9]*alpha[26]+alpha[5]*fUpwind_r[22]+fUpwind_r[5]*alpha[22]+alpha[20]*fUpwind_r[21]+fUpwind_r[20]*alpha[21]+alpha[4]*fUpwind_r[19]+fUpwind_r[4]*alpha[19]+alpha[17]*fUpwind_r[18]+fUpwind_r[17]*alpha[18]+alpha[1]*fUpwind_r[16]+fUpwind_r[1]*alpha[16]+alpha[13]*fUpwind_r[14]+fUpwind_r[13]*alpha[14]+alpha[10]*fUpwind_r[11]+fUpwind_r[10]*alpha[11]+alpha[0]*fUpwind_r[8]+fUpwind_r[0]*alpha[8]+alpha[6]*fUpwind_r[7]+fUpwind_r[6]*alpha[7]+alpha[2]*fUpwind_r[3]+fUpwind_r[2]*alpha[3]); 
  Ghat_r[9] = 0.01178511301977579*(13.41640786499874*alpha[19]*fUpwind_r[43]+13.41640786499874*(alpha[26]*fUpwind_r[39]+alpha[11]*fUpwind_r[38]+alpha[10]*fUpwind_r[37])+13.41640786499874*(alpha[18]*fUpwind_r[35]+alpha[17]*fUpwind_r[34]+alpha[4]*fUpwind_r[33])+13.41640786499874*alpha[9]*fUpwind_r[32]+15.0*(alpha[22]*fUpwind_r[31]+alpha[27]*fUpwind_r[30]+alpha[14]*fUpwind_r[29]+alpha[13]*fUpwind_r[28]+alpha[8]*fUpwind_r[26]+fUpwind_r[8]*alpha[26]+alpha[21]*fUpwind_r[25]+alpha[20]*fUpwind_r[24]+alpha[5]*fUpwind_r[23]+alpha[16]*fUpwind_r[19]+fUpwind_r[16]*alpha[19]+alpha[3]*fUpwind_r[18]+fUpwind_r[3]*alpha[18]+alpha[2]*fUpwind_r[17]+fUpwind_r[2]*alpha[17]+alpha[12]*fUpwind_r[15]+alpha[7]*fUpwind_r[11]+fUpwind_r[7]*alpha[11]+alpha[6]*fUpwind_r[10]+fUpwind_r[6]*alpha[10]+alpha[0]*fUpwind_r[9]+fUpwind_r[0]*alpha[9]+alpha[1]*fUpwind_r[4]+fUpwind_r[1]*alpha[4])); 
  Ghat_r[10] = 0.01178511301977579*(13.41640786499874*alpha[18]*fUpwind_r[43]+13.41640786499874*(alpha[11]*fUpwind_r[39]+alpha[26]*fUpwind_r[38]+alpha[9]*fUpwind_r[37])+13.41640786499874*(alpha[19]*fUpwind_r[35]+alpha[4]*fUpwind_r[34]+alpha[17]*fUpwind_r[33])+13.41640786499874*alpha[10]*fUpwind_r[32]+15.0*(alpha[21]*fUpwind_r[31]+alpha[14]*fUpwind_r[30]+alpha[27]*fUpwind_r[29]+alpha[12]*fUpwind_r[28]+alpha[7]*fUpwind_r[26]+fUpwind_r[7]*alpha[26]+alpha[22]*fUpwind_r[25]+alpha[5]*fUpwind_r[24]+alpha[20]*fUpwind_r[23]+alpha[3]*fUpwind_r[19]+fUpwind_r[3]*alpha[19]+alpha[16]*fUpwind_r[18]+fUpwind_r[16]*alpha[18]+alpha[1]*fUpwind_r[17]+fUpwind_r[1]*alpha[17]+alpha[13]*fUpwind_r[15]+alpha[8]*fUpwind_r[11]+fUpwind_r[8]*alpha[11]+alpha[0]*fUpwind_r[10]+fUpwind_r[0]*alpha[10]+alpha[6]*fUpwind_r[9]+fUpwind_r[6]*alpha[9]+alpha[2]*fUpwind_r[4]+fUpwind_r[2]*alpha[4])); 
  Ghat_r[11] = 0.01178511301977579*(13.41640786499874*alpha[17]*fUpwind_r[43]+13.41640786499874*(alpha[10]*fUpwind_r[39]+alpha[9]*fUpwind_r[38]+alpha[26]*fUpwind_r[37])+13.41640786499874*(alpha[4]*fUpwind_r[35]+alpha[19]*fUpwind_r[34]+alpha[18]*fUpwind_r[33])+13.41640786499874*alpha[11]*fUpwind_r[32]+15.0*(alpha[20]*fUpwind_r[31]+alpha[13]*fUpwind_r[30]+alpha[12]*fUpwind_r[29]+alpha[27]*fUpwind_r[28]+alpha[6]*fUpwind_r[26]+fUpwind_r[6]*alpha[26]+alpha[5]*fUpwind_r[25]+alpha[22]*fUpwind_r[24]+alpha[21]*fUpwind_r[23]+alpha[2]*fUpwind_r[19]+fUpwind_r[2]*alpha[19]+alpha[1]*fUpwind_r[18]+fUpwind_r[1]*alpha[18]+alpha[16]*fUpwind_r[17]+fUpwind_r[16]*alpha[17]+alpha[14]*fUpwind_r[15]+alpha[0]*fUpwind_r[11]+fUpwind_r[0]*alpha[11]+alpha[8]*fUpwind_r[10]+fUpwind_r[8]*alpha[10]+alpha[7]*fUpwind_r[9]+fUpwind_r[7]*alpha[9]+alpha[3]*fUpwind_r[4]+fUpwind_r[3]*alpha[4])); 
  Ghat_r[12] = 0.01178511301977579*(13.41640786499874*alpha[22]*fUpwind_r[59]+13.41640786499874*(alpha[27]*fUpwind_r[55]+alpha[14]*fUpwind_r[54]+alpha[13]*fUpwind_r[53])+13.41640786499874*(alpha[21]*fUpwind_r[51]+alpha[20]*fUpwind_r[50]+alpha[5]*fUpwind_r[49])+13.41640786499874*alpha[12]*fUpwind_r[48]+15.0*(alpha[19]*fUpwind_r[31]+alpha[26]*fUpwind_r[30]+alpha[11]*fUpwind_r[29]+alpha[10]*fUpwind_r[28]+alpha[8]*fUpwind_r[27]+fUpwind_r[8]*alpha[27]+alpha[18]*fUpwind_r[25]+alpha[17]*fUpwind_r[24]+alpha[4]*fUpwind_r[23]+alpha[16]*fUpwind_r[22]+fUpwind_r[16]*alpha[22]+alpha[3]*fUpwind_r[21]+fUpwind_r[3]*alpha[21]+alpha[2]*fUpwind_r[20]+fUpwind_r[2]*alpha[20]+alpha[9]*fUpwind_r[15]+alpha[7]*fUpwind_r[14]+fUpwind_r[7]*alpha[14]+alpha[6]*fUpwind_r[13]+fUpwind_r[6]*alpha[13]+alpha[0]*fUpwind_r[12]+fUpwind_r[0]*alpha[12]+alpha[1]*fUpwind_r[5]+fUpwind_r[1]*alpha[5])); 
  Ghat_r[13] = 0.01178511301977579*(13.41640786499874*alpha[21]*fUpwind_r[59]+13.41640786499874*(alpha[14]*fUpwind_r[55]+alpha[27]*fUpwind_r[54]+alpha[12]*fUpwind_r[53])+13.41640786499874*(alpha[22]*fUpwind_r[51]+alpha[5]*fUpwind_r[50]+alpha[20]*fUpwind_r[49])+13.41640786499874*alpha[13]*fUpwind_r[48]+15.0*(alpha[18]*fUpwind_r[31]+alpha[11]*fUpwind_r[30]+alpha[26]*fUpwind_r[29]+alpha[9]*fUpwind_r[28]+alpha[7]*fUpwind_r[27]+fUpwind_r[7]*alpha[27]+alpha[19]*fUpwind_r[25]+alpha[4]*fUpwind_r[24]+alpha[17]*fUpwind_r[23]+alpha[3]*fUpwind_r[22]+fUpwind_r[3]*alpha[22]+alpha[16]*fUpwind_r[21]+fUpwind_r[16]*alpha[21]+alpha[1]*fUpwind_r[20]+fUpwind_r[1]*alpha[20]+alpha[10]*fUpwind_r[15]+alpha[8]*fUpwind_r[14]+fUpwind_r[8]*alpha[14]+alpha[0]*fUpwind_r[13]+fUpwind_r[0]*alpha[13]+alpha[6]*fUpwind_r[12]+fUpwind_r[6]*alpha[12]+alpha[2]*fUpwind_r[5]+fUpwind_r[2]*alpha[5])); 
  Ghat_r[14] = 0.01178511301977579*(13.41640786499874*alpha[20]*fUpwind_r[59]+13.41640786499874*(alpha[13]*fUpwind_r[55]+alpha[12]*fUpwind_r[54]+alpha[27]*fUpwind_r[53])+13.41640786499874*(alpha[5]*fUpwind_r[51]+alpha[22]*fUpwind_r[50]+alpha[21]*fUpwind_r[49])+13.41640786499874*alpha[14]*fUpwind_r[48]+15.0*(alpha[17]*fUpwind_r[31]+alpha[10]*fUpwind_r[30]+alpha[9]*fUpwind_r[29]+alpha[26]*fUpwind_r[28]+alpha[6]*fUpwind_r[27]+fUpwind_r[6]*alpha[27]+alpha[4]*fUpwind_r[25]+alpha[19]*fUpwind_r[24]+alpha[18]*fUpwind_r[23]+alpha[2]*fUpwind_r[22]+fUpwind_r[2]*alpha[22]+alpha[1]*fUpwind_r[21]+fUpwind_r[1]*alpha[21]+alpha[16]*fUpwind_r[20]+fUpwind_r[16]*alpha[20]+alpha[11]*fUpwind_r[15]+alpha[0]*fUpwind_r[14]+fUpwind_r[0]*alpha[14]+alpha[8]*fUpwind_r[13]+fUpwind_r[8]*alpha[13]+alpha[7]*fUpwind_r[12]+fUpwind_r[7]*alpha[12]+alpha[3]*fUpwind_r[5]+fUpwind_r[3]*alpha[5])); 
  Ghat_r[15] = 0.01178511301977579*(13.41640786499874*alpha[27]*fUpwind_r[63]+13.41640786499874*(alpha[22]*fUpwind_r[62]+alpha[21]*fUpwind_r[61]+alpha[20]*fUpwind_r[60])+13.41640786499874*(alpha[14]*fUpwind_r[58]+alpha[13]*fUpwind_r[57]+alpha[12]*fUpwind_r[56])+13.41640786499874*alpha[5]*fUpwind_r[52]+13.41640786499874*alpha[26]*fUpwind_r[47]+13.41640786499874*(alpha[19]*fUpwind_r[46]+alpha[18]*fUpwind_r[45]+alpha[17]*fUpwind_r[44])+13.41640786499874*(alpha[11]*fUpwind_r[42]+alpha[10]*fUpwind_r[41]+alpha[9]*fUpwind_r[40])+13.41640786499874*alpha[4]*fUpwind_r[36]+15.0*(alpha[16]*fUpwind_r[31]+alpha[8]*fUpwind_r[30]+alpha[7]*fUpwind_r[29]+alpha[6]*fUpwind_r[28]+alpha[26]*fUpwind_r[27]+fUpwind_r[26]*alpha[27]+alpha[3]*fUpwind_r[25]+alpha[2]*fUpwind_r[24]+alpha[1]*fUpwind_r[23]+alpha[19]*fUpwind_r[22]+fUpwind_r[19]*alpha[22]+alpha[18]*fUpwind_r[21]+fUpwind_r[18]*alpha[21]+alpha[17]*fUpwind_r[20]+fUpwind_r[17]*alpha[20]+alpha[0]*fUpwind_r[15]+alpha[11]*fUpwind_r[14]+fUpwind_r[11]*alpha[14]+alpha[10]*fUpwind_r[13]+fUpwind_r[10]*alpha[13]+alpha[9]*fUpwind_r[12]+fUpwind_r[9]*alpha[12]+alpha[4]*fUpwind_r[5]+fUpwind_r[4]*alpha[5])); 
  Ghat_r[16] = 0.1767766952966368*(alpha[5]*fUpwind_r[27]+fUpwind_r[5]*alpha[27]+alpha[4]*fUpwind_r[26]+fUpwind_r[4]*alpha[26]+alpha[12]*fUpwind_r[22]+fUpwind_r[12]*alpha[22]+alpha[13]*fUpwind_r[21]+fUpwind_r[13]*alpha[21]+alpha[14]*fUpwind_r[20]+fUpwind_r[14]*alpha[20]+alpha[9]*fUpwind_r[19]+fUpwind_r[9]*alpha[19]+alpha[10]*fUpwind_r[18]+fUpwind_r[10]*alpha[18]+alpha[11]*fUpwind_r[17]+fUpwind_r[11]*alpha[17]+alpha[0]*fUpwind_r[16]+fUpwind_r[0]*alpha[16]+alpha[1]*fUpwind_r[8]+fUpwind_r[1]*alpha[8]+alpha[2]*fUpwind_r[7]+fUpwind_r[2]*alpha[7]+alpha[3]*fUpwind_r[6]+fUpwind_r[3]*alpha[6]); 
  Ghat_r[17] = 0.01178511301977579*(13.41640786499874*alpha[11]*fUpwind_r[43]+13.41640786499874*(alpha[18]*fUpwind_r[39]+alpha[19]*fUpwind_r[38]+alpha[4]*fUpwind_r[37])+13.41640786499874*(alpha[26]*fUpwind_r[35]+alpha[9]*fUpwind_r[34]+alpha[10]*fUpwind_r[33])+13.41640786499874*alpha[17]*fUpwind_r[32]+15.0*(alpha[14]*fUpwind_r[31]+alpha[21]*fUpwind_r[30]+alpha[22]*fUpwind_r[29]+alpha[5]*fUpwind_r[28]+fUpwind_r[25]*alpha[27]+alpha[3]*fUpwind_r[26]+fUpwind_r[3]*alpha[26]+alpha[12]*fUpwind_r[24]+alpha[13]*fUpwind_r[23]+fUpwind_r[15]*alpha[20]+alpha[7]*fUpwind_r[19]+fUpwind_r[7]*alpha[19]+alpha[8]*fUpwind_r[18]+fUpwind_r[8]*alpha[18]+alpha[0]*fUpwind_r[17]+fUpwind_r[0]*alpha[17]+alpha[11]*fUpwind_r[16]+fUpwind_r[11]*alpha[16]+alpha[1]*fUpwind_r[10]+fUpwind_r[1]*alpha[10]+alpha[2]*fUpwind_r[9]+fUpwind_r[2]*alpha[9]+alpha[4]*fUpwind_r[6]+fUpwind_r[4]*alpha[6])); 
  Ghat_r[18] = 0.01178511301977579*(13.41640786499874*alpha[10]*fUpwind_r[43]+13.41640786499874*(alpha[17]*fUpwind_r[39]+alpha[4]*fUpwind_r[38]+alpha[19]*fUpwind_r[37])+13.41640786499874*(alpha[9]*fUpwind_r[35]+alpha[26]*fUpwind_r[34]+alpha[11]*fUpwind_r[33])+13.41640786499874*alpha[18]*fUpwind_r[32]+15.0*(alpha[13]*fUpwind_r[31]+alpha[20]*fUpwind_r[30]+alpha[5]*fUpwind_r[29]+alpha[22]*fUpwind_r[28]+fUpwind_r[24]*alpha[27]+alpha[2]*fUpwind_r[26]+fUpwind_r[2]*alpha[26]+alpha[12]*fUpwind_r[25]+alpha[14]*fUpwind_r[23]+fUpwind_r[15]*alpha[21]+alpha[6]*fUpwind_r[19]+fUpwind_r[6]*alpha[19]+alpha[0]*fUpwind_r[18]+fUpwind_r[0]*alpha[18]+alpha[8]*fUpwind_r[17]+fUpwind_r[8]*alpha[17]+alpha[10]*fUpwind_r[16]+fUpwind_r[10]*alpha[16]+alpha[1]*fUpwind_r[11]+fUpwind_r[1]*alpha[11]+alpha[3]*fUpwind_r[9]+fUpwind_r[3]*alpha[9]+alpha[4]*fUpwind_r[7]+fUpwind_r[4]*alpha[7])); 
  Ghat_r[19] = 0.01178511301977579*(13.41640786499874*alpha[9]*fUpwind_r[43]+13.41640786499874*(alpha[4]*fUpwind_r[39]+alpha[17]*fUpwind_r[38]+alpha[18]*fUpwind_r[37])+13.41640786499874*(alpha[10]*fUpwind_r[35]+alpha[11]*fUpwind_r[34]+alpha[26]*fUpwind_r[33])+13.41640786499874*alpha[19]*fUpwind_r[32]+15.0*(alpha[12]*fUpwind_r[31]+alpha[5]*fUpwind_r[30]+alpha[20]*fUpwind_r[29]+alpha[21]*fUpwind_r[28]+fUpwind_r[23]*alpha[27]+alpha[1]*fUpwind_r[26]+fUpwind_r[1]*alpha[26]+alpha[13]*fUpwind_r[25]+alpha[14]*fUpwind_r[24]+fUpwind_r[15]*alpha[22]+alpha[0]*fUpwind_r[19]+fUpwind_r[0]*alpha[19]+alpha[6]*fUpwind_r[18]+fUpwind_r[6]*alpha[18]+alpha[7]*fUpwind_r[17]+fUpwind_r[7]*alpha[17]+alpha[9]*fUpwind_r[16]+fUpwind_r[9]*alpha[16]+alpha[2]*fUpwind_r[11]+fUpwind_r[2]*alpha[11]+alpha[3]*fUpwind_r[10]+fUpwind_r[3]*alpha[10]+alpha[4]*fUpwind_r[8]+fUpwind_r[4]*alpha[8])); 
  Ghat_r[20] = 0.01178511301977579*(13.41640786499874*alpha[14]*fUpwind_r[59]+13.41640786499874*(alpha[21]*fUpwind_r[55]+alpha[22]*fUpwind_r[54]+alpha[5]*fUpwind_r[53])+13.41640786499874*(alpha[27]*fUpwind_r[51]+alpha[12]*fUpwind_r[50]+alpha[13]*fUpwind_r[49])+13.41640786499874*alpha[20]*fUpwind_r[48]+15.0*(alpha[11]*fUpwind_r[31]+alpha[18]*fUpwind_r[30]+alpha[19]*fUpwind_r[29]+alpha[4]*fUpwind_r[28]+alpha[3]*fUpwind_r[27]+fUpwind_r[3]*alpha[27]+fUpwind_r[25]*alpha[26]+alpha[9]*fUpwind_r[24]+alpha[10]*fUpwind_r[23]+alpha[7]*fUpwind_r[22]+fUpwind_r[7]*alpha[22]+alpha[8]*fUpwind_r[21]+fUpwind_r[8]*alpha[21]+alpha[0]*fUpwind_r[20]+fUpwind_r[0]*alpha[20]+fUpwind_r[15]*alpha[17]+alpha[14]*fUpwind_r[16]+fUpwind_r[14]*alpha[16]+alpha[1]*fUpwind_r[13]+fUpwind_r[1]*alpha[13]+alpha[2]*fUpwind_r[12]+fUpwind_r[2]*alpha[12]+alpha[5]*fUpwind_r[6]+fUpwind_r[5]*alpha[6])); 
  Ghat_r[21] = 0.01178511301977579*(13.41640786499874*alpha[13]*fUpwind_r[59]+13.41640786499874*(alpha[20]*fUpwind_r[55]+alpha[5]*fUpwind_r[54]+alpha[22]*fUpwind_r[53])+13.41640786499874*(alpha[12]*fUpwind_r[51]+alpha[27]*fUpwind_r[50]+alpha[14]*fUpwind_r[49])+13.41640786499874*alpha[21]*fUpwind_r[48]+15.0*(alpha[10]*fUpwind_r[31]+alpha[17]*fUpwind_r[30]+alpha[4]*fUpwind_r[29]+alpha[19]*fUpwind_r[28]+alpha[2]*fUpwind_r[27]+fUpwind_r[2]*alpha[27]+fUpwind_r[24]*alpha[26]+alpha[9]*fUpwind_r[25]+alpha[11]*fUpwind_r[23]+alpha[6]*fUpwind_r[22]+fUpwind_r[6]*alpha[22]+alpha[0]*fUpwind_r[21]+fUpwind_r[0]*alpha[21]+alpha[8]*fUpwind_r[20]+fUpwind_r[8]*alpha[20]+fUpwind_r[15]*alpha[18]+alpha[13]*fUpwind_r[16]+fUpwind_r[13]*alpha[16]+alpha[1]*fUpwind_r[14]+fUpwind_r[1]*alpha[14]+alpha[3]*fUpwind_r[12]+fUpwind_r[3]*alpha[12]+alpha[5]*fUpwind_r[7]+fUpwind_r[5]*alpha[7])); 
  Ghat_r[22] = 0.01178511301977579*(13.41640786499874*alpha[12]*fUpwind_r[59]+13.41640786499874*(alpha[5]*fUpwind_r[55]+alpha[20]*fUpwind_r[54]+alpha[21]*fUpwind_r[53])+13.41640786499874*(alpha[13]*fUpwind_r[51]+alpha[14]*fUpwind_r[50]+alpha[27]*fUpwind_r[49])+13.41640786499874*alpha[22]*fUpwind_r[48]+15.0*(alpha[9]*fUpwind_r[31]+alpha[4]*fUpwind_r[30]+alpha[17]*fUpwind_r[29]+alpha[18]*fUpwind_r[28]+alpha[1]*fUpwind_r[27]+fUpwind_r[1]*alpha[27]+fUpwind_r[23]*alpha[26]+alpha[10]*fUpwind_r[25]+alpha[11]*fUpwind_r[24]+alpha[0]*fUpwind_r[22]+fUpwind_r[0]*alpha[22]+alpha[6]*fUpwind_r[21]+fUpwind_r[6]*alpha[21]+alpha[7]*fUpwind_r[20]+fUpwind_r[7]*alpha[20]+fUpwind_r[15]*alpha[19]+alpha[12]*fUpwind_r[16]+fUpwind_r[12]*alpha[16]+alpha[2]*fUpwind_r[14]+fUpwind_r[2]*alpha[14]+alpha[3]*fUpwind_r[13]+fUpwind_r[3]*alpha[13]+alpha[5]*fUpwind_r[8]+fUpwind_r[5]*alpha[8])); 
  Ghat_r[23] = 0.01178511301977579*(13.41640786499874*alpha[22]*fUpwind_r[63]+13.41640786499874*(alpha[27]*fUpwind_r[62]+alpha[14]*fUpwind_r[61]+alpha[13]*fUpwind_r[60])+13.41640786499874*(alpha[21]*fUpwind_r[58]+alpha[20]*fUpwind_r[57]+alpha[5]*fUpwind_r[56])+13.41640786499874*alpha[12]*fUpwind_r[52]+13.41640786499874*alpha[19]*fUpwind_r[47]+13.41640786499874*(alpha[26]*fUpwind_r[46]+alpha[11]*fUpwind_r[45]+alpha[10]*fUpwind_r[44])+13.41640786499874*(alpha[18]*fUpwind_r[42]+alpha[17]*fUpwind_r[41]+alpha[4]*fUpwind_r[40])+13.41640786499874*alpha[9]*fUpwind_r[36]+15.0*(alpha[8]*fUpwind_r[31]+alpha[16]*fUpwind_r[30]+alpha[3]*fUpwind_r[29]+alpha[2]*fUpwind_r[28]+alpha[19]*fUpwind_r[27]+fUpwind_r[19]*alpha[27]+alpha[22]*fUpwind_r[26]+fUpwind_r[22]*alpha[26]+alpha[7]*fUpwind_r[25]+alpha[6]*fUpwind_r[24]+alpha[0]*fUpwind_r[23]+alpha[11]*fUpwind_r[21]+fUpwind_r[11]*alpha[21]+alpha[10]*fUpwind_r[20]+fUpwind_r[10]*alpha[20]+alpha[14]*fUpwind_r[18]+fUpwind_r[14]*alpha[18]+alpha[13]*fUpwind_r[17]+fUpwind_r[13]*alpha[17]+alpha[1]*fUpwind_r[15]+alpha[4]*fUpwind_r[12]+fUpwind_r[4]*alpha[12]+alpha[5]*fUpwind_r[9]+fUpwind_r[5]*alpha[9])); 
  Ghat_r[24] = 0.01178511301977579*(13.41640786499874*alpha[21]*fUpwind_r[63]+13.41640786499874*(alpha[14]*fUpwind_r[62]+alpha[27]*fUpwind_r[61]+alpha[12]*fUpwind_r[60])+13.41640786499874*(alpha[22]*fUpwind_r[58]+alpha[5]*fUpwind_r[57]+alpha[20]*fUpwind_r[56])+13.41640786499874*alpha[13]*fUpwind_r[52]+13.41640786499874*alpha[18]*fUpwind_r[47]+13.41640786499874*(alpha[11]*fUpwind_r[46]+alpha[26]*fUpwind_r[45]+alpha[9]*fUpwind_r[44])+13.41640786499874*(alpha[19]*fUpwind_r[42]+alpha[4]*fUpwind_r[41]+alpha[17]*fUpwind_r[40])+13.41640786499874*alpha[10]*fUpwind_r[36]+15.0*(alpha[7]*fUpwind_r[31]+alpha[3]*fUpwind_r[30]+alpha[16]*fUpwind_r[29]+alpha[1]*fUpwind_r[28]+alpha[18]*fUpwind_r[27]+fUpwind_r[18]*alpha[27]+alpha[21]*fUpwind_r[26]+fUpwind_r[21]*alpha[26]+alpha[8]*fUpwind_r[25]+alpha[0]*fUpwind_r[24]+alpha[6]*fUpwind_r[23]+alpha[11]*fUpwind_r[22]+fUpwind_r[11]*alpha[22]+alpha[9]*fUpwind_r[20]+fUpwind_r[9]*alpha[20]+alpha[14]*fUpwind_r[19]+fUpwind_r[14]*alpha[19]+alpha[12]*fUpwind_r[17]+fUpwind_r[12]*alpha[17]+alpha[2]*fUpwind_r[15]+alpha[4]*fUpwind_r[13]+fUpwind_r[4]*alpha[13]+alpha[5]*fUpwind_r[10]+fUpwind_r[5]*alpha[10])); 
  Ghat_r[25] = 0.01178511301977579*(13.41640786499874*alpha[20]*fUpwind_r[63]+13.41640786499874*(alpha[13]*fUpwind_r[62]+alpha[12]*fUpwind_r[61]+alpha[27]*fUpwind_r[60])+13.41640786499874*(alpha[5]*fUpwind_r[58]+alpha[22]*fUpwind_r[57]+alpha[21]*fUpwind_r[56])+13.41640786499874*alpha[14]*fUpwind_r[52]+13.41640786499874*alpha[17]*fUpwind_r[47]+13.41640786499874*(alpha[10]*fUpwind_r[46]+alpha[9]*fUpwind_r[45]+alpha[26]*fUpwind_r[44])+13.41640786499874*(alpha[4]*fUpwind_r[42]+alpha[19]*fUpwind_r[41]+alpha[18]*fUpwind_r[40])+13.41640786499874*alpha[11]*fUpwind_r[36]+15.0*(alpha[6]*fUpwind_r[31]+alpha[2]*fUpwind_r[30]+alpha[1]*fUpwind_r[29]+alpha[16]*fUpwind_r[28]+alpha[17]*fUpwind_r[27]+fUpwind_r[17]*alpha[27]+alpha[20]*fUpwind_r[26]+fUpwind_r[20]*alpha[26]+alpha[0]*fUpwind_r[25]+alpha[8]*fUpwind_r[24]+alpha[7]*fUpwind_r[23]+alpha[10]*fUpwind_r[22]+fUpwind_r[10]*alpha[22]+alpha[9]*fUpwind_r[21]+fUpwind_r[9]*alpha[21]+alpha[13]*fUpwind_r[19]+fUpwind_r[13]*alpha[19]+alpha[12]*fUpwind_r[18]+fUpwind_r[12]*alpha[18]+alpha[3]*fUpwind_r[15]+alpha[4]*fUpwind_r[14]+fUpwind_r[4]*alpha[14]+alpha[5]*fUpwind_r[11]+fUpwind_r[5]*alpha[11])); 
  Ghat_r[26] = 0.01178511301977579*(13.41640786499874*alpha[4]*fUpwind_r[43]+13.41640786499874*(alpha[9]*fUpwind_r[39]+alpha[10]*fUpwind_r[38]+alpha[11]*fUpwind_r[37])+13.41640786499874*(alpha[17]*fUpwind_r[35]+alpha[18]*fUpwind_r[34]+alpha[19]*fUpwind_r[33])+13.41640786499874*alpha[26]*fUpwind_r[32]+15.0*(alpha[5]*fUpwind_r[31]+alpha[12]*fUpwind_r[30]+alpha[13]*fUpwind_r[29]+alpha[14]*fUpwind_r[28]+fUpwind_r[15]*alpha[27]+alpha[0]*fUpwind_r[26]+fUpwind_r[0]*alpha[26]+alpha[20]*fUpwind_r[25]+alpha[21]*fUpwind_r[24]+alpha[22]*fUpwind_r[23]+alpha[1]*fUpwind_r[19]+fUpwind_r[1]*alpha[19]+alpha[2]*fUpwind_r[18]+fUpwind_r[2]*alpha[18]+alpha[3]*fUpwind_r[17]+fUpwind_r[3]*alpha[17]+alpha[4]*fUpwind_r[16]+fUpwind_r[4]*alpha[16]+alpha[6]*fUpwind_r[11]+fUpwind_r[6]*alpha[11]+alpha[7]*fUpwind_r[10]+fUpwind_r[7]*alpha[10]+alpha[8]*fUpwind_r[9]+fUpwind_r[8]*alpha[9])); 
  Ghat_r[27] = 0.01178511301977579*(13.41640786499874*alpha[5]*fUpwind_r[59]+13.41640786499874*(alpha[12]*fUpwind_r[55]+alpha[13]*fUpwind_r[54]+alpha[14]*fUpwind_r[53])+13.41640786499874*(alpha[20]*fUpwind_r[51]+alpha[21]*fUpwind_r[50]+alpha[22]*fUpwind_r[49])+13.41640786499874*alpha[27]*fUpwind_r[48]+15.0*(alpha[4]*fUpwind_r[31]+alpha[9]*fUpwind_r[30]+alpha[10]*fUpwind_r[29]+alpha[11]*fUpwind_r[28]+alpha[0]*fUpwind_r[27]+fUpwind_r[0]*alpha[27]+fUpwind_r[15]*alpha[26]+alpha[17]*fUpwind_r[25]+alpha[18]*fUpwind_r[24]+alpha[19]*fUpwind_r[23]+alpha[1]*fUpwind_r[22]+fUpwind_r[1]*alpha[22]+alpha[2]*fUpwind_r[21]+fUpwind_r[2]*alpha[21]+alpha[3]*fUpwind_r[20]+fUpwind_r[3]*alpha[20]+alpha[5]*fUpwind_r[16]+fUpwind_r[5]*alpha[16]+alpha[6]*fUpwind_r[14]+fUpwind_r[6]*alpha[14]+alpha[7]*fUpwind_r[13]+fUpwind_r[7]*alpha[13]+alpha[8]*fUpwind_r[12]+fUpwind_r[8]*alpha[12])); 
  Ghat_r[28] = 0.01178511301977579*(13.41640786499874*alpha[14]*fUpwind_r[63]+13.41640786499874*(alpha[21]*fUpwind_r[62]+alpha[22]*fUpwind_r[61]+alpha[5]*fUpwind_r[60])+13.41640786499874*(alpha[27]*fUpwind_r[58]+alpha[12]*fUpwind_r[57]+alpha[13]*fUpwind_r[56])+13.41640786499874*alpha[20]*fUpwind_r[52]+13.41640786499874*alpha[11]*fUpwind_r[47]+13.41640786499874*(alpha[18]*fUpwind_r[46]+alpha[19]*fUpwind_r[45]+alpha[4]*fUpwind_r[44])+13.41640786499874*(alpha[26]*fUpwind_r[42]+alpha[9]*fUpwind_r[41]+alpha[10]*fUpwind_r[40])+13.41640786499874*alpha[17]*fUpwind_r[36]+15.0*(alpha[3]*fUpwind_r[31]+alpha[7]*fUpwind_r[30]+alpha[8]*fUpwind_r[29]+alpha[0]*fUpwind_r[28]+alpha[11]*fUpwind_r[27]+fUpwind_r[11]*alpha[27]+alpha[14]*fUpwind_r[26]+fUpwind_r[14]*alpha[26]+alpha[16]*fUpwind_r[25]+alpha[1]*fUpwind_r[24]+alpha[2]*fUpwind_r[23]+alpha[18]*fUpwind_r[22]+fUpwind_r[18]*alpha[22]+alpha[19]*fUpwind_r[21]+fUpwind_r[19]*alpha[21]+alpha[4]*fUpwind_r[20]+fUpwind_r[4]*alpha[20]+alpha[5]*fUpwind_r[17]+fUpwind_r[5]*alpha[17]+alpha[6]*fUpwind_r[15]+alpha[9]*fUpwind_r[13]+fUpwind_r[9]*alpha[13]+alpha[10]*fUpwind_r[12]+fUpwind_r[10]*alpha[12])); 
  Ghat_r[29] = 0.01178511301977579*(13.41640786499874*alpha[13]*fUpwind_r[63]+13.41640786499874*(alpha[20]*fUpwind_r[62]+alpha[5]*fUpwind_r[61]+alpha[22]*fUpwind_r[60])+13.41640786499874*(alpha[12]*fUpwind_r[58]+alpha[27]*fUpwind_r[57]+alpha[14]*fUpwind_r[56])+13.41640786499874*alpha[21]*fUpwind_r[52]+13.41640786499874*alpha[10]*fUpwind_r[47]+13.41640786499874*(alpha[17]*fUpwind_r[46]+alpha[4]*fUpwind_r[45]+alpha[19]*fUpwind_r[44])+13.41640786499874*(alpha[9]*fUpwind_r[42]+alpha[26]*fUpwind_r[41]+alpha[11]*fUpwind_r[40])+13.41640786499874*alpha[18]*fUpwind_r[36]+15.0*(alpha[2]*fUpwind_r[31]+alpha[6]*fUpwind_r[30]+alpha[0]*fUpwind_r[29]+alpha[8]*fUpwind_r[28]+alpha[10]*fUpwind_r[27]+fUpwind_r[10]*alpha[27]+alpha[13]*fUpwind_r[26]+fUpwind_r[13]*alpha[26]+alpha[1]*fUpwind_r[25]+alpha[16]*fUpwind_r[24]+alpha[3]*fUpwind_r[23]+alpha[17]*fUpwind_r[22]+fUpwind_r[17]*alpha[22]+alpha[4]*fUpwind_r[21]+fUpwind_r[4]*alpha[21]+alpha[19]*fUpwind_r[20]+fUpwind_r[19]*alpha[20]+alpha[5]*fUpwind_r[18]+fUpwind_r[5]*alpha[18]+alpha[7]*fUpwind_r[15]+alpha[9]*fUpwind_r[14]+fUpwind_r[9]*alpha[14]+alpha[11]*fUpwind_r[12]+fUpwind_r[11]*alpha[12])); 
  Ghat_r[30] = 0.01178511301977579*(13.41640786499874*alpha[12]*fUpwind_r[63]+13.41640786499874*(alpha[5]*fUpwind_r[62]+alpha[20]*fUpwind_r[61]+alpha[21]*fUpwind_r[60])+13.41640786499874*(alpha[13]*fUpwind_r[58]+alpha[14]*fUpwind_r[57]+alpha[27]*fUpwind_r[56])+13.41640786499874*alpha[22]*fUpwind_r[52]+13.41640786499874*alpha[9]*fUpwind_r[47]+13.41640786499874*(alpha[4]*fUpwind_r[46]+alpha[17]*fUpwind_r[45]+alpha[18]*fUpwind_r[44])+13.41640786499874*(alpha[10]*fUpwind_r[42]+alpha[11]*fUpwind_r[41]+alpha[26]*fUpwind_r[40])+13.41640786499874*alpha[19]*fUpwind_r[36]+15.0*(alpha[1]*fUpwind_r[31]+alpha[0]*fUpwind_r[30]+alpha[6]*fUpwind_r[29]+alpha[7]*fUpwind_r[28]+alpha[9]*fUpwind_r[27]+fUpwind_r[9]*alpha[27]+alpha[12]*fUpwind_r[26]+fUpwind_r[12]*alpha[26]+alpha[2]*fUpwind_r[25]+alpha[3]*fUpwind_r[24]+alpha[16]*fUpwind_r[23]+alpha[4]*fUpwind_r[22]+fUpwind_r[4]*alpha[22]+alpha[17]*fUpwind_r[21]+fUpwind_r[17]*alpha[21]+alpha[18]*fUpwind_r[20]+fUpwind_r[18]*alpha[20]+alpha[5]*fUpwind_r[19]+fUpwind_r[5]*alpha[19]+alpha[8]*fUpwind_r[15]+alpha[10]*fUpwind_r[14]+fUpwind_r[10]*alpha[14]+alpha[11]*fUpwind_r[13]+fUpwind_r[11]*alpha[13])); 
  Ghat_r[31] = 0.01178511301977579*(13.41640786499874*alpha[5]*fUpwind_r[63]+13.41640786499874*(alpha[12]*fUpwind_r[62]+alpha[13]*fUpwind_r[61]+alpha[14]*fUpwind_r[60])+13.41640786499874*(alpha[20]*fUpwind_r[58]+alpha[21]*fUpwind_r[57]+alpha[22]*fUpwind_r[56])+13.41640786499874*alpha[27]*fUpwind_r[52]+13.41640786499874*alpha[4]*fUpwind_r[47]+13.41640786499874*(alpha[9]*fUpwind_r[46]+alpha[10]*fUpwind_r[45]+alpha[11]*fUpwind_r[44])+13.41640786499874*(alpha[17]*fUpwind_r[42]+alpha[18]*fUpwind_r[41]+alpha[19]*fUpwind_r[40])+13.41640786499874*alpha[26]*fUpwind_r[36]+15.0*(alpha[0]*fUpwind_r[31]+alpha[1]*fUpwind_r[30]+alpha[2]*fUpwind_r[29]+alpha[3]*fUpwind_r[28]+alpha[4]*fUpwind_r[27]+fUpwind_r[4]*alpha[27]+alpha[5]*fUpwind_r[26]+fUpwind_r[5]*alpha[26]+alpha[6]*fUpwind_r[25]+alpha[7]*fUpwind_r[24]+alpha[8]*fUpwind_r[23]+alpha[9]*fUpwind_r[22]+fUpwind_r[9]*alpha[22]+alpha[10]*fUpwind_r[21]+fUpwind_r[10]*alpha[21]+alpha[11]*fUpwind_r[20]+fUpwind_r[11]*alpha[20]+alpha[12]*fUpwind_r[19]+fUpwind_r[12]*alpha[19]+alpha[13]*fUpwind_r[18]+fUpwind_r[13]*alpha[18]+alpha[14]*fUpwind_r[17]+fUpwind_r[14]*alpha[17]+fUpwind_r[15]*alpha[16])); 
  Ghat_r[32] = 0.01178511301977579*(15.0*alpha[27]*fUpwind_r[47]+15.0*(alpha[22]*fUpwind_r[46]+alpha[21]*fUpwind_r[45]+alpha[20]*fUpwind_r[44]+alpha[16]*fUpwind_r[43])+15.0*(alpha[14]*fUpwind_r[42]+alpha[13]*fUpwind_r[41]+alpha[12]*fUpwind_r[40]+alpha[8]*fUpwind_r[39]+alpha[7]*fUpwind_r[38]+alpha[6]*fUpwind_r[37])+15.0*(alpha[5]*fUpwind_r[36]+alpha[3]*fUpwind_r[35]+alpha[2]*fUpwind_r[34]+alpha[1]*fUpwind_r[33])+15.0*alpha[0]*fUpwind_r[32]+13.41640786499874*(alpha[26]*fUpwind_r[26]+alpha[19]*fUpwind_r[19]+alpha[18]*fUpwind_r[18]+alpha[17]*fUpwind_r[17]+alpha[11]*fUpwind_r[11]+alpha[10]*fUpwind_r[10]+alpha[9]*fUpwind_r[9]+alpha[4]*fUpwind_r[4])); 
  Ghat_r[33] = 0.01178511301977579*(15.0*alpha[22]*fUpwind_r[47]+15.0*(alpha[27]*fUpwind_r[46]+alpha[14]*fUpwind_r[45]+alpha[13]*fUpwind_r[44]+alpha[8]*fUpwind_r[43])+15.0*(alpha[21]*fUpwind_r[42]+alpha[20]*fUpwind_r[41]+alpha[5]*fUpwind_r[40]+alpha[16]*fUpwind_r[39]+alpha[3]*fUpwind_r[38]+alpha[2]*fUpwind_r[37])+15.0*(alpha[12]*fUpwind_r[36]+alpha[7]*fUpwind_r[35]+alpha[6]*fUpwind_r[34]+alpha[0]*fUpwind_r[33])+15.0*alpha[1]*fUpwind_r[32]+13.41640786499874*(alpha[19]*fUpwind_r[26]+fUpwind_r[19]*alpha[26]+alpha[11]*fUpwind_r[18]+fUpwind_r[11]*alpha[18]+alpha[10]*fUpwind_r[17]+fUpwind_r[10]*alpha[17]+alpha[4]*fUpwind_r[9]+fUpwind_r[4]*alpha[9])); 
  Ghat_r[34] = 0.01178511301977579*(15.0*alpha[21]*fUpwind_r[47]+15.0*(alpha[14]*fUpwind_r[46]+alpha[27]*fUpwind_r[45]+alpha[12]*fUpwind_r[44]+alpha[7]*fUpwind_r[43])+15.0*(alpha[22]*fUpwind_r[42]+alpha[5]*fUpwind_r[41]+alpha[20]*fUpwind_r[40]+alpha[3]*fUpwind_r[39]+alpha[16]*fUpwind_r[38]+alpha[1]*fUpwind_r[37])+15.0*(alpha[13]*fUpwind_r[36]+alpha[8]*fUpwind_r[35]+alpha[0]*fUpwind_r[34]+alpha[6]*fUpwind_r[33])+15.0*alpha[2]*fUpwind_r[32]+13.41640786499874*(alpha[18]*fUpwind_r[26]+fUpwind_r[18]*alpha[26]+alpha[11]*fUpwind_r[19]+fUpwind_r[11]*alpha[19]+alpha[9]*fUpwind_r[17]+fUpwind_r[9]*alpha[17]+alpha[4]*fUpwind_r[10]+fUpwind_r[4]*alpha[10])); 
  Ghat_r[35] = 0.01178511301977579*(15.0*alpha[20]*fUpwind_r[47]+15.0*(alpha[13]*fUpwind_r[46]+alpha[12]*fUpwind_r[45]+alpha[27]*fUpwind_r[44]+alpha[6]*fUpwind_r[43])+15.0*(alpha[5]*fUpwind_r[42]+alpha[22]*fUpwind_r[41]+alpha[21]*fUpwind_r[40]+alpha[2]*fUpwind_r[39]+alpha[1]*fUpwind_r[38]+alpha[16]*fUpwind_r[37])+15.0*(alpha[14]*fUpwind_r[36]+alpha[0]*fUpwind_r[35]+alpha[8]*fUpwind_r[34]+alpha[7]*fUpwind_r[33])+15.0*alpha[3]*fUpwind_r[32]+13.41640786499874*(alpha[17]*fUpwind_r[26]+fUpwind_r[17]*alpha[26]+alpha[10]*fUpwind_r[19]+fUpwind_r[10]*alpha[19]+alpha[9]*fUpwind_r[18]+fUpwind_r[9]*alpha[18]+alpha[4]*fUpwind_r[11]+fUpwind_r[4]*alpha[11])); 
  Ghat_r[36] = 0.01178511301977579*(15.0*alpha[16]*fUpwind_r[47]+15.0*(alpha[8]*fUpwind_r[46]+alpha[7]*fUpwind_r[45]+alpha[6]*fUpwind_r[44]+alpha[27]*fUpwind_r[43])+15.0*(alpha[3]*fUpwind_r[42]+alpha[2]*fUpwind_r[41]+alpha[1]*fUpwind_r[40]+alpha[22]*fUpwind_r[39]+alpha[21]*fUpwind_r[38]+alpha[20]*fUpwind_r[37])+15.0*(alpha[0]*fUpwind_r[36]+alpha[14]*fUpwind_r[35]+alpha[13]*fUpwind_r[34]+alpha[12]*fUpwind_r[33])+15.0*alpha[5]*fUpwind_r[32]+13.41640786499874*(alpha[26]*fUpwind_r[31]+alpha[19]*fUpwind_r[30]+alpha[18]*fUpwind_r[29]+alpha[17]*fUpwind_r[28]+alpha[11]*fUpwind_r[25]+alpha[10]*fUpwind_r[24]+alpha[9]*fUpwind_r[23]+alpha[4]*fUpwind_r[15])); 
  Ghat_r[37] = 0.01178511301977579*(15.0*alpha[14]*fUpwind_r[47]+15.0*(alpha[21]*fUpwind_r[46]+alpha[22]*fUpwind_r[45]+alpha[5]*fUpwind_r[44]+alpha[3]*fUpwind_r[43])+15.0*(alpha[27]*fUpwind_r[42]+alpha[12]*fUpwind_r[41]+alpha[13]*fUpwind_r[40]+alpha[7]*fUpwind_r[39]+alpha[8]*fUpwind_r[38]+alpha[0]*fUpwind_r[37])+15.0*(alpha[20]*fUpwind_r[36]+alpha[16]*fUpwind_r[35]+alpha[1]*fUpwind_r[34]+alpha[2]*fUpwind_r[33])+15.0*alpha[6]*fUpwind_r[32]+13.41640786499874*(alpha[11]*fUpwind_r[26]+fUpwind_r[11]*alpha[26]+alpha[18]*fUpwind_r[19]+fUpwind_r[18]*alpha[19]+alpha[4]*fUpwind_r[17]+fUpwind_r[4]*alpha[17]+alpha[9]*fUpwind_r[10]+fUpwind_r[9]*alpha[10])); 
  Ghat_r[38] = 0.01178511301977579*(15.0*alpha[13]*fUpwind_r[47]+15.0*(alpha[20]*fUpwind_r[46]+alpha[5]*fUpwind_r[45]+alpha[22]*fUpwind_r[44]+alpha[2]*fUpwind_r[43])+15.0*(alpha[12]*fUpwind_r[42]+alpha[27]*fUpwind_r[41]+alpha[14]*fUpwind_r[40]+alpha[6]*fUpwind_r[39]+alpha[0]*fUpwind_r[38]+alpha[8]*fUpwind_r[37])+15.0*(alpha[21]*fUpwind_r[36]+alpha[1]*fUpwind_r[35]+alpha[16]*fUpwind_r[34]+alpha[3]*fUpwind_r[33])+15.0*alpha[7]*fUpwind_r[32]+13.41640786499874*(alpha[10]*fUpwind_r[26]+fUpwind_r[10]*alpha[26]+alpha[17]*fUpwind_r[19]+fUpwind_r[17]*alpha[19]+alpha[4]*fUpwind_r[18]+fUpwind_r[4]*alpha[18]+alpha[9]*fUpwind_r[11]+fUpwind_r[9]*alpha[11])); 
  Ghat_r[39] = 0.01178511301977579*(15.0*alpha[12]*fUpwind_r[47]+15.0*(alpha[5]*fUpwind_r[46]+alpha[20]*fUpwind_r[45]+alpha[21]*fUpwind_r[44]+alpha[1]*fUpwind_r[43])+15.0*(alpha[13]*fUpwind_r[42]+alpha[14]*fUpwind_r[41]+alpha[27]*fUpwind_r[40]+alpha[0]*fUpwind_r[39]+alpha[6]*fUpwind_r[38]+alpha[7]*fUpwind_r[37])+15.0*(alpha[22]*fUpwind_r[36]+alpha[2]*fUpwind_r[35]+alpha[3]*fUpwind_r[34]+alpha[16]*fUpwind_r[33])+15.0*alpha[8]*fUpwind_r[32]+13.41640786499874*(alpha[9]*fUpwind_r[26]+fUpwind_r[9]*alpha[26]+alpha[4]*fUpwind_r[19]+fUpwind_r[4]*alpha[19]+alpha[17]*fUpwind_r[18]+fUpwind_r[17]*alpha[18]+alpha[10]*fUpwind_r[11]+fUpwind_r[10]*alpha[11])); 
  Ghat_r[40] = 0.01178511301977579*(15.0*alpha[8]*fUpwind_r[47]+15.0*(alpha[16]*fUpwind_r[46]+alpha[3]*fUpwind_r[45]+alpha[2]*fUpwind_r[44]+alpha[22]*fUpwind_r[43])+15.0*(alpha[7]*fUpwind_r[42]+alpha[6]*fUpwind_r[41]+alpha[0]*fUpwind_r[40]+alpha[27]*fUpwind_r[39]+alpha[14]*fUpwind_r[38]+alpha[13]*fUpwind_r[37])+15.0*(alpha[1]*fUpwind_r[36]+alpha[21]*fUpwind_r[35]+alpha[20]*fUpwind_r[34]+alpha[5]*fUpwind_r[33])+15.0*alpha[12]*fUpwind_r[32]+13.41640786499874*(alpha[19]*fUpwind_r[31]+alpha[26]*fUpwind_r[30]+alpha[11]*fUpwind_r[29]+alpha[10]*fUpwind_r[28]+alpha[18]*fUpwind_r[25]+alpha[17]*fUpwind_r[24]+alpha[4]*fUpwind_r[23]+alpha[9]*fUpwind_r[15])); 
  Ghat_r[41] = 0.01178511301977579*(15.0*alpha[7]*fUpwind_r[47]+15.0*(alpha[3]*fUpwind_r[46]+alpha[16]*fUpwind_r[45]+alpha[1]*fUpwind_r[44]+alpha[21]*fUpwind_r[43])+15.0*(alpha[8]*fUpwind_r[42]+alpha[0]*fUpwind_r[41]+alpha[6]*fUpwind_r[40]+alpha[14]*fUpwind_r[39]+alpha[27]*fUpwind_r[38]+alpha[12]*fUpwind_r[37])+15.0*(alpha[2]*fUpwind_r[36]+alpha[22]*fUpwind_r[35]+alpha[5]*fUpwind_r[34]+alpha[20]*fUpwind_r[33])+15.0*alpha[13]*fUpwind_r[32]+13.41640786499874*(alpha[18]*fUpwind_r[31]+alpha[11]*fUpwind_r[30]+alpha[26]*fUpwind_r[29]+alpha[9]*fUpwind_r[28]+alpha[19]*fUpwind_r[25]+alpha[4]*fUpwind_r[24]+alpha[17]*fUpwind_r[23]+alpha[10]*fUpwind_r[15])); 
  Ghat_r[42] = 0.01178511301977579*(15.0*alpha[6]*fUpwind_r[47]+15.0*(alpha[2]*fUpwind_r[46]+alpha[1]*fUpwind_r[45]+alpha[16]*fUpwind_r[44]+alpha[20]*fUpwind_r[43])+15.0*(alpha[0]*fUpwind_r[42]+alpha[8]*fUpwind_r[41]+alpha[7]*fUpwind_r[40]+alpha[13]*fUpwind_r[39]+alpha[12]*fUpwind_r[38]+alpha[27]*fUpwind_r[37])+15.0*(alpha[3]*fUpwind_r[36]+alpha[5]*fUpwind_r[35]+alpha[22]*fUpwind_r[34]+alpha[21]*fUpwind_r[33])+15.0*alpha[14]*fUpwind_r[32]+13.41640786499874*(alpha[17]*fUpwind_r[31]+alpha[10]*fUpwind_r[30]+alpha[9]*fUpwind_r[29]+alpha[26]*fUpwind_r[28]+alpha[4]*fUpwind_r[25]+alpha[19]*fUpwind_r[24]+alpha[18]*fUpwind_r[23]+alpha[11]*fUpwind_r[15])); 
  Ghat_r[43] = 0.01178511301977579*(15.0*alpha[5]*fUpwind_r[47]+15.0*(alpha[12]*fUpwind_r[46]+alpha[13]*fUpwind_r[45]+alpha[14]*fUpwind_r[44]+alpha[0]*fUpwind_r[43])+15.0*(alpha[20]*fUpwind_r[42]+alpha[21]*fUpwind_r[41]+alpha[22]*fUpwind_r[40]+alpha[1]*fUpwind_r[39]+alpha[2]*fUpwind_r[38]+alpha[3]*fUpwind_r[37])+15.0*(alpha[27]*fUpwind_r[36]+alpha[6]*fUpwind_r[35]+alpha[7]*fUpwind_r[34]+alpha[8]*fUpwind_r[33])+15.0*alpha[16]*fUpwind_r[32]+13.41640786499874*(alpha[4]*fUpwind_r[26]+fUpwind_r[4]*alpha[26]+alpha[9]*fUpwind_r[19]+fUpwind_r[9]*alpha[19]+alpha[10]*fUpwind_r[18]+fUpwind_r[10]*alpha[18]+alpha[11]*fUpwind_r[17]+fUpwind_r[11]*alpha[17])); 
  Ghat_r[44] = 0.01178511301977579*(15.0*alpha[3]*fUpwind_r[47]+15.0*(alpha[7]*fUpwind_r[46]+alpha[8]*fUpwind_r[45]+alpha[0]*fUpwind_r[44]+alpha[14]*fUpwind_r[43])+15.0*(alpha[16]*fUpwind_r[42]+alpha[1]*fUpwind_r[41]+alpha[2]*fUpwind_r[40]+alpha[21]*fUpwind_r[39]+alpha[22]*fUpwind_r[38]+alpha[5]*fUpwind_r[37])+15.0*(alpha[6]*fUpwind_r[36]+alpha[27]*fUpwind_r[35]+alpha[12]*fUpwind_r[34]+alpha[13]*fUpwind_r[33])+15.0*alpha[20]*fUpwind_r[32]+13.41640786499874*(alpha[11]*fUpwind_r[31]+alpha[18]*fUpwind_r[30]+alpha[19]*fUpwind_r[29]+alpha[4]*fUpwind_r[28]+fUpwind_r[25]*alpha[26]+alpha[9]*fUpwind_r[24]+alpha[10]*fUpwind_r[23]+fUpwind_r[15]*alpha[17])); 
  Ghat_r[45] = 0.01178511301977579*(15.0*alpha[2]*fUpwind_r[47]+15.0*(alpha[6]*fUpwind_r[46]+alpha[0]*fUpwind_r[45]+alpha[8]*fUpwind_r[44]+alpha[13]*fUpwind_r[43])+15.0*(alpha[1]*fUpwind_r[42]+alpha[16]*fUpwind_r[41]+alpha[3]*fUpwind_r[40]+alpha[20]*fUpwind_r[39]+alpha[5]*fUpwind_r[38]+alpha[22]*fUpwind_r[37])+15.0*(alpha[7]*fUpwind_r[36]+alpha[12]*fUpwind_r[35]+alpha[27]*fUpwind_r[34]+alpha[14]*fUpwind_r[33])+15.0*alpha[21]*fUpwind_r[32]+13.41640786499874*(alpha[10]*fUpwind_r[31]+alpha[17]*fUpwind_r[30]+alpha[4]*fUpwind_r[29]+alpha[19]*fUpwind_r[28]+fUpwind_r[24]*alpha[26]+alpha[9]*fUpwind_r[25]+alpha[11]*fUpwind_r[23]+fUpwind_r[15]*alpha[18])); 
  Ghat_r[46] = 0.01178511301977579*(15.0*alpha[1]*fUpwind_r[47]+15.0*(alpha[0]*fUpwind_r[46]+alpha[6]*fUpwind_r[45]+alpha[7]*fUpwind_r[44]+alpha[12]*fUpwind_r[43])+15.0*(alpha[2]*fUpwind_r[42]+alpha[3]*fUpwind_r[41]+alpha[16]*fUpwind_r[40]+alpha[5]*fUpwind_r[39]+alpha[20]*fUpwind_r[38]+alpha[21]*fUpwind_r[37])+15.0*(alpha[8]*fUpwind_r[36]+alpha[13]*fUpwind_r[35]+alpha[14]*fUpwind_r[34]+alpha[27]*fUpwind_r[33])+15.0*alpha[22]*fUpwind_r[32]+13.41640786499874*(alpha[9]*fUpwind_r[31]+alpha[4]*fUpwind_r[30]+alpha[17]*fUpwind_r[29]+alpha[18]*fUpwind_r[28]+fUpwind_r[23]*alpha[26]+alpha[10]*fUpwind_r[25]+alpha[11]*fUpwind_r[24]+fUpwind_r[15]*alpha[19])); 
  Ghat_r[47] = 0.01178511301977579*(15.0*alpha[0]*fUpwind_r[47]+15.0*(alpha[1]*fUpwind_r[46]+alpha[2]*fUpwind_r[45]+alpha[3]*fUpwind_r[44]+alpha[5]*fUpwind_r[43])+15.0*(alpha[6]*fUpwind_r[42]+alpha[7]*fUpwind_r[41]+alpha[8]*fUpwind_r[40]+alpha[12]*fUpwind_r[39]+alpha[13]*fUpwind_r[38]+alpha[14]*fUpwind_r[37])+15.0*(alpha[16]*fUpwind_r[36]+alpha[20]*fUpwind_r[35]+alpha[21]*fUpwind_r[34]+alpha[22]*fUpwind_r[33])+15.0*alpha[27]*fUpwind_r[32]+13.41640786499874*(alpha[4]*fUpwind_r[31]+alpha[9]*fUpwind_r[30]+alpha[10]*fUpwind_r[29]+alpha[11]*fUpwind_r[28]+fUpwind_r[15]*alpha[26]+alpha[17]*fUpwind_r[25]+alpha[18]*fUpwind_r[24]+alpha[19]*fUpwind_r[23])); 
  Ghat_r[48] = 0.01178511301977579*(15.0*alpha[26]*fUpwind_r[63]+15.0*(alpha[19]*fUpwind_r[62]+alpha[18]*fUpwind_r[61]+alpha[17]*fUpwind_r[60]+alpha[16]*fUpwind_r[59])+15.0*(alpha[11]*fUpwind_r[58]+alpha[10]*fUpwind_r[57]+alpha[9]*fUpwind_r[56]+alpha[8]*fUpwind_r[55]+alpha[7]*fUpwind_r[54]+alpha[6]*fUpwind_r[53])+15.0*(alpha[4]*fUpwind_r[52]+alpha[3]*fUpwind_r[51]+alpha[2]*fUpwind_r[50]+alpha[1]*fUpwind_r[49])+15.0*alpha[0]*fUpwind_r[48]+13.41640786499874*(alpha[27]*fUpwind_r[27]+alpha[22]*fUpwind_r[22]+alpha[21]*fUpwind_r[21]+alpha[20]*fUpwind_r[20]+alpha[14]*fUpwind_r[14]+alpha[13]*fUpwind_r[13]+alpha[12]*fUpwind_r[12]+alpha[5]*fUpwind_r[5])); 
  Ghat_r[49] = 0.01178511301977579*(15.0*alpha[19]*fUpwind_r[63]+15.0*(alpha[26]*fUpwind_r[62]+alpha[11]*fUpwind_r[61]+alpha[10]*fUpwind_r[60]+alpha[8]*fUpwind_r[59])+15.0*(alpha[18]*fUpwind_r[58]+alpha[17]*fUpwind_r[57]+alpha[4]*fUpwind_r[56]+alpha[16]*fUpwind_r[55]+alpha[3]*fUpwind_r[54]+alpha[2]*fUpwind_r[53])+15.0*(alpha[9]*fUpwind_r[52]+alpha[7]*fUpwind_r[51]+alpha[6]*fUpwind_r[50]+alpha[0]*fUpwind_r[49])+15.0*alpha[1]*fUpwind_r[48]+13.41640786499874*(alpha[22]*fUpwind_r[27]+fUpwind_r[22]*alpha[27]+alpha[14]*fUpwind_r[21]+fUpwind_r[14]*alpha[21]+alpha[13]*fUpwind_r[20]+fUpwind_r[13]*alpha[20]+alpha[5]*fUpwind_r[12]+fUpwind_r[5]*alpha[12])); 
  Ghat_r[50] = 0.01178511301977579*(15.0*alpha[18]*fUpwind_r[63]+15.0*(alpha[11]*fUpwind_r[62]+alpha[26]*fUpwind_r[61]+alpha[9]*fUpwind_r[60]+alpha[7]*fUpwind_r[59])+15.0*(alpha[19]*fUpwind_r[58]+alpha[4]*fUpwind_r[57]+alpha[17]*fUpwind_r[56]+alpha[3]*fUpwind_r[55]+alpha[16]*fUpwind_r[54]+alpha[1]*fUpwind_r[53])+15.0*(alpha[10]*fUpwind_r[52]+alpha[8]*fUpwind_r[51]+alpha[0]*fUpwind_r[50]+alpha[6]*fUpwind_r[49])+15.0*alpha[2]*fUpwind_r[48]+13.41640786499874*(alpha[21]*fUpwind_r[27]+fUpwind_r[21]*alpha[27]+alpha[14]*fUpwind_r[22]+fUpwind_r[14]*alpha[22]+alpha[12]*fUpwind_r[20]+fUpwind_r[12]*alpha[20]+alpha[5]*fUpwind_r[13]+fUpwind_r[5]*alpha[13])); 
  Ghat_r[51] = 0.01178511301977579*(15.0*alpha[17]*fUpwind_r[63]+15.0*(alpha[10]*fUpwind_r[62]+alpha[9]*fUpwind_r[61]+alpha[26]*fUpwind_r[60]+alpha[6]*fUpwind_r[59])+15.0*(alpha[4]*fUpwind_r[58]+alpha[19]*fUpwind_r[57]+alpha[18]*fUpwind_r[56]+alpha[2]*fUpwind_r[55]+alpha[1]*fUpwind_r[54]+alpha[16]*fUpwind_r[53])+15.0*(alpha[11]*fUpwind_r[52]+alpha[0]*fUpwind_r[51]+alpha[8]*fUpwind_r[50]+alpha[7]*fUpwind_r[49])+15.0*alpha[3]*fUpwind_r[48]+13.41640786499874*(alpha[20]*fUpwind_r[27]+fUpwind_r[20]*alpha[27]+alpha[13]*fUpwind_r[22]+fUpwind_r[13]*alpha[22]+alpha[12]*fUpwind_r[21]+fUpwind_r[12]*alpha[21]+alpha[5]*fUpwind_r[14]+fUpwind_r[5]*alpha[14])); 
  Ghat_r[52] = 0.01178511301977579*(15.0*alpha[16]*fUpwind_r[63]+15.0*(alpha[8]*fUpwind_r[62]+alpha[7]*fUpwind_r[61]+alpha[6]*fUpwind_r[60]+alpha[26]*fUpwind_r[59])+15.0*(alpha[3]*fUpwind_r[58]+alpha[2]*fUpwind_r[57]+alpha[1]*fUpwind_r[56]+alpha[19]*fUpwind_r[55]+alpha[18]*fUpwind_r[54]+alpha[17]*fUpwind_r[53])+15.0*(alpha[0]*fUpwind_r[52]+alpha[11]*fUpwind_r[51]+alpha[10]*fUpwind_r[50]+alpha[9]*fUpwind_r[49])+15.0*alpha[4]*fUpwind_r[48]+13.41640786499874*(alpha[27]*fUpwind_r[31]+alpha[22]*fUpwind_r[30]+alpha[21]*fUpwind_r[29]+alpha[20]*fUpwind_r[28]+alpha[14]*fUpwind_r[25]+alpha[13]*fUpwind_r[24]+alpha[12]*fUpwind_r[23]+alpha[5]*fUpwind_r[15])); 
  Ghat_r[53] = 0.01178511301977579*(15.0*alpha[11]*fUpwind_r[63]+15.0*(alpha[18]*fUpwind_r[62]+alpha[19]*fUpwind_r[61]+alpha[4]*fUpwind_r[60]+alpha[3]*fUpwind_r[59])+15.0*(alpha[26]*fUpwind_r[58]+alpha[9]*fUpwind_r[57]+alpha[10]*fUpwind_r[56]+alpha[7]*fUpwind_r[55]+alpha[8]*fUpwind_r[54]+alpha[0]*fUpwind_r[53])+15.0*(alpha[17]*fUpwind_r[52]+alpha[16]*fUpwind_r[51]+alpha[1]*fUpwind_r[50]+alpha[2]*fUpwind_r[49])+15.0*alpha[6]*fUpwind_r[48]+13.41640786499874*(alpha[14]*fUpwind_r[27]+fUpwind_r[14]*alpha[27]+alpha[21]*fUpwind_r[22]+fUpwind_r[21]*alpha[22]+alpha[5]*fUpwind_r[20]+fUpwind_r[5]*alpha[20]+alpha[12]*fUpwind_r[13]+fUpwind_r[12]*alpha[13])); 
  Ghat_r[54] = 0.01178511301977579*(15.0*alpha[10]*fUpwind_r[63]+15.0*(alpha[17]*fUpwind_r[62]+alpha[4]*fUpwind_r[61]+alpha[19]*fUpwind_r[60]+alpha[2]*fUpwind_r[59])+15.0*(alpha[9]*fUpwind_r[58]+alpha[26]*fUpwind_r[57]+alpha[11]*fUpwind_r[56]+alpha[6]*fUpwind_r[55]+alpha[0]*fUpwind_r[54]+alpha[8]*fUpwind_r[53])+15.0*(alpha[18]*fUpwind_r[52]+alpha[1]*fUpwind_r[51]+alpha[16]*fUpwind_r[50]+alpha[3]*fUpwind_r[49])+15.0*alpha[7]*fUpwind_r[48]+13.41640786499874*(alpha[13]*fUpwind_r[27]+fUpwind_r[13]*alpha[27]+alpha[20]*fUpwind_r[22]+fUpwind_r[20]*alpha[22]+alpha[5]*fUpwind_r[21]+fUpwind_r[5]*alpha[21]+alpha[12]*fUpwind_r[14]+fUpwind_r[12]*alpha[14])); 
  Ghat_r[55] = 0.01178511301977579*(15.0*alpha[9]*fUpwind_r[63]+15.0*(alpha[4]*fUpwind_r[62]+alpha[17]*fUpwind_r[61]+alpha[18]*fUpwind_r[60]+alpha[1]*fUpwind_r[59])+15.0*(alpha[10]*fUpwind_r[58]+alpha[11]*fUpwind_r[57]+alpha[26]*fUpwind_r[56]+alpha[0]*fUpwind_r[55]+alpha[6]*fUpwind_r[54]+alpha[7]*fUpwind_r[53])+15.0*(alpha[19]*fUpwind_r[52]+alpha[2]*fUpwind_r[51]+alpha[3]*fUpwind_r[50]+alpha[16]*fUpwind_r[49])+15.0*alpha[8]*fUpwind_r[48]+13.41640786499874*(alpha[12]*fUpwind_r[27]+fUpwind_r[12]*alpha[27]+alpha[5]*fUpwind_r[22]+fUpwind_r[5]*alpha[22]+alpha[20]*fUpwind_r[21]+fUpwind_r[20]*alpha[21]+alpha[13]*fUpwind_r[14]+fUpwind_r[13]*alpha[14])); 
  Ghat_r[56] = 0.01178511301977579*(15.0*alpha[8]*fUpwind_r[63]+15.0*(alpha[16]*fUpwind_r[62]+alpha[3]*fUpwind_r[61]+alpha[2]*fUpwind_r[60]+alpha[19]*fUpwind_r[59])+15.0*(alpha[7]*fUpwind_r[58]+alpha[6]*fUpwind_r[57]+alpha[0]*fUpwind_r[56]+alpha[26]*fUpwind_r[55]+alpha[11]*fUpwind_r[54]+alpha[10]*fUpwind_r[53])+15.0*(alpha[1]*fUpwind_r[52]+alpha[18]*fUpwind_r[51]+alpha[17]*fUpwind_r[50]+alpha[4]*fUpwind_r[49])+15.0*alpha[9]*fUpwind_r[48]+13.41640786499874*(alpha[22]*fUpwind_r[31]+alpha[27]*fUpwind_r[30]+alpha[14]*fUpwind_r[29]+alpha[13]*fUpwind_r[28]+alpha[21]*fUpwind_r[25]+alpha[20]*fUpwind_r[24]+alpha[5]*fUpwind_r[23]+alpha[12]*fUpwind_r[15])); 
  Ghat_r[57] = 0.01178511301977579*(15.0*alpha[7]*fUpwind_r[63]+15.0*(alpha[3]*fUpwind_r[62]+alpha[16]*fUpwind_r[61]+alpha[1]*fUpwind_r[60]+alpha[18]*fUpwind_r[59])+15.0*(alpha[8]*fUpwind_r[58]+alpha[0]*fUpwind_r[57]+alpha[6]*fUpwind_r[56]+alpha[11]*fUpwind_r[55]+alpha[26]*fUpwind_r[54]+alpha[9]*fUpwind_r[53])+15.0*(alpha[2]*fUpwind_r[52]+alpha[19]*fUpwind_r[51]+alpha[4]*fUpwind_r[50]+alpha[17]*fUpwind_r[49])+15.0*alpha[10]*fUpwind_r[48]+13.41640786499874*(alpha[21]*fUpwind_r[31]+alpha[14]*fUpwind_r[30]+alpha[27]*fUpwind_r[29]+alpha[12]*fUpwind_r[28]+alpha[22]*fUpwind_r[25]+alpha[5]*fUpwind_r[24]+alpha[20]*fUpwind_r[23]+alpha[13]*fUpwind_r[15])); 
  Ghat_r[58] = 0.01178511301977579*(15.0*alpha[6]*fUpwind_r[63]+15.0*(alpha[2]*fUpwind_r[62]+alpha[1]*fUpwind_r[61]+alpha[16]*fUpwind_r[60]+alpha[17]*fUpwind_r[59])+15.0*(alpha[0]*fUpwind_r[58]+alpha[8]*fUpwind_r[57]+alpha[7]*fUpwind_r[56]+alpha[10]*fUpwind_r[55]+alpha[9]*fUpwind_r[54]+alpha[26]*fUpwind_r[53])+15.0*(alpha[3]*fUpwind_r[52]+alpha[4]*fUpwind_r[51]+alpha[19]*fUpwind_r[50]+alpha[18]*fUpwind_r[49])+15.0*alpha[11]*fUpwind_r[48]+13.41640786499874*(alpha[20]*fUpwind_r[31]+alpha[13]*fUpwind_r[30]+alpha[12]*fUpwind_r[29]+alpha[27]*fUpwind_r[28]+alpha[5]*fUpwind_r[25]+alpha[22]*fUpwind_r[24]+alpha[21]*fUpwind_r[23]+alpha[14]*fUpwind_r[15])); 
  Ghat_r[59] = 0.01178511301977579*(15.0*alpha[4]*fUpwind_r[63]+15.0*(alpha[9]*fUpwind_r[62]+alpha[10]*fUpwind_r[61]+alpha[11]*fUpwind_r[60]+alpha[0]*fUpwind_r[59])+15.0*(alpha[17]*fUpwind_r[58]+alpha[18]*fUpwind_r[57]+alpha[19]*fUpwind_r[56]+alpha[1]*fUpwind_r[55]+alpha[2]*fUpwind_r[54]+alpha[3]*fUpwind_r[53])+15.0*(alpha[26]*fUpwind_r[52]+alpha[6]*fUpwind_r[51]+alpha[7]*fUpwind_r[50]+alpha[8]*fUpwind_r[49])+15.0*alpha[16]*fUpwind_r[48]+13.41640786499874*(alpha[5]*fUpwind_r[27]+fUpwind_r[5]*alpha[27]+alpha[12]*fUpwind_r[22]+fUpwind_r[12]*alpha[22]+alpha[13]*fUpwind_r[21]+fUpwind_r[13]*alpha[21]+alpha[14]*fUpwind_r[20]+fUpwind_r[14]*alpha[20])); 
  Ghat_r[60] = 0.01178511301977579*(15.0*alpha[3]*fUpwind_r[63]+15.0*(alpha[7]*fUpwind_r[62]+alpha[8]*fUpwind_r[61]+alpha[0]*fUpwind_r[60]+alpha[11]*fUpwind_r[59])+15.0*(alpha[16]*fUpwind_r[58]+alpha[1]*fUpwind_r[57]+alpha[2]*fUpwind_r[56]+alpha[18]*fUpwind_r[55]+alpha[19]*fUpwind_r[54]+alpha[4]*fUpwind_r[53])+15.0*(alpha[6]*fUpwind_r[52]+alpha[26]*fUpwind_r[51]+alpha[9]*fUpwind_r[50]+alpha[10]*fUpwind_r[49])+15.0*alpha[17]*fUpwind_r[48]+13.41640786499874*(alpha[14]*fUpwind_r[31]+alpha[21]*fUpwind_r[30]+alpha[22]*fUpwind_r[29]+alpha[5]*fUpwind_r[28]+fUpwind_r[25]*alpha[27]+alpha[12]*fUpwind_r[24]+alpha[13]*fUpwind_r[23]+fUpwind_r[15]*alpha[20])); 
  Ghat_r[61] = 0.01178511301977579*(15.0*alpha[2]*fUpwind_r[63]+15.0*(alpha[6]*fUpwind_r[62]+alpha[0]*fUpwind_r[61]+alpha[8]*fUpwind_r[60]+alpha[10]*fUpwind_r[59])+15.0*(alpha[1]*fUpwind_r[58]+alpha[16]*fUpwind_r[57]+alpha[3]*fUpwind_r[56]+alpha[17]*fUpwind_r[55]+alpha[4]*fUpwind_r[54]+alpha[19]*fUpwind_r[53])+15.0*(alpha[7]*fUpwind_r[52]+alpha[9]*fUpwind_r[51]+alpha[26]*fUpwind_r[50]+alpha[11]*fUpwind_r[49])+15.0*alpha[18]*fUpwind_r[48]+13.41640786499874*(alpha[13]*fUpwind_r[31]+alpha[20]*fUpwind_r[30]+alpha[5]*fUpwind_r[29]+alpha[22]*fUpwind_r[28]+fUpwind_r[24]*alpha[27]+alpha[12]*fUpwind_r[25]+alpha[14]*fUpwind_r[23]+fUpwind_r[15]*alpha[21])); 
  Ghat_r[62] = 0.01178511301977579*(15.0*alpha[1]*fUpwind_r[63]+15.0*(alpha[0]*fUpwind_r[62]+alpha[6]*fUpwind_r[61]+alpha[7]*fUpwind_r[60]+alpha[9]*fUpwind_r[59])+15.0*(alpha[2]*fUpwind_r[58]+alpha[3]*fUpwind_r[57]+alpha[16]*fUpwind_r[56]+alpha[4]*fUpwind_r[55]+alpha[17]*fUpwind_r[54]+alpha[18]*fUpwind_r[53])+15.0*(alpha[8]*fUpwind_r[52]+alpha[10]*fUpwind_r[51]+alpha[11]*fUpwind_r[50]+alpha[26]*fUpwind_r[49])+15.0*alpha[19]*fUpwind_r[48]+13.41640786499874*(alpha[12]*fUpwind_r[31]+alpha[5]*fUpwind_r[30]+alpha[20]*fUpwind_r[29]+alpha[21]*fUpwind_r[28]+fUpwind_r[23]*alpha[27]+alpha[13]*fUpwind_r[25]+alpha[14]*fUpwind_r[24]+fUpwind_r[15]*alpha[22])); 
  Ghat_r[63] = 0.01178511301977579*(15.0*alpha[0]*fUpwind_r[63]+15.0*(alpha[1]*fUpwind_r[62]+alpha[2]*fUpwind_r[61]+alpha[3]*fUpwind_r[60]+alpha[4]*fUpwind_r[59])+15.0*(alpha[6]*fUpwind_r[58]+alpha[7]*fUpwind_r[57]+alpha[8]*fUpwind_r[56]+alpha[9]*fUpwind_r[55]+alpha[10]*fUpwind_r[54]+alpha[11]*fUpwind_r[53])+15.0*(alpha[16]*fUpwind_r[52]+alpha[17]*fUpwind_r[51]+alpha[18]*fUpwind_r[50]+alpha[19]*fUpwind_r[49])+15.0*alpha[26]*fUpwind_r[48]+13.41640786499874*(alpha[5]*fUpwind_r[31]+alpha[12]*fUpwind_r[30]+alpha[13]*fUpwind_r[29]+alpha[14]*fUpwind_r[28]+fUpwind_r[15]*alpha[27]+alpha[20]*fUpwind_r[25]+alpha[21]*fUpwind_r[24]+alpha[22]*fUpwind_r[23])); 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv12; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv12; 
  out[2] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv12; 
  out[3] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv12; 
  out[4] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dv12; 
  out[5] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dv12; 
  out[6] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv12; 
  out[7] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dv12; 
  out[8] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dv12; 
  out[9] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*dv12; 
  out[10] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*dv12; 
  out[11] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*dv12; 
  out[12] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*dv12; 
  out[13] += (0.7071067811865475*Ghat_l[12]-0.7071067811865475*Ghat_r[12])*dv12; 
  out[14] += (0.7071067811865475*Ghat_l[13]-0.7071067811865475*Ghat_r[13])*dv12; 
  out[15] += (0.7071067811865475*Ghat_l[14]-0.7071067811865475*Ghat_r[14])*dv12; 
  out[16] += (0.7071067811865475*Ghat_l[15]-0.7071067811865475*Ghat_r[15])*dv12; 
  out[17] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv12; 
  out[18] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv12; 
  out[19] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv12; 
  out[20] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dv12; 
  out[21] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dv12; 
  out[22] += (0.7071067811865475*Ghat_l[16]-0.7071067811865475*Ghat_r[16])*dv12; 
  out[23] += (0.7071067811865475*Ghat_l[17]-0.7071067811865475*Ghat_r[17])*dv12; 
  out[24] += (0.7071067811865475*Ghat_l[18]-0.7071067811865475*Ghat_r[18])*dv12; 
  out[25] += (0.7071067811865475*Ghat_l[19]-0.7071067811865475*Ghat_r[19])*dv12; 
  out[26] += (0.7071067811865475*Ghat_l[20]-0.7071067811865475*Ghat_r[20])*dv12; 
  out[27] += (0.7071067811865475*Ghat_l[21]-0.7071067811865475*Ghat_r[21])*dv12; 
  out[28] += (0.7071067811865475*Ghat_l[22]-0.7071067811865475*Ghat_r[22])*dv12; 
  out[29] += (0.7071067811865475*Ghat_l[23]-0.7071067811865475*Ghat_r[23])*dv12; 
  out[30] += (0.7071067811865475*Ghat_l[24]-0.7071067811865475*Ghat_r[24])*dv12; 
  out[31] += (0.7071067811865475*Ghat_l[25]-0.7071067811865475*Ghat_r[25])*dv12; 
  out[32] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dv12; 
  out[33] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dv12; 
  out[34] += -1.224744871391589*(Ghat_r[8]+Ghat_l[8])*dv12; 
  out[35] += -1.224744871391589*(Ghat_r[9]+Ghat_l[9])*dv12; 
  out[36] += -1.224744871391589*(Ghat_r[10]+Ghat_l[10])*dv12; 
  out[37] += -1.224744871391589*(Ghat_r[11]+Ghat_l[11])*dv12; 
  out[38] += -1.224744871391589*(Ghat_r[12]+Ghat_l[12])*dv12; 
  out[39] += -1.224744871391589*(Ghat_r[13]+Ghat_l[13])*dv12; 
  out[40] += -1.224744871391589*(Ghat_r[14]+Ghat_l[14])*dv12; 
  out[41] += -1.224744871391589*(Ghat_r[15]+Ghat_l[15])*dv12; 
  out[42] += (0.7071067811865475*Ghat_l[26]-0.7071067811865475*Ghat_r[26])*dv12; 
  out[43] += (0.7071067811865475*Ghat_l[27]-0.7071067811865475*Ghat_r[27])*dv12; 
  out[44] += (0.7071067811865475*Ghat_l[28]-0.7071067811865475*Ghat_r[28])*dv12; 
  out[45] += (0.7071067811865475*Ghat_l[29]-0.7071067811865475*Ghat_r[29])*dv12; 
  out[46] += (0.7071067811865475*Ghat_l[30]-0.7071067811865475*Ghat_r[30])*dv12; 
  out[47] += -1.224744871391589*(Ghat_r[16]+Ghat_l[16])*dv12; 
  out[48] += -1.224744871391589*(Ghat_r[17]+Ghat_l[17])*dv12; 
  out[49] += -1.224744871391589*(Ghat_r[18]+Ghat_l[18])*dv12; 
  out[50] += -1.224744871391589*(Ghat_r[19]+Ghat_l[19])*dv12; 
  out[51] += -1.224744871391589*(Ghat_r[20]+Ghat_l[20])*dv12; 
  out[52] += -1.224744871391589*(Ghat_r[21]+Ghat_l[21])*dv12; 
  out[53] += -1.224744871391589*(Ghat_r[22]+Ghat_l[22])*dv12; 
  out[54] += -1.224744871391589*(Ghat_r[23]+Ghat_l[23])*dv12; 
  out[55] += -1.224744871391589*(Ghat_r[24]+Ghat_l[24])*dv12; 
  out[56] += -1.224744871391589*(Ghat_r[25]+Ghat_l[25])*dv12; 
  out[57] += (0.7071067811865475*Ghat_l[31]-0.7071067811865475*Ghat_r[31])*dv12; 
  out[58] += -1.224744871391589*(Ghat_r[26]+Ghat_l[26])*dv12; 
  out[59] += -1.224744871391589*(Ghat_r[27]+Ghat_l[27])*dv12; 
  out[60] += -1.224744871391589*(Ghat_r[28]+Ghat_l[28])*dv12; 
  out[61] += -1.224744871391589*(Ghat_r[29]+Ghat_l[29])*dv12; 
  out[62] += -1.224744871391589*(Ghat_r[30]+Ghat_l[30])*dv12; 
  out[63] += -1.224744871391589*(Ghat_r[31]+Ghat_l[31])*dv12; 
  out[64] += (0.7071067811865475*Ghat_l[32]-0.7071067811865475*Ghat_r[32])*dv12; 
  out[65] += (0.7071067811865475*Ghat_l[33]-0.7071067811865475*Ghat_r[33])*dv12; 
  out[66] += (0.7071067811865475*Ghat_l[34]-0.7071067811865475*Ghat_r[34])*dv12; 
  out[67] += (0.7071067811865475*Ghat_l[35]-0.7071067811865475*Ghat_r[35])*dv12; 
  out[68] += (0.7071067811865475*Ghat_l[36]-0.7071067811865475*Ghat_r[36])*dv12; 
  out[69] += -1.224744871391589*(Ghat_r[32]+Ghat_l[32])*dv12; 
  out[70] += (0.7071067811865475*Ghat_l[37]-0.7071067811865475*Ghat_r[37])*dv12; 
  out[71] += (0.7071067811865475*Ghat_l[38]-0.7071067811865475*Ghat_r[38])*dv12; 
  out[72] += (0.7071067811865475*Ghat_l[39]-0.7071067811865475*Ghat_r[39])*dv12; 
  out[73] += (0.7071067811865475*Ghat_l[40]-0.7071067811865475*Ghat_r[40])*dv12; 
  out[74] += (0.7071067811865475*Ghat_l[41]-0.7071067811865475*Ghat_r[41])*dv12; 
  out[75] += (0.7071067811865475*Ghat_l[42]-0.7071067811865475*Ghat_r[42])*dv12; 
  out[76] += -1.224744871391589*(Ghat_r[33]+Ghat_l[33])*dv12; 
  out[77] += -1.224744871391589*(Ghat_r[34]+Ghat_l[34])*dv12; 
  out[78] += -1.224744871391589*(Ghat_r[35]+Ghat_l[35])*dv12; 
  out[79] += -1.224744871391589*(Ghat_r[36]+Ghat_l[36])*dv12; 
  out[80] += (0.7071067811865475*Ghat_l[43]-0.7071067811865475*Ghat_r[43])*dv12; 
  out[81] += (0.7071067811865475*Ghat_l[44]-0.7071067811865475*Ghat_r[44])*dv12; 
  out[82] += (0.7071067811865475*Ghat_l[45]-0.7071067811865475*Ghat_r[45])*dv12; 
  out[83] += (0.7071067811865475*Ghat_l[46]-0.7071067811865475*Ghat_r[46])*dv12; 
  out[84] += -1.224744871391589*(Ghat_r[37]+Ghat_l[37])*dv12; 
  out[85] += -1.224744871391589*(Ghat_r[38]+Ghat_l[38])*dv12; 
  out[86] += -1.224744871391589*(Ghat_r[39]+Ghat_l[39])*dv12; 
  out[87] += -1.224744871391589*(Ghat_r[40]+Ghat_l[40])*dv12; 
  out[88] += -1.224744871391589*(Ghat_r[41]+Ghat_l[41])*dv12; 
  out[89] += -1.224744871391589*(Ghat_r[42]+Ghat_l[42])*dv12; 
  out[90] += (0.7071067811865475*Ghat_l[47]-0.7071067811865475*Ghat_r[47])*dv12; 
  out[91] += -1.224744871391589*(Ghat_r[43]+Ghat_l[43])*dv12; 
  out[92] += -1.224744871391589*(Ghat_r[44]+Ghat_l[44])*dv12; 
  out[93] += -1.224744871391589*(Ghat_r[45]+Ghat_l[45])*dv12; 
  out[94] += -1.224744871391589*(Ghat_r[46]+Ghat_l[46])*dv12; 
  out[95] += -1.224744871391589*(Ghat_r[47]+Ghat_l[47])*dv12; 
  out[96] += (0.7071067811865475*Ghat_l[48]-0.7071067811865475*Ghat_r[48])*dv12; 
  out[97] += (0.7071067811865475*Ghat_l[49]-0.7071067811865475*Ghat_r[49])*dv12; 
  out[98] += (0.7071067811865475*Ghat_l[50]-0.7071067811865475*Ghat_r[50])*dv12; 
  out[99] += (0.7071067811865475*Ghat_l[51]-0.7071067811865475*Ghat_r[51])*dv12; 
  out[100] += (0.7071067811865475*Ghat_l[52]-0.7071067811865475*Ghat_r[52])*dv12; 
  out[101] += -1.224744871391589*(Ghat_r[48]+Ghat_l[48])*dv12; 
  out[102] += (0.7071067811865475*Ghat_l[53]-0.7071067811865475*Ghat_r[53])*dv12; 
  out[103] += (0.7071067811865475*Ghat_l[54]-0.7071067811865475*Ghat_r[54])*dv12; 
  out[104] += (0.7071067811865475*Ghat_l[55]-0.7071067811865475*Ghat_r[55])*dv12; 
  out[105] += (0.7071067811865475*Ghat_l[56]-0.7071067811865475*Ghat_r[56])*dv12; 
  out[106] += (0.7071067811865475*Ghat_l[57]-0.7071067811865475*Ghat_r[57])*dv12; 
  out[107] += (0.7071067811865475*Ghat_l[58]-0.7071067811865475*Ghat_r[58])*dv12; 
  out[108] += -1.224744871391589*(Ghat_r[49]+Ghat_l[49])*dv12; 
  out[109] += -1.224744871391589*(Ghat_r[50]+Ghat_l[50])*dv12; 
  out[110] += -1.224744871391589*(Ghat_r[51]+Ghat_l[51])*dv12; 
  out[111] += -1.224744871391589*(Ghat_r[52]+Ghat_l[52])*dv12; 
  out[112] += (0.7071067811865475*Ghat_l[59]-0.7071067811865475*Ghat_r[59])*dv12; 
  out[113] += (0.7071067811865475*Ghat_l[60]-0.7071067811865475*Ghat_r[60])*dv12; 
  out[114] += (0.7071067811865475*Ghat_l[61]-0.7071067811865475*Ghat_r[61])*dv12; 
  out[115] += (0.7071067811865475*Ghat_l[62]-0.7071067811865475*Ghat_r[62])*dv12; 
  out[116] += -1.224744871391589*(Ghat_r[53]+Ghat_l[53])*dv12; 
  out[117] += -1.224744871391589*(Ghat_r[54]+Ghat_l[54])*dv12; 
  out[118] += -1.224744871391589*(Ghat_r[55]+Ghat_l[55])*dv12; 
  out[119] += -1.224744871391589*(Ghat_r[56]+Ghat_l[56])*dv12; 
  out[120] += -1.224744871391589*(Ghat_r[57]+Ghat_l[57])*dv12; 
  out[121] += -1.224744871391589*(Ghat_r[58]+Ghat_l[58])*dv12; 
  out[122] += (0.7071067811865475*Ghat_l[63]-0.7071067811865475*Ghat_r[63])*dv12; 
  out[123] += -1.224744871391589*(Ghat_r[59]+Ghat_l[59])*dv12; 
  out[124] += -1.224744871391589*(Ghat_r[60]+Ghat_l[60])*dv12; 
  out[125] += -1.224744871391589*(Ghat_r[61]+Ghat_l[61])*dv12; 
  out[126] += -1.224744871391589*(Ghat_r[62]+Ghat_l[62])*dv12; 
  out[127] += -1.224744871391589*(Ghat_r[63]+Ghat_l[63])*dv12; 
  out[128] += (1.58113883008419*Ghat_l[0]-1.58113883008419*Ghat_r[0])*dv12; 
  out[129] += (1.58113883008419*Ghat_l[1]-1.58113883008419*Ghat_r[1])*dv12; 
  out[130] += (1.58113883008419*Ghat_l[2]-1.58113883008419*Ghat_r[2])*dv12; 
  out[131] += (1.58113883008419*Ghat_l[3]-1.58113883008419*Ghat_r[3])*dv12; 
  out[132] += (1.58113883008419*Ghat_l[4]-1.58113883008419*Ghat_r[4])*dv12; 
  out[133] += (1.58113883008419*Ghat_l[5]-1.58113883008419*Ghat_r[5])*dv12; 
  out[134] += (1.58113883008419*Ghat_l[6]-1.58113883008419*Ghat_r[6])*dv12; 
  out[135] += (1.58113883008419*Ghat_l[7]-1.58113883008419*Ghat_r[7])*dv12; 
  out[136] += (1.58113883008419*Ghat_l[8]-1.58113883008419*Ghat_r[8])*dv12; 
  out[137] += (1.58113883008419*Ghat_l[9]-1.58113883008419*Ghat_r[9])*dv12; 
  out[138] += (1.58113883008419*Ghat_l[10]-1.58113883008419*Ghat_r[10])*dv12; 
  out[139] += (1.58113883008419*Ghat_l[11]-1.58113883008419*Ghat_r[11])*dv12; 
  out[140] += (1.58113883008419*Ghat_l[12]-1.58113883008419*Ghat_r[12])*dv12; 
  out[141] += (1.58113883008419*Ghat_l[13]-1.58113883008419*Ghat_r[13])*dv12; 
  out[142] += (1.58113883008419*Ghat_l[14]-1.58113883008419*Ghat_r[14])*dv12; 
  out[143] += (1.58113883008419*Ghat_l[15]-1.58113883008419*Ghat_r[15])*dv12; 
  out[144] += (1.58113883008419*Ghat_l[16]-1.58113883008419*Ghat_r[16])*dv12; 
  out[145] += (1.58113883008419*Ghat_l[17]-1.58113883008419*Ghat_r[17])*dv12; 
  out[146] += (1.58113883008419*Ghat_l[18]-1.58113883008419*Ghat_r[18])*dv12; 
  out[147] += (1.58113883008419*Ghat_l[19]-1.58113883008419*Ghat_r[19])*dv12; 
  out[148] += (1.58113883008419*Ghat_l[20]-1.58113883008419*Ghat_r[20])*dv12; 
  out[149] += (1.58113883008419*Ghat_l[21]-1.58113883008419*Ghat_r[21])*dv12; 
  out[150] += (1.58113883008419*Ghat_l[22]-1.58113883008419*Ghat_r[22])*dv12; 
  out[151] += (1.58113883008419*Ghat_l[23]-1.58113883008419*Ghat_r[23])*dv12; 
  out[152] += (1.58113883008419*Ghat_l[24]-1.58113883008419*Ghat_r[24])*dv12; 
  out[153] += (1.58113883008419*Ghat_l[25]-1.58113883008419*Ghat_r[25])*dv12; 
  out[154] += (1.58113883008419*Ghat_l[26]-1.58113883008419*Ghat_r[26])*dv12; 
  out[155] += (1.58113883008419*Ghat_l[27]-1.58113883008419*Ghat_r[27])*dv12; 
  out[156] += (1.58113883008419*Ghat_l[28]-1.58113883008419*Ghat_r[28])*dv12; 
  out[157] += (1.58113883008419*Ghat_l[29]-1.58113883008419*Ghat_r[29])*dv12; 
  out[158] += (1.58113883008419*Ghat_l[30]-1.58113883008419*Ghat_r[30])*dv12; 
  out[159] += (1.58113883008419*Ghat_l[31]-1.58113883008419*Ghat_r[31])*dv12; 

} 
