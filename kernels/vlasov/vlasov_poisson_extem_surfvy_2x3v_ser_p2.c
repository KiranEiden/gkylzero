#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_5x_p2_surfx4_eval_quad.h> 
#include <gkyl_basis_ser_5x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_poisson_extem_surfvy_2x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // field:     potential (scaled by appropriate factors).
  // ext_field: vector potential (scaled by appropriate factors). 
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 
  const double dv11 = 2/dxv[3]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double dv3 = dxv[4], wv3 = w[4]; 
  const double *phi = &field[0]; 
  const double dx10 = 2/dxv[0]; 
  const double dx11 = 2/dxv[1]; 
  const double *A0 = &ext_field[0]; 
  const double *A1 = &ext_field[8]; 
  const double *A2 = &ext_field[16]; 
  double alpha[48] = {0.0}; 

  alpha[0] = 3.464101615137754*A2[2]*dx11*wv3+3.464101615137754*A0[2]*dx11*wv1-3.464101615137754*A1[1]*dx10*wv1-3.464101615137754*phi[2]*dx11; 
  alpha[1] = 3.464101615137754*A2[3]*dx11*wv3+3.464101615137754*A0[3]*dx11*wv1-7.745966692414834*A1[4]*dx10*wv1-3.464101615137754*phi[3]*dx11; 
  alpha[2] = 7.745966692414834*A2[5]*dx11*wv3+7.745966692414834*A0[5]*dx11*wv1-3.464101615137754*A1[3]*dx10*wv1-7.745966692414834*phi[5]*dx11; 
  alpha[3] = A0[2]*dv1*dx11-1.0*A1[1]*dv1*dx10; 
  alpha[4] = A2[2]*dv3*dx11; 
  alpha[5] = 7.745966692414834*A2[7]*dx11*wv3+7.745966692414834*A0[7]*dx11*wv1-7.745966692414834*A1[6]*dx10*wv1-7.745966692414834*phi[7]*dx11; 
  alpha[6] = A0[3]*dv1*dx11-2.23606797749979*A1[4]*dv1*dx10; 
  alpha[7] = 2.23606797749979*A0[5]*dv1*dx11-1.0*A1[3]*dv1*dx10; 
  alpha[8] = A2[3]*dv3*dx11; 
  alpha[9] = 2.23606797749979*A2[5]*dv3*dx11; 
  alpha[11] = 3.464101615137755*A2[6]*dx11*wv3+3.464101615137755*A0[6]*dx11*wv1-3.464101615137755*phi[6]*dx11; 
  alpha[12] = -3.464101615137755*A1[7]*dx10*wv1; 
  alpha[15] = 2.23606797749979*A0[7]*dv1*dx11-2.23606797749979*A1[6]*dv1*dx10; 
  alpha[16] = 2.23606797749979*A2[7]*dv3*dx11; 
  alpha[21] = A0[6]*dv1*dx11; 
  alpha[22] = -1.0*A1[7]*dv1*dx10; 
  alpha[25] = A2[6]*dv3*dx11; 

  double fUpwindQuad_l[81] = {0.0};
  double fUpwindQuad_r[81] = {0.0};
  double fUpwind_l[48] = {0.0};;
  double fUpwind_r[48] = {0.0};
  double Ghat_l[48] = {0.0}; 
  double Ghat_r[48] = {0.0}; 

  if ((-0.2999999999999999*(alpha[25]+alpha[22]+alpha[21]))-0.6037383539249431*(alpha[16]+alpha[15])+0.223606797749979*(alpha[12]+alpha[11])+0.45*(alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5])-0.3354101966249685*(alpha[4]+alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[0] = ser_5x_p2_surfx4_eval_quad_node_0_r(fl); 
    fUpwindQuad_r[0] = ser_5x_p2_surfx4_eval_quad_node_0_r(fc); 
  } else { 
    fUpwindQuad_l[0] = ser_5x_p2_surfx4_eval_quad_node_0_l(fc); 
    fUpwindQuad_r[0] = ser_5x_p2_surfx4_eval_quad_node_0_l(fr); 
  } 
  if ((-0.2999999999999999*(alpha[22]+alpha[21]))-0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])+0.45*(alpha[7]+alpha[6]+alpha[5])-0.3354101966249685*(alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[1] = ser_5x_p2_surfx4_eval_quad_node_1_r(fl); 
    fUpwindQuad_r[1] = ser_5x_p2_surfx4_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_l[1] = ser_5x_p2_surfx4_eval_quad_node_1_l(fc); 
    fUpwindQuad_r[1] = ser_5x_p2_surfx4_eval_quad_node_1_l(fr); 
  } 
  if (0.2999999999999999*alpha[25]-0.2999999999999999*(alpha[22]+alpha[21])+0.6037383539249431*alpha[16]-0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])-0.45*(alpha[9]+alpha[8])+0.45*(alpha[7]+alpha[6]+alpha[5])+0.3354101966249685*alpha[4]-0.3354101966249685*(alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[2] = ser_5x_p2_surfx4_eval_quad_node_2_r(fl); 
    fUpwindQuad_r[2] = ser_5x_p2_surfx4_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_l[2] = ser_5x_p2_surfx4_eval_quad_node_2_l(fc); 
    fUpwindQuad_r[2] = ser_5x_p2_surfx4_eval_quad_node_2_l(fr); 
  } 
  if ((-0.2999999999999999*alpha[25])-0.6037383539249431*alpha[16]+0.223606797749979*(alpha[12]+alpha[11])+0.45*(alpha[9]+alpha[8]+alpha[5])-0.3354101966249685*(alpha[4]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[3] = ser_5x_p2_surfx4_eval_quad_node_3_r(fl); 
    fUpwindQuad_r[3] = ser_5x_p2_surfx4_eval_quad_node_3_r(fc); 
  } else { 
    fUpwindQuad_l[3] = ser_5x_p2_surfx4_eval_quad_node_3_l(fc); 
    fUpwindQuad_r[3] = ser_5x_p2_surfx4_eval_quad_node_3_l(fr); 
  } 
  if (0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[5]-0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[4] = ser_5x_p2_surfx4_eval_quad_node_4_r(fl); 
    fUpwindQuad_r[4] = ser_5x_p2_surfx4_eval_quad_node_4_r(fc); 
  } else { 
    fUpwindQuad_l[4] = ser_5x_p2_surfx4_eval_quad_node_4_l(fc); 
    fUpwindQuad_r[4] = ser_5x_p2_surfx4_eval_quad_node_4_l(fr); 
  } 
  if (0.2999999999999999*alpha[25]+0.6037383539249431*alpha[16]+0.223606797749979*(alpha[12]+alpha[11])-0.45*(alpha[9]+alpha[8])+0.45*alpha[5]+0.3354101966249685*alpha[4]-0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[5] = ser_5x_p2_surfx4_eval_quad_node_5_r(fl); 
    fUpwindQuad_r[5] = ser_5x_p2_surfx4_eval_quad_node_5_r(fc); 
  } else { 
    fUpwindQuad_l[5] = ser_5x_p2_surfx4_eval_quad_node_5_l(fc); 
    fUpwindQuad_r[5] = ser_5x_p2_surfx4_eval_quad_node_5_l(fr); 
  } 
  if ((-0.2999999999999999*alpha[25])+0.2999999999999999*(alpha[22]+alpha[21])-0.6037383539249431*alpha[16]+0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])+0.45*(alpha[9]+alpha[8])-0.45*(alpha[7]+alpha[6])+0.45*alpha[5]-0.3354101966249685*alpha[4]+0.3354101966249685*alpha[3]-0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[6] = ser_5x_p2_surfx4_eval_quad_node_6_r(fl); 
    fUpwindQuad_r[6] = ser_5x_p2_surfx4_eval_quad_node_6_r(fc); 
  } else { 
    fUpwindQuad_l[6] = ser_5x_p2_surfx4_eval_quad_node_6_l(fc); 
    fUpwindQuad_r[6] = ser_5x_p2_surfx4_eval_quad_node_6_l(fr); 
  } 
  if (0.2999999999999999*(alpha[22]+alpha[21])+0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])-0.45*(alpha[7]+alpha[6])+0.45*alpha[5]+0.3354101966249685*alpha[3]-0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[7] = ser_5x_p2_surfx4_eval_quad_node_7_r(fl); 
    fUpwindQuad_r[7] = ser_5x_p2_surfx4_eval_quad_node_7_r(fc); 
  } else { 
    fUpwindQuad_l[7] = ser_5x_p2_surfx4_eval_quad_node_7_l(fc); 
    fUpwindQuad_r[7] = ser_5x_p2_surfx4_eval_quad_node_7_l(fr); 
  } 
  if (0.2999999999999999*(alpha[25]+alpha[22]+alpha[21])+0.6037383539249431*(alpha[16]+alpha[15])+0.223606797749979*(alpha[12]+alpha[11])-0.45*(alpha[9]+alpha[8]+alpha[7]+alpha[6])+0.45*alpha[5]+0.3354101966249685*(alpha[4]+alpha[3])-0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[8] = ser_5x_p2_surfx4_eval_quad_node_8_r(fl); 
    fUpwindQuad_r[8] = ser_5x_p2_surfx4_eval_quad_node_8_r(fc); 
  } else { 
    fUpwindQuad_l[8] = ser_5x_p2_surfx4_eval_quad_node_8_l(fc); 
    fUpwindQuad_r[8] = ser_5x_p2_surfx4_eval_quad_node_8_l(fr); 
  } 
  if ((-0.2999999999999999*alpha[25])+0.375*alpha[22]-0.2999999999999999*alpha[21]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]+0.45*(alpha[8]+alpha[6])-0.3354101966249685*(alpha[4]+alpha[3]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[9] = ser_5x_p2_surfx4_eval_quad_node_9_r(fl); 
    fUpwindQuad_r[9] = ser_5x_p2_surfx4_eval_quad_node_9_r(fc); 
  } else { 
    fUpwindQuad_l[9] = ser_5x_p2_surfx4_eval_quad_node_9_l(fc); 
    fUpwindQuad_r[9] = ser_5x_p2_surfx4_eval_quad_node_9_l(fr); 
  } 
  if (0.375*alpha[22]-0.2999999999999999*alpha[21]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]+0.45*alpha[6]-0.3354101966249685*(alpha[3]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[10] = ser_5x_p2_surfx4_eval_quad_node_10_r(fl); 
    fUpwindQuad_r[10] = ser_5x_p2_surfx4_eval_quad_node_10_r(fc); 
  } else { 
    fUpwindQuad_l[10] = ser_5x_p2_surfx4_eval_quad_node_10_l(fc); 
    fUpwindQuad_r[10] = ser_5x_p2_surfx4_eval_quad_node_10_l(fr); 
  } 
  if (0.2999999999999999*alpha[25]+0.375*alpha[22]-0.2999999999999999*alpha[21]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]-0.45*alpha[8]+0.45*alpha[6]+0.3354101966249685*alpha[4]-0.3354101966249685*(alpha[3]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[11] = ser_5x_p2_surfx4_eval_quad_node_11_r(fl); 
    fUpwindQuad_r[11] = ser_5x_p2_surfx4_eval_quad_node_11_r(fc); 
  } else { 
    fUpwindQuad_l[11] = ser_5x_p2_surfx4_eval_quad_node_11_l(fc); 
    fUpwindQuad_r[11] = ser_5x_p2_surfx4_eval_quad_node_11_l(fr); 
  } 
  if ((-0.2999999999999999*alpha[25])-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]+0.45*alpha[8]-0.3354101966249685*(alpha[4]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[12] = ser_5x_p2_surfx4_eval_quad_node_12_r(fl); 
    fUpwindQuad_r[12] = ser_5x_p2_surfx4_eval_quad_node_12_r(fc); 
  } else { 
    fUpwindQuad_l[12] = ser_5x_p2_surfx4_eval_quad_node_12_l(fc); 
    fUpwindQuad_r[12] = ser_5x_p2_surfx4_eval_quad_node_12_l(fr); 
  } 
  if ((-0.2795084971874737*alpha[12])+0.223606797749979*alpha[11]-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[13] = ser_5x_p2_surfx4_eval_quad_node_13_r(fl); 
    fUpwindQuad_r[13] = ser_5x_p2_surfx4_eval_quad_node_13_r(fc); 
  } else { 
    fUpwindQuad_l[13] = ser_5x_p2_surfx4_eval_quad_node_13_l(fc); 
    fUpwindQuad_r[13] = ser_5x_p2_surfx4_eval_quad_node_13_l(fr); 
  } 
  if (0.2999999999999999*alpha[25]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]-0.45*alpha[8]+0.3354101966249685*alpha[4]-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[14] = ser_5x_p2_surfx4_eval_quad_node_14_r(fl); 
    fUpwindQuad_r[14] = ser_5x_p2_surfx4_eval_quad_node_14_r(fc); 
  } else { 
    fUpwindQuad_l[14] = ser_5x_p2_surfx4_eval_quad_node_14_l(fc); 
    fUpwindQuad_r[14] = ser_5x_p2_surfx4_eval_quad_node_14_l(fr); 
  } 
  if ((-0.2999999999999999*alpha[25])-0.375*alpha[22]+0.2999999999999999*alpha[21]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]+0.45*alpha[8]-0.45*alpha[6]-0.3354101966249685*alpha[4]+0.3354101966249685*alpha[3]-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[15] = ser_5x_p2_surfx4_eval_quad_node_15_r(fl); 
    fUpwindQuad_r[15] = ser_5x_p2_surfx4_eval_quad_node_15_r(fc); 
  } else { 
    fUpwindQuad_l[15] = ser_5x_p2_surfx4_eval_quad_node_15_l(fc); 
    fUpwindQuad_r[15] = ser_5x_p2_surfx4_eval_quad_node_15_l(fr); 
  } 
  if ((-0.375*alpha[22])+0.2999999999999999*alpha[21]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]-0.45*alpha[6]+0.3354101966249685*alpha[3]-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[16] = ser_5x_p2_surfx4_eval_quad_node_16_r(fl); 
    fUpwindQuad_r[16] = ser_5x_p2_surfx4_eval_quad_node_16_r(fc); 
  } else { 
    fUpwindQuad_l[16] = ser_5x_p2_surfx4_eval_quad_node_16_l(fc); 
    fUpwindQuad_r[16] = ser_5x_p2_surfx4_eval_quad_node_16_l(fr); 
  } 
  if (0.2999999999999999*alpha[25]-0.375*alpha[22]+0.2999999999999999*alpha[21]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]-0.45*(alpha[8]+alpha[6])+0.3354101966249685*(alpha[4]+alpha[3])-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[17] = ser_5x_p2_surfx4_eval_quad_node_17_r(fl); 
    fUpwindQuad_r[17] = ser_5x_p2_surfx4_eval_quad_node_17_r(fc); 
  } else { 
    fUpwindQuad_l[17] = ser_5x_p2_surfx4_eval_quad_node_17_l(fc); 
    fUpwindQuad_r[17] = ser_5x_p2_surfx4_eval_quad_node_17_l(fr); 
  } 
  if ((-0.2999999999999999*(alpha[25]+alpha[22]+alpha[21]))+0.6037383539249431*(alpha[16]+alpha[15])+0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[9]+0.45*alpha[8]-0.45*alpha[7]+0.45*alpha[6]-0.45*alpha[5]-0.3354101966249685*(alpha[4]+alpha[3])+0.3354101966249685*alpha[2]-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[18] = ser_5x_p2_surfx4_eval_quad_node_18_r(fl); 
    fUpwindQuad_r[18] = ser_5x_p2_surfx4_eval_quad_node_18_r(fc); 
  } else { 
    fUpwindQuad_l[18] = ser_5x_p2_surfx4_eval_quad_node_18_l(fc); 
    fUpwindQuad_r[18] = ser_5x_p2_surfx4_eval_quad_node_18_l(fr); 
  } 
  if ((-0.2999999999999999*(alpha[22]+alpha[21]))+0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[7]+0.45*alpha[6]-0.45*alpha[5]-0.3354101966249685*alpha[3]+0.3354101966249685*alpha[2]-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[19] = ser_5x_p2_surfx4_eval_quad_node_19_r(fl); 
    fUpwindQuad_r[19] = ser_5x_p2_surfx4_eval_quad_node_19_r(fc); 
  } else { 
    fUpwindQuad_l[19] = ser_5x_p2_surfx4_eval_quad_node_19_l(fc); 
    fUpwindQuad_r[19] = ser_5x_p2_surfx4_eval_quad_node_19_l(fr); 
  } 
  if (0.2999999999999999*alpha[25]-0.2999999999999999*(alpha[22]+alpha[21])-0.6037383539249431*alpha[16]+0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[9]-0.45*(alpha[8]+alpha[7])+0.45*alpha[6]-0.45*alpha[5]+0.3354101966249685*alpha[4]-0.3354101966249685*alpha[3]+0.3354101966249685*alpha[2]-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[20] = ser_5x_p2_surfx4_eval_quad_node_20_r(fl); 
    fUpwindQuad_r[20] = ser_5x_p2_surfx4_eval_quad_node_20_r(fc); 
  } else { 
    fUpwindQuad_l[20] = ser_5x_p2_surfx4_eval_quad_node_20_l(fc); 
    fUpwindQuad_r[20] = ser_5x_p2_surfx4_eval_quad_node_20_l(fr); 
  } 
  if ((-0.2999999999999999*alpha[25])+0.6037383539249431*alpha[16]+0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[9]+0.45*alpha[8]-0.45*alpha[5]-0.3354101966249685*alpha[4]+0.3354101966249685*alpha[2]-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[21] = ser_5x_p2_surfx4_eval_quad_node_21_r(fl); 
    fUpwindQuad_r[21] = ser_5x_p2_surfx4_eval_quad_node_21_r(fc); 
  } else { 
    fUpwindQuad_l[21] = ser_5x_p2_surfx4_eval_quad_node_21_l(fc); 
    fUpwindQuad_r[21] = ser_5x_p2_surfx4_eval_quad_node_21_l(fr); 
  } 
  if (0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[5]+0.3354101966249685*alpha[2]-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[22] = ser_5x_p2_surfx4_eval_quad_node_22_r(fl); 
    fUpwindQuad_r[22] = ser_5x_p2_surfx4_eval_quad_node_22_r(fc); 
  } else { 
    fUpwindQuad_l[22] = ser_5x_p2_surfx4_eval_quad_node_22_l(fc); 
    fUpwindQuad_r[22] = ser_5x_p2_surfx4_eval_quad_node_22_l(fr); 
  } 
  if (0.2999999999999999*alpha[25]-0.6037383539249431*alpha[16]+0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[9]-0.45*(alpha[8]+alpha[5])+0.3354101966249685*(alpha[4]+alpha[2])-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[23] = ser_5x_p2_surfx4_eval_quad_node_23_r(fl); 
    fUpwindQuad_r[23] = ser_5x_p2_surfx4_eval_quad_node_23_r(fc); 
  } else { 
    fUpwindQuad_l[23] = ser_5x_p2_surfx4_eval_quad_node_23_l(fc); 
    fUpwindQuad_r[23] = ser_5x_p2_surfx4_eval_quad_node_23_l(fr); 
  } 
  if ((-0.2999999999999999*alpha[25])+0.2999999999999999*(alpha[22]+alpha[21])+0.6037383539249431*alpha[16]-0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[9]+0.45*(alpha[8]+alpha[7])-0.45*(alpha[6]+alpha[5])-0.3354101966249685*alpha[4]+0.3354101966249685*(alpha[3]+alpha[2])-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[24] = ser_5x_p2_surfx4_eval_quad_node_24_r(fl); 
    fUpwindQuad_r[24] = ser_5x_p2_surfx4_eval_quad_node_24_r(fc); 
  } else { 
    fUpwindQuad_l[24] = ser_5x_p2_surfx4_eval_quad_node_24_l(fc); 
    fUpwindQuad_r[24] = ser_5x_p2_surfx4_eval_quad_node_24_l(fr); 
  } 
  if (0.2999999999999999*(alpha[22]+alpha[21])-0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[7]-0.45*(alpha[6]+alpha[5])+0.3354101966249685*(alpha[3]+alpha[2])-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[25] = ser_5x_p2_surfx4_eval_quad_node_25_r(fl); 
    fUpwindQuad_r[25] = ser_5x_p2_surfx4_eval_quad_node_25_r(fc); 
  } else { 
    fUpwindQuad_l[25] = ser_5x_p2_surfx4_eval_quad_node_25_l(fc); 
    fUpwindQuad_r[25] = ser_5x_p2_surfx4_eval_quad_node_25_l(fr); 
  } 
  if (0.2999999999999999*(alpha[25]+alpha[22]+alpha[21])-0.6037383539249431*(alpha[16]+alpha[15])+0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[9]-0.45*alpha[8]+0.45*alpha[7]-0.45*(alpha[6]+alpha[5])+0.3354101966249685*(alpha[4]+alpha[3]+alpha[2])-0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[26] = ser_5x_p2_surfx4_eval_quad_node_26_r(fl); 
    fUpwindQuad_r[26] = ser_5x_p2_surfx4_eval_quad_node_26_r(fc); 
  } else { 
    fUpwindQuad_l[26] = ser_5x_p2_surfx4_eval_quad_node_26_l(fc); 
    fUpwindQuad_r[26] = ser_5x_p2_surfx4_eval_quad_node_26_l(fr); 
  } 
  if (0.375*alpha[25]-0.2999999999999999*alpha[22]+0.375*alpha[21]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]+0.45*(alpha[9]+alpha[7])-0.3354101966249685*(alpha[4]+alpha[3]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[27] = ser_5x_p2_surfx4_eval_quad_node_27_r(fl); 
    fUpwindQuad_r[27] = ser_5x_p2_surfx4_eval_quad_node_27_r(fc); 
  } else { 
    fUpwindQuad_l[27] = ser_5x_p2_surfx4_eval_quad_node_27_l(fc); 
    fUpwindQuad_r[27] = ser_5x_p2_surfx4_eval_quad_node_27_l(fr); 
  } 
  if ((-0.2999999999999999*alpha[22])+0.375*alpha[21]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]+0.45*alpha[7]-0.3354101966249685*(alpha[3]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[28] = ser_5x_p2_surfx4_eval_quad_node_28_r(fl); 
    fUpwindQuad_r[28] = ser_5x_p2_surfx4_eval_quad_node_28_r(fc); 
  } else { 
    fUpwindQuad_l[28] = ser_5x_p2_surfx4_eval_quad_node_28_l(fc); 
    fUpwindQuad_r[28] = ser_5x_p2_surfx4_eval_quad_node_28_l(fr); 
  } 
  if ((-0.375*alpha[25])-0.2999999999999999*alpha[22]+0.375*alpha[21]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]-0.45*alpha[9]+0.45*alpha[7]+0.3354101966249685*alpha[4]-0.3354101966249685*(alpha[3]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[29] = ser_5x_p2_surfx4_eval_quad_node_29_r(fl); 
    fUpwindQuad_r[29] = ser_5x_p2_surfx4_eval_quad_node_29_r(fc); 
  } else { 
    fUpwindQuad_l[29] = ser_5x_p2_surfx4_eval_quad_node_29_l(fc); 
    fUpwindQuad_r[29] = ser_5x_p2_surfx4_eval_quad_node_29_l(fr); 
  } 
  if (0.375*alpha[25]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]+0.45*alpha[9]-0.3354101966249685*(alpha[4]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[30] = ser_5x_p2_surfx4_eval_quad_node_30_r(fl); 
    fUpwindQuad_r[30] = ser_5x_p2_surfx4_eval_quad_node_30_r(fc); 
  } else { 
    fUpwindQuad_l[30] = ser_5x_p2_surfx4_eval_quad_node_30_l(fc); 
    fUpwindQuad_r[30] = ser_5x_p2_surfx4_eval_quad_node_30_l(fr); 
  } 
  if (0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]-0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[31] = ser_5x_p2_surfx4_eval_quad_node_31_r(fl); 
    fUpwindQuad_r[31] = ser_5x_p2_surfx4_eval_quad_node_31_r(fc); 
  } else { 
    fUpwindQuad_l[31] = ser_5x_p2_surfx4_eval_quad_node_31_l(fc); 
    fUpwindQuad_r[31] = ser_5x_p2_surfx4_eval_quad_node_31_l(fr); 
  } 
  if ((-0.375*alpha[25])+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]-0.45*alpha[9]+0.3354101966249685*alpha[4]-0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[32] = ser_5x_p2_surfx4_eval_quad_node_32_r(fl); 
    fUpwindQuad_r[32] = ser_5x_p2_surfx4_eval_quad_node_32_r(fc); 
  } else { 
    fUpwindQuad_l[32] = ser_5x_p2_surfx4_eval_quad_node_32_l(fc); 
    fUpwindQuad_r[32] = ser_5x_p2_surfx4_eval_quad_node_32_l(fr); 
  } 
  if (0.375*alpha[25]+0.2999999999999999*alpha[22]-0.375*alpha[21]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]+0.45*alpha[9]-0.45*alpha[7]-0.3354101966249685*alpha[4]+0.3354101966249685*alpha[3]-0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[33] = ser_5x_p2_surfx4_eval_quad_node_33_r(fl); 
    fUpwindQuad_r[33] = ser_5x_p2_surfx4_eval_quad_node_33_r(fc); 
  } else { 
    fUpwindQuad_l[33] = ser_5x_p2_surfx4_eval_quad_node_33_l(fc); 
    fUpwindQuad_r[33] = ser_5x_p2_surfx4_eval_quad_node_33_l(fr); 
  } 
  if (0.2999999999999999*alpha[22]-0.375*alpha[21]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]-0.45*alpha[7]+0.3354101966249685*alpha[3]-0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[34] = ser_5x_p2_surfx4_eval_quad_node_34_r(fl); 
    fUpwindQuad_r[34] = ser_5x_p2_surfx4_eval_quad_node_34_r(fc); 
  } else { 
    fUpwindQuad_l[34] = ser_5x_p2_surfx4_eval_quad_node_34_l(fc); 
    fUpwindQuad_r[34] = ser_5x_p2_surfx4_eval_quad_node_34_l(fr); 
  } 
  if ((-0.375*alpha[25])+0.2999999999999999*alpha[22]-0.375*alpha[21]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]-0.45*(alpha[9]+alpha[7])+0.3354101966249685*(alpha[4]+alpha[3])-0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[35] = ser_5x_p2_surfx4_eval_quad_node_35_r(fl); 
    fUpwindQuad_r[35] = ser_5x_p2_surfx4_eval_quad_node_35_r(fc); 
  } else { 
    fUpwindQuad_l[35] = ser_5x_p2_surfx4_eval_quad_node_35_l(fc); 
    fUpwindQuad_r[35] = ser_5x_p2_surfx4_eval_quad_node_35_l(fr); 
  } 
  if (0.375*(alpha[25]+alpha[22]+alpha[21])-0.2795084971874737*(alpha[12]+alpha[11])-0.3354101966249685*(alpha[4]+alpha[3])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[36] = ser_5x_p2_surfx4_eval_quad_node_36_r(fl); 
    fUpwindQuad_r[36] = ser_5x_p2_surfx4_eval_quad_node_36_r(fc); 
  } else { 
    fUpwindQuad_l[36] = ser_5x_p2_surfx4_eval_quad_node_36_l(fc); 
    fUpwindQuad_r[36] = ser_5x_p2_surfx4_eval_quad_node_36_l(fr); 
  } 
  if (0.375*(alpha[22]+alpha[21])-0.2795084971874737*(alpha[12]+alpha[11])-0.3354101966249685*alpha[3]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[37] = ser_5x_p2_surfx4_eval_quad_node_37_r(fl); 
    fUpwindQuad_r[37] = ser_5x_p2_surfx4_eval_quad_node_37_r(fc); 
  } else { 
    fUpwindQuad_l[37] = ser_5x_p2_surfx4_eval_quad_node_37_l(fc); 
    fUpwindQuad_r[37] = ser_5x_p2_surfx4_eval_quad_node_37_l(fr); 
  } 
  if ((-0.375*alpha[25])+0.375*(alpha[22]+alpha[21])-0.2795084971874737*(alpha[12]+alpha[11])+0.3354101966249685*alpha[4]-0.3354101966249685*alpha[3]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[38] = ser_5x_p2_surfx4_eval_quad_node_38_r(fl); 
    fUpwindQuad_r[38] = ser_5x_p2_surfx4_eval_quad_node_38_r(fc); 
  } else { 
    fUpwindQuad_l[38] = ser_5x_p2_surfx4_eval_quad_node_38_l(fc); 
    fUpwindQuad_r[38] = ser_5x_p2_surfx4_eval_quad_node_38_l(fr); 
  } 
  if (0.375*alpha[25]-0.2795084971874737*(alpha[12]+alpha[11])-0.3354101966249685*alpha[4]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[39] = ser_5x_p2_surfx4_eval_quad_node_39_r(fl); 
    fUpwindQuad_r[39] = ser_5x_p2_surfx4_eval_quad_node_39_r(fc); 
  } else { 
    fUpwindQuad_l[39] = ser_5x_p2_surfx4_eval_quad_node_39_l(fc); 
    fUpwindQuad_r[39] = ser_5x_p2_surfx4_eval_quad_node_39_l(fr); 
  } 
  if (0.25*alpha[0]-0.2795084971874737*(alpha[12]+alpha[11]) > 0) { 
    fUpwindQuad_l[40] = ser_5x_p2_surfx4_eval_quad_node_40_r(fl); 
    fUpwindQuad_r[40] = ser_5x_p2_surfx4_eval_quad_node_40_r(fc); 
  } else { 
    fUpwindQuad_l[40] = ser_5x_p2_surfx4_eval_quad_node_40_l(fc); 
    fUpwindQuad_r[40] = ser_5x_p2_surfx4_eval_quad_node_40_l(fr); 
  } 
  if ((-0.375*alpha[25])-0.2795084971874737*(alpha[12]+alpha[11])+0.3354101966249685*alpha[4]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[41] = ser_5x_p2_surfx4_eval_quad_node_41_r(fl); 
    fUpwindQuad_r[41] = ser_5x_p2_surfx4_eval_quad_node_41_r(fc); 
  } else { 
    fUpwindQuad_l[41] = ser_5x_p2_surfx4_eval_quad_node_41_l(fc); 
    fUpwindQuad_r[41] = ser_5x_p2_surfx4_eval_quad_node_41_l(fr); 
  } 
  if (0.375*alpha[25]-0.375*(alpha[22]+alpha[21])-0.2795084971874737*(alpha[12]+alpha[11])-0.3354101966249685*alpha[4]+0.3354101966249685*alpha[3]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[42] = ser_5x_p2_surfx4_eval_quad_node_42_r(fl); 
    fUpwindQuad_r[42] = ser_5x_p2_surfx4_eval_quad_node_42_r(fc); 
  } else { 
    fUpwindQuad_l[42] = ser_5x_p2_surfx4_eval_quad_node_42_l(fc); 
    fUpwindQuad_r[42] = ser_5x_p2_surfx4_eval_quad_node_42_l(fr); 
  } 
  if ((-0.375*(alpha[22]+alpha[21]))-0.2795084971874737*(alpha[12]+alpha[11])+0.3354101966249685*alpha[3]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[43] = ser_5x_p2_surfx4_eval_quad_node_43_r(fl); 
    fUpwindQuad_r[43] = ser_5x_p2_surfx4_eval_quad_node_43_r(fc); 
  } else { 
    fUpwindQuad_l[43] = ser_5x_p2_surfx4_eval_quad_node_43_l(fc); 
    fUpwindQuad_r[43] = ser_5x_p2_surfx4_eval_quad_node_43_l(fr); 
  } 
  if ((-0.375*(alpha[25]+alpha[22]+alpha[21]))-0.2795084971874737*(alpha[12]+alpha[11])+0.3354101966249685*(alpha[4]+alpha[3])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[44] = ser_5x_p2_surfx4_eval_quad_node_44_r(fl); 
    fUpwindQuad_r[44] = ser_5x_p2_surfx4_eval_quad_node_44_r(fc); 
  } else { 
    fUpwindQuad_l[44] = ser_5x_p2_surfx4_eval_quad_node_44_l(fc); 
    fUpwindQuad_r[44] = ser_5x_p2_surfx4_eval_quad_node_44_l(fr); 
  } 
  if (0.375*alpha[25]-0.2999999999999999*alpha[22]+0.375*alpha[21]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]-0.45*(alpha[9]+alpha[7])-0.3354101966249685*(alpha[4]+alpha[3])+0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[45] = ser_5x_p2_surfx4_eval_quad_node_45_r(fl); 
    fUpwindQuad_r[45] = ser_5x_p2_surfx4_eval_quad_node_45_r(fc); 
  } else { 
    fUpwindQuad_l[45] = ser_5x_p2_surfx4_eval_quad_node_45_l(fc); 
    fUpwindQuad_r[45] = ser_5x_p2_surfx4_eval_quad_node_45_l(fr); 
  } 
  if ((-0.2999999999999999*alpha[22])+0.375*alpha[21]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]-0.45*alpha[7]-0.3354101966249685*alpha[3]+0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[46] = ser_5x_p2_surfx4_eval_quad_node_46_r(fl); 
    fUpwindQuad_r[46] = ser_5x_p2_surfx4_eval_quad_node_46_r(fc); 
  } else { 
    fUpwindQuad_l[46] = ser_5x_p2_surfx4_eval_quad_node_46_l(fc); 
    fUpwindQuad_r[46] = ser_5x_p2_surfx4_eval_quad_node_46_l(fr); 
  } 
  if ((-0.375*alpha[25])-0.2999999999999999*alpha[22]+0.375*alpha[21]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]+0.45*alpha[9]-0.45*alpha[7]+0.3354101966249685*alpha[4]-0.3354101966249685*alpha[3]+0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[47] = ser_5x_p2_surfx4_eval_quad_node_47_r(fl); 
    fUpwindQuad_r[47] = ser_5x_p2_surfx4_eval_quad_node_47_r(fc); 
  } else { 
    fUpwindQuad_l[47] = ser_5x_p2_surfx4_eval_quad_node_47_l(fc); 
    fUpwindQuad_r[47] = ser_5x_p2_surfx4_eval_quad_node_47_l(fr); 
  } 
  if (0.375*alpha[25]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]-0.45*alpha[9]-0.3354101966249685*alpha[4]+0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[48] = ser_5x_p2_surfx4_eval_quad_node_48_r(fl); 
    fUpwindQuad_r[48] = ser_5x_p2_surfx4_eval_quad_node_48_r(fc); 
  } else { 
    fUpwindQuad_l[48] = ser_5x_p2_surfx4_eval_quad_node_48_l(fc); 
    fUpwindQuad_r[48] = ser_5x_p2_surfx4_eval_quad_node_48_l(fr); 
  } 
  if (0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]+0.3354101966249685*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[49] = ser_5x_p2_surfx4_eval_quad_node_49_r(fl); 
    fUpwindQuad_r[49] = ser_5x_p2_surfx4_eval_quad_node_49_r(fc); 
  } else { 
    fUpwindQuad_l[49] = ser_5x_p2_surfx4_eval_quad_node_49_l(fc); 
    fUpwindQuad_r[49] = ser_5x_p2_surfx4_eval_quad_node_49_l(fr); 
  } 
  if ((-0.375*alpha[25])+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]+0.45*alpha[9]+0.3354101966249685*(alpha[4]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[50] = ser_5x_p2_surfx4_eval_quad_node_50_r(fl); 
    fUpwindQuad_r[50] = ser_5x_p2_surfx4_eval_quad_node_50_r(fc); 
  } else { 
    fUpwindQuad_l[50] = ser_5x_p2_surfx4_eval_quad_node_50_l(fc); 
    fUpwindQuad_r[50] = ser_5x_p2_surfx4_eval_quad_node_50_l(fr); 
  } 
  if (0.375*alpha[25]+0.2999999999999999*alpha[22]-0.375*alpha[21]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]-0.45*alpha[9]+0.45*alpha[7]-0.3354101966249685*alpha[4]+0.3354101966249685*(alpha[3]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[51] = ser_5x_p2_surfx4_eval_quad_node_51_r(fl); 
    fUpwindQuad_r[51] = ser_5x_p2_surfx4_eval_quad_node_51_r(fc); 
  } else { 
    fUpwindQuad_l[51] = ser_5x_p2_surfx4_eval_quad_node_51_l(fc); 
    fUpwindQuad_r[51] = ser_5x_p2_surfx4_eval_quad_node_51_l(fr); 
  } 
  if (0.2999999999999999*alpha[22]-0.375*alpha[21]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]+0.45*alpha[7]+0.3354101966249685*(alpha[3]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[52] = ser_5x_p2_surfx4_eval_quad_node_52_r(fl); 
    fUpwindQuad_r[52] = ser_5x_p2_surfx4_eval_quad_node_52_r(fc); 
  } else { 
    fUpwindQuad_l[52] = ser_5x_p2_surfx4_eval_quad_node_52_l(fc); 
    fUpwindQuad_r[52] = ser_5x_p2_surfx4_eval_quad_node_52_l(fr); 
  } 
  if ((-0.375*alpha[25])+0.2999999999999999*alpha[22]-0.375*alpha[21]+0.223606797749979*alpha[12]-0.2795084971874737*alpha[11]+0.45*(alpha[9]+alpha[7])+0.3354101966249685*(alpha[4]+alpha[3]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[53] = ser_5x_p2_surfx4_eval_quad_node_53_r(fl); 
    fUpwindQuad_r[53] = ser_5x_p2_surfx4_eval_quad_node_53_r(fc); 
  } else { 
    fUpwindQuad_l[53] = ser_5x_p2_surfx4_eval_quad_node_53_l(fc); 
    fUpwindQuad_r[53] = ser_5x_p2_surfx4_eval_quad_node_53_l(fr); 
  } 
  if ((-0.2999999999999999*(alpha[25]+alpha[22]+alpha[21]))+0.6037383539249431*(alpha[16]+alpha[15])+0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[9]-0.45*alpha[8]+0.45*alpha[7]-0.45*(alpha[6]+alpha[5])-0.3354101966249685*(alpha[4]+alpha[3]+alpha[2])+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[54] = ser_5x_p2_surfx4_eval_quad_node_54_r(fl); 
    fUpwindQuad_r[54] = ser_5x_p2_surfx4_eval_quad_node_54_r(fc); 
  } else { 
    fUpwindQuad_l[54] = ser_5x_p2_surfx4_eval_quad_node_54_l(fc); 
    fUpwindQuad_r[54] = ser_5x_p2_surfx4_eval_quad_node_54_l(fr); 
  } 
  if ((-0.2999999999999999*(alpha[22]+alpha[21]))+0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[7]-0.45*(alpha[6]+alpha[5])-0.3354101966249685*(alpha[3]+alpha[2])+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[55] = ser_5x_p2_surfx4_eval_quad_node_55_r(fl); 
    fUpwindQuad_r[55] = ser_5x_p2_surfx4_eval_quad_node_55_r(fc); 
  } else { 
    fUpwindQuad_l[55] = ser_5x_p2_surfx4_eval_quad_node_55_l(fc); 
    fUpwindQuad_r[55] = ser_5x_p2_surfx4_eval_quad_node_55_l(fr); 
  } 
  if (0.2999999999999999*alpha[25]-0.2999999999999999*(alpha[22]+alpha[21])-0.6037383539249431*alpha[16]+0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[9]+0.45*(alpha[8]+alpha[7])-0.45*(alpha[6]+alpha[5])+0.3354101966249685*alpha[4]-0.3354101966249685*(alpha[3]+alpha[2])+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[56] = ser_5x_p2_surfx4_eval_quad_node_56_r(fl); 
    fUpwindQuad_r[56] = ser_5x_p2_surfx4_eval_quad_node_56_r(fc); 
  } else { 
    fUpwindQuad_l[56] = ser_5x_p2_surfx4_eval_quad_node_56_l(fc); 
    fUpwindQuad_r[56] = ser_5x_p2_surfx4_eval_quad_node_56_l(fr); 
  } 
  if ((-0.2999999999999999*alpha[25])+0.6037383539249431*alpha[16]+0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[9]-0.45*(alpha[8]+alpha[5])-0.3354101966249685*(alpha[4]+alpha[2])+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[57] = ser_5x_p2_surfx4_eval_quad_node_57_r(fl); 
    fUpwindQuad_r[57] = ser_5x_p2_surfx4_eval_quad_node_57_r(fc); 
  } else { 
    fUpwindQuad_l[57] = ser_5x_p2_surfx4_eval_quad_node_57_l(fc); 
    fUpwindQuad_r[57] = ser_5x_p2_surfx4_eval_quad_node_57_l(fr); 
  } 
  if (0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[5]-0.3354101966249685*alpha[2]+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[58] = ser_5x_p2_surfx4_eval_quad_node_58_r(fl); 
    fUpwindQuad_r[58] = ser_5x_p2_surfx4_eval_quad_node_58_r(fc); 
  } else { 
    fUpwindQuad_l[58] = ser_5x_p2_surfx4_eval_quad_node_58_l(fc); 
    fUpwindQuad_r[58] = ser_5x_p2_surfx4_eval_quad_node_58_l(fr); 
  } 
  if (0.2999999999999999*alpha[25]-0.6037383539249431*alpha[16]+0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[9]+0.45*alpha[8]-0.45*alpha[5]+0.3354101966249685*alpha[4]-0.3354101966249685*alpha[2]+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[59] = ser_5x_p2_surfx4_eval_quad_node_59_r(fl); 
    fUpwindQuad_r[59] = ser_5x_p2_surfx4_eval_quad_node_59_r(fc); 
  } else { 
    fUpwindQuad_l[59] = ser_5x_p2_surfx4_eval_quad_node_59_l(fc); 
    fUpwindQuad_r[59] = ser_5x_p2_surfx4_eval_quad_node_59_l(fr); 
  } 
  if ((-0.2999999999999999*alpha[25])+0.2999999999999999*(alpha[22]+alpha[21])+0.6037383539249431*alpha[16]-0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[9]-0.45*(alpha[8]+alpha[7])+0.45*alpha[6]-0.45*alpha[5]-0.3354101966249685*alpha[4]+0.3354101966249685*alpha[3]-0.3354101966249685*alpha[2]+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[60] = ser_5x_p2_surfx4_eval_quad_node_60_r(fl); 
    fUpwindQuad_r[60] = ser_5x_p2_surfx4_eval_quad_node_60_r(fc); 
  } else { 
    fUpwindQuad_l[60] = ser_5x_p2_surfx4_eval_quad_node_60_l(fc); 
    fUpwindQuad_r[60] = ser_5x_p2_surfx4_eval_quad_node_60_l(fr); 
  } 
  if (0.2999999999999999*(alpha[22]+alpha[21])-0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[7]+0.45*alpha[6]-0.45*alpha[5]+0.3354101966249685*alpha[3]-0.3354101966249685*alpha[2]+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[61] = ser_5x_p2_surfx4_eval_quad_node_61_r(fl); 
    fUpwindQuad_r[61] = ser_5x_p2_surfx4_eval_quad_node_61_r(fc); 
  } else { 
    fUpwindQuad_l[61] = ser_5x_p2_surfx4_eval_quad_node_61_l(fc); 
    fUpwindQuad_r[61] = ser_5x_p2_surfx4_eval_quad_node_61_l(fr); 
  } 
  if (0.2999999999999999*(alpha[25]+alpha[22]+alpha[21])-0.6037383539249431*(alpha[16]+alpha[15])+0.223606797749979*(alpha[12]+alpha[11])-0.45*alpha[9]+0.45*alpha[8]-0.45*alpha[7]+0.45*alpha[6]-0.45*alpha[5]+0.3354101966249685*(alpha[4]+alpha[3])-0.3354101966249685*alpha[2]+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[62] = ser_5x_p2_surfx4_eval_quad_node_62_r(fl); 
    fUpwindQuad_r[62] = ser_5x_p2_surfx4_eval_quad_node_62_r(fc); 
  } else { 
    fUpwindQuad_l[62] = ser_5x_p2_surfx4_eval_quad_node_62_l(fc); 
    fUpwindQuad_r[62] = ser_5x_p2_surfx4_eval_quad_node_62_l(fr); 
  } 
  if ((-0.2999999999999999*alpha[25])+0.375*alpha[22]-0.2999999999999999*alpha[21]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]-0.45*(alpha[8]+alpha[6])-0.3354101966249685*(alpha[4]+alpha[3])+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[63] = ser_5x_p2_surfx4_eval_quad_node_63_r(fl); 
    fUpwindQuad_r[63] = ser_5x_p2_surfx4_eval_quad_node_63_r(fc); 
  } else { 
    fUpwindQuad_l[63] = ser_5x_p2_surfx4_eval_quad_node_63_l(fc); 
    fUpwindQuad_r[63] = ser_5x_p2_surfx4_eval_quad_node_63_l(fr); 
  } 
  if (0.375*alpha[22]-0.2999999999999999*alpha[21]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]-0.45*alpha[6]-0.3354101966249685*alpha[3]+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[64] = ser_5x_p2_surfx4_eval_quad_node_64_r(fl); 
    fUpwindQuad_r[64] = ser_5x_p2_surfx4_eval_quad_node_64_r(fc); 
  } else { 
    fUpwindQuad_l[64] = ser_5x_p2_surfx4_eval_quad_node_64_l(fc); 
    fUpwindQuad_r[64] = ser_5x_p2_surfx4_eval_quad_node_64_l(fr); 
  } 
  if (0.2999999999999999*alpha[25]+0.375*alpha[22]-0.2999999999999999*alpha[21]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]+0.45*alpha[8]-0.45*alpha[6]+0.3354101966249685*alpha[4]-0.3354101966249685*alpha[3]+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[65] = ser_5x_p2_surfx4_eval_quad_node_65_r(fl); 
    fUpwindQuad_r[65] = ser_5x_p2_surfx4_eval_quad_node_65_r(fc); 
  } else { 
    fUpwindQuad_l[65] = ser_5x_p2_surfx4_eval_quad_node_65_l(fc); 
    fUpwindQuad_r[65] = ser_5x_p2_surfx4_eval_quad_node_65_l(fr); 
  } 
  if ((-0.2999999999999999*alpha[25])-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]-0.45*alpha[8]-0.3354101966249685*alpha[4]+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[66] = ser_5x_p2_surfx4_eval_quad_node_66_r(fl); 
    fUpwindQuad_r[66] = ser_5x_p2_surfx4_eval_quad_node_66_r(fc); 
  } else { 
    fUpwindQuad_l[66] = ser_5x_p2_surfx4_eval_quad_node_66_l(fc); 
    fUpwindQuad_r[66] = ser_5x_p2_surfx4_eval_quad_node_66_l(fr); 
  } 
  if ((-0.2795084971874737*alpha[12])+0.223606797749979*alpha[11]+0.3354101966249685*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[67] = ser_5x_p2_surfx4_eval_quad_node_67_r(fl); 
    fUpwindQuad_r[67] = ser_5x_p2_surfx4_eval_quad_node_67_r(fc); 
  } else { 
    fUpwindQuad_l[67] = ser_5x_p2_surfx4_eval_quad_node_67_l(fc); 
    fUpwindQuad_r[67] = ser_5x_p2_surfx4_eval_quad_node_67_l(fr); 
  } 
  if (0.2999999999999999*alpha[25]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]+0.45*alpha[8]+0.3354101966249685*(alpha[4]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[68] = ser_5x_p2_surfx4_eval_quad_node_68_r(fl); 
    fUpwindQuad_r[68] = ser_5x_p2_surfx4_eval_quad_node_68_r(fc); 
  } else { 
    fUpwindQuad_l[68] = ser_5x_p2_surfx4_eval_quad_node_68_l(fc); 
    fUpwindQuad_r[68] = ser_5x_p2_surfx4_eval_quad_node_68_l(fr); 
  } 
  if ((-0.2999999999999999*alpha[25])-0.375*alpha[22]+0.2999999999999999*alpha[21]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]-0.45*alpha[8]+0.45*alpha[6]-0.3354101966249685*alpha[4]+0.3354101966249685*(alpha[3]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[69] = ser_5x_p2_surfx4_eval_quad_node_69_r(fl); 
    fUpwindQuad_r[69] = ser_5x_p2_surfx4_eval_quad_node_69_r(fc); 
  } else { 
    fUpwindQuad_l[69] = ser_5x_p2_surfx4_eval_quad_node_69_l(fc); 
    fUpwindQuad_r[69] = ser_5x_p2_surfx4_eval_quad_node_69_l(fr); 
  } 
  if ((-0.375*alpha[22])+0.2999999999999999*alpha[21]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]+0.45*alpha[6]+0.3354101966249685*(alpha[3]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[70] = ser_5x_p2_surfx4_eval_quad_node_70_r(fl); 
    fUpwindQuad_r[70] = ser_5x_p2_surfx4_eval_quad_node_70_r(fc); 
  } else { 
    fUpwindQuad_l[70] = ser_5x_p2_surfx4_eval_quad_node_70_l(fc); 
    fUpwindQuad_r[70] = ser_5x_p2_surfx4_eval_quad_node_70_l(fr); 
  } 
  if (0.2999999999999999*alpha[25]-0.375*alpha[22]+0.2999999999999999*alpha[21]-0.2795084971874737*alpha[12]+0.223606797749979*alpha[11]+0.45*(alpha[8]+alpha[6])+0.3354101966249685*(alpha[4]+alpha[3]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[71] = ser_5x_p2_surfx4_eval_quad_node_71_r(fl); 
    fUpwindQuad_r[71] = ser_5x_p2_surfx4_eval_quad_node_71_r(fc); 
  } else { 
    fUpwindQuad_l[71] = ser_5x_p2_surfx4_eval_quad_node_71_l(fc); 
    fUpwindQuad_r[71] = ser_5x_p2_surfx4_eval_quad_node_71_l(fr); 
  } 
  if ((-0.2999999999999999*(alpha[25]+alpha[22]+alpha[21]))-0.6037383539249431*(alpha[16]+alpha[15])+0.223606797749979*(alpha[12]+alpha[11])-0.45*(alpha[9]+alpha[8]+alpha[7]+alpha[6])+0.45*alpha[5]-0.3354101966249685*(alpha[4]+alpha[3])+0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[72] = ser_5x_p2_surfx4_eval_quad_node_72_r(fl); 
    fUpwindQuad_r[72] = ser_5x_p2_surfx4_eval_quad_node_72_r(fc); 
  } else { 
    fUpwindQuad_l[72] = ser_5x_p2_surfx4_eval_quad_node_72_l(fc); 
    fUpwindQuad_r[72] = ser_5x_p2_surfx4_eval_quad_node_72_l(fr); 
  } 
  if ((-0.2999999999999999*(alpha[22]+alpha[21]))-0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])-0.45*(alpha[7]+alpha[6])+0.45*alpha[5]-0.3354101966249685*alpha[3]+0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[73] = ser_5x_p2_surfx4_eval_quad_node_73_r(fl); 
    fUpwindQuad_r[73] = ser_5x_p2_surfx4_eval_quad_node_73_r(fc); 
  } else { 
    fUpwindQuad_l[73] = ser_5x_p2_surfx4_eval_quad_node_73_l(fc); 
    fUpwindQuad_r[73] = ser_5x_p2_surfx4_eval_quad_node_73_l(fr); 
  } 
  if (0.2999999999999999*alpha[25]-0.2999999999999999*(alpha[22]+alpha[21])+0.6037383539249431*alpha[16]-0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])+0.45*(alpha[9]+alpha[8])-0.45*(alpha[7]+alpha[6])+0.45*alpha[5]+0.3354101966249685*alpha[4]-0.3354101966249685*alpha[3]+0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[74] = ser_5x_p2_surfx4_eval_quad_node_74_r(fl); 
    fUpwindQuad_r[74] = ser_5x_p2_surfx4_eval_quad_node_74_r(fc); 
  } else { 
    fUpwindQuad_l[74] = ser_5x_p2_surfx4_eval_quad_node_74_l(fc); 
    fUpwindQuad_r[74] = ser_5x_p2_surfx4_eval_quad_node_74_l(fr); 
  } 
  if ((-0.2999999999999999*alpha[25])-0.6037383539249431*alpha[16]+0.223606797749979*(alpha[12]+alpha[11])-0.45*(alpha[9]+alpha[8])+0.45*alpha[5]-0.3354101966249685*alpha[4]+0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[75] = ser_5x_p2_surfx4_eval_quad_node_75_r(fl); 
    fUpwindQuad_r[75] = ser_5x_p2_surfx4_eval_quad_node_75_r(fc); 
  } else { 
    fUpwindQuad_l[75] = ser_5x_p2_surfx4_eval_quad_node_75_l(fc); 
    fUpwindQuad_r[75] = ser_5x_p2_surfx4_eval_quad_node_75_l(fr); 
  } 
  if (0.223606797749979*(alpha[12]+alpha[11])+0.45*alpha[5]+0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[76] = ser_5x_p2_surfx4_eval_quad_node_76_r(fl); 
    fUpwindQuad_r[76] = ser_5x_p2_surfx4_eval_quad_node_76_r(fc); 
  } else { 
    fUpwindQuad_l[76] = ser_5x_p2_surfx4_eval_quad_node_76_l(fc); 
    fUpwindQuad_r[76] = ser_5x_p2_surfx4_eval_quad_node_76_l(fr); 
  } 
  if (0.2999999999999999*alpha[25]+0.6037383539249431*alpha[16]+0.223606797749979*(alpha[12]+alpha[11])+0.45*(alpha[9]+alpha[8]+alpha[5])+0.3354101966249685*(alpha[4]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[77] = ser_5x_p2_surfx4_eval_quad_node_77_r(fl); 
    fUpwindQuad_r[77] = ser_5x_p2_surfx4_eval_quad_node_77_r(fc); 
  } else { 
    fUpwindQuad_l[77] = ser_5x_p2_surfx4_eval_quad_node_77_l(fc); 
    fUpwindQuad_r[77] = ser_5x_p2_surfx4_eval_quad_node_77_l(fr); 
  } 
  if ((-0.2999999999999999*alpha[25])+0.2999999999999999*(alpha[22]+alpha[21])-0.6037383539249431*alpha[16]+0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])-0.45*(alpha[9]+alpha[8])+0.45*(alpha[7]+alpha[6]+alpha[5])-0.3354101966249685*alpha[4]+0.3354101966249685*(alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[78] = ser_5x_p2_surfx4_eval_quad_node_78_r(fl); 
    fUpwindQuad_r[78] = ser_5x_p2_surfx4_eval_quad_node_78_r(fc); 
  } else { 
    fUpwindQuad_l[78] = ser_5x_p2_surfx4_eval_quad_node_78_l(fc); 
    fUpwindQuad_r[78] = ser_5x_p2_surfx4_eval_quad_node_78_l(fr); 
  } 
  if (0.2999999999999999*(alpha[22]+alpha[21])+0.6037383539249431*alpha[15]+0.223606797749979*(alpha[12]+alpha[11])+0.45*(alpha[7]+alpha[6]+alpha[5])+0.3354101966249685*(alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[79] = ser_5x_p2_surfx4_eval_quad_node_79_r(fl); 
    fUpwindQuad_r[79] = ser_5x_p2_surfx4_eval_quad_node_79_r(fc); 
  } else { 
    fUpwindQuad_l[79] = ser_5x_p2_surfx4_eval_quad_node_79_l(fc); 
    fUpwindQuad_r[79] = ser_5x_p2_surfx4_eval_quad_node_79_l(fr); 
  } 
  if (0.2999999999999999*(alpha[25]+alpha[22]+alpha[21])+0.6037383539249431*(alpha[16]+alpha[15])+0.223606797749979*(alpha[12]+alpha[11])+0.45*(alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5])+0.3354101966249685*(alpha[4]+alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[80] = ser_5x_p2_surfx4_eval_quad_node_80_r(fl); 
    fUpwindQuad_r[80] = ser_5x_p2_surfx4_eval_quad_node_80_r(fc); 
  } else { 
    fUpwindQuad_l[80] = ser_5x_p2_surfx4_eval_quad_node_80_l(fc); 
    fUpwindQuad_r[80] = ser_5x_p2_surfx4_eval_quad_node_80_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_5x_p2_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  ser_5x_p2_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.25*alpha[25]*fUpwind_l[25]+0.25*alpha[22]*fUpwind_l[22]+0.25*alpha[21]*fUpwind_l[21]+0.25*alpha[16]*fUpwind_l[16]+0.25*alpha[15]*fUpwind_l[15]+0.25*alpha[12]*fUpwind_l[12]+0.25*alpha[11]*fUpwind_l[11]+0.25*alpha[9]*fUpwind_l[9]+0.25*alpha[8]*fUpwind_l[8]+0.25*alpha[7]*fUpwind_l[7]+0.25*alpha[6]*fUpwind_l[6]+0.25*alpha[5]*fUpwind_l[5]+0.25*alpha[4]*fUpwind_l[4]+0.25*alpha[3]*fUpwind_l[3]+0.25*alpha[2]*fUpwind_l[2]+0.25*alpha[1]*fUpwind_l[1]+0.25*alpha[0]*fUpwind_l[0]; 
  Ghat_l[1] = 0.223606797749979*alpha[16]*fUpwind_l[35]+0.2500000000000001*alpha[22]*fUpwind_l[33]+0.223606797749979*alpha[15]*fUpwind_l[32]+0.223606797749979*alpha[8]*fUpwind_l[25]+0.223606797749979*fUpwind_l[8]*alpha[25]+0.223606797749979*alpha[6]*fUpwind_l[21]+0.223606797749979*fUpwind_l[6]*alpha[21]+0.2500000000000001*alpha[12]*fUpwind_l[20]+0.223606797749979*alpha[5]*fUpwind_l[19]+0.25*alpha[9]*fUpwind_l[16]+0.25*fUpwind_l[9]*alpha[16]+0.25*alpha[7]*fUpwind_l[15]+0.25*fUpwind_l[7]*alpha[15]+0.223606797749979*alpha[1]*fUpwind_l[11]+0.223606797749979*fUpwind_l[1]*alpha[11]+0.25*alpha[4]*fUpwind_l[8]+0.25*fUpwind_l[4]*alpha[8]+0.25*alpha[3]*fUpwind_l[6]+0.25*fUpwind_l[3]*alpha[6]+0.25*alpha[2]*fUpwind_l[5]+0.25*fUpwind_l[2]*alpha[5]+0.25*alpha[0]*fUpwind_l[1]+0.25*fUpwind_l[0]*alpha[1]; 
  Ghat_l[2] = 0.223606797749979*alpha[16]*fUpwind_l[36]+0.2500000000000001*alpha[25]*fUpwind_l[35]+0.223606797749979*alpha[15]*fUpwind_l[33]+0.2500000000000001*alpha[21]*fUpwind_l[32]+0.223606797749979*alpha[9]*fUpwind_l[26]+0.223606797749979*alpha[7]*fUpwind_l[22]+0.223606797749979*fUpwind_l[7]*alpha[22]+0.223606797749979*alpha[5]*fUpwind_l[20]+0.2500000000000001*alpha[11]*fUpwind_l[19]+0.25*alpha[8]*fUpwind_l[16]+0.25*fUpwind_l[8]*alpha[16]+0.25*alpha[6]*fUpwind_l[15]+0.25*fUpwind_l[6]*alpha[15]+0.223606797749979*alpha[2]*fUpwind_l[12]+0.223606797749979*fUpwind_l[2]*alpha[12]+0.25*alpha[4]*fUpwind_l[9]+0.25*fUpwind_l[4]*alpha[9]+0.25*alpha[3]*fUpwind_l[7]+0.25*fUpwind_l[3]*alpha[7]+0.25*alpha[1]*fUpwind_l[5]+0.25*fUpwind_l[1]*alpha[5]+0.25*alpha[0]*fUpwind_l[2]+0.25*fUpwind_l[0]*alpha[2]; 
  Ghat_l[3] = 0.2500000000000001*alpha[25]*fUpwind_l[37]+0.223606797749979*alpha[15]*fUpwind_l[34]+0.25*alpha[16]*fUpwind_l[31]+0.223606797749979*alpha[7]*fUpwind_l[24]+0.223606797749979*alpha[6]*fUpwind_l[23]+0.2500000000000001*alpha[12]*fUpwind_l[22]+0.2500000000000001*fUpwind_l[12]*alpha[22]+0.2500000000000001*alpha[11]*fUpwind_l[21]+0.2500000000000001*fUpwind_l[11]*alpha[21]+0.25*alpha[9]*fUpwind_l[18]+0.25*alpha[8]*fUpwind_l[17]+0.25*alpha[5]*fUpwind_l[15]+0.25*fUpwind_l[5]*alpha[15]+0.223606797749979*alpha[3]*fUpwind_l[13]+0.25*alpha[4]*fUpwind_l[10]+0.25*alpha[2]*fUpwind_l[7]+0.25*fUpwind_l[2]*alpha[7]+0.25*alpha[1]*fUpwind_l[6]+0.25*fUpwind_l[1]*alpha[6]+0.25*alpha[0]*fUpwind_l[3]+0.25*fUpwind_l[0]*alpha[3]; 
  Ghat_l[4] = 0.223606797749979*alpha[16]*fUpwind_l[41]+0.2500000000000001*alpha[22]*fUpwind_l[38]+0.2500000000000001*alpha[21]*fUpwind_l[37]+0.25*alpha[15]*fUpwind_l[31]+0.223606797749979*alpha[9]*fUpwind_l[29]+0.223606797749979*alpha[8]*fUpwind_l[28]+0.2500000000000001*alpha[12]*fUpwind_l[26]+0.2500000000000001*alpha[11]*fUpwind_l[25]+0.2500000000000001*fUpwind_l[11]*alpha[25]+0.25*alpha[7]*fUpwind_l[18]+0.25*alpha[6]*fUpwind_l[17]+0.25*alpha[5]*fUpwind_l[16]+0.25*fUpwind_l[5]*alpha[16]+0.223606797749979*alpha[4]*fUpwind_l[14]+0.25*alpha[3]*fUpwind_l[10]+0.25*alpha[2]*fUpwind_l[9]+0.25*fUpwind_l[2]*alpha[9]+0.25*alpha[1]*fUpwind_l[8]+0.25*fUpwind_l[1]*alpha[8]+0.25*alpha[0]*fUpwind_l[4]+0.25*fUpwind_l[0]*alpha[4]; 
  Ghat_l[5] = 0.223606797749979*alpha[9]*fUpwind_l[36]+0.223606797749979*alpha[8]*fUpwind_l[35]+0.223606797749979*alpha[7]*fUpwind_l[33]+0.223606797749979*alpha[6]*fUpwind_l[32]+0.223606797749979*alpha[16]*fUpwind_l[26]+0.223606797749979*alpha[16]*fUpwind_l[25]+0.223606797749979*fUpwind_l[16]*alpha[25]+0.223606797749979*alpha[15]*fUpwind_l[22]+0.223606797749979*fUpwind_l[15]*alpha[22]+0.223606797749979*alpha[15]*fUpwind_l[21]+0.223606797749979*fUpwind_l[15]*alpha[21]+0.223606797749979*alpha[2]*fUpwind_l[20]+0.223606797749979*alpha[1]*fUpwind_l[19]+0.25*alpha[4]*fUpwind_l[16]+0.25*fUpwind_l[4]*alpha[16]+0.25*alpha[3]*fUpwind_l[15]+0.25*fUpwind_l[3]*alpha[15]+0.223606797749979*alpha[5]*fUpwind_l[12]+0.223606797749979*fUpwind_l[5]*alpha[12]+0.223606797749979*alpha[5]*fUpwind_l[11]+0.223606797749979*fUpwind_l[5]*alpha[11]+0.25*alpha[8]*fUpwind_l[9]+0.25*fUpwind_l[8]*alpha[9]+0.25*alpha[6]*fUpwind_l[7]+0.25*fUpwind_l[6]*alpha[7]+0.25*alpha[0]*fUpwind_l[5]+0.25*fUpwind_l[0]*alpha[5]+0.25*alpha[1]*fUpwind_l[2]+0.25*fUpwind_l[1]*alpha[2]; 
  Ghat_l[6] = 0.223606797749979*alpha[16]*fUpwind_l[44]+0.223606797749979*alpha[8]*fUpwind_l[37]+0.223606797749979*alpha[7]*fUpwind_l[34]+0.25*alpha[12]*fUpwind_l[33]+0.223606797749979*alpha[5]*fUpwind_l[32]+0.25*alpha[9]*fUpwind_l[31]+0.223606797749979*fUpwind_l[17]*alpha[25]+0.223606797749979*alpha[15]*fUpwind_l[24]+0.2*alpha[21]*fUpwind_l[23]+0.223606797749979*alpha[3]*fUpwind_l[23]+0.25*fUpwind_l[20]*alpha[22]+0.223606797749979*alpha[1]*fUpwind_l[21]+0.223606797749979*fUpwind_l[1]*alpha[21]+0.223606797749979*alpha[15]*fUpwind_l[19]+0.25*alpha[16]*fUpwind_l[18]+0.25*alpha[4]*fUpwind_l[17]+0.25*alpha[2]*fUpwind_l[15]+0.25*fUpwind_l[2]*alpha[15]+0.223606797749979*alpha[6]*fUpwind_l[13]+0.223606797749979*alpha[6]*fUpwind_l[11]+0.223606797749979*fUpwind_l[6]*alpha[11]+0.25*alpha[8]*fUpwind_l[10]+0.25*alpha[5]*fUpwind_l[7]+0.25*fUpwind_l[5]*alpha[7]+0.25*alpha[0]*fUpwind_l[6]+0.25*fUpwind_l[0]*alpha[6]+0.25*alpha[1]*fUpwind_l[3]+0.25*fUpwind_l[1]*alpha[3]; 
  Ghat_l[7] = 0.223606797749979*alpha[16]*fUpwind_l[45]+0.25*alpha[25]*fUpwind_l[44]+0.223606797749979*alpha[9]*fUpwind_l[38]+0.223606797749979*alpha[6]*fUpwind_l[34]+0.223606797749979*alpha[5]*fUpwind_l[33]+0.25*alpha[11]*fUpwind_l[32]+0.25*alpha[8]*fUpwind_l[31]+0.2*alpha[22]*fUpwind_l[24]+0.223606797749979*alpha[3]*fUpwind_l[24]+0.223606797749979*alpha[15]*fUpwind_l[23]+0.223606797749979*alpha[2]*fUpwind_l[22]+0.223606797749979*fUpwind_l[2]*alpha[22]+0.25*fUpwind_l[19]*alpha[21]+0.223606797749979*alpha[15]*fUpwind_l[20]+0.25*alpha[4]*fUpwind_l[18]+0.25*alpha[16]*fUpwind_l[17]+0.25*alpha[1]*fUpwind_l[15]+0.25*fUpwind_l[1]*alpha[15]+0.223606797749979*alpha[7]*fUpwind_l[13]+0.223606797749979*alpha[7]*fUpwind_l[12]+0.223606797749979*fUpwind_l[7]*alpha[12]+0.25*alpha[9]*fUpwind_l[10]+0.25*alpha[0]*fUpwind_l[7]+0.25*fUpwind_l[0]*alpha[7]+0.25*alpha[5]*fUpwind_l[6]+0.25*fUpwind_l[5]*alpha[6]+0.25*alpha[2]*fUpwind_l[3]+0.25*fUpwind_l[2]*alpha[3]; 
  Ghat_l[8] = 0.25*alpha[22]*fUpwind_l[45]+0.223606797749979*alpha[15]*fUpwind_l[44]+0.223606797749979*alpha[9]*fUpwind_l[41]+0.223606797749979*alpha[6]*fUpwind_l[37]+0.25*alpha[12]*fUpwind_l[36]+0.223606797749979*alpha[5]*fUpwind_l[35]+0.25*alpha[7]*fUpwind_l[31]+0.223606797749979*alpha[16]*fUpwind_l[29]+0.2*alpha[25]*fUpwind_l[28]+0.223606797749979*alpha[4]*fUpwind_l[28]+0.223606797749979*alpha[1]*fUpwind_l[25]+0.223606797749979*fUpwind_l[1]*alpha[25]+0.223606797749979*fUpwind_l[17]*alpha[21]+0.223606797749979*alpha[16]*fUpwind_l[19]+0.25*alpha[15]*fUpwind_l[18]+0.25*alpha[3]*fUpwind_l[17]+0.25*alpha[2]*fUpwind_l[16]+0.25*fUpwind_l[2]*alpha[16]+0.223606797749979*alpha[8]*fUpwind_l[14]+0.223606797749979*alpha[8]*fUpwind_l[11]+0.223606797749979*fUpwind_l[8]*alpha[11]+0.25*alpha[6]*fUpwind_l[10]+0.25*alpha[5]*fUpwind_l[9]+0.25*fUpwind_l[5]*alpha[9]+0.25*alpha[0]*fUpwind_l[8]+0.25*fUpwind_l[0]*alpha[8]+0.25*alpha[1]*fUpwind_l[4]+0.25*fUpwind_l[1]*alpha[4]; 
  Ghat_l[9] = 0.223606797749979*alpha[15]*fUpwind_l[45]+0.25*alpha[21]*fUpwind_l[44]+0.223606797749979*alpha[8]*fUpwind_l[41]+0.223606797749979*alpha[7]*fUpwind_l[38]+0.223606797749979*alpha[5]*fUpwind_l[36]+0.25*alpha[11]*fUpwind_l[35]+0.25*alpha[6]*fUpwind_l[31]+0.223606797749979*alpha[4]*fUpwind_l[29]+0.223606797749979*alpha[16]*fUpwind_l[28]+0.223606797749979*alpha[2]*fUpwind_l[26]+0.25*fUpwind_l[19]*alpha[25]+0.223606797749979*fUpwind_l[18]*alpha[22]+0.223606797749979*alpha[16]*fUpwind_l[20]+0.25*alpha[3]*fUpwind_l[18]+0.25*alpha[15]*fUpwind_l[17]+0.25*alpha[1]*fUpwind_l[16]+0.25*fUpwind_l[1]*alpha[16]+0.223606797749979*alpha[9]*fUpwind_l[14]+0.223606797749979*alpha[9]*fUpwind_l[12]+0.223606797749979*fUpwind_l[9]*alpha[12]+0.25*alpha[7]*fUpwind_l[10]+0.25*alpha[0]*fUpwind_l[9]+0.25*fUpwind_l[0]*alpha[9]+0.25*alpha[5]*fUpwind_l[8]+0.25*fUpwind_l[5]*alpha[8]+0.25*alpha[2]*fUpwind_l[4]+0.25*fUpwind_l[2]*alpha[4]; 
  Ghat_l[10] = 0.223606797749979*alpha[16]*fUpwind_l[47]+0.223606797749979*alpha[15]*fUpwind_l[46]+0.223606797749979*alpha[9]*fUpwind_l[43]+0.223606797749979*alpha[8]*fUpwind_l[42]+0.223606797749979*alpha[7]*fUpwind_l[40]+0.223606797749979*alpha[6]*fUpwind_l[39]+0.25*alpha[12]*fUpwind_l[38]+0.25*alpha[11]*fUpwind_l[37]+0.25*alpha[5]*fUpwind_l[31]+0.223606797749979*alpha[4]*fUpwind_l[30]+0.223606797749979*alpha[3]*fUpwind_l[27]+0.25*alpha[22]*fUpwind_l[26]+0.25*alpha[21]*fUpwind_l[25]+0.25*fUpwind_l[21]*alpha[25]+0.25*alpha[2]*fUpwind_l[18]+0.25*alpha[1]*fUpwind_l[17]+0.25*alpha[15]*fUpwind_l[16]+0.25*fUpwind_l[15]*alpha[16]+0.25*alpha[0]*fUpwind_l[10]+0.25*alpha[7]*fUpwind_l[9]+0.25*fUpwind_l[7]*alpha[9]+0.25*alpha[6]*fUpwind_l[8]+0.25*fUpwind_l[6]*alpha[8]+0.25*alpha[3]*fUpwind_l[4]+0.25*fUpwind_l[3]*alpha[4]; 
  Ghat_l[11] = 0.25*alpha[9]*fUpwind_l[35]+0.25*alpha[7]*fUpwind_l[32]+0.159719141249985*alpha[25]*fUpwind_l[25]+0.2500000000000001*alpha[4]*fUpwind_l[25]+0.2500000000000001*fUpwind_l[4]*alpha[25]+0.159719141249985*alpha[21]*fUpwind_l[21]+0.2500000000000001*alpha[3]*fUpwind_l[21]+0.2500000000000001*fUpwind_l[3]*alpha[21]+0.2500000000000001*alpha[2]*fUpwind_l[19]+0.223606797749979*alpha[16]*fUpwind_l[16]+0.223606797749979*alpha[15]*fUpwind_l[15]+0.159719141249985*alpha[11]*fUpwind_l[11]+0.25*alpha[0]*fUpwind_l[11]+0.25*fUpwind_l[0]*alpha[11]+0.223606797749979*alpha[8]*fUpwind_l[8]+0.223606797749979*alpha[6]*fUpwind_l[6]+0.223606797749979*alpha[5]*fUpwind_l[5]+0.223606797749979*alpha[1]*fUpwind_l[1]; 
  Ghat_l[12] = 0.25*alpha[8]*fUpwind_l[36]+0.25*alpha[6]*fUpwind_l[33]+0.2500000000000001*alpha[4]*fUpwind_l[26]+0.159719141249985*alpha[22]*fUpwind_l[22]+0.2500000000000001*alpha[3]*fUpwind_l[22]+0.2500000000000001*fUpwind_l[3]*alpha[22]+0.2500000000000001*alpha[1]*fUpwind_l[20]+0.223606797749979*alpha[16]*fUpwind_l[16]+0.223606797749979*alpha[15]*fUpwind_l[15]+0.159719141249985*alpha[12]*fUpwind_l[12]+0.25*alpha[0]*fUpwind_l[12]+0.25*fUpwind_l[0]*alpha[12]+0.223606797749979*alpha[9]*fUpwind_l[9]+0.223606797749979*alpha[7]*fUpwind_l[7]+0.223606797749979*alpha[5]*fUpwind_l[5]+0.223606797749979*alpha[2]*fUpwind_l[2]; 
  Ghat_l[13] = 0.2500000000000001*alpha[16]*fUpwind_l[46]+0.25*alpha[9]*fUpwind_l[40]+0.25*alpha[8]*fUpwind_l[39]+0.25*alpha[5]*fUpwind_l[34]+0.2500000000000001*alpha[4]*fUpwind_l[27]+0.2500000000000001*alpha[2]*fUpwind_l[24]+0.2500000000000001*alpha[1]*fUpwind_l[23]+0.223606797749979*alpha[22]*fUpwind_l[22]+0.223606797749979*alpha[21]*fUpwind_l[21]+0.223606797749979*alpha[15]*fUpwind_l[15]+0.25*alpha[0]*fUpwind_l[13]+0.223606797749979*alpha[7]*fUpwind_l[7]+0.223606797749979*alpha[6]*fUpwind_l[6]+0.223606797749979*alpha[3]*fUpwind_l[3]; 
  Ghat_l[14] = 0.2500000000000001*alpha[15]*fUpwind_l[47]+0.25*alpha[7]*fUpwind_l[43]+0.25*alpha[6]*fUpwind_l[42]+0.25*alpha[5]*fUpwind_l[41]+0.2500000000000001*alpha[3]*fUpwind_l[30]+0.2500000000000001*alpha[2]*fUpwind_l[29]+0.2500000000000001*alpha[1]*fUpwind_l[28]+0.223606797749979*alpha[25]*fUpwind_l[25]+0.223606797749979*alpha[16]*fUpwind_l[16]+0.25*alpha[0]*fUpwind_l[14]+0.223606797749979*alpha[9]*fUpwind_l[9]+0.223606797749979*alpha[8]*fUpwind_l[8]+0.223606797749979*alpha[4]*fUpwind_l[4]; 
  Ghat_l[15] = 0.223606797749979*alpha[9]*fUpwind_l[45]+0.223606797749979*alpha[8]*fUpwind_l[44]+0.223606797749979*alpha[16]*fUpwind_l[38]+0.223606797749979*alpha[16]*fUpwind_l[37]+0.2*alpha[22]*fUpwind_l[34]+0.2*alpha[21]*fUpwind_l[34]+0.223606797749979*alpha[3]*fUpwind_l[34]+0.223606797749979*alpha[2]*fUpwind_l[33]+0.223606797749979*alpha[1]*fUpwind_l[32]+0.223606797749979*alpha[25]*fUpwind_l[31]+0.25*alpha[4]*fUpwind_l[31]+0.223606797749979*alpha[6]*fUpwind_l[24]+0.223606797749979*alpha[7]*fUpwind_l[23]+0.223606797749979*alpha[5]*fUpwind_l[22]+0.223606797749979*fUpwind_l[5]*alpha[22]+0.223606797749979*alpha[5]*fUpwind_l[21]+0.223606797749979*fUpwind_l[5]*alpha[21]+0.223606797749979*alpha[7]*fUpwind_l[20]+0.223606797749979*alpha[6]*fUpwind_l[19]+0.25*alpha[8]*fUpwind_l[18]+0.25*alpha[9]*fUpwind_l[17]+0.25*fUpwind_l[10]*alpha[16]+0.223606797749979*alpha[12]*fUpwind_l[15]+0.223606797749979*alpha[11]*fUpwind_l[15]+0.25*alpha[0]*fUpwind_l[15]+0.223606797749979*fUpwind_l[13]*alpha[15]+0.223606797749979*fUpwind_l[12]*alpha[15]+0.223606797749979*fUpwind_l[11]*alpha[15]+0.25*fUpwind_l[0]*alpha[15]+0.25*alpha[1]*fUpwind_l[7]+0.25*fUpwind_l[1]*alpha[7]+0.25*alpha[2]*fUpwind_l[6]+0.25*fUpwind_l[2]*alpha[6]+0.25*alpha[3]*fUpwind_l[5]+0.25*fUpwind_l[3]*alpha[5]; 
  Ghat_l[16] = 0.223606797749979*alpha[7]*fUpwind_l[45]+0.223606797749979*alpha[6]*fUpwind_l[44]+0.2*alpha[25]*fUpwind_l[41]+0.223606797749979*alpha[4]*fUpwind_l[41]+0.223606797749979*alpha[15]*fUpwind_l[38]+0.223606797749979*alpha[15]*fUpwind_l[37]+0.223606797749979*alpha[2]*fUpwind_l[36]+0.223606797749979*alpha[1]*fUpwind_l[35]+0.223606797749979*alpha[22]*fUpwind_l[31]+0.223606797749979*alpha[21]*fUpwind_l[31]+0.25*alpha[3]*fUpwind_l[31]+0.223606797749979*alpha[8]*fUpwind_l[29]+0.223606797749979*alpha[9]*fUpwind_l[28]+0.223606797749979*alpha[5]*fUpwind_l[26]+0.223606797749979*alpha[5]*fUpwind_l[25]+0.223606797749979*fUpwind_l[5]*alpha[25]+0.223606797749979*alpha[9]*fUpwind_l[20]+0.223606797749979*alpha[8]*fUpwind_l[19]+0.25*alpha[6]*fUpwind_l[18]+0.25*alpha[7]*fUpwind_l[17]+0.223606797749979*alpha[12]*fUpwind_l[16]+0.223606797749979*alpha[11]*fUpwind_l[16]+0.25*alpha[0]*fUpwind_l[16]+0.223606797749979*fUpwind_l[14]*alpha[16]+0.223606797749979*fUpwind_l[12]*alpha[16]+0.223606797749979*fUpwind_l[11]*alpha[16]+0.25*fUpwind_l[0]*alpha[16]+0.25*fUpwind_l[10]*alpha[15]+0.25*alpha[1]*fUpwind_l[9]+0.25*fUpwind_l[1]*alpha[9]+0.25*alpha[2]*fUpwind_l[8]+0.25*fUpwind_l[2]*alpha[8]+0.25*alpha[4]*fUpwind_l[5]+0.25*fUpwind_l[4]*alpha[5]; 
  Ghat_l[17] = 0.223606797749979*alpha[9]*fUpwind_l[47]+0.223606797749979*alpha[7]*fUpwind_l[46]+0.2500000000000001*alpha[12]*fUpwind_l[45]+0.223606797749979*alpha[5]*fUpwind_l[44]+0.223606797749979*alpha[16]*fUpwind_l[43]+0.2*alpha[25]*fUpwind_l[42]+0.223606797749979*alpha[4]*fUpwind_l[42]+0.223606797749979*alpha[15]*fUpwind_l[40]+0.2*alpha[21]*fUpwind_l[39]+0.223606797749979*alpha[3]*fUpwind_l[39]+0.223606797749979*alpha[1]*fUpwind_l[37]+0.2500000000000001*alpha[22]*fUpwind_l[36]+0.223606797749979*alpha[15]*fUpwind_l[35]+0.223606797749979*alpha[16]*fUpwind_l[32]+0.25*alpha[2]*fUpwind_l[31]+0.223606797749979*alpha[8]*fUpwind_l[30]+0.223606797749979*alpha[6]*fUpwind_l[27]+0.223606797749979*alpha[6]*fUpwind_l[25]+0.223606797749979*fUpwind_l[6]*alpha[25]+0.223606797749979*alpha[8]*fUpwind_l[21]+0.223606797749979*fUpwind_l[8]*alpha[21]+0.25*alpha[5]*fUpwind_l[18]+0.223606797749979*alpha[11]*fUpwind_l[17]+0.25*alpha[0]*fUpwind_l[17]+0.25*alpha[7]*fUpwind_l[16]+0.25*fUpwind_l[7]*alpha[16]+0.25*alpha[9]*fUpwind_l[15]+0.25*fUpwind_l[9]*alpha[15]+0.25*alpha[1]*fUpwind_l[10]+0.25*alpha[3]*fUpwind_l[8]+0.25*fUpwind_l[3]*alpha[8]+0.25*alpha[4]*fUpwind_l[6]+0.25*fUpwind_l[4]*alpha[6]; 
  Ghat_l[18] = 0.223606797749979*alpha[8]*fUpwind_l[47]+0.223606797749979*alpha[6]*fUpwind_l[46]+0.223606797749979*alpha[5]*fUpwind_l[45]+0.2500000000000001*alpha[11]*fUpwind_l[44]+0.223606797749979*alpha[4]*fUpwind_l[43]+0.223606797749979*alpha[16]*fUpwind_l[42]+0.2*alpha[22]*fUpwind_l[40]+0.223606797749979*alpha[3]*fUpwind_l[40]+0.223606797749979*alpha[15]*fUpwind_l[39]+0.223606797749979*alpha[2]*fUpwind_l[38]+0.223606797749979*alpha[15]*fUpwind_l[36]+0.2500000000000001*alpha[21]*fUpwind_l[35]+0.223606797749979*alpha[16]*fUpwind_l[33]+0.2500000000000001*alpha[25]*fUpwind_l[32]+0.25*alpha[1]*fUpwind_l[31]+0.223606797749979*alpha[9]*fUpwind_l[30]+0.223606797749979*alpha[7]*fUpwind_l[27]+0.223606797749979*alpha[7]*fUpwind_l[26]+0.223606797749979*alpha[9]*fUpwind_l[22]+0.223606797749979*fUpwind_l[9]*alpha[22]+0.223606797749979*alpha[12]*fUpwind_l[18]+0.25*alpha[0]*fUpwind_l[18]+0.25*alpha[5]*fUpwind_l[17]+0.25*alpha[6]*fUpwind_l[16]+0.25*fUpwind_l[6]*alpha[16]+0.25*alpha[8]*fUpwind_l[15]+0.25*fUpwind_l[8]*alpha[15]+0.25*alpha[2]*fUpwind_l[10]+0.25*alpha[3]*fUpwind_l[9]+0.25*fUpwind_l[3]*alpha[9]+0.25*alpha[4]*fUpwind_l[7]+0.25*fUpwind_l[4]*alpha[7]; 
  Ghat_l[19] = 0.2*alpha[16]*fUpwind_l[36]+0.159719141249985*alpha[25]*fUpwind_l[35]+0.2500000000000001*alpha[4]*fUpwind_l[35]+0.2*alpha[15]*fUpwind_l[33]+0.223606797749979*alpha[22]*fUpwind_l[32]+0.159719141249985*alpha[21]*fUpwind_l[32]+0.2500000000000001*alpha[3]*fUpwind_l[32]+0.25*alpha[9]*fUpwind_l[25]+0.25*fUpwind_l[9]*alpha[25]+0.25*alpha[7]*fUpwind_l[21]+0.25*fUpwind_l[7]*alpha[21]+0.2*alpha[5]*fUpwind_l[20]+0.223606797749979*alpha[12]*fUpwind_l[19]+0.159719141249985*alpha[11]*fUpwind_l[19]+0.25*alpha[0]*fUpwind_l[19]+0.223606797749979*alpha[8]*fUpwind_l[16]+0.223606797749979*fUpwind_l[8]*alpha[16]+0.223606797749979*alpha[6]*fUpwind_l[15]+0.223606797749979*fUpwind_l[6]*alpha[15]+0.2500000000000001*alpha[2]*fUpwind_l[11]+0.2500000000000001*fUpwind_l[2]*alpha[11]+0.223606797749979*alpha[1]*fUpwind_l[5]+0.223606797749979*fUpwind_l[1]*alpha[5]; 
  Ghat_l[20] = 0.223606797749979*alpha[25]*fUpwind_l[36]+0.2500000000000001*alpha[4]*fUpwind_l[36]+0.2*alpha[16]*fUpwind_l[35]+0.159719141249985*alpha[22]*fUpwind_l[33]+0.223606797749979*alpha[21]*fUpwind_l[33]+0.2500000000000001*alpha[3]*fUpwind_l[33]+0.2*alpha[15]*fUpwind_l[32]+0.25*alpha[8]*fUpwind_l[26]+0.25*alpha[6]*fUpwind_l[22]+0.25*fUpwind_l[6]*alpha[22]+0.159719141249985*alpha[12]*fUpwind_l[20]+0.223606797749979*alpha[11]*fUpwind_l[20]+0.25*alpha[0]*fUpwind_l[20]+0.2*alpha[5]*fUpwind_l[19]+0.223606797749979*alpha[9]*fUpwind_l[16]+0.223606797749979*fUpwind_l[9]*alpha[16]+0.223606797749979*alpha[7]*fUpwind_l[15]+0.223606797749979*fUpwind_l[7]*alpha[15]+0.2500000000000001*alpha[1]*fUpwind_l[12]+0.2500000000000001*fUpwind_l[1]*alpha[12]+0.223606797749979*alpha[2]*fUpwind_l[5]+0.223606797749979*fUpwind_l[2]*alpha[5]; 
  Ghat_l[21] = 0.25*alpha[9]*fUpwind_l[44]+0.159719141249985*alpha[25]*fUpwind_l[37]+0.2500000000000001*alpha[4]*fUpwind_l[37]+0.2*alpha[15]*fUpwind_l[34]+0.2500000000000001*alpha[2]*fUpwind_l[32]+0.223606797749979*alpha[16]*fUpwind_l[31]+0.25*fUpwind_l[10]*alpha[25]+0.2*alpha[6]*fUpwind_l[23]+0.159719141249985*alpha[11]*fUpwind_l[21]+0.25*alpha[0]*fUpwind_l[21]+0.223606797749979*fUpwind_l[13]*alpha[21]+0.159719141249985*fUpwind_l[11]*alpha[21]+0.25*fUpwind_l[0]*alpha[21]+0.25*alpha[7]*fUpwind_l[19]+0.223606797749979*alpha[8]*fUpwind_l[17]+0.223606797749979*alpha[5]*fUpwind_l[15]+0.223606797749979*fUpwind_l[5]*alpha[15]+0.2500000000000001*alpha[3]*fUpwind_l[11]+0.2500000000000001*fUpwind_l[3]*alpha[11]+0.223606797749979*alpha[1]*fUpwind_l[6]+0.223606797749979*fUpwind_l[1]*alpha[6]; 
  Ghat_l[22] = 0.25*alpha[8]*fUpwind_l[45]+0.2500000000000001*alpha[4]*fUpwind_l[38]+0.2*alpha[15]*fUpwind_l[34]+0.2500000000000001*alpha[1]*fUpwind_l[33]+0.223606797749979*alpha[16]*fUpwind_l[31]+0.2*alpha[7]*fUpwind_l[24]+0.159719141249985*alpha[12]*fUpwind_l[22]+0.25*alpha[0]*fUpwind_l[22]+0.223606797749979*fUpwind_l[13]*alpha[22]+0.159719141249985*fUpwind_l[12]*alpha[22]+0.25*fUpwind_l[0]*alpha[22]+0.25*alpha[6]*fUpwind_l[20]+0.223606797749979*alpha[9]*fUpwind_l[18]+0.223606797749979*alpha[5]*fUpwind_l[15]+0.223606797749979*fUpwind_l[5]*alpha[15]+0.2500000000000001*alpha[3]*fUpwind_l[12]+0.2500000000000001*fUpwind_l[3]*alpha[12]+0.223606797749979*alpha[2]*fUpwind_l[7]+0.223606797749979*fUpwind_l[2]*alpha[7]; 
  Ghat_l[23] = 0.25*alpha[9]*fUpwind_l[46]+0.2500000000000001*alpha[16]*fUpwind_l[40]+0.223606797749979*alpha[25]*fUpwind_l[39]+0.2500000000000001*alpha[4]*fUpwind_l[39]+0.2500000000000001*alpha[2]*fUpwind_l[34]+0.223606797749979*alpha[22]*fUpwind_l[33]+0.2*alpha[15]*fUpwind_l[32]+0.25*alpha[8]*fUpwind_l[27]+0.25*alpha[5]*fUpwind_l[24]+0.223606797749979*alpha[11]*fUpwind_l[23]+0.25*alpha[0]*fUpwind_l[23]+0.2*alpha[6]*fUpwind_l[21]+0.2*fUpwind_l[6]*alpha[21]+0.223606797749979*alpha[7]*fUpwind_l[15]+0.223606797749979*fUpwind_l[7]*alpha[15]+0.2500000000000001*alpha[1]*fUpwind_l[13]+0.223606797749979*alpha[3]*fUpwind_l[6]+0.223606797749979*fUpwind_l[3]*alpha[6]; 
  Ghat_l[24] = 0.25*alpha[8]*fUpwind_l[46]+0.2500000000000001*alpha[4]*fUpwind_l[40]+0.2500000000000001*alpha[16]*fUpwind_l[39]+0.2500000000000001*alpha[1]*fUpwind_l[34]+0.2*alpha[15]*fUpwind_l[33]+0.223606797749979*alpha[21]*fUpwind_l[32]+0.25*alpha[9]*fUpwind_l[27]+0.223606797749979*alpha[12]*fUpwind_l[24]+0.25*alpha[0]*fUpwind_l[24]+0.25*alpha[5]*fUpwind_l[23]+0.2*alpha[7]*fUpwind_l[22]+0.2*fUpwind_l[7]*alpha[22]+0.223606797749979*alpha[6]*fUpwind_l[15]+0.223606797749979*fUpwind_l[6]*alpha[15]+0.2500000000000001*alpha[2]*fUpwind_l[13]+0.223606797749979*alpha[3]*fUpwind_l[7]+0.223606797749979*fUpwind_l[3]*alpha[7]; 
  Ghat_l[25] = 0.25*alpha[7]*fUpwind_l[44]+0.2*alpha[16]*fUpwind_l[41]+0.159719141249985*alpha[21]*fUpwind_l[37]+0.2500000000000001*alpha[3]*fUpwind_l[37]+0.2500000000000001*alpha[2]*fUpwind_l[35]+0.223606797749979*alpha[15]*fUpwind_l[31]+0.2*alpha[8]*fUpwind_l[28]+0.159719141249985*alpha[11]*fUpwind_l[25]+0.25*alpha[0]*fUpwind_l[25]+0.223606797749979*fUpwind_l[14]*alpha[25]+0.159719141249985*fUpwind_l[11]*alpha[25]+0.25*fUpwind_l[0]*alpha[25]+0.25*fUpwind_l[10]*alpha[21]+0.25*alpha[9]*fUpwind_l[19]+0.223606797749979*alpha[6]*fUpwind_l[17]+0.223606797749979*alpha[5]*fUpwind_l[16]+0.223606797749979*fUpwind_l[5]*alpha[16]+0.2500000000000001*alpha[4]*fUpwind_l[11]+0.2500000000000001*fUpwind_l[4]*alpha[11]+0.223606797749979*alpha[1]*fUpwind_l[8]+0.223606797749979*fUpwind_l[1]*alpha[8]; 
  Ghat_l[26] = 0.25*alpha[6]*fUpwind_l[45]+0.2*alpha[16]*fUpwind_l[41]+0.159719141249985*alpha[22]*fUpwind_l[38]+0.2500000000000001*alpha[3]*fUpwind_l[38]+0.2500000000000001*alpha[1]*fUpwind_l[36]+0.223606797749979*alpha[15]*fUpwind_l[31]+0.2*alpha[9]*fUpwind_l[29]+0.159719141249985*alpha[12]*fUpwind_l[26]+0.25*alpha[0]*fUpwind_l[26]+0.25*fUpwind_l[10]*alpha[22]+0.25*alpha[8]*fUpwind_l[20]+0.223606797749979*alpha[7]*fUpwind_l[18]+0.223606797749979*alpha[5]*fUpwind_l[16]+0.223606797749979*fUpwind_l[5]*alpha[16]+0.2500000000000001*alpha[4]*fUpwind_l[12]+0.2500000000000001*fUpwind_l[4]*alpha[12]+0.223606797749979*alpha[2]*fUpwind_l[9]+0.223606797749979*fUpwind_l[2]*alpha[9]; 
  Ghat_l[27] = 0.25*alpha[5]*fUpwind_l[46]+0.2500000000000001*alpha[2]*fUpwind_l[40]+0.2500000000000001*alpha[1]*fUpwind_l[39]+0.223606797749979*alpha[22]*fUpwind_l[38]+0.223606797749979*alpha[21]*fUpwind_l[37]+0.2500000000000001*alpha[16]*fUpwind_l[34]+0.223606797749979*alpha[15]*fUpwind_l[31]+0.25*alpha[0]*fUpwind_l[27]+0.25*alpha[9]*fUpwind_l[24]+0.25*alpha[8]*fUpwind_l[23]+0.223606797749979*alpha[7]*fUpwind_l[18]+0.223606797749979*alpha[6]*fUpwind_l[17]+0.2500000000000001*alpha[4]*fUpwind_l[13]+0.223606797749979*alpha[3]*fUpwind_l[10]; 
  Ghat_l[28] = 0.25*alpha[7]*fUpwind_l[47]+0.2500000000000001*alpha[15]*fUpwind_l[43]+0.223606797749979*alpha[21]*fUpwind_l[42]+0.2500000000000001*alpha[3]*fUpwind_l[42]+0.2500000000000001*alpha[2]*fUpwind_l[41]+0.2*alpha[16]*fUpwind_l[35]+0.25*alpha[6]*fUpwind_l[30]+0.25*alpha[5]*fUpwind_l[29]+0.223606797749979*alpha[11]*fUpwind_l[28]+0.25*alpha[0]*fUpwind_l[28]+0.2*alpha[8]*fUpwind_l[25]+0.2*fUpwind_l[8]*alpha[25]+0.223606797749979*alpha[9]*fUpwind_l[16]+0.223606797749979*fUpwind_l[9]*alpha[16]+0.2500000000000001*alpha[1]*fUpwind_l[14]+0.223606797749979*alpha[4]*fUpwind_l[8]+0.223606797749979*fUpwind_l[4]*alpha[8]; 
  Ghat_l[29] = 0.25*alpha[6]*fUpwind_l[47]+0.223606797749979*alpha[22]*fUpwind_l[43]+0.2500000000000001*alpha[3]*fUpwind_l[43]+0.2500000000000001*alpha[15]*fUpwind_l[42]+0.2500000000000001*alpha[1]*fUpwind_l[41]+0.2*alpha[16]*fUpwind_l[36]+0.223606797749979*alpha[25]*fUpwind_l[35]+0.25*alpha[7]*fUpwind_l[30]+0.223606797749979*alpha[12]*fUpwind_l[29]+0.25*alpha[0]*fUpwind_l[29]+0.25*alpha[5]*fUpwind_l[28]+0.2*alpha[9]*fUpwind_l[26]+0.223606797749979*alpha[8]*fUpwind_l[16]+0.223606797749979*fUpwind_l[8]*alpha[16]+0.2500000000000001*alpha[2]*fUpwind_l[14]+0.223606797749979*alpha[4]*fUpwind_l[9]+0.223606797749979*fUpwind_l[4]*alpha[9]; 
  Ghat_l[30] = 0.25*alpha[5]*fUpwind_l[47]+0.2500000000000001*alpha[2]*fUpwind_l[43]+0.2500000000000001*alpha[1]*fUpwind_l[42]+0.2500000000000001*alpha[15]*fUpwind_l[41]+0.223606797749979*alpha[25]*fUpwind_l[37]+0.223606797749979*alpha[16]*fUpwind_l[31]+0.25*alpha[0]*fUpwind_l[30]+0.25*alpha[7]*fUpwind_l[29]+0.25*alpha[6]*fUpwind_l[28]+0.223606797749979*alpha[9]*fUpwind_l[18]+0.223606797749979*alpha[8]*fUpwind_l[17]+0.2500000000000001*alpha[3]*fUpwind_l[14]+0.223606797749979*alpha[4]*fUpwind_l[10]; 
  Ghat_l[31] = 0.2*alpha[25]*fUpwind_l[47]+0.223606797749979*alpha[4]*fUpwind_l[47]+0.2*alpha[22]*fUpwind_l[46]+0.2*alpha[21]*fUpwind_l[46]+0.223606797749979*alpha[3]*fUpwind_l[46]+0.223606797749979*alpha[2]*fUpwind_l[45]+0.223606797749979*alpha[1]*fUpwind_l[44]+0.223606797749979*alpha[8]*fUpwind_l[43]+0.223606797749979*alpha[9]*fUpwind_l[42]+0.223606797749979*alpha[6]*fUpwind_l[40]+0.223606797749979*alpha[7]*fUpwind_l[39]+0.223606797749979*alpha[5]*fUpwind_l[38]+0.223606797749979*alpha[5]*fUpwind_l[37]+0.223606797749979*alpha[7]*fUpwind_l[36]+0.223606797749979*alpha[6]*fUpwind_l[35]+0.223606797749979*alpha[9]*fUpwind_l[33]+0.223606797749979*alpha[8]*fUpwind_l[32]+0.223606797749979*alpha[12]*fUpwind_l[31]+0.223606797749979*alpha[11]*fUpwind_l[31]+0.25*alpha[0]*fUpwind_l[31]+0.223606797749979*alpha[16]*fUpwind_l[30]+0.223606797749979*alpha[15]*fUpwind_l[27]+0.223606797749979*alpha[15]*fUpwind_l[26]+0.223606797749979*alpha[15]*fUpwind_l[25]+0.223606797749979*fUpwind_l[15]*alpha[25]+0.223606797749979*alpha[16]*fUpwind_l[22]+0.223606797749979*fUpwind_l[16]*alpha[22]+0.223606797749979*alpha[16]*fUpwind_l[21]+0.223606797749979*fUpwind_l[16]*alpha[21]+0.25*alpha[1]*fUpwind_l[18]+0.25*alpha[2]*fUpwind_l[17]+0.25*alpha[3]*fUpwind_l[16]+0.25*fUpwind_l[3]*alpha[16]+0.25*alpha[4]*fUpwind_l[15]+0.25*fUpwind_l[4]*alpha[15]+0.25*alpha[5]*fUpwind_l[10]+0.25*alpha[6]*fUpwind_l[9]+0.25*fUpwind_l[6]*alpha[9]+0.25*alpha[7]*fUpwind_l[8]+0.25*fUpwind_l[7]*alpha[8]; 
  Ghat_l[32] = 0.2*alpha[16]*fUpwind_l[45]+0.159719141249985*alpha[25]*fUpwind_l[44]+0.2500000000000001*alpha[4]*fUpwind_l[44]+0.25*alpha[9]*fUpwind_l[37]+0.2*alpha[6]*fUpwind_l[34]+0.2*alpha[5]*fUpwind_l[33]+0.223606797749979*alpha[12]*fUpwind_l[32]+0.159719141249985*alpha[11]*fUpwind_l[32]+0.25*alpha[0]*fUpwind_l[32]+0.223606797749979*alpha[8]*fUpwind_l[31]+0.2500000000000001*fUpwind_l[18]*alpha[25]+0.223606797749979*alpha[21]*fUpwind_l[24]+0.2*alpha[15]*fUpwind_l[23]+0.223606797749979*fUpwind_l[19]*alpha[22]+0.2500000000000001*alpha[2]*fUpwind_l[21]+0.159719141249985*fUpwind_l[19]*alpha[21]+0.2500000000000001*fUpwind_l[2]*alpha[21]+0.2*alpha[15]*fUpwind_l[20]+0.2500000000000001*alpha[3]*fUpwind_l[19]+0.223606797749979*alpha[16]*fUpwind_l[17]+0.223606797749979*alpha[1]*fUpwind_l[15]+0.223606797749979*fUpwind_l[1]*alpha[15]+0.25*alpha[7]*fUpwind_l[11]+0.25*fUpwind_l[7]*alpha[11]+0.223606797749979*alpha[5]*fUpwind_l[6]+0.223606797749979*fUpwind_l[5]*alpha[6]; 
  Ghat_l[33] = 0.223606797749979*alpha[25]*fUpwind_l[45]+0.2500000000000001*alpha[4]*fUpwind_l[45]+0.2*alpha[16]*fUpwind_l[44]+0.25*alpha[8]*fUpwind_l[38]+0.2*alpha[7]*fUpwind_l[34]+0.159719141249985*alpha[12]*fUpwind_l[33]+0.223606797749979*alpha[11]*fUpwind_l[33]+0.25*alpha[0]*fUpwind_l[33]+0.2*alpha[5]*fUpwind_l[32]+0.223606797749979*alpha[9]*fUpwind_l[31]+0.2*alpha[15]*fUpwind_l[24]+0.223606797749979*alpha[22]*fUpwind_l[23]+0.2500000000000001*alpha[1]*fUpwind_l[22]+0.159719141249985*fUpwind_l[20]*alpha[22]+0.2500000000000001*fUpwind_l[1]*alpha[22]+0.223606797749979*fUpwind_l[20]*alpha[21]+0.2500000000000001*alpha[3]*fUpwind_l[20]+0.2*alpha[15]*fUpwind_l[19]+0.223606797749979*alpha[16]*fUpwind_l[18]+0.223606797749979*alpha[2]*fUpwind_l[15]+0.223606797749979*fUpwind_l[2]*alpha[15]+0.25*alpha[6]*fUpwind_l[12]+0.25*fUpwind_l[6]*alpha[12]+0.223606797749979*alpha[5]*fUpwind_l[7]+0.223606797749979*fUpwind_l[5]*alpha[7]; 
  Ghat_l[34] = 0.223606797749979*alpha[25]*fUpwind_l[46]+0.2500000000000001*alpha[4]*fUpwind_l[46]+0.25*alpha[8]*fUpwind_l[40]+0.25*alpha[9]*fUpwind_l[39]+0.223606797749979*alpha[12]*fUpwind_l[34]+0.223606797749979*alpha[11]*fUpwind_l[34]+0.25*alpha[0]*fUpwind_l[34]+0.2*alpha[7]*fUpwind_l[33]+0.2*alpha[6]*fUpwind_l[32]+0.2500000000000001*alpha[16]*fUpwind_l[27]+0.2500000000000001*alpha[1]*fUpwind_l[24]+0.2500000000000001*alpha[2]*fUpwind_l[23]+0.2*alpha[15]*fUpwind_l[22]+0.2*fUpwind_l[15]*alpha[22]+0.2*alpha[15]*fUpwind_l[21]+0.2*fUpwind_l[15]*alpha[21]+0.223606797749979*alpha[3]*fUpwind_l[15]+0.223606797749979*fUpwind_l[3]*alpha[15]+0.25*alpha[5]*fUpwind_l[13]+0.223606797749979*alpha[6]*fUpwind_l[7]+0.223606797749979*fUpwind_l[6]*alpha[7]; 
  Ghat_l[35] = 0.2*alpha[15]*fUpwind_l[45]+0.223606797749979*alpha[22]*fUpwind_l[44]+0.159719141249985*alpha[21]*fUpwind_l[44]+0.2500000000000001*alpha[3]*fUpwind_l[44]+0.2*alpha[8]*fUpwind_l[41]+0.25*alpha[7]*fUpwind_l[37]+0.2*alpha[5]*fUpwind_l[36]+0.223606797749979*alpha[12]*fUpwind_l[35]+0.159719141249985*alpha[11]*fUpwind_l[35]+0.25*alpha[0]*fUpwind_l[35]+0.223606797749979*alpha[6]*fUpwind_l[31]+0.223606797749979*alpha[25]*fUpwind_l[29]+0.2*alpha[16]*fUpwind_l[28]+0.2500000000000001*alpha[2]*fUpwind_l[25]+0.159719141249985*fUpwind_l[19]*alpha[25]+0.2500000000000001*fUpwind_l[2]*alpha[25]+0.2500000000000001*fUpwind_l[18]*alpha[21]+0.2*alpha[16]*fUpwind_l[20]+0.2500000000000001*alpha[4]*fUpwind_l[19]+0.223606797749979*alpha[15]*fUpwind_l[17]+0.223606797749979*alpha[1]*fUpwind_l[16]+0.223606797749979*fUpwind_l[1]*alpha[16]+0.25*alpha[9]*fUpwind_l[11]+0.25*fUpwind_l[9]*alpha[11]+0.223606797749979*alpha[5]*fUpwind_l[8]+0.223606797749979*fUpwind_l[5]*alpha[8]; 
  Ghat_l[36] = 0.159719141249985*alpha[22]*fUpwind_l[45]+0.223606797749979*alpha[21]*fUpwind_l[45]+0.2500000000000001*alpha[3]*fUpwind_l[45]+0.2*alpha[15]*fUpwind_l[44]+0.2*alpha[9]*fUpwind_l[41]+0.25*alpha[6]*fUpwind_l[38]+0.159719141249985*alpha[12]*fUpwind_l[36]+0.223606797749979*alpha[11]*fUpwind_l[36]+0.25*alpha[0]*fUpwind_l[36]+0.2*alpha[5]*fUpwind_l[35]+0.223606797749979*alpha[7]*fUpwind_l[31]+0.2*alpha[16]*fUpwind_l[29]+0.2500000000000001*alpha[1]*fUpwind_l[26]+0.223606797749979*fUpwind_l[20]*alpha[25]+0.2500000000000001*fUpwind_l[17]*alpha[22]+0.2500000000000001*alpha[4]*fUpwind_l[20]+0.2*alpha[16]*fUpwind_l[19]+0.223606797749979*alpha[15]*fUpwind_l[18]+0.223606797749979*alpha[2]*fUpwind_l[16]+0.223606797749979*fUpwind_l[2]*alpha[16]+0.25*alpha[8]*fUpwind_l[12]+0.25*fUpwind_l[8]*alpha[12]+0.223606797749979*alpha[5]*fUpwind_l[9]+0.223606797749979*fUpwind_l[5]*alpha[9]; 
  Ghat_l[37] = 0.2*alpha[16]*fUpwind_l[47]+0.2*alpha[15]*fUpwind_l[46]+0.2500000000000001*alpha[2]*fUpwind_l[44]+0.2*alpha[8]*fUpwind_l[42]+0.2*alpha[6]*fUpwind_l[39]+0.159719141249985*alpha[11]*fUpwind_l[37]+0.25*alpha[0]*fUpwind_l[37]+0.25*alpha[7]*fUpwind_l[35]+0.25*alpha[9]*fUpwind_l[32]+0.223606797749979*alpha[5]*fUpwind_l[31]+0.223606797749979*alpha[25]*fUpwind_l[30]+0.223606797749979*alpha[21]*fUpwind_l[27]+0.159719141249985*alpha[21]*fUpwind_l[25]+0.2500000000000001*alpha[3]*fUpwind_l[25]+0.159719141249985*fUpwind_l[21]*alpha[25]+0.2500000000000001*fUpwind_l[3]*alpha[25]+0.2500000000000001*alpha[4]*fUpwind_l[21]+0.2500000000000001*fUpwind_l[4]*alpha[21]+0.223606797749979*alpha[1]*fUpwind_l[17]+0.223606797749979*alpha[15]*fUpwind_l[16]+0.223606797749979*fUpwind_l[15]*alpha[16]+0.25*fUpwind_l[10]*alpha[11]+0.223606797749979*alpha[6]*fUpwind_l[8]+0.223606797749979*fUpwind_l[6]*alpha[8]; 
  Ghat_l[38] = 0.2*alpha[16]*fUpwind_l[47]+0.2*alpha[15]*fUpwind_l[46]+0.2500000000000001*alpha[1]*fUpwind_l[45]+0.2*alpha[9]*fUpwind_l[43]+0.2*alpha[7]*fUpwind_l[40]+0.159719141249985*alpha[12]*fUpwind_l[38]+0.25*alpha[0]*fUpwind_l[38]+0.25*alpha[6]*fUpwind_l[36]+0.25*alpha[8]*fUpwind_l[33]+0.223606797749979*alpha[5]*fUpwind_l[31]+0.223606797749979*alpha[22]*fUpwind_l[27]+0.159719141249985*alpha[22]*fUpwind_l[26]+0.2500000000000001*alpha[3]*fUpwind_l[26]+0.2500000000000001*alpha[4]*fUpwind_l[22]+0.2500000000000001*fUpwind_l[4]*alpha[22]+0.223606797749979*alpha[2]*fUpwind_l[18]+0.223606797749979*alpha[15]*fUpwind_l[16]+0.223606797749979*fUpwind_l[15]*alpha[16]+0.25*fUpwind_l[10]*alpha[12]+0.223606797749979*alpha[7]*fUpwind_l[9]+0.223606797749979*fUpwind_l[7]*alpha[9]; 
  Ghat_l[39] = 0.2500000000000001*alpha[2]*fUpwind_l[46]+0.223606797749979*alpha[22]*fUpwind_l[45]+0.2*alpha[15]*fUpwind_l[44]+0.25*alpha[5]*fUpwind_l[40]+0.223606797749979*alpha[11]*fUpwind_l[39]+0.25*alpha[0]*fUpwind_l[39]+0.2*alpha[6]*fUpwind_l[37]+0.25*alpha[9]*fUpwind_l[34]+0.223606797749979*alpha[7]*fUpwind_l[31]+0.2500000000000001*alpha[1]*fUpwind_l[27]+0.223606797749979*fUpwind_l[23]*alpha[25]+0.2500000000000001*alpha[16]*fUpwind_l[24]+0.2500000000000001*alpha[4]*fUpwind_l[23]+0.2*fUpwind_l[17]*alpha[21]+0.223606797749979*alpha[15]*fUpwind_l[18]+0.223606797749979*alpha[3]*fUpwind_l[17]+0.25*alpha[8]*fUpwind_l[13]+0.223606797749979*alpha[6]*fUpwind_l[10]; 
  Ghat_l[40] = 0.2500000000000001*alpha[1]*fUpwind_l[46]+0.2*alpha[15]*fUpwind_l[45]+0.223606797749979*alpha[21]*fUpwind_l[44]+0.223606797749979*alpha[12]*fUpwind_l[40]+0.25*alpha[0]*fUpwind_l[40]+0.25*alpha[5]*fUpwind_l[39]+0.2*alpha[7]*fUpwind_l[38]+0.25*alpha[8]*fUpwind_l[34]+0.223606797749979*alpha[6]*fUpwind_l[31]+0.2500000000000001*alpha[2]*fUpwind_l[27]+0.2500000000000001*alpha[4]*fUpwind_l[24]+0.2500000000000001*alpha[16]*fUpwind_l[23]+0.2*fUpwind_l[18]*alpha[22]+0.223606797749979*alpha[3]*fUpwind_l[18]+0.223606797749979*alpha[15]*fUpwind_l[17]+0.25*alpha[9]*fUpwind_l[13]+0.223606797749979*alpha[7]*fUpwind_l[10]; 
  Ghat_l[41] = 0.223606797749979*alpha[22]*fUpwind_l[47]+0.223606797749979*alpha[21]*fUpwind_l[47]+0.2500000000000001*alpha[3]*fUpwind_l[47]+0.25*alpha[6]*fUpwind_l[43]+0.25*alpha[7]*fUpwind_l[42]+0.223606797749979*alpha[12]*fUpwind_l[41]+0.223606797749979*alpha[11]*fUpwind_l[41]+0.25*alpha[0]*fUpwind_l[41]+0.2*alpha[9]*fUpwind_l[36]+0.2*alpha[8]*fUpwind_l[35]+0.2500000000000001*alpha[15]*fUpwind_l[30]+0.2500000000000001*alpha[1]*fUpwind_l[29]+0.2500000000000001*alpha[2]*fUpwind_l[28]+0.2*alpha[16]*fUpwind_l[26]+0.2*alpha[16]*fUpwind_l[25]+0.2*fUpwind_l[16]*alpha[25]+0.223606797749979*alpha[4]*fUpwind_l[16]+0.223606797749979*fUpwind_l[4]*alpha[16]+0.25*alpha[5]*fUpwind_l[14]+0.223606797749979*alpha[8]*fUpwind_l[9]+0.223606797749979*fUpwind_l[8]*alpha[9]; 
  Ghat_l[42] = 0.2500000000000001*alpha[2]*fUpwind_l[47]+0.2*alpha[16]*fUpwind_l[44]+0.25*alpha[5]*fUpwind_l[43]+0.223606797749979*alpha[11]*fUpwind_l[42]+0.25*alpha[0]*fUpwind_l[42]+0.25*alpha[7]*fUpwind_l[41]+0.2*alpha[8]*fUpwind_l[37]+0.223606797749979*alpha[9]*fUpwind_l[31]+0.2500000000000001*alpha[1]*fUpwind_l[30]+0.2500000000000001*alpha[15]*fUpwind_l[29]+0.223606797749979*alpha[21]*fUpwind_l[28]+0.2500000000000001*alpha[3]*fUpwind_l[28]+0.2*fUpwind_l[17]*alpha[25]+0.223606797749979*alpha[16]*fUpwind_l[18]+0.223606797749979*alpha[4]*fUpwind_l[17]+0.25*alpha[6]*fUpwind_l[14]+0.223606797749979*alpha[8]*fUpwind_l[10]; 
  Ghat_l[43] = 0.2500000000000001*alpha[1]*fUpwind_l[47]+0.2*alpha[16]*fUpwind_l[45]+0.223606797749979*alpha[25]*fUpwind_l[44]+0.223606797749979*alpha[12]*fUpwind_l[43]+0.25*alpha[0]*fUpwind_l[43]+0.25*alpha[5]*fUpwind_l[42]+0.25*alpha[6]*fUpwind_l[41]+0.2*alpha[9]*fUpwind_l[38]+0.223606797749979*alpha[8]*fUpwind_l[31]+0.2500000000000001*alpha[2]*fUpwind_l[30]+0.223606797749979*alpha[22]*fUpwind_l[29]+0.2500000000000001*alpha[3]*fUpwind_l[29]+0.2500000000000001*alpha[15]*fUpwind_l[28]+0.223606797749979*alpha[4]*fUpwind_l[18]+0.223606797749979*alpha[16]*fUpwind_l[17]+0.25*alpha[7]*fUpwind_l[14]+0.223606797749979*alpha[9]*fUpwind_l[10]; 
  Ghat_l[44] = 0.2*alpha[8]*fUpwind_l[47]+0.2*alpha[6]*fUpwind_l[46]+0.2*alpha[5]*fUpwind_l[45]+0.223606797749979*alpha[12]*fUpwind_l[44]+0.159719141249985*alpha[11]*fUpwind_l[44]+0.25*alpha[0]*fUpwind_l[44]+0.223606797749979*alpha[25]*fUpwind_l[43]+0.2*alpha[16]*fUpwind_l[42]+0.223606797749979*alpha[21]*fUpwind_l[40]+0.2*alpha[15]*fUpwind_l[39]+0.2500000000000001*alpha[2]*fUpwind_l[37]+0.2*alpha[15]*fUpwind_l[36]+0.223606797749979*alpha[22]*fUpwind_l[35]+0.159719141249985*alpha[21]*fUpwind_l[35]+0.2500000000000001*alpha[3]*fUpwind_l[35]+0.2*alpha[16]*fUpwind_l[33]+0.159719141249985*alpha[25]*fUpwind_l[32]+0.2500000000000001*alpha[4]*fUpwind_l[32]+0.223606797749979*alpha[1]*fUpwind_l[31]+0.25*alpha[7]*fUpwind_l[25]+0.25*fUpwind_l[7]*alpha[25]+0.25*alpha[9]*fUpwind_l[21]+0.25*fUpwind_l[9]*alpha[21]+0.2500000000000001*alpha[11]*fUpwind_l[18]+0.223606797749979*alpha[5]*fUpwind_l[17]+0.223606797749979*alpha[6]*fUpwind_l[16]+0.223606797749979*fUpwind_l[6]*alpha[16]+0.223606797749979*alpha[8]*fUpwind_l[15]+0.223606797749979*fUpwind_l[8]*alpha[15]; 
  Ghat_l[45] = 0.2*alpha[9]*fUpwind_l[47]+0.2*alpha[7]*fUpwind_l[46]+0.159719141249985*alpha[12]*fUpwind_l[45]+0.223606797749979*alpha[11]*fUpwind_l[45]+0.25*alpha[0]*fUpwind_l[45]+0.2*alpha[5]*fUpwind_l[44]+0.2*alpha[16]*fUpwind_l[43]+0.2*alpha[15]*fUpwind_l[40]+0.223606797749979*alpha[22]*fUpwind_l[39]+0.2500000000000001*alpha[1]*fUpwind_l[38]+0.159719141249985*alpha[22]*fUpwind_l[36]+0.223606797749979*alpha[21]*fUpwind_l[36]+0.2500000000000001*alpha[3]*fUpwind_l[36]+0.2*alpha[15]*fUpwind_l[35]+0.223606797749979*alpha[25]*fUpwind_l[33]+0.2500000000000001*alpha[4]*fUpwind_l[33]+0.2*alpha[16]*fUpwind_l[32]+0.223606797749979*alpha[2]*fUpwind_l[31]+0.25*alpha[6]*fUpwind_l[26]+0.25*alpha[8]*fUpwind_l[22]+0.25*fUpwind_l[8]*alpha[22]+0.223606797749979*alpha[5]*fUpwind_l[18]+0.2500000000000001*alpha[12]*fUpwind_l[17]+0.223606797749979*alpha[7]*fUpwind_l[16]+0.223606797749979*fUpwind_l[7]*alpha[16]+0.223606797749979*alpha[9]*fUpwind_l[15]+0.223606797749979*fUpwind_l[9]*alpha[15]; 
  Ghat_l[46] = 0.223606797749979*alpha[12]*fUpwind_l[46]+0.223606797749979*alpha[11]*fUpwind_l[46]+0.25*alpha[0]*fUpwind_l[46]+0.2*alpha[7]*fUpwind_l[45]+0.2*alpha[6]*fUpwind_l[44]+0.2500000000000001*alpha[1]*fUpwind_l[40]+0.2500000000000001*alpha[2]*fUpwind_l[39]+0.2*alpha[15]*fUpwind_l[38]+0.2*alpha[15]*fUpwind_l[37]+0.223606797749979*alpha[25]*fUpwind_l[34]+0.2500000000000001*alpha[4]*fUpwind_l[34]+0.2*alpha[22]*fUpwind_l[31]+0.2*alpha[21]*fUpwind_l[31]+0.223606797749979*alpha[3]*fUpwind_l[31]+0.25*alpha[5]*fUpwind_l[27]+0.25*alpha[8]*fUpwind_l[24]+0.25*alpha[9]*fUpwind_l[23]+0.223606797749979*alpha[6]*fUpwind_l[18]+0.223606797749979*alpha[7]*fUpwind_l[17]+0.2500000000000001*fUpwind_l[13]*alpha[16]+0.223606797749979*fUpwind_l[10]*alpha[15]; 
  Ghat_l[47] = 0.223606797749979*alpha[12]*fUpwind_l[47]+0.223606797749979*alpha[11]*fUpwind_l[47]+0.25*alpha[0]*fUpwind_l[47]+0.2*alpha[9]*fUpwind_l[45]+0.2*alpha[8]*fUpwind_l[44]+0.2500000000000001*alpha[1]*fUpwind_l[43]+0.2500000000000001*alpha[2]*fUpwind_l[42]+0.223606797749979*alpha[22]*fUpwind_l[41]+0.223606797749979*alpha[21]*fUpwind_l[41]+0.2500000000000001*alpha[3]*fUpwind_l[41]+0.2*alpha[16]*fUpwind_l[38]+0.2*alpha[16]*fUpwind_l[37]+0.2*alpha[25]*fUpwind_l[31]+0.223606797749979*alpha[4]*fUpwind_l[31]+0.25*alpha[5]*fUpwind_l[30]+0.25*alpha[6]*fUpwind_l[29]+0.25*alpha[7]*fUpwind_l[28]+0.223606797749979*alpha[8]*fUpwind_l[18]+0.223606797749979*alpha[9]*fUpwind_l[17]+0.223606797749979*fUpwind_l[10]*alpha[16]+0.2500000000000001*fUpwind_l[14]*alpha[15]; 

  Ghat_r[0] = 0.25*alpha[25]*fUpwind_r[25]+0.25*alpha[22]*fUpwind_r[22]+0.25*alpha[21]*fUpwind_r[21]+0.25*alpha[16]*fUpwind_r[16]+0.25*alpha[15]*fUpwind_r[15]+0.25*alpha[12]*fUpwind_r[12]+0.25*alpha[11]*fUpwind_r[11]+0.25*alpha[9]*fUpwind_r[9]+0.25*alpha[8]*fUpwind_r[8]+0.25*alpha[7]*fUpwind_r[7]+0.25*alpha[6]*fUpwind_r[6]+0.25*alpha[5]*fUpwind_r[5]+0.25*alpha[4]*fUpwind_r[4]+0.25*alpha[3]*fUpwind_r[3]+0.25*alpha[2]*fUpwind_r[2]+0.25*alpha[1]*fUpwind_r[1]+0.25*alpha[0]*fUpwind_r[0]; 
  Ghat_r[1] = 0.223606797749979*alpha[16]*fUpwind_r[35]+0.2500000000000001*alpha[22]*fUpwind_r[33]+0.223606797749979*alpha[15]*fUpwind_r[32]+0.223606797749979*alpha[8]*fUpwind_r[25]+0.223606797749979*fUpwind_r[8]*alpha[25]+0.223606797749979*alpha[6]*fUpwind_r[21]+0.223606797749979*fUpwind_r[6]*alpha[21]+0.2500000000000001*alpha[12]*fUpwind_r[20]+0.223606797749979*alpha[5]*fUpwind_r[19]+0.25*alpha[9]*fUpwind_r[16]+0.25*fUpwind_r[9]*alpha[16]+0.25*alpha[7]*fUpwind_r[15]+0.25*fUpwind_r[7]*alpha[15]+0.223606797749979*alpha[1]*fUpwind_r[11]+0.223606797749979*fUpwind_r[1]*alpha[11]+0.25*alpha[4]*fUpwind_r[8]+0.25*fUpwind_r[4]*alpha[8]+0.25*alpha[3]*fUpwind_r[6]+0.25*fUpwind_r[3]*alpha[6]+0.25*alpha[2]*fUpwind_r[5]+0.25*fUpwind_r[2]*alpha[5]+0.25*alpha[0]*fUpwind_r[1]+0.25*fUpwind_r[0]*alpha[1]; 
  Ghat_r[2] = 0.223606797749979*alpha[16]*fUpwind_r[36]+0.2500000000000001*alpha[25]*fUpwind_r[35]+0.223606797749979*alpha[15]*fUpwind_r[33]+0.2500000000000001*alpha[21]*fUpwind_r[32]+0.223606797749979*alpha[9]*fUpwind_r[26]+0.223606797749979*alpha[7]*fUpwind_r[22]+0.223606797749979*fUpwind_r[7]*alpha[22]+0.223606797749979*alpha[5]*fUpwind_r[20]+0.2500000000000001*alpha[11]*fUpwind_r[19]+0.25*alpha[8]*fUpwind_r[16]+0.25*fUpwind_r[8]*alpha[16]+0.25*alpha[6]*fUpwind_r[15]+0.25*fUpwind_r[6]*alpha[15]+0.223606797749979*alpha[2]*fUpwind_r[12]+0.223606797749979*fUpwind_r[2]*alpha[12]+0.25*alpha[4]*fUpwind_r[9]+0.25*fUpwind_r[4]*alpha[9]+0.25*alpha[3]*fUpwind_r[7]+0.25*fUpwind_r[3]*alpha[7]+0.25*alpha[1]*fUpwind_r[5]+0.25*fUpwind_r[1]*alpha[5]+0.25*alpha[0]*fUpwind_r[2]+0.25*fUpwind_r[0]*alpha[2]; 
  Ghat_r[3] = 0.2500000000000001*alpha[25]*fUpwind_r[37]+0.223606797749979*alpha[15]*fUpwind_r[34]+0.25*alpha[16]*fUpwind_r[31]+0.223606797749979*alpha[7]*fUpwind_r[24]+0.223606797749979*alpha[6]*fUpwind_r[23]+0.2500000000000001*alpha[12]*fUpwind_r[22]+0.2500000000000001*fUpwind_r[12]*alpha[22]+0.2500000000000001*alpha[11]*fUpwind_r[21]+0.2500000000000001*fUpwind_r[11]*alpha[21]+0.25*alpha[9]*fUpwind_r[18]+0.25*alpha[8]*fUpwind_r[17]+0.25*alpha[5]*fUpwind_r[15]+0.25*fUpwind_r[5]*alpha[15]+0.223606797749979*alpha[3]*fUpwind_r[13]+0.25*alpha[4]*fUpwind_r[10]+0.25*alpha[2]*fUpwind_r[7]+0.25*fUpwind_r[2]*alpha[7]+0.25*alpha[1]*fUpwind_r[6]+0.25*fUpwind_r[1]*alpha[6]+0.25*alpha[0]*fUpwind_r[3]+0.25*fUpwind_r[0]*alpha[3]; 
  Ghat_r[4] = 0.223606797749979*alpha[16]*fUpwind_r[41]+0.2500000000000001*alpha[22]*fUpwind_r[38]+0.2500000000000001*alpha[21]*fUpwind_r[37]+0.25*alpha[15]*fUpwind_r[31]+0.223606797749979*alpha[9]*fUpwind_r[29]+0.223606797749979*alpha[8]*fUpwind_r[28]+0.2500000000000001*alpha[12]*fUpwind_r[26]+0.2500000000000001*alpha[11]*fUpwind_r[25]+0.2500000000000001*fUpwind_r[11]*alpha[25]+0.25*alpha[7]*fUpwind_r[18]+0.25*alpha[6]*fUpwind_r[17]+0.25*alpha[5]*fUpwind_r[16]+0.25*fUpwind_r[5]*alpha[16]+0.223606797749979*alpha[4]*fUpwind_r[14]+0.25*alpha[3]*fUpwind_r[10]+0.25*alpha[2]*fUpwind_r[9]+0.25*fUpwind_r[2]*alpha[9]+0.25*alpha[1]*fUpwind_r[8]+0.25*fUpwind_r[1]*alpha[8]+0.25*alpha[0]*fUpwind_r[4]+0.25*fUpwind_r[0]*alpha[4]; 
  Ghat_r[5] = 0.223606797749979*alpha[9]*fUpwind_r[36]+0.223606797749979*alpha[8]*fUpwind_r[35]+0.223606797749979*alpha[7]*fUpwind_r[33]+0.223606797749979*alpha[6]*fUpwind_r[32]+0.223606797749979*alpha[16]*fUpwind_r[26]+0.223606797749979*alpha[16]*fUpwind_r[25]+0.223606797749979*fUpwind_r[16]*alpha[25]+0.223606797749979*alpha[15]*fUpwind_r[22]+0.223606797749979*fUpwind_r[15]*alpha[22]+0.223606797749979*alpha[15]*fUpwind_r[21]+0.223606797749979*fUpwind_r[15]*alpha[21]+0.223606797749979*alpha[2]*fUpwind_r[20]+0.223606797749979*alpha[1]*fUpwind_r[19]+0.25*alpha[4]*fUpwind_r[16]+0.25*fUpwind_r[4]*alpha[16]+0.25*alpha[3]*fUpwind_r[15]+0.25*fUpwind_r[3]*alpha[15]+0.223606797749979*alpha[5]*fUpwind_r[12]+0.223606797749979*fUpwind_r[5]*alpha[12]+0.223606797749979*alpha[5]*fUpwind_r[11]+0.223606797749979*fUpwind_r[5]*alpha[11]+0.25*alpha[8]*fUpwind_r[9]+0.25*fUpwind_r[8]*alpha[9]+0.25*alpha[6]*fUpwind_r[7]+0.25*fUpwind_r[6]*alpha[7]+0.25*alpha[0]*fUpwind_r[5]+0.25*fUpwind_r[0]*alpha[5]+0.25*alpha[1]*fUpwind_r[2]+0.25*fUpwind_r[1]*alpha[2]; 
  Ghat_r[6] = 0.223606797749979*alpha[16]*fUpwind_r[44]+0.223606797749979*alpha[8]*fUpwind_r[37]+0.223606797749979*alpha[7]*fUpwind_r[34]+0.25*alpha[12]*fUpwind_r[33]+0.223606797749979*alpha[5]*fUpwind_r[32]+0.25*alpha[9]*fUpwind_r[31]+0.223606797749979*fUpwind_r[17]*alpha[25]+0.223606797749979*alpha[15]*fUpwind_r[24]+0.2*alpha[21]*fUpwind_r[23]+0.223606797749979*alpha[3]*fUpwind_r[23]+0.25*fUpwind_r[20]*alpha[22]+0.223606797749979*alpha[1]*fUpwind_r[21]+0.223606797749979*fUpwind_r[1]*alpha[21]+0.223606797749979*alpha[15]*fUpwind_r[19]+0.25*alpha[16]*fUpwind_r[18]+0.25*alpha[4]*fUpwind_r[17]+0.25*alpha[2]*fUpwind_r[15]+0.25*fUpwind_r[2]*alpha[15]+0.223606797749979*alpha[6]*fUpwind_r[13]+0.223606797749979*alpha[6]*fUpwind_r[11]+0.223606797749979*fUpwind_r[6]*alpha[11]+0.25*alpha[8]*fUpwind_r[10]+0.25*alpha[5]*fUpwind_r[7]+0.25*fUpwind_r[5]*alpha[7]+0.25*alpha[0]*fUpwind_r[6]+0.25*fUpwind_r[0]*alpha[6]+0.25*alpha[1]*fUpwind_r[3]+0.25*fUpwind_r[1]*alpha[3]; 
  Ghat_r[7] = 0.223606797749979*alpha[16]*fUpwind_r[45]+0.25*alpha[25]*fUpwind_r[44]+0.223606797749979*alpha[9]*fUpwind_r[38]+0.223606797749979*alpha[6]*fUpwind_r[34]+0.223606797749979*alpha[5]*fUpwind_r[33]+0.25*alpha[11]*fUpwind_r[32]+0.25*alpha[8]*fUpwind_r[31]+0.2*alpha[22]*fUpwind_r[24]+0.223606797749979*alpha[3]*fUpwind_r[24]+0.223606797749979*alpha[15]*fUpwind_r[23]+0.223606797749979*alpha[2]*fUpwind_r[22]+0.223606797749979*fUpwind_r[2]*alpha[22]+0.25*fUpwind_r[19]*alpha[21]+0.223606797749979*alpha[15]*fUpwind_r[20]+0.25*alpha[4]*fUpwind_r[18]+0.25*alpha[16]*fUpwind_r[17]+0.25*alpha[1]*fUpwind_r[15]+0.25*fUpwind_r[1]*alpha[15]+0.223606797749979*alpha[7]*fUpwind_r[13]+0.223606797749979*alpha[7]*fUpwind_r[12]+0.223606797749979*fUpwind_r[7]*alpha[12]+0.25*alpha[9]*fUpwind_r[10]+0.25*alpha[0]*fUpwind_r[7]+0.25*fUpwind_r[0]*alpha[7]+0.25*alpha[5]*fUpwind_r[6]+0.25*fUpwind_r[5]*alpha[6]+0.25*alpha[2]*fUpwind_r[3]+0.25*fUpwind_r[2]*alpha[3]; 
  Ghat_r[8] = 0.25*alpha[22]*fUpwind_r[45]+0.223606797749979*alpha[15]*fUpwind_r[44]+0.223606797749979*alpha[9]*fUpwind_r[41]+0.223606797749979*alpha[6]*fUpwind_r[37]+0.25*alpha[12]*fUpwind_r[36]+0.223606797749979*alpha[5]*fUpwind_r[35]+0.25*alpha[7]*fUpwind_r[31]+0.223606797749979*alpha[16]*fUpwind_r[29]+0.2*alpha[25]*fUpwind_r[28]+0.223606797749979*alpha[4]*fUpwind_r[28]+0.223606797749979*alpha[1]*fUpwind_r[25]+0.223606797749979*fUpwind_r[1]*alpha[25]+0.223606797749979*fUpwind_r[17]*alpha[21]+0.223606797749979*alpha[16]*fUpwind_r[19]+0.25*alpha[15]*fUpwind_r[18]+0.25*alpha[3]*fUpwind_r[17]+0.25*alpha[2]*fUpwind_r[16]+0.25*fUpwind_r[2]*alpha[16]+0.223606797749979*alpha[8]*fUpwind_r[14]+0.223606797749979*alpha[8]*fUpwind_r[11]+0.223606797749979*fUpwind_r[8]*alpha[11]+0.25*alpha[6]*fUpwind_r[10]+0.25*alpha[5]*fUpwind_r[9]+0.25*fUpwind_r[5]*alpha[9]+0.25*alpha[0]*fUpwind_r[8]+0.25*fUpwind_r[0]*alpha[8]+0.25*alpha[1]*fUpwind_r[4]+0.25*fUpwind_r[1]*alpha[4]; 
  Ghat_r[9] = 0.223606797749979*alpha[15]*fUpwind_r[45]+0.25*alpha[21]*fUpwind_r[44]+0.223606797749979*alpha[8]*fUpwind_r[41]+0.223606797749979*alpha[7]*fUpwind_r[38]+0.223606797749979*alpha[5]*fUpwind_r[36]+0.25*alpha[11]*fUpwind_r[35]+0.25*alpha[6]*fUpwind_r[31]+0.223606797749979*alpha[4]*fUpwind_r[29]+0.223606797749979*alpha[16]*fUpwind_r[28]+0.223606797749979*alpha[2]*fUpwind_r[26]+0.25*fUpwind_r[19]*alpha[25]+0.223606797749979*fUpwind_r[18]*alpha[22]+0.223606797749979*alpha[16]*fUpwind_r[20]+0.25*alpha[3]*fUpwind_r[18]+0.25*alpha[15]*fUpwind_r[17]+0.25*alpha[1]*fUpwind_r[16]+0.25*fUpwind_r[1]*alpha[16]+0.223606797749979*alpha[9]*fUpwind_r[14]+0.223606797749979*alpha[9]*fUpwind_r[12]+0.223606797749979*fUpwind_r[9]*alpha[12]+0.25*alpha[7]*fUpwind_r[10]+0.25*alpha[0]*fUpwind_r[9]+0.25*fUpwind_r[0]*alpha[9]+0.25*alpha[5]*fUpwind_r[8]+0.25*fUpwind_r[5]*alpha[8]+0.25*alpha[2]*fUpwind_r[4]+0.25*fUpwind_r[2]*alpha[4]; 
  Ghat_r[10] = 0.223606797749979*alpha[16]*fUpwind_r[47]+0.223606797749979*alpha[15]*fUpwind_r[46]+0.223606797749979*alpha[9]*fUpwind_r[43]+0.223606797749979*alpha[8]*fUpwind_r[42]+0.223606797749979*alpha[7]*fUpwind_r[40]+0.223606797749979*alpha[6]*fUpwind_r[39]+0.25*alpha[12]*fUpwind_r[38]+0.25*alpha[11]*fUpwind_r[37]+0.25*alpha[5]*fUpwind_r[31]+0.223606797749979*alpha[4]*fUpwind_r[30]+0.223606797749979*alpha[3]*fUpwind_r[27]+0.25*alpha[22]*fUpwind_r[26]+0.25*alpha[21]*fUpwind_r[25]+0.25*fUpwind_r[21]*alpha[25]+0.25*alpha[2]*fUpwind_r[18]+0.25*alpha[1]*fUpwind_r[17]+0.25*alpha[15]*fUpwind_r[16]+0.25*fUpwind_r[15]*alpha[16]+0.25*alpha[0]*fUpwind_r[10]+0.25*alpha[7]*fUpwind_r[9]+0.25*fUpwind_r[7]*alpha[9]+0.25*alpha[6]*fUpwind_r[8]+0.25*fUpwind_r[6]*alpha[8]+0.25*alpha[3]*fUpwind_r[4]+0.25*fUpwind_r[3]*alpha[4]; 
  Ghat_r[11] = 0.25*alpha[9]*fUpwind_r[35]+0.25*alpha[7]*fUpwind_r[32]+0.159719141249985*alpha[25]*fUpwind_r[25]+0.2500000000000001*alpha[4]*fUpwind_r[25]+0.2500000000000001*fUpwind_r[4]*alpha[25]+0.159719141249985*alpha[21]*fUpwind_r[21]+0.2500000000000001*alpha[3]*fUpwind_r[21]+0.2500000000000001*fUpwind_r[3]*alpha[21]+0.2500000000000001*alpha[2]*fUpwind_r[19]+0.223606797749979*alpha[16]*fUpwind_r[16]+0.223606797749979*alpha[15]*fUpwind_r[15]+0.159719141249985*alpha[11]*fUpwind_r[11]+0.25*alpha[0]*fUpwind_r[11]+0.25*fUpwind_r[0]*alpha[11]+0.223606797749979*alpha[8]*fUpwind_r[8]+0.223606797749979*alpha[6]*fUpwind_r[6]+0.223606797749979*alpha[5]*fUpwind_r[5]+0.223606797749979*alpha[1]*fUpwind_r[1]; 
  Ghat_r[12] = 0.25*alpha[8]*fUpwind_r[36]+0.25*alpha[6]*fUpwind_r[33]+0.2500000000000001*alpha[4]*fUpwind_r[26]+0.159719141249985*alpha[22]*fUpwind_r[22]+0.2500000000000001*alpha[3]*fUpwind_r[22]+0.2500000000000001*fUpwind_r[3]*alpha[22]+0.2500000000000001*alpha[1]*fUpwind_r[20]+0.223606797749979*alpha[16]*fUpwind_r[16]+0.223606797749979*alpha[15]*fUpwind_r[15]+0.159719141249985*alpha[12]*fUpwind_r[12]+0.25*alpha[0]*fUpwind_r[12]+0.25*fUpwind_r[0]*alpha[12]+0.223606797749979*alpha[9]*fUpwind_r[9]+0.223606797749979*alpha[7]*fUpwind_r[7]+0.223606797749979*alpha[5]*fUpwind_r[5]+0.223606797749979*alpha[2]*fUpwind_r[2]; 
  Ghat_r[13] = 0.2500000000000001*alpha[16]*fUpwind_r[46]+0.25*alpha[9]*fUpwind_r[40]+0.25*alpha[8]*fUpwind_r[39]+0.25*alpha[5]*fUpwind_r[34]+0.2500000000000001*alpha[4]*fUpwind_r[27]+0.2500000000000001*alpha[2]*fUpwind_r[24]+0.2500000000000001*alpha[1]*fUpwind_r[23]+0.223606797749979*alpha[22]*fUpwind_r[22]+0.223606797749979*alpha[21]*fUpwind_r[21]+0.223606797749979*alpha[15]*fUpwind_r[15]+0.25*alpha[0]*fUpwind_r[13]+0.223606797749979*alpha[7]*fUpwind_r[7]+0.223606797749979*alpha[6]*fUpwind_r[6]+0.223606797749979*alpha[3]*fUpwind_r[3]; 
  Ghat_r[14] = 0.2500000000000001*alpha[15]*fUpwind_r[47]+0.25*alpha[7]*fUpwind_r[43]+0.25*alpha[6]*fUpwind_r[42]+0.25*alpha[5]*fUpwind_r[41]+0.2500000000000001*alpha[3]*fUpwind_r[30]+0.2500000000000001*alpha[2]*fUpwind_r[29]+0.2500000000000001*alpha[1]*fUpwind_r[28]+0.223606797749979*alpha[25]*fUpwind_r[25]+0.223606797749979*alpha[16]*fUpwind_r[16]+0.25*alpha[0]*fUpwind_r[14]+0.223606797749979*alpha[9]*fUpwind_r[9]+0.223606797749979*alpha[8]*fUpwind_r[8]+0.223606797749979*alpha[4]*fUpwind_r[4]; 
  Ghat_r[15] = 0.223606797749979*alpha[9]*fUpwind_r[45]+0.223606797749979*alpha[8]*fUpwind_r[44]+0.223606797749979*alpha[16]*fUpwind_r[38]+0.223606797749979*alpha[16]*fUpwind_r[37]+0.2*alpha[22]*fUpwind_r[34]+0.2*alpha[21]*fUpwind_r[34]+0.223606797749979*alpha[3]*fUpwind_r[34]+0.223606797749979*alpha[2]*fUpwind_r[33]+0.223606797749979*alpha[1]*fUpwind_r[32]+0.223606797749979*alpha[25]*fUpwind_r[31]+0.25*alpha[4]*fUpwind_r[31]+0.223606797749979*alpha[6]*fUpwind_r[24]+0.223606797749979*alpha[7]*fUpwind_r[23]+0.223606797749979*alpha[5]*fUpwind_r[22]+0.223606797749979*fUpwind_r[5]*alpha[22]+0.223606797749979*alpha[5]*fUpwind_r[21]+0.223606797749979*fUpwind_r[5]*alpha[21]+0.223606797749979*alpha[7]*fUpwind_r[20]+0.223606797749979*alpha[6]*fUpwind_r[19]+0.25*alpha[8]*fUpwind_r[18]+0.25*alpha[9]*fUpwind_r[17]+0.25*fUpwind_r[10]*alpha[16]+0.223606797749979*alpha[12]*fUpwind_r[15]+0.223606797749979*alpha[11]*fUpwind_r[15]+0.25*alpha[0]*fUpwind_r[15]+0.223606797749979*fUpwind_r[13]*alpha[15]+0.223606797749979*fUpwind_r[12]*alpha[15]+0.223606797749979*fUpwind_r[11]*alpha[15]+0.25*fUpwind_r[0]*alpha[15]+0.25*alpha[1]*fUpwind_r[7]+0.25*fUpwind_r[1]*alpha[7]+0.25*alpha[2]*fUpwind_r[6]+0.25*fUpwind_r[2]*alpha[6]+0.25*alpha[3]*fUpwind_r[5]+0.25*fUpwind_r[3]*alpha[5]; 
  Ghat_r[16] = 0.223606797749979*alpha[7]*fUpwind_r[45]+0.223606797749979*alpha[6]*fUpwind_r[44]+0.2*alpha[25]*fUpwind_r[41]+0.223606797749979*alpha[4]*fUpwind_r[41]+0.223606797749979*alpha[15]*fUpwind_r[38]+0.223606797749979*alpha[15]*fUpwind_r[37]+0.223606797749979*alpha[2]*fUpwind_r[36]+0.223606797749979*alpha[1]*fUpwind_r[35]+0.223606797749979*alpha[22]*fUpwind_r[31]+0.223606797749979*alpha[21]*fUpwind_r[31]+0.25*alpha[3]*fUpwind_r[31]+0.223606797749979*alpha[8]*fUpwind_r[29]+0.223606797749979*alpha[9]*fUpwind_r[28]+0.223606797749979*alpha[5]*fUpwind_r[26]+0.223606797749979*alpha[5]*fUpwind_r[25]+0.223606797749979*fUpwind_r[5]*alpha[25]+0.223606797749979*alpha[9]*fUpwind_r[20]+0.223606797749979*alpha[8]*fUpwind_r[19]+0.25*alpha[6]*fUpwind_r[18]+0.25*alpha[7]*fUpwind_r[17]+0.223606797749979*alpha[12]*fUpwind_r[16]+0.223606797749979*alpha[11]*fUpwind_r[16]+0.25*alpha[0]*fUpwind_r[16]+0.223606797749979*fUpwind_r[14]*alpha[16]+0.223606797749979*fUpwind_r[12]*alpha[16]+0.223606797749979*fUpwind_r[11]*alpha[16]+0.25*fUpwind_r[0]*alpha[16]+0.25*fUpwind_r[10]*alpha[15]+0.25*alpha[1]*fUpwind_r[9]+0.25*fUpwind_r[1]*alpha[9]+0.25*alpha[2]*fUpwind_r[8]+0.25*fUpwind_r[2]*alpha[8]+0.25*alpha[4]*fUpwind_r[5]+0.25*fUpwind_r[4]*alpha[5]; 
  Ghat_r[17] = 0.223606797749979*alpha[9]*fUpwind_r[47]+0.223606797749979*alpha[7]*fUpwind_r[46]+0.2500000000000001*alpha[12]*fUpwind_r[45]+0.223606797749979*alpha[5]*fUpwind_r[44]+0.223606797749979*alpha[16]*fUpwind_r[43]+0.2*alpha[25]*fUpwind_r[42]+0.223606797749979*alpha[4]*fUpwind_r[42]+0.223606797749979*alpha[15]*fUpwind_r[40]+0.2*alpha[21]*fUpwind_r[39]+0.223606797749979*alpha[3]*fUpwind_r[39]+0.223606797749979*alpha[1]*fUpwind_r[37]+0.2500000000000001*alpha[22]*fUpwind_r[36]+0.223606797749979*alpha[15]*fUpwind_r[35]+0.223606797749979*alpha[16]*fUpwind_r[32]+0.25*alpha[2]*fUpwind_r[31]+0.223606797749979*alpha[8]*fUpwind_r[30]+0.223606797749979*alpha[6]*fUpwind_r[27]+0.223606797749979*alpha[6]*fUpwind_r[25]+0.223606797749979*fUpwind_r[6]*alpha[25]+0.223606797749979*alpha[8]*fUpwind_r[21]+0.223606797749979*fUpwind_r[8]*alpha[21]+0.25*alpha[5]*fUpwind_r[18]+0.223606797749979*alpha[11]*fUpwind_r[17]+0.25*alpha[0]*fUpwind_r[17]+0.25*alpha[7]*fUpwind_r[16]+0.25*fUpwind_r[7]*alpha[16]+0.25*alpha[9]*fUpwind_r[15]+0.25*fUpwind_r[9]*alpha[15]+0.25*alpha[1]*fUpwind_r[10]+0.25*alpha[3]*fUpwind_r[8]+0.25*fUpwind_r[3]*alpha[8]+0.25*alpha[4]*fUpwind_r[6]+0.25*fUpwind_r[4]*alpha[6]; 
  Ghat_r[18] = 0.223606797749979*alpha[8]*fUpwind_r[47]+0.223606797749979*alpha[6]*fUpwind_r[46]+0.223606797749979*alpha[5]*fUpwind_r[45]+0.2500000000000001*alpha[11]*fUpwind_r[44]+0.223606797749979*alpha[4]*fUpwind_r[43]+0.223606797749979*alpha[16]*fUpwind_r[42]+0.2*alpha[22]*fUpwind_r[40]+0.223606797749979*alpha[3]*fUpwind_r[40]+0.223606797749979*alpha[15]*fUpwind_r[39]+0.223606797749979*alpha[2]*fUpwind_r[38]+0.223606797749979*alpha[15]*fUpwind_r[36]+0.2500000000000001*alpha[21]*fUpwind_r[35]+0.223606797749979*alpha[16]*fUpwind_r[33]+0.2500000000000001*alpha[25]*fUpwind_r[32]+0.25*alpha[1]*fUpwind_r[31]+0.223606797749979*alpha[9]*fUpwind_r[30]+0.223606797749979*alpha[7]*fUpwind_r[27]+0.223606797749979*alpha[7]*fUpwind_r[26]+0.223606797749979*alpha[9]*fUpwind_r[22]+0.223606797749979*fUpwind_r[9]*alpha[22]+0.223606797749979*alpha[12]*fUpwind_r[18]+0.25*alpha[0]*fUpwind_r[18]+0.25*alpha[5]*fUpwind_r[17]+0.25*alpha[6]*fUpwind_r[16]+0.25*fUpwind_r[6]*alpha[16]+0.25*alpha[8]*fUpwind_r[15]+0.25*fUpwind_r[8]*alpha[15]+0.25*alpha[2]*fUpwind_r[10]+0.25*alpha[3]*fUpwind_r[9]+0.25*fUpwind_r[3]*alpha[9]+0.25*alpha[4]*fUpwind_r[7]+0.25*fUpwind_r[4]*alpha[7]; 
  Ghat_r[19] = 0.2*alpha[16]*fUpwind_r[36]+0.159719141249985*alpha[25]*fUpwind_r[35]+0.2500000000000001*alpha[4]*fUpwind_r[35]+0.2*alpha[15]*fUpwind_r[33]+0.223606797749979*alpha[22]*fUpwind_r[32]+0.159719141249985*alpha[21]*fUpwind_r[32]+0.2500000000000001*alpha[3]*fUpwind_r[32]+0.25*alpha[9]*fUpwind_r[25]+0.25*fUpwind_r[9]*alpha[25]+0.25*alpha[7]*fUpwind_r[21]+0.25*fUpwind_r[7]*alpha[21]+0.2*alpha[5]*fUpwind_r[20]+0.223606797749979*alpha[12]*fUpwind_r[19]+0.159719141249985*alpha[11]*fUpwind_r[19]+0.25*alpha[0]*fUpwind_r[19]+0.223606797749979*alpha[8]*fUpwind_r[16]+0.223606797749979*fUpwind_r[8]*alpha[16]+0.223606797749979*alpha[6]*fUpwind_r[15]+0.223606797749979*fUpwind_r[6]*alpha[15]+0.2500000000000001*alpha[2]*fUpwind_r[11]+0.2500000000000001*fUpwind_r[2]*alpha[11]+0.223606797749979*alpha[1]*fUpwind_r[5]+0.223606797749979*fUpwind_r[1]*alpha[5]; 
  Ghat_r[20] = 0.223606797749979*alpha[25]*fUpwind_r[36]+0.2500000000000001*alpha[4]*fUpwind_r[36]+0.2*alpha[16]*fUpwind_r[35]+0.159719141249985*alpha[22]*fUpwind_r[33]+0.223606797749979*alpha[21]*fUpwind_r[33]+0.2500000000000001*alpha[3]*fUpwind_r[33]+0.2*alpha[15]*fUpwind_r[32]+0.25*alpha[8]*fUpwind_r[26]+0.25*alpha[6]*fUpwind_r[22]+0.25*fUpwind_r[6]*alpha[22]+0.159719141249985*alpha[12]*fUpwind_r[20]+0.223606797749979*alpha[11]*fUpwind_r[20]+0.25*alpha[0]*fUpwind_r[20]+0.2*alpha[5]*fUpwind_r[19]+0.223606797749979*alpha[9]*fUpwind_r[16]+0.223606797749979*fUpwind_r[9]*alpha[16]+0.223606797749979*alpha[7]*fUpwind_r[15]+0.223606797749979*fUpwind_r[7]*alpha[15]+0.2500000000000001*alpha[1]*fUpwind_r[12]+0.2500000000000001*fUpwind_r[1]*alpha[12]+0.223606797749979*alpha[2]*fUpwind_r[5]+0.223606797749979*fUpwind_r[2]*alpha[5]; 
  Ghat_r[21] = 0.25*alpha[9]*fUpwind_r[44]+0.159719141249985*alpha[25]*fUpwind_r[37]+0.2500000000000001*alpha[4]*fUpwind_r[37]+0.2*alpha[15]*fUpwind_r[34]+0.2500000000000001*alpha[2]*fUpwind_r[32]+0.223606797749979*alpha[16]*fUpwind_r[31]+0.25*fUpwind_r[10]*alpha[25]+0.2*alpha[6]*fUpwind_r[23]+0.159719141249985*alpha[11]*fUpwind_r[21]+0.25*alpha[0]*fUpwind_r[21]+0.223606797749979*fUpwind_r[13]*alpha[21]+0.159719141249985*fUpwind_r[11]*alpha[21]+0.25*fUpwind_r[0]*alpha[21]+0.25*alpha[7]*fUpwind_r[19]+0.223606797749979*alpha[8]*fUpwind_r[17]+0.223606797749979*alpha[5]*fUpwind_r[15]+0.223606797749979*fUpwind_r[5]*alpha[15]+0.2500000000000001*alpha[3]*fUpwind_r[11]+0.2500000000000001*fUpwind_r[3]*alpha[11]+0.223606797749979*alpha[1]*fUpwind_r[6]+0.223606797749979*fUpwind_r[1]*alpha[6]; 
  Ghat_r[22] = 0.25*alpha[8]*fUpwind_r[45]+0.2500000000000001*alpha[4]*fUpwind_r[38]+0.2*alpha[15]*fUpwind_r[34]+0.2500000000000001*alpha[1]*fUpwind_r[33]+0.223606797749979*alpha[16]*fUpwind_r[31]+0.2*alpha[7]*fUpwind_r[24]+0.159719141249985*alpha[12]*fUpwind_r[22]+0.25*alpha[0]*fUpwind_r[22]+0.223606797749979*fUpwind_r[13]*alpha[22]+0.159719141249985*fUpwind_r[12]*alpha[22]+0.25*fUpwind_r[0]*alpha[22]+0.25*alpha[6]*fUpwind_r[20]+0.223606797749979*alpha[9]*fUpwind_r[18]+0.223606797749979*alpha[5]*fUpwind_r[15]+0.223606797749979*fUpwind_r[5]*alpha[15]+0.2500000000000001*alpha[3]*fUpwind_r[12]+0.2500000000000001*fUpwind_r[3]*alpha[12]+0.223606797749979*alpha[2]*fUpwind_r[7]+0.223606797749979*fUpwind_r[2]*alpha[7]; 
  Ghat_r[23] = 0.25*alpha[9]*fUpwind_r[46]+0.2500000000000001*alpha[16]*fUpwind_r[40]+0.223606797749979*alpha[25]*fUpwind_r[39]+0.2500000000000001*alpha[4]*fUpwind_r[39]+0.2500000000000001*alpha[2]*fUpwind_r[34]+0.223606797749979*alpha[22]*fUpwind_r[33]+0.2*alpha[15]*fUpwind_r[32]+0.25*alpha[8]*fUpwind_r[27]+0.25*alpha[5]*fUpwind_r[24]+0.223606797749979*alpha[11]*fUpwind_r[23]+0.25*alpha[0]*fUpwind_r[23]+0.2*alpha[6]*fUpwind_r[21]+0.2*fUpwind_r[6]*alpha[21]+0.223606797749979*alpha[7]*fUpwind_r[15]+0.223606797749979*fUpwind_r[7]*alpha[15]+0.2500000000000001*alpha[1]*fUpwind_r[13]+0.223606797749979*alpha[3]*fUpwind_r[6]+0.223606797749979*fUpwind_r[3]*alpha[6]; 
  Ghat_r[24] = 0.25*alpha[8]*fUpwind_r[46]+0.2500000000000001*alpha[4]*fUpwind_r[40]+0.2500000000000001*alpha[16]*fUpwind_r[39]+0.2500000000000001*alpha[1]*fUpwind_r[34]+0.2*alpha[15]*fUpwind_r[33]+0.223606797749979*alpha[21]*fUpwind_r[32]+0.25*alpha[9]*fUpwind_r[27]+0.223606797749979*alpha[12]*fUpwind_r[24]+0.25*alpha[0]*fUpwind_r[24]+0.25*alpha[5]*fUpwind_r[23]+0.2*alpha[7]*fUpwind_r[22]+0.2*fUpwind_r[7]*alpha[22]+0.223606797749979*alpha[6]*fUpwind_r[15]+0.223606797749979*fUpwind_r[6]*alpha[15]+0.2500000000000001*alpha[2]*fUpwind_r[13]+0.223606797749979*alpha[3]*fUpwind_r[7]+0.223606797749979*fUpwind_r[3]*alpha[7]; 
  Ghat_r[25] = 0.25*alpha[7]*fUpwind_r[44]+0.2*alpha[16]*fUpwind_r[41]+0.159719141249985*alpha[21]*fUpwind_r[37]+0.2500000000000001*alpha[3]*fUpwind_r[37]+0.2500000000000001*alpha[2]*fUpwind_r[35]+0.223606797749979*alpha[15]*fUpwind_r[31]+0.2*alpha[8]*fUpwind_r[28]+0.159719141249985*alpha[11]*fUpwind_r[25]+0.25*alpha[0]*fUpwind_r[25]+0.223606797749979*fUpwind_r[14]*alpha[25]+0.159719141249985*fUpwind_r[11]*alpha[25]+0.25*fUpwind_r[0]*alpha[25]+0.25*fUpwind_r[10]*alpha[21]+0.25*alpha[9]*fUpwind_r[19]+0.223606797749979*alpha[6]*fUpwind_r[17]+0.223606797749979*alpha[5]*fUpwind_r[16]+0.223606797749979*fUpwind_r[5]*alpha[16]+0.2500000000000001*alpha[4]*fUpwind_r[11]+0.2500000000000001*fUpwind_r[4]*alpha[11]+0.223606797749979*alpha[1]*fUpwind_r[8]+0.223606797749979*fUpwind_r[1]*alpha[8]; 
  Ghat_r[26] = 0.25*alpha[6]*fUpwind_r[45]+0.2*alpha[16]*fUpwind_r[41]+0.159719141249985*alpha[22]*fUpwind_r[38]+0.2500000000000001*alpha[3]*fUpwind_r[38]+0.2500000000000001*alpha[1]*fUpwind_r[36]+0.223606797749979*alpha[15]*fUpwind_r[31]+0.2*alpha[9]*fUpwind_r[29]+0.159719141249985*alpha[12]*fUpwind_r[26]+0.25*alpha[0]*fUpwind_r[26]+0.25*fUpwind_r[10]*alpha[22]+0.25*alpha[8]*fUpwind_r[20]+0.223606797749979*alpha[7]*fUpwind_r[18]+0.223606797749979*alpha[5]*fUpwind_r[16]+0.223606797749979*fUpwind_r[5]*alpha[16]+0.2500000000000001*alpha[4]*fUpwind_r[12]+0.2500000000000001*fUpwind_r[4]*alpha[12]+0.223606797749979*alpha[2]*fUpwind_r[9]+0.223606797749979*fUpwind_r[2]*alpha[9]; 
  Ghat_r[27] = 0.25*alpha[5]*fUpwind_r[46]+0.2500000000000001*alpha[2]*fUpwind_r[40]+0.2500000000000001*alpha[1]*fUpwind_r[39]+0.223606797749979*alpha[22]*fUpwind_r[38]+0.223606797749979*alpha[21]*fUpwind_r[37]+0.2500000000000001*alpha[16]*fUpwind_r[34]+0.223606797749979*alpha[15]*fUpwind_r[31]+0.25*alpha[0]*fUpwind_r[27]+0.25*alpha[9]*fUpwind_r[24]+0.25*alpha[8]*fUpwind_r[23]+0.223606797749979*alpha[7]*fUpwind_r[18]+0.223606797749979*alpha[6]*fUpwind_r[17]+0.2500000000000001*alpha[4]*fUpwind_r[13]+0.223606797749979*alpha[3]*fUpwind_r[10]; 
  Ghat_r[28] = 0.25*alpha[7]*fUpwind_r[47]+0.2500000000000001*alpha[15]*fUpwind_r[43]+0.223606797749979*alpha[21]*fUpwind_r[42]+0.2500000000000001*alpha[3]*fUpwind_r[42]+0.2500000000000001*alpha[2]*fUpwind_r[41]+0.2*alpha[16]*fUpwind_r[35]+0.25*alpha[6]*fUpwind_r[30]+0.25*alpha[5]*fUpwind_r[29]+0.223606797749979*alpha[11]*fUpwind_r[28]+0.25*alpha[0]*fUpwind_r[28]+0.2*alpha[8]*fUpwind_r[25]+0.2*fUpwind_r[8]*alpha[25]+0.223606797749979*alpha[9]*fUpwind_r[16]+0.223606797749979*fUpwind_r[9]*alpha[16]+0.2500000000000001*alpha[1]*fUpwind_r[14]+0.223606797749979*alpha[4]*fUpwind_r[8]+0.223606797749979*fUpwind_r[4]*alpha[8]; 
  Ghat_r[29] = 0.25*alpha[6]*fUpwind_r[47]+0.223606797749979*alpha[22]*fUpwind_r[43]+0.2500000000000001*alpha[3]*fUpwind_r[43]+0.2500000000000001*alpha[15]*fUpwind_r[42]+0.2500000000000001*alpha[1]*fUpwind_r[41]+0.2*alpha[16]*fUpwind_r[36]+0.223606797749979*alpha[25]*fUpwind_r[35]+0.25*alpha[7]*fUpwind_r[30]+0.223606797749979*alpha[12]*fUpwind_r[29]+0.25*alpha[0]*fUpwind_r[29]+0.25*alpha[5]*fUpwind_r[28]+0.2*alpha[9]*fUpwind_r[26]+0.223606797749979*alpha[8]*fUpwind_r[16]+0.223606797749979*fUpwind_r[8]*alpha[16]+0.2500000000000001*alpha[2]*fUpwind_r[14]+0.223606797749979*alpha[4]*fUpwind_r[9]+0.223606797749979*fUpwind_r[4]*alpha[9]; 
  Ghat_r[30] = 0.25*alpha[5]*fUpwind_r[47]+0.2500000000000001*alpha[2]*fUpwind_r[43]+0.2500000000000001*alpha[1]*fUpwind_r[42]+0.2500000000000001*alpha[15]*fUpwind_r[41]+0.223606797749979*alpha[25]*fUpwind_r[37]+0.223606797749979*alpha[16]*fUpwind_r[31]+0.25*alpha[0]*fUpwind_r[30]+0.25*alpha[7]*fUpwind_r[29]+0.25*alpha[6]*fUpwind_r[28]+0.223606797749979*alpha[9]*fUpwind_r[18]+0.223606797749979*alpha[8]*fUpwind_r[17]+0.2500000000000001*alpha[3]*fUpwind_r[14]+0.223606797749979*alpha[4]*fUpwind_r[10]; 
  Ghat_r[31] = 0.2*alpha[25]*fUpwind_r[47]+0.223606797749979*alpha[4]*fUpwind_r[47]+0.2*alpha[22]*fUpwind_r[46]+0.2*alpha[21]*fUpwind_r[46]+0.223606797749979*alpha[3]*fUpwind_r[46]+0.223606797749979*alpha[2]*fUpwind_r[45]+0.223606797749979*alpha[1]*fUpwind_r[44]+0.223606797749979*alpha[8]*fUpwind_r[43]+0.223606797749979*alpha[9]*fUpwind_r[42]+0.223606797749979*alpha[6]*fUpwind_r[40]+0.223606797749979*alpha[7]*fUpwind_r[39]+0.223606797749979*alpha[5]*fUpwind_r[38]+0.223606797749979*alpha[5]*fUpwind_r[37]+0.223606797749979*alpha[7]*fUpwind_r[36]+0.223606797749979*alpha[6]*fUpwind_r[35]+0.223606797749979*alpha[9]*fUpwind_r[33]+0.223606797749979*alpha[8]*fUpwind_r[32]+0.223606797749979*alpha[12]*fUpwind_r[31]+0.223606797749979*alpha[11]*fUpwind_r[31]+0.25*alpha[0]*fUpwind_r[31]+0.223606797749979*alpha[16]*fUpwind_r[30]+0.223606797749979*alpha[15]*fUpwind_r[27]+0.223606797749979*alpha[15]*fUpwind_r[26]+0.223606797749979*alpha[15]*fUpwind_r[25]+0.223606797749979*fUpwind_r[15]*alpha[25]+0.223606797749979*alpha[16]*fUpwind_r[22]+0.223606797749979*fUpwind_r[16]*alpha[22]+0.223606797749979*alpha[16]*fUpwind_r[21]+0.223606797749979*fUpwind_r[16]*alpha[21]+0.25*alpha[1]*fUpwind_r[18]+0.25*alpha[2]*fUpwind_r[17]+0.25*alpha[3]*fUpwind_r[16]+0.25*fUpwind_r[3]*alpha[16]+0.25*alpha[4]*fUpwind_r[15]+0.25*fUpwind_r[4]*alpha[15]+0.25*alpha[5]*fUpwind_r[10]+0.25*alpha[6]*fUpwind_r[9]+0.25*fUpwind_r[6]*alpha[9]+0.25*alpha[7]*fUpwind_r[8]+0.25*fUpwind_r[7]*alpha[8]; 
  Ghat_r[32] = 0.2*alpha[16]*fUpwind_r[45]+0.159719141249985*alpha[25]*fUpwind_r[44]+0.2500000000000001*alpha[4]*fUpwind_r[44]+0.25*alpha[9]*fUpwind_r[37]+0.2*alpha[6]*fUpwind_r[34]+0.2*alpha[5]*fUpwind_r[33]+0.223606797749979*alpha[12]*fUpwind_r[32]+0.159719141249985*alpha[11]*fUpwind_r[32]+0.25*alpha[0]*fUpwind_r[32]+0.223606797749979*alpha[8]*fUpwind_r[31]+0.2500000000000001*fUpwind_r[18]*alpha[25]+0.223606797749979*alpha[21]*fUpwind_r[24]+0.2*alpha[15]*fUpwind_r[23]+0.223606797749979*fUpwind_r[19]*alpha[22]+0.2500000000000001*alpha[2]*fUpwind_r[21]+0.159719141249985*fUpwind_r[19]*alpha[21]+0.2500000000000001*fUpwind_r[2]*alpha[21]+0.2*alpha[15]*fUpwind_r[20]+0.2500000000000001*alpha[3]*fUpwind_r[19]+0.223606797749979*alpha[16]*fUpwind_r[17]+0.223606797749979*alpha[1]*fUpwind_r[15]+0.223606797749979*fUpwind_r[1]*alpha[15]+0.25*alpha[7]*fUpwind_r[11]+0.25*fUpwind_r[7]*alpha[11]+0.223606797749979*alpha[5]*fUpwind_r[6]+0.223606797749979*fUpwind_r[5]*alpha[6]; 
  Ghat_r[33] = 0.223606797749979*alpha[25]*fUpwind_r[45]+0.2500000000000001*alpha[4]*fUpwind_r[45]+0.2*alpha[16]*fUpwind_r[44]+0.25*alpha[8]*fUpwind_r[38]+0.2*alpha[7]*fUpwind_r[34]+0.159719141249985*alpha[12]*fUpwind_r[33]+0.223606797749979*alpha[11]*fUpwind_r[33]+0.25*alpha[0]*fUpwind_r[33]+0.2*alpha[5]*fUpwind_r[32]+0.223606797749979*alpha[9]*fUpwind_r[31]+0.2*alpha[15]*fUpwind_r[24]+0.223606797749979*alpha[22]*fUpwind_r[23]+0.2500000000000001*alpha[1]*fUpwind_r[22]+0.159719141249985*fUpwind_r[20]*alpha[22]+0.2500000000000001*fUpwind_r[1]*alpha[22]+0.223606797749979*fUpwind_r[20]*alpha[21]+0.2500000000000001*alpha[3]*fUpwind_r[20]+0.2*alpha[15]*fUpwind_r[19]+0.223606797749979*alpha[16]*fUpwind_r[18]+0.223606797749979*alpha[2]*fUpwind_r[15]+0.223606797749979*fUpwind_r[2]*alpha[15]+0.25*alpha[6]*fUpwind_r[12]+0.25*fUpwind_r[6]*alpha[12]+0.223606797749979*alpha[5]*fUpwind_r[7]+0.223606797749979*fUpwind_r[5]*alpha[7]; 
  Ghat_r[34] = 0.223606797749979*alpha[25]*fUpwind_r[46]+0.2500000000000001*alpha[4]*fUpwind_r[46]+0.25*alpha[8]*fUpwind_r[40]+0.25*alpha[9]*fUpwind_r[39]+0.223606797749979*alpha[12]*fUpwind_r[34]+0.223606797749979*alpha[11]*fUpwind_r[34]+0.25*alpha[0]*fUpwind_r[34]+0.2*alpha[7]*fUpwind_r[33]+0.2*alpha[6]*fUpwind_r[32]+0.2500000000000001*alpha[16]*fUpwind_r[27]+0.2500000000000001*alpha[1]*fUpwind_r[24]+0.2500000000000001*alpha[2]*fUpwind_r[23]+0.2*alpha[15]*fUpwind_r[22]+0.2*fUpwind_r[15]*alpha[22]+0.2*alpha[15]*fUpwind_r[21]+0.2*fUpwind_r[15]*alpha[21]+0.223606797749979*alpha[3]*fUpwind_r[15]+0.223606797749979*fUpwind_r[3]*alpha[15]+0.25*alpha[5]*fUpwind_r[13]+0.223606797749979*alpha[6]*fUpwind_r[7]+0.223606797749979*fUpwind_r[6]*alpha[7]; 
  Ghat_r[35] = 0.2*alpha[15]*fUpwind_r[45]+0.223606797749979*alpha[22]*fUpwind_r[44]+0.159719141249985*alpha[21]*fUpwind_r[44]+0.2500000000000001*alpha[3]*fUpwind_r[44]+0.2*alpha[8]*fUpwind_r[41]+0.25*alpha[7]*fUpwind_r[37]+0.2*alpha[5]*fUpwind_r[36]+0.223606797749979*alpha[12]*fUpwind_r[35]+0.159719141249985*alpha[11]*fUpwind_r[35]+0.25*alpha[0]*fUpwind_r[35]+0.223606797749979*alpha[6]*fUpwind_r[31]+0.223606797749979*alpha[25]*fUpwind_r[29]+0.2*alpha[16]*fUpwind_r[28]+0.2500000000000001*alpha[2]*fUpwind_r[25]+0.159719141249985*fUpwind_r[19]*alpha[25]+0.2500000000000001*fUpwind_r[2]*alpha[25]+0.2500000000000001*fUpwind_r[18]*alpha[21]+0.2*alpha[16]*fUpwind_r[20]+0.2500000000000001*alpha[4]*fUpwind_r[19]+0.223606797749979*alpha[15]*fUpwind_r[17]+0.223606797749979*alpha[1]*fUpwind_r[16]+0.223606797749979*fUpwind_r[1]*alpha[16]+0.25*alpha[9]*fUpwind_r[11]+0.25*fUpwind_r[9]*alpha[11]+0.223606797749979*alpha[5]*fUpwind_r[8]+0.223606797749979*fUpwind_r[5]*alpha[8]; 
  Ghat_r[36] = 0.159719141249985*alpha[22]*fUpwind_r[45]+0.223606797749979*alpha[21]*fUpwind_r[45]+0.2500000000000001*alpha[3]*fUpwind_r[45]+0.2*alpha[15]*fUpwind_r[44]+0.2*alpha[9]*fUpwind_r[41]+0.25*alpha[6]*fUpwind_r[38]+0.159719141249985*alpha[12]*fUpwind_r[36]+0.223606797749979*alpha[11]*fUpwind_r[36]+0.25*alpha[0]*fUpwind_r[36]+0.2*alpha[5]*fUpwind_r[35]+0.223606797749979*alpha[7]*fUpwind_r[31]+0.2*alpha[16]*fUpwind_r[29]+0.2500000000000001*alpha[1]*fUpwind_r[26]+0.223606797749979*fUpwind_r[20]*alpha[25]+0.2500000000000001*fUpwind_r[17]*alpha[22]+0.2500000000000001*alpha[4]*fUpwind_r[20]+0.2*alpha[16]*fUpwind_r[19]+0.223606797749979*alpha[15]*fUpwind_r[18]+0.223606797749979*alpha[2]*fUpwind_r[16]+0.223606797749979*fUpwind_r[2]*alpha[16]+0.25*alpha[8]*fUpwind_r[12]+0.25*fUpwind_r[8]*alpha[12]+0.223606797749979*alpha[5]*fUpwind_r[9]+0.223606797749979*fUpwind_r[5]*alpha[9]; 
  Ghat_r[37] = 0.2*alpha[16]*fUpwind_r[47]+0.2*alpha[15]*fUpwind_r[46]+0.2500000000000001*alpha[2]*fUpwind_r[44]+0.2*alpha[8]*fUpwind_r[42]+0.2*alpha[6]*fUpwind_r[39]+0.159719141249985*alpha[11]*fUpwind_r[37]+0.25*alpha[0]*fUpwind_r[37]+0.25*alpha[7]*fUpwind_r[35]+0.25*alpha[9]*fUpwind_r[32]+0.223606797749979*alpha[5]*fUpwind_r[31]+0.223606797749979*alpha[25]*fUpwind_r[30]+0.223606797749979*alpha[21]*fUpwind_r[27]+0.159719141249985*alpha[21]*fUpwind_r[25]+0.2500000000000001*alpha[3]*fUpwind_r[25]+0.159719141249985*fUpwind_r[21]*alpha[25]+0.2500000000000001*fUpwind_r[3]*alpha[25]+0.2500000000000001*alpha[4]*fUpwind_r[21]+0.2500000000000001*fUpwind_r[4]*alpha[21]+0.223606797749979*alpha[1]*fUpwind_r[17]+0.223606797749979*alpha[15]*fUpwind_r[16]+0.223606797749979*fUpwind_r[15]*alpha[16]+0.25*fUpwind_r[10]*alpha[11]+0.223606797749979*alpha[6]*fUpwind_r[8]+0.223606797749979*fUpwind_r[6]*alpha[8]; 
  Ghat_r[38] = 0.2*alpha[16]*fUpwind_r[47]+0.2*alpha[15]*fUpwind_r[46]+0.2500000000000001*alpha[1]*fUpwind_r[45]+0.2*alpha[9]*fUpwind_r[43]+0.2*alpha[7]*fUpwind_r[40]+0.159719141249985*alpha[12]*fUpwind_r[38]+0.25*alpha[0]*fUpwind_r[38]+0.25*alpha[6]*fUpwind_r[36]+0.25*alpha[8]*fUpwind_r[33]+0.223606797749979*alpha[5]*fUpwind_r[31]+0.223606797749979*alpha[22]*fUpwind_r[27]+0.159719141249985*alpha[22]*fUpwind_r[26]+0.2500000000000001*alpha[3]*fUpwind_r[26]+0.2500000000000001*alpha[4]*fUpwind_r[22]+0.2500000000000001*fUpwind_r[4]*alpha[22]+0.223606797749979*alpha[2]*fUpwind_r[18]+0.223606797749979*alpha[15]*fUpwind_r[16]+0.223606797749979*fUpwind_r[15]*alpha[16]+0.25*fUpwind_r[10]*alpha[12]+0.223606797749979*alpha[7]*fUpwind_r[9]+0.223606797749979*fUpwind_r[7]*alpha[9]; 
  Ghat_r[39] = 0.2500000000000001*alpha[2]*fUpwind_r[46]+0.223606797749979*alpha[22]*fUpwind_r[45]+0.2*alpha[15]*fUpwind_r[44]+0.25*alpha[5]*fUpwind_r[40]+0.223606797749979*alpha[11]*fUpwind_r[39]+0.25*alpha[0]*fUpwind_r[39]+0.2*alpha[6]*fUpwind_r[37]+0.25*alpha[9]*fUpwind_r[34]+0.223606797749979*alpha[7]*fUpwind_r[31]+0.2500000000000001*alpha[1]*fUpwind_r[27]+0.223606797749979*fUpwind_r[23]*alpha[25]+0.2500000000000001*alpha[16]*fUpwind_r[24]+0.2500000000000001*alpha[4]*fUpwind_r[23]+0.2*fUpwind_r[17]*alpha[21]+0.223606797749979*alpha[15]*fUpwind_r[18]+0.223606797749979*alpha[3]*fUpwind_r[17]+0.25*alpha[8]*fUpwind_r[13]+0.223606797749979*alpha[6]*fUpwind_r[10]; 
  Ghat_r[40] = 0.2500000000000001*alpha[1]*fUpwind_r[46]+0.2*alpha[15]*fUpwind_r[45]+0.223606797749979*alpha[21]*fUpwind_r[44]+0.223606797749979*alpha[12]*fUpwind_r[40]+0.25*alpha[0]*fUpwind_r[40]+0.25*alpha[5]*fUpwind_r[39]+0.2*alpha[7]*fUpwind_r[38]+0.25*alpha[8]*fUpwind_r[34]+0.223606797749979*alpha[6]*fUpwind_r[31]+0.2500000000000001*alpha[2]*fUpwind_r[27]+0.2500000000000001*alpha[4]*fUpwind_r[24]+0.2500000000000001*alpha[16]*fUpwind_r[23]+0.2*fUpwind_r[18]*alpha[22]+0.223606797749979*alpha[3]*fUpwind_r[18]+0.223606797749979*alpha[15]*fUpwind_r[17]+0.25*alpha[9]*fUpwind_r[13]+0.223606797749979*alpha[7]*fUpwind_r[10]; 
  Ghat_r[41] = 0.223606797749979*alpha[22]*fUpwind_r[47]+0.223606797749979*alpha[21]*fUpwind_r[47]+0.2500000000000001*alpha[3]*fUpwind_r[47]+0.25*alpha[6]*fUpwind_r[43]+0.25*alpha[7]*fUpwind_r[42]+0.223606797749979*alpha[12]*fUpwind_r[41]+0.223606797749979*alpha[11]*fUpwind_r[41]+0.25*alpha[0]*fUpwind_r[41]+0.2*alpha[9]*fUpwind_r[36]+0.2*alpha[8]*fUpwind_r[35]+0.2500000000000001*alpha[15]*fUpwind_r[30]+0.2500000000000001*alpha[1]*fUpwind_r[29]+0.2500000000000001*alpha[2]*fUpwind_r[28]+0.2*alpha[16]*fUpwind_r[26]+0.2*alpha[16]*fUpwind_r[25]+0.2*fUpwind_r[16]*alpha[25]+0.223606797749979*alpha[4]*fUpwind_r[16]+0.223606797749979*fUpwind_r[4]*alpha[16]+0.25*alpha[5]*fUpwind_r[14]+0.223606797749979*alpha[8]*fUpwind_r[9]+0.223606797749979*fUpwind_r[8]*alpha[9]; 
  Ghat_r[42] = 0.2500000000000001*alpha[2]*fUpwind_r[47]+0.2*alpha[16]*fUpwind_r[44]+0.25*alpha[5]*fUpwind_r[43]+0.223606797749979*alpha[11]*fUpwind_r[42]+0.25*alpha[0]*fUpwind_r[42]+0.25*alpha[7]*fUpwind_r[41]+0.2*alpha[8]*fUpwind_r[37]+0.223606797749979*alpha[9]*fUpwind_r[31]+0.2500000000000001*alpha[1]*fUpwind_r[30]+0.2500000000000001*alpha[15]*fUpwind_r[29]+0.223606797749979*alpha[21]*fUpwind_r[28]+0.2500000000000001*alpha[3]*fUpwind_r[28]+0.2*fUpwind_r[17]*alpha[25]+0.223606797749979*alpha[16]*fUpwind_r[18]+0.223606797749979*alpha[4]*fUpwind_r[17]+0.25*alpha[6]*fUpwind_r[14]+0.223606797749979*alpha[8]*fUpwind_r[10]; 
  Ghat_r[43] = 0.2500000000000001*alpha[1]*fUpwind_r[47]+0.2*alpha[16]*fUpwind_r[45]+0.223606797749979*alpha[25]*fUpwind_r[44]+0.223606797749979*alpha[12]*fUpwind_r[43]+0.25*alpha[0]*fUpwind_r[43]+0.25*alpha[5]*fUpwind_r[42]+0.25*alpha[6]*fUpwind_r[41]+0.2*alpha[9]*fUpwind_r[38]+0.223606797749979*alpha[8]*fUpwind_r[31]+0.2500000000000001*alpha[2]*fUpwind_r[30]+0.223606797749979*alpha[22]*fUpwind_r[29]+0.2500000000000001*alpha[3]*fUpwind_r[29]+0.2500000000000001*alpha[15]*fUpwind_r[28]+0.223606797749979*alpha[4]*fUpwind_r[18]+0.223606797749979*alpha[16]*fUpwind_r[17]+0.25*alpha[7]*fUpwind_r[14]+0.223606797749979*alpha[9]*fUpwind_r[10]; 
  Ghat_r[44] = 0.2*alpha[8]*fUpwind_r[47]+0.2*alpha[6]*fUpwind_r[46]+0.2*alpha[5]*fUpwind_r[45]+0.223606797749979*alpha[12]*fUpwind_r[44]+0.159719141249985*alpha[11]*fUpwind_r[44]+0.25*alpha[0]*fUpwind_r[44]+0.223606797749979*alpha[25]*fUpwind_r[43]+0.2*alpha[16]*fUpwind_r[42]+0.223606797749979*alpha[21]*fUpwind_r[40]+0.2*alpha[15]*fUpwind_r[39]+0.2500000000000001*alpha[2]*fUpwind_r[37]+0.2*alpha[15]*fUpwind_r[36]+0.223606797749979*alpha[22]*fUpwind_r[35]+0.159719141249985*alpha[21]*fUpwind_r[35]+0.2500000000000001*alpha[3]*fUpwind_r[35]+0.2*alpha[16]*fUpwind_r[33]+0.159719141249985*alpha[25]*fUpwind_r[32]+0.2500000000000001*alpha[4]*fUpwind_r[32]+0.223606797749979*alpha[1]*fUpwind_r[31]+0.25*alpha[7]*fUpwind_r[25]+0.25*fUpwind_r[7]*alpha[25]+0.25*alpha[9]*fUpwind_r[21]+0.25*fUpwind_r[9]*alpha[21]+0.2500000000000001*alpha[11]*fUpwind_r[18]+0.223606797749979*alpha[5]*fUpwind_r[17]+0.223606797749979*alpha[6]*fUpwind_r[16]+0.223606797749979*fUpwind_r[6]*alpha[16]+0.223606797749979*alpha[8]*fUpwind_r[15]+0.223606797749979*fUpwind_r[8]*alpha[15]; 
  Ghat_r[45] = 0.2*alpha[9]*fUpwind_r[47]+0.2*alpha[7]*fUpwind_r[46]+0.159719141249985*alpha[12]*fUpwind_r[45]+0.223606797749979*alpha[11]*fUpwind_r[45]+0.25*alpha[0]*fUpwind_r[45]+0.2*alpha[5]*fUpwind_r[44]+0.2*alpha[16]*fUpwind_r[43]+0.2*alpha[15]*fUpwind_r[40]+0.223606797749979*alpha[22]*fUpwind_r[39]+0.2500000000000001*alpha[1]*fUpwind_r[38]+0.159719141249985*alpha[22]*fUpwind_r[36]+0.223606797749979*alpha[21]*fUpwind_r[36]+0.2500000000000001*alpha[3]*fUpwind_r[36]+0.2*alpha[15]*fUpwind_r[35]+0.223606797749979*alpha[25]*fUpwind_r[33]+0.2500000000000001*alpha[4]*fUpwind_r[33]+0.2*alpha[16]*fUpwind_r[32]+0.223606797749979*alpha[2]*fUpwind_r[31]+0.25*alpha[6]*fUpwind_r[26]+0.25*alpha[8]*fUpwind_r[22]+0.25*fUpwind_r[8]*alpha[22]+0.223606797749979*alpha[5]*fUpwind_r[18]+0.2500000000000001*alpha[12]*fUpwind_r[17]+0.223606797749979*alpha[7]*fUpwind_r[16]+0.223606797749979*fUpwind_r[7]*alpha[16]+0.223606797749979*alpha[9]*fUpwind_r[15]+0.223606797749979*fUpwind_r[9]*alpha[15]; 
  Ghat_r[46] = 0.223606797749979*alpha[12]*fUpwind_r[46]+0.223606797749979*alpha[11]*fUpwind_r[46]+0.25*alpha[0]*fUpwind_r[46]+0.2*alpha[7]*fUpwind_r[45]+0.2*alpha[6]*fUpwind_r[44]+0.2500000000000001*alpha[1]*fUpwind_r[40]+0.2500000000000001*alpha[2]*fUpwind_r[39]+0.2*alpha[15]*fUpwind_r[38]+0.2*alpha[15]*fUpwind_r[37]+0.223606797749979*alpha[25]*fUpwind_r[34]+0.2500000000000001*alpha[4]*fUpwind_r[34]+0.2*alpha[22]*fUpwind_r[31]+0.2*alpha[21]*fUpwind_r[31]+0.223606797749979*alpha[3]*fUpwind_r[31]+0.25*alpha[5]*fUpwind_r[27]+0.25*alpha[8]*fUpwind_r[24]+0.25*alpha[9]*fUpwind_r[23]+0.223606797749979*alpha[6]*fUpwind_r[18]+0.223606797749979*alpha[7]*fUpwind_r[17]+0.2500000000000001*fUpwind_r[13]*alpha[16]+0.223606797749979*fUpwind_r[10]*alpha[15]; 
  Ghat_r[47] = 0.223606797749979*alpha[12]*fUpwind_r[47]+0.223606797749979*alpha[11]*fUpwind_r[47]+0.25*alpha[0]*fUpwind_r[47]+0.2*alpha[9]*fUpwind_r[45]+0.2*alpha[8]*fUpwind_r[44]+0.2500000000000001*alpha[1]*fUpwind_r[43]+0.2500000000000001*alpha[2]*fUpwind_r[42]+0.223606797749979*alpha[22]*fUpwind_r[41]+0.223606797749979*alpha[21]*fUpwind_r[41]+0.2500000000000001*alpha[3]*fUpwind_r[41]+0.2*alpha[16]*fUpwind_r[38]+0.2*alpha[16]*fUpwind_r[37]+0.2*alpha[25]*fUpwind_r[31]+0.223606797749979*alpha[4]*fUpwind_r[31]+0.25*alpha[5]*fUpwind_r[30]+0.25*alpha[6]*fUpwind_r[29]+0.25*alpha[7]*fUpwind_r[28]+0.223606797749979*alpha[8]*fUpwind_r[18]+0.223606797749979*alpha[9]*fUpwind_r[17]+0.223606797749979*fUpwind_r[10]*alpha[16]+0.2500000000000001*fUpwind_r[14]*alpha[15]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv11; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv11; 
  out[2] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv11; 
  out[3] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv11; 
  out[4] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv11; 
  out[5] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dv11; 
  out[6] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dv11; 
  out[7] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dv11; 
  out[8] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dv11; 
  out[9] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv11; 
  out[10] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv11; 
  out[11] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv11; 
  out[12] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*dv11; 
  out[13] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*dv11; 
  out[14] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*dv11; 
  out[15] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dv11; 
  out[16] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*dv11; 
  out[17] += (0.7071067811865475*Ghat_l[12]-0.7071067811865475*Ghat_r[12])*dv11; 
  out[18] += (0.7071067811865475*Ghat_l[13]-0.7071067811865475*Ghat_r[13])*dv11; 
  out[19] += (1.58113883008419*Ghat_l[0]-1.58113883008419*Ghat_r[0])*dv11; 
  out[20] += (0.7071067811865475*Ghat_l[14]-0.7071067811865475*Ghat_r[14])*dv11; 
  out[21] += (0.7071067811865475*Ghat_l[15]-0.7071067811865475*Ghat_r[15])*dv11; 
  out[22] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dv11; 
  out[23] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dv11; 
  out[24] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dv11; 
  out[25] += (0.7071067811865475*Ghat_l[16]-0.7071067811865475*Ghat_r[16])*dv11; 
  out[26] += (0.7071067811865475*Ghat_l[17]-0.7071067811865475*Ghat_r[17])*dv11; 
  out[27] += (0.7071067811865475*Ghat_l[18]-0.7071067811865475*Ghat_r[18])*dv11; 
  out[28] += -1.224744871391589*(Ghat_r[8]+Ghat_l[8])*dv11; 
  out[29] += -1.224744871391589*(Ghat_r[9]+Ghat_l[9])*dv11; 
  out[30] += -1.224744871391589*(Ghat_r[10]+Ghat_l[10])*dv11; 
  out[31] += (0.7071067811865475*Ghat_l[19]-0.7071067811865475*Ghat_r[19])*dv11; 
  out[32] += (0.7071067811865475*Ghat_l[20]-0.7071067811865475*Ghat_r[20])*dv11; 
  out[33] += (0.7071067811865475*Ghat_l[21]-0.7071067811865475*Ghat_r[21])*dv11; 
  out[34] += (0.7071067811865475*Ghat_l[22]-0.7071067811865475*Ghat_r[22])*dv11; 
  out[35] += (0.7071067811865475*Ghat_l[23]-0.7071067811865475*Ghat_r[23])*dv11; 
  out[36] += (0.7071067811865475*Ghat_l[24]-0.7071067811865475*Ghat_r[24])*dv11; 
  out[37] += -1.224744871391589*(Ghat_r[11]+Ghat_l[11])*dv11; 
  out[38] += -1.224744871391589*(Ghat_r[12]+Ghat_l[12])*dv11; 
  out[39] += -1.224744871391589*(Ghat_r[13]+Ghat_l[13])*dv11; 
  out[40] += (1.58113883008419*Ghat_l[1]-1.58113883008419*Ghat_r[1])*dv11; 
  out[41] += (1.58113883008419*Ghat_l[2]-1.58113883008419*Ghat_r[2])*dv11; 
  out[42] += (1.58113883008419*Ghat_l[3]-1.58113883008419*Ghat_r[3])*dv11; 
  out[43] += (0.7071067811865475*Ghat_l[25]-0.7071067811865475*Ghat_r[25])*dv11; 
  out[44] += (0.7071067811865475*Ghat_l[26]-0.7071067811865475*Ghat_r[26])*dv11; 
  out[45] += (0.7071067811865475*Ghat_l[27]-0.7071067811865475*Ghat_r[27])*dv11; 
  out[46] += (1.58113883008419*Ghat_l[4]-1.58113883008419*Ghat_r[4])*dv11; 
  out[47] += (0.7071067811865475*Ghat_l[28]-0.7071067811865475*Ghat_r[28])*dv11; 
  out[48] += (0.7071067811865475*Ghat_l[29]-0.7071067811865475*Ghat_r[29])*dv11; 
  out[49] += (0.7071067811865475*Ghat_l[30]-0.7071067811865475*Ghat_r[30])*dv11; 
  out[50] += -1.224744871391589*(Ghat_r[14]+Ghat_l[14])*dv11; 
  out[51] += -1.224744871391589*(Ghat_r[15]+Ghat_l[15])*dv11; 
  out[52] += (0.7071067811865475*Ghat_l[31]-0.7071067811865475*Ghat_r[31])*dv11; 
  out[53] += -1.224744871391589*(Ghat_r[16]+Ghat_l[16])*dv11; 
  out[54] += -1.224744871391589*(Ghat_r[17]+Ghat_l[17])*dv11; 
  out[55] += -1.224744871391589*(Ghat_r[18]+Ghat_l[18])*dv11; 
  out[56] += (0.7071067811865475*Ghat_l[32]-0.7071067811865475*Ghat_r[32])*dv11; 
  out[57] += (0.7071067811865475*Ghat_l[33]-0.7071067811865475*Ghat_r[33])*dv11; 
  out[58] += (0.7071067811865475*Ghat_l[34]-0.7071067811865475*Ghat_r[34])*dv11; 
  out[59] += -1.224744871391589*(Ghat_r[19]+Ghat_l[19])*dv11; 
  out[60] += -1.224744871391589*(Ghat_r[20]+Ghat_l[20])*dv11; 
  out[61] += -1.224744871391589*(Ghat_r[21]+Ghat_l[21])*dv11; 
  out[62] += -1.224744871391589*(Ghat_r[22]+Ghat_l[22])*dv11; 
  out[63] += -1.224744871391589*(Ghat_r[23]+Ghat_l[23])*dv11; 
  out[64] += -1.224744871391589*(Ghat_r[24]+Ghat_l[24])*dv11; 
  out[65] += (1.58113883008419*Ghat_l[5]-1.58113883008419*Ghat_r[5])*dv11; 
  out[66] += (1.58113883008419*Ghat_l[6]-1.58113883008419*Ghat_r[6])*dv11; 
  out[67] += (1.58113883008419*Ghat_l[7]-1.58113883008419*Ghat_r[7])*dv11; 
  out[68] += (0.7071067811865475*Ghat_l[35]-0.7071067811865475*Ghat_r[35])*dv11; 
  out[69] += (0.7071067811865475*Ghat_l[36]-0.7071067811865475*Ghat_r[36])*dv11; 
  out[70] += (0.7071067811865475*Ghat_l[37]-0.7071067811865475*Ghat_r[37])*dv11; 
  out[71] += (0.7071067811865475*Ghat_l[38]-0.7071067811865475*Ghat_r[38])*dv11; 
  out[72] += (0.7071067811865475*Ghat_l[39]-0.7071067811865475*Ghat_r[39])*dv11; 
  out[73] += (0.7071067811865475*Ghat_l[40]-0.7071067811865475*Ghat_r[40])*dv11; 
  out[74] += -1.224744871391589*(Ghat_r[25]+Ghat_l[25])*dv11; 
  out[75] += -1.224744871391589*(Ghat_r[26]+Ghat_l[26])*dv11; 
  out[76] += -1.224744871391589*(Ghat_r[27]+Ghat_l[27])*dv11; 
  out[77] += (1.58113883008419*Ghat_l[8]-1.58113883008419*Ghat_r[8])*dv11; 
  out[78] += (1.58113883008419*Ghat_l[9]-1.58113883008419*Ghat_r[9])*dv11; 
  out[79] += (1.58113883008419*Ghat_l[10]-1.58113883008419*Ghat_r[10])*dv11; 
  out[80] += (0.7071067811865475*Ghat_l[41]-0.7071067811865475*Ghat_r[41])*dv11; 
  out[81] += (0.7071067811865475*Ghat_l[42]-0.7071067811865475*Ghat_r[42])*dv11; 
  out[82] += (0.7071067811865475*Ghat_l[43]-0.7071067811865475*Ghat_r[43])*dv11; 
  out[83] += -1.224744871391589*(Ghat_r[28]+Ghat_l[28])*dv11; 
  out[84] += -1.224744871391589*(Ghat_r[29]+Ghat_l[29])*dv11; 
  out[85] += -1.224744871391589*(Ghat_r[30]+Ghat_l[30])*dv11; 
  out[86] += -1.224744871391589*(Ghat_r[31]+Ghat_l[31])*dv11; 
  out[87] += -1.224744871391589*(Ghat_r[32]+Ghat_l[32])*dv11; 
  out[88] += -1.224744871391589*(Ghat_r[33]+Ghat_l[33])*dv11; 
  out[89] += -1.224744871391589*(Ghat_r[34]+Ghat_l[34])*dv11; 
  out[90] += (1.58113883008419*Ghat_l[15]-1.58113883008419*Ghat_r[15])*dv11; 
  out[91] += (0.7071067811865475*Ghat_l[44]-0.7071067811865475*Ghat_r[44])*dv11; 
  out[92] += (0.7071067811865475*Ghat_l[45]-0.7071067811865475*Ghat_r[45])*dv11; 
  out[93] += (0.7071067811865475*Ghat_l[46]-0.7071067811865475*Ghat_r[46])*dv11; 
  out[94] += -1.224744871391589*(Ghat_r[35]+Ghat_l[35])*dv11; 
  out[95] += -1.224744871391589*(Ghat_r[36]+Ghat_l[36])*dv11; 
  out[96] += -1.224744871391589*(Ghat_r[37]+Ghat_l[37])*dv11; 
  out[97] += -1.224744871391589*(Ghat_r[38]+Ghat_l[38])*dv11; 
  out[98] += -1.224744871391589*(Ghat_r[39]+Ghat_l[39])*dv11; 
  out[99] += -1.224744871391589*(Ghat_r[40]+Ghat_l[40])*dv11; 
  out[100] += (1.58113883008419*Ghat_l[16]-1.58113883008419*Ghat_r[16])*dv11; 
  out[101] += (1.58113883008419*Ghat_l[17]-1.58113883008419*Ghat_r[17])*dv11; 
  out[102] += (1.58113883008419*Ghat_l[18]-1.58113883008419*Ghat_r[18])*dv11; 
  out[103] += (0.7071067811865475*Ghat_l[47]-0.7071067811865475*Ghat_r[47])*dv11; 
  out[104] += -1.224744871391589*(Ghat_r[41]+Ghat_l[41])*dv11; 
  out[105] += -1.224744871391589*(Ghat_r[42]+Ghat_l[42])*dv11; 
  out[106] += -1.224744871391589*(Ghat_r[43]+Ghat_l[43])*dv11; 
  out[107] += -1.224744871391589*(Ghat_r[44]+Ghat_l[44])*dv11; 
  out[108] += -1.224744871391589*(Ghat_r[45]+Ghat_l[45])*dv11; 
  out[109] += -1.224744871391589*(Ghat_r[46]+Ghat_l[46])*dv11; 
  out[110] += (1.58113883008419*Ghat_l[31]-1.58113883008419*Ghat_r[31])*dv11; 
  out[111] += -1.224744871391589*(Ghat_r[47]+Ghat_l[47])*dv11; 

  return 0.;

} 
