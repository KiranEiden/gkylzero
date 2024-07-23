#include <gkyl_vlasov_sr_kernels.h> 
#include <gkyl_basis_hyb_1x3v_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_hyb_1x3v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_sr_surfvx_1x3v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // gamma:     Particle Lorentz boost factor sqrt(1 + p^2).
  // qmem:      q/m*EM fields.
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  const double dv12 = 2.0/dxv[3]; 
  const double *E0 = &qmem[0]; 
  double p1_over_gamma_l[8] = {0.0}; 
  double p1_over_gamma_r[8] = {0.0}; 
  p1_over_gamma_l[0] = 2.738612787525831*gamma[11]*dv11-2.121320343559642*gamma[4]*dv11+1.224744871391589*gamma[2]*dv11; 
  p1_over_gamma_l[1] = 2.738612787525831*gamma[8]*dv11-4.743416490252569*gamma[12]*dv11; 
  p1_over_gamma_l[2] = 2.738612787525831*gamma[17]*dv11-2.121320343559642*gamma[10]*dv11+1.224744871391589*gamma[6]*dv11; 
  p1_over_gamma_l[3] = 2.738612787525831*gamma[14]*dv11-4.743416490252569*gamma[18]*dv11; 
  p1_over_gamma_l[5] = 1.224744871391589*gamma[16]*dv11-2.121320343559642*gamma[19]*dv11; 
  p1_over_gamma_r[0] = 2.738612787525831*gamma[11]*dv11+2.121320343559642*gamma[4]*dv11+1.224744871391589*gamma[2]*dv11; 
  p1_over_gamma_r[1] = 4.743416490252569*gamma[12]*dv11+2.738612787525831*gamma[8]*dv11; 
  p1_over_gamma_r[2] = 2.738612787525831*gamma[17]*dv11+2.121320343559642*gamma[10]*dv11+1.224744871391589*gamma[6]*dv11; 
  p1_over_gamma_r[3] = 4.743416490252569*gamma[18]*dv11+2.738612787525831*gamma[14]*dv11; 
  p1_over_gamma_r[5] = 2.121320343559642*gamma[19]*dv11+1.224744871391589*gamma[16]*dv11; 

  double p2_over_gamma_l[8] = {0.0}; 
  double p2_over_gamma_r[8] = {0.0}; 
  p2_over_gamma_l[0] = 2.738612787525831*gamma[13]*dv12-2.121320343559642*gamma[5]*dv12+1.224744871391589*gamma[3]*dv12; 
  p2_over_gamma_l[1] = 2.738612787525831*gamma[17]*dv12-2.121320343559642*gamma[10]*dv12+1.224744871391589*gamma[6]*dv12; 
  p2_over_gamma_l[2] = 2.738612787525831*gamma[9]*dv12-4.743416490252569*gamma[15]*dv12; 
  p2_over_gamma_l[3] = 2.738612787525831*gamma[16]*dv12-4.743416490252569*gamma[19]*dv12; 
  p2_over_gamma_l[4] = 1.224744871391589*gamma[14]*dv12-2.121320343559642*gamma[18]*dv12; 
  p2_over_gamma_r[0] = 2.738612787525831*gamma[13]*dv12+2.121320343559642*gamma[5]*dv12+1.224744871391589*gamma[3]*dv12; 
  p2_over_gamma_r[1] = 2.738612787525831*gamma[17]*dv12+2.121320343559642*gamma[10]*dv12+1.224744871391589*gamma[6]*dv12; 
  p2_over_gamma_r[2] = 4.743416490252569*gamma[15]*dv12+2.738612787525831*gamma[9]*dv12; 
  p2_over_gamma_r[3] = 4.743416490252569*gamma[19]*dv12+2.738612787525831*gamma[16]*dv12; 
  p2_over_gamma_r[4] = 2.121320343559642*gamma[18]*dv12+1.224744871391589*gamma[14]*dv12; 

  const double *B1 = &qmem[8]; 
  const double *B2 = &qmem[10]; 

  double alpha_l[16] = {0.0}; 
  double alpha_r[16] = {0.0}; 

  alpha_l[0] = (-1.0*B1[0]*p2_over_gamma_l[0])+B2[0]*p1_over_gamma_l[0]+2.0*E0[0]; 
  alpha_l[1] = 2.0*E0[1]+p1_over_gamma_l[0]*B2[1]-1.0*p2_over_gamma_l[0]*B1[1]; 
  alpha_l[2] = B2[0]*p1_over_gamma_l[1]-1.0*B1[0]*p2_over_gamma_l[1]; 
  alpha_l[3] = B2[0]*p1_over_gamma_l[2]-1.0*B1[0]*p2_over_gamma_l[2]; 
  alpha_l[4] = B2[1]*p1_over_gamma_l[1]-1.0*B1[1]*p2_over_gamma_l[1]; 
  alpha_l[5] = B2[1]*p1_over_gamma_l[2]-1.0*B1[1]*p2_over_gamma_l[2]; 
  alpha_l[6] = B2[0]*p1_over_gamma_l[3]-1.0*B1[0]*p2_over_gamma_l[3]; 
  alpha_l[7] = B2[1]*p1_over_gamma_l[3]-1.0*B1[1]*p2_over_gamma_l[3]; 
  alpha_l[8] = -1.0*B1[0]*p2_over_gamma_l[4]; 
  alpha_l[9] = -1.0*B1[1]*p2_over_gamma_l[4]; 
  alpha_l[12] = B2[0]*p1_over_gamma_l[5]; 
  alpha_l[13] = 1.0*B2[1]*p1_over_gamma_l[5]; 

  alpha_r[0] = (-1.0*B1[0]*p2_over_gamma_r[0])+B2[0]*p1_over_gamma_r[0]+2.0*E0[0]; 
  alpha_r[1] = 2.0*E0[1]+p1_over_gamma_r[0]*B2[1]-1.0*p2_over_gamma_r[0]*B1[1]; 
  alpha_r[2] = B2[0]*p1_over_gamma_r[1]-1.0*B1[0]*p2_over_gamma_r[1]; 
  alpha_r[3] = B2[0]*p1_over_gamma_r[2]-1.0*B1[0]*p2_over_gamma_r[2]; 
  alpha_r[4] = B2[1]*p1_over_gamma_r[1]-1.0*B1[1]*p2_over_gamma_r[1]; 
  alpha_r[5] = B2[1]*p1_over_gamma_r[2]-1.0*B1[1]*p2_over_gamma_r[2]; 
  alpha_r[6] = B2[0]*p1_over_gamma_r[3]-1.0*B1[0]*p2_over_gamma_r[3]; 
  alpha_r[7] = B2[1]*p1_over_gamma_r[3]-1.0*B1[1]*p2_over_gamma_r[3]; 
  alpha_r[8] = -1.0*B1[0]*p2_over_gamma_r[4]; 
  alpha_r[9] = -1.0*B1[1]*p2_over_gamma_r[4]; 
  alpha_r[12] = B2[0]*p1_over_gamma_r[5]; 
  alpha_r[13] = 1.0*B2[1]*p1_over_gamma_r[5]; 

  double fUpwindQuad_l[18] = {0.0};
  double fUpwindQuad_r[18] = {0.0};
  double fUpwind_l[16] = {0.0};;
  double fUpwind_r[16] = {0.0};
  double Ghat_l[16] = {0.0}; 
  double Ghat_r[16] = {0.0}; 

  if ((-0.3162277660168379*alpha_l[13])+0.3162277660168379*alpha_l[12]-0.3162277660168379*alpha_l[9]+0.3162277660168379*alpha_l[8]-0.6363961030678926*alpha_l[7]+0.6363961030678926*alpha_l[6]+0.4743416490252568*(alpha_l[5]+alpha_l[4])-0.4743416490252568*(alpha_l[3]+alpha_l[2])-0.3535533905932737*alpha_l[1]+0.3535533905932737*alpha_l[0] > 0) { 
    fUpwindQuad_l[0] = hyb_1x3v_p1_surfx2_eval_quad_node_0_r(fl); 
  } else { 
    fUpwindQuad_l[0] = hyb_1x3v_p1_surfx2_eval_quad_node_0_l(fc); 
  } 
  if ((-0.3162277660168379*alpha_r[13])+0.3162277660168379*alpha_r[12]-0.3162277660168379*alpha_r[9]+0.3162277660168379*alpha_r[8]-0.6363961030678926*alpha_r[7]+0.6363961030678926*alpha_r[6]+0.4743416490252568*(alpha_r[5]+alpha_r[4])-0.4743416490252568*(alpha_r[3]+alpha_r[2])-0.3535533905932737*alpha_r[1]+0.3535533905932737*alpha_r[0] > 0) { 
    fUpwindQuad_r[0] = hyb_1x3v_p1_surfx2_eval_quad_node_0_r(fc); 
  } else { 
    fUpwindQuad_r[0] = hyb_1x3v_p1_surfx2_eval_quad_node_0_l(fr); 
  } 
  if (0.3952847075210473*alpha_l[13]-0.3952847075210473*alpha_l[12]-0.3162277660168379*alpha_l[9]+0.3162277660168379*alpha_l[8]+0.4743416490252568*alpha_l[4]-0.4743416490252568*alpha_l[2]-0.3535533905932737*alpha_l[1]+0.3535533905932737*alpha_l[0] > 0) { 
    fUpwindQuad_l[1] = hyb_1x3v_p1_surfx2_eval_quad_node_1_r(fl); 
  } else { 
    fUpwindQuad_l[1] = hyb_1x3v_p1_surfx2_eval_quad_node_1_l(fc); 
  } 
  if (0.3952847075210473*alpha_r[13]-0.3952847075210473*alpha_r[12]-0.3162277660168379*alpha_r[9]+0.3162277660168379*alpha_r[8]+0.4743416490252568*alpha_r[4]-0.4743416490252568*alpha_r[2]-0.3535533905932737*alpha_r[1]+0.3535533905932737*alpha_r[0] > 0) { 
    fUpwindQuad_r[1] = hyb_1x3v_p1_surfx2_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_r[1] = hyb_1x3v_p1_surfx2_eval_quad_node_1_l(fr); 
  } 
  if ((-0.3162277660168379*alpha_l[13])+0.3162277660168379*alpha_l[12]-0.3162277660168379*alpha_l[9]+0.3162277660168379*alpha_l[8]+0.6363961030678926*alpha_l[7]-0.6363961030678926*alpha_l[6]-0.4743416490252568*alpha_l[5]+0.4743416490252568*(alpha_l[4]+alpha_l[3])-0.4743416490252568*alpha_l[2]-0.3535533905932737*alpha_l[1]+0.3535533905932737*alpha_l[0] > 0) { 
    fUpwindQuad_l[2] = hyb_1x3v_p1_surfx2_eval_quad_node_2_r(fl); 
  } else { 
    fUpwindQuad_l[2] = hyb_1x3v_p1_surfx2_eval_quad_node_2_l(fc); 
  } 
  if ((-0.3162277660168379*alpha_r[13])+0.3162277660168379*alpha_r[12]-0.3162277660168379*alpha_r[9]+0.3162277660168379*alpha_r[8]+0.6363961030678926*alpha_r[7]-0.6363961030678926*alpha_r[6]-0.4743416490252568*alpha_r[5]+0.4743416490252568*(alpha_r[4]+alpha_r[3])-0.4743416490252568*alpha_r[2]-0.3535533905932737*alpha_r[1]+0.3535533905932737*alpha_r[0] > 0) { 
    fUpwindQuad_r[2] = hyb_1x3v_p1_surfx2_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_r[2] = hyb_1x3v_p1_surfx2_eval_quad_node_2_l(fr); 
  } 
  if ((-0.3162277660168379*alpha_l[13])+0.3162277660168379*alpha_l[12]+0.3952847075210473*alpha_l[9]-0.3952847075210473*alpha_l[8]+0.4743416490252568*alpha_l[5]-0.4743416490252568*alpha_l[3]-0.3535533905932737*alpha_l[1]+0.3535533905932737*alpha_l[0] > 0) { 
    fUpwindQuad_l[3] = hyb_1x3v_p1_surfx2_eval_quad_node_3_r(fl); 
  } else { 
    fUpwindQuad_l[3] = hyb_1x3v_p1_surfx2_eval_quad_node_3_l(fc); 
  } 
  if ((-0.3162277660168379*alpha_r[13])+0.3162277660168379*alpha_r[12]+0.3952847075210473*alpha_r[9]-0.3952847075210473*alpha_r[8]+0.4743416490252568*alpha_r[5]-0.4743416490252568*alpha_r[3]-0.3535533905932737*alpha_r[1]+0.3535533905932737*alpha_r[0] > 0) { 
    fUpwindQuad_r[3] = hyb_1x3v_p1_surfx2_eval_quad_node_3_r(fc); 
  } else { 
    fUpwindQuad_r[3] = hyb_1x3v_p1_surfx2_eval_quad_node_3_l(fr); 
  } 
  if (0.3952847075210473*alpha_l[13]-0.3952847075210473*alpha_l[12]+0.3952847075210473*alpha_l[9]-0.3952847075210473*alpha_l[8]-0.3535533905932737*alpha_l[1]+0.3535533905932737*alpha_l[0] > 0) { 
    fUpwindQuad_l[4] = hyb_1x3v_p1_surfx2_eval_quad_node_4_r(fl); 
  } else { 
    fUpwindQuad_l[4] = hyb_1x3v_p1_surfx2_eval_quad_node_4_l(fc); 
  } 
  if (0.3952847075210473*alpha_r[13]-0.3952847075210473*alpha_r[12]+0.3952847075210473*alpha_r[9]-0.3952847075210473*alpha_r[8]-0.3535533905932737*alpha_r[1]+0.3535533905932737*alpha_r[0] > 0) { 
    fUpwindQuad_r[4] = hyb_1x3v_p1_surfx2_eval_quad_node_4_r(fc); 
  } else { 
    fUpwindQuad_r[4] = hyb_1x3v_p1_surfx2_eval_quad_node_4_l(fr); 
  } 
  if ((-0.3162277660168379*alpha_l[13])+0.3162277660168379*alpha_l[12]+0.3952847075210473*alpha_l[9]-0.3952847075210473*alpha_l[8]-0.4743416490252568*alpha_l[5]+0.4743416490252568*alpha_l[3]-0.3535533905932737*alpha_l[1]+0.3535533905932737*alpha_l[0] > 0) { 
    fUpwindQuad_l[5] = hyb_1x3v_p1_surfx2_eval_quad_node_5_r(fl); 
  } else { 
    fUpwindQuad_l[5] = hyb_1x3v_p1_surfx2_eval_quad_node_5_l(fc); 
  } 
  if ((-0.3162277660168379*alpha_r[13])+0.3162277660168379*alpha_r[12]+0.3952847075210473*alpha_r[9]-0.3952847075210473*alpha_r[8]-0.4743416490252568*alpha_r[5]+0.4743416490252568*alpha_r[3]-0.3535533905932737*alpha_r[1]+0.3535533905932737*alpha_r[0] > 0) { 
    fUpwindQuad_r[5] = hyb_1x3v_p1_surfx2_eval_quad_node_5_r(fc); 
  } else { 
    fUpwindQuad_r[5] = hyb_1x3v_p1_surfx2_eval_quad_node_5_l(fr); 
  } 
  if ((-0.3162277660168379*alpha_l[13])+0.3162277660168379*alpha_l[12]-0.3162277660168379*alpha_l[9]+0.3162277660168379*alpha_l[8]+0.6363961030678926*alpha_l[7]-0.6363961030678926*alpha_l[6]+0.4743416490252568*alpha_l[5]-0.4743416490252568*(alpha_l[4]+alpha_l[3])+0.4743416490252568*alpha_l[2]-0.3535533905932737*alpha_l[1]+0.3535533905932737*alpha_l[0] > 0) { 
    fUpwindQuad_l[6] = hyb_1x3v_p1_surfx2_eval_quad_node_6_r(fl); 
  } else { 
    fUpwindQuad_l[6] = hyb_1x3v_p1_surfx2_eval_quad_node_6_l(fc); 
  } 
  if ((-0.3162277660168379*alpha_r[13])+0.3162277660168379*alpha_r[12]-0.3162277660168379*alpha_r[9]+0.3162277660168379*alpha_r[8]+0.6363961030678926*alpha_r[7]-0.6363961030678926*alpha_r[6]+0.4743416490252568*alpha_r[5]-0.4743416490252568*(alpha_r[4]+alpha_r[3])+0.4743416490252568*alpha_r[2]-0.3535533905932737*alpha_r[1]+0.3535533905932737*alpha_r[0] > 0) { 
    fUpwindQuad_r[6] = hyb_1x3v_p1_surfx2_eval_quad_node_6_r(fc); 
  } else { 
    fUpwindQuad_r[6] = hyb_1x3v_p1_surfx2_eval_quad_node_6_l(fr); 
  } 
  if (0.3952847075210473*alpha_l[13]-0.3952847075210473*alpha_l[12]-0.3162277660168379*alpha_l[9]+0.3162277660168379*alpha_l[8]-0.4743416490252568*alpha_l[4]+0.4743416490252568*alpha_l[2]-0.3535533905932737*alpha_l[1]+0.3535533905932737*alpha_l[0] > 0) { 
    fUpwindQuad_l[7] = hyb_1x3v_p1_surfx2_eval_quad_node_7_r(fl); 
  } else { 
    fUpwindQuad_l[7] = hyb_1x3v_p1_surfx2_eval_quad_node_7_l(fc); 
  } 
  if (0.3952847075210473*alpha_r[13]-0.3952847075210473*alpha_r[12]-0.3162277660168379*alpha_r[9]+0.3162277660168379*alpha_r[8]-0.4743416490252568*alpha_r[4]+0.4743416490252568*alpha_r[2]-0.3535533905932737*alpha_r[1]+0.3535533905932737*alpha_r[0] > 0) { 
    fUpwindQuad_r[7] = hyb_1x3v_p1_surfx2_eval_quad_node_7_r(fc); 
  } else { 
    fUpwindQuad_r[7] = hyb_1x3v_p1_surfx2_eval_quad_node_7_l(fr); 
  } 
  if ((-0.3162277660168379*alpha_l[13])+0.3162277660168379*alpha_l[12]-0.3162277660168379*alpha_l[9]+0.3162277660168379*alpha_l[8]-0.6363961030678926*alpha_l[7]+0.6363961030678926*alpha_l[6]-0.4743416490252568*(alpha_l[5]+alpha_l[4])+0.4743416490252568*(alpha_l[3]+alpha_l[2])-0.3535533905932737*alpha_l[1]+0.3535533905932737*alpha_l[0] > 0) { 
    fUpwindQuad_l[8] = hyb_1x3v_p1_surfx2_eval_quad_node_8_r(fl); 
  } else { 
    fUpwindQuad_l[8] = hyb_1x3v_p1_surfx2_eval_quad_node_8_l(fc); 
  } 
  if ((-0.3162277660168379*alpha_r[13])+0.3162277660168379*alpha_r[12]-0.3162277660168379*alpha_r[9]+0.3162277660168379*alpha_r[8]-0.6363961030678926*alpha_r[7]+0.6363961030678926*alpha_r[6]-0.4743416490252568*(alpha_r[5]+alpha_r[4])+0.4743416490252568*(alpha_r[3]+alpha_r[2])-0.3535533905932737*alpha_r[1]+0.3535533905932737*alpha_r[0] > 0) { 
    fUpwindQuad_r[8] = hyb_1x3v_p1_surfx2_eval_quad_node_8_r(fc); 
  } else { 
    fUpwindQuad_r[8] = hyb_1x3v_p1_surfx2_eval_quad_node_8_l(fr); 
  } 
  if (0.3162277660168379*alpha_l[13]+0.3162277660168379*alpha_l[12]+0.3162277660168379*alpha_l[9]+0.3162277660168379*alpha_l[8]+0.6363961030678926*(alpha_l[7]+alpha_l[6])-0.4743416490252568*(alpha_l[5]+alpha_l[4]+alpha_l[3]+alpha_l[2])+0.3535533905932737*(alpha_l[1]+alpha_l[0]) > 0) { 
    fUpwindQuad_l[9] = hyb_1x3v_p1_surfx2_eval_quad_node_9_r(fl); 
  } else { 
    fUpwindQuad_l[9] = hyb_1x3v_p1_surfx2_eval_quad_node_9_l(fc); 
  } 
  if (0.3162277660168379*alpha_r[13]+0.3162277660168379*alpha_r[12]+0.3162277660168379*alpha_r[9]+0.3162277660168379*alpha_r[8]+0.6363961030678926*(alpha_r[7]+alpha_r[6])-0.4743416490252568*(alpha_r[5]+alpha_r[4]+alpha_r[3]+alpha_r[2])+0.3535533905932737*(alpha_r[1]+alpha_r[0]) > 0) { 
    fUpwindQuad_r[9] = hyb_1x3v_p1_surfx2_eval_quad_node_9_r(fc); 
  } else { 
    fUpwindQuad_r[9] = hyb_1x3v_p1_surfx2_eval_quad_node_9_l(fr); 
  } 
  if ((-0.3952847075210473*alpha_l[13])-0.3952847075210473*alpha_l[12]+0.3162277660168379*alpha_l[9]+0.3162277660168379*alpha_l[8]-0.4743416490252568*(alpha_l[4]+alpha_l[2])+0.3535533905932737*(alpha_l[1]+alpha_l[0]) > 0) { 
    fUpwindQuad_l[10] = hyb_1x3v_p1_surfx2_eval_quad_node_10_r(fl); 
  } else { 
    fUpwindQuad_l[10] = hyb_1x3v_p1_surfx2_eval_quad_node_10_l(fc); 
  } 
  if ((-0.3952847075210473*alpha_r[13])-0.3952847075210473*alpha_r[12]+0.3162277660168379*alpha_r[9]+0.3162277660168379*alpha_r[8]-0.4743416490252568*(alpha_r[4]+alpha_r[2])+0.3535533905932737*(alpha_r[1]+alpha_r[0]) > 0) { 
    fUpwindQuad_r[10] = hyb_1x3v_p1_surfx2_eval_quad_node_10_r(fc); 
  } else { 
    fUpwindQuad_r[10] = hyb_1x3v_p1_surfx2_eval_quad_node_10_l(fr); 
  } 
  if (0.3162277660168379*alpha_l[13]+0.3162277660168379*alpha_l[12]+0.3162277660168379*alpha_l[9]+0.3162277660168379*alpha_l[8]-0.6363961030678926*(alpha_l[7]+alpha_l[6])+0.4743416490252568*alpha_l[5]-0.4743416490252568*alpha_l[4]+0.4743416490252568*alpha_l[3]-0.4743416490252568*alpha_l[2]+0.3535533905932737*(alpha_l[1]+alpha_l[0]) > 0) { 
    fUpwindQuad_l[11] = hyb_1x3v_p1_surfx2_eval_quad_node_11_r(fl); 
  } else { 
    fUpwindQuad_l[11] = hyb_1x3v_p1_surfx2_eval_quad_node_11_l(fc); 
  } 
  if (0.3162277660168379*alpha_r[13]+0.3162277660168379*alpha_r[12]+0.3162277660168379*alpha_r[9]+0.3162277660168379*alpha_r[8]-0.6363961030678926*(alpha_r[7]+alpha_r[6])+0.4743416490252568*alpha_r[5]-0.4743416490252568*alpha_r[4]+0.4743416490252568*alpha_r[3]-0.4743416490252568*alpha_r[2]+0.3535533905932737*(alpha_r[1]+alpha_r[0]) > 0) { 
    fUpwindQuad_r[11] = hyb_1x3v_p1_surfx2_eval_quad_node_11_r(fc); 
  } else { 
    fUpwindQuad_r[11] = hyb_1x3v_p1_surfx2_eval_quad_node_11_l(fr); 
  } 
  if (0.3162277660168379*alpha_l[13]+0.3162277660168379*alpha_l[12]-0.3952847075210473*alpha_l[9]-0.3952847075210473*alpha_l[8]-0.4743416490252568*(alpha_l[5]+alpha_l[3])+0.3535533905932737*(alpha_l[1]+alpha_l[0]) > 0) { 
    fUpwindQuad_l[12] = hyb_1x3v_p1_surfx2_eval_quad_node_12_r(fl); 
  } else { 
    fUpwindQuad_l[12] = hyb_1x3v_p1_surfx2_eval_quad_node_12_l(fc); 
  } 
  if (0.3162277660168379*alpha_r[13]+0.3162277660168379*alpha_r[12]-0.3952847075210473*alpha_r[9]-0.3952847075210473*alpha_r[8]-0.4743416490252568*(alpha_r[5]+alpha_r[3])+0.3535533905932737*(alpha_r[1]+alpha_r[0]) > 0) { 
    fUpwindQuad_r[12] = hyb_1x3v_p1_surfx2_eval_quad_node_12_r(fc); 
  } else { 
    fUpwindQuad_r[12] = hyb_1x3v_p1_surfx2_eval_quad_node_12_l(fr); 
  } 
  if ((-0.3952847075210473*alpha_l[13])-0.3952847075210473*alpha_l[12]-0.3952847075210473*alpha_l[9]-0.3952847075210473*alpha_l[8]+0.3535533905932737*(alpha_l[1]+alpha_l[0]) > 0) { 
    fUpwindQuad_l[13] = hyb_1x3v_p1_surfx2_eval_quad_node_13_r(fl); 
  } else { 
    fUpwindQuad_l[13] = hyb_1x3v_p1_surfx2_eval_quad_node_13_l(fc); 
  } 
  if ((-0.3952847075210473*alpha_r[13])-0.3952847075210473*alpha_r[12]-0.3952847075210473*alpha_r[9]-0.3952847075210473*alpha_r[8]+0.3535533905932737*(alpha_r[1]+alpha_r[0]) > 0) { 
    fUpwindQuad_r[13] = hyb_1x3v_p1_surfx2_eval_quad_node_13_r(fc); 
  } else { 
    fUpwindQuad_r[13] = hyb_1x3v_p1_surfx2_eval_quad_node_13_l(fr); 
  } 
  if (0.3162277660168379*alpha_l[13]+0.3162277660168379*alpha_l[12]-0.3952847075210473*alpha_l[9]-0.3952847075210473*alpha_l[8]+0.4743416490252568*(alpha_l[5]+alpha_l[3])+0.3535533905932737*(alpha_l[1]+alpha_l[0]) > 0) { 
    fUpwindQuad_l[14] = hyb_1x3v_p1_surfx2_eval_quad_node_14_r(fl); 
  } else { 
    fUpwindQuad_l[14] = hyb_1x3v_p1_surfx2_eval_quad_node_14_l(fc); 
  } 
  if (0.3162277660168379*alpha_r[13]+0.3162277660168379*alpha_r[12]-0.3952847075210473*alpha_r[9]-0.3952847075210473*alpha_r[8]+0.4743416490252568*(alpha_r[5]+alpha_r[3])+0.3535533905932737*(alpha_r[1]+alpha_r[0]) > 0) { 
    fUpwindQuad_r[14] = hyb_1x3v_p1_surfx2_eval_quad_node_14_r(fc); 
  } else { 
    fUpwindQuad_r[14] = hyb_1x3v_p1_surfx2_eval_quad_node_14_l(fr); 
  } 
  if (0.3162277660168379*alpha_l[13]+0.3162277660168379*alpha_l[12]+0.3162277660168379*alpha_l[9]+0.3162277660168379*alpha_l[8]-0.6363961030678926*(alpha_l[7]+alpha_l[6])-0.4743416490252568*alpha_l[5]+0.4743416490252568*alpha_l[4]-0.4743416490252568*alpha_l[3]+0.4743416490252568*alpha_l[2]+0.3535533905932737*(alpha_l[1]+alpha_l[0]) > 0) { 
    fUpwindQuad_l[15] = hyb_1x3v_p1_surfx2_eval_quad_node_15_r(fl); 
  } else { 
    fUpwindQuad_l[15] = hyb_1x3v_p1_surfx2_eval_quad_node_15_l(fc); 
  } 
  if (0.3162277660168379*alpha_r[13]+0.3162277660168379*alpha_r[12]+0.3162277660168379*alpha_r[9]+0.3162277660168379*alpha_r[8]-0.6363961030678926*(alpha_r[7]+alpha_r[6])-0.4743416490252568*alpha_r[5]+0.4743416490252568*alpha_r[4]-0.4743416490252568*alpha_r[3]+0.4743416490252568*alpha_r[2]+0.3535533905932737*(alpha_r[1]+alpha_r[0]) > 0) { 
    fUpwindQuad_r[15] = hyb_1x3v_p1_surfx2_eval_quad_node_15_r(fc); 
  } else { 
    fUpwindQuad_r[15] = hyb_1x3v_p1_surfx2_eval_quad_node_15_l(fr); 
  } 
  if ((-0.3952847075210473*alpha_l[13])-0.3952847075210473*alpha_l[12]+0.3162277660168379*alpha_l[9]+0.3162277660168379*alpha_l[8]+0.4743416490252568*(alpha_l[4]+alpha_l[2])+0.3535533905932737*(alpha_l[1]+alpha_l[0]) > 0) { 
    fUpwindQuad_l[16] = hyb_1x3v_p1_surfx2_eval_quad_node_16_r(fl); 
  } else { 
    fUpwindQuad_l[16] = hyb_1x3v_p1_surfx2_eval_quad_node_16_l(fc); 
  } 
  if ((-0.3952847075210473*alpha_r[13])-0.3952847075210473*alpha_r[12]+0.3162277660168379*alpha_r[9]+0.3162277660168379*alpha_r[8]+0.4743416490252568*(alpha_r[4]+alpha_r[2])+0.3535533905932737*(alpha_r[1]+alpha_r[0]) > 0) { 
    fUpwindQuad_r[16] = hyb_1x3v_p1_surfx2_eval_quad_node_16_r(fc); 
  } else { 
    fUpwindQuad_r[16] = hyb_1x3v_p1_surfx2_eval_quad_node_16_l(fr); 
  } 
  if (0.3162277660168379*alpha_l[13]+0.3162277660168379*alpha_l[12]+0.3162277660168379*alpha_l[9]+0.3162277660168379*alpha_l[8]+0.6363961030678926*(alpha_l[7]+alpha_l[6])+0.4743416490252568*(alpha_l[5]+alpha_l[4]+alpha_l[3]+alpha_l[2])+0.3535533905932737*(alpha_l[1]+alpha_l[0]) > 0) { 
    fUpwindQuad_l[17] = hyb_1x3v_p1_surfx2_eval_quad_node_17_r(fl); 
  } else { 
    fUpwindQuad_l[17] = hyb_1x3v_p1_surfx2_eval_quad_node_17_l(fc); 
  } 
  if (0.3162277660168379*alpha_r[13]+0.3162277660168379*alpha_r[12]+0.3162277660168379*alpha_r[9]+0.3162277660168379*alpha_r[8]+0.6363961030678926*(alpha_r[7]+alpha_r[6])+0.4743416490252568*(alpha_r[5]+alpha_r[4]+alpha_r[3]+alpha_r[2])+0.3535533905932737*(alpha_r[1]+alpha_r[0]) > 0) { 
    fUpwindQuad_r[17] = hyb_1x3v_p1_surfx2_eval_quad_node_17_r(fc); 
  } else { 
    fUpwindQuad_r[17] = hyb_1x3v_p1_surfx2_eval_quad_node_17_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_1x3v_p1_vdir_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  hyb_1x3v_p1_vdir_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 
  Ghat_l[0] = 0.3535533905932737*alpha_l[13]*fUpwind_l[13]+0.3535533905932737*alpha_l[12]*fUpwind_l[12]+0.3535533905932737*alpha_l[9]*fUpwind_l[9]+0.3535533905932737*alpha_l[8]*fUpwind_l[8]+0.3535533905932737*alpha_l[7]*fUpwind_l[7]+0.3535533905932737*alpha_l[6]*fUpwind_l[6]+0.3535533905932737*alpha_l[5]*fUpwind_l[5]+0.3535533905932737*alpha_l[4]*fUpwind_l[4]+0.3535533905932737*alpha_l[3]*fUpwind_l[3]+0.3535533905932737*alpha_l[2]*fUpwind_l[2]+0.3535533905932737*alpha_l[1]*fUpwind_l[1]+0.3535533905932737*alpha_l[0]*fUpwind_l[0]; 
  Ghat_l[1] = 0.3535533905932737*alpha_l[12]*fUpwind_l[13]+0.3535533905932737*fUpwind_l[12]*alpha_l[13]+0.3535533905932737*alpha_l[8]*fUpwind_l[9]+0.3535533905932737*fUpwind_l[8]*alpha_l[9]+0.3535533905932737*alpha_l[6]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[6]*alpha_l[7]+0.3535533905932737*alpha_l[3]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[3]*alpha_l[5]+0.3535533905932737*alpha_l[2]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[2]*alpha_l[4]+0.3535533905932737*alpha_l[0]*fUpwind_l[1]+0.3535533905932737*fUpwind_l[0]*alpha_l[1]; 
  Ghat_l[2] = 0.3535533905932737*alpha_l[13]*fUpwind_l[15]+0.3535533905932737*alpha_l[12]*fUpwind_l[14]+0.3162277660168379*alpha_l[7]*fUpwind_l[11]+0.3162277660168379*alpha_l[6]*fUpwind_l[10]+0.3162277660168379*alpha_l[4]*fUpwind_l[9]+0.3162277660168379*fUpwind_l[4]*alpha_l[9]+0.3162277660168379*alpha_l[2]*fUpwind_l[8]+0.3162277660168379*fUpwind_l[2]*alpha_l[8]+0.3535533905932737*alpha_l[5]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[5]*alpha_l[7]+0.3535533905932737*alpha_l[3]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[3]*alpha_l[6]+0.3535533905932737*alpha_l[1]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[1]*alpha_l[4]+0.3535533905932737*alpha_l[0]*fUpwind_l[2]+0.3535533905932737*fUpwind_l[0]*alpha_l[2]; 
  Ghat_l[3] = 0.3162277660168379*alpha_l[7]*fUpwind_l[15]+0.3162277660168379*alpha_l[6]*fUpwind_l[14]+0.3162277660168379*alpha_l[5]*fUpwind_l[13]+0.3162277660168379*fUpwind_l[5]*alpha_l[13]+0.3162277660168379*alpha_l[3]*fUpwind_l[12]+0.3162277660168379*fUpwind_l[3]*alpha_l[12]+0.3535533905932737*alpha_l[9]*fUpwind_l[11]+0.3535533905932737*alpha_l[8]*fUpwind_l[10]+0.3535533905932737*alpha_l[4]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[4]*alpha_l[7]+0.3535533905932737*alpha_l[2]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[2]*alpha_l[6]+0.3535533905932737*alpha_l[1]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[1]*alpha_l[5]+0.3535533905932737*alpha_l[0]*fUpwind_l[3]+0.3535533905932737*fUpwind_l[0]*alpha_l[3]; 
  Ghat_l[4] = 0.3535533905932737*alpha_l[12]*fUpwind_l[15]+0.3535533905932737*alpha_l[13]*fUpwind_l[14]+0.3162277660168379*alpha_l[6]*fUpwind_l[11]+0.3162277660168379*alpha_l[7]*fUpwind_l[10]+0.3162277660168379*alpha_l[2]*fUpwind_l[9]+0.3162277660168379*fUpwind_l[2]*alpha_l[9]+0.3162277660168379*alpha_l[4]*fUpwind_l[8]+0.3162277660168379*fUpwind_l[4]*alpha_l[8]+0.3535533905932737*alpha_l[3]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[3]*alpha_l[7]+0.3535533905932737*alpha_l[5]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[5]*alpha_l[6]+0.3535533905932737*alpha_l[0]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[0]*alpha_l[4]+0.3535533905932737*alpha_l[1]*fUpwind_l[2]+0.3535533905932737*fUpwind_l[1]*alpha_l[2]; 
  Ghat_l[5] = 0.3162277660168379*alpha_l[6]*fUpwind_l[15]+0.3162277660168379*alpha_l[7]*fUpwind_l[14]+0.3162277660168379*alpha_l[3]*fUpwind_l[13]+0.3162277660168379*fUpwind_l[3]*alpha_l[13]+0.3162277660168379*alpha_l[5]*fUpwind_l[12]+0.3162277660168379*fUpwind_l[5]*alpha_l[12]+0.3535533905932737*alpha_l[8]*fUpwind_l[11]+0.3535533905932737*alpha_l[9]*fUpwind_l[10]+0.3535533905932737*alpha_l[2]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[2]*alpha_l[7]+0.3535533905932737*alpha_l[4]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[4]*alpha_l[6]+0.3535533905932737*alpha_l[0]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[0]*alpha_l[5]+0.3535533905932737*alpha_l[1]*fUpwind_l[3]+0.3535533905932737*fUpwind_l[1]*alpha_l[3]; 
  Ghat_l[6] = 0.3162277660168379*alpha_l[5]*fUpwind_l[15]+0.3162277660168379*alpha_l[3]*fUpwind_l[14]+0.3162277660168379*alpha_l[7]*fUpwind_l[13]+0.3162277660168379*fUpwind_l[7]*alpha_l[13]+0.3162277660168379*alpha_l[6]*fUpwind_l[12]+0.3162277660168379*fUpwind_l[6]*alpha_l[12]+0.3162277660168379*alpha_l[4]*fUpwind_l[11]+0.3162277660168379*alpha_l[2]*fUpwind_l[10]+0.3162277660168379*alpha_l[7]*fUpwind_l[9]+0.3162277660168379*fUpwind_l[7]*alpha_l[9]+0.3162277660168379*alpha_l[6]*fUpwind_l[8]+0.3162277660168379*fUpwind_l[6]*alpha_l[8]+0.3535533905932737*alpha_l[1]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[1]*alpha_l[7]+0.3535533905932737*alpha_l[0]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[0]*alpha_l[6]+0.3535533905932737*alpha_l[4]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[4]*alpha_l[5]+0.3535533905932737*alpha_l[2]*fUpwind_l[3]+0.3535533905932737*fUpwind_l[2]*alpha_l[3]; 
  Ghat_l[7] = 0.3162277660168379*alpha_l[3]*fUpwind_l[15]+0.3162277660168379*alpha_l[5]*fUpwind_l[14]+0.3162277660168379*alpha_l[6]*fUpwind_l[13]+0.3162277660168379*fUpwind_l[6]*alpha_l[13]+0.3162277660168379*alpha_l[7]*fUpwind_l[12]+0.3162277660168379*fUpwind_l[7]*alpha_l[12]+0.3162277660168379*alpha_l[2]*fUpwind_l[11]+0.3162277660168379*alpha_l[4]*fUpwind_l[10]+0.3162277660168379*alpha_l[6]*fUpwind_l[9]+0.3162277660168379*fUpwind_l[6]*alpha_l[9]+0.3162277660168379*alpha_l[7]*fUpwind_l[8]+0.3162277660168379*fUpwind_l[7]*alpha_l[8]+0.3535533905932737*alpha_l[0]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[0]*alpha_l[7]+0.3535533905932737*alpha_l[1]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[1]*alpha_l[6]+0.3535533905932737*alpha_l[2]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[2]*alpha_l[5]+0.3535533905932737*alpha_l[3]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[3]*alpha_l[4]; 
  Ghat_l[8] = 0.3535533905932737*alpha_l[5]*fUpwind_l[11]+0.3535533905932737*alpha_l[3]*fUpwind_l[10]+0.2258769757263128*alpha_l[9]*fUpwind_l[9]+0.3535533905932737*alpha_l[1]*fUpwind_l[9]+0.3535533905932737*fUpwind_l[1]*alpha_l[9]+0.2258769757263128*alpha_l[8]*fUpwind_l[8]+0.3535533905932737*alpha_l[0]*fUpwind_l[8]+0.3535533905932737*fUpwind_l[0]*alpha_l[8]+0.3162277660168379*alpha_l[7]*fUpwind_l[7]+0.3162277660168379*alpha_l[6]*fUpwind_l[6]+0.3162277660168379*alpha_l[4]*fUpwind_l[4]+0.3162277660168379*alpha_l[2]*fUpwind_l[2]; 
  Ghat_l[9] = 0.3535533905932737*alpha_l[3]*fUpwind_l[11]+0.3535533905932737*alpha_l[5]*fUpwind_l[10]+0.2258769757263128*alpha_l[8]*fUpwind_l[9]+0.3535533905932737*alpha_l[0]*fUpwind_l[9]+0.2258769757263128*fUpwind_l[8]*alpha_l[9]+0.3535533905932737*fUpwind_l[0]*alpha_l[9]+0.3535533905932737*alpha_l[1]*fUpwind_l[8]+0.3535533905932737*fUpwind_l[1]*alpha_l[8]+0.3162277660168379*alpha_l[6]*fUpwind_l[7]+0.3162277660168379*fUpwind_l[6]*alpha_l[7]+0.3162277660168379*alpha_l[2]*fUpwind_l[4]+0.3162277660168379*fUpwind_l[2]*alpha_l[4]; 
  Ghat_l[10] = 0.282842712474619*alpha_l[7]*fUpwind_l[15]+0.2828427124746191*alpha_l[6]*fUpwind_l[14]+0.3162277660168379*fUpwind_l[11]*alpha_l[13]+0.3162277660168379*fUpwind_l[10]*alpha_l[12]+0.2258769757263128*alpha_l[9]*fUpwind_l[11]+0.3535533905932737*alpha_l[1]*fUpwind_l[11]+0.2258769757263128*alpha_l[8]*fUpwind_l[10]+0.3535533905932737*alpha_l[0]*fUpwind_l[10]+0.3535533905932737*alpha_l[5]*fUpwind_l[9]+0.3535533905932737*fUpwind_l[5]*alpha_l[9]+0.3535533905932737*alpha_l[3]*fUpwind_l[8]+0.3535533905932737*fUpwind_l[3]*alpha_l[8]+0.3162277660168379*alpha_l[4]*fUpwind_l[7]+0.3162277660168379*fUpwind_l[4]*alpha_l[7]+0.3162277660168379*alpha_l[2]*fUpwind_l[6]+0.3162277660168379*fUpwind_l[2]*alpha_l[6]; 
  Ghat_l[11] = 0.2828427124746191*alpha_l[6]*fUpwind_l[15]+0.282842712474619*alpha_l[7]*fUpwind_l[14]+0.3162277660168379*fUpwind_l[10]*alpha_l[13]+0.3162277660168379*fUpwind_l[11]*alpha_l[12]+0.2258769757263128*alpha_l[8]*fUpwind_l[11]+0.3535533905932737*alpha_l[0]*fUpwind_l[11]+0.2258769757263128*alpha_l[9]*fUpwind_l[10]+0.3535533905932737*alpha_l[1]*fUpwind_l[10]+0.3535533905932737*alpha_l[3]*fUpwind_l[9]+0.3535533905932737*fUpwind_l[3]*alpha_l[9]+0.3535533905932737*alpha_l[5]*fUpwind_l[8]+0.3535533905932737*fUpwind_l[5]*alpha_l[8]+0.3162277660168379*alpha_l[2]*fUpwind_l[7]+0.3162277660168379*fUpwind_l[2]*alpha_l[7]+0.3162277660168379*alpha_l[4]*fUpwind_l[6]+0.3162277660168379*fUpwind_l[4]*alpha_l[6]; 
  Ghat_l[12] = 0.3535533905932737*alpha_l[4]*fUpwind_l[15]+0.3535533905932737*alpha_l[2]*fUpwind_l[14]+0.2258769757263128*alpha_l[13]*fUpwind_l[13]+0.3535533905932737*alpha_l[1]*fUpwind_l[13]+0.3535533905932737*fUpwind_l[1]*alpha_l[13]+0.2258769757263128*alpha_l[12]*fUpwind_l[12]+0.3535533905932737*alpha_l[0]*fUpwind_l[12]+0.3535533905932737*fUpwind_l[0]*alpha_l[12]+0.3162277660168379*alpha_l[7]*fUpwind_l[7]+0.3162277660168379*alpha_l[6]*fUpwind_l[6]+0.3162277660168379*alpha_l[5]*fUpwind_l[5]+0.3162277660168379*alpha_l[3]*fUpwind_l[3]; 
  Ghat_l[13] = 0.3535533905932737*alpha_l[2]*fUpwind_l[15]+0.3535533905932737*alpha_l[4]*fUpwind_l[14]+0.2258769757263128*alpha_l[12]*fUpwind_l[13]+0.3535533905932737*alpha_l[0]*fUpwind_l[13]+0.2258769757263128*fUpwind_l[12]*alpha_l[13]+0.3535533905932737*fUpwind_l[0]*alpha_l[13]+0.3535533905932737*alpha_l[1]*fUpwind_l[12]+0.3535533905932737*fUpwind_l[1]*alpha_l[12]+0.3162277660168379*alpha_l[6]*fUpwind_l[7]+0.3162277660168379*fUpwind_l[6]*alpha_l[7]+0.3162277660168379*alpha_l[3]*fUpwind_l[5]+0.3162277660168379*fUpwind_l[3]*alpha_l[5]; 
  Ghat_l[14] = 0.2258769757263128*alpha_l[13]*fUpwind_l[15]+0.3162277660168379*alpha_l[9]*fUpwind_l[15]+0.3535533905932737*alpha_l[1]*fUpwind_l[15]+0.2258769757263128*alpha_l[12]*fUpwind_l[14]+0.3162277660168379*alpha_l[8]*fUpwind_l[14]+0.3535533905932737*alpha_l[0]*fUpwind_l[14]+0.3535533905932737*alpha_l[4]*fUpwind_l[13]+0.3535533905932737*fUpwind_l[4]*alpha_l[13]+0.3535533905932737*alpha_l[2]*fUpwind_l[12]+0.3535533905932737*fUpwind_l[2]*alpha_l[12]+0.282842712474619*alpha_l[7]*fUpwind_l[11]+0.2828427124746191*alpha_l[6]*fUpwind_l[10]+0.3162277660168379*alpha_l[5]*fUpwind_l[7]+0.3162277660168379*fUpwind_l[5]*alpha_l[7]+0.3162277660168379*alpha_l[3]*fUpwind_l[6]+0.3162277660168379*fUpwind_l[3]*alpha_l[6]; 
  Ghat_l[15] = 0.2258769757263128*alpha_l[12]*fUpwind_l[15]+0.3162277660168379*alpha_l[8]*fUpwind_l[15]+0.3535533905932737*alpha_l[0]*fUpwind_l[15]+0.2258769757263128*alpha_l[13]*fUpwind_l[14]+0.3162277660168379*alpha_l[9]*fUpwind_l[14]+0.3535533905932737*alpha_l[1]*fUpwind_l[14]+0.3535533905932737*alpha_l[2]*fUpwind_l[13]+0.3535533905932737*fUpwind_l[2]*alpha_l[13]+0.3535533905932737*alpha_l[4]*fUpwind_l[12]+0.3535533905932737*fUpwind_l[4]*alpha_l[12]+0.2828427124746191*alpha_l[6]*fUpwind_l[11]+0.282842712474619*alpha_l[7]*fUpwind_l[10]+0.3162277660168379*alpha_l[3]*fUpwind_l[7]+0.3162277660168379*fUpwind_l[3]*alpha_l[7]+0.3162277660168379*alpha_l[5]*fUpwind_l[6]+0.3162277660168379*fUpwind_l[5]*alpha_l[6]; 

  Ghat_r[0] = 0.3535533905932737*alpha_r[13]*fUpwind_r[13]+0.3535533905932737*alpha_r[12]*fUpwind_r[12]+0.3535533905932737*alpha_r[9]*fUpwind_r[9]+0.3535533905932737*alpha_r[8]*fUpwind_r[8]+0.3535533905932737*alpha_r[7]*fUpwind_r[7]+0.3535533905932737*alpha_r[6]*fUpwind_r[6]+0.3535533905932737*alpha_r[5]*fUpwind_r[5]+0.3535533905932737*alpha_r[4]*fUpwind_r[4]+0.3535533905932737*alpha_r[3]*fUpwind_r[3]+0.3535533905932737*alpha_r[2]*fUpwind_r[2]+0.3535533905932737*alpha_r[1]*fUpwind_r[1]+0.3535533905932737*alpha_r[0]*fUpwind_r[0]; 
  Ghat_r[1] = 0.3535533905932737*alpha_r[12]*fUpwind_r[13]+0.3535533905932737*fUpwind_r[12]*alpha_r[13]+0.3535533905932737*alpha_r[8]*fUpwind_r[9]+0.3535533905932737*fUpwind_r[8]*alpha_r[9]+0.3535533905932737*alpha_r[6]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[6]*alpha_r[7]+0.3535533905932737*alpha_r[3]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[3]*alpha_r[5]+0.3535533905932737*alpha_r[2]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[2]*alpha_r[4]+0.3535533905932737*alpha_r[0]*fUpwind_r[1]+0.3535533905932737*fUpwind_r[0]*alpha_r[1]; 
  Ghat_r[2] = 0.3535533905932737*alpha_r[13]*fUpwind_r[15]+0.3535533905932737*alpha_r[12]*fUpwind_r[14]+0.3162277660168379*alpha_r[7]*fUpwind_r[11]+0.3162277660168379*alpha_r[6]*fUpwind_r[10]+0.3162277660168379*alpha_r[4]*fUpwind_r[9]+0.3162277660168379*fUpwind_r[4]*alpha_r[9]+0.3162277660168379*alpha_r[2]*fUpwind_r[8]+0.3162277660168379*fUpwind_r[2]*alpha_r[8]+0.3535533905932737*alpha_r[5]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[5]*alpha_r[7]+0.3535533905932737*alpha_r[3]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[3]*alpha_r[6]+0.3535533905932737*alpha_r[1]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[1]*alpha_r[4]+0.3535533905932737*alpha_r[0]*fUpwind_r[2]+0.3535533905932737*fUpwind_r[0]*alpha_r[2]; 
  Ghat_r[3] = 0.3162277660168379*alpha_r[7]*fUpwind_r[15]+0.3162277660168379*alpha_r[6]*fUpwind_r[14]+0.3162277660168379*alpha_r[5]*fUpwind_r[13]+0.3162277660168379*fUpwind_r[5]*alpha_r[13]+0.3162277660168379*alpha_r[3]*fUpwind_r[12]+0.3162277660168379*fUpwind_r[3]*alpha_r[12]+0.3535533905932737*alpha_r[9]*fUpwind_r[11]+0.3535533905932737*alpha_r[8]*fUpwind_r[10]+0.3535533905932737*alpha_r[4]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[4]*alpha_r[7]+0.3535533905932737*alpha_r[2]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[2]*alpha_r[6]+0.3535533905932737*alpha_r[1]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[1]*alpha_r[5]+0.3535533905932737*alpha_r[0]*fUpwind_r[3]+0.3535533905932737*fUpwind_r[0]*alpha_r[3]; 
  Ghat_r[4] = 0.3535533905932737*alpha_r[12]*fUpwind_r[15]+0.3535533905932737*alpha_r[13]*fUpwind_r[14]+0.3162277660168379*alpha_r[6]*fUpwind_r[11]+0.3162277660168379*alpha_r[7]*fUpwind_r[10]+0.3162277660168379*alpha_r[2]*fUpwind_r[9]+0.3162277660168379*fUpwind_r[2]*alpha_r[9]+0.3162277660168379*alpha_r[4]*fUpwind_r[8]+0.3162277660168379*fUpwind_r[4]*alpha_r[8]+0.3535533905932737*alpha_r[3]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[3]*alpha_r[7]+0.3535533905932737*alpha_r[5]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[5]*alpha_r[6]+0.3535533905932737*alpha_r[0]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[0]*alpha_r[4]+0.3535533905932737*alpha_r[1]*fUpwind_r[2]+0.3535533905932737*fUpwind_r[1]*alpha_r[2]; 
  Ghat_r[5] = 0.3162277660168379*alpha_r[6]*fUpwind_r[15]+0.3162277660168379*alpha_r[7]*fUpwind_r[14]+0.3162277660168379*alpha_r[3]*fUpwind_r[13]+0.3162277660168379*fUpwind_r[3]*alpha_r[13]+0.3162277660168379*alpha_r[5]*fUpwind_r[12]+0.3162277660168379*fUpwind_r[5]*alpha_r[12]+0.3535533905932737*alpha_r[8]*fUpwind_r[11]+0.3535533905932737*alpha_r[9]*fUpwind_r[10]+0.3535533905932737*alpha_r[2]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[2]*alpha_r[7]+0.3535533905932737*alpha_r[4]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[4]*alpha_r[6]+0.3535533905932737*alpha_r[0]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[0]*alpha_r[5]+0.3535533905932737*alpha_r[1]*fUpwind_r[3]+0.3535533905932737*fUpwind_r[1]*alpha_r[3]; 
  Ghat_r[6] = 0.3162277660168379*alpha_r[5]*fUpwind_r[15]+0.3162277660168379*alpha_r[3]*fUpwind_r[14]+0.3162277660168379*alpha_r[7]*fUpwind_r[13]+0.3162277660168379*fUpwind_r[7]*alpha_r[13]+0.3162277660168379*alpha_r[6]*fUpwind_r[12]+0.3162277660168379*fUpwind_r[6]*alpha_r[12]+0.3162277660168379*alpha_r[4]*fUpwind_r[11]+0.3162277660168379*alpha_r[2]*fUpwind_r[10]+0.3162277660168379*alpha_r[7]*fUpwind_r[9]+0.3162277660168379*fUpwind_r[7]*alpha_r[9]+0.3162277660168379*alpha_r[6]*fUpwind_r[8]+0.3162277660168379*fUpwind_r[6]*alpha_r[8]+0.3535533905932737*alpha_r[1]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[1]*alpha_r[7]+0.3535533905932737*alpha_r[0]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[0]*alpha_r[6]+0.3535533905932737*alpha_r[4]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[4]*alpha_r[5]+0.3535533905932737*alpha_r[2]*fUpwind_r[3]+0.3535533905932737*fUpwind_r[2]*alpha_r[3]; 
  Ghat_r[7] = 0.3162277660168379*alpha_r[3]*fUpwind_r[15]+0.3162277660168379*alpha_r[5]*fUpwind_r[14]+0.3162277660168379*alpha_r[6]*fUpwind_r[13]+0.3162277660168379*fUpwind_r[6]*alpha_r[13]+0.3162277660168379*alpha_r[7]*fUpwind_r[12]+0.3162277660168379*fUpwind_r[7]*alpha_r[12]+0.3162277660168379*alpha_r[2]*fUpwind_r[11]+0.3162277660168379*alpha_r[4]*fUpwind_r[10]+0.3162277660168379*alpha_r[6]*fUpwind_r[9]+0.3162277660168379*fUpwind_r[6]*alpha_r[9]+0.3162277660168379*alpha_r[7]*fUpwind_r[8]+0.3162277660168379*fUpwind_r[7]*alpha_r[8]+0.3535533905932737*alpha_r[0]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[0]*alpha_r[7]+0.3535533905932737*alpha_r[1]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[1]*alpha_r[6]+0.3535533905932737*alpha_r[2]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[2]*alpha_r[5]+0.3535533905932737*alpha_r[3]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[3]*alpha_r[4]; 
  Ghat_r[8] = 0.3535533905932737*alpha_r[5]*fUpwind_r[11]+0.3535533905932737*alpha_r[3]*fUpwind_r[10]+0.2258769757263128*alpha_r[9]*fUpwind_r[9]+0.3535533905932737*alpha_r[1]*fUpwind_r[9]+0.3535533905932737*fUpwind_r[1]*alpha_r[9]+0.2258769757263128*alpha_r[8]*fUpwind_r[8]+0.3535533905932737*alpha_r[0]*fUpwind_r[8]+0.3535533905932737*fUpwind_r[0]*alpha_r[8]+0.3162277660168379*alpha_r[7]*fUpwind_r[7]+0.3162277660168379*alpha_r[6]*fUpwind_r[6]+0.3162277660168379*alpha_r[4]*fUpwind_r[4]+0.3162277660168379*alpha_r[2]*fUpwind_r[2]; 
  Ghat_r[9] = 0.3535533905932737*alpha_r[3]*fUpwind_r[11]+0.3535533905932737*alpha_r[5]*fUpwind_r[10]+0.2258769757263128*alpha_r[8]*fUpwind_r[9]+0.3535533905932737*alpha_r[0]*fUpwind_r[9]+0.2258769757263128*fUpwind_r[8]*alpha_r[9]+0.3535533905932737*fUpwind_r[0]*alpha_r[9]+0.3535533905932737*alpha_r[1]*fUpwind_r[8]+0.3535533905932737*fUpwind_r[1]*alpha_r[8]+0.3162277660168379*alpha_r[6]*fUpwind_r[7]+0.3162277660168379*fUpwind_r[6]*alpha_r[7]+0.3162277660168379*alpha_r[2]*fUpwind_r[4]+0.3162277660168379*fUpwind_r[2]*alpha_r[4]; 
  Ghat_r[10] = 0.282842712474619*alpha_r[7]*fUpwind_r[15]+0.2828427124746191*alpha_r[6]*fUpwind_r[14]+0.3162277660168379*fUpwind_r[11]*alpha_r[13]+0.3162277660168379*fUpwind_r[10]*alpha_r[12]+0.2258769757263128*alpha_r[9]*fUpwind_r[11]+0.3535533905932737*alpha_r[1]*fUpwind_r[11]+0.2258769757263128*alpha_r[8]*fUpwind_r[10]+0.3535533905932737*alpha_r[0]*fUpwind_r[10]+0.3535533905932737*alpha_r[5]*fUpwind_r[9]+0.3535533905932737*fUpwind_r[5]*alpha_r[9]+0.3535533905932737*alpha_r[3]*fUpwind_r[8]+0.3535533905932737*fUpwind_r[3]*alpha_r[8]+0.3162277660168379*alpha_r[4]*fUpwind_r[7]+0.3162277660168379*fUpwind_r[4]*alpha_r[7]+0.3162277660168379*alpha_r[2]*fUpwind_r[6]+0.3162277660168379*fUpwind_r[2]*alpha_r[6]; 
  Ghat_r[11] = 0.2828427124746191*alpha_r[6]*fUpwind_r[15]+0.282842712474619*alpha_r[7]*fUpwind_r[14]+0.3162277660168379*fUpwind_r[10]*alpha_r[13]+0.3162277660168379*fUpwind_r[11]*alpha_r[12]+0.2258769757263128*alpha_r[8]*fUpwind_r[11]+0.3535533905932737*alpha_r[0]*fUpwind_r[11]+0.2258769757263128*alpha_r[9]*fUpwind_r[10]+0.3535533905932737*alpha_r[1]*fUpwind_r[10]+0.3535533905932737*alpha_r[3]*fUpwind_r[9]+0.3535533905932737*fUpwind_r[3]*alpha_r[9]+0.3535533905932737*alpha_r[5]*fUpwind_r[8]+0.3535533905932737*fUpwind_r[5]*alpha_r[8]+0.3162277660168379*alpha_r[2]*fUpwind_r[7]+0.3162277660168379*fUpwind_r[2]*alpha_r[7]+0.3162277660168379*alpha_r[4]*fUpwind_r[6]+0.3162277660168379*fUpwind_r[4]*alpha_r[6]; 
  Ghat_r[12] = 0.3535533905932737*alpha_r[4]*fUpwind_r[15]+0.3535533905932737*alpha_r[2]*fUpwind_r[14]+0.2258769757263128*alpha_r[13]*fUpwind_r[13]+0.3535533905932737*alpha_r[1]*fUpwind_r[13]+0.3535533905932737*fUpwind_r[1]*alpha_r[13]+0.2258769757263128*alpha_r[12]*fUpwind_r[12]+0.3535533905932737*alpha_r[0]*fUpwind_r[12]+0.3535533905932737*fUpwind_r[0]*alpha_r[12]+0.3162277660168379*alpha_r[7]*fUpwind_r[7]+0.3162277660168379*alpha_r[6]*fUpwind_r[6]+0.3162277660168379*alpha_r[5]*fUpwind_r[5]+0.3162277660168379*alpha_r[3]*fUpwind_r[3]; 
  Ghat_r[13] = 0.3535533905932737*alpha_r[2]*fUpwind_r[15]+0.3535533905932737*alpha_r[4]*fUpwind_r[14]+0.2258769757263128*alpha_r[12]*fUpwind_r[13]+0.3535533905932737*alpha_r[0]*fUpwind_r[13]+0.2258769757263128*fUpwind_r[12]*alpha_r[13]+0.3535533905932737*fUpwind_r[0]*alpha_r[13]+0.3535533905932737*alpha_r[1]*fUpwind_r[12]+0.3535533905932737*fUpwind_r[1]*alpha_r[12]+0.3162277660168379*alpha_r[6]*fUpwind_r[7]+0.3162277660168379*fUpwind_r[6]*alpha_r[7]+0.3162277660168379*alpha_r[3]*fUpwind_r[5]+0.3162277660168379*fUpwind_r[3]*alpha_r[5]; 
  Ghat_r[14] = 0.2258769757263128*alpha_r[13]*fUpwind_r[15]+0.3162277660168379*alpha_r[9]*fUpwind_r[15]+0.3535533905932737*alpha_r[1]*fUpwind_r[15]+0.2258769757263128*alpha_r[12]*fUpwind_r[14]+0.3162277660168379*alpha_r[8]*fUpwind_r[14]+0.3535533905932737*alpha_r[0]*fUpwind_r[14]+0.3535533905932737*alpha_r[4]*fUpwind_r[13]+0.3535533905932737*fUpwind_r[4]*alpha_r[13]+0.3535533905932737*alpha_r[2]*fUpwind_r[12]+0.3535533905932737*fUpwind_r[2]*alpha_r[12]+0.282842712474619*alpha_r[7]*fUpwind_r[11]+0.2828427124746191*alpha_r[6]*fUpwind_r[10]+0.3162277660168379*alpha_r[5]*fUpwind_r[7]+0.3162277660168379*fUpwind_r[5]*alpha_r[7]+0.3162277660168379*alpha_r[3]*fUpwind_r[6]+0.3162277660168379*fUpwind_r[3]*alpha_r[6]; 
  Ghat_r[15] = 0.2258769757263128*alpha_r[12]*fUpwind_r[15]+0.3162277660168379*alpha_r[8]*fUpwind_r[15]+0.3535533905932737*alpha_r[0]*fUpwind_r[15]+0.2258769757263128*alpha_r[13]*fUpwind_r[14]+0.3162277660168379*alpha_r[9]*fUpwind_r[14]+0.3535533905932737*alpha_r[1]*fUpwind_r[14]+0.3535533905932737*alpha_r[2]*fUpwind_r[13]+0.3535533905932737*fUpwind_r[2]*alpha_r[13]+0.3535533905932737*alpha_r[4]*fUpwind_r[12]+0.3535533905932737*fUpwind_r[4]*alpha_r[12]+0.2828427124746191*alpha_r[6]*fUpwind_r[11]+0.282842712474619*alpha_r[7]*fUpwind_r[10]+0.3162277660168379*alpha_r[3]*fUpwind_r[7]+0.3162277660168379*fUpwind_r[3]*alpha_r[7]+0.3162277660168379*alpha_r[5]*fUpwind_r[6]+0.3162277660168379*fUpwind_r[5]*alpha_r[6]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv10; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv10; 
  out[2] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv10; 
  out[3] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv10; 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv10; 
  out[5] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv10; 
  out[6] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dv10; 
  out[7] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv10; 
  out[8] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dv10; 
  out[9] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv10; 
  out[10] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dv10; 
  out[11] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dv10; 
  out[12] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dv10; 
  out[13] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dv10; 
  out[14] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dv10; 
  out[15] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dv10; 
  out[16] += (1.58113883008419*Ghat_l[0]-1.58113883008419*Ghat_r[0])*dv10; 
  out[17] += (1.58113883008419*Ghat_l[1]-1.58113883008419*Ghat_r[1])*dv10; 
  out[18] += (1.58113883008419*Ghat_l[2]-1.58113883008419*Ghat_r[2])*dv10; 
  out[19] += (1.58113883008419*Ghat_l[3]-1.58113883008419*Ghat_r[3])*dv10; 
  out[20] += (1.58113883008419*Ghat_l[4]-1.58113883008419*Ghat_r[4])*dv10; 
  out[21] += (1.58113883008419*Ghat_l[5]-1.58113883008419*Ghat_r[5])*dv10; 
  out[22] += (1.58113883008419*Ghat_l[6]-1.58113883008419*Ghat_r[6])*dv10; 
  out[23] += (1.58113883008419*Ghat_l[7]-1.58113883008419*Ghat_r[7])*dv10; 
  out[24] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*dv10; 
  out[25] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*dv10; 
  out[26] += -1.224744871391589*(Ghat_r[8]+Ghat_l[8])*dv10; 
  out[27] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*dv10; 
  out[28] += -1.224744871391589*(Ghat_r[9]+Ghat_l[9])*dv10; 
  out[29] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*dv10; 
  out[30] += -1.224744871391589*(Ghat_r[10]+Ghat_l[10])*dv10; 
  out[31] += -1.224744871391589*(Ghat_r[11]+Ghat_l[11])*dv10; 
  out[32] += (0.7071067811865475*Ghat_l[12]-0.7071067811865475*Ghat_r[12])*dv10; 
  out[33] += (0.7071067811865475*Ghat_l[13]-0.7071067811865475*Ghat_r[13])*dv10; 
  out[34] += -1.224744871391589*(Ghat_r[12]+Ghat_l[12])*dv10; 
  out[35] += (0.7071067811865475*Ghat_l[14]-0.7071067811865475*Ghat_r[14])*dv10; 
  out[36] += -1.224744871391589*(Ghat_r[13]+Ghat_l[13])*dv10; 
  out[37] += (0.7071067811865475*Ghat_l[15]-0.7071067811865475*Ghat_r[15])*dv10; 
  out[38] += -1.224744871391589*(Ghat_r[14]+Ghat_l[14])*dv10; 
  out[39] += -1.224744871391589*(Ghat_r[15]+Ghat_l[15])*dv10; 

  return 0.;

} 
