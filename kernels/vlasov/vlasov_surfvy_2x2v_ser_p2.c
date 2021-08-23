#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_2x2v_p2_surfvy_quad.h> 
GKYL_CU_DH void vlasov_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // qmem:      q/m*EM fields.
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 
  const double dv11 = 2/dxv[3]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double *E1 = &qmem[8]; 
  const double *B2 = &qmem[40]; 

  double alpha[20] = {0.0}; 

  alpha[0] = 1.414213562373095*E1[0]-1.414213562373095*B2[0]*wv1; 
  alpha[1] = 1.414213562373095*E1[1]-1.414213562373095*B2[1]*wv1; 
  alpha[2] = 1.414213562373095*E1[2]-1.414213562373095*B2[2]*wv1; 
  alpha[3] = -0.408248290463863*B2[0]*dv1; 
  alpha[4] = 1.414213562373095*E1[3]-1.414213562373095*B2[3]*wv1; 
  alpha[5] = -0.408248290463863*B2[1]*dv1; 
  alpha[6] = -0.408248290463863*B2[2]*dv1; 
  alpha[7] = 1.414213562373095*E1[4]-1.414213562373095*B2[4]*wv1; 
  alpha[8] = 1.414213562373095*E1[5]-1.414213562373095*B2[5]*wv1; 
  alpha[10] = -0.408248290463863*B2[3]*dv1; 
  alpha[11] = 1.414213562373095*E1[6]-1.414213562373095*B2[6]*wv1; 
  alpha[12] = 1.414213562373095*E1[7]-1.414213562373095*B2[7]*wv1; 
  alpha[13] = -0.408248290463863*B2[4]*dv1; 
  alpha[14] = -0.408248290463863*B2[5]*dv1; 
  alpha[17] = -0.408248290463863*B2[6]*dv1; 
  alpha[18] = -0.408248290463863*B2[7]*dv1; 

  double fUpwindQuad_l[27] = {0.0};
  double fUpwindQuad_r[27] = {0.0};
  double fUpwind_l[20] = {0.0};;
  double fUpwind_r[20] = {0.0};
  double Ghat_l[20] = {0.0}; 
  double Ghat_r[20] = {0.0}; 

  if (0.5692099788303082*(alpha[18]+alpha[17])-0.4242640687119285*(alpha[14]+alpha[13])-0.4242640687119281*alpha[12]-0.4242640687119285*alpha[11]-0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*(alpha[6]+alpha[5]+alpha[4])-0.4743416490252568*(alpha[3]+alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[0] = ser_2x2v_p2_surfvy_quad_0(1, fl); 
    fUpwindQuad_r[0] = ser_2x2v_p2_surfvy_quad_0(1, fc); 
  } else { 

    fUpwindQuad_l[0] = ser_2x2v_p2_surfvy_quad_0(-1, fc); 
    fUpwindQuad_r[0] = ser_2x2v_p2_surfvy_quad_0(-1, fr); 
  } 
  if ((-0.711512473537885*alpha[17])-0.4242640687119285*alpha[14]+0.5303300858899104*(alpha[13]+alpha[11])+0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]+0.6363961030678926*alpha[6]-0.4743416490252568*(alpha[3]+alpha[2])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[1] = ser_2x2v_p2_surfvy_quad_1(1, fl); 
    fUpwindQuad_r[1] = ser_2x2v_p2_surfvy_quad_1(1, fc); 
  } else { 

    fUpwindQuad_l[1] = ser_2x2v_p2_surfvy_quad_1(-1, fc); 
    fUpwindQuad_r[1] = ser_2x2v_p2_surfvy_quad_1(-1, fr); 
  } 
  if ((-0.5692099788303082*alpha[18])+0.5692099788303082*alpha[17]-0.4242640687119285*(alpha[14]+alpha[13])+0.4242640687119281*alpha[12]-0.4242640687119285*alpha[11]+0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*alpha[6]-0.6363961030678926*(alpha[5]+alpha[4])-0.4743416490252568*(alpha[3]+alpha[2])+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[2] = ser_2x2v_p2_surfvy_quad_2(1, fl); 
    fUpwindQuad_r[2] = ser_2x2v_p2_surfvy_quad_2(1, fc); 
  } else { 

    fUpwindQuad_l[2] = ser_2x2v_p2_surfvy_quad_2(-1, fc); 
    fUpwindQuad_r[2] = ser_2x2v_p2_surfvy_quad_2(-1, fr); 
  } 
  if ((-0.711512473537885*alpha[18])+0.5303300858899104*alpha[14]-0.4242640687119285*alpha[13]+0.5303300858899104*alpha[12]-0.3952847075210473*alpha[8]+0.3162277660168379*alpha[7]+0.6363961030678926*alpha[5]-0.4743416490252568*(alpha[3]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[3] = ser_2x2v_p2_surfvy_quad_3(1, fl); 
    fUpwindQuad_r[3] = ser_2x2v_p2_surfvy_quad_3(1, fc); 
  } else { 

    fUpwindQuad_l[3] = ser_2x2v_p2_surfvy_quad_3(-1, fc); 
    fUpwindQuad_r[3] = ser_2x2v_p2_surfvy_quad_3(-1, fr); 
  } 
  if (0.5303300858899104*(alpha[14]+alpha[13])-0.3952847075210473*(alpha[8]+alpha[7])-0.4743416490252568*alpha[3]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[4] = ser_2x2v_p2_surfvy_quad_4(1, fl); 
    fUpwindQuad_r[4] = ser_2x2v_p2_surfvy_quad_4(1, fc); 
  } else { 

    fUpwindQuad_l[4] = ser_2x2v_p2_surfvy_quad_4(-1, fc); 
    fUpwindQuad_r[4] = ser_2x2v_p2_surfvy_quad_4(-1, fr); 
  } 
  if (0.711512473537885*alpha[18]+0.5303300858899104*alpha[14]-0.4242640687119285*alpha[13]-0.5303300858899104*alpha[12]-0.3952847075210473*alpha[8]+0.3162277660168379*alpha[7]-0.6363961030678926*alpha[5]-0.4743416490252568*alpha[3]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[5] = ser_2x2v_p2_surfvy_quad_5(1, fl); 
    fUpwindQuad_r[5] = ser_2x2v_p2_surfvy_quad_5(1, fc); 
  } else { 

    fUpwindQuad_l[5] = ser_2x2v_p2_surfvy_quad_5(-1, fc); 
    fUpwindQuad_r[5] = ser_2x2v_p2_surfvy_quad_5(-1, fr); 
  } 
  if (0.5692099788303082*alpha[18]-0.5692099788303082*alpha[17]-0.4242640687119285*(alpha[14]+alpha[13])-0.4242640687119281*alpha[12]+0.4242640687119285*alpha[11]+0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*alpha[6]+0.6363961030678926*alpha[5]-0.6363961030678926*alpha[4]-0.4743416490252568*alpha[3]+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[6] = ser_2x2v_p2_surfvy_quad_6(1, fl); 
    fUpwindQuad_r[6] = ser_2x2v_p2_surfvy_quad_6(1, fc); 
  } else { 

    fUpwindQuad_l[6] = ser_2x2v_p2_surfvy_quad_6(-1, fc); 
    fUpwindQuad_r[6] = ser_2x2v_p2_surfvy_quad_6(-1, fr); 
  } 
  if (0.711512473537885*alpha[17]-0.4242640687119285*alpha[14]+0.5303300858899104*alpha[13]-0.5303300858899104*alpha[11]+0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]-0.6363961030678926*alpha[6]-0.4743416490252568*alpha[3]+0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[7] = ser_2x2v_p2_surfvy_quad_7(1, fl); 
    fUpwindQuad_r[7] = ser_2x2v_p2_surfvy_quad_7(1, fc); 
  } else { 

    fUpwindQuad_l[7] = ser_2x2v_p2_surfvy_quad_7(-1, fc); 
    fUpwindQuad_r[7] = ser_2x2v_p2_surfvy_quad_7(-1, fr); 
  } 
  if ((-0.5692099788303082*(alpha[18]+alpha[17]))-0.4242640687119285*(alpha[14]+alpha[13])+0.4242640687119281*alpha[12]+0.4242640687119285*alpha[11]-0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*(alpha[6]+alpha[5])+0.6363961030678926*alpha[4]-0.4743416490252568*alpha[3]+0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[8] = ser_2x2v_p2_surfvy_quad_8(1, fl); 
    fUpwindQuad_r[8] = ser_2x2v_p2_surfvy_quad_8(1, fc); 
  } else { 

    fUpwindQuad_l[8] = ser_2x2v_p2_surfvy_quad_8(-1, fc); 
    fUpwindQuad_r[8] = ser_2x2v_p2_surfvy_quad_8(-1, fr); 
  } 
  if ((-0.4242640687119281*alpha[12])-0.4242640687119285*alpha[11]+0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*alpha[4]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[9] = ser_2x2v_p2_surfvy_quad_9(1, fl); 
    fUpwindQuad_r[9] = ser_2x2v_p2_surfvy_quad_9(1, fc); 
  } else { 

    fUpwindQuad_l[9] = ser_2x2v_p2_surfvy_quad_9(-1, fc); 
    fUpwindQuad_r[9] = ser_2x2v_p2_surfvy_quad_9(-1, fr); 
  } 
  if (0.5303300858899104*alpha[11]+0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]-0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[10] = ser_2x2v_p2_surfvy_quad_10(1, fl); 
    fUpwindQuad_r[10] = ser_2x2v_p2_surfvy_quad_10(1, fc); 
  } else { 

    fUpwindQuad_l[10] = ser_2x2v_p2_surfvy_quad_10(-1, fc); 
    fUpwindQuad_r[10] = ser_2x2v_p2_surfvy_quad_10(-1, fr); 
  } 
  if (0.4242640687119281*alpha[12]-0.4242640687119285*alpha[11]+0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*alpha[4]-0.4743416490252568*alpha[2]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[11] = ser_2x2v_p2_surfvy_quad_11(1, fl); 
    fUpwindQuad_r[11] = ser_2x2v_p2_surfvy_quad_11(1, fc); 
  } else { 

    fUpwindQuad_l[11] = ser_2x2v_p2_surfvy_quad_11(-1, fc); 
    fUpwindQuad_r[11] = ser_2x2v_p2_surfvy_quad_11(-1, fr); 
  } 
  if (0.5303300858899104*alpha[12]-0.3952847075210473*alpha[8]+0.3162277660168379*alpha[7]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[12] = ser_2x2v_p2_surfvy_quad_12(1, fl); 
    fUpwindQuad_r[12] = ser_2x2v_p2_surfvy_quad_12(1, fc); 
  } else { 

    fUpwindQuad_l[12] = ser_2x2v_p2_surfvy_quad_12(-1, fc); 
    fUpwindQuad_r[12] = ser_2x2v_p2_surfvy_quad_12(-1, fr); 
  } 
  if (0.3535533905932737*alpha[0]-0.3952847075210473*(alpha[8]+alpha[7]) > 0) { 

    fUpwindQuad_l[13] = ser_2x2v_p2_surfvy_quad_13(1, fl); 
    fUpwindQuad_r[13] = ser_2x2v_p2_surfvy_quad_13(1, fc); 
  } else { 

    fUpwindQuad_l[13] = ser_2x2v_p2_surfvy_quad_13(-1, fc); 
    fUpwindQuad_r[13] = ser_2x2v_p2_surfvy_quad_13(-1, fr); 
  } 
  if ((-0.5303300858899104*alpha[12])-0.3952847075210473*alpha[8]+0.3162277660168379*alpha[7]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[14] = ser_2x2v_p2_surfvy_quad_14(1, fl); 
    fUpwindQuad_r[14] = ser_2x2v_p2_surfvy_quad_14(1, fc); 
  } else { 

    fUpwindQuad_l[14] = ser_2x2v_p2_surfvy_quad_14(-1, fc); 
    fUpwindQuad_r[14] = ser_2x2v_p2_surfvy_quad_14(-1, fr); 
  } 
  if ((-0.4242640687119281*alpha[12])+0.4242640687119285*alpha[11]+0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*alpha[4]+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[15] = ser_2x2v_p2_surfvy_quad_15(1, fl); 
    fUpwindQuad_r[15] = ser_2x2v_p2_surfvy_quad_15(1, fc); 
  } else { 

    fUpwindQuad_l[15] = ser_2x2v_p2_surfvy_quad_15(-1, fc); 
    fUpwindQuad_r[15] = ser_2x2v_p2_surfvy_quad_15(-1, fr); 
  } 
  if ((-0.5303300858899104*alpha[11])+0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]+0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[16] = ser_2x2v_p2_surfvy_quad_16(1, fl); 
    fUpwindQuad_r[16] = ser_2x2v_p2_surfvy_quad_16(1, fc); 
  } else { 

    fUpwindQuad_l[16] = ser_2x2v_p2_surfvy_quad_16(-1, fc); 
    fUpwindQuad_r[16] = ser_2x2v_p2_surfvy_quad_16(-1, fr); 
  } 
  if (0.4242640687119281*alpha[12]+0.4242640687119285*alpha[11]+0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*alpha[4]+0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[17] = ser_2x2v_p2_surfvy_quad_17(1, fl); 
    fUpwindQuad_r[17] = ser_2x2v_p2_surfvy_quad_17(1, fc); 
  } else { 

    fUpwindQuad_l[17] = ser_2x2v_p2_surfvy_quad_17(-1, fc); 
    fUpwindQuad_r[17] = ser_2x2v_p2_surfvy_quad_17(-1, fr); 
  } 
  if ((-0.5692099788303082*(alpha[18]+alpha[17]))+0.4242640687119285*(alpha[14]+alpha[13])-0.4242640687119281*alpha[12]-0.4242640687119285*alpha[11]+0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*(alpha[6]+alpha[5])+0.6363961030678926*alpha[4]+0.4743416490252568*alpha[3]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[18] = ser_2x2v_p2_surfvy_quad_18(1, fl); 
    fUpwindQuad_r[18] = ser_2x2v_p2_surfvy_quad_18(1, fc); 
  } else { 

    fUpwindQuad_l[18] = ser_2x2v_p2_surfvy_quad_18(-1, fc); 
    fUpwindQuad_r[18] = ser_2x2v_p2_surfvy_quad_18(-1, fr); 
  } 
  if (0.711512473537885*alpha[17]+0.4242640687119285*alpha[14]-0.5303300858899104*alpha[13]+0.5303300858899104*alpha[11]+0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]-0.6363961030678926*alpha[6]+0.4743416490252568*alpha[3]-0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[19] = ser_2x2v_p2_surfvy_quad_19(1, fl); 
    fUpwindQuad_r[19] = ser_2x2v_p2_surfvy_quad_19(1, fc); 
  } else { 

    fUpwindQuad_l[19] = ser_2x2v_p2_surfvy_quad_19(-1, fc); 
    fUpwindQuad_r[19] = ser_2x2v_p2_surfvy_quad_19(-1, fr); 
  } 
  if (0.5692099788303082*alpha[18]-0.5692099788303082*alpha[17]+0.4242640687119285*(alpha[14]+alpha[13])+0.4242640687119281*alpha[12]-0.4242640687119285*alpha[11]-0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*alpha[6]+0.6363961030678926*alpha[5]-0.6363961030678926*alpha[4]+0.4743416490252568*alpha[3]-0.4743416490252568*alpha[2]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[20] = ser_2x2v_p2_surfvy_quad_20(1, fl); 
    fUpwindQuad_r[20] = ser_2x2v_p2_surfvy_quad_20(1, fc); 
  } else { 

    fUpwindQuad_l[20] = ser_2x2v_p2_surfvy_quad_20(-1, fc); 
    fUpwindQuad_r[20] = ser_2x2v_p2_surfvy_quad_20(-1, fr); 
  } 
  if (0.711512473537885*alpha[18]-0.5303300858899104*alpha[14]+0.4242640687119285*alpha[13]+0.5303300858899104*alpha[12]-0.3952847075210473*alpha[8]+0.3162277660168379*alpha[7]-0.6363961030678926*alpha[5]+0.4743416490252568*alpha[3]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[21] = ser_2x2v_p2_surfvy_quad_21(1, fl); 
    fUpwindQuad_r[21] = ser_2x2v_p2_surfvy_quad_21(1, fc); 
  } else { 

    fUpwindQuad_l[21] = ser_2x2v_p2_surfvy_quad_21(-1, fc); 
    fUpwindQuad_r[21] = ser_2x2v_p2_surfvy_quad_21(-1, fr); 
  } 
  if ((-0.5303300858899104*(alpha[14]+alpha[13]))-0.3952847075210473*(alpha[8]+alpha[7])+0.4743416490252568*alpha[3]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[22] = ser_2x2v_p2_surfvy_quad_22(1, fl); 
    fUpwindQuad_r[22] = ser_2x2v_p2_surfvy_quad_22(1, fc); 
  } else { 

    fUpwindQuad_l[22] = ser_2x2v_p2_surfvy_quad_22(-1, fc); 
    fUpwindQuad_r[22] = ser_2x2v_p2_surfvy_quad_22(-1, fr); 
  } 
  if ((-0.711512473537885*alpha[18])-0.5303300858899104*alpha[14]+0.4242640687119285*alpha[13]-0.5303300858899104*alpha[12]-0.3952847075210473*alpha[8]+0.3162277660168379*alpha[7]+0.6363961030678926*alpha[5]+0.4743416490252568*(alpha[3]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[23] = ser_2x2v_p2_surfvy_quad_23(1, fl); 
    fUpwindQuad_r[23] = ser_2x2v_p2_surfvy_quad_23(1, fc); 
  } else { 

    fUpwindQuad_l[23] = ser_2x2v_p2_surfvy_quad_23(-1, fc); 
    fUpwindQuad_r[23] = ser_2x2v_p2_surfvy_quad_23(-1, fr); 
  } 
  if ((-0.5692099788303082*alpha[18])+0.5692099788303082*alpha[17]+0.4242640687119285*(alpha[14]+alpha[13])-0.4242640687119281*alpha[12]+0.4242640687119285*alpha[11]-0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*alpha[6]-0.6363961030678926*(alpha[5]+alpha[4])+0.4743416490252568*(alpha[3]+alpha[2])-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[24] = ser_2x2v_p2_surfvy_quad_24(1, fl); 
    fUpwindQuad_r[24] = ser_2x2v_p2_surfvy_quad_24(1, fc); 
  } else { 

    fUpwindQuad_l[24] = ser_2x2v_p2_surfvy_quad_24(-1, fc); 
    fUpwindQuad_r[24] = ser_2x2v_p2_surfvy_quad_24(-1, fr); 
  } 
  if ((-0.711512473537885*alpha[17])+0.4242640687119285*alpha[14]-0.5303300858899104*(alpha[13]+alpha[11])+0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]+0.6363961030678926*alpha[6]+0.4743416490252568*(alpha[3]+alpha[2])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[25] = ser_2x2v_p2_surfvy_quad_25(1, fl); 
    fUpwindQuad_r[25] = ser_2x2v_p2_surfvy_quad_25(1, fc); 
  } else { 

    fUpwindQuad_l[25] = ser_2x2v_p2_surfvy_quad_25(-1, fc); 
    fUpwindQuad_r[25] = ser_2x2v_p2_surfvy_quad_25(-1, fr); 
  } 
  if (0.5692099788303082*(alpha[18]+alpha[17])+0.4242640687119285*(alpha[14]+alpha[13])+0.4242640687119281*alpha[12]+0.4242640687119285*alpha[11]+0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*(alpha[6]+alpha[5]+alpha[4])+0.4743416490252568*(alpha[3]+alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[26] = ser_2x2v_p2_surfvy_quad_26(1, fl); 
    fUpwindQuad_r[26] = ser_2x2v_p2_surfvy_quad_26(1, fc); 
  } else { 

    fUpwindQuad_l[26] = ser_2x2v_p2_surfvy_quad_26(-1, fc); 
    fUpwindQuad_r[26] = ser_2x2v_p2_surfvy_quad_26(-1, fr); 
  } 
  fUpwind_l[0] = 0.06062300936098657*fUpwindQuad_l[26]+0.09699681497757856*fUpwindQuad_l[25]+0.06062300936098657*fUpwindQuad_l[24]+0.09699681497757856*fUpwindQuad_l[23]+0.1551949039641257*fUpwindQuad_l[22]+0.09699681497757856*fUpwindQuad_l[21]+0.06062300936098657*fUpwindQuad_l[20]+0.09699681497757856*fUpwindQuad_l[19]+0.06062300936098657*fUpwindQuad_l[18]+0.09699681497757856*fUpwindQuad_l[17]+0.1551949039641257*fUpwindQuad_l[16]+0.09699681497757856*fUpwindQuad_l[15]+0.1551949039641257*fUpwindQuad_l[14]+0.2483118463426013*fUpwindQuad_l[13]+0.1551949039641257*fUpwindQuad_l[12]+0.09699681497757856*fUpwindQuad_l[11]+0.1551949039641257*fUpwindQuad_l[10]+0.09699681497757856*fUpwindQuad_l[9]+0.06062300936098657*fUpwindQuad_l[8]+0.09699681497757856*fUpwindQuad_l[7]+0.06062300936098657*fUpwindQuad_l[6]+0.09699681497757856*fUpwindQuad_l[5]+0.1551949039641257*fUpwindQuad_l[4]+0.09699681497757856*fUpwindQuad_l[3]+0.06062300936098657*fUpwindQuad_l[2]+0.09699681497757856*fUpwindQuad_l[1]+0.06062300936098657*fUpwindQuad_l[0]; 
  fUpwind_l[1] = 0.08133430195906327*fUpwindQuad_l[26]-0.08133430195906327*fUpwindQuad_l[24]+0.1301348831345013*fUpwindQuad_l[23]-0.1301348831345013*fUpwindQuad_l[21]+0.08133430195906327*fUpwindQuad_l[20]-0.08133430195906327*fUpwindQuad_l[18]+0.1301348831345013*fUpwindQuad_l[17]-0.1301348831345013*fUpwindQuad_l[15]+0.2082158130152021*fUpwindQuad_l[14]-0.2082158130152021*fUpwindQuad_l[12]+0.1301348831345013*fUpwindQuad_l[11]-0.1301348831345013*fUpwindQuad_l[9]+0.08133430195906327*fUpwindQuad_l[8]-0.08133430195906327*fUpwindQuad_l[6]+0.1301348831345013*fUpwindQuad_l[5]-0.1301348831345013*fUpwindQuad_l[3]+0.08133430195906327*fUpwindQuad_l[2]-0.08133430195906327*fUpwindQuad_l[0]; 
  fUpwind_l[2] = 0.08133430195906327*fUpwindQuad_l[26]+0.1301348831345013*fUpwindQuad_l[25]+0.08133430195906327*fUpwindQuad_l[24]-0.08133430195906327*fUpwindQuad_l[20]-0.1301348831345013*fUpwindQuad_l[19]-0.08133430195906327*fUpwindQuad_l[18]+0.1301348831345013*fUpwindQuad_l[17]+0.2082158130152021*fUpwindQuad_l[16]+0.1301348831345013*fUpwindQuad_l[15]-0.1301348831345013*fUpwindQuad_l[11]-0.2082158130152021*fUpwindQuad_l[10]-0.1301348831345013*fUpwindQuad_l[9]+0.08133430195906327*fUpwindQuad_l[8]+0.1301348831345013*fUpwindQuad_l[7]+0.08133430195906327*fUpwindQuad_l[6]-0.08133430195906327*fUpwindQuad_l[2]-0.1301348831345013*fUpwindQuad_l[1]-0.08133430195906327*fUpwindQuad_l[0]; 
  fUpwind_l[3] = 0.08133430195906327*fUpwindQuad_l[26]+0.1301348831345013*fUpwindQuad_l[25]+0.08133430195906327*fUpwindQuad_l[24]+0.1301348831345013*fUpwindQuad_l[23]+0.2082158130152021*fUpwindQuad_l[22]+0.1301348831345013*fUpwindQuad_l[21]+0.08133430195906327*fUpwindQuad_l[20]+0.1301348831345013*fUpwindQuad_l[19]+0.08133430195906327*fUpwindQuad_l[18]-0.08133430195906327*fUpwindQuad_l[8]-0.1301348831345013*fUpwindQuad_l[7]-0.08133430195906327*fUpwindQuad_l[6]-0.1301348831345013*fUpwindQuad_l[5]-0.2082158130152021*fUpwindQuad_l[4]-0.1301348831345013*fUpwindQuad_l[3]-0.08133430195906327*fUpwindQuad_l[2]-0.1301348831345013*fUpwindQuad_l[1]-0.08133430195906327*fUpwindQuad_l[0]; 
  fUpwind_l[4] = 0.1091214168497758*fUpwindQuad_l[26]-0.1091214168497758*fUpwindQuad_l[24]-0.1091214168497758*fUpwindQuad_l[20]+0.1091214168497758*fUpwindQuad_l[18]+0.1745942669596414*fUpwindQuad_l[17]-0.1745942669596414*fUpwindQuad_l[15]-0.1745942669596414*fUpwindQuad_l[11]+0.1745942669596414*fUpwindQuad_l[9]+0.1091214168497758*fUpwindQuad_l[8]-0.1091214168497758*fUpwindQuad_l[6]-0.1091214168497758*fUpwindQuad_l[2]+0.1091214168497758*fUpwindQuad_l[0]; 
  fUpwind_l[5] = 0.1091214168497758*fUpwindQuad_l[26]-0.1091214168497758*fUpwindQuad_l[24]+0.1745942669596414*fUpwindQuad_l[23]-0.1745942669596414*fUpwindQuad_l[21]+0.1091214168497758*fUpwindQuad_l[20]-0.1091214168497758*fUpwindQuad_l[18]-0.1091214168497758*fUpwindQuad_l[8]+0.1091214168497758*fUpwindQuad_l[6]-0.1745942669596414*fUpwindQuad_l[5]+0.1745942669596414*fUpwindQuad_l[3]-0.1091214168497758*fUpwindQuad_l[2]+0.1091214168497758*fUpwindQuad_l[0]; 
  fUpwind_l[6] = 0.1091214168497758*fUpwindQuad_l[26]+0.1745942669596414*fUpwindQuad_l[25]+0.1091214168497758*fUpwindQuad_l[24]-0.1091214168497758*fUpwindQuad_l[20]-0.1745942669596414*fUpwindQuad_l[19]-0.1091214168497758*fUpwindQuad_l[18]-0.1091214168497758*fUpwindQuad_l[8]-0.1745942669596414*fUpwindQuad_l[7]-0.1091214168497758*fUpwindQuad_l[6]+0.1091214168497758*fUpwindQuad_l[2]+0.1745942669596414*fUpwindQuad_l[1]+0.1091214168497758*fUpwindQuad_l[0]; 
  fUpwind_l[7] = 0.05422286797270884*fUpwindQuad_l[26]-0.1084457359454177*fUpwindQuad_l[25]+0.05422286797270884*fUpwindQuad_l[24]+0.08675658875633419*fUpwindQuad_l[23]-0.1735131775126684*fUpwindQuad_l[22]+0.08675658875633419*fUpwindQuad_l[21]+0.05422286797270884*fUpwindQuad_l[20]-0.1084457359454177*fUpwindQuad_l[19]+0.05422286797270884*fUpwindQuad_l[18]+0.08675658875633419*fUpwindQuad_l[17]-0.1735131775126684*fUpwindQuad_l[16]+0.08675658875633419*fUpwindQuad_l[15]+0.1388105420101347*fUpwindQuad_l[14]-0.2776210840202695*fUpwindQuad_l[13]+0.1388105420101347*fUpwindQuad_l[12]+0.08675658875633419*fUpwindQuad_l[11]-0.1735131775126684*fUpwindQuad_l[10]+0.08675658875633419*fUpwindQuad_l[9]+0.05422286797270884*fUpwindQuad_l[8]-0.1084457359454177*fUpwindQuad_l[7]+0.05422286797270884*fUpwindQuad_l[6]+0.08675658875633419*fUpwindQuad_l[5]-0.1735131775126684*fUpwindQuad_l[4]+0.08675658875633419*fUpwindQuad_l[3]+0.05422286797270884*fUpwindQuad_l[2]-0.1084457359454177*fUpwindQuad_l[1]+0.05422286797270884*fUpwindQuad_l[0]; 
  fUpwind_l[8] = 0.05422286797270884*fUpwindQuad_l[26]+0.08675658875633419*fUpwindQuad_l[25]+0.05422286797270884*fUpwindQuad_l[24]-0.1084457359454177*fUpwindQuad_l[23]-0.1735131775126684*fUpwindQuad_l[22]-0.1084457359454177*fUpwindQuad_l[21]+0.05422286797270884*fUpwindQuad_l[20]+0.08675658875633419*fUpwindQuad_l[19]+0.05422286797270884*fUpwindQuad_l[18]+0.08675658875633419*fUpwindQuad_l[17]+0.1388105420101347*fUpwindQuad_l[16]+0.08675658875633419*fUpwindQuad_l[15]-0.1735131775126684*fUpwindQuad_l[14]-0.2776210840202695*fUpwindQuad_l[13]-0.1735131775126684*fUpwindQuad_l[12]+0.08675658875633419*fUpwindQuad_l[11]+0.1388105420101347*fUpwindQuad_l[10]+0.08675658875633419*fUpwindQuad_l[9]+0.05422286797270884*fUpwindQuad_l[8]+0.08675658875633419*fUpwindQuad_l[7]+0.05422286797270884*fUpwindQuad_l[6]-0.1084457359454177*fUpwindQuad_l[5]-0.1735131775126684*fUpwindQuad_l[4]-0.1084457359454177*fUpwindQuad_l[3]+0.05422286797270884*fUpwindQuad_l[2]+0.08675658875633419*fUpwindQuad_l[1]+0.05422286797270884*fUpwindQuad_l[0]; 
  fUpwind_l[9] = 0.05422286797270884*fUpwindQuad_l[26]+0.08675658875633419*fUpwindQuad_l[25]+0.05422286797270884*fUpwindQuad_l[24]+0.08675658875633419*fUpwindQuad_l[23]+0.1388105420101347*fUpwindQuad_l[22]+0.08675658875633419*fUpwindQuad_l[21]+0.05422286797270884*fUpwindQuad_l[20]+0.08675658875633419*fUpwindQuad_l[19]+0.05422286797270884*fUpwindQuad_l[18]-0.1084457359454177*fUpwindQuad_l[17]-0.1735131775126684*fUpwindQuad_l[16]-0.1084457359454177*fUpwindQuad_l[15]-0.1735131775126684*fUpwindQuad_l[14]-0.2776210840202695*fUpwindQuad_l[13]-0.1735131775126684*fUpwindQuad_l[12]-0.1084457359454177*fUpwindQuad_l[11]-0.1735131775126684*fUpwindQuad_l[10]-0.1084457359454177*fUpwindQuad_l[9]+0.05422286797270884*fUpwindQuad_l[8]+0.08675658875633419*fUpwindQuad_l[7]+0.05422286797270884*fUpwindQuad_l[6]+0.08675658875633419*fUpwindQuad_l[5]+0.1388105420101347*fUpwindQuad_l[4]+0.08675658875633419*fUpwindQuad_l[3]+0.05422286797270884*fUpwindQuad_l[2]+0.08675658875633419*fUpwindQuad_l[1]+0.05422286797270884*fUpwindQuad_l[0]; 
  fUpwind_l[10] = 0.1464017435263139*fUpwindQuad_l[26]-0.1464017435263139*fUpwindQuad_l[24]-0.1464017435263139*fUpwindQuad_l[20]+0.1464017435263139*fUpwindQuad_l[18]-0.1464017435263139*fUpwindQuad_l[8]+0.1464017435263139*fUpwindQuad_l[6]+0.1464017435263139*fUpwindQuad_l[2]-0.1464017435263139*fUpwindQuad_l[0]; 
  fUpwind_l[11] = 0.07274761123318395*fUpwindQuad_l[26]-0.1454952224663679*fUpwindQuad_l[25]+0.07274761123318395*fUpwindQuad_l[24]-0.07274761123318395*fUpwindQuad_l[20]+0.1454952224663679*fUpwindQuad_l[19]-0.07274761123318395*fUpwindQuad_l[18]+0.1163961779730944*fUpwindQuad_l[17]-0.2327923559461888*fUpwindQuad_l[16]+0.1163961779730944*fUpwindQuad_l[15]-0.1163961779730944*fUpwindQuad_l[11]+0.2327923559461888*fUpwindQuad_l[10]-0.1163961779730944*fUpwindQuad_l[9]+0.07274761123318395*fUpwindQuad_l[8]-0.1454952224663679*fUpwindQuad_l[7]+0.07274761123318395*fUpwindQuad_l[6]-0.07274761123318395*fUpwindQuad_l[2]+0.1454952224663679*fUpwindQuad_l[1]-0.07274761123318395*fUpwindQuad_l[0]; 
  fUpwind_l[12] = 0.07274761123318395*fUpwindQuad_l[26]-0.07274761123318395*fUpwindQuad_l[24]-0.1454952224663679*fUpwindQuad_l[23]+0.1454952224663679*fUpwindQuad_l[21]+0.07274761123318395*fUpwindQuad_l[20]-0.07274761123318395*fUpwindQuad_l[18]+0.1163961779730944*fUpwindQuad_l[17]-0.1163961779730944*fUpwindQuad_l[15]-0.2327923559461888*fUpwindQuad_l[14]+0.2327923559461888*fUpwindQuad_l[12]+0.1163961779730944*fUpwindQuad_l[11]-0.1163961779730944*fUpwindQuad_l[9]+0.07274761123318395*fUpwindQuad_l[8]-0.07274761123318395*fUpwindQuad_l[6]-0.1454952224663679*fUpwindQuad_l[5]+0.1454952224663679*fUpwindQuad_l[3]+0.07274761123318395*fUpwindQuad_l[2]-0.07274761123318395*fUpwindQuad_l[0]; 
  fUpwind_l[13] = 0.07274761123318395*fUpwindQuad_l[26]-0.1454952224663679*fUpwindQuad_l[25]+0.07274761123318395*fUpwindQuad_l[24]+0.1163961779730944*fUpwindQuad_l[23]-0.2327923559461888*fUpwindQuad_l[22]+0.1163961779730944*fUpwindQuad_l[21]+0.07274761123318395*fUpwindQuad_l[20]-0.1454952224663679*fUpwindQuad_l[19]+0.07274761123318395*fUpwindQuad_l[18]-0.07274761123318395*fUpwindQuad_l[8]+0.1454952224663679*fUpwindQuad_l[7]-0.07274761123318395*fUpwindQuad_l[6]-0.1163961779730944*fUpwindQuad_l[5]+0.2327923559461888*fUpwindQuad_l[4]-0.1163961779730944*fUpwindQuad_l[3]-0.07274761123318395*fUpwindQuad_l[2]+0.1454952224663679*fUpwindQuad_l[1]-0.07274761123318395*fUpwindQuad_l[0]; 
  fUpwind_l[14] = 0.07274761123318395*fUpwindQuad_l[26]+0.1163961779730944*fUpwindQuad_l[25]+0.07274761123318395*fUpwindQuad_l[24]-0.1454952224663679*fUpwindQuad_l[23]-0.2327923559461888*fUpwindQuad_l[22]-0.1454952224663679*fUpwindQuad_l[21]+0.07274761123318395*fUpwindQuad_l[20]+0.1163961779730944*fUpwindQuad_l[19]+0.07274761123318395*fUpwindQuad_l[18]-0.07274761123318395*fUpwindQuad_l[8]-0.1163961779730944*fUpwindQuad_l[7]-0.07274761123318395*fUpwindQuad_l[6]+0.1454952224663679*fUpwindQuad_l[5]+0.2327923559461888*fUpwindQuad_l[4]+0.1454952224663679*fUpwindQuad_l[3]-0.07274761123318395*fUpwindQuad_l[2]-0.1163961779730944*fUpwindQuad_l[1]-0.07274761123318395*fUpwindQuad_l[0]; 
  fUpwind_l[15] = 0.07274761123318395*fUpwindQuad_l[26]-0.07274761123318395*fUpwindQuad_l[24]+0.1163961779730944*fUpwindQuad_l[23]-0.1163961779730944*fUpwindQuad_l[21]+0.07274761123318395*fUpwindQuad_l[20]-0.07274761123318395*fUpwindQuad_l[18]-0.1454952224663679*fUpwindQuad_l[17]+0.1454952224663679*fUpwindQuad_l[15]-0.2327923559461888*fUpwindQuad_l[14]+0.2327923559461888*fUpwindQuad_l[12]-0.1454952224663679*fUpwindQuad_l[11]+0.1454952224663679*fUpwindQuad_l[9]+0.07274761123318395*fUpwindQuad_l[8]-0.07274761123318395*fUpwindQuad_l[6]+0.1163961779730944*fUpwindQuad_l[5]-0.1163961779730944*fUpwindQuad_l[3]+0.07274761123318395*fUpwindQuad_l[2]-0.07274761123318395*fUpwindQuad_l[0]; 
  fUpwind_l[16] = 0.07274761123318395*fUpwindQuad_l[26]+0.1163961779730944*fUpwindQuad_l[25]+0.07274761123318395*fUpwindQuad_l[24]-0.07274761123318395*fUpwindQuad_l[20]-0.1163961779730944*fUpwindQuad_l[19]-0.07274761123318395*fUpwindQuad_l[18]-0.1454952224663679*fUpwindQuad_l[17]-0.2327923559461888*fUpwindQuad_l[16]-0.1454952224663679*fUpwindQuad_l[15]+0.1454952224663679*fUpwindQuad_l[11]+0.2327923559461888*fUpwindQuad_l[10]+0.1454952224663679*fUpwindQuad_l[9]+0.07274761123318395*fUpwindQuad_l[8]+0.1163961779730944*fUpwindQuad_l[7]+0.07274761123318395*fUpwindQuad_l[6]-0.07274761123318395*fUpwindQuad_l[2]-0.1163961779730944*fUpwindQuad_l[1]-0.07274761123318395*fUpwindQuad_l[0]; 
  fUpwind_l[17] = 0.09760116235087592*fUpwindQuad_l[26]-0.1952023247017519*fUpwindQuad_l[25]+0.09760116235087592*fUpwindQuad_l[24]-0.09760116235087592*fUpwindQuad_l[20]+0.1952023247017519*fUpwindQuad_l[19]-0.09760116235087592*fUpwindQuad_l[18]-0.09760116235087592*fUpwindQuad_l[8]+0.1952023247017519*fUpwindQuad_l[7]-0.09760116235087592*fUpwindQuad_l[6]+0.09760116235087592*fUpwindQuad_l[2]-0.1952023247017519*fUpwindQuad_l[1]+0.09760116235087592*fUpwindQuad_l[0]; 
  fUpwind_l[18] = 0.09760116235087592*fUpwindQuad_l[26]-0.09760116235087592*fUpwindQuad_l[24]-0.1952023247017519*fUpwindQuad_l[23]+0.1952023247017519*fUpwindQuad_l[21]+0.09760116235087592*fUpwindQuad_l[20]-0.09760116235087592*fUpwindQuad_l[18]-0.09760116235087592*fUpwindQuad_l[8]+0.09760116235087592*fUpwindQuad_l[6]+0.1952023247017519*fUpwindQuad_l[5]-0.1952023247017519*fUpwindQuad_l[3]-0.09760116235087592*fUpwindQuad_l[2]+0.09760116235087592*fUpwindQuad_l[0]; 
  fUpwind_l[19] = 0.09760116235087592*fUpwindQuad_l[26]-0.09760116235087592*fUpwindQuad_l[24]-0.09760116235087592*fUpwindQuad_l[20]+0.09760116235087592*fUpwindQuad_l[18]-0.1952023247017519*fUpwindQuad_l[17]+0.1952023247017519*fUpwindQuad_l[15]+0.1952023247017519*fUpwindQuad_l[11]-0.1952023247017519*fUpwindQuad_l[9]+0.09760116235087592*fUpwindQuad_l[8]-0.09760116235087592*fUpwindQuad_l[6]-0.09760116235087592*fUpwindQuad_l[2]+0.09760116235087592*fUpwindQuad_l[0]; 
  fUpwind_r[0] = 0.06062300936098657*fUpwindQuad_r[26]+0.09699681497757856*fUpwindQuad_r[25]+0.06062300936098657*fUpwindQuad_r[24]+0.09699681497757856*fUpwindQuad_r[23]+0.1551949039641257*fUpwindQuad_r[22]+0.09699681497757856*fUpwindQuad_r[21]+0.06062300936098657*fUpwindQuad_r[20]+0.09699681497757856*fUpwindQuad_r[19]+0.06062300936098657*fUpwindQuad_r[18]+0.09699681497757856*fUpwindQuad_r[17]+0.1551949039641257*fUpwindQuad_r[16]+0.09699681497757856*fUpwindQuad_r[15]+0.1551949039641257*fUpwindQuad_r[14]+0.2483118463426013*fUpwindQuad_r[13]+0.1551949039641257*fUpwindQuad_r[12]+0.09699681497757856*fUpwindQuad_r[11]+0.1551949039641257*fUpwindQuad_r[10]+0.09699681497757856*fUpwindQuad_r[9]+0.06062300936098657*fUpwindQuad_r[8]+0.09699681497757856*fUpwindQuad_r[7]+0.06062300936098657*fUpwindQuad_r[6]+0.09699681497757856*fUpwindQuad_r[5]+0.1551949039641257*fUpwindQuad_r[4]+0.09699681497757856*fUpwindQuad_r[3]+0.06062300936098657*fUpwindQuad_r[2]+0.09699681497757856*fUpwindQuad_r[1]+0.06062300936098657*fUpwindQuad_r[0]; 
  fUpwind_r[1] = 0.08133430195906327*fUpwindQuad_r[26]-0.08133430195906327*fUpwindQuad_r[24]+0.1301348831345013*fUpwindQuad_r[23]-0.1301348831345013*fUpwindQuad_r[21]+0.08133430195906327*fUpwindQuad_r[20]-0.08133430195906327*fUpwindQuad_r[18]+0.1301348831345013*fUpwindQuad_r[17]-0.1301348831345013*fUpwindQuad_r[15]+0.2082158130152021*fUpwindQuad_r[14]-0.2082158130152021*fUpwindQuad_r[12]+0.1301348831345013*fUpwindQuad_r[11]-0.1301348831345013*fUpwindQuad_r[9]+0.08133430195906327*fUpwindQuad_r[8]-0.08133430195906327*fUpwindQuad_r[6]+0.1301348831345013*fUpwindQuad_r[5]-0.1301348831345013*fUpwindQuad_r[3]+0.08133430195906327*fUpwindQuad_r[2]-0.08133430195906327*fUpwindQuad_r[0]; 
  fUpwind_r[2] = 0.08133430195906327*fUpwindQuad_r[26]+0.1301348831345013*fUpwindQuad_r[25]+0.08133430195906327*fUpwindQuad_r[24]-0.08133430195906327*fUpwindQuad_r[20]-0.1301348831345013*fUpwindQuad_r[19]-0.08133430195906327*fUpwindQuad_r[18]+0.1301348831345013*fUpwindQuad_r[17]+0.2082158130152021*fUpwindQuad_r[16]+0.1301348831345013*fUpwindQuad_r[15]-0.1301348831345013*fUpwindQuad_r[11]-0.2082158130152021*fUpwindQuad_r[10]-0.1301348831345013*fUpwindQuad_r[9]+0.08133430195906327*fUpwindQuad_r[8]+0.1301348831345013*fUpwindQuad_r[7]+0.08133430195906327*fUpwindQuad_r[6]-0.08133430195906327*fUpwindQuad_r[2]-0.1301348831345013*fUpwindQuad_r[1]-0.08133430195906327*fUpwindQuad_r[0]; 
  fUpwind_r[3] = 0.08133430195906327*fUpwindQuad_r[26]+0.1301348831345013*fUpwindQuad_r[25]+0.08133430195906327*fUpwindQuad_r[24]+0.1301348831345013*fUpwindQuad_r[23]+0.2082158130152021*fUpwindQuad_r[22]+0.1301348831345013*fUpwindQuad_r[21]+0.08133430195906327*fUpwindQuad_r[20]+0.1301348831345013*fUpwindQuad_r[19]+0.08133430195906327*fUpwindQuad_r[18]-0.08133430195906327*fUpwindQuad_r[8]-0.1301348831345013*fUpwindQuad_r[7]-0.08133430195906327*fUpwindQuad_r[6]-0.1301348831345013*fUpwindQuad_r[5]-0.2082158130152021*fUpwindQuad_r[4]-0.1301348831345013*fUpwindQuad_r[3]-0.08133430195906327*fUpwindQuad_r[2]-0.1301348831345013*fUpwindQuad_r[1]-0.08133430195906327*fUpwindQuad_r[0]; 
  fUpwind_r[4] = 0.1091214168497758*fUpwindQuad_r[26]-0.1091214168497758*fUpwindQuad_r[24]-0.1091214168497758*fUpwindQuad_r[20]+0.1091214168497758*fUpwindQuad_r[18]+0.1745942669596414*fUpwindQuad_r[17]-0.1745942669596414*fUpwindQuad_r[15]-0.1745942669596414*fUpwindQuad_r[11]+0.1745942669596414*fUpwindQuad_r[9]+0.1091214168497758*fUpwindQuad_r[8]-0.1091214168497758*fUpwindQuad_r[6]-0.1091214168497758*fUpwindQuad_r[2]+0.1091214168497758*fUpwindQuad_r[0]; 
  fUpwind_r[5] = 0.1091214168497758*fUpwindQuad_r[26]-0.1091214168497758*fUpwindQuad_r[24]+0.1745942669596414*fUpwindQuad_r[23]-0.1745942669596414*fUpwindQuad_r[21]+0.1091214168497758*fUpwindQuad_r[20]-0.1091214168497758*fUpwindQuad_r[18]-0.1091214168497758*fUpwindQuad_r[8]+0.1091214168497758*fUpwindQuad_r[6]-0.1745942669596414*fUpwindQuad_r[5]+0.1745942669596414*fUpwindQuad_r[3]-0.1091214168497758*fUpwindQuad_r[2]+0.1091214168497758*fUpwindQuad_r[0]; 
  fUpwind_r[6] = 0.1091214168497758*fUpwindQuad_r[26]+0.1745942669596414*fUpwindQuad_r[25]+0.1091214168497758*fUpwindQuad_r[24]-0.1091214168497758*fUpwindQuad_r[20]-0.1745942669596414*fUpwindQuad_r[19]-0.1091214168497758*fUpwindQuad_r[18]-0.1091214168497758*fUpwindQuad_r[8]-0.1745942669596414*fUpwindQuad_r[7]-0.1091214168497758*fUpwindQuad_r[6]+0.1091214168497758*fUpwindQuad_r[2]+0.1745942669596414*fUpwindQuad_r[1]+0.1091214168497758*fUpwindQuad_r[0]; 
  fUpwind_r[7] = 0.05422286797270884*fUpwindQuad_r[26]-0.1084457359454177*fUpwindQuad_r[25]+0.05422286797270884*fUpwindQuad_r[24]+0.08675658875633419*fUpwindQuad_r[23]-0.1735131775126684*fUpwindQuad_r[22]+0.08675658875633419*fUpwindQuad_r[21]+0.05422286797270884*fUpwindQuad_r[20]-0.1084457359454177*fUpwindQuad_r[19]+0.05422286797270884*fUpwindQuad_r[18]+0.08675658875633419*fUpwindQuad_r[17]-0.1735131775126684*fUpwindQuad_r[16]+0.08675658875633419*fUpwindQuad_r[15]+0.1388105420101347*fUpwindQuad_r[14]-0.2776210840202695*fUpwindQuad_r[13]+0.1388105420101347*fUpwindQuad_r[12]+0.08675658875633419*fUpwindQuad_r[11]-0.1735131775126684*fUpwindQuad_r[10]+0.08675658875633419*fUpwindQuad_r[9]+0.05422286797270884*fUpwindQuad_r[8]-0.1084457359454177*fUpwindQuad_r[7]+0.05422286797270884*fUpwindQuad_r[6]+0.08675658875633419*fUpwindQuad_r[5]-0.1735131775126684*fUpwindQuad_r[4]+0.08675658875633419*fUpwindQuad_r[3]+0.05422286797270884*fUpwindQuad_r[2]-0.1084457359454177*fUpwindQuad_r[1]+0.05422286797270884*fUpwindQuad_r[0]; 
  fUpwind_r[8] = 0.05422286797270884*fUpwindQuad_r[26]+0.08675658875633419*fUpwindQuad_r[25]+0.05422286797270884*fUpwindQuad_r[24]-0.1084457359454177*fUpwindQuad_r[23]-0.1735131775126684*fUpwindQuad_r[22]-0.1084457359454177*fUpwindQuad_r[21]+0.05422286797270884*fUpwindQuad_r[20]+0.08675658875633419*fUpwindQuad_r[19]+0.05422286797270884*fUpwindQuad_r[18]+0.08675658875633419*fUpwindQuad_r[17]+0.1388105420101347*fUpwindQuad_r[16]+0.08675658875633419*fUpwindQuad_r[15]-0.1735131775126684*fUpwindQuad_r[14]-0.2776210840202695*fUpwindQuad_r[13]-0.1735131775126684*fUpwindQuad_r[12]+0.08675658875633419*fUpwindQuad_r[11]+0.1388105420101347*fUpwindQuad_r[10]+0.08675658875633419*fUpwindQuad_r[9]+0.05422286797270884*fUpwindQuad_r[8]+0.08675658875633419*fUpwindQuad_r[7]+0.05422286797270884*fUpwindQuad_r[6]-0.1084457359454177*fUpwindQuad_r[5]-0.1735131775126684*fUpwindQuad_r[4]-0.1084457359454177*fUpwindQuad_r[3]+0.05422286797270884*fUpwindQuad_r[2]+0.08675658875633419*fUpwindQuad_r[1]+0.05422286797270884*fUpwindQuad_r[0]; 
  fUpwind_r[9] = 0.05422286797270884*fUpwindQuad_r[26]+0.08675658875633419*fUpwindQuad_r[25]+0.05422286797270884*fUpwindQuad_r[24]+0.08675658875633419*fUpwindQuad_r[23]+0.1388105420101347*fUpwindQuad_r[22]+0.08675658875633419*fUpwindQuad_r[21]+0.05422286797270884*fUpwindQuad_r[20]+0.08675658875633419*fUpwindQuad_r[19]+0.05422286797270884*fUpwindQuad_r[18]-0.1084457359454177*fUpwindQuad_r[17]-0.1735131775126684*fUpwindQuad_r[16]-0.1084457359454177*fUpwindQuad_r[15]-0.1735131775126684*fUpwindQuad_r[14]-0.2776210840202695*fUpwindQuad_r[13]-0.1735131775126684*fUpwindQuad_r[12]-0.1084457359454177*fUpwindQuad_r[11]-0.1735131775126684*fUpwindQuad_r[10]-0.1084457359454177*fUpwindQuad_r[9]+0.05422286797270884*fUpwindQuad_r[8]+0.08675658875633419*fUpwindQuad_r[7]+0.05422286797270884*fUpwindQuad_r[6]+0.08675658875633419*fUpwindQuad_r[5]+0.1388105420101347*fUpwindQuad_r[4]+0.08675658875633419*fUpwindQuad_r[3]+0.05422286797270884*fUpwindQuad_r[2]+0.08675658875633419*fUpwindQuad_r[1]+0.05422286797270884*fUpwindQuad_r[0]; 
  fUpwind_r[10] = 0.1464017435263139*fUpwindQuad_r[26]-0.1464017435263139*fUpwindQuad_r[24]-0.1464017435263139*fUpwindQuad_r[20]+0.1464017435263139*fUpwindQuad_r[18]-0.1464017435263139*fUpwindQuad_r[8]+0.1464017435263139*fUpwindQuad_r[6]+0.1464017435263139*fUpwindQuad_r[2]-0.1464017435263139*fUpwindQuad_r[0]; 
  fUpwind_r[11] = 0.07274761123318395*fUpwindQuad_r[26]-0.1454952224663679*fUpwindQuad_r[25]+0.07274761123318395*fUpwindQuad_r[24]-0.07274761123318395*fUpwindQuad_r[20]+0.1454952224663679*fUpwindQuad_r[19]-0.07274761123318395*fUpwindQuad_r[18]+0.1163961779730944*fUpwindQuad_r[17]-0.2327923559461888*fUpwindQuad_r[16]+0.1163961779730944*fUpwindQuad_r[15]-0.1163961779730944*fUpwindQuad_r[11]+0.2327923559461888*fUpwindQuad_r[10]-0.1163961779730944*fUpwindQuad_r[9]+0.07274761123318395*fUpwindQuad_r[8]-0.1454952224663679*fUpwindQuad_r[7]+0.07274761123318395*fUpwindQuad_r[6]-0.07274761123318395*fUpwindQuad_r[2]+0.1454952224663679*fUpwindQuad_r[1]-0.07274761123318395*fUpwindQuad_r[0]; 
  fUpwind_r[12] = 0.07274761123318395*fUpwindQuad_r[26]-0.07274761123318395*fUpwindQuad_r[24]-0.1454952224663679*fUpwindQuad_r[23]+0.1454952224663679*fUpwindQuad_r[21]+0.07274761123318395*fUpwindQuad_r[20]-0.07274761123318395*fUpwindQuad_r[18]+0.1163961779730944*fUpwindQuad_r[17]-0.1163961779730944*fUpwindQuad_r[15]-0.2327923559461888*fUpwindQuad_r[14]+0.2327923559461888*fUpwindQuad_r[12]+0.1163961779730944*fUpwindQuad_r[11]-0.1163961779730944*fUpwindQuad_r[9]+0.07274761123318395*fUpwindQuad_r[8]-0.07274761123318395*fUpwindQuad_r[6]-0.1454952224663679*fUpwindQuad_r[5]+0.1454952224663679*fUpwindQuad_r[3]+0.07274761123318395*fUpwindQuad_r[2]-0.07274761123318395*fUpwindQuad_r[0]; 
  fUpwind_r[13] = 0.07274761123318395*fUpwindQuad_r[26]-0.1454952224663679*fUpwindQuad_r[25]+0.07274761123318395*fUpwindQuad_r[24]+0.1163961779730944*fUpwindQuad_r[23]-0.2327923559461888*fUpwindQuad_r[22]+0.1163961779730944*fUpwindQuad_r[21]+0.07274761123318395*fUpwindQuad_r[20]-0.1454952224663679*fUpwindQuad_r[19]+0.07274761123318395*fUpwindQuad_r[18]-0.07274761123318395*fUpwindQuad_r[8]+0.1454952224663679*fUpwindQuad_r[7]-0.07274761123318395*fUpwindQuad_r[6]-0.1163961779730944*fUpwindQuad_r[5]+0.2327923559461888*fUpwindQuad_r[4]-0.1163961779730944*fUpwindQuad_r[3]-0.07274761123318395*fUpwindQuad_r[2]+0.1454952224663679*fUpwindQuad_r[1]-0.07274761123318395*fUpwindQuad_r[0]; 
  fUpwind_r[14] = 0.07274761123318395*fUpwindQuad_r[26]+0.1163961779730944*fUpwindQuad_r[25]+0.07274761123318395*fUpwindQuad_r[24]-0.1454952224663679*fUpwindQuad_r[23]-0.2327923559461888*fUpwindQuad_r[22]-0.1454952224663679*fUpwindQuad_r[21]+0.07274761123318395*fUpwindQuad_r[20]+0.1163961779730944*fUpwindQuad_r[19]+0.07274761123318395*fUpwindQuad_r[18]-0.07274761123318395*fUpwindQuad_r[8]-0.1163961779730944*fUpwindQuad_r[7]-0.07274761123318395*fUpwindQuad_r[6]+0.1454952224663679*fUpwindQuad_r[5]+0.2327923559461888*fUpwindQuad_r[4]+0.1454952224663679*fUpwindQuad_r[3]-0.07274761123318395*fUpwindQuad_r[2]-0.1163961779730944*fUpwindQuad_r[1]-0.07274761123318395*fUpwindQuad_r[0]; 
  fUpwind_r[15] = 0.07274761123318395*fUpwindQuad_r[26]-0.07274761123318395*fUpwindQuad_r[24]+0.1163961779730944*fUpwindQuad_r[23]-0.1163961779730944*fUpwindQuad_r[21]+0.07274761123318395*fUpwindQuad_r[20]-0.07274761123318395*fUpwindQuad_r[18]-0.1454952224663679*fUpwindQuad_r[17]+0.1454952224663679*fUpwindQuad_r[15]-0.2327923559461888*fUpwindQuad_r[14]+0.2327923559461888*fUpwindQuad_r[12]-0.1454952224663679*fUpwindQuad_r[11]+0.1454952224663679*fUpwindQuad_r[9]+0.07274761123318395*fUpwindQuad_r[8]-0.07274761123318395*fUpwindQuad_r[6]+0.1163961779730944*fUpwindQuad_r[5]-0.1163961779730944*fUpwindQuad_r[3]+0.07274761123318395*fUpwindQuad_r[2]-0.07274761123318395*fUpwindQuad_r[0]; 
  fUpwind_r[16] = 0.07274761123318395*fUpwindQuad_r[26]+0.1163961779730944*fUpwindQuad_r[25]+0.07274761123318395*fUpwindQuad_r[24]-0.07274761123318395*fUpwindQuad_r[20]-0.1163961779730944*fUpwindQuad_r[19]-0.07274761123318395*fUpwindQuad_r[18]-0.1454952224663679*fUpwindQuad_r[17]-0.2327923559461888*fUpwindQuad_r[16]-0.1454952224663679*fUpwindQuad_r[15]+0.1454952224663679*fUpwindQuad_r[11]+0.2327923559461888*fUpwindQuad_r[10]+0.1454952224663679*fUpwindQuad_r[9]+0.07274761123318395*fUpwindQuad_r[8]+0.1163961779730944*fUpwindQuad_r[7]+0.07274761123318395*fUpwindQuad_r[6]-0.07274761123318395*fUpwindQuad_r[2]-0.1163961779730944*fUpwindQuad_r[1]-0.07274761123318395*fUpwindQuad_r[0]; 
  fUpwind_r[17] = 0.09760116235087592*fUpwindQuad_r[26]-0.1952023247017519*fUpwindQuad_r[25]+0.09760116235087592*fUpwindQuad_r[24]-0.09760116235087592*fUpwindQuad_r[20]+0.1952023247017519*fUpwindQuad_r[19]-0.09760116235087592*fUpwindQuad_r[18]-0.09760116235087592*fUpwindQuad_r[8]+0.1952023247017519*fUpwindQuad_r[7]-0.09760116235087592*fUpwindQuad_r[6]+0.09760116235087592*fUpwindQuad_r[2]-0.1952023247017519*fUpwindQuad_r[1]+0.09760116235087592*fUpwindQuad_r[0]; 
  fUpwind_r[18] = 0.09760116235087592*fUpwindQuad_r[26]-0.09760116235087592*fUpwindQuad_r[24]-0.1952023247017519*fUpwindQuad_r[23]+0.1952023247017519*fUpwindQuad_r[21]+0.09760116235087592*fUpwindQuad_r[20]-0.09760116235087592*fUpwindQuad_r[18]-0.09760116235087592*fUpwindQuad_r[8]+0.09760116235087592*fUpwindQuad_r[6]+0.1952023247017519*fUpwindQuad_r[5]-0.1952023247017519*fUpwindQuad_r[3]-0.09760116235087592*fUpwindQuad_r[2]+0.09760116235087592*fUpwindQuad_r[0]; 
  fUpwind_r[19] = 0.09760116235087592*fUpwindQuad_r[26]-0.09760116235087592*fUpwindQuad_r[24]-0.09760116235087592*fUpwindQuad_r[20]+0.09760116235087592*fUpwindQuad_r[18]-0.1952023247017519*fUpwindQuad_r[17]+0.1952023247017519*fUpwindQuad_r[15]+0.1952023247017519*fUpwindQuad_r[11]-0.1952023247017519*fUpwindQuad_r[9]+0.09760116235087592*fUpwindQuad_r[8]-0.09760116235087592*fUpwindQuad_r[6]-0.09760116235087592*fUpwindQuad_r[2]+0.09760116235087592*fUpwindQuad_r[0]; 

  Ghat_l[0] += 0.3535533905932737*(alpha[18]*fUpwind_l[18]+alpha[17]*fUpwind_l[17]+alpha[14]*fUpwind_l[14]+alpha[13]*fUpwind_l[13]+alpha[12]*fUpwind_l[12]+alpha[11]*fUpwind_l[11]+alpha[10]*fUpwind_l[10]+alpha[8]*fUpwind_l[8]+alpha[7]*fUpwind_l[7]+alpha[6]*fUpwind_l[6]+alpha[5]*fUpwind_l[5]+alpha[4]*fUpwind_l[4]+alpha[3]*fUpwind_l[3]+alpha[2]*fUpwind_l[2]+alpha[1]*fUpwind_l[1]+alpha[0]*fUpwind_l[0]); 
  Ghat_l[1] += 0.3535533905932737*(alpha[14]*fUpwind_l[18]+fUpwind_l[14]*alpha[18])+0.3162277660168379*(alpha[10]*fUpwind_l[17]+fUpwind_l[10]*alpha[17])+0.3162277660168379*(alpha[5]*fUpwind_l[13]+fUpwind_l[5]*alpha[13])+0.3535533905932737*(alpha[8]*fUpwind_l[12]+fUpwind_l[8]*alpha[12])+0.3162277660168379*(alpha[4]*fUpwind_l[11]+fUpwind_l[4]*alpha[11])+0.3535533905932737*(alpha[6]*fUpwind_l[10]+fUpwind_l[6]*alpha[10])+0.3162277660168379*(alpha[1]*fUpwind_l[7]+fUpwind_l[1]*alpha[7])+0.3535533905932737*(alpha[3]*fUpwind_l[5]+fUpwind_l[3]*alpha[5]+alpha[2]*fUpwind_l[4]+fUpwind_l[2]*alpha[4]+alpha[0]*fUpwind_l[1]+fUpwind_l[0]*alpha[1]); 
  Ghat_l[2] += 0.3162277660168379*(alpha[10]*fUpwind_l[18]+fUpwind_l[10]*alpha[18])+0.3535533905932737*(alpha[13]*fUpwind_l[17]+fUpwind_l[13]*alpha[17])+0.3162277660168379*(alpha[6]*fUpwind_l[14]+fUpwind_l[6]*alpha[14]+alpha[4]*fUpwind_l[12]+fUpwind_l[4]*alpha[12])+0.3535533905932737*(alpha[7]*fUpwind_l[11]+fUpwind_l[7]*alpha[11]+alpha[5]*fUpwind_l[10]+fUpwind_l[5]*alpha[10])+0.3162277660168379*(alpha[2]*fUpwind_l[8]+fUpwind_l[2]*alpha[8])+0.3535533905932737*(alpha[3]*fUpwind_l[6]+fUpwind_l[3]*alpha[6]+alpha[1]*fUpwind_l[4]+fUpwind_l[1]*alpha[4]+alpha[0]*fUpwind_l[2]+fUpwind_l[0]*alpha[2]); 
  Ghat_l[3] += 0.3162277660168379*alpha[10]*fUpwind_l[19]+0.3535533905932737*(alpha[12]*fUpwind_l[18]+fUpwind_l[12]*alpha[18]+alpha[11]*fUpwind_l[17]+fUpwind_l[11]*alpha[17])+0.3162277660168379*(alpha[6]*fUpwind_l[16]+alpha[5]*fUpwind_l[15])+0.3535533905932737*(alpha[8]*fUpwind_l[14]+fUpwind_l[8]*alpha[14]+alpha[7]*fUpwind_l[13]+fUpwind_l[7]*alpha[13]+alpha[4]*fUpwind_l[10]+fUpwind_l[4]*alpha[10])+0.3162277660168379*alpha[3]*fUpwind_l[9]+0.3535533905932737*(alpha[2]*fUpwind_l[6]+fUpwind_l[2]*alpha[6]+alpha[1]*fUpwind_l[5]+fUpwind_l[1]*alpha[5]+alpha[0]*fUpwind_l[3]+fUpwind_l[0]*alpha[3]); 
  Ghat_l[4] += (0.2828427124746191*alpha[17]+0.3162277660168379*alpha[6])*fUpwind_l[18]+0.2828427124746191*fUpwind_l[17]*alpha[18]+0.3162277660168379*(fUpwind_l[6]*alpha[18]+alpha[5]*fUpwind_l[17]+fUpwind_l[5]*alpha[17])+0.3162277660168379*(alpha[10]*fUpwind_l[14]+fUpwind_l[10]*alpha[14]+alpha[10]*fUpwind_l[13]+fUpwind_l[10]*alpha[13])+(0.2828427124746191*alpha[11]+0.3162277660168379*alpha[2])*fUpwind_l[12]+0.2828427124746191*fUpwind_l[11]*alpha[12]+0.3162277660168379*(fUpwind_l[2]*alpha[12]+alpha[1]*fUpwind_l[11]+fUpwind_l[1]*alpha[11])+0.3535533905932737*(alpha[3]*fUpwind_l[10]+fUpwind_l[3]*alpha[10])+0.3162277660168379*(alpha[4]*fUpwind_l[8]+fUpwind_l[4]*alpha[8]+alpha[4]*fUpwind_l[7]+fUpwind_l[4]*alpha[7])+0.3535533905932737*(alpha[5]*fUpwind_l[6]+fUpwind_l[5]*alpha[6]+alpha[0]*fUpwind_l[4]+fUpwind_l[0]*alpha[4]+alpha[1]*fUpwind_l[2]+fUpwind_l[1]*alpha[2]); 
  Ghat_l[5] += (0.2828427124746191*alpha[17]+0.3162277660168379*alpha[6])*fUpwind_l[19]+0.3535533905932737*(alpha[8]*fUpwind_l[18]+fUpwind_l[8]*alpha[18])+0.3162277660168379*(alpha[4]*fUpwind_l[17]+fUpwind_l[4]*alpha[17])+0.3162277660168379*alpha[10]*fUpwind_l[16]+(0.2828427124746191*alpha[13]+0.3162277660168379*alpha[3])*fUpwind_l[15]+0.3535533905932737*(alpha[12]*fUpwind_l[14]+fUpwind_l[12]*alpha[14])+0.3162277660168379*(alpha[1]*fUpwind_l[13]+fUpwind_l[1]*alpha[13]+alpha[10]*fUpwind_l[11]+fUpwind_l[10]*alpha[11])+0.3535533905932737*(alpha[2]*fUpwind_l[10]+fUpwind_l[2]*alpha[10])+0.3162277660168379*(alpha[5]*(fUpwind_l[9]+fUpwind_l[7])+fUpwind_l[5]*alpha[7])+0.3535533905932737*(alpha[4]*fUpwind_l[6]+fUpwind_l[4]*alpha[6]+alpha[0]*fUpwind_l[5]+fUpwind_l[0]*alpha[5]+alpha[1]*fUpwind_l[3]+fUpwind_l[1]*alpha[3]); 
  Ghat_l[6] += 0.2828427124746191*alpha[18]*fUpwind_l[19]+0.3162277660168379*(alpha[5]*fUpwind_l[19]+alpha[4]*fUpwind_l[18]+fUpwind_l[4]*alpha[18])+0.3535533905932737*(alpha[7]*fUpwind_l[17]+fUpwind_l[7]*alpha[17])+0.2828427124746191*alpha[14]*fUpwind_l[16]+0.3162277660168379*(alpha[3]*fUpwind_l[16]+alpha[10]*fUpwind_l[15]+alpha[2]*fUpwind_l[14]+fUpwind_l[2]*alpha[14])+0.3535533905932737*(alpha[11]*fUpwind_l[13]+fUpwind_l[11]*alpha[13])+0.3162277660168379*(alpha[10]*fUpwind_l[12]+fUpwind_l[10]*alpha[12])+0.3535533905932737*(alpha[1]*fUpwind_l[10]+fUpwind_l[1]*alpha[10])+0.3162277660168379*(alpha[6]*(fUpwind_l[9]+fUpwind_l[8])+fUpwind_l[6]*alpha[8])+0.3535533905932737*(alpha[0]*fUpwind_l[6]+fUpwind_l[0]*alpha[6]+alpha[4]*fUpwind_l[5]+fUpwind_l[4]*alpha[5]+alpha[2]*fUpwind_l[3]+fUpwind_l[2]*alpha[3]); 
  Ghat_l[7] += 0.3162277660168379*alpha[18]*fUpwind_l[18]+0.2258769757263128*alpha[17]*fUpwind_l[17]+0.3535533905932737*(alpha[6]*fUpwind_l[17]+fUpwind_l[6]*alpha[17])+0.2258769757263128*alpha[13]*fUpwind_l[13]+0.3535533905932737*(alpha[3]*fUpwind_l[13]+fUpwind_l[3]*alpha[13])+0.3162277660168379*alpha[12]*fUpwind_l[12]+0.2258769757263128*alpha[11]*fUpwind_l[11]+0.3535533905932737*(alpha[2]*fUpwind_l[11]+fUpwind_l[2]*alpha[11])+0.3162277660168379*alpha[10]*fUpwind_l[10]+0.2258769757263128*alpha[7]*fUpwind_l[7]+0.3535533905932737*(alpha[0]*fUpwind_l[7]+fUpwind_l[0]*alpha[7])+0.3162277660168379*(alpha[5]*fUpwind_l[5]+alpha[4]*fUpwind_l[4]+alpha[1]*fUpwind_l[1]); 
  Ghat_l[8] += 0.2258769757263128*alpha[18]*fUpwind_l[18]+0.3535533905932737*(alpha[5]*fUpwind_l[18]+fUpwind_l[5]*alpha[18])+0.3162277660168379*alpha[17]*fUpwind_l[17]+0.2258769757263128*alpha[14]*fUpwind_l[14]+0.3535533905932737*(alpha[3]*fUpwind_l[14]+fUpwind_l[3]*alpha[14])+0.2258769757263128*alpha[12]*fUpwind_l[12]+0.3535533905932737*(alpha[1]*fUpwind_l[12]+fUpwind_l[1]*alpha[12])+0.3162277660168379*(alpha[11]*fUpwind_l[11]+alpha[10]*fUpwind_l[10])+0.2258769757263128*alpha[8]*fUpwind_l[8]+0.3535533905932737*(alpha[0]*fUpwind_l[8]+fUpwind_l[0]*alpha[8])+0.3162277660168379*(alpha[6]*fUpwind_l[6]+alpha[4]*fUpwind_l[4]+alpha[2]*fUpwind_l[2]); 
  Ghat_l[9] += 0.3535533905932737*alpha[4]*fUpwind_l[19]+0.3162277660168379*(alpha[18]*fUpwind_l[18]+alpha[17]*fUpwind_l[17])+0.3535533905932737*(alpha[2]*fUpwind_l[16]+alpha[1]*fUpwind_l[15])+0.3162277660168379*(alpha[14]*fUpwind_l[14]+alpha[13]*fUpwind_l[13]+alpha[10]*fUpwind_l[10])+0.3535533905932737*alpha[0]*fUpwind_l[9]+0.3162277660168379*(alpha[6]*fUpwind_l[6]+alpha[5]*fUpwind_l[5]+alpha[3]*fUpwind_l[3]); 
  Ghat_l[10] += (0.282842712474619*(alpha[14]+alpha[13])+0.3162277660168379*alpha[3])*fUpwind_l[19]+(0.282842712474619*alpha[11]+0.3162277660168379*alpha[2])*fUpwind_l[18]+(0.282842712474619*(fUpwind_l[16]+fUpwind_l[11])+0.3162277660168379*fUpwind_l[2])*alpha[18]+(0.282842712474619*alpha[12]+0.3162277660168379*alpha[1])*fUpwind_l[17]+(0.282842712474619*(fUpwind_l[15]+fUpwind_l[12])+0.3162277660168379*fUpwind_l[1])*alpha[17]+0.3162277660168379*(alpha[5]*fUpwind_l[16]+alpha[6]*fUpwind_l[15]+alpha[4]*fUpwind_l[14]+fUpwind_l[4]*alpha[14]+alpha[4]*fUpwind_l[13]+fUpwind_l[4]*alpha[13]+alpha[6]*fUpwind_l[12]+fUpwind_l[6]*alpha[12]+alpha[5]*fUpwind_l[11]+fUpwind_l[5]*alpha[11])+(0.3162277660168379*(alpha[8]+alpha[7])+0.3535533905932737*alpha[0])*fUpwind_l[10]+0.3162277660168379*(fUpwind_l[9]+fUpwind_l[8]+fUpwind_l[7])*alpha[10]+0.3535533905932737*(fUpwind_l[0]*alpha[10]+alpha[1]*fUpwind_l[6]+fUpwind_l[1]*alpha[6]+alpha[2]*fUpwind_l[5]+fUpwind_l[2]*alpha[5]+alpha[3]*fUpwind_l[4]+fUpwind_l[3]*alpha[4]); 
  Ghat_l[11] += 0.282842712474619*(alpha[10]*fUpwind_l[18]+fUpwind_l[10]*alpha[18])+(0.3162277660168379*alpha[14]+0.2258769757263128*alpha[13]+0.3535533905932737*alpha[3])*fUpwind_l[17]+(0.3162277660168379*fUpwind_l[14]+0.2258769757263128*fUpwind_l[13])*alpha[17]+0.3535533905932737*(fUpwind_l[3]*alpha[17]+alpha[6]*fUpwind_l[13]+fUpwind_l[6]*alpha[13])+0.2828427124746191*(alpha[4]*fUpwind_l[12]+fUpwind_l[4]*alpha[12])+(0.3162277660168379*alpha[8]+0.2258769757263128*alpha[7]+0.3535533905932737*alpha[0])*fUpwind_l[11]+(0.3162277660168379*fUpwind_l[8]+0.2258769757263128*fUpwind_l[7]+0.3535533905932737*fUpwind_l[0])*alpha[11]+0.3162277660168379*(alpha[5]*fUpwind_l[10]+fUpwind_l[5]*alpha[10])+0.3535533905932737*(alpha[2]*fUpwind_l[7]+fUpwind_l[2]*alpha[7])+0.3162277660168379*(alpha[1]*fUpwind_l[4]+fUpwind_l[1]*alpha[4]); 
  Ghat_l[12] += (0.2258769757263128*alpha[14]+0.3162277660168379*alpha[13]+0.3535533905932737*alpha[3])*fUpwind_l[18]+(0.2258769757263128*fUpwind_l[14]+0.3162277660168379*fUpwind_l[13]+0.3535533905932737*fUpwind_l[3])*alpha[18]+0.282842712474619*(alpha[10]*fUpwind_l[17]+fUpwind_l[10]*alpha[17])+0.3535533905932737*(alpha[5]*fUpwind_l[14]+fUpwind_l[5]*alpha[14])+(0.2258769757263128*alpha[8]+0.3162277660168379*alpha[7]+0.3535533905932737*alpha[0])*fUpwind_l[12]+(0.2258769757263128*fUpwind_l[8]+0.3162277660168379*fUpwind_l[7]+0.3535533905932737*fUpwind_l[0])*alpha[12]+0.2828427124746191*(alpha[4]*fUpwind_l[11]+fUpwind_l[4]*alpha[11])+0.3162277660168379*(alpha[6]*fUpwind_l[10]+fUpwind_l[6]*alpha[10])+0.3535533905932737*(alpha[1]*fUpwind_l[8]+fUpwind_l[1]*alpha[8])+0.3162277660168379*(alpha[2]*fUpwind_l[4]+fUpwind_l[2]*alpha[4]); 
  Ghat_l[13] += 0.282842712474619*alpha[10]*fUpwind_l[19]+0.3162277660168379*(alpha[12]*fUpwind_l[18]+fUpwind_l[12]*alpha[18])+(0.2258769757263128*alpha[11]+0.3535533905932737*alpha[2])*fUpwind_l[17]+(0.3162277660168379*fUpwind_l[16]+0.2258769757263128*fUpwind_l[11]+0.3535533905932737*fUpwind_l[2])*alpha[17]+0.2828427124746191*alpha[5]*fUpwind_l[15]+(0.2258769757263128*alpha[7]+0.3535533905932737*alpha[0])*fUpwind_l[13]+(0.3162277660168379*fUpwind_l[9]+0.2258769757263128*fUpwind_l[7])*alpha[13]+0.3535533905932737*(fUpwind_l[0]*alpha[13]+alpha[6]*fUpwind_l[11]+fUpwind_l[6]*alpha[11])+0.3162277660168379*(alpha[4]*fUpwind_l[10]+fUpwind_l[4]*alpha[10])+0.3535533905932737*(alpha[3]*fUpwind_l[7]+fUpwind_l[3]*alpha[7])+0.3162277660168379*(alpha[1]*fUpwind_l[5]+fUpwind_l[1]*alpha[5]); 
  Ghat_l[14] += 0.282842712474619*alpha[10]*fUpwind_l[19]+(0.2258769757263128*alpha[12]+0.3535533905932737*alpha[1])*fUpwind_l[18]+(0.3162277660168379*fUpwind_l[15]+0.2258769757263128*fUpwind_l[12]+0.3535533905932737*fUpwind_l[1])*alpha[18]+0.3162277660168379*(alpha[11]*fUpwind_l[17]+fUpwind_l[11]*alpha[17])+0.2828427124746191*alpha[6]*fUpwind_l[16]+(0.2258769757263128*alpha[8]+0.3535533905932737*alpha[0])*fUpwind_l[14]+(0.3162277660168379*fUpwind_l[9]+0.2258769757263128*fUpwind_l[8])*alpha[14]+0.3535533905932737*(fUpwind_l[0]*alpha[14]+alpha[5]*fUpwind_l[12]+fUpwind_l[5]*alpha[12])+0.3162277660168379*(alpha[4]*fUpwind_l[10]+fUpwind_l[4]*alpha[10])+0.3535533905932737*(alpha[3]*fUpwind_l[8]+fUpwind_l[3]*alpha[8])+0.3162277660168379*(alpha[2]*fUpwind_l[6]+fUpwind_l[2]*alpha[6]); 
  Ghat_l[15] += (0.3162277660168379*alpha[11]+0.3535533905932737*alpha[2])*fUpwind_l[19]+0.3162277660168379*(alpha[14]*fUpwind_l[18]+fUpwind_l[14]*alpha[18])+0.282842712474619*(alpha[10]*fUpwind_l[17]+fUpwind_l[10]*alpha[17])+0.3535533905932737*alpha[4]*fUpwind_l[16]+(0.3162277660168379*alpha[7]+0.3535533905932737*alpha[0])*fUpwind_l[15]+0.2828427124746191*(alpha[5]*fUpwind_l[13]+fUpwind_l[5]*alpha[13])+0.3162277660168379*(alpha[6]*fUpwind_l[10]+fUpwind_l[6]*alpha[10])+0.3535533905932737*alpha[1]*fUpwind_l[9]+0.3162277660168379*(alpha[3]*fUpwind_l[5]+fUpwind_l[3]*alpha[5]); 
  Ghat_l[16] += (0.3162277660168379*alpha[12]+0.3535533905932737*alpha[1])*fUpwind_l[19]+0.282842712474619*(alpha[10]*fUpwind_l[18]+fUpwind_l[10]*alpha[18])+0.3162277660168379*(alpha[13]*fUpwind_l[17]+fUpwind_l[13]*alpha[17]+alpha[8]*fUpwind_l[16])+0.3535533905932737*(alpha[0]*fUpwind_l[16]+alpha[4]*fUpwind_l[15])+0.2828427124746191*(alpha[6]*fUpwind_l[14]+fUpwind_l[6]*alpha[14])+0.3162277660168379*(alpha[5]*fUpwind_l[10]+fUpwind_l[5]*alpha[10])+0.3535533905932737*alpha[2]*fUpwind_l[9]+0.3162277660168379*(alpha[3]*fUpwind_l[6]+fUpwind_l[3]*alpha[6]); 
  Ghat_l[17] += 0.2529822128134704*alpha[18]*fUpwind_l[19]+0.2828427124746191*(alpha[5]*fUpwind_l[19]+alpha[4]*fUpwind_l[18]+fUpwind_l[4]*alpha[18])+(0.3162277660168379*alpha[8]+0.2258769757263128*alpha[7]+0.3535533905932737*alpha[0])*fUpwind_l[17]+(0.3162277660168379*(fUpwind_l[9]+fUpwind_l[8])+0.2258769757263128*fUpwind_l[7]+0.3535533905932737*fUpwind_l[0])*alpha[17]+0.3162277660168379*alpha[13]*fUpwind_l[16]+0.282842712474619*alpha[10]*fUpwind_l[15]+0.3162277660168379*(alpha[11]*fUpwind_l[14]+fUpwind_l[11]*alpha[14])+(0.2258769757263128*alpha[11]+0.3535533905932737*alpha[2])*fUpwind_l[13]+(0.2258769757263128*fUpwind_l[11]+0.3535533905932737*fUpwind_l[2])*alpha[13]+0.282842712474619*(alpha[10]*fUpwind_l[12]+fUpwind_l[10]*alpha[12])+0.3535533905932737*(alpha[3]*fUpwind_l[11]+fUpwind_l[3]*alpha[11])+0.3162277660168379*(alpha[1]*fUpwind_l[10]+fUpwind_l[1]*alpha[10])+0.3535533905932737*(alpha[6]*fUpwind_l[7]+fUpwind_l[6]*alpha[7])+0.3162277660168379*(alpha[4]*fUpwind_l[5]+fUpwind_l[4]*alpha[5]); 
  Ghat_l[18] += (0.2529822128134704*alpha[17]+0.2828427124746191*alpha[6])*fUpwind_l[19]+(0.2258769757263128*alpha[8]+0.3162277660168379*alpha[7]+0.3535533905932737*alpha[0])*fUpwind_l[18]+(0.3162277660168379*fUpwind_l[9]+0.2258769757263128*fUpwind_l[8]+0.3162277660168379*fUpwind_l[7]+0.3535533905932737*fUpwind_l[0])*alpha[18]+0.2828427124746191*(alpha[4]*fUpwind_l[17]+fUpwind_l[4]*alpha[17])+0.282842712474619*alpha[10]*fUpwind_l[16]+0.3162277660168379*alpha[14]*fUpwind_l[15]+(0.2258769757263128*alpha[12]+0.3535533905932737*alpha[1])*fUpwind_l[14]+(0.2258769757263128*fUpwind_l[12]+0.3535533905932737*fUpwind_l[1])*alpha[14]+0.3162277660168379*(alpha[12]*fUpwind_l[13]+fUpwind_l[12]*alpha[13])+0.3535533905932737*(alpha[3]*fUpwind_l[12]+fUpwind_l[3]*alpha[12])+0.282842712474619*(alpha[10]*fUpwind_l[11]+fUpwind_l[10]*alpha[11])+0.3162277660168379*(alpha[2]*fUpwind_l[10]+fUpwind_l[2]*alpha[10])+0.3535533905932737*(alpha[5]*fUpwind_l[8]+fUpwind_l[5]*alpha[8])+0.3162277660168379*(alpha[4]*fUpwind_l[6]+fUpwind_l[4]*alpha[6]); 
  Ghat_l[19] += (0.3162277660168379*(alpha[8]+alpha[7])+0.3535533905932737*alpha[0])*fUpwind_l[19]+(0.2529822128134704*alpha[17]+0.2828427124746191*alpha[6])*fUpwind_l[18]+0.2529822128134704*fUpwind_l[17]*alpha[18]+0.2828427124746191*(fUpwind_l[6]*alpha[18]+alpha[5]*fUpwind_l[17]+fUpwind_l[5]*alpha[17])+(0.3162277660168379*alpha[12]+0.3535533905932737*alpha[1])*fUpwind_l[16]+(0.3162277660168379*alpha[11]+0.3535533905932737*alpha[2])*fUpwind_l[15]+0.282842712474619*(alpha[10]*fUpwind_l[14]+fUpwind_l[10]*alpha[14]+alpha[10]*fUpwind_l[13]+fUpwind_l[10]*alpha[13])+0.3162277660168379*(alpha[3]*fUpwind_l[10]+fUpwind_l[3]*alpha[10])+0.3535533905932737*alpha[4]*fUpwind_l[9]+0.3162277660168379*(alpha[5]*fUpwind_l[6]+fUpwind_l[5]*alpha[6]); 

  Ghat_r[0] += 0.3535533905932737*(alpha[18]*fUpwind_r[18]+alpha[17]*fUpwind_r[17]+alpha[14]*fUpwind_r[14]+alpha[13]*fUpwind_r[13]+alpha[12]*fUpwind_r[12]+alpha[11]*fUpwind_r[11]+alpha[10]*fUpwind_r[10]+alpha[8]*fUpwind_r[8]+alpha[7]*fUpwind_r[7]+alpha[6]*fUpwind_r[6]+alpha[5]*fUpwind_r[5]+alpha[4]*fUpwind_r[4]+alpha[3]*fUpwind_r[3]+alpha[2]*fUpwind_r[2]+alpha[1]*fUpwind_r[1]+alpha[0]*fUpwind_r[0]); 
  Ghat_r[1] += 0.3535533905932737*(alpha[14]*fUpwind_r[18]+fUpwind_r[14]*alpha[18])+0.3162277660168379*(alpha[10]*fUpwind_r[17]+fUpwind_r[10]*alpha[17])+0.3162277660168379*(alpha[5]*fUpwind_r[13]+fUpwind_r[5]*alpha[13])+0.3535533905932737*(alpha[8]*fUpwind_r[12]+fUpwind_r[8]*alpha[12])+0.3162277660168379*(alpha[4]*fUpwind_r[11]+fUpwind_r[4]*alpha[11])+0.3535533905932737*(alpha[6]*fUpwind_r[10]+fUpwind_r[6]*alpha[10])+0.3162277660168379*(alpha[1]*fUpwind_r[7]+fUpwind_r[1]*alpha[7])+0.3535533905932737*(alpha[3]*fUpwind_r[5]+fUpwind_r[3]*alpha[5]+alpha[2]*fUpwind_r[4]+fUpwind_r[2]*alpha[4]+alpha[0]*fUpwind_r[1]+fUpwind_r[0]*alpha[1]); 
  Ghat_r[2] += 0.3162277660168379*(alpha[10]*fUpwind_r[18]+fUpwind_r[10]*alpha[18])+0.3535533905932737*(alpha[13]*fUpwind_r[17]+fUpwind_r[13]*alpha[17])+0.3162277660168379*(alpha[6]*fUpwind_r[14]+fUpwind_r[6]*alpha[14]+alpha[4]*fUpwind_r[12]+fUpwind_r[4]*alpha[12])+0.3535533905932737*(alpha[7]*fUpwind_r[11]+fUpwind_r[7]*alpha[11]+alpha[5]*fUpwind_r[10]+fUpwind_r[5]*alpha[10])+0.3162277660168379*(alpha[2]*fUpwind_r[8]+fUpwind_r[2]*alpha[8])+0.3535533905932737*(alpha[3]*fUpwind_r[6]+fUpwind_r[3]*alpha[6]+alpha[1]*fUpwind_r[4]+fUpwind_r[1]*alpha[4]+alpha[0]*fUpwind_r[2]+fUpwind_r[0]*alpha[2]); 
  Ghat_r[3] += 0.3162277660168379*alpha[10]*fUpwind_r[19]+0.3535533905932737*(alpha[12]*fUpwind_r[18]+fUpwind_r[12]*alpha[18]+alpha[11]*fUpwind_r[17]+fUpwind_r[11]*alpha[17])+0.3162277660168379*(alpha[6]*fUpwind_r[16]+alpha[5]*fUpwind_r[15])+0.3535533905932737*(alpha[8]*fUpwind_r[14]+fUpwind_r[8]*alpha[14]+alpha[7]*fUpwind_r[13]+fUpwind_r[7]*alpha[13]+alpha[4]*fUpwind_r[10]+fUpwind_r[4]*alpha[10])+0.3162277660168379*alpha[3]*fUpwind_r[9]+0.3535533905932737*(alpha[2]*fUpwind_r[6]+fUpwind_r[2]*alpha[6]+alpha[1]*fUpwind_r[5]+fUpwind_r[1]*alpha[5]+alpha[0]*fUpwind_r[3]+fUpwind_r[0]*alpha[3]); 
  Ghat_r[4] += (0.2828427124746191*alpha[17]+0.3162277660168379*alpha[6])*fUpwind_r[18]+0.2828427124746191*fUpwind_r[17]*alpha[18]+0.3162277660168379*(fUpwind_r[6]*alpha[18]+alpha[5]*fUpwind_r[17]+fUpwind_r[5]*alpha[17])+0.3162277660168379*(alpha[10]*fUpwind_r[14]+fUpwind_r[10]*alpha[14]+alpha[10]*fUpwind_r[13]+fUpwind_r[10]*alpha[13])+(0.2828427124746191*alpha[11]+0.3162277660168379*alpha[2])*fUpwind_r[12]+0.2828427124746191*fUpwind_r[11]*alpha[12]+0.3162277660168379*(fUpwind_r[2]*alpha[12]+alpha[1]*fUpwind_r[11]+fUpwind_r[1]*alpha[11])+0.3535533905932737*(alpha[3]*fUpwind_r[10]+fUpwind_r[3]*alpha[10])+0.3162277660168379*(alpha[4]*fUpwind_r[8]+fUpwind_r[4]*alpha[8]+alpha[4]*fUpwind_r[7]+fUpwind_r[4]*alpha[7])+0.3535533905932737*(alpha[5]*fUpwind_r[6]+fUpwind_r[5]*alpha[6]+alpha[0]*fUpwind_r[4]+fUpwind_r[0]*alpha[4]+alpha[1]*fUpwind_r[2]+fUpwind_r[1]*alpha[2]); 
  Ghat_r[5] += (0.2828427124746191*alpha[17]+0.3162277660168379*alpha[6])*fUpwind_r[19]+0.3535533905932737*(alpha[8]*fUpwind_r[18]+fUpwind_r[8]*alpha[18])+0.3162277660168379*(alpha[4]*fUpwind_r[17]+fUpwind_r[4]*alpha[17])+0.3162277660168379*alpha[10]*fUpwind_r[16]+(0.2828427124746191*alpha[13]+0.3162277660168379*alpha[3])*fUpwind_r[15]+0.3535533905932737*(alpha[12]*fUpwind_r[14]+fUpwind_r[12]*alpha[14])+0.3162277660168379*(alpha[1]*fUpwind_r[13]+fUpwind_r[1]*alpha[13]+alpha[10]*fUpwind_r[11]+fUpwind_r[10]*alpha[11])+0.3535533905932737*(alpha[2]*fUpwind_r[10]+fUpwind_r[2]*alpha[10])+0.3162277660168379*(alpha[5]*(fUpwind_r[9]+fUpwind_r[7])+fUpwind_r[5]*alpha[7])+0.3535533905932737*(alpha[4]*fUpwind_r[6]+fUpwind_r[4]*alpha[6]+alpha[0]*fUpwind_r[5]+fUpwind_r[0]*alpha[5]+alpha[1]*fUpwind_r[3]+fUpwind_r[1]*alpha[3]); 
  Ghat_r[6] += 0.2828427124746191*alpha[18]*fUpwind_r[19]+0.3162277660168379*(alpha[5]*fUpwind_r[19]+alpha[4]*fUpwind_r[18]+fUpwind_r[4]*alpha[18])+0.3535533905932737*(alpha[7]*fUpwind_r[17]+fUpwind_r[7]*alpha[17])+0.2828427124746191*alpha[14]*fUpwind_r[16]+0.3162277660168379*(alpha[3]*fUpwind_r[16]+alpha[10]*fUpwind_r[15]+alpha[2]*fUpwind_r[14]+fUpwind_r[2]*alpha[14])+0.3535533905932737*(alpha[11]*fUpwind_r[13]+fUpwind_r[11]*alpha[13])+0.3162277660168379*(alpha[10]*fUpwind_r[12]+fUpwind_r[10]*alpha[12])+0.3535533905932737*(alpha[1]*fUpwind_r[10]+fUpwind_r[1]*alpha[10])+0.3162277660168379*(alpha[6]*(fUpwind_r[9]+fUpwind_r[8])+fUpwind_r[6]*alpha[8])+0.3535533905932737*(alpha[0]*fUpwind_r[6]+fUpwind_r[0]*alpha[6]+alpha[4]*fUpwind_r[5]+fUpwind_r[4]*alpha[5]+alpha[2]*fUpwind_r[3]+fUpwind_r[2]*alpha[3]); 
  Ghat_r[7] += 0.3162277660168379*alpha[18]*fUpwind_r[18]+0.2258769757263128*alpha[17]*fUpwind_r[17]+0.3535533905932737*(alpha[6]*fUpwind_r[17]+fUpwind_r[6]*alpha[17])+0.2258769757263128*alpha[13]*fUpwind_r[13]+0.3535533905932737*(alpha[3]*fUpwind_r[13]+fUpwind_r[3]*alpha[13])+0.3162277660168379*alpha[12]*fUpwind_r[12]+0.2258769757263128*alpha[11]*fUpwind_r[11]+0.3535533905932737*(alpha[2]*fUpwind_r[11]+fUpwind_r[2]*alpha[11])+0.3162277660168379*alpha[10]*fUpwind_r[10]+0.2258769757263128*alpha[7]*fUpwind_r[7]+0.3535533905932737*(alpha[0]*fUpwind_r[7]+fUpwind_r[0]*alpha[7])+0.3162277660168379*(alpha[5]*fUpwind_r[5]+alpha[4]*fUpwind_r[4]+alpha[1]*fUpwind_r[1]); 
  Ghat_r[8] += 0.2258769757263128*alpha[18]*fUpwind_r[18]+0.3535533905932737*(alpha[5]*fUpwind_r[18]+fUpwind_r[5]*alpha[18])+0.3162277660168379*alpha[17]*fUpwind_r[17]+0.2258769757263128*alpha[14]*fUpwind_r[14]+0.3535533905932737*(alpha[3]*fUpwind_r[14]+fUpwind_r[3]*alpha[14])+0.2258769757263128*alpha[12]*fUpwind_r[12]+0.3535533905932737*(alpha[1]*fUpwind_r[12]+fUpwind_r[1]*alpha[12])+0.3162277660168379*(alpha[11]*fUpwind_r[11]+alpha[10]*fUpwind_r[10])+0.2258769757263128*alpha[8]*fUpwind_r[8]+0.3535533905932737*(alpha[0]*fUpwind_r[8]+fUpwind_r[0]*alpha[8])+0.3162277660168379*(alpha[6]*fUpwind_r[6]+alpha[4]*fUpwind_r[4]+alpha[2]*fUpwind_r[2]); 
  Ghat_r[9] += 0.3535533905932737*alpha[4]*fUpwind_r[19]+0.3162277660168379*(alpha[18]*fUpwind_r[18]+alpha[17]*fUpwind_r[17])+0.3535533905932737*(alpha[2]*fUpwind_r[16]+alpha[1]*fUpwind_r[15])+0.3162277660168379*(alpha[14]*fUpwind_r[14]+alpha[13]*fUpwind_r[13]+alpha[10]*fUpwind_r[10])+0.3535533905932737*alpha[0]*fUpwind_r[9]+0.3162277660168379*(alpha[6]*fUpwind_r[6]+alpha[5]*fUpwind_r[5]+alpha[3]*fUpwind_r[3]); 
  Ghat_r[10] += (0.282842712474619*(alpha[14]+alpha[13])+0.3162277660168379*alpha[3])*fUpwind_r[19]+(0.282842712474619*alpha[11]+0.3162277660168379*alpha[2])*fUpwind_r[18]+(0.282842712474619*(fUpwind_r[16]+fUpwind_r[11])+0.3162277660168379*fUpwind_r[2])*alpha[18]+(0.282842712474619*alpha[12]+0.3162277660168379*alpha[1])*fUpwind_r[17]+(0.282842712474619*(fUpwind_r[15]+fUpwind_r[12])+0.3162277660168379*fUpwind_r[1])*alpha[17]+0.3162277660168379*(alpha[5]*fUpwind_r[16]+alpha[6]*fUpwind_r[15]+alpha[4]*fUpwind_r[14]+fUpwind_r[4]*alpha[14]+alpha[4]*fUpwind_r[13]+fUpwind_r[4]*alpha[13]+alpha[6]*fUpwind_r[12]+fUpwind_r[6]*alpha[12]+alpha[5]*fUpwind_r[11]+fUpwind_r[5]*alpha[11])+(0.3162277660168379*(alpha[8]+alpha[7])+0.3535533905932737*alpha[0])*fUpwind_r[10]+0.3162277660168379*(fUpwind_r[9]+fUpwind_r[8]+fUpwind_r[7])*alpha[10]+0.3535533905932737*(fUpwind_r[0]*alpha[10]+alpha[1]*fUpwind_r[6]+fUpwind_r[1]*alpha[6]+alpha[2]*fUpwind_r[5]+fUpwind_r[2]*alpha[5]+alpha[3]*fUpwind_r[4]+fUpwind_r[3]*alpha[4]); 
  Ghat_r[11] += 0.282842712474619*(alpha[10]*fUpwind_r[18]+fUpwind_r[10]*alpha[18])+(0.3162277660168379*alpha[14]+0.2258769757263128*alpha[13]+0.3535533905932737*alpha[3])*fUpwind_r[17]+(0.3162277660168379*fUpwind_r[14]+0.2258769757263128*fUpwind_r[13])*alpha[17]+0.3535533905932737*(fUpwind_r[3]*alpha[17]+alpha[6]*fUpwind_r[13]+fUpwind_r[6]*alpha[13])+0.2828427124746191*(alpha[4]*fUpwind_r[12]+fUpwind_r[4]*alpha[12])+(0.3162277660168379*alpha[8]+0.2258769757263128*alpha[7]+0.3535533905932737*alpha[0])*fUpwind_r[11]+(0.3162277660168379*fUpwind_r[8]+0.2258769757263128*fUpwind_r[7]+0.3535533905932737*fUpwind_r[0])*alpha[11]+0.3162277660168379*(alpha[5]*fUpwind_r[10]+fUpwind_r[5]*alpha[10])+0.3535533905932737*(alpha[2]*fUpwind_r[7]+fUpwind_r[2]*alpha[7])+0.3162277660168379*(alpha[1]*fUpwind_r[4]+fUpwind_r[1]*alpha[4]); 
  Ghat_r[12] += (0.2258769757263128*alpha[14]+0.3162277660168379*alpha[13]+0.3535533905932737*alpha[3])*fUpwind_r[18]+(0.2258769757263128*fUpwind_r[14]+0.3162277660168379*fUpwind_r[13]+0.3535533905932737*fUpwind_r[3])*alpha[18]+0.282842712474619*(alpha[10]*fUpwind_r[17]+fUpwind_r[10]*alpha[17])+0.3535533905932737*(alpha[5]*fUpwind_r[14]+fUpwind_r[5]*alpha[14])+(0.2258769757263128*alpha[8]+0.3162277660168379*alpha[7]+0.3535533905932737*alpha[0])*fUpwind_r[12]+(0.2258769757263128*fUpwind_r[8]+0.3162277660168379*fUpwind_r[7]+0.3535533905932737*fUpwind_r[0])*alpha[12]+0.2828427124746191*(alpha[4]*fUpwind_r[11]+fUpwind_r[4]*alpha[11])+0.3162277660168379*(alpha[6]*fUpwind_r[10]+fUpwind_r[6]*alpha[10])+0.3535533905932737*(alpha[1]*fUpwind_r[8]+fUpwind_r[1]*alpha[8])+0.3162277660168379*(alpha[2]*fUpwind_r[4]+fUpwind_r[2]*alpha[4]); 
  Ghat_r[13] += 0.282842712474619*alpha[10]*fUpwind_r[19]+0.3162277660168379*(alpha[12]*fUpwind_r[18]+fUpwind_r[12]*alpha[18])+(0.2258769757263128*alpha[11]+0.3535533905932737*alpha[2])*fUpwind_r[17]+(0.3162277660168379*fUpwind_r[16]+0.2258769757263128*fUpwind_r[11]+0.3535533905932737*fUpwind_r[2])*alpha[17]+0.2828427124746191*alpha[5]*fUpwind_r[15]+(0.2258769757263128*alpha[7]+0.3535533905932737*alpha[0])*fUpwind_r[13]+(0.3162277660168379*fUpwind_r[9]+0.2258769757263128*fUpwind_r[7])*alpha[13]+0.3535533905932737*(fUpwind_r[0]*alpha[13]+alpha[6]*fUpwind_r[11]+fUpwind_r[6]*alpha[11])+0.3162277660168379*(alpha[4]*fUpwind_r[10]+fUpwind_r[4]*alpha[10])+0.3535533905932737*(alpha[3]*fUpwind_r[7]+fUpwind_r[3]*alpha[7])+0.3162277660168379*(alpha[1]*fUpwind_r[5]+fUpwind_r[1]*alpha[5]); 
  Ghat_r[14] += 0.282842712474619*alpha[10]*fUpwind_r[19]+(0.2258769757263128*alpha[12]+0.3535533905932737*alpha[1])*fUpwind_r[18]+(0.3162277660168379*fUpwind_r[15]+0.2258769757263128*fUpwind_r[12]+0.3535533905932737*fUpwind_r[1])*alpha[18]+0.3162277660168379*(alpha[11]*fUpwind_r[17]+fUpwind_r[11]*alpha[17])+0.2828427124746191*alpha[6]*fUpwind_r[16]+(0.2258769757263128*alpha[8]+0.3535533905932737*alpha[0])*fUpwind_r[14]+(0.3162277660168379*fUpwind_r[9]+0.2258769757263128*fUpwind_r[8])*alpha[14]+0.3535533905932737*(fUpwind_r[0]*alpha[14]+alpha[5]*fUpwind_r[12]+fUpwind_r[5]*alpha[12])+0.3162277660168379*(alpha[4]*fUpwind_r[10]+fUpwind_r[4]*alpha[10])+0.3535533905932737*(alpha[3]*fUpwind_r[8]+fUpwind_r[3]*alpha[8])+0.3162277660168379*(alpha[2]*fUpwind_r[6]+fUpwind_r[2]*alpha[6]); 
  Ghat_r[15] += (0.3162277660168379*alpha[11]+0.3535533905932737*alpha[2])*fUpwind_r[19]+0.3162277660168379*(alpha[14]*fUpwind_r[18]+fUpwind_r[14]*alpha[18])+0.282842712474619*(alpha[10]*fUpwind_r[17]+fUpwind_r[10]*alpha[17])+0.3535533905932737*alpha[4]*fUpwind_r[16]+(0.3162277660168379*alpha[7]+0.3535533905932737*alpha[0])*fUpwind_r[15]+0.2828427124746191*(alpha[5]*fUpwind_r[13]+fUpwind_r[5]*alpha[13])+0.3162277660168379*(alpha[6]*fUpwind_r[10]+fUpwind_r[6]*alpha[10])+0.3535533905932737*alpha[1]*fUpwind_r[9]+0.3162277660168379*(alpha[3]*fUpwind_r[5]+fUpwind_r[3]*alpha[5]); 
  Ghat_r[16] += (0.3162277660168379*alpha[12]+0.3535533905932737*alpha[1])*fUpwind_r[19]+0.282842712474619*(alpha[10]*fUpwind_r[18]+fUpwind_r[10]*alpha[18])+0.3162277660168379*(alpha[13]*fUpwind_r[17]+fUpwind_r[13]*alpha[17]+alpha[8]*fUpwind_r[16])+0.3535533905932737*(alpha[0]*fUpwind_r[16]+alpha[4]*fUpwind_r[15])+0.2828427124746191*(alpha[6]*fUpwind_r[14]+fUpwind_r[6]*alpha[14])+0.3162277660168379*(alpha[5]*fUpwind_r[10]+fUpwind_r[5]*alpha[10])+0.3535533905932737*alpha[2]*fUpwind_r[9]+0.3162277660168379*(alpha[3]*fUpwind_r[6]+fUpwind_r[3]*alpha[6]); 
  Ghat_r[17] += 0.2529822128134704*alpha[18]*fUpwind_r[19]+0.2828427124746191*(alpha[5]*fUpwind_r[19]+alpha[4]*fUpwind_r[18]+fUpwind_r[4]*alpha[18])+(0.3162277660168379*alpha[8]+0.2258769757263128*alpha[7]+0.3535533905932737*alpha[0])*fUpwind_r[17]+(0.3162277660168379*(fUpwind_r[9]+fUpwind_r[8])+0.2258769757263128*fUpwind_r[7]+0.3535533905932737*fUpwind_r[0])*alpha[17]+0.3162277660168379*alpha[13]*fUpwind_r[16]+0.282842712474619*alpha[10]*fUpwind_r[15]+0.3162277660168379*(alpha[11]*fUpwind_r[14]+fUpwind_r[11]*alpha[14])+(0.2258769757263128*alpha[11]+0.3535533905932737*alpha[2])*fUpwind_r[13]+(0.2258769757263128*fUpwind_r[11]+0.3535533905932737*fUpwind_r[2])*alpha[13]+0.282842712474619*(alpha[10]*fUpwind_r[12]+fUpwind_r[10]*alpha[12])+0.3535533905932737*(alpha[3]*fUpwind_r[11]+fUpwind_r[3]*alpha[11])+0.3162277660168379*(alpha[1]*fUpwind_r[10]+fUpwind_r[1]*alpha[10])+0.3535533905932737*(alpha[6]*fUpwind_r[7]+fUpwind_r[6]*alpha[7])+0.3162277660168379*(alpha[4]*fUpwind_r[5]+fUpwind_r[4]*alpha[5]); 
  Ghat_r[18] += (0.2529822128134704*alpha[17]+0.2828427124746191*alpha[6])*fUpwind_r[19]+(0.2258769757263128*alpha[8]+0.3162277660168379*alpha[7]+0.3535533905932737*alpha[0])*fUpwind_r[18]+(0.3162277660168379*fUpwind_r[9]+0.2258769757263128*fUpwind_r[8]+0.3162277660168379*fUpwind_r[7]+0.3535533905932737*fUpwind_r[0])*alpha[18]+0.2828427124746191*(alpha[4]*fUpwind_r[17]+fUpwind_r[4]*alpha[17])+0.282842712474619*alpha[10]*fUpwind_r[16]+0.3162277660168379*alpha[14]*fUpwind_r[15]+(0.2258769757263128*alpha[12]+0.3535533905932737*alpha[1])*fUpwind_r[14]+(0.2258769757263128*fUpwind_r[12]+0.3535533905932737*fUpwind_r[1])*alpha[14]+0.3162277660168379*(alpha[12]*fUpwind_r[13]+fUpwind_r[12]*alpha[13])+0.3535533905932737*(alpha[3]*fUpwind_r[12]+fUpwind_r[3]*alpha[12])+0.282842712474619*(alpha[10]*fUpwind_r[11]+fUpwind_r[10]*alpha[11])+0.3162277660168379*(alpha[2]*fUpwind_r[10]+fUpwind_r[2]*alpha[10])+0.3535533905932737*(alpha[5]*fUpwind_r[8]+fUpwind_r[5]*alpha[8])+0.3162277660168379*(alpha[4]*fUpwind_r[6]+fUpwind_r[4]*alpha[6]); 
  Ghat_r[19] += (0.3162277660168379*(alpha[8]+alpha[7])+0.3535533905932737*alpha[0])*fUpwind_r[19]+(0.2529822128134704*alpha[17]+0.2828427124746191*alpha[6])*fUpwind_r[18]+0.2529822128134704*fUpwind_r[17]*alpha[18]+0.2828427124746191*(fUpwind_r[6]*alpha[18]+alpha[5]*fUpwind_r[17]+fUpwind_r[5]*alpha[17])+(0.3162277660168379*alpha[12]+0.3535533905932737*alpha[1])*fUpwind_r[16]+(0.3162277660168379*alpha[11]+0.3535533905932737*alpha[2])*fUpwind_r[15]+0.282842712474619*(alpha[10]*fUpwind_r[14]+fUpwind_r[10]*alpha[14]+alpha[10]*fUpwind_r[13]+fUpwind_r[10]*alpha[13])+0.3162277660168379*(alpha[3]*fUpwind_r[10]+fUpwind_r[3]*alpha[10])+0.3535533905932737*alpha[4]*fUpwind_r[9]+0.3162277660168379*(alpha[5]*fUpwind_r[6]+fUpwind_r[5]*alpha[6]); 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv11; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv11; 
  out[2] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv11; 
  out[3] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv11; 
  out[4] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv11; 
  out[5] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dv11; 
  out[6] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dv11; 
  out[7] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dv11; 
  out[8] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv11; 
  out[9] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv11; 
  out[10] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv11; 
  out[11] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dv11; 
  out[12] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*dv11; 
  out[13] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*dv11; 
  out[14] += (1.58113883008419*Ghat_l[0]-1.58113883008419*Ghat_r[0])*dv11; 
  out[15] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*dv11; 
  out[16] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dv11; 
  out[17] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dv11; 
  out[18] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dv11; 
  out[19] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*dv11; 
  out[20] += (0.7071067811865475*Ghat_l[12]-0.7071067811865475*Ghat_r[12])*dv11; 
  out[21] += (0.7071067811865475*Ghat_l[13]-0.7071067811865475*Ghat_r[13])*dv11; 
  out[22] += (0.7071067811865475*Ghat_l[14]-0.7071067811865475*Ghat_r[14])*dv11; 
  out[23] += (0.7071067811865475*Ghat_l[15]-0.7071067811865475*Ghat_r[15])*dv11; 
  out[24] += (0.7071067811865475*Ghat_l[16]-0.7071067811865475*Ghat_r[16])*dv11; 
  out[25] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dv11; 
  out[26] += -1.224744871391589*(Ghat_r[8]+Ghat_l[8])*dv11; 
  out[27] += -1.224744871391589*(Ghat_r[9]+Ghat_l[9])*dv11; 
  out[28] += (1.58113883008419*Ghat_l[1]-1.58113883008419*Ghat_r[1])*dv11; 
  out[29] += (1.58113883008419*Ghat_l[2]-1.58113883008419*Ghat_r[2])*dv11; 
  out[30] += (1.58113883008419*Ghat_l[3]-1.58113883008419*Ghat_r[3])*dv11; 
  out[31] += -1.224744871391589*(Ghat_r[10]+Ghat_l[10])*dv11; 
  out[32] += (0.7071067811865475*Ghat_l[17]-0.7071067811865475*Ghat_r[17])*dv11; 
  out[33] += (0.7071067811865475*Ghat_l[18]-0.7071067811865475*Ghat_r[18])*dv11; 
  out[34] += (0.7071067811865475*Ghat_l[19]-0.7071067811865475*Ghat_r[19])*dv11; 
  out[35] += -1.224744871391589*(Ghat_r[11]+Ghat_l[11])*dv11; 
  out[36] += -1.224744871391589*(Ghat_r[12]+Ghat_l[12])*dv11; 
  out[37] += -1.224744871391589*(Ghat_r[13]+Ghat_l[13])*dv11; 
  out[38] += -1.224744871391589*(Ghat_r[14]+Ghat_l[14])*dv11; 
  out[39] += -1.224744871391589*(Ghat_r[15]+Ghat_l[15])*dv11; 
  out[40] += -1.224744871391589*(Ghat_r[16]+Ghat_l[16])*dv11; 
  out[41] += (1.58113883008419*Ghat_l[4]-1.58113883008419*Ghat_r[4])*dv11; 
  out[42] += (1.58113883008419*Ghat_l[5]-1.58113883008419*Ghat_r[5])*dv11; 
  out[43] += (1.58113883008419*Ghat_l[6]-1.58113883008419*Ghat_r[6])*dv11; 
  out[44] += -1.224744871391589*(Ghat_r[17]+Ghat_l[17])*dv11; 
  out[45] += -1.224744871391589*(Ghat_r[18]+Ghat_l[18])*dv11; 
  out[46] += -1.224744871391589*(Ghat_r[19]+Ghat_l[19])*dv11; 
  out[47] += (1.58113883008419*Ghat_l[10]-1.58113883008419*Ghat_r[10])*dv11; 

} 
