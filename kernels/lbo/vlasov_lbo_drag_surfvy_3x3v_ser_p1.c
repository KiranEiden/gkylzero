#include <gkyl_vlasov_lbo_kernels.h> 
#include <gkyl_basis_ser_3x3v_p1_surfvy_quad.h> 
GKYL_CU_DH void vlasov_lbo_drag_surfvy_3x3v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[6]:         cell-center coordinates. 
  // dxv[6]:       cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[24]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[8]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      distribution function in cells 
  // out:           incremented distribution function in cell 
  double rdv2 = 2.0/dxv[4]; 
  double rdvSq4 = 4.0/(dxv[4]*dxv[4]); 

  const double *sumNuUy = &nuUSum[8]; 

  double alphaDrSurf_l[32] = {0.0}; 
  alphaDrSurf_l[0] = 2.0*nuSum[0]*w[4]-1.0*nuSum[0]*dxv[4]-2.0*sumNuUy[0]; 
  alphaDrSurf_l[1] = 2.0*nuSum[1]*w[4]-1.0*nuSum[1]*dxv[4]-2.0*sumNuUy[1]; 
  alphaDrSurf_l[2] = 2.0*nuSum[2]*w[4]-1.0*nuSum[2]*dxv[4]-2.0*sumNuUy[2]; 
  alphaDrSurf_l[3] = 2.0*nuSum[3]*w[4]-1.0*nuSum[3]*dxv[4]-2.0*sumNuUy[3]; 
  alphaDrSurf_l[6] = 2.0*nuSum[4]*w[4]-2.0*sumNuUy[4]-1.0*dxv[4]*nuSum[4]; 
  alphaDrSurf_l[7] = (-2.0*sumNuUy[5])+2.0*w[4]*nuSum[5]-1.0*dxv[4]*nuSum[5]; 
  alphaDrSurf_l[8] = (-2.0*sumNuUy[6])+2.0*w[4]*nuSum[6]-1.0*dxv[4]*nuSum[6]; 
  alphaDrSurf_l[16] = (-2.0*sumNuUy[7])+2.0*w[4]*nuSum[7]-1.0*dxv[4]*nuSum[7]; 

  double alphaDrSurf_r[32] = {0.0}; 
  alphaDrSurf_r[0] = 2.0*nuSum[0]*w[4]+nuSum[0]*dxv[4]-2.0*sumNuUy[0]; 
  alphaDrSurf_r[1] = 2.0*nuSum[1]*w[4]+nuSum[1]*dxv[4]-2.0*sumNuUy[1]; 
  alphaDrSurf_r[2] = 2.0*nuSum[2]*w[4]+nuSum[2]*dxv[4]-2.0*sumNuUy[2]; 
  alphaDrSurf_r[3] = 2.0*nuSum[3]*w[4]+nuSum[3]*dxv[4]-2.0*sumNuUy[3]; 
  alphaDrSurf_r[6] = 2.0*nuSum[4]*w[4]-2.0*sumNuUy[4]+dxv[4]*nuSum[4]; 
  alphaDrSurf_r[7] = (-2.0*sumNuUy[5])+2.0*w[4]*nuSum[5]+dxv[4]*nuSum[5]; 
  alphaDrSurf_r[8] = (-2.0*sumNuUy[6])+2.0*w[4]*nuSum[6]+dxv[4]*nuSum[6]; 
  alphaDrSurf_r[16] = (-2.0*sumNuUy[7])+2.0*w[4]*nuSum[7]+dxv[4]*nuSum[7]; 

  double fUpwindQuad_l[32] = {0.0};
  double fUpwindQuad_r[32] = {0.0};
  double fUpwind_l[32] = {0.0};;
  double fUpwind_r[32] = {0.0};
  double Gdrag_l[32] = {0.0}; 
  double Gdrag_r[32] = {0.0}; 

  if ((-alphaDrSurf_l[16])+alphaDrSurf_l[8]+alphaDrSurf_l[7]+alphaDrSurf_l[6]-alphaDrSurf_l[3]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[0] = ser_3x3v_p1_surfvy_quad_0(1, fl); 
  } else { 
    fUpwindQuad_l[0] = ser_3x3v_p1_surfvy_quad_0(-1, fc); 
  } 
  if (alphaDrSurf_l[16]+alphaDrSurf_l[8]-alphaDrSurf_l[7]-alphaDrSurf_l[6]-alphaDrSurf_l[3]-alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[1] = ser_3x3v_p1_surfvy_quad_1(1, fl); 
  } else { 
    fUpwindQuad_l[1] = ser_3x3v_p1_surfvy_quad_1(-1, fc); 
  } 
  if (alphaDrSurf_l[16]-alphaDrSurf_l[8]+alphaDrSurf_l[7]-alphaDrSurf_l[6]-alphaDrSurf_l[3]+alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[2] = ser_3x3v_p1_surfvy_quad_2(1, fl); 
  } else { 
    fUpwindQuad_l[2] = ser_3x3v_p1_surfvy_quad_2(-1, fc); 
  } 
  if ((-alphaDrSurf_l[16])-alphaDrSurf_l[8]-alphaDrSurf_l[7]+alphaDrSurf_l[6]-alphaDrSurf_l[3]+alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[3] = ser_3x3v_p1_surfvy_quad_3(1, fl); 
  } else { 
    fUpwindQuad_l[3] = ser_3x3v_p1_surfvy_quad_3(-1, fc); 
  } 
  if (alphaDrSurf_l[16]-alphaDrSurf_l[8]-alphaDrSurf_l[7]+alphaDrSurf_l[6]+alphaDrSurf_l[3]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[4] = ser_3x3v_p1_surfvy_quad_4(1, fl); 
  } else { 
    fUpwindQuad_l[4] = ser_3x3v_p1_surfvy_quad_4(-1, fc); 
  } 
  if ((-alphaDrSurf_l[16])-alphaDrSurf_l[8]+alphaDrSurf_l[7]-alphaDrSurf_l[6]+alphaDrSurf_l[3]-alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[5] = ser_3x3v_p1_surfvy_quad_5(1, fl); 
  } else { 
    fUpwindQuad_l[5] = ser_3x3v_p1_surfvy_quad_5(-1, fc); 
  } 
  if ((-alphaDrSurf_l[16])+alphaDrSurf_l[8]-alphaDrSurf_l[7]-alphaDrSurf_l[6]+alphaDrSurf_l[3]+alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[6] = ser_3x3v_p1_surfvy_quad_6(1, fl); 
  } else { 
    fUpwindQuad_l[6] = ser_3x3v_p1_surfvy_quad_6(-1, fc); 
  } 
  if (alphaDrSurf_l[16]+alphaDrSurf_l[8]+alphaDrSurf_l[7]+alphaDrSurf_l[6]+alphaDrSurf_l[3]+alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[7] = ser_3x3v_p1_surfvy_quad_7(1, fl); 
  } else { 
    fUpwindQuad_l[7] = ser_3x3v_p1_surfvy_quad_7(-1, fc); 
  } 
  if ((-alphaDrSurf_l[16])+alphaDrSurf_l[8]+alphaDrSurf_l[7]+alphaDrSurf_l[6]-alphaDrSurf_l[3]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[8] = ser_3x3v_p1_surfvy_quad_8(1, fl); 
  } else { 
    fUpwindQuad_l[8] = ser_3x3v_p1_surfvy_quad_8(-1, fc); 
  } 
  if (alphaDrSurf_l[16]+alphaDrSurf_l[8]-alphaDrSurf_l[7]-alphaDrSurf_l[6]-alphaDrSurf_l[3]-alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[9] = ser_3x3v_p1_surfvy_quad_9(1, fl); 
  } else { 
    fUpwindQuad_l[9] = ser_3x3v_p1_surfvy_quad_9(-1, fc); 
  } 
  if (alphaDrSurf_l[16]-alphaDrSurf_l[8]+alphaDrSurf_l[7]-alphaDrSurf_l[6]-alphaDrSurf_l[3]+alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[10] = ser_3x3v_p1_surfvy_quad_10(1, fl); 
  } else { 
    fUpwindQuad_l[10] = ser_3x3v_p1_surfvy_quad_10(-1, fc); 
  } 
  if ((-alphaDrSurf_l[16])-alphaDrSurf_l[8]-alphaDrSurf_l[7]+alphaDrSurf_l[6]-alphaDrSurf_l[3]+alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[11] = ser_3x3v_p1_surfvy_quad_11(1, fl); 
  } else { 
    fUpwindQuad_l[11] = ser_3x3v_p1_surfvy_quad_11(-1, fc); 
  } 
  if (alphaDrSurf_l[16]-alphaDrSurf_l[8]-alphaDrSurf_l[7]+alphaDrSurf_l[6]+alphaDrSurf_l[3]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[12] = ser_3x3v_p1_surfvy_quad_12(1, fl); 
  } else { 
    fUpwindQuad_l[12] = ser_3x3v_p1_surfvy_quad_12(-1, fc); 
  } 
  if ((-alphaDrSurf_l[16])-alphaDrSurf_l[8]+alphaDrSurf_l[7]-alphaDrSurf_l[6]+alphaDrSurf_l[3]-alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[13] = ser_3x3v_p1_surfvy_quad_13(1, fl); 
  } else { 
    fUpwindQuad_l[13] = ser_3x3v_p1_surfvy_quad_13(-1, fc); 
  } 
  if ((-alphaDrSurf_l[16])+alphaDrSurf_l[8]-alphaDrSurf_l[7]-alphaDrSurf_l[6]+alphaDrSurf_l[3]+alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[14] = ser_3x3v_p1_surfvy_quad_14(1, fl); 
  } else { 
    fUpwindQuad_l[14] = ser_3x3v_p1_surfvy_quad_14(-1, fc); 
  } 
  if (alphaDrSurf_l[16]+alphaDrSurf_l[8]+alphaDrSurf_l[7]+alphaDrSurf_l[6]+alphaDrSurf_l[3]+alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[15] = ser_3x3v_p1_surfvy_quad_15(1, fl); 
  } else { 
    fUpwindQuad_l[15] = ser_3x3v_p1_surfvy_quad_15(-1, fc); 
  } 
  if ((-alphaDrSurf_l[16])+alphaDrSurf_l[8]+alphaDrSurf_l[7]+alphaDrSurf_l[6]-alphaDrSurf_l[3]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[16] = ser_3x3v_p1_surfvy_quad_16(1, fl); 
  } else { 
    fUpwindQuad_l[16] = ser_3x3v_p1_surfvy_quad_16(-1, fc); 
  } 
  if (alphaDrSurf_l[16]+alphaDrSurf_l[8]-alphaDrSurf_l[7]-alphaDrSurf_l[6]-alphaDrSurf_l[3]-alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[17] = ser_3x3v_p1_surfvy_quad_17(1, fl); 
  } else { 
    fUpwindQuad_l[17] = ser_3x3v_p1_surfvy_quad_17(-1, fc); 
  } 
  if (alphaDrSurf_l[16]-alphaDrSurf_l[8]+alphaDrSurf_l[7]-alphaDrSurf_l[6]-alphaDrSurf_l[3]+alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[18] = ser_3x3v_p1_surfvy_quad_18(1, fl); 
  } else { 
    fUpwindQuad_l[18] = ser_3x3v_p1_surfvy_quad_18(-1, fc); 
  } 
  if ((-alphaDrSurf_l[16])-alphaDrSurf_l[8]-alphaDrSurf_l[7]+alphaDrSurf_l[6]-alphaDrSurf_l[3]+alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[19] = ser_3x3v_p1_surfvy_quad_19(1, fl); 
  } else { 
    fUpwindQuad_l[19] = ser_3x3v_p1_surfvy_quad_19(-1, fc); 
  } 
  if (alphaDrSurf_l[16]-alphaDrSurf_l[8]-alphaDrSurf_l[7]+alphaDrSurf_l[6]+alphaDrSurf_l[3]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[20] = ser_3x3v_p1_surfvy_quad_20(1, fl); 
  } else { 
    fUpwindQuad_l[20] = ser_3x3v_p1_surfvy_quad_20(-1, fc); 
  } 
  if ((-alphaDrSurf_l[16])-alphaDrSurf_l[8]+alphaDrSurf_l[7]-alphaDrSurf_l[6]+alphaDrSurf_l[3]-alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[21] = ser_3x3v_p1_surfvy_quad_21(1, fl); 
  } else { 
    fUpwindQuad_l[21] = ser_3x3v_p1_surfvy_quad_21(-1, fc); 
  } 
  if ((-alphaDrSurf_l[16])+alphaDrSurf_l[8]-alphaDrSurf_l[7]-alphaDrSurf_l[6]+alphaDrSurf_l[3]+alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[22] = ser_3x3v_p1_surfvy_quad_22(1, fl); 
  } else { 
    fUpwindQuad_l[22] = ser_3x3v_p1_surfvy_quad_22(-1, fc); 
  } 
  if (alphaDrSurf_l[16]+alphaDrSurf_l[8]+alphaDrSurf_l[7]+alphaDrSurf_l[6]+alphaDrSurf_l[3]+alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[23] = ser_3x3v_p1_surfvy_quad_23(1, fl); 
  } else { 
    fUpwindQuad_l[23] = ser_3x3v_p1_surfvy_quad_23(-1, fc); 
  } 
  if ((-alphaDrSurf_l[16])+alphaDrSurf_l[8]+alphaDrSurf_l[7]+alphaDrSurf_l[6]-alphaDrSurf_l[3]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[24] = ser_3x3v_p1_surfvy_quad_24(1, fl); 
  } else { 
    fUpwindQuad_l[24] = ser_3x3v_p1_surfvy_quad_24(-1, fc); 
  } 
  if (alphaDrSurf_l[16]+alphaDrSurf_l[8]-alphaDrSurf_l[7]-alphaDrSurf_l[6]-alphaDrSurf_l[3]-alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[25] = ser_3x3v_p1_surfvy_quad_25(1, fl); 
  } else { 
    fUpwindQuad_l[25] = ser_3x3v_p1_surfvy_quad_25(-1, fc); 
  } 
  if (alphaDrSurf_l[16]-alphaDrSurf_l[8]+alphaDrSurf_l[7]-alphaDrSurf_l[6]-alphaDrSurf_l[3]+alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[26] = ser_3x3v_p1_surfvy_quad_26(1, fl); 
  } else { 
    fUpwindQuad_l[26] = ser_3x3v_p1_surfvy_quad_26(-1, fc); 
  } 
  if ((-alphaDrSurf_l[16])-alphaDrSurf_l[8]-alphaDrSurf_l[7]+alphaDrSurf_l[6]-alphaDrSurf_l[3]+alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[27] = ser_3x3v_p1_surfvy_quad_27(1, fl); 
  } else { 
    fUpwindQuad_l[27] = ser_3x3v_p1_surfvy_quad_27(-1, fc); 
  } 
  if (alphaDrSurf_l[16]-alphaDrSurf_l[8]-alphaDrSurf_l[7]+alphaDrSurf_l[6]+alphaDrSurf_l[3]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[28] = ser_3x3v_p1_surfvy_quad_28(1, fl); 
  } else { 
    fUpwindQuad_l[28] = ser_3x3v_p1_surfvy_quad_28(-1, fc); 
  } 
  if ((-alphaDrSurf_l[16])-alphaDrSurf_l[8]+alphaDrSurf_l[7]-alphaDrSurf_l[6]+alphaDrSurf_l[3]-alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[29] = ser_3x3v_p1_surfvy_quad_29(1, fl); 
  } else { 
    fUpwindQuad_l[29] = ser_3x3v_p1_surfvy_quad_29(-1, fc); 
  } 
  if ((-alphaDrSurf_l[16])+alphaDrSurf_l[8]-alphaDrSurf_l[7]-alphaDrSurf_l[6]+alphaDrSurf_l[3]+alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[30] = ser_3x3v_p1_surfvy_quad_30(1, fl); 
  } else { 
    fUpwindQuad_l[30] = ser_3x3v_p1_surfvy_quad_30(-1, fc); 
  } 
  if (alphaDrSurf_l[16]+alphaDrSurf_l[8]+alphaDrSurf_l[7]+alphaDrSurf_l[6]+alphaDrSurf_l[3]+alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[31] = ser_3x3v_p1_surfvy_quad_31(1, fl); 
  } else { 
    fUpwindQuad_l[31] = ser_3x3v_p1_surfvy_quad_31(-1, fc); 
  } 
  if ((-alphaDrSurf_r[16])+alphaDrSurf_r[8]+alphaDrSurf_r[7]+alphaDrSurf_r[6]-alphaDrSurf_r[3]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[0] = ser_3x3v_p1_surfvy_quad_0(1, fc); 
  } else { 
    fUpwindQuad_r[0] = ser_3x3v_p1_surfvy_quad_0(-1, fr); 
  } 
  if (alphaDrSurf_r[16]+alphaDrSurf_r[8]-alphaDrSurf_r[7]-alphaDrSurf_r[6]-alphaDrSurf_r[3]-alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[1] = ser_3x3v_p1_surfvy_quad_1(1, fc); 
  } else { 
    fUpwindQuad_r[1] = ser_3x3v_p1_surfvy_quad_1(-1, fr); 
  } 
  if (alphaDrSurf_r[16]-alphaDrSurf_r[8]+alphaDrSurf_r[7]-alphaDrSurf_r[6]-alphaDrSurf_r[3]+alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[2] = ser_3x3v_p1_surfvy_quad_2(1, fc); 
  } else { 
    fUpwindQuad_r[2] = ser_3x3v_p1_surfvy_quad_2(-1, fr); 
  } 
  if ((-alphaDrSurf_r[16])-alphaDrSurf_r[8]-alphaDrSurf_r[7]+alphaDrSurf_r[6]-alphaDrSurf_r[3]+alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[3] = ser_3x3v_p1_surfvy_quad_3(1, fc); 
  } else { 
    fUpwindQuad_r[3] = ser_3x3v_p1_surfvy_quad_3(-1, fr); 
  } 
  if (alphaDrSurf_r[16]-alphaDrSurf_r[8]-alphaDrSurf_r[7]+alphaDrSurf_r[6]+alphaDrSurf_r[3]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[4] = ser_3x3v_p1_surfvy_quad_4(1, fc); 
  } else { 
    fUpwindQuad_r[4] = ser_3x3v_p1_surfvy_quad_4(-1, fr); 
  } 
  if ((-alphaDrSurf_r[16])-alphaDrSurf_r[8]+alphaDrSurf_r[7]-alphaDrSurf_r[6]+alphaDrSurf_r[3]-alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[5] = ser_3x3v_p1_surfvy_quad_5(1, fc); 
  } else { 
    fUpwindQuad_r[5] = ser_3x3v_p1_surfvy_quad_5(-1, fr); 
  } 
  if ((-alphaDrSurf_r[16])+alphaDrSurf_r[8]-alphaDrSurf_r[7]-alphaDrSurf_r[6]+alphaDrSurf_r[3]+alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[6] = ser_3x3v_p1_surfvy_quad_6(1, fc); 
  } else { 
    fUpwindQuad_r[6] = ser_3x3v_p1_surfvy_quad_6(-1, fr); 
  } 
  if (alphaDrSurf_r[16]+alphaDrSurf_r[8]+alphaDrSurf_r[7]+alphaDrSurf_r[6]+alphaDrSurf_r[3]+alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[7] = ser_3x3v_p1_surfvy_quad_7(1, fc); 
  } else { 
    fUpwindQuad_r[7] = ser_3x3v_p1_surfvy_quad_7(-1, fr); 
  } 
  if ((-alphaDrSurf_r[16])+alphaDrSurf_r[8]+alphaDrSurf_r[7]+alphaDrSurf_r[6]-alphaDrSurf_r[3]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[8] = ser_3x3v_p1_surfvy_quad_8(1, fc); 
  } else { 
    fUpwindQuad_r[8] = ser_3x3v_p1_surfvy_quad_8(-1, fr); 
  } 
  if (alphaDrSurf_r[16]+alphaDrSurf_r[8]-alphaDrSurf_r[7]-alphaDrSurf_r[6]-alphaDrSurf_r[3]-alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[9] = ser_3x3v_p1_surfvy_quad_9(1, fc); 
  } else { 
    fUpwindQuad_r[9] = ser_3x3v_p1_surfvy_quad_9(-1, fr); 
  } 
  if (alphaDrSurf_r[16]-alphaDrSurf_r[8]+alphaDrSurf_r[7]-alphaDrSurf_r[6]-alphaDrSurf_r[3]+alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[10] = ser_3x3v_p1_surfvy_quad_10(1, fc); 
  } else { 
    fUpwindQuad_r[10] = ser_3x3v_p1_surfvy_quad_10(-1, fr); 
  } 
  if ((-alphaDrSurf_r[16])-alphaDrSurf_r[8]-alphaDrSurf_r[7]+alphaDrSurf_r[6]-alphaDrSurf_r[3]+alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[11] = ser_3x3v_p1_surfvy_quad_11(1, fc); 
  } else { 
    fUpwindQuad_r[11] = ser_3x3v_p1_surfvy_quad_11(-1, fr); 
  } 
  if (alphaDrSurf_r[16]-alphaDrSurf_r[8]-alphaDrSurf_r[7]+alphaDrSurf_r[6]+alphaDrSurf_r[3]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[12] = ser_3x3v_p1_surfvy_quad_12(1, fc); 
  } else { 
    fUpwindQuad_r[12] = ser_3x3v_p1_surfvy_quad_12(-1, fr); 
  } 
  if ((-alphaDrSurf_r[16])-alphaDrSurf_r[8]+alphaDrSurf_r[7]-alphaDrSurf_r[6]+alphaDrSurf_r[3]-alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[13] = ser_3x3v_p1_surfvy_quad_13(1, fc); 
  } else { 
    fUpwindQuad_r[13] = ser_3x3v_p1_surfvy_quad_13(-1, fr); 
  } 
  if ((-alphaDrSurf_r[16])+alphaDrSurf_r[8]-alphaDrSurf_r[7]-alphaDrSurf_r[6]+alphaDrSurf_r[3]+alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[14] = ser_3x3v_p1_surfvy_quad_14(1, fc); 
  } else { 
    fUpwindQuad_r[14] = ser_3x3v_p1_surfvy_quad_14(-1, fr); 
  } 
  if (alphaDrSurf_r[16]+alphaDrSurf_r[8]+alphaDrSurf_r[7]+alphaDrSurf_r[6]+alphaDrSurf_r[3]+alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[15] = ser_3x3v_p1_surfvy_quad_15(1, fc); 
  } else { 
    fUpwindQuad_r[15] = ser_3x3v_p1_surfvy_quad_15(-1, fr); 
  } 
  if ((-alphaDrSurf_r[16])+alphaDrSurf_r[8]+alphaDrSurf_r[7]+alphaDrSurf_r[6]-alphaDrSurf_r[3]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[16] = ser_3x3v_p1_surfvy_quad_16(1, fc); 
  } else { 
    fUpwindQuad_r[16] = ser_3x3v_p1_surfvy_quad_16(-1, fr); 
  } 
  if (alphaDrSurf_r[16]+alphaDrSurf_r[8]-alphaDrSurf_r[7]-alphaDrSurf_r[6]-alphaDrSurf_r[3]-alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[17] = ser_3x3v_p1_surfvy_quad_17(1, fc); 
  } else { 
    fUpwindQuad_r[17] = ser_3x3v_p1_surfvy_quad_17(-1, fr); 
  } 
  if (alphaDrSurf_r[16]-alphaDrSurf_r[8]+alphaDrSurf_r[7]-alphaDrSurf_r[6]-alphaDrSurf_r[3]+alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[18] = ser_3x3v_p1_surfvy_quad_18(1, fc); 
  } else { 
    fUpwindQuad_r[18] = ser_3x3v_p1_surfvy_quad_18(-1, fr); 
  } 
  if ((-alphaDrSurf_r[16])-alphaDrSurf_r[8]-alphaDrSurf_r[7]+alphaDrSurf_r[6]-alphaDrSurf_r[3]+alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[19] = ser_3x3v_p1_surfvy_quad_19(1, fc); 
  } else { 
    fUpwindQuad_r[19] = ser_3x3v_p1_surfvy_quad_19(-1, fr); 
  } 
  if (alphaDrSurf_r[16]-alphaDrSurf_r[8]-alphaDrSurf_r[7]+alphaDrSurf_r[6]+alphaDrSurf_r[3]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[20] = ser_3x3v_p1_surfvy_quad_20(1, fc); 
  } else { 
    fUpwindQuad_r[20] = ser_3x3v_p1_surfvy_quad_20(-1, fr); 
  } 
  if ((-alphaDrSurf_r[16])-alphaDrSurf_r[8]+alphaDrSurf_r[7]-alphaDrSurf_r[6]+alphaDrSurf_r[3]-alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[21] = ser_3x3v_p1_surfvy_quad_21(1, fc); 
  } else { 
    fUpwindQuad_r[21] = ser_3x3v_p1_surfvy_quad_21(-1, fr); 
  } 
  if ((-alphaDrSurf_r[16])+alphaDrSurf_r[8]-alphaDrSurf_r[7]-alphaDrSurf_r[6]+alphaDrSurf_r[3]+alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[22] = ser_3x3v_p1_surfvy_quad_22(1, fc); 
  } else { 
    fUpwindQuad_r[22] = ser_3x3v_p1_surfvy_quad_22(-1, fr); 
  } 
  if (alphaDrSurf_r[16]+alphaDrSurf_r[8]+alphaDrSurf_r[7]+alphaDrSurf_r[6]+alphaDrSurf_r[3]+alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[23] = ser_3x3v_p1_surfvy_quad_23(1, fc); 
  } else { 
    fUpwindQuad_r[23] = ser_3x3v_p1_surfvy_quad_23(-1, fr); 
  } 
  if ((-alphaDrSurf_r[16])+alphaDrSurf_r[8]+alphaDrSurf_r[7]+alphaDrSurf_r[6]-alphaDrSurf_r[3]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[24] = ser_3x3v_p1_surfvy_quad_24(1, fc); 
  } else { 
    fUpwindQuad_r[24] = ser_3x3v_p1_surfvy_quad_24(-1, fr); 
  } 
  if (alphaDrSurf_r[16]+alphaDrSurf_r[8]-alphaDrSurf_r[7]-alphaDrSurf_r[6]-alphaDrSurf_r[3]-alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[25] = ser_3x3v_p1_surfvy_quad_25(1, fc); 
  } else { 
    fUpwindQuad_r[25] = ser_3x3v_p1_surfvy_quad_25(-1, fr); 
  } 
  if (alphaDrSurf_r[16]-alphaDrSurf_r[8]+alphaDrSurf_r[7]-alphaDrSurf_r[6]-alphaDrSurf_r[3]+alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[26] = ser_3x3v_p1_surfvy_quad_26(1, fc); 
  } else { 
    fUpwindQuad_r[26] = ser_3x3v_p1_surfvy_quad_26(-1, fr); 
  } 
  if ((-alphaDrSurf_r[16])-alphaDrSurf_r[8]-alphaDrSurf_r[7]+alphaDrSurf_r[6]-alphaDrSurf_r[3]+alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[27] = ser_3x3v_p1_surfvy_quad_27(1, fc); 
  } else { 
    fUpwindQuad_r[27] = ser_3x3v_p1_surfvy_quad_27(-1, fr); 
  } 
  if (alphaDrSurf_r[16]-alphaDrSurf_r[8]-alphaDrSurf_r[7]+alphaDrSurf_r[6]+alphaDrSurf_r[3]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[28] = ser_3x3v_p1_surfvy_quad_28(1, fc); 
  } else { 
    fUpwindQuad_r[28] = ser_3x3v_p1_surfvy_quad_28(-1, fr); 
  } 
  if ((-alphaDrSurf_r[16])-alphaDrSurf_r[8]+alphaDrSurf_r[7]-alphaDrSurf_r[6]+alphaDrSurf_r[3]-alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[29] = ser_3x3v_p1_surfvy_quad_29(1, fc); 
  } else { 
    fUpwindQuad_r[29] = ser_3x3v_p1_surfvy_quad_29(-1, fr); 
  } 
  if ((-alphaDrSurf_r[16])+alphaDrSurf_r[8]-alphaDrSurf_r[7]-alphaDrSurf_r[6]+alphaDrSurf_r[3]+alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[30] = ser_3x3v_p1_surfvy_quad_30(1, fc); 
  } else { 
    fUpwindQuad_r[30] = ser_3x3v_p1_surfvy_quad_30(-1, fr); 
  } 
  if (alphaDrSurf_r[16]+alphaDrSurf_r[8]+alphaDrSurf_r[7]+alphaDrSurf_r[6]+alphaDrSurf_r[3]+alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[31] = ser_3x3v_p1_surfvy_quad_31(1, fc); 
  } else { 
    fUpwindQuad_r[31] = ser_3x3v_p1_surfvy_quad_31(-1, fr); 
  } 
  fUpwind_l[0] = 0.1767766952966368*(fUpwindQuad_l[31]+fUpwindQuad_l[30]+fUpwindQuad_l[29]+fUpwindQuad_l[28]+fUpwindQuad_l[27]+fUpwindQuad_l[26]+fUpwindQuad_l[25]+fUpwindQuad_l[24]+fUpwindQuad_l[23]+fUpwindQuad_l[22]+fUpwindQuad_l[21]+fUpwindQuad_l[20]+fUpwindQuad_l[19]+fUpwindQuad_l[18]+fUpwindQuad_l[17]+fUpwindQuad_l[16]+fUpwindQuad_l[15]+fUpwindQuad_l[14]+fUpwindQuad_l[13]+fUpwindQuad_l[12]+fUpwindQuad_l[11]+fUpwindQuad_l[10]+fUpwindQuad_l[9]+fUpwindQuad_l[8]+fUpwindQuad_l[7]+fUpwindQuad_l[6]+fUpwindQuad_l[5]+fUpwindQuad_l[4]+fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[1] = 0.1767766952966368*(fUpwindQuad_l[31]-1.0*fUpwindQuad_l[30]+fUpwindQuad_l[29]-1.0*fUpwindQuad_l[28]+fUpwindQuad_l[27]-1.0*fUpwindQuad_l[26]+fUpwindQuad_l[25]-1.0*fUpwindQuad_l[24]+fUpwindQuad_l[23]-1.0*fUpwindQuad_l[22]+fUpwindQuad_l[21]-1.0*fUpwindQuad_l[20]+fUpwindQuad_l[19]-1.0*fUpwindQuad_l[18]+fUpwindQuad_l[17]-1.0*fUpwindQuad_l[16]+fUpwindQuad_l[15]-1.0*fUpwindQuad_l[14]+fUpwindQuad_l[13]-1.0*fUpwindQuad_l[12]+fUpwindQuad_l[11]-1.0*fUpwindQuad_l[10]+fUpwindQuad_l[9]-1.0*fUpwindQuad_l[8]+fUpwindQuad_l[7]-1.0*fUpwindQuad_l[6]+fUpwindQuad_l[5]-1.0*fUpwindQuad_l[4]+fUpwindQuad_l[3]-1.0*fUpwindQuad_l[2]+fUpwindQuad_l[1]-1.0*fUpwindQuad_l[0]); 
  fUpwind_l[2] = 0.1767766952966368*(fUpwindQuad_l[31]+fUpwindQuad_l[30]-1.0*(fUpwindQuad_l[29]+fUpwindQuad_l[28])+fUpwindQuad_l[27]+fUpwindQuad_l[26]-1.0*(fUpwindQuad_l[25]+fUpwindQuad_l[24])+fUpwindQuad_l[23]+fUpwindQuad_l[22]-1.0*(fUpwindQuad_l[21]+fUpwindQuad_l[20])+fUpwindQuad_l[19]+fUpwindQuad_l[18]-1.0*(fUpwindQuad_l[17]+fUpwindQuad_l[16])+fUpwindQuad_l[15]+fUpwindQuad_l[14]-1.0*(fUpwindQuad_l[13]+fUpwindQuad_l[12])+fUpwindQuad_l[11]+fUpwindQuad_l[10]-1.0*(fUpwindQuad_l[9]+fUpwindQuad_l[8])+fUpwindQuad_l[7]+fUpwindQuad_l[6]-1.0*(fUpwindQuad_l[5]+fUpwindQuad_l[4])+fUpwindQuad_l[3]+fUpwindQuad_l[2]-1.0*(fUpwindQuad_l[1]+fUpwindQuad_l[0])); 
  fUpwind_l[3] = 0.1767766952966368*(fUpwindQuad_l[31]+fUpwindQuad_l[30]+fUpwindQuad_l[29]+fUpwindQuad_l[28]-1.0*(fUpwindQuad_l[27]+fUpwindQuad_l[26]+fUpwindQuad_l[25]+fUpwindQuad_l[24])+fUpwindQuad_l[23]+fUpwindQuad_l[22]+fUpwindQuad_l[21]+fUpwindQuad_l[20]-1.0*(fUpwindQuad_l[19]+fUpwindQuad_l[18]+fUpwindQuad_l[17]+fUpwindQuad_l[16])+fUpwindQuad_l[15]+fUpwindQuad_l[14]+fUpwindQuad_l[13]+fUpwindQuad_l[12]-1.0*(fUpwindQuad_l[11]+fUpwindQuad_l[10]+fUpwindQuad_l[9]+fUpwindQuad_l[8])+fUpwindQuad_l[7]+fUpwindQuad_l[6]+fUpwindQuad_l[5]+fUpwindQuad_l[4]-1.0*(fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]+fUpwindQuad_l[0])); 
  fUpwind_l[4] = 0.1767766952966368*(fUpwindQuad_l[31]+fUpwindQuad_l[30]+fUpwindQuad_l[29]+fUpwindQuad_l[28]+fUpwindQuad_l[27]+fUpwindQuad_l[26]+fUpwindQuad_l[25]+fUpwindQuad_l[24]-1.0*(fUpwindQuad_l[23]+fUpwindQuad_l[22]+fUpwindQuad_l[21]+fUpwindQuad_l[20]+fUpwindQuad_l[19]+fUpwindQuad_l[18]+fUpwindQuad_l[17]+fUpwindQuad_l[16])+fUpwindQuad_l[15]+fUpwindQuad_l[14]+fUpwindQuad_l[13]+fUpwindQuad_l[12]+fUpwindQuad_l[11]+fUpwindQuad_l[10]+fUpwindQuad_l[9]+fUpwindQuad_l[8]-1.0*(fUpwindQuad_l[7]+fUpwindQuad_l[6]+fUpwindQuad_l[5]+fUpwindQuad_l[4]+fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]+fUpwindQuad_l[0])); 
  fUpwind_l[5] = 0.1767766952966368*(fUpwindQuad_l[31]+fUpwindQuad_l[30]+fUpwindQuad_l[29]+fUpwindQuad_l[28]+fUpwindQuad_l[27]+fUpwindQuad_l[26]+fUpwindQuad_l[25]+fUpwindQuad_l[24]+fUpwindQuad_l[23]+fUpwindQuad_l[22]+fUpwindQuad_l[21]+fUpwindQuad_l[20]+fUpwindQuad_l[19]+fUpwindQuad_l[18]+fUpwindQuad_l[17]+fUpwindQuad_l[16]-1.0*(fUpwindQuad_l[15]+fUpwindQuad_l[14]+fUpwindQuad_l[13]+fUpwindQuad_l[12]+fUpwindQuad_l[11]+fUpwindQuad_l[10]+fUpwindQuad_l[9]+fUpwindQuad_l[8]+fUpwindQuad_l[7]+fUpwindQuad_l[6]+fUpwindQuad_l[5]+fUpwindQuad_l[4]+fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]+fUpwindQuad_l[0])); 
  fUpwind_l[6] = 0.1767766952966368*(fUpwindQuad_l[31]-1.0*(fUpwindQuad_l[30]+fUpwindQuad_l[29])+fUpwindQuad_l[28]+fUpwindQuad_l[27]-1.0*(fUpwindQuad_l[26]+fUpwindQuad_l[25])+fUpwindQuad_l[24]+fUpwindQuad_l[23]-1.0*(fUpwindQuad_l[22]+fUpwindQuad_l[21])+fUpwindQuad_l[20]+fUpwindQuad_l[19]-1.0*(fUpwindQuad_l[18]+fUpwindQuad_l[17])+fUpwindQuad_l[16]+fUpwindQuad_l[15]-1.0*(fUpwindQuad_l[14]+fUpwindQuad_l[13])+fUpwindQuad_l[12]+fUpwindQuad_l[11]-1.0*(fUpwindQuad_l[10]+fUpwindQuad_l[9])+fUpwindQuad_l[8]+fUpwindQuad_l[7]-1.0*(fUpwindQuad_l[6]+fUpwindQuad_l[5])+fUpwindQuad_l[4]+fUpwindQuad_l[3]-1.0*(fUpwindQuad_l[2]+fUpwindQuad_l[1])+fUpwindQuad_l[0]); 
  fUpwind_l[7] = 0.1767766952966368*(fUpwindQuad_l[31]-1.0*fUpwindQuad_l[30]+fUpwindQuad_l[29]-1.0*(fUpwindQuad_l[28]+fUpwindQuad_l[27])+fUpwindQuad_l[26]-1.0*fUpwindQuad_l[25]+fUpwindQuad_l[24]+fUpwindQuad_l[23]-1.0*fUpwindQuad_l[22]+fUpwindQuad_l[21]-1.0*(fUpwindQuad_l[20]+fUpwindQuad_l[19])+fUpwindQuad_l[18]-1.0*fUpwindQuad_l[17]+fUpwindQuad_l[16]+fUpwindQuad_l[15]-1.0*fUpwindQuad_l[14]+fUpwindQuad_l[13]-1.0*(fUpwindQuad_l[12]+fUpwindQuad_l[11])+fUpwindQuad_l[10]-1.0*fUpwindQuad_l[9]+fUpwindQuad_l[8]+fUpwindQuad_l[7]-1.0*fUpwindQuad_l[6]+fUpwindQuad_l[5]-1.0*(fUpwindQuad_l[4]+fUpwindQuad_l[3])+fUpwindQuad_l[2]-1.0*fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[8] = 0.1767766952966368*(fUpwindQuad_l[31]+fUpwindQuad_l[30]-1.0*(fUpwindQuad_l[29]+fUpwindQuad_l[28]+fUpwindQuad_l[27]+fUpwindQuad_l[26])+fUpwindQuad_l[25]+fUpwindQuad_l[24]+fUpwindQuad_l[23]+fUpwindQuad_l[22]-1.0*(fUpwindQuad_l[21]+fUpwindQuad_l[20]+fUpwindQuad_l[19]+fUpwindQuad_l[18])+fUpwindQuad_l[17]+fUpwindQuad_l[16]+fUpwindQuad_l[15]+fUpwindQuad_l[14]-1.0*(fUpwindQuad_l[13]+fUpwindQuad_l[12]+fUpwindQuad_l[11]+fUpwindQuad_l[10])+fUpwindQuad_l[9]+fUpwindQuad_l[8]+fUpwindQuad_l[7]+fUpwindQuad_l[6]-1.0*(fUpwindQuad_l[5]+fUpwindQuad_l[4]+fUpwindQuad_l[3]+fUpwindQuad_l[2])+fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[9] = 0.1767766952966368*(fUpwindQuad_l[31]-1.0*fUpwindQuad_l[30]+fUpwindQuad_l[29]-1.0*fUpwindQuad_l[28]+fUpwindQuad_l[27]-1.0*fUpwindQuad_l[26]+fUpwindQuad_l[25]-1.0*(fUpwindQuad_l[24]+fUpwindQuad_l[23])+fUpwindQuad_l[22]-1.0*fUpwindQuad_l[21]+fUpwindQuad_l[20]-1.0*fUpwindQuad_l[19]+fUpwindQuad_l[18]-1.0*fUpwindQuad_l[17]+fUpwindQuad_l[16]+fUpwindQuad_l[15]-1.0*fUpwindQuad_l[14]+fUpwindQuad_l[13]-1.0*fUpwindQuad_l[12]+fUpwindQuad_l[11]-1.0*fUpwindQuad_l[10]+fUpwindQuad_l[9]-1.0*(fUpwindQuad_l[8]+fUpwindQuad_l[7])+fUpwindQuad_l[6]-1.0*fUpwindQuad_l[5]+fUpwindQuad_l[4]-1.0*fUpwindQuad_l[3]+fUpwindQuad_l[2]-1.0*fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[10] = 0.1767766952966368*(fUpwindQuad_l[31]+fUpwindQuad_l[30]-1.0*(fUpwindQuad_l[29]+fUpwindQuad_l[28])+fUpwindQuad_l[27]+fUpwindQuad_l[26]-1.0*(fUpwindQuad_l[25]+fUpwindQuad_l[24]+fUpwindQuad_l[23]+fUpwindQuad_l[22])+fUpwindQuad_l[21]+fUpwindQuad_l[20]-1.0*(fUpwindQuad_l[19]+fUpwindQuad_l[18])+fUpwindQuad_l[17]+fUpwindQuad_l[16]+fUpwindQuad_l[15]+fUpwindQuad_l[14]-1.0*(fUpwindQuad_l[13]+fUpwindQuad_l[12])+fUpwindQuad_l[11]+fUpwindQuad_l[10]-1.0*(fUpwindQuad_l[9]+fUpwindQuad_l[8]+fUpwindQuad_l[7]+fUpwindQuad_l[6])+fUpwindQuad_l[5]+fUpwindQuad_l[4]-1.0*(fUpwindQuad_l[3]+fUpwindQuad_l[2])+fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[11] = 0.1767766952966368*(fUpwindQuad_l[31]+fUpwindQuad_l[30]+fUpwindQuad_l[29]+fUpwindQuad_l[28]-1.0*(fUpwindQuad_l[27]+fUpwindQuad_l[26]+fUpwindQuad_l[25]+fUpwindQuad_l[24]+fUpwindQuad_l[23]+fUpwindQuad_l[22]+fUpwindQuad_l[21]+fUpwindQuad_l[20])+fUpwindQuad_l[19]+fUpwindQuad_l[18]+fUpwindQuad_l[17]+fUpwindQuad_l[16]+fUpwindQuad_l[15]+fUpwindQuad_l[14]+fUpwindQuad_l[13]+fUpwindQuad_l[12]-1.0*(fUpwindQuad_l[11]+fUpwindQuad_l[10]+fUpwindQuad_l[9]+fUpwindQuad_l[8]+fUpwindQuad_l[7]+fUpwindQuad_l[6]+fUpwindQuad_l[5]+fUpwindQuad_l[4])+fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[12] = 0.1767766952966368*(fUpwindQuad_l[31]-1.0*fUpwindQuad_l[30]+fUpwindQuad_l[29]-1.0*fUpwindQuad_l[28]+fUpwindQuad_l[27]-1.0*fUpwindQuad_l[26]+fUpwindQuad_l[25]-1.0*fUpwindQuad_l[24]+fUpwindQuad_l[23]-1.0*fUpwindQuad_l[22]+fUpwindQuad_l[21]-1.0*fUpwindQuad_l[20]+fUpwindQuad_l[19]-1.0*fUpwindQuad_l[18]+fUpwindQuad_l[17]-1.0*(fUpwindQuad_l[16]+fUpwindQuad_l[15])+fUpwindQuad_l[14]-1.0*fUpwindQuad_l[13]+fUpwindQuad_l[12]-1.0*fUpwindQuad_l[11]+fUpwindQuad_l[10]-1.0*fUpwindQuad_l[9]+fUpwindQuad_l[8]-1.0*fUpwindQuad_l[7]+fUpwindQuad_l[6]-1.0*fUpwindQuad_l[5]+fUpwindQuad_l[4]-1.0*fUpwindQuad_l[3]+fUpwindQuad_l[2]-1.0*fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[13] = 0.1767766952966368*(fUpwindQuad_l[31]+fUpwindQuad_l[30]-1.0*(fUpwindQuad_l[29]+fUpwindQuad_l[28])+fUpwindQuad_l[27]+fUpwindQuad_l[26]-1.0*(fUpwindQuad_l[25]+fUpwindQuad_l[24])+fUpwindQuad_l[23]+fUpwindQuad_l[22]-1.0*(fUpwindQuad_l[21]+fUpwindQuad_l[20])+fUpwindQuad_l[19]+fUpwindQuad_l[18]-1.0*(fUpwindQuad_l[17]+fUpwindQuad_l[16]+fUpwindQuad_l[15]+fUpwindQuad_l[14])+fUpwindQuad_l[13]+fUpwindQuad_l[12]-1.0*(fUpwindQuad_l[11]+fUpwindQuad_l[10])+fUpwindQuad_l[9]+fUpwindQuad_l[8]-1.0*(fUpwindQuad_l[7]+fUpwindQuad_l[6])+fUpwindQuad_l[5]+fUpwindQuad_l[4]-1.0*(fUpwindQuad_l[3]+fUpwindQuad_l[2])+fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[14] = 0.1767766952966368*(fUpwindQuad_l[31]+fUpwindQuad_l[30]+fUpwindQuad_l[29]+fUpwindQuad_l[28]-1.0*(fUpwindQuad_l[27]+fUpwindQuad_l[26]+fUpwindQuad_l[25]+fUpwindQuad_l[24])+fUpwindQuad_l[23]+fUpwindQuad_l[22]+fUpwindQuad_l[21]+fUpwindQuad_l[20]-1.0*(fUpwindQuad_l[19]+fUpwindQuad_l[18]+fUpwindQuad_l[17]+fUpwindQuad_l[16]+fUpwindQuad_l[15]+fUpwindQuad_l[14]+fUpwindQuad_l[13]+fUpwindQuad_l[12])+fUpwindQuad_l[11]+fUpwindQuad_l[10]+fUpwindQuad_l[9]+fUpwindQuad_l[8]-1.0*(fUpwindQuad_l[7]+fUpwindQuad_l[6]+fUpwindQuad_l[5]+fUpwindQuad_l[4])+fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[15] = 0.1767766952966368*(fUpwindQuad_l[31]+fUpwindQuad_l[30]+fUpwindQuad_l[29]+fUpwindQuad_l[28]+fUpwindQuad_l[27]+fUpwindQuad_l[26]+fUpwindQuad_l[25]+fUpwindQuad_l[24]-1.0*(fUpwindQuad_l[23]+fUpwindQuad_l[22]+fUpwindQuad_l[21]+fUpwindQuad_l[20]+fUpwindQuad_l[19]+fUpwindQuad_l[18]+fUpwindQuad_l[17]+fUpwindQuad_l[16]+fUpwindQuad_l[15]+fUpwindQuad_l[14]+fUpwindQuad_l[13]+fUpwindQuad_l[12]+fUpwindQuad_l[11]+fUpwindQuad_l[10]+fUpwindQuad_l[9]+fUpwindQuad_l[8])+fUpwindQuad_l[7]+fUpwindQuad_l[6]+fUpwindQuad_l[5]+fUpwindQuad_l[4]+fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[16] = 0.1767766952966368*(fUpwindQuad_l[31]-1.0*(fUpwindQuad_l[30]+fUpwindQuad_l[29])+fUpwindQuad_l[28]-1.0*fUpwindQuad_l[27]+fUpwindQuad_l[26]+fUpwindQuad_l[25]-1.0*fUpwindQuad_l[24]+fUpwindQuad_l[23]-1.0*(fUpwindQuad_l[22]+fUpwindQuad_l[21])+fUpwindQuad_l[20]-1.0*fUpwindQuad_l[19]+fUpwindQuad_l[18]+fUpwindQuad_l[17]-1.0*fUpwindQuad_l[16]+fUpwindQuad_l[15]-1.0*(fUpwindQuad_l[14]+fUpwindQuad_l[13])+fUpwindQuad_l[12]-1.0*fUpwindQuad_l[11]+fUpwindQuad_l[10]+fUpwindQuad_l[9]-1.0*fUpwindQuad_l[8]+fUpwindQuad_l[7]-1.0*(fUpwindQuad_l[6]+fUpwindQuad_l[5])+fUpwindQuad_l[4]-1.0*fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]-1.0*fUpwindQuad_l[0]); 
  fUpwind_l[17] = 0.1767766952966368*(fUpwindQuad_l[31]-1.0*(fUpwindQuad_l[30]+fUpwindQuad_l[29])+fUpwindQuad_l[28]+fUpwindQuad_l[27]-1.0*(fUpwindQuad_l[26]+fUpwindQuad_l[25])+fUpwindQuad_l[24]-1.0*fUpwindQuad_l[23]+fUpwindQuad_l[22]+fUpwindQuad_l[21]-1.0*(fUpwindQuad_l[20]+fUpwindQuad_l[19])+fUpwindQuad_l[18]+fUpwindQuad_l[17]-1.0*fUpwindQuad_l[16]+fUpwindQuad_l[15]-1.0*(fUpwindQuad_l[14]+fUpwindQuad_l[13])+fUpwindQuad_l[12]+fUpwindQuad_l[11]-1.0*(fUpwindQuad_l[10]+fUpwindQuad_l[9])+fUpwindQuad_l[8]-1.0*fUpwindQuad_l[7]+fUpwindQuad_l[6]+fUpwindQuad_l[5]-1.0*(fUpwindQuad_l[4]+fUpwindQuad_l[3])+fUpwindQuad_l[2]+fUpwindQuad_l[1]-1.0*fUpwindQuad_l[0]); 
  fUpwind_l[18] = 0.1767766952966368*(fUpwindQuad_l[31]-1.0*fUpwindQuad_l[30]+fUpwindQuad_l[29]-1.0*(fUpwindQuad_l[28]+fUpwindQuad_l[27])+fUpwindQuad_l[26]-1.0*fUpwindQuad_l[25]+fUpwindQuad_l[24]-1.0*fUpwindQuad_l[23]+fUpwindQuad_l[22]-1.0*fUpwindQuad_l[21]+fUpwindQuad_l[20]+fUpwindQuad_l[19]-1.0*fUpwindQuad_l[18]+fUpwindQuad_l[17]-1.0*fUpwindQuad_l[16]+fUpwindQuad_l[15]-1.0*fUpwindQuad_l[14]+fUpwindQuad_l[13]-1.0*(fUpwindQuad_l[12]+fUpwindQuad_l[11])+fUpwindQuad_l[10]-1.0*fUpwindQuad_l[9]+fUpwindQuad_l[8]-1.0*fUpwindQuad_l[7]+fUpwindQuad_l[6]-1.0*fUpwindQuad_l[5]+fUpwindQuad_l[4]+fUpwindQuad_l[3]-1.0*fUpwindQuad_l[2]+fUpwindQuad_l[1]-1.0*fUpwindQuad_l[0]); 
  fUpwind_l[19] = 0.1767766952966368*(fUpwindQuad_l[31]+fUpwindQuad_l[30]-1.0*(fUpwindQuad_l[29]+fUpwindQuad_l[28]+fUpwindQuad_l[27]+fUpwindQuad_l[26])+fUpwindQuad_l[25]+fUpwindQuad_l[24]-1.0*(fUpwindQuad_l[23]+fUpwindQuad_l[22])+fUpwindQuad_l[21]+fUpwindQuad_l[20]+fUpwindQuad_l[19]+fUpwindQuad_l[18]-1.0*(fUpwindQuad_l[17]+fUpwindQuad_l[16])+fUpwindQuad_l[15]+fUpwindQuad_l[14]-1.0*(fUpwindQuad_l[13]+fUpwindQuad_l[12]+fUpwindQuad_l[11]+fUpwindQuad_l[10])+fUpwindQuad_l[9]+fUpwindQuad_l[8]-1.0*(fUpwindQuad_l[7]+fUpwindQuad_l[6])+fUpwindQuad_l[5]+fUpwindQuad_l[4]+fUpwindQuad_l[3]+fUpwindQuad_l[2]-1.0*(fUpwindQuad_l[1]+fUpwindQuad_l[0])); 
  fUpwind_l[20] = 0.1767766952966368*(fUpwindQuad_l[31]-1.0*(fUpwindQuad_l[30]+fUpwindQuad_l[29])+fUpwindQuad_l[28]+fUpwindQuad_l[27]-1.0*(fUpwindQuad_l[26]+fUpwindQuad_l[25])+fUpwindQuad_l[24]+fUpwindQuad_l[23]-1.0*(fUpwindQuad_l[22]+fUpwindQuad_l[21])+fUpwindQuad_l[20]+fUpwindQuad_l[19]-1.0*(fUpwindQuad_l[18]+fUpwindQuad_l[17])+fUpwindQuad_l[16]-1.0*fUpwindQuad_l[15]+fUpwindQuad_l[14]+fUpwindQuad_l[13]-1.0*(fUpwindQuad_l[12]+fUpwindQuad_l[11])+fUpwindQuad_l[10]+fUpwindQuad_l[9]-1.0*(fUpwindQuad_l[8]+fUpwindQuad_l[7])+fUpwindQuad_l[6]+fUpwindQuad_l[5]-1.0*(fUpwindQuad_l[4]+fUpwindQuad_l[3])+fUpwindQuad_l[2]+fUpwindQuad_l[1]-1.0*fUpwindQuad_l[0]); 
  fUpwind_l[21] = 0.1767766952966368*(fUpwindQuad_l[31]-1.0*fUpwindQuad_l[30]+fUpwindQuad_l[29]-1.0*(fUpwindQuad_l[28]+fUpwindQuad_l[27])+fUpwindQuad_l[26]-1.0*fUpwindQuad_l[25]+fUpwindQuad_l[24]+fUpwindQuad_l[23]-1.0*fUpwindQuad_l[22]+fUpwindQuad_l[21]-1.0*(fUpwindQuad_l[20]+fUpwindQuad_l[19])+fUpwindQuad_l[18]-1.0*fUpwindQuad_l[17]+fUpwindQuad_l[16]-1.0*fUpwindQuad_l[15]+fUpwindQuad_l[14]-1.0*fUpwindQuad_l[13]+fUpwindQuad_l[12]+fUpwindQuad_l[11]-1.0*fUpwindQuad_l[10]+fUpwindQuad_l[9]-1.0*(fUpwindQuad_l[8]+fUpwindQuad_l[7])+fUpwindQuad_l[6]-1.0*fUpwindQuad_l[5]+fUpwindQuad_l[4]+fUpwindQuad_l[3]-1.0*fUpwindQuad_l[2]+fUpwindQuad_l[1]-1.0*fUpwindQuad_l[0]); 
  fUpwind_l[22] = 0.1767766952966368*(fUpwindQuad_l[31]+fUpwindQuad_l[30]-1.0*(fUpwindQuad_l[29]+fUpwindQuad_l[28]+fUpwindQuad_l[27]+fUpwindQuad_l[26])+fUpwindQuad_l[25]+fUpwindQuad_l[24]+fUpwindQuad_l[23]+fUpwindQuad_l[22]-1.0*(fUpwindQuad_l[21]+fUpwindQuad_l[20]+fUpwindQuad_l[19]+fUpwindQuad_l[18])+fUpwindQuad_l[17]+fUpwindQuad_l[16]-1.0*(fUpwindQuad_l[15]+fUpwindQuad_l[14])+fUpwindQuad_l[13]+fUpwindQuad_l[12]+fUpwindQuad_l[11]+fUpwindQuad_l[10]-1.0*(fUpwindQuad_l[9]+fUpwindQuad_l[8]+fUpwindQuad_l[7]+fUpwindQuad_l[6])+fUpwindQuad_l[5]+fUpwindQuad_l[4]+fUpwindQuad_l[3]+fUpwindQuad_l[2]-1.0*(fUpwindQuad_l[1]+fUpwindQuad_l[0])); 
  fUpwind_l[23] = 0.1767766952966368*(fUpwindQuad_l[31]-1.0*fUpwindQuad_l[30]+fUpwindQuad_l[29]-1.0*fUpwindQuad_l[28]+fUpwindQuad_l[27]-1.0*fUpwindQuad_l[26]+fUpwindQuad_l[25]-1.0*(fUpwindQuad_l[24]+fUpwindQuad_l[23])+fUpwindQuad_l[22]-1.0*fUpwindQuad_l[21]+fUpwindQuad_l[20]-1.0*fUpwindQuad_l[19]+fUpwindQuad_l[18]-1.0*fUpwindQuad_l[17]+fUpwindQuad_l[16]-1.0*fUpwindQuad_l[15]+fUpwindQuad_l[14]-1.0*fUpwindQuad_l[13]+fUpwindQuad_l[12]-1.0*fUpwindQuad_l[11]+fUpwindQuad_l[10]-1.0*fUpwindQuad_l[9]+fUpwindQuad_l[8]+fUpwindQuad_l[7]-1.0*fUpwindQuad_l[6]+fUpwindQuad_l[5]-1.0*fUpwindQuad_l[4]+fUpwindQuad_l[3]-1.0*fUpwindQuad_l[2]+fUpwindQuad_l[1]-1.0*fUpwindQuad_l[0]); 
  fUpwind_l[24] = 0.1767766952966368*(fUpwindQuad_l[31]+fUpwindQuad_l[30]-1.0*(fUpwindQuad_l[29]+fUpwindQuad_l[28])+fUpwindQuad_l[27]+fUpwindQuad_l[26]-1.0*(fUpwindQuad_l[25]+fUpwindQuad_l[24]+fUpwindQuad_l[23]+fUpwindQuad_l[22])+fUpwindQuad_l[21]+fUpwindQuad_l[20]-1.0*(fUpwindQuad_l[19]+fUpwindQuad_l[18])+fUpwindQuad_l[17]+fUpwindQuad_l[16]-1.0*(fUpwindQuad_l[15]+fUpwindQuad_l[14])+fUpwindQuad_l[13]+fUpwindQuad_l[12]-1.0*(fUpwindQuad_l[11]+fUpwindQuad_l[10])+fUpwindQuad_l[9]+fUpwindQuad_l[8]+fUpwindQuad_l[7]+fUpwindQuad_l[6]-1.0*(fUpwindQuad_l[5]+fUpwindQuad_l[4])+fUpwindQuad_l[3]+fUpwindQuad_l[2]-1.0*(fUpwindQuad_l[1]+fUpwindQuad_l[0])); 
  fUpwind_l[25] = 0.1767766952966368*(fUpwindQuad_l[31]+fUpwindQuad_l[30]+fUpwindQuad_l[29]+fUpwindQuad_l[28]-1.0*(fUpwindQuad_l[27]+fUpwindQuad_l[26]+fUpwindQuad_l[25]+fUpwindQuad_l[24]+fUpwindQuad_l[23]+fUpwindQuad_l[22]+fUpwindQuad_l[21]+fUpwindQuad_l[20])+fUpwindQuad_l[19]+fUpwindQuad_l[18]+fUpwindQuad_l[17]+fUpwindQuad_l[16]-1.0*(fUpwindQuad_l[15]+fUpwindQuad_l[14]+fUpwindQuad_l[13]+fUpwindQuad_l[12])+fUpwindQuad_l[11]+fUpwindQuad_l[10]+fUpwindQuad_l[9]+fUpwindQuad_l[8]+fUpwindQuad_l[7]+fUpwindQuad_l[6]+fUpwindQuad_l[5]+fUpwindQuad_l[4]-1.0*(fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]+fUpwindQuad_l[0])); 
  fUpwind_l[26] = 0.1767766952966368*(fUpwindQuad_l[31]-1.0*(fUpwindQuad_l[30]+fUpwindQuad_l[29])+fUpwindQuad_l[28]-1.0*fUpwindQuad_l[27]+fUpwindQuad_l[26]+fUpwindQuad_l[25]-1.0*(fUpwindQuad_l[24]+fUpwindQuad_l[23])+fUpwindQuad_l[22]+fUpwindQuad_l[21]-1.0*fUpwindQuad_l[20]+fUpwindQuad_l[19]-1.0*(fUpwindQuad_l[18]+fUpwindQuad_l[17])+fUpwindQuad_l[16]+fUpwindQuad_l[15]-1.0*(fUpwindQuad_l[14]+fUpwindQuad_l[13])+fUpwindQuad_l[12]-1.0*fUpwindQuad_l[11]+fUpwindQuad_l[10]+fUpwindQuad_l[9]-1.0*(fUpwindQuad_l[8]+fUpwindQuad_l[7])+fUpwindQuad_l[6]+fUpwindQuad_l[5]-1.0*fUpwindQuad_l[4]+fUpwindQuad_l[3]-1.0*(fUpwindQuad_l[2]+fUpwindQuad_l[1])+fUpwindQuad_l[0]); 
  fUpwind_l[27] = 0.1767766952966368*(fUpwindQuad_l[31]-1.0*(fUpwindQuad_l[30]+fUpwindQuad_l[29])+fUpwindQuad_l[28]-1.0*fUpwindQuad_l[27]+fUpwindQuad_l[26]+fUpwindQuad_l[25]-1.0*fUpwindQuad_l[24]+fUpwindQuad_l[23]-1.0*(fUpwindQuad_l[22]+fUpwindQuad_l[21])+fUpwindQuad_l[20]-1.0*fUpwindQuad_l[19]+fUpwindQuad_l[18]+fUpwindQuad_l[17]-1.0*(fUpwindQuad_l[16]+fUpwindQuad_l[15])+fUpwindQuad_l[14]+fUpwindQuad_l[13]-1.0*fUpwindQuad_l[12]+fUpwindQuad_l[11]-1.0*(fUpwindQuad_l[10]+fUpwindQuad_l[9])+fUpwindQuad_l[8]-1.0*fUpwindQuad_l[7]+fUpwindQuad_l[6]+fUpwindQuad_l[5]-1.0*fUpwindQuad_l[4]+fUpwindQuad_l[3]-1.0*(fUpwindQuad_l[2]+fUpwindQuad_l[1])+fUpwindQuad_l[0]); 
  fUpwind_l[28] = 0.1767766952966368*(fUpwindQuad_l[31]-1.0*(fUpwindQuad_l[30]+fUpwindQuad_l[29])+fUpwindQuad_l[28]+fUpwindQuad_l[27]-1.0*(fUpwindQuad_l[26]+fUpwindQuad_l[25])+fUpwindQuad_l[24]-1.0*fUpwindQuad_l[23]+fUpwindQuad_l[22]+fUpwindQuad_l[21]-1.0*(fUpwindQuad_l[20]+fUpwindQuad_l[19])+fUpwindQuad_l[18]+fUpwindQuad_l[17]-1.0*(fUpwindQuad_l[16]+fUpwindQuad_l[15])+fUpwindQuad_l[14]+fUpwindQuad_l[13]-1.0*(fUpwindQuad_l[12]+fUpwindQuad_l[11])+fUpwindQuad_l[10]+fUpwindQuad_l[9]-1.0*fUpwindQuad_l[8]+fUpwindQuad_l[7]-1.0*(fUpwindQuad_l[6]+fUpwindQuad_l[5])+fUpwindQuad_l[4]+fUpwindQuad_l[3]-1.0*(fUpwindQuad_l[2]+fUpwindQuad_l[1])+fUpwindQuad_l[0]); 
  fUpwind_l[29] = 0.1767766952966368*(fUpwindQuad_l[31]-1.0*fUpwindQuad_l[30]+fUpwindQuad_l[29]-1.0*(fUpwindQuad_l[28]+fUpwindQuad_l[27])+fUpwindQuad_l[26]-1.0*fUpwindQuad_l[25]+fUpwindQuad_l[24]-1.0*fUpwindQuad_l[23]+fUpwindQuad_l[22]-1.0*fUpwindQuad_l[21]+fUpwindQuad_l[20]+fUpwindQuad_l[19]-1.0*fUpwindQuad_l[18]+fUpwindQuad_l[17]-1.0*(fUpwindQuad_l[16]+fUpwindQuad_l[15])+fUpwindQuad_l[14]-1.0*fUpwindQuad_l[13]+fUpwindQuad_l[12]+fUpwindQuad_l[11]-1.0*fUpwindQuad_l[10]+fUpwindQuad_l[9]-1.0*fUpwindQuad_l[8]+fUpwindQuad_l[7]-1.0*fUpwindQuad_l[6]+fUpwindQuad_l[5]-1.0*(fUpwindQuad_l[4]+fUpwindQuad_l[3])+fUpwindQuad_l[2]-1.0*fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[30] = 0.1767766952966368*(fUpwindQuad_l[31]+fUpwindQuad_l[30]-1.0*(fUpwindQuad_l[29]+fUpwindQuad_l[28]+fUpwindQuad_l[27]+fUpwindQuad_l[26])+fUpwindQuad_l[25]+fUpwindQuad_l[24]-1.0*(fUpwindQuad_l[23]+fUpwindQuad_l[22])+fUpwindQuad_l[21]+fUpwindQuad_l[20]+fUpwindQuad_l[19]+fUpwindQuad_l[18]-1.0*(fUpwindQuad_l[17]+fUpwindQuad_l[16]+fUpwindQuad_l[15]+fUpwindQuad_l[14])+fUpwindQuad_l[13]+fUpwindQuad_l[12]+fUpwindQuad_l[11]+fUpwindQuad_l[10]-1.0*(fUpwindQuad_l[9]+fUpwindQuad_l[8])+fUpwindQuad_l[7]+fUpwindQuad_l[6]-1.0*(fUpwindQuad_l[5]+fUpwindQuad_l[4]+fUpwindQuad_l[3]+fUpwindQuad_l[2])+fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[31] = 0.1767766952966368*(fUpwindQuad_l[31]-1.0*(fUpwindQuad_l[30]+fUpwindQuad_l[29])+fUpwindQuad_l[28]-1.0*fUpwindQuad_l[27]+fUpwindQuad_l[26]+fUpwindQuad_l[25]-1.0*(fUpwindQuad_l[24]+fUpwindQuad_l[23])+fUpwindQuad_l[22]+fUpwindQuad_l[21]-1.0*fUpwindQuad_l[20]+fUpwindQuad_l[19]-1.0*(fUpwindQuad_l[18]+fUpwindQuad_l[17])+fUpwindQuad_l[16]-1.0*fUpwindQuad_l[15]+fUpwindQuad_l[14]+fUpwindQuad_l[13]-1.0*fUpwindQuad_l[12]+fUpwindQuad_l[11]-1.0*(fUpwindQuad_l[10]+fUpwindQuad_l[9])+fUpwindQuad_l[8]+fUpwindQuad_l[7]-1.0*(fUpwindQuad_l[6]+fUpwindQuad_l[5])+fUpwindQuad_l[4]-1.0*fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]-1.0*fUpwindQuad_l[0]); 

  fUpwind_r[0] = 0.1767766952966368*(fUpwindQuad_r[31]+fUpwindQuad_r[30]+fUpwindQuad_r[29]+fUpwindQuad_r[28]+fUpwindQuad_r[27]+fUpwindQuad_r[26]+fUpwindQuad_r[25]+fUpwindQuad_r[24]+fUpwindQuad_r[23]+fUpwindQuad_r[22]+fUpwindQuad_r[21]+fUpwindQuad_r[20]+fUpwindQuad_r[19]+fUpwindQuad_r[18]+fUpwindQuad_r[17]+fUpwindQuad_r[16]+fUpwindQuad_r[15]+fUpwindQuad_r[14]+fUpwindQuad_r[13]+fUpwindQuad_r[12]+fUpwindQuad_r[11]+fUpwindQuad_r[10]+fUpwindQuad_r[9]+fUpwindQuad_r[8]+fUpwindQuad_r[7]+fUpwindQuad_r[6]+fUpwindQuad_r[5]+fUpwindQuad_r[4]+fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[1] = 0.1767766952966368*(fUpwindQuad_r[31]-1.0*fUpwindQuad_r[30]+fUpwindQuad_r[29]-1.0*fUpwindQuad_r[28]+fUpwindQuad_r[27]-1.0*fUpwindQuad_r[26]+fUpwindQuad_r[25]-1.0*fUpwindQuad_r[24]+fUpwindQuad_r[23]-1.0*fUpwindQuad_r[22]+fUpwindQuad_r[21]-1.0*fUpwindQuad_r[20]+fUpwindQuad_r[19]-1.0*fUpwindQuad_r[18]+fUpwindQuad_r[17]-1.0*fUpwindQuad_r[16]+fUpwindQuad_r[15]-1.0*fUpwindQuad_r[14]+fUpwindQuad_r[13]-1.0*fUpwindQuad_r[12]+fUpwindQuad_r[11]-1.0*fUpwindQuad_r[10]+fUpwindQuad_r[9]-1.0*fUpwindQuad_r[8]+fUpwindQuad_r[7]-1.0*fUpwindQuad_r[6]+fUpwindQuad_r[5]-1.0*fUpwindQuad_r[4]+fUpwindQuad_r[3]-1.0*fUpwindQuad_r[2]+fUpwindQuad_r[1]-1.0*fUpwindQuad_r[0]); 
  fUpwind_r[2] = 0.1767766952966368*(fUpwindQuad_r[31]+fUpwindQuad_r[30]-1.0*(fUpwindQuad_r[29]+fUpwindQuad_r[28])+fUpwindQuad_r[27]+fUpwindQuad_r[26]-1.0*(fUpwindQuad_r[25]+fUpwindQuad_r[24])+fUpwindQuad_r[23]+fUpwindQuad_r[22]-1.0*(fUpwindQuad_r[21]+fUpwindQuad_r[20])+fUpwindQuad_r[19]+fUpwindQuad_r[18]-1.0*(fUpwindQuad_r[17]+fUpwindQuad_r[16])+fUpwindQuad_r[15]+fUpwindQuad_r[14]-1.0*(fUpwindQuad_r[13]+fUpwindQuad_r[12])+fUpwindQuad_r[11]+fUpwindQuad_r[10]-1.0*(fUpwindQuad_r[9]+fUpwindQuad_r[8])+fUpwindQuad_r[7]+fUpwindQuad_r[6]-1.0*(fUpwindQuad_r[5]+fUpwindQuad_r[4])+fUpwindQuad_r[3]+fUpwindQuad_r[2]-1.0*(fUpwindQuad_r[1]+fUpwindQuad_r[0])); 
  fUpwind_r[3] = 0.1767766952966368*(fUpwindQuad_r[31]+fUpwindQuad_r[30]+fUpwindQuad_r[29]+fUpwindQuad_r[28]-1.0*(fUpwindQuad_r[27]+fUpwindQuad_r[26]+fUpwindQuad_r[25]+fUpwindQuad_r[24])+fUpwindQuad_r[23]+fUpwindQuad_r[22]+fUpwindQuad_r[21]+fUpwindQuad_r[20]-1.0*(fUpwindQuad_r[19]+fUpwindQuad_r[18]+fUpwindQuad_r[17]+fUpwindQuad_r[16])+fUpwindQuad_r[15]+fUpwindQuad_r[14]+fUpwindQuad_r[13]+fUpwindQuad_r[12]-1.0*(fUpwindQuad_r[11]+fUpwindQuad_r[10]+fUpwindQuad_r[9]+fUpwindQuad_r[8])+fUpwindQuad_r[7]+fUpwindQuad_r[6]+fUpwindQuad_r[5]+fUpwindQuad_r[4]-1.0*(fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]+fUpwindQuad_r[0])); 
  fUpwind_r[4] = 0.1767766952966368*(fUpwindQuad_r[31]+fUpwindQuad_r[30]+fUpwindQuad_r[29]+fUpwindQuad_r[28]+fUpwindQuad_r[27]+fUpwindQuad_r[26]+fUpwindQuad_r[25]+fUpwindQuad_r[24]-1.0*(fUpwindQuad_r[23]+fUpwindQuad_r[22]+fUpwindQuad_r[21]+fUpwindQuad_r[20]+fUpwindQuad_r[19]+fUpwindQuad_r[18]+fUpwindQuad_r[17]+fUpwindQuad_r[16])+fUpwindQuad_r[15]+fUpwindQuad_r[14]+fUpwindQuad_r[13]+fUpwindQuad_r[12]+fUpwindQuad_r[11]+fUpwindQuad_r[10]+fUpwindQuad_r[9]+fUpwindQuad_r[8]-1.0*(fUpwindQuad_r[7]+fUpwindQuad_r[6]+fUpwindQuad_r[5]+fUpwindQuad_r[4]+fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]+fUpwindQuad_r[0])); 
  fUpwind_r[5] = 0.1767766952966368*(fUpwindQuad_r[31]+fUpwindQuad_r[30]+fUpwindQuad_r[29]+fUpwindQuad_r[28]+fUpwindQuad_r[27]+fUpwindQuad_r[26]+fUpwindQuad_r[25]+fUpwindQuad_r[24]+fUpwindQuad_r[23]+fUpwindQuad_r[22]+fUpwindQuad_r[21]+fUpwindQuad_r[20]+fUpwindQuad_r[19]+fUpwindQuad_r[18]+fUpwindQuad_r[17]+fUpwindQuad_r[16]-1.0*(fUpwindQuad_r[15]+fUpwindQuad_r[14]+fUpwindQuad_r[13]+fUpwindQuad_r[12]+fUpwindQuad_r[11]+fUpwindQuad_r[10]+fUpwindQuad_r[9]+fUpwindQuad_r[8]+fUpwindQuad_r[7]+fUpwindQuad_r[6]+fUpwindQuad_r[5]+fUpwindQuad_r[4]+fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]+fUpwindQuad_r[0])); 
  fUpwind_r[6] = 0.1767766952966368*(fUpwindQuad_r[31]-1.0*(fUpwindQuad_r[30]+fUpwindQuad_r[29])+fUpwindQuad_r[28]+fUpwindQuad_r[27]-1.0*(fUpwindQuad_r[26]+fUpwindQuad_r[25])+fUpwindQuad_r[24]+fUpwindQuad_r[23]-1.0*(fUpwindQuad_r[22]+fUpwindQuad_r[21])+fUpwindQuad_r[20]+fUpwindQuad_r[19]-1.0*(fUpwindQuad_r[18]+fUpwindQuad_r[17])+fUpwindQuad_r[16]+fUpwindQuad_r[15]-1.0*(fUpwindQuad_r[14]+fUpwindQuad_r[13])+fUpwindQuad_r[12]+fUpwindQuad_r[11]-1.0*(fUpwindQuad_r[10]+fUpwindQuad_r[9])+fUpwindQuad_r[8]+fUpwindQuad_r[7]-1.0*(fUpwindQuad_r[6]+fUpwindQuad_r[5])+fUpwindQuad_r[4]+fUpwindQuad_r[3]-1.0*(fUpwindQuad_r[2]+fUpwindQuad_r[1])+fUpwindQuad_r[0]); 
  fUpwind_r[7] = 0.1767766952966368*(fUpwindQuad_r[31]-1.0*fUpwindQuad_r[30]+fUpwindQuad_r[29]-1.0*(fUpwindQuad_r[28]+fUpwindQuad_r[27])+fUpwindQuad_r[26]-1.0*fUpwindQuad_r[25]+fUpwindQuad_r[24]+fUpwindQuad_r[23]-1.0*fUpwindQuad_r[22]+fUpwindQuad_r[21]-1.0*(fUpwindQuad_r[20]+fUpwindQuad_r[19])+fUpwindQuad_r[18]-1.0*fUpwindQuad_r[17]+fUpwindQuad_r[16]+fUpwindQuad_r[15]-1.0*fUpwindQuad_r[14]+fUpwindQuad_r[13]-1.0*(fUpwindQuad_r[12]+fUpwindQuad_r[11])+fUpwindQuad_r[10]-1.0*fUpwindQuad_r[9]+fUpwindQuad_r[8]+fUpwindQuad_r[7]-1.0*fUpwindQuad_r[6]+fUpwindQuad_r[5]-1.0*(fUpwindQuad_r[4]+fUpwindQuad_r[3])+fUpwindQuad_r[2]-1.0*fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[8] = 0.1767766952966368*(fUpwindQuad_r[31]+fUpwindQuad_r[30]-1.0*(fUpwindQuad_r[29]+fUpwindQuad_r[28]+fUpwindQuad_r[27]+fUpwindQuad_r[26])+fUpwindQuad_r[25]+fUpwindQuad_r[24]+fUpwindQuad_r[23]+fUpwindQuad_r[22]-1.0*(fUpwindQuad_r[21]+fUpwindQuad_r[20]+fUpwindQuad_r[19]+fUpwindQuad_r[18])+fUpwindQuad_r[17]+fUpwindQuad_r[16]+fUpwindQuad_r[15]+fUpwindQuad_r[14]-1.0*(fUpwindQuad_r[13]+fUpwindQuad_r[12]+fUpwindQuad_r[11]+fUpwindQuad_r[10])+fUpwindQuad_r[9]+fUpwindQuad_r[8]+fUpwindQuad_r[7]+fUpwindQuad_r[6]-1.0*(fUpwindQuad_r[5]+fUpwindQuad_r[4]+fUpwindQuad_r[3]+fUpwindQuad_r[2])+fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[9] = 0.1767766952966368*(fUpwindQuad_r[31]-1.0*fUpwindQuad_r[30]+fUpwindQuad_r[29]-1.0*fUpwindQuad_r[28]+fUpwindQuad_r[27]-1.0*fUpwindQuad_r[26]+fUpwindQuad_r[25]-1.0*(fUpwindQuad_r[24]+fUpwindQuad_r[23])+fUpwindQuad_r[22]-1.0*fUpwindQuad_r[21]+fUpwindQuad_r[20]-1.0*fUpwindQuad_r[19]+fUpwindQuad_r[18]-1.0*fUpwindQuad_r[17]+fUpwindQuad_r[16]+fUpwindQuad_r[15]-1.0*fUpwindQuad_r[14]+fUpwindQuad_r[13]-1.0*fUpwindQuad_r[12]+fUpwindQuad_r[11]-1.0*fUpwindQuad_r[10]+fUpwindQuad_r[9]-1.0*(fUpwindQuad_r[8]+fUpwindQuad_r[7])+fUpwindQuad_r[6]-1.0*fUpwindQuad_r[5]+fUpwindQuad_r[4]-1.0*fUpwindQuad_r[3]+fUpwindQuad_r[2]-1.0*fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[10] = 0.1767766952966368*(fUpwindQuad_r[31]+fUpwindQuad_r[30]-1.0*(fUpwindQuad_r[29]+fUpwindQuad_r[28])+fUpwindQuad_r[27]+fUpwindQuad_r[26]-1.0*(fUpwindQuad_r[25]+fUpwindQuad_r[24]+fUpwindQuad_r[23]+fUpwindQuad_r[22])+fUpwindQuad_r[21]+fUpwindQuad_r[20]-1.0*(fUpwindQuad_r[19]+fUpwindQuad_r[18])+fUpwindQuad_r[17]+fUpwindQuad_r[16]+fUpwindQuad_r[15]+fUpwindQuad_r[14]-1.0*(fUpwindQuad_r[13]+fUpwindQuad_r[12])+fUpwindQuad_r[11]+fUpwindQuad_r[10]-1.0*(fUpwindQuad_r[9]+fUpwindQuad_r[8]+fUpwindQuad_r[7]+fUpwindQuad_r[6])+fUpwindQuad_r[5]+fUpwindQuad_r[4]-1.0*(fUpwindQuad_r[3]+fUpwindQuad_r[2])+fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[11] = 0.1767766952966368*(fUpwindQuad_r[31]+fUpwindQuad_r[30]+fUpwindQuad_r[29]+fUpwindQuad_r[28]-1.0*(fUpwindQuad_r[27]+fUpwindQuad_r[26]+fUpwindQuad_r[25]+fUpwindQuad_r[24]+fUpwindQuad_r[23]+fUpwindQuad_r[22]+fUpwindQuad_r[21]+fUpwindQuad_r[20])+fUpwindQuad_r[19]+fUpwindQuad_r[18]+fUpwindQuad_r[17]+fUpwindQuad_r[16]+fUpwindQuad_r[15]+fUpwindQuad_r[14]+fUpwindQuad_r[13]+fUpwindQuad_r[12]-1.0*(fUpwindQuad_r[11]+fUpwindQuad_r[10]+fUpwindQuad_r[9]+fUpwindQuad_r[8]+fUpwindQuad_r[7]+fUpwindQuad_r[6]+fUpwindQuad_r[5]+fUpwindQuad_r[4])+fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[12] = 0.1767766952966368*(fUpwindQuad_r[31]-1.0*fUpwindQuad_r[30]+fUpwindQuad_r[29]-1.0*fUpwindQuad_r[28]+fUpwindQuad_r[27]-1.0*fUpwindQuad_r[26]+fUpwindQuad_r[25]-1.0*fUpwindQuad_r[24]+fUpwindQuad_r[23]-1.0*fUpwindQuad_r[22]+fUpwindQuad_r[21]-1.0*fUpwindQuad_r[20]+fUpwindQuad_r[19]-1.0*fUpwindQuad_r[18]+fUpwindQuad_r[17]-1.0*(fUpwindQuad_r[16]+fUpwindQuad_r[15])+fUpwindQuad_r[14]-1.0*fUpwindQuad_r[13]+fUpwindQuad_r[12]-1.0*fUpwindQuad_r[11]+fUpwindQuad_r[10]-1.0*fUpwindQuad_r[9]+fUpwindQuad_r[8]-1.0*fUpwindQuad_r[7]+fUpwindQuad_r[6]-1.0*fUpwindQuad_r[5]+fUpwindQuad_r[4]-1.0*fUpwindQuad_r[3]+fUpwindQuad_r[2]-1.0*fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[13] = 0.1767766952966368*(fUpwindQuad_r[31]+fUpwindQuad_r[30]-1.0*(fUpwindQuad_r[29]+fUpwindQuad_r[28])+fUpwindQuad_r[27]+fUpwindQuad_r[26]-1.0*(fUpwindQuad_r[25]+fUpwindQuad_r[24])+fUpwindQuad_r[23]+fUpwindQuad_r[22]-1.0*(fUpwindQuad_r[21]+fUpwindQuad_r[20])+fUpwindQuad_r[19]+fUpwindQuad_r[18]-1.0*(fUpwindQuad_r[17]+fUpwindQuad_r[16]+fUpwindQuad_r[15]+fUpwindQuad_r[14])+fUpwindQuad_r[13]+fUpwindQuad_r[12]-1.0*(fUpwindQuad_r[11]+fUpwindQuad_r[10])+fUpwindQuad_r[9]+fUpwindQuad_r[8]-1.0*(fUpwindQuad_r[7]+fUpwindQuad_r[6])+fUpwindQuad_r[5]+fUpwindQuad_r[4]-1.0*(fUpwindQuad_r[3]+fUpwindQuad_r[2])+fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[14] = 0.1767766952966368*(fUpwindQuad_r[31]+fUpwindQuad_r[30]+fUpwindQuad_r[29]+fUpwindQuad_r[28]-1.0*(fUpwindQuad_r[27]+fUpwindQuad_r[26]+fUpwindQuad_r[25]+fUpwindQuad_r[24])+fUpwindQuad_r[23]+fUpwindQuad_r[22]+fUpwindQuad_r[21]+fUpwindQuad_r[20]-1.0*(fUpwindQuad_r[19]+fUpwindQuad_r[18]+fUpwindQuad_r[17]+fUpwindQuad_r[16]+fUpwindQuad_r[15]+fUpwindQuad_r[14]+fUpwindQuad_r[13]+fUpwindQuad_r[12])+fUpwindQuad_r[11]+fUpwindQuad_r[10]+fUpwindQuad_r[9]+fUpwindQuad_r[8]-1.0*(fUpwindQuad_r[7]+fUpwindQuad_r[6]+fUpwindQuad_r[5]+fUpwindQuad_r[4])+fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[15] = 0.1767766952966368*(fUpwindQuad_r[31]+fUpwindQuad_r[30]+fUpwindQuad_r[29]+fUpwindQuad_r[28]+fUpwindQuad_r[27]+fUpwindQuad_r[26]+fUpwindQuad_r[25]+fUpwindQuad_r[24]-1.0*(fUpwindQuad_r[23]+fUpwindQuad_r[22]+fUpwindQuad_r[21]+fUpwindQuad_r[20]+fUpwindQuad_r[19]+fUpwindQuad_r[18]+fUpwindQuad_r[17]+fUpwindQuad_r[16]+fUpwindQuad_r[15]+fUpwindQuad_r[14]+fUpwindQuad_r[13]+fUpwindQuad_r[12]+fUpwindQuad_r[11]+fUpwindQuad_r[10]+fUpwindQuad_r[9]+fUpwindQuad_r[8])+fUpwindQuad_r[7]+fUpwindQuad_r[6]+fUpwindQuad_r[5]+fUpwindQuad_r[4]+fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[16] = 0.1767766952966368*(fUpwindQuad_r[31]-1.0*(fUpwindQuad_r[30]+fUpwindQuad_r[29])+fUpwindQuad_r[28]-1.0*fUpwindQuad_r[27]+fUpwindQuad_r[26]+fUpwindQuad_r[25]-1.0*fUpwindQuad_r[24]+fUpwindQuad_r[23]-1.0*(fUpwindQuad_r[22]+fUpwindQuad_r[21])+fUpwindQuad_r[20]-1.0*fUpwindQuad_r[19]+fUpwindQuad_r[18]+fUpwindQuad_r[17]-1.0*fUpwindQuad_r[16]+fUpwindQuad_r[15]-1.0*(fUpwindQuad_r[14]+fUpwindQuad_r[13])+fUpwindQuad_r[12]-1.0*fUpwindQuad_r[11]+fUpwindQuad_r[10]+fUpwindQuad_r[9]-1.0*fUpwindQuad_r[8]+fUpwindQuad_r[7]-1.0*(fUpwindQuad_r[6]+fUpwindQuad_r[5])+fUpwindQuad_r[4]-1.0*fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]-1.0*fUpwindQuad_r[0]); 
  fUpwind_r[17] = 0.1767766952966368*(fUpwindQuad_r[31]-1.0*(fUpwindQuad_r[30]+fUpwindQuad_r[29])+fUpwindQuad_r[28]+fUpwindQuad_r[27]-1.0*(fUpwindQuad_r[26]+fUpwindQuad_r[25])+fUpwindQuad_r[24]-1.0*fUpwindQuad_r[23]+fUpwindQuad_r[22]+fUpwindQuad_r[21]-1.0*(fUpwindQuad_r[20]+fUpwindQuad_r[19])+fUpwindQuad_r[18]+fUpwindQuad_r[17]-1.0*fUpwindQuad_r[16]+fUpwindQuad_r[15]-1.0*(fUpwindQuad_r[14]+fUpwindQuad_r[13])+fUpwindQuad_r[12]+fUpwindQuad_r[11]-1.0*(fUpwindQuad_r[10]+fUpwindQuad_r[9])+fUpwindQuad_r[8]-1.0*fUpwindQuad_r[7]+fUpwindQuad_r[6]+fUpwindQuad_r[5]-1.0*(fUpwindQuad_r[4]+fUpwindQuad_r[3])+fUpwindQuad_r[2]+fUpwindQuad_r[1]-1.0*fUpwindQuad_r[0]); 
  fUpwind_r[18] = 0.1767766952966368*(fUpwindQuad_r[31]-1.0*fUpwindQuad_r[30]+fUpwindQuad_r[29]-1.0*(fUpwindQuad_r[28]+fUpwindQuad_r[27])+fUpwindQuad_r[26]-1.0*fUpwindQuad_r[25]+fUpwindQuad_r[24]-1.0*fUpwindQuad_r[23]+fUpwindQuad_r[22]-1.0*fUpwindQuad_r[21]+fUpwindQuad_r[20]+fUpwindQuad_r[19]-1.0*fUpwindQuad_r[18]+fUpwindQuad_r[17]-1.0*fUpwindQuad_r[16]+fUpwindQuad_r[15]-1.0*fUpwindQuad_r[14]+fUpwindQuad_r[13]-1.0*(fUpwindQuad_r[12]+fUpwindQuad_r[11])+fUpwindQuad_r[10]-1.0*fUpwindQuad_r[9]+fUpwindQuad_r[8]-1.0*fUpwindQuad_r[7]+fUpwindQuad_r[6]-1.0*fUpwindQuad_r[5]+fUpwindQuad_r[4]+fUpwindQuad_r[3]-1.0*fUpwindQuad_r[2]+fUpwindQuad_r[1]-1.0*fUpwindQuad_r[0]); 
  fUpwind_r[19] = 0.1767766952966368*(fUpwindQuad_r[31]+fUpwindQuad_r[30]-1.0*(fUpwindQuad_r[29]+fUpwindQuad_r[28]+fUpwindQuad_r[27]+fUpwindQuad_r[26])+fUpwindQuad_r[25]+fUpwindQuad_r[24]-1.0*(fUpwindQuad_r[23]+fUpwindQuad_r[22])+fUpwindQuad_r[21]+fUpwindQuad_r[20]+fUpwindQuad_r[19]+fUpwindQuad_r[18]-1.0*(fUpwindQuad_r[17]+fUpwindQuad_r[16])+fUpwindQuad_r[15]+fUpwindQuad_r[14]-1.0*(fUpwindQuad_r[13]+fUpwindQuad_r[12]+fUpwindQuad_r[11]+fUpwindQuad_r[10])+fUpwindQuad_r[9]+fUpwindQuad_r[8]-1.0*(fUpwindQuad_r[7]+fUpwindQuad_r[6])+fUpwindQuad_r[5]+fUpwindQuad_r[4]+fUpwindQuad_r[3]+fUpwindQuad_r[2]-1.0*(fUpwindQuad_r[1]+fUpwindQuad_r[0])); 
  fUpwind_r[20] = 0.1767766952966368*(fUpwindQuad_r[31]-1.0*(fUpwindQuad_r[30]+fUpwindQuad_r[29])+fUpwindQuad_r[28]+fUpwindQuad_r[27]-1.0*(fUpwindQuad_r[26]+fUpwindQuad_r[25])+fUpwindQuad_r[24]+fUpwindQuad_r[23]-1.0*(fUpwindQuad_r[22]+fUpwindQuad_r[21])+fUpwindQuad_r[20]+fUpwindQuad_r[19]-1.0*(fUpwindQuad_r[18]+fUpwindQuad_r[17])+fUpwindQuad_r[16]-1.0*fUpwindQuad_r[15]+fUpwindQuad_r[14]+fUpwindQuad_r[13]-1.0*(fUpwindQuad_r[12]+fUpwindQuad_r[11])+fUpwindQuad_r[10]+fUpwindQuad_r[9]-1.0*(fUpwindQuad_r[8]+fUpwindQuad_r[7])+fUpwindQuad_r[6]+fUpwindQuad_r[5]-1.0*(fUpwindQuad_r[4]+fUpwindQuad_r[3])+fUpwindQuad_r[2]+fUpwindQuad_r[1]-1.0*fUpwindQuad_r[0]); 
  fUpwind_r[21] = 0.1767766952966368*(fUpwindQuad_r[31]-1.0*fUpwindQuad_r[30]+fUpwindQuad_r[29]-1.0*(fUpwindQuad_r[28]+fUpwindQuad_r[27])+fUpwindQuad_r[26]-1.0*fUpwindQuad_r[25]+fUpwindQuad_r[24]+fUpwindQuad_r[23]-1.0*fUpwindQuad_r[22]+fUpwindQuad_r[21]-1.0*(fUpwindQuad_r[20]+fUpwindQuad_r[19])+fUpwindQuad_r[18]-1.0*fUpwindQuad_r[17]+fUpwindQuad_r[16]-1.0*fUpwindQuad_r[15]+fUpwindQuad_r[14]-1.0*fUpwindQuad_r[13]+fUpwindQuad_r[12]+fUpwindQuad_r[11]-1.0*fUpwindQuad_r[10]+fUpwindQuad_r[9]-1.0*(fUpwindQuad_r[8]+fUpwindQuad_r[7])+fUpwindQuad_r[6]-1.0*fUpwindQuad_r[5]+fUpwindQuad_r[4]+fUpwindQuad_r[3]-1.0*fUpwindQuad_r[2]+fUpwindQuad_r[1]-1.0*fUpwindQuad_r[0]); 
  fUpwind_r[22] = 0.1767766952966368*(fUpwindQuad_r[31]+fUpwindQuad_r[30]-1.0*(fUpwindQuad_r[29]+fUpwindQuad_r[28]+fUpwindQuad_r[27]+fUpwindQuad_r[26])+fUpwindQuad_r[25]+fUpwindQuad_r[24]+fUpwindQuad_r[23]+fUpwindQuad_r[22]-1.0*(fUpwindQuad_r[21]+fUpwindQuad_r[20]+fUpwindQuad_r[19]+fUpwindQuad_r[18])+fUpwindQuad_r[17]+fUpwindQuad_r[16]-1.0*(fUpwindQuad_r[15]+fUpwindQuad_r[14])+fUpwindQuad_r[13]+fUpwindQuad_r[12]+fUpwindQuad_r[11]+fUpwindQuad_r[10]-1.0*(fUpwindQuad_r[9]+fUpwindQuad_r[8]+fUpwindQuad_r[7]+fUpwindQuad_r[6])+fUpwindQuad_r[5]+fUpwindQuad_r[4]+fUpwindQuad_r[3]+fUpwindQuad_r[2]-1.0*(fUpwindQuad_r[1]+fUpwindQuad_r[0])); 
  fUpwind_r[23] = 0.1767766952966368*(fUpwindQuad_r[31]-1.0*fUpwindQuad_r[30]+fUpwindQuad_r[29]-1.0*fUpwindQuad_r[28]+fUpwindQuad_r[27]-1.0*fUpwindQuad_r[26]+fUpwindQuad_r[25]-1.0*(fUpwindQuad_r[24]+fUpwindQuad_r[23])+fUpwindQuad_r[22]-1.0*fUpwindQuad_r[21]+fUpwindQuad_r[20]-1.0*fUpwindQuad_r[19]+fUpwindQuad_r[18]-1.0*fUpwindQuad_r[17]+fUpwindQuad_r[16]-1.0*fUpwindQuad_r[15]+fUpwindQuad_r[14]-1.0*fUpwindQuad_r[13]+fUpwindQuad_r[12]-1.0*fUpwindQuad_r[11]+fUpwindQuad_r[10]-1.0*fUpwindQuad_r[9]+fUpwindQuad_r[8]+fUpwindQuad_r[7]-1.0*fUpwindQuad_r[6]+fUpwindQuad_r[5]-1.0*fUpwindQuad_r[4]+fUpwindQuad_r[3]-1.0*fUpwindQuad_r[2]+fUpwindQuad_r[1]-1.0*fUpwindQuad_r[0]); 
  fUpwind_r[24] = 0.1767766952966368*(fUpwindQuad_r[31]+fUpwindQuad_r[30]-1.0*(fUpwindQuad_r[29]+fUpwindQuad_r[28])+fUpwindQuad_r[27]+fUpwindQuad_r[26]-1.0*(fUpwindQuad_r[25]+fUpwindQuad_r[24]+fUpwindQuad_r[23]+fUpwindQuad_r[22])+fUpwindQuad_r[21]+fUpwindQuad_r[20]-1.0*(fUpwindQuad_r[19]+fUpwindQuad_r[18])+fUpwindQuad_r[17]+fUpwindQuad_r[16]-1.0*(fUpwindQuad_r[15]+fUpwindQuad_r[14])+fUpwindQuad_r[13]+fUpwindQuad_r[12]-1.0*(fUpwindQuad_r[11]+fUpwindQuad_r[10])+fUpwindQuad_r[9]+fUpwindQuad_r[8]+fUpwindQuad_r[7]+fUpwindQuad_r[6]-1.0*(fUpwindQuad_r[5]+fUpwindQuad_r[4])+fUpwindQuad_r[3]+fUpwindQuad_r[2]-1.0*(fUpwindQuad_r[1]+fUpwindQuad_r[0])); 
  fUpwind_r[25] = 0.1767766952966368*(fUpwindQuad_r[31]+fUpwindQuad_r[30]+fUpwindQuad_r[29]+fUpwindQuad_r[28]-1.0*(fUpwindQuad_r[27]+fUpwindQuad_r[26]+fUpwindQuad_r[25]+fUpwindQuad_r[24]+fUpwindQuad_r[23]+fUpwindQuad_r[22]+fUpwindQuad_r[21]+fUpwindQuad_r[20])+fUpwindQuad_r[19]+fUpwindQuad_r[18]+fUpwindQuad_r[17]+fUpwindQuad_r[16]-1.0*(fUpwindQuad_r[15]+fUpwindQuad_r[14]+fUpwindQuad_r[13]+fUpwindQuad_r[12])+fUpwindQuad_r[11]+fUpwindQuad_r[10]+fUpwindQuad_r[9]+fUpwindQuad_r[8]+fUpwindQuad_r[7]+fUpwindQuad_r[6]+fUpwindQuad_r[5]+fUpwindQuad_r[4]-1.0*(fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]+fUpwindQuad_r[0])); 
  fUpwind_r[26] = 0.1767766952966368*(fUpwindQuad_r[31]-1.0*(fUpwindQuad_r[30]+fUpwindQuad_r[29])+fUpwindQuad_r[28]-1.0*fUpwindQuad_r[27]+fUpwindQuad_r[26]+fUpwindQuad_r[25]-1.0*(fUpwindQuad_r[24]+fUpwindQuad_r[23])+fUpwindQuad_r[22]+fUpwindQuad_r[21]-1.0*fUpwindQuad_r[20]+fUpwindQuad_r[19]-1.0*(fUpwindQuad_r[18]+fUpwindQuad_r[17])+fUpwindQuad_r[16]+fUpwindQuad_r[15]-1.0*(fUpwindQuad_r[14]+fUpwindQuad_r[13])+fUpwindQuad_r[12]-1.0*fUpwindQuad_r[11]+fUpwindQuad_r[10]+fUpwindQuad_r[9]-1.0*(fUpwindQuad_r[8]+fUpwindQuad_r[7])+fUpwindQuad_r[6]+fUpwindQuad_r[5]-1.0*fUpwindQuad_r[4]+fUpwindQuad_r[3]-1.0*(fUpwindQuad_r[2]+fUpwindQuad_r[1])+fUpwindQuad_r[0]); 
  fUpwind_r[27] = 0.1767766952966368*(fUpwindQuad_r[31]-1.0*(fUpwindQuad_r[30]+fUpwindQuad_r[29])+fUpwindQuad_r[28]-1.0*fUpwindQuad_r[27]+fUpwindQuad_r[26]+fUpwindQuad_r[25]-1.0*fUpwindQuad_r[24]+fUpwindQuad_r[23]-1.0*(fUpwindQuad_r[22]+fUpwindQuad_r[21])+fUpwindQuad_r[20]-1.0*fUpwindQuad_r[19]+fUpwindQuad_r[18]+fUpwindQuad_r[17]-1.0*(fUpwindQuad_r[16]+fUpwindQuad_r[15])+fUpwindQuad_r[14]+fUpwindQuad_r[13]-1.0*fUpwindQuad_r[12]+fUpwindQuad_r[11]-1.0*(fUpwindQuad_r[10]+fUpwindQuad_r[9])+fUpwindQuad_r[8]-1.0*fUpwindQuad_r[7]+fUpwindQuad_r[6]+fUpwindQuad_r[5]-1.0*fUpwindQuad_r[4]+fUpwindQuad_r[3]-1.0*(fUpwindQuad_r[2]+fUpwindQuad_r[1])+fUpwindQuad_r[0]); 
  fUpwind_r[28] = 0.1767766952966368*(fUpwindQuad_r[31]-1.0*(fUpwindQuad_r[30]+fUpwindQuad_r[29])+fUpwindQuad_r[28]+fUpwindQuad_r[27]-1.0*(fUpwindQuad_r[26]+fUpwindQuad_r[25])+fUpwindQuad_r[24]-1.0*fUpwindQuad_r[23]+fUpwindQuad_r[22]+fUpwindQuad_r[21]-1.0*(fUpwindQuad_r[20]+fUpwindQuad_r[19])+fUpwindQuad_r[18]+fUpwindQuad_r[17]-1.0*(fUpwindQuad_r[16]+fUpwindQuad_r[15])+fUpwindQuad_r[14]+fUpwindQuad_r[13]-1.0*(fUpwindQuad_r[12]+fUpwindQuad_r[11])+fUpwindQuad_r[10]+fUpwindQuad_r[9]-1.0*fUpwindQuad_r[8]+fUpwindQuad_r[7]-1.0*(fUpwindQuad_r[6]+fUpwindQuad_r[5])+fUpwindQuad_r[4]+fUpwindQuad_r[3]-1.0*(fUpwindQuad_r[2]+fUpwindQuad_r[1])+fUpwindQuad_r[0]); 
  fUpwind_r[29] = 0.1767766952966368*(fUpwindQuad_r[31]-1.0*fUpwindQuad_r[30]+fUpwindQuad_r[29]-1.0*(fUpwindQuad_r[28]+fUpwindQuad_r[27])+fUpwindQuad_r[26]-1.0*fUpwindQuad_r[25]+fUpwindQuad_r[24]-1.0*fUpwindQuad_r[23]+fUpwindQuad_r[22]-1.0*fUpwindQuad_r[21]+fUpwindQuad_r[20]+fUpwindQuad_r[19]-1.0*fUpwindQuad_r[18]+fUpwindQuad_r[17]-1.0*(fUpwindQuad_r[16]+fUpwindQuad_r[15])+fUpwindQuad_r[14]-1.0*fUpwindQuad_r[13]+fUpwindQuad_r[12]+fUpwindQuad_r[11]-1.0*fUpwindQuad_r[10]+fUpwindQuad_r[9]-1.0*fUpwindQuad_r[8]+fUpwindQuad_r[7]-1.0*fUpwindQuad_r[6]+fUpwindQuad_r[5]-1.0*(fUpwindQuad_r[4]+fUpwindQuad_r[3])+fUpwindQuad_r[2]-1.0*fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[30] = 0.1767766952966368*(fUpwindQuad_r[31]+fUpwindQuad_r[30]-1.0*(fUpwindQuad_r[29]+fUpwindQuad_r[28]+fUpwindQuad_r[27]+fUpwindQuad_r[26])+fUpwindQuad_r[25]+fUpwindQuad_r[24]-1.0*(fUpwindQuad_r[23]+fUpwindQuad_r[22])+fUpwindQuad_r[21]+fUpwindQuad_r[20]+fUpwindQuad_r[19]+fUpwindQuad_r[18]-1.0*(fUpwindQuad_r[17]+fUpwindQuad_r[16]+fUpwindQuad_r[15]+fUpwindQuad_r[14])+fUpwindQuad_r[13]+fUpwindQuad_r[12]+fUpwindQuad_r[11]+fUpwindQuad_r[10]-1.0*(fUpwindQuad_r[9]+fUpwindQuad_r[8])+fUpwindQuad_r[7]+fUpwindQuad_r[6]-1.0*(fUpwindQuad_r[5]+fUpwindQuad_r[4]+fUpwindQuad_r[3]+fUpwindQuad_r[2])+fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[31] = 0.1767766952966368*(fUpwindQuad_r[31]-1.0*(fUpwindQuad_r[30]+fUpwindQuad_r[29])+fUpwindQuad_r[28]-1.0*fUpwindQuad_r[27]+fUpwindQuad_r[26]+fUpwindQuad_r[25]-1.0*(fUpwindQuad_r[24]+fUpwindQuad_r[23])+fUpwindQuad_r[22]+fUpwindQuad_r[21]-1.0*fUpwindQuad_r[20]+fUpwindQuad_r[19]-1.0*(fUpwindQuad_r[18]+fUpwindQuad_r[17])+fUpwindQuad_r[16]-1.0*fUpwindQuad_r[15]+fUpwindQuad_r[14]+fUpwindQuad_r[13]-1.0*fUpwindQuad_r[12]+fUpwindQuad_r[11]-1.0*(fUpwindQuad_r[10]+fUpwindQuad_r[9])+fUpwindQuad_r[8]+fUpwindQuad_r[7]-1.0*(fUpwindQuad_r[6]+fUpwindQuad_r[5])+fUpwindQuad_r[4]-1.0*fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]-1.0*fUpwindQuad_r[0]); 

  Gdrag_l[0] = 0.1767766952966368*alphaDrSurf_l[16]*fUpwind_l[16]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[8]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[7]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[6]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[3]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[2]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[1]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[0]; 
  Gdrag_l[1] = 0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[16]+0.1767766952966368*fUpwind_l[8]*alphaDrSurf_l[16]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[7]+0.1767766952966368*fUpwind_l[3]*alphaDrSurf_l[7]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[6]+0.1767766952966368*fUpwind_l[2]*alphaDrSurf_l[6]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[1]+0.1767766952966368*fUpwind_l[0]*alphaDrSurf_l[1]; 
  Gdrag_l[2] = 0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[16]+0.1767766952966368*fUpwind_l[7]*alphaDrSurf_l[16]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[8]+0.1767766952966368*fUpwind_l[3]*alphaDrSurf_l[8]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[6]+0.1767766952966368*fUpwind_l[1]*alphaDrSurf_l[6]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[2]+0.1767766952966368*fUpwind_l[0]*alphaDrSurf_l[2]; 
  Gdrag_l[3] = 0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[16]+0.1767766952966368*fUpwind_l[6]*alphaDrSurf_l[16]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[8]+0.1767766952966368*fUpwind_l[2]*alphaDrSurf_l[8]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[7]+0.1767766952966368*fUpwind_l[1]*alphaDrSurf_l[7]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[3]+0.1767766952966368*fUpwind_l[0]*alphaDrSurf_l[3]; 
  Gdrag_l[4] = 0.1767766952966368*alphaDrSurf_l[16]*fUpwind_l[26]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[19]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[18]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[17]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[11]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[10]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[9]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[4]; 
  Gdrag_l[5] = 0.1767766952966368*alphaDrSurf_l[16]*fUpwind_l[27]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[22]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[21]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[20]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[14]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[13]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[12]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[5]; 
  Gdrag_l[6] = 0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[16]+0.1767766952966368*fUpwind_l[3]*alphaDrSurf_l[16]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[8]+0.1767766952966368*fUpwind_l[7]*alphaDrSurf_l[8]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[6]+0.1767766952966368*fUpwind_l[0]*alphaDrSurf_l[6]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[2]+0.1767766952966368*fUpwind_l[1]*alphaDrSurf_l[2]; 
  Gdrag_l[7] = 0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[16]+0.1767766952966368*fUpwind_l[2]*alphaDrSurf_l[16]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[8]+0.1767766952966368*fUpwind_l[6]*alphaDrSurf_l[8]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[7]+0.1767766952966368*fUpwind_l[0]*alphaDrSurf_l[7]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[3]+0.1767766952966368*fUpwind_l[1]*alphaDrSurf_l[3]; 
  Gdrag_l[8] = 0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[16]+0.1767766952966368*fUpwind_l[1]*alphaDrSurf_l[16]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[8]+0.1767766952966368*fUpwind_l[0]*alphaDrSurf_l[8]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[7]+0.1767766952966368*fUpwind_l[6]*alphaDrSurf_l[7]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[3]+0.1767766952966368*fUpwind_l[2]*alphaDrSurf_l[3]; 
  Gdrag_l[9] = 0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[26]+0.1767766952966368*alphaDrSurf_l[16]*fUpwind_l[19]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[18]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[17]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[11]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[10]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[9]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[4]; 
  Gdrag_l[10] = 0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[26]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[19]+0.1767766952966368*alphaDrSurf_l[16]*fUpwind_l[18]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[17]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[11]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[10]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[9]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[4]; 
  Gdrag_l[11] = 0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[26]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[19]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[18]+0.1767766952966368*alphaDrSurf_l[16]*fUpwind_l[17]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[11]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[10]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[9]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[4]; 
  Gdrag_l[12] = 0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[27]+0.1767766952966368*alphaDrSurf_l[16]*fUpwind_l[22]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[21]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[20]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[14]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[13]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[12]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[5]; 
  Gdrag_l[13] = 0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[27]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[22]+0.1767766952966368*alphaDrSurf_l[16]*fUpwind_l[21]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[20]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[14]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[13]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[12]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[5]; 
  Gdrag_l[14] = 0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[27]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[22]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[21]+0.1767766952966368*alphaDrSurf_l[16]*fUpwind_l[20]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[14]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[13]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[12]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[5]; 
  Gdrag_l[15] = 0.1767766952966368*alphaDrSurf_l[16]*fUpwind_l[31]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[30]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[29]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[28]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[25]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[24]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[23]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[15]; 
  Gdrag_l[16] = 0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[16]+0.1767766952966368*fUpwind_l[0]*alphaDrSurf_l[16]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[8]+0.1767766952966368*fUpwind_l[1]*alphaDrSurf_l[8]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[7]+0.1767766952966368*fUpwind_l[2]*alphaDrSurf_l[7]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[6]+0.1767766952966368*fUpwind_l[3]*alphaDrSurf_l[6]; 
  Gdrag_l[17] = 0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[26]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[19]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[18]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[17]+0.1767766952966368*fUpwind_l[11]*alphaDrSurf_l[16]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[10]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[9]+0.1767766952966368*fUpwind_l[4]*alphaDrSurf_l[6]; 
  Gdrag_l[18] = 0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[26]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[19]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[18]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[17]+0.1767766952966368*fUpwind_l[10]*alphaDrSurf_l[16]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[11]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[9]+0.1767766952966368*fUpwind_l[4]*alphaDrSurf_l[7]; 
  Gdrag_l[19] = 0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[26]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[19]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[18]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[17]+0.1767766952966368*fUpwind_l[9]*alphaDrSurf_l[16]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[11]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[10]+0.1767766952966368*fUpwind_l[4]*alphaDrSurf_l[8]; 
  Gdrag_l[20] = 0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[27]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[22]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[21]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[20]+0.1767766952966368*fUpwind_l[14]*alphaDrSurf_l[16]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[13]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[12]+0.1767766952966368*fUpwind_l[5]*alphaDrSurf_l[6]; 
  Gdrag_l[21] = 0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[27]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[22]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[21]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[20]+0.1767766952966368*fUpwind_l[13]*alphaDrSurf_l[16]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[14]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[12]+0.1767766952966368*fUpwind_l[5]*alphaDrSurf_l[7]; 
  Gdrag_l[22] = 0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[27]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[22]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[21]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[20]+0.1767766952966368*fUpwind_l[12]*alphaDrSurf_l[16]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[14]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[13]+0.1767766952966368*fUpwind_l[5]*alphaDrSurf_l[8]; 
  Gdrag_l[23] = 0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[31]+0.1767766952966368*alphaDrSurf_l[16]*fUpwind_l[30]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[29]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[28]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[25]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[24]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[23]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[15]; 
  Gdrag_l[24] = 0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[31]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[30]+0.1767766952966368*alphaDrSurf_l[16]*fUpwind_l[29]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[28]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[25]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[24]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[23]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[15]; 
  Gdrag_l[25] = 0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[31]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[30]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[29]+0.1767766952966368*alphaDrSurf_l[16]*fUpwind_l[28]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[25]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[24]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[23]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[15]; 
  Gdrag_l[26] = 0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[26]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[19]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[18]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[17]+0.1767766952966368*fUpwind_l[4]*alphaDrSurf_l[16]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[11]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[10]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[9]; 
  Gdrag_l[27] = 0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[27]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[22]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[21]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[20]+0.1767766952966368*fUpwind_l[5]*alphaDrSurf_l[16]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[14]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[13]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[12]; 
  Gdrag_l[28] = 0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[31]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[30]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[29]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[28]+0.1767766952966368*alphaDrSurf_l[16]*fUpwind_l[25]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[24]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[23]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[15]; 
  Gdrag_l[29] = 0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[31]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[30]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[29]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[28]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[25]+0.1767766952966368*alphaDrSurf_l[16]*fUpwind_l[24]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[23]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[15]; 
  Gdrag_l[30] = 0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[31]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[30]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[29]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[28]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[25]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[24]+0.1767766952966368*alphaDrSurf_l[16]*fUpwind_l[23]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[15]; 
  Gdrag_l[31] = 0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[31]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[30]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[29]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[28]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[25]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[24]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[23]+0.1767766952966368*fUpwind_l[15]*alphaDrSurf_l[16]; 

  Gdrag_r[0] = 0.1767766952966368*alphaDrSurf_r[16]*fUpwind_r[16]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[8]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[7]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[6]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[3]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[2]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[1]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[0]; 
  Gdrag_r[1] = 0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[16]+0.1767766952966368*fUpwind_r[8]*alphaDrSurf_r[16]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[7]+0.1767766952966368*fUpwind_r[3]*alphaDrSurf_r[7]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[6]+0.1767766952966368*fUpwind_r[2]*alphaDrSurf_r[6]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[1]+0.1767766952966368*fUpwind_r[0]*alphaDrSurf_r[1]; 
  Gdrag_r[2] = 0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[16]+0.1767766952966368*fUpwind_r[7]*alphaDrSurf_r[16]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[8]+0.1767766952966368*fUpwind_r[3]*alphaDrSurf_r[8]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[6]+0.1767766952966368*fUpwind_r[1]*alphaDrSurf_r[6]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[2]+0.1767766952966368*fUpwind_r[0]*alphaDrSurf_r[2]; 
  Gdrag_r[3] = 0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[16]+0.1767766952966368*fUpwind_r[6]*alphaDrSurf_r[16]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[8]+0.1767766952966368*fUpwind_r[2]*alphaDrSurf_r[8]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[7]+0.1767766952966368*fUpwind_r[1]*alphaDrSurf_r[7]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[3]+0.1767766952966368*fUpwind_r[0]*alphaDrSurf_r[3]; 
  Gdrag_r[4] = 0.1767766952966368*alphaDrSurf_r[16]*fUpwind_r[26]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[19]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[18]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[17]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[11]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[10]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[9]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[4]; 
  Gdrag_r[5] = 0.1767766952966368*alphaDrSurf_r[16]*fUpwind_r[27]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[22]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[21]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[20]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[14]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[13]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[12]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[5]; 
  Gdrag_r[6] = 0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[16]+0.1767766952966368*fUpwind_r[3]*alphaDrSurf_r[16]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[8]+0.1767766952966368*fUpwind_r[7]*alphaDrSurf_r[8]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[6]+0.1767766952966368*fUpwind_r[0]*alphaDrSurf_r[6]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[2]+0.1767766952966368*fUpwind_r[1]*alphaDrSurf_r[2]; 
  Gdrag_r[7] = 0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[16]+0.1767766952966368*fUpwind_r[2]*alphaDrSurf_r[16]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[8]+0.1767766952966368*fUpwind_r[6]*alphaDrSurf_r[8]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[7]+0.1767766952966368*fUpwind_r[0]*alphaDrSurf_r[7]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[3]+0.1767766952966368*fUpwind_r[1]*alphaDrSurf_r[3]; 
  Gdrag_r[8] = 0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[16]+0.1767766952966368*fUpwind_r[1]*alphaDrSurf_r[16]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[8]+0.1767766952966368*fUpwind_r[0]*alphaDrSurf_r[8]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[7]+0.1767766952966368*fUpwind_r[6]*alphaDrSurf_r[7]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[3]+0.1767766952966368*fUpwind_r[2]*alphaDrSurf_r[3]; 
  Gdrag_r[9] = 0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[26]+0.1767766952966368*alphaDrSurf_r[16]*fUpwind_r[19]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[18]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[17]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[11]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[10]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[9]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[4]; 
  Gdrag_r[10] = 0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[26]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[19]+0.1767766952966368*alphaDrSurf_r[16]*fUpwind_r[18]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[17]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[11]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[10]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[9]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[4]; 
  Gdrag_r[11] = 0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[26]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[19]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[18]+0.1767766952966368*alphaDrSurf_r[16]*fUpwind_r[17]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[11]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[10]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[9]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[4]; 
  Gdrag_r[12] = 0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[27]+0.1767766952966368*alphaDrSurf_r[16]*fUpwind_r[22]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[21]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[20]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[14]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[13]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[12]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[5]; 
  Gdrag_r[13] = 0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[27]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[22]+0.1767766952966368*alphaDrSurf_r[16]*fUpwind_r[21]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[20]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[14]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[13]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[12]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[5]; 
  Gdrag_r[14] = 0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[27]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[22]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[21]+0.1767766952966368*alphaDrSurf_r[16]*fUpwind_r[20]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[14]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[13]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[12]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[5]; 
  Gdrag_r[15] = 0.1767766952966368*alphaDrSurf_r[16]*fUpwind_r[31]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[30]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[29]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[28]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[25]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[24]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[23]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[15]; 
  Gdrag_r[16] = 0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[16]+0.1767766952966368*fUpwind_r[0]*alphaDrSurf_r[16]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[8]+0.1767766952966368*fUpwind_r[1]*alphaDrSurf_r[8]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[7]+0.1767766952966368*fUpwind_r[2]*alphaDrSurf_r[7]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[6]+0.1767766952966368*fUpwind_r[3]*alphaDrSurf_r[6]; 
  Gdrag_r[17] = 0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[26]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[19]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[18]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[17]+0.1767766952966368*fUpwind_r[11]*alphaDrSurf_r[16]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[10]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[9]+0.1767766952966368*fUpwind_r[4]*alphaDrSurf_r[6]; 
  Gdrag_r[18] = 0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[26]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[19]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[18]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[17]+0.1767766952966368*fUpwind_r[10]*alphaDrSurf_r[16]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[11]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[9]+0.1767766952966368*fUpwind_r[4]*alphaDrSurf_r[7]; 
  Gdrag_r[19] = 0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[26]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[19]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[18]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[17]+0.1767766952966368*fUpwind_r[9]*alphaDrSurf_r[16]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[11]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[10]+0.1767766952966368*fUpwind_r[4]*alphaDrSurf_r[8]; 
  Gdrag_r[20] = 0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[27]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[22]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[21]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[20]+0.1767766952966368*fUpwind_r[14]*alphaDrSurf_r[16]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[13]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[12]+0.1767766952966368*fUpwind_r[5]*alphaDrSurf_r[6]; 
  Gdrag_r[21] = 0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[27]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[22]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[21]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[20]+0.1767766952966368*fUpwind_r[13]*alphaDrSurf_r[16]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[14]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[12]+0.1767766952966368*fUpwind_r[5]*alphaDrSurf_r[7]; 
  Gdrag_r[22] = 0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[27]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[22]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[21]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[20]+0.1767766952966368*fUpwind_r[12]*alphaDrSurf_r[16]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[14]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[13]+0.1767766952966368*fUpwind_r[5]*alphaDrSurf_r[8]; 
  Gdrag_r[23] = 0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[31]+0.1767766952966368*alphaDrSurf_r[16]*fUpwind_r[30]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[29]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[28]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[25]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[24]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[23]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[15]; 
  Gdrag_r[24] = 0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[31]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[30]+0.1767766952966368*alphaDrSurf_r[16]*fUpwind_r[29]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[28]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[25]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[24]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[23]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[15]; 
  Gdrag_r[25] = 0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[31]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[30]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[29]+0.1767766952966368*alphaDrSurf_r[16]*fUpwind_r[28]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[25]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[24]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[23]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[15]; 
  Gdrag_r[26] = 0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[26]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[19]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[18]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[17]+0.1767766952966368*fUpwind_r[4]*alphaDrSurf_r[16]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[11]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[10]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[9]; 
  Gdrag_r[27] = 0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[27]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[22]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[21]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[20]+0.1767766952966368*fUpwind_r[5]*alphaDrSurf_r[16]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[14]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[13]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[12]; 
  Gdrag_r[28] = 0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[31]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[30]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[29]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[28]+0.1767766952966368*alphaDrSurf_r[16]*fUpwind_r[25]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[24]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[23]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[15]; 
  Gdrag_r[29] = 0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[31]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[30]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[29]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[28]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[25]+0.1767766952966368*alphaDrSurf_r[16]*fUpwind_r[24]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[23]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[15]; 
  Gdrag_r[30] = 0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[31]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[30]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[29]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[28]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[25]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[24]+0.1767766952966368*alphaDrSurf_r[16]*fUpwind_r[23]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[15]; 
  Gdrag_r[31] = 0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[31]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[30]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[29]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[28]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[25]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[24]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[23]+0.1767766952966368*fUpwind_r[15]*alphaDrSurf_r[16]; 

  out[0] += 0.7071067811865475*Gdrag_r[0]*rdv2-0.7071067811865475*Gdrag_l[0]*rdv2; 
  out[1] += 0.7071067811865475*Gdrag_r[1]*rdv2-0.7071067811865475*Gdrag_l[1]*rdv2; 
  out[2] += 0.7071067811865475*Gdrag_r[2]*rdv2-0.7071067811865475*Gdrag_l[2]*rdv2; 
  out[3] += 0.7071067811865475*Gdrag_r[3]*rdv2-0.7071067811865475*Gdrag_l[3]*rdv2; 
  out[4] += 0.7071067811865475*Gdrag_r[4]*rdv2-0.7071067811865475*Gdrag_l[4]*rdv2; 
  out[5] += 1.224744871391589*Gdrag_r[0]*rdv2+1.224744871391589*Gdrag_l[0]*rdv2; 
  out[6] += 0.7071067811865475*Gdrag_r[5]*rdv2-0.7071067811865475*Gdrag_l[5]*rdv2; 
  out[7] += 0.7071067811865475*Gdrag_r[6]*rdv2-0.7071067811865475*Gdrag_l[6]*rdv2; 
  out[8] += 0.7071067811865475*Gdrag_r[7]*rdv2-0.7071067811865475*Gdrag_l[7]*rdv2; 
  out[9] += 0.7071067811865475*Gdrag_r[8]*rdv2-0.7071067811865475*Gdrag_l[8]*rdv2; 
  out[10] += 0.7071067811865475*Gdrag_r[9]*rdv2-0.7071067811865475*Gdrag_l[9]*rdv2; 
  out[11] += 0.7071067811865475*Gdrag_r[10]*rdv2-0.7071067811865475*Gdrag_l[10]*rdv2; 
  out[12] += 0.7071067811865475*Gdrag_r[11]*rdv2-0.7071067811865475*Gdrag_l[11]*rdv2; 
  out[13] += 1.224744871391589*Gdrag_r[1]*rdv2+1.224744871391589*Gdrag_l[1]*rdv2; 
  out[14] += 1.224744871391589*Gdrag_r[2]*rdv2+1.224744871391589*Gdrag_l[2]*rdv2; 
  out[15] += 1.224744871391589*Gdrag_r[3]*rdv2+1.224744871391589*Gdrag_l[3]*rdv2; 
  out[16] += 1.224744871391589*Gdrag_r[4]*rdv2+1.224744871391589*Gdrag_l[4]*rdv2; 
  out[17] += 0.7071067811865475*Gdrag_r[12]*rdv2-0.7071067811865475*Gdrag_l[12]*rdv2; 
  out[18] += 0.7071067811865475*Gdrag_r[13]*rdv2-0.7071067811865475*Gdrag_l[13]*rdv2; 
  out[19] += 0.7071067811865475*Gdrag_r[14]*rdv2-0.7071067811865475*Gdrag_l[14]*rdv2; 
  out[20] += 0.7071067811865475*Gdrag_r[15]*rdv2-0.7071067811865475*Gdrag_l[15]*rdv2; 
  out[21] += 1.224744871391589*Gdrag_r[5]*rdv2+1.224744871391589*Gdrag_l[5]*rdv2; 
  out[22] += 0.7071067811865475*Gdrag_r[16]*rdv2-0.7071067811865475*Gdrag_l[16]*rdv2; 
  out[23] += 0.7071067811865475*Gdrag_r[17]*rdv2-0.7071067811865475*Gdrag_l[17]*rdv2; 
  out[24] += 0.7071067811865475*Gdrag_r[18]*rdv2-0.7071067811865475*Gdrag_l[18]*rdv2; 
  out[25] += 0.7071067811865475*Gdrag_r[19]*rdv2-0.7071067811865475*Gdrag_l[19]*rdv2; 
  out[26] += 1.224744871391589*Gdrag_r[6]*rdv2+1.224744871391589*Gdrag_l[6]*rdv2; 
  out[27] += 1.224744871391589*Gdrag_r[7]*rdv2+1.224744871391589*Gdrag_l[7]*rdv2; 
  out[28] += 1.224744871391589*Gdrag_r[8]*rdv2+1.224744871391589*Gdrag_l[8]*rdv2; 
  out[29] += 1.224744871391589*Gdrag_r[9]*rdv2+1.224744871391589*Gdrag_l[9]*rdv2; 
  out[30] += 1.224744871391589*Gdrag_r[10]*rdv2+1.224744871391589*Gdrag_l[10]*rdv2; 
  out[31] += 1.224744871391589*Gdrag_r[11]*rdv2+1.224744871391589*Gdrag_l[11]*rdv2; 
  out[32] += 0.7071067811865475*Gdrag_r[20]*rdv2-0.7071067811865475*Gdrag_l[20]*rdv2; 
  out[33] += 0.7071067811865475*Gdrag_r[21]*rdv2-0.7071067811865475*Gdrag_l[21]*rdv2; 
  out[34] += 0.7071067811865475*Gdrag_r[22]*rdv2-0.7071067811865475*Gdrag_l[22]*rdv2; 
  out[35] += 0.7071067811865475*Gdrag_r[23]*rdv2-0.7071067811865475*Gdrag_l[23]*rdv2; 
  out[36] += 0.7071067811865475*Gdrag_r[24]*rdv2-0.7071067811865475*Gdrag_l[24]*rdv2; 
  out[37] += 0.7071067811865475*Gdrag_r[25]*rdv2-0.7071067811865475*Gdrag_l[25]*rdv2; 
  out[38] += 1.224744871391589*Gdrag_r[12]*rdv2+1.224744871391589*Gdrag_l[12]*rdv2; 
  out[39] += 1.224744871391589*Gdrag_r[13]*rdv2+1.224744871391589*Gdrag_l[13]*rdv2; 
  out[40] += 1.224744871391589*Gdrag_r[14]*rdv2+1.224744871391589*Gdrag_l[14]*rdv2; 
  out[41] += 1.224744871391589*Gdrag_r[15]*rdv2+1.224744871391589*Gdrag_l[15]*rdv2; 
  out[42] += 0.7071067811865475*Gdrag_r[26]*rdv2-0.7071067811865475*Gdrag_l[26]*rdv2; 
  out[43] += 1.224744871391589*Gdrag_r[16]*rdv2+1.224744871391589*Gdrag_l[16]*rdv2; 
  out[44] += 1.224744871391589*Gdrag_r[17]*rdv2+1.224744871391589*Gdrag_l[17]*rdv2; 
  out[45] += 1.224744871391589*Gdrag_r[18]*rdv2+1.224744871391589*Gdrag_l[18]*rdv2; 
  out[46] += 1.224744871391589*Gdrag_r[19]*rdv2+1.224744871391589*Gdrag_l[19]*rdv2; 
  out[47] += 0.7071067811865475*Gdrag_r[27]*rdv2-0.7071067811865475*Gdrag_l[27]*rdv2; 
  out[48] += 0.7071067811865475*Gdrag_r[28]*rdv2-0.7071067811865475*Gdrag_l[28]*rdv2; 
  out[49] += 0.7071067811865475*Gdrag_r[29]*rdv2-0.7071067811865475*Gdrag_l[29]*rdv2; 
  out[50] += 0.7071067811865475*Gdrag_r[30]*rdv2-0.7071067811865475*Gdrag_l[30]*rdv2; 
  out[51] += 1.224744871391589*Gdrag_r[20]*rdv2+1.224744871391589*Gdrag_l[20]*rdv2; 
  out[52] += 1.224744871391589*Gdrag_r[21]*rdv2+1.224744871391589*Gdrag_l[21]*rdv2; 
  out[53] += 1.224744871391589*Gdrag_r[22]*rdv2+1.224744871391589*Gdrag_l[22]*rdv2; 
  out[54] += 1.224744871391589*Gdrag_r[23]*rdv2+1.224744871391589*Gdrag_l[23]*rdv2; 
  out[55] += 1.224744871391589*Gdrag_r[24]*rdv2+1.224744871391589*Gdrag_l[24]*rdv2; 
  out[56] += 1.224744871391589*Gdrag_r[25]*rdv2+1.224744871391589*Gdrag_l[25]*rdv2; 
  out[57] += 1.224744871391589*Gdrag_r[26]*rdv2+1.224744871391589*Gdrag_l[26]*rdv2; 
  out[58] += 0.7071067811865475*Gdrag_r[31]*rdv2-0.7071067811865475*Gdrag_l[31]*rdv2; 
  out[59] += 1.224744871391589*Gdrag_r[27]*rdv2+1.224744871391589*Gdrag_l[27]*rdv2; 
  out[60] += 1.224744871391589*Gdrag_r[28]*rdv2+1.224744871391589*Gdrag_l[28]*rdv2; 
  out[61] += 1.224744871391589*Gdrag_r[29]*rdv2+1.224744871391589*Gdrag_l[29]*rdv2; 
  out[62] += 1.224744871391589*Gdrag_r[30]*rdv2+1.224744871391589*Gdrag_l[30]*rdv2; 
  out[63] += 1.224744871391589*Gdrag_r[31]*rdv2+1.224744871391589*Gdrag_l[31]*rdv2; 
} 