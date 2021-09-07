#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_tensor_1x3v_p2_surfvz_quad.h> 
GKYL_CU_DH void vlasov_poisson_extem_boundary_surfvz_1x3v_tensor_p2(const double *w, const double *dxv, const double *fac_phi, const double *vecA, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // fac_phi:     potential (scaled by appropriate factors).
  // vecA:        vector potential (scaled by appropriate factors). Unused in pure Vlasov-Poisson. 
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv12 = 2/dxv[3]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double dv3 = dxv[3], wv3 = w[3]; 
  const double *phi = &fac_phi[0]; 
  const double dx10 = 2/dxv[0]; 
  const double *A0 = &vecA[0]; 
  const double *A1 = &vecA[3]; 
  const double *A2 = &vecA[6]; 
  double alpha[27] = {0.0}; 

  alpha[0] = -3.464101615137754*A2[1]*dx10*wv1; 
  alpha[1] = -7.745966692414834*A2[2]*dx10*wv1; 
  alpha[2] = -1.0*A2[1]*dv1*dx10; 
  alpha[4] = -2.23606797749979*A2[2]*dv1*dx10; 

  double fUpwindQuad[27] = {0.0};
  double fUpwind[27] = {0.0};
  double Ghat[27] = {0.0}; 

  if (edge == -1) { 

  if (0.6363961030678926*alpha[4]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[0] = tensor_1x3v_p2_surfvz_quad_0(1, fSkin); 
  } else { 

    fUpwindQuad[0] = tensor_1x3v_p2_surfvz_quad_0(-1, fEdge); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[2] > 0) { 

    fUpwindQuad[1] = tensor_1x3v_p2_surfvz_quad_1(1, fSkin); 
  } else { 

    fUpwindQuad[1] = tensor_1x3v_p2_surfvz_quad_1(-1, fEdge); 
  } 
  if ((-0.6363961030678926*alpha[4])-0.4743416490252568*alpha[2]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[2] = tensor_1x3v_p2_surfvz_quad_2(1, fSkin); 
  } else { 

    fUpwindQuad[2] = tensor_1x3v_p2_surfvz_quad_2(-1, fEdge); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[1] > 0) { 

    fUpwindQuad[3] = tensor_1x3v_p2_surfvz_quad_3(1, fSkin); 
  } else { 

    fUpwindQuad[3] = tensor_1x3v_p2_surfvz_quad_3(-1, fEdge); 
  } 
  if (0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[4] = tensor_1x3v_p2_surfvz_quad_4(1, fSkin); 
  } else { 

    fUpwindQuad[4] = tensor_1x3v_p2_surfvz_quad_4(-1, fEdge); 
  } 
  if (0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[5] = tensor_1x3v_p2_surfvz_quad_5(1, fSkin); 
  } else { 

    fUpwindQuad[5] = tensor_1x3v_p2_surfvz_quad_5(-1, fEdge); 
  } 
  if ((-0.6363961030678926*alpha[4])+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[6] = tensor_1x3v_p2_surfvz_quad_6(1, fSkin); 
  } else { 

    fUpwindQuad[6] = tensor_1x3v_p2_surfvz_quad_6(-1, fEdge); 
  } 
  if (0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[7] = tensor_1x3v_p2_surfvz_quad_7(1, fSkin); 
  } else { 

    fUpwindQuad[7] = tensor_1x3v_p2_surfvz_quad_7(-1, fEdge); 
  } 
  if (0.6363961030678926*alpha[4]+0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[8] = tensor_1x3v_p2_surfvz_quad_8(1, fSkin); 
  } else { 

    fUpwindQuad[8] = tensor_1x3v_p2_surfvz_quad_8(-1, fEdge); 
  } 
  if (0.6363961030678926*alpha[4]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[9] = tensor_1x3v_p2_surfvz_quad_9(1, fSkin); 
  } else { 

    fUpwindQuad[9] = tensor_1x3v_p2_surfvz_quad_9(-1, fEdge); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[2] > 0) { 

    fUpwindQuad[10] = tensor_1x3v_p2_surfvz_quad_10(1, fSkin); 
  } else { 

    fUpwindQuad[10] = tensor_1x3v_p2_surfvz_quad_10(-1, fEdge); 
  } 
  if ((-0.6363961030678926*alpha[4])-0.4743416490252568*alpha[2]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[11] = tensor_1x3v_p2_surfvz_quad_11(1, fSkin); 
  } else { 

    fUpwindQuad[11] = tensor_1x3v_p2_surfvz_quad_11(-1, fEdge); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[1] > 0) { 

    fUpwindQuad[12] = tensor_1x3v_p2_surfvz_quad_12(1, fSkin); 
  } else { 

    fUpwindQuad[12] = tensor_1x3v_p2_surfvz_quad_12(-1, fEdge); 
  } 
  if (0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[13] = tensor_1x3v_p2_surfvz_quad_13(1, fSkin); 
  } else { 

    fUpwindQuad[13] = tensor_1x3v_p2_surfvz_quad_13(-1, fEdge); 
  } 
  if (0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[14] = tensor_1x3v_p2_surfvz_quad_14(1, fSkin); 
  } else { 

    fUpwindQuad[14] = tensor_1x3v_p2_surfvz_quad_14(-1, fEdge); 
  } 
  if ((-0.6363961030678926*alpha[4])+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[15] = tensor_1x3v_p2_surfvz_quad_15(1, fSkin); 
  } else { 

    fUpwindQuad[15] = tensor_1x3v_p2_surfvz_quad_15(-1, fEdge); 
  } 
  if (0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[16] = tensor_1x3v_p2_surfvz_quad_16(1, fSkin); 
  } else { 

    fUpwindQuad[16] = tensor_1x3v_p2_surfvz_quad_16(-1, fEdge); 
  } 
  if (0.6363961030678926*alpha[4]+0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[17] = tensor_1x3v_p2_surfvz_quad_17(1, fSkin); 
  } else { 

    fUpwindQuad[17] = tensor_1x3v_p2_surfvz_quad_17(-1, fEdge); 
  } 
  if (0.6363961030678926*alpha[4]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[18] = tensor_1x3v_p2_surfvz_quad_18(1, fSkin); 
  } else { 

    fUpwindQuad[18] = tensor_1x3v_p2_surfvz_quad_18(-1, fEdge); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[2] > 0) { 

    fUpwindQuad[19] = tensor_1x3v_p2_surfvz_quad_19(1, fSkin); 
  } else { 

    fUpwindQuad[19] = tensor_1x3v_p2_surfvz_quad_19(-1, fEdge); 
  } 
  if ((-0.6363961030678926*alpha[4])-0.4743416490252568*alpha[2]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[20] = tensor_1x3v_p2_surfvz_quad_20(1, fSkin); 
  } else { 

    fUpwindQuad[20] = tensor_1x3v_p2_surfvz_quad_20(-1, fEdge); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[1] > 0) { 

    fUpwindQuad[21] = tensor_1x3v_p2_surfvz_quad_21(1, fSkin); 
  } else { 

    fUpwindQuad[21] = tensor_1x3v_p2_surfvz_quad_21(-1, fEdge); 
  } 
  if (0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[22] = tensor_1x3v_p2_surfvz_quad_22(1, fSkin); 
  } else { 

    fUpwindQuad[22] = tensor_1x3v_p2_surfvz_quad_22(-1, fEdge); 
  } 
  if (0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[23] = tensor_1x3v_p2_surfvz_quad_23(1, fSkin); 
  } else { 

    fUpwindQuad[23] = tensor_1x3v_p2_surfvz_quad_23(-1, fEdge); 
  } 
  if ((-0.6363961030678926*alpha[4])+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[24] = tensor_1x3v_p2_surfvz_quad_24(1, fSkin); 
  } else { 

    fUpwindQuad[24] = tensor_1x3v_p2_surfvz_quad_24(-1, fEdge); 
  } 
  if (0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[25] = tensor_1x3v_p2_surfvz_quad_25(1, fSkin); 
  } else { 

    fUpwindQuad[25] = tensor_1x3v_p2_surfvz_quad_25(-1, fEdge); 
  } 
  if (0.6363961030678926*alpha[4]+0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[26] = tensor_1x3v_p2_surfvz_quad_26(1, fSkin); 
  } else { 

    fUpwindQuad[26] = tensor_1x3v_p2_surfvz_quad_26(-1, fEdge); 
  } 

  fUpwind[0] = 0.06062300936098657*fUpwindQuad[26]+0.09699681497757856*fUpwindQuad[25]+0.06062300936098657*fUpwindQuad[24]+0.09699681497757856*fUpwindQuad[23]+0.1551949039641257*fUpwindQuad[22]+0.09699681497757856*fUpwindQuad[21]+0.06062300936098657*fUpwindQuad[20]+0.09699681497757856*fUpwindQuad[19]+0.06062300936098657*fUpwindQuad[18]+0.09699681497757856*fUpwindQuad[17]+0.1551949039641257*fUpwindQuad[16]+0.09699681497757856*fUpwindQuad[15]+0.1551949039641257*fUpwindQuad[14]+0.2483118463426013*fUpwindQuad[13]+0.1551949039641257*fUpwindQuad[12]+0.09699681497757856*fUpwindQuad[11]+0.1551949039641257*fUpwindQuad[10]+0.09699681497757856*fUpwindQuad[9]+0.06062300936098657*fUpwindQuad[8]+0.09699681497757856*fUpwindQuad[7]+0.06062300936098657*fUpwindQuad[6]+0.09699681497757856*fUpwindQuad[5]+0.1551949039641257*fUpwindQuad[4]+0.09699681497757856*fUpwindQuad[3]+0.06062300936098657*fUpwindQuad[2]+0.09699681497757856*fUpwindQuad[1]+0.06062300936098657*fUpwindQuad[0]; 
  fUpwind[1] = 0.08133430195906327*fUpwindQuad[26]-0.08133430195906327*fUpwindQuad[24]+0.1301348831345013*fUpwindQuad[23]-0.1301348831345013*fUpwindQuad[21]+0.08133430195906327*fUpwindQuad[20]-0.08133430195906327*fUpwindQuad[18]+0.1301348831345013*fUpwindQuad[17]-0.1301348831345013*fUpwindQuad[15]+0.2082158130152021*fUpwindQuad[14]-0.2082158130152021*fUpwindQuad[12]+0.1301348831345013*fUpwindQuad[11]-0.1301348831345013*fUpwindQuad[9]+0.08133430195906327*fUpwindQuad[8]-0.08133430195906327*fUpwindQuad[6]+0.1301348831345013*fUpwindQuad[5]-0.1301348831345013*fUpwindQuad[3]+0.08133430195906327*fUpwindQuad[2]-0.08133430195906327*fUpwindQuad[0]; 
  fUpwind[2] = 0.08133430195906327*fUpwindQuad[26]+0.1301348831345013*fUpwindQuad[25]+0.08133430195906327*fUpwindQuad[24]-0.08133430195906327*fUpwindQuad[20]-0.1301348831345013*fUpwindQuad[19]-0.08133430195906327*fUpwindQuad[18]+0.1301348831345013*fUpwindQuad[17]+0.2082158130152021*fUpwindQuad[16]+0.1301348831345013*fUpwindQuad[15]-0.1301348831345013*fUpwindQuad[11]-0.2082158130152021*fUpwindQuad[10]-0.1301348831345013*fUpwindQuad[9]+0.08133430195906327*fUpwindQuad[8]+0.1301348831345013*fUpwindQuad[7]+0.08133430195906327*fUpwindQuad[6]-0.08133430195906327*fUpwindQuad[2]-0.1301348831345013*fUpwindQuad[1]-0.08133430195906327*fUpwindQuad[0]; 
  fUpwind[3] = 0.08133430195906327*fUpwindQuad[26]+0.1301348831345013*fUpwindQuad[25]+0.08133430195906327*fUpwindQuad[24]+0.1301348831345013*fUpwindQuad[23]+0.2082158130152021*fUpwindQuad[22]+0.1301348831345013*fUpwindQuad[21]+0.08133430195906327*fUpwindQuad[20]+0.1301348831345013*fUpwindQuad[19]+0.08133430195906327*fUpwindQuad[18]-0.08133430195906327*fUpwindQuad[8]-0.1301348831345013*fUpwindQuad[7]-0.08133430195906327*fUpwindQuad[6]-0.1301348831345013*fUpwindQuad[5]-0.2082158130152021*fUpwindQuad[4]-0.1301348831345013*fUpwindQuad[3]-0.08133430195906327*fUpwindQuad[2]-0.1301348831345013*fUpwindQuad[1]-0.08133430195906327*fUpwindQuad[0]; 
  fUpwind[4] = 0.1091214168497758*fUpwindQuad[26]-0.1091214168497758*fUpwindQuad[24]-0.1091214168497758*fUpwindQuad[20]+0.1091214168497758*fUpwindQuad[18]+0.1745942669596414*fUpwindQuad[17]-0.1745942669596414*fUpwindQuad[15]-0.1745942669596414*fUpwindQuad[11]+0.1745942669596414*fUpwindQuad[9]+0.1091214168497758*fUpwindQuad[8]-0.1091214168497758*fUpwindQuad[6]-0.1091214168497758*fUpwindQuad[2]+0.1091214168497758*fUpwindQuad[0]; 
  fUpwind[5] = 0.1091214168497758*fUpwindQuad[26]-0.1091214168497758*fUpwindQuad[24]+0.1745942669596414*fUpwindQuad[23]-0.1745942669596414*fUpwindQuad[21]+0.1091214168497758*fUpwindQuad[20]-0.1091214168497758*fUpwindQuad[18]-0.1091214168497758*fUpwindQuad[8]+0.1091214168497758*fUpwindQuad[6]-0.1745942669596414*fUpwindQuad[5]+0.1745942669596414*fUpwindQuad[3]-0.1091214168497758*fUpwindQuad[2]+0.1091214168497758*fUpwindQuad[0]; 
  fUpwind[6] = 0.1091214168497758*fUpwindQuad[26]+0.1745942669596414*fUpwindQuad[25]+0.1091214168497758*fUpwindQuad[24]-0.1091214168497758*fUpwindQuad[20]-0.1745942669596414*fUpwindQuad[19]-0.1091214168497758*fUpwindQuad[18]-0.1091214168497758*fUpwindQuad[8]-0.1745942669596414*fUpwindQuad[7]-0.1091214168497758*fUpwindQuad[6]+0.1091214168497758*fUpwindQuad[2]+0.1745942669596414*fUpwindQuad[1]+0.1091214168497758*fUpwindQuad[0]; 
  fUpwind[7] = 0.05422286797270884*fUpwindQuad[26]-0.1084457359454177*fUpwindQuad[25]+0.05422286797270884*fUpwindQuad[24]+0.08675658875633419*fUpwindQuad[23]-0.1735131775126684*fUpwindQuad[22]+0.08675658875633419*fUpwindQuad[21]+0.05422286797270884*fUpwindQuad[20]-0.1084457359454177*fUpwindQuad[19]+0.05422286797270884*fUpwindQuad[18]+0.08675658875633419*fUpwindQuad[17]-0.1735131775126684*fUpwindQuad[16]+0.08675658875633419*fUpwindQuad[15]+0.1388105420101347*fUpwindQuad[14]-0.2776210840202695*fUpwindQuad[13]+0.1388105420101347*fUpwindQuad[12]+0.08675658875633419*fUpwindQuad[11]-0.1735131775126684*fUpwindQuad[10]+0.08675658875633419*fUpwindQuad[9]+0.05422286797270884*fUpwindQuad[8]-0.1084457359454177*fUpwindQuad[7]+0.05422286797270884*fUpwindQuad[6]+0.08675658875633419*fUpwindQuad[5]-0.1735131775126684*fUpwindQuad[4]+0.08675658875633419*fUpwindQuad[3]+0.05422286797270884*fUpwindQuad[2]-0.1084457359454177*fUpwindQuad[1]+0.05422286797270884*fUpwindQuad[0]; 
  fUpwind[8] = 0.05422286797270884*fUpwindQuad[26]+0.08675658875633419*fUpwindQuad[25]+0.05422286797270884*fUpwindQuad[24]-0.1084457359454177*fUpwindQuad[23]-0.1735131775126684*fUpwindQuad[22]-0.1084457359454177*fUpwindQuad[21]+0.05422286797270884*fUpwindQuad[20]+0.08675658875633419*fUpwindQuad[19]+0.05422286797270884*fUpwindQuad[18]+0.08675658875633419*fUpwindQuad[17]+0.1388105420101347*fUpwindQuad[16]+0.08675658875633419*fUpwindQuad[15]-0.1735131775126684*fUpwindQuad[14]-0.2776210840202695*fUpwindQuad[13]-0.1735131775126684*fUpwindQuad[12]+0.08675658875633419*fUpwindQuad[11]+0.1388105420101347*fUpwindQuad[10]+0.08675658875633419*fUpwindQuad[9]+0.05422286797270884*fUpwindQuad[8]+0.08675658875633419*fUpwindQuad[7]+0.05422286797270884*fUpwindQuad[6]-0.1084457359454177*fUpwindQuad[5]-0.1735131775126684*fUpwindQuad[4]-0.1084457359454177*fUpwindQuad[3]+0.05422286797270884*fUpwindQuad[2]+0.08675658875633419*fUpwindQuad[1]+0.05422286797270884*fUpwindQuad[0]; 
  fUpwind[9] = 0.05422286797270884*fUpwindQuad[26]+0.08675658875633419*fUpwindQuad[25]+0.05422286797270884*fUpwindQuad[24]+0.08675658875633419*fUpwindQuad[23]+0.1388105420101347*fUpwindQuad[22]+0.08675658875633419*fUpwindQuad[21]+0.05422286797270884*fUpwindQuad[20]+0.08675658875633419*fUpwindQuad[19]+0.05422286797270884*fUpwindQuad[18]-0.1084457359454177*fUpwindQuad[17]-0.1735131775126684*fUpwindQuad[16]-0.1084457359454177*fUpwindQuad[15]-0.1735131775126684*fUpwindQuad[14]-0.2776210840202695*fUpwindQuad[13]-0.1735131775126684*fUpwindQuad[12]-0.1084457359454177*fUpwindQuad[11]-0.1735131775126684*fUpwindQuad[10]-0.1084457359454177*fUpwindQuad[9]+0.05422286797270884*fUpwindQuad[8]+0.08675658875633419*fUpwindQuad[7]+0.05422286797270884*fUpwindQuad[6]+0.08675658875633419*fUpwindQuad[5]+0.1388105420101347*fUpwindQuad[4]+0.08675658875633419*fUpwindQuad[3]+0.05422286797270884*fUpwindQuad[2]+0.08675658875633419*fUpwindQuad[1]+0.05422286797270884*fUpwindQuad[0]; 
  fUpwind[10] = 0.1464017435263139*fUpwindQuad[26]-0.1464017435263139*fUpwindQuad[24]-0.1464017435263139*fUpwindQuad[20]+0.1464017435263139*fUpwindQuad[18]-0.1464017435263139*fUpwindQuad[8]+0.1464017435263139*fUpwindQuad[6]+0.1464017435263139*fUpwindQuad[2]-0.1464017435263139*fUpwindQuad[0]; 
  fUpwind[11] = 0.07274761123318395*fUpwindQuad[26]-0.1454952224663679*fUpwindQuad[25]+0.07274761123318395*fUpwindQuad[24]-0.07274761123318395*fUpwindQuad[20]+0.1454952224663679*fUpwindQuad[19]-0.07274761123318395*fUpwindQuad[18]+0.1163961779730944*fUpwindQuad[17]-0.2327923559461888*fUpwindQuad[16]+0.1163961779730944*fUpwindQuad[15]-0.1163961779730944*fUpwindQuad[11]+0.2327923559461888*fUpwindQuad[10]-0.1163961779730944*fUpwindQuad[9]+0.07274761123318395*fUpwindQuad[8]-0.1454952224663679*fUpwindQuad[7]+0.07274761123318395*fUpwindQuad[6]-0.07274761123318395*fUpwindQuad[2]+0.1454952224663679*fUpwindQuad[1]-0.07274761123318395*fUpwindQuad[0]; 
  fUpwind[12] = 0.07274761123318395*fUpwindQuad[26]-0.07274761123318395*fUpwindQuad[24]-0.1454952224663679*fUpwindQuad[23]+0.1454952224663679*fUpwindQuad[21]+0.07274761123318395*fUpwindQuad[20]-0.07274761123318395*fUpwindQuad[18]+0.1163961779730944*fUpwindQuad[17]-0.1163961779730944*fUpwindQuad[15]-0.2327923559461888*fUpwindQuad[14]+0.2327923559461888*fUpwindQuad[12]+0.1163961779730944*fUpwindQuad[11]-0.1163961779730944*fUpwindQuad[9]+0.07274761123318395*fUpwindQuad[8]-0.07274761123318395*fUpwindQuad[6]-0.1454952224663679*fUpwindQuad[5]+0.1454952224663679*fUpwindQuad[3]+0.07274761123318395*fUpwindQuad[2]-0.07274761123318395*fUpwindQuad[0]; 
  fUpwind[13] = 0.07274761123318395*fUpwindQuad[26]-0.1454952224663679*fUpwindQuad[25]+0.07274761123318395*fUpwindQuad[24]+0.1163961779730944*fUpwindQuad[23]-0.2327923559461888*fUpwindQuad[22]+0.1163961779730944*fUpwindQuad[21]+0.07274761123318395*fUpwindQuad[20]-0.1454952224663679*fUpwindQuad[19]+0.07274761123318395*fUpwindQuad[18]-0.07274761123318395*fUpwindQuad[8]+0.1454952224663679*fUpwindQuad[7]-0.07274761123318395*fUpwindQuad[6]-0.1163961779730944*fUpwindQuad[5]+0.2327923559461888*fUpwindQuad[4]-0.1163961779730944*fUpwindQuad[3]-0.07274761123318395*fUpwindQuad[2]+0.1454952224663679*fUpwindQuad[1]-0.07274761123318395*fUpwindQuad[0]; 
  fUpwind[14] = 0.07274761123318395*fUpwindQuad[26]+0.1163961779730944*fUpwindQuad[25]+0.07274761123318395*fUpwindQuad[24]-0.1454952224663679*fUpwindQuad[23]-0.2327923559461888*fUpwindQuad[22]-0.1454952224663679*fUpwindQuad[21]+0.07274761123318395*fUpwindQuad[20]+0.1163961779730944*fUpwindQuad[19]+0.07274761123318395*fUpwindQuad[18]-0.07274761123318395*fUpwindQuad[8]-0.1163961779730944*fUpwindQuad[7]-0.07274761123318395*fUpwindQuad[6]+0.1454952224663679*fUpwindQuad[5]+0.2327923559461888*fUpwindQuad[4]+0.1454952224663679*fUpwindQuad[3]-0.07274761123318395*fUpwindQuad[2]-0.1163961779730944*fUpwindQuad[1]-0.07274761123318395*fUpwindQuad[0]; 
  fUpwind[15] = 0.07274761123318395*fUpwindQuad[26]-0.07274761123318395*fUpwindQuad[24]+0.1163961779730944*fUpwindQuad[23]-0.1163961779730944*fUpwindQuad[21]+0.07274761123318395*fUpwindQuad[20]-0.07274761123318395*fUpwindQuad[18]-0.1454952224663679*fUpwindQuad[17]+0.1454952224663679*fUpwindQuad[15]-0.2327923559461888*fUpwindQuad[14]+0.2327923559461888*fUpwindQuad[12]-0.1454952224663679*fUpwindQuad[11]+0.1454952224663679*fUpwindQuad[9]+0.07274761123318395*fUpwindQuad[8]-0.07274761123318395*fUpwindQuad[6]+0.1163961779730944*fUpwindQuad[5]-0.1163961779730944*fUpwindQuad[3]+0.07274761123318395*fUpwindQuad[2]-0.07274761123318395*fUpwindQuad[0]; 
  fUpwind[16] = 0.07274761123318395*fUpwindQuad[26]+0.1163961779730944*fUpwindQuad[25]+0.07274761123318395*fUpwindQuad[24]-0.07274761123318395*fUpwindQuad[20]-0.1163961779730944*fUpwindQuad[19]-0.07274761123318395*fUpwindQuad[18]-0.1454952224663679*fUpwindQuad[17]-0.2327923559461888*fUpwindQuad[16]-0.1454952224663679*fUpwindQuad[15]+0.1454952224663679*fUpwindQuad[11]+0.2327923559461888*fUpwindQuad[10]+0.1454952224663679*fUpwindQuad[9]+0.07274761123318395*fUpwindQuad[8]+0.1163961779730944*fUpwindQuad[7]+0.07274761123318395*fUpwindQuad[6]-0.07274761123318395*fUpwindQuad[2]-0.1163961779730944*fUpwindQuad[1]-0.07274761123318395*fUpwindQuad[0]; 
  fUpwind[17] = 0.09760116235087592*fUpwindQuad[26]-0.1952023247017519*fUpwindQuad[25]+0.09760116235087592*fUpwindQuad[24]-0.09760116235087592*fUpwindQuad[20]+0.1952023247017519*fUpwindQuad[19]-0.09760116235087592*fUpwindQuad[18]-0.09760116235087592*fUpwindQuad[8]+0.1952023247017519*fUpwindQuad[7]-0.09760116235087592*fUpwindQuad[6]+0.09760116235087592*fUpwindQuad[2]-0.1952023247017519*fUpwindQuad[1]+0.09760116235087592*fUpwindQuad[0]; 
  fUpwind[18] = 0.09760116235087592*fUpwindQuad[26]-0.09760116235087592*fUpwindQuad[24]-0.1952023247017519*fUpwindQuad[23]+0.1952023247017519*fUpwindQuad[21]+0.09760116235087592*fUpwindQuad[20]-0.09760116235087592*fUpwindQuad[18]-0.09760116235087592*fUpwindQuad[8]+0.09760116235087592*fUpwindQuad[6]+0.1952023247017519*fUpwindQuad[5]-0.1952023247017519*fUpwindQuad[3]-0.09760116235087592*fUpwindQuad[2]+0.09760116235087592*fUpwindQuad[0]; 
  fUpwind[19] = 0.09760116235087592*fUpwindQuad[26]-0.09760116235087592*fUpwindQuad[24]-0.09760116235087592*fUpwindQuad[20]+0.09760116235087592*fUpwindQuad[18]-0.1952023247017519*fUpwindQuad[17]+0.1952023247017519*fUpwindQuad[15]+0.1952023247017519*fUpwindQuad[11]-0.1952023247017519*fUpwindQuad[9]+0.09760116235087592*fUpwindQuad[8]-0.09760116235087592*fUpwindQuad[6]-0.09760116235087592*fUpwindQuad[2]+0.09760116235087592*fUpwindQuad[0]; 
  fUpwind[20] = 0.04849840748878927*fUpwindQuad[26]-0.09699681497757856*fUpwindQuad[25]+0.04849840748878927*fUpwindQuad[24]-0.09699681497757856*fUpwindQuad[23]+0.1939936299551571*fUpwindQuad[22]-0.09699681497757856*fUpwindQuad[21]+0.04849840748878927*fUpwindQuad[20]-0.09699681497757856*fUpwindQuad[19]+0.04849840748878927*fUpwindQuad[18]+0.07759745198206287*fUpwindQuad[17]-0.1551949039641257*fUpwindQuad[16]+0.07759745198206287*fUpwindQuad[15]-0.1551949039641257*fUpwindQuad[14]+0.3103898079282515*fUpwindQuad[13]-0.1551949039641257*fUpwindQuad[12]+0.07759745198206287*fUpwindQuad[11]-0.1551949039641257*fUpwindQuad[10]+0.07759745198206287*fUpwindQuad[9]+0.04849840748878927*fUpwindQuad[8]-0.09699681497757856*fUpwindQuad[7]+0.04849840748878927*fUpwindQuad[6]-0.09699681497757856*fUpwindQuad[5]+0.1939936299551571*fUpwindQuad[4]-0.09699681497757856*fUpwindQuad[3]+0.04849840748878927*fUpwindQuad[2]-0.09699681497757856*fUpwindQuad[1]+0.04849840748878927*fUpwindQuad[0]; 
  fUpwind[21] = 0.04849840748878927*fUpwindQuad[26]-0.09699681497757856*fUpwindQuad[25]+0.04849840748878927*fUpwindQuad[24]+0.07759745198206287*fUpwindQuad[23]-0.1551949039641257*fUpwindQuad[22]+0.07759745198206287*fUpwindQuad[21]+0.04849840748878927*fUpwindQuad[20]-0.09699681497757856*fUpwindQuad[19]+0.04849840748878927*fUpwindQuad[18]-0.09699681497757856*fUpwindQuad[17]+0.1939936299551571*fUpwindQuad[16]-0.09699681497757856*fUpwindQuad[15]-0.1551949039641257*fUpwindQuad[14]+0.3103898079282515*fUpwindQuad[13]-0.1551949039641257*fUpwindQuad[12]-0.09699681497757856*fUpwindQuad[11]+0.1939936299551571*fUpwindQuad[10]-0.09699681497757856*fUpwindQuad[9]+0.04849840748878927*fUpwindQuad[8]-0.09699681497757856*fUpwindQuad[7]+0.04849840748878927*fUpwindQuad[6]+0.07759745198206287*fUpwindQuad[5]-0.1551949039641257*fUpwindQuad[4]+0.07759745198206287*fUpwindQuad[3]+0.04849840748878927*fUpwindQuad[2]-0.09699681497757856*fUpwindQuad[1]+0.04849840748878927*fUpwindQuad[0]; 
  fUpwind[22] = 0.04849840748878927*fUpwindQuad[26]+0.07759745198206287*fUpwindQuad[25]+0.04849840748878927*fUpwindQuad[24]-0.09699681497757856*fUpwindQuad[23]-0.1551949039641257*fUpwindQuad[22]-0.09699681497757856*fUpwindQuad[21]+0.04849840748878927*fUpwindQuad[20]+0.07759745198206287*fUpwindQuad[19]+0.04849840748878927*fUpwindQuad[18]-0.09699681497757856*fUpwindQuad[17]-0.1551949039641257*fUpwindQuad[16]-0.09699681497757856*fUpwindQuad[15]+0.1939936299551571*fUpwindQuad[14]+0.3103898079282515*fUpwindQuad[13]+0.1939936299551571*fUpwindQuad[12]-0.09699681497757856*fUpwindQuad[11]-0.1551949039641257*fUpwindQuad[10]-0.09699681497757856*fUpwindQuad[9]+0.04849840748878927*fUpwindQuad[8]+0.07759745198206287*fUpwindQuad[7]+0.04849840748878927*fUpwindQuad[6]-0.09699681497757856*fUpwindQuad[5]-0.1551949039641257*fUpwindQuad[4]-0.09699681497757856*fUpwindQuad[3]+0.04849840748878927*fUpwindQuad[2]+0.07759745198206287*fUpwindQuad[1]+0.04849840748878927*fUpwindQuad[0]; 
  fUpwind[23] = 0.06506744156725063*fUpwindQuad[26]-0.1301348831345013*fUpwindQuad[25]+0.06506744156725063*fUpwindQuad[24]-0.1301348831345013*fUpwindQuad[23]+0.2602697662690026*fUpwindQuad[22]-0.1301348831345013*fUpwindQuad[21]+0.06506744156725063*fUpwindQuad[20]-0.1301348831345013*fUpwindQuad[19]+0.06506744156725063*fUpwindQuad[18]-0.06506744156725063*fUpwindQuad[8]+0.1301348831345013*fUpwindQuad[7]-0.06506744156725063*fUpwindQuad[6]+0.1301348831345013*fUpwindQuad[5]-0.2602697662690026*fUpwindQuad[4]+0.1301348831345013*fUpwindQuad[3]-0.06506744156725063*fUpwindQuad[2]+0.1301348831345013*fUpwindQuad[1]-0.06506744156725063*fUpwindQuad[0]; 
  fUpwind[24] = 0.06506744156725063*fUpwindQuad[26]-0.1301348831345013*fUpwindQuad[25]+0.06506744156725063*fUpwindQuad[24]-0.06506744156725063*fUpwindQuad[20]+0.1301348831345013*fUpwindQuad[19]-0.06506744156725063*fUpwindQuad[18]-0.1301348831345013*fUpwindQuad[17]+0.2602697662690026*fUpwindQuad[16]-0.1301348831345013*fUpwindQuad[15]+0.1301348831345013*fUpwindQuad[11]-0.2602697662690026*fUpwindQuad[10]+0.1301348831345013*fUpwindQuad[9]+0.06506744156725063*fUpwindQuad[8]-0.1301348831345013*fUpwindQuad[7]+0.06506744156725063*fUpwindQuad[6]-0.06506744156725063*fUpwindQuad[2]+0.1301348831345013*fUpwindQuad[1]-0.06506744156725063*fUpwindQuad[0]; 
  fUpwind[25] = 0.06506744156725063*fUpwindQuad[26]-0.06506744156725063*fUpwindQuad[24]-0.1301348831345013*fUpwindQuad[23]+0.1301348831345013*fUpwindQuad[21]+0.06506744156725063*fUpwindQuad[20]-0.06506744156725063*fUpwindQuad[18]-0.1301348831345013*fUpwindQuad[17]+0.1301348831345013*fUpwindQuad[15]+0.2602697662690026*fUpwindQuad[14]-0.2602697662690026*fUpwindQuad[12]-0.1301348831345013*fUpwindQuad[11]+0.1301348831345013*fUpwindQuad[9]+0.06506744156725063*fUpwindQuad[8]-0.06506744156725063*fUpwindQuad[6]-0.1301348831345013*fUpwindQuad[5]+0.1301348831345013*fUpwindQuad[3]+0.06506744156725063*fUpwindQuad[2]-0.06506744156725063*fUpwindQuad[0]; 
  fUpwind[26] = 0.04337829437816709*fUpwindQuad[26]-0.08675658875633419*fUpwindQuad[25]+0.04337829437816709*fUpwindQuad[24]-0.08675658875633419*fUpwindQuad[23]+0.1735131775126684*fUpwindQuad[22]-0.08675658875633419*fUpwindQuad[21]+0.04337829437816709*fUpwindQuad[20]-0.08675658875633419*fUpwindQuad[19]+0.04337829437816709*fUpwindQuad[18]-0.08675658875633419*fUpwindQuad[17]+0.1735131775126684*fUpwindQuad[16]-0.08675658875633419*fUpwindQuad[15]+0.1735131775126684*fUpwindQuad[14]-0.3470263550253368*fUpwindQuad[13]+0.1735131775126684*fUpwindQuad[12]-0.08675658875633419*fUpwindQuad[11]+0.1735131775126684*fUpwindQuad[10]-0.08675658875633419*fUpwindQuad[9]+0.04337829437816709*fUpwindQuad[8]-0.08675658875633419*fUpwindQuad[7]+0.04337829437816709*fUpwindQuad[6]-0.08675658875633419*fUpwindQuad[5]+0.1735131775126684*fUpwindQuad[4]-0.08675658875633419*fUpwindQuad[3]+0.04337829437816709*fUpwindQuad[2]-0.08675658875633419*fUpwindQuad[1]+0.04337829437816709*fUpwindQuad[0]; 

  Ghat[0] += 0.3535533905932737*(alpha[4]*fUpwind[4]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] += 0.3162277660168379*alpha[4]*fUpwind[11]+0.3162277660168379*alpha[1]*fUpwind[7]+0.3535533905932737*(alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] += 0.3162277660168379*alpha[4]*fUpwind[12]+0.3162277660168379*alpha[2]*fUpwind[8]+0.3535533905932737*(alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] += 0.3535533905932737*(alpha[4]*fUpwind[10]+alpha[2]*fUpwind[6]+alpha[1]*fUpwind[5]+alpha[0]*fUpwind[3]); 
  Ghat[4] += 0.2828427124746191*alpha[4]*fUpwind[20]+0.3162277660168379*(alpha[2]*fUpwind[12]+alpha[1]*fUpwind[11])+0.3162277660168379*alpha[4]*(fUpwind[8]+fUpwind[7])+0.3535533905932737*(alpha[0]*fUpwind[4]+fUpwind[0]*alpha[4]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[5] += 0.3162277660168379*alpha[4]*fUpwind[17]+0.3162277660168379*alpha[1]*fUpwind[13]+0.3535533905932737*(alpha[2]*fUpwind[10]+alpha[4]*fUpwind[6]+alpha[0]*fUpwind[5]+alpha[1]*fUpwind[3]); 
  Ghat[6] += 0.3162277660168379*alpha[4]*fUpwind[18]+0.3162277660168379*alpha[2]*fUpwind[14]+0.3535533905932737*(alpha[1]*fUpwind[10]+alpha[0]*fUpwind[6]+alpha[4]*fUpwind[5]+alpha[2]*fUpwind[3]); 
  Ghat[7] += 0.3535533905932737*(alpha[2]*fUpwind[11]+alpha[0]*fUpwind[7])+0.3162277660168379*(alpha[4]*fUpwind[4]+alpha[1]*fUpwind[1]); 
  Ghat[8] += 0.3535533905932737*(alpha[1]*fUpwind[12]+alpha[0]*fUpwind[8])+0.3162277660168379*(alpha[4]*fUpwind[4]+alpha[2]*fUpwind[2]); 
  Ghat[9] += 0.3535533905932737*(alpha[4]*fUpwind[19]+alpha[2]*fUpwind[16]+alpha[1]*fUpwind[15]+alpha[0]*fUpwind[9]); 
  Ghat[10] += 0.2828427124746191*alpha[4]*fUpwind[23]+0.3162277660168379*(alpha[2]*fUpwind[18]+alpha[1]*fUpwind[17])+0.3162277660168379*alpha[4]*(fUpwind[14]+fUpwind[13])+0.3535533905932737*(alpha[0]*fUpwind[10]+alpha[1]*fUpwind[6]+alpha[2]*fUpwind[5]+fUpwind[3]*alpha[4]); 
  Ghat[11] += 0.3162277660168379*alpha[2]*fUpwind[20]+0.2828427124746191*alpha[4]*fUpwind[12]+0.3535533905932737*(alpha[0]*fUpwind[11]+alpha[2]*fUpwind[7])+0.3162277660168379*(alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]); 
  Ghat[12] += 0.3162277660168379*alpha[1]*fUpwind[20]+0.3535533905932737*alpha[0]*fUpwind[12]+0.2828427124746191*alpha[4]*fUpwind[11]+0.3535533905932737*alpha[1]*fUpwind[8]+0.3162277660168379*(alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]); 
  Ghat[13] += 0.3535533905932737*(alpha[2]*fUpwind[17]+alpha[0]*fUpwind[13])+0.3162277660168379*(alpha[4]*fUpwind[10]+alpha[1]*fUpwind[5]); 
  Ghat[14] += 0.3535533905932737*(alpha[1]*fUpwind[18]+alpha[0]*fUpwind[14])+0.3162277660168379*(alpha[4]*fUpwind[10]+alpha[2]*fUpwind[6]); 
  Ghat[15] += 0.3162277660168379*(alpha[4]*fUpwind[24]+alpha[1]*fUpwind[21])+0.3535533905932737*(alpha[2]*fUpwind[19]+alpha[4]*fUpwind[16]+alpha[0]*fUpwind[15]+alpha[1]*fUpwind[9]); 
  Ghat[16] += 0.3162277660168379*(alpha[4]*fUpwind[25]+alpha[2]*fUpwind[22])+0.3535533905932737*(alpha[1]*fUpwind[19]+alpha[0]*fUpwind[16]+alpha[4]*fUpwind[15]+alpha[2]*fUpwind[9]); 
  Ghat[17] += 0.3162277660168379*alpha[2]*fUpwind[23]+0.2828427124746191*alpha[4]*fUpwind[18]+0.3535533905932737*(alpha[0]*fUpwind[17]+alpha[2]*fUpwind[13])+0.3162277660168379*(alpha[1]*fUpwind[10]+alpha[4]*fUpwind[5]); 
  Ghat[18] += 0.3162277660168379*alpha[1]*fUpwind[23]+0.3535533905932737*alpha[0]*fUpwind[18]+0.2828427124746191*alpha[4]*fUpwind[17]+0.3535533905932737*alpha[1]*fUpwind[14]+0.3162277660168379*(alpha[2]*fUpwind[10]+alpha[4]*fUpwind[6]); 
  Ghat[19] += 0.2828427124746191*alpha[4]*fUpwind[26]+0.3162277660168379*(alpha[2]*fUpwind[25]+alpha[1]*fUpwind[24]+alpha[4]*(fUpwind[22]+fUpwind[21]))+0.3535533905932737*(alpha[0]*fUpwind[19]+alpha[1]*fUpwind[16]+alpha[2]*fUpwind[15]+alpha[4]*fUpwind[9]); 
  Ghat[20] += 0.3535533905932737*alpha[0]*fUpwind[20]+0.3162277660168379*(alpha[1]*fUpwind[12]+alpha[2]*fUpwind[11])+0.2828427124746191*alpha[4]*fUpwind[4]; 
  Ghat[21] += 0.3535533905932737*(alpha[2]*fUpwind[24]+alpha[0]*fUpwind[21])+0.3162277660168379*alpha[4]*fUpwind[19]+0.3162277660168379*alpha[1]*fUpwind[15]; 
  Ghat[22] += 0.3535533905932737*(alpha[1]*fUpwind[25]+alpha[0]*fUpwind[22])+0.3162277660168379*alpha[4]*fUpwind[19]+0.3162277660168379*alpha[2]*fUpwind[16]; 
  Ghat[23] += 0.3535533905932737*alpha[0]*fUpwind[23]+0.3162277660168379*(alpha[1]*fUpwind[18]+alpha[2]*fUpwind[17])+0.2828427124746191*alpha[4]*fUpwind[10]; 
  Ghat[24] += 0.3162277660168379*alpha[2]*fUpwind[26]+0.2828427124746191*alpha[4]*fUpwind[25]+0.3535533905932737*(alpha[0]*fUpwind[24]+alpha[2]*fUpwind[21])+0.3162277660168379*alpha[1]*fUpwind[19]+0.3162277660168379*alpha[4]*fUpwind[15]; 
  Ghat[25] += 0.3162277660168379*alpha[1]*fUpwind[26]+0.3535533905932737*alpha[0]*fUpwind[25]+0.2828427124746191*alpha[4]*fUpwind[24]+0.3535533905932737*alpha[1]*fUpwind[22]+0.3162277660168379*alpha[2]*fUpwind[19]+0.3162277660168379*alpha[4]*fUpwind[16]; 
  Ghat[26] += 0.3535533905932737*alpha[0]*fUpwind[26]+0.3162277660168379*(alpha[1]*fUpwind[25]+alpha[2]*fUpwind[24])+0.2828427124746191*alpha[4]*fUpwind[19]; 

  out[0] += -0.7071067811865475*Ghat[0]*dv12; 
  out[1] += -0.7071067811865475*Ghat[1]*dv12; 
  out[2] += -0.7071067811865475*Ghat[2]*dv12; 
  out[3] += -0.7071067811865475*Ghat[3]*dv12; 
  out[4] += -1.224744871391589*Ghat[0]*dv12; 
  out[5] += -0.7071067811865475*Ghat[4]*dv12; 
  out[6] += -0.7071067811865475*Ghat[5]*dv12; 
  out[7] += -0.7071067811865475*Ghat[6]*dv12; 
  out[8] += -1.224744871391589*Ghat[1]*dv12; 
  out[9] += -1.224744871391589*Ghat[2]*dv12; 
  out[10] += -1.224744871391589*Ghat[3]*dv12; 
  out[11] += -0.7071067811865475*Ghat[7]*dv12; 
  out[12] += -0.7071067811865475*Ghat[8]*dv12; 
  out[13] += -0.7071067811865475*Ghat[9]*dv12; 
  out[14] += -1.58113883008419*Ghat[0]*dv12; 
  out[15] += -0.7071067811865475*Ghat[10]*dv12; 
  out[16] += -1.224744871391589*Ghat[4]*dv12; 
  out[17] += -1.224744871391589*Ghat[5]*dv12; 
  out[18] += -1.224744871391589*Ghat[6]*dv12; 
  out[19] += -0.7071067811865475*Ghat[11]*dv12; 
  out[20] += -0.7071067811865475*Ghat[12]*dv12; 
  out[21] += -0.7071067811865475*Ghat[13]*dv12; 
  out[22] += -0.7071067811865475*Ghat[14]*dv12; 
  out[23] += -0.7071067811865475*Ghat[15]*dv12; 
  out[24] += -0.7071067811865475*Ghat[16]*dv12; 
  out[25] += -1.224744871391589*Ghat[7]*dv12; 
  out[26] += -1.224744871391589*Ghat[8]*dv12; 
  out[27] += -1.224744871391589*Ghat[9]*dv12; 
  out[28] += -1.58113883008419*Ghat[1]*dv12; 
  out[29] += -1.58113883008419*Ghat[2]*dv12; 
  out[30] += -1.58113883008419*Ghat[3]*dv12; 
  out[31] += -1.224744871391589*Ghat[10]*dv12; 
  out[32] += -0.7071067811865475*Ghat[17]*dv12; 
  out[33] += -0.7071067811865475*Ghat[18]*dv12; 
  out[34] += -0.7071067811865475*Ghat[19]*dv12; 
  out[35] += -1.224744871391589*Ghat[11]*dv12; 
  out[36] += -1.224744871391589*Ghat[12]*dv12; 
  out[37] += -1.224744871391589*Ghat[13]*dv12; 
  out[38] += -1.224744871391589*Ghat[14]*dv12; 
  out[39] += -1.224744871391589*Ghat[15]*dv12; 
  out[40] += -1.224744871391589*Ghat[16]*dv12; 
  out[41] += -1.58113883008419*Ghat[4]*dv12; 
  out[42] += -1.58113883008419*Ghat[5]*dv12; 
  out[43] += -1.58113883008419*Ghat[6]*dv12; 
  out[44] += -0.7071067811865475*Ghat[20]*dv12; 
  out[45] += -0.7071067811865475*Ghat[21]*dv12; 
  out[46] += -0.7071067811865475*Ghat[22]*dv12; 
  out[47] += -1.58113883008419*Ghat[7]*dv12; 
  out[48] += -1.58113883008419*Ghat[8]*dv12; 
  out[49] += -1.58113883008419*Ghat[9]*dv12; 
  out[50] += -1.224744871391589*Ghat[17]*dv12; 
  out[51] += -1.224744871391589*Ghat[18]*dv12; 
  out[52] += -1.224744871391589*Ghat[19]*dv12; 
  out[53] += -1.58113883008419*Ghat[10]*dv12; 
  out[54] += -0.7071067811865475*Ghat[23]*dv12; 
  out[55] += -0.7071067811865475*Ghat[24]*dv12; 
  out[56] += -0.7071067811865475*Ghat[25]*dv12; 
  out[57] += -1.224744871391589*Ghat[20]*dv12; 
  out[58] += -1.224744871391589*Ghat[21]*dv12; 
  out[59] += -1.224744871391589*Ghat[22]*dv12; 
  out[60] += -1.58113883008419*Ghat[11]*dv12; 
  out[61] += -1.58113883008419*Ghat[12]*dv12; 
  out[62] += -1.58113883008419*Ghat[13]*dv12; 
  out[63] += -1.58113883008419*Ghat[14]*dv12; 
  out[64] += -1.58113883008419*Ghat[15]*dv12; 
  out[65] += -1.58113883008419*Ghat[16]*dv12; 
  out[66] += -1.224744871391589*Ghat[23]*dv12; 
  out[67] += -1.224744871391589*Ghat[24]*dv12; 
  out[68] += -1.224744871391589*Ghat[25]*dv12; 
  out[69] += -1.58113883008419*Ghat[17]*dv12; 
  out[70] += -1.58113883008419*Ghat[18]*dv12; 
  out[71] += -1.58113883008419*Ghat[19]*dv12; 
  out[72] += -0.7071067811865475*Ghat[26]*dv12; 
  out[73] += -1.58113883008419*Ghat[20]*dv12; 
  out[74] += -1.58113883008419*Ghat[21]*dv12; 
  out[75] += -1.58113883008419*Ghat[22]*dv12; 
  out[76] += -1.224744871391589*Ghat[26]*dv12; 
  out[77] += -1.58113883008419*Ghat[23]*dv12; 
  out[78] += -1.58113883008419*Ghat[24]*dv12; 
  out[79] += -1.58113883008419*Ghat[25]*dv12; 
  out[80] += -1.58113883008419*Ghat[26]*dv12; 

  } else { 

  if (0.6363961030678926*alpha[4]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[0] = tensor_1x3v_p2_surfvz_quad_0(1, fEdge); 
  } else { 

    fUpwindQuad[0] = tensor_1x3v_p2_surfvz_quad_0(-1, fSkin); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[2] > 0) { 

    fUpwindQuad[1] = tensor_1x3v_p2_surfvz_quad_1(1, fEdge); 
  } else { 

    fUpwindQuad[1] = tensor_1x3v_p2_surfvz_quad_1(-1, fSkin); 
  } 
  if ((-0.6363961030678926*alpha[4])-0.4743416490252568*alpha[2]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[2] = tensor_1x3v_p2_surfvz_quad_2(1, fEdge); 
  } else { 

    fUpwindQuad[2] = tensor_1x3v_p2_surfvz_quad_2(-1, fSkin); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[1] > 0) { 

    fUpwindQuad[3] = tensor_1x3v_p2_surfvz_quad_3(1, fEdge); 
  } else { 

    fUpwindQuad[3] = tensor_1x3v_p2_surfvz_quad_3(-1, fSkin); 
  } 
  if (0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[4] = tensor_1x3v_p2_surfvz_quad_4(1, fEdge); 
  } else { 

    fUpwindQuad[4] = tensor_1x3v_p2_surfvz_quad_4(-1, fSkin); 
  } 
  if (0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[5] = tensor_1x3v_p2_surfvz_quad_5(1, fEdge); 
  } else { 

    fUpwindQuad[5] = tensor_1x3v_p2_surfvz_quad_5(-1, fSkin); 
  } 
  if ((-0.6363961030678926*alpha[4])+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[6] = tensor_1x3v_p2_surfvz_quad_6(1, fEdge); 
  } else { 

    fUpwindQuad[6] = tensor_1x3v_p2_surfvz_quad_6(-1, fSkin); 
  } 
  if (0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[7] = tensor_1x3v_p2_surfvz_quad_7(1, fEdge); 
  } else { 

    fUpwindQuad[7] = tensor_1x3v_p2_surfvz_quad_7(-1, fSkin); 
  } 
  if (0.6363961030678926*alpha[4]+0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[8] = tensor_1x3v_p2_surfvz_quad_8(1, fEdge); 
  } else { 

    fUpwindQuad[8] = tensor_1x3v_p2_surfvz_quad_8(-1, fSkin); 
  } 
  if (0.6363961030678926*alpha[4]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[9] = tensor_1x3v_p2_surfvz_quad_9(1, fEdge); 
  } else { 

    fUpwindQuad[9] = tensor_1x3v_p2_surfvz_quad_9(-1, fSkin); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[2] > 0) { 

    fUpwindQuad[10] = tensor_1x3v_p2_surfvz_quad_10(1, fEdge); 
  } else { 

    fUpwindQuad[10] = tensor_1x3v_p2_surfvz_quad_10(-1, fSkin); 
  } 
  if ((-0.6363961030678926*alpha[4])-0.4743416490252568*alpha[2]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[11] = tensor_1x3v_p2_surfvz_quad_11(1, fEdge); 
  } else { 

    fUpwindQuad[11] = tensor_1x3v_p2_surfvz_quad_11(-1, fSkin); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[1] > 0) { 

    fUpwindQuad[12] = tensor_1x3v_p2_surfvz_quad_12(1, fEdge); 
  } else { 

    fUpwindQuad[12] = tensor_1x3v_p2_surfvz_quad_12(-1, fSkin); 
  } 
  if (0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[13] = tensor_1x3v_p2_surfvz_quad_13(1, fEdge); 
  } else { 

    fUpwindQuad[13] = tensor_1x3v_p2_surfvz_quad_13(-1, fSkin); 
  } 
  if (0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[14] = tensor_1x3v_p2_surfvz_quad_14(1, fEdge); 
  } else { 

    fUpwindQuad[14] = tensor_1x3v_p2_surfvz_quad_14(-1, fSkin); 
  } 
  if ((-0.6363961030678926*alpha[4])+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[15] = tensor_1x3v_p2_surfvz_quad_15(1, fEdge); 
  } else { 

    fUpwindQuad[15] = tensor_1x3v_p2_surfvz_quad_15(-1, fSkin); 
  } 
  if (0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[16] = tensor_1x3v_p2_surfvz_quad_16(1, fEdge); 
  } else { 

    fUpwindQuad[16] = tensor_1x3v_p2_surfvz_quad_16(-1, fSkin); 
  } 
  if (0.6363961030678926*alpha[4]+0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[17] = tensor_1x3v_p2_surfvz_quad_17(1, fEdge); 
  } else { 

    fUpwindQuad[17] = tensor_1x3v_p2_surfvz_quad_17(-1, fSkin); 
  } 
  if (0.6363961030678926*alpha[4]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[18] = tensor_1x3v_p2_surfvz_quad_18(1, fEdge); 
  } else { 

    fUpwindQuad[18] = tensor_1x3v_p2_surfvz_quad_18(-1, fSkin); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[2] > 0) { 

    fUpwindQuad[19] = tensor_1x3v_p2_surfvz_quad_19(1, fEdge); 
  } else { 

    fUpwindQuad[19] = tensor_1x3v_p2_surfvz_quad_19(-1, fSkin); 
  } 
  if ((-0.6363961030678926*alpha[4])-0.4743416490252568*alpha[2]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[20] = tensor_1x3v_p2_surfvz_quad_20(1, fEdge); 
  } else { 

    fUpwindQuad[20] = tensor_1x3v_p2_surfvz_quad_20(-1, fSkin); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[1] > 0) { 

    fUpwindQuad[21] = tensor_1x3v_p2_surfvz_quad_21(1, fEdge); 
  } else { 

    fUpwindQuad[21] = tensor_1x3v_p2_surfvz_quad_21(-1, fSkin); 
  } 
  if (0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[22] = tensor_1x3v_p2_surfvz_quad_22(1, fEdge); 
  } else { 

    fUpwindQuad[22] = tensor_1x3v_p2_surfvz_quad_22(-1, fSkin); 
  } 
  if (0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[23] = tensor_1x3v_p2_surfvz_quad_23(1, fEdge); 
  } else { 

    fUpwindQuad[23] = tensor_1x3v_p2_surfvz_quad_23(-1, fSkin); 
  } 
  if ((-0.6363961030678926*alpha[4])+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[24] = tensor_1x3v_p2_surfvz_quad_24(1, fEdge); 
  } else { 

    fUpwindQuad[24] = tensor_1x3v_p2_surfvz_quad_24(-1, fSkin); 
  } 
  if (0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[25] = tensor_1x3v_p2_surfvz_quad_25(1, fEdge); 
  } else { 

    fUpwindQuad[25] = tensor_1x3v_p2_surfvz_quad_25(-1, fSkin); 
  } 
  if (0.6363961030678926*alpha[4]+0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[26] = tensor_1x3v_p2_surfvz_quad_26(1, fEdge); 
  } else { 

    fUpwindQuad[26] = tensor_1x3v_p2_surfvz_quad_26(-1, fSkin); 
  } 

  fUpwind[0] = 0.06062300936098657*fUpwindQuad[26]+0.09699681497757856*fUpwindQuad[25]+0.06062300936098657*fUpwindQuad[24]+0.09699681497757856*fUpwindQuad[23]+0.1551949039641257*fUpwindQuad[22]+0.09699681497757856*fUpwindQuad[21]+0.06062300936098657*fUpwindQuad[20]+0.09699681497757856*fUpwindQuad[19]+0.06062300936098657*fUpwindQuad[18]+0.09699681497757856*fUpwindQuad[17]+0.1551949039641257*fUpwindQuad[16]+0.09699681497757856*fUpwindQuad[15]+0.1551949039641257*fUpwindQuad[14]+0.2483118463426013*fUpwindQuad[13]+0.1551949039641257*fUpwindQuad[12]+0.09699681497757856*fUpwindQuad[11]+0.1551949039641257*fUpwindQuad[10]+0.09699681497757856*fUpwindQuad[9]+0.06062300936098657*fUpwindQuad[8]+0.09699681497757856*fUpwindQuad[7]+0.06062300936098657*fUpwindQuad[6]+0.09699681497757856*fUpwindQuad[5]+0.1551949039641257*fUpwindQuad[4]+0.09699681497757856*fUpwindQuad[3]+0.06062300936098657*fUpwindQuad[2]+0.09699681497757856*fUpwindQuad[1]+0.06062300936098657*fUpwindQuad[0]; 
  fUpwind[1] = 0.08133430195906327*fUpwindQuad[26]-0.08133430195906327*fUpwindQuad[24]+0.1301348831345013*fUpwindQuad[23]-0.1301348831345013*fUpwindQuad[21]+0.08133430195906327*fUpwindQuad[20]-0.08133430195906327*fUpwindQuad[18]+0.1301348831345013*fUpwindQuad[17]-0.1301348831345013*fUpwindQuad[15]+0.2082158130152021*fUpwindQuad[14]-0.2082158130152021*fUpwindQuad[12]+0.1301348831345013*fUpwindQuad[11]-0.1301348831345013*fUpwindQuad[9]+0.08133430195906327*fUpwindQuad[8]-0.08133430195906327*fUpwindQuad[6]+0.1301348831345013*fUpwindQuad[5]-0.1301348831345013*fUpwindQuad[3]+0.08133430195906327*fUpwindQuad[2]-0.08133430195906327*fUpwindQuad[0]; 
  fUpwind[2] = 0.08133430195906327*fUpwindQuad[26]+0.1301348831345013*fUpwindQuad[25]+0.08133430195906327*fUpwindQuad[24]-0.08133430195906327*fUpwindQuad[20]-0.1301348831345013*fUpwindQuad[19]-0.08133430195906327*fUpwindQuad[18]+0.1301348831345013*fUpwindQuad[17]+0.2082158130152021*fUpwindQuad[16]+0.1301348831345013*fUpwindQuad[15]-0.1301348831345013*fUpwindQuad[11]-0.2082158130152021*fUpwindQuad[10]-0.1301348831345013*fUpwindQuad[9]+0.08133430195906327*fUpwindQuad[8]+0.1301348831345013*fUpwindQuad[7]+0.08133430195906327*fUpwindQuad[6]-0.08133430195906327*fUpwindQuad[2]-0.1301348831345013*fUpwindQuad[1]-0.08133430195906327*fUpwindQuad[0]; 
  fUpwind[3] = 0.08133430195906327*fUpwindQuad[26]+0.1301348831345013*fUpwindQuad[25]+0.08133430195906327*fUpwindQuad[24]+0.1301348831345013*fUpwindQuad[23]+0.2082158130152021*fUpwindQuad[22]+0.1301348831345013*fUpwindQuad[21]+0.08133430195906327*fUpwindQuad[20]+0.1301348831345013*fUpwindQuad[19]+0.08133430195906327*fUpwindQuad[18]-0.08133430195906327*fUpwindQuad[8]-0.1301348831345013*fUpwindQuad[7]-0.08133430195906327*fUpwindQuad[6]-0.1301348831345013*fUpwindQuad[5]-0.2082158130152021*fUpwindQuad[4]-0.1301348831345013*fUpwindQuad[3]-0.08133430195906327*fUpwindQuad[2]-0.1301348831345013*fUpwindQuad[1]-0.08133430195906327*fUpwindQuad[0]; 
  fUpwind[4] = 0.1091214168497758*fUpwindQuad[26]-0.1091214168497758*fUpwindQuad[24]-0.1091214168497758*fUpwindQuad[20]+0.1091214168497758*fUpwindQuad[18]+0.1745942669596414*fUpwindQuad[17]-0.1745942669596414*fUpwindQuad[15]-0.1745942669596414*fUpwindQuad[11]+0.1745942669596414*fUpwindQuad[9]+0.1091214168497758*fUpwindQuad[8]-0.1091214168497758*fUpwindQuad[6]-0.1091214168497758*fUpwindQuad[2]+0.1091214168497758*fUpwindQuad[0]; 
  fUpwind[5] = 0.1091214168497758*fUpwindQuad[26]-0.1091214168497758*fUpwindQuad[24]+0.1745942669596414*fUpwindQuad[23]-0.1745942669596414*fUpwindQuad[21]+0.1091214168497758*fUpwindQuad[20]-0.1091214168497758*fUpwindQuad[18]-0.1091214168497758*fUpwindQuad[8]+0.1091214168497758*fUpwindQuad[6]-0.1745942669596414*fUpwindQuad[5]+0.1745942669596414*fUpwindQuad[3]-0.1091214168497758*fUpwindQuad[2]+0.1091214168497758*fUpwindQuad[0]; 
  fUpwind[6] = 0.1091214168497758*fUpwindQuad[26]+0.1745942669596414*fUpwindQuad[25]+0.1091214168497758*fUpwindQuad[24]-0.1091214168497758*fUpwindQuad[20]-0.1745942669596414*fUpwindQuad[19]-0.1091214168497758*fUpwindQuad[18]-0.1091214168497758*fUpwindQuad[8]-0.1745942669596414*fUpwindQuad[7]-0.1091214168497758*fUpwindQuad[6]+0.1091214168497758*fUpwindQuad[2]+0.1745942669596414*fUpwindQuad[1]+0.1091214168497758*fUpwindQuad[0]; 
  fUpwind[7] = 0.05422286797270884*fUpwindQuad[26]-0.1084457359454177*fUpwindQuad[25]+0.05422286797270884*fUpwindQuad[24]+0.08675658875633419*fUpwindQuad[23]-0.1735131775126684*fUpwindQuad[22]+0.08675658875633419*fUpwindQuad[21]+0.05422286797270884*fUpwindQuad[20]-0.1084457359454177*fUpwindQuad[19]+0.05422286797270884*fUpwindQuad[18]+0.08675658875633419*fUpwindQuad[17]-0.1735131775126684*fUpwindQuad[16]+0.08675658875633419*fUpwindQuad[15]+0.1388105420101347*fUpwindQuad[14]-0.2776210840202695*fUpwindQuad[13]+0.1388105420101347*fUpwindQuad[12]+0.08675658875633419*fUpwindQuad[11]-0.1735131775126684*fUpwindQuad[10]+0.08675658875633419*fUpwindQuad[9]+0.05422286797270884*fUpwindQuad[8]-0.1084457359454177*fUpwindQuad[7]+0.05422286797270884*fUpwindQuad[6]+0.08675658875633419*fUpwindQuad[5]-0.1735131775126684*fUpwindQuad[4]+0.08675658875633419*fUpwindQuad[3]+0.05422286797270884*fUpwindQuad[2]-0.1084457359454177*fUpwindQuad[1]+0.05422286797270884*fUpwindQuad[0]; 
  fUpwind[8] = 0.05422286797270884*fUpwindQuad[26]+0.08675658875633419*fUpwindQuad[25]+0.05422286797270884*fUpwindQuad[24]-0.1084457359454177*fUpwindQuad[23]-0.1735131775126684*fUpwindQuad[22]-0.1084457359454177*fUpwindQuad[21]+0.05422286797270884*fUpwindQuad[20]+0.08675658875633419*fUpwindQuad[19]+0.05422286797270884*fUpwindQuad[18]+0.08675658875633419*fUpwindQuad[17]+0.1388105420101347*fUpwindQuad[16]+0.08675658875633419*fUpwindQuad[15]-0.1735131775126684*fUpwindQuad[14]-0.2776210840202695*fUpwindQuad[13]-0.1735131775126684*fUpwindQuad[12]+0.08675658875633419*fUpwindQuad[11]+0.1388105420101347*fUpwindQuad[10]+0.08675658875633419*fUpwindQuad[9]+0.05422286797270884*fUpwindQuad[8]+0.08675658875633419*fUpwindQuad[7]+0.05422286797270884*fUpwindQuad[6]-0.1084457359454177*fUpwindQuad[5]-0.1735131775126684*fUpwindQuad[4]-0.1084457359454177*fUpwindQuad[3]+0.05422286797270884*fUpwindQuad[2]+0.08675658875633419*fUpwindQuad[1]+0.05422286797270884*fUpwindQuad[0]; 
  fUpwind[9] = 0.05422286797270884*fUpwindQuad[26]+0.08675658875633419*fUpwindQuad[25]+0.05422286797270884*fUpwindQuad[24]+0.08675658875633419*fUpwindQuad[23]+0.1388105420101347*fUpwindQuad[22]+0.08675658875633419*fUpwindQuad[21]+0.05422286797270884*fUpwindQuad[20]+0.08675658875633419*fUpwindQuad[19]+0.05422286797270884*fUpwindQuad[18]-0.1084457359454177*fUpwindQuad[17]-0.1735131775126684*fUpwindQuad[16]-0.1084457359454177*fUpwindQuad[15]-0.1735131775126684*fUpwindQuad[14]-0.2776210840202695*fUpwindQuad[13]-0.1735131775126684*fUpwindQuad[12]-0.1084457359454177*fUpwindQuad[11]-0.1735131775126684*fUpwindQuad[10]-0.1084457359454177*fUpwindQuad[9]+0.05422286797270884*fUpwindQuad[8]+0.08675658875633419*fUpwindQuad[7]+0.05422286797270884*fUpwindQuad[6]+0.08675658875633419*fUpwindQuad[5]+0.1388105420101347*fUpwindQuad[4]+0.08675658875633419*fUpwindQuad[3]+0.05422286797270884*fUpwindQuad[2]+0.08675658875633419*fUpwindQuad[1]+0.05422286797270884*fUpwindQuad[0]; 
  fUpwind[10] = 0.1464017435263139*fUpwindQuad[26]-0.1464017435263139*fUpwindQuad[24]-0.1464017435263139*fUpwindQuad[20]+0.1464017435263139*fUpwindQuad[18]-0.1464017435263139*fUpwindQuad[8]+0.1464017435263139*fUpwindQuad[6]+0.1464017435263139*fUpwindQuad[2]-0.1464017435263139*fUpwindQuad[0]; 
  fUpwind[11] = 0.07274761123318395*fUpwindQuad[26]-0.1454952224663679*fUpwindQuad[25]+0.07274761123318395*fUpwindQuad[24]-0.07274761123318395*fUpwindQuad[20]+0.1454952224663679*fUpwindQuad[19]-0.07274761123318395*fUpwindQuad[18]+0.1163961779730944*fUpwindQuad[17]-0.2327923559461888*fUpwindQuad[16]+0.1163961779730944*fUpwindQuad[15]-0.1163961779730944*fUpwindQuad[11]+0.2327923559461888*fUpwindQuad[10]-0.1163961779730944*fUpwindQuad[9]+0.07274761123318395*fUpwindQuad[8]-0.1454952224663679*fUpwindQuad[7]+0.07274761123318395*fUpwindQuad[6]-0.07274761123318395*fUpwindQuad[2]+0.1454952224663679*fUpwindQuad[1]-0.07274761123318395*fUpwindQuad[0]; 
  fUpwind[12] = 0.07274761123318395*fUpwindQuad[26]-0.07274761123318395*fUpwindQuad[24]-0.1454952224663679*fUpwindQuad[23]+0.1454952224663679*fUpwindQuad[21]+0.07274761123318395*fUpwindQuad[20]-0.07274761123318395*fUpwindQuad[18]+0.1163961779730944*fUpwindQuad[17]-0.1163961779730944*fUpwindQuad[15]-0.2327923559461888*fUpwindQuad[14]+0.2327923559461888*fUpwindQuad[12]+0.1163961779730944*fUpwindQuad[11]-0.1163961779730944*fUpwindQuad[9]+0.07274761123318395*fUpwindQuad[8]-0.07274761123318395*fUpwindQuad[6]-0.1454952224663679*fUpwindQuad[5]+0.1454952224663679*fUpwindQuad[3]+0.07274761123318395*fUpwindQuad[2]-0.07274761123318395*fUpwindQuad[0]; 
  fUpwind[13] = 0.07274761123318395*fUpwindQuad[26]-0.1454952224663679*fUpwindQuad[25]+0.07274761123318395*fUpwindQuad[24]+0.1163961779730944*fUpwindQuad[23]-0.2327923559461888*fUpwindQuad[22]+0.1163961779730944*fUpwindQuad[21]+0.07274761123318395*fUpwindQuad[20]-0.1454952224663679*fUpwindQuad[19]+0.07274761123318395*fUpwindQuad[18]-0.07274761123318395*fUpwindQuad[8]+0.1454952224663679*fUpwindQuad[7]-0.07274761123318395*fUpwindQuad[6]-0.1163961779730944*fUpwindQuad[5]+0.2327923559461888*fUpwindQuad[4]-0.1163961779730944*fUpwindQuad[3]-0.07274761123318395*fUpwindQuad[2]+0.1454952224663679*fUpwindQuad[1]-0.07274761123318395*fUpwindQuad[0]; 
  fUpwind[14] = 0.07274761123318395*fUpwindQuad[26]+0.1163961779730944*fUpwindQuad[25]+0.07274761123318395*fUpwindQuad[24]-0.1454952224663679*fUpwindQuad[23]-0.2327923559461888*fUpwindQuad[22]-0.1454952224663679*fUpwindQuad[21]+0.07274761123318395*fUpwindQuad[20]+0.1163961779730944*fUpwindQuad[19]+0.07274761123318395*fUpwindQuad[18]-0.07274761123318395*fUpwindQuad[8]-0.1163961779730944*fUpwindQuad[7]-0.07274761123318395*fUpwindQuad[6]+0.1454952224663679*fUpwindQuad[5]+0.2327923559461888*fUpwindQuad[4]+0.1454952224663679*fUpwindQuad[3]-0.07274761123318395*fUpwindQuad[2]-0.1163961779730944*fUpwindQuad[1]-0.07274761123318395*fUpwindQuad[0]; 
  fUpwind[15] = 0.07274761123318395*fUpwindQuad[26]-0.07274761123318395*fUpwindQuad[24]+0.1163961779730944*fUpwindQuad[23]-0.1163961779730944*fUpwindQuad[21]+0.07274761123318395*fUpwindQuad[20]-0.07274761123318395*fUpwindQuad[18]-0.1454952224663679*fUpwindQuad[17]+0.1454952224663679*fUpwindQuad[15]-0.2327923559461888*fUpwindQuad[14]+0.2327923559461888*fUpwindQuad[12]-0.1454952224663679*fUpwindQuad[11]+0.1454952224663679*fUpwindQuad[9]+0.07274761123318395*fUpwindQuad[8]-0.07274761123318395*fUpwindQuad[6]+0.1163961779730944*fUpwindQuad[5]-0.1163961779730944*fUpwindQuad[3]+0.07274761123318395*fUpwindQuad[2]-0.07274761123318395*fUpwindQuad[0]; 
  fUpwind[16] = 0.07274761123318395*fUpwindQuad[26]+0.1163961779730944*fUpwindQuad[25]+0.07274761123318395*fUpwindQuad[24]-0.07274761123318395*fUpwindQuad[20]-0.1163961779730944*fUpwindQuad[19]-0.07274761123318395*fUpwindQuad[18]-0.1454952224663679*fUpwindQuad[17]-0.2327923559461888*fUpwindQuad[16]-0.1454952224663679*fUpwindQuad[15]+0.1454952224663679*fUpwindQuad[11]+0.2327923559461888*fUpwindQuad[10]+0.1454952224663679*fUpwindQuad[9]+0.07274761123318395*fUpwindQuad[8]+0.1163961779730944*fUpwindQuad[7]+0.07274761123318395*fUpwindQuad[6]-0.07274761123318395*fUpwindQuad[2]-0.1163961779730944*fUpwindQuad[1]-0.07274761123318395*fUpwindQuad[0]; 
  fUpwind[17] = 0.09760116235087592*fUpwindQuad[26]-0.1952023247017519*fUpwindQuad[25]+0.09760116235087592*fUpwindQuad[24]-0.09760116235087592*fUpwindQuad[20]+0.1952023247017519*fUpwindQuad[19]-0.09760116235087592*fUpwindQuad[18]-0.09760116235087592*fUpwindQuad[8]+0.1952023247017519*fUpwindQuad[7]-0.09760116235087592*fUpwindQuad[6]+0.09760116235087592*fUpwindQuad[2]-0.1952023247017519*fUpwindQuad[1]+0.09760116235087592*fUpwindQuad[0]; 
  fUpwind[18] = 0.09760116235087592*fUpwindQuad[26]-0.09760116235087592*fUpwindQuad[24]-0.1952023247017519*fUpwindQuad[23]+0.1952023247017519*fUpwindQuad[21]+0.09760116235087592*fUpwindQuad[20]-0.09760116235087592*fUpwindQuad[18]-0.09760116235087592*fUpwindQuad[8]+0.09760116235087592*fUpwindQuad[6]+0.1952023247017519*fUpwindQuad[5]-0.1952023247017519*fUpwindQuad[3]-0.09760116235087592*fUpwindQuad[2]+0.09760116235087592*fUpwindQuad[0]; 
  fUpwind[19] = 0.09760116235087592*fUpwindQuad[26]-0.09760116235087592*fUpwindQuad[24]-0.09760116235087592*fUpwindQuad[20]+0.09760116235087592*fUpwindQuad[18]-0.1952023247017519*fUpwindQuad[17]+0.1952023247017519*fUpwindQuad[15]+0.1952023247017519*fUpwindQuad[11]-0.1952023247017519*fUpwindQuad[9]+0.09760116235087592*fUpwindQuad[8]-0.09760116235087592*fUpwindQuad[6]-0.09760116235087592*fUpwindQuad[2]+0.09760116235087592*fUpwindQuad[0]; 
  fUpwind[20] = 0.04849840748878927*fUpwindQuad[26]-0.09699681497757856*fUpwindQuad[25]+0.04849840748878927*fUpwindQuad[24]-0.09699681497757856*fUpwindQuad[23]+0.1939936299551571*fUpwindQuad[22]-0.09699681497757856*fUpwindQuad[21]+0.04849840748878927*fUpwindQuad[20]-0.09699681497757856*fUpwindQuad[19]+0.04849840748878927*fUpwindQuad[18]+0.07759745198206287*fUpwindQuad[17]-0.1551949039641257*fUpwindQuad[16]+0.07759745198206287*fUpwindQuad[15]-0.1551949039641257*fUpwindQuad[14]+0.3103898079282515*fUpwindQuad[13]-0.1551949039641257*fUpwindQuad[12]+0.07759745198206287*fUpwindQuad[11]-0.1551949039641257*fUpwindQuad[10]+0.07759745198206287*fUpwindQuad[9]+0.04849840748878927*fUpwindQuad[8]-0.09699681497757856*fUpwindQuad[7]+0.04849840748878927*fUpwindQuad[6]-0.09699681497757856*fUpwindQuad[5]+0.1939936299551571*fUpwindQuad[4]-0.09699681497757856*fUpwindQuad[3]+0.04849840748878927*fUpwindQuad[2]-0.09699681497757856*fUpwindQuad[1]+0.04849840748878927*fUpwindQuad[0]; 
  fUpwind[21] = 0.04849840748878927*fUpwindQuad[26]-0.09699681497757856*fUpwindQuad[25]+0.04849840748878927*fUpwindQuad[24]+0.07759745198206287*fUpwindQuad[23]-0.1551949039641257*fUpwindQuad[22]+0.07759745198206287*fUpwindQuad[21]+0.04849840748878927*fUpwindQuad[20]-0.09699681497757856*fUpwindQuad[19]+0.04849840748878927*fUpwindQuad[18]-0.09699681497757856*fUpwindQuad[17]+0.1939936299551571*fUpwindQuad[16]-0.09699681497757856*fUpwindQuad[15]-0.1551949039641257*fUpwindQuad[14]+0.3103898079282515*fUpwindQuad[13]-0.1551949039641257*fUpwindQuad[12]-0.09699681497757856*fUpwindQuad[11]+0.1939936299551571*fUpwindQuad[10]-0.09699681497757856*fUpwindQuad[9]+0.04849840748878927*fUpwindQuad[8]-0.09699681497757856*fUpwindQuad[7]+0.04849840748878927*fUpwindQuad[6]+0.07759745198206287*fUpwindQuad[5]-0.1551949039641257*fUpwindQuad[4]+0.07759745198206287*fUpwindQuad[3]+0.04849840748878927*fUpwindQuad[2]-0.09699681497757856*fUpwindQuad[1]+0.04849840748878927*fUpwindQuad[0]; 
  fUpwind[22] = 0.04849840748878927*fUpwindQuad[26]+0.07759745198206287*fUpwindQuad[25]+0.04849840748878927*fUpwindQuad[24]-0.09699681497757856*fUpwindQuad[23]-0.1551949039641257*fUpwindQuad[22]-0.09699681497757856*fUpwindQuad[21]+0.04849840748878927*fUpwindQuad[20]+0.07759745198206287*fUpwindQuad[19]+0.04849840748878927*fUpwindQuad[18]-0.09699681497757856*fUpwindQuad[17]-0.1551949039641257*fUpwindQuad[16]-0.09699681497757856*fUpwindQuad[15]+0.1939936299551571*fUpwindQuad[14]+0.3103898079282515*fUpwindQuad[13]+0.1939936299551571*fUpwindQuad[12]-0.09699681497757856*fUpwindQuad[11]-0.1551949039641257*fUpwindQuad[10]-0.09699681497757856*fUpwindQuad[9]+0.04849840748878927*fUpwindQuad[8]+0.07759745198206287*fUpwindQuad[7]+0.04849840748878927*fUpwindQuad[6]-0.09699681497757856*fUpwindQuad[5]-0.1551949039641257*fUpwindQuad[4]-0.09699681497757856*fUpwindQuad[3]+0.04849840748878927*fUpwindQuad[2]+0.07759745198206287*fUpwindQuad[1]+0.04849840748878927*fUpwindQuad[0]; 
  fUpwind[23] = 0.06506744156725063*fUpwindQuad[26]-0.1301348831345013*fUpwindQuad[25]+0.06506744156725063*fUpwindQuad[24]-0.1301348831345013*fUpwindQuad[23]+0.2602697662690026*fUpwindQuad[22]-0.1301348831345013*fUpwindQuad[21]+0.06506744156725063*fUpwindQuad[20]-0.1301348831345013*fUpwindQuad[19]+0.06506744156725063*fUpwindQuad[18]-0.06506744156725063*fUpwindQuad[8]+0.1301348831345013*fUpwindQuad[7]-0.06506744156725063*fUpwindQuad[6]+0.1301348831345013*fUpwindQuad[5]-0.2602697662690026*fUpwindQuad[4]+0.1301348831345013*fUpwindQuad[3]-0.06506744156725063*fUpwindQuad[2]+0.1301348831345013*fUpwindQuad[1]-0.06506744156725063*fUpwindQuad[0]; 
  fUpwind[24] = 0.06506744156725063*fUpwindQuad[26]-0.1301348831345013*fUpwindQuad[25]+0.06506744156725063*fUpwindQuad[24]-0.06506744156725063*fUpwindQuad[20]+0.1301348831345013*fUpwindQuad[19]-0.06506744156725063*fUpwindQuad[18]-0.1301348831345013*fUpwindQuad[17]+0.2602697662690026*fUpwindQuad[16]-0.1301348831345013*fUpwindQuad[15]+0.1301348831345013*fUpwindQuad[11]-0.2602697662690026*fUpwindQuad[10]+0.1301348831345013*fUpwindQuad[9]+0.06506744156725063*fUpwindQuad[8]-0.1301348831345013*fUpwindQuad[7]+0.06506744156725063*fUpwindQuad[6]-0.06506744156725063*fUpwindQuad[2]+0.1301348831345013*fUpwindQuad[1]-0.06506744156725063*fUpwindQuad[0]; 
  fUpwind[25] = 0.06506744156725063*fUpwindQuad[26]-0.06506744156725063*fUpwindQuad[24]-0.1301348831345013*fUpwindQuad[23]+0.1301348831345013*fUpwindQuad[21]+0.06506744156725063*fUpwindQuad[20]-0.06506744156725063*fUpwindQuad[18]-0.1301348831345013*fUpwindQuad[17]+0.1301348831345013*fUpwindQuad[15]+0.2602697662690026*fUpwindQuad[14]-0.2602697662690026*fUpwindQuad[12]-0.1301348831345013*fUpwindQuad[11]+0.1301348831345013*fUpwindQuad[9]+0.06506744156725063*fUpwindQuad[8]-0.06506744156725063*fUpwindQuad[6]-0.1301348831345013*fUpwindQuad[5]+0.1301348831345013*fUpwindQuad[3]+0.06506744156725063*fUpwindQuad[2]-0.06506744156725063*fUpwindQuad[0]; 
  fUpwind[26] = 0.04337829437816709*fUpwindQuad[26]-0.08675658875633419*fUpwindQuad[25]+0.04337829437816709*fUpwindQuad[24]-0.08675658875633419*fUpwindQuad[23]+0.1735131775126684*fUpwindQuad[22]-0.08675658875633419*fUpwindQuad[21]+0.04337829437816709*fUpwindQuad[20]-0.08675658875633419*fUpwindQuad[19]+0.04337829437816709*fUpwindQuad[18]-0.08675658875633419*fUpwindQuad[17]+0.1735131775126684*fUpwindQuad[16]-0.08675658875633419*fUpwindQuad[15]+0.1735131775126684*fUpwindQuad[14]-0.3470263550253368*fUpwindQuad[13]+0.1735131775126684*fUpwindQuad[12]-0.08675658875633419*fUpwindQuad[11]+0.1735131775126684*fUpwindQuad[10]-0.08675658875633419*fUpwindQuad[9]+0.04337829437816709*fUpwindQuad[8]-0.08675658875633419*fUpwindQuad[7]+0.04337829437816709*fUpwindQuad[6]-0.08675658875633419*fUpwindQuad[5]+0.1735131775126684*fUpwindQuad[4]-0.08675658875633419*fUpwindQuad[3]+0.04337829437816709*fUpwindQuad[2]-0.08675658875633419*fUpwindQuad[1]+0.04337829437816709*fUpwindQuad[0]; 

  Ghat[0] += 0.3535533905932737*(alpha[4]*fUpwind[4]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] += 0.3162277660168379*alpha[4]*fUpwind[11]+0.3162277660168379*alpha[1]*fUpwind[7]+0.3535533905932737*(alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] += 0.3162277660168379*alpha[4]*fUpwind[12]+0.3162277660168379*alpha[2]*fUpwind[8]+0.3535533905932737*(alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] += 0.3535533905932737*(alpha[4]*fUpwind[10]+alpha[2]*fUpwind[6]+alpha[1]*fUpwind[5]+alpha[0]*fUpwind[3]); 
  Ghat[4] += 0.2828427124746191*alpha[4]*fUpwind[20]+0.3162277660168379*(alpha[2]*fUpwind[12]+alpha[1]*fUpwind[11])+0.3162277660168379*alpha[4]*(fUpwind[8]+fUpwind[7])+0.3535533905932737*(alpha[0]*fUpwind[4]+fUpwind[0]*alpha[4]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[5] += 0.3162277660168379*alpha[4]*fUpwind[17]+0.3162277660168379*alpha[1]*fUpwind[13]+0.3535533905932737*(alpha[2]*fUpwind[10]+alpha[4]*fUpwind[6]+alpha[0]*fUpwind[5]+alpha[1]*fUpwind[3]); 
  Ghat[6] += 0.3162277660168379*alpha[4]*fUpwind[18]+0.3162277660168379*alpha[2]*fUpwind[14]+0.3535533905932737*(alpha[1]*fUpwind[10]+alpha[0]*fUpwind[6]+alpha[4]*fUpwind[5]+alpha[2]*fUpwind[3]); 
  Ghat[7] += 0.3535533905932737*(alpha[2]*fUpwind[11]+alpha[0]*fUpwind[7])+0.3162277660168379*(alpha[4]*fUpwind[4]+alpha[1]*fUpwind[1]); 
  Ghat[8] += 0.3535533905932737*(alpha[1]*fUpwind[12]+alpha[0]*fUpwind[8])+0.3162277660168379*(alpha[4]*fUpwind[4]+alpha[2]*fUpwind[2]); 
  Ghat[9] += 0.3535533905932737*(alpha[4]*fUpwind[19]+alpha[2]*fUpwind[16]+alpha[1]*fUpwind[15]+alpha[0]*fUpwind[9]); 
  Ghat[10] += 0.2828427124746191*alpha[4]*fUpwind[23]+0.3162277660168379*(alpha[2]*fUpwind[18]+alpha[1]*fUpwind[17])+0.3162277660168379*alpha[4]*(fUpwind[14]+fUpwind[13])+0.3535533905932737*(alpha[0]*fUpwind[10]+alpha[1]*fUpwind[6]+alpha[2]*fUpwind[5]+fUpwind[3]*alpha[4]); 
  Ghat[11] += 0.3162277660168379*alpha[2]*fUpwind[20]+0.2828427124746191*alpha[4]*fUpwind[12]+0.3535533905932737*(alpha[0]*fUpwind[11]+alpha[2]*fUpwind[7])+0.3162277660168379*(alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]); 
  Ghat[12] += 0.3162277660168379*alpha[1]*fUpwind[20]+0.3535533905932737*alpha[0]*fUpwind[12]+0.2828427124746191*alpha[4]*fUpwind[11]+0.3535533905932737*alpha[1]*fUpwind[8]+0.3162277660168379*(alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]); 
  Ghat[13] += 0.3535533905932737*(alpha[2]*fUpwind[17]+alpha[0]*fUpwind[13])+0.3162277660168379*(alpha[4]*fUpwind[10]+alpha[1]*fUpwind[5]); 
  Ghat[14] += 0.3535533905932737*(alpha[1]*fUpwind[18]+alpha[0]*fUpwind[14])+0.3162277660168379*(alpha[4]*fUpwind[10]+alpha[2]*fUpwind[6]); 
  Ghat[15] += 0.3162277660168379*(alpha[4]*fUpwind[24]+alpha[1]*fUpwind[21])+0.3535533905932737*(alpha[2]*fUpwind[19]+alpha[4]*fUpwind[16]+alpha[0]*fUpwind[15]+alpha[1]*fUpwind[9]); 
  Ghat[16] += 0.3162277660168379*(alpha[4]*fUpwind[25]+alpha[2]*fUpwind[22])+0.3535533905932737*(alpha[1]*fUpwind[19]+alpha[0]*fUpwind[16]+alpha[4]*fUpwind[15]+alpha[2]*fUpwind[9]); 
  Ghat[17] += 0.3162277660168379*alpha[2]*fUpwind[23]+0.2828427124746191*alpha[4]*fUpwind[18]+0.3535533905932737*(alpha[0]*fUpwind[17]+alpha[2]*fUpwind[13])+0.3162277660168379*(alpha[1]*fUpwind[10]+alpha[4]*fUpwind[5]); 
  Ghat[18] += 0.3162277660168379*alpha[1]*fUpwind[23]+0.3535533905932737*alpha[0]*fUpwind[18]+0.2828427124746191*alpha[4]*fUpwind[17]+0.3535533905932737*alpha[1]*fUpwind[14]+0.3162277660168379*(alpha[2]*fUpwind[10]+alpha[4]*fUpwind[6]); 
  Ghat[19] += 0.2828427124746191*alpha[4]*fUpwind[26]+0.3162277660168379*(alpha[2]*fUpwind[25]+alpha[1]*fUpwind[24]+alpha[4]*(fUpwind[22]+fUpwind[21]))+0.3535533905932737*(alpha[0]*fUpwind[19]+alpha[1]*fUpwind[16]+alpha[2]*fUpwind[15]+alpha[4]*fUpwind[9]); 
  Ghat[20] += 0.3535533905932737*alpha[0]*fUpwind[20]+0.3162277660168379*(alpha[1]*fUpwind[12]+alpha[2]*fUpwind[11])+0.2828427124746191*alpha[4]*fUpwind[4]; 
  Ghat[21] += 0.3535533905932737*(alpha[2]*fUpwind[24]+alpha[0]*fUpwind[21])+0.3162277660168379*alpha[4]*fUpwind[19]+0.3162277660168379*alpha[1]*fUpwind[15]; 
  Ghat[22] += 0.3535533905932737*(alpha[1]*fUpwind[25]+alpha[0]*fUpwind[22])+0.3162277660168379*alpha[4]*fUpwind[19]+0.3162277660168379*alpha[2]*fUpwind[16]; 
  Ghat[23] += 0.3535533905932737*alpha[0]*fUpwind[23]+0.3162277660168379*(alpha[1]*fUpwind[18]+alpha[2]*fUpwind[17])+0.2828427124746191*alpha[4]*fUpwind[10]; 
  Ghat[24] += 0.3162277660168379*alpha[2]*fUpwind[26]+0.2828427124746191*alpha[4]*fUpwind[25]+0.3535533905932737*(alpha[0]*fUpwind[24]+alpha[2]*fUpwind[21])+0.3162277660168379*alpha[1]*fUpwind[19]+0.3162277660168379*alpha[4]*fUpwind[15]; 
  Ghat[25] += 0.3162277660168379*alpha[1]*fUpwind[26]+0.3535533905932737*alpha[0]*fUpwind[25]+0.2828427124746191*alpha[4]*fUpwind[24]+0.3535533905932737*alpha[1]*fUpwind[22]+0.3162277660168379*alpha[2]*fUpwind[19]+0.3162277660168379*alpha[4]*fUpwind[16]; 
  Ghat[26] += 0.3535533905932737*alpha[0]*fUpwind[26]+0.3162277660168379*(alpha[1]*fUpwind[25]+alpha[2]*fUpwind[24])+0.2828427124746191*alpha[4]*fUpwind[19]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv12; 
  out[1] += 0.7071067811865475*Ghat[1]*dv12; 
  out[2] += 0.7071067811865475*Ghat[2]*dv12; 
  out[3] += 0.7071067811865475*Ghat[3]*dv12; 
  out[4] += -1.224744871391589*Ghat[0]*dv12; 
  out[5] += 0.7071067811865475*Ghat[4]*dv12; 
  out[6] += 0.7071067811865475*Ghat[5]*dv12; 
  out[7] += 0.7071067811865475*Ghat[6]*dv12; 
  out[8] += -1.224744871391589*Ghat[1]*dv12; 
  out[9] += -1.224744871391589*Ghat[2]*dv12; 
  out[10] += -1.224744871391589*Ghat[3]*dv12; 
  out[11] += 0.7071067811865475*Ghat[7]*dv12; 
  out[12] += 0.7071067811865475*Ghat[8]*dv12; 
  out[13] += 0.7071067811865475*Ghat[9]*dv12; 
  out[14] += 1.58113883008419*Ghat[0]*dv12; 
  out[15] += 0.7071067811865475*Ghat[10]*dv12; 
  out[16] += -1.224744871391589*Ghat[4]*dv12; 
  out[17] += -1.224744871391589*Ghat[5]*dv12; 
  out[18] += -1.224744871391589*Ghat[6]*dv12; 
  out[19] += 0.7071067811865475*Ghat[11]*dv12; 
  out[20] += 0.7071067811865475*Ghat[12]*dv12; 
  out[21] += 0.7071067811865475*Ghat[13]*dv12; 
  out[22] += 0.7071067811865475*Ghat[14]*dv12; 
  out[23] += 0.7071067811865475*Ghat[15]*dv12; 
  out[24] += 0.7071067811865475*Ghat[16]*dv12; 
  out[25] += -1.224744871391589*Ghat[7]*dv12; 
  out[26] += -1.224744871391589*Ghat[8]*dv12; 
  out[27] += -1.224744871391589*Ghat[9]*dv12; 
  out[28] += 1.58113883008419*Ghat[1]*dv12; 
  out[29] += 1.58113883008419*Ghat[2]*dv12; 
  out[30] += 1.58113883008419*Ghat[3]*dv12; 
  out[31] += -1.224744871391589*Ghat[10]*dv12; 
  out[32] += 0.7071067811865475*Ghat[17]*dv12; 
  out[33] += 0.7071067811865475*Ghat[18]*dv12; 
  out[34] += 0.7071067811865475*Ghat[19]*dv12; 
  out[35] += -1.224744871391589*Ghat[11]*dv12; 
  out[36] += -1.224744871391589*Ghat[12]*dv12; 
  out[37] += -1.224744871391589*Ghat[13]*dv12; 
  out[38] += -1.224744871391589*Ghat[14]*dv12; 
  out[39] += -1.224744871391589*Ghat[15]*dv12; 
  out[40] += -1.224744871391589*Ghat[16]*dv12; 
  out[41] += 1.58113883008419*Ghat[4]*dv12; 
  out[42] += 1.58113883008419*Ghat[5]*dv12; 
  out[43] += 1.58113883008419*Ghat[6]*dv12; 
  out[44] += 0.7071067811865475*Ghat[20]*dv12; 
  out[45] += 0.7071067811865475*Ghat[21]*dv12; 
  out[46] += 0.7071067811865475*Ghat[22]*dv12; 
  out[47] += 1.58113883008419*Ghat[7]*dv12; 
  out[48] += 1.58113883008419*Ghat[8]*dv12; 
  out[49] += 1.58113883008419*Ghat[9]*dv12; 
  out[50] += -1.224744871391589*Ghat[17]*dv12; 
  out[51] += -1.224744871391589*Ghat[18]*dv12; 
  out[52] += -1.224744871391589*Ghat[19]*dv12; 
  out[53] += 1.58113883008419*Ghat[10]*dv12; 
  out[54] += 0.7071067811865475*Ghat[23]*dv12; 
  out[55] += 0.7071067811865475*Ghat[24]*dv12; 
  out[56] += 0.7071067811865475*Ghat[25]*dv12; 
  out[57] += -1.224744871391589*Ghat[20]*dv12; 
  out[58] += -1.224744871391589*Ghat[21]*dv12; 
  out[59] += -1.224744871391589*Ghat[22]*dv12; 
  out[60] += 1.58113883008419*Ghat[11]*dv12; 
  out[61] += 1.58113883008419*Ghat[12]*dv12; 
  out[62] += 1.58113883008419*Ghat[13]*dv12; 
  out[63] += 1.58113883008419*Ghat[14]*dv12; 
  out[64] += 1.58113883008419*Ghat[15]*dv12; 
  out[65] += 1.58113883008419*Ghat[16]*dv12; 
  out[66] += -1.224744871391589*Ghat[23]*dv12; 
  out[67] += -1.224744871391589*Ghat[24]*dv12; 
  out[68] += -1.224744871391589*Ghat[25]*dv12; 
  out[69] += 1.58113883008419*Ghat[17]*dv12; 
  out[70] += 1.58113883008419*Ghat[18]*dv12; 
  out[71] += 1.58113883008419*Ghat[19]*dv12; 
  out[72] += 0.7071067811865475*Ghat[26]*dv12; 
  out[73] += 1.58113883008419*Ghat[20]*dv12; 
  out[74] += 1.58113883008419*Ghat[21]*dv12; 
  out[75] += 1.58113883008419*Ghat[22]*dv12; 
  out[76] += -1.224744871391589*Ghat[26]*dv12; 
  out[77] += 1.58113883008419*Ghat[23]*dv12; 
  out[78] += 1.58113883008419*Ghat[24]*dv12; 
  out[79] += 1.58113883008419*Ghat[25]*dv12; 
  out[80] += 1.58113883008419*Ghat[26]*dv12; 

  } 
} 