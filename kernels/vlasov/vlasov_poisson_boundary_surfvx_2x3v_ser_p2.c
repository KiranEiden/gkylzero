#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_5x_p2_surfx3_eval_quad.h> 
#include <gkyl_basis_ser_5x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_poisson_boundary_surfvx_2x3v_ser_p2(const double *w, const double *dxv, const double *fac_phi, const double *vecA, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // fac_phi:     potential (scaled by appropriate factors).
  // vecA:        vector potential (scaled by appropriate factors). Unused in pure Vlasov-Poisson. 
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv10 = 2/dxv[2]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double dv3 = dxv[4], wv3 = w[4]; 
  const double *phi = &fac_phi[0]; 
  const double dx10 = 2/dxv[0]; 
  const double dx11 = 2/dxv[1]; 
  double alpha[48] = {0.0}; 

  alpha[0] = -3.464101615137754*phi[1]*dx10; 
  alpha[1] = -7.745966692414834*phi[4]*dx10; 
  alpha[2] = -3.464101615137754*phi[3]*dx10; 
  alpha[5] = -7.745966692414834*phi[6]*dx10; 
  alpha[12] = -3.464101615137755*phi[7]*dx10; 

  double fUpwindQuad[81] = {0.0};
  double fUpwind[48] = {0.0};
  double Ghat[48] = {0.0}; 

  if (edge == -1) { 

  if (0.223606797749979*alpha[12]+0.45*alpha[5]-0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[0] = ser_5x_p2_surfx3_eval_quad_node_0_r(fSkin); 
    fUpwindQuad[1] = ser_5x_p2_surfx3_eval_quad_node_1_r(fSkin); 
    fUpwindQuad[2] = ser_5x_p2_surfx3_eval_quad_node_2_r(fSkin); 
    fUpwindQuad[3] = ser_5x_p2_surfx3_eval_quad_node_3_r(fSkin); 
    fUpwindQuad[4] = ser_5x_p2_surfx3_eval_quad_node_4_r(fSkin); 
    fUpwindQuad[5] = ser_5x_p2_surfx3_eval_quad_node_5_r(fSkin); 
    fUpwindQuad[6] = ser_5x_p2_surfx3_eval_quad_node_6_r(fSkin); 
    fUpwindQuad[7] = ser_5x_p2_surfx3_eval_quad_node_7_r(fSkin); 
    fUpwindQuad[8] = ser_5x_p2_surfx3_eval_quad_node_8_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_5x_p2_surfx3_eval_quad_node_0_l(fEdge); 
    fUpwindQuad[1] = ser_5x_p2_surfx3_eval_quad_node_1_l(fEdge); 
    fUpwindQuad[2] = ser_5x_p2_surfx3_eval_quad_node_2_l(fEdge); 
    fUpwindQuad[3] = ser_5x_p2_surfx3_eval_quad_node_3_l(fEdge); 
    fUpwindQuad[4] = ser_5x_p2_surfx3_eval_quad_node_4_l(fEdge); 
    fUpwindQuad[5] = ser_5x_p2_surfx3_eval_quad_node_5_l(fEdge); 
    fUpwindQuad[6] = ser_5x_p2_surfx3_eval_quad_node_6_l(fEdge); 
    fUpwindQuad[7] = ser_5x_p2_surfx3_eval_quad_node_7_l(fEdge); 
    fUpwindQuad[8] = ser_5x_p2_surfx3_eval_quad_node_8_l(fEdge); 
  } 
  if (0.223606797749979*alpha[12]+0.45*alpha[5]-0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[9] = ser_5x_p2_surfx3_eval_quad_node_9_r(fSkin); 
    fUpwindQuad[10] = ser_5x_p2_surfx3_eval_quad_node_10_r(fSkin); 
    fUpwindQuad[11] = ser_5x_p2_surfx3_eval_quad_node_11_r(fSkin); 
    fUpwindQuad[12] = ser_5x_p2_surfx3_eval_quad_node_12_r(fSkin); 
    fUpwindQuad[13] = ser_5x_p2_surfx3_eval_quad_node_13_r(fSkin); 
    fUpwindQuad[14] = ser_5x_p2_surfx3_eval_quad_node_14_r(fSkin); 
    fUpwindQuad[15] = ser_5x_p2_surfx3_eval_quad_node_15_r(fSkin); 
    fUpwindQuad[16] = ser_5x_p2_surfx3_eval_quad_node_16_r(fSkin); 
    fUpwindQuad[17] = ser_5x_p2_surfx3_eval_quad_node_17_r(fSkin); 
  } else { 
    fUpwindQuad[9] = ser_5x_p2_surfx3_eval_quad_node_9_l(fEdge); 
    fUpwindQuad[10] = ser_5x_p2_surfx3_eval_quad_node_10_l(fEdge); 
    fUpwindQuad[11] = ser_5x_p2_surfx3_eval_quad_node_11_l(fEdge); 
    fUpwindQuad[12] = ser_5x_p2_surfx3_eval_quad_node_12_l(fEdge); 
    fUpwindQuad[13] = ser_5x_p2_surfx3_eval_quad_node_13_l(fEdge); 
    fUpwindQuad[14] = ser_5x_p2_surfx3_eval_quad_node_14_l(fEdge); 
    fUpwindQuad[15] = ser_5x_p2_surfx3_eval_quad_node_15_l(fEdge); 
    fUpwindQuad[16] = ser_5x_p2_surfx3_eval_quad_node_16_l(fEdge); 
    fUpwindQuad[17] = ser_5x_p2_surfx3_eval_quad_node_17_l(fEdge); 
  } 
  if (0.223606797749979*alpha[12]+0.45*alpha[5]-0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[18] = ser_5x_p2_surfx3_eval_quad_node_18_r(fSkin); 
    fUpwindQuad[19] = ser_5x_p2_surfx3_eval_quad_node_19_r(fSkin); 
    fUpwindQuad[20] = ser_5x_p2_surfx3_eval_quad_node_20_r(fSkin); 
    fUpwindQuad[21] = ser_5x_p2_surfx3_eval_quad_node_21_r(fSkin); 
    fUpwindQuad[22] = ser_5x_p2_surfx3_eval_quad_node_22_r(fSkin); 
    fUpwindQuad[23] = ser_5x_p2_surfx3_eval_quad_node_23_r(fSkin); 
    fUpwindQuad[24] = ser_5x_p2_surfx3_eval_quad_node_24_r(fSkin); 
    fUpwindQuad[25] = ser_5x_p2_surfx3_eval_quad_node_25_r(fSkin); 
    fUpwindQuad[26] = ser_5x_p2_surfx3_eval_quad_node_26_r(fSkin); 
  } else { 
    fUpwindQuad[18] = ser_5x_p2_surfx3_eval_quad_node_18_l(fEdge); 
    fUpwindQuad[19] = ser_5x_p2_surfx3_eval_quad_node_19_l(fEdge); 
    fUpwindQuad[20] = ser_5x_p2_surfx3_eval_quad_node_20_l(fEdge); 
    fUpwindQuad[21] = ser_5x_p2_surfx3_eval_quad_node_21_l(fEdge); 
    fUpwindQuad[22] = ser_5x_p2_surfx3_eval_quad_node_22_l(fEdge); 
    fUpwindQuad[23] = ser_5x_p2_surfx3_eval_quad_node_23_l(fEdge); 
    fUpwindQuad[24] = ser_5x_p2_surfx3_eval_quad_node_24_l(fEdge); 
    fUpwindQuad[25] = ser_5x_p2_surfx3_eval_quad_node_25_l(fEdge); 
    fUpwindQuad[26] = ser_5x_p2_surfx3_eval_quad_node_26_l(fEdge); 
  } 
  if (0.223606797749979*alpha[12]+0.45*alpha[5]-0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[27] = ser_5x_p2_surfx3_eval_quad_node_27_r(fSkin); 
    fUpwindQuad[28] = ser_5x_p2_surfx3_eval_quad_node_28_r(fSkin); 
    fUpwindQuad[29] = ser_5x_p2_surfx3_eval_quad_node_29_r(fSkin); 
    fUpwindQuad[30] = ser_5x_p2_surfx3_eval_quad_node_30_r(fSkin); 
    fUpwindQuad[31] = ser_5x_p2_surfx3_eval_quad_node_31_r(fSkin); 
    fUpwindQuad[32] = ser_5x_p2_surfx3_eval_quad_node_32_r(fSkin); 
    fUpwindQuad[33] = ser_5x_p2_surfx3_eval_quad_node_33_r(fSkin); 
    fUpwindQuad[34] = ser_5x_p2_surfx3_eval_quad_node_34_r(fSkin); 
    fUpwindQuad[35] = ser_5x_p2_surfx3_eval_quad_node_35_r(fSkin); 
  } else { 
    fUpwindQuad[27] = ser_5x_p2_surfx3_eval_quad_node_27_l(fEdge); 
    fUpwindQuad[28] = ser_5x_p2_surfx3_eval_quad_node_28_l(fEdge); 
    fUpwindQuad[29] = ser_5x_p2_surfx3_eval_quad_node_29_l(fEdge); 
    fUpwindQuad[30] = ser_5x_p2_surfx3_eval_quad_node_30_l(fEdge); 
    fUpwindQuad[31] = ser_5x_p2_surfx3_eval_quad_node_31_l(fEdge); 
    fUpwindQuad[32] = ser_5x_p2_surfx3_eval_quad_node_32_l(fEdge); 
    fUpwindQuad[33] = ser_5x_p2_surfx3_eval_quad_node_33_l(fEdge); 
    fUpwindQuad[34] = ser_5x_p2_surfx3_eval_quad_node_34_l(fEdge); 
    fUpwindQuad[35] = ser_5x_p2_surfx3_eval_quad_node_35_l(fEdge); 
  } 
  if (0.223606797749979*alpha[12]+0.45*alpha[5]-0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[36] = ser_5x_p2_surfx3_eval_quad_node_36_r(fSkin); 
    fUpwindQuad[37] = ser_5x_p2_surfx3_eval_quad_node_37_r(fSkin); 
    fUpwindQuad[38] = ser_5x_p2_surfx3_eval_quad_node_38_r(fSkin); 
    fUpwindQuad[39] = ser_5x_p2_surfx3_eval_quad_node_39_r(fSkin); 
    fUpwindQuad[40] = ser_5x_p2_surfx3_eval_quad_node_40_r(fSkin); 
    fUpwindQuad[41] = ser_5x_p2_surfx3_eval_quad_node_41_r(fSkin); 
    fUpwindQuad[42] = ser_5x_p2_surfx3_eval_quad_node_42_r(fSkin); 
    fUpwindQuad[43] = ser_5x_p2_surfx3_eval_quad_node_43_r(fSkin); 
    fUpwindQuad[44] = ser_5x_p2_surfx3_eval_quad_node_44_r(fSkin); 
  } else { 
    fUpwindQuad[36] = ser_5x_p2_surfx3_eval_quad_node_36_l(fEdge); 
    fUpwindQuad[37] = ser_5x_p2_surfx3_eval_quad_node_37_l(fEdge); 
    fUpwindQuad[38] = ser_5x_p2_surfx3_eval_quad_node_38_l(fEdge); 
    fUpwindQuad[39] = ser_5x_p2_surfx3_eval_quad_node_39_l(fEdge); 
    fUpwindQuad[40] = ser_5x_p2_surfx3_eval_quad_node_40_l(fEdge); 
    fUpwindQuad[41] = ser_5x_p2_surfx3_eval_quad_node_41_l(fEdge); 
    fUpwindQuad[42] = ser_5x_p2_surfx3_eval_quad_node_42_l(fEdge); 
    fUpwindQuad[43] = ser_5x_p2_surfx3_eval_quad_node_43_l(fEdge); 
    fUpwindQuad[44] = ser_5x_p2_surfx3_eval_quad_node_44_l(fEdge); 
  } 
  if (0.223606797749979*alpha[12]+0.45*alpha[5]-0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[45] = ser_5x_p2_surfx3_eval_quad_node_45_r(fSkin); 
    fUpwindQuad[46] = ser_5x_p2_surfx3_eval_quad_node_46_r(fSkin); 
    fUpwindQuad[47] = ser_5x_p2_surfx3_eval_quad_node_47_r(fSkin); 
    fUpwindQuad[48] = ser_5x_p2_surfx3_eval_quad_node_48_r(fSkin); 
    fUpwindQuad[49] = ser_5x_p2_surfx3_eval_quad_node_49_r(fSkin); 
    fUpwindQuad[50] = ser_5x_p2_surfx3_eval_quad_node_50_r(fSkin); 
    fUpwindQuad[51] = ser_5x_p2_surfx3_eval_quad_node_51_r(fSkin); 
    fUpwindQuad[52] = ser_5x_p2_surfx3_eval_quad_node_52_r(fSkin); 
    fUpwindQuad[53] = ser_5x_p2_surfx3_eval_quad_node_53_r(fSkin); 
  } else { 
    fUpwindQuad[45] = ser_5x_p2_surfx3_eval_quad_node_45_l(fEdge); 
    fUpwindQuad[46] = ser_5x_p2_surfx3_eval_quad_node_46_l(fEdge); 
    fUpwindQuad[47] = ser_5x_p2_surfx3_eval_quad_node_47_l(fEdge); 
    fUpwindQuad[48] = ser_5x_p2_surfx3_eval_quad_node_48_l(fEdge); 
    fUpwindQuad[49] = ser_5x_p2_surfx3_eval_quad_node_49_l(fEdge); 
    fUpwindQuad[50] = ser_5x_p2_surfx3_eval_quad_node_50_l(fEdge); 
    fUpwindQuad[51] = ser_5x_p2_surfx3_eval_quad_node_51_l(fEdge); 
    fUpwindQuad[52] = ser_5x_p2_surfx3_eval_quad_node_52_l(fEdge); 
    fUpwindQuad[53] = ser_5x_p2_surfx3_eval_quad_node_53_l(fEdge); 
  } 
  if (0.223606797749979*alpha[12]+0.45*alpha[5]-0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[54] = ser_5x_p2_surfx3_eval_quad_node_54_r(fSkin); 
    fUpwindQuad[55] = ser_5x_p2_surfx3_eval_quad_node_55_r(fSkin); 
    fUpwindQuad[56] = ser_5x_p2_surfx3_eval_quad_node_56_r(fSkin); 
    fUpwindQuad[57] = ser_5x_p2_surfx3_eval_quad_node_57_r(fSkin); 
    fUpwindQuad[58] = ser_5x_p2_surfx3_eval_quad_node_58_r(fSkin); 
    fUpwindQuad[59] = ser_5x_p2_surfx3_eval_quad_node_59_r(fSkin); 
    fUpwindQuad[60] = ser_5x_p2_surfx3_eval_quad_node_60_r(fSkin); 
    fUpwindQuad[61] = ser_5x_p2_surfx3_eval_quad_node_61_r(fSkin); 
    fUpwindQuad[62] = ser_5x_p2_surfx3_eval_quad_node_62_r(fSkin); 
  } else { 
    fUpwindQuad[54] = ser_5x_p2_surfx3_eval_quad_node_54_l(fEdge); 
    fUpwindQuad[55] = ser_5x_p2_surfx3_eval_quad_node_55_l(fEdge); 
    fUpwindQuad[56] = ser_5x_p2_surfx3_eval_quad_node_56_l(fEdge); 
    fUpwindQuad[57] = ser_5x_p2_surfx3_eval_quad_node_57_l(fEdge); 
    fUpwindQuad[58] = ser_5x_p2_surfx3_eval_quad_node_58_l(fEdge); 
    fUpwindQuad[59] = ser_5x_p2_surfx3_eval_quad_node_59_l(fEdge); 
    fUpwindQuad[60] = ser_5x_p2_surfx3_eval_quad_node_60_l(fEdge); 
    fUpwindQuad[61] = ser_5x_p2_surfx3_eval_quad_node_61_l(fEdge); 
    fUpwindQuad[62] = ser_5x_p2_surfx3_eval_quad_node_62_l(fEdge); 
  } 
  if (0.223606797749979*alpha[12]+0.45*alpha[5]-0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[63] = ser_5x_p2_surfx3_eval_quad_node_63_r(fSkin); 
    fUpwindQuad[64] = ser_5x_p2_surfx3_eval_quad_node_64_r(fSkin); 
    fUpwindQuad[65] = ser_5x_p2_surfx3_eval_quad_node_65_r(fSkin); 
    fUpwindQuad[66] = ser_5x_p2_surfx3_eval_quad_node_66_r(fSkin); 
    fUpwindQuad[67] = ser_5x_p2_surfx3_eval_quad_node_67_r(fSkin); 
    fUpwindQuad[68] = ser_5x_p2_surfx3_eval_quad_node_68_r(fSkin); 
    fUpwindQuad[69] = ser_5x_p2_surfx3_eval_quad_node_69_r(fSkin); 
    fUpwindQuad[70] = ser_5x_p2_surfx3_eval_quad_node_70_r(fSkin); 
    fUpwindQuad[71] = ser_5x_p2_surfx3_eval_quad_node_71_r(fSkin); 
  } else { 
    fUpwindQuad[63] = ser_5x_p2_surfx3_eval_quad_node_63_l(fEdge); 
    fUpwindQuad[64] = ser_5x_p2_surfx3_eval_quad_node_64_l(fEdge); 
    fUpwindQuad[65] = ser_5x_p2_surfx3_eval_quad_node_65_l(fEdge); 
    fUpwindQuad[66] = ser_5x_p2_surfx3_eval_quad_node_66_l(fEdge); 
    fUpwindQuad[67] = ser_5x_p2_surfx3_eval_quad_node_67_l(fEdge); 
    fUpwindQuad[68] = ser_5x_p2_surfx3_eval_quad_node_68_l(fEdge); 
    fUpwindQuad[69] = ser_5x_p2_surfx3_eval_quad_node_69_l(fEdge); 
    fUpwindQuad[70] = ser_5x_p2_surfx3_eval_quad_node_70_l(fEdge); 
    fUpwindQuad[71] = ser_5x_p2_surfx3_eval_quad_node_71_l(fEdge); 
  } 
  if (0.223606797749979*alpha[12]+0.45*alpha[5]-0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[72] = ser_5x_p2_surfx3_eval_quad_node_72_r(fSkin); 
    fUpwindQuad[73] = ser_5x_p2_surfx3_eval_quad_node_73_r(fSkin); 
    fUpwindQuad[74] = ser_5x_p2_surfx3_eval_quad_node_74_r(fSkin); 
    fUpwindQuad[75] = ser_5x_p2_surfx3_eval_quad_node_75_r(fSkin); 
    fUpwindQuad[76] = ser_5x_p2_surfx3_eval_quad_node_76_r(fSkin); 
    fUpwindQuad[77] = ser_5x_p2_surfx3_eval_quad_node_77_r(fSkin); 
    fUpwindQuad[78] = ser_5x_p2_surfx3_eval_quad_node_78_r(fSkin); 
    fUpwindQuad[79] = ser_5x_p2_surfx3_eval_quad_node_79_r(fSkin); 
    fUpwindQuad[80] = ser_5x_p2_surfx3_eval_quad_node_80_r(fSkin); 
  } else { 
    fUpwindQuad[72] = ser_5x_p2_surfx3_eval_quad_node_72_l(fEdge); 
    fUpwindQuad[73] = ser_5x_p2_surfx3_eval_quad_node_73_l(fEdge); 
    fUpwindQuad[74] = ser_5x_p2_surfx3_eval_quad_node_74_l(fEdge); 
    fUpwindQuad[75] = ser_5x_p2_surfx3_eval_quad_node_75_l(fEdge); 
    fUpwindQuad[76] = ser_5x_p2_surfx3_eval_quad_node_76_l(fEdge); 
    fUpwindQuad[77] = ser_5x_p2_surfx3_eval_quad_node_77_l(fEdge); 
    fUpwindQuad[78] = ser_5x_p2_surfx3_eval_quad_node_78_l(fEdge); 
    fUpwindQuad[79] = ser_5x_p2_surfx3_eval_quad_node_79_l(fEdge); 
    fUpwindQuad[80] = ser_5x_p2_surfx3_eval_quad_node_80_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_5x_p2_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.25*alpha[12]*fUpwind[12]+0.25*alpha[5]*fUpwind[5]+0.25*alpha[2]*fUpwind[2]+0.25*alpha[1]*fUpwind[1]+0.25*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.2500000000000001*alpha[12]*fUpwind[20]+0.223606797749979*alpha[5]*fUpwind[19]+0.223606797749979*alpha[1]*fUpwind[11]+0.25*alpha[2]*fUpwind[5]+0.25*fUpwind[2]*alpha[5]+0.25*alpha[0]*fUpwind[1]+0.25*fUpwind[0]*alpha[1]; 
  Ghat[2] = 0.223606797749979*alpha[5]*fUpwind[20]+0.223606797749979*alpha[2]*fUpwind[12]+0.223606797749979*fUpwind[2]*alpha[12]+0.25*alpha[1]*fUpwind[5]+0.25*fUpwind[1]*alpha[5]+0.25*alpha[0]*fUpwind[2]+0.25*fUpwind[0]*alpha[2]; 
  Ghat[3] = 0.2500000000000001*alpha[12]*fUpwind[22]+0.25*alpha[5]*fUpwind[15]+0.25*alpha[2]*fUpwind[7]+0.25*alpha[1]*fUpwind[6]+0.25*alpha[0]*fUpwind[3]; 
  Ghat[4] = 0.2500000000000001*alpha[12]*fUpwind[26]+0.25*alpha[5]*fUpwind[16]+0.25*alpha[2]*fUpwind[9]+0.25*alpha[1]*fUpwind[8]+0.25*alpha[0]*fUpwind[4]; 
  Ghat[5] = 0.223606797749979*alpha[2]*fUpwind[20]+0.223606797749979*alpha[1]*fUpwind[19]+0.223606797749979*alpha[5]*fUpwind[12]+0.223606797749979*fUpwind[5]*alpha[12]+0.223606797749979*alpha[5]*fUpwind[11]+0.25*alpha[0]*fUpwind[5]+0.25*fUpwind[0]*alpha[5]+0.25*alpha[1]*fUpwind[2]+0.25*fUpwind[1]*alpha[2]; 
  Ghat[6] = 0.25*alpha[12]*fUpwind[33]+0.223606797749979*alpha[5]*fUpwind[32]+0.223606797749979*alpha[1]*fUpwind[21]+0.25*alpha[2]*fUpwind[15]+0.25*alpha[5]*fUpwind[7]+0.25*alpha[0]*fUpwind[6]+0.25*alpha[1]*fUpwind[3]; 
  Ghat[7] = 0.223606797749979*alpha[5]*fUpwind[33]+0.223606797749979*alpha[2]*fUpwind[22]+0.25*alpha[1]*fUpwind[15]+0.223606797749979*fUpwind[7]*alpha[12]+0.25*alpha[0]*fUpwind[7]+0.25*alpha[5]*fUpwind[6]+0.25*alpha[2]*fUpwind[3]; 
  Ghat[8] = 0.25*alpha[12]*fUpwind[36]+0.223606797749979*alpha[5]*fUpwind[35]+0.223606797749979*alpha[1]*fUpwind[25]+0.25*alpha[2]*fUpwind[16]+0.25*alpha[5]*fUpwind[9]+0.25*alpha[0]*fUpwind[8]+0.25*alpha[1]*fUpwind[4]; 
  Ghat[9] = 0.223606797749979*alpha[5]*fUpwind[36]+0.223606797749979*alpha[2]*fUpwind[26]+0.25*alpha[1]*fUpwind[16]+0.223606797749979*fUpwind[9]*alpha[12]+0.25*alpha[0]*fUpwind[9]+0.25*alpha[5]*fUpwind[8]+0.25*alpha[2]*fUpwind[4]; 
  Ghat[10] = 0.25*alpha[12]*fUpwind[38]+0.25*alpha[5]*fUpwind[31]+0.25*alpha[2]*fUpwind[18]+0.25*alpha[1]*fUpwind[17]+0.25*alpha[0]*fUpwind[10]; 
  Ghat[11] = 0.2500000000000001*alpha[2]*fUpwind[19]+0.25*alpha[0]*fUpwind[11]+0.223606797749979*alpha[5]*fUpwind[5]+0.223606797749979*alpha[1]*fUpwind[1]; 
  Ghat[12] = 0.2500000000000001*alpha[1]*fUpwind[20]+0.159719141249985*alpha[12]*fUpwind[12]+0.25*alpha[0]*fUpwind[12]+0.25*fUpwind[0]*alpha[12]+0.223606797749979*alpha[5]*fUpwind[5]+0.223606797749979*alpha[2]*fUpwind[2]; 
  Ghat[13] = 0.25*alpha[5]*fUpwind[34]+0.2500000000000001*alpha[2]*fUpwind[24]+0.2500000000000001*alpha[1]*fUpwind[23]+0.25*alpha[0]*fUpwind[13]; 
  Ghat[14] = 0.25*alpha[5]*fUpwind[41]+0.2500000000000001*alpha[2]*fUpwind[29]+0.2500000000000001*alpha[1]*fUpwind[28]+0.25*alpha[0]*fUpwind[14]; 
  Ghat[15] = 0.223606797749979*alpha[2]*fUpwind[33]+0.223606797749979*alpha[1]*fUpwind[32]+0.223606797749979*alpha[5]*fUpwind[22]+0.223606797749979*alpha[5]*fUpwind[21]+0.223606797749979*alpha[12]*fUpwind[15]+0.25*alpha[0]*fUpwind[15]+0.25*alpha[1]*fUpwind[7]+0.25*alpha[2]*fUpwind[6]+0.25*fUpwind[3]*alpha[5]; 
  Ghat[16] = 0.223606797749979*alpha[2]*fUpwind[36]+0.223606797749979*alpha[1]*fUpwind[35]+0.223606797749979*alpha[5]*fUpwind[26]+0.223606797749979*alpha[5]*fUpwind[25]+0.223606797749979*alpha[12]*fUpwind[16]+0.25*alpha[0]*fUpwind[16]+0.25*alpha[1]*fUpwind[9]+0.25*alpha[2]*fUpwind[8]+0.25*fUpwind[4]*alpha[5]; 
  Ghat[17] = 0.2500000000000001*alpha[12]*fUpwind[45]+0.223606797749979*alpha[5]*fUpwind[44]+0.223606797749979*alpha[1]*fUpwind[37]+0.25*alpha[2]*fUpwind[31]+0.25*alpha[5]*fUpwind[18]+0.25*alpha[0]*fUpwind[17]+0.25*alpha[1]*fUpwind[10]; 
  Ghat[18] = 0.223606797749979*alpha[5]*fUpwind[45]+0.223606797749979*alpha[2]*fUpwind[38]+0.25*alpha[1]*fUpwind[31]+0.223606797749979*alpha[12]*fUpwind[18]+0.25*alpha[0]*fUpwind[18]+0.25*alpha[5]*fUpwind[17]+0.25*alpha[2]*fUpwind[10]; 
  Ghat[19] = 0.2*alpha[5]*fUpwind[20]+0.223606797749979*alpha[12]*fUpwind[19]+0.25*alpha[0]*fUpwind[19]+0.2500000000000001*alpha[2]*fUpwind[11]+0.223606797749979*alpha[1]*fUpwind[5]+0.223606797749979*fUpwind[1]*alpha[5]; 
  Ghat[20] = 0.159719141249985*alpha[12]*fUpwind[20]+0.25*alpha[0]*fUpwind[20]+0.2*alpha[5]*fUpwind[19]+0.2500000000000001*alpha[1]*fUpwind[12]+0.2500000000000001*fUpwind[1]*alpha[12]+0.223606797749979*alpha[2]*fUpwind[5]+0.223606797749979*fUpwind[2]*alpha[5]; 
  Ghat[21] = 0.2500000000000001*alpha[2]*fUpwind[32]+0.25*alpha[0]*fUpwind[21]+0.223606797749979*alpha[5]*fUpwind[15]+0.223606797749979*alpha[1]*fUpwind[6]; 
  Ghat[22] = 0.2500000000000001*alpha[1]*fUpwind[33]+0.159719141249985*alpha[12]*fUpwind[22]+0.25*alpha[0]*fUpwind[22]+0.223606797749979*alpha[5]*fUpwind[15]+0.2500000000000001*fUpwind[3]*alpha[12]+0.223606797749979*alpha[2]*fUpwind[7]; 
  Ghat[23] = 0.2500000000000001*alpha[2]*fUpwind[34]+0.25*alpha[5]*fUpwind[24]+0.25*alpha[0]*fUpwind[23]+0.2500000000000001*alpha[1]*fUpwind[13]; 
  Ghat[24] = 0.2500000000000001*alpha[1]*fUpwind[34]+0.223606797749979*alpha[12]*fUpwind[24]+0.25*alpha[0]*fUpwind[24]+0.25*alpha[5]*fUpwind[23]+0.2500000000000001*alpha[2]*fUpwind[13]; 
  Ghat[25] = 0.2500000000000001*alpha[2]*fUpwind[35]+0.25*alpha[0]*fUpwind[25]+0.223606797749979*alpha[5]*fUpwind[16]+0.223606797749979*alpha[1]*fUpwind[8]; 
  Ghat[26] = 0.2500000000000001*alpha[1]*fUpwind[36]+0.159719141249985*alpha[12]*fUpwind[26]+0.25*alpha[0]*fUpwind[26]+0.223606797749979*alpha[5]*fUpwind[16]+0.2500000000000001*fUpwind[4]*alpha[12]+0.223606797749979*alpha[2]*fUpwind[9]; 
  Ghat[27] = 0.25*alpha[5]*fUpwind[46]+0.2500000000000001*alpha[2]*fUpwind[40]+0.2500000000000001*alpha[1]*fUpwind[39]+0.25*alpha[0]*fUpwind[27]; 
  Ghat[28] = 0.2500000000000001*alpha[2]*fUpwind[41]+0.25*alpha[5]*fUpwind[29]+0.25*alpha[0]*fUpwind[28]+0.2500000000000001*alpha[1]*fUpwind[14]; 
  Ghat[29] = 0.2500000000000001*alpha[1]*fUpwind[41]+0.223606797749979*alpha[12]*fUpwind[29]+0.25*alpha[0]*fUpwind[29]+0.25*alpha[5]*fUpwind[28]+0.2500000000000001*alpha[2]*fUpwind[14]; 
  Ghat[30] = 0.25*alpha[5]*fUpwind[47]+0.2500000000000001*alpha[2]*fUpwind[43]+0.2500000000000001*alpha[1]*fUpwind[42]+0.25*alpha[0]*fUpwind[30]; 
  Ghat[31] = 0.223606797749979*alpha[2]*fUpwind[45]+0.223606797749979*alpha[1]*fUpwind[44]+0.223606797749979*alpha[5]*fUpwind[38]+0.223606797749979*alpha[5]*fUpwind[37]+0.223606797749979*alpha[12]*fUpwind[31]+0.25*alpha[0]*fUpwind[31]+0.25*alpha[1]*fUpwind[18]+0.25*alpha[2]*fUpwind[17]+0.25*alpha[5]*fUpwind[10]; 
  Ghat[32] = 0.2*alpha[5]*fUpwind[33]+0.223606797749979*alpha[12]*fUpwind[32]+0.25*alpha[0]*fUpwind[32]+0.2500000000000001*alpha[2]*fUpwind[21]+0.223606797749979*alpha[1]*fUpwind[15]+0.223606797749979*alpha[5]*fUpwind[6]; 
  Ghat[33] = 0.159719141249985*alpha[12]*fUpwind[33]+0.25*alpha[0]*fUpwind[33]+0.2*alpha[5]*fUpwind[32]+0.2500000000000001*alpha[1]*fUpwind[22]+0.223606797749979*alpha[2]*fUpwind[15]+0.25*fUpwind[6]*alpha[12]+0.223606797749979*alpha[5]*fUpwind[7]; 
  Ghat[34] = 0.223606797749979*alpha[12]*fUpwind[34]+0.25*alpha[0]*fUpwind[34]+0.2500000000000001*alpha[1]*fUpwind[24]+0.2500000000000001*alpha[2]*fUpwind[23]+0.25*alpha[5]*fUpwind[13]; 
  Ghat[35] = 0.2*alpha[5]*fUpwind[36]+0.223606797749979*alpha[12]*fUpwind[35]+0.25*alpha[0]*fUpwind[35]+0.2500000000000001*alpha[2]*fUpwind[25]+0.223606797749979*alpha[1]*fUpwind[16]+0.223606797749979*alpha[5]*fUpwind[8]; 
  Ghat[36] = 0.159719141249985*alpha[12]*fUpwind[36]+0.25*alpha[0]*fUpwind[36]+0.2*alpha[5]*fUpwind[35]+0.2500000000000001*alpha[1]*fUpwind[26]+0.223606797749979*alpha[2]*fUpwind[16]+0.25*fUpwind[8]*alpha[12]+0.223606797749979*alpha[5]*fUpwind[9]; 
  Ghat[37] = 0.2500000000000001*alpha[2]*fUpwind[44]+0.25*alpha[0]*fUpwind[37]+0.223606797749979*alpha[5]*fUpwind[31]+0.223606797749979*alpha[1]*fUpwind[17]; 
  Ghat[38] = 0.2500000000000001*alpha[1]*fUpwind[45]+0.159719141249985*alpha[12]*fUpwind[38]+0.25*alpha[0]*fUpwind[38]+0.223606797749979*alpha[5]*fUpwind[31]+0.223606797749979*alpha[2]*fUpwind[18]+0.25*fUpwind[10]*alpha[12]; 
  Ghat[39] = 0.2500000000000001*alpha[2]*fUpwind[46]+0.25*alpha[5]*fUpwind[40]+0.25*alpha[0]*fUpwind[39]+0.2500000000000001*alpha[1]*fUpwind[27]; 
  Ghat[40] = 0.2500000000000001*alpha[1]*fUpwind[46]+0.223606797749979*alpha[12]*fUpwind[40]+0.25*alpha[0]*fUpwind[40]+0.25*alpha[5]*fUpwind[39]+0.2500000000000001*alpha[2]*fUpwind[27]; 
  Ghat[41] = 0.223606797749979*alpha[12]*fUpwind[41]+0.25*alpha[0]*fUpwind[41]+0.2500000000000001*alpha[1]*fUpwind[29]+0.2500000000000001*alpha[2]*fUpwind[28]+0.25*alpha[5]*fUpwind[14]; 
  Ghat[42] = 0.2500000000000001*alpha[2]*fUpwind[47]+0.25*alpha[5]*fUpwind[43]+0.25*alpha[0]*fUpwind[42]+0.2500000000000001*alpha[1]*fUpwind[30]; 
  Ghat[43] = 0.2500000000000001*alpha[1]*fUpwind[47]+0.223606797749979*alpha[12]*fUpwind[43]+0.25*alpha[0]*fUpwind[43]+0.25*alpha[5]*fUpwind[42]+0.2500000000000001*alpha[2]*fUpwind[30]; 
  Ghat[44] = 0.2*alpha[5]*fUpwind[45]+0.223606797749979*alpha[12]*fUpwind[44]+0.25*alpha[0]*fUpwind[44]+0.2500000000000001*alpha[2]*fUpwind[37]+0.223606797749979*alpha[1]*fUpwind[31]+0.223606797749979*alpha[5]*fUpwind[17]; 
  Ghat[45] = 0.159719141249985*alpha[12]*fUpwind[45]+0.25*alpha[0]*fUpwind[45]+0.2*alpha[5]*fUpwind[44]+0.2500000000000001*alpha[1]*fUpwind[38]+0.223606797749979*alpha[2]*fUpwind[31]+0.223606797749979*alpha[5]*fUpwind[18]+0.2500000000000001*alpha[12]*fUpwind[17]; 
  Ghat[46] = 0.223606797749979*alpha[12]*fUpwind[46]+0.25*alpha[0]*fUpwind[46]+0.2500000000000001*alpha[1]*fUpwind[40]+0.2500000000000001*alpha[2]*fUpwind[39]+0.25*alpha[5]*fUpwind[27]; 
  Ghat[47] = 0.223606797749979*alpha[12]*fUpwind[47]+0.25*alpha[0]*fUpwind[47]+0.2500000000000001*alpha[1]*fUpwind[43]+0.2500000000000001*alpha[2]*fUpwind[42]+0.25*alpha[5]*fUpwind[30]; 

  out[0] += -0.7071067811865475*Ghat[0]*dv10; 
  out[1] += -0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -0.7071067811865475*Ghat[2]*dv10; 
  out[3] += -1.224744871391589*Ghat[0]*dv10; 
  out[4] += -0.7071067811865475*Ghat[3]*dv10; 
  out[5] += -0.7071067811865475*Ghat[4]*dv10; 
  out[6] += -0.7071067811865475*Ghat[5]*dv10; 
  out[7] += -1.224744871391589*Ghat[1]*dv10; 
  out[8] += -1.224744871391589*Ghat[2]*dv10; 
  out[9] += -0.7071067811865475*Ghat[6]*dv10; 
  out[10] += -0.7071067811865475*Ghat[7]*dv10; 
  out[11] += -1.224744871391589*Ghat[3]*dv10; 
  out[12] += -0.7071067811865475*Ghat[8]*dv10; 
  out[13] += -0.7071067811865475*Ghat[9]*dv10; 
  out[14] += -1.224744871391589*Ghat[4]*dv10; 
  out[15] += -0.7071067811865475*Ghat[10]*dv10; 
  out[16] += -0.7071067811865475*Ghat[11]*dv10; 
  out[17] += -0.7071067811865475*Ghat[12]*dv10; 
  out[18] += -1.58113883008419*Ghat[0]*dv10; 
  out[19] += -0.7071067811865475*Ghat[13]*dv10; 
  out[20] += -0.7071067811865475*Ghat[14]*dv10; 
  out[21] += -1.224744871391589*Ghat[5]*dv10; 
  out[22] += -0.7071067811865475*Ghat[15]*dv10; 
  out[23] += -1.224744871391589*Ghat[6]*dv10; 
  out[24] += -1.224744871391589*Ghat[7]*dv10; 
  out[25] += -0.7071067811865475*Ghat[16]*dv10; 
  out[26] += -1.224744871391589*Ghat[8]*dv10; 
  out[27] += -1.224744871391589*Ghat[9]*dv10; 
  out[28] += -0.7071067811865475*Ghat[17]*dv10; 
  out[29] += -0.7071067811865475*Ghat[18]*dv10; 
  out[30] += -1.224744871391589*Ghat[10]*dv10; 
  out[31] += -0.7071067811865475*Ghat[19]*dv10; 
  out[32] += -0.7071067811865475*Ghat[20]*dv10; 
  out[33] += -1.224744871391589*Ghat[11]*dv10; 
  out[34] += -1.224744871391589*Ghat[12]*dv10; 
  out[35] += -1.58113883008419*Ghat[1]*dv10; 
  out[36] += -1.58113883008419*Ghat[2]*dv10; 
  out[37] += -0.7071067811865475*Ghat[21]*dv10; 
  out[38] += -0.7071067811865475*Ghat[22]*dv10; 
  out[39] += -1.58113883008419*Ghat[3]*dv10; 
  out[40] += -0.7071067811865475*Ghat[23]*dv10; 
  out[41] += -0.7071067811865475*Ghat[24]*dv10; 
  out[42] += -1.224744871391589*Ghat[13]*dv10; 
  out[43] += -0.7071067811865475*Ghat[25]*dv10; 
  out[44] += -0.7071067811865475*Ghat[26]*dv10; 
  out[45] += -1.58113883008419*Ghat[4]*dv10; 
  out[46] += -0.7071067811865475*Ghat[27]*dv10; 
  out[47] += -0.7071067811865475*Ghat[28]*dv10; 
  out[48] += -0.7071067811865475*Ghat[29]*dv10; 
  out[49] += -1.224744871391589*Ghat[14]*dv10; 
  out[50] += -0.7071067811865475*Ghat[30]*dv10; 
  out[51] += -1.224744871391589*Ghat[15]*dv10; 
  out[52] += -1.224744871391589*Ghat[16]*dv10; 
  out[53] += -0.7071067811865475*Ghat[31]*dv10; 
  out[54] += -1.224744871391589*Ghat[17]*dv10; 
  out[55] += -1.224744871391589*Ghat[18]*dv10; 
  out[56] += -1.224744871391589*Ghat[19]*dv10; 
  out[57] += -1.224744871391589*Ghat[20]*dv10; 
  out[58] += -1.58113883008419*Ghat[5]*dv10; 
  out[59] += -0.7071067811865475*Ghat[32]*dv10; 
  out[60] += -0.7071067811865475*Ghat[33]*dv10; 
  out[61] += -1.224744871391589*Ghat[21]*dv10; 
  out[62] += -1.224744871391589*Ghat[22]*dv10; 
  out[63] += -1.58113883008419*Ghat[6]*dv10; 
  out[64] += -1.58113883008419*Ghat[7]*dv10; 
  out[65] += -0.7071067811865475*Ghat[34]*dv10; 
  out[66] += -1.224744871391589*Ghat[23]*dv10; 
  out[67] += -1.224744871391589*Ghat[24]*dv10; 
  out[68] += -0.7071067811865475*Ghat[35]*dv10; 
  out[69] += -0.7071067811865475*Ghat[36]*dv10; 
  out[70] += -1.224744871391589*Ghat[25]*dv10; 
  out[71] += -1.224744871391589*Ghat[26]*dv10; 
  out[72] += -1.58113883008419*Ghat[8]*dv10; 
  out[73] += -1.58113883008419*Ghat[9]*dv10; 
  out[74] += -0.7071067811865475*Ghat[37]*dv10; 
  out[75] += -0.7071067811865475*Ghat[38]*dv10; 
  out[76] += -1.58113883008419*Ghat[10]*dv10; 
  out[77] += -0.7071067811865475*Ghat[39]*dv10; 
  out[78] += -0.7071067811865475*Ghat[40]*dv10; 
  out[79] += -1.224744871391589*Ghat[27]*dv10; 
  out[80] += -0.7071067811865475*Ghat[41]*dv10; 
  out[81] += -1.224744871391589*Ghat[28]*dv10; 
  out[82] += -1.224744871391589*Ghat[29]*dv10; 
  out[83] += -0.7071067811865475*Ghat[42]*dv10; 
  out[84] += -0.7071067811865475*Ghat[43]*dv10; 
  out[85] += -1.224744871391589*Ghat[30]*dv10; 
  out[86] += -1.224744871391589*Ghat[31]*dv10; 
  out[87] += -1.224744871391589*Ghat[32]*dv10; 
  out[88] += -1.224744871391589*Ghat[33]*dv10; 
  out[89] += -1.58113883008419*Ghat[15]*dv10; 
  out[90] += -1.224744871391589*Ghat[34]*dv10; 
  out[91] += -1.224744871391589*Ghat[35]*dv10; 
  out[92] += -1.224744871391589*Ghat[36]*dv10; 
  out[93] += -1.58113883008419*Ghat[16]*dv10; 
  out[94] += -0.7071067811865475*Ghat[44]*dv10; 
  out[95] += -0.7071067811865475*Ghat[45]*dv10; 
  out[96] += -1.224744871391589*Ghat[37]*dv10; 
  out[97] += -1.224744871391589*Ghat[38]*dv10; 
  out[98] += -1.58113883008419*Ghat[17]*dv10; 
  out[99] += -1.58113883008419*Ghat[18]*dv10; 
  out[100] += -0.7071067811865475*Ghat[46]*dv10; 
  out[101] += -1.224744871391589*Ghat[39]*dv10; 
  out[102] += -1.224744871391589*Ghat[40]*dv10; 
  out[103] += -1.224744871391589*Ghat[41]*dv10; 
  out[104] += -0.7071067811865475*Ghat[47]*dv10; 
  out[105] += -1.224744871391589*Ghat[42]*dv10; 
  out[106] += -1.224744871391589*Ghat[43]*dv10; 
  out[107] += -1.224744871391589*Ghat[44]*dv10; 
  out[108] += -1.224744871391589*Ghat[45]*dv10; 
  out[109] += -1.58113883008419*Ghat[31]*dv10; 
  out[110] += -1.224744871391589*Ghat[46]*dv10; 
  out[111] += -1.224744871391589*Ghat[47]*dv10; 

  } else { 

  if (0.223606797749979*alpha[12]+0.45*alpha[5]-0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[0] = ser_5x_p2_surfx3_eval_quad_node_0_r(fEdge); 
    fUpwindQuad[1] = ser_5x_p2_surfx3_eval_quad_node_1_r(fEdge); 
    fUpwindQuad[2] = ser_5x_p2_surfx3_eval_quad_node_2_r(fEdge); 
    fUpwindQuad[3] = ser_5x_p2_surfx3_eval_quad_node_3_r(fEdge); 
    fUpwindQuad[4] = ser_5x_p2_surfx3_eval_quad_node_4_r(fEdge); 
    fUpwindQuad[5] = ser_5x_p2_surfx3_eval_quad_node_5_r(fEdge); 
    fUpwindQuad[6] = ser_5x_p2_surfx3_eval_quad_node_6_r(fEdge); 
    fUpwindQuad[7] = ser_5x_p2_surfx3_eval_quad_node_7_r(fEdge); 
    fUpwindQuad[8] = ser_5x_p2_surfx3_eval_quad_node_8_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_5x_p2_surfx3_eval_quad_node_0_l(fSkin); 
    fUpwindQuad[1] = ser_5x_p2_surfx3_eval_quad_node_1_l(fSkin); 
    fUpwindQuad[2] = ser_5x_p2_surfx3_eval_quad_node_2_l(fSkin); 
    fUpwindQuad[3] = ser_5x_p2_surfx3_eval_quad_node_3_l(fSkin); 
    fUpwindQuad[4] = ser_5x_p2_surfx3_eval_quad_node_4_l(fSkin); 
    fUpwindQuad[5] = ser_5x_p2_surfx3_eval_quad_node_5_l(fSkin); 
    fUpwindQuad[6] = ser_5x_p2_surfx3_eval_quad_node_6_l(fSkin); 
    fUpwindQuad[7] = ser_5x_p2_surfx3_eval_quad_node_7_l(fSkin); 
    fUpwindQuad[8] = ser_5x_p2_surfx3_eval_quad_node_8_l(fSkin); 
  } 
  if (0.223606797749979*alpha[12]+0.45*alpha[5]-0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[9] = ser_5x_p2_surfx3_eval_quad_node_9_r(fEdge); 
    fUpwindQuad[10] = ser_5x_p2_surfx3_eval_quad_node_10_r(fEdge); 
    fUpwindQuad[11] = ser_5x_p2_surfx3_eval_quad_node_11_r(fEdge); 
    fUpwindQuad[12] = ser_5x_p2_surfx3_eval_quad_node_12_r(fEdge); 
    fUpwindQuad[13] = ser_5x_p2_surfx3_eval_quad_node_13_r(fEdge); 
    fUpwindQuad[14] = ser_5x_p2_surfx3_eval_quad_node_14_r(fEdge); 
    fUpwindQuad[15] = ser_5x_p2_surfx3_eval_quad_node_15_r(fEdge); 
    fUpwindQuad[16] = ser_5x_p2_surfx3_eval_quad_node_16_r(fEdge); 
    fUpwindQuad[17] = ser_5x_p2_surfx3_eval_quad_node_17_r(fEdge); 
  } else { 
    fUpwindQuad[9] = ser_5x_p2_surfx3_eval_quad_node_9_l(fSkin); 
    fUpwindQuad[10] = ser_5x_p2_surfx3_eval_quad_node_10_l(fSkin); 
    fUpwindQuad[11] = ser_5x_p2_surfx3_eval_quad_node_11_l(fSkin); 
    fUpwindQuad[12] = ser_5x_p2_surfx3_eval_quad_node_12_l(fSkin); 
    fUpwindQuad[13] = ser_5x_p2_surfx3_eval_quad_node_13_l(fSkin); 
    fUpwindQuad[14] = ser_5x_p2_surfx3_eval_quad_node_14_l(fSkin); 
    fUpwindQuad[15] = ser_5x_p2_surfx3_eval_quad_node_15_l(fSkin); 
    fUpwindQuad[16] = ser_5x_p2_surfx3_eval_quad_node_16_l(fSkin); 
    fUpwindQuad[17] = ser_5x_p2_surfx3_eval_quad_node_17_l(fSkin); 
  } 
  if (0.223606797749979*alpha[12]+0.45*alpha[5]-0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[18] = ser_5x_p2_surfx3_eval_quad_node_18_r(fEdge); 
    fUpwindQuad[19] = ser_5x_p2_surfx3_eval_quad_node_19_r(fEdge); 
    fUpwindQuad[20] = ser_5x_p2_surfx3_eval_quad_node_20_r(fEdge); 
    fUpwindQuad[21] = ser_5x_p2_surfx3_eval_quad_node_21_r(fEdge); 
    fUpwindQuad[22] = ser_5x_p2_surfx3_eval_quad_node_22_r(fEdge); 
    fUpwindQuad[23] = ser_5x_p2_surfx3_eval_quad_node_23_r(fEdge); 
    fUpwindQuad[24] = ser_5x_p2_surfx3_eval_quad_node_24_r(fEdge); 
    fUpwindQuad[25] = ser_5x_p2_surfx3_eval_quad_node_25_r(fEdge); 
    fUpwindQuad[26] = ser_5x_p2_surfx3_eval_quad_node_26_r(fEdge); 
  } else { 
    fUpwindQuad[18] = ser_5x_p2_surfx3_eval_quad_node_18_l(fSkin); 
    fUpwindQuad[19] = ser_5x_p2_surfx3_eval_quad_node_19_l(fSkin); 
    fUpwindQuad[20] = ser_5x_p2_surfx3_eval_quad_node_20_l(fSkin); 
    fUpwindQuad[21] = ser_5x_p2_surfx3_eval_quad_node_21_l(fSkin); 
    fUpwindQuad[22] = ser_5x_p2_surfx3_eval_quad_node_22_l(fSkin); 
    fUpwindQuad[23] = ser_5x_p2_surfx3_eval_quad_node_23_l(fSkin); 
    fUpwindQuad[24] = ser_5x_p2_surfx3_eval_quad_node_24_l(fSkin); 
    fUpwindQuad[25] = ser_5x_p2_surfx3_eval_quad_node_25_l(fSkin); 
    fUpwindQuad[26] = ser_5x_p2_surfx3_eval_quad_node_26_l(fSkin); 
  } 
  if (0.223606797749979*alpha[12]+0.45*alpha[5]-0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[27] = ser_5x_p2_surfx3_eval_quad_node_27_r(fEdge); 
    fUpwindQuad[28] = ser_5x_p2_surfx3_eval_quad_node_28_r(fEdge); 
    fUpwindQuad[29] = ser_5x_p2_surfx3_eval_quad_node_29_r(fEdge); 
    fUpwindQuad[30] = ser_5x_p2_surfx3_eval_quad_node_30_r(fEdge); 
    fUpwindQuad[31] = ser_5x_p2_surfx3_eval_quad_node_31_r(fEdge); 
    fUpwindQuad[32] = ser_5x_p2_surfx3_eval_quad_node_32_r(fEdge); 
    fUpwindQuad[33] = ser_5x_p2_surfx3_eval_quad_node_33_r(fEdge); 
    fUpwindQuad[34] = ser_5x_p2_surfx3_eval_quad_node_34_r(fEdge); 
    fUpwindQuad[35] = ser_5x_p2_surfx3_eval_quad_node_35_r(fEdge); 
  } else { 
    fUpwindQuad[27] = ser_5x_p2_surfx3_eval_quad_node_27_l(fSkin); 
    fUpwindQuad[28] = ser_5x_p2_surfx3_eval_quad_node_28_l(fSkin); 
    fUpwindQuad[29] = ser_5x_p2_surfx3_eval_quad_node_29_l(fSkin); 
    fUpwindQuad[30] = ser_5x_p2_surfx3_eval_quad_node_30_l(fSkin); 
    fUpwindQuad[31] = ser_5x_p2_surfx3_eval_quad_node_31_l(fSkin); 
    fUpwindQuad[32] = ser_5x_p2_surfx3_eval_quad_node_32_l(fSkin); 
    fUpwindQuad[33] = ser_5x_p2_surfx3_eval_quad_node_33_l(fSkin); 
    fUpwindQuad[34] = ser_5x_p2_surfx3_eval_quad_node_34_l(fSkin); 
    fUpwindQuad[35] = ser_5x_p2_surfx3_eval_quad_node_35_l(fSkin); 
  } 
  if (0.223606797749979*alpha[12]+0.45*alpha[5]-0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[36] = ser_5x_p2_surfx3_eval_quad_node_36_r(fEdge); 
    fUpwindQuad[37] = ser_5x_p2_surfx3_eval_quad_node_37_r(fEdge); 
    fUpwindQuad[38] = ser_5x_p2_surfx3_eval_quad_node_38_r(fEdge); 
    fUpwindQuad[39] = ser_5x_p2_surfx3_eval_quad_node_39_r(fEdge); 
    fUpwindQuad[40] = ser_5x_p2_surfx3_eval_quad_node_40_r(fEdge); 
    fUpwindQuad[41] = ser_5x_p2_surfx3_eval_quad_node_41_r(fEdge); 
    fUpwindQuad[42] = ser_5x_p2_surfx3_eval_quad_node_42_r(fEdge); 
    fUpwindQuad[43] = ser_5x_p2_surfx3_eval_quad_node_43_r(fEdge); 
    fUpwindQuad[44] = ser_5x_p2_surfx3_eval_quad_node_44_r(fEdge); 
  } else { 
    fUpwindQuad[36] = ser_5x_p2_surfx3_eval_quad_node_36_l(fSkin); 
    fUpwindQuad[37] = ser_5x_p2_surfx3_eval_quad_node_37_l(fSkin); 
    fUpwindQuad[38] = ser_5x_p2_surfx3_eval_quad_node_38_l(fSkin); 
    fUpwindQuad[39] = ser_5x_p2_surfx3_eval_quad_node_39_l(fSkin); 
    fUpwindQuad[40] = ser_5x_p2_surfx3_eval_quad_node_40_l(fSkin); 
    fUpwindQuad[41] = ser_5x_p2_surfx3_eval_quad_node_41_l(fSkin); 
    fUpwindQuad[42] = ser_5x_p2_surfx3_eval_quad_node_42_l(fSkin); 
    fUpwindQuad[43] = ser_5x_p2_surfx3_eval_quad_node_43_l(fSkin); 
    fUpwindQuad[44] = ser_5x_p2_surfx3_eval_quad_node_44_l(fSkin); 
  } 
  if (0.223606797749979*alpha[12]+0.45*alpha[5]-0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[45] = ser_5x_p2_surfx3_eval_quad_node_45_r(fEdge); 
    fUpwindQuad[46] = ser_5x_p2_surfx3_eval_quad_node_46_r(fEdge); 
    fUpwindQuad[47] = ser_5x_p2_surfx3_eval_quad_node_47_r(fEdge); 
    fUpwindQuad[48] = ser_5x_p2_surfx3_eval_quad_node_48_r(fEdge); 
    fUpwindQuad[49] = ser_5x_p2_surfx3_eval_quad_node_49_r(fEdge); 
    fUpwindQuad[50] = ser_5x_p2_surfx3_eval_quad_node_50_r(fEdge); 
    fUpwindQuad[51] = ser_5x_p2_surfx3_eval_quad_node_51_r(fEdge); 
    fUpwindQuad[52] = ser_5x_p2_surfx3_eval_quad_node_52_r(fEdge); 
    fUpwindQuad[53] = ser_5x_p2_surfx3_eval_quad_node_53_r(fEdge); 
  } else { 
    fUpwindQuad[45] = ser_5x_p2_surfx3_eval_quad_node_45_l(fSkin); 
    fUpwindQuad[46] = ser_5x_p2_surfx3_eval_quad_node_46_l(fSkin); 
    fUpwindQuad[47] = ser_5x_p2_surfx3_eval_quad_node_47_l(fSkin); 
    fUpwindQuad[48] = ser_5x_p2_surfx3_eval_quad_node_48_l(fSkin); 
    fUpwindQuad[49] = ser_5x_p2_surfx3_eval_quad_node_49_l(fSkin); 
    fUpwindQuad[50] = ser_5x_p2_surfx3_eval_quad_node_50_l(fSkin); 
    fUpwindQuad[51] = ser_5x_p2_surfx3_eval_quad_node_51_l(fSkin); 
    fUpwindQuad[52] = ser_5x_p2_surfx3_eval_quad_node_52_l(fSkin); 
    fUpwindQuad[53] = ser_5x_p2_surfx3_eval_quad_node_53_l(fSkin); 
  } 
  if (0.223606797749979*alpha[12]+0.45*alpha[5]-0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[54] = ser_5x_p2_surfx3_eval_quad_node_54_r(fEdge); 
    fUpwindQuad[55] = ser_5x_p2_surfx3_eval_quad_node_55_r(fEdge); 
    fUpwindQuad[56] = ser_5x_p2_surfx3_eval_quad_node_56_r(fEdge); 
    fUpwindQuad[57] = ser_5x_p2_surfx3_eval_quad_node_57_r(fEdge); 
    fUpwindQuad[58] = ser_5x_p2_surfx3_eval_quad_node_58_r(fEdge); 
    fUpwindQuad[59] = ser_5x_p2_surfx3_eval_quad_node_59_r(fEdge); 
    fUpwindQuad[60] = ser_5x_p2_surfx3_eval_quad_node_60_r(fEdge); 
    fUpwindQuad[61] = ser_5x_p2_surfx3_eval_quad_node_61_r(fEdge); 
    fUpwindQuad[62] = ser_5x_p2_surfx3_eval_quad_node_62_r(fEdge); 
  } else { 
    fUpwindQuad[54] = ser_5x_p2_surfx3_eval_quad_node_54_l(fSkin); 
    fUpwindQuad[55] = ser_5x_p2_surfx3_eval_quad_node_55_l(fSkin); 
    fUpwindQuad[56] = ser_5x_p2_surfx3_eval_quad_node_56_l(fSkin); 
    fUpwindQuad[57] = ser_5x_p2_surfx3_eval_quad_node_57_l(fSkin); 
    fUpwindQuad[58] = ser_5x_p2_surfx3_eval_quad_node_58_l(fSkin); 
    fUpwindQuad[59] = ser_5x_p2_surfx3_eval_quad_node_59_l(fSkin); 
    fUpwindQuad[60] = ser_5x_p2_surfx3_eval_quad_node_60_l(fSkin); 
    fUpwindQuad[61] = ser_5x_p2_surfx3_eval_quad_node_61_l(fSkin); 
    fUpwindQuad[62] = ser_5x_p2_surfx3_eval_quad_node_62_l(fSkin); 
  } 
  if (0.223606797749979*alpha[12]+0.45*alpha[5]-0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[63] = ser_5x_p2_surfx3_eval_quad_node_63_r(fEdge); 
    fUpwindQuad[64] = ser_5x_p2_surfx3_eval_quad_node_64_r(fEdge); 
    fUpwindQuad[65] = ser_5x_p2_surfx3_eval_quad_node_65_r(fEdge); 
    fUpwindQuad[66] = ser_5x_p2_surfx3_eval_quad_node_66_r(fEdge); 
    fUpwindQuad[67] = ser_5x_p2_surfx3_eval_quad_node_67_r(fEdge); 
    fUpwindQuad[68] = ser_5x_p2_surfx3_eval_quad_node_68_r(fEdge); 
    fUpwindQuad[69] = ser_5x_p2_surfx3_eval_quad_node_69_r(fEdge); 
    fUpwindQuad[70] = ser_5x_p2_surfx3_eval_quad_node_70_r(fEdge); 
    fUpwindQuad[71] = ser_5x_p2_surfx3_eval_quad_node_71_r(fEdge); 
  } else { 
    fUpwindQuad[63] = ser_5x_p2_surfx3_eval_quad_node_63_l(fSkin); 
    fUpwindQuad[64] = ser_5x_p2_surfx3_eval_quad_node_64_l(fSkin); 
    fUpwindQuad[65] = ser_5x_p2_surfx3_eval_quad_node_65_l(fSkin); 
    fUpwindQuad[66] = ser_5x_p2_surfx3_eval_quad_node_66_l(fSkin); 
    fUpwindQuad[67] = ser_5x_p2_surfx3_eval_quad_node_67_l(fSkin); 
    fUpwindQuad[68] = ser_5x_p2_surfx3_eval_quad_node_68_l(fSkin); 
    fUpwindQuad[69] = ser_5x_p2_surfx3_eval_quad_node_69_l(fSkin); 
    fUpwindQuad[70] = ser_5x_p2_surfx3_eval_quad_node_70_l(fSkin); 
    fUpwindQuad[71] = ser_5x_p2_surfx3_eval_quad_node_71_l(fSkin); 
  } 
  if (0.223606797749979*alpha[12]+0.45*alpha[5]-0.3354101966249685*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[72] = ser_5x_p2_surfx3_eval_quad_node_72_r(fEdge); 
    fUpwindQuad[73] = ser_5x_p2_surfx3_eval_quad_node_73_r(fEdge); 
    fUpwindQuad[74] = ser_5x_p2_surfx3_eval_quad_node_74_r(fEdge); 
    fUpwindQuad[75] = ser_5x_p2_surfx3_eval_quad_node_75_r(fEdge); 
    fUpwindQuad[76] = ser_5x_p2_surfx3_eval_quad_node_76_r(fEdge); 
    fUpwindQuad[77] = ser_5x_p2_surfx3_eval_quad_node_77_r(fEdge); 
    fUpwindQuad[78] = ser_5x_p2_surfx3_eval_quad_node_78_r(fEdge); 
    fUpwindQuad[79] = ser_5x_p2_surfx3_eval_quad_node_79_r(fEdge); 
    fUpwindQuad[80] = ser_5x_p2_surfx3_eval_quad_node_80_r(fEdge); 
  } else { 
    fUpwindQuad[72] = ser_5x_p2_surfx3_eval_quad_node_72_l(fSkin); 
    fUpwindQuad[73] = ser_5x_p2_surfx3_eval_quad_node_73_l(fSkin); 
    fUpwindQuad[74] = ser_5x_p2_surfx3_eval_quad_node_74_l(fSkin); 
    fUpwindQuad[75] = ser_5x_p2_surfx3_eval_quad_node_75_l(fSkin); 
    fUpwindQuad[76] = ser_5x_p2_surfx3_eval_quad_node_76_l(fSkin); 
    fUpwindQuad[77] = ser_5x_p2_surfx3_eval_quad_node_77_l(fSkin); 
    fUpwindQuad[78] = ser_5x_p2_surfx3_eval_quad_node_78_l(fSkin); 
    fUpwindQuad[79] = ser_5x_p2_surfx3_eval_quad_node_79_l(fSkin); 
    fUpwindQuad[80] = ser_5x_p2_surfx3_eval_quad_node_80_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_5x_p2_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.25*alpha[12]*fUpwind[12]+0.25*alpha[5]*fUpwind[5]+0.25*alpha[2]*fUpwind[2]+0.25*alpha[1]*fUpwind[1]+0.25*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.2500000000000001*alpha[12]*fUpwind[20]+0.223606797749979*alpha[5]*fUpwind[19]+0.223606797749979*alpha[1]*fUpwind[11]+0.25*alpha[2]*fUpwind[5]+0.25*fUpwind[2]*alpha[5]+0.25*alpha[0]*fUpwind[1]+0.25*fUpwind[0]*alpha[1]; 
  Ghat[2] = 0.223606797749979*alpha[5]*fUpwind[20]+0.223606797749979*alpha[2]*fUpwind[12]+0.223606797749979*fUpwind[2]*alpha[12]+0.25*alpha[1]*fUpwind[5]+0.25*fUpwind[1]*alpha[5]+0.25*alpha[0]*fUpwind[2]+0.25*fUpwind[0]*alpha[2]; 
  Ghat[3] = 0.2500000000000001*alpha[12]*fUpwind[22]+0.25*alpha[5]*fUpwind[15]+0.25*alpha[2]*fUpwind[7]+0.25*alpha[1]*fUpwind[6]+0.25*alpha[0]*fUpwind[3]; 
  Ghat[4] = 0.2500000000000001*alpha[12]*fUpwind[26]+0.25*alpha[5]*fUpwind[16]+0.25*alpha[2]*fUpwind[9]+0.25*alpha[1]*fUpwind[8]+0.25*alpha[0]*fUpwind[4]; 
  Ghat[5] = 0.223606797749979*alpha[2]*fUpwind[20]+0.223606797749979*alpha[1]*fUpwind[19]+0.223606797749979*alpha[5]*fUpwind[12]+0.223606797749979*fUpwind[5]*alpha[12]+0.223606797749979*alpha[5]*fUpwind[11]+0.25*alpha[0]*fUpwind[5]+0.25*fUpwind[0]*alpha[5]+0.25*alpha[1]*fUpwind[2]+0.25*fUpwind[1]*alpha[2]; 
  Ghat[6] = 0.25*alpha[12]*fUpwind[33]+0.223606797749979*alpha[5]*fUpwind[32]+0.223606797749979*alpha[1]*fUpwind[21]+0.25*alpha[2]*fUpwind[15]+0.25*alpha[5]*fUpwind[7]+0.25*alpha[0]*fUpwind[6]+0.25*alpha[1]*fUpwind[3]; 
  Ghat[7] = 0.223606797749979*alpha[5]*fUpwind[33]+0.223606797749979*alpha[2]*fUpwind[22]+0.25*alpha[1]*fUpwind[15]+0.223606797749979*fUpwind[7]*alpha[12]+0.25*alpha[0]*fUpwind[7]+0.25*alpha[5]*fUpwind[6]+0.25*alpha[2]*fUpwind[3]; 
  Ghat[8] = 0.25*alpha[12]*fUpwind[36]+0.223606797749979*alpha[5]*fUpwind[35]+0.223606797749979*alpha[1]*fUpwind[25]+0.25*alpha[2]*fUpwind[16]+0.25*alpha[5]*fUpwind[9]+0.25*alpha[0]*fUpwind[8]+0.25*alpha[1]*fUpwind[4]; 
  Ghat[9] = 0.223606797749979*alpha[5]*fUpwind[36]+0.223606797749979*alpha[2]*fUpwind[26]+0.25*alpha[1]*fUpwind[16]+0.223606797749979*fUpwind[9]*alpha[12]+0.25*alpha[0]*fUpwind[9]+0.25*alpha[5]*fUpwind[8]+0.25*alpha[2]*fUpwind[4]; 
  Ghat[10] = 0.25*alpha[12]*fUpwind[38]+0.25*alpha[5]*fUpwind[31]+0.25*alpha[2]*fUpwind[18]+0.25*alpha[1]*fUpwind[17]+0.25*alpha[0]*fUpwind[10]; 
  Ghat[11] = 0.2500000000000001*alpha[2]*fUpwind[19]+0.25*alpha[0]*fUpwind[11]+0.223606797749979*alpha[5]*fUpwind[5]+0.223606797749979*alpha[1]*fUpwind[1]; 
  Ghat[12] = 0.2500000000000001*alpha[1]*fUpwind[20]+0.159719141249985*alpha[12]*fUpwind[12]+0.25*alpha[0]*fUpwind[12]+0.25*fUpwind[0]*alpha[12]+0.223606797749979*alpha[5]*fUpwind[5]+0.223606797749979*alpha[2]*fUpwind[2]; 
  Ghat[13] = 0.25*alpha[5]*fUpwind[34]+0.2500000000000001*alpha[2]*fUpwind[24]+0.2500000000000001*alpha[1]*fUpwind[23]+0.25*alpha[0]*fUpwind[13]; 
  Ghat[14] = 0.25*alpha[5]*fUpwind[41]+0.2500000000000001*alpha[2]*fUpwind[29]+0.2500000000000001*alpha[1]*fUpwind[28]+0.25*alpha[0]*fUpwind[14]; 
  Ghat[15] = 0.223606797749979*alpha[2]*fUpwind[33]+0.223606797749979*alpha[1]*fUpwind[32]+0.223606797749979*alpha[5]*fUpwind[22]+0.223606797749979*alpha[5]*fUpwind[21]+0.223606797749979*alpha[12]*fUpwind[15]+0.25*alpha[0]*fUpwind[15]+0.25*alpha[1]*fUpwind[7]+0.25*alpha[2]*fUpwind[6]+0.25*fUpwind[3]*alpha[5]; 
  Ghat[16] = 0.223606797749979*alpha[2]*fUpwind[36]+0.223606797749979*alpha[1]*fUpwind[35]+0.223606797749979*alpha[5]*fUpwind[26]+0.223606797749979*alpha[5]*fUpwind[25]+0.223606797749979*alpha[12]*fUpwind[16]+0.25*alpha[0]*fUpwind[16]+0.25*alpha[1]*fUpwind[9]+0.25*alpha[2]*fUpwind[8]+0.25*fUpwind[4]*alpha[5]; 
  Ghat[17] = 0.2500000000000001*alpha[12]*fUpwind[45]+0.223606797749979*alpha[5]*fUpwind[44]+0.223606797749979*alpha[1]*fUpwind[37]+0.25*alpha[2]*fUpwind[31]+0.25*alpha[5]*fUpwind[18]+0.25*alpha[0]*fUpwind[17]+0.25*alpha[1]*fUpwind[10]; 
  Ghat[18] = 0.223606797749979*alpha[5]*fUpwind[45]+0.223606797749979*alpha[2]*fUpwind[38]+0.25*alpha[1]*fUpwind[31]+0.223606797749979*alpha[12]*fUpwind[18]+0.25*alpha[0]*fUpwind[18]+0.25*alpha[5]*fUpwind[17]+0.25*alpha[2]*fUpwind[10]; 
  Ghat[19] = 0.2*alpha[5]*fUpwind[20]+0.223606797749979*alpha[12]*fUpwind[19]+0.25*alpha[0]*fUpwind[19]+0.2500000000000001*alpha[2]*fUpwind[11]+0.223606797749979*alpha[1]*fUpwind[5]+0.223606797749979*fUpwind[1]*alpha[5]; 
  Ghat[20] = 0.159719141249985*alpha[12]*fUpwind[20]+0.25*alpha[0]*fUpwind[20]+0.2*alpha[5]*fUpwind[19]+0.2500000000000001*alpha[1]*fUpwind[12]+0.2500000000000001*fUpwind[1]*alpha[12]+0.223606797749979*alpha[2]*fUpwind[5]+0.223606797749979*fUpwind[2]*alpha[5]; 
  Ghat[21] = 0.2500000000000001*alpha[2]*fUpwind[32]+0.25*alpha[0]*fUpwind[21]+0.223606797749979*alpha[5]*fUpwind[15]+0.223606797749979*alpha[1]*fUpwind[6]; 
  Ghat[22] = 0.2500000000000001*alpha[1]*fUpwind[33]+0.159719141249985*alpha[12]*fUpwind[22]+0.25*alpha[0]*fUpwind[22]+0.223606797749979*alpha[5]*fUpwind[15]+0.2500000000000001*fUpwind[3]*alpha[12]+0.223606797749979*alpha[2]*fUpwind[7]; 
  Ghat[23] = 0.2500000000000001*alpha[2]*fUpwind[34]+0.25*alpha[5]*fUpwind[24]+0.25*alpha[0]*fUpwind[23]+0.2500000000000001*alpha[1]*fUpwind[13]; 
  Ghat[24] = 0.2500000000000001*alpha[1]*fUpwind[34]+0.223606797749979*alpha[12]*fUpwind[24]+0.25*alpha[0]*fUpwind[24]+0.25*alpha[5]*fUpwind[23]+0.2500000000000001*alpha[2]*fUpwind[13]; 
  Ghat[25] = 0.2500000000000001*alpha[2]*fUpwind[35]+0.25*alpha[0]*fUpwind[25]+0.223606797749979*alpha[5]*fUpwind[16]+0.223606797749979*alpha[1]*fUpwind[8]; 
  Ghat[26] = 0.2500000000000001*alpha[1]*fUpwind[36]+0.159719141249985*alpha[12]*fUpwind[26]+0.25*alpha[0]*fUpwind[26]+0.223606797749979*alpha[5]*fUpwind[16]+0.2500000000000001*fUpwind[4]*alpha[12]+0.223606797749979*alpha[2]*fUpwind[9]; 
  Ghat[27] = 0.25*alpha[5]*fUpwind[46]+0.2500000000000001*alpha[2]*fUpwind[40]+0.2500000000000001*alpha[1]*fUpwind[39]+0.25*alpha[0]*fUpwind[27]; 
  Ghat[28] = 0.2500000000000001*alpha[2]*fUpwind[41]+0.25*alpha[5]*fUpwind[29]+0.25*alpha[0]*fUpwind[28]+0.2500000000000001*alpha[1]*fUpwind[14]; 
  Ghat[29] = 0.2500000000000001*alpha[1]*fUpwind[41]+0.223606797749979*alpha[12]*fUpwind[29]+0.25*alpha[0]*fUpwind[29]+0.25*alpha[5]*fUpwind[28]+0.2500000000000001*alpha[2]*fUpwind[14]; 
  Ghat[30] = 0.25*alpha[5]*fUpwind[47]+0.2500000000000001*alpha[2]*fUpwind[43]+0.2500000000000001*alpha[1]*fUpwind[42]+0.25*alpha[0]*fUpwind[30]; 
  Ghat[31] = 0.223606797749979*alpha[2]*fUpwind[45]+0.223606797749979*alpha[1]*fUpwind[44]+0.223606797749979*alpha[5]*fUpwind[38]+0.223606797749979*alpha[5]*fUpwind[37]+0.223606797749979*alpha[12]*fUpwind[31]+0.25*alpha[0]*fUpwind[31]+0.25*alpha[1]*fUpwind[18]+0.25*alpha[2]*fUpwind[17]+0.25*alpha[5]*fUpwind[10]; 
  Ghat[32] = 0.2*alpha[5]*fUpwind[33]+0.223606797749979*alpha[12]*fUpwind[32]+0.25*alpha[0]*fUpwind[32]+0.2500000000000001*alpha[2]*fUpwind[21]+0.223606797749979*alpha[1]*fUpwind[15]+0.223606797749979*alpha[5]*fUpwind[6]; 
  Ghat[33] = 0.159719141249985*alpha[12]*fUpwind[33]+0.25*alpha[0]*fUpwind[33]+0.2*alpha[5]*fUpwind[32]+0.2500000000000001*alpha[1]*fUpwind[22]+0.223606797749979*alpha[2]*fUpwind[15]+0.25*fUpwind[6]*alpha[12]+0.223606797749979*alpha[5]*fUpwind[7]; 
  Ghat[34] = 0.223606797749979*alpha[12]*fUpwind[34]+0.25*alpha[0]*fUpwind[34]+0.2500000000000001*alpha[1]*fUpwind[24]+0.2500000000000001*alpha[2]*fUpwind[23]+0.25*alpha[5]*fUpwind[13]; 
  Ghat[35] = 0.2*alpha[5]*fUpwind[36]+0.223606797749979*alpha[12]*fUpwind[35]+0.25*alpha[0]*fUpwind[35]+0.2500000000000001*alpha[2]*fUpwind[25]+0.223606797749979*alpha[1]*fUpwind[16]+0.223606797749979*alpha[5]*fUpwind[8]; 
  Ghat[36] = 0.159719141249985*alpha[12]*fUpwind[36]+0.25*alpha[0]*fUpwind[36]+0.2*alpha[5]*fUpwind[35]+0.2500000000000001*alpha[1]*fUpwind[26]+0.223606797749979*alpha[2]*fUpwind[16]+0.25*fUpwind[8]*alpha[12]+0.223606797749979*alpha[5]*fUpwind[9]; 
  Ghat[37] = 0.2500000000000001*alpha[2]*fUpwind[44]+0.25*alpha[0]*fUpwind[37]+0.223606797749979*alpha[5]*fUpwind[31]+0.223606797749979*alpha[1]*fUpwind[17]; 
  Ghat[38] = 0.2500000000000001*alpha[1]*fUpwind[45]+0.159719141249985*alpha[12]*fUpwind[38]+0.25*alpha[0]*fUpwind[38]+0.223606797749979*alpha[5]*fUpwind[31]+0.223606797749979*alpha[2]*fUpwind[18]+0.25*fUpwind[10]*alpha[12]; 
  Ghat[39] = 0.2500000000000001*alpha[2]*fUpwind[46]+0.25*alpha[5]*fUpwind[40]+0.25*alpha[0]*fUpwind[39]+0.2500000000000001*alpha[1]*fUpwind[27]; 
  Ghat[40] = 0.2500000000000001*alpha[1]*fUpwind[46]+0.223606797749979*alpha[12]*fUpwind[40]+0.25*alpha[0]*fUpwind[40]+0.25*alpha[5]*fUpwind[39]+0.2500000000000001*alpha[2]*fUpwind[27]; 
  Ghat[41] = 0.223606797749979*alpha[12]*fUpwind[41]+0.25*alpha[0]*fUpwind[41]+0.2500000000000001*alpha[1]*fUpwind[29]+0.2500000000000001*alpha[2]*fUpwind[28]+0.25*alpha[5]*fUpwind[14]; 
  Ghat[42] = 0.2500000000000001*alpha[2]*fUpwind[47]+0.25*alpha[5]*fUpwind[43]+0.25*alpha[0]*fUpwind[42]+0.2500000000000001*alpha[1]*fUpwind[30]; 
  Ghat[43] = 0.2500000000000001*alpha[1]*fUpwind[47]+0.223606797749979*alpha[12]*fUpwind[43]+0.25*alpha[0]*fUpwind[43]+0.25*alpha[5]*fUpwind[42]+0.2500000000000001*alpha[2]*fUpwind[30]; 
  Ghat[44] = 0.2*alpha[5]*fUpwind[45]+0.223606797749979*alpha[12]*fUpwind[44]+0.25*alpha[0]*fUpwind[44]+0.2500000000000001*alpha[2]*fUpwind[37]+0.223606797749979*alpha[1]*fUpwind[31]+0.223606797749979*alpha[5]*fUpwind[17]; 
  Ghat[45] = 0.159719141249985*alpha[12]*fUpwind[45]+0.25*alpha[0]*fUpwind[45]+0.2*alpha[5]*fUpwind[44]+0.2500000000000001*alpha[1]*fUpwind[38]+0.223606797749979*alpha[2]*fUpwind[31]+0.223606797749979*alpha[5]*fUpwind[18]+0.2500000000000001*alpha[12]*fUpwind[17]; 
  Ghat[46] = 0.223606797749979*alpha[12]*fUpwind[46]+0.25*alpha[0]*fUpwind[46]+0.2500000000000001*alpha[1]*fUpwind[40]+0.2500000000000001*alpha[2]*fUpwind[39]+0.25*alpha[5]*fUpwind[27]; 
  Ghat[47] = 0.223606797749979*alpha[12]*fUpwind[47]+0.25*alpha[0]*fUpwind[47]+0.2500000000000001*alpha[1]*fUpwind[43]+0.2500000000000001*alpha[2]*fUpwind[42]+0.25*alpha[5]*fUpwind[30]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv10; 
  out[1] += 0.7071067811865475*Ghat[1]*dv10; 
  out[2] += 0.7071067811865475*Ghat[2]*dv10; 
  out[3] += -1.224744871391589*Ghat[0]*dv10; 
  out[4] += 0.7071067811865475*Ghat[3]*dv10; 
  out[5] += 0.7071067811865475*Ghat[4]*dv10; 
  out[6] += 0.7071067811865475*Ghat[5]*dv10; 
  out[7] += -1.224744871391589*Ghat[1]*dv10; 
  out[8] += -1.224744871391589*Ghat[2]*dv10; 
  out[9] += 0.7071067811865475*Ghat[6]*dv10; 
  out[10] += 0.7071067811865475*Ghat[7]*dv10; 
  out[11] += -1.224744871391589*Ghat[3]*dv10; 
  out[12] += 0.7071067811865475*Ghat[8]*dv10; 
  out[13] += 0.7071067811865475*Ghat[9]*dv10; 
  out[14] += -1.224744871391589*Ghat[4]*dv10; 
  out[15] += 0.7071067811865475*Ghat[10]*dv10; 
  out[16] += 0.7071067811865475*Ghat[11]*dv10; 
  out[17] += 0.7071067811865475*Ghat[12]*dv10; 
  out[18] += 1.58113883008419*Ghat[0]*dv10; 
  out[19] += 0.7071067811865475*Ghat[13]*dv10; 
  out[20] += 0.7071067811865475*Ghat[14]*dv10; 
  out[21] += -1.224744871391589*Ghat[5]*dv10; 
  out[22] += 0.7071067811865475*Ghat[15]*dv10; 
  out[23] += -1.224744871391589*Ghat[6]*dv10; 
  out[24] += -1.224744871391589*Ghat[7]*dv10; 
  out[25] += 0.7071067811865475*Ghat[16]*dv10; 
  out[26] += -1.224744871391589*Ghat[8]*dv10; 
  out[27] += -1.224744871391589*Ghat[9]*dv10; 
  out[28] += 0.7071067811865475*Ghat[17]*dv10; 
  out[29] += 0.7071067811865475*Ghat[18]*dv10; 
  out[30] += -1.224744871391589*Ghat[10]*dv10; 
  out[31] += 0.7071067811865475*Ghat[19]*dv10; 
  out[32] += 0.7071067811865475*Ghat[20]*dv10; 
  out[33] += -1.224744871391589*Ghat[11]*dv10; 
  out[34] += -1.224744871391589*Ghat[12]*dv10; 
  out[35] += 1.58113883008419*Ghat[1]*dv10; 
  out[36] += 1.58113883008419*Ghat[2]*dv10; 
  out[37] += 0.7071067811865475*Ghat[21]*dv10; 
  out[38] += 0.7071067811865475*Ghat[22]*dv10; 
  out[39] += 1.58113883008419*Ghat[3]*dv10; 
  out[40] += 0.7071067811865475*Ghat[23]*dv10; 
  out[41] += 0.7071067811865475*Ghat[24]*dv10; 
  out[42] += -1.224744871391589*Ghat[13]*dv10; 
  out[43] += 0.7071067811865475*Ghat[25]*dv10; 
  out[44] += 0.7071067811865475*Ghat[26]*dv10; 
  out[45] += 1.58113883008419*Ghat[4]*dv10; 
  out[46] += 0.7071067811865475*Ghat[27]*dv10; 
  out[47] += 0.7071067811865475*Ghat[28]*dv10; 
  out[48] += 0.7071067811865475*Ghat[29]*dv10; 
  out[49] += -1.224744871391589*Ghat[14]*dv10; 
  out[50] += 0.7071067811865475*Ghat[30]*dv10; 
  out[51] += -1.224744871391589*Ghat[15]*dv10; 
  out[52] += -1.224744871391589*Ghat[16]*dv10; 
  out[53] += 0.7071067811865475*Ghat[31]*dv10; 
  out[54] += -1.224744871391589*Ghat[17]*dv10; 
  out[55] += -1.224744871391589*Ghat[18]*dv10; 
  out[56] += -1.224744871391589*Ghat[19]*dv10; 
  out[57] += -1.224744871391589*Ghat[20]*dv10; 
  out[58] += 1.58113883008419*Ghat[5]*dv10; 
  out[59] += 0.7071067811865475*Ghat[32]*dv10; 
  out[60] += 0.7071067811865475*Ghat[33]*dv10; 
  out[61] += -1.224744871391589*Ghat[21]*dv10; 
  out[62] += -1.224744871391589*Ghat[22]*dv10; 
  out[63] += 1.58113883008419*Ghat[6]*dv10; 
  out[64] += 1.58113883008419*Ghat[7]*dv10; 
  out[65] += 0.7071067811865475*Ghat[34]*dv10; 
  out[66] += -1.224744871391589*Ghat[23]*dv10; 
  out[67] += -1.224744871391589*Ghat[24]*dv10; 
  out[68] += 0.7071067811865475*Ghat[35]*dv10; 
  out[69] += 0.7071067811865475*Ghat[36]*dv10; 
  out[70] += -1.224744871391589*Ghat[25]*dv10; 
  out[71] += -1.224744871391589*Ghat[26]*dv10; 
  out[72] += 1.58113883008419*Ghat[8]*dv10; 
  out[73] += 1.58113883008419*Ghat[9]*dv10; 
  out[74] += 0.7071067811865475*Ghat[37]*dv10; 
  out[75] += 0.7071067811865475*Ghat[38]*dv10; 
  out[76] += 1.58113883008419*Ghat[10]*dv10; 
  out[77] += 0.7071067811865475*Ghat[39]*dv10; 
  out[78] += 0.7071067811865475*Ghat[40]*dv10; 
  out[79] += -1.224744871391589*Ghat[27]*dv10; 
  out[80] += 0.7071067811865475*Ghat[41]*dv10; 
  out[81] += -1.224744871391589*Ghat[28]*dv10; 
  out[82] += -1.224744871391589*Ghat[29]*dv10; 
  out[83] += 0.7071067811865475*Ghat[42]*dv10; 
  out[84] += 0.7071067811865475*Ghat[43]*dv10; 
  out[85] += -1.224744871391589*Ghat[30]*dv10; 
  out[86] += -1.224744871391589*Ghat[31]*dv10; 
  out[87] += -1.224744871391589*Ghat[32]*dv10; 
  out[88] += -1.224744871391589*Ghat[33]*dv10; 
  out[89] += 1.58113883008419*Ghat[15]*dv10; 
  out[90] += -1.224744871391589*Ghat[34]*dv10; 
  out[91] += -1.224744871391589*Ghat[35]*dv10; 
  out[92] += -1.224744871391589*Ghat[36]*dv10; 
  out[93] += 1.58113883008419*Ghat[16]*dv10; 
  out[94] += 0.7071067811865475*Ghat[44]*dv10; 
  out[95] += 0.7071067811865475*Ghat[45]*dv10; 
  out[96] += -1.224744871391589*Ghat[37]*dv10; 
  out[97] += -1.224744871391589*Ghat[38]*dv10; 
  out[98] += 1.58113883008419*Ghat[17]*dv10; 
  out[99] += 1.58113883008419*Ghat[18]*dv10; 
  out[100] += 0.7071067811865475*Ghat[46]*dv10; 
  out[101] += -1.224744871391589*Ghat[39]*dv10; 
  out[102] += -1.224744871391589*Ghat[40]*dv10; 
  out[103] += -1.224744871391589*Ghat[41]*dv10; 
  out[104] += 0.7071067811865475*Ghat[47]*dv10; 
  out[105] += -1.224744871391589*Ghat[42]*dv10; 
  out[106] += -1.224744871391589*Ghat[43]*dv10; 
  out[107] += -1.224744871391589*Ghat[44]*dv10; 
  out[108] += -1.224744871391589*Ghat[45]*dv10; 
  out[109] += 1.58113883008419*Ghat[31]*dv10; 
  out[110] += -1.224744871391589*Ghat[46]*dv10; 
  out[111] += -1.224744871391589*Ghat[47]*dv10; 

  } 

  return 0.;

} 
