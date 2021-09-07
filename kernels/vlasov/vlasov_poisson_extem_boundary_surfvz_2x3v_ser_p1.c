#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_2x3v_p1_surfvz_quad.h> 
GKYL_CU_DH void vlasov_poisson_extem_boundary_surfvz_2x3v_ser_p1(const double *w, const double *dxv, const double *fac_phi, const double *vecA, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // fac_phi:     potential (scaled by appropriate factors).
  // vecA:        vector potential (scaled by appropriate factors). Unused in pure Vlasov-Poisson. 
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv12 = 2/dxv[4]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double dv3 = dxv[4], wv3 = w[4]; 
  const double *phi = &fac_phi[0]; 
  const double dx10 = 2/dxv[0]; 
  const double dx11 = 2/dxv[1]; 
  const double *A0 = &vecA[0]; 
  const double *A1 = &vecA[4]; 
  const double *A2 = &vecA[8]; 
  double alpha[16] = {0.0}; 

  alpha[0] = (-3.464101615137754*A2[2]*dx11*wv2)-3.464101615137754*A2[1]*dx10*wv1; 
  alpha[1] = -3.464101615137754*A2[3]*dx11*wv2; 
  alpha[2] = -3.464101615137754*A2[3]*dx10*wv1; 
  alpha[3] = -1.0*A2[1]*dv1*dx10; 
  alpha[4] = -1.0*A2[2]*dv2*dx11; 
  alpha[7] = -1.0*A2[3]*dv1*dx10; 
  alpha[8] = -1.0*A2[3]*dv2*dx11; 

  double fUpwindQuad[16] = {0.0};
  double fUpwind[16] = {0.0};
  double Ghat[16] = {0.0}; 

  if (edge == -1) { 

  if (alpha[8]+alpha[7]-alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[0] = ser_2x3v_p1_surfvz_quad_0(1, fSkin); 
  } else { 

    fUpwindQuad[0] = ser_2x3v_p1_surfvz_quad_0(-1, fEdge); 
  } 
  if ((-alpha[8])+alpha[7]-alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[1] = ser_2x3v_p1_surfvz_quad_1(1, fSkin); 
  } else { 

    fUpwindQuad[1] = ser_2x3v_p1_surfvz_quad_1(-1, fEdge); 
  } 
  if (alpha[8]-alpha[7]-alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[2] = ser_2x3v_p1_surfvz_quad_2(1, fSkin); 
  } else { 

    fUpwindQuad[2] = ser_2x3v_p1_surfvz_quad_2(-1, fEdge); 
  } 
  if ((-alpha[8])-alpha[7]-alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[3] = ser_2x3v_p1_surfvz_quad_3(1, fSkin); 
  } else { 

    fUpwindQuad[3] = ser_2x3v_p1_surfvz_quad_3(-1, fEdge); 
  } 
  if (alpha[8]-alpha[7]-alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[4] = ser_2x3v_p1_surfvz_quad_4(1, fSkin); 
  } else { 

    fUpwindQuad[4] = ser_2x3v_p1_surfvz_quad_4(-1, fEdge); 
  } 
  if ((-alpha[8])-alpha[7]-alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[5] = ser_2x3v_p1_surfvz_quad_5(1, fSkin); 
  } else { 

    fUpwindQuad[5] = ser_2x3v_p1_surfvz_quad_5(-1, fEdge); 
  } 
  if (alpha[8]+alpha[7]-alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[6] = ser_2x3v_p1_surfvz_quad_6(1, fSkin); 
  } else { 

    fUpwindQuad[6] = ser_2x3v_p1_surfvz_quad_6(-1, fEdge); 
  } 
  if ((-alpha[8])+alpha[7]-alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[7] = ser_2x3v_p1_surfvz_quad_7(1, fSkin); 
  } else { 

    fUpwindQuad[7] = ser_2x3v_p1_surfvz_quad_7(-1, fEdge); 
  } 
  if ((-alpha[8])+alpha[7]+alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[8] = ser_2x3v_p1_surfvz_quad_8(1, fSkin); 
  } else { 

    fUpwindQuad[8] = ser_2x3v_p1_surfvz_quad_8(-1, fEdge); 
  } 
  if (alpha[8]+alpha[7]+alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[9] = ser_2x3v_p1_surfvz_quad_9(1, fSkin); 
  } else { 

    fUpwindQuad[9] = ser_2x3v_p1_surfvz_quad_9(-1, fEdge); 
  } 
  if ((-alpha[8])-alpha[7]+alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[10] = ser_2x3v_p1_surfvz_quad_10(1, fSkin); 
  } else { 

    fUpwindQuad[10] = ser_2x3v_p1_surfvz_quad_10(-1, fEdge); 
  } 
  if (alpha[8]-alpha[7]+alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[11] = ser_2x3v_p1_surfvz_quad_11(1, fSkin); 
  } else { 

    fUpwindQuad[11] = ser_2x3v_p1_surfvz_quad_11(-1, fEdge); 
  } 
  if ((-alpha[8])-alpha[7]+alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[12] = ser_2x3v_p1_surfvz_quad_12(1, fSkin); 
  } else { 

    fUpwindQuad[12] = ser_2x3v_p1_surfvz_quad_12(-1, fEdge); 
  } 
  if (alpha[8]-alpha[7]+alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[13] = ser_2x3v_p1_surfvz_quad_13(1, fSkin); 
  } else { 

    fUpwindQuad[13] = ser_2x3v_p1_surfvz_quad_13(-1, fEdge); 
  } 
  if ((-alpha[8])+alpha[7]+alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[14] = ser_2x3v_p1_surfvz_quad_14(1, fSkin); 
  } else { 

    fUpwindQuad[14] = ser_2x3v_p1_surfvz_quad_14(-1, fEdge); 
  } 
  if (alpha[8]+alpha[7]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[15] = ser_2x3v_p1_surfvz_quad_15(1, fSkin); 
  } else { 

    fUpwindQuad[15] = ser_2x3v_p1_surfvz_quad_15(-1, fEdge); 
  } 

  fUpwind[0] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.25*(fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*fUpwindQuad[12]+fUpwindQuad[11]-1.0*fUpwindQuad[10]+fUpwindQuad[9]-1.0*fUpwindQuad[8]+fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[2] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12])+fUpwindQuad[11]+fUpwindQuad[10]-1.0*(fUpwindQuad[9]+fUpwindQuad[8])+fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[3] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]-1.0*(fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8])+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[4] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]-1.0*(fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[5] = 0.25*(fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]+fUpwindQuad[11]-1.0*(fUpwindQuad[10]+fUpwindQuad[9])+fUpwindQuad[8]+fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 
  fUpwind[6] = 0.25*(fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*(fUpwindQuad[12]+fUpwindQuad[11])+fUpwindQuad[10]-1.0*fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[7] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10])+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[8] = 0.25*(fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*fUpwindQuad[12]+fUpwindQuad[11]-1.0*fUpwindQuad[10]+fUpwindQuad[9]-1.0*(fUpwindQuad[8]+fUpwindQuad[7])+fUpwindQuad[6]-1.0*fUpwindQuad[5]+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[9] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12])+fUpwindQuad[11]+fUpwindQuad[10]-1.0*(fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6])+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[10] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]-1.0*(fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[11] = 0.25*(fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]-1.0*fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]-1.0*fUpwindQuad[8]+fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[12] = 0.25*(fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]+fUpwindQuad[11]-1.0*(fUpwindQuad[10]+fUpwindQuad[9])+fUpwindQuad[8]-1.0*fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[13] = 0.25*(fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*(fUpwindQuad[12]+fUpwindQuad[11])+fUpwindQuad[10]-1.0*fUpwindQuad[9]+fUpwindQuad[8]-1.0*fUpwindQuad[7]+fUpwindQuad[6]-1.0*fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[14] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10])+fUpwindQuad[9]+fUpwindQuad[8]-1.0*(fUpwindQuad[7]+fUpwindQuad[6])+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[15] = 0.25*(fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]-1.0*fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]-1.0*(fUpwindQuad[8]+fUpwindQuad[7])+fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 

  Ghat[0] += 0.25*(alpha[8]*fUpwind[8]+alpha[7]*fUpwind[7]+alpha[4]*fUpwind[4]+alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] += 0.25*(alpha[7]*fUpwind[11]+alpha[4]*fUpwind[8]+fUpwind[4]*alpha[8]+alpha[3]*fUpwind[6]+alpha[2]*fUpwind[5]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] += 0.25*(alpha[8]*fUpwind[12]+alpha[4]*fUpwind[9]+alpha[3]*fUpwind[7]+fUpwind[3]*alpha[7]+alpha[1]*fUpwind[5]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] += 0.25*(alpha[8]*fUpwind[13]+alpha[4]*fUpwind[10]+alpha[2]*fUpwind[7]+fUpwind[2]*alpha[7]+alpha[1]*fUpwind[6]+alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]); 
  Ghat[4] += 0.25*(alpha[7]*fUpwind[14]+alpha[3]*fUpwind[10]+alpha[2]*fUpwind[9]+alpha[1]*fUpwind[8]+fUpwind[1]*alpha[8]+alpha[0]*fUpwind[4]+fUpwind[0]*alpha[4]); 
  Ghat[5] += 0.25*(alpha[4]*fUpwind[12]+alpha[3]*fUpwind[11]+alpha[8]*fUpwind[9]+fUpwind[6]*alpha[7]+alpha[0]*fUpwind[5]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[6] += 0.25*(alpha[4]*fUpwind[13]+alpha[2]*fUpwind[11]+alpha[8]*fUpwind[10]+fUpwind[5]*alpha[7]+alpha[0]*fUpwind[6]+alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]); 
  Ghat[7] += 0.25*(alpha[8]*fUpwind[15]+alpha[4]*fUpwind[14]+alpha[1]*fUpwind[11]+alpha[0]*fUpwind[7]+fUpwind[0]*alpha[7]+alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]); 
  Ghat[8] += 0.25*(alpha[7]*fUpwind[15]+alpha[3]*fUpwind[13]+alpha[2]*fUpwind[12]+alpha[0]*fUpwind[8]+fUpwind[0]*alpha[8]+alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]); 
  Ghat[9] += 0.25*(alpha[3]*fUpwind[14]+alpha[1]*fUpwind[12]+alpha[7]*fUpwind[10]+alpha[0]*fUpwind[9]+fUpwind[5]*alpha[8]+alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]); 
  Ghat[10] += 0.25*(alpha[2]*fUpwind[14]+alpha[1]*fUpwind[13]+alpha[0]*fUpwind[10]+alpha[7]*fUpwind[9]+fUpwind[6]*alpha[8]+alpha[3]*fUpwind[4]+fUpwind[3]*alpha[4]); 
  Ghat[11] += 0.25*(alpha[4]*fUpwind[15]+alpha[8]*fUpwind[14]+alpha[0]*fUpwind[11]+alpha[1]*fUpwind[7]+fUpwind[1]*alpha[7]+alpha[2]*fUpwind[6]+alpha[3]*fUpwind[5]); 
  Ghat[12] += 0.25*(alpha[3]*fUpwind[15]+alpha[7]*fUpwind[13]+alpha[0]*fUpwind[12]+alpha[1]*fUpwind[9]+alpha[2]*fUpwind[8]+fUpwind[2]*alpha[8]+alpha[4]*fUpwind[5]); 
  Ghat[13] += 0.25*(alpha[2]*fUpwind[15]+alpha[0]*fUpwind[13]+alpha[7]*fUpwind[12]+alpha[1]*fUpwind[10]+alpha[3]*fUpwind[8]+fUpwind[3]*alpha[8]+alpha[4]*fUpwind[6]); 
  Ghat[14] += 0.25*(alpha[1]*fUpwind[15]+alpha[0]*fUpwind[14]+alpha[8]*fUpwind[11]+alpha[2]*fUpwind[10]+alpha[3]*fUpwind[9]+alpha[4]*fUpwind[7]+fUpwind[4]*alpha[7]); 
  Ghat[15] += 0.25*(alpha[0]*fUpwind[15]+alpha[1]*fUpwind[14]+alpha[2]*fUpwind[13]+alpha[3]*fUpwind[12]+alpha[4]*fUpwind[11]+alpha[7]*fUpwind[8]+fUpwind[7]*alpha[8]); 

  out[0] += -0.7071067811865475*Ghat[0]*dv12; 
  out[1] += -0.7071067811865475*Ghat[1]*dv12; 
  out[2] += -0.7071067811865475*Ghat[2]*dv12; 
  out[3] += -0.7071067811865475*Ghat[3]*dv12; 
  out[4] += -0.7071067811865475*Ghat[4]*dv12; 
  out[5] += -1.224744871391589*Ghat[0]*dv12; 
  out[6] += -0.7071067811865475*Ghat[5]*dv12; 
  out[7] += -0.7071067811865475*Ghat[6]*dv12; 
  out[8] += -0.7071067811865475*Ghat[7]*dv12; 
  out[9] += -0.7071067811865475*Ghat[8]*dv12; 
  out[10] += -0.7071067811865475*Ghat[9]*dv12; 
  out[11] += -0.7071067811865475*Ghat[10]*dv12; 
  out[12] += -1.224744871391589*Ghat[1]*dv12; 
  out[13] += -1.224744871391589*Ghat[2]*dv12; 
  out[14] += -1.224744871391589*Ghat[3]*dv12; 
  out[15] += -1.224744871391589*Ghat[4]*dv12; 
  out[16] += -0.7071067811865475*Ghat[11]*dv12; 
  out[17] += -0.7071067811865475*Ghat[12]*dv12; 
  out[18] += -0.7071067811865475*Ghat[13]*dv12; 
  out[19] += -0.7071067811865475*Ghat[14]*dv12; 
  out[20] += -1.224744871391589*Ghat[5]*dv12; 
  out[21] += -1.224744871391589*Ghat[6]*dv12; 
  out[22] += -1.224744871391589*Ghat[7]*dv12; 
  out[23] += -1.224744871391589*Ghat[8]*dv12; 
  out[24] += -1.224744871391589*Ghat[9]*dv12; 
  out[25] += -1.224744871391589*Ghat[10]*dv12; 
  out[26] += -0.7071067811865475*Ghat[15]*dv12; 
  out[27] += -1.224744871391589*Ghat[11]*dv12; 
  out[28] += -1.224744871391589*Ghat[12]*dv12; 
  out[29] += -1.224744871391589*Ghat[13]*dv12; 
  out[30] += -1.224744871391589*Ghat[14]*dv12; 
  out[31] += -1.224744871391589*Ghat[15]*dv12; 

  } else { 

  if (alpha[8]+alpha[7]-alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[0] = ser_2x3v_p1_surfvz_quad_0(1, fEdge); 
  } else { 

    fUpwindQuad[0] = ser_2x3v_p1_surfvz_quad_0(-1, fSkin); 
  } 
  if ((-alpha[8])+alpha[7]-alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[1] = ser_2x3v_p1_surfvz_quad_1(1, fEdge); 
  } else { 

    fUpwindQuad[1] = ser_2x3v_p1_surfvz_quad_1(-1, fSkin); 
  } 
  if (alpha[8]-alpha[7]-alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[2] = ser_2x3v_p1_surfvz_quad_2(1, fEdge); 
  } else { 

    fUpwindQuad[2] = ser_2x3v_p1_surfvz_quad_2(-1, fSkin); 
  } 
  if ((-alpha[8])-alpha[7]-alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[3] = ser_2x3v_p1_surfvz_quad_3(1, fEdge); 
  } else { 

    fUpwindQuad[3] = ser_2x3v_p1_surfvz_quad_3(-1, fSkin); 
  } 
  if (alpha[8]-alpha[7]-alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[4] = ser_2x3v_p1_surfvz_quad_4(1, fEdge); 
  } else { 

    fUpwindQuad[4] = ser_2x3v_p1_surfvz_quad_4(-1, fSkin); 
  } 
  if ((-alpha[8])-alpha[7]-alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[5] = ser_2x3v_p1_surfvz_quad_5(1, fEdge); 
  } else { 

    fUpwindQuad[5] = ser_2x3v_p1_surfvz_quad_5(-1, fSkin); 
  } 
  if (alpha[8]+alpha[7]-alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[6] = ser_2x3v_p1_surfvz_quad_6(1, fEdge); 
  } else { 

    fUpwindQuad[6] = ser_2x3v_p1_surfvz_quad_6(-1, fSkin); 
  } 
  if ((-alpha[8])+alpha[7]-alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[7] = ser_2x3v_p1_surfvz_quad_7(1, fEdge); 
  } else { 

    fUpwindQuad[7] = ser_2x3v_p1_surfvz_quad_7(-1, fSkin); 
  } 
  if ((-alpha[8])+alpha[7]+alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[8] = ser_2x3v_p1_surfvz_quad_8(1, fEdge); 
  } else { 

    fUpwindQuad[8] = ser_2x3v_p1_surfvz_quad_8(-1, fSkin); 
  } 
  if (alpha[8]+alpha[7]+alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[9] = ser_2x3v_p1_surfvz_quad_9(1, fEdge); 
  } else { 

    fUpwindQuad[9] = ser_2x3v_p1_surfvz_quad_9(-1, fSkin); 
  } 
  if ((-alpha[8])-alpha[7]+alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[10] = ser_2x3v_p1_surfvz_quad_10(1, fEdge); 
  } else { 

    fUpwindQuad[10] = ser_2x3v_p1_surfvz_quad_10(-1, fSkin); 
  } 
  if (alpha[8]-alpha[7]+alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[11] = ser_2x3v_p1_surfvz_quad_11(1, fEdge); 
  } else { 

    fUpwindQuad[11] = ser_2x3v_p1_surfvz_quad_11(-1, fSkin); 
  } 
  if ((-alpha[8])-alpha[7]+alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[12] = ser_2x3v_p1_surfvz_quad_12(1, fEdge); 
  } else { 

    fUpwindQuad[12] = ser_2x3v_p1_surfvz_quad_12(-1, fSkin); 
  } 
  if (alpha[8]-alpha[7]+alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[13] = ser_2x3v_p1_surfvz_quad_13(1, fEdge); 
  } else { 

    fUpwindQuad[13] = ser_2x3v_p1_surfvz_quad_13(-1, fSkin); 
  } 
  if ((-alpha[8])+alpha[7]+alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[14] = ser_2x3v_p1_surfvz_quad_14(1, fEdge); 
  } else { 

    fUpwindQuad[14] = ser_2x3v_p1_surfvz_quad_14(-1, fSkin); 
  } 
  if (alpha[8]+alpha[7]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[15] = ser_2x3v_p1_surfvz_quad_15(1, fEdge); 
  } else { 

    fUpwindQuad[15] = ser_2x3v_p1_surfvz_quad_15(-1, fSkin); 
  } 

  fUpwind[0] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.25*(fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*fUpwindQuad[12]+fUpwindQuad[11]-1.0*fUpwindQuad[10]+fUpwindQuad[9]-1.0*fUpwindQuad[8]+fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[2] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12])+fUpwindQuad[11]+fUpwindQuad[10]-1.0*(fUpwindQuad[9]+fUpwindQuad[8])+fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[3] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]-1.0*(fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8])+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[4] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]-1.0*(fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[5] = 0.25*(fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]+fUpwindQuad[11]-1.0*(fUpwindQuad[10]+fUpwindQuad[9])+fUpwindQuad[8]+fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 
  fUpwind[6] = 0.25*(fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*(fUpwindQuad[12]+fUpwindQuad[11])+fUpwindQuad[10]-1.0*fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[7] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10])+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[8] = 0.25*(fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*fUpwindQuad[12]+fUpwindQuad[11]-1.0*fUpwindQuad[10]+fUpwindQuad[9]-1.0*(fUpwindQuad[8]+fUpwindQuad[7])+fUpwindQuad[6]-1.0*fUpwindQuad[5]+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[9] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12])+fUpwindQuad[11]+fUpwindQuad[10]-1.0*(fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6])+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[10] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]-1.0*(fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[11] = 0.25*(fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]-1.0*fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]-1.0*fUpwindQuad[8]+fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[12] = 0.25*(fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]+fUpwindQuad[11]-1.0*(fUpwindQuad[10]+fUpwindQuad[9])+fUpwindQuad[8]-1.0*fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[13] = 0.25*(fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*(fUpwindQuad[12]+fUpwindQuad[11])+fUpwindQuad[10]-1.0*fUpwindQuad[9]+fUpwindQuad[8]-1.0*fUpwindQuad[7]+fUpwindQuad[6]-1.0*fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[14] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10])+fUpwindQuad[9]+fUpwindQuad[8]-1.0*(fUpwindQuad[7]+fUpwindQuad[6])+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[15] = 0.25*(fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]-1.0*fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]-1.0*(fUpwindQuad[8]+fUpwindQuad[7])+fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 

  Ghat[0] += 0.25*(alpha[8]*fUpwind[8]+alpha[7]*fUpwind[7]+alpha[4]*fUpwind[4]+alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] += 0.25*(alpha[7]*fUpwind[11]+alpha[4]*fUpwind[8]+fUpwind[4]*alpha[8]+alpha[3]*fUpwind[6]+alpha[2]*fUpwind[5]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] += 0.25*(alpha[8]*fUpwind[12]+alpha[4]*fUpwind[9]+alpha[3]*fUpwind[7]+fUpwind[3]*alpha[7]+alpha[1]*fUpwind[5]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] += 0.25*(alpha[8]*fUpwind[13]+alpha[4]*fUpwind[10]+alpha[2]*fUpwind[7]+fUpwind[2]*alpha[7]+alpha[1]*fUpwind[6]+alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]); 
  Ghat[4] += 0.25*(alpha[7]*fUpwind[14]+alpha[3]*fUpwind[10]+alpha[2]*fUpwind[9]+alpha[1]*fUpwind[8]+fUpwind[1]*alpha[8]+alpha[0]*fUpwind[4]+fUpwind[0]*alpha[4]); 
  Ghat[5] += 0.25*(alpha[4]*fUpwind[12]+alpha[3]*fUpwind[11]+alpha[8]*fUpwind[9]+fUpwind[6]*alpha[7]+alpha[0]*fUpwind[5]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[6] += 0.25*(alpha[4]*fUpwind[13]+alpha[2]*fUpwind[11]+alpha[8]*fUpwind[10]+fUpwind[5]*alpha[7]+alpha[0]*fUpwind[6]+alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]); 
  Ghat[7] += 0.25*(alpha[8]*fUpwind[15]+alpha[4]*fUpwind[14]+alpha[1]*fUpwind[11]+alpha[0]*fUpwind[7]+fUpwind[0]*alpha[7]+alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]); 
  Ghat[8] += 0.25*(alpha[7]*fUpwind[15]+alpha[3]*fUpwind[13]+alpha[2]*fUpwind[12]+alpha[0]*fUpwind[8]+fUpwind[0]*alpha[8]+alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]); 
  Ghat[9] += 0.25*(alpha[3]*fUpwind[14]+alpha[1]*fUpwind[12]+alpha[7]*fUpwind[10]+alpha[0]*fUpwind[9]+fUpwind[5]*alpha[8]+alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]); 
  Ghat[10] += 0.25*(alpha[2]*fUpwind[14]+alpha[1]*fUpwind[13]+alpha[0]*fUpwind[10]+alpha[7]*fUpwind[9]+fUpwind[6]*alpha[8]+alpha[3]*fUpwind[4]+fUpwind[3]*alpha[4]); 
  Ghat[11] += 0.25*(alpha[4]*fUpwind[15]+alpha[8]*fUpwind[14]+alpha[0]*fUpwind[11]+alpha[1]*fUpwind[7]+fUpwind[1]*alpha[7]+alpha[2]*fUpwind[6]+alpha[3]*fUpwind[5]); 
  Ghat[12] += 0.25*(alpha[3]*fUpwind[15]+alpha[7]*fUpwind[13]+alpha[0]*fUpwind[12]+alpha[1]*fUpwind[9]+alpha[2]*fUpwind[8]+fUpwind[2]*alpha[8]+alpha[4]*fUpwind[5]); 
  Ghat[13] += 0.25*(alpha[2]*fUpwind[15]+alpha[0]*fUpwind[13]+alpha[7]*fUpwind[12]+alpha[1]*fUpwind[10]+alpha[3]*fUpwind[8]+fUpwind[3]*alpha[8]+alpha[4]*fUpwind[6]); 
  Ghat[14] += 0.25*(alpha[1]*fUpwind[15]+alpha[0]*fUpwind[14]+alpha[8]*fUpwind[11]+alpha[2]*fUpwind[10]+alpha[3]*fUpwind[9]+alpha[4]*fUpwind[7]+fUpwind[4]*alpha[7]); 
  Ghat[15] += 0.25*(alpha[0]*fUpwind[15]+alpha[1]*fUpwind[14]+alpha[2]*fUpwind[13]+alpha[3]*fUpwind[12]+alpha[4]*fUpwind[11]+alpha[7]*fUpwind[8]+fUpwind[7]*alpha[8]); 

  out[0] += 0.7071067811865475*Ghat[0]*dv12; 
  out[1] += 0.7071067811865475*Ghat[1]*dv12; 
  out[2] += 0.7071067811865475*Ghat[2]*dv12; 
  out[3] += 0.7071067811865475*Ghat[3]*dv12; 
  out[4] += 0.7071067811865475*Ghat[4]*dv12; 
  out[5] += -1.224744871391589*Ghat[0]*dv12; 
  out[6] += 0.7071067811865475*Ghat[5]*dv12; 
  out[7] += 0.7071067811865475*Ghat[6]*dv12; 
  out[8] += 0.7071067811865475*Ghat[7]*dv12; 
  out[9] += 0.7071067811865475*Ghat[8]*dv12; 
  out[10] += 0.7071067811865475*Ghat[9]*dv12; 
  out[11] += 0.7071067811865475*Ghat[10]*dv12; 
  out[12] += -1.224744871391589*Ghat[1]*dv12; 
  out[13] += -1.224744871391589*Ghat[2]*dv12; 
  out[14] += -1.224744871391589*Ghat[3]*dv12; 
  out[15] += -1.224744871391589*Ghat[4]*dv12; 
  out[16] += 0.7071067811865475*Ghat[11]*dv12; 
  out[17] += 0.7071067811865475*Ghat[12]*dv12; 
  out[18] += 0.7071067811865475*Ghat[13]*dv12; 
  out[19] += 0.7071067811865475*Ghat[14]*dv12; 
  out[20] += -1.224744871391589*Ghat[5]*dv12; 
  out[21] += -1.224744871391589*Ghat[6]*dv12; 
  out[22] += -1.224744871391589*Ghat[7]*dv12; 
  out[23] += -1.224744871391589*Ghat[8]*dv12; 
  out[24] += -1.224744871391589*Ghat[9]*dv12; 
  out[25] += -1.224744871391589*Ghat[10]*dv12; 
  out[26] += 0.7071067811865475*Ghat[15]*dv12; 
  out[27] += -1.224744871391589*Ghat[11]*dv12; 
  out[28] += -1.224744871391589*Ghat[12]*dv12; 
  out[29] += -1.224744871391589*Ghat[13]*dv12; 
  out[30] += -1.224744871391589*Ghat[14]*dv12; 
  out[31] += -1.224744871391589*Ghat[15]*dv12; 

  } 
} 