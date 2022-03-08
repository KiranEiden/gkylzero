#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_4x_p1_surfx2_quad.h> 
#include <gkyl_basis_ser_4x_p1_upwind.h> 
GKYL_CU_DH void vlasov_boundary_surfvx_1x3v_ser_p1(const double *w, const double *dxv, const double *qmem, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // qmem:        q/m*EM fields.
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv10 = 2/dxv[1]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double dv3 = dxv[3], wv3 = w[3]; 
  const double *E0 = &qmem[0]; 
  const double *B0 = &qmem[6]; 
  const double *B1 = &qmem[8]; 
  const double *B2 = &qmem[10]; 

  double alpha[8] = {0.0}; 

  alpha[0] = 2.0*(B2[0]*wv2+E0[0])-2.0*B1[0]*wv3; 
  alpha[1] = 2.0*(B2[1]*wv2+E0[1])-2.0*B1[1]*wv3; 
  alpha[2] = 0.5773502691896258*B2[0]*dv2; 
  alpha[3] = -0.5773502691896258*B1[0]*dv3; 
  alpha[4] = 0.5773502691896258*B2[1]*dv2; 
  alpha[5] = -0.5773502691896258*B1[1]*dv3; 

  double fUpwindQuad[8] = {0.0};
  double fUpwind[8] = {0.0};
  double Ghat[8] = {0.0}; 

  if (edge == -1) { 

  if (alpha[5]+alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[0] = ser_4x_p1_surfx2_quad_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_4x_p1_surfx2_quad_0_l(fEdge); 
  } 
  if ((-alpha[5])-alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[1] = ser_4x_p1_surfx2_quad_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = ser_4x_p1_surfx2_quad_1_l(fEdge); 
  } 
  if (alpha[5]-alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[2] = ser_4x_p1_surfx2_quad_2_r(fSkin); 
  } else { 
    fUpwindQuad[2] = ser_4x_p1_surfx2_quad_2_l(fEdge); 
  } 
  if ((-alpha[5])+alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[3] = ser_4x_p1_surfx2_quad_3_r(fSkin); 
  } else { 
    fUpwindQuad[3] = ser_4x_p1_surfx2_quad_3_l(fEdge); 
  } 
  if ((-alpha[5])+alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[4] = ser_4x_p1_surfx2_quad_4_r(fSkin); 
  } else { 
    fUpwindQuad[4] = ser_4x_p1_surfx2_quad_4_l(fEdge); 
  } 
  if (alpha[5]-alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[5] = ser_4x_p1_surfx2_quad_5_r(fSkin); 
  } else { 
    fUpwindQuad[5] = ser_4x_p1_surfx2_quad_5_l(fEdge); 
  } 
  if ((-alpha[5])-alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[6] = ser_4x_p1_surfx2_quad_6_r(fSkin); 
  } else { 
    fUpwindQuad[6] = ser_4x_p1_surfx2_quad_6_l(fEdge); 
  } 
  if (alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[7] = ser_4x_p1_surfx2_quad_7_r(fSkin); 
  } else { 
    fUpwindQuad[7] = ser_4x_p1_surfx2_quad_7_l(fEdge); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_4x_p1_upwind(fUpwindQuad, fUpwind); 

  Ghat[0] += 0.3535533905932737*(alpha[5]*fUpwind[5]+alpha[4]*fUpwind[4]+alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] += 0.3535533905932737*(alpha[3]*fUpwind[5]+fUpwind[3]*alpha[5]+alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] += 0.3535533905932737*(alpha[5]*fUpwind[7]+alpha[3]*fUpwind[6]+alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] += 0.3535533905932737*(alpha[4]*fUpwind[7]+alpha[2]*fUpwind[6]+alpha[1]*fUpwind[5]+fUpwind[1]*alpha[5]+alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]); 
  Ghat[4] += 0.3535533905932737*(alpha[3]*fUpwind[7]+alpha[5]*fUpwind[6]+alpha[0]*fUpwind[4]+fUpwind[0]*alpha[4]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[5] += 0.3535533905932737*(alpha[2]*fUpwind[7]+alpha[4]*fUpwind[6]+alpha[0]*fUpwind[5]+fUpwind[0]*alpha[5]+alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]); 
  Ghat[6] += 0.3535533905932737*(alpha[1]*fUpwind[7]+alpha[0]*fUpwind[6]+alpha[4]*fUpwind[5]+fUpwind[4]*alpha[5]+alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]); 
  Ghat[7] += 0.3535533905932737*(alpha[0]*fUpwind[7]+alpha[1]*fUpwind[6]+alpha[2]*fUpwind[5]+fUpwind[2]*alpha[5]+alpha[3]*fUpwind[4]+fUpwind[3]*alpha[4]); 

  out[0] += -0.7071067811865475*Ghat[0]*dv10; 
  out[1] += -0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -1.224744871391589*Ghat[0]*dv10; 
  out[3] += -0.7071067811865475*Ghat[2]*dv10; 
  out[4] += -0.7071067811865475*Ghat[3]*dv10; 
  out[5] += -1.224744871391589*Ghat[1]*dv10; 
  out[6] += -0.7071067811865475*Ghat[4]*dv10; 
  out[7] += -1.224744871391589*Ghat[2]*dv10; 
  out[8] += -0.7071067811865475*Ghat[5]*dv10; 
  out[9] += -1.224744871391589*Ghat[3]*dv10; 
  out[10] += -0.7071067811865475*Ghat[6]*dv10; 
  out[11] += -1.224744871391589*Ghat[4]*dv10; 
  out[12] += -1.224744871391589*Ghat[5]*dv10; 
  out[13] += -0.7071067811865475*Ghat[7]*dv10; 
  out[14] += -1.224744871391589*Ghat[6]*dv10; 
  out[15] += -1.224744871391589*Ghat[7]*dv10; 

  } else { 

  if (alpha[5]+alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[0] = ser_4x_p1_surfx2_quad_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_4x_p1_surfx2_quad_0_l(fSkin); 
  } 
  if ((-alpha[5])-alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[1] = ser_4x_p1_surfx2_quad_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = ser_4x_p1_surfx2_quad_1_l(fSkin); 
  } 
  if (alpha[5]-alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[2] = ser_4x_p1_surfx2_quad_2_r(fEdge); 
  } else { 
    fUpwindQuad[2] = ser_4x_p1_surfx2_quad_2_l(fSkin); 
  } 
  if ((-alpha[5])+alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[3] = ser_4x_p1_surfx2_quad_3_r(fEdge); 
  } else { 
    fUpwindQuad[3] = ser_4x_p1_surfx2_quad_3_l(fSkin); 
  } 
  if ((-alpha[5])+alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[4] = ser_4x_p1_surfx2_quad_4_r(fEdge); 
  } else { 
    fUpwindQuad[4] = ser_4x_p1_surfx2_quad_4_l(fSkin); 
  } 
  if (alpha[5]-alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[5] = ser_4x_p1_surfx2_quad_5_r(fEdge); 
  } else { 
    fUpwindQuad[5] = ser_4x_p1_surfx2_quad_5_l(fSkin); 
  } 
  if ((-alpha[5])-alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[6] = ser_4x_p1_surfx2_quad_6_r(fEdge); 
  } else { 
    fUpwindQuad[6] = ser_4x_p1_surfx2_quad_6_l(fSkin); 
  } 
  if (alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[7] = ser_4x_p1_surfx2_quad_7_r(fEdge); 
  } else { 
    fUpwindQuad[7] = ser_4x_p1_surfx2_quad_7_l(fSkin); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_4x_p1_upwind(fUpwindQuad, fUpwind); 

  Ghat[0] += 0.3535533905932737*(alpha[5]*fUpwind[5]+alpha[4]*fUpwind[4]+alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] += 0.3535533905932737*(alpha[3]*fUpwind[5]+fUpwind[3]*alpha[5]+alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] += 0.3535533905932737*(alpha[5]*fUpwind[7]+alpha[3]*fUpwind[6]+alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] += 0.3535533905932737*(alpha[4]*fUpwind[7]+alpha[2]*fUpwind[6]+alpha[1]*fUpwind[5]+fUpwind[1]*alpha[5]+alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]); 
  Ghat[4] += 0.3535533905932737*(alpha[3]*fUpwind[7]+alpha[5]*fUpwind[6]+alpha[0]*fUpwind[4]+fUpwind[0]*alpha[4]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[5] += 0.3535533905932737*(alpha[2]*fUpwind[7]+alpha[4]*fUpwind[6]+alpha[0]*fUpwind[5]+fUpwind[0]*alpha[5]+alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]); 
  Ghat[6] += 0.3535533905932737*(alpha[1]*fUpwind[7]+alpha[0]*fUpwind[6]+alpha[4]*fUpwind[5]+fUpwind[4]*alpha[5]+alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]); 
  Ghat[7] += 0.3535533905932737*(alpha[0]*fUpwind[7]+alpha[1]*fUpwind[6]+alpha[2]*fUpwind[5]+fUpwind[2]*alpha[5]+alpha[3]*fUpwind[4]+fUpwind[3]*alpha[4]); 

  out[0] += 0.7071067811865475*Ghat[0]*dv10; 
  out[1] += 0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -1.224744871391589*Ghat[0]*dv10; 
  out[3] += 0.7071067811865475*Ghat[2]*dv10; 
  out[4] += 0.7071067811865475*Ghat[3]*dv10; 
  out[5] += -1.224744871391589*Ghat[1]*dv10; 
  out[6] += 0.7071067811865475*Ghat[4]*dv10; 
  out[7] += -1.224744871391589*Ghat[2]*dv10; 
  out[8] += 0.7071067811865475*Ghat[5]*dv10; 
  out[9] += -1.224744871391589*Ghat[3]*dv10; 
  out[10] += 0.7071067811865475*Ghat[6]*dv10; 
  out[11] += -1.224744871391589*Ghat[4]*dv10; 
  out[12] += -1.224744871391589*Ghat[5]*dv10; 
  out[13] += 0.7071067811865475*Ghat[7]*dv10; 
  out[14] += -1.224744871391589*Ghat[6]*dv10; 
  out[15] += -1.224744871391589*Ghat[7]*dv10; 

  } 
} 
