#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_ser_5x_p1_surfx4_quad.h> 
#include <gkyl_basis_ser_5x_p1_upwind.h> 
GKYL_CU_DH void lbo_vlasov_drag_surfvy_2x3v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[5]:         cell-center coordinates. 
  // dxv[5]:       cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[12]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[4]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      distribution function in cells 
  // out:           incremented distribution function in cell 
  double rdv2 = 2.0/dxv[3]; 

  const double *sumNuUy = &nuUSum[4]; 

  double alphaDrSurf_l[16] = {0.0}; 
  alphaDrSurf_l[0] = 2.0*nuSum[0]*w[3]-1.0*nuSum[0]*dxv[3]-2.0*sumNuUy[0]; 
  alphaDrSurf_l[1] = 2.0*nuSum[1]*w[3]-1.0*nuSum[1]*dxv[3]-2.0*sumNuUy[1]; 
  alphaDrSurf_l[2] = 2.0*nuSum[2]*w[3]-1.0*nuSum[2]*dxv[3]-2.0*sumNuUy[2]; 
  alphaDrSurf_l[5] = 2.0*nuSum[3]*w[3]-2.0*sumNuUy[3]-1.0*dxv[3]*nuSum[3]; 

  double alphaDrSurf_r[16] = {0.0}; 
  alphaDrSurf_r[0] = 2.0*nuSum[0]*w[3]+nuSum[0]*dxv[3]-2.0*sumNuUy[0]; 
  alphaDrSurf_r[1] = 2.0*nuSum[1]*w[3]+nuSum[1]*dxv[3]-2.0*sumNuUy[1]; 
  alphaDrSurf_r[2] = 2.0*nuSum[2]*w[3]+nuSum[2]*dxv[3]-2.0*sumNuUy[2]; 
  alphaDrSurf_r[5] = 2.0*nuSum[3]*w[3]-2.0*sumNuUy[3]+dxv[3]*nuSum[3]; 

  double fUpwindQuad_l[16] = {0.0};
  double fUpwindQuad_r[16] = {0.0};
  double fUpwind_l[16] = {0.0};
  double fUpwind_r[16] = {0.0};
  double Gdrag_l[16] = {0.0}; 
  double Gdrag_r[16] = {0.0}; 

  if (alphaDrSurf_l[5]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[0] = ser_5x_p1_surfx4_quad_0_r(fl); 
    fUpwindQuad_l[4] = ser_5x_p1_surfx4_quad_4_r(fl); 
    fUpwindQuad_l[8] = ser_5x_p1_surfx4_quad_8_r(fl); 
    fUpwindQuad_l[12] = ser_5x_p1_surfx4_quad_12_r(fl); 
  } else { 
    fUpwindQuad_l[0] = ser_5x_p1_surfx4_quad_0_l(fc); 
    fUpwindQuad_l[4] = ser_5x_p1_surfx4_quad_4_l(fc); 
    fUpwindQuad_l[8] = ser_5x_p1_surfx4_quad_8_l(fc); 
    fUpwindQuad_l[12] = ser_5x_p1_surfx4_quad_12_l(fc); 
  } 
  if (alphaDrSurf_r[5]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[0] = ser_5x_p1_surfx4_quad_0_r(fc); 
    fUpwindQuad_r[4] = ser_5x_p1_surfx4_quad_4_r(fc); 
    fUpwindQuad_r[8] = ser_5x_p1_surfx4_quad_8_r(fc); 
    fUpwindQuad_r[12] = ser_5x_p1_surfx4_quad_12_r(fc); 
  } else { 
    fUpwindQuad_r[0] = ser_5x_p1_surfx4_quad_0_l(fr); 
    fUpwindQuad_r[4] = ser_5x_p1_surfx4_quad_4_l(fr); 
    fUpwindQuad_r[8] = ser_5x_p1_surfx4_quad_8_l(fr); 
    fUpwindQuad_r[12] = ser_5x_p1_surfx4_quad_12_l(fr); 
  } 
  if ((-alphaDrSurf_l[5])-alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[1] = ser_5x_p1_surfx4_quad_1_r(fl); 
    fUpwindQuad_l[5] = ser_5x_p1_surfx4_quad_5_r(fl); 
    fUpwindQuad_l[9] = ser_5x_p1_surfx4_quad_9_r(fl); 
    fUpwindQuad_l[13] = ser_5x_p1_surfx4_quad_13_r(fl); 
  } else { 
    fUpwindQuad_l[1] = ser_5x_p1_surfx4_quad_1_l(fc); 
    fUpwindQuad_l[5] = ser_5x_p1_surfx4_quad_5_l(fc); 
    fUpwindQuad_l[9] = ser_5x_p1_surfx4_quad_9_l(fc); 
    fUpwindQuad_l[13] = ser_5x_p1_surfx4_quad_13_l(fc); 
  } 
  if ((-alphaDrSurf_r[5])-alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[1] = ser_5x_p1_surfx4_quad_1_r(fc); 
    fUpwindQuad_r[5] = ser_5x_p1_surfx4_quad_5_r(fc); 
    fUpwindQuad_r[9] = ser_5x_p1_surfx4_quad_9_r(fc); 
    fUpwindQuad_r[13] = ser_5x_p1_surfx4_quad_13_r(fc); 
  } else { 
    fUpwindQuad_r[1] = ser_5x_p1_surfx4_quad_1_l(fr); 
    fUpwindQuad_r[5] = ser_5x_p1_surfx4_quad_5_l(fr); 
    fUpwindQuad_r[9] = ser_5x_p1_surfx4_quad_9_l(fr); 
    fUpwindQuad_r[13] = ser_5x_p1_surfx4_quad_13_l(fr); 
  } 
  if ((-alphaDrSurf_l[5])+alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[2] = ser_5x_p1_surfx4_quad_2_r(fl); 
    fUpwindQuad_l[6] = ser_5x_p1_surfx4_quad_6_r(fl); 
    fUpwindQuad_l[10] = ser_5x_p1_surfx4_quad_10_r(fl); 
    fUpwindQuad_l[14] = ser_5x_p1_surfx4_quad_14_r(fl); 
  } else { 
    fUpwindQuad_l[2] = ser_5x_p1_surfx4_quad_2_l(fc); 
    fUpwindQuad_l[6] = ser_5x_p1_surfx4_quad_6_l(fc); 
    fUpwindQuad_l[10] = ser_5x_p1_surfx4_quad_10_l(fc); 
    fUpwindQuad_l[14] = ser_5x_p1_surfx4_quad_14_l(fc); 
  } 
  if ((-alphaDrSurf_r[5])+alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[2] = ser_5x_p1_surfx4_quad_2_r(fc); 
    fUpwindQuad_r[6] = ser_5x_p1_surfx4_quad_6_r(fc); 
    fUpwindQuad_r[10] = ser_5x_p1_surfx4_quad_10_r(fc); 
    fUpwindQuad_r[14] = ser_5x_p1_surfx4_quad_14_r(fc); 
  } else { 
    fUpwindQuad_r[2] = ser_5x_p1_surfx4_quad_2_l(fr); 
    fUpwindQuad_r[6] = ser_5x_p1_surfx4_quad_6_l(fr); 
    fUpwindQuad_r[10] = ser_5x_p1_surfx4_quad_10_l(fr); 
    fUpwindQuad_r[14] = ser_5x_p1_surfx4_quad_14_l(fr); 
  } 
  if (alphaDrSurf_l[5]+alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[3] = ser_5x_p1_surfx4_quad_3_r(fl); 
    fUpwindQuad_l[7] = ser_5x_p1_surfx4_quad_7_r(fl); 
    fUpwindQuad_l[11] = ser_5x_p1_surfx4_quad_11_r(fl); 
    fUpwindQuad_l[15] = ser_5x_p1_surfx4_quad_15_r(fl); 
  } else { 
    fUpwindQuad_l[3] = ser_5x_p1_surfx4_quad_3_l(fc); 
    fUpwindQuad_l[7] = ser_5x_p1_surfx4_quad_7_l(fc); 
    fUpwindQuad_l[11] = ser_5x_p1_surfx4_quad_11_l(fc); 
    fUpwindQuad_l[15] = ser_5x_p1_surfx4_quad_15_l(fc); 
  } 
  if (alphaDrSurf_r[5]+alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[3] = ser_5x_p1_surfx4_quad_3_r(fc); 
    fUpwindQuad_r[7] = ser_5x_p1_surfx4_quad_7_r(fc); 
    fUpwindQuad_r[11] = ser_5x_p1_surfx4_quad_11_r(fc); 
    fUpwindQuad_r[15] = ser_5x_p1_surfx4_quad_15_r(fc); 
  } else { 
    fUpwindQuad_r[3] = ser_5x_p1_surfx4_quad_3_l(fr); 
    fUpwindQuad_r[7] = ser_5x_p1_surfx4_quad_7_l(fr); 
    fUpwindQuad_r[11] = ser_5x_p1_surfx4_quad_11_l(fr); 
    fUpwindQuad_r[15] = ser_5x_p1_surfx4_quad_15_l(fr); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_5x_p1_upwind(fUpwindQuad_l, fUpwind_l); 
  ser_5x_p1_upwind(fUpwindQuad_r, fUpwind_r); 

  Gdrag_l[0] = 0.25*alphaDrSurf_l[5]*fUpwind_l[5]+0.25*alphaDrSurf_l[2]*fUpwind_l[2]+0.25*alphaDrSurf_l[1]*fUpwind_l[1]+0.25*alphaDrSurf_l[0]*fUpwind_l[0]; 
  Gdrag_l[1] = 0.25*alphaDrSurf_l[2]*fUpwind_l[5]+0.25*fUpwind_l[2]*alphaDrSurf_l[5]+0.25*alphaDrSurf_l[0]*fUpwind_l[1]+0.25*fUpwind_l[0]*alphaDrSurf_l[1]; 
  Gdrag_l[2] = 0.25*alphaDrSurf_l[1]*fUpwind_l[5]+0.25*fUpwind_l[1]*alphaDrSurf_l[5]+0.25*alphaDrSurf_l[0]*fUpwind_l[2]+0.25*fUpwind_l[0]*alphaDrSurf_l[2]; 
  Gdrag_l[3] = 0.25*alphaDrSurf_l[5]*fUpwind_l[11]+0.25*alphaDrSurf_l[2]*fUpwind_l[7]+0.25*alphaDrSurf_l[1]*fUpwind_l[6]+0.25*alphaDrSurf_l[0]*fUpwind_l[3]; 
  Gdrag_l[4] = 0.25*alphaDrSurf_l[5]*fUpwind_l[12]+0.25*alphaDrSurf_l[2]*fUpwind_l[9]+0.25*alphaDrSurf_l[1]*fUpwind_l[8]+0.25*alphaDrSurf_l[0]*fUpwind_l[4]; 
  Gdrag_l[5] = 0.25*alphaDrSurf_l[0]*fUpwind_l[5]+0.25*fUpwind_l[0]*alphaDrSurf_l[5]+0.25*alphaDrSurf_l[1]*fUpwind_l[2]+0.25*fUpwind_l[1]*alphaDrSurf_l[2]; 
  Gdrag_l[6] = 0.25*alphaDrSurf_l[2]*fUpwind_l[11]+0.25*alphaDrSurf_l[5]*fUpwind_l[7]+0.25*alphaDrSurf_l[0]*fUpwind_l[6]+0.25*alphaDrSurf_l[1]*fUpwind_l[3]; 
  Gdrag_l[7] = 0.25*alphaDrSurf_l[1]*fUpwind_l[11]+0.25*alphaDrSurf_l[0]*fUpwind_l[7]+0.25*alphaDrSurf_l[5]*fUpwind_l[6]+0.25*alphaDrSurf_l[2]*fUpwind_l[3]; 
  Gdrag_l[8] = 0.25*alphaDrSurf_l[2]*fUpwind_l[12]+0.25*alphaDrSurf_l[5]*fUpwind_l[9]+0.25*alphaDrSurf_l[0]*fUpwind_l[8]+0.25*alphaDrSurf_l[1]*fUpwind_l[4]; 
  Gdrag_l[9] = 0.25*alphaDrSurf_l[1]*fUpwind_l[12]+0.25*alphaDrSurf_l[0]*fUpwind_l[9]+0.25*alphaDrSurf_l[5]*fUpwind_l[8]+0.25*alphaDrSurf_l[2]*fUpwind_l[4]; 
  Gdrag_l[10] = 0.25*alphaDrSurf_l[5]*fUpwind_l[15]+0.25*alphaDrSurf_l[2]*fUpwind_l[14]+0.25*alphaDrSurf_l[1]*fUpwind_l[13]+0.25*alphaDrSurf_l[0]*fUpwind_l[10]; 
  Gdrag_l[11] = 0.25*alphaDrSurf_l[0]*fUpwind_l[11]+0.25*alphaDrSurf_l[1]*fUpwind_l[7]+0.25*alphaDrSurf_l[2]*fUpwind_l[6]+0.25*fUpwind_l[3]*alphaDrSurf_l[5]; 
  Gdrag_l[12] = 0.25*alphaDrSurf_l[0]*fUpwind_l[12]+0.25*alphaDrSurf_l[1]*fUpwind_l[9]+0.25*alphaDrSurf_l[2]*fUpwind_l[8]+0.25*fUpwind_l[4]*alphaDrSurf_l[5]; 
  Gdrag_l[13] = 0.25*alphaDrSurf_l[2]*fUpwind_l[15]+0.25*alphaDrSurf_l[5]*fUpwind_l[14]+0.25*alphaDrSurf_l[0]*fUpwind_l[13]+0.25*alphaDrSurf_l[1]*fUpwind_l[10]; 
  Gdrag_l[14] = 0.25*alphaDrSurf_l[1]*fUpwind_l[15]+0.25*alphaDrSurf_l[0]*fUpwind_l[14]+0.25*alphaDrSurf_l[5]*fUpwind_l[13]+0.25*alphaDrSurf_l[2]*fUpwind_l[10]; 
  Gdrag_l[15] = 0.25*alphaDrSurf_l[0]*fUpwind_l[15]+0.25*alphaDrSurf_l[1]*fUpwind_l[14]+0.25*alphaDrSurf_l[2]*fUpwind_l[13]+0.25*alphaDrSurf_l[5]*fUpwind_l[10]; 

  Gdrag_r[0] = 0.25*alphaDrSurf_r[5]*fUpwind_r[5]+0.25*alphaDrSurf_r[2]*fUpwind_r[2]+0.25*alphaDrSurf_r[1]*fUpwind_r[1]+0.25*alphaDrSurf_r[0]*fUpwind_r[0]; 
  Gdrag_r[1] = 0.25*alphaDrSurf_r[2]*fUpwind_r[5]+0.25*fUpwind_r[2]*alphaDrSurf_r[5]+0.25*alphaDrSurf_r[0]*fUpwind_r[1]+0.25*fUpwind_r[0]*alphaDrSurf_r[1]; 
  Gdrag_r[2] = 0.25*alphaDrSurf_r[1]*fUpwind_r[5]+0.25*fUpwind_r[1]*alphaDrSurf_r[5]+0.25*alphaDrSurf_r[0]*fUpwind_r[2]+0.25*fUpwind_r[0]*alphaDrSurf_r[2]; 
  Gdrag_r[3] = 0.25*alphaDrSurf_r[5]*fUpwind_r[11]+0.25*alphaDrSurf_r[2]*fUpwind_r[7]+0.25*alphaDrSurf_r[1]*fUpwind_r[6]+0.25*alphaDrSurf_r[0]*fUpwind_r[3]; 
  Gdrag_r[4] = 0.25*alphaDrSurf_r[5]*fUpwind_r[12]+0.25*alphaDrSurf_r[2]*fUpwind_r[9]+0.25*alphaDrSurf_r[1]*fUpwind_r[8]+0.25*alphaDrSurf_r[0]*fUpwind_r[4]; 
  Gdrag_r[5] = 0.25*alphaDrSurf_r[0]*fUpwind_r[5]+0.25*fUpwind_r[0]*alphaDrSurf_r[5]+0.25*alphaDrSurf_r[1]*fUpwind_r[2]+0.25*fUpwind_r[1]*alphaDrSurf_r[2]; 
  Gdrag_r[6] = 0.25*alphaDrSurf_r[2]*fUpwind_r[11]+0.25*alphaDrSurf_r[5]*fUpwind_r[7]+0.25*alphaDrSurf_r[0]*fUpwind_r[6]+0.25*alphaDrSurf_r[1]*fUpwind_r[3]; 
  Gdrag_r[7] = 0.25*alphaDrSurf_r[1]*fUpwind_r[11]+0.25*alphaDrSurf_r[0]*fUpwind_r[7]+0.25*alphaDrSurf_r[5]*fUpwind_r[6]+0.25*alphaDrSurf_r[2]*fUpwind_r[3]; 
  Gdrag_r[8] = 0.25*alphaDrSurf_r[2]*fUpwind_r[12]+0.25*alphaDrSurf_r[5]*fUpwind_r[9]+0.25*alphaDrSurf_r[0]*fUpwind_r[8]+0.25*alphaDrSurf_r[1]*fUpwind_r[4]; 
  Gdrag_r[9] = 0.25*alphaDrSurf_r[1]*fUpwind_r[12]+0.25*alphaDrSurf_r[0]*fUpwind_r[9]+0.25*alphaDrSurf_r[5]*fUpwind_r[8]+0.25*alphaDrSurf_r[2]*fUpwind_r[4]; 
  Gdrag_r[10] = 0.25*alphaDrSurf_r[5]*fUpwind_r[15]+0.25*alphaDrSurf_r[2]*fUpwind_r[14]+0.25*alphaDrSurf_r[1]*fUpwind_r[13]+0.25*alphaDrSurf_r[0]*fUpwind_r[10]; 
  Gdrag_r[11] = 0.25*alphaDrSurf_r[0]*fUpwind_r[11]+0.25*alphaDrSurf_r[1]*fUpwind_r[7]+0.25*alphaDrSurf_r[2]*fUpwind_r[6]+0.25*fUpwind_r[3]*alphaDrSurf_r[5]; 
  Gdrag_r[12] = 0.25*alphaDrSurf_r[0]*fUpwind_r[12]+0.25*alphaDrSurf_r[1]*fUpwind_r[9]+0.25*alphaDrSurf_r[2]*fUpwind_r[8]+0.25*fUpwind_r[4]*alphaDrSurf_r[5]; 
  Gdrag_r[13] = 0.25*alphaDrSurf_r[2]*fUpwind_r[15]+0.25*alphaDrSurf_r[5]*fUpwind_r[14]+0.25*alphaDrSurf_r[0]*fUpwind_r[13]+0.25*alphaDrSurf_r[1]*fUpwind_r[10]; 
  Gdrag_r[14] = 0.25*alphaDrSurf_r[1]*fUpwind_r[15]+0.25*alphaDrSurf_r[0]*fUpwind_r[14]+0.25*alphaDrSurf_r[5]*fUpwind_r[13]+0.25*alphaDrSurf_r[2]*fUpwind_r[10]; 
  Gdrag_r[15] = 0.25*alphaDrSurf_r[0]*fUpwind_r[15]+0.25*alphaDrSurf_r[1]*fUpwind_r[14]+0.25*alphaDrSurf_r[2]*fUpwind_r[13]+0.25*alphaDrSurf_r[5]*fUpwind_r[10]; 

  out[0] += 0.7071067811865475*Gdrag_r[0]*rdv2-0.7071067811865475*Gdrag_l[0]*rdv2; 
  out[1] += 0.7071067811865475*Gdrag_r[1]*rdv2-0.7071067811865475*Gdrag_l[1]*rdv2; 
  out[2] += 0.7071067811865475*Gdrag_r[2]*rdv2-0.7071067811865475*Gdrag_l[2]*rdv2; 
  out[3] += 0.7071067811865475*Gdrag_r[3]*rdv2-0.7071067811865475*Gdrag_l[3]*rdv2; 
  out[4] += 1.224744871391589*Gdrag_r[0]*rdv2+1.224744871391589*Gdrag_l[0]*rdv2; 
  out[5] += 0.7071067811865475*Gdrag_r[4]*rdv2-0.7071067811865475*Gdrag_l[4]*rdv2; 
  out[6] += 0.7071067811865475*Gdrag_r[5]*rdv2-0.7071067811865475*Gdrag_l[5]*rdv2; 
  out[7] += 0.7071067811865475*Gdrag_r[6]*rdv2-0.7071067811865475*Gdrag_l[6]*rdv2; 
  out[8] += 0.7071067811865475*Gdrag_r[7]*rdv2-0.7071067811865475*Gdrag_l[7]*rdv2; 
  out[9] += 1.224744871391589*Gdrag_r[1]*rdv2+1.224744871391589*Gdrag_l[1]*rdv2; 
  out[10] += 1.224744871391589*Gdrag_r[2]*rdv2+1.224744871391589*Gdrag_l[2]*rdv2; 
  out[11] += 1.224744871391589*Gdrag_r[3]*rdv2+1.224744871391589*Gdrag_l[3]*rdv2; 
  out[12] += 0.7071067811865475*Gdrag_r[8]*rdv2-0.7071067811865475*Gdrag_l[8]*rdv2; 
  out[13] += 0.7071067811865475*Gdrag_r[9]*rdv2-0.7071067811865475*Gdrag_l[9]*rdv2; 
  out[14] += 0.7071067811865475*Gdrag_r[10]*rdv2-0.7071067811865475*Gdrag_l[10]*rdv2; 
  out[15] += 1.224744871391589*Gdrag_r[4]*rdv2+1.224744871391589*Gdrag_l[4]*rdv2; 
  out[16] += 0.7071067811865475*Gdrag_r[11]*rdv2-0.7071067811865475*Gdrag_l[11]*rdv2; 
  out[17] += 1.224744871391589*Gdrag_r[5]*rdv2+1.224744871391589*Gdrag_l[5]*rdv2; 
  out[18] += 1.224744871391589*Gdrag_r[6]*rdv2+1.224744871391589*Gdrag_l[6]*rdv2; 
  out[19] += 1.224744871391589*Gdrag_r[7]*rdv2+1.224744871391589*Gdrag_l[7]*rdv2; 
  out[20] += 0.7071067811865475*Gdrag_r[12]*rdv2-0.7071067811865475*Gdrag_l[12]*rdv2; 
  out[21] += 0.7071067811865475*Gdrag_r[13]*rdv2-0.7071067811865475*Gdrag_l[13]*rdv2; 
  out[22] += 0.7071067811865475*Gdrag_r[14]*rdv2-0.7071067811865475*Gdrag_l[14]*rdv2; 
  out[23] += 1.224744871391589*Gdrag_r[8]*rdv2+1.224744871391589*Gdrag_l[8]*rdv2; 
  out[24] += 1.224744871391589*Gdrag_r[9]*rdv2+1.224744871391589*Gdrag_l[9]*rdv2; 
  out[25] += 1.224744871391589*Gdrag_r[10]*rdv2+1.224744871391589*Gdrag_l[10]*rdv2; 
  out[26] += 1.224744871391589*Gdrag_r[11]*rdv2+1.224744871391589*Gdrag_l[11]*rdv2; 
  out[27] += 0.7071067811865475*Gdrag_r[15]*rdv2-0.7071067811865475*Gdrag_l[15]*rdv2; 
  out[28] += 1.224744871391589*Gdrag_r[12]*rdv2+1.224744871391589*Gdrag_l[12]*rdv2; 
  out[29] += 1.224744871391589*Gdrag_r[13]*rdv2+1.224744871391589*Gdrag_l[13]*rdv2; 
  out[30] += 1.224744871391589*Gdrag_r[14]*rdv2+1.224744871391589*Gdrag_l[14]*rdv2; 
  out[31] += 1.224744871391589*Gdrag_r[15]*rdv2+1.224744871391589*Gdrag_l[15]*rdv2; 
} 
