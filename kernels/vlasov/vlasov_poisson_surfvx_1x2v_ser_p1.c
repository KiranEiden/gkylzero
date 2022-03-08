#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_3x_p1_surfx2_quad.h> 
#include <gkyl_basis_ser_3x_p1_upwind.h> 
GKYL_CU_DH void vlasov_poisson_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double *fac_phi, const double *vecA, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fac_phi:   potential (scaled by appropriate factors).
  // vecA:      vector potential (scaled by appropriate factors). Unused in pure Vlasov-Poisson. 
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 
  const double dv10 = 2/dxv[1]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double *phi = &fac_phi[0]; 
  const double dx10 = 2/dxv[0]; 
  double alpha[4] = {0.0}; 

  alpha[0] = -2.449489742783178*phi[1]*dx10; 

  double fUpwindQuad_l[4] = {0.0};
  double fUpwindQuad_r[4] = {0.0};
  double fUpwind_l[4] = {0.0};
  double fUpwind_r[4] = {0.0};
  double Ghat_l[4] = {0.0}; 
  double Ghat_r[4] = {0.0}; 

  if (alpha[0] > 0) { 
    fUpwindQuad_l[0] = ser_3x_p1_surfx2_quad_0_r(fl); 
    fUpwindQuad_r[0] = ser_3x_p1_surfx2_quad_0_r(fc); 
    fUpwindQuad_l[2] = ser_3x_p1_surfx2_quad_2_r(fl); 
    fUpwindQuad_r[2] = ser_3x_p1_surfx2_quad_2_r(fc); 
  } else { 
    fUpwindQuad_l[0] = ser_3x_p1_surfx2_quad_0_l(fc); 
    fUpwindQuad_r[0] = ser_3x_p1_surfx2_quad_0_l(fr); 
    fUpwindQuad_l[2] = ser_3x_p1_surfx2_quad_2_l(fc); 
    fUpwindQuad_r[2] = ser_3x_p1_surfx2_quad_2_l(fr); 
  } 
  if (alpha[0] > 0) { 
    fUpwindQuad_l[1] = ser_3x_p1_surfx2_quad_1_r(fl); 
    fUpwindQuad_r[1] = ser_3x_p1_surfx2_quad_1_r(fc); 
    fUpwindQuad_l[3] = ser_3x_p1_surfx2_quad_3_r(fl); 
    fUpwindQuad_r[3] = ser_3x_p1_surfx2_quad_3_r(fc); 
  } else { 
    fUpwindQuad_l[1] = ser_3x_p1_surfx2_quad_1_l(fc); 
    fUpwindQuad_r[1] = ser_3x_p1_surfx2_quad_1_l(fr); 
    fUpwindQuad_l[3] = ser_3x_p1_surfx2_quad_3_l(fc); 
    fUpwindQuad_r[3] = ser_3x_p1_surfx2_quad_3_l(fr); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_3x_p1_upwind(fUpwindQuad_l, fUpwind_l); 
  ser_3x_p1_upwind(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] += 0.5*alpha[0]*fUpwind_l[0]; 
  Ghat_l[1] += 0.5*alpha[0]*fUpwind_l[1]; 
  Ghat_l[2] += 0.5*alpha[0]*fUpwind_l[2]; 
  Ghat_l[3] += 0.5*alpha[0]*fUpwind_l[3]; 

  Ghat_r[0] += 0.5*alpha[0]*fUpwind_r[0]; 
  Ghat_r[1] += 0.5*alpha[0]*fUpwind_r[1]; 
  Ghat_r[2] += 0.5*alpha[0]*fUpwind_r[2]; 
  Ghat_r[3] += 0.5*alpha[0]*fUpwind_r[3]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv10; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv10; 
  out[2] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv10; 
  out[3] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv10; 
  out[4] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv10; 
  out[5] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv10; 
  out[6] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv10; 
  out[7] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv10; 

} 
