#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_1x1v_p2_surfvx_quad.h> 
GKYL_CU_DH void vlasov_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // qmem:      q/m*EM fields.
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 
  const double dv10 = 2/dxv[1]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double *E0 = &qmem[0]; 

  double alpha[3] = {0.0}; 

  alpha[0] = E0[0]; 
  alpha[1] = E0[1]; 
  alpha[2] = E0[2]; 

  double fUpwindQuad_l[3] = {0.0};
  double fUpwindQuad_r[3] = {0.0};
  double fUpwind_l[3] = {0.0};;
  double fUpwind_r[3] = {0.0};
  double Ghat_l[3] = {0.0}; 
  double Ghat_r[3] = {0.0}; 

  if (0.6324555320336759*alpha[2]-0.9486832980505137*alpha[1]+0.7071067811865475*alpha[0] > 0) { 

    fUpwindQuad_l[0] = ser_1x1v_p2_surfvx_quad_0(1, fl); 
    fUpwindQuad_r[0] = ser_1x1v_p2_surfvx_quad_0(1, fc); 
  } else { 

    fUpwindQuad_l[0] = ser_1x1v_p2_surfvx_quad_0(-1, fc); 
    fUpwindQuad_r[0] = ser_1x1v_p2_surfvx_quad_0(-1, fr); 
  } 
  if (0.7071067811865475*alpha[0]-0.7905694150420947*alpha[2] > 0) { 

    fUpwindQuad_l[1] = ser_1x1v_p2_surfvx_quad_1(1, fl); 
    fUpwindQuad_r[1] = ser_1x1v_p2_surfvx_quad_1(1, fc); 
  } else { 

    fUpwindQuad_l[1] = ser_1x1v_p2_surfvx_quad_1(-1, fc); 
    fUpwindQuad_r[1] = ser_1x1v_p2_surfvx_quad_1(-1, fr); 
  } 
  if (0.6324555320336759*alpha[2]+0.9486832980505137*alpha[1]+0.7071067811865475*alpha[0] > 0) { 

    fUpwindQuad_l[2] = ser_1x1v_p2_surfvx_quad_2(1, fl); 
    fUpwindQuad_r[2] = ser_1x1v_p2_surfvx_quad_2(1, fc); 
  } else { 

    fUpwindQuad_l[2] = ser_1x1v_p2_surfvx_quad_2(-1, fc); 
    fUpwindQuad_r[2] = ser_1x1v_p2_surfvx_quad_2(-1, fr); 
  } 
  fUpwind_l[0] = 0.392837100659193*fUpwindQuad_l[2]+0.6285393610547091*fUpwindQuad_l[1]+0.392837100659193*fUpwindQuad_l[0]; 
  fUpwind_l[1] = 0.5270462766947298*fUpwindQuad_l[2]-0.5270462766947298*fUpwindQuad_l[0]; 
  fUpwind_l[2] = 0.3513641844631533*fUpwindQuad_l[2]-0.7027283689263066*fUpwindQuad_l[1]+0.3513641844631533*fUpwindQuad_l[0]; 
  fUpwind_r[0] = 0.392837100659193*fUpwindQuad_r[2]+0.6285393610547091*fUpwindQuad_r[1]+0.392837100659193*fUpwindQuad_r[0]; 
  fUpwind_r[1] = 0.5270462766947298*fUpwindQuad_r[2]-0.5270462766947298*fUpwindQuad_r[0]; 
  fUpwind_r[2] = 0.3513641844631533*fUpwindQuad_r[2]-0.7027283689263066*fUpwindQuad_r[1]+0.3513641844631533*fUpwindQuad_r[0]; 

  Ghat_l[0] += 0.7071067811865475*(alpha[2]*fUpwind_l[2]+alpha[1]*fUpwind_l[1]+alpha[0]*fUpwind_l[0]); 
  Ghat_l[1] += 0.6324555320336759*(alpha[1]*fUpwind_l[2]+fUpwind_l[1]*alpha[2])+0.7071067811865475*(alpha[0]*fUpwind_l[1]+fUpwind_l[0]*alpha[1]); 
  Ghat_l[2] += 0.4517539514526256*alpha[2]*fUpwind_l[2]+0.7071067811865475*(alpha[0]*fUpwind_l[2]+fUpwind_l[0]*alpha[2])+0.6324555320336759*alpha[1]*fUpwind_l[1]; 

  Ghat_r[0] += 0.7071067811865475*(alpha[2]*fUpwind_r[2]+alpha[1]*fUpwind_r[1]+alpha[0]*fUpwind_r[0]); 
  Ghat_r[1] += 0.6324555320336759*(alpha[1]*fUpwind_r[2]+fUpwind_r[1]*alpha[2])+0.7071067811865475*(alpha[0]*fUpwind_r[1]+fUpwind_r[0]*alpha[1]); 
  Ghat_r[2] += 0.4517539514526256*alpha[2]*fUpwind_r[2]+0.7071067811865475*(alpha[0]*fUpwind_r[2]+fUpwind_r[0]*alpha[2])+0.6324555320336759*alpha[1]*fUpwind_r[1]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv10; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv10; 
  out[2] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv10; 
  out[3] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv10; 
  out[4] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv10; 
  out[5] += (1.58113883008419*Ghat_l[0]-1.58113883008419*Ghat_r[0])*dv10; 
  out[6] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv10; 
  out[7] += (1.58113883008419*Ghat_l[1]-1.58113883008419*Ghat_r[1])*dv10; 

} 
