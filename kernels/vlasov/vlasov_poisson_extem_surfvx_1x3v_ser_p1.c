#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_1x3v_p1_surfvx_quad.h> 
GKYL_CU_DH void vlasov_poisson_extem_surfvx_1x3v_ser_p1(const double *w, const double *dxv, const double *fac_phi, const double *vecA, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
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
  const double dv3 = dxv[3], wv3 = w[3]; 
  const double *phi = &fac_phi[0]; 
  const double dx10 = 2/dxv[0]; 
  const double *A0 = &vecA[0]; 
  const double *A1 = &vecA[2]; 
  const double *A2 = &vecA[4]; 
  double alpha[8] = {0.0}; 

  alpha[0] = 3.464101615137754*A2[1]*dx10*wv3+3.464101615137754*A1[1]*dx10*wv2-3.464101615137754*phi[1]*dx10; 
  alpha[2] = A1[1]*dv2*dx10; 
  alpha[3] = A2[1]*dv3*dx10; 

  double fUpwindQuad_l[8] = {0.0};
  double fUpwindQuad_r[8] = {0.0};
  double fUpwind_l[8] = {0.0};;
  double fUpwind_r[8] = {0.0};
  double Ghat_l[8] = {0.0}; 
  double Ghat_r[8] = {0.0}; 

  if ((-alpha[3])-alpha[2]+alpha[0] > 0) { 

    fUpwindQuad_l[0] = ser_1x3v_p1_surfvx_quad_0(1, fl); 
    fUpwindQuad_r[0] = ser_1x3v_p1_surfvx_quad_0(1, fc); 
  } else { 

    fUpwindQuad_l[0] = ser_1x3v_p1_surfvx_quad_0(-1, fc); 
    fUpwindQuad_r[0] = ser_1x3v_p1_surfvx_quad_0(-1, fr); 
  } 
  if ((-alpha[3])-alpha[2]+alpha[0] > 0) { 

    fUpwindQuad_l[1] = ser_1x3v_p1_surfvx_quad_1(1, fl); 
    fUpwindQuad_r[1] = ser_1x3v_p1_surfvx_quad_1(1, fc); 
  } else { 

    fUpwindQuad_l[1] = ser_1x3v_p1_surfvx_quad_1(-1, fc); 
    fUpwindQuad_r[1] = ser_1x3v_p1_surfvx_quad_1(-1, fr); 
  } 
  if ((-alpha[3])+alpha[2]+alpha[0] > 0) { 

    fUpwindQuad_l[2] = ser_1x3v_p1_surfvx_quad_2(1, fl); 
    fUpwindQuad_r[2] = ser_1x3v_p1_surfvx_quad_2(1, fc); 
  } else { 

    fUpwindQuad_l[2] = ser_1x3v_p1_surfvx_quad_2(-1, fc); 
    fUpwindQuad_r[2] = ser_1x3v_p1_surfvx_quad_2(-1, fr); 
  } 
  if ((-alpha[3])+alpha[2]+alpha[0] > 0) { 

    fUpwindQuad_l[3] = ser_1x3v_p1_surfvx_quad_3(1, fl); 
    fUpwindQuad_r[3] = ser_1x3v_p1_surfvx_quad_3(1, fc); 
  } else { 

    fUpwindQuad_l[3] = ser_1x3v_p1_surfvx_quad_3(-1, fc); 
    fUpwindQuad_r[3] = ser_1x3v_p1_surfvx_quad_3(-1, fr); 
  } 
  if (alpha[3]-alpha[2]+alpha[0] > 0) { 

    fUpwindQuad_l[4] = ser_1x3v_p1_surfvx_quad_4(1, fl); 
    fUpwindQuad_r[4] = ser_1x3v_p1_surfvx_quad_4(1, fc); 
  } else { 

    fUpwindQuad_l[4] = ser_1x3v_p1_surfvx_quad_4(-1, fc); 
    fUpwindQuad_r[4] = ser_1x3v_p1_surfvx_quad_4(-1, fr); 
  } 
  if (alpha[3]-alpha[2]+alpha[0] > 0) { 

    fUpwindQuad_l[5] = ser_1x3v_p1_surfvx_quad_5(1, fl); 
    fUpwindQuad_r[5] = ser_1x3v_p1_surfvx_quad_5(1, fc); 
  } else { 

    fUpwindQuad_l[5] = ser_1x3v_p1_surfvx_quad_5(-1, fc); 
    fUpwindQuad_r[5] = ser_1x3v_p1_surfvx_quad_5(-1, fr); 
  } 
  if (alpha[3]+alpha[2]+alpha[0] > 0) { 

    fUpwindQuad_l[6] = ser_1x3v_p1_surfvx_quad_6(1, fl); 
    fUpwindQuad_r[6] = ser_1x3v_p1_surfvx_quad_6(1, fc); 
  } else { 

    fUpwindQuad_l[6] = ser_1x3v_p1_surfvx_quad_6(-1, fc); 
    fUpwindQuad_r[6] = ser_1x3v_p1_surfvx_quad_6(-1, fr); 
  } 
  if (alpha[3]+alpha[2]+alpha[0] > 0) { 

    fUpwindQuad_l[7] = ser_1x3v_p1_surfvx_quad_7(1, fl); 
    fUpwindQuad_r[7] = ser_1x3v_p1_surfvx_quad_7(1, fc); 
  } else { 

    fUpwindQuad_l[7] = ser_1x3v_p1_surfvx_quad_7(-1, fc); 
    fUpwindQuad_r[7] = ser_1x3v_p1_surfvx_quad_7(-1, fr); 
  } 
  fUpwind_l[0] = 0.3535533905932737*(fUpwindQuad_l[7]+fUpwindQuad_l[6]+fUpwindQuad_l[5]+fUpwindQuad_l[4]+fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[1] = 0.3535533905932737*(fUpwindQuad_l[7]-1.0*fUpwindQuad_l[6]+fUpwindQuad_l[5]-1.0*fUpwindQuad_l[4]+fUpwindQuad_l[3]-1.0*fUpwindQuad_l[2]+fUpwindQuad_l[1]-1.0*fUpwindQuad_l[0]); 
  fUpwind_l[2] = 0.3535533905932737*(fUpwindQuad_l[7]+fUpwindQuad_l[6]-1.0*(fUpwindQuad_l[5]+fUpwindQuad_l[4])+fUpwindQuad_l[3]+fUpwindQuad_l[2]-1.0*(fUpwindQuad_l[1]+fUpwindQuad_l[0])); 
  fUpwind_l[3] = 0.3535533905932737*(fUpwindQuad_l[7]+fUpwindQuad_l[6]+fUpwindQuad_l[5]+fUpwindQuad_l[4]-1.0*(fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]+fUpwindQuad_l[0])); 
  fUpwind_l[4] = 0.3535533905932737*(fUpwindQuad_l[7]-1.0*(fUpwindQuad_l[6]+fUpwindQuad_l[5])+fUpwindQuad_l[4]+fUpwindQuad_l[3]-1.0*(fUpwindQuad_l[2]+fUpwindQuad_l[1])+fUpwindQuad_l[0]); 
  fUpwind_l[5] = 0.3535533905932737*(fUpwindQuad_l[7]-1.0*fUpwindQuad_l[6]+fUpwindQuad_l[5]-1.0*(fUpwindQuad_l[4]+fUpwindQuad_l[3])+fUpwindQuad_l[2]-1.0*fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[6] = 0.3535533905932737*(fUpwindQuad_l[7]+fUpwindQuad_l[6]-1.0*(fUpwindQuad_l[5]+fUpwindQuad_l[4]+fUpwindQuad_l[3]+fUpwindQuad_l[2])+fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[7] = 0.3535533905932737*(fUpwindQuad_l[7]-1.0*(fUpwindQuad_l[6]+fUpwindQuad_l[5])+fUpwindQuad_l[4]-1.0*fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]-1.0*fUpwindQuad_l[0]); 

  fUpwind_r[0] = 0.3535533905932737*(fUpwindQuad_r[7]+fUpwindQuad_r[6]+fUpwindQuad_r[5]+fUpwindQuad_r[4]+fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[1] = 0.3535533905932737*(fUpwindQuad_r[7]-1.0*fUpwindQuad_r[6]+fUpwindQuad_r[5]-1.0*fUpwindQuad_r[4]+fUpwindQuad_r[3]-1.0*fUpwindQuad_r[2]+fUpwindQuad_r[1]-1.0*fUpwindQuad_r[0]); 
  fUpwind_r[2] = 0.3535533905932737*(fUpwindQuad_r[7]+fUpwindQuad_r[6]-1.0*(fUpwindQuad_r[5]+fUpwindQuad_r[4])+fUpwindQuad_r[3]+fUpwindQuad_r[2]-1.0*(fUpwindQuad_r[1]+fUpwindQuad_r[0])); 
  fUpwind_r[3] = 0.3535533905932737*(fUpwindQuad_r[7]+fUpwindQuad_r[6]+fUpwindQuad_r[5]+fUpwindQuad_r[4]-1.0*(fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]+fUpwindQuad_r[0])); 
  fUpwind_r[4] = 0.3535533905932737*(fUpwindQuad_r[7]-1.0*(fUpwindQuad_r[6]+fUpwindQuad_r[5])+fUpwindQuad_r[4]+fUpwindQuad_r[3]-1.0*(fUpwindQuad_r[2]+fUpwindQuad_r[1])+fUpwindQuad_r[0]); 
  fUpwind_r[5] = 0.3535533905932737*(fUpwindQuad_r[7]-1.0*fUpwindQuad_r[6]+fUpwindQuad_r[5]-1.0*(fUpwindQuad_r[4]+fUpwindQuad_r[3])+fUpwindQuad_r[2]-1.0*fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[6] = 0.3535533905932737*(fUpwindQuad_r[7]+fUpwindQuad_r[6]-1.0*(fUpwindQuad_r[5]+fUpwindQuad_r[4]+fUpwindQuad_r[3]+fUpwindQuad_r[2])+fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[7] = 0.3535533905932737*(fUpwindQuad_r[7]-1.0*(fUpwindQuad_r[6]+fUpwindQuad_r[5])+fUpwindQuad_r[4]-1.0*fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]-1.0*fUpwindQuad_r[0]); 

  Ghat_l[0] += 0.3535533905932737*(alpha[3]*fUpwind_l[3]+alpha[2]*fUpwind_l[2]+alpha[0]*fUpwind_l[0]); 
  Ghat_l[1] += 0.3535533905932737*(alpha[3]*fUpwind_l[5]+alpha[2]*fUpwind_l[4]+alpha[0]*fUpwind_l[1]); 
  Ghat_l[2] += 0.3535533905932737*(alpha[3]*fUpwind_l[6]+alpha[0]*fUpwind_l[2]+fUpwind_l[0]*alpha[2]); 
  Ghat_l[3] += 0.3535533905932737*(alpha[2]*fUpwind_l[6]+alpha[0]*fUpwind_l[3]+fUpwind_l[0]*alpha[3]); 
  Ghat_l[4] += 0.3535533905932737*(alpha[3]*fUpwind_l[7]+alpha[0]*fUpwind_l[4]+fUpwind_l[1]*alpha[2]); 
  Ghat_l[5] += 0.3535533905932737*(alpha[2]*fUpwind_l[7]+alpha[0]*fUpwind_l[5]+fUpwind_l[1]*alpha[3]); 
  Ghat_l[6] += 0.3535533905932737*(alpha[0]*fUpwind_l[6]+alpha[2]*fUpwind_l[3]+fUpwind_l[2]*alpha[3]); 
  Ghat_l[7] += 0.3535533905932737*(alpha[0]*fUpwind_l[7]+alpha[2]*fUpwind_l[5]+alpha[3]*fUpwind_l[4]); 

  Ghat_r[0] += 0.3535533905932737*(alpha[3]*fUpwind_r[3]+alpha[2]*fUpwind_r[2]+alpha[0]*fUpwind_r[0]); 
  Ghat_r[1] += 0.3535533905932737*(alpha[3]*fUpwind_r[5]+alpha[2]*fUpwind_r[4]+alpha[0]*fUpwind_r[1]); 
  Ghat_r[2] += 0.3535533905932737*(alpha[3]*fUpwind_r[6]+alpha[0]*fUpwind_r[2]+fUpwind_r[0]*alpha[2]); 
  Ghat_r[3] += 0.3535533905932737*(alpha[2]*fUpwind_r[6]+alpha[0]*fUpwind_r[3]+fUpwind_r[0]*alpha[3]); 
  Ghat_r[4] += 0.3535533905932737*(alpha[3]*fUpwind_r[7]+alpha[0]*fUpwind_r[4]+fUpwind_r[1]*alpha[2]); 
  Ghat_r[5] += 0.3535533905932737*(alpha[2]*fUpwind_r[7]+alpha[0]*fUpwind_r[5]+fUpwind_r[1]*alpha[3]); 
  Ghat_r[6] += 0.3535533905932737*(alpha[0]*fUpwind_r[6]+alpha[2]*fUpwind_r[3]+fUpwind_r[2]*alpha[3]); 
  Ghat_r[7] += 0.3535533905932737*(alpha[0]*fUpwind_r[7]+alpha[2]*fUpwind_r[5]+alpha[3]*fUpwind_r[4]); 

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

} 