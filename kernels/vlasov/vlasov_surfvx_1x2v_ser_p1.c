#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_1x2v_p1_surfvx_quad.h> 
GKYL_CU_DH void vlasov_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // qmem:      q/m*EM fields.
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 
  const double dv10 = 2/dxv[1]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double *E0 = &qmem[0]; 
  const double *B2 = &qmem[10]; 

  double alpha[4] = {0.0}; 

  alpha[0] = 1.414213562373095*B2[0]*wv2+1.414213562373095*E0[0]; 
  alpha[1] = 1.414213562373095*B2[1]*wv2+1.414213562373095*E0[1]; 
  alpha[2] = 0.408248290463863*B2[0]*dv2; 
  alpha[3] = 0.408248290463863*B2[1]*dv2; 

  double fUpwindQuad_l[4] = {0.0};
  double fUpwindQuad_r[4] = {0.0};
  double fUpwind_l[4] = {0.0};;
  double fUpwind_r[4] = {0.0};
  double Ghat_l[4] = {0.0}; 
  double Ghat_r[4] = {0.0}; 

  if (alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[0] = ser_1x2v_p1_surfvx_quad_0(1, fl); 
    fUpwindQuad_r[0] = ser_1x2v_p1_surfvx_quad_0(1, fc); 
  } else { 

    fUpwindQuad_l[0] = ser_1x2v_p1_surfvx_quad_0(-1, fc); 
    fUpwindQuad_r[0] = ser_1x2v_p1_surfvx_quad_0(-1, fr); 
  } 
  if ((-alpha[3])-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[1] = ser_1x2v_p1_surfvx_quad_1(1, fl); 
    fUpwindQuad_r[1] = ser_1x2v_p1_surfvx_quad_1(1, fc); 
  } else { 

    fUpwindQuad_l[1] = ser_1x2v_p1_surfvx_quad_1(-1, fc); 
    fUpwindQuad_r[1] = ser_1x2v_p1_surfvx_quad_1(-1, fr); 
  } 
  if ((-alpha[3])+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[2] = ser_1x2v_p1_surfvx_quad_2(1, fl); 
    fUpwindQuad_r[2] = ser_1x2v_p1_surfvx_quad_2(1, fc); 
  } else { 

    fUpwindQuad_l[2] = ser_1x2v_p1_surfvx_quad_2(-1, fc); 
    fUpwindQuad_r[2] = ser_1x2v_p1_surfvx_quad_2(-1, fr); 
  } 
  if (alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[3] = ser_1x2v_p1_surfvx_quad_3(1, fl); 
    fUpwindQuad_r[3] = ser_1x2v_p1_surfvx_quad_3(1, fc); 
  } else { 

    fUpwindQuad_l[3] = ser_1x2v_p1_surfvx_quad_3(-1, fc); 
    fUpwindQuad_r[3] = ser_1x2v_p1_surfvx_quad_3(-1, fr); 
  } 
  fUpwind_l[0] = 0.5*(fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[1] = 0.5*(fUpwindQuad_l[3]-1.0*fUpwindQuad_l[2]+fUpwindQuad_l[1]-1.0*fUpwindQuad_l[0]); 
  fUpwind_l[2] = 0.5*(fUpwindQuad_l[3]+fUpwindQuad_l[2]-1.0*(fUpwindQuad_l[1]+fUpwindQuad_l[0])); 
  fUpwind_l[3] = 0.5*(fUpwindQuad_l[3]-1.0*(fUpwindQuad_l[2]+fUpwindQuad_l[1])+fUpwindQuad_l[0]); 

  fUpwind_r[0] = 0.5*(fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[1] = 0.5*(fUpwindQuad_r[3]-1.0*fUpwindQuad_r[2]+fUpwindQuad_r[1]-1.0*fUpwindQuad_r[0]); 
  fUpwind_r[2] = 0.5*(fUpwindQuad_r[3]+fUpwindQuad_r[2]-1.0*(fUpwindQuad_r[1]+fUpwindQuad_r[0])); 
  fUpwind_r[3] = 0.5*(fUpwindQuad_r[3]-1.0*(fUpwindQuad_r[2]+fUpwindQuad_r[1])+fUpwindQuad_r[0]); 

  Ghat_l[0] += 0.5*(alpha[3]*fUpwind_l[3]+alpha[2]*fUpwind_l[2]+alpha[1]*fUpwind_l[1]+alpha[0]*fUpwind_l[0]); 
  Ghat_l[1] += 0.5*(alpha[2]*fUpwind_l[3]+fUpwind_l[2]*alpha[3]+alpha[0]*fUpwind_l[1]+fUpwind_l[0]*alpha[1]); 
  Ghat_l[2] += 0.5*(alpha[1]*fUpwind_l[3]+fUpwind_l[1]*alpha[3]+alpha[0]*fUpwind_l[2]+fUpwind_l[0]*alpha[2]); 
  Ghat_l[3] += 0.5*(alpha[0]*fUpwind_l[3]+fUpwind_l[0]*alpha[3]+alpha[1]*fUpwind_l[2]+fUpwind_l[1]*alpha[2]); 

  Ghat_r[0] += 0.5*(alpha[3]*fUpwind_r[3]+alpha[2]*fUpwind_r[2]+alpha[1]*fUpwind_r[1]+alpha[0]*fUpwind_r[0]); 
  Ghat_r[1] += 0.5*(alpha[2]*fUpwind_r[3]+fUpwind_r[2]*alpha[3]+alpha[0]*fUpwind_r[1]+fUpwind_r[0]*alpha[1]); 
  Ghat_r[2] += 0.5*(alpha[1]*fUpwind_r[3]+fUpwind_r[1]*alpha[3]+alpha[0]*fUpwind_r[2]+fUpwind_r[0]*alpha[2]); 
  Ghat_r[3] += 0.5*(alpha[0]*fUpwind_r[3]+fUpwind_r[0]*alpha[3]+alpha[1]*fUpwind_r[2]+fUpwind_r[1]*alpha[2]); 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv10; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv10; 
  out[2] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv10; 
  out[3] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv10; 
  out[4] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv10; 
  out[5] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv10; 
  out[6] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv10; 
  out[7] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv10; 

} 
