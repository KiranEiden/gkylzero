#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_4x_p1_surfx4_quad.h> 
#include <gkyl_basis_ser_4x_p1_upwind.h> 
GKYL_CU_DH void vlasov_surfvz_1x3v_ser_p1(const double *w, const double *dxv, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // qmem:      q/m*EM fields.
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 
  const double dv12 = 2/dxv[3]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double dv3 = dxv[3], wv3 = w[3]; 
  const double *E2 = &qmem[4]; 
  const double *B0 = &qmem[6]; 
  const double *B1 = &qmem[8]; 
  const double *B2 = &qmem[10]; 

  double alpha[8] = {0.0}; 

  alpha[0] = (-2.0*B0[0]*wv2)+2.0*B1[0]*wv1+2.0*E2[0]; 
  alpha[1] = (-2.0*B0[1]*wv2)+2.0*B1[1]*wv1+2.0*E2[1]; 
  alpha[2] = 0.5773502691896258*B1[0]*dv1; 
  alpha[3] = -0.5773502691896258*B0[0]*dv2; 
  alpha[4] = 0.5773502691896258*B1[1]*dv1; 
  alpha[5] = -0.5773502691896258*B0[1]*dv2; 

  double fUpwindQuad_l[8] = {0.0};
  double fUpwindQuad_r[8] = {0.0};
  double fUpwind_l[8] = {0.0};;
  double fUpwind_r[8] = {0.0};
  double Ghat_l[8] = {0.0}; 
  double Ghat_r[8] = {0.0}; 

  if (alpha[5]+alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad_l[0] = ser_4x_p1_surfx4_quad_0_r(fl); 
    fUpwindQuad_r[0] = ser_4x_p1_surfx4_quad_0_r(fc); 
  } else { 
    fUpwindQuad_l[0] = ser_4x_p1_surfx4_quad_0_l(fc); 
    fUpwindQuad_r[0] = ser_4x_p1_surfx4_quad_0_l(fr); 
  } 
  if ((-alpha[5])-alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad_l[1] = ser_4x_p1_surfx4_quad_1_r(fl); 
    fUpwindQuad_r[1] = ser_4x_p1_surfx4_quad_1_r(fc); 
  } else { 
    fUpwindQuad_l[1] = ser_4x_p1_surfx4_quad_1_l(fc); 
    fUpwindQuad_r[1] = ser_4x_p1_surfx4_quad_1_l(fr); 
  } 
  if (alpha[5]-alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad_l[2] = ser_4x_p1_surfx4_quad_2_r(fl); 
    fUpwindQuad_r[2] = ser_4x_p1_surfx4_quad_2_r(fc); 
  } else { 
    fUpwindQuad_l[2] = ser_4x_p1_surfx4_quad_2_l(fc); 
    fUpwindQuad_r[2] = ser_4x_p1_surfx4_quad_2_l(fr); 
  } 
  if ((-alpha[5])+alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad_l[3] = ser_4x_p1_surfx4_quad_3_r(fl); 
    fUpwindQuad_r[3] = ser_4x_p1_surfx4_quad_3_r(fc); 
  } else { 
    fUpwindQuad_l[3] = ser_4x_p1_surfx4_quad_3_l(fc); 
    fUpwindQuad_r[3] = ser_4x_p1_surfx4_quad_3_l(fr); 
  } 
  if ((-alpha[5])+alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad_l[4] = ser_4x_p1_surfx4_quad_4_r(fl); 
    fUpwindQuad_r[4] = ser_4x_p1_surfx4_quad_4_r(fc); 
  } else { 
    fUpwindQuad_l[4] = ser_4x_p1_surfx4_quad_4_l(fc); 
    fUpwindQuad_r[4] = ser_4x_p1_surfx4_quad_4_l(fr); 
  } 
  if (alpha[5]-alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad_l[5] = ser_4x_p1_surfx4_quad_5_r(fl); 
    fUpwindQuad_r[5] = ser_4x_p1_surfx4_quad_5_r(fc); 
  } else { 
    fUpwindQuad_l[5] = ser_4x_p1_surfx4_quad_5_l(fc); 
    fUpwindQuad_r[5] = ser_4x_p1_surfx4_quad_5_l(fr); 
  } 
  if ((-alpha[5])-alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad_l[6] = ser_4x_p1_surfx4_quad_6_r(fl); 
    fUpwindQuad_r[6] = ser_4x_p1_surfx4_quad_6_r(fc); 
  } else { 
    fUpwindQuad_l[6] = ser_4x_p1_surfx4_quad_6_l(fc); 
    fUpwindQuad_r[6] = ser_4x_p1_surfx4_quad_6_l(fr); 
  } 
  if (alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad_l[7] = ser_4x_p1_surfx4_quad_7_r(fl); 
    fUpwindQuad_r[7] = ser_4x_p1_surfx4_quad_7_r(fc); 
  } else { 
    fUpwindQuad_l[7] = ser_4x_p1_surfx4_quad_7_l(fc); 
    fUpwindQuad_r[7] = ser_4x_p1_surfx4_quad_7_l(fr); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_4x_p1_upwind(fUpwindQuad_l, fUpwind_l); 
  ser_4x_p1_upwind(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] += 0.3535533905932737*(alpha[5]*fUpwind_l[5]+alpha[4]*fUpwind_l[4]+alpha[3]*fUpwind_l[3]+alpha[2]*fUpwind_l[2]+alpha[1]*fUpwind_l[1]+alpha[0]*fUpwind_l[0]); 
  Ghat_l[1] += 0.3535533905932737*(alpha[3]*fUpwind_l[5]+fUpwind_l[3]*alpha[5]+alpha[2]*fUpwind_l[4]+fUpwind_l[2]*alpha[4]+alpha[0]*fUpwind_l[1]+fUpwind_l[0]*alpha[1]); 
  Ghat_l[2] += 0.3535533905932737*(alpha[5]*fUpwind_l[7]+alpha[3]*fUpwind_l[6]+alpha[1]*fUpwind_l[4]+fUpwind_l[1]*alpha[4]+alpha[0]*fUpwind_l[2]+fUpwind_l[0]*alpha[2]); 
  Ghat_l[3] += 0.3535533905932737*(alpha[4]*fUpwind_l[7]+alpha[2]*fUpwind_l[6]+alpha[1]*fUpwind_l[5]+fUpwind_l[1]*alpha[5]+alpha[0]*fUpwind_l[3]+fUpwind_l[0]*alpha[3]); 
  Ghat_l[4] += 0.3535533905932737*(alpha[3]*fUpwind_l[7]+alpha[5]*fUpwind_l[6]+alpha[0]*fUpwind_l[4]+fUpwind_l[0]*alpha[4]+alpha[1]*fUpwind_l[2]+fUpwind_l[1]*alpha[2]); 
  Ghat_l[5] += 0.3535533905932737*(alpha[2]*fUpwind_l[7]+alpha[4]*fUpwind_l[6]+alpha[0]*fUpwind_l[5]+fUpwind_l[0]*alpha[5]+alpha[1]*fUpwind_l[3]+fUpwind_l[1]*alpha[3]); 
  Ghat_l[6] += 0.3535533905932737*(alpha[1]*fUpwind_l[7]+alpha[0]*fUpwind_l[6]+alpha[4]*fUpwind_l[5]+fUpwind_l[4]*alpha[5]+alpha[2]*fUpwind_l[3]+fUpwind_l[2]*alpha[3]); 
  Ghat_l[7] += 0.3535533905932737*(alpha[0]*fUpwind_l[7]+alpha[1]*fUpwind_l[6]+alpha[2]*fUpwind_l[5]+fUpwind_l[2]*alpha[5]+alpha[3]*fUpwind_l[4]+fUpwind_l[3]*alpha[4]); 

  Ghat_r[0] += 0.3535533905932737*(alpha[5]*fUpwind_r[5]+alpha[4]*fUpwind_r[4]+alpha[3]*fUpwind_r[3]+alpha[2]*fUpwind_r[2]+alpha[1]*fUpwind_r[1]+alpha[0]*fUpwind_r[0]); 
  Ghat_r[1] += 0.3535533905932737*(alpha[3]*fUpwind_r[5]+fUpwind_r[3]*alpha[5]+alpha[2]*fUpwind_r[4]+fUpwind_r[2]*alpha[4]+alpha[0]*fUpwind_r[1]+fUpwind_r[0]*alpha[1]); 
  Ghat_r[2] += 0.3535533905932737*(alpha[5]*fUpwind_r[7]+alpha[3]*fUpwind_r[6]+alpha[1]*fUpwind_r[4]+fUpwind_r[1]*alpha[4]+alpha[0]*fUpwind_r[2]+fUpwind_r[0]*alpha[2]); 
  Ghat_r[3] += 0.3535533905932737*(alpha[4]*fUpwind_r[7]+alpha[2]*fUpwind_r[6]+alpha[1]*fUpwind_r[5]+fUpwind_r[1]*alpha[5]+alpha[0]*fUpwind_r[3]+fUpwind_r[0]*alpha[3]); 
  Ghat_r[4] += 0.3535533905932737*(alpha[3]*fUpwind_r[7]+alpha[5]*fUpwind_r[6]+alpha[0]*fUpwind_r[4]+fUpwind_r[0]*alpha[4]+alpha[1]*fUpwind_r[2]+fUpwind_r[1]*alpha[2]); 
  Ghat_r[5] += 0.3535533905932737*(alpha[2]*fUpwind_r[7]+alpha[4]*fUpwind_r[6]+alpha[0]*fUpwind_r[5]+fUpwind_r[0]*alpha[5]+alpha[1]*fUpwind_r[3]+fUpwind_r[1]*alpha[3]); 
  Ghat_r[6] += 0.3535533905932737*(alpha[1]*fUpwind_r[7]+alpha[0]*fUpwind_r[6]+alpha[4]*fUpwind_r[5]+fUpwind_r[4]*alpha[5]+alpha[2]*fUpwind_r[3]+fUpwind_r[2]*alpha[3]); 
  Ghat_r[7] += 0.3535533905932737*(alpha[0]*fUpwind_r[7]+alpha[1]*fUpwind_r[6]+alpha[2]*fUpwind_r[5]+fUpwind_r[2]*alpha[5]+alpha[3]*fUpwind_r[4]+fUpwind_r[3]*alpha[4]); 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv12; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv12; 
  out[2] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv12; 
  out[3] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv12; 
  out[4] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv12; 
  out[5] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dv12; 
  out[6] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dv12; 
  out[7] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dv12; 
  out[8] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv12; 
  out[9] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv12; 
  out[10] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv12; 
  out[11] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dv12; 
  out[12] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dv12; 
  out[13] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dv12; 
  out[14] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dv12; 
  out[15] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dv12; 

} 
