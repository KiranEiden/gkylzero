#include <gkyl_vlasov_kernels.h> 
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

  double alpha[3]; 

  alpha[0] = E0[0]; 
  alpha[1] = E0[1]; 
  alpha[2] = E0[2]; 

  double fUpwindQuad_l[3];
  double fUpwindQuad_r[3];
  if (0.6324555320336759*alpha[2]-0.9486832980505137*alpha[1]+0.7071067811865475*alpha[0] > 0) { 

    fUpwindQuad_l[0] = (-1.5*fl[7])+0.7745966692414833*fl[6]+1.118033988749895*fl[5]+0.4472135954999579*fl[4]-1.161895003862225*fl[3]+0.8660254037844386*fl[2]-0.6708203932499369*fl[1]+0.5*fl[0]; 
    fUpwindQuad_r[0] = (-1.5*fc[7])+0.7745966692414833*fc[6]+1.118033988749895*fc[5]+0.4472135954999579*fc[4]-1.161895003862225*fc[3]+0.8660254037844386*fc[2]-0.6708203932499369*fc[1]+0.5*fc[0]; 
  } else { 

    fUpwindQuad_l[0] = (-1.5*fc[7])-0.7745966692414833*fc[6]+1.118033988749895*fc[5]+0.4472135954999579*fc[4]+1.161895003862225*fc[3]-0.8660254037844386*fc[2]-0.6708203932499369*fc[1]+0.5*fc[0]; 
    fUpwindQuad_r[0] = (-1.5*fr[7])-0.7745966692414833*fr[6]+1.118033988749895*fr[5]+0.4472135954999579*fr[4]+1.161895003862225*fr[3]-0.8660254037844386*fr[2]-0.6708203932499369*fr[1]+0.5*fr[0]; 
  } 
  if (0.7071067811865475*alpha[0]-0.7905694150420947*alpha[2] > 0) { 

    fUpwindQuad_l[1] = (-0.9682458365518543*fl[6])+1.118033988749895*fl[5]-0.5590169943749475*fl[4]+0.8660254037844386*fl[2]+0.5*fl[0]; 
    fUpwindQuad_r[1] = (-0.9682458365518543*fc[6])+1.118033988749895*fc[5]-0.5590169943749475*fc[4]+0.8660254037844386*fc[2]+0.5*fc[0]; 
  } else { 

    fUpwindQuad_l[1] = 0.9682458365518543*fc[6]+1.118033988749895*fc[5]-0.5590169943749475*fc[4]-0.8660254037844386*fc[2]+0.5*fc[0]; 
    fUpwindQuad_r[1] = 0.9682458365518543*fr[6]+1.118033988749895*fr[5]-0.5590169943749475*fr[4]-0.8660254037844386*fr[2]+0.5*fr[0]; 
  } 
  if (0.6324555320336759*alpha[2]+0.9486832980505137*alpha[1]+0.7071067811865475*alpha[0] > 0) { 

    fUpwindQuad_l[2] = 1.5*fl[7]+0.7745966692414833*fl[6]+1.118033988749895*fl[5]+0.4472135954999579*fl[4]+1.161895003862225*fl[3]+0.8660254037844386*fl[2]+0.6708203932499369*fl[1]+0.5*fl[0]; 
    fUpwindQuad_r[2] = 1.5*fc[7]+0.7745966692414833*fc[6]+1.118033988749895*fc[5]+0.4472135954999579*fc[4]+1.161895003862225*fc[3]+0.8660254037844386*fc[2]+0.6708203932499369*fc[1]+0.5*fc[0]; 
  } else { 

    fUpwindQuad_l[2] = 1.5*fc[7]-0.7745966692414833*fc[6]+1.118033988749895*fc[5]+0.4472135954999579*fc[4]-1.161895003862225*fc[3]-0.8660254037844386*fc[2]+0.6708203932499369*fc[1]+0.5*fc[0]; 
    fUpwindQuad_r[2] = 1.5*fr[7]-0.7745966692414833*fr[6]+1.118033988749895*fr[5]+0.4472135954999579*fr[4]-1.161895003862225*fr[3]-0.8660254037844386*fr[2]+0.6708203932499369*fr[1]+0.5*fr[0]; 
  } 
  double fUpwind_l[3];
  fUpwind_l[0] = 0.07856742013183861*(5.0*fUpwindQuad_l[2]+8.0*fUpwindQuad_l[1]+5.0*fUpwindQuad_l[0]); 
  fUpwind_l[1] = 0.5270462766947298*(fUpwindQuad_l[2]-1.0*fUpwindQuad_l[0]); 
  fUpwind_l[2] = 0.3513641844631533*(fUpwindQuad_l[2]-2.0*fUpwindQuad_l[1]+fUpwindQuad_l[0]); 

  double fUpwind_r[3];
  fUpwind_r[0] = 0.07856742013183861*(5.0*fUpwindQuad_r[2]+8.0*fUpwindQuad_r[1]+5.0*fUpwindQuad_r[0]); 
  fUpwind_r[1] = 0.5270462766947298*(fUpwindQuad_r[2]-1.0*fUpwindQuad_r[0]); 
  fUpwind_r[2] = 0.3513641844631533*(fUpwindQuad_r[2]-2.0*fUpwindQuad_r[1]+fUpwindQuad_r[0]); 

  double Ghat_l[3]; 
  double Ghat_r[3]; 
  Ghat_l[0] = 0.7071067811865475*(alpha[2]*fUpwind_l[2]+alpha[1]*fUpwind_l[1]+alpha[0]*fUpwind_l[0]); 
  Ghat_l[1] = 0.1414213562373095*(4.47213595499958*(alpha[1]*fUpwind_l[2]+fUpwind_l[1]*alpha[2])+5.0*(alpha[0]*fUpwind_l[1]+fUpwind_l[0]*alpha[1])); 
  Ghat_l[2] = 0.02020305089104421*((22.3606797749979*alpha[2]+35.0*alpha[0])*fUpwind_l[2]+35.0*fUpwind_l[0]*alpha[2]+31.30495168499706*alpha[1]*fUpwind_l[1]); 

  Ghat_r[0] = 0.7071067811865475*(alpha[2]*fUpwind_r[2]+alpha[1]*fUpwind_r[1]+alpha[0]*fUpwind_r[0]); 
  Ghat_r[1] = 0.1414213562373095*(4.47213595499958*(alpha[1]*fUpwind_r[2]+fUpwind_r[1]*alpha[2])+5.0*(alpha[0]*fUpwind_r[1]+fUpwind_r[0]*alpha[1])); 
  Ghat_r[2] = 0.02020305089104421*((22.3606797749979*alpha[2]+35.0*alpha[0])*fUpwind_r[2]+35.0*fUpwind_r[0]*alpha[2]+31.30495168499706*alpha[1]*fUpwind_r[1]); 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv10; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv10; 
  out[2] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv10; 
  out[3] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv10; 
  out[4] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv10; 
  out[5] += (1.58113883008419*Ghat_l[0]-1.58113883008419*Ghat_r[0])*dv10; 
  out[6] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv10; 
  out[7] += (1.58113883008419*Ghat_l[1]-1.58113883008419*Ghat_r[1])*dv10; 

} 
