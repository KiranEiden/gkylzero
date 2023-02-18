#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH void lbo_vlasov_pkpm_drag_boundary_surfvpar_1x1v_ser_p2(const double *w, const double *dxv, const double *nu, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[2]:       Cell-center coordinates. 
  // dxv[2]:     Cell spacing. 
  // nu:          Collisionality. 
  // fSkin/fEdge: Input Distribution function [F_0, T_perp G = T_perp (F_1 - F_0)] in skin cell/last edge cell 
  // out:         Incremented distribution function in cell 
  const double dv1par = 2.0/dxv[1]; 
  const double dvpar = dxv[1], wvpar = w[1]; 
  const double *F_0Skin = &fSkin[0]; 
  const double *G_1Skin = &fSkin[8]; 
  const double *F_0Edge = &fEdge[0]; 
  const double *G_1Edge = &fEdge[8]; 
  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[8]; 

  double alphaDrSurf[3] = {0.0}; 
  double Ghat_F_0[3] = {0.0}; 
  double Ghat_G_1[3] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = nu[0]*wvpar+0.5*nu[0]*dvpar; 
  alphaDrSurf[1] = nu[1]*wvpar+0.5*nu[1]*dvpar; 
  alphaDrSurf[2] = nu[2]*wvpar+0.5*nu[2]*dvpar; 

  Ghat_F_0[0] = 1.118033988749895*alphaDrSurf[1]*F_0Skin[7]+0.8660254037844387*alphaDrSurf[2]*F_0Skin[6]+1.118033988749895*alphaDrSurf[0]*F_0Skin[5]+0.5*alphaDrSurf[2]*F_0Skin[4]+0.8660254037844386*alphaDrSurf[1]*F_0Skin[3]+0.8660254037844386*alphaDrSurf[0]*F_0Skin[2]+0.5*F_0Skin[1]*alphaDrSurf[1]+0.5*F_0Skin[0]*alphaDrSurf[0]; 
  Ghat_F_0[1] = 1.0*alphaDrSurf[2]*F_0Skin[7]+1.118033988749895*alphaDrSurf[0]*F_0Skin[7]+0.7745966692414834*alphaDrSurf[1]*F_0Skin[6]+1.118033988749895*alphaDrSurf[1]*F_0Skin[5]+0.4472135954999579*alphaDrSurf[1]*F_0Skin[4]+0.7745966692414833*alphaDrSurf[2]*F_0Skin[3]+0.8660254037844386*alphaDrSurf[0]*F_0Skin[3]+0.4472135954999579*F_0Skin[1]*alphaDrSurf[2]+0.8660254037844386*alphaDrSurf[1]*F_0Skin[2]+0.5*F_0Skin[0]*alphaDrSurf[1]+0.5*alphaDrSurf[0]*F_0Skin[1]; 
  Ghat_F_0[2] = 1.0*alphaDrSurf[1]*F_0Skin[7]+0.5532833351724881*alphaDrSurf[2]*F_0Skin[6]+0.8660254037844387*alphaDrSurf[0]*F_0Skin[6]+1.118033988749895*alphaDrSurf[2]*F_0Skin[5]+0.31943828249997*alphaDrSurf[2]*F_0Skin[4]+0.5*alphaDrSurf[0]*F_0Skin[4]+0.7745966692414833*alphaDrSurf[1]*F_0Skin[3]+0.8660254037844386*F_0Skin[2]*alphaDrSurf[2]+0.5*F_0Skin[0]*alphaDrSurf[2]+0.4472135954999579*F_0Skin[1]*alphaDrSurf[1]; 
  Ghat_G_1[0] = 1.118033988749895*alphaDrSurf[1]*G_1Skin[7]+0.8660254037844387*alphaDrSurf[2]*G_1Skin[6]+1.118033988749895*alphaDrSurf[0]*G_1Skin[5]+0.5*alphaDrSurf[2]*G_1Skin[4]+0.8660254037844386*alphaDrSurf[1]*G_1Skin[3]+0.8660254037844386*alphaDrSurf[0]*G_1Skin[2]+0.5*G_1Skin[1]*alphaDrSurf[1]+0.5*G_1Skin[0]*alphaDrSurf[0]; 
  Ghat_G_1[1] = 1.0*alphaDrSurf[2]*G_1Skin[7]+1.118033988749895*alphaDrSurf[0]*G_1Skin[7]+0.7745966692414834*alphaDrSurf[1]*G_1Skin[6]+1.118033988749895*alphaDrSurf[1]*G_1Skin[5]+0.4472135954999579*alphaDrSurf[1]*G_1Skin[4]+0.7745966692414833*alphaDrSurf[2]*G_1Skin[3]+0.8660254037844386*alphaDrSurf[0]*G_1Skin[3]+0.4472135954999579*G_1Skin[1]*alphaDrSurf[2]+0.8660254037844386*alphaDrSurf[1]*G_1Skin[2]+0.5*G_1Skin[0]*alphaDrSurf[1]+0.5*alphaDrSurf[0]*G_1Skin[1]; 
  Ghat_G_1[2] = 1.0*alphaDrSurf[1]*G_1Skin[7]+0.5532833351724881*alphaDrSurf[2]*G_1Skin[6]+0.8660254037844387*alphaDrSurf[0]*G_1Skin[6]+1.118033988749895*alphaDrSurf[2]*G_1Skin[5]+0.31943828249997*alphaDrSurf[2]*G_1Skin[4]+0.5*alphaDrSurf[0]*G_1Skin[4]+0.7745966692414833*alphaDrSurf[1]*G_1Skin[3]+0.8660254037844386*G_1Skin[2]*alphaDrSurf[2]+0.5*G_1Skin[0]*alphaDrSurf[2]+0.4472135954999579*G_1Skin[1]*alphaDrSurf[1]; 

  out_F_0[0] += 0.7071067811865475*Ghat_F_0[0]*dv1par; 
  out_F_0[1] += 0.7071067811865475*Ghat_F_0[1]*dv1par; 
  out_F_0[2] += 1.224744871391589*Ghat_F_0[0]*dv1par; 
  out_F_0[3] += 1.224744871391589*Ghat_F_0[1]*dv1par; 
  out_F_0[4] += 0.7071067811865475*Ghat_F_0[2]*dv1par; 
  out_F_0[5] += 1.58113883008419*Ghat_F_0[0]*dv1par; 
  out_F_0[6] += 1.224744871391589*Ghat_F_0[2]*dv1par; 
  out_F_0[7] += 1.58113883008419*Ghat_F_0[1]*dv1par; 
  out_G_1[0] += 0.7071067811865475*Ghat_G_1[0]*dv1par; 
  out_G_1[1] += 0.7071067811865475*Ghat_G_1[1]*dv1par; 
  out_G_1[2] += 1.224744871391589*Ghat_G_1[0]*dv1par; 
  out_G_1[3] += 1.224744871391589*Ghat_G_1[1]*dv1par; 
  out_G_1[4] += 0.7071067811865475*Ghat_G_1[2]*dv1par; 
  out_G_1[5] += 1.58113883008419*Ghat_G_1[0]*dv1par; 
  out_G_1[6] += 1.224744871391589*Ghat_G_1[2]*dv1par; 
  out_G_1[7] += 1.58113883008419*Ghat_G_1[1]*dv1par; 

  } else { 

  alphaDrSurf[0] = nu[0]*wvpar-0.5*nu[0]*dvpar; 
  alphaDrSurf[1] = nu[1]*wvpar-0.5*nu[1]*dvpar; 
  alphaDrSurf[2] = nu[2]*wvpar-0.5*nu[2]*dvpar; 

  Ghat_F_0[0] = 1.118033988749895*alphaDrSurf[1]*F_0Skin[7]-0.8660254037844387*alphaDrSurf[2]*F_0Skin[6]+1.118033988749895*alphaDrSurf[0]*F_0Skin[5]+0.5*alphaDrSurf[2]*F_0Skin[4]-0.8660254037844386*alphaDrSurf[1]*F_0Skin[3]-0.8660254037844386*alphaDrSurf[0]*F_0Skin[2]+0.5*F_0Skin[1]*alphaDrSurf[1]+0.5*F_0Skin[0]*alphaDrSurf[0]; 
  Ghat_F_0[1] = 1.0*alphaDrSurf[2]*F_0Skin[7]+1.118033988749895*alphaDrSurf[0]*F_0Skin[7]-0.7745966692414834*alphaDrSurf[1]*F_0Skin[6]+1.118033988749895*alphaDrSurf[1]*F_0Skin[5]+0.4472135954999579*alphaDrSurf[1]*F_0Skin[4]-0.7745966692414833*alphaDrSurf[2]*F_0Skin[3]-0.8660254037844386*alphaDrSurf[0]*F_0Skin[3]+0.4472135954999579*F_0Skin[1]*alphaDrSurf[2]-0.8660254037844386*alphaDrSurf[1]*F_0Skin[2]+0.5*F_0Skin[0]*alphaDrSurf[1]+0.5*alphaDrSurf[0]*F_0Skin[1]; 
  Ghat_F_0[2] = 1.0*alphaDrSurf[1]*F_0Skin[7]-0.5532833351724881*alphaDrSurf[2]*F_0Skin[6]-0.8660254037844387*alphaDrSurf[0]*F_0Skin[6]+1.118033988749895*alphaDrSurf[2]*F_0Skin[5]+0.31943828249997*alphaDrSurf[2]*F_0Skin[4]+0.5*alphaDrSurf[0]*F_0Skin[4]-0.7745966692414833*alphaDrSurf[1]*F_0Skin[3]-0.8660254037844386*F_0Skin[2]*alphaDrSurf[2]+0.5*F_0Skin[0]*alphaDrSurf[2]+0.4472135954999579*F_0Skin[1]*alphaDrSurf[1]; 
  Ghat_G_1[0] = 1.118033988749895*alphaDrSurf[1]*G_1Skin[7]-0.8660254037844387*alphaDrSurf[2]*G_1Skin[6]+1.118033988749895*alphaDrSurf[0]*G_1Skin[5]+0.5*alphaDrSurf[2]*G_1Skin[4]-0.8660254037844386*alphaDrSurf[1]*G_1Skin[3]-0.8660254037844386*alphaDrSurf[0]*G_1Skin[2]+0.5*G_1Skin[1]*alphaDrSurf[1]+0.5*G_1Skin[0]*alphaDrSurf[0]; 
  Ghat_G_1[1] = 1.0*alphaDrSurf[2]*G_1Skin[7]+1.118033988749895*alphaDrSurf[0]*G_1Skin[7]-0.7745966692414834*alphaDrSurf[1]*G_1Skin[6]+1.118033988749895*alphaDrSurf[1]*G_1Skin[5]+0.4472135954999579*alphaDrSurf[1]*G_1Skin[4]-0.7745966692414833*alphaDrSurf[2]*G_1Skin[3]-0.8660254037844386*alphaDrSurf[0]*G_1Skin[3]+0.4472135954999579*G_1Skin[1]*alphaDrSurf[2]-0.8660254037844386*alphaDrSurf[1]*G_1Skin[2]+0.5*G_1Skin[0]*alphaDrSurf[1]+0.5*alphaDrSurf[0]*G_1Skin[1]; 
  Ghat_G_1[2] = 1.0*alphaDrSurf[1]*G_1Skin[7]-0.5532833351724881*alphaDrSurf[2]*G_1Skin[6]-0.8660254037844387*alphaDrSurf[0]*G_1Skin[6]+1.118033988749895*alphaDrSurf[2]*G_1Skin[5]+0.31943828249997*alphaDrSurf[2]*G_1Skin[4]+0.5*alphaDrSurf[0]*G_1Skin[4]-0.7745966692414833*alphaDrSurf[1]*G_1Skin[3]-0.8660254037844386*G_1Skin[2]*alphaDrSurf[2]+0.5*G_1Skin[0]*alphaDrSurf[2]+0.4472135954999579*G_1Skin[1]*alphaDrSurf[1]; 

  out_F_0[0] += -0.7071067811865475*Ghat_F_0[0]*dv1par; 
  out_F_0[1] += -0.7071067811865475*Ghat_F_0[1]*dv1par; 
  out_F_0[2] += 1.224744871391589*Ghat_F_0[0]*dv1par; 
  out_F_0[3] += 1.224744871391589*Ghat_F_0[1]*dv1par; 
  out_F_0[4] += -0.7071067811865475*Ghat_F_0[2]*dv1par; 
  out_F_0[5] += -1.58113883008419*Ghat_F_0[0]*dv1par; 
  out_F_0[6] += 1.224744871391589*Ghat_F_0[2]*dv1par; 
  out_F_0[7] += -1.58113883008419*Ghat_F_0[1]*dv1par; 
  out_G_1[0] += -0.7071067811865475*Ghat_G_1[0]*dv1par; 
  out_G_1[1] += -0.7071067811865475*Ghat_G_1[1]*dv1par; 
  out_G_1[2] += 1.224744871391589*Ghat_G_1[0]*dv1par; 
  out_G_1[3] += 1.224744871391589*Ghat_G_1[1]*dv1par; 
  out_G_1[4] += -0.7071067811865475*Ghat_G_1[2]*dv1par; 
  out_G_1[5] += -1.58113883008419*Ghat_G_1[0]*dv1par; 
  out_G_1[6] += 1.224744871391589*Ghat_G_1[2]*dv1par; 
  out_G_1[7] += -1.58113883008419*Ghat_G_1[1]*dv1par; 

  } 
} 
