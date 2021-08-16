#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_surfx_1x1v_tensor_p2(const double *w, const double *dxv, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fl/fc/fr:  Input Distribution function in left/center/right cells.
  // out:       Incremented distribution function in center cell.
  const double dx10 = 2/dxv[0]; 
  const double dv = dxv[1], wv = w[1]; 
  double Ghat_r[3]; 
  double Ghat_l[3]; 
  if (wv>0) { 

  Ghat_r[0] = (1.58113883008419*fc[4]+1.224744871391589*fc[1]+0.7071067811865475*fc[0])*wv+(0.4564354645876384*fc[6]+0.3535533905932737*fc[3]+0.2041241452319315*fc[2])*dv; 
  Ghat_r[1] = (1.58113883008419*fc[6]+1.224744871391589*fc[3]+0.7071067811865475*fc[2])*wv+(0.408248290463863*fc[8]+0.3162277660168379*fc[7]+0.1825741858350554*fc[5]+0.4564354645876384*fc[4]+0.3535533905932737*fc[1]+0.2041241452319315*fc[0])*dv; 
  Ghat_r[2] = (1.58113883008419*fc[8]+1.224744871391589*fc[7]+0.7071067811865475*fc[5])*wv+(0.408248290463863*fc[6]+0.3162277660168379*fc[3]+0.1825741858350554*fc[2])*dv; 

  Ghat_l[0] = (1.58113883008419*fl[4]+1.224744871391589*fl[1]+0.7071067811865475*fl[0])*wv+(0.4564354645876384*fl[6]+0.3535533905932737*fl[3]+0.2041241452319315*fl[2])*dv; 
  Ghat_l[1] = (1.58113883008419*fl[6]+1.224744871391589*fl[3]+0.7071067811865475*fl[2])*wv+(0.408248290463863*fl[8]+0.3162277660168379*fl[7]+0.1825741858350554*fl[5]+0.4564354645876384*fl[4]+0.3535533905932737*fl[1]+0.2041241452319315*fl[0])*dv; 
  Ghat_l[2] = (1.58113883008419*fl[8]+1.224744871391589*fl[7]+0.7071067811865475*fl[5])*wv+(0.408248290463863*fl[6]+0.3162277660168379*fl[3]+0.1825741858350554*fl[2])*dv; 

  } else { 

  Ghat_r[0] = (1.58113883008419*fr[4]-1.224744871391589*fr[1]+0.7071067811865475*fr[0])*wv+(0.4564354645876384*fr[6]-0.3535533905932737*fr[3]+0.2041241452319315*fr[2])*dv; 
  Ghat_r[1] = (1.58113883008419*fr[6]-1.224744871391589*fr[3]+0.7071067811865475*fr[2])*wv+(0.408248290463863*fr[8]-0.3162277660168379*fr[7]+0.1825741858350554*fr[5]+0.4564354645876384*fr[4]-0.3535533905932737*fr[1]+0.2041241452319315*fr[0])*dv; 
  Ghat_r[2] = (1.58113883008419*fr[8]-1.224744871391589*fr[7]+0.7071067811865475*fr[5])*wv+(0.408248290463863*fr[6]-0.3162277660168379*fr[3]+0.1825741858350554*fr[2])*dv; 

  Ghat_l[0] = (1.58113883008419*fc[4]-1.224744871391589*fc[1]+0.7071067811865475*fc[0])*wv+(0.4564354645876384*fc[6]-0.3535533905932737*fc[3]+0.2041241452319315*fc[2])*dv; 
  Ghat_l[1] = (1.58113883008419*fc[6]-1.224744871391589*fc[3]+0.7071067811865475*fc[2])*wv+(0.408248290463863*fc[8]-0.3162277660168379*fc[7]+0.1825741858350554*fc[5]+0.4564354645876384*fc[4]-0.3535533905932737*fc[1]+0.2041241452319315*fc[0])*dv; 
  Ghat_l[2] = (1.58113883008419*fc[8]-1.224744871391589*fc[7]+0.7071067811865475*fc[5])*wv+(0.408248290463863*fc[6]-0.3162277660168379*fc[3]+0.1825741858350554*fc[2])*dv; 

  } 
  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dx10; 
  out[1] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dx10; 
  out[2] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dx10; 
  out[3] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dx10; 
  out[4] += (1.58113883008419*Ghat_l[0]-1.58113883008419*Ghat_r[0])*dx10; 
  out[5] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dx10; 
  out[6] += (1.58113883008419*Ghat_l[1]-1.58113883008419*Ghat_r[1])*dx10; 
  out[7] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dx10; 
  out[8] += (1.58113883008419*Ghat_l[2]-1.58113883008419*Ghat_r[2])*dx10; 
} 