#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_boundary_surfx_1x2v_ser_p2(const double *w, const double *dxv, const double *alpha_geo, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:     Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // alpha_geo:   Fields used only for general geometry.
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Incremented distribution function in center cell.
  const double dx10 = 2/dxv[0]; 
  const double dv = dxv[1], wv = w[1]; 
  double Ghat[8]; 

  if (edge == -1) { 

  if (wv>0) { 

  Ghat[0] = (1.58113883008419*fSkin[7]+1.224744871391589*fSkin[1]+0.7071067811865475*fSkin[0])*wv+(0.4564354645876384*fSkin[11]+0.3535533905932737*fSkin[4]+0.2041241452319315*fSkin[2])*dv; 
  Ghat[1] = (1.58113883008419*fSkin[11]+1.224744871391589*fSkin[4]+0.7071067811865475*fSkin[2])*wv+(0.3162277660168379*fSkin[12]+0.1825741858350554*fSkin[8]+0.4564354645876384*fSkin[7]+0.3535533905932737*fSkin[1]+0.2041241452319315*fSkin[0])*dv; 
  Ghat[2] = (1.58113883008419*fSkin[13]+1.224744871391589*fSkin[5]+0.7071067811865475*fSkin[3])*wv+(0.4564354645876384*fSkin[17]+0.3535533905932737*fSkin[10]+0.2041241452319315*fSkin[6])*dv; 
  Ghat[3] = (1.58113883008419*fSkin[17]+1.224744871391589*fSkin[10]+0.7071067811865475*fSkin[6])*wv+(0.3162277660168379*fSkin[18]+0.1825741858350554*fSkin[14]+0.4564354645876384*fSkin[13]+0.3535533905932737*fSkin[5]+0.2041241452319315*fSkin[3])*dv; 
  Ghat[4] = (1.224744871391589*fSkin[12]+0.7071067811865475*fSkin[8])*wv+(0.408248290463863*fSkin[11]+0.3162277660168379*fSkin[4]+0.1825741858350554*fSkin[2])*dv; 
  Ghat[5] = (1.224744871391589*fSkin[15]+0.7071067811865475*fSkin[9])*wv+(0.3535533905932737*fSkin[19]+0.2041241452319315*fSkin[16])*dv; 
  Ghat[6] = (1.224744871391589*fSkin[18]+0.7071067811865475*fSkin[14])*wv+(0.408248290463863*fSkin[17]+0.3162277660168379*fSkin[10]+0.1825741858350554*fSkin[6])*dv; 
  Ghat[7] = (1.224744871391589*fSkin[19]+0.7071067811865475*fSkin[16])*wv+(0.3535533905932737*fSkin[15]+0.2041241452319315*fSkin[9])*dv; 

  } else { 

  Ghat[0] = 1.58113883008419*fEdge[7]*wv-1.224744871391589*fEdge[1]*wv+0.7071067811865475*fEdge[0]*wv+0.4564354645876383*fEdge[11]*dv-0.3535533905932737*fEdge[4]*dv+0.2041241452319315*fEdge[2]*dv; 
  Ghat[1] = 1.58113883008419*fEdge[11]*wv-1.224744871391589*fEdge[4]*wv+0.7071067811865475*fEdge[2]*wv-0.3162277660168379*fEdge[12]*dv+0.1825741858350554*fEdge[8]*dv+0.4564354645876384*fEdge[7]*dv-0.3535533905932737*fEdge[1]*dv+0.2041241452319315*fEdge[0]*dv; 
  Ghat[2] = 1.58113883008419*fEdge[13]*wv-1.224744871391589*fEdge[5]*wv+0.7071067811865475*fEdge[3]*wv+0.4564354645876384*fEdge[17]*dv-0.3535533905932737*fEdge[10]*dv+0.2041241452319315*fEdge[6]*dv; 
  Ghat[3] = 1.58113883008419*fEdge[17]*wv-1.224744871391589*fEdge[10]*wv+0.7071067811865475*fEdge[6]*wv-0.3162277660168379*fEdge[18]*dv+0.1825741858350553*fEdge[14]*dv+0.4564354645876383*fEdge[13]*dv-0.3535533905932737*fEdge[5]*dv+0.2041241452319315*fEdge[3]*dv; 
  Ghat[4] = (-1.224744871391589*fEdge[12]*wv)+0.7071067811865475*fEdge[8]*wv+0.408248290463863*fEdge[11]*dv-0.3162277660168379*fEdge[4]*dv+0.1825741858350554*fEdge[2]*dv; 
  Ghat[5] = (-1.224744871391589*fEdge[15]*wv)+0.7071067811865475*fEdge[9]*wv-0.3535533905932737*fEdge[19]*dv+0.2041241452319315*fEdge[16]*dv; 
  Ghat[6] = (-1.224744871391589*fEdge[18]*wv)+0.7071067811865475*fEdge[14]*wv+0.408248290463863*fEdge[17]*dv-0.3162277660168379*fEdge[10]*dv+0.1825741858350553*fEdge[6]*dv; 
  Ghat[7] = (-1.224744871391589*fEdge[19]*wv)+0.7071067811865475*fEdge[16]*wv-0.3535533905932737*fEdge[15]*dv+0.2041241452319315*fEdge[9]*dv; 

  } 

  out[0] += 0.7071067811865475*Ghat[0]*dx10; 
  out[1] += 1.224744871391589*Ghat[0]*dx10; 
  out[2] += 0.7071067811865475*Ghat[1]*dx10; 
  out[3] += 0.7071067811865475*Ghat[2]*dx10; 
  out[4] += 1.224744871391589*Ghat[1]*dx10; 
  out[5] += 1.224744871391589*Ghat[2]*dx10; 
  out[6] += 0.7071067811865475*Ghat[3]*dx10; 
  out[7] += 1.58113883008419*Ghat[0]*dx10; 
  out[8] += 0.7071067811865475*Ghat[4]*dx10; 
  out[9] += 0.7071067811865475*Ghat[5]*dx10; 
  out[10] += 1.224744871391589*Ghat[3]*dx10; 
  out[11] += 1.58113883008419*Ghat[1]*dx10; 
  out[12] += 1.224744871391589*Ghat[4]*dx10; 
  out[13] += 1.58113883008419*Ghat[2]*dx10; 
  out[14] += 0.7071067811865475*Ghat[6]*dx10; 
  out[15] += 1.224744871391589*Ghat[5]*dx10; 
  out[16] += 0.7071067811865475*Ghat[7]*dx10; 
  out[17] += 1.58113883008419*Ghat[3]*dx10; 
  out[18] += 1.224744871391589*Ghat[6]*dx10; 
  out[19] += 1.224744871391589*Ghat[7]*dx10; 

  } else { 

  if (wv>0) { 

  Ghat[0] = (1.58113883008419*fEdge[7]+1.224744871391589*fEdge[1]+0.7071067811865475*fEdge[0])*wv+(0.4564354645876384*fEdge[11]+0.3535533905932737*fEdge[4]+0.2041241452319315*fEdge[2])*dv; 
  Ghat[1] = (1.58113883008419*fEdge[11]+1.224744871391589*fEdge[4]+0.7071067811865475*fEdge[2])*wv+(0.3162277660168379*fEdge[12]+0.1825741858350554*fEdge[8]+0.4564354645876384*fEdge[7]+0.3535533905932737*fEdge[1]+0.2041241452319315*fEdge[0])*dv; 
  Ghat[2] = (1.58113883008419*fEdge[13]+1.224744871391589*fEdge[5]+0.7071067811865475*fEdge[3])*wv+(0.4564354645876384*fEdge[17]+0.3535533905932737*fEdge[10]+0.2041241452319315*fEdge[6])*dv; 
  Ghat[3] = (1.58113883008419*fEdge[17]+1.224744871391589*fEdge[10]+0.7071067811865475*fEdge[6])*wv+(0.3162277660168379*fEdge[18]+0.1825741858350554*fEdge[14]+0.4564354645876384*fEdge[13]+0.3535533905932737*fEdge[5]+0.2041241452319315*fEdge[3])*dv; 
  Ghat[4] = (1.224744871391589*fEdge[12]+0.7071067811865475*fEdge[8])*wv+(0.408248290463863*fEdge[11]+0.3162277660168379*fEdge[4]+0.1825741858350554*fEdge[2])*dv; 
  Ghat[5] = (1.224744871391589*fEdge[15]+0.7071067811865475*fEdge[9])*wv+(0.3535533905932737*fEdge[19]+0.2041241452319315*fEdge[16])*dv; 
  Ghat[6] = (1.224744871391589*fEdge[18]+0.7071067811865475*fEdge[14])*wv+(0.408248290463863*fEdge[17]+0.3162277660168379*fEdge[10]+0.1825741858350554*fEdge[6])*dv; 
  Ghat[7] = (1.224744871391589*fEdge[19]+0.7071067811865475*fEdge[16])*wv+(0.3535533905932737*fEdge[15]+0.2041241452319315*fEdge[9])*dv; 

  } else { 

  Ghat[0] = 1.58113883008419*fSkin[7]*wv-1.224744871391589*fSkin[1]*wv+0.7071067811865475*fSkin[0]*wv+0.4564354645876383*fSkin[11]*dv-0.3535533905932737*fSkin[4]*dv+0.2041241452319315*fSkin[2]*dv; 
  Ghat[1] = 1.58113883008419*fSkin[11]*wv-1.224744871391589*fSkin[4]*wv+0.7071067811865475*fSkin[2]*wv-0.3162277660168379*fSkin[12]*dv+0.1825741858350554*fSkin[8]*dv+0.4564354645876384*fSkin[7]*dv-0.3535533905932737*fSkin[1]*dv+0.2041241452319315*fSkin[0]*dv; 
  Ghat[2] = 1.58113883008419*fSkin[13]*wv-1.224744871391589*fSkin[5]*wv+0.7071067811865475*fSkin[3]*wv+0.4564354645876384*fSkin[17]*dv-0.3535533905932737*fSkin[10]*dv+0.2041241452319315*fSkin[6]*dv; 
  Ghat[3] = 1.58113883008419*fSkin[17]*wv-1.224744871391589*fSkin[10]*wv+0.7071067811865475*fSkin[6]*wv-0.3162277660168379*fSkin[18]*dv+0.1825741858350553*fSkin[14]*dv+0.4564354645876383*fSkin[13]*dv-0.3535533905932737*fSkin[5]*dv+0.2041241452319315*fSkin[3]*dv; 
  Ghat[4] = (-1.224744871391589*fSkin[12]*wv)+0.7071067811865475*fSkin[8]*wv+0.408248290463863*fSkin[11]*dv-0.3162277660168379*fSkin[4]*dv+0.1825741858350554*fSkin[2]*dv; 
  Ghat[5] = (-1.224744871391589*fSkin[15]*wv)+0.7071067811865475*fSkin[9]*wv-0.3535533905932737*fSkin[19]*dv+0.2041241452319315*fSkin[16]*dv; 
  Ghat[6] = (-1.224744871391589*fSkin[18]*wv)+0.7071067811865475*fSkin[14]*wv+0.408248290463863*fSkin[17]*dv-0.3162277660168379*fSkin[10]*dv+0.1825741858350553*fSkin[6]*dv; 
  Ghat[7] = (-1.224744871391589*fSkin[19]*wv)+0.7071067811865475*fSkin[16]*wv-0.3535533905932737*fSkin[15]*dv+0.2041241452319315*fSkin[9]*dv; 

  } 

  out[0] += 0.7071067811865475*Ghat[0]*dx10; 
  out[1] += -1.224744871391589*Ghat[0]*dx10; 
  out[2] += 0.7071067811865475*Ghat[1]*dx10; 
  out[3] += 0.7071067811865475*Ghat[2]*dx10; 
  out[4] += -1.224744871391589*Ghat[1]*dx10; 
  out[5] += -1.224744871391589*Ghat[2]*dx10; 
  out[6] += 0.7071067811865475*Ghat[3]*dx10; 
  out[7] += 1.58113883008419*Ghat[0]*dx10; 
  out[8] += 0.7071067811865475*Ghat[4]*dx10; 
  out[9] += 0.7071067811865475*Ghat[5]*dx10; 
  out[10] += -1.224744871391589*Ghat[3]*dx10; 
  out[11] += 1.58113883008419*Ghat[1]*dx10; 
  out[12] += -1.224744871391589*Ghat[4]*dx10; 
  out[13] += 1.58113883008419*Ghat[2]*dx10; 
  out[14] += 0.7071067811865475*Ghat[6]*dx10; 
  out[15] += -1.224744871391589*Ghat[5]*dx10; 
  out[16] += 0.7071067811865475*Ghat[7]*dx10; 
  out[17] += 1.58113883008419*Ghat[3]*dx10; 
  out[18] += -1.224744871391589*Ghat[6]*dx10; 
  out[19] += -1.224744871391589*Ghat[7]*dx10; 

  } 
} 
