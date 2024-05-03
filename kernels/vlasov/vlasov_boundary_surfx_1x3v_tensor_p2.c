#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_boundary_surfx_1x3v_tensor_p2(const double *w, const double *dxv, const double *alpha_geo, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:     Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // alpha_geo:   Fields used only for general geometry.
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Incremented distribution function in center cell.
  const double dx10 = 2/dxv[0]; 
  const double dv = dxv[1], wv = w[1]; 
  double Ghat[27]; 

  if (edge == -1) { 

  if (wv>0) { 

  Ghat[0] = (1.58113883008419*fSkin[11]+1.224744871391589*fSkin[1]+0.7071067811865475*fSkin[0])*wv+(0.4564354645876384*fSkin[19]+0.3535533905932737*fSkin[5]+0.2041241452319315*fSkin[2])*dv; 
  Ghat[1] = (1.58113883008419*fSkin[19]+1.224744871391589*fSkin[5]+0.7071067811865475*fSkin[2])*wv+(0.408248290463863*fSkin[44]+0.3162277660168379*fSkin[20]+0.1825741858350554*fSkin[12]+0.4564354645876384*fSkin[11]+0.3535533905932737*fSkin[1]+0.2041241452319315*fSkin[0])*dv; 
  Ghat[2] = (1.58113883008419*fSkin[21]+1.224744871391589*fSkin[6]+0.7071067811865475*fSkin[3])*wv+(0.4564354645876384*fSkin[32]+0.3535533905932737*fSkin[15]+0.2041241452319315*fSkin[7])*dv; 
  Ghat[3] = (1.58113883008419*fSkin[25]+1.224744871391589*fSkin[8]+0.7071067811865475*fSkin[4])*wv+(0.4564354645876384*fSkin[35]+0.3535533905932737*fSkin[16]+0.2041241452319315*fSkin[9])*dv; 
  Ghat[4] = (1.58113883008419*fSkin[32]+1.224744871391589*fSkin[15]+0.7071067811865475*fSkin[7])*wv+(0.408248290463863*fSkin[54]+0.3162277660168379*fSkin[33]+0.1825741858350554*fSkin[22]+0.4564354645876384*fSkin[21]+0.3535533905932737*fSkin[6]+0.2041241452319315*fSkin[3])*dv; 
  Ghat[5] = (1.58113883008419*fSkin[35]+1.224744871391589*fSkin[16]+0.7071067811865475*fSkin[9])*wv+(0.408248290463863*fSkin[57]+0.3162277660168379*fSkin[36]+0.1825741858350554*fSkin[26]+0.4564354645876384*fSkin[25]+0.3535533905932737*fSkin[8]+0.2041241452319315*fSkin[4])*dv; 
  Ghat[6] = (1.58113883008419*fSkin[37]+1.224744871391589*fSkin[17]+0.7071067811865475*fSkin[10])*wv+(0.4564354645876384*fSkin[50]+0.3535533905932737*fSkin[31]+0.2041241452319315*fSkin[18])*dv; 
  Ghat[7] = (1.58113883008419*fSkin[44]+1.224744871391589*fSkin[20]+0.7071067811865475*fSkin[12])*wv+(0.408248290463863*fSkin[19]+0.3162277660168379*fSkin[5]+0.1825741858350554*fSkin[2])*dv; 
  Ghat[8] = (1.58113883008419*fSkin[45]+1.224744871391589*fSkin[23]+0.7071067811865475*fSkin[13])*wv+(0.4564354645876384*fSkin[55]+0.3535533905932737*fSkin[34]+0.2041241452319315*fSkin[24])*dv; 
  Ghat[9] = (1.58113883008419*fSkin[47]+1.224744871391589*fSkin[28]+0.7071067811865475*fSkin[14])*wv+(0.4564354645876384*fSkin[60]+0.3535533905932737*fSkin[41]+0.2041241452319315*fSkin[29])*dv; 
  Ghat[10] = (1.58113883008419*fSkin[50]+1.224744871391589*fSkin[31]+0.7071067811865475*fSkin[18])*wv+(0.408248290463863*fSkin[66]+0.3162277660168379*fSkin[51]+0.1825741858350554*fSkin[38]+0.4564354645876384*fSkin[37]+0.3535533905932737*fSkin[17]+0.2041241452319315*fSkin[10])*dv; 
  Ghat[11] = (1.58113883008419*fSkin[54]+1.224744871391589*fSkin[33]+0.7071067811865475*fSkin[22])*wv+(0.408248290463863*fSkin[32]+0.3162277660168379*fSkin[15]+0.1825741858350554*fSkin[7])*dv; 
  Ghat[12] = (1.58113883008419*fSkin[55]+1.224744871391589*fSkin[34]+0.7071067811865475*fSkin[24])*wv+(0.408248290463863*fSkin[72]+0.3162277660168379*fSkin[56]+0.1825741858350554*fSkin[46]+0.4564354645876384*fSkin[45]+0.3535533905932737*fSkin[23]+0.2041241452319315*fSkin[13])*dv; 
  Ghat[13] = (1.58113883008419*fSkin[57]+1.224744871391589*fSkin[36]+0.7071067811865475*fSkin[26])*wv+(0.408248290463863*fSkin[35]+0.3162277660168379*fSkin[16]+0.1825741858350554*fSkin[9])*dv; 
  Ghat[14] = (1.58113883008419*fSkin[58]+1.224744871391589*fSkin[39]+0.7071067811865475*fSkin[27])*wv+(0.4564354645876384*fSkin[67]+0.3535533905932737*fSkin[52]+0.2041241452319315*fSkin[40])*dv; 
  Ghat[15] = (1.58113883008419*fSkin[60]+1.224744871391589*fSkin[41]+0.7071067811865475*fSkin[29])*wv+(0.408248290463863*fSkin[73]+0.3162277660168379*fSkin[61]+0.1825741858350554*fSkin[48]+0.4564354645876384*fSkin[47]+0.3535533905932737*fSkin[28]+0.2041241452319315*fSkin[14])*dv; 
  Ghat[16] = (1.58113883008419*fSkin[62]+1.224744871391589*fSkin[42]+0.7071067811865475*fSkin[30])*wv+(0.4564354645876384*fSkin[69]+0.3535533905932737*fSkin[53]+0.2041241452319315*fSkin[43])*dv; 
  Ghat[17] = (1.58113883008419*fSkin[66]+1.224744871391589*fSkin[51]+0.7071067811865475*fSkin[38])*wv+(0.408248290463863*fSkin[50]+0.3162277660168379*fSkin[31]+0.1825741858350554*fSkin[18])*dv; 
  Ghat[18] = (1.58113883008419*fSkin[67]+1.224744871391589*fSkin[52]+0.7071067811865475*fSkin[40])*wv+(0.408248290463863*fSkin[76]+0.3162277660168379*fSkin[68]+0.1825741858350554*fSkin[59]+0.4564354645876384*fSkin[58]+0.3535533905932737*fSkin[39]+0.2041241452319315*fSkin[27])*dv; 
  Ghat[19] = (1.58113883008419*fSkin[69]+1.224744871391589*fSkin[53]+0.7071067811865475*fSkin[43])*wv+(0.408248290463863*fSkin[77]+0.3162277660168379*fSkin[70]+0.1825741858350554*fSkin[63]+0.4564354645876384*fSkin[62]+0.3535533905932737*fSkin[42]+0.2041241452319315*fSkin[30])*dv; 
  Ghat[20] = (1.58113883008419*fSkin[72]+1.224744871391589*fSkin[56]+0.7071067811865475*fSkin[46])*wv+(0.408248290463863*fSkin[55]+0.3162277660168379*fSkin[34]+0.1825741858350554*fSkin[24])*dv; 
  Ghat[21] = (1.58113883008419*fSkin[73]+1.224744871391589*fSkin[61]+0.7071067811865475*fSkin[48])*wv+(0.408248290463863*fSkin[60]+0.3162277660168379*fSkin[41]+0.1825741858350554*fSkin[29])*dv; 
  Ghat[22] = (1.58113883008419*fSkin[74]+1.224744871391589*fSkin[64]+0.7071067811865475*fSkin[49])*wv+(0.4564354645876384*fSkin[78]+0.3535533905932737*fSkin[71]+0.2041241452319315*fSkin[65])*dv; 
  Ghat[23] = (1.58113883008419*fSkin[76]+1.224744871391589*fSkin[68]+0.7071067811865475*fSkin[59])*wv+(0.408248290463863*fSkin[67]+0.3162277660168379*fSkin[52]+0.1825741858350554*fSkin[40])*dv; 
  Ghat[24] = (1.58113883008419*fSkin[77]+1.224744871391589*fSkin[70]+0.7071067811865475*fSkin[63])*wv+(0.408248290463863*fSkin[69]+0.3162277660168379*fSkin[53]+0.1825741858350554*fSkin[43])*dv; 
  Ghat[25] = (1.58113883008419*fSkin[78]+1.224744871391589*fSkin[71]+0.7071067811865475*fSkin[65])*wv+(0.408248290463863*fSkin[80]+0.3162277660168379*fSkin[79]+0.1825741858350554*fSkin[75]+0.4564354645876384*fSkin[74]+0.3535533905932737*fSkin[64]+0.2041241452319315*fSkin[49])*dv; 
  Ghat[26] = (1.58113883008419*fSkin[80]+1.224744871391589*fSkin[79]+0.7071067811865475*fSkin[75])*wv+(0.408248290463863*fSkin[78]+0.3162277660168379*fSkin[71]+0.1825741858350554*fSkin[65])*dv; 

  } else { 

  Ghat[0] = 1.58113883008419*fEdge[11]*wv-1.224744871391589*fEdge[1]*wv+0.7071067811865475*fEdge[0]*wv+0.4564354645876383*fEdge[19]*dv-0.3535533905932737*fEdge[5]*dv+0.2041241452319315*fEdge[2]*dv; 
  Ghat[1] = 1.58113883008419*fEdge[19]*wv-1.224744871391589*fEdge[5]*wv+0.7071067811865475*fEdge[2]*wv+0.408248290463863*fEdge[44]*dv-0.3162277660168379*fEdge[20]*dv+0.1825741858350554*fEdge[12]*dv+0.4564354645876384*fEdge[11]*dv-0.3535533905932737*fEdge[1]*dv+0.2041241452319315*fEdge[0]*dv; 
  Ghat[2] = 1.58113883008419*fEdge[21]*wv-1.224744871391589*fEdge[6]*wv+0.7071067811865475*fEdge[3]*wv+0.4564354645876384*fEdge[32]*dv-0.3535533905932737*fEdge[15]*dv+0.2041241452319315*fEdge[7]*dv; 
  Ghat[3] = 1.58113883008419*fEdge[25]*wv-1.224744871391589*fEdge[8]*wv+0.7071067811865475*fEdge[4]*wv+0.4564354645876384*fEdge[35]*dv-0.3535533905932737*fEdge[16]*dv+0.2041241452319315*fEdge[9]*dv; 
  Ghat[4] = 1.58113883008419*fEdge[32]*wv-1.224744871391589*fEdge[15]*wv+0.7071067811865475*fEdge[7]*wv+0.408248290463863*fEdge[54]*dv-0.3162277660168379*fEdge[33]*dv+0.1825741858350553*fEdge[22]*dv+0.4564354645876383*fEdge[21]*dv-0.3535533905932737*fEdge[6]*dv+0.2041241452319315*fEdge[3]*dv; 
  Ghat[5] = 1.58113883008419*fEdge[35]*wv-1.224744871391589*fEdge[16]*wv+0.7071067811865475*fEdge[9]*wv+0.408248290463863*fEdge[57]*dv-0.3162277660168379*fEdge[36]*dv+0.1825741858350553*fEdge[26]*dv+0.4564354645876383*fEdge[25]*dv-0.3535533905932737*fEdge[8]*dv+0.2041241452319315*fEdge[4]*dv; 
  Ghat[6] = 1.58113883008419*fEdge[37]*wv-1.224744871391589*fEdge[17]*wv+0.7071067811865475*fEdge[10]*wv+0.4564354645876383*fEdge[50]*dv-0.3535533905932737*fEdge[31]*dv+0.2041241452319315*fEdge[18]*dv; 
  Ghat[7] = 1.58113883008419*fEdge[44]*wv-1.224744871391589*fEdge[20]*wv+0.7071067811865475*fEdge[12]*wv+0.408248290463863*fEdge[19]*dv-0.3162277660168379*fEdge[5]*dv+0.1825741858350554*fEdge[2]*dv; 
  Ghat[8] = 1.58113883008419*fEdge[45]*wv-1.224744871391589*fEdge[23]*wv+0.7071067811865475*fEdge[13]*wv+0.4564354645876384*fEdge[55]*dv-0.3535533905932737*fEdge[34]*dv+0.2041241452319315*fEdge[24]*dv; 
  Ghat[9] = 1.58113883008419*fEdge[47]*wv-1.224744871391589*fEdge[28]*wv+0.7071067811865475*fEdge[14]*wv+0.4564354645876384*fEdge[60]*dv-0.3535533905932737*fEdge[41]*dv+0.2041241452319315*fEdge[29]*dv; 
  Ghat[10] = 1.58113883008419*fEdge[50]*wv-1.224744871391589*fEdge[31]*wv+0.7071067811865475*fEdge[18]*wv+0.408248290463863*fEdge[66]*dv-0.3162277660168379*fEdge[51]*dv+0.1825741858350554*fEdge[38]*dv+0.4564354645876384*fEdge[37]*dv-0.3535533905932737*fEdge[17]*dv+0.2041241452319315*fEdge[10]*dv; 
  Ghat[11] = 1.58113883008419*fEdge[54]*wv-1.224744871391589*fEdge[33]*wv+0.7071067811865475*fEdge[22]*wv+0.408248290463863*fEdge[32]*dv-0.3162277660168379*fEdge[15]*dv+0.1825741858350553*fEdge[7]*dv; 
  Ghat[12] = 1.58113883008419*fEdge[55]*wv-1.224744871391589*fEdge[34]*wv+0.7071067811865475*fEdge[24]*wv+0.408248290463863*fEdge[72]*dv-0.3162277660168379*fEdge[56]*dv+0.1825741858350553*fEdge[46]*dv+0.4564354645876383*fEdge[45]*dv-0.3535533905932737*fEdge[23]*dv+0.2041241452319315*fEdge[13]*dv; 
  Ghat[13] = 1.58113883008419*fEdge[57]*wv-1.224744871391589*fEdge[36]*wv+0.7071067811865475*fEdge[26]*wv+0.408248290463863*fEdge[35]*dv-0.3162277660168379*fEdge[16]*dv+0.1825741858350553*fEdge[9]*dv; 
  Ghat[14] = 1.58113883008419*fEdge[58]*wv-1.224744871391589*fEdge[39]*wv+0.7071067811865475*fEdge[27]*wv+0.4564354645876383*fEdge[67]*dv-0.3535533905932737*fEdge[52]*dv+0.2041241452319315*fEdge[40]*dv; 
  Ghat[15] = 1.58113883008419*fEdge[60]*wv-1.224744871391589*fEdge[41]*wv+0.7071067811865475*fEdge[29]*wv+0.408248290463863*fEdge[73]*dv-0.3162277660168379*fEdge[61]*dv+0.1825741858350553*fEdge[48]*dv+0.4564354645876383*fEdge[47]*dv-0.3535533905932737*fEdge[28]*dv+0.2041241452319315*fEdge[14]*dv; 
  Ghat[16] = 1.58113883008419*fEdge[62]*wv-1.224744871391589*fEdge[42]*wv+0.7071067811865475*fEdge[30]*wv+0.4564354645876383*fEdge[69]*dv-0.3535533905932737*fEdge[53]*dv+0.2041241452319315*fEdge[43]*dv; 
  Ghat[17] = 1.58113883008419*fEdge[66]*wv-1.224744871391589*fEdge[51]*wv+0.7071067811865475*fEdge[38]*wv+0.408248290463863*fEdge[50]*dv-0.3162277660168379*fEdge[31]*dv+0.1825741858350554*fEdge[18]*dv; 
  Ghat[18] = 1.58113883008419*fEdge[67]*wv-1.224744871391589*fEdge[52]*wv+0.7071067811865475*fEdge[40]*wv+0.408248290463863*fEdge[76]*dv-0.3162277660168379*fEdge[68]*dv+0.1825741858350554*fEdge[59]*dv+0.4564354645876384*fEdge[58]*dv-0.3535533905932737*fEdge[39]*dv+0.2041241452319315*fEdge[27]*dv; 
  Ghat[19] = 1.58113883008419*fEdge[69]*wv-1.224744871391589*fEdge[53]*wv+0.7071067811865475*fEdge[43]*wv+0.408248290463863*fEdge[77]*dv-0.3162277660168379*fEdge[70]*dv+0.1825741858350554*fEdge[63]*dv+0.4564354645876384*fEdge[62]*dv-0.3535533905932737*fEdge[42]*dv+0.2041241452319315*fEdge[30]*dv; 
  Ghat[20] = 1.58113883008419*fEdge[72]*wv-1.224744871391589*fEdge[56]*wv+0.7071067811865475*fEdge[46]*wv+0.408248290463863*fEdge[55]*dv-0.3162277660168379*fEdge[34]*dv+0.1825741858350553*fEdge[24]*dv; 
  Ghat[21] = 1.58113883008419*fEdge[73]*wv-1.224744871391589*fEdge[61]*wv+0.7071067811865475*fEdge[48]*wv+0.408248290463863*fEdge[60]*dv-0.3162277660168379*fEdge[41]*dv+0.1825741858350553*fEdge[29]*dv; 
  Ghat[22] = 1.58113883008419*fEdge[74]*wv-1.224744871391589*fEdge[64]*wv+0.7071067811865475*fEdge[49]*wv+0.4564354645876383*fEdge[78]*dv-0.3535533905932737*fEdge[71]*dv+0.2041241452319315*fEdge[65]*dv; 
  Ghat[23] = 1.58113883008419*fEdge[76]*wv-1.224744871391589*fEdge[68]*wv+0.7071067811865475*fEdge[59]*wv+0.408248290463863*fEdge[67]*dv-0.3162277660168379*fEdge[52]*dv+0.1825741858350554*fEdge[40]*dv; 
  Ghat[24] = 1.58113883008419*fEdge[77]*wv-1.224744871391589*fEdge[70]*wv+0.7071067811865475*fEdge[63]*wv+0.408248290463863*fEdge[69]*dv-0.3162277660168379*fEdge[53]*dv+0.1825741858350554*fEdge[43]*dv; 
  Ghat[25] = 1.58113883008419*fEdge[78]*wv-1.224744871391589*fEdge[71]*wv+0.7071067811865475*fEdge[65]*wv+0.408248290463863*fEdge[80]*dv-0.3162277660168379*fEdge[79]*dv+0.1825741858350554*fEdge[75]*dv+0.4564354645876384*fEdge[74]*dv-0.3535533905932737*fEdge[64]*dv+0.2041241452319315*fEdge[49]*dv; 
  Ghat[26] = 1.58113883008419*fEdge[80]*wv-1.224744871391589*fEdge[79]*wv+0.7071067811865475*fEdge[75]*wv+0.408248290463863*fEdge[78]*dv-0.3162277660168379*fEdge[71]*dv+0.1825741858350554*fEdge[65]*dv; 

  } 

  out[0] += -0.7071067811865475*Ghat[0]*dx10; 
  out[1] += -1.224744871391589*Ghat[0]*dx10; 
  out[2] += -0.7071067811865475*Ghat[1]*dx10; 
  out[3] += -0.7071067811865475*Ghat[2]*dx10; 
  out[4] += -0.7071067811865475*Ghat[3]*dx10; 
  out[5] += -1.224744871391589*Ghat[1]*dx10; 
  out[6] += -1.224744871391589*Ghat[2]*dx10; 
  out[7] += -0.7071067811865475*Ghat[4]*dx10; 
  out[8] += -1.224744871391589*Ghat[3]*dx10; 
  out[9] += -0.7071067811865475*Ghat[5]*dx10; 
  out[10] += -0.7071067811865475*Ghat[6]*dx10; 
  out[11] += -1.58113883008419*Ghat[0]*dx10; 
  out[12] += -0.7071067811865475*Ghat[7]*dx10; 
  out[13] += -0.7071067811865475*Ghat[8]*dx10; 
  out[14] += -0.7071067811865475*Ghat[9]*dx10; 
  out[15] += -1.224744871391589*Ghat[4]*dx10; 
  out[16] += -1.224744871391589*Ghat[5]*dx10; 
  out[17] += -1.224744871391589*Ghat[6]*dx10; 
  out[18] += -0.7071067811865475*Ghat[10]*dx10; 
  out[19] += -1.58113883008419*Ghat[1]*dx10; 
  out[20] += -1.224744871391589*Ghat[7]*dx10; 
  out[21] += -1.58113883008419*Ghat[2]*dx10; 
  out[22] += -0.7071067811865475*Ghat[11]*dx10; 
  out[23] += -1.224744871391589*Ghat[8]*dx10; 
  out[24] += -0.7071067811865475*Ghat[12]*dx10; 
  out[25] += -1.58113883008419*Ghat[3]*dx10; 
  out[26] += -0.7071067811865475*Ghat[13]*dx10; 
  out[27] += -0.7071067811865475*Ghat[14]*dx10; 
  out[28] += -1.224744871391589*Ghat[9]*dx10; 
  out[29] += -0.7071067811865475*Ghat[15]*dx10; 
  out[30] += -0.7071067811865475*Ghat[16]*dx10; 
  out[31] += -1.224744871391589*Ghat[10]*dx10; 
  out[32] += -1.58113883008419*Ghat[4]*dx10; 
  out[33] += -1.224744871391589*Ghat[11]*dx10; 
  out[34] += -1.224744871391589*Ghat[12]*dx10; 
  out[35] += -1.58113883008419*Ghat[5]*dx10; 
  out[36] += -1.224744871391589*Ghat[13]*dx10; 
  out[37] += -1.58113883008419*Ghat[6]*dx10; 
  out[38] += -0.7071067811865475*Ghat[17]*dx10; 
  out[39] += -1.224744871391589*Ghat[14]*dx10; 
  out[40] += -0.7071067811865475*Ghat[18]*dx10; 
  out[41] += -1.224744871391589*Ghat[15]*dx10; 
  out[42] += -1.224744871391589*Ghat[16]*dx10; 
  out[43] += -0.7071067811865475*Ghat[19]*dx10; 
  out[44] += -1.58113883008419*Ghat[7]*dx10; 
  out[45] += -1.58113883008419*Ghat[8]*dx10; 
  out[46] += -0.7071067811865475*Ghat[20]*dx10; 
  out[47] += -1.58113883008419*Ghat[9]*dx10; 
  out[48] += -0.7071067811865475*Ghat[21]*dx10; 
  out[49] += -0.7071067811865475*Ghat[22]*dx10; 
  out[50] += -1.58113883008419*Ghat[10]*dx10; 
  out[51] += -1.224744871391589*Ghat[17]*dx10; 
  out[52] += -1.224744871391589*Ghat[18]*dx10; 
  out[53] += -1.224744871391589*Ghat[19]*dx10; 
  out[54] += -1.58113883008419*Ghat[11]*dx10; 
  out[55] += -1.58113883008419*Ghat[12]*dx10; 
  out[56] += -1.224744871391589*Ghat[20]*dx10; 
  out[57] += -1.58113883008419*Ghat[13]*dx10; 
  out[58] += -1.58113883008419*Ghat[14]*dx10; 
  out[59] += -0.7071067811865475*Ghat[23]*dx10; 
  out[60] += -1.58113883008419*Ghat[15]*dx10; 
  out[61] += -1.224744871391589*Ghat[21]*dx10; 
  out[62] += -1.58113883008419*Ghat[16]*dx10; 
  out[63] += -0.7071067811865475*Ghat[24]*dx10; 
  out[64] += -1.224744871391589*Ghat[22]*dx10; 
  out[65] += -0.7071067811865475*Ghat[25]*dx10; 
  out[66] += -1.58113883008419*Ghat[17]*dx10; 
  out[67] += -1.58113883008419*Ghat[18]*dx10; 
  out[68] += -1.224744871391589*Ghat[23]*dx10; 
  out[69] += -1.58113883008419*Ghat[19]*dx10; 
  out[70] += -1.224744871391589*Ghat[24]*dx10; 
  out[71] += -1.224744871391589*Ghat[25]*dx10; 
  out[72] += -1.58113883008419*Ghat[20]*dx10; 
  out[73] += -1.58113883008419*Ghat[21]*dx10; 
  out[74] += -1.58113883008419*Ghat[22]*dx10; 
  out[75] += -0.7071067811865475*Ghat[26]*dx10; 
  out[76] += -1.58113883008419*Ghat[23]*dx10; 
  out[77] += -1.58113883008419*Ghat[24]*dx10; 
  out[78] += -1.58113883008419*Ghat[25]*dx10; 
  out[79] += -1.224744871391589*Ghat[26]*dx10; 
  out[80] += -1.58113883008419*Ghat[26]*dx10; 

  } else { 

  if (wv>0) { 

  Ghat[0] = (1.58113883008419*fEdge[11]+1.224744871391589*fEdge[1]+0.7071067811865475*fEdge[0])*wv+(0.4564354645876384*fEdge[19]+0.3535533905932737*fEdge[5]+0.2041241452319315*fEdge[2])*dv; 
  Ghat[1] = (1.58113883008419*fEdge[19]+1.224744871391589*fEdge[5]+0.7071067811865475*fEdge[2])*wv+(0.408248290463863*fEdge[44]+0.3162277660168379*fEdge[20]+0.1825741858350554*fEdge[12]+0.4564354645876384*fEdge[11]+0.3535533905932737*fEdge[1]+0.2041241452319315*fEdge[0])*dv; 
  Ghat[2] = (1.58113883008419*fEdge[21]+1.224744871391589*fEdge[6]+0.7071067811865475*fEdge[3])*wv+(0.4564354645876384*fEdge[32]+0.3535533905932737*fEdge[15]+0.2041241452319315*fEdge[7])*dv; 
  Ghat[3] = (1.58113883008419*fEdge[25]+1.224744871391589*fEdge[8]+0.7071067811865475*fEdge[4])*wv+(0.4564354645876384*fEdge[35]+0.3535533905932737*fEdge[16]+0.2041241452319315*fEdge[9])*dv; 
  Ghat[4] = (1.58113883008419*fEdge[32]+1.224744871391589*fEdge[15]+0.7071067811865475*fEdge[7])*wv+(0.408248290463863*fEdge[54]+0.3162277660168379*fEdge[33]+0.1825741858350554*fEdge[22]+0.4564354645876384*fEdge[21]+0.3535533905932737*fEdge[6]+0.2041241452319315*fEdge[3])*dv; 
  Ghat[5] = (1.58113883008419*fEdge[35]+1.224744871391589*fEdge[16]+0.7071067811865475*fEdge[9])*wv+(0.408248290463863*fEdge[57]+0.3162277660168379*fEdge[36]+0.1825741858350554*fEdge[26]+0.4564354645876384*fEdge[25]+0.3535533905932737*fEdge[8]+0.2041241452319315*fEdge[4])*dv; 
  Ghat[6] = (1.58113883008419*fEdge[37]+1.224744871391589*fEdge[17]+0.7071067811865475*fEdge[10])*wv+(0.4564354645876384*fEdge[50]+0.3535533905932737*fEdge[31]+0.2041241452319315*fEdge[18])*dv; 
  Ghat[7] = (1.58113883008419*fEdge[44]+1.224744871391589*fEdge[20]+0.7071067811865475*fEdge[12])*wv+(0.408248290463863*fEdge[19]+0.3162277660168379*fEdge[5]+0.1825741858350554*fEdge[2])*dv; 
  Ghat[8] = (1.58113883008419*fEdge[45]+1.224744871391589*fEdge[23]+0.7071067811865475*fEdge[13])*wv+(0.4564354645876384*fEdge[55]+0.3535533905932737*fEdge[34]+0.2041241452319315*fEdge[24])*dv; 
  Ghat[9] = (1.58113883008419*fEdge[47]+1.224744871391589*fEdge[28]+0.7071067811865475*fEdge[14])*wv+(0.4564354645876384*fEdge[60]+0.3535533905932737*fEdge[41]+0.2041241452319315*fEdge[29])*dv; 
  Ghat[10] = (1.58113883008419*fEdge[50]+1.224744871391589*fEdge[31]+0.7071067811865475*fEdge[18])*wv+(0.408248290463863*fEdge[66]+0.3162277660168379*fEdge[51]+0.1825741858350554*fEdge[38]+0.4564354645876384*fEdge[37]+0.3535533905932737*fEdge[17]+0.2041241452319315*fEdge[10])*dv; 
  Ghat[11] = (1.58113883008419*fEdge[54]+1.224744871391589*fEdge[33]+0.7071067811865475*fEdge[22])*wv+(0.408248290463863*fEdge[32]+0.3162277660168379*fEdge[15]+0.1825741858350554*fEdge[7])*dv; 
  Ghat[12] = (1.58113883008419*fEdge[55]+1.224744871391589*fEdge[34]+0.7071067811865475*fEdge[24])*wv+(0.408248290463863*fEdge[72]+0.3162277660168379*fEdge[56]+0.1825741858350554*fEdge[46]+0.4564354645876384*fEdge[45]+0.3535533905932737*fEdge[23]+0.2041241452319315*fEdge[13])*dv; 
  Ghat[13] = (1.58113883008419*fEdge[57]+1.224744871391589*fEdge[36]+0.7071067811865475*fEdge[26])*wv+(0.408248290463863*fEdge[35]+0.3162277660168379*fEdge[16]+0.1825741858350554*fEdge[9])*dv; 
  Ghat[14] = (1.58113883008419*fEdge[58]+1.224744871391589*fEdge[39]+0.7071067811865475*fEdge[27])*wv+(0.4564354645876384*fEdge[67]+0.3535533905932737*fEdge[52]+0.2041241452319315*fEdge[40])*dv; 
  Ghat[15] = (1.58113883008419*fEdge[60]+1.224744871391589*fEdge[41]+0.7071067811865475*fEdge[29])*wv+(0.408248290463863*fEdge[73]+0.3162277660168379*fEdge[61]+0.1825741858350554*fEdge[48]+0.4564354645876384*fEdge[47]+0.3535533905932737*fEdge[28]+0.2041241452319315*fEdge[14])*dv; 
  Ghat[16] = (1.58113883008419*fEdge[62]+1.224744871391589*fEdge[42]+0.7071067811865475*fEdge[30])*wv+(0.4564354645876384*fEdge[69]+0.3535533905932737*fEdge[53]+0.2041241452319315*fEdge[43])*dv; 
  Ghat[17] = (1.58113883008419*fEdge[66]+1.224744871391589*fEdge[51]+0.7071067811865475*fEdge[38])*wv+(0.408248290463863*fEdge[50]+0.3162277660168379*fEdge[31]+0.1825741858350554*fEdge[18])*dv; 
  Ghat[18] = (1.58113883008419*fEdge[67]+1.224744871391589*fEdge[52]+0.7071067811865475*fEdge[40])*wv+(0.408248290463863*fEdge[76]+0.3162277660168379*fEdge[68]+0.1825741858350554*fEdge[59]+0.4564354645876384*fEdge[58]+0.3535533905932737*fEdge[39]+0.2041241452319315*fEdge[27])*dv; 
  Ghat[19] = (1.58113883008419*fEdge[69]+1.224744871391589*fEdge[53]+0.7071067811865475*fEdge[43])*wv+(0.408248290463863*fEdge[77]+0.3162277660168379*fEdge[70]+0.1825741858350554*fEdge[63]+0.4564354645876384*fEdge[62]+0.3535533905932737*fEdge[42]+0.2041241452319315*fEdge[30])*dv; 
  Ghat[20] = (1.58113883008419*fEdge[72]+1.224744871391589*fEdge[56]+0.7071067811865475*fEdge[46])*wv+(0.408248290463863*fEdge[55]+0.3162277660168379*fEdge[34]+0.1825741858350554*fEdge[24])*dv; 
  Ghat[21] = (1.58113883008419*fEdge[73]+1.224744871391589*fEdge[61]+0.7071067811865475*fEdge[48])*wv+(0.408248290463863*fEdge[60]+0.3162277660168379*fEdge[41]+0.1825741858350554*fEdge[29])*dv; 
  Ghat[22] = (1.58113883008419*fEdge[74]+1.224744871391589*fEdge[64]+0.7071067811865475*fEdge[49])*wv+(0.4564354645876384*fEdge[78]+0.3535533905932737*fEdge[71]+0.2041241452319315*fEdge[65])*dv; 
  Ghat[23] = (1.58113883008419*fEdge[76]+1.224744871391589*fEdge[68]+0.7071067811865475*fEdge[59])*wv+(0.408248290463863*fEdge[67]+0.3162277660168379*fEdge[52]+0.1825741858350554*fEdge[40])*dv; 
  Ghat[24] = (1.58113883008419*fEdge[77]+1.224744871391589*fEdge[70]+0.7071067811865475*fEdge[63])*wv+(0.408248290463863*fEdge[69]+0.3162277660168379*fEdge[53]+0.1825741858350554*fEdge[43])*dv; 
  Ghat[25] = (1.58113883008419*fEdge[78]+1.224744871391589*fEdge[71]+0.7071067811865475*fEdge[65])*wv+(0.408248290463863*fEdge[80]+0.3162277660168379*fEdge[79]+0.1825741858350554*fEdge[75]+0.4564354645876384*fEdge[74]+0.3535533905932737*fEdge[64]+0.2041241452319315*fEdge[49])*dv; 
  Ghat[26] = (1.58113883008419*fEdge[80]+1.224744871391589*fEdge[79]+0.7071067811865475*fEdge[75])*wv+(0.408248290463863*fEdge[78]+0.3162277660168379*fEdge[71]+0.1825741858350554*fEdge[65])*dv; 

  } else { 

  Ghat[0] = 1.58113883008419*fSkin[11]*wv-1.224744871391589*fSkin[1]*wv+0.7071067811865475*fSkin[0]*wv+0.4564354645876383*fSkin[19]*dv-0.3535533905932737*fSkin[5]*dv+0.2041241452319315*fSkin[2]*dv; 
  Ghat[1] = 1.58113883008419*fSkin[19]*wv-1.224744871391589*fSkin[5]*wv+0.7071067811865475*fSkin[2]*wv+0.408248290463863*fSkin[44]*dv-0.3162277660168379*fSkin[20]*dv+0.1825741858350554*fSkin[12]*dv+0.4564354645876384*fSkin[11]*dv-0.3535533905932737*fSkin[1]*dv+0.2041241452319315*fSkin[0]*dv; 
  Ghat[2] = 1.58113883008419*fSkin[21]*wv-1.224744871391589*fSkin[6]*wv+0.7071067811865475*fSkin[3]*wv+0.4564354645876384*fSkin[32]*dv-0.3535533905932737*fSkin[15]*dv+0.2041241452319315*fSkin[7]*dv; 
  Ghat[3] = 1.58113883008419*fSkin[25]*wv-1.224744871391589*fSkin[8]*wv+0.7071067811865475*fSkin[4]*wv+0.4564354645876384*fSkin[35]*dv-0.3535533905932737*fSkin[16]*dv+0.2041241452319315*fSkin[9]*dv; 
  Ghat[4] = 1.58113883008419*fSkin[32]*wv-1.224744871391589*fSkin[15]*wv+0.7071067811865475*fSkin[7]*wv+0.408248290463863*fSkin[54]*dv-0.3162277660168379*fSkin[33]*dv+0.1825741858350553*fSkin[22]*dv+0.4564354645876383*fSkin[21]*dv-0.3535533905932737*fSkin[6]*dv+0.2041241452319315*fSkin[3]*dv; 
  Ghat[5] = 1.58113883008419*fSkin[35]*wv-1.224744871391589*fSkin[16]*wv+0.7071067811865475*fSkin[9]*wv+0.408248290463863*fSkin[57]*dv-0.3162277660168379*fSkin[36]*dv+0.1825741858350553*fSkin[26]*dv+0.4564354645876383*fSkin[25]*dv-0.3535533905932737*fSkin[8]*dv+0.2041241452319315*fSkin[4]*dv; 
  Ghat[6] = 1.58113883008419*fSkin[37]*wv-1.224744871391589*fSkin[17]*wv+0.7071067811865475*fSkin[10]*wv+0.4564354645876383*fSkin[50]*dv-0.3535533905932737*fSkin[31]*dv+0.2041241452319315*fSkin[18]*dv; 
  Ghat[7] = 1.58113883008419*fSkin[44]*wv-1.224744871391589*fSkin[20]*wv+0.7071067811865475*fSkin[12]*wv+0.408248290463863*fSkin[19]*dv-0.3162277660168379*fSkin[5]*dv+0.1825741858350554*fSkin[2]*dv; 
  Ghat[8] = 1.58113883008419*fSkin[45]*wv-1.224744871391589*fSkin[23]*wv+0.7071067811865475*fSkin[13]*wv+0.4564354645876384*fSkin[55]*dv-0.3535533905932737*fSkin[34]*dv+0.2041241452319315*fSkin[24]*dv; 
  Ghat[9] = 1.58113883008419*fSkin[47]*wv-1.224744871391589*fSkin[28]*wv+0.7071067811865475*fSkin[14]*wv+0.4564354645876384*fSkin[60]*dv-0.3535533905932737*fSkin[41]*dv+0.2041241452319315*fSkin[29]*dv; 
  Ghat[10] = 1.58113883008419*fSkin[50]*wv-1.224744871391589*fSkin[31]*wv+0.7071067811865475*fSkin[18]*wv+0.408248290463863*fSkin[66]*dv-0.3162277660168379*fSkin[51]*dv+0.1825741858350554*fSkin[38]*dv+0.4564354645876384*fSkin[37]*dv-0.3535533905932737*fSkin[17]*dv+0.2041241452319315*fSkin[10]*dv; 
  Ghat[11] = 1.58113883008419*fSkin[54]*wv-1.224744871391589*fSkin[33]*wv+0.7071067811865475*fSkin[22]*wv+0.408248290463863*fSkin[32]*dv-0.3162277660168379*fSkin[15]*dv+0.1825741858350553*fSkin[7]*dv; 
  Ghat[12] = 1.58113883008419*fSkin[55]*wv-1.224744871391589*fSkin[34]*wv+0.7071067811865475*fSkin[24]*wv+0.408248290463863*fSkin[72]*dv-0.3162277660168379*fSkin[56]*dv+0.1825741858350553*fSkin[46]*dv+0.4564354645876383*fSkin[45]*dv-0.3535533905932737*fSkin[23]*dv+0.2041241452319315*fSkin[13]*dv; 
  Ghat[13] = 1.58113883008419*fSkin[57]*wv-1.224744871391589*fSkin[36]*wv+0.7071067811865475*fSkin[26]*wv+0.408248290463863*fSkin[35]*dv-0.3162277660168379*fSkin[16]*dv+0.1825741858350553*fSkin[9]*dv; 
  Ghat[14] = 1.58113883008419*fSkin[58]*wv-1.224744871391589*fSkin[39]*wv+0.7071067811865475*fSkin[27]*wv+0.4564354645876383*fSkin[67]*dv-0.3535533905932737*fSkin[52]*dv+0.2041241452319315*fSkin[40]*dv; 
  Ghat[15] = 1.58113883008419*fSkin[60]*wv-1.224744871391589*fSkin[41]*wv+0.7071067811865475*fSkin[29]*wv+0.408248290463863*fSkin[73]*dv-0.3162277660168379*fSkin[61]*dv+0.1825741858350553*fSkin[48]*dv+0.4564354645876383*fSkin[47]*dv-0.3535533905932737*fSkin[28]*dv+0.2041241452319315*fSkin[14]*dv; 
  Ghat[16] = 1.58113883008419*fSkin[62]*wv-1.224744871391589*fSkin[42]*wv+0.7071067811865475*fSkin[30]*wv+0.4564354645876383*fSkin[69]*dv-0.3535533905932737*fSkin[53]*dv+0.2041241452319315*fSkin[43]*dv; 
  Ghat[17] = 1.58113883008419*fSkin[66]*wv-1.224744871391589*fSkin[51]*wv+0.7071067811865475*fSkin[38]*wv+0.408248290463863*fSkin[50]*dv-0.3162277660168379*fSkin[31]*dv+0.1825741858350554*fSkin[18]*dv; 
  Ghat[18] = 1.58113883008419*fSkin[67]*wv-1.224744871391589*fSkin[52]*wv+0.7071067811865475*fSkin[40]*wv+0.408248290463863*fSkin[76]*dv-0.3162277660168379*fSkin[68]*dv+0.1825741858350554*fSkin[59]*dv+0.4564354645876384*fSkin[58]*dv-0.3535533905932737*fSkin[39]*dv+0.2041241452319315*fSkin[27]*dv; 
  Ghat[19] = 1.58113883008419*fSkin[69]*wv-1.224744871391589*fSkin[53]*wv+0.7071067811865475*fSkin[43]*wv+0.408248290463863*fSkin[77]*dv-0.3162277660168379*fSkin[70]*dv+0.1825741858350554*fSkin[63]*dv+0.4564354645876384*fSkin[62]*dv-0.3535533905932737*fSkin[42]*dv+0.2041241452319315*fSkin[30]*dv; 
  Ghat[20] = 1.58113883008419*fSkin[72]*wv-1.224744871391589*fSkin[56]*wv+0.7071067811865475*fSkin[46]*wv+0.408248290463863*fSkin[55]*dv-0.3162277660168379*fSkin[34]*dv+0.1825741858350553*fSkin[24]*dv; 
  Ghat[21] = 1.58113883008419*fSkin[73]*wv-1.224744871391589*fSkin[61]*wv+0.7071067811865475*fSkin[48]*wv+0.408248290463863*fSkin[60]*dv-0.3162277660168379*fSkin[41]*dv+0.1825741858350553*fSkin[29]*dv; 
  Ghat[22] = 1.58113883008419*fSkin[74]*wv-1.224744871391589*fSkin[64]*wv+0.7071067811865475*fSkin[49]*wv+0.4564354645876383*fSkin[78]*dv-0.3535533905932737*fSkin[71]*dv+0.2041241452319315*fSkin[65]*dv; 
  Ghat[23] = 1.58113883008419*fSkin[76]*wv-1.224744871391589*fSkin[68]*wv+0.7071067811865475*fSkin[59]*wv+0.408248290463863*fSkin[67]*dv-0.3162277660168379*fSkin[52]*dv+0.1825741858350554*fSkin[40]*dv; 
  Ghat[24] = 1.58113883008419*fSkin[77]*wv-1.224744871391589*fSkin[70]*wv+0.7071067811865475*fSkin[63]*wv+0.408248290463863*fSkin[69]*dv-0.3162277660168379*fSkin[53]*dv+0.1825741858350554*fSkin[43]*dv; 
  Ghat[25] = 1.58113883008419*fSkin[78]*wv-1.224744871391589*fSkin[71]*wv+0.7071067811865475*fSkin[65]*wv+0.408248290463863*fSkin[80]*dv-0.3162277660168379*fSkin[79]*dv+0.1825741858350554*fSkin[75]*dv+0.4564354645876384*fSkin[74]*dv-0.3535533905932737*fSkin[64]*dv+0.2041241452319315*fSkin[49]*dv; 
  Ghat[26] = 1.58113883008419*fSkin[80]*wv-1.224744871391589*fSkin[79]*wv+0.7071067811865475*fSkin[75]*wv+0.408248290463863*fSkin[78]*dv-0.3162277660168379*fSkin[71]*dv+0.1825741858350554*fSkin[65]*dv; 

  } 

  out[0] += 0.7071067811865475*Ghat[0]*dx10; 
  out[1] += -1.224744871391589*Ghat[0]*dx10; 
  out[2] += 0.7071067811865475*Ghat[1]*dx10; 
  out[3] += 0.7071067811865475*Ghat[2]*dx10; 
  out[4] += 0.7071067811865475*Ghat[3]*dx10; 
  out[5] += -1.224744871391589*Ghat[1]*dx10; 
  out[6] += -1.224744871391589*Ghat[2]*dx10; 
  out[7] += 0.7071067811865475*Ghat[4]*dx10; 
  out[8] += -1.224744871391589*Ghat[3]*dx10; 
  out[9] += 0.7071067811865475*Ghat[5]*dx10; 
  out[10] += 0.7071067811865475*Ghat[6]*dx10; 
  out[11] += 1.58113883008419*Ghat[0]*dx10; 
  out[12] += 0.7071067811865475*Ghat[7]*dx10; 
  out[13] += 0.7071067811865475*Ghat[8]*dx10; 
  out[14] += 0.7071067811865475*Ghat[9]*dx10; 
  out[15] += -1.224744871391589*Ghat[4]*dx10; 
  out[16] += -1.224744871391589*Ghat[5]*dx10; 
  out[17] += -1.224744871391589*Ghat[6]*dx10; 
  out[18] += 0.7071067811865475*Ghat[10]*dx10; 
  out[19] += 1.58113883008419*Ghat[1]*dx10; 
  out[20] += -1.224744871391589*Ghat[7]*dx10; 
  out[21] += 1.58113883008419*Ghat[2]*dx10; 
  out[22] += 0.7071067811865475*Ghat[11]*dx10; 
  out[23] += -1.224744871391589*Ghat[8]*dx10; 
  out[24] += 0.7071067811865475*Ghat[12]*dx10; 
  out[25] += 1.58113883008419*Ghat[3]*dx10; 
  out[26] += 0.7071067811865475*Ghat[13]*dx10; 
  out[27] += 0.7071067811865475*Ghat[14]*dx10; 
  out[28] += -1.224744871391589*Ghat[9]*dx10; 
  out[29] += 0.7071067811865475*Ghat[15]*dx10; 
  out[30] += 0.7071067811865475*Ghat[16]*dx10; 
  out[31] += -1.224744871391589*Ghat[10]*dx10; 
  out[32] += 1.58113883008419*Ghat[4]*dx10; 
  out[33] += -1.224744871391589*Ghat[11]*dx10; 
  out[34] += -1.224744871391589*Ghat[12]*dx10; 
  out[35] += 1.58113883008419*Ghat[5]*dx10; 
  out[36] += -1.224744871391589*Ghat[13]*dx10; 
  out[37] += 1.58113883008419*Ghat[6]*dx10; 
  out[38] += 0.7071067811865475*Ghat[17]*dx10; 
  out[39] += -1.224744871391589*Ghat[14]*dx10; 
  out[40] += 0.7071067811865475*Ghat[18]*dx10; 
  out[41] += -1.224744871391589*Ghat[15]*dx10; 
  out[42] += -1.224744871391589*Ghat[16]*dx10; 
  out[43] += 0.7071067811865475*Ghat[19]*dx10; 
  out[44] += 1.58113883008419*Ghat[7]*dx10; 
  out[45] += 1.58113883008419*Ghat[8]*dx10; 
  out[46] += 0.7071067811865475*Ghat[20]*dx10; 
  out[47] += 1.58113883008419*Ghat[9]*dx10; 
  out[48] += 0.7071067811865475*Ghat[21]*dx10; 
  out[49] += 0.7071067811865475*Ghat[22]*dx10; 
  out[50] += 1.58113883008419*Ghat[10]*dx10; 
  out[51] += -1.224744871391589*Ghat[17]*dx10; 
  out[52] += -1.224744871391589*Ghat[18]*dx10; 
  out[53] += -1.224744871391589*Ghat[19]*dx10; 
  out[54] += 1.58113883008419*Ghat[11]*dx10; 
  out[55] += 1.58113883008419*Ghat[12]*dx10; 
  out[56] += -1.224744871391589*Ghat[20]*dx10; 
  out[57] += 1.58113883008419*Ghat[13]*dx10; 
  out[58] += 1.58113883008419*Ghat[14]*dx10; 
  out[59] += 0.7071067811865475*Ghat[23]*dx10; 
  out[60] += 1.58113883008419*Ghat[15]*dx10; 
  out[61] += -1.224744871391589*Ghat[21]*dx10; 
  out[62] += 1.58113883008419*Ghat[16]*dx10; 
  out[63] += 0.7071067811865475*Ghat[24]*dx10; 
  out[64] += -1.224744871391589*Ghat[22]*dx10; 
  out[65] += 0.7071067811865475*Ghat[25]*dx10; 
  out[66] += 1.58113883008419*Ghat[17]*dx10; 
  out[67] += 1.58113883008419*Ghat[18]*dx10; 
  out[68] += -1.224744871391589*Ghat[23]*dx10; 
  out[69] += 1.58113883008419*Ghat[19]*dx10; 
  out[70] += -1.224744871391589*Ghat[24]*dx10; 
  out[71] += -1.224744871391589*Ghat[25]*dx10; 
  out[72] += 1.58113883008419*Ghat[20]*dx10; 
  out[73] += 1.58113883008419*Ghat[21]*dx10; 
  out[74] += 1.58113883008419*Ghat[22]*dx10; 
  out[75] += 0.7071067811865475*Ghat[26]*dx10; 
  out[76] += 1.58113883008419*Ghat[23]*dx10; 
  out[77] += 1.58113883008419*Ghat[24]*dx10; 
  out[78] += 1.58113883008419*Ghat[25]*dx10; 
  out[79] += -1.224744871391589*Ghat[26]*dx10; 
  out[80] += 1.58113883008419*Ghat[26]*dx10; 

  } 
  return 0.;

} 