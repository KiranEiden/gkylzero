#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_boundary_surfx_3x3v_tensor_p1(const double *w, const double *dxv, const double *alpha_geo, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:     Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // alpha_geo:   Fields used only for general geometry.
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Incremented distribution function in center cell.
  const double dx10 = 2/dxv[0]; 
  const double dv = dxv[3], wv = w[3]; 
  double Ghat[80]; 

  if (edge == -1) { 

  if (wv>0) { 

  Ghat[0] = (1.224744871391589*fSkin[1]+0.7071067811865475*fSkin[0])*wv+(0.3535533905932737*fSkin[10]+0.2041241452319315*fSkin[4])*dv; 
  Ghat[1] = (1.224744871391589*fSkin[7]+0.7071067811865475*fSkin[2])*wv+(0.3535533905932737*fSkin[23]+0.2041241452319315*fSkin[11])*dv; 
  Ghat[2] = (1.224744871391589*fSkin[8]+0.7071067811865475*fSkin[3])*wv+(0.3535533905932737*fSkin[24]+0.2041241452319315*fSkin[12])*dv; 
  Ghat[3] = (1.224744871391589*fSkin[10]+0.7071067811865475*fSkin[4])*wv+(0.3535533905932737*fSkin[1]+0.2041241452319315*fSkin[0])*dv; 
  Ghat[4] = (1.224744871391589*fSkin[13]+0.7071067811865475*fSkin[5])*wv+(0.3535533905932737*fSkin[29]+0.2041241452319315*fSkin[16])*dv; 
  Ghat[5] = (1.224744871391589*fSkin[17]+0.7071067811865475*fSkin[6])*wv+(0.3535533905932737*fSkin[35]+0.2041241452319315*fSkin[20])*dv; 
  Ghat[6] = (1.224744871391589*fSkin[22]+0.7071067811865475*fSkin[9])*wv+(0.3535533905932737*fSkin[42]+0.2041241452319315*fSkin[25])*dv; 
  Ghat[7] = (1.224744871391589*fSkin[23]+0.7071067811865475*fSkin[11])*wv+(0.3535533905932737*fSkin[7]+0.2041241452319315*fSkin[2])*dv; 
  Ghat[8] = (1.224744871391589*fSkin[24]+0.7071067811865475*fSkin[12])*wv+(0.3535533905932737*fSkin[8]+0.2041241452319315*fSkin[3])*dv; 
  Ghat[9] = (1.224744871391589*fSkin[26]+0.7071067811865475*fSkin[14])*wv+(0.3535533905932737*fSkin[44]+0.2041241452319315*fSkin[30])*dv; 
  Ghat[10] = (1.224744871391589*fSkin[27]+0.7071067811865475*fSkin[15])*wv+(0.3535533905932737*fSkin[45]+0.2041241452319315*fSkin[31])*dv; 
  Ghat[11] = (1.224744871391589*fSkin[29]+0.7071067811865475*fSkin[16])*wv+(0.3535533905932737*fSkin[13]+0.2041241452319315*fSkin[5])*dv; 
  Ghat[12] = (1.224744871391589*fSkin[32]+0.7071067811865475*fSkin[18])*wv+(0.3535533905932737*fSkin[48]+0.2041241452319315*fSkin[36])*dv; 
  Ghat[13] = (1.224744871391589*fSkin[33]+0.7071067811865475*fSkin[19])*wv+(0.3535533905932737*fSkin[49]+0.2041241452319315*fSkin[37])*dv; 
  Ghat[14] = (1.224744871391589*fSkin[35]+0.7071067811865475*fSkin[20])*wv+(0.3535533905932737*fSkin[17]+0.2041241452319315*fSkin[6])*dv; 
  Ghat[15] = (1.224744871391589*fSkin[38]+0.7071067811865475*fSkin[21])*wv+(0.3535533905932737*fSkin[54]+0.2041241452319315*fSkin[41])*dv; 
  Ghat[16] = (1.224744871391589*fSkin[42]+0.7071067811865475*fSkin[25])*wv+(0.3535533905932737*fSkin[22]+0.2041241452319315*fSkin[9])*dv; 
  Ghat[17] = (1.224744871391589*fSkin[43]+0.7071067811865475*fSkin[28])*wv+(0.3535533905932737*fSkin[57]+0.2041241452319315*fSkin[46])*dv; 
  Ghat[18] = (1.224744871391589*fSkin[44]+0.7071067811865475*fSkin[30])*wv+(0.3535533905932737*fSkin[26]+0.2041241452319315*fSkin[14])*dv; 
  Ghat[19] = (1.224744871391589*fSkin[45]+0.7071067811865475*fSkin[31])*wv+(0.3535533905932737*fSkin[27]+0.2041241452319315*fSkin[15])*dv; 
  Ghat[20] = (1.224744871391589*fSkin[47]+0.7071067811865475*fSkin[34])*wv+(0.3535533905932737*fSkin[58]+0.2041241452319315*fSkin[50])*dv; 
  Ghat[21] = (1.224744871391589*fSkin[48]+0.7071067811865475*fSkin[36])*wv+(0.3535533905932737*fSkin[32]+0.2041241452319315*fSkin[18])*dv; 
  Ghat[22] = (1.224744871391589*fSkin[49]+0.7071067811865475*fSkin[37])*wv+(0.3535533905932737*fSkin[33]+0.2041241452319315*fSkin[19])*dv; 
  Ghat[23] = (1.224744871391589*fSkin[51]+0.7071067811865475*fSkin[39])*wv+(0.3535533905932737*fSkin[60]+0.2041241452319315*fSkin[55])*dv; 
  Ghat[24] = (1.224744871391589*fSkin[52]+0.7071067811865475*fSkin[40])*wv+(0.3535533905932737*fSkin[61]+0.2041241452319315*fSkin[56])*dv; 
  Ghat[25] = (1.224744871391589*fSkin[54]+0.7071067811865475*fSkin[41])*wv+(0.3535533905932737*fSkin[38]+0.2041241452319315*fSkin[21])*dv; 
  Ghat[26] = (1.224744871391589*fSkin[57]+0.7071067811865475*fSkin[46])*wv+(0.3535533905932737*fSkin[43]+0.2041241452319315*fSkin[28])*dv; 
  Ghat[27] = (1.224744871391589*fSkin[58]+0.7071067811865475*fSkin[50])*wv+(0.3535533905932737*fSkin[47]+0.2041241452319315*fSkin[34])*dv; 
  Ghat[28] = (1.224744871391589*fSkin[59]+0.7071067811865475*fSkin[53])*wv+(0.3535533905932737*fSkin[63]+0.2041241452319315*fSkin[62])*dv; 
  Ghat[29] = (1.224744871391589*fSkin[60]+0.7071067811865475*fSkin[55])*wv+(0.3535533905932737*fSkin[51]+0.2041241452319315*fSkin[39])*dv; 
  Ghat[30] = (1.224744871391589*fSkin[61]+0.7071067811865475*fSkin[56])*wv+(0.3535533905932737*fSkin[52]+0.2041241452319315*fSkin[40])*dv; 
  Ghat[31] = (1.224744871391589*fSkin[63]+0.7071067811865475*fSkin[62])*wv+(0.3535533905932737*fSkin[59]+0.2041241452319315*fSkin[53])*dv; 
  Ghat[32] = (0.3162277660168379*fSkin[10]+0.1825741858350554*fSkin[4])*dv; 
  Ghat[33] = (0.3162277660168379*fSkin[23]+0.1825741858350554*fSkin[11])*dv; 
  Ghat[34] = (0.3162277660168379*fSkin[24]+0.1825741858350554*fSkin[12])*dv; 
  Ghat[35] = (0.3162277660168379*fSkin[29]+0.1825741858350554*fSkin[16])*dv; 
  Ghat[36] = (0.3162277660168379*fSkin[35]+0.1825741858350554*fSkin[20])*dv; 
  Ghat[37] = (0.3162277660168379*fSkin[42]+0.1825741858350554*fSkin[25])*dv; 
  Ghat[38] = (0.3162277660168379*fSkin[44]+0.1825741858350554*fSkin[30])*dv; 
  Ghat[39] = (0.3162277660168379*fSkin[45]+0.1825741858350554*fSkin[31])*dv; 
  Ghat[40] = (0.3162277660168379*fSkin[48]+0.1825741858350554*fSkin[36])*dv; 
  Ghat[41] = (0.3162277660168379*fSkin[49]+0.1825741858350554*fSkin[37])*dv; 
  Ghat[42] = (0.3162277660168379*fSkin[54]+0.1825741858350554*fSkin[41])*dv; 
  Ghat[43] = (0.3162277660168379*fSkin[57]+0.1825741858350554*fSkin[46])*dv; 
  Ghat[44] = (0.3162277660168379*fSkin[58]+0.1825741858350554*fSkin[50])*dv; 
  Ghat[45] = (0.3162277660168379*fSkin[60]+0.1825741858350554*fSkin[55])*dv; 
  Ghat[46] = (0.3162277660168379*fSkin[61]+0.1825741858350554*fSkin[56])*dv; 
  Ghat[47] = (0.3162277660168379*fSkin[63]+0.1825741858350554*fSkin[62])*dv; 

  } else { 

  Ghat[0] = -0.1178511301977579*((10.39230484541326*fEdge[1]-6.0*fEdge[0])*wv+(3.0*fEdge[10]-1.732050807568877*fEdge[4])*dv); 
  Ghat[1] = -0.1178511301977579*((10.39230484541326*fEdge[7]-6.0*fEdge[2])*wv+(3.0*fEdge[23]-1.732050807568877*fEdge[11])*dv); 
  Ghat[2] = -0.1178511301977579*((10.39230484541326*fEdge[8]-6.0*fEdge[3])*wv+(3.0*fEdge[24]-1.732050807568877*fEdge[12])*dv); 
  Ghat[3] = -0.1178511301977579*((10.39230484541326*fEdge[10]-6.0*fEdge[4])*wv+(3.0*fEdge[1]-1.732050807568877*fEdge[0])*dv); 
  Ghat[4] = -0.1178511301977579*((10.39230484541326*fEdge[13]-6.0*fEdge[5])*wv+(3.0*fEdge[29]-1.732050807568877*fEdge[16])*dv); 
  Ghat[5] = -0.1178511301977579*((10.39230484541326*fEdge[17]-6.0*fEdge[6])*wv+(3.0*fEdge[35]-1.732050807568877*fEdge[20])*dv); 
  Ghat[6] = -0.1178511301977579*((10.39230484541326*fEdge[22]-6.0*fEdge[9])*wv+(3.0*fEdge[42]-1.732050807568877*fEdge[25])*dv); 
  Ghat[7] = -0.1178511301977579*((10.39230484541326*fEdge[23]-6.0*fEdge[11])*wv+(3.0*fEdge[7]-1.732050807568877*fEdge[2])*dv); 
  Ghat[8] = -0.1178511301977579*((10.39230484541326*fEdge[24]-6.0*fEdge[12])*wv+(3.0*fEdge[8]-1.732050807568877*fEdge[3])*dv); 
  Ghat[9] = -0.1178511301977579*((10.39230484541326*fEdge[26]-6.0*fEdge[14])*wv+(3.0*fEdge[44]-1.732050807568877*fEdge[30])*dv); 
  Ghat[10] = -0.1178511301977579*((10.39230484541326*fEdge[27]-6.0*fEdge[15])*wv+(3.0*fEdge[45]-1.732050807568877*fEdge[31])*dv); 
  Ghat[11] = -0.1178511301977579*((10.39230484541326*fEdge[29]-6.0*fEdge[16])*wv+(3.0*fEdge[13]-1.732050807568877*fEdge[5])*dv); 
  Ghat[12] = -0.1178511301977579*((10.39230484541326*fEdge[32]-6.0*fEdge[18])*wv+(3.0*fEdge[48]-1.732050807568877*fEdge[36])*dv); 
  Ghat[13] = -0.1178511301977579*((10.39230484541326*fEdge[33]-6.0*fEdge[19])*wv+(3.0*fEdge[49]-1.732050807568877*fEdge[37])*dv); 
  Ghat[14] = -0.1178511301977579*((10.39230484541326*fEdge[35]-6.0*fEdge[20])*wv+(3.0*fEdge[17]-1.732050807568877*fEdge[6])*dv); 
  Ghat[15] = -0.1178511301977579*((10.39230484541326*fEdge[38]-6.0*fEdge[21])*wv+(3.0*fEdge[54]-1.732050807568877*fEdge[41])*dv); 
  Ghat[16] = -0.1178511301977579*((10.39230484541326*fEdge[42]-6.0*fEdge[25])*wv+(3.0*fEdge[22]-1.732050807568877*fEdge[9])*dv); 
  Ghat[17] = -0.1178511301977579*((10.39230484541326*fEdge[43]-6.0*fEdge[28])*wv+(3.0*fEdge[57]-1.732050807568877*fEdge[46])*dv); 
  Ghat[18] = -0.1178511301977579*((10.39230484541326*fEdge[44]-6.0*fEdge[30])*wv+(3.0*fEdge[26]-1.732050807568877*fEdge[14])*dv); 
  Ghat[19] = -0.1178511301977579*((10.39230484541326*fEdge[45]-6.0*fEdge[31])*wv+(3.0*fEdge[27]-1.732050807568877*fEdge[15])*dv); 
  Ghat[20] = -0.1178511301977579*((10.39230484541326*fEdge[47]-6.0*fEdge[34])*wv+(3.0*fEdge[58]-1.732050807568877*fEdge[50])*dv); 
  Ghat[21] = -0.1178511301977579*((10.39230484541326*fEdge[48]-6.0*fEdge[36])*wv+(3.0*fEdge[32]-1.732050807568877*fEdge[18])*dv); 
  Ghat[22] = -0.1178511301977579*((10.39230484541326*fEdge[49]-6.0*fEdge[37])*wv+(3.0*fEdge[33]-1.732050807568877*fEdge[19])*dv); 
  Ghat[23] = -0.1178511301977579*((10.39230484541326*fEdge[51]-6.0*fEdge[39])*wv+(3.0*fEdge[60]-1.732050807568877*fEdge[55])*dv); 
  Ghat[24] = -0.1178511301977579*((10.39230484541326*fEdge[52]-6.0*fEdge[40])*wv+(3.0*fEdge[61]-1.732050807568877*fEdge[56])*dv); 
  Ghat[25] = -0.1178511301977579*((10.39230484541326*fEdge[54]-6.0*fEdge[41])*wv+(3.0*fEdge[38]-1.732050807568877*fEdge[21])*dv); 
  Ghat[26] = -0.1178511301977579*((10.39230484541326*fEdge[57]-6.0*fEdge[46])*wv+(3.0*fEdge[43]-1.732050807568877*fEdge[28])*dv); 
  Ghat[27] = -0.1178511301977579*((10.39230484541326*fEdge[58]-6.0*fEdge[50])*wv+(3.0*fEdge[47]-1.732050807568877*fEdge[34])*dv); 
  Ghat[28] = -0.1178511301977579*((10.39230484541326*fEdge[59]-6.0*fEdge[53])*wv+(3.0*fEdge[63]-1.732050807568877*fEdge[62])*dv); 
  Ghat[29] = -0.1178511301977579*((10.39230484541326*fEdge[60]-6.0*fEdge[55])*wv+(3.0*fEdge[51]-1.732050807568877*fEdge[39])*dv); 
  Ghat[30] = -0.1178511301977579*((10.39230484541326*fEdge[61]-6.0*fEdge[56])*wv+(3.0*fEdge[52]-1.732050807568877*fEdge[40])*dv); 
  Ghat[31] = -0.1178511301977579*((10.39230484541326*fEdge[63]-6.0*fEdge[62])*wv+(3.0*fEdge[59]-1.732050807568877*fEdge[53])*dv); 
  Ghat[32] = -0.04714045207910316*(6.708203932499369*fEdge[10]-3.872983346207417*fEdge[4])*dv; 
  Ghat[33] = -0.04714045207910316*(6.708203932499369*fEdge[23]-3.872983346207417*fEdge[11])*dv; 
  Ghat[34] = -0.04714045207910316*(6.708203932499369*fEdge[24]-3.872983346207417*fEdge[12])*dv; 
  Ghat[35] = -0.04714045207910316*(6.708203932499369*fEdge[29]-3.872983346207417*fEdge[16])*dv; 
  Ghat[36] = -0.04714045207910316*(6.708203932499369*fEdge[35]-3.872983346207417*fEdge[20])*dv; 
  Ghat[37] = -0.04714045207910316*(6.708203932499369*fEdge[42]-3.872983346207417*fEdge[25])*dv; 
  Ghat[38] = -0.04714045207910316*(6.708203932499369*fEdge[44]-3.872983346207417*fEdge[30])*dv; 
  Ghat[39] = -0.04714045207910316*(6.708203932499369*fEdge[45]-3.872983346207417*fEdge[31])*dv; 
  Ghat[40] = -0.04714045207910316*(6.708203932499369*fEdge[48]-3.872983346207417*fEdge[36])*dv; 
  Ghat[41] = -0.04714045207910316*(6.708203932499369*fEdge[49]-3.872983346207417*fEdge[37])*dv; 
  Ghat[42] = -0.04714045207910316*(6.708203932499369*fEdge[54]-3.872983346207417*fEdge[41])*dv; 
  Ghat[43] = -0.04714045207910316*(6.708203932499369*fEdge[57]-3.872983346207417*fEdge[46])*dv; 
  Ghat[44] = -0.04714045207910316*(6.708203932499369*fEdge[58]-3.872983346207417*fEdge[50])*dv; 
  Ghat[45] = -0.04714045207910316*(6.708203932499369*fEdge[60]-3.872983346207417*fEdge[55])*dv; 
  Ghat[46] = -0.04714045207910316*(6.708203932499369*fEdge[61]-3.872983346207417*fEdge[56])*dv; 
  Ghat[47] = -0.04714045207910316*(6.708203932499369*fEdge[63]-3.872983346207417*fEdge[62])*dv; 

  } 

  out[0] += -0.7071067811865475*Ghat[0]*dx10; 
  out[1] += -1.224744871391589*Ghat[0]*dx10; 
  out[2] += -0.7071067811865475*Ghat[1]*dx10; 
  out[3] += -0.7071067811865475*Ghat[2]*dx10; 
  out[4] += -0.7071067811865475*Ghat[3]*dx10; 
  out[5] += -0.7071067811865475*Ghat[4]*dx10; 
  out[6] += -0.7071067811865475*Ghat[5]*dx10; 
  out[7] += -1.224744871391589*Ghat[1]*dx10; 
  out[8] += -1.224744871391589*Ghat[2]*dx10; 
  out[9] += -0.7071067811865475*Ghat[6]*dx10; 
  out[10] += -1.224744871391589*Ghat[3]*dx10; 
  out[11] += -0.7071067811865475*Ghat[7]*dx10; 
  out[12] += -0.7071067811865475*Ghat[8]*dx10; 
  out[13] += -1.224744871391589*Ghat[4]*dx10; 
  out[14] += -0.7071067811865475*Ghat[9]*dx10; 
  out[15] += -0.7071067811865475*Ghat[10]*dx10; 
  out[16] += -0.7071067811865475*Ghat[11]*dx10; 
  out[17] += -1.224744871391589*Ghat[5]*dx10; 
  out[18] += -0.7071067811865475*Ghat[12]*dx10; 
  out[19] += -0.7071067811865475*Ghat[13]*dx10; 
  out[20] += -0.7071067811865475*Ghat[14]*dx10; 
  out[21] += -0.7071067811865475*Ghat[15]*dx10; 
  out[22] += -1.224744871391589*Ghat[6]*dx10; 
  out[23] += -1.224744871391589*Ghat[7]*dx10; 
  out[24] += -1.224744871391589*Ghat[8]*dx10; 
  out[25] += -0.7071067811865475*Ghat[16]*dx10; 
  out[26] += -1.224744871391589*Ghat[9]*dx10; 
  out[27] += -1.224744871391589*Ghat[10]*dx10; 
  out[28] += -0.7071067811865475*Ghat[17]*dx10; 
  out[29] += -1.224744871391589*Ghat[11]*dx10; 
  out[30] += -0.7071067811865475*Ghat[18]*dx10; 
  out[31] += -0.7071067811865475*Ghat[19]*dx10; 
  out[32] += -1.224744871391589*Ghat[12]*dx10; 
  out[33] += -1.224744871391589*Ghat[13]*dx10; 
  out[34] += -0.7071067811865475*Ghat[20]*dx10; 
  out[35] += -1.224744871391589*Ghat[14]*dx10; 
  out[36] += -0.7071067811865475*Ghat[21]*dx10; 
  out[37] += -0.7071067811865475*Ghat[22]*dx10; 
  out[38] += -1.224744871391589*Ghat[15]*dx10; 
  out[39] += -0.7071067811865475*Ghat[23]*dx10; 
  out[40] += -0.7071067811865475*Ghat[24]*dx10; 
  out[41] += -0.7071067811865475*Ghat[25]*dx10; 
  out[42] += -1.224744871391589*Ghat[16]*dx10; 
  out[43] += -1.224744871391589*Ghat[17]*dx10; 
  out[44] += -1.224744871391589*Ghat[18]*dx10; 
  out[45] += -1.224744871391589*Ghat[19]*dx10; 
  out[46] += -0.7071067811865475*Ghat[26]*dx10; 
  out[47] += -1.224744871391589*Ghat[20]*dx10; 
  out[48] += -1.224744871391589*Ghat[21]*dx10; 
  out[49] += -1.224744871391589*Ghat[22]*dx10; 
  out[50] += -0.7071067811865475*Ghat[27]*dx10; 
  out[51] += -1.224744871391589*Ghat[23]*dx10; 
  out[52] += -1.224744871391589*Ghat[24]*dx10; 
  out[53] += -0.7071067811865475*Ghat[28]*dx10; 
  out[54] += -1.224744871391589*Ghat[25]*dx10; 
  out[55] += -0.7071067811865475*Ghat[29]*dx10; 
  out[56] += -0.7071067811865475*Ghat[30]*dx10; 
  out[57] += -1.224744871391589*Ghat[26]*dx10; 
  out[58] += -1.224744871391589*Ghat[27]*dx10; 
  out[59] += -1.224744871391589*Ghat[28]*dx10; 
  out[60] += -1.224744871391589*Ghat[29]*dx10; 
  out[61] += -1.224744871391589*Ghat[30]*dx10; 
  out[62] += -0.7071067811865475*Ghat[31]*dx10; 
  out[63] += -1.224744871391589*Ghat[31]*dx10; 

  } else { 

  if (wv>0) { 

  Ghat[0] = (1.224744871391589*fEdge[1]+0.7071067811865475*fEdge[0])*wv+(0.3535533905932737*fEdge[10]+0.2041241452319315*fEdge[4])*dv; 
  Ghat[1] = (1.224744871391589*fEdge[7]+0.7071067811865475*fEdge[2])*wv+(0.3535533905932737*fEdge[23]+0.2041241452319315*fEdge[11])*dv; 
  Ghat[2] = (1.224744871391589*fEdge[8]+0.7071067811865475*fEdge[3])*wv+(0.3535533905932737*fEdge[24]+0.2041241452319315*fEdge[12])*dv; 
  Ghat[3] = (1.224744871391589*fEdge[10]+0.7071067811865475*fEdge[4])*wv+(0.3535533905932737*fEdge[1]+0.2041241452319315*fEdge[0])*dv; 
  Ghat[4] = (1.224744871391589*fEdge[13]+0.7071067811865475*fEdge[5])*wv+(0.3535533905932737*fEdge[29]+0.2041241452319315*fEdge[16])*dv; 
  Ghat[5] = (1.224744871391589*fEdge[17]+0.7071067811865475*fEdge[6])*wv+(0.3535533905932737*fEdge[35]+0.2041241452319315*fEdge[20])*dv; 
  Ghat[6] = (1.224744871391589*fEdge[22]+0.7071067811865475*fEdge[9])*wv+(0.3535533905932737*fEdge[42]+0.2041241452319315*fEdge[25])*dv; 
  Ghat[7] = (1.224744871391589*fEdge[23]+0.7071067811865475*fEdge[11])*wv+(0.3535533905932737*fEdge[7]+0.2041241452319315*fEdge[2])*dv; 
  Ghat[8] = (1.224744871391589*fEdge[24]+0.7071067811865475*fEdge[12])*wv+(0.3535533905932737*fEdge[8]+0.2041241452319315*fEdge[3])*dv; 
  Ghat[9] = (1.224744871391589*fEdge[26]+0.7071067811865475*fEdge[14])*wv+(0.3535533905932737*fEdge[44]+0.2041241452319315*fEdge[30])*dv; 
  Ghat[10] = (1.224744871391589*fEdge[27]+0.7071067811865475*fEdge[15])*wv+(0.3535533905932737*fEdge[45]+0.2041241452319315*fEdge[31])*dv; 
  Ghat[11] = (1.224744871391589*fEdge[29]+0.7071067811865475*fEdge[16])*wv+(0.3535533905932737*fEdge[13]+0.2041241452319315*fEdge[5])*dv; 
  Ghat[12] = (1.224744871391589*fEdge[32]+0.7071067811865475*fEdge[18])*wv+(0.3535533905932737*fEdge[48]+0.2041241452319315*fEdge[36])*dv; 
  Ghat[13] = (1.224744871391589*fEdge[33]+0.7071067811865475*fEdge[19])*wv+(0.3535533905932737*fEdge[49]+0.2041241452319315*fEdge[37])*dv; 
  Ghat[14] = (1.224744871391589*fEdge[35]+0.7071067811865475*fEdge[20])*wv+(0.3535533905932737*fEdge[17]+0.2041241452319315*fEdge[6])*dv; 
  Ghat[15] = (1.224744871391589*fEdge[38]+0.7071067811865475*fEdge[21])*wv+(0.3535533905932737*fEdge[54]+0.2041241452319315*fEdge[41])*dv; 
  Ghat[16] = (1.224744871391589*fEdge[42]+0.7071067811865475*fEdge[25])*wv+(0.3535533905932737*fEdge[22]+0.2041241452319315*fEdge[9])*dv; 
  Ghat[17] = (1.224744871391589*fEdge[43]+0.7071067811865475*fEdge[28])*wv+(0.3535533905932737*fEdge[57]+0.2041241452319315*fEdge[46])*dv; 
  Ghat[18] = (1.224744871391589*fEdge[44]+0.7071067811865475*fEdge[30])*wv+(0.3535533905932737*fEdge[26]+0.2041241452319315*fEdge[14])*dv; 
  Ghat[19] = (1.224744871391589*fEdge[45]+0.7071067811865475*fEdge[31])*wv+(0.3535533905932737*fEdge[27]+0.2041241452319315*fEdge[15])*dv; 
  Ghat[20] = (1.224744871391589*fEdge[47]+0.7071067811865475*fEdge[34])*wv+(0.3535533905932737*fEdge[58]+0.2041241452319315*fEdge[50])*dv; 
  Ghat[21] = (1.224744871391589*fEdge[48]+0.7071067811865475*fEdge[36])*wv+(0.3535533905932737*fEdge[32]+0.2041241452319315*fEdge[18])*dv; 
  Ghat[22] = (1.224744871391589*fEdge[49]+0.7071067811865475*fEdge[37])*wv+(0.3535533905932737*fEdge[33]+0.2041241452319315*fEdge[19])*dv; 
  Ghat[23] = (1.224744871391589*fEdge[51]+0.7071067811865475*fEdge[39])*wv+(0.3535533905932737*fEdge[60]+0.2041241452319315*fEdge[55])*dv; 
  Ghat[24] = (1.224744871391589*fEdge[52]+0.7071067811865475*fEdge[40])*wv+(0.3535533905932737*fEdge[61]+0.2041241452319315*fEdge[56])*dv; 
  Ghat[25] = (1.224744871391589*fEdge[54]+0.7071067811865475*fEdge[41])*wv+(0.3535533905932737*fEdge[38]+0.2041241452319315*fEdge[21])*dv; 
  Ghat[26] = (1.224744871391589*fEdge[57]+0.7071067811865475*fEdge[46])*wv+(0.3535533905932737*fEdge[43]+0.2041241452319315*fEdge[28])*dv; 
  Ghat[27] = (1.224744871391589*fEdge[58]+0.7071067811865475*fEdge[50])*wv+(0.3535533905932737*fEdge[47]+0.2041241452319315*fEdge[34])*dv; 
  Ghat[28] = (1.224744871391589*fEdge[59]+0.7071067811865475*fEdge[53])*wv+(0.3535533905932737*fEdge[63]+0.2041241452319315*fEdge[62])*dv; 
  Ghat[29] = (1.224744871391589*fEdge[60]+0.7071067811865475*fEdge[55])*wv+(0.3535533905932737*fEdge[51]+0.2041241452319315*fEdge[39])*dv; 
  Ghat[30] = (1.224744871391589*fEdge[61]+0.7071067811865475*fEdge[56])*wv+(0.3535533905932737*fEdge[52]+0.2041241452319315*fEdge[40])*dv; 
  Ghat[31] = (1.224744871391589*fEdge[63]+0.7071067811865475*fEdge[62])*wv+(0.3535533905932737*fEdge[59]+0.2041241452319315*fEdge[53])*dv; 
  Ghat[32] = (0.3162277660168379*fEdge[10]+0.1825741858350554*fEdge[4])*dv; 
  Ghat[33] = (0.3162277660168379*fEdge[23]+0.1825741858350554*fEdge[11])*dv; 
  Ghat[34] = (0.3162277660168379*fEdge[24]+0.1825741858350554*fEdge[12])*dv; 
  Ghat[35] = (0.3162277660168379*fEdge[29]+0.1825741858350554*fEdge[16])*dv; 
  Ghat[36] = (0.3162277660168379*fEdge[35]+0.1825741858350554*fEdge[20])*dv; 
  Ghat[37] = (0.3162277660168379*fEdge[42]+0.1825741858350554*fEdge[25])*dv; 
  Ghat[38] = (0.3162277660168379*fEdge[44]+0.1825741858350554*fEdge[30])*dv; 
  Ghat[39] = (0.3162277660168379*fEdge[45]+0.1825741858350554*fEdge[31])*dv; 
  Ghat[40] = (0.3162277660168379*fEdge[48]+0.1825741858350554*fEdge[36])*dv; 
  Ghat[41] = (0.3162277660168379*fEdge[49]+0.1825741858350554*fEdge[37])*dv; 
  Ghat[42] = (0.3162277660168379*fEdge[54]+0.1825741858350554*fEdge[41])*dv; 
  Ghat[43] = (0.3162277660168379*fEdge[57]+0.1825741858350554*fEdge[46])*dv; 
  Ghat[44] = (0.3162277660168379*fEdge[58]+0.1825741858350554*fEdge[50])*dv; 
  Ghat[45] = (0.3162277660168379*fEdge[60]+0.1825741858350554*fEdge[55])*dv; 
  Ghat[46] = (0.3162277660168379*fEdge[61]+0.1825741858350554*fEdge[56])*dv; 
  Ghat[47] = (0.3162277660168379*fEdge[63]+0.1825741858350554*fEdge[62])*dv; 

  } else { 

  Ghat[0] = -0.1178511301977579*((10.39230484541326*fSkin[1]-6.0*fSkin[0])*wv+(3.0*fSkin[10]-1.732050807568877*fSkin[4])*dv); 
  Ghat[1] = -0.1178511301977579*((10.39230484541326*fSkin[7]-6.0*fSkin[2])*wv+(3.0*fSkin[23]-1.732050807568877*fSkin[11])*dv); 
  Ghat[2] = -0.1178511301977579*((10.39230484541326*fSkin[8]-6.0*fSkin[3])*wv+(3.0*fSkin[24]-1.732050807568877*fSkin[12])*dv); 
  Ghat[3] = -0.1178511301977579*((10.39230484541326*fSkin[10]-6.0*fSkin[4])*wv+(3.0*fSkin[1]-1.732050807568877*fSkin[0])*dv); 
  Ghat[4] = -0.1178511301977579*((10.39230484541326*fSkin[13]-6.0*fSkin[5])*wv+(3.0*fSkin[29]-1.732050807568877*fSkin[16])*dv); 
  Ghat[5] = -0.1178511301977579*((10.39230484541326*fSkin[17]-6.0*fSkin[6])*wv+(3.0*fSkin[35]-1.732050807568877*fSkin[20])*dv); 
  Ghat[6] = -0.1178511301977579*((10.39230484541326*fSkin[22]-6.0*fSkin[9])*wv+(3.0*fSkin[42]-1.732050807568877*fSkin[25])*dv); 
  Ghat[7] = -0.1178511301977579*((10.39230484541326*fSkin[23]-6.0*fSkin[11])*wv+(3.0*fSkin[7]-1.732050807568877*fSkin[2])*dv); 
  Ghat[8] = -0.1178511301977579*((10.39230484541326*fSkin[24]-6.0*fSkin[12])*wv+(3.0*fSkin[8]-1.732050807568877*fSkin[3])*dv); 
  Ghat[9] = -0.1178511301977579*((10.39230484541326*fSkin[26]-6.0*fSkin[14])*wv+(3.0*fSkin[44]-1.732050807568877*fSkin[30])*dv); 
  Ghat[10] = -0.1178511301977579*((10.39230484541326*fSkin[27]-6.0*fSkin[15])*wv+(3.0*fSkin[45]-1.732050807568877*fSkin[31])*dv); 
  Ghat[11] = -0.1178511301977579*((10.39230484541326*fSkin[29]-6.0*fSkin[16])*wv+(3.0*fSkin[13]-1.732050807568877*fSkin[5])*dv); 
  Ghat[12] = -0.1178511301977579*((10.39230484541326*fSkin[32]-6.0*fSkin[18])*wv+(3.0*fSkin[48]-1.732050807568877*fSkin[36])*dv); 
  Ghat[13] = -0.1178511301977579*((10.39230484541326*fSkin[33]-6.0*fSkin[19])*wv+(3.0*fSkin[49]-1.732050807568877*fSkin[37])*dv); 
  Ghat[14] = -0.1178511301977579*((10.39230484541326*fSkin[35]-6.0*fSkin[20])*wv+(3.0*fSkin[17]-1.732050807568877*fSkin[6])*dv); 
  Ghat[15] = -0.1178511301977579*((10.39230484541326*fSkin[38]-6.0*fSkin[21])*wv+(3.0*fSkin[54]-1.732050807568877*fSkin[41])*dv); 
  Ghat[16] = -0.1178511301977579*((10.39230484541326*fSkin[42]-6.0*fSkin[25])*wv+(3.0*fSkin[22]-1.732050807568877*fSkin[9])*dv); 
  Ghat[17] = -0.1178511301977579*((10.39230484541326*fSkin[43]-6.0*fSkin[28])*wv+(3.0*fSkin[57]-1.732050807568877*fSkin[46])*dv); 
  Ghat[18] = -0.1178511301977579*((10.39230484541326*fSkin[44]-6.0*fSkin[30])*wv+(3.0*fSkin[26]-1.732050807568877*fSkin[14])*dv); 
  Ghat[19] = -0.1178511301977579*((10.39230484541326*fSkin[45]-6.0*fSkin[31])*wv+(3.0*fSkin[27]-1.732050807568877*fSkin[15])*dv); 
  Ghat[20] = -0.1178511301977579*((10.39230484541326*fSkin[47]-6.0*fSkin[34])*wv+(3.0*fSkin[58]-1.732050807568877*fSkin[50])*dv); 
  Ghat[21] = -0.1178511301977579*((10.39230484541326*fSkin[48]-6.0*fSkin[36])*wv+(3.0*fSkin[32]-1.732050807568877*fSkin[18])*dv); 
  Ghat[22] = -0.1178511301977579*((10.39230484541326*fSkin[49]-6.0*fSkin[37])*wv+(3.0*fSkin[33]-1.732050807568877*fSkin[19])*dv); 
  Ghat[23] = -0.1178511301977579*((10.39230484541326*fSkin[51]-6.0*fSkin[39])*wv+(3.0*fSkin[60]-1.732050807568877*fSkin[55])*dv); 
  Ghat[24] = -0.1178511301977579*((10.39230484541326*fSkin[52]-6.0*fSkin[40])*wv+(3.0*fSkin[61]-1.732050807568877*fSkin[56])*dv); 
  Ghat[25] = -0.1178511301977579*((10.39230484541326*fSkin[54]-6.0*fSkin[41])*wv+(3.0*fSkin[38]-1.732050807568877*fSkin[21])*dv); 
  Ghat[26] = -0.1178511301977579*((10.39230484541326*fSkin[57]-6.0*fSkin[46])*wv+(3.0*fSkin[43]-1.732050807568877*fSkin[28])*dv); 
  Ghat[27] = -0.1178511301977579*((10.39230484541326*fSkin[58]-6.0*fSkin[50])*wv+(3.0*fSkin[47]-1.732050807568877*fSkin[34])*dv); 
  Ghat[28] = -0.1178511301977579*((10.39230484541326*fSkin[59]-6.0*fSkin[53])*wv+(3.0*fSkin[63]-1.732050807568877*fSkin[62])*dv); 
  Ghat[29] = -0.1178511301977579*((10.39230484541326*fSkin[60]-6.0*fSkin[55])*wv+(3.0*fSkin[51]-1.732050807568877*fSkin[39])*dv); 
  Ghat[30] = -0.1178511301977579*((10.39230484541326*fSkin[61]-6.0*fSkin[56])*wv+(3.0*fSkin[52]-1.732050807568877*fSkin[40])*dv); 
  Ghat[31] = -0.1178511301977579*((10.39230484541326*fSkin[63]-6.0*fSkin[62])*wv+(3.0*fSkin[59]-1.732050807568877*fSkin[53])*dv); 
  Ghat[32] = -0.04714045207910316*(6.708203932499369*fSkin[10]-3.872983346207417*fSkin[4])*dv; 
  Ghat[33] = -0.04714045207910316*(6.708203932499369*fSkin[23]-3.872983346207417*fSkin[11])*dv; 
  Ghat[34] = -0.04714045207910316*(6.708203932499369*fSkin[24]-3.872983346207417*fSkin[12])*dv; 
  Ghat[35] = -0.04714045207910316*(6.708203932499369*fSkin[29]-3.872983346207417*fSkin[16])*dv; 
  Ghat[36] = -0.04714045207910316*(6.708203932499369*fSkin[35]-3.872983346207417*fSkin[20])*dv; 
  Ghat[37] = -0.04714045207910316*(6.708203932499369*fSkin[42]-3.872983346207417*fSkin[25])*dv; 
  Ghat[38] = -0.04714045207910316*(6.708203932499369*fSkin[44]-3.872983346207417*fSkin[30])*dv; 
  Ghat[39] = -0.04714045207910316*(6.708203932499369*fSkin[45]-3.872983346207417*fSkin[31])*dv; 
  Ghat[40] = -0.04714045207910316*(6.708203932499369*fSkin[48]-3.872983346207417*fSkin[36])*dv; 
  Ghat[41] = -0.04714045207910316*(6.708203932499369*fSkin[49]-3.872983346207417*fSkin[37])*dv; 
  Ghat[42] = -0.04714045207910316*(6.708203932499369*fSkin[54]-3.872983346207417*fSkin[41])*dv; 
  Ghat[43] = -0.04714045207910316*(6.708203932499369*fSkin[57]-3.872983346207417*fSkin[46])*dv; 
  Ghat[44] = -0.04714045207910316*(6.708203932499369*fSkin[58]-3.872983346207417*fSkin[50])*dv; 
  Ghat[45] = -0.04714045207910316*(6.708203932499369*fSkin[60]-3.872983346207417*fSkin[55])*dv; 
  Ghat[46] = -0.04714045207910316*(6.708203932499369*fSkin[61]-3.872983346207417*fSkin[56])*dv; 
  Ghat[47] = -0.04714045207910316*(6.708203932499369*fSkin[63]-3.872983346207417*fSkin[62])*dv; 

  } 

  out[0] += 0.7071067811865475*Ghat[0]*dx10; 
  out[1] += -1.224744871391589*Ghat[0]*dx10; 
  out[2] += 0.7071067811865475*Ghat[1]*dx10; 
  out[3] += 0.7071067811865475*Ghat[2]*dx10; 
  out[4] += 0.7071067811865475*Ghat[3]*dx10; 
  out[5] += 0.7071067811865475*Ghat[4]*dx10; 
  out[6] += 0.7071067811865475*Ghat[5]*dx10; 
  out[7] += -1.224744871391589*Ghat[1]*dx10; 
  out[8] += -1.224744871391589*Ghat[2]*dx10; 
  out[9] += 0.7071067811865475*Ghat[6]*dx10; 
  out[10] += -1.224744871391589*Ghat[3]*dx10; 
  out[11] += 0.7071067811865475*Ghat[7]*dx10; 
  out[12] += 0.7071067811865475*Ghat[8]*dx10; 
  out[13] += -1.224744871391589*Ghat[4]*dx10; 
  out[14] += 0.7071067811865475*Ghat[9]*dx10; 
  out[15] += 0.7071067811865475*Ghat[10]*dx10; 
  out[16] += 0.7071067811865475*Ghat[11]*dx10; 
  out[17] += -1.224744871391589*Ghat[5]*dx10; 
  out[18] += 0.7071067811865475*Ghat[12]*dx10; 
  out[19] += 0.7071067811865475*Ghat[13]*dx10; 
  out[20] += 0.7071067811865475*Ghat[14]*dx10; 
  out[21] += 0.7071067811865475*Ghat[15]*dx10; 
  out[22] += -1.224744871391589*Ghat[6]*dx10; 
  out[23] += -1.224744871391589*Ghat[7]*dx10; 
  out[24] += -1.224744871391589*Ghat[8]*dx10; 
  out[25] += 0.7071067811865475*Ghat[16]*dx10; 
  out[26] += -1.224744871391589*Ghat[9]*dx10; 
  out[27] += -1.224744871391589*Ghat[10]*dx10; 
  out[28] += 0.7071067811865475*Ghat[17]*dx10; 
  out[29] += -1.224744871391589*Ghat[11]*dx10; 
  out[30] += 0.7071067811865475*Ghat[18]*dx10; 
  out[31] += 0.7071067811865475*Ghat[19]*dx10; 
  out[32] += -1.224744871391589*Ghat[12]*dx10; 
  out[33] += -1.224744871391589*Ghat[13]*dx10; 
  out[34] += 0.7071067811865475*Ghat[20]*dx10; 
  out[35] += -1.224744871391589*Ghat[14]*dx10; 
  out[36] += 0.7071067811865475*Ghat[21]*dx10; 
  out[37] += 0.7071067811865475*Ghat[22]*dx10; 
  out[38] += -1.224744871391589*Ghat[15]*dx10; 
  out[39] += 0.7071067811865475*Ghat[23]*dx10; 
  out[40] += 0.7071067811865475*Ghat[24]*dx10; 
  out[41] += 0.7071067811865475*Ghat[25]*dx10; 
  out[42] += -1.224744871391589*Ghat[16]*dx10; 
  out[43] += -1.224744871391589*Ghat[17]*dx10; 
  out[44] += -1.224744871391589*Ghat[18]*dx10; 
  out[45] += -1.224744871391589*Ghat[19]*dx10; 
  out[46] += 0.7071067811865475*Ghat[26]*dx10; 
  out[47] += -1.224744871391589*Ghat[20]*dx10; 
  out[48] += -1.224744871391589*Ghat[21]*dx10; 
  out[49] += -1.224744871391589*Ghat[22]*dx10; 
  out[50] += 0.7071067811865475*Ghat[27]*dx10; 
  out[51] += -1.224744871391589*Ghat[23]*dx10; 
  out[52] += -1.224744871391589*Ghat[24]*dx10; 
  out[53] += 0.7071067811865475*Ghat[28]*dx10; 
  out[54] += -1.224744871391589*Ghat[25]*dx10; 
  out[55] += 0.7071067811865475*Ghat[29]*dx10; 
  out[56] += 0.7071067811865475*Ghat[30]*dx10; 
  out[57] += -1.224744871391589*Ghat[26]*dx10; 
  out[58] += -1.224744871391589*Ghat[27]*dx10; 
  out[59] += -1.224744871391589*Ghat[28]*dx10; 
  out[60] += -1.224744871391589*Ghat[29]*dx10; 
  out[61] += -1.224744871391589*Ghat[30]*dx10; 
  out[62] += 0.7071067811865475*Ghat[31]*dx10; 
  out[63] += -1.224744871391589*Ghat[31]*dx10; 

  } 
  return 0.;

} 
