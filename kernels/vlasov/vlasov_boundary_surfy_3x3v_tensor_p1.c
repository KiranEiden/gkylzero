#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_boundary_surfy_3x3v_tensor_p1(const double *w, const double *dxv, const double *alpha_geo, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:     Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // alpha_geo:   Fields used only for general geometry.
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Incremented distribution function in center cell.
  const double dx11 = 2/dxv[1]; 
  const double dv = dxv[4], wv = w[4]; 
  double Ghat[80]; 

  if (edge == -1) { 

  if (wv>0) { 

  Ghat[0] = (1.224744871391589*fSkin[2]+0.7071067811865475*fSkin[0])*wv+(0.3535533905932737*fSkin[14]+0.2041241452319315*fSkin[5])*dv; 
  Ghat[1] = (1.224744871391589*fSkin[7]+0.7071067811865475*fSkin[1])*wv+(0.3535533905932737*fSkin[26]+0.2041241452319315*fSkin[13])*dv; 
  Ghat[2] = (1.224744871391589*fSkin[9]+0.7071067811865475*fSkin[3])*wv+(0.3535533905932737*fSkin[28]+0.2041241452319315*fSkin[15])*dv; 
  Ghat[3] = (1.224744871391589*fSkin[11]+0.7071067811865475*fSkin[4])*wv+(0.3535533905932737*fSkin[30]+0.2041241452319315*fSkin[16])*dv; 
  Ghat[4] = (1.224744871391589*fSkin[14]+0.7071067811865475*fSkin[5])*wv+(0.3535533905932737*fSkin[2]+0.2041241452319315*fSkin[0])*dv; 
  Ghat[5] = (1.224744871391589*fSkin[18]+0.7071067811865475*fSkin[6])*wv+(0.3535533905932737*fSkin[39]+0.2041241452319315*fSkin[21])*dv; 
  Ghat[6] = (1.224744871391589*fSkin[22]+0.7071067811865475*fSkin[8])*wv+(0.3535533905932737*fSkin[43]+0.2041241452319315*fSkin[27])*dv; 
  Ghat[7] = (1.224744871391589*fSkin[23]+0.7071067811865475*fSkin[10])*wv+(0.3535533905932737*fSkin[44]+0.2041241452319315*fSkin[29])*dv; 
  Ghat[8] = (1.224744871391589*fSkin[25]+0.7071067811865475*fSkin[12])*wv+(0.3535533905932737*fSkin[46]+0.2041241452319315*fSkin[31])*dv; 
  Ghat[9] = (1.224744871391589*fSkin[26]+0.7071067811865475*fSkin[13])*wv+(0.3535533905932737*fSkin[7]+0.2041241452319315*fSkin[1])*dv; 
  Ghat[10] = (1.224744871391589*fSkin[28]+0.7071067811865475*fSkin[15])*wv+(0.3535533905932737*fSkin[9]+0.2041241452319315*fSkin[3])*dv; 
  Ghat[11] = (1.224744871391589*fSkin[30]+0.7071067811865475*fSkin[16])*wv+(0.3535533905932737*fSkin[11]+0.2041241452319315*fSkin[4])*dv; 
  Ghat[12] = (1.224744871391589*fSkin[32]+0.7071067811865475*fSkin[17])*wv+(0.3535533905932737*fSkin[51]+0.2041241452319315*fSkin[38])*dv; 
  Ghat[13] = (1.224744871391589*fSkin[34]+0.7071067811865475*fSkin[19])*wv+(0.3535533905932737*fSkin[53]+0.2041241452319315*fSkin[40])*dv; 
  Ghat[14] = (1.224744871391589*fSkin[36]+0.7071067811865475*fSkin[20])*wv+(0.3535533905932737*fSkin[55]+0.2041241452319315*fSkin[41])*dv; 
  Ghat[15] = (1.224744871391589*fSkin[39]+0.7071067811865475*fSkin[21])*wv+(0.3535533905932737*fSkin[18]+0.2041241452319315*fSkin[6])*dv; 
  Ghat[16] = (1.224744871391589*fSkin[42]+0.7071067811865475*fSkin[24])*wv+(0.3535533905932737*fSkin[57]+0.2041241452319315*fSkin[45])*dv; 
  Ghat[17] = (1.224744871391589*fSkin[43]+0.7071067811865475*fSkin[27])*wv+(0.3535533905932737*fSkin[22]+0.2041241452319315*fSkin[8])*dv; 
  Ghat[18] = (1.224744871391589*fSkin[44]+0.7071067811865475*fSkin[29])*wv+(0.3535533905932737*fSkin[23]+0.2041241452319315*fSkin[10])*dv; 
  Ghat[19] = (1.224744871391589*fSkin[46]+0.7071067811865475*fSkin[31])*wv+(0.3535533905932737*fSkin[25]+0.2041241452319315*fSkin[12])*dv; 
  Ghat[20] = (1.224744871391589*fSkin[47]+0.7071067811865475*fSkin[33])*wv+(0.3535533905932737*fSkin[59]+0.2041241452319315*fSkin[52])*dv; 
  Ghat[21] = (1.224744871391589*fSkin[48]+0.7071067811865475*fSkin[35])*wv+(0.3535533905932737*fSkin[60]+0.2041241452319315*fSkin[54])*dv; 
  Ghat[22] = (1.224744871391589*fSkin[50]+0.7071067811865475*fSkin[37])*wv+(0.3535533905932737*fSkin[62]+0.2041241452319315*fSkin[56])*dv; 
  Ghat[23] = (1.224744871391589*fSkin[51]+0.7071067811865475*fSkin[38])*wv+(0.3535533905932737*fSkin[32]+0.2041241452319315*fSkin[17])*dv; 
  Ghat[24] = (1.224744871391589*fSkin[53]+0.7071067811865475*fSkin[40])*wv+(0.3535533905932737*fSkin[34]+0.2041241452319315*fSkin[19])*dv; 
  Ghat[25] = (1.224744871391589*fSkin[55]+0.7071067811865475*fSkin[41])*wv+(0.3535533905932737*fSkin[36]+0.2041241452319315*fSkin[20])*dv; 
  Ghat[26] = (1.224744871391589*fSkin[57]+0.7071067811865475*fSkin[45])*wv+(0.3535533905932737*fSkin[42]+0.2041241452319315*fSkin[24])*dv; 
  Ghat[27] = (1.224744871391589*fSkin[58]+0.7071067811865475*fSkin[49])*wv+(0.3535533905932737*fSkin[63]+0.2041241452319315*fSkin[61])*dv; 
  Ghat[28] = (1.224744871391589*fSkin[59]+0.7071067811865475*fSkin[52])*wv+(0.3535533905932737*fSkin[47]+0.2041241452319315*fSkin[33])*dv; 
  Ghat[29] = (1.224744871391589*fSkin[60]+0.7071067811865475*fSkin[54])*wv+(0.3535533905932737*fSkin[48]+0.2041241452319315*fSkin[35])*dv; 
  Ghat[30] = (1.224744871391589*fSkin[62]+0.7071067811865475*fSkin[56])*wv+(0.3535533905932737*fSkin[50]+0.2041241452319315*fSkin[37])*dv; 
  Ghat[31] = (1.224744871391589*fSkin[63]+0.7071067811865475*fSkin[61])*wv+(0.3535533905932737*fSkin[58]+0.2041241452319315*fSkin[49])*dv; 
  Ghat[48] = (0.3162277660168379*fSkin[14]+0.1825741858350554*fSkin[5])*dv; 
  Ghat[49] = (0.3162277660168379*fSkin[26]+0.1825741858350554*fSkin[13])*dv; 
  Ghat[50] = (0.3162277660168379*fSkin[28]+0.1825741858350554*fSkin[15])*dv; 
  Ghat[51] = (0.3162277660168379*fSkin[30]+0.1825741858350554*fSkin[16])*dv; 
  Ghat[52] = (0.3162277660168379*fSkin[39]+0.1825741858350554*fSkin[21])*dv; 
  Ghat[53] = (0.3162277660168379*fSkin[43]+0.1825741858350554*fSkin[27])*dv; 
  Ghat[54] = (0.3162277660168379*fSkin[44]+0.1825741858350554*fSkin[29])*dv; 
  Ghat[55] = (0.3162277660168379*fSkin[46]+0.1825741858350554*fSkin[31])*dv; 
  Ghat[56] = (0.3162277660168379*fSkin[51]+0.1825741858350554*fSkin[38])*dv; 
  Ghat[57] = (0.3162277660168379*fSkin[53]+0.1825741858350554*fSkin[40])*dv; 
  Ghat[58] = (0.3162277660168379*fSkin[55]+0.1825741858350554*fSkin[41])*dv; 
  Ghat[59] = (0.3162277660168379*fSkin[57]+0.1825741858350554*fSkin[45])*dv; 
  Ghat[60] = (0.3162277660168379*fSkin[59]+0.1825741858350554*fSkin[52])*dv; 
  Ghat[61] = (0.3162277660168379*fSkin[60]+0.1825741858350554*fSkin[54])*dv; 
  Ghat[62] = (0.3162277660168379*fSkin[62]+0.1825741858350554*fSkin[56])*dv; 
  Ghat[63] = (0.3162277660168379*fSkin[63]+0.1825741858350554*fSkin[61])*dv; 

  } else { 

  Ghat[0] = -0.1178511301977579*((10.39230484541326*fEdge[2]-6.0*fEdge[0])*wv+(3.0*fEdge[14]-1.732050807568877*fEdge[5])*dv); 
  Ghat[1] = -0.1178511301977579*((10.39230484541326*fEdge[7]-6.0*fEdge[1])*wv+(3.0*fEdge[26]-1.732050807568877*fEdge[13])*dv); 
  Ghat[2] = -0.1178511301977579*((10.39230484541326*fEdge[9]-6.0*fEdge[3])*wv+(3.0*fEdge[28]-1.732050807568877*fEdge[15])*dv); 
  Ghat[3] = -0.1178511301977579*((10.39230484541326*fEdge[11]-6.0*fEdge[4])*wv+(3.0*fEdge[30]-1.732050807568877*fEdge[16])*dv); 
  Ghat[4] = -0.1178511301977579*((10.39230484541326*fEdge[14]-6.0*fEdge[5])*wv+(3.0*fEdge[2]-1.732050807568877*fEdge[0])*dv); 
  Ghat[5] = -0.1178511301977579*((10.39230484541326*fEdge[18]-6.0*fEdge[6])*wv+(3.0*fEdge[39]-1.732050807568877*fEdge[21])*dv); 
  Ghat[6] = -0.1178511301977579*((10.39230484541326*fEdge[22]-6.0*fEdge[8])*wv+(3.0*fEdge[43]-1.732050807568877*fEdge[27])*dv); 
  Ghat[7] = -0.1178511301977579*((10.39230484541326*fEdge[23]-6.0*fEdge[10])*wv+(3.0*fEdge[44]-1.732050807568877*fEdge[29])*dv); 
  Ghat[8] = -0.1178511301977579*((10.39230484541326*fEdge[25]-6.0*fEdge[12])*wv+(3.0*fEdge[46]-1.732050807568877*fEdge[31])*dv); 
  Ghat[9] = -0.1178511301977579*((10.39230484541326*fEdge[26]-6.0*fEdge[13])*wv+(3.0*fEdge[7]-1.732050807568877*fEdge[1])*dv); 
  Ghat[10] = -0.1178511301977579*((10.39230484541326*fEdge[28]-6.0*fEdge[15])*wv+(3.0*fEdge[9]-1.732050807568877*fEdge[3])*dv); 
  Ghat[11] = -0.1178511301977579*((10.39230484541326*fEdge[30]-6.0*fEdge[16])*wv+(3.0*fEdge[11]-1.732050807568877*fEdge[4])*dv); 
  Ghat[12] = -0.1178511301977579*((10.39230484541326*fEdge[32]-6.0*fEdge[17])*wv+(3.0*fEdge[51]-1.732050807568877*fEdge[38])*dv); 
  Ghat[13] = -0.1178511301977579*((10.39230484541326*fEdge[34]-6.0*fEdge[19])*wv+(3.0*fEdge[53]-1.732050807568877*fEdge[40])*dv); 
  Ghat[14] = -0.1178511301977579*((10.39230484541326*fEdge[36]-6.0*fEdge[20])*wv+(3.0*fEdge[55]-1.732050807568877*fEdge[41])*dv); 
  Ghat[15] = -0.1178511301977579*((10.39230484541326*fEdge[39]-6.0*fEdge[21])*wv+(3.0*fEdge[18]-1.732050807568877*fEdge[6])*dv); 
  Ghat[16] = -0.1178511301977579*((10.39230484541326*fEdge[42]-6.0*fEdge[24])*wv+(3.0*fEdge[57]-1.732050807568877*fEdge[45])*dv); 
  Ghat[17] = -0.1178511301977579*((10.39230484541326*fEdge[43]-6.0*fEdge[27])*wv+(3.0*fEdge[22]-1.732050807568877*fEdge[8])*dv); 
  Ghat[18] = -0.1178511301977579*((10.39230484541326*fEdge[44]-6.0*fEdge[29])*wv+(3.0*fEdge[23]-1.732050807568877*fEdge[10])*dv); 
  Ghat[19] = -0.1178511301977579*((10.39230484541326*fEdge[46]-6.0*fEdge[31])*wv+(3.0*fEdge[25]-1.732050807568877*fEdge[12])*dv); 
  Ghat[20] = -0.1178511301977579*((10.39230484541326*fEdge[47]-6.0*fEdge[33])*wv+(3.0*fEdge[59]-1.732050807568877*fEdge[52])*dv); 
  Ghat[21] = -0.1178511301977579*((10.39230484541326*fEdge[48]-6.0*fEdge[35])*wv+(3.0*fEdge[60]-1.732050807568877*fEdge[54])*dv); 
  Ghat[22] = -0.1178511301977579*((10.39230484541326*fEdge[50]-6.0*fEdge[37])*wv+(3.0*fEdge[62]-1.732050807568877*fEdge[56])*dv); 
  Ghat[23] = -0.1178511301977579*((10.39230484541326*fEdge[51]-6.0*fEdge[38])*wv+(3.0*fEdge[32]-1.732050807568877*fEdge[17])*dv); 
  Ghat[24] = -0.1178511301977579*((10.39230484541326*fEdge[53]-6.0*fEdge[40])*wv+(3.0*fEdge[34]-1.732050807568877*fEdge[19])*dv); 
  Ghat[25] = -0.1178511301977579*((10.39230484541326*fEdge[55]-6.0*fEdge[41])*wv+(3.0*fEdge[36]-1.732050807568877*fEdge[20])*dv); 
  Ghat[26] = -0.1178511301977579*((10.39230484541326*fEdge[57]-6.0*fEdge[45])*wv+(3.0*fEdge[42]-1.732050807568877*fEdge[24])*dv); 
  Ghat[27] = -0.1178511301977579*((10.39230484541326*fEdge[58]-6.0*fEdge[49])*wv+(3.0*fEdge[63]-1.732050807568877*fEdge[61])*dv); 
  Ghat[28] = -0.1178511301977579*((10.39230484541326*fEdge[59]-6.0*fEdge[52])*wv+(3.0*fEdge[47]-1.732050807568877*fEdge[33])*dv); 
  Ghat[29] = -0.1178511301977579*((10.39230484541326*fEdge[60]-6.0*fEdge[54])*wv+(3.0*fEdge[48]-1.732050807568877*fEdge[35])*dv); 
  Ghat[30] = -0.1178511301977579*((10.39230484541326*fEdge[62]-6.0*fEdge[56])*wv+(3.0*fEdge[50]-1.732050807568877*fEdge[37])*dv); 
  Ghat[31] = -0.1178511301977579*((10.39230484541326*fEdge[63]-6.0*fEdge[61])*wv+(3.0*fEdge[58]-1.732050807568877*fEdge[49])*dv); 
  Ghat[48] = -0.04714045207910316*(6.708203932499369*fEdge[14]-3.872983346207417*fEdge[5])*dv; 
  Ghat[49] = -0.04714045207910316*(6.708203932499369*fEdge[26]-3.872983346207417*fEdge[13])*dv; 
  Ghat[50] = -0.04714045207910316*(6.708203932499369*fEdge[28]-3.872983346207417*fEdge[15])*dv; 
  Ghat[51] = -0.04714045207910316*(6.708203932499369*fEdge[30]-3.872983346207417*fEdge[16])*dv; 
  Ghat[52] = -0.04714045207910316*(6.708203932499369*fEdge[39]-3.872983346207417*fEdge[21])*dv; 
  Ghat[53] = -0.04714045207910316*(6.708203932499369*fEdge[43]-3.872983346207417*fEdge[27])*dv; 
  Ghat[54] = -0.04714045207910316*(6.708203932499369*fEdge[44]-3.872983346207417*fEdge[29])*dv; 
  Ghat[55] = -0.04714045207910316*(6.708203932499369*fEdge[46]-3.872983346207417*fEdge[31])*dv; 
  Ghat[56] = -0.04714045207910316*(6.708203932499369*fEdge[51]-3.872983346207417*fEdge[38])*dv; 
  Ghat[57] = -0.04714045207910316*(6.708203932499369*fEdge[53]-3.872983346207417*fEdge[40])*dv; 
  Ghat[58] = -0.04714045207910316*(6.708203932499369*fEdge[55]-3.872983346207417*fEdge[41])*dv; 
  Ghat[59] = -0.04714045207910316*(6.708203932499369*fEdge[57]-3.872983346207417*fEdge[45])*dv; 
  Ghat[60] = -0.04714045207910316*(6.708203932499369*fEdge[59]-3.872983346207417*fEdge[52])*dv; 
  Ghat[61] = -0.04714045207910316*(6.708203932499369*fEdge[60]-3.872983346207417*fEdge[54])*dv; 
  Ghat[62] = -0.04714045207910316*(6.708203932499369*fEdge[62]-3.872983346207417*fEdge[56])*dv; 
  Ghat[63] = -0.04714045207910316*(6.708203932499369*fEdge[63]-3.872983346207417*fEdge[61])*dv; 

  } 

  out[0] += -0.7071067811865475*Ghat[0]*dx11; 
  out[1] += -0.7071067811865475*Ghat[1]*dx11; 
  out[2] += -1.224744871391589*Ghat[0]*dx11; 
  out[3] += -0.7071067811865475*Ghat[2]*dx11; 
  out[4] += -0.7071067811865475*Ghat[3]*dx11; 
  out[5] += -0.7071067811865475*Ghat[4]*dx11; 
  out[6] += -0.7071067811865475*Ghat[5]*dx11; 
  out[7] += -1.224744871391589*Ghat[1]*dx11; 
  out[8] += -0.7071067811865475*Ghat[6]*dx11; 
  out[9] += -1.224744871391589*Ghat[2]*dx11; 
  out[10] += -0.7071067811865475*Ghat[7]*dx11; 
  out[11] += -1.224744871391589*Ghat[3]*dx11; 
  out[12] += -0.7071067811865475*Ghat[8]*dx11; 
  out[13] += -0.7071067811865475*Ghat[9]*dx11; 
  out[14] += -1.224744871391589*Ghat[4]*dx11; 
  out[15] += -0.7071067811865475*Ghat[10]*dx11; 
  out[16] += -0.7071067811865475*Ghat[11]*dx11; 
  out[17] += -0.7071067811865475*Ghat[12]*dx11; 
  out[18] += -1.224744871391589*Ghat[5]*dx11; 
  out[19] += -0.7071067811865475*Ghat[13]*dx11; 
  out[20] += -0.7071067811865475*Ghat[14]*dx11; 
  out[21] += -0.7071067811865475*Ghat[15]*dx11; 
  out[22] += -1.224744871391589*Ghat[6]*dx11; 
  out[23] += -1.224744871391589*Ghat[7]*dx11; 
  out[24] += -0.7071067811865475*Ghat[16]*dx11; 
  out[25] += -1.224744871391589*Ghat[8]*dx11; 
  out[26] += -1.224744871391589*Ghat[9]*dx11; 
  out[27] += -0.7071067811865475*Ghat[17]*dx11; 
  out[28] += -1.224744871391589*Ghat[10]*dx11; 
  out[29] += -0.7071067811865475*Ghat[18]*dx11; 
  out[30] += -1.224744871391589*Ghat[11]*dx11; 
  out[31] += -0.7071067811865475*Ghat[19]*dx11; 
  out[32] += -1.224744871391589*Ghat[12]*dx11; 
  out[33] += -0.7071067811865475*Ghat[20]*dx11; 
  out[34] += -1.224744871391589*Ghat[13]*dx11; 
  out[35] += -0.7071067811865475*Ghat[21]*dx11; 
  out[36] += -1.224744871391589*Ghat[14]*dx11; 
  out[37] += -0.7071067811865475*Ghat[22]*dx11; 
  out[38] += -0.7071067811865475*Ghat[23]*dx11; 
  out[39] += -1.224744871391589*Ghat[15]*dx11; 
  out[40] += -0.7071067811865475*Ghat[24]*dx11; 
  out[41] += -0.7071067811865475*Ghat[25]*dx11; 
  out[42] += -1.224744871391589*Ghat[16]*dx11; 
  out[43] += -1.224744871391589*Ghat[17]*dx11; 
  out[44] += -1.224744871391589*Ghat[18]*dx11; 
  out[45] += -0.7071067811865475*Ghat[26]*dx11; 
  out[46] += -1.224744871391589*Ghat[19]*dx11; 
  out[47] += -1.224744871391589*Ghat[20]*dx11; 
  out[48] += -1.224744871391589*Ghat[21]*dx11; 
  out[49] += -0.7071067811865475*Ghat[27]*dx11; 
  out[50] += -1.224744871391589*Ghat[22]*dx11; 
  out[51] += -1.224744871391589*Ghat[23]*dx11; 
  out[52] += -0.7071067811865475*Ghat[28]*dx11; 
  out[53] += -1.224744871391589*Ghat[24]*dx11; 
  out[54] += -0.7071067811865475*Ghat[29]*dx11; 
  out[55] += -1.224744871391589*Ghat[25]*dx11; 
  out[56] += -0.7071067811865475*Ghat[30]*dx11; 
  out[57] += -1.224744871391589*Ghat[26]*dx11; 
  out[58] += -1.224744871391589*Ghat[27]*dx11; 
  out[59] += -1.224744871391589*Ghat[28]*dx11; 
  out[60] += -1.224744871391589*Ghat[29]*dx11; 
  out[61] += -0.7071067811865475*Ghat[31]*dx11; 
  out[62] += -1.224744871391589*Ghat[30]*dx11; 
  out[63] += -1.224744871391589*Ghat[31]*dx11; 

  } else { 

  if (wv>0) { 

  Ghat[0] = (1.224744871391589*fEdge[2]+0.7071067811865475*fEdge[0])*wv+(0.3535533905932737*fEdge[14]+0.2041241452319315*fEdge[5])*dv; 
  Ghat[1] = (1.224744871391589*fEdge[7]+0.7071067811865475*fEdge[1])*wv+(0.3535533905932737*fEdge[26]+0.2041241452319315*fEdge[13])*dv; 
  Ghat[2] = (1.224744871391589*fEdge[9]+0.7071067811865475*fEdge[3])*wv+(0.3535533905932737*fEdge[28]+0.2041241452319315*fEdge[15])*dv; 
  Ghat[3] = (1.224744871391589*fEdge[11]+0.7071067811865475*fEdge[4])*wv+(0.3535533905932737*fEdge[30]+0.2041241452319315*fEdge[16])*dv; 
  Ghat[4] = (1.224744871391589*fEdge[14]+0.7071067811865475*fEdge[5])*wv+(0.3535533905932737*fEdge[2]+0.2041241452319315*fEdge[0])*dv; 
  Ghat[5] = (1.224744871391589*fEdge[18]+0.7071067811865475*fEdge[6])*wv+(0.3535533905932737*fEdge[39]+0.2041241452319315*fEdge[21])*dv; 
  Ghat[6] = (1.224744871391589*fEdge[22]+0.7071067811865475*fEdge[8])*wv+(0.3535533905932737*fEdge[43]+0.2041241452319315*fEdge[27])*dv; 
  Ghat[7] = (1.224744871391589*fEdge[23]+0.7071067811865475*fEdge[10])*wv+(0.3535533905932737*fEdge[44]+0.2041241452319315*fEdge[29])*dv; 
  Ghat[8] = (1.224744871391589*fEdge[25]+0.7071067811865475*fEdge[12])*wv+(0.3535533905932737*fEdge[46]+0.2041241452319315*fEdge[31])*dv; 
  Ghat[9] = (1.224744871391589*fEdge[26]+0.7071067811865475*fEdge[13])*wv+(0.3535533905932737*fEdge[7]+0.2041241452319315*fEdge[1])*dv; 
  Ghat[10] = (1.224744871391589*fEdge[28]+0.7071067811865475*fEdge[15])*wv+(0.3535533905932737*fEdge[9]+0.2041241452319315*fEdge[3])*dv; 
  Ghat[11] = (1.224744871391589*fEdge[30]+0.7071067811865475*fEdge[16])*wv+(0.3535533905932737*fEdge[11]+0.2041241452319315*fEdge[4])*dv; 
  Ghat[12] = (1.224744871391589*fEdge[32]+0.7071067811865475*fEdge[17])*wv+(0.3535533905932737*fEdge[51]+0.2041241452319315*fEdge[38])*dv; 
  Ghat[13] = (1.224744871391589*fEdge[34]+0.7071067811865475*fEdge[19])*wv+(0.3535533905932737*fEdge[53]+0.2041241452319315*fEdge[40])*dv; 
  Ghat[14] = (1.224744871391589*fEdge[36]+0.7071067811865475*fEdge[20])*wv+(0.3535533905932737*fEdge[55]+0.2041241452319315*fEdge[41])*dv; 
  Ghat[15] = (1.224744871391589*fEdge[39]+0.7071067811865475*fEdge[21])*wv+(0.3535533905932737*fEdge[18]+0.2041241452319315*fEdge[6])*dv; 
  Ghat[16] = (1.224744871391589*fEdge[42]+0.7071067811865475*fEdge[24])*wv+(0.3535533905932737*fEdge[57]+0.2041241452319315*fEdge[45])*dv; 
  Ghat[17] = (1.224744871391589*fEdge[43]+0.7071067811865475*fEdge[27])*wv+(0.3535533905932737*fEdge[22]+0.2041241452319315*fEdge[8])*dv; 
  Ghat[18] = (1.224744871391589*fEdge[44]+0.7071067811865475*fEdge[29])*wv+(0.3535533905932737*fEdge[23]+0.2041241452319315*fEdge[10])*dv; 
  Ghat[19] = (1.224744871391589*fEdge[46]+0.7071067811865475*fEdge[31])*wv+(0.3535533905932737*fEdge[25]+0.2041241452319315*fEdge[12])*dv; 
  Ghat[20] = (1.224744871391589*fEdge[47]+0.7071067811865475*fEdge[33])*wv+(0.3535533905932737*fEdge[59]+0.2041241452319315*fEdge[52])*dv; 
  Ghat[21] = (1.224744871391589*fEdge[48]+0.7071067811865475*fEdge[35])*wv+(0.3535533905932737*fEdge[60]+0.2041241452319315*fEdge[54])*dv; 
  Ghat[22] = (1.224744871391589*fEdge[50]+0.7071067811865475*fEdge[37])*wv+(0.3535533905932737*fEdge[62]+0.2041241452319315*fEdge[56])*dv; 
  Ghat[23] = (1.224744871391589*fEdge[51]+0.7071067811865475*fEdge[38])*wv+(0.3535533905932737*fEdge[32]+0.2041241452319315*fEdge[17])*dv; 
  Ghat[24] = (1.224744871391589*fEdge[53]+0.7071067811865475*fEdge[40])*wv+(0.3535533905932737*fEdge[34]+0.2041241452319315*fEdge[19])*dv; 
  Ghat[25] = (1.224744871391589*fEdge[55]+0.7071067811865475*fEdge[41])*wv+(0.3535533905932737*fEdge[36]+0.2041241452319315*fEdge[20])*dv; 
  Ghat[26] = (1.224744871391589*fEdge[57]+0.7071067811865475*fEdge[45])*wv+(0.3535533905932737*fEdge[42]+0.2041241452319315*fEdge[24])*dv; 
  Ghat[27] = (1.224744871391589*fEdge[58]+0.7071067811865475*fEdge[49])*wv+(0.3535533905932737*fEdge[63]+0.2041241452319315*fEdge[61])*dv; 
  Ghat[28] = (1.224744871391589*fEdge[59]+0.7071067811865475*fEdge[52])*wv+(0.3535533905932737*fEdge[47]+0.2041241452319315*fEdge[33])*dv; 
  Ghat[29] = (1.224744871391589*fEdge[60]+0.7071067811865475*fEdge[54])*wv+(0.3535533905932737*fEdge[48]+0.2041241452319315*fEdge[35])*dv; 
  Ghat[30] = (1.224744871391589*fEdge[62]+0.7071067811865475*fEdge[56])*wv+(0.3535533905932737*fEdge[50]+0.2041241452319315*fEdge[37])*dv; 
  Ghat[31] = (1.224744871391589*fEdge[63]+0.7071067811865475*fEdge[61])*wv+(0.3535533905932737*fEdge[58]+0.2041241452319315*fEdge[49])*dv; 
  Ghat[48] = (0.3162277660168379*fEdge[14]+0.1825741858350554*fEdge[5])*dv; 
  Ghat[49] = (0.3162277660168379*fEdge[26]+0.1825741858350554*fEdge[13])*dv; 
  Ghat[50] = (0.3162277660168379*fEdge[28]+0.1825741858350554*fEdge[15])*dv; 
  Ghat[51] = (0.3162277660168379*fEdge[30]+0.1825741858350554*fEdge[16])*dv; 
  Ghat[52] = (0.3162277660168379*fEdge[39]+0.1825741858350554*fEdge[21])*dv; 
  Ghat[53] = (0.3162277660168379*fEdge[43]+0.1825741858350554*fEdge[27])*dv; 
  Ghat[54] = (0.3162277660168379*fEdge[44]+0.1825741858350554*fEdge[29])*dv; 
  Ghat[55] = (0.3162277660168379*fEdge[46]+0.1825741858350554*fEdge[31])*dv; 
  Ghat[56] = (0.3162277660168379*fEdge[51]+0.1825741858350554*fEdge[38])*dv; 
  Ghat[57] = (0.3162277660168379*fEdge[53]+0.1825741858350554*fEdge[40])*dv; 
  Ghat[58] = (0.3162277660168379*fEdge[55]+0.1825741858350554*fEdge[41])*dv; 
  Ghat[59] = (0.3162277660168379*fEdge[57]+0.1825741858350554*fEdge[45])*dv; 
  Ghat[60] = (0.3162277660168379*fEdge[59]+0.1825741858350554*fEdge[52])*dv; 
  Ghat[61] = (0.3162277660168379*fEdge[60]+0.1825741858350554*fEdge[54])*dv; 
  Ghat[62] = (0.3162277660168379*fEdge[62]+0.1825741858350554*fEdge[56])*dv; 
  Ghat[63] = (0.3162277660168379*fEdge[63]+0.1825741858350554*fEdge[61])*dv; 

  } else { 

  Ghat[0] = -0.1178511301977579*((10.39230484541326*fSkin[2]-6.0*fSkin[0])*wv+(3.0*fSkin[14]-1.732050807568877*fSkin[5])*dv); 
  Ghat[1] = -0.1178511301977579*((10.39230484541326*fSkin[7]-6.0*fSkin[1])*wv+(3.0*fSkin[26]-1.732050807568877*fSkin[13])*dv); 
  Ghat[2] = -0.1178511301977579*((10.39230484541326*fSkin[9]-6.0*fSkin[3])*wv+(3.0*fSkin[28]-1.732050807568877*fSkin[15])*dv); 
  Ghat[3] = -0.1178511301977579*((10.39230484541326*fSkin[11]-6.0*fSkin[4])*wv+(3.0*fSkin[30]-1.732050807568877*fSkin[16])*dv); 
  Ghat[4] = -0.1178511301977579*((10.39230484541326*fSkin[14]-6.0*fSkin[5])*wv+(3.0*fSkin[2]-1.732050807568877*fSkin[0])*dv); 
  Ghat[5] = -0.1178511301977579*((10.39230484541326*fSkin[18]-6.0*fSkin[6])*wv+(3.0*fSkin[39]-1.732050807568877*fSkin[21])*dv); 
  Ghat[6] = -0.1178511301977579*((10.39230484541326*fSkin[22]-6.0*fSkin[8])*wv+(3.0*fSkin[43]-1.732050807568877*fSkin[27])*dv); 
  Ghat[7] = -0.1178511301977579*((10.39230484541326*fSkin[23]-6.0*fSkin[10])*wv+(3.0*fSkin[44]-1.732050807568877*fSkin[29])*dv); 
  Ghat[8] = -0.1178511301977579*((10.39230484541326*fSkin[25]-6.0*fSkin[12])*wv+(3.0*fSkin[46]-1.732050807568877*fSkin[31])*dv); 
  Ghat[9] = -0.1178511301977579*((10.39230484541326*fSkin[26]-6.0*fSkin[13])*wv+(3.0*fSkin[7]-1.732050807568877*fSkin[1])*dv); 
  Ghat[10] = -0.1178511301977579*((10.39230484541326*fSkin[28]-6.0*fSkin[15])*wv+(3.0*fSkin[9]-1.732050807568877*fSkin[3])*dv); 
  Ghat[11] = -0.1178511301977579*((10.39230484541326*fSkin[30]-6.0*fSkin[16])*wv+(3.0*fSkin[11]-1.732050807568877*fSkin[4])*dv); 
  Ghat[12] = -0.1178511301977579*((10.39230484541326*fSkin[32]-6.0*fSkin[17])*wv+(3.0*fSkin[51]-1.732050807568877*fSkin[38])*dv); 
  Ghat[13] = -0.1178511301977579*((10.39230484541326*fSkin[34]-6.0*fSkin[19])*wv+(3.0*fSkin[53]-1.732050807568877*fSkin[40])*dv); 
  Ghat[14] = -0.1178511301977579*((10.39230484541326*fSkin[36]-6.0*fSkin[20])*wv+(3.0*fSkin[55]-1.732050807568877*fSkin[41])*dv); 
  Ghat[15] = -0.1178511301977579*((10.39230484541326*fSkin[39]-6.0*fSkin[21])*wv+(3.0*fSkin[18]-1.732050807568877*fSkin[6])*dv); 
  Ghat[16] = -0.1178511301977579*((10.39230484541326*fSkin[42]-6.0*fSkin[24])*wv+(3.0*fSkin[57]-1.732050807568877*fSkin[45])*dv); 
  Ghat[17] = -0.1178511301977579*((10.39230484541326*fSkin[43]-6.0*fSkin[27])*wv+(3.0*fSkin[22]-1.732050807568877*fSkin[8])*dv); 
  Ghat[18] = -0.1178511301977579*((10.39230484541326*fSkin[44]-6.0*fSkin[29])*wv+(3.0*fSkin[23]-1.732050807568877*fSkin[10])*dv); 
  Ghat[19] = -0.1178511301977579*((10.39230484541326*fSkin[46]-6.0*fSkin[31])*wv+(3.0*fSkin[25]-1.732050807568877*fSkin[12])*dv); 
  Ghat[20] = -0.1178511301977579*((10.39230484541326*fSkin[47]-6.0*fSkin[33])*wv+(3.0*fSkin[59]-1.732050807568877*fSkin[52])*dv); 
  Ghat[21] = -0.1178511301977579*((10.39230484541326*fSkin[48]-6.0*fSkin[35])*wv+(3.0*fSkin[60]-1.732050807568877*fSkin[54])*dv); 
  Ghat[22] = -0.1178511301977579*((10.39230484541326*fSkin[50]-6.0*fSkin[37])*wv+(3.0*fSkin[62]-1.732050807568877*fSkin[56])*dv); 
  Ghat[23] = -0.1178511301977579*((10.39230484541326*fSkin[51]-6.0*fSkin[38])*wv+(3.0*fSkin[32]-1.732050807568877*fSkin[17])*dv); 
  Ghat[24] = -0.1178511301977579*((10.39230484541326*fSkin[53]-6.0*fSkin[40])*wv+(3.0*fSkin[34]-1.732050807568877*fSkin[19])*dv); 
  Ghat[25] = -0.1178511301977579*((10.39230484541326*fSkin[55]-6.0*fSkin[41])*wv+(3.0*fSkin[36]-1.732050807568877*fSkin[20])*dv); 
  Ghat[26] = -0.1178511301977579*((10.39230484541326*fSkin[57]-6.0*fSkin[45])*wv+(3.0*fSkin[42]-1.732050807568877*fSkin[24])*dv); 
  Ghat[27] = -0.1178511301977579*((10.39230484541326*fSkin[58]-6.0*fSkin[49])*wv+(3.0*fSkin[63]-1.732050807568877*fSkin[61])*dv); 
  Ghat[28] = -0.1178511301977579*((10.39230484541326*fSkin[59]-6.0*fSkin[52])*wv+(3.0*fSkin[47]-1.732050807568877*fSkin[33])*dv); 
  Ghat[29] = -0.1178511301977579*((10.39230484541326*fSkin[60]-6.0*fSkin[54])*wv+(3.0*fSkin[48]-1.732050807568877*fSkin[35])*dv); 
  Ghat[30] = -0.1178511301977579*((10.39230484541326*fSkin[62]-6.0*fSkin[56])*wv+(3.0*fSkin[50]-1.732050807568877*fSkin[37])*dv); 
  Ghat[31] = -0.1178511301977579*((10.39230484541326*fSkin[63]-6.0*fSkin[61])*wv+(3.0*fSkin[58]-1.732050807568877*fSkin[49])*dv); 
  Ghat[48] = -0.04714045207910316*(6.708203932499369*fSkin[14]-3.872983346207417*fSkin[5])*dv; 
  Ghat[49] = -0.04714045207910316*(6.708203932499369*fSkin[26]-3.872983346207417*fSkin[13])*dv; 
  Ghat[50] = -0.04714045207910316*(6.708203932499369*fSkin[28]-3.872983346207417*fSkin[15])*dv; 
  Ghat[51] = -0.04714045207910316*(6.708203932499369*fSkin[30]-3.872983346207417*fSkin[16])*dv; 
  Ghat[52] = -0.04714045207910316*(6.708203932499369*fSkin[39]-3.872983346207417*fSkin[21])*dv; 
  Ghat[53] = -0.04714045207910316*(6.708203932499369*fSkin[43]-3.872983346207417*fSkin[27])*dv; 
  Ghat[54] = -0.04714045207910316*(6.708203932499369*fSkin[44]-3.872983346207417*fSkin[29])*dv; 
  Ghat[55] = -0.04714045207910316*(6.708203932499369*fSkin[46]-3.872983346207417*fSkin[31])*dv; 
  Ghat[56] = -0.04714045207910316*(6.708203932499369*fSkin[51]-3.872983346207417*fSkin[38])*dv; 
  Ghat[57] = -0.04714045207910316*(6.708203932499369*fSkin[53]-3.872983346207417*fSkin[40])*dv; 
  Ghat[58] = -0.04714045207910316*(6.708203932499369*fSkin[55]-3.872983346207417*fSkin[41])*dv; 
  Ghat[59] = -0.04714045207910316*(6.708203932499369*fSkin[57]-3.872983346207417*fSkin[45])*dv; 
  Ghat[60] = -0.04714045207910316*(6.708203932499369*fSkin[59]-3.872983346207417*fSkin[52])*dv; 
  Ghat[61] = -0.04714045207910316*(6.708203932499369*fSkin[60]-3.872983346207417*fSkin[54])*dv; 
  Ghat[62] = -0.04714045207910316*(6.708203932499369*fSkin[62]-3.872983346207417*fSkin[56])*dv; 
  Ghat[63] = -0.04714045207910316*(6.708203932499369*fSkin[63]-3.872983346207417*fSkin[61])*dv; 

  } 

  out[0] += 0.7071067811865475*Ghat[0]*dx11; 
  out[1] += 0.7071067811865475*Ghat[1]*dx11; 
  out[2] += -1.224744871391589*Ghat[0]*dx11; 
  out[3] += 0.7071067811865475*Ghat[2]*dx11; 
  out[4] += 0.7071067811865475*Ghat[3]*dx11; 
  out[5] += 0.7071067811865475*Ghat[4]*dx11; 
  out[6] += 0.7071067811865475*Ghat[5]*dx11; 
  out[7] += -1.224744871391589*Ghat[1]*dx11; 
  out[8] += 0.7071067811865475*Ghat[6]*dx11; 
  out[9] += -1.224744871391589*Ghat[2]*dx11; 
  out[10] += 0.7071067811865475*Ghat[7]*dx11; 
  out[11] += -1.224744871391589*Ghat[3]*dx11; 
  out[12] += 0.7071067811865475*Ghat[8]*dx11; 
  out[13] += 0.7071067811865475*Ghat[9]*dx11; 
  out[14] += -1.224744871391589*Ghat[4]*dx11; 
  out[15] += 0.7071067811865475*Ghat[10]*dx11; 
  out[16] += 0.7071067811865475*Ghat[11]*dx11; 
  out[17] += 0.7071067811865475*Ghat[12]*dx11; 
  out[18] += -1.224744871391589*Ghat[5]*dx11; 
  out[19] += 0.7071067811865475*Ghat[13]*dx11; 
  out[20] += 0.7071067811865475*Ghat[14]*dx11; 
  out[21] += 0.7071067811865475*Ghat[15]*dx11; 
  out[22] += -1.224744871391589*Ghat[6]*dx11; 
  out[23] += -1.224744871391589*Ghat[7]*dx11; 
  out[24] += 0.7071067811865475*Ghat[16]*dx11; 
  out[25] += -1.224744871391589*Ghat[8]*dx11; 
  out[26] += -1.224744871391589*Ghat[9]*dx11; 
  out[27] += 0.7071067811865475*Ghat[17]*dx11; 
  out[28] += -1.224744871391589*Ghat[10]*dx11; 
  out[29] += 0.7071067811865475*Ghat[18]*dx11; 
  out[30] += -1.224744871391589*Ghat[11]*dx11; 
  out[31] += 0.7071067811865475*Ghat[19]*dx11; 
  out[32] += -1.224744871391589*Ghat[12]*dx11; 
  out[33] += 0.7071067811865475*Ghat[20]*dx11; 
  out[34] += -1.224744871391589*Ghat[13]*dx11; 
  out[35] += 0.7071067811865475*Ghat[21]*dx11; 
  out[36] += -1.224744871391589*Ghat[14]*dx11; 
  out[37] += 0.7071067811865475*Ghat[22]*dx11; 
  out[38] += 0.7071067811865475*Ghat[23]*dx11; 
  out[39] += -1.224744871391589*Ghat[15]*dx11; 
  out[40] += 0.7071067811865475*Ghat[24]*dx11; 
  out[41] += 0.7071067811865475*Ghat[25]*dx11; 
  out[42] += -1.224744871391589*Ghat[16]*dx11; 
  out[43] += -1.224744871391589*Ghat[17]*dx11; 
  out[44] += -1.224744871391589*Ghat[18]*dx11; 
  out[45] += 0.7071067811865475*Ghat[26]*dx11; 
  out[46] += -1.224744871391589*Ghat[19]*dx11; 
  out[47] += -1.224744871391589*Ghat[20]*dx11; 
  out[48] += -1.224744871391589*Ghat[21]*dx11; 
  out[49] += 0.7071067811865475*Ghat[27]*dx11; 
  out[50] += -1.224744871391589*Ghat[22]*dx11; 
  out[51] += -1.224744871391589*Ghat[23]*dx11; 
  out[52] += 0.7071067811865475*Ghat[28]*dx11; 
  out[53] += -1.224744871391589*Ghat[24]*dx11; 
  out[54] += 0.7071067811865475*Ghat[29]*dx11; 
  out[55] += -1.224744871391589*Ghat[25]*dx11; 
  out[56] += 0.7071067811865475*Ghat[30]*dx11; 
  out[57] += -1.224744871391589*Ghat[26]*dx11; 
  out[58] += -1.224744871391589*Ghat[27]*dx11; 
  out[59] += -1.224744871391589*Ghat[28]*dx11; 
  out[60] += -1.224744871391589*Ghat[29]*dx11; 
  out[61] += 0.7071067811865475*Ghat[31]*dx11; 
  out[62] += -1.224744871391589*Ghat[30]*dx11; 
  out[63] += -1.224744871391589*Ghat[31]*dx11; 

  } 
  return 0.;

} 