#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_surfz_3x3v_tensor_p1(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // alpha_geo: Fields used only for general geometry.
  // fl/fc/fr:  Input Distribution function in left/center/right cells.
  // out:       Incremented distribution function in center cell.
  const double dx12 = 2/dxv[2]; 
  const double dv = dxv[5], wv = w[5]; 
  double Ghat_r[80]; 
  double Ghat_l[80]; 
  if (wv>0) { 

  Ghat_r[0] = (1.224744871391589*fc[3]+0.7071067811865475*fc[0])*wv+(0.3535533905932737*fc[19]+0.2041241452319315*fc[6])*dv; 
  Ghat_r[1] = (1.224744871391589*fc[8]+0.7071067811865475*fc[1])*wv+(0.3535533905932737*fc[33]+0.2041241452319315*fc[17])*dv; 
  Ghat_r[2] = (1.224744871391589*fc[9]+0.7071067811865475*fc[2])*wv+(0.3535533905932737*fc[34]+0.2041241452319315*fc[18])*dv; 
  Ghat_r[3] = (1.224744871391589*fc[12]+0.7071067811865475*fc[4])*wv+(0.3535533905932737*fc[37]+0.2041241452319315*fc[20])*dv; 
  Ghat_r[4] = (1.224744871391589*fc[15]+0.7071067811865475*fc[5])*wv+(0.3535533905932737*fc[40]+0.2041241452319315*fc[21])*dv; 
  Ghat_r[5] = (1.224744871391589*fc[19]+0.7071067811865475*fc[6])*wv+(0.3535533905932737*fc[3]+0.2041241452319315*fc[0])*dv; 
  Ghat_r[6] = (1.224744871391589*fc[22]+0.7071067811865475*fc[7])*wv+(0.3535533905932737*fc[47]+0.2041241452319315*fc[32])*dv; 
  Ghat_r[7] = (1.224744871391589*fc[24]+0.7071067811865475*fc[10])*wv+(0.3535533905932737*fc[49]+0.2041241452319315*fc[35])*dv; 
  Ghat_r[8] = (1.224744871391589*fc[25]+0.7071067811865475*fc[11])*wv+(0.3535533905932737*fc[50]+0.2041241452319315*fc[36])*dv; 
  Ghat_r[9] = (1.224744871391589*fc[27]+0.7071067811865475*fc[13])*wv+(0.3535533905932737*fc[52]+0.2041241452319315*fc[38])*dv; 
  Ghat_r[10] = (1.224744871391589*fc[28]+0.7071067811865475*fc[14])*wv+(0.3535533905932737*fc[53]+0.2041241452319315*fc[39])*dv; 
  Ghat_r[11] = (1.224744871391589*fc[31]+0.7071067811865475*fc[16])*wv+(0.3535533905932737*fc[56]+0.2041241452319315*fc[41])*dv; 
  Ghat_r[12] = (1.224744871391589*fc[33]+0.7071067811865475*fc[17])*wv+(0.3535533905932737*fc[8]+0.2041241452319315*fc[1])*dv; 
  Ghat_r[13] = (1.224744871391589*fc[34]+0.7071067811865475*fc[18])*wv+(0.3535533905932737*fc[9]+0.2041241452319315*fc[2])*dv; 
  Ghat_r[14] = (1.224744871391589*fc[37]+0.7071067811865475*fc[20])*wv+(0.3535533905932737*fc[12]+0.2041241452319315*fc[4])*dv; 
  Ghat_r[15] = (1.224744871391589*fc[40]+0.7071067811865475*fc[21])*wv+(0.3535533905932737*fc[15]+0.2041241452319315*fc[5])*dv; 
  Ghat_r[16] = (1.224744871391589*fc[42]+0.7071067811865475*fc[23])*wv+(0.3535533905932737*fc[58]+0.2041241452319315*fc[48])*dv; 
  Ghat_r[17] = (1.224744871391589*fc[43]+0.7071067811865475*fc[26])*wv+(0.3535533905932737*fc[59]+0.2041241452319315*fc[51])*dv; 
  Ghat_r[18] = (1.224744871391589*fc[45]+0.7071067811865475*fc[29])*wv+(0.3535533905932737*fc[61]+0.2041241452319315*fc[54])*dv; 
  Ghat_r[19] = (1.224744871391589*fc[46]+0.7071067811865475*fc[30])*wv+(0.3535533905932737*fc[62]+0.2041241452319315*fc[55])*dv; 
  Ghat_r[20] = (1.224744871391589*fc[47]+0.7071067811865475*fc[32])*wv+(0.3535533905932737*fc[22]+0.2041241452319315*fc[7])*dv; 
  Ghat_r[21] = (1.224744871391589*fc[49]+0.7071067811865475*fc[35])*wv+(0.3535533905932737*fc[24]+0.2041241452319315*fc[10])*dv; 
  Ghat_r[22] = (1.224744871391589*fc[50]+0.7071067811865475*fc[36])*wv+(0.3535533905932737*fc[25]+0.2041241452319315*fc[11])*dv; 
  Ghat_r[23] = (1.224744871391589*fc[52]+0.7071067811865475*fc[38])*wv+(0.3535533905932737*fc[27]+0.2041241452319315*fc[13])*dv; 
  Ghat_r[24] = (1.224744871391589*fc[53]+0.7071067811865475*fc[39])*wv+(0.3535533905932737*fc[28]+0.2041241452319315*fc[14])*dv; 
  Ghat_r[25] = (1.224744871391589*fc[56]+0.7071067811865475*fc[41])*wv+(0.3535533905932737*fc[31]+0.2041241452319315*fc[16])*dv; 
  Ghat_r[26] = (1.224744871391589*fc[57]+0.7071067811865475*fc[44])*wv+(0.3535533905932737*fc[63]+0.2041241452319315*fc[60])*dv; 
  Ghat_r[27] = (1.224744871391589*fc[58]+0.7071067811865475*fc[48])*wv+(0.3535533905932737*fc[42]+0.2041241452319315*fc[23])*dv; 
  Ghat_r[28] = (1.224744871391589*fc[59]+0.7071067811865475*fc[51])*wv+(0.3535533905932737*fc[43]+0.2041241452319315*fc[26])*dv; 
  Ghat_r[29] = (1.224744871391589*fc[61]+0.7071067811865475*fc[54])*wv+(0.3535533905932737*fc[45]+0.2041241452319315*fc[29])*dv; 
  Ghat_r[30] = (1.224744871391589*fc[62]+0.7071067811865475*fc[55])*wv+(0.3535533905932737*fc[46]+0.2041241452319315*fc[30])*dv; 
  Ghat_r[31] = (1.224744871391589*fc[63]+0.7071067811865475*fc[60])*wv+(0.3535533905932737*fc[57]+0.2041241452319315*fc[44])*dv; 
  Ghat_r[64] = (0.3162277660168379*fc[19]+0.1825741858350554*fc[6])*dv; 
  Ghat_r[65] = (0.3162277660168379*fc[33]+0.1825741858350554*fc[17])*dv; 
  Ghat_r[66] = (0.3162277660168379*fc[34]+0.1825741858350554*fc[18])*dv; 
  Ghat_r[67] = (0.3162277660168379*fc[37]+0.1825741858350554*fc[20])*dv; 
  Ghat_r[68] = (0.3162277660168379*fc[40]+0.1825741858350554*fc[21])*dv; 
  Ghat_r[69] = (0.3162277660168379*fc[47]+0.1825741858350554*fc[32])*dv; 
  Ghat_r[70] = (0.3162277660168379*fc[49]+0.1825741858350554*fc[35])*dv; 
  Ghat_r[71] = (0.3162277660168379*fc[50]+0.1825741858350554*fc[36])*dv; 
  Ghat_r[72] = (0.3162277660168379*fc[52]+0.1825741858350554*fc[38])*dv; 
  Ghat_r[73] = (0.3162277660168379*fc[53]+0.1825741858350554*fc[39])*dv; 
  Ghat_r[74] = (0.3162277660168379*fc[56]+0.1825741858350554*fc[41])*dv; 
  Ghat_r[75] = (0.3162277660168379*fc[58]+0.1825741858350554*fc[48])*dv; 
  Ghat_r[76] = (0.3162277660168379*fc[59]+0.1825741858350554*fc[51])*dv; 
  Ghat_r[77] = (0.3162277660168379*fc[61]+0.1825741858350554*fc[54])*dv; 
  Ghat_r[78] = (0.3162277660168379*fc[62]+0.1825741858350554*fc[55])*dv; 
  Ghat_r[79] = (0.3162277660168379*fc[63]+0.1825741858350554*fc[60])*dv; 

  Ghat_l[0] = (1.224744871391589*fl[3]+0.7071067811865475*fl[0])*wv+(0.3535533905932737*fl[19]+0.2041241452319315*fl[6])*dv; 
  Ghat_l[1] = (1.224744871391589*fl[8]+0.7071067811865475*fl[1])*wv+(0.3535533905932737*fl[33]+0.2041241452319315*fl[17])*dv; 
  Ghat_l[2] = (1.224744871391589*fl[9]+0.7071067811865475*fl[2])*wv+(0.3535533905932737*fl[34]+0.2041241452319315*fl[18])*dv; 
  Ghat_l[3] = (1.224744871391589*fl[12]+0.7071067811865475*fl[4])*wv+(0.3535533905932737*fl[37]+0.2041241452319315*fl[20])*dv; 
  Ghat_l[4] = (1.224744871391589*fl[15]+0.7071067811865475*fl[5])*wv+(0.3535533905932737*fl[40]+0.2041241452319315*fl[21])*dv; 
  Ghat_l[5] = (1.224744871391589*fl[19]+0.7071067811865475*fl[6])*wv+(0.3535533905932737*fl[3]+0.2041241452319315*fl[0])*dv; 
  Ghat_l[6] = (1.224744871391589*fl[22]+0.7071067811865475*fl[7])*wv+(0.3535533905932737*fl[47]+0.2041241452319315*fl[32])*dv; 
  Ghat_l[7] = (1.224744871391589*fl[24]+0.7071067811865475*fl[10])*wv+(0.3535533905932737*fl[49]+0.2041241452319315*fl[35])*dv; 
  Ghat_l[8] = (1.224744871391589*fl[25]+0.7071067811865475*fl[11])*wv+(0.3535533905932737*fl[50]+0.2041241452319315*fl[36])*dv; 
  Ghat_l[9] = (1.224744871391589*fl[27]+0.7071067811865475*fl[13])*wv+(0.3535533905932737*fl[52]+0.2041241452319315*fl[38])*dv; 
  Ghat_l[10] = (1.224744871391589*fl[28]+0.7071067811865475*fl[14])*wv+(0.3535533905932737*fl[53]+0.2041241452319315*fl[39])*dv; 
  Ghat_l[11] = (1.224744871391589*fl[31]+0.7071067811865475*fl[16])*wv+(0.3535533905932737*fl[56]+0.2041241452319315*fl[41])*dv; 
  Ghat_l[12] = (1.224744871391589*fl[33]+0.7071067811865475*fl[17])*wv+(0.3535533905932737*fl[8]+0.2041241452319315*fl[1])*dv; 
  Ghat_l[13] = (1.224744871391589*fl[34]+0.7071067811865475*fl[18])*wv+(0.3535533905932737*fl[9]+0.2041241452319315*fl[2])*dv; 
  Ghat_l[14] = (1.224744871391589*fl[37]+0.7071067811865475*fl[20])*wv+(0.3535533905932737*fl[12]+0.2041241452319315*fl[4])*dv; 
  Ghat_l[15] = (1.224744871391589*fl[40]+0.7071067811865475*fl[21])*wv+(0.3535533905932737*fl[15]+0.2041241452319315*fl[5])*dv; 
  Ghat_l[16] = (1.224744871391589*fl[42]+0.7071067811865475*fl[23])*wv+(0.3535533905932737*fl[58]+0.2041241452319315*fl[48])*dv; 
  Ghat_l[17] = (1.224744871391589*fl[43]+0.7071067811865475*fl[26])*wv+(0.3535533905932737*fl[59]+0.2041241452319315*fl[51])*dv; 
  Ghat_l[18] = (1.224744871391589*fl[45]+0.7071067811865475*fl[29])*wv+(0.3535533905932737*fl[61]+0.2041241452319315*fl[54])*dv; 
  Ghat_l[19] = (1.224744871391589*fl[46]+0.7071067811865475*fl[30])*wv+(0.3535533905932737*fl[62]+0.2041241452319315*fl[55])*dv; 
  Ghat_l[20] = (1.224744871391589*fl[47]+0.7071067811865475*fl[32])*wv+(0.3535533905932737*fl[22]+0.2041241452319315*fl[7])*dv; 
  Ghat_l[21] = (1.224744871391589*fl[49]+0.7071067811865475*fl[35])*wv+(0.3535533905932737*fl[24]+0.2041241452319315*fl[10])*dv; 
  Ghat_l[22] = (1.224744871391589*fl[50]+0.7071067811865475*fl[36])*wv+(0.3535533905932737*fl[25]+0.2041241452319315*fl[11])*dv; 
  Ghat_l[23] = (1.224744871391589*fl[52]+0.7071067811865475*fl[38])*wv+(0.3535533905932737*fl[27]+0.2041241452319315*fl[13])*dv; 
  Ghat_l[24] = (1.224744871391589*fl[53]+0.7071067811865475*fl[39])*wv+(0.3535533905932737*fl[28]+0.2041241452319315*fl[14])*dv; 
  Ghat_l[25] = (1.224744871391589*fl[56]+0.7071067811865475*fl[41])*wv+(0.3535533905932737*fl[31]+0.2041241452319315*fl[16])*dv; 
  Ghat_l[26] = (1.224744871391589*fl[57]+0.7071067811865475*fl[44])*wv+(0.3535533905932737*fl[63]+0.2041241452319315*fl[60])*dv; 
  Ghat_l[27] = (1.224744871391589*fl[58]+0.7071067811865475*fl[48])*wv+(0.3535533905932737*fl[42]+0.2041241452319315*fl[23])*dv; 
  Ghat_l[28] = (1.224744871391589*fl[59]+0.7071067811865475*fl[51])*wv+(0.3535533905932737*fl[43]+0.2041241452319315*fl[26])*dv; 
  Ghat_l[29] = (1.224744871391589*fl[61]+0.7071067811865475*fl[54])*wv+(0.3535533905932737*fl[45]+0.2041241452319315*fl[29])*dv; 
  Ghat_l[30] = (1.224744871391589*fl[62]+0.7071067811865475*fl[55])*wv+(0.3535533905932737*fl[46]+0.2041241452319315*fl[30])*dv; 
  Ghat_l[31] = (1.224744871391589*fl[63]+0.7071067811865475*fl[60])*wv+(0.3535533905932737*fl[57]+0.2041241452319315*fl[44])*dv; 
  Ghat_l[64] = (0.3162277660168379*fl[19]+0.1825741858350554*fl[6])*dv; 
  Ghat_l[65] = (0.3162277660168379*fl[33]+0.1825741858350554*fl[17])*dv; 
  Ghat_l[66] = (0.3162277660168379*fl[34]+0.1825741858350554*fl[18])*dv; 
  Ghat_l[67] = (0.3162277660168379*fl[37]+0.1825741858350554*fl[20])*dv; 
  Ghat_l[68] = (0.3162277660168379*fl[40]+0.1825741858350554*fl[21])*dv; 
  Ghat_l[69] = (0.3162277660168379*fl[47]+0.1825741858350554*fl[32])*dv; 
  Ghat_l[70] = (0.3162277660168379*fl[49]+0.1825741858350554*fl[35])*dv; 
  Ghat_l[71] = (0.3162277660168379*fl[50]+0.1825741858350554*fl[36])*dv; 
  Ghat_l[72] = (0.3162277660168379*fl[52]+0.1825741858350554*fl[38])*dv; 
  Ghat_l[73] = (0.3162277660168379*fl[53]+0.1825741858350554*fl[39])*dv; 
  Ghat_l[74] = (0.3162277660168379*fl[56]+0.1825741858350554*fl[41])*dv; 
  Ghat_l[75] = (0.3162277660168379*fl[58]+0.1825741858350554*fl[48])*dv; 
  Ghat_l[76] = (0.3162277660168379*fl[59]+0.1825741858350554*fl[51])*dv; 
  Ghat_l[77] = (0.3162277660168379*fl[61]+0.1825741858350554*fl[54])*dv; 
  Ghat_l[78] = (0.3162277660168379*fl[62]+0.1825741858350554*fl[55])*dv; 
  Ghat_l[79] = (0.3162277660168379*fl[63]+0.1825741858350554*fl[60])*dv; 

  } else { 

  Ghat_r[0] = -0.1178511301977579*((10.39230484541326*fr[3]-6.0*fr[0])*wv+(3.0*fr[19]-1.732050807568877*fr[6])*dv); 
  Ghat_r[1] = -0.1178511301977579*((10.39230484541326*fr[8]-6.0*fr[1])*wv+(3.0*fr[33]-1.732050807568877*fr[17])*dv); 
  Ghat_r[2] = -0.1178511301977579*((10.39230484541326*fr[9]-6.0*fr[2])*wv+(3.0*fr[34]-1.732050807568877*fr[18])*dv); 
  Ghat_r[3] = -0.1178511301977579*((10.39230484541326*fr[12]-6.0*fr[4])*wv+(3.0*fr[37]-1.732050807568877*fr[20])*dv); 
  Ghat_r[4] = -0.1178511301977579*((10.39230484541326*fr[15]-6.0*fr[5])*wv+(3.0*fr[40]-1.732050807568877*fr[21])*dv); 
  Ghat_r[5] = -0.1178511301977579*((10.39230484541326*fr[19]-6.0*fr[6])*wv+(3.0*fr[3]-1.732050807568877*fr[0])*dv); 
  Ghat_r[6] = -0.1178511301977579*((10.39230484541326*fr[22]-6.0*fr[7])*wv+(3.0*fr[47]-1.732050807568877*fr[32])*dv); 
  Ghat_r[7] = -0.1178511301977579*((10.39230484541326*fr[24]-6.0*fr[10])*wv+(3.0*fr[49]-1.732050807568877*fr[35])*dv); 
  Ghat_r[8] = -0.1178511301977579*((10.39230484541326*fr[25]-6.0*fr[11])*wv+(3.0*fr[50]-1.732050807568877*fr[36])*dv); 
  Ghat_r[9] = -0.1178511301977579*((10.39230484541326*fr[27]-6.0*fr[13])*wv+(3.0*fr[52]-1.732050807568877*fr[38])*dv); 
  Ghat_r[10] = -0.1178511301977579*((10.39230484541326*fr[28]-6.0*fr[14])*wv+(3.0*fr[53]-1.732050807568877*fr[39])*dv); 
  Ghat_r[11] = -0.1178511301977579*((10.39230484541326*fr[31]-6.0*fr[16])*wv+(3.0*fr[56]-1.732050807568877*fr[41])*dv); 
  Ghat_r[12] = -0.1178511301977579*((10.39230484541326*fr[33]-6.0*fr[17])*wv+(3.0*fr[8]-1.732050807568877*fr[1])*dv); 
  Ghat_r[13] = -0.1178511301977579*((10.39230484541326*fr[34]-6.0*fr[18])*wv+(3.0*fr[9]-1.732050807568877*fr[2])*dv); 
  Ghat_r[14] = -0.1178511301977579*((10.39230484541326*fr[37]-6.0*fr[20])*wv+(3.0*fr[12]-1.732050807568877*fr[4])*dv); 
  Ghat_r[15] = -0.1178511301977579*((10.39230484541326*fr[40]-6.0*fr[21])*wv+(3.0*fr[15]-1.732050807568877*fr[5])*dv); 
  Ghat_r[16] = -0.1178511301977579*((10.39230484541326*fr[42]-6.0*fr[23])*wv+(3.0*fr[58]-1.732050807568877*fr[48])*dv); 
  Ghat_r[17] = -0.1178511301977579*((10.39230484541326*fr[43]-6.0*fr[26])*wv+(3.0*fr[59]-1.732050807568877*fr[51])*dv); 
  Ghat_r[18] = -0.1178511301977579*((10.39230484541326*fr[45]-6.0*fr[29])*wv+(3.0*fr[61]-1.732050807568877*fr[54])*dv); 
  Ghat_r[19] = -0.1178511301977579*((10.39230484541326*fr[46]-6.0*fr[30])*wv+(3.0*fr[62]-1.732050807568877*fr[55])*dv); 
  Ghat_r[20] = -0.1178511301977579*((10.39230484541326*fr[47]-6.0*fr[32])*wv+(3.0*fr[22]-1.732050807568877*fr[7])*dv); 
  Ghat_r[21] = -0.1178511301977579*((10.39230484541326*fr[49]-6.0*fr[35])*wv+(3.0*fr[24]-1.732050807568877*fr[10])*dv); 
  Ghat_r[22] = -0.1178511301977579*((10.39230484541326*fr[50]-6.0*fr[36])*wv+(3.0*fr[25]-1.732050807568877*fr[11])*dv); 
  Ghat_r[23] = -0.1178511301977579*((10.39230484541326*fr[52]-6.0*fr[38])*wv+(3.0*fr[27]-1.732050807568877*fr[13])*dv); 
  Ghat_r[24] = -0.1178511301977579*((10.39230484541326*fr[53]-6.0*fr[39])*wv+(3.0*fr[28]-1.732050807568877*fr[14])*dv); 
  Ghat_r[25] = -0.1178511301977579*((10.39230484541326*fr[56]-6.0*fr[41])*wv+(3.0*fr[31]-1.732050807568877*fr[16])*dv); 
  Ghat_r[26] = -0.1178511301977579*((10.39230484541326*fr[57]-6.0*fr[44])*wv+(3.0*fr[63]-1.732050807568877*fr[60])*dv); 
  Ghat_r[27] = -0.1178511301977579*((10.39230484541326*fr[58]-6.0*fr[48])*wv+(3.0*fr[42]-1.732050807568877*fr[23])*dv); 
  Ghat_r[28] = -0.1178511301977579*((10.39230484541326*fr[59]-6.0*fr[51])*wv+(3.0*fr[43]-1.732050807568877*fr[26])*dv); 
  Ghat_r[29] = -0.1178511301977579*((10.39230484541326*fr[61]-6.0*fr[54])*wv+(3.0*fr[45]-1.732050807568877*fr[29])*dv); 
  Ghat_r[30] = -0.1178511301977579*((10.39230484541326*fr[62]-6.0*fr[55])*wv+(3.0*fr[46]-1.732050807568877*fr[30])*dv); 
  Ghat_r[31] = -0.1178511301977579*((10.39230484541326*fr[63]-6.0*fr[60])*wv+(3.0*fr[57]-1.732050807568877*fr[44])*dv); 
  Ghat_r[64] = -0.04714045207910316*(6.708203932499369*fr[19]-3.872983346207417*fr[6])*dv; 
  Ghat_r[65] = -0.04714045207910316*(6.708203932499369*fr[33]-3.872983346207417*fr[17])*dv; 
  Ghat_r[66] = -0.04714045207910316*(6.708203932499369*fr[34]-3.872983346207417*fr[18])*dv; 
  Ghat_r[67] = -0.04714045207910316*(6.708203932499369*fr[37]-3.872983346207417*fr[20])*dv; 
  Ghat_r[68] = -0.04714045207910316*(6.708203932499369*fr[40]-3.872983346207417*fr[21])*dv; 
  Ghat_r[69] = -0.04714045207910316*(6.708203932499369*fr[47]-3.872983346207417*fr[32])*dv; 
  Ghat_r[70] = -0.04714045207910316*(6.708203932499369*fr[49]-3.872983346207417*fr[35])*dv; 
  Ghat_r[71] = -0.04714045207910316*(6.708203932499369*fr[50]-3.872983346207417*fr[36])*dv; 
  Ghat_r[72] = -0.04714045207910316*(6.708203932499369*fr[52]-3.872983346207417*fr[38])*dv; 
  Ghat_r[73] = -0.04714045207910316*(6.708203932499369*fr[53]-3.872983346207417*fr[39])*dv; 
  Ghat_r[74] = -0.04714045207910316*(6.708203932499369*fr[56]-3.872983346207417*fr[41])*dv; 
  Ghat_r[75] = -0.04714045207910316*(6.708203932499369*fr[58]-3.872983346207417*fr[48])*dv; 
  Ghat_r[76] = -0.04714045207910316*(6.708203932499369*fr[59]-3.872983346207417*fr[51])*dv; 
  Ghat_r[77] = -0.04714045207910316*(6.708203932499369*fr[61]-3.872983346207417*fr[54])*dv; 
  Ghat_r[78] = -0.04714045207910316*(6.708203932499369*fr[62]-3.872983346207417*fr[55])*dv; 
  Ghat_r[79] = -0.04714045207910316*(6.708203932499369*fr[63]-3.872983346207417*fr[60])*dv; 

  Ghat_l[0] = -0.1178511301977579*((10.39230484541326*fc[3]-6.0*fc[0])*wv+(3.0*fc[19]-1.732050807568877*fc[6])*dv); 
  Ghat_l[1] = -0.1178511301977579*((10.39230484541326*fc[8]-6.0*fc[1])*wv+(3.0*fc[33]-1.732050807568877*fc[17])*dv); 
  Ghat_l[2] = -0.1178511301977579*((10.39230484541326*fc[9]-6.0*fc[2])*wv+(3.0*fc[34]-1.732050807568877*fc[18])*dv); 
  Ghat_l[3] = -0.1178511301977579*((10.39230484541326*fc[12]-6.0*fc[4])*wv+(3.0*fc[37]-1.732050807568877*fc[20])*dv); 
  Ghat_l[4] = -0.1178511301977579*((10.39230484541326*fc[15]-6.0*fc[5])*wv+(3.0*fc[40]-1.732050807568877*fc[21])*dv); 
  Ghat_l[5] = -0.1178511301977579*((10.39230484541326*fc[19]-6.0*fc[6])*wv+(3.0*fc[3]-1.732050807568877*fc[0])*dv); 
  Ghat_l[6] = -0.1178511301977579*((10.39230484541326*fc[22]-6.0*fc[7])*wv+(3.0*fc[47]-1.732050807568877*fc[32])*dv); 
  Ghat_l[7] = -0.1178511301977579*((10.39230484541326*fc[24]-6.0*fc[10])*wv+(3.0*fc[49]-1.732050807568877*fc[35])*dv); 
  Ghat_l[8] = -0.1178511301977579*((10.39230484541326*fc[25]-6.0*fc[11])*wv+(3.0*fc[50]-1.732050807568877*fc[36])*dv); 
  Ghat_l[9] = -0.1178511301977579*((10.39230484541326*fc[27]-6.0*fc[13])*wv+(3.0*fc[52]-1.732050807568877*fc[38])*dv); 
  Ghat_l[10] = -0.1178511301977579*((10.39230484541326*fc[28]-6.0*fc[14])*wv+(3.0*fc[53]-1.732050807568877*fc[39])*dv); 
  Ghat_l[11] = -0.1178511301977579*((10.39230484541326*fc[31]-6.0*fc[16])*wv+(3.0*fc[56]-1.732050807568877*fc[41])*dv); 
  Ghat_l[12] = -0.1178511301977579*((10.39230484541326*fc[33]-6.0*fc[17])*wv+(3.0*fc[8]-1.732050807568877*fc[1])*dv); 
  Ghat_l[13] = -0.1178511301977579*((10.39230484541326*fc[34]-6.0*fc[18])*wv+(3.0*fc[9]-1.732050807568877*fc[2])*dv); 
  Ghat_l[14] = -0.1178511301977579*((10.39230484541326*fc[37]-6.0*fc[20])*wv+(3.0*fc[12]-1.732050807568877*fc[4])*dv); 
  Ghat_l[15] = -0.1178511301977579*((10.39230484541326*fc[40]-6.0*fc[21])*wv+(3.0*fc[15]-1.732050807568877*fc[5])*dv); 
  Ghat_l[16] = -0.1178511301977579*((10.39230484541326*fc[42]-6.0*fc[23])*wv+(3.0*fc[58]-1.732050807568877*fc[48])*dv); 
  Ghat_l[17] = -0.1178511301977579*((10.39230484541326*fc[43]-6.0*fc[26])*wv+(3.0*fc[59]-1.732050807568877*fc[51])*dv); 
  Ghat_l[18] = -0.1178511301977579*((10.39230484541326*fc[45]-6.0*fc[29])*wv+(3.0*fc[61]-1.732050807568877*fc[54])*dv); 
  Ghat_l[19] = -0.1178511301977579*((10.39230484541326*fc[46]-6.0*fc[30])*wv+(3.0*fc[62]-1.732050807568877*fc[55])*dv); 
  Ghat_l[20] = -0.1178511301977579*((10.39230484541326*fc[47]-6.0*fc[32])*wv+(3.0*fc[22]-1.732050807568877*fc[7])*dv); 
  Ghat_l[21] = -0.1178511301977579*((10.39230484541326*fc[49]-6.0*fc[35])*wv+(3.0*fc[24]-1.732050807568877*fc[10])*dv); 
  Ghat_l[22] = -0.1178511301977579*((10.39230484541326*fc[50]-6.0*fc[36])*wv+(3.0*fc[25]-1.732050807568877*fc[11])*dv); 
  Ghat_l[23] = -0.1178511301977579*((10.39230484541326*fc[52]-6.0*fc[38])*wv+(3.0*fc[27]-1.732050807568877*fc[13])*dv); 
  Ghat_l[24] = -0.1178511301977579*((10.39230484541326*fc[53]-6.0*fc[39])*wv+(3.0*fc[28]-1.732050807568877*fc[14])*dv); 
  Ghat_l[25] = -0.1178511301977579*((10.39230484541326*fc[56]-6.0*fc[41])*wv+(3.0*fc[31]-1.732050807568877*fc[16])*dv); 
  Ghat_l[26] = -0.1178511301977579*((10.39230484541326*fc[57]-6.0*fc[44])*wv+(3.0*fc[63]-1.732050807568877*fc[60])*dv); 
  Ghat_l[27] = -0.1178511301977579*((10.39230484541326*fc[58]-6.0*fc[48])*wv+(3.0*fc[42]-1.732050807568877*fc[23])*dv); 
  Ghat_l[28] = -0.1178511301977579*((10.39230484541326*fc[59]-6.0*fc[51])*wv+(3.0*fc[43]-1.732050807568877*fc[26])*dv); 
  Ghat_l[29] = -0.1178511301977579*((10.39230484541326*fc[61]-6.0*fc[54])*wv+(3.0*fc[45]-1.732050807568877*fc[29])*dv); 
  Ghat_l[30] = -0.1178511301977579*((10.39230484541326*fc[62]-6.0*fc[55])*wv+(3.0*fc[46]-1.732050807568877*fc[30])*dv); 
  Ghat_l[31] = -0.1178511301977579*((10.39230484541326*fc[63]-6.0*fc[60])*wv+(3.0*fc[57]-1.732050807568877*fc[44])*dv); 
  Ghat_l[64] = -0.04714045207910316*(6.708203932499369*fc[19]-3.872983346207417*fc[6])*dv; 
  Ghat_l[65] = -0.04714045207910316*(6.708203932499369*fc[33]-3.872983346207417*fc[17])*dv; 
  Ghat_l[66] = -0.04714045207910316*(6.708203932499369*fc[34]-3.872983346207417*fc[18])*dv; 
  Ghat_l[67] = -0.04714045207910316*(6.708203932499369*fc[37]-3.872983346207417*fc[20])*dv; 
  Ghat_l[68] = -0.04714045207910316*(6.708203932499369*fc[40]-3.872983346207417*fc[21])*dv; 
  Ghat_l[69] = -0.04714045207910316*(6.708203932499369*fc[47]-3.872983346207417*fc[32])*dv; 
  Ghat_l[70] = -0.04714045207910316*(6.708203932499369*fc[49]-3.872983346207417*fc[35])*dv; 
  Ghat_l[71] = -0.04714045207910316*(6.708203932499369*fc[50]-3.872983346207417*fc[36])*dv; 
  Ghat_l[72] = -0.04714045207910316*(6.708203932499369*fc[52]-3.872983346207417*fc[38])*dv; 
  Ghat_l[73] = -0.04714045207910316*(6.708203932499369*fc[53]-3.872983346207417*fc[39])*dv; 
  Ghat_l[74] = -0.04714045207910316*(6.708203932499369*fc[56]-3.872983346207417*fc[41])*dv; 
  Ghat_l[75] = -0.04714045207910316*(6.708203932499369*fc[58]-3.872983346207417*fc[48])*dv; 
  Ghat_l[76] = -0.04714045207910316*(6.708203932499369*fc[59]-3.872983346207417*fc[51])*dv; 
  Ghat_l[77] = -0.04714045207910316*(6.708203932499369*fc[61]-3.872983346207417*fc[54])*dv; 
  Ghat_l[78] = -0.04714045207910316*(6.708203932499369*fc[62]-3.872983346207417*fc[55])*dv; 
  Ghat_l[79] = -0.04714045207910316*(6.708203932499369*fc[63]-3.872983346207417*fc[60])*dv; 

  } 
  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dx12; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dx12; 
  out[2] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dx12; 
  out[3] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dx12; 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dx12; 
  out[5] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dx12; 
  out[6] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dx12; 
  out[7] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dx12; 
  out[8] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dx12; 
  out[9] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dx12; 
  out[10] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dx12; 
  out[11] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*dx12; 
  out[12] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dx12; 
  out[13] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*dx12; 
  out[14] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*dx12; 
  out[15] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dx12; 
  out[16] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*dx12; 
  out[17] += (0.7071067811865475*Ghat_l[12]-0.7071067811865475*Ghat_r[12])*dx12; 
  out[18] += (0.7071067811865475*Ghat_l[13]-0.7071067811865475*Ghat_r[13])*dx12; 
  out[19] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dx12; 
  out[20] += (0.7071067811865475*Ghat_l[14]-0.7071067811865475*Ghat_r[14])*dx12; 
  out[21] += (0.7071067811865475*Ghat_l[15]-0.7071067811865475*Ghat_r[15])*dx12; 
  out[22] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dx12; 
  out[23] += (0.7071067811865475*Ghat_l[16]-0.7071067811865475*Ghat_r[16])*dx12; 
  out[24] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dx12; 
  out[25] += -1.224744871391589*(Ghat_r[8]+Ghat_l[8])*dx12; 
  out[26] += (0.7071067811865475*Ghat_l[17]-0.7071067811865475*Ghat_r[17])*dx12; 
  out[27] += -1.224744871391589*(Ghat_r[9]+Ghat_l[9])*dx12; 
  out[28] += -1.224744871391589*(Ghat_r[10]+Ghat_l[10])*dx12; 
  out[29] += (0.7071067811865475*Ghat_l[18]-0.7071067811865475*Ghat_r[18])*dx12; 
  out[30] += (0.7071067811865475*Ghat_l[19]-0.7071067811865475*Ghat_r[19])*dx12; 
  out[31] += -1.224744871391589*(Ghat_r[11]+Ghat_l[11])*dx12; 
  out[32] += (0.7071067811865475*Ghat_l[20]-0.7071067811865475*Ghat_r[20])*dx12; 
  out[33] += -1.224744871391589*(Ghat_r[12]+Ghat_l[12])*dx12; 
  out[34] += -1.224744871391589*(Ghat_r[13]+Ghat_l[13])*dx12; 
  out[35] += (0.7071067811865475*Ghat_l[21]-0.7071067811865475*Ghat_r[21])*dx12; 
  out[36] += (0.7071067811865475*Ghat_l[22]-0.7071067811865475*Ghat_r[22])*dx12; 
  out[37] += -1.224744871391589*(Ghat_r[14]+Ghat_l[14])*dx12; 
  out[38] += (0.7071067811865475*Ghat_l[23]-0.7071067811865475*Ghat_r[23])*dx12; 
  out[39] += (0.7071067811865475*Ghat_l[24]-0.7071067811865475*Ghat_r[24])*dx12; 
  out[40] += -1.224744871391589*(Ghat_r[15]+Ghat_l[15])*dx12; 
  out[41] += (0.7071067811865475*Ghat_l[25]-0.7071067811865475*Ghat_r[25])*dx12; 
  out[42] += -1.224744871391589*(Ghat_r[16]+Ghat_l[16])*dx12; 
  out[43] += -1.224744871391589*(Ghat_r[17]+Ghat_l[17])*dx12; 
  out[44] += (0.7071067811865475*Ghat_l[26]-0.7071067811865475*Ghat_r[26])*dx12; 
  out[45] += -1.224744871391589*(Ghat_r[18]+Ghat_l[18])*dx12; 
  out[46] += -1.224744871391589*(Ghat_r[19]+Ghat_l[19])*dx12; 
  out[47] += -1.224744871391589*(Ghat_r[20]+Ghat_l[20])*dx12; 
  out[48] += (0.7071067811865475*Ghat_l[27]-0.7071067811865475*Ghat_r[27])*dx12; 
  out[49] += -1.224744871391589*(Ghat_r[21]+Ghat_l[21])*dx12; 
  out[50] += -1.224744871391589*(Ghat_r[22]+Ghat_l[22])*dx12; 
  out[51] += (0.7071067811865475*Ghat_l[28]-0.7071067811865475*Ghat_r[28])*dx12; 
  out[52] += -1.224744871391589*(Ghat_r[23]+Ghat_l[23])*dx12; 
  out[53] += -1.224744871391589*(Ghat_r[24]+Ghat_l[24])*dx12; 
  out[54] += (0.7071067811865475*Ghat_l[29]-0.7071067811865475*Ghat_r[29])*dx12; 
  out[55] += (0.7071067811865475*Ghat_l[30]-0.7071067811865475*Ghat_r[30])*dx12; 
  out[56] += -1.224744871391589*(Ghat_r[25]+Ghat_l[25])*dx12; 
  out[57] += -1.224744871391589*(Ghat_r[26]+Ghat_l[26])*dx12; 
  out[58] += -1.224744871391589*(Ghat_r[27]+Ghat_l[27])*dx12; 
  out[59] += -1.224744871391589*(Ghat_r[28]+Ghat_l[28])*dx12; 
  out[60] += (0.7071067811865475*Ghat_l[31]-0.7071067811865475*Ghat_r[31])*dx12; 
  out[61] += -1.224744871391589*(Ghat_r[29]+Ghat_l[29])*dx12; 
  out[62] += -1.224744871391589*(Ghat_r[30]+Ghat_l[30])*dx12; 
  out[63] += -1.224744871391589*(Ghat_r[31]+Ghat_l[31])*dx12; 

  return 0.;

} 
