#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_surfx_2x2v_ser_p2(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // alpha_geo: Fields used only for general geometry.
  // fl/fc/fr:  Input Distribution function in left/center/right cells.
  // out:       Incremented distribution function in center cell.
  const double dx10 = 2/dxv[0]; 
  const double dv = dxv[2], wv = w[2]; 
  double Ghat_r[20]; 
  double Ghat_l[20]; 
  if (wv>0) { 

  Ghat_r[0] = (1.58113883008419*fc[11]+1.224744871391589*fc[1]+0.7071067811865475*fc[0])*wv+(0.4564354645876384*fc[21]+0.3535533905932737*fc[6]+0.2041241452319315*fc[3])*dv; 
  Ghat_r[1] = (1.58113883008419*fc[19]+1.224744871391589*fc[5]+0.7071067811865475*fc[2])*wv+(0.4564354645876384*fc[32]+0.3535533905932737*fc[15]+0.2041241452319315*fc[7])*dv; 
  Ghat_r[2] = (1.58113883008419*fc[21]+1.224744871391589*fc[6]+0.7071067811865475*fc[3])*wv+(0.3162277660168379*fc[23]+0.1825741858350554*fc[13]+0.4564354645876384*fc[11]+0.3535533905932737*fc[1]+0.2041241452319315*fc[0])*dv; 
  Ghat_r[3] = (1.58113883008419*fc[25]+1.224744871391589*fc[8]+0.7071067811865475*fc[4])*wv+(0.4564354645876384*fc[37]+0.3535533905932737*fc[17]+0.2041241452319315*fc[10])*dv; 
  Ghat_r[4] = (1.58113883008419*fc[32]+1.224744871391589*fc[15]+0.7071067811865475*fc[7])*wv+(0.3162277660168379*fc[34]+0.1825741858350554*fc[24]+0.4564354645876384*fc[19]+0.3535533905932737*fc[5]+0.2041241452319315*fc[2])*dv; 
  Ghat_r[5] = (1.58113883008419*fc[35]+1.224744871391589*fc[16]+0.7071067811865475*fc[9])*wv+(0.4564354645876384*fc[44]+0.3535533905932737*fc[31]+0.2041241452319315*fc[18])*dv; 
  Ghat_r[6] = (1.58113883008419*fc[37]+1.224744871391589*fc[17]+0.7071067811865475*fc[10])*wv+(0.3162277660168379*fc[39]+0.1825741858350554*fc[27]+0.4564354645876384*fc[25]+0.3535533905932737*fc[8]+0.2041241452319315*fc[4])*dv; 
  Ghat_r[7] = (1.224744871391589*fc[20]+0.7071067811865475*fc[12])*wv+(0.3535533905932737*fc[33]+0.2041241452319315*fc[22])*dv; 
  Ghat_r[8] = (1.224744871391589*fc[23]+0.7071067811865475*fc[13])*wv+(0.408248290463863*fc[21]+0.3162277660168379*fc[6]+0.1825741858350554*fc[3])*dv; 
  Ghat_r[9] = (1.224744871391589*fc[28]+0.7071067811865475*fc[14])*wv+(0.3535533905932737*fc[42]+0.2041241452319315*fc[30])*dv; 
  Ghat_r[10] = (1.58113883008419*fc[44]+1.224744871391589*fc[31]+0.7071067811865475*fc[18])*wv+(0.3162277660168379*fc[46]+0.1825741858350554*fc[40]+0.4564354645876384*fc[35]+0.3535533905932737*fc[16]+0.2041241452319315*fc[9])*dv; 
  Ghat_r[11] = (1.224744871391589*fc[33]+0.7071067811865475*fc[22])*wv+(0.3535533905932737*fc[20]+0.2041241452319315*fc[12])*dv; 
  Ghat_r[12] = (1.224744871391589*fc[34]+0.7071067811865475*fc[24])*wv+(0.408248290463863*fc[32]+0.3162277660168379*fc[15]+0.1825741858350554*fc[7])*dv; 
  Ghat_r[13] = (1.224744871391589*fc[36]+0.7071067811865475*fc[26])*wv+(0.3535533905932737*fc[45]+0.2041241452319315*fc[38])*dv; 
  Ghat_r[14] = (1.224744871391589*fc[39]+0.7071067811865475*fc[27])*wv+(0.408248290463863*fc[37]+0.3162277660168379*fc[17]+0.1825741858350554*fc[10])*dv; 
  Ghat_r[15] = (1.224744871391589*fc[41]+0.7071067811865475*fc[29])*wv+(0.3535533905932737*fc[47]+0.2041241452319315*fc[43])*dv; 
  Ghat_r[16] = (1.224744871391589*fc[42]+0.7071067811865475*fc[30])*wv+(0.3535533905932737*fc[28]+0.2041241452319315*fc[14])*dv; 
  Ghat_r[17] = (1.224744871391589*fc[45]+0.7071067811865475*fc[38])*wv+(0.3535533905932737*fc[36]+0.2041241452319315*fc[26])*dv; 
  Ghat_r[18] = (1.224744871391589*fc[46]+0.7071067811865475*fc[40])*wv+(0.408248290463863*fc[44]+0.3162277660168379*fc[31]+0.1825741858350554*fc[18])*dv; 
  Ghat_r[19] = (1.224744871391589*fc[47]+0.7071067811865475*fc[43])*wv+(0.3535533905932737*fc[41]+0.2041241452319315*fc[29])*dv; 

  Ghat_l[0] = (1.58113883008419*fl[11]+1.224744871391589*fl[1]+0.7071067811865475*fl[0])*wv+(0.4564354645876384*fl[21]+0.3535533905932737*fl[6]+0.2041241452319315*fl[3])*dv; 
  Ghat_l[1] = (1.58113883008419*fl[19]+1.224744871391589*fl[5]+0.7071067811865475*fl[2])*wv+(0.4564354645876384*fl[32]+0.3535533905932737*fl[15]+0.2041241452319315*fl[7])*dv; 
  Ghat_l[2] = (1.58113883008419*fl[21]+1.224744871391589*fl[6]+0.7071067811865475*fl[3])*wv+(0.3162277660168379*fl[23]+0.1825741858350554*fl[13]+0.4564354645876384*fl[11]+0.3535533905932737*fl[1]+0.2041241452319315*fl[0])*dv; 
  Ghat_l[3] = (1.58113883008419*fl[25]+1.224744871391589*fl[8]+0.7071067811865475*fl[4])*wv+(0.4564354645876384*fl[37]+0.3535533905932737*fl[17]+0.2041241452319315*fl[10])*dv; 
  Ghat_l[4] = (1.58113883008419*fl[32]+1.224744871391589*fl[15]+0.7071067811865475*fl[7])*wv+(0.3162277660168379*fl[34]+0.1825741858350554*fl[24]+0.4564354645876384*fl[19]+0.3535533905932737*fl[5]+0.2041241452319315*fl[2])*dv; 
  Ghat_l[5] = (1.58113883008419*fl[35]+1.224744871391589*fl[16]+0.7071067811865475*fl[9])*wv+(0.4564354645876384*fl[44]+0.3535533905932737*fl[31]+0.2041241452319315*fl[18])*dv; 
  Ghat_l[6] = (1.58113883008419*fl[37]+1.224744871391589*fl[17]+0.7071067811865475*fl[10])*wv+(0.3162277660168379*fl[39]+0.1825741858350554*fl[27]+0.4564354645876384*fl[25]+0.3535533905932737*fl[8]+0.2041241452319315*fl[4])*dv; 
  Ghat_l[7] = (1.224744871391589*fl[20]+0.7071067811865475*fl[12])*wv+(0.3535533905932737*fl[33]+0.2041241452319315*fl[22])*dv; 
  Ghat_l[8] = (1.224744871391589*fl[23]+0.7071067811865475*fl[13])*wv+(0.408248290463863*fl[21]+0.3162277660168379*fl[6]+0.1825741858350554*fl[3])*dv; 
  Ghat_l[9] = (1.224744871391589*fl[28]+0.7071067811865475*fl[14])*wv+(0.3535533905932737*fl[42]+0.2041241452319315*fl[30])*dv; 
  Ghat_l[10] = (1.58113883008419*fl[44]+1.224744871391589*fl[31]+0.7071067811865475*fl[18])*wv+(0.3162277660168379*fl[46]+0.1825741858350554*fl[40]+0.4564354645876384*fl[35]+0.3535533905932737*fl[16]+0.2041241452319315*fl[9])*dv; 
  Ghat_l[11] = (1.224744871391589*fl[33]+0.7071067811865475*fl[22])*wv+(0.3535533905932737*fl[20]+0.2041241452319315*fl[12])*dv; 
  Ghat_l[12] = (1.224744871391589*fl[34]+0.7071067811865475*fl[24])*wv+(0.408248290463863*fl[32]+0.3162277660168379*fl[15]+0.1825741858350554*fl[7])*dv; 
  Ghat_l[13] = (1.224744871391589*fl[36]+0.7071067811865475*fl[26])*wv+(0.3535533905932737*fl[45]+0.2041241452319315*fl[38])*dv; 
  Ghat_l[14] = (1.224744871391589*fl[39]+0.7071067811865475*fl[27])*wv+(0.408248290463863*fl[37]+0.3162277660168379*fl[17]+0.1825741858350554*fl[10])*dv; 
  Ghat_l[15] = (1.224744871391589*fl[41]+0.7071067811865475*fl[29])*wv+(0.3535533905932737*fl[47]+0.2041241452319315*fl[43])*dv; 
  Ghat_l[16] = (1.224744871391589*fl[42]+0.7071067811865475*fl[30])*wv+(0.3535533905932737*fl[28]+0.2041241452319315*fl[14])*dv; 
  Ghat_l[17] = (1.224744871391589*fl[45]+0.7071067811865475*fl[38])*wv+(0.3535533905932737*fl[36]+0.2041241452319315*fl[26])*dv; 
  Ghat_l[18] = (1.224744871391589*fl[46]+0.7071067811865475*fl[40])*wv+(0.408248290463863*fl[44]+0.3162277660168379*fl[31]+0.1825741858350554*fl[18])*dv; 
  Ghat_l[19] = (1.224744871391589*fl[47]+0.7071067811865475*fl[43])*wv+(0.3535533905932737*fl[41]+0.2041241452319315*fl[29])*dv; 

  } else { 

  Ghat_r[0] = 1.58113883008419*fr[11]*wv-1.224744871391589*fr[1]*wv+0.7071067811865475*fr[0]*wv+0.4564354645876383*fr[21]*dv-0.3535533905932737*fr[6]*dv+0.2041241452319315*fr[3]*dv; 
  Ghat_r[1] = 1.58113883008419*fr[19]*wv-1.224744871391589*fr[5]*wv+0.7071067811865475*fr[2]*wv+0.4564354645876384*fr[32]*dv-0.3535533905932737*fr[15]*dv+0.2041241452319315*fr[7]*dv; 
  Ghat_r[2] = 1.58113883008419*fr[21]*wv-1.224744871391589*fr[6]*wv+0.7071067811865475*fr[3]*wv-0.3162277660168379*fr[23]*dv+0.1825741858350554*fr[13]*dv+0.4564354645876384*fr[11]*dv-0.3535533905932737*fr[1]*dv+0.2041241452319315*fr[0]*dv; 
  Ghat_r[3] = 1.58113883008419*fr[25]*wv-1.224744871391589*fr[8]*wv+0.7071067811865475*fr[4]*wv+0.4564354645876384*fr[37]*dv-0.3535533905932737*fr[17]*dv+0.2041241452319315*fr[10]*dv; 
  Ghat_r[4] = 1.58113883008419*fr[32]*wv-1.224744871391589*fr[15]*wv+0.7071067811865475*fr[7]*wv-0.3162277660168379*fr[34]*dv+0.1825741858350553*fr[24]*dv+0.4564354645876383*fr[19]*dv-0.3535533905932737*fr[5]*dv+0.2041241452319315*fr[2]*dv; 
  Ghat_r[5] = 1.58113883008419*fr[35]*wv-1.224744871391589*fr[16]*wv+0.7071067811865475*fr[9]*wv+0.4564354645876383*fr[44]*dv-0.3535533905932737*fr[31]*dv+0.2041241452319315*fr[18]*dv; 
  Ghat_r[6] = 1.58113883008419*fr[37]*wv-1.224744871391589*fr[17]*wv+0.7071067811865475*fr[10]*wv-0.3162277660168379*fr[39]*dv+0.1825741858350553*fr[27]*dv+0.4564354645876383*fr[25]*dv-0.3535533905932737*fr[8]*dv+0.2041241452319315*fr[4]*dv; 
  Ghat_r[7] = (-1.224744871391589*fr[20]*wv)+0.7071067811865475*fr[12]*wv-0.3535533905932737*fr[33]*dv+0.2041241452319315*fr[22]*dv; 
  Ghat_r[8] = (-1.224744871391589*fr[23]*wv)+0.7071067811865475*fr[13]*wv+0.408248290463863*fr[21]*dv-0.3162277660168379*fr[6]*dv+0.1825741858350554*fr[3]*dv; 
  Ghat_r[9] = (-1.224744871391589*fr[28]*wv)+0.7071067811865475*fr[14]*wv-0.3535533905932737*fr[42]*dv+0.2041241452319315*fr[30]*dv; 
  Ghat_r[10] = 1.58113883008419*fr[44]*wv-1.224744871391589*fr[31]*wv+0.7071067811865475*fr[18]*wv-0.3162277660168379*fr[46]*dv+0.1825741858350554*fr[40]*dv+0.4564354645876384*fr[35]*dv-0.3535533905932737*fr[16]*dv+0.2041241452319315*fr[9]*dv; 
  Ghat_r[11] = (-1.224744871391589*fr[33]*wv)+0.7071067811865475*fr[22]*wv-0.3535533905932737*fr[20]*dv+0.2041241452319315*fr[12]*dv; 
  Ghat_r[12] = (-1.224744871391589*fr[34]*wv)+0.7071067811865475*fr[24]*wv+0.408248290463863*fr[32]*dv-0.3162277660168379*fr[15]*dv+0.1825741858350553*fr[7]*dv; 
  Ghat_r[13] = (-1.224744871391589*fr[36]*wv)+0.7071067811865475*fr[26]*wv-0.3535533905932737*fr[45]*dv+0.2041241452319315*fr[38]*dv; 
  Ghat_r[14] = (-1.224744871391589*fr[39]*wv)+0.7071067811865475*fr[27]*wv+0.408248290463863*fr[37]*dv-0.3162277660168379*fr[17]*dv+0.1825741858350553*fr[10]*dv; 
  Ghat_r[15] = (-1.224744871391589*fr[41]*wv)+0.7071067811865475*fr[29]*wv-0.3535533905932737*fr[47]*dv+0.2041241452319315*fr[43]*dv; 
  Ghat_r[16] = (-1.224744871391589*fr[42]*wv)+0.7071067811865475*fr[30]*wv-0.3535533905932737*fr[28]*dv+0.2041241452319315*fr[14]*dv; 
  Ghat_r[17] = (-1.224744871391589*fr[45]*wv)+0.7071067811865475*fr[38]*wv-0.3535533905932737*fr[36]*dv+0.2041241452319315*fr[26]*dv; 
  Ghat_r[18] = (-1.224744871391589*fr[46]*wv)+0.7071067811865475*fr[40]*wv+0.408248290463863*fr[44]*dv-0.3162277660168379*fr[31]*dv+0.1825741858350554*fr[18]*dv; 
  Ghat_r[19] = (-1.224744871391589*fr[47]*wv)+0.7071067811865475*fr[43]*wv-0.3535533905932737*fr[41]*dv+0.2041241452319315*fr[29]*dv; 

  Ghat_l[0] = 1.58113883008419*fc[11]*wv-1.224744871391589*fc[1]*wv+0.7071067811865475*fc[0]*wv+0.4564354645876383*fc[21]*dv-0.3535533905932737*fc[6]*dv+0.2041241452319315*fc[3]*dv; 
  Ghat_l[1] = 1.58113883008419*fc[19]*wv-1.224744871391589*fc[5]*wv+0.7071067811865475*fc[2]*wv+0.4564354645876384*fc[32]*dv-0.3535533905932737*fc[15]*dv+0.2041241452319315*fc[7]*dv; 
  Ghat_l[2] = 1.58113883008419*fc[21]*wv-1.224744871391589*fc[6]*wv+0.7071067811865475*fc[3]*wv-0.3162277660168379*fc[23]*dv+0.1825741858350554*fc[13]*dv+0.4564354645876384*fc[11]*dv-0.3535533905932737*fc[1]*dv+0.2041241452319315*fc[0]*dv; 
  Ghat_l[3] = 1.58113883008419*fc[25]*wv-1.224744871391589*fc[8]*wv+0.7071067811865475*fc[4]*wv+0.4564354645876384*fc[37]*dv-0.3535533905932737*fc[17]*dv+0.2041241452319315*fc[10]*dv; 
  Ghat_l[4] = 1.58113883008419*fc[32]*wv-1.224744871391589*fc[15]*wv+0.7071067811865475*fc[7]*wv-0.3162277660168379*fc[34]*dv+0.1825741858350553*fc[24]*dv+0.4564354645876383*fc[19]*dv-0.3535533905932737*fc[5]*dv+0.2041241452319315*fc[2]*dv; 
  Ghat_l[5] = 1.58113883008419*fc[35]*wv-1.224744871391589*fc[16]*wv+0.7071067811865475*fc[9]*wv+0.4564354645876383*fc[44]*dv-0.3535533905932737*fc[31]*dv+0.2041241452319315*fc[18]*dv; 
  Ghat_l[6] = 1.58113883008419*fc[37]*wv-1.224744871391589*fc[17]*wv+0.7071067811865475*fc[10]*wv-0.3162277660168379*fc[39]*dv+0.1825741858350553*fc[27]*dv+0.4564354645876383*fc[25]*dv-0.3535533905932737*fc[8]*dv+0.2041241452319315*fc[4]*dv; 
  Ghat_l[7] = (-1.224744871391589*fc[20]*wv)+0.7071067811865475*fc[12]*wv-0.3535533905932737*fc[33]*dv+0.2041241452319315*fc[22]*dv; 
  Ghat_l[8] = (-1.224744871391589*fc[23]*wv)+0.7071067811865475*fc[13]*wv+0.408248290463863*fc[21]*dv-0.3162277660168379*fc[6]*dv+0.1825741858350554*fc[3]*dv; 
  Ghat_l[9] = (-1.224744871391589*fc[28]*wv)+0.7071067811865475*fc[14]*wv-0.3535533905932737*fc[42]*dv+0.2041241452319315*fc[30]*dv; 
  Ghat_l[10] = 1.58113883008419*fc[44]*wv-1.224744871391589*fc[31]*wv+0.7071067811865475*fc[18]*wv-0.3162277660168379*fc[46]*dv+0.1825741858350554*fc[40]*dv+0.4564354645876384*fc[35]*dv-0.3535533905932737*fc[16]*dv+0.2041241452319315*fc[9]*dv; 
  Ghat_l[11] = (-1.224744871391589*fc[33]*wv)+0.7071067811865475*fc[22]*wv-0.3535533905932737*fc[20]*dv+0.2041241452319315*fc[12]*dv; 
  Ghat_l[12] = (-1.224744871391589*fc[34]*wv)+0.7071067811865475*fc[24]*wv+0.408248290463863*fc[32]*dv-0.3162277660168379*fc[15]*dv+0.1825741858350553*fc[7]*dv; 
  Ghat_l[13] = (-1.224744871391589*fc[36]*wv)+0.7071067811865475*fc[26]*wv-0.3535533905932737*fc[45]*dv+0.2041241452319315*fc[38]*dv; 
  Ghat_l[14] = (-1.224744871391589*fc[39]*wv)+0.7071067811865475*fc[27]*wv+0.408248290463863*fc[37]*dv-0.3162277660168379*fc[17]*dv+0.1825741858350553*fc[10]*dv; 
  Ghat_l[15] = (-1.224744871391589*fc[41]*wv)+0.7071067811865475*fc[29]*wv-0.3535533905932737*fc[47]*dv+0.2041241452319315*fc[43]*dv; 
  Ghat_l[16] = (-1.224744871391589*fc[42]*wv)+0.7071067811865475*fc[30]*wv-0.3535533905932737*fc[28]*dv+0.2041241452319315*fc[14]*dv; 
  Ghat_l[17] = (-1.224744871391589*fc[45]*wv)+0.7071067811865475*fc[38]*wv-0.3535533905932737*fc[36]*dv+0.2041241452319315*fc[26]*dv; 
  Ghat_l[18] = (-1.224744871391589*fc[46]*wv)+0.7071067811865475*fc[40]*wv+0.408248290463863*fc[44]*dv-0.3162277660168379*fc[31]*dv+0.1825741858350554*fc[18]*dv; 
  Ghat_l[19] = (-1.224744871391589*fc[47]*wv)+0.7071067811865475*fc[43]*wv-0.3535533905932737*fc[41]*dv+0.2041241452319315*fc[29]*dv; 

  } 
  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dx10; 
  out[1] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dx10; 
  out[2] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dx10; 
  out[3] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dx10; 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dx10; 
  out[5] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dx10; 
  out[6] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dx10; 
  out[7] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dx10; 
  out[8] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dx10; 
  out[9] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dx10; 
  out[10] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dx10; 
  out[11] += (1.58113883008419*Ghat_l[0]-1.58113883008419*Ghat_r[0])*dx10; 
  out[12] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dx10; 
  out[13] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*dx10; 
  out[14] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*dx10; 
  out[15] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dx10; 
  out[16] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dx10; 
  out[17] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dx10; 
  out[18] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*dx10; 
  out[19] += (1.58113883008419*Ghat_l[1]-1.58113883008419*Ghat_r[1])*dx10; 
  out[20] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dx10; 
  out[21] += (1.58113883008419*Ghat_l[2]-1.58113883008419*Ghat_r[2])*dx10; 
  out[22] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*dx10; 
  out[23] += -1.224744871391589*(Ghat_r[8]+Ghat_l[8])*dx10; 
  out[24] += (0.7071067811865475*Ghat_l[12]-0.7071067811865475*Ghat_r[12])*dx10; 
  out[25] += (1.58113883008419*Ghat_l[3]-1.58113883008419*Ghat_r[3])*dx10; 
  out[26] += (0.7071067811865475*Ghat_l[13]-0.7071067811865475*Ghat_r[13])*dx10; 
  out[27] += (0.7071067811865475*Ghat_l[14]-0.7071067811865475*Ghat_r[14])*dx10; 
  out[28] += -1.224744871391589*(Ghat_r[9]+Ghat_l[9])*dx10; 
  out[29] += (0.7071067811865475*Ghat_l[15]-0.7071067811865475*Ghat_r[15])*dx10; 
  out[30] += (0.7071067811865475*Ghat_l[16]-0.7071067811865475*Ghat_r[16])*dx10; 
  out[31] += -1.224744871391589*(Ghat_r[10]+Ghat_l[10])*dx10; 
  out[32] += (1.58113883008419*Ghat_l[4]-1.58113883008419*Ghat_r[4])*dx10; 
  out[33] += -1.224744871391589*(Ghat_r[11]+Ghat_l[11])*dx10; 
  out[34] += -1.224744871391589*(Ghat_r[12]+Ghat_l[12])*dx10; 
  out[35] += (1.58113883008419*Ghat_l[5]-1.58113883008419*Ghat_r[5])*dx10; 
  out[36] += -1.224744871391589*(Ghat_r[13]+Ghat_l[13])*dx10; 
  out[37] += (1.58113883008419*Ghat_l[6]-1.58113883008419*Ghat_r[6])*dx10; 
  out[38] += (0.7071067811865475*Ghat_l[17]-0.7071067811865475*Ghat_r[17])*dx10; 
  out[39] += -1.224744871391589*(Ghat_r[14]+Ghat_l[14])*dx10; 
  out[40] += (0.7071067811865475*Ghat_l[18]-0.7071067811865475*Ghat_r[18])*dx10; 
  out[41] += -1.224744871391589*(Ghat_r[15]+Ghat_l[15])*dx10; 
  out[42] += -1.224744871391589*(Ghat_r[16]+Ghat_l[16])*dx10; 
  out[43] += (0.7071067811865475*Ghat_l[19]-0.7071067811865475*Ghat_r[19])*dx10; 
  out[44] += (1.58113883008419*Ghat_l[10]-1.58113883008419*Ghat_r[10])*dx10; 
  out[45] += -1.224744871391589*(Ghat_r[17]+Ghat_l[17])*dx10; 
  out[46] += -1.224744871391589*(Ghat_r[18]+Ghat_l[18])*dx10; 
  out[47] += -1.224744871391589*(Ghat_r[19]+Ghat_l[19])*dx10; 

  return 0.;

} 
