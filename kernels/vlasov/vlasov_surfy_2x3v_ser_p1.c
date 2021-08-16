#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_surfy_2x3v_ser_p1(const double *w, const double *dxv, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fl/fc/fr:  Input Distribution function in left/center/right cells.
  // out:       Incremented distribution function in center cell.
  const double dx11 = 2/dxv[1]; 
  const double dv = dxv[3], wv = w[3]; 
  double Ghat_r[16]; 
  double Ghat_l[16]; 
  if (wv>0) { 

  Ghat_r[0] = (1.224744871391589*fc[2]+0.7071067811865475*fc[0])*wv+(0.3535533905932737*fc[10]+0.2041241452319315*fc[4])*dv; 
  Ghat_r[1] = (1.224744871391589*fc[6]+0.7071067811865475*fc[1])*wv+(0.3535533905932737*fc[17]+0.2041241452319315*fc[9])*dv; 
  Ghat_r[2] = (1.224744871391589*fc[8]+0.7071067811865475*fc[3])*wv+(0.3535533905932737*fc[19]+0.2041241452319315*fc[11])*dv; 
  Ghat_r[3] = (1.224744871391589*fc[10]+0.7071067811865475*fc[4])*wv+(0.3535533905932737*fc[2]+0.2041241452319315*fc[0])*dv; 
  Ghat_r[4] = (1.224744871391589*fc[13]+0.7071067811865475*fc[5])*wv+(0.3535533905932737*fc[24]+0.2041241452319315*fc[15])*dv; 
  Ghat_r[5] = (1.224744871391589*fc[16]+0.7071067811865475*fc[7])*wv+(0.3535533905932737*fc[26]+0.2041241452319315*fc[18])*dv; 
  Ghat_r[6] = (1.224744871391589*fc[17]+0.7071067811865475*fc[9])*wv+(0.3535533905932737*fc[6]+0.2041241452319315*fc[1])*dv; 
  Ghat_r[7] = (1.224744871391589*fc[19]+0.7071067811865475*fc[11])*wv+(0.3535533905932737*fc[8]+0.2041241452319315*fc[3])*dv; 
  Ghat_r[8] = (1.224744871391589*fc[20]+0.7071067811865475*fc[12])*wv+(0.3535533905932737*fc[28]+0.2041241452319315*fc[23])*dv; 
  Ghat_r[9] = (1.224744871391589*fc[22]+0.7071067811865475*fc[14])*wv+(0.3535533905932737*fc[30]+0.2041241452319315*fc[25])*dv; 
  Ghat_r[10] = (1.224744871391589*fc[24]+0.7071067811865475*fc[15])*wv+(0.3535533905932737*fc[13]+0.2041241452319315*fc[5])*dv; 
  Ghat_r[11] = (1.224744871391589*fc[26]+0.7071067811865475*fc[18])*wv+(0.3535533905932737*fc[16]+0.2041241452319315*fc[7])*dv; 
  Ghat_r[12] = (1.224744871391589*fc[27]+0.7071067811865475*fc[21])*wv+(0.3535533905932737*fc[31]+0.2041241452319315*fc[29])*dv; 
  Ghat_r[13] = (1.224744871391589*fc[28]+0.7071067811865475*fc[23])*wv+(0.3535533905932737*fc[20]+0.2041241452319315*fc[12])*dv; 
  Ghat_r[14] = (1.224744871391589*fc[30]+0.7071067811865475*fc[25])*wv+(0.3535533905932737*fc[22]+0.2041241452319315*fc[14])*dv; 
  Ghat_r[15] = (1.224744871391589*fc[31]+0.7071067811865475*fc[29])*wv+(0.3535533905932737*fc[27]+0.2041241452319315*fc[21])*dv; 

  Ghat_l[0] = (1.224744871391589*fl[2]+0.7071067811865475*fl[0])*wv+(0.3535533905932737*fl[10]+0.2041241452319315*fl[4])*dv; 
  Ghat_l[1] = (1.224744871391589*fl[6]+0.7071067811865475*fl[1])*wv+(0.3535533905932737*fl[17]+0.2041241452319315*fl[9])*dv; 
  Ghat_l[2] = (1.224744871391589*fl[8]+0.7071067811865475*fl[3])*wv+(0.3535533905932737*fl[19]+0.2041241452319315*fl[11])*dv; 
  Ghat_l[3] = (1.224744871391589*fl[10]+0.7071067811865475*fl[4])*wv+(0.3535533905932737*fl[2]+0.2041241452319315*fl[0])*dv; 
  Ghat_l[4] = (1.224744871391589*fl[13]+0.7071067811865475*fl[5])*wv+(0.3535533905932737*fl[24]+0.2041241452319315*fl[15])*dv; 
  Ghat_l[5] = (1.224744871391589*fl[16]+0.7071067811865475*fl[7])*wv+(0.3535533905932737*fl[26]+0.2041241452319315*fl[18])*dv; 
  Ghat_l[6] = (1.224744871391589*fl[17]+0.7071067811865475*fl[9])*wv+(0.3535533905932737*fl[6]+0.2041241452319315*fl[1])*dv; 
  Ghat_l[7] = (1.224744871391589*fl[19]+0.7071067811865475*fl[11])*wv+(0.3535533905932737*fl[8]+0.2041241452319315*fl[3])*dv; 
  Ghat_l[8] = (1.224744871391589*fl[20]+0.7071067811865475*fl[12])*wv+(0.3535533905932737*fl[28]+0.2041241452319315*fl[23])*dv; 
  Ghat_l[9] = (1.224744871391589*fl[22]+0.7071067811865475*fl[14])*wv+(0.3535533905932737*fl[30]+0.2041241452319315*fl[25])*dv; 
  Ghat_l[10] = (1.224744871391589*fl[24]+0.7071067811865475*fl[15])*wv+(0.3535533905932737*fl[13]+0.2041241452319315*fl[5])*dv; 
  Ghat_l[11] = (1.224744871391589*fl[26]+0.7071067811865475*fl[18])*wv+(0.3535533905932737*fl[16]+0.2041241452319315*fl[7])*dv; 
  Ghat_l[12] = (1.224744871391589*fl[27]+0.7071067811865475*fl[21])*wv+(0.3535533905932737*fl[31]+0.2041241452319315*fl[29])*dv; 
  Ghat_l[13] = (1.224744871391589*fl[28]+0.7071067811865475*fl[23])*wv+(0.3535533905932737*fl[20]+0.2041241452319315*fl[12])*dv; 
  Ghat_l[14] = (1.224744871391589*fl[30]+0.7071067811865475*fl[25])*wv+(0.3535533905932737*fl[22]+0.2041241452319315*fl[14])*dv; 
  Ghat_l[15] = (1.224744871391589*fl[31]+0.7071067811865475*fl[29])*wv+(0.3535533905932737*fl[27]+0.2041241452319315*fl[21])*dv; 

  } else { 

  Ghat_r[0] = (0.7071067811865475*fr[0]-1.224744871391589*fr[2])*wv+(0.2041241452319315*fr[4]-0.3535533905932737*fr[10])*dv; 
  Ghat_r[1] = (0.7071067811865475*fr[1]-1.224744871391589*fr[6])*wv+(0.2041241452319315*fr[9]-0.3535533905932737*fr[17])*dv; 
  Ghat_r[2] = (0.7071067811865475*fr[3]-1.224744871391589*fr[8])*wv+(0.2041241452319315*fr[11]-0.3535533905932737*fr[19])*dv; 
  Ghat_r[3] = (0.7071067811865475*fr[4]-1.224744871391589*fr[10])*wv+(0.2041241452319315*fr[0]-0.3535533905932737*fr[2])*dv; 
  Ghat_r[4] = (0.7071067811865475*fr[5]-1.224744871391589*fr[13])*wv+(0.2041241452319315*fr[15]-0.3535533905932737*fr[24])*dv; 
  Ghat_r[5] = (0.7071067811865475*fr[7]-1.224744871391589*fr[16])*wv+(0.2041241452319315*fr[18]-0.3535533905932737*fr[26])*dv; 
  Ghat_r[6] = (0.7071067811865475*fr[9]-1.224744871391589*fr[17])*wv+(0.2041241452319315*fr[1]-0.3535533905932737*fr[6])*dv; 
  Ghat_r[7] = (0.7071067811865475*fr[11]-1.224744871391589*fr[19])*wv+(0.2041241452319315*fr[3]-0.3535533905932737*fr[8])*dv; 
  Ghat_r[8] = (0.7071067811865475*fr[12]-1.224744871391589*fr[20])*wv+(0.2041241452319315*fr[23]-0.3535533905932737*fr[28])*dv; 
  Ghat_r[9] = (0.7071067811865475*fr[14]-1.224744871391589*fr[22])*wv+(0.2041241452319315*fr[25]-0.3535533905932737*fr[30])*dv; 
  Ghat_r[10] = (0.7071067811865475*fr[15]-1.224744871391589*fr[24])*wv+(0.2041241452319315*fr[5]-0.3535533905932737*fr[13])*dv; 
  Ghat_r[11] = (0.7071067811865475*fr[18]-1.224744871391589*fr[26])*wv+(0.2041241452319315*fr[7]-0.3535533905932737*fr[16])*dv; 
  Ghat_r[12] = (0.7071067811865475*fr[21]-1.224744871391589*fr[27])*wv+(0.2041241452319315*fr[29]-0.3535533905932737*fr[31])*dv; 
  Ghat_r[13] = (0.7071067811865475*fr[23]-1.224744871391589*fr[28])*wv+(0.2041241452319315*fr[12]-0.3535533905932737*fr[20])*dv; 
  Ghat_r[14] = (0.7071067811865475*fr[25]-1.224744871391589*fr[30])*wv+(0.2041241452319315*fr[14]-0.3535533905932737*fr[22])*dv; 
  Ghat_r[15] = (0.7071067811865475*fr[29]-1.224744871391589*fr[31])*wv+(0.2041241452319315*fr[21]-0.3535533905932737*fr[27])*dv; 

  Ghat_l[0] = (0.7071067811865475*fc[0]-1.224744871391589*fc[2])*wv+(0.2041241452319315*fc[4]-0.3535533905932737*fc[10])*dv; 
  Ghat_l[1] = (0.7071067811865475*fc[1]-1.224744871391589*fc[6])*wv+(0.2041241452319315*fc[9]-0.3535533905932737*fc[17])*dv; 
  Ghat_l[2] = (0.7071067811865475*fc[3]-1.224744871391589*fc[8])*wv+(0.2041241452319315*fc[11]-0.3535533905932737*fc[19])*dv; 
  Ghat_l[3] = (0.7071067811865475*fc[4]-1.224744871391589*fc[10])*wv+(0.2041241452319315*fc[0]-0.3535533905932737*fc[2])*dv; 
  Ghat_l[4] = (0.7071067811865475*fc[5]-1.224744871391589*fc[13])*wv+(0.2041241452319315*fc[15]-0.3535533905932737*fc[24])*dv; 
  Ghat_l[5] = (0.7071067811865475*fc[7]-1.224744871391589*fc[16])*wv+(0.2041241452319315*fc[18]-0.3535533905932737*fc[26])*dv; 
  Ghat_l[6] = (0.7071067811865475*fc[9]-1.224744871391589*fc[17])*wv+(0.2041241452319315*fc[1]-0.3535533905932737*fc[6])*dv; 
  Ghat_l[7] = (0.7071067811865475*fc[11]-1.224744871391589*fc[19])*wv+(0.2041241452319315*fc[3]-0.3535533905932737*fc[8])*dv; 
  Ghat_l[8] = (0.7071067811865475*fc[12]-1.224744871391589*fc[20])*wv+(0.2041241452319315*fc[23]-0.3535533905932737*fc[28])*dv; 
  Ghat_l[9] = (0.7071067811865475*fc[14]-1.224744871391589*fc[22])*wv+(0.2041241452319315*fc[25]-0.3535533905932737*fc[30])*dv; 
  Ghat_l[10] = (0.7071067811865475*fc[15]-1.224744871391589*fc[24])*wv+(0.2041241452319315*fc[5]-0.3535533905932737*fc[13])*dv; 
  Ghat_l[11] = (0.7071067811865475*fc[18]-1.224744871391589*fc[26])*wv+(0.2041241452319315*fc[7]-0.3535533905932737*fc[16])*dv; 
  Ghat_l[12] = (0.7071067811865475*fc[21]-1.224744871391589*fc[27])*wv+(0.2041241452319315*fc[29]-0.3535533905932737*fc[31])*dv; 
  Ghat_l[13] = (0.7071067811865475*fc[23]-1.224744871391589*fc[28])*wv+(0.2041241452319315*fc[12]-0.3535533905932737*fc[20])*dv; 
  Ghat_l[14] = (0.7071067811865475*fc[25]-1.224744871391589*fc[30])*wv+(0.2041241452319315*fc[14]-0.3535533905932737*fc[22])*dv; 
  Ghat_l[15] = (0.7071067811865475*fc[29]-1.224744871391589*fc[31])*wv+(0.2041241452319315*fc[21]-0.3535533905932737*fc[27])*dv; 

  } 
  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dx11; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dx11; 
  out[2] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dx11; 
  out[3] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dx11; 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dx11; 
  out[5] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dx11; 
  out[6] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dx11; 
  out[7] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dx11; 
  out[8] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dx11; 
  out[9] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dx11; 
  out[10] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dx11; 
  out[11] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dx11; 
  out[12] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*dx11; 
  out[13] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dx11; 
  out[14] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*dx11; 
  out[15] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*dx11; 
  out[16] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dx11; 
  out[17] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dx11; 
  out[18] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*dx11; 
  out[19] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dx11; 
  out[20] += -1.224744871391589*(Ghat_r[8]+Ghat_l[8])*dx11; 
  out[21] += (0.7071067811865475*Ghat_l[12]-0.7071067811865475*Ghat_r[12])*dx11; 
  out[22] += -1.224744871391589*(Ghat_r[9]+Ghat_l[9])*dx11; 
  out[23] += (0.7071067811865475*Ghat_l[13]-0.7071067811865475*Ghat_r[13])*dx11; 
  out[24] += -1.224744871391589*(Ghat_r[10]+Ghat_l[10])*dx11; 
  out[25] += (0.7071067811865475*Ghat_l[14]-0.7071067811865475*Ghat_r[14])*dx11; 
  out[26] += -1.224744871391589*(Ghat_r[11]+Ghat_l[11])*dx11; 
  out[27] += -1.224744871391589*(Ghat_r[12]+Ghat_l[12])*dx11; 
  out[28] += -1.224744871391589*(Ghat_r[13]+Ghat_l[13])*dx11; 
  out[29] += (0.7071067811865475*Ghat_l[15]-0.7071067811865475*Ghat_r[15])*dx11; 
  out[30] += -1.224744871391589*(Ghat_r[14]+Ghat_l[14])*dx11; 
  out[31] += -1.224744871391589*(Ghat_r[15]+Ghat_l[15])*dx11; 
} 