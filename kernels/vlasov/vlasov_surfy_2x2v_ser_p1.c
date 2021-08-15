#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_surfy_2x2v_ser_p1(const double *w, const double *dxv, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fl/fc/fr:  Input Distribution function in left/center/right cells.
  // out:       Incremented distribution function in center cell.
  const double dx11 = 2/dxv[1]; 
  const double dv = dxv[3], wv = w[3]; 
  double Ghat_r[8]; 
  double Ghat_l[8]; 
  if (wv>0) { 

  Ghat_r[0] = (1.224744871391589*fc[2]+0.7071067811865475*fc[0])*wv+(0.3535533905932737*fc[9]+0.2041241452319315*fc[4])*dv; 
  Ghat_r[1] = (1.224744871391589*fc[5]+0.7071067811865475*fc[1])*wv+(0.3535533905932737*fc[12]+0.2041241452319315*fc[8])*dv; 
  Ghat_r[2] = (1.224744871391589*fc[7]+0.7071067811865475*fc[3])*wv+(0.3535533905932737*fc[14]+0.2041241452319315*fc[10])*dv; 
  Ghat_r[3] = (1.224744871391589*fc[9]+0.7071067811865475*fc[4])*wv+(0.3535533905932737*fc[2]+0.2041241452319315*fc[0])*dv; 
  Ghat_r[4] = (1.224744871391589*fc[11]+0.7071067811865475*fc[6])*wv+(0.3535533905932737*fc[15]+0.2041241452319315*fc[13])*dv; 
  Ghat_r[5] = (1.224744871391589*fc[12]+0.7071067811865475*fc[8])*wv+(0.3535533905932737*fc[5]+0.2041241452319315*fc[1])*dv; 
  Ghat_r[6] = (1.224744871391589*fc[14]+0.7071067811865475*fc[10])*wv+(0.3535533905932737*fc[7]+0.2041241452319315*fc[3])*dv; 
  Ghat_r[7] = (1.224744871391589*fc[15]+0.7071067811865475*fc[13])*wv+(0.3535533905932737*fc[11]+0.2041241452319315*fc[6])*dv; 

  Ghat_l[0] = (1.224744871391589*fl[2]+0.7071067811865475*fl[0])*wv+(0.3535533905932737*fl[9]+0.2041241452319315*fl[4])*dv; 
  Ghat_l[1] = (1.224744871391589*fl[5]+0.7071067811865475*fl[1])*wv+(0.3535533905932737*fl[12]+0.2041241452319315*fl[8])*dv; 
  Ghat_l[2] = (1.224744871391589*fl[7]+0.7071067811865475*fl[3])*wv+(0.3535533905932737*fl[14]+0.2041241452319315*fl[10])*dv; 
  Ghat_l[3] = (1.224744871391589*fl[9]+0.7071067811865475*fl[4])*wv+(0.3535533905932737*fl[2]+0.2041241452319315*fl[0])*dv; 
  Ghat_l[4] = (1.224744871391589*fl[11]+0.7071067811865475*fl[6])*wv+(0.3535533905932737*fl[15]+0.2041241452319315*fl[13])*dv; 
  Ghat_l[5] = (1.224744871391589*fl[12]+0.7071067811865475*fl[8])*wv+(0.3535533905932737*fl[5]+0.2041241452319315*fl[1])*dv; 
  Ghat_l[6] = (1.224744871391589*fl[14]+0.7071067811865475*fl[10])*wv+(0.3535533905932737*fl[7]+0.2041241452319315*fl[3])*dv; 
  Ghat_l[7] = (1.224744871391589*fl[15]+0.7071067811865475*fl[13])*wv+(0.3535533905932737*fl[11]+0.2041241452319315*fl[6])*dv; 

  } else { 

  Ghat_r[0] = (0.7071067811865475*fr[0]-1.224744871391589*fr[2])*wv+(0.2041241452319315*fr[4]-0.3535533905932737*fr[9])*dv; 
  Ghat_r[1] = (0.7071067811865475*fr[1]-1.224744871391589*fr[5])*wv+(0.2041241452319315*fr[8]-0.3535533905932737*fr[12])*dv; 
  Ghat_r[2] = (0.7071067811865475*fr[3]-1.224744871391589*fr[7])*wv+(0.2041241452319315*fr[10]-0.3535533905932737*fr[14])*dv; 
  Ghat_r[3] = (0.7071067811865475*fr[4]-1.224744871391589*fr[9])*wv+(0.2041241452319315*fr[0]-0.3535533905932737*fr[2])*dv; 
  Ghat_r[4] = (0.7071067811865475*fr[6]-1.224744871391589*fr[11])*wv+(0.2041241452319315*fr[13]-0.3535533905932737*fr[15])*dv; 
  Ghat_r[5] = (0.7071067811865475*fr[8]-1.224744871391589*fr[12])*wv+(0.2041241452319315*fr[1]-0.3535533905932737*fr[5])*dv; 
  Ghat_r[6] = (0.7071067811865475*fr[10]-1.224744871391589*fr[14])*wv+(0.2041241452319315*fr[3]-0.3535533905932737*fr[7])*dv; 
  Ghat_r[7] = (0.7071067811865475*fr[13]-1.224744871391589*fr[15])*wv+(0.2041241452319315*fr[6]-0.3535533905932737*fr[11])*dv; 

  Ghat_l[0] = (0.7071067811865475*fc[0]-1.224744871391589*fc[2])*wv+(0.2041241452319315*fc[4]-0.3535533905932737*fc[9])*dv; 
  Ghat_l[1] = (0.7071067811865475*fc[1]-1.224744871391589*fc[5])*wv+(0.2041241452319315*fc[8]-0.3535533905932737*fc[12])*dv; 
  Ghat_l[2] = (0.7071067811865475*fc[3]-1.224744871391589*fc[7])*wv+(0.2041241452319315*fc[10]-0.3535533905932737*fc[14])*dv; 
  Ghat_l[3] = (0.7071067811865475*fc[4]-1.224744871391589*fc[9])*wv+(0.2041241452319315*fc[0]-0.3535533905932737*fc[2])*dv; 
  Ghat_l[4] = (0.7071067811865475*fc[6]-1.224744871391589*fc[11])*wv+(0.2041241452319315*fc[13]-0.3535533905932737*fc[15])*dv; 
  Ghat_l[5] = (0.7071067811865475*fc[8]-1.224744871391589*fc[12])*wv+(0.2041241452319315*fc[1]-0.3535533905932737*fc[5])*dv; 
  Ghat_l[6] = (0.7071067811865475*fc[10]-1.224744871391589*fc[14])*wv+(0.2041241452319315*fc[3]-0.3535533905932737*fc[7])*dv; 
  Ghat_l[7] = (0.7071067811865475*fc[13]-1.224744871391589*fc[15])*wv+(0.2041241452319315*fc[6]-0.3535533905932737*fc[11])*dv; 

  } 
  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dx11; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dx11; 
  out[2] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dx11; 
  out[3] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dx11; 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dx11; 
  out[5] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dx11; 
  out[6] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dx11; 
  out[7] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dx11; 
  out[8] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dx11; 
  out[9] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dx11; 
  out[10] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dx11; 
  out[11] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dx11; 
  out[12] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dx11; 
  out[13] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dx11; 
  out[14] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dx11; 
  out[15] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dx11; 
} 
