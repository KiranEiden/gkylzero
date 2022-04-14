#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_surfx_1x2v_ser_p1(const double *w, const double *dxv, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fl/fc/fr:  Input Distribution function in left/center/right cells.
  // out:       Incremented distribution function in center cell.
  const double dx10 = 2/dxv[0]; 
  const double dv = dxv[1], wv = w[1]; 
  double Ghat_r[4]; 
  double Ghat_l[4]; 
  if (wv>0) { 

  Ghat_r[0] = (1.224744871391589*fc[1]+0.7071067811865475*fc[0])*wv+(0.3535533905932737*fc[4]+0.2041241452319315*fc[2])*dv; 
  Ghat_r[1] = (1.224744871391589*fc[4]+0.7071067811865475*fc[2])*wv+(0.3535533905932737*fc[1]+0.2041241452319315*fc[0])*dv; 
  Ghat_r[2] = (1.224744871391589*fc[5]+0.7071067811865475*fc[3])*wv+(0.3535533905932737*fc[7]+0.2041241452319315*fc[6])*dv; 
  Ghat_r[3] = (1.224744871391589*fc[7]+0.7071067811865475*fc[6])*wv+(0.3535533905932737*fc[5]+0.2041241452319315*fc[3])*dv; 

  Ghat_l[0] = (1.224744871391589*fl[1]+0.7071067811865475*fl[0])*wv+(0.3535533905932737*fl[4]+0.2041241452319315*fl[2])*dv; 
  Ghat_l[1] = (1.224744871391589*fl[4]+0.7071067811865475*fl[2])*wv+(0.3535533905932737*fl[1]+0.2041241452319315*fl[0])*dv; 
  Ghat_l[2] = (1.224744871391589*fl[5]+0.7071067811865475*fl[3])*wv+(0.3535533905932737*fl[7]+0.2041241452319315*fl[6])*dv; 
  Ghat_l[3] = (1.224744871391589*fl[7]+0.7071067811865475*fl[6])*wv+(0.3535533905932737*fl[5]+0.2041241452319315*fl[3])*dv; 

  } else { 

  Ghat_r[0] = -0.08333333333333333*((14.69693845669907*fr[1]-8.485281374238571*fr[0])*wv+(4.242640687119286*fr[4]-2.449489742783178*fr[2])*dv); 
  Ghat_r[1] = -0.08333333333333333*((14.69693845669907*fr[4]-8.485281374238571*fr[2])*wv+(4.242640687119286*fr[1]-2.449489742783178*fr[0])*dv); 
  Ghat_r[2] = -0.08333333333333333*((14.69693845669907*fr[5]-8.485281374238571*fr[3])*wv+(4.242640687119286*fr[7]-2.449489742783178*fr[6])*dv); 
  Ghat_r[3] = -0.08333333333333333*((14.69693845669907*fr[7]-8.485281374238571*fr[6])*wv+(4.242640687119286*fr[5]-2.449489742783178*fr[3])*dv); 

  Ghat_l[0] = -0.08333333333333333*((14.69693845669907*fc[1]-8.485281374238571*fc[0])*wv+(4.242640687119286*fc[4]-2.449489742783178*fc[2])*dv); 
  Ghat_l[1] = -0.08333333333333333*((14.69693845669907*fc[4]-8.485281374238571*fc[2])*wv+(4.242640687119286*fc[1]-2.449489742783178*fc[0])*dv); 
  Ghat_l[2] = -0.08333333333333333*((14.69693845669907*fc[5]-8.485281374238571*fc[3])*wv+(4.242640687119286*fc[7]-2.449489742783178*fc[6])*dv); 
  Ghat_l[3] = -0.08333333333333333*((14.69693845669907*fc[7]-8.485281374238571*fc[6])*wv+(4.242640687119286*fc[5]-2.449489742783178*fc[3])*dv); 

  } 
  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dx10; 
  out[1] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dx10; 
  out[2] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dx10; 
  out[3] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dx10; 
  out[4] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dx10; 
  out[5] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dx10; 
  out[6] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dx10; 
  out[7] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dx10; 
} 
