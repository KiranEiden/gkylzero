#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_surfvx_2x3v_ser_p1(const double *w, const double *dxv, const double amax, const double *qmem, const double *fl, const double *fc, const double *fr, double* restrict out) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // amax:      amax in global lax flux.
  // qmem:      q/m*EM fields.
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  const double dv10 = 2/dxv[2]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double dv3 = dxv[4], wv3 = w[4]; 
  const double *E0 = &qmem[0]; 
  const double *B0 = &qmem[12]; 
  const double *B1 = &qmem[16]; 
  const double *B2 = &qmem[20]; 

  double Ghat_r[16]; 
  double favg_r[16]; 
  double Ghat_l[16]; 
  double favg_l[16]; 
  double alpha[16]; 

  favg_r[0] = (-1.224744871391589*fr[3])+1.224744871391589*fc[3]+0.7071067811865475*(fr[0]+fc[0]); 
  favg_r[1] = (-1.224744871391589*fr[7])+1.224744871391589*fc[7]+0.7071067811865475*(fr[1]+fc[1]); 
  favg_r[2] = (-1.224744871391589*fr[8])+1.224744871391589*fc[8]+0.7071067811865475*(fr[2]+fc[2]); 
  favg_r[3] = (-1.224744871391589*fr[11])+1.224744871391589*fc[11]+0.7071067811865475*(fr[4]+fc[4]); 
  favg_r[4] = (-1.224744871391589*fr[14])+1.224744871391589*fc[14]+0.7071067811865475*(fr[5]+fc[5]); 
  favg_r[5] = (-1.224744871391589*fr[16])+1.224744871391589*fc[16]+0.7071067811865475*(fr[6]+fc[6]); 
  favg_r[6] = (-1.224744871391589*fr[18])+1.224744871391589*fc[18]+0.7071067811865475*(fr[9]+fc[9]); 
  favg_r[7] = (-1.224744871391589*fr[19])+1.224744871391589*fc[19]+0.7071067811865475*(fr[10]+fc[10]); 
  favg_r[8] = (-1.224744871391589*fr[21])+1.224744871391589*fc[21]+0.7071067811865475*(fr[12]+fc[12]); 
  favg_r[9] = (-1.224744871391589*fr[22])+1.224744871391589*fc[22]+0.7071067811865475*(fr[13]+fc[13]); 
  favg_r[10] = (-1.224744871391589*fr[25])+1.224744871391589*fc[25]+0.7071067811865475*(fr[15]+fc[15]); 
  favg_r[11] = (-1.224744871391589*fr[26])+1.224744871391589*fc[26]+0.7071067811865475*(fr[17]+fc[17]); 
  favg_r[12] = (-1.224744871391589*fr[27])+1.224744871391589*fc[27]+0.7071067811865475*(fr[20]+fc[20]); 
  favg_r[13] = (-1.224744871391589*fr[29])+1.224744871391589*fc[29]+0.7071067811865475*(fr[23]+fc[23]); 
  favg_r[14] = (-1.224744871391589*fr[30])+1.224744871391589*fc[30]+0.7071067811865475*(fr[24]+fc[24]); 
  favg_r[15] = (-1.224744871391589*fr[31])+1.224744871391589*fc[31]+0.7071067811865475*(fr[28]+fc[28]); 

  favg_l[0] = 1.224744871391589*fl[3]-1.224744871391589*fc[3]+0.7071067811865475*(fl[0]+fc[0]); 
  favg_l[1] = 1.224744871391589*fl[7]-1.224744871391589*fc[7]+0.7071067811865475*(fl[1]+fc[1]); 
  favg_l[2] = 1.224744871391589*fl[8]-1.224744871391589*fc[8]+0.7071067811865475*(fl[2]+fc[2]); 
  favg_l[3] = 1.224744871391589*fl[11]-1.224744871391589*fc[11]+0.7071067811865475*(fl[4]+fc[4]); 
  favg_l[4] = 1.224744871391589*fl[14]-1.224744871391589*fc[14]+0.7071067811865475*(fl[5]+fc[5]); 
  favg_l[5] = 1.224744871391589*fl[16]-1.224744871391589*fc[16]+0.7071067811865475*(fl[6]+fc[6]); 
  favg_l[6] = 1.224744871391589*fl[18]-1.224744871391589*fc[18]+0.7071067811865475*(fl[9]+fc[9]); 
  favg_l[7] = 1.224744871391589*fl[19]-1.224744871391589*fc[19]+0.7071067811865475*(fl[10]+fc[10]); 
  favg_l[8] = 1.224744871391589*fl[21]-1.224744871391589*fc[21]+0.7071067811865475*(fl[12]+fc[12]); 
  favg_l[9] = 1.224744871391589*fl[22]-1.224744871391589*fc[22]+0.7071067811865475*(fl[13]+fc[13]); 
  favg_l[10] = 1.224744871391589*fl[25]-1.224744871391589*fc[25]+0.7071067811865475*(fl[15]+fc[15]); 
  favg_l[11] = 1.224744871391589*fl[26]-1.224744871391589*fc[26]+0.7071067811865475*(fl[17]+fc[17]); 
  favg_l[12] = 1.224744871391589*fl[27]-1.224744871391589*fc[27]+0.7071067811865475*(fl[20]+fc[20]); 
  favg_l[13] = 1.224744871391589*fl[29]-1.224744871391589*fc[29]+0.7071067811865475*(fl[23]+fc[23]); 
  favg_l[14] = 1.224744871391589*fl[30]-1.224744871391589*fc[30]+0.7071067811865475*(fl[24]+fc[24]); 
  favg_l[15] = 1.224744871391589*fl[31]-1.224744871391589*fc[31]+0.7071067811865475*(fl[28]+fc[28]); 

  alpha[0] = 2.0*(B2[0]*wv2+E0[0])-2.0*B1[0]*wv3; 
  alpha[1] = 2.0*(B2[1]*wv2+E0[1])-2.0*B1[1]*wv3; 
  alpha[2] = 2.0*(B2[2]*wv2+E0[2])-2.0*B1[2]*wv3; 
  alpha[3] = 0.5773502691896258*B2[0]*dv2; 
  alpha[4] = -0.5773502691896258*B1[0]*dv3; 
  alpha[5] = 2.0*(B2[3]*wv2+E0[3])-2.0*B1[3]*wv3; 
  alpha[6] = 0.5773502691896258*B2[1]*dv2; 
  alpha[7] = 0.5773502691896258*B2[2]*dv2; 
  alpha[8] = -0.5773502691896258*B1[1]*dv3; 
  alpha[9] = -0.5773502691896258*B1[2]*dv3; 
  alpha[11] = 0.5773502691896258*B2[3]*dv2; 
  alpha[12] = -0.5773502691896258*B1[3]*dv3; 

  double amid = 0.25*alpha[0]; 

  Ghat_r[0] = (0.6123724356957944*(fr[3]+fc[3])-0.3535533905932737*fr[0]+0.3535533905932737*fc[0])*amax+0.125*(alpha[12]*favg_r[12]+alpha[11]*favg_r[11]+alpha[9]*favg_r[9]+alpha[8]*favg_r[8]+alpha[7]*favg_r[7]+alpha[6]*favg_r[6]+alpha[5]*favg_r[5]+alpha[4]*favg_r[4]+alpha[3]*favg_r[3]+alpha[2]*favg_r[2]+alpha[1]*favg_r[1]+alpha[0]*favg_r[0]); 
  Ghat_r[1] = (0.6123724356957944*(fr[7]+fc[7])-0.3535533905932737*fr[1]+0.3535533905932737*fc[1])*amax+0.125*(alpha[9]*favg_r[12]+favg_r[9]*alpha[12]+alpha[7]*favg_r[11]+favg_r[7]*alpha[11]+alpha[4]*favg_r[8]+favg_r[4]*alpha[8]+alpha[3]*favg_r[6]+favg_r[3]*alpha[6]+alpha[2]*favg_r[5]+favg_r[2]*alpha[5]+alpha[0]*favg_r[1]+favg_r[0]*alpha[1]); 
  Ghat_r[2] = (0.6123724356957944*(fr[8]+fc[8])-0.3535533905932737*fr[2]+0.3535533905932737*fc[2])*amax+0.125*(alpha[8]*favg_r[12]+favg_r[8]*alpha[12]+alpha[6]*favg_r[11]+favg_r[6]*alpha[11]+alpha[4]*favg_r[9]+favg_r[4]*alpha[9]+alpha[3]*favg_r[7]+favg_r[3]*alpha[7]+alpha[1]*favg_r[5]+favg_r[1]*alpha[5]+alpha[0]*favg_r[2]+favg_r[0]*alpha[2]); 
  Ghat_r[3] = (0.6123724356957944*(fr[11]+fc[11])-0.3535533905932737*fr[4]+0.3535533905932737*fc[4])*amax+0.125*(alpha[12]*favg_r[15]+alpha[9]*favg_r[14]+alpha[8]*favg_r[13]+alpha[5]*favg_r[11]+favg_r[5]*alpha[11]+alpha[4]*favg_r[10]+alpha[2]*favg_r[7]+favg_r[2]*alpha[7]+alpha[1]*favg_r[6]+favg_r[1]*alpha[6]+alpha[0]*favg_r[3]+favg_r[0]*alpha[3]); 
  Ghat_r[4] = (0.6123724356957944*(fr[14]+fc[14])-0.3535533905932737*fr[5]+0.3535533905932737*fc[5])*amax+0.125*(alpha[11]*favg_r[15]+alpha[7]*favg_r[14]+alpha[6]*favg_r[13]+alpha[5]*favg_r[12]+favg_r[5]*alpha[12]+alpha[3]*favg_r[10]+alpha[2]*favg_r[9]+favg_r[2]*alpha[9]+alpha[1]*favg_r[8]+favg_r[1]*alpha[8]+alpha[0]*favg_r[4]+favg_r[0]*alpha[4]); 
  Ghat_r[5] = (0.6123724356957944*(fr[16]+fc[16])-0.3535533905932737*fr[6]+0.3535533905932737*fc[6])*amax+0.125*(alpha[4]*favg_r[12]+favg_r[4]*alpha[12]+alpha[3]*favg_r[11]+favg_r[3]*alpha[11]+alpha[8]*favg_r[9]+favg_r[8]*alpha[9]+alpha[6]*favg_r[7]+favg_r[6]*alpha[7]+alpha[0]*favg_r[5]+favg_r[0]*alpha[5]+alpha[1]*favg_r[2]+favg_r[1]*alpha[2]); 
  Ghat_r[6] = (0.6123724356957944*(fr[18]+fc[18])-0.3535533905932737*fr[9]+0.3535533905932737*fc[9])*amax+0.125*(alpha[9]*favg_r[15]+alpha[12]*favg_r[14]+alpha[4]*favg_r[13]+alpha[2]*favg_r[11]+favg_r[2]*alpha[11]+alpha[8]*favg_r[10]+alpha[5]*favg_r[7]+favg_r[5]*alpha[7]+alpha[0]*favg_r[6]+favg_r[0]*alpha[6]+alpha[1]*favg_r[3]+favg_r[1]*alpha[3]); 
  Ghat_r[7] = (0.6123724356957944*(fr[19]+fc[19])-0.3535533905932737*fr[10]+0.3535533905932737*fc[10])*amax+0.125*(alpha[8]*favg_r[15]+alpha[4]*favg_r[14]+alpha[12]*favg_r[13]+alpha[1]*favg_r[11]+favg_r[1]*alpha[11]+alpha[9]*favg_r[10]+alpha[0]*favg_r[7]+favg_r[0]*alpha[7]+alpha[5]*favg_r[6]+favg_r[5]*alpha[6]+alpha[2]*favg_r[3]+favg_r[2]*alpha[3]); 
  Ghat_r[8] = (0.6123724356957944*(fr[21]+fc[21])-0.3535533905932737*fr[12]+0.3535533905932737*fc[12])*amax+0.125*(alpha[7]*favg_r[15]+alpha[11]*favg_r[14]+alpha[3]*favg_r[13]+alpha[2]*favg_r[12]+favg_r[2]*alpha[12]+alpha[6]*favg_r[10]+alpha[5]*favg_r[9]+favg_r[5]*alpha[9]+alpha[0]*favg_r[8]+favg_r[0]*alpha[8]+alpha[1]*favg_r[4]+favg_r[1]*alpha[4]); 
  Ghat_r[9] = (0.6123724356957944*(fr[22]+fc[22])-0.3535533905932737*fr[13]+0.3535533905932737*fc[13])*amax+0.125*(alpha[6]*favg_r[15]+alpha[3]*favg_r[14]+alpha[11]*favg_r[13]+alpha[1]*favg_r[12]+favg_r[1]*alpha[12]+alpha[7]*favg_r[10]+alpha[0]*favg_r[9]+favg_r[0]*alpha[9]+alpha[5]*favg_r[8]+favg_r[5]*alpha[8]+alpha[2]*favg_r[4]+favg_r[2]*alpha[4]); 
  Ghat_r[10] = (0.6123724356957944*(fr[25]+fc[25])-0.3535533905932737*fr[15]+0.3535533905932737*fc[15])*amax+0.125*(alpha[5]*favg_r[15]+alpha[2]*favg_r[14]+alpha[1]*favg_r[13]+alpha[11]*favg_r[12]+favg_r[11]*alpha[12]+alpha[0]*favg_r[10]+alpha[7]*favg_r[9]+favg_r[7]*alpha[9]+alpha[6]*favg_r[8]+favg_r[6]*alpha[8]+alpha[3]*favg_r[4]+favg_r[3]*alpha[4]); 
  Ghat_r[11] = (0.6123724356957944*(fr[26]+fc[26])-0.3535533905932737*fr[17]+0.3535533905932737*fc[17])*amax+0.125*(alpha[4]*favg_r[15]+alpha[8]*favg_r[14]+alpha[9]*favg_r[13]+favg_r[10]*alpha[12]+alpha[0]*favg_r[11]+favg_r[0]*alpha[11]+alpha[1]*favg_r[7]+favg_r[1]*alpha[7]+alpha[2]*favg_r[6]+favg_r[2]*alpha[6]+alpha[3]*favg_r[5]+favg_r[3]*alpha[5]); 
  Ghat_r[12] = (0.6123724356957944*(fr[27]+fc[27])-0.3535533905932737*fr[20]+0.3535533905932737*fc[20])*amax+0.125*(alpha[3]*favg_r[15]+alpha[6]*favg_r[14]+alpha[7]*favg_r[13]+alpha[0]*favg_r[12]+favg_r[0]*alpha[12]+favg_r[10]*alpha[11]+alpha[1]*favg_r[9]+favg_r[1]*alpha[9]+alpha[2]*favg_r[8]+favg_r[2]*alpha[8]+alpha[4]*favg_r[5]+favg_r[4]*alpha[5]); 
  Ghat_r[13] = (0.6123724356957944*(fr[29]+fc[29])-0.3535533905932737*fr[23]+0.3535533905932737*fc[23])*amax+0.125*(alpha[2]*favg_r[15]+alpha[5]*favg_r[14]+alpha[0]*favg_r[13]+alpha[7]*favg_r[12]+favg_r[7]*alpha[12]+alpha[9]*favg_r[11]+favg_r[9]*alpha[11]+alpha[1]*favg_r[10]+alpha[3]*favg_r[8]+favg_r[3]*alpha[8]+alpha[4]*favg_r[6]+favg_r[4]*alpha[6]); 
  Ghat_r[14] = (0.6123724356957944*(fr[30]+fc[30])-0.3535533905932737*fr[24]+0.3535533905932737*fc[24])*amax+0.125*(alpha[1]*favg_r[15]+alpha[0]*favg_r[14]+alpha[5]*favg_r[13]+alpha[6]*favg_r[12]+favg_r[6]*alpha[12]+alpha[8]*favg_r[11]+favg_r[8]*alpha[11]+alpha[2]*favg_r[10]+alpha[3]*favg_r[9]+favg_r[3]*alpha[9]+alpha[4]*favg_r[7]+favg_r[4]*alpha[7]); 
  Ghat_r[15] = (0.6123724356957944*(fr[31]+fc[31])-0.3535533905932737*fr[28]+0.3535533905932737*fc[28])*amax+0.125*(alpha[0]*favg_r[15]+alpha[1]*favg_r[14]+alpha[2]*favg_r[13]+alpha[3]*favg_r[12]+favg_r[3]*alpha[12]+alpha[4]*favg_r[11]+favg_r[4]*alpha[11]+alpha[5]*favg_r[10]+alpha[6]*favg_r[9]+favg_r[6]*alpha[9]+alpha[7]*favg_r[8]+favg_r[7]*alpha[8]); 

  Ghat_l[0] = (0.6123724356957944*(fl[3]+fc[3])+0.3535533905932737*fl[0]-0.3535533905932737*fc[0])*amax+0.125*(alpha[12]*favg_l[12]+alpha[11]*favg_l[11]+alpha[9]*favg_l[9]+alpha[8]*favg_l[8]+alpha[7]*favg_l[7]+alpha[6]*favg_l[6]+alpha[5]*favg_l[5]+alpha[4]*favg_l[4]+alpha[3]*favg_l[3]+alpha[2]*favg_l[2]+alpha[1]*favg_l[1]+alpha[0]*favg_l[0]); 
  Ghat_l[1] = (0.6123724356957944*(fl[7]+fc[7])+0.3535533905932737*fl[1]-0.3535533905932737*fc[1])*amax+0.125*(alpha[9]*favg_l[12]+favg_l[9]*alpha[12]+alpha[7]*favg_l[11]+favg_l[7]*alpha[11]+alpha[4]*favg_l[8]+favg_l[4]*alpha[8]+alpha[3]*favg_l[6]+favg_l[3]*alpha[6]+alpha[2]*favg_l[5]+favg_l[2]*alpha[5]+alpha[0]*favg_l[1]+favg_l[0]*alpha[1]); 
  Ghat_l[2] = (0.6123724356957944*(fl[8]+fc[8])+0.3535533905932737*fl[2]-0.3535533905932737*fc[2])*amax+0.125*(alpha[8]*favg_l[12]+favg_l[8]*alpha[12]+alpha[6]*favg_l[11]+favg_l[6]*alpha[11]+alpha[4]*favg_l[9]+favg_l[4]*alpha[9]+alpha[3]*favg_l[7]+favg_l[3]*alpha[7]+alpha[1]*favg_l[5]+favg_l[1]*alpha[5]+alpha[0]*favg_l[2]+favg_l[0]*alpha[2]); 
  Ghat_l[3] = (0.6123724356957944*(fl[11]+fc[11])+0.3535533905932737*fl[4]-0.3535533905932737*fc[4])*amax+0.125*(alpha[12]*favg_l[15]+alpha[9]*favg_l[14]+alpha[8]*favg_l[13]+alpha[5]*favg_l[11]+favg_l[5]*alpha[11]+alpha[4]*favg_l[10]+alpha[2]*favg_l[7]+favg_l[2]*alpha[7]+alpha[1]*favg_l[6]+favg_l[1]*alpha[6]+alpha[0]*favg_l[3]+favg_l[0]*alpha[3]); 
  Ghat_l[4] = (0.6123724356957944*(fl[14]+fc[14])+0.3535533905932737*fl[5]-0.3535533905932737*fc[5])*amax+0.125*(alpha[11]*favg_l[15]+alpha[7]*favg_l[14]+alpha[6]*favg_l[13]+alpha[5]*favg_l[12]+favg_l[5]*alpha[12]+alpha[3]*favg_l[10]+alpha[2]*favg_l[9]+favg_l[2]*alpha[9]+alpha[1]*favg_l[8]+favg_l[1]*alpha[8]+alpha[0]*favg_l[4]+favg_l[0]*alpha[4]); 
  Ghat_l[5] = (0.6123724356957944*(fl[16]+fc[16])+0.3535533905932737*fl[6]-0.3535533905932737*fc[6])*amax+0.125*(alpha[4]*favg_l[12]+favg_l[4]*alpha[12]+alpha[3]*favg_l[11]+favg_l[3]*alpha[11]+alpha[8]*favg_l[9]+favg_l[8]*alpha[9]+alpha[6]*favg_l[7]+favg_l[6]*alpha[7]+alpha[0]*favg_l[5]+favg_l[0]*alpha[5]+alpha[1]*favg_l[2]+favg_l[1]*alpha[2]); 
  Ghat_l[6] = (0.6123724356957944*(fl[18]+fc[18])+0.3535533905932737*fl[9]-0.3535533905932737*fc[9])*amax+0.125*(alpha[9]*favg_l[15]+alpha[12]*favg_l[14]+alpha[4]*favg_l[13]+alpha[2]*favg_l[11]+favg_l[2]*alpha[11]+alpha[8]*favg_l[10]+alpha[5]*favg_l[7]+favg_l[5]*alpha[7]+alpha[0]*favg_l[6]+favg_l[0]*alpha[6]+alpha[1]*favg_l[3]+favg_l[1]*alpha[3]); 
  Ghat_l[7] = (0.6123724356957944*(fl[19]+fc[19])+0.3535533905932737*fl[10]-0.3535533905932737*fc[10])*amax+0.125*(alpha[8]*favg_l[15]+alpha[4]*favg_l[14]+alpha[12]*favg_l[13]+alpha[1]*favg_l[11]+favg_l[1]*alpha[11]+alpha[9]*favg_l[10]+alpha[0]*favg_l[7]+favg_l[0]*alpha[7]+alpha[5]*favg_l[6]+favg_l[5]*alpha[6]+alpha[2]*favg_l[3]+favg_l[2]*alpha[3]); 
  Ghat_l[8] = (0.6123724356957944*(fl[21]+fc[21])+0.3535533905932737*fl[12]-0.3535533905932737*fc[12])*amax+0.125*(alpha[7]*favg_l[15]+alpha[11]*favg_l[14]+alpha[3]*favg_l[13]+alpha[2]*favg_l[12]+favg_l[2]*alpha[12]+alpha[6]*favg_l[10]+alpha[5]*favg_l[9]+favg_l[5]*alpha[9]+alpha[0]*favg_l[8]+favg_l[0]*alpha[8]+alpha[1]*favg_l[4]+favg_l[1]*alpha[4]); 
  Ghat_l[9] = (0.6123724356957944*(fl[22]+fc[22])+0.3535533905932737*fl[13]-0.3535533905932737*fc[13])*amax+0.125*(alpha[6]*favg_l[15]+alpha[3]*favg_l[14]+alpha[11]*favg_l[13]+alpha[1]*favg_l[12]+favg_l[1]*alpha[12]+alpha[7]*favg_l[10]+alpha[0]*favg_l[9]+favg_l[0]*alpha[9]+alpha[5]*favg_l[8]+favg_l[5]*alpha[8]+alpha[2]*favg_l[4]+favg_l[2]*alpha[4]); 
  Ghat_l[10] = (0.6123724356957944*(fl[25]+fc[25])+0.3535533905932737*fl[15]-0.3535533905932737*fc[15])*amax+0.125*(alpha[5]*favg_l[15]+alpha[2]*favg_l[14]+alpha[1]*favg_l[13]+alpha[11]*favg_l[12]+favg_l[11]*alpha[12]+alpha[0]*favg_l[10]+alpha[7]*favg_l[9]+favg_l[7]*alpha[9]+alpha[6]*favg_l[8]+favg_l[6]*alpha[8]+alpha[3]*favg_l[4]+favg_l[3]*alpha[4]); 
  Ghat_l[11] = (0.6123724356957944*(fl[26]+fc[26])+0.3535533905932737*fl[17]-0.3535533905932737*fc[17])*amax+0.125*(alpha[4]*favg_l[15]+alpha[8]*favg_l[14]+alpha[9]*favg_l[13]+favg_l[10]*alpha[12]+alpha[0]*favg_l[11]+favg_l[0]*alpha[11]+alpha[1]*favg_l[7]+favg_l[1]*alpha[7]+alpha[2]*favg_l[6]+favg_l[2]*alpha[6]+alpha[3]*favg_l[5]+favg_l[3]*alpha[5]); 
  Ghat_l[12] = (0.6123724356957944*(fl[27]+fc[27])+0.3535533905932737*fl[20]-0.3535533905932737*fc[20])*amax+0.125*(alpha[3]*favg_l[15]+alpha[6]*favg_l[14]+alpha[7]*favg_l[13]+alpha[0]*favg_l[12]+favg_l[0]*alpha[12]+favg_l[10]*alpha[11]+alpha[1]*favg_l[9]+favg_l[1]*alpha[9]+alpha[2]*favg_l[8]+favg_l[2]*alpha[8]+alpha[4]*favg_l[5]+favg_l[4]*alpha[5]); 
  Ghat_l[13] = (0.6123724356957944*(fl[29]+fc[29])+0.3535533905932737*fl[23]-0.3535533905932737*fc[23])*amax+0.125*(alpha[2]*favg_l[15]+alpha[5]*favg_l[14]+alpha[0]*favg_l[13]+alpha[7]*favg_l[12]+favg_l[7]*alpha[12]+alpha[9]*favg_l[11]+favg_l[9]*alpha[11]+alpha[1]*favg_l[10]+alpha[3]*favg_l[8]+favg_l[3]*alpha[8]+alpha[4]*favg_l[6]+favg_l[4]*alpha[6]); 
  Ghat_l[14] = (0.6123724356957944*(fl[30]+fc[30])+0.3535533905932737*fl[24]-0.3535533905932737*fc[24])*amax+0.125*(alpha[1]*favg_l[15]+alpha[0]*favg_l[14]+alpha[5]*favg_l[13]+alpha[6]*favg_l[12]+favg_l[6]*alpha[12]+alpha[8]*favg_l[11]+favg_l[8]*alpha[11]+alpha[2]*favg_l[10]+alpha[3]*favg_l[9]+favg_l[3]*alpha[9]+alpha[4]*favg_l[7]+favg_l[4]*alpha[7]); 
  Ghat_l[15] = (0.6123724356957944*(fl[31]+fc[31])+0.3535533905932737*fl[28]-0.3535533905932737*fc[28])*amax+0.125*(alpha[0]*favg_l[15]+alpha[1]*favg_l[14]+alpha[2]*favg_l[13]+alpha[3]*favg_l[12]+favg_l[3]*alpha[12]+alpha[4]*favg_l[11]+favg_l[4]*alpha[11]+alpha[5]*favg_l[10]+alpha[6]*favg_l[9]+favg_l[6]*alpha[9]+alpha[7]*favg_l[8]+favg_l[7]*alpha[8]); 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv10; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv10; 
  out[2] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv10; 
  out[3] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv10; 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv10; 
  out[5] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dv10; 
  out[6] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dv10; 
  out[7] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv10; 
  out[8] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv10; 
  out[9] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dv10; 
  out[10] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dv10; 
  out[11] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv10; 
  out[12] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*dv10; 
  out[13] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*dv10; 
  out[14] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dv10; 
  out[15] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*dv10; 
  out[16] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dv10; 
  out[17] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*dv10; 
  out[18] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dv10; 
  out[19] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dv10; 
  out[20] += (0.7071067811865475*Ghat_l[12]-0.7071067811865475*Ghat_r[12])*dv10; 
  out[21] += -1.224744871391589*(Ghat_r[8]+Ghat_l[8])*dv10; 
  out[22] += -1.224744871391589*(Ghat_r[9]+Ghat_l[9])*dv10; 
  out[23] += (0.7071067811865475*Ghat_l[13]-0.7071067811865475*Ghat_r[13])*dv10; 
  out[24] += (0.7071067811865475*Ghat_l[14]-0.7071067811865475*Ghat_r[14])*dv10; 
  out[25] += -1.224744871391589*(Ghat_r[10]+Ghat_l[10])*dv10; 
  out[26] += -1.224744871391589*(Ghat_r[11]+Ghat_l[11])*dv10; 
  out[27] += -1.224744871391589*(Ghat_r[12]+Ghat_l[12])*dv10; 
  out[28] += (0.7071067811865475*Ghat_l[15]-0.7071067811865475*Ghat_r[15])*dv10; 
  out[29] += -1.224744871391589*(Ghat_r[13]+Ghat_l[13])*dv10; 
  out[30] += -1.224744871391589*(Ghat_r[14]+Ghat_l[14])*dv10; 
  out[31] += -1.224744871391589*(Ghat_r[15]+Ghat_l[15])*dv10; 

  return fabs(amid); 
} 
