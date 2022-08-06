#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_surfx_2x3v_ser_p1(const double *w, const double *dxv, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fl/fc/fr:  Input Distribution function in left/center/right cells.
  // out:       Incremented distribution function in center cell.
  const double dx10 = 2/dxv[0]; 
  const double dv = dxv[2], wv = w[2]; 
  double Ghat_r[40]; 
  double Ghat_l[40]; 
  if (wv>0) { 

  Ghat_r[0] = (1.224744871391589*fc[1]+0.7071067811865475*fc[0])*wv+(0.3535533905932737*fc[7]+0.2041241452319315*fc[3])*dv; 
  Ghat_r[1] = (1.224744871391589*fc[6]+0.7071067811865475*fc[2])*wv+(0.3535533905932737*fc[16]+0.2041241452319315*fc[8])*dv; 
  Ghat_r[2] = (1.224744871391589*fc[7]+0.7071067811865475*fc[3])*wv+(0.3162277660168379*fc[33]+0.1825741858350554*fc[32]+0.3535533905932737*fc[1]+0.2041241452319315*fc[0])*dv; 
  Ghat_r[3] = (1.224744871391589*fc[9]+0.7071067811865475*fc[4])*wv+(0.3535533905932737*fc[18]+0.2041241452319315*fc[11])*dv; 
  Ghat_r[4] = (1.224744871391589*fc[12]+0.7071067811865475*fc[5])*wv+(0.3535533905932737*fc[21]+0.2041241452319315*fc[14])*dv; 
  Ghat_r[5] = (1.224744871391589*fc[16]+0.7071067811865475*fc[8])*wv+(0.3162277660168379*fc[37]+0.1825741858350554*fc[34]+0.3535533905932737*fc[6]+0.2041241452319315*fc[2])*dv; 
  Ghat_r[6] = (1.224744871391589*fc[17]+0.7071067811865475*fc[10])*wv+(0.3535533905932737*fc[26]+0.2041241452319315*fc[19])*dv; 
  Ghat_r[7] = (1.224744871391589*fc[18]+0.7071067811865475*fc[11])*wv+(0.3162277660168379*fc[38]+0.1825741858350554*fc[35]+0.3535533905932737*fc[9]+0.2041241452319315*fc[4])*dv; 
  Ghat_r[8] = (1.224744871391589*fc[20]+0.7071067811865475*fc[13])*wv+(0.3535533905932737*fc[27]+0.2041241452319315*fc[22])*dv; 
  Ghat_r[9] = (1.224744871391589*fc[21]+0.7071067811865475*fc[14])*wv+(0.3162277660168379*fc[40]+0.1825741858350554*fc[36]+0.3535533905932737*fc[12]+0.2041241452319315*fc[5])*dv; 
  Ghat_r[10] = (1.224744871391589*fc[23]+0.7071067811865475*fc[15])*wv+(0.3535533905932737*fc[29]+0.2041241452319315*fc[25])*dv; 
  Ghat_r[11] = (1.224744871391589*fc[26]+0.7071067811865475*fc[19])*wv+(0.3162277660168379*fc[43]+0.1825741858350554*fc[39]+0.3535533905932737*fc[17]+0.2041241452319315*fc[10])*dv; 
  Ghat_r[12] = (1.224744871391589*fc[27]+0.7071067811865475*fc[22])*wv+(0.3162277660168379*fc[44]+0.1825741858350554*fc[41]+0.3535533905932737*fc[20]+0.2041241452319315*fc[13])*dv; 
  Ghat_r[13] = (1.224744871391589*fc[28]+0.7071067811865475*fc[24])*wv+(0.3535533905932737*fc[31]+0.2041241452319315*fc[30])*dv; 
  Ghat_r[14] = (1.224744871391589*fc[29]+0.7071067811865475*fc[25])*wv+(0.3162277660168379*fc[45]+0.1825741858350554*fc[42]+0.3535533905932737*fc[23]+0.2041241452319315*fc[15])*dv; 
  Ghat_r[15] = (1.224744871391589*fc[31]+0.7071067811865475*fc[30])*wv+(0.3162277660168379*fc[47]+0.1825741858350554*fc[46]+0.3535533905932737*fc[28]+0.2041241452319315*fc[24])*dv; 
  Ghat_r[16] = (1.224744871391589*fc[33]+0.7071067811865475*fc[32])*wv+(0.3162277660168379*fc[7]+0.1825741858350554*fc[3])*dv; 
  Ghat_r[17] = (1.224744871391589*fc[37]+0.7071067811865475*fc[34])*wv+(0.3162277660168379*fc[16]+0.1825741858350554*fc[8])*dv; 
  Ghat_r[18] = (1.224744871391589*fc[38]+0.7071067811865475*fc[35])*wv+(0.3162277660168379*fc[18]+0.1825741858350554*fc[11])*dv; 
  Ghat_r[19] = (1.224744871391589*fc[40]+0.7071067811865475*fc[36])*wv+(0.3162277660168379*fc[21]+0.1825741858350554*fc[14])*dv; 
  Ghat_r[20] = (1.224744871391589*fc[43]+0.7071067811865475*fc[39])*wv+(0.3162277660168379*fc[26]+0.1825741858350554*fc[19])*dv; 
  Ghat_r[21] = (1.224744871391589*fc[44]+0.7071067811865475*fc[41])*wv+(0.3162277660168379*fc[27]+0.1825741858350554*fc[22])*dv; 
  Ghat_r[22] = (1.224744871391589*fc[45]+0.7071067811865475*fc[42])*wv+(0.3162277660168379*fc[29]+0.1825741858350554*fc[25])*dv; 
  Ghat_r[23] = (1.224744871391589*fc[47]+0.7071067811865475*fc[46])*wv+(0.3162277660168379*fc[31]+0.1825741858350554*fc[30])*dv; 
  Ghat_r[24] = (1.224744871391589*fc[49]+0.7071067811865475*fc[48])*wv+(0.3535533905932737*fc[54]+0.2041241452319315*fc[51])*dv; 
  Ghat_r[25] = (1.224744871391589*fc[53]+0.7071067811865475*fc[50])*wv+(0.3535533905932737*fc[59]+0.2041241452319315*fc[55])*dv; 
  Ghat_r[26] = (1.224744871391589*fc[54]+0.7071067811865475*fc[51])*wv+(0.3535533905932737*fc[49]+0.2041241452319315*fc[48])*dv; 
  Ghat_r[27] = (1.224744871391589*fc[56]+0.7071067811865475*fc[52])*wv+(0.3535533905932737*fc[61]+0.2041241452319315*fc[58])*dv; 
  Ghat_r[28] = (1.224744871391589*fc[59]+0.7071067811865475*fc[55])*wv+(0.3535533905932737*fc[53]+0.2041241452319315*fc[50])*dv; 
  Ghat_r[29] = (1.224744871391589*fc[60]+0.7071067811865475*fc[57])*wv+(0.3535533905932737*fc[63]+0.2041241452319315*fc[62])*dv; 
  Ghat_r[30] = (1.224744871391589*fc[61]+0.7071067811865475*fc[58])*wv+(0.3535533905932737*fc[56]+0.2041241452319315*fc[52])*dv; 
  Ghat_r[31] = (1.224744871391589*fc[63]+0.7071067811865475*fc[62])*wv+(0.3535533905932737*fc[60]+0.2041241452319315*fc[57])*dv; 
  Ghat_r[32] = (1.224744871391589*fc[65]+0.7071067811865475*fc[64])*wv+(0.3535533905932737*fc[70]+0.2041241452319315*fc[67])*dv; 
  Ghat_r[33] = (1.224744871391589*fc[69]+0.7071067811865475*fc[66])*wv+(0.3535533905932737*fc[75]+0.2041241452319315*fc[71])*dv; 
  Ghat_r[34] = (1.224744871391589*fc[70]+0.7071067811865475*fc[67])*wv+(0.3535533905932737*fc[65]+0.2041241452319315*fc[64])*dv; 
  Ghat_r[35] = (1.224744871391589*fc[72]+0.7071067811865475*fc[68])*wv+(0.3535533905932737*fc[77]+0.2041241452319315*fc[74])*dv; 
  Ghat_r[36] = (1.224744871391589*fc[75]+0.7071067811865475*fc[71])*wv+(0.3535533905932737*fc[69]+0.2041241452319315*fc[66])*dv; 
  Ghat_r[37] = (1.224744871391589*fc[76]+0.7071067811865475*fc[73])*wv+(0.3535533905932737*fc[79]+0.2041241452319315*fc[78])*dv; 
  Ghat_r[38] = (1.224744871391589*fc[77]+0.7071067811865475*fc[74])*wv+(0.3535533905932737*fc[72]+0.2041241452319315*fc[68])*dv; 
  Ghat_r[39] = (1.224744871391589*fc[79]+0.7071067811865475*fc[78])*wv+(0.3535533905932737*fc[76]+0.2041241452319315*fc[73])*dv; 

  Ghat_l[0] = (1.224744871391589*fl[1]+0.7071067811865475*fl[0])*wv+(0.3535533905932737*fl[7]+0.2041241452319315*fl[3])*dv; 
  Ghat_l[1] = (1.224744871391589*fl[6]+0.7071067811865475*fl[2])*wv+(0.3535533905932737*fl[16]+0.2041241452319315*fl[8])*dv; 
  Ghat_l[2] = (1.224744871391589*fl[7]+0.7071067811865475*fl[3])*wv+(0.3162277660168379*fl[33]+0.1825741858350554*fl[32]+0.3535533905932737*fl[1]+0.2041241452319315*fl[0])*dv; 
  Ghat_l[3] = (1.224744871391589*fl[9]+0.7071067811865475*fl[4])*wv+(0.3535533905932737*fl[18]+0.2041241452319315*fl[11])*dv; 
  Ghat_l[4] = (1.224744871391589*fl[12]+0.7071067811865475*fl[5])*wv+(0.3535533905932737*fl[21]+0.2041241452319315*fl[14])*dv; 
  Ghat_l[5] = (1.224744871391589*fl[16]+0.7071067811865475*fl[8])*wv+(0.3162277660168379*fl[37]+0.1825741858350554*fl[34]+0.3535533905932737*fl[6]+0.2041241452319315*fl[2])*dv; 
  Ghat_l[6] = (1.224744871391589*fl[17]+0.7071067811865475*fl[10])*wv+(0.3535533905932737*fl[26]+0.2041241452319315*fl[19])*dv; 
  Ghat_l[7] = (1.224744871391589*fl[18]+0.7071067811865475*fl[11])*wv+(0.3162277660168379*fl[38]+0.1825741858350554*fl[35]+0.3535533905932737*fl[9]+0.2041241452319315*fl[4])*dv; 
  Ghat_l[8] = (1.224744871391589*fl[20]+0.7071067811865475*fl[13])*wv+(0.3535533905932737*fl[27]+0.2041241452319315*fl[22])*dv; 
  Ghat_l[9] = (1.224744871391589*fl[21]+0.7071067811865475*fl[14])*wv+(0.3162277660168379*fl[40]+0.1825741858350554*fl[36]+0.3535533905932737*fl[12]+0.2041241452319315*fl[5])*dv; 
  Ghat_l[10] = (1.224744871391589*fl[23]+0.7071067811865475*fl[15])*wv+(0.3535533905932737*fl[29]+0.2041241452319315*fl[25])*dv; 
  Ghat_l[11] = (1.224744871391589*fl[26]+0.7071067811865475*fl[19])*wv+(0.3162277660168379*fl[43]+0.1825741858350554*fl[39]+0.3535533905932737*fl[17]+0.2041241452319315*fl[10])*dv; 
  Ghat_l[12] = (1.224744871391589*fl[27]+0.7071067811865475*fl[22])*wv+(0.3162277660168379*fl[44]+0.1825741858350554*fl[41]+0.3535533905932737*fl[20]+0.2041241452319315*fl[13])*dv; 
  Ghat_l[13] = (1.224744871391589*fl[28]+0.7071067811865475*fl[24])*wv+(0.3535533905932737*fl[31]+0.2041241452319315*fl[30])*dv; 
  Ghat_l[14] = (1.224744871391589*fl[29]+0.7071067811865475*fl[25])*wv+(0.3162277660168379*fl[45]+0.1825741858350554*fl[42]+0.3535533905932737*fl[23]+0.2041241452319315*fl[15])*dv; 
  Ghat_l[15] = (1.224744871391589*fl[31]+0.7071067811865475*fl[30])*wv+(0.3162277660168379*fl[47]+0.1825741858350554*fl[46]+0.3535533905932737*fl[28]+0.2041241452319315*fl[24])*dv; 
  Ghat_l[16] = (1.224744871391589*fl[33]+0.7071067811865475*fl[32])*wv+(0.3162277660168379*fl[7]+0.1825741858350554*fl[3])*dv; 
  Ghat_l[17] = (1.224744871391589*fl[37]+0.7071067811865475*fl[34])*wv+(0.3162277660168379*fl[16]+0.1825741858350554*fl[8])*dv; 
  Ghat_l[18] = (1.224744871391589*fl[38]+0.7071067811865475*fl[35])*wv+(0.3162277660168379*fl[18]+0.1825741858350554*fl[11])*dv; 
  Ghat_l[19] = (1.224744871391589*fl[40]+0.7071067811865475*fl[36])*wv+(0.3162277660168379*fl[21]+0.1825741858350554*fl[14])*dv; 
  Ghat_l[20] = (1.224744871391589*fl[43]+0.7071067811865475*fl[39])*wv+(0.3162277660168379*fl[26]+0.1825741858350554*fl[19])*dv; 
  Ghat_l[21] = (1.224744871391589*fl[44]+0.7071067811865475*fl[41])*wv+(0.3162277660168379*fl[27]+0.1825741858350554*fl[22])*dv; 
  Ghat_l[22] = (1.224744871391589*fl[45]+0.7071067811865475*fl[42])*wv+(0.3162277660168379*fl[29]+0.1825741858350554*fl[25])*dv; 
  Ghat_l[23] = (1.224744871391589*fl[47]+0.7071067811865475*fl[46])*wv+(0.3162277660168379*fl[31]+0.1825741858350554*fl[30])*dv; 
  Ghat_l[24] = (1.224744871391589*fl[49]+0.7071067811865475*fl[48])*wv+(0.3535533905932737*fl[54]+0.2041241452319315*fl[51])*dv; 
  Ghat_l[25] = (1.224744871391589*fl[53]+0.7071067811865475*fl[50])*wv+(0.3535533905932737*fl[59]+0.2041241452319315*fl[55])*dv; 
  Ghat_l[26] = (1.224744871391589*fl[54]+0.7071067811865475*fl[51])*wv+(0.3535533905932737*fl[49]+0.2041241452319315*fl[48])*dv; 
  Ghat_l[27] = (1.224744871391589*fl[56]+0.7071067811865475*fl[52])*wv+(0.3535533905932737*fl[61]+0.2041241452319315*fl[58])*dv; 
  Ghat_l[28] = (1.224744871391589*fl[59]+0.7071067811865475*fl[55])*wv+(0.3535533905932737*fl[53]+0.2041241452319315*fl[50])*dv; 
  Ghat_l[29] = (1.224744871391589*fl[60]+0.7071067811865475*fl[57])*wv+(0.3535533905932737*fl[63]+0.2041241452319315*fl[62])*dv; 
  Ghat_l[30] = (1.224744871391589*fl[61]+0.7071067811865475*fl[58])*wv+(0.3535533905932737*fl[56]+0.2041241452319315*fl[52])*dv; 
  Ghat_l[31] = (1.224744871391589*fl[63]+0.7071067811865475*fl[62])*wv+(0.3535533905932737*fl[60]+0.2041241452319315*fl[57])*dv; 
  Ghat_l[32] = (1.224744871391589*fl[65]+0.7071067811865475*fl[64])*wv+(0.3535533905932737*fl[70]+0.2041241452319315*fl[67])*dv; 
  Ghat_l[33] = (1.224744871391589*fl[69]+0.7071067811865475*fl[66])*wv+(0.3535533905932737*fl[75]+0.2041241452319315*fl[71])*dv; 
  Ghat_l[34] = (1.224744871391589*fl[70]+0.7071067811865475*fl[67])*wv+(0.3535533905932737*fl[65]+0.2041241452319315*fl[64])*dv; 
  Ghat_l[35] = (1.224744871391589*fl[72]+0.7071067811865475*fl[68])*wv+(0.3535533905932737*fl[77]+0.2041241452319315*fl[74])*dv; 
  Ghat_l[36] = (1.224744871391589*fl[75]+0.7071067811865475*fl[71])*wv+(0.3535533905932737*fl[69]+0.2041241452319315*fl[66])*dv; 
  Ghat_l[37] = (1.224744871391589*fl[76]+0.7071067811865475*fl[73])*wv+(0.3535533905932737*fl[79]+0.2041241452319315*fl[78])*dv; 
  Ghat_l[38] = (1.224744871391589*fl[77]+0.7071067811865475*fl[74])*wv+(0.3535533905932737*fl[72]+0.2041241452319315*fl[68])*dv; 
  Ghat_l[39] = (1.224744871391589*fl[79]+0.7071067811865475*fl[78])*wv+(0.3535533905932737*fl[76]+0.2041241452319315*fl[73])*dv; 

  } else { 

  Ghat_r[0] = -0.08333333333333333*((14.69693845669907*fr[1]-8.485281374238571*fr[0])*wv+(4.242640687119286*fr[7]-2.449489742783178*fr[3])*dv); 
  Ghat_r[1] = -0.08333333333333333*((14.69693845669907*fr[6]-8.485281374238571*fr[2])*wv+(4.242640687119286*fr[16]-2.449489742783178*fr[8])*dv); 
  Ghat_r[2] = -0.01666666666666667*((73.48469228349535*fr[7]-42.42640687119286*fr[3])*wv+(18.97366596101028*fr[33]-10.95445115010332*fr[32]+21.21320343559643*fr[1]-12.24744871391589*fr[0])*dv); 
  Ghat_r[3] = -0.08333333333333333*((14.69693845669907*fr[9]-8.485281374238571*fr[4])*wv+(4.242640687119286*fr[18]-2.449489742783178*fr[11])*dv); 
  Ghat_r[4] = -0.08333333333333333*((14.69693845669907*fr[12]-8.485281374238571*fr[5])*wv+(4.242640687119286*fr[21]-2.449489742783178*fr[14])*dv); 
  Ghat_r[5] = -0.01666666666666667*((73.48469228349535*fr[16]-42.42640687119286*fr[8])*wv+(18.97366596101028*fr[37]-10.95445115010333*fr[34]+21.21320343559643*fr[6]-12.24744871391589*fr[2])*dv); 
  Ghat_r[6] = -0.08333333333333333*((14.69693845669907*fr[17]-8.485281374238571*fr[10])*wv+(4.242640687119286*fr[26]-2.449489742783178*fr[19])*dv); 
  Ghat_r[7] = -0.01666666666666667*((73.48469228349535*fr[18]-42.42640687119286*fr[11])*wv+(18.97366596101028*fr[38]-10.95445115010333*fr[35]+21.21320343559643*fr[9]-12.24744871391589*fr[4])*dv); 
  Ghat_r[8] = -0.08333333333333333*((14.69693845669907*fr[20]-8.485281374238571*fr[13])*wv+(4.242640687119286*fr[27]-2.449489742783178*fr[22])*dv); 
  Ghat_r[9] = -0.01666666666666667*((73.48469228349535*fr[21]-42.42640687119286*fr[14])*wv+(18.97366596101028*fr[40]-10.95445115010333*fr[36]+21.21320343559643*fr[12]-12.24744871391589*fr[5])*dv); 
  Ghat_r[10] = -0.08333333333333333*((14.69693845669907*fr[23]-8.485281374238571*fr[15])*wv+(4.242640687119286*fr[29]-2.449489742783178*fr[25])*dv); 
  Ghat_r[11] = -0.01666666666666667*((73.48469228349535*fr[26]-42.42640687119286*fr[19])*wv+(18.97366596101028*fr[43]-10.95445115010332*fr[39]+21.21320343559643*fr[17]-12.24744871391589*fr[10])*dv); 
  Ghat_r[12] = -0.01666666666666667*((73.48469228349535*fr[27]-42.42640687119286*fr[22])*wv+(18.97366596101028*fr[44]-10.95445115010332*fr[41]+21.21320343559643*fr[20]-12.24744871391589*fr[13])*dv); 
  Ghat_r[13] = -0.08333333333333333*((14.69693845669907*fr[28]-8.485281374238571*fr[24])*wv+(4.242640687119286*fr[31]-2.449489742783178*fr[30])*dv); 
  Ghat_r[14] = -0.01666666666666667*((73.48469228349535*fr[29]-42.42640687119286*fr[25])*wv+(18.97366596101028*fr[45]-10.95445115010332*fr[42]+21.21320343559643*fr[23]-12.24744871391589*fr[15])*dv); 
  Ghat_r[15] = -0.01666666666666667*((73.48469228349535*fr[31]-42.42640687119286*fr[30])*wv+(18.97366596101028*fr[47]-10.95445115010333*fr[46]+21.21320343559643*fr[28]-12.24744871391589*fr[24])*dv); 
  Ghat_r[16] = -0.03333333333333333*((36.74234614174768*fr[33]-21.21320343559643*fr[32])*wv+(9.48683298050514*fr[7]-5.477225575051662*fr[3])*dv); 
  Ghat_r[17] = -0.03333333333333333*((36.74234614174768*fr[37]-21.21320343559643*fr[34])*wv+(9.48683298050514*fr[16]-5.477225575051662*fr[8])*dv); 
  Ghat_r[18] = -0.03333333333333333*((36.74234614174768*fr[38]-21.21320343559643*fr[35])*wv+(9.48683298050514*fr[18]-5.477225575051662*fr[11])*dv); 
  Ghat_r[19] = -0.03333333333333333*((36.74234614174768*fr[40]-21.21320343559643*fr[36])*wv+(9.48683298050514*fr[21]-5.477225575051662*fr[14])*dv); 
  Ghat_r[20] = -0.03333333333333333*((36.74234614174768*fr[43]-21.21320343559643*fr[39])*wv+(9.48683298050514*fr[26]-5.477225575051662*fr[19])*dv); 
  Ghat_r[21] = -0.03333333333333333*((36.74234614174768*fr[44]-21.21320343559643*fr[41])*wv+(9.48683298050514*fr[27]-5.477225575051662*fr[22])*dv); 
  Ghat_r[22] = -0.03333333333333333*((36.74234614174768*fr[45]-21.21320343559643*fr[42])*wv+(9.48683298050514*fr[29]-5.477225575051662*fr[25])*dv); 
  Ghat_r[23] = -0.03333333333333333*((36.74234614174768*fr[47]-21.21320343559643*fr[46])*wv+(9.48683298050514*fr[31]-5.477225575051662*fr[30])*dv); 
  Ghat_r[24] = -0.01666666666666667*((73.48469228349536*fr[49]-42.42640687119286*fr[48])*wv+(21.21320343559643*fr[54]-12.24744871391589*fr[51])*dv); 
  Ghat_r[25] = -0.01666666666666667*((73.48469228349536*fr[53]-42.42640687119286*fr[50])*wv+(21.21320343559643*fr[59]-12.24744871391589*fr[55])*dv); 
  Ghat_r[26] = -0.01666666666666667*((73.48469228349536*fr[54]-42.42640687119286*fr[51])*wv+(21.21320343559643*fr[49]-12.24744871391589*fr[48])*dv); 
  Ghat_r[27] = -0.01666666666666667*((73.48469228349536*fr[56]-42.42640687119286*fr[52])*wv+(21.21320343559643*fr[61]-12.24744871391589*fr[58])*dv); 
  Ghat_r[28] = -0.01666666666666667*((73.48469228349536*fr[59]-42.42640687119286*fr[55])*wv+(21.21320343559643*fr[53]-12.24744871391589*fr[50])*dv); 
  Ghat_r[29] = -0.01666666666666667*((73.48469228349536*fr[60]-42.42640687119286*fr[57])*wv+(21.21320343559643*fr[63]-12.24744871391589*fr[62])*dv); 
  Ghat_r[30] = -0.01666666666666667*((73.48469228349536*fr[61]-42.42640687119286*fr[58])*wv+(21.21320343559643*fr[56]-12.24744871391589*fr[52])*dv); 
  Ghat_r[31] = -0.01666666666666667*((73.48469228349536*fr[63]-42.42640687119286*fr[62])*wv+(21.21320343559643*fr[60]-12.24744871391589*fr[57])*dv); 
  Ghat_r[32] = -0.01666666666666667*((73.48469228349536*fr[65]-42.42640687119286*fr[64])*wv+(21.21320343559643*fr[70]-12.24744871391589*fr[67])*dv); 
  Ghat_r[33] = -0.01666666666666667*((73.48469228349536*fr[69]-42.42640687119286*fr[66])*wv+(21.21320343559643*fr[75]-12.24744871391589*fr[71])*dv); 
  Ghat_r[34] = -0.01666666666666667*((73.48469228349536*fr[70]-42.42640687119286*fr[67])*wv+(21.21320343559643*fr[65]-12.24744871391589*fr[64])*dv); 
  Ghat_r[35] = -0.01666666666666667*((73.48469228349536*fr[72]-42.42640687119286*fr[68])*wv+(21.21320343559643*fr[77]-12.24744871391589*fr[74])*dv); 
  Ghat_r[36] = -0.01666666666666667*((73.48469228349536*fr[75]-42.42640687119286*fr[71])*wv+(21.21320343559643*fr[69]-12.24744871391589*fr[66])*dv); 
  Ghat_r[37] = -0.01666666666666667*((73.48469228349536*fr[76]-42.42640687119286*fr[73])*wv+(21.21320343559643*fr[79]-12.24744871391589*fr[78])*dv); 
  Ghat_r[38] = -0.01666666666666667*((73.48469228349536*fr[77]-42.42640687119286*fr[74])*wv+(21.21320343559643*fr[72]-12.24744871391589*fr[68])*dv); 
  Ghat_r[39] = -0.01666666666666667*((73.48469228349536*fr[79]-42.42640687119286*fr[78])*wv+(21.21320343559643*fr[76]-12.24744871391589*fr[73])*dv); 

  Ghat_l[0] = -0.08333333333333333*((14.69693845669907*fc[1]-8.485281374238571*fc[0])*wv+(4.242640687119286*fc[7]-2.449489742783178*fc[3])*dv); 
  Ghat_l[1] = -0.08333333333333333*((14.69693845669907*fc[6]-8.485281374238571*fc[2])*wv+(4.242640687119286*fc[16]-2.449489742783178*fc[8])*dv); 
  Ghat_l[2] = -0.01666666666666667*((73.48469228349535*fc[7]-42.42640687119286*fc[3])*wv+(18.97366596101028*fc[33]-10.95445115010332*fc[32]+21.21320343559643*fc[1]-12.24744871391589*fc[0])*dv); 
  Ghat_l[3] = -0.08333333333333333*((14.69693845669907*fc[9]-8.485281374238571*fc[4])*wv+(4.242640687119286*fc[18]-2.449489742783178*fc[11])*dv); 
  Ghat_l[4] = -0.08333333333333333*((14.69693845669907*fc[12]-8.485281374238571*fc[5])*wv+(4.242640687119286*fc[21]-2.449489742783178*fc[14])*dv); 
  Ghat_l[5] = -0.01666666666666667*((73.48469228349535*fc[16]-42.42640687119286*fc[8])*wv+(18.97366596101028*fc[37]-10.95445115010333*fc[34]+21.21320343559643*fc[6]-12.24744871391589*fc[2])*dv); 
  Ghat_l[6] = -0.08333333333333333*((14.69693845669907*fc[17]-8.485281374238571*fc[10])*wv+(4.242640687119286*fc[26]-2.449489742783178*fc[19])*dv); 
  Ghat_l[7] = -0.01666666666666667*((73.48469228349535*fc[18]-42.42640687119286*fc[11])*wv+(18.97366596101028*fc[38]-10.95445115010333*fc[35]+21.21320343559643*fc[9]-12.24744871391589*fc[4])*dv); 
  Ghat_l[8] = -0.08333333333333333*((14.69693845669907*fc[20]-8.485281374238571*fc[13])*wv+(4.242640687119286*fc[27]-2.449489742783178*fc[22])*dv); 
  Ghat_l[9] = -0.01666666666666667*((73.48469228349535*fc[21]-42.42640687119286*fc[14])*wv+(18.97366596101028*fc[40]-10.95445115010333*fc[36]+21.21320343559643*fc[12]-12.24744871391589*fc[5])*dv); 
  Ghat_l[10] = -0.08333333333333333*((14.69693845669907*fc[23]-8.485281374238571*fc[15])*wv+(4.242640687119286*fc[29]-2.449489742783178*fc[25])*dv); 
  Ghat_l[11] = -0.01666666666666667*((73.48469228349535*fc[26]-42.42640687119286*fc[19])*wv+(18.97366596101028*fc[43]-10.95445115010332*fc[39]+21.21320343559643*fc[17]-12.24744871391589*fc[10])*dv); 
  Ghat_l[12] = -0.01666666666666667*((73.48469228349535*fc[27]-42.42640687119286*fc[22])*wv+(18.97366596101028*fc[44]-10.95445115010332*fc[41]+21.21320343559643*fc[20]-12.24744871391589*fc[13])*dv); 
  Ghat_l[13] = -0.08333333333333333*((14.69693845669907*fc[28]-8.485281374238571*fc[24])*wv+(4.242640687119286*fc[31]-2.449489742783178*fc[30])*dv); 
  Ghat_l[14] = -0.01666666666666667*((73.48469228349535*fc[29]-42.42640687119286*fc[25])*wv+(18.97366596101028*fc[45]-10.95445115010332*fc[42]+21.21320343559643*fc[23]-12.24744871391589*fc[15])*dv); 
  Ghat_l[15] = -0.01666666666666667*((73.48469228349535*fc[31]-42.42640687119286*fc[30])*wv+(18.97366596101028*fc[47]-10.95445115010333*fc[46]+21.21320343559643*fc[28]-12.24744871391589*fc[24])*dv); 
  Ghat_l[16] = -0.03333333333333333*((36.74234614174768*fc[33]-21.21320343559643*fc[32])*wv+(9.48683298050514*fc[7]-5.477225575051662*fc[3])*dv); 
  Ghat_l[17] = -0.03333333333333333*((36.74234614174768*fc[37]-21.21320343559643*fc[34])*wv+(9.48683298050514*fc[16]-5.477225575051662*fc[8])*dv); 
  Ghat_l[18] = -0.03333333333333333*((36.74234614174768*fc[38]-21.21320343559643*fc[35])*wv+(9.48683298050514*fc[18]-5.477225575051662*fc[11])*dv); 
  Ghat_l[19] = -0.03333333333333333*((36.74234614174768*fc[40]-21.21320343559643*fc[36])*wv+(9.48683298050514*fc[21]-5.477225575051662*fc[14])*dv); 
  Ghat_l[20] = -0.03333333333333333*((36.74234614174768*fc[43]-21.21320343559643*fc[39])*wv+(9.48683298050514*fc[26]-5.477225575051662*fc[19])*dv); 
  Ghat_l[21] = -0.03333333333333333*((36.74234614174768*fc[44]-21.21320343559643*fc[41])*wv+(9.48683298050514*fc[27]-5.477225575051662*fc[22])*dv); 
  Ghat_l[22] = -0.03333333333333333*((36.74234614174768*fc[45]-21.21320343559643*fc[42])*wv+(9.48683298050514*fc[29]-5.477225575051662*fc[25])*dv); 
  Ghat_l[23] = -0.03333333333333333*((36.74234614174768*fc[47]-21.21320343559643*fc[46])*wv+(9.48683298050514*fc[31]-5.477225575051662*fc[30])*dv); 
  Ghat_l[24] = -0.01666666666666667*((73.48469228349536*fc[49]-42.42640687119286*fc[48])*wv+(21.21320343559643*fc[54]-12.24744871391589*fc[51])*dv); 
  Ghat_l[25] = -0.01666666666666667*((73.48469228349536*fc[53]-42.42640687119286*fc[50])*wv+(21.21320343559643*fc[59]-12.24744871391589*fc[55])*dv); 
  Ghat_l[26] = -0.01666666666666667*((73.48469228349536*fc[54]-42.42640687119286*fc[51])*wv+(21.21320343559643*fc[49]-12.24744871391589*fc[48])*dv); 
  Ghat_l[27] = -0.01666666666666667*((73.48469228349536*fc[56]-42.42640687119286*fc[52])*wv+(21.21320343559643*fc[61]-12.24744871391589*fc[58])*dv); 
  Ghat_l[28] = -0.01666666666666667*((73.48469228349536*fc[59]-42.42640687119286*fc[55])*wv+(21.21320343559643*fc[53]-12.24744871391589*fc[50])*dv); 
  Ghat_l[29] = -0.01666666666666667*((73.48469228349536*fc[60]-42.42640687119286*fc[57])*wv+(21.21320343559643*fc[63]-12.24744871391589*fc[62])*dv); 
  Ghat_l[30] = -0.01666666666666667*((73.48469228349536*fc[61]-42.42640687119286*fc[58])*wv+(21.21320343559643*fc[56]-12.24744871391589*fc[52])*dv); 
  Ghat_l[31] = -0.01666666666666667*((73.48469228349536*fc[63]-42.42640687119286*fc[62])*wv+(21.21320343559643*fc[60]-12.24744871391589*fc[57])*dv); 
  Ghat_l[32] = -0.01666666666666667*((73.48469228349536*fc[65]-42.42640687119286*fc[64])*wv+(21.21320343559643*fc[70]-12.24744871391589*fc[67])*dv); 
  Ghat_l[33] = -0.01666666666666667*((73.48469228349536*fc[69]-42.42640687119286*fc[66])*wv+(21.21320343559643*fc[75]-12.24744871391589*fc[71])*dv); 
  Ghat_l[34] = -0.01666666666666667*((73.48469228349536*fc[70]-42.42640687119286*fc[67])*wv+(21.21320343559643*fc[65]-12.24744871391589*fc[64])*dv); 
  Ghat_l[35] = -0.01666666666666667*((73.48469228349536*fc[72]-42.42640687119286*fc[68])*wv+(21.21320343559643*fc[77]-12.24744871391589*fc[74])*dv); 
  Ghat_l[36] = -0.01666666666666667*((73.48469228349536*fc[75]-42.42640687119286*fc[71])*wv+(21.21320343559643*fc[69]-12.24744871391589*fc[66])*dv); 
  Ghat_l[37] = -0.01666666666666667*((73.48469228349536*fc[76]-42.42640687119286*fc[73])*wv+(21.21320343559643*fc[79]-12.24744871391589*fc[78])*dv); 
  Ghat_l[38] = -0.01666666666666667*((73.48469228349536*fc[77]-42.42640687119286*fc[74])*wv+(21.21320343559643*fc[72]-12.24744871391589*fc[68])*dv); 
  Ghat_l[39] = -0.01666666666666667*((73.48469228349536*fc[79]-42.42640687119286*fc[78])*wv+(21.21320343559643*fc[76]-12.24744871391589*fc[73])*dv); 

  } 
  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dx10; 
  out[1] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dx10; 
  out[2] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dx10; 
  out[3] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dx10; 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dx10; 
  out[5] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dx10; 
  out[6] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dx10; 
  out[7] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dx10; 
  out[8] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dx10; 
  out[9] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dx10; 
  out[10] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dx10; 
  out[11] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dx10; 
  out[12] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dx10; 
  out[13] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*dx10; 
  out[14] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*dx10; 
  out[15] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*dx10; 
  out[16] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dx10; 
  out[17] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dx10; 
  out[18] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dx10; 
  out[19] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*dx10; 
  out[20] += -1.224744871391589*(Ghat_r[8]+Ghat_l[8])*dx10; 
  out[21] += -1.224744871391589*(Ghat_r[9]+Ghat_l[9])*dx10; 
  out[22] += (0.7071067811865475*Ghat_l[12]-0.7071067811865475*Ghat_r[12])*dx10; 
  out[23] += -1.224744871391589*(Ghat_r[10]+Ghat_l[10])*dx10; 
  out[24] += (0.7071067811865475*Ghat_l[13]-0.7071067811865475*Ghat_r[13])*dx10; 
  out[25] += (0.7071067811865475*Ghat_l[14]-0.7071067811865475*Ghat_r[14])*dx10; 
  out[26] += -1.224744871391589*(Ghat_r[11]+Ghat_l[11])*dx10; 
  out[27] += -1.224744871391589*(Ghat_r[12]+Ghat_l[12])*dx10; 
  out[28] += -1.224744871391589*(Ghat_r[13]+Ghat_l[13])*dx10; 
  out[29] += -1.224744871391589*(Ghat_r[14]+Ghat_l[14])*dx10; 
  out[30] += (0.7071067811865475*Ghat_l[15]-0.7071067811865475*Ghat_r[15])*dx10; 
  out[31] += -1.224744871391589*(Ghat_r[15]+Ghat_l[15])*dx10; 
  out[32] += (0.7071067811865475*Ghat_l[16]-0.7071067811865475*Ghat_r[16])*dx10; 
  out[33] += -1.224744871391589*(Ghat_r[16]+Ghat_l[16])*dx10; 
  out[34] += (0.7071067811865475*Ghat_l[17]-0.7071067811865475*Ghat_r[17])*dx10; 
  out[35] += (0.7071067811865475*Ghat_l[18]-0.7071067811865475*Ghat_r[18])*dx10; 
  out[36] += (0.7071067811865475*Ghat_l[19]-0.7071067811865475*Ghat_r[19])*dx10; 
  out[37] += -1.224744871391589*(Ghat_r[17]+Ghat_l[17])*dx10; 
  out[38] += -1.224744871391589*(Ghat_r[18]+Ghat_l[18])*dx10; 
  out[39] += (0.7071067811865475*Ghat_l[20]-0.7071067811865475*Ghat_r[20])*dx10; 
  out[40] += -1.224744871391589*(Ghat_r[19]+Ghat_l[19])*dx10; 
  out[41] += (0.7071067811865475*Ghat_l[21]-0.7071067811865475*Ghat_r[21])*dx10; 
  out[42] += (0.7071067811865475*Ghat_l[22]-0.7071067811865475*Ghat_r[22])*dx10; 
  out[43] += -1.224744871391589*(Ghat_r[20]+Ghat_l[20])*dx10; 
  out[44] += -1.224744871391589*(Ghat_r[21]+Ghat_l[21])*dx10; 
  out[45] += -1.224744871391589*(Ghat_r[22]+Ghat_l[22])*dx10; 
  out[46] += (0.7071067811865475*Ghat_l[23]-0.7071067811865475*Ghat_r[23])*dx10; 
  out[47] += -1.224744871391589*(Ghat_r[23]+Ghat_l[23])*dx10; 
  out[48] += (0.7071067811865475*Ghat_l[24]-0.7071067811865475*Ghat_r[24])*dx10; 
  out[49] += -1.224744871391589*(Ghat_r[24]+Ghat_l[24])*dx10; 
  out[50] += (0.7071067811865475*Ghat_l[25]-0.7071067811865475*Ghat_r[25])*dx10; 
  out[51] += (0.7071067811865475*Ghat_l[26]-0.7071067811865475*Ghat_r[26])*dx10; 
  out[52] += (0.7071067811865475*Ghat_l[27]-0.7071067811865475*Ghat_r[27])*dx10; 
  out[53] += -1.224744871391589*(Ghat_r[25]+Ghat_l[25])*dx10; 
  out[54] += -1.224744871391589*(Ghat_r[26]+Ghat_l[26])*dx10; 
  out[55] += (0.7071067811865475*Ghat_l[28]-0.7071067811865475*Ghat_r[28])*dx10; 
  out[56] += -1.224744871391589*(Ghat_r[27]+Ghat_l[27])*dx10; 
  out[57] += (0.7071067811865475*Ghat_l[29]-0.7071067811865475*Ghat_r[29])*dx10; 
  out[58] += (0.7071067811865475*Ghat_l[30]-0.7071067811865475*Ghat_r[30])*dx10; 
  out[59] += -1.224744871391589*(Ghat_r[28]+Ghat_l[28])*dx10; 
  out[60] += -1.224744871391589*(Ghat_r[29]+Ghat_l[29])*dx10; 
  out[61] += -1.224744871391589*(Ghat_r[30]+Ghat_l[30])*dx10; 
  out[62] += (0.7071067811865475*Ghat_l[31]-0.7071067811865475*Ghat_r[31])*dx10; 
  out[63] += -1.224744871391589*(Ghat_r[31]+Ghat_l[31])*dx10; 
  out[64] += (0.7071067811865475*Ghat_l[32]-0.7071067811865475*Ghat_r[32])*dx10; 
  out[65] += -1.224744871391589*(Ghat_r[32]+Ghat_l[32])*dx10; 
  out[66] += (0.7071067811865475*Ghat_l[33]-0.7071067811865475*Ghat_r[33])*dx10; 
  out[67] += (0.7071067811865475*Ghat_l[34]-0.7071067811865475*Ghat_r[34])*dx10; 
  out[68] += (0.7071067811865475*Ghat_l[35]-0.7071067811865475*Ghat_r[35])*dx10; 
  out[69] += -1.224744871391589*(Ghat_r[33]+Ghat_l[33])*dx10; 
  out[70] += -1.224744871391589*(Ghat_r[34]+Ghat_l[34])*dx10; 
  out[71] += (0.7071067811865475*Ghat_l[36]-0.7071067811865475*Ghat_r[36])*dx10; 
  out[72] += -1.224744871391589*(Ghat_r[35]+Ghat_l[35])*dx10; 
  out[73] += (0.7071067811865475*Ghat_l[37]-0.7071067811865475*Ghat_r[37])*dx10; 
  out[74] += (0.7071067811865475*Ghat_l[38]-0.7071067811865475*Ghat_r[38])*dx10; 
  out[75] += -1.224744871391589*(Ghat_r[36]+Ghat_l[36])*dx10; 
  out[76] += -1.224744871391589*(Ghat_r[37]+Ghat_l[37])*dx10; 
  out[77] += -1.224744871391589*(Ghat_r[38]+Ghat_l[38])*dx10; 
  out[78] += (0.7071067811865475*Ghat_l[39]-0.7071067811865475*Ghat_r[39])*dx10; 
  out[79] += -1.224744871391589*(Ghat_r[39]+Ghat_l[39])*dx10; 
} 
