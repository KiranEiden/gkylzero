#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_boundary_surfx_2x3v_tensor_p1(const double *w, const double *dxv, const double *alpha_geo, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:     Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // alpha_geo:   Fields used only for general geometry.
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Incremented distribution function in center cell.
  const double dx10 = 2/dxv[0]; 
  const double dv = dxv[2], wv = w[2]; 
  double Ghat[40]; 

  if (edge == -1) { 

  if (wv>0) { 

  Ghat[0] = (1.224744871391589*fSkin[1]+0.7071067811865475*fSkin[0])*wv+(0.3535533905932737*fSkin[7]+0.2041241452319315*fSkin[3])*dv; 
  Ghat[1] = (1.224744871391589*fSkin[6]+0.7071067811865475*fSkin[2])*wv+(0.3535533905932737*fSkin[16]+0.2041241452319315*fSkin[8])*dv; 
  Ghat[2] = (1.224744871391589*fSkin[7]+0.7071067811865475*fSkin[3])*wv+(0.3535533905932737*fSkin[1]+0.2041241452319315*fSkin[0])*dv; 
  Ghat[3] = (1.224744871391589*fSkin[9]+0.7071067811865475*fSkin[4])*wv+(0.3535533905932737*fSkin[18]+0.2041241452319315*fSkin[11])*dv; 
  Ghat[4] = (1.224744871391589*fSkin[12]+0.7071067811865475*fSkin[5])*wv+(0.3535533905932737*fSkin[21]+0.2041241452319315*fSkin[14])*dv; 
  Ghat[5] = (1.224744871391589*fSkin[16]+0.7071067811865475*fSkin[8])*wv+(0.3535533905932737*fSkin[6]+0.2041241452319315*fSkin[2])*dv; 
  Ghat[6] = (1.224744871391589*fSkin[17]+0.7071067811865475*fSkin[10])*wv+(0.3535533905932737*fSkin[26]+0.2041241452319315*fSkin[19])*dv; 
  Ghat[7] = (1.224744871391589*fSkin[18]+0.7071067811865475*fSkin[11])*wv+(0.3535533905932737*fSkin[9]+0.2041241452319315*fSkin[4])*dv; 
  Ghat[8] = (1.224744871391589*fSkin[20]+0.7071067811865475*fSkin[13])*wv+(0.3535533905932737*fSkin[27]+0.2041241452319315*fSkin[22])*dv; 
  Ghat[9] = (1.224744871391589*fSkin[21]+0.7071067811865475*fSkin[14])*wv+(0.3535533905932737*fSkin[12]+0.2041241452319315*fSkin[5])*dv; 
  Ghat[10] = (1.224744871391589*fSkin[23]+0.7071067811865475*fSkin[15])*wv+(0.3535533905932737*fSkin[29]+0.2041241452319315*fSkin[25])*dv; 
  Ghat[11] = (1.224744871391589*fSkin[26]+0.7071067811865475*fSkin[19])*wv+(0.3535533905932737*fSkin[17]+0.2041241452319315*fSkin[10])*dv; 
  Ghat[12] = (1.224744871391589*fSkin[27]+0.7071067811865475*fSkin[22])*wv+(0.3535533905932737*fSkin[20]+0.2041241452319315*fSkin[13])*dv; 
  Ghat[13] = (1.224744871391589*fSkin[28]+0.7071067811865475*fSkin[24])*wv+(0.3535533905932737*fSkin[31]+0.2041241452319315*fSkin[30])*dv; 
  Ghat[14] = (1.224744871391589*fSkin[29]+0.7071067811865475*fSkin[25])*wv+(0.3535533905932737*fSkin[23]+0.2041241452319315*fSkin[15])*dv; 
  Ghat[15] = (1.224744871391589*fSkin[31]+0.7071067811865475*fSkin[30])*wv+(0.3535533905932737*fSkin[28]+0.2041241452319315*fSkin[24])*dv; 
  Ghat[16] = (0.3162277660168379*fSkin[7]+0.1825741858350554*fSkin[3])*dv; 
  Ghat[17] = (0.3162277660168379*fSkin[16]+0.1825741858350554*fSkin[8])*dv; 
  Ghat[18] = (0.3162277660168379*fSkin[18]+0.1825741858350554*fSkin[11])*dv; 
  Ghat[19] = (0.3162277660168379*fSkin[21]+0.1825741858350554*fSkin[14])*dv; 
  Ghat[20] = (0.3162277660168379*fSkin[26]+0.1825741858350554*fSkin[19])*dv; 
  Ghat[21] = (0.3162277660168379*fSkin[27]+0.1825741858350554*fSkin[22])*dv; 
  Ghat[22] = (0.3162277660168379*fSkin[29]+0.1825741858350554*fSkin[25])*dv; 
  Ghat[23] = (0.3162277660168379*fSkin[31]+0.1825741858350554*fSkin[30])*dv; 

  } else { 

  Ghat[0] = -0.08333333333333333*((14.69693845669907*fEdge[1]-8.485281374238571*fEdge[0])*wv+(4.242640687119286*fEdge[7]-2.449489742783178*fEdge[3])*dv); 
  Ghat[1] = -0.08333333333333333*((14.69693845669907*fEdge[6]-8.485281374238571*fEdge[2])*wv+(4.242640687119286*fEdge[16]-2.449489742783178*fEdge[8])*dv); 
  Ghat[2] = -0.08333333333333333*((14.69693845669907*fEdge[7]-8.485281374238571*fEdge[3])*wv+(4.242640687119286*fEdge[1]-2.449489742783178*fEdge[0])*dv); 
  Ghat[3] = -0.08333333333333333*((14.69693845669907*fEdge[9]-8.485281374238571*fEdge[4])*wv+(4.242640687119286*fEdge[18]-2.449489742783178*fEdge[11])*dv); 
  Ghat[4] = -0.08333333333333333*((14.69693845669907*fEdge[12]-8.485281374238571*fEdge[5])*wv+(4.242640687119286*fEdge[21]-2.449489742783178*fEdge[14])*dv); 
  Ghat[5] = -0.08333333333333333*((14.69693845669907*fEdge[16]-8.485281374238571*fEdge[8])*wv+(4.242640687119286*fEdge[6]-2.449489742783178*fEdge[2])*dv); 
  Ghat[6] = -0.08333333333333333*((14.69693845669907*fEdge[17]-8.485281374238571*fEdge[10])*wv+(4.242640687119286*fEdge[26]-2.449489742783178*fEdge[19])*dv); 
  Ghat[7] = -0.08333333333333333*((14.69693845669907*fEdge[18]-8.485281374238571*fEdge[11])*wv+(4.242640687119286*fEdge[9]-2.449489742783178*fEdge[4])*dv); 
  Ghat[8] = -0.08333333333333333*((14.69693845669907*fEdge[20]-8.485281374238571*fEdge[13])*wv+(4.242640687119286*fEdge[27]-2.449489742783178*fEdge[22])*dv); 
  Ghat[9] = -0.08333333333333333*((14.69693845669907*fEdge[21]-8.485281374238571*fEdge[14])*wv+(4.242640687119286*fEdge[12]-2.449489742783178*fEdge[5])*dv); 
  Ghat[10] = -0.08333333333333333*((14.69693845669907*fEdge[23]-8.485281374238571*fEdge[15])*wv+(4.242640687119286*fEdge[29]-2.449489742783178*fEdge[25])*dv); 
  Ghat[11] = -0.08333333333333333*((14.69693845669907*fEdge[26]-8.485281374238571*fEdge[19])*wv+(4.242640687119286*fEdge[17]-2.449489742783178*fEdge[10])*dv); 
  Ghat[12] = -0.08333333333333333*((14.69693845669907*fEdge[27]-8.485281374238571*fEdge[22])*wv+(4.242640687119286*fEdge[20]-2.449489742783178*fEdge[13])*dv); 
  Ghat[13] = -0.08333333333333333*((14.69693845669907*fEdge[28]-8.485281374238571*fEdge[24])*wv+(4.242640687119286*fEdge[31]-2.449489742783178*fEdge[30])*dv); 
  Ghat[14] = -0.08333333333333333*((14.69693845669907*fEdge[29]-8.485281374238571*fEdge[25])*wv+(4.242640687119286*fEdge[23]-2.449489742783178*fEdge[15])*dv); 
  Ghat[15] = -0.08333333333333333*((14.69693845669907*fEdge[31]-8.485281374238571*fEdge[30])*wv+(4.242640687119286*fEdge[28]-2.449489742783178*fEdge[24])*dv); 
  Ghat[16] = -0.03333333333333333*(9.48683298050514*fEdge[7]-5.477225575051662*fEdge[3])*dv; 
  Ghat[17] = -0.03333333333333333*(9.48683298050514*fEdge[16]-5.477225575051662*fEdge[8])*dv; 
  Ghat[18] = -0.03333333333333333*(9.48683298050514*fEdge[18]-5.477225575051662*fEdge[11])*dv; 
  Ghat[19] = -0.03333333333333333*(9.48683298050514*fEdge[21]-5.477225575051662*fEdge[14])*dv; 
  Ghat[20] = -0.03333333333333333*(9.48683298050514*fEdge[26]-5.477225575051662*fEdge[19])*dv; 
  Ghat[21] = -0.03333333333333333*(9.48683298050514*fEdge[27]-5.477225575051662*fEdge[22])*dv; 
  Ghat[22] = -0.03333333333333333*(9.48683298050514*fEdge[29]-5.477225575051662*fEdge[25])*dv; 
  Ghat[23] = -0.03333333333333333*(9.48683298050514*fEdge[31]-5.477225575051662*fEdge[30])*dv; 

  } 

  out[0] += -0.7071067811865475*Ghat[0]*dx10; 
  out[1] += -1.224744871391589*Ghat[0]*dx10; 
  out[2] += -0.7071067811865475*Ghat[1]*dx10; 
  out[3] += -0.7071067811865475*Ghat[2]*dx10; 
  out[4] += -0.7071067811865475*Ghat[3]*dx10; 
  out[5] += -0.7071067811865475*Ghat[4]*dx10; 
  out[6] += -1.224744871391589*Ghat[1]*dx10; 
  out[7] += -1.224744871391589*Ghat[2]*dx10; 
  out[8] += -0.7071067811865475*Ghat[5]*dx10; 
  out[9] += -1.224744871391589*Ghat[3]*dx10; 
  out[10] += -0.7071067811865475*Ghat[6]*dx10; 
  out[11] += -0.7071067811865475*Ghat[7]*dx10; 
  out[12] += -1.224744871391589*Ghat[4]*dx10; 
  out[13] += -0.7071067811865475*Ghat[8]*dx10; 
  out[14] += -0.7071067811865475*Ghat[9]*dx10; 
  out[15] += -0.7071067811865475*Ghat[10]*dx10; 
  out[16] += -1.224744871391589*Ghat[5]*dx10; 
  out[17] += -1.224744871391589*Ghat[6]*dx10; 
  out[18] += -1.224744871391589*Ghat[7]*dx10; 
  out[19] += -0.7071067811865475*Ghat[11]*dx10; 
  out[20] += -1.224744871391589*Ghat[8]*dx10; 
  out[21] += -1.224744871391589*Ghat[9]*dx10; 
  out[22] += -0.7071067811865475*Ghat[12]*dx10; 
  out[23] += -1.224744871391589*Ghat[10]*dx10; 
  out[24] += -0.7071067811865475*Ghat[13]*dx10; 
  out[25] += -0.7071067811865475*Ghat[14]*dx10; 
  out[26] += -1.224744871391589*Ghat[11]*dx10; 
  out[27] += -1.224744871391589*Ghat[12]*dx10; 
  out[28] += -1.224744871391589*Ghat[13]*dx10; 
  out[29] += -1.224744871391589*Ghat[14]*dx10; 
  out[30] += -0.7071067811865475*Ghat[15]*dx10; 
  out[31] += -1.224744871391589*Ghat[15]*dx10; 

  } else { 

  if (wv>0) { 

  Ghat[0] = (1.224744871391589*fEdge[1]+0.7071067811865475*fEdge[0])*wv+(0.3535533905932737*fEdge[7]+0.2041241452319315*fEdge[3])*dv; 
  Ghat[1] = (1.224744871391589*fEdge[6]+0.7071067811865475*fEdge[2])*wv+(0.3535533905932737*fEdge[16]+0.2041241452319315*fEdge[8])*dv; 
  Ghat[2] = (1.224744871391589*fEdge[7]+0.7071067811865475*fEdge[3])*wv+(0.3535533905932737*fEdge[1]+0.2041241452319315*fEdge[0])*dv; 
  Ghat[3] = (1.224744871391589*fEdge[9]+0.7071067811865475*fEdge[4])*wv+(0.3535533905932737*fEdge[18]+0.2041241452319315*fEdge[11])*dv; 
  Ghat[4] = (1.224744871391589*fEdge[12]+0.7071067811865475*fEdge[5])*wv+(0.3535533905932737*fEdge[21]+0.2041241452319315*fEdge[14])*dv; 
  Ghat[5] = (1.224744871391589*fEdge[16]+0.7071067811865475*fEdge[8])*wv+(0.3535533905932737*fEdge[6]+0.2041241452319315*fEdge[2])*dv; 
  Ghat[6] = (1.224744871391589*fEdge[17]+0.7071067811865475*fEdge[10])*wv+(0.3535533905932737*fEdge[26]+0.2041241452319315*fEdge[19])*dv; 
  Ghat[7] = (1.224744871391589*fEdge[18]+0.7071067811865475*fEdge[11])*wv+(0.3535533905932737*fEdge[9]+0.2041241452319315*fEdge[4])*dv; 
  Ghat[8] = (1.224744871391589*fEdge[20]+0.7071067811865475*fEdge[13])*wv+(0.3535533905932737*fEdge[27]+0.2041241452319315*fEdge[22])*dv; 
  Ghat[9] = (1.224744871391589*fEdge[21]+0.7071067811865475*fEdge[14])*wv+(0.3535533905932737*fEdge[12]+0.2041241452319315*fEdge[5])*dv; 
  Ghat[10] = (1.224744871391589*fEdge[23]+0.7071067811865475*fEdge[15])*wv+(0.3535533905932737*fEdge[29]+0.2041241452319315*fEdge[25])*dv; 
  Ghat[11] = (1.224744871391589*fEdge[26]+0.7071067811865475*fEdge[19])*wv+(0.3535533905932737*fEdge[17]+0.2041241452319315*fEdge[10])*dv; 
  Ghat[12] = (1.224744871391589*fEdge[27]+0.7071067811865475*fEdge[22])*wv+(0.3535533905932737*fEdge[20]+0.2041241452319315*fEdge[13])*dv; 
  Ghat[13] = (1.224744871391589*fEdge[28]+0.7071067811865475*fEdge[24])*wv+(0.3535533905932737*fEdge[31]+0.2041241452319315*fEdge[30])*dv; 
  Ghat[14] = (1.224744871391589*fEdge[29]+0.7071067811865475*fEdge[25])*wv+(0.3535533905932737*fEdge[23]+0.2041241452319315*fEdge[15])*dv; 
  Ghat[15] = (1.224744871391589*fEdge[31]+0.7071067811865475*fEdge[30])*wv+(0.3535533905932737*fEdge[28]+0.2041241452319315*fEdge[24])*dv; 
  Ghat[16] = (0.3162277660168379*fEdge[7]+0.1825741858350554*fEdge[3])*dv; 
  Ghat[17] = (0.3162277660168379*fEdge[16]+0.1825741858350554*fEdge[8])*dv; 
  Ghat[18] = (0.3162277660168379*fEdge[18]+0.1825741858350554*fEdge[11])*dv; 
  Ghat[19] = (0.3162277660168379*fEdge[21]+0.1825741858350554*fEdge[14])*dv; 
  Ghat[20] = (0.3162277660168379*fEdge[26]+0.1825741858350554*fEdge[19])*dv; 
  Ghat[21] = (0.3162277660168379*fEdge[27]+0.1825741858350554*fEdge[22])*dv; 
  Ghat[22] = (0.3162277660168379*fEdge[29]+0.1825741858350554*fEdge[25])*dv; 
  Ghat[23] = (0.3162277660168379*fEdge[31]+0.1825741858350554*fEdge[30])*dv; 

  } else { 

  Ghat[0] = -0.08333333333333333*((14.69693845669907*fSkin[1]-8.485281374238571*fSkin[0])*wv+(4.242640687119286*fSkin[7]-2.449489742783178*fSkin[3])*dv); 
  Ghat[1] = -0.08333333333333333*((14.69693845669907*fSkin[6]-8.485281374238571*fSkin[2])*wv+(4.242640687119286*fSkin[16]-2.449489742783178*fSkin[8])*dv); 
  Ghat[2] = -0.08333333333333333*((14.69693845669907*fSkin[7]-8.485281374238571*fSkin[3])*wv+(4.242640687119286*fSkin[1]-2.449489742783178*fSkin[0])*dv); 
  Ghat[3] = -0.08333333333333333*((14.69693845669907*fSkin[9]-8.485281374238571*fSkin[4])*wv+(4.242640687119286*fSkin[18]-2.449489742783178*fSkin[11])*dv); 
  Ghat[4] = -0.08333333333333333*((14.69693845669907*fSkin[12]-8.485281374238571*fSkin[5])*wv+(4.242640687119286*fSkin[21]-2.449489742783178*fSkin[14])*dv); 
  Ghat[5] = -0.08333333333333333*((14.69693845669907*fSkin[16]-8.485281374238571*fSkin[8])*wv+(4.242640687119286*fSkin[6]-2.449489742783178*fSkin[2])*dv); 
  Ghat[6] = -0.08333333333333333*((14.69693845669907*fSkin[17]-8.485281374238571*fSkin[10])*wv+(4.242640687119286*fSkin[26]-2.449489742783178*fSkin[19])*dv); 
  Ghat[7] = -0.08333333333333333*((14.69693845669907*fSkin[18]-8.485281374238571*fSkin[11])*wv+(4.242640687119286*fSkin[9]-2.449489742783178*fSkin[4])*dv); 
  Ghat[8] = -0.08333333333333333*((14.69693845669907*fSkin[20]-8.485281374238571*fSkin[13])*wv+(4.242640687119286*fSkin[27]-2.449489742783178*fSkin[22])*dv); 
  Ghat[9] = -0.08333333333333333*((14.69693845669907*fSkin[21]-8.485281374238571*fSkin[14])*wv+(4.242640687119286*fSkin[12]-2.449489742783178*fSkin[5])*dv); 
  Ghat[10] = -0.08333333333333333*((14.69693845669907*fSkin[23]-8.485281374238571*fSkin[15])*wv+(4.242640687119286*fSkin[29]-2.449489742783178*fSkin[25])*dv); 
  Ghat[11] = -0.08333333333333333*((14.69693845669907*fSkin[26]-8.485281374238571*fSkin[19])*wv+(4.242640687119286*fSkin[17]-2.449489742783178*fSkin[10])*dv); 
  Ghat[12] = -0.08333333333333333*((14.69693845669907*fSkin[27]-8.485281374238571*fSkin[22])*wv+(4.242640687119286*fSkin[20]-2.449489742783178*fSkin[13])*dv); 
  Ghat[13] = -0.08333333333333333*((14.69693845669907*fSkin[28]-8.485281374238571*fSkin[24])*wv+(4.242640687119286*fSkin[31]-2.449489742783178*fSkin[30])*dv); 
  Ghat[14] = -0.08333333333333333*((14.69693845669907*fSkin[29]-8.485281374238571*fSkin[25])*wv+(4.242640687119286*fSkin[23]-2.449489742783178*fSkin[15])*dv); 
  Ghat[15] = -0.08333333333333333*((14.69693845669907*fSkin[31]-8.485281374238571*fSkin[30])*wv+(4.242640687119286*fSkin[28]-2.449489742783178*fSkin[24])*dv); 
  Ghat[16] = -0.03333333333333333*(9.48683298050514*fSkin[7]-5.477225575051662*fSkin[3])*dv; 
  Ghat[17] = -0.03333333333333333*(9.48683298050514*fSkin[16]-5.477225575051662*fSkin[8])*dv; 
  Ghat[18] = -0.03333333333333333*(9.48683298050514*fSkin[18]-5.477225575051662*fSkin[11])*dv; 
  Ghat[19] = -0.03333333333333333*(9.48683298050514*fSkin[21]-5.477225575051662*fSkin[14])*dv; 
  Ghat[20] = -0.03333333333333333*(9.48683298050514*fSkin[26]-5.477225575051662*fSkin[19])*dv; 
  Ghat[21] = -0.03333333333333333*(9.48683298050514*fSkin[27]-5.477225575051662*fSkin[22])*dv; 
  Ghat[22] = -0.03333333333333333*(9.48683298050514*fSkin[29]-5.477225575051662*fSkin[25])*dv; 
  Ghat[23] = -0.03333333333333333*(9.48683298050514*fSkin[31]-5.477225575051662*fSkin[30])*dv; 

  } 

  out[0] += 0.7071067811865475*Ghat[0]*dx10; 
  out[1] += -1.224744871391589*Ghat[0]*dx10; 
  out[2] += 0.7071067811865475*Ghat[1]*dx10; 
  out[3] += 0.7071067811865475*Ghat[2]*dx10; 
  out[4] += 0.7071067811865475*Ghat[3]*dx10; 
  out[5] += 0.7071067811865475*Ghat[4]*dx10; 
  out[6] += -1.224744871391589*Ghat[1]*dx10; 
  out[7] += -1.224744871391589*Ghat[2]*dx10; 
  out[8] += 0.7071067811865475*Ghat[5]*dx10; 
  out[9] += -1.224744871391589*Ghat[3]*dx10; 
  out[10] += 0.7071067811865475*Ghat[6]*dx10; 
  out[11] += 0.7071067811865475*Ghat[7]*dx10; 
  out[12] += -1.224744871391589*Ghat[4]*dx10; 
  out[13] += 0.7071067811865475*Ghat[8]*dx10; 
  out[14] += 0.7071067811865475*Ghat[9]*dx10; 
  out[15] += 0.7071067811865475*Ghat[10]*dx10; 
  out[16] += -1.224744871391589*Ghat[5]*dx10; 
  out[17] += -1.224744871391589*Ghat[6]*dx10; 
  out[18] += -1.224744871391589*Ghat[7]*dx10; 
  out[19] += 0.7071067811865475*Ghat[11]*dx10; 
  out[20] += -1.224744871391589*Ghat[8]*dx10; 
  out[21] += -1.224744871391589*Ghat[9]*dx10; 
  out[22] += 0.7071067811865475*Ghat[12]*dx10; 
  out[23] += -1.224744871391589*Ghat[10]*dx10; 
  out[24] += 0.7071067811865475*Ghat[13]*dx10; 
  out[25] += 0.7071067811865475*Ghat[14]*dx10; 
  out[26] += -1.224744871391589*Ghat[11]*dx10; 
  out[27] += -1.224744871391589*Ghat[12]*dx10; 
  out[28] += -1.224744871391589*Ghat[13]*dx10; 
  out[29] += -1.224744871391589*Ghat[14]*dx10; 
  out[30] += 0.7071067811865475*Ghat[15]*dx10; 
  out[31] += -1.224744871391589*Ghat[15]*dx10; 

  } 
  return 0.;

} 