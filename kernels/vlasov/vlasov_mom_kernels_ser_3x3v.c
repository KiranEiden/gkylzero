#include <gkyl_vlasov_mom_kernels.h>

void vlasov_mom_3x3v_m0_p1(const double *w, const double *dxv, const long *idx, const double *f, double* restrict out) 
{ 
  const double volFact = dxv[3]*dxv[4]*dxv[5]/8; 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact; 
  out[2] += 2.828427124746191*f[2]*volFact; 
  out[3] += 2.828427124746191*f[3]*volFact; 
  out[4] += 2.828427124746191*f[7]*volFact; 
  out[5] += 2.828427124746191*f[8]*volFact; 
  out[6] += 2.828427124746191*f[9]*volFact; 
  out[7] += 2.828427124746191*f[22]*volFact; 
} 
void vlasov_mom_3x3v_m1i_p1(const double *w, const double *dxv, const long *idx, const double *f, double* restrict out) 
{ 
  const double volFact = dxv[3]*dxv[4]*dxv[5]/8; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx3 = w[5], dv3 = dxv[5]; 
  out[0] += volFact*(2.828427124746191*f[0]*wx1+0.8164965809277261*f[4]*dv1); 
  out[1] += volFact*(2.828427124746191*f[1]*wx1+0.8164965809277261*f[10]*dv1); 
  out[2] += volFact*(2.828427124746191*f[2]*wx1+0.8164965809277261*f[11]*dv1); 
  out[3] += volFact*(2.828427124746191*f[3]*wx1+0.8164965809277261*f[12]*dv1); 
  out[4] += volFact*(2.828427124746191*f[7]*wx1+0.8164965809277261*f[23]*dv1); 
  out[5] += volFact*(2.828427124746191*f[8]*wx1+0.8164965809277261*f[24]*dv1); 
  out[6] += volFact*(2.828427124746191*f[9]*wx1+0.8164965809277261*f[25]*dv1); 
  out[7] += volFact*(2.828427124746191*f[22]*wx1+0.8164965809277261*f[42]*dv1); 
  out[8] += volFact*(2.828427124746191*f[0]*wx2+0.8164965809277261*f[5]*dv2); 
  out[9] += volFact*(2.828427124746191*f[1]*wx2+0.8164965809277261*f[13]*dv2); 
  out[10] += volFact*(2.828427124746191*f[2]*wx2+0.8164965809277261*f[14]*dv2); 
  out[11] += volFact*(2.828427124746191*f[3]*wx2+0.8164965809277261*f[15]*dv2); 
  out[12] += volFact*(2.828427124746191*f[7]*wx2+0.8164965809277261*f[26]*dv2); 
  out[13] += volFact*(2.828427124746191*f[8]*wx2+0.8164965809277261*f[27]*dv2); 
  out[14] += volFact*(2.828427124746191*f[9]*wx2+0.8164965809277261*f[28]*dv2); 
  out[15] += volFact*(2.828427124746191*f[22]*wx2+0.8164965809277261*f[43]*dv2); 
  out[16] += volFact*(2.828427124746191*f[0]*wx3+0.8164965809277261*f[6]*dv3); 
  out[17] += volFact*(2.828427124746191*f[1]*wx3+0.8164965809277261*f[17]*dv3); 
  out[18] += volFact*(2.828427124746191*f[2]*wx3+0.8164965809277261*f[18]*dv3); 
  out[19] += volFact*(2.828427124746191*f[3]*wx3+0.8164965809277261*f[19]*dv3); 
  out[20] += volFact*(2.828427124746191*f[7]*wx3+0.8164965809277261*f[32]*dv3); 
  out[21] += volFact*(2.828427124746191*f[8]*wx3+0.8164965809277261*f[33]*dv3); 
  out[22] += volFact*(2.828427124746191*f[9]*wx3+0.8164965809277261*f[34]*dv3); 
  out[23] += volFact*(2.828427124746191*f[22]*wx3+0.8164965809277261*f[47]*dv3); 
} 
void vlasov_mom_3x3v_m2ij_p1(const double *w, const double *dxv, const long *idx, const double *f, double* restrict out) 
{ 
  const double volFact = dxv[3]*dxv[4]*dxv[5]/8; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[5], dv3 = dxv[5]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  out[0] += volFact*(2.828427124746191*f[0]*wx1_sq+1.632993161855453*f[4]*dv1*wx1+0.2357022603955158*f[0]*dv1_sq); 
  out[1] += volFact*(2.828427124746191*f[1]*wx1_sq+1.632993161855453*f[10]*dv1*wx1+0.2357022603955158*f[1]*dv1_sq); 
  out[2] += volFact*(2.828427124746191*f[2]*wx1_sq+1.632993161855453*f[11]*dv1*wx1+0.2357022603955158*f[2]*dv1_sq); 
  out[3] += volFact*(2.828427124746191*f[3]*wx1_sq+1.632993161855453*f[12]*dv1*wx1+0.2357022603955158*f[3]*dv1_sq); 
  out[4] += volFact*(2.828427124746191*f[7]*wx1_sq+1.632993161855453*f[23]*dv1*wx1+0.2357022603955158*f[7]*dv1_sq); 
  out[5] += volFact*(2.828427124746191*f[8]*wx1_sq+1.632993161855453*f[24]*dv1*wx1+0.2357022603955158*f[8]*dv1_sq); 
  out[6] += volFact*(2.828427124746191*f[9]*wx1_sq+1.632993161855453*f[25]*dv1*wx1+0.2357022603955158*f[9]*dv1_sq); 
  out[7] += volFact*(2.828427124746191*f[22]*wx1_sq+1.632993161855453*f[42]*dv1*wx1+0.2357022603955158*f[22]*dv1_sq); 
  out[8] += volFact*(2.828427124746191*f[0]*wx1*wx2+0.8164965809277261*f[4]*dv1*wx2+0.8164965809277261*f[5]*dv2*wx1+0.2357022603955158*f[16]*dv1*dv2); 
  out[9] += volFact*(2.828427124746191*f[1]*wx1*wx2+0.8164965809277261*f[10]*dv1*wx2+0.8164965809277261*f[13]*dv2*wx1+0.2357022603955158*f[29]*dv1*dv2); 
  out[10] += volFact*(2.828427124746191*f[2]*wx1*wx2+0.8164965809277261*f[11]*dv1*wx2+0.8164965809277261*f[14]*dv2*wx1+0.2357022603955158*f[30]*dv1*dv2); 
  out[11] += volFact*(2.828427124746191*f[3]*wx1*wx2+0.8164965809277261*f[12]*dv1*wx2+0.8164965809277261*f[15]*dv2*wx1+0.2357022603955158*f[31]*dv1*dv2); 
  out[12] += volFact*(2.828427124746191*f[7]*wx1*wx2+0.8164965809277261*f[23]*dv1*wx2+0.8164965809277261*f[26]*dv2*wx1+0.2357022603955158*f[44]*dv1*dv2); 
  out[13] += volFact*(2.828427124746191*f[8]*wx1*wx2+0.8164965809277261*f[24]*dv1*wx2+0.8164965809277261*f[27]*dv2*wx1+0.2357022603955158*f[45]*dv1*dv2); 
  out[14] += volFact*(2.828427124746191*f[9]*wx1*wx2+0.8164965809277261*f[25]*dv1*wx2+0.8164965809277261*f[28]*dv2*wx1+0.2357022603955158*f[46]*dv1*dv2); 
  out[15] += volFact*(2.828427124746191*f[22]*wx1*wx2+0.8164965809277261*f[42]*dv1*wx2+0.8164965809277261*f[43]*dv2*wx1+0.2357022603955158*f[57]*dv1*dv2); 
  out[16] += volFact*(2.828427124746191*f[0]*wx1*wx3+0.8164965809277261*f[4]*dv1*wx3+0.8164965809277261*f[6]*dv3*wx1+0.2357022603955158*f[20]*dv1*dv3); 
  out[17] += volFact*(2.828427124746191*f[1]*wx1*wx3+0.8164965809277261*f[10]*dv1*wx3+0.8164965809277261*f[17]*dv3*wx1+0.2357022603955158*f[35]*dv1*dv3); 
  out[18] += volFact*(2.828427124746191*f[2]*wx1*wx3+0.8164965809277261*f[11]*dv1*wx3+0.8164965809277261*f[18]*dv3*wx1+0.2357022603955158*f[36]*dv1*dv3); 
  out[19] += volFact*(2.828427124746191*f[3]*wx1*wx3+0.8164965809277261*f[12]*dv1*wx3+0.8164965809277261*f[19]*dv3*wx1+0.2357022603955158*f[37]*dv1*dv3); 
  out[20] += volFact*(2.828427124746191*f[7]*wx1*wx3+0.8164965809277261*f[23]*dv1*wx3+0.8164965809277261*f[32]*dv3*wx1+0.2357022603955158*f[48]*dv1*dv3); 
  out[21] += volFact*(2.828427124746191*f[8]*wx1*wx3+0.8164965809277261*f[24]*dv1*wx3+0.8164965809277261*f[33]*dv3*wx1+0.2357022603955158*f[49]*dv1*dv3); 
  out[22] += volFact*(2.828427124746191*f[9]*wx1*wx3+0.8164965809277261*f[25]*dv1*wx3+0.8164965809277261*f[34]*dv3*wx1+0.2357022603955158*f[50]*dv1*dv3); 
  out[23] += volFact*(2.828427124746191*f[22]*wx1*wx3+0.8164965809277261*f[42]*dv1*wx3+0.8164965809277261*f[47]*dv3*wx1+0.2357022603955158*f[58]*dv1*dv3); 
  out[24] += volFact*(2.828427124746191*f[0]*wx2_sq+1.632993161855453*f[5]*dv2*wx2+0.2357022603955158*f[0]*dv2_sq); 
  out[25] += volFact*(2.828427124746191*f[1]*wx2_sq+1.632993161855453*f[13]*dv2*wx2+0.2357022603955158*f[1]*dv2_sq); 
  out[26] += volFact*(2.828427124746191*f[2]*wx2_sq+1.632993161855453*f[14]*dv2*wx2+0.2357022603955158*f[2]*dv2_sq); 
  out[27] += volFact*(2.828427124746191*f[3]*wx2_sq+1.632993161855453*f[15]*dv2*wx2+0.2357022603955158*f[3]*dv2_sq); 
  out[28] += volFact*(2.828427124746191*f[7]*wx2_sq+1.632993161855453*f[26]*dv2*wx2+0.2357022603955158*f[7]*dv2_sq); 
  out[29] += volFact*(2.828427124746191*f[8]*wx2_sq+1.632993161855453*f[27]*dv2*wx2+0.2357022603955158*f[8]*dv2_sq); 
  out[30] += volFact*(2.828427124746191*f[9]*wx2_sq+1.632993161855453*f[28]*dv2*wx2+0.2357022603955158*f[9]*dv2_sq); 
  out[31] += volFact*(2.828427124746191*f[22]*wx2_sq+1.632993161855453*f[43]*dv2*wx2+0.2357022603955158*f[22]*dv2_sq); 
  out[32] += volFact*(2.828427124746191*f[0]*wx2*wx3+0.8164965809277261*f[5]*dv2*wx3+0.8164965809277261*f[6]*dv3*wx2+0.2357022603955158*f[21]*dv2*dv3); 
  out[33] += volFact*(2.828427124746191*f[1]*wx2*wx3+0.8164965809277261*f[13]*dv2*wx3+0.8164965809277261*f[17]*dv3*wx2+0.2357022603955158*f[38]*dv2*dv3); 
  out[34] += volFact*(2.828427124746191*f[2]*wx2*wx3+0.8164965809277261*f[14]*dv2*wx3+0.8164965809277261*f[18]*dv3*wx2+0.2357022603955158*f[39]*dv2*dv3); 
  out[35] += volFact*(2.828427124746191*f[3]*wx2*wx3+0.8164965809277261*f[15]*dv2*wx3+0.8164965809277261*f[19]*dv3*wx2+0.2357022603955158*f[40]*dv2*dv3); 
  out[36] += volFact*(2.828427124746191*f[7]*wx2*wx3+0.8164965809277261*f[26]*dv2*wx3+0.8164965809277261*f[32]*dv3*wx2+0.2357022603955158*f[51]*dv2*dv3); 
  out[37] += volFact*(2.828427124746191*f[8]*wx2*wx3+0.8164965809277261*f[27]*dv2*wx3+0.8164965809277261*f[33]*dv3*wx2+0.2357022603955158*f[52]*dv2*dv3); 
  out[38] += volFact*(2.828427124746191*f[9]*wx2*wx3+0.8164965809277261*f[28]*dv2*wx3+0.8164965809277261*f[34]*dv3*wx2+0.2357022603955158*f[53]*dv2*dv3); 
  out[39] += volFact*(2.828427124746191*f[22]*wx2*wx3+0.8164965809277261*f[43]*dv2*wx3+0.8164965809277261*f[47]*dv3*wx2+0.2357022603955158*f[59]*dv2*dv3); 
  out[40] += volFact*(2.828427124746191*f[0]*wx3_sq+1.632993161855453*f[6]*dv3*wx3+0.2357022603955158*f[0]*dv3_sq); 
  out[41] += volFact*(2.828427124746191*f[1]*wx3_sq+1.632993161855453*f[17]*dv3*wx3+0.2357022603955158*f[1]*dv3_sq); 
  out[42] += volFact*(2.828427124746191*f[2]*wx3_sq+1.632993161855453*f[18]*dv3*wx3+0.2357022603955158*f[2]*dv3_sq); 
  out[43] += volFact*(2.828427124746191*f[3]*wx3_sq+1.632993161855453*f[19]*dv3*wx3+0.2357022603955158*f[3]*dv3_sq); 
  out[44] += volFact*(2.828427124746191*f[7]*wx3_sq+1.632993161855453*f[32]*dv3*wx3+0.2357022603955158*f[7]*dv3_sq); 
  out[45] += volFact*(2.828427124746191*f[8]*wx3_sq+1.632993161855453*f[33]*dv3*wx3+0.2357022603955158*f[8]*dv3_sq); 
  out[46] += volFact*(2.828427124746191*f[9]*wx3_sq+1.632993161855453*f[34]*dv3*wx3+0.2357022603955158*f[9]*dv3_sq); 
  out[47] += volFact*(2.828427124746191*f[22]*wx3_sq+1.632993161855453*f[47]*dv3*wx3+0.2357022603955158*f[22]*dv3_sq); 
} 
void vlasov_mom_3x3v_m2_p1(const double *w, const double *dxv, const long *idx, const double *f, double* restrict out) 
{ 
  const double volFact = dxv[3]*dxv[4]*dxv[5]/8; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[5], dv3 = dxv[5]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  out[0] += volFact*(2.828427124746191*f[0]*wx3_sq+1.632993161855453*f[6]*dv3*wx3+2.828427124746191*f[0]*wx2_sq+1.632993161855453*f[5]*dv2*wx2+2.828427124746191*f[0]*wx1_sq+1.632993161855453*f[4]*dv1*wx1+0.2357022603955158*f[0]*dv3_sq+0.2357022603955158*f[0]*dv2_sq+0.2357022603955158*f[0]*dv1_sq); 
  out[1] += volFact*(2.828427124746191*f[1]*wx3_sq+1.632993161855453*f[17]*dv3*wx3+2.828427124746191*f[1]*wx2_sq+1.632993161855453*f[13]*dv2*wx2+2.828427124746191*f[1]*wx1_sq+1.632993161855453*f[10]*dv1*wx1+0.2357022603955158*f[1]*dv3_sq+0.2357022603955158*f[1]*dv2_sq+0.2357022603955158*f[1]*dv1_sq); 
  out[2] += volFact*(2.828427124746191*f[2]*wx3_sq+1.632993161855453*f[18]*dv3*wx3+2.828427124746191*f[2]*wx2_sq+1.632993161855453*f[14]*dv2*wx2+2.828427124746191*f[2]*wx1_sq+1.632993161855453*f[11]*dv1*wx1+0.2357022603955158*f[2]*dv3_sq+0.2357022603955158*f[2]*dv2_sq+0.2357022603955158*f[2]*dv1_sq); 
  out[3] += volFact*(2.828427124746191*f[3]*wx3_sq+1.632993161855453*f[19]*dv3*wx3+2.828427124746191*f[3]*wx2_sq+1.632993161855453*f[15]*dv2*wx2+2.828427124746191*f[3]*wx1_sq+1.632993161855453*f[12]*dv1*wx1+0.2357022603955158*f[3]*dv3_sq+0.2357022603955158*f[3]*dv2_sq+0.2357022603955158*f[3]*dv1_sq); 
  out[4] += volFact*(2.828427124746191*f[7]*wx3_sq+1.632993161855453*f[32]*dv3*wx3+2.828427124746191*f[7]*wx2_sq+1.632993161855453*f[26]*dv2*wx2+2.828427124746191*f[7]*wx1_sq+1.632993161855453*f[23]*dv1*wx1+0.2357022603955158*f[7]*dv3_sq+0.2357022603955158*f[7]*dv2_sq+0.2357022603955158*f[7]*dv1_sq); 
  out[5] += volFact*(2.828427124746191*f[8]*wx3_sq+1.632993161855453*f[33]*dv3*wx3+2.828427124746191*f[8]*wx2_sq+1.632993161855453*f[27]*dv2*wx2+2.828427124746191*f[8]*wx1_sq+1.632993161855453*f[24]*dv1*wx1+0.2357022603955158*f[8]*dv3_sq+0.2357022603955158*f[8]*dv2_sq+0.2357022603955158*f[8]*dv1_sq); 
  out[6] += volFact*(2.828427124746191*f[9]*wx3_sq+1.632993161855453*f[34]*dv3*wx3+2.828427124746191*f[9]*wx2_sq+1.632993161855453*f[28]*dv2*wx2+2.828427124746191*f[9]*wx1_sq+1.632993161855453*f[25]*dv1*wx1+0.2357022603955158*f[9]*dv3_sq+0.2357022603955158*f[9]*dv2_sq+0.2357022603955158*f[9]*dv1_sq); 
  out[7] += volFact*(2.828427124746191*f[22]*wx3_sq+1.632993161855453*f[47]*dv3*wx3+2.828427124746191*f[22]*wx2_sq+1.632993161855453*f[43]*dv2*wx2+2.828427124746191*f[22]*wx1_sq+1.632993161855453*f[42]*dv1*wx1+0.2357022603955158*f[22]*dv3_sq+0.2357022603955158*f[22]*dv2_sq+0.2357022603955158*f[22]*dv1_sq); 
} 
void vlasov_mom_3x3v_m3i_p1(const double *w, const double *dxv, const long *idx, const double *f, double* restrict out) 
{ 
  const double volFact = dxv[3]*dxv[4]*dxv[5]/8; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  const double wx3 = w[5], dv3 = dxv[5]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  const double wx3_cu = wx3*wx3*wx3, dv3_cu = dv3*dv3*dv3; 
  out[0] += volFact*(2.828427124746191*f[0]*wx1*wx3_sq+0.8164965809277261*f[4]*dv1*wx3_sq+1.632993161855453*f[6]*dv3*wx1*wx3+0.4714045207910317*f[20]*dv1*dv3*wx3+2.828427124746191*f[0]*wx1*wx2_sq+0.8164965809277261*f[4]*dv1*wx2_sq+1.632993161855453*f[5]*dv2*wx1*wx2+0.4714045207910317*f[16]*dv1*dv2*wx2+2.828427124746191*f[0]*wx1*wx1_sq+2.449489742783178*f[4]*dv1*wx1_sq+0.2357022603955158*f[0]*dv3_sq*wx1+0.2357022603955158*f[0]*dv2_sq*wx1+0.7071067811865475*f[0]*dv1_sq*wx1+0.06804138174397717*f[4]*dv1*dv3_sq+0.06804138174397717*f[4]*dv1*dv2_sq+0.1224744871391589*f[4]*dv1*dv1_sq); 
  out[1] += volFact*(2.828427124746191*f[1]*wx1*wx3_sq+0.8164965809277261*f[10]*dv1*wx3_sq+1.632993161855453*f[17]*dv3*wx1*wx3+0.4714045207910317*f[35]*dv1*dv3*wx3+2.828427124746191*f[1]*wx1*wx2_sq+0.8164965809277261*f[10]*dv1*wx2_sq+1.632993161855453*f[13]*dv2*wx1*wx2+0.4714045207910317*f[29]*dv1*dv2*wx2+2.828427124746191*f[1]*wx1*wx1_sq+2.449489742783178*f[10]*dv1*wx1_sq+0.2357022603955158*f[1]*dv3_sq*wx1+0.2357022603955158*f[1]*dv2_sq*wx1+0.7071067811865475*f[1]*dv1_sq*wx1+0.06804138174397717*f[10]*dv1*dv3_sq+0.06804138174397717*f[10]*dv1*dv2_sq+0.1224744871391589*f[10]*dv1*dv1_sq); 
  out[2] += volFact*(2.828427124746191*f[2]*wx1*wx3_sq+0.8164965809277261*f[11]*dv1*wx3_sq+1.632993161855453*f[18]*dv3*wx1*wx3+0.4714045207910317*f[36]*dv1*dv3*wx3+2.828427124746191*f[2]*wx1*wx2_sq+0.8164965809277261*f[11]*dv1*wx2_sq+1.632993161855453*f[14]*dv2*wx1*wx2+0.4714045207910317*f[30]*dv1*dv2*wx2+2.828427124746191*f[2]*wx1*wx1_sq+2.449489742783178*f[11]*dv1*wx1_sq+0.2357022603955158*f[2]*dv3_sq*wx1+0.2357022603955158*f[2]*dv2_sq*wx1+0.7071067811865475*f[2]*dv1_sq*wx1+0.06804138174397717*f[11]*dv1*dv3_sq+0.06804138174397717*f[11]*dv1*dv2_sq+0.1224744871391589*f[11]*dv1*dv1_sq); 
  out[3] += volFact*(2.828427124746191*f[3]*wx1*wx3_sq+0.8164965809277261*f[12]*dv1*wx3_sq+1.632993161855453*f[19]*dv3*wx1*wx3+0.4714045207910317*f[37]*dv1*dv3*wx3+2.828427124746191*f[3]*wx1*wx2_sq+0.8164965809277261*f[12]*dv1*wx2_sq+1.632993161855453*f[15]*dv2*wx1*wx2+0.4714045207910317*f[31]*dv1*dv2*wx2+2.828427124746191*f[3]*wx1*wx1_sq+2.449489742783178*f[12]*dv1*wx1_sq+0.2357022603955158*f[3]*dv3_sq*wx1+0.2357022603955158*f[3]*dv2_sq*wx1+0.7071067811865475*f[3]*dv1_sq*wx1+0.06804138174397717*f[12]*dv1*dv3_sq+0.06804138174397717*f[12]*dv1*dv2_sq+0.1224744871391589*f[12]*dv1*dv1_sq); 
  out[4] += volFact*(2.828427124746191*f[7]*wx1*wx3_sq+0.8164965809277261*f[23]*dv1*wx3_sq+1.632993161855453*f[32]*dv3*wx1*wx3+0.4714045207910317*f[48]*dv1*dv3*wx3+2.828427124746191*f[7]*wx1*wx2_sq+0.8164965809277261*f[23]*dv1*wx2_sq+1.632993161855453*f[26]*dv2*wx1*wx2+0.4714045207910317*f[44]*dv1*dv2*wx2+2.828427124746191*f[7]*wx1*wx1_sq+2.449489742783178*f[23]*dv1*wx1_sq+0.2357022603955158*f[7]*dv3_sq*wx1+0.2357022603955158*f[7]*dv2_sq*wx1+0.7071067811865475*f[7]*dv1_sq*wx1+0.06804138174397717*f[23]*dv1*dv3_sq+0.06804138174397717*f[23]*dv1*dv2_sq+0.1224744871391589*f[23]*dv1*dv1_sq); 
  out[5] += volFact*(2.828427124746191*f[8]*wx1*wx3_sq+0.8164965809277261*f[24]*dv1*wx3_sq+1.632993161855453*f[33]*dv3*wx1*wx3+0.4714045207910317*f[49]*dv1*dv3*wx3+2.828427124746191*f[8]*wx1*wx2_sq+0.8164965809277261*f[24]*dv1*wx2_sq+1.632993161855453*f[27]*dv2*wx1*wx2+0.4714045207910317*f[45]*dv1*dv2*wx2+2.828427124746191*f[8]*wx1*wx1_sq+2.449489742783178*f[24]*dv1*wx1_sq+0.2357022603955158*f[8]*dv3_sq*wx1+0.2357022603955158*f[8]*dv2_sq*wx1+0.7071067811865475*f[8]*dv1_sq*wx1+0.06804138174397717*f[24]*dv1*dv3_sq+0.06804138174397717*f[24]*dv1*dv2_sq+0.1224744871391589*f[24]*dv1*dv1_sq); 
  out[6] += volFact*(2.828427124746191*f[9]*wx1*wx3_sq+0.8164965809277261*f[25]*dv1*wx3_sq+1.632993161855453*f[34]*dv3*wx1*wx3+0.4714045207910317*f[50]*dv1*dv3*wx3+2.828427124746191*f[9]*wx1*wx2_sq+0.8164965809277261*f[25]*dv1*wx2_sq+1.632993161855453*f[28]*dv2*wx1*wx2+0.4714045207910317*f[46]*dv1*dv2*wx2+2.828427124746191*f[9]*wx1*wx1_sq+2.449489742783178*f[25]*dv1*wx1_sq+0.2357022603955158*f[9]*dv3_sq*wx1+0.2357022603955158*f[9]*dv2_sq*wx1+0.7071067811865475*f[9]*dv1_sq*wx1+0.06804138174397717*f[25]*dv1*dv3_sq+0.06804138174397717*f[25]*dv1*dv2_sq+0.1224744871391589*f[25]*dv1*dv1_sq); 
  out[7] += volFact*(2.828427124746191*f[22]*wx1*wx3_sq+0.8164965809277261*f[42]*dv1*wx3_sq+1.632993161855453*f[47]*dv3*wx1*wx3+0.4714045207910317*f[58]*dv1*dv3*wx3+2.828427124746191*f[22]*wx1*wx2_sq+0.8164965809277261*f[42]*dv1*wx2_sq+1.632993161855453*f[43]*dv2*wx1*wx2+0.4714045207910317*f[57]*dv1*dv2*wx2+2.828427124746191*f[22]*wx1*wx1_sq+2.449489742783178*f[42]*dv1*wx1_sq+0.2357022603955158*f[22]*dv3_sq*wx1+0.2357022603955158*f[22]*dv2_sq*wx1+0.7071067811865475*f[22]*dv1_sq*wx1+0.06804138174397717*f[42]*dv1*dv3_sq+0.06804138174397717*f[42]*dv1*dv2_sq+0.1224744871391589*f[42]*dv1*dv1_sq); 
  out[8] += volFact*(2.828427124746191*f[0]*wx2*wx3_sq+0.8164965809277261*f[5]*dv2*wx3_sq+1.632993161855453*f[6]*dv3*wx2*wx3+0.4714045207910317*f[21]*dv2*dv3*wx3+2.828427124746191*f[0]*wx2*wx2_sq+2.449489742783178*f[5]*dv2*wx2_sq+2.828427124746191*f[0]*wx1_sq*wx2+1.632993161855453*f[4]*dv1*wx1*wx2+0.2357022603955158*f[0]*dv3_sq*wx2+0.7071067811865475*f[0]*dv2_sq*wx2+0.2357022603955158*f[0]*dv1_sq*wx2+0.8164965809277261*f[5]*dv2*wx1_sq+0.4714045207910317*f[16]*dv1*dv2*wx1+0.06804138174397717*f[5]*dv2*dv3_sq+0.1224744871391589*f[5]*dv2*dv2_sq+0.06804138174397717*f[5]*dv1_sq*dv2); 
  out[9] += volFact*(2.828427124746191*f[1]*wx2*wx3_sq+0.8164965809277261*f[13]*dv2*wx3_sq+1.632993161855453*f[17]*dv3*wx2*wx3+0.4714045207910317*f[38]*dv2*dv3*wx3+2.828427124746191*f[1]*wx2*wx2_sq+2.449489742783178*f[13]*dv2*wx2_sq+2.828427124746191*f[1]*wx1_sq*wx2+1.632993161855453*f[10]*dv1*wx1*wx2+0.2357022603955158*f[1]*dv3_sq*wx2+0.7071067811865475*f[1]*dv2_sq*wx2+0.2357022603955158*f[1]*dv1_sq*wx2+0.8164965809277261*f[13]*dv2*wx1_sq+0.4714045207910317*f[29]*dv1*dv2*wx1+0.06804138174397717*f[13]*dv2*dv3_sq+0.1224744871391589*f[13]*dv2*dv2_sq+0.06804138174397717*f[13]*dv1_sq*dv2); 
  out[10] += volFact*(2.828427124746191*f[2]*wx2*wx3_sq+0.8164965809277261*f[14]*dv2*wx3_sq+1.632993161855453*f[18]*dv3*wx2*wx3+0.4714045207910317*f[39]*dv2*dv3*wx3+2.828427124746191*f[2]*wx2*wx2_sq+2.449489742783178*f[14]*dv2*wx2_sq+2.828427124746191*f[2]*wx1_sq*wx2+1.632993161855453*f[11]*dv1*wx1*wx2+0.2357022603955158*f[2]*dv3_sq*wx2+0.7071067811865475*f[2]*dv2_sq*wx2+0.2357022603955158*f[2]*dv1_sq*wx2+0.8164965809277261*f[14]*dv2*wx1_sq+0.4714045207910317*f[30]*dv1*dv2*wx1+0.06804138174397717*f[14]*dv2*dv3_sq+0.1224744871391589*f[14]*dv2*dv2_sq+0.06804138174397717*f[14]*dv1_sq*dv2); 
  out[11] += volFact*(2.828427124746191*f[3]*wx2*wx3_sq+0.8164965809277261*f[15]*dv2*wx3_sq+1.632993161855453*f[19]*dv3*wx2*wx3+0.4714045207910317*f[40]*dv2*dv3*wx3+2.828427124746191*f[3]*wx2*wx2_sq+2.449489742783178*f[15]*dv2*wx2_sq+2.828427124746191*f[3]*wx1_sq*wx2+1.632993161855453*f[12]*dv1*wx1*wx2+0.2357022603955158*f[3]*dv3_sq*wx2+0.7071067811865475*f[3]*dv2_sq*wx2+0.2357022603955158*f[3]*dv1_sq*wx2+0.8164965809277261*f[15]*dv2*wx1_sq+0.4714045207910317*f[31]*dv1*dv2*wx1+0.06804138174397717*f[15]*dv2*dv3_sq+0.1224744871391589*f[15]*dv2*dv2_sq+0.06804138174397717*f[15]*dv1_sq*dv2); 
  out[12] += volFact*(2.828427124746191*f[7]*wx2*wx3_sq+0.8164965809277261*f[26]*dv2*wx3_sq+1.632993161855453*f[32]*dv3*wx2*wx3+0.4714045207910317*f[51]*dv2*dv3*wx3+2.828427124746191*f[7]*wx2*wx2_sq+2.449489742783178*f[26]*dv2*wx2_sq+2.828427124746191*f[7]*wx1_sq*wx2+1.632993161855453*f[23]*dv1*wx1*wx2+0.2357022603955158*f[7]*dv3_sq*wx2+0.7071067811865475*f[7]*dv2_sq*wx2+0.2357022603955158*f[7]*dv1_sq*wx2+0.8164965809277261*f[26]*dv2*wx1_sq+0.4714045207910317*f[44]*dv1*dv2*wx1+0.06804138174397717*f[26]*dv2*dv3_sq+0.1224744871391589*f[26]*dv2*dv2_sq+0.06804138174397717*f[26]*dv1_sq*dv2); 
  out[13] += volFact*(2.828427124746191*f[8]*wx2*wx3_sq+0.8164965809277261*f[27]*dv2*wx3_sq+1.632993161855453*f[33]*dv3*wx2*wx3+0.4714045207910317*f[52]*dv2*dv3*wx3+2.828427124746191*f[8]*wx2*wx2_sq+2.449489742783178*f[27]*dv2*wx2_sq+2.828427124746191*f[8]*wx1_sq*wx2+1.632993161855453*f[24]*dv1*wx1*wx2+0.2357022603955158*f[8]*dv3_sq*wx2+0.7071067811865475*f[8]*dv2_sq*wx2+0.2357022603955158*f[8]*dv1_sq*wx2+0.8164965809277261*f[27]*dv2*wx1_sq+0.4714045207910317*f[45]*dv1*dv2*wx1+0.06804138174397717*f[27]*dv2*dv3_sq+0.1224744871391589*f[27]*dv2*dv2_sq+0.06804138174397717*f[27]*dv1_sq*dv2); 
  out[14] += volFact*(2.828427124746191*f[9]*wx2*wx3_sq+0.8164965809277261*f[28]*dv2*wx3_sq+1.632993161855453*f[34]*dv3*wx2*wx3+0.4714045207910317*f[53]*dv2*dv3*wx3+2.828427124746191*f[9]*wx2*wx2_sq+2.449489742783178*f[28]*dv2*wx2_sq+2.828427124746191*f[9]*wx1_sq*wx2+1.632993161855453*f[25]*dv1*wx1*wx2+0.2357022603955158*f[9]*dv3_sq*wx2+0.7071067811865475*f[9]*dv2_sq*wx2+0.2357022603955158*f[9]*dv1_sq*wx2+0.8164965809277261*f[28]*dv2*wx1_sq+0.4714045207910317*f[46]*dv1*dv2*wx1+0.06804138174397717*f[28]*dv2*dv3_sq+0.1224744871391589*f[28]*dv2*dv2_sq+0.06804138174397717*f[28]*dv1_sq*dv2); 
  out[15] += volFact*(2.828427124746191*f[22]*wx2*wx3_sq+0.8164965809277261*f[43]*dv2*wx3_sq+1.632993161855453*f[47]*dv3*wx2*wx3+0.4714045207910317*f[59]*dv2*dv3*wx3+2.828427124746191*f[22]*wx2*wx2_sq+2.449489742783178*f[43]*dv2*wx2_sq+2.828427124746191*f[22]*wx1_sq*wx2+1.632993161855453*f[42]*dv1*wx1*wx2+0.2357022603955158*f[22]*dv3_sq*wx2+0.7071067811865475*f[22]*dv2_sq*wx2+0.2357022603955158*f[22]*dv1_sq*wx2+0.8164965809277261*f[43]*dv2*wx1_sq+0.4714045207910317*f[57]*dv1*dv2*wx1+0.06804138174397717*f[43]*dv2*dv3_sq+0.1224744871391589*f[43]*dv2*dv2_sq+0.06804138174397717*f[43]*dv1_sq*dv2); 
  out[16] += volFact*(2.828427124746191*f[0]*wx3*wx3_sq+2.449489742783178*f[6]*dv3*wx3_sq+2.828427124746191*f[0]*wx2_sq*wx3+1.632993161855453*f[5]*dv2*wx2*wx3+2.828427124746191*f[0]*wx1_sq*wx3+1.632993161855453*f[4]*dv1*wx1*wx3+0.7071067811865475*f[0]*dv3_sq*wx3+0.2357022603955158*f[0]*dv2_sq*wx3+0.2357022603955158*f[0]*dv1_sq*wx3+0.8164965809277261*f[6]*dv3*wx2_sq+0.4714045207910317*f[21]*dv2*dv3*wx2+0.8164965809277261*f[6]*dv3*wx1_sq+0.4714045207910317*f[20]*dv1*dv3*wx1+0.1224744871391589*f[6]*dv3*dv3_sq+0.06804138174397717*f[6]*dv2_sq*dv3+0.06804138174397717*f[6]*dv1_sq*dv3); 
  out[17] += volFact*(2.828427124746191*f[1]*wx3*wx3_sq+2.449489742783178*f[17]*dv3*wx3_sq+2.828427124746191*f[1]*wx2_sq*wx3+1.632993161855453*f[13]*dv2*wx2*wx3+2.828427124746191*f[1]*wx1_sq*wx3+1.632993161855453*f[10]*dv1*wx1*wx3+0.7071067811865475*f[1]*dv3_sq*wx3+0.2357022603955158*f[1]*dv2_sq*wx3+0.2357022603955158*f[1]*dv1_sq*wx3+0.8164965809277261*f[17]*dv3*wx2_sq+0.4714045207910317*f[38]*dv2*dv3*wx2+0.8164965809277261*f[17]*dv3*wx1_sq+0.4714045207910317*f[35]*dv1*dv3*wx1+0.1224744871391589*f[17]*dv3*dv3_sq+0.06804138174397717*f[17]*dv2_sq*dv3+0.06804138174397717*f[17]*dv1_sq*dv3); 
  out[18] += volFact*(2.828427124746191*f[2]*wx3*wx3_sq+2.449489742783178*f[18]*dv3*wx3_sq+2.828427124746191*f[2]*wx2_sq*wx3+1.632993161855453*f[14]*dv2*wx2*wx3+2.828427124746191*f[2]*wx1_sq*wx3+1.632993161855453*f[11]*dv1*wx1*wx3+0.7071067811865475*f[2]*dv3_sq*wx3+0.2357022603955158*f[2]*dv2_sq*wx3+0.2357022603955158*f[2]*dv1_sq*wx3+0.8164965809277261*f[18]*dv3*wx2_sq+0.4714045207910317*f[39]*dv2*dv3*wx2+0.8164965809277261*f[18]*dv3*wx1_sq+0.4714045207910317*f[36]*dv1*dv3*wx1+0.1224744871391589*f[18]*dv3*dv3_sq+0.06804138174397717*f[18]*dv2_sq*dv3+0.06804138174397717*f[18]*dv1_sq*dv3); 
  out[19] += volFact*(2.828427124746191*f[3]*wx3*wx3_sq+2.449489742783178*f[19]*dv3*wx3_sq+2.828427124746191*f[3]*wx2_sq*wx3+1.632993161855453*f[15]*dv2*wx2*wx3+2.828427124746191*f[3]*wx1_sq*wx3+1.632993161855453*f[12]*dv1*wx1*wx3+0.7071067811865475*f[3]*dv3_sq*wx3+0.2357022603955158*f[3]*dv2_sq*wx3+0.2357022603955158*f[3]*dv1_sq*wx3+0.8164965809277261*f[19]*dv3*wx2_sq+0.4714045207910317*f[40]*dv2*dv3*wx2+0.8164965809277261*f[19]*dv3*wx1_sq+0.4714045207910317*f[37]*dv1*dv3*wx1+0.1224744871391589*f[19]*dv3*dv3_sq+0.06804138174397717*f[19]*dv2_sq*dv3+0.06804138174397717*f[19]*dv1_sq*dv3); 
  out[20] += volFact*(2.828427124746191*f[7]*wx3*wx3_sq+2.449489742783178*f[32]*dv3*wx3_sq+2.828427124746191*f[7]*wx2_sq*wx3+1.632993161855453*f[26]*dv2*wx2*wx3+2.828427124746191*f[7]*wx1_sq*wx3+1.632993161855453*f[23]*dv1*wx1*wx3+0.7071067811865475*f[7]*dv3_sq*wx3+0.2357022603955158*f[7]*dv2_sq*wx3+0.2357022603955158*f[7]*dv1_sq*wx3+0.8164965809277261*f[32]*dv3*wx2_sq+0.4714045207910317*f[51]*dv2*dv3*wx2+0.8164965809277261*f[32]*dv3*wx1_sq+0.4714045207910317*f[48]*dv1*dv3*wx1+0.1224744871391589*f[32]*dv3*dv3_sq+0.06804138174397717*f[32]*dv2_sq*dv3+0.06804138174397717*f[32]*dv1_sq*dv3); 
  out[21] += volFact*(2.828427124746191*f[8]*wx3*wx3_sq+2.449489742783178*f[33]*dv3*wx3_sq+2.828427124746191*f[8]*wx2_sq*wx3+1.632993161855453*f[27]*dv2*wx2*wx3+2.828427124746191*f[8]*wx1_sq*wx3+1.632993161855453*f[24]*dv1*wx1*wx3+0.7071067811865475*f[8]*dv3_sq*wx3+0.2357022603955158*f[8]*dv2_sq*wx3+0.2357022603955158*f[8]*dv1_sq*wx3+0.8164965809277261*f[33]*dv3*wx2_sq+0.4714045207910317*f[52]*dv2*dv3*wx2+0.8164965809277261*f[33]*dv3*wx1_sq+0.4714045207910317*f[49]*dv1*dv3*wx1+0.1224744871391589*f[33]*dv3*dv3_sq+0.06804138174397717*f[33]*dv2_sq*dv3+0.06804138174397717*f[33]*dv1_sq*dv3); 
  out[22] += volFact*(2.828427124746191*f[9]*wx3*wx3_sq+2.449489742783178*f[34]*dv3*wx3_sq+2.828427124746191*f[9]*wx2_sq*wx3+1.632993161855453*f[28]*dv2*wx2*wx3+2.828427124746191*f[9]*wx1_sq*wx3+1.632993161855453*f[25]*dv1*wx1*wx3+0.7071067811865475*f[9]*dv3_sq*wx3+0.2357022603955158*f[9]*dv2_sq*wx3+0.2357022603955158*f[9]*dv1_sq*wx3+0.8164965809277261*f[34]*dv3*wx2_sq+0.4714045207910317*f[53]*dv2*dv3*wx2+0.8164965809277261*f[34]*dv3*wx1_sq+0.4714045207910317*f[50]*dv1*dv3*wx1+0.1224744871391589*f[34]*dv3*dv3_sq+0.06804138174397717*f[34]*dv2_sq*dv3+0.06804138174397717*f[34]*dv1_sq*dv3); 
  out[23] += volFact*(2.828427124746191*f[22]*wx3*wx3_sq+2.449489742783178*f[47]*dv3*wx3_sq+2.828427124746191*f[22]*wx2_sq*wx3+1.632993161855453*f[43]*dv2*wx2*wx3+2.828427124746191*f[22]*wx1_sq*wx3+1.632993161855453*f[42]*dv1*wx1*wx3+0.7071067811865475*f[22]*dv3_sq*wx3+0.2357022603955158*f[22]*dv2_sq*wx3+0.2357022603955158*f[22]*dv1_sq*wx3+0.8164965809277261*f[47]*dv3*wx2_sq+0.4714045207910317*f[59]*dv2*dv3*wx2+0.8164965809277261*f[47]*dv3*wx1_sq+0.4714045207910317*f[58]*dv1*dv3*wx1+0.1224744871391589*f[47]*dv3*dv3_sq+0.06804138174397717*f[47]*dv2_sq*dv3+0.06804138174397717*f[47]*dv1_sq*dv3); 
} 
