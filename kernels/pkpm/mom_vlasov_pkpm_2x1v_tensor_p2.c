#include <gkyl_mom_vlasov_pkpm_kernels.h> 
GKYL_CU_DH void mom_vlasov_pkpm_2x1v_tensor_p2(const double *w, const double *dxv, const int *idx, double mass, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]/2.0; 
  const double wvpar = w[2], dvpar = dxv[2]; 
  const double wvpar_sq = wvpar*wvpar, dvpar_sq = dvpar*dvpar; 

  const double *F_0 = &f[0]; 
  const double *G_1 = &f[27]; 
  out[0] += 1.414213562373095*F_0[0]*mass*volFact; 
  out[1] += 1.414213562373095*F_0[1]*mass*volFact; 
  out[2] += 1.414213562373095*F_0[2]*mass*volFact; 
  out[3] += 1.414213562373095*F_0[4]*mass*volFact; 
  out[4] += 1.414213562373095*F_0[7]*mass*volFact; 
  out[5] += 1.414213562373095*F_0[8]*mass*volFact; 
  out[6] += 1.414213562373095*F_0[11]*mass*volFact; 
  out[7] += 1.414213562373095*F_0[12]*mass*volFact; 
  out[8] += 1.414213562373095*F_0[20]*mass*volFact; 
  out[9] += mass*volFact*(1.414213562373095*F_0[0]*wvpar_sq+0.8164965809277261*F_0[3]*dvpar*wvpar+0.105409255338946*F_0[9]*dvpar_sq+0.1178511301977579*F_0[0]*dvpar_sq); 
  out[10] += mass*volFact*(1.414213562373095*F_0[1]*wvpar_sq+0.8164965809277261*F_0[5]*dvpar*wvpar+0.105409255338946*F_0[15]*dvpar_sq+0.1178511301977579*F_0[1]*dvpar_sq); 
  out[11] += mass*volFact*(1.414213562373095*F_0[2]*wvpar_sq+0.8164965809277261*F_0[6]*dvpar*wvpar+0.105409255338946*F_0[16]*dvpar_sq+0.1178511301977579*F_0[2]*dvpar_sq); 
  out[12] += mass*volFact*(1.414213562373095*F_0[4]*wvpar_sq+0.8164965809277261*F_0[10]*dvpar*wvpar+0.105409255338946*F_0[19]*dvpar_sq+0.1178511301977579*F_0[4]*dvpar_sq); 
  out[13] += mass*volFact*(1.414213562373095*F_0[7]*wvpar_sq+0.816496580927726*F_0[13]*dvpar*wvpar+0.105409255338946*F_0[21]*dvpar_sq+0.1178511301977579*F_0[7]*dvpar_sq); 
  out[14] += mass*volFact*(1.414213562373095*F_0[8]*wvpar_sq+0.816496580927726*F_0[14]*dvpar*wvpar+0.105409255338946*F_0[22]*dvpar_sq+0.1178511301977579*F_0[8]*dvpar_sq); 
  out[15] += mass*volFact*(1.414213562373095*F_0[11]*wvpar_sq+0.816496580927726*F_0[17]*dvpar*wvpar+0.105409255338946*F_0[24]*dvpar_sq+0.1178511301977579*F_0[11]*dvpar_sq); 
  out[16] += mass*volFact*(1.414213562373095*F_0[12]*wvpar_sq+0.816496580927726*F_0[18]*dvpar*wvpar+0.105409255338946*F_0[25]*dvpar_sq+0.1178511301977579*F_0[12]*dvpar_sq); 
  out[17] += mass*volFact*(1.414213562373095*F_0[20]*wvpar_sq+0.8164965809277261*F_0[23]*dvpar*wvpar+0.105409255338946*F_0[26]*dvpar_sq+0.1178511301977579*F_0[20]*dvpar_sq); 
  out[18] += 1.414213562373095*G_1[0]*mass*volFact; 
  out[19] += 1.414213562373095*G_1[1]*mass*volFact; 
  out[20] += 1.414213562373095*G_1[2]*mass*volFact; 
  out[21] += 1.414213562373095*G_1[4]*mass*volFact; 
  out[22] += 1.414213562373095*G_1[7]*mass*volFact; 
  out[23] += 1.414213562373095*G_1[8]*mass*volFact; 
  out[24] += 1.414213562373095*G_1[11]*mass*volFact; 
  out[25] += 1.414213562373095*G_1[12]*mass*volFact; 
  out[26] += 1.414213562373095*G_1[20]*mass*volFact; 
} 
GKYL_CU_DH void mom_vlasov_pkpm_diag_2x1v_tensor_p2(const double *w, const double *dxv, const int *idx, double mass, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]/2.0; 
  const double wvpar = w[2], dvpar = dxv[2]; 
  const double wvpar_sq = wvpar*wvpar, dvpar_sq = dvpar*dvpar; 
  const double wvpar_cu = wvpar*wvpar*wvpar, dvpar_cu = dvpar*dvpar*dvpar; 
  const double wvpar_qu = wvpar*wvpar*wvpar*wvpar, dvpar_qu = dvpar*dvpar*dvpar*dvpar; 

  const double *F_0 = &f[0]; 
  const double *G_1 = &f[27]; 
  out[0] += 1.414213562373095*F_0[0]*mass*volFact; 
  out[1] += 1.414213562373095*F_0[1]*mass*volFact; 
  out[2] += 1.414213562373095*F_0[2]*mass*volFact; 
  out[3] += 1.414213562373095*F_0[4]*mass*volFact; 
  out[4] += 1.414213562373095*F_0[7]*mass*volFact; 
  out[5] += 1.414213562373095*F_0[8]*mass*volFact; 
  out[6] += 1.414213562373095*F_0[11]*mass*volFact; 
  out[7] += 1.414213562373095*F_0[12]*mass*volFact; 
  out[8] += 1.414213562373095*F_0[20]*mass*volFact; 
  out[9] += mass*volFact*(1.414213562373095*F_0[0]*wvpar+0.408248290463863*F_0[3]*dvpar); 
  out[10] += mass*volFact*(1.414213562373095*F_0[1]*wvpar+0.408248290463863*F_0[5]*dvpar); 
  out[11] += mass*volFact*(1.414213562373095*F_0[2]*wvpar+0.408248290463863*F_0[6]*dvpar); 
  out[12] += mass*volFact*(1.414213562373095*F_0[4]*wvpar+0.408248290463863*F_0[10]*dvpar); 
  out[13] += mass*volFact*(1.414213562373095*F_0[7]*wvpar+0.408248290463863*F_0[13]*dvpar); 
  out[14] += mass*volFact*(1.414213562373095*F_0[8]*wvpar+0.408248290463863*F_0[14]*dvpar); 
  out[15] += mass*volFact*(1.414213562373095*F_0[11]*wvpar+0.408248290463863*F_0[17]*dvpar); 
  out[16] += mass*volFact*(1.414213562373095*F_0[12]*wvpar+0.408248290463863*F_0[18]*dvpar); 
  out[17] += mass*volFact*(1.414213562373095*F_0[20]*wvpar+0.408248290463863*F_0[23]*dvpar); 
  out[18] += mass*volFact*(1.414213562373095*F_0[0]*wvpar_sq+0.8164965809277261*F_0[3]*dvpar*wvpar+0.105409255338946*F_0[9]*dvpar_sq+0.1178511301977579*F_0[0]*dvpar_sq); 
  out[19] += mass*volFact*(1.414213562373095*F_0[1]*wvpar_sq+0.8164965809277261*F_0[5]*dvpar*wvpar+0.105409255338946*F_0[15]*dvpar_sq+0.1178511301977579*F_0[1]*dvpar_sq); 
  out[20] += mass*volFact*(1.414213562373095*F_0[2]*wvpar_sq+0.8164965809277261*F_0[6]*dvpar*wvpar+0.105409255338946*F_0[16]*dvpar_sq+0.1178511301977579*F_0[2]*dvpar_sq); 
  out[21] += mass*volFact*(1.414213562373095*F_0[4]*wvpar_sq+0.8164965809277261*F_0[10]*dvpar*wvpar+0.105409255338946*F_0[19]*dvpar_sq+0.1178511301977579*F_0[4]*dvpar_sq); 
  out[22] += mass*volFact*(1.414213562373095*F_0[7]*wvpar_sq+0.816496580927726*F_0[13]*dvpar*wvpar+0.105409255338946*F_0[21]*dvpar_sq+0.1178511301977579*F_0[7]*dvpar_sq); 
  out[23] += mass*volFact*(1.414213562373095*F_0[8]*wvpar_sq+0.816496580927726*F_0[14]*dvpar*wvpar+0.105409255338946*F_0[22]*dvpar_sq+0.1178511301977579*F_0[8]*dvpar_sq); 
  out[24] += mass*volFact*(1.414213562373095*F_0[11]*wvpar_sq+0.816496580927726*F_0[17]*dvpar*wvpar+0.105409255338946*F_0[24]*dvpar_sq+0.1178511301977579*F_0[11]*dvpar_sq); 
  out[25] += mass*volFact*(1.414213562373095*F_0[12]*wvpar_sq+0.816496580927726*F_0[18]*dvpar*wvpar+0.105409255338946*F_0[25]*dvpar_sq+0.1178511301977579*F_0[12]*dvpar_sq); 
  out[26] += mass*volFact*(1.414213562373095*F_0[20]*wvpar_sq+0.8164965809277261*F_0[23]*dvpar*wvpar+0.105409255338946*F_0[26]*dvpar_sq+0.1178511301977579*F_0[20]*dvpar_sq); 
  out[27] += 1.414213562373095*G_1[0]*mass*volFact; 
  out[28] += 1.414213562373095*G_1[1]*mass*volFact; 
  out[29] += 1.414213562373095*G_1[2]*mass*volFact; 
  out[30] += 1.414213562373095*G_1[4]*mass*volFact; 
  out[31] += 1.414213562373095*G_1[7]*mass*volFact; 
  out[32] += 1.414213562373095*G_1[8]*mass*volFact; 
  out[33] += 1.414213562373095*G_1[11]*mass*volFact; 
  out[34] += 1.414213562373095*G_1[12]*mass*volFact; 
  out[35] += 1.414213562373095*G_1[20]*mass*volFact; 
  out[36] += mass*volFact*(1.224744871391589*F_0[3]*dvpar*wvpar_sq+1.414213562373095*F_0[0]*wvpar_cu+0.3162277660168379*F_0[9]*dvpar_sq*wvpar+0.3535533905932737*F_0[0]*dvpar_sq*wvpar+0.06123724356957942*F_0[3]*dvpar_cu); 
  out[37] += mass*volFact*(1.224744871391589*F_0[5]*dvpar*wvpar_sq+1.414213562373095*F_0[1]*wvpar_cu+0.3162277660168379*F_0[15]*dvpar_sq*wvpar+0.3535533905932737*F_0[1]*dvpar_sq*wvpar+0.06123724356957942*F_0[5]*dvpar_cu); 
  out[38] += mass*volFact*(1.224744871391589*F_0[6]*dvpar*wvpar_sq+1.414213562373095*F_0[2]*wvpar_cu+0.3162277660168379*F_0[16]*dvpar_sq*wvpar+0.3535533905932737*F_0[2]*dvpar_sq*wvpar+0.06123724356957942*F_0[6]*dvpar_cu); 
  out[39] += mass*volFact*(1.224744871391589*F_0[10]*dvpar*wvpar_sq+1.414213562373095*F_0[4]*wvpar_cu+0.3162277660168379*F_0[19]*dvpar_sq*wvpar+0.3535533905932737*F_0[4]*dvpar_sq*wvpar+0.06123724356957942*F_0[10]*dvpar_cu); 
  out[40] += mass*volFact*(1.224744871391589*F_0[13]*dvpar*wvpar_sq+1.414213562373095*F_0[7]*wvpar_cu+0.3162277660168379*F_0[21]*dvpar_sq*wvpar+0.3535533905932737*F_0[7]*dvpar_sq*wvpar+0.06123724356957942*F_0[13]*dvpar_cu); 
  out[41] += mass*volFact*(1.224744871391589*F_0[14]*dvpar*wvpar_sq+1.414213562373095*F_0[8]*wvpar_cu+0.3162277660168379*F_0[22]*dvpar_sq*wvpar+0.3535533905932737*F_0[8]*dvpar_sq*wvpar+0.06123724356957942*F_0[14]*dvpar_cu); 
  out[42] += mass*volFact*(1.224744871391589*F_0[17]*dvpar*wvpar_sq+1.414213562373095*F_0[11]*wvpar_cu+0.3162277660168379*F_0[24]*dvpar_sq*wvpar+0.3535533905932737*F_0[11]*dvpar_sq*wvpar+0.06123724356957942*F_0[17]*dvpar_cu); 
  out[43] += mass*volFact*(1.224744871391589*F_0[18]*dvpar*wvpar_sq+1.414213562373095*F_0[12]*wvpar_cu+0.3162277660168379*F_0[25]*dvpar_sq*wvpar+0.3535533905932737*F_0[12]*dvpar_sq*wvpar+0.06123724356957942*F_0[18]*dvpar_cu); 
  out[44] += mass*volFact*(1.224744871391589*F_0[23]*dvpar*wvpar_sq+1.414213562373095*F_0[20]*wvpar_cu+0.3162277660168379*F_0[26]*dvpar_sq*wvpar+0.3535533905932737*F_0[20]*dvpar_sq*wvpar+0.06123724356957942*F_0[23]*dvpar_cu); 
  out[45] += mass*volFact*(1.414213562373095*G_1[0]*wvpar+0.408248290463863*G_1[3]*dvpar); 
  out[46] += mass*volFact*(1.414213562373095*G_1[1]*wvpar+0.408248290463863*G_1[5]*dvpar); 
  out[47] += mass*volFact*(1.414213562373095*G_1[2]*wvpar+0.408248290463863*G_1[6]*dvpar); 
  out[48] += mass*volFact*(1.414213562373095*G_1[4]*wvpar+0.408248290463863*G_1[10]*dvpar); 
  out[49] += mass*volFact*(1.414213562373095*G_1[7]*wvpar+0.408248290463863*G_1[13]*dvpar); 
  out[50] += mass*volFact*(1.414213562373095*G_1[8]*wvpar+0.408248290463863*G_1[14]*dvpar); 
  out[51] += mass*volFact*(1.414213562373095*G_1[11]*wvpar+0.408248290463863*G_1[17]*dvpar); 
  out[52] += mass*volFact*(1.414213562373095*G_1[12]*wvpar+0.408248290463863*G_1[18]*dvpar); 
  out[53] += mass*volFact*(1.414213562373095*G_1[20]*wvpar+0.408248290463863*G_1[23]*dvpar); 
  out[54] += mass*volFact*(0.6324555320336759*F_0[9]*dvpar_sq*wvpar_sq+0.7071067811865475*F_0[0]*dvpar_sq*wvpar_sq+1.414213562373095*F_0[0]*wvpar_qu+1.632993161855453*F_0[3]*dvpar*wvpar_cu+0.2449489742783178*F_0[3]*dvpar_cu*wvpar+0.02258769757263127*F_0[9]*dvpar_qu+0.01767766952966368*F_0[0]*dvpar_qu); 
  out[55] += mass*volFact*(0.632455532033676*F_0[15]*dvpar_sq*wvpar_sq+0.7071067811865475*F_0[1]*dvpar_sq*wvpar_sq+1.414213562373095*F_0[1]*wvpar_qu+1.632993161855453*F_0[5]*dvpar*wvpar_cu+0.2449489742783178*F_0[5]*dvpar_cu*wvpar+0.02258769757263128*F_0[15]*dvpar_qu+0.01767766952966368*F_0[1]*dvpar_qu); 
  out[56] += mass*volFact*(0.632455532033676*F_0[16]*dvpar_sq*wvpar_sq+0.7071067811865475*F_0[2]*dvpar_sq*wvpar_sq+1.414213562373095*F_0[2]*wvpar_qu+1.632993161855453*F_0[6]*dvpar*wvpar_cu+0.2449489742783178*F_0[6]*dvpar_cu*wvpar+0.02258769757263128*F_0[16]*dvpar_qu+0.01767766952966368*F_0[2]*dvpar_qu); 
  out[57] += mass*volFact*(0.6324555320336759*F_0[19]*dvpar_sq*wvpar_sq+0.7071067811865475*F_0[4]*dvpar_sq*wvpar_sq+1.414213562373095*F_0[4]*wvpar_qu+1.632993161855453*F_0[10]*dvpar*wvpar_cu+0.2449489742783178*F_0[10]*dvpar_cu*wvpar+0.02258769757263127*F_0[19]*dvpar_qu+0.01767766952966368*F_0[4]*dvpar_qu); 
  out[58] += mass*volFact*(0.6324555320336759*F_0[21]*dvpar_sq*wvpar_sq+0.7071067811865475*F_0[7]*dvpar_sq*wvpar_sq+1.414213562373095*F_0[7]*wvpar_qu+1.632993161855453*F_0[13]*dvpar*wvpar_cu+0.2449489742783177*F_0[13]*dvpar_cu*wvpar+0.02258769757263127*F_0[21]*dvpar_qu+0.01767766952966368*F_0[7]*dvpar_qu); 
  out[59] += mass*volFact*(0.6324555320336759*F_0[22]*dvpar_sq*wvpar_sq+0.7071067811865475*F_0[8]*dvpar_sq*wvpar_sq+1.414213562373095*F_0[8]*wvpar_qu+1.632993161855453*F_0[14]*dvpar*wvpar_cu+0.2449489742783177*F_0[14]*dvpar_cu*wvpar+0.02258769757263127*F_0[22]*dvpar_qu+0.01767766952966368*F_0[8]*dvpar_qu); 
  out[60] += mass*volFact*(0.632455532033676*F_0[24]*dvpar_sq*wvpar_sq+0.7071067811865475*F_0[11]*dvpar_sq*wvpar_sq+1.414213562373095*F_0[11]*wvpar_qu+1.632993161855453*F_0[17]*dvpar*wvpar_cu+0.2449489742783177*F_0[17]*dvpar_cu*wvpar+0.02258769757263128*F_0[24]*dvpar_qu+0.01767766952966368*F_0[11]*dvpar_qu); 
  out[61] += mass*volFact*(0.632455532033676*F_0[25]*dvpar_sq*wvpar_sq+0.7071067811865475*F_0[12]*dvpar_sq*wvpar_sq+1.414213562373095*F_0[12]*wvpar_qu+1.632993161855453*F_0[18]*dvpar*wvpar_cu+0.2449489742783177*F_0[18]*dvpar_cu*wvpar+0.02258769757263128*F_0[25]*dvpar_qu+0.01767766952966368*F_0[12]*dvpar_qu); 
  out[62] += mass*volFact*(0.6324555320336759*F_0[26]*dvpar_sq*wvpar_sq+0.7071067811865475*F_0[20]*dvpar_sq*wvpar_sq+1.414213562373095*F_0[20]*wvpar_qu+1.632993161855453*F_0[23]*dvpar*wvpar_cu+0.2449489742783178*F_0[23]*dvpar_cu*wvpar+0.02258769757263127*F_0[26]*dvpar_qu+0.01767766952966368*F_0[20]*dvpar_qu); 
  out[63] += mass*volFact*(1.414213562373095*G_1[0]*wvpar_sq+0.8164965809277261*G_1[3]*dvpar*wvpar+0.105409255338946*G_1[9]*dvpar_sq+0.1178511301977579*G_1[0]*dvpar_sq); 
  out[64] += mass*volFact*(1.414213562373095*G_1[1]*wvpar_sq+0.8164965809277261*G_1[5]*dvpar*wvpar+0.105409255338946*G_1[15]*dvpar_sq+0.1178511301977579*G_1[1]*dvpar_sq); 
  out[65] += mass*volFact*(1.414213562373095*G_1[2]*wvpar_sq+0.8164965809277261*G_1[6]*dvpar*wvpar+0.105409255338946*G_1[16]*dvpar_sq+0.1178511301977579*G_1[2]*dvpar_sq); 
  out[66] += mass*volFact*(1.414213562373095*G_1[4]*wvpar_sq+0.8164965809277261*G_1[10]*dvpar*wvpar+0.105409255338946*G_1[19]*dvpar_sq+0.1178511301977579*G_1[4]*dvpar_sq); 
  out[67] += mass*volFact*(1.414213562373095*G_1[7]*wvpar_sq+0.816496580927726*G_1[13]*dvpar*wvpar+0.105409255338946*G_1[21]*dvpar_sq+0.1178511301977579*G_1[7]*dvpar_sq); 
  out[68] += mass*volFact*(1.414213562373095*G_1[8]*wvpar_sq+0.816496580927726*G_1[14]*dvpar*wvpar+0.105409255338946*G_1[22]*dvpar_sq+0.1178511301977579*G_1[8]*dvpar_sq); 
  out[69] += mass*volFact*(1.414213562373095*G_1[11]*wvpar_sq+0.816496580927726*G_1[17]*dvpar*wvpar+0.105409255338946*G_1[24]*dvpar_sq+0.1178511301977579*G_1[11]*dvpar_sq); 
  out[70] += mass*volFact*(1.414213562373095*G_1[12]*wvpar_sq+0.816496580927726*G_1[18]*dvpar*wvpar+0.105409255338946*G_1[25]*dvpar_sq+0.1178511301977579*G_1[12]*dvpar_sq); 
  out[71] += mass*volFact*(1.414213562373095*G_1[20]*wvpar_sq+0.8164965809277261*G_1[23]*dvpar*wvpar+0.105409255338946*G_1[26]*dvpar_sq+0.1178511301977579*G_1[20]*dvpar_sq); 
} 
