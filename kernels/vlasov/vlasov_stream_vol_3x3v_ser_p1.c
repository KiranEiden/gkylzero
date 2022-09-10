#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_stream_vol_3x3v_ser_p1(const double *w, const double *dxv, const double *qmem, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // qmem:      q/m*EM fields (unused in streaming-only update).
  // f:         Input distribution function.
  // out:       Incremented output.
  double w3Ddx0  = w[3]/dxv[0]; 
  double dv3Ddx0 = dxv[3]/dxv[0]; 
  double w4Ddx1  = w[4]/dxv[1]; 
  double dv4Ddx1 = dxv[4]/dxv[1]; 
  double w5Ddx2  = w[5]/dxv[2]; 
  double dv5Ddx2 = dxv[5]/dxv[2]; 

  out[1] += 3.464101615137754*f[0]*w3Ddx0+f[4]*dv3Ddx0; 
  out[2] += 3.464101615137754*f[0]*w4Ddx1+f[5]*dv4Ddx1; 
  out[3] += 3.464101615137754*f[0]*w5Ddx2+f[6]*dv5Ddx2; 
  out[7] += 3.464101615137754*f[1]*w4Ddx1+3.464101615137754*f[2]*w3Ddx0+f[13]*dv4Ddx1+f[11]*dv3Ddx0; 
  out[8] += 3.464101615137754*f[1]*w5Ddx2+3.464101615137754*f[3]*w3Ddx0+f[17]*dv5Ddx2+f[12]*dv3Ddx0; 
  out[9] += 3.464101615137754*f[2]*w5Ddx2+3.464101615137754*f[3]*w4Ddx1+f[18]*dv5Ddx2+f[15]*dv4Ddx1; 
  out[10] += 3.464101615137754*f[4]*w3Ddx0+(0.8944271909999159*f[64]+f[0])*dv3Ddx0; 
  out[11] += 3.464101615137754*f[4]*w4Ddx1+f[16]*dv4Ddx1; 
  out[12] += 3.464101615137754*f[4]*w5Ddx2+f[20]*dv5Ddx2; 
  out[13] += 3.464101615137754*f[5]*w3Ddx0+f[16]*dv3Ddx0; 
  out[14] += 3.464101615137754*f[5]*w4Ddx1+(0.8944271909999159*f[96]+f[0])*dv4Ddx1; 
  out[15] += 3.464101615137754*f[5]*w5Ddx2+f[21]*dv5Ddx2; 
  out[17] += 3.464101615137754*f[6]*w3Ddx0+f[20]*dv3Ddx0; 
  out[18] += 3.464101615137754*f[6]*w4Ddx1+f[21]*dv4Ddx1; 
  out[19] += 3.464101615137754*f[6]*w5Ddx2+(0.8944271909999159*f[128]+f[0])*dv5Ddx2; 
  out[22] += 3.464101615137754*f[7]*w5Ddx2+3.464101615137754*f[8]*w4Ddx1+3.464101615137754*f[9]*w3Ddx0+f[32]*dv5Ddx2+f[27]*dv4Ddx1+f[25]*dv3Ddx0; 
  out[23] += 3.464101615137754*f[10]*w4Ddx1+3.464101615137754*f[11]*w3Ddx0+f[29]*dv4Ddx1+(0.8944271909999161*f[66]+f[2])*dv3Ddx0; 
  out[24] += 3.464101615137754*f[10]*w5Ddx2+3.464101615137754*f[12]*w3Ddx0+f[35]*dv5Ddx2+(0.8944271909999161*f[67]+f[3])*dv3Ddx0; 
  out[25] += 3.464101615137754*f[11]*w5Ddx2+3.464101615137754*f[12]*w4Ddx1+f[36]*dv5Ddx2+f[31]*dv4Ddx1; 
  out[26] += 3.464101615137754*f[13]*w4Ddx1+3.464101615137754*f[14]*w3Ddx0+(0.8944271909999161*f[97]+f[1])*dv4Ddx1+f[30]*dv3Ddx0; 
  out[27] += 3.464101615137754*f[13]*w5Ddx2+3.464101615137754*f[15]*w3Ddx0+f[38]*dv5Ddx2+f[31]*dv3Ddx0; 
  out[28] += 3.464101615137754*f[14]*w5Ddx2+3.464101615137754*f[15]*w4Ddx1+f[39]*dv5Ddx2+(0.8944271909999161*f[99]+f[3])*dv4Ddx1; 
  out[29] += 3.464101615137754*f[16]*w3Ddx0+(0.8944271909999161*f[68]+f[5])*dv3Ddx0; 
  out[30] += 3.464101615137754*f[16]*w4Ddx1+(0.8944271909999161*f[100]+f[4])*dv4Ddx1; 
  out[31] += 3.464101615137754*f[16]*w5Ddx2+f[41]*dv5Ddx2; 
  out[32] += 3.464101615137754*f[17]*w4Ddx1+3.464101615137754*f[18]*w3Ddx0+f[38]*dv4Ddx1+f[36]*dv3Ddx0; 
  out[33] += 3.464101615137754*f[17]*w5Ddx2+3.464101615137754*f[19]*w3Ddx0+(0.8944271909999161*f[129]+f[1])*dv5Ddx2+f[37]*dv3Ddx0; 
  out[34] += 3.464101615137754*f[18]*w5Ddx2+3.464101615137754*f[19]*w4Ddx1+(0.8944271909999161*f[130]+f[2])*dv5Ddx2+f[40]*dv4Ddx1; 
  out[35] += 3.464101615137754*f[20]*w3Ddx0+(0.8944271909999161*f[69]+f[6])*dv3Ddx0; 
  out[36] += 3.464101615137754*f[20]*w4Ddx1+f[41]*dv4Ddx1; 
  out[37] += 3.464101615137754*f[20]*w5Ddx2+(0.8944271909999161*f[132]+f[4])*dv5Ddx2; 
  out[38] += 3.464101615137754*f[21]*w3Ddx0+f[41]*dv3Ddx0; 
  out[39] += 3.464101615137754*f[21]*w4Ddx1+(0.8944271909999161*f[101]+f[6])*dv4Ddx1; 
  out[40] += 3.464101615137754*f[21]*w5Ddx2+(0.8944271909999161*f[133]+f[5])*dv5Ddx2; 
  out[42] += 3.464101615137754*f[23]*w5Ddx2+3.464101615137754*f[24]*w4Ddx1+3.464101615137754*f[25]*w3Ddx0+f[48]*dv5Ddx2+f[45]*dv4Ddx1+(0.8944271909999159*f[72]+f[9])*dv3Ddx0; 
  out[43] += 3.464101615137754*f[26]*w5Ddx2+3.464101615137754*f[27]*w4Ddx1+3.464101615137754*f[28]*w3Ddx0+f[51]*dv5Ddx2+(0.8944271909999159*f[103]+f[8])*dv4Ddx1+f[46]*dv3Ddx0; 
  out[44] += 3.464101615137754*f[29]*w4Ddx1+3.464101615137754*f[30]*w3Ddx0+(0.8944271909999159*f[105]+f[10])*dv4Ddx1+(0.8944271909999159*f[74]+f[14])*dv3Ddx0; 
  out[45] += 3.464101615137754*f[29]*w5Ddx2+3.464101615137754*f[31]*w3Ddx0+f[54]*dv5Ddx2+(0.8944271909999159*f[75]+f[15])*dv3Ddx0; 
  out[46] += 3.464101615137754*f[30]*w5Ddx2+3.464101615137754*f[31]*w4Ddx1+f[55]*dv5Ddx2+(0.8944271909999159*f[107]+f[12])*dv4Ddx1; 
  out[47] += 3.464101615137754*f[32]*w5Ddx2+3.464101615137754*f[33]*w4Ddx1+3.464101615137754*f[34]*w3Ddx0+(0.8944271909999159*f[134]+f[7])*dv5Ddx2+f[52]*dv4Ddx1+f[50]*dv3Ddx0; 
  out[48] += 3.464101615137754*f[35]*w4Ddx1+3.464101615137754*f[36]*w3Ddx0+f[54]*dv4Ddx1+(0.8944271909999159*f[77]+f[18])*dv3Ddx0; 
  out[49] += 3.464101615137754*f[35]*w5Ddx2+3.464101615137754*f[37]*w3Ddx0+(0.8944271909999159*f[137]+f[10])*dv5Ddx2+(0.8944271909999159*f[78]+f[19])*dv3Ddx0; 
  out[50] += 3.464101615137754*f[36]*w5Ddx2+3.464101615137754*f[37]*w4Ddx1+(0.8944271909999159*f[138]+f[11])*dv5Ddx2+f[56]*dv4Ddx1; 
  out[51] += 3.464101615137754*f[38]*w4Ddx1+3.464101615137754*f[39]*w3Ddx0+(0.8944271909999159*f[108]+f[17])*dv4Ddx1+f[55]*dv3Ddx0; 
  out[52] += 3.464101615137754*f[38]*w5Ddx2+3.464101615137754*f[40]*w3Ddx0+(0.8944271909999159*f[140]+f[13])*dv5Ddx2+f[56]*dv3Ddx0; 
  out[53] += 3.464101615137754*f[39]*w5Ddx2+3.464101615137754*f[40]*w4Ddx1+(0.8944271909999159*f[141]+f[14])*dv5Ddx2+(0.8944271909999159*f[110]+f[19])*dv4Ddx1; 
  out[54] += 3.464101615137754*f[41]*w3Ddx0+(0.8944271909999159*f[79]+f[21])*dv3Ddx0; 
  out[55] += 3.464101615137754*f[41]*w4Ddx1+(0.8944271909999159*f[111]+f[20])*dv4Ddx1; 
  out[56] += 3.464101615137754*f[41]*w5Ddx2+(0.8944271909999159*f[143]+f[16])*dv5Ddx2; 
  out[57] += 3.464101615137754*f[44]*w5Ddx2+3.464101615137754*f[45]*w4Ddx1+3.464101615137754*f[46]*w3Ddx0+f[60]*dv5Ddx2+(0.8944271909999161*f[114]+f[24])*dv4Ddx1+(0.8944271909999161*f[83]+f[28])*dv3Ddx0; 
  out[58] += 3.464101615137754*f[48]*w5Ddx2+3.464101615137754*f[49]*w4Ddx1+3.464101615137754*f[50]*w3Ddx0+(0.8944271909999161*f[145]+f[23])*dv5Ddx2+f[61]*dv4Ddx1+(0.8944271909999161*f[86]+f[34])*dv3Ddx0; 
  out[59] += 3.464101615137754*f[51]*w5Ddx2+3.464101615137754*f[52]*w4Ddx1+3.464101615137754*f[53]*w3Ddx0+(0.8944271909999161*f[148]+f[26])*dv5Ddx2+(0.8944271909999161*f[117]+f[33])*dv4Ddx1+f[62]*dv3Ddx0; 
  out[60] += 3.464101615137754*f[54]*w4Ddx1+3.464101615137754*f[55]*w3Ddx0+(0.8944271909999161*f[119]+f[35])*dv4Ddx1+(0.8944271909999161*f[88]+f[39])*dv3Ddx0; 
  out[61] += 3.464101615137754*f[54]*w5Ddx2+3.464101615137754*f[56]*w3Ddx0+(0.8944271909999161*f[151]+f[29])*dv5Ddx2+(0.8944271909999161*f[89]+f[40])*dv3Ddx0; 
  out[62] += 3.464101615137754*f[55]*w5Ddx2+3.464101615137754*f[56]*w4Ddx1+(0.8944271909999161*f[152]+f[30])*dv5Ddx2+(0.8944271909999161*f[121]+f[37])*dv4Ddx1; 
  out[63] += 3.464101615137754*f[60]*w5Ddx2+3.464101615137754*f[61]*w4Ddx1+3.464101615137754*f[62]*w3Ddx0+(0.8944271909999159*f[156]+f[44])*dv5Ddx2+(0.8944271909999159*f[125]+f[49])*dv4Ddx1+(0.8944271909999159*f[94]+f[53])*dv3Ddx0; 
  out[65] += 3.464101615137755*f[64]*w3Ddx0+0.8944271909999161*f[4]*dv3Ddx0; 
  out[66] += 3.464101615137755*f[64]*w4Ddx1+f[68]*dv4Ddx1; 
  out[67] += 3.464101615137755*f[64]*w5Ddx2+f[69]*dv5Ddx2; 
  out[70] += 3.464101615137755*f[65]*w4Ddx1+3.464101615137755*f[66]*w3Ddx0+f[73]*dv4Ddx1+0.8944271909999159*f[11]*dv3Ddx0; 
  out[71] += 3.464101615137755*f[65]*w5Ddx2+3.464101615137755*f[67]*w3Ddx0+f[76]*dv5Ddx2+0.8944271909999159*f[12]*dv3Ddx0; 
  out[72] += 3.464101615137755*f[66]*w5Ddx2+3.464101615137755*f[67]*w4Ddx1+f[77]*dv5Ddx2+f[75]*dv4Ddx1; 
  out[73] += 3.464101615137755*f[68]*w3Ddx0+0.8944271909999159*f[16]*dv3Ddx0; 
  out[74] += 3.464101615137755*f[68]*w4Ddx1+f[64]*dv4Ddx1; 
  out[75] += 3.464101615137755*f[68]*w5Ddx2+f[79]*dv5Ddx2; 
  out[76] += 3.464101615137755*f[69]*w3Ddx0+0.8944271909999159*f[20]*dv3Ddx0; 
  out[77] += 3.464101615137755*f[69]*w4Ddx1+f[79]*dv4Ddx1; 
  out[78] += 3.464101615137755*f[69]*w5Ddx2+f[64]*dv5Ddx2; 
  out[80] += 3.464101615137755*f[70]*w5Ddx2+3.464101615137755*f[71]*w4Ddx1+3.464101615137755*f[72]*w3Ddx0+f[84]*dv5Ddx2+f[82]*dv4Ddx1+0.8944271909999161*f[25]*dv3Ddx0; 
  out[81] += 3.464101615137755*f[73]*w4Ddx1+3.464101615137755*f[74]*w3Ddx0+f[65]*dv4Ddx1+0.8944271909999161*f[30]*dv3Ddx0; 
  out[82] += 3.464101615137755*f[73]*w5Ddx2+3.464101615137755*f[75]*w3Ddx0+f[87]*dv5Ddx2+0.8944271909999161*f[31]*dv3Ddx0; 
  out[83] += 3.464101615137755*f[74]*w5Ddx2+3.464101615137755*f[75]*w4Ddx1+f[88]*dv5Ddx2+f[67]*dv4Ddx1; 
  out[84] += 3.464101615137755*f[76]*w4Ddx1+3.464101615137755*f[77]*w3Ddx0+f[87]*dv4Ddx1+0.8944271909999161*f[36]*dv3Ddx0; 
  out[85] += 3.464101615137755*f[76]*w5Ddx2+3.464101615137755*f[78]*w3Ddx0+f[65]*dv5Ddx2+0.8944271909999161*f[37]*dv3Ddx0; 
  out[86] += 3.464101615137755*f[77]*w5Ddx2+3.464101615137755*f[78]*w4Ddx1+f[66]*dv5Ddx2+f[89]*dv4Ddx1; 
  out[87] += 3.464101615137755*f[79]*w3Ddx0+0.8944271909999161*f[41]*dv3Ddx0; 
  out[88] += 3.464101615137755*f[79]*w4Ddx1+f[69]*dv4Ddx1; 
  out[89] += 3.464101615137755*f[79]*w5Ddx2+f[68]*dv5Ddx2; 
  out[90] += 3.464101615137755*f[81]*w5Ddx2+3.464101615137755*f[82]*w4Ddx1+3.464101615137755*f[83]*w3Ddx0+f[92]*dv5Ddx2+f[71]*dv4Ddx1+0.8944271909999159*f[46]*dv3Ddx0; 
  out[91] += 3.464101615137755*f[84]*w5Ddx2+3.464101615137755*f[85]*w4Ddx1+3.464101615137755*f[86]*w3Ddx0+f[70]*dv5Ddx2+f[93]*dv4Ddx1+0.8944271909999159*f[50]*dv3Ddx0; 
  out[92] += 3.464101615137755*f[87]*w4Ddx1+3.464101615137755*f[88]*w3Ddx0+f[76]*dv4Ddx1+0.8944271909999159*f[55]*dv3Ddx0; 
  out[93] += 3.464101615137755*f[87]*w5Ddx2+3.464101615137755*f[89]*w3Ddx0+f[73]*dv5Ddx2+0.8944271909999159*f[56]*dv3Ddx0; 
  out[94] += 3.464101615137755*f[88]*w5Ddx2+3.464101615137755*f[89]*w4Ddx1+f[74]*dv5Ddx2+f[78]*dv4Ddx1; 
  out[95] += 3.464101615137755*f[92]*w5Ddx2+3.464101615137755*f[93]*w4Ddx1+3.464101615137755*f[94]*w3Ddx0+f[81]*dv5Ddx2+f[85]*dv4Ddx1+0.8944271909999161*f[62]*dv3Ddx0; 
  out[97] += 3.464101615137755*f[96]*w3Ddx0+f[100]*dv3Ddx0; 
  out[98] += 3.464101615137755*f[96]*w4Ddx1+0.8944271909999161*f[5]*dv4Ddx1; 
  out[99] += 3.464101615137755*f[96]*w5Ddx2+f[101]*dv5Ddx2; 
  out[102] += 3.464101615137755*f[97]*w4Ddx1+3.464101615137755*f[98]*w3Ddx0+0.8944271909999159*f[13]*dv4Ddx1+f[106]*dv3Ddx0; 
  out[103] += 3.464101615137755*f[97]*w5Ddx2+3.464101615137755*f[99]*w3Ddx0+f[108]*dv5Ddx2+f[107]*dv3Ddx0; 
  out[104] += 3.464101615137755*f[98]*w5Ddx2+3.464101615137755*f[99]*w4Ddx1+f[109]*dv5Ddx2+0.8944271909999159*f[15]*dv4Ddx1; 
  out[105] += 3.464101615137755*f[100]*w3Ddx0+f[96]*dv3Ddx0; 
  out[106] += 3.464101615137755*f[100]*w4Ddx1+0.8944271909999159*f[16]*dv4Ddx1; 
  out[107] += 3.464101615137755*f[100]*w5Ddx2+f[111]*dv5Ddx2; 
  out[108] += 3.464101615137755*f[101]*w3Ddx0+f[111]*dv3Ddx0; 
  out[109] += 3.464101615137755*f[101]*w4Ddx1+0.8944271909999159*f[21]*dv4Ddx1; 
  out[110] += 3.464101615137755*f[101]*w5Ddx2+f[96]*dv5Ddx2; 
  out[112] += 3.464101615137755*f[102]*w5Ddx2+3.464101615137755*f[103]*w4Ddx1+3.464101615137755*f[104]*w3Ddx0+f[116]*dv5Ddx2+0.8944271909999161*f[27]*dv4Ddx1+f[115]*dv3Ddx0; 
  out[113] += 3.464101615137755*f[105]*w4Ddx1+3.464101615137755*f[106]*w3Ddx0+0.8944271909999161*f[29]*dv4Ddx1+f[98]*dv3Ddx0; 
  out[114] += 3.464101615137755*f[105]*w5Ddx2+3.464101615137755*f[107]*w3Ddx0+f[119]*dv5Ddx2+f[99]*dv3Ddx0; 
  out[115] += 3.464101615137755*f[106]*w5Ddx2+3.464101615137755*f[107]*w4Ddx1+f[120]*dv5Ddx2+0.8944271909999161*f[31]*dv4Ddx1; 
  out[116] += 3.464101615137755*f[108]*w4Ddx1+3.464101615137755*f[109]*w3Ddx0+0.8944271909999161*f[38]*dv4Ddx1+f[120]*dv3Ddx0; 
  out[117] += 3.464101615137755*f[108]*w5Ddx2+3.464101615137755*f[110]*w3Ddx0+f[97]*dv5Ddx2+f[121]*dv3Ddx0; 
  out[118] += 3.464101615137755*f[109]*w5Ddx2+3.464101615137755*f[110]*w4Ddx1+f[98]*dv5Ddx2+0.8944271909999161*f[40]*dv4Ddx1; 
  out[119] += 3.464101615137755*f[111]*w3Ddx0+f[101]*dv3Ddx0; 
  out[120] += 3.464101615137755*f[111]*w4Ddx1+0.8944271909999161*f[41]*dv4Ddx1; 
  out[121] += 3.464101615137755*f[111]*w5Ddx2+f[100]*dv5Ddx2; 
  out[122] += 3.464101615137755*f[113]*w5Ddx2+3.464101615137755*f[114]*w4Ddx1+3.464101615137755*f[115]*w3Ddx0+f[124]*dv5Ddx2+0.8944271909999159*f[45]*dv4Ddx1+f[104]*dv3Ddx0; 
  out[123] += 3.464101615137755*f[116]*w5Ddx2+3.464101615137755*f[117]*w4Ddx1+3.464101615137755*f[118]*w3Ddx0+f[102]*dv5Ddx2+0.8944271909999159*f[52]*dv4Ddx1+f[126]*dv3Ddx0; 
  out[124] += 3.464101615137755*f[119]*w4Ddx1+3.464101615137755*f[120]*w3Ddx0+0.8944271909999159*f[54]*dv4Ddx1+f[109]*dv3Ddx0; 
  out[125] += 3.464101615137755*f[119]*w5Ddx2+3.464101615137755*f[121]*w3Ddx0+f[105]*dv5Ddx2+f[110]*dv3Ddx0; 
  out[126] += 3.464101615137755*f[120]*w5Ddx2+3.464101615137755*f[121]*w4Ddx1+f[106]*dv5Ddx2+0.8944271909999159*f[56]*dv4Ddx1; 
  out[127] += 3.464101615137755*f[124]*w5Ddx2+3.464101615137755*f[125]*w4Ddx1+3.464101615137755*f[126]*w3Ddx0+f[113]*dv5Ddx2+0.8944271909999161*f[61]*dv4Ddx1+f[118]*dv3Ddx0; 
  out[129] += 3.464101615137755*f[128]*w3Ddx0+f[132]*dv3Ddx0; 
  out[130] += 3.464101615137755*f[128]*w4Ddx1+f[133]*dv4Ddx1; 
  out[131] += 3.464101615137755*f[128]*w5Ddx2+0.8944271909999161*f[6]*dv5Ddx2; 
  out[134] += 3.464101615137755*f[129]*w4Ddx1+3.464101615137755*f[130]*w3Ddx0+f[140]*dv4Ddx1+f[138]*dv3Ddx0; 
  out[135] += 3.464101615137755*f[129]*w5Ddx2+3.464101615137755*f[131]*w3Ddx0+0.8944271909999159*f[17]*dv5Ddx2+f[139]*dv3Ddx0; 
  out[136] += 3.464101615137755*f[130]*w5Ddx2+3.464101615137755*f[131]*w4Ddx1+0.8944271909999159*f[18]*dv5Ddx2+f[142]*dv4Ddx1; 
  out[137] += 3.464101615137755*f[132]*w3Ddx0+f[128]*dv3Ddx0; 
  out[138] += 3.464101615137755*f[132]*w4Ddx1+f[143]*dv4Ddx1; 
  out[139] += 3.464101615137755*f[132]*w5Ddx2+0.8944271909999159*f[20]*dv5Ddx2; 
  out[140] += 3.464101615137755*f[133]*w3Ddx0+f[143]*dv3Ddx0; 
  out[141] += 3.464101615137755*f[133]*w4Ddx1+f[128]*dv4Ddx1; 
  out[142] += 3.464101615137755*f[133]*w5Ddx2+0.8944271909999159*f[21]*dv5Ddx2; 
  out[144] += 3.464101615137755*f[134]*w5Ddx2+3.464101615137755*f[135]*w4Ddx1+3.464101615137755*f[136]*w3Ddx0+0.8944271909999161*f[32]*dv5Ddx2+f[149]*dv4Ddx1+f[147]*dv3Ddx0; 
  out[145] += 3.464101615137755*f[137]*w4Ddx1+3.464101615137755*f[138]*w3Ddx0+f[151]*dv4Ddx1+f[130]*dv3Ddx0; 
  out[146] += 3.464101615137755*f[137]*w5Ddx2+3.464101615137755*f[139]*w3Ddx0+0.8944271909999161*f[35]*dv5Ddx2+f[131]*dv3Ddx0; 
  out[147] += 3.464101615137755*f[138]*w5Ddx2+3.464101615137755*f[139]*w4Ddx1+0.8944271909999161*f[36]*dv5Ddx2+f[153]*dv4Ddx1; 
  out[148] += 3.464101615137755*f[140]*w4Ddx1+3.464101615137755*f[141]*w3Ddx0+f[129]*dv4Ddx1+f[152]*dv3Ddx0; 
  out[149] += 3.464101615137755*f[140]*w5Ddx2+3.464101615137755*f[142]*w3Ddx0+0.8944271909999161*f[38]*dv5Ddx2+f[153]*dv3Ddx0; 
  out[150] += 3.464101615137755*f[141]*w5Ddx2+3.464101615137755*f[142]*w4Ddx1+0.8944271909999161*f[39]*dv5Ddx2+f[131]*dv4Ddx1; 
  out[151] += 3.464101615137755*f[143]*w3Ddx0+f[133]*dv3Ddx0; 
  out[152] += 3.464101615137755*f[143]*w4Ddx1+f[132]*dv4Ddx1; 
  out[153] += 3.464101615137755*f[143]*w5Ddx2+0.8944271909999161*f[41]*dv5Ddx2; 
  out[154] += 3.464101615137755*f[145]*w5Ddx2+3.464101615137755*f[146]*w4Ddx1+3.464101615137755*f[147]*w3Ddx0+0.8944271909999159*f[48]*dv5Ddx2+f[157]*dv4Ddx1+f[136]*dv3Ddx0; 
  out[155] += 3.464101615137755*f[148]*w5Ddx2+3.464101615137755*f[149]*w4Ddx1+3.464101615137755*f[150]*w3Ddx0+0.8944271909999159*f[51]*dv5Ddx2+f[135]*dv4Ddx1+f[158]*dv3Ddx0; 
  out[156] += 3.464101615137755*f[151]*w4Ddx1+3.464101615137755*f[152]*w3Ddx0+f[137]*dv4Ddx1+f[141]*dv3Ddx0; 
  out[157] += 3.464101615137755*f[151]*w5Ddx2+3.464101615137755*f[153]*w3Ddx0+0.8944271909999159*f[54]*dv5Ddx2+f[142]*dv3Ddx0; 
  out[158] += 3.464101615137755*f[152]*w5Ddx2+3.464101615137755*f[153]*w4Ddx1+0.8944271909999159*f[55]*dv5Ddx2+f[139]*dv4Ddx1; 
  out[159] += 3.464101615137755*f[156]*w5Ddx2+3.464101615137755*f[157]*w4Ddx1+3.464101615137755*f[158]*w3Ddx0+0.8944271909999161*f[60]*dv5Ddx2+f[146]*dv4Ddx1+f[150]*dv3Ddx0; 

  return 3.0*(fabs(w3Ddx0)+0.5*dv3Ddx0+fabs(w4Ddx1)+0.5*dv4Ddx1+fabs(w5Ddx2)+0.5*dv5Ddx2);
} 
