#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_stream_vol_2x3v_tensor_p2(const double *w, const double *dxv, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // f:         Input distribution function.
  // out:       Incremented output.
  double w2Ddx0  = w[2]/dxv[0]; 
  double dv2Ddx0 = dxv[2]/dxv[0]; 
  double w3Ddx1  = w[3]/dxv[1]; 
  double dv3Ddx1 = dxv[3]/dxv[1]; 

  out[1] += 3.464101615137754*f[0]*w2Ddx0+f[3]*dv2Ddx0; 
  out[2] += 3.464101615137754*f[0]*w3Ddx1+f[4]*dv3Ddx1; 
  out[6] += 3.464101615137754*f[1]*w3Ddx1+3.464101615137754*f[2]*w2Ddx0+f[9]*dv3Ddx1+f[8]*dv2Ddx0; 
  out[7] += 3.464101615137754*f[3]*w2Ddx0+(0.8944271909999159*f[18]+f[0])*dv2Ddx0; 
  out[8] += 3.464101615137754*f[3]*w3Ddx1+f[11]*dv3Ddx1; 
  out[9] += 3.464101615137754*f[4]*w2Ddx0+f[11]*dv2Ddx0; 
  out[10] += 3.464101615137754*f[4]*w3Ddx1+(0.8944271909999159*f[19]+f[0])*dv3Ddx1; 
  out[12] += 3.464101615137754*f[5]*w2Ddx0+f[14]*dv2Ddx0; 
  out[13] += 3.464101615137754*f[5]*w3Ddx1+f[15]*dv3Ddx1; 
  out[16] += 7.745966692414834*f[1]*w2Ddx0+2.23606797749979*f[7]*dv2Ddx0; 
  out[17] += 7.745966692414834*f[2]*w3Ddx1+2.23606797749979*f[10]*dv3Ddx1; 
  out[21] += 3.464101615137754*f[7]*w3Ddx1+3.464101615137754*f[8]*w2Ddx0+f[23]*dv3Ddx1+(0.8944271909999161*f[36]+f[2])*dv2Ddx0; 
  out[22] += 3.464101615137754*f[9]*w3Ddx1+3.464101615137754*f[10]*w2Ddx0+(0.8944271909999161*f[40]+f[1])*dv3Ddx1+f[24]*dv2Ddx0; 
  out[23] += 3.464101615137754*f[11]*w2Ddx0+(0.8944271909999161*f[39]+f[4])*dv2Ddx0; 
  out[24] += 3.464101615137754*f[11]*w3Ddx1+(0.8944271909999161*f[42]+f[3])*dv3Ddx1; 
  out[25] += 3.464101615137754*f[12]*w3Ddx1+3.464101615137754*f[13]*w2Ddx0+f[28]*dv3Ddx1+f[27]*dv2Ddx0; 
  out[26] += 3.464101615137754*f[14]*w2Ddx0+(0.8944271909999161*f[45]+f[5])*dv2Ddx0; 
  out[27] += 3.464101615137754*f[14]*w3Ddx1+f[30]*dv3Ddx1; 
  out[28] += 3.464101615137754*f[15]*w2Ddx0+f[30]*dv2Ddx0; 
  out[29] += 3.464101615137754*f[15]*w3Ddx1+(0.8944271909999161*f[46]+f[5])*dv3Ddx1; 
  out[31] += 3.464101615137755*f[16]*w3Ddx1+7.745966692414834*f[6]*w2Ddx0+f[37]*dv3Ddx1+2.23606797749979*f[21]*dv2Ddx0; 
  out[32] += 7.745966692414834*f[6]*w3Ddx1+3.464101615137755*f[17]*w2Ddx0+2.23606797749979*f[22]*dv3Ddx1+f[34]*dv2Ddx0; 
  out[33] += 7.745966692414834*f[7]*w2Ddx0+(2.0*f[35]+2.23606797749979*f[1])*dv2Ddx0; 
  out[34] += 7.745966692414834*f[8]*w3Ddx1+2.23606797749979*f[24]*dv3Ddx1; 
  out[35] += 3.464101615137755*f[18]*w2Ddx0+0.8944271909999161*f[3]*dv2Ddx0; 
  out[36] += 3.464101615137755*f[18]*w3Ddx1+f[39]*dv3Ddx1; 
  out[37] += 7.745966692414834*f[9]*w2Ddx0+2.23606797749979*f[23]*dv2Ddx0; 
  out[38] += 7.745966692414834*f[10]*w3Ddx1+(2.0*f[41]+2.23606797749979*f[2])*dv3Ddx1; 
  out[40] += 3.464101615137755*f[19]*w2Ddx0+f[42]*dv2Ddx0; 
  out[41] += 3.464101615137755*f[19]*w3Ddx1+0.8944271909999161*f[4]*dv3Ddx1; 
  out[43] += 7.745966692414834*f[12]*w2Ddx0+2.23606797749979*f[26]*dv2Ddx0; 
  out[44] += 7.745966692414834*f[13]*w3Ddx1+2.23606797749979*f[29]*dv3Ddx1; 
  out[47] += 3.464101615137755*f[20]*w2Ddx0+f[49]*dv2Ddx0; 
  out[48] += 3.464101615137755*f[20]*w3Ddx1+f[50]*dv3Ddx1; 
  out[51] += 3.464101615137754*f[23]*w3Ddx1+3.464101615137754*f[24]*w2Ddx0+(0.8944271909999159*f[66]+f[7])*dv3Ddx1+(0.8944271909999159*f[64]+f[10])*dv2Ddx0; 
  out[52] += 3.464101615137754*f[26]*w3Ddx1+3.464101615137754*f[27]*w2Ddx0+f[54]*dv3Ddx1+(0.8944271909999159*f[73]+f[13])*dv2Ddx0; 
  out[53] += 3.464101615137754*f[28]*w3Ddx1+3.464101615137754*f[29]*w2Ddx0+(0.8944271909999159*f[77]+f[12])*dv3Ddx1+f[55]*dv2Ddx0; 
  out[54] += 3.464101615137754*f[30]*w2Ddx0+(0.8944271909999159*f[76]+f[15])*dv2Ddx0; 
  out[55] += 3.464101615137754*f[30]*w3Ddx1+(0.8944271909999159*f[79]+f[14])*dv3Ddx1; 
  out[56] += 3.464101615137755*f[33]*w3Ddx1+7.745966692414834*f[21]*w2Ddx0+f[61]*dv3Ddx1+(2.0*f[58]+2.23606797749979*f[6])*dv2Ddx0; 
  out[57] += 7.745966692414834*f[21]*w3Ddx1+3.464101615137755*f[34]*w2Ddx0+2.23606797749979*f[51]*dv3Ddx1+(0.8944271909999159*f[88]+f[17])*dv2Ddx0; 
  out[58] += 3.464101615137755*f[35]*w3Ddx1+3.464101615137755*f[36]*w2Ddx0+f[63]*dv3Ddx1+0.8944271909999159*f[8]*dv2Ddx0; 
  out[59] += 3.464101615137755*f[37]*w3Ddx1+7.745966692414834*f[22]*w2Ddx0+(0.8944271909999159*f[89]+f[16])*dv3Ddx1+2.23606797749979*f[51]*dv2Ddx0; 
  out[60] += 7.745966692414834*f[22]*w3Ddx1+3.464101615137755*f[38]*w2Ddx0+(2.0*f[65]+2.23606797749979*f[6])*dv3Ddx1+f[62]*dv2Ddx0; 
  out[61] += 7.745966692414834*f[23]*w2Ddx0+(2.0*f[63]+2.23606797749979*f[9])*dv2Ddx0; 
  out[62] += 7.745966692414834*f[24]*w3Ddx1+(2.0*f[67]+2.23606797749979*f[8])*dv3Ddx1; 
  out[63] += 3.464101615137755*f[39]*w2Ddx0+0.8944271909999159*f[11]*dv2Ddx0; 
  out[64] += 3.464101615137755*f[39]*w3Ddx1+(0.8944271909999159*f[91]+f[18])*dv3Ddx1; 
  out[65] += 3.464101615137755*f[40]*w3Ddx1+3.464101615137755*f[41]*w2Ddx0+0.8944271909999159*f[9]*dv3Ddx1+f[67]*dv2Ddx0; 
  out[66] += 3.464101615137755*f[42]*w2Ddx0+(0.8944271909999159*f[91]+f[19])*dv2Ddx0; 
  out[67] += 3.464101615137755*f[42]*w3Ddx1+0.8944271909999159*f[11]*dv3Ddx1; 
  out[68] += 3.464101615137755*f[43]*w3Ddx1+7.745966692414834*f[25]*w2Ddx0+f[74]*dv3Ddx1+2.23606797749979*f[52]*dv2Ddx0; 
  out[69] += 7.745966692414834*f[25]*w3Ddx1+3.464101615137755*f[44]*w2Ddx0+2.23606797749979*f[53]*dv3Ddx1+f[71]*dv2Ddx0; 
  out[70] += 7.745966692414834*f[26]*w2Ddx0+(2.0*f[72]+2.23606797749979*f[12])*dv2Ddx0; 
  out[71] += 7.745966692414834*f[27]*w3Ddx1+2.23606797749979*f[55]*dv3Ddx1; 
  out[72] += 3.464101615137755*f[45]*w2Ddx0+0.8944271909999159*f[14]*dv2Ddx0; 
  out[73] += 3.464101615137755*f[45]*w3Ddx1+f[76]*dv3Ddx1; 
  out[74] += 7.745966692414834*f[28]*w2Ddx0+2.23606797749979*f[54]*dv2Ddx0; 
  out[75] += 7.745966692414834*f[29]*w3Ddx1+(2.0*f[78]+2.23606797749979*f[13])*dv3Ddx1; 
  out[77] += 3.464101615137755*f[46]*w2Ddx0+f[79]*dv2Ddx0; 
  out[78] += 3.464101615137755*f[46]*w3Ddx1+0.8944271909999159*f[15]*dv3Ddx1; 
  out[80] += 3.464101615137755*f[47]*w3Ddx1+3.464101615137755*f[48]*w2Ddx0+f[83]*dv3Ddx1+f[82]*dv2Ddx0; 
  out[81] += 3.464101615137755*f[49]*w2Ddx0+(0.8944271909999159*f[94]+f[20])*dv2Ddx0; 
  out[82] += 3.464101615137755*f[49]*w3Ddx1+f[85]*dv3Ddx1; 
  out[83] += 3.464101615137755*f[50]*w2Ddx0+f[85]*dv2Ddx0; 
  out[84] += 3.464101615137755*f[50]*w3Ddx1+(0.8944271909999159*f[95]+f[20])*dv3Ddx1; 
  out[86] += 7.745966692414834*f[31]*w3Ddx1+7.745966692414834*f[32]*w2Ddx0+2.23606797749979*f[59]*dv3Ddx1+2.23606797749979*f[57]*dv2Ddx0; 
  out[87] += 7.745966692414834*f[35]*w2Ddx0+2.0*f[7]*dv2Ddx0; 
  out[88] += 7.745966692414834*f[36]*w3Ddx1+2.23606797749979*f[64]*dv3Ddx1; 
  out[89] += 7.745966692414834*f[40]*w2Ddx0+2.23606797749979*f[66]*dv2Ddx0; 
  out[90] += 7.745966692414834*f[41]*w3Ddx1+2.0*f[10]*dv3Ddx1; 
  out[92] += 7.745966692414834*f[47]*w2Ddx0+2.23606797749979*f[81]*dv2Ddx0; 
  out[93] += 7.745966692414834*f[48]*w3Ddx1+2.23606797749979*f[84]*dv3Ddx1; 
  out[96] += 3.464101615137754*f[54]*w3Ddx1+3.464101615137754*f[55]*w2Ddx0+(0.8944271909999161*f[111]+f[26])*dv3Ddx1+(0.8944271909999161*f[109]+f[29])*dv2Ddx0; 
  out[97] += 3.464101615137755*f[61]*w3Ddx1+7.745966692414834*f[51]*w2Ddx0+(0.8944271909999161*f[125]+f[33])*dv3Ddx1+(2.0*f[99]+2.23606797749979*f[22])*dv2Ddx0; 
  out[98] += 7.745966692414834*f[51]*w3Ddx1+3.464101615137755*f[62]*w2Ddx0+(2.0*f[100]+2.23606797749979*f[21])*dv3Ddx1+(0.8944271909999161*f[122]+f[38])*dv2Ddx0; 
  out[99] += 3.464101615137755*f[63]*w3Ddx1+3.464101615137755*f[64]*w2Ddx0+(0.8944271909999161*f[127]+f[35])*dv3Ddx1+0.8944271909999161*f[24]*dv2Ddx0; 
  out[100] += 3.464101615137755*f[66]*w3Ddx1+3.464101615137755*f[67]*w2Ddx0+0.8944271909999161*f[23]*dv3Ddx1+(0.8944271909999161*f[128]+f[41])*dv2Ddx0; 
  out[101] += 3.464101615137755*f[70]*w3Ddx1+7.745966692414834*f[52]*w2Ddx0+f[106]*dv3Ddx1+(2.0*f[103]+2.23606797749979*f[25])*dv2Ddx0; 
  out[102] += 7.745966692414834*f[52]*w3Ddx1+3.464101615137755*f[71]*w2Ddx0+2.23606797749979*f[96]*dv3Ddx1+(0.8944271909999161*f[131]+f[44])*dv2Ddx0; 
  out[103] += 3.464101615137755*f[72]*w3Ddx1+3.464101615137755*f[73]*w2Ddx0+f[108]*dv3Ddx1+0.8944271909999161*f[27]*dv2Ddx0; 
  out[104] += 3.464101615137755*f[74]*w3Ddx1+7.745966692414834*f[53]*w2Ddx0+(0.8944271909999161*f[132]+f[43])*dv3Ddx1+2.23606797749979*f[96]*dv2Ddx0; 
  out[105] += 7.745966692414834*f[53]*w3Ddx1+3.464101615137755*f[75]*w2Ddx0+(2.0*f[110]+2.23606797749979*f[25])*dv3Ddx1+f[107]*dv2Ddx0; 
  out[106] += 7.745966692414834*f[54]*w2Ddx0+(2.0*f[108]+2.23606797749979*f[28])*dv2Ddx0; 
  out[107] += 7.745966692414834*f[55]*w3Ddx1+(2.0*f[112]+2.23606797749979*f[27])*dv3Ddx1; 
  out[108] += 3.464101615137755*f[76]*w2Ddx0+0.8944271909999161*f[30]*dv2Ddx0; 
  out[109] += 3.464101615137755*f[76]*w3Ddx1+(0.8944271909999161*f[134]+f[45])*dv3Ddx1; 
  out[110] += 3.464101615137755*f[77]*w3Ddx1+3.464101615137755*f[78]*w2Ddx0+0.8944271909999161*f[28]*dv3Ddx1+f[112]*dv2Ddx0; 
  out[111] += 3.464101615137755*f[79]*w2Ddx0+(0.8944271909999161*f[134]+f[46])*dv2Ddx0; 
  out[112] += 3.464101615137755*f[79]*w3Ddx1+0.8944271909999161*f[30]*dv3Ddx1; 
  out[113] += 3.464101615137755*f[81]*w3Ddx1+3.464101615137755*f[82]*w2Ddx0+f[115]*dv3Ddx1+(0.8944271909999161*f[140]+f[48])*dv2Ddx0; 
  out[114] += 3.464101615137755*f[83]*w3Ddx1+3.464101615137755*f[84]*w2Ddx0+(0.8944271909999161*f[144]+f[47])*dv3Ddx1+f[116]*dv2Ddx0; 
  out[115] += 3.464101615137755*f[85]*w2Ddx0+(0.8944271909999161*f[143]+f[50])*dv2Ddx0; 
  out[116] += 3.464101615137755*f[85]*w3Ddx1+(0.8944271909999161*f[146]+f[49])*dv3Ddx1; 
  out[117] += 7.745966692414834*f[56]*w3Ddx1+7.745966692414834*f[57]*w2Ddx0+2.23606797749979*f[97]*dv3Ddx1+(2.0*f[119]+2.23606797749979*f[32])*dv2Ddx0; 
  out[118] += 3.464101615137754*f[87]*w3Ddx1+7.745966692414834*f[58]*w2Ddx0+f[121]*dv3Ddx1+2.0*f[21]*dv2Ddx0; 
  out[119] += 7.745966692414834*f[58]*w3Ddx1+3.464101615137754*f[88]*w2Ddx0+2.23606797749979*f[99]*dv3Ddx1+0.8944271909999161*f[34]*dv2Ddx0; 
  out[120] += 7.745966692414834*f[59]*w3Ddx1+7.745966692414834*f[60]*w2Ddx0+(2.0*f[123]+2.23606797749979*f[31])*dv3Ddx1+2.23606797749979*f[98]*dv2Ddx0; 
  out[121] += 7.745966692414834*f[63]*w2Ddx0+2.0*f[23]*dv2Ddx0; 
  out[122] += 7.745966692414834*f[64]*w3Ddx1+(2.0*f[128]+2.23606797749979*f[36])*dv3Ddx1; 
  out[123] += 3.464101615137754*f[89]*w3Ddx1+7.745966692414834*f[65]*w2Ddx0+0.8944271909999161*f[37]*dv3Ddx1+2.23606797749979*f[100]*dv2Ddx0; 
  out[124] += 7.745966692414834*f[65]*w3Ddx1+3.464101615137754*f[90]*w2Ddx0+2.0*f[22]*dv3Ddx1+f[126]*dv2Ddx0; 
  out[125] += 7.745966692414834*f[66]*w2Ddx0+(2.0*f[127]+2.23606797749979*f[40])*dv2Ddx0; 
  out[126] += 7.745966692414834*f[67]*w3Ddx1+2.0*f[24]*dv3Ddx1; 
  out[127] += 3.464101615137754*f[91]*w2Ddx0+0.8944271909999161*f[42]*dv2Ddx0; 
  out[128] += 3.464101615137754*f[91]*w3Ddx1+0.8944271909999161*f[39]*dv3Ddx1; 
  out[129] += 7.745966692414834*f[68]*w3Ddx1+7.745966692414834*f[69]*w2Ddx0+2.23606797749979*f[104]*dv3Ddx1+2.23606797749979*f[102]*dv2Ddx0; 
  out[130] += 7.745966692414834*f[72]*w2Ddx0+2.0*f[26]*dv2Ddx0; 
  out[131] += 7.745966692414834*f[73]*w3Ddx1+2.23606797749979*f[109]*dv3Ddx1; 
  out[132] += 7.745966692414834*f[77]*w2Ddx0+2.23606797749979*f[111]*dv2Ddx0; 
  out[133] += 7.745966692414834*f[78]*w3Ddx1+2.0*f[29]*dv3Ddx1; 
  out[135] += 3.464101615137754*f[92]*w3Ddx1+7.745966692414834*f[80]*w2Ddx0+f[141]*dv3Ddx1+2.23606797749979*f[113]*dv2Ddx0; 
  out[136] += 7.745966692414834*f[80]*w3Ddx1+3.464101615137754*f[93]*w2Ddx0+2.23606797749979*f[114]*dv3Ddx1+f[138]*dv2Ddx0; 
  out[137] += 7.745966692414834*f[81]*w2Ddx0+(2.0*f[139]+2.23606797749979*f[47])*dv2Ddx0; 
  out[138] += 7.745966692414834*f[82]*w3Ddx1+2.23606797749979*f[116]*dv3Ddx1; 
  out[139] += 3.464101615137754*f[94]*w2Ddx0+0.8944271909999161*f[49]*dv2Ddx0; 
  out[140] += 3.464101615137754*f[94]*w3Ddx1+f[143]*dv3Ddx1; 
  out[141] += 7.745966692414834*f[83]*w2Ddx0+2.23606797749979*f[115]*dv2Ddx0; 
  out[142] += 7.745966692414834*f[84]*w3Ddx1+(2.0*f[145]+2.23606797749979*f[48])*dv3Ddx1; 
  out[144] += 3.464101615137754*f[95]*w2Ddx0+f[146]*dv2Ddx0; 
  out[145] += 3.464101615137754*f[95]*w3Ddx1+0.8944271909999161*f[50]*dv3Ddx1; 
  out[147] += 3.464101615137755*f[106]*w3Ddx1+7.745966692414834*f[96]*w2Ddx0+(0.8944271909999159*f[166]+f[70])*dv3Ddx1+(2.0*f[149]+2.23606797749979*f[53])*dv2Ddx0; 
  out[148] += 7.745966692414834*f[96]*w3Ddx1+3.464101615137755*f[107]*w2Ddx0+(2.0*f[150]+2.23606797749979*f[52])*dv3Ddx1+(0.8944271909999159*f[163]+f[75])*dv2Ddx0; 
  out[149] += 3.464101615137755*f[108]*w3Ddx1+3.464101615137755*f[109]*w2Ddx0+(0.8944271909999159*f[168]+f[72])*dv3Ddx1+0.8944271909999159*f[55]*dv2Ddx0; 
  out[150] += 3.464101615137755*f[111]*w3Ddx1+3.464101615137755*f[112]*w2Ddx0+0.8944271909999159*f[54]*dv3Ddx1+(0.8944271909999159*f[169]+f[78])*dv2Ddx0; 
  out[151] += 3.464101615137755*f[115]*w3Ddx1+3.464101615137755*f[116]*w2Ddx0+(0.8944271909999159*f[180]+f[81])*dv3Ddx1+(0.8944271909999159*f[178]+f[84])*dv2Ddx0; 
  out[152] += 7.745966692414834*f[97]*w3Ddx1+7.745966692414834*f[98]*w2Ddx0+(2.0*f[155]+2.23606797749979*f[56])*dv3Ddx1+(2.0*f[154]+2.23606797749979*f[60])*dv2Ddx0; 
  out[153] += 3.464101615137754*f[121]*w3Ddx1+7.745966692414834*f[99]*w2Ddx0+(0.8944271909999159*f[184]+f[87])*dv3Ddx1+2.0*f[51]*dv2Ddx0; 
  out[154] += 7.745966692414834*f[99]*w3Ddx1+3.464101615137754*f[122]*w2Ddx0+(2.0*f[157]+2.23606797749979*f[58])*dv3Ddx1+0.8944271909999159*f[62]*dv2Ddx0; 
  out[155] += 3.464101615137754*f[125]*w3Ddx1+7.745966692414834*f[100]*w2Ddx0+0.8944271909999159*f[61]*dv3Ddx1+(2.0*f[157]+2.23606797749979*f[65])*dv2Ddx0; 
  out[156] += 7.745966692414834*f[100]*w3Ddx1+3.464101615137754*f[126]*w2Ddx0+2.0*f[51]*dv3Ddx1+(0.8944271909999159*f[185]+f[90])*dv2Ddx0; 
  out[157] += 3.464101615137754*f[127]*w3Ddx1+3.464101615137754*f[128]*w2Ddx0+0.8944271909999159*f[63]*dv3Ddx1+0.8944271909999159*f[67]*dv2Ddx0; 
  out[158] += 7.745966692414834*f[101]*w3Ddx1+7.745966692414834*f[102]*w2Ddx0+2.23606797749979*f[147]*dv3Ddx1+(2.0*f[160]+2.23606797749979*f[69])*dv2Ddx0; 
  out[159] += 3.464101615137754*f[130]*w3Ddx1+7.745966692414834*f[103]*w2Ddx0+f[162]*dv3Ddx1+2.0*f[52]*dv2Ddx0; 
  out[160] += 7.745966692414834*f[103]*w3Ddx1+3.464101615137754*f[131]*w2Ddx0+2.23606797749979*f[149]*dv3Ddx1+0.8944271909999159*f[71]*dv2Ddx0; 
  out[161] += 7.745966692414834*f[104]*w3Ddx1+7.745966692414834*f[105]*w2Ddx0+(2.0*f[164]+2.23606797749979*f[68])*dv3Ddx1+2.23606797749979*f[148]*dv2Ddx0; 
  out[162] += 7.745966692414834*f[108]*w2Ddx0+2.0*f[54]*dv2Ddx0; 
  out[163] += 7.745966692414834*f[109]*w3Ddx1+(2.0*f[169]+2.23606797749979*f[73])*dv3Ddx1; 
  out[164] += 3.464101615137754*f[132]*w3Ddx1+7.745966692414834*f[110]*w2Ddx0+0.8944271909999159*f[74]*dv3Ddx1+2.23606797749979*f[150]*dv2Ddx0; 
  out[165] += 7.745966692414834*f[110]*w3Ddx1+3.464101615137754*f[133]*w2Ddx0+2.0*f[53]*dv3Ddx1+f[167]*dv2Ddx0; 
  out[166] += 7.745966692414834*f[111]*w2Ddx0+(2.0*f[168]+2.23606797749979*f[77])*dv2Ddx0; 
  out[167] += 7.745966692414834*f[112]*w3Ddx1+2.0*f[55]*dv3Ddx1; 
  out[168] += 3.464101615137754*f[134]*w2Ddx0+0.8944271909999159*f[79]*dv2Ddx0; 
  out[169] += 3.464101615137754*f[134]*w3Ddx1+0.8944271909999159*f[76]*dv3Ddx1; 
  out[170] += 3.464101615137754*f[137]*w3Ddx1+7.745966692414834*f[113]*w2Ddx0+f[175]*dv3Ddx1+(2.0*f[172]+2.23606797749979*f[80])*dv2Ddx0; 
  out[171] += 7.745966692414834*f[113]*w3Ddx1+3.464101615137754*f[138]*w2Ddx0+2.23606797749979*f[151]*dv3Ddx1+(0.8944271909999159*f[188]+f[93])*dv2Ddx0; 
  out[172] += 3.464101615137754*f[139]*w3Ddx1+3.464101615137754*f[140]*w2Ddx0+f[177]*dv3Ddx1+0.8944271909999159*f[82]*dv2Ddx0; 
  out[173] += 3.464101615137754*f[141]*w3Ddx1+7.745966692414834*f[114]*w2Ddx0+(0.8944271909999159*f[189]+f[92])*dv3Ddx1+2.23606797749979*f[151]*dv2Ddx0; 
  out[174] += 7.745966692414834*f[114]*w3Ddx1+3.464101615137754*f[142]*w2Ddx0+(2.0*f[179]+2.23606797749979*f[80])*dv3Ddx1+f[176]*dv2Ddx0; 
  out[175] += 7.745966692414834*f[115]*w2Ddx0+(2.0*f[177]+2.23606797749979*f[83])*dv2Ddx0; 
  out[176] += 7.745966692414834*f[116]*w3Ddx1+(2.0*f[181]+2.23606797749979*f[82])*dv3Ddx1; 
  out[177] += 3.464101615137754*f[143]*w2Ddx0+0.8944271909999159*f[85]*dv2Ddx0; 
  out[178] += 3.464101615137754*f[143]*w3Ddx1+(0.8944271909999159*f[191]+f[94])*dv3Ddx1; 
  out[179] += 3.464101615137754*f[144]*w3Ddx1+3.464101615137754*f[145]*w2Ddx0+0.8944271909999159*f[83]*dv3Ddx1+f[181]*dv2Ddx0; 
  out[180] += 3.464101615137754*f[146]*w2Ddx0+(0.8944271909999159*f[191]+f[95])*dv2Ddx0; 
  out[181] += 3.464101615137754*f[146]*w3Ddx1+0.8944271909999159*f[85]*dv3Ddx1; 
  out[182] += 7.745966692414834*f[118]*w3Ddx1+7.745966692414834*f[119]*w2Ddx0+2.23606797749979*f[153]*dv3Ddx1+2.0*f[57]*dv2Ddx0; 
  out[183] += 7.745966692414834*f[123]*w3Ddx1+7.745966692414834*f[124]*w2Ddx0+2.0*f[59]*dv3Ddx1+2.23606797749979*f[156]*dv2Ddx0; 
  out[184] += 7.745966692414834*f[127]*w2Ddx0+2.0*f[66]*dv2Ddx0; 
  out[185] += 7.745966692414834*f[128]*w3Ddx1+2.0*f[64]*dv3Ddx1; 
  out[186] += 7.745966692414834*f[135]*w3Ddx1+7.745966692414834*f[136]*w2Ddx0+2.23606797749979*f[173]*dv3Ddx1+2.23606797749979*f[171]*dv2Ddx0; 
  out[187] += 7.745966692414834*f[139]*w2Ddx0+2.0*f[81]*dv2Ddx0; 
  out[188] += 7.745966692414834*f[140]*w3Ddx1+2.23606797749979*f[178]*dv3Ddx1; 
  out[189] += 7.745966692414834*f[144]*w2Ddx0+2.23606797749979*f[180]*dv2Ddx0; 
  out[190] += 7.745966692414834*f[145]*w3Ddx1+2.0*f[84]*dv3Ddx1; 
  out[192] += 7.745966692414834*f[147]*w3Ddx1+7.745966692414834*f[148]*w2Ddx0+(2.0*f[195]+2.23606797749979*f[101])*dv3Ddx1+(2.0*f[194]+2.23606797749979*f[105])*dv2Ddx0; 
  out[193] += 3.464101615137754*f[162]*w3Ddx1+7.745966692414834*f[149]*w2Ddx0+(0.8944271909999161*f[208]+f[130])*dv3Ddx1+2.0*f[96]*dv2Ddx0; 
  out[194] += 7.745966692414834*f[149]*w3Ddx1+3.464101615137754*f[163]*w2Ddx0+(2.0*f[197]+2.23606797749979*f[103])*dv3Ddx1+0.8944271909999161*f[107]*dv2Ddx0; 
  out[195] += 3.464101615137754*f[166]*w3Ddx1+7.745966692414834*f[150]*w2Ddx0+0.8944271909999161*f[106]*dv3Ddx1+(2.0*f[197]+2.23606797749979*f[110])*dv2Ddx0; 
  out[196] += 7.745966692414834*f[150]*w3Ddx1+3.464101615137754*f[167]*w2Ddx0+2.0*f[96]*dv3Ddx1+(0.8944271909999161*f[209]+f[133])*dv2Ddx0; 
  out[197] += 3.464101615137754*f[168]*w3Ddx1+3.464101615137754*f[169]*w2Ddx0+0.8944271909999161*f[108]*dv3Ddx1+0.8944271909999161*f[112]*dv2Ddx0; 
  out[198] += 3.464101615137754*f[175]*w3Ddx1+7.745966692414834*f[151]*w2Ddx0+(0.8944271909999161*f[218]+f[137])*dv3Ddx1+(2.0*f[200]+2.23606797749979*f[114])*dv2Ddx0; 
  out[199] += 7.745966692414834*f[151]*w3Ddx1+3.464101615137754*f[176]*w2Ddx0+(2.0*f[201]+2.23606797749979*f[113])*dv3Ddx1+(0.8944271909999161*f[215]+f[142])*dv2Ddx0; 
  out[200] += 3.464101615137754*f[177]*w3Ddx1+3.464101615137754*f[178]*w2Ddx0+(0.8944271909999161*f[220]+f[139])*dv3Ddx1+0.8944271909999161*f[116]*dv2Ddx0; 
  out[201] += 3.464101615137754*f[180]*w3Ddx1+3.464101615137754*f[181]*w2Ddx0+0.8944271909999161*f[115]*dv3Ddx1+(0.8944271909999161*f[221]+f[145])*dv2Ddx0; 
  out[202] += 7.745966692414834*f[153]*w3Ddx1+7.745966692414834*f[154]*w2Ddx0+(2.0*f[204]+2.23606797749979*f[118])*dv3Ddx1+2.0*f[98]*dv2Ddx0; 
  out[203] += 7.745966692414834*f[155]*w3Ddx1+7.745966692414834*f[156]*w2Ddx0+2.0*f[97]*dv3Ddx1+(2.0*f[205]+2.23606797749979*f[124])*dv2Ddx0; 
  out[204] += 3.464101615137755*f[184]*w3Ddx1+7.745966692414834*f[157]*w2Ddx0+0.8944271909999161*f[121]*dv3Ddx1+2.0*f[100]*dv2Ddx0; 
  out[205] += 7.745966692414834*f[157]*w3Ddx1+3.464101615137755*f[185]*w2Ddx0+2.0*f[99]*dv3Ddx1+0.8944271909999161*f[126]*dv2Ddx0; 
  out[206] += 7.745966692414834*f[159]*w3Ddx1+7.745966692414834*f[160]*w2Ddx0+2.23606797749979*f[193]*dv3Ddx1+2.0*f[102]*dv2Ddx0; 
  out[207] += 7.745966692414834*f[164]*w3Ddx1+7.745966692414834*f[165]*w2Ddx0+2.0*f[104]*dv3Ddx1+2.23606797749979*f[196]*dv2Ddx0; 
  out[208] += 7.745966692414834*f[168]*w2Ddx0+2.0*f[111]*dv2Ddx0; 
  out[209] += 7.745966692414834*f[169]*w3Ddx1+2.0*f[109]*dv3Ddx1; 
  out[210] += 7.745966692414834*f[170]*w3Ddx1+7.745966692414834*f[171]*w2Ddx0+2.23606797749979*f[198]*dv3Ddx1+(2.0*f[212]+2.23606797749979*f[136])*dv2Ddx0; 
  out[211] += 3.464101615137755*f[187]*w3Ddx1+7.745966692414834*f[172]*w2Ddx0+f[214]*dv3Ddx1+2.0*f[113]*dv2Ddx0; 
  out[212] += 7.745966692414834*f[172]*w3Ddx1+3.464101615137755*f[188]*w2Ddx0+2.23606797749979*f[200]*dv3Ddx1+0.8944271909999161*f[138]*dv2Ddx0; 
  out[213] += 7.745966692414834*f[173]*w3Ddx1+7.745966692414834*f[174]*w2Ddx0+(2.0*f[216]+2.23606797749979*f[135])*dv3Ddx1+2.23606797749979*f[199]*dv2Ddx0; 
  out[214] += 7.745966692414834*f[177]*w2Ddx0+2.0*f[115]*dv2Ddx0; 
  out[215] += 7.745966692414834*f[178]*w3Ddx1+(2.0*f[221]+2.23606797749979*f[140])*dv3Ddx1; 
  out[216] += 3.464101615137755*f[189]*w3Ddx1+7.745966692414834*f[179]*w2Ddx0+0.8944271909999161*f[141]*dv3Ddx1+2.23606797749979*f[201]*dv2Ddx0; 
  out[217] += 7.745966692414834*f[179]*w3Ddx1+3.464101615137755*f[190]*w2Ddx0+2.0*f[114]*dv3Ddx1+f[219]*dv2Ddx0; 
  out[218] += 7.745966692414834*f[180]*w2Ddx0+(2.0*f[220]+2.23606797749979*f[144])*dv2Ddx0; 
  out[219] += 7.745966692414834*f[181]*w3Ddx1+2.0*f[116]*dv3Ddx1; 
  out[220] += 3.464101615137755*f[191]*w2Ddx0+0.8944271909999161*f[146]*dv2Ddx0; 
  out[221] += 3.464101615137755*f[191]*w3Ddx1+0.8944271909999161*f[143]*dv3Ddx1; 
  out[222] += 7.745966692414834*f[193]*w3Ddx1+7.745966692414834*f[194]*w2Ddx0+(2.0*f[224]+2.23606797749979*f[159])*dv3Ddx1+2.0*f[148]*dv2Ddx0; 
  out[223] += 7.745966692414834*f[195]*w3Ddx1+7.745966692414834*f[196]*w2Ddx0+2.0*f[147]*dv3Ddx1+(2.0*f[225]+2.23606797749979*f[165])*dv2Ddx0; 
  out[224] += 3.464101615137755*f[208]*w3Ddx1+7.745966692414834*f[197]*w2Ddx0+0.8944271909999159*f[162]*dv3Ddx1+2.0*f[150]*dv2Ddx0; 
  out[225] += 7.745966692414834*f[197]*w3Ddx1+3.464101615137755*f[209]*w2Ddx0+2.0*f[149]*dv3Ddx1+0.8944271909999159*f[167]*dv2Ddx0; 
  out[226] += 7.745966692414834*f[198]*w3Ddx1+7.745966692414834*f[199]*w2Ddx0+(2.0*f[229]+2.23606797749979*f[170])*dv3Ddx1+(2.0*f[228]+2.23606797749979*f[174])*dv2Ddx0; 
  out[227] += 3.464101615137755*f[214]*w3Ddx1+7.745966692414834*f[200]*w2Ddx0+(0.8944271909999159*f[235]+f[187])*dv3Ddx1+2.0*f[151]*dv2Ddx0; 
  out[228] += 7.745966692414834*f[200]*w3Ddx1+3.464101615137755*f[215]*w2Ddx0+(2.0*f[231]+2.23606797749979*f[172])*dv3Ddx1+0.8944271909999159*f[176]*dv2Ddx0; 
  out[229] += 3.464101615137755*f[218]*w3Ddx1+7.745966692414834*f[201]*w2Ddx0+0.8944271909999159*f[175]*dv3Ddx1+(2.0*f[231]+2.23606797749979*f[179])*dv2Ddx0; 
  out[230] += 7.745966692414834*f[201]*w3Ddx1+3.464101615137755*f[219]*w2Ddx0+2.0*f[151]*dv3Ddx1+(0.8944271909999159*f[236]+f[190])*dv2Ddx0; 
  out[231] += 3.464101615137755*f[220]*w3Ddx1+3.464101615137755*f[221]*w2Ddx0+0.8944271909999159*f[177]*dv3Ddx1+0.8944271909999159*f[181]*dv2Ddx0; 
  out[232] += 7.745966692414834*f[204]*w3Ddx1+7.745966692414834*f[205]*w2Ddx0+2.0*f[153]*dv3Ddx1+2.0*f[156]*dv2Ddx0; 
  out[233] += 7.745966692414834*f[211]*w3Ddx1+7.745966692414834*f[212]*w2Ddx0+2.23606797749979*f[227]*dv3Ddx1+2.0*f[171]*dv2Ddx0; 
  out[234] += 7.745966692414834*f[216]*w3Ddx1+7.745966692414834*f[217]*w2Ddx0+2.0*f[173]*dv3Ddx1+2.23606797749979*f[230]*dv2Ddx0; 
  out[235] += 7.745966692414834*f[220]*w2Ddx0+2.0*f[180]*dv2Ddx0; 
  out[236] += 7.745966692414834*f[221]*w3Ddx1+2.0*f[178]*dv3Ddx1; 
  out[237] += 7.745966692414834*f[224]*w3Ddx1+7.745966692414834*f[225]*w2Ddx0+2.0*f[193]*dv3Ddx1+2.0*f[196]*dv2Ddx0; 
  out[238] += 7.745966692414834*f[227]*w3Ddx1+7.745966692414834*f[228]*w2Ddx0+(2.0*f[240]+2.23606797749979*f[211])*dv3Ddx1+2.0*f[199]*dv2Ddx0; 
  out[239] += 7.745966692414834*f[229]*w3Ddx1+7.745966692414834*f[230]*w2Ddx0+2.0*f[198]*dv3Ddx1+(2.0*f[241]+2.23606797749979*f[217])*dv2Ddx0; 
  out[240] += 3.464101615137754*f[235]*w3Ddx1+7.745966692414834*f[231]*w2Ddx0+0.8944271909999161*f[214]*dv3Ddx1+2.0*f[201]*dv2Ddx0; 
  out[241] += 7.745966692414834*f[231]*w3Ddx1+3.464101615137754*f[236]*w2Ddx0+2.0*f[200]*dv3Ddx1+0.8944271909999161*f[219]*dv2Ddx0; 
  out[242] += 7.745966692414834*f[240]*w3Ddx1+7.745966692414834*f[241]*w2Ddx0+2.0*f[227]*dv3Ddx1+2.0*f[230]*dv2Ddx0; 

  return 5.0*(fabs(w2Ddx0)+0.5*dv2Ddx0+fabs(w3Ddx1)+0.5*dv3Ddx1);
} 