#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_stream_vol_2x3v_ser_p2(const double *w, const double *dxv, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // f: Input distribution function.
  // out: Incremented output.
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
  out[57] += 7.745966692414834*f[21]*w3Ddx1+3.464101615137755*f[34]*w2Ddx0+2.23606797749979*f[51]*dv3Ddx1+f[17]*dv2Ddx0; 
  out[58] += 3.464101615137755*f[35]*w3Ddx1+3.464101615137755*f[36]*w2Ddx0+f[63]*dv3Ddx1+0.8944271909999159*f[8]*dv2Ddx0; 
  out[59] += 3.464101615137755*f[37]*w3Ddx1+7.745966692414834*f[22]*w2Ddx0+f[16]*dv3Ddx1+2.23606797749979*f[51]*dv2Ddx0; 
  out[60] += 7.745966692414834*f[22]*w3Ddx1+3.464101615137755*f[38]*w2Ddx0+(2.0*f[65]+2.23606797749979*f[6])*dv3Ddx1+f[62]*dv2Ddx0; 
  out[61] += 7.745966692414834*f[23]*w2Ddx0+(2.0*f[63]+2.23606797749979*f[9])*dv2Ddx0; 
  out[62] += 7.745966692414834*f[24]*w3Ddx1+(2.0*f[67]+2.23606797749979*f[8])*dv3Ddx1; 
  out[63] += 3.464101615137755*f[39]*w2Ddx0+0.8944271909999159*f[11]*dv2Ddx0; 
  out[64] += 3.464101615137755*f[39]*w3Ddx1+f[18]*dv3Ddx1; 
  out[65] += 3.464101615137755*f[40]*w3Ddx1+3.464101615137755*f[41]*w2Ddx0+0.8944271909999159*f[9]*dv3Ddx1+f[67]*dv2Ddx0; 
  out[66] += 3.464101615137755*f[42]*w2Ddx0+f[19]*dv2Ddx0; 
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
  out[81] += 3.464101615137755*f[49]*w2Ddx0+f[20]*dv2Ddx0; 
  out[82] += 3.464101615137755*f[49]*w3Ddx1+f[85]*dv3Ddx1; 
  out[83] += 3.464101615137755*f[50]*w2Ddx0+f[85]*dv2Ddx0; 
  out[84] += 3.464101615137755*f[50]*w3Ddx1+f[20]*dv3Ddx1; 
  out[86] += 3.464101615137754*f[54]*w3Ddx1+3.464101615137754*f[55]*w2Ddx0+(0.8944271909999161*f[101]+f[26])*dv3Ddx1+(0.8944271909999161*f[99]+f[29])*dv2Ddx0; 
  out[87] += 3.464101615137755*f[61]*w3Ddx1+7.745966692414834*f[51]*w2Ddx0+f[33]*dv3Ddx1+(2.0*f[89]+2.23606797749979*f[22])*dv2Ddx0; 
  out[88] += 7.745966692414834*f[51]*w3Ddx1+3.464101615137755*f[62]*w2Ddx0+(2.0*f[90]+2.23606797749979*f[21])*dv3Ddx1+f[38]*dv2Ddx0; 
  out[89] += 3.464101615137755*f[63]*w3Ddx1+3.464101615137755*f[64]*w2Ddx0+f[35]*dv3Ddx1+0.8944271909999161*f[24]*dv2Ddx0; 
  out[90] += 3.464101615137755*f[66]*w3Ddx1+3.464101615137755*f[67]*w2Ddx0+0.8944271909999161*f[23]*dv3Ddx1+f[41]*dv2Ddx0; 
  out[91] += 3.464101615137755*f[70]*w3Ddx1+7.745966692414834*f[52]*w2Ddx0+f[96]*dv3Ddx1+(2.0*f[93]+2.23606797749979*f[25])*dv2Ddx0; 
  out[92] += 7.745966692414834*f[52]*w3Ddx1+3.464101615137755*f[71]*w2Ddx0+2.23606797749979*f[86]*dv3Ddx1+f[44]*dv2Ddx0; 
  out[93] += 3.464101615137755*f[72]*w3Ddx1+3.464101615137755*f[73]*w2Ddx0+f[98]*dv3Ddx1+0.8944271909999161*f[27]*dv2Ddx0; 
  out[94] += 3.464101615137755*f[74]*w3Ddx1+7.745966692414834*f[53]*w2Ddx0+f[43]*dv3Ddx1+2.23606797749979*f[86]*dv2Ddx0; 
  out[95] += 7.745966692414834*f[53]*w3Ddx1+3.464101615137755*f[75]*w2Ddx0+(2.0*f[100]+2.23606797749979*f[25])*dv3Ddx1+f[97]*dv2Ddx0; 
  out[96] += 7.745966692414834*f[54]*w2Ddx0+(2.0*f[98]+2.23606797749979*f[28])*dv2Ddx0; 
  out[97] += 7.745966692414834*f[55]*w3Ddx1+(2.0*f[102]+2.23606797749979*f[27])*dv3Ddx1; 
  out[98] += 3.464101615137755*f[76]*w2Ddx0+0.8944271909999161*f[30]*dv2Ddx0; 
  out[99] += 3.464101615137755*f[76]*w3Ddx1+f[45]*dv3Ddx1; 
  out[100] += 3.464101615137755*f[77]*w3Ddx1+3.464101615137755*f[78]*w2Ddx0+0.8944271909999161*f[28]*dv3Ddx1+f[102]*dv2Ddx0; 
  out[101] += 3.464101615137755*f[79]*w2Ddx0+f[46]*dv2Ddx0; 
  out[102] += 3.464101615137755*f[79]*w3Ddx1+0.8944271909999161*f[30]*dv3Ddx1; 
  out[103] += 3.464101615137755*f[81]*w3Ddx1+3.464101615137755*f[82]*w2Ddx0+f[105]*dv3Ddx1+f[48]*dv2Ddx0; 
  out[104] += 3.464101615137755*f[83]*w3Ddx1+3.464101615137755*f[84]*w2Ddx0+f[47]*dv3Ddx1+f[106]*dv2Ddx0; 
  out[105] += 3.464101615137755*f[85]*w2Ddx0+f[50]*dv2Ddx0; 
  out[106] += 3.464101615137755*f[85]*w3Ddx1+f[49]*dv3Ddx1; 
  out[107] += 3.464101615137755*f[96]*w3Ddx1+7.745966692414834*f[86]*w2Ddx0+f[70]*dv3Ddx1+(2.0*f[109]+2.23606797749979*f[53])*dv2Ddx0; 
  out[108] += 7.745966692414834*f[86]*w3Ddx1+3.464101615137755*f[97]*w2Ddx0+(2.0*f[110]+2.23606797749979*f[52])*dv3Ddx1+f[75]*dv2Ddx0; 
  out[109] += 3.464101615137755*f[98]*w3Ddx1+3.464101615137755*f[99]*w2Ddx0+f[72]*dv3Ddx1+0.8944271909999159*f[55]*dv2Ddx0; 
  out[110] += 3.464101615137755*f[101]*w3Ddx1+3.464101615137755*f[102]*w2Ddx0+0.8944271909999159*f[54]*dv3Ddx1+f[78]*dv2Ddx0; 
  out[111] += 3.464101615137755*f[105]*w3Ddx1+3.464101615137755*f[106]*w2Ddx0+f[81]*dv3Ddx1+f[84]*dv2Ddx0; 

  return fabs(w2Ddx0)+0.5*dv2Ddx0+fabs(w3Ddx1)+0.5*dv3Ddx1;
} 
