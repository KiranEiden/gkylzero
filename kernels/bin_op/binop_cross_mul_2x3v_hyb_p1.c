// Thu Jul 28 13:07:03 2022
#include <gkyl_binop_cross_mul_hyb.h>
GKYL_CU_DH
void
binop_cross_mul_2x3v_hyb_p1(const double *f, const double *g, double *fg )
{
  double tmp[80] = {0.};
  tmp[0] =  5.0000000000000000e-01*f[3]*g[6]+5.0000000000000000e-01*g[0]*f[0]+5.0000000000000000e-01*g[2]*f[2]+5.0000000000000000e-01*g[1]*f[1];
  tmp[1] =  5.0000000000000000e-01*f[2]*g[6]+5.0000000000000000e-01*g[2]*f[3]+5.0000000000000000e-01*g[0]*f[1]+5.0000000000000000e-01*g[1]*f[0];
  tmp[2] =  5.0000000000000000e-01*g[2]*f[0]+5.0000000000000000e-01*g[1]*f[3]+5.0000000000000000e-01*g[0]*f[2]+5.0000000000000000e-01*f[1]*g[6];
  tmp[3] =  5.0000000000000000e-01*g[8]*f[2]+5.0000000000000000e-01*g[16]*f[3]+5.0000000000000000e-01*g[3]*f[0]+5.0000000000000000e-01*f[1]*g[7];
  tmp[4] =  5.0000000000000000e-01*g[9]*f[1]+5.0000000000000000e-01*g[10]*f[2]+5.0000000000000000e-01*f[3]*g[17]+5.0000000000000000e-01*f[0]*g[4];
  tmp[5] =  5.0000000000000000e-01*g[20]*f[3]+5.0000000000000000e-01*f[1]*g[12]+5.0000000000000000e-01*g[13]*f[2]+5.0000000000000000e-01*f[0]*g[5];
  tmp[6] =  5.0000000000000000e-01*g[0]*f[3]+5.0000000000000000e-01*g[1]*f[2]+5.0000000000000000e-01*g[2]*f[1]+5.0000000000000000e-01*f[0]*g[6];
  tmp[7] =  5.0000000000000000e-01*g[16]*f[2]+5.0000000000000000e-01*g[8]*f[3]+5.0000000000000000e-01*g[3]*f[1]+5.0000000000000000e-01*f[0]*g[7];
  tmp[8] =  5.0000000000000000e-01*g[3]*f[2]+5.0000000000000000e-01*g[8]*f[0]+5.0000000000000000e-01*f[3]*g[7]+5.0000000000000000e-01*g[16]*f[1];
  tmp[9] =  5.0000000000000000e-01*g[9]*f[0]+5.0000000000000000e-01*g[17]*f[2]+5.0000000000000000e-01*f[1]*g[4]+5.0000000000000000e-01*g[10]*f[3];
  tmp[10] =  5.0000000000000000e-01*f[0]*g[10]+5.0000000000000000e-01*g[4]*f[2]+5.0000000000000000e-01*f[1]*g[17]+5.0000000000000000e-01*g[9]*f[3];
  tmp[11] =  5.0000000000000000e-01*g[11]*f[0]+5.0000000000000000e-01*f[1]*g[18]+5.0000000000000000e-01*g[19]*f[2]+5.0000000000000000e-01*g[26]*f[3];
  tmp[12] =  5.0000000000000000e-01*g[12]*f[0]+5.0000000000000000e-01*g[20]*f[2]+5.0000000000000000e-01*g[13]*f[3]+5.0000000000000000e-01*f[1]*g[5];
  tmp[13] =  5.0000000000000000e-01*g[12]*f[3]+5.0000000000000000e-01*g[20]*f[1]+5.0000000000000000e-01*g[5]*f[2]+5.0000000000000000e-01*g[13]*f[0];
  tmp[14] =  5.0000000000000000e-01*f[0]*g[14]+5.0000000000000000e-01*g[27]*f[3]+5.0000000000000000e-01*g[22]*f[2]+5.0000000000000000e-01*f[1]*g[21];
  tmp[15] =  5.0000000000000000e-01*f[1]*g[23]+5.0000000000000000e-01*g[28]*f[3]+5.0000000000000000e-01*g[24]*f[2]+5.0000000000000000e-01*g[15]*f[0];
  tmp[16] =  5.0000000000000000e-01*g[3]*f[3]+5.0000000000000000e-01*g[7]*f[2]+5.0000000000000000e-01*g[8]*f[1]+5.0000000000000000e-01*g[16]*f[0];
  tmp[17] =  5.0000000000000000e-01*f[1]*g[10]+5.0000000000000000e-01*f[0]*g[17]+5.0000000000000000e-01*g[9]*f[2]+5.0000000000000000e-01*g[4]*f[3];
  tmp[18] =  5.0000000000000000e-01*g[19]*f[3]+5.0000000000000000e-01*g[11]*f[1]+5.0000000000000000e-01*f[0]*g[18]+5.0000000000000000e-01*g[26]*f[2];
  tmp[19] =  5.0000000000000000e-01*g[18]*f[3]+5.0000000000000000e-01*f[1]*g[26]+5.0000000000000000e-01*g[19]*f[0]+5.0000000000000000e-01*g[11]*f[2];
  tmp[20] =  5.0000000000000000e-01*g[5]*f[3]+5.0000000000000000e-01*g[20]*f[0]+5.0000000000000000e-01*g[12]*f[2]+5.0000000000000000e-01*g[13]*f[1];
  tmp[21] =  5.0000000000000000e-01*f[1]*g[14]+5.0000000000000000e-01*f[0]*g[21]+5.0000000000000000e-01*g[22]*f[3]+5.0000000000000000e-01*g[27]*f[2];
  tmp[22] =  5.0000000000000000e-01*f[1]*g[27]+5.0000000000000000e-01*g[14]*f[2]+5.0000000000000000e-01*f[3]*g[21]+5.0000000000000000e-01*f[0]*g[22];
  tmp[23] =  5.0000000000000000e-01*g[23]*f[0]+5.0000000000000000e-01*g[24]*f[3]+5.0000000000000000e-01*g[28]*f[2]+5.0000000000000000e-01*g[15]*f[1];
  tmp[24] =  5.0000000000000000e-01*g[15]*f[2]+5.0000000000000000e-01*g[28]*f[1]+5.0000000000000000e-01*g[24]*f[0]+5.0000000000000000e-01*g[23]*f[3];
  tmp[25] =  5.0000000000000000e-01*g[29]*f[1]+5.0000000000000000e-01*g[25]*f[0]+5.0000000000000000e-01*g[31]*f[3]+5.0000000000000000e-01*g[30]*f[2];
  tmp[26] =  5.0000000000000000e-01*g[18]*f[2]+5.0000000000000000e-01*g[11]*f[3]+5.0000000000000000e-01*g[19]*f[1]+5.0000000000000000e-01*f[0]*g[26];
  tmp[27] =  5.0000000000000000e-01*g[27]*f[0]+5.0000000000000000e-01*g[21]*f[2]+5.0000000000000000e-01*g[14]*f[3]+5.0000000000000000e-01*f[1]*g[22];
  tmp[28] =  5.0000000000000000e-01*g[28]*f[0]+5.0000000000000000e-01*g[15]*f[3]+5.0000000000000000e-01*g[24]*f[1]+5.0000000000000000e-01*g[23]*f[2];
  tmp[29] =  5.0000000000000000e-01*g[29]*f[0]+5.0000000000000000e-01*g[31]*f[2]+5.0000000000000000e-01*g[25]*f[1]+5.0000000000000000e-01*g[30]*f[3];
  tmp[30] =  5.0000000000000000e-01*f[1]*g[31]+5.0000000000000000e-01*g[29]*f[3]+5.0000000000000000e-01*g[25]*f[2]+5.0000000000000000e-01*f[0]*g[30];
  tmp[31] =  5.0000000000000000e-01*f[0]*g[31]+5.0000000000000000e-01*f[1]*g[30]+5.0000000000000000e-01*g[25]*f[3]+5.0000000000000000e-01*g[29]*f[2];
  tmp[32] =  5.0000000000000000e-01*g[34]*f[2]+5.0000000000000000e-01*g[33]*f[1]+5.0000000000000000e-01*f[0]*g[32]+5.0000000000000000e-01*g[37]*f[3];
  tmp[33] =  5.0000000000000000e-01*g[34]*f[3]+5.0000000000000000e-01*g[33]*f[0]+5.0000000000000000e-01*f[1]*g[32]+5.0000000000000000e-01*g[37]*f[2];
  tmp[34] =  5.0000000000000000e-01*f[0]*g[34]+5.0000000000000000e-01*f[1]*g[37]+5.0000000000000000e-01*g[32]*f[2]+5.0000000000000000e-01*g[33]*f[3];
  tmp[35] =  5.0000000000000000e-01*g[35]*f[0]+5.0000000000000000e-01*g[39]*f[2]+5.0000000000000000e-01*g[43]*f[3]+5.0000000000000000e-01*g[38]*f[1];
  tmp[36] =  5.0000000000000000e-01*g[44]*f[3]+5.0000000000000000e-01*g[41]*f[2]+5.0000000000000000e-01*g[36]*f[0]+5.0000000000000000e-01*f[1]*g[40];
  tmp[37] =  5.0000000000000000e-01*f[1]*g[34]+5.0000000000000000e-01*g[32]*f[3]+5.0000000000000000e-01*g[33]*f[2]+5.0000000000000000e-01*f[0]*g[37];
  tmp[38] =  5.0000000000000000e-01*g[35]*f[1]+5.0000000000000000e-01*g[39]*f[3]+5.0000000000000000e-01*g[38]*f[0]+5.0000000000000000e-01*g[43]*f[2];
  tmp[39] =  5.0000000000000000e-01*g[43]*f[1]+5.0000000000000000e-01*g[38]*f[3]+5.0000000000000000e-01*f[0]*g[39]+5.0000000000000000e-01*g[35]*f[2];
  tmp[40] =  5.0000000000000000e-01*g[44]*f[2]+5.0000000000000000e-01*f[0]*g[40]+5.0000000000000000e-01*g[36]*f[1]+5.0000000000000000e-01*g[41]*f[3];
  tmp[41] =  5.0000000000000000e-01*g[36]*f[2]+5.0000000000000000e-01*g[41]*f[0]+5.0000000000000000e-01*f[3]*g[40]+5.0000000000000000e-01*g[44]*f[1];
  tmp[42] =  5.0000000000000000e-01*f[0]*g[42]+5.0000000000000000e-01*g[46]*f[2]+5.0000000000000000e-01*f[1]*g[45]+5.0000000000000000e-01*g[47]*f[3];
  tmp[43] =  5.0000000000000000e-01*g[38]*f[2]+5.0000000000000000e-01*f[1]*g[39]+5.0000000000000000e-01*g[43]*f[0]+5.0000000000000000e-01*g[35]*f[3];
  tmp[44] =  5.0000000000000000e-01*f[2]*g[40]+5.0000000000000000e-01*g[41]*f[1]+5.0000000000000000e-01*g[36]*f[3]+5.0000000000000000e-01*g[44]*f[0];
  tmp[45] =  5.0000000000000000e-01*g[46]*f[3]+5.0000000000000000e-01*f[1]*g[42]+5.0000000000000000e-01*f[0]*g[45]+5.0000000000000000e-01*g[47]*f[2];
  tmp[46] =  5.0000000000000000e-01*f[1]*g[47]+5.0000000000000000e-01*g[42]*f[2]+5.0000000000000000e-01*g[46]*f[0]+5.0000000000000000e-01*g[45]*f[3];
  tmp[47] =  5.0000000000000000e-01*g[47]*f[0]+5.0000000000000000e-01*g[42]*f[3]+5.0000000000000000e-01*g[46]*f[1]+5.0000000000000000e-01*g[45]*f[2];
  tmp[48] =  5.0000000000000000e-01*g[49]*f[1]+5.0000000000000000e-01*g[50]*f[2]+5.0000000000000000e-01*f[0]*g[48]+5.0000000000000000e-01*f[3]*g[53];
  tmp[49] =  5.0000000000000000e-01*g[49]*f[0]+5.0000000000000000e-01*g[50]*f[3]+5.0000000000000000e-01*f[2]*g[53]+5.0000000000000000e-01*f[1]*g[48];
  tmp[50] =  5.0000000000000000e-01*f[0]*g[50]+5.0000000000000000e-01*g[48]*f[2]+5.0000000000000000e-01*g[49]*f[3]+5.0000000000000000e-01*f[1]*g[53];
  tmp[51] =  5.0000000000000000e-01*g[55]*f[2]+5.0000000000000000e-01*g[54]*f[1]+5.0000000000000000e-01*g[59]*f[3]+5.0000000000000000e-01*g[51]*f[0];
  tmp[52] =  5.0000000000000000e-01*g[57]*f[2]+5.0000000000000000e-01*g[52]*f[0]+5.0000000000000000e-01*g[56]*f[1]+5.0000000000000000e-01*g[60]*f[3];
  tmp[53] =  5.0000000000000000e-01*f[1]*g[50]+5.0000000000000000e-01*f[3]*g[48]+5.0000000000000000e-01*f[0]*g[53]+5.0000000000000000e-01*g[49]*f[2];
  tmp[54] =  5.0000000000000000e-01*g[59]*f[2]+5.0000000000000000e-01*g[54]*f[0]+5.0000000000000000e-01*g[55]*f[3]+5.0000000000000000e-01*g[51]*f[1];
  tmp[55] =  5.0000000000000000e-01*g[51]*f[2]+5.0000000000000000e-01*f[0]*g[55]+5.0000000000000000e-01*g[59]*f[1]+5.0000000000000000e-01*g[54]*f[3];
  tmp[56] =  5.0000000000000000e-01*f[1]*g[52]+5.0000000000000000e-01*g[57]*f[3]+5.0000000000000000e-01*g[56]*f[0]+5.0000000000000000e-01*g[60]*f[2];
  tmp[57] =  5.0000000000000000e-01*f[1]*g[60]+5.0000000000000000e-01*g[52]*f[2]+5.0000000000000000e-01*g[56]*f[3]+5.0000000000000000e-01*g[57]*f[0];
  tmp[58] =  5.0000000000000000e-01*g[63]*f[3]+5.0000000000000000e-01*f[1]*g[61]+5.0000000000000000e-01*f[0]*g[58]+5.0000000000000000e-01*g[62]*f[2];
  tmp[59] =  5.0000000000000000e-01*f[1]*g[55]+5.0000000000000000e-01*g[51]*f[3]+5.0000000000000000e-01*g[54]*f[2]+5.0000000000000000e-01*g[59]*f[0];
  tmp[60] =  5.0000000000000000e-01*g[60]*f[0]+5.0000000000000000e-01*g[56]*f[2]+5.0000000000000000e-01*g[52]*f[3]+5.0000000000000000e-01*g[57]*f[1];
  tmp[61] =  5.0000000000000000e-01*g[63]*f[2]+5.0000000000000000e-01*f[0]*g[61]+5.0000000000000000e-01*g[62]*f[3]+5.0000000000000000e-01*f[1]*g[58];
  tmp[62] =  5.0000000000000000e-01*g[62]*f[0]+5.0000000000000000e-01*f[1]*g[63]+5.0000000000000000e-01*g[58]*f[2]+5.0000000000000000e-01*f[3]*g[61];
  tmp[63] =  5.0000000000000000e-01*g[62]*f[1]+5.0000000000000000e-01*f[0]*g[63]+5.0000000000000000e-01*g[61]*f[2]+5.0000000000000000e-01*g[58]*f[3];
  tmp[64] =  5.0000000000000000e-01*g[65]*f[1]+5.0000000000000000e-01*f[0]*g[64]+5.0000000000000000e-01*g[66]*f[2]+5.0000000000000000e-01*g[69]*f[3];
  tmp[65] =  5.0000000000000000e-01*f[1]*g[64]+5.0000000000000000e-01*g[65]*f[0]+5.0000000000000000e-01*g[69]*f[2]+5.0000000000000000e-01*f[3]*g[66];
  tmp[66] =  5.0000000000000000e-01*g[69]*f[1]+5.0000000000000000e-01*g[64]*f[2]+5.0000000000000000e-01*f[0]*g[66]+5.0000000000000000e-01*g[65]*f[3];
  tmp[67] =  5.0000000000000000e-01*g[71]*f[2]+5.0000000000000000e-01*f[0]*g[67]+5.0000000000000000e-01*g[75]*f[3]+5.0000000000000000e-01*f[1]*g[70];
  tmp[68] =  5.0000000000000000e-01*g[72]*f[1]+5.0000000000000000e-01*g[68]*f[0]+5.0000000000000000e-01*g[76]*f[3]+5.0000000000000000e-01*g[73]*f[2];
  tmp[69] =  5.0000000000000000e-01*g[69]*f[0]+5.0000000000000000e-01*g[65]*f[2]+5.0000000000000000e-01*f[1]*g[66]+5.0000000000000000e-01*g[64]*f[3];
  tmp[70] =  5.0000000000000000e-01*g[75]*f[2]+5.0000000000000000e-01*f[0]*g[70]+5.0000000000000000e-01*f[1]*g[67]+5.0000000000000000e-01*g[71]*f[3];
  tmp[71] =  5.0000000000000000e-01*g[70]*f[3]+5.0000000000000000e-01*g[67]*f[2]+5.0000000000000000e-01*g[71]*f[0]+5.0000000000000000e-01*g[75]*f[1];
  tmp[72] =  5.0000000000000000e-01*g[68]*f[1]+5.0000000000000000e-01*g[72]*f[0]+5.0000000000000000e-01*g[76]*f[2]+5.0000000000000000e-01*g[73]*f[3];
  tmp[73] =  5.0000000000000000e-01*f[1]*g[76]+5.0000000000000000e-01*g[72]*f[3]+5.0000000000000000e-01*g[68]*f[2]+5.0000000000000000e-01*f[0]*g[73];
  tmp[74] =  5.0000000000000000e-01*f[1]*g[77]+5.0000000000000000e-01*g[78]*f[2]+5.0000000000000000e-01*f[3]*g[79]+5.0000000000000000e-01*g[74]*f[0];
  tmp[75] =  5.0000000000000000e-01*g[67]*f[3]+5.0000000000000000e-01*g[70]*f[2]+5.0000000000000000e-01*g[71]*f[1]+5.0000000000000000e-01*g[75]*f[0];
  tmp[76] =  5.0000000000000000e-01*f[1]*g[73]+5.0000000000000000e-01*g[68]*f[3]+5.0000000000000000e-01*f[0]*g[76]+5.0000000000000000e-01*g[72]*f[2];
  tmp[77] =  5.0000000000000000e-01*g[78]*f[3]+5.0000000000000000e-01*f[2]*g[79]+5.0000000000000000e-01*f[0]*g[77]+5.0000000000000000e-01*g[74]*f[1];
  tmp[78] =  5.0000000000000000e-01*g[74]*f[2]+5.0000000000000000e-01*g[78]*f[0]+5.0000000000000000e-01*g[77]*f[3]+5.0000000000000000e-01*f[1]*g[79];
  tmp[79] =  5.0000000000000000e-01*g[74]*f[3]+5.0000000000000000e-01*g[77]*f[2]+5.0000000000000000e-01*f[0]*g[79]+5.0000000000000000e-01*g[78]*f[1];
 
  fg[0] = tmp[0];
  fg[1] = tmp[1];
  fg[2] = tmp[2];
  fg[3] = tmp[3];
  fg[4] = tmp[4];
  fg[5] = tmp[5];
  fg[6] = tmp[6];
  fg[7] = tmp[7];
  fg[8] = tmp[8];
  fg[9] = tmp[9];
  fg[10] = tmp[10];
  fg[11] = tmp[11];
  fg[12] = tmp[12];
  fg[13] = tmp[13];
  fg[14] = tmp[14];
  fg[15] = tmp[15];
  fg[16] = tmp[16];
  fg[17] = tmp[17];
  fg[18] = tmp[18];
  fg[19] = tmp[19];
  fg[20] = tmp[20];
  fg[21] = tmp[21];
  fg[22] = tmp[22];
  fg[23] = tmp[23];
  fg[24] = tmp[24];
  fg[25] = tmp[25];
  fg[26] = tmp[26];
  fg[27] = tmp[27];
  fg[28] = tmp[28];
  fg[29] = tmp[29];
  fg[30] = tmp[30];
  fg[31] = tmp[31];
  fg[32] = tmp[32];
  fg[33] = tmp[33];
  fg[34] = tmp[34];
  fg[35] = tmp[35];
  fg[36] = tmp[36];
  fg[37] = tmp[37];
  fg[38] = tmp[38];
  fg[39] = tmp[39];
  fg[40] = tmp[40];
  fg[41] = tmp[41];
  fg[42] = tmp[42];
  fg[43] = tmp[43];
  fg[44] = tmp[44];
  fg[45] = tmp[45];
  fg[46] = tmp[46];
  fg[47] = tmp[47];
  fg[48] = tmp[48];
  fg[49] = tmp[49];
  fg[50] = tmp[50];
  fg[51] = tmp[51];
  fg[52] = tmp[52];
  fg[53] = tmp[53];
  fg[54] = tmp[54];
  fg[55] = tmp[55];
  fg[56] = tmp[56];
  fg[57] = tmp[57];
  fg[58] = tmp[58];
  fg[59] = tmp[59];
  fg[60] = tmp[60];
  fg[61] = tmp[61];
  fg[62] = tmp[62];
  fg[63] = tmp[63];
  fg[64] = tmp[64];
  fg[65] = tmp[65];
  fg[66] = tmp[66];
  fg[67] = tmp[67];
  fg[68] = tmp[68];
  fg[69] = tmp[69];
  fg[70] = tmp[70];
  fg[71] = tmp[71];
  fg[72] = tmp[72];
  fg[73] = tmp[73];
  fg[74] = tmp[74];
  fg[75] = tmp[75];
  fg[76] = tmp[76];
  fg[77] = tmp[77];
  fg[78] = tmp[78];
  fg[79] = tmp[79];
}

