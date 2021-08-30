// Thu Aug 26 15:51:57 2021
#include <gkyl_binop_cross_mul_ser.h>
GKYL_CU_DH
void
binop_cross_mul_1d_4d_ser_p3(const double *f, const double *g, double *fg )
{
  fg[0] =  7.0710678118654757e-01*f[3]*g[31]+7.0710678118654757e-01*f[1]*g[1]+7.0710678118654757e-01*f[2]*g[11]+7.0710678118654757e-01*f[0]*g[0];
  fg[1] =  6.2105900340811881e-01*f[3]*g[11]+6.3245553203367588e-01*f[2]*g[1]+6.3245553203367588e-01*f[1]*g[11]+7.0710678118654757e-01*f[1]*g[0]+6.2105900340811881e-01*f[2]*g[31]+7.0710678118654757e-01*f[0]*g[1];
  fg[2] =  7.0710678118654757e-01*g[5]*f[1]+7.0710678118654757e-01*g[48]*f[3]+7.0710678118654757e-01*f[0]*g[2]+7.0710678118654757e-01*f[2]*g[19];
  fg[3] =  7.0710678118654757e-01*f[3]*g[50]+7.0710678118654757e-01*f[0]*g[3]+7.0710678118654757e-01*g[21]*f[2]+7.0710678118654757e-01*g[6]*f[1];
  fg[4] =  7.0710678118654757e-01*g[4]*f[0]+7.0710678118654757e-01*f[2]*g[25]+7.0710678118654757e-01*f[3]*g[54]+7.0710678118654757e-01*f[1]*g[8];
  fg[5] =  7.0710678118654757e-01*g[5]*f[0]+7.0710678118654757e-01*f[1]*g[2]+6.2105900340811881e-01*g[48]*f[2]+6.3245553203367588e-01*f[1]*g[19]+6.2105900340811881e-01*f[3]*g[19]+6.3245553203367588e-01*g[5]*f[2];
  fg[6] =  6.2105900340811881e-01*g[21]*f[3]+6.3245553203367588e-01*g[6]*f[2]+7.0710678118654757e-01*f[1]*g[3]+6.3245553203367588e-01*g[21]*f[1]+6.2105900340811881e-01*f[2]*g[50]+7.0710678118654757e-01*g[6]*f[0];
  fg[7] =  7.0710678118654757e-01*f[1]*g[15]+7.0710678118654757e-01*f[2]*g[36]+7.0710678118654757e-01*g[64]*f[3]+7.0710678118654757e-01*f[0]*g[7];
  fg[8] =  6.2105900340811881e-01*f[2]*g[54]+6.2105900340811881e-01*f[3]*g[25]+6.3245553203367588e-01*f[1]*g[25]+7.0710678118654757e-01*g[4]*f[1]+6.3245553203367588e-01*f[2]*g[8]+7.0710678118654757e-01*f[0]*g[8];
  fg[9] =  7.0710678118654757e-01*f[1]*g[16]+7.0710678118654757e-01*f[0]*g[9]+7.0710678118654757e-01*f[3]*g[67]+7.0710678118654757e-01*g[39]*f[2];
  fg[10] =  7.0710678118654757e-01*f[2]*g[41]+7.0710678118654757e-01*f[0]*g[10]+7.0710678118654757e-01*f[3]*g[69]+7.0710678118654757e-01*f[1]*g[17];
  fg[11] =  6.2105900340811881e-01*f[1]*g[31]+7.0710678118654757e-01*f[0]*g[11]+7.0710678118654757e-01*f[2]*g[0]+4.2163702135578390e-01*f[3]*g[31]+6.3245553203367588e-01*f[1]*g[1]+4.5175395145262565e-01*f[2]*g[11]+6.2105900340811881e-01*f[3]*g[1];
  fg[12] =  7.0710678118654757e-01*f[1]*g[20]+7.0710678118654757e-01*f[0]*g[12];
  fg[13] =  7.0710678118654757e-01*f[0]*g[13]+7.0710678118654757e-01*f[1]*g[23];
  fg[14] =  7.0710678118654757e-01*f[0]*g[14]+7.0710678118654757e-01*f[1]*g[28];
  fg[15] =  7.0710678118654757e-01*f[0]*g[15]+6.3245553203367588e-01*f[1]*g[36]+6.2105900340811881e-01*g[64]*f[2]+7.0710678118654757e-01*f[1]*g[7]+6.3245553203367588e-01*f[2]*g[15]+6.2105900340811881e-01*f[3]*g[36];
  fg[16] =  6.3245553203367588e-01*f[1]*g[39]+7.0710678118654757e-01*f[0]*g[16]+7.0710678118654757e-01*f[1]*g[9]+6.3245553203367588e-01*f[2]*g[16]+6.2105900340811881e-01*f[2]*g[67]+6.2105900340811881e-01*g[39]*f[3];
  fg[17] =  6.3245553203367588e-01*g[17]*f[2]+7.0710678118654757e-01*f[1]*g[10]+6.2105900340811881e-01*f[2]*g[69]+6.2105900340811881e-01*f[3]*g[41]+6.3245553203367588e-01*f[1]*g[41]+7.0710678118654757e-01*f[0]*g[17];
  fg[18] =  7.0710678118654757e-01*f[1]*g[35]+7.0710678118654757e-01*f[0]*g[18]+7.0710678118654757e-01*g[76]*f[3]+7.0710678118654757e-01*f[2]*g[60];
  fg[19] =  6.3245553203367588e-01*g[5]*f[1]+4.2163702135578390e-01*g[48]*f[3]+4.5175395145262565e-01*f[2]*g[19]+6.2105900340811881e-01*g[48]*f[1]+7.0710678118654757e-01*f[2]*g[2]+6.2105900340811881e-01*g[5]*f[3]+7.0710678118654757e-01*f[0]*g[19];
  fg[20] =  7.0710678118654757e-01*f[1]*g[12]+7.0710678118654757e-01*f[0]*g[20]+6.3245553203367588e-01*f[2]*g[20];
  fg[21] =  6.2105900340811881e-01*f[1]*g[50]+4.2163702135578390e-01*f[3]*g[50]+4.5175395145262565e-01*g[21]*f[2]+6.2105900340811881e-01*g[6]*f[3]+6.3245553203367588e-01*g[6]*f[1]+7.0710678118654757e-01*f[2]*g[3]+7.0710678118654757e-01*f[0]*g[21];
  fg[22] =  7.0710678118654757e-01*g[22]*f[0]+7.0710678118654757e-01*g[37]*f[1];
  fg[23] =  7.0710678118654757e-01*f[0]*g[23]+7.0710678118654757e-01*f[1]*g[13]+6.3245553203367588e-01*f[2]*g[23];
  fg[24] =  7.0710678118654757e-01*f[1]*g[38]+7.0710678118654757e-01*f[0]*g[24];
  fg[25] =  7.0710678118654757e-01*f[0]*g[25]+4.5175395145262565e-01*f[2]*g[25]+4.2163702135578390e-01*f[3]*g[54]+7.0710678118654757e-01*g[4]*f[2]+6.3245553203367588e-01*f[1]*g[8]+6.2105900340811881e-01*f[3]*g[8]+6.2105900340811881e-01*f[1]*g[54];
  fg[26] =  7.0710678118654757e-01*g[26]*f[0]+7.0710678118654757e-01*g[40]*f[1];
  fg[27] =  7.0710678118654757e-01*f[0]*g[27]+7.0710678118654757e-01*f[1]*g[43];
  fg[28] =  7.0710678118654757e-01*f[1]*g[14]+6.3245553203367588e-01*f[2]*g[28]+7.0710678118654757e-01*f[0]*g[28];
  fg[29] =  7.0710678118654757e-01*g[45]*f[1]+7.0710678118654757e-01*f[0]*g[29];
  fg[30] =  7.0710678118654757e-01*g[30]*f[0]+7.0710678118654757e-01*f[1]*g[46];
  fg[31] =  4.2163702135578390e-01*f[3]*g[11]+7.0710678118654757e-01*f[0]*g[31]+6.2105900340811881e-01*f[2]*g[1]+6.2105900340811881e-01*f[1]*g[11]+7.0710678118654757e-01*f[3]*g[0]+4.2163702135578390e-01*f[2]*g[31];
  fg[32] =  7.0710678118654757e-01*g[32]*f[0]+7.0710678118654757e-01*f[1]*g[49];
  fg[33] =  7.0710678118654757e-01*f[1]*g[52]+7.0710678118654757e-01*f[0]*g[33];
  fg[34] =  7.0710678118654757e-01*f[0]*g[34]+7.0710678118654757e-01*f[1]*g[57];
  fg[35] =  7.0710678118654757e-01*f[0]*g[35]+7.0710678118654757e-01*f[1]*g[18]+6.3245553203367588e-01*f[1]*g[60]+6.3245553203367588e-01*f[2]*g[35]+6.2105900340811881e-01*g[60]*f[3]+6.2105900340811881e-01*g[76]*f[2];
  fg[36] =  6.3245553203367588e-01*f[1]*g[15]+6.2105900340811881e-01*f[1]*g[64]+4.5175395145262565e-01*f[2]*g[36]+7.0710678118654757e-01*g[7]*f[2]+6.2105900340811881e-01*f[3]*g[15]+7.0710678118654757e-01*f[0]*g[36]+4.2163702135578390e-01*g[64]*f[3];
  fg[37] =  7.0710678118654757e-01*g[37]*f[0]+7.0710678118654757e-01*g[22]*f[1]+6.3245553203367588e-01*g[37]*f[2];
  fg[38] =  7.0710678118654757e-01*f[0]*g[38]+7.0710678118654757e-01*f[1]*g[24]+6.3245553203367588e-01*f[2]*g[38];
  fg[39] =  6.3245553203367588e-01*f[1]*g[16]+6.2105900340811881e-01*f[1]*g[67]+7.0710678118654757e-01*f[0]*g[39]+4.2163702135578390e-01*f[3]*g[67]+6.2105900340811881e-01*f[3]*g[16]+7.0710678118654757e-01*f[2]*g[9]+4.5175395145262565e-01*g[39]*f[2];
  fg[40] =  7.0710678118654757e-01*g[26]*f[1]+6.3245553203367588e-01*g[40]*f[2]+7.0710678118654757e-01*g[40]*f[0];
  fg[41] =  4.5175395145262565e-01*f[2]*g[41]+6.2105900340811881e-01*g[17]*f[3]+4.2163702135578390e-01*f[3]*g[69]+6.2105900340811881e-01*f[1]*g[69]+7.0710678118654757e-01*f[0]*g[41]+7.0710678118654757e-01*g[10]*f[2]+6.3245553203367588e-01*f[1]*g[17];
  fg[42] =  7.0710678118654757e-01*f[0]*g[42]+7.0710678118654757e-01*f[1]*g[61];
  fg[43] =  7.0710678118654757e-01*f[0]*g[43]+7.0710678118654757e-01*f[1]*g[27]+6.3245553203367588e-01*f[2]*g[43];
  fg[44] =  7.0710678118654757e-01*f[0]*g[44]+7.0710678118654757e-01*f[1]*g[62];
  fg[45] =  7.0710678118654757e-01*g[45]*f[0]+7.0710678118654757e-01*f[1]*g[29]+6.3245553203367588e-01*g[45]*f[2];
  fg[46] =  7.0710678118654757e-01*g[30]*f[1]+7.0710678118654757e-01*f[0]*g[46]+6.3245553203367588e-01*f[2]*g[46];
  fg[47] =  7.0710678118654757e-01*f[0]*g[47]+7.0710678118654757e-01*f[1]*g[63];
  fg[48] =  7.0710678118654757e-01*f[3]*g[2]+4.2163702135578390e-01*g[48]*f[2]+6.2105900340811881e-01*f[1]*g[19]+4.2163702135578390e-01*f[3]*g[19]+6.2105900340811881e-01*g[5]*f[2]+7.0710678118654757e-01*f[0]*g[48];
  fg[49] =  7.0710678118654757e-01*g[32]*f[1]+6.3245553203367588e-01*f[2]*g[49]+7.0710678118654757e-01*f[0]*g[49];
  fg[50] =  4.2163702135578390e-01*g[21]*f[3]+6.2105900340811881e-01*g[6]*f[2]+7.0710678118654757e-01*f[3]*g[3]+7.0710678118654757e-01*f[0]*g[50]+6.2105900340811881e-01*g[21]*f[1]+4.2163702135578390e-01*f[2]*g[50];
  fg[51] =  7.0710678118654757e-01*f[0]*g[51]+7.0710678118654757e-01*f[1]*g[65];
  fg[52] =  7.0710678118654757e-01*f[0]*g[52]+6.3245553203367588e-01*f[2]*g[52]+7.0710678118654757e-01*f[1]*g[33];
  fg[53] =  7.0710678118654757e-01*f[1]*g[66]+7.0710678118654757e-01*f[0]*g[53];
  fg[54] =  4.2163702135578390e-01*f[2]*g[54]+4.2163702135578390e-01*f[3]*g[25]+6.2105900340811881e-01*f[1]*g[25]+6.2105900340811881e-01*f[2]*g[8]+7.0710678118654757e-01*f[0]*g[54]+7.0710678118654757e-01*g[4]*f[3];
  fg[55] =  7.0710678118654757e-01*f[0]*g[55]+7.0710678118654757e-01*f[1]*g[68];
  fg[56] =  7.0710678118654757e-01*f[0]*g[56]+7.0710678118654757e-01*f[1]*g[71];
  fg[57] =  7.0710678118654757e-01*f[1]*g[34]+6.3245553203367588e-01*f[2]*g[57]+7.0710678118654757e-01*f[0]*g[57];
  fg[58] =  7.0710678118654757e-01*g[73]*f[1]+7.0710678118654757e-01*g[58]*f[0];
  fg[59] =  7.0710678118654757e-01*f[0]*g[59]+7.0710678118654757e-01*f[1]*g[74];
  fg[60] =  6.2105900340811881e-01*g[76]*f[1]+7.0710678118654757e-01*f[0]*g[60]+6.2105900340811881e-01*f[3]*g[35]+6.3245553203367588e-01*f[1]*g[35]+4.2163702135578390e-01*g[76]*f[3]+4.5175395145262565e-01*f[2]*g[60]+7.0710678118654757e-01*f[2]*g[18];
  fg[61] =  7.0710678118654757e-01*f[1]*g[42]+6.3245553203367588e-01*g[61]*f[2]+7.0710678118654757e-01*f[0]*g[61];
  fg[62] =  6.3245553203367588e-01*f[2]*g[62]+7.0710678118654757e-01*f[0]*g[62]+7.0710678118654757e-01*f[1]*g[44];
  fg[63] =  7.0710678118654757e-01*f[1]*g[47]+7.0710678118654757e-01*f[0]*g[63]+6.3245553203367588e-01*f[2]*g[63];
  fg[64] =  7.0710678118654757e-01*g[7]*f[3]+7.0710678118654757e-01*f[0]*g[64]+6.2105900340811881e-01*f[1]*g[36]+4.2163702135578390e-01*g[64]*f[2]+6.2105900340811881e-01*f[2]*g[15]+4.2163702135578390e-01*f[3]*g[36];
  fg[65] =  7.0710678118654757e-01*f[1]*g[51]+6.3245553203367588e-01*f[2]*g[65]+7.0710678118654757e-01*f[0]*g[65];
  fg[66] =  6.3245553203367588e-01*g[66]*f[2]+7.0710678118654757e-01*f[0]*g[66]+7.0710678118654757e-01*g[53]*f[1];
  fg[67] =  7.0710678118654757e-01*f[3]*g[9]+6.2105900340811881e-01*f[1]*g[39]+7.0710678118654757e-01*f[0]*g[67]+6.2105900340811881e-01*f[2]*g[16]+4.2163702135578390e-01*f[2]*g[67]+4.2163702135578390e-01*g[39]*f[3];
  fg[68] =  7.0710678118654757e-01*f[1]*g[55]+6.3245553203367588e-01*f[2]*g[68]+7.0710678118654757e-01*f[0]*g[68];
  fg[69] =  7.0710678118654757e-01*f[0]*g[69]+6.2105900340811881e-01*g[17]*f[2]+4.2163702135578390e-01*f[2]*g[69]+4.2163702135578390e-01*f[3]*g[41]+6.2105900340811881e-01*f[1]*g[41]+7.0710678118654757e-01*g[10]*f[3];
  fg[70] =  7.0710678118654757e-01*f[1]*g[77]+7.0710678118654757e-01*f[0]*g[70];
  fg[71] =  7.0710678118654757e-01*f[1]*g[56]+6.3245553203367588e-01*f[2]*g[71]+7.0710678118654757e-01*f[0]*g[71];
  fg[72] =  7.0710678118654757e-01*f[0]*g[72]+7.0710678118654757e-01*f[1]*g[78];
  fg[73] =  7.0710678118654757e-01*g[73]*f[0]+7.0710678118654757e-01*g[58]*f[1]+6.3245553203367588e-01*g[73]*f[2];
  fg[74] =  7.0710678118654757e-01*f[1]*g[59]+6.3245553203367588e-01*f[2]*g[74]+7.0710678118654757e-01*f[0]*g[74];
  fg[75] =  7.0710678118654757e-01*f[0]*g[75]+7.0710678118654757e-01*g[79]*f[1];
  fg[76] =  7.0710678118654757e-01*g[76]*f[0]+7.0710678118654757e-01*f[3]*g[18]+6.2105900340811881e-01*f[1]*g[60]+6.2105900340811881e-01*f[2]*g[35]+4.2163702135578390e-01*g[60]*f[3]+4.2163702135578390e-01*g[76]*f[2];
  fg[77] =  7.0710678118654757e-01*f[0]*g[77]+7.0710678118654757e-01*f[1]*g[70]+6.3245553203367588e-01*g[77]*f[2];
  fg[78] =  7.0710678118654757e-01*f[1]*g[72]+7.0710678118654757e-01*f[0]*g[78]+6.3245553203367588e-01*f[2]*g[78];
  fg[79] =  7.0710678118654757e-01*f[1]*g[75]+6.3245553203367588e-01*g[79]*f[2]+7.0710678118654757e-01*g[79]*f[0];
}

