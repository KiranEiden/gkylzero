// Thu Aug 26 15:51:57 2021
#include <gkyl_binop_cross_mul_ser.h>
GKYL_CU_DH
void
binop_cross_mul_2d_4d_ser_p2(const double *f, const double *g, double *fg )
{
  fg[0] =  5.0000000000000000e-01*f[0]*g[0]+5.0000000000000000e-01*g[11]*f[4]+5.0000000000000000e-01*g[12]*f[5]+5.0000000000000000e-01*f[2]*g[2]+5.0000000000000000e-01*f[3]*g[5]+5.0000000000000000e-01*f[1]*g[1]+5.0000000000000000e-01*g[20]*f[7]+5.0000000000000000e-01*g[19]*f[6];
  fg[1] =  5.0000000000000000e-01*f[0]*g[1]+4.4721359549995793e-01*f[3]*g[19]+5.0000000000000000e-01*g[12]*f[7]+5.0000000000000000e-01*f[1]*g[0]+4.4721359549995793e-01*g[5]*f[6]+4.4721359549995793e-01*f[1]*g[11]+5.0000000000000000e-01*g[20]*f[5]+5.0000000000000000e-01*f[2]*g[5]+5.0000000000000000e-01*f[3]*g[2]+4.4721359549995793e-01*f[4]*g[1];
  fg[2] =  5.0000000000000000e-01*f[3]*g[1]+4.4721359549995793e-01*g[2]*f[5]+4.4721359549995793e-01*g[5]*f[7]+5.0000000000000000e-01*g[11]*f[6]+4.4721359549995793e-01*f[2]*g[12]+5.0000000000000000e-01*g[5]*f[1]+5.0000000000000000e-01*f[2]*g[0]+4.4721359549995793e-01*f[3]*g[20]+5.0000000000000000e-01*g[19]*f[4]+5.0000000000000000e-01*f[0]*g[2];
  fg[3] =  5.0000000000000000e-01*f[7]*g[33]+5.0000000000000000e-01*g[32]*f[6]+5.0000000000000000e-01*g[22]*f[5]+5.0000000000000000e-01*f[2]*g[7]+5.0000000000000000e-01*f[3]*g[15]+5.0000000000000000e-01*f[0]*g[3]+5.0000000000000000e-01*g[21]*f[4]+5.0000000000000000e-01*g[6]*f[1];
  fg[4] =  5.0000000000000000e-01*g[35]*f[6]+5.0000000000000000e-01*g[26]*f[5]+5.0000000000000000e-01*f[1]*g[8]+5.0000000000000000e-01*f[2]*g[9]+5.0000000000000000e-01*g[25]*f[4]+5.0000000000000000e-01*f[0]*g[4]+5.0000000000000000e-01*f[3]*g[16]+5.0000000000000000e-01*g[36]*f[7];
  fg[5] =  4.4721359549995793e-01*g[5]*f[5]+4.4721359549995793e-01*f[2]*g[20]+5.0000000000000000e-01*f[3]*g[0]+4.4721359549995793e-01*g[5]*f[4]+4.4721359549995793e-01*g[2]*f[7]+4.4721359549995793e-01*f[1]*g[19]+5.0000000000000000e-01*g[5]*f[0]+5.0000000000000000e-01*f[1]*g[2]+5.0000000000000000e-01*f[2]*g[1]+4.4721359549995793e-01*g[1]*f[6]+4.0000000000000002e-01*g[19]*f[7]+4.4721359549995793e-01*f[3]*g[12]+4.0000000000000002e-01*g[20]*f[6]+4.4721359549995793e-01*f[3]*g[11];
  fg[6] =  5.0000000000000000e-01*g[22]*f[7]+5.0000000000000000e-01*f[2]*g[15]+5.0000000000000000e-01*f[3]*g[7]+4.4721359549995793e-01*f[3]*g[32]+5.0000000000000000e-01*g[6]*f[0]+4.4721359549995793e-01*g[15]*f[6]+5.0000000000000000e-01*f[1]*g[3]+4.4721359549995793e-01*g[21]*f[1]+5.0000000000000000e-01*f[5]*g[33]+4.4721359549995793e-01*g[6]*f[4];
  fg[7] =  5.0000000000000000e-01*f[2]*g[3]+5.0000000000000000e-01*f[3]*g[6]+5.0000000000000000e-01*g[32]*f[4]+4.4721359549995793e-01*f[3]*g[33]+4.4721359549995793e-01*g[15]*f[7]+4.4721359549995793e-01*g[7]*f[5]+5.0000000000000000e-01*g[21]*f[6]+4.4721359549995793e-01*f[2]*g[22]+5.0000000000000000e-01*g[7]*f[0]+5.0000000000000000e-01*f[1]*g[15];
  fg[8] =  5.0000000000000000e-01*f[0]*g[8]+5.0000000000000000e-01*f[2]*g[16]+5.0000000000000000e-01*f[1]*g[4]+5.0000000000000000e-01*g[26]*f[7]+5.0000000000000000e-01*f[5]*g[36]+4.4721359549995793e-01*g[16]*f[6]+4.4721359549995793e-01*g[8]*f[4]+4.4721359549995793e-01*f[1]*g[25]+5.0000000000000000e-01*f[3]*g[9]+4.4721359549995793e-01*f[3]*g[35];
  fg[9] =  4.4721359549995793e-01*f[3]*g[36]+5.0000000000000000e-01*f[3]*g[8]+5.0000000000000000e-01*g[25]*f[6]+5.0000000000000000e-01*g[35]*f[4]+4.4721359549995793e-01*g[9]*f[5]+4.4721359549995793e-01*g[16]*f[7]+4.4721359549995793e-01*f[2]*g[26]+5.0000000000000000e-01*f[2]*g[4]+5.0000000000000000e-01*f[1]*g[16]+5.0000000000000000e-01*f[0]*g[9];
  fg[10] =  5.0000000000000000e-01*g[10]*f[0]+5.0000000000000000e-01*g[31]*f[3]+5.0000000000000000e-01*g[18]*f[2]+5.0000000000000000e-01*g[37]*f[4]+5.0000000000000000e-01*g[45]*f[7]+5.0000000000000000e-01*g[44]*f[6]+5.0000000000000000e-01*g[38]*f[5]+5.0000000000000000e-01*g[17]*f[1];
  fg[11] =  3.1943828249996997e-01*g[11]*f[4]+4.4721359549995793e-01*f[3]*g[5]+4.4721359549995793e-01*f[1]*g[1]+5.0000000000000000e-01*g[2]*f[6]+5.0000000000000000e-01*f[2]*g[19]+5.0000000000000000e-01*g[0]*f[4]+4.4721359549995793e-01*g[20]*f[7]+5.0000000000000000e-01*f[0]*g[11]+3.1943828249996997e-01*g[19]*f[6];
  fg[12] =  3.1943828249996997e-01*g[12]*f[5]+5.0000000000000000e-01*g[12]*f[0]+4.4721359549995793e-01*f[2]*g[2]+4.4721359549995793e-01*f[3]*g[5]+5.0000000000000000e-01*g[1]*f[7]+5.0000000000000000e-01*f[5]*g[0]+3.1943828249996997e-01*g[20]*f[7]+5.0000000000000000e-01*f[1]*g[20]+4.4721359549995793e-01*g[19]*f[6];
  fg[13] =  5.0000000000000000e-01*g[23]*f[1]+5.0000000000000000e-01*g[34]*f[3]+5.0000000000000000e-01*f[0]*g[13]+5.0000000000000000e-01*f[2]*g[24];
  fg[14] =  5.0000000000000000e-01*f[1]*g[28]+5.0000000000000000e-01*f[3]*g[41]+5.0000000000000000e-01*f[2]*g[29]+5.0000000000000000e-01*g[14]*f[0];
  fg[15] =  4.0000000000000002e-01*g[33]*f[6]+4.4721359549995793e-01*f[3]*g[21]+4.4721359549995793e-01*f[3]*g[22]+4.0000000000000002e-01*g[32]*f[7]+4.4721359549995793e-01*g[32]*f[1]+5.0000000000000000e-01*f[0]*g[15]+5.0000000000000000e-01*f[3]*g[3]+4.4721359549995793e-01*g[6]*f[6]+5.0000000000000000e-01*f[2]*g[6]+5.0000000000000000e-01*g[7]*f[1]+4.4721359549995793e-01*g[15]*f[5]+4.4721359549995793e-01*g[15]*f[4]+4.4721359549995793e-01*f[2]*g[33]+4.4721359549995793e-01*g[7]*f[7];
  fg[16] =  4.4721359549995793e-01*g[9]*f[7]+4.0000000000000002e-01*g[35]*f[7]+5.0000000000000000e-01*f[3]*g[4]+4.4721359549995793e-01*f[3]*g[26]+4.0000000000000002e-01*g[36]*f[6]+5.0000000000000000e-01*f[2]*g[8]+5.0000000000000000e-01*f[0]*g[16]+5.0000000000000000e-01*f[1]*g[9]+4.4721359549995793e-01*f[1]*g[35]+4.4721359549995793e-01*f[2]*g[36]+4.4721359549995793e-01*f[3]*g[25]+4.4721359549995793e-01*g[16]*f[5]+4.4721359549995793e-01*g[8]*f[6]+4.4721359549995793e-01*g[16]*f[4];
  fg[17] =  4.4721359549995793e-01*g[37]*f[1]+5.0000000000000000e-01*g[45]*f[5]+5.0000000000000000e-01*f[2]*g[31]+5.0000000000000000e-01*g[10]*f[1]+4.4721359549995793e-01*f[3]*g[44]+5.0000000000000000e-01*g[18]*f[3]+5.0000000000000000e-01*g[38]*f[7]+4.4721359549995793e-01*g[31]*f[6]+4.4721359549995793e-01*g[17]*f[4]+5.0000000000000000e-01*g[17]*f[0];
  fg[18] =  5.0000000000000000e-01*g[37]*f[6]+5.0000000000000000e-01*f[3]*g[17]+4.4721359549995793e-01*f[2]*g[38]+4.4721359549995793e-01*f[3]*g[45]+5.0000000000000000e-01*f[2]*g[10]+5.0000000000000000e-01*g[44]*f[4]+5.0000000000000000e-01*g[31]*f[1]+4.4721359549995793e-01*g[31]*f[7]+4.4721359549995793e-01*g[18]*f[5]+5.0000000000000000e-01*g[18]*f[0];
  fg[19] =  4.4721359549995793e-01*f[3]*g[1]+5.0000000000000000e-01*f[0]*g[19]+5.0000000000000000e-01*g[2]*f[4]+5.0000000000000000e-01*f[2]*g[11]+4.0000000000000002e-01*g[5]*f[7]+3.1943828249996997e-01*g[11]*f[6]+4.4721359549995793e-01*g[5]*f[1]+4.4721359549995793e-01*g[12]*f[6]+4.4721359549995793e-01*g[19]*f[5]+4.0000000000000002e-01*f[3]*g[20]+3.1943828249996997e-01*g[19]*f[4]+5.0000000000000000e-01*g[0]*f[6];
  fg[20] =  4.0000000000000002e-01*f[3]*g[19]+3.1943828249996997e-01*g[12]*f[7]+4.4721359549995793e-01*g[11]*f[7]+4.0000000000000002e-01*g[5]*f[6]+5.0000000000000000e-01*g[12]*f[1]+5.0000000000000000e-01*f[0]*g[20]+3.1943828249996997e-01*g[20]*f[5]+4.4721359549995793e-01*g[20]*f[4]+5.0000000000000000e-01*f[5]*g[1]+5.0000000000000000e-01*g[0]*f[7]+4.4721359549995793e-01*f[2]*g[5]+4.4721359549995793e-01*f[3]*g[2];
  fg[21] =  4.4721359549995793e-01*f[7]*g[33]+3.1943828249996997e-01*g[32]*f[6]+5.0000000000000000e-01*g[3]*f[4]+4.4721359549995793e-01*f[3]*g[15]+3.1943828249996997e-01*g[21]*f[4]+4.4721359549995793e-01*g[6]*f[1]+5.0000000000000000e-01*g[7]*f[6]+5.0000000000000000e-01*g[21]*f[0]+5.0000000000000000e-01*f[2]*g[32];
  fg[22] =  5.0000000000000000e-01*f[1]*g[33]+3.1943828249996997e-01*f[7]*g[33]+4.4721359549995793e-01*g[32]*f[6]+3.1943828249996997e-01*g[22]*f[5]+5.0000000000000000e-01*f[0]*g[22]+4.4721359549995793e-01*f[2]*g[7]+4.4721359549995793e-01*f[3]*g[15]+5.0000000000000000e-01*g[6]*f[7]+5.0000000000000000e-01*g[3]*f[5];
  fg[23] =  5.0000000000000000e-01*f[2]*g[34]+5.0000000000000000e-01*g[23]*f[0]+4.4721359549995793e-01*g[23]*f[4]+5.0000000000000000e-01*f[3]*g[24]+4.4721359549995793e-01*g[34]*f[6]+5.0000000000000000e-01*f[1]*g[13];
  fg[24] =  5.0000000000000000e-01*f[2]*g[13]+5.0000000000000000e-01*f[0]*g[24]+5.0000000000000000e-01*g[34]*f[1]+4.4721359549995793e-01*g[34]*f[7]+5.0000000000000000e-01*f[3]*g[23]+4.4721359549995793e-01*f[5]*g[24];
  fg[25] =  5.0000000000000000e-01*g[9]*f[6]+3.1943828249996997e-01*g[35]*f[6]+4.4721359549995793e-01*f[1]*g[8]+5.0000000000000000e-01*f[2]*g[35]+3.1943828249996997e-01*g[25]*f[4]+5.0000000000000000e-01*g[4]*f[4]+4.4721359549995793e-01*f[3]*g[16]+4.4721359549995793e-01*g[36]*f[7]+5.0000000000000000e-01*f[0]*g[25];
  fg[26] =  4.4721359549995793e-01*g[35]*f[6]+3.1943828249996997e-01*g[26]*f[5]+5.0000000000000000e-01*g[4]*f[5]+5.0000000000000000e-01*f[1]*g[36]+4.4721359549995793e-01*f[2]*g[9]+5.0000000000000000e-01*g[26]*f[0]+4.4721359549995793e-01*f[3]*g[16]+3.1943828249996997e-01*g[36]*f[7]+5.0000000000000000e-01*g[8]*f[7];
  fg[27] =  5.0000000000000000e-01*g[39]*f[1]+5.0000000000000000e-01*g[27]*f[0]+5.0000000000000000e-01*f[2]*g[40]+5.0000000000000000e-01*f[3]*g[46];
  fg[28] =  4.4721359549995793e-01*f[6]*g[41]+5.0000000000000000e-01*f[0]*g[28]+5.0000000000000000e-01*g[14]*f[1]+4.4721359549995793e-01*f[4]*g[28]+5.0000000000000000e-01*f[3]*g[29]+5.0000000000000000e-01*f[2]*g[41];
  fg[29] =  4.4721359549995793e-01*g[29]*f[5]+5.0000000000000000e-01*f[2]*g[14]+5.0000000000000000e-01*f[1]*g[41]+4.4721359549995793e-01*f[7]*g[41]+5.0000000000000000e-01*f[3]*g[28]+5.0000000000000000e-01*f[0]*g[29];
  fg[30] =  5.0000000000000000e-01*g[30]*f[0]+5.0000000000000000e-01*f[2]*g[43]+5.0000000000000000e-01*g[42]*f[1]+5.0000000000000000e-01*f[3]*g[47];
  fg[31] =  4.4721359549995793e-01*f[1]*g[44]+4.0000000000000002e-01*g[45]*f[6]+4.4721359549995793e-01*f[3]*g[37]+4.4721359549995793e-01*f[3]*g[38]+5.0000000000000000e-01*f[2]*g[17]+5.0000000000000000e-01*g[31]*f[0]+4.4721359549995793e-01*f[2]*g[45]+5.0000000000000000e-01*g[18]*f[1]+5.0000000000000000e-01*f[3]*g[10]+4.4721359549995793e-01*g[31]*f[4]+4.4721359549995793e-01*g[17]*f[6]+4.0000000000000002e-01*g[44]*f[7]+4.4721359549995793e-01*g[31]*f[5]+4.4721359549995793e-01*g[18]*f[7];
  fg[32] =  4.4721359549995793e-01*g[32]*f[5]+4.4721359549995793e-01*f[3]*g[6]+4.4721359549995793e-01*g[22]*f[6]+3.1943828249996997e-01*g[32]*f[4]+5.0000000000000000e-01*g[32]*f[0]+4.0000000000000002e-01*f[3]*g[33]+5.0000000000000000e-01*g[7]*f[4]+5.0000000000000000e-01*f[2]*g[21]+4.0000000000000002e-01*g[15]*f[7]+5.0000000000000000e-01*g[3]*f[6]+3.1943828249996997e-01*g[21]*f[6]+4.4721359549995793e-01*f[1]*g[15];
  fg[33] =  5.0000000000000000e-01*g[22]*f[1]+3.1943828249996997e-01*g[22]*f[7]+4.4721359549995793e-01*f[2]*g[15]+4.4721359549995793e-01*f[3]*g[7]+4.0000000000000002e-01*f[3]*g[32]+5.0000000000000000e-01*f[0]*g[33]+4.0000000000000002e-01*g[15]*f[6]+4.4721359549995793e-01*f[4]*g[33]+4.4721359549995793e-01*g[21]*f[7]+5.0000000000000000e-01*g[6]*f[5]+3.1943828249996997e-01*f[5]*g[33]+5.0000000000000000e-01*g[3]*f[7];
  fg[34] =  4.4721359549995793e-01*g[23]*f[6]+5.0000000000000000e-01*f[1]*g[24]+4.4721359549995793e-01*f[7]*g[24]+5.0000000000000000e-01*f[2]*g[23]+5.0000000000000000e-01*f[3]*g[13]+5.0000000000000000e-01*g[34]*f[0]+4.4721359549995793e-01*g[34]*f[4]+4.4721359549995793e-01*g[34]*f[5];
  fg[35] =  4.0000000000000002e-01*f[3]*g[36]+4.4721359549995793e-01*f[3]*g[8]+3.1943828249996997e-01*g[25]*f[6]+5.0000000000000000e-01*g[4]*f[6]+5.0000000000000000e-01*g[9]*f[4]+3.1943828249996997e-01*g[35]*f[4]+4.4721359549995793e-01*g[26]*f[6]+4.4721359549995793e-01*g[35]*f[5]+5.0000000000000000e-01*f[2]*g[25]+4.0000000000000002e-01*g[16]*f[7]+4.4721359549995793e-01*f[1]*g[16]+5.0000000000000000e-01*f[0]*g[35];
  fg[36] =  5.0000000000000000e-01*g[4]*f[7]+4.4721359549995793e-01*g[25]*f[7]+5.0000000000000000e-01*f[0]*g[36]+4.4721359549995793e-01*f[2]*g[16]+5.0000000000000000e-01*g[26]*f[1]+3.1943828249996997e-01*g[26]*f[7]+5.0000000000000000e-01*f[5]*g[8]+3.1943828249996997e-01*f[5]*g[36]+4.0000000000000002e-01*g[16]*f[6]+4.4721359549995793e-01*g[36]*f[4]+4.4721359549995793e-01*f[3]*g[9]+4.0000000000000002e-01*f[3]*g[35];
  fg[37] =  4.4721359549995793e-01*g[31]*f[3]+3.1943828249996997e-01*g[37]*f[4]+4.4721359549995793e-01*g[45]*f[7]+5.0000000000000000e-01*f[0]*g[37]+5.0000000000000000e-01*g[10]*f[4]+5.0000000000000000e-01*f[2]*g[44]+5.0000000000000000e-01*g[18]*f[6]+3.1943828249996997e-01*g[44]*f[6]+4.4721359549995793e-01*g[17]*f[1];
  fg[38] =  4.4721359549995793e-01*g[31]*f[3]+4.4721359549995793e-01*g[18]*f[2]+5.0000000000000000e-01*g[45]*f[1]+5.0000000000000000e-01*g[10]*f[5]+3.1943828249996997e-01*g[45]*f[7]+5.0000000000000000e-01*f[0]*g[38]+4.4721359549995793e-01*g[44]*f[6]+3.1943828249996997e-01*g[38]*f[5]+5.0000000000000000e-01*g[17]*f[7];
  fg[39] =  5.0000000000000000e-01*f[3]*g[40]+5.0000000000000000e-01*f[2]*g[46]+5.0000000000000000e-01*g[39]*f[0]+5.0000000000000000e-01*g[27]*f[1]+4.4721359549995793e-01*g[39]*f[4]+4.4721359549995793e-01*g[46]*f[6];
  fg[40] =  5.0000000000000000e-01*f[2]*g[27]+5.0000000000000000e-01*g[40]*f[0]+5.0000000000000000e-01*f[3]*g[39]+4.4721359549995793e-01*g[40]*f[5]+5.0000000000000000e-01*f[1]*g[46]+4.4721359549995793e-01*g[46]*f[7];
  fg[41] =  5.0000000000000000e-01*g[14]*f[3]+4.4721359549995793e-01*g[29]*f[7]+5.0000000000000000e-01*f[0]*g[41]+4.4721359549995793e-01*g[28]*f[6]+5.0000000000000000e-01*f[1]*g[29]+4.4721359549995793e-01*f[4]*g[41]+5.0000000000000000e-01*f[2]*g[28]+4.4721359549995793e-01*f[5]*g[41];
  fg[42] =  5.0000000000000000e-01*f[2]*g[47]+5.0000000000000000e-01*g[30]*f[1]+4.4721359549995793e-01*g[47]*f[6]+4.4721359549995793e-01*g[42]*f[4]+5.0000000000000000e-01*g[42]*f[0]+5.0000000000000000e-01*f[3]*g[43];
  fg[43] =  4.4721359549995793e-01*g[43]*f[5]+4.4721359549995793e-01*g[47]*f[7]+5.0000000000000000e-01*g[42]*f[3]+5.0000000000000000e-01*g[47]*f[1]+5.0000000000000000e-01*f[2]*g[30]+5.0000000000000000e-01*f[0]*g[43];
  fg[44] =  3.1943828249996997e-01*g[37]*f[6]+5.0000000000000000e-01*f[0]*g[44]+5.0000000000000000e-01*g[10]*f[6]+4.4721359549995793e-01*f[3]*g[17]+4.0000000000000002e-01*f[3]*g[45]+5.0000000000000000e-01*g[18]*f[4]+3.1943828249996997e-01*g[44]*f[4]+4.4721359549995793e-01*g[31]*f[1]+4.0000000000000002e-01*g[31]*f[7]+4.4721359549995793e-01*g[44]*f[5]+4.4721359549995793e-01*g[38]*f[6]+5.0000000000000000e-01*f[2]*g[37];
  fg[45] =  4.4721359549995793e-01*g[45]*f[4]+5.0000000000000000e-01*g[10]*f[7]+3.1943828249996997e-01*g[45]*f[5]+4.4721359549995793e-01*f[2]*g[31]+4.4721359549995793e-01*g[37]*f[7]+4.0000000000000002e-01*f[3]*g[44]+5.0000000000000000e-01*g[45]*f[0]+4.4721359549995793e-01*g[18]*f[3]+3.1943828249996997e-01*g[38]*f[7]+4.0000000000000002e-01*g[31]*f[6]+5.0000000000000000e-01*g[17]*f[5]+5.0000000000000000e-01*f[1]*g[38];
  fg[46] =  4.4721359549995793e-01*g[39]*f[6]+5.0000000000000000e-01*f[2]*g[39]+4.4721359549995793e-01*g[46]*f[4]+5.0000000000000000e-01*f[0]*g[46]+4.4721359549995793e-01*g[40]*f[7]+4.4721359549995793e-01*g[46]*f[5]+5.0000000000000000e-01*g[40]*f[1]+5.0000000000000000e-01*g[27]*f[3];
  fg[47] =  4.4721359549995793e-01*g[47]*f[5]+4.4721359549995793e-01*g[47]*f[4]+5.0000000000000000e-01*g[42]*f[2]+5.0000000000000000e-01*g[47]*f[0]+5.0000000000000000e-01*f[3]*g[30]+4.4721359549995793e-01*g[43]*f[7]+5.0000000000000000e-01*f[1]*g[43]+4.4721359549995793e-01*g[42]*f[6];
}
