// Thu Jul 28 13:04:52 2022
#include <gkyl_binop_cross_mul_ser.h>
GKYL_CU_DH
void
binop_cross_mul_1d_3d_ser_p3(const double *f, const double *g, double *fg )
{
  double tmp[32] = {0.};
  tmp[0] =  7.0710678118654757e-01*f[1]*g[1]+7.0710678118654757e-01*f[0]*g[0]+7.0710678118654757e-01*g[17]*f[3]+7.0710678118654757e-01*f[2]*g[7];
  tmp[1] =  7.0710678118654757e-01*f[1]*g[0]+6.3245553203367588e-01*f[1]*g[7]+7.0710678118654757e-01*f[0]*g[1]+6.2105900340811881e-01*g[7]*f[3]+6.3245553203367588e-01*f[2]*g[1]+6.2105900340811881e-01*f[2]*g[17];
  tmp[2] =  7.0710678118654757e-01*f[2]*g[11]+7.0710678118654757e-01*f[1]*g[4]+7.0710678118654757e-01*g[23]*f[3]+7.0710678118654757e-01*f[0]*g[2];
  tmp[3] =  7.0710678118654757e-01*f[2]*g[13]+7.0710678118654757e-01*f[1]*g[5]+7.0710678118654757e-01*f[3]*g[25]+7.0710678118654757e-01*f[0]*g[3];
  tmp[4] =  6.2105900340811881e-01*g[23]*f[2]+7.0710678118654757e-01*f[0]*g[4]+6.2105900340811881e-01*f[3]*g[11]+6.3245553203367588e-01*f[2]*g[4]+7.0710678118654757e-01*f[1]*g[2]+6.3245553203367588e-01*f[1]*g[11];
  tmp[5] =  6.2105900340811881e-01*f[2]*g[25]+7.0710678118654757e-01*f[0]*g[5]+6.3245553203367588e-01*f[1]*g[13]+6.3245553203367588e-01*f[2]*g[5]+7.0710678118654757e-01*f[1]*g[3]+6.2105900340811881e-01*f[3]*g[13];
  tmp[6] =  7.0710678118654757e-01*f[0]*g[6]+7.0710678118654757e-01*f[2]*g[20]+7.0710678118654757e-01*f[3]*g[29]+7.0710678118654757e-01*f[1]*g[10];
  tmp[7] =  6.3245553203367588e-01*f[1]*g[1]+6.2105900340811881e-01*f[1]*g[17]+7.0710678118654757e-01*f[0]*g[7]+6.2105900340811881e-01*f[3]*g[1]+4.2163702135578390e-01*g[17]*f[3]+7.0710678118654757e-01*f[2]*g[0]+4.5175395145262565e-01*f[2]*g[7];
  tmp[8] =  7.0710678118654757e-01*f[0]*g[8]+7.0710678118654757e-01*g[12]*f[1];
  tmp[9] =  7.0710678118654757e-01*f[0]*g[9]+7.0710678118654757e-01*f[1]*g[15];
  tmp[10] =  7.0710678118654757e-01*f[1]*g[6]+6.2105900340811881e-01*f[2]*g[29]+6.3245553203367588e-01*g[10]*f[2]+7.0710678118654757e-01*f[0]*g[10]+6.3245553203367588e-01*f[1]*g[20]+6.2105900340811881e-01*f[3]*g[20];
  tmp[11] =  7.0710678118654757e-01*f[2]*g[2]+6.2105900340811881e-01*f[3]*g[4]+4.5175395145262565e-01*f[2]*g[11]+6.3245553203367588e-01*f[1]*g[4]+4.2163702135578390e-01*g[23]*f[3]+7.0710678118654757e-01*f[0]*g[11]+6.2105900340811881e-01*f[1]*g[23];
  tmp[12] =  7.0710678118654757e-01*f[1]*g[8]+6.3245553203367588e-01*g[12]*f[2]+7.0710678118654757e-01*g[12]*f[0];
  tmp[13] =  6.2105900340811881e-01*g[5]*f[3]+4.5175395145262565e-01*f[2]*g[13]+6.3245553203367588e-01*f[1]*g[5]+7.0710678118654757e-01*f[2]*g[3]+4.2163702135578390e-01*f[3]*g[25]+6.2105900340811881e-01*f[1]*g[25]+7.0710678118654757e-01*f[0]*g[13];
  tmp[14] =  7.0710678118654757e-01*f[1]*g[21]+7.0710678118654757e-01*g[14]*f[0];
  tmp[15] =  6.3245553203367588e-01*f[2]*g[15]+7.0710678118654757e-01*f[1]*g[9]+7.0710678118654757e-01*f[0]*g[15];
  tmp[16] =  7.0710678118654757e-01*f[1]*g[22]+7.0710678118654757e-01*f[0]*g[16];
  tmp[17] =  6.2105900340811881e-01*f[1]*g[7]+7.0710678118654757e-01*f[3]*g[0]+7.0710678118654757e-01*f[0]*g[17]+4.2163702135578390e-01*g[7]*f[3]+6.2105900340811881e-01*f[2]*g[1]+4.2163702135578390e-01*f[2]*g[17];
  tmp[18] =  7.0710678118654757e-01*g[18]*f[0]+7.0710678118654757e-01*f[1]*g[24];
  tmp[19] =  7.0710678118654757e-01*f[0]*g[19]+7.0710678118654757e-01*g[27]*f[1];
  tmp[20] =  4.5175395145262565e-01*f[2]*g[20]+6.2105900340811881e-01*g[10]*f[3]+4.2163702135578390e-01*f[3]*g[29]+7.0710678118654757e-01*f[2]*g[6]+6.3245553203367588e-01*f[1]*g[10]+6.2105900340811881e-01*f[1]*g[29]+7.0710678118654757e-01*f[0]*g[20];
  tmp[21] =  7.0710678118654757e-01*g[14]*f[1]+7.0710678118654757e-01*f[0]*g[21]+6.3245553203367588e-01*f[2]*g[21];
  tmp[22] =  7.0710678118654757e-01*f[0]*g[22]+7.0710678118654757e-01*f[1]*g[16]+6.3245553203367588e-01*f[2]*g[22];
  tmp[23] =  4.2163702135578390e-01*g[23]*f[2]+4.2163702135578390e-01*f[3]*g[11]+6.2105900340811881e-01*f[2]*g[4]+7.0710678118654757e-01*f[0]*g[23]+7.0710678118654757e-01*f[3]*g[2]+6.2105900340811881e-01*f[1]*g[11];
  tmp[24] =  7.0710678118654757e-01*f[0]*g[24]+7.0710678118654757e-01*g[18]*f[1]+6.3245553203367588e-01*f[2]*g[24];
  tmp[25] =  4.2163702135578390e-01*f[2]*g[25]+7.0710678118654757e-01*f[0]*g[25]+6.2105900340811881e-01*f[1]*g[13]+6.2105900340811881e-01*f[2]*g[5]+7.0710678118654757e-01*f[3]*g[3]+4.2163702135578390e-01*f[3]*g[13];
  tmp[26] =  7.0710678118654757e-01*f[0]*g[26]+7.0710678118654757e-01*f[1]*g[30];
  tmp[27] =  7.0710678118654757e-01*f[1]*g[19]+6.3245553203367588e-01*g[27]*f[2]+7.0710678118654757e-01*f[0]*g[27];
  tmp[28] =  7.0710678118654757e-01*f[0]*g[28]+7.0710678118654757e-01*g[31]*f[1];
  tmp[29] =  4.2163702135578390e-01*f[2]*g[29]+6.2105900340811881e-01*g[10]*f[2]+6.2105900340811881e-01*f[1]*g[20]+7.0710678118654757e-01*f[0]*g[29]+4.2163702135578390e-01*f[3]*g[20]+7.0710678118654757e-01*g[6]*f[3];
  tmp[30] =  7.0710678118654757e-01*f[0]*g[30]+7.0710678118654757e-01*f[1]*g[26]+6.3245553203367588e-01*f[2]*g[30];
  tmp[31] =  7.0710678118654757e-01*g[31]*f[0]+7.0710678118654757e-01*f[1]*g[28]+6.3245553203367588e-01*g[31]*f[2];
 
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
}

