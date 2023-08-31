#include <gkyl_dg_diffusion_kernels.h> 
GKYL_CU_DH double dg_diffusion4_iso_euler_surfz_3x_tensor_p2(const double* w, const double* dx, double D, 
  const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: Cell-center coordinates
  // dxv[NDIM]: Cell spacing
  // D: Diffusion coefficient in the center cell
  // ql: Input field in the left cell
  // qc: Input field in the center cell
  // qr: Input field in the right cell
  // out: Incremented output

  const double dx1 = 2.0/dx[2]; 
  const double J = -1.0*pow(dx1, 4.0);

  const double *q0l = &ql[0]; 
  const double *q0c = &qc[0]; 
  const double *q0r = &qr[0]; 
  double *out0 = &out[0]; 

  out0[0] += J*D*((-6.708203932499369*q0r[9])-6.708203932499369*q0l[9]+13.41640786499874*q0c[9]+8.11898816047911*q0r[3]-8.11898816047911*q0l[3]-4.6875*q0r[0]-4.6875*q0l[0]+9.375*q0c[0]); 
  out0[1] += J*D*((-6.708203932499369*q0r[15])-6.708203932499369*q0l[15]+13.41640786499874*q0c[15]+8.11898816047911*q0r[5]-8.11898816047911*q0l[5]-4.6875*q0r[1]-4.6875*q0l[1]+9.375*q0c[1]); 
  out0[2] += J*D*((-6.708203932499369*q0r[16])-6.708203932499369*q0l[16]+13.41640786499874*q0c[16]+8.11898816047911*q0r[6]-8.11898816047911*q0l[6]-4.6875*q0r[2]-4.6875*q0l[2]+9.375*q0c[2]); 
  out0[3] += J*D*((-9.077304717673634*q0r[9])+9.077304717673634*q0l[9]+12.65625*q0r[3]+12.65625*q0l[3]+30.9375*q0c[3]-8.11898816047911*q0r[0]+8.11898816047911*q0l[0]); 
  out0[4] += J*D*((-6.708203932499369*q0r[19])-6.708203932499369*q0l[19]+13.41640786499874*q0c[19]+8.11898816047911*q0r[10]-8.11898816047911*q0l[10]-4.6875*q0r[4]-4.6875*q0l[4]+9.375*q0c[4]); 
  out0[5] += J*D*((-9.077304717673634*q0r[15])+9.077304717673634*q0l[15]+12.65625*q0r[5]+12.65625*q0l[5]+30.9375*q0c[5]-8.11898816047911*q0r[1]+8.11898816047911*q0l[1]); 
  out0[6] += J*D*((-9.077304717673634*q0r[16])+9.077304717673634*q0l[16]+12.65625*q0r[6]+12.65625*q0l[6]+30.9375*q0c[6]-8.11898816047911*q0r[2]+8.11898816047911*q0l[2]); 
  out0[7] += J*D*((-6.708203932499369*q0r[21])-6.708203932499369*q0l[21]+13.41640786499874*q0c[21]+8.118988160479114*q0r[13]-8.118988160479114*q0l[13]-4.6875*q0r[7]-4.6875*q0l[7]+9.375*q0c[7]); 
  out0[8] += J*D*((-6.708203932499369*q0r[22])-6.708203932499369*q0l[22]+13.41640786499874*q0c[22]+8.118988160479114*q0r[14]-8.118988160479114*q0l[14]-4.6875*q0r[8]-4.6875*q0l[8]+9.375*q0c[8]); 
  out0[9] += J*D*((-0.65625*q0r[9])-0.65625*q0l[9]+40.6875*q0c[9]+4.720198453190289*q0r[3]-4.720198453190289*q0l[3]-4.192627457812106*q0r[0]-4.192627457812106*q0l[0]+8.385254915624213*q0c[0]); 
  out0[10] += J*D*((-9.077304717673634*q0r[19])+9.077304717673634*q0l[19]+12.65625*q0r[10]+12.65625*q0l[10]+30.9375*q0c[10]-8.11898816047911*q0r[4]+8.11898816047911*q0l[4]); 
  out0[11] += J*D*((-6.708203932499369*q0r[24])-6.708203932499369*q0l[24]+13.41640786499874*q0c[24]+8.118988160479114*q0r[17]-8.118988160479114*q0l[17]-4.6875*q0r[11]-4.6875*q0l[11]+9.375*q0c[11]); 
  out0[12] += J*D*((-6.708203932499369*q0r[25])-6.708203932499369*q0l[25]+13.41640786499874*q0c[25]+8.118988160479114*q0r[18]-8.118988160479114*q0l[18]-4.6875*q0r[12]-4.6875*q0l[12]+9.375*q0c[12]); 
  out0[13] += J*D*((-9.077304717673634*q0r[21])+9.077304717673634*q0l[21]+12.65625*q0r[13]+12.65625*q0l[13]+30.9375*q0c[13]-8.118988160479114*q0r[7]+8.118988160479114*q0l[7]); 
  out0[14] += J*D*((-9.077304717673634*q0r[22])+9.077304717673634*q0l[22]+12.65625*q0r[14]+12.65625*q0l[14]+30.9375*q0c[14]-8.118988160479114*q0r[8]+8.118988160479114*q0l[8]); 
  out0[15] += J*D*((-0.65625*q0r[15])-0.65625*q0l[15]+40.6875*q0c[15]+4.72019845319029*q0r[5]-4.72019845319029*q0l[5]-4.192627457812105*q0r[1]-4.192627457812105*q0l[1]+8.38525491562421*q0c[1]); 
  out0[16] += J*D*((-0.65625*q0r[16])-0.65625*q0l[16]+40.6875*q0c[16]+4.72019845319029*q0r[6]-4.72019845319029*q0l[6]-4.192627457812105*q0r[2]-4.192627457812105*q0l[2]+8.38525491562421*q0c[2]); 
  out0[17] += J*D*((-9.077304717673634*q0r[24])+9.077304717673634*q0l[24]+12.65625*q0r[17]+12.65625*q0l[17]+30.9375*q0c[17]-8.118988160479114*q0r[11]+8.118988160479114*q0l[11]); 
  out0[18] += J*D*((-9.077304717673634*q0r[25])+9.077304717673634*q0l[25]+12.65625*q0r[18]+12.65625*q0l[18]+30.9375*q0c[18]-8.118988160479114*q0r[12]+8.118988160479114*q0l[12]); 
  out0[19] += J*D*((-0.65625*q0r[19])-0.65625*q0l[19]+40.6875*q0c[19]+4.720198453190289*q0r[10]-4.720198453190289*q0l[10]-4.192627457812106*q0r[4]-4.192627457812106*q0l[4]+8.385254915624213*q0c[4]); 
  out0[20] += J*D*((-6.708203932499369*q0r[26])-6.708203932499369*q0l[26]+13.41640786499874*q0c[26]+8.11898816047911*q0r[23]-8.11898816047911*q0l[23]-4.6875*q0r[20]-4.6875*q0l[20]+9.375*q0c[20]); 
  out0[21] += J*D*((-0.65625*q0r[21])-0.65625*q0l[21]+40.6875*q0c[21]+4.72019845319029*q0r[13]-4.72019845319029*q0l[13]-4.192627457812106*q0r[7]-4.192627457812106*q0l[7]+8.385254915624213*q0c[7]); 
  out0[22] += J*D*((-0.65625*q0r[22])-0.65625*q0l[22]+40.6875*q0c[22]+4.72019845319029*q0r[14]-4.72019845319029*q0l[14]-4.192627457812106*q0r[8]-4.192627457812106*q0l[8]+8.385254915624213*q0c[8]); 
  out0[23] += J*D*((-9.077304717673634*q0r[26])+9.077304717673634*q0l[26]+12.65625*q0r[23]+12.65625*q0l[23]+30.9375*q0c[23]-8.11898816047911*q0r[20]+8.11898816047911*q0l[20]); 
  out0[24] += J*D*((-0.65625*q0r[24])-0.65625*q0l[24]+40.6875*q0c[24]+4.720198453190289*q0r[17]-4.720198453190289*q0l[17]-4.192627457812105*q0r[11]-4.192627457812105*q0l[11]+8.38525491562421*q0c[11]); 
  out0[25] += J*D*((-0.65625*q0r[25])-0.65625*q0l[25]+40.6875*q0c[25]+4.720198453190289*q0r[18]-4.720198453190289*q0l[18]-4.192627457812105*q0r[12]-4.192627457812105*q0l[12]+8.38525491562421*q0c[12]); 
  out0[26] += J*D*((-0.65625*q0r[26])-0.65625*q0l[26]+40.6875*q0c[26]+4.720198453190289*q0r[23]-4.720198453190289*q0l[23]-4.192627457812106*q0r[20]-4.192627457812106*q0l[20]+8.385254915624213*q0c[20]); 

  const double *q1l = &ql[27]; 
  const double *q1c = &qc[27]; 
  const double *q1r = &qr[27]; 
  double *out1 = &out[27]; 

  out1[0] += J*D*((-6.708203932499369*q1r[9])-6.708203932499369*q1l[9]+13.41640786499874*q1c[9]+8.11898816047911*q1r[3]-8.11898816047911*q1l[3]-4.6875*q1r[0]-4.6875*q1l[0]+9.375*q1c[0]); 
  out1[1] += J*D*((-6.708203932499369*q1r[15])-6.708203932499369*q1l[15]+13.41640786499874*q1c[15]+8.11898816047911*q1r[5]-8.11898816047911*q1l[5]-4.6875*q1r[1]-4.6875*q1l[1]+9.375*q1c[1]); 
  out1[2] += J*D*((-6.708203932499369*q1r[16])-6.708203932499369*q1l[16]+13.41640786499874*q1c[16]+8.11898816047911*q1r[6]-8.11898816047911*q1l[6]-4.6875*q1r[2]-4.6875*q1l[2]+9.375*q1c[2]); 
  out1[3] += J*D*((-9.077304717673634*q1r[9])+9.077304717673634*q1l[9]+12.65625*q1r[3]+12.65625*q1l[3]+30.9375*q1c[3]-8.11898816047911*q1r[0]+8.11898816047911*q1l[0]); 
  out1[4] += J*D*((-6.708203932499369*q1r[19])-6.708203932499369*q1l[19]+13.41640786499874*q1c[19]+8.11898816047911*q1r[10]-8.11898816047911*q1l[10]-4.6875*q1r[4]-4.6875*q1l[4]+9.375*q1c[4]); 
  out1[5] += J*D*((-9.077304717673634*q1r[15])+9.077304717673634*q1l[15]+12.65625*q1r[5]+12.65625*q1l[5]+30.9375*q1c[5]-8.11898816047911*q1r[1]+8.11898816047911*q1l[1]); 
  out1[6] += J*D*((-9.077304717673634*q1r[16])+9.077304717673634*q1l[16]+12.65625*q1r[6]+12.65625*q1l[6]+30.9375*q1c[6]-8.11898816047911*q1r[2]+8.11898816047911*q1l[2]); 
  out1[7] += J*D*((-6.708203932499369*q1r[21])-6.708203932499369*q1l[21]+13.41640786499874*q1c[21]+8.118988160479114*q1r[13]-8.118988160479114*q1l[13]-4.6875*q1r[7]-4.6875*q1l[7]+9.375*q1c[7]); 
  out1[8] += J*D*((-6.708203932499369*q1r[22])-6.708203932499369*q1l[22]+13.41640786499874*q1c[22]+8.118988160479114*q1r[14]-8.118988160479114*q1l[14]-4.6875*q1r[8]-4.6875*q1l[8]+9.375*q1c[8]); 
  out1[9] += J*D*((-0.65625*q1r[9])-0.65625*q1l[9]+40.6875*q1c[9]+4.720198453190289*q1r[3]-4.720198453190289*q1l[3]-4.192627457812106*q1r[0]-4.192627457812106*q1l[0]+8.385254915624213*q1c[0]); 
  out1[10] += J*D*((-9.077304717673634*q1r[19])+9.077304717673634*q1l[19]+12.65625*q1r[10]+12.65625*q1l[10]+30.9375*q1c[10]-8.11898816047911*q1r[4]+8.11898816047911*q1l[4]); 
  out1[11] += J*D*((-6.708203932499369*q1r[24])-6.708203932499369*q1l[24]+13.41640786499874*q1c[24]+8.118988160479114*q1r[17]-8.118988160479114*q1l[17]-4.6875*q1r[11]-4.6875*q1l[11]+9.375*q1c[11]); 
  out1[12] += J*D*((-6.708203932499369*q1r[25])-6.708203932499369*q1l[25]+13.41640786499874*q1c[25]+8.118988160479114*q1r[18]-8.118988160479114*q1l[18]-4.6875*q1r[12]-4.6875*q1l[12]+9.375*q1c[12]); 
  out1[13] += J*D*((-9.077304717673634*q1r[21])+9.077304717673634*q1l[21]+12.65625*q1r[13]+12.65625*q1l[13]+30.9375*q1c[13]-8.118988160479114*q1r[7]+8.118988160479114*q1l[7]); 
  out1[14] += J*D*((-9.077304717673634*q1r[22])+9.077304717673634*q1l[22]+12.65625*q1r[14]+12.65625*q1l[14]+30.9375*q1c[14]-8.118988160479114*q1r[8]+8.118988160479114*q1l[8]); 
  out1[15] += J*D*((-0.65625*q1r[15])-0.65625*q1l[15]+40.6875*q1c[15]+4.72019845319029*q1r[5]-4.72019845319029*q1l[5]-4.192627457812105*q1r[1]-4.192627457812105*q1l[1]+8.38525491562421*q1c[1]); 
  out1[16] += J*D*((-0.65625*q1r[16])-0.65625*q1l[16]+40.6875*q1c[16]+4.72019845319029*q1r[6]-4.72019845319029*q1l[6]-4.192627457812105*q1r[2]-4.192627457812105*q1l[2]+8.38525491562421*q1c[2]); 
  out1[17] += J*D*((-9.077304717673634*q1r[24])+9.077304717673634*q1l[24]+12.65625*q1r[17]+12.65625*q1l[17]+30.9375*q1c[17]-8.118988160479114*q1r[11]+8.118988160479114*q1l[11]); 
  out1[18] += J*D*((-9.077304717673634*q1r[25])+9.077304717673634*q1l[25]+12.65625*q1r[18]+12.65625*q1l[18]+30.9375*q1c[18]-8.118988160479114*q1r[12]+8.118988160479114*q1l[12]); 
  out1[19] += J*D*((-0.65625*q1r[19])-0.65625*q1l[19]+40.6875*q1c[19]+4.720198453190289*q1r[10]-4.720198453190289*q1l[10]-4.192627457812106*q1r[4]-4.192627457812106*q1l[4]+8.385254915624213*q1c[4]); 
  out1[20] += J*D*((-6.708203932499369*q1r[26])-6.708203932499369*q1l[26]+13.41640786499874*q1c[26]+8.11898816047911*q1r[23]-8.11898816047911*q1l[23]-4.6875*q1r[20]-4.6875*q1l[20]+9.375*q1c[20]); 
  out1[21] += J*D*((-0.65625*q1r[21])-0.65625*q1l[21]+40.6875*q1c[21]+4.72019845319029*q1r[13]-4.72019845319029*q1l[13]-4.192627457812106*q1r[7]-4.192627457812106*q1l[7]+8.385254915624213*q1c[7]); 
  out1[22] += J*D*((-0.65625*q1r[22])-0.65625*q1l[22]+40.6875*q1c[22]+4.72019845319029*q1r[14]-4.72019845319029*q1l[14]-4.192627457812106*q1r[8]-4.192627457812106*q1l[8]+8.385254915624213*q1c[8]); 
  out1[23] += J*D*((-9.077304717673634*q1r[26])+9.077304717673634*q1l[26]+12.65625*q1r[23]+12.65625*q1l[23]+30.9375*q1c[23]-8.11898816047911*q1r[20]+8.11898816047911*q1l[20]); 
  out1[24] += J*D*((-0.65625*q1r[24])-0.65625*q1l[24]+40.6875*q1c[24]+4.720198453190289*q1r[17]-4.720198453190289*q1l[17]-4.192627457812105*q1r[11]-4.192627457812105*q1l[11]+8.38525491562421*q1c[11]); 
  out1[25] += J*D*((-0.65625*q1r[25])-0.65625*q1l[25]+40.6875*q1c[25]+4.720198453190289*q1r[18]-4.720198453190289*q1l[18]-4.192627457812105*q1r[12]-4.192627457812105*q1l[12]+8.38525491562421*q1c[12]); 
  out1[26] += J*D*((-0.65625*q1r[26])-0.65625*q1l[26]+40.6875*q1c[26]+4.720198453190289*q1r[23]-4.720198453190289*q1l[23]-4.192627457812106*q1r[20]-4.192627457812106*q1l[20]+8.385254915624213*q1c[20]); 

  const double *q2l = &ql[54]; 
  const double *q2c = &qc[54]; 
  const double *q2r = &qr[54]; 
  double *out2 = &out[54]; 

  out2[0] += J*D*((-6.708203932499369*q2r[9])-6.708203932499369*q2l[9]+13.41640786499874*q2c[9]+8.11898816047911*q2r[3]-8.11898816047911*q2l[3]-4.6875*q2r[0]-4.6875*q2l[0]+9.375*q2c[0]); 
  out2[1] += J*D*((-6.708203932499369*q2r[15])-6.708203932499369*q2l[15]+13.41640786499874*q2c[15]+8.11898816047911*q2r[5]-8.11898816047911*q2l[5]-4.6875*q2r[1]-4.6875*q2l[1]+9.375*q2c[1]); 
  out2[2] += J*D*((-6.708203932499369*q2r[16])-6.708203932499369*q2l[16]+13.41640786499874*q2c[16]+8.11898816047911*q2r[6]-8.11898816047911*q2l[6]-4.6875*q2r[2]-4.6875*q2l[2]+9.375*q2c[2]); 
  out2[3] += J*D*((-9.077304717673634*q2r[9])+9.077304717673634*q2l[9]+12.65625*q2r[3]+12.65625*q2l[3]+30.9375*q2c[3]-8.11898816047911*q2r[0]+8.11898816047911*q2l[0]); 
  out2[4] += J*D*((-6.708203932499369*q2r[19])-6.708203932499369*q2l[19]+13.41640786499874*q2c[19]+8.11898816047911*q2r[10]-8.11898816047911*q2l[10]-4.6875*q2r[4]-4.6875*q2l[4]+9.375*q2c[4]); 
  out2[5] += J*D*((-9.077304717673634*q2r[15])+9.077304717673634*q2l[15]+12.65625*q2r[5]+12.65625*q2l[5]+30.9375*q2c[5]-8.11898816047911*q2r[1]+8.11898816047911*q2l[1]); 
  out2[6] += J*D*((-9.077304717673634*q2r[16])+9.077304717673634*q2l[16]+12.65625*q2r[6]+12.65625*q2l[6]+30.9375*q2c[6]-8.11898816047911*q2r[2]+8.11898816047911*q2l[2]); 
  out2[7] += J*D*((-6.708203932499369*q2r[21])-6.708203932499369*q2l[21]+13.41640786499874*q2c[21]+8.118988160479114*q2r[13]-8.118988160479114*q2l[13]-4.6875*q2r[7]-4.6875*q2l[7]+9.375*q2c[7]); 
  out2[8] += J*D*((-6.708203932499369*q2r[22])-6.708203932499369*q2l[22]+13.41640786499874*q2c[22]+8.118988160479114*q2r[14]-8.118988160479114*q2l[14]-4.6875*q2r[8]-4.6875*q2l[8]+9.375*q2c[8]); 
  out2[9] += J*D*((-0.65625*q2r[9])-0.65625*q2l[9]+40.6875*q2c[9]+4.720198453190289*q2r[3]-4.720198453190289*q2l[3]-4.192627457812106*q2r[0]-4.192627457812106*q2l[0]+8.385254915624213*q2c[0]); 
  out2[10] += J*D*((-9.077304717673634*q2r[19])+9.077304717673634*q2l[19]+12.65625*q2r[10]+12.65625*q2l[10]+30.9375*q2c[10]-8.11898816047911*q2r[4]+8.11898816047911*q2l[4]); 
  out2[11] += J*D*((-6.708203932499369*q2r[24])-6.708203932499369*q2l[24]+13.41640786499874*q2c[24]+8.118988160479114*q2r[17]-8.118988160479114*q2l[17]-4.6875*q2r[11]-4.6875*q2l[11]+9.375*q2c[11]); 
  out2[12] += J*D*((-6.708203932499369*q2r[25])-6.708203932499369*q2l[25]+13.41640786499874*q2c[25]+8.118988160479114*q2r[18]-8.118988160479114*q2l[18]-4.6875*q2r[12]-4.6875*q2l[12]+9.375*q2c[12]); 
  out2[13] += J*D*((-9.077304717673634*q2r[21])+9.077304717673634*q2l[21]+12.65625*q2r[13]+12.65625*q2l[13]+30.9375*q2c[13]-8.118988160479114*q2r[7]+8.118988160479114*q2l[7]); 
  out2[14] += J*D*((-9.077304717673634*q2r[22])+9.077304717673634*q2l[22]+12.65625*q2r[14]+12.65625*q2l[14]+30.9375*q2c[14]-8.118988160479114*q2r[8]+8.118988160479114*q2l[8]); 
  out2[15] += J*D*((-0.65625*q2r[15])-0.65625*q2l[15]+40.6875*q2c[15]+4.72019845319029*q2r[5]-4.72019845319029*q2l[5]-4.192627457812105*q2r[1]-4.192627457812105*q2l[1]+8.38525491562421*q2c[1]); 
  out2[16] += J*D*((-0.65625*q2r[16])-0.65625*q2l[16]+40.6875*q2c[16]+4.72019845319029*q2r[6]-4.72019845319029*q2l[6]-4.192627457812105*q2r[2]-4.192627457812105*q2l[2]+8.38525491562421*q2c[2]); 
  out2[17] += J*D*((-9.077304717673634*q2r[24])+9.077304717673634*q2l[24]+12.65625*q2r[17]+12.65625*q2l[17]+30.9375*q2c[17]-8.118988160479114*q2r[11]+8.118988160479114*q2l[11]); 
  out2[18] += J*D*((-9.077304717673634*q2r[25])+9.077304717673634*q2l[25]+12.65625*q2r[18]+12.65625*q2l[18]+30.9375*q2c[18]-8.118988160479114*q2r[12]+8.118988160479114*q2l[12]); 
  out2[19] += J*D*((-0.65625*q2r[19])-0.65625*q2l[19]+40.6875*q2c[19]+4.720198453190289*q2r[10]-4.720198453190289*q2l[10]-4.192627457812106*q2r[4]-4.192627457812106*q2l[4]+8.385254915624213*q2c[4]); 
  out2[20] += J*D*((-6.708203932499369*q2r[26])-6.708203932499369*q2l[26]+13.41640786499874*q2c[26]+8.11898816047911*q2r[23]-8.11898816047911*q2l[23]-4.6875*q2r[20]-4.6875*q2l[20]+9.375*q2c[20]); 
  out2[21] += J*D*((-0.65625*q2r[21])-0.65625*q2l[21]+40.6875*q2c[21]+4.72019845319029*q2r[13]-4.72019845319029*q2l[13]-4.192627457812106*q2r[7]-4.192627457812106*q2l[7]+8.385254915624213*q2c[7]); 
  out2[22] += J*D*((-0.65625*q2r[22])-0.65625*q2l[22]+40.6875*q2c[22]+4.72019845319029*q2r[14]-4.72019845319029*q2l[14]-4.192627457812106*q2r[8]-4.192627457812106*q2l[8]+8.385254915624213*q2c[8]); 
  out2[23] += J*D*((-9.077304717673634*q2r[26])+9.077304717673634*q2l[26]+12.65625*q2r[23]+12.65625*q2l[23]+30.9375*q2c[23]-8.11898816047911*q2r[20]+8.11898816047911*q2l[20]); 
  out2[24] += J*D*((-0.65625*q2r[24])-0.65625*q2l[24]+40.6875*q2c[24]+4.720198453190289*q2r[17]-4.720198453190289*q2l[17]-4.192627457812105*q2r[11]-4.192627457812105*q2l[11]+8.38525491562421*q2c[11]); 
  out2[25] += J*D*((-0.65625*q2r[25])-0.65625*q2l[25]+40.6875*q2c[25]+4.720198453190289*q2r[18]-4.720198453190289*q2l[18]-4.192627457812105*q2r[12]-4.192627457812105*q2l[12]+8.38525491562421*q2c[12]); 
  out2[26] += J*D*((-0.65625*q2r[26])-0.65625*q2l[26]+40.6875*q2c[26]+4.720198453190289*q2r[23]-4.720198453190289*q2l[23]-4.192627457812106*q2r[20]-4.192627457812106*q2l[20]+8.385254915624213*q2c[20]); 

  const double *q3l = &ql[81]; 
  const double *q3c = &qc[81]; 
  const double *q3r = &qr[81]; 
  double *out3 = &out[81]; 

  out3[0] += J*D*((-6.708203932499369*q3r[9])-6.708203932499369*q3l[9]+13.41640786499874*q3c[9]+8.11898816047911*q3r[3]-8.11898816047911*q3l[3]-4.6875*q3r[0]-4.6875*q3l[0]+9.375*q3c[0]); 
  out3[1] += J*D*((-6.708203932499369*q3r[15])-6.708203932499369*q3l[15]+13.41640786499874*q3c[15]+8.11898816047911*q3r[5]-8.11898816047911*q3l[5]-4.6875*q3r[1]-4.6875*q3l[1]+9.375*q3c[1]); 
  out3[2] += J*D*((-6.708203932499369*q3r[16])-6.708203932499369*q3l[16]+13.41640786499874*q3c[16]+8.11898816047911*q3r[6]-8.11898816047911*q3l[6]-4.6875*q3r[2]-4.6875*q3l[2]+9.375*q3c[2]); 
  out3[3] += J*D*((-9.077304717673634*q3r[9])+9.077304717673634*q3l[9]+12.65625*q3r[3]+12.65625*q3l[3]+30.9375*q3c[3]-8.11898816047911*q3r[0]+8.11898816047911*q3l[0]); 
  out3[4] += J*D*((-6.708203932499369*q3r[19])-6.708203932499369*q3l[19]+13.41640786499874*q3c[19]+8.11898816047911*q3r[10]-8.11898816047911*q3l[10]-4.6875*q3r[4]-4.6875*q3l[4]+9.375*q3c[4]); 
  out3[5] += J*D*((-9.077304717673634*q3r[15])+9.077304717673634*q3l[15]+12.65625*q3r[5]+12.65625*q3l[5]+30.9375*q3c[5]-8.11898816047911*q3r[1]+8.11898816047911*q3l[1]); 
  out3[6] += J*D*((-9.077304717673634*q3r[16])+9.077304717673634*q3l[16]+12.65625*q3r[6]+12.65625*q3l[6]+30.9375*q3c[6]-8.11898816047911*q3r[2]+8.11898816047911*q3l[2]); 
  out3[7] += J*D*((-6.708203932499369*q3r[21])-6.708203932499369*q3l[21]+13.41640786499874*q3c[21]+8.118988160479114*q3r[13]-8.118988160479114*q3l[13]-4.6875*q3r[7]-4.6875*q3l[7]+9.375*q3c[7]); 
  out3[8] += J*D*((-6.708203932499369*q3r[22])-6.708203932499369*q3l[22]+13.41640786499874*q3c[22]+8.118988160479114*q3r[14]-8.118988160479114*q3l[14]-4.6875*q3r[8]-4.6875*q3l[8]+9.375*q3c[8]); 
  out3[9] += J*D*((-0.65625*q3r[9])-0.65625*q3l[9]+40.6875*q3c[9]+4.720198453190289*q3r[3]-4.720198453190289*q3l[3]-4.192627457812106*q3r[0]-4.192627457812106*q3l[0]+8.385254915624213*q3c[0]); 
  out3[10] += J*D*((-9.077304717673634*q3r[19])+9.077304717673634*q3l[19]+12.65625*q3r[10]+12.65625*q3l[10]+30.9375*q3c[10]-8.11898816047911*q3r[4]+8.11898816047911*q3l[4]); 
  out3[11] += J*D*((-6.708203932499369*q3r[24])-6.708203932499369*q3l[24]+13.41640786499874*q3c[24]+8.118988160479114*q3r[17]-8.118988160479114*q3l[17]-4.6875*q3r[11]-4.6875*q3l[11]+9.375*q3c[11]); 
  out3[12] += J*D*((-6.708203932499369*q3r[25])-6.708203932499369*q3l[25]+13.41640786499874*q3c[25]+8.118988160479114*q3r[18]-8.118988160479114*q3l[18]-4.6875*q3r[12]-4.6875*q3l[12]+9.375*q3c[12]); 
  out3[13] += J*D*((-9.077304717673634*q3r[21])+9.077304717673634*q3l[21]+12.65625*q3r[13]+12.65625*q3l[13]+30.9375*q3c[13]-8.118988160479114*q3r[7]+8.118988160479114*q3l[7]); 
  out3[14] += J*D*((-9.077304717673634*q3r[22])+9.077304717673634*q3l[22]+12.65625*q3r[14]+12.65625*q3l[14]+30.9375*q3c[14]-8.118988160479114*q3r[8]+8.118988160479114*q3l[8]); 
  out3[15] += J*D*((-0.65625*q3r[15])-0.65625*q3l[15]+40.6875*q3c[15]+4.72019845319029*q3r[5]-4.72019845319029*q3l[5]-4.192627457812105*q3r[1]-4.192627457812105*q3l[1]+8.38525491562421*q3c[1]); 
  out3[16] += J*D*((-0.65625*q3r[16])-0.65625*q3l[16]+40.6875*q3c[16]+4.72019845319029*q3r[6]-4.72019845319029*q3l[6]-4.192627457812105*q3r[2]-4.192627457812105*q3l[2]+8.38525491562421*q3c[2]); 
  out3[17] += J*D*((-9.077304717673634*q3r[24])+9.077304717673634*q3l[24]+12.65625*q3r[17]+12.65625*q3l[17]+30.9375*q3c[17]-8.118988160479114*q3r[11]+8.118988160479114*q3l[11]); 
  out3[18] += J*D*((-9.077304717673634*q3r[25])+9.077304717673634*q3l[25]+12.65625*q3r[18]+12.65625*q3l[18]+30.9375*q3c[18]-8.118988160479114*q3r[12]+8.118988160479114*q3l[12]); 
  out3[19] += J*D*((-0.65625*q3r[19])-0.65625*q3l[19]+40.6875*q3c[19]+4.720198453190289*q3r[10]-4.720198453190289*q3l[10]-4.192627457812106*q3r[4]-4.192627457812106*q3l[4]+8.385254915624213*q3c[4]); 
  out3[20] += J*D*((-6.708203932499369*q3r[26])-6.708203932499369*q3l[26]+13.41640786499874*q3c[26]+8.11898816047911*q3r[23]-8.11898816047911*q3l[23]-4.6875*q3r[20]-4.6875*q3l[20]+9.375*q3c[20]); 
  out3[21] += J*D*((-0.65625*q3r[21])-0.65625*q3l[21]+40.6875*q3c[21]+4.72019845319029*q3r[13]-4.72019845319029*q3l[13]-4.192627457812106*q3r[7]-4.192627457812106*q3l[7]+8.385254915624213*q3c[7]); 
  out3[22] += J*D*((-0.65625*q3r[22])-0.65625*q3l[22]+40.6875*q3c[22]+4.72019845319029*q3r[14]-4.72019845319029*q3l[14]-4.192627457812106*q3r[8]-4.192627457812106*q3l[8]+8.385254915624213*q3c[8]); 
  out3[23] += J*D*((-9.077304717673634*q3r[26])+9.077304717673634*q3l[26]+12.65625*q3r[23]+12.65625*q3l[23]+30.9375*q3c[23]-8.11898816047911*q3r[20]+8.11898816047911*q3l[20]); 
  out3[24] += J*D*((-0.65625*q3r[24])-0.65625*q3l[24]+40.6875*q3c[24]+4.720198453190289*q3r[17]-4.720198453190289*q3l[17]-4.192627457812105*q3r[11]-4.192627457812105*q3l[11]+8.38525491562421*q3c[11]); 
  out3[25] += J*D*((-0.65625*q3r[25])-0.65625*q3l[25]+40.6875*q3c[25]+4.720198453190289*q3r[18]-4.720198453190289*q3l[18]-4.192627457812105*q3r[12]-4.192627457812105*q3l[12]+8.38525491562421*q3c[12]); 
  out3[26] += J*D*((-0.65625*q3r[26])-0.65625*q3l[26]+40.6875*q3c[26]+4.720198453190289*q3r[23]-4.720198453190289*q3l[23]-4.192627457812106*q3r[20]-4.192627457812106*q3l[20]+8.385254915624213*q3c[20]); 

  return 0.;

} 
