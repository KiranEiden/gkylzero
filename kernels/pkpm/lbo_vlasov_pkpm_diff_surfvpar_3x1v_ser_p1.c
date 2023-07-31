#include <gkyl_lbo_vlasov_pkpm_kernels.h> 
GKYL_CU_DH double lbo_vlasov_pkpm_diff_surfvpar_3x1v_ser_p1(const double *w, const double *dxv, const double *nuVtSq, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[4]:         cell-center coordinates. 
  // dxv[4]:       cell spacing. 
  // nuVtSqSum[8]: Sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      Input Distribution function [F_0, T_perp G = T_perp (F_1 - F_0)] in left/center/right cells 
  // out:           Incremented distribution function in cell 
  double rdvSq4 = 4.0/(dxv[3]*dxv[3]); 
  const double *F_0l = &fl[0]; 
  const double *G_1l = &fl[24]; 
  const double *F_0c = &fc[0]; 
  const double *G_1c = &fc[24]; 
  const double *F_0r = &fr[0]; 
  const double *G_1r = &fr[24]; 
  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[24]; 
  double incr_F_0[24] = {0.0}; 
  double incr_G_1[24] = {0.0}; 

  double F_0_xx[24] = {0.0}; 
  double G_1_xx[24] = {0.0}; 
  F_0_xx[0] = 0.6708203932499369*F_0r[16]+0.6708203932499369*F_0l[16]-1.341640786499874*F_0c[16]-1.190784930203603*F_0r[4]+1.190784930203603*F_0l[4]+0.9375*F_0r[0]+0.9375*F_0l[0]-1.875*F_0c[0]; 
  F_0_xx[1] = 0.6708203932499369*F_0r[17]+0.6708203932499369*F_0l[17]-1.341640786499874*F_0c[17]-1.190784930203603*F_0r[8]+1.190784930203603*F_0l[8]+0.9375*F_0r[1]+0.9375*F_0l[1]-1.875*F_0c[1]; 
  F_0_xx[2] = 0.6708203932499369*F_0r[18]+0.6708203932499369*F_0l[18]-1.341640786499874*F_0c[18]-1.190784930203603*F_0r[9]+1.190784930203603*F_0l[9]+0.9375*F_0r[2]+0.9375*F_0l[2]-1.875*F_0c[2]; 
  F_0_xx[3] = 0.6708203932499369*F_0r[19]+0.6708203932499369*F_0l[19]-1.341640786499874*F_0c[19]-1.190784930203603*F_0r[10]+1.190784930203603*F_0l[10]+0.9375*F_0r[3]+0.9375*F_0l[3]-1.875*F_0c[3]; 
  F_0_xx[4] = 0.7382874503707888*F_0r[16]-0.7382874503707888*F_0l[16]-1.453125*F_0r[4]-1.453125*F_0l[4]-5.34375*F_0c[4]+1.190784930203603*F_0r[0]-1.190784930203603*F_0l[0]; 
  F_0_xx[5] = 0.6708203932499369*F_0r[20]+0.6708203932499369*F_0l[20]-1.341640786499874*F_0c[20]-1.190784930203603*F_0r[12]+1.190784930203603*F_0l[12]+0.9375*F_0r[5]+0.9375*F_0l[5]-1.875*F_0c[5]; 
  F_0_xx[6] = 0.6708203932499369*F_0r[21]+0.6708203932499369*F_0l[21]-1.341640786499874*F_0c[21]-1.190784930203603*F_0r[13]+1.190784930203603*F_0l[13]+0.9375*F_0r[6]+0.9375*F_0l[6]-1.875*F_0c[6]; 
  F_0_xx[7] = 0.6708203932499369*F_0r[22]+0.6708203932499369*F_0l[22]-1.341640786499874*F_0c[22]-1.190784930203603*F_0r[14]+1.190784930203603*F_0l[14]+0.9375*F_0r[7]+0.9375*F_0l[7]-1.875*F_0c[7]; 
  F_0_xx[8] = 0.7382874503707888*F_0r[17]-0.7382874503707888*F_0l[17]-1.453125*F_0r[8]-1.453125*F_0l[8]-5.34375*F_0c[8]+1.190784930203603*F_0r[1]-1.190784930203603*F_0l[1]; 
  F_0_xx[9] = 0.7382874503707888*F_0r[18]-0.7382874503707888*F_0l[18]-1.453125*F_0r[9]-1.453125*F_0l[9]-5.34375*F_0c[9]+1.190784930203603*F_0r[2]-1.190784930203603*F_0l[2]; 
  F_0_xx[10] = 0.7382874503707888*F_0r[19]-0.7382874503707888*F_0l[19]-1.453125*F_0r[10]-1.453125*F_0l[10]-5.34375*F_0c[10]+1.190784930203603*F_0r[3]-1.190784930203603*F_0l[3]; 
  F_0_xx[11] = 0.6708203932499369*F_0r[23]+0.6708203932499369*F_0l[23]-1.341640786499874*F_0c[23]-1.190784930203603*F_0r[15]+1.190784930203603*F_0l[15]+0.9375*F_0r[11]+0.9375*F_0l[11]-1.875*F_0c[11]; 
  F_0_xx[12] = 0.7382874503707888*F_0r[20]-0.7382874503707888*F_0l[20]-1.453125*F_0r[12]-1.453125*F_0l[12]-5.34375*F_0c[12]+1.190784930203603*F_0r[5]-1.190784930203603*F_0l[5]; 
  F_0_xx[13] = 0.7382874503707888*F_0r[21]-0.7382874503707888*F_0l[21]-1.453125*F_0r[13]-1.453125*F_0l[13]-5.34375*F_0c[13]+1.190784930203603*F_0r[6]-1.190784930203603*F_0l[6]; 
  F_0_xx[14] = 0.7382874503707888*F_0r[22]-0.7382874503707888*F_0l[22]-1.453125*F_0r[14]-1.453125*F_0l[14]-5.34375*F_0c[14]+1.190784930203603*F_0r[7]-1.190784930203603*F_0l[7]; 
  F_0_xx[15] = 0.7382874503707888*F_0r[23]-0.7382874503707888*F_0l[23]-1.453125*F_0r[15]-1.453125*F_0l[15]-5.34375*F_0c[15]+1.190784930203603*F_0r[11]-1.190784930203603*F_0l[11]; 
  F_0_xx[16] = (-0.140625*F_0r[16])-0.140625*F_0l[16]-6.28125*F_0c[16]-0.3025768239224545*F_0r[4]+0.3025768239224545*F_0l[4]+0.4192627457812106*F_0r[0]+0.4192627457812106*F_0l[0]-0.8385254915624212*F_0c[0]; 
  F_0_xx[17] = (-0.140625*F_0r[17])-0.140625*F_0l[17]-6.28125*F_0c[17]-0.3025768239224544*F_0r[8]+0.3025768239224544*F_0l[8]+0.4192627457812105*F_0r[1]+0.4192627457812105*F_0l[1]-0.8385254915624211*F_0c[1]; 
  F_0_xx[18] = (-0.140625*F_0r[18])-0.140625*F_0l[18]-6.28125*F_0c[18]-0.3025768239224544*F_0r[9]+0.3025768239224544*F_0l[9]+0.4192627457812105*F_0r[2]+0.4192627457812105*F_0l[2]-0.8385254915624211*F_0c[2]; 
  F_0_xx[19] = (-0.140625*F_0r[19])-0.140625*F_0l[19]-6.28125*F_0c[19]-0.3025768239224544*F_0r[10]+0.3025768239224544*F_0l[10]+0.4192627457812105*F_0r[3]+0.4192627457812105*F_0l[3]-0.8385254915624211*F_0c[3]; 
  F_0_xx[20] = (-0.140625*F_0r[20])-0.140625*F_0l[20]-6.28125*F_0c[20]-0.3025768239224545*F_0r[12]+0.3025768239224545*F_0l[12]+0.4192627457812106*F_0r[5]+0.4192627457812106*F_0l[5]-0.8385254915624212*F_0c[5]; 
  F_0_xx[21] = (-0.140625*F_0r[21])-0.140625*F_0l[21]-6.28125*F_0c[21]-0.3025768239224545*F_0r[13]+0.3025768239224545*F_0l[13]+0.4192627457812106*F_0r[6]+0.4192627457812106*F_0l[6]-0.8385254915624212*F_0c[6]; 
  F_0_xx[22] = (-0.140625*F_0r[22])-0.140625*F_0l[22]-6.28125*F_0c[22]-0.3025768239224545*F_0r[14]+0.3025768239224545*F_0l[14]+0.4192627457812106*F_0r[7]+0.4192627457812106*F_0l[7]-0.8385254915624212*F_0c[7]; 
  F_0_xx[23] = (-0.140625*F_0r[23])-0.140625*F_0l[23]-6.28125*F_0c[23]-0.3025768239224544*F_0r[15]+0.3025768239224544*F_0l[15]+0.4192627457812105*F_0r[11]+0.4192627457812105*F_0l[11]-0.8385254915624211*F_0c[11]; 
  G_1_xx[0] = 0.6708203932499369*G_1r[16]+0.6708203932499369*G_1l[16]-1.341640786499874*G_1c[16]-1.190784930203603*G_1r[4]+1.190784930203603*G_1l[4]+0.9375*G_1r[0]+0.9375*G_1l[0]-1.875*G_1c[0]; 
  G_1_xx[1] = 0.6708203932499369*G_1r[17]+0.6708203932499369*G_1l[17]-1.341640786499874*G_1c[17]-1.190784930203603*G_1r[8]+1.190784930203603*G_1l[8]+0.9375*G_1r[1]+0.9375*G_1l[1]-1.875*G_1c[1]; 
  G_1_xx[2] = 0.6708203932499369*G_1r[18]+0.6708203932499369*G_1l[18]-1.341640786499874*G_1c[18]-1.190784930203603*G_1r[9]+1.190784930203603*G_1l[9]+0.9375*G_1r[2]+0.9375*G_1l[2]-1.875*G_1c[2]; 
  G_1_xx[3] = 0.6708203932499369*G_1r[19]+0.6708203932499369*G_1l[19]-1.341640786499874*G_1c[19]-1.190784930203603*G_1r[10]+1.190784930203603*G_1l[10]+0.9375*G_1r[3]+0.9375*G_1l[3]-1.875*G_1c[3]; 
  G_1_xx[4] = 0.7382874503707888*G_1r[16]-0.7382874503707888*G_1l[16]-1.453125*G_1r[4]-1.453125*G_1l[4]-5.34375*G_1c[4]+1.190784930203603*G_1r[0]-1.190784930203603*G_1l[0]; 
  G_1_xx[5] = 0.6708203932499369*G_1r[20]+0.6708203932499369*G_1l[20]-1.341640786499874*G_1c[20]-1.190784930203603*G_1r[12]+1.190784930203603*G_1l[12]+0.9375*G_1r[5]+0.9375*G_1l[5]-1.875*G_1c[5]; 
  G_1_xx[6] = 0.6708203932499369*G_1r[21]+0.6708203932499369*G_1l[21]-1.341640786499874*G_1c[21]-1.190784930203603*G_1r[13]+1.190784930203603*G_1l[13]+0.9375*G_1r[6]+0.9375*G_1l[6]-1.875*G_1c[6]; 
  G_1_xx[7] = 0.6708203932499369*G_1r[22]+0.6708203932499369*G_1l[22]-1.341640786499874*G_1c[22]-1.190784930203603*G_1r[14]+1.190784930203603*G_1l[14]+0.9375*G_1r[7]+0.9375*G_1l[7]-1.875*G_1c[7]; 
  G_1_xx[8] = 0.7382874503707888*G_1r[17]-0.7382874503707888*G_1l[17]-1.453125*G_1r[8]-1.453125*G_1l[8]-5.34375*G_1c[8]+1.190784930203603*G_1r[1]-1.190784930203603*G_1l[1]; 
  G_1_xx[9] = 0.7382874503707888*G_1r[18]-0.7382874503707888*G_1l[18]-1.453125*G_1r[9]-1.453125*G_1l[9]-5.34375*G_1c[9]+1.190784930203603*G_1r[2]-1.190784930203603*G_1l[2]; 
  G_1_xx[10] = 0.7382874503707888*G_1r[19]-0.7382874503707888*G_1l[19]-1.453125*G_1r[10]-1.453125*G_1l[10]-5.34375*G_1c[10]+1.190784930203603*G_1r[3]-1.190784930203603*G_1l[3]; 
  G_1_xx[11] = 0.6708203932499369*G_1r[23]+0.6708203932499369*G_1l[23]-1.341640786499874*G_1c[23]-1.190784930203603*G_1r[15]+1.190784930203603*G_1l[15]+0.9375*G_1r[11]+0.9375*G_1l[11]-1.875*G_1c[11]; 
  G_1_xx[12] = 0.7382874503707888*G_1r[20]-0.7382874503707888*G_1l[20]-1.453125*G_1r[12]-1.453125*G_1l[12]-5.34375*G_1c[12]+1.190784930203603*G_1r[5]-1.190784930203603*G_1l[5]; 
  G_1_xx[13] = 0.7382874503707888*G_1r[21]-0.7382874503707888*G_1l[21]-1.453125*G_1r[13]-1.453125*G_1l[13]-5.34375*G_1c[13]+1.190784930203603*G_1r[6]-1.190784930203603*G_1l[6]; 
  G_1_xx[14] = 0.7382874503707888*G_1r[22]-0.7382874503707888*G_1l[22]-1.453125*G_1r[14]-1.453125*G_1l[14]-5.34375*G_1c[14]+1.190784930203603*G_1r[7]-1.190784930203603*G_1l[7]; 
  G_1_xx[15] = 0.7382874503707888*G_1r[23]-0.7382874503707888*G_1l[23]-1.453125*G_1r[15]-1.453125*G_1l[15]-5.34375*G_1c[15]+1.190784930203603*G_1r[11]-1.190784930203603*G_1l[11]; 
  G_1_xx[16] = (-0.140625*G_1r[16])-0.140625*G_1l[16]-6.28125*G_1c[16]-0.3025768239224545*G_1r[4]+0.3025768239224545*G_1l[4]+0.4192627457812106*G_1r[0]+0.4192627457812106*G_1l[0]-0.8385254915624212*G_1c[0]; 
  G_1_xx[17] = (-0.140625*G_1r[17])-0.140625*G_1l[17]-6.28125*G_1c[17]-0.3025768239224544*G_1r[8]+0.3025768239224544*G_1l[8]+0.4192627457812105*G_1r[1]+0.4192627457812105*G_1l[1]-0.8385254915624211*G_1c[1]; 
  G_1_xx[18] = (-0.140625*G_1r[18])-0.140625*G_1l[18]-6.28125*G_1c[18]-0.3025768239224544*G_1r[9]+0.3025768239224544*G_1l[9]+0.4192627457812105*G_1r[2]+0.4192627457812105*G_1l[2]-0.8385254915624211*G_1c[2]; 
  G_1_xx[19] = (-0.140625*G_1r[19])-0.140625*G_1l[19]-6.28125*G_1c[19]-0.3025768239224544*G_1r[10]+0.3025768239224544*G_1l[10]+0.4192627457812105*G_1r[3]+0.4192627457812105*G_1l[3]-0.8385254915624211*G_1c[3]; 
  G_1_xx[20] = (-0.140625*G_1r[20])-0.140625*G_1l[20]-6.28125*G_1c[20]-0.3025768239224545*G_1r[12]+0.3025768239224545*G_1l[12]+0.4192627457812106*G_1r[5]+0.4192627457812106*G_1l[5]-0.8385254915624212*G_1c[5]; 
  G_1_xx[21] = (-0.140625*G_1r[21])-0.140625*G_1l[21]-6.28125*G_1c[21]-0.3025768239224545*G_1r[13]+0.3025768239224545*G_1l[13]+0.4192627457812106*G_1r[6]+0.4192627457812106*G_1l[6]-0.8385254915624212*G_1c[6]; 
  G_1_xx[22] = (-0.140625*G_1r[22])-0.140625*G_1l[22]-6.28125*G_1c[22]-0.3025768239224545*G_1r[14]+0.3025768239224545*G_1l[14]+0.4192627457812106*G_1r[7]+0.4192627457812106*G_1l[7]-0.8385254915624212*G_1c[7]; 
  G_1_xx[23] = (-0.140625*G_1r[23])-0.140625*G_1l[23]-6.28125*G_1c[23]-0.3025768239224544*G_1r[15]+0.3025768239224544*G_1l[15]+0.4192627457812105*G_1r[11]+0.4192627457812105*G_1l[11]-0.8385254915624211*G_1c[11]; 

  incr_F_0[0] = 0.3535533905932737*nuVtSq[7]*F_0_xx[11]+0.3535533905932737*nuVtSq[6]*F_0_xx[7]+0.3535533905932737*nuVtSq[5]*F_0_xx[6]+0.3535533905932737*nuVtSq[4]*F_0_xx[5]+0.3535533905932737*F_0_xx[3]*nuVtSq[3]+0.3535533905932737*F_0_xx[2]*nuVtSq[2]+0.3535533905932737*F_0_xx[1]*nuVtSq[1]+0.3535533905932737*F_0_xx[0]*nuVtSq[0]; 
  incr_F_0[1] = 0.3535533905932737*nuVtSq[6]*F_0_xx[11]+0.3535533905932737*F_0_xx[7]*nuVtSq[7]+0.3535533905932737*nuVtSq[3]*F_0_xx[6]+0.3535533905932737*F_0_xx[3]*nuVtSq[5]+0.3535533905932737*nuVtSq[2]*F_0_xx[5]+0.3535533905932737*F_0_xx[2]*nuVtSq[4]+0.3535533905932737*F_0_xx[0]*nuVtSq[1]+0.3535533905932737*nuVtSq[0]*F_0_xx[1]; 
  incr_F_0[2] = 0.3535533905932737*nuVtSq[5]*F_0_xx[11]+0.3535533905932737*F_0_xx[6]*nuVtSq[7]+0.3535533905932737*nuVtSq[3]*F_0_xx[7]+0.3535533905932737*F_0_xx[3]*nuVtSq[6]+0.3535533905932737*nuVtSq[1]*F_0_xx[5]+0.3535533905932737*F_0_xx[1]*nuVtSq[4]+0.3535533905932737*F_0_xx[0]*nuVtSq[2]+0.3535533905932737*nuVtSq[0]*F_0_xx[2]; 
  incr_F_0[3] = 0.3535533905932737*nuVtSq[4]*F_0_xx[11]+0.3535533905932737*F_0_xx[5]*nuVtSq[7]+0.3535533905932737*nuVtSq[2]*F_0_xx[7]+0.3535533905932737*F_0_xx[2]*nuVtSq[6]+0.3535533905932737*nuVtSq[1]*F_0_xx[6]+0.3535533905932737*F_0_xx[1]*nuVtSq[5]+0.3535533905932737*F_0_xx[0]*nuVtSq[3]+0.3535533905932737*nuVtSq[0]*F_0_xx[3]; 
  incr_F_0[4] = 0.3535533905932737*nuVtSq[7]*F_0_xx[15]+0.3535533905932737*nuVtSq[6]*F_0_xx[14]+0.3535533905932737*nuVtSq[5]*F_0_xx[13]+0.3535533905932737*nuVtSq[4]*F_0_xx[12]+0.3535533905932737*nuVtSq[3]*F_0_xx[10]+0.3535533905932737*nuVtSq[2]*F_0_xx[9]+0.3535533905932737*nuVtSq[1]*F_0_xx[8]+0.3535533905932737*nuVtSq[0]*F_0_xx[4]; 
  incr_F_0[5] = 0.3535533905932737*nuVtSq[3]*F_0_xx[11]+0.3535533905932737*F_0_xx[3]*nuVtSq[7]+0.3535533905932737*nuVtSq[5]*F_0_xx[7]+0.3535533905932737*F_0_xx[6]*nuVtSq[6]+0.3535533905932737*nuVtSq[0]*F_0_xx[5]+0.3535533905932737*F_0_xx[0]*nuVtSq[4]+0.3535533905932737*F_0_xx[1]*nuVtSq[2]+0.3535533905932737*nuVtSq[1]*F_0_xx[2]; 
  incr_F_0[6] = 0.3535533905932737*nuVtSq[2]*F_0_xx[11]+0.3535533905932737*F_0_xx[2]*nuVtSq[7]+0.3535533905932737*nuVtSq[4]*F_0_xx[7]+0.3535533905932737*F_0_xx[5]*nuVtSq[6]+0.3535533905932737*nuVtSq[0]*F_0_xx[6]+0.3535533905932737*F_0_xx[0]*nuVtSq[5]+0.3535533905932737*F_0_xx[1]*nuVtSq[3]+0.3535533905932737*nuVtSq[1]*F_0_xx[3]; 
  incr_F_0[7] = 0.3535533905932737*nuVtSq[1]*F_0_xx[11]+0.3535533905932737*F_0_xx[1]*nuVtSq[7]+0.3535533905932737*nuVtSq[0]*F_0_xx[7]+0.3535533905932737*F_0_xx[0]*nuVtSq[6]+0.3535533905932737*nuVtSq[4]*F_0_xx[6]+0.3535533905932737*F_0_xx[5]*nuVtSq[5]+0.3535533905932737*F_0_xx[2]*nuVtSq[3]+0.3535533905932737*nuVtSq[2]*F_0_xx[3]; 
  incr_F_0[8] = 0.3535533905932737*nuVtSq[6]*F_0_xx[15]+0.3535533905932737*nuVtSq[7]*F_0_xx[14]+0.3535533905932737*nuVtSq[3]*F_0_xx[13]+0.3535533905932737*nuVtSq[2]*F_0_xx[12]+0.3535533905932737*nuVtSq[5]*F_0_xx[10]+0.3535533905932737*nuVtSq[4]*F_0_xx[9]+0.3535533905932737*nuVtSq[0]*F_0_xx[8]+0.3535533905932737*nuVtSq[1]*F_0_xx[4]; 
  incr_F_0[9] = 0.3535533905932737*nuVtSq[5]*F_0_xx[15]+0.3535533905932737*nuVtSq[3]*F_0_xx[14]+0.3535533905932737*nuVtSq[7]*F_0_xx[13]+0.3535533905932737*nuVtSq[1]*F_0_xx[12]+0.3535533905932737*nuVtSq[6]*F_0_xx[10]+0.3535533905932737*nuVtSq[0]*F_0_xx[9]+0.3535533905932737*nuVtSq[4]*F_0_xx[8]+0.3535533905932737*nuVtSq[2]*F_0_xx[4]; 
  incr_F_0[10] = 0.3535533905932737*nuVtSq[4]*F_0_xx[15]+0.3535533905932737*nuVtSq[2]*F_0_xx[14]+0.3535533905932737*nuVtSq[1]*F_0_xx[13]+0.3535533905932737*nuVtSq[7]*F_0_xx[12]+0.3535533905932737*nuVtSq[0]*F_0_xx[10]+0.3535533905932737*nuVtSq[6]*F_0_xx[9]+0.3535533905932737*nuVtSq[5]*F_0_xx[8]+0.3535533905932737*nuVtSq[3]*F_0_xx[4]; 
  incr_F_0[11] = 0.3535533905932737*nuVtSq[0]*F_0_xx[11]+0.3535533905932737*F_0_xx[0]*nuVtSq[7]+0.3535533905932737*nuVtSq[1]*F_0_xx[7]+0.3535533905932737*F_0_xx[1]*nuVtSq[6]+0.3535533905932737*nuVtSq[2]*F_0_xx[6]+0.3535533905932737*F_0_xx[2]*nuVtSq[5]+0.3535533905932737*nuVtSq[3]*F_0_xx[5]+0.3535533905932737*F_0_xx[3]*nuVtSq[4]; 
  incr_F_0[12] = 0.3535533905932737*nuVtSq[3]*F_0_xx[15]+0.3535533905932737*nuVtSq[5]*F_0_xx[14]+0.3535533905932737*nuVtSq[6]*F_0_xx[13]+0.3535533905932737*nuVtSq[0]*F_0_xx[12]+0.3535533905932737*nuVtSq[7]*F_0_xx[10]+0.3535533905932737*nuVtSq[1]*F_0_xx[9]+0.3535533905932737*nuVtSq[2]*F_0_xx[8]+0.3535533905932737*F_0_xx[4]*nuVtSq[4]; 
  incr_F_0[13] = 0.3535533905932737*nuVtSq[2]*F_0_xx[15]+0.3535533905932737*nuVtSq[4]*F_0_xx[14]+0.3535533905932737*nuVtSq[0]*F_0_xx[13]+0.3535533905932737*nuVtSq[6]*F_0_xx[12]+0.3535533905932737*nuVtSq[1]*F_0_xx[10]+0.3535533905932737*nuVtSq[7]*F_0_xx[9]+0.3535533905932737*nuVtSq[3]*F_0_xx[8]+0.3535533905932737*F_0_xx[4]*nuVtSq[5]; 
  incr_F_0[14] = 0.3535533905932737*nuVtSq[1]*F_0_xx[15]+0.3535533905932737*nuVtSq[0]*F_0_xx[14]+0.3535533905932737*nuVtSq[4]*F_0_xx[13]+0.3535533905932737*nuVtSq[5]*F_0_xx[12]+0.3535533905932737*nuVtSq[2]*F_0_xx[10]+0.3535533905932737*nuVtSq[3]*F_0_xx[9]+0.3535533905932737*nuVtSq[7]*F_0_xx[8]+0.3535533905932737*F_0_xx[4]*nuVtSq[6]; 
  incr_F_0[15] = 0.3535533905932737*nuVtSq[0]*F_0_xx[15]+0.3535533905932737*nuVtSq[1]*F_0_xx[14]+0.3535533905932737*nuVtSq[2]*F_0_xx[13]+0.3535533905932737*nuVtSq[3]*F_0_xx[12]+0.3535533905932737*nuVtSq[4]*F_0_xx[10]+0.3535533905932737*nuVtSq[5]*F_0_xx[9]+0.3535533905932737*nuVtSq[6]*F_0_xx[8]+0.3535533905932737*F_0_xx[4]*nuVtSq[7]; 
  incr_F_0[16] = 0.3535533905932737*nuVtSq[7]*F_0_xx[23]+0.3535533905932737*nuVtSq[6]*F_0_xx[22]+0.3535533905932737*nuVtSq[5]*F_0_xx[21]+0.3535533905932737*nuVtSq[4]*F_0_xx[20]+0.3535533905932737*nuVtSq[3]*F_0_xx[19]+0.3535533905932737*nuVtSq[2]*F_0_xx[18]+0.3535533905932737*nuVtSq[1]*F_0_xx[17]+0.3535533905932737*nuVtSq[0]*F_0_xx[16]; 
  incr_F_0[17] = 0.3535533905932737*nuVtSq[6]*F_0_xx[23]+0.3535533905932737*nuVtSq[7]*F_0_xx[22]+0.3535533905932737*nuVtSq[3]*F_0_xx[21]+0.3535533905932737*nuVtSq[2]*F_0_xx[20]+0.3535533905932737*nuVtSq[5]*F_0_xx[19]+0.3535533905932737*nuVtSq[4]*F_0_xx[18]+0.3535533905932737*nuVtSq[0]*F_0_xx[17]+0.3535533905932737*nuVtSq[1]*F_0_xx[16]; 
  incr_F_0[18] = 0.3535533905932737*nuVtSq[5]*F_0_xx[23]+0.3535533905932737*nuVtSq[3]*F_0_xx[22]+0.3535533905932737*nuVtSq[7]*F_0_xx[21]+0.3535533905932737*nuVtSq[1]*F_0_xx[20]+0.3535533905932737*nuVtSq[6]*F_0_xx[19]+0.3535533905932737*nuVtSq[0]*F_0_xx[18]+0.3535533905932737*nuVtSq[4]*F_0_xx[17]+0.3535533905932737*nuVtSq[2]*F_0_xx[16]; 
  incr_F_0[19] = 0.3535533905932737*nuVtSq[4]*F_0_xx[23]+0.3535533905932737*nuVtSq[2]*F_0_xx[22]+0.3535533905932737*nuVtSq[1]*F_0_xx[21]+0.3535533905932737*nuVtSq[7]*F_0_xx[20]+0.3535533905932737*nuVtSq[0]*F_0_xx[19]+0.3535533905932737*nuVtSq[6]*F_0_xx[18]+0.3535533905932737*nuVtSq[5]*F_0_xx[17]+0.3535533905932737*nuVtSq[3]*F_0_xx[16]; 
  incr_F_0[20] = 0.3535533905932737*nuVtSq[3]*F_0_xx[23]+0.3535533905932737*nuVtSq[5]*F_0_xx[22]+0.3535533905932737*nuVtSq[6]*F_0_xx[21]+0.3535533905932737*nuVtSq[0]*F_0_xx[20]+0.3535533905932737*nuVtSq[7]*F_0_xx[19]+0.3535533905932737*nuVtSq[1]*F_0_xx[18]+0.3535533905932737*nuVtSq[2]*F_0_xx[17]+0.3535533905932737*nuVtSq[4]*F_0_xx[16]; 
  incr_F_0[21] = 0.3535533905932737*nuVtSq[2]*F_0_xx[23]+0.3535533905932737*nuVtSq[4]*F_0_xx[22]+0.3535533905932737*nuVtSq[0]*F_0_xx[21]+0.3535533905932737*nuVtSq[6]*F_0_xx[20]+0.3535533905932737*nuVtSq[1]*F_0_xx[19]+0.3535533905932737*nuVtSq[7]*F_0_xx[18]+0.3535533905932737*nuVtSq[3]*F_0_xx[17]+0.3535533905932737*nuVtSq[5]*F_0_xx[16]; 
  incr_F_0[22] = 0.3535533905932737*nuVtSq[1]*F_0_xx[23]+0.3535533905932737*nuVtSq[0]*F_0_xx[22]+0.3535533905932737*nuVtSq[4]*F_0_xx[21]+0.3535533905932737*nuVtSq[5]*F_0_xx[20]+0.3535533905932737*nuVtSq[2]*F_0_xx[19]+0.3535533905932737*nuVtSq[3]*F_0_xx[18]+0.3535533905932737*nuVtSq[7]*F_0_xx[17]+0.3535533905932737*nuVtSq[6]*F_0_xx[16]; 
  incr_F_0[23] = 0.3535533905932737*nuVtSq[0]*F_0_xx[23]+0.3535533905932737*nuVtSq[1]*F_0_xx[22]+0.3535533905932737*nuVtSq[2]*F_0_xx[21]+0.3535533905932737*nuVtSq[3]*F_0_xx[20]+0.3535533905932737*nuVtSq[4]*F_0_xx[19]+0.3535533905932737*nuVtSq[5]*F_0_xx[18]+0.3535533905932737*nuVtSq[6]*F_0_xx[17]+0.3535533905932737*nuVtSq[7]*F_0_xx[16]; 
  incr_G_1[0] = 0.3535533905932737*nuVtSq[7]*G_1_xx[11]+0.3535533905932737*nuVtSq[6]*G_1_xx[7]+0.3535533905932737*nuVtSq[5]*G_1_xx[6]+0.3535533905932737*nuVtSq[4]*G_1_xx[5]+0.3535533905932737*G_1_xx[3]*nuVtSq[3]+0.3535533905932737*G_1_xx[2]*nuVtSq[2]+0.3535533905932737*G_1_xx[1]*nuVtSq[1]+0.3535533905932737*G_1_xx[0]*nuVtSq[0]; 
  incr_G_1[1] = 0.3535533905932737*nuVtSq[6]*G_1_xx[11]+0.3535533905932737*G_1_xx[7]*nuVtSq[7]+0.3535533905932737*nuVtSq[3]*G_1_xx[6]+0.3535533905932737*G_1_xx[3]*nuVtSq[5]+0.3535533905932737*nuVtSq[2]*G_1_xx[5]+0.3535533905932737*G_1_xx[2]*nuVtSq[4]+0.3535533905932737*G_1_xx[0]*nuVtSq[1]+0.3535533905932737*nuVtSq[0]*G_1_xx[1]; 
  incr_G_1[2] = 0.3535533905932737*nuVtSq[5]*G_1_xx[11]+0.3535533905932737*G_1_xx[6]*nuVtSq[7]+0.3535533905932737*nuVtSq[3]*G_1_xx[7]+0.3535533905932737*G_1_xx[3]*nuVtSq[6]+0.3535533905932737*nuVtSq[1]*G_1_xx[5]+0.3535533905932737*G_1_xx[1]*nuVtSq[4]+0.3535533905932737*G_1_xx[0]*nuVtSq[2]+0.3535533905932737*nuVtSq[0]*G_1_xx[2]; 
  incr_G_1[3] = 0.3535533905932737*nuVtSq[4]*G_1_xx[11]+0.3535533905932737*G_1_xx[5]*nuVtSq[7]+0.3535533905932737*nuVtSq[2]*G_1_xx[7]+0.3535533905932737*G_1_xx[2]*nuVtSq[6]+0.3535533905932737*nuVtSq[1]*G_1_xx[6]+0.3535533905932737*G_1_xx[1]*nuVtSq[5]+0.3535533905932737*G_1_xx[0]*nuVtSq[3]+0.3535533905932737*nuVtSq[0]*G_1_xx[3]; 
  incr_G_1[4] = 0.3535533905932737*nuVtSq[7]*G_1_xx[15]+0.3535533905932737*nuVtSq[6]*G_1_xx[14]+0.3535533905932737*nuVtSq[5]*G_1_xx[13]+0.3535533905932737*nuVtSq[4]*G_1_xx[12]+0.3535533905932737*nuVtSq[3]*G_1_xx[10]+0.3535533905932737*nuVtSq[2]*G_1_xx[9]+0.3535533905932737*nuVtSq[1]*G_1_xx[8]+0.3535533905932737*nuVtSq[0]*G_1_xx[4]; 
  incr_G_1[5] = 0.3535533905932737*nuVtSq[3]*G_1_xx[11]+0.3535533905932737*G_1_xx[3]*nuVtSq[7]+0.3535533905932737*nuVtSq[5]*G_1_xx[7]+0.3535533905932737*G_1_xx[6]*nuVtSq[6]+0.3535533905932737*nuVtSq[0]*G_1_xx[5]+0.3535533905932737*G_1_xx[0]*nuVtSq[4]+0.3535533905932737*G_1_xx[1]*nuVtSq[2]+0.3535533905932737*nuVtSq[1]*G_1_xx[2]; 
  incr_G_1[6] = 0.3535533905932737*nuVtSq[2]*G_1_xx[11]+0.3535533905932737*G_1_xx[2]*nuVtSq[7]+0.3535533905932737*nuVtSq[4]*G_1_xx[7]+0.3535533905932737*G_1_xx[5]*nuVtSq[6]+0.3535533905932737*nuVtSq[0]*G_1_xx[6]+0.3535533905932737*G_1_xx[0]*nuVtSq[5]+0.3535533905932737*G_1_xx[1]*nuVtSq[3]+0.3535533905932737*nuVtSq[1]*G_1_xx[3]; 
  incr_G_1[7] = 0.3535533905932737*nuVtSq[1]*G_1_xx[11]+0.3535533905932737*G_1_xx[1]*nuVtSq[7]+0.3535533905932737*nuVtSq[0]*G_1_xx[7]+0.3535533905932737*G_1_xx[0]*nuVtSq[6]+0.3535533905932737*nuVtSq[4]*G_1_xx[6]+0.3535533905932737*G_1_xx[5]*nuVtSq[5]+0.3535533905932737*G_1_xx[2]*nuVtSq[3]+0.3535533905932737*nuVtSq[2]*G_1_xx[3]; 
  incr_G_1[8] = 0.3535533905932737*nuVtSq[6]*G_1_xx[15]+0.3535533905932737*nuVtSq[7]*G_1_xx[14]+0.3535533905932737*nuVtSq[3]*G_1_xx[13]+0.3535533905932737*nuVtSq[2]*G_1_xx[12]+0.3535533905932737*nuVtSq[5]*G_1_xx[10]+0.3535533905932737*nuVtSq[4]*G_1_xx[9]+0.3535533905932737*nuVtSq[0]*G_1_xx[8]+0.3535533905932737*nuVtSq[1]*G_1_xx[4]; 
  incr_G_1[9] = 0.3535533905932737*nuVtSq[5]*G_1_xx[15]+0.3535533905932737*nuVtSq[3]*G_1_xx[14]+0.3535533905932737*nuVtSq[7]*G_1_xx[13]+0.3535533905932737*nuVtSq[1]*G_1_xx[12]+0.3535533905932737*nuVtSq[6]*G_1_xx[10]+0.3535533905932737*nuVtSq[0]*G_1_xx[9]+0.3535533905932737*nuVtSq[4]*G_1_xx[8]+0.3535533905932737*nuVtSq[2]*G_1_xx[4]; 
  incr_G_1[10] = 0.3535533905932737*nuVtSq[4]*G_1_xx[15]+0.3535533905932737*nuVtSq[2]*G_1_xx[14]+0.3535533905932737*nuVtSq[1]*G_1_xx[13]+0.3535533905932737*nuVtSq[7]*G_1_xx[12]+0.3535533905932737*nuVtSq[0]*G_1_xx[10]+0.3535533905932737*nuVtSq[6]*G_1_xx[9]+0.3535533905932737*nuVtSq[5]*G_1_xx[8]+0.3535533905932737*nuVtSq[3]*G_1_xx[4]; 
  incr_G_1[11] = 0.3535533905932737*nuVtSq[0]*G_1_xx[11]+0.3535533905932737*G_1_xx[0]*nuVtSq[7]+0.3535533905932737*nuVtSq[1]*G_1_xx[7]+0.3535533905932737*G_1_xx[1]*nuVtSq[6]+0.3535533905932737*nuVtSq[2]*G_1_xx[6]+0.3535533905932737*G_1_xx[2]*nuVtSq[5]+0.3535533905932737*nuVtSq[3]*G_1_xx[5]+0.3535533905932737*G_1_xx[3]*nuVtSq[4]; 
  incr_G_1[12] = 0.3535533905932737*nuVtSq[3]*G_1_xx[15]+0.3535533905932737*nuVtSq[5]*G_1_xx[14]+0.3535533905932737*nuVtSq[6]*G_1_xx[13]+0.3535533905932737*nuVtSq[0]*G_1_xx[12]+0.3535533905932737*nuVtSq[7]*G_1_xx[10]+0.3535533905932737*nuVtSq[1]*G_1_xx[9]+0.3535533905932737*nuVtSq[2]*G_1_xx[8]+0.3535533905932737*G_1_xx[4]*nuVtSq[4]; 
  incr_G_1[13] = 0.3535533905932737*nuVtSq[2]*G_1_xx[15]+0.3535533905932737*nuVtSq[4]*G_1_xx[14]+0.3535533905932737*nuVtSq[0]*G_1_xx[13]+0.3535533905932737*nuVtSq[6]*G_1_xx[12]+0.3535533905932737*nuVtSq[1]*G_1_xx[10]+0.3535533905932737*nuVtSq[7]*G_1_xx[9]+0.3535533905932737*nuVtSq[3]*G_1_xx[8]+0.3535533905932737*G_1_xx[4]*nuVtSq[5]; 
  incr_G_1[14] = 0.3535533905932737*nuVtSq[1]*G_1_xx[15]+0.3535533905932737*nuVtSq[0]*G_1_xx[14]+0.3535533905932737*nuVtSq[4]*G_1_xx[13]+0.3535533905932737*nuVtSq[5]*G_1_xx[12]+0.3535533905932737*nuVtSq[2]*G_1_xx[10]+0.3535533905932737*nuVtSq[3]*G_1_xx[9]+0.3535533905932737*nuVtSq[7]*G_1_xx[8]+0.3535533905932737*G_1_xx[4]*nuVtSq[6]; 
  incr_G_1[15] = 0.3535533905932737*nuVtSq[0]*G_1_xx[15]+0.3535533905932737*nuVtSq[1]*G_1_xx[14]+0.3535533905932737*nuVtSq[2]*G_1_xx[13]+0.3535533905932737*nuVtSq[3]*G_1_xx[12]+0.3535533905932737*nuVtSq[4]*G_1_xx[10]+0.3535533905932737*nuVtSq[5]*G_1_xx[9]+0.3535533905932737*nuVtSq[6]*G_1_xx[8]+0.3535533905932737*G_1_xx[4]*nuVtSq[7]; 
  incr_G_1[16] = 0.3535533905932737*nuVtSq[7]*G_1_xx[23]+0.3535533905932737*nuVtSq[6]*G_1_xx[22]+0.3535533905932737*nuVtSq[5]*G_1_xx[21]+0.3535533905932737*nuVtSq[4]*G_1_xx[20]+0.3535533905932737*nuVtSq[3]*G_1_xx[19]+0.3535533905932737*nuVtSq[2]*G_1_xx[18]+0.3535533905932737*nuVtSq[1]*G_1_xx[17]+0.3535533905932737*nuVtSq[0]*G_1_xx[16]; 
  incr_G_1[17] = 0.3535533905932737*nuVtSq[6]*G_1_xx[23]+0.3535533905932737*nuVtSq[7]*G_1_xx[22]+0.3535533905932737*nuVtSq[3]*G_1_xx[21]+0.3535533905932737*nuVtSq[2]*G_1_xx[20]+0.3535533905932737*nuVtSq[5]*G_1_xx[19]+0.3535533905932737*nuVtSq[4]*G_1_xx[18]+0.3535533905932737*nuVtSq[0]*G_1_xx[17]+0.3535533905932737*nuVtSq[1]*G_1_xx[16]; 
  incr_G_1[18] = 0.3535533905932737*nuVtSq[5]*G_1_xx[23]+0.3535533905932737*nuVtSq[3]*G_1_xx[22]+0.3535533905932737*nuVtSq[7]*G_1_xx[21]+0.3535533905932737*nuVtSq[1]*G_1_xx[20]+0.3535533905932737*nuVtSq[6]*G_1_xx[19]+0.3535533905932737*nuVtSq[0]*G_1_xx[18]+0.3535533905932737*nuVtSq[4]*G_1_xx[17]+0.3535533905932737*nuVtSq[2]*G_1_xx[16]; 
  incr_G_1[19] = 0.3535533905932737*nuVtSq[4]*G_1_xx[23]+0.3535533905932737*nuVtSq[2]*G_1_xx[22]+0.3535533905932737*nuVtSq[1]*G_1_xx[21]+0.3535533905932737*nuVtSq[7]*G_1_xx[20]+0.3535533905932737*nuVtSq[0]*G_1_xx[19]+0.3535533905932737*nuVtSq[6]*G_1_xx[18]+0.3535533905932737*nuVtSq[5]*G_1_xx[17]+0.3535533905932737*nuVtSq[3]*G_1_xx[16]; 
  incr_G_1[20] = 0.3535533905932737*nuVtSq[3]*G_1_xx[23]+0.3535533905932737*nuVtSq[5]*G_1_xx[22]+0.3535533905932737*nuVtSq[6]*G_1_xx[21]+0.3535533905932737*nuVtSq[0]*G_1_xx[20]+0.3535533905932737*nuVtSq[7]*G_1_xx[19]+0.3535533905932737*nuVtSq[1]*G_1_xx[18]+0.3535533905932737*nuVtSq[2]*G_1_xx[17]+0.3535533905932737*nuVtSq[4]*G_1_xx[16]; 
  incr_G_1[21] = 0.3535533905932737*nuVtSq[2]*G_1_xx[23]+0.3535533905932737*nuVtSq[4]*G_1_xx[22]+0.3535533905932737*nuVtSq[0]*G_1_xx[21]+0.3535533905932737*nuVtSq[6]*G_1_xx[20]+0.3535533905932737*nuVtSq[1]*G_1_xx[19]+0.3535533905932737*nuVtSq[7]*G_1_xx[18]+0.3535533905932737*nuVtSq[3]*G_1_xx[17]+0.3535533905932737*nuVtSq[5]*G_1_xx[16]; 
  incr_G_1[22] = 0.3535533905932737*nuVtSq[1]*G_1_xx[23]+0.3535533905932737*nuVtSq[0]*G_1_xx[22]+0.3535533905932737*nuVtSq[4]*G_1_xx[21]+0.3535533905932737*nuVtSq[5]*G_1_xx[20]+0.3535533905932737*nuVtSq[2]*G_1_xx[19]+0.3535533905932737*nuVtSq[3]*G_1_xx[18]+0.3535533905932737*nuVtSq[7]*G_1_xx[17]+0.3535533905932737*nuVtSq[6]*G_1_xx[16]; 
  incr_G_1[23] = 0.3535533905932737*nuVtSq[0]*G_1_xx[23]+0.3535533905932737*nuVtSq[1]*G_1_xx[22]+0.3535533905932737*nuVtSq[2]*G_1_xx[21]+0.3535533905932737*nuVtSq[3]*G_1_xx[20]+0.3535533905932737*nuVtSq[4]*G_1_xx[19]+0.3535533905932737*nuVtSq[5]*G_1_xx[18]+0.3535533905932737*nuVtSq[6]*G_1_xx[17]+0.3535533905932737*nuVtSq[7]*G_1_xx[16]; 

  out_F_0[0] += incr_F_0[0]*rdvSq4; 
  out_F_0[1] += incr_F_0[1]*rdvSq4; 
  out_F_0[2] += incr_F_0[2]*rdvSq4; 
  out_F_0[3] += incr_F_0[3]*rdvSq4; 
  out_F_0[4] += incr_F_0[4]*rdvSq4; 
  out_F_0[5] += incr_F_0[5]*rdvSq4; 
  out_F_0[6] += incr_F_0[6]*rdvSq4; 
  out_F_0[7] += incr_F_0[7]*rdvSq4; 
  out_F_0[8] += incr_F_0[8]*rdvSq4; 
  out_F_0[9] += incr_F_0[9]*rdvSq4; 
  out_F_0[10] += incr_F_0[10]*rdvSq4; 
  out_F_0[11] += incr_F_0[11]*rdvSq4; 
  out_F_0[12] += incr_F_0[12]*rdvSq4; 
  out_F_0[13] += incr_F_0[13]*rdvSq4; 
  out_F_0[14] += incr_F_0[14]*rdvSq4; 
  out_F_0[15] += incr_F_0[15]*rdvSq4; 
  out_F_0[16] += incr_F_0[16]*rdvSq4; 
  out_F_0[17] += incr_F_0[17]*rdvSq4; 
  out_F_0[18] += incr_F_0[18]*rdvSq4; 
  out_F_0[19] += incr_F_0[19]*rdvSq4; 
  out_F_0[20] += incr_F_0[20]*rdvSq4; 
  out_F_0[21] += incr_F_0[21]*rdvSq4; 
  out_F_0[22] += incr_F_0[22]*rdvSq4; 
  out_F_0[23] += incr_F_0[23]*rdvSq4; 
  out_G_1[0] += incr_G_1[0]*rdvSq4; 
  out_G_1[1] += incr_G_1[1]*rdvSq4; 
  out_G_1[2] += incr_G_1[2]*rdvSq4; 
  out_G_1[3] += incr_G_1[3]*rdvSq4; 
  out_G_1[4] += incr_G_1[4]*rdvSq4; 
  out_G_1[5] += incr_G_1[5]*rdvSq4; 
  out_G_1[6] += incr_G_1[6]*rdvSq4; 
  out_G_1[7] += incr_G_1[7]*rdvSq4; 
  out_G_1[8] += incr_G_1[8]*rdvSq4; 
  out_G_1[9] += incr_G_1[9]*rdvSq4; 
  out_G_1[10] += incr_G_1[10]*rdvSq4; 
  out_G_1[11] += incr_G_1[11]*rdvSq4; 
  out_G_1[12] += incr_G_1[12]*rdvSq4; 
  out_G_1[13] += incr_G_1[13]*rdvSq4; 
  out_G_1[14] += incr_G_1[14]*rdvSq4; 
  out_G_1[15] += incr_G_1[15]*rdvSq4; 
  out_G_1[16] += incr_G_1[16]*rdvSq4; 
  out_G_1[17] += incr_G_1[17]*rdvSq4; 
  out_G_1[18] += incr_G_1[18]*rdvSq4; 
  out_G_1[19] += incr_G_1[19]*rdvSq4; 
  out_G_1[20] += incr_G_1[20]*rdvSq4; 
  out_G_1[21] += incr_G_1[21]*rdvSq4; 
  out_G_1[22] += incr_G_1[22]*rdvSq4; 
  out_G_1[23] += incr_G_1[23]*rdvSq4; 

  return 0.;

} 
