#include <gkyl_dg_diffusion_kernels.h> 
GKYL_CU_DH void dg_diffusion_iso_euler_surfy_2x_tensor_p2(const double* w, const double* dx, double D, 
  const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: Cell-center coordinates
  // dxv[NDIM]: Cell spacing
  // D: Diffusion coefficient in the center cell
  // ql: Input field in the left cell
  // qc: Input field in the center cell
  // qr: Input field in the right cell
  // out: Incremented output

  const double dx1 = 2.0/dx[1]; 
  const double J = pow(dx1, 2.0);

  const double *q0l = &ql[0]; 
  const double *q0c = &qc[0]; 
  const double *q0r = &qr[0]; 
  double *out0= &out[0]; 

  out0[0] += J*D*(0.6708203932499369*q0r[5]+0.6708203932499369*q0l[5]-1.341640786499874*q0c[5]-1.190784930203603*q0r[2]+1.190784930203603*q0l[2]+0.9375*q0r[0]+0.9375*q0l[0]-1.875*q0c[0]); 
  out0[1] += J*D*(0.6708203932499369*q0r[7]+0.6708203932499369*q0l[7]-1.341640786499874*q0c[7]-1.190784930203603*q0r[3]+1.190784930203603*q0l[3]+0.9375*q0r[1]+0.9375*q0l[1]-1.875*q0c[1]); 
  out0[2] += J*D*(0.7382874503707888*q0r[5]-0.7382874503707888*q0l[5]-1.453125*q0r[2]-1.453125*q0l[2]-5.34375*q0c[2]+1.190784930203603*q0r[0]-1.190784930203603*q0l[0]); 
  out0[3] += J*D*(0.7382874503707888*q0r[7]-0.7382874503707888*q0l[7]-1.453125*q0r[3]-1.453125*q0l[3]-5.34375*q0c[3]+1.190784930203603*q0r[1]-1.190784930203603*q0l[1]); 
  out0[4] += J*D*(0.6708203932499369*q0r[8]+0.6708203932499369*q0l[8]-1.341640786499874*q0c[8]-1.190784930203603*q0r[6]+1.190784930203603*q0l[6]+0.9375*q0r[4]+0.9375*q0l[4]-1.875*q0c[4]); 
  out0[5] += J*D*((-0.140625*q0r[5])-0.140625*q0l[5]-6.28125*q0c[5]-0.3025768239224545*q0r[2]+0.3025768239224545*q0l[2]+0.4192627457812106*q0r[0]+0.4192627457812106*q0l[0]-0.8385254915624212*q0c[0]); 
  out0[6] += J*D*(0.7382874503707888*q0r[8]-0.7382874503707888*q0l[8]-1.453125*q0r[6]-1.453125*q0l[6]-5.34375*q0c[6]+1.190784930203603*q0r[4]-1.190784930203603*q0l[4]); 
  out0[7] += J*D*((-0.140625*q0r[7])-0.140625*q0l[7]-6.28125*q0c[7]-0.3025768239224544*q0r[3]+0.3025768239224544*q0l[3]+0.4192627457812105*q0r[1]+0.4192627457812105*q0l[1]-0.8385254915624211*q0c[1]); 
  out0[8] += J*D*((-0.140625*q0r[8])-0.140625*q0l[8]-6.28125*q0c[8]-0.3025768239224544*q0r[6]+0.3025768239224544*q0l[6]+0.4192627457812106*q0r[4]+0.4192627457812106*q0l[4]-0.8385254915624212*q0c[4]); 

  const double *q1l = &ql[9]; 
  const double *q1c = &qc[9]; 
  const double *q1r = &qr[9]; 
  double *out1= &out[9]; 

  out1[0] += J*D*(0.6708203932499369*q1r[5]+0.6708203932499369*q1l[5]-1.341640786499874*q1c[5]-1.190784930203603*q1r[2]+1.190784930203603*q1l[2]+0.9375*q1r[0]+0.9375*q1l[0]-1.875*q1c[0]); 
  out1[1] += J*D*(0.6708203932499369*q1r[7]+0.6708203932499369*q1l[7]-1.341640786499874*q1c[7]-1.190784930203603*q1r[3]+1.190784930203603*q1l[3]+0.9375*q1r[1]+0.9375*q1l[1]-1.875*q1c[1]); 
  out1[2] += J*D*(0.7382874503707888*q1r[5]-0.7382874503707888*q1l[5]-1.453125*q1r[2]-1.453125*q1l[2]-5.34375*q1c[2]+1.190784930203603*q1r[0]-1.190784930203603*q1l[0]); 
  out1[3] += J*D*(0.7382874503707888*q1r[7]-0.7382874503707888*q1l[7]-1.453125*q1r[3]-1.453125*q1l[3]-5.34375*q1c[3]+1.190784930203603*q1r[1]-1.190784930203603*q1l[1]); 
  out1[4] += J*D*(0.6708203932499369*q1r[8]+0.6708203932499369*q1l[8]-1.341640786499874*q1c[8]-1.190784930203603*q1r[6]+1.190784930203603*q1l[6]+0.9375*q1r[4]+0.9375*q1l[4]-1.875*q1c[4]); 
  out1[5] += J*D*((-0.140625*q1r[5])-0.140625*q1l[5]-6.28125*q1c[5]-0.3025768239224545*q1r[2]+0.3025768239224545*q1l[2]+0.4192627457812106*q1r[0]+0.4192627457812106*q1l[0]-0.8385254915624212*q1c[0]); 
  out1[6] += J*D*(0.7382874503707888*q1r[8]-0.7382874503707888*q1l[8]-1.453125*q1r[6]-1.453125*q1l[6]-5.34375*q1c[6]+1.190784930203603*q1r[4]-1.190784930203603*q1l[4]); 
  out1[7] += J*D*((-0.140625*q1r[7])-0.140625*q1l[7]-6.28125*q1c[7]-0.3025768239224544*q1r[3]+0.3025768239224544*q1l[3]+0.4192627457812105*q1r[1]+0.4192627457812105*q1l[1]-0.8385254915624211*q1c[1]); 
  out1[8] += J*D*((-0.140625*q1r[8])-0.140625*q1l[8]-6.28125*q1c[8]-0.3025768239224544*q1r[6]+0.3025768239224544*q1l[6]+0.4192627457812106*q1r[4]+0.4192627457812106*q1l[4]-0.8385254915624212*q1c[4]); 

  const double *q2l = &ql[18]; 
  const double *q2c = &qc[18]; 
  const double *q2r = &qr[18]; 
  double *out2= &out[18]; 

  out2[0] += J*D*(0.6708203932499369*q2r[5]+0.6708203932499369*q2l[5]-1.341640786499874*q2c[5]-1.190784930203603*q2r[2]+1.190784930203603*q2l[2]+0.9375*q2r[0]+0.9375*q2l[0]-1.875*q2c[0]); 
  out2[1] += J*D*(0.6708203932499369*q2r[7]+0.6708203932499369*q2l[7]-1.341640786499874*q2c[7]-1.190784930203603*q2r[3]+1.190784930203603*q2l[3]+0.9375*q2r[1]+0.9375*q2l[1]-1.875*q2c[1]); 
  out2[2] += J*D*(0.7382874503707888*q2r[5]-0.7382874503707888*q2l[5]-1.453125*q2r[2]-1.453125*q2l[2]-5.34375*q2c[2]+1.190784930203603*q2r[0]-1.190784930203603*q2l[0]); 
  out2[3] += J*D*(0.7382874503707888*q2r[7]-0.7382874503707888*q2l[7]-1.453125*q2r[3]-1.453125*q2l[3]-5.34375*q2c[3]+1.190784930203603*q2r[1]-1.190784930203603*q2l[1]); 
  out2[4] += J*D*(0.6708203932499369*q2r[8]+0.6708203932499369*q2l[8]-1.341640786499874*q2c[8]-1.190784930203603*q2r[6]+1.190784930203603*q2l[6]+0.9375*q2r[4]+0.9375*q2l[4]-1.875*q2c[4]); 
  out2[5] += J*D*((-0.140625*q2r[5])-0.140625*q2l[5]-6.28125*q2c[5]-0.3025768239224545*q2r[2]+0.3025768239224545*q2l[2]+0.4192627457812106*q2r[0]+0.4192627457812106*q2l[0]-0.8385254915624212*q2c[0]); 
  out2[6] += J*D*(0.7382874503707888*q2r[8]-0.7382874503707888*q2l[8]-1.453125*q2r[6]-1.453125*q2l[6]-5.34375*q2c[6]+1.190784930203603*q2r[4]-1.190784930203603*q2l[4]); 
  out2[7] += J*D*((-0.140625*q2r[7])-0.140625*q2l[7]-6.28125*q2c[7]-0.3025768239224544*q2r[3]+0.3025768239224544*q2l[3]+0.4192627457812105*q2r[1]+0.4192627457812105*q2l[1]-0.8385254915624211*q2c[1]); 
  out2[8] += J*D*((-0.140625*q2r[8])-0.140625*q2l[8]-6.28125*q2c[8]-0.3025768239224544*q2r[6]+0.3025768239224544*q2l[6]+0.4192627457812106*q2r[4]+0.4192627457812106*q2l[4]-0.8385254915624212*q2c[4]); 

  const double *q3l = &ql[27]; 
  const double *q3c = &qc[27]; 
  const double *q3r = &qr[27]; 
  double *out3= &out[27]; 

  out3[0] += J*D*(0.6708203932499369*q3r[5]+0.6708203932499369*q3l[5]-1.341640786499874*q3c[5]-1.190784930203603*q3r[2]+1.190784930203603*q3l[2]+0.9375*q3r[0]+0.9375*q3l[0]-1.875*q3c[0]); 
  out3[1] += J*D*(0.6708203932499369*q3r[7]+0.6708203932499369*q3l[7]-1.341640786499874*q3c[7]-1.190784930203603*q3r[3]+1.190784930203603*q3l[3]+0.9375*q3r[1]+0.9375*q3l[1]-1.875*q3c[1]); 
  out3[2] += J*D*(0.7382874503707888*q3r[5]-0.7382874503707888*q3l[5]-1.453125*q3r[2]-1.453125*q3l[2]-5.34375*q3c[2]+1.190784930203603*q3r[0]-1.190784930203603*q3l[0]); 
  out3[3] += J*D*(0.7382874503707888*q3r[7]-0.7382874503707888*q3l[7]-1.453125*q3r[3]-1.453125*q3l[3]-5.34375*q3c[3]+1.190784930203603*q3r[1]-1.190784930203603*q3l[1]); 
  out3[4] += J*D*(0.6708203932499369*q3r[8]+0.6708203932499369*q3l[8]-1.341640786499874*q3c[8]-1.190784930203603*q3r[6]+1.190784930203603*q3l[6]+0.9375*q3r[4]+0.9375*q3l[4]-1.875*q3c[4]); 
  out3[5] += J*D*((-0.140625*q3r[5])-0.140625*q3l[5]-6.28125*q3c[5]-0.3025768239224545*q3r[2]+0.3025768239224545*q3l[2]+0.4192627457812106*q3r[0]+0.4192627457812106*q3l[0]-0.8385254915624212*q3c[0]); 
  out3[6] += J*D*(0.7382874503707888*q3r[8]-0.7382874503707888*q3l[8]-1.453125*q3r[6]-1.453125*q3l[6]-5.34375*q3c[6]+1.190784930203603*q3r[4]-1.190784930203603*q3l[4]); 
  out3[7] += J*D*((-0.140625*q3r[7])-0.140625*q3l[7]-6.28125*q3c[7]-0.3025768239224544*q3r[3]+0.3025768239224544*q3l[3]+0.4192627457812105*q3r[1]+0.4192627457812105*q3l[1]-0.8385254915624211*q3c[1]); 
  out3[8] += J*D*((-0.140625*q3r[8])-0.140625*q3l[8]-6.28125*q3c[8]-0.3025768239224544*q3r[6]+0.3025768239224544*q3l[6]+0.4192627457812106*q3r[4]+0.4192627457812106*q3l[4]-0.8385254915624212*q3c[4]); 

} 
