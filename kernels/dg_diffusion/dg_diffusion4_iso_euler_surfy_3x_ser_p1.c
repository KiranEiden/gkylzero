#include <gkyl_dg_diffusion_kernels.h> 
GKYL_CU_DH void dg_diffusion4_iso_euler_surfy_3x_ser_p1(const double* w, const double* dx, double D, 
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
  const double J = -1.0*pow(dx1, 4.0);

  const double *q0l = &ql[0]; 
  const double *q0c = &qc[0]; 
  const double *q0r = &qr[0]; 
  double *out0= &out[0]; 

  out0[0] += J*D*(1.623797632095822*q0r[2]-1.623797632095822*q0l[2]-0.9375*q0r[0]-0.9375*q0l[0]+1.875*q0c[0]); 
  out0[1] += J*D*(1.623797632095822*q0r[4]-1.623797632095822*q0l[4]-0.9375*q0r[1]-0.9375*q0l[1]+1.875*q0c[1]); 
  out0[2] += J*D*(2.0625*q0r[2]+2.0625*q0l[2]+7.125*q0c[2]-1.623797632095822*q0r[0]+1.623797632095822*q0l[0]); 
  out0[3] += J*D*(1.623797632095822*q0r[6]-1.623797632095822*q0l[6]-0.9375*q0r[3]-0.9375*q0l[3]+1.875*q0c[3]); 
  out0[4] += J*D*(2.0625*q0r[4]+2.0625*q0l[4]+7.125*q0c[4]-1.623797632095822*q0r[1]+1.623797632095822*q0l[1]); 
  out0[5] += J*D*(1.623797632095822*q0r[7]-1.623797632095822*q0l[7]-0.9375*q0r[5]-0.9375*q0l[5]+1.875*q0c[5]); 
  out0[6] += J*D*(2.0625*q0r[6]+2.0625*q0l[6]+7.125*q0c[6]-1.623797632095822*q0r[3]+1.623797632095822*q0l[3]); 
  out0[7] += J*D*(2.0625*q0r[7]+2.0625*q0l[7]+7.125*q0c[7]-1.623797632095822*q0r[5]+1.623797632095822*q0l[5]); 

  const double *q1l = &ql[8]; 
  const double *q1c = &qc[8]; 
  const double *q1r = &qr[8]; 
  double *out1= &out[8]; 

  out1[0] += J*D*(1.623797632095822*q1r[2]-1.623797632095822*q1l[2]-0.9375*q1r[0]-0.9375*q1l[0]+1.875*q1c[0]); 
  out1[1] += J*D*(1.623797632095822*q1r[4]-1.623797632095822*q1l[4]-0.9375*q1r[1]-0.9375*q1l[1]+1.875*q1c[1]); 
  out1[2] += J*D*(2.0625*q1r[2]+2.0625*q1l[2]+7.125*q1c[2]-1.623797632095822*q1r[0]+1.623797632095822*q1l[0]); 
  out1[3] += J*D*(1.623797632095822*q1r[6]-1.623797632095822*q1l[6]-0.9375*q1r[3]-0.9375*q1l[3]+1.875*q1c[3]); 
  out1[4] += J*D*(2.0625*q1r[4]+2.0625*q1l[4]+7.125*q1c[4]-1.623797632095822*q1r[1]+1.623797632095822*q1l[1]); 
  out1[5] += J*D*(1.623797632095822*q1r[7]-1.623797632095822*q1l[7]-0.9375*q1r[5]-0.9375*q1l[5]+1.875*q1c[5]); 
  out1[6] += J*D*(2.0625*q1r[6]+2.0625*q1l[6]+7.125*q1c[6]-1.623797632095822*q1r[3]+1.623797632095822*q1l[3]); 
  out1[7] += J*D*(2.0625*q1r[7]+2.0625*q1l[7]+7.125*q1c[7]-1.623797632095822*q1r[5]+1.623797632095822*q1l[5]); 

  const double *q2l = &ql[16]; 
  const double *q2c = &qc[16]; 
  const double *q2r = &qr[16]; 
  double *out2= &out[16]; 

  out2[0] += J*D*(1.623797632095822*q2r[2]-1.623797632095822*q2l[2]-0.9375*q2r[0]-0.9375*q2l[0]+1.875*q2c[0]); 
  out2[1] += J*D*(1.623797632095822*q2r[4]-1.623797632095822*q2l[4]-0.9375*q2r[1]-0.9375*q2l[1]+1.875*q2c[1]); 
  out2[2] += J*D*(2.0625*q2r[2]+2.0625*q2l[2]+7.125*q2c[2]-1.623797632095822*q2r[0]+1.623797632095822*q2l[0]); 
  out2[3] += J*D*(1.623797632095822*q2r[6]-1.623797632095822*q2l[6]-0.9375*q2r[3]-0.9375*q2l[3]+1.875*q2c[3]); 
  out2[4] += J*D*(2.0625*q2r[4]+2.0625*q2l[4]+7.125*q2c[4]-1.623797632095822*q2r[1]+1.623797632095822*q2l[1]); 
  out2[5] += J*D*(1.623797632095822*q2r[7]-1.623797632095822*q2l[7]-0.9375*q2r[5]-0.9375*q2l[5]+1.875*q2c[5]); 
  out2[6] += J*D*(2.0625*q2r[6]+2.0625*q2l[6]+7.125*q2c[6]-1.623797632095822*q2r[3]+1.623797632095822*q2l[3]); 
  out2[7] += J*D*(2.0625*q2r[7]+2.0625*q2l[7]+7.125*q2c[7]-1.623797632095822*q2r[5]+1.623797632095822*q2l[5]); 

  const double *q3l = &ql[24]; 
  const double *q3c = &qc[24]; 
  const double *q3r = &qr[24]; 
  double *out3= &out[24]; 

  out3[0] += J*D*(1.623797632095822*q3r[2]-1.623797632095822*q3l[2]-0.9375*q3r[0]-0.9375*q3l[0]+1.875*q3c[0]); 
  out3[1] += J*D*(1.623797632095822*q3r[4]-1.623797632095822*q3l[4]-0.9375*q3r[1]-0.9375*q3l[1]+1.875*q3c[1]); 
  out3[2] += J*D*(2.0625*q3r[2]+2.0625*q3l[2]+7.125*q3c[2]-1.623797632095822*q3r[0]+1.623797632095822*q3l[0]); 
  out3[3] += J*D*(1.623797632095822*q3r[6]-1.623797632095822*q3l[6]-0.9375*q3r[3]-0.9375*q3l[3]+1.875*q3c[3]); 
  out3[4] += J*D*(2.0625*q3r[4]+2.0625*q3l[4]+7.125*q3c[4]-1.623797632095822*q3r[1]+1.623797632095822*q3l[1]); 
  out3[5] += J*D*(1.623797632095822*q3r[7]-1.623797632095822*q3l[7]-0.9375*q3r[5]-0.9375*q3l[5]+1.875*q3c[5]); 
  out3[6] += J*D*(2.0625*q3r[6]+2.0625*q3l[6]+7.125*q3c[6]-1.623797632095822*q3r[3]+1.623797632095822*q3l[3]); 
  out3[7] += J*D*(2.0625*q3r[7]+2.0625*q3l[7]+7.125*q3c[7]-1.623797632095822*q3r[5]+1.623797632095822*q3l[5]); 

} 
