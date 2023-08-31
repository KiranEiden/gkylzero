#include <gkyl_dg_diffusion_kernels.h> 
GKYL_CU_DH double dg_diffusion4_pkpm_surfx_2x_ser_p1(const double* w, const double* dx, double D, 
  const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: Cell-center coordinates
  // dxv[NDIM]: Cell spacing
  // D: Diffusion coefficient in the center cell
  // ql: Input field in the left cell
  // qc: Input field in the center cell
  // qr: Input field in the right cell
  // out: Incremented output

  const double dx1 = 2.0/dx[0]; 
  const double J = -1.0*pow(dx1, 4.0);

  const double *q0l = &ql[0]; 
  const double *q0c = &qc[0]; 
  const double *q0r = &qr[0]; 
  double *out0 = &out[0]; 

  out0[0] += J*D*(1.623797632095822*q0r[1]-1.623797632095822*q0l[1]-0.9375*q0r[0]-0.9375*q0l[0]+1.875*q0c[0]); 
  out0[1] += J*D*(2.0625*q0r[1]+2.0625*q0l[1]+7.125*q0c[1]-1.623797632095822*q0r[0]+1.623797632095822*q0l[0]); 
  out0[2] += J*D*(1.623797632095822*q0r[3]-1.623797632095822*q0l[3]-0.9375*q0r[2]-0.9375*q0l[2]+1.875*q0c[2]); 
  out0[3] += J*D*(2.0625*q0r[3]+2.0625*q0l[3]+7.125*q0c[3]-1.623797632095822*q0r[2]+1.623797632095822*q0l[2]); 

  const double *q1l = &ql[4]; 
  const double *q1c = &qc[4]; 
  const double *q1r = &qr[4]; 
  double *out1 = &out[4]; 

  out1[0] += J*D*(1.623797632095822*q1r[1]-1.623797632095822*q1l[1]-0.9375*q1r[0]-0.9375*q1l[0]+1.875*q1c[0]); 
  out1[1] += J*D*(2.0625*q1r[1]+2.0625*q1l[1]+7.125*q1c[1]-1.623797632095822*q1r[0]+1.623797632095822*q1l[0]); 
  out1[2] += J*D*(1.623797632095822*q1r[3]-1.623797632095822*q1l[3]-0.9375*q1r[2]-0.9375*q1l[2]+1.875*q1c[2]); 
  out1[3] += J*D*(2.0625*q1r[3]+2.0625*q1l[3]+7.125*q1c[3]-1.623797632095822*q1r[2]+1.623797632095822*q1l[2]); 

  const double *q2l = &ql[8]; 
  const double *q2c = &qc[8]; 
  const double *q2r = &qr[8]; 
  double *out2 = &out[8]; 

  out2[0] += J*D*(1.623797632095822*q2r[1]-1.623797632095822*q2l[1]-0.9375*q2r[0]-0.9375*q2l[0]+1.875*q2c[0]); 
  out2[1] += J*D*(2.0625*q2r[1]+2.0625*q2l[1]+7.125*q2c[1]-1.623797632095822*q2r[0]+1.623797632095822*q2l[0]); 
  out2[2] += J*D*(1.623797632095822*q2r[3]-1.623797632095822*q2l[3]-0.9375*q2r[2]-0.9375*q2l[2]+1.875*q2c[2]); 
  out2[3] += J*D*(2.0625*q2r[3]+2.0625*q2l[3]+7.125*q2c[3]-1.623797632095822*q2r[2]+1.623797632095822*q2l[2]); 

  return 0.;

} 
