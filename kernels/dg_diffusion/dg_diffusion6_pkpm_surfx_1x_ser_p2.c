#include <gkyl_dg_diffusion_kernels.h> 
GKYL_CU_DH double dg_diffusion6_pkpm_surfx_1x_ser_p2(const double* w, const double* dx, double D, 
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
  const double J = pow(dx1, 6.0);

  const double *q0l = &ql[0]; 
  const double *q0c = &qc[0]; 
  const double *q0r = &qr[0]; 
  double *out0 = &out[0]; 

  out0[0] += J*D*(35.21807064562169*q0r[2]+35.21807064562169*q0l[2]-70.43614129124337*q0c[2]-34.09975027401226*q0r[1]+34.09975027401226*q0l[1]+19.6875*q0r[0]+19.6875*q0l[0]-39.375*q0c[0]); 
  out0[1] += J*D*(51.46831774920947*q0r[2]-51.46831774920947*q0l[2]-56.6015625*q0r[1]-56.6015625*q0l[1]-123.046875*q0c[1]+34.09975027401226*q0r[0]-34.09975027401226*q0l[0]); 
  out0[2] += J*D*((-3.1640625*q0r[2])-3.1640625*q0l[2]-141.328125*q0c[2]-12.2543613688594*q0r[1]+12.2543613688594*q0l[1]+12.57788237343632*q0r[0]+12.57788237343632*q0l[0]-25.15576474687264*q0c[0]); 

  const double *q1l = &ql[3]; 
  const double *q1c = &qc[3]; 
  const double *q1r = &qr[3]; 
  double *out1 = &out[3]; 

  out1[0] += J*D*(35.21807064562169*q1r[2]+35.21807064562169*q1l[2]-70.43614129124337*q1c[2]-34.09975027401226*q1r[1]+34.09975027401226*q1l[1]+19.6875*q1r[0]+19.6875*q1l[0]-39.375*q1c[0]); 
  out1[1] += J*D*(51.46831774920947*q1r[2]-51.46831774920947*q1l[2]-56.6015625*q1r[1]-56.6015625*q1l[1]-123.046875*q1c[1]+34.09975027401226*q1r[0]-34.09975027401226*q1l[0]); 
  out1[2] += J*D*((-3.1640625*q1r[2])-3.1640625*q1l[2]-141.328125*q1c[2]-12.2543613688594*q1r[1]+12.2543613688594*q1l[1]+12.57788237343632*q1r[0]+12.57788237343632*q1l[0]-25.15576474687264*q1c[0]); 

  const double *q2l = &ql[6]; 
  const double *q2c = &qc[6]; 
  const double *q2r = &qr[6]; 
  double *out2 = &out[6]; 

  out2[0] += J*D*(35.21807064562169*q2r[2]+35.21807064562169*q2l[2]-70.43614129124337*q2c[2]-34.09975027401226*q2r[1]+34.09975027401226*q2l[1]+19.6875*q2r[0]+19.6875*q2l[0]-39.375*q2c[0]); 
  out2[1] += J*D*(51.46831774920947*q2r[2]-51.46831774920947*q2l[2]-56.6015625*q2r[1]-56.6015625*q2l[1]-123.046875*q2c[1]+34.09975027401226*q2r[0]-34.09975027401226*q2l[0]); 
  out2[2] += J*D*((-3.1640625*q2r[2])-3.1640625*q2l[2]-141.328125*q2c[2]-12.2543613688594*q2r[1]+12.2543613688594*q2l[1]+12.57788237343632*q2r[0]+12.57788237343632*q2l[0]-25.15576474687264*q2c[0]); 

  return 0.;

} 
