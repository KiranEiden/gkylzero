#include <gkyl_dg_diffusion_kernels.h> 
GKYL_CU_DH void dg_diffusion6_pkpm_surfx_2x_ser_p2(const double* w, const double* dx, double D, 
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
  double *out0= &out[0]; 

  out0[0] += J*D*(35.21807064562169*q0r[4]+35.21807064562169*q0l[4]-70.43614129124337*q0c[4]-34.09975027401226*q0r[1]+34.09975027401226*q0l[1]+19.6875*q0r[0]+19.6875*q0l[0]-39.375*q0c[0]); 
  out0[1] += J*D*(51.46831774920947*q0r[4]-51.46831774920947*q0l[4]-56.6015625*q0r[1]-56.6015625*q0l[1]-123.046875*q0c[1]+34.09975027401226*q0r[0]-34.09975027401226*q0l[0]); 
  out0[2] += J*D*(35.21807064562168*q0r[6]+35.21807064562168*q0l[6]-70.43614129124336*q0c[6]-34.09975027401226*q0r[3]+34.09975027401226*q0l[3]+19.6875*q0r[2]+19.6875*q0l[2]-39.375*q0c[2]); 
  out0[3] += J*D*(51.4683177492095*q0r[6]-51.4683177492095*q0l[6]-56.6015625*q0r[3]-56.6015625*q0l[3]-123.046875*q0c[3]+34.09975027401226*q0r[2]-34.09975027401226*q0l[2]); 
  out0[4] += J*D*((-3.1640625*q0r[4])-3.1640625*q0l[4]-141.328125*q0c[4]-12.2543613688594*q0r[1]+12.2543613688594*q0l[1]+12.57788237343632*q0r[0]+12.57788237343632*q0l[0]-25.15576474687264*q0c[0]); 
  out0[5] += J*D*((-34.09975027401227*q0r[7])+34.09975027401227*q0l[7]+19.6875*q0r[5]+19.6875*q0l[5]-39.375*q0c[5]); 
  out0[6] += J*D*((-3.1640625*q0r[6])-3.1640625*q0l[6]-141.328125*q0c[6]-12.25436136885941*q0r[3]+12.25436136885941*q0l[3]+12.57788237343632*q0r[2]+12.57788237343632*q0l[2]-25.15576474687263*q0c[2]); 
  out0[7] += J*D*((-56.6015625*q0r[7])-56.6015625*q0l[7]-123.046875*q0c[7]+34.09975027401227*q0r[5]-34.09975027401227*q0l[5]); 

  const double *q1l = &ql[8]; 
  const double *q1c = &qc[8]; 
  const double *q1r = &qr[8]; 
  double *out1= &out[8]; 

  out1[0] += J*D*(35.21807064562169*q1r[4]+35.21807064562169*q1l[4]-70.43614129124337*q1c[4]-34.09975027401226*q1r[1]+34.09975027401226*q1l[1]+19.6875*q1r[0]+19.6875*q1l[0]-39.375*q1c[0]); 
  out1[1] += J*D*(51.46831774920947*q1r[4]-51.46831774920947*q1l[4]-56.6015625*q1r[1]-56.6015625*q1l[1]-123.046875*q1c[1]+34.09975027401226*q1r[0]-34.09975027401226*q1l[0]); 
  out1[2] += J*D*(35.21807064562168*q1r[6]+35.21807064562168*q1l[6]-70.43614129124336*q1c[6]-34.09975027401226*q1r[3]+34.09975027401226*q1l[3]+19.6875*q1r[2]+19.6875*q1l[2]-39.375*q1c[2]); 
  out1[3] += J*D*(51.4683177492095*q1r[6]-51.4683177492095*q1l[6]-56.6015625*q1r[3]-56.6015625*q1l[3]-123.046875*q1c[3]+34.09975027401226*q1r[2]-34.09975027401226*q1l[2]); 
  out1[4] += J*D*((-3.1640625*q1r[4])-3.1640625*q1l[4]-141.328125*q1c[4]-12.2543613688594*q1r[1]+12.2543613688594*q1l[1]+12.57788237343632*q1r[0]+12.57788237343632*q1l[0]-25.15576474687264*q1c[0]); 
  out1[5] += J*D*((-34.09975027401227*q1r[7])+34.09975027401227*q1l[7]+19.6875*q1r[5]+19.6875*q1l[5]-39.375*q1c[5]); 
  out1[6] += J*D*((-3.1640625*q1r[6])-3.1640625*q1l[6]-141.328125*q1c[6]-12.25436136885941*q1r[3]+12.25436136885941*q1l[3]+12.57788237343632*q1r[2]+12.57788237343632*q1l[2]-25.15576474687263*q1c[2]); 
  out1[7] += J*D*((-56.6015625*q1r[7])-56.6015625*q1l[7]-123.046875*q1c[7]+34.09975027401227*q1r[5]-34.09975027401227*q1l[5]); 

  const double *q2l = &ql[16]; 
  const double *q2c = &qc[16]; 
  const double *q2r = &qr[16]; 
  double *out2= &out[16]; 

  out2[0] += J*D*(35.21807064562169*q2r[4]+35.21807064562169*q2l[4]-70.43614129124337*q2c[4]-34.09975027401226*q2r[1]+34.09975027401226*q2l[1]+19.6875*q2r[0]+19.6875*q2l[0]-39.375*q2c[0]); 
  out2[1] += J*D*(51.46831774920947*q2r[4]-51.46831774920947*q2l[4]-56.6015625*q2r[1]-56.6015625*q2l[1]-123.046875*q2c[1]+34.09975027401226*q2r[0]-34.09975027401226*q2l[0]); 
  out2[2] += J*D*(35.21807064562168*q2r[6]+35.21807064562168*q2l[6]-70.43614129124336*q2c[6]-34.09975027401226*q2r[3]+34.09975027401226*q2l[3]+19.6875*q2r[2]+19.6875*q2l[2]-39.375*q2c[2]); 
  out2[3] += J*D*(51.4683177492095*q2r[6]-51.4683177492095*q2l[6]-56.6015625*q2r[3]-56.6015625*q2l[3]-123.046875*q2c[3]+34.09975027401226*q2r[2]-34.09975027401226*q2l[2]); 
  out2[4] += J*D*((-3.1640625*q2r[4])-3.1640625*q2l[4]-141.328125*q2c[4]-12.2543613688594*q2r[1]+12.2543613688594*q2l[1]+12.57788237343632*q2r[0]+12.57788237343632*q2l[0]-25.15576474687264*q2c[0]); 
  out2[5] += J*D*((-34.09975027401227*q2r[7])+34.09975027401227*q2l[7]+19.6875*q2r[5]+19.6875*q2l[5]-39.375*q2c[5]); 
  out2[6] += J*D*((-3.1640625*q2r[6])-3.1640625*q2l[6]-141.328125*q2c[6]-12.25436136885941*q2r[3]+12.25436136885941*q2l[3]+12.57788237343632*q2r[2]+12.57788237343632*q2l[2]-25.15576474687263*q2c[2]); 
  out2[7] += J*D*((-56.6015625*q2r[7])-56.6015625*q2l[7]-123.046875*q2c[7]+34.09975027401227*q2r[5]-34.09975027401227*q2l[5]); 

} 
