#include <gkyl_dg_diffusion_kernels.h> 
GKYL_CU_DH double dg_diffusion6_pkpm_surfx_3x_ser_p2(const double* w, const double* dx, double D, 
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

  out0[0] += J*D*(35.21807064562169*q0r[7]+35.21807064562169*q0l[7]-70.43614129124337*q0c[7]-34.09975027401226*q0r[1]+34.09975027401226*q0l[1]+19.6875*q0r[0]+19.6875*q0l[0]-39.375*q0c[0]); 
  out0[1] += J*D*(51.46831774920947*q0r[7]-51.46831774920947*q0l[7]-56.6015625*q0r[1]-56.6015625*q0l[1]-123.046875*q0c[1]+34.09975027401226*q0r[0]-34.09975027401226*q0l[0]); 
  out0[2] += J*D*(35.21807064562168*q0r[11]+35.21807064562168*q0l[11]-70.43614129124336*q0c[11]-34.09975027401226*q0r[4]+34.09975027401226*q0l[4]+19.6875*q0r[2]+19.6875*q0l[2]-39.375*q0c[2]); 
  out0[3] += J*D*(35.21807064562168*q0r[13]+35.21807064562168*q0l[13]-70.43614129124336*q0c[13]-34.09975027401226*q0r[5]+34.09975027401226*q0l[5]+19.6875*q0r[3]+19.6875*q0l[3]-39.375*q0c[3]); 
  out0[4] += J*D*(51.4683177492095*q0r[11]-51.4683177492095*q0l[11]-56.6015625*q0r[4]-56.6015625*q0l[4]-123.046875*q0c[4]+34.09975027401226*q0r[2]-34.09975027401226*q0l[2]); 
  out0[5] += J*D*(51.4683177492095*q0r[13]-51.4683177492095*q0l[13]-56.6015625*q0r[5]-56.6015625*q0l[5]-123.046875*q0c[5]+34.09975027401226*q0r[3]-34.09975027401226*q0l[3]); 
  out0[6] += J*D*(35.21807064562169*q0r[17]+35.21807064562169*q0l[17]-70.43614129124337*q0c[17]-34.09975027401226*q0r[10]+34.09975027401226*q0l[10]+19.6875*q0r[6]+19.6875*q0l[6]-39.375*q0c[6]); 
  out0[7] += J*D*((-3.1640625*q0r[7])-3.1640625*q0l[7]-141.328125*q0c[7]-12.2543613688594*q0r[1]+12.2543613688594*q0l[1]+12.57788237343632*q0r[0]+12.57788237343632*q0l[0]-25.15576474687264*q0c[0]); 
  out0[8] += J*D*((-34.09975027401227*q0r[12])+34.09975027401227*q0l[12]+19.6875*q0r[8]+19.6875*q0l[8]-39.375*q0c[8]); 
  out0[9] += J*D*((-34.09975027401227*q0r[15])+34.09975027401227*q0l[15]+19.6875*q0r[9]+19.6875*q0l[9]-39.375*q0c[9]); 
  out0[10] += J*D*(51.46831774920947*q0r[17]-51.46831774920947*q0l[17]-56.6015625*q0r[10]-56.6015625*q0l[10]-123.046875*q0c[10]+34.09975027401226*q0r[6]-34.09975027401226*q0l[6]); 
  out0[11] += J*D*((-3.1640625*q0r[11])-3.1640625*q0l[11]-141.328125*q0c[11]-12.25436136885941*q0r[4]+12.25436136885941*q0l[4]+12.57788237343632*q0r[2]+12.57788237343632*q0l[2]-25.15576474687263*q0c[2]); 
  out0[12] += J*D*((-56.6015625*q0r[12])-56.6015625*q0l[12]-123.046875*q0c[12]+34.09975027401227*q0r[8]-34.09975027401227*q0l[8]); 
  out0[13] += J*D*((-3.1640625*q0r[13])-3.1640625*q0l[13]-141.328125*q0c[13]-12.25436136885941*q0r[5]+12.25436136885941*q0l[5]+12.57788237343632*q0r[3]+12.57788237343632*q0l[3]-25.15576474687263*q0c[3]); 
  out0[14] += J*D*((-34.09975027401227*q0r[18])+34.09975027401227*q0l[18]+19.6875*q0r[14]+19.6875*q0l[14]-39.375*q0c[14]); 
  out0[15] += J*D*((-56.6015625*q0r[15])-56.6015625*q0l[15]-123.046875*q0c[15]+34.09975027401227*q0r[9]-34.09975027401227*q0l[9]); 
  out0[16] += J*D*((-34.09975027401227*q0r[19])+34.09975027401227*q0l[19]+19.6875*q0r[16]+19.6875*q0l[16]-39.375*q0c[16]); 
  out0[17] += J*D*((-3.1640625*q0r[17])-3.1640625*q0l[17]-141.328125*q0c[17]-12.2543613688594*q0r[10]+12.2543613688594*q0l[10]+12.57788237343632*q0r[6]+12.57788237343632*q0l[6]-25.15576474687264*q0c[6]); 
  out0[18] += J*D*((-56.6015625*q0r[18])-56.6015625*q0l[18]-123.046875*q0c[18]+34.09975027401227*q0r[14]-34.09975027401227*q0l[14]); 
  out0[19] += J*D*((-56.6015625*q0r[19])-56.6015625*q0l[19]-123.046875*q0c[19]+34.09975027401227*q0r[16]-34.09975027401227*q0l[16]); 

  const double *q1l = &ql[20]; 
  const double *q1c = &qc[20]; 
  const double *q1r = &qr[20]; 
  double *out1 = &out[20]; 

  out1[0] += J*D*(35.21807064562169*q1r[7]+35.21807064562169*q1l[7]-70.43614129124337*q1c[7]-34.09975027401226*q1r[1]+34.09975027401226*q1l[1]+19.6875*q1r[0]+19.6875*q1l[0]-39.375*q1c[0]); 
  out1[1] += J*D*(51.46831774920947*q1r[7]-51.46831774920947*q1l[7]-56.6015625*q1r[1]-56.6015625*q1l[1]-123.046875*q1c[1]+34.09975027401226*q1r[0]-34.09975027401226*q1l[0]); 
  out1[2] += J*D*(35.21807064562168*q1r[11]+35.21807064562168*q1l[11]-70.43614129124336*q1c[11]-34.09975027401226*q1r[4]+34.09975027401226*q1l[4]+19.6875*q1r[2]+19.6875*q1l[2]-39.375*q1c[2]); 
  out1[3] += J*D*(35.21807064562168*q1r[13]+35.21807064562168*q1l[13]-70.43614129124336*q1c[13]-34.09975027401226*q1r[5]+34.09975027401226*q1l[5]+19.6875*q1r[3]+19.6875*q1l[3]-39.375*q1c[3]); 
  out1[4] += J*D*(51.4683177492095*q1r[11]-51.4683177492095*q1l[11]-56.6015625*q1r[4]-56.6015625*q1l[4]-123.046875*q1c[4]+34.09975027401226*q1r[2]-34.09975027401226*q1l[2]); 
  out1[5] += J*D*(51.4683177492095*q1r[13]-51.4683177492095*q1l[13]-56.6015625*q1r[5]-56.6015625*q1l[5]-123.046875*q1c[5]+34.09975027401226*q1r[3]-34.09975027401226*q1l[3]); 
  out1[6] += J*D*(35.21807064562169*q1r[17]+35.21807064562169*q1l[17]-70.43614129124337*q1c[17]-34.09975027401226*q1r[10]+34.09975027401226*q1l[10]+19.6875*q1r[6]+19.6875*q1l[6]-39.375*q1c[6]); 
  out1[7] += J*D*((-3.1640625*q1r[7])-3.1640625*q1l[7]-141.328125*q1c[7]-12.2543613688594*q1r[1]+12.2543613688594*q1l[1]+12.57788237343632*q1r[0]+12.57788237343632*q1l[0]-25.15576474687264*q1c[0]); 
  out1[8] += J*D*((-34.09975027401227*q1r[12])+34.09975027401227*q1l[12]+19.6875*q1r[8]+19.6875*q1l[8]-39.375*q1c[8]); 
  out1[9] += J*D*((-34.09975027401227*q1r[15])+34.09975027401227*q1l[15]+19.6875*q1r[9]+19.6875*q1l[9]-39.375*q1c[9]); 
  out1[10] += J*D*(51.46831774920947*q1r[17]-51.46831774920947*q1l[17]-56.6015625*q1r[10]-56.6015625*q1l[10]-123.046875*q1c[10]+34.09975027401226*q1r[6]-34.09975027401226*q1l[6]); 
  out1[11] += J*D*((-3.1640625*q1r[11])-3.1640625*q1l[11]-141.328125*q1c[11]-12.25436136885941*q1r[4]+12.25436136885941*q1l[4]+12.57788237343632*q1r[2]+12.57788237343632*q1l[2]-25.15576474687263*q1c[2]); 
  out1[12] += J*D*((-56.6015625*q1r[12])-56.6015625*q1l[12]-123.046875*q1c[12]+34.09975027401227*q1r[8]-34.09975027401227*q1l[8]); 
  out1[13] += J*D*((-3.1640625*q1r[13])-3.1640625*q1l[13]-141.328125*q1c[13]-12.25436136885941*q1r[5]+12.25436136885941*q1l[5]+12.57788237343632*q1r[3]+12.57788237343632*q1l[3]-25.15576474687263*q1c[3]); 
  out1[14] += J*D*((-34.09975027401227*q1r[18])+34.09975027401227*q1l[18]+19.6875*q1r[14]+19.6875*q1l[14]-39.375*q1c[14]); 
  out1[15] += J*D*((-56.6015625*q1r[15])-56.6015625*q1l[15]-123.046875*q1c[15]+34.09975027401227*q1r[9]-34.09975027401227*q1l[9]); 
  out1[16] += J*D*((-34.09975027401227*q1r[19])+34.09975027401227*q1l[19]+19.6875*q1r[16]+19.6875*q1l[16]-39.375*q1c[16]); 
  out1[17] += J*D*((-3.1640625*q1r[17])-3.1640625*q1l[17]-141.328125*q1c[17]-12.2543613688594*q1r[10]+12.2543613688594*q1l[10]+12.57788237343632*q1r[6]+12.57788237343632*q1l[6]-25.15576474687264*q1c[6]); 
  out1[18] += J*D*((-56.6015625*q1r[18])-56.6015625*q1l[18]-123.046875*q1c[18]+34.09975027401227*q1r[14]-34.09975027401227*q1l[14]); 
  out1[19] += J*D*((-56.6015625*q1r[19])-56.6015625*q1l[19]-123.046875*q1c[19]+34.09975027401227*q1r[16]-34.09975027401227*q1l[16]); 

  const double *q2l = &ql[40]; 
  const double *q2c = &qc[40]; 
  const double *q2r = &qr[40]; 
  double *out2 = &out[40]; 

  out2[0] += J*D*(35.21807064562169*q2r[7]+35.21807064562169*q2l[7]-70.43614129124337*q2c[7]-34.09975027401226*q2r[1]+34.09975027401226*q2l[1]+19.6875*q2r[0]+19.6875*q2l[0]-39.375*q2c[0]); 
  out2[1] += J*D*(51.46831774920947*q2r[7]-51.46831774920947*q2l[7]-56.6015625*q2r[1]-56.6015625*q2l[1]-123.046875*q2c[1]+34.09975027401226*q2r[0]-34.09975027401226*q2l[0]); 
  out2[2] += J*D*(35.21807064562168*q2r[11]+35.21807064562168*q2l[11]-70.43614129124336*q2c[11]-34.09975027401226*q2r[4]+34.09975027401226*q2l[4]+19.6875*q2r[2]+19.6875*q2l[2]-39.375*q2c[2]); 
  out2[3] += J*D*(35.21807064562168*q2r[13]+35.21807064562168*q2l[13]-70.43614129124336*q2c[13]-34.09975027401226*q2r[5]+34.09975027401226*q2l[5]+19.6875*q2r[3]+19.6875*q2l[3]-39.375*q2c[3]); 
  out2[4] += J*D*(51.4683177492095*q2r[11]-51.4683177492095*q2l[11]-56.6015625*q2r[4]-56.6015625*q2l[4]-123.046875*q2c[4]+34.09975027401226*q2r[2]-34.09975027401226*q2l[2]); 
  out2[5] += J*D*(51.4683177492095*q2r[13]-51.4683177492095*q2l[13]-56.6015625*q2r[5]-56.6015625*q2l[5]-123.046875*q2c[5]+34.09975027401226*q2r[3]-34.09975027401226*q2l[3]); 
  out2[6] += J*D*(35.21807064562169*q2r[17]+35.21807064562169*q2l[17]-70.43614129124337*q2c[17]-34.09975027401226*q2r[10]+34.09975027401226*q2l[10]+19.6875*q2r[6]+19.6875*q2l[6]-39.375*q2c[6]); 
  out2[7] += J*D*((-3.1640625*q2r[7])-3.1640625*q2l[7]-141.328125*q2c[7]-12.2543613688594*q2r[1]+12.2543613688594*q2l[1]+12.57788237343632*q2r[0]+12.57788237343632*q2l[0]-25.15576474687264*q2c[0]); 
  out2[8] += J*D*((-34.09975027401227*q2r[12])+34.09975027401227*q2l[12]+19.6875*q2r[8]+19.6875*q2l[8]-39.375*q2c[8]); 
  out2[9] += J*D*((-34.09975027401227*q2r[15])+34.09975027401227*q2l[15]+19.6875*q2r[9]+19.6875*q2l[9]-39.375*q2c[9]); 
  out2[10] += J*D*(51.46831774920947*q2r[17]-51.46831774920947*q2l[17]-56.6015625*q2r[10]-56.6015625*q2l[10]-123.046875*q2c[10]+34.09975027401226*q2r[6]-34.09975027401226*q2l[6]); 
  out2[11] += J*D*((-3.1640625*q2r[11])-3.1640625*q2l[11]-141.328125*q2c[11]-12.25436136885941*q2r[4]+12.25436136885941*q2l[4]+12.57788237343632*q2r[2]+12.57788237343632*q2l[2]-25.15576474687263*q2c[2]); 
  out2[12] += J*D*((-56.6015625*q2r[12])-56.6015625*q2l[12]-123.046875*q2c[12]+34.09975027401227*q2r[8]-34.09975027401227*q2l[8]); 
  out2[13] += J*D*((-3.1640625*q2r[13])-3.1640625*q2l[13]-141.328125*q2c[13]-12.25436136885941*q2r[5]+12.25436136885941*q2l[5]+12.57788237343632*q2r[3]+12.57788237343632*q2l[3]-25.15576474687263*q2c[3]); 
  out2[14] += J*D*((-34.09975027401227*q2r[18])+34.09975027401227*q2l[18]+19.6875*q2r[14]+19.6875*q2l[14]-39.375*q2c[14]); 
  out2[15] += J*D*((-56.6015625*q2r[15])-56.6015625*q2l[15]-123.046875*q2c[15]+34.09975027401227*q2r[9]-34.09975027401227*q2l[9]); 
  out2[16] += J*D*((-34.09975027401227*q2r[19])+34.09975027401227*q2l[19]+19.6875*q2r[16]+19.6875*q2l[16]-39.375*q2c[16]); 
  out2[17] += J*D*((-3.1640625*q2r[17])-3.1640625*q2l[17]-141.328125*q2c[17]-12.2543613688594*q2r[10]+12.2543613688594*q2l[10]+12.57788237343632*q2r[6]+12.57788237343632*q2l[6]-25.15576474687264*q2c[6]); 
  out2[18] += J*D*((-56.6015625*q2r[18])-56.6015625*q2l[18]-123.046875*q2c[18]+34.09975027401227*q2r[14]-34.09975027401227*q2l[14]); 
  out2[19] += J*D*((-56.6015625*q2r[19])-56.6015625*q2l[19]-123.046875*q2c[19]+34.09975027401227*q2r[16]-34.09975027401227*q2l[16]); 

  return 0.;

} 
