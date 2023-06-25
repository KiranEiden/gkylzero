#include <gkyl_dg_diffusion_kernels.h> 
GKYL_CU_DH void dg_diffusion6_surfz_3x_tensor_p2(const double* w, const double* dx, double D, 
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
  const double J = pow(dx1, 6.0);

  const double *q0l = &ql[0]; 
  const double *q0c = &qc[0]; 
  const double *q0r = &qr[0]; 
  double *out0= &out[0]; 

  out0[0] += J*D*(35.21807064562169*q0r[9]+35.21807064562169*q0l[9]-70.43614129124337*q0c[9]-34.09975027401226*q0r[3]+34.09975027401226*q0l[3]+19.6875*q0r[0]+19.6875*q0l[0]-39.375*q0c[0]); 
  out0[1] += J*D*(35.21807064562168*q0r[15]+35.21807064562168*q0l[15]-70.43614129124336*q0c[15]-34.09975027401226*q0r[5]+34.09975027401226*q0l[5]+19.6875*q0r[1]+19.6875*q0l[1]-39.375*q0c[1]); 
  out0[2] += J*D*(35.21807064562168*q0r[16]+35.21807064562168*q0l[16]-70.43614129124336*q0c[16]-34.09975027401226*q0r[6]+34.09975027401226*q0l[6]+19.6875*q0r[2]+19.6875*q0l[2]-39.375*q0c[2]); 
  out0[3] += J*D*(51.46831774920947*q0r[9]-51.46831774920947*q0l[9]-56.6015625*q0r[3]-56.6015625*q0l[3]-123.046875*q0c[3]+34.09975027401226*q0r[0]-34.09975027401226*q0l[0]); 
  out0[4] += J*D*(35.21807064562169*q0r[19]+35.21807064562169*q0l[19]-70.43614129124337*q0c[19]-34.09975027401226*q0r[10]+34.09975027401226*q0l[10]+19.6875*q0r[4]+19.6875*q0l[4]-39.375*q0c[4]); 
  out0[5] += J*D*(51.4683177492095*q0r[15]-51.4683177492095*q0l[15]-56.6015625*q0r[5]-56.6015625*q0l[5]-123.046875*q0c[5]+34.09975027401226*q0r[1]-34.09975027401226*q0l[1]); 
  out0[6] += J*D*(51.4683177492095*q0r[16]-51.4683177492095*q0l[16]-56.6015625*q0r[6]-56.6015625*q0l[6]-123.046875*q0c[6]+34.09975027401226*q0r[2]-34.09975027401226*q0l[2]); 
  out0[7] += J*D*(35.21807064562169*q0r[21]+35.21807064562169*q0l[21]-70.43614129124337*q0c[21]-34.09975027401227*q0r[13]+34.09975027401227*q0l[13]+19.6875*q0r[7]+19.6875*q0l[7]-39.375*q0c[7]); 
  out0[8] += J*D*(35.21807064562169*q0r[22]+35.21807064562169*q0l[22]-70.43614129124337*q0c[22]-34.09975027401227*q0r[14]+34.09975027401227*q0l[14]+19.6875*q0r[8]+19.6875*q0l[8]-39.375*q0c[8]); 
  out0[9] += J*D*((-3.1640625*q0r[9])-3.1640625*q0l[9]-141.328125*q0c[9]-12.2543613688594*q0r[3]+12.2543613688594*q0l[3]+12.57788237343632*q0r[0]+12.57788237343632*q0l[0]-25.15576474687264*q0c[0]); 
  out0[10] += J*D*(51.46831774920947*q0r[19]-51.46831774920947*q0l[19]-56.6015625*q0r[10]-56.6015625*q0l[10]-123.046875*q0c[10]+34.09975027401226*q0r[4]-34.09975027401226*q0l[4]); 
  out0[11] += J*D*(35.21807064562168*q0r[24]+35.21807064562168*q0l[24]-70.43614129124336*q0c[24]-34.09975027401227*q0r[17]+34.09975027401227*q0l[17]+19.6875*q0r[11]+19.6875*q0l[11]-39.375*q0c[11]); 
  out0[12] += J*D*(35.21807064562168*q0r[25]+35.21807064562168*q0l[25]-70.43614129124336*q0c[25]-34.09975027401227*q0r[18]+34.09975027401227*q0l[18]+19.6875*q0r[12]+19.6875*q0l[12]-39.375*q0c[12]); 
  out0[13] += J*D*(51.4683177492095*q0r[21]-51.4683177492095*q0l[21]-56.6015625*q0r[13]-56.6015625*q0l[13]-123.046875*q0c[13]+34.09975027401227*q0r[7]-34.09975027401227*q0l[7]); 
  out0[14] += J*D*(51.4683177492095*q0r[22]-51.4683177492095*q0l[22]-56.6015625*q0r[14]-56.6015625*q0l[14]-123.046875*q0c[14]+34.09975027401227*q0r[8]-34.09975027401227*q0l[8]); 
  out0[15] += J*D*((-3.1640625*q0r[15])-3.1640625*q0l[15]-141.328125*q0c[15]-12.25436136885941*q0r[5]+12.25436136885941*q0l[5]+12.57788237343632*q0r[1]+12.57788237343632*q0l[1]-25.15576474687263*q0c[1]); 
  out0[16] += J*D*((-3.1640625*q0r[16])-3.1640625*q0l[16]-141.328125*q0c[16]-12.25436136885941*q0r[6]+12.25436136885941*q0l[6]+12.57788237343632*q0r[2]+12.57788237343632*q0l[2]-25.15576474687263*q0c[2]); 
  out0[17] += J*D*(51.46831774920947*q0r[24]-51.46831774920947*q0l[24]-56.6015625*q0r[17]-56.6015625*q0l[17]-123.046875*q0c[17]+34.09975027401227*q0r[11]-34.09975027401227*q0l[11]); 
  out0[18] += J*D*(51.46831774920947*q0r[25]-51.46831774920947*q0l[25]-56.6015625*q0r[18]-56.6015625*q0l[18]-123.046875*q0c[18]+34.09975027401227*q0r[12]-34.09975027401227*q0l[12]); 
  out0[19] += J*D*((-3.1640625*q0r[19])-3.1640625*q0l[19]-141.328125*q0c[19]-12.2543613688594*q0r[10]+12.2543613688594*q0l[10]+12.57788237343632*q0r[4]+12.57788237343632*q0l[4]-25.15576474687264*q0c[4]); 
  out0[20] += J*D*(35.21807064562169*q0r[26]+35.21807064562169*q0l[26]-70.43614129124337*q0c[26]-34.09975027401226*q0r[23]+34.09975027401226*q0l[23]+19.6875*q0r[20]+19.6875*q0l[20]-39.375*q0c[20]); 
  out0[21] += J*D*((-3.1640625*q0r[21])-3.1640625*q0l[21]-141.328125*q0c[21]-12.25436136885941*q0r[13]+12.25436136885941*q0l[13]+12.57788237343632*q0r[7]+12.57788237343632*q0l[7]-25.15576474687264*q0c[7]); 
  out0[22] += J*D*((-3.1640625*q0r[22])-3.1640625*q0l[22]-141.328125*q0c[22]-12.25436136885941*q0r[14]+12.25436136885941*q0l[14]+12.57788237343632*q0r[8]+12.57788237343632*q0l[8]-25.15576474687264*q0c[8]); 
  out0[23] += J*D*(51.46831774920947*q0r[26]-51.46831774920947*q0l[26]-56.6015625*q0r[23]-56.6015625*q0l[23]-123.046875*q0c[23]+34.09975027401226*q0r[20]-34.09975027401226*q0l[20]); 
  out0[24] += J*D*((-3.1640625*q0r[24])-3.1640625*q0l[24]-141.328125*q0c[24]-12.2543613688594*q0r[17]+12.2543613688594*q0l[17]+12.57788237343632*q0r[11]+12.57788237343632*q0l[11]-25.15576474687263*q0c[11]); 
  out0[25] += J*D*((-3.1640625*q0r[25])-3.1640625*q0l[25]-141.328125*q0c[25]-12.2543613688594*q0r[18]+12.2543613688594*q0l[18]+12.57788237343632*q0r[12]+12.57788237343632*q0l[12]-25.15576474687263*q0c[12]); 
  out0[26] += J*D*((-3.1640625*q0r[26])-3.1640625*q0l[26]-141.328125*q0c[26]-12.2543613688594*q0r[23]+12.2543613688594*q0l[23]+12.57788237343632*q0r[20]+12.57788237343632*q0l[20]-25.15576474687264*q0c[20]); 

} 
