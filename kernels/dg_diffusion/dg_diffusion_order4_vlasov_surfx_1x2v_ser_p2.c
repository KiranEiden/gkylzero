#include <gkyl_dg_diffusion_kernels.h>

GKYL_CU_DH double dg_diffusion_order4_vlasov_surfx_1x2v_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],4.);

  out[0] += 0.0625*(107.3312629199899*coeff[0]*qr[7]+107.3312629199899*coeff[0]*ql[7]-214.6625258399798*coeff[0]*qc[7]-129.9038105676658*coeff[0]*qr[1]+129.9038105676658*coeff[0]*ql[1]+75.0*coeff[0]*qr[0]+75.0*coeff[0]*ql[0]-150.0*coeff[0]*qc[0])*Jfac; 
  out[1] += 0.03125*(290.4737509655563*coeff[0]*qr[7]-290.4737509655563*coeff[0]*ql[7]-405.0*coeff[0]*qr[1]-405.0*coeff[0]*ql[1]-990.0*coeff[0]*qc[1]+259.8076211353315*coeff[0]*qr[0]-259.8076211353315*coeff[0]*ql[0])*Jfac; 
  out[2] += 0.0625*(107.3312629199899*coeff[0]*qr[11]+107.3312629199899*coeff[0]*ql[11]-214.6625258399798*coeff[0]*qc[11]-129.9038105676658*coeff[0]*qr[4]+129.9038105676658*coeff[0]*ql[4]+75.0*coeff[0]*qr[2]+75.0*coeff[0]*ql[2]-150.0*coeff[0]*qc[2])*Jfac; 
  out[3] += 0.0625*(107.3312629199899*coeff[0]*qr[13]+107.3312629199899*coeff[0]*ql[13]-214.6625258399798*coeff[0]*qc[13]-129.9038105676658*coeff[0]*qr[5]+129.9038105676658*coeff[0]*ql[5]+75.0*coeff[0]*qr[3]+75.0*coeff[0]*ql[3]-150.0*coeff[0]*qc[3])*Jfac; 
  out[4] += 0.03125*(290.4737509655563*coeff[0]*qr[11]-290.4737509655563*coeff[0]*ql[11]-405.0*coeff[0]*qr[4]-405.0*coeff[0]*ql[4]-990.0*coeff[0]*qc[4]+259.8076211353315*coeff[0]*qr[2]-259.8076211353315*coeff[0]*ql[2])*Jfac; 
  out[5] += 0.03125*(290.4737509655563*coeff[0]*qr[13]-290.4737509655563*coeff[0]*ql[13]-405.0*coeff[0]*qr[5]-405.0*coeff[0]*ql[5]-990.0*coeff[0]*qc[5]+259.8076211353315*coeff[0]*qr[3]-259.8076211353315*coeff[0]*ql[3])*Jfac; 
  out[6] += 0.0625*(107.3312629199899*coeff[0]*qr[17]+107.3312629199899*coeff[0]*ql[17]-214.6625258399798*coeff[0]*qc[17]-129.9038105676658*coeff[0]*qr[10]+129.9038105676658*coeff[0]*ql[10]+75.0*coeff[0]*qr[6]+75.0*coeff[0]*ql[6]-150.0*coeff[0]*qc[6])*Jfac; 
  out[7] += 0.03125*(21.0*coeff[0]*qr[7]+21.0*coeff[0]*ql[7]-1302.0*coeff[0]*qc[7]-151.0463505020892*coeff[0]*qr[1]+151.0463505020892*coeff[0]*ql[1]+134.1640786499874*coeff[0]*qr[0]+134.1640786499874*coeff[0]*ql[0]-268.3281572999748*coeff[0]*qc[0])*Jfac; 
  out[8] += -0.0625*(129.9038105676658*coeff[0]*qr[12]-129.9038105676658*coeff[0]*ql[12]-75.0*coeff[0]*qr[8]-75.0*coeff[0]*ql[8]+150.0*coeff[0]*qc[8])*Jfac; 
  out[9] += -0.0625*(129.9038105676658*coeff[0]*qr[15]-129.9038105676658*coeff[0]*ql[15]-75.0*coeff[0]*qr[9]-75.0*coeff[0]*ql[9]+150.0*coeff[0]*qc[9])*Jfac; 
  out[10] += 0.03125*(290.4737509655563*coeff[0]*qr[17]-290.4737509655563*coeff[0]*ql[17]-405.0*coeff[0]*qr[10]-405.0*coeff[0]*ql[10]-990.0*coeff[0]*qc[10]+259.8076211353315*coeff[0]*qr[6]-259.8076211353315*coeff[0]*ql[6])*Jfac; 
  out[11] += 0.03125*(21.0*coeff[0]*qr[11]+21.0*coeff[0]*ql[11]-1302.0*coeff[0]*qc[11]-151.0463505020893*coeff[0]*qr[4]+151.0463505020893*coeff[0]*ql[4]+134.1640786499874*coeff[0]*qr[2]+134.1640786499874*coeff[0]*ql[2]-268.3281572999747*coeff[0]*qc[2])*Jfac; 
  out[12] += -0.03125*(405.0*coeff[0]*qr[12]+405.0*coeff[0]*ql[12]+990.0*coeff[0]*qc[12]-259.8076211353317*coeff[0]*qr[8]+259.8076211353317*coeff[0]*ql[8])*Jfac; 
  out[13] += 0.03125*(21.0*coeff[0]*qr[13]+21.0*coeff[0]*ql[13]-1302.0*coeff[0]*qc[13]-151.0463505020893*coeff[0]*qr[5]+151.0463505020893*coeff[0]*ql[5]+134.1640786499874*coeff[0]*qr[3]+134.1640786499874*coeff[0]*ql[3]-268.3281572999747*coeff[0]*qc[3])*Jfac; 
  out[14] += -0.0625*(129.9038105676658*coeff[0]*qr[18]-129.9038105676658*coeff[0]*ql[18]-75.0*coeff[0]*qr[14]-75.0*coeff[0]*ql[14]+150.0*coeff[0]*qc[14])*Jfac; 
  out[15] += -0.03125*(405.0*coeff[0]*qr[15]+405.0*coeff[0]*ql[15]+990.0*coeff[0]*qc[15]-259.8076211353317*coeff[0]*qr[9]+259.8076211353317*coeff[0]*ql[9])*Jfac; 
  out[16] += -0.0625*(129.9038105676658*coeff[0]*qr[19]-129.9038105676658*coeff[0]*ql[19]-75.0*coeff[0]*qr[16]-75.0*coeff[0]*ql[16]+150.0*coeff[0]*qc[16])*Jfac; 
  out[17] += 0.03125*(21.0*coeff[0]*qr[17]+21.0*coeff[0]*ql[17]-1302.0*coeff[0]*qc[17]-151.0463505020892*coeff[0]*qr[10]+151.0463505020892*coeff[0]*ql[10]+134.1640786499874*coeff[0]*qr[6]+134.1640786499874*coeff[0]*ql[6]-268.3281572999748*coeff[0]*qc[6])*Jfac; 
  out[18] += -0.03125*(405.0*coeff[0]*qr[18]+405.0*coeff[0]*ql[18]+990.0*coeff[0]*qc[18]-259.8076211353317*coeff[0]*qr[14]+259.8076211353317*coeff[0]*ql[14])*Jfac; 
  out[19] += -0.03125*(405.0*coeff[0]*qr[19]+405.0*coeff[0]*ql[19]+990.0*coeff[0]*qc[19]-259.8076211353317*coeff[0]*qr[16]+259.8076211353317*coeff[0]*ql[16])*Jfac; 

  return 0.;

}

GKYL_CU_DH double dg_diffusion_order4_vlasov_surfx_1x2v_ser_p2_varcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],4.);

  out[0] += 0.02209708691207959*((795.0*coeff[2]+453.1390515062677*coeff[1]+214.6625258399798*coeff[0])*qr[7]+(795.0*coeff[2]-453.1390515062677*coeff[1]+214.6625258399798*coeff[0])*ql[7]+((-330.0*coeff[2])-429.3250516799596*coeff[0])*qc[7]+((-755.2317525104463*qr[1])+755.2317525104463*ql[1]+335.4101966249685*qr[0]+335.4101966249685*ql[0]-670.8203932499371*qc[0])*coeff[2]+((-495.0*coeff[1])-259.8076211353315*coeff[0])*qr[1]+(259.8076211353315*coeff[0]-495.0*coeff[1])*ql[1]-810.0*coeff[1]*qc[1]+(259.8076211353315*qr[0]-259.8076211353315*ql[0])*coeff[1]+150.0*coeff[0]*qr[0]+150.0*coeff[0]*ql[0]-300.0*coeff[0]*qc[0])*Jfac; 
  out[1] += 0.02209708691207959*((1195.115057222525*coeff[2]+643.9875775199395*coeff[1]+290.4737509655563*coeff[0])*qr[7]+((-1195.115057222525*coeff[2])+643.9875775199395*coeff[1]-290.4737509655563*coeff[0])*ql[7]-1287.975155039879*coeff[1]*qc[7]+((-1207.476707849887*qr[1])-1207.476707849887*ql[1]-1609.968943799849*qc[1]+580.9475019311126*qr[0]-580.9475019311126*ql[0])*coeff[2]+((-779.4228634059946*coeff[1])-405.0*coeff[0])*qr[1]+(779.4228634059946*coeff[1]-405.0*coeff[0])*ql[1]-990.0*coeff[0]*qc[1]+(450.0*qr[0]+450.0*ql[0]-900.0*qc[0])*coeff[1]+259.8076211353315*coeff[0]*qr[0]-259.8076211353315*coeff[0]*ql[0])*Jfac; 
  out[2] += 0.02209708691207959*((795.0*coeff[2]+453.1390515062678*coeff[1]+214.6625258399798*coeff[0])*qr[11]+(795.0*coeff[2]-453.1390515062678*coeff[1]+214.6625258399798*coeff[0])*ql[11]+((-330.0000000000001*coeff[2])-429.3250516799596*coeff[0])*qc[11]+((-755.2317525104463*coeff[2])-495.0*coeff[1]-259.8076211353315*coeff[0])*qr[4]+(755.2317525104463*coeff[2]-495.0*coeff[1]+259.8076211353315*coeff[0])*ql[4]-810.0*coeff[1]*qc[4]+(335.4101966249685*coeff[2]+259.8076211353315*coeff[1]+150.0*coeff[0])*qr[2]+(335.4101966249685*coeff[2]-259.8076211353315*coeff[1]+150.0*coeff[0])*ql[2]+((-670.8203932499371*coeff[2])-300.0*coeff[0])*qc[2])*Jfac; 
  out[3] += 0.02209708691207959*((795.0*coeff[2]+453.1390515062678*coeff[1]+214.6625258399798*coeff[0])*qr[13]+(795.0*coeff[2]-453.1390515062678*coeff[1]+214.6625258399798*coeff[0])*ql[13]+((-330.0000000000001*coeff[2])-429.3250516799596*coeff[0])*qc[13]+((-755.2317525104463*coeff[2])-495.0*coeff[1]-259.8076211353315*coeff[0])*qr[5]+(755.2317525104463*coeff[2]-495.0*coeff[1]+259.8076211353315*coeff[0])*ql[5]-810.0*coeff[1]*qc[5]+(335.4101966249685*coeff[2]+259.8076211353315*coeff[1]+150.0*coeff[0])*qr[3]+(335.4101966249685*coeff[2]-259.8076211353315*coeff[1]+150.0*coeff[0])*ql[3]+((-670.8203932499371*coeff[2])-300.0*coeff[0])*qc[3])*Jfac; 
  out[4] += 0.02209708691207959*((1195.115057222525*coeff[2]+643.9875775199394*coeff[1]+290.4737509655563*coeff[0])*qr[11]+((-1195.115057222525*coeff[2])+643.9875775199394*coeff[1]-290.4737509655563*coeff[0])*ql[11]-1287.975155039879*coeff[1]*qc[11]+((-1207.476707849887*coeff[2])-779.4228634059946*coeff[1]-405.0*coeff[0])*qr[4]+((-1207.476707849887*coeff[2])+779.4228634059946*coeff[1]-405.0*coeff[0])*ql[4]+((-1609.968943799849*coeff[2])-990.0*coeff[0])*qc[4]+(580.9475019311126*coeff[2]+450.0*coeff[1]+259.8076211353315*coeff[0])*qr[2]+((-580.9475019311126*coeff[2])+450.0*coeff[1]-259.8076211353315*coeff[0])*ql[2]-900.0*coeff[1]*qc[2])*Jfac; 
  out[5] += 0.02209708691207959*((1195.115057222525*coeff[2]+643.9875775199394*coeff[1]+290.4737509655563*coeff[0])*qr[13]+((-1195.115057222525*coeff[2])+643.9875775199394*coeff[1]-290.4737509655563*coeff[0])*ql[13]-1287.975155039879*coeff[1]*qc[13]+((-1207.476707849887*coeff[2])-779.4228634059946*coeff[1]-405.0*coeff[0])*qr[5]+((-1207.476707849887*coeff[2])+779.4228634059946*coeff[1]-405.0*coeff[0])*ql[5]+((-1609.968943799849*coeff[2])-990.0*coeff[0])*qc[5]+(580.9475019311126*coeff[2]+450.0*coeff[1]+259.8076211353315*coeff[0])*qr[3]+((-580.9475019311126*coeff[2])+450.0*coeff[1]-259.8076211353315*coeff[0])*ql[3]-900.0*coeff[1]*qc[3])*Jfac; 
  out[6] += 0.02209708691207959*((795.0*coeff[2]+453.1390515062677*coeff[1]+214.6625258399798*coeff[0])*qr[17]+(795.0*coeff[2]-453.1390515062677*coeff[1]+214.6625258399798*coeff[0])*ql[17]+((-330.0*coeff[2])-429.3250516799596*coeff[0])*qc[17]+((-755.2317525104463*coeff[2])-495.0*coeff[1]-259.8076211353315*coeff[0])*qr[10]+(755.2317525104463*coeff[2]-495.0*coeff[1]+259.8076211353315*coeff[0])*ql[10]-810.0*coeff[1]*qc[10]+(335.4101966249685*coeff[2]+259.8076211353315*coeff[1]+150.0*coeff[0])*qr[6]+(335.4101966249685*coeff[2]-259.8076211353315*coeff[1]+150.0*coeff[0])*ql[6]+((-670.8203932499371*coeff[2])-300.0*coeff[0])*qc[6])*Jfac; 
  out[7] += 0.0110485434560398*((2206.999093792293*coeff[2]+618.3421383020891*coeff[1]+42.0*coeff[0])*qr[7]+(2206.999093792293*coeff[2]-618.3421383020891*coeff[1]+42.0*coeff[0])*ql[7]+((-1596.55253593485*coeff[2])-2604.0*coeff[0])*qc[7]+((-2468.17240078565*qr[1])+2468.17240078565*ql[1]+1320.0*qr[0]+1320.0*ql[0]-2640.0*qc[0])*coeff[2]+((-986.1059780774073*coeff[1])-302.0927010041785*coeff[0])*qr[1]+(302.0927010041785*coeff[0]-986.1059780774073*coeff[1])*ql[1]-2535.701086484762*coeff[1]*qc[1]+(650.6612021628459*qr[0]-650.6612021628459*ql[0])*coeff[1]+268.3281572999748*coeff[0]*qr[0]+268.3281572999748*coeff[0]*ql[0]-536.6563145999496*coeff[0]*qc[0])*Jfac; 
  out[8] += -0.02209708691207959*((755.2317525104464*coeff[2]+495.0*coeff[1]+259.8076211353317*coeff[0])*qr[12]+((-755.2317525104464*coeff[2])+495.0*coeff[1]-259.8076211353317*coeff[0])*ql[12]+809.9999999999998*coeff[1]*qc[12]+((-335.4101966249685*coeff[2])-259.8076211353315*coeff[1]-150.0*coeff[0])*qr[8]+((-335.4101966249685*coeff[2])+259.8076211353315*coeff[1]-150.0*coeff[0])*ql[8]+(670.8203932499371*coeff[2]+300.0*coeff[0])*qc[8])*Jfac; 
  out[9] += -0.02209708691207959*((755.2317525104464*coeff[2]+495.0*coeff[1]+259.8076211353317*coeff[0])*qr[15]+((-755.2317525104464*coeff[2])+495.0*coeff[1]-259.8076211353317*coeff[0])*ql[15]+809.9999999999998*coeff[1]*qc[15]+((-335.4101966249685*coeff[2])-259.8076211353315*coeff[1]-150.0*coeff[0])*qr[9]+((-335.4101966249685*coeff[2])+259.8076211353315*coeff[1]-150.0*coeff[0])*ql[9]+(670.8203932499371*coeff[2]+300.0*coeff[0])*qc[9])*Jfac; 
  out[10] += 0.02209708691207959*((1195.115057222525*coeff[2]+643.9875775199395*coeff[1]+290.4737509655563*coeff[0])*qr[17]+((-1195.115057222525*coeff[2])+643.9875775199395*coeff[1]-290.4737509655563*coeff[0])*ql[17]-1287.975155039879*coeff[1]*qc[17]+((-1207.476707849887*coeff[2])-779.4228634059946*coeff[1]-405.0*coeff[0])*qr[10]+((-1207.476707849887*coeff[2])+779.4228634059946*coeff[1]-405.0*coeff[0])*ql[10]+((-1609.968943799849*coeff[2])-990.0*coeff[0])*qc[10]+(580.9475019311126*coeff[2]+450.0*coeff[1]+259.8076211353315*coeff[0])*qr[6]+((-580.9475019311126*coeff[2])+450.0*coeff[1]-259.8076211353315*coeff[0])*ql[6]-900.0*coeff[1]*qc[6])*Jfac; 
  out[11] += 0.0110485434560398*((2206.999093792293*coeff[2]+618.3421383020891*coeff[1]+42.0*coeff[0])*qr[11]+(2206.999093792293*coeff[2]-618.3421383020891*coeff[1]+42.0*coeff[0])*ql[11]+((-1596.55253593485*coeff[2])-2604.0*coeff[0])*qc[11]+((-2468.172400785651*coeff[2])-986.1059780774071*coeff[1]-302.0927010041785*coeff[0])*qr[4]+(2468.172400785651*coeff[2]-986.1059780774071*coeff[1]+302.0927010041785*coeff[0])*ql[4]-2535.70108648476*coeff[1]*qc[4]+(1320.0*coeff[2]+650.661202162846*coeff[1]+268.3281572999747*coeff[0])*qr[2]+(1320.0*coeff[2]-650.661202162846*coeff[1]+268.3281572999747*coeff[0])*ql[2]+((-2640.0*coeff[2])-536.6563145999495*coeff[0])*qc[2])*Jfac; 
  out[12] += -0.02209708691207959*((1207.476707849887*coeff[2]+779.4228634059946*coeff[1]+405.0*coeff[0])*qr[12]+(1207.476707849887*coeff[2]-779.4228634059946*coeff[1]+405.0*coeff[0])*ql[12]+(1609.968943799849*coeff[2]+990.0*coeff[0])*qc[12]+((-580.9475019311126*coeff[2])-450.0000000000001*coeff[1]-259.8076211353317*coeff[0])*qr[8]+(580.9475019311126*coeff[2]-450.0000000000001*coeff[1]+259.8076211353317*coeff[0])*ql[8]+900.0000000000001*coeff[1]*qc[8])*Jfac; 
  out[13] += 0.0110485434560398*((2206.999093792293*coeff[2]+618.3421383020891*coeff[1]+42.0*coeff[0])*qr[13]+(2206.999093792293*coeff[2]-618.3421383020891*coeff[1]+42.0*coeff[0])*ql[13]+((-1596.55253593485*coeff[2])-2604.0*coeff[0])*qc[13]+((-2468.172400785651*coeff[2])-986.1059780774071*coeff[1]-302.0927010041785*coeff[0])*qr[5]+(2468.172400785651*coeff[2]-986.1059780774071*coeff[1]+302.0927010041785*coeff[0])*ql[5]-2535.70108648476*coeff[1]*qc[5]+(1320.0*coeff[2]+650.661202162846*coeff[1]+268.3281572999747*coeff[0])*qr[3]+(1320.0*coeff[2]-650.661202162846*coeff[1]+268.3281572999747*coeff[0])*ql[3]+((-2640.0*coeff[2])-536.6563145999495*coeff[0])*qc[3])*Jfac; 
  out[14] += -0.02209708691207959*((755.2317525104464*coeff[2]+495.0*coeff[1]+259.8076211353317*coeff[0])*qr[18]+((-755.2317525104464*coeff[2])+495.0*coeff[1]-259.8076211353317*coeff[0])*ql[18]+809.9999999999998*coeff[1]*qc[18]+((-335.4101966249685*coeff[2])-259.8076211353315*coeff[1]-150.0*coeff[0])*qr[14]+((-335.4101966249685*coeff[2])+259.8076211353315*coeff[1]-150.0*coeff[0])*ql[14]+(670.8203932499371*coeff[2]+300.0*coeff[0])*qc[14])*Jfac; 
  out[15] += -0.02209708691207959*((1207.476707849887*coeff[2]+779.4228634059946*coeff[1]+405.0*coeff[0])*qr[15]+(1207.476707849887*coeff[2]-779.4228634059946*coeff[1]+405.0*coeff[0])*ql[15]+(1609.968943799849*coeff[2]+990.0*coeff[0])*qc[15]+((-580.9475019311126*coeff[2])-450.0000000000001*coeff[1]-259.8076211353317*coeff[0])*qr[9]+(580.9475019311126*coeff[2]-450.0000000000001*coeff[1]+259.8076211353317*coeff[0])*ql[9]+900.0000000000001*coeff[1]*qc[9])*Jfac; 
  out[16] += -0.02209708691207959*((755.2317525104464*coeff[2]+495.0*coeff[1]+259.8076211353317*coeff[0])*qr[19]+((-755.2317525104464*coeff[2])+495.0*coeff[1]-259.8076211353317*coeff[0])*ql[19]+809.9999999999998*coeff[1]*qc[19]+((-335.4101966249685*coeff[2])-259.8076211353315*coeff[1]-150.0*coeff[0])*qr[16]+((-335.4101966249685*coeff[2])+259.8076211353315*coeff[1]-150.0*coeff[0])*ql[16]+(670.8203932499371*coeff[2]+300.0*coeff[0])*qc[16])*Jfac; 
  out[17] += 0.0110485434560398*((2206.999093792293*coeff[2]+618.3421383020891*coeff[1]+42.0*coeff[0])*qr[17]+(2206.999093792293*coeff[2]-618.3421383020891*coeff[1]+42.0*coeff[0])*ql[17]+((-1596.55253593485*coeff[2])-2604.0*coeff[0])*qc[17]+((-2468.17240078565*coeff[2])-986.1059780774073*coeff[1]-302.0927010041785*coeff[0])*qr[10]+(2468.17240078565*coeff[2]-986.1059780774073*coeff[1]+302.0927010041785*coeff[0])*ql[10]-2535.701086484762*coeff[1]*qc[10]+(1320.0*coeff[2]+650.6612021628459*coeff[1]+268.3281572999748*coeff[0])*qr[6]+(1320.0*coeff[2]-650.6612021628459*coeff[1]+268.3281572999748*coeff[0])*ql[6]+((-2640.0*coeff[2])-536.6563145999496*coeff[0])*qc[6])*Jfac; 
  out[18] += -0.02209708691207959*((1207.476707849887*coeff[2]+779.4228634059946*coeff[1]+405.0*coeff[0])*qr[18]+(1207.476707849887*coeff[2]-779.4228634059946*coeff[1]+405.0*coeff[0])*ql[18]+(1609.968943799849*coeff[2]+990.0*coeff[0])*qc[18]+((-580.9475019311126*coeff[2])-450.0000000000001*coeff[1]-259.8076211353317*coeff[0])*qr[14]+(580.9475019311126*coeff[2]-450.0000000000001*coeff[1]+259.8076211353317*coeff[0])*ql[14]+900.0000000000001*coeff[1]*qc[14])*Jfac; 
  out[19] += -0.02209708691207959*((1207.476707849887*coeff[2]+779.4228634059946*coeff[1]+405.0*coeff[0])*qr[19]+(1207.476707849887*coeff[2]-779.4228634059946*coeff[1]+405.0*coeff[0])*ql[19]+(1609.968943799849*coeff[2]+990.0*coeff[0])*qc[19]+((-580.9475019311126*coeff[2])-450.0000000000001*coeff[1]-259.8076211353317*coeff[0])*qr[16]+(580.9475019311126*coeff[2]-450.0000000000001*coeff[1]+259.8076211353317*coeff[0])*ql[16]+900.0000000000001*coeff[1]*qc[16])*Jfac; 

  return 0.;

}

