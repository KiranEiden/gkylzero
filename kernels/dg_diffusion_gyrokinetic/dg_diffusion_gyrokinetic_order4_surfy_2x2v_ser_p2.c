#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_surfy_2x2v_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[1],4.);

  out[0] += 0.0625*(107.3312629199899*coeff[1]*qr[12]+107.3312629199899*coeff[1]*ql[12]-214.6625258399798*coeff[1]*qc[12]-129.9038105676658*coeff[1]*qr[2]+129.9038105676658*coeff[1]*ql[2]+(75.0*qr[0]+75.0*ql[0]-150.0*qc[0])*coeff[1])*Jfac; 
  out[1] += 0.0625*(107.3312629199899*coeff[1]*qr[20]+107.3312629199899*coeff[1]*ql[20]-214.6625258399798*coeff[1]*qc[20]-129.9038105676658*coeff[1]*qr[5]+129.9038105676658*coeff[1]*ql[5]+75.0*coeff[1]*qr[1]+75.0*coeff[1]*ql[1]-150.0*coeff[1]*qc[1])*Jfac; 
  out[2] += 0.03125*(290.4737509655563*coeff[1]*qr[12]-290.4737509655563*coeff[1]*ql[12]-405.0*coeff[1]*qr[2]-405.0*coeff[1]*ql[2]-990.0*coeff[1]*qc[2]+(259.8076211353315*qr[0]-259.8076211353315*ql[0])*coeff[1])*Jfac; 
  out[3] += 0.0625*(107.3312629199899*coeff[1]*qr[22]+107.3312629199899*coeff[1]*ql[22]-214.6625258399798*coeff[1]*qc[22]-129.9038105676658*coeff[1]*qr[7]+129.9038105676658*coeff[1]*ql[7]+75.0*coeff[1]*qr[3]+75.0*coeff[1]*ql[3]-150.0*coeff[1]*qc[3])*Jfac; 
  out[4] += 0.0625*(107.3312629199899*coeff[1]*qr[26]+107.3312629199899*coeff[1]*ql[26]-214.6625258399798*coeff[1]*qc[26]-129.9038105676658*coeff[1]*qr[9]+129.9038105676658*coeff[1]*ql[9]+75.0*coeff[1]*qr[4]+75.0*coeff[1]*ql[4]-150.0*coeff[1]*qc[4])*Jfac; 
  out[5] += 0.03125*(290.4737509655563*coeff[1]*qr[20]-290.4737509655563*coeff[1]*ql[20]-405.0*coeff[1]*qr[5]-405.0*coeff[1]*ql[5]-990.0*coeff[1]*qc[5]+259.8076211353315*coeff[1]*qr[1]-259.8076211353315*coeff[1]*ql[1])*Jfac; 
  out[6] += 0.0625*(107.3312629199899*coeff[1]*qr[33]+107.3312629199899*coeff[1]*ql[33]-214.6625258399798*coeff[1]*qc[33]-129.9038105676658*coeff[1]*qr[15]+129.9038105676658*coeff[1]*ql[15]+75.0*coeff[1]*qr[6]+75.0*coeff[1]*ql[6]-150.0*coeff[1]*qc[6])*Jfac; 
  out[7] += 0.03125*(290.4737509655563*coeff[1]*qr[22]-290.4737509655563*coeff[1]*ql[22]-405.0*coeff[1]*qr[7]-405.0*coeff[1]*ql[7]-990.0*coeff[1]*qc[7]+259.8076211353315*coeff[1]*qr[3]-259.8076211353315*coeff[1]*ql[3])*Jfac; 
  out[8] += 0.0625*(107.3312629199899*coeff[1]*qr[36]+107.3312629199899*coeff[1]*ql[36]-214.6625258399798*coeff[1]*qc[36]-129.9038105676658*coeff[1]*qr[16]+129.9038105676658*coeff[1]*ql[16]+75.0*coeff[1]*qr[8]+75.0*coeff[1]*ql[8]-150.0*coeff[1]*qc[8])*Jfac; 
  out[9] += 0.03125*(290.4737509655563*coeff[1]*qr[26]-290.4737509655563*coeff[1]*ql[26]-405.0*coeff[1]*qr[9]-405.0*coeff[1]*ql[9]-990.0*coeff[1]*qc[9]+259.8076211353315*coeff[1]*qr[4]-259.8076211353315*coeff[1]*ql[4])*Jfac; 
  out[10] += 0.0625*(107.3312629199899*coeff[1]*qr[38]+107.3312629199899*coeff[1]*ql[38]-214.6625258399798*coeff[1]*qc[38]-129.9038105676658*coeff[1]*qr[18]+129.9038105676658*coeff[1]*ql[18]+75.0*coeff[1]*qr[10]+75.0*coeff[1]*ql[10]-150.0*coeff[1]*qc[10])*Jfac; 
  out[11] += -0.0625*(129.9038105676658*coeff[1]*qr[19]-129.9038105676658*coeff[1]*ql[19]-75.0*coeff[1]*qr[11]-75.0*coeff[1]*ql[11]+150.0*coeff[1]*qc[11])*Jfac; 
  out[12] += 0.03125*(21.0*coeff[1]*qr[12]+21.0*coeff[1]*ql[12]-1302.0*coeff[1]*qc[12]-151.0463505020892*coeff[1]*qr[2]+151.0463505020892*coeff[1]*ql[2]+(134.1640786499874*qr[0]+134.1640786499874*ql[0]-268.3281572999748*qc[0])*coeff[1])*Jfac; 
  out[13] += -0.0625*(129.9038105676658*coeff[1]*qr[24]-129.9038105676658*coeff[1]*ql[24]-75.0*coeff[1]*qr[13]-75.0*coeff[1]*ql[13]+150.0*coeff[1]*qc[13])*Jfac; 
  out[14] += -0.0625*(129.9038105676658*coeff[1]*qr[29]-129.9038105676658*coeff[1]*ql[29]-75.0*coeff[1]*qr[14]-75.0*coeff[1]*ql[14]+150.0*coeff[1]*qc[14])*Jfac; 
  out[15] += 0.03125*(290.4737509655563*coeff[1]*qr[33]-290.4737509655563*coeff[1]*ql[33]-405.0*coeff[1]*qr[15]-405.0*coeff[1]*ql[15]-990.0*coeff[1]*qc[15]+259.8076211353315*coeff[1]*qr[6]-259.8076211353315*coeff[1]*ql[6])*Jfac; 
  out[16] += 0.03125*(290.4737509655563*coeff[1]*qr[36]-290.4737509655563*coeff[1]*ql[36]-405.0*coeff[1]*qr[16]-405.0*coeff[1]*ql[16]-990.0*coeff[1]*qc[16]+259.8076211353315*coeff[1]*qr[8]-259.8076211353315*coeff[1]*ql[8])*Jfac; 
  out[17] += 0.0625*(107.3312629199899*coeff[1]*qr[45]+107.3312629199899*coeff[1]*ql[45]-214.6625258399798*coeff[1]*qc[45]-129.9038105676658*coeff[1]*qr[31]+129.9038105676658*coeff[1]*ql[31]+75.0*coeff[1]*qr[17]+75.0*coeff[1]*ql[17]-150.0*coeff[1]*qc[17])*Jfac; 
  out[18] += 0.03125*(290.4737509655563*coeff[1]*qr[38]-290.4737509655563*coeff[1]*ql[38]-405.0*coeff[1]*qr[18]-405.0*coeff[1]*ql[18]-990.0*coeff[1]*qc[18]+259.8076211353315*coeff[1]*qr[10]-259.8076211353315*coeff[1]*ql[10])*Jfac; 
  out[19] += -0.03125*(405.0*coeff[1]*qr[19]+405.0*coeff[1]*ql[19]+990.0*coeff[1]*qc[19]-259.8076211353317*coeff[1]*qr[11]+259.8076211353317*coeff[1]*ql[11])*Jfac; 
  out[20] += 0.03125*(21.0*coeff[1]*qr[20]+21.0*coeff[1]*ql[20]-1302.0*coeff[1]*qc[20]-151.0463505020893*coeff[1]*qr[5]+151.0463505020893*coeff[1]*ql[5]+134.1640786499874*coeff[1]*qr[1]+134.1640786499874*coeff[1]*ql[1]-268.3281572999747*coeff[1]*qc[1])*Jfac; 
  out[21] += -0.0625*(129.9038105676658*coeff[1]*qr[32]-129.9038105676658*coeff[1]*ql[32]-75.0*coeff[1]*qr[21]-75.0*coeff[1]*ql[21]+150.0*coeff[1]*qc[21])*Jfac; 
  out[22] += 0.03125*(21.0*coeff[1]*qr[22]+21.0*coeff[1]*ql[22]-1302.0*coeff[1]*qc[22]-151.0463505020893*coeff[1]*qr[7]+151.0463505020893*coeff[1]*ql[7]+134.1640786499874*coeff[1]*qr[3]+134.1640786499874*coeff[1]*ql[3]-268.3281572999747*coeff[1]*qc[3])*Jfac; 
  out[23] += -0.0625*(129.9038105676658*coeff[1]*qr[34]-129.9038105676658*coeff[1]*ql[34]-75.0*coeff[1]*qr[23]-75.0*coeff[1]*ql[23]+150.0*coeff[1]*qc[23])*Jfac; 
  out[24] += -0.03125*(405.0*coeff[1]*qr[24]+405.0*coeff[1]*ql[24]+990.0*coeff[1]*qc[24]-259.8076211353317*coeff[1]*qr[13]+259.8076211353317*coeff[1]*ql[13])*Jfac; 
  out[25] += -0.0625*(129.9038105676658*coeff[1]*qr[35]-129.9038105676658*coeff[1]*ql[35]-75.0*coeff[1]*qr[25]-75.0*coeff[1]*ql[25]+150.0*coeff[1]*qc[25])*Jfac; 
  out[26] += 0.03125*(21.0*coeff[1]*qr[26]+21.0*coeff[1]*ql[26]-1302.0*coeff[1]*qc[26]-151.0463505020893*coeff[1]*qr[9]+151.0463505020893*coeff[1]*ql[9]+134.1640786499874*coeff[1]*qr[4]+134.1640786499874*coeff[1]*ql[4]-268.3281572999747*coeff[1]*qc[4])*Jfac; 
  out[27] += -0.0625*(129.9038105676658*coeff[1]*qr[40]-129.9038105676658*coeff[1]*ql[40]-75.0*coeff[1]*qr[27]-75.0*coeff[1]*ql[27]+150.0*coeff[1]*qc[27])*Jfac; 
  out[28] += -0.0625*(129.9038105676658*coeff[1]*qr[41]-129.9038105676658*coeff[1]*ql[41]-75.0*coeff[1]*qr[28]-75.0*coeff[1]*ql[28]+150.0*coeff[1]*qc[28])*Jfac; 
  out[29] += -0.03125*(405.0*coeff[1]*qr[29]+405.0*coeff[1]*ql[29]+990.0*coeff[1]*qc[29]-259.8076211353317*coeff[1]*qr[14]+259.8076211353317*coeff[1]*ql[14])*Jfac; 
  out[30] += -0.0625*(129.9038105676658*coeff[1]*qr[43]-129.9038105676658*coeff[1]*ql[43]-75.0*coeff[1]*qr[30]-75.0*coeff[1]*ql[30]+150.0*coeff[1]*qc[30])*Jfac; 
  out[31] += 0.03125*(290.4737509655563*coeff[1]*qr[45]-290.4737509655563*coeff[1]*ql[45]-405.0*coeff[1]*qr[31]-405.0*coeff[1]*ql[31]-990.0*coeff[1]*qc[31]+259.8076211353315*coeff[1]*qr[17]-259.8076211353315*coeff[1]*ql[17])*Jfac; 
  out[32] += -0.03125*(405.0*coeff[1]*qr[32]+405.0*coeff[1]*ql[32]+990.0*coeff[1]*qc[32]-259.8076211353317*coeff[1]*qr[21]+259.8076211353317*coeff[1]*ql[21])*Jfac; 
  out[33] += 0.03125*(21.0*coeff[1]*qr[33]+21.0*coeff[1]*ql[33]-1302.0*coeff[1]*qc[33]-151.0463505020892*coeff[1]*qr[15]+151.0463505020892*coeff[1]*ql[15]+134.1640786499874*coeff[1]*qr[6]+134.1640786499874*coeff[1]*ql[6]-268.3281572999748*coeff[1]*qc[6])*Jfac; 
  out[34] += -0.03125*(405.0*coeff[1]*qr[34]+405.0*coeff[1]*ql[34]+990.0*coeff[1]*qc[34]-259.8076211353317*coeff[1]*qr[23]+259.8076211353317*coeff[1]*ql[23])*Jfac; 
  out[35] += -0.03125*(405.0*coeff[1]*qr[35]+405.0*coeff[1]*ql[35]+990.0*coeff[1]*qc[35]-259.8076211353317*coeff[1]*qr[25]+259.8076211353317*coeff[1]*ql[25])*Jfac; 
  out[36] += 0.03125*(21.0*coeff[1]*qr[36]+21.0*coeff[1]*ql[36]-1302.0*coeff[1]*qc[36]-151.0463505020892*coeff[1]*qr[16]+151.0463505020892*coeff[1]*ql[16]+134.1640786499874*coeff[1]*qr[8]+134.1640786499874*coeff[1]*ql[8]-268.3281572999748*coeff[1]*qc[8])*Jfac; 
  out[37] += -0.0625*(129.9038105676658*coeff[1]*qr[44]-129.9038105676658*coeff[1]*ql[44]-75.0*coeff[1]*qr[37]-75.0*coeff[1]*ql[37]+150.0*coeff[1]*qc[37])*Jfac; 
  out[38] += 0.03125*(21.0*coeff[1]*qr[38]+21.0*coeff[1]*ql[38]-1302.0*coeff[1]*qc[38]-151.0463505020892*coeff[1]*qr[18]+151.0463505020892*coeff[1]*ql[18]+134.1640786499874*coeff[1]*qr[10]+134.1640786499874*coeff[1]*ql[10]-268.3281572999748*coeff[1]*qc[10])*Jfac; 
  out[39] += -0.0625*(129.9038105676658*coeff[1]*qr[46]-129.9038105676658*coeff[1]*ql[46]-75.0*coeff[1]*qr[39]-75.0*coeff[1]*ql[39]+150.0*coeff[1]*qc[39])*Jfac; 
  out[40] += -0.03125*(405.0*coeff[1]*qr[40]+405.0*coeff[1]*ql[40]+990.0*coeff[1]*qc[40]-259.8076211353317*coeff[1]*qr[27]+259.8076211353317*coeff[1]*ql[27])*Jfac; 
  out[41] += -0.03125*(405.0*coeff[1]*qr[41]+405.0*coeff[1]*ql[41]+990.0*coeff[1]*qc[41]-259.8076211353317*coeff[1]*qr[28]+259.8076211353317*coeff[1]*ql[28])*Jfac; 
  out[42] += -0.0625*(129.9038105676658*coeff[1]*qr[47]-129.9038105676658*coeff[1]*ql[47]-75.0*coeff[1]*qr[42]-75.0*coeff[1]*ql[42]+150.0*coeff[1]*qc[42])*Jfac; 
  out[43] += -0.03125*(405.0*coeff[1]*qr[43]+405.0*coeff[1]*ql[43]+990.0*coeff[1]*qc[43]-259.8076211353317*coeff[1]*qr[30]+259.8076211353317*coeff[1]*ql[30])*Jfac; 
  out[44] += -0.03125*(405.0*coeff[1]*qr[44]+405.0*coeff[1]*ql[44]+990.0*coeff[1]*qc[44]-259.8076211353317*coeff[1]*qr[37]+259.8076211353317*coeff[1]*ql[37])*Jfac; 
  out[45] += 0.03125*(21.0*coeff[1]*qr[45]+21.0*coeff[1]*ql[45]-1302.0*coeff[1]*qc[45]-151.0463505020893*coeff[1]*qr[31]+151.0463505020893*coeff[1]*ql[31]+134.1640786499874*coeff[1]*qr[17]+134.1640786499874*coeff[1]*ql[17]-268.3281572999747*coeff[1]*qc[17])*Jfac; 
  out[46] += -0.03125*(405.0*coeff[1]*qr[46]+405.0*coeff[1]*ql[46]+990.0*coeff[1]*qc[46]-259.8076211353317*coeff[1]*qr[39]+259.8076211353317*coeff[1]*ql[39])*Jfac; 
  out[47] += -0.03125*(405.0*coeff[1]*qr[47]+405.0*coeff[1]*ql[47]+990.0*coeff[1]*qc[47]-259.8076211353317*coeff[1]*qr[42]+259.8076211353317*coeff[1]*ql[42])*Jfac; 

  return 0.;

}

