#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order4_surfz_3x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[2],4.);

  out[0] += 0.0625*(107.3312629199899*coeff[2]*qr[9]+107.3312629199899*coeff[2]*ql[9]-214.6625258399798*coeff[2]*qc[9]-129.9038105676658*coeff[2]*qr[3]+129.9038105676658*coeff[2]*ql[3]+(75.0*qr[0]+75.0*ql[0]-150.0*qc[0])*coeff[2])*Jfac; 
  out[1] += 0.0625*(107.3312629199899*coeff[2]*qr[15]+107.3312629199899*coeff[2]*ql[15]-214.6625258399798*coeff[2]*qc[15]-129.9038105676658*coeff[2]*qr[5]+129.9038105676658*coeff[2]*ql[5]+(75.0*qr[1]+75.0*ql[1]-150.0*qc[1])*coeff[2])*Jfac; 
  out[2] += 0.0625*(107.3312629199899*coeff[2]*qr[16]+107.3312629199899*coeff[2]*ql[16]-214.6625258399798*coeff[2]*qc[16]-129.9038105676658*coeff[2]*qr[6]+129.9038105676658*coeff[2]*ql[6]+75.0*coeff[2]*qr[2]+75.0*coeff[2]*ql[2]-150.0*coeff[2]*qc[2])*Jfac; 
  out[3] += 0.03125*(290.4737509655563*coeff[2]*qr[9]-290.4737509655563*coeff[2]*ql[9]-405.0*coeff[2]*qr[3]-405.0*coeff[2]*ql[3]-990.0*coeff[2]*qc[3]+(259.8076211353315*qr[0]-259.8076211353315*ql[0])*coeff[2])*Jfac; 
  out[4] += 0.0625*(107.3312629199899*coeff[2]*qr[19]+107.3312629199899*coeff[2]*ql[19]-214.6625258399798*coeff[2]*qc[19]-129.9038105676658*coeff[2]*qr[10]+129.9038105676658*coeff[2]*ql[10]+75.0*coeff[2]*qr[4]+75.0*coeff[2]*ql[4]-150.0*coeff[2]*qc[4])*Jfac; 
  out[5] += 0.03125*(290.4737509655563*coeff[2]*qr[15]-290.4737509655563*coeff[2]*ql[15]-405.0*coeff[2]*qr[5]-405.0*coeff[2]*ql[5]-990.0*coeff[2]*qc[5]+(259.8076211353315*qr[1]-259.8076211353315*ql[1])*coeff[2])*Jfac; 
  out[6] += 0.03125*(290.4737509655563*coeff[2]*qr[16]-290.4737509655563*coeff[2]*ql[16]-405.0*coeff[2]*qr[6]-405.0*coeff[2]*ql[6]-990.0*coeff[2]*qc[6]+259.8076211353315*coeff[2]*qr[2]-259.8076211353315*coeff[2]*ql[2])*Jfac; 
  out[7] += -0.0625*(129.9038105676658*coeff[2]*qr[13]-129.9038105676658*coeff[2]*ql[13]-75.0*coeff[2]*qr[7]-75.0*coeff[2]*ql[7]+150.0*coeff[2]*qc[7])*Jfac; 
  out[8] += -0.0625*(129.9038105676658*coeff[2]*qr[14]-129.9038105676658*coeff[2]*ql[14]-75.0*coeff[2]*qr[8]-75.0*coeff[2]*ql[8]+150.0*coeff[2]*qc[8])*Jfac; 
  out[9] += 0.03125*(21.0*coeff[2]*qr[9]+21.0*coeff[2]*ql[9]-1302.0*coeff[2]*qc[9]-151.0463505020892*coeff[2]*qr[3]+151.0463505020892*coeff[2]*ql[3]+(134.1640786499874*qr[0]+134.1640786499874*ql[0]-268.3281572999748*qc[0])*coeff[2])*Jfac; 
  out[10] += 0.03125*(290.4737509655563*coeff[2]*qr[19]-290.4737509655563*coeff[2]*ql[19]-405.0*coeff[2]*qr[10]-405.0*coeff[2]*ql[10]-990.0*coeff[2]*qc[10]+259.8076211353315*coeff[2]*qr[4]-259.8076211353315*coeff[2]*ql[4])*Jfac; 
  out[11] += -0.0625*(129.9038105676658*coeff[2]*qr[17]-129.9038105676658*coeff[2]*ql[17]-75.0*coeff[2]*qr[11]-75.0*coeff[2]*ql[11]+150.0*coeff[2]*qc[11])*Jfac; 
  out[12] += -0.0625*(129.9038105676658*coeff[2]*qr[18]-129.9038105676658*coeff[2]*ql[18]-75.0*coeff[2]*qr[12]-75.0*coeff[2]*ql[12]+150.0*coeff[2]*qc[12])*Jfac; 
  out[13] += -0.03125*(405.0*coeff[2]*qr[13]+405.0*coeff[2]*ql[13]+990.0*coeff[2]*qc[13]-259.8076211353317*coeff[2]*qr[7]+259.8076211353317*coeff[2]*ql[7])*Jfac; 
  out[14] += -0.03125*(405.0*coeff[2]*qr[14]+405.0*coeff[2]*ql[14]+990.0*coeff[2]*qc[14]-259.8076211353317*coeff[2]*qr[8]+259.8076211353317*coeff[2]*ql[8])*Jfac; 
  out[15] += 0.03125*(21.0*coeff[2]*qr[15]+21.0*coeff[2]*ql[15]-1302.0*coeff[2]*qc[15]-151.0463505020893*coeff[2]*qr[5]+151.0463505020893*coeff[2]*ql[5]+(134.1640786499874*qr[1]+134.1640786499874*ql[1]-268.3281572999747*qc[1])*coeff[2])*Jfac; 
  out[16] += 0.03125*(21.0*coeff[2]*qr[16]+21.0*coeff[2]*ql[16]-1302.0*coeff[2]*qc[16]-151.0463505020893*coeff[2]*qr[6]+151.0463505020893*coeff[2]*ql[6]+134.1640786499874*coeff[2]*qr[2]+134.1640786499874*coeff[2]*ql[2]-268.3281572999747*coeff[2]*qc[2])*Jfac; 
  out[17] += -0.03125*(405.0*coeff[2]*qr[17]+405.0*coeff[2]*ql[17]+990.0*coeff[2]*qc[17]-259.8076211353317*coeff[2]*qr[11]+259.8076211353317*coeff[2]*ql[11])*Jfac; 
  out[18] += -0.03125*(405.0*coeff[2]*qr[18]+405.0*coeff[2]*ql[18]+990.0*coeff[2]*qc[18]-259.8076211353317*coeff[2]*qr[12]+259.8076211353317*coeff[2]*ql[12])*Jfac; 
  out[19] += 0.03125*(21.0*coeff[2]*qr[19]+21.0*coeff[2]*ql[19]-1302.0*coeff[2]*qc[19]-151.0463505020892*coeff[2]*qr[10]+151.0463505020892*coeff[2]*ql[10]+134.1640786499874*coeff[2]*qr[4]+134.1640786499874*coeff[2]*ql[4]-268.3281572999748*coeff[2]*qc[4])*Jfac; 

  return 0.;

}

