#include <gkyl_vlasov_sr_kernels.h> 
GKYL_CU_DH double vlasov_sr_surfy_2x2v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // gamma:     Particle Lorentz boost factor sqrt(1 + p^2).
  // fl/fc/fr:  Input Distribution function in left/center/right cells.
  // out:       Incremented distribution function in center cell.
  const double dx11 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[3]; 
  const double dv = dxv[3], wv = w[3]; 

  double p_over_gamma[8] = {0.0}; 
  p_over_gamma[0] = 1.732050807568877*gamma[2]*dv11; 
  p_over_gamma[1] = 1.732050807568877*gamma[3]*dv11; 
  p_over_gamma[2] = 3.872983346207417*gamma[5]*dv11; 
  p_over_gamma[3] = 3.872983346207417*gamma[7]*dv11; 
  p_over_gamma[4] = 1.732050807568877*gamma[6]*dv11; 
  double alpha[16] = {0.0}; 
  alpha[0] = 1.414213562373095*p_over_gamma[0]; 
  alpha[2] = 1.414213562373095*p_over_gamma[1]; 
  alpha[3] = 1.414213562373095*p_over_gamma[2]; 
  alpha[6] = 1.414213562373095*p_over_gamma[3]; 
  alpha[8] = 1.414213562373095*p_over_gamma[4]; 

  double Ghat_r[16]; 
  double Ghat_l[16]; 
  if (wv>0) { 

  Ghat_r[0] = alpha[8]*(0.4330127018922194*fc[18]+0.25*fc[16])+alpha[6]*(0.4330127018922193*fc[14]+0.25*fc[10])+0.4330127018922193*(alpha[3]*fc[9]+alpha[2]*fc[7])+0.25*(alpha[3]*fc[4]+alpha[2]*fc[3])+alpha[0]*(0.4330127018922193*fc[2]+0.25*fc[0]); 
  Ghat_r[1] = alpha[8]*(0.4330127018922193*fc[20]+0.2500000000000001*fc[17])+alpha[6]*(0.4330127018922193*fc[15]+0.25*fc[13])+0.4330127018922193*(alpha[3]*fc[12]+alpha[2]*fc[11])+0.25*(alpha[3]*fc[8]+alpha[2]*fc[6])+alpha[0]*(0.4330127018922193*fc[5]+0.25*fc[1]); 
  Ghat_r[2] = alpha[6]*(0.3872983346207416*fc[22]+0.223606797749979*fc[19])+alpha[2]*(0.3872983346207417*fc[18]+0.223606797749979*fc[16])+alpha[3]*(0.4330127018922193*fc[14]+0.25*fc[10])+0.4330127018922193*alpha[6]*fc[9]+(0.3872983346207416*fc[7]+0.223606797749979*fc[3])*alpha[8]+0.4330127018922193*alpha[0]*fc[7]+0.25*(fc[4]*alpha[6]+alpha[0]*fc[3])+alpha[2]*(0.4330127018922193*fc[2]+0.25*fc[0]); 
  Ghat_r[3] = alpha[6]*(0.3872983346207416*fc[30]+0.223606797749979*fc[27])+alpha[3]*(0.3872983346207417*fc[26]+0.223606797749979*fc[24])+alpha[8]*(0.4330127018922193*fc[22]+0.2500000000000001*fc[19])+alpha[2]*(0.4330127018922193*fc[14]+0.25*fc[10])+0.4330127018922193*(alpha[0]*fc[9]+alpha[6]*fc[7])+0.25*(fc[3]*alpha[6]+alpha[0]*fc[4])+(0.4330127018922193*fc[2]+0.25*fc[0])*alpha[3]; 
  Ghat_r[4] = alpha[6]*(0.3872983346207417*fc[23]+0.223606797749979*fc[21])+alpha[2]*(0.3872983346207416*fc[20]+0.223606797749979*fc[17])+alpha[3]*(0.4330127018922193*fc[15]+0.25*fc[13])+0.4330127018922193*alpha[6]*fc[12]+(0.3872983346207416*alpha[8]+0.4330127018922193*alpha[0])*fc[11]+0.25*alpha[6]*fc[8]+fc[6]*(0.223606797749979*alpha[8]+0.25*alpha[0])+alpha[2]*(0.4330127018922193*fc[5]+0.25*fc[1]); 
  Ghat_r[5] = alpha[6]*(0.3872983346207417*fc[31]+0.223606797749979*fc[29])+alpha[3]*(0.3872983346207416*fc[28]+0.223606797749979*fc[25])+alpha[8]*(0.4330127018922194*fc[23]+0.25*fc[21])+alpha[2]*(0.4330127018922193*fc[15]+0.25*fc[13])+0.4330127018922193*(alpha[0]*fc[12]+alpha[6]*fc[11])+0.25*(alpha[0]*fc[8]+alpha[6]*fc[6])+alpha[3]*(0.4330127018922193*fc[5]+0.25*fc[1]); 
  Ghat_r[6] = alpha[3]*(0.3872983346207416*fc[30]+0.223606797749979*fc[27])+alpha[6]*(0.3872983346207417*fc[26]+0.223606797749979*fc[24])+alpha[2]*(0.3872983346207416*fc[22]+0.223606797749979*fc[19])+alpha[6]*(0.3872983346207417*fc[18]+0.223606797749979*fc[16])+(0.3872983346207416*alpha[8]+0.4330127018922193*alpha[0])*fc[14]+(0.223606797749979*alpha[8]+0.25*alpha[0])*fc[10]+0.4330127018922193*(alpha[2]*fc[9]+alpha[3]*fc[7]+fc[2]*alpha[6])+0.25*(fc[0]*alpha[6]+alpha[2]*fc[4]+alpha[3]*fc[3]); 
  Ghat_r[7] = alpha[3]*(0.3872983346207417*fc[31]+0.223606797749979*fc[29])+alpha[6]*(0.3872983346207416*fc[28]+0.223606797749979*fc[25])+alpha[2]*(0.3872983346207417*fc[23]+0.223606797749979*fc[21])+alpha[6]*(0.3872983346207416*fc[20]+0.223606797749979*fc[17])+(0.3872983346207416*alpha[8]+0.4330127018922193*alpha[0])*fc[15]+(0.223606797749979*alpha[8]+0.25*alpha[0])*fc[13]+0.4330127018922193*(alpha[2]*fc[12]+alpha[3]*fc[11])+0.25*(alpha[2]*fc[8]+alpha[3]*fc[6])+(0.4330127018922193*fc[5]+0.25*fc[1])*alpha[6]; 
  Ghat_r[8] = alpha[3]*(0.4330127018922193*fc[22]+0.2500000000000001*fc[19])+(0.276641667586244*alpha[8]+0.4330127018922194*alpha[0])*fc[18]+(0.159719141249985*alpha[8]+0.25*alpha[0])*fc[16]+alpha[6]*(0.3872983346207416*fc[14]+0.223606797749979*fc[10])+(0.4330127018922193*fc[2]+0.25*fc[0])*alpha[8]+alpha[2]*(0.3872983346207416*fc[7]+0.223606797749979*fc[3]); 
  Ghat_r[9] = alpha[3]*(0.4330127018922193*fc[23]+0.2500000000000001*fc[21])+(0.276641667586244*alpha[8]+0.4330127018922194*alpha[0])*fc[20]+(0.159719141249985*alpha[8]+0.25*alpha[0])*fc[17]+alpha[6]*(0.3872983346207417*fc[15]+0.223606797749979*fc[13])+0.3872983346207417*alpha[2]*fc[11]+(0.4330127018922194*fc[5]+0.2500000000000001*fc[1])*alpha[8]+0.223606797749979*alpha[2]*fc[6]; 
  Ghat_r[10] = alpha[6]*(0.3464101615137754*fc[30]+0.2*fc[27])+(0.276641667586244*alpha[8]+0.4330127018922194*alpha[0])*fc[22]+(0.159719141249985*alpha[8]+0.25*alpha[0])*fc[19]+alpha[3]*(0.4330127018922193*fc[18]+0.2500000000000001*fc[16])+alpha[2]*(0.3872983346207417*fc[14]+0.223606797749979*fc[10])+alpha[8]*(0.4330127018922194*fc[9]+0.2500000000000001*fc[4])+alpha[6]*(0.3872983346207417*fc[7]+0.223606797749979*fc[3]); 
  Ghat_r[11] = alpha[6]*(0.3464101615137754*fc[31]+0.2*fc[29])+(0.276641667586244*alpha[8]+0.4330127018922194*alpha[0])*fc[23]+(0.159719141249985*alpha[8]+0.25*alpha[0])*fc[21]+alpha[3]*(0.4330127018922193*fc[20]+0.2500000000000001*fc[17])+alpha[2]*(0.3872983346207416*fc[15]+0.223606797749979*fc[13])+0.4330127018922193*alpha[8]*fc[12]+0.3872983346207416*alpha[6]*fc[11]+0.25*alpha[8]*fc[8]+0.223606797749979*alpha[6]*fc[6]; 
  Ghat_r[12] = alpha[2]*(0.4330127018922193*fc[30]+0.2500000000000001*fc[27])+alpha[0]*(0.4330127018922194*fc[26]+0.25*fc[24])+alpha[6]*(0.3872983346207416*fc[14]+0.223606797749979*fc[10])+alpha[3]*(0.3872983346207416*fc[9]+0.223606797749979*fc[4]); 
  Ghat_r[13] = alpha[2]*(0.4330127018922193*fc[31]+0.2500000000000001*fc[29])+alpha[0]*(0.4330127018922194*fc[28]+0.25*fc[25])+alpha[6]*(0.3872983346207417*fc[15]+0.223606797749979*fc[13])+alpha[3]*(0.3872983346207417*fc[12]+0.223606797749979*fc[8]); 
  Ghat_r[14] = (0.3872983346207417*alpha[8]+0.4330127018922194*alpha[0])*fc[30]+(0.223606797749979*alpha[8]+0.25*alpha[0])*fc[27]+alpha[2]*(0.4330127018922193*fc[26]+0.2500000000000001*fc[24])+alpha[6]*(0.3464101615137754*fc[22]+0.2*fc[19])+alpha[3]*(0.3872983346207417*fc[14]+0.223606797749979*fc[10])+alpha[6]*(0.3872983346207417*fc[9]+0.223606797749979*fc[4]); 
  Ghat_r[15] = (0.3872983346207417*alpha[8]+0.4330127018922194*alpha[0])*fc[31]+(0.223606797749979*alpha[8]+0.25*alpha[0])*fc[29]+alpha[2]*(0.4330127018922193*fc[28]+0.2500000000000001*fc[25])+alpha[6]*(0.3464101615137754*fc[23]+0.2*fc[21])+alpha[3]*(0.3872983346207416*fc[15]+0.223606797749979*fc[13])+alpha[6]*(0.3872983346207416*fc[12]+0.223606797749979*fc[8]); 

  Ghat_l[0] = alpha[8]*(0.4330127018922194*fl[18]+0.25*fl[16])+alpha[6]*(0.4330127018922193*fl[14]+0.25*fl[10])+0.4330127018922193*(alpha[3]*fl[9]+alpha[2]*fl[7])+0.25*(alpha[3]*fl[4]+alpha[2]*fl[3])+alpha[0]*(0.4330127018922193*fl[2]+0.25*fl[0]); 
  Ghat_l[1] = alpha[8]*(0.4330127018922193*fl[20]+0.2500000000000001*fl[17])+alpha[6]*(0.4330127018922193*fl[15]+0.25*fl[13])+0.4330127018922193*(alpha[3]*fl[12]+alpha[2]*fl[11])+0.25*(alpha[3]*fl[8]+alpha[2]*fl[6])+alpha[0]*(0.4330127018922193*fl[5]+0.25*fl[1]); 
  Ghat_l[2] = alpha[6]*(0.3872983346207416*fl[22]+0.223606797749979*fl[19])+alpha[2]*(0.3872983346207417*fl[18]+0.223606797749979*fl[16])+alpha[3]*(0.4330127018922193*fl[14]+0.25*fl[10])+0.4330127018922193*alpha[6]*fl[9]+(0.3872983346207416*fl[7]+0.223606797749979*fl[3])*alpha[8]+0.4330127018922193*alpha[0]*fl[7]+0.25*(fl[4]*alpha[6]+alpha[0]*fl[3])+alpha[2]*(0.4330127018922193*fl[2]+0.25*fl[0]); 
  Ghat_l[3] = alpha[6]*(0.3872983346207416*fl[30]+0.223606797749979*fl[27])+alpha[3]*(0.3872983346207417*fl[26]+0.223606797749979*fl[24])+alpha[8]*(0.4330127018922193*fl[22]+0.2500000000000001*fl[19])+alpha[2]*(0.4330127018922193*fl[14]+0.25*fl[10])+0.4330127018922193*(alpha[0]*fl[9]+alpha[6]*fl[7])+0.25*(fl[3]*alpha[6]+alpha[0]*fl[4])+(0.4330127018922193*fl[2]+0.25*fl[0])*alpha[3]; 
  Ghat_l[4] = alpha[6]*(0.3872983346207417*fl[23]+0.223606797749979*fl[21])+alpha[2]*(0.3872983346207416*fl[20]+0.223606797749979*fl[17])+alpha[3]*(0.4330127018922193*fl[15]+0.25*fl[13])+0.4330127018922193*alpha[6]*fl[12]+(0.3872983346207416*alpha[8]+0.4330127018922193*alpha[0])*fl[11]+0.25*alpha[6]*fl[8]+fl[6]*(0.223606797749979*alpha[8]+0.25*alpha[0])+alpha[2]*(0.4330127018922193*fl[5]+0.25*fl[1]); 
  Ghat_l[5] = alpha[6]*(0.3872983346207417*fl[31]+0.223606797749979*fl[29])+alpha[3]*(0.3872983346207416*fl[28]+0.223606797749979*fl[25])+alpha[8]*(0.4330127018922194*fl[23]+0.25*fl[21])+alpha[2]*(0.4330127018922193*fl[15]+0.25*fl[13])+0.4330127018922193*(alpha[0]*fl[12]+alpha[6]*fl[11])+0.25*(alpha[0]*fl[8]+alpha[6]*fl[6])+alpha[3]*(0.4330127018922193*fl[5]+0.25*fl[1]); 
  Ghat_l[6] = alpha[3]*(0.3872983346207416*fl[30]+0.223606797749979*fl[27])+alpha[6]*(0.3872983346207417*fl[26]+0.223606797749979*fl[24])+alpha[2]*(0.3872983346207416*fl[22]+0.223606797749979*fl[19])+alpha[6]*(0.3872983346207417*fl[18]+0.223606797749979*fl[16])+(0.3872983346207416*alpha[8]+0.4330127018922193*alpha[0])*fl[14]+(0.223606797749979*alpha[8]+0.25*alpha[0])*fl[10]+0.4330127018922193*(alpha[2]*fl[9]+alpha[3]*fl[7]+fl[2]*alpha[6])+0.25*(fl[0]*alpha[6]+alpha[2]*fl[4]+alpha[3]*fl[3]); 
  Ghat_l[7] = alpha[3]*(0.3872983346207417*fl[31]+0.223606797749979*fl[29])+alpha[6]*(0.3872983346207416*fl[28]+0.223606797749979*fl[25])+alpha[2]*(0.3872983346207417*fl[23]+0.223606797749979*fl[21])+alpha[6]*(0.3872983346207416*fl[20]+0.223606797749979*fl[17])+(0.3872983346207416*alpha[8]+0.4330127018922193*alpha[0])*fl[15]+(0.223606797749979*alpha[8]+0.25*alpha[0])*fl[13]+0.4330127018922193*(alpha[2]*fl[12]+alpha[3]*fl[11])+0.25*(alpha[2]*fl[8]+alpha[3]*fl[6])+(0.4330127018922193*fl[5]+0.25*fl[1])*alpha[6]; 
  Ghat_l[8] = alpha[3]*(0.4330127018922193*fl[22]+0.2500000000000001*fl[19])+(0.276641667586244*alpha[8]+0.4330127018922194*alpha[0])*fl[18]+(0.159719141249985*alpha[8]+0.25*alpha[0])*fl[16]+alpha[6]*(0.3872983346207416*fl[14]+0.223606797749979*fl[10])+(0.4330127018922193*fl[2]+0.25*fl[0])*alpha[8]+alpha[2]*(0.3872983346207416*fl[7]+0.223606797749979*fl[3]); 
  Ghat_l[9] = alpha[3]*(0.4330127018922193*fl[23]+0.2500000000000001*fl[21])+(0.276641667586244*alpha[8]+0.4330127018922194*alpha[0])*fl[20]+(0.159719141249985*alpha[8]+0.25*alpha[0])*fl[17]+alpha[6]*(0.3872983346207417*fl[15]+0.223606797749979*fl[13])+0.3872983346207417*alpha[2]*fl[11]+(0.4330127018922194*fl[5]+0.2500000000000001*fl[1])*alpha[8]+0.223606797749979*alpha[2]*fl[6]; 
  Ghat_l[10] = alpha[6]*(0.3464101615137754*fl[30]+0.2*fl[27])+(0.276641667586244*alpha[8]+0.4330127018922194*alpha[0])*fl[22]+(0.159719141249985*alpha[8]+0.25*alpha[0])*fl[19]+alpha[3]*(0.4330127018922193*fl[18]+0.2500000000000001*fl[16])+alpha[2]*(0.3872983346207417*fl[14]+0.223606797749979*fl[10])+alpha[8]*(0.4330127018922194*fl[9]+0.2500000000000001*fl[4])+alpha[6]*(0.3872983346207417*fl[7]+0.223606797749979*fl[3]); 
  Ghat_l[11] = alpha[6]*(0.3464101615137754*fl[31]+0.2*fl[29])+(0.276641667586244*alpha[8]+0.4330127018922194*alpha[0])*fl[23]+(0.159719141249985*alpha[8]+0.25*alpha[0])*fl[21]+alpha[3]*(0.4330127018922193*fl[20]+0.2500000000000001*fl[17])+alpha[2]*(0.3872983346207416*fl[15]+0.223606797749979*fl[13])+0.4330127018922193*alpha[8]*fl[12]+0.3872983346207416*alpha[6]*fl[11]+0.25*alpha[8]*fl[8]+0.223606797749979*alpha[6]*fl[6]; 
  Ghat_l[12] = alpha[2]*(0.4330127018922193*fl[30]+0.2500000000000001*fl[27])+alpha[0]*(0.4330127018922194*fl[26]+0.25*fl[24])+alpha[6]*(0.3872983346207416*fl[14]+0.223606797749979*fl[10])+alpha[3]*(0.3872983346207416*fl[9]+0.223606797749979*fl[4]); 
  Ghat_l[13] = alpha[2]*(0.4330127018922193*fl[31]+0.2500000000000001*fl[29])+alpha[0]*(0.4330127018922194*fl[28]+0.25*fl[25])+alpha[6]*(0.3872983346207417*fl[15]+0.223606797749979*fl[13])+alpha[3]*(0.3872983346207417*fl[12]+0.223606797749979*fl[8]); 
  Ghat_l[14] = (0.3872983346207417*alpha[8]+0.4330127018922194*alpha[0])*fl[30]+(0.223606797749979*alpha[8]+0.25*alpha[0])*fl[27]+alpha[2]*(0.4330127018922193*fl[26]+0.2500000000000001*fl[24])+alpha[6]*(0.3464101615137754*fl[22]+0.2*fl[19])+alpha[3]*(0.3872983346207417*fl[14]+0.223606797749979*fl[10])+alpha[6]*(0.3872983346207417*fl[9]+0.223606797749979*fl[4]); 
  Ghat_l[15] = (0.3872983346207417*alpha[8]+0.4330127018922194*alpha[0])*fl[31]+(0.223606797749979*alpha[8]+0.25*alpha[0])*fl[29]+alpha[2]*(0.4330127018922193*fl[28]+0.2500000000000001*fl[25])+alpha[6]*(0.3464101615137754*fl[23]+0.2*fl[21])+alpha[3]*(0.3872983346207416*fl[15]+0.223606797749979*fl[13])+alpha[6]*(0.3872983346207416*fl[12]+0.223606797749979*fl[8]); 

  } else { 

  Ghat_r[0] = (-0.4330127018922194*alpha[8]*fr[18])+0.25*alpha[8]*fr[16]-0.4330127018922193*alpha[6]*fr[14]+0.25*alpha[6]*fr[10]-0.4330127018922193*alpha[3]*fr[9]-0.4330127018922193*alpha[2]*fr[7]+0.25*alpha[3]*fr[4]+0.25*alpha[2]*fr[3]-0.4330127018922193*alpha[0]*fr[2]+0.25*alpha[0]*fr[0]; 
  Ghat_r[1] = (-0.4330127018922193*alpha[8]*fr[20])+0.2500000000000001*alpha[8]*fr[17]-0.4330127018922193*alpha[6]*fr[15]+0.25*alpha[6]*fr[13]-0.4330127018922193*alpha[3]*fr[12]-0.4330127018922193*alpha[2]*fr[11]+0.25*alpha[3]*fr[8]+0.25*alpha[2]*fr[6]-0.4330127018922193*alpha[0]*fr[5]+0.25*alpha[0]*fr[1]; 
  Ghat_r[2] = (-0.3872983346207416*alpha[6]*fr[22])+0.223606797749979*alpha[6]*fr[19]-0.3872983346207417*alpha[2]*fr[18]+0.223606797749979*alpha[2]*fr[16]-0.4330127018922193*alpha[3]*fr[14]+0.25*alpha[3]*fr[10]-0.4330127018922193*alpha[6]*fr[9]-0.3872983346207416*fr[7]*alpha[8]+0.223606797749979*fr[3]*alpha[8]-0.4330127018922193*alpha[0]*fr[7]+0.25*fr[4]*alpha[6]+0.25*alpha[0]*fr[3]-0.4330127018922193*alpha[2]*fr[2]+0.25*fr[0]*alpha[2]; 
  Ghat_r[3] = (-0.3872983346207416*alpha[6]*fr[30])+0.223606797749979*alpha[6]*fr[27]-0.3872983346207417*alpha[3]*fr[26]+0.223606797749979*alpha[3]*fr[24]-0.4330127018922193*alpha[8]*fr[22]+0.2500000000000001*alpha[8]*fr[19]-0.4330127018922193*alpha[2]*fr[14]+0.25*alpha[2]*fr[10]-0.4330127018922193*alpha[0]*fr[9]-0.4330127018922193*alpha[6]*fr[7]+0.25*fr[3]*alpha[6]+0.25*alpha[0]*fr[4]-0.4330127018922193*fr[2]*alpha[3]+0.25*fr[0]*alpha[3]; 
  Ghat_r[4] = (-0.3872983346207417*alpha[6]*fr[23])+0.223606797749979*alpha[6]*fr[21]-0.3872983346207416*alpha[2]*fr[20]+0.223606797749979*alpha[2]*fr[17]-0.4330127018922193*alpha[3]*fr[15]+0.25*alpha[3]*fr[13]-0.4330127018922193*alpha[6]*fr[12]-0.3872983346207416*alpha[8]*fr[11]-0.4330127018922193*alpha[0]*fr[11]+0.25*alpha[6]*fr[8]+0.223606797749979*fr[6]*alpha[8]+0.25*alpha[0]*fr[6]-0.4330127018922193*alpha[2]*fr[5]+0.25*fr[1]*alpha[2]; 
  Ghat_r[5] = (-0.3872983346207417*alpha[6]*fr[31])+0.223606797749979*alpha[6]*fr[29]-0.3872983346207416*alpha[3]*fr[28]+0.223606797749979*alpha[3]*fr[25]-0.4330127018922194*alpha[8]*fr[23]+0.25*alpha[8]*fr[21]-0.4330127018922193*alpha[2]*fr[15]+0.25*alpha[2]*fr[13]-0.4330127018922193*alpha[0]*fr[12]-0.4330127018922193*alpha[6]*fr[11]+0.25*alpha[0]*fr[8]+0.25*alpha[6]*fr[6]-0.4330127018922193*alpha[3]*fr[5]+0.25*fr[1]*alpha[3]; 
  Ghat_r[6] = (-0.3872983346207416*alpha[3]*fr[30])+0.223606797749979*alpha[3]*fr[27]-0.3872983346207417*alpha[6]*fr[26]+0.223606797749979*alpha[6]*fr[24]-0.3872983346207416*alpha[2]*fr[22]+0.223606797749979*alpha[2]*fr[19]-0.3872983346207417*alpha[6]*fr[18]+0.223606797749979*alpha[6]*fr[16]-0.3872983346207416*alpha[8]*fr[14]-0.4330127018922193*alpha[0]*fr[14]+0.223606797749979*alpha[8]*fr[10]+0.25*alpha[0]*fr[10]-0.4330127018922193*alpha[2]*fr[9]-0.4330127018922193*alpha[3]*fr[7]-0.4330127018922193*fr[2]*alpha[6]+0.25*fr[0]*alpha[6]+0.25*alpha[2]*fr[4]+0.25*alpha[3]*fr[3]; 
  Ghat_r[7] = (-0.3872983346207417*alpha[3]*fr[31])+0.223606797749979*alpha[3]*fr[29]-0.3872983346207416*alpha[6]*fr[28]+0.223606797749979*alpha[6]*fr[25]-0.3872983346207417*alpha[2]*fr[23]+0.223606797749979*alpha[2]*fr[21]-0.3872983346207416*alpha[6]*fr[20]+0.223606797749979*alpha[6]*fr[17]-0.3872983346207416*alpha[8]*fr[15]-0.4330127018922193*alpha[0]*fr[15]+0.223606797749979*alpha[8]*fr[13]+0.25*alpha[0]*fr[13]-0.4330127018922193*alpha[2]*fr[12]-0.4330127018922193*alpha[3]*fr[11]+0.25*alpha[2]*fr[8]+0.25*alpha[3]*fr[6]-0.4330127018922193*fr[5]*alpha[6]+0.25*fr[1]*alpha[6]; 
  Ghat_r[8] = (-0.4330127018922193*alpha[3]*fr[22])+0.2500000000000001*alpha[3]*fr[19]-0.276641667586244*alpha[8]*fr[18]-0.4330127018922194*alpha[0]*fr[18]+0.159719141249985*alpha[8]*fr[16]+0.25*alpha[0]*fr[16]-0.3872983346207416*alpha[6]*fr[14]+0.223606797749979*alpha[6]*fr[10]-0.4330127018922193*fr[2]*alpha[8]+0.25*fr[0]*alpha[8]-0.3872983346207416*alpha[2]*fr[7]+0.223606797749979*alpha[2]*fr[3]; 
  Ghat_r[9] = (-0.4330127018922193*alpha[3]*fr[23])+0.2500000000000001*alpha[3]*fr[21]-0.276641667586244*alpha[8]*fr[20]-0.4330127018922194*alpha[0]*fr[20]+0.159719141249985*alpha[8]*fr[17]+0.25*alpha[0]*fr[17]-0.3872983346207417*alpha[6]*fr[15]+0.223606797749979*alpha[6]*fr[13]-0.3872983346207417*alpha[2]*fr[11]-0.4330127018922194*fr[5]*alpha[8]+0.2500000000000001*fr[1]*alpha[8]+0.223606797749979*alpha[2]*fr[6]; 
  Ghat_r[10] = (-0.3464101615137754*alpha[6]*fr[30])+0.2*alpha[6]*fr[27]-0.276641667586244*alpha[8]*fr[22]-0.4330127018922194*alpha[0]*fr[22]+0.159719141249985*alpha[8]*fr[19]+0.25*alpha[0]*fr[19]-0.4330127018922193*alpha[3]*fr[18]+0.2500000000000001*alpha[3]*fr[16]-0.3872983346207417*alpha[2]*fr[14]+0.223606797749979*alpha[2]*fr[10]-0.4330127018922194*alpha[8]*fr[9]+0.2500000000000001*fr[4]*alpha[8]-0.3872983346207417*alpha[6]*fr[7]+0.223606797749979*fr[3]*alpha[6]; 
  Ghat_r[11] = (-0.3464101615137754*alpha[6]*fr[31])+0.2*alpha[6]*fr[29]-0.276641667586244*alpha[8]*fr[23]-0.4330127018922194*alpha[0]*fr[23]+0.159719141249985*alpha[8]*fr[21]+0.25*alpha[0]*fr[21]-0.4330127018922193*alpha[3]*fr[20]+0.2500000000000001*alpha[3]*fr[17]-0.3872983346207416*alpha[2]*fr[15]+0.223606797749979*alpha[2]*fr[13]-0.4330127018922193*alpha[8]*fr[12]-0.3872983346207416*alpha[6]*fr[11]+0.25*alpha[8]*fr[8]+0.223606797749979*alpha[6]*fr[6]; 
  Ghat_r[12] = (-0.4330127018922193*alpha[2]*fr[30])+0.2500000000000001*alpha[2]*fr[27]-0.4330127018922194*alpha[0]*fr[26]+0.25*alpha[0]*fr[24]-0.3872983346207416*alpha[6]*fr[14]+0.223606797749979*alpha[6]*fr[10]-0.3872983346207416*alpha[3]*fr[9]+0.223606797749979*alpha[3]*fr[4]; 
  Ghat_r[13] = (-0.4330127018922193*alpha[2]*fr[31])+0.2500000000000001*alpha[2]*fr[29]-0.4330127018922194*alpha[0]*fr[28]+0.25*alpha[0]*fr[25]-0.3872983346207417*alpha[6]*fr[15]+0.223606797749979*alpha[6]*fr[13]-0.3872983346207417*alpha[3]*fr[12]+0.223606797749979*alpha[3]*fr[8]; 
  Ghat_r[14] = (-0.3872983346207417*alpha[8]*fr[30])-0.4330127018922194*alpha[0]*fr[30]+0.223606797749979*alpha[8]*fr[27]+0.25*alpha[0]*fr[27]-0.4330127018922193*alpha[2]*fr[26]+0.2500000000000001*alpha[2]*fr[24]-0.3464101615137754*alpha[6]*fr[22]+0.2*alpha[6]*fr[19]-0.3872983346207417*alpha[3]*fr[14]+0.223606797749979*alpha[3]*fr[10]-0.3872983346207417*alpha[6]*fr[9]+0.223606797749979*fr[4]*alpha[6]; 
  Ghat_r[15] = (-0.3872983346207417*alpha[8]*fr[31])-0.4330127018922194*alpha[0]*fr[31]+0.223606797749979*alpha[8]*fr[29]+0.25*alpha[0]*fr[29]-0.4330127018922193*alpha[2]*fr[28]+0.2500000000000001*alpha[2]*fr[25]-0.3464101615137754*alpha[6]*fr[23]+0.2*alpha[6]*fr[21]-0.3872983346207416*alpha[3]*fr[15]+0.223606797749979*alpha[3]*fr[13]-0.3872983346207416*alpha[6]*fr[12]+0.223606797749979*alpha[6]*fr[8]; 

  Ghat_l[0] = (-0.4330127018922194*alpha[8]*fc[18])+0.25*alpha[8]*fc[16]-0.4330127018922193*alpha[6]*fc[14]+0.25*alpha[6]*fc[10]-0.4330127018922193*alpha[3]*fc[9]-0.4330127018922193*alpha[2]*fc[7]+0.25*alpha[3]*fc[4]+0.25*alpha[2]*fc[3]-0.4330127018922193*alpha[0]*fc[2]+0.25*alpha[0]*fc[0]; 
  Ghat_l[1] = (-0.4330127018922193*alpha[8]*fc[20])+0.2500000000000001*alpha[8]*fc[17]-0.4330127018922193*alpha[6]*fc[15]+0.25*alpha[6]*fc[13]-0.4330127018922193*alpha[3]*fc[12]-0.4330127018922193*alpha[2]*fc[11]+0.25*alpha[3]*fc[8]+0.25*alpha[2]*fc[6]-0.4330127018922193*alpha[0]*fc[5]+0.25*alpha[0]*fc[1]; 
  Ghat_l[2] = (-0.3872983346207416*alpha[6]*fc[22])+0.223606797749979*alpha[6]*fc[19]-0.3872983346207417*alpha[2]*fc[18]+0.223606797749979*alpha[2]*fc[16]-0.4330127018922193*alpha[3]*fc[14]+0.25*alpha[3]*fc[10]-0.4330127018922193*alpha[6]*fc[9]-0.3872983346207416*fc[7]*alpha[8]+0.223606797749979*fc[3]*alpha[8]-0.4330127018922193*alpha[0]*fc[7]+0.25*fc[4]*alpha[6]+0.25*alpha[0]*fc[3]-0.4330127018922193*alpha[2]*fc[2]+0.25*fc[0]*alpha[2]; 
  Ghat_l[3] = (-0.3872983346207416*alpha[6]*fc[30])+0.223606797749979*alpha[6]*fc[27]-0.3872983346207417*alpha[3]*fc[26]+0.223606797749979*alpha[3]*fc[24]-0.4330127018922193*alpha[8]*fc[22]+0.2500000000000001*alpha[8]*fc[19]-0.4330127018922193*alpha[2]*fc[14]+0.25*alpha[2]*fc[10]-0.4330127018922193*alpha[0]*fc[9]-0.4330127018922193*alpha[6]*fc[7]+0.25*fc[3]*alpha[6]+0.25*alpha[0]*fc[4]-0.4330127018922193*fc[2]*alpha[3]+0.25*fc[0]*alpha[3]; 
  Ghat_l[4] = (-0.3872983346207417*alpha[6]*fc[23])+0.223606797749979*alpha[6]*fc[21]-0.3872983346207416*alpha[2]*fc[20]+0.223606797749979*alpha[2]*fc[17]-0.4330127018922193*alpha[3]*fc[15]+0.25*alpha[3]*fc[13]-0.4330127018922193*alpha[6]*fc[12]-0.3872983346207416*alpha[8]*fc[11]-0.4330127018922193*alpha[0]*fc[11]+0.25*alpha[6]*fc[8]+0.223606797749979*fc[6]*alpha[8]+0.25*alpha[0]*fc[6]-0.4330127018922193*alpha[2]*fc[5]+0.25*fc[1]*alpha[2]; 
  Ghat_l[5] = (-0.3872983346207417*alpha[6]*fc[31])+0.223606797749979*alpha[6]*fc[29]-0.3872983346207416*alpha[3]*fc[28]+0.223606797749979*alpha[3]*fc[25]-0.4330127018922194*alpha[8]*fc[23]+0.25*alpha[8]*fc[21]-0.4330127018922193*alpha[2]*fc[15]+0.25*alpha[2]*fc[13]-0.4330127018922193*alpha[0]*fc[12]-0.4330127018922193*alpha[6]*fc[11]+0.25*alpha[0]*fc[8]+0.25*alpha[6]*fc[6]-0.4330127018922193*alpha[3]*fc[5]+0.25*fc[1]*alpha[3]; 
  Ghat_l[6] = (-0.3872983346207416*alpha[3]*fc[30])+0.223606797749979*alpha[3]*fc[27]-0.3872983346207417*alpha[6]*fc[26]+0.223606797749979*alpha[6]*fc[24]-0.3872983346207416*alpha[2]*fc[22]+0.223606797749979*alpha[2]*fc[19]-0.3872983346207417*alpha[6]*fc[18]+0.223606797749979*alpha[6]*fc[16]-0.3872983346207416*alpha[8]*fc[14]-0.4330127018922193*alpha[0]*fc[14]+0.223606797749979*alpha[8]*fc[10]+0.25*alpha[0]*fc[10]-0.4330127018922193*alpha[2]*fc[9]-0.4330127018922193*alpha[3]*fc[7]-0.4330127018922193*fc[2]*alpha[6]+0.25*fc[0]*alpha[6]+0.25*alpha[2]*fc[4]+0.25*alpha[3]*fc[3]; 
  Ghat_l[7] = (-0.3872983346207417*alpha[3]*fc[31])+0.223606797749979*alpha[3]*fc[29]-0.3872983346207416*alpha[6]*fc[28]+0.223606797749979*alpha[6]*fc[25]-0.3872983346207417*alpha[2]*fc[23]+0.223606797749979*alpha[2]*fc[21]-0.3872983346207416*alpha[6]*fc[20]+0.223606797749979*alpha[6]*fc[17]-0.3872983346207416*alpha[8]*fc[15]-0.4330127018922193*alpha[0]*fc[15]+0.223606797749979*alpha[8]*fc[13]+0.25*alpha[0]*fc[13]-0.4330127018922193*alpha[2]*fc[12]-0.4330127018922193*alpha[3]*fc[11]+0.25*alpha[2]*fc[8]+0.25*alpha[3]*fc[6]-0.4330127018922193*fc[5]*alpha[6]+0.25*fc[1]*alpha[6]; 
  Ghat_l[8] = (-0.4330127018922193*alpha[3]*fc[22])+0.2500000000000001*alpha[3]*fc[19]-0.276641667586244*alpha[8]*fc[18]-0.4330127018922194*alpha[0]*fc[18]+0.159719141249985*alpha[8]*fc[16]+0.25*alpha[0]*fc[16]-0.3872983346207416*alpha[6]*fc[14]+0.223606797749979*alpha[6]*fc[10]-0.4330127018922193*fc[2]*alpha[8]+0.25*fc[0]*alpha[8]-0.3872983346207416*alpha[2]*fc[7]+0.223606797749979*alpha[2]*fc[3]; 
  Ghat_l[9] = (-0.4330127018922193*alpha[3]*fc[23])+0.2500000000000001*alpha[3]*fc[21]-0.276641667586244*alpha[8]*fc[20]-0.4330127018922194*alpha[0]*fc[20]+0.159719141249985*alpha[8]*fc[17]+0.25*alpha[0]*fc[17]-0.3872983346207417*alpha[6]*fc[15]+0.223606797749979*alpha[6]*fc[13]-0.3872983346207417*alpha[2]*fc[11]-0.4330127018922194*fc[5]*alpha[8]+0.2500000000000001*fc[1]*alpha[8]+0.223606797749979*alpha[2]*fc[6]; 
  Ghat_l[10] = (-0.3464101615137754*alpha[6]*fc[30])+0.2*alpha[6]*fc[27]-0.276641667586244*alpha[8]*fc[22]-0.4330127018922194*alpha[0]*fc[22]+0.159719141249985*alpha[8]*fc[19]+0.25*alpha[0]*fc[19]-0.4330127018922193*alpha[3]*fc[18]+0.2500000000000001*alpha[3]*fc[16]-0.3872983346207417*alpha[2]*fc[14]+0.223606797749979*alpha[2]*fc[10]-0.4330127018922194*alpha[8]*fc[9]+0.2500000000000001*fc[4]*alpha[8]-0.3872983346207417*alpha[6]*fc[7]+0.223606797749979*fc[3]*alpha[6]; 
  Ghat_l[11] = (-0.3464101615137754*alpha[6]*fc[31])+0.2*alpha[6]*fc[29]-0.276641667586244*alpha[8]*fc[23]-0.4330127018922194*alpha[0]*fc[23]+0.159719141249985*alpha[8]*fc[21]+0.25*alpha[0]*fc[21]-0.4330127018922193*alpha[3]*fc[20]+0.2500000000000001*alpha[3]*fc[17]-0.3872983346207416*alpha[2]*fc[15]+0.223606797749979*alpha[2]*fc[13]-0.4330127018922193*alpha[8]*fc[12]-0.3872983346207416*alpha[6]*fc[11]+0.25*alpha[8]*fc[8]+0.223606797749979*alpha[6]*fc[6]; 
  Ghat_l[12] = (-0.4330127018922193*alpha[2]*fc[30])+0.2500000000000001*alpha[2]*fc[27]-0.4330127018922194*alpha[0]*fc[26]+0.25*alpha[0]*fc[24]-0.3872983346207416*alpha[6]*fc[14]+0.223606797749979*alpha[6]*fc[10]-0.3872983346207416*alpha[3]*fc[9]+0.223606797749979*alpha[3]*fc[4]; 
  Ghat_l[13] = (-0.4330127018922193*alpha[2]*fc[31])+0.2500000000000001*alpha[2]*fc[29]-0.4330127018922194*alpha[0]*fc[28]+0.25*alpha[0]*fc[25]-0.3872983346207417*alpha[6]*fc[15]+0.223606797749979*alpha[6]*fc[13]-0.3872983346207417*alpha[3]*fc[12]+0.223606797749979*alpha[3]*fc[8]; 
  Ghat_l[14] = (-0.3872983346207417*alpha[8]*fc[30])-0.4330127018922194*alpha[0]*fc[30]+0.223606797749979*alpha[8]*fc[27]+0.25*alpha[0]*fc[27]-0.4330127018922193*alpha[2]*fc[26]+0.2500000000000001*alpha[2]*fc[24]-0.3464101615137754*alpha[6]*fc[22]+0.2*alpha[6]*fc[19]-0.3872983346207417*alpha[3]*fc[14]+0.223606797749979*alpha[3]*fc[10]-0.3872983346207417*alpha[6]*fc[9]+0.223606797749979*fc[4]*alpha[6]; 
  Ghat_l[15] = (-0.3872983346207417*alpha[8]*fc[31])-0.4330127018922194*alpha[0]*fc[31]+0.223606797749979*alpha[8]*fc[29]+0.25*alpha[0]*fc[29]-0.4330127018922193*alpha[2]*fc[28]+0.2500000000000001*alpha[2]*fc[25]-0.3464101615137754*alpha[6]*fc[23]+0.2*alpha[6]*fc[21]-0.3872983346207416*alpha[3]*fc[15]+0.223606797749979*alpha[3]*fc[13]-0.3872983346207416*alpha[6]*fc[12]+0.223606797749979*alpha[6]*fc[8]; 

  } 
  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dx11; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dx11; 
  out[2] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dx11; 
  out[3] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dx11; 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dx11; 
  out[5] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dx11; 
  out[6] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dx11; 
  out[7] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dx11; 
  out[8] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dx11; 
  out[9] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dx11; 
  out[10] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dx11; 
  out[11] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dx11; 
  out[12] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dx11; 
  out[13] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dx11; 
  out[14] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dx11; 
  out[15] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dx11; 
  out[16] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*dx11; 
  out[17] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*dx11; 
  out[18] += -1.224744871391589*(Ghat_r[8]+Ghat_l[8])*dx11; 
  out[19] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*dx11; 
  out[20] += -1.224744871391589*(Ghat_r[9]+Ghat_l[9])*dx11; 
  out[21] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*dx11; 
  out[22] += -1.224744871391589*(Ghat_r[10]+Ghat_l[10])*dx11; 
  out[23] += -1.224744871391589*(Ghat_r[11]+Ghat_l[11])*dx11; 
  out[24] += (0.7071067811865475*Ghat_l[12]-0.7071067811865475*Ghat_r[12])*dx11; 
  out[25] += (0.7071067811865475*Ghat_l[13]-0.7071067811865475*Ghat_r[13])*dx11; 
  out[26] += -1.224744871391589*(Ghat_r[12]+Ghat_l[12])*dx11; 
  out[27] += (0.7071067811865475*Ghat_l[14]-0.7071067811865475*Ghat_r[14])*dx11; 
  out[28] += -1.224744871391589*(Ghat_r[13]+Ghat_l[13])*dx11; 
  out[29] += (0.7071067811865475*Ghat_l[15]-0.7071067811865475*Ghat_r[15])*dx11; 
  out[30] += -1.224744871391589*(Ghat_r[14]+Ghat_l[14])*dx11; 
  out[31] += -1.224744871391589*(Ghat_r[15]+Ghat_l[15])*dx11; 

  return 0.;

} 
