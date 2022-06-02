#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_sr_surfy_2x2v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // p_over_gamma: p/gamma (velocity).
  // fl/fc/fr:  Input Distribution function in left/center/right cells.
  // out:       Incremented distribution function in center cell.
  const double *p1_over_gamma = &p_over_gamma[4]; 
  double alpha[8] = {0.0}; 
  alpha[0] = 1.414213562373095*p1_over_gamma[0]; 
  alpha[2] = 1.414213562373095*p1_over_gamma[1]; 
  alpha[3] = 1.414213562373095*p1_over_gamma[2]; 
  alpha[6] = 1.414213562373095*p1_over_gamma[3]; 

  const double dx11 = 2/dxv[1]; 
  const double dv = dxv[3], wv = w[3]; 
  double Ghat_r[8]; 
  double Ghat_l[8]; 
  if (wv>0) { 

  Ghat_r[0] = alpha[6]*(0.4330127018922193*fc[14]+0.25*fc[10])+0.4330127018922193*(alpha[3]*fc[9]+alpha[2]*fc[7])+0.25*(alpha[3]*fc[4]+alpha[2]*fc[3])+alpha[0]*(0.4330127018922193*fc[2]+0.25*fc[0]); 
  Ghat_r[1] = alpha[6]*(0.4330127018922193*fc[15]+0.25*fc[13])+0.4330127018922193*(alpha[3]*fc[12]+alpha[2]*fc[11])+0.25*(alpha[3]*fc[8]+alpha[2]*fc[6])+alpha[0]*(0.4330127018922193*fc[5]+0.25*fc[1]); 
  Ghat_r[2] = alpha[3]*(0.4330127018922193*fc[14]+0.25*fc[10])+0.4330127018922193*(alpha[6]*fc[9]+alpha[0]*fc[7])+0.25*(fc[4]*alpha[6]+alpha[0]*fc[3])+alpha[2]*(0.4330127018922193*fc[2]+0.25*fc[0]); 
  Ghat_r[3] = alpha[2]*(0.4330127018922193*fc[14]+0.25*fc[10])+0.4330127018922193*(alpha[0]*fc[9]+alpha[6]*fc[7])+0.25*(fc[3]*alpha[6]+alpha[0]*fc[4])+(0.4330127018922193*fc[2]+0.25*fc[0])*alpha[3]; 
  Ghat_r[4] = alpha[3]*(0.4330127018922193*fc[15]+0.25*fc[13])+0.4330127018922193*(alpha[6]*fc[12]+alpha[0]*fc[11])+0.25*(alpha[6]*fc[8]+alpha[0]*fc[6])+alpha[2]*(0.4330127018922193*fc[5]+0.25*fc[1]); 
  Ghat_r[5] = alpha[2]*(0.4330127018922193*fc[15]+0.25*fc[13])+0.4330127018922193*(alpha[0]*fc[12]+alpha[6]*fc[11])+0.25*(alpha[0]*fc[8]+alpha[6]*fc[6])+alpha[3]*(0.4330127018922193*fc[5]+0.25*fc[1]); 
  Ghat_r[6] = alpha[0]*(0.4330127018922193*fc[14]+0.25*fc[10])+0.4330127018922193*(alpha[2]*fc[9]+alpha[3]*fc[7]+fc[2]*alpha[6])+0.25*(fc[0]*alpha[6]+alpha[2]*fc[4]+alpha[3]*fc[3]); 
  Ghat_r[7] = alpha[0]*(0.4330127018922193*fc[15]+0.25*fc[13])+0.4330127018922193*(alpha[2]*fc[12]+alpha[3]*fc[11])+0.25*(alpha[2]*fc[8]+alpha[3]*fc[6])+(0.4330127018922193*fc[5]+0.25*fc[1])*alpha[6]; 

  Ghat_l[0] = alpha[6]*(0.4330127018922193*fl[14]+0.25*fl[10])+0.4330127018922193*(alpha[3]*fl[9]+alpha[2]*fl[7])+0.25*(alpha[3]*fl[4]+alpha[2]*fl[3])+alpha[0]*(0.4330127018922193*fl[2]+0.25*fl[0]); 
  Ghat_l[1] = alpha[6]*(0.4330127018922193*fl[15]+0.25*fl[13])+0.4330127018922193*(alpha[3]*fl[12]+alpha[2]*fl[11])+0.25*(alpha[3]*fl[8]+alpha[2]*fl[6])+alpha[0]*(0.4330127018922193*fl[5]+0.25*fl[1]); 
  Ghat_l[2] = alpha[3]*(0.4330127018922193*fl[14]+0.25*fl[10])+0.4330127018922193*(alpha[6]*fl[9]+alpha[0]*fl[7])+0.25*(fl[4]*alpha[6]+alpha[0]*fl[3])+alpha[2]*(0.4330127018922193*fl[2]+0.25*fl[0]); 
  Ghat_l[3] = alpha[2]*(0.4330127018922193*fl[14]+0.25*fl[10])+0.4330127018922193*(alpha[0]*fl[9]+alpha[6]*fl[7])+0.25*(fl[3]*alpha[6]+alpha[0]*fl[4])+(0.4330127018922193*fl[2]+0.25*fl[0])*alpha[3]; 
  Ghat_l[4] = alpha[3]*(0.4330127018922193*fl[15]+0.25*fl[13])+0.4330127018922193*(alpha[6]*fl[12]+alpha[0]*fl[11])+0.25*(alpha[6]*fl[8]+alpha[0]*fl[6])+alpha[2]*(0.4330127018922193*fl[5]+0.25*fl[1]); 
  Ghat_l[5] = alpha[2]*(0.4330127018922193*fl[15]+0.25*fl[13])+0.4330127018922193*(alpha[0]*fl[12]+alpha[6]*fl[11])+0.25*(alpha[0]*fl[8]+alpha[6]*fl[6])+alpha[3]*(0.4330127018922193*fl[5]+0.25*fl[1]); 
  Ghat_l[6] = alpha[0]*(0.4330127018922193*fl[14]+0.25*fl[10])+0.4330127018922193*(alpha[2]*fl[9]+alpha[3]*fl[7]+fl[2]*alpha[6])+0.25*(fl[0]*alpha[6]+alpha[2]*fl[4]+alpha[3]*fl[3]); 
  Ghat_l[7] = alpha[0]*(0.4330127018922193*fl[15]+0.25*fl[13])+0.4330127018922193*(alpha[2]*fl[12]+alpha[3]*fl[11])+0.25*(alpha[2]*fl[8]+alpha[3]*fl[6])+(0.4330127018922193*fl[5]+0.25*fl[1])*alpha[6]; 

  } else { 

  Ghat_r[0] = -0.25*(alpha[6]*(1.732050807568877*fr[14]-1.0*fr[10])+1.732050807568877*(alpha[3]*fr[9]+alpha[2]*fr[7])-1.0*(alpha[3]*fr[4]+alpha[2]*fr[3])+alpha[0]*(1.732050807568877*fr[2]-1.0*fr[0])); 
  Ghat_r[1] = -0.25*(alpha[6]*(1.732050807568877*fr[15]-1.0*fr[13])+1.732050807568877*(alpha[3]*fr[12]+alpha[2]*fr[11])-1.0*(alpha[3]*fr[8]+alpha[2]*fr[6])+alpha[0]*(1.732050807568877*fr[5]-1.0*fr[1])); 
  Ghat_r[2] = -0.25*(alpha[3]*(1.732050807568877*fr[14]-1.0*fr[10])+1.732050807568877*(alpha[6]*fr[9]+alpha[0]*fr[7])-1.0*(fr[4]*alpha[6]+alpha[0]*fr[3])+alpha[2]*(1.732050807568877*fr[2]-1.0*fr[0])); 
  Ghat_r[3] = -0.25*(alpha[2]*(1.732050807568877*fr[14]-1.0*fr[10])+1.732050807568877*(alpha[0]*fr[9]+alpha[6]*fr[7])-1.0*(fr[3]*alpha[6]+alpha[0]*fr[4])+(1.732050807568877*fr[2]-1.0*fr[0])*alpha[3]); 
  Ghat_r[4] = -0.25*(alpha[3]*(1.732050807568877*fr[15]-1.0*fr[13])+1.732050807568877*(alpha[6]*fr[12]+alpha[0]*fr[11])-1.0*(alpha[6]*fr[8]+alpha[0]*fr[6])+alpha[2]*(1.732050807568877*fr[5]-1.0*fr[1])); 
  Ghat_r[5] = -0.25*(alpha[2]*(1.732050807568877*fr[15]-1.0*fr[13])+1.732050807568877*(alpha[0]*fr[12]+alpha[6]*fr[11])-1.0*(alpha[0]*fr[8]+alpha[6]*fr[6])+alpha[3]*(1.732050807568877*fr[5]-1.0*fr[1])); 
  Ghat_r[6] = -0.25*(alpha[0]*(1.732050807568877*fr[14]-1.0*fr[10])+1.732050807568877*(alpha[2]*fr[9]+alpha[3]*fr[7])+(1.732050807568877*fr[2]-1.0*fr[0])*alpha[6]-1.0*(alpha[2]*fr[4]+alpha[3]*fr[3])); 
  Ghat_r[7] = -0.25*(alpha[0]*(1.732050807568877*fr[15]-1.0*fr[13])+1.732050807568877*(alpha[2]*fr[12]+alpha[3]*fr[11])-1.0*(alpha[2]*fr[8]+alpha[3]*fr[6])+(1.732050807568877*fr[5]-1.0*fr[1])*alpha[6]); 

  Ghat_l[0] = -0.25*(alpha[6]*(1.732050807568877*fc[14]-1.0*fc[10])+1.732050807568877*(alpha[3]*fc[9]+alpha[2]*fc[7])-1.0*(alpha[3]*fc[4]+alpha[2]*fc[3])+alpha[0]*(1.732050807568877*fc[2]-1.0*fc[0])); 
  Ghat_l[1] = -0.25*(alpha[6]*(1.732050807568877*fc[15]-1.0*fc[13])+1.732050807568877*(alpha[3]*fc[12]+alpha[2]*fc[11])-1.0*(alpha[3]*fc[8]+alpha[2]*fc[6])+alpha[0]*(1.732050807568877*fc[5]-1.0*fc[1])); 
  Ghat_l[2] = -0.25*(alpha[3]*(1.732050807568877*fc[14]-1.0*fc[10])+1.732050807568877*(alpha[6]*fc[9]+alpha[0]*fc[7])-1.0*(fc[4]*alpha[6]+alpha[0]*fc[3])+alpha[2]*(1.732050807568877*fc[2]-1.0*fc[0])); 
  Ghat_l[3] = -0.25*(alpha[2]*(1.732050807568877*fc[14]-1.0*fc[10])+1.732050807568877*(alpha[0]*fc[9]+alpha[6]*fc[7])-1.0*(fc[3]*alpha[6]+alpha[0]*fc[4])+(1.732050807568877*fc[2]-1.0*fc[0])*alpha[3]); 
  Ghat_l[4] = -0.25*(alpha[3]*(1.732050807568877*fc[15]-1.0*fc[13])+1.732050807568877*(alpha[6]*fc[12]+alpha[0]*fc[11])-1.0*(alpha[6]*fc[8]+alpha[0]*fc[6])+alpha[2]*(1.732050807568877*fc[5]-1.0*fc[1])); 
  Ghat_l[5] = -0.25*(alpha[2]*(1.732050807568877*fc[15]-1.0*fc[13])+1.732050807568877*(alpha[0]*fc[12]+alpha[6]*fc[11])-1.0*(alpha[0]*fc[8]+alpha[6]*fc[6])+alpha[3]*(1.732050807568877*fc[5]-1.0*fc[1])); 
  Ghat_l[6] = -0.25*(alpha[0]*(1.732050807568877*fc[14]-1.0*fc[10])+1.732050807568877*(alpha[2]*fc[9]+alpha[3]*fc[7])+(1.732050807568877*fc[2]-1.0*fc[0])*alpha[6]-1.0*(alpha[2]*fc[4]+alpha[3]*fc[3])); 
  Ghat_l[7] = -0.25*(alpha[0]*(1.732050807568877*fc[15]-1.0*fc[13])+1.732050807568877*(alpha[2]*fc[12]+alpha[3]*fc[11])-1.0*(alpha[2]*fc[8]+alpha[3]*fc[6])+(1.732050807568877*fc[5]-1.0*fc[1])*alpha[6]); 

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
} 