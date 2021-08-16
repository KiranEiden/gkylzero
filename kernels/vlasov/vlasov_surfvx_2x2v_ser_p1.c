#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_surfvx_2x2v_ser_p1(const double *w, const double *dxv, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // qmem:      q/m*EM fields.
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 
  const double dv10 = 2/dxv[2]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double *E0 = &qmem[0]; 
  const double *B2 = &qmem[20]; 

  double alpha[8]; 

  alpha[0] = 1.414213562373095*(B2[0]*wv2+E0[0]); 
  alpha[1] = 1.414213562373095*(B2[1]*wv2+E0[1]); 
  alpha[2] = 1.414213562373095*(B2[2]*wv2+E0[2]); 
  alpha[3] = 0.408248290463863*B2[0]*dv2; 
  alpha[4] = 1.414213562373095*(B2[3]*wv2+E0[3]); 
  alpha[5] = 0.408248290463863*B2[1]*dv2; 
  alpha[6] = 0.408248290463863*B2[2]*dv2; 
  alpha[7] = 0.408248290463863*B2[3]*dv2; 

  double fUpwindQuad_l[8];
  double fUpwindQuad_r[8];
  if ((-alpha[7])+alpha[6]+alpha[5]+alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[0] = (-0.4330127018922193*fl[15])+0.4330127018922193*(fl[14]+fl[13])-0.25*fl[12]+0.4330127018922193*fl[11]-0.4330127018922193*fl[10]+0.25*(fl[9]+fl[8])-0.4330127018922193*(fl[7]+fl[6])+0.25*fl[5]-0.25*fl[4]+0.4330127018922193*fl[3]-0.25*(fl[2]+fl[1])+0.25*fl[0]; 
    fUpwindQuad_r[0] = (-0.4330127018922193*fc[15])+0.4330127018922193*(fc[14]+fc[13])-0.25*fc[12]+0.4330127018922193*fc[11]-0.4330127018922193*fc[10]+0.25*(fc[9]+fc[8])-0.4330127018922193*(fc[7]+fc[6])+0.25*fc[5]-0.25*fc[4]+0.4330127018922193*fc[3]-0.25*(fc[2]+fc[1])+0.25*fc[0]; 
  } else { 

    fUpwindQuad_l[0] = 0.4330127018922193*fc[15]-0.4330127018922193*(fc[14]+fc[13])-0.25*fc[12]-0.4330127018922193*fc[11]+0.4330127018922193*fc[10]+0.25*(fc[9]+fc[8])+0.4330127018922193*(fc[7]+fc[6])+0.25*fc[5]-0.25*fc[4]-0.4330127018922193*fc[3]-0.25*(fc[2]+fc[1])+0.25*fc[0]; 
    fUpwindQuad_r[0] = 0.4330127018922193*fr[15]-0.4330127018922193*(fr[14]+fr[13])-0.25*fr[12]-0.4330127018922193*fr[11]+0.4330127018922193*fr[10]+0.25*(fr[9]+fr[8])+0.4330127018922193*(fr[7]+fr[6])+0.25*fr[5]-0.25*fr[4]-0.4330127018922193*fr[3]-0.25*(fr[2]+fr[1])+0.25*fr[0]; 
  } 
  if (alpha[7]+alpha[6]-alpha[5]-alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[1] = 0.4330127018922193*(fl[15]+fl[14])-0.4330127018922193*fl[13]+0.25*fl[12]-0.4330127018922193*(fl[11]+fl[10])+0.25*fl[9]-0.25*fl[8]-0.4330127018922193*fl[7]+0.4330127018922193*fl[6]-0.25*(fl[5]+fl[4])+0.4330127018922193*fl[3]-0.25*fl[2]+0.25*(fl[1]+fl[0]); 
    fUpwindQuad_r[1] = 0.4330127018922193*(fc[15]+fc[14])-0.4330127018922193*fc[13]+0.25*fc[12]-0.4330127018922193*(fc[11]+fc[10])+0.25*fc[9]-0.25*fc[8]-0.4330127018922193*fc[7]+0.4330127018922193*fc[6]-0.25*(fc[5]+fc[4])+0.4330127018922193*fc[3]-0.25*fc[2]+0.25*(fc[1]+fc[0]); 
  } else { 

    fUpwindQuad_l[1] = (-0.4330127018922193*(fc[15]+fc[14]))+0.4330127018922193*fc[13]+0.25*fc[12]+0.4330127018922193*(fc[11]+fc[10])+0.25*fc[9]-0.25*fc[8]+0.4330127018922193*fc[7]-0.4330127018922193*fc[6]-0.25*(fc[5]+fc[4])-0.4330127018922193*fc[3]-0.25*fc[2]+0.25*(fc[1]+fc[0]); 
    fUpwindQuad_r[1] = (-0.4330127018922193*(fr[15]+fr[14]))+0.4330127018922193*fr[13]+0.25*fr[12]+0.4330127018922193*(fr[11]+fr[10])+0.25*fr[9]-0.25*fr[8]+0.4330127018922193*fr[7]-0.4330127018922193*fr[6]-0.25*(fr[5]+fr[4])-0.4330127018922193*fr[3]-0.25*fr[2]+0.25*(fr[1]+fr[0]); 
  } 
  if (alpha[7]-alpha[6]+alpha[5]-alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[2] = 0.4330127018922193*fl[15]-0.4330127018922193*fl[14]+0.4330127018922193*fl[13]+0.25*fl[12]-0.4330127018922193*(fl[11]+fl[10])-0.25*fl[9]+0.25*fl[8]+0.4330127018922193*fl[7]-0.4330127018922193*fl[6]-0.25*(fl[5]+fl[4])+0.4330127018922193*fl[3]+0.25*fl[2]-0.25*fl[1]+0.25*fl[0]; 
    fUpwindQuad_r[2] = 0.4330127018922193*fc[15]-0.4330127018922193*fc[14]+0.4330127018922193*fc[13]+0.25*fc[12]-0.4330127018922193*(fc[11]+fc[10])-0.25*fc[9]+0.25*fc[8]+0.4330127018922193*fc[7]-0.4330127018922193*fc[6]-0.25*(fc[5]+fc[4])+0.4330127018922193*fc[3]+0.25*fc[2]-0.25*fc[1]+0.25*fc[0]; 
  } else { 

    fUpwindQuad_l[2] = (-0.4330127018922193*fc[15])+0.4330127018922193*fc[14]-0.4330127018922193*fc[13]+0.25*fc[12]+0.4330127018922193*(fc[11]+fc[10])-0.25*fc[9]+0.25*fc[8]-0.4330127018922193*fc[7]+0.4330127018922193*fc[6]-0.25*(fc[5]+fc[4])-0.4330127018922193*fc[3]+0.25*fc[2]-0.25*fc[1]+0.25*fc[0]; 
    fUpwindQuad_r[2] = (-0.4330127018922193*fr[15])+0.4330127018922193*fr[14]-0.4330127018922193*fr[13]+0.25*fr[12]+0.4330127018922193*(fr[11]+fr[10])-0.25*fr[9]+0.25*fr[8]-0.4330127018922193*fr[7]+0.4330127018922193*fr[6]-0.25*(fr[5]+fr[4])-0.4330127018922193*fr[3]+0.25*fr[2]-0.25*fr[1]+0.25*fr[0]; 
  } 
  if ((-alpha[7])-alpha[6]-alpha[5]+alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[3] = (-0.4330127018922193*(fl[15]+fl[14]+fl[13]))-0.25*fl[12]+0.4330127018922193*fl[11]-0.4330127018922193*fl[10]-0.25*(fl[9]+fl[8])+0.4330127018922193*(fl[7]+fl[6])+0.25*fl[5]-0.25*fl[4]+0.4330127018922193*fl[3]+0.25*(fl[2]+fl[1]+fl[0]); 
    fUpwindQuad_r[3] = (-0.4330127018922193*(fc[15]+fc[14]+fc[13]))-0.25*fc[12]+0.4330127018922193*fc[11]-0.4330127018922193*fc[10]-0.25*(fc[9]+fc[8])+0.4330127018922193*(fc[7]+fc[6])+0.25*fc[5]-0.25*fc[4]+0.4330127018922193*fc[3]+0.25*(fc[2]+fc[1]+fc[0]); 
  } else { 

    fUpwindQuad_l[3] = 0.4330127018922193*(fc[15]+fc[14]+fc[13])-0.25*fc[12]-0.4330127018922193*fc[11]+0.4330127018922193*fc[10]-0.25*(fc[9]+fc[8])-0.4330127018922193*(fc[7]+fc[6])+0.25*fc[5]-0.25*fc[4]-0.4330127018922193*fc[3]+0.25*(fc[2]+fc[1]+fc[0]); 
    fUpwindQuad_r[3] = 0.4330127018922193*(fr[15]+fr[14]+fr[13])-0.25*fr[12]-0.4330127018922193*fr[11]+0.4330127018922193*fr[10]-0.25*(fr[9]+fr[8])-0.4330127018922193*(fr[7]+fr[6])+0.25*fr[5]-0.25*fr[4]-0.4330127018922193*fr[3]+0.25*(fr[2]+fr[1]+fr[0]); 
  } 
  if (alpha[7]-alpha[6]-alpha[5]+alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[4] = 0.4330127018922193*fl[15]-0.4330127018922193*(fl[14]+fl[13])+0.25*fl[12]+0.4330127018922193*(fl[11]+fl[10])-0.25*(fl[9]+fl[8])-0.4330127018922193*(fl[7]+fl[6])+0.25*(fl[5]+fl[4])+0.4330127018922193*fl[3]-0.25*(fl[2]+fl[1])+0.25*fl[0]; 
    fUpwindQuad_r[4] = 0.4330127018922193*fc[15]-0.4330127018922193*(fc[14]+fc[13])+0.25*fc[12]+0.4330127018922193*(fc[11]+fc[10])-0.25*(fc[9]+fc[8])-0.4330127018922193*(fc[7]+fc[6])+0.25*(fc[5]+fc[4])+0.4330127018922193*fc[3]-0.25*(fc[2]+fc[1])+0.25*fc[0]; 
  } else { 

    fUpwindQuad_l[4] = (-0.4330127018922193*fc[15])+0.4330127018922193*(fc[14]+fc[13])+0.25*fc[12]-0.4330127018922193*(fc[11]+fc[10])-0.25*(fc[9]+fc[8])+0.4330127018922193*(fc[7]+fc[6])+0.25*(fc[5]+fc[4])-0.4330127018922193*fc[3]-0.25*(fc[2]+fc[1])+0.25*fc[0]; 
    fUpwindQuad_r[4] = (-0.4330127018922193*fr[15])+0.4330127018922193*(fr[14]+fr[13])+0.25*fr[12]-0.4330127018922193*(fr[11]+fr[10])-0.25*(fr[9]+fr[8])+0.4330127018922193*(fr[7]+fr[6])+0.25*(fr[5]+fr[4])-0.4330127018922193*fr[3]-0.25*(fr[2]+fr[1])+0.25*fr[0]; 
  } 
  if ((-alpha[7])-alpha[6]+alpha[5]-alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[5] = (-0.4330127018922193*(fl[15]+fl[14]))+0.4330127018922193*fl[13]-0.25*fl[12]-0.4330127018922193*fl[11]+0.4330127018922193*fl[10]-0.25*fl[9]+0.25*fl[8]-0.4330127018922193*fl[7]+0.4330127018922193*fl[6]-0.25*fl[5]+0.25*fl[4]+0.4330127018922193*fl[3]-0.25*fl[2]+0.25*(fl[1]+fl[0]); 
    fUpwindQuad_r[5] = (-0.4330127018922193*(fc[15]+fc[14]))+0.4330127018922193*fc[13]-0.25*fc[12]-0.4330127018922193*fc[11]+0.4330127018922193*fc[10]-0.25*fc[9]+0.25*fc[8]-0.4330127018922193*fc[7]+0.4330127018922193*fc[6]-0.25*fc[5]+0.25*fc[4]+0.4330127018922193*fc[3]-0.25*fc[2]+0.25*(fc[1]+fc[0]); 
  } else { 

    fUpwindQuad_l[5] = 0.4330127018922193*(fc[15]+fc[14])-0.4330127018922193*fc[13]-0.25*fc[12]+0.4330127018922193*fc[11]-0.4330127018922193*fc[10]-0.25*fc[9]+0.25*fc[8]+0.4330127018922193*fc[7]-0.4330127018922193*fc[6]-0.25*fc[5]+0.25*fc[4]-0.4330127018922193*fc[3]-0.25*fc[2]+0.25*(fc[1]+fc[0]); 
    fUpwindQuad_r[5] = 0.4330127018922193*(fr[15]+fr[14])-0.4330127018922193*fr[13]-0.25*fr[12]+0.4330127018922193*fr[11]-0.4330127018922193*fr[10]-0.25*fr[9]+0.25*fr[8]+0.4330127018922193*fr[7]-0.4330127018922193*fr[6]-0.25*fr[5]+0.25*fr[4]-0.4330127018922193*fr[3]-0.25*fr[2]+0.25*(fr[1]+fr[0]); 
  } 
  if ((-alpha[7])+alpha[6]-alpha[5]-alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[6] = (-0.4330127018922193*fl[15])+0.4330127018922193*fl[14]-0.4330127018922193*fl[13]-0.25*fl[12]-0.4330127018922193*fl[11]+0.4330127018922193*fl[10]+0.25*fl[9]-0.25*fl[8]+0.4330127018922193*fl[7]-0.4330127018922193*fl[6]-0.25*fl[5]+0.25*fl[4]+0.4330127018922193*fl[3]+0.25*fl[2]-0.25*fl[1]+0.25*fl[0]; 
    fUpwindQuad_r[6] = (-0.4330127018922193*fc[15])+0.4330127018922193*fc[14]-0.4330127018922193*fc[13]-0.25*fc[12]-0.4330127018922193*fc[11]+0.4330127018922193*fc[10]+0.25*fc[9]-0.25*fc[8]+0.4330127018922193*fc[7]-0.4330127018922193*fc[6]-0.25*fc[5]+0.25*fc[4]+0.4330127018922193*fc[3]+0.25*fc[2]-0.25*fc[1]+0.25*fc[0]; 
  } else { 

    fUpwindQuad_l[6] = 0.4330127018922193*fc[15]-0.4330127018922193*fc[14]+0.4330127018922193*fc[13]-0.25*fc[12]+0.4330127018922193*fc[11]-0.4330127018922193*fc[10]+0.25*fc[9]-0.25*fc[8]-0.4330127018922193*fc[7]+0.4330127018922193*fc[6]-0.25*fc[5]+0.25*fc[4]-0.4330127018922193*fc[3]+0.25*fc[2]-0.25*fc[1]+0.25*fc[0]; 
    fUpwindQuad_r[6] = 0.4330127018922193*fr[15]-0.4330127018922193*fr[14]+0.4330127018922193*fr[13]-0.25*fr[12]+0.4330127018922193*fr[11]-0.4330127018922193*fr[10]+0.25*fr[9]-0.25*fr[8]-0.4330127018922193*fr[7]+0.4330127018922193*fr[6]-0.25*fr[5]+0.25*fr[4]-0.4330127018922193*fr[3]+0.25*fr[2]-0.25*fr[1]+0.25*fr[0]; 
  } 
  if (alpha[7]+alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[7] = 0.4330127018922193*(fl[15]+fl[14]+fl[13])+0.25*fl[12]+0.4330127018922193*(fl[11]+fl[10])+0.25*(fl[9]+fl[8])+0.4330127018922193*(fl[7]+fl[6])+0.25*(fl[5]+fl[4])+0.4330127018922193*fl[3]+0.25*(fl[2]+fl[1]+fl[0]); 
    fUpwindQuad_r[7] = 0.4330127018922193*(fc[15]+fc[14]+fc[13])+0.25*fc[12]+0.4330127018922193*(fc[11]+fc[10])+0.25*(fc[9]+fc[8])+0.4330127018922193*(fc[7]+fc[6])+0.25*(fc[5]+fc[4])+0.4330127018922193*fc[3]+0.25*(fc[2]+fc[1]+fc[0]); 
  } else { 

    fUpwindQuad_l[7] = (-0.4330127018922193*(fc[15]+fc[14]+fc[13]))+0.25*fc[12]-0.4330127018922193*(fc[11]+fc[10])+0.25*(fc[9]+fc[8])-0.4330127018922193*(fc[7]+fc[6])+0.25*(fc[5]+fc[4])-0.4330127018922193*fc[3]+0.25*(fc[2]+fc[1]+fc[0]); 
    fUpwindQuad_r[7] = (-0.4330127018922193*(fr[15]+fr[14]+fr[13]))+0.25*fr[12]-0.4330127018922193*(fr[11]+fr[10])+0.25*(fr[9]+fr[8])-0.4330127018922193*(fr[7]+fr[6])+0.25*(fr[5]+fr[4])-0.4330127018922193*fr[3]+0.25*(fr[2]+fr[1]+fr[0]); 
  } 
  double fUpwind_l[8];
  fUpwind_l[0] = 0.3535533905932737*(fUpwindQuad_l[7]+fUpwindQuad_l[6]+fUpwindQuad_l[5]+fUpwindQuad_l[4]+fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[1] = 0.3535533905932737*(fUpwindQuad_l[7]-1.0*fUpwindQuad_l[6]+fUpwindQuad_l[5]-1.0*fUpwindQuad_l[4]+fUpwindQuad_l[3]-1.0*fUpwindQuad_l[2]+fUpwindQuad_l[1]-1.0*fUpwindQuad_l[0]); 
  fUpwind_l[2] = 0.3535533905932737*(fUpwindQuad_l[7]+fUpwindQuad_l[6]-1.0*(fUpwindQuad_l[5]+fUpwindQuad_l[4])+fUpwindQuad_l[3]+fUpwindQuad_l[2]-1.0*(fUpwindQuad_l[1]+fUpwindQuad_l[0])); 
  fUpwind_l[3] = 0.3535533905932737*(fUpwindQuad_l[7]+fUpwindQuad_l[6]+fUpwindQuad_l[5]+fUpwindQuad_l[4]-1.0*(fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]+fUpwindQuad_l[0])); 
  fUpwind_l[4] = 0.3535533905932737*(fUpwindQuad_l[7]-1.0*(fUpwindQuad_l[6]+fUpwindQuad_l[5])+fUpwindQuad_l[4]+fUpwindQuad_l[3]-1.0*(fUpwindQuad_l[2]+fUpwindQuad_l[1])+fUpwindQuad_l[0]); 
  fUpwind_l[5] = 0.3535533905932737*(fUpwindQuad_l[7]-1.0*fUpwindQuad_l[6]+fUpwindQuad_l[5]-1.0*(fUpwindQuad_l[4]+fUpwindQuad_l[3])+fUpwindQuad_l[2]-1.0*fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[6] = 0.3535533905932737*(fUpwindQuad_l[7]+fUpwindQuad_l[6]-1.0*(fUpwindQuad_l[5]+fUpwindQuad_l[4]+fUpwindQuad_l[3]+fUpwindQuad_l[2])+fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[7] = 0.3535533905932737*(fUpwindQuad_l[7]-1.0*(fUpwindQuad_l[6]+fUpwindQuad_l[5])+fUpwindQuad_l[4]-1.0*fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]-1.0*fUpwindQuad_l[0]); 

  double fUpwind_r[8];
  fUpwind_r[0] = 0.3535533905932737*(fUpwindQuad_r[7]+fUpwindQuad_r[6]+fUpwindQuad_r[5]+fUpwindQuad_r[4]+fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[1] = 0.3535533905932737*(fUpwindQuad_r[7]-1.0*fUpwindQuad_r[6]+fUpwindQuad_r[5]-1.0*fUpwindQuad_r[4]+fUpwindQuad_r[3]-1.0*fUpwindQuad_r[2]+fUpwindQuad_r[1]-1.0*fUpwindQuad_r[0]); 
  fUpwind_r[2] = 0.3535533905932737*(fUpwindQuad_r[7]+fUpwindQuad_r[6]-1.0*(fUpwindQuad_r[5]+fUpwindQuad_r[4])+fUpwindQuad_r[3]+fUpwindQuad_r[2]-1.0*(fUpwindQuad_r[1]+fUpwindQuad_r[0])); 
  fUpwind_r[3] = 0.3535533905932737*(fUpwindQuad_r[7]+fUpwindQuad_r[6]+fUpwindQuad_r[5]+fUpwindQuad_r[4]-1.0*(fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]+fUpwindQuad_r[0])); 
  fUpwind_r[4] = 0.3535533905932737*(fUpwindQuad_r[7]-1.0*(fUpwindQuad_r[6]+fUpwindQuad_r[5])+fUpwindQuad_r[4]+fUpwindQuad_r[3]-1.0*(fUpwindQuad_r[2]+fUpwindQuad_r[1])+fUpwindQuad_r[0]); 
  fUpwind_r[5] = 0.3535533905932737*(fUpwindQuad_r[7]-1.0*fUpwindQuad_r[6]+fUpwindQuad_r[5]-1.0*(fUpwindQuad_r[4]+fUpwindQuad_r[3])+fUpwindQuad_r[2]-1.0*fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[6] = 0.3535533905932737*(fUpwindQuad_r[7]+fUpwindQuad_r[6]-1.0*(fUpwindQuad_r[5]+fUpwindQuad_r[4]+fUpwindQuad_r[3]+fUpwindQuad_r[2])+fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[7] = 0.3535533905932737*(fUpwindQuad_r[7]-1.0*(fUpwindQuad_r[6]+fUpwindQuad_r[5])+fUpwindQuad_r[4]-1.0*fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]-1.0*fUpwindQuad_r[0]); 

  double Ghat_l[8]; 
  double Ghat_r[8]; 
  Ghat_l[0] = 0.3535533905932737*(alpha[7]*fUpwind_l[7]+alpha[6]*fUpwind_l[6]+alpha[5]*fUpwind_l[5]+alpha[4]*fUpwind_l[4]+alpha[3]*fUpwind_l[3]+alpha[2]*fUpwind_l[2]+alpha[1]*fUpwind_l[1]+alpha[0]*fUpwind_l[0]); 
  Ghat_l[1] = 0.3535533905932737*(alpha[6]*fUpwind_l[7]+fUpwind_l[6]*alpha[7]+alpha[3]*fUpwind_l[5]+fUpwind_l[3]*alpha[5]+alpha[2]*fUpwind_l[4]+fUpwind_l[2]*alpha[4]+alpha[0]*fUpwind_l[1]+fUpwind_l[0]*alpha[1]); 
  Ghat_l[2] = 0.3535533905932737*(alpha[5]*fUpwind_l[7]+fUpwind_l[5]*alpha[7]+alpha[3]*fUpwind_l[6]+fUpwind_l[3]*alpha[6]+alpha[1]*fUpwind_l[4]+fUpwind_l[1]*alpha[4]+alpha[0]*fUpwind_l[2]+fUpwind_l[0]*alpha[2]); 
  Ghat_l[3] = 0.3535533905932737*(alpha[4]*fUpwind_l[7]+fUpwind_l[4]*alpha[7]+alpha[2]*fUpwind_l[6]+fUpwind_l[2]*alpha[6]+alpha[1]*fUpwind_l[5]+fUpwind_l[1]*alpha[5]+alpha[0]*fUpwind_l[3]+fUpwind_l[0]*alpha[3]); 
  Ghat_l[4] = 0.3535533905932737*(alpha[3]*fUpwind_l[7]+fUpwind_l[3]*alpha[7]+alpha[5]*fUpwind_l[6]+fUpwind_l[5]*alpha[6]+alpha[0]*fUpwind_l[4]+fUpwind_l[0]*alpha[4]+alpha[1]*fUpwind_l[2]+fUpwind_l[1]*alpha[2]); 
  Ghat_l[5] = 0.3535533905932737*(alpha[2]*fUpwind_l[7]+fUpwind_l[2]*alpha[7]+alpha[4]*fUpwind_l[6]+fUpwind_l[4]*alpha[6]+alpha[0]*fUpwind_l[5]+fUpwind_l[0]*alpha[5]+alpha[1]*fUpwind_l[3]+fUpwind_l[1]*alpha[3]); 
  Ghat_l[6] = 0.3535533905932737*(alpha[1]*fUpwind_l[7]+fUpwind_l[1]*alpha[7]+alpha[0]*fUpwind_l[6]+fUpwind_l[0]*alpha[6]+alpha[4]*fUpwind_l[5]+fUpwind_l[4]*alpha[5]+alpha[2]*fUpwind_l[3]+fUpwind_l[2]*alpha[3]); 
  Ghat_l[7] = 0.3535533905932737*(alpha[0]*fUpwind_l[7]+fUpwind_l[0]*alpha[7]+alpha[1]*fUpwind_l[6]+fUpwind_l[1]*alpha[6]+alpha[2]*fUpwind_l[5]+fUpwind_l[2]*alpha[5]+alpha[3]*fUpwind_l[4]+fUpwind_l[3]*alpha[4]); 

  Ghat_r[0] = 0.3535533905932737*(alpha[7]*fUpwind_r[7]+alpha[6]*fUpwind_r[6]+alpha[5]*fUpwind_r[5]+alpha[4]*fUpwind_r[4]+alpha[3]*fUpwind_r[3]+alpha[2]*fUpwind_r[2]+alpha[1]*fUpwind_r[1]+alpha[0]*fUpwind_r[0]); 
  Ghat_r[1] = 0.3535533905932737*(alpha[6]*fUpwind_r[7]+fUpwind_r[6]*alpha[7]+alpha[3]*fUpwind_r[5]+fUpwind_r[3]*alpha[5]+alpha[2]*fUpwind_r[4]+fUpwind_r[2]*alpha[4]+alpha[0]*fUpwind_r[1]+fUpwind_r[0]*alpha[1]); 
  Ghat_r[2] = 0.3535533905932737*(alpha[5]*fUpwind_r[7]+fUpwind_r[5]*alpha[7]+alpha[3]*fUpwind_r[6]+fUpwind_r[3]*alpha[6]+alpha[1]*fUpwind_r[4]+fUpwind_r[1]*alpha[4]+alpha[0]*fUpwind_r[2]+fUpwind_r[0]*alpha[2]); 
  Ghat_r[3] = 0.3535533905932737*(alpha[4]*fUpwind_r[7]+fUpwind_r[4]*alpha[7]+alpha[2]*fUpwind_r[6]+fUpwind_r[2]*alpha[6]+alpha[1]*fUpwind_r[5]+fUpwind_r[1]*alpha[5]+alpha[0]*fUpwind_r[3]+fUpwind_r[0]*alpha[3]); 
  Ghat_r[4] = 0.3535533905932737*(alpha[3]*fUpwind_r[7]+fUpwind_r[3]*alpha[7]+alpha[5]*fUpwind_r[6]+fUpwind_r[5]*alpha[6]+alpha[0]*fUpwind_r[4]+fUpwind_r[0]*alpha[4]+alpha[1]*fUpwind_r[2]+fUpwind_r[1]*alpha[2]); 
  Ghat_r[5] = 0.3535533905932737*(alpha[2]*fUpwind_r[7]+fUpwind_r[2]*alpha[7]+alpha[4]*fUpwind_r[6]+fUpwind_r[4]*alpha[6]+alpha[0]*fUpwind_r[5]+fUpwind_r[0]*alpha[5]+alpha[1]*fUpwind_r[3]+fUpwind_r[1]*alpha[3]); 
  Ghat_r[6] = 0.3535533905932737*(alpha[1]*fUpwind_r[7]+fUpwind_r[1]*alpha[7]+alpha[0]*fUpwind_r[6]+fUpwind_r[0]*alpha[6]+alpha[4]*fUpwind_r[5]+fUpwind_r[4]*alpha[5]+alpha[2]*fUpwind_r[3]+fUpwind_r[2]*alpha[3]); 
  Ghat_r[7] = 0.3535533905932737*(alpha[0]*fUpwind_r[7]+fUpwind_r[0]*alpha[7]+alpha[1]*fUpwind_r[6]+fUpwind_r[1]*alpha[6]+alpha[2]*fUpwind_r[5]+fUpwind_r[2]*alpha[5]+alpha[3]*fUpwind_r[4]+fUpwind_r[3]*alpha[4]); 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv10; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv10; 
  out[2] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv10; 
  out[3] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv10; 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv10; 
  out[5] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dv10; 
  out[6] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv10; 
  out[7] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv10; 
  out[8] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dv10; 
  out[9] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dv10; 
  out[10] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv10; 
  out[11] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dv10; 
  out[12] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dv10; 
  out[13] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dv10; 
  out[14] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dv10; 
  out[15] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dv10; 

} 
