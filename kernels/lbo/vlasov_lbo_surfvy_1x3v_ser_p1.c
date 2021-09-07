#include <gkyl_vlasov_lbo_kernels.h> 
#include <gkyl_basis_ser_1x3v_p1_surfvy_quad.h> 
GKYL_CU_DH void vlasov_lbo_surfvy_1x3v_ser_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[4]:         cell-center coordinates. 
  // dxv[4]:       cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[6]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      distribution function in cells 
  // out:           incremented distribution function in cell 
  double rdv2 = 2.0/dxv[2]; 
  double rdvSq4 = 4.0/(dxv[2]*dxv[2]); 

  const double *sumNuUy = &nuUSum[2]; 

  double alphaDrSurf_l[8] = {0.0}; 
  alphaDrSurf_l[0] = 2.828427124746191*w[2]*nuSum-1.414213562373095*dxv[2]*nuSum-2.0*sumNuUy[0]; 
  alphaDrSurf_l[1] = -2.0*sumNuUy[1]; 

  double alphaDrSurf_r[8] = {0.0}; 
  alphaDrSurf_r[0] = 2.828427124746191*w[2]*nuSum+1.414213562373095*dxv[2]*nuSum-2.0*sumNuUy[0]; 
  alphaDrSurf_r[1] = -2.0*sumNuUy[1]; 

  double fUpwindQuad_l[8] = {0.0};
  double fUpwindQuad_r[8] = {0.0};
  double fUpwind_l[8] = {0.0};;
  double fUpwind_r[8] = {0.0};
  double Ghat_l[8] = {0.0}; 
  double Ghat_r[8] = {0.0}; 
  double Gdiff_l[8] = {0.0}; 
  double Gdiff_r[8] = {0.0}; 
  double Gdiff2_l[8] = {0.0}; 
  double Gdiff2_r[8] = {0.0}; 

  if (alphaDrSurf_l[0]-alphaDrSurf_l[1] > 0) { 
    fUpwindQuad_l[0] = ser_1x3v_p1_surfvy_quad_0(1, fl); 
  } else { 
    fUpwindQuad_l[0] = ser_1x3v_p1_surfvy_quad_0(-1, fc); 
  } 
  if (alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[1] = ser_1x3v_p1_surfvy_quad_1(1, fl); 
  } else { 
    fUpwindQuad_l[1] = ser_1x3v_p1_surfvy_quad_1(-1, fc); 
  } 
  if (alphaDrSurf_l[0]-alphaDrSurf_l[1] > 0) { 
    fUpwindQuad_l[2] = ser_1x3v_p1_surfvy_quad_2(1, fl); 
  } else { 
    fUpwindQuad_l[2] = ser_1x3v_p1_surfvy_quad_2(-1, fc); 
  } 
  if (alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[3] = ser_1x3v_p1_surfvy_quad_3(1, fl); 
  } else { 
    fUpwindQuad_l[3] = ser_1x3v_p1_surfvy_quad_3(-1, fc); 
  } 
  if (alphaDrSurf_l[0]-alphaDrSurf_l[1] > 0) { 
    fUpwindQuad_l[4] = ser_1x3v_p1_surfvy_quad_4(1, fl); 
  } else { 
    fUpwindQuad_l[4] = ser_1x3v_p1_surfvy_quad_4(-1, fc); 
  } 
  if (alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[5] = ser_1x3v_p1_surfvy_quad_5(1, fl); 
  } else { 
    fUpwindQuad_l[5] = ser_1x3v_p1_surfvy_quad_5(-1, fc); 
  } 
  if (alphaDrSurf_l[0]-alphaDrSurf_l[1] > 0) { 
    fUpwindQuad_l[6] = ser_1x3v_p1_surfvy_quad_6(1, fl); 
  } else { 
    fUpwindQuad_l[6] = ser_1x3v_p1_surfvy_quad_6(-1, fc); 
  } 
  if (alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[7] = ser_1x3v_p1_surfvy_quad_7(1, fl); 
  } else { 
    fUpwindQuad_l[7] = ser_1x3v_p1_surfvy_quad_7(-1, fc); 
  } 
  if (alphaDrSurf_r[0]-alphaDrSurf_r[1] > 0) { 
    fUpwindQuad_r[0] = ser_1x3v_p1_surfvy_quad_0(1, fc); 
  } else { 
    fUpwindQuad_r[0] = ser_1x3v_p1_surfvy_quad_0(-1, fr); 
  } 
  if (alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[1] = ser_1x3v_p1_surfvy_quad_1(1, fc); 
  } else { 
    fUpwindQuad_r[1] = ser_1x3v_p1_surfvy_quad_1(-1, fr); 
  } 
  if (alphaDrSurf_r[0]-alphaDrSurf_r[1] > 0) { 
    fUpwindQuad_r[2] = ser_1x3v_p1_surfvy_quad_2(1, fc); 
  } else { 
    fUpwindQuad_r[2] = ser_1x3v_p1_surfvy_quad_2(-1, fr); 
  } 
  if (alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[3] = ser_1x3v_p1_surfvy_quad_3(1, fc); 
  } else { 
    fUpwindQuad_r[3] = ser_1x3v_p1_surfvy_quad_3(-1, fr); 
  } 
  if (alphaDrSurf_r[0]-alphaDrSurf_r[1] > 0) { 
    fUpwindQuad_r[4] = ser_1x3v_p1_surfvy_quad_4(1, fc); 
  } else { 
    fUpwindQuad_r[4] = ser_1x3v_p1_surfvy_quad_4(-1, fr); 
  } 
  if (alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[5] = ser_1x3v_p1_surfvy_quad_5(1, fc); 
  } else { 
    fUpwindQuad_r[5] = ser_1x3v_p1_surfvy_quad_5(-1, fr); 
  } 
  if (alphaDrSurf_r[0]-alphaDrSurf_r[1] > 0) { 
    fUpwindQuad_r[6] = ser_1x3v_p1_surfvy_quad_6(1, fc); 
  } else { 
    fUpwindQuad_r[6] = ser_1x3v_p1_surfvy_quad_6(-1, fr); 
  } 
  if (alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[7] = ser_1x3v_p1_surfvy_quad_7(1, fc); 
  } else { 
    fUpwindQuad_r[7] = ser_1x3v_p1_surfvy_quad_7(-1, fr); 
  } 
  fUpwind_l[0] = 0.3535533905932737*(fUpwindQuad_l[7]+fUpwindQuad_l[6]+fUpwindQuad_l[5]+fUpwindQuad_l[4]+fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[1] = 0.3535533905932737*(fUpwindQuad_l[7]-1.0*fUpwindQuad_l[6]+fUpwindQuad_l[5]-1.0*fUpwindQuad_l[4]+fUpwindQuad_l[3]-1.0*fUpwindQuad_l[2]+fUpwindQuad_l[1]-1.0*fUpwindQuad_l[0]); 
  fUpwind_l[2] = 0.3535533905932737*(fUpwindQuad_l[7]+fUpwindQuad_l[6]-1.0*(fUpwindQuad_l[5]+fUpwindQuad_l[4])+fUpwindQuad_l[3]+fUpwindQuad_l[2]-1.0*(fUpwindQuad_l[1]+fUpwindQuad_l[0])); 
  fUpwind_l[3] = 0.3535533905932737*(fUpwindQuad_l[7]+fUpwindQuad_l[6]+fUpwindQuad_l[5]+fUpwindQuad_l[4]-1.0*(fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]+fUpwindQuad_l[0])); 
  fUpwind_l[4] = 0.3535533905932737*(fUpwindQuad_l[7]-1.0*(fUpwindQuad_l[6]+fUpwindQuad_l[5])+fUpwindQuad_l[4]+fUpwindQuad_l[3]-1.0*(fUpwindQuad_l[2]+fUpwindQuad_l[1])+fUpwindQuad_l[0]); 
  fUpwind_l[5] = 0.3535533905932737*(fUpwindQuad_l[7]-1.0*fUpwindQuad_l[6]+fUpwindQuad_l[5]-1.0*(fUpwindQuad_l[4]+fUpwindQuad_l[3])+fUpwindQuad_l[2]-1.0*fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[6] = 0.3535533905932737*(fUpwindQuad_l[7]+fUpwindQuad_l[6]-1.0*(fUpwindQuad_l[5]+fUpwindQuad_l[4]+fUpwindQuad_l[3]+fUpwindQuad_l[2])+fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[7] = 0.3535533905932737*(fUpwindQuad_l[7]-1.0*(fUpwindQuad_l[6]+fUpwindQuad_l[5])+fUpwindQuad_l[4]-1.0*fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]-1.0*fUpwindQuad_l[0]); 

  fUpwind_r[0] = 0.3535533905932737*(fUpwindQuad_r[7]+fUpwindQuad_r[6]+fUpwindQuad_r[5]+fUpwindQuad_r[4]+fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[1] = 0.3535533905932737*(fUpwindQuad_r[7]-1.0*fUpwindQuad_r[6]+fUpwindQuad_r[5]-1.0*fUpwindQuad_r[4]+fUpwindQuad_r[3]-1.0*fUpwindQuad_r[2]+fUpwindQuad_r[1]-1.0*fUpwindQuad_r[0]); 
  fUpwind_r[2] = 0.3535533905932737*(fUpwindQuad_r[7]+fUpwindQuad_r[6]-1.0*(fUpwindQuad_r[5]+fUpwindQuad_r[4])+fUpwindQuad_r[3]+fUpwindQuad_r[2]-1.0*(fUpwindQuad_r[1]+fUpwindQuad_r[0])); 
  fUpwind_r[3] = 0.3535533905932737*(fUpwindQuad_r[7]+fUpwindQuad_r[6]+fUpwindQuad_r[5]+fUpwindQuad_r[4]-1.0*(fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]+fUpwindQuad_r[0])); 
  fUpwind_r[4] = 0.3535533905932737*(fUpwindQuad_r[7]-1.0*(fUpwindQuad_r[6]+fUpwindQuad_r[5])+fUpwindQuad_r[4]+fUpwindQuad_r[3]-1.0*(fUpwindQuad_r[2]+fUpwindQuad_r[1])+fUpwindQuad_r[0]); 
  fUpwind_r[5] = 0.3535533905932737*(fUpwindQuad_r[7]-1.0*fUpwindQuad_r[6]+fUpwindQuad_r[5]-1.0*(fUpwindQuad_r[4]+fUpwindQuad_r[3])+fUpwindQuad_r[2]-1.0*fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[6] = 0.3535533905932737*(fUpwindQuad_r[7]+fUpwindQuad_r[6]-1.0*(fUpwindQuad_r[5]+fUpwindQuad_r[4]+fUpwindQuad_r[3]+fUpwindQuad_r[2])+fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[7] = 0.3535533905932737*(fUpwindQuad_r[7]-1.0*(fUpwindQuad_r[6]+fUpwindQuad_r[5])+fUpwindQuad_r[4]-1.0*fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]-1.0*fUpwindQuad_r[0]); 

  Gdiff2_l[0] = 0.2886751345948129*nuVtSqSum[1]*fl[6]-0.2886751345948129*nuVtSqSum[1]*fc[6]+0.2886751345948129*nuVtSqSum[0]*fl[3]-0.2886751345948129*nuVtSqSum[0]*fc[3]+0.25*fl[1]*nuVtSqSum[1]+0.25*fc[1]*nuVtSqSum[1]+0.25*fl[0]*nuVtSqSum[0]+0.25*fc[0]*nuVtSqSum[0]; 
  Gdiff2_l[1] = 0.2886751345948129*nuVtSqSum[0]*fl[6]-0.2886751345948129*nuVtSqSum[0]*fc[6]+0.2886751345948129*nuVtSqSum[1]*fl[3]-0.2886751345948129*nuVtSqSum[1]*fc[3]+0.25*fl[0]*nuVtSqSum[1]+0.25*fc[0]*nuVtSqSum[1]+0.25*nuVtSqSum[0]*fl[1]+0.25*nuVtSqSum[0]*fc[1]; 
  Gdiff2_l[2] = 0.2886751345948129*nuVtSqSum[1]*fl[11]-0.2886751345948129*nuVtSqSum[1]*fc[11]+0.2886751345948129*nuVtSqSum[0]*fl[7]-0.2886751345948129*nuVtSqSum[0]*fc[7]+0.25*nuVtSqSum[1]*fl[5]+0.25*nuVtSqSum[1]*fc[5]+0.25*nuVtSqSum[0]*fl[2]+0.25*nuVtSqSum[0]*fc[2]; 
  Gdiff2_l[3] = 0.2886751345948129*nuVtSqSum[1]*fl[13]-0.2886751345948129*nuVtSqSum[1]*fc[13]+0.2886751345948129*nuVtSqSum[0]*fl[10]-0.2886751345948129*nuVtSqSum[0]*fc[10]+0.25*nuVtSqSum[1]*fl[8]+0.25*nuVtSqSum[1]*fc[8]+0.25*nuVtSqSum[0]*fl[4]+0.25*nuVtSqSum[0]*fc[4]; 
  Gdiff2_l[4] = 0.2886751345948129*nuVtSqSum[0]*fl[11]-0.2886751345948129*nuVtSqSum[0]*fc[11]+0.2886751345948129*nuVtSqSum[1]*fl[7]-0.2886751345948129*nuVtSqSum[1]*fc[7]+0.25*nuVtSqSum[0]*fl[5]+0.25*nuVtSqSum[0]*fc[5]+0.25*nuVtSqSum[1]*fl[2]+0.25*nuVtSqSum[1]*fc[2]; 
  Gdiff2_l[5] = 0.2886751345948129*nuVtSqSum[0]*fl[13]-0.2886751345948129*nuVtSqSum[0]*fc[13]+0.2886751345948129*nuVtSqSum[1]*fl[10]-0.2886751345948129*nuVtSqSum[1]*fc[10]+0.25*nuVtSqSum[0]*fl[8]+0.25*nuVtSqSum[0]*fc[8]+0.25*nuVtSqSum[1]*fl[4]+0.25*nuVtSqSum[1]*fc[4]; 
  Gdiff2_l[6] = 0.2886751345948129*nuVtSqSum[1]*fl[15]-0.2886751345948129*nuVtSqSum[1]*fc[15]+0.2886751345948129*nuVtSqSum[0]*fl[14]-0.2886751345948129*nuVtSqSum[0]*fc[14]+0.25*nuVtSqSum[1]*fl[12]+0.25*nuVtSqSum[1]*fc[12]+0.25*nuVtSqSum[0]*fl[9]+0.25*nuVtSqSum[0]*fc[9]; 
  Gdiff2_l[7] = 0.2886751345948129*nuVtSqSum[0]*fl[15]-0.2886751345948129*nuVtSqSum[0]*fc[15]+0.2886751345948129*nuVtSqSum[1]*fl[14]-0.2886751345948129*nuVtSqSum[1]*fc[14]+0.25*nuVtSqSum[0]*fl[12]+0.25*nuVtSqSum[0]*fc[12]+0.25*nuVtSqSum[1]*fl[9]+0.25*nuVtSqSum[1]*fc[9]; 

  Gdiff2_r[0] = (-0.2886751345948129*nuVtSqSum[1]*fr[6])+0.2886751345948129*nuVtSqSum[1]*fc[6]-0.2886751345948129*nuVtSqSum[0]*fr[3]+0.2886751345948129*nuVtSqSum[0]*fc[3]+0.25*fr[1]*nuVtSqSum[1]+0.25*fc[1]*nuVtSqSum[1]+0.25*fr[0]*nuVtSqSum[0]+0.25*fc[0]*nuVtSqSum[0]; 
  Gdiff2_r[1] = (-0.2886751345948129*nuVtSqSum[0]*fr[6])+0.2886751345948129*nuVtSqSum[0]*fc[6]-0.2886751345948129*nuVtSqSum[1]*fr[3]+0.2886751345948129*nuVtSqSum[1]*fc[3]+0.25*fr[0]*nuVtSqSum[1]+0.25*fc[0]*nuVtSqSum[1]+0.25*nuVtSqSum[0]*fr[1]+0.25*nuVtSqSum[0]*fc[1]; 
  Gdiff2_r[2] = (-0.2886751345948129*nuVtSqSum[1]*fr[11])+0.2886751345948129*nuVtSqSum[1]*fc[11]-0.2886751345948129*nuVtSqSum[0]*fr[7]+0.2886751345948129*nuVtSqSum[0]*fc[7]+0.25*nuVtSqSum[1]*fr[5]+0.25*nuVtSqSum[1]*fc[5]+0.25*nuVtSqSum[0]*fr[2]+0.25*nuVtSqSum[0]*fc[2]; 
  Gdiff2_r[3] = (-0.2886751345948129*nuVtSqSum[1]*fr[13])+0.2886751345948129*nuVtSqSum[1]*fc[13]-0.2886751345948129*nuVtSqSum[0]*fr[10]+0.2886751345948129*nuVtSqSum[0]*fc[10]+0.25*nuVtSqSum[1]*fr[8]+0.25*nuVtSqSum[1]*fc[8]+0.25*nuVtSqSum[0]*fr[4]+0.25*nuVtSqSum[0]*fc[4]; 
  Gdiff2_r[4] = (-0.2886751345948129*nuVtSqSum[0]*fr[11])+0.2886751345948129*nuVtSqSum[0]*fc[11]-0.2886751345948129*nuVtSqSum[1]*fr[7]+0.2886751345948129*nuVtSqSum[1]*fc[7]+0.25*nuVtSqSum[0]*fr[5]+0.25*nuVtSqSum[0]*fc[5]+0.25*nuVtSqSum[1]*fr[2]+0.25*nuVtSqSum[1]*fc[2]; 
  Gdiff2_r[5] = (-0.2886751345948129*nuVtSqSum[0]*fr[13])+0.2886751345948129*nuVtSqSum[0]*fc[13]-0.2886751345948129*nuVtSqSum[1]*fr[10]+0.2886751345948129*nuVtSqSum[1]*fc[10]+0.25*nuVtSqSum[0]*fr[8]+0.25*nuVtSqSum[0]*fc[8]+0.25*nuVtSqSum[1]*fr[4]+0.25*nuVtSqSum[1]*fc[4]; 
  Gdiff2_r[6] = (-0.2886751345948129*nuVtSqSum[1]*fr[15])+0.2886751345948129*nuVtSqSum[1]*fc[15]-0.2886751345948129*nuVtSqSum[0]*fr[14]+0.2886751345948129*nuVtSqSum[0]*fc[14]+0.25*nuVtSqSum[1]*fr[12]+0.25*nuVtSqSum[1]*fc[12]+0.25*nuVtSqSum[0]*fr[9]+0.25*nuVtSqSum[0]*fc[9]; 
  Gdiff2_r[7] = (-0.2886751345948129*nuVtSqSum[0]*fr[15])+0.2886751345948129*nuVtSqSum[0]*fc[15]-0.2886751345948129*nuVtSqSum[1]*fr[14]+0.2886751345948129*nuVtSqSum[1]*fc[14]+0.25*nuVtSqSum[0]*fr[12]+0.25*nuVtSqSum[0]*fc[12]+0.25*nuVtSqSum[1]*fr[9]+0.25*nuVtSqSum[1]*fc[9]; 

  Gdiff_l[0] = (-0.5412658773652741*nuVtSqSum[1]*fl[6])-0.5412658773652741*nuVtSqSum[1]*fc[6]-0.5412658773652741*nuVtSqSum[0]*fl[3]-0.5412658773652741*nuVtSqSum[0]*fc[3]-0.5625*fl[1]*nuVtSqSum[1]+0.5625*fc[1]*nuVtSqSum[1]-0.5625*fl[0]*nuVtSqSum[0]+0.5625*fc[0]*nuVtSqSum[0]; 
  Gdiff_l[1] = (-0.5412658773652741*nuVtSqSum[0]*fl[6])-0.5412658773652741*nuVtSqSum[0]*fc[6]-0.5412658773652741*nuVtSqSum[1]*fl[3]-0.5412658773652741*nuVtSqSum[1]*fc[3]-0.5625*fl[0]*nuVtSqSum[1]+0.5625*fc[0]*nuVtSqSum[1]-0.5625*nuVtSqSum[0]*fl[1]+0.5625*nuVtSqSum[0]*fc[1]; 
  Gdiff_l[2] = (-0.5412658773652741*nuVtSqSum[1]*fl[11])-0.5412658773652741*nuVtSqSum[1]*fc[11]-0.5412658773652741*nuVtSqSum[0]*fl[7]-0.5412658773652741*nuVtSqSum[0]*fc[7]-0.5625*nuVtSqSum[1]*fl[5]+0.5625*nuVtSqSum[1]*fc[5]-0.5625*nuVtSqSum[0]*fl[2]+0.5625*nuVtSqSum[0]*fc[2]; 
  Gdiff_l[3] = (-0.5412658773652741*nuVtSqSum[1]*fl[13])-0.5412658773652741*nuVtSqSum[1]*fc[13]-0.5412658773652741*nuVtSqSum[0]*fl[10]-0.5412658773652741*nuVtSqSum[0]*fc[10]-0.5625*nuVtSqSum[1]*fl[8]+0.5625*nuVtSqSum[1]*fc[8]-0.5625*nuVtSqSum[0]*fl[4]+0.5625*nuVtSqSum[0]*fc[4]; 
  Gdiff_l[4] = (-0.5412658773652741*nuVtSqSum[0]*fl[11])-0.5412658773652741*nuVtSqSum[0]*fc[11]-0.5412658773652741*nuVtSqSum[1]*fl[7]-0.5412658773652741*nuVtSqSum[1]*fc[7]-0.5625*nuVtSqSum[0]*fl[5]+0.5625*nuVtSqSum[0]*fc[5]-0.5625*nuVtSqSum[1]*fl[2]+0.5625*nuVtSqSum[1]*fc[2]; 
  Gdiff_l[5] = (-0.5412658773652741*nuVtSqSum[0]*fl[13])-0.5412658773652741*nuVtSqSum[0]*fc[13]-0.5412658773652741*nuVtSqSum[1]*fl[10]-0.5412658773652741*nuVtSqSum[1]*fc[10]-0.5625*nuVtSqSum[0]*fl[8]+0.5625*nuVtSqSum[0]*fc[8]-0.5625*nuVtSqSum[1]*fl[4]+0.5625*nuVtSqSum[1]*fc[4]; 
  Gdiff_l[6] = (-0.5412658773652741*nuVtSqSum[1]*fl[15])-0.5412658773652741*nuVtSqSum[1]*fc[15]-0.5412658773652741*nuVtSqSum[0]*fl[14]-0.5412658773652741*nuVtSqSum[0]*fc[14]-0.5625*nuVtSqSum[1]*fl[12]+0.5625*nuVtSqSum[1]*fc[12]-0.5625*nuVtSqSum[0]*fl[9]+0.5625*nuVtSqSum[0]*fc[9]; 
  Gdiff_l[7] = (-0.5412658773652741*nuVtSqSum[0]*fl[15])-0.5412658773652741*nuVtSqSum[0]*fc[15]-0.5412658773652741*nuVtSqSum[1]*fl[14]-0.5412658773652741*nuVtSqSum[1]*fc[14]-0.5625*nuVtSqSum[0]*fl[12]+0.5625*nuVtSqSum[0]*fc[12]-0.5625*nuVtSqSum[1]*fl[9]+0.5625*nuVtSqSum[1]*fc[9]; 

  Gdiff_r[0] = (-0.5412658773652741*nuVtSqSum[1]*fr[6])-0.5412658773652741*nuVtSqSum[1]*fc[6]-0.5412658773652741*nuVtSqSum[0]*fr[3]-0.5412658773652741*nuVtSqSum[0]*fc[3]+0.5625*fr[1]*nuVtSqSum[1]-0.5625*fc[1]*nuVtSqSum[1]+0.5625*fr[0]*nuVtSqSum[0]-0.5625*fc[0]*nuVtSqSum[0]; 
  Gdiff_r[1] = (-0.5412658773652741*nuVtSqSum[0]*fr[6])-0.5412658773652741*nuVtSqSum[0]*fc[6]-0.5412658773652741*nuVtSqSum[1]*fr[3]-0.5412658773652741*nuVtSqSum[1]*fc[3]+0.5625*fr[0]*nuVtSqSum[1]-0.5625*fc[0]*nuVtSqSum[1]+0.5625*nuVtSqSum[0]*fr[1]-0.5625*nuVtSqSum[0]*fc[1]; 
  Gdiff_r[2] = (-0.5412658773652741*nuVtSqSum[1]*fr[11])-0.5412658773652741*nuVtSqSum[1]*fc[11]-0.5412658773652741*nuVtSqSum[0]*fr[7]-0.5412658773652741*nuVtSqSum[0]*fc[7]+0.5625*nuVtSqSum[1]*fr[5]-0.5625*nuVtSqSum[1]*fc[5]+0.5625*nuVtSqSum[0]*fr[2]-0.5625*nuVtSqSum[0]*fc[2]; 
  Gdiff_r[3] = (-0.5412658773652741*nuVtSqSum[1]*fr[13])-0.5412658773652741*nuVtSqSum[1]*fc[13]-0.5412658773652741*nuVtSqSum[0]*fr[10]-0.5412658773652741*nuVtSqSum[0]*fc[10]+0.5625*nuVtSqSum[1]*fr[8]-0.5625*nuVtSqSum[1]*fc[8]+0.5625*nuVtSqSum[0]*fr[4]-0.5625*nuVtSqSum[0]*fc[4]; 
  Gdiff_r[4] = (-0.5412658773652741*nuVtSqSum[0]*fr[11])-0.5412658773652741*nuVtSqSum[0]*fc[11]-0.5412658773652741*nuVtSqSum[1]*fr[7]-0.5412658773652741*nuVtSqSum[1]*fc[7]+0.5625*nuVtSqSum[0]*fr[5]-0.5625*nuVtSqSum[0]*fc[5]+0.5625*nuVtSqSum[1]*fr[2]-0.5625*nuVtSqSum[1]*fc[2]; 
  Gdiff_r[5] = (-0.5412658773652741*nuVtSqSum[0]*fr[13])-0.5412658773652741*nuVtSqSum[0]*fc[13]-0.5412658773652741*nuVtSqSum[1]*fr[10]-0.5412658773652741*nuVtSqSum[1]*fc[10]+0.5625*nuVtSqSum[0]*fr[8]-0.5625*nuVtSqSum[0]*fc[8]+0.5625*nuVtSqSum[1]*fr[4]-0.5625*nuVtSqSum[1]*fc[4]; 
  Gdiff_r[6] = (-0.5412658773652741*nuVtSqSum[1]*fr[15])-0.5412658773652741*nuVtSqSum[1]*fc[15]-0.5412658773652741*nuVtSqSum[0]*fr[14]-0.5412658773652741*nuVtSqSum[0]*fc[14]+0.5625*nuVtSqSum[1]*fr[12]-0.5625*nuVtSqSum[1]*fc[12]+0.5625*nuVtSqSum[0]*fr[9]-0.5625*nuVtSqSum[0]*fc[9]; 
  Gdiff_r[7] = (-0.5412658773652741*nuVtSqSum[0]*fr[15])-0.5412658773652741*nuVtSqSum[0]*fc[15]-0.5412658773652741*nuVtSqSum[1]*fr[14]-0.5412658773652741*nuVtSqSum[1]*fc[14]+0.5625*nuVtSqSum[0]*fr[12]-0.5625*nuVtSqSum[0]*fc[12]+0.5625*nuVtSqSum[1]*fr[9]-0.5625*nuVtSqSum[1]*fc[9]; 

  Ghat_l[0] = Gdiff_l[0]*rdv2+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[1]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[0]; 
  Ghat_l[1] = Gdiff_l[1]*rdv2+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[1]+0.3535533905932737*fUpwind_l[0]*alphaDrSurf_l[1]; 
  Ghat_l[2] = Gdiff_l[2]*rdv2+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[4]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[2]; 
  Ghat_l[3] = Gdiff_l[3]*rdv2+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[5]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[3]; 
  Ghat_l[4] = Gdiff_l[4]*rdv2+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[4]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[2]; 
  Ghat_l[5] = Gdiff_l[5]*rdv2+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[5]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[3]; 
  Ghat_l[6] = Gdiff_l[6]*rdv2+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[7]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[6]; 
  Ghat_l[7] = Gdiff_l[7]*rdv2+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[7]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[6]; 

  Ghat_r[0] = (-1.0*Gdiff_r[0]*rdv2)+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[1]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[0]; 
  Ghat_r[1] = (-1.0*Gdiff_r[1]*rdv2)+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[1]+0.3535533905932737*fUpwind_r[0]*alphaDrSurf_r[1]; 
  Ghat_r[2] = (-1.0*Gdiff_r[2]*rdv2)+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[4]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[2]; 
  Ghat_r[3] = (-1.0*Gdiff_r[3]*rdv2)+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[5]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[3]; 
  Ghat_r[4] = (-1.0*Gdiff_r[4]*rdv2)+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[4]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[2]; 
  Ghat_r[5] = (-1.0*Gdiff_r[5]*rdv2)+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[5]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[3]; 
  Ghat_r[6] = (-1.0*Gdiff_r[6]*rdv2)+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[7]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[6]; 
  Ghat_r[7] = (-1.0*Gdiff_r[7]*rdv2)+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[7]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[6]; 

  out[0] += 0.7071067811865475*Ghat_l[0]*rdv2-0.7071067811865475*Ghat_r[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat_l[1]*rdv2-0.7071067811865475*Ghat_r[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat_l[2]*rdv2-0.7071067811865475*Ghat_r[2]*rdv2; 
  out[3] += 1.224744871391589*Gdiff2_r[0]*rdvSq4-1.224744871391589*Gdiff2_l[0]*rdvSq4-1.224744871391589*Ghat_r[0]*rdv2-1.224744871391589*Ghat_l[0]*rdv2; 
  out[4] += 0.7071067811865475*Ghat_l[3]*rdv2-0.7071067811865475*Ghat_r[3]*rdv2; 
  out[5] += 0.7071067811865475*Ghat_l[4]*rdv2-0.7071067811865475*Ghat_r[4]*rdv2; 
  out[6] += 1.224744871391589*Gdiff2_r[1]*rdvSq4-1.224744871391589*Gdiff2_l[1]*rdvSq4-1.224744871391589*Ghat_r[1]*rdv2-1.224744871391589*Ghat_l[1]*rdv2; 
  out[7] += 1.224744871391589*Gdiff2_r[2]*rdvSq4-1.224744871391589*Gdiff2_l[2]*rdvSq4-1.224744871391589*Ghat_r[2]*rdv2-1.224744871391589*Ghat_l[2]*rdv2; 
  out[8] += 0.7071067811865475*Ghat_l[5]*rdv2-0.7071067811865475*Ghat_r[5]*rdv2; 
  out[9] += 0.7071067811865475*Ghat_l[6]*rdv2-0.7071067811865475*Ghat_r[6]*rdv2; 
  out[10] += 1.224744871391589*Gdiff2_r[3]*rdvSq4-1.224744871391589*Gdiff2_l[3]*rdvSq4-1.224744871391589*Ghat_r[3]*rdv2-1.224744871391589*Ghat_l[3]*rdv2; 
  out[11] += 1.224744871391589*Gdiff2_r[4]*rdvSq4-1.224744871391589*Gdiff2_l[4]*rdvSq4-1.224744871391589*Ghat_r[4]*rdv2-1.224744871391589*Ghat_l[4]*rdv2; 
  out[12] += 0.7071067811865475*Ghat_l[7]*rdv2-0.7071067811865475*Ghat_r[7]*rdv2; 
  out[13] += 1.224744871391589*Gdiff2_r[5]*rdvSq4-1.224744871391589*Gdiff2_l[5]*rdvSq4-1.224744871391589*Ghat_r[5]*rdv2-1.224744871391589*Ghat_l[5]*rdv2; 
  out[14] += 1.224744871391589*Gdiff2_r[6]*rdvSq4-1.224744871391589*Gdiff2_l[6]*rdvSq4-1.224744871391589*Ghat_r[6]*rdv2-1.224744871391589*Ghat_l[6]*rdv2; 
  out[15] += 1.224744871391589*Gdiff2_r[7]*rdvSq4-1.224744871391589*Gdiff2_l[7]*rdvSq4-1.224744871391589*Ghat_r[7]*rdv2-1.224744871391589*Ghat_l[7]*rdv2; 
} 