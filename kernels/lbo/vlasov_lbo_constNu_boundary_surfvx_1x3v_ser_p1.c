#include <gkyl_lbo_kernels.h> 
GKYL_CU_DH void vlasov_lbo_constNu_boundary_surfvx_1x3v_ser_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSum, const double *nuVtSqSum_l, const double *nuVtSqSum_r, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[4]:          Cell-center coordinates. 
  // dxv[4]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[6]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]:  sum of thermal speeds squared time their respective collisionalities. 
  // fSkin/Edge:      Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[1]; 
  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 

  double Gdiff2_l[8]; 
  double Gdiff2_r[8]; 


  Gdiff2_l[0] = -0.08333333333333333*(3.464101615137754*nuVtSqSum_l[1]*fSkin[5]-3.464101615137754*nuVtSqSum_l[1]*fEdge[5]+3.464101615137754*nuVtSqSum_l[0]*fSkin[2]-3.464101615137754*nuVtSqSum_l[0]*fEdge[2]+((-3.0*fSkin[1])-3.0*fEdge[1])*nuVtSqSum_l[1]+((-3.0*fSkin[0])-3.0*fEdge[0])*nuVtSqSum_l[0]); 
  Gdiff2_l[1] = -0.08333333333333333*(3.464101615137754*nuVtSqSum_l[0]*fSkin[5]-3.464101615137754*nuVtSqSum_l[0]*fEdge[5]+3.464101615137754*nuVtSqSum_l[1]*fSkin[2]-3.464101615137754*nuVtSqSum_l[1]*fEdge[2]+((-3.0*fSkin[0])-3.0*fEdge[0])*nuVtSqSum_l[1]-3.0*nuVtSqSum_l[0]*fSkin[1]-3.0*nuVtSqSum_l[0]*fEdge[1]); 
  Gdiff2_l[2] = -0.08333333333333333*(3.464101615137754*nuVtSqSum_l[1]*fSkin[11]-3.464101615137754*nuVtSqSum_l[1]*fEdge[11]+3.464101615137754*nuVtSqSum_l[0]*fSkin[7]-3.464101615137754*nuVtSqSum_l[0]*fEdge[7]-3.0*nuVtSqSum_l[1]*fSkin[6]-3.0*nuVtSqSum_l[1]*fEdge[6]-3.0*nuVtSqSum_l[0]*fSkin[3]-3.0*nuVtSqSum_l[0]*fEdge[3]); 
  Gdiff2_l[3] = -0.08333333333333333*(3.464101615137754*nuVtSqSum_l[1]*fSkin[12]-3.464101615137754*nuVtSqSum_l[1]*fEdge[12]+3.464101615137754*nuVtSqSum_l[0]*fSkin[9]-3.464101615137754*nuVtSqSum_l[0]*fEdge[9]-3.0*nuVtSqSum_l[1]*fSkin[8]-3.0*nuVtSqSum_l[1]*fEdge[8]-3.0*nuVtSqSum_l[0]*fSkin[4]-3.0*nuVtSqSum_l[0]*fEdge[4]); 
  Gdiff2_l[4] = -0.08333333333333333*(3.464101615137754*nuVtSqSum_l[0]*fSkin[11]-3.464101615137754*nuVtSqSum_l[0]*fEdge[11]+3.464101615137754*nuVtSqSum_l[1]*fSkin[7]-3.464101615137754*nuVtSqSum_l[1]*fEdge[7]-3.0*nuVtSqSum_l[0]*fSkin[6]-3.0*nuVtSqSum_l[0]*fEdge[6]-3.0*nuVtSqSum_l[1]*fSkin[3]-3.0*nuVtSqSum_l[1]*fEdge[3]); 
  Gdiff2_l[5] = -0.08333333333333333*(3.464101615137754*nuVtSqSum_l[0]*fSkin[12]-3.464101615137754*nuVtSqSum_l[0]*fEdge[12]+3.464101615137754*nuVtSqSum_l[1]*fSkin[9]-3.464101615137754*nuVtSqSum_l[1]*fEdge[9]-3.0*nuVtSqSum_l[0]*fSkin[8]-3.0*nuVtSqSum_l[0]*fEdge[8]-3.0*nuVtSqSum_l[1]*fSkin[4]-3.0*nuVtSqSum_l[1]*fEdge[4]); 
  Gdiff2_l[6] = -0.08333333333333333*(3.464101615137754*nuVtSqSum_l[1]*fSkin[15]-3.464101615137754*nuVtSqSum_l[1]*fEdge[15]+3.464101615137754*nuVtSqSum_l[0]*fSkin[14]-3.464101615137754*nuVtSqSum_l[0]*fEdge[14]-3.0*nuVtSqSum_l[1]*fSkin[13]-3.0*nuVtSqSum_l[1]*fEdge[13]-3.0*nuVtSqSum_l[0]*fSkin[10]-3.0*nuVtSqSum_l[0]*fEdge[10]); 
  Gdiff2_l[7] = -0.08333333333333333*(3.464101615137754*nuVtSqSum_l[0]*fSkin[15]-3.464101615137754*nuVtSqSum_l[0]*fEdge[15]+3.464101615137754*nuVtSqSum_l[1]*fSkin[14]-3.464101615137754*nuVtSqSum_l[1]*fEdge[14]-3.0*nuVtSqSum_l[0]*fSkin[13]-3.0*nuVtSqSum_l[0]*fEdge[13]-3.0*nuVtSqSum_l[1]*fSkin[10]-3.0*nuVtSqSum_l[1]*fEdge[10]); 


  Gdiff2_r[0] = 0.08333333333333333*(3.464101615137754*nuVtSqSum_r[1]*fSkin[5]-3.464101615137754*nuVtSqSum_r[1]*fEdge[5]+3.464101615137754*nuVtSqSum_r[0]*fSkin[2]-3.464101615137754*nuVtSqSum_r[0]*fEdge[2]+(3.0*fSkin[1]+3.0*fEdge[1])*nuVtSqSum_r[1]+(3.0*fSkin[0]+3.0*fEdge[0])*nuVtSqSum_r[0]); 
  Gdiff2_r[1] = 0.08333333333333333*(3.464101615137754*nuVtSqSum_r[0]*fSkin[5]-3.464101615137754*nuVtSqSum_r[0]*fEdge[5]+3.464101615137754*nuVtSqSum_r[1]*fSkin[2]-3.464101615137754*nuVtSqSum_r[1]*fEdge[2]+(3.0*fSkin[0]+3.0*fEdge[0])*nuVtSqSum_r[1]+3.0*nuVtSqSum_r[0]*fSkin[1]+3.0*nuVtSqSum_r[0]*fEdge[1]); 
  Gdiff2_r[2] = 0.08333333333333333*(3.464101615137754*nuVtSqSum_r[1]*fSkin[11]-3.464101615137754*nuVtSqSum_r[1]*fEdge[11]+3.464101615137754*nuVtSqSum_r[0]*fSkin[7]-3.464101615137754*nuVtSqSum_r[0]*fEdge[7]+3.0*nuVtSqSum_r[1]*fSkin[6]+3.0*nuVtSqSum_r[1]*fEdge[6]+3.0*nuVtSqSum_r[0]*fSkin[3]+3.0*nuVtSqSum_r[0]*fEdge[3]); 
  Gdiff2_r[3] = 0.08333333333333333*(3.464101615137754*nuVtSqSum_r[1]*fSkin[12]-3.464101615137754*nuVtSqSum_r[1]*fEdge[12]+3.464101615137754*nuVtSqSum_r[0]*fSkin[9]-3.464101615137754*nuVtSqSum_r[0]*fEdge[9]+3.0*nuVtSqSum_r[1]*fSkin[8]+3.0*nuVtSqSum_r[1]*fEdge[8]+3.0*nuVtSqSum_r[0]*fSkin[4]+3.0*nuVtSqSum_r[0]*fEdge[4]); 
  Gdiff2_r[4] = 0.08333333333333333*(3.464101615137754*nuVtSqSum_r[0]*fSkin[11]-3.464101615137754*nuVtSqSum_r[0]*fEdge[11]+3.464101615137754*nuVtSqSum_r[1]*fSkin[7]-3.464101615137754*nuVtSqSum_r[1]*fEdge[7]+3.0*nuVtSqSum_r[0]*fSkin[6]+3.0*nuVtSqSum_r[0]*fEdge[6]+3.0*nuVtSqSum_r[1]*fSkin[3]+3.0*nuVtSqSum_r[1]*fEdge[3]); 
  Gdiff2_r[5] = 0.08333333333333333*(3.464101615137754*nuVtSqSum_r[0]*fSkin[12]-3.464101615137754*nuVtSqSum_r[0]*fEdge[12]+3.464101615137754*nuVtSqSum_r[1]*fSkin[9]-3.464101615137754*nuVtSqSum_r[1]*fEdge[9]+3.0*nuVtSqSum_r[0]*fSkin[8]+3.0*nuVtSqSum_r[0]*fEdge[8]+3.0*nuVtSqSum_r[1]*fSkin[4]+3.0*nuVtSqSum_r[1]*fEdge[4]); 
  Gdiff2_r[6] = 0.08333333333333333*(3.464101615137754*nuVtSqSum_r[1]*fSkin[15]-3.464101615137754*nuVtSqSum_r[1]*fEdge[15]+3.464101615137754*nuVtSqSum_r[0]*fSkin[14]-3.464101615137754*nuVtSqSum_r[0]*fEdge[14]+3.0*nuVtSqSum_r[1]*fSkin[13]+3.0*nuVtSqSum_r[1]*fEdge[13]+3.0*nuVtSqSum_r[0]*fSkin[10]+3.0*nuVtSqSum_r[0]*fEdge[10]); 
  Gdiff2_r[7] = 0.08333333333333333*(3.464101615137754*nuVtSqSum_r[0]*fSkin[15]-3.464101615137754*nuVtSqSum_r[0]*fEdge[15]+3.464101615137754*nuVtSqSum_r[1]*fSkin[14]-3.464101615137754*nuVtSqSum_r[1]*fEdge[14]+3.0*nuVtSqSum_r[0]*fSkin[13]+3.0*nuVtSqSum_r[0]*fEdge[13]+3.0*nuVtSqSum_r[1]*fSkin[10]+3.0*nuVtSqSum_r[1]*fEdge[10]); 

  if (edge == -1) { 

  const double *sumNuUx_r = &nuUSum[0]; 

  double alphaDrSurf_r[8]; 
  alphaDrSurf_r[0] = ((-2.82842712474619*w[1])-1.414213562373095*dxv[1])*nuSum+2.0*sumNuUx_r[0]; 
  alphaDrSurf_r[1] = 2.0*sumNuUx_r[1]; 

  double fUpwindQuad_r[8];
  if (alphaDrSurf_r[0]-alphaDrSurf_r[1] > 0) { 
    fUpwindQuad_r[0] = (-0.4330127018922193*fSkin[15])+0.4330127018922193*fSkin[14]-0.25*fSkin[13]+0.4330127018922193*(fSkin[12]+fSkin[11])+0.25*fSkin[10]-0.4330127018922193*fSkin[9]+0.25*fSkin[8]-0.4330127018922193*fSkin[7]+0.25*fSkin[6]-0.4330127018922193*fSkin[5]-0.25*(fSkin[4]+fSkin[3])+0.4330127018922193*fSkin[2]-0.25*fSkin[1]+0.25*fSkin[0]; 
  } else { 
    fUpwindQuad_r[0] = 0.4330127018922193*fEdge[15]-0.4330127018922193*fEdge[14]-0.25*fEdge[13]-0.4330127018922193*(fEdge[12]+fEdge[11])+0.25*fEdge[10]+0.4330127018922193*fEdge[9]+0.25*fEdge[8]+0.4330127018922193*fEdge[7]+0.25*fEdge[6]+0.4330127018922193*fEdge[5]-0.25*(fEdge[4]+fEdge[3])-0.4330127018922193*fEdge[2]-0.25*fEdge[1]+0.25*fEdge[0]; 
  } 
  if (alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[1] = 0.4330127018922193*(fSkin[15]+fSkin[14])+0.25*fSkin[13]-0.4330127018922193*(fSkin[12]+fSkin[11])+0.25*fSkin[10]-0.4330127018922193*fSkin[9]-0.25*fSkin[8]-0.4330127018922193*fSkin[7]-0.25*fSkin[6]+0.4330127018922193*fSkin[5]-0.25*(fSkin[4]+fSkin[3])+0.4330127018922193*fSkin[2]+0.25*(fSkin[1]+fSkin[0]); 
  } else { 
    fUpwindQuad_r[1] = (-0.4330127018922193*(fEdge[15]+fEdge[14]))+0.25*fEdge[13]+0.4330127018922193*(fEdge[12]+fEdge[11])+0.25*fEdge[10]+0.4330127018922193*fEdge[9]-0.25*fEdge[8]+0.4330127018922193*fEdge[7]-0.25*fEdge[6]-0.4330127018922193*fEdge[5]-0.25*(fEdge[4]+fEdge[3])-0.4330127018922193*fEdge[2]+0.25*(fEdge[1]+fEdge[0]); 
  } 
  if (alphaDrSurf_r[0]-alphaDrSurf_r[1] > 0) { 
    fUpwindQuad_r[2] = 0.4330127018922193*fSkin[15]-0.4330127018922193*fSkin[14]+0.25*fSkin[13]+0.4330127018922193*fSkin[12]-0.4330127018922193*fSkin[11]-0.25*fSkin[10]-0.4330127018922193*fSkin[9]+0.25*fSkin[8]+0.4330127018922193*fSkin[7]-0.25*fSkin[6]-0.4330127018922193*fSkin[5]-0.25*fSkin[4]+0.25*fSkin[3]+0.4330127018922193*fSkin[2]-0.25*fSkin[1]+0.25*fSkin[0]; 
  } else { 
    fUpwindQuad_r[2] = (-0.4330127018922193*fEdge[15])+0.4330127018922193*fEdge[14]+0.25*fEdge[13]-0.4330127018922193*fEdge[12]+0.4330127018922193*fEdge[11]-0.25*fEdge[10]+0.4330127018922193*fEdge[9]+0.25*fEdge[8]-0.4330127018922193*fEdge[7]-0.25*fEdge[6]+0.4330127018922193*fEdge[5]-0.25*fEdge[4]+0.25*fEdge[3]-0.4330127018922193*fEdge[2]-0.25*fEdge[1]+0.25*fEdge[0]; 
  } 
  if (alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[3] = (-0.4330127018922193*(fSkin[15]+fSkin[14]))-0.25*fSkin[13]-0.4330127018922193*fSkin[12]+0.4330127018922193*fSkin[11]-0.25*fSkin[10]-0.4330127018922193*fSkin[9]-0.25*fSkin[8]+0.4330127018922193*fSkin[7]+0.25*fSkin[6]+0.4330127018922193*fSkin[5]-0.25*fSkin[4]+0.25*fSkin[3]+0.4330127018922193*fSkin[2]+0.25*(fSkin[1]+fSkin[0]); 
  } else { 
    fUpwindQuad_r[3] = 0.4330127018922193*(fEdge[15]+fEdge[14])-0.25*fEdge[13]+0.4330127018922193*fEdge[12]-0.4330127018922193*fEdge[11]-0.25*fEdge[10]+0.4330127018922193*fEdge[9]-0.25*fEdge[8]-0.4330127018922193*fEdge[7]+0.25*fEdge[6]-0.4330127018922193*fEdge[5]-0.25*fEdge[4]+0.25*fEdge[3]-0.4330127018922193*fEdge[2]+0.25*(fEdge[1]+fEdge[0]); 
  } 
  if (alphaDrSurf_r[0]-alphaDrSurf_r[1] > 0) { 
    fUpwindQuad_r[4] = 0.4330127018922193*fSkin[15]-0.4330127018922193*fSkin[14]+0.25*fSkin[13]-0.4330127018922193*fSkin[12]+0.4330127018922193*fSkin[11]-0.25*fSkin[10]+0.4330127018922193*fSkin[9]-0.25*fSkin[8]-0.4330127018922193*fSkin[7]+0.25*fSkin[6]-0.4330127018922193*fSkin[5]+0.25*fSkin[4]-0.25*fSkin[3]+0.4330127018922193*fSkin[2]-0.25*fSkin[1]+0.25*fSkin[0]; 
  } else { 
    fUpwindQuad_r[4] = (-0.4330127018922193*fEdge[15])+0.4330127018922193*fEdge[14]+0.25*fEdge[13]+0.4330127018922193*fEdge[12]-0.4330127018922193*fEdge[11]-0.25*fEdge[10]-0.4330127018922193*fEdge[9]-0.25*fEdge[8]+0.4330127018922193*fEdge[7]+0.25*fEdge[6]+0.4330127018922193*fEdge[5]+0.25*fEdge[4]-0.25*fEdge[3]-0.4330127018922193*fEdge[2]-0.25*fEdge[1]+0.25*fEdge[0]; 
  } 
  if (alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[5] = (-0.4330127018922193*(fSkin[15]+fSkin[14]))-0.25*fSkin[13]+0.4330127018922193*fSkin[12]-0.4330127018922193*fSkin[11]-0.25*fSkin[10]+0.4330127018922193*fSkin[9]+0.25*fSkin[8]-0.4330127018922193*fSkin[7]-0.25*fSkin[6]+0.4330127018922193*fSkin[5]+0.25*fSkin[4]-0.25*fSkin[3]+0.4330127018922193*fSkin[2]+0.25*(fSkin[1]+fSkin[0]); 
  } else { 
    fUpwindQuad_r[5] = 0.4330127018922193*(fEdge[15]+fEdge[14])-0.25*fEdge[13]-0.4330127018922193*fEdge[12]+0.4330127018922193*fEdge[11]-0.25*fEdge[10]-0.4330127018922193*fEdge[9]+0.25*fEdge[8]+0.4330127018922193*fEdge[7]-0.25*fEdge[6]-0.4330127018922193*fEdge[5]+0.25*fEdge[4]-0.25*fEdge[3]-0.4330127018922193*fEdge[2]+0.25*(fEdge[1]+fEdge[0]); 
  } 
  if (alphaDrSurf_r[0]-alphaDrSurf_r[1] > 0) { 
    fUpwindQuad_r[6] = (-0.4330127018922193*fSkin[15])+0.4330127018922193*fSkin[14]-0.25*fSkin[13]-0.4330127018922193*(fSkin[12]+fSkin[11])+0.25*fSkin[10]+0.4330127018922193*fSkin[9]-0.25*fSkin[8]+0.4330127018922193*fSkin[7]-0.25*fSkin[6]-0.4330127018922193*fSkin[5]+0.25*(fSkin[4]+fSkin[3])+0.4330127018922193*fSkin[2]-0.25*fSkin[1]+0.25*fSkin[0]; 
  } else { 
    fUpwindQuad_r[6] = 0.4330127018922193*fEdge[15]-0.4330127018922193*fEdge[14]-0.25*fEdge[13]+0.4330127018922193*(fEdge[12]+fEdge[11])+0.25*fEdge[10]-0.4330127018922193*fEdge[9]-0.25*fEdge[8]-0.4330127018922193*fEdge[7]-0.25*fEdge[6]+0.4330127018922193*fEdge[5]+0.25*(fEdge[4]+fEdge[3])-0.4330127018922193*fEdge[2]-0.25*fEdge[1]+0.25*fEdge[0]; 
  } 
  if (alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[7] = 0.4330127018922193*(fSkin[15]+fSkin[14])+0.25*fSkin[13]+0.4330127018922193*(fSkin[12]+fSkin[11])+0.25*fSkin[10]+0.4330127018922193*fSkin[9]+0.25*fSkin[8]+0.4330127018922193*fSkin[7]+0.25*fSkin[6]+0.4330127018922193*fSkin[5]+0.25*(fSkin[4]+fSkin[3])+0.4330127018922193*fSkin[2]+0.25*(fSkin[1]+fSkin[0]); 
  } else { 
    fUpwindQuad_r[7] = (-0.4330127018922193*(fEdge[15]+fEdge[14]))+0.25*fEdge[13]-0.4330127018922193*(fEdge[12]+fEdge[11])+0.25*fEdge[10]-0.4330127018922193*fEdge[9]+0.25*fEdge[8]-0.4330127018922193*fEdge[7]+0.25*fEdge[6]-0.4330127018922193*fEdge[5]+0.25*(fEdge[4]+fEdge[3])-0.4330127018922193*fEdge[2]+0.25*(fEdge[1]+fEdge[0]); 
  } 

  double fUpwind_r[8];
  fUpwind_r[0] = 0.3535533905932737*(fUpwindQuad_r[7]+fUpwindQuad_r[6]+fUpwindQuad_r[5]+fUpwindQuad_r[4]+fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[1] = 0.3535533905932737*(fUpwindQuad_r[7]-1.0*fUpwindQuad_r[6]+fUpwindQuad_r[5]-1.0*fUpwindQuad_r[4]+fUpwindQuad_r[3]-1.0*fUpwindQuad_r[2]+fUpwindQuad_r[1]-1.0*fUpwindQuad_r[0]); 
  fUpwind_r[2] = 0.3535533905932737*(fUpwindQuad_r[7]+fUpwindQuad_r[6]-1.0*(fUpwindQuad_r[5]+fUpwindQuad_r[4])+fUpwindQuad_r[3]+fUpwindQuad_r[2]-1.0*(fUpwindQuad_r[1]+fUpwindQuad_r[0])); 
  fUpwind_r[3] = 0.3535533905932737*(fUpwindQuad_r[7]+fUpwindQuad_r[6]+fUpwindQuad_r[5]+fUpwindQuad_r[4]-1.0*(fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]+fUpwindQuad_r[0])); 
  fUpwind_r[4] = 0.3535533905932737*(fUpwindQuad_r[7]-1.0*(fUpwindQuad_r[6]+fUpwindQuad_r[5])+fUpwindQuad_r[4]+fUpwindQuad_r[3]-1.0*(fUpwindQuad_r[2]+fUpwindQuad_r[1])+fUpwindQuad_r[0]); 
  fUpwind_r[5] = 0.3535533905932737*(fUpwindQuad_r[7]-1.0*fUpwindQuad_r[6]+fUpwindQuad_r[5]-1.0*(fUpwindQuad_r[4]+fUpwindQuad_r[3])+fUpwindQuad_r[2]-1.0*fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[6] = 0.3535533905932737*(fUpwindQuad_r[7]+fUpwindQuad_r[6]-1.0*(fUpwindQuad_r[5]+fUpwindQuad_r[4]+fUpwindQuad_r[3]+fUpwindQuad_r[2])+fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[7] = 0.3535533905932737*(fUpwindQuad_r[7]-1.0*(fUpwindQuad_r[6]+fUpwindQuad_r[5])+fUpwindQuad_r[4]-1.0*fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]-1.0*fUpwindQuad_r[0]); 

  double Gdiff_r[8]; 
  double Ghat_r[8]; 


  Gdiff_r[0] = -0.0625*(8.660254037844386*nuVtSqSum_r[1]*fSkin[5]+8.660254037844386*nuVtSqSum_r[1]*fEdge[5]+8.660254037844386*nuVtSqSum_r[0]*fSkin[2]+8.660254037844386*nuVtSqSum_r[0]*fEdge[2]+(9.0*fSkin[1]-9.0*fEdge[1])*nuVtSqSum_r[1]+(9.0*fSkin[0]-9.0*fEdge[0])*nuVtSqSum_r[0]); 
  Gdiff_r[1] = -0.0625*(8.660254037844386*nuVtSqSum_r[0]*fSkin[5]+8.660254037844386*nuVtSqSum_r[0]*fEdge[5]+8.660254037844386*nuVtSqSum_r[1]*fSkin[2]+8.660254037844386*nuVtSqSum_r[1]*fEdge[2]+(9.0*fSkin[0]-9.0*fEdge[0])*nuVtSqSum_r[1]+9.0*nuVtSqSum_r[0]*fSkin[1]-9.0*nuVtSqSum_r[0]*fEdge[1]); 
  Gdiff_r[2] = -0.0625*(8.660254037844386*nuVtSqSum_r[1]*fSkin[11]+8.660254037844386*nuVtSqSum_r[1]*fEdge[11]+8.660254037844386*nuVtSqSum_r[0]*fSkin[7]+8.660254037844386*nuVtSqSum_r[0]*fEdge[7]+9.0*nuVtSqSum_r[1]*fSkin[6]-9.0*nuVtSqSum_r[1]*fEdge[6]+9.0*nuVtSqSum_r[0]*fSkin[3]-9.0*nuVtSqSum_r[0]*fEdge[3]); 
  Gdiff_r[3] = -0.0625*(8.660254037844386*nuVtSqSum_r[1]*fSkin[12]+8.660254037844386*nuVtSqSum_r[1]*fEdge[12]+8.660254037844386*nuVtSqSum_r[0]*fSkin[9]+8.660254037844386*nuVtSqSum_r[0]*fEdge[9]+9.0*nuVtSqSum_r[1]*fSkin[8]-9.0*nuVtSqSum_r[1]*fEdge[8]+9.0*nuVtSqSum_r[0]*fSkin[4]-9.0*nuVtSqSum_r[0]*fEdge[4]); 
  Gdiff_r[4] = -0.0625*(8.660254037844386*nuVtSqSum_r[0]*fSkin[11]+8.660254037844386*nuVtSqSum_r[0]*fEdge[11]+8.660254037844386*nuVtSqSum_r[1]*fSkin[7]+8.660254037844386*nuVtSqSum_r[1]*fEdge[7]+9.0*nuVtSqSum_r[0]*fSkin[6]-9.0*nuVtSqSum_r[0]*fEdge[6]+9.0*nuVtSqSum_r[1]*fSkin[3]-9.0*nuVtSqSum_r[1]*fEdge[3]); 
  Gdiff_r[5] = -0.0625*(8.660254037844386*nuVtSqSum_r[0]*fSkin[12]+8.660254037844386*nuVtSqSum_r[0]*fEdge[12]+8.660254037844386*nuVtSqSum_r[1]*fSkin[9]+8.660254037844386*nuVtSqSum_r[1]*fEdge[9]+9.0*nuVtSqSum_r[0]*fSkin[8]-9.0*nuVtSqSum_r[0]*fEdge[8]+9.0*nuVtSqSum_r[1]*fSkin[4]-9.0*nuVtSqSum_r[1]*fEdge[4]); 
  Gdiff_r[6] = -0.0625*(8.660254037844386*nuVtSqSum_r[1]*fSkin[15]+8.660254037844386*nuVtSqSum_r[1]*fEdge[15]+8.660254037844386*nuVtSqSum_r[0]*fSkin[14]+8.660254037844386*nuVtSqSum_r[0]*fEdge[14]+9.0*nuVtSqSum_r[1]*fSkin[13]-9.0*nuVtSqSum_r[1]*fEdge[13]+9.0*nuVtSqSum_r[0]*fSkin[10]-9.0*nuVtSqSum_r[0]*fEdge[10]); 
  Gdiff_r[7] = -0.0625*(8.660254037844386*nuVtSqSum_r[0]*fSkin[15]+8.660254037844386*nuVtSqSum_r[0]*fEdge[15]+8.660254037844386*nuVtSqSum_r[1]*fSkin[14]+8.660254037844386*nuVtSqSum_r[1]*fEdge[14]+9.0*nuVtSqSum_r[0]*fSkin[13]-9.0*nuVtSqSum_r[0]*fEdge[13]+9.0*nuVtSqSum_r[1]*fSkin[10]-9.0*nuVtSqSum_r[1]*fEdge[10]); 

  Ghat_r[0] = (-1.0*Gdiff_r[0]*rdv2)+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[1]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[0]; 
  Ghat_r[1] = (-1.0*Gdiff_r[1]*rdv2)+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[1]+0.3535533905932737*fUpwind_r[0]*alphaDrSurf_r[1]; 
  Ghat_r[2] = (-1.0*Gdiff_r[2]*rdv2)+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[4]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[2]; 
  Ghat_r[3] = (-1.0*Gdiff_r[3]*rdv2)+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[5]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[3]; 
  Ghat_r[4] = (-1.0*Gdiff_r[4]*rdv2)+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[4]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[2]; 
  Ghat_r[5] = (-1.0*Gdiff_r[5]*rdv2)+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[5]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[3]; 
  Ghat_r[6] = (-1.0*Gdiff_r[6]*rdv2)+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[7]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[6]; 
  Ghat_r[7] = (-1.0*Gdiff_r[7]*rdv2)+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[7]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[6]; 

  out[0] += -0.7071067811865475*Ghat_r[0]*rdv2; 
  out[1] += -0.7071067811865475*Ghat_r[1]*rdv2; 
  out[2] += 1.060660171779821*nuVtSqSum_l[1]*fSkin[5]*rdvSq4+1.060660171779821*nuVtSqSum_l[0]*fSkin[2]*rdvSq4-0.6123724356957944*fSkin[1]*nuVtSqSum_l[1]*rdvSq4-0.6123724356957944*fSkin[0]*nuVtSqSum_l[0]*rdvSq4+1.224744871391589*Gdiff2_r[0]*rdvSq4-1.224744871391589*Ghat_r[0]*rdv2; 
  out[3] += -0.7071067811865475*Ghat_r[2]*rdv2; 
  out[4] += -0.7071067811865475*Ghat_r[3]*rdv2; 
  out[5] += 1.060660171779821*nuVtSqSum_l[0]*fSkin[5]*rdvSq4+1.060660171779821*nuVtSqSum_l[1]*fSkin[2]*rdvSq4-0.6123724356957944*fSkin[0]*nuVtSqSum_l[1]*rdvSq4-0.6123724356957944*nuVtSqSum_l[0]*fSkin[1]*rdvSq4+1.224744871391589*Gdiff2_r[1]*rdvSq4-1.224744871391589*Ghat_r[1]*rdv2; 
  out[6] += -0.7071067811865475*Ghat_r[4]*rdv2; 
  out[7] += 1.060660171779821*nuVtSqSum_l[1]*fSkin[11]*rdvSq4+1.060660171779821*nuVtSqSum_l[0]*fSkin[7]*rdvSq4-0.6123724356957944*nuVtSqSum_l[1]*fSkin[6]*rdvSq4-0.6123724356957944*nuVtSqSum_l[0]*fSkin[3]*rdvSq4+1.224744871391589*Gdiff2_r[2]*rdvSq4-1.224744871391589*Ghat_r[2]*rdv2; 
  out[8] += -0.7071067811865475*Ghat_r[5]*rdv2; 
  out[9] += 1.060660171779821*nuVtSqSum_l[1]*fSkin[12]*rdvSq4+1.060660171779821*nuVtSqSum_l[0]*fSkin[9]*rdvSq4-0.6123724356957944*nuVtSqSum_l[1]*fSkin[8]*rdvSq4-0.6123724356957944*nuVtSqSum_l[0]*fSkin[4]*rdvSq4+1.224744871391589*Gdiff2_r[3]*rdvSq4-1.224744871391589*Ghat_r[3]*rdv2; 
  out[10] += -0.7071067811865475*Ghat_r[6]*rdv2; 
  out[11] += 1.060660171779821*nuVtSqSum_l[0]*fSkin[11]*rdvSq4+1.060660171779821*nuVtSqSum_l[1]*fSkin[7]*rdvSq4-0.6123724356957944*nuVtSqSum_l[0]*fSkin[6]*rdvSq4+1.224744871391589*Gdiff2_r[4]*rdvSq4-0.6123724356957944*nuVtSqSum_l[1]*fSkin[3]*rdvSq4-1.224744871391589*Ghat_r[4]*rdv2; 
  out[12] += 1.060660171779821*nuVtSqSum_l[0]*fSkin[12]*rdvSq4+1.060660171779821*nuVtSqSum_l[1]*fSkin[9]*rdvSq4-0.6123724356957944*nuVtSqSum_l[0]*fSkin[8]*rdvSq4+1.224744871391589*Gdiff2_r[5]*rdvSq4-0.6123724356957944*nuVtSqSum_l[1]*fSkin[4]*rdvSq4-1.224744871391589*Ghat_r[5]*rdv2; 
  out[13] += -0.7071067811865475*Ghat_r[7]*rdv2; 
  out[14] += 1.060660171779821*nuVtSqSum_l[1]*fSkin[15]*rdvSq4+1.060660171779821*nuVtSqSum_l[0]*fSkin[14]*rdvSq4-0.6123724356957944*nuVtSqSum_l[1]*fSkin[13]*rdvSq4-0.6123724356957944*nuVtSqSum_l[0]*fSkin[10]*rdvSq4+1.224744871391589*Gdiff2_r[6]*rdvSq4-1.224744871391589*Ghat_r[6]*rdv2; 
  out[15] += 1.060660171779821*nuVtSqSum_l[0]*fSkin[15]*rdvSq4+1.060660171779821*nuVtSqSum_l[1]*fSkin[14]*rdvSq4-0.6123724356957944*nuVtSqSum_l[0]*fSkin[13]*rdvSq4-0.6123724356957944*nuVtSqSum_l[1]*fSkin[10]*rdvSq4+1.224744871391589*Gdiff2_r[7]*rdvSq4-1.224744871391589*Ghat_r[7]*rdv2; 

  } else {

  const double *sumNuUx_l = &nuUSum[0]; 

  double alphaDrSurf_l[8]; 
  alphaDrSurf_l[0] = (2.82842712474619*w[1]+1.414213562373095*dxv[1])*nuSum-2.0*sumNuUx_l[0]; 
  alphaDrSurf_l[1] = -2.0*sumNuUx_l[1]; 

  double fUpwindQuad_l[8];
  if (alphaDrSurf_l[0]-alphaDrSurf_l[1] > 0) { 
    fUpwindQuad_l[0] = (-0.4330127018922193*fEdge[15])+0.4330127018922193*fEdge[14]-0.25*fEdge[13]+0.4330127018922193*(fEdge[12]+fEdge[11])+0.25*fEdge[10]-0.4330127018922193*fEdge[9]+0.25*fEdge[8]-0.4330127018922193*fEdge[7]+0.25*fEdge[6]-0.4330127018922193*fEdge[5]-0.25*(fEdge[4]+fEdge[3])+0.4330127018922193*fEdge[2]-0.25*fEdge[1]+0.25*fEdge[0]; 
  } else { 
    fUpwindQuad_l[0] = 0.4330127018922193*fSkin[15]-0.4330127018922193*fSkin[14]-0.25*fSkin[13]-0.4330127018922193*(fSkin[12]+fSkin[11])+0.25*fSkin[10]+0.4330127018922193*fSkin[9]+0.25*fSkin[8]+0.4330127018922193*fSkin[7]+0.25*fSkin[6]+0.4330127018922193*fSkin[5]-0.25*(fSkin[4]+fSkin[3])-0.4330127018922193*fSkin[2]-0.25*fSkin[1]+0.25*fSkin[0]; 
  } 
  if (alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[1] = 0.4330127018922193*(fEdge[15]+fEdge[14])+0.25*fEdge[13]-0.4330127018922193*(fEdge[12]+fEdge[11])+0.25*fEdge[10]-0.4330127018922193*fEdge[9]-0.25*fEdge[8]-0.4330127018922193*fEdge[7]-0.25*fEdge[6]+0.4330127018922193*fEdge[5]-0.25*(fEdge[4]+fEdge[3])+0.4330127018922193*fEdge[2]+0.25*(fEdge[1]+fEdge[0]); 
  } else { 
    fUpwindQuad_l[1] = (-0.4330127018922193*(fSkin[15]+fSkin[14]))+0.25*fSkin[13]+0.4330127018922193*(fSkin[12]+fSkin[11])+0.25*fSkin[10]+0.4330127018922193*fSkin[9]-0.25*fSkin[8]+0.4330127018922193*fSkin[7]-0.25*fSkin[6]-0.4330127018922193*fSkin[5]-0.25*(fSkin[4]+fSkin[3])-0.4330127018922193*fSkin[2]+0.25*(fSkin[1]+fSkin[0]); 
  } 
  if (alphaDrSurf_l[0]-alphaDrSurf_l[1] > 0) { 
    fUpwindQuad_l[2] = 0.4330127018922193*fEdge[15]-0.4330127018922193*fEdge[14]+0.25*fEdge[13]+0.4330127018922193*fEdge[12]-0.4330127018922193*fEdge[11]-0.25*fEdge[10]-0.4330127018922193*fEdge[9]+0.25*fEdge[8]+0.4330127018922193*fEdge[7]-0.25*fEdge[6]-0.4330127018922193*fEdge[5]-0.25*fEdge[4]+0.25*fEdge[3]+0.4330127018922193*fEdge[2]-0.25*fEdge[1]+0.25*fEdge[0]; 
  } else { 
    fUpwindQuad_l[2] = (-0.4330127018922193*fSkin[15])+0.4330127018922193*fSkin[14]+0.25*fSkin[13]-0.4330127018922193*fSkin[12]+0.4330127018922193*fSkin[11]-0.25*fSkin[10]+0.4330127018922193*fSkin[9]+0.25*fSkin[8]-0.4330127018922193*fSkin[7]-0.25*fSkin[6]+0.4330127018922193*fSkin[5]-0.25*fSkin[4]+0.25*fSkin[3]-0.4330127018922193*fSkin[2]-0.25*fSkin[1]+0.25*fSkin[0]; 
  } 
  if (alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[3] = (-0.4330127018922193*(fEdge[15]+fEdge[14]))-0.25*fEdge[13]-0.4330127018922193*fEdge[12]+0.4330127018922193*fEdge[11]-0.25*fEdge[10]-0.4330127018922193*fEdge[9]-0.25*fEdge[8]+0.4330127018922193*fEdge[7]+0.25*fEdge[6]+0.4330127018922193*fEdge[5]-0.25*fEdge[4]+0.25*fEdge[3]+0.4330127018922193*fEdge[2]+0.25*(fEdge[1]+fEdge[0]); 
  } else { 
    fUpwindQuad_l[3] = 0.4330127018922193*(fSkin[15]+fSkin[14])-0.25*fSkin[13]+0.4330127018922193*fSkin[12]-0.4330127018922193*fSkin[11]-0.25*fSkin[10]+0.4330127018922193*fSkin[9]-0.25*fSkin[8]-0.4330127018922193*fSkin[7]+0.25*fSkin[6]-0.4330127018922193*fSkin[5]-0.25*fSkin[4]+0.25*fSkin[3]-0.4330127018922193*fSkin[2]+0.25*(fSkin[1]+fSkin[0]); 
  } 
  if (alphaDrSurf_l[0]-alphaDrSurf_l[1] > 0) { 
    fUpwindQuad_l[4] = 0.4330127018922193*fEdge[15]-0.4330127018922193*fEdge[14]+0.25*fEdge[13]-0.4330127018922193*fEdge[12]+0.4330127018922193*fEdge[11]-0.25*fEdge[10]+0.4330127018922193*fEdge[9]-0.25*fEdge[8]-0.4330127018922193*fEdge[7]+0.25*fEdge[6]-0.4330127018922193*fEdge[5]+0.25*fEdge[4]-0.25*fEdge[3]+0.4330127018922193*fEdge[2]-0.25*fEdge[1]+0.25*fEdge[0]; 
  } else { 
    fUpwindQuad_l[4] = (-0.4330127018922193*fSkin[15])+0.4330127018922193*fSkin[14]+0.25*fSkin[13]+0.4330127018922193*fSkin[12]-0.4330127018922193*fSkin[11]-0.25*fSkin[10]-0.4330127018922193*fSkin[9]-0.25*fSkin[8]+0.4330127018922193*fSkin[7]+0.25*fSkin[6]+0.4330127018922193*fSkin[5]+0.25*fSkin[4]-0.25*fSkin[3]-0.4330127018922193*fSkin[2]-0.25*fSkin[1]+0.25*fSkin[0]; 
  } 
  if (alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[5] = (-0.4330127018922193*(fEdge[15]+fEdge[14]))-0.25*fEdge[13]+0.4330127018922193*fEdge[12]-0.4330127018922193*fEdge[11]-0.25*fEdge[10]+0.4330127018922193*fEdge[9]+0.25*fEdge[8]-0.4330127018922193*fEdge[7]-0.25*fEdge[6]+0.4330127018922193*fEdge[5]+0.25*fEdge[4]-0.25*fEdge[3]+0.4330127018922193*fEdge[2]+0.25*(fEdge[1]+fEdge[0]); 
  } else { 
    fUpwindQuad_l[5] = 0.4330127018922193*(fSkin[15]+fSkin[14])-0.25*fSkin[13]-0.4330127018922193*fSkin[12]+0.4330127018922193*fSkin[11]-0.25*fSkin[10]-0.4330127018922193*fSkin[9]+0.25*fSkin[8]+0.4330127018922193*fSkin[7]-0.25*fSkin[6]-0.4330127018922193*fSkin[5]+0.25*fSkin[4]-0.25*fSkin[3]-0.4330127018922193*fSkin[2]+0.25*(fSkin[1]+fSkin[0]); 
  } 
  if (alphaDrSurf_l[0]-alphaDrSurf_l[1] > 0) { 
    fUpwindQuad_l[6] = (-0.4330127018922193*fEdge[15])+0.4330127018922193*fEdge[14]-0.25*fEdge[13]-0.4330127018922193*(fEdge[12]+fEdge[11])+0.25*fEdge[10]+0.4330127018922193*fEdge[9]-0.25*fEdge[8]+0.4330127018922193*fEdge[7]-0.25*fEdge[6]-0.4330127018922193*fEdge[5]+0.25*(fEdge[4]+fEdge[3])+0.4330127018922193*fEdge[2]-0.25*fEdge[1]+0.25*fEdge[0]; 
  } else { 
    fUpwindQuad_l[6] = 0.4330127018922193*fSkin[15]-0.4330127018922193*fSkin[14]-0.25*fSkin[13]+0.4330127018922193*(fSkin[12]+fSkin[11])+0.25*fSkin[10]-0.4330127018922193*fSkin[9]-0.25*fSkin[8]-0.4330127018922193*fSkin[7]-0.25*fSkin[6]+0.4330127018922193*fSkin[5]+0.25*(fSkin[4]+fSkin[3])-0.4330127018922193*fSkin[2]-0.25*fSkin[1]+0.25*fSkin[0]; 
  } 
  if (alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[7] = 0.4330127018922193*(fEdge[15]+fEdge[14])+0.25*fEdge[13]+0.4330127018922193*(fEdge[12]+fEdge[11])+0.25*fEdge[10]+0.4330127018922193*fEdge[9]+0.25*fEdge[8]+0.4330127018922193*fEdge[7]+0.25*fEdge[6]+0.4330127018922193*fEdge[5]+0.25*(fEdge[4]+fEdge[3])+0.4330127018922193*fEdge[2]+0.25*(fEdge[1]+fEdge[0]); 
  } else { 
    fUpwindQuad_l[7] = (-0.4330127018922193*(fSkin[15]+fSkin[14]))+0.25*fSkin[13]-0.4330127018922193*(fSkin[12]+fSkin[11])+0.25*fSkin[10]-0.4330127018922193*fSkin[9]+0.25*fSkin[8]-0.4330127018922193*fSkin[7]+0.25*fSkin[6]-0.4330127018922193*fSkin[5]+0.25*(fSkin[4]+fSkin[3])-0.4330127018922193*fSkin[2]+0.25*(fSkin[1]+fSkin[0]); 
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

  double Gdiff_l[8]; 
  double Ghat_l[8]; 


  Gdiff_l[0] = -0.0625*(8.660254037844386*nuVtSqSum_l[1]*fSkin[5]+8.660254037844386*nuVtSqSum_l[1]*fEdge[5]+8.660254037844386*nuVtSqSum_l[0]*fSkin[2]+8.660254037844386*nuVtSqSum_l[0]*fEdge[2]+(9.0*fEdge[1]-9.0*fSkin[1])*nuVtSqSum_l[1]+(9.0*fEdge[0]-9.0*fSkin[0])*nuVtSqSum_l[0]); 
  Gdiff_l[1] = -0.0625*(8.660254037844386*nuVtSqSum_l[0]*fSkin[5]+8.660254037844386*nuVtSqSum_l[0]*fEdge[5]+8.660254037844386*nuVtSqSum_l[1]*fSkin[2]+8.660254037844386*nuVtSqSum_l[1]*fEdge[2]+(9.0*fEdge[0]-9.0*fSkin[0])*nuVtSqSum_l[1]-9.0*nuVtSqSum_l[0]*fSkin[1]+9.0*nuVtSqSum_l[0]*fEdge[1]); 
  Gdiff_l[2] = -0.0625*(8.660254037844386*nuVtSqSum_l[1]*fSkin[11]+8.660254037844386*nuVtSqSum_l[1]*fEdge[11]+8.660254037844386*nuVtSqSum_l[0]*fSkin[7]+8.660254037844386*nuVtSqSum_l[0]*fEdge[7]-9.0*nuVtSqSum_l[1]*fSkin[6]+9.0*nuVtSqSum_l[1]*fEdge[6]-9.0*nuVtSqSum_l[0]*fSkin[3]+9.0*nuVtSqSum_l[0]*fEdge[3]); 
  Gdiff_l[3] = -0.0625*(8.660254037844386*nuVtSqSum_l[1]*fSkin[12]+8.660254037844386*nuVtSqSum_l[1]*fEdge[12]+8.660254037844386*nuVtSqSum_l[0]*fSkin[9]+8.660254037844386*nuVtSqSum_l[0]*fEdge[9]-9.0*nuVtSqSum_l[1]*fSkin[8]+9.0*nuVtSqSum_l[1]*fEdge[8]-9.0*nuVtSqSum_l[0]*fSkin[4]+9.0*nuVtSqSum_l[0]*fEdge[4]); 
  Gdiff_l[4] = -0.0625*(8.660254037844386*nuVtSqSum_l[0]*fSkin[11]+8.660254037844386*nuVtSqSum_l[0]*fEdge[11]+8.660254037844386*nuVtSqSum_l[1]*fSkin[7]+8.660254037844386*nuVtSqSum_l[1]*fEdge[7]-9.0*nuVtSqSum_l[0]*fSkin[6]+9.0*nuVtSqSum_l[0]*fEdge[6]-9.0*nuVtSqSum_l[1]*fSkin[3]+9.0*nuVtSqSum_l[1]*fEdge[3]); 
  Gdiff_l[5] = -0.0625*(8.660254037844386*nuVtSqSum_l[0]*fSkin[12]+8.660254037844386*nuVtSqSum_l[0]*fEdge[12]+8.660254037844386*nuVtSqSum_l[1]*fSkin[9]+8.660254037844386*nuVtSqSum_l[1]*fEdge[9]-9.0*nuVtSqSum_l[0]*fSkin[8]+9.0*nuVtSqSum_l[0]*fEdge[8]-9.0*nuVtSqSum_l[1]*fSkin[4]+9.0*nuVtSqSum_l[1]*fEdge[4]); 
  Gdiff_l[6] = -0.0625*(8.660254037844386*nuVtSqSum_l[1]*fSkin[15]+8.660254037844386*nuVtSqSum_l[1]*fEdge[15]+8.660254037844386*nuVtSqSum_l[0]*fSkin[14]+8.660254037844386*nuVtSqSum_l[0]*fEdge[14]-9.0*nuVtSqSum_l[1]*fSkin[13]+9.0*nuVtSqSum_l[1]*fEdge[13]-9.0*nuVtSqSum_l[0]*fSkin[10]+9.0*nuVtSqSum_l[0]*fEdge[10]); 
  Gdiff_l[7] = -0.0625*(8.660254037844386*nuVtSqSum_l[0]*fSkin[15]+8.660254037844386*nuVtSqSum_l[0]*fEdge[15]+8.660254037844386*nuVtSqSum_l[1]*fSkin[14]+8.660254037844386*nuVtSqSum_l[1]*fEdge[14]-9.0*nuVtSqSum_l[0]*fSkin[13]+9.0*nuVtSqSum_l[0]*fEdge[13]-9.0*nuVtSqSum_l[1]*fSkin[10]+9.0*nuVtSqSum_l[1]*fEdge[10]); 

  Ghat_l[0] = Gdiff_l[0]*rdv2+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[1]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[0]; 
  Ghat_l[1] = Gdiff_l[1]*rdv2+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[1]+0.3535533905932737*fUpwind_l[0]*alphaDrSurf_l[1]; 
  Ghat_l[2] = Gdiff_l[2]*rdv2+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[4]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[2]; 
  Ghat_l[3] = Gdiff_l[3]*rdv2+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[5]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[3]; 
  Ghat_l[4] = Gdiff_l[4]*rdv2+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[4]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[2]; 
  Ghat_l[5] = Gdiff_l[5]*rdv2+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[5]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[3]; 
  Ghat_l[6] = Gdiff_l[6]*rdv2+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[7]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[6]; 
  Ghat_l[7] = Gdiff_l[7]*rdv2+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[7]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[6]; 

  out[0] += 0.7071067811865475*Ghat_l[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat_l[1]*rdv2; 
  out[2] += 1.060660171779821*nuVtSqSum_r[1]*fSkin[5]*rdvSq4+1.060660171779821*nuVtSqSum_r[0]*fSkin[2]*rdvSq4+0.6123724356957944*fSkin[1]*nuVtSqSum_r[1]*rdvSq4+0.6123724356957944*fSkin[0]*nuVtSqSum_r[0]*rdvSq4-1.224744871391589*Gdiff2_l[0]*rdvSq4-1.224744871391589*Ghat_l[0]*rdv2; 
  out[3] += 0.7071067811865475*Ghat_l[2]*rdv2; 
  out[4] += 0.7071067811865475*Ghat_l[3]*rdv2; 
  out[5] += 1.060660171779821*nuVtSqSum_r[0]*fSkin[5]*rdvSq4+1.060660171779821*nuVtSqSum_r[1]*fSkin[2]*rdvSq4+0.6123724356957944*fSkin[0]*nuVtSqSum_r[1]*rdvSq4+0.6123724356957944*nuVtSqSum_r[0]*fSkin[1]*rdvSq4-1.224744871391589*Gdiff2_l[1]*rdvSq4-1.224744871391589*Ghat_l[1]*rdv2; 
  out[6] += 0.7071067811865475*Ghat_l[4]*rdv2; 
  out[7] += 1.060660171779821*nuVtSqSum_r[1]*fSkin[11]*rdvSq4+1.060660171779821*nuVtSqSum_r[0]*fSkin[7]*rdvSq4+0.6123724356957944*nuVtSqSum_r[1]*fSkin[6]*rdvSq4+0.6123724356957944*nuVtSqSum_r[0]*fSkin[3]*rdvSq4-1.224744871391589*Gdiff2_l[2]*rdvSq4-1.224744871391589*Ghat_l[2]*rdv2; 
  out[8] += 0.7071067811865475*Ghat_l[5]*rdv2; 
  out[9] += 1.060660171779821*nuVtSqSum_r[1]*fSkin[12]*rdvSq4+1.060660171779821*nuVtSqSum_r[0]*fSkin[9]*rdvSq4+0.6123724356957944*nuVtSqSum_r[1]*fSkin[8]*rdvSq4+0.6123724356957944*nuVtSqSum_r[0]*fSkin[4]*rdvSq4-1.224744871391589*Gdiff2_l[3]*rdvSq4-1.224744871391589*Ghat_l[3]*rdv2; 
  out[10] += 0.7071067811865475*Ghat_l[6]*rdv2; 
  out[11] += 1.060660171779821*nuVtSqSum_r[0]*fSkin[11]*rdvSq4+1.060660171779821*nuVtSqSum_r[1]*fSkin[7]*rdvSq4+0.6123724356957944*nuVtSqSum_r[0]*fSkin[6]*rdvSq4-1.224744871391589*Gdiff2_l[4]*rdvSq4+0.6123724356957944*nuVtSqSum_r[1]*fSkin[3]*rdvSq4-1.224744871391589*Ghat_l[4]*rdv2; 
  out[12] += 1.060660171779821*nuVtSqSum_r[0]*fSkin[12]*rdvSq4+1.060660171779821*nuVtSqSum_r[1]*fSkin[9]*rdvSq4+0.6123724356957944*nuVtSqSum_r[0]*fSkin[8]*rdvSq4-1.224744871391589*Gdiff2_l[5]*rdvSq4+0.6123724356957944*nuVtSqSum_r[1]*fSkin[4]*rdvSq4-1.224744871391589*Ghat_l[5]*rdv2; 
  out[13] += 0.7071067811865475*Ghat_l[7]*rdv2; 
  out[14] += 1.060660171779821*nuVtSqSum_r[1]*fSkin[15]*rdvSq4+1.060660171779821*nuVtSqSum_r[0]*fSkin[14]*rdvSq4+0.6123724356957944*nuVtSqSum_r[1]*fSkin[13]*rdvSq4+0.6123724356957944*nuVtSqSum_r[0]*fSkin[10]*rdvSq4-1.224744871391589*Gdiff2_l[6]*rdvSq4-1.224744871391589*Ghat_l[6]*rdv2; 
  out[15] += 1.060660171779821*nuVtSqSum_r[0]*fSkin[15]*rdvSq4+1.060660171779821*nuVtSqSum_r[1]*fSkin[14]*rdvSq4+0.6123724356957944*nuVtSqSum_r[0]*fSkin[13]*rdvSq4+0.6123724356957944*nuVtSqSum_r[1]*fSkin[10]*rdvSq4-1.224744871391589*Gdiff2_l[7]*rdvSq4-1.224744871391589*Ghat_l[7]*rdv2; 
  } 
} 
