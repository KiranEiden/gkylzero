#include <gkyl_vlasov_lbo_kernels.h> 
#include <gkyl_basis_ser_1x2v_p2_surfvy_quad.h> 
GKYL_CU_DH void vlasov_lbo_boundary_surfvy_1x2v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[3]:         Cell-center coordinates. 
  // dxv[3]:       Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[6]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]: sum of thermal speeds squared time their respective collisionalities. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[2]; 
  double rdvSq4 = 4.0/(dxv[2]*dxv[2]); 

  const double *sumNuUy = &nuUSum[3]; 

  double alphaDrSurf[8] = {0.0}; 
  double fUpwindQuad[9] = {0.0};
  double fUpwind[8] = {0.0};;
  double Ghat[8] = {0.0}; 
  double Gdiff[8] = {0.0}; 
  double Gdiff2[8] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*(2.828427124746191*w[2]+1.414213562373095*dxv[2])-2.828427124746191*sumNuUy[0]); 
  alphaDrSurf[1] = 0.5*(nuSum[1]*(2.828427124746191*w[2]+1.414213562373095*dxv[2])-2.828427124746191*sumNuUy[1]); 
  alphaDrSurf[4] = 0.5*(2.828427124746191*nuSum[2]*w[2]-2.828427124746191*sumNuUy[2]+1.414213562373095*dxv[2]*nuSum[2]); 

  if (0.4472135954999579*alphaDrSurf[4]-0.6708203932499369*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = ser_1x2v_p2_surfvy_quad_0(1, fSkin); 
  } else { 
    fUpwindQuad[0] = ser_1x2v_p2_surfvy_quad_0(-1, fEdge); 
  } 
  if (0.5*alphaDrSurf[0]-0.5590169943749475*alphaDrSurf[4] < 0) { 
    fUpwindQuad[1] = ser_1x2v_p2_surfvy_quad_1(1, fSkin); 
  } else { 
    fUpwindQuad[1] = ser_1x2v_p2_surfvy_quad_1(-1, fEdge); 
  } 
  if (0.4472135954999579*alphaDrSurf[4]+0.6708203932499369*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[2] = ser_1x2v_p2_surfvy_quad_2(1, fSkin); 
  } else { 
    fUpwindQuad[2] = ser_1x2v_p2_surfvy_quad_2(-1, fEdge); 
  } 
  if (0.4472135954999579*alphaDrSurf[4]-0.6708203932499369*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = ser_1x2v_p2_surfvy_quad_3(1, fSkin); 
  } else { 
    fUpwindQuad[3] = ser_1x2v_p2_surfvy_quad_3(-1, fEdge); 
  } 
  if (0.5*alphaDrSurf[0]-0.5590169943749475*alphaDrSurf[4] < 0) { 
    fUpwindQuad[4] = ser_1x2v_p2_surfvy_quad_4(1, fSkin); 
  } else { 
    fUpwindQuad[4] = ser_1x2v_p2_surfvy_quad_4(-1, fEdge); 
  } 
  if (0.4472135954999579*alphaDrSurf[4]+0.6708203932499369*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[5] = ser_1x2v_p2_surfvy_quad_5(1, fSkin); 
  } else { 
    fUpwindQuad[5] = ser_1x2v_p2_surfvy_quad_5(-1, fEdge); 
  } 
  if (0.4472135954999579*alphaDrSurf[4]-0.6708203932499369*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[6] = ser_1x2v_p2_surfvy_quad_6(1, fSkin); 
  } else { 
    fUpwindQuad[6] = ser_1x2v_p2_surfvy_quad_6(-1, fEdge); 
  } 
  if (0.5*alphaDrSurf[0]-0.5590169943749475*alphaDrSurf[4] < 0) { 
    fUpwindQuad[7] = ser_1x2v_p2_surfvy_quad_7(1, fSkin); 
  } else { 
    fUpwindQuad[7] = ser_1x2v_p2_surfvy_quad_7(-1, fEdge); 
  } 
  if (0.4472135954999579*alphaDrSurf[4]+0.6708203932499369*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[8] = ser_1x2v_p2_surfvy_quad_8(1, fSkin); 
  } else { 
    fUpwindQuad[8] = ser_1x2v_p2_surfvy_quad_8(-1, fEdge); 
  } 

  fUpwind[0] = 0.154320987654321*fUpwindQuad[8]+0.2469135802469136*fUpwindQuad[7]+0.154320987654321*fUpwindQuad[6]+0.2469135802469136*fUpwindQuad[5]+0.3950617283950617*fUpwindQuad[4]+0.2469135802469136*fUpwindQuad[3]+0.154320987654321*fUpwindQuad[2]+0.2469135802469136*fUpwindQuad[1]+0.154320987654321*fUpwindQuad[0]; 
  fUpwind[1] = 0.2070433312499806*fUpwindQuad[8]-0.2070433312499806*fUpwindQuad[6]+0.3312693299999688*fUpwindQuad[5]-0.3312693299999688*fUpwindQuad[3]+0.2070433312499806*fUpwindQuad[2]-0.2070433312499806*fUpwindQuad[0]; 
  fUpwind[2] = 0.2070433312499806*fUpwindQuad[8]+0.3312693299999688*fUpwindQuad[7]+0.2070433312499806*fUpwindQuad[6]-0.2070433312499806*fUpwindQuad[2]-0.3312693299999688*fUpwindQuad[1]-0.2070433312499806*fUpwindQuad[0]; 
  fUpwind[3] = 0.2777777777777778*fUpwindQuad[8]-0.2777777777777778*fUpwindQuad[6]-0.2777777777777778*fUpwindQuad[2]+0.2777777777777778*fUpwindQuad[0]; 
  fUpwind[4] = 0.138028887499987*fUpwindQuad[8]-0.2760577749999741*fUpwindQuad[7]+0.138028887499987*fUpwindQuad[6]+0.2208462199999792*fUpwindQuad[5]-0.4416924399999584*fUpwindQuad[4]+0.2208462199999792*fUpwindQuad[3]+0.138028887499987*fUpwindQuad[2]-0.2760577749999741*fUpwindQuad[1]+0.138028887499987*fUpwindQuad[0]; 
  fUpwind[5] = 0.138028887499987*fUpwindQuad[8]+0.2208462199999792*fUpwindQuad[7]+0.138028887499987*fUpwindQuad[6]-0.2760577749999741*fUpwindQuad[5]-0.4416924399999584*fUpwindQuad[4]-0.2760577749999741*fUpwindQuad[3]+0.138028887499987*fUpwindQuad[2]+0.2208462199999792*fUpwindQuad[1]+0.138028887499987*fUpwindQuad[0]; 
  fUpwind[6] = 0.1851851851851853*fUpwindQuad[8]-0.3703703703703705*fUpwindQuad[7]+0.1851851851851853*fUpwindQuad[6]-0.1851851851851853*fUpwindQuad[2]+0.3703703703703705*fUpwindQuad[1]-0.1851851851851853*fUpwindQuad[0]; 
  fUpwind[7] = 0.1851851851851853*fUpwindQuad[8]-0.1851851851851853*fUpwindQuad[6]-0.3703703703703705*fUpwindQuad[5]+0.3703703703703705*fUpwindQuad[3]+0.1851851851851853*fUpwindQuad[2]-0.1851851851851853*fUpwindQuad[0]; 

  Gdiff2[0] = 0.2445699350390395*nuVtSqSum[1]*fSkin[15]+0.2445699350390395*nuVtSqSum[1]*fEdge[15]+0.3518228202874282*nuVtSqSum[2]*fSkin[13]-0.3518228202874282*nuVtSqSum[2]*fEdge[13]+0.2445699350390395*nuVtSqSum[0]*fSkin[9]+0.2445699350390395*nuVtSqSum[0]*fEdge[9]+0.25*nuVtSqSum[2]*fSkin[7]+0.25*nuVtSqSum[2]*fEdge[7]+0.3518228202874282*nuVtSqSum[1]*fSkin[5]-0.3518228202874282*nuVtSqSum[1]*fEdge[5]+0.3518228202874282*nuVtSqSum[0]*fSkin[3]-0.3518228202874282*nuVtSqSum[0]*fEdge[3]+0.25*fSkin[1]*nuVtSqSum[1]+0.25*fEdge[1]*nuVtSqSum[1]+0.25*fSkin[0]*nuVtSqSum[0]+0.25*fEdge[0]*nuVtSqSum[0]; 
  Gdiff2[1] = 0.21875*nuVtSqSum[2]*fSkin[15]+0.2445699350390395*nuVtSqSum[0]*fSkin[15]+0.21875*nuVtSqSum[2]*fEdge[15]+0.2445699350390395*nuVtSqSum[0]*fEdge[15]+0.3146798968793527*nuVtSqSum[1]*fSkin[13]-0.3146798968793527*nuVtSqSum[1]*fEdge[13]+0.2445699350390395*nuVtSqSum[1]*fSkin[9]+0.2445699350390395*nuVtSqSum[1]*fEdge[9]+0.223606797749979*nuVtSqSum[1]*fSkin[7]+0.223606797749979*nuVtSqSum[1]*fEdge[7]+0.3146798968793526*nuVtSqSum[2]*fSkin[5]+0.3518228202874282*nuVtSqSum[0]*fSkin[5]-0.3146798968793526*nuVtSqSum[2]*fEdge[5]-0.3518228202874282*nuVtSqSum[0]*fEdge[5]+0.3518228202874282*nuVtSqSum[1]*fSkin[3]-0.3518228202874282*nuVtSqSum[1]*fEdge[3]+0.223606797749979*fSkin[1]*nuVtSqSum[2]+0.223606797749979*fEdge[1]*nuVtSqSum[2]+0.25*fSkin[0]*nuVtSqSum[1]+0.25*fEdge[0]*nuVtSqSum[1]+0.25*nuVtSqSum[0]*fSkin[1]+0.25*nuVtSqSum[0]*fEdge[1]; 
  Gdiff2[2] = 0.2445699350390395*nuVtSqSum[1]*fSkin[19]+0.2445699350390395*nuVtSqSum[1]*fEdge[19]+0.3518228202874282*nuVtSqSum[2]*fSkin[17]-0.3518228202874282*nuVtSqSum[2]*fEdge[17]+0.2445699350390395*nuVtSqSum[0]*fSkin[16]+0.2445699350390395*nuVtSqSum[0]*fEdge[16]+0.2500000000000001*nuVtSqSum[2]*fSkin[11]+0.2500000000000001*nuVtSqSum[2]*fEdge[11]+0.3518228202874282*nuVtSqSum[1]*fSkin[10]-0.3518228202874282*nuVtSqSum[1]*fEdge[10]+0.3518228202874282*nuVtSqSum[0]*fSkin[6]-0.3518228202874282*nuVtSqSum[0]*fEdge[6]+0.25*nuVtSqSum[1]*fSkin[4]+0.25*nuVtSqSum[1]*fEdge[4]+0.25*nuVtSqSum[0]*fSkin[2]+0.25*nuVtSqSum[0]*fEdge[2]; 
  Gdiff2[3] = 0.21875*nuVtSqSum[2]*fSkin[19]+0.2445699350390395*nuVtSqSum[0]*fSkin[19]+0.21875*nuVtSqSum[2]*fEdge[19]+0.2445699350390395*nuVtSqSum[0]*fEdge[19]+0.3146798968793526*nuVtSqSum[1]*fSkin[17]-0.3146798968793526*nuVtSqSum[1]*fEdge[17]+0.2445699350390395*nuVtSqSum[1]*fSkin[16]+0.2445699350390395*nuVtSqSum[1]*fEdge[16]+0.223606797749979*nuVtSqSum[1]*fSkin[11]+0.223606797749979*nuVtSqSum[1]*fEdge[11]+0.3146798968793526*nuVtSqSum[2]*fSkin[10]+0.3518228202874282*nuVtSqSum[0]*fSkin[10]-0.3146798968793526*nuVtSqSum[2]*fEdge[10]-0.3518228202874282*nuVtSqSum[0]*fEdge[10]+0.3518228202874282*nuVtSqSum[1]*fSkin[6]-0.3518228202874282*nuVtSqSum[1]*fEdge[6]+0.223606797749979*nuVtSqSum[2]*fSkin[4]+0.25*nuVtSqSum[0]*fSkin[4]+0.223606797749979*nuVtSqSum[2]*fEdge[4]+0.25*nuVtSqSum[0]*fEdge[4]+0.25*nuVtSqSum[1]*fSkin[2]+0.25*nuVtSqSum[1]*fEdge[2]; 
  Gdiff2[4] = 0.21875*nuVtSqSum[1]*fSkin[15]+0.21875*nuVtSqSum[1]*fEdge[15]+0.2247713549138233*nuVtSqSum[2]*fSkin[13]+0.3518228202874282*nuVtSqSum[0]*fSkin[13]-0.2247713549138233*nuVtSqSum[2]*fEdge[13]-0.3518228202874282*nuVtSqSum[0]*fEdge[13]+0.2445699350390395*nuVtSqSum[2]*fSkin[9]+0.2445699350390395*nuVtSqSum[2]*fEdge[9]+0.159719141249985*nuVtSqSum[2]*fSkin[7]+0.25*nuVtSqSum[0]*fSkin[7]+0.159719141249985*nuVtSqSum[2]*fEdge[7]+0.25*nuVtSqSum[0]*fEdge[7]+0.3146798968793526*nuVtSqSum[1]*fSkin[5]-0.3146798968793526*nuVtSqSum[1]*fEdge[5]+0.3518228202874282*nuVtSqSum[2]*fSkin[3]-0.3518228202874282*nuVtSqSum[2]*fEdge[3]+0.25*fSkin[0]*nuVtSqSum[2]+0.25*fEdge[0]*nuVtSqSum[2]+0.223606797749979*fSkin[1]*nuVtSqSum[1]+0.223606797749979*fEdge[1]*nuVtSqSum[1]; 
  Gdiff2[5] = 0.3518228202874282*nuVtSqSum[1]*fSkin[18]-0.3518228202874282*nuVtSqSum[1]*fEdge[18]+0.3518228202874282*nuVtSqSum[0]*fSkin[14]-0.3518228202874282*nuVtSqSum[0]*fEdge[14]+0.2500000000000001*nuVtSqSum[1]*fSkin[12]+0.2500000000000001*nuVtSqSum[1]*fEdge[12]+0.25*nuVtSqSum[0]*fSkin[8]+0.25*nuVtSqSum[0]*fEdge[8]; 
  Gdiff2[6] = 0.21875*nuVtSqSum[1]*fSkin[19]+0.21875*nuVtSqSum[1]*fEdge[19]+0.2247713549138233*nuVtSqSum[2]*fSkin[17]+0.3518228202874282*nuVtSqSum[0]*fSkin[17]-0.2247713549138233*nuVtSqSum[2]*fEdge[17]-0.3518228202874282*nuVtSqSum[0]*fEdge[17]+0.2445699350390395*nuVtSqSum[2]*fSkin[16]+0.2445699350390395*nuVtSqSum[2]*fEdge[16]+0.159719141249985*nuVtSqSum[2]*fSkin[11]+0.25*nuVtSqSum[0]*fSkin[11]+0.159719141249985*nuVtSqSum[2]*fEdge[11]+0.25*nuVtSqSum[0]*fEdge[11]+0.3146798968793527*nuVtSqSum[1]*fSkin[10]-0.3146798968793527*nuVtSqSum[1]*fEdge[10]+0.3518228202874282*nuVtSqSum[2]*fSkin[6]-0.3518228202874282*nuVtSqSum[2]*fEdge[6]+0.223606797749979*nuVtSqSum[1]*fSkin[4]+0.223606797749979*nuVtSqSum[1]*fEdge[4]+0.2500000000000001*fSkin[2]*nuVtSqSum[2]+0.2500000000000001*fEdge[2]*nuVtSqSum[2]; 
  Gdiff2[7] = 0.3146798968793527*nuVtSqSum[2]*fSkin[18]+0.3518228202874282*nuVtSqSum[0]*fSkin[18]-0.3146798968793527*nuVtSqSum[2]*fEdge[18]-0.3518228202874282*nuVtSqSum[0]*fEdge[18]+0.3518228202874282*nuVtSqSum[1]*fSkin[14]-0.3518228202874282*nuVtSqSum[1]*fEdge[14]+0.223606797749979*nuVtSqSum[2]*fSkin[12]+0.25*nuVtSqSum[0]*fSkin[12]+0.223606797749979*nuVtSqSum[2]*fEdge[12]+0.25*nuVtSqSum[0]*fEdge[12]+0.2500000000000001*nuVtSqSum[1]*fSkin[8]+0.2500000000000001*nuVtSqSum[1]*fEdge[8]; 

  Gdiff[0] = (-0.6708203932499369*nuVtSqSum[1]*fSkin[15])+0.6708203932499369*nuVtSqSum[1]*fEdge[15]-1.190784930203603*nuVtSqSum[2]*fSkin[13]-1.190784930203603*nuVtSqSum[2]*fEdge[13]-0.6708203932499369*nuVtSqSum[0]*fSkin[9]+0.6708203932499369*nuVtSqSum[0]*fEdge[9]-0.9375*nuVtSqSum[2]*fSkin[7]+0.9375*nuVtSqSum[2]*fEdge[7]-1.190784930203603*nuVtSqSum[1]*fSkin[5]-1.190784930203603*nuVtSqSum[1]*fEdge[5]-1.190784930203603*nuVtSqSum[0]*fSkin[3]-1.190784930203603*nuVtSqSum[0]*fEdge[3]-0.9375*fSkin[1]*nuVtSqSum[1]+0.9375*fEdge[1]*nuVtSqSum[1]-0.9375*fSkin[0]*nuVtSqSum[0]+0.9375*fEdge[0]*nuVtSqSum[0]; 
  Gdiff[1] = (-0.5999999999999999*nuVtSqSum[2]*fSkin[15])-0.6708203932499369*nuVtSqSum[0]*fSkin[15]+0.5999999999999999*nuVtSqSum[2]*fEdge[15]+0.6708203932499369*nuVtSqSum[0]*fEdge[15]-1.06507042020704*nuVtSqSum[1]*fSkin[13]-1.06507042020704*nuVtSqSum[1]*fEdge[13]-0.6708203932499369*nuVtSqSum[1]*fSkin[9]+0.6708203932499369*nuVtSqSum[1]*fEdge[9]-0.8385254915624212*nuVtSqSum[1]*fSkin[7]+0.8385254915624212*nuVtSqSum[1]*fEdge[7]-1.06507042020704*nuVtSqSum[2]*fSkin[5]-1.190784930203603*nuVtSqSum[0]*fSkin[5]-1.06507042020704*nuVtSqSum[2]*fEdge[5]-1.190784930203603*nuVtSqSum[0]*fEdge[5]-1.190784930203603*nuVtSqSum[1]*fSkin[3]-1.190784930203603*nuVtSqSum[1]*fEdge[3]-0.8385254915624212*fSkin[1]*nuVtSqSum[2]+0.8385254915624212*fEdge[1]*nuVtSqSum[2]-0.9375*fSkin[0]*nuVtSqSum[1]+0.9375*fEdge[0]*nuVtSqSum[1]-0.9375*nuVtSqSum[0]*fSkin[1]+0.9375*nuVtSqSum[0]*fEdge[1]; 
  Gdiff[2] = (-0.6708203932499369*nuVtSqSum[1]*fSkin[19])+0.6708203932499369*nuVtSqSum[1]*fEdge[19]-1.190784930203603*nuVtSqSum[2]*fSkin[17]-1.190784930203603*nuVtSqSum[2]*fEdge[17]-0.6708203932499369*nuVtSqSum[0]*fSkin[16]+0.6708203932499369*nuVtSqSum[0]*fEdge[16]-0.9375000000000001*nuVtSqSum[2]*fSkin[11]+0.9375000000000001*nuVtSqSum[2]*fEdge[11]-1.190784930203603*nuVtSqSum[1]*fSkin[10]-1.190784930203603*nuVtSqSum[1]*fEdge[10]-1.190784930203603*nuVtSqSum[0]*fSkin[6]-1.190784930203603*nuVtSqSum[0]*fEdge[6]-0.9375*nuVtSqSum[1]*fSkin[4]+0.9375*nuVtSqSum[1]*fEdge[4]-0.9375*nuVtSqSum[0]*fSkin[2]+0.9375*nuVtSqSum[0]*fEdge[2]; 
  Gdiff[3] = (-0.6*nuVtSqSum[2]*fSkin[19])-0.6708203932499369*nuVtSqSum[0]*fSkin[19]+0.6*nuVtSqSum[2]*fEdge[19]+0.6708203932499369*nuVtSqSum[0]*fEdge[19]-1.06507042020704*nuVtSqSum[1]*fSkin[17]-1.06507042020704*nuVtSqSum[1]*fEdge[17]-0.6708203932499369*nuVtSqSum[1]*fSkin[16]+0.6708203932499369*nuVtSqSum[1]*fEdge[16]-0.8385254915624211*nuVtSqSum[1]*fSkin[11]+0.8385254915624211*nuVtSqSum[1]*fEdge[11]-1.06507042020704*nuVtSqSum[2]*fSkin[10]-1.190784930203603*nuVtSqSum[0]*fSkin[10]-1.06507042020704*nuVtSqSum[2]*fEdge[10]-1.190784930203603*nuVtSqSum[0]*fEdge[10]-1.190784930203603*nuVtSqSum[1]*fSkin[6]-1.190784930203603*nuVtSqSum[1]*fEdge[6]-0.8385254915624212*nuVtSqSum[2]*fSkin[4]-0.9375*nuVtSqSum[0]*fSkin[4]+0.8385254915624212*nuVtSqSum[2]*fEdge[4]+0.9375*nuVtSqSum[0]*fEdge[4]-0.9375*nuVtSqSum[1]*fSkin[2]+0.9375*nuVtSqSum[1]*fEdge[2]; 
  Gdiff[4] = (-0.5999999999999999*nuVtSqSum[1]*fSkin[15])+0.5999999999999999*nuVtSqSum[1]*fEdge[15]-0.7607645858621712*nuVtSqSum[2]*fSkin[13]-1.190784930203603*nuVtSqSum[0]*fSkin[13]-0.7607645858621712*nuVtSqSum[2]*fEdge[13]-1.190784930203603*nuVtSqSum[0]*fEdge[13]-0.6708203932499369*nuVtSqSum[2]*fSkin[9]+0.6708203932499369*nuVtSqSum[2]*fEdge[9]-0.5989467796874438*nuVtSqSum[2]*fSkin[7]-0.9375*nuVtSqSum[0]*fSkin[7]+0.5989467796874438*nuVtSqSum[2]*fEdge[7]+0.9375*nuVtSqSum[0]*fEdge[7]-1.06507042020704*nuVtSqSum[1]*fSkin[5]-1.06507042020704*nuVtSqSum[1]*fEdge[5]-1.190784930203603*nuVtSqSum[2]*fSkin[3]-1.190784930203603*nuVtSqSum[2]*fEdge[3]-0.9375*fSkin[0]*nuVtSqSum[2]+0.9375*fEdge[0]*nuVtSqSum[2]-0.8385254915624212*fSkin[1]*nuVtSqSum[1]+0.8385254915624212*fEdge[1]*nuVtSqSum[1]; 
  Gdiff[5] = (-1.190784930203603*nuVtSqSum[1]*fSkin[18])-1.190784930203603*nuVtSqSum[1]*fEdge[18]-1.190784930203603*nuVtSqSum[0]*fSkin[14]-1.190784930203603*nuVtSqSum[0]*fEdge[14]-0.9375000000000001*nuVtSqSum[1]*fSkin[12]+0.9375000000000001*nuVtSqSum[1]*fEdge[12]-0.9375*nuVtSqSum[0]*fSkin[8]+0.9375*nuVtSqSum[0]*fEdge[8]; 
  Gdiff[6] = (-0.5999999999999999*nuVtSqSum[1]*fSkin[19])+0.5999999999999999*nuVtSqSum[1]*fEdge[19]-0.7607645858621712*nuVtSqSum[2]*fSkin[17]-1.190784930203603*nuVtSqSum[0]*fSkin[17]-0.7607645858621712*nuVtSqSum[2]*fEdge[17]-1.190784930203603*nuVtSqSum[0]*fEdge[17]-0.6708203932499369*nuVtSqSum[2]*fSkin[16]+0.6708203932499369*nuVtSqSum[2]*fEdge[16]-0.5989467796874438*nuVtSqSum[2]*fSkin[11]-0.9375*nuVtSqSum[0]*fSkin[11]+0.5989467796874438*nuVtSqSum[2]*fEdge[11]+0.9375*nuVtSqSum[0]*fEdge[11]-1.06507042020704*nuVtSqSum[1]*fSkin[10]-1.06507042020704*nuVtSqSum[1]*fEdge[10]-1.190784930203603*nuVtSqSum[2]*fSkin[6]-1.190784930203603*nuVtSqSum[2]*fEdge[6]-0.8385254915624211*nuVtSqSum[1]*fSkin[4]+0.8385254915624211*nuVtSqSum[1]*fEdge[4]-0.9375000000000001*fSkin[2]*nuVtSqSum[2]+0.9375000000000001*fEdge[2]*nuVtSqSum[2]; 
  Gdiff[7] = (-1.06507042020704*nuVtSqSum[2]*fSkin[18])-1.190784930203603*nuVtSqSum[0]*fSkin[18]-1.06507042020704*nuVtSqSum[2]*fEdge[18]-1.190784930203603*nuVtSqSum[0]*fEdge[18]-1.190784930203603*nuVtSqSum[1]*fSkin[14]-1.190784930203603*nuVtSqSum[1]*fEdge[14]-0.8385254915624212*nuVtSqSum[2]*fSkin[12]-0.9375*nuVtSqSum[0]*fSkin[12]+0.8385254915624212*nuVtSqSum[2]*fEdge[12]+0.9375*nuVtSqSum[0]*fEdge[12]-0.9375000000000001*nuVtSqSum[1]*fSkin[8]+0.9375000000000001*nuVtSqSum[1]*fEdge[8]; 

  Ghat[0] = Gdiff[0]*rdv2+0.5*alphaDrSurf[4]*fUpwind[4]+0.5*alphaDrSurf[1]*fUpwind[1]+0.5*alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = Gdiff[1]*rdv2+0.4472135954999579*alphaDrSurf[1]*fUpwind[4]+0.4472135954999579*fUpwind[1]*alphaDrSurf[4]+0.5*alphaDrSurf[0]*fUpwind[1]+0.5*fUpwind[0]*alphaDrSurf[1]; 
  Ghat[2] = Gdiff[2]*rdv2+0.5000000000000001*alphaDrSurf[4]*fUpwind[6]+0.5*alphaDrSurf[1]*fUpwind[3]+0.5*alphaDrSurf[0]*fUpwind[2]; 
  Ghat[3] = Gdiff[3]*rdv2+0.447213595499958*alphaDrSurf[1]*fUpwind[6]+0.4472135954999579*fUpwind[3]*alphaDrSurf[4]+0.5*alphaDrSurf[0]*fUpwind[3]+0.5*alphaDrSurf[1]*fUpwind[2]; 
  Ghat[4] = Gdiff[4]*rdv2+0.31943828249997*alphaDrSurf[4]*fUpwind[4]+0.5*alphaDrSurf[0]*fUpwind[4]+0.5*fUpwind[0]*alphaDrSurf[4]+0.4472135954999579*alphaDrSurf[1]*fUpwind[1]; 
  Ghat[5] = Gdiff[5]*rdv2+0.5000000000000001*alphaDrSurf[1]*fUpwind[7]+0.5*alphaDrSurf[0]*fUpwind[5]; 
  Ghat[6] = Gdiff[6]*rdv2+0.31943828249997*alphaDrSurf[4]*fUpwind[6]+0.5*alphaDrSurf[0]*fUpwind[6]+0.5000000000000001*fUpwind[2]*alphaDrSurf[4]+0.447213595499958*alphaDrSurf[1]*fUpwind[3]; 
  Ghat[7] = Gdiff[7]*rdv2+0.4472135954999579*alphaDrSurf[4]*fUpwind[7]+0.5*alphaDrSurf[0]*fUpwind[7]+0.5000000000000001*alphaDrSurf[1]*fUpwind[5]; 

  out[0] += 0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += 1.369306393762915*nuVtSqSum[1]*fSkin[15]*rdvSq4-1.060660171779821*nuVtSqSum[2]*fSkin[13]*rdvSq4+1.369306393762915*nuVtSqSum[0]*fSkin[9]*rdvSq4+0.6123724356957944*nuVtSqSum[2]*fSkin[7]*rdvSq4-1.060660171779821*nuVtSqSum[1]*fSkin[5]*rdvSq4-1.060660171779821*nuVtSqSum[0]*fSkin[3]*rdvSq4+0.6123724356957944*fSkin[1]*nuVtSqSum[1]*rdvSq4+0.6123724356957944*fSkin[0]*nuVtSqSum[0]*rdvSq4-1.224744871391589*Gdiff2[0]*rdvSq4+1.224744871391589*Ghat[0]*rdv2; 
  out[4] += 0.7071067811865475*Ghat[3]*rdv2; 
  out[5] += 1.224744871391589*nuVtSqSum[2]*fSkin[15]*rdvSq4+1.369306393762915*nuVtSqSum[0]*fSkin[15]*rdvSq4-0.9486832980505138*nuVtSqSum[1]*fSkin[13]*rdvSq4+1.369306393762915*nuVtSqSum[1]*fSkin[9]*rdvSq4+0.5477225575051661*nuVtSqSum[1]*fSkin[7]*rdvSq4-0.9486832980505137*nuVtSqSum[2]*fSkin[5]*rdvSq4-1.060660171779821*nuVtSqSum[0]*fSkin[5]*rdvSq4-1.060660171779821*nuVtSqSum[1]*fSkin[3]*rdvSq4+0.5477225575051661*fSkin[1]*nuVtSqSum[2]*rdvSq4+0.6123724356957944*fSkin[0]*nuVtSqSum[1]*rdvSq4+0.6123724356957944*nuVtSqSum[0]*fSkin[1]*rdvSq4-1.224744871391589*Gdiff2[1]*rdvSq4+1.224744871391589*Ghat[1]*rdv2; 
  out[6] += 1.369306393762915*nuVtSqSum[1]*fSkin[19]*rdvSq4-1.060660171779821*nuVtSqSum[2]*fSkin[17]*rdvSq4+1.369306393762915*nuVtSqSum[0]*fSkin[16]*rdvSq4+0.6123724356957944*nuVtSqSum[2]*fSkin[11]*rdvSq4-1.060660171779821*nuVtSqSum[1]*fSkin[10]*rdvSq4-1.060660171779821*nuVtSqSum[0]*fSkin[6]*rdvSq4+0.6123724356957944*nuVtSqSum[1]*fSkin[4]*rdvSq4+0.6123724356957944*nuVtSqSum[0]*fSkin[2]*rdvSq4-1.224744871391589*Gdiff2[2]*rdvSq4+1.224744871391589*Ghat[2]*rdv2; 
  out[7] += 0.7071067811865475*Ghat[4]*rdv2; 
  out[8] += 0.7071067811865475*Ghat[5]*rdv2; 
  out[9] += (-5.303300858899106*nuVtSqSum[1]*fSkin[15]*rdvSq4)+4.107919181288745*nuVtSqSum[2]*fSkin[13]*rdvSq4-5.303300858899105*nuVtSqSum[0]*fSkin[9]*rdvSq4-2.371708245126284*nuVtSqSum[2]*fSkin[7]*rdvSq4+4.107919181288745*nuVtSqSum[1]*fSkin[5]*rdvSq4+4.107919181288745*nuVtSqSum[0]*fSkin[3]*rdvSq4-2.371708245126284*fSkin[1]*nuVtSqSum[1]*rdvSq4-2.371708245126284*fSkin[0]*nuVtSqSum[0]*rdvSq4-4.743416490252569*Gdiff2[0]*rdvSq4+1.58113883008419*Ghat[0]*rdv2; 
  out[10] += 1.224744871391589*nuVtSqSum[2]*fSkin[19]*rdvSq4+1.369306393762915*nuVtSqSum[0]*fSkin[19]*rdvSq4-0.9486832980505137*nuVtSqSum[1]*fSkin[17]*rdvSq4+1.369306393762915*nuVtSqSum[1]*fSkin[16]*rdvSq4+0.5477225575051661*nuVtSqSum[1]*fSkin[11]*rdvSq4-0.9486832980505137*nuVtSqSum[2]*fSkin[10]*rdvSq4-1.060660171779821*nuVtSqSum[0]*fSkin[10]*rdvSq4-1.060660171779821*nuVtSqSum[1]*fSkin[6]*rdvSq4+0.5477225575051661*nuVtSqSum[2]*fSkin[4]*rdvSq4+0.6123724356957944*nuVtSqSum[0]*fSkin[4]*rdvSq4-1.224744871391589*Gdiff2[3]*rdvSq4+0.6123724356957944*nuVtSqSum[1]*fSkin[2]*rdvSq4+1.224744871391589*Ghat[3]*rdv2; 
  out[11] += 0.7071067811865475*Ghat[6]*rdv2; 
  out[12] += 0.7071067811865475*Ghat[7]*rdv2; 
  out[13] += 1.224744871391589*nuVtSqSum[1]*fSkin[15]*rdvSq4-0.6776309271789384*nuVtSqSum[2]*fSkin[13]*rdvSq4-1.060660171779821*nuVtSqSum[0]*fSkin[13]*rdvSq4+1.369306393762915*nuVtSqSum[2]*fSkin[9]*rdvSq4+0.3912303982179757*nuVtSqSum[2]*fSkin[7]*rdvSq4+0.6123724356957944*nuVtSqSum[0]*fSkin[7]*rdvSq4-0.9486832980505138*nuVtSqSum[1]*fSkin[5]*rdvSq4-1.224744871391589*Gdiff2[4]*rdvSq4-1.060660171779821*nuVtSqSum[2]*fSkin[3]*rdvSq4+0.6123724356957944*fSkin[0]*nuVtSqSum[2]*rdvSq4+0.5477225575051661*fSkin[1]*nuVtSqSum[1]*rdvSq4+1.224744871391589*Ghat[4]*rdv2; 
  out[14] += (-1.060660171779821*nuVtSqSum[1]*fSkin[18]*rdvSq4)-1.060660171779821*nuVtSqSum[0]*fSkin[14]*rdvSq4+0.6123724356957944*nuVtSqSum[1]*fSkin[12]*rdvSq4+0.6123724356957944*nuVtSqSum[0]*fSkin[8]*rdvSq4-1.224744871391589*Gdiff2[5]*rdvSq4+1.224744871391589*Ghat[5]*rdv2; 
  out[15] += (-4.743416490252569*nuVtSqSum[2]*fSkin[15]*rdvSq4)-5.303300858899105*nuVtSqSum[0]*fSkin[15]*rdvSq4+3.674234614174766*nuVtSqSum[1]*fSkin[13]*rdvSq4-5.303300858899106*nuVtSqSum[1]*fSkin[9]*rdvSq4-2.121320343559642*nuVtSqSum[1]*fSkin[7]*rdvSq4+3.674234614174767*nuVtSqSum[2]*fSkin[5]*rdvSq4+4.107919181288746*nuVtSqSum[0]*fSkin[5]*rdvSq4+4.107919181288746*nuVtSqSum[1]*fSkin[3]*rdvSq4-2.121320343559642*fSkin[1]*nuVtSqSum[2]*rdvSq4-2.371708245126284*fSkin[0]*nuVtSqSum[1]*rdvSq4-2.371708245126284*nuVtSqSum[0]*fSkin[1]*rdvSq4-4.743416490252569*Gdiff2[1]*rdvSq4+1.58113883008419*Ghat[1]*rdv2; 
  out[16] += (-5.303300858899106*nuVtSqSum[1]*fSkin[19]*rdvSq4)+4.107919181288745*nuVtSqSum[2]*fSkin[17]*rdvSq4-5.303300858899105*nuVtSqSum[0]*fSkin[16]*rdvSq4-2.371708245126284*nuVtSqSum[2]*fSkin[11]*rdvSq4+4.107919181288745*nuVtSqSum[1]*fSkin[10]*rdvSq4+4.107919181288745*nuVtSqSum[0]*fSkin[6]*rdvSq4-2.371708245126284*nuVtSqSum[1]*fSkin[4]*rdvSq4-2.371708245126284*nuVtSqSum[0]*fSkin[2]*rdvSq4-4.743416490252569*Gdiff2[2]*rdvSq4+1.58113883008419*Ghat[2]*rdv2; 
  out[17] += 1.224744871391589*nuVtSqSum[1]*fSkin[19]*rdvSq4-0.6776309271789384*nuVtSqSum[2]*fSkin[17]*rdvSq4-1.060660171779821*nuVtSqSum[0]*fSkin[17]*rdvSq4+1.369306393762915*nuVtSqSum[2]*fSkin[16]*rdvSq4+0.3912303982179757*nuVtSqSum[2]*fSkin[11]*rdvSq4+0.6123724356957944*nuVtSqSum[0]*fSkin[11]*rdvSq4-0.9486832980505137*nuVtSqSum[1]*fSkin[10]*rdvSq4-1.060660171779821*nuVtSqSum[2]*fSkin[6]*rdvSq4-1.224744871391589*Gdiff2[6]*rdvSq4+0.5477225575051661*nuVtSqSum[1]*fSkin[4]*rdvSq4+0.6123724356957944*fSkin[2]*nuVtSqSum[2]*rdvSq4+1.224744871391589*Ghat[6]*rdv2; 
  out[18] += (-0.9486832980505137*nuVtSqSum[2]*fSkin[18]*rdvSq4)-1.060660171779821*nuVtSqSum[0]*fSkin[18]*rdvSq4-1.060660171779821*nuVtSqSum[1]*fSkin[14]*rdvSq4+0.5477225575051661*nuVtSqSum[2]*fSkin[12]*rdvSq4+0.6123724356957944*nuVtSqSum[0]*fSkin[12]*rdvSq4+0.6123724356957944*nuVtSqSum[1]*fSkin[8]*rdvSq4-1.224744871391589*Gdiff2[7]*rdvSq4+1.224744871391589*Ghat[7]*rdv2; 
  out[19] += (-4.743416490252569*nuVtSqSum[2]*fSkin[19]*rdvSq4)-5.303300858899105*nuVtSqSum[0]*fSkin[19]*rdvSq4+3.674234614174766*nuVtSqSum[1]*fSkin[17]*rdvSq4-5.303300858899106*nuVtSqSum[1]*fSkin[16]*rdvSq4-2.121320343559642*nuVtSqSum[1]*fSkin[11]*rdvSq4+3.674234614174766*nuVtSqSum[2]*fSkin[10]*rdvSq4+4.107919181288745*nuVtSqSum[0]*fSkin[10]*rdvSq4+4.107919181288745*nuVtSqSum[1]*fSkin[6]*rdvSq4-2.121320343559642*nuVtSqSum[2]*fSkin[4]*rdvSq4-2.371708245126284*nuVtSqSum[0]*fSkin[4]*rdvSq4-4.743416490252569*Gdiff2[3]*rdvSq4-2.371708245126284*nuVtSqSum[1]*fSkin[2]*rdvSq4+1.58113883008419*Ghat[3]*rdv2; 

  } else { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*(2.828427124746191*w[2]-1.414213562373095*dxv[2])-2.828427124746191*sumNuUy[0]); 
  alphaDrSurf[1] = 0.5*(nuSum[1]*(2.828427124746191*w[2]-1.414213562373095*dxv[2])-2.828427124746191*sumNuUy[1]); 
  alphaDrSurf[4] = 0.5*(2.828427124746191*nuSum[2]*w[2]-2.828427124746191*sumNuUy[2]-1.414213562373095*dxv[2]*nuSum[2]); 

  if (0.4472135954999579*alphaDrSurf[4]-0.6708203932499369*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = ser_1x2v_p2_surfvy_quad_0(1, fEdge); 
  } else { 
    fUpwindQuad[0] = ser_1x2v_p2_surfvy_quad_0(-1, fSkin); 
  } 
  if (0.5*alphaDrSurf[0]-0.5590169943749475*alphaDrSurf[4] < 0) { 
    fUpwindQuad[1] = ser_1x2v_p2_surfvy_quad_1(1, fEdge); 
  } else { 
    fUpwindQuad[1] = ser_1x2v_p2_surfvy_quad_1(-1, fSkin); 
  } 
  if (0.4472135954999579*alphaDrSurf[4]+0.6708203932499369*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[2] = ser_1x2v_p2_surfvy_quad_2(1, fEdge); 
  } else { 
    fUpwindQuad[2] = ser_1x2v_p2_surfvy_quad_2(-1, fSkin); 
  } 
  if (0.4472135954999579*alphaDrSurf[4]-0.6708203932499369*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = ser_1x2v_p2_surfvy_quad_3(1, fEdge); 
  } else { 
    fUpwindQuad[3] = ser_1x2v_p2_surfvy_quad_3(-1, fSkin); 
  } 
  if (0.5*alphaDrSurf[0]-0.5590169943749475*alphaDrSurf[4] < 0) { 
    fUpwindQuad[4] = ser_1x2v_p2_surfvy_quad_4(1, fEdge); 
  } else { 
    fUpwindQuad[4] = ser_1x2v_p2_surfvy_quad_4(-1, fSkin); 
  } 
  if (0.4472135954999579*alphaDrSurf[4]+0.6708203932499369*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[5] = ser_1x2v_p2_surfvy_quad_5(1, fEdge); 
  } else { 
    fUpwindQuad[5] = ser_1x2v_p2_surfvy_quad_5(-1, fSkin); 
  } 
  if (0.4472135954999579*alphaDrSurf[4]-0.6708203932499369*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[6] = ser_1x2v_p2_surfvy_quad_6(1, fEdge); 
  } else { 
    fUpwindQuad[6] = ser_1x2v_p2_surfvy_quad_6(-1, fSkin); 
  } 
  if (0.5*alphaDrSurf[0]-0.5590169943749475*alphaDrSurf[4] < 0) { 
    fUpwindQuad[7] = ser_1x2v_p2_surfvy_quad_7(1, fEdge); 
  } else { 
    fUpwindQuad[7] = ser_1x2v_p2_surfvy_quad_7(-1, fSkin); 
  } 
  if (0.4472135954999579*alphaDrSurf[4]+0.6708203932499369*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[8] = ser_1x2v_p2_surfvy_quad_8(1, fEdge); 
  } else { 
    fUpwindQuad[8] = ser_1x2v_p2_surfvy_quad_8(-1, fSkin); 
  } 

  fUpwind[0] = 0.154320987654321*fUpwindQuad[8]+0.2469135802469136*fUpwindQuad[7]+0.154320987654321*fUpwindQuad[6]+0.2469135802469136*fUpwindQuad[5]+0.3950617283950617*fUpwindQuad[4]+0.2469135802469136*fUpwindQuad[3]+0.154320987654321*fUpwindQuad[2]+0.2469135802469136*fUpwindQuad[1]+0.154320987654321*fUpwindQuad[0]; 
  fUpwind[1] = 0.2070433312499806*fUpwindQuad[8]-0.2070433312499806*fUpwindQuad[6]+0.3312693299999688*fUpwindQuad[5]-0.3312693299999688*fUpwindQuad[3]+0.2070433312499806*fUpwindQuad[2]-0.2070433312499806*fUpwindQuad[0]; 
  fUpwind[2] = 0.2070433312499806*fUpwindQuad[8]+0.3312693299999688*fUpwindQuad[7]+0.2070433312499806*fUpwindQuad[6]-0.2070433312499806*fUpwindQuad[2]-0.3312693299999688*fUpwindQuad[1]-0.2070433312499806*fUpwindQuad[0]; 
  fUpwind[3] = 0.2777777777777778*fUpwindQuad[8]-0.2777777777777778*fUpwindQuad[6]-0.2777777777777778*fUpwindQuad[2]+0.2777777777777778*fUpwindQuad[0]; 
  fUpwind[4] = 0.138028887499987*fUpwindQuad[8]-0.2760577749999741*fUpwindQuad[7]+0.138028887499987*fUpwindQuad[6]+0.2208462199999792*fUpwindQuad[5]-0.4416924399999584*fUpwindQuad[4]+0.2208462199999792*fUpwindQuad[3]+0.138028887499987*fUpwindQuad[2]-0.2760577749999741*fUpwindQuad[1]+0.138028887499987*fUpwindQuad[0]; 
  fUpwind[5] = 0.138028887499987*fUpwindQuad[8]+0.2208462199999792*fUpwindQuad[7]+0.138028887499987*fUpwindQuad[6]-0.2760577749999741*fUpwindQuad[5]-0.4416924399999584*fUpwindQuad[4]-0.2760577749999741*fUpwindQuad[3]+0.138028887499987*fUpwindQuad[2]+0.2208462199999792*fUpwindQuad[1]+0.138028887499987*fUpwindQuad[0]; 
  fUpwind[6] = 0.1851851851851853*fUpwindQuad[8]-0.3703703703703705*fUpwindQuad[7]+0.1851851851851853*fUpwindQuad[6]-0.1851851851851853*fUpwindQuad[2]+0.3703703703703705*fUpwindQuad[1]-0.1851851851851853*fUpwindQuad[0]; 
  fUpwind[7] = 0.1851851851851853*fUpwindQuad[8]-0.1851851851851853*fUpwindQuad[6]-0.3703703703703705*fUpwindQuad[5]+0.3703703703703705*fUpwindQuad[3]+0.1851851851851853*fUpwindQuad[2]-0.1851851851851853*fUpwindQuad[0]; 

  Gdiff2[0] = 0.2445699350390395*nuVtSqSum[1]*fSkin[15]+0.2445699350390395*nuVtSqSum[1]*fEdge[15]-0.3518228202874282*nuVtSqSum[2]*fSkin[13]+0.3518228202874282*nuVtSqSum[2]*fEdge[13]+0.2445699350390395*nuVtSqSum[0]*fSkin[9]+0.2445699350390395*nuVtSqSum[0]*fEdge[9]+0.25*nuVtSqSum[2]*fSkin[7]+0.25*nuVtSqSum[2]*fEdge[7]-0.3518228202874282*nuVtSqSum[1]*fSkin[5]+0.3518228202874282*nuVtSqSum[1]*fEdge[5]-0.3518228202874282*nuVtSqSum[0]*fSkin[3]+0.3518228202874282*nuVtSqSum[0]*fEdge[3]+0.25*fSkin[1]*nuVtSqSum[1]+0.25*fEdge[1]*nuVtSqSum[1]+0.25*fSkin[0]*nuVtSqSum[0]+0.25*fEdge[0]*nuVtSqSum[0]; 
  Gdiff2[1] = 0.21875*nuVtSqSum[2]*fSkin[15]+0.2445699350390395*nuVtSqSum[0]*fSkin[15]+0.21875*nuVtSqSum[2]*fEdge[15]+0.2445699350390395*nuVtSqSum[0]*fEdge[15]-0.3146798968793527*nuVtSqSum[1]*fSkin[13]+0.3146798968793527*nuVtSqSum[1]*fEdge[13]+0.2445699350390395*nuVtSqSum[1]*fSkin[9]+0.2445699350390395*nuVtSqSum[1]*fEdge[9]+0.223606797749979*nuVtSqSum[1]*fSkin[7]+0.223606797749979*nuVtSqSum[1]*fEdge[7]-0.3146798968793526*nuVtSqSum[2]*fSkin[5]-0.3518228202874282*nuVtSqSum[0]*fSkin[5]+0.3146798968793526*nuVtSqSum[2]*fEdge[5]+0.3518228202874282*nuVtSqSum[0]*fEdge[5]-0.3518228202874282*nuVtSqSum[1]*fSkin[3]+0.3518228202874282*nuVtSqSum[1]*fEdge[3]+0.223606797749979*fSkin[1]*nuVtSqSum[2]+0.223606797749979*fEdge[1]*nuVtSqSum[2]+0.25*fSkin[0]*nuVtSqSum[1]+0.25*fEdge[0]*nuVtSqSum[1]+0.25*nuVtSqSum[0]*fSkin[1]+0.25*nuVtSqSum[0]*fEdge[1]; 
  Gdiff2[2] = 0.2445699350390395*nuVtSqSum[1]*fSkin[19]+0.2445699350390395*nuVtSqSum[1]*fEdge[19]-0.3518228202874282*nuVtSqSum[2]*fSkin[17]+0.3518228202874282*nuVtSqSum[2]*fEdge[17]+0.2445699350390395*nuVtSqSum[0]*fSkin[16]+0.2445699350390395*nuVtSqSum[0]*fEdge[16]+0.2500000000000001*nuVtSqSum[2]*fSkin[11]+0.2500000000000001*nuVtSqSum[2]*fEdge[11]-0.3518228202874282*nuVtSqSum[1]*fSkin[10]+0.3518228202874282*nuVtSqSum[1]*fEdge[10]-0.3518228202874282*nuVtSqSum[0]*fSkin[6]+0.3518228202874282*nuVtSqSum[0]*fEdge[6]+0.25*nuVtSqSum[1]*fSkin[4]+0.25*nuVtSqSum[1]*fEdge[4]+0.25*nuVtSqSum[0]*fSkin[2]+0.25*nuVtSqSum[0]*fEdge[2]; 
  Gdiff2[3] = 0.21875*nuVtSqSum[2]*fSkin[19]+0.2445699350390395*nuVtSqSum[0]*fSkin[19]+0.21875*nuVtSqSum[2]*fEdge[19]+0.2445699350390395*nuVtSqSum[0]*fEdge[19]-0.3146798968793526*nuVtSqSum[1]*fSkin[17]+0.3146798968793526*nuVtSqSum[1]*fEdge[17]+0.2445699350390395*nuVtSqSum[1]*fSkin[16]+0.2445699350390395*nuVtSqSum[1]*fEdge[16]+0.223606797749979*nuVtSqSum[1]*fSkin[11]+0.223606797749979*nuVtSqSum[1]*fEdge[11]-0.3146798968793526*nuVtSqSum[2]*fSkin[10]-0.3518228202874282*nuVtSqSum[0]*fSkin[10]+0.3146798968793526*nuVtSqSum[2]*fEdge[10]+0.3518228202874282*nuVtSqSum[0]*fEdge[10]-0.3518228202874282*nuVtSqSum[1]*fSkin[6]+0.3518228202874282*nuVtSqSum[1]*fEdge[6]+0.223606797749979*nuVtSqSum[2]*fSkin[4]+0.25*nuVtSqSum[0]*fSkin[4]+0.223606797749979*nuVtSqSum[2]*fEdge[4]+0.25*nuVtSqSum[0]*fEdge[4]+0.25*nuVtSqSum[1]*fSkin[2]+0.25*nuVtSqSum[1]*fEdge[2]; 
  Gdiff2[4] = 0.21875*nuVtSqSum[1]*fSkin[15]+0.21875*nuVtSqSum[1]*fEdge[15]-0.2247713549138233*nuVtSqSum[2]*fSkin[13]-0.3518228202874282*nuVtSqSum[0]*fSkin[13]+0.2247713549138233*nuVtSqSum[2]*fEdge[13]+0.3518228202874282*nuVtSqSum[0]*fEdge[13]+0.2445699350390395*nuVtSqSum[2]*fSkin[9]+0.2445699350390395*nuVtSqSum[2]*fEdge[9]+0.159719141249985*nuVtSqSum[2]*fSkin[7]+0.25*nuVtSqSum[0]*fSkin[7]+0.159719141249985*nuVtSqSum[2]*fEdge[7]+0.25*nuVtSqSum[0]*fEdge[7]-0.3146798968793526*nuVtSqSum[1]*fSkin[5]+0.3146798968793526*nuVtSqSum[1]*fEdge[5]-0.3518228202874282*nuVtSqSum[2]*fSkin[3]+0.3518228202874282*nuVtSqSum[2]*fEdge[3]+0.25*fSkin[0]*nuVtSqSum[2]+0.25*fEdge[0]*nuVtSqSum[2]+0.223606797749979*fSkin[1]*nuVtSqSum[1]+0.223606797749979*fEdge[1]*nuVtSqSum[1]; 
  Gdiff2[5] = (-0.3518228202874282*nuVtSqSum[1]*fSkin[18])+0.3518228202874282*nuVtSqSum[1]*fEdge[18]-0.3518228202874282*nuVtSqSum[0]*fSkin[14]+0.3518228202874282*nuVtSqSum[0]*fEdge[14]+0.2500000000000001*nuVtSqSum[1]*fSkin[12]+0.2500000000000001*nuVtSqSum[1]*fEdge[12]+0.25*nuVtSqSum[0]*fSkin[8]+0.25*nuVtSqSum[0]*fEdge[8]; 
  Gdiff2[6] = 0.21875*nuVtSqSum[1]*fSkin[19]+0.21875*nuVtSqSum[1]*fEdge[19]-0.2247713549138233*nuVtSqSum[2]*fSkin[17]-0.3518228202874282*nuVtSqSum[0]*fSkin[17]+0.2247713549138233*nuVtSqSum[2]*fEdge[17]+0.3518228202874282*nuVtSqSum[0]*fEdge[17]+0.2445699350390395*nuVtSqSum[2]*fSkin[16]+0.2445699350390395*nuVtSqSum[2]*fEdge[16]+0.159719141249985*nuVtSqSum[2]*fSkin[11]+0.25*nuVtSqSum[0]*fSkin[11]+0.159719141249985*nuVtSqSum[2]*fEdge[11]+0.25*nuVtSqSum[0]*fEdge[11]-0.3146798968793527*nuVtSqSum[1]*fSkin[10]+0.3146798968793527*nuVtSqSum[1]*fEdge[10]-0.3518228202874282*nuVtSqSum[2]*fSkin[6]+0.3518228202874282*nuVtSqSum[2]*fEdge[6]+0.223606797749979*nuVtSqSum[1]*fSkin[4]+0.223606797749979*nuVtSqSum[1]*fEdge[4]+0.2500000000000001*fSkin[2]*nuVtSqSum[2]+0.2500000000000001*fEdge[2]*nuVtSqSum[2]; 
  Gdiff2[7] = (-0.3146798968793527*nuVtSqSum[2]*fSkin[18])-0.3518228202874282*nuVtSqSum[0]*fSkin[18]+0.3146798968793527*nuVtSqSum[2]*fEdge[18]+0.3518228202874282*nuVtSqSum[0]*fEdge[18]-0.3518228202874282*nuVtSqSum[1]*fSkin[14]+0.3518228202874282*nuVtSqSum[1]*fEdge[14]+0.223606797749979*nuVtSqSum[2]*fSkin[12]+0.25*nuVtSqSum[0]*fSkin[12]+0.223606797749979*nuVtSqSum[2]*fEdge[12]+0.25*nuVtSqSum[0]*fEdge[12]+0.2500000000000001*nuVtSqSum[1]*fSkin[8]+0.2500000000000001*nuVtSqSum[1]*fEdge[8]; 

  Gdiff[0] = 0.6708203932499369*nuVtSqSum[1]*fSkin[15]-0.6708203932499369*nuVtSqSum[1]*fEdge[15]-1.190784930203603*nuVtSqSum[2]*fSkin[13]-1.190784930203603*nuVtSqSum[2]*fEdge[13]+0.6708203932499369*nuVtSqSum[0]*fSkin[9]-0.6708203932499369*nuVtSqSum[0]*fEdge[9]+0.9375*nuVtSqSum[2]*fSkin[7]-0.9375*nuVtSqSum[2]*fEdge[7]-1.190784930203603*nuVtSqSum[1]*fSkin[5]-1.190784930203603*nuVtSqSum[1]*fEdge[5]-1.190784930203603*nuVtSqSum[0]*fSkin[3]-1.190784930203603*nuVtSqSum[0]*fEdge[3]+0.9375*fSkin[1]*nuVtSqSum[1]-0.9375*fEdge[1]*nuVtSqSum[1]+0.9375*fSkin[0]*nuVtSqSum[0]-0.9375*fEdge[0]*nuVtSqSum[0]; 
  Gdiff[1] = 0.5999999999999999*nuVtSqSum[2]*fSkin[15]+0.6708203932499369*nuVtSqSum[0]*fSkin[15]-0.5999999999999999*nuVtSqSum[2]*fEdge[15]-0.6708203932499369*nuVtSqSum[0]*fEdge[15]-1.06507042020704*nuVtSqSum[1]*fSkin[13]-1.06507042020704*nuVtSqSum[1]*fEdge[13]+0.6708203932499369*nuVtSqSum[1]*fSkin[9]-0.6708203932499369*nuVtSqSum[1]*fEdge[9]+0.8385254915624212*nuVtSqSum[1]*fSkin[7]-0.8385254915624212*nuVtSqSum[1]*fEdge[7]-1.06507042020704*nuVtSqSum[2]*fSkin[5]-1.190784930203603*nuVtSqSum[0]*fSkin[5]-1.06507042020704*nuVtSqSum[2]*fEdge[5]-1.190784930203603*nuVtSqSum[0]*fEdge[5]-1.190784930203603*nuVtSqSum[1]*fSkin[3]-1.190784930203603*nuVtSqSum[1]*fEdge[3]+0.8385254915624212*fSkin[1]*nuVtSqSum[2]-0.8385254915624212*fEdge[1]*nuVtSqSum[2]+0.9375*fSkin[0]*nuVtSqSum[1]-0.9375*fEdge[0]*nuVtSqSum[1]+0.9375*nuVtSqSum[0]*fSkin[1]-0.9375*nuVtSqSum[0]*fEdge[1]; 
  Gdiff[2] = 0.6708203932499369*nuVtSqSum[1]*fSkin[19]-0.6708203932499369*nuVtSqSum[1]*fEdge[19]-1.190784930203603*nuVtSqSum[2]*fSkin[17]-1.190784930203603*nuVtSqSum[2]*fEdge[17]+0.6708203932499369*nuVtSqSum[0]*fSkin[16]-0.6708203932499369*nuVtSqSum[0]*fEdge[16]+0.9375000000000001*nuVtSqSum[2]*fSkin[11]-0.9375000000000001*nuVtSqSum[2]*fEdge[11]-1.190784930203603*nuVtSqSum[1]*fSkin[10]-1.190784930203603*nuVtSqSum[1]*fEdge[10]-1.190784930203603*nuVtSqSum[0]*fSkin[6]-1.190784930203603*nuVtSqSum[0]*fEdge[6]+0.9375*nuVtSqSum[1]*fSkin[4]-0.9375*nuVtSqSum[1]*fEdge[4]+0.9375*nuVtSqSum[0]*fSkin[2]-0.9375*nuVtSqSum[0]*fEdge[2]; 
  Gdiff[3] = 0.6*nuVtSqSum[2]*fSkin[19]+0.6708203932499369*nuVtSqSum[0]*fSkin[19]-0.6*nuVtSqSum[2]*fEdge[19]-0.6708203932499369*nuVtSqSum[0]*fEdge[19]-1.06507042020704*nuVtSqSum[1]*fSkin[17]-1.06507042020704*nuVtSqSum[1]*fEdge[17]+0.6708203932499369*nuVtSqSum[1]*fSkin[16]-0.6708203932499369*nuVtSqSum[1]*fEdge[16]+0.8385254915624211*nuVtSqSum[1]*fSkin[11]-0.8385254915624211*nuVtSqSum[1]*fEdge[11]-1.06507042020704*nuVtSqSum[2]*fSkin[10]-1.190784930203603*nuVtSqSum[0]*fSkin[10]-1.06507042020704*nuVtSqSum[2]*fEdge[10]-1.190784930203603*nuVtSqSum[0]*fEdge[10]-1.190784930203603*nuVtSqSum[1]*fSkin[6]-1.190784930203603*nuVtSqSum[1]*fEdge[6]+0.8385254915624212*nuVtSqSum[2]*fSkin[4]+0.9375*nuVtSqSum[0]*fSkin[4]-0.8385254915624212*nuVtSqSum[2]*fEdge[4]-0.9375*nuVtSqSum[0]*fEdge[4]+0.9375*nuVtSqSum[1]*fSkin[2]-0.9375*nuVtSqSum[1]*fEdge[2]; 
  Gdiff[4] = 0.5999999999999999*nuVtSqSum[1]*fSkin[15]-0.5999999999999999*nuVtSqSum[1]*fEdge[15]-0.7607645858621712*nuVtSqSum[2]*fSkin[13]-1.190784930203603*nuVtSqSum[0]*fSkin[13]-0.7607645858621712*nuVtSqSum[2]*fEdge[13]-1.190784930203603*nuVtSqSum[0]*fEdge[13]+0.6708203932499369*nuVtSqSum[2]*fSkin[9]-0.6708203932499369*nuVtSqSum[2]*fEdge[9]+0.5989467796874438*nuVtSqSum[2]*fSkin[7]+0.9375*nuVtSqSum[0]*fSkin[7]-0.5989467796874438*nuVtSqSum[2]*fEdge[7]-0.9375*nuVtSqSum[0]*fEdge[7]-1.06507042020704*nuVtSqSum[1]*fSkin[5]-1.06507042020704*nuVtSqSum[1]*fEdge[5]-1.190784930203603*nuVtSqSum[2]*fSkin[3]-1.190784930203603*nuVtSqSum[2]*fEdge[3]+0.9375*fSkin[0]*nuVtSqSum[2]-0.9375*fEdge[0]*nuVtSqSum[2]+0.8385254915624212*fSkin[1]*nuVtSqSum[1]-0.8385254915624212*fEdge[1]*nuVtSqSum[1]; 
  Gdiff[5] = (-1.190784930203603*nuVtSqSum[1]*fSkin[18])-1.190784930203603*nuVtSqSum[1]*fEdge[18]-1.190784930203603*nuVtSqSum[0]*fSkin[14]-1.190784930203603*nuVtSqSum[0]*fEdge[14]+0.9375000000000001*nuVtSqSum[1]*fSkin[12]-0.9375000000000001*nuVtSqSum[1]*fEdge[12]+0.9375*nuVtSqSum[0]*fSkin[8]-0.9375*nuVtSqSum[0]*fEdge[8]; 
  Gdiff[6] = 0.5999999999999999*nuVtSqSum[1]*fSkin[19]-0.5999999999999999*nuVtSqSum[1]*fEdge[19]-0.7607645858621712*nuVtSqSum[2]*fSkin[17]-1.190784930203603*nuVtSqSum[0]*fSkin[17]-0.7607645858621712*nuVtSqSum[2]*fEdge[17]-1.190784930203603*nuVtSqSum[0]*fEdge[17]+0.6708203932499369*nuVtSqSum[2]*fSkin[16]-0.6708203932499369*nuVtSqSum[2]*fEdge[16]+0.5989467796874438*nuVtSqSum[2]*fSkin[11]+0.9375*nuVtSqSum[0]*fSkin[11]-0.5989467796874438*nuVtSqSum[2]*fEdge[11]-0.9375*nuVtSqSum[0]*fEdge[11]-1.06507042020704*nuVtSqSum[1]*fSkin[10]-1.06507042020704*nuVtSqSum[1]*fEdge[10]-1.190784930203603*nuVtSqSum[2]*fSkin[6]-1.190784930203603*nuVtSqSum[2]*fEdge[6]+0.8385254915624211*nuVtSqSum[1]*fSkin[4]-0.8385254915624211*nuVtSqSum[1]*fEdge[4]+0.9375000000000001*fSkin[2]*nuVtSqSum[2]-0.9375000000000001*fEdge[2]*nuVtSqSum[2]; 
  Gdiff[7] = (-1.06507042020704*nuVtSqSum[2]*fSkin[18])-1.190784930203603*nuVtSqSum[0]*fSkin[18]-1.06507042020704*nuVtSqSum[2]*fEdge[18]-1.190784930203603*nuVtSqSum[0]*fEdge[18]-1.190784930203603*nuVtSqSum[1]*fSkin[14]-1.190784930203603*nuVtSqSum[1]*fEdge[14]+0.8385254915624212*nuVtSqSum[2]*fSkin[12]+0.9375*nuVtSqSum[0]*fSkin[12]-0.8385254915624212*nuVtSqSum[2]*fEdge[12]-0.9375*nuVtSqSum[0]*fEdge[12]+0.9375000000000001*nuVtSqSum[1]*fSkin[8]-0.9375000000000001*nuVtSqSum[1]*fEdge[8]; 

  Ghat[0] = Gdiff[0]*rdv2+0.5*alphaDrSurf[4]*fUpwind[4]+0.5*alphaDrSurf[1]*fUpwind[1]+0.5*alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = Gdiff[1]*rdv2+0.4472135954999579*alphaDrSurf[1]*fUpwind[4]+0.4472135954999579*fUpwind[1]*alphaDrSurf[4]+0.5*alphaDrSurf[0]*fUpwind[1]+0.5*fUpwind[0]*alphaDrSurf[1]; 
  Ghat[2] = Gdiff[2]*rdv2+0.5000000000000001*alphaDrSurf[4]*fUpwind[6]+0.5*alphaDrSurf[1]*fUpwind[3]+0.5*alphaDrSurf[0]*fUpwind[2]; 
  Ghat[3] = Gdiff[3]*rdv2+0.447213595499958*alphaDrSurf[1]*fUpwind[6]+0.4472135954999579*fUpwind[3]*alphaDrSurf[4]+0.5*alphaDrSurf[0]*fUpwind[3]+0.5*alphaDrSurf[1]*fUpwind[2]; 
  Ghat[4] = Gdiff[4]*rdv2+0.31943828249997*alphaDrSurf[4]*fUpwind[4]+0.5*alphaDrSurf[0]*fUpwind[4]+0.5*fUpwind[0]*alphaDrSurf[4]+0.4472135954999579*alphaDrSurf[1]*fUpwind[1]; 
  Ghat[5] = Gdiff[5]*rdv2+0.5000000000000001*alphaDrSurf[1]*fUpwind[7]+0.5*alphaDrSurf[0]*fUpwind[5]; 
  Ghat[6] = Gdiff[6]*rdv2+0.31943828249997*alphaDrSurf[4]*fUpwind[6]+0.5*alphaDrSurf[0]*fUpwind[6]+0.5000000000000001*fUpwind[2]*alphaDrSurf[4]+0.447213595499958*alphaDrSurf[1]*fUpwind[3]; 
  Ghat[7] = Gdiff[7]*rdv2+0.4472135954999579*alphaDrSurf[4]*fUpwind[7]+0.5*alphaDrSurf[0]*fUpwind[7]+0.5000000000000001*alphaDrSurf[1]*fUpwind[5]; 

  out[0] += -0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += -0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += -0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += (-1.369306393762915*nuVtSqSum[1]*fSkin[15]*rdvSq4)-1.060660171779821*nuVtSqSum[2]*fSkin[13]*rdvSq4-1.369306393762915*nuVtSqSum[0]*fSkin[9]*rdvSq4-0.6123724356957944*nuVtSqSum[2]*fSkin[7]*rdvSq4-1.060660171779821*nuVtSqSum[1]*fSkin[5]*rdvSq4-1.060660171779821*nuVtSqSum[0]*fSkin[3]*rdvSq4-0.6123724356957944*fSkin[1]*nuVtSqSum[1]*rdvSq4-0.6123724356957944*fSkin[0]*nuVtSqSum[0]*rdvSq4+1.224744871391589*Gdiff2[0]*rdvSq4+1.224744871391589*Ghat[0]*rdv2; 
  out[4] += -0.7071067811865475*Ghat[3]*rdv2; 
  out[5] += (-1.224744871391589*nuVtSqSum[2]*fSkin[15]*rdvSq4)-1.369306393762915*nuVtSqSum[0]*fSkin[15]*rdvSq4-0.9486832980505138*nuVtSqSum[1]*fSkin[13]*rdvSq4-1.369306393762915*nuVtSqSum[1]*fSkin[9]*rdvSq4-0.5477225575051661*nuVtSqSum[1]*fSkin[7]*rdvSq4-0.9486832980505137*nuVtSqSum[2]*fSkin[5]*rdvSq4-1.060660171779821*nuVtSqSum[0]*fSkin[5]*rdvSq4-1.060660171779821*nuVtSqSum[1]*fSkin[3]*rdvSq4-0.5477225575051661*fSkin[1]*nuVtSqSum[2]*rdvSq4-0.6123724356957944*fSkin[0]*nuVtSqSum[1]*rdvSq4-0.6123724356957944*nuVtSqSum[0]*fSkin[1]*rdvSq4+1.224744871391589*Gdiff2[1]*rdvSq4+1.224744871391589*Ghat[1]*rdv2; 
  out[6] += (-1.369306393762915*nuVtSqSum[1]*fSkin[19]*rdvSq4)-1.060660171779821*nuVtSqSum[2]*fSkin[17]*rdvSq4-1.369306393762915*nuVtSqSum[0]*fSkin[16]*rdvSq4-0.6123724356957944*nuVtSqSum[2]*fSkin[11]*rdvSq4-1.060660171779821*nuVtSqSum[1]*fSkin[10]*rdvSq4-1.060660171779821*nuVtSqSum[0]*fSkin[6]*rdvSq4-0.6123724356957944*nuVtSqSum[1]*fSkin[4]*rdvSq4-0.6123724356957944*nuVtSqSum[0]*fSkin[2]*rdvSq4+1.224744871391589*Gdiff2[2]*rdvSq4+1.224744871391589*Ghat[2]*rdv2; 
  out[7] += -0.7071067811865475*Ghat[4]*rdv2; 
  out[8] += -0.7071067811865475*Ghat[5]*rdv2; 
  out[9] += (-5.303300858899106*nuVtSqSum[1]*fSkin[15]*rdvSq4)-4.107919181288745*nuVtSqSum[2]*fSkin[13]*rdvSq4-5.303300858899105*nuVtSqSum[0]*fSkin[9]*rdvSq4-2.371708245126284*nuVtSqSum[2]*fSkin[7]*rdvSq4-4.107919181288745*nuVtSqSum[1]*fSkin[5]*rdvSq4-4.107919181288745*nuVtSqSum[0]*fSkin[3]*rdvSq4-2.371708245126284*fSkin[1]*nuVtSqSum[1]*rdvSq4-2.371708245126284*fSkin[0]*nuVtSqSum[0]*rdvSq4-4.743416490252569*Gdiff2[0]*rdvSq4-1.58113883008419*Ghat[0]*rdv2; 
  out[10] += (-1.224744871391589*nuVtSqSum[2]*fSkin[19]*rdvSq4)-1.369306393762915*nuVtSqSum[0]*fSkin[19]*rdvSq4-0.9486832980505137*nuVtSqSum[1]*fSkin[17]*rdvSq4-1.369306393762915*nuVtSqSum[1]*fSkin[16]*rdvSq4-0.5477225575051661*nuVtSqSum[1]*fSkin[11]*rdvSq4-0.9486832980505137*nuVtSqSum[2]*fSkin[10]*rdvSq4-1.060660171779821*nuVtSqSum[0]*fSkin[10]*rdvSq4-1.060660171779821*nuVtSqSum[1]*fSkin[6]*rdvSq4-0.5477225575051661*nuVtSqSum[2]*fSkin[4]*rdvSq4-0.6123724356957944*nuVtSqSum[0]*fSkin[4]*rdvSq4+1.224744871391589*Gdiff2[3]*rdvSq4-0.6123724356957944*nuVtSqSum[1]*fSkin[2]*rdvSq4+1.224744871391589*Ghat[3]*rdv2; 
  out[11] += -0.7071067811865475*Ghat[6]*rdv2; 
  out[12] += -0.7071067811865475*Ghat[7]*rdv2; 
  out[13] += (-1.224744871391589*nuVtSqSum[1]*fSkin[15]*rdvSq4)-0.6776309271789384*nuVtSqSum[2]*fSkin[13]*rdvSq4-1.060660171779821*nuVtSqSum[0]*fSkin[13]*rdvSq4-1.369306393762915*nuVtSqSum[2]*fSkin[9]*rdvSq4-0.3912303982179757*nuVtSqSum[2]*fSkin[7]*rdvSq4-0.6123724356957944*nuVtSqSum[0]*fSkin[7]*rdvSq4-0.9486832980505138*nuVtSqSum[1]*fSkin[5]*rdvSq4+1.224744871391589*Gdiff2[4]*rdvSq4-1.060660171779821*nuVtSqSum[2]*fSkin[3]*rdvSq4-0.6123724356957944*fSkin[0]*nuVtSqSum[2]*rdvSq4-0.5477225575051661*fSkin[1]*nuVtSqSum[1]*rdvSq4+1.224744871391589*Ghat[4]*rdv2; 
  out[14] += (-1.060660171779821*nuVtSqSum[1]*fSkin[18]*rdvSq4)-1.060660171779821*nuVtSqSum[0]*fSkin[14]*rdvSq4-0.6123724356957944*nuVtSqSum[1]*fSkin[12]*rdvSq4-0.6123724356957944*nuVtSqSum[0]*fSkin[8]*rdvSq4+1.224744871391589*Gdiff2[5]*rdvSq4+1.224744871391589*Ghat[5]*rdv2; 
  out[15] += (-4.743416490252569*nuVtSqSum[2]*fSkin[15]*rdvSq4)-5.303300858899105*nuVtSqSum[0]*fSkin[15]*rdvSq4-3.674234614174766*nuVtSqSum[1]*fSkin[13]*rdvSq4-5.303300858899106*nuVtSqSum[1]*fSkin[9]*rdvSq4-2.121320343559642*nuVtSqSum[1]*fSkin[7]*rdvSq4-3.674234614174767*nuVtSqSum[2]*fSkin[5]*rdvSq4-4.107919181288746*nuVtSqSum[0]*fSkin[5]*rdvSq4-4.107919181288746*nuVtSqSum[1]*fSkin[3]*rdvSq4-2.121320343559642*fSkin[1]*nuVtSqSum[2]*rdvSq4-2.371708245126284*fSkin[0]*nuVtSqSum[1]*rdvSq4-2.371708245126284*nuVtSqSum[0]*fSkin[1]*rdvSq4-4.743416490252569*Gdiff2[1]*rdvSq4-1.58113883008419*Ghat[1]*rdv2; 
  out[16] += (-5.303300858899106*nuVtSqSum[1]*fSkin[19]*rdvSq4)-4.107919181288745*nuVtSqSum[2]*fSkin[17]*rdvSq4-5.303300858899105*nuVtSqSum[0]*fSkin[16]*rdvSq4-2.371708245126284*nuVtSqSum[2]*fSkin[11]*rdvSq4-4.107919181288745*nuVtSqSum[1]*fSkin[10]*rdvSq4-4.107919181288745*nuVtSqSum[0]*fSkin[6]*rdvSq4-2.371708245126284*nuVtSqSum[1]*fSkin[4]*rdvSq4-2.371708245126284*nuVtSqSum[0]*fSkin[2]*rdvSq4-4.743416490252569*Gdiff2[2]*rdvSq4-1.58113883008419*Ghat[2]*rdv2; 
  out[17] += (-1.224744871391589*nuVtSqSum[1]*fSkin[19]*rdvSq4)-0.6776309271789384*nuVtSqSum[2]*fSkin[17]*rdvSq4-1.060660171779821*nuVtSqSum[0]*fSkin[17]*rdvSq4-1.369306393762915*nuVtSqSum[2]*fSkin[16]*rdvSq4-0.3912303982179757*nuVtSqSum[2]*fSkin[11]*rdvSq4-0.6123724356957944*nuVtSqSum[0]*fSkin[11]*rdvSq4-0.9486832980505137*nuVtSqSum[1]*fSkin[10]*rdvSq4-1.060660171779821*nuVtSqSum[2]*fSkin[6]*rdvSq4+1.224744871391589*Gdiff2[6]*rdvSq4-0.5477225575051661*nuVtSqSum[1]*fSkin[4]*rdvSq4-0.6123724356957944*fSkin[2]*nuVtSqSum[2]*rdvSq4+1.224744871391589*Ghat[6]*rdv2; 
  out[18] += (-0.9486832980505137*nuVtSqSum[2]*fSkin[18]*rdvSq4)-1.060660171779821*nuVtSqSum[0]*fSkin[18]*rdvSq4-1.060660171779821*nuVtSqSum[1]*fSkin[14]*rdvSq4-0.5477225575051661*nuVtSqSum[2]*fSkin[12]*rdvSq4-0.6123724356957944*nuVtSqSum[0]*fSkin[12]*rdvSq4-0.6123724356957944*nuVtSqSum[1]*fSkin[8]*rdvSq4+1.224744871391589*Gdiff2[7]*rdvSq4+1.224744871391589*Ghat[7]*rdv2; 
  out[19] += (-4.743416490252569*nuVtSqSum[2]*fSkin[19]*rdvSq4)-5.303300858899105*nuVtSqSum[0]*fSkin[19]*rdvSq4-3.674234614174766*nuVtSqSum[1]*fSkin[17]*rdvSq4-5.303300858899106*nuVtSqSum[1]*fSkin[16]*rdvSq4-2.121320343559642*nuVtSqSum[1]*fSkin[11]*rdvSq4-3.674234614174766*nuVtSqSum[2]*fSkin[10]*rdvSq4-4.107919181288745*nuVtSqSum[0]*fSkin[10]*rdvSq4-4.107919181288745*nuVtSqSum[1]*fSkin[6]*rdvSq4-2.121320343559642*nuVtSqSum[2]*fSkin[4]*rdvSq4-2.371708245126284*nuVtSqSum[0]*fSkin[4]*rdvSq4-4.743416490252569*Gdiff2[3]*rdvSq4-2.371708245126284*nuVtSqSum[1]*fSkin[2]*rdvSq4-1.58113883008419*Ghat[3]*rdv2; 

  } 
} 
