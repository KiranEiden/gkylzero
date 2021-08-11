#include <gkyl_lbo_kernels.h> 
GKYL_CU_DH void lbo_constNu_surf_1x1v_ser_vx_p2(const double *w, const double *dxv, const double nuSum, const double *nuUSumL, const double *nuUSumR, const double *nuVtSqSumL, const double *nuVtSqSumR, const double *fl, const double *fc, const double *fr, double *out) 
{ 
  // w[2]:          Cell-center coordinates. 
  // dxv[2]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[3]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[1]; 
  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 

  const double *sumNuUxL = &nuUSumL[0]; 
  const double *sumNuUxR = &nuUSumR[0]; 

  double alphaDrSurfL[3]; 
  alphaDrSurfL[0] = 0.7071067811865475*((2.0*w[1]+dxv[1])*nuSum-1.414213562373095*sumNuUxL[0]); 
  alphaDrSurfL[1] = -1.0*sumNuUxL[1]; 
  alphaDrSurfL[2] = -1.0*sumNuUxL[2]; 

  double alphaDrSurfR[3]; 
  alphaDrSurfR[0] = -0.7071067811865475*((2.0*w[1]+dxv[1])*nuSum-1.414213562373095*sumNuUxR[0]); 
  alphaDrSurfR[1] = sumNuUxR[1]; 
  alphaDrSurfR[2] = sumNuUxR[2]; 

  double fUpwindQuadL[3];
  fUpwindQuadL[0] = 0.5*((-1.5*(fl[7]+fc[7]))+0.7745966692414833*fl[6]-0.7745966692414833*fc[6]+1.118033988749895*(fl[5]+fc[5])+0.4472135954999579*(fl[4]+fc[4])-1.161895003862224*fl[3]+1.161895003862224*fc[3]+0.8660254037844386*fl[2]-0.8660254037844386*fc[2]-0.6708203932499369*(fl[1]+fc[1])+0.5*(fl[0]+fc[0]))-(0.5*(0.6324555320336759*alphaDrSurfL[2]-0.9486832980505137*alphaDrSurfL[1]+0.7071067811865475*alphaDrSurfL[0])*((-1.5*fl[7])+1.5*fc[7]+0.7745966692414833*(fl[6]+fc[6])+1.118033988749895*fl[5]-1.118033988749895*fc[5]+0.4472135954999579*fl[4]-0.4472135954999579*fc[4]-1.161895003862224*(fl[3]+fc[3])+0.8660254037844386*(fl[2]+fc[2])-0.6708203932499369*fl[1]+0.6708203932499369*fc[1]+0.5*fl[0]-0.5*fc[0]))/fabs(0.6324555320336759*alphaDrSurfL[2]-0.9486832980505137*alphaDrSurfL[1]+0.7071067811865475*alphaDrSurfL[0]); 
  fUpwindQuadL[1] = 0.5*((-0.9682458365518543*fl[6])+0.9682458365518543*fc[6]+1.118033988749895*(fl[5]+fc[5])-0.5590169943749475*(fl[4]+fc[4])+0.8660254037844386*fl[2]-0.8660254037844386*fc[2]+0.5*(fl[0]+fc[0]))-(0.5*(0.7071067811865475*alphaDrSurfL[0]-0.7905694150420947*alphaDrSurfL[2])*((-0.9682458365518543*(fl[6]+fc[6]))+1.118033988749895*fl[5]-1.118033988749895*fc[5]-0.5590169943749475*fl[4]+0.5590169943749475*fc[4]+0.8660254037844386*(fl[2]+fc[2])+0.5*fl[0]-0.5*fc[0]))/fabs(0.7071067811865475*alphaDrSurfL[0]-0.7905694150420947*alphaDrSurfL[2]); 
  fUpwindQuadL[2] = 0.5*(1.5*(fl[7]+fc[7])+0.7745966692414833*fl[6]-0.7745966692414833*fc[6]+1.118033988749895*(fl[5]+fc[5])+0.4472135954999579*(fl[4]+fc[4])+1.161895003862224*fl[3]-1.161895003862224*fc[3]+0.8660254037844386*fl[2]-0.8660254037844386*fc[2]+0.6708203932499369*(fl[1]+fc[1])+0.5*(fl[0]+fc[0]))-(0.5*(0.6324555320336759*alphaDrSurfL[2]+0.9486832980505137*alphaDrSurfL[1]+0.7071067811865475*alphaDrSurfL[0])*(1.5*fl[7]-1.5*fc[7]+0.7745966692414833*(fl[6]+fc[6])+1.118033988749895*fl[5]-1.118033988749895*fc[5]+0.4472135954999579*fl[4]-0.4472135954999579*fc[4]+1.161895003862224*(fl[3]+fc[3])+0.8660254037844386*(fl[2]+fc[2])+0.6708203932499369*fl[1]-0.6708203932499369*fc[1]+0.5*fl[0]-0.5*fc[0]))/fabs(0.6324555320336759*alphaDrSurfL[2]+0.9486832980505137*alphaDrSurfL[1]+0.7071067811865475*alphaDrSurfL[0]); 

  double fUpwindL[3];
  fUpwindL[0] = 0.07856742013183861*(5.0*fUpwindQuadL[2]+8.0*fUpwindQuadL[1]+5.0*fUpwindQuadL[0]); 
  fUpwindL[1] = 0.5270462766947298*(fUpwindQuadL[2]-1.0*fUpwindQuadL[0]); 
  fUpwindL[2] = 0.3513641844631532*(fUpwindQuadL[2]-2.0*fUpwindQuadL[1]+fUpwindQuadL[0]); 

  double fUpwindQuadR[3];
  fUpwindQuadR[0] = 0.5*((-1.5*(fr[7]+fc[7]))-0.7745966692414833*fr[6]+0.7745966692414833*fc[6]+1.118033988749895*(fr[5]+fc[5])+0.4472135954999579*(fr[4]+fc[4])+1.161895003862224*fr[3]-1.161895003862224*fc[3]-0.8660254037844386*fr[2]+0.8660254037844386*fc[2]-0.6708203932499369*(fr[1]+fc[1])+0.5*(fr[0]+fc[0]))-(0.5*(0.6324555320336759*alphaDrSurfR[2]-0.9486832980505137*alphaDrSurfR[1]+0.7071067811865475*alphaDrSurfR[0])*(1.5*fr[7]-1.5*fc[7]+0.7745966692414833*(fr[6]+fc[6])-1.118033988749895*fr[5]+1.118033988749895*fc[5]-0.4472135954999579*fr[4]+0.4472135954999579*fc[4]-1.161895003862224*(fr[3]+fc[3])+0.8660254037844386*(fr[2]+fc[2])+0.6708203932499369*fr[1]-0.6708203932499369*fc[1]-0.5*fr[0]+0.5*fc[0]))/fabs(0.6324555320336759*alphaDrSurfR[2]-0.9486832980505137*alphaDrSurfR[1]+0.7071067811865475*alphaDrSurfR[0]); 
  fUpwindQuadR[1] = 0.5*(0.9682458365518543*fr[6]-0.9682458365518543*fc[6]+1.118033988749895*(fr[5]+fc[5])-0.5590169943749475*(fr[4]+fc[4])-0.8660254037844386*fr[2]+0.8660254037844386*fc[2]+0.5*(fr[0]+fc[0]))-(0.5*(0.7071067811865475*alphaDrSurfR[0]-0.7905694150420947*alphaDrSurfR[2])*((-0.9682458365518543*(fr[6]+fc[6]))-1.118033988749895*fr[5]+1.118033988749895*fc[5]+0.5590169943749475*fr[4]-0.5590169943749475*fc[4]+0.8660254037844386*(fr[2]+fc[2])-0.5*fr[0]+0.5*fc[0]))/fabs(0.7071067811865475*alphaDrSurfR[0]-0.7905694150420947*alphaDrSurfR[2]); 
  fUpwindQuadR[2] = 0.5*(1.5*(fr[7]+fc[7])-0.7745966692414833*fr[6]+0.7745966692414833*fc[6]+1.118033988749895*(fr[5]+fc[5])+0.4472135954999579*(fr[4]+fc[4])-1.161895003862224*fr[3]+1.161895003862224*fc[3]-0.8660254037844386*fr[2]+0.8660254037844386*fc[2]+0.6708203932499369*(fr[1]+fc[1])+0.5*(fr[0]+fc[0]))-(0.5*(0.6324555320336759*alphaDrSurfR[2]+0.9486832980505137*alphaDrSurfR[1]+0.7071067811865475*alphaDrSurfR[0])*((-1.5*fr[7])+1.5*fc[7]+0.7745966692414833*(fr[6]+fc[6])-1.118033988749895*fr[5]+1.118033988749895*fc[5]-0.4472135954999579*fr[4]+0.4472135954999579*fc[4]+1.161895003862224*(fr[3]+fc[3])+0.8660254037844386*(fr[2]+fc[2])-0.6708203932499369*fr[1]+0.6708203932499369*fc[1]-0.5*fr[0]+0.5*fc[0]))/fabs(0.6324555320336759*alphaDrSurfR[2]+0.9486832980505137*alphaDrSurfR[1]+0.7071067811865475*alphaDrSurfR[0]); 

  double fUpwindR[3];
  fUpwindR[0] = 0.07856742013183861*(5.0*fUpwindQuadR[2]+8.0*fUpwindQuadR[1]+5.0*fUpwindQuadR[0]); 
  fUpwindR[1] = 0.5270462766947298*(fUpwindQuadR[2]-1.0*fUpwindQuadR[0]); 
  fUpwindR[2] = 0.3513641844631532*(fUpwindQuadR[2]-2.0*fUpwindQuadR[1]+fUpwindQuadR[0]); 

  double GdiffL[8]; 
  double GdiffR[8]; 
  double GhatL[8]; 
  double GhatR[8]; 
  double incr2_l[8]; 
  double incr2_r[8]; 


  incr2_l[2] = -0.0078125*(38.34057902536163*nuVtSqSumL[1]*fl[7]+38.34057902536163*nuVtSqSumL[1]*fc[7]+55.15432893255071*nuVtSqSumL[2]*fl[6]-55.15432893255071*nuVtSqSumL[2]*fc[6]+38.34057902536163*nuVtSqSumL[0]*fl[5]+38.34057902536163*nuVtSqSumL[0]*fc[5]+39.19183588453089*nuVtSqSumL[2]*fl[4]+39.19183588453089*nuVtSqSumL[2]*fc[4]+55.15432893255071*nuVtSqSumL[1]*fl[3]-55.15432893255071*nuVtSqSumL[1]*fc[3]+55.15432893255071*nuVtSqSumL[0]*fl[2]-55.15432893255071*nuVtSqSumL[0]*fc[2]+(39.19183588453089*fl[1]+39.19183588453089*fc[1])*nuVtSqSumL[1]+(39.19183588453089*fl[0]+39.19183588453089*fc[0])*nuVtSqSumL[0]); 
  incr2_l[3] = -0.0015625*((171.4642819948225*nuVtSqSumL[2]+191.7028951268081*nuVtSqSumL[0])*fl[7]+(171.4642819948225*nuVtSqSumL[2]+191.7028951268081*nuVtSqSumL[0])*fc[7]+246.6576574931336*nuVtSqSumL[1]*fl[6]-246.6576574931336*nuVtSqSumL[1]*fc[6]+191.7028951268081*nuVtSqSumL[1]*fl[5]+191.7028951268081*nuVtSqSumL[1]*fc[5]+175.2712184016533*nuVtSqSumL[1]*fl[4]+175.2712184016533*nuVtSqSumL[1]*fc[4]+(246.6576574931336*nuVtSqSumL[2]+275.7716446627535*nuVtSqSumL[0])*fl[3]+((-246.6576574931336*nuVtSqSumL[2])-275.7716446627535*nuVtSqSumL[0])*fc[3]+(175.2712184016533*fl[1]+175.2712184016533*fc[1])*nuVtSqSumL[2]+275.7716446627535*nuVtSqSumL[1]*fl[2]-275.7716446627535*nuVtSqSumL[1]*fc[2]+(195.9591794226544*fl[0]+195.9591794226544*fc[0])*nuVtSqSumL[1]+195.9591794226544*nuVtSqSumL[0]*fl[1]+195.9591794226544*nuVtSqSumL[0]*fc[1]); 
  incr2_l[5] = 0.0078125*(148.4924240491749*nuVtSqSumL[1]*fl[7]+148.4924240491749*nuVtSqSumL[1]*fc[7]+213.6117974270148*nuVtSqSumL[2]*fl[6]-213.6117974270148*nuVtSqSumL[2]*fc[6]+148.4924240491749*nuVtSqSumL[0]*fl[5]+148.4924240491749*nuVtSqSumL[0]*fc[5]+151.7893276880823*nuVtSqSumL[2]*fl[4]+151.7893276880823*nuVtSqSumL[2]*fc[4]+213.6117974270148*nuVtSqSumL[1]*fl[3]-213.6117974270148*nuVtSqSumL[1]*fc[3]+213.6117974270148*nuVtSqSumL[0]*fl[2]-213.6117974270148*nuVtSqSumL[0]*fc[2]+(151.7893276880823*fl[1]+151.7893276880823*fc[1])*nuVtSqSumL[1]+(151.7893276880823*fl[0]+151.7893276880823*fc[0])*nuVtSqSumL[0]); 
  incr2_l[6] = -2.232142857142857E-4*(1200.249973963757*nuVtSqSumL[1]*fl[7]+1200.249973963757*nuVtSqSumL[1]*fc[7]+(1233.288287465668*nuVtSqSumL[2]+1930.401512639274*nuVtSqSumL[0])*fl[6]+((-1233.288287465668*nuVtSqSumL[2])-1930.401512639274*nuVtSqSumL[0])*fc[6]+1341.920265887657*nuVtSqSumL[2]*fl[5]+1341.920265887657*nuVtSqSumL[2]*fc[5]+(876.356092008267*nuVtSqSumL[2]+1371.714255958581*nuVtSqSumL[0])*fl[4]+(876.356092008267*nuVtSqSumL[2]+1371.714255958581*nuVtSqSumL[0])*fc[4]+1726.603602451935*nuVtSqSumL[1]*fl[3]-1726.603602451935*nuVtSqSumL[1]*fc[3]+(1930.401512639274*fl[2]-1930.401512639274*fc[2]+1371.714255958581*fl[0]+1371.714255958581*fc[0])*nuVtSqSumL[2]+(1226.898528811573*fl[1]+1226.898528811573*fc[1])*nuVtSqSumL[1]); 
  incr2_l[7] = 0.0078125*((132.8156617270719*nuVtSqSumL[2]+148.4924240491749*nuVtSqSumL[0])*fl[7]+(132.8156617270719*nuVtSqSumL[2]+148.4924240491749*nuVtSqSumL[0])*fc[7]+191.0601999370879*nuVtSqSumL[1]*fl[6]-191.0601999370879*nuVtSqSumL[1]*fc[6]+148.4924240491749*nuVtSqSumL[1]*fl[5]+148.4924240491749*nuVtSqSumL[1]*fc[5]+135.7645019878172*nuVtSqSumL[1]*fl[4]+135.7645019878172*nuVtSqSumL[1]*fc[4]+(191.0601999370879*nuVtSqSumL[2]+213.6117974270148*nuVtSqSumL[0])*fl[3]+((-191.0601999370879*nuVtSqSumL[2])-213.6117974270148*nuVtSqSumL[0])*fc[3]+(135.7645019878172*fl[1]+135.7645019878172*fc[1])*nuVtSqSumL[2]+213.6117974270148*nuVtSqSumL[1]*fl[2]-213.6117974270148*nuVtSqSumL[1]*fc[2]+(151.7893276880823*fl[0]+151.7893276880823*fc[0])*nuVtSqSumL[1]+151.7893276880823*nuVtSqSumL[0]*fl[1]+151.7893276880823*nuVtSqSumL[0]*fc[1]); 


  incr2_r[2] = 0.0078125*(38.34057902536163*nuVtSqSumR[1]*fr[7]+38.34057902536163*nuVtSqSumR[1]*fc[7]-55.15432893255071*nuVtSqSumR[2]*fr[6]+55.15432893255071*nuVtSqSumR[2]*fc[6]+38.34057902536163*nuVtSqSumR[0]*fr[5]+38.34057902536163*nuVtSqSumR[0]*fc[5]+39.19183588453089*nuVtSqSumR[2]*fr[4]+39.19183588453089*nuVtSqSumR[2]*fc[4]-55.15432893255071*nuVtSqSumR[1]*fr[3]+55.15432893255071*nuVtSqSumR[1]*fc[3]-55.15432893255071*nuVtSqSumR[0]*fr[2]+55.15432893255071*nuVtSqSumR[0]*fc[2]+(39.19183588453089*fr[1]+39.19183588453089*fc[1])*nuVtSqSumR[1]+(39.19183588453089*fr[0]+39.19183588453089*fc[0])*nuVtSqSumR[0]); 
  incr2_r[3] = 0.0015625*((171.4642819948225*nuVtSqSumR[2]+191.7028951268081*nuVtSqSumR[0])*fr[7]+(171.4642819948225*nuVtSqSumR[2]+191.7028951268081*nuVtSqSumR[0])*fc[7]-246.6576574931336*nuVtSqSumR[1]*fr[6]+246.6576574931336*nuVtSqSumR[1]*fc[6]+191.7028951268081*nuVtSqSumR[1]*fr[5]+191.7028951268081*nuVtSqSumR[1]*fc[5]+175.2712184016533*nuVtSqSumR[1]*fr[4]+175.2712184016533*nuVtSqSumR[1]*fc[4]+((-246.6576574931336*nuVtSqSumR[2])-275.7716446627535*nuVtSqSumR[0])*fr[3]+(246.6576574931336*nuVtSqSumR[2]+275.7716446627535*nuVtSqSumR[0])*fc[3]+(175.2712184016533*fr[1]+175.2712184016533*fc[1])*nuVtSqSumR[2]-275.7716446627535*nuVtSqSumR[1]*fr[2]+275.7716446627535*nuVtSqSumR[1]*fc[2]+(195.9591794226544*fr[0]+195.9591794226544*fc[0])*nuVtSqSumR[1]+195.9591794226544*nuVtSqSumR[0]*fr[1]+195.9591794226544*nuVtSqSumR[0]*fc[1]); 
  incr2_r[5] = 0.0078125*(148.4924240491749*nuVtSqSumR[1]*fr[7]+148.4924240491749*nuVtSqSumR[1]*fc[7]-213.6117974270148*nuVtSqSumR[2]*fr[6]+213.6117974270148*nuVtSqSumR[2]*fc[6]+148.4924240491749*nuVtSqSumR[0]*fr[5]+148.4924240491749*nuVtSqSumR[0]*fc[5]+151.7893276880823*nuVtSqSumR[2]*fr[4]+151.7893276880823*nuVtSqSumR[2]*fc[4]-213.6117974270148*nuVtSqSumR[1]*fr[3]+213.6117974270148*nuVtSqSumR[1]*fc[3]-213.6117974270148*nuVtSqSumR[0]*fr[2]+213.6117974270148*nuVtSqSumR[0]*fc[2]+(151.7893276880823*fr[1]+151.7893276880823*fc[1])*nuVtSqSumR[1]+(151.7893276880823*fr[0]+151.7893276880823*fc[0])*nuVtSqSumR[0]); 
  incr2_r[6] = 2.232142857142857E-4*(1200.249973963757*nuVtSqSumR[1]*fr[7]+1200.249973963757*nuVtSqSumR[1]*fc[7]+((-1233.288287465668*nuVtSqSumR[2])-1930.401512639274*nuVtSqSumR[0])*fr[6]+(1233.288287465668*nuVtSqSumR[2]+1930.401512639274*nuVtSqSumR[0])*fc[6]+1341.920265887657*nuVtSqSumR[2]*fr[5]+1341.920265887657*nuVtSqSumR[2]*fc[5]+(876.356092008267*nuVtSqSumR[2]+1371.714255958581*nuVtSqSumR[0])*fr[4]+(876.356092008267*nuVtSqSumR[2]+1371.714255958581*nuVtSqSumR[0])*fc[4]-1726.603602451935*nuVtSqSumR[1]*fr[3]+1726.603602451935*nuVtSqSumR[1]*fc[3]+((-1930.401512639274*fr[2])+1930.401512639274*fc[2]+1371.714255958581*fr[0]+1371.714255958581*fc[0])*nuVtSqSumR[2]+(1226.898528811573*fr[1]+1226.898528811573*fc[1])*nuVtSqSumR[1]); 
  incr2_r[7] = 0.0078125*((132.8156617270719*nuVtSqSumR[2]+148.4924240491749*nuVtSqSumR[0])*fr[7]+(132.8156617270719*nuVtSqSumR[2]+148.4924240491749*nuVtSqSumR[0])*fc[7]-191.0601999370879*nuVtSqSumR[1]*fr[6]+191.0601999370879*nuVtSqSumR[1]*fc[6]+148.4924240491749*nuVtSqSumR[1]*fr[5]+148.4924240491749*nuVtSqSumR[1]*fc[5]+135.7645019878172*nuVtSqSumR[1]*fr[4]+135.7645019878172*nuVtSqSumR[1]*fc[4]+((-191.0601999370879*nuVtSqSumR[2])-213.6117974270148*nuVtSqSumR[0])*fr[3]+(191.0601999370879*nuVtSqSumR[2]+213.6117974270148*nuVtSqSumR[0])*fc[3]+(135.7645019878172*fr[1]+135.7645019878172*fc[1])*nuVtSqSumR[2]-213.6117974270148*nuVtSqSumR[1]*fr[2]+213.6117974270148*nuVtSqSumR[1]*fc[2]+(151.7893276880823*fr[0]+151.7893276880823*fc[0])*nuVtSqSumR[1]+151.7893276880823*nuVtSqSumR[0]*fr[1]+151.7893276880823*nuVtSqSumR[0]*fc[1]); 


  GdiffL[0] = -0.0125*(75.89466384404116*nuVtSqSumL[1]*fl[7]-75.89466384404116*nuVtSqSumL[1]*fc[7]+134.7219358530748*nuVtSqSumL[2]*fl[6]+134.7219358530748*nuVtSqSumL[2]*fc[6]+75.89466384404116*nuVtSqSumL[0]*fl[5]-75.89466384404116*nuVtSqSumL[0]*fc[5]+106.0660171779821*nuVtSqSumL[2]*fl[4]-106.0660171779821*nuVtSqSumL[2]*fc[4]+134.7219358530748*nuVtSqSumL[1]*fl[3]+134.7219358530748*nuVtSqSumL[1]*fc[3]+134.7219358530748*nuVtSqSumL[0]*fl[2]+134.7219358530748*nuVtSqSumL[0]*fc[2]+(106.0660171779821*fl[1]-106.0660171779821*fc[1])*nuVtSqSumL[1]+(106.0660171779821*fl[0]-106.0660171779821*fc[0])*nuVtSqSumL[0]); 
  GdiffL[1] = -0.0025*((339.4112549695431*nuVtSqSumL[2]+379.4733192202057*nuVtSqSumL[0])*fl[7]+((-339.4112549695431*nuVtSqSumL[2])-379.4733192202057*nuVtSqSumL[0])*fc[7]+602.4948132556829*nuVtSqSumL[1]*fl[6]+602.4948132556829*nuVtSqSumL[1]*fc[6]+379.4733192202059*nuVtSqSumL[1]*fl[5]-379.4733192202059*nuVtSqSumL[1]*fc[5]+474.3416490252572*nuVtSqSumL[1]*fl[4]-474.3416490252572*nuVtSqSumL[1]*fc[4]+(602.494813255683*nuVtSqSumL[2]+673.609679265374*nuVtSqSumL[0])*fl[3]+(602.494813255683*nuVtSqSumL[2]+673.609679265374*nuVtSqSumL[0])*fc[3]+(474.3416490252572*fl[1]-474.3416490252572*fc[1])*nuVtSqSumL[2]+673.609679265374*nuVtSqSumL[1]*fl[2]+673.609679265374*nuVtSqSumL[1]*fc[2]+(530.3300858899107*fl[0]-530.3300858899107*fc[0])*nuVtSqSumL[1]+530.3300858899107*nuVtSqSumL[0]*fl[1]-530.3300858899107*nuVtSqSumL[0]*fc[1]); 
  GdiffL[4] = -3.571428571428571E-4*(2375.878784786802*nuVtSqSumL[1]*fl[7]-2375.878784786802*nuVtSqSumL[1]*fc[7]+(3012.474066278414*nuVtSqSumL[2]+4715.267754857619*nuVtSqSumL[0])*fl[6]+(3012.474066278414*nuVtSqSumL[2]+4715.267754857619*nuVtSqSumL[0])*fc[6]+2656.313234541441*nuVtSqSumL[2]*fl[5]-2656.313234541441*nuVtSqSumL[2]*fc[5]+(2371.708245126286*nuVtSqSumL[2]+3712.310601229374*nuVtSqSumL[0])*fl[4]+((-2371.708245126286*nuVtSqSumL[2])-3712.310601229374*nuVtSqSumL[0])*fc[4]+4217.463692789781*nuVtSqSumL[1]*fl[3]+4217.463692789781*nuVtSqSumL[1]*fc[3]+(4715.267754857618*fl[2]+4715.267754857618*fc[2]+3712.310601229374*fl[0]-3712.310601229374*fc[0])*nuVtSqSumL[2]+(3320.3915431768*fl[1]-3320.3915431768*fc[1])*nuVtSqSumL[1]); 


  GdiffR[0] = 0.0125*(75.89466384404116*nuVtSqSumR[1]*fr[7]-75.89466384404116*nuVtSqSumR[1]*fc[7]-134.7219358530748*nuVtSqSumR[2]*fr[6]-134.7219358530748*nuVtSqSumR[2]*fc[6]+75.89466384404116*nuVtSqSumR[0]*fr[5]-75.89466384404116*nuVtSqSumR[0]*fc[5]+106.0660171779821*nuVtSqSumR[2]*fr[4]-106.0660171779821*nuVtSqSumR[2]*fc[4]-134.7219358530748*nuVtSqSumR[1]*fr[3]-134.7219358530748*nuVtSqSumR[1]*fc[3]-134.7219358530748*nuVtSqSumR[0]*fr[2]-134.7219358530748*nuVtSqSumR[0]*fc[2]+(106.0660171779821*fr[1]-106.0660171779821*fc[1])*nuVtSqSumR[1]+(106.0660171779821*fr[0]-106.0660171779821*fc[0])*nuVtSqSumR[0]); 
  GdiffR[1] = 0.0025*((339.4112549695431*nuVtSqSumR[2]+379.4733192202057*nuVtSqSumR[0])*fr[7]+((-339.4112549695431*nuVtSqSumR[2])-379.4733192202057*nuVtSqSumR[0])*fc[7]-602.4948132556829*nuVtSqSumR[1]*fr[6]-602.4948132556829*nuVtSqSumR[1]*fc[6]+379.4733192202059*nuVtSqSumR[1]*fr[5]-379.4733192202059*nuVtSqSumR[1]*fc[5]+474.3416490252572*nuVtSqSumR[1]*fr[4]-474.3416490252572*nuVtSqSumR[1]*fc[4]+((-602.494813255683*nuVtSqSumR[2])-673.609679265374*nuVtSqSumR[0])*fr[3]+((-602.494813255683*nuVtSqSumR[2])-673.609679265374*nuVtSqSumR[0])*fc[3]+(474.3416490252572*fr[1]-474.3416490252572*fc[1])*nuVtSqSumR[2]-673.609679265374*nuVtSqSumR[1]*fr[2]-673.609679265374*nuVtSqSumR[1]*fc[2]+(530.3300858899107*fr[0]-530.3300858899107*fc[0])*nuVtSqSumR[1]+530.3300858899107*nuVtSqSumR[0]*fr[1]-530.3300858899107*nuVtSqSumR[0]*fc[1]); 
  GdiffR[4] = 3.571428571428571E-4*(2375.878784786802*nuVtSqSumR[1]*fr[7]-2375.878784786802*nuVtSqSumR[1]*fc[7]+((-3012.474066278414*nuVtSqSumR[2])-4715.267754857619*nuVtSqSumR[0])*fr[6]+((-3012.474066278414*nuVtSqSumR[2])-4715.267754857619*nuVtSqSumR[0])*fc[6]+2656.313234541441*nuVtSqSumR[2]*fr[5]-2656.313234541441*nuVtSqSumR[2]*fc[5]+(2371.708245126286*nuVtSqSumR[2]+3712.310601229374*nuVtSqSumR[0])*fr[4]+((-2371.708245126286*nuVtSqSumR[2])-3712.310601229374*nuVtSqSumR[0])*fc[4]-4217.463692789781*nuVtSqSumR[1]*fr[3]-4217.463692789781*nuVtSqSumR[1]*fc[3]+((-4715.267754857618*fr[2])-4715.267754857618*fc[2]+3712.310601229374*fr[0]-3712.310601229374*fc[0])*nuVtSqSumR[2]+(3320.3915431768*fr[1]-3320.3915431768*fc[1])*nuVtSqSumR[1]); 

  GhatL[0] = GdiffL[0]*rdv2+alphaDrSurfL[2]*fUpwindL[2]+alphaDrSurfL[1]*fUpwindL[1]+alphaDrSurfL[0]*fUpwindL[0]; 
  GhatL[1] = GdiffL[1]*rdv2+0.8944271909999159*alphaDrSurfL[1]*fUpwindL[2]+0.8944271909999159*fUpwindL[1]*alphaDrSurfL[2]+alphaDrSurfL[0]*fUpwindL[1]+fUpwindL[0]*alphaDrSurfL[1]; 
  GhatL[4] = GdiffL[4]*rdv2+0.6388765649999399*alphaDrSurfL[2]*fUpwindL[2]+alphaDrSurfL[0]*fUpwindL[2]+fUpwindL[0]*alphaDrSurfL[2]+0.8944271909999159*alphaDrSurfL[1]*fUpwindL[1]; 

  GhatR[0] = (-1.0*GdiffR[0]*rdv2)+alphaDrSurfR[2]*fUpwindR[2]+alphaDrSurfR[1]*fUpwindR[1]+alphaDrSurfR[0]*fUpwindR[0]; 
  GhatR[1] = (-1.0*GdiffR[1]*rdv2)+0.8944271909999159*alphaDrSurfR[1]*fUpwindR[2]+0.8944271909999159*fUpwindR[1]*alphaDrSurfR[2]+alphaDrSurfR[0]*fUpwindR[1]+fUpwindR[0]*alphaDrSurfR[1]; 
  GhatR[4] = (-1.0*GdiffR[4]*rdv2)+0.6388765649999399*alphaDrSurfR[2]*fUpwindR[2]+alphaDrSurfR[0]*fUpwindR[2]+fUpwindR[0]*alphaDrSurfR[2]+0.8944271909999159*alphaDrSurfR[1]*fUpwindR[1]; 

  out[0] += 0.5*GhatL[0]*rdv2-0.5*GhatR[0]*rdv2; 
  out[1] += 0.5*GhatL[1]*rdv2-0.5*GhatR[1]*rdv2; 
  out[2] += incr2_r[2]*rdvSq4+incr2_l[2]*rdvSq4-0.8660254037844386*GhatR[0]*rdv2-0.8660254037844386*GhatL[0]*rdv2; 
  out[3] += incr2_r[3]*rdvSq4+incr2_l[3]*rdvSq4-0.8660254037844386*GhatR[1]*rdv2-0.8660254037844386*GhatL[1]*rdv2; 
  out[4] += 0.5*GhatL[4]*rdv2-0.5*GhatR[4]*rdv2; 
  out[5] += incr2_r[5]*rdvSq4+incr2_l[5]*rdvSq4-1.118033988749895*GhatR[0]*rdv2+1.118033988749895*GhatL[0]*rdv2; 
  out[6] += incr2_r[6]*rdvSq4+incr2_l[6]*rdvSq4-0.8660254037844387*GhatR[4]*rdv2-0.8660254037844387*GhatL[4]*rdv2; 
  out[7] += incr2_r[7]*rdvSq4+incr2_l[7]*rdvSq4-1.118033988749895*GhatR[1]*rdv2+1.118033988749895*GhatL[1]*rdv2; 
} 
