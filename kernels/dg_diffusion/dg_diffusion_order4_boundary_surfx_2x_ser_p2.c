#include <gkyl_dg_diffusion_kernels.h>

GKYL_CU_DH double dg_diffusion_order4_boundary_surfx_2x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],4.);

  double vol_incr[8] = {0.0}; 

  double edgeSurf_incr[8] = {0.0}; 
  double boundSurf_incr[8] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = 6.708203932499369*coeff[0]*fSkin[4]-6.708203932499369*coeff[0]*fEdge[4]+8.11898816047911*coeff[0]*fSkin[1]+8.11898816047911*coeff[0]*fEdge[1]+4.6875*coeff[0]*fSkin[0]-4.6875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 14.16059535957086*coeff[0]*fSkin[4]-9.077304717673634*coeff[0]*fEdge[4]+15.46875*coeff[0]*fSkin[1]+12.65625*coeff[0]*fEdge[1]+8.11898816047911*coeff[0]*fSkin[0]-8.11898816047911*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 6.708203932499369*coeff[0]*fSkin[6]-6.708203932499369*coeff[0]*fEdge[6]+8.11898816047911*coeff[0]*fSkin[3]+8.11898816047911*coeff[0]*fEdge[3]+4.6875*coeff[0]*fSkin[2]-4.6875*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 14.16059535957087*coeff[0]*fSkin[6]-9.077304717673634*coeff[0]*fEdge[6]+15.46875*coeff[0]*fSkin[3]+12.65625*coeff[0]*fEdge[3]+8.11898816047911*coeff[0]*fSkin[2]-8.11898816047911*coeff[0]*fEdge[2]; 
  edgeSurf_incr[4] = 20.34375*coeff[0]*fSkin[4]-0.65625*coeff[0]*fEdge[4]+15.61296411439865*coeff[0]*fSkin[1]+4.720198453190292*coeff[0]*fEdge[1]+4.192627457812107*coeff[0]*fSkin[0]-4.192627457812107*coeff[0]*fEdge[0]; 
  edgeSurf_incr[5] = 8.118988160479114*coeff[0]*fSkin[7]+8.118988160479114*coeff[0]*fEdge[7]+4.6875*coeff[0]*fSkin[5]-4.6875*coeff[0]*fEdge[5]; 
  edgeSurf_incr[6] = 20.34375*coeff[0]*fSkin[6]-0.65625*coeff[0]*fEdge[6]+15.61296411439865*coeff[0]*fSkin[3]+4.72019845319029*coeff[0]*fEdge[3]+4.192627457812107*coeff[0]*fSkin[2]-4.192627457812107*coeff[0]*fEdge[2]; 
  edgeSurf_incr[7] = 15.46875*coeff[0]*fSkin[7]+12.65625*coeff[0]*fEdge[7]+8.118988160479114*coeff[0]*fSkin[5]-8.118988160479114*coeff[0]*fEdge[5]; 

  boundSurf_incr[1] = 2.8125*coeff[0]*fSkin[1]-5.083290641897234*coeff[0]*fSkin[4]; 
  boundSurf_incr[3] = 2.8125*coeff[0]*fSkin[3]-5.083290641897235*coeff[0]*fSkin[6]; 
  boundSurf_incr[4] = 19.6875*coeff[0]*fSkin[4]-10.89276566120836*coeff[0]*fSkin[1]; 
  boundSurf_incr[6] = 19.6875*coeff[0]*fSkin[6]-10.89276566120836*coeff[0]*fSkin[3]; 
  boundSurf_incr[7] = 2.8125*coeff[0]*fSkin[7]; 

  } else { 

  edgeSurf_incr[0] = 6.708203932499369*coeff[0]*fSkin[4]-6.708203932499369*coeff[0]*fEdge[4]-8.11898816047911*coeff[0]*fSkin[1]-8.11898816047911*coeff[0]*fEdge[1]+4.6875*coeff[0]*fSkin[0]-4.6875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = (-14.16059535957086*coeff[0]*fSkin[4])+9.077304717673634*coeff[0]*fEdge[4]+15.46875*coeff[0]*fSkin[1]+12.65625*coeff[0]*fEdge[1]-8.11898816047911*coeff[0]*fSkin[0]+8.11898816047911*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 6.708203932499369*coeff[0]*fSkin[6]-6.708203932499369*coeff[0]*fEdge[6]-8.11898816047911*coeff[0]*fSkin[3]-8.11898816047911*coeff[0]*fEdge[3]+4.6875*coeff[0]*fSkin[2]-4.6875*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = (-14.16059535957087*coeff[0]*fSkin[6])+9.077304717673634*coeff[0]*fEdge[6]+15.46875*coeff[0]*fSkin[3]+12.65625*coeff[0]*fEdge[3]-8.11898816047911*coeff[0]*fSkin[2]+8.11898816047911*coeff[0]*fEdge[2]; 
  edgeSurf_incr[4] = 20.34375*coeff[0]*fSkin[4]-0.65625*coeff[0]*fEdge[4]-15.61296411439865*coeff[0]*fSkin[1]-4.720198453190292*coeff[0]*fEdge[1]+4.192627457812107*coeff[0]*fSkin[0]-4.192627457812107*coeff[0]*fEdge[0]; 
  edgeSurf_incr[5] = (-8.118988160479114*coeff[0]*fSkin[7])-8.118988160479114*coeff[0]*fEdge[7]+4.6875*coeff[0]*fSkin[5]-4.6875*coeff[0]*fEdge[5]; 
  edgeSurf_incr[6] = 20.34375*coeff[0]*fSkin[6]-0.65625*coeff[0]*fEdge[6]-15.61296411439865*coeff[0]*fSkin[3]-4.72019845319029*coeff[0]*fEdge[3]+4.192627457812107*coeff[0]*fSkin[2]-4.192627457812107*coeff[0]*fEdge[2]; 
  edgeSurf_incr[7] = 15.46875*coeff[0]*fSkin[7]+12.65625*coeff[0]*fEdge[7]-8.118988160479114*coeff[0]*fSkin[5]+8.118988160479114*coeff[0]*fEdge[5]; 

  boundSurf_incr[1] = 5.083290641897234*coeff[0]*fSkin[4]+2.8125*coeff[0]*fSkin[1]; 
  boundSurf_incr[3] = 5.083290641897235*coeff[0]*fSkin[6]+2.8125*coeff[0]*fSkin[3]; 
  boundSurf_incr[4] = 19.6875*coeff[0]*fSkin[4]+10.89276566120836*coeff[0]*fSkin[1]; 
  boundSurf_incr[6] = 19.6875*coeff[0]*fSkin[6]+10.89276566120836*coeff[0]*fSkin[3]; 
  boundSurf_incr[7] = 2.8125*coeff[0]*fSkin[7]; 

  out[0] += -1.0*(vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += -1.0*(vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 
  out[2] += -1.0*(vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac; 
  out[3] += -1.0*(vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*Jfac; 
  out[4] += -1.0*(vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*Jfac; 
  out[5] += -1.0*(vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*Jfac; 
  out[6] += -1.0*(vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*Jfac; 
  out[7] += -1.0*(vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*Jfac; 

  }

  return 0.;
}

GKYL_CU_DH double dg_diffusion_order4_boundary_surfx_2x_ser_p2_varcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],4.);

  double vol_incr[8] = {0.0}; 
  vol_incr[4] = (-22.5*fSkin[2]*coeff[6])-22.5*fSkin[0]*coeff[4]; 
  vol_incr[6] = (-20.12461179749811*fSkin[5]*coeff[6])-22.5*fSkin[0]*coeff[6]-22.5*fSkin[2]*coeff[4]; 

  double edgeSurf_incr[8] = {0.0}; 
  double boundSurf_incr[8] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = (-0.703125*coeff[7]*fSkin[7])+4.059494080239557*coeff[5]*fSkin[7]+0.703125*coeff[7]*fEdge[7]+4.059494080239557*coeff[5]*fEdge[7]-3.75*coeff[6]*fSkin[6]-1.270822660474309*coeff[3]*fSkin[6]+3.354101966249684*coeff[2]*fSkin[6]+3.75*coeff[6]*fEdge[6]-1.270822660474309*coeff[3]*fEdge[6]-3.354101966249684*coeff[2]*fEdge[6]-4.538652358836817*fSkin[3]*coeff[6]-4.538652358836817*fEdge[3]*coeff[6]-2.620392161132566*fSkin[2]*coeff[6]+2.620392161132566*fEdge[2]*coeff[6]+2.34375*coeff[5]*fSkin[5]-2.34375*coeff[5]*fEdge[5]-3.75*coeff[4]*fSkin[4]-1.270822660474308*coeff[1]*fSkin[4]+3.354101966249685*coeff[0]*fSkin[4]+3.75*coeff[4]*fEdge[4]-1.270822660474308*coeff[1]*fEdge[4]-3.354101966249685*coeff[0]*fEdge[4]-4.538652358836817*fSkin[1]*coeff[4]-4.538652358836817*fEdge[1]*coeff[4]-2.620392161132567*fSkin[0]*coeff[4]+2.620392161132567*fEdge[0]*coeff[4]-0.703125*coeff[3]*fSkin[3]+4.059494080239555*coeff[2]*fSkin[3]+0.703125*coeff[3]*fEdge[3]+4.059494080239555*coeff[2]*fEdge[3]+2.34375*coeff[2]*fSkin[2]-2.34375*coeff[2]*fEdge[2]-0.703125*coeff[1]*fSkin[1]+4.059494080239555*coeff[0]*fSkin[1]+0.703125*coeff[1]*fEdge[1]+4.059494080239555*coeff[0]*fEdge[1]+2.34375*coeff[0]*fSkin[0]-2.34375*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = (-1.217848224071866*coeff[7]*fSkin[7])+7.734375000000001*coeff[5]*fSkin[7]+1.217848224071866*coeff[7]*fEdge[7]+6.328125000000001*coeff[5]*fEdge[7]-7.916013456467132*coeff[6]*fSkin[6]-2.201129415351355*coeff[3]*fSkin[6]+7.080297679785434*coeff[2]*fSkin[6]+5.074367600299444*coeff[6]*fEdge[6]-2.201129415351355*coeff[3]*fEdge[6]-4.538652358836817*coeff[2]*fEdge[6]-8.647294131737468*fSkin[3]*coeff[6]-7.075058835057925*fEdge[3]*coeff[6]-4.538652358836817*fSkin[2]*coeff[6]+4.538652358836817*fEdge[2]*coeff[6]+4.059494080239555*coeff[5]*fSkin[5]-4.059494080239555*coeff[5]*fEdge[5]-7.916013456467132*coeff[4]*fSkin[4]-2.201129415351355*coeff[1]*fSkin[4]+7.080297679785432*coeff[0]*fSkin[4]+5.074367600299444*coeff[4]*fEdge[4]-2.201129415351355*coeff[1]*fEdge[4]-4.538652358836817*coeff[0]*fEdge[4]-8.64729413173747*fSkin[1]*coeff[4]-7.07505883505793*fEdge[1]*coeff[4]-4.538652358836817*fSkin[0]*coeff[4]+4.538652358836817*fEdge[0]*coeff[4]-1.217848224071866*coeff[3]*fSkin[3]+7.734375*coeff[2]*fSkin[3]+1.217848224071866*coeff[3]*fEdge[3]+6.328125*coeff[2]*fEdge[3]+4.059494080239555*coeff[2]*fSkin[2]-4.059494080239555*coeff[2]*fEdge[2]-1.217848224071866*coeff[1]*fSkin[1]+7.734375*coeff[0]*fSkin[1]+1.217848224071866*coeff[1]*fEdge[1]+6.328125*coeff[0]*fEdge[1]+4.059494080239555*coeff[0]*fSkin[0]-4.059494080239555*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = (-4.059494080239555*coeff[6]*fSkin[7])-0.6288941186718158*coeff[3]*fSkin[7]+3.630921887069454*coeff[2]*fSkin[7]-4.059494080239555*coeff[6]*fEdge[7]+0.6288941186718158*coeff[3]*fEdge[7]+3.630921887069454*coeff[2]*fEdge[7]-1.136658342467076*fSkin[6]*coeff[7]-1.136658342467076*fEdge[6]*coeff[7]-0.6288941186718158*fSkin[3]*coeff[7]+0.6288941186718158*fEdge[3]*coeff[7]+3.0*coeff[5]*fSkin[6]-3.75*coeff[4]*fSkin[6]-1.270822660474309*coeff[1]*fSkin[6]+3.354101966249684*coeff[0]*fSkin[6]-3.0*coeff[5]*fEdge[6]+3.75*coeff[4]*fEdge[6]-1.270822660474309*coeff[1]*fEdge[6]-3.354101966249684*coeff[0]*fEdge[6]-2.34375*fSkin[5]*coeff[6]+2.34375*fEdge[5]*coeff[6]-3.75*fSkin[4]*coeff[6]+3.75*fEdge[4]*coeff[6]-4.538652358836817*fSkin[1]*coeff[6]-4.538652358836817*fEdge[1]*coeff[6]-2.620392161132566*fSkin[0]*coeff[6]+2.620392161132566*fEdge[0]*coeff[6]+2.096313728906053*coeff[2]*fSkin[5]-2.096313728906053*coeff[2]*fEdge[5]+3.630921887069454*fSkin[3]*coeff[5]+3.630921887069454*fEdge[3]*coeff[5]+2.096313728906053*fSkin[2]*coeff[5]-2.096313728906053*fEdge[2]*coeff[5]-1.270822660474308*coeff[3]*fSkin[4]+3.354101966249685*coeff[2]*fSkin[4]-1.270822660474308*coeff[3]*fEdge[4]-3.354101966249685*coeff[2]*fEdge[4]-4.538652358836817*fSkin[3]*coeff[4]-4.538652358836817*fEdge[3]*coeff[4]-2.620392161132567*fSkin[2]*coeff[4]+2.620392161132567*fEdge[2]*coeff[4]-0.703125*coeff[1]*fSkin[3]+4.059494080239555*coeff[0]*fSkin[3]+0.703125*coeff[1]*fEdge[3]+4.059494080239555*coeff[0]*fEdge[3]-0.703125*fSkin[1]*coeff[3]+0.703125*fEdge[1]*coeff[3]+2.34375*coeff[0]*fSkin[2]-2.34375*coeff[0]*fEdge[2]+4.059494080239555*fSkin[1]*coeff[2]+4.059494080239555*fEdge[1]*coeff[2]+2.34375*fSkin[0]*coeff[2]-2.34375*fEdge[0]*coeff[2]; 
  edgeSurf_incr[3] = (-7.734375*coeff[6]*fSkin[7])-1.089276566120836*coeff[3]*fSkin[7]+6.917835305389974*coeff[2]*fSkin[7]-6.328125*coeff[6]*fEdge[7]+1.089276566120836*coeff[3]*fEdge[7]+5.660047068046341*coeff[2]*fEdge[7]-1.96875*fSkin[6]*coeff[7]-1.96875*fEdge[6]*coeff[7]-1.089276566120836*fSkin[3]*coeff[7]+1.089276566120836*fEdge[3]*coeff[7]+6.332810765173708*coeff[5]*fSkin[6]-7.916013456467136*coeff[4]*fSkin[6]-2.201129415351355*coeff[1]*fSkin[6]+7.080297679785435*coeff[0]*fSkin[6]-4.059494080239557*coeff[5]*fEdge[6]+5.074367600299446*coeff[4]*fEdge[6]-2.201129415351355*coeff[1]*fEdge[6]-4.538652358836817*coeff[0]*fEdge[6]-4.059494080239557*fSkin[5]*coeff[6]+4.059494080239557*fEdge[5]*coeff[6]-7.916013456467136*fSkin[4]*coeff[6]+5.074367600299446*fEdge[4]*coeff[6]-8.647294131737468*fSkin[1]*coeff[6]-7.075058835057925*fEdge[1]*coeff[6]-4.538652358836817*fSkin[0]*coeff[6]+4.538652358836817*fEdge[0]*coeff[6]+3.630921887069454*coeff[2]*fSkin[5]-3.630921887069454*coeff[2]*fEdge[5]+6.917835305389976*fSkin[3]*coeff[5]+5.660047068046343*fEdge[3]*coeff[5]+3.630921887069454*fSkin[2]*coeff[5]-3.630921887069454*fEdge[2]*coeff[5]-2.201129415351355*coeff[3]*fSkin[4]+7.080297679785432*coeff[2]*fSkin[4]-2.201129415351355*coeff[3]*fEdge[4]-4.538652358836817*coeff[2]*fEdge[4]-8.64729413173747*fSkin[3]*coeff[4]-7.07505883505793*fEdge[3]*coeff[4]-4.538652358836817*fSkin[2]*coeff[4]+4.538652358836817*fEdge[2]*coeff[4]-1.217848224071866*coeff[1]*fSkin[3]+7.734375*coeff[0]*fSkin[3]+1.217848224071866*coeff[1]*fEdge[3]+6.328125*coeff[0]*fEdge[3]-1.217848224071866*fSkin[1]*coeff[3]+1.217848224071866*fEdge[1]*coeff[3]+4.059494080239555*coeff[0]*fSkin[2]-4.059494080239555*coeff[0]*fEdge[2]+7.734375*fSkin[1]*coeff[2]+6.328125*fEdge[1]*coeff[2]+4.059494080239555*fSkin[0]*coeff[2]-4.059494080239555*fEdge[0]*coeff[2]; 
  edgeSurf_incr[4] = (-10.53397648775292*coeff[7]*fSkin[7])+7.806482057199325*coeff[5]*fSkin[7]-3.301694123027033*coeff[7]*fEdge[7]+2.360099226595145*coeff[5]*fEdge[7]-6.898751585431963*fSkin[5]*coeff[7]+3.9940140757764*fEdge[5]*coeff[7]-24.42205494175552*coeff[6]*fSkin[6]-8.159583101281509*coeff[3]*fSkin[6]+10.171875*coeff[2]*fSkin[6]+2.41076078824196*coeff[6]*fEdge[6]-0.3653544672215601*coeff[3]*fEdge[6]-0.3281250000000027*coeff[2]*fEdge[6]-30.04025619377272*fSkin[3]*coeff[6]-8.118988160479114*fEdge[3]*coeff[6]-18.515625*fSkin[2]*coeff[6]+7.265625000000002*fEdge[2]*coeff[6]+2.096313728906054*coeff[5]*fSkin[5]-2.096313728906054*coeff[5]*fEdge[5]-24.42205494175552*coeff[4]*fSkin[4]-8.159583101281505*coeff[1]*fSkin[4]+10.171875*coeff[0]*fSkin[4]+2.41076078824196*coeff[4]*fEdge[4]-0.3653544672215599*coeff[1]*fEdge[4]-0.328125*coeff[0]*fEdge[4]-30.04025619377271*fSkin[1]*coeff[4]-8.118988160479118*fEdge[1]*coeff[4]-18.515625*fSkin[0]*coeff[4]+7.265625*fEdge[0]*coeff[4]-10.53397648775292*coeff[3]*fSkin[3]+7.806482057199325*coeff[2]*fSkin[3]-3.301694123027033*coeff[3]*fEdge[3]+2.360099226595146*coeff[2]*fEdge[3]-6.898751585431961*fSkin[2]*coeff[3]+3.994014075776398*fEdge[2]*coeff[3]+2.096313728906054*coeff[2]*fSkin[2]-2.096313728906054*coeff[2]*fEdge[2]-10.53397648775292*coeff[1]*fSkin[1]+7.806482057199325*coeff[0]*fSkin[1]-3.301694123027033*coeff[1]*fEdge[1]+2.360099226595146*coeff[0]*fEdge[1]-6.898751585431961*fSkin[0]*coeff[1]+3.994014075776398*fEdge[0]*coeff[1]+2.096313728906054*coeff[0]*fSkin[0]-2.096313728906054*coeff[0]*fEdge[0]; 
  edgeSurf_incr[5] = (-0.4492100847655829*coeff[7]*fSkin[7])+2.593515633621038*coeff[5]*fSkin[7]-4.538652358836817*coeff[4]*fSkin[7]-0.703125*coeff[1]*fSkin[7]+4.059494080239557*coeff[0]*fSkin[7]+0.4492100847655829*coeff[7]*fEdge[7]+2.593515633621038*coeff[5]*fEdge[7]-4.538652358836817*coeff[4]*fEdge[7]+0.703125*coeff[1]*fEdge[7]+4.059494080239557*coeff[0]*fEdge[7]-1.270822660474309*fSkin[4]*coeff[7]-1.270822660474309*fEdge[4]*coeff[7]-0.703125*fSkin[1]*coeff[7]+0.703125*fEdge[1]*coeff[7]-3.354101966249685*coeff[6]*fSkin[6]-1.136658342467076*coeff[3]*fSkin[6]+3.0*coeff[2]*fSkin[6]+3.354101966249685*coeff[6]*fEdge[6]-1.136658342467076*coeff[3]*fEdge[6]-3.0*coeff[2]*fEdge[6]-4.059494080239557*fSkin[3]*coeff[6]-4.059494080239557*fEdge[3]*coeff[6]-2.34375*fSkin[2]*coeff[6]+2.34375*fEdge[2]*coeff[6]+1.497366949218609*coeff[5]*fSkin[5]-2.620392161132567*coeff[4]*fSkin[5]+2.34375*coeff[0]*fSkin[5]-1.497366949218609*coeff[5]*fEdge[5]+2.620392161132567*coeff[4]*fEdge[5]-2.34375*coeff[0]*fEdge[5]+3.354101966249685*fSkin[4]*coeff[5]-3.354101966249685*fEdge[4]*coeff[5]+4.059494080239555*fSkin[1]*coeff[5]+4.059494080239555*fEdge[1]*coeff[5]+2.34375*fSkin[0]*coeff[5]-2.34375*fEdge[0]*coeff[5]-0.6288941186718159*coeff[3]*fSkin[3]+3.630921887069454*coeff[2]*fSkin[3]+0.6288941186718159*coeff[3]*fEdge[3]+3.630921887069454*coeff[2]*fEdge[3]+2.096313728906053*coeff[2]*fSkin[2]-2.096313728906053*coeff[2]*fEdge[2]; 
  edgeSurf_incr[6] = (-26.86882196431396*coeff[6]*fSkin[7])-9.421875*coeff[3]*fSkin[7]+6.982329818012035*coeff[2]*fSkin[7]-7.261843774138907*coeff[6]*fEdge[7]-2.953125*coeff[3]*fEdge[7]+2.110936921724569*coeff[2]*fEdge[7]-7.298152993009602*fSkin[6]*coeff[7]-0.3267829698362508*fEdge[6]*coeff[7]-9.421875*fSkin[3]*coeff[7]-2.953125*fEdge[3]*coeff[7]-6.170431001964124*fSkin[2]*coeff[7]+3.572354790610808*fEdge[2]*coeff[7]+9.09800158345227*coeff[5]*fSkin[6]-24.42205494175552*coeff[4]*fSkin[6]-8.159583101281505*coeff[1]*fSkin[6]+10.171875*coeff[0]*fSkin[6]-0.2934839220468475*coeff[5]*fEdge[6]+2.41076078824196*coeff[4]*fEdge[6]-0.3653544672215599*coeff[1]*fEdge[6]-0.328125*coeff[0]*fEdge[6]-16.56087845835782*fSkin[5]*coeff[6]+6.498572559608766*fEdge[5]*coeff[6]-24.42205494175552*fSkin[4]*coeff[6]+2.41076078824196*fEdge[4]*coeff[6]-30.04025619377271*fSkin[1]*coeff[6]-8.118988160479118*fEdge[1]*coeff[6]-18.515625*fSkin[0]*coeff[6]+7.265625*fEdge[0]*coeff[6]-6.170431001964126*coeff[3]*fSkin[5]+1.875000000000001*coeff[2]*fSkin[5]+3.57235479061081*coeff[3]*fEdge[5]-1.875000000000001*coeff[2]*fEdge[5]+6.982329818012039*fSkin[3]*coeff[5]+2.110936921724571*fEdge[3]*coeff[5]+1.875000000000001*fSkin[2]*coeff[5]-1.875000000000001*fEdge[2]*coeff[5]-8.159583101281509*coeff[3]*fSkin[4]+10.171875*coeff[2]*fSkin[4]-0.3653544672215601*coeff[3]*fEdge[4]-0.3281250000000027*coeff[2]*fEdge[4]-30.04025619377272*fSkin[3]*coeff[4]-8.118988160479114*fEdge[3]*coeff[4]-18.515625*fSkin[2]*coeff[4]+7.265625000000002*fEdge[2]*coeff[4]-10.53397648775291*coeff[1]*fSkin[3]+7.806482057199325*coeff[0]*fSkin[3]-3.301694123027032*coeff[1]*fEdge[3]+2.360099226595145*coeff[0]*fEdge[3]-10.53397648775291*fSkin[1]*coeff[3]-3.301694123027032*fEdge[1]*coeff[3]-6.898751585431963*fSkin[0]*coeff[3]+3.9940140757764*fEdge[0]*coeff[3]-6.898751585431963*coeff[1]*fSkin[2]+2.096313728906054*coeff[0]*fSkin[2]+3.9940140757764*coeff[1]*fEdge[2]-2.096313728906054*coeff[0]*fEdge[2]+7.806482057199325*fSkin[1]*coeff[2]+2.360099226595145*fEdge[1]*coeff[2]+2.096313728906054*fSkin[0]*coeff[2]-2.096313728906054*fEdge[0]*coeff[2]; 
  edgeSurf_incr[7] = (-0.7780546900863115*coeff[7]*fSkin[7])+4.941310932421412*coeff[5]*fSkin[7]-8.64729413173747*coeff[4]*fSkin[7]-1.217848224071866*coeff[1]*fSkin[7]+7.734375*coeff[0]*fSkin[7]+0.7780546900863115*coeff[7]*fEdge[7]+4.042890762890246*coeff[5]*fEdge[7]-7.07505883505793*coeff[4]*fEdge[7]+1.217848224071866*coeff[1]*fEdge[7]+6.328125*coeff[0]*fEdge[7]-2.201129415351355*fSkin[4]*coeff[7]-2.201129415351355*fEdge[4]*coeff[7]-1.217848224071866*fSkin[1]*coeff[7]+1.217848224071866*fEdge[1]*coeff[7]-7.080297679785434*coeff[6]*fSkin[6]-1.96875*coeff[3]*fSkin[6]+6.332810765173706*coeff[2]*fSkin[6]+4.538652358836817*coeff[6]*fEdge[6]-1.96875*coeff[3]*fEdge[6]-4.059494080239555*coeff[2]*fEdge[6]-7.734375*fSkin[3]*coeff[6]-6.328125*fEdge[3]*coeff[6]-4.059494080239555*fSkin[2]*coeff[6]+4.059494080239555*fEdge[2]*coeff[6]+2.593515633621038*coeff[5]*fSkin[5]-4.538652358836817*coeff[4]*fSkin[5]+4.059494080239557*coeff[0]*fSkin[5]-2.593515633621038*coeff[5]*fEdge[5]+4.538652358836817*coeff[4]*fEdge[5]-4.059494080239557*coeff[0]*fEdge[5]+7.080297679785434*fSkin[4]*coeff[5]-4.538652358836817*fEdge[4]*coeff[5]+7.734375000000001*fSkin[1]*coeff[5]+6.328125000000001*fEdge[1]*coeff[5]+4.059494080239557*fSkin[0]*coeff[5]-4.059494080239557*fEdge[0]*coeff[5]-1.089276566120836*coeff[3]*fSkin[3]+6.917835305389974*coeff[2]*fSkin[3]+1.089276566120836*coeff[3]*fEdge[3]+5.660047068046341*coeff[2]*fEdge[3]+3.630921887069454*coeff[2]*fSkin[2]-3.630921887069454*coeff[2]*fEdge[2]; 

  boundSurf_incr[0] = (-1.40625*coeff[7]*fSkin[7])+2.541645320948617*coeff[3]*fSkin[6]+2.541645320948617*coeff[1]*fSkin[4]-1.40625*coeff[3]*fSkin[3]-1.40625*coeff[1]*fSkin[1]; 
  boundSurf_incr[1] = 2.435696448143733*coeff[7]*fSkin[7]+1.40625*coeff[5]*fSkin[7]+2.841645856167689*coeff[6]*fSkin[6]-4.40225883070271*coeff[3]*fSkin[6]-2.541645320948617*coeff[2]*fSkin[6]-1.572235296679539*fSkin[3]*coeff[6]+2.841645856167689*coeff[4]*fSkin[4]-4.402258830702711*coeff[1]*fSkin[4]-2.541645320948617*coeff[0]*fSkin[4]-1.57223529667954*fSkin[1]*coeff[4]+2.435696448143733*coeff[3]*fSkin[3]+1.40625*coeff[2]*fSkin[3]+2.435696448143733*coeff[1]*fSkin[1]+1.40625*coeff[0]*fSkin[1]; 
  boundSurf_incr[2] = (-1.257788237343632*coeff[3]*fSkin[7])+2.273316684934151*fSkin[6]*coeff[7]-1.257788237343632*fSkin[3]*coeff[7]+2.541645320948617*coeff[1]*fSkin[6]+2.541645320948617*coeff[3]*fSkin[4]-1.40625*coeff[1]*fSkin[3]-1.40625*fSkin[1]*coeff[3]; 
  boundSurf_incr[3] = (-1.40625*coeff[6]*fSkin[7])+2.178553132241672*coeff[3]*fSkin[7]+1.257788237343632*coeff[2]*fSkin[7]-3.9375*fSkin[6]*coeff[7]+2.178553132241672*fSkin[3]*coeff[7]-2.273316684934152*coeff[5]*fSkin[6]+2.841645856167689*coeff[4]*fSkin[6]-4.40225883070271*coeff[1]*fSkin[6]-2.541645320948617*coeff[0]*fSkin[6]+2.841645856167689*fSkin[4]*coeff[6]-1.572235296679539*fSkin[1]*coeff[6]+1.257788237343632*fSkin[3]*coeff[5]-4.402258830702711*coeff[3]*fSkin[4]-2.541645320948617*coeff[2]*fSkin[4]-1.57223529667954*fSkin[3]*coeff[4]+2.435696448143733*coeff[1]*fSkin[3]+1.40625*coeff[0]*fSkin[3]+2.435696448143733*fSkin[1]*coeff[3]+1.40625*fSkin[1]*coeff[2]; 
  boundSurf_incr[4] = (-7.232282364725883*coeff[7]*fSkin[7])-5.446382830604181*coeff[5]*fSkin[7]+2.904737509655563*fSkin[5]*coeff[7]-22.01129415351356*coeff[6]*fSkin[6]+8.524937568503068*coeff[3]*fSkin[6]+9.843749999999998*coeff[2]*fSkin[6]+21.92126803329361*fSkin[3]*coeff[6]-11.25*fSkin[2]*coeff[6]-22.01129415351356*coeff[4]*fSkin[4]+8.524937568503065*coeff[1]*fSkin[4]+9.84375*coeff[0]*fSkin[4]+21.92126803329359*fSkin[1]*coeff[4]-11.25*fSkin[0]*coeff[4]-7.232282364725883*coeff[3]*fSkin[3]-5.446382830604179*coeff[2]*fSkin[3]+2.904737509655563*fSkin[2]*coeff[3]-7.232282364725883*coeff[1]*fSkin[1]-5.446382830604179*coeff[0]*fSkin[1]+2.904737509655563*fSkin[0]*coeff[1]; 
  boundSurf_incr[5] = (-0.8984201695311658*coeff[7]*fSkin[7])-1.40625*coeff[1]*fSkin[7]+2.541645320948617*fSkin[4]*coeff[7]-1.40625*fSkin[1]*coeff[7]+2.273316684934152*coeff[3]*fSkin[6]-1.257788237343632*coeff[3]*fSkin[3]; 
  boundSurf_incr[6] = 19.60697819017505*coeff[6]*fSkin[7]-6.46875*coeff[3]*fSkin[7]-4.871392896287466*coeff[2]*fSkin[7]+7.624935962845852*fSkin[6]*coeff[7]-6.46875*fSkin[3]*coeff[7]+2.598076211353316*fSkin[2]*coeff[7]+8.804517661405422*coeff[5]*fSkin[6]-22.01129415351356*coeff[4]*fSkin[6]+8.524937568503065*coeff[1]*fSkin[6]+9.84375*coeff[0]*fSkin[6]-10.06230589874905*fSkin[5]*coeff[6]-22.01129415351356*fSkin[4]*coeff[6]+21.92126803329359*fSkin[1]*coeff[6]-11.25*fSkin[0]*coeff[6]+2.598076211353316*coeff[3]*fSkin[5]-4.871392896287468*fSkin[3]*coeff[5]+8.524937568503068*coeff[3]*fSkin[4]+9.843749999999998*coeff[2]*fSkin[4]+21.92126803329361*fSkin[3]*coeff[4]-11.25*fSkin[2]*coeff[4]-7.232282364725882*coeff[1]*fSkin[3]-5.446382830604181*coeff[0]*fSkin[3]-7.232282364725882*fSkin[1]*coeff[3]+2.904737509655563*fSkin[0]*coeff[3]+2.904737509655563*coeff[1]*fSkin[2]-5.446382830604181*fSkin[1]*coeff[2]; 
  boundSurf_incr[7] = 1.556109380172623*coeff[7]*fSkin[7]+0.8984201695311658*coeff[5]*fSkin[7]-1.57223529667954*coeff[4]*fSkin[7]+2.435696448143733*coeff[1]*fSkin[7]+1.40625*coeff[0]*fSkin[7]-4.402258830702711*fSkin[4]*coeff[7]+2.435696448143733*fSkin[1]*coeff[7]+2.541645320948617*coeff[6]*fSkin[6]-3.9375*coeff[3]*fSkin[6]-2.273316684934151*coeff[2]*fSkin[6]-1.40625*fSkin[3]*coeff[6]-2.541645320948617*fSkin[4]*coeff[5]+1.40625*fSkin[1]*coeff[5]+2.178553132241672*coeff[3]*fSkin[3]+1.257788237343632*coeff[2]*fSkin[3]; 

  } else { 

  edgeSurf_incr[0] = (-0.703125*coeff[7]*fSkin[7])-4.059494080239557*coeff[5]*fSkin[7]+0.703125*coeff[7]*fEdge[7]-4.059494080239557*coeff[5]*fEdge[7]-3.75*coeff[6]*fSkin[6]+1.270822660474309*coeff[3]*fSkin[6]+3.354101966249684*coeff[2]*fSkin[6]+3.75*coeff[6]*fEdge[6]+1.270822660474309*coeff[3]*fEdge[6]-3.354101966249684*coeff[2]*fEdge[6]+4.538652358836817*fSkin[3]*coeff[6]+4.538652358836817*fEdge[3]*coeff[6]-2.620392161132566*fSkin[2]*coeff[6]+2.620392161132566*fEdge[2]*coeff[6]+2.34375*coeff[5]*fSkin[5]-2.34375*coeff[5]*fEdge[5]-3.75*coeff[4]*fSkin[4]+1.270822660474308*coeff[1]*fSkin[4]+3.354101966249685*coeff[0]*fSkin[4]+3.75*coeff[4]*fEdge[4]+1.270822660474308*coeff[1]*fEdge[4]-3.354101966249685*coeff[0]*fEdge[4]+4.538652358836817*fSkin[1]*coeff[4]+4.538652358836817*fEdge[1]*coeff[4]-2.620392161132567*fSkin[0]*coeff[4]+2.620392161132567*fEdge[0]*coeff[4]-0.703125*coeff[3]*fSkin[3]-4.059494080239555*coeff[2]*fSkin[3]+0.703125*coeff[3]*fEdge[3]-4.059494080239555*coeff[2]*fEdge[3]+2.34375*coeff[2]*fSkin[2]-2.34375*coeff[2]*fEdge[2]-0.703125*coeff[1]*fSkin[1]-4.059494080239555*coeff[0]*fSkin[1]+0.703125*coeff[1]*fEdge[1]-4.059494080239555*coeff[0]*fEdge[1]+2.34375*coeff[0]*fSkin[0]-2.34375*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 1.217848224071866*coeff[7]*fSkin[7]+7.734375000000001*coeff[5]*fSkin[7]-1.217848224071866*coeff[7]*fEdge[7]+6.328125000000001*coeff[5]*fEdge[7]+7.916013456467132*coeff[6]*fSkin[6]-2.201129415351355*coeff[3]*fSkin[6]-7.080297679785434*coeff[2]*fSkin[6]-5.074367600299444*coeff[6]*fEdge[6]-2.201129415351355*coeff[3]*fEdge[6]+4.538652358836817*coeff[2]*fEdge[6]-8.647294131737468*fSkin[3]*coeff[6]-7.075058835057925*fEdge[3]*coeff[6]+4.538652358836817*fSkin[2]*coeff[6]-4.538652358836817*fEdge[2]*coeff[6]-4.059494080239555*coeff[5]*fSkin[5]+4.059494080239555*coeff[5]*fEdge[5]+7.916013456467132*coeff[4]*fSkin[4]-2.201129415351355*coeff[1]*fSkin[4]-7.080297679785432*coeff[0]*fSkin[4]-5.074367600299444*coeff[4]*fEdge[4]-2.201129415351355*coeff[1]*fEdge[4]+4.538652358836817*coeff[0]*fEdge[4]-8.64729413173747*fSkin[1]*coeff[4]-7.07505883505793*fEdge[1]*coeff[4]+4.538652358836817*fSkin[0]*coeff[4]-4.538652358836817*fEdge[0]*coeff[4]+1.217848224071866*coeff[3]*fSkin[3]+7.734375*coeff[2]*fSkin[3]-1.217848224071866*coeff[3]*fEdge[3]+6.328125*coeff[2]*fEdge[3]-4.059494080239555*coeff[2]*fSkin[2]+4.059494080239555*coeff[2]*fEdge[2]+1.217848224071866*coeff[1]*fSkin[1]+7.734375*coeff[0]*fSkin[1]-1.217848224071866*coeff[1]*fEdge[1]+6.328125*coeff[0]*fEdge[1]-4.059494080239555*coeff[0]*fSkin[0]+4.059494080239555*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 4.059494080239555*coeff[6]*fSkin[7]-0.6288941186718158*coeff[3]*fSkin[7]-3.630921887069454*coeff[2]*fSkin[7]+4.059494080239555*coeff[6]*fEdge[7]+0.6288941186718158*coeff[3]*fEdge[7]-3.630921887069454*coeff[2]*fEdge[7]+1.136658342467076*fSkin[6]*coeff[7]+1.136658342467076*fEdge[6]*coeff[7]-0.6288941186718158*fSkin[3]*coeff[7]+0.6288941186718158*fEdge[3]*coeff[7]+3.0*coeff[5]*fSkin[6]-3.75*coeff[4]*fSkin[6]+1.270822660474309*coeff[1]*fSkin[6]+3.354101966249684*coeff[0]*fSkin[6]-3.0*coeff[5]*fEdge[6]+3.75*coeff[4]*fEdge[6]+1.270822660474309*coeff[1]*fEdge[6]-3.354101966249684*coeff[0]*fEdge[6]-2.34375*fSkin[5]*coeff[6]+2.34375*fEdge[5]*coeff[6]-3.75*fSkin[4]*coeff[6]+3.75*fEdge[4]*coeff[6]+4.538652358836817*fSkin[1]*coeff[6]+4.538652358836817*fEdge[1]*coeff[6]-2.620392161132566*fSkin[0]*coeff[6]+2.620392161132566*fEdge[0]*coeff[6]+2.096313728906053*coeff[2]*fSkin[5]-2.096313728906053*coeff[2]*fEdge[5]-3.630921887069454*fSkin[3]*coeff[5]-3.630921887069454*fEdge[3]*coeff[5]+2.096313728906053*fSkin[2]*coeff[5]-2.096313728906053*fEdge[2]*coeff[5]+1.270822660474308*coeff[3]*fSkin[4]+3.354101966249685*coeff[2]*fSkin[4]+1.270822660474308*coeff[3]*fEdge[4]-3.354101966249685*coeff[2]*fEdge[4]+4.538652358836817*fSkin[3]*coeff[4]+4.538652358836817*fEdge[3]*coeff[4]-2.620392161132567*fSkin[2]*coeff[4]+2.620392161132567*fEdge[2]*coeff[4]-0.703125*coeff[1]*fSkin[3]-4.059494080239555*coeff[0]*fSkin[3]+0.703125*coeff[1]*fEdge[3]-4.059494080239555*coeff[0]*fEdge[3]-0.703125*fSkin[1]*coeff[3]+0.703125*fEdge[1]*coeff[3]+2.34375*coeff[0]*fSkin[2]-2.34375*coeff[0]*fEdge[2]-4.059494080239555*fSkin[1]*coeff[2]-4.059494080239555*fEdge[1]*coeff[2]+2.34375*fSkin[0]*coeff[2]-2.34375*fEdge[0]*coeff[2]; 
  edgeSurf_incr[3] = (-7.734375*coeff[6]*fSkin[7])+1.089276566120836*coeff[3]*fSkin[7]+6.917835305389974*coeff[2]*fSkin[7]-6.328125*coeff[6]*fEdge[7]-1.089276566120836*coeff[3]*fEdge[7]+5.660047068046341*coeff[2]*fEdge[7]-1.96875*fSkin[6]*coeff[7]-1.96875*fEdge[6]*coeff[7]+1.089276566120836*fSkin[3]*coeff[7]-1.089276566120836*fEdge[3]*coeff[7]-6.332810765173708*coeff[5]*fSkin[6]+7.916013456467136*coeff[4]*fSkin[6]-2.201129415351355*coeff[1]*fSkin[6]-7.080297679785435*coeff[0]*fSkin[6]+4.059494080239557*coeff[5]*fEdge[6]-5.074367600299446*coeff[4]*fEdge[6]-2.201129415351355*coeff[1]*fEdge[6]+4.538652358836817*coeff[0]*fEdge[6]+4.059494080239557*fSkin[5]*coeff[6]-4.059494080239557*fEdge[5]*coeff[6]+7.916013456467136*fSkin[4]*coeff[6]-5.074367600299446*fEdge[4]*coeff[6]-8.647294131737468*fSkin[1]*coeff[6]-7.075058835057925*fEdge[1]*coeff[6]+4.538652358836817*fSkin[0]*coeff[6]-4.538652358836817*fEdge[0]*coeff[6]-3.630921887069454*coeff[2]*fSkin[5]+3.630921887069454*coeff[2]*fEdge[5]+6.917835305389976*fSkin[3]*coeff[5]+5.660047068046343*fEdge[3]*coeff[5]-3.630921887069454*fSkin[2]*coeff[5]+3.630921887069454*fEdge[2]*coeff[5]-2.201129415351355*coeff[3]*fSkin[4]-7.080297679785432*coeff[2]*fSkin[4]-2.201129415351355*coeff[3]*fEdge[4]+4.538652358836817*coeff[2]*fEdge[4]-8.64729413173747*fSkin[3]*coeff[4]-7.07505883505793*fEdge[3]*coeff[4]+4.538652358836817*fSkin[2]*coeff[4]-4.538652358836817*fEdge[2]*coeff[4]+1.217848224071866*coeff[1]*fSkin[3]+7.734375*coeff[0]*fSkin[3]-1.217848224071866*coeff[1]*fEdge[3]+6.328125*coeff[0]*fEdge[3]+1.217848224071866*fSkin[1]*coeff[3]-1.217848224071866*fEdge[1]*coeff[3]-4.059494080239555*coeff[0]*fSkin[2]+4.059494080239555*coeff[0]*fEdge[2]+7.734375*fSkin[1]*coeff[2]+6.328125*fEdge[1]*coeff[2]-4.059494080239555*fSkin[0]*coeff[2]+4.059494080239555*fEdge[0]*coeff[2]; 
  edgeSurf_incr[4] = (-10.53397648775292*coeff[7]*fSkin[7])-7.806482057199325*coeff[5]*fSkin[7]-3.301694123027033*coeff[7]*fEdge[7]-2.360099226595145*coeff[5]*fEdge[7]+6.898751585431963*fSkin[5]*coeff[7]-3.9940140757764*fEdge[5]*coeff[7]-24.42205494175552*coeff[6]*fSkin[6]+8.159583101281509*coeff[3]*fSkin[6]+10.171875*coeff[2]*fSkin[6]+2.41076078824196*coeff[6]*fEdge[6]+0.3653544672215601*coeff[3]*fEdge[6]-0.3281250000000027*coeff[2]*fEdge[6]+30.04025619377272*fSkin[3]*coeff[6]+8.118988160479114*fEdge[3]*coeff[6]-18.515625*fSkin[2]*coeff[6]+7.265625000000002*fEdge[2]*coeff[6]+2.096313728906054*coeff[5]*fSkin[5]-2.096313728906054*coeff[5]*fEdge[5]-24.42205494175552*coeff[4]*fSkin[4]+8.159583101281505*coeff[1]*fSkin[4]+10.171875*coeff[0]*fSkin[4]+2.41076078824196*coeff[4]*fEdge[4]+0.3653544672215599*coeff[1]*fEdge[4]-0.328125*coeff[0]*fEdge[4]+30.04025619377271*fSkin[1]*coeff[4]+8.118988160479118*fEdge[1]*coeff[4]-18.515625*fSkin[0]*coeff[4]+7.265625*fEdge[0]*coeff[4]-10.53397648775292*coeff[3]*fSkin[3]-7.806482057199325*coeff[2]*fSkin[3]-3.301694123027033*coeff[3]*fEdge[3]-2.360099226595146*coeff[2]*fEdge[3]+6.898751585431961*fSkin[2]*coeff[3]-3.994014075776398*fEdge[2]*coeff[3]+2.096313728906054*coeff[2]*fSkin[2]-2.096313728906054*coeff[2]*fEdge[2]-10.53397648775292*coeff[1]*fSkin[1]-7.806482057199325*coeff[0]*fSkin[1]-3.301694123027033*coeff[1]*fEdge[1]-2.360099226595146*coeff[0]*fEdge[1]+6.898751585431961*fSkin[0]*coeff[1]-3.994014075776398*fEdge[0]*coeff[1]+2.096313728906054*coeff[0]*fSkin[0]-2.096313728906054*coeff[0]*fEdge[0]; 
  edgeSurf_incr[5] = (-0.4492100847655829*coeff[7]*fSkin[7])-2.593515633621038*coeff[5]*fSkin[7]+4.538652358836817*coeff[4]*fSkin[7]-0.703125*coeff[1]*fSkin[7]-4.059494080239557*coeff[0]*fSkin[7]+0.4492100847655829*coeff[7]*fEdge[7]-2.593515633621038*coeff[5]*fEdge[7]+4.538652358836817*coeff[4]*fEdge[7]+0.703125*coeff[1]*fEdge[7]-4.059494080239557*coeff[0]*fEdge[7]+1.270822660474309*fSkin[4]*coeff[7]+1.270822660474309*fEdge[4]*coeff[7]-0.703125*fSkin[1]*coeff[7]+0.703125*fEdge[1]*coeff[7]-3.354101966249685*coeff[6]*fSkin[6]+1.136658342467076*coeff[3]*fSkin[6]+3.0*coeff[2]*fSkin[6]+3.354101966249685*coeff[6]*fEdge[6]+1.136658342467076*coeff[3]*fEdge[6]-3.0*coeff[2]*fEdge[6]+4.059494080239557*fSkin[3]*coeff[6]+4.059494080239557*fEdge[3]*coeff[6]-2.34375*fSkin[2]*coeff[6]+2.34375*fEdge[2]*coeff[6]+1.497366949218609*coeff[5]*fSkin[5]-2.620392161132567*coeff[4]*fSkin[5]+2.34375*coeff[0]*fSkin[5]-1.497366949218609*coeff[5]*fEdge[5]+2.620392161132567*coeff[4]*fEdge[5]-2.34375*coeff[0]*fEdge[5]+3.354101966249685*fSkin[4]*coeff[5]-3.354101966249685*fEdge[4]*coeff[5]-4.059494080239555*fSkin[1]*coeff[5]-4.059494080239555*fEdge[1]*coeff[5]+2.34375*fSkin[0]*coeff[5]-2.34375*fEdge[0]*coeff[5]-0.6288941186718159*coeff[3]*fSkin[3]-3.630921887069454*coeff[2]*fSkin[3]+0.6288941186718159*coeff[3]*fEdge[3]-3.630921887069454*coeff[2]*fEdge[3]+2.096313728906053*coeff[2]*fSkin[2]-2.096313728906053*coeff[2]*fEdge[2]; 
  edgeSurf_incr[6] = 26.86882196431396*coeff[6]*fSkin[7]-9.421875*coeff[3]*fSkin[7]-6.982329818012035*coeff[2]*fSkin[7]+7.261843774138907*coeff[6]*fEdge[7]-2.953125*coeff[3]*fEdge[7]-2.110936921724569*coeff[2]*fEdge[7]+7.298152993009602*fSkin[6]*coeff[7]+0.3267829698362508*fEdge[6]*coeff[7]-9.421875*fSkin[3]*coeff[7]-2.953125*fEdge[3]*coeff[7]+6.170431001964124*fSkin[2]*coeff[7]-3.572354790610808*fEdge[2]*coeff[7]+9.09800158345227*coeff[5]*fSkin[6]-24.42205494175552*coeff[4]*fSkin[6]+8.159583101281505*coeff[1]*fSkin[6]+10.171875*coeff[0]*fSkin[6]-0.2934839220468475*coeff[5]*fEdge[6]+2.41076078824196*coeff[4]*fEdge[6]+0.3653544672215599*coeff[1]*fEdge[6]-0.328125*coeff[0]*fEdge[6]-16.56087845835782*fSkin[5]*coeff[6]+6.498572559608766*fEdge[5]*coeff[6]-24.42205494175552*fSkin[4]*coeff[6]+2.41076078824196*fEdge[4]*coeff[6]+30.04025619377271*fSkin[1]*coeff[6]+8.118988160479118*fEdge[1]*coeff[6]-18.515625*fSkin[0]*coeff[6]+7.265625*fEdge[0]*coeff[6]+6.170431001964126*coeff[3]*fSkin[5]+1.875000000000001*coeff[2]*fSkin[5]-3.57235479061081*coeff[3]*fEdge[5]-1.875000000000001*coeff[2]*fEdge[5]-6.982329818012039*fSkin[3]*coeff[5]-2.110936921724571*fEdge[3]*coeff[5]+1.875000000000001*fSkin[2]*coeff[5]-1.875000000000001*fEdge[2]*coeff[5]+8.159583101281509*coeff[3]*fSkin[4]+10.171875*coeff[2]*fSkin[4]+0.3653544672215601*coeff[3]*fEdge[4]-0.3281250000000027*coeff[2]*fEdge[4]+30.04025619377272*fSkin[3]*coeff[4]+8.118988160479114*fEdge[3]*coeff[4]-18.515625*fSkin[2]*coeff[4]+7.265625000000002*fEdge[2]*coeff[4]-10.53397648775291*coeff[1]*fSkin[3]-7.806482057199325*coeff[0]*fSkin[3]-3.301694123027032*coeff[1]*fEdge[3]-2.360099226595145*coeff[0]*fEdge[3]-10.53397648775291*fSkin[1]*coeff[3]-3.301694123027032*fEdge[1]*coeff[3]+6.898751585431963*fSkin[0]*coeff[3]-3.9940140757764*fEdge[0]*coeff[3]+6.898751585431963*coeff[1]*fSkin[2]+2.096313728906054*coeff[0]*fSkin[2]-3.9940140757764*coeff[1]*fEdge[2]-2.096313728906054*coeff[0]*fEdge[2]-7.806482057199325*fSkin[1]*coeff[2]-2.360099226595145*fEdge[1]*coeff[2]+2.096313728906054*fSkin[0]*coeff[2]-2.096313728906054*fEdge[0]*coeff[2]; 
  edgeSurf_incr[7] = 0.7780546900863115*coeff[7]*fSkin[7]+4.941310932421412*coeff[5]*fSkin[7]-8.64729413173747*coeff[4]*fSkin[7]+1.217848224071866*coeff[1]*fSkin[7]+7.734375*coeff[0]*fSkin[7]-0.7780546900863115*coeff[7]*fEdge[7]+4.042890762890246*coeff[5]*fEdge[7]-7.07505883505793*coeff[4]*fEdge[7]-1.217848224071866*coeff[1]*fEdge[7]+6.328125*coeff[0]*fEdge[7]-2.201129415351355*fSkin[4]*coeff[7]-2.201129415351355*fEdge[4]*coeff[7]+1.217848224071866*fSkin[1]*coeff[7]-1.217848224071866*fEdge[1]*coeff[7]+7.080297679785434*coeff[6]*fSkin[6]-1.96875*coeff[3]*fSkin[6]-6.332810765173706*coeff[2]*fSkin[6]-4.538652358836817*coeff[6]*fEdge[6]-1.96875*coeff[3]*fEdge[6]+4.059494080239555*coeff[2]*fEdge[6]-7.734375*fSkin[3]*coeff[6]-6.328125*fEdge[3]*coeff[6]+4.059494080239555*fSkin[2]*coeff[6]-4.059494080239555*fEdge[2]*coeff[6]-2.593515633621038*coeff[5]*fSkin[5]+4.538652358836817*coeff[4]*fSkin[5]-4.059494080239557*coeff[0]*fSkin[5]+2.593515633621038*coeff[5]*fEdge[5]-4.538652358836817*coeff[4]*fEdge[5]+4.059494080239557*coeff[0]*fEdge[5]-7.080297679785434*fSkin[4]*coeff[5]+4.538652358836817*fEdge[4]*coeff[5]+7.734375000000001*fSkin[1]*coeff[5]+6.328125000000001*fEdge[1]*coeff[5]-4.059494080239557*fSkin[0]*coeff[5]+4.059494080239557*fEdge[0]*coeff[5]+1.089276566120836*coeff[3]*fSkin[3]+6.917835305389974*coeff[2]*fSkin[3]-1.089276566120836*coeff[3]*fEdge[3]+5.660047068046341*coeff[2]*fEdge[3]-3.630921887069454*coeff[2]*fSkin[2]+3.630921887069454*coeff[2]*fEdge[2]; 

  boundSurf_incr[0] = (-1.40625*coeff[7]*fSkin[7])-2.541645320948617*coeff[3]*fSkin[6]-2.541645320948617*coeff[1]*fSkin[4]-1.40625*coeff[3]*fSkin[3]-1.40625*coeff[1]*fSkin[1]; 
  boundSurf_incr[1] = (-2.435696448143733*coeff[7]*fSkin[7])+1.40625*coeff[5]*fSkin[7]-2.841645856167689*coeff[6]*fSkin[6]-4.40225883070271*coeff[3]*fSkin[6]+2.541645320948617*coeff[2]*fSkin[6]-1.572235296679539*fSkin[3]*coeff[6]-2.841645856167689*coeff[4]*fSkin[4]-4.402258830702711*coeff[1]*fSkin[4]+2.541645320948617*coeff[0]*fSkin[4]-1.57223529667954*fSkin[1]*coeff[4]-2.435696448143733*coeff[3]*fSkin[3]+1.40625*coeff[2]*fSkin[3]-2.435696448143733*coeff[1]*fSkin[1]+1.40625*coeff[0]*fSkin[1]; 
  boundSurf_incr[2] = (-1.257788237343632*coeff[3]*fSkin[7])-2.273316684934151*fSkin[6]*coeff[7]-1.257788237343632*fSkin[3]*coeff[7]-2.541645320948617*coeff[1]*fSkin[6]-2.541645320948617*coeff[3]*fSkin[4]-1.40625*coeff[1]*fSkin[3]-1.40625*fSkin[1]*coeff[3]; 
  boundSurf_incr[3] = (-1.40625*coeff[6]*fSkin[7])-2.178553132241672*coeff[3]*fSkin[7]+1.257788237343632*coeff[2]*fSkin[7]-3.9375*fSkin[6]*coeff[7]-2.178553132241672*fSkin[3]*coeff[7]+2.273316684934152*coeff[5]*fSkin[6]-2.841645856167689*coeff[4]*fSkin[6]-4.40225883070271*coeff[1]*fSkin[6]+2.541645320948617*coeff[0]*fSkin[6]-2.841645856167689*fSkin[4]*coeff[6]-1.572235296679539*fSkin[1]*coeff[6]+1.257788237343632*fSkin[3]*coeff[5]-4.402258830702711*coeff[3]*fSkin[4]+2.541645320948617*coeff[2]*fSkin[4]-1.57223529667954*fSkin[3]*coeff[4]-2.435696448143733*coeff[1]*fSkin[3]+1.40625*coeff[0]*fSkin[3]-2.435696448143733*fSkin[1]*coeff[3]+1.40625*fSkin[1]*coeff[2]; 
  boundSurf_incr[4] = (-7.232282364725883*coeff[7]*fSkin[7])+5.446382830604181*coeff[5]*fSkin[7]-2.904737509655563*fSkin[5]*coeff[7]-22.01129415351356*coeff[6]*fSkin[6]-8.524937568503068*coeff[3]*fSkin[6]+9.843749999999998*coeff[2]*fSkin[6]-21.92126803329361*fSkin[3]*coeff[6]-11.25*fSkin[2]*coeff[6]-22.01129415351356*coeff[4]*fSkin[4]-8.524937568503065*coeff[1]*fSkin[4]+9.84375*coeff[0]*fSkin[4]-21.92126803329359*fSkin[1]*coeff[4]-11.25*fSkin[0]*coeff[4]-7.232282364725883*coeff[3]*fSkin[3]+5.446382830604179*coeff[2]*fSkin[3]-2.904737509655563*fSkin[2]*coeff[3]-7.232282364725883*coeff[1]*fSkin[1]+5.446382830604179*coeff[0]*fSkin[1]-2.904737509655563*fSkin[0]*coeff[1]; 
  boundSurf_incr[5] = (-0.8984201695311658*coeff[7]*fSkin[7])-1.40625*coeff[1]*fSkin[7]-2.541645320948617*fSkin[4]*coeff[7]-1.40625*fSkin[1]*coeff[7]-2.273316684934152*coeff[3]*fSkin[6]-1.257788237343632*coeff[3]*fSkin[3]; 
  boundSurf_incr[6] = (-19.60697819017505*coeff[6]*fSkin[7])-6.46875*coeff[3]*fSkin[7]+4.871392896287466*coeff[2]*fSkin[7]-7.624935962845852*fSkin[6]*coeff[7]-6.46875*fSkin[3]*coeff[7]-2.598076211353316*fSkin[2]*coeff[7]+8.804517661405422*coeff[5]*fSkin[6]-22.01129415351356*coeff[4]*fSkin[6]-8.524937568503065*coeff[1]*fSkin[6]+9.84375*coeff[0]*fSkin[6]-10.06230589874905*fSkin[5]*coeff[6]-22.01129415351356*fSkin[4]*coeff[6]-21.92126803329359*fSkin[1]*coeff[6]-11.25*fSkin[0]*coeff[6]-2.598076211353316*coeff[3]*fSkin[5]+4.871392896287468*fSkin[3]*coeff[5]-8.524937568503068*coeff[3]*fSkin[4]+9.843749999999998*coeff[2]*fSkin[4]-21.92126803329361*fSkin[3]*coeff[4]-11.25*fSkin[2]*coeff[4]-7.232282364725882*coeff[1]*fSkin[3]+5.446382830604181*coeff[0]*fSkin[3]-7.232282364725882*fSkin[1]*coeff[3]-2.904737509655563*fSkin[0]*coeff[3]-2.904737509655563*coeff[1]*fSkin[2]+5.446382830604181*fSkin[1]*coeff[2]; 
  boundSurf_incr[7] = (-1.556109380172623*coeff[7]*fSkin[7])+0.8984201695311658*coeff[5]*fSkin[7]-1.57223529667954*coeff[4]*fSkin[7]-2.435696448143733*coeff[1]*fSkin[7]+1.40625*coeff[0]*fSkin[7]-4.402258830702711*fSkin[4]*coeff[7]-2.435696448143733*fSkin[1]*coeff[7]-2.541645320948617*coeff[6]*fSkin[6]-3.9375*coeff[3]*fSkin[6]+2.273316684934151*coeff[2]*fSkin[6]-1.40625*fSkin[3]*coeff[6]+2.541645320948617*fSkin[4]*coeff[5]+1.40625*fSkin[1]*coeff[5]-2.178553132241672*coeff[3]*fSkin[3]+1.257788237343632*coeff[2]*fSkin[3]; 

  out[0] += -1.0*(vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += -1.0*(vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 
  out[2] += -1.0*(vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac; 
  out[3] += -1.0*(vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*Jfac; 
  out[4] += -1.0*(vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*Jfac; 
  out[5] += -1.0*(vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*Jfac; 
  out[6] += -1.0*(vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*Jfac; 
  out[7] += -1.0*(vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*Jfac; 

  }

  return 0.;
}

