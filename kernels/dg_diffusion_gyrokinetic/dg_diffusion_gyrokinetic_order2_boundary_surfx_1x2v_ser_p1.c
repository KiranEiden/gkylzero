#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_boundary_surfx_1x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],2.);

  double vol_incr[12] = {0.0}; 

  double edgeSurf_incr[12] = {0.0}; 
  double boundSurf_incr[12] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = (-0.5412658773652741*coeff[0]*fSkin[1])-0.5412658773652741*coeff[0]*fEdge[1]-0.5625*coeff[0]*fSkin[0]+0.5625*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = (-1.4375*coeff[0]*fSkin[1])-0.4375*coeff[0]*fEdge[1]-1.407291281149712*coeff[0]*fSkin[0]+0.5412658773652739*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = (-0.5412658773652741*coeff[0]*fSkin[4])-0.5412658773652741*coeff[0]*fEdge[4]-0.5625*coeff[0]*fSkin[2]+0.5625*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = (-0.5412658773652741*coeff[0]*fSkin[5])-0.5412658773652741*coeff[0]*fEdge[5]-0.5625*coeff[0]*fSkin[3]+0.5625*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = (-1.4375*coeff[0]*fSkin[4])-0.4375*coeff[0]*fEdge[4]-1.407291281149712*coeff[0]*fSkin[2]+0.5412658773652739*coeff[0]*fEdge[2]; 
  edgeSurf_incr[5] = (-1.4375*coeff[0]*fSkin[5])-0.4375*coeff[0]*fEdge[5]-1.407291281149712*coeff[0]*fSkin[3]+0.5412658773652739*coeff[0]*fEdge[3]; 
  edgeSurf_incr[6] = (-0.5412658773652741*coeff[0]*fSkin[7])-0.5412658773652741*coeff[0]*fEdge[7]-0.5625*coeff[0]*fSkin[6]+0.5625*coeff[0]*fEdge[6]; 
  edgeSurf_incr[7] = (-1.4375*coeff[0]*fSkin[7])-0.4375*coeff[0]*fEdge[7]-1.407291281149712*coeff[0]*fSkin[6]+0.5412658773652739*coeff[0]*fEdge[6]; 
  edgeSurf_incr[8] = (-0.5412658773652742*coeff[0]*fSkin[9])-0.5412658773652742*coeff[0]*fEdge[9]-0.5625*coeff[0]*fSkin[8]+0.5625*coeff[0]*fEdge[8]; 
  edgeSurf_incr[9] = (-1.4375*coeff[0]*fSkin[9])-0.4375*coeff[0]*fEdge[9]-1.407291281149713*coeff[0]*fSkin[8]+0.5412658773652742*coeff[0]*fEdge[8]; 
  edgeSurf_incr[10] = (-0.5412658773652742*coeff[0]*fSkin[11])-0.5412658773652742*coeff[0]*fEdge[11]-0.5625*coeff[0]*fSkin[10]+0.5625*coeff[0]*fEdge[10]; 
  edgeSurf_incr[11] = (-1.4375*coeff[0]*fSkin[11])-0.4375*coeff[0]*fEdge[11]-1.407291281149713*coeff[0]*fSkin[10]+0.5412658773652742*coeff[0]*fEdge[10]; 

  boundSurf_incr[1] = 0.8660254037844386*coeff[0]*fSkin[0]-1.0*coeff[0]*fSkin[1]; 
  boundSurf_incr[4] = 0.8660254037844386*coeff[0]*fSkin[2]-1.0*coeff[0]*fSkin[4]; 
  boundSurf_incr[5] = 0.8660254037844386*coeff[0]*fSkin[3]-1.0*coeff[0]*fSkin[5]; 
  boundSurf_incr[7] = 0.8660254037844386*coeff[0]*fSkin[6]-1.0*coeff[0]*fSkin[7]; 
  boundSurf_incr[9] = 0.8660254037844387*coeff[0]*fSkin[8]-1.0*coeff[0]*fSkin[9]; 
  boundSurf_incr[11] = 0.8660254037844387*coeff[0]*fSkin[10]-1.0*coeff[0]*fSkin[11]; 

  } else { 

  edgeSurf_incr[0] = 0.5412658773652741*coeff[0]*fSkin[1]+0.5412658773652741*coeff[0]*fEdge[1]-0.5625*coeff[0]*fSkin[0]+0.5625*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = (-1.4375*coeff[0]*fSkin[1])-0.4375*coeff[0]*fEdge[1]+1.407291281149712*coeff[0]*fSkin[0]-0.5412658773652739*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 0.5412658773652741*coeff[0]*fSkin[4]+0.5412658773652741*coeff[0]*fEdge[4]-0.5625*coeff[0]*fSkin[2]+0.5625*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 0.5412658773652741*coeff[0]*fSkin[5]+0.5412658773652741*coeff[0]*fEdge[5]-0.5625*coeff[0]*fSkin[3]+0.5625*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = (-1.4375*coeff[0]*fSkin[4])-0.4375*coeff[0]*fEdge[4]+1.407291281149712*coeff[0]*fSkin[2]-0.5412658773652739*coeff[0]*fEdge[2]; 
  edgeSurf_incr[5] = (-1.4375*coeff[0]*fSkin[5])-0.4375*coeff[0]*fEdge[5]+1.407291281149712*coeff[0]*fSkin[3]-0.5412658773652739*coeff[0]*fEdge[3]; 
  edgeSurf_incr[6] = 0.5412658773652741*coeff[0]*fSkin[7]+0.5412658773652741*coeff[0]*fEdge[7]-0.5625*coeff[0]*fSkin[6]+0.5625*coeff[0]*fEdge[6]; 
  edgeSurf_incr[7] = (-1.4375*coeff[0]*fSkin[7])-0.4375*coeff[0]*fEdge[7]+1.407291281149712*coeff[0]*fSkin[6]-0.5412658773652739*coeff[0]*fEdge[6]; 
  edgeSurf_incr[8] = 0.5412658773652742*coeff[0]*fSkin[9]+0.5412658773652742*coeff[0]*fEdge[9]-0.5625*coeff[0]*fSkin[8]+0.5625*coeff[0]*fEdge[8]; 
  edgeSurf_incr[9] = (-1.4375*coeff[0]*fSkin[9])-0.4375*coeff[0]*fEdge[9]+1.407291281149713*coeff[0]*fSkin[8]-0.5412658773652742*coeff[0]*fEdge[8]; 
  edgeSurf_incr[10] = 0.5412658773652742*coeff[0]*fSkin[11]+0.5412658773652742*coeff[0]*fEdge[11]-0.5625*coeff[0]*fSkin[10]+0.5625*coeff[0]*fEdge[10]; 
  edgeSurf_incr[11] = (-1.4375*coeff[0]*fSkin[11])-0.4375*coeff[0]*fEdge[11]+1.407291281149713*coeff[0]*fSkin[10]-0.5412658773652742*coeff[0]*fEdge[10]; 

  boundSurf_incr[1] = (-1.0*coeff[0]*fSkin[1])-0.8660254037844386*coeff[0]*fSkin[0]; 
  boundSurf_incr[4] = (-1.0*coeff[0]*fSkin[4])-0.8660254037844386*coeff[0]*fSkin[2]; 
  boundSurf_incr[5] = (-1.0*coeff[0]*fSkin[5])-0.8660254037844386*coeff[0]*fSkin[3]; 
  boundSurf_incr[7] = (-1.0*coeff[0]*fSkin[7])-0.8660254037844386*coeff[0]*fSkin[6]; 
  boundSurf_incr[9] = (-1.0*coeff[0]*fSkin[9])-0.8660254037844387*coeff[0]*fSkin[8]; 
  boundSurf_incr[11] = (-1.0*coeff[0]*fSkin[11])-0.8660254037844387*coeff[0]*fSkin[10]; 

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*Jfac; 
  out[4] += (vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*Jfac; 
  out[5] += (vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*Jfac; 
  out[6] += (vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*Jfac; 
  out[7] += (vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*Jfac; 
  out[8] += (vol_incr[8]+edgeSurf_incr[8]+boundSurf_incr[8])*Jfac; 
  out[9] += (vol_incr[9]+edgeSurf_incr[9]+boundSurf_incr[9])*Jfac; 
  out[10] += (vol_incr[10]+edgeSurf_incr[10]+boundSurf_incr[10])*Jfac; 
  out[11] += (vol_incr[11]+edgeSurf_incr[11]+boundSurf_incr[11])*Jfac; 

  }

  return 0.;
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_boundary_surfx_1x2v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],2.);

  double vol_incr[12] = {0.0}; 
  vol_incr[1] = 2.121320343559642*fSkin[0]*coeff[1]; 
  vol_incr[4] = 2.121320343559642*coeff[1]*fSkin[2]; 
  vol_incr[5] = 2.121320343559642*coeff[1]*fSkin[3]; 
  vol_incr[7] = 2.121320343559642*coeff[1]*fSkin[6]; 
  vol_incr[9] = 2.121320343559642*coeff[1]*fSkin[8]; 
  vol_incr[11] = 2.121320343559642*coeff[1]*fSkin[10]; 

  double edgeSurf_incr[12] = {0.0}; 
  double boundSurf_incr[12] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = (-0.3827327723098713*coeff[0]*fSkin[1])-0.3827327723098713*coeff[0]*fEdge[1]-0.3977475644174328*coeff[0]*fSkin[0]+0.3977475644174328*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = (-0.6123724356957944*coeff[1]*fSkin[1])-1.016465997955662*coeff[0]*fSkin[1]+0.6123724356957944*coeff[1]*fEdge[1]-0.3093592167691142*coeff[0]*fEdge[1]-0.5303300858899105*fSkin[0]*coeff[1]-0.5303300858899105*fEdge[0]*coeff[1]-0.9951052080056654*coeff[0]*fSkin[0]+0.3827327723098711*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = (-0.3827327723098713*coeff[0]*fSkin[4])-0.3827327723098713*coeff[0]*fEdge[4]-0.3977475644174328*coeff[0]*fSkin[2]+0.3977475644174328*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = (-0.3827327723098713*coeff[0]*fSkin[5])-0.3827327723098713*coeff[0]*fEdge[5]-0.3977475644174328*coeff[0]*fSkin[3]+0.3977475644174328*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = (-0.6123724356957944*coeff[1]*fSkin[4])-1.016465997955662*coeff[0]*fSkin[4]+0.6123724356957944*coeff[1]*fEdge[4]-0.3093592167691142*coeff[0]*fEdge[4]-0.5303300858899105*coeff[1]*fSkin[2]-0.9951052080056654*coeff[0]*fSkin[2]-0.5303300858899105*coeff[1]*fEdge[2]+0.3827327723098711*coeff[0]*fEdge[2]; 
  edgeSurf_incr[5] = (-0.6123724356957944*coeff[1]*fSkin[5])-1.016465997955662*coeff[0]*fSkin[5]+0.6123724356957944*coeff[1]*fEdge[5]-0.3093592167691142*coeff[0]*fEdge[5]-0.5303300858899105*coeff[1]*fSkin[3]-0.9951052080056654*coeff[0]*fSkin[3]-0.5303300858899105*coeff[1]*fEdge[3]+0.3827327723098711*coeff[0]*fEdge[3]; 
  edgeSurf_incr[6] = (-0.3827327723098713*coeff[0]*fSkin[7])-0.3827327723098713*coeff[0]*fEdge[7]-0.3977475644174328*coeff[0]*fSkin[6]+0.3977475644174328*coeff[0]*fEdge[6]; 
  edgeSurf_incr[7] = (-0.6123724356957944*coeff[1]*fSkin[7])-1.016465997955662*coeff[0]*fSkin[7]+0.6123724356957944*coeff[1]*fEdge[7]-0.3093592167691142*coeff[0]*fEdge[7]-0.5303300858899105*coeff[1]*fSkin[6]-0.9951052080056654*coeff[0]*fSkin[6]-0.5303300858899105*coeff[1]*fEdge[6]+0.3827327723098711*coeff[0]*fEdge[6]; 
  edgeSurf_incr[8] = (-0.3827327723098714*coeff[0]*fSkin[9])-0.3827327723098714*coeff[0]*fEdge[9]-0.3977475644174328*coeff[0]*fSkin[8]+0.3977475644174328*coeff[0]*fEdge[8]; 
  edgeSurf_incr[9] = (-0.6123724356957944*coeff[1]*fSkin[9])-1.016465997955662*coeff[0]*fSkin[9]+0.6123724356957944*coeff[1]*fEdge[9]-0.3093592167691142*coeff[0]*fEdge[9]-0.5303300858899104*coeff[1]*fSkin[8]-0.9951052080056656*coeff[0]*fSkin[8]-0.5303300858899104*coeff[1]*fEdge[8]+0.3827327723098713*coeff[0]*fEdge[8]; 
  edgeSurf_incr[10] = (-0.3827327723098714*coeff[0]*fSkin[11])-0.3827327723098714*coeff[0]*fEdge[11]-0.3977475644174328*coeff[0]*fSkin[10]+0.3977475644174328*coeff[0]*fEdge[10]; 
  edgeSurf_incr[11] = (-0.6123724356957944*coeff[1]*fSkin[11])-1.016465997955662*coeff[0]*fSkin[11]+0.6123724356957944*coeff[1]*fEdge[11]-0.3093592167691142*coeff[0]*fEdge[11]-0.5303300858899104*coeff[1]*fSkin[10]-0.9951052080056656*coeff[0]*fSkin[10]-0.5303300858899104*coeff[1]*fEdge[10]+0.3827327723098713*coeff[0]*fEdge[10]; 

  boundSurf_incr[1] = 1.224744871391589*coeff[1]*fSkin[1]-0.7071067811865475*coeff[0]*fSkin[1]-1.060660171779821*fSkin[0]*coeff[1]+0.6123724356957944*coeff[0]*fSkin[0]; 
  boundSurf_incr[4] = 1.224744871391589*coeff[1]*fSkin[4]-0.7071067811865475*coeff[0]*fSkin[4]-1.060660171779821*coeff[1]*fSkin[2]+0.6123724356957944*coeff[0]*fSkin[2]; 
  boundSurf_incr[5] = 1.224744871391589*coeff[1]*fSkin[5]-0.7071067811865475*coeff[0]*fSkin[5]-1.060660171779821*coeff[1]*fSkin[3]+0.6123724356957944*coeff[0]*fSkin[3]; 
  boundSurf_incr[7] = 1.224744871391589*coeff[1]*fSkin[7]-0.7071067811865475*coeff[0]*fSkin[7]-1.060660171779821*coeff[1]*fSkin[6]+0.6123724356957944*coeff[0]*fSkin[6]; 
  boundSurf_incr[9] = 1.224744871391589*coeff[1]*fSkin[9]-0.7071067811865475*coeff[0]*fSkin[9]-1.060660171779821*coeff[1]*fSkin[8]+0.6123724356957944*coeff[0]*fSkin[8]; 
  boundSurf_incr[11] = 1.224744871391589*coeff[1]*fSkin[11]-0.7071067811865475*coeff[0]*fSkin[11]-1.060660171779821*coeff[1]*fSkin[10]+0.6123724356957944*coeff[0]*fSkin[10]; 

  } else { 

  edgeSurf_incr[0] = 0.3827327723098713*coeff[0]*fSkin[1]+0.3827327723098713*coeff[0]*fEdge[1]-0.3977475644174328*coeff[0]*fSkin[0]+0.3977475644174328*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 0.6123724356957944*coeff[1]*fSkin[1]-1.016465997955662*coeff[0]*fSkin[1]-0.6123724356957944*coeff[1]*fEdge[1]-0.3093592167691142*coeff[0]*fEdge[1]-0.5303300858899105*fSkin[0]*coeff[1]-0.5303300858899105*fEdge[0]*coeff[1]+0.9951052080056654*coeff[0]*fSkin[0]-0.3827327723098711*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 0.3827327723098713*coeff[0]*fSkin[4]+0.3827327723098713*coeff[0]*fEdge[4]-0.3977475644174328*coeff[0]*fSkin[2]+0.3977475644174328*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 0.3827327723098713*coeff[0]*fSkin[5]+0.3827327723098713*coeff[0]*fEdge[5]-0.3977475644174328*coeff[0]*fSkin[3]+0.3977475644174328*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = 0.6123724356957944*coeff[1]*fSkin[4]-1.016465997955662*coeff[0]*fSkin[4]-0.6123724356957944*coeff[1]*fEdge[4]-0.3093592167691142*coeff[0]*fEdge[4]-0.5303300858899105*coeff[1]*fSkin[2]+0.9951052080056654*coeff[0]*fSkin[2]-0.5303300858899105*coeff[1]*fEdge[2]-0.3827327723098711*coeff[0]*fEdge[2]; 
  edgeSurf_incr[5] = 0.6123724356957944*coeff[1]*fSkin[5]-1.016465997955662*coeff[0]*fSkin[5]-0.6123724356957944*coeff[1]*fEdge[5]-0.3093592167691142*coeff[0]*fEdge[5]-0.5303300858899105*coeff[1]*fSkin[3]+0.9951052080056654*coeff[0]*fSkin[3]-0.5303300858899105*coeff[1]*fEdge[3]-0.3827327723098711*coeff[0]*fEdge[3]; 
  edgeSurf_incr[6] = 0.3827327723098713*coeff[0]*fSkin[7]+0.3827327723098713*coeff[0]*fEdge[7]-0.3977475644174328*coeff[0]*fSkin[6]+0.3977475644174328*coeff[0]*fEdge[6]; 
  edgeSurf_incr[7] = 0.6123724356957944*coeff[1]*fSkin[7]-1.016465997955662*coeff[0]*fSkin[7]-0.6123724356957944*coeff[1]*fEdge[7]-0.3093592167691142*coeff[0]*fEdge[7]-0.5303300858899105*coeff[1]*fSkin[6]+0.9951052080056654*coeff[0]*fSkin[6]-0.5303300858899105*coeff[1]*fEdge[6]-0.3827327723098711*coeff[0]*fEdge[6]; 
  edgeSurf_incr[8] = 0.3827327723098714*coeff[0]*fSkin[9]+0.3827327723098714*coeff[0]*fEdge[9]-0.3977475644174328*coeff[0]*fSkin[8]+0.3977475644174328*coeff[0]*fEdge[8]; 
  edgeSurf_incr[9] = 0.6123724356957944*coeff[1]*fSkin[9]-1.016465997955662*coeff[0]*fSkin[9]-0.6123724356957944*coeff[1]*fEdge[9]-0.3093592167691142*coeff[0]*fEdge[9]-0.5303300858899104*coeff[1]*fSkin[8]+0.9951052080056656*coeff[0]*fSkin[8]-0.5303300858899104*coeff[1]*fEdge[8]-0.3827327723098713*coeff[0]*fEdge[8]; 
  edgeSurf_incr[10] = 0.3827327723098714*coeff[0]*fSkin[11]+0.3827327723098714*coeff[0]*fEdge[11]-0.3977475644174328*coeff[0]*fSkin[10]+0.3977475644174328*coeff[0]*fEdge[10]; 
  edgeSurf_incr[11] = 0.6123724356957944*coeff[1]*fSkin[11]-1.016465997955662*coeff[0]*fSkin[11]-0.6123724356957944*coeff[1]*fEdge[11]-0.3093592167691142*coeff[0]*fEdge[11]-0.5303300858899104*coeff[1]*fSkin[10]+0.9951052080056656*coeff[0]*fSkin[10]-0.5303300858899104*coeff[1]*fEdge[10]-0.3827327723098713*coeff[0]*fEdge[10]; 

  boundSurf_incr[1] = (-1.224744871391589*coeff[1]*fSkin[1])-0.7071067811865475*coeff[0]*fSkin[1]-1.060660171779821*fSkin[0]*coeff[1]-0.6123724356957944*coeff[0]*fSkin[0]; 
  boundSurf_incr[4] = (-1.224744871391589*coeff[1]*fSkin[4])-0.7071067811865475*coeff[0]*fSkin[4]-1.060660171779821*coeff[1]*fSkin[2]-0.6123724356957944*coeff[0]*fSkin[2]; 
  boundSurf_incr[5] = (-1.224744871391589*coeff[1]*fSkin[5])-0.7071067811865475*coeff[0]*fSkin[5]-1.060660171779821*coeff[1]*fSkin[3]-0.6123724356957944*coeff[0]*fSkin[3]; 
  boundSurf_incr[7] = (-1.224744871391589*coeff[1]*fSkin[7])-0.7071067811865475*coeff[0]*fSkin[7]-1.060660171779821*coeff[1]*fSkin[6]-0.6123724356957944*coeff[0]*fSkin[6]; 
  boundSurf_incr[9] = (-1.224744871391589*coeff[1]*fSkin[9])-0.7071067811865475*coeff[0]*fSkin[9]-1.060660171779821*coeff[1]*fSkin[8]-0.6123724356957944*coeff[0]*fSkin[8]; 
  boundSurf_incr[11] = (-1.224744871391589*coeff[1]*fSkin[11])-0.7071067811865475*coeff[0]*fSkin[11]-1.060660171779821*coeff[1]*fSkin[10]-0.6123724356957944*coeff[0]*fSkin[10]; 

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*Jfac; 
  out[4] += (vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*Jfac; 
  out[5] += (vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*Jfac; 
  out[6] += (vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*Jfac; 
  out[7] += (vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*Jfac; 
  out[8] += (vol_incr[8]+edgeSurf_incr[8]+boundSurf_incr[8])*Jfac; 
  out[9] += (vol_incr[9]+edgeSurf_incr[9]+boundSurf_incr[9])*Jfac; 
  out[10] += (vol_incr[10]+edgeSurf_incr[10]+boundSurf_incr[10])*Jfac; 
  out[11] += (vol_incr[11]+edgeSurf_incr[11]+boundSurf_incr[11])*Jfac; 

  }

  return 0.;
}

