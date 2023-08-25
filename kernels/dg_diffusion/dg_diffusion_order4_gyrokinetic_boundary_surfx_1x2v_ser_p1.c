#include <gkyl_dg_diffusion_kernels.h>

GKYL_CU_DH double dg_diffusion_order4_gyrokinetic_boundary_surfx_1x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],4.);

  double vol_incr[12] = {0.0}; 

  double edgeSurf_incr[12] = {0.0}; 
  double boundSurf_incr[12] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = 1.623797632095822*coeff[0]*fSkin[1]+1.623797632095822*coeff[0]*fEdge[1]+0.9375*coeff[0]*fSkin[0]-0.9375*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 3.5625*coeff[0]*fSkin[1]+2.0625*coeff[0]*fEdge[1]+1.623797632095822*coeff[0]*fSkin[0]-1.623797632095822*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 1.623797632095822*coeff[0]*fSkin[4]+1.623797632095822*coeff[0]*fEdge[4]+0.9375*coeff[0]*fSkin[2]-0.9375*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 1.623797632095822*coeff[0]*fSkin[5]+1.623797632095822*coeff[0]*fEdge[5]+0.9375*coeff[0]*fSkin[3]-0.9375*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = 3.5625*coeff[0]*fSkin[4]+2.0625*coeff[0]*fEdge[4]+1.623797632095822*coeff[0]*fSkin[2]-1.623797632095822*coeff[0]*fEdge[2]; 
  edgeSurf_incr[5] = 3.5625*coeff[0]*fSkin[5]+2.0625*coeff[0]*fEdge[5]+1.623797632095822*coeff[0]*fSkin[3]-1.623797632095822*coeff[0]*fEdge[3]; 
  edgeSurf_incr[6] = 1.623797632095822*coeff[0]*fSkin[7]+1.623797632095822*coeff[0]*fEdge[7]+0.9375*coeff[0]*fSkin[6]-0.9375*coeff[0]*fEdge[6]; 
  edgeSurf_incr[7] = 3.5625*coeff[0]*fSkin[7]+2.0625*coeff[0]*fEdge[7]+1.623797632095822*coeff[0]*fSkin[6]-1.623797632095822*coeff[0]*fEdge[6]; 
  edgeSurf_incr[8] = 1.623797632095823*coeff[0]*fSkin[9]+1.623797632095823*coeff[0]*fEdge[9]+0.9375*coeff[0]*fSkin[8]-0.9375*coeff[0]*fEdge[8]; 
  edgeSurf_incr[9] = 3.5625*coeff[0]*fSkin[9]+2.0625*coeff[0]*fEdge[9]+1.623797632095823*coeff[0]*fSkin[8]-1.623797632095823*coeff[0]*fEdge[8]; 
  edgeSurf_incr[10] = 1.623797632095823*coeff[0]*fSkin[11]+1.623797632095823*coeff[0]*fEdge[11]+0.9375*coeff[0]*fSkin[10]-0.9375*coeff[0]*fEdge[10]; 
  edgeSurf_incr[11] = 3.5625*coeff[0]*fSkin[11]+2.0625*coeff[0]*fEdge[11]+1.623797632095823*coeff[0]*fSkin[10]-1.623797632095823*coeff[0]*fEdge[10]; 

  boundSurf_incr[1] = 1.5*coeff[0]*fSkin[1]; 
  boundSurf_incr[4] = 1.5*coeff[0]*fSkin[4]; 
  boundSurf_incr[5] = 1.5*coeff[0]*fSkin[5]; 
  boundSurf_incr[7] = 1.5*coeff[0]*fSkin[7]; 
  boundSurf_incr[9] = 1.5*coeff[0]*fSkin[9]; 
  boundSurf_incr[11] = 1.5*coeff[0]*fSkin[11]; 

  } else { 

  edgeSurf_incr[0] = (-1.623797632095822*coeff[0]*fSkin[1])-1.623797632095822*coeff[0]*fEdge[1]+0.9375*coeff[0]*fSkin[0]-0.9375*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 3.5625*coeff[0]*fSkin[1]+2.0625*coeff[0]*fEdge[1]-1.623797632095822*coeff[0]*fSkin[0]+1.623797632095822*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = (-1.623797632095822*coeff[0]*fSkin[4])-1.623797632095822*coeff[0]*fEdge[4]+0.9375*coeff[0]*fSkin[2]-0.9375*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = (-1.623797632095822*coeff[0]*fSkin[5])-1.623797632095822*coeff[0]*fEdge[5]+0.9375*coeff[0]*fSkin[3]-0.9375*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = 3.5625*coeff[0]*fSkin[4]+2.0625*coeff[0]*fEdge[4]-1.623797632095822*coeff[0]*fSkin[2]+1.623797632095822*coeff[0]*fEdge[2]; 
  edgeSurf_incr[5] = 3.5625*coeff[0]*fSkin[5]+2.0625*coeff[0]*fEdge[5]-1.623797632095822*coeff[0]*fSkin[3]+1.623797632095822*coeff[0]*fEdge[3]; 
  edgeSurf_incr[6] = (-1.623797632095822*coeff[0]*fSkin[7])-1.623797632095822*coeff[0]*fEdge[7]+0.9375*coeff[0]*fSkin[6]-0.9375*coeff[0]*fEdge[6]; 
  edgeSurf_incr[7] = 3.5625*coeff[0]*fSkin[7]+2.0625*coeff[0]*fEdge[7]-1.623797632095822*coeff[0]*fSkin[6]+1.623797632095822*coeff[0]*fEdge[6]; 
  edgeSurf_incr[8] = (-1.623797632095823*coeff[0]*fSkin[9])-1.623797632095823*coeff[0]*fEdge[9]+0.9375*coeff[0]*fSkin[8]-0.9375*coeff[0]*fEdge[8]; 
  edgeSurf_incr[9] = 3.5625*coeff[0]*fSkin[9]+2.0625*coeff[0]*fEdge[9]-1.623797632095823*coeff[0]*fSkin[8]+1.623797632095823*coeff[0]*fEdge[8]; 
  edgeSurf_incr[10] = (-1.623797632095823*coeff[0]*fSkin[11])-1.623797632095823*coeff[0]*fEdge[11]+0.9375*coeff[0]*fSkin[10]-0.9375*coeff[0]*fEdge[10]; 
  edgeSurf_incr[11] = 3.5625*coeff[0]*fSkin[11]+2.0625*coeff[0]*fEdge[11]-1.623797632095823*coeff[0]*fSkin[10]+1.623797632095823*coeff[0]*fEdge[10]; 

  boundSurf_incr[1] = 1.5*coeff[0]*fSkin[1]; 
  boundSurf_incr[4] = 1.5*coeff[0]*fSkin[4]; 
  boundSurf_incr[5] = 1.5*coeff[0]*fSkin[5]; 
  boundSurf_incr[7] = 1.5*coeff[0]*fSkin[7]; 
  boundSurf_incr[9] = 1.5*coeff[0]*fSkin[9]; 
  boundSurf_incr[11] = 1.5*coeff[0]*fSkin[11]; 

  out[0] += -1.0*(vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += -1.0*(vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 
  out[2] += -1.0*(vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac; 
  out[3] += -1.0*(vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*Jfac; 
  out[4] += -1.0*(vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*Jfac; 
  out[5] += -1.0*(vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*Jfac; 
  out[6] += -1.0*(vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*Jfac; 
  out[7] += -1.0*(vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*Jfac; 
  out[8] += -1.0*(vol_incr[8]+edgeSurf_incr[8]+boundSurf_incr[8])*Jfac; 
  out[9] += -1.0*(vol_incr[9]+edgeSurf_incr[9]+boundSurf_incr[9])*Jfac; 
  out[10] += -1.0*(vol_incr[10]+edgeSurf_incr[10]+boundSurf_incr[10])*Jfac; 
  out[11] += -1.0*(vol_incr[11]+edgeSurf_incr[11]+boundSurf_incr[11])*Jfac; 

  }

  return 0.;
}

GKYL_CU_DH double dg_diffusion_order4_gyrokinetic_boundary_surfx_1x2v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],4.);

  double vol_incr[12] = {0.0}; 

  double edgeSurf_incr[12] = {0.0}; 
  double boundSurf_incr[12] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = (-0.5303300858899105*coeff[1]*fSkin[1])+1.148198316929614*coeff[0]*fSkin[1]+0.5303300858899105*coeff[1]*fEdge[1]+1.148198316929614*coeff[0]*fEdge[1]+0.6629126073623879*coeff[0]*fSkin[0]-0.6629126073623879*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = (-0.9185586535436913*coeff[1]*fSkin[1])+2.519067907977075*coeff[0]*fSkin[1]+0.9185586535436913*coeff[1]*fEdge[1]+1.458407736197253*coeff[0]*fEdge[1]+1.148198316929614*coeff[0]*fSkin[0]-1.148198316929614*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = (-0.5303300858899105*coeff[1]*fSkin[4])+1.148198316929614*coeff[0]*fSkin[4]+0.5303300858899105*coeff[1]*fEdge[4]+1.148198316929614*coeff[0]*fEdge[4]+0.6629126073623879*coeff[0]*fSkin[2]-0.6629126073623879*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = (-0.5303300858899105*coeff[1]*fSkin[5])+1.148198316929614*coeff[0]*fSkin[5]+0.5303300858899105*coeff[1]*fEdge[5]+1.148198316929614*coeff[0]*fEdge[5]+0.6629126073623879*coeff[0]*fSkin[3]-0.6629126073623879*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = (-0.9185586535436913*coeff[1]*fSkin[4])+2.519067907977075*coeff[0]*fSkin[4]+0.9185586535436913*coeff[1]*fEdge[4]+1.458407736197253*coeff[0]*fEdge[4]+1.148198316929614*coeff[0]*fSkin[2]-1.148198316929614*coeff[0]*fEdge[2]; 
  edgeSurf_incr[5] = (-0.9185586535436913*coeff[1]*fSkin[5])+2.519067907977075*coeff[0]*fSkin[5]+0.9185586535436913*coeff[1]*fEdge[5]+1.458407736197253*coeff[0]*fEdge[5]+1.148198316929614*coeff[0]*fSkin[3]-1.148198316929614*coeff[0]*fEdge[3]; 
  edgeSurf_incr[6] = (-0.5303300858899105*coeff[1]*fSkin[7])+1.148198316929614*coeff[0]*fSkin[7]+0.5303300858899105*coeff[1]*fEdge[7]+1.148198316929614*coeff[0]*fEdge[7]+0.6629126073623879*coeff[0]*fSkin[6]-0.6629126073623879*coeff[0]*fEdge[6]; 
  edgeSurf_incr[7] = (-0.9185586535436913*coeff[1]*fSkin[7])+2.519067907977075*coeff[0]*fSkin[7]+0.9185586535436913*coeff[1]*fEdge[7]+1.458407736197253*coeff[0]*fEdge[7]+1.148198316929614*coeff[0]*fSkin[6]-1.148198316929614*coeff[0]*fEdge[6]; 
  edgeSurf_incr[8] = (-0.5303300858899104*coeff[1]*fSkin[9])+1.148198316929614*coeff[0]*fSkin[9]+0.5303300858899104*coeff[1]*fEdge[9]+1.148198316929614*coeff[0]*fEdge[9]+0.6629126073623879*coeff[0]*fSkin[8]-0.6629126073623879*coeff[0]*fEdge[8]; 
  edgeSurf_incr[9] = (-0.9185586535436913*coeff[1]*fSkin[9])+2.519067907977075*coeff[0]*fSkin[9]+0.9185586535436913*coeff[1]*fEdge[9]+1.458407736197253*coeff[0]*fEdge[9]+1.148198316929614*coeff[0]*fSkin[8]-1.148198316929614*coeff[0]*fEdge[8]; 
  edgeSurf_incr[10] = (-0.5303300858899104*coeff[1]*fSkin[11])+1.148198316929614*coeff[0]*fSkin[11]+0.5303300858899104*coeff[1]*fEdge[11]+1.148198316929614*coeff[0]*fEdge[11]+0.6629126073623879*coeff[0]*fSkin[10]-0.6629126073623879*coeff[0]*fEdge[10]; 
  edgeSurf_incr[11] = (-0.9185586535436913*coeff[1]*fSkin[11])+2.519067907977075*coeff[0]*fSkin[11]+0.9185586535436913*coeff[1]*fEdge[11]+1.458407736197253*coeff[0]*fEdge[11]+1.148198316929614*coeff[0]*fSkin[10]-1.148198316929614*coeff[0]*fEdge[10]; 

  boundSurf_incr[0] = -1.060660171779821*coeff[1]*fSkin[1]; 
  boundSurf_incr[1] = 1.837117307087383*coeff[1]*fSkin[1]+1.060660171779821*coeff[0]*fSkin[1]; 
  boundSurf_incr[2] = -1.060660171779821*coeff[1]*fSkin[4]; 
  boundSurf_incr[3] = -1.060660171779821*coeff[1]*fSkin[5]; 
  boundSurf_incr[4] = 1.837117307087383*coeff[1]*fSkin[4]+1.060660171779821*coeff[0]*fSkin[4]; 
  boundSurf_incr[5] = 1.837117307087383*coeff[1]*fSkin[5]+1.060660171779821*coeff[0]*fSkin[5]; 
  boundSurf_incr[6] = -1.060660171779821*coeff[1]*fSkin[7]; 
  boundSurf_incr[7] = 1.837117307087383*coeff[1]*fSkin[7]+1.060660171779821*coeff[0]*fSkin[7]; 
  boundSurf_incr[8] = -1.060660171779821*coeff[1]*fSkin[9]; 
  boundSurf_incr[9] = 1.837117307087383*coeff[1]*fSkin[9]+1.060660171779821*coeff[0]*fSkin[9]; 
  boundSurf_incr[10] = -1.060660171779821*coeff[1]*fSkin[11]; 
  boundSurf_incr[11] = 1.837117307087383*coeff[1]*fSkin[11]+1.060660171779821*coeff[0]*fSkin[11]; 

  } else { 

  edgeSurf_incr[0] = (-0.5303300858899105*coeff[1]*fSkin[1])-1.148198316929614*coeff[0]*fSkin[1]+0.5303300858899105*coeff[1]*fEdge[1]-1.148198316929614*coeff[0]*fEdge[1]+0.6629126073623879*coeff[0]*fSkin[0]-0.6629126073623879*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 0.9185586535436913*coeff[1]*fSkin[1]+2.519067907977075*coeff[0]*fSkin[1]-0.9185586535436913*coeff[1]*fEdge[1]+1.458407736197253*coeff[0]*fEdge[1]-1.148198316929614*coeff[0]*fSkin[0]+1.148198316929614*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = (-0.5303300858899105*coeff[1]*fSkin[4])-1.148198316929614*coeff[0]*fSkin[4]+0.5303300858899105*coeff[1]*fEdge[4]-1.148198316929614*coeff[0]*fEdge[4]+0.6629126073623879*coeff[0]*fSkin[2]-0.6629126073623879*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = (-0.5303300858899105*coeff[1]*fSkin[5])-1.148198316929614*coeff[0]*fSkin[5]+0.5303300858899105*coeff[1]*fEdge[5]-1.148198316929614*coeff[0]*fEdge[5]+0.6629126073623879*coeff[0]*fSkin[3]-0.6629126073623879*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = 0.9185586535436913*coeff[1]*fSkin[4]+2.519067907977075*coeff[0]*fSkin[4]-0.9185586535436913*coeff[1]*fEdge[4]+1.458407736197253*coeff[0]*fEdge[4]-1.148198316929614*coeff[0]*fSkin[2]+1.148198316929614*coeff[0]*fEdge[2]; 
  edgeSurf_incr[5] = 0.9185586535436913*coeff[1]*fSkin[5]+2.519067907977075*coeff[0]*fSkin[5]-0.9185586535436913*coeff[1]*fEdge[5]+1.458407736197253*coeff[0]*fEdge[5]-1.148198316929614*coeff[0]*fSkin[3]+1.148198316929614*coeff[0]*fEdge[3]; 
  edgeSurf_incr[6] = (-0.5303300858899105*coeff[1]*fSkin[7])-1.148198316929614*coeff[0]*fSkin[7]+0.5303300858899105*coeff[1]*fEdge[7]-1.148198316929614*coeff[0]*fEdge[7]+0.6629126073623879*coeff[0]*fSkin[6]-0.6629126073623879*coeff[0]*fEdge[6]; 
  edgeSurf_incr[7] = 0.9185586535436913*coeff[1]*fSkin[7]+2.519067907977075*coeff[0]*fSkin[7]-0.9185586535436913*coeff[1]*fEdge[7]+1.458407736197253*coeff[0]*fEdge[7]-1.148198316929614*coeff[0]*fSkin[6]+1.148198316929614*coeff[0]*fEdge[6]; 
  edgeSurf_incr[8] = (-0.5303300858899104*coeff[1]*fSkin[9])-1.148198316929614*coeff[0]*fSkin[9]+0.5303300858899104*coeff[1]*fEdge[9]-1.148198316929614*coeff[0]*fEdge[9]+0.6629126073623879*coeff[0]*fSkin[8]-0.6629126073623879*coeff[0]*fEdge[8]; 
  edgeSurf_incr[9] = 0.9185586535436913*coeff[1]*fSkin[9]+2.519067907977075*coeff[0]*fSkin[9]-0.9185586535436913*coeff[1]*fEdge[9]+1.458407736197253*coeff[0]*fEdge[9]-1.148198316929614*coeff[0]*fSkin[8]+1.148198316929614*coeff[0]*fEdge[8]; 
  edgeSurf_incr[10] = (-0.5303300858899104*coeff[1]*fSkin[11])-1.148198316929614*coeff[0]*fSkin[11]+0.5303300858899104*coeff[1]*fEdge[11]-1.148198316929614*coeff[0]*fEdge[11]+0.6629126073623879*coeff[0]*fSkin[10]-0.6629126073623879*coeff[0]*fEdge[10]; 
  edgeSurf_incr[11] = 0.9185586535436913*coeff[1]*fSkin[11]+2.519067907977075*coeff[0]*fSkin[11]-0.9185586535436913*coeff[1]*fEdge[11]+1.458407736197253*coeff[0]*fEdge[11]-1.148198316929614*coeff[0]*fSkin[10]+1.148198316929614*coeff[0]*fEdge[10]; 

  boundSurf_incr[0] = -1.060660171779821*coeff[1]*fSkin[1]; 
  boundSurf_incr[1] = 1.060660171779821*coeff[0]*fSkin[1]-1.837117307087383*coeff[1]*fSkin[1]; 
  boundSurf_incr[2] = -1.060660171779821*coeff[1]*fSkin[4]; 
  boundSurf_incr[3] = -1.060660171779821*coeff[1]*fSkin[5]; 
  boundSurf_incr[4] = 1.060660171779821*coeff[0]*fSkin[4]-1.837117307087383*coeff[1]*fSkin[4]; 
  boundSurf_incr[5] = 1.060660171779821*coeff[0]*fSkin[5]-1.837117307087383*coeff[1]*fSkin[5]; 
  boundSurf_incr[6] = -1.060660171779821*coeff[1]*fSkin[7]; 
  boundSurf_incr[7] = 1.060660171779821*coeff[0]*fSkin[7]-1.837117307087383*coeff[1]*fSkin[7]; 
  boundSurf_incr[8] = -1.060660171779821*coeff[1]*fSkin[9]; 
  boundSurf_incr[9] = 1.060660171779821*coeff[0]*fSkin[9]-1.837117307087383*coeff[1]*fSkin[9]; 
  boundSurf_incr[10] = -1.060660171779821*coeff[1]*fSkin[11]; 
  boundSurf_incr[11] = 1.060660171779821*coeff[0]*fSkin[11]-1.837117307087383*coeff[1]*fSkin[11]; 

  out[0] += -1.0*(vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += -1.0*(vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 
  out[2] += -1.0*(vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac; 
  out[3] += -1.0*(vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*Jfac; 
  out[4] += -1.0*(vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*Jfac; 
  out[5] += -1.0*(vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*Jfac; 
  out[6] += -1.0*(vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*Jfac; 
  out[7] += -1.0*(vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*Jfac; 
  out[8] += -1.0*(vol_incr[8]+edgeSurf_incr[8]+boundSurf_incr[8])*Jfac; 
  out[9] += -1.0*(vol_incr[9]+edgeSurf_incr[9]+boundSurf_incr[9])*Jfac; 
  out[10] += -1.0*(vol_incr[10]+edgeSurf_incr[10]+boundSurf_incr[10])*Jfac; 
  out[11] += -1.0*(vol_incr[11]+edgeSurf_incr[11]+boundSurf_incr[11])*Jfac; 

  }

  return 0.;
}

