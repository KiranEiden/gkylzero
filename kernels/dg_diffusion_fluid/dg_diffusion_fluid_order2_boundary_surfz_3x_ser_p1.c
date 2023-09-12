#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfz_3x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[2],2.);

  double vol_incr[8] = {0.0}; 

  double edgeSurf_incr[8] = {0.0}; 
  double boundSurf_incr[8] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = (-0.5412658773652741*coeff[2]*fSkin[3])-0.5412658773652741*coeff[2]*fEdge[3]-0.5625*fSkin[0]*coeff[2]+0.5625*fEdge[0]*coeff[2]; 
  edgeSurf_incr[1] = (-0.5412658773652741*coeff[2]*fSkin[5])-0.5412658773652741*coeff[2]*fEdge[5]-0.5625*fSkin[1]*coeff[2]+0.5625*fEdge[1]*coeff[2]; 
  edgeSurf_incr[2] = (-0.5412658773652741*coeff[2]*fSkin[6])-0.5412658773652741*coeff[2]*fEdge[6]-0.5625*coeff[2]*fSkin[2]+0.5625*coeff[2]*fEdge[2]; 
  edgeSurf_incr[3] = (-1.4375*coeff[2]*fSkin[3])-0.4375*coeff[2]*fEdge[3]-1.407291281149712*fSkin[0]*coeff[2]+0.5412658773652739*fEdge[0]*coeff[2]; 
  edgeSurf_incr[4] = (-0.5412658773652741*coeff[2]*fSkin[7])-0.5412658773652741*coeff[2]*fEdge[7]-0.5625*coeff[2]*fSkin[4]+0.5625*coeff[2]*fEdge[4]; 
  edgeSurf_incr[5] = (-1.4375*coeff[2]*fSkin[5])-0.4375*coeff[2]*fEdge[5]-1.407291281149712*fSkin[1]*coeff[2]+0.5412658773652739*fEdge[1]*coeff[2]; 
  edgeSurf_incr[6] = (-1.4375*coeff[2]*fSkin[6])-0.4375*coeff[2]*fEdge[6]-1.407291281149712*coeff[2]*fSkin[2]+0.5412658773652739*coeff[2]*fEdge[2]; 
  edgeSurf_incr[7] = (-1.4375*coeff[2]*fSkin[7])-0.4375*coeff[2]*fEdge[7]-1.407291281149712*coeff[2]*fSkin[4]+0.5412658773652739*coeff[2]*fEdge[4]; 

  boundSurf_incr[3] = 0.8660254037844386*fSkin[0]*coeff[2]-1.0*coeff[2]*fSkin[3]; 
  boundSurf_incr[5] = 0.8660254037844386*fSkin[1]*coeff[2]-1.0*coeff[2]*fSkin[5]; 
  boundSurf_incr[6] = 0.8660254037844386*coeff[2]*fSkin[2]-1.0*coeff[2]*fSkin[6]; 
  boundSurf_incr[7] = 0.8660254037844386*coeff[2]*fSkin[4]-1.0*coeff[2]*fSkin[7]; 

  } else { 

  edgeSurf_incr[0] = 0.5412658773652741*coeff[2]*fSkin[3]+0.5412658773652741*coeff[2]*fEdge[3]-0.5625*fSkin[0]*coeff[2]+0.5625*fEdge[0]*coeff[2]; 
  edgeSurf_incr[1] = 0.5412658773652741*coeff[2]*fSkin[5]+0.5412658773652741*coeff[2]*fEdge[5]-0.5625*fSkin[1]*coeff[2]+0.5625*fEdge[1]*coeff[2]; 
  edgeSurf_incr[2] = 0.5412658773652741*coeff[2]*fSkin[6]+0.5412658773652741*coeff[2]*fEdge[6]-0.5625*coeff[2]*fSkin[2]+0.5625*coeff[2]*fEdge[2]; 
  edgeSurf_incr[3] = (-1.4375*coeff[2]*fSkin[3])-0.4375*coeff[2]*fEdge[3]+1.407291281149712*fSkin[0]*coeff[2]-0.5412658773652739*fEdge[0]*coeff[2]; 
  edgeSurf_incr[4] = 0.5412658773652741*coeff[2]*fSkin[7]+0.5412658773652741*coeff[2]*fEdge[7]-0.5625*coeff[2]*fSkin[4]+0.5625*coeff[2]*fEdge[4]; 
  edgeSurf_incr[5] = (-1.4375*coeff[2]*fSkin[5])-0.4375*coeff[2]*fEdge[5]+1.407291281149712*fSkin[1]*coeff[2]-0.5412658773652739*fEdge[1]*coeff[2]; 
  edgeSurf_incr[6] = (-1.4375*coeff[2]*fSkin[6])-0.4375*coeff[2]*fEdge[6]+1.407291281149712*coeff[2]*fSkin[2]-0.5412658773652739*coeff[2]*fEdge[2]; 
  edgeSurf_incr[7] = (-1.4375*coeff[2]*fSkin[7])-0.4375*coeff[2]*fEdge[7]+1.407291281149712*coeff[2]*fSkin[4]-0.5412658773652739*coeff[2]*fEdge[4]; 

  boundSurf_incr[3] = (-1.0*coeff[2]*fSkin[3])-0.8660254037844386*fSkin[0]*coeff[2]; 
  boundSurf_incr[5] = (-1.0*coeff[2]*fSkin[5])-0.8660254037844386*fSkin[1]*coeff[2]; 
  boundSurf_incr[6] = (-1.0*coeff[2]*fSkin[6])-0.8660254037844386*coeff[2]*fSkin[2]; 
  boundSurf_incr[7] = (-1.0*coeff[2]*fSkin[7])-0.8660254037844386*coeff[2]*fSkin[4]; 

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*Jfac; 
  out[4] += (vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*Jfac; 
  out[5] += (vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*Jfac; 
  out[6] += (vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*Jfac; 
  out[7] += (vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*Jfac; 

  }

  return 0.;
}

GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfz_3x_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[2],2.);

  double vol_incr[8] = {0.0}; 
  vol_incr[3] = 1.060660171779821*fSkin[4]*coeff[23]+1.060660171779821*fSkin[2]*coeff[22]+1.060660171779821*fSkin[1]*coeff[21]+1.060660171779821*fSkin[0]*coeff[19]; 
  vol_incr[5] = 1.060660171779821*fSkin[2]*coeff[23]+1.060660171779821*fSkin[4]*coeff[22]+1.060660171779821*fSkin[0]*coeff[21]+1.060660171779821*fSkin[1]*coeff[19]; 
  vol_incr[6] = 1.060660171779821*fSkin[1]*coeff[23]+1.060660171779821*fSkin[0]*coeff[22]+1.060660171779821*fSkin[4]*coeff[21]+1.060660171779821*fSkin[2]*coeff[19]; 
  vol_incr[7] = 1.060660171779821*fSkin[0]*coeff[23]+1.060660171779821*fSkin[1]*coeff[22]+1.060660171779821*fSkin[2]*coeff[21]+1.060660171779821*fSkin[4]*coeff[19]; 

  double edgeSurf_incr[8] = {0.0}; 
  double boundSurf_incr[8] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = (-0.1913663861549356*fSkin[7]*coeff[20])-0.1913663861549356*fEdge[7]*coeff[20]-0.1988737822087163*fSkin[4]*coeff[20]+0.1988737822087163*fEdge[4]*coeff[20]-0.1913663861549356*fSkin[6]*coeff[18]-0.1913663861549356*fEdge[6]*coeff[18]-0.1988737822087163*fSkin[2]*coeff[18]+0.1988737822087163*fEdge[2]*coeff[18]-0.1913663861549356*fSkin[5]*coeff[17]-0.1913663861549356*fEdge[5]*coeff[17]-0.1988737822087163*fSkin[1]*coeff[17]+0.1988737822087163*fEdge[1]*coeff[17]-0.1913663861549356*fSkin[3]*coeff[16]-0.1913663861549356*fEdge[3]*coeff[16]-0.1988737822087163*fSkin[0]*coeff[16]+0.1988737822087163*fEdge[0]*coeff[16]; 
  edgeSurf_incr[1] = (-0.1913663861549356*fSkin[6]*coeff[20])-0.1913663861549356*fEdge[6]*coeff[20]-0.1988737822087163*fSkin[2]*coeff[20]+0.1988737822087163*fEdge[2]*coeff[20]-0.1913663861549356*fSkin[7]*coeff[18]-0.1913663861549356*fEdge[7]*coeff[18]-0.1988737822087163*fSkin[4]*coeff[18]+0.1988737822087163*fEdge[4]*coeff[18]-0.1913663861549356*fSkin[3]*coeff[17]-0.1913663861549356*fEdge[3]*coeff[17]-0.1988737822087163*fSkin[0]*coeff[17]+0.1988737822087163*fEdge[0]*coeff[17]-0.1913663861549356*fSkin[5]*coeff[16]-0.1913663861549356*fEdge[5]*coeff[16]-0.1988737822087163*fSkin[1]*coeff[16]+0.1988737822087163*fEdge[1]*coeff[16]; 
  edgeSurf_incr[2] = (-0.1913663861549356*fSkin[5]*coeff[20])-0.1913663861549356*fEdge[5]*coeff[20]-0.1988737822087163*fSkin[1]*coeff[20]+0.1988737822087163*fEdge[1]*coeff[20]-0.1913663861549356*fSkin[3]*coeff[18]-0.1913663861549356*fEdge[3]*coeff[18]-0.1988737822087163*fSkin[0]*coeff[18]+0.1988737822087163*fEdge[0]*coeff[18]-0.1913663861549356*fSkin[7]*coeff[17]-0.1913663861549356*fEdge[7]*coeff[17]-0.1988737822087163*fSkin[4]*coeff[17]+0.1988737822087163*fEdge[4]*coeff[17]-0.1913663861549356*fSkin[6]*coeff[16]-0.1913663861549356*fEdge[6]*coeff[16]-0.1988737822087163*fSkin[2]*coeff[16]+0.1988737822087163*fEdge[2]*coeff[16]; 
  edgeSurf_incr[3] = (-0.3061862178478971*fSkin[7]*coeff[23])+0.3061862178478971*fEdge[7]*coeff[23]-0.2651650429449552*fSkin[4]*coeff[23]-0.2651650429449552*fEdge[4]*coeff[23]-0.3061862178478971*fSkin[6]*coeff[22]+0.3061862178478971*fEdge[6]*coeff[22]-0.2651650429449552*fSkin[2]*coeff[22]-0.2651650429449552*fEdge[2]*coeff[22]-0.3061862178478971*fSkin[5]*coeff[21]+0.3061862178478971*fEdge[5]*coeff[21]-0.2651650429449552*fSkin[1]*coeff[21]-0.2651650429449552*fEdge[1]*coeff[21]-0.5082329989778307*fSkin[7]*coeff[20]-0.1546796083845571*fEdge[7]*coeff[20]-0.4975526040028326*fSkin[4]*coeff[20]+0.1913663861549355*fEdge[4]*coeff[20]-0.3061862178478971*fSkin[3]*coeff[19]+0.3061862178478971*fEdge[3]*coeff[19]-0.2651650429449552*fSkin[0]*coeff[19]-0.2651650429449552*fEdge[0]*coeff[19]-0.5082329989778307*fSkin[6]*coeff[18]-0.1546796083845571*fEdge[6]*coeff[18]-0.4975526040028326*fSkin[2]*coeff[18]+0.1913663861549355*fEdge[2]*coeff[18]-0.5082329989778307*fSkin[5]*coeff[17]-0.1546796083845571*fEdge[5]*coeff[17]-0.4975526040028326*fSkin[1]*coeff[17]+0.1913663861549355*fEdge[1]*coeff[17]-0.5082329989778307*fSkin[3]*coeff[16]-0.1546796083845571*fEdge[3]*coeff[16]-0.4975526040028326*fSkin[0]*coeff[16]+0.1913663861549355*fEdge[0]*coeff[16]; 
  edgeSurf_incr[4] = (-0.1913663861549356*fSkin[3]*coeff[20])-0.1913663861549356*fEdge[3]*coeff[20]-0.1988737822087163*fSkin[0]*coeff[20]+0.1988737822087163*fEdge[0]*coeff[20]-0.1913663861549356*fSkin[5]*coeff[18]-0.1913663861549356*fEdge[5]*coeff[18]-0.1988737822087163*fSkin[1]*coeff[18]+0.1988737822087163*fEdge[1]*coeff[18]-0.1913663861549356*fSkin[6]*coeff[17]-0.1913663861549356*fEdge[6]*coeff[17]-0.1988737822087163*fSkin[2]*coeff[17]+0.1988737822087163*fEdge[2]*coeff[17]-0.1913663861549356*fSkin[7]*coeff[16]-0.1913663861549356*fEdge[7]*coeff[16]-0.1988737822087163*fSkin[4]*coeff[16]+0.1988737822087163*fEdge[4]*coeff[16]; 
  edgeSurf_incr[5] = (-0.3061862178478971*fSkin[6]*coeff[23])+0.3061862178478971*fEdge[6]*coeff[23]-0.2651650429449552*fSkin[2]*coeff[23]-0.2651650429449552*fEdge[2]*coeff[23]-0.3061862178478971*fSkin[7]*coeff[22]+0.3061862178478971*fEdge[7]*coeff[22]-0.2651650429449552*fSkin[4]*coeff[22]-0.2651650429449552*fEdge[4]*coeff[22]-0.3061862178478971*fSkin[3]*coeff[21]+0.3061862178478971*fEdge[3]*coeff[21]-0.2651650429449552*fSkin[0]*coeff[21]-0.2651650429449552*fEdge[0]*coeff[21]-0.5082329989778307*fSkin[6]*coeff[20]-0.1546796083845571*fEdge[6]*coeff[20]-0.4975526040028326*fSkin[2]*coeff[20]+0.1913663861549355*fEdge[2]*coeff[20]-0.3061862178478971*fSkin[5]*coeff[19]+0.3061862178478971*fEdge[5]*coeff[19]-0.2651650429449552*fSkin[1]*coeff[19]-0.2651650429449552*fEdge[1]*coeff[19]-0.5082329989778307*fSkin[7]*coeff[18]-0.1546796083845571*fEdge[7]*coeff[18]-0.4975526040028326*fSkin[4]*coeff[18]+0.1913663861549355*fEdge[4]*coeff[18]-0.5082329989778307*fSkin[3]*coeff[17]-0.1546796083845571*fEdge[3]*coeff[17]-0.4975526040028326*fSkin[0]*coeff[17]+0.1913663861549355*fEdge[0]*coeff[17]-0.5082329989778307*fSkin[5]*coeff[16]-0.1546796083845571*fEdge[5]*coeff[16]-0.4975526040028326*fSkin[1]*coeff[16]+0.1913663861549355*fEdge[1]*coeff[16]; 
  edgeSurf_incr[6] = (-0.3061862178478971*fSkin[5]*coeff[23])+0.3061862178478971*fEdge[5]*coeff[23]-0.2651650429449552*fSkin[1]*coeff[23]-0.2651650429449552*fEdge[1]*coeff[23]-0.3061862178478971*fSkin[3]*coeff[22]+0.3061862178478971*fEdge[3]*coeff[22]-0.2651650429449552*fSkin[0]*coeff[22]-0.2651650429449552*fEdge[0]*coeff[22]-0.3061862178478971*fSkin[7]*coeff[21]+0.3061862178478971*fEdge[7]*coeff[21]-0.2651650429449552*fSkin[4]*coeff[21]-0.2651650429449552*fEdge[4]*coeff[21]-0.5082329989778307*fSkin[5]*coeff[20]-0.1546796083845571*fEdge[5]*coeff[20]-0.4975526040028326*fSkin[1]*coeff[20]+0.1913663861549355*fEdge[1]*coeff[20]-0.3061862178478971*fSkin[6]*coeff[19]+0.3061862178478971*fEdge[6]*coeff[19]-0.2651650429449552*fSkin[2]*coeff[19]-0.2651650429449552*fEdge[2]*coeff[19]-0.5082329989778307*fSkin[3]*coeff[18]-0.1546796083845571*fEdge[3]*coeff[18]-0.4975526040028326*fSkin[0]*coeff[18]+0.1913663861549355*fEdge[0]*coeff[18]-0.5082329989778307*fSkin[7]*coeff[17]-0.1546796083845571*fEdge[7]*coeff[17]-0.4975526040028326*fSkin[4]*coeff[17]+0.1913663861549355*fEdge[4]*coeff[17]-0.5082329989778307*fSkin[6]*coeff[16]-0.1546796083845571*fEdge[6]*coeff[16]-0.4975526040028326*fSkin[2]*coeff[16]+0.1913663861549355*fEdge[2]*coeff[16]; 
  edgeSurf_incr[7] = (-0.3061862178478971*fSkin[3]*coeff[23])+0.3061862178478971*fEdge[3]*coeff[23]-0.2651650429449552*fSkin[0]*coeff[23]-0.2651650429449552*fEdge[0]*coeff[23]-0.3061862178478971*fSkin[5]*coeff[22]+0.3061862178478971*fEdge[5]*coeff[22]-0.2651650429449552*fSkin[1]*coeff[22]-0.2651650429449552*fEdge[1]*coeff[22]-0.3061862178478971*fSkin[6]*coeff[21]+0.3061862178478971*fEdge[6]*coeff[21]-0.2651650429449552*fSkin[2]*coeff[21]-0.2651650429449552*fEdge[2]*coeff[21]-0.5082329989778307*fSkin[3]*coeff[20]-0.1546796083845571*fEdge[3]*coeff[20]-0.4975526040028326*fSkin[0]*coeff[20]+0.1913663861549355*fEdge[0]*coeff[20]-0.3061862178478971*fSkin[7]*coeff[19]+0.3061862178478971*fEdge[7]*coeff[19]-0.2651650429449552*fSkin[4]*coeff[19]-0.2651650429449552*fEdge[4]*coeff[19]-0.5082329989778307*fSkin[5]*coeff[18]-0.1546796083845571*fEdge[5]*coeff[18]-0.4975526040028326*fSkin[1]*coeff[18]+0.1913663861549355*fEdge[1]*coeff[18]-0.5082329989778307*fSkin[6]*coeff[17]-0.1546796083845571*fEdge[6]*coeff[17]-0.4975526040028326*fSkin[2]*coeff[17]+0.1913663861549355*fEdge[2]*coeff[17]-0.5082329989778307*fSkin[7]*coeff[16]-0.1546796083845571*fEdge[7]*coeff[16]-0.4975526040028326*fSkin[4]*coeff[16]+0.1913663861549355*fEdge[4]*coeff[16]; 

  boundSurf_incr[3] = 0.6123724356957944*fSkin[7]*coeff[23]-0.5303300858899105*fSkin[4]*coeff[23]+0.6123724356957944*fSkin[6]*coeff[22]-0.5303300858899105*fSkin[2]*coeff[22]+0.6123724356957944*fSkin[5]*coeff[21]-0.5303300858899105*fSkin[1]*coeff[21]-0.3535533905932737*fSkin[7]*coeff[20]+0.3061862178478971*fSkin[4]*coeff[20]+0.6123724356957944*fSkin[3]*coeff[19]-0.5303300858899105*fSkin[0]*coeff[19]-0.3535533905932737*fSkin[6]*coeff[18]+0.3061862178478971*fSkin[2]*coeff[18]-0.3535533905932737*fSkin[5]*coeff[17]+0.3061862178478971*fSkin[1]*coeff[17]-0.3535533905932737*fSkin[3]*coeff[16]+0.3061862178478971*fSkin[0]*coeff[16]; 
  boundSurf_incr[5] = 0.6123724356957944*fSkin[6]*coeff[23]-0.5303300858899105*fSkin[2]*coeff[23]+0.6123724356957944*fSkin[7]*coeff[22]-0.5303300858899105*fSkin[4]*coeff[22]+0.6123724356957944*fSkin[3]*coeff[21]-0.5303300858899105*fSkin[0]*coeff[21]-0.3535533905932737*fSkin[6]*coeff[20]+0.3061862178478971*fSkin[2]*coeff[20]+0.6123724356957944*fSkin[5]*coeff[19]-0.5303300858899105*fSkin[1]*coeff[19]-0.3535533905932737*fSkin[7]*coeff[18]+0.3061862178478971*fSkin[4]*coeff[18]-0.3535533905932737*fSkin[3]*coeff[17]+0.3061862178478971*fSkin[0]*coeff[17]-0.3535533905932737*fSkin[5]*coeff[16]+0.3061862178478971*fSkin[1]*coeff[16]; 
  boundSurf_incr[6] = 0.6123724356957944*fSkin[5]*coeff[23]-0.5303300858899105*fSkin[1]*coeff[23]+0.6123724356957944*fSkin[3]*coeff[22]-0.5303300858899105*fSkin[0]*coeff[22]+0.6123724356957944*fSkin[7]*coeff[21]-0.5303300858899105*fSkin[4]*coeff[21]-0.3535533905932737*fSkin[5]*coeff[20]+0.3061862178478971*fSkin[1]*coeff[20]+0.6123724356957944*fSkin[6]*coeff[19]-0.5303300858899105*fSkin[2]*coeff[19]-0.3535533905932737*fSkin[3]*coeff[18]+0.3061862178478971*fSkin[0]*coeff[18]-0.3535533905932737*fSkin[7]*coeff[17]+0.3061862178478971*fSkin[4]*coeff[17]-0.3535533905932737*fSkin[6]*coeff[16]+0.3061862178478971*fSkin[2]*coeff[16]; 
  boundSurf_incr[7] = 0.6123724356957944*fSkin[3]*coeff[23]-0.5303300858899105*fSkin[0]*coeff[23]+0.6123724356957944*fSkin[5]*coeff[22]-0.5303300858899105*fSkin[1]*coeff[22]+0.6123724356957944*fSkin[6]*coeff[21]-0.5303300858899105*fSkin[2]*coeff[21]-0.3535533905932737*fSkin[3]*coeff[20]+0.3061862178478971*fSkin[0]*coeff[20]+0.6123724356957944*fSkin[7]*coeff[19]-0.5303300858899105*fSkin[4]*coeff[19]-0.3535533905932737*fSkin[5]*coeff[18]+0.3061862178478971*fSkin[1]*coeff[18]-0.3535533905932737*fSkin[6]*coeff[17]+0.3061862178478971*fSkin[2]*coeff[17]-0.3535533905932737*fSkin[7]*coeff[16]+0.3061862178478971*fSkin[4]*coeff[16]; 

  } else { 

  edgeSurf_incr[0] = 0.1913663861549356*fSkin[7]*coeff[20]+0.1913663861549356*fEdge[7]*coeff[20]-0.1988737822087163*fSkin[4]*coeff[20]+0.1988737822087163*fEdge[4]*coeff[20]+0.1913663861549356*fSkin[6]*coeff[18]+0.1913663861549356*fEdge[6]*coeff[18]-0.1988737822087163*fSkin[2]*coeff[18]+0.1988737822087163*fEdge[2]*coeff[18]+0.1913663861549356*fSkin[5]*coeff[17]+0.1913663861549356*fEdge[5]*coeff[17]-0.1988737822087163*fSkin[1]*coeff[17]+0.1988737822087163*fEdge[1]*coeff[17]+0.1913663861549356*fSkin[3]*coeff[16]+0.1913663861549356*fEdge[3]*coeff[16]-0.1988737822087163*fSkin[0]*coeff[16]+0.1988737822087163*fEdge[0]*coeff[16]; 
  edgeSurf_incr[1] = 0.1913663861549356*fSkin[6]*coeff[20]+0.1913663861549356*fEdge[6]*coeff[20]-0.1988737822087163*fSkin[2]*coeff[20]+0.1988737822087163*fEdge[2]*coeff[20]+0.1913663861549356*fSkin[7]*coeff[18]+0.1913663861549356*fEdge[7]*coeff[18]-0.1988737822087163*fSkin[4]*coeff[18]+0.1988737822087163*fEdge[4]*coeff[18]+0.1913663861549356*fSkin[3]*coeff[17]+0.1913663861549356*fEdge[3]*coeff[17]-0.1988737822087163*fSkin[0]*coeff[17]+0.1988737822087163*fEdge[0]*coeff[17]+0.1913663861549356*fSkin[5]*coeff[16]+0.1913663861549356*fEdge[5]*coeff[16]-0.1988737822087163*fSkin[1]*coeff[16]+0.1988737822087163*fEdge[1]*coeff[16]; 
  edgeSurf_incr[2] = 0.1913663861549356*fSkin[5]*coeff[20]+0.1913663861549356*fEdge[5]*coeff[20]-0.1988737822087163*fSkin[1]*coeff[20]+0.1988737822087163*fEdge[1]*coeff[20]+0.1913663861549356*fSkin[3]*coeff[18]+0.1913663861549356*fEdge[3]*coeff[18]-0.1988737822087163*fSkin[0]*coeff[18]+0.1988737822087163*fEdge[0]*coeff[18]+0.1913663861549356*fSkin[7]*coeff[17]+0.1913663861549356*fEdge[7]*coeff[17]-0.1988737822087163*fSkin[4]*coeff[17]+0.1988737822087163*fEdge[4]*coeff[17]+0.1913663861549356*fSkin[6]*coeff[16]+0.1913663861549356*fEdge[6]*coeff[16]-0.1988737822087163*fSkin[2]*coeff[16]+0.1988737822087163*fEdge[2]*coeff[16]; 
  edgeSurf_incr[3] = 0.3061862178478971*fSkin[7]*coeff[23]-0.3061862178478971*fEdge[7]*coeff[23]-0.2651650429449552*fSkin[4]*coeff[23]-0.2651650429449552*fEdge[4]*coeff[23]+0.3061862178478971*fSkin[6]*coeff[22]-0.3061862178478971*fEdge[6]*coeff[22]-0.2651650429449552*fSkin[2]*coeff[22]-0.2651650429449552*fEdge[2]*coeff[22]+0.3061862178478971*fSkin[5]*coeff[21]-0.3061862178478971*fEdge[5]*coeff[21]-0.2651650429449552*fSkin[1]*coeff[21]-0.2651650429449552*fEdge[1]*coeff[21]-0.5082329989778307*fSkin[7]*coeff[20]-0.1546796083845571*fEdge[7]*coeff[20]+0.4975526040028326*fSkin[4]*coeff[20]-0.1913663861549355*fEdge[4]*coeff[20]+0.3061862178478971*fSkin[3]*coeff[19]-0.3061862178478971*fEdge[3]*coeff[19]-0.2651650429449552*fSkin[0]*coeff[19]-0.2651650429449552*fEdge[0]*coeff[19]-0.5082329989778307*fSkin[6]*coeff[18]-0.1546796083845571*fEdge[6]*coeff[18]+0.4975526040028326*fSkin[2]*coeff[18]-0.1913663861549355*fEdge[2]*coeff[18]-0.5082329989778307*fSkin[5]*coeff[17]-0.1546796083845571*fEdge[5]*coeff[17]+0.4975526040028326*fSkin[1]*coeff[17]-0.1913663861549355*fEdge[1]*coeff[17]-0.5082329989778307*fSkin[3]*coeff[16]-0.1546796083845571*fEdge[3]*coeff[16]+0.4975526040028326*fSkin[0]*coeff[16]-0.1913663861549355*fEdge[0]*coeff[16]; 
  edgeSurf_incr[4] = 0.1913663861549356*fSkin[3]*coeff[20]+0.1913663861549356*fEdge[3]*coeff[20]-0.1988737822087163*fSkin[0]*coeff[20]+0.1988737822087163*fEdge[0]*coeff[20]+0.1913663861549356*fSkin[5]*coeff[18]+0.1913663861549356*fEdge[5]*coeff[18]-0.1988737822087163*fSkin[1]*coeff[18]+0.1988737822087163*fEdge[1]*coeff[18]+0.1913663861549356*fSkin[6]*coeff[17]+0.1913663861549356*fEdge[6]*coeff[17]-0.1988737822087163*fSkin[2]*coeff[17]+0.1988737822087163*fEdge[2]*coeff[17]+0.1913663861549356*fSkin[7]*coeff[16]+0.1913663861549356*fEdge[7]*coeff[16]-0.1988737822087163*fSkin[4]*coeff[16]+0.1988737822087163*fEdge[4]*coeff[16]; 
  edgeSurf_incr[5] = 0.3061862178478971*fSkin[6]*coeff[23]-0.3061862178478971*fEdge[6]*coeff[23]-0.2651650429449552*fSkin[2]*coeff[23]-0.2651650429449552*fEdge[2]*coeff[23]+0.3061862178478971*fSkin[7]*coeff[22]-0.3061862178478971*fEdge[7]*coeff[22]-0.2651650429449552*fSkin[4]*coeff[22]-0.2651650429449552*fEdge[4]*coeff[22]+0.3061862178478971*fSkin[3]*coeff[21]-0.3061862178478971*fEdge[3]*coeff[21]-0.2651650429449552*fSkin[0]*coeff[21]-0.2651650429449552*fEdge[0]*coeff[21]-0.5082329989778307*fSkin[6]*coeff[20]-0.1546796083845571*fEdge[6]*coeff[20]+0.4975526040028326*fSkin[2]*coeff[20]-0.1913663861549355*fEdge[2]*coeff[20]+0.3061862178478971*fSkin[5]*coeff[19]-0.3061862178478971*fEdge[5]*coeff[19]-0.2651650429449552*fSkin[1]*coeff[19]-0.2651650429449552*fEdge[1]*coeff[19]-0.5082329989778307*fSkin[7]*coeff[18]-0.1546796083845571*fEdge[7]*coeff[18]+0.4975526040028326*fSkin[4]*coeff[18]-0.1913663861549355*fEdge[4]*coeff[18]-0.5082329989778307*fSkin[3]*coeff[17]-0.1546796083845571*fEdge[3]*coeff[17]+0.4975526040028326*fSkin[0]*coeff[17]-0.1913663861549355*fEdge[0]*coeff[17]-0.5082329989778307*fSkin[5]*coeff[16]-0.1546796083845571*fEdge[5]*coeff[16]+0.4975526040028326*fSkin[1]*coeff[16]-0.1913663861549355*fEdge[1]*coeff[16]; 
  edgeSurf_incr[6] = 0.3061862178478971*fSkin[5]*coeff[23]-0.3061862178478971*fEdge[5]*coeff[23]-0.2651650429449552*fSkin[1]*coeff[23]-0.2651650429449552*fEdge[1]*coeff[23]+0.3061862178478971*fSkin[3]*coeff[22]-0.3061862178478971*fEdge[3]*coeff[22]-0.2651650429449552*fSkin[0]*coeff[22]-0.2651650429449552*fEdge[0]*coeff[22]+0.3061862178478971*fSkin[7]*coeff[21]-0.3061862178478971*fEdge[7]*coeff[21]-0.2651650429449552*fSkin[4]*coeff[21]-0.2651650429449552*fEdge[4]*coeff[21]-0.5082329989778307*fSkin[5]*coeff[20]-0.1546796083845571*fEdge[5]*coeff[20]+0.4975526040028326*fSkin[1]*coeff[20]-0.1913663861549355*fEdge[1]*coeff[20]+0.3061862178478971*fSkin[6]*coeff[19]-0.3061862178478971*fEdge[6]*coeff[19]-0.2651650429449552*fSkin[2]*coeff[19]-0.2651650429449552*fEdge[2]*coeff[19]-0.5082329989778307*fSkin[3]*coeff[18]-0.1546796083845571*fEdge[3]*coeff[18]+0.4975526040028326*fSkin[0]*coeff[18]-0.1913663861549355*fEdge[0]*coeff[18]-0.5082329989778307*fSkin[7]*coeff[17]-0.1546796083845571*fEdge[7]*coeff[17]+0.4975526040028326*fSkin[4]*coeff[17]-0.1913663861549355*fEdge[4]*coeff[17]-0.5082329989778307*fSkin[6]*coeff[16]-0.1546796083845571*fEdge[6]*coeff[16]+0.4975526040028326*fSkin[2]*coeff[16]-0.1913663861549355*fEdge[2]*coeff[16]; 
  edgeSurf_incr[7] = 0.3061862178478971*fSkin[3]*coeff[23]-0.3061862178478971*fEdge[3]*coeff[23]-0.2651650429449552*fSkin[0]*coeff[23]-0.2651650429449552*fEdge[0]*coeff[23]+0.3061862178478971*fSkin[5]*coeff[22]-0.3061862178478971*fEdge[5]*coeff[22]-0.2651650429449552*fSkin[1]*coeff[22]-0.2651650429449552*fEdge[1]*coeff[22]+0.3061862178478971*fSkin[6]*coeff[21]-0.3061862178478971*fEdge[6]*coeff[21]-0.2651650429449552*fSkin[2]*coeff[21]-0.2651650429449552*fEdge[2]*coeff[21]-0.5082329989778307*fSkin[3]*coeff[20]-0.1546796083845571*fEdge[3]*coeff[20]+0.4975526040028326*fSkin[0]*coeff[20]-0.1913663861549355*fEdge[0]*coeff[20]+0.3061862178478971*fSkin[7]*coeff[19]-0.3061862178478971*fEdge[7]*coeff[19]-0.2651650429449552*fSkin[4]*coeff[19]-0.2651650429449552*fEdge[4]*coeff[19]-0.5082329989778307*fSkin[5]*coeff[18]-0.1546796083845571*fEdge[5]*coeff[18]+0.4975526040028326*fSkin[1]*coeff[18]-0.1913663861549355*fEdge[1]*coeff[18]-0.5082329989778307*fSkin[6]*coeff[17]-0.1546796083845571*fEdge[6]*coeff[17]+0.4975526040028326*fSkin[2]*coeff[17]-0.1913663861549355*fEdge[2]*coeff[17]-0.5082329989778307*fSkin[7]*coeff[16]-0.1546796083845571*fEdge[7]*coeff[16]+0.4975526040028326*fSkin[4]*coeff[16]-0.1913663861549355*fEdge[4]*coeff[16]; 

  boundSurf_incr[3] = (-0.6123724356957944*fSkin[7]*coeff[23])-0.5303300858899105*fSkin[4]*coeff[23]-0.6123724356957944*fSkin[6]*coeff[22]-0.5303300858899105*fSkin[2]*coeff[22]-0.6123724356957944*fSkin[5]*coeff[21]-0.5303300858899105*fSkin[1]*coeff[21]-0.3535533905932737*fSkin[7]*coeff[20]-0.3061862178478971*fSkin[4]*coeff[20]-0.6123724356957944*fSkin[3]*coeff[19]-0.5303300858899105*fSkin[0]*coeff[19]-0.3535533905932737*fSkin[6]*coeff[18]-0.3061862178478971*fSkin[2]*coeff[18]-0.3535533905932737*fSkin[5]*coeff[17]-0.3061862178478971*fSkin[1]*coeff[17]-0.3535533905932737*fSkin[3]*coeff[16]-0.3061862178478971*fSkin[0]*coeff[16]; 
  boundSurf_incr[5] = (-0.6123724356957944*fSkin[6]*coeff[23])-0.5303300858899105*fSkin[2]*coeff[23]-0.6123724356957944*fSkin[7]*coeff[22]-0.5303300858899105*fSkin[4]*coeff[22]-0.6123724356957944*fSkin[3]*coeff[21]-0.5303300858899105*fSkin[0]*coeff[21]-0.3535533905932737*fSkin[6]*coeff[20]-0.3061862178478971*fSkin[2]*coeff[20]-0.6123724356957944*fSkin[5]*coeff[19]-0.5303300858899105*fSkin[1]*coeff[19]-0.3535533905932737*fSkin[7]*coeff[18]-0.3061862178478971*fSkin[4]*coeff[18]-0.3535533905932737*fSkin[3]*coeff[17]-0.3061862178478971*fSkin[0]*coeff[17]-0.3535533905932737*fSkin[5]*coeff[16]-0.3061862178478971*fSkin[1]*coeff[16]; 
  boundSurf_incr[6] = (-0.6123724356957944*fSkin[5]*coeff[23])-0.5303300858899105*fSkin[1]*coeff[23]-0.6123724356957944*fSkin[3]*coeff[22]-0.5303300858899105*fSkin[0]*coeff[22]-0.6123724356957944*fSkin[7]*coeff[21]-0.5303300858899105*fSkin[4]*coeff[21]-0.3535533905932737*fSkin[5]*coeff[20]-0.3061862178478971*fSkin[1]*coeff[20]-0.6123724356957944*fSkin[6]*coeff[19]-0.5303300858899105*fSkin[2]*coeff[19]-0.3535533905932737*fSkin[3]*coeff[18]-0.3061862178478971*fSkin[0]*coeff[18]-0.3535533905932737*fSkin[7]*coeff[17]-0.3061862178478971*fSkin[4]*coeff[17]-0.3535533905932737*fSkin[6]*coeff[16]-0.3061862178478971*fSkin[2]*coeff[16]; 
  boundSurf_incr[7] = (-0.6123724356957944*fSkin[3]*coeff[23])-0.5303300858899105*fSkin[0]*coeff[23]-0.6123724356957944*fSkin[5]*coeff[22]-0.5303300858899105*fSkin[1]*coeff[22]-0.6123724356957944*fSkin[6]*coeff[21]-0.5303300858899105*fSkin[2]*coeff[21]-0.3535533905932737*fSkin[3]*coeff[20]-0.3061862178478971*fSkin[0]*coeff[20]-0.6123724356957944*fSkin[7]*coeff[19]-0.5303300858899105*fSkin[4]*coeff[19]-0.3535533905932737*fSkin[5]*coeff[18]-0.3061862178478971*fSkin[1]*coeff[18]-0.3535533905932737*fSkin[6]*coeff[17]-0.3061862178478971*fSkin[2]*coeff[17]-0.3535533905932737*fSkin[7]*coeff[16]-0.3061862178478971*fSkin[4]*coeff[16]; 

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*Jfac; 
  out[4] += (vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*Jfac; 
  out[5] += (vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*Jfac; 
  out[6] += (vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*Jfac; 
  out[7] += (vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*Jfac; 

  }

  return 0.;
}
