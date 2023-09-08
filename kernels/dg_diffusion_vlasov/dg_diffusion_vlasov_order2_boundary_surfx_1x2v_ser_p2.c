#include <gkyl_dg_diffusion_vlasov_kernels.h>

GKYL_CU_DH double dg_diffusion_vlasov_order2_boundary_surfx_1x2v_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],2.);

  double vol_incr[20] = {0.0}; 
  vol_incr[7] = 6.708203932499369*coeff[0]*fSkin[0]; 
  vol_incr[11] = 6.708203932499369*coeff[0]*fSkin[2]; 
  vol_incr[13] = 6.708203932499369*coeff[0]*fSkin[3]; 
  vol_incr[17] = 6.708203932499369*coeff[0]*fSkin[6]; 

  double edgeSurf_incr[20] = {0.0}; 
  double boundSurf_incr[20] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = (-0.6708203932499369*coeff[0]*fSkin[7])+0.6708203932499369*coeff[0]*fEdge[7]-1.190784930203603*coeff[0]*fSkin[1]-1.190784930203603*coeff[0]*fEdge[1]-0.9375*coeff[0]*fSkin[0]+0.9375*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = (-1.585502557353661*coeff[0]*fSkin[7])+0.7382874503707886*coeff[0]*fEdge[7]-2.671875*coeff[0]*fSkin[1]-1.453125*coeff[0]*fEdge[1]-2.056810333988042*coeff[0]*fSkin[0]+1.190784930203603*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = (-0.6708203932499369*coeff[0]*fSkin[11])+0.6708203932499369*coeff[0]*fEdge[11]-1.190784930203603*coeff[0]*fSkin[4]-1.190784930203603*coeff[0]*fEdge[4]-0.9375*coeff[0]*fSkin[2]+0.9375*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = (-0.6708203932499369*coeff[0]*fSkin[13])+0.6708203932499369*coeff[0]*fEdge[13]-1.190784930203603*coeff[0]*fSkin[5]-1.190784930203603*coeff[0]*fEdge[5]-0.9375*coeff[0]*fSkin[3]+0.9375*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = (-1.585502557353661*coeff[0]*fSkin[11])+0.7382874503707888*coeff[0]*fEdge[11]-2.671875*coeff[0]*fSkin[4]-1.453125*coeff[0]*fEdge[4]-2.056810333988042*coeff[0]*fSkin[2]+1.190784930203603*coeff[0]*fEdge[2]; 
  edgeSurf_incr[5] = (-1.585502557353661*coeff[0]*fSkin[13])+0.7382874503707888*coeff[0]*fEdge[13]-2.671875*coeff[0]*fSkin[5]-1.453125*coeff[0]*fEdge[5]-2.056810333988042*coeff[0]*fSkin[3]+1.190784930203603*coeff[0]*fEdge[3]; 
  edgeSurf_incr[6] = (-0.6708203932499369*coeff[0]*fSkin[17])+0.6708203932499369*coeff[0]*fEdge[17]-1.190784930203603*coeff[0]*fSkin[10]-1.190784930203603*coeff[0]*fEdge[10]-0.9375*coeff[0]*fSkin[6]+0.9375*coeff[0]*fEdge[6]; 
  edgeSurf_incr[7] = (-3.140625*coeff[0]*fSkin[7])-0.140625*coeff[0]*fEdge[7]-5.022775277112744*coeff[0]*fSkin[1]-0.3025768239224549*coeff[0]*fEdge[1]-3.773364712030896*coeff[0]*fSkin[0]+0.4192627457812108*coeff[0]*fEdge[0]; 
  edgeSurf_incr[8] = (-1.190784930203603*coeff[0]*fSkin[12])-1.190784930203603*coeff[0]*fEdge[12]-0.9375*coeff[0]*fSkin[8]+0.9375*coeff[0]*fEdge[8]; 
  edgeSurf_incr[9] = (-1.190784930203603*coeff[0]*fSkin[15])-1.190784930203603*coeff[0]*fEdge[15]-0.9375*coeff[0]*fSkin[9]+0.9375*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = (-1.585502557353661*coeff[0]*fSkin[17])+0.7382874503707886*coeff[0]*fEdge[17]-2.671875*coeff[0]*fSkin[10]-1.453125*coeff[0]*fEdge[10]-2.056810333988042*coeff[0]*fSkin[6]+1.190784930203603*coeff[0]*fEdge[6]; 
  edgeSurf_incr[11] = (-3.140625*coeff[0]*fSkin[11])-0.140625*coeff[0]*fEdge[11]-5.022775277112744*coeff[0]*fSkin[4]-0.3025768239224544*coeff[0]*fEdge[4]-3.773364712030894*coeff[0]*fSkin[2]+0.4192627457812105*coeff[0]*fEdge[2]; 
  edgeSurf_incr[12] = (-2.671875*coeff[0]*fSkin[12])-1.453125*coeff[0]*fEdge[12]-2.056810333988042*coeff[0]*fSkin[8]+1.190784930203603*coeff[0]*fEdge[8]; 
  edgeSurf_incr[13] = (-3.140625*coeff[0]*fSkin[13])-0.140625*coeff[0]*fEdge[13]-5.022775277112744*coeff[0]*fSkin[5]-0.3025768239224544*coeff[0]*fEdge[5]-3.773364712030894*coeff[0]*fSkin[3]+0.4192627457812105*coeff[0]*fEdge[3]; 
  edgeSurf_incr[14] = (-1.190784930203603*coeff[0]*fSkin[18])-1.190784930203603*coeff[0]*fEdge[18]-0.9375*coeff[0]*fSkin[14]+0.9375*coeff[0]*fEdge[14]; 
  edgeSurf_incr[15] = (-2.671875*coeff[0]*fSkin[15])-1.453125*coeff[0]*fEdge[15]-2.056810333988042*coeff[0]*fSkin[9]+1.190784930203603*coeff[0]*fEdge[9]; 
  edgeSurf_incr[16] = (-1.190784930203603*coeff[0]*fSkin[19])-1.190784930203603*coeff[0]*fEdge[19]-0.9375*coeff[0]*fSkin[16]+0.9375*coeff[0]*fEdge[16]; 
  edgeSurf_incr[17] = (-3.140625*coeff[0]*fSkin[17])-0.140625*coeff[0]*fEdge[17]-5.022775277112744*coeff[0]*fSkin[10]-0.3025768239224549*coeff[0]*fEdge[10]-3.773364712030896*coeff[0]*fSkin[6]+0.4192627457812108*coeff[0]*fEdge[6]; 
  edgeSurf_incr[18] = (-2.671875*coeff[0]*fSkin[18])-1.453125*coeff[0]*fEdge[18]-2.056810333988042*coeff[0]*fSkin[14]+1.190784930203603*coeff[0]*fEdge[14]; 
  edgeSurf_incr[19] = (-2.671875*coeff[0]*fSkin[19])-1.453125*coeff[0]*fEdge[19]-2.056810333988042*coeff[0]*fSkin[16]+1.190784930203603*coeff[0]*fEdge[16]; 

  boundSurf_incr[1] = 0.9682458365518543*coeff[0]*fSkin[7]-1.25*coeff[0]*fSkin[1]+0.8660254037844386*coeff[0]*fSkin[0]; 
  boundSurf_incr[4] = 0.9682458365518543*coeff[0]*fSkin[11]-1.25*coeff[0]*fSkin[4]+0.8660254037844386*coeff[0]*fSkin[2]; 
  boundSurf_incr[5] = 0.9682458365518543*coeff[0]*fSkin[13]-1.25*coeff[0]*fSkin[5]+0.8660254037844386*coeff[0]*fSkin[3]; 
  boundSurf_incr[7] = (-3.75*coeff[0]*fSkin[7])+4.841229182759272*coeff[0]*fSkin[1]-3.354101966249685*coeff[0]*fSkin[0]; 
  boundSurf_incr[10] = 0.9682458365518543*coeff[0]*fSkin[17]-1.25*coeff[0]*fSkin[10]+0.8660254037844386*coeff[0]*fSkin[6]; 
  boundSurf_incr[11] = (-3.75*coeff[0]*fSkin[11])+4.841229182759271*coeff[0]*fSkin[4]-3.354101966249684*coeff[0]*fSkin[2]; 
  boundSurf_incr[12] = 0.8660254037844387*coeff[0]*fSkin[8]-1.25*coeff[0]*fSkin[12]; 
  boundSurf_incr[13] = (-3.75*coeff[0]*fSkin[13])+4.841229182759271*coeff[0]*fSkin[5]-3.354101966249684*coeff[0]*fSkin[3]; 
  boundSurf_incr[15] = 0.8660254037844387*coeff[0]*fSkin[9]-1.25*coeff[0]*fSkin[15]; 
  boundSurf_incr[17] = (-3.75*coeff[0]*fSkin[17])+4.841229182759272*coeff[0]*fSkin[10]-3.354101966249685*coeff[0]*fSkin[6]; 
  boundSurf_incr[18] = 0.8660254037844387*coeff[0]*fSkin[14]-1.25*coeff[0]*fSkin[18]; 
  boundSurf_incr[19] = 0.8660254037844387*coeff[0]*fSkin[16]-1.25*coeff[0]*fSkin[19]; 

  } else { 

  edgeSurf_incr[0] = (-0.6708203932499369*coeff[0]*fSkin[7])+0.6708203932499369*coeff[0]*fEdge[7]+1.190784930203603*coeff[0]*fSkin[1]+1.190784930203603*coeff[0]*fEdge[1]-0.9375*coeff[0]*fSkin[0]+0.9375*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 1.585502557353661*coeff[0]*fSkin[7]-0.7382874503707886*coeff[0]*fEdge[7]-2.671875*coeff[0]*fSkin[1]-1.453125*coeff[0]*fEdge[1]+2.056810333988042*coeff[0]*fSkin[0]-1.190784930203603*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = (-0.6708203932499369*coeff[0]*fSkin[11])+0.6708203932499369*coeff[0]*fEdge[11]+1.190784930203603*coeff[0]*fSkin[4]+1.190784930203603*coeff[0]*fEdge[4]-0.9375*coeff[0]*fSkin[2]+0.9375*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = (-0.6708203932499369*coeff[0]*fSkin[13])+0.6708203932499369*coeff[0]*fEdge[13]+1.190784930203603*coeff[0]*fSkin[5]+1.190784930203603*coeff[0]*fEdge[5]-0.9375*coeff[0]*fSkin[3]+0.9375*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = 1.585502557353661*coeff[0]*fSkin[11]-0.7382874503707888*coeff[0]*fEdge[11]-2.671875*coeff[0]*fSkin[4]-1.453125*coeff[0]*fEdge[4]+2.056810333988042*coeff[0]*fSkin[2]-1.190784930203603*coeff[0]*fEdge[2]; 
  edgeSurf_incr[5] = 1.585502557353661*coeff[0]*fSkin[13]-0.7382874503707888*coeff[0]*fEdge[13]-2.671875*coeff[0]*fSkin[5]-1.453125*coeff[0]*fEdge[5]+2.056810333988042*coeff[0]*fSkin[3]-1.190784930203603*coeff[0]*fEdge[3]; 
  edgeSurf_incr[6] = (-0.6708203932499369*coeff[0]*fSkin[17])+0.6708203932499369*coeff[0]*fEdge[17]+1.190784930203603*coeff[0]*fSkin[10]+1.190784930203603*coeff[0]*fEdge[10]-0.9375*coeff[0]*fSkin[6]+0.9375*coeff[0]*fEdge[6]; 
  edgeSurf_incr[7] = (-3.140625*coeff[0]*fSkin[7])-0.140625*coeff[0]*fEdge[7]+5.022775277112744*coeff[0]*fSkin[1]+0.3025768239224549*coeff[0]*fEdge[1]-3.773364712030896*coeff[0]*fSkin[0]+0.4192627457812108*coeff[0]*fEdge[0]; 
  edgeSurf_incr[8] = 1.190784930203603*coeff[0]*fSkin[12]+1.190784930203603*coeff[0]*fEdge[12]-0.9375*coeff[0]*fSkin[8]+0.9375*coeff[0]*fEdge[8]; 
  edgeSurf_incr[9] = 1.190784930203603*coeff[0]*fSkin[15]+1.190784930203603*coeff[0]*fEdge[15]-0.9375*coeff[0]*fSkin[9]+0.9375*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = 1.585502557353661*coeff[0]*fSkin[17]-0.7382874503707886*coeff[0]*fEdge[17]-2.671875*coeff[0]*fSkin[10]-1.453125*coeff[0]*fEdge[10]+2.056810333988042*coeff[0]*fSkin[6]-1.190784930203603*coeff[0]*fEdge[6]; 
  edgeSurf_incr[11] = (-3.140625*coeff[0]*fSkin[11])-0.140625*coeff[0]*fEdge[11]+5.022775277112744*coeff[0]*fSkin[4]+0.3025768239224544*coeff[0]*fEdge[4]-3.773364712030894*coeff[0]*fSkin[2]+0.4192627457812105*coeff[0]*fEdge[2]; 
  edgeSurf_incr[12] = (-2.671875*coeff[0]*fSkin[12])-1.453125*coeff[0]*fEdge[12]+2.056810333988042*coeff[0]*fSkin[8]-1.190784930203603*coeff[0]*fEdge[8]; 
  edgeSurf_incr[13] = (-3.140625*coeff[0]*fSkin[13])-0.140625*coeff[0]*fEdge[13]+5.022775277112744*coeff[0]*fSkin[5]+0.3025768239224544*coeff[0]*fEdge[5]-3.773364712030894*coeff[0]*fSkin[3]+0.4192627457812105*coeff[0]*fEdge[3]; 
  edgeSurf_incr[14] = 1.190784930203603*coeff[0]*fSkin[18]+1.190784930203603*coeff[0]*fEdge[18]-0.9375*coeff[0]*fSkin[14]+0.9375*coeff[0]*fEdge[14]; 
  edgeSurf_incr[15] = (-2.671875*coeff[0]*fSkin[15])-1.453125*coeff[0]*fEdge[15]+2.056810333988042*coeff[0]*fSkin[9]-1.190784930203603*coeff[0]*fEdge[9]; 
  edgeSurf_incr[16] = 1.190784930203603*coeff[0]*fSkin[19]+1.190784930203603*coeff[0]*fEdge[19]-0.9375*coeff[0]*fSkin[16]+0.9375*coeff[0]*fEdge[16]; 
  edgeSurf_incr[17] = (-3.140625*coeff[0]*fSkin[17])-0.140625*coeff[0]*fEdge[17]+5.022775277112744*coeff[0]*fSkin[10]+0.3025768239224549*coeff[0]*fEdge[10]-3.773364712030896*coeff[0]*fSkin[6]+0.4192627457812108*coeff[0]*fEdge[6]; 
  edgeSurf_incr[18] = (-2.671875*coeff[0]*fSkin[18])-1.453125*coeff[0]*fEdge[18]+2.056810333988042*coeff[0]*fSkin[14]-1.190784930203603*coeff[0]*fEdge[14]; 
  edgeSurf_incr[19] = (-2.671875*coeff[0]*fSkin[19])-1.453125*coeff[0]*fEdge[19]+2.056810333988042*coeff[0]*fSkin[16]-1.190784930203603*coeff[0]*fEdge[16]; 

  boundSurf_incr[1] = (-0.9682458365518543*coeff[0]*fSkin[7])-1.25*coeff[0]*fSkin[1]-0.8660254037844386*coeff[0]*fSkin[0]; 
  boundSurf_incr[4] = (-0.9682458365518543*coeff[0]*fSkin[11])-1.25*coeff[0]*fSkin[4]-0.8660254037844386*coeff[0]*fSkin[2]; 
  boundSurf_incr[5] = (-0.9682458365518543*coeff[0]*fSkin[13])-1.25*coeff[0]*fSkin[5]-0.8660254037844386*coeff[0]*fSkin[3]; 
  boundSurf_incr[7] = (-3.75*coeff[0]*fSkin[7])-4.841229182759272*coeff[0]*fSkin[1]-3.354101966249685*coeff[0]*fSkin[0]; 
  boundSurf_incr[10] = (-0.9682458365518543*coeff[0]*fSkin[17])-1.25*coeff[0]*fSkin[10]-0.8660254037844386*coeff[0]*fSkin[6]; 
  boundSurf_incr[11] = (-3.75*coeff[0]*fSkin[11])-4.841229182759271*coeff[0]*fSkin[4]-3.354101966249684*coeff[0]*fSkin[2]; 
  boundSurf_incr[12] = (-1.25*coeff[0]*fSkin[12])-0.8660254037844387*coeff[0]*fSkin[8]; 
  boundSurf_incr[13] = (-3.75*coeff[0]*fSkin[13])-4.841229182759271*coeff[0]*fSkin[5]-3.354101966249684*coeff[0]*fSkin[3]; 
  boundSurf_incr[15] = (-1.25*coeff[0]*fSkin[15])-0.8660254037844387*coeff[0]*fSkin[9]; 
  boundSurf_incr[17] = (-3.75*coeff[0]*fSkin[17])-4.841229182759272*coeff[0]*fSkin[10]-3.354101966249685*coeff[0]*fSkin[6]; 
  boundSurf_incr[18] = (-1.25*coeff[0]*fSkin[18])-0.8660254037844387*coeff[0]*fSkin[14]; 
  boundSurf_incr[19] = (-1.25*coeff[0]*fSkin[19])-0.8660254037844387*coeff[0]*fSkin[16]; 

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
  out[12] += (vol_incr[12]+edgeSurf_incr[12]+boundSurf_incr[12])*Jfac; 
  out[13] += (vol_incr[13]+edgeSurf_incr[13]+boundSurf_incr[13])*Jfac; 
  out[14] += (vol_incr[14]+edgeSurf_incr[14]+boundSurf_incr[14])*Jfac; 
  out[15] += (vol_incr[15]+edgeSurf_incr[15]+boundSurf_incr[15])*Jfac; 
  out[16] += (vol_incr[16]+edgeSurf_incr[16]+boundSurf_incr[16])*Jfac; 
  out[17] += (vol_incr[17]+edgeSurf_incr[17]+boundSurf_incr[17])*Jfac; 
  out[18] += (vol_incr[18]+edgeSurf_incr[18]+boundSurf_incr[18])*Jfac; 
  out[19] += (vol_incr[19]+edgeSurf_incr[19]+boundSurf_incr[19])*Jfac; 

  }

  return 0.;
}

GKYL_CU_DH double dg_diffusion_vlasov_order2_boundary_surfx_1x2v_ser_p2_varcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],2.);

  double vol_incr[20] = {0.0}; 
  vol_incr[1] = 4.743416490252569*fSkin[1]*coeff[2]+2.121320343559642*fSkin[0]*coeff[1]; 
  vol_incr[4] = 4.743416490252569*coeff[2]*fSkin[4]+2.121320343559642*coeff[1]*fSkin[2]; 
  vol_incr[5] = 4.743416490252569*coeff[2]*fSkin[5]+2.121320343559642*coeff[1]*fSkin[3]; 
  vol_incr[7] = 14.23024947075771*coeff[2]*fSkin[7]+10.60660171779821*fSkin[0]*coeff[2]+9.48683298050514*coeff[1]*fSkin[1]+4.743416490252569*coeff[0]*fSkin[0]; 
  vol_incr[10] = 4.743416490252569*coeff[2]*fSkin[10]+2.121320343559642*coeff[1]*fSkin[6]; 
  vol_incr[11] = 14.23024947075771*coeff[2]*fSkin[11]+9.48683298050514*coeff[1]*fSkin[4]+10.60660171779821*coeff[2]*fSkin[2]+4.743416490252569*coeff[0]*fSkin[2]; 
  vol_incr[12] = 4.743416490252569*coeff[2]*fSkin[12]+2.121320343559642*coeff[1]*fSkin[8]; 
  vol_incr[13] = 14.23024947075771*coeff[2]*fSkin[13]+9.48683298050514*coeff[1]*fSkin[5]+10.60660171779821*coeff[2]*fSkin[3]+4.743416490252569*coeff[0]*fSkin[3]; 
  vol_incr[15] = 4.743416490252569*coeff[2]*fSkin[15]+2.121320343559642*coeff[1]*fSkin[9]; 
  vol_incr[17] = 14.23024947075771*coeff[2]*fSkin[17]+9.48683298050514*coeff[1]*fSkin[10]+10.60660171779821*coeff[2]*fSkin[6]+4.743416490252569*coeff[0]*fSkin[6]; 
  vol_incr[18] = 4.743416490252569*coeff[2]*fSkin[18]+2.121320343559642*coeff[1]*fSkin[14]; 
  vol_incr[19] = 4.743416490252569*coeff[2]*fSkin[19]+2.121320343559642*coeff[1]*fSkin[16]; 

  double edgeSurf_incr[20] = {0.0}; 
  double boundSurf_incr[20] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = 0.5303300858899105*coeff[2]*fSkin[7]-0.4743416490252568*coeff[0]*fSkin[7]-0.5303300858899105*coeff[2]*fEdge[7]+0.4743416490252568*coeff[0]*fEdge[7]+0.9413981457120035*fSkin[1]*coeff[2]+0.9413981457120035*fEdge[1]*coeff[2]+0.7411588266019635*fSkin[0]*coeff[2]-0.7411588266019635*fEdge[0]*coeff[2]-0.8420120990817169*coeff[0]*fSkin[1]-0.8420120990817169*coeff[0]*fEdge[1]-0.6629126073623879*coeff[0]*fSkin[0]+0.6629126073623879*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 0.2487763020014165*coeff[2]*fSkin[7]-0.5188111786213743*coeff[1]*fSkin[7]-1.121119609893386*coeff[0]*fSkin[7]-1.588341005085966*coeff[2]*fEdge[7]-0.5188111786213743*coeff[1]*fEdge[7]+0.5220480626221115*coeff[0]*fEdge[7]+0.6670429439417671*fSkin[1]*coeff[2]+2.594055893106872*fEdge[1]*coeff[2]+0.5990715472712751*fSkin[0]*coeff[2]-1.96837794103419*fEdge[0]*coeff[2]-0.7463289060042488*coeff[1]*fSkin[1]-1.889300930982805*coeff[0]*fSkin[1]+0.7463289060042488*coeff[1]*fEdge[1]-1.027514541411701*coeff[0]*fEdge[1]-0.5303300858899105*fSkin[0]*coeff[1]-0.5303300858899105*fEdge[0]*coeff[1]-1.454384534777511*coeff[0]*fSkin[0]+0.8420120990817168*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 0.5303300858899104*coeff[2]*fSkin[11]-0.4743416490252568*coeff[0]*fSkin[11]-0.5303300858899104*coeff[2]*fEdge[11]+0.4743416490252568*coeff[0]*fEdge[11]+0.9413981457120035*coeff[2]*fSkin[4]-0.8420120990817169*coeff[0]*fSkin[4]+0.9413981457120035*coeff[2]*fEdge[4]-0.8420120990817169*coeff[0]*fEdge[4]+0.7411588266019635*coeff[2]*fSkin[2]-0.6629126073623879*coeff[0]*fSkin[2]-0.7411588266019635*coeff[2]*fEdge[2]+0.6629126073623879*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 0.5303300858899104*coeff[2]*fSkin[13]-0.4743416490252568*coeff[0]*fSkin[13]-0.5303300858899104*coeff[2]*fEdge[13]+0.4743416490252568*coeff[0]*fEdge[13]+0.9413981457120035*coeff[2]*fSkin[5]-0.8420120990817169*coeff[0]*fSkin[5]+0.9413981457120035*coeff[2]*fEdge[5]-0.8420120990817169*coeff[0]*fEdge[5]+0.7411588266019635*coeff[2]*fSkin[3]-0.6629126073623879*coeff[0]*fSkin[3]-0.7411588266019635*coeff[2]*fEdge[3]+0.6629126073623879*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = 0.2487763020014169*coeff[2]*fSkin[11]-0.5188111786213743*coeff[1]*fSkin[11]-1.121119609893386*coeff[0]*fSkin[11]-1.588341005085966*coeff[2]*fEdge[11]-0.5188111786213743*coeff[1]*fEdge[11]+0.5220480626221116*coeff[0]*fEdge[11]+0.6670429439417671*coeff[2]*fSkin[4]-0.7463289060042488*coeff[1]*fSkin[4]-1.889300930982805*coeff[0]*fSkin[4]+2.594055893106872*coeff[2]*fEdge[4]+0.7463289060042488*coeff[1]*fEdge[4]-1.027514541411701*coeff[0]*fEdge[4]+0.5990715472712751*coeff[2]*fSkin[2]-0.5303300858899105*coeff[1]*fSkin[2]-1.454384534777511*coeff[0]*fSkin[2]-1.96837794103419*coeff[2]*fEdge[2]-0.5303300858899105*coeff[1]*fEdge[2]+0.8420120990817168*coeff[0]*fEdge[2]; 
  edgeSurf_incr[5] = 0.2487763020014169*coeff[2]*fSkin[13]-0.5188111786213743*coeff[1]*fSkin[13]-1.121119609893386*coeff[0]*fSkin[13]-1.588341005085966*coeff[2]*fEdge[13]-0.5188111786213743*coeff[1]*fEdge[13]+0.5220480626221116*coeff[0]*fEdge[13]+0.6670429439417671*coeff[2]*fSkin[5]-0.7463289060042488*coeff[1]*fSkin[5]-1.889300930982805*coeff[0]*fSkin[5]+2.594055893106872*coeff[2]*fEdge[5]+0.7463289060042488*coeff[1]*fEdge[5]-1.027514541411701*coeff[0]*fEdge[5]+0.5990715472712751*coeff[2]*fSkin[3]-0.5303300858899105*coeff[1]*fSkin[3]-1.454384534777511*coeff[0]*fSkin[3]-1.96837794103419*coeff[2]*fEdge[3]-0.5303300858899105*coeff[1]*fEdge[3]+0.8420120990817168*coeff[0]*fEdge[3]; 
  edgeSurf_incr[6] = 0.5303300858899105*coeff[2]*fSkin[17]-0.4743416490252568*coeff[0]*fSkin[17]-0.5303300858899105*coeff[2]*fEdge[17]+0.4743416490252568*coeff[0]*fEdge[17]+0.9413981457120035*coeff[2]*fSkin[10]-0.8420120990817169*coeff[0]*fSkin[10]+0.9413981457120035*coeff[2]*fEdge[10]-0.8420120990817169*coeff[0]*fEdge[10]+0.7411588266019635*coeff[2]*fSkin[6]-0.6629126073623879*coeff[0]*fSkin[6]-0.7411588266019635*coeff[2]*fEdge[6]+0.6629126073623879*coeff[0]*fEdge[6]; 
  edgeSurf_incr[7] = (-1.40820177054373*coeff[2]*fSkin[7])-2.009347054626824*coeff[1]*fSkin[7]-2.220757234663999*coeff[0]*fSkin[7]-3.779910015670014*coeff[2]*fEdge[7]-2.009347054626824*coeff[1]*fEdge[7]-0.0994368911043575*coeff[0]*fEdge[7]-1.626614282316952*fSkin[1]*coeff[2]+5.836674777725536*fEdge[1]*coeff[2]-0.9943689110435825*fSkin[0]*coeff[2]-4.308931947855521*fEdge[0]*coeff[2]-2.890519423747657*coeff[1]*fSkin[1]-3.551638458822559*coeff[0]*fSkin[1]+2.890519423747657*coeff[1]*fEdge[1]-0.2139541240254559*coeff[0]*fEdge[1]-2.053959590644372*fSkin[0]*coeff[1]-2.053959590644372*fEdge[0]*coeff[1]-2.668171775767069*coeff[0]*fSkin[0]+0.2964635306407852*coeff[0]*fEdge[0]; 
  edgeSurf_incr[8] = 0.9413981457120036*coeff[2]*fSkin[12]-0.842012099081717*coeff[0]*fSkin[12]+0.9413981457120036*coeff[2]*fEdge[12]-0.842012099081717*coeff[0]*fEdge[12]+0.7411588266019635*coeff[2]*fSkin[8]-0.6629126073623879*coeff[0]*fSkin[8]-0.7411588266019635*coeff[2]*fEdge[8]+0.6629126073623879*coeff[0]*fEdge[8]; 
  edgeSurf_incr[9] = 0.9413981457120036*coeff[2]*fSkin[15]-0.842012099081717*coeff[0]*fSkin[15]+0.9413981457120036*coeff[2]*fEdge[15]-0.842012099081717*coeff[0]*fEdge[15]+0.7411588266019635*coeff[2]*fSkin[9]-0.6629126073623879*coeff[0]*fSkin[9]-0.7411588266019635*coeff[2]*fEdge[9]+0.6629126073623879*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = 0.2487763020014165*coeff[2]*fSkin[17]-0.5188111786213743*coeff[1]*fSkin[17]-1.121119609893386*coeff[0]*fSkin[17]-1.588341005085966*coeff[2]*fEdge[17]-0.5188111786213743*coeff[1]*fEdge[17]+0.5220480626221115*coeff[0]*fEdge[17]+0.6670429439417671*coeff[2]*fSkin[10]-0.7463289060042488*coeff[1]*fSkin[10]-1.889300930982805*coeff[0]*fSkin[10]+2.594055893106872*coeff[2]*fEdge[10]+0.7463289060042488*coeff[1]*fEdge[10]-1.027514541411701*coeff[0]*fEdge[10]+0.5990715472712751*coeff[2]*fSkin[6]-0.5303300858899105*coeff[1]*fSkin[6]-1.454384534777511*coeff[0]*fSkin[6]-1.96837794103419*coeff[2]*fEdge[6]-0.5303300858899105*coeff[1]*fEdge[6]+0.8420120990817168*coeff[0]*fEdge[6]; 
  edgeSurf_incr[11] = (-1.40820177054373*coeff[2]*fSkin[11])-2.009347054626824*coeff[1]*fSkin[11]-2.220757234663999*coeff[0]*fSkin[11]-3.779910015670014*coeff[2]*fEdge[11]-2.009347054626824*coeff[1]*fEdge[11]-0.0994368911043575*coeff[0]*fEdge[11]-1.626614282316953*coeff[2]*fSkin[4]-2.890519423747656*coeff[1]*fSkin[4]-3.551638458822559*coeff[0]*fSkin[4]+5.836674777725538*coeff[2]*fEdge[4]+2.890519423747656*coeff[1]*fEdge[4]-0.2139541240254559*coeff[0]*fEdge[4]-0.9943689110435823*coeff[2]*fSkin[2]-2.053959590644372*coeff[1]*fSkin[2]-2.668171775767069*coeff[0]*fSkin[2]-4.308931947855521*coeff[2]*fEdge[2]-2.053959590644372*coeff[1]*fEdge[2]+0.2964635306407852*coeff[0]*fEdge[2]; 
  edgeSurf_incr[12] = 0.6670429439417671*coeff[2]*fSkin[12]-0.7463289060042488*coeff[1]*fSkin[12]-1.889300930982805*coeff[0]*fSkin[12]+2.594055893106872*coeff[2]*fEdge[12]+0.7463289060042488*coeff[1]*fEdge[12]-1.027514541411701*coeff[0]*fEdge[12]+0.599071547271275*coeff[2]*fSkin[8]-0.5303300858899104*coeff[1]*fSkin[8]-1.454384534777511*coeff[0]*fSkin[8]-1.96837794103419*coeff[2]*fEdge[8]-0.5303300858899104*coeff[1]*fEdge[8]+0.842012099081717*coeff[0]*fEdge[8]; 
  edgeSurf_incr[13] = (-1.40820177054373*coeff[2]*fSkin[13])-2.009347054626824*coeff[1]*fSkin[13]-2.220757234663999*coeff[0]*fSkin[13]-3.779910015670014*coeff[2]*fEdge[13]-2.009347054626824*coeff[1]*fEdge[13]-0.0994368911043575*coeff[0]*fEdge[13]-1.626614282316953*coeff[2]*fSkin[5]-2.890519423747656*coeff[1]*fSkin[5]-3.551638458822559*coeff[0]*fSkin[5]+5.836674777725538*coeff[2]*fEdge[5]+2.890519423747656*coeff[1]*fEdge[5]-0.2139541240254559*coeff[0]*fEdge[5]-0.9943689110435823*coeff[2]*fSkin[3]-2.053959590644372*coeff[1]*fSkin[3]-2.668171775767069*coeff[0]*fSkin[3]-4.308931947855521*coeff[2]*fEdge[3]-2.053959590644372*coeff[1]*fEdge[3]+0.2964635306407852*coeff[0]*fEdge[3]; 
  edgeSurf_incr[14] = 0.9413981457120036*coeff[2]*fSkin[18]-0.842012099081717*coeff[0]*fSkin[18]+0.9413981457120036*coeff[2]*fEdge[18]-0.842012099081717*coeff[0]*fEdge[18]+0.7411588266019635*coeff[2]*fSkin[14]-0.6629126073623879*coeff[0]*fSkin[14]-0.7411588266019635*coeff[2]*fEdge[14]+0.6629126073623879*coeff[0]*fEdge[14]; 
  edgeSurf_incr[15] = 0.6670429439417671*coeff[2]*fSkin[15]-0.7463289060042488*coeff[1]*fSkin[15]-1.889300930982805*coeff[0]*fSkin[15]+2.594055893106872*coeff[2]*fEdge[15]+0.7463289060042488*coeff[1]*fEdge[15]-1.027514541411701*coeff[0]*fEdge[15]+0.599071547271275*coeff[2]*fSkin[9]-0.5303300858899104*coeff[1]*fSkin[9]-1.454384534777511*coeff[0]*fSkin[9]-1.96837794103419*coeff[2]*fEdge[9]-0.5303300858899104*coeff[1]*fEdge[9]+0.842012099081717*coeff[0]*fEdge[9]; 
  edgeSurf_incr[16] = 0.9413981457120036*coeff[2]*fSkin[19]-0.842012099081717*coeff[0]*fSkin[19]+0.9413981457120036*coeff[2]*fEdge[19]-0.842012099081717*coeff[0]*fEdge[19]+0.7411588266019635*coeff[2]*fSkin[16]-0.6629126073623879*coeff[0]*fSkin[16]-0.7411588266019635*coeff[2]*fEdge[16]+0.6629126073623879*coeff[0]*fEdge[16]; 
  edgeSurf_incr[17] = (-1.40820177054373*coeff[2]*fSkin[17])-2.009347054626824*coeff[1]*fSkin[17]-2.220757234663999*coeff[0]*fSkin[17]-3.779910015670014*coeff[2]*fEdge[17]-2.009347054626824*coeff[1]*fEdge[17]-0.0994368911043575*coeff[0]*fEdge[17]-1.626614282316952*coeff[2]*fSkin[10]-2.890519423747657*coeff[1]*fSkin[10]-3.551638458822559*coeff[0]*fSkin[10]+5.836674777725536*coeff[2]*fEdge[10]+2.890519423747657*coeff[1]*fEdge[10]-0.2139541240254559*coeff[0]*fEdge[10]-0.9943689110435825*coeff[2]*fSkin[6]-2.053959590644372*coeff[1]*fSkin[6]-2.668171775767069*coeff[0]*fSkin[6]-4.308931947855521*coeff[2]*fEdge[6]-2.053959590644372*coeff[1]*fEdge[6]+0.2964635306407852*coeff[0]*fEdge[6]; 
  edgeSurf_incr[18] = 0.6670429439417671*coeff[2]*fSkin[18]-0.7463289060042488*coeff[1]*fSkin[18]-1.889300930982805*coeff[0]*fSkin[18]+2.594055893106872*coeff[2]*fEdge[18]+0.7463289060042488*coeff[1]*fEdge[18]-1.027514541411701*coeff[0]*fEdge[18]+0.599071547271275*coeff[2]*fSkin[14]-0.5303300858899104*coeff[1]*fSkin[14]-1.454384534777511*coeff[0]*fSkin[14]-1.96837794103419*coeff[2]*fEdge[14]-0.5303300858899104*coeff[1]*fEdge[14]+0.842012099081717*coeff[0]*fEdge[14]; 
  edgeSurf_incr[19] = 0.6670429439417671*coeff[2]*fSkin[19]-0.7463289060042488*coeff[1]*fSkin[19]-1.889300930982805*coeff[0]*fSkin[19]+2.594055893106872*coeff[2]*fEdge[19]+0.7463289060042488*coeff[1]*fEdge[19]-1.027514541411701*coeff[0]*fEdge[19]+0.599071547271275*coeff[2]*fSkin[16]-0.5303300858899104*coeff[1]*fSkin[16]-1.454384534777511*coeff[0]*fSkin[16]-1.96837794103419*coeff[2]*fEdge[16]-0.5303300858899104*coeff[1]*fEdge[16]+0.842012099081717*coeff[0]*fEdge[16]; 

  boundSurf_incr[1] = 1.530931089239486*coeff[2]*fSkin[7]-1.185854122563142*coeff[1]*fSkin[7]+0.6846531968814573*coeff[0]*fSkin[7]-1.976423537605237*fSkin[1]*coeff[2]+1.369306393762915*fSkin[0]*coeff[2]+1.530931089239486*coeff[1]*fSkin[1]-0.883883476483184*coeff[0]*fSkin[1]-1.060660171779821*fSkin[0]*coeff[1]+0.6123724356957944*coeff[0]*fSkin[0]; 
  boundSurf_incr[4] = 1.530931089239486*coeff[2]*fSkin[11]-1.185854122563142*coeff[1]*fSkin[11]+0.6846531968814574*coeff[0]*fSkin[11]-1.976423537605237*coeff[2]*fSkin[4]+1.530931089239486*coeff[1]*fSkin[4]-0.883883476483184*coeff[0]*fSkin[4]+1.369306393762915*coeff[2]*fSkin[2]-1.060660171779821*coeff[1]*fSkin[2]+0.6123724356957944*coeff[0]*fSkin[2]; 
  boundSurf_incr[5] = 1.530931089239486*coeff[2]*fSkin[13]-1.185854122563142*coeff[1]*fSkin[13]+0.6846531968814574*coeff[0]*fSkin[13]-1.976423537605237*coeff[2]*fSkin[5]+1.530931089239486*coeff[1]*fSkin[5]-0.883883476483184*coeff[0]*fSkin[5]+1.369306393762915*coeff[2]*fSkin[3]-1.060660171779821*coeff[1]*fSkin[3]+0.6123724356957944*coeff[0]*fSkin[3]; 
  boundSurf_incr[7] = (-5.929270612815711*coeff[2]*fSkin[7])+4.592793267718456*coeff[1]*fSkin[7]-2.651650429449552*coeff[0]*fSkin[7]+7.654655446197428*fSkin[1]*coeff[2]-5.303300858899105*fSkin[0]*coeff[2]-5.929270612815711*coeff[1]*fSkin[1]+3.423265984407287*coeff[0]*fSkin[1]+4.107919181288745*fSkin[0]*coeff[1]-2.371708245126284*coeff[0]*fSkin[0]; 
  boundSurf_incr[10] = 1.530931089239486*coeff[2]*fSkin[17]-1.185854122563142*coeff[1]*fSkin[17]+0.6846531968814573*coeff[0]*fSkin[17]-1.976423537605237*coeff[2]*fSkin[10]+1.530931089239486*coeff[1]*fSkin[10]-0.883883476483184*coeff[0]*fSkin[10]+1.369306393762915*coeff[2]*fSkin[6]-1.060660171779821*coeff[1]*fSkin[6]+0.6123724356957944*coeff[0]*fSkin[6]; 
  boundSurf_incr[11] = (-5.929270612815711*coeff[2]*fSkin[11])+4.592793267718456*coeff[1]*fSkin[11]-2.651650429449552*coeff[0]*fSkin[11]+7.65465544619743*coeff[2]*fSkin[4]-5.929270612815709*coeff[1]*fSkin[4]+3.423265984407287*coeff[0]*fSkin[4]-5.303300858899106*coeff[2]*fSkin[2]+4.107919181288745*coeff[1]*fSkin[2]-2.371708245126284*coeff[0]*fSkin[2]; 
  boundSurf_incr[12] = (-1.976423537605237*coeff[2]*fSkin[12])+1.530931089239486*coeff[1]*fSkin[12]-0.883883476483184*coeff[0]*fSkin[12]+1.369306393762915*coeff[2]*fSkin[8]-1.060660171779821*coeff[1]*fSkin[8]+0.6123724356957944*coeff[0]*fSkin[8]; 
  boundSurf_incr[13] = (-5.929270612815711*coeff[2]*fSkin[13])+4.592793267718456*coeff[1]*fSkin[13]-2.651650429449552*coeff[0]*fSkin[13]+7.65465544619743*coeff[2]*fSkin[5]-5.929270612815709*coeff[1]*fSkin[5]+3.423265984407287*coeff[0]*fSkin[5]-5.303300858899106*coeff[2]*fSkin[3]+4.107919181288745*coeff[1]*fSkin[3]-2.371708245126284*coeff[0]*fSkin[3]; 
  boundSurf_incr[15] = (-1.976423537605237*coeff[2]*fSkin[15])+1.530931089239486*coeff[1]*fSkin[15]-0.883883476483184*coeff[0]*fSkin[15]+1.369306393762915*coeff[2]*fSkin[9]-1.060660171779821*coeff[1]*fSkin[9]+0.6123724356957944*coeff[0]*fSkin[9]; 
  boundSurf_incr[17] = (-5.929270612815711*coeff[2]*fSkin[17])+4.592793267718456*coeff[1]*fSkin[17]-2.651650429449552*coeff[0]*fSkin[17]+7.654655446197428*coeff[2]*fSkin[10]-5.929270612815711*coeff[1]*fSkin[10]+3.423265984407287*coeff[0]*fSkin[10]-5.303300858899105*coeff[2]*fSkin[6]+4.107919181288745*coeff[1]*fSkin[6]-2.371708245126284*coeff[0]*fSkin[6]; 
  boundSurf_incr[18] = (-1.976423537605237*coeff[2]*fSkin[18])+1.530931089239486*coeff[1]*fSkin[18]-0.883883476483184*coeff[0]*fSkin[18]+1.369306393762915*coeff[2]*fSkin[14]-1.060660171779821*coeff[1]*fSkin[14]+0.6123724356957944*coeff[0]*fSkin[14]; 
  boundSurf_incr[19] = (-1.976423537605237*coeff[2]*fSkin[19])+1.530931089239486*coeff[1]*fSkin[19]-0.883883476483184*coeff[0]*fSkin[19]+1.369306393762915*coeff[2]*fSkin[16]-1.060660171779821*coeff[1]*fSkin[16]+0.6123724356957944*coeff[0]*fSkin[16]; 

  } else { 

  edgeSurf_incr[0] = 0.5303300858899105*coeff[2]*fSkin[7]-0.4743416490252568*coeff[0]*fSkin[7]-0.5303300858899105*coeff[2]*fEdge[7]+0.4743416490252568*coeff[0]*fEdge[7]-0.9413981457120035*fSkin[1]*coeff[2]-0.9413981457120035*fEdge[1]*coeff[2]+0.7411588266019635*fSkin[0]*coeff[2]-0.7411588266019635*fEdge[0]*coeff[2]+0.8420120990817169*coeff[0]*fSkin[1]+0.8420120990817169*coeff[0]*fEdge[1]-0.6629126073623879*coeff[0]*fSkin[0]+0.6629126073623879*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = (-0.2487763020014165*coeff[2]*fSkin[7])-0.5188111786213743*coeff[1]*fSkin[7]+1.121119609893386*coeff[0]*fSkin[7]+1.588341005085966*coeff[2]*fEdge[7]-0.5188111786213743*coeff[1]*fEdge[7]-0.5220480626221115*coeff[0]*fEdge[7]+0.6670429439417671*fSkin[1]*coeff[2]+2.594055893106872*fEdge[1]*coeff[2]-0.5990715472712751*fSkin[0]*coeff[2]+1.96837794103419*fEdge[0]*coeff[2]+0.7463289060042488*coeff[1]*fSkin[1]-1.889300930982805*coeff[0]*fSkin[1]-0.7463289060042488*coeff[1]*fEdge[1]-1.027514541411701*coeff[0]*fEdge[1]-0.5303300858899105*fSkin[0]*coeff[1]-0.5303300858899105*fEdge[0]*coeff[1]+1.454384534777511*coeff[0]*fSkin[0]-0.8420120990817168*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 0.5303300858899104*coeff[2]*fSkin[11]-0.4743416490252568*coeff[0]*fSkin[11]-0.5303300858899104*coeff[2]*fEdge[11]+0.4743416490252568*coeff[0]*fEdge[11]-0.9413981457120035*coeff[2]*fSkin[4]+0.8420120990817169*coeff[0]*fSkin[4]-0.9413981457120035*coeff[2]*fEdge[4]+0.8420120990817169*coeff[0]*fEdge[4]+0.7411588266019635*coeff[2]*fSkin[2]-0.6629126073623879*coeff[0]*fSkin[2]-0.7411588266019635*coeff[2]*fEdge[2]+0.6629126073623879*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 0.5303300858899104*coeff[2]*fSkin[13]-0.4743416490252568*coeff[0]*fSkin[13]-0.5303300858899104*coeff[2]*fEdge[13]+0.4743416490252568*coeff[0]*fEdge[13]-0.9413981457120035*coeff[2]*fSkin[5]+0.8420120990817169*coeff[0]*fSkin[5]-0.9413981457120035*coeff[2]*fEdge[5]+0.8420120990817169*coeff[0]*fEdge[5]+0.7411588266019635*coeff[2]*fSkin[3]-0.6629126073623879*coeff[0]*fSkin[3]-0.7411588266019635*coeff[2]*fEdge[3]+0.6629126073623879*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = (-0.2487763020014169*coeff[2]*fSkin[11])-0.5188111786213743*coeff[1]*fSkin[11]+1.121119609893386*coeff[0]*fSkin[11]+1.588341005085966*coeff[2]*fEdge[11]-0.5188111786213743*coeff[1]*fEdge[11]-0.5220480626221116*coeff[0]*fEdge[11]+0.6670429439417671*coeff[2]*fSkin[4]+0.7463289060042488*coeff[1]*fSkin[4]-1.889300930982805*coeff[0]*fSkin[4]+2.594055893106872*coeff[2]*fEdge[4]-0.7463289060042488*coeff[1]*fEdge[4]-1.027514541411701*coeff[0]*fEdge[4]-0.5990715472712751*coeff[2]*fSkin[2]-0.5303300858899105*coeff[1]*fSkin[2]+1.454384534777511*coeff[0]*fSkin[2]+1.96837794103419*coeff[2]*fEdge[2]-0.5303300858899105*coeff[1]*fEdge[2]-0.8420120990817168*coeff[0]*fEdge[2]; 
  edgeSurf_incr[5] = (-0.2487763020014169*coeff[2]*fSkin[13])-0.5188111786213743*coeff[1]*fSkin[13]+1.121119609893386*coeff[0]*fSkin[13]+1.588341005085966*coeff[2]*fEdge[13]-0.5188111786213743*coeff[1]*fEdge[13]-0.5220480626221116*coeff[0]*fEdge[13]+0.6670429439417671*coeff[2]*fSkin[5]+0.7463289060042488*coeff[1]*fSkin[5]-1.889300930982805*coeff[0]*fSkin[5]+2.594055893106872*coeff[2]*fEdge[5]-0.7463289060042488*coeff[1]*fEdge[5]-1.027514541411701*coeff[0]*fEdge[5]-0.5990715472712751*coeff[2]*fSkin[3]-0.5303300858899105*coeff[1]*fSkin[3]+1.454384534777511*coeff[0]*fSkin[3]+1.96837794103419*coeff[2]*fEdge[3]-0.5303300858899105*coeff[1]*fEdge[3]-0.8420120990817168*coeff[0]*fEdge[3]; 
  edgeSurf_incr[6] = 0.5303300858899105*coeff[2]*fSkin[17]-0.4743416490252568*coeff[0]*fSkin[17]-0.5303300858899105*coeff[2]*fEdge[17]+0.4743416490252568*coeff[0]*fEdge[17]-0.9413981457120035*coeff[2]*fSkin[10]+0.8420120990817169*coeff[0]*fSkin[10]-0.9413981457120035*coeff[2]*fEdge[10]+0.8420120990817169*coeff[0]*fEdge[10]+0.7411588266019635*coeff[2]*fSkin[6]-0.6629126073623879*coeff[0]*fSkin[6]-0.7411588266019635*coeff[2]*fEdge[6]+0.6629126073623879*coeff[0]*fEdge[6]; 
  edgeSurf_incr[7] = (-1.40820177054373*coeff[2]*fSkin[7])+2.009347054626824*coeff[1]*fSkin[7]-2.220757234663999*coeff[0]*fSkin[7]-3.779910015670014*coeff[2]*fEdge[7]+2.009347054626824*coeff[1]*fEdge[7]-0.0994368911043575*coeff[0]*fEdge[7]+1.626614282316952*fSkin[1]*coeff[2]-5.836674777725536*fEdge[1]*coeff[2]-0.9943689110435825*fSkin[0]*coeff[2]-4.308931947855521*fEdge[0]*coeff[2]-2.890519423747657*coeff[1]*fSkin[1]+3.551638458822559*coeff[0]*fSkin[1]+2.890519423747657*coeff[1]*fEdge[1]+0.2139541240254559*coeff[0]*fEdge[1]+2.053959590644372*fSkin[0]*coeff[1]+2.053959590644372*fEdge[0]*coeff[1]-2.668171775767069*coeff[0]*fSkin[0]+0.2964635306407852*coeff[0]*fEdge[0]; 
  edgeSurf_incr[8] = (-0.9413981457120036*coeff[2]*fSkin[12])+0.842012099081717*coeff[0]*fSkin[12]-0.9413981457120036*coeff[2]*fEdge[12]+0.842012099081717*coeff[0]*fEdge[12]+0.7411588266019635*coeff[2]*fSkin[8]-0.6629126073623879*coeff[0]*fSkin[8]-0.7411588266019635*coeff[2]*fEdge[8]+0.6629126073623879*coeff[0]*fEdge[8]; 
  edgeSurf_incr[9] = (-0.9413981457120036*coeff[2]*fSkin[15])+0.842012099081717*coeff[0]*fSkin[15]-0.9413981457120036*coeff[2]*fEdge[15]+0.842012099081717*coeff[0]*fEdge[15]+0.7411588266019635*coeff[2]*fSkin[9]-0.6629126073623879*coeff[0]*fSkin[9]-0.7411588266019635*coeff[2]*fEdge[9]+0.6629126073623879*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = (-0.2487763020014165*coeff[2]*fSkin[17])-0.5188111786213743*coeff[1]*fSkin[17]+1.121119609893386*coeff[0]*fSkin[17]+1.588341005085966*coeff[2]*fEdge[17]-0.5188111786213743*coeff[1]*fEdge[17]-0.5220480626221115*coeff[0]*fEdge[17]+0.6670429439417671*coeff[2]*fSkin[10]+0.7463289060042488*coeff[1]*fSkin[10]-1.889300930982805*coeff[0]*fSkin[10]+2.594055893106872*coeff[2]*fEdge[10]-0.7463289060042488*coeff[1]*fEdge[10]-1.027514541411701*coeff[0]*fEdge[10]-0.5990715472712751*coeff[2]*fSkin[6]-0.5303300858899105*coeff[1]*fSkin[6]+1.454384534777511*coeff[0]*fSkin[6]+1.96837794103419*coeff[2]*fEdge[6]-0.5303300858899105*coeff[1]*fEdge[6]-0.8420120990817168*coeff[0]*fEdge[6]; 
  edgeSurf_incr[11] = (-1.40820177054373*coeff[2]*fSkin[11])+2.009347054626824*coeff[1]*fSkin[11]-2.220757234663999*coeff[0]*fSkin[11]-3.779910015670014*coeff[2]*fEdge[11]+2.009347054626824*coeff[1]*fEdge[11]-0.0994368911043575*coeff[0]*fEdge[11]+1.626614282316953*coeff[2]*fSkin[4]-2.890519423747656*coeff[1]*fSkin[4]+3.551638458822559*coeff[0]*fSkin[4]-5.836674777725538*coeff[2]*fEdge[4]+2.890519423747656*coeff[1]*fEdge[4]+0.2139541240254559*coeff[0]*fEdge[4]-0.9943689110435823*coeff[2]*fSkin[2]+2.053959590644372*coeff[1]*fSkin[2]-2.668171775767069*coeff[0]*fSkin[2]-4.308931947855521*coeff[2]*fEdge[2]+2.053959590644372*coeff[1]*fEdge[2]+0.2964635306407852*coeff[0]*fEdge[2]; 
  edgeSurf_incr[12] = 0.6670429439417671*coeff[2]*fSkin[12]+0.7463289060042488*coeff[1]*fSkin[12]-1.889300930982805*coeff[0]*fSkin[12]+2.594055893106872*coeff[2]*fEdge[12]-0.7463289060042488*coeff[1]*fEdge[12]-1.027514541411701*coeff[0]*fEdge[12]-0.599071547271275*coeff[2]*fSkin[8]-0.5303300858899104*coeff[1]*fSkin[8]+1.454384534777511*coeff[0]*fSkin[8]+1.96837794103419*coeff[2]*fEdge[8]-0.5303300858899104*coeff[1]*fEdge[8]-0.842012099081717*coeff[0]*fEdge[8]; 
  edgeSurf_incr[13] = (-1.40820177054373*coeff[2]*fSkin[13])+2.009347054626824*coeff[1]*fSkin[13]-2.220757234663999*coeff[0]*fSkin[13]-3.779910015670014*coeff[2]*fEdge[13]+2.009347054626824*coeff[1]*fEdge[13]-0.0994368911043575*coeff[0]*fEdge[13]+1.626614282316953*coeff[2]*fSkin[5]-2.890519423747656*coeff[1]*fSkin[5]+3.551638458822559*coeff[0]*fSkin[5]-5.836674777725538*coeff[2]*fEdge[5]+2.890519423747656*coeff[1]*fEdge[5]+0.2139541240254559*coeff[0]*fEdge[5]-0.9943689110435823*coeff[2]*fSkin[3]+2.053959590644372*coeff[1]*fSkin[3]-2.668171775767069*coeff[0]*fSkin[3]-4.308931947855521*coeff[2]*fEdge[3]+2.053959590644372*coeff[1]*fEdge[3]+0.2964635306407852*coeff[0]*fEdge[3]; 
  edgeSurf_incr[14] = (-0.9413981457120036*coeff[2]*fSkin[18])+0.842012099081717*coeff[0]*fSkin[18]-0.9413981457120036*coeff[2]*fEdge[18]+0.842012099081717*coeff[0]*fEdge[18]+0.7411588266019635*coeff[2]*fSkin[14]-0.6629126073623879*coeff[0]*fSkin[14]-0.7411588266019635*coeff[2]*fEdge[14]+0.6629126073623879*coeff[0]*fEdge[14]; 
  edgeSurf_incr[15] = 0.6670429439417671*coeff[2]*fSkin[15]+0.7463289060042488*coeff[1]*fSkin[15]-1.889300930982805*coeff[0]*fSkin[15]+2.594055893106872*coeff[2]*fEdge[15]-0.7463289060042488*coeff[1]*fEdge[15]-1.027514541411701*coeff[0]*fEdge[15]-0.599071547271275*coeff[2]*fSkin[9]-0.5303300858899104*coeff[1]*fSkin[9]+1.454384534777511*coeff[0]*fSkin[9]+1.96837794103419*coeff[2]*fEdge[9]-0.5303300858899104*coeff[1]*fEdge[9]-0.842012099081717*coeff[0]*fEdge[9]; 
  edgeSurf_incr[16] = (-0.9413981457120036*coeff[2]*fSkin[19])+0.842012099081717*coeff[0]*fSkin[19]-0.9413981457120036*coeff[2]*fEdge[19]+0.842012099081717*coeff[0]*fEdge[19]+0.7411588266019635*coeff[2]*fSkin[16]-0.6629126073623879*coeff[0]*fSkin[16]-0.7411588266019635*coeff[2]*fEdge[16]+0.6629126073623879*coeff[0]*fEdge[16]; 
  edgeSurf_incr[17] = (-1.40820177054373*coeff[2]*fSkin[17])+2.009347054626824*coeff[1]*fSkin[17]-2.220757234663999*coeff[0]*fSkin[17]-3.779910015670014*coeff[2]*fEdge[17]+2.009347054626824*coeff[1]*fEdge[17]-0.0994368911043575*coeff[0]*fEdge[17]+1.626614282316952*coeff[2]*fSkin[10]-2.890519423747657*coeff[1]*fSkin[10]+3.551638458822559*coeff[0]*fSkin[10]-5.836674777725536*coeff[2]*fEdge[10]+2.890519423747657*coeff[1]*fEdge[10]+0.2139541240254559*coeff[0]*fEdge[10]-0.9943689110435825*coeff[2]*fSkin[6]+2.053959590644372*coeff[1]*fSkin[6]-2.668171775767069*coeff[0]*fSkin[6]-4.308931947855521*coeff[2]*fEdge[6]+2.053959590644372*coeff[1]*fEdge[6]+0.2964635306407852*coeff[0]*fEdge[6]; 
  edgeSurf_incr[18] = 0.6670429439417671*coeff[2]*fSkin[18]+0.7463289060042488*coeff[1]*fSkin[18]-1.889300930982805*coeff[0]*fSkin[18]+2.594055893106872*coeff[2]*fEdge[18]-0.7463289060042488*coeff[1]*fEdge[18]-1.027514541411701*coeff[0]*fEdge[18]-0.599071547271275*coeff[2]*fSkin[14]-0.5303300858899104*coeff[1]*fSkin[14]+1.454384534777511*coeff[0]*fSkin[14]+1.96837794103419*coeff[2]*fEdge[14]-0.5303300858899104*coeff[1]*fEdge[14]-0.842012099081717*coeff[0]*fEdge[14]; 
  edgeSurf_incr[19] = 0.6670429439417671*coeff[2]*fSkin[19]+0.7463289060042488*coeff[1]*fSkin[19]-1.889300930982805*coeff[0]*fSkin[19]+2.594055893106872*coeff[2]*fEdge[19]-0.7463289060042488*coeff[1]*fEdge[19]-1.027514541411701*coeff[0]*fEdge[19]-0.599071547271275*coeff[2]*fSkin[16]-0.5303300858899104*coeff[1]*fSkin[16]+1.454384534777511*coeff[0]*fSkin[16]+1.96837794103419*coeff[2]*fEdge[16]-0.5303300858899104*coeff[1]*fEdge[16]-0.842012099081717*coeff[0]*fEdge[16]; 

  boundSurf_incr[1] = (-1.530931089239486*coeff[2]*fSkin[7])-1.185854122563142*coeff[1]*fSkin[7]-0.6846531968814573*coeff[0]*fSkin[7]-1.976423537605237*fSkin[1]*coeff[2]-1.369306393762915*fSkin[0]*coeff[2]-1.530931089239486*coeff[1]*fSkin[1]-0.883883476483184*coeff[0]*fSkin[1]-1.060660171779821*fSkin[0]*coeff[1]-0.6123724356957944*coeff[0]*fSkin[0]; 
  boundSurf_incr[4] = (-1.530931089239486*coeff[2]*fSkin[11])-1.185854122563142*coeff[1]*fSkin[11]-0.6846531968814574*coeff[0]*fSkin[11]-1.976423537605237*coeff[2]*fSkin[4]-1.530931089239486*coeff[1]*fSkin[4]-0.883883476483184*coeff[0]*fSkin[4]-1.369306393762915*coeff[2]*fSkin[2]-1.060660171779821*coeff[1]*fSkin[2]-0.6123724356957944*coeff[0]*fSkin[2]; 
  boundSurf_incr[5] = (-1.530931089239486*coeff[2]*fSkin[13])-1.185854122563142*coeff[1]*fSkin[13]-0.6846531968814574*coeff[0]*fSkin[13]-1.976423537605237*coeff[2]*fSkin[5]-1.530931089239486*coeff[1]*fSkin[5]-0.883883476483184*coeff[0]*fSkin[5]-1.369306393762915*coeff[2]*fSkin[3]-1.060660171779821*coeff[1]*fSkin[3]-0.6123724356957944*coeff[0]*fSkin[3]; 
  boundSurf_incr[7] = (-5.929270612815711*coeff[2]*fSkin[7])-4.592793267718456*coeff[1]*fSkin[7]-2.651650429449552*coeff[0]*fSkin[7]-7.654655446197428*fSkin[1]*coeff[2]-5.303300858899105*fSkin[0]*coeff[2]-5.929270612815711*coeff[1]*fSkin[1]-3.423265984407287*coeff[0]*fSkin[1]-4.107919181288745*fSkin[0]*coeff[1]-2.371708245126284*coeff[0]*fSkin[0]; 
  boundSurf_incr[10] = (-1.530931089239486*coeff[2]*fSkin[17])-1.185854122563142*coeff[1]*fSkin[17]-0.6846531968814573*coeff[0]*fSkin[17]-1.976423537605237*coeff[2]*fSkin[10]-1.530931089239486*coeff[1]*fSkin[10]-0.883883476483184*coeff[0]*fSkin[10]-1.369306393762915*coeff[2]*fSkin[6]-1.060660171779821*coeff[1]*fSkin[6]-0.6123724356957944*coeff[0]*fSkin[6]; 
  boundSurf_incr[11] = (-5.929270612815711*coeff[2]*fSkin[11])-4.592793267718456*coeff[1]*fSkin[11]-2.651650429449552*coeff[0]*fSkin[11]-7.65465544619743*coeff[2]*fSkin[4]-5.929270612815709*coeff[1]*fSkin[4]-3.423265984407287*coeff[0]*fSkin[4]-5.303300858899106*coeff[2]*fSkin[2]-4.107919181288745*coeff[1]*fSkin[2]-2.371708245126284*coeff[0]*fSkin[2]; 
  boundSurf_incr[12] = (-1.976423537605237*coeff[2]*fSkin[12])-1.530931089239486*coeff[1]*fSkin[12]-0.883883476483184*coeff[0]*fSkin[12]-1.369306393762915*coeff[2]*fSkin[8]-1.060660171779821*coeff[1]*fSkin[8]-0.6123724356957944*coeff[0]*fSkin[8]; 
  boundSurf_incr[13] = (-5.929270612815711*coeff[2]*fSkin[13])-4.592793267718456*coeff[1]*fSkin[13]-2.651650429449552*coeff[0]*fSkin[13]-7.65465544619743*coeff[2]*fSkin[5]-5.929270612815709*coeff[1]*fSkin[5]-3.423265984407287*coeff[0]*fSkin[5]-5.303300858899106*coeff[2]*fSkin[3]-4.107919181288745*coeff[1]*fSkin[3]-2.371708245126284*coeff[0]*fSkin[3]; 
  boundSurf_incr[15] = (-1.976423537605237*coeff[2]*fSkin[15])-1.530931089239486*coeff[1]*fSkin[15]-0.883883476483184*coeff[0]*fSkin[15]-1.369306393762915*coeff[2]*fSkin[9]-1.060660171779821*coeff[1]*fSkin[9]-0.6123724356957944*coeff[0]*fSkin[9]; 
  boundSurf_incr[17] = (-5.929270612815711*coeff[2]*fSkin[17])-4.592793267718456*coeff[1]*fSkin[17]-2.651650429449552*coeff[0]*fSkin[17]-7.654655446197428*coeff[2]*fSkin[10]-5.929270612815711*coeff[1]*fSkin[10]-3.423265984407287*coeff[0]*fSkin[10]-5.303300858899105*coeff[2]*fSkin[6]-4.107919181288745*coeff[1]*fSkin[6]-2.371708245126284*coeff[0]*fSkin[6]; 
  boundSurf_incr[18] = (-1.976423537605237*coeff[2]*fSkin[18])-1.530931089239486*coeff[1]*fSkin[18]-0.883883476483184*coeff[0]*fSkin[18]-1.369306393762915*coeff[2]*fSkin[14]-1.060660171779821*coeff[1]*fSkin[14]-0.6123724356957944*coeff[0]*fSkin[14]; 
  boundSurf_incr[19] = (-1.976423537605237*coeff[2]*fSkin[19])-1.530931089239486*coeff[1]*fSkin[19]-0.883883476483184*coeff[0]*fSkin[19]-1.369306393762915*coeff[2]*fSkin[16]-1.060660171779821*coeff[1]*fSkin[16]-0.6123724356957944*coeff[0]*fSkin[16]; 

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
  out[12] += (vol_incr[12]+edgeSurf_incr[12]+boundSurf_incr[12])*Jfac; 
  out[13] += (vol_incr[13]+edgeSurf_incr[13]+boundSurf_incr[13])*Jfac; 
  out[14] += (vol_incr[14]+edgeSurf_incr[14]+boundSurf_incr[14])*Jfac; 
  out[15] += (vol_incr[15]+edgeSurf_incr[15]+boundSurf_incr[15])*Jfac; 
  out[16] += (vol_incr[16]+edgeSurf_incr[16]+boundSurf_incr[16])*Jfac; 
  out[17] += (vol_incr[17]+edgeSurf_incr[17]+boundSurf_incr[17])*Jfac; 
  out[18] += (vol_incr[18]+edgeSurf_incr[18]+boundSurf_incr[18])*Jfac; 
  out[19] += (vol_incr[19]+edgeSurf_incr[19]+boundSurf_incr[19])*Jfac; 

  }

  return 0.;
}

