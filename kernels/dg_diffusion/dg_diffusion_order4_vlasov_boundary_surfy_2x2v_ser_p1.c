#include <gkyl_dg_diffusion_kernels.h>

GKYL_CU_DH double dg_diffusion_order4_vlasov_boundary_surfy_2x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[1],4.);

  double vol_incr[32] = {0.0}; 

  double edgeSurf_incr[32] = {0.0}; 
  double boundSurf_incr[32] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = 1.623797632095822*coeff[1]*fSkin[2]+1.623797632095822*coeff[1]*fEdge[2]+0.9375*fSkin[0]*coeff[1]-0.9375*fEdge[0]*coeff[1]; 
  edgeSurf_incr[1] = 1.623797632095822*coeff[1]*fSkin[5]+1.623797632095822*coeff[1]*fEdge[5]+0.9375*coeff[1]*fSkin[1]-0.9375*coeff[1]*fEdge[1]; 
  edgeSurf_incr[2] = 3.5625*coeff[1]*fSkin[2]+2.0625*coeff[1]*fEdge[2]+1.623797632095822*fSkin[0]*coeff[1]-1.623797632095822*fEdge[0]*coeff[1]; 
  edgeSurf_incr[3] = 1.623797632095822*coeff[1]*fSkin[7]+1.623797632095822*coeff[1]*fEdge[7]+0.9375*coeff[1]*fSkin[3]-0.9375*coeff[1]*fEdge[3]; 
  edgeSurf_incr[4] = 1.623797632095822*coeff[1]*fSkin[9]+1.623797632095822*coeff[1]*fEdge[9]+0.9375*coeff[1]*fSkin[4]-0.9375*coeff[1]*fEdge[4]; 
  edgeSurf_incr[5] = 3.5625*coeff[1]*fSkin[5]+2.0625*coeff[1]*fEdge[5]+1.623797632095822*coeff[1]*fSkin[1]-1.623797632095822*coeff[1]*fEdge[1]; 
  edgeSurf_incr[6] = 1.623797632095822*coeff[1]*fSkin[11]+1.623797632095822*coeff[1]*fEdge[11]+0.9375*coeff[1]*fSkin[6]-0.9375*coeff[1]*fEdge[6]; 
  edgeSurf_incr[7] = 3.5625*coeff[1]*fSkin[7]+2.0625*coeff[1]*fEdge[7]+1.623797632095822*coeff[1]*fSkin[3]-1.623797632095822*coeff[1]*fEdge[3]; 
  edgeSurf_incr[8] = 1.623797632095822*coeff[1]*fSkin[12]+1.623797632095822*coeff[1]*fEdge[12]+0.9375*coeff[1]*fSkin[8]-0.9375*coeff[1]*fEdge[8]; 
  edgeSurf_incr[9] = 3.5625*coeff[1]*fSkin[9]+2.0625*coeff[1]*fEdge[9]+1.623797632095822*coeff[1]*fSkin[4]-1.623797632095822*coeff[1]*fEdge[4]; 
  edgeSurf_incr[10] = 1.623797632095822*coeff[1]*fSkin[14]+1.623797632095822*coeff[1]*fEdge[14]+0.9375*coeff[1]*fSkin[10]-0.9375*coeff[1]*fEdge[10]; 
  edgeSurf_incr[11] = 3.5625*coeff[1]*fSkin[11]+2.0625*coeff[1]*fEdge[11]+1.623797632095822*coeff[1]*fSkin[6]-1.623797632095822*coeff[1]*fEdge[6]; 
  edgeSurf_incr[12] = 3.5625*coeff[1]*fSkin[12]+2.0625*coeff[1]*fEdge[12]+1.623797632095822*coeff[1]*fSkin[8]-1.623797632095822*coeff[1]*fEdge[8]; 
  edgeSurf_incr[13] = 1.623797632095822*coeff[1]*fSkin[15]+1.623797632095822*coeff[1]*fEdge[15]+0.9375*coeff[1]*fSkin[13]-0.9375*coeff[1]*fEdge[13]; 
  edgeSurf_incr[14] = 3.5625*coeff[1]*fSkin[14]+2.0625*coeff[1]*fEdge[14]+1.623797632095822*coeff[1]*fSkin[10]-1.623797632095822*coeff[1]*fEdge[10]; 
  edgeSurf_incr[15] = 3.5625*coeff[1]*fSkin[15]+2.0625*coeff[1]*fEdge[15]+1.623797632095822*coeff[1]*fSkin[13]-1.623797632095822*coeff[1]*fEdge[13]; 
  edgeSurf_incr[16] = 1.623797632095823*coeff[1]*fSkin[18]+1.623797632095823*coeff[1]*fEdge[18]+0.9375*coeff[1]*fSkin[16]-0.9375*coeff[1]*fEdge[16]; 
  edgeSurf_incr[17] = 1.623797632095823*coeff[1]*fSkin[20]+1.623797632095823*coeff[1]*fEdge[20]+0.9375*coeff[1]*fSkin[17]-0.9375*coeff[1]*fEdge[17]; 
  edgeSurf_incr[18] = 3.5625*coeff[1]*fSkin[18]+2.0625*coeff[1]*fEdge[18]+1.623797632095823*coeff[1]*fSkin[16]-1.623797632095823*coeff[1]*fEdge[16]; 
  edgeSurf_incr[19] = 1.623797632095823*coeff[1]*fSkin[22]+1.623797632095823*coeff[1]*fEdge[22]+0.9375*coeff[1]*fSkin[19]-0.9375*coeff[1]*fEdge[19]; 
  edgeSurf_incr[20] = 3.5625*coeff[1]*fSkin[20]+2.0625*coeff[1]*fEdge[20]+1.623797632095823*coeff[1]*fSkin[17]-1.623797632095823*coeff[1]*fEdge[17]; 
  edgeSurf_incr[21] = 1.623797632095823*coeff[1]*fSkin[23]+1.623797632095823*coeff[1]*fEdge[23]+0.9375*coeff[1]*fSkin[21]-0.9375*coeff[1]*fEdge[21]; 
  edgeSurf_incr[22] = 3.5625*coeff[1]*fSkin[22]+2.0625*coeff[1]*fEdge[22]+1.623797632095823*coeff[1]*fSkin[19]-1.623797632095823*coeff[1]*fEdge[19]; 
  edgeSurf_incr[23] = 3.5625*coeff[1]*fSkin[23]+2.0625*coeff[1]*fEdge[23]+1.623797632095823*coeff[1]*fSkin[21]-1.623797632095823*coeff[1]*fEdge[21]; 
  edgeSurf_incr[24] = 1.623797632095823*coeff[1]*fSkin[26]+1.623797632095823*coeff[1]*fEdge[26]+0.9375*coeff[1]*fSkin[24]-0.9375*coeff[1]*fEdge[24]; 
  edgeSurf_incr[25] = 1.623797632095823*coeff[1]*fSkin[28]+1.623797632095823*coeff[1]*fEdge[28]+0.9375*coeff[1]*fSkin[25]-0.9375*coeff[1]*fEdge[25]; 
  edgeSurf_incr[26] = 3.5625*coeff[1]*fSkin[26]+2.0625*coeff[1]*fEdge[26]+1.623797632095823*coeff[1]*fSkin[24]-1.623797632095823*coeff[1]*fEdge[24]; 
  edgeSurf_incr[27] = 1.623797632095823*coeff[1]*fSkin[30]+1.623797632095823*coeff[1]*fEdge[30]+0.9375*coeff[1]*fSkin[27]-0.9375*coeff[1]*fEdge[27]; 
  edgeSurf_incr[28] = 3.5625*coeff[1]*fSkin[28]+2.0625*coeff[1]*fEdge[28]+1.623797632095823*coeff[1]*fSkin[25]-1.623797632095823*coeff[1]*fEdge[25]; 
  edgeSurf_incr[29] = 1.623797632095823*coeff[1]*fSkin[31]+1.623797632095823*coeff[1]*fEdge[31]+0.9375*coeff[1]*fSkin[29]-0.9375*coeff[1]*fEdge[29]; 
  edgeSurf_incr[30] = 3.5625*coeff[1]*fSkin[30]+2.0625*coeff[1]*fEdge[30]+1.623797632095823*coeff[1]*fSkin[27]-1.623797632095823*coeff[1]*fEdge[27]; 
  edgeSurf_incr[31] = 3.5625*coeff[1]*fSkin[31]+2.0625*coeff[1]*fEdge[31]+1.623797632095823*coeff[1]*fSkin[29]-1.623797632095823*coeff[1]*fEdge[29]; 

  boundSurf_incr[2] = 1.5*coeff[1]*fSkin[2]; 
  boundSurf_incr[5] = 1.5*coeff[1]*fSkin[5]; 
  boundSurf_incr[7] = 1.5*coeff[1]*fSkin[7]; 
  boundSurf_incr[9] = 1.5*coeff[1]*fSkin[9]; 
  boundSurf_incr[11] = 1.5*coeff[1]*fSkin[11]; 
  boundSurf_incr[12] = 1.5*coeff[1]*fSkin[12]; 
  boundSurf_incr[14] = 1.5*coeff[1]*fSkin[14]; 
  boundSurf_incr[15] = 1.5*coeff[1]*fSkin[15]; 
  boundSurf_incr[18] = 1.5*coeff[1]*fSkin[18]; 
  boundSurf_incr[20] = 1.5*coeff[1]*fSkin[20]; 
  boundSurf_incr[22] = 1.5*coeff[1]*fSkin[22]; 
  boundSurf_incr[23] = 1.5*coeff[1]*fSkin[23]; 
  boundSurf_incr[26] = 1.5*coeff[1]*fSkin[26]; 
  boundSurf_incr[28] = 1.5*coeff[1]*fSkin[28]; 
  boundSurf_incr[30] = 1.5*coeff[1]*fSkin[30]; 
  boundSurf_incr[31] = 1.5*coeff[1]*fSkin[31]; 

  } else { 

  edgeSurf_incr[0] = (-1.623797632095822*coeff[1]*fSkin[2])-1.623797632095822*coeff[1]*fEdge[2]+0.9375*fSkin[0]*coeff[1]-0.9375*fEdge[0]*coeff[1]; 
  edgeSurf_incr[1] = (-1.623797632095822*coeff[1]*fSkin[5])-1.623797632095822*coeff[1]*fEdge[5]+0.9375*coeff[1]*fSkin[1]-0.9375*coeff[1]*fEdge[1]; 
  edgeSurf_incr[2] = 3.5625*coeff[1]*fSkin[2]+2.0625*coeff[1]*fEdge[2]-1.623797632095822*fSkin[0]*coeff[1]+1.623797632095822*fEdge[0]*coeff[1]; 
  edgeSurf_incr[3] = (-1.623797632095822*coeff[1]*fSkin[7])-1.623797632095822*coeff[1]*fEdge[7]+0.9375*coeff[1]*fSkin[3]-0.9375*coeff[1]*fEdge[3]; 
  edgeSurf_incr[4] = (-1.623797632095822*coeff[1]*fSkin[9])-1.623797632095822*coeff[1]*fEdge[9]+0.9375*coeff[1]*fSkin[4]-0.9375*coeff[1]*fEdge[4]; 
  edgeSurf_incr[5] = 3.5625*coeff[1]*fSkin[5]+2.0625*coeff[1]*fEdge[5]-1.623797632095822*coeff[1]*fSkin[1]+1.623797632095822*coeff[1]*fEdge[1]; 
  edgeSurf_incr[6] = (-1.623797632095822*coeff[1]*fSkin[11])-1.623797632095822*coeff[1]*fEdge[11]+0.9375*coeff[1]*fSkin[6]-0.9375*coeff[1]*fEdge[6]; 
  edgeSurf_incr[7] = 3.5625*coeff[1]*fSkin[7]+2.0625*coeff[1]*fEdge[7]-1.623797632095822*coeff[1]*fSkin[3]+1.623797632095822*coeff[1]*fEdge[3]; 
  edgeSurf_incr[8] = (-1.623797632095822*coeff[1]*fSkin[12])-1.623797632095822*coeff[1]*fEdge[12]+0.9375*coeff[1]*fSkin[8]-0.9375*coeff[1]*fEdge[8]; 
  edgeSurf_incr[9] = 3.5625*coeff[1]*fSkin[9]+2.0625*coeff[1]*fEdge[9]-1.623797632095822*coeff[1]*fSkin[4]+1.623797632095822*coeff[1]*fEdge[4]; 
  edgeSurf_incr[10] = (-1.623797632095822*coeff[1]*fSkin[14])-1.623797632095822*coeff[1]*fEdge[14]+0.9375*coeff[1]*fSkin[10]-0.9375*coeff[1]*fEdge[10]; 
  edgeSurf_incr[11] = 3.5625*coeff[1]*fSkin[11]+2.0625*coeff[1]*fEdge[11]-1.623797632095822*coeff[1]*fSkin[6]+1.623797632095822*coeff[1]*fEdge[6]; 
  edgeSurf_incr[12] = 3.5625*coeff[1]*fSkin[12]+2.0625*coeff[1]*fEdge[12]-1.623797632095822*coeff[1]*fSkin[8]+1.623797632095822*coeff[1]*fEdge[8]; 
  edgeSurf_incr[13] = (-1.623797632095822*coeff[1]*fSkin[15])-1.623797632095822*coeff[1]*fEdge[15]+0.9375*coeff[1]*fSkin[13]-0.9375*coeff[1]*fEdge[13]; 
  edgeSurf_incr[14] = 3.5625*coeff[1]*fSkin[14]+2.0625*coeff[1]*fEdge[14]-1.623797632095822*coeff[1]*fSkin[10]+1.623797632095822*coeff[1]*fEdge[10]; 
  edgeSurf_incr[15] = 3.5625*coeff[1]*fSkin[15]+2.0625*coeff[1]*fEdge[15]-1.623797632095822*coeff[1]*fSkin[13]+1.623797632095822*coeff[1]*fEdge[13]; 
  edgeSurf_incr[16] = (-1.623797632095823*coeff[1]*fSkin[18])-1.623797632095823*coeff[1]*fEdge[18]+0.9375*coeff[1]*fSkin[16]-0.9375*coeff[1]*fEdge[16]; 
  edgeSurf_incr[17] = (-1.623797632095823*coeff[1]*fSkin[20])-1.623797632095823*coeff[1]*fEdge[20]+0.9375*coeff[1]*fSkin[17]-0.9375*coeff[1]*fEdge[17]; 
  edgeSurf_incr[18] = 3.5625*coeff[1]*fSkin[18]+2.0625*coeff[1]*fEdge[18]-1.623797632095823*coeff[1]*fSkin[16]+1.623797632095823*coeff[1]*fEdge[16]; 
  edgeSurf_incr[19] = (-1.623797632095823*coeff[1]*fSkin[22])-1.623797632095823*coeff[1]*fEdge[22]+0.9375*coeff[1]*fSkin[19]-0.9375*coeff[1]*fEdge[19]; 
  edgeSurf_incr[20] = 3.5625*coeff[1]*fSkin[20]+2.0625*coeff[1]*fEdge[20]-1.623797632095823*coeff[1]*fSkin[17]+1.623797632095823*coeff[1]*fEdge[17]; 
  edgeSurf_incr[21] = (-1.623797632095823*coeff[1]*fSkin[23])-1.623797632095823*coeff[1]*fEdge[23]+0.9375*coeff[1]*fSkin[21]-0.9375*coeff[1]*fEdge[21]; 
  edgeSurf_incr[22] = 3.5625*coeff[1]*fSkin[22]+2.0625*coeff[1]*fEdge[22]-1.623797632095823*coeff[1]*fSkin[19]+1.623797632095823*coeff[1]*fEdge[19]; 
  edgeSurf_incr[23] = 3.5625*coeff[1]*fSkin[23]+2.0625*coeff[1]*fEdge[23]-1.623797632095823*coeff[1]*fSkin[21]+1.623797632095823*coeff[1]*fEdge[21]; 
  edgeSurf_incr[24] = (-1.623797632095823*coeff[1]*fSkin[26])-1.623797632095823*coeff[1]*fEdge[26]+0.9375*coeff[1]*fSkin[24]-0.9375*coeff[1]*fEdge[24]; 
  edgeSurf_incr[25] = (-1.623797632095823*coeff[1]*fSkin[28])-1.623797632095823*coeff[1]*fEdge[28]+0.9375*coeff[1]*fSkin[25]-0.9375*coeff[1]*fEdge[25]; 
  edgeSurf_incr[26] = 3.5625*coeff[1]*fSkin[26]+2.0625*coeff[1]*fEdge[26]-1.623797632095823*coeff[1]*fSkin[24]+1.623797632095823*coeff[1]*fEdge[24]; 
  edgeSurf_incr[27] = (-1.623797632095823*coeff[1]*fSkin[30])-1.623797632095823*coeff[1]*fEdge[30]+0.9375*coeff[1]*fSkin[27]-0.9375*coeff[1]*fEdge[27]; 
  edgeSurf_incr[28] = 3.5625*coeff[1]*fSkin[28]+2.0625*coeff[1]*fEdge[28]-1.623797632095823*coeff[1]*fSkin[25]+1.623797632095823*coeff[1]*fEdge[25]; 
  edgeSurf_incr[29] = (-1.623797632095823*coeff[1]*fSkin[31])-1.623797632095823*coeff[1]*fEdge[31]+0.9375*coeff[1]*fSkin[29]-0.9375*coeff[1]*fEdge[29]; 
  edgeSurf_incr[30] = 3.5625*coeff[1]*fSkin[30]+2.0625*coeff[1]*fEdge[30]-1.623797632095823*coeff[1]*fSkin[27]+1.623797632095823*coeff[1]*fEdge[27]; 
  edgeSurf_incr[31] = 3.5625*coeff[1]*fSkin[31]+2.0625*coeff[1]*fEdge[31]-1.623797632095823*coeff[1]*fSkin[29]+1.623797632095823*coeff[1]*fEdge[29]; 

  boundSurf_incr[2] = 1.5*coeff[1]*fSkin[2]; 
  boundSurf_incr[5] = 1.5*coeff[1]*fSkin[5]; 
  boundSurf_incr[7] = 1.5*coeff[1]*fSkin[7]; 
  boundSurf_incr[9] = 1.5*coeff[1]*fSkin[9]; 
  boundSurf_incr[11] = 1.5*coeff[1]*fSkin[11]; 
  boundSurf_incr[12] = 1.5*coeff[1]*fSkin[12]; 
  boundSurf_incr[14] = 1.5*coeff[1]*fSkin[14]; 
  boundSurf_incr[15] = 1.5*coeff[1]*fSkin[15]; 
  boundSurf_incr[18] = 1.5*coeff[1]*fSkin[18]; 
  boundSurf_incr[20] = 1.5*coeff[1]*fSkin[20]; 
  boundSurf_incr[22] = 1.5*coeff[1]*fSkin[22]; 
  boundSurf_incr[23] = 1.5*coeff[1]*fSkin[23]; 
  boundSurf_incr[26] = 1.5*coeff[1]*fSkin[26]; 
  boundSurf_incr[28] = 1.5*coeff[1]*fSkin[28]; 
  boundSurf_incr[30] = 1.5*coeff[1]*fSkin[30]; 
  boundSurf_incr[31] = 1.5*coeff[1]*fSkin[31]; 

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
  out[12] += -1.0*(vol_incr[12]+edgeSurf_incr[12]+boundSurf_incr[12])*Jfac; 
  out[13] += -1.0*(vol_incr[13]+edgeSurf_incr[13]+boundSurf_incr[13])*Jfac; 
  out[14] += -1.0*(vol_incr[14]+edgeSurf_incr[14]+boundSurf_incr[14])*Jfac; 
  out[15] += -1.0*(vol_incr[15]+edgeSurf_incr[15]+boundSurf_incr[15])*Jfac; 
  out[16] += -1.0*(vol_incr[16]+edgeSurf_incr[16]+boundSurf_incr[16])*Jfac; 
  out[17] += -1.0*(vol_incr[17]+edgeSurf_incr[17]+boundSurf_incr[17])*Jfac; 
  out[18] += -1.0*(vol_incr[18]+edgeSurf_incr[18]+boundSurf_incr[18])*Jfac; 
  out[19] += -1.0*(vol_incr[19]+edgeSurf_incr[19]+boundSurf_incr[19])*Jfac; 
  out[20] += -1.0*(vol_incr[20]+edgeSurf_incr[20]+boundSurf_incr[20])*Jfac; 
  out[21] += -1.0*(vol_incr[21]+edgeSurf_incr[21]+boundSurf_incr[21])*Jfac; 
  out[22] += -1.0*(vol_incr[22]+edgeSurf_incr[22]+boundSurf_incr[22])*Jfac; 
  out[23] += -1.0*(vol_incr[23]+edgeSurf_incr[23]+boundSurf_incr[23])*Jfac; 
  out[24] += -1.0*(vol_incr[24]+edgeSurf_incr[24]+boundSurf_incr[24])*Jfac; 
  out[25] += -1.0*(vol_incr[25]+edgeSurf_incr[25]+boundSurf_incr[25])*Jfac; 
  out[26] += -1.0*(vol_incr[26]+edgeSurf_incr[26]+boundSurf_incr[26])*Jfac; 
  out[27] += -1.0*(vol_incr[27]+edgeSurf_incr[27]+boundSurf_incr[27])*Jfac; 
  out[28] += -1.0*(vol_incr[28]+edgeSurf_incr[28]+boundSurf_incr[28])*Jfac; 
  out[29] += -1.0*(vol_incr[29]+edgeSurf_incr[29]+boundSurf_incr[29])*Jfac; 
  out[30] += -1.0*(vol_incr[30]+edgeSurf_incr[30]+boundSurf_incr[30])*Jfac; 
  out[31] += -1.0*(vol_incr[31]+edgeSurf_incr[31]+boundSurf_incr[31])*Jfac; 

  }

  return 0.;
}

GKYL_CU_DH double dg_diffusion_order4_vlasov_boundary_surfy_2x2v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[1],4.);

  double vol_incr[32] = {0.0}; 

  double edgeSurf_incr[32] = {0.0}; 
  double boundSurf_incr[32] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = (-0.375*fSkin[5]*coeff[7])+0.375*fEdge[5]*coeff[7]-0.375*fSkin[2]*coeff[6]+0.375*fEdge[2]*coeff[6]+0.8118988160479111*coeff[5]*fSkin[5]+0.8118988160479111*coeff[5]*fEdge[5]+0.46875*fSkin[1]*coeff[5]-0.46875*fEdge[1]*coeff[5]+0.8118988160479111*fSkin[2]*coeff[4]+0.8118988160479111*fEdge[2]*coeff[4]+0.46875*fSkin[0]*coeff[4]-0.46875*fEdge[0]*coeff[4]; 
  edgeSurf_incr[1] = (-0.375*fSkin[2]*coeff[7])+0.375*fEdge[2]*coeff[7]-0.375*fSkin[5]*coeff[6]+0.375*fEdge[5]*coeff[6]+0.8118988160479111*coeff[4]*fSkin[5]+0.8118988160479111*coeff[4]*fEdge[5]+0.8118988160479111*fSkin[2]*coeff[5]+0.8118988160479111*fEdge[2]*coeff[5]+0.46875*fSkin[0]*coeff[5]-0.46875*fEdge[0]*coeff[5]+0.46875*fSkin[1]*coeff[4]-0.46875*fEdge[1]*coeff[4]; 
  edgeSurf_incr[2] = (-0.6495190528383289*fSkin[5]*coeff[7])+0.6495190528383289*fEdge[5]*coeff[7]-0.6495190528383289*fSkin[2]*coeff[6]+0.6495190528383289*fEdge[2]*coeff[6]+1.78125*coeff[5]*fSkin[5]+1.03125*coeff[5]*fEdge[5]+0.8118988160479111*fSkin[1]*coeff[5]-0.8118988160479111*fEdge[1]*coeff[5]+1.78125*fSkin[2]*coeff[4]+1.03125*fEdge[2]*coeff[4]+0.8118988160479111*fSkin[0]*coeff[4]-0.8118988160479111*fEdge[0]*coeff[4]; 
  edgeSurf_incr[3] = (-0.375*coeff[7]*fSkin[11])+0.8118988160479111*coeff[5]*fSkin[11]+0.375*coeff[7]*fEdge[11]+0.8118988160479111*coeff[5]*fEdge[11]-0.375*coeff[6]*fSkin[7]+0.8118988160479111*coeff[4]*fSkin[7]+0.375*coeff[6]*fEdge[7]+0.8118988160479111*coeff[4]*fEdge[7]+0.46875*coeff[5]*fSkin[6]-0.46875*coeff[5]*fEdge[6]+0.46875*fSkin[3]*coeff[4]-0.46875*fEdge[3]*coeff[4]; 
  edgeSurf_incr[4] = (-0.375*coeff[7]*fSkin[12])+0.8118988160479111*coeff[5]*fSkin[12]+0.375*coeff[7]*fEdge[12]+0.8118988160479111*coeff[5]*fEdge[12]-0.375*coeff[6]*fSkin[9]+0.8118988160479111*coeff[4]*fSkin[9]+0.375*coeff[6]*fEdge[9]+0.8118988160479111*coeff[4]*fEdge[9]+0.46875*coeff[5]*fSkin[8]-0.46875*coeff[5]*fEdge[8]+0.46875*coeff[4]*fSkin[4]-0.46875*coeff[4]*fEdge[4]; 
  edgeSurf_incr[5] = (-0.6495190528383289*fSkin[2]*coeff[7])+0.6495190528383289*fEdge[2]*coeff[7]-0.6495190528383289*fSkin[5]*coeff[6]+0.6495190528383289*fEdge[5]*coeff[6]+1.78125*coeff[4]*fSkin[5]+1.03125*coeff[4]*fEdge[5]+1.78125*fSkin[2]*coeff[5]+1.03125*fEdge[2]*coeff[5]+0.8118988160479111*fSkin[0]*coeff[5]-0.8118988160479111*fEdge[0]*coeff[5]+0.8118988160479111*fSkin[1]*coeff[4]-0.8118988160479111*fEdge[1]*coeff[4]; 
  edgeSurf_incr[6] = (-0.375*coeff[6]*fSkin[11])+0.8118988160479111*coeff[4]*fSkin[11]+0.375*coeff[6]*fEdge[11]+0.8118988160479111*coeff[4]*fEdge[11]-0.375*coeff[7]*fSkin[7]+0.8118988160479111*coeff[5]*fSkin[7]+0.375*coeff[7]*fEdge[7]+0.8118988160479111*coeff[5]*fEdge[7]+0.46875*coeff[4]*fSkin[6]-0.46875*coeff[4]*fEdge[6]+0.46875*fSkin[3]*coeff[5]-0.46875*fEdge[3]*coeff[5]; 
  edgeSurf_incr[7] = (-0.6495190528383289*coeff[7]*fSkin[11])+1.78125*coeff[5]*fSkin[11]+0.6495190528383289*coeff[7]*fEdge[11]+1.03125*coeff[5]*fEdge[11]-0.6495190528383289*coeff[6]*fSkin[7]+1.78125*coeff[4]*fSkin[7]+0.6495190528383289*coeff[6]*fEdge[7]+1.03125*coeff[4]*fEdge[7]+0.8118988160479111*coeff[5]*fSkin[6]-0.8118988160479111*coeff[5]*fEdge[6]+0.8118988160479111*fSkin[3]*coeff[4]-0.8118988160479111*fEdge[3]*coeff[4]; 
  edgeSurf_incr[8] = (-0.375*coeff[6]*fSkin[12])+0.8118988160479111*coeff[4]*fSkin[12]+0.375*coeff[6]*fEdge[12]+0.8118988160479111*coeff[4]*fEdge[12]-0.375*coeff[7]*fSkin[9]+0.8118988160479111*coeff[5]*fSkin[9]+0.375*coeff[7]*fEdge[9]+0.8118988160479111*coeff[5]*fEdge[9]+0.46875*coeff[4]*fSkin[8]-0.46875*coeff[4]*fEdge[8]+0.46875*fSkin[4]*coeff[5]-0.46875*fEdge[4]*coeff[5]; 
  edgeSurf_incr[9] = (-0.6495190528383289*coeff[7]*fSkin[12])+1.78125*coeff[5]*fSkin[12]+0.6495190528383289*coeff[7]*fEdge[12]+1.03125*coeff[5]*fEdge[12]-0.6495190528383289*coeff[6]*fSkin[9]+1.78125*coeff[4]*fSkin[9]+0.6495190528383289*coeff[6]*fEdge[9]+1.03125*coeff[4]*fEdge[9]+0.8118988160479111*coeff[5]*fSkin[8]-0.8118988160479111*coeff[5]*fEdge[8]+0.8118988160479111*coeff[4]*fSkin[4]-0.8118988160479111*coeff[4]*fEdge[4]; 
  edgeSurf_incr[10] = (-0.375*coeff[7]*fSkin[15])+0.8118988160479111*coeff[5]*fSkin[15]+0.375*coeff[7]*fEdge[15]+0.8118988160479111*coeff[5]*fEdge[15]-0.375*coeff[6]*fSkin[14]+0.8118988160479111*coeff[4]*fSkin[14]+0.375*coeff[6]*fEdge[14]+0.8118988160479111*coeff[4]*fEdge[14]+0.46875*coeff[5]*fSkin[13]-0.46875*coeff[5]*fEdge[13]+0.46875*coeff[4]*fSkin[10]-0.46875*coeff[4]*fEdge[10]; 
  edgeSurf_incr[11] = (-0.6495190528383289*coeff[6]*fSkin[11])+1.78125*coeff[4]*fSkin[11]+0.6495190528383289*coeff[6]*fEdge[11]+1.03125*coeff[4]*fEdge[11]-0.6495190528383289*coeff[7]*fSkin[7]+1.78125*coeff[5]*fSkin[7]+0.6495190528383289*coeff[7]*fEdge[7]+1.03125*coeff[5]*fEdge[7]+0.8118988160479111*coeff[4]*fSkin[6]-0.8118988160479111*coeff[4]*fEdge[6]+0.8118988160479111*fSkin[3]*coeff[5]-0.8118988160479111*fEdge[3]*coeff[5]; 
  edgeSurf_incr[12] = (-0.6495190528383289*coeff[6]*fSkin[12])+1.78125*coeff[4]*fSkin[12]+0.6495190528383289*coeff[6]*fEdge[12]+1.03125*coeff[4]*fEdge[12]-0.6495190528383289*coeff[7]*fSkin[9]+1.78125*coeff[5]*fSkin[9]+0.6495190528383289*coeff[7]*fEdge[9]+1.03125*coeff[5]*fEdge[9]+0.8118988160479111*coeff[4]*fSkin[8]-0.8118988160479111*coeff[4]*fEdge[8]+0.8118988160479111*fSkin[4]*coeff[5]-0.8118988160479111*fEdge[4]*coeff[5]; 
  edgeSurf_incr[13] = (-0.375*coeff[6]*fSkin[15])+0.8118988160479111*coeff[4]*fSkin[15]+0.375*coeff[6]*fEdge[15]+0.8118988160479111*coeff[4]*fEdge[15]-0.375*coeff[7]*fSkin[14]+0.8118988160479111*coeff[5]*fSkin[14]+0.375*coeff[7]*fEdge[14]+0.8118988160479111*coeff[5]*fEdge[14]+0.46875*coeff[4]*fSkin[13]-0.46875*coeff[4]*fEdge[13]+0.46875*coeff[5]*fSkin[10]-0.46875*coeff[5]*fEdge[10]; 
  edgeSurf_incr[14] = (-0.6495190528383289*coeff[7]*fSkin[15])+1.78125*coeff[5]*fSkin[15]+0.6495190528383289*coeff[7]*fEdge[15]+1.03125*coeff[5]*fEdge[15]-0.6495190528383289*coeff[6]*fSkin[14]+1.78125*coeff[4]*fSkin[14]+0.6495190528383289*coeff[6]*fEdge[14]+1.03125*coeff[4]*fEdge[14]+0.8118988160479111*coeff[5]*fSkin[13]-0.8118988160479111*coeff[5]*fEdge[13]+0.8118988160479111*coeff[4]*fSkin[10]-0.8118988160479111*coeff[4]*fEdge[10]; 
  edgeSurf_incr[15] = (-0.6495190528383289*coeff[6]*fSkin[15])+1.78125*coeff[4]*fSkin[15]+0.6495190528383289*coeff[6]*fEdge[15]+1.03125*coeff[4]*fEdge[15]-0.6495190528383289*coeff[7]*fSkin[14]+1.78125*coeff[5]*fSkin[14]+0.6495190528383289*coeff[7]*fEdge[14]+1.03125*coeff[5]*fEdge[14]+0.8118988160479111*coeff[4]*fSkin[13]-0.8118988160479111*coeff[4]*fEdge[13]+0.8118988160479111*coeff[5]*fSkin[10]-0.8118988160479111*coeff[5]*fEdge[10]; 
  edgeSurf_incr[16] = (-0.375*coeff[7]*fSkin[20])+0.8118988160479111*coeff[5]*fSkin[20]+0.375*coeff[7]*fEdge[20]+0.8118988160479111*coeff[5]*fEdge[20]-0.375*coeff[6]*fSkin[18]+0.8118988160479113*coeff[4]*fSkin[18]+0.375*coeff[6]*fEdge[18]+0.8118988160479113*coeff[4]*fEdge[18]+0.4687500000000001*coeff[5]*fSkin[17]-0.4687500000000001*coeff[5]*fEdge[17]+0.46875*coeff[4]*fSkin[16]-0.46875*coeff[4]*fEdge[16]; 
  edgeSurf_incr[17] = (-0.375*coeff[6]*fSkin[20])+0.8118988160479113*coeff[4]*fSkin[20]+0.375*coeff[6]*fEdge[20]+0.8118988160479113*coeff[4]*fEdge[20]-0.375*coeff[7]*fSkin[18]+0.8118988160479111*coeff[5]*fSkin[18]+0.375*coeff[7]*fEdge[18]+0.8118988160479111*coeff[5]*fEdge[18]+0.46875*coeff[4]*fSkin[17]-0.46875*coeff[4]*fEdge[17]+0.4687500000000001*coeff[5]*fSkin[16]-0.4687500000000001*coeff[5]*fEdge[16]; 
  edgeSurf_incr[18] = (-0.649519052838329*coeff[7]*fSkin[20])+1.78125*coeff[5]*fSkin[20]+0.649519052838329*coeff[7]*fEdge[20]+1.03125*coeff[5]*fEdge[20]-0.6495190528383289*coeff[6]*fSkin[18]+1.78125*coeff[4]*fSkin[18]+0.6495190528383289*coeff[6]*fEdge[18]+1.03125*coeff[4]*fEdge[18]+0.8118988160479111*coeff[5]*fSkin[17]-0.8118988160479111*coeff[5]*fEdge[17]+0.8118988160479113*coeff[4]*fSkin[16]-0.8118988160479113*coeff[4]*fEdge[16]; 
  edgeSurf_incr[19] = (-0.375*coeff[7]*fSkin[23])+0.8118988160479111*coeff[5]*fSkin[23]+0.375*coeff[7]*fEdge[23]+0.8118988160479111*coeff[5]*fEdge[23]-0.375*coeff[6]*fSkin[22]+0.8118988160479113*coeff[4]*fSkin[22]+0.375*coeff[6]*fEdge[22]+0.8118988160479113*coeff[4]*fEdge[22]+0.4687500000000001*coeff[5]*fSkin[21]-0.4687500000000001*coeff[5]*fEdge[21]+0.46875*coeff[4]*fSkin[19]-0.46875*coeff[4]*fEdge[19]; 
  edgeSurf_incr[20] = (-0.6495190528383289*coeff[6]*fSkin[20])+1.78125*coeff[4]*fSkin[20]+0.6495190528383289*coeff[6]*fEdge[20]+1.03125*coeff[4]*fEdge[20]-0.649519052838329*coeff[7]*fSkin[18]+1.78125*coeff[5]*fSkin[18]+0.649519052838329*coeff[7]*fEdge[18]+1.03125*coeff[5]*fEdge[18]+0.8118988160479113*coeff[4]*fSkin[17]-0.8118988160479113*coeff[4]*fEdge[17]+0.8118988160479111*coeff[5]*fSkin[16]-0.8118988160479111*coeff[5]*fEdge[16]; 
  edgeSurf_incr[21] = (-0.375*coeff[6]*fSkin[23])+0.8118988160479113*coeff[4]*fSkin[23]+0.375*coeff[6]*fEdge[23]+0.8118988160479113*coeff[4]*fEdge[23]-0.375*coeff[7]*fSkin[22]+0.8118988160479111*coeff[5]*fSkin[22]+0.375*coeff[7]*fEdge[22]+0.8118988160479111*coeff[5]*fEdge[22]+0.46875*coeff[4]*fSkin[21]-0.46875*coeff[4]*fEdge[21]+0.4687500000000001*coeff[5]*fSkin[19]-0.4687500000000001*coeff[5]*fEdge[19]; 
  edgeSurf_incr[22] = (-0.649519052838329*coeff[7]*fSkin[23])+1.78125*coeff[5]*fSkin[23]+0.649519052838329*coeff[7]*fEdge[23]+1.03125*coeff[5]*fEdge[23]-0.6495190528383289*coeff[6]*fSkin[22]+1.78125*coeff[4]*fSkin[22]+0.6495190528383289*coeff[6]*fEdge[22]+1.03125*coeff[4]*fEdge[22]+0.8118988160479111*coeff[5]*fSkin[21]-0.8118988160479111*coeff[5]*fEdge[21]+0.8118988160479113*coeff[4]*fSkin[19]-0.8118988160479113*coeff[4]*fEdge[19]; 
  edgeSurf_incr[23] = (-0.6495190528383289*coeff[6]*fSkin[23])+1.78125*coeff[4]*fSkin[23]+0.6495190528383289*coeff[6]*fEdge[23]+1.03125*coeff[4]*fEdge[23]-0.649519052838329*coeff[7]*fSkin[22]+1.78125*coeff[5]*fSkin[22]+0.649519052838329*coeff[7]*fEdge[22]+1.03125*coeff[5]*fEdge[22]+0.8118988160479113*coeff[4]*fSkin[21]-0.8118988160479113*coeff[4]*fEdge[21]+0.8118988160479111*coeff[5]*fSkin[19]-0.8118988160479111*coeff[5]*fEdge[19]; 
  edgeSurf_incr[24] = (-0.375*coeff[7]*fSkin[28])+0.8118988160479111*coeff[5]*fSkin[28]+0.375*coeff[7]*fEdge[28]+0.8118988160479111*coeff[5]*fEdge[28]-0.375*coeff[6]*fSkin[26]+0.8118988160479113*coeff[4]*fSkin[26]+0.375*coeff[6]*fEdge[26]+0.8118988160479113*coeff[4]*fEdge[26]+0.4687500000000001*coeff[5]*fSkin[25]-0.4687500000000001*coeff[5]*fEdge[25]+0.46875*coeff[4]*fSkin[24]-0.46875*coeff[4]*fEdge[24]; 
  edgeSurf_incr[25] = (-0.375*coeff[6]*fSkin[28])+0.8118988160479113*coeff[4]*fSkin[28]+0.375*coeff[6]*fEdge[28]+0.8118988160479113*coeff[4]*fEdge[28]-0.375*coeff[7]*fSkin[26]+0.8118988160479111*coeff[5]*fSkin[26]+0.375*coeff[7]*fEdge[26]+0.8118988160479111*coeff[5]*fEdge[26]+0.46875*coeff[4]*fSkin[25]-0.46875*coeff[4]*fEdge[25]+0.4687500000000001*coeff[5]*fSkin[24]-0.4687500000000001*coeff[5]*fEdge[24]; 
  edgeSurf_incr[26] = (-0.649519052838329*coeff[7]*fSkin[28])+1.78125*coeff[5]*fSkin[28]+0.649519052838329*coeff[7]*fEdge[28]+1.03125*coeff[5]*fEdge[28]-0.6495190528383289*coeff[6]*fSkin[26]+1.78125*coeff[4]*fSkin[26]+0.6495190528383289*coeff[6]*fEdge[26]+1.03125*coeff[4]*fEdge[26]+0.8118988160479111*coeff[5]*fSkin[25]-0.8118988160479111*coeff[5]*fEdge[25]+0.8118988160479113*coeff[4]*fSkin[24]-0.8118988160479113*coeff[4]*fEdge[24]; 
  edgeSurf_incr[27] = (-0.375*coeff[7]*fSkin[31])+0.8118988160479111*coeff[5]*fSkin[31]+0.375*coeff[7]*fEdge[31]+0.8118988160479111*coeff[5]*fEdge[31]-0.375*coeff[6]*fSkin[30]+0.8118988160479113*coeff[4]*fSkin[30]+0.375*coeff[6]*fEdge[30]+0.8118988160479113*coeff[4]*fEdge[30]+0.4687500000000001*coeff[5]*fSkin[29]-0.4687500000000001*coeff[5]*fEdge[29]+0.46875*coeff[4]*fSkin[27]-0.46875*coeff[4]*fEdge[27]; 
  edgeSurf_incr[28] = (-0.6495190528383289*coeff[6]*fSkin[28])+1.78125*coeff[4]*fSkin[28]+0.6495190528383289*coeff[6]*fEdge[28]+1.03125*coeff[4]*fEdge[28]-0.649519052838329*coeff[7]*fSkin[26]+1.78125*coeff[5]*fSkin[26]+0.649519052838329*coeff[7]*fEdge[26]+1.03125*coeff[5]*fEdge[26]+0.8118988160479113*coeff[4]*fSkin[25]-0.8118988160479113*coeff[4]*fEdge[25]+0.8118988160479111*coeff[5]*fSkin[24]-0.8118988160479111*coeff[5]*fEdge[24]; 
  edgeSurf_incr[29] = (-0.375*coeff[6]*fSkin[31])+0.8118988160479113*coeff[4]*fSkin[31]+0.375*coeff[6]*fEdge[31]+0.8118988160479113*coeff[4]*fEdge[31]-0.375*coeff[7]*fSkin[30]+0.8118988160479111*coeff[5]*fSkin[30]+0.375*coeff[7]*fEdge[30]+0.8118988160479111*coeff[5]*fEdge[30]+0.46875*coeff[4]*fSkin[29]-0.46875*coeff[4]*fEdge[29]+0.4687500000000001*coeff[5]*fSkin[27]-0.4687500000000001*coeff[5]*fEdge[27]; 
  edgeSurf_incr[30] = (-0.649519052838329*coeff[7]*fSkin[31])+1.78125*coeff[5]*fSkin[31]+0.649519052838329*coeff[7]*fEdge[31]+1.03125*coeff[5]*fEdge[31]-0.6495190528383289*coeff[6]*fSkin[30]+1.78125*coeff[4]*fSkin[30]+0.6495190528383289*coeff[6]*fEdge[30]+1.03125*coeff[4]*fEdge[30]+0.8118988160479111*coeff[5]*fSkin[29]-0.8118988160479111*coeff[5]*fEdge[29]+0.8118988160479113*coeff[4]*fSkin[27]-0.8118988160479113*coeff[4]*fEdge[27]; 
  edgeSurf_incr[31] = (-0.6495190528383289*coeff[6]*fSkin[31])+1.78125*coeff[4]*fSkin[31]+0.6495190528383289*coeff[6]*fEdge[31]+1.03125*coeff[4]*fEdge[31]-0.649519052838329*coeff[7]*fSkin[30]+1.78125*coeff[5]*fSkin[30]+0.649519052838329*coeff[7]*fEdge[30]+1.03125*coeff[5]*fEdge[30]+0.8118988160479113*coeff[4]*fSkin[29]-0.8118988160479113*coeff[4]*fEdge[29]+0.8118988160479111*coeff[5]*fSkin[27]-0.8118988160479111*coeff[5]*fEdge[27]; 

  boundSurf_incr[0] = (-0.75*fSkin[5]*coeff[7])-0.75*fSkin[2]*coeff[6]; 
  boundSurf_incr[1] = (-0.75*fSkin[2]*coeff[7])-0.75*fSkin[5]*coeff[6]; 
  boundSurf_incr[2] = 1.299038105676658*fSkin[5]*coeff[7]+1.299038105676658*fSkin[2]*coeff[6]+0.75*coeff[5]*fSkin[5]+0.75*fSkin[2]*coeff[4]; 
  boundSurf_incr[3] = (-0.75*coeff[7]*fSkin[11])-0.75*coeff[6]*fSkin[7]; 
  boundSurf_incr[4] = (-0.75*coeff[7]*fSkin[12])-0.75*coeff[6]*fSkin[9]; 
  boundSurf_incr[5] = 1.299038105676658*fSkin[2]*coeff[7]+1.299038105676658*fSkin[5]*coeff[6]+0.75*coeff[4]*fSkin[5]+0.75*fSkin[2]*coeff[5]; 
  boundSurf_incr[6] = (-0.75*coeff[6]*fSkin[11])-0.75*coeff[7]*fSkin[7]; 
  boundSurf_incr[7] = 1.299038105676658*coeff[7]*fSkin[11]+0.75*coeff[5]*fSkin[11]+1.299038105676658*coeff[6]*fSkin[7]+0.75*coeff[4]*fSkin[7]; 
  boundSurf_incr[8] = (-0.75*coeff[6]*fSkin[12])-0.75*coeff[7]*fSkin[9]; 
  boundSurf_incr[9] = 1.299038105676658*coeff[7]*fSkin[12]+0.75*coeff[5]*fSkin[12]+1.299038105676658*coeff[6]*fSkin[9]+0.75*coeff[4]*fSkin[9]; 
  boundSurf_incr[10] = (-0.75*coeff[7]*fSkin[15])-0.75*coeff[6]*fSkin[14]; 
  boundSurf_incr[11] = 1.299038105676658*coeff[6]*fSkin[11]+0.75*coeff[4]*fSkin[11]+1.299038105676658*coeff[7]*fSkin[7]+0.75*coeff[5]*fSkin[7]; 
  boundSurf_incr[12] = 1.299038105676658*coeff[6]*fSkin[12]+0.75*coeff[4]*fSkin[12]+1.299038105676658*coeff[7]*fSkin[9]+0.75*coeff[5]*fSkin[9]; 
  boundSurf_incr[13] = (-0.75*coeff[6]*fSkin[15])-0.75*coeff[7]*fSkin[14]; 
  boundSurf_incr[14] = 1.299038105676658*coeff[7]*fSkin[15]+0.75*coeff[5]*fSkin[15]+1.299038105676658*coeff[6]*fSkin[14]+0.75*coeff[4]*fSkin[14]; 
  boundSurf_incr[15] = 1.299038105676658*coeff[6]*fSkin[15]+0.75*coeff[4]*fSkin[15]+1.299038105676658*coeff[7]*fSkin[14]+0.75*coeff[5]*fSkin[14]; 
  boundSurf_incr[16] = (-0.75*coeff[7]*fSkin[20])-0.75*coeff[6]*fSkin[18]; 
  boundSurf_incr[17] = (-0.75*coeff[6]*fSkin[20])-0.75*coeff[7]*fSkin[18]; 
  boundSurf_incr[18] = 1.299038105676658*coeff[7]*fSkin[20]+0.75*coeff[5]*fSkin[20]+1.299038105676658*coeff[6]*fSkin[18]+0.75*coeff[4]*fSkin[18]; 
  boundSurf_incr[19] = (-0.75*coeff[7]*fSkin[23])-0.75*coeff[6]*fSkin[22]; 
  boundSurf_incr[20] = 1.299038105676658*coeff[6]*fSkin[20]+0.75*coeff[4]*fSkin[20]+1.299038105676658*coeff[7]*fSkin[18]+0.75*coeff[5]*fSkin[18]; 
  boundSurf_incr[21] = (-0.75*coeff[6]*fSkin[23])-0.75*coeff[7]*fSkin[22]; 
  boundSurf_incr[22] = 1.299038105676658*coeff[7]*fSkin[23]+0.75*coeff[5]*fSkin[23]+1.299038105676658*coeff[6]*fSkin[22]+0.75*coeff[4]*fSkin[22]; 
  boundSurf_incr[23] = 1.299038105676658*coeff[6]*fSkin[23]+0.75*coeff[4]*fSkin[23]+1.299038105676658*coeff[7]*fSkin[22]+0.75*coeff[5]*fSkin[22]; 
  boundSurf_incr[24] = (-0.75*coeff[7]*fSkin[28])-0.75*coeff[6]*fSkin[26]; 
  boundSurf_incr[25] = (-0.75*coeff[6]*fSkin[28])-0.75*coeff[7]*fSkin[26]; 
  boundSurf_incr[26] = 1.299038105676658*coeff[7]*fSkin[28]+0.75*coeff[5]*fSkin[28]+1.299038105676658*coeff[6]*fSkin[26]+0.75*coeff[4]*fSkin[26]; 
  boundSurf_incr[27] = (-0.75*coeff[7]*fSkin[31])-0.75*coeff[6]*fSkin[30]; 
  boundSurf_incr[28] = 1.299038105676658*coeff[6]*fSkin[28]+0.75*coeff[4]*fSkin[28]+1.299038105676658*coeff[7]*fSkin[26]+0.75*coeff[5]*fSkin[26]; 
  boundSurf_incr[29] = (-0.75*coeff[6]*fSkin[31])-0.75*coeff[7]*fSkin[30]; 
  boundSurf_incr[30] = 1.299038105676658*coeff[7]*fSkin[31]+0.75*coeff[5]*fSkin[31]+1.299038105676658*coeff[6]*fSkin[30]+0.75*coeff[4]*fSkin[30]; 
  boundSurf_incr[31] = 1.299038105676658*coeff[6]*fSkin[31]+0.75*coeff[4]*fSkin[31]+1.299038105676658*coeff[7]*fSkin[30]+0.75*coeff[5]*fSkin[30]; 

  } else { 

  edgeSurf_incr[0] = (-0.375*fSkin[5]*coeff[7])+0.375*fEdge[5]*coeff[7]-0.375*fSkin[2]*coeff[6]+0.375*fEdge[2]*coeff[6]-0.8118988160479111*coeff[5]*fSkin[5]-0.8118988160479111*coeff[5]*fEdge[5]+0.46875*fSkin[1]*coeff[5]-0.46875*fEdge[1]*coeff[5]-0.8118988160479111*fSkin[2]*coeff[4]-0.8118988160479111*fEdge[2]*coeff[4]+0.46875*fSkin[0]*coeff[4]-0.46875*fEdge[0]*coeff[4]; 
  edgeSurf_incr[1] = (-0.375*fSkin[2]*coeff[7])+0.375*fEdge[2]*coeff[7]-0.375*fSkin[5]*coeff[6]+0.375*fEdge[5]*coeff[6]-0.8118988160479111*coeff[4]*fSkin[5]-0.8118988160479111*coeff[4]*fEdge[5]-0.8118988160479111*fSkin[2]*coeff[5]-0.8118988160479111*fEdge[2]*coeff[5]+0.46875*fSkin[0]*coeff[5]-0.46875*fEdge[0]*coeff[5]+0.46875*fSkin[1]*coeff[4]-0.46875*fEdge[1]*coeff[4]; 
  edgeSurf_incr[2] = 0.6495190528383289*fSkin[5]*coeff[7]-0.6495190528383289*fEdge[5]*coeff[7]+0.6495190528383289*fSkin[2]*coeff[6]-0.6495190528383289*fEdge[2]*coeff[6]+1.78125*coeff[5]*fSkin[5]+1.03125*coeff[5]*fEdge[5]-0.8118988160479111*fSkin[1]*coeff[5]+0.8118988160479111*fEdge[1]*coeff[5]+1.78125*fSkin[2]*coeff[4]+1.03125*fEdge[2]*coeff[4]-0.8118988160479111*fSkin[0]*coeff[4]+0.8118988160479111*fEdge[0]*coeff[4]; 
  edgeSurf_incr[3] = (-0.375*coeff[7]*fSkin[11])-0.8118988160479111*coeff[5]*fSkin[11]+0.375*coeff[7]*fEdge[11]-0.8118988160479111*coeff[5]*fEdge[11]-0.375*coeff[6]*fSkin[7]-0.8118988160479111*coeff[4]*fSkin[7]+0.375*coeff[6]*fEdge[7]-0.8118988160479111*coeff[4]*fEdge[7]+0.46875*coeff[5]*fSkin[6]-0.46875*coeff[5]*fEdge[6]+0.46875*fSkin[3]*coeff[4]-0.46875*fEdge[3]*coeff[4]; 
  edgeSurf_incr[4] = (-0.375*coeff[7]*fSkin[12])-0.8118988160479111*coeff[5]*fSkin[12]+0.375*coeff[7]*fEdge[12]-0.8118988160479111*coeff[5]*fEdge[12]-0.375*coeff[6]*fSkin[9]-0.8118988160479111*coeff[4]*fSkin[9]+0.375*coeff[6]*fEdge[9]-0.8118988160479111*coeff[4]*fEdge[9]+0.46875*coeff[5]*fSkin[8]-0.46875*coeff[5]*fEdge[8]+0.46875*coeff[4]*fSkin[4]-0.46875*coeff[4]*fEdge[4]; 
  edgeSurf_incr[5] = 0.6495190528383289*fSkin[2]*coeff[7]-0.6495190528383289*fEdge[2]*coeff[7]+0.6495190528383289*fSkin[5]*coeff[6]-0.6495190528383289*fEdge[5]*coeff[6]+1.78125*coeff[4]*fSkin[5]+1.03125*coeff[4]*fEdge[5]+1.78125*fSkin[2]*coeff[5]+1.03125*fEdge[2]*coeff[5]-0.8118988160479111*fSkin[0]*coeff[5]+0.8118988160479111*fEdge[0]*coeff[5]-0.8118988160479111*fSkin[1]*coeff[4]+0.8118988160479111*fEdge[1]*coeff[4]; 
  edgeSurf_incr[6] = (-0.375*coeff[6]*fSkin[11])-0.8118988160479111*coeff[4]*fSkin[11]+0.375*coeff[6]*fEdge[11]-0.8118988160479111*coeff[4]*fEdge[11]-0.375*coeff[7]*fSkin[7]-0.8118988160479111*coeff[5]*fSkin[7]+0.375*coeff[7]*fEdge[7]-0.8118988160479111*coeff[5]*fEdge[7]+0.46875*coeff[4]*fSkin[6]-0.46875*coeff[4]*fEdge[6]+0.46875*fSkin[3]*coeff[5]-0.46875*fEdge[3]*coeff[5]; 
  edgeSurf_incr[7] = 0.6495190528383289*coeff[7]*fSkin[11]+1.78125*coeff[5]*fSkin[11]-0.6495190528383289*coeff[7]*fEdge[11]+1.03125*coeff[5]*fEdge[11]+0.6495190528383289*coeff[6]*fSkin[7]+1.78125*coeff[4]*fSkin[7]-0.6495190528383289*coeff[6]*fEdge[7]+1.03125*coeff[4]*fEdge[7]-0.8118988160479111*coeff[5]*fSkin[6]+0.8118988160479111*coeff[5]*fEdge[6]-0.8118988160479111*fSkin[3]*coeff[4]+0.8118988160479111*fEdge[3]*coeff[4]; 
  edgeSurf_incr[8] = (-0.375*coeff[6]*fSkin[12])-0.8118988160479111*coeff[4]*fSkin[12]+0.375*coeff[6]*fEdge[12]-0.8118988160479111*coeff[4]*fEdge[12]-0.375*coeff[7]*fSkin[9]-0.8118988160479111*coeff[5]*fSkin[9]+0.375*coeff[7]*fEdge[9]-0.8118988160479111*coeff[5]*fEdge[9]+0.46875*coeff[4]*fSkin[8]-0.46875*coeff[4]*fEdge[8]+0.46875*fSkin[4]*coeff[5]-0.46875*fEdge[4]*coeff[5]; 
  edgeSurf_incr[9] = 0.6495190528383289*coeff[7]*fSkin[12]+1.78125*coeff[5]*fSkin[12]-0.6495190528383289*coeff[7]*fEdge[12]+1.03125*coeff[5]*fEdge[12]+0.6495190528383289*coeff[6]*fSkin[9]+1.78125*coeff[4]*fSkin[9]-0.6495190528383289*coeff[6]*fEdge[9]+1.03125*coeff[4]*fEdge[9]-0.8118988160479111*coeff[5]*fSkin[8]+0.8118988160479111*coeff[5]*fEdge[8]-0.8118988160479111*coeff[4]*fSkin[4]+0.8118988160479111*coeff[4]*fEdge[4]; 
  edgeSurf_incr[10] = (-0.375*coeff[7]*fSkin[15])-0.8118988160479111*coeff[5]*fSkin[15]+0.375*coeff[7]*fEdge[15]-0.8118988160479111*coeff[5]*fEdge[15]-0.375*coeff[6]*fSkin[14]-0.8118988160479111*coeff[4]*fSkin[14]+0.375*coeff[6]*fEdge[14]-0.8118988160479111*coeff[4]*fEdge[14]+0.46875*coeff[5]*fSkin[13]-0.46875*coeff[5]*fEdge[13]+0.46875*coeff[4]*fSkin[10]-0.46875*coeff[4]*fEdge[10]; 
  edgeSurf_incr[11] = 0.6495190528383289*coeff[6]*fSkin[11]+1.78125*coeff[4]*fSkin[11]-0.6495190528383289*coeff[6]*fEdge[11]+1.03125*coeff[4]*fEdge[11]+0.6495190528383289*coeff[7]*fSkin[7]+1.78125*coeff[5]*fSkin[7]-0.6495190528383289*coeff[7]*fEdge[7]+1.03125*coeff[5]*fEdge[7]-0.8118988160479111*coeff[4]*fSkin[6]+0.8118988160479111*coeff[4]*fEdge[6]-0.8118988160479111*fSkin[3]*coeff[5]+0.8118988160479111*fEdge[3]*coeff[5]; 
  edgeSurf_incr[12] = 0.6495190528383289*coeff[6]*fSkin[12]+1.78125*coeff[4]*fSkin[12]-0.6495190528383289*coeff[6]*fEdge[12]+1.03125*coeff[4]*fEdge[12]+0.6495190528383289*coeff[7]*fSkin[9]+1.78125*coeff[5]*fSkin[9]-0.6495190528383289*coeff[7]*fEdge[9]+1.03125*coeff[5]*fEdge[9]-0.8118988160479111*coeff[4]*fSkin[8]+0.8118988160479111*coeff[4]*fEdge[8]-0.8118988160479111*fSkin[4]*coeff[5]+0.8118988160479111*fEdge[4]*coeff[5]; 
  edgeSurf_incr[13] = (-0.375*coeff[6]*fSkin[15])-0.8118988160479111*coeff[4]*fSkin[15]+0.375*coeff[6]*fEdge[15]-0.8118988160479111*coeff[4]*fEdge[15]-0.375*coeff[7]*fSkin[14]-0.8118988160479111*coeff[5]*fSkin[14]+0.375*coeff[7]*fEdge[14]-0.8118988160479111*coeff[5]*fEdge[14]+0.46875*coeff[4]*fSkin[13]-0.46875*coeff[4]*fEdge[13]+0.46875*coeff[5]*fSkin[10]-0.46875*coeff[5]*fEdge[10]; 
  edgeSurf_incr[14] = 0.6495190528383289*coeff[7]*fSkin[15]+1.78125*coeff[5]*fSkin[15]-0.6495190528383289*coeff[7]*fEdge[15]+1.03125*coeff[5]*fEdge[15]+0.6495190528383289*coeff[6]*fSkin[14]+1.78125*coeff[4]*fSkin[14]-0.6495190528383289*coeff[6]*fEdge[14]+1.03125*coeff[4]*fEdge[14]-0.8118988160479111*coeff[5]*fSkin[13]+0.8118988160479111*coeff[5]*fEdge[13]-0.8118988160479111*coeff[4]*fSkin[10]+0.8118988160479111*coeff[4]*fEdge[10]; 
  edgeSurf_incr[15] = 0.6495190528383289*coeff[6]*fSkin[15]+1.78125*coeff[4]*fSkin[15]-0.6495190528383289*coeff[6]*fEdge[15]+1.03125*coeff[4]*fEdge[15]+0.6495190528383289*coeff[7]*fSkin[14]+1.78125*coeff[5]*fSkin[14]-0.6495190528383289*coeff[7]*fEdge[14]+1.03125*coeff[5]*fEdge[14]-0.8118988160479111*coeff[4]*fSkin[13]+0.8118988160479111*coeff[4]*fEdge[13]-0.8118988160479111*coeff[5]*fSkin[10]+0.8118988160479111*coeff[5]*fEdge[10]; 
  edgeSurf_incr[16] = (-0.375*coeff[7]*fSkin[20])-0.8118988160479111*coeff[5]*fSkin[20]+0.375*coeff[7]*fEdge[20]-0.8118988160479111*coeff[5]*fEdge[20]-0.375*coeff[6]*fSkin[18]-0.8118988160479113*coeff[4]*fSkin[18]+0.375*coeff[6]*fEdge[18]-0.8118988160479113*coeff[4]*fEdge[18]+0.4687500000000001*coeff[5]*fSkin[17]-0.4687500000000001*coeff[5]*fEdge[17]+0.46875*coeff[4]*fSkin[16]-0.46875*coeff[4]*fEdge[16]; 
  edgeSurf_incr[17] = (-0.375*coeff[6]*fSkin[20])-0.8118988160479113*coeff[4]*fSkin[20]+0.375*coeff[6]*fEdge[20]-0.8118988160479113*coeff[4]*fEdge[20]-0.375*coeff[7]*fSkin[18]-0.8118988160479111*coeff[5]*fSkin[18]+0.375*coeff[7]*fEdge[18]-0.8118988160479111*coeff[5]*fEdge[18]+0.46875*coeff[4]*fSkin[17]-0.46875*coeff[4]*fEdge[17]+0.4687500000000001*coeff[5]*fSkin[16]-0.4687500000000001*coeff[5]*fEdge[16]; 
  edgeSurf_incr[18] = 0.649519052838329*coeff[7]*fSkin[20]+1.78125*coeff[5]*fSkin[20]-0.649519052838329*coeff[7]*fEdge[20]+1.03125*coeff[5]*fEdge[20]+0.6495190528383289*coeff[6]*fSkin[18]+1.78125*coeff[4]*fSkin[18]-0.6495190528383289*coeff[6]*fEdge[18]+1.03125*coeff[4]*fEdge[18]-0.8118988160479111*coeff[5]*fSkin[17]+0.8118988160479111*coeff[5]*fEdge[17]-0.8118988160479113*coeff[4]*fSkin[16]+0.8118988160479113*coeff[4]*fEdge[16]; 
  edgeSurf_incr[19] = (-0.375*coeff[7]*fSkin[23])-0.8118988160479111*coeff[5]*fSkin[23]+0.375*coeff[7]*fEdge[23]-0.8118988160479111*coeff[5]*fEdge[23]-0.375*coeff[6]*fSkin[22]-0.8118988160479113*coeff[4]*fSkin[22]+0.375*coeff[6]*fEdge[22]-0.8118988160479113*coeff[4]*fEdge[22]+0.4687500000000001*coeff[5]*fSkin[21]-0.4687500000000001*coeff[5]*fEdge[21]+0.46875*coeff[4]*fSkin[19]-0.46875*coeff[4]*fEdge[19]; 
  edgeSurf_incr[20] = 0.6495190528383289*coeff[6]*fSkin[20]+1.78125*coeff[4]*fSkin[20]-0.6495190528383289*coeff[6]*fEdge[20]+1.03125*coeff[4]*fEdge[20]+0.649519052838329*coeff[7]*fSkin[18]+1.78125*coeff[5]*fSkin[18]-0.649519052838329*coeff[7]*fEdge[18]+1.03125*coeff[5]*fEdge[18]-0.8118988160479113*coeff[4]*fSkin[17]+0.8118988160479113*coeff[4]*fEdge[17]-0.8118988160479111*coeff[5]*fSkin[16]+0.8118988160479111*coeff[5]*fEdge[16]; 
  edgeSurf_incr[21] = (-0.375*coeff[6]*fSkin[23])-0.8118988160479113*coeff[4]*fSkin[23]+0.375*coeff[6]*fEdge[23]-0.8118988160479113*coeff[4]*fEdge[23]-0.375*coeff[7]*fSkin[22]-0.8118988160479111*coeff[5]*fSkin[22]+0.375*coeff[7]*fEdge[22]-0.8118988160479111*coeff[5]*fEdge[22]+0.46875*coeff[4]*fSkin[21]-0.46875*coeff[4]*fEdge[21]+0.4687500000000001*coeff[5]*fSkin[19]-0.4687500000000001*coeff[5]*fEdge[19]; 
  edgeSurf_incr[22] = 0.649519052838329*coeff[7]*fSkin[23]+1.78125*coeff[5]*fSkin[23]-0.649519052838329*coeff[7]*fEdge[23]+1.03125*coeff[5]*fEdge[23]+0.6495190528383289*coeff[6]*fSkin[22]+1.78125*coeff[4]*fSkin[22]-0.6495190528383289*coeff[6]*fEdge[22]+1.03125*coeff[4]*fEdge[22]-0.8118988160479111*coeff[5]*fSkin[21]+0.8118988160479111*coeff[5]*fEdge[21]-0.8118988160479113*coeff[4]*fSkin[19]+0.8118988160479113*coeff[4]*fEdge[19]; 
  edgeSurf_incr[23] = 0.6495190528383289*coeff[6]*fSkin[23]+1.78125*coeff[4]*fSkin[23]-0.6495190528383289*coeff[6]*fEdge[23]+1.03125*coeff[4]*fEdge[23]+0.649519052838329*coeff[7]*fSkin[22]+1.78125*coeff[5]*fSkin[22]-0.649519052838329*coeff[7]*fEdge[22]+1.03125*coeff[5]*fEdge[22]-0.8118988160479113*coeff[4]*fSkin[21]+0.8118988160479113*coeff[4]*fEdge[21]-0.8118988160479111*coeff[5]*fSkin[19]+0.8118988160479111*coeff[5]*fEdge[19]; 
  edgeSurf_incr[24] = (-0.375*coeff[7]*fSkin[28])-0.8118988160479111*coeff[5]*fSkin[28]+0.375*coeff[7]*fEdge[28]-0.8118988160479111*coeff[5]*fEdge[28]-0.375*coeff[6]*fSkin[26]-0.8118988160479113*coeff[4]*fSkin[26]+0.375*coeff[6]*fEdge[26]-0.8118988160479113*coeff[4]*fEdge[26]+0.4687500000000001*coeff[5]*fSkin[25]-0.4687500000000001*coeff[5]*fEdge[25]+0.46875*coeff[4]*fSkin[24]-0.46875*coeff[4]*fEdge[24]; 
  edgeSurf_incr[25] = (-0.375*coeff[6]*fSkin[28])-0.8118988160479113*coeff[4]*fSkin[28]+0.375*coeff[6]*fEdge[28]-0.8118988160479113*coeff[4]*fEdge[28]-0.375*coeff[7]*fSkin[26]-0.8118988160479111*coeff[5]*fSkin[26]+0.375*coeff[7]*fEdge[26]-0.8118988160479111*coeff[5]*fEdge[26]+0.46875*coeff[4]*fSkin[25]-0.46875*coeff[4]*fEdge[25]+0.4687500000000001*coeff[5]*fSkin[24]-0.4687500000000001*coeff[5]*fEdge[24]; 
  edgeSurf_incr[26] = 0.649519052838329*coeff[7]*fSkin[28]+1.78125*coeff[5]*fSkin[28]-0.649519052838329*coeff[7]*fEdge[28]+1.03125*coeff[5]*fEdge[28]+0.6495190528383289*coeff[6]*fSkin[26]+1.78125*coeff[4]*fSkin[26]-0.6495190528383289*coeff[6]*fEdge[26]+1.03125*coeff[4]*fEdge[26]-0.8118988160479111*coeff[5]*fSkin[25]+0.8118988160479111*coeff[5]*fEdge[25]-0.8118988160479113*coeff[4]*fSkin[24]+0.8118988160479113*coeff[4]*fEdge[24]; 
  edgeSurf_incr[27] = (-0.375*coeff[7]*fSkin[31])-0.8118988160479111*coeff[5]*fSkin[31]+0.375*coeff[7]*fEdge[31]-0.8118988160479111*coeff[5]*fEdge[31]-0.375*coeff[6]*fSkin[30]-0.8118988160479113*coeff[4]*fSkin[30]+0.375*coeff[6]*fEdge[30]-0.8118988160479113*coeff[4]*fEdge[30]+0.4687500000000001*coeff[5]*fSkin[29]-0.4687500000000001*coeff[5]*fEdge[29]+0.46875*coeff[4]*fSkin[27]-0.46875*coeff[4]*fEdge[27]; 
  edgeSurf_incr[28] = 0.6495190528383289*coeff[6]*fSkin[28]+1.78125*coeff[4]*fSkin[28]-0.6495190528383289*coeff[6]*fEdge[28]+1.03125*coeff[4]*fEdge[28]+0.649519052838329*coeff[7]*fSkin[26]+1.78125*coeff[5]*fSkin[26]-0.649519052838329*coeff[7]*fEdge[26]+1.03125*coeff[5]*fEdge[26]-0.8118988160479113*coeff[4]*fSkin[25]+0.8118988160479113*coeff[4]*fEdge[25]-0.8118988160479111*coeff[5]*fSkin[24]+0.8118988160479111*coeff[5]*fEdge[24]; 
  edgeSurf_incr[29] = (-0.375*coeff[6]*fSkin[31])-0.8118988160479113*coeff[4]*fSkin[31]+0.375*coeff[6]*fEdge[31]-0.8118988160479113*coeff[4]*fEdge[31]-0.375*coeff[7]*fSkin[30]-0.8118988160479111*coeff[5]*fSkin[30]+0.375*coeff[7]*fEdge[30]-0.8118988160479111*coeff[5]*fEdge[30]+0.46875*coeff[4]*fSkin[29]-0.46875*coeff[4]*fEdge[29]+0.4687500000000001*coeff[5]*fSkin[27]-0.4687500000000001*coeff[5]*fEdge[27]; 
  edgeSurf_incr[30] = 0.649519052838329*coeff[7]*fSkin[31]+1.78125*coeff[5]*fSkin[31]-0.649519052838329*coeff[7]*fEdge[31]+1.03125*coeff[5]*fEdge[31]+0.6495190528383289*coeff[6]*fSkin[30]+1.78125*coeff[4]*fSkin[30]-0.6495190528383289*coeff[6]*fEdge[30]+1.03125*coeff[4]*fEdge[30]-0.8118988160479111*coeff[5]*fSkin[29]+0.8118988160479111*coeff[5]*fEdge[29]-0.8118988160479113*coeff[4]*fSkin[27]+0.8118988160479113*coeff[4]*fEdge[27]; 
  edgeSurf_incr[31] = 0.6495190528383289*coeff[6]*fSkin[31]+1.78125*coeff[4]*fSkin[31]-0.6495190528383289*coeff[6]*fEdge[31]+1.03125*coeff[4]*fEdge[31]+0.649519052838329*coeff[7]*fSkin[30]+1.78125*coeff[5]*fSkin[30]-0.649519052838329*coeff[7]*fEdge[30]+1.03125*coeff[5]*fEdge[30]-0.8118988160479113*coeff[4]*fSkin[29]+0.8118988160479113*coeff[4]*fEdge[29]-0.8118988160479111*coeff[5]*fSkin[27]+0.8118988160479111*coeff[5]*fEdge[27]; 

  boundSurf_incr[0] = (-0.75*fSkin[5]*coeff[7])-0.75*fSkin[2]*coeff[6]; 
  boundSurf_incr[1] = (-0.75*fSkin[2]*coeff[7])-0.75*fSkin[5]*coeff[6]; 
  boundSurf_incr[2] = (-1.299038105676658*fSkin[5]*coeff[7])-1.299038105676658*fSkin[2]*coeff[6]+0.75*coeff[5]*fSkin[5]+0.75*fSkin[2]*coeff[4]; 
  boundSurf_incr[3] = (-0.75*coeff[7]*fSkin[11])-0.75*coeff[6]*fSkin[7]; 
  boundSurf_incr[4] = (-0.75*coeff[7]*fSkin[12])-0.75*coeff[6]*fSkin[9]; 
  boundSurf_incr[5] = (-1.299038105676658*fSkin[2]*coeff[7])-1.299038105676658*fSkin[5]*coeff[6]+0.75*coeff[4]*fSkin[5]+0.75*fSkin[2]*coeff[5]; 
  boundSurf_incr[6] = (-0.75*coeff[6]*fSkin[11])-0.75*coeff[7]*fSkin[7]; 
  boundSurf_incr[7] = (-1.299038105676658*coeff[7]*fSkin[11])+0.75*coeff[5]*fSkin[11]-1.299038105676658*coeff[6]*fSkin[7]+0.75*coeff[4]*fSkin[7]; 
  boundSurf_incr[8] = (-0.75*coeff[6]*fSkin[12])-0.75*coeff[7]*fSkin[9]; 
  boundSurf_incr[9] = (-1.299038105676658*coeff[7]*fSkin[12])+0.75*coeff[5]*fSkin[12]-1.299038105676658*coeff[6]*fSkin[9]+0.75*coeff[4]*fSkin[9]; 
  boundSurf_incr[10] = (-0.75*coeff[7]*fSkin[15])-0.75*coeff[6]*fSkin[14]; 
  boundSurf_incr[11] = (-1.299038105676658*coeff[6]*fSkin[11])+0.75*coeff[4]*fSkin[11]-1.299038105676658*coeff[7]*fSkin[7]+0.75*coeff[5]*fSkin[7]; 
  boundSurf_incr[12] = (-1.299038105676658*coeff[6]*fSkin[12])+0.75*coeff[4]*fSkin[12]-1.299038105676658*coeff[7]*fSkin[9]+0.75*coeff[5]*fSkin[9]; 
  boundSurf_incr[13] = (-0.75*coeff[6]*fSkin[15])-0.75*coeff[7]*fSkin[14]; 
  boundSurf_incr[14] = (-1.299038105676658*coeff[7]*fSkin[15])+0.75*coeff[5]*fSkin[15]-1.299038105676658*coeff[6]*fSkin[14]+0.75*coeff[4]*fSkin[14]; 
  boundSurf_incr[15] = (-1.299038105676658*coeff[6]*fSkin[15])+0.75*coeff[4]*fSkin[15]-1.299038105676658*coeff[7]*fSkin[14]+0.75*coeff[5]*fSkin[14]; 
  boundSurf_incr[16] = (-0.75*coeff[7]*fSkin[20])-0.75*coeff[6]*fSkin[18]; 
  boundSurf_incr[17] = (-0.75*coeff[6]*fSkin[20])-0.75*coeff[7]*fSkin[18]; 
  boundSurf_incr[18] = (-1.299038105676658*coeff[7]*fSkin[20])+0.75*coeff[5]*fSkin[20]-1.299038105676658*coeff[6]*fSkin[18]+0.75*coeff[4]*fSkin[18]; 
  boundSurf_incr[19] = (-0.75*coeff[7]*fSkin[23])-0.75*coeff[6]*fSkin[22]; 
  boundSurf_incr[20] = (-1.299038105676658*coeff[6]*fSkin[20])+0.75*coeff[4]*fSkin[20]-1.299038105676658*coeff[7]*fSkin[18]+0.75*coeff[5]*fSkin[18]; 
  boundSurf_incr[21] = (-0.75*coeff[6]*fSkin[23])-0.75*coeff[7]*fSkin[22]; 
  boundSurf_incr[22] = (-1.299038105676658*coeff[7]*fSkin[23])+0.75*coeff[5]*fSkin[23]-1.299038105676658*coeff[6]*fSkin[22]+0.75*coeff[4]*fSkin[22]; 
  boundSurf_incr[23] = (-1.299038105676658*coeff[6]*fSkin[23])+0.75*coeff[4]*fSkin[23]-1.299038105676658*coeff[7]*fSkin[22]+0.75*coeff[5]*fSkin[22]; 
  boundSurf_incr[24] = (-0.75*coeff[7]*fSkin[28])-0.75*coeff[6]*fSkin[26]; 
  boundSurf_incr[25] = (-0.75*coeff[6]*fSkin[28])-0.75*coeff[7]*fSkin[26]; 
  boundSurf_incr[26] = (-1.299038105676658*coeff[7]*fSkin[28])+0.75*coeff[5]*fSkin[28]-1.299038105676658*coeff[6]*fSkin[26]+0.75*coeff[4]*fSkin[26]; 
  boundSurf_incr[27] = (-0.75*coeff[7]*fSkin[31])-0.75*coeff[6]*fSkin[30]; 
  boundSurf_incr[28] = (-1.299038105676658*coeff[6]*fSkin[28])+0.75*coeff[4]*fSkin[28]-1.299038105676658*coeff[7]*fSkin[26]+0.75*coeff[5]*fSkin[26]; 
  boundSurf_incr[29] = (-0.75*coeff[6]*fSkin[31])-0.75*coeff[7]*fSkin[30]; 
  boundSurf_incr[30] = (-1.299038105676658*coeff[7]*fSkin[31])+0.75*coeff[5]*fSkin[31]-1.299038105676658*coeff[6]*fSkin[30]+0.75*coeff[4]*fSkin[30]; 
  boundSurf_incr[31] = (-1.299038105676658*coeff[6]*fSkin[31])+0.75*coeff[4]*fSkin[31]-1.299038105676658*coeff[7]*fSkin[30]+0.75*coeff[5]*fSkin[30]; 

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
  out[12] += -1.0*(vol_incr[12]+edgeSurf_incr[12]+boundSurf_incr[12])*Jfac; 
  out[13] += -1.0*(vol_incr[13]+edgeSurf_incr[13]+boundSurf_incr[13])*Jfac; 
  out[14] += -1.0*(vol_incr[14]+edgeSurf_incr[14]+boundSurf_incr[14])*Jfac; 
  out[15] += -1.0*(vol_incr[15]+edgeSurf_incr[15]+boundSurf_incr[15])*Jfac; 
  out[16] += -1.0*(vol_incr[16]+edgeSurf_incr[16]+boundSurf_incr[16])*Jfac; 
  out[17] += -1.0*(vol_incr[17]+edgeSurf_incr[17]+boundSurf_incr[17])*Jfac; 
  out[18] += -1.0*(vol_incr[18]+edgeSurf_incr[18]+boundSurf_incr[18])*Jfac; 
  out[19] += -1.0*(vol_incr[19]+edgeSurf_incr[19]+boundSurf_incr[19])*Jfac; 
  out[20] += -1.0*(vol_incr[20]+edgeSurf_incr[20]+boundSurf_incr[20])*Jfac; 
  out[21] += -1.0*(vol_incr[21]+edgeSurf_incr[21]+boundSurf_incr[21])*Jfac; 
  out[22] += -1.0*(vol_incr[22]+edgeSurf_incr[22]+boundSurf_incr[22])*Jfac; 
  out[23] += -1.0*(vol_incr[23]+edgeSurf_incr[23]+boundSurf_incr[23])*Jfac; 
  out[24] += -1.0*(vol_incr[24]+edgeSurf_incr[24]+boundSurf_incr[24])*Jfac; 
  out[25] += -1.0*(vol_incr[25]+edgeSurf_incr[25]+boundSurf_incr[25])*Jfac; 
  out[26] += -1.0*(vol_incr[26]+edgeSurf_incr[26]+boundSurf_incr[26])*Jfac; 
  out[27] += -1.0*(vol_incr[27]+edgeSurf_incr[27]+boundSurf_incr[27])*Jfac; 
  out[28] += -1.0*(vol_incr[28]+edgeSurf_incr[28]+boundSurf_incr[28])*Jfac; 
  out[29] += -1.0*(vol_incr[29]+edgeSurf_incr[29]+boundSurf_incr[29])*Jfac; 
  out[30] += -1.0*(vol_incr[30]+edgeSurf_incr[30]+boundSurf_incr[30])*Jfac; 
  out[31] += -1.0*(vol_incr[31]+edgeSurf_incr[31]+boundSurf_incr[31])*Jfac; 

  }

  return 0.;
}

