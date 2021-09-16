#include <gkyl_vlasov_kernels.h> 

GKYL_CU_DH double vlasov_poisson_vol_1x3v_ser_p1(const double *w, const double *dxv, const double *fac_phi, const double *vecA, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fac_phi:   potential (scaled by appropriate factors).
  // vecA:      vector potential (scaled by appropriate factors). Unused in pure Vlasov-Poisson. 
  // f:         Input distribution function.
  // out:       Incremented output.
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  const double dx10 = 2/dxv[0]; 
  const double *phi = &fac_phi[0]; 
  const double dv10 = 2/dxv[1]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv11 = 2/dxv[2]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double dv12 = 2/dxv[3]; 
  const double dv3 = dxv[3], wv3 = w[3]; 

  double alpha_mid = 0.0; 
  double alpha_cdim[16]; 
  double alpha_vdim[48]; 

  alpha_cdim[0] = 8.0*w0dx0; 
  alpha_cdim[2] = 2.309401076758503*dv0dx0; 
  alpha_mid += fabs(w0dx0)+0.5*dv0dx0; 

  alpha_vdim[0] = -4.898979485566357*phi[1]*dv10*dx10; 
  alpha_mid += fabs(0.125*alpha_vdim[0]); 

  alpha_mid += fabs(0.0); 

  alpha_mid += fabs(0.0); 

  out[1] += 0.4330127018922193*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.4330127018922193*alpha_vdim[0]*f[0]; 
  out[5] += 0.4330127018922193*(alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]); 
  out[6] += 0.4330127018922193*(alpha_cdim[2]*f[7]+alpha_cdim[0]*f[3]); 
  out[7] += 0.4330127018922193*alpha_vdim[0]*f[3]; 
  out[8] += 0.4330127018922193*(alpha_cdim[2]*f[9]+alpha_cdim[0]*f[4]); 
  out[9] += 0.4330127018922193*alpha_vdim[0]*f[4]; 
  out[11] += 0.4330127018922193*(alpha_cdim[0]*f[7]+alpha_vdim[0]*f[6]+alpha_cdim[2]*f[3]); 
  out[12] += 0.4330127018922193*(alpha_cdim[0]*f[9]+alpha_vdim[0]*f[8]+alpha_cdim[2]*f[4]); 
  out[13] += 0.4330127018922193*(alpha_cdim[2]*f[14]+alpha_cdim[0]*f[10]); 
  out[14] += 0.4330127018922193*alpha_vdim[0]*f[10]; 
  out[15] += 0.4330127018922193*(alpha_cdim[0]*f[14]+alpha_vdim[0]*f[13]+alpha_cdim[2]*f[10]); 

  return alpha_mid; 
} 

