#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_pkpm_drag_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *nu, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[2]:      Cell-center coordinates. 
  // dxv[2]:    Cell spacing. 
  // nu:     collisionality. 
  // f:         Input distribution function.
  // out:       Incremented output 
  const double rdvpar = 2.0/dxv[1]; 
  const double dvpar = dxv[1], wvpar = w[1]; 
  double alphaDrag[6]; 
  // Expand rdvpar*(nu*vx) in phase basis.
  alphaDrag[0] = -1.414213562373095*nu[0]*rdvpar*wvpar; 
  alphaDrag[1] = -1.414213562373095*nu[1]*rdvpar*wvpar; 
  alphaDrag[2] = -0.408248290463863*nu[0]*dvpar*rdvpar; 
  alphaDrag[3] = -0.408248290463863*nu[1]*dvpar*rdvpar; 

  out[2] += 0.8660254037844386*(alphaDrag[3]*f[3]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.8660254037844386*(alphaDrag[2]*f[3]+f[2]*alphaDrag[3]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[4] += 1.732050807568877*(alphaDrag[3]*f[5]+alphaDrag[2]*f[4])+1.936491673103709*(alphaDrag[1]*f[3]+f[1]*alphaDrag[3]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[5] += 1.732050807568877*(alphaDrag[2]*f[5]+alphaDrag[3]*f[4])+1.936491673103709*(alphaDrag[0]*f[3]+f[0]*alphaDrag[3]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 

  return fabs(1.25*alphaDrag[0]); 

} 
