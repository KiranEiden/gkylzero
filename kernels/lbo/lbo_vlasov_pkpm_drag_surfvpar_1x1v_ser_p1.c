#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH void lbo_vlasov_pkpm_drag_surfvpar_1x1v_ser_p1(const double *w, const double *dxv, const double *nu, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[2]:         cell-center coordinates. 
  // dxv[2]:       cell spacing. 
  // nu:         collisionalities added (self and cross species collisionalities). 
  // fl/fc/fr:      distribution function in cells 
  // out:           incremented distribution function in cell 

  const double dv1par = 2.0/dxv[1]; 
  const double dvpar = dxv[1], wvpar = w[1]; 

  double alphaDrSurf_l[2] = {0.0}; 
  alphaDrSurf_l[0] = nu[0]*wvpar-0.5*nu[0]*dvpar; 
  alphaDrSurf_l[1] = nu[1]*wvpar-0.5*nu[1]*dvpar; 

  double alphaDrSurf_r[2] = {0.0}; 
  alphaDrSurf_r[0] = nu[0]*wvpar+0.5*nu[0]*dvpar; 
  alphaDrSurf_r[1] = nu[1]*wvpar+0.5*nu[1]*dvpar; 

  double Ghat_r[2]; 
  double Ghat_l[2]; 
  if (wvpar>0) { 

  Ghat_r[0] = 1.118033988749895*(alphaDrSurf_r[1]*fr[5]+alphaDrSurf_r[0]*fr[4])-0.8660254037844386*(alphaDrSurf_r[1]*fr[3]+alphaDrSurf_r[0]*fr[2])+0.5*(alphaDrSurf_r[1]*fr[1]+alphaDrSurf_r[0]*fr[0]); 
  Ghat_r[1] = 1.118033988749895*(alphaDrSurf_r[0]*fr[5]+alphaDrSurf_r[1]*fr[4])-0.8660254037844386*(alphaDrSurf_r[0]*fr[3]+alphaDrSurf_r[1]*fr[2])+0.5*(alphaDrSurf_r[0]*fr[1]+fr[0]*alphaDrSurf_r[1]); 

  Ghat_l[0] = 1.118033988749895*(alphaDrSurf_l[1]*fc[5]+alphaDrSurf_l[0]*fc[4])-0.8660254037844386*(alphaDrSurf_l[1]*fc[3]+alphaDrSurf_l[0]*fc[2])+0.5*(alphaDrSurf_l[1]*fc[1]+alphaDrSurf_l[0]*fc[0]); 
  Ghat_l[1] = 1.118033988749895*(alphaDrSurf_l[0]*fc[5]+alphaDrSurf_l[1]*fc[4])-0.8660254037844386*(alphaDrSurf_l[0]*fc[3]+alphaDrSurf_l[1]*fc[2])+0.5*(alphaDrSurf_l[0]*fc[1]+fc[0]*alphaDrSurf_l[1]); 

  } else { 

  Ghat_r[0] = 0.1666666666666667*(6.708203932499369*alphaDrSurf_r[1]*fc[5]+6.708203932499369*alphaDrSurf_r[0]*fc[4]+5.196152422706631*(alphaDrSurf_r[1]*fc[3]+alphaDrSurf_r[0]*fc[2])+3.0*(alphaDrSurf_r[1]*fc[1]+alphaDrSurf_r[0]*fc[0])); 
  Ghat_r[1] = 0.1666666666666667*(6.708203932499369*alphaDrSurf_r[0]*fc[5]+6.708203932499369*alphaDrSurf_r[1]*fc[4]+5.196152422706631*(alphaDrSurf_r[0]*fc[3]+alphaDrSurf_r[1]*fc[2])+3.0*(alphaDrSurf_r[0]*fc[1]+fc[0]*alphaDrSurf_r[1])); 

  Ghat_l[0] = 0.1666666666666667*(6.708203932499369*alphaDrSurf_l[1]*fl[5]+6.708203932499369*alphaDrSurf_l[0]*fl[4]+5.196152422706631*(alphaDrSurf_l[1]*fl[3]+alphaDrSurf_l[0]*fl[2])+3.0*(alphaDrSurf_l[1]*fl[1]+alphaDrSurf_l[0]*fl[0])); 
  Ghat_l[1] = 0.1666666666666667*(6.708203932499369*alphaDrSurf_l[0]*fl[5]+6.708203932499369*alphaDrSurf_l[1]*fl[4]+5.196152422706631*(alphaDrSurf_l[0]*fl[3]+alphaDrSurf_l[1]*fl[2])+3.0*(alphaDrSurf_l[0]*fl[1]+fl[0]*alphaDrSurf_l[1])); 

  } 
  out[0] += (0.7071067811865475*Ghat_r[0]-0.7071067811865475*Ghat_l[0])*dv1par; 
  out[1] += (0.7071067811865475*Ghat_r[1]-0.7071067811865475*Ghat_l[1])*dv1par; 
  out[2] += 1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv1par; 
  out[3] += 1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv1par; 
  out[4] += (1.58113883008419*Ghat_r[0]-1.58113883008419*Ghat_l[0])*dv1par; 
  out[5] += (1.58113883008419*Ghat_r[1]-1.58113883008419*Ghat_l[1])*dv1par; 
} 
