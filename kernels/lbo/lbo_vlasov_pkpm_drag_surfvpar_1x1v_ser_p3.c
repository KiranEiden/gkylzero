#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH void lbo_vlasov_pkpm_drag_surfvpar_1x1v_ser_p3(const double *w, const double *dxv, const double *nu, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[2]:         cell-center coordinates. 
  // dxv[2]:       cell spacing. 
  // nu:         collisionalities added (self and cross species collisionalities). 
  // fl/fc/fr:      distribution function in cells 
  // out:           incremented distribution function in cell 

  const double dv1par = 2.0/dxv[1]; 
  const double dvpar = dxv[1], wvpar = w[1]; 

  double alphaDrSurf_l[4] = {0.0}; 
  alphaDrSurf_l[0] = nu[0]*wvpar-0.5*nu[0]*dvpar; 
  alphaDrSurf_l[1] = nu[1]*wvpar-0.5*nu[1]*dvpar; 
  alphaDrSurf_l[2] = nu[2]*wvpar-0.5*nu[2]*dvpar; 
  alphaDrSurf_l[3] = nu[3]*wvpar-0.5*nu[3]*dvpar; 

  double alphaDrSurf_r[4] = {0.0}; 
  alphaDrSurf_r[0] = nu[0]*wvpar+0.5*nu[0]*dvpar; 
  alphaDrSurf_r[1] = nu[1]*wvpar+0.5*nu[1]*dvpar; 
  alphaDrSurf_r[2] = nu[2]*wvpar+0.5*nu[2]*dvpar; 
  alphaDrSurf_r[3] = nu[3]*wvpar+0.5*nu[3]*dvpar; 

  double Ghat_r[4]; 
  double Ghat_l[4]; 
  if (wvpar>0) { 

  Ghat_r[0] = (-1.322875655532295*alphaDrSurf_r[1]*fr[11])-0.8660254037844386*alphaDrSurf_r[3]*fr[10]-1.322875655532295*alphaDrSurf_r[0]*fr[9]+0.5*alphaDrSurf_r[3]*fr[8]+1.118033988749895*alphaDrSurf_r[1]*fr[7]-0.8660254037844386*alphaDrSurf_r[2]*fr[6]+1.118033988749895*alphaDrSurf_r[0]*fr[5]+0.5*alphaDrSurf_r[2]*fr[4]-0.8660254037844386*(alphaDrSurf_r[1]*fr[3]+alphaDrSurf_r[0]*fr[2])+0.5*(alphaDrSurf_r[1]*fr[1]+alphaDrSurf_r[0]*fr[0]); 
  Ghat_r[1] = ((-1.183215956619923*alphaDrSurf_r[2])-1.322875655532295*alphaDrSurf_r[0])*fr[11]-0.7606388292556648*alphaDrSurf_r[2]*fr[10]-1.322875655532295*alphaDrSurf_r[1]*fr[9]+0.4391550328268398*alphaDrSurf_r[2]*fr[8]+(alphaDrSurf_r[2]+1.118033988749895*alphaDrSurf_r[0])*fr[7]-0.7606388292556648*alphaDrSurf_r[3]*fr[6]+alphaDrSurf_r[1]*(1.118033988749895*fr[5]-0.7745966692414833*fr[6])+(0.4391550328268398*alphaDrSurf_r[3]+0.4472135954999579*alphaDrSurf_r[1])*fr[4]-0.7745966692414833*alphaDrSurf_r[2]*fr[3]-0.8660254037844386*(alphaDrSurf_r[0]*fr[3]+alphaDrSurf_r[1]*fr[2])+0.4472135954999579*fr[1]*alphaDrSurf_r[2]+0.5*(alphaDrSurf_r[0]*fr[1]+fr[0]*alphaDrSurf_r[1]); 
  Ghat_r[2] = ((-1.161895003862225*alphaDrSurf_r[3])-1.183215956619923*alphaDrSurf_r[1])*fr[11]+((-0.5163977794943223*alphaDrSurf_r[3])-0.7606388292556648*alphaDrSurf_r[1])*fr[10]-1.322875655532295*alphaDrSurf_r[2]*fr[9]+(0.2981423969999719*alphaDrSurf_r[3]+0.4391550328268398*alphaDrSurf_r[1])*fr[8]+(0.9819805060619655*alphaDrSurf_r[3]+alphaDrSurf_r[1])*fr[7]+((-0.5532833351724881*alphaDrSurf_r[2])-0.8660254037844386*alphaDrSurf_r[0])*fr[6]+1.118033988749895*alphaDrSurf_r[2]*fr[5]+(0.31943828249997*alphaDrSurf_r[2]+0.5*alphaDrSurf_r[0])*fr[4]+((-0.7606388292556648*alphaDrSurf_r[3])-0.7745966692414833*alphaDrSurf_r[1])*fr[3]+0.4391550328268398*fr[1]*alphaDrSurf_r[3]+alphaDrSurf_r[2]*(0.5*fr[0]-0.8660254037844386*fr[2])+0.4472135954999579*alphaDrSurf_r[1]*fr[1]; 
  Ghat_r[3] = (-1.161895003862225*alphaDrSurf_r[2]*fr[11])+((-0.5163977794943223*alphaDrSurf_r[2])-0.8660254037844386*alphaDrSurf_r[0])*fr[10]-1.322875655532295*alphaDrSurf_r[3]*fr[9]+(0.2981423969999719*alphaDrSurf_r[2]+0.5*alphaDrSurf_r[0])*fr[8]+0.9819805060619655*alphaDrSurf_r[2]*fr[7]+((-0.5163977794943223*alphaDrSurf_r[3])-0.7606388292556648*alphaDrSurf_r[1])*fr[6]+1.118033988749895*alphaDrSurf_r[3]*fr[5]+(0.2981423969999719*alphaDrSurf_r[3]+0.4391550328268398*alphaDrSurf_r[1])*fr[4]-0.7606388292556648*alphaDrSurf_r[2]*fr[3]+(0.5*fr[0]-0.8660254037844386*fr[2])*alphaDrSurf_r[3]+0.4391550328268398*fr[1]*alphaDrSurf_r[2]; 

  Ghat_l[0] = (-1.322875655532295*alphaDrSurf_l[1]*fc[11])-0.8660254037844386*alphaDrSurf_l[3]*fc[10]-1.322875655532295*alphaDrSurf_l[0]*fc[9]+0.5*alphaDrSurf_l[3]*fc[8]+1.118033988749895*alphaDrSurf_l[1]*fc[7]-0.8660254037844386*alphaDrSurf_l[2]*fc[6]+1.118033988749895*alphaDrSurf_l[0]*fc[5]+0.5*alphaDrSurf_l[2]*fc[4]-0.8660254037844386*(alphaDrSurf_l[1]*fc[3]+alphaDrSurf_l[0]*fc[2])+0.5*(alphaDrSurf_l[1]*fc[1]+alphaDrSurf_l[0]*fc[0]); 
  Ghat_l[1] = ((-1.183215956619923*alphaDrSurf_l[2])-1.322875655532295*alphaDrSurf_l[0])*fc[11]-0.7606388292556648*alphaDrSurf_l[2]*fc[10]-1.322875655532295*alphaDrSurf_l[1]*fc[9]+0.4391550328268398*alphaDrSurf_l[2]*fc[8]+(alphaDrSurf_l[2]+1.118033988749895*alphaDrSurf_l[0])*fc[7]-0.7606388292556648*alphaDrSurf_l[3]*fc[6]+alphaDrSurf_l[1]*(1.118033988749895*fc[5]-0.7745966692414833*fc[6])+(0.4391550328268398*alphaDrSurf_l[3]+0.4472135954999579*alphaDrSurf_l[1])*fc[4]-0.7745966692414833*alphaDrSurf_l[2]*fc[3]-0.8660254037844386*(alphaDrSurf_l[0]*fc[3]+alphaDrSurf_l[1]*fc[2])+0.4472135954999579*fc[1]*alphaDrSurf_l[2]+0.5*(alphaDrSurf_l[0]*fc[1]+fc[0]*alphaDrSurf_l[1]); 
  Ghat_l[2] = ((-1.161895003862225*alphaDrSurf_l[3])-1.183215956619923*alphaDrSurf_l[1])*fc[11]+((-0.5163977794943223*alphaDrSurf_l[3])-0.7606388292556648*alphaDrSurf_l[1])*fc[10]-1.322875655532295*alphaDrSurf_l[2]*fc[9]+(0.2981423969999719*alphaDrSurf_l[3]+0.4391550328268398*alphaDrSurf_l[1])*fc[8]+(0.9819805060619655*alphaDrSurf_l[3]+alphaDrSurf_l[1])*fc[7]+((-0.5532833351724881*alphaDrSurf_l[2])-0.8660254037844386*alphaDrSurf_l[0])*fc[6]+1.118033988749895*alphaDrSurf_l[2]*fc[5]+(0.31943828249997*alphaDrSurf_l[2]+0.5*alphaDrSurf_l[0])*fc[4]+((-0.7606388292556648*alphaDrSurf_l[3])-0.7745966692414833*alphaDrSurf_l[1])*fc[3]+0.4391550328268398*fc[1]*alphaDrSurf_l[3]+alphaDrSurf_l[2]*(0.5*fc[0]-0.8660254037844386*fc[2])+0.4472135954999579*alphaDrSurf_l[1]*fc[1]; 
  Ghat_l[3] = (-1.161895003862225*alphaDrSurf_l[2]*fc[11])+((-0.5163977794943223*alphaDrSurf_l[2])-0.8660254037844386*alphaDrSurf_l[0])*fc[10]-1.322875655532295*alphaDrSurf_l[3]*fc[9]+(0.2981423969999719*alphaDrSurf_l[2]+0.5*alphaDrSurf_l[0])*fc[8]+0.9819805060619655*alphaDrSurf_l[2]*fc[7]+((-0.5163977794943223*alphaDrSurf_l[3])-0.7606388292556648*alphaDrSurf_l[1])*fc[6]+1.118033988749895*alphaDrSurf_l[3]*fc[5]+(0.2981423969999719*alphaDrSurf_l[3]+0.4391550328268398*alphaDrSurf_l[1])*fc[4]-0.7606388292556648*alphaDrSurf_l[2]*fc[3]+(0.5*fc[0]-0.8660254037844386*fc[2])*alphaDrSurf_l[3]+0.4391550328268398*fc[1]*alphaDrSurf_l[2]; 

  } else { 

  Ghat_r[0] = 1.322875655532295*alphaDrSurf_r[1]*fc[11]+0.8660254037844386*alphaDrSurf_r[3]*fc[10]+1.322875655532295*alphaDrSurf_r[0]*fc[9]+0.5*alphaDrSurf_r[3]*fc[8]+1.118033988749895*alphaDrSurf_r[1]*fc[7]+0.8660254037844387*alphaDrSurf_r[2]*fc[6]+1.118033988749895*alphaDrSurf_r[0]*fc[5]+0.5*alphaDrSurf_r[2]*fc[4]+0.8660254037844386*alphaDrSurf_r[1]*fc[3]+0.8660254037844386*alphaDrSurf_r[0]*fc[2]+0.5*alphaDrSurf_r[1]*fc[1]+0.5*alphaDrSurf_r[0]*fc[0]; 
  Ghat_r[1] = 1.183215956619923*alphaDrSurf_r[2]*fc[11]+1.322875655532295*alphaDrSurf_r[0]*fc[11]+0.7606388292556647*alphaDrSurf_r[2]*fc[10]+1.322875655532295*alphaDrSurf_r[1]*fc[9]+0.4391550328268398*alphaDrSurf_r[2]*fc[8]+1.0*alphaDrSurf_r[2]*fc[7]+1.118033988749895*alphaDrSurf_r[0]*fc[7]+0.7606388292556648*alphaDrSurf_r[3]*fc[6]+0.7745966692414834*alphaDrSurf_r[1]*fc[6]+1.118033988749895*alphaDrSurf_r[1]*fc[5]+0.4391550328268398*alphaDrSurf_r[3]*fc[4]+0.4472135954999579*alphaDrSurf_r[1]*fc[4]+0.7745966692414833*alphaDrSurf_r[2]*fc[3]+0.8660254037844386*alphaDrSurf_r[0]*fc[3]+0.8660254037844386*alphaDrSurf_r[1]*fc[2]+0.4472135954999579*fc[1]*alphaDrSurf_r[2]+0.5*alphaDrSurf_r[0]*fc[1]+0.5*fc[0]*alphaDrSurf_r[1]; 
  Ghat_r[2] = 1.161895003862225*alphaDrSurf_r[3]*fc[11]+1.183215956619923*alphaDrSurf_r[1]*fc[11]+0.5163977794943221*alphaDrSurf_r[3]*fc[10]+0.7606388292556647*alphaDrSurf_r[1]*fc[10]+1.322875655532295*alphaDrSurf_r[2]*fc[9]+0.2981423969999719*alphaDrSurf_r[3]*fc[8]+0.4391550328268398*alphaDrSurf_r[1]*fc[8]+0.9819805060619657*alphaDrSurf_r[3]*fc[7]+1.0*alphaDrSurf_r[1]*fc[7]+0.5532833351724881*alphaDrSurf_r[2]*fc[6]+0.8660254037844387*alphaDrSurf_r[0]*fc[6]+1.118033988749895*alphaDrSurf_r[2]*fc[5]+0.31943828249997*alphaDrSurf_r[2]*fc[4]+0.5*alphaDrSurf_r[0]*fc[4]+0.7606388292556648*alphaDrSurf_r[3]*fc[3]+0.7745966692414833*alphaDrSurf_r[1]*fc[3]+0.4391550328268398*fc[1]*alphaDrSurf_r[3]+0.8660254037844386*alphaDrSurf_r[2]*fc[2]+0.5*fc[0]*alphaDrSurf_r[2]+0.4472135954999579*alphaDrSurf_r[1]*fc[1]; 
  Ghat_r[3] = 1.161895003862225*alphaDrSurf_r[2]*fc[11]+0.5163977794943221*alphaDrSurf_r[2]*fc[10]+0.8660254037844386*alphaDrSurf_r[0]*fc[10]+1.322875655532295*alphaDrSurf_r[3]*fc[9]+0.2981423969999719*alphaDrSurf_r[2]*fc[8]+0.5*alphaDrSurf_r[0]*fc[8]+0.9819805060619657*alphaDrSurf_r[2]*fc[7]+0.5163977794943222*alphaDrSurf_r[3]*fc[6]+0.7606388292556648*alphaDrSurf_r[1]*fc[6]+1.118033988749895*alphaDrSurf_r[3]*fc[5]+0.2981423969999719*alphaDrSurf_r[3]*fc[4]+0.4391550328268398*alphaDrSurf_r[1]*fc[4]+0.7606388292556648*alphaDrSurf_r[2]*fc[3]+0.8660254037844386*fc[2]*alphaDrSurf_r[3]+0.5*fc[0]*alphaDrSurf_r[3]+0.4391550328268398*fc[1]*alphaDrSurf_r[2]; 

  Ghat_l[0] = 1.322875655532295*alphaDrSurf_l[1]*fl[11]+0.8660254037844386*alphaDrSurf_l[3]*fl[10]+1.322875655532295*alphaDrSurf_l[0]*fl[9]+0.5*alphaDrSurf_l[3]*fl[8]+1.118033988749895*alphaDrSurf_l[1]*fl[7]+0.8660254037844387*alphaDrSurf_l[2]*fl[6]+1.118033988749895*alphaDrSurf_l[0]*fl[5]+0.5*alphaDrSurf_l[2]*fl[4]+0.8660254037844386*alphaDrSurf_l[1]*fl[3]+0.8660254037844386*alphaDrSurf_l[0]*fl[2]+0.5*alphaDrSurf_l[1]*fl[1]+0.5*alphaDrSurf_l[0]*fl[0]; 
  Ghat_l[1] = 1.183215956619923*alphaDrSurf_l[2]*fl[11]+1.322875655532295*alphaDrSurf_l[0]*fl[11]+0.7606388292556647*alphaDrSurf_l[2]*fl[10]+1.322875655532295*alphaDrSurf_l[1]*fl[9]+0.4391550328268398*alphaDrSurf_l[2]*fl[8]+1.0*alphaDrSurf_l[2]*fl[7]+1.118033988749895*alphaDrSurf_l[0]*fl[7]+0.7606388292556648*alphaDrSurf_l[3]*fl[6]+0.7745966692414834*alphaDrSurf_l[1]*fl[6]+1.118033988749895*alphaDrSurf_l[1]*fl[5]+0.4391550328268398*alphaDrSurf_l[3]*fl[4]+0.4472135954999579*alphaDrSurf_l[1]*fl[4]+0.7745966692414833*alphaDrSurf_l[2]*fl[3]+0.8660254037844386*alphaDrSurf_l[0]*fl[3]+0.8660254037844386*alphaDrSurf_l[1]*fl[2]+0.4472135954999579*fl[1]*alphaDrSurf_l[2]+0.5*alphaDrSurf_l[0]*fl[1]+0.5*fl[0]*alphaDrSurf_l[1]; 
  Ghat_l[2] = 1.161895003862225*alphaDrSurf_l[3]*fl[11]+1.183215956619923*alphaDrSurf_l[1]*fl[11]+0.5163977794943221*alphaDrSurf_l[3]*fl[10]+0.7606388292556647*alphaDrSurf_l[1]*fl[10]+1.322875655532295*alphaDrSurf_l[2]*fl[9]+0.2981423969999719*alphaDrSurf_l[3]*fl[8]+0.4391550328268398*alphaDrSurf_l[1]*fl[8]+0.9819805060619657*alphaDrSurf_l[3]*fl[7]+1.0*alphaDrSurf_l[1]*fl[7]+0.5532833351724881*alphaDrSurf_l[2]*fl[6]+0.8660254037844387*alphaDrSurf_l[0]*fl[6]+1.118033988749895*alphaDrSurf_l[2]*fl[5]+0.31943828249997*alphaDrSurf_l[2]*fl[4]+0.5*alphaDrSurf_l[0]*fl[4]+0.7606388292556648*alphaDrSurf_l[3]*fl[3]+0.7745966692414833*alphaDrSurf_l[1]*fl[3]+0.4391550328268398*fl[1]*alphaDrSurf_l[3]+0.8660254037844386*alphaDrSurf_l[2]*fl[2]+0.5*fl[0]*alphaDrSurf_l[2]+0.4472135954999579*alphaDrSurf_l[1]*fl[1]; 
  Ghat_l[3] = 1.161895003862225*alphaDrSurf_l[2]*fl[11]+0.5163977794943221*alphaDrSurf_l[2]*fl[10]+0.8660254037844386*alphaDrSurf_l[0]*fl[10]+1.322875655532295*alphaDrSurf_l[3]*fl[9]+0.2981423969999719*alphaDrSurf_l[2]*fl[8]+0.5*alphaDrSurf_l[0]*fl[8]+0.9819805060619657*alphaDrSurf_l[2]*fl[7]+0.5163977794943222*alphaDrSurf_l[3]*fl[6]+0.7606388292556648*alphaDrSurf_l[1]*fl[6]+1.118033988749895*alphaDrSurf_l[3]*fl[5]+0.2981423969999719*alphaDrSurf_l[3]*fl[4]+0.4391550328268398*alphaDrSurf_l[1]*fl[4]+0.7606388292556648*alphaDrSurf_l[2]*fl[3]+0.8660254037844386*fl[2]*alphaDrSurf_l[3]+0.5*fl[0]*alphaDrSurf_l[3]+0.4391550328268398*fl[1]*alphaDrSurf_l[2]; 

  } 
  out[0] += (0.7071067811865475*Ghat_r[0]-0.7071067811865475*Ghat_l[0])*dv1par; 
  out[1] += (0.7071067811865475*Ghat_r[1]-0.7071067811865475*Ghat_l[1])*dv1par; 
  out[2] += 1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv1par; 
  out[3] += 1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv1par; 
  out[4] += (0.7071067811865475*Ghat_r[2]-0.7071067811865475*Ghat_l[2])*dv1par; 
  out[5] += (1.58113883008419*Ghat_r[0]-1.58113883008419*Ghat_l[0])*dv1par; 
  out[6] += 1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv1par; 
  out[7] += (1.58113883008419*Ghat_r[1]-1.58113883008419*Ghat_l[1])*dv1par; 
  out[8] += (0.7071067811865475*Ghat_r[3]-0.7071067811865475*Ghat_l[3])*dv1par; 
  out[9] += 1.870828693386971*(Ghat_r[0]+Ghat_l[0])*dv1par; 
  out[10] += 1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv1par; 
  out[11] += 1.870828693386971*(Ghat_r[1]+Ghat_l[1])*dv1par; 
} 
