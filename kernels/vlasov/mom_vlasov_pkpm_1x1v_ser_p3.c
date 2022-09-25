#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void mom_vlasov_pkpm_1x1v_ser_p3(const double *w, const double *dxv, const int *idx, double mass, const double *bvar, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]/2.0; 
  const double wvpar = w[1], dvpar = dxv[1]; 
  const double wvpar_sq = wvpar*wvpar, dvpar_sq = dvpar*dvpar; 
  const double wvpar_cu = wvpar*wvpar*wvpar, dvpar_cu = dvpar*dvpar*dvpar; 

  const double *bx = &bvar[0]; 
  const double *by = &bvar[4]; 
  const double *bz = &bvar[8]; 

  out[0] += 1.414213562373095*f[0]*mass*volFact; 
  out[1] += 1.414213562373095*f[1]*mass*volFact; 
  out[2] += 1.414213562373095*f[4]*mass*volFact; 
  out[3] += 1.414213562373095*f[8]*mass*volFact; 
  out[4] += mass*volFact*(1.414213562373095*f[0]*wvpar_sq+0.8164965809277261*f[2]*dvpar*wvpar+0.105409255338946*f[5]*dvpar_sq+0.1178511301977579*f[0]*dvpar_sq); 
  out[5] += mass*volFact*(1.414213562373095*f[1]*wvpar_sq+0.8164965809277261*f[3]*dvpar*wvpar+0.105409255338946*f[7]*dvpar_sq+0.1178511301977579*f[1]*dvpar_sq); 
  out[6] += mass*volFact*(1.414213562373095*f[4]*wvpar_sq+0.816496580927726*f[6]*dvpar*wvpar+0.1178511301977579*f[4]*dvpar_sq); 
  out[7] += mass*volFact*(1.414213562373095*f[8]*wvpar_sq+0.8164965809277258*f[10]*dvpar*wvpar+0.1178511301977579*f[8]*dvpar_sq); 
  out[8] += mass*volFact*(0.8660254037844386*bx[3]*f[10]*dvpar*wvpar_sq+0.8660254037844387*bx[2]*f[6]*dvpar*wvpar_sq+0.8660254037844386*bx[1]*f[3]*dvpar*wvpar_sq+0.8660254037844386*bx[0]*f[2]*dvpar*wvpar_sq+bx[3]*f[8]*wvpar_cu+bx[2]*f[4]*wvpar_cu+bx[1]*f[1]*wvpar_cu+bx[0]*f[0]*wvpar_cu+0.25*bx[3]*f[8]*dvpar_sq*wvpar+0.223606797749979*bx[1]*f[7]*dvpar_sq*wvpar+0.223606797749979*bx[0]*f[5]*dvpar_sq*wvpar+0.25*bx[2]*f[4]*dvpar_sq*wvpar+0.25*bx[1]*f[1]*dvpar_sq*wvpar+0.25*bx[0]*f[0]*dvpar_sq*wvpar+0.01889822365046136*bx[1]*f[11]*dvpar_cu+0.04330127018922193*bx[3]*f[10]*dvpar_cu+0.01889822365046136*bx[0]*f[9]*dvpar_cu+0.04330127018922193*bx[2]*f[6]*dvpar_cu+0.04330127018922193*bx[1]*f[3]*dvpar_cu+0.04330127018922193*bx[0]*f[2]*dvpar_cu); 
  out[9] += mass*volFact*(0.7606388292556647*bx[2]*f[10]*dvpar*wvpar_sq+0.7606388292556648*bx[3]*f[6]*dvpar*wvpar_sq+0.7745966692414834*bx[1]*f[6]*dvpar*wvpar_sq+0.7745966692414833*bx[2]*f[3]*dvpar*wvpar_sq+0.8660254037844386*bx[0]*f[3]*dvpar*wvpar_sq+0.8660254037844386*bx[1]*f[2]*dvpar*wvpar_sq+0.8783100656536796*bx[2]*f[8]*wvpar_cu+0.8783100656536796*bx[3]*f[4]*wvpar_cu+0.8944271909999159*bx[1]*f[4]*wvpar_cu+0.8944271909999159*f[1]*bx[2]*wvpar_cu+bx[0]*f[1]*wvpar_cu+f[0]*bx[1]*wvpar_cu+0.2195775164134199*bx[2]*f[8]*dvpar_sq*wvpar+0.2*bx[2]*f[7]*dvpar_sq*wvpar+0.223606797749979*bx[0]*f[7]*dvpar_sq*wvpar+0.223606797749979*bx[1]*f[5]*dvpar_sq*wvpar+0.2195775164134199*bx[3]*f[4]*dvpar_sq*wvpar+0.223606797749979*bx[1]*f[4]*dvpar_sq*wvpar+0.223606797749979*f[1]*bx[2]*dvpar_sq*wvpar+0.25*bx[0]*f[1]*dvpar_sq*wvpar+0.25*f[0]*bx[1]*dvpar_sq*wvpar+0.01690308509457033*bx[2]*f[11]*dvpar_cu+0.01889822365046136*bx[0]*f[11]*dvpar_cu+0.03803194146278323*bx[2]*f[10]*dvpar_cu+0.01889822365046136*bx[1]*f[9]*dvpar_cu+0.03803194146278324*bx[3]*f[6]*dvpar_cu+0.03872983346207417*bx[1]*f[6]*dvpar_cu+0.03872983346207416*bx[2]*f[3]*dvpar_cu+0.04330127018922193*bx[0]*f[3]*dvpar_cu+0.04330127018922193*bx[1]*f[2]*dvpar_cu); 
  out[10] += mass*volFact*(0.5163977794943221*bx[3]*f[10]*dvpar*wvpar_sq+0.7606388292556647*bx[1]*f[10]*dvpar*wvpar_sq+0.5532833351724881*bx[2]*f[6]*dvpar*wvpar_sq+0.8660254037844387*bx[0]*f[6]*dvpar*wvpar_sq+0.7606388292556648*bx[3]*f[3]*dvpar*wvpar_sq+0.7745966692414833*bx[1]*f[3]*dvpar*wvpar_sq+0.8660254037844386*bx[2]*f[2]*dvpar*wvpar_sq+0.5962847939999438*bx[3]*f[8]*wvpar_cu+0.8783100656536796*bx[1]*f[8]*wvpar_cu+0.6388765649999399*bx[2]*f[4]*wvpar_cu+bx[0]*f[4]*wvpar_cu+0.8783100656536796*f[1]*bx[3]*wvpar_cu+f[0]*bx[2]*wvpar_cu+0.8944271909999159*bx[1]*f[1]*wvpar_cu+0.149071198499986*bx[3]*f[8]*dvpar_sq*wvpar+0.2195775164134199*bx[1]*f[8]*dvpar_sq*wvpar+0.1963961012123931*bx[3]*f[7]*dvpar_sq*wvpar+0.2*bx[1]*f[7]*dvpar_sq*wvpar+0.223606797749979*bx[2]*f[5]*dvpar_sq*wvpar+0.159719141249985*bx[2]*f[4]*dvpar_sq*wvpar+0.25*bx[0]*f[4]*dvpar_sq*wvpar+0.2195775164134199*f[1]*bx[3]*dvpar_sq*wvpar+0.25*f[0]*bx[2]*dvpar_sq*wvpar+0.223606797749979*bx[1]*f[1]*dvpar_sq*wvpar+0.01659850005517464*bx[3]*f[11]*dvpar_cu+0.01690308509457033*bx[1]*f[11]*dvpar_cu+0.02581988897471611*bx[3]*f[10]*dvpar_cu+0.03803194146278323*bx[1]*f[10]*dvpar_cu+0.01889822365046136*bx[2]*f[9]*dvpar_cu+0.02766416675862441*bx[2]*f[6]*dvpar_cu+0.04330127018922193*bx[0]*f[6]*dvpar_cu+0.03803194146278324*bx[3]*f[3]*dvpar_cu+0.03872983346207416*bx[1]*f[3]*dvpar_cu+0.04330127018922193*bx[2]*f[2]*dvpar_cu); 
  out[11] += mass*volFact*(0.5163977794943221*bx[2]*f[10]*dvpar*wvpar_sq+0.8660254037844386*bx[0]*f[10]*dvpar*wvpar_sq+0.5163977794943222*bx[3]*f[6]*dvpar*wvpar_sq+0.7606388292556648*bx[1]*f[6]*dvpar*wvpar_sq+0.7606388292556648*bx[2]*f[3]*dvpar*wvpar_sq+0.8660254037844386*f[2]*bx[3]*dvpar*wvpar_sq+0.5962847939999438*bx[2]*f[8]*wvpar_cu+bx[0]*f[8]*wvpar_cu+0.5962847939999438*bx[3]*f[4]*wvpar_cu+0.8783100656536796*bx[1]*f[4]*wvpar_cu+f[0]*bx[3]*wvpar_cu+0.8783100656536796*f[1]*bx[2]*wvpar_cu+0.149071198499986*bx[2]*f[8]*dvpar_sq*wvpar+0.25*bx[0]*f[8]*dvpar_sq*wvpar+0.1963961012123931*bx[2]*f[7]*dvpar_sq*wvpar+0.223606797749979*bx[3]*f[5]*dvpar_sq*wvpar+0.149071198499986*bx[3]*f[4]*dvpar_sq*wvpar+0.2195775164134199*bx[1]*f[4]*dvpar_sq*wvpar+0.25*f[0]*bx[3]*dvpar_sq*wvpar+0.2195775164134199*f[1]*bx[2]*dvpar_sq*wvpar+0.01659850005517464*bx[2]*f[11]*dvpar_cu+0.02581988897471611*bx[2]*f[10]*dvpar_cu+0.04330127018922193*bx[0]*f[10]*dvpar_cu+0.01889822365046136*bx[3]*f[9]*dvpar_cu+0.02581988897471611*bx[3]*f[6]*dvpar_cu+0.03803194146278324*bx[1]*f[6]*dvpar_cu+0.03803194146278324*bx[2]*f[3]*dvpar_cu+0.04330127018922193*f[2]*bx[3]*dvpar_cu); 
} 
