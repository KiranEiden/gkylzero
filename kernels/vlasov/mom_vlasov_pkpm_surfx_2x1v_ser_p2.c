#include <gkyl_mom_vlasov_kernels.h> 
#include <gkyl_basis_ser_3x_p2_surfx1_eval_quad.h> 
#include <gkyl_basis_ser_3x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH void mom_vlasov_pkpm_surfx_2x1v_ser_p2(const double *w, const double *dxv, const int *idx, double mass, const double *u_il, const double *u_ir, const double *bvarl, const double *bvarr, const double *fl, const double *fr, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]/2.0; 
  const double wvpar = w[0], dvpar = dxv[0]; 
  const double wvpar_sq = wvpar*wvpar, dvpar_sq = dvpar*dvpar; 
  const double wvpar_cu = wvpar*wvpar*wvpar, dvpar_cu = dvpar*dvpar*dvpar; 

  const double *ul_0 = &u_il[0]; 
  const double *ur_0 = &u_ir[0]; 
  const double *bl_0 = &bvarl[0]; 
  const double *br_0 = &bvarr[0]; 
  double *rho_flux = &out[0]; 
  double *heat_flux = &out[6]; 
  double u_l[8] = {0.0}; 
  double u_r[8] = {0.0}; 
  double q_l[8] = {0.0}; 
  double q_r[8] = {0.0}; 

  u_l[0] = 2.23606797749979*ul_0[4]+1.732050807568877*ul_0[1]+ul_0[0]; 
  u_l[1] = 2.23606797749979*ul_0[6]+1.732050807568877*ul_0[3]+ul_0[2]; 
  u_l[4] = 1.732050807568877*ul_0[7]+ul_0[5]; 

  u_r[0] = 2.23606797749979*ur_0[4]-1.732050807568877*ur_0[1]+ur_0[0]; 
  u_r[1] = 2.23606797749979*ur_0[6]-1.732050807568877*ur_0[3]+ur_0[2]; 
  u_r[4] = ur_0[5]-1.732050807568877*ur_0[7]; 

  q_l[0] = 1.118033988749895*bl_0[4]*wvpar_cu+0.8660254037844386*bl_0[1]*wvpar_cu+0.5*bl_0[0]*wvpar_cu+0.2795084971874737*bl_0[4]*dvpar_sq*wvpar+0.2165063509461096*bl_0[1]*dvpar_sq*wvpar+0.125*bl_0[0]*dvpar_sq*wvpar; 
  q_l[1] = 1.118033988749895*bl_0[6]*wvpar_cu+0.8660254037844386*bl_0[3]*wvpar_cu+0.5*bl_0[2]*wvpar_cu+0.2795084971874738*bl_0[6]*dvpar_sq*wvpar+0.2165063509461096*bl_0[3]*dvpar_sq*wvpar+0.125*bl_0[2]*dvpar_sq*wvpar; 
  q_l[2] = 0.9682458365518543*bl_0[4]*dvpar*wvpar_sq+0.75*bl_0[1]*dvpar*wvpar_sq+0.4330127018922193*bl_0[0]*dvpar*wvpar_sq+0.04841229182759271*bl_0[4]*dvpar_cu+0.0375*bl_0[1]*dvpar_cu+0.02165063509461097*bl_0[0]*dvpar_cu; 
  q_l[3] = 0.9682458365518543*bl_0[6]*dvpar*wvpar_sq+0.75*bl_0[3]*dvpar*wvpar_sq+0.4330127018922193*bl_0[2]*dvpar*wvpar_sq+0.04841229182759271*bl_0[6]*dvpar_cu+0.0375*bl_0[3]*dvpar_cu+0.02165063509461097*bl_0[2]*dvpar_cu; 
  q_l[4] = 0.8660254037844387*bl_0[7]*wvpar_cu+0.5*bl_0[5]*wvpar_cu+0.2165063509461097*bl_0[7]*dvpar_sq*wvpar+0.125*bl_0[5]*dvpar_sq*wvpar; 
  q_l[5] = 0.25*bl_0[4]*dvpar_sq*wvpar+0.1936491673103708*bl_0[1]*dvpar_sq*wvpar+0.1118033988749895*bl_0[0]*dvpar_sq*wvpar; 
  q_l[6] = 0.75*bl_0[7]*dvpar*wvpar_sq+0.4330127018922194*bl_0[5]*dvpar*wvpar_sq+0.0375*bl_0[7]*dvpar_cu+0.02165063509461096*bl_0[5]*dvpar_cu; 
  q_l[7] = 0.25*bl_0[6]*dvpar_sq*wvpar+0.1936491673103709*bl_0[3]*dvpar_sq*wvpar+0.1118033988749895*bl_0[2]*dvpar_sq*wvpar; 

  q_r[0] = 1.118033988749895*br_0[4]*wvpar_cu-0.8660254037844386*br_0[1]*wvpar_cu+0.5*br_0[0]*wvpar_cu+0.2795084971874737*br_0[4]*dvpar_sq*wvpar-0.2165063509461096*br_0[1]*dvpar_sq*wvpar+0.125*br_0[0]*dvpar_sq*wvpar; 
  q_r[1] = 1.118033988749895*br_0[6]*wvpar_cu-0.8660254037844386*br_0[3]*wvpar_cu+0.5*br_0[2]*wvpar_cu+0.2795084971874738*br_0[6]*dvpar_sq*wvpar-0.2165063509461096*br_0[3]*dvpar_sq*wvpar+0.125*br_0[2]*dvpar_sq*wvpar; 
  q_r[2] = 0.9682458365518543*br_0[4]*dvpar*wvpar_sq-0.75*br_0[1]*dvpar*wvpar_sq+0.4330127018922193*br_0[0]*dvpar*wvpar_sq+0.04841229182759271*br_0[4]*dvpar_cu-0.0375*br_0[1]*dvpar_cu+0.02165063509461097*br_0[0]*dvpar_cu; 
  q_r[3] = 0.9682458365518543*br_0[6]*dvpar*wvpar_sq-0.75*br_0[3]*dvpar*wvpar_sq+0.4330127018922193*br_0[2]*dvpar*wvpar_sq+0.04841229182759271*br_0[6]*dvpar_cu-0.0375*br_0[3]*dvpar_cu+0.02165063509461097*br_0[2]*dvpar_cu; 
  q_r[4] = (-0.8660254037844387*br_0[7]*wvpar_cu)+0.5*br_0[5]*wvpar_cu-0.2165063509461097*br_0[7]*dvpar_sq*wvpar+0.125*br_0[5]*dvpar_sq*wvpar; 
  q_r[5] = 0.25*br_0[4]*dvpar_sq*wvpar-0.1936491673103708*br_0[1]*dvpar_sq*wvpar+0.1118033988749895*br_0[0]*dvpar_sq*wvpar; 
  q_r[6] = (-0.75*br_0[7]*dvpar*wvpar_sq)+0.4330127018922194*br_0[5]*dvpar*wvpar_sq-0.0375*br_0[7]*dvpar_cu+0.02165063509461096*br_0[5]*dvpar_cu; 
  q_r[7] = 0.25*br_0[6]*dvpar_sq*wvpar-0.1936491673103709*br_0[3]*dvpar_sq*wvpar+0.1118033988749895*br_0[2]*dvpar_sq*wvpar; 

  double uQuad[9] = {0.0};
  double qQuad[9] = {0.0};
  double uMax[8] = {0.0};
  double qMax[8] = {0.0};
  double ev_u_l = 0.0; 
  double ev_u_r = 0.0; 
  double ev_q_l = 0.0; 
  double ev_q_r = 0.0; 

  ev_u_l = 0.4472135954999579*u_l[4]-0.6708203932499369*u_l[1]+0.5*u_l[0]; 
  ev_u_r = 0.4472135954999579*u_r[4]-0.6708203932499369*u_r[1]+0.5*u_r[0]; 
  ev_q_l = (-0.5999999999999995*q_l[7])-0.5999999999999999*q_l[6]+0.4472135954999579*(q_l[5]+q_l[4])+0.9*q_l[3]-0.6708203932499369*(q_l[2]+q_l[1])+0.5*q_l[0]; 
  ev_q_r = (-0.5999999999999995*q_l[7])-0.5999999999999999*q_l[6]+0.4472135954999579*(q_l[5]+q_l[4])+0.9*q_l[3]-0.6708203932499369*(q_l[2]+q_l[1])+0.5*q_l[0]; 
  uQuad[0] = fmax(fabs(ev_u_l), fabs(ev_u_r)); 
  qQuad[0] = fmax(fabs(ev_q_l), fabs(ev_q_r)); 
  ev_u_l = 0.4472135954999579*u_l[4]-0.6708203932499369*u_l[1]+0.5*u_l[0]; 
  ev_u_r = 0.4472135954999579*u_r[4]-0.6708203932499369*u_r[1]+0.5*u_r[0]; 
  ev_q_l = 0.75*q_l[7]-0.5590169943749475*q_l[5]+0.4472135954999579*q_l[4]-0.6708203932499369*q_l[1]+0.5*q_l[0]; 
  ev_q_r = 0.75*q_l[7]-0.5590169943749475*q_l[5]+0.4472135954999579*q_l[4]-0.6708203932499369*q_l[1]+0.5*q_l[0]; 
  uQuad[1] = fmax(fabs(ev_u_l), fabs(ev_u_r)); 
  qQuad[1] = fmax(fabs(ev_q_l), fabs(ev_q_r)); 
  ev_u_l = 0.4472135954999579*u_l[4]-0.6708203932499369*u_l[1]+0.5*u_l[0]; 
  ev_u_r = 0.4472135954999579*u_r[4]-0.6708203932499369*u_r[1]+0.5*u_r[0]; 
  ev_q_l = (-0.5999999999999995*q_l[7])+0.5999999999999999*q_l[6]+0.4472135954999579*(q_l[5]+q_l[4])-0.9*q_l[3]+0.6708203932499369*q_l[2]-0.6708203932499369*q_l[1]+0.5*q_l[0]; 
  ev_q_r = (-0.5999999999999995*q_l[7])+0.5999999999999999*q_l[6]+0.4472135954999579*(q_l[5]+q_l[4])-0.9*q_l[3]+0.6708203932499369*q_l[2]-0.6708203932499369*q_l[1]+0.5*q_l[0]; 
  uQuad[2] = fmax(fabs(ev_u_l), fabs(ev_u_r)); 
  qQuad[2] = fmax(fabs(ev_q_l), fabs(ev_q_r)); 
  ev_u_l = 0.5*u_l[0]-0.5590169943749475*u_l[4]; 
  ev_u_r = 0.5*u_r[0]-0.5590169943749475*u_r[4]; 
  ev_q_l = 0.75*q_l[6]+0.4472135954999579*q_l[5]-0.5590169943749475*q_l[4]-0.6708203932499369*q_l[2]+0.5*q_l[0]; 
  ev_q_r = 0.75*q_l[6]+0.4472135954999579*q_l[5]-0.5590169943749475*q_l[4]-0.6708203932499369*q_l[2]+0.5*q_l[0]; 
  uQuad[3] = fmax(fabs(ev_u_l), fabs(ev_u_r)); 
  qQuad[3] = fmax(fabs(ev_q_l), fabs(ev_q_r)); 
  ev_u_l = 0.5*u_l[0]-0.5590169943749475*u_l[4]; 
  ev_u_r = 0.5*u_r[0]-0.5590169943749475*u_r[4]; 
  ev_q_l = 0.5*q_l[0]-0.5590169943749475*(q_l[5]+q_l[4]); 
  ev_q_r = 0.5*q_l[0]-0.5590169943749475*(q_l[5]+q_l[4]); 
  uQuad[4] = fmax(fabs(ev_u_l), fabs(ev_u_r)); 
  qQuad[4] = fmax(fabs(ev_q_l), fabs(ev_q_r)); 
  ev_u_l = 0.5*u_l[0]-0.5590169943749475*u_l[4]; 
  ev_u_r = 0.5*u_r[0]-0.5590169943749475*u_r[4]; 
  ev_q_l = (-0.75*q_l[6])+0.4472135954999579*q_l[5]-0.5590169943749475*q_l[4]+0.6708203932499369*q_l[2]+0.5*q_l[0]; 
  ev_q_r = (-0.75*q_l[6])+0.4472135954999579*q_l[5]-0.5590169943749475*q_l[4]+0.6708203932499369*q_l[2]+0.5*q_l[0]; 
  uQuad[5] = fmax(fabs(ev_u_l), fabs(ev_u_r)); 
  qQuad[5] = fmax(fabs(ev_q_l), fabs(ev_q_r)); 
  ev_u_l = 0.4472135954999579*u_l[4]+0.6708203932499369*u_l[1]+0.5*u_l[0]; 
  ev_u_r = 0.4472135954999579*u_r[4]+0.6708203932499369*u_r[1]+0.5*u_r[0]; 
  ev_q_l = 0.5999999999999995*q_l[7]-0.5999999999999999*q_l[6]+0.4472135954999579*(q_l[5]+q_l[4])-0.9*q_l[3]-0.6708203932499369*q_l[2]+0.6708203932499369*q_l[1]+0.5*q_l[0]; 
  ev_q_r = 0.5999999999999995*q_l[7]-0.5999999999999999*q_l[6]+0.4472135954999579*(q_l[5]+q_l[4])-0.9*q_l[3]-0.6708203932499369*q_l[2]+0.6708203932499369*q_l[1]+0.5*q_l[0]; 
  uQuad[6] = fmax(fabs(ev_u_l), fabs(ev_u_r)); 
  qQuad[6] = fmax(fabs(ev_q_l), fabs(ev_q_r)); 
  ev_u_l = 0.4472135954999579*u_l[4]+0.6708203932499369*u_l[1]+0.5*u_l[0]; 
  ev_u_r = 0.4472135954999579*u_r[4]+0.6708203932499369*u_r[1]+0.5*u_r[0]; 
  ev_q_l = (-0.75*q_l[7])-0.5590169943749475*q_l[5]+0.4472135954999579*q_l[4]+0.6708203932499369*q_l[1]+0.5*q_l[0]; 
  ev_q_r = (-0.75*q_l[7])-0.5590169943749475*q_l[5]+0.4472135954999579*q_l[4]+0.6708203932499369*q_l[1]+0.5*q_l[0]; 
  uQuad[7] = fmax(fabs(ev_u_l), fabs(ev_u_r)); 
  qQuad[7] = fmax(fabs(ev_q_l), fabs(ev_q_r)); 
  ev_u_l = 0.4472135954999579*u_l[4]+0.6708203932499369*u_l[1]+0.5*u_l[0]; 
  ev_u_r = 0.4472135954999579*u_r[4]+0.6708203932499369*u_r[1]+0.5*u_r[0]; 
  ev_q_l = 0.5999999999999995*q_l[7]+0.5999999999999999*q_l[6]+0.4472135954999579*(q_l[5]+q_l[4])+0.9*q_l[3]+0.6708203932499369*(q_l[2]+q_l[1])+0.5*q_l[0]; 
  ev_q_r = 0.5999999999999995*q_l[7]+0.5999999999999999*q_l[6]+0.4472135954999579*(q_l[5]+q_l[4])+0.9*q_l[3]+0.6708203932499369*(q_l[2]+q_l[1])+0.5*q_l[0]; 
  uQuad[8] = fmax(fabs(ev_u_l), fabs(ev_u_r)); 
  qQuad[8] = fmax(fabs(ev_q_l), fabs(ev_q_r)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_3x_p2_upwind_quad_to_modal(uQuad, uMax); 
  ser_3x_p2_upwind_quad_to_modal(qQuad, qMax); 
  rho_flux[0] += (0.4330127018922194*uMax[7]*fr[19]+0.4330127018922194*uMax[7]*fl[19]+0.4330127018922194*uMax[6]*fr[18]+0.4330127018922194*uMax[6]*fl[18]-0.5590169943749475*uMax[3]*fr[17]+0.5590169943749475*uMax[3]*fl[17]-0.25*uMax[7]*fr[16]+0.25*uMax[7]*fl[16]+0.4330127018922194*uMax[5]*fr[15]+0.4330127018922194*uMax[5]*fl[15]-0.25*uMax[6]*fr[14]+0.25*uMax[6]*fl[14]-0.5590169943749476*uMax[2]*fr[13]+0.5590169943749476*uMax[2]*fl[13]-0.4330127018922194*u_r[4]*fr[12]+0.4330127018922194*uMax[4]*fr[12]+0.4330127018922194*u_l[4]*fl[12]+0.4330127018922194*uMax[4]*fl[12]+0.5590169943749476*u_r[1]*fr[11]-0.5590169943749476*uMax[1]*fr[11]+0.5590169943749476*u_l[1]*fl[11]+0.5590169943749476*uMax[1]*fl[11]+0.4330127018922193*uMax[3]*fr[10]+0.4330127018922193*uMax[3]*fl[10]-0.25*uMax[5]*fr[9]+0.25*uMax[5]*fl[9]+0.25*u_r[4]*fr[8]-0.25*uMax[4]*fr[8]+0.25*u_l[4]*fl[8]+0.25*uMax[4]*fl[8]+0.5590169943749475*u_r[0]*fr[7]-0.5590169943749475*uMax[0]*fr[7]+0.5590169943749475*u_l[0]*fl[7]+0.5590169943749475*uMax[0]*fl[7]-0.25*uMax[3]*fr[6]+0.25*uMax[3]*fl[6]+0.4330127018922193*uMax[2]*fr[5]+0.4330127018922193*uMax[2]*fl[5]-0.4330127018922193*u_r[1]*fr[4]+0.4330127018922193*uMax[1]*fr[4]+0.4330127018922193*u_l[1]*fl[4]+0.4330127018922193*uMax[1]*fl[4]-0.25*uMax[2]*fr[3]+0.25*uMax[2]*fl[3]+0.25*u_r[1]*fr[2]-0.25*uMax[1]*fr[2]+0.25*u_l[1]*fl[2]+0.25*uMax[1]*fl[2]-0.4330127018922193*u_r[0]*fr[1]+0.4330127018922193*uMax[0]*fr[1]+0.4330127018922193*u_l[0]*fl[1]+0.4330127018922193*uMax[0]*fl[1]+0.25*fr[0]*u_r[0]+0.25*fl[0]*u_l[0]-0.25*fr[0]*uMax[0]+0.25*fl[0]*uMax[0])*mass*volFact; 
  rho_flux[1] += (0.4330127018922193*uMax[5]*fr[19]+0.4330127018922193*uMax[5]*fl[19]+0.3872983346207416*uMax[3]*fr[18]+0.3872983346207416*uMax[3]*fl[18]-0.5000000000000001*uMax[6]*fr[17]-0.5590169943749475*uMax[2]*fr[17]+0.5000000000000001*uMax[6]*fl[17]+0.5590169943749475*uMax[2]*fl[17]-0.2500000000000001*uMax[5]*fr[16]+0.2500000000000001*uMax[5]*fl[16]+0.4330127018922193*uMax[7]*fr[15]+0.4330127018922193*uMax[7]*fl[15]-0.223606797749979*uMax[3]*fr[14]+0.223606797749979*uMax[3]*fl[14]-0.5590169943749476*uMax[3]*fr[13]+0.5590169943749476*uMax[3]*fl[13]-0.3872983346207417*u_r[1]*fr[12]+0.3872983346207417*uMax[1]*fr[12]+0.3872983346207417*u_l[1]*fl[12]+0.3872983346207417*uMax[1]*fl[12]+0.5000000000000001*u_r[4]*fr[11]-0.5000000000000001*uMax[4]*fr[11]+0.5590169943749476*u_r[0]*fr[11]-0.5590169943749476*uMax[0]*fr[11]+0.5000000000000001*u_l[4]*fl[11]+0.5000000000000001*uMax[4]*fl[11]+0.5590169943749476*u_l[0]*fl[11]+0.5590169943749476*uMax[0]*fl[11]+0.3872983346207417*uMax[6]*fr[10]+0.4330127018922193*uMax[2]*fr[10]+0.3872983346207417*uMax[6]*fl[10]+0.4330127018922193*uMax[2]*fl[10]-0.2500000000000001*uMax[7]*fr[9]+0.2500000000000001*uMax[7]*fl[9]+0.223606797749979*u_r[1]*fr[8]-0.223606797749979*uMax[1]*fr[8]+0.223606797749979*u_l[1]*fl[8]+0.223606797749979*uMax[1]*fl[8]+0.5590169943749475*u_r[1]*fr[7]-0.5590169943749475*uMax[1]*fr[7]+0.5590169943749475*u_l[1]*fl[7]+0.5590169943749475*uMax[1]*fl[7]-0.223606797749979*fr[6]*uMax[6]+0.223606797749979*fl[6]*uMax[6]-0.25*uMax[2]*fr[6]+0.25*uMax[2]*fl[6]+0.4330127018922193*uMax[3]*fr[5]+0.4330127018922193*uMax[3]*fl[5]-0.3872983346207416*fr[4]*u_r[4]+0.223606797749979*fr[2]*u_r[4]+0.3872983346207416*fl[4]*u_l[4]+0.223606797749979*fl[2]*u_l[4]+0.3872983346207416*fr[4]*uMax[4]+0.3872983346207416*fl[4]*uMax[4]-0.223606797749979*fr[2]*uMax[4]+0.223606797749979*fl[2]*uMax[4]-0.4330127018922193*u_r[0]*fr[4]+0.4330127018922193*uMax[0]*fr[4]+0.4330127018922193*u_l[0]*fl[4]+0.4330127018922193*uMax[0]*fl[4]-0.25*fr[3]*uMax[3]+0.25*fl[3]*uMax[3]+0.25*u_r[0]*fr[2]-0.25*uMax[0]*fr[2]+0.25*u_l[0]*fl[2]+0.25*uMax[0]*fl[2]-0.4330127018922193*fr[1]*u_r[1]+0.25*fr[0]*u_r[1]+0.4330127018922193*fl[1]*u_l[1]+0.25*fl[0]*u_l[1]+0.4330127018922193*fr[1]*uMax[1]+0.4330127018922193*fl[1]*uMax[1]-0.25*fr[0]*uMax[1]+0.25*fl[0]*uMax[1])*mass*volFact; 
  rho_flux[2] += (0.3872983346207417*uMax[7]*fr[19]+0.3872983346207417*uMax[7]*fl[19]+0.276641667586244*uMax[6]*fr[18]+0.4330127018922193*uMax[2]*fr[18]+0.276641667586244*uMax[6]*fl[18]+0.4330127018922193*uMax[2]*fl[18]-0.5*uMax[3]*fr[17]+0.5*uMax[3]*fl[17]-0.223606797749979*uMax[7]*fr[16]+0.223606797749979*uMax[7]*fl[16]-0.159719141249985*uMax[6]*fr[14]-0.2500000000000001*uMax[2]*fr[14]+0.159719141249985*uMax[6]*fl[14]+0.2500000000000001*uMax[2]*fl[14]-0.5590169943749475*uMax[6]*fr[13]+0.5590169943749475*uMax[6]*fl[13]-0.276641667586244*u_r[4]*fr[12]+0.276641667586244*uMax[4]*fr[12]-0.4330127018922194*u_r[0]*fr[12]+0.4330127018922194*uMax[0]*fr[12]+0.276641667586244*u_l[4]*fl[12]+0.276641667586244*uMax[4]*fl[12]+0.4330127018922194*u_l[0]*fl[12]+0.4330127018922194*uMax[0]*fl[12]+0.5000000000000001*u_r[1]*fr[11]-0.5000000000000001*uMax[1]*fr[11]+0.5000000000000001*u_l[1]*fl[11]+0.5000000000000001*uMax[1]*fl[11]+0.3872983346207416*uMax[3]*fr[10]+0.3872983346207416*uMax[3]*fl[10]+0.159719141249985*u_r[4]*fr[8]-0.159719141249985*uMax[4]*fr[8]+0.25*u_r[0]*fr[8]-0.25*uMax[0]*fr[8]+0.159719141249985*u_l[4]*fl[8]+0.159719141249985*uMax[4]*fl[8]+0.25*u_l[0]*fl[8]+0.25*uMax[0]*fl[8]+0.5590169943749475*u_r[4]*fr[7]-0.5590169943749475*uMax[4]*fr[7]+0.5590169943749475*u_l[4]*fl[7]+0.5590169943749475*uMax[4]*fl[7]+0.4330127018922194*fr[5]*uMax[6]+0.4330127018922194*fl[5]*uMax[6]-0.2500000000000001*fr[3]*uMax[6]+0.2500000000000001*fl[3]*uMax[6]-0.223606797749979*uMax[3]*fr[6]+0.223606797749979*uMax[3]*fl[6]-0.4330127018922193*fr[1]*u_r[4]+0.25*fr[0]*u_r[4]+0.4330127018922193*fl[1]*u_l[4]+0.25*fl[0]*u_l[4]+0.4330127018922193*fr[1]*uMax[4]+0.4330127018922193*fl[1]*uMax[4]-0.25*fr[0]*uMax[4]+0.25*fl[0]*uMax[4]-0.3872983346207416*u_r[1]*fr[4]+0.3872983346207416*uMax[1]*fr[4]+0.3872983346207416*u_l[1]*fl[4]+0.3872983346207416*uMax[1]*fl[4]+0.223606797749979*u_r[1]*fr[2]-0.223606797749979*uMax[1]*fr[2]+0.223606797749979*u_l[1]*fl[2]+0.223606797749979*uMax[1]*fl[2])*mass*volFact; 
  heat_flux[0] += ((-0.4330127018922194*q_r[7]*fr[19])+0.4330127018922194*qMax[7]*fr[19]+0.4330127018922194*q_l[7]*fl[19]+0.4330127018922194*qMax[7]*fl[19]-0.4330127018922194*q_r[6]*fr[18]+0.4330127018922194*qMax[6]*fr[18]+0.4330127018922194*q_l[6]*fl[18]+0.4330127018922194*qMax[6]*fl[18]+0.5590169943749475*q_r[3]*fr[17]-0.5590169943749475*qMax[3]*fr[17]+0.5590169943749475*q_l[3]*fl[17]+0.5590169943749475*qMax[3]*fl[17]+0.25*q_r[7]*fr[16]-0.25*qMax[7]*fr[16]+0.25*q_l[7]*fl[16]+0.25*qMax[7]*fl[16]-0.4330127018922194*q_r[5]*fr[15]+0.4330127018922194*qMax[5]*fr[15]+0.4330127018922194*q_l[5]*fl[15]+0.4330127018922194*qMax[5]*fl[15]+0.25*q_r[6]*fr[14]-0.25*qMax[6]*fr[14]+0.25*q_l[6]*fl[14]+0.25*qMax[6]*fl[14]+0.5590169943749476*q_r[2]*fr[13]-0.5590169943749476*qMax[2]*fr[13]+0.5590169943749476*q_l[2]*fl[13]+0.5590169943749476*qMax[2]*fl[13]-0.4330127018922194*q_r[4]*fr[12]+0.4330127018922194*qMax[4]*fr[12]+0.4330127018922194*q_l[4]*fl[12]+0.4330127018922194*qMax[4]*fl[12]+0.5590169943749476*q_r[1]*fr[11]-0.5590169943749476*qMax[1]*fr[11]+0.5590169943749476*q_l[1]*fl[11]+0.5590169943749476*qMax[1]*fl[11]-0.4330127018922193*q_r[3]*fr[10]+0.4330127018922193*qMax[3]*fr[10]+0.4330127018922193*q_l[3]*fl[10]+0.4330127018922193*qMax[3]*fl[10]+0.25*q_r[5]*fr[9]-0.25*qMax[5]*fr[9]+0.25*q_l[5]*fl[9]+0.25*qMax[5]*fl[9]+0.25*q_r[4]*fr[8]-0.25*qMax[4]*fr[8]+0.25*q_l[4]*fl[8]+0.25*qMax[4]*fl[8]+0.5590169943749475*q_r[0]*fr[7]-0.5590169943749475*qMax[0]*fr[7]+0.5590169943749475*q_l[0]*fl[7]+0.5590169943749475*qMax[0]*fl[7]+0.25*q_r[3]*fr[6]-0.25*qMax[3]*fr[6]+0.25*q_l[3]*fl[6]+0.25*qMax[3]*fl[6]-0.4330127018922193*q_r[2]*fr[5]+0.4330127018922193*qMax[2]*fr[5]+0.4330127018922193*q_l[2]*fl[5]+0.4330127018922193*qMax[2]*fl[5]-0.4330127018922193*q_r[1]*fr[4]+0.4330127018922193*qMax[1]*fr[4]+0.4330127018922193*q_l[1]*fl[4]+0.4330127018922193*qMax[1]*fl[4]+0.25*q_r[2]*fr[3]-0.25*qMax[2]*fr[3]+0.25*q_l[2]*fl[3]+0.25*qMax[2]*fl[3]+0.25*q_r[1]*fr[2]-0.25*qMax[1]*fr[2]+0.25*q_l[1]*fl[2]+0.25*qMax[1]*fl[2]-0.4330127018922193*q_r[0]*fr[1]+0.4330127018922193*qMax[0]*fr[1]+0.4330127018922193*q_l[0]*fl[1]+0.4330127018922193*qMax[0]*fl[1]+0.25*fr[0]*q_r[0]+0.25*fl[0]*q_l[0]-0.25*fr[0]*qMax[0]+0.25*fl[0]*qMax[0])*mass*volFact; 
  heat_flux[1] += ((-0.4330127018922193*q_r[5]*fr[19])+0.4330127018922193*qMax[5]*fr[19]+0.4330127018922193*q_l[5]*fl[19]+0.4330127018922193*qMax[5]*fl[19]-0.3872983346207416*q_r[3]*fr[18]+0.3872983346207416*qMax[3]*fr[18]+0.3872983346207416*q_l[3]*fl[18]+0.3872983346207416*qMax[3]*fl[18]+0.5000000000000001*q_r[6]*fr[17]-0.5000000000000001*qMax[6]*fr[17]+0.5590169943749475*q_r[2]*fr[17]-0.5590169943749475*qMax[2]*fr[17]+0.5000000000000001*q_l[6]*fl[17]+0.5000000000000001*qMax[6]*fl[17]+0.5590169943749475*q_l[2]*fl[17]+0.5590169943749475*qMax[2]*fl[17]+0.2500000000000001*q_r[5]*fr[16]-0.2500000000000001*qMax[5]*fr[16]+0.2500000000000001*q_l[5]*fl[16]+0.2500000000000001*qMax[5]*fl[16]-0.4330127018922193*q_r[7]*fr[15]+0.4330127018922193*qMax[7]*fr[15]+0.4330127018922193*q_l[7]*fl[15]+0.4330127018922193*qMax[7]*fl[15]+0.223606797749979*q_r[3]*fr[14]-0.223606797749979*qMax[3]*fr[14]+0.223606797749979*q_l[3]*fl[14]+0.223606797749979*qMax[3]*fl[14]+0.5590169943749476*q_r[3]*fr[13]-0.5590169943749476*qMax[3]*fr[13]+0.5590169943749476*q_l[3]*fl[13]+0.5590169943749476*qMax[3]*fl[13]-0.3872983346207417*q_r[1]*fr[12]+0.3872983346207417*qMax[1]*fr[12]+0.3872983346207417*q_l[1]*fl[12]+0.3872983346207417*qMax[1]*fl[12]+0.5000000000000001*q_r[4]*fr[11]-0.5000000000000001*qMax[4]*fr[11]+0.5590169943749476*q_r[0]*fr[11]-0.5590169943749476*qMax[0]*fr[11]+0.5000000000000001*q_l[4]*fl[11]+0.5000000000000001*qMax[4]*fl[11]+0.5590169943749476*q_l[0]*fl[11]+0.5590169943749476*qMax[0]*fl[11]-0.3872983346207417*q_r[6]*fr[10]+0.3872983346207417*qMax[6]*fr[10]-0.4330127018922193*q_r[2]*fr[10]+0.4330127018922193*qMax[2]*fr[10]+0.3872983346207417*q_l[6]*fl[10]+0.3872983346207417*qMax[6]*fl[10]+0.4330127018922193*q_l[2]*fl[10]+0.4330127018922193*qMax[2]*fl[10]+0.2500000000000001*q_r[7]*fr[9]-0.2500000000000001*qMax[7]*fr[9]+0.2500000000000001*q_l[7]*fl[9]+0.2500000000000001*qMax[7]*fl[9]+0.223606797749979*q_r[1]*fr[8]-0.223606797749979*qMax[1]*fr[8]+0.223606797749979*q_l[1]*fl[8]+0.223606797749979*qMax[1]*fl[8]+0.5590169943749475*q_r[1]*fr[7]-0.5590169943749475*qMax[1]*fr[7]+0.5590169943749475*q_l[1]*fl[7]+0.5590169943749475*qMax[1]*fl[7]+0.223606797749979*fr[6]*q_r[6]+0.223606797749979*fl[6]*q_l[6]-0.223606797749979*fr[6]*qMax[6]+0.223606797749979*fl[6]*qMax[6]+0.25*q_r[2]*fr[6]-0.25*qMax[2]*fr[6]+0.25*q_l[2]*fl[6]+0.25*qMax[2]*fl[6]-0.4330127018922193*q_r[3]*fr[5]+0.4330127018922193*qMax[3]*fr[5]+0.4330127018922193*q_l[3]*fl[5]+0.4330127018922193*qMax[3]*fl[5]-0.3872983346207416*fr[4]*q_r[4]+0.223606797749979*fr[2]*q_r[4]+0.3872983346207416*fl[4]*q_l[4]+0.223606797749979*fl[2]*q_l[4]+0.3872983346207416*fr[4]*qMax[4]+0.3872983346207416*fl[4]*qMax[4]-0.223606797749979*fr[2]*qMax[4]+0.223606797749979*fl[2]*qMax[4]-0.4330127018922193*q_r[0]*fr[4]+0.4330127018922193*qMax[0]*fr[4]+0.4330127018922193*q_l[0]*fl[4]+0.4330127018922193*qMax[0]*fl[4]+0.25*fr[3]*q_r[3]+0.25*fl[3]*q_l[3]-0.25*fr[3]*qMax[3]+0.25*fl[3]*qMax[3]+0.25*q_r[0]*fr[2]-0.25*qMax[0]*fr[2]+0.25*q_l[0]*fl[2]+0.25*qMax[0]*fl[2]-0.4330127018922193*fr[1]*q_r[1]+0.25*fr[0]*q_r[1]+0.4330127018922193*fl[1]*q_l[1]+0.25*fl[0]*q_l[1]+0.4330127018922193*fr[1]*qMax[1]+0.4330127018922193*fl[1]*qMax[1]-0.25*fr[0]*qMax[1]+0.25*fl[0]*qMax[1])*mass*volFact; 
  heat_flux[2] += ((-0.3872983346207417*q_r[7]*fr[19])+0.3872983346207417*qMax[7]*fr[19]+0.3872983346207417*q_l[7]*fl[19]+0.3872983346207417*qMax[7]*fl[19]-0.276641667586244*q_r[6]*fr[18]+0.276641667586244*qMax[6]*fr[18]-0.4330127018922193*q_r[2]*fr[18]+0.4330127018922193*qMax[2]*fr[18]+0.276641667586244*q_l[6]*fl[18]+0.276641667586244*qMax[6]*fl[18]+0.4330127018922193*q_l[2]*fl[18]+0.4330127018922193*qMax[2]*fl[18]+0.5*q_r[3]*fr[17]-0.5*qMax[3]*fr[17]+0.5*q_l[3]*fl[17]+0.5*qMax[3]*fl[17]+0.223606797749979*q_r[7]*fr[16]-0.223606797749979*qMax[7]*fr[16]+0.223606797749979*q_l[7]*fl[16]+0.223606797749979*qMax[7]*fl[16]+0.159719141249985*q_r[6]*fr[14]-0.159719141249985*qMax[6]*fr[14]+0.2500000000000001*q_r[2]*fr[14]-0.2500000000000001*qMax[2]*fr[14]+0.159719141249985*q_l[6]*fl[14]+0.159719141249985*qMax[6]*fl[14]+0.2500000000000001*q_l[2]*fl[14]+0.2500000000000001*qMax[2]*fl[14]+0.5590169943749475*q_r[6]*fr[13]-0.5590169943749475*qMax[6]*fr[13]+0.5590169943749475*q_l[6]*fl[13]+0.5590169943749475*qMax[6]*fl[13]-0.276641667586244*q_r[4]*fr[12]+0.276641667586244*qMax[4]*fr[12]-0.4330127018922194*q_r[0]*fr[12]+0.4330127018922194*qMax[0]*fr[12]+0.276641667586244*q_l[4]*fl[12]+0.276641667586244*qMax[4]*fl[12]+0.4330127018922194*q_l[0]*fl[12]+0.4330127018922194*qMax[0]*fl[12]+0.5000000000000001*q_r[1]*fr[11]-0.5000000000000001*qMax[1]*fr[11]+0.5000000000000001*q_l[1]*fl[11]+0.5000000000000001*qMax[1]*fl[11]-0.3872983346207416*q_r[3]*fr[10]+0.3872983346207416*qMax[3]*fr[10]+0.3872983346207416*q_l[3]*fl[10]+0.3872983346207416*qMax[3]*fl[10]+0.159719141249985*q_r[4]*fr[8]-0.159719141249985*qMax[4]*fr[8]+0.25*q_r[0]*fr[8]-0.25*qMax[0]*fr[8]+0.159719141249985*q_l[4]*fl[8]+0.159719141249985*qMax[4]*fl[8]+0.25*q_l[0]*fl[8]+0.25*qMax[0]*fl[8]+0.5590169943749475*q_r[4]*fr[7]-0.5590169943749475*qMax[4]*fr[7]+0.5590169943749475*q_l[4]*fl[7]+0.5590169943749475*qMax[4]*fl[7]-0.4330127018922194*fr[5]*q_r[6]+0.2500000000000001*fr[3]*q_r[6]+0.4330127018922194*fl[5]*q_l[6]+0.2500000000000001*fl[3]*q_l[6]+0.4330127018922194*fr[5]*qMax[6]+0.4330127018922194*fl[5]*qMax[6]-0.2500000000000001*fr[3]*qMax[6]+0.2500000000000001*fl[3]*qMax[6]+0.223606797749979*q_r[3]*fr[6]-0.223606797749979*qMax[3]*fr[6]+0.223606797749979*q_l[3]*fl[6]+0.223606797749979*qMax[3]*fl[6]-0.4330127018922193*fr[1]*q_r[4]+0.25*fr[0]*q_r[4]+0.4330127018922193*fl[1]*q_l[4]+0.25*fl[0]*q_l[4]+0.4330127018922193*fr[1]*qMax[4]+0.4330127018922193*fl[1]*qMax[4]-0.25*fr[0]*qMax[4]+0.25*fl[0]*qMax[4]-0.3872983346207416*q_r[1]*fr[4]+0.3872983346207416*qMax[1]*fr[4]+0.3872983346207416*q_l[1]*fl[4]+0.3872983346207416*qMax[1]*fl[4]+0.223606797749979*q_r[1]*fr[2]-0.223606797749979*qMax[1]*fr[2]+0.223606797749979*q_l[1]*fl[2]+0.223606797749979*qMax[1]*fl[2])*mass*volFact; 
} 
