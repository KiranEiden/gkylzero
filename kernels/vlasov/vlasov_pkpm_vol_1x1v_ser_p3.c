#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_pkpm_vol_1x1v_ser_p3(const double *w, const double *dxv, const double *bvar, const double *u_i, const double *bb_grad_u, const double *p_force, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // bvar:      magnetic field unit vector (nine components; first three components, b_i, other six components, b_i b_j.) 
  // u_i:       flow velocity  // p_force:   total pressure force = 1/rho (b . div(P) + p_perp div(b)) for Euler PKPM.
  // bb_grad_u: bb : grad(u).
  // f:         Input distribution function.
  // out:       Incremented output.
  const double dx0 = 2.0/dxv[0]; 
  const double dv1par = 2.0/dxv[1]; 
  const double dvpar = dxv[1], wvpar = w[1]; 
  const double *ux = &u_i[0]; 
  const double *uy = &u_i[4]; 
  const double *uz = &u_i[8]; 
  const double *bx = &bvar[0]; 
  const double *by = &bvar[4]; 
  const double *bz = &bvar[8]; 

  double cflFreq_mid = 0.0; 
  double alpha_cdim[12] = {0.0}; 
  double alpha_vdim[12] = {0.0}; 

  alpha_cdim[0] = 1.414213562373095*dx0*(bx[0]*wvpar+ux[0]); 
  alpha_cdim[1] = 1.414213562373095*dx0*(bx[1]*wvpar+ux[1]); 
  alpha_cdim[2] = 0.408248290463863*bx[0]*dvpar*dx0; 
  alpha_cdim[3] = 0.408248290463863*bx[1]*dvpar*dx0; 
  alpha_cdim[4] = 1.414213562373095*dx0*(bx[2]*wvpar+ux[2]); 
  alpha_cdim[6] = 0.408248290463863*bx[2]*dvpar*dx0; 
  alpha_cdim[8] = 1.414213562373095*dx0*(bx[3]*wvpar+ux[3]); 
  alpha_cdim[10] = 0.408248290463863*bx[3]*dvpar*dx0; 
  cflFreq_mid += 7.0*fabs(0.25*alpha_cdim[0]-0.2795084971874737*alpha_cdim[4]); 

  alpha_vdim[0] = 1.414213562373095*p_force[0]*dv1par-1.414213562373095*bb_grad_u[0]*dv1par*wvpar; 
  alpha_vdim[1] = 1.414213562373095*p_force[1]*dv1par-1.414213562373095*bb_grad_u[1]*dv1par*wvpar; 
  alpha_vdim[2] = -0.408248290463863*bb_grad_u[0]*dv1par*dvpar; 
  alpha_vdim[3] = -0.408248290463863*bb_grad_u[1]*dv1par*dvpar; 
  alpha_vdim[4] = 1.414213562373095*p_force[2]*dv1par-1.414213562373095*bb_grad_u[2]*dv1par*wvpar; 
  alpha_vdim[6] = -0.408248290463863*bb_grad_u[2]*dv1par*dvpar; 
  alpha_vdim[8] = 1.414213562373095*p_force[3]*dv1par-1.414213562373095*bb_grad_u[3]*dv1par*wvpar; 
  alpha_vdim[10] = -0.4082482904638629*bb_grad_u[3]*dv1par*dvpar; 
  cflFreq_mid += 7.0*fabs(0.25*alpha_vdim[0]-0.2795084971874737*alpha_vdim[4]); 

  out[1] += 0.8660254037844386*(alpha_cdim[10]*f[10]+alpha_cdim[8]*f[8]+alpha_cdim[6]*f[6]+alpha_cdim[4]*f[4]+alpha_cdim[3]*f[3]+alpha_cdim[2]*f[2]+alpha_cdim[1]*f[1]+alpha_cdim[0]*f[0]); 
  out[2] += 0.8660254037844386*(alpha_vdim[10]*f[10]+alpha_vdim[8]*f[8]+alpha_vdim[6]*f[6]+alpha_vdim[4]*f[4]+alpha_vdim[3]*f[3]+alpha_vdim[2]*f[2]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.8660254037844386*alpha_cdim[8]*f[10]+0.7606388292556648*(alpha_vdim[6]*f[10]+f[6]*alpha_vdim[10])+0.8660254037844386*f[8]*alpha_cdim[10]+0.7606388292556648*(alpha_vdim[4]*f[8]+f[4]*alpha_vdim[8])+0.7745966692414833*alpha_cdim[3]*f[7]+0.8660254037844386*alpha_cdim[4]*f[6]+0.7745966692414833*(alpha_vdim[3]*f[6]+f[3]*alpha_vdim[6])+0.8660254037844386*f[4]*alpha_cdim[6]+0.7745966692414833*(alpha_cdim[2]*f[5]+alpha_vdim[1]*f[4]+f[1]*alpha_vdim[4])+0.8660254037844386*((alpha_vdim[2]+alpha_cdim[1])*f[3]+f[2]*alpha_vdim[3]+f[1]*alpha_cdim[3]+alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[4] += 1.700840128541522*(alpha_cdim[6]*f[10]+f[6]*alpha_cdim[10]+alpha_cdim[4]*f[8]+f[4]*alpha_cdim[8])+1.732050807568877*(alpha_cdim[3]*f[6]+f[3]*alpha_cdim[6]+alpha_cdim[1]*f[4]+f[1]*alpha_cdim[4])+1.936491673103709*(alpha_cdim[2]*f[3]+f[2]*alpha_cdim[3]+alpha_cdim[0]*f[1]+f[0]*alpha_cdim[1]); 
  out[5] += 1.936491673103709*(alpha_vdim[8]*f[10]+f[8]*alpha_vdim[10])+1.732050807568877*alpha_vdim[3]*f[7]+1.936491673103709*(alpha_vdim[4]*f[6]+f[4]*alpha_vdim[6])+1.732050807568877*alpha_vdim[2]*f[5]+1.936491673103709*(alpha_vdim[1]*f[3]+f[1]*alpha_vdim[3]+alpha_vdim[0]*f[2]+f[0]*alpha_vdim[2]); 
  out[6] += (0.5163977794943223*alpha_vdim[10]+1.700840128541522*alpha_cdim[4])*f[10]+0.7606388292556648*(alpha_vdim[3]*f[10]+f[3]*alpha_vdim[10])+1.700840128541522*f[4]*alpha_cdim[10]+(0.5163977794943223*alpha_vdim[8]+1.700840128541522*alpha_cdim[6])*f[8]+0.7606388292556648*(alpha_vdim[1]*f[8]+f[1]*alpha_vdim[8])+1.700840128541522*f[6]*alpha_cdim[8]+(1.549193338482967*alpha_cdim[6]+1.732050807568877*alpha_cdim[2])*f[7]+(0.5532833351724881*alpha_vdim[6]+0.8660254037844386*alpha_vdim[2]+1.732050807568877*alpha_cdim[1])*f[6]+0.8660254037844386*f[2]*alpha_vdim[6]+1.732050807568877*(f[1]*alpha_cdim[6]+alpha_cdim[3]*f[5])+(0.5532833351724881*alpha_vdim[4]+1.732050807568877*alpha_cdim[3])*f[4]+0.8660254037844386*(alpha_vdim[0]*f[4]+f[0]*alpha_vdim[4])+f[3]*(1.732050807568877*alpha_cdim[4]+0.7745966692414833*alpha_vdim[3])+1.936491673103709*(alpha_cdim[0]*f[3]+f[0]*alpha_cdim[3]+alpha_cdim[1]*f[2])+f[1]*(1.936491673103709*alpha_cdim[2]+0.7745966692414833*alpha_vdim[1]); 
  out[7] += 0.7606388292556648*alpha_cdim[3]*f[11]+0.7745966692414833*alpha_cdim[10]*f[10]+1.700840128541522*(alpha_vdim[4]*f[10]+f[4]*alpha_vdim[10])+0.7606388292556648*alpha_cdim[2]*f[9]+1.700840128541522*(alpha_vdim[6]*f[8]+f[6]*alpha_vdim[8])+(1.549193338482967*alpha_vdim[6]+1.732050807568877*alpha_vdim[2]+0.8660254037844386*alpha_cdim[1])*f[7]+0.7745966692414833*alpha_cdim[6]*f[6]+1.732050807568877*(alpha_vdim[1]*f[6]+f[1]*alpha_vdim[6])+(1.732050807568877*alpha_vdim[3]+0.8660254037844386*alpha_cdim[0])*f[5]+1.732050807568877*alpha_vdim[3]*f[4]+f[3]*(1.732050807568877*alpha_vdim[4]+0.7745966692414833*alpha_cdim[3])+1.936491673103709*(alpha_vdim[0]*f[3]+f[0]*alpha_vdim[3])+0.7745966692414833*alpha_cdim[2]*f[2]+1.936491673103709*(alpha_vdim[1]*f[2]+f[1]*alpha_vdim[2]); 
  out[8] += 3.086709862908689*alpha_cdim[10]*f[10]+2.598076211353316*(alpha_cdim[3]*f[10]+f[3]*alpha_cdim[10])+3.086709862908689*alpha_cdim[8]*f[8]+2.598076211353316*(alpha_cdim[1]*f[8]+f[1]*alpha_cdim[8])+3.212698020578431*alpha_cdim[6]*f[6]+2.958039891549809*(alpha_cdim[2]*f[6]+f[2]*alpha_cdim[6])+3.212698020578431*alpha_cdim[4]*f[4]+2.958039891549809*(alpha_cdim[0]*f[4]+f[0]*alpha_cdim[4])+3.968626966596886*alpha_cdim[3]*f[3]+1.322875655532295*alpha_cdim[2]*f[2]+3.968626966596886*alpha_cdim[1]*f[1]+1.322875655532295*alpha_cdim[0]*f[0]; 
  out[9] += 2.598076211353316*alpha_vdim[3]*f[11]+3.968626966596886*alpha_vdim[10]*f[10]+2.598076211353316*alpha_vdim[2]*f[9]+1.322875655532295*alpha_vdim[8]*f[8]+2.958039891549809*alpha_vdim[1]*f[7]+3.968626966596886*alpha_vdim[6]*f[6]+2.958039891549809*alpha_vdim[0]*f[5]+1.322875655532295*alpha_vdim[4]*f[4]+3.968626966596886*(alpha_vdim[3]*f[3]+alpha_vdim[2]*f[2])+1.322875655532295*(alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[10] += (3.086709862908689*alpha_cdim[8]+0.5163977794943223*alpha_vdim[6]+0.8660254037844386*alpha_vdim[2]+2.598076211353316*alpha_cdim[1])*f[10]+(0.5163977794943223*f[6]+0.8660254037844386*f[2])*alpha_vdim[10]+(3.086709862908689*f[8]+2.32379000772445*f[7]+2.598076211353316*f[1])*alpha_cdim[10]+(0.5163977794943223*alpha_vdim[4]+2.598076211353316*alpha_cdim[3]+0.8660254037844386*alpha_vdim[0])*f[8]+(0.5163977794943223*f[4]+0.8660254037844386*f[0])*alpha_vdim[8]+2.598076211353316*f[3]*alpha_cdim[8]+3.54964786985977*alpha_cdim[3]*f[7]+(3.212698020578431*alpha_cdim[4]+0.7606388292556648*alpha_vdim[3]+2.958039891549809*alpha_cdim[0])*f[6]+0.7606388292556648*f[3]*alpha_vdim[6]+(2.645751311064591*f[5]+3.212698020578431*f[4]+2.958039891549809*f[0])*alpha_cdim[6]+alpha_cdim[2]*(1.183215956619923*f[5]+2.958039891549809*f[4])+0.7606388292556648*(alpha_vdim[1]*f[4]+f[1]*alpha_vdim[4])+2.958039891549809*f[2]*alpha_cdim[4]+3.968626966596886*(alpha_cdim[1]*f[3]+f[1]*alpha_cdim[3])+1.322875655532295*(alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]); 
  out[11] += (2.32379000772445*alpha_vdim[6]+2.598076211353316*alpha_vdim[2]+0.8660254037844386*alpha_cdim[1])*f[11]+3.485685011586674*(alpha_vdim[6]*f[10]+f[6]*alpha_vdim[10])+(2.598076211353316*alpha_vdim[3]+0.8660254037844386*alpha_cdim[0])*f[9]+1.161895003862225*(alpha_vdim[4]*f[8]+f[4]*alpha_vdim[8])+(2.645751311064591*alpha_vdim[4]+0.7606388292556648*alpha_cdim[3]+2.958039891549809*alpha_vdim[0])*f[7]+3.54964786985977*(alpha_vdim[3]*f[6]+f[3]*alpha_vdim[6])+(0.7606388292556648*alpha_cdim[2]+2.958039891549809*alpha_vdim[1])*f[5]+1.183215956619923*(alpha_vdim[1]*f[4]+f[1]*alpha_vdim[4])+3.968626966596886*(alpha_vdim[2]*f[3]+f[2]*alpha_vdim[3])+1.322875655532295*(alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 

  return cflFreq_mid; 
} 
