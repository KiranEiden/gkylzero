#include <gkyl_isoeuler_kernels.h> 
GKYL_CU_DH double isoeuler_vol_2x_ser_p2(const double *w, const double *dxv, const double vth, const double *uvar, const double *statevec, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // statevec: [rho, rho ux, rho uy, rho uz].
  // uvar: [ux, uy, uz].
  // out: Incremented output.

  const double *rho = &statevec[0]; 
  const double *rhou0 = &statevec[8]; 
  const double *rhou1 = &statevec[16]; 
  const double *rhou2 = &statevec[24]; 
  const double *uvar0 = &uvar[0]; 
  const double *uvar1 = &uvar[8]; 
  const double *uvar2 = &uvar[16]; 
  double *outrho = &out[0]; 
  double *outrhou0 = &out[8]; 
  double *outrhou1 = &out[16]; 
  double *outrhou2 = &out[24]; 
  double dx10 = 2./dxv[0]; 
  double dx11 = 2./dxv[1]; 

  double vthsq = vth*vth; 
  double alpha_mid = 0.0; 
  alpha_mid += 0.5*dx10*(fabs(0.5*uvar0[0]-0.5590169943749475*(uvar0[5]+uvar0[4]))+sqrt((5*vthsq)/(3*(0.5*rho[0]-0.5590169943749475*(rho[5]+rho[4]))))); 
  alpha_mid += 0.5*dx11*(fabs(0.5*uvar1[0]-0.5590169943749475*(uvar1[5]+uvar1[4]))+sqrt((5*vthsq)/(3*(0.5*rho[0]-0.5590169943749475*(rho[5]+rho[4]))))); 

  outrho[1] += 0.8660254037844386*rho[7]*rhou0[7]*dx10+0.8660254037844386*rho[6]*rhou0[6]*dx10+0.8660254037844386*rho[5]*rhou0[5]*dx10+0.8660254037844386*rho[4]*rhou0[4]*dx10+0.8660254037844386*rho[3]*rhou0[3]*dx10+0.8660254037844386*rho[2]*rhou0[2]*dx10+0.8660254037844386*rho[1]*rhou0[1]*dx10+0.8660254037844386*rho[0]*rhou0[0]*dx10; 
  outrho[2] += 0.8660254037844386*rho[7]*rhou1[7]*dx11+0.8660254037844386*rho[6]*rhou1[6]*dx11+0.8660254037844386*rho[5]*rhou1[5]*dx11+0.8660254037844386*rho[4]*rhou1[4]*dx11+0.8660254037844386*rho[3]*rhou1[3]*dx11+0.8660254037844386*rho[2]*rhou1[2]*dx11+0.8660254037844386*rho[1]*rhou1[1]*dx11+0.8660254037844386*rho[0]*rhou1[0]*dx11; 
  outrho[3] += 0.8660254037844387*rho[5]*rhou1[7]*dx11+0.8660254037844387*rhou1[5]*rho[7]*dx11+0.7745966692414834*rho[3]*rhou1[6]*dx11+0.7745966692414834*rhou1[3]*rho[6]*dx11+0.7745966692414833*rho[1]*rhou1[4]*dx11+0.7745966692414833*rhou1[1]*rho[4]*dx11+0.8660254037844386*rho[2]*rhou1[3]*dx11+0.8660254037844386*rhou1[2]*rho[3]*dx11+0.8660254037844386*rho[0]*rhou1[1]*dx11+0.8660254037844386*rhou1[0]*rho[1]*dx11+0.7745966692414834*rho[3]*rhou0[7]*dx10+0.7745966692414834*rhou0[3]*rho[7]*dx10+0.8660254037844387*rho[4]*rhou0[6]*dx10+0.8660254037844387*rhou0[4]*rho[6]*dx10+0.7745966692414833*rho[2]*rhou0[5]*dx10+0.7745966692414833*rhou0[2]*rho[5]*dx10+0.8660254037844386*rho[1]*rhou0[3]*dx10+0.8660254037844386*rhou0[1]*rho[3]*dx10+0.8660254037844386*rho[0]*rhou0[2]*dx10+0.8660254037844386*rhou0[0]*rho[2]*dx10; 
  outrho[4] += 1.936491673103709*rho[5]*rhou0[7]*dx10+1.936491673103709*rhou0[5]*rho[7]*dx10+1.732050807568877*rho[3]*rhou0[6]*dx10+1.732050807568877*rhou0[3]*rho[6]*dx10+1.732050807568877*rho[1]*rhou0[4]*dx10+1.732050807568877*rhou0[1]*rho[4]*dx10+1.936491673103709*rho[2]*rhou0[3]*dx10+1.936491673103709*rhou0[2]*rho[3]*dx10+1.936491673103709*rho[0]*rhou0[1]*dx10+1.936491673103709*rhou0[0]*rho[1]*dx10; 
  outrho[5] += 1.732050807568877*rho[3]*rhou1[7]*dx11+1.732050807568877*rhou1[3]*rho[7]*dx11+1.936491673103709*rho[4]*rhou1[6]*dx11+1.936491673103709*rhou1[4]*rho[6]*dx11+1.732050807568877*rho[2]*rhou1[5]*dx11+1.732050807568877*rhou1[2]*rho[5]*dx11+1.936491673103709*rho[1]*rhou1[3]*dx11+1.936491673103709*rhou1[1]*rho[3]*dx11+1.936491673103709*rho[0]*rhou1[2]*dx11+1.936491673103709*rhou1[0]*rho[2]*dx11; 
  outrho[6] += 0.7745966692414834*rho[7]*rhou1[7]*dx11+0.5532833351724881*rho[6]*rhou1[6]*dx11+0.8660254037844386*rho[2]*rhou1[6]*dx11+0.8660254037844386*rhou1[2]*rho[6]*dx11+0.5532833351724881*rho[4]*rhou1[4]*dx11+0.8660254037844387*rho[0]*rhou1[4]*dx11+0.8660254037844387*rhou1[0]*rho[4]*dx11+0.7745966692414834*rho[3]*rhou1[3]*dx11+0.7745966692414834*rho[1]*rhou1[1]*dx11+1.549193338482967*rho[6]*rhou0[7]*dx10+1.732050807568877*rho[2]*rhou0[7]*dx10+1.549193338482967*rhou0[6]*rho[7]*dx10+1.732050807568877*rhou0[2]*rho[7]*dx10+1.732050807568877*rho[1]*rhou0[6]*dx10+1.732050807568877*rhou0[1]*rho[6]*dx10+1.732050807568877*rho[3]*rhou0[5]*dx10+1.732050807568877*rhou0[3]*rho[5]*dx10+1.732050807568877*rho[3]*rhou0[4]*dx10+1.732050807568877*rhou0[3]*rho[4]*dx10+1.936491673103709*rho[0]*rhou0[3]*dx10+1.936491673103709*rhou0[0]*rho[3]*dx10+1.936491673103709*rho[1]*rhou0[2]*dx10+1.936491673103709*rhou0[1]*rho[2]*dx10; 
  outrho[7] += 1.549193338482967*rho[6]*rhou1[7]*dx11+1.732050807568877*rho[2]*rhou1[7]*dx11+1.549193338482967*rhou1[6]*rho[7]*dx11+1.732050807568877*rhou1[2]*rho[7]*dx11+1.732050807568877*rho[1]*rhou1[6]*dx11+1.732050807568877*rhou1[1]*rho[6]*dx11+1.732050807568877*rho[3]*rhou1[5]*dx11+1.732050807568877*rhou1[3]*rho[5]*dx11+1.732050807568877*rho[3]*rhou1[4]*dx11+1.732050807568877*rhou1[3]*rho[4]*dx11+1.936491673103709*rho[0]*rhou1[3]*dx11+1.936491673103709*rhou1[0]*rho[3]*dx11+1.936491673103709*rho[1]*rhou1[2]*dx11+1.936491673103709*rhou1[1]*rho[2]*dx11+0.5532833351724881*rho[7]*rhou0[7]*dx10+0.8660254037844386*rho[1]*rhou0[7]*dx10+0.8660254037844386*rhou0[1]*rho[7]*dx10+0.7745966692414834*rho[6]*rhou0[6]*dx10+0.5532833351724881*rho[5]*rhou0[5]*dx10+0.8660254037844387*rho[0]*rhou0[5]*dx10+0.8660254037844387*rhou0[0]*rho[5]*dx10+0.7745966692414834*rho[3]*rhou0[3]*dx10+0.7745966692414834*rho[2]*rhou0[2]*dx10; 

  outrhou0[1] += 3.0*rho[1]*dx10*vthsq+0.8660254037844386*rhou0[7]*uvar0[7]*dx10+0.8660254037844386*rhou0[6]*uvar0[6]*dx10+0.8660254037844386*rhou0[5]*uvar0[5]*dx10+0.8660254037844386*rhou0[4]*uvar0[4]*dx10+0.8660254037844386*rhou0[3]*uvar0[3]*dx10+0.8660254037844386*rhou0[2]*uvar0[2]*dx10+0.8660254037844386*rhou0[1]*uvar0[1]*dx10+0.8660254037844386*rhou0[0]*uvar0[0]*dx10; 
  outrhou0[2] += 0.8660254037844386*rhou0[7]*uvar1[7]*dx11+0.8660254037844386*rhou0[6]*uvar1[6]*dx11+0.8660254037844386*rhou0[5]*uvar1[5]*dx11+0.8660254037844386*rhou0[4]*uvar1[4]*dx11+0.8660254037844386*rhou0[3]*uvar1[3]*dx11+0.8660254037844386*rhou0[2]*uvar1[2]*dx11+0.8660254037844386*rhou0[1]*uvar1[1]*dx11+0.8660254037844386*rhou0[0]*uvar1[0]*dx11; 
  outrhou0[3] += 3.0*rho[3]*dx10*vthsq+0.8660254037844387*rhou0[5]*uvar1[7]*dx11+0.8660254037844387*uvar1[5]*rhou0[7]*dx11+0.7745966692414834*rhou0[3]*uvar1[6]*dx11+0.7745966692414834*uvar1[3]*rhou0[6]*dx11+0.7745966692414833*rhou0[1]*uvar1[4]*dx11+0.7745966692414833*uvar1[1]*rhou0[4]*dx11+0.8660254037844386*rhou0[2]*uvar1[3]*dx11+0.8660254037844386*uvar1[2]*rhou0[3]*dx11+0.8660254037844386*rhou0[0]*uvar1[1]*dx11+0.8660254037844386*uvar1[0]*rhou0[1]*dx11+0.7745966692414834*rhou0[3]*uvar0[7]*dx10+0.7745966692414834*uvar0[3]*rhou0[7]*dx10+0.8660254037844387*rhou0[4]*uvar0[6]*dx10+0.8660254037844387*uvar0[4]*rhou0[6]*dx10+0.7745966692414833*rhou0[2]*uvar0[5]*dx10+0.7745966692414833*uvar0[2]*rhou0[5]*dx10+0.8660254037844386*rhou0[1]*uvar0[3]*dx10+0.8660254037844386*uvar0[1]*rhou0[3]*dx10+0.8660254037844386*rhou0[0]*uvar0[2]*dx10+0.8660254037844386*uvar0[0]*rhou0[2]*dx10; 
  outrhou0[4] += 15.0*rho[4]*dx10*vthsq+1.936491673103709*rhou0[5]*uvar0[7]*dx10+1.936491673103709*uvar0[5]*rhou0[7]*dx10+1.732050807568877*rhou0[3]*uvar0[6]*dx10+1.732050807568877*uvar0[3]*rhou0[6]*dx10+1.732050807568877*rhou0[1]*uvar0[4]*dx10+1.732050807568877*uvar0[1]*rhou0[4]*dx10+1.936491673103709*rhou0[2]*uvar0[3]*dx10+1.936491673103709*uvar0[2]*rhou0[3]*dx10+1.936491673103709*rhou0[0]*uvar0[1]*dx10+1.936491673103709*uvar0[0]*rhou0[1]*dx10; 
  outrhou0[5] += 1.732050807568877*rhou0[3]*uvar1[7]*dx11+1.732050807568877*uvar1[3]*rhou0[7]*dx11+1.936491673103709*rhou0[4]*uvar1[6]*dx11+1.936491673103709*uvar1[4]*rhou0[6]*dx11+1.732050807568877*rhou0[2]*uvar1[5]*dx11+1.732050807568877*uvar1[2]*rhou0[5]*dx11+1.936491673103709*rhou0[1]*uvar1[3]*dx11+1.936491673103709*uvar1[1]*rhou0[3]*dx11+1.936491673103709*rhou0[0]*uvar1[2]*dx11+1.936491673103709*uvar1[0]*rhou0[2]*dx11; 
  outrhou0[6] += 15.0*rho[6]*dx10*vthsq+0.7745966692414834*rhou0[7]*uvar1[7]*dx11+0.5532833351724881*rhou0[6]*uvar1[6]*dx11+0.8660254037844386*rhou0[2]*uvar1[6]*dx11+0.8660254037844386*uvar1[2]*rhou0[6]*dx11+0.5532833351724881*rhou0[4]*uvar1[4]*dx11+0.8660254037844387*rhou0[0]*uvar1[4]*dx11+0.8660254037844387*uvar1[0]*rhou0[4]*dx11+0.7745966692414834*rhou0[3]*uvar1[3]*dx11+0.7745966692414834*rhou0[1]*uvar1[1]*dx11+1.549193338482967*rhou0[6]*uvar0[7]*dx10+1.732050807568877*rhou0[2]*uvar0[7]*dx10+1.549193338482967*uvar0[6]*rhou0[7]*dx10+1.732050807568877*uvar0[2]*rhou0[7]*dx10+1.732050807568877*rhou0[1]*uvar0[6]*dx10+1.732050807568877*uvar0[1]*rhou0[6]*dx10+1.732050807568877*rhou0[3]*uvar0[5]*dx10+1.732050807568877*uvar0[3]*rhou0[5]*dx10+1.732050807568877*rhou0[3]*uvar0[4]*dx10+1.732050807568877*uvar0[3]*rhou0[4]*dx10+1.936491673103709*rhou0[0]*uvar0[3]*dx10+1.936491673103709*uvar0[0]*rhou0[3]*dx10+1.936491673103709*rhou0[1]*uvar0[2]*dx10+1.936491673103709*uvar0[1]*rhou0[2]*dx10; 
  outrhou0[7] += 3.0*rho[7]*dx10*vthsq+1.549193338482967*rhou0[6]*uvar1[7]*dx11+1.732050807568877*rhou0[2]*uvar1[7]*dx11+1.549193338482967*uvar1[6]*rhou0[7]*dx11+1.732050807568877*uvar1[2]*rhou0[7]*dx11+1.732050807568877*rhou0[1]*uvar1[6]*dx11+1.732050807568877*uvar1[1]*rhou0[6]*dx11+1.732050807568877*rhou0[3]*uvar1[5]*dx11+1.732050807568877*uvar1[3]*rhou0[5]*dx11+1.732050807568877*rhou0[3]*uvar1[4]*dx11+1.732050807568877*uvar1[3]*rhou0[4]*dx11+1.936491673103709*rhou0[0]*uvar1[3]*dx11+1.936491673103709*uvar1[0]*rhou0[3]*dx11+1.936491673103709*rhou0[1]*uvar1[2]*dx11+1.936491673103709*uvar1[1]*rhou0[2]*dx11+0.5532833351724881*rhou0[7]*uvar0[7]*dx10+0.8660254037844386*rhou0[1]*uvar0[7]*dx10+0.8660254037844386*uvar0[1]*rhou0[7]*dx10+0.7745966692414834*rhou0[6]*uvar0[6]*dx10+0.5532833351724881*rhou0[5]*uvar0[5]*dx10+0.8660254037844387*rhou0[0]*uvar0[5]*dx10+0.8660254037844387*uvar0[0]*rhou0[5]*dx10+0.7745966692414834*rhou0[3]*uvar0[3]*dx10+0.7745966692414834*rhou0[2]*uvar0[2]*dx10; 

  outrhou1[1] += 0.8660254037844386*rhou1[7]*uvar0[7]*dx10+0.8660254037844386*rhou1[6]*uvar0[6]*dx10+0.8660254037844386*rhou1[5]*uvar0[5]*dx10+0.8660254037844386*rhou1[4]*uvar0[4]*dx10+0.8660254037844386*rhou1[3]*uvar0[3]*dx10+0.8660254037844386*rhou1[2]*uvar0[2]*dx10+0.8660254037844386*rhou1[1]*uvar0[1]*dx10+0.8660254037844386*rhou1[0]*uvar0[0]*dx10; 
  outrhou1[2] += 3.0*rho[2]*dx11*vthsq+0.8660254037844386*rhou1[7]*uvar1[7]*dx11+0.8660254037844386*rhou1[6]*uvar1[6]*dx11+0.8660254037844386*rhou1[5]*uvar1[5]*dx11+0.8660254037844386*rhou1[4]*uvar1[4]*dx11+0.8660254037844386*rhou1[3]*uvar1[3]*dx11+0.8660254037844386*rhou1[2]*uvar1[2]*dx11+0.8660254037844386*rhou1[1]*uvar1[1]*dx11+0.8660254037844386*rhou1[0]*uvar1[0]*dx11; 
  outrhou1[3] += 3.0*rho[3]*dx11*vthsq+0.8660254037844387*rhou1[5]*uvar1[7]*dx11+0.8660254037844387*uvar1[5]*rhou1[7]*dx11+0.7745966692414834*rhou1[3]*uvar1[6]*dx11+0.7745966692414834*uvar1[3]*rhou1[6]*dx11+0.7745966692414833*rhou1[1]*uvar1[4]*dx11+0.7745966692414833*uvar1[1]*rhou1[4]*dx11+0.8660254037844386*rhou1[2]*uvar1[3]*dx11+0.8660254037844386*uvar1[2]*rhou1[3]*dx11+0.8660254037844386*rhou1[0]*uvar1[1]*dx11+0.8660254037844386*uvar1[0]*rhou1[1]*dx11+0.7745966692414834*rhou1[3]*uvar0[7]*dx10+0.7745966692414834*uvar0[3]*rhou1[7]*dx10+0.8660254037844387*rhou1[4]*uvar0[6]*dx10+0.8660254037844387*uvar0[4]*rhou1[6]*dx10+0.7745966692414833*rhou1[2]*uvar0[5]*dx10+0.7745966692414833*uvar0[2]*rhou1[5]*dx10+0.8660254037844386*rhou1[1]*uvar0[3]*dx10+0.8660254037844386*uvar0[1]*rhou1[3]*dx10+0.8660254037844386*rhou1[0]*uvar0[2]*dx10+0.8660254037844386*uvar0[0]*rhou1[2]*dx10; 
  outrhou1[4] += 1.936491673103709*rhou1[5]*uvar0[7]*dx10+1.936491673103709*uvar0[5]*rhou1[7]*dx10+1.732050807568877*rhou1[3]*uvar0[6]*dx10+1.732050807568877*uvar0[3]*rhou1[6]*dx10+1.732050807568877*rhou1[1]*uvar0[4]*dx10+1.732050807568877*uvar0[1]*rhou1[4]*dx10+1.936491673103709*rhou1[2]*uvar0[3]*dx10+1.936491673103709*uvar0[2]*rhou1[3]*dx10+1.936491673103709*rhou1[0]*uvar0[1]*dx10+1.936491673103709*uvar0[0]*rhou1[1]*dx10; 
  outrhou1[5] += 15.0*rho[5]*dx11*vthsq+1.732050807568877*rhou1[3]*uvar1[7]*dx11+1.732050807568877*uvar1[3]*rhou1[7]*dx11+1.936491673103709*rhou1[4]*uvar1[6]*dx11+1.936491673103709*uvar1[4]*rhou1[6]*dx11+1.732050807568877*rhou1[2]*uvar1[5]*dx11+1.732050807568877*uvar1[2]*rhou1[5]*dx11+1.936491673103709*rhou1[1]*uvar1[3]*dx11+1.936491673103709*uvar1[1]*rhou1[3]*dx11+1.936491673103709*rhou1[0]*uvar1[2]*dx11+1.936491673103709*uvar1[0]*rhou1[2]*dx11; 
  outrhou1[6] += 3.0*rho[6]*dx11*vthsq+0.7745966692414834*rhou1[7]*uvar1[7]*dx11+0.5532833351724881*rhou1[6]*uvar1[6]*dx11+0.8660254037844386*rhou1[2]*uvar1[6]*dx11+0.8660254037844386*uvar1[2]*rhou1[6]*dx11+0.5532833351724881*rhou1[4]*uvar1[4]*dx11+0.8660254037844387*rhou1[0]*uvar1[4]*dx11+0.8660254037844387*uvar1[0]*rhou1[4]*dx11+0.7745966692414834*rhou1[3]*uvar1[3]*dx11+0.7745966692414834*rhou1[1]*uvar1[1]*dx11+1.549193338482967*rhou1[6]*uvar0[7]*dx10+1.732050807568877*rhou1[2]*uvar0[7]*dx10+1.549193338482967*uvar0[6]*rhou1[7]*dx10+1.732050807568877*uvar0[2]*rhou1[7]*dx10+1.732050807568877*rhou1[1]*uvar0[6]*dx10+1.732050807568877*uvar0[1]*rhou1[6]*dx10+1.732050807568877*rhou1[3]*uvar0[5]*dx10+1.732050807568877*uvar0[3]*rhou1[5]*dx10+1.732050807568877*rhou1[3]*uvar0[4]*dx10+1.732050807568877*uvar0[3]*rhou1[4]*dx10+1.936491673103709*rhou1[0]*uvar0[3]*dx10+1.936491673103709*uvar0[0]*rhou1[3]*dx10+1.936491673103709*rhou1[1]*uvar0[2]*dx10+1.936491673103709*uvar0[1]*rhou1[2]*dx10; 
  outrhou1[7] += 15.0*rho[7]*dx11*vthsq+1.549193338482967*rhou1[6]*uvar1[7]*dx11+1.732050807568877*rhou1[2]*uvar1[7]*dx11+1.549193338482967*uvar1[6]*rhou1[7]*dx11+1.732050807568877*uvar1[2]*rhou1[7]*dx11+1.732050807568877*rhou1[1]*uvar1[6]*dx11+1.732050807568877*uvar1[1]*rhou1[6]*dx11+1.732050807568877*rhou1[3]*uvar1[5]*dx11+1.732050807568877*uvar1[3]*rhou1[5]*dx11+1.732050807568877*rhou1[3]*uvar1[4]*dx11+1.732050807568877*uvar1[3]*rhou1[4]*dx11+1.936491673103709*rhou1[0]*uvar1[3]*dx11+1.936491673103709*uvar1[0]*rhou1[3]*dx11+1.936491673103709*rhou1[1]*uvar1[2]*dx11+1.936491673103709*uvar1[1]*rhou1[2]*dx11+0.5532833351724881*rhou1[7]*uvar0[7]*dx10+0.8660254037844386*rhou1[1]*uvar0[7]*dx10+0.8660254037844386*uvar0[1]*rhou1[7]*dx10+0.7745966692414834*rhou1[6]*uvar0[6]*dx10+0.5532833351724881*rhou1[5]*uvar0[5]*dx10+0.8660254037844387*rhou1[0]*uvar0[5]*dx10+0.8660254037844387*uvar0[0]*rhou1[5]*dx10+0.7745966692414834*rhou1[3]*uvar0[3]*dx10+0.7745966692414834*rhou1[2]*uvar0[2]*dx10; 

  return alpha_mid; 
} 
