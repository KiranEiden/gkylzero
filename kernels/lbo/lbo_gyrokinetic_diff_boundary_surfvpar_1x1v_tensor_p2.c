#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH void lbo_gyrokinetic_diff_boundary_surfvpar_1x1v_tensor_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  // w[2]: Cell-center coordinates. 
  // dxv[2]: Cell spacing. 
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[6]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fskin/edge: Distribution function in cells 
  // out: Incremented distribution function in cell 

  const double *nuVtSqSum = &nuPrimMomsSum[3];

  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 

  double bFacFskin[9] = {0.0}; 
  bFacFskin[5] = 6.708203932499369*fskin[0]; 
  bFacFskin[7] = 6.708203932499369*fskin[1]; 
  bFacFskin[8] = 6.708203932499369*fskin[4]; 

  double vol_incr[9] = {0.0};
  vol_incr[5] = 0.105409255338946*(6.708203932499369*nuVtSqSum[2]*bFacFskin[8]+6.708203932499369*nuVtSqSum[1]*bFacFskin[7]+6.708203932499369*nuVtSqSum[0]*bFacFskin[5]); 
  vol_incr[7] = 0.02581988897471611*(24.49489742783179*nuVtSqSum[1]*bFacFskin[8]+(24.49489742783179*nuVtSqSum[2]+27.38612787525831*nuVtSqSum[0])*bFacFskin[7]+27.38612787525831*nuVtSqSum[1]*bFacFskin[5]); 
  vol_incr[8] = 0.004761904761904762*((94.86832980505142*nuVtSqSum[2]+148.492424049175*nuVtSqSum[0])*bFacFskin[8]+132.815661727072*nuVtSqSum[1]*bFacFskin[7]+148.492424049175*nuVtSqSum[2]*bFacFskin[5]); 

  double edgeSurf_incr[9] = {0.0}; 
  double boundSurf_incr[9] = {0.0}; 

  if (edge == -1) { 

  double edgeSurf[9] = {0.0}; 
  edgeSurf[0] = (-0.6708203932499369*fskin[5])+0.6708203932499369*fedge[5]-1.190784930203603*fskin[2]-1.190784930203603*fedge[2]-0.9375*fskin[0]+0.9375*fedge[0]; 
  edgeSurf[1] = (-0.6708203932499369*fskin[7])+0.6708203932499369*fedge[7]-1.190784930203603*fskin[3]-1.190784930203603*fedge[3]-0.9375*fskin[1]+0.9375*fedge[1]; 
  edgeSurf[2] = (-1.585502557353661*fskin[5])+0.7382874503707888*fedge[5]-2.671875*fskin[2]-1.453125*fedge[2]-2.056810333988042*fskin[0]+1.190784930203603*fedge[0]; 
  edgeSurf[3] = (-1.585502557353661*fskin[7])+0.7382874503707888*fedge[7]-2.671875*fskin[3]-1.453125*fedge[3]-2.056810333988042*fskin[1]+1.190784930203603*fedge[1]; 
  edgeSurf[4] = (-0.6708203932499369*fskin[8])+0.6708203932499369*fedge[8]-1.190784930203603*fskin[6]-1.190784930203603*fedge[6]-0.9375*fskin[4]+0.9375*fedge[4]; 
  edgeSurf[5] = (-3.140625*fskin[5])-0.140625*fedge[5]-5.022775277112744*fskin[2]-0.3025768239224545*fedge[2]-3.773364712030896*fskin[0]+0.4192627457812106*fedge[0]; 
  edgeSurf[6] = (-1.585502557353661*fskin[8])+0.7382874503707888*fedge[8]-2.671875*fskin[6]-1.453125*fedge[6]-2.056810333988042*fskin[4]+1.190784930203603*fedge[4]; 
  edgeSurf[7] = (-3.140625*fskin[7])-0.140625*fedge[7]-5.022775277112744*fskin[3]-0.3025768239224544*fedge[3]-3.773364712030894*fskin[1]+0.4192627457812105*fedge[1]; 
  edgeSurf[8] = (-3.140625*fskin[8])-0.140625*fedge[8]-5.022775277112744*fskin[6]-0.3025768239224544*fedge[6]-3.773364712030896*fskin[4]+0.4192627457812106*fedge[4]; 

  double boundSurf[9] = {0.0}; 
  boundSurf[2] = 1.936491673103709*fskin[5]-1.5*fskin[2]+0.8660254037844386*fskin[0]; 
  boundSurf[3] = 1.936491673103709*fskin[7]-1.5*fskin[3]+0.8660254037844386*fskin[1]; 
  boundSurf[5] = (-7.5*fskin[5])+5.809475019311125*fskin[2]-3.354101966249685*fskin[0]; 
  boundSurf[6] = 1.936491673103709*fskin[8]-1.5*fskin[6]+0.8660254037844387*fskin[4]; 
  boundSurf[7] = (-7.5*fskin[7])+5.809475019311126*fskin[3]-3.354101966249684*fskin[1]; 
  boundSurf[8] = (-7.5*fskin[8])+5.809475019311126*fskin[6]-3.354101966249685*fskin[4]; 

  edgeSurf_incr[0] = 0.7071067811865476*(nuVtSqSum[2]*edgeSurf[4]+edgeSurf[1]*nuVtSqSum[1]+edgeSurf[0]*nuVtSqSum[0]); 
  edgeSurf_incr[1] = 0.1*(6.324555320336761*(nuVtSqSum[1]*edgeSurf[4]+edgeSurf[1]*nuVtSqSum[2])+7.071067811865476*(edgeSurf[0]*nuVtSqSum[1]+nuVtSqSum[0]*edgeSurf[1])); 
  edgeSurf_incr[2] = 0.03333333333333333*(21.21320343559643*nuVtSqSum[2]*edgeSurf[6]+21.21320343559643*(nuVtSqSum[1]*edgeSurf[3]+nuVtSqSum[0]*edgeSurf[2])); 
  edgeSurf_incr[3] = 0.03333333333333333*(18.97366596101028*nuVtSqSum[1]*edgeSurf[6]+(18.97366596101028*nuVtSqSum[2]+21.21320343559643*nuVtSqSum[0])*edgeSurf[3]+21.21320343559643*nuVtSqSum[1]*edgeSurf[2]); 
  edgeSurf_incr[4] = 0.01428571428571429*((31.62277660168381*nuVtSqSum[2]+49.49747468305833*nuVtSqSum[0])*edgeSurf[4]+49.49747468305833*edgeSurf[0]*nuVtSqSum[2]+44.27188724235732*edgeSurf[1]*nuVtSqSum[1]); 
  edgeSurf_incr[5] = 0.03333333333333333*(21.21320343559643*nuVtSqSum[2]*edgeSurf[8]+21.21320343559643*nuVtSqSum[1]*edgeSurf[7]+21.21320343559643*nuVtSqSum[0]*edgeSurf[5]); 
  edgeSurf_incr[6] = 0.004761904761904762*((94.86832980505142*nuVtSqSum[2]+148.492424049175*nuVtSqSum[0])*edgeSurf[6]+132.815661727072*nuVtSqSum[1]*edgeSurf[3]+148.492424049175*edgeSurf[2]*nuVtSqSum[2]); 
  edgeSurf_incr[7] = 0.03333333333333333*(18.97366596101028*nuVtSqSum[1]*edgeSurf[8]+(18.97366596101028*nuVtSqSum[2]+21.21320343559643*nuVtSqSum[0])*edgeSurf[7]+21.21320343559643*nuVtSqSum[1]*edgeSurf[5]); 
  edgeSurf_incr[8] = 0.004761904761904762*((94.86832980505142*nuVtSqSum[2]+148.492424049175*nuVtSqSum[0])*edgeSurf[8]+132.815661727072*nuVtSqSum[1]*edgeSurf[7]+148.492424049175*nuVtSqSum[2]*edgeSurf[5]); 

  boundSurf_incr[0] = 0.7071067811865476*(nuVtSqSum[2]*boundSurf[4]+boundSurf[1]*nuVtSqSum[1]+boundSurf[0]*nuVtSqSum[0]); 
  boundSurf_incr[1] = 0.1*(6.324555320336761*(nuVtSqSum[1]*boundSurf[4]+boundSurf[1]*nuVtSqSum[2])+7.071067811865476*(boundSurf[0]*nuVtSqSum[1]+nuVtSqSum[0]*boundSurf[1])); 
  boundSurf_incr[2] = 0.03333333333333333*(21.21320343559643*nuVtSqSum[2]*boundSurf[6]+21.21320343559643*(nuVtSqSum[1]*boundSurf[3]+nuVtSqSum[0]*boundSurf[2])); 
  boundSurf_incr[3] = 0.03333333333333333*(18.97366596101028*nuVtSqSum[1]*boundSurf[6]+(18.97366596101028*nuVtSqSum[2]+21.21320343559643*nuVtSqSum[0])*boundSurf[3]+21.21320343559643*nuVtSqSum[1]*boundSurf[2]); 
  boundSurf_incr[4] = 0.01428571428571429*((31.62277660168381*nuVtSqSum[2]+49.49747468305833*nuVtSqSum[0])*boundSurf[4]+49.49747468305833*boundSurf[0]*nuVtSqSum[2]+44.27188724235732*boundSurf[1]*nuVtSqSum[1]); 
  boundSurf_incr[5] = 0.03333333333333333*(21.21320343559643*nuVtSqSum[2]*boundSurf[8]+21.21320343559643*nuVtSqSum[1]*boundSurf[7]+21.21320343559643*nuVtSqSum[0]*boundSurf[5]); 
  boundSurf_incr[6] = 0.004761904761904762*((94.86832980505142*nuVtSqSum[2]+148.492424049175*nuVtSqSum[0])*boundSurf[6]+132.815661727072*nuVtSqSum[1]*boundSurf[3]+148.492424049175*boundSurf[2]*nuVtSqSum[2]); 
  boundSurf_incr[7] = 0.03333333333333333*(18.97366596101028*nuVtSqSum[1]*boundSurf[8]+(18.97366596101028*nuVtSqSum[2]+21.21320343559643*nuVtSqSum[0])*boundSurf[7]+21.21320343559643*nuVtSqSum[1]*boundSurf[5]); 
  boundSurf_incr[8] = 0.004761904761904762*((94.86832980505142*nuVtSqSum[2]+148.492424049175*nuVtSqSum[0])*boundSurf[8]+132.815661727072*nuVtSqSum[1]*boundSurf[7]+148.492424049175*nuVtSqSum[2]*boundSurf[5]); 


  } else { 

  double edgeSurf[9] = {0.0}; 
  edgeSurf[0] = -0.0125*(53.66563145999496*fskin[5]-53.66563145999496*fedge[5]-95.26279441628824*(fskin[2]+fedge[2])+75.0*fskin[0]-75.0*fedge[0]); 
  edgeSurf[1] = -0.0125*(53.66563145999495*fskin[7]-53.66563145999495*fedge[7]-95.26279441628824*(fskin[3]+fedge[3])+75.0*fskin[1]-75.0*fedge[1]); 
  edgeSurf[2] = 0.003125*(507.3608183531716*fskin[5]-236.2519841186524*fedge[5]-855.0*fskin[2]-465.0*fedge[2]+658.1793068761733*fskin[0]-381.051177665153*fedge[0]); 
  edgeSurf[3] = 0.003125*(507.3608183531716*fskin[7]-236.2519841186524*fedge[7]-855.0*fskin[3]-465.0*fedge[3]+658.1793068761733*fskin[1]-381.051177665153*fedge[1]); 
  edgeSurf[4] = -0.0125*(53.66563145999496*fskin[8]-53.66563145999496*fedge[8]-95.26279441628826*(fskin[6]+fedge[6])+75.0*fskin[4]-75.0*fedge[4]); 
  edgeSurf[5] = -0.015625*(201.0*fskin[5]+9.0*fedge[5]-321.4576177352156*fskin[2]-19.36491673103709*fedge[2]+241.4953415699773*fskin[0]-26.83281572999748*fedge[0]); 
  edgeSurf[6] = 0.003125*(507.3608183531716*fskin[8]-236.2519841186524*fedge[8]-855.0*fskin[6]-465.0*fedge[6]+658.1793068761734*fskin[4]-381.051177665153*fedge[4]); 
  edgeSurf[7] = -0.015625*(201.0*fskin[7]+9.0*fedge[7]-321.4576177352156*fskin[3]-19.36491673103708*fedge[3]+241.4953415699772*fskin[1]-26.83281572999747*fedge[1]); 
  edgeSurf[8] = -0.015625*(201.0*fskin[8]+9.0*fedge[8]-321.4576177352156*fskin[6]-19.36491673103708*fedge[6]+241.4953415699773*fskin[4]-26.83281572999748*fedge[4]); 

  double boundSurf[9] = {0.0}; 
  boundSurf[2] = -0.5*(3.872983346207417*fskin[5]+3.0*fskin[2]+1.732050807568877*fskin[0]); 
  boundSurf[3] = -0.5*(3.872983346207417*fskin[7]+3.0*fskin[3]+1.732050807568877*fskin[1]); 
  boundSurf[5] = -0.5*(15.0*fskin[5]+11.61895003862225*fskin[2]+6.708203932499369*fskin[0]); 
  boundSurf[6] = -0.1*(19.36491673103708*fskin[8]+15.0*fskin[6]+8.660254037844387*fskin[4]); 
  boundSurf[7] = -0.5*(15.0*fskin[7]+11.61895003862225*fskin[3]+6.708203932499369*fskin[1]); 
  boundSurf[8] = -0.5*(15.0*fskin[8]+11.61895003862225*fskin[6]+6.708203932499369*fskin[4]); 

  edgeSurf_incr[0] = 0.7071067811865476*(nuVtSqSum[2]*edgeSurf[4]+edgeSurf[1]*nuVtSqSum[1]+edgeSurf[0]*nuVtSqSum[0]); 
  edgeSurf_incr[1] = 0.1*(6.324555320336761*(nuVtSqSum[1]*edgeSurf[4]+edgeSurf[1]*nuVtSqSum[2])+7.071067811865476*(edgeSurf[0]*nuVtSqSum[1]+nuVtSqSum[0]*edgeSurf[1])); 
  edgeSurf_incr[2] = 0.03333333333333333*(21.21320343559643*nuVtSqSum[2]*edgeSurf[6]+21.21320343559643*(nuVtSqSum[1]*edgeSurf[3]+nuVtSqSum[0]*edgeSurf[2])); 
  edgeSurf_incr[3] = 0.03333333333333333*(18.97366596101028*nuVtSqSum[1]*edgeSurf[6]+(18.97366596101028*nuVtSqSum[2]+21.21320343559643*nuVtSqSum[0])*edgeSurf[3]+21.21320343559643*nuVtSqSum[1]*edgeSurf[2]); 
  edgeSurf_incr[4] = 0.01428571428571429*((31.62277660168381*nuVtSqSum[2]+49.49747468305833*nuVtSqSum[0])*edgeSurf[4]+49.49747468305833*edgeSurf[0]*nuVtSqSum[2]+44.27188724235732*edgeSurf[1]*nuVtSqSum[1]); 
  edgeSurf_incr[5] = 0.03333333333333333*(21.21320343559643*nuVtSqSum[2]*edgeSurf[8]+21.21320343559643*nuVtSqSum[1]*edgeSurf[7]+21.21320343559643*nuVtSqSum[0]*edgeSurf[5]); 
  edgeSurf_incr[6] = 0.004761904761904762*((94.86832980505142*nuVtSqSum[2]+148.492424049175*nuVtSqSum[0])*edgeSurf[6]+132.815661727072*nuVtSqSum[1]*edgeSurf[3]+148.492424049175*edgeSurf[2]*nuVtSqSum[2]); 
  edgeSurf_incr[7] = 0.03333333333333333*(18.97366596101028*nuVtSqSum[1]*edgeSurf[8]+(18.97366596101028*nuVtSqSum[2]+21.21320343559643*nuVtSqSum[0])*edgeSurf[7]+21.21320343559643*nuVtSqSum[1]*edgeSurf[5]); 
  edgeSurf_incr[8] = 0.004761904761904762*((94.86832980505142*nuVtSqSum[2]+148.492424049175*nuVtSqSum[0])*edgeSurf[8]+132.815661727072*nuVtSqSum[1]*edgeSurf[7]+148.492424049175*nuVtSqSum[2]*edgeSurf[5]); 

  boundSurf_incr[0] = 0.7071067811865476*(nuVtSqSum[2]*boundSurf[4]+boundSurf[1]*nuVtSqSum[1]+boundSurf[0]*nuVtSqSum[0]); 
  boundSurf_incr[1] = 0.1*(6.324555320336761*(nuVtSqSum[1]*boundSurf[4]+boundSurf[1]*nuVtSqSum[2])+7.071067811865476*(boundSurf[0]*nuVtSqSum[1]+nuVtSqSum[0]*boundSurf[1])); 
  boundSurf_incr[2] = 0.03333333333333333*(21.21320343559643*nuVtSqSum[2]*boundSurf[6]+21.21320343559643*(nuVtSqSum[1]*boundSurf[3]+nuVtSqSum[0]*boundSurf[2])); 
  boundSurf_incr[3] = 0.03333333333333333*(18.97366596101028*nuVtSqSum[1]*boundSurf[6]+(18.97366596101028*nuVtSqSum[2]+21.21320343559643*nuVtSqSum[0])*boundSurf[3]+21.21320343559643*nuVtSqSum[1]*boundSurf[2]); 
  boundSurf_incr[4] = 0.01428571428571429*((31.62277660168381*nuVtSqSum[2]+49.49747468305833*nuVtSqSum[0])*boundSurf[4]+49.49747468305833*boundSurf[0]*nuVtSqSum[2]+44.27188724235732*boundSurf[1]*nuVtSqSum[1]); 
  boundSurf_incr[5] = 0.03333333333333333*(21.21320343559643*nuVtSqSum[2]*boundSurf[8]+21.21320343559643*nuVtSqSum[1]*boundSurf[7]+21.21320343559643*nuVtSqSum[0]*boundSurf[5]); 
  boundSurf_incr[6] = 0.004761904761904762*((94.86832980505142*nuVtSqSum[2]+148.492424049175*nuVtSqSum[0])*boundSurf[6]+132.815661727072*nuVtSqSum[1]*boundSurf[3]+148.492424049175*boundSurf[2]*nuVtSqSum[2]); 
  boundSurf_incr[7] = 0.03333333333333333*(18.97366596101028*nuVtSqSum[1]*boundSurf[8]+(18.97366596101028*nuVtSqSum[2]+21.21320343559643*nuVtSqSum[0])*boundSurf[7]+21.21320343559643*nuVtSqSum[1]*boundSurf[5]); 
  boundSurf_incr[8] = 0.004761904761904762*((94.86832980505142*nuVtSqSum[2]+148.492424049175*nuVtSqSum[0])*boundSurf[8]+132.815661727072*nuVtSqSum[1]*boundSurf[7]+148.492424049175*nuVtSqSum[2]*boundSurf[5]); 

  } 

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*rdvSq4; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*rdvSq4; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*rdvSq4; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*rdvSq4; 
  out[4] += (vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*rdvSq4; 
  out[5] += (vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*rdvSq4; 
  out[6] += (vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*rdvSq4; 
  out[7] += (vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*rdvSq4; 
  out[8] += (vol_incr[8]+edgeSurf_incr[8]+boundSurf_incr[8])*rdvSq4; 
} 
