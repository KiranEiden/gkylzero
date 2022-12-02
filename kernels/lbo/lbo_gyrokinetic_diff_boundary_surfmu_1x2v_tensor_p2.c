#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH void lbo_gyrokinetic_diff_boundary_surfmu_1x2v_tensor_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  // w[3]: Cell-center coordinates. 
  // dxv[3]: Cell spacing. 
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[6]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fskin/edge: Distribution function in cells 
  // out: Incremented distribution function in cell 

  const double *nuVtSqSum = &nuPrimMomsSum[3];

  double rdvSq4 = 4.0/(dxv[2]*dxv[2]); 

  double surfVar_l = w[2]-0.5*dxv[2];
  double surfVar_r = w[2]+0.5*dxv[2];

  double facDiff[3]; 
  // Expand diffusion coefficient in conf basis.
  facDiff[0] = 1.414213562373095*(bmag_inv[2]*nuVtSqSum[2]+bmag_inv[1]*nuVtSqSum[1]+bmag_inv[0]*nuVtSqSum[0])*m_; 
  facDiff[1] = (1.264911064067352*(bmag_inv[1]*nuVtSqSum[2]+nuVtSqSum[1]*bmag_inv[2])+1.414213562373095*(bmag_inv[0]*nuVtSqSum[1]+nuVtSqSum[0]*bmag_inv[1]))*m_; 
  facDiff[2] = (0.9035079029052515*bmag_inv[2]*nuVtSqSum[2]+1.414213562373095*(bmag_inv[0]*nuVtSqSum[2]+nuVtSqSum[0]*bmag_inv[2])+1.264911064067352*bmag_inv[1]*nuVtSqSum[1])*m_; 

  double bFacFskin[27] = {0.0}; 
  bFacFskin[3] = 0.8660254037844386*fskin[0]*dxv[2]; 
  bFacFskin[5] = 0.8660254037844386*fskin[1]*dxv[2]; 
  bFacFskin[6] = 0.8660254037844386*dxv[2]*fskin[2]; 
  bFacFskin[9] = 3.872983346207417*dxv[2]*fskin[3]+6.708203932499369*fskin[0]*w[2]; 
  bFacFskin[10] = 0.8660254037844386*dxv[2]*fskin[4]; 
  bFacFskin[13] = 0.8660254037844387*dxv[2]*fskin[7]; 
  bFacFskin[14] = 0.8660254037844387*dxv[2]*fskin[8]; 
  bFacFskin[15] = 3.872983346207417*dxv[2]*fskin[5]+6.708203932499369*fskin[1]*w[2]; 
  bFacFskin[16] = 3.872983346207417*dxv[2]*fskin[6]+6.708203932499369*fskin[2]*w[2]; 
  bFacFskin[17] = 0.8660254037844387*dxv[2]*fskin[11]; 
  bFacFskin[18] = 0.8660254037844387*dxv[2]*fskin[12]; 
  bFacFskin[19] = 3.872983346207417*dxv[2]*fskin[10]+6.708203932499369*w[2]*fskin[4]; 
  bFacFskin[21] = 3.872983346207417*dxv[2]*fskin[13]+6.708203932499369*w[2]*fskin[7]; 
  bFacFskin[22] = 3.872983346207417*dxv[2]*fskin[14]+6.708203932499369*w[2]*fskin[8]; 
  bFacFskin[23] = 0.8660254037844386*dxv[2]*fskin[20]; 
  bFacFskin[24] = 3.872983346207417*dxv[2]*fskin[17]+6.708203932499369*w[2]*fskin[11]; 
  bFacFskin[25] = 3.872983346207417*dxv[2]*fskin[18]+6.708203932499369*w[2]*fskin[12]; 
  bFacFskin[26] = 3.872983346207417*dxv[2]*fskin[23]+6.708203932499369*w[2]*fskin[20]; 

  double vol_incr[27] = {0.0};
  vol_incr[3] = 0.04714045207910316*(15.0*facDiff[2]*bFacFskin[13]+15.0*(facDiff[1]*bFacFskin[5]+facDiff[0]*bFacFskin[3])); 
  vol_incr[5] = 0.04714045207910316*(13.41640786499874*facDiff[1]*bFacFskin[13]+(13.41640786499874*facDiff[2]+15.0*facDiff[0])*bFacFskin[5]+15.0*facDiff[1]*bFacFskin[3]); 
  vol_incr[6] = 0.7071067811865475*(facDiff[2]*bFacFskin[17]+facDiff[1]*bFacFskin[10]+facDiff[0]*bFacFskin[6]); 
  vol_incr[9] = 0.04714045207910316*(15.0*facDiff[2]*bFacFskin[21]+15.0*facDiff[1]*bFacFskin[15]+15.0*facDiff[0]*bFacFskin[9]); 
  vol_incr[10] = 0.1414213562373095*(4.47213595499958*facDiff[1]*bFacFskin[17]+(4.47213595499958*facDiff[2]+5.0*facDiff[0])*bFacFskin[10]+5.0*facDiff[1]*bFacFskin[6]); 
  vol_incr[13] = 0.006734350297014738*((67.0820393249937*facDiff[2]+105.0*facDiff[0])*bFacFskin[13]+93.91485505499116*facDiff[1]*bFacFskin[5]+105.0*facDiff[2]*bFacFskin[3]); 
  vol_incr[14] = 0.04714045207910316*(15.0*(facDiff[2]*bFacFskin[23]+facDiff[1]*bFacFskin[18])+15.0*facDiff[0]*bFacFskin[14]); 
  vol_incr[15] = 0.04714045207910316*(13.41640786499874*facDiff[1]*bFacFskin[21]+(13.41640786499874*facDiff[2]+15.0*facDiff[0])*bFacFskin[15]+15.0*facDiff[1]*bFacFskin[9]); 
  vol_incr[16] = 0.04714045207910316*(15.0*(facDiff[2]*bFacFskin[24]+facDiff[1]*bFacFskin[19])+15.0*facDiff[0]*bFacFskin[16]); 
  vol_incr[17] = 0.02020305089104421*((22.3606797749979*facDiff[2]+35.0*facDiff[0])*bFacFskin[17]+31.30495168499706*facDiff[1]*bFacFskin[10]+35.0*facDiff[2]*bFacFskin[6]); 
  vol_incr[18] = 0.04714045207910316*(13.41640786499874*facDiff[1]*bFacFskin[23]+(13.41640786499874*facDiff[2]+15.0*facDiff[0])*bFacFskin[18]+15.0*facDiff[1]*bFacFskin[14]); 
  vol_incr[19] = 0.04714045207910316*(13.41640786499874*facDiff[1]*bFacFskin[24]+(13.41640786499874*facDiff[2]+15.0*facDiff[0])*bFacFskin[19]+15.0*facDiff[1]*bFacFskin[16]); 
  vol_incr[21] = 0.006734350297014738*((67.0820393249937*facDiff[2]+105.0*facDiff[0])*bFacFskin[21]+93.91485505499116*facDiff[1]*bFacFskin[15]+105.0*facDiff[2]*bFacFskin[9]); 
  vol_incr[22] = 0.7071067811865475*(facDiff[2]*bFacFskin[26]+facDiff[1]*bFacFskin[25]+facDiff[0]*bFacFskin[22]); 
  vol_incr[23] = 0.006734350297014738*((67.0820393249937*facDiff[2]+105.0*facDiff[0])*bFacFskin[23]+93.91485505499116*facDiff[1]*bFacFskin[18]+105.0*facDiff[2]*bFacFskin[14]); 
  vol_incr[24] = 0.006734350297014738*((67.0820393249937*facDiff[2]+105.0*facDiff[0])*bFacFskin[24]+93.91485505499116*facDiff[1]*bFacFskin[19]+105.0*facDiff[2]*bFacFskin[16]); 
  vol_incr[25] = 0.1414213562373095*(4.47213595499958*facDiff[1]*bFacFskin[26]+(4.47213595499958*facDiff[2]+5.0*facDiff[0])*bFacFskin[25]+5.0*facDiff[1]*bFacFskin[22]); 
  vol_incr[26] = 0.02020305089104421*((22.3606797749979*facDiff[2]+35.0*facDiff[0])*bFacFskin[26]+31.30495168499706*facDiff[1]*bFacFskin[25]+35.0*facDiff[2]*bFacFskin[22]); 

  double edgeSurf_incr[27] = {0.0}; 
  double boundSurf_incr[27] = {0.0}; 

  if (edge == -1) { 

  double edgeSurf[27] = {0.0}; 
  edgeSurf[0] = (-0.6708203932499369*fskin[9]*surfVar_r)+0.6708203932499369*fedge[9]*surfVar_r-1.190784930203603*fskin[3]*surfVar_r-1.190784930203603*fedge[3]*surfVar_r-0.9375*fskin[0]*surfVar_r+0.9375*fedge[0]*surfVar_r; 
  edgeSurf[1] = (-0.6708203932499369*fskin[15]*surfVar_r)+0.6708203932499369*fedge[15]*surfVar_r-1.190784930203603*fskin[5]*surfVar_r-1.190784930203603*fedge[5]*surfVar_r-0.9375*fskin[1]*surfVar_r+0.9375*fedge[1]*surfVar_r; 
  edgeSurf[2] = (-0.6708203932499369*fskin[16]*surfVar_r)+0.6708203932499369*fedge[16]*surfVar_r-1.190784930203603*fskin[6]*surfVar_r-1.190784930203603*fedge[6]*surfVar_r-0.9375*fskin[2]*surfVar_r+0.9375*fedge[2]*surfVar_r; 
  edgeSurf[3] = (-1.585502557353661*fskin[9]*surfVar_r)+0.7382874503707888*fedge[9]*surfVar_r-2.671875*fskin[3]*surfVar_r-1.453125*fedge[3]*surfVar_r-2.056810333988042*fskin[0]*surfVar_r+1.190784930203603*fedge[0]*surfVar_r; 
  edgeSurf[4] = (-0.6708203932499369*fskin[19]*surfVar_r)+0.6708203932499369*fedge[19]*surfVar_r-1.190784930203603*fskin[10]*surfVar_r-1.190784930203603*fedge[10]*surfVar_r-0.9375*fskin[4]*surfVar_r+0.9375*fedge[4]*surfVar_r; 
  edgeSurf[5] = (-1.585502557353661*fskin[15]*surfVar_r)+0.7382874503707888*fedge[15]*surfVar_r-2.671875*fskin[5]*surfVar_r-1.453125*fedge[5]*surfVar_r-2.056810333988042*fskin[1]*surfVar_r+1.190784930203603*fedge[1]*surfVar_r; 
  edgeSurf[6] = (-1.585502557353661*fskin[16]*surfVar_r)+0.7382874503707888*fedge[16]*surfVar_r-2.671875*fskin[6]*surfVar_r-1.453125*fedge[6]*surfVar_r-2.056810333988042*fskin[2]*surfVar_r+1.190784930203603*fedge[2]*surfVar_r; 
  edgeSurf[7] = (-0.6708203932499369*fskin[21]*surfVar_r)+0.6708203932499369*fedge[21]*surfVar_r-1.190784930203603*fskin[13]*surfVar_r-1.190784930203603*fedge[13]*surfVar_r-0.9375*fskin[7]*surfVar_r+0.9375*fedge[7]*surfVar_r; 
  edgeSurf[8] = (-0.6708203932499369*fskin[22]*surfVar_r)+0.6708203932499369*fedge[22]*surfVar_r-1.190784930203603*fskin[14]*surfVar_r-1.190784930203603*fedge[14]*surfVar_r-0.9375*fskin[8]*surfVar_r+0.9375*fedge[8]*surfVar_r; 
  edgeSurf[9] = (-3.140625*fskin[9]*surfVar_r)-0.140625*fedge[9]*surfVar_r-5.022775277112744*fskin[3]*surfVar_r-0.3025768239224545*fedge[3]*surfVar_r-3.773364712030896*fskin[0]*surfVar_r+0.4192627457812106*fedge[0]*surfVar_r; 
  edgeSurf[10] = (-1.585502557353661*fskin[19]*surfVar_r)+0.7382874503707888*fedge[19]*surfVar_r-2.671875*fskin[10]*surfVar_r-1.453125*fedge[10]*surfVar_r-2.056810333988042*fskin[4]*surfVar_r+1.190784930203603*fedge[4]*surfVar_r; 
  edgeSurf[11] = (-0.6708203932499369*fskin[24]*surfVar_r)+0.6708203932499369*fedge[24]*surfVar_r-1.190784930203603*fskin[17]*surfVar_r-1.190784930203603*fedge[17]*surfVar_r-0.9375*fskin[11]*surfVar_r+0.9375*fedge[11]*surfVar_r; 
  edgeSurf[12] = (-0.6708203932499369*fskin[25]*surfVar_r)+0.6708203932499369*fedge[25]*surfVar_r-1.190784930203603*fskin[18]*surfVar_r-1.190784930203603*fedge[18]*surfVar_r-0.9375*fskin[12]*surfVar_r+0.9375*fedge[12]*surfVar_r; 
  edgeSurf[13] = (-1.585502557353661*fskin[21]*surfVar_r)+0.7382874503707888*fedge[21]*surfVar_r-2.671875*fskin[13]*surfVar_r-1.453125*fedge[13]*surfVar_r-2.056810333988042*fskin[7]*surfVar_r+1.190784930203603*fedge[7]*surfVar_r; 
  edgeSurf[14] = (-1.585502557353661*fskin[22]*surfVar_r)+0.7382874503707888*fedge[22]*surfVar_r-2.671875*fskin[14]*surfVar_r-1.453125*fedge[14]*surfVar_r-2.056810333988042*fskin[8]*surfVar_r+1.190784930203603*fedge[8]*surfVar_r; 
  edgeSurf[15] = (-3.140625*fskin[15]*surfVar_r)-0.140625*fedge[15]*surfVar_r-5.022775277112744*fskin[5]*surfVar_r-0.3025768239224544*fedge[5]*surfVar_r-3.773364712030894*fskin[1]*surfVar_r+0.4192627457812105*fedge[1]*surfVar_r; 
  edgeSurf[16] = (-3.140625*fskin[16]*surfVar_r)-0.140625*fedge[16]*surfVar_r-5.022775277112744*fskin[6]*surfVar_r-0.3025768239224544*fedge[6]*surfVar_r-3.773364712030894*fskin[2]*surfVar_r+0.4192627457812105*fedge[2]*surfVar_r; 
  edgeSurf[17] = (-1.585502557353661*fskin[24]*surfVar_r)+0.7382874503707888*fedge[24]*surfVar_r-2.671875*fskin[17]*surfVar_r-1.453125*fedge[17]*surfVar_r-2.056810333988042*fskin[11]*surfVar_r+1.190784930203603*fedge[11]*surfVar_r; 
  edgeSurf[18] = (-1.585502557353661*fskin[25]*surfVar_r)+0.7382874503707888*fedge[25]*surfVar_r-2.671875*fskin[18]*surfVar_r-1.453125*fedge[18]*surfVar_r-2.056810333988042*fskin[12]*surfVar_r+1.190784930203603*fedge[12]*surfVar_r; 
  edgeSurf[19] = (-3.140625*fskin[19]*surfVar_r)-0.140625*fedge[19]*surfVar_r-5.022775277112744*fskin[10]*surfVar_r-0.3025768239224545*fedge[10]*surfVar_r-3.773364712030896*fskin[4]*surfVar_r+0.4192627457812106*fedge[4]*surfVar_r; 
  edgeSurf[20] = (-0.6708203932499369*fskin[26]*surfVar_r)+0.6708203932499369*fedge[26]*surfVar_r-1.190784930203603*fskin[23]*surfVar_r-1.190784930203603*fedge[23]*surfVar_r-0.9375*fskin[20]*surfVar_r+0.9375*fedge[20]*surfVar_r; 
  edgeSurf[21] = (-3.140625*fskin[21]*surfVar_r)-0.140625*fedge[21]*surfVar_r-5.022775277112744*fskin[13]*surfVar_r-0.3025768239224544*fedge[13]*surfVar_r-3.773364712030896*fskin[7]*surfVar_r+0.4192627457812106*fedge[7]*surfVar_r; 
  edgeSurf[22] = (-3.140625*fskin[22]*surfVar_r)-0.140625*fedge[22]*surfVar_r-5.022775277112744*fskin[14]*surfVar_r-0.3025768239224544*fedge[14]*surfVar_r-3.773364712030896*fskin[8]*surfVar_r+0.4192627457812106*fedge[8]*surfVar_r; 
  edgeSurf[23] = (-1.585502557353661*fskin[26]*surfVar_r)+0.7382874503707888*fedge[26]*surfVar_r-2.671875*fskin[23]*surfVar_r-1.453125*fedge[23]*surfVar_r-2.056810333988042*fskin[20]*surfVar_r+1.190784930203603*fedge[20]*surfVar_r; 
  edgeSurf[24] = (-3.140625*fskin[24]*surfVar_r)-0.140625*fedge[24]*surfVar_r-5.022775277112744*fskin[17]*surfVar_r-0.3025768239224545*fedge[17]*surfVar_r-3.773364712030894*fskin[11]*surfVar_r+0.4192627457812105*fedge[11]*surfVar_r; 
  edgeSurf[25] = (-3.140625*fskin[25]*surfVar_r)-0.140625*fedge[25]*surfVar_r-5.022775277112744*fskin[18]*surfVar_r-0.3025768239224545*fedge[18]*surfVar_r-3.773364712030894*fskin[12]*surfVar_r+0.4192627457812105*fedge[12]*surfVar_r; 
  edgeSurf[26] = (-3.140625*fskin[26]*surfVar_r)-0.140625*fedge[26]*surfVar_r-5.022775277112744*fskin[23]*surfVar_r-0.3025768239224545*fedge[23]*surfVar_r-3.773364712030896*fskin[20]*surfVar_r+0.4192627457812106*fedge[20]*surfVar_r; 

  double boundSurf[27] = {0.0}; 
  boundSurf[3] = 1.936491673103709*fskin[9]*surfVar_l-1.5*fskin[3]*surfVar_l+0.8660254037844386*fskin[0]*surfVar_l; 
  boundSurf[5] = 1.936491673103709*fskin[15]*surfVar_l-1.5*fskin[5]*surfVar_l+0.8660254037844386*fskin[1]*surfVar_l; 
  boundSurf[6] = 1.936491673103709*fskin[16]*surfVar_l-1.5*fskin[6]*surfVar_l+0.8660254037844386*fskin[2]*surfVar_l; 
  boundSurf[9] = (-7.5*fskin[9]*surfVar_l)+5.809475019311125*fskin[3]*surfVar_l-3.354101966249685*fskin[0]*surfVar_l; 
  boundSurf[10] = 1.936491673103709*fskin[19]*surfVar_l-1.5*fskin[10]*surfVar_l+0.8660254037844386*fskin[4]*surfVar_l; 
  boundSurf[13] = 1.936491673103709*fskin[21]*surfVar_l-1.5*fskin[13]*surfVar_l+0.8660254037844387*fskin[7]*surfVar_l; 
  boundSurf[14] = 1.936491673103709*fskin[22]*surfVar_l-1.5*fskin[14]*surfVar_l+0.8660254037844387*fskin[8]*surfVar_l; 
  boundSurf[15] = (-7.5*fskin[15]*surfVar_l)+5.809475019311126*fskin[5]*surfVar_l-3.354101966249684*fskin[1]*surfVar_l; 
  boundSurf[16] = (-7.5*fskin[16]*surfVar_l)+5.809475019311126*fskin[6]*surfVar_l-3.354101966249684*fskin[2]*surfVar_l; 
  boundSurf[17] = 1.936491673103709*fskin[24]*surfVar_l-1.5*fskin[17]*surfVar_l+0.8660254037844387*fskin[11]*surfVar_l; 
  boundSurf[18] = 1.936491673103709*fskin[25]*surfVar_l-1.5*fskin[18]*surfVar_l+0.8660254037844387*fskin[12]*surfVar_l; 
  boundSurf[19] = (-7.5*fskin[19]*surfVar_l)+5.809475019311125*fskin[10]*surfVar_l-3.354101966249685*fskin[4]*surfVar_l; 
  boundSurf[21] = (-7.5*fskin[21]*surfVar_l)+5.809475019311126*fskin[13]*surfVar_l-3.354101966249685*fskin[7]*surfVar_l; 
  boundSurf[22] = (-7.5*fskin[22]*surfVar_l)+5.809475019311126*fskin[14]*surfVar_l-3.354101966249685*fskin[8]*surfVar_l; 
  boundSurf[23] = 1.936491673103709*fskin[26]*surfVar_l-1.5*fskin[23]*surfVar_l+0.8660254037844386*fskin[20]*surfVar_l; 
  boundSurf[24] = (-7.5*fskin[24]*surfVar_l)+5.809475019311125*fskin[17]*surfVar_l-3.354101966249684*fskin[11]*surfVar_l; 
  boundSurf[25] = (-7.5*fskin[25]*surfVar_l)+5.809475019311125*fskin[18]*surfVar_l-3.354101966249684*fskin[12]*surfVar_l; 
  boundSurf[26] = (-7.5*fskin[26]*surfVar_l)+5.809475019311125*fskin[23]*surfVar_l-3.354101966249685*fskin[20]*surfVar_l; 

  edgeSurf_incr[0] = 0.7071067811865475*(facDiff[2]*edgeSurf[7]+edgeSurf[1]*facDiff[1]+edgeSurf[0]*facDiff[0]); 
  edgeSurf_incr[1] = 0.1414213562373095*(4.47213595499958*(facDiff[1]*edgeSurf[7]+edgeSurf[1]*facDiff[2])+5.0*(edgeSurf[0]*facDiff[1]+facDiff[0]*edgeSurf[1])); 
  edgeSurf_incr[2] = 0.04714045207910316*(15.0*facDiff[2]*edgeSurf[11]+15.0*(facDiff[1]*edgeSurf[4]+facDiff[0]*edgeSurf[2])); 
  edgeSurf_incr[3] = 0.04714045207910316*(15.0*facDiff[2]*edgeSurf[13]+15.0*(facDiff[1]*edgeSurf[5]+facDiff[0]*edgeSurf[3])); 
  edgeSurf_incr[4] = 0.04714045207910316*(13.41640786499874*facDiff[1]*edgeSurf[11]+(13.41640786499874*facDiff[2]+15.0*facDiff[0])*edgeSurf[4]+15.0*facDiff[1]*edgeSurf[2]); 
  edgeSurf_incr[5] = 0.04714045207910316*(13.41640786499874*facDiff[1]*edgeSurf[13]+(13.41640786499874*facDiff[2]+15.0*facDiff[0])*edgeSurf[5]+15.0*facDiff[1]*edgeSurf[3]); 
  edgeSurf_incr[6] = 0.7071067811865475*(facDiff[2]*edgeSurf[17]+facDiff[1]*edgeSurf[10]+facDiff[0]*edgeSurf[6]); 
  edgeSurf_incr[7] = 0.02020305089104421*((22.3606797749979*facDiff[2]+35.0*facDiff[0])*edgeSurf[7]+35.0*edgeSurf[0]*facDiff[2]+31.30495168499706*edgeSurf[1]*facDiff[1]); 
  edgeSurf_incr[8] = 0.04714045207910316*(15.0*facDiff[2]*edgeSurf[20]+15.0*facDiff[1]*edgeSurf[12]+15.0*facDiff[0]*edgeSurf[8]); 
  edgeSurf_incr[9] = 0.04714045207910316*(15.0*facDiff[2]*edgeSurf[21]+15.0*facDiff[1]*edgeSurf[15]+15.0*facDiff[0]*edgeSurf[9]); 
  edgeSurf_incr[10] = 0.1414213562373095*(4.47213595499958*facDiff[1]*edgeSurf[17]+(4.47213595499958*facDiff[2]+5.0*facDiff[0])*edgeSurf[10]+5.0*facDiff[1]*edgeSurf[6]); 
  edgeSurf_incr[11] = 0.006734350297014738*((67.0820393249937*facDiff[2]+105.0*facDiff[0])*edgeSurf[11]+93.91485505499116*facDiff[1]*edgeSurf[4]+105.0*edgeSurf[2]*facDiff[2]); 
  edgeSurf_incr[12] = 0.04714045207910316*(13.41640786499874*facDiff[1]*edgeSurf[20]+(13.41640786499874*facDiff[2]+15.0*facDiff[0])*edgeSurf[12]+15.0*facDiff[1]*edgeSurf[8]); 
  edgeSurf_incr[13] = 0.006734350297014738*((67.0820393249937*facDiff[2]+105.0*facDiff[0])*edgeSurf[13]+93.91485505499116*facDiff[1]*edgeSurf[5]+105.0*facDiff[2]*edgeSurf[3]); 
  edgeSurf_incr[14] = 0.04714045207910316*(15.0*(facDiff[2]*edgeSurf[23]+facDiff[1]*edgeSurf[18])+15.0*facDiff[0]*edgeSurf[14]); 
  edgeSurf_incr[15] = 0.04714045207910316*(13.41640786499874*facDiff[1]*edgeSurf[21]+(13.41640786499874*facDiff[2]+15.0*facDiff[0])*edgeSurf[15]+15.0*facDiff[1]*edgeSurf[9]); 
  edgeSurf_incr[16] = 0.04714045207910316*(15.0*(facDiff[2]*edgeSurf[24]+facDiff[1]*edgeSurf[19])+15.0*facDiff[0]*edgeSurf[16]); 
  edgeSurf_incr[17] = 0.02020305089104421*((22.3606797749979*facDiff[2]+35.0*facDiff[0])*edgeSurf[17]+31.30495168499706*facDiff[1]*edgeSurf[10]+35.0*facDiff[2]*edgeSurf[6]); 
  edgeSurf_incr[18] = 0.04714045207910316*(13.41640786499874*facDiff[1]*edgeSurf[23]+(13.41640786499874*facDiff[2]+15.0*facDiff[0])*edgeSurf[18]+15.0*facDiff[1]*edgeSurf[14]); 
  edgeSurf_incr[19] = 0.04714045207910316*(13.41640786499874*facDiff[1]*edgeSurf[24]+(13.41640786499874*facDiff[2]+15.0*facDiff[0])*edgeSurf[19]+15.0*facDiff[1]*edgeSurf[16]); 
  edgeSurf_incr[20] = 0.006734350297014738*((67.0820393249937*facDiff[2]+105.0*facDiff[0])*edgeSurf[20]+93.91485505499116*facDiff[1]*edgeSurf[12]+105.0*facDiff[2]*edgeSurf[8]); 
  edgeSurf_incr[21] = 0.006734350297014738*((67.0820393249937*facDiff[2]+105.0*facDiff[0])*edgeSurf[21]+93.91485505499116*facDiff[1]*edgeSurf[15]+105.0*facDiff[2]*edgeSurf[9]); 
  edgeSurf_incr[22] = 0.7071067811865475*(facDiff[2]*edgeSurf[26]+facDiff[1]*edgeSurf[25]+facDiff[0]*edgeSurf[22]); 
  edgeSurf_incr[23] = 0.006734350297014738*((67.0820393249937*facDiff[2]+105.0*facDiff[0])*edgeSurf[23]+93.91485505499116*facDiff[1]*edgeSurf[18]+105.0*facDiff[2]*edgeSurf[14]); 
  edgeSurf_incr[24] = 0.006734350297014738*((67.0820393249937*facDiff[2]+105.0*facDiff[0])*edgeSurf[24]+93.91485505499116*facDiff[1]*edgeSurf[19]+105.0*facDiff[2]*edgeSurf[16]); 
  edgeSurf_incr[25] = 0.1414213562373095*(4.47213595499958*facDiff[1]*edgeSurf[26]+(4.47213595499958*facDiff[2]+5.0*facDiff[0])*edgeSurf[25]+5.0*facDiff[1]*edgeSurf[22]); 
  edgeSurf_incr[26] = 0.02020305089104421*((22.3606797749979*facDiff[2]+35.0*facDiff[0])*edgeSurf[26]+31.30495168499706*facDiff[1]*edgeSurf[25]+35.0*facDiff[2]*edgeSurf[22]); 

  boundSurf_incr[0] = 0.7071067811865475*(facDiff[2]*boundSurf[7]+boundSurf[1]*facDiff[1]+boundSurf[0]*facDiff[0]); 
  boundSurf_incr[1] = 0.1414213562373095*(4.47213595499958*(facDiff[1]*boundSurf[7]+boundSurf[1]*facDiff[2])+5.0*(boundSurf[0]*facDiff[1]+facDiff[0]*boundSurf[1])); 
  boundSurf_incr[2] = 0.04714045207910316*(15.0*facDiff[2]*boundSurf[11]+15.0*(facDiff[1]*boundSurf[4]+facDiff[0]*boundSurf[2])); 
  boundSurf_incr[3] = 0.04714045207910316*(15.0*facDiff[2]*boundSurf[13]+15.0*(facDiff[1]*boundSurf[5]+facDiff[0]*boundSurf[3])); 
  boundSurf_incr[4] = 0.04714045207910316*(13.41640786499874*facDiff[1]*boundSurf[11]+(13.41640786499874*facDiff[2]+15.0*facDiff[0])*boundSurf[4]+15.0*facDiff[1]*boundSurf[2]); 
  boundSurf_incr[5] = 0.04714045207910316*(13.41640786499874*facDiff[1]*boundSurf[13]+(13.41640786499874*facDiff[2]+15.0*facDiff[0])*boundSurf[5]+15.0*facDiff[1]*boundSurf[3]); 
  boundSurf_incr[6] = 0.7071067811865475*(facDiff[2]*boundSurf[17]+facDiff[1]*boundSurf[10]+facDiff[0]*boundSurf[6]); 
  boundSurf_incr[7] = 0.02020305089104421*((22.3606797749979*facDiff[2]+35.0*facDiff[0])*boundSurf[7]+35.0*boundSurf[0]*facDiff[2]+31.30495168499706*boundSurf[1]*facDiff[1]); 
  boundSurf_incr[8] = 0.04714045207910316*(15.0*facDiff[2]*boundSurf[20]+15.0*facDiff[1]*boundSurf[12]+15.0*facDiff[0]*boundSurf[8]); 
  boundSurf_incr[9] = 0.04714045207910316*(15.0*facDiff[2]*boundSurf[21]+15.0*facDiff[1]*boundSurf[15]+15.0*facDiff[0]*boundSurf[9]); 
  boundSurf_incr[10] = 0.1414213562373095*(4.47213595499958*facDiff[1]*boundSurf[17]+(4.47213595499958*facDiff[2]+5.0*facDiff[0])*boundSurf[10]+5.0*facDiff[1]*boundSurf[6]); 
  boundSurf_incr[11] = 0.006734350297014738*((67.0820393249937*facDiff[2]+105.0*facDiff[0])*boundSurf[11]+93.91485505499116*facDiff[1]*boundSurf[4]+105.0*boundSurf[2]*facDiff[2]); 
  boundSurf_incr[12] = 0.04714045207910316*(13.41640786499874*facDiff[1]*boundSurf[20]+(13.41640786499874*facDiff[2]+15.0*facDiff[0])*boundSurf[12]+15.0*facDiff[1]*boundSurf[8]); 
  boundSurf_incr[13] = 0.006734350297014738*((67.0820393249937*facDiff[2]+105.0*facDiff[0])*boundSurf[13]+93.91485505499116*facDiff[1]*boundSurf[5]+105.0*facDiff[2]*boundSurf[3]); 
  boundSurf_incr[14] = 0.04714045207910316*(15.0*(facDiff[2]*boundSurf[23]+facDiff[1]*boundSurf[18])+15.0*facDiff[0]*boundSurf[14]); 
  boundSurf_incr[15] = 0.04714045207910316*(13.41640786499874*facDiff[1]*boundSurf[21]+(13.41640786499874*facDiff[2]+15.0*facDiff[0])*boundSurf[15]+15.0*facDiff[1]*boundSurf[9]); 
  boundSurf_incr[16] = 0.04714045207910316*(15.0*(facDiff[2]*boundSurf[24]+facDiff[1]*boundSurf[19])+15.0*facDiff[0]*boundSurf[16]); 
  boundSurf_incr[17] = 0.02020305089104421*((22.3606797749979*facDiff[2]+35.0*facDiff[0])*boundSurf[17]+31.30495168499706*facDiff[1]*boundSurf[10]+35.0*facDiff[2]*boundSurf[6]); 
  boundSurf_incr[18] = 0.04714045207910316*(13.41640786499874*facDiff[1]*boundSurf[23]+(13.41640786499874*facDiff[2]+15.0*facDiff[0])*boundSurf[18]+15.0*facDiff[1]*boundSurf[14]); 
  boundSurf_incr[19] = 0.04714045207910316*(13.41640786499874*facDiff[1]*boundSurf[24]+(13.41640786499874*facDiff[2]+15.0*facDiff[0])*boundSurf[19]+15.0*facDiff[1]*boundSurf[16]); 
  boundSurf_incr[20] = 0.006734350297014738*((67.0820393249937*facDiff[2]+105.0*facDiff[0])*boundSurf[20]+93.91485505499116*facDiff[1]*boundSurf[12]+105.0*facDiff[2]*boundSurf[8]); 
  boundSurf_incr[21] = 0.006734350297014738*((67.0820393249937*facDiff[2]+105.0*facDiff[0])*boundSurf[21]+93.91485505499116*facDiff[1]*boundSurf[15]+105.0*facDiff[2]*boundSurf[9]); 
  boundSurf_incr[22] = 0.7071067811865475*(facDiff[2]*boundSurf[26]+facDiff[1]*boundSurf[25]+facDiff[0]*boundSurf[22]); 
  boundSurf_incr[23] = 0.006734350297014738*((67.0820393249937*facDiff[2]+105.0*facDiff[0])*boundSurf[23]+93.91485505499116*facDiff[1]*boundSurf[18]+105.0*facDiff[2]*boundSurf[14]); 
  boundSurf_incr[24] = 0.006734350297014738*((67.0820393249937*facDiff[2]+105.0*facDiff[0])*boundSurf[24]+93.91485505499116*facDiff[1]*boundSurf[19]+105.0*facDiff[2]*boundSurf[16]); 
  boundSurf_incr[25] = 0.1414213562373095*(4.47213595499958*facDiff[1]*boundSurf[26]+(4.47213595499958*facDiff[2]+5.0*facDiff[0])*boundSurf[25]+5.0*facDiff[1]*boundSurf[22]); 
  boundSurf_incr[26] = 0.02020305089104421*((22.3606797749979*facDiff[2]+35.0*facDiff[0])*boundSurf[26]+31.30495168499706*facDiff[1]*boundSurf[25]+35.0*facDiff[2]*boundSurf[22]); 


  } else { 

  double edgeSurf[27] = {0.0}; 
  edgeSurf[0] = -0.0125*(53.66563145999496*fskin[9]-53.66563145999496*fedge[9]-95.26279441628824*(fskin[3]+fedge[3])+75.0*fskin[0]-75.0*fedge[0])*surfVar_l; 
  edgeSurf[1] = -0.0125*(53.66563145999495*fskin[15]-53.66563145999495*fedge[15]-95.26279441628824*(fskin[5]+fedge[5])+75.0*fskin[1]-75.0*fedge[1])*surfVar_l; 
  edgeSurf[2] = -0.0125*(53.66563145999495*fskin[16]-53.66563145999495*fedge[16]-95.26279441628824*(fskin[6]+fedge[6])+75.0*fskin[2]-75.0*fedge[2])*surfVar_l; 
  edgeSurf[3] = 0.003125*(507.3608183531716*fskin[9]-236.2519841186524*fedge[9]-855.0*fskin[3]-465.0*fedge[3]+658.1793068761733*fskin[0]-381.051177665153*fedge[0])*surfVar_l; 
  edgeSurf[4] = -0.0125*(53.66563145999496*fskin[19]-53.66563145999496*fedge[19]-95.26279441628824*(fskin[10]+fedge[10])+75.0*fskin[4]-75.0*fedge[4])*surfVar_l; 
  edgeSurf[5] = 0.003125*(507.3608183531716*fskin[15]-236.2519841186524*fedge[15]-855.0*fskin[5]-465.0*fedge[5]+658.1793068761733*fskin[1]-381.051177665153*fedge[1])*surfVar_l; 
  edgeSurf[6] = 0.003125*(507.3608183531716*fskin[16]-236.2519841186524*fedge[16]-855.0*fskin[6]-465.0*fedge[6]+658.1793068761733*fskin[2]-381.051177665153*fedge[2])*surfVar_l; 
  edgeSurf[7] = -0.0125*(53.66563145999496*fskin[21]-53.66563145999496*fedge[21]-95.26279441628826*(fskin[13]+fedge[13])+75.0*fskin[7]-75.0*fedge[7])*surfVar_l; 
  edgeSurf[8] = -0.0125*(53.66563145999496*fskin[22]-53.66563145999496*fedge[22]-95.26279441628826*(fskin[14]+fedge[14])+75.0*fskin[8]-75.0*fedge[8])*surfVar_l; 
  edgeSurf[9] = -0.015625*(201.0*fskin[9]+9.0*fedge[9]-321.4576177352156*fskin[3]-19.36491673103709*fedge[3]+241.4953415699773*fskin[0]-26.83281572999748*fedge[0])*surfVar_l; 
  edgeSurf[10] = 0.003125*(507.3608183531716*fskin[19]-236.2519841186524*fedge[19]-855.0*fskin[10]-465.0*fedge[10]+658.1793068761733*fskin[4]-381.051177665153*fedge[4])*surfVar_l; 
  edgeSurf[11] = -0.0125*(53.66563145999495*fskin[24]-53.66563145999495*fedge[24]-95.26279441628826*(fskin[17]+fedge[17])+75.0*fskin[11]-75.0*fedge[11])*surfVar_l; 
  edgeSurf[12] = -0.0125*(53.66563145999495*fskin[25]-53.66563145999495*fedge[25]-95.26279441628826*(fskin[18]+fedge[18])+75.0*fskin[12]-75.0*fedge[12])*surfVar_l; 
  edgeSurf[13] = 0.003125*(507.3608183531716*fskin[21]-236.2519841186524*fedge[21]-855.0*fskin[13]-465.0*fedge[13]+658.1793068761734*fskin[7]-381.051177665153*fedge[7])*surfVar_l; 
  edgeSurf[14] = 0.003125*(507.3608183531716*fskin[22]-236.2519841186524*fedge[22]-855.0*fskin[14]-465.0*fedge[14]+658.1793068761734*fskin[8]-381.051177665153*fedge[8])*surfVar_l; 
  edgeSurf[15] = -0.015625*(201.0*fskin[15]+9.0*fedge[15]-321.4576177352156*fskin[5]-19.36491673103708*fedge[5]+241.4953415699772*fskin[1]-26.83281572999747*fedge[1])*surfVar_l; 
  edgeSurf[16] = -0.015625*(201.0*fskin[16]+9.0*fedge[16]-321.4576177352156*fskin[6]-19.36491673103708*fedge[6]+241.4953415699772*fskin[2]-26.83281572999747*fedge[2])*surfVar_l; 
  edgeSurf[17] = 0.003125*(507.3608183531716*fskin[24]-236.2519841186524*fedge[24]-855.0*fskin[17]-465.0*fedge[17]+658.1793068761734*fskin[11]-381.051177665153*fedge[11])*surfVar_l; 
  edgeSurf[18] = 0.003125*(507.3608183531716*fskin[25]-236.2519841186524*fedge[25]-855.0*fskin[18]-465.0*fedge[18]+658.1793068761734*fskin[12]-381.051177665153*fedge[12])*surfVar_l; 
  edgeSurf[19] = -0.015625*(201.0*fskin[19]+9.0*fedge[19]-321.4576177352156*fskin[10]-19.36491673103709*fedge[10]+241.4953415699773*fskin[4]-26.83281572999748*fedge[4])*surfVar_l; 
  edgeSurf[20] = -0.0125*(53.66563145999496*fskin[26]-53.66563145999496*fedge[26]-95.26279441628824*(fskin[23]+fedge[23])+75.0*fskin[20]-75.0*fedge[20])*surfVar_l; 
  edgeSurf[21] = -0.015625*(201.0*fskin[21]+9.0*fedge[21]-321.4576177352156*fskin[13]-19.36491673103708*fedge[13]+241.4953415699773*fskin[7]-26.83281572999748*fedge[7])*surfVar_l; 
  edgeSurf[22] = -0.015625*(201.0*fskin[22]+9.0*fedge[22]-321.4576177352156*fskin[14]-19.36491673103708*fedge[14]+241.4953415699773*fskin[8]-26.83281572999748*fedge[8])*surfVar_l; 
  edgeSurf[23] = 0.003125*(507.3608183531716*fskin[26]-236.2519841186524*fedge[26]-855.0*fskin[23]-465.0*fedge[23]+658.1793068761733*fskin[20]-381.051177665153*fedge[20])*surfVar_l; 
  edgeSurf[24] = -0.015625*(201.0*fskin[24]+9.0*fedge[24]-321.4576177352156*fskin[17]-19.36491673103709*fedge[17]+241.4953415699772*fskin[11]-26.83281572999747*fedge[11])*surfVar_l; 
  edgeSurf[25] = -0.015625*(201.0*fskin[25]+9.0*fedge[25]-321.4576177352156*fskin[18]-19.36491673103709*fedge[18]+241.4953415699772*fskin[12]-26.83281572999747*fedge[12])*surfVar_l; 
  edgeSurf[26] = -0.015625*(201.0*fskin[26]+9.0*fedge[26]-321.4576177352156*fskin[23]-19.36491673103709*fedge[23]+241.4953415699773*fskin[20]-26.83281572999748*fedge[20])*surfVar_l; 

  double boundSurf[27] = {0.0}; 
  boundSurf[3] = -0.5*(3.872983346207417*fskin[9]+3.0*fskin[3]+1.732050807568877*fskin[0])*surfVar_r; 
  boundSurf[5] = -0.5*(3.872983346207417*fskin[15]+3.0*fskin[5]+1.732050807568877*fskin[1])*surfVar_r; 
  boundSurf[6] = -0.5*(3.872983346207417*fskin[16]+3.0*fskin[6]+1.732050807568877*fskin[2])*surfVar_r; 
  boundSurf[9] = -0.5*(15.0*fskin[9]+11.61895003862225*fskin[3]+6.708203932499369*fskin[0])*surfVar_r; 
  boundSurf[10] = -0.5*(3.872983346207417*fskin[19]+3.0*fskin[10]+1.732050807568877*fskin[4])*surfVar_r; 
  boundSurf[13] = -0.1*(19.36491673103708*fskin[21]+15.0*fskin[13]+8.660254037844387*fskin[7])*surfVar_r; 
  boundSurf[14] = -0.1*(19.36491673103708*fskin[22]+15.0*fskin[14]+8.660254037844387*fskin[8])*surfVar_r; 
  boundSurf[15] = -0.5*(15.0*fskin[15]+11.61895003862225*fskin[5]+6.708203932499369*fskin[1])*surfVar_r; 
  boundSurf[16] = -0.5*(15.0*fskin[16]+11.61895003862225*fskin[6]+6.708203932499369*fskin[2])*surfVar_r; 
  boundSurf[17] = -0.1*(19.36491673103709*fskin[24]+15.0*fskin[17]+8.660254037844387*fskin[11])*surfVar_r; 
  boundSurf[18] = -0.1*(19.36491673103709*fskin[25]+15.0*fskin[18]+8.660254037844387*fskin[12])*surfVar_r; 
  boundSurf[19] = -0.5*(15.0*fskin[19]+11.61895003862225*fskin[10]+6.708203932499369*fskin[4])*surfVar_r; 
  boundSurf[21] = -0.5*(15.0*fskin[21]+11.61895003862225*fskin[13]+6.708203932499369*fskin[7])*surfVar_r; 
  boundSurf[22] = -0.5*(15.0*fskin[22]+11.61895003862225*fskin[14]+6.708203932499369*fskin[8])*surfVar_r; 
  boundSurf[23] = -0.5*(3.872983346207417*fskin[26]+3.0*fskin[23]+1.732050807568877*fskin[20])*surfVar_r; 
  boundSurf[24] = -0.5*(15.0*fskin[24]+11.61895003862225*fskin[17]+6.708203932499369*fskin[11])*surfVar_r; 
  boundSurf[25] = -0.5*(15.0*fskin[25]+11.61895003862225*fskin[18]+6.708203932499369*fskin[12])*surfVar_r; 
  boundSurf[26] = -0.5*(15.0*fskin[26]+11.61895003862225*fskin[23]+6.708203932499369*fskin[20])*surfVar_r; 

  edgeSurf_incr[0] = 0.7071067811865475*(facDiff[2]*edgeSurf[7]+edgeSurf[1]*facDiff[1]+edgeSurf[0]*facDiff[0]); 
  edgeSurf_incr[1] = 0.1414213562373095*(4.47213595499958*(facDiff[1]*edgeSurf[7]+edgeSurf[1]*facDiff[2])+5.0*(edgeSurf[0]*facDiff[1]+facDiff[0]*edgeSurf[1])); 
  edgeSurf_incr[2] = 0.04714045207910316*(15.0*facDiff[2]*edgeSurf[11]+15.0*(facDiff[1]*edgeSurf[4]+facDiff[0]*edgeSurf[2])); 
  edgeSurf_incr[3] = 0.04714045207910316*(15.0*facDiff[2]*edgeSurf[13]+15.0*(facDiff[1]*edgeSurf[5]+facDiff[0]*edgeSurf[3])); 
  edgeSurf_incr[4] = 0.04714045207910316*(13.41640786499874*facDiff[1]*edgeSurf[11]+(13.41640786499874*facDiff[2]+15.0*facDiff[0])*edgeSurf[4]+15.0*facDiff[1]*edgeSurf[2]); 
  edgeSurf_incr[5] = 0.04714045207910316*(13.41640786499874*facDiff[1]*edgeSurf[13]+(13.41640786499874*facDiff[2]+15.0*facDiff[0])*edgeSurf[5]+15.0*facDiff[1]*edgeSurf[3]); 
  edgeSurf_incr[6] = 0.7071067811865475*(facDiff[2]*edgeSurf[17]+facDiff[1]*edgeSurf[10]+facDiff[0]*edgeSurf[6]); 
  edgeSurf_incr[7] = 0.02020305089104421*((22.3606797749979*facDiff[2]+35.0*facDiff[0])*edgeSurf[7]+35.0*edgeSurf[0]*facDiff[2]+31.30495168499706*edgeSurf[1]*facDiff[1]); 
  edgeSurf_incr[8] = 0.04714045207910316*(15.0*facDiff[2]*edgeSurf[20]+15.0*facDiff[1]*edgeSurf[12]+15.0*facDiff[0]*edgeSurf[8]); 
  edgeSurf_incr[9] = 0.04714045207910316*(15.0*facDiff[2]*edgeSurf[21]+15.0*facDiff[1]*edgeSurf[15]+15.0*facDiff[0]*edgeSurf[9]); 
  edgeSurf_incr[10] = 0.1414213562373095*(4.47213595499958*facDiff[1]*edgeSurf[17]+(4.47213595499958*facDiff[2]+5.0*facDiff[0])*edgeSurf[10]+5.0*facDiff[1]*edgeSurf[6]); 
  edgeSurf_incr[11] = 0.006734350297014738*((67.0820393249937*facDiff[2]+105.0*facDiff[0])*edgeSurf[11]+93.91485505499116*facDiff[1]*edgeSurf[4]+105.0*edgeSurf[2]*facDiff[2]); 
  edgeSurf_incr[12] = 0.04714045207910316*(13.41640786499874*facDiff[1]*edgeSurf[20]+(13.41640786499874*facDiff[2]+15.0*facDiff[0])*edgeSurf[12]+15.0*facDiff[1]*edgeSurf[8]); 
  edgeSurf_incr[13] = 0.006734350297014738*((67.0820393249937*facDiff[2]+105.0*facDiff[0])*edgeSurf[13]+93.91485505499116*facDiff[1]*edgeSurf[5]+105.0*facDiff[2]*edgeSurf[3]); 
  edgeSurf_incr[14] = 0.04714045207910316*(15.0*(facDiff[2]*edgeSurf[23]+facDiff[1]*edgeSurf[18])+15.0*facDiff[0]*edgeSurf[14]); 
  edgeSurf_incr[15] = 0.04714045207910316*(13.41640786499874*facDiff[1]*edgeSurf[21]+(13.41640786499874*facDiff[2]+15.0*facDiff[0])*edgeSurf[15]+15.0*facDiff[1]*edgeSurf[9]); 
  edgeSurf_incr[16] = 0.04714045207910316*(15.0*(facDiff[2]*edgeSurf[24]+facDiff[1]*edgeSurf[19])+15.0*facDiff[0]*edgeSurf[16]); 
  edgeSurf_incr[17] = 0.02020305089104421*((22.3606797749979*facDiff[2]+35.0*facDiff[0])*edgeSurf[17]+31.30495168499706*facDiff[1]*edgeSurf[10]+35.0*facDiff[2]*edgeSurf[6]); 
  edgeSurf_incr[18] = 0.04714045207910316*(13.41640786499874*facDiff[1]*edgeSurf[23]+(13.41640786499874*facDiff[2]+15.0*facDiff[0])*edgeSurf[18]+15.0*facDiff[1]*edgeSurf[14]); 
  edgeSurf_incr[19] = 0.04714045207910316*(13.41640786499874*facDiff[1]*edgeSurf[24]+(13.41640786499874*facDiff[2]+15.0*facDiff[0])*edgeSurf[19]+15.0*facDiff[1]*edgeSurf[16]); 
  edgeSurf_incr[20] = 0.006734350297014738*((67.0820393249937*facDiff[2]+105.0*facDiff[0])*edgeSurf[20]+93.91485505499116*facDiff[1]*edgeSurf[12]+105.0*facDiff[2]*edgeSurf[8]); 
  edgeSurf_incr[21] = 0.006734350297014738*((67.0820393249937*facDiff[2]+105.0*facDiff[0])*edgeSurf[21]+93.91485505499116*facDiff[1]*edgeSurf[15]+105.0*facDiff[2]*edgeSurf[9]); 
  edgeSurf_incr[22] = 0.7071067811865475*(facDiff[2]*edgeSurf[26]+facDiff[1]*edgeSurf[25]+facDiff[0]*edgeSurf[22]); 
  edgeSurf_incr[23] = 0.006734350297014738*((67.0820393249937*facDiff[2]+105.0*facDiff[0])*edgeSurf[23]+93.91485505499116*facDiff[1]*edgeSurf[18]+105.0*facDiff[2]*edgeSurf[14]); 
  edgeSurf_incr[24] = 0.006734350297014738*((67.0820393249937*facDiff[2]+105.0*facDiff[0])*edgeSurf[24]+93.91485505499116*facDiff[1]*edgeSurf[19]+105.0*facDiff[2]*edgeSurf[16]); 
  edgeSurf_incr[25] = 0.1414213562373095*(4.47213595499958*facDiff[1]*edgeSurf[26]+(4.47213595499958*facDiff[2]+5.0*facDiff[0])*edgeSurf[25]+5.0*facDiff[1]*edgeSurf[22]); 
  edgeSurf_incr[26] = 0.02020305089104421*((22.3606797749979*facDiff[2]+35.0*facDiff[0])*edgeSurf[26]+31.30495168499706*facDiff[1]*edgeSurf[25]+35.0*facDiff[2]*edgeSurf[22]); 

  boundSurf_incr[0] = 0.7071067811865475*(facDiff[2]*boundSurf[7]+boundSurf[1]*facDiff[1]+boundSurf[0]*facDiff[0]); 
  boundSurf_incr[1] = 0.1414213562373095*(4.47213595499958*(facDiff[1]*boundSurf[7]+boundSurf[1]*facDiff[2])+5.0*(boundSurf[0]*facDiff[1]+facDiff[0]*boundSurf[1])); 
  boundSurf_incr[2] = 0.04714045207910316*(15.0*facDiff[2]*boundSurf[11]+15.0*(facDiff[1]*boundSurf[4]+facDiff[0]*boundSurf[2])); 
  boundSurf_incr[3] = 0.04714045207910316*(15.0*facDiff[2]*boundSurf[13]+15.0*(facDiff[1]*boundSurf[5]+facDiff[0]*boundSurf[3])); 
  boundSurf_incr[4] = 0.04714045207910316*(13.41640786499874*facDiff[1]*boundSurf[11]+(13.41640786499874*facDiff[2]+15.0*facDiff[0])*boundSurf[4]+15.0*facDiff[1]*boundSurf[2]); 
  boundSurf_incr[5] = 0.04714045207910316*(13.41640786499874*facDiff[1]*boundSurf[13]+(13.41640786499874*facDiff[2]+15.0*facDiff[0])*boundSurf[5]+15.0*facDiff[1]*boundSurf[3]); 
  boundSurf_incr[6] = 0.7071067811865475*(facDiff[2]*boundSurf[17]+facDiff[1]*boundSurf[10]+facDiff[0]*boundSurf[6]); 
  boundSurf_incr[7] = 0.02020305089104421*((22.3606797749979*facDiff[2]+35.0*facDiff[0])*boundSurf[7]+35.0*boundSurf[0]*facDiff[2]+31.30495168499706*boundSurf[1]*facDiff[1]); 
  boundSurf_incr[8] = 0.04714045207910316*(15.0*facDiff[2]*boundSurf[20]+15.0*facDiff[1]*boundSurf[12]+15.0*facDiff[0]*boundSurf[8]); 
  boundSurf_incr[9] = 0.04714045207910316*(15.0*facDiff[2]*boundSurf[21]+15.0*facDiff[1]*boundSurf[15]+15.0*facDiff[0]*boundSurf[9]); 
  boundSurf_incr[10] = 0.1414213562373095*(4.47213595499958*facDiff[1]*boundSurf[17]+(4.47213595499958*facDiff[2]+5.0*facDiff[0])*boundSurf[10]+5.0*facDiff[1]*boundSurf[6]); 
  boundSurf_incr[11] = 0.006734350297014738*((67.0820393249937*facDiff[2]+105.0*facDiff[0])*boundSurf[11]+93.91485505499116*facDiff[1]*boundSurf[4]+105.0*boundSurf[2]*facDiff[2]); 
  boundSurf_incr[12] = 0.04714045207910316*(13.41640786499874*facDiff[1]*boundSurf[20]+(13.41640786499874*facDiff[2]+15.0*facDiff[0])*boundSurf[12]+15.0*facDiff[1]*boundSurf[8]); 
  boundSurf_incr[13] = 0.006734350297014738*((67.0820393249937*facDiff[2]+105.0*facDiff[0])*boundSurf[13]+93.91485505499116*facDiff[1]*boundSurf[5]+105.0*facDiff[2]*boundSurf[3]); 
  boundSurf_incr[14] = 0.04714045207910316*(15.0*(facDiff[2]*boundSurf[23]+facDiff[1]*boundSurf[18])+15.0*facDiff[0]*boundSurf[14]); 
  boundSurf_incr[15] = 0.04714045207910316*(13.41640786499874*facDiff[1]*boundSurf[21]+(13.41640786499874*facDiff[2]+15.0*facDiff[0])*boundSurf[15]+15.0*facDiff[1]*boundSurf[9]); 
  boundSurf_incr[16] = 0.04714045207910316*(15.0*(facDiff[2]*boundSurf[24]+facDiff[1]*boundSurf[19])+15.0*facDiff[0]*boundSurf[16]); 
  boundSurf_incr[17] = 0.02020305089104421*((22.3606797749979*facDiff[2]+35.0*facDiff[0])*boundSurf[17]+31.30495168499706*facDiff[1]*boundSurf[10]+35.0*facDiff[2]*boundSurf[6]); 
  boundSurf_incr[18] = 0.04714045207910316*(13.41640786499874*facDiff[1]*boundSurf[23]+(13.41640786499874*facDiff[2]+15.0*facDiff[0])*boundSurf[18]+15.0*facDiff[1]*boundSurf[14]); 
  boundSurf_incr[19] = 0.04714045207910316*(13.41640786499874*facDiff[1]*boundSurf[24]+(13.41640786499874*facDiff[2]+15.0*facDiff[0])*boundSurf[19]+15.0*facDiff[1]*boundSurf[16]); 
  boundSurf_incr[20] = 0.006734350297014738*((67.0820393249937*facDiff[2]+105.0*facDiff[0])*boundSurf[20]+93.91485505499116*facDiff[1]*boundSurf[12]+105.0*facDiff[2]*boundSurf[8]); 
  boundSurf_incr[21] = 0.006734350297014738*((67.0820393249937*facDiff[2]+105.0*facDiff[0])*boundSurf[21]+93.91485505499116*facDiff[1]*boundSurf[15]+105.0*facDiff[2]*boundSurf[9]); 
  boundSurf_incr[22] = 0.7071067811865475*(facDiff[2]*boundSurf[26]+facDiff[1]*boundSurf[25]+facDiff[0]*boundSurf[22]); 
  boundSurf_incr[23] = 0.006734350297014738*((67.0820393249937*facDiff[2]+105.0*facDiff[0])*boundSurf[23]+93.91485505499116*facDiff[1]*boundSurf[18]+105.0*facDiff[2]*boundSurf[14]); 
  boundSurf_incr[24] = 0.006734350297014738*((67.0820393249937*facDiff[2]+105.0*facDiff[0])*boundSurf[24]+93.91485505499116*facDiff[1]*boundSurf[19]+105.0*facDiff[2]*boundSurf[16]); 
  boundSurf_incr[25] = 0.1414213562373095*(4.47213595499958*facDiff[1]*boundSurf[26]+(4.47213595499958*facDiff[2]+5.0*facDiff[0])*boundSurf[25]+5.0*facDiff[1]*boundSurf[22]); 
  boundSurf_incr[26] = 0.02020305089104421*((22.3606797749979*facDiff[2]+35.0*facDiff[0])*boundSurf[26]+31.30495168499706*facDiff[1]*boundSurf[25]+35.0*facDiff[2]*boundSurf[22]); 

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
  out[9] += (vol_incr[9]+edgeSurf_incr[9]+boundSurf_incr[9])*rdvSq4; 
  out[10] += (vol_incr[10]+edgeSurf_incr[10]+boundSurf_incr[10])*rdvSq4; 
  out[11] += (vol_incr[11]+edgeSurf_incr[11]+boundSurf_incr[11])*rdvSq4; 
  out[12] += (vol_incr[12]+edgeSurf_incr[12]+boundSurf_incr[12])*rdvSq4; 
  out[13] += (vol_incr[13]+edgeSurf_incr[13]+boundSurf_incr[13])*rdvSq4; 
  out[14] += (vol_incr[14]+edgeSurf_incr[14]+boundSurf_incr[14])*rdvSq4; 
  out[15] += (vol_incr[15]+edgeSurf_incr[15]+boundSurf_incr[15])*rdvSq4; 
  out[16] += (vol_incr[16]+edgeSurf_incr[16]+boundSurf_incr[16])*rdvSq4; 
  out[17] += (vol_incr[17]+edgeSurf_incr[17]+boundSurf_incr[17])*rdvSq4; 
  out[18] += (vol_incr[18]+edgeSurf_incr[18]+boundSurf_incr[18])*rdvSq4; 
  out[19] += (vol_incr[19]+edgeSurf_incr[19]+boundSurf_incr[19])*rdvSq4; 
  out[20] += (vol_incr[20]+edgeSurf_incr[20]+boundSurf_incr[20])*rdvSq4; 
  out[21] += (vol_incr[21]+edgeSurf_incr[21]+boundSurf_incr[21])*rdvSq4; 
  out[22] += (vol_incr[22]+edgeSurf_incr[22]+boundSurf_incr[22])*rdvSq4; 
  out[23] += (vol_incr[23]+edgeSurf_incr[23]+boundSurf_incr[23])*rdvSq4; 
  out[24] += (vol_incr[24]+edgeSurf_incr[24]+boundSurf_incr[24])*rdvSq4; 
  out[25] += (vol_incr[25]+edgeSurf_incr[25]+boundSurf_incr[25])*rdvSq4; 
  out[26] += (vol_incr[26]+edgeSurf_incr[26]+boundSurf_incr[26])*rdvSq4; 
} 
