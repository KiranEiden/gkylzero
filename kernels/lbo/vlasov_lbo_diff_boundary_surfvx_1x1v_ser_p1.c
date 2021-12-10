#include <gkyl_vlasov_lbo_kernels.h> 
GKYL_CU_DH void vlasov_lbo_diff_boundary_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[2]:         Cell-center coordinates. 
  // dxv[2]:       Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[2]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]: sum of thermal speeds squared time their respective collisionalities. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 

  double diff_incr[4] = {0.0}; 

  if (edge == -1) { 

  diff_incr[0] = 0.2296396633859228*nuVtSqSum[1]*fSkin[3]-0.3827327723098713*nuVtSqSum[1]*fEdge[3]+0.2296396633859228*nuVtSqSum[0]*fSkin[2]-0.3827327723098713*nuVtSqSum[0]*fEdge[2]-0.7513009550107064*fSkin[1]*nuVtSqSum[1]+0.3977475644174328*fEdge[1]*nuVtSqSum[1]-0.7513009550107064*fSkin[0]*nuVtSqSum[0]+0.3977475644174328*fEdge[0]*nuVtSqSum[0]; 
  diff_incr[1] = 0.2296396633859228*nuVtSqSum[0]*fSkin[3]-0.3827327723098713*nuVtSqSum[0]*fEdge[3]+0.2296396633859228*nuVtSqSum[1]*fSkin[2]-0.3827327723098713*nuVtSqSum[1]*fEdge[2]-0.7513009550107064*fSkin[0]*nuVtSqSum[1]+0.3977475644174328*fEdge[0]*nuVtSqSum[1]-0.7513009550107064*nuVtSqSum[0]*fSkin[1]+0.3977475644174328*nuVtSqSum[0]*fEdge[1]; 
  diff_incr[2] = (-2.077126169735482*nuVtSqSum[1]*fSkin[3])-0.3093592167691144*nuVtSqSum[1]*fEdge[3]-2.077126169735482*nuVtSqSum[0]*fSkin[2]-0.3093592167691144*nuVtSqSum[0]*fEdge[2]-0.3827327723098713*fSkin[1]*nuVtSqSum[1]+0.3827327723098713*fEdge[1]*nuVtSqSum[1]-0.3827327723098713*fSkin[0]*nuVtSqSum[0]+0.3827327723098713*fEdge[0]*nuVtSqSum[0]; 
  diff_incr[3] = (-2.077126169735482*nuVtSqSum[0]*fSkin[3])-0.3093592167691144*nuVtSqSum[0]*fEdge[3]-2.077126169735482*nuVtSqSum[1]*fSkin[2]-0.3093592167691144*nuVtSqSum[1]*fEdge[2]-0.3827327723098713*fSkin[0]*nuVtSqSum[1]+0.3827327723098713*fEdge[0]*nuVtSqSum[1]-0.3827327723098713*nuVtSqSum[0]*fSkin[1]+0.3827327723098713*nuVtSqSum[0]*fEdge[1]; 


  } else { 

  diff_incr[0] = 0.9951052080056654*nuVtSqSum[1]*fSkin[3]+0.3827327723098713*nuVtSqSum[1]*fEdge[3]+0.9951052080056654*nuVtSqSum[0]*fSkin[2]+0.3827327723098713*nuVtSqSum[0]*fEdge[2]-0.0441941738241592*fSkin[1]*nuVtSqSum[1]+0.3977475644174328*fEdge[1]*nuVtSqSum[1]-0.0441941738241592*fSkin[0]*nuVtSqSum[0]+0.3977475644174328*fEdge[0]*nuVtSqSum[0]; 
  diff_incr[1] = 0.9951052080056654*nuVtSqSum[0]*fSkin[3]+0.3827327723098713*nuVtSqSum[0]*fEdge[3]+0.9951052080056654*nuVtSqSum[1]*fSkin[2]+0.3827327723098713*nuVtSqSum[1]*fEdge[2]-0.0441941738241592*fSkin[0]*nuVtSqSum[1]+0.3977475644174328*fEdge[0]*nuVtSqSum[1]-0.0441941738241592*nuVtSqSum[0]*fSkin[1]+0.3977475644174328*nuVtSqSum[0]*fEdge[1]; 
  diff_incr[2] = 0.0441941738241592*nuVtSqSum[1]*fSkin[3]-0.3093592167691144*nuVtSqSum[1]*fEdge[3]+0.0441941738241592*nuVtSqSum[0]*fSkin[2]-0.3093592167691144*nuVtSqSum[0]*fEdge[2]+1.60747764370146*fSkin[1]*nuVtSqSum[1]-0.3827327723098713*fEdge[1]*nuVtSqSum[1]+1.60747764370146*fSkin[0]*nuVtSqSum[0]-0.3827327723098713*fEdge[0]*nuVtSqSum[0]; 
  diff_incr[3] = 0.0441941738241592*nuVtSqSum[0]*fSkin[3]-0.3093592167691144*nuVtSqSum[0]*fEdge[3]+0.0441941738241592*nuVtSqSum[1]*fSkin[2]-0.3093592167691144*nuVtSqSum[1]*fEdge[2]+1.60747764370146*fSkin[0]*nuVtSqSum[1]-0.3827327723098713*fEdge[0]*nuVtSqSum[1]+1.60747764370146*nuVtSqSum[0]*fSkin[1]-0.3827327723098713*nuVtSqSum[0]*fEdge[1]; 

  } 

  out[0] += diff_incr[0]*rdvSq4; 
  out[1] += diff_incr[1]*rdvSq4; 
  out[2] += diff_incr[2]*rdvSq4; 
  out[3] += diff_incr[3]*rdvSq4; 
} 
