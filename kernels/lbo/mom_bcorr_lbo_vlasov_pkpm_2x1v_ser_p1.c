#include <gkyl_mom_bcorr_lbo_vlasov_kernels.h> 
GKYL_CU_DH void mom_bcorr_lbo_vlasov_pkpm_2x1v_ser_p1(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, double mass, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // edge:      indicator of which velocity grid edge is being considered (VX/VPAR). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:    cell length in each direction. 
  // mass:      species mass. 
  // fIn[12]:    distribution function at lower/upper velocity boundaries. 
  // out:       int dS of vf^(vmax)_(vmin). 
 
  if (edge == GKYL_VX_LOWER) {
    out[0] += vBoundary[0]*((-1.58113883008419*fIn[8])+1.224744871391589*fIn[3]-0.7071067811865475*fIn[0])*mass; 
    out[1] += vBoundary[0]*((-1.58113883008419*fIn[9])+1.224744871391589*fIn[5]-0.7071067811865475*fIn[1])*mass; 
    out[2] += vBoundary[0]*((-1.58113883008419*fIn[10])+1.224744871391589*fIn[6]-0.7071067811865475*fIn[2])*mass; 
    out[3] += vBoundary[0]*((-1.58113883008419*fIn[11])+1.224744871391589*fIn[7]-0.7071067811865475*fIn[4])*mass; 
  } 
  else if (edge == GKYL_VX_UPPER) {
    out[0] += vBoundary[1]*(1.58113883008419*fIn[8]+1.224744871391589*fIn[3]+0.7071067811865475*fIn[0])*mass; 
    out[1] += vBoundary[1]*(1.58113883008419*fIn[9]+1.224744871391589*fIn[5]+0.7071067811865475*fIn[1])*mass; 
    out[2] += vBoundary[1]*(1.58113883008419*fIn[10]+1.224744871391589*fIn[6]+0.7071067811865475*fIn[2])*mass; 
    out[3] += vBoundary[1]*(1.58113883008419*fIn[11]+1.224744871391589*fIn[7]+0.7071067811865475*fIn[4])*mass; 
  }
 
} 
 
