#include <gkyl_mom_bcorr_lbo_vlasov_kernels.h> 
GKYL_CU_DH void mom_bcorr_lbo_vlasov_1x3v_ser_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // edge:      indicator of which velocity grid edge is being considered (VX/VPAR, VY/MU, VZ). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:   cell length in each direction. 
  // fIn[48]:   distribution function at lower/upper velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) [first vdim components] and vf^(vmax)_(vmin) [last component]. 
 
  double dS = 0.0; 
 
  if (edge == GKYL_VX_LOWER) {
 
    dS = 0.25*dxv[2]*dxv[3]; 
    out[0] += ((-3.16227766016838*fIn[12])+2.449489742783178*fIn[2]-1.414213562373095*fIn[0])*dS; 
    out[1] += ((-3.16227766016838*fIn[20])+2.449489742783178*fIn[5]-1.414213562373095*fIn[1])*dS; 
    out[2] += (2.449489742783178*fIn[19]-1.414213562373095*fIn[11])*dS; 
    out[9] += vBoundary[0]*((-3.16227766016838*fIn[12])+2.449489742783178*fIn[2]-1.414213562373095*fIn[0])*dS; 
    out[10] += vBoundary[0]*((-3.16227766016838*fIn[20])+2.449489742783178*fIn[5]-1.414213562373095*fIn[1])*dS; 
    out[11] += vBoundary[0]*(2.449489742783178*fIn[19]-1.414213562373095*fIn[11])*dS; 
 
  } else if (edge == GKYL_VX_UPPER) {
 
    dS = 0.25*dxv[2]*dxv[3]; 
    out[0] += (3.16227766016838*fIn[12]+2.449489742783178*fIn[2]+1.414213562373095*fIn[0])*dS; 
    out[1] += (3.16227766016838*fIn[20]+2.449489742783178*fIn[5]+1.414213562373095*fIn[1])*dS; 
    out[2] += (2.449489742783178*fIn[19]+1.414213562373095*fIn[11])*dS; 
    out[9] += vBoundary[3]*(3.16227766016838*fIn[12]+2.449489742783178*fIn[2]+1.414213562373095*fIn[0])*dS; 
    out[10] += vBoundary[3]*(3.16227766016838*fIn[20]+2.449489742783178*fIn[5]+1.414213562373095*fIn[1])*dS; 
    out[11] += vBoundary[3]*(2.449489742783178*fIn[19]+1.414213562373095*fIn[11])*dS; 
  }
 
  if (edge == GKYL_VY_LOWER) {
 
    dS = 0.25*dxv[1]*dxv[3]; 
    out[3] += ((-3.16227766016838*fIn[13])+2.449489742783178*fIn[3]-1.414213562373095*fIn[0])*dS; 
    out[4] += ((-3.16227766016838*fIn[23])+2.449489742783178*fIn[6]-1.414213562373095*fIn[1])*dS; 
    out[5] += (2.449489742783178*fIn[21]-1.414213562373095*fIn[11])*dS; 
    out[9] += vBoundary[1]*((-3.16227766016838*fIn[13])+2.449489742783178*fIn[3]-1.414213562373095*fIn[0])*dS; 
    out[10] += vBoundary[1]*((-3.16227766016838*fIn[23])+2.449489742783178*fIn[6]-1.414213562373095*fIn[1])*dS; 
    out[11] += vBoundary[1]*(2.449489742783178*fIn[21]-1.414213562373095*fIn[11])*dS; 
 
  } else if (edge == GKYL_VY_UPPER) {
 
    dS = 0.25*dxv[1]*dxv[3]; 
    out[3] += (3.16227766016838*fIn[13]+2.449489742783178*fIn[3]+1.414213562373095*fIn[0])*dS; 
    out[4] += (3.16227766016838*fIn[23]+2.449489742783178*fIn[6]+1.414213562373095*fIn[1])*dS; 
    out[5] += (2.449489742783178*fIn[21]+1.414213562373095*fIn[11])*dS; 
    out[9] += vBoundary[4]*(3.16227766016838*fIn[13]+2.449489742783178*fIn[3]+1.414213562373095*fIn[0])*dS; 
    out[10] += vBoundary[4]*(3.16227766016838*fIn[23]+2.449489742783178*fIn[6]+1.414213562373095*fIn[1])*dS; 
    out[11] += vBoundary[4]*(2.449489742783178*fIn[21]+1.414213562373095*fIn[11])*dS; 
  }
 
  if (edge == GKYL_VZ_LOWER) {
 
    dS = 0.25*dxv[1]*dxv[2]; 
    out[6] += ((-3.16227766016838*fIn[14])+2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS; 
    out[7] += ((-3.16227766016838*fIn[28])+2.449489742783178*fIn[8]-1.414213562373095*fIn[1])*dS; 
    out[8] += (2.449489742783178*fIn[25]-1.414213562373095*fIn[11])*dS; 
    out[9] += vBoundary[2]*((-3.16227766016838*fIn[14])+2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS; 
    out[10] += vBoundary[2]*((-3.16227766016838*fIn[28])+2.449489742783178*fIn[8]-1.414213562373095*fIn[1])*dS; 
    out[11] += vBoundary[2]*(2.449489742783178*fIn[25]-1.414213562373095*fIn[11])*dS; 
 
  } else if (edge == GKYL_VZ_UPPER) {
 
    dS = 0.25*dxv[1]*dxv[2]; 
    out[6] += (3.16227766016838*fIn[14]+2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS; 
    out[7] += (3.16227766016838*fIn[28]+2.449489742783178*fIn[8]+1.414213562373095*fIn[1])*dS; 
    out[8] += (2.449489742783178*fIn[25]+1.414213562373095*fIn[11])*dS; 
    out[9] += vBoundary[5]*(3.16227766016838*fIn[14]+2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS; 
    out[10] += vBoundary[5]*(3.16227766016838*fIn[28]+2.449489742783178*fIn[8]+1.414213562373095*fIn[1])*dS; 
    out[11] += vBoundary[5]*(2.449489742783178*fIn[25]+1.414213562373095*fIn[11])*dS; 
  }
 
} 
 
