#include <gkyl_vlasov_lbo_mom_kernels.h> 
GKYL_CU_DH void vlasov_boundary_integral_1x3v_ser_f_p1(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // edge:   indicator of which velocity grid edge is being considered. 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:   cell length in each direction. 
  // fIn[NP]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  double dS = 0.0; 
 
  if (edge == GKYL_VX_LOWER) {
 
    dS = 0.25*dxv[2]*dxv[3]; 
    out[0] += (2.449489742783178*fIn[2]-1.414213562373095*fIn[0])*dS; 
    out[1] += (2.449489742783178*fIn[5]-1.414213562373095*fIn[1])*dS; 
 
  } else if (edge == GKYL_VX_UPPER) {
 
    dS = 0.25*dxv[2]*dxv[3]; 
    out[0] += (2.449489742783178*fIn[2]+1.414213562373095*fIn[0])*dS; 
    out[1] += (2.449489742783178*fIn[5]+1.414213562373095*fIn[1])*dS; 
 
  } 
 
  if (edge == GKYL_VY_LOWER) {
 
    dS = 0.25*dxv[1]*dxv[3]; 
    out[2] += (2.449489742783178*fIn[3]-1.414213562373095*fIn[0])*dS; 
    out[3] += (2.449489742783178*fIn[6]-1.414213562373095*fIn[1])*dS; 
 
  } else if (edge == GKYL_VY_UPPER) {
 
    dS = 0.25*dxv[1]*dxv[3]; 
    out[2] += (2.449489742783178*fIn[3]+1.414213562373095*fIn[0])*dS; 
    out[3] += (2.449489742783178*fIn[6]+1.414213562373095*fIn[1])*dS; 
 
  } 
 
  if (edge == GKYL_VZ_LOWER) {
 
    dS = 0.25*dxv[1]*dxv[2]; 
    out[4] += (2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS; 
    out[5] += (2.449489742783178*fIn[8]-1.414213562373095*fIn[1])*dS; 
 
  } else if (edge == GKYL_VZ_UPPER) {
 
    dS = 0.25*dxv[1]*dxv[2]; 
    out[4] += (2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS; 
    out[5] += (2.449489742783178*fIn[8]+1.414213562373095*fIn[1])*dS; 
 
  } 
 
}
 
GKYL_CU_DH void vlasov_boundary_integral_1x3v_ser_vf_p1(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // edge:   indicator of which velocity grid edge is being considered. 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:   cell length in each direction. 
  // fIn[NP]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  double dS = 0.0; 
 
  if (edge == GKYL_VX_LOWER) {
 
    dS = 0.25*dxv[2]*dxv[3]; 
      out[0] += (1.224744871391589*dxv[1]*fIn[2]-0.7071067811865475*fIn[0]*dxv[1])*dS+vBoundary[0]*(2.449489742783178*fIn[2]-1.414213562373095*fIn[0])*dS; 
      out[1] += (1.224744871391589*dxv[1]*fIn[5]-0.7071067811865475*dxv[1]*fIn[1])*dS+vBoundary[0]*(2.449489742783178*fIn[5]-1.414213562373095*fIn[1])*dS; 
 
  } else if (edge == GKYL_VX_UPPER) {
 
    dS = 0.25*dxv[2]*dxv[3]; 
      out[0] += (2.449489742783178*fIn[2]+1.414213562373095*fIn[0])*vBoundary[3]*dS+((-1.224744871391589*dxv[1]*fIn[2])-0.7071067811865475*fIn[0]*dxv[1])*dS; 
      out[1] += ((-1.224744871391589*dxv[1]*fIn[5])-0.7071067811865475*dxv[1]*fIn[1])*dS+vBoundary[3]*(2.449489742783178*fIn[5]+1.414213562373095*fIn[1])*dS; 
 
  } 
 
  if (edge == GKYL_VY_LOWER) {
 
    dS = 0.25*dxv[1]*dxv[3]; 
      out[0] += (1.224744871391589*dxv[2]*fIn[3]-0.7071067811865475*fIn[0]*dxv[2])*dS+vBoundary[1]*(2.449489742783178*fIn[3]-1.414213562373095*fIn[0])*dS; 
      out[1] += (1.224744871391589*dxv[2]*fIn[6]-0.7071067811865475*fIn[1]*dxv[2])*dS+vBoundary[1]*(2.449489742783178*fIn[6]-1.414213562373095*fIn[1])*dS; 
 
  } else if (edge == GKYL_VY_UPPER) {
 
    dS = 0.25*dxv[1]*dxv[3]; 
      out[0] += (2.449489742783178*fIn[3]+1.414213562373095*fIn[0])*vBoundary[4]*dS+((-1.224744871391589*dxv[2]*fIn[3])-0.7071067811865475*fIn[0]*dxv[2])*dS; 
      out[1] += ((-1.224744871391589*dxv[2]*fIn[6])-0.7071067811865475*fIn[1]*dxv[2])*dS+vBoundary[4]*(2.449489742783178*fIn[6]+1.414213562373095*fIn[1])*dS; 
 
  } 
 
  if (edge == GKYL_VZ_LOWER) {
 
    dS = 0.25*dxv[1]*dxv[2]; 
      out[0] += (1.224744871391589*dxv[3]*fIn[4]-0.7071067811865475*fIn[0]*dxv[3])*dS+vBoundary[2]*(2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS; 
      out[1] += (1.224744871391589*dxv[3]*fIn[8]-0.7071067811865475*fIn[1]*dxv[3])*dS+vBoundary[2]*(2.449489742783178*fIn[8]-1.414213562373095*fIn[1])*dS; 
 
  } else if (edge == GKYL_VZ_UPPER) {
 
    dS = 0.25*dxv[1]*dxv[2]; 
      out[0] += (2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*vBoundary[5]*dS+((-1.224744871391589*dxv[3]*fIn[4])-0.7071067811865475*fIn[0]*dxv[3])*dS; 
      out[1] += ((-1.224744871391589*dxv[3]*fIn[8])-0.7071067811865475*fIn[1]*dxv[3])*dS+vBoundary[5]*(2.449489742783178*fIn[8]+1.414213562373095*fIn[1])*dS; 
 
  } 
 
}
 
GKYL_CU_DH void vlasov_boundary_integral_1x3v_ser_f_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // edge:   indicator of which velocity grid edge is being considered. 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:   cell length in each direction. 
  // fIn[NP]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  double dS = 0.0; 
 
  if (edge == GKYL_VX_LOWER) {
 
    dS = 0.25*dxv[2]*dxv[3]; 
    out[0] += ((-3.16227766016838*fIn[12])+2.449489742783178*fIn[2]-1.414213562373095*fIn[0])*dS; 
    out[1] += ((-3.16227766016838*fIn[20])+2.449489742783178*fIn[5]-1.414213562373095*fIn[1])*dS; 
    out[2] += (2.449489742783178*fIn[19]-1.414213562373095*fIn[11])*dS; 
 
  } else if (edge == GKYL_VX_UPPER) {
 
    dS = 0.25*dxv[2]*dxv[3]; 
    out[0] += (3.16227766016838*fIn[12]+2.449489742783178*fIn[2]+1.414213562373095*fIn[0])*dS; 
    out[1] += (3.16227766016838*fIn[20]+2.449489742783178*fIn[5]+1.414213562373095*fIn[1])*dS; 
    out[2] += (2.449489742783178*fIn[19]+1.414213562373095*fIn[11])*dS; 
 
  } 
 
  if (edge == GKYL_VY_LOWER) {
 
    dS = 0.25*dxv[1]*dxv[3]; 
    out[3] += ((-3.16227766016838*fIn[13])+2.449489742783178*fIn[3]-1.414213562373095*fIn[0])*dS; 
    out[4] += ((-3.16227766016838*fIn[23])+2.449489742783178*fIn[6]-1.414213562373095*fIn[1])*dS; 
    out[5] += (2.449489742783178*fIn[21]-1.414213562373095*fIn[11])*dS; 
 
  } else if (edge == GKYL_VY_UPPER) {
 
    dS = 0.25*dxv[1]*dxv[3]; 
    out[3] += (3.16227766016838*fIn[13]+2.449489742783178*fIn[3]+1.414213562373095*fIn[0])*dS; 
    out[4] += (3.16227766016838*fIn[23]+2.449489742783178*fIn[6]+1.414213562373095*fIn[1])*dS; 
    out[5] += (2.449489742783178*fIn[21]+1.414213562373095*fIn[11])*dS; 
 
  } 
 
  if (edge == GKYL_VZ_LOWER) {
 
    dS = 0.25*dxv[1]*dxv[2]; 
    out[6] += ((-3.16227766016838*fIn[14])+2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS; 
    out[7] += ((-3.16227766016838*fIn[28])+2.449489742783178*fIn[8]-1.414213562373095*fIn[1])*dS; 
    out[8] += (2.449489742783178*fIn[25]-1.414213562373095*fIn[11])*dS; 
 
  } else if (edge == GKYL_VZ_UPPER) {
 
    dS = 0.25*dxv[1]*dxv[2]; 
    out[6] += (3.16227766016838*fIn[14]+2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS; 
    out[7] += (3.16227766016838*fIn[28]+2.449489742783178*fIn[8]+1.414213562373095*fIn[1])*dS; 
    out[8] += (2.449489742783178*fIn[25]+1.414213562373095*fIn[11])*dS; 
 
  } 
 
}
 
GKYL_CU_DH void vlasov_boundary_integral_1x3v_ser_vf_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // edge:   indicator of which velocity grid edge is being considered. 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:   cell length in each direction. 
  // fIn[NP]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  double dS = 0.0; 
 
  if (edge == GKYL_VX_LOWER) {
 
    dS = 0.25*dxv[2]*dxv[3]; 
      out[0] += vBoundary[0]*((-3.16227766016838*fIn[12])+2.449489742783178*fIn[2]-1.414213562373095*fIn[0])*dS; 
      out[1] += vBoundary[0]*((-3.16227766016838*fIn[20])+2.449489742783178*fIn[5]-1.414213562373095*fIn[1])*dS; 
      out[2] += vBoundary[0]*(2.449489742783178*fIn[19]-1.414213562373095*fIn[11])*dS; 
 
  } else if (edge == GKYL_VX_UPPER) {
 
    dS = 0.25*dxv[2]*dxv[3]; 
      out[0] += vBoundary[3]*(3.16227766016838*fIn[12]+2.449489742783178*fIn[2]+1.414213562373095*fIn[0])*dS; 
      out[1] += vBoundary[3]*(3.16227766016838*fIn[20]+2.449489742783178*fIn[5]+1.414213562373095*fIn[1])*dS; 
      out[2] += vBoundary[3]*(2.449489742783178*fIn[19]+1.414213562373095*fIn[11])*dS; 
 
  } 
 
  if (edge == GKYL_VY_LOWER) {
 
    dS = 0.25*dxv[1]*dxv[3]; 
      out[0] += vBoundary[1]*((-3.16227766016838*fIn[13])+2.449489742783178*fIn[3]-1.414213562373095*fIn[0])*dS; 
      out[1] += vBoundary[1]*((-3.16227766016838*fIn[23])+2.449489742783178*fIn[6]-1.414213562373095*fIn[1])*dS; 
      out[2] += vBoundary[1]*(2.449489742783178*fIn[21]-1.414213562373095*fIn[11])*dS; 
 
  } else if (edge == GKYL_VY_UPPER) {
 
    dS = 0.25*dxv[1]*dxv[3]; 
      out[0] += vBoundary[4]*(3.16227766016838*fIn[13]+2.449489742783178*fIn[3]+1.414213562373095*fIn[0])*dS; 
      out[1] += vBoundary[4]*(3.16227766016838*fIn[23]+2.449489742783178*fIn[6]+1.414213562373095*fIn[1])*dS; 
      out[2] += vBoundary[4]*(2.449489742783178*fIn[21]+1.414213562373095*fIn[11])*dS; 
 
  } 
 
  if (edge == GKYL_VZ_LOWER) {
 
    dS = 0.25*dxv[1]*dxv[2]; 
      out[0] += vBoundary[2]*((-3.16227766016838*fIn[14])+2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS; 
      out[1] += vBoundary[2]*((-3.16227766016838*fIn[28])+2.449489742783178*fIn[8]-1.414213562373095*fIn[1])*dS; 
      out[2] += vBoundary[2]*(2.449489742783178*fIn[25]-1.414213562373095*fIn[11])*dS; 
 
  } else if (edge == GKYL_VZ_UPPER) {
 
    dS = 0.25*dxv[1]*dxv[2]; 
      out[0] += vBoundary[5]*(3.16227766016838*fIn[14]+2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS; 
      out[1] += vBoundary[5]*(3.16227766016838*fIn[28]+2.449489742783178*fIn[8]+1.414213562373095*fIn[1])*dS; 
      out[2] += vBoundary[5]*(2.449489742783178*fIn[25]+1.414213562373095*fIn[11])*dS; 
 
  } 
 
}
 
