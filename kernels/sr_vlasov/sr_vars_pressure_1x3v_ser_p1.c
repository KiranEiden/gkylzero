#include <gkyl_sr_Gamma_kernels.h> 
GKYL_CU_DH void sr_vars_pressure_1x3v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure) 
{ 
  // gamma:       Particle Lorentz boost factor sqrt(1 + p^2).
  // gamma_inv:   Inverse particle Lorentz boost factor 1/sqrt(1 + p^2).
  // u_i:         Spatial components of bulk four-velocity = GammaV*V_drift. 
  // u_i_sq:      Squared spatial components of bulk four-velocity = u_i^2. 
  // GammaV:      Bulk four-velocity Lorentz factor = sqrt(1 + |u_i|^2). 
  // GammaV_sq:   Squared bulk four-velocity Lorentz factor = 1 + |u_i|^2. 
  // f:           Input distribution function.
  // sr_pressure: Output relativistic pressure.
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double *V_0 = &u_i[0]; 
  const double *V_0_sq = &u_i_sq[0]; 
 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double *V_1 = &u_i[2]; 
  const double *V_1_sq = &u_i_sq[2]; 
 
  const double wx3 = w[3], dv3 = dxv[3]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  const double *V_2 = &u_i[4]; 
  const double *V_2_sq = &u_i_sq[4]; 
 
  double temp[40] = {0.0}; 
  double temp_sq[40] = {0.0}; 
  double p_fac[40] = {0.0}; 
  temp[0] = 2.828427124746191*V_2[0]*wx3+2.828427124746191*V_1[0]*wx2+2.828427124746191*V_0[0]*wx1; 
  temp[1] = 2.828427124746191*V_2[1]*wx3+2.828427124746191*V_1[1]*wx2+2.828427124746191*V_0[1]*wx1; 
  temp[2] = 0.8164965809277261*V_0[0]*dv1; 
  temp[3] = 0.8164965809277261*V_1[0]*dv2; 
  temp[4] = 0.8164965809277261*V_2[0]*dv3; 
  temp[5] = 0.8164965809277261*V_0[1]*dv1; 
  temp[6] = 0.8164965809277261*V_1[1]*dv2; 
  temp[8] = 0.8164965809277261*V_2[1]*dv3; 

  temp_sq[0] = 2.828427124746191*V_2_sq[0]*wx3_sq+4.0*V_1[1]*V_2[1]*wx2*wx3+4.0*V_1[0]*V_2[0]*wx2*wx3+4.0*V_0[1]*V_2[1]*wx1*wx3+4.0*V_0[0]*V_2[0]*wx1*wx3+2.828427124746191*V_1_sq[0]*wx2_sq+4.0*V_0[1]*V_1[1]*wx1*wx2+4.0*V_0[0]*V_1[0]*wx1*wx2+2.828427124746191*V_0_sq[0]*wx1_sq+0.2357022603955158*V_2_sq[0]*dv3_sq+0.2357022603955158*V_1_sq[0]*dv2_sq+0.2357022603955158*V_0_sq[0]*dv1_sq; 
  temp_sq[1] = 2.828427124746191*V_2_sq[1]*wx3_sq+4.0*V_1[0]*V_2[1]*wx2*wx3+4.0*V_2[0]*V_1[1]*wx2*wx3+4.0*V_0[0]*V_2[1]*wx1*wx3+4.0*V_2[0]*V_0[1]*wx1*wx3+2.828427124746191*V_1_sq[1]*wx2_sq+4.0*V_0[0]*V_1[1]*wx1*wx2+4.0*V_1[0]*V_0[1]*wx1*wx2+2.828427124746191*V_0_sq[1]*wx1_sq+0.2357022603955158*V_2_sq[1]*dv3_sq+0.2357022603955158*V_1_sq[1]*dv2_sq+0.2357022603955158*V_0_sq[1]*dv1_sq; 
  temp_sq[2] = 1.154700538379252*V_0[1]*V_2[1]*dv1*wx3+1.154700538379252*V_0[0]*V_2[0]*dv1*wx3+1.154700538379252*V_0[1]*V_1[1]*dv1*wx2+1.154700538379252*V_0[0]*V_1[0]*dv1*wx2+1.632993161855453*V_0_sq[0]*dv1*wx1; 
  temp_sq[3] = 1.154700538379252*V_1[1]*V_2[1]*dv2*wx3+1.154700538379252*V_1[0]*V_2[0]*dv2*wx3+1.632993161855453*V_1_sq[0]*dv2*wx2+1.154700538379252*V_0[1]*V_1[1]*dv2*wx1+1.154700538379252*V_0[0]*V_1[0]*dv2*wx1; 
  temp_sq[4] = 1.632993161855453*V_2_sq[0]*dv3*wx3+1.154700538379252*V_1[1]*V_2[1]*dv3*wx2+1.154700538379252*V_1[0]*V_2[0]*dv3*wx2+1.154700538379252*V_0[1]*V_2[1]*dv3*wx1+1.154700538379252*V_0[0]*V_2[0]*dv3*wx1; 
  temp_sq[5] = 1.154700538379252*V_0[0]*V_2[1]*dv1*wx3+1.154700538379252*V_2[0]*V_0[1]*dv1*wx3+1.154700538379252*V_0[0]*V_1[1]*dv1*wx2+1.154700538379252*V_1[0]*V_0[1]*dv1*wx2+1.632993161855453*V_0_sq[1]*dv1*wx1; 
  temp_sq[6] = 1.154700538379252*V_1[0]*V_2[1]*dv2*wx3+1.154700538379252*V_2[0]*V_1[1]*dv2*wx3+1.632993161855453*V_1_sq[1]*dv2*wx2+1.154700538379252*V_0[0]*V_1[1]*dv2*wx1+1.154700538379252*V_1[0]*V_0[1]*dv2*wx1; 
  temp_sq[7] = 0.3333333333333333*V_0[1]*V_1[1]*dv1*dv2+0.3333333333333333*V_0[0]*V_1[0]*dv1*dv2; 
  temp_sq[8] = 1.632993161855453*V_2_sq[1]*dv3*wx3+1.154700538379252*V_1[0]*V_2[1]*dv3*wx2+1.154700538379252*V_2[0]*V_1[1]*dv3*wx2+1.154700538379252*V_0[0]*V_2[1]*dv3*wx1+1.154700538379252*V_2[0]*V_0[1]*dv3*wx1; 
  temp_sq[9] = 0.3333333333333333*V_0[1]*V_2[1]*dv1*dv3+0.3333333333333333*V_0[0]*V_2[0]*dv1*dv3; 
  temp_sq[10] = 0.3333333333333333*V_1[1]*V_2[1]*dv2*dv3+0.3333333333333333*V_1[0]*V_2[0]*dv2*dv3; 
  temp_sq[11] = 0.3333333333333333*V_0[0]*V_1[1]*dv1*dv2+0.3333333333333333*V_1[0]*V_0[1]*dv1*dv2; 
  temp_sq[12] = 0.3333333333333333*V_0[0]*V_2[1]*dv1*dv3+0.3333333333333333*V_2[0]*V_0[1]*dv1*dv3; 
  temp_sq[13] = 0.3333333333333333*V_1[0]*V_2[1]*dv2*dv3+0.3333333333333333*V_2[0]*V_1[1]*dv2*dv3; 
  temp_sq[16] = 0.210818510677892*V_0_sq[0]*dv1_sq; 
  temp_sq[17] = 0.2108185106778921*V_0_sq[1]*dv1_sq; 
  temp_sq[24] = 0.210818510677892*V_1_sq[0]*dv2_sq; 
  temp_sq[25] = 0.2108185106778921*V_1_sq[1]*dv2_sq; 
  temp_sq[32] = 0.210818510677892*V_2_sq[0]*dv3_sq; 
  temp_sq[33] = 0.2108185106778921*V_2_sq[1]*dv3_sq; 

  p_fac[0] = 0.3535533905932737*gamma_inv[9]*temp_sq[32]+0.3535533905932737*gamma_inv[8]*temp_sq[24]+0.3535533905932737*gamma_inv[7]*temp_sq[16]+0.3535533905932737*gamma_inv[6]*temp_sq[10]+0.3535533905932737*gamma_inv[5]*temp_sq[9]+0.3535533905932737*gamma_inv[4]*temp_sq[7]+0.3535533905932737*gamma_inv[3]*temp_sq[4]+0.3535533905932737*gamma_inv[2]*temp_sq[3]+0.3535533905932737*gamma_inv[1]*temp_sq[2]-1.414213562373095*GammaV[1]*temp[1]+GammaV_sq[0]*gamma[0]+0.3535533905932737*gamma_inv[0]*temp_sq[0]-1.414213562373095*GammaV[0]*temp[0]-1.414213562373095*gamma_inv[0]; 
  p_fac[1] = 0.3535533905932737*gamma_inv[9]*temp_sq[33]+0.3535533905932737*gamma_inv[8]*temp_sq[25]+0.3535533905932737*gamma_inv[7]*temp_sq[17]+0.3535533905932737*gamma_inv[6]*temp_sq[13]+0.3535533905932737*gamma_inv[5]*temp_sq[12]+0.3535533905932737*gamma_inv[4]*temp_sq[11]+0.3535533905932737*gamma_inv[3]*temp_sq[8]+0.3535533905932737*gamma_inv[2]*temp_sq[6]+0.3535533905932737*gamma_inv[1]*temp_sq[5]+0.3535533905932737*gamma_inv[0]*temp_sq[1]-1.414213562373095*GammaV[0]*temp[1]+gamma[0]*GammaV_sq[1]-1.414213562373095*temp[0]*GammaV[1]; 
  p_fac[2] = 0.3535533905932737*gamma_inv[15]*temp_sq[32]+0.3535533905932737*gamma_inv[12]*temp_sq[24]+0.3162277660168379*gamma_inv[1]*temp_sq[16]+0.3162277660168379*temp_sq[9]*gamma_inv[13]+0.3162277660168379*temp_sq[7]*gamma_inv[11]+0.3535533905932737*gamma_inv[10]*temp_sq[10]+0.3535533905932737*gamma_inv[3]*temp_sq[9]+0.3535533905932737*gamma_inv[2]*temp_sq[7]+0.3162277660168379*temp_sq[2]*gamma_inv[7]-1.414213562373095*GammaV[1]*temp[5]+0.3535533905932737*temp_sq[4]*gamma_inv[5]+0.3535533905932737*temp_sq[3]*gamma_inv[4]+0.3535533905932737*gamma_inv[0]*temp_sq[2]-1.414213562373095*GammaV[0]*temp[2]+GammaV_sq[0]*gamma[1]+0.3535533905932737*temp_sq[0]*gamma_inv[1]-1.414213562373095*gamma_inv[1]; 
  p_fac[3] = 0.3535533905932737*gamma_inv[16]*temp_sq[32]+0.3162277660168379*gamma_inv[2]*temp_sq[24]+0.3535533905932737*gamma_inv[11]*temp_sq[16]+0.3162277660168379*temp_sq[10]*gamma_inv[14]+0.3162277660168379*temp_sq[7]*gamma_inv[12]+0.3535533905932737*gamma_inv[3]*temp_sq[10]+0.3535533905932737*temp_sq[9]*gamma_inv[10]+0.3162277660168379*temp_sq[3]*gamma_inv[8]+0.3535533905932737*gamma_inv[1]*temp_sq[7]-1.414213562373095*GammaV[1]*temp[6]+0.3535533905932737*temp_sq[4]*gamma_inv[6]+0.3535533905932737*temp_sq[2]*gamma_inv[4]+0.3535533905932737*gamma_inv[0]*temp_sq[3]-1.414213562373095*GammaV[0]*temp[3]+GammaV_sq[0]*gamma[2]+0.3535533905932737*temp_sq[0]*gamma_inv[2]-1.414213562373095*gamma_inv[2]; 
  p_fac[4] = 0.3162277660168379*gamma_inv[3]*temp_sq[32]+0.3535533905932737*gamma_inv[14]*temp_sq[24]+0.3535533905932737*gamma_inv[13]*temp_sq[16]+0.3162277660168379*temp_sq[10]*gamma_inv[16]+0.3162277660168379*temp_sq[9]*gamma_inv[15]+0.3535533905932737*gamma_inv[2]*temp_sq[10]+0.3535533905932737*temp_sq[7]*gamma_inv[10]+0.3535533905932737*gamma_inv[1]*temp_sq[9]+0.3162277660168379*temp_sq[4]*gamma_inv[9]-1.414213562373095*GammaV[1]*temp[8]+0.3535533905932737*temp_sq[3]*gamma_inv[6]+0.3535533905932737*temp_sq[2]*gamma_inv[5]+0.3535533905932737*gamma_inv[0]*temp_sq[4]-1.414213562373095*GammaV[0]*temp[4]+GammaV_sq[0]*gamma[3]+0.3535533905932737*temp_sq[0]*gamma_inv[3]-1.414213562373095*gamma_inv[3]; 
  p_fac[5] = 0.3535533905932737*gamma_inv[15]*temp_sq[33]+0.3535533905932737*gamma_inv[12]*temp_sq[25]+0.3162277660168379*gamma_inv[1]*temp_sq[17]+0.3535533905932737*gamma_inv[10]*temp_sq[13]+0.3162277660168379*temp_sq[12]*gamma_inv[13]+0.3535533905932737*gamma_inv[3]*temp_sq[12]+0.3162277660168379*gamma_inv[11]*temp_sq[11]+0.3535533905932737*gamma_inv[2]*temp_sq[11]+0.3535533905932737*gamma_inv[5]*temp_sq[8]+0.3162277660168379*temp_sq[5]*gamma_inv[7]+0.3535533905932737*gamma_inv[4]*temp_sq[6]+0.3535533905932737*gamma_inv[0]*temp_sq[5]-1.414213562373095*GammaV[0]*temp[5]-1.414213562373095*GammaV[1]*temp[2]+GammaV_sq[1]*gamma[1]+0.3535533905932737*gamma_inv[1]*temp_sq[1]; 
  p_fac[6] = 0.3535533905932737*gamma_inv[16]*temp_sq[33]+0.3162277660168379*gamma_inv[2]*temp_sq[25]+0.3535533905932737*gamma_inv[11]*temp_sq[17]+0.3162277660168379*temp_sq[13]*gamma_inv[14]+0.3535533905932737*gamma_inv[3]*temp_sq[13]+0.3535533905932737*gamma_inv[10]*temp_sq[12]+0.3162277660168379*temp_sq[11]*gamma_inv[12]+0.3535533905932737*gamma_inv[1]*temp_sq[11]+0.3535533905932737*gamma_inv[6]*temp_sq[8]+0.3162277660168379*temp_sq[6]*gamma_inv[8]+0.3535533905932737*gamma_inv[0]*temp_sq[6]-1.414213562373095*GammaV[0]*temp[6]+0.3535533905932737*gamma_inv[4]*temp_sq[5]-1.414213562373095*GammaV[1]*temp[3]+GammaV_sq[1]*gamma[2]+0.3535533905932737*temp_sq[1]*gamma_inv[2]; 
  p_fac[7] = 0.3535533905932737*gamma_inv[19]*temp_sq[32]+0.3162277660168379*gamma_inv[4]*temp_sq[24]+0.3162277660168379*temp_sq[10]*gamma_inv[18]+0.3162277660168379*temp_sq[9]*gamma_inv[17]+0.3162277660168379*gamma_inv[4]*temp_sq[16]+0.3162277660168379*temp_sq[3]*gamma_inv[12]+0.3162277660168379*temp_sq[2]*gamma_inv[11]+0.3535533905932737*gamma_inv[5]*temp_sq[10]+0.3535533905932737*temp_sq[4]*gamma_inv[10]+0.3535533905932737*gamma_inv[6]*temp_sq[9]+0.3162277660168379*temp_sq[7]*gamma_inv[8]+0.3162277660168379*gamma_inv[7]*temp_sq[7]+0.3535533905932737*gamma_inv[0]*temp_sq[7]+GammaV_sq[0]*gamma[4]+0.3535533905932737*temp_sq[0]*gamma_inv[4]-1.414213562373095*gamma_inv[4]+0.3535533905932737*gamma_inv[1]*temp_sq[3]+0.3535533905932737*gamma_inv[2]*temp_sq[2]; 
  p_fac[8] = 0.3162277660168379*gamma_inv[3]*temp_sq[33]+0.3535533905932737*gamma_inv[14]*temp_sq[25]+0.3535533905932737*gamma_inv[13]*temp_sq[17]+0.3162277660168379*temp_sq[13]*gamma_inv[16]+0.3162277660168379*temp_sq[12]*gamma_inv[15]+0.3535533905932737*gamma_inv[2]*temp_sq[13]+0.3535533905932737*gamma_inv[1]*temp_sq[12]+0.3535533905932737*gamma_inv[10]*temp_sq[11]+0.3162277660168379*temp_sq[8]*gamma_inv[9]+0.3535533905932737*gamma_inv[0]*temp_sq[8]-1.414213562373095*GammaV[0]*temp[8]+0.3535533905932737*gamma_inv[6]*temp_sq[6]+0.3535533905932737*gamma_inv[5]*temp_sq[5]-1.414213562373095*GammaV[1]*temp[4]+GammaV_sq[1]*gamma[3]+0.3535533905932737*temp_sq[1]*gamma_inv[3]; 
  p_fac[9] = 0.3162277660168379*gamma_inv[5]*temp_sq[32]+0.3535533905932737*gamma_inv[18]*temp_sq[24]+0.3162277660168379*temp_sq[10]*gamma_inv[19]+0.3162277660168379*temp_sq[7]*gamma_inv[17]+0.3162277660168379*gamma_inv[5]*temp_sq[16]+0.3162277660168379*temp_sq[4]*gamma_inv[15]+0.3162277660168379*temp_sq[2]*gamma_inv[13]+0.3535533905932737*gamma_inv[4]*temp_sq[10]+0.3535533905932737*temp_sq[3]*gamma_inv[10]+0.3162277660168379*gamma_inv[9]*temp_sq[9]+0.3162277660168379*gamma_inv[7]*temp_sq[9]+0.3535533905932737*gamma_inv[0]*temp_sq[9]+0.3535533905932737*gamma_inv[6]*temp_sq[7]+GammaV_sq[0]*gamma[5]+0.3535533905932737*temp_sq[0]*gamma_inv[5]-1.414213562373095*gamma_inv[5]+0.3535533905932737*gamma_inv[1]*temp_sq[4]+0.3535533905932737*temp_sq[2]*gamma_inv[3]; 
  p_fac[10] = 0.3162277660168379*gamma_inv[6]*temp_sq[32]+0.3162277660168379*gamma_inv[6]*temp_sq[24]+0.3162277660168379*temp_sq[9]*gamma_inv[19]+0.3162277660168379*temp_sq[7]*gamma_inv[18]+0.3535533905932737*temp_sq[16]*gamma_inv[17]+0.3162277660168379*temp_sq[4]*gamma_inv[16]+0.3162277660168379*temp_sq[3]*gamma_inv[14]+0.3162277660168379*gamma_inv[9]*temp_sq[10]+0.3162277660168379*gamma_inv[8]*temp_sq[10]+0.3535533905932737*gamma_inv[0]*temp_sq[10]+0.3535533905932737*temp_sq[2]*gamma_inv[10]+0.3535533905932737*gamma_inv[4]*temp_sq[9]+0.3535533905932737*gamma_inv[5]*temp_sq[7]+GammaV_sq[0]*gamma[6]+0.3535533905932737*temp_sq[0]*gamma_inv[6]-1.414213562373095*gamma_inv[6]+0.3535533905932737*gamma_inv[2]*temp_sq[4]+0.3535533905932737*gamma_inv[3]*temp_sq[3]; 
  p_fac[11] = 0.3535533905932737*gamma_inv[19]*temp_sq[33]+0.3162277660168379*gamma_inv[4]*temp_sq[25]+0.3162277660168379*temp_sq[13]*gamma_inv[18]+0.3162277660168379*gamma_inv[4]*temp_sq[17]+0.3162277660168379*temp_sq[12]*gamma_inv[17]+0.3535533905932737*gamma_inv[5]*temp_sq[13]+0.3535533905932737*gamma_inv[6]*temp_sq[12]+0.3162277660168379*temp_sq[6]*gamma_inv[12]+0.3162277660168379*gamma_inv[8]*temp_sq[11]+0.3162277660168379*gamma_inv[7]*temp_sq[11]+0.3535533905932737*gamma_inv[0]*temp_sq[11]+0.3162277660168379*temp_sq[5]*gamma_inv[11]+0.3535533905932737*temp_sq[8]*gamma_inv[10]+0.3535533905932737*gamma_inv[1]*temp_sq[6]+0.3535533905932737*gamma_inv[2]*temp_sq[5]+GammaV_sq[1]*gamma[4]+0.3535533905932737*temp_sq[1]*gamma_inv[4]; 
  p_fac[12] = 0.3162277660168379*gamma_inv[5]*temp_sq[33]+0.3535533905932737*gamma_inv[18]*temp_sq[25]+0.3162277660168379*temp_sq[13]*gamma_inv[19]+0.3162277660168379*gamma_inv[5]*temp_sq[17]+0.3162277660168379*temp_sq[11]*gamma_inv[17]+0.3162277660168379*temp_sq[8]*gamma_inv[15]+0.3535533905932737*gamma_inv[4]*temp_sq[13]+0.3162277660168379*temp_sq[5]*gamma_inv[13]+0.3162277660168379*gamma_inv[9]*temp_sq[12]+0.3162277660168379*gamma_inv[7]*temp_sq[12]+0.3535533905932737*gamma_inv[0]*temp_sq[12]+0.3535533905932737*gamma_inv[6]*temp_sq[11]+0.3535533905932737*temp_sq[6]*gamma_inv[10]+0.3535533905932737*gamma_inv[1]*temp_sq[8]+GammaV_sq[1]*gamma[5]+0.3535533905932737*gamma_inv[3]*temp_sq[5]+0.3535533905932737*temp_sq[1]*gamma_inv[5]; 
  p_fac[13] = 0.3162277660168379*gamma_inv[6]*temp_sq[33]+0.3162277660168379*gamma_inv[6]*temp_sq[25]+0.3162277660168379*temp_sq[12]*gamma_inv[19]+0.3162277660168379*temp_sq[11]*gamma_inv[18]+0.3535533905932737*gamma_inv[17]*temp_sq[17]+0.3162277660168379*temp_sq[8]*gamma_inv[16]+0.3162277660168379*temp_sq[6]*gamma_inv[14]+0.3162277660168379*gamma_inv[9]*temp_sq[13]+0.3162277660168379*gamma_inv[8]*temp_sq[13]+0.3535533905932737*gamma_inv[0]*temp_sq[13]+0.3535533905932737*gamma_inv[4]*temp_sq[12]+0.3535533905932737*gamma_inv[5]*temp_sq[11]+0.3535533905932737*temp_sq[5]*gamma_inv[10]+0.3535533905932737*gamma_inv[2]*temp_sq[8]+GammaV_sq[1]*gamma[6]+0.3535533905932737*gamma_inv[3]*temp_sq[6]+0.3535533905932737*temp_sq[1]*gamma_inv[6]; 
  p_fac[14] = 0.3162277660168379*gamma_inv[10]*temp_sq[32]+0.3162277660168379*gamma_inv[10]*temp_sq[24]+0.3162277660168379*temp_sq[4]*gamma_inv[19]+0.3162277660168379*temp_sq[3]*gamma_inv[18]+0.3162277660168379*temp_sq[2]*gamma_inv[17]+0.3162277660168379*gamma_inv[10]*temp_sq[16]+0.3162277660168379*temp_sq[9]*gamma_inv[16]+0.3162277660168379*temp_sq[10]*gamma_inv[15]+0.3162277660168379*temp_sq[7]*gamma_inv[14]+0.3162277660168379*temp_sq[7]*gamma_inv[13]+0.3162277660168379*temp_sq[10]*gamma_inv[12]+0.3162277660168379*temp_sq[9]*gamma_inv[11]+GammaV_sq[0]*gamma[10]+0.3535533905932737*gamma_inv[1]*temp_sq[10]+0.3535533905932737*temp_sq[0]*gamma_inv[10]-1.414213562373095*gamma_inv[10]+0.3535533905932737*gamma_inv[2]*temp_sq[9]+0.3535533905932737*gamma_inv[3]*temp_sq[7]+0.3535533905932737*temp_sq[2]*gamma_inv[6]+0.3535533905932737*temp_sq[3]*gamma_inv[5]+0.3535533905932737*gamma_inv[4]*temp_sq[4]; 
  p_fac[15] = 0.3162277660168379*gamma_inv[10]*temp_sq[33]+0.3162277660168379*gamma_inv[10]*temp_sq[25]+0.3162277660168379*temp_sq[8]*gamma_inv[19]+0.3162277660168379*temp_sq[6]*gamma_inv[18]+0.3162277660168379*gamma_inv[10]*temp_sq[17]+0.3162277660168379*temp_sq[5]*gamma_inv[17]+0.3162277660168379*temp_sq[12]*gamma_inv[16]+0.3162277660168379*temp_sq[13]*gamma_inv[15]+0.3162277660168379*temp_sq[11]*gamma_inv[14]+0.3162277660168379*gamma_inv[12]*temp_sq[13]+0.3535533905932737*gamma_inv[1]*temp_sq[13]+0.3162277660168379*temp_sq[11]*gamma_inv[13]+0.3162277660168379*gamma_inv[11]*temp_sq[12]+0.3535533905932737*gamma_inv[2]*temp_sq[12]+0.3535533905932737*gamma_inv[3]*temp_sq[11]+GammaV_sq[1]*gamma[10]+0.3535533905932737*temp_sq[1]*gamma_inv[10]+0.3535533905932737*gamma_inv[4]*temp_sq[8]+0.3535533905932737*gamma_inv[5]*temp_sq[6]+0.3535533905932737*temp_sq[5]*gamma_inv[6]; 
  p_fac[16] = 0.3535533905932737*temp_sq[10]*gamma_inv[17]+0.2258769757263128*gamma_inv[7]*temp_sq[16]+0.3535533905932737*gamma_inv[0]*temp_sq[16]+0.3535533905932737*temp_sq[4]*gamma_inv[13]+0.3535533905932737*temp_sq[3]*gamma_inv[11]+0.3162277660168379*gamma_inv[5]*temp_sq[9]+GammaV_sq[0]*gamma[7]+0.3162277660168379*gamma_inv[4]*temp_sq[7]+0.3535533905932737*temp_sq[0]*gamma_inv[7]-1.414213562373095*gamma_inv[7]+0.3162277660168379*gamma_inv[1]*temp_sq[2]; 
  p_fac[17] = 0.2258769757263128*gamma_inv[7]*temp_sq[17]+0.3535533905932737*gamma_inv[0]*temp_sq[17]+0.3535533905932737*temp_sq[13]*gamma_inv[17]+0.3535533905932737*temp_sq[8]*gamma_inv[13]+0.3162277660168379*gamma_inv[5]*temp_sq[12]+0.3162277660168379*gamma_inv[4]*temp_sq[11]+0.3535533905932737*temp_sq[6]*gamma_inv[11]+1.0*GammaV_sq[1]*gamma[7]+0.3535533905932737*temp_sq[1]*gamma_inv[7]+0.3162277660168379*gamma_inv[1]*temp_sq[5]; 
  p_fac[18] = 0.3162277660168379*gamma_inv[11]*temp_sq[24]+0.3535533905932737*temp_sq[4]*gamma_inv[17]+0.2258769757263128*gamma_inv[11]*temp_sq[16]+0.3535533905932737*gamma_inv[2]*temp_sq[16]+0.3535533905932737*temp_sq[10]*gamma_inv[13]+0.2828427124746191*temp_sq[7]*gamma_inv[12]+GammaV_sq[0]*gamma[11]+0.3535533905932737*temp_sq[0]*gamma_inv[11]-1.414213562373095*gamma_inv[11]+0.3162277660168379*temp_sq[9]*gamma_inv[10]+0.3162277660168379*gamma_inv[1]*temp_sq[7]+0.3535533905932737*temp_sq[3]*gamma_inv[7]+0.3162277660168379*temp_sq[2]*gamma_inv[4]; 
  p_fac[19] = 0.3162277660168379*gamma_inv[13]*temp_sq[32]+0.3535533905932737*temp_sq[3]*gamma_inv[17]+0.2258769757263128*gamma_inv[13]*temp_sq[16]+0.3535533905932737*gamma_inv[3]*temp_sq[16]+0.2828427124746191*temp_sq[9]*gamma_inv[15]+GammaV_sq[0]*gamma[13]+0.3535533905932737*temp_sq[0]*gamma_inv[13]-1.414213562373095*gamma_inv[13]+0.3535533905932737*temp_sq[10]*gamma_inv[11]+0.3162277660168379*temp_sq[7]*gamma_inv[10]+0.3162277660168379*gamma_inv[1]*temp_sq[9]+0.3535533905932737*temp_sq[4]*gamma_inv[7]+0.3162277660168379*temp_sq[2]*gamma_inv[5]; 
  p_fac[20] = 0.3162277660168379*gamma_inv[11]*temp_sq[25]+0.2258769757263128*gamma_inv[11]*temp_sq[17]+0.3535533905932737*gamma_inv[2]*temp_sq[17]+0.3535533905932737*temp_sq[8]*gamma_inv[17]+0.3535533905932737*gamma_inv[13]*temp_sq[13]+0.3162277660168379*gamma_inv[10]*temp_sq[12]+0.282842712474619*temp_sq[11]*gamma_inv[12]+1.0*GammaV_sq[1]*gamma[11]+0.3162277660168379*gamma_inv[1]*temp_sq[11]+0.3535533905932737*temp_sq[1]*gamma_inv[11]+0.3535533905932737*temp_sq[6]*gamma_inv[7]+0.3162277660168379*gamma_inv[4]*temp_sq[5]; 
  p_fac[21] = 0.3162277660168379*gamma_inv[13]*temp_sq[33]+0.2258769757263128*gamma_inv[13]*temp_sq[17]+0.3535533905932737*gamma_inv[3]*temp_sq[17]+0.3535533905932737*temp_sq[6]*gamma_inv[17]+0.282842712474619*temp_sq[12]*gamma_inv[15]+1.0*GammaV_sq[1]*gamma[13]+0.3535533905932737*gamma_inv[11]*temp_sq[13]+0.3535533905932737*temp_sq[1]*gamma_inv[13]+0.3162277660168379*gamma_inv[1]*temp_sq[12]+0.3162277660168379*gamma_inv[10]*temp_sq[11]+0.3535533905932737*gamma_inv[7]*temp_sq[8]+0.3162277660168379*gamma_inv[5]*temp_sq[5]; 
  p_fac[22] = 0.3162277660168379*gamma_inv[17]*temp_sq[32]+0.3162277660168379*gamma_inv[17]*temp_sq[24]+0.2828427124746191*temp_sq[9]*gamma_inv[19]+0.2828427124746191*temp_sq[7]*gamma_inv[18]+GammaV_sq[0]*gamma[17]+0.2258769757263128*temp_sq[16]*gamma_inv[17]+0.3535533905932737*temp_sq[0]*gamma_inv[17]-1.414213562373095*gamma_inv[17]+0.3535533905932737*gamma_inv[6]*temp_sq[16]+0.3535533905932737*temp_sq[3]*gamma_inv[13]+0.3535533905932737*temp_sq[4]*gamma_inv[11]+0.3535533905932737*gamma_inv[7]*temp_sq[10]+0.3162277660168379*temp_sq[2]*gamma_inv[10]+0.3162277660168379*gamma_inv[4]*temp_sq[9]+0.3162277660168379*gamma_inv[5]*temp_sq[7]; 
  p_fac[23] = 0.3162277660168379*gamma_inv[17]*temp_sq[33]+0.3162277660168379*gamma_inv[17]*temp_sq[25]+0.282842712474619*temp_sq[12]*gamma_inv[19]+0.282842712474619*temp_sq[11]*gamma_inv[18]+1.0*GammaV_sq[1]*gamma[17]+0.2258769757263128*gamma_inv[17]*temp_sq[17]+0.3535533905932737*gamma_inv[6]*temp_sq[17]+0.3535533905932737*temp_sq[1]*gamma_inv[17]+0.3535533905932737*gamma_inv[7]*temp_sq[13]+0.3535533905932737*temp_sq[6]*gamma_inv[13]+0.3162277660168379*gamma_inv[4]*temp_sq[12]+0.3162277660168379*gamma_inv[5]*temp_sq[11]+0.3535533905932737*temp_sq[8]*gamma_inv[11]+0.3162277660168379*temp_sq[5]*gamma_inv[10]; 
  p_fac[24] = 0.2258769757263128*gamma_inv[8]*temp_sq[24]+0.3535533905932737*gamma_inv[0]*temp_sq[24]+0.3535533905932737*temp_sq[9]*gamma_inv[18]+0.3535533905932737*temp_sq[4]*gamma_inv[14]+0.3535533905932737*temp_sq[2]*gamma_inv[12]+0.3162277660168379*gamma_inv[6]*temp_sq[10]+GammaV_sq[0]*gamma[8]+0.3535533905932737*temp_sq[0]*gamma_inv[8]-1.414213562373095*gamma_inv[8]+0.3162277660168379*gamma_inv[4]*temp_sq[7]+0.3162277660168379*gamma_inv[2]*temp_sq[3]; 
  p_fac[25] = 0.2258769757263128*gamma_inv[8]*temp_sq[25]+0.3535533905932737*gamma_inv[0]*temp_sq[25]+0.3535533905932737*temp_sq[12]*gamma_inv[18]+0.3535533905932737*temp_sq[8]*gamma_inv[14]+0.3162277660168379*gamma_inv[6]*temp_sq[13]+0.3535533905932737*temp_sq[5]*gamma_inv[12]+0.3162277660168379*gamma_inv[4]*temp_sq[11]+1.0*GammaV_sq[1]*gamma[8]+0.3535533905932737*temp_sq[1]*gamma_inv[8]+0.3162277660168379*gamma_inv[2]*temp_sq[6]; 
  p_fac[26] = 0.2258769757263128*gamma_inv[12]*temp_sq[24]+0.3535533905932737*gamma_inv[1]*temp_sq[24]+0.3535533905932737*temp_sq[4]*gamma_inv[18]+0.3162277660168379*gamma_inv[12]*temp_sq[16]+0.3535533905932737*temp_sq[9]*gamma_inv[14]+GammaV_sq[0]*gamma[12]+0.3535533905932737*temp_sq[0]*gamma_inv[12]-1.414213562373095*gamma_inv[12]+0.2828427124746191*temp_sq[7]*gamma_inv[11]+0.3162277660168379*gamma_inv[10]*temp_sq[10]+0.3535533905932737*temp_sq[2]*gamma_inv[8]+0.3162277660168379*gamma_inv[2]*temp_sq[7]+0.3162277660168379*temp_sq[3]*gamma_inv[4]; 
  p_fac[27] = 0.3162277660168379*gamma_inv[14]*temp_sq[32]+0.2258769757263128*gamma_inv[14]*temp_sq[24]+0.3535533905932737*gamma_inv[3]*temp_sq[24]+0.3535533905932737*temp_sq[2]*gamma_inv[18]+0.2828427124746191*temp_sq[10]*gamma_inv[16]+GammaV_sq[0]*gamma[14]+0.3535533905932737*temp_sq[0]*gamma_inv[14]-1.414213562373095*gamma_inv[14]+0.3535533905932737*temp_sq[9]*gamma_inv[12]+0.3162277660168379*gamma_inv[2]*temp_sq[10]+0.3162277660168379*temp_sq[7]*gamma_inv[10]+0.3535533905932737*temp_sq[4]*gamma_inv[8]+0.3162277660168379*temp_sq[3]*gamma_inv[6]; 
  p_fac[28] = 0.2258769757263128*gamma_inv[12]*temp_sq[25]+0.3535533905932737*gamma_inv[1]*temp_sq[25]+0.3535533905932737*temp_sq[8]*gamma_inv[18]+0.3162277660168379*gamma_inv[12]*temp_sq[17]+0.3535533905932737*temp_sq[12]*gamma_inv[14]+0.3162277660168379*gamma_inv[10]*temp_sq[13]+1.0*GammaV_sq[1]*gamma[12]+0.3535533905932737*temp_sq[1]*gamma_inv[12]+0.282842712474619*gamma_inv[11]*temp_sq[11]+0.3162277660168379*gamma_inv[2]*temp_sq[11]+0.3535533905932737*temp_sq[5]*gamma_inv[8]+0.3162277660168379*gamma_inv[4]*temp_sq[6]; 
  p_fac[29] = 0.3162277660168379*gamma_inv[14]*temp_sq[33]+0.2258769757263128*gamma_inv[14]*temp_sq[25]+0.3535533905932737*gamma_inv[3]*temp_sq[25]+0.3535533905932737*temp_sq[5]*gamma_inv[18]+0.282842712474619*temp_sq[13]*gamma_inv[16]+1.0*GammaV_sq[1]*gamma[14]+0.3535533905932737*temp_sq[1]*gamma_inv[14]+0.3162277660168379*gamma_inv[2]*temp_sq[13]+0.3535533905932737*gamma_inv[12]*temp_sq[12]+0.3162277660168379*gamma_inv[10]*temp_sq[11]+0.3535533905932737*gamma_inv[8]*temp_sq[8]+0.3162277660168379*gamma_inv[6]*temp_sq[6]; 
  p_fac[30] = 0.3162277660168379*gamma_inv[18]*temp_sq[32]+0.2258769757263128*gamma_inv[18]*temp_sq[24]+0.3535533905932737*gamma_inv[5]*temp_sq[24]+0.2828427124746191*temp_sq[10]*gamma_inv[19]+GammaV_sq[0]*gamma[18]+0.3162277660168379*temp_sq[16]*gamma_inv[18]+0.3535533905932737*temp_sq[0]*gamma_inv[18]-1.414213562373095*gamma_inv[18]+0.2828427124746191*temp_sq[7]*gamma_inv[17]+0.3535533905932737*temp_sq[2]*gamma_inv[14]+0.3535533905932737*temp_sq[4]*gamma_inv[12]+0.3162277660168379*gamma_inv[4]*temp_sq[10]+0.3162277660168379*temp_sq[3]*gamma_inv[10]+0.3535533905932737*gamma_inv[8]*temp_sq[9]+0.3162277660168379*gamma_inv[6]*temp_sq[7]; 
  p_fac[31] = 0.3162277660168379*gamma_inv[18]*temp_sq[33]+0.2258769757263128*gamma_inv[18]*temp_sq[25]+0.3535533905932737*gamma_inv[5]*temp_sq[25]+0.282842712474619*temp_sq[13]*gamma_inv[19]+1.0*GammaV_sq[1]*gamma[18]+0.3162277660168379*temp_sq[17]*gamma_inv[18]+0.3535533905932737*temp_sq[1]*gamma_inv[18]+0.282842712474619*temp_sq[11]*gamma_inv[17]+0.3535533905932737*temp_sq[5]*gamma_inv[14]+0.3162277660168379*gamma_inv[4]*temp_sq[13]+0.3535533905932737*gamma_inv[8]*temp_sq[12]+0.3535533905932737*temp_sq[8]*gamma_inv[12]+0.3162277660168379*gamma_inv[6]*temp_sq[11]+0.3162277660168379*temp_sq[6]*gamma_inv[10]; 
  p_fac[32] = 0.2258769757263128*gamma_inv[9]*temp_sq[32]+0.3535533905932737*gamma_inv[0]*temp_sq[32]+0.3535533905932737*temp_sq[7]*gamma_inv[19]+0.3535533905932737*temp_sq[3]*gamma_inv[16]+0.3535533905932737*temp_sq[2]*gamma_inv[15]+0.3162277660168379*gamma_inv[6]*temp_sq[10]+GammaV_sq[0]*gamma[9]+0.3162277660168379*gamma_inv[5]*temp_sq[9]+0.3535533905932737*temp_sq[0]*gamma_inv[9]-1.414213562373095*gamma_inv[9]+0.3162277660168379*gamma_inv[3]*temp_sq[4]; 
  p_fac[33] = 0.2258769757263128*gamma_inv[9]*temp_sq[33]+0.3535533905932737*gamma_inv[0]*temp_sq[33]+0.3535533905932737*temp_sq[11]*gamma_inv[19]+0.3535533905932737*temp_sq[6]*gamma_inv[16]+0.3535533905932737*temp_sq[5]*gamma_inv[15]+0.3162277660168379*gamma_inv[6]*temp_sq[13]+0.3162277660168379*gamma_inv[5]*temp_sq[12]+1.0*GammaV_sq[1]*gamma[9]+0.3535533905932737*temp_sq[1]*gamma_inv[9]+0.3162277660168379*gamma_inv[3]*temp_sq[8]; 
  p_fac[34] = 0.2258769757263128*gamma_inv[15]*temp_sq[32]+0.3535533905932737*gamma_inv[1]*temp_sq[32]+0.3535533905932737*temp_sq[3]*gamma_inv[19]+0.3162277660168379*gamma_inv[15]*temp_sq[16]+0.3535533905932737*temp_sq[7]*gamma_inv[16]+GammaV_sq[0]*gamma[15]+0.3535533905932737*temp_sq[0]*gamma_inv[15]-1.414213562373095*gamma_inv[15]+0.2828427124746191*temp_sq[9]*gamma_inv[13]+0.3162277660168379*gamma_inv[10]*temp_sq[10]+0.3162277660168379*gamma_inv[3]*temp_sq[9]+0.3535533905932737*temp_sq[2]*gamma_inv[9]+0.3162277660168379*temp_sq[4]*gamma_inv[5]; 
  p_fac[35] = 0.2258769757263128*gamma_inv[16]*temp_sq[32]+0.3535533905932737*gamma_inv[2]*temp_sq[32]+0.3162277660168379*gamma_inv[16]*temp_sq[24]+0.3535533905932737*temp_sq[2]*gamma_inv[19]+GammaV_sq[0]*gamma[16]+0.3535533905932737*temp_sq[0]*gamma_inv[16]-1.414213562373095*gamma_inv[16]+0.3535533905932737*temp_sq[7]*gamma_inv[15]+0.2828427124746191*temp_sq[10]*gamma_inv[14]+0.3162277660168379*gamma_inv[3]*temp_sq[10]+0.3162277660168379*temp_sq[9]*gamma_inv[10]+0.3535533905932737*temp_sq[3]*gamma_inv[9]+0.3162277660168379*temp_sq[4]*gamma_inv[6]; 
  p_fac[36] = 0.2258769757263128*gamma_inv[15]*temp_sq[33]+0.3535533905932737*gamma_inv[1]*temp_sq[33]+0.3535533905932737*temp_sq[6]*gamma_inv[19]+0.3162277660168379*gamma_inv[15]*temp_sq[17]+0.3535533905932737*temp_sq[11]*gamma_inv[16]+1.0*GammaV_sq[1]*gamma[15]+0.3535533905932737*temp_sq[1]*gamma_inv[15]+0.3162277660168379*gamma_inv[10]*temp_sq[13]+0.282842712474619*temp_sq[12]*gamma_inv[13]+0.3162277660168379*gamma_inv[3]*temp_sq[12]+0.3535533905932737*temp_sq[5]*gamma_inv[9]+0.3162277660168379*gamma_inv[5]*temp_sq[8]; 
  p_fac[37] = 0.2258769757263128*gamma_inv[16]*temp_sq[33]+0.3535533905932737*gamma_inv[2]*temp_sq[33]+0.3162277660168379*gamma_inv[16]*temp_sq[25]+0.3535533905932737*temp_sq[5]*gamma_inv[19]+1.0*GammaV_sq[1]*gamma[16]+0.3535533905932737*temp_sq[1]*gamma_inv[16]+0.3535533905932737*temp_sq[11]*gamma_inv[15]+0.282842712474619*temp_sq[13]*gamma_inv[14]+0.3162277660168379*gamma_inv[3]*temp_sq[13]+0.3162277660168379*gamma_inv[10]*temp_sq[12]+0.3535533905932737*temp_sq[6]*gamma_inv[9]+0.3162277660168379*gamma_inv[6]*temp_sq[8]; 
  p_fac[38] = 0.2258769757263128*gamma_inv[19]*temp_sq[32]+0.3535533905932737*gamma_inv[4]*temp_sq[32]+0.3162277660168379*gamma_inv[19]*temp_sq[24]+GammaV_sq[0]*gamma[19]+0.3162277660168379*temp_sq[16]*gamma_inv[19]+0.3535533905932737*temp_sq[0]*gamma_inv[19]-1.414213562373095*gamma_inv[19]+0.2828427124746191*temp_sq[10]*gamma_inv[18]+0.2828427124746191*temp_sq[9]*gamma_inv[17]+0.3535533905932737*temp_sq[2]*gamma_inv[16]+0.3535533905932737*temp_sq[3]*gamma_inv[15]+0.3162277660168379*gamma_inv[5]*temp_sq[10]+0.3162277660168379*temp_sq[4]*gamma_inv[10]+0.3162277660168379*gamma_inv[6]*temp_sq[9]+0.3535533905932737*temp_sq[7]*gamma_inv[9]; 
  p_fac[39] = 0.2258769757263128*gamma_inv[19]*temp_sq[33]+0.3535533905932737*gamma_inv[4]*temp_sq[33]+0.3162277660168379*gamma_inv[19]*temp_sq[25]+1.0*GammaV_sq[1]*gamma[19]+0.3162277660168379*temp_sq[17]*gamma_inv[19]+0.3535533905932737*temp_sq[1]*gamma_inv[19]+0.282842712474619*temp_sq[13]*gamma_inv[18]+0.282842712474619*temp_sq[12]*gamma_inv[17]+0.3535533905932737*temp_sq[5]*gamma_inv[16]+0.3535533905932737*temp_sq[6]*gamma_inv[15]+0.3162277660168379*gamma_inv[5]*temp_sq[13]+0.3162277660168379*gamma_inv[6]*temp_sq[12]+0.3535533905932737*gamma_inv[9]*temp_sq[11]+0.3162277660168379*temp_sq[8]*gamma_inv[10]; 

  sr_pressure[0] += (0.2357022603955158*f[39]*p_fac[39]+0.2357022603955158*f[38]*p_fac[38]+0.2357022603955158*f[37]*p_fac[37]+0.2357022603955158*f[36]*p_fac[36]+0.2357022603955158*f[35]*p_fac[35]+0.2357022603955158*f[34]*p_fac[34]+0.2357022603955158*f[33]*p_fac[33]+0.2357022603955158*f[32]*p_fac[32]+0.2357022603955158*f[31]*p_fac[31]+0.2357022603955158*f[30]*p_fac[30]+0.2357022603955158*f[29]*p_fac[29]+0.2357022603955158*f[28]*p_fac[28]+0.2357022603955158*f[27]*p_fac[27]+0.2357022603955158*f[26]*p_fac[26]+0.2357022603955158*f[25]*p_fac[25]+0.2357022603955158*f[24]*p_fac[24]+0.2357022603955158*f[23]*p_fac[23]+0.2357022603955158*f[22]*p_fac[22]+0.2357022603955158*f[21]*p_fac[21]+0.2357022603955158*f[20]*p_fac[20]+0.2357022603955158*f[19]*p_fac[19]+0.2357022603955158*f[18]*p_fac[18]+0.2357022603955158*f[17]*p_fac[17]+0.2357022603955158*f[16]*p_fac[16]+0.2357022603955158*f[15]*p_fac[15]+0.2357022603955158*f[14]*p_fac[14]+0.2357022603955158*f[13]*p_fac[13]+0.2357022603955158*f[12]*p_fac[12]+0.2357022603955158*f[11]*p_fac[11]+0.2357022603955158*f[10]*p_fac[10]+0.2357022603955158*f[9]*p_fac[9]+0.2357022603955158*f[8]*p_fac[8]+0.2357022603955158*f[7]*p_fac[7]+0.2357022603955158*f[6]*p_fac[6]+0.2357022603955158*f[5]*p_fac[5]+0.2357022603955158*f[4]*p_fac[4]+0.2357022603955158*f[3]*p_fac[3]+0.2357022603955158*f[2]*p_fac[2]+0.2357022603955158*f[1]*p_fac[1]+0.2357022603955158*f[0]*p_fac[0])*volFact; 
  sr_pressure[1] += (0.2357022603955158*f[38]*p_fac[39]+0.2357022603955158*p_fac[38]*f[39]+0.2357022603955158*f[35]*p_fac[37]+0.2357022603955158*p_fac[35]*f[37]+0.2357022603955158*f[34]*p_fac[36]+0.2357022603955158*p_fac[34]*f[36]+0.2357022603955158*f[32]*p_fac[33]+0.2357022603955158*p_fac[32]*f[33]+0.2357022603955158*f[30]*p_fac[31]+0.2357022603955158*p_fac[30]*f[31]+0.2357022603955158*f[27]*p_fac[29]+0.2357022603955158*p_fac[27]*f[29]+0.2357022603955158*f[26]*p_fac[28]+0.2357022603955158*p_fac[26]*f[28]+0.2357022603955158*f[24]*p_fac[25]+0.2357022603955158*p_fac[24]*f[25]+0.2357022603955158*f[22]*p_fac[23]+0.2357022603955158*p_fac[22]*f[23]+0.2357022603955158*f[19]*p_fac[21]+0.2357022603955158*p_fac[19]*f[21]+0.2357022603955158*f[18]*p_fac[20]+0.2357022603955158*p_fac[18]*f[20]+0.2357022603955158*f[16]*p_fac[17]+0.2357022603955158*p_fac[16]*f[17]+0.2357022603955158*f[14]*p_fac[15]+0.2357022603955158*p_fac[14]*f[15]+0.2357022603955158*f[10]*p_fac[13]+0.2357022603955158*p_fac[10]*f[13]+0.2357022603955158*f[9]*p_fac[12]+0.2357022603955158*p_fac[9]*f[12]+0.2357022603955158*f[7]*p_fac[11]+0.2357022603955158*p_fac[7]*f[11]+0.2357022603955158*f[4]*p_fac[8]+0.2357022603955158*p_fac[4]*f[8]+0.2357022603955158*f[3]*p_fac[6]+0.2357022603955158*p_fac[3]*f[6]+0.2357022603955158*f[2]*p_fac[5]+0.2357022603955158*p_fac[2]*f[5]+0.2357022603955158*f[0]*p_fac[1]+0.2357022603955158*p_fac[0]*f[1])*volFact; 
} 
