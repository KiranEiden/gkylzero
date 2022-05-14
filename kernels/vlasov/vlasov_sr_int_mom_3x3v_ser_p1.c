#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_sr_int_mom_3x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]*dxv[5]*0.015625; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double *p0_over_gamma = &p_over_gamma[0]; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double *p1_over_gamma = &p_over_gamma[8]; 
  const double wx3 = w[5], dv3 = dxv[5]; 
  const double *p2_over_gamma = &p_over_gamma[16]; 
 
  out[0] += 8.0*f[0]*volFact; 
  out[1] += (2.828427124746191*p0_over_gamma[7]*f[41]+2.828427124746191*p0_over_gamma[6]*f[21]+2.828427124746191*p0_over_gamma[5]*f[20]+2.828427124746191*p0_over_gamma[4]*f[16]+2.828427124746191*p0_over_gamma[3]*f[6]+2.828427124746191*p0_over_gamma[2]*f[5]+2.828427124746191*p0_over_gamma[1]*f[4]+2.828427124746191*f[0]*p0_over_gamma[0])*volFact; 
  out[2] += (2.828427124746191*p1_over_gamma[7]*f[41]+2.828427124746191*p1_over_gamma[6]*f[21]+2.828427124746191*p1_over_gamma[5]*f[20]+2.828427124746191*p1_over_gamma[4]*f[16]+2.828427124746191*p1_over_gamma[3]*f[6]+2.828427124746191*p1_over_gamma[2]*f[5]+2.828427124746191*p1_over_gamma[1]*f[4]+2.828427124746191*f[0]*p1_over_gamma[0])*volFact; 
  out[3] += (2.828427124746191*p2_over_gamma[7]*f[41]+2.828427124746191*p2_over_gamma[6]*f[21]+2.828427124746191*p2_over_gamma[5]*f[20]+2.828427124746191*p2_over_gamma[4]*f[16]+2.828427124746191*p2_over_gamma[3]*f[6]+2.828427124746191*p2_over_gamma[2]*f[5]+2.828427124746191*p2_over_gamma[1]*f[4]+2.828427124746191*f[0]*p2_over_gamma[0])*volFact; 
  out[4] += volFact*(2.828427124746191*p2_over_gamma[7]*f[41]*wx3+2.828427124746191*p2_over_gamma[6]*f[21]*wx3+2.828427124746191*p2_over_gamma[5]*f[20]*wx3+2.828427124746191*p2_over_gamma[4]*f[16]*wx3+2.828427124746191*p2_over_gamma[3]*f[6]*wx3+2.828427124746191*p2_over_gamma[2]*f[5]*wx3+2.828427124746191*p2_over_gamma[1]*f[4]*wx3+2.828427124746191*f[0]*p2_over_gamma[0]*wx3+2.828427124746191*p1_over_gamma[7]*f[41]*wx2+2.828427124746191*p1_over_gamma[6]*f[21]*wx2+2.828427124746191*p1_over_gamma[5]*f[20]*wx2+2.828427124746191*p1_over_gamma[4]*f[16]*wx2+2.828427124746191*p1_over_gamma[3]*f[6]*wx2+2.828427124746191*p1_over_gamma[2]*f[5]*wx2+2.828427124746191*p1_over_gamma[1]*f[4]*wx2+2.828427124746191*f[0]*p1_over_gamma[0]*wx2+2.828427124746191*p0_over_gamma[7]*f[41]*wx1+2.828427124746191*p0_over_gamma[6]*f[21]*wx1+2.828427124746191*p0_over_gamma[5]*f[20]*wx1+2.828427124746191*p0_over_gamma[4]*f[16]*wx1+2.828427124746191*p0_over_gamma[3]*f[6]*wx1+2.828427124746191*p0_over_gamma[2]*f[5]*wx1+2.828427124746191*p0_over_gamma[1]*f[4]*wx1+2.828427124746191*f[0]*p0_over_gamma[0]*wx1+0.8164965809277261*p2_over_gamma[4]*f[41]*dv3+0.8164965809277261*p2_over_gamma[2]*f[21]*dv3+0.8164965809277261*p2_over_gamma[1]*f[20]*dv3+0.8164965809277261*p2_over_gamma[7]*f[16]*dv3+0.8164965809277261*f[5]*p2_over_gamma[6]*dv3+0.8164965809277261*p2_over_gamma[0]*f[6]*dv3+0.8164965809277261*f[4]*p2_over_gamma[5]*dv3+0.8164965809277261*f[0]*p2_over_gamma[3]*dv3+0.8164965809277261*p1_over_gamma[5]*f[41]*dv2+0.8164965809277261*p1_over_gamma[3]*f[21]*dv2+0.8164965809277261*p1_over_gamma[7]*f[20]*dv2+0.8164965809277261*p1_over_gamma[1]*f[16]*dv2+0.8164965809277261*f[6]*p1_over_gamma[6]*dv2+0.8164965809277261*p1_over_gamma[0]*f[5]*dv2+0.8164965809277261*f[4]*p1_over_gamma[4]*dv2+0.8164965809277261*f[0]*p1_over_gamma[2]*dv2+0.8164965809277261*p0_over_gamma[6]*f[41]*dv1+0.8164965809277261*p0_over_gamma[7]*f[21]*dv1+0.8164965809277261*p0_over_gamma[3]*f[20]*dv1+0.8164965809277261*p0_over_gamma[2]*f[16]*dv1+0.8164965809277261*p0_over_gamma[5]*f[6]*dv1+0.8164965809277261*p0_over_gamma[4]*f[5]*dv1+0.8164965809277261*p0_over_gamma[0]*f[4]*dv1+0.8164965809277261*f[0]*p0_over_gamma[1]*dv1); 
} 
