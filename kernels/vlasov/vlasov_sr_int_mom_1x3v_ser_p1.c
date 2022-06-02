#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_sr_int_mom_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]*0.0625; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double *p0_over_gamma = &p_over_gamma[0]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double *p1_over_gamma = &p_over_gamma[8]; 
  const double wx3 = w[3], dv3 = dxv[3]; 
  const double *p2_over_gamma = &p_over_gamma[16]; 
 
  out[0] += 4.0*f[0]*volFact; 
  out[1] += (1.414213562373095*p0_over_gamma[7]*f[14]+1.414213562373095*p0_over_gamma[6]*f[10]+1.414213562373095*p0_over_gamma[5]*f[9]+1.414213562373095*p0_over_gamma[4]*f[7]+1.414213562373095*p0_over_gamma[3]*f[4]+1.414213562373095*p0_over_gamma[2]*f[3]+1.414213562373095*p0_over_gamma[1]*f[2]+1.414213562373095*f[0]*p0_over_gamma[0])*volFact; 
  out[2] += (1.414213562373095*p1_over_gamma[7]*f[14]+1.414213562373095*p1_over_gamma[6]*f[10]+1.414213562373095*p1_over_gamma[5]*f[9]+1.414213562373095*p1_over_gamma[4]*f[7]+1.414213562373095*p1_over_gamma[3]*f[4]+1.414213562373095*p1_over_gamma[2]*f[3]+1.414213562373095*p1_over_gamma[1]*f[2]+1.414213562373095*f[0]*p1_over_gamma[0])*volFact; 
  out[3] += (1.414213562373095*p2_over_gamma[7]*f[14]+1.414213562373095*p2_over_gamma[6]*f[10]+1.414213562373095*p2_over_gamma[5]*f[9]+1.414213562373095*p2_over_gamma[4]*f[7]+1.414213562373095*p2_over_gamma[3]*f[4]+1.414213562373095*p2_over_gamma[2]*f[3]+1.414213562373095*p2_over_gamma[1]*f[2]+1.414213562373095*f[0]*p2_over_gamma[0])*volFact; 
  out[4] += volFact*(1.414213562373095*p2_over_gamma[7]*f[14]*wx3+1.414213562373095*p2_over_gamma[6]*f[10]*wx3+1.414213562373095*p2_over_gamma[5]*f[9]*wx3+1.414213562373095*p2_over_gamma[4]*f[7]*wx3+1.414213562373095*p2_over_gamma[3]*f[4]*wx3+1.414213562373095*p2_over_gamma[2]*f[3]*wx3+1.414213562373095*p2_over_gamma[1]*f[2]*wx3+1.414213562373095*f[0]*p2_over_gamma[0]*wx3+1.414213562373095*p1_over_gamma[7]*f[14]*wx2+1.414213562373095*p1_over_gamma[6]*f[10]*wx2+1.414213562373095*p1_over_gamma[5]*f[9]*wx2+1.414213562373095*p1_over_gamma[4]*f[7]*wx2+1.414213562373095*p1_over_gamma[3]*f[4]*wx2+1.414213562373095*p1_over_gamma[2]*f[3]*wx2+1.414213562373095*p1_over_gamma[1]*f[2]*wx2+1.414213562373095*f[0]*p1_over_gamma[0]*wx2+1.414213562373095*p0_over_gamma[7]*f[14]*wx1+1.414213562373095*p0_over_gamma[6]*f[10]*wx1+1.414213562373095*p0_over_gamma[5]*f[9]*wx1+1.414213562373095*p0_over_gamma[4]*f[7]*wx1+1.414213562373095*p0_over_gamma[3]*f[4]*wx1+1.414213562373095*p0_over_gamma[2]*f[3]*wx1+1.414213562373095*p0_over_gamma[1]*f[2]*wx1+1.414213562373095*f[0]*p0_over_gamma[0]*wx1+0.408248290463863*p2_over_gamma[4]*f[14]*dv3+0.408248290463863*p2_over_gamma[2]*f[10]*dv3+0.408248290463863*p2_over_gamma[1]*f[9]*dv3+0.408248290463863*f[7]*p2_over_gamma[7]*dv3+0.408248290463863*f[3]*p2_over_gamma[6]*dv3+0.408248290463863*f[2]*p2_over_gamma[5]*dv3+0.408248290463863*p2_over_gamma[0]*f[4]*dv3+0.408248290463863*f[0]*p2_over_gamma[3]*dv3+0.408248290463863*p1_over_gamma[5]*f[14]*dv2+0.408248290463863*p1_over_gamma[3]*f[10]*dv2+0.408248290463863*p1_over_gamma[7]*f[9]*dv2+0.408248290463863*p1_over_gamma[1]*f[7]*dv2+0.408248290463863*f[4]*p1_over_gamma[6]*dv2+0.408248290463863*f[2]*p1_over_gamma[4]*dv2+0.408248290463863*p1_over_gamma[0]*f[3]*dv2+0.408248290463863*f[0]*p1_over_gamma[2]*dv2+0.408248290463863*p0_over_gamma[6]*f[14]*dv1+0.408248290463863*p0_over_gamma[7]*f[10]*dv1+0.408248290463863*p0_over_gamma[3]*f[9]*dv1+0.408248290463863*p0_over_gamma[2]*f[7]*dv1+0.408248290463863*f[4]*p0_over_gamma[5]*dv1+0.408248290463863*f[3]*p0_over_gamma[4]*dv1+0.408248290463863*p0_over_gamma[0]*f[2]*dv1+0.408248290463863*f[0]*p0_over_gamma[1]*dv1); 
} 