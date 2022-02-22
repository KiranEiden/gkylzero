#include <gkyl_mom_gyrokinetic_kernels.h> 
GKYL_CU_DH void gyrokinetic_M0_3x2v_ser_p1(const double *w, const double *dxv, const int *idx, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[3]*dxv[4]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[3]*volFact; 
  out[4] += 2.0*f[6]*volFact; 
  out[5] += 2.0*f[7]*volFact; 
  out[6] += 2.0*f[8]*volFact; 
  out[7] += 2.0*f[16]*volFact; 
} 
GKYL_CU_DH void gyrokinetic_M1_3x2v_ser_p1(const double *w, const double *dxv, const int *idx, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  out[0] += volFact*(2.0*f[0]*wx1+0.5773502691896258*f[4]*dv1); 
  out[1] += volFact*(2.0*f[1]*wx1+0.5773502691896258*f[9]*dv1); 
  out[2] += volFact*(2.0*f[2]*wx1+0.5773502691896258*f[10]*dv1); 
  out[3] += volFact*(2.0*f[3]*wx1+0.5773502691896258*f[11]*dv1); 
  out[4] += volFact*(2.0*f[6]*wx1+0.5773502691896258*f[17]*dv1); 
  out[5] += volFact*(2.0*f[7]*wx1+0.5773502691896258*f[18]*dv1); 
  out[6] += volFact*(2.0*f[8]*wx1+0.5773502691896258*f[19]*dv1); 
  out[7] += volFact*(2.0*f[16]*wx1+0.5773502691896258*f[26]*dv1); 
} 
GKYL_CU_DH void gyrokinetic_M2_3x2v_ser_p1(const double *w, const double *dxv, const int *idx, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += volFact*(2.0*f[0]*wx1_sq+1.154700538379252*f[4]*dv1*wx1+0.1666666666666667*f[0]*dv1_sq); 
  out[1] += volFact*(2.0*f[1]*wx1_sq+1.154700538379252*f[9]*dv1*wx1+0.1666666666666667*f[1]*dv1_sq); 
  out[2] += volFact*(2.0*f[2]*wx1_sq+1.154700538379252*f[10]*dv1*wx1+0.1666666666666667*f[2]*dv1_sq); 
  out[3] += volFact*(2.0*f[3]*wx1_sq+1.154700538379252*f[11]*dv1*wx1+0.1666666666666667*f[3]*dv1_sq); 
  out[4] += volFact*(2.0*f[6]*wx1_sq+1.154700538379252*f[17]*dv1*wx1+0.1666666666666667*f[6]*dv1_sq); 
  out[5] += volFact*(2.0*f[7]*wx1_sq+1.154700538379252*f[18]*dv1*wx1+0.1666666666666667*f[7]*dv1_sq); 
  out[6] += volFact*(2.0*f[8]*wx1_sq+1.154700538379252*f[19]*dv1*wx1+0.1666666666666667*f[8]*dv1_sq); 
  out[7] += volFact*(2.0*f[16]*wx1_sq+1.154700538379252*f[26]*dv1*wx1+0.1666666666666667*f[16]*dv1_sq); 
  double tmp[8]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[5]*dv2; 
  tmp[1] = 2.0*f[1]*wx2+0.5773502691896258*f[12]*dv2; 
  tmp[2] = 2.0*f[2]*wx2+0.5773502691896258*f[13]*dv2; 
  tmp[3] = 2.0*f[3]*wx2+0.5773502691896258*f[14]*dv2; 
  tmp[4] = 2.0*f[6]*wx2+0.5773502691896258*f[20]*dv2; 
  tmp[5] = 2.0*f[7]*wx2+0.5773502691896258*f[21]*dv2; 
  tmp[6] = 2.0*f[8]*wx2+0.5773502691896258*f[22]*dv2; 
  tmp[7] = 2.0*f[16]*wx2+0.5773502691896258*f[27]*dv2; 
  out[0] += (2.0*(0.3535533905932737*bmag[7]*tmp[7]+0.3535533905932737*bmag[6]*tmp[6]+0.3535533905932737*bmag[5]*tmp[5]+0.3535533905932737*bmag[4]*tmp[4]+0.3535533905932737*bmag[3]*tmp[3]+0.3535533905932737*bmag[2]*tmp[2]+0.3535533905932737*bmag[1]*tmp[1]+0.3535533905932737*bmag[0]*tmp[0])*volFact)/m_; 
  out[1] += (2.0*(0.3535533905932737*bmag[6]*tmp[7]+0.3535533905932737*tmp[6]*bmag[7]+0.3535533905932737*bmag[3]*tmp[5]+0.3535533905932737*tmp[3]*bmag[5]+0.3535533905932737*bmag[2]*tmp[4]+0.3535533905932737*tmp[2]*bmag[4]+0.3535533905932737*bmag[0]*tmp[1]+0.3535533905932737*tmp[0]*bmag[1])*volFact)/m_; 
  out[2] += (2.0*(0.3535533905932737*bmag[5]*tmp[7]+0.3535533905932737*tmp[5]*bmag[7]+0.3535533905932737*bmag[3]*tmp[6]+0.3535533905932737*tmp[3]*bmag[6]+0.3535533905932737*bmag[1]*tmp[4]+0.3535533905932737*tmp[1]*bmag[4]+0.3535533905932737*bmag[0]*tmp[2]+0.3535533905932737*tmp[0]*bmag[2])*volFact)/m_; 
  out[3] += (2.0*(0.3535533905932737*bmag[4]*tmp[7]+0.3535533905932737*tmp[4]*bmag[7]+0.3535533905932737*bmag[2]*tmp[6]+0.3535533905932737*tmp[2]*bmag[6]+0.3535533905932737*bmag[1]*tmp[5]+0.3535533905932737*tmp[1]*bmag[5]+0.3535533905932737*bmag[0]*tmp[3]+0.3535533905932737*tmp[0]*bmag[3])*volFact)/m_; 
  out[4] += (2.0*(0.3535533905932737*bmag[3]*tmp[7]+0.3535533905932737*tmp[3]*bmag[7]+0.3535533905932737*bmag[5]*tmp[6]+0.3535533905932737*tmp[5]*bmag[6]+0.3535533905932737*bmag[0]*tmp[4]+0.3535533905932737*tmp[0]*bmag[4]+0.3535533905932737*bmag[1]*tmp[2]+0.3535533905932737*tmp[1]*bmag[2])*volFact)/m_; 
  out[5] += (2.0*(0.3535533905932737*bmag[2]*tmp[7]+0.3535533905932737*tmp[2]*bmag[7]+0.3535533905932737*bmag[4]*tmp[6]+0.3535533905932737*tmp[4]*bmag[6]+0.3535533905932737*bmag[0]*tmp[5]+0.3535533905932737*tmp[0]*bmag[5]+0.3535533905932737*bmag[1]*tmp[3]+0.3535533905932737*tmp[1]*bmag[3])*volFact)/m_; 
  out[6] += (2.0*(0.3535533905932737*bmag[1]*tmp[7]+0.3535533905932737*tmp[1]*bmag[7]+0.3535533905932737*bmag[0]*tmp[6]+0.3535533905932737*tmp[0]*bmag[6]+0.3535533905932737*bmag[4]*tmp[5]+0.3535533905932737*tmp[4]*bmag[5]+0.3535533905932737*bmag[2]*tmp[3]+0.3535533905932737*tmp[2]*bmag[3])*volFact)/m_; 
  out[7] += (2.0*(0.3535533905932737*bmag[0]*tmp[7]+0.3535533905932737*tmp[0]*bmag[7]+0.3535533905932737*bmag[1]*tmp[6]+0.3535533905932737*tmp[1]*bmag[6]+0.3535533905932737*bmag[2]*tmp[5]+0.3535533905932737*tmp[2]*bmag[5]+0.3535533905932737*bmag[3]*tmp[4]+0.3535533905932737*tmp[3]*bmag[4])*volFact)/m_; 
} 
GKYL_CU_DH void gyrokinetic_M2_par_3x2v_ser_p1(const double *w, const double *dxv, const int *idx, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += volFact*(2.0*f[0]*wx1_sq+1.154700538379252*f[4]*dv1*wx1+0.1666666666666667*f[0]*dv1_sq); 
  out[1] += volFact*(2.0*f[1]*wx1_sq+1.154700538379252*f[9]*dv1*wx1+0.1666666666666667*f[1]*dv1_sq); 
  out[2] += volFact*(2.0*f[2]*wx1_sq+1.154700538379252*f[10]*dv1*wx1+0.1666666666666667*f[2]*dv1_sq); 
  out[3] += volFact*(2.0*f[3]*wx1_sq+1.154700538379252*f[11]*dv1*wx1+0.1666666666666667*f[3]*dv1_sq); 
  out[4] += volFact*(2.0*f[6]*wx1_sq+1.154700538379252*f[17]*dv1*wx1+0.1666666666666667*f[6]*dv1_sq); 
  out[5] += volFact*(2.0*f[7]*wx1_sq+1.154700538379252*f[18]*dv1*wx1+0.1666666666666667*f[7]*dv1_sq); 
  out[6] += volFact*(2.0*f[8]*wx1_sq+1.154700538379252*f[19]*dv1*wx1+0.1666666666666667*f[8]*dv1_sq); 
  out[7] += volFact*(2.0*f[16]*wx1_sq+1.154700538379252*f[26]*dv1*wx1+0.1666666666666667*f[16]*dv1_sq); 
} 
GKYL_CU_DH void gyrokinetic_M2_perp_3x2v_ser_p1(const double *w, const double *dxv, const int *idx, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  double tmp[8]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[5]*dv2; 
  tmp[1] = 2.0*f[1]*wx2+0.5773502691896258*f[12]*dv2; 
  tmp[2] = 2.0*f[2]*wx2+0.5773502691896258*f[13]*dv2; 
  tmp[3] = 2.0*f[3]*wx2+0.5773502691896258*f[14]*dv2; 
  tmp[4] = 2.0*f[6]*wx2+0.5773502691896258*f[20]*dv2; 
  tmp[5] = 2.0*f[7]*wx2+0.5773502691896258*f[21]*dv2; 
  tmp[6] = 2.0*f[8]*wx2+0.5773502691896258*f[22]*dv2; 
  tmp[7] = 2.0*f[16]*wx2+0.5773502691896258*f[27]*dv2; 
  out[0] += ((0.3535533905932737*bmag[7]*tmp[7]+0.3535533905932737*bmag[6]*tmp[6]+0.3535533905932737*bmag[5]*tmp[5]+0.3535533905932737*bmag[4]*tmp[4]+0.3535533905932737*bmag[3]*tmp[3]+0.3535533905932737*bmag[2]*tmp[2]+0.3535533905932737*bmag[1]*tmp[1]+0.3535533905932737*bmag[0]*tmp[0])*volFact)/m_; 
  out[1] += ((0.3535533905932737*bmag[6]*tmp[7]+0.3535533905932737*tmp[6]*bmag[7]+0.3535533905932737*bmag[3]*tmp[5]+0.3535533905932737*tmp[3]*bmag[5]+0.3535533905932737*bmag[2]*tmp[4]+0.3535533905932737*tmp[2]*bmag[4]+0.3535533905932737*bmag[0]*tmp[1]+0.3535533905932737*tmp[0]*bmag[1])*volFact)/m_; 
  out[2] += ((0.3535533905932737*bmag[5]*tmp[7]+0.3535533905932737*tmp[5]*bmag[7]+0.3535533905932737*bmag[3]*tmp[6]+0.3535533905932737*tmp[3]*bmag[6]+0.3535533905932737*bmag[1]*tmp[4]+0.3535533905932737*tmp[1]*bmag[4]+0.3535533905932737*bmag[0]*tmp[2]+0.3535533905932737*tmp[0]*bmag[2])*volFact)/m_; 
  out[3] += ((0.3535533905932737*bmag[4]*tmp[7]+0.3535533905932737*tmp[4]*bmag[7]+0.3535533905932737*bmag[2]*tmp[6]+0.3535533905932737*tmp[2]*bmag[6]+0.3535533905932737*bmag[1]*tmp[5]+0.3535533905932737*tmp[1]*bmag[5]+0.3535533905932737*bmag[0]*tmp[3]+0.3535533905932737*tmp[0]*bmag[3])*volFact)/m_; 
  out[4] += ((0.3535533905932737*bmag[3]*tmp[7]+0.3535533905932737*tmp[3]*bmag[7]+0.3535533905932737*bmag[5]*tmp[6]+0.3535533905932737*tmp[5]*bmag[6]+0.3535533905932737*bmag[0]*tmp[4]+0.3535533905932737*tmp[0]*bmag[4]+0.3535533905932737*bmag[1]*tmp[2]+0.3535533905932737*tmp[1]*bmag[2])*volFact)/m_; 
  out[5] += ((0.3535533905932737*bmag[2]*tmp[7]+0.3535533905932737*tmp[2]*bmag[7]+0.3535533905932737*bmag[4]*tmp[6]+0.3535533905932737*tmp[4]*bmag[6]+0.3535533905932737*bmag[0]*tmp[5]+0.3535533905932737*tmp[0]*bmag[5]+0.3535533905932737*bmag[1]*tmp[3]+0.3535533905932737*tmp[1]*bmag[3])*volFact)/m_; 
  out[6] += ((0.3535533905932737*bmag[1]*tmp[7]+0.3535533905932737*tmp[1]*bmag[7]+0.3535533905932737*bmag[0]*tmp[6]+0.3535533905932737*tmp[0]*bmag[6]+0.3535533905932737*bmag[4]*tmp[5]+0.3535533905932737*tmp[4]*bmag[5]+0.3535533905932737*bmag[2]*tmp[3]+0.3535533905932737*tmp[2]*bmag[3])*volFact)/m_; 
  out[7] += ((0.3535533905932737*bmag[0]*tmp[7]+0.3535533905932737*tmp[0]*bmag[7]+0.3535533905932737*bmag[1]*tmp[6]+0.3535533905932737*tmp[1]*bmag[6]+0.3535533905932737*bmag[2]*tmp[5]+0.3535533905932737*tmp[2]*bmag[5]+0.3535533905932737*bmag[3]*tmp[4]+0.3535533905932737*tmp[3]*bmag[4])*volFact)/m_; 
} 
GKYL_CU_DH void gyrokinetic_M3_par_3x2v_ser_p1(const double *w, const double *dxv, const int *idx, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  out[0] += volFact*(2.0*f[0]*wx1*wx1_sq+1.732050807568877*f[4]*dv1*wx1_sq+0.5*f[0]*dv1_sq*wx1+0.08660254037844387*f[4]*dv1*dv1_sq); 
  out[1] += volFact*(2.0*f[1]*wx1*wx1_sq+1.732050807568877*f[9]*dv1*wx1_sq+0.5*f[1]*dv1_sq*wx1+0.08660254037844387*f[9]*dv1*dv1_sq); 
  out[2] += volFact*(2.0*f[2]*wx1*wx1_sq+1.732050807568877*f[10]*dv1*wx1_sq+0.5*f[2]*dv1_sq*wx1+0.08660254037844387*f[10]*dv1*dv1_sq); 
  out[3] += volFact*(2.0*f[3]*wx1*wx1_sq+1.732050807568877*f[11]*dv1*wx1_sq+0.5*f[3]*dv1_sq*wx1+0.08660254037844387*f[11]*dv1*dv1_sq); 
  out[4] += volFact*(2.0*f[6]*wx1*wx1_sq+1.732050807568877*f[17]*dv1*wx1_sq+0.5*f[6]*dv1_sq*wx1+0.08660254037844387*f[17]*dv1*dv1_sq); 
  out[5] += volFact*(2.0*f[7]*wx1*wx1_sq+1.732050807568877*f[18]*dv1*wx1_sq+0.5*f[7]*dv1_sq*wx1+0.08660254037844387*f[18]*dv1*dv1_sq); 
  out[6] += volFact*(2.0*f[8]*wx1*wx1_sq+1.732050807568877*f[19]*dv1*wx1_sq+0.5*f[8]*dv1_sq*wx1+0.08660254037844387*f[19]*dv1*dv1_sq); 
  out[7] += volFact*(2.0*f[16]*wx1*wx1_sq+1.732050807568877*f[26]*dv1*wx1_sq+0.5*f[16]*dv1_sq*wx1+0.08660254037844387*f[26]*dv1*dv1_sq); 
} 
GKYL_CU_DH void gyrokinetic_M3_perp_3x2v_ser_p1(const double *w, const double *dxv, const int *idx, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  out[0] += (volFact*(0.7071067811865475*bmag[7]*f[16]*wx1*wx2+0.7071067811865475*bmag[6]*f[8]*wx1*wx2+0.7071067811865475*bmag[5]*f[7]*wx1*wx2+0.7071067811865475*bmag[4]*f[6]*wx1*wx2+0.7071067811865475*bmag[3]*f[3]*wx1*wx2+0.7071067811865475*bmag[2]*f[2]*wx1*wx2+0.7071067811865475*bmag[1]*f[1]*wx1*wx2+0.7071067811865475*bmag[0]*f[0]*wx1*wx2+0.2041241452319315*bmag[7]*f[26]*dv1*wx2+0.2041241452319315*bmag[6]*f[19]*dv1*wx2+0.2041241452319315*bmag[5]*f[18]*dv1*wx2+0.2041241452319315*bmag[4]*f[17]*dv1*wx2+0.2041241452319315*bmag[3]*f[11]*dv1*wx2+0.2041241452319315*bmag[2]*f[10]*dv1*wx2+0.2041241452319315*bmag[1]*f[9]*dv1*wx2+0.2041241452319315*bmag[0]*f[4]*dv1*wx2+0.2041241452319315*bmag[7]*f[27]*dv2*wx1+0.2041241452319315*bmag[6]*f[22]*dv2*wx1+0.2041241452319315*bmag[5]*f[21]*dv2*wx1+0.2041241452319315*bmag[4]*f[20]*dv2*wx1+0.2041241452319315*bmag[3]*f[14]*dv2*wx1+0.2041241452319315*bmag[2]*f[13]*dv2*wx1+0.2041241452319315*bmag[1]*f[12]*dv2*wx1+0.2041241452319315*bmag[0]*f[5]*dv2*wx1+0.05892556509887893*bmag[7]*f[31]*dv1*dv2+0.05892556509887893*bmag[6]*f[30]*dv1*dv2+0.05892556509887893*bmag[5]*f[29]*dv1*dv2+0.05892556509887893*bmag[4]*f[28]*dv1*dv2+0.05892556509887893*bmag[3]*f[25]*dv1*dv2+0.05892556509887893*bmag[2]*f[24]*dv1*dv2+0.05892556509887893*bmag[1]*f[23]*dv1*dv2+0.05892556509887893*bmag[0]*f[15]*dv1*dv2))/m_; 
  out[1] += (volFact*(0.7071067811865475*bmag[6]*f[16]*wx1*wx2+0.7071067811865475*bmag[7]*f[8]*wx1*wx2+0.7071067811865475*bmag[3]*f[7]*wx1*wx2+0.7071067811865475*bmag[2]*f[6]*wx1*wx2+0.7071067811865475*f[3]*bmag[5]*wx1*wx2+0.7071067811865475*f[2]*bmag[4]*wx1*wx2+0.7071067811865475*bmag[0]*f[1]*wx1*wx2+0.7071067811865475*f[0]*bmag[1]*wx1*wx2+0.2041241452319315*bmag[6]*f[26]*dv1*wx2+0.2041241452319315*bmag[7]*f[19]*dv1*wx2+0.2041241452319315*bmag[3]*f[18]*dv1*wx2+0.2041241452319315*bmag[2]*f[17]*dv1*wx2+0.2041241452319315*bmag[5]*f[11]*dv1*wx2+0.2041241452319315*bmag[4]*f[10]*dv1*wx2+0.2041241452319315*bmag[0]*f[9]*dv1*wx2+0.2041241452319315*bmag[1]*f[4]*dv1*wx2+0.2041241452319315*bmag[6]*f[27]*dv2*wx1+0.2041241452319315*bmag[7]*f[22]*dv2*wx1+0.2041241452319315*bmag[3]*f[21]*dv2*wx1+0.2041241452319315*bmag[2]*f[20]*dv2*wx1+0.2041241452319315*bmag[5]*f[14]*dv2*wx1+0.2041241452319315*bmag[4]*f[13]*dv2*wx1+0.2041241452319315*bmag[0]*f[12]*dv2*wx1+0.2041241452319315*bmag[1]*f[5]*dv2*wx1+0.05892556509887893*bmag[6]*f[31]*dv1*dv2+0.05892556509887893*bmag[7]*f[30]*dv1*dv2+0.05892556509887893*bmag[3]*f[29]*dv1*dv2+0.05892556509887893*bmag[2]*f[28]*dv1*dv2+0.05892556509887893*bmag[5]*f[25]*dv1*dv2+0.05892556509887893*bmag[4]*f[24]*dv1*dv2+0.05892556509887893*bmag[0]*f[23]*dv1*dv2+0.05892556509887893*bmag[1]*f[15]*dv1*dv2))/m_; 
  out[2] += (volFact*(0.7071067811865475*bmag[5]*f[16]*wx1*wx2+0.7071067811865475*bmag[3]*f[8]*wx1*wx2+0.7071067811865475*bmag[7]*f[7]*wx1*wx2+0.7071067811865475*bmag[1]*f[6]*wx1*wx2+0.7071067811865475*f[3]*bmag[6]*wx1*wx2+0.7071067811865475*f[1]*bmag[4]*wx1*wx2+0.7071067811865475*bmag[0]*f[2]*wx1*wx2+0.7071067811865475*f[0]*bmag[2]*wx1*wx2+0.2041241452319315*bmag[5]*f[26]*dv1*wx2+0.2041241452319315*bmag[3]*f[19]*dv1*wx2+0.2041241452319315*bmag[7]*f[18]*dv1*wx2+0.2041241452319315*bmag[1]*f[17]*dv1*wx2+0.2041241452319315*bmag[6]*f[11]*dv1*wx2+0.2041241452319315*bmag[0]*f[10]*dv1*wx2+0.2041241452319315*bmag[4]*f[9]*dv1*wx2+0.2041241452319315*bmag[2]*f[4]*dv1*wx2+0.2041241452319315*bmag[5]*f[27]*dv2*wx1+0.2041241452319315*bmag[3]*f[22]*dv2*wx1+0.2041241452319315*bmag[7]*f[21]*dv2*wx1+0.2041241452319315*bmag[1]*f[20]*dv2*wx1+0.2041241452319315*bmag[6]*f[14]*dv2*wx1+0.2041241452319315*bmag[0]*f[13]*dv2*wx1+0.2041241452319315*bmag[4]*f[12]*dv2*wx1+0.2041241452319315*bmag[2]*f[5]*dv2*wx1+0.05892556509887893*bmag[5]*f[31]*dv1*dv2+0.05892556509887893*bmag[3]*f[30]*dv1*dv2+0.05892556509887893*bmag[7]*f[29]*dv1*dv2+0.05892556509887893*bmag[1]*f[28]*dv1*dv2+0.05892556509887893*bmag[6]*f[25]*dv1*dv2+0.05892556509887893*bmag[0]*f[24]*dv1*dv2+0.05892556509887893*bmag[4]*f[23]*dv1*dv2+0.05892556509887893*bmag[2]*f[15]*dv1*dv2))/m_; 
  out[3] += (volFact*(0.7071067811865475*bmag[4]*f[16]*wx1*wx2+0.7071067811865475*bmag[2]*f[8]*wx1*wx2+0.7071067811865475*bmag[1]*f[7]*wx1*wx2+0.7071067811865475*f[6]*bmag[7]*wx1*wx2+0.7071067811865475*f[2]*bmag[6]*wx1*wx2+0.7071067811865475*f[1]*bmag[5]*wx1*wx2+0.7071067811865475*bmag[0]*f[3]*wx1*wx2+0.7071067811865475*f[0]*bmag[3]*wx1*wx2+0.2041241452319315*bmag[4]*f[26]*dv1*wx2+0.2041241452319315*bmag[2]*f[19]*dv1*wx2+0.2041241452319315*bmag[1]*f[18]*dv1*wx2+0.2041241452319315*bmag[7]*f[17]*dv1*wx2+0.2041241452319315*bmag[0]*f[11]*dv1*wx2+0.2041241452319315*bmag[6]*f[10]*dv1*wx2+0.2041241452319315*bmag[5]*f[9]*dv1*wx2+0.2041241452319315*bmag[3]*f[4]*dv1*wx2+0.2041241452319315*bmag[4]*f[27]*dv2*wx1+0.2041241452319315*bmag[2]*f[22]*dv2*wx1+0.2041241452319315*bmag[1]*f[21]*dv2*wx1+0.2041241452319315*bmag[7]*f[20]*dv2*wx1+0.2041241452319315*bmag[0]*f[14]*dv2*wx1+0.2041241452319315*bmag[6]*f[13]*dv2*wx1+0.2041241452319315*bmag[5]*f[12]*dv2*wx1+0.2041241452319315*bmag[3]*f[5]*dv2*wx1+0.05892556509887893*bmag[4]*f[31]*dv1*dv2+0.05892556509887893*bmag[2]*f[30]*dv1*dv2+0.05892556509887893*bmag[1]*f[29]*dv1*dv2+0.05892556509887893*bmag[7]*f[28]*dv1*dv2+0.05892556509887893*bmag[0]*f[25]*dv1*dv2+0.05892556509887893*bmag[6]*f[24]*dv1*dv2+0.05892556509887893*bmag[5]*f[23]*dv1*dv2+0.05892556509887893*bmag[3]*f[15]*dv1*dv2))/m_; 
  out[4] += (volFact*(0.7071067811865475*bmag[3]*f[16]*wx1*wx2+0.7071067811865475*bmag[5]*f[8]*wx1*wx2+0.7071067811865475*bmag[6]*f[7]*wx1*wx2+0.7071067811865475*f[3]*bmag[7]*wx1*wx2+0.7071067811865475*bmag[0]*f[6]*wx1*wx2+0.7071067811865475*f[0]*bmag[4]*wx1*wx2+0.7071067811865475*bmag[1]*f[2]*wx1*wx2+0.7071067811865475*f[1]*bmag[2]*wx1*wx2+0.2041241452319315*bmag[3]*f[26]*dv1*wx2+0.2041241452319315*bmag[5]*f[19]*dv1*wx2+0.2041241452319315*bmag[6]*f[18]*dv1*wx2+0.2041241452319315*bmag[0]*f[17]*dv1*wx2+0.2041241452319315*bmag[7]*f[11]*dv1*wx2+0.2041241452319315*bmag[1]*f[10]*dv1*wx2+0.2041241452319315*bmag[2]*f[9]*dv1*wx2+0.2041241452319315*bmag[4]*f[4]*dv1*wx2+0.2041241452319315*bmag[3]*f[27]*dv2*wx1+0.2041241452319315*bmag[5]*f[22]*dv2*wx1+0.2041241452319315*bmag[6]*f[21]*dv2*wx1+0.2041241452319315*bmag[0]*f[20]*dv2*wx1+0.2041241452319315*bmag[7]*f[14]*dv2*wx1+0.2041241452319315*bmag[1]*f[13]*dv2*wx1+0.2041241452319315*bmag[2]*f[12]*dv2*wx1+0.2041241452319315*bmag[4]*f[5]*dv2*wx1+0.05892556509887893*bmag[3]*f[31]*dv1*dv2+0.05892556509887893*bmag[5]*f[30]*dv1*dv2+0.05892556509887893*bmag[6]*f[29]*dv1*dv2+0.05892556509887893*bmag[0]*f[28]*dv1*dv2+0.05892556509887893*bmag[7]*f[25]*dv1*dv2+0.05892556509887893*bmag[1]*f[24]*dv1*dv2+0.05892556509887893*bmag[2]*f[23]*dv1*dv2+0.05892556509887893*bmag[4]*f[15]*dv1*dv2))/m_; 
  out[5] += (volFact*(0.7071067811865475*bmag[2]*f[16]*wx1*wx2+0.7071067811865475*bmag[4]*f[8]*wx1*wx2+0.7071067811865475*bmag[0]*f[7]*wx1*wx2+0.7071067811865475*f[2]*bmag[7]*wx1*wx2+0.7071067811865475*bmag[6]*f[6]*wx1*wx2+0.7071067811865475*f[0]*bmag[5]*wx1*wx2+0.7071067811865475*bmag[1]*f[3]*wx1*wx2+0.7071067811865475*f[1]*bmag[3]*wx1*wx2+0.2041241452319315*bmag[2]*f[26]*dv1*wx2+0.2041241452319315*bmag[4]*f[19]*dv1*wx2+0.2041241452319315*bmag[0]*f[18]*dv1*wx2+0.2041241452319315*bmag[6]*f[17]*dv1*wx2+0.2041241452319315*bmag[1]*f[11]*dv1*wx2+0.2041241452319315*bmag[7]*f[10]*dv1*wx2+0.2041241452319315*bmag[3]*f[9]*dv1*wx2+0.2041241452319315*f[4]*bmag[5]*dv1*wx2+0.2041241452319315*bmag[2]*f[27]*dv2*wx1+0.2041241452319315*bmag[4]*f[22]*dv2*wx1+0.2041241452319315*bmag[0]*f[21]*dv2*wx1+0.2041241452319315*bmag[6]*f[20]*dv2*wx1+0.2041241452319315*bmag[1]*f[14]*dv2*wx1+0.2041241452319315*bmag[7]*f[13]*dv2*wx1+0.2041241452319315*bmag[3]*f[12]*dv2*wx1+0.2041241452319315*bmag[5]*f[5]*dv2*wx1+0.05892556509887893*bmag[2]*f[31]*dv1*dv2+0.05892556509887893*bmag[4]*f[30]*dv1*dv2+0.05892556509887893*bmag[0]*f[29]*dv1*dv2+0.05892556509887893*bmag[6]*f[28]*dv1*dv2+0.05892556509887893*bmag[1]*f[25]*dv1*dv2+0.05892556509887893*bmag[7]*f[24]*dv1*dv2+0.05892556509887893*bmag[3]*f[23]*dv1*dv2+0.05892556509887893*bmag[5]*f[15]*dv1*dv2))/m_; 
  out[6] += (volFact*(0.7071067811865475*bmag[1]*f[16]*wx1*wx2+0.7071067811865475*bmag[0]*f[8]*wx1*wx2+0.7071067811865475*bmag[4]*f[7]*wx1*wx2+0.7071067811865475*f[1]*bmag[7]*wx1*wx2+0.7071067811865475*bmag[5]*f[6]*wx1*wx2+0.7071067811865475*f[0]*bmag[6]*wx1*wx2+0.7071067811865475*bmag[2]*f[3]*wx1*wx2+0.7071067811865475*f[2]*bmag[3]*wx1*wx2+0.2041241452319315*bmag[1]*f[26]*dv1*wx2+0.2041241452319315*bmag[0]*f[19]*dv1*wx2+0.2041241452319315*bmag[4]*f[18]*dv1*wx2+0.2041241452319315*bmag[5]*f[17]*dv1*wx2+0.2041241452319315*bmag[2]*f[11]*dv1*wx2+0.2041241452319315*bmag[3]*f[10]*dv1*wx2+0.2041241452319315*bmag[7]*f[9]*dv1*wx2+0.2041241452319315*f[4]*bmag[6]*dv1*wx2+0.2041241452319315*bmag[1]*f[27]*dv2*wx1+0.2041241452319315*bmag[0]*f[22]*dv2*wx1+0.2041241452319315*bmag[4]*f[21]*dv2*wx1+0.2041241452319315*bmag[5]*f[20]*dv2*wx1+0.2041241452319315*bmag[2]*f[14]*dv2*wx1+0.2041241452319315*bmag[3]*f[13]*dv2*wx1+0.2041241452319315*bmag[7]*f[12]*dv2*wx1+0.2041241452319315*f[5]*bmag[6]*dv2*wx1+0.05892556509887893*bmag[1]*f[31]*dv1*dv2+0.05892556509887893*bmag[0]*f[30]*dv1*dv2+0.05892556509887893*bmag[4]*f[29]*dv1*dv2+0.05892556509887893*bmag[5]*f[28]*dv1*dv2+0.05892556509887893*bmag[2]*f[25]*dv1*dv2+0.05892556509887893*bmag[3]*f[24]*dv1*dv2+0.05892556509887893*bmag[7]*f[23]*dv1*dv2+0.05892556509887893*bmag[6]*f[15]*dv1*dv2))/m_; 
  out[7] += (volFact*(0.7071067811865475*bmag[0]*f[16]*wx1*wx2+0.7071067811865475*bmag[1]*f[8]*wx1*wx2+0.7071067811865475*bmag[2]*f[7]*wx1*wx2+0.7071067811865475*f[0]*bmag[7]*wx1*wx2+0.7071067811865475*bmag[3]*f[6]*wx1*wx2+0.7071067811865475*f[1]*bmag[6]*wx1*wx2+0.7071067811865475*f[2]*bmag[5]*wx1*wx2+0.7071067811865475*f[3]*bmag[4]*wx1*wx2+0.2041241452319315*bmag[0]*f[26]*dv1*wx2+0.2041241452319315*bmag[1]*f[19]*dv1*wx2+0.2041241452319315*bmag[2]*f[18]*dv1*wx2+0.2041241452319315*bmag[3]*f[17]*dv1*wx2+0.2041241452319315*bmag[4]*f[11]*dv1*wx2+0.2041241452319315*bmag[5]*f[10]*dv1*wx2+0.2041241452319315*bmag[6]*f[9]*dv1*wx2+0.2041241452319315*f[4]*bmag[7]*dv1*wx2+0.2041241452319315*bmag[0]*f[27]*dv2*wx1+0.2041241452319315*bmag[1]*f[22]*dv2*wx1+0.2041241452319315*bmag[2]*f[21]*dv2*wx1+0.2041241452319315*bmag[3]*f[20]*dv2*wx1+0.2041241452319315*bmag[4]*f[14]*dv2*wx1+0.2041241452319315*bmag[5]*f[13]*dv2*wx1+0.2041241452319315*bmag[6]*f[12]*dv2*wx1+0.2041241452319315*f[5]*bmag[7]*dv2*wx1+0.05892556509887893*bmag[0]*f[31]*dv1*dv2+0.05892556509887893*bmag[1]*f[30]*dv1*dv2+0.05892556509887893*bmag[2]*f[29]*dv1*dv2+0.05892556509887893*bmag[3]*f[28]*dv1*dv2+0.05892556509887893*bmag[4]*f[25]*dv1*dv2+0.05892556509887893*bmag[5]*f[24]*dv1*dv2+0.05892556509887893*bmag[6]*f[23]*dv1*dv2+0.05892556509887893*bmag[7]*f[15]*dv1*dv2))/m_; 
} 
GKYL_CU_DH void gyrokinetic_ThreeMoments_3x2v_ser_p1(const double *w, const double *dxv, const int *idx, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT outM0, double* GKYL_RESTRICT outM1, double* GKYL_RESTRICT outM2) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  outM0[0] += 2.0*f[0]*volFact; 
  outM0[1] += 2.0*f[1]*volFact; 
  outM0[2] += 2.0*f[2]*volFact; 
  outM0[3] += 2.0*f[3]*volFact; 
  outM0[4] += 2.0*f[6]*volFact; 
  outM0[5] += 2.0*f[7]*volFact; 
  outM0[6] += 2.0*f[8]*volFact; 
  outM0[7] += 2.0*f[16]*volFact; 
  outM1[0] += 0.3333333333333333*volFact*(6.0*f[0]*wx1+1.732050807568877*f[4]*dv1); 
  outM1[1] += 0.3333333333333333*volFact*(6.0*f[1]*wx1+1.732050807568877*f[9]*dv1); 
  outM1[2] += 0.3333333333333333*volFact*(6.0*f[2]*wx1+1.732050807568877*f[10]*dv1); 
  outM1[3] += 0.3333333333333333*volFact*(6.0*f[3]*wx1+1.732050807568877*f[11]*dv1); 
  outM1[4] += 0.3333333333333333*volFact*(6.0*f[6]*wx1+1.732050807568877*f[17]*dv1); 
  outM1[5] += 0.3333333333333333*volFact*(6.0*f[7]*wx1+1.732050807568877*f[18]*dv1); 
  outM1[6] += 0.3333333333333333*volFact*(6.0*f[8]*wx1+1.732050807568877*f[19]*dv1); 
  outM1[7] += 0.3333333333333333*volFact*(6.0*f[16]*wx1+1.732050807568877*f[26]*dv1); 
  outM2[0] += volFact*(2.0*f[0]*wx1_sq+1.154700538379252*f[4]*dv1*wx1+0.1666666666666667*f[0]*dv1_sq); 
  outM2[1] += volFact*(2.0*f[1]*wx1_sq+1.154700538379252*f[9]*dv1*wx1+0.1666666666666667*f[1]*dv1_sq); 
  outM2[2] += volFact*(2.0*f[2]*wx1_sq+1.154700538379252*f[10]*dv1*wx1+0.1666666666666667*f[2]*dv1_sq); 
  outM2[3] += volFact*(2.0*f[3]*wx1_sq+1.154700538379252*f[11]*dv1*wx1+0.1666666666666667*f[3]*dv1_sq); 
  outM2[4] += volFact*(2.0*f[6]*wx1_sq+1.154700538379252*f[17]*dv1*wx1+0.1666666666666667*f[6]*dv1_sq); 
  outM2[5] += volFact*(2.0*f[7]*wx1_sq+1.154700538379252*f[18]*dv1*wx1+0.1666666666666667*f[7]*dv1_sq); 
  outM2[6] += volFact*(2.0*f[8]*wx1_sq+1.154700538379252*f[19]*dv1*wx1+0.1666666666666667*f[8]*dv1_sq); 
  outM2[7] += volFact*(2.0*f[16]*wx1_sq+1.154700538379252*f[26]*dv1*wx1+0.1666666666666667*f[16]*dv1_sq); 
  double tmp[8]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[5]*dv2; 
  tmp[1] = 2.0*f[1]*wx2+0.5773502691896258*f[12]*dv2; 
  tmp[2] = 2.0*f[2]*wx2+0.5773502691896258*f[13]*dv2; 
  tmp[3] = 2.0*f[3]*wx2+0.5773502691896258*f[14]*dv2; 
  tmp[4] = 2.0*f[6]*wx2+0.5773502691896258*f[20]*dv2; 
  tmp[5] = 2.0*f[7]*wx2+0.5773502691896258*f[21]*dv2; 
  tmp[6] = 2.0*f[8]*wx2+0.5773502691896258*f[22]*dv2; 
  tmp[7] = 2.0*f[16]*wx2+0.5773502691896258*f[27]*dv2; 
  outM2[0] += (2.0*(0.3535533905932737*bmag[7]*tmp[7]+0.3535533905932737*bmag[6]*tmp[6]+0.3535533905932737*bmag[5]*tmp[5]+0.3535533905932737*bmag[4]*tmp[4]+0.3535533905932737*bmag[3]*tmp[3]+0.3535533905932737*bmag[2]*tmp[2]+0.3535533905932737*bmag[1]*tmp[1]+0.3535533905932737*bmag[0]*tmp[0])*volFact)/m_; 
  outM2[1] += (2.0*(0.3535533905932737*bmag[6]*tmp[7]+0.3535533905932737*tmp[6]*bmag[7]+0.3535533905932737*bmag[3]*tmp[5]+0.3535533905932737*tmp[3]*bmag[5]+0.3535533905932737*bmag[2]*tmp[4]+0.3535533905932737*tmp[2]*bmag[4]+0.3535533905932737*bmag[0]*tmp[1]+0.3535533905932737*tmp[0]*bmag[1])*volFact)/m_; 
  outM2[2] += (2.0*(0.3535533905932737*bmag[5]*tmp[7]+0.3535533905932737*tmp[5]*bmag[7]+0.3535533905932737*bmag[3]*tmp[6]+0.3535533905932737*tmp[3]*bmag[6]+0.3535533905932737*bmag[1]*tmp[4]+0.3535533905932737*tmp[1]*bmag[4]+0.3535533905932737*bmag[0]*tmp[2]+0.3535533905932737*tmp[0]*bmag[2])*volFact)/m_; 
  outM2[3] += (2.0*(0.3535533905932737*bmag[4]*tmp[7]+0.3535533905932737*tmp[4]*bmag[7]+0.3535533905932737*bmag[2]*tmp[6]+0.3535533905932737*tmp[2]*bmag[6]+0.3535533905932737*bmag[1]*tmp[5]+0.3535533905932737*tmp[1]*bmag[5]+0.3535533905932737*bmag[0]*tmp[3]+0.3535533905932737*tmp[0]*bmag[3])*volFact)/m_; 
  outM2[4] += (2.0*(0.3535533905932737*bmag[3]*tmp[7]+0.3535533905932737*tmp[3]*bmag[7]+0.3535533905932737*bmag[5]*tmp[6]+0.3535533905932737*tmp[5]*bmag[6]+0.3535533905932737*bmag[0]*tmp[4]+0.3535533905932737*tmp[0]*bmag[4]+0.3535533905932737*bmag[1]*tmp[2]+0.3535533905932737*tmp[1]*bmag[2])*volFact)/m_; 
  outM2[5] += (2.0*(0.3535533905932737*bmag[2]*tmp[7]+0.3535533905932737*tmp[2]*bmag[7]+0.3535533905932737*bmag[4]*tmp[6]+0.3535533905932737*tmp[4]*bmag[6]+0.3535533905932737*bmag[0]*tmp[5]+0.3535533905932737*tmp[0]*bmag[5]+0.3535533905932737*bmag[1]*tmp[3]+0.3535533905932737*tmp[1]*bmag[3])*volFact)/m_; 
  outM2[6] += (2.0*(0.3535533905932737*bmag[1]*tmp[7]+0.3535533905932737*tmp[1]*bmag[7]+0.3535533905932737*bmag[0]*tmp[6]+0.3535533905932737*tmp[0]*bmag[6]+0.3535533905932737*bmag[4]*tmp[5]+0.3535533905932737*tmp[4]*bmag[5]+0.3535533905932737*bmag[2]*tmp[3]+0.3535533905932737*tmp[2]*bmag[3])*volFact)/m_; 
  outM2[7] += (2.0*(0.3535533905932737*bmag[0]*tmp[7]+0.3535533905932737*tmp[0]*bmag[7]+0.3535533905932737*bmag[1]*tmp[6]+0.3535533905932737*tmp[1]*bmag[6]+0.3535533905932737*bmag[2]*tmp[5]+0.3535533905932737*tmp[2]*bmag[5]+0.3535533905932737*bmag[3]*tmp[4]+0.3535533905932737*tmp[3]*bmag[4])*volFact)/m_; 
} 
GKYL_CU_DH void gyrokinetic_M0_step1_3x2v_ser_p1(const double *w, const double *dxv, const int *idx, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[3]/2; 
  out[0] += 1.414213562373095*f[0]*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact; 
  out[2] += 1.414213562373095*f[2]*volFact; 
  out[3] += 1.414213562373095*f[3]*volFact; 
  out[4] += 1.414213562373095*f[5]*volFact; 
  out[5] += 1.414213562373095*f[6]*volFact; 
  out[6] += 1.414213562373095*f[7]*volFact; 
  out[7] += 1.414213562373095*f[8]*volFact; 
  out[8] += 1.414213562373095*f[12]*volFact; 
  out[9] += 1.414213562373095*f[13]*volFact; 
  out[10] += 1.414213562373095*f[14]*volFact; 
  out[11] += 1.414213562373095*f[16]*volFact; 
  out[12] += 1.414213562373095*f[20]*volFact; 
  out[13] += 1.414213562373095*f[21]*volFact; 
  out[14] += 1.414213562373095*f[22]*volFact; 
  out[15] += 1.414213562373095*f[27]*volFact; 
} 
GKYL_CU_DH void gyrokinetic_M0_step2_3x2v_ser_p1(const double *w, const double *dxv, const int *idx, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[4]/2; 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact; 
  out[2] += 2.828427124746191*f[2]*volFact; 
  out[3] += 2.828427124746191*f[3]*volFact; 
  out[4] += 2.828427124746191*f[5]*volFact; 
  out[5] += 2.828427124746191*f[6]*volFact; 
  out[6] += 2.828427124746191*f[7]*volFact; 
  out[7] += 2.828427124746191*f[11]*volFact; 
} 
