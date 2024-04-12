#include <gkyl_canonical_pb_kernels.h> 
GKYL_CU_DH int canonical_pb_alpha_surfvx_2x2v_ser_p1(const double *w, const double *dxv, const double *hamil,
   double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // hamil: hamiltonian.
  // alpha_surf: output surface phase space flux in each direction (cdim + 1 components).
  //             Note: Each cell owns their *lower* edge surface evaluation.
  // sgn_alpha_surf: output sign(alpha_surf) in each direction at quadrature points (cdim + 1 components).
  //                 Note: Each cell owns their *lower* edge sign(alpha_surf).
  // returns int const_sgn_alpha (true if sign(alpha_surf) is only one sign, either +1 or -1).

  double wx = w[0];
  double rdx2 = 2.0/dxv[0];
  double wy = w[1];
  double rdy2 = 2.0/dxv[1];
  double wvx = w[2];
  double rdvx2 = 2.0/dxv[2];
  double wvy = w[3];
  double rdvy2 = 2.0/dxv[3];

  double *alphaL = &alpha_surf[16];
  double *sgn_alpha_surfL = &sgn_alpha_surf[16];
  alphaL[0] = (-2.738612787525831*hamil[17]*rdvx2*rdx2)+2.121320343559642*hamil[6]*rdvx2*rdx2-1.224744871391589*hamil[1]*rdvx2*rdx2; 
  alphaL[2] = (-2.738612787525831*hamil[20]*rdvx2*rdx2)+2.121320343559642*hamil[11]*rdvx2*rdx2-1.224744871391589*hamil[5]*rdvx2*rdx2; 
  alphaL[3] = (-2.738612787525831*hamil[21]*rdvx2*rdx2)+2.121320343559642*hamil[13]*rdvx2*rdx2-1.224744871391589*hamil[8]*rdvx2*rdx2; 
  alphaL[6] = (-2.738612787525831*hamil[23]*rdvx2*rdx2)+2.121320343559642*hamil[15]*rdvx2*rdx2-1.224744871391589*hamil[12]*rdvx2*rdx2; 

  int const_sgn_alpha_surf = 1;  
  
  if (0.3535533905932734*alphaL[6]-0.3535533905932734*(alphaL[3]+alphaL[2])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[0] = 1.0; 
  else  
    sgn_alpha_surfL[0] = -1.0; 
  
  if ((-0.3535533905932734*alphaL[6])+0.3535533905932734*alphaL[3]-0.3535533905932734*alphaL[2]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[1] = 1.0; 
  else  
    sgn_alpha_surfL[1] = -1.0; 
  
  if (sgn_alpha_surfL[1] == sgn_alpha_surfL[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*(alphaL[2]+alphaL[0])-0.3535533905932734*(alphaL[6]+alphaL[3]) > 0.) 
    sgn_alpha_surfL[2] = 1.0; 
  else  
    sgn_alpha_surfL[2] = -1.0; 
  
  if (sgn_alpha_surfL[2] == sgn_alpha_surfL[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*(alphaL[6]+alphaL[3]+alphaL[2]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[3] = 1.0; 
  else  
    sgn_alpha_surfL[3] = -1.0; 
  
  if (sgn_alpha_surfL[3] == sgn_alpha_surfL[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*alphaL[6]-0.3535533905932734*(alphaL[3]+alphaL[2])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[4] = 1.0; 
  else  
    sgn_alpha_surfL[4] = -1.0; 
  
  if (sgn_alpha_surfL[4] == sgn_alpha_surfL[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3535533905932734*alphaL[6])+0.3535533905932734*alphaL[3]-0.3535533905932734*alphaL[2]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[5] = 1.0; 
  else  
    sgn_alpha_surfL[5] = -1.0; 
  
  if (sgn_alpha_surfL[5] == sgn_alpha_surfL[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*(alphaL[2]+alphaL[0])-0.3535533905932734*(alphaL[6]+alphaL[3]) > 0.) 
    sgn_alpha_surfL[6] = 1.0; 
  else  
    sgn_alpha_surfL[6] = -1.0; 
  
  if (sgn_alpha_surfL[6] == sgn_alpha_surfL[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*(alphaL[6]+alphaL[3]+alphaL[2]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[7] = 1.0; 
  else  
    sgn_alpha_surfL[7] = -1.0; 
  
  if (sgn_alpha_surfL[7] == sgn_alpha_surfL[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
