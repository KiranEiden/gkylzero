#include <gkyl_canonical_pb_kernels.h> 
GKYL_CU_DH int canonical_pb_alpha_edge_surfy_3x3v_ser_p1(const double *w, const double *dxv, const double *hamil,
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
  double wz = w[2];
  double rdz2 = 2.0/dxv[2];
  double wvx = w[3];
  double rdvx2 = 2.0/dxv[3];
  double wvy = w[4];
  double rdvy2 = 2.0/dxv[4];
  double wvz = w[5];
  double rdvz2 = 2.0/dxv[5];

  double *alphaR = &alpha_surf[32];
  double *sgn_alpha_surfR = &sgn_alpha_surf[32];
  alphaR[0] = 2.121320343559642*hamil[14]*rdvy2*rdy2+1.224744871391589*hamil[5]*rdvy2*rdy2; 
  alphaR[1] = 2.121320343559642*hamil[26]*rdvy2*rdy2+1.224744871391589*hamil[13]*rdvy2*rdy2; 
  alphaR[2] = 2.121320343559642*hamil[28]*rdvy2*rdy2+1.224744871391589*hamil[15]*rdvy2*rdy2; 
  alphaR[3] = 2.121320343559642*hamil[30]*rdvy2*rdy2+1.224744871391589*hamil[16]*rdvy2*rdy2; 
  alphaR[4] = 4.743416490252569*hamil[98]*rdvy2*rdy2+2.738612787525831*hamil[96]*rdvy2*rdy2; 
  alphaR[5] = 2.121320343559642*hamil[39]*rdvy2*rdy2+1.224744871391589*hamil[21]*rdvy2*rdy2; 
  alphaR[6] = 2.121320343559642*hamil[43]*rdvy2*rdy2+1.224744871391589*hamil[27]*rdvy2*rdy2; 
  alphaR[7] = 2.121320343559642*hamil[44]*rdvy2*rdy2+1.224744871391589*hamil[29]*rdvy2*rdy2; 
  alphaR[8] = 2.121320343559642*hamil[46]*rdvy2*rdy2+1.224744871391589*hamil[31]*rdvy2*rdy2; 
  alphaR[9] = 4.743416490252569*hamil[102]*rdvy2*rdy2+2.738612787525831*hamil[97]*rdvy2*rdy2; 
  alphaR[10] = 4.743416490252569*hamil[104]*rdvy2*rdy2+2.738612787525831*hamil[99]*rdvy2*rdy2; 
  alphaR[11] = 4.743416490252569*hamil[106]*rdvy2*rdy2+2.738612787525831*hamil[100]*rdvy2*rdy2; 
  alphaR[12] = 2.121320343559642*hamil[51]*rdvy2*rdy2+1.224744871391589*hamil[38]*rdvy2*rdy2; 
  alphaR[13] = 2.121320343559642*hamil[53]*rdvy2*rdy2+1.224744871391589*hamil[40]*rdvy2*rdy2; 
  alphaR[14] = 2.121320343559642*hamil[55]*rdvy2*rdy2+1.224744871391589*hamil[41]*rdvy2*rdy2; 
  alphaR[15] = 4.743416490252569*hamil[109]*rdvy2*rdy2+2.738612787525831*hamil[101]*rdvy2*rdy2; 
  alphaR[16] = 2.121320343559642*hamil[57]*rdvy2*rdy2+1.224744871391589*hamil[45]*rdvy2*rdy2; 
  alphaR[17] = 4.743416490252569*hamil[112]*rdvy2*rdy2+2.738612787525831*hamil[103]*rdvy2*rdy2; 
  alphaR[18] = 4.743416490252569*hamil[113]*rdvy2*rdy2+2.738612787525831*hamil[105]*rdvy2*rdy2; 
  alphaR[19] = 4.743416490252569*hamil[115]*rdvy2*rdy2+2.738612787525831*hamil[107]*rdvy2*rdy2; 
  alphaR[20] = 2.121320343559642*hamil[59]*rdvy2*rdy2+1.224744871391589*hamil[52]*rdvy2*rdy2; 
  alphaR[21] = 2.121320343559642*hamil[60]*rdvy2*rdy2+1.224744871391589*hamil[54]*rdvy2*rdy2; 
  alphaR[22] = 2.121320343559642*hamil[62]*rdvy2*rdy2+1.224744871391589*hamil[56]*rdvy2*rdy2; 
  alphaR[23] = 4.743416490252569*hamil[116]*rdvy2*rdy2+2.738612787525831*hamil[108]*rdvy2*rdy2; 
  alphaR[24] = 4.743416490252569*hamil[118]*rdvy2*rdy2+2.738612787525831*hamil[110]*rdvy2*rdy2; 
  alphaR[25] = 4.743416490252569*hamil[120]*rdvy2*rdy2+2.738612787525831*hamil[111]*rdvy2*rdy2; 
  alphaR[26] = 4.743416490252569*hamil[122]*rdvy2*rdy2+2.738612787525831*hamil[114]*rdvy2*rdy2; 
  alphaR[27] = 2.121320343559642*hamil[63]*rdvy2*rdy2+1.224744871391589*hamil[61]*rdvy2*rdy2; 
  alphaR[28] = 4.743416490252569*hamil[123]*rdvy2*rdy2+2.738612787525831*hamil[117]*rdvy2*rdy2; 
  alphaR[29] = 4.743416490252569*hamil[124]*rdvy2*rdy2+2.738612787525831*hamil[119]*rdvy2*rdy2; 
  alphaR[30] = 4.743416490252569*hamil[126]*rdvy2*rdy2+2.738612787525831*hamil[121]*rdvy2*rdy2; 
  alphaR[31] = 4.743416490252569*hamil[127]*rdvy2*rdy2+2.738612787525831*hamil[125]*rdvy2*rdy2; 

  int const_sgn_alpha_surf = 1;  
  
  if ((-0.1767766952966367*alphaR[31])+0.1767766952966367*(alphaR[30]+alphaR[29]+alphaR[28]+alphaR[27]+alphaR[26])-0.1767766952966367*(alphaR[25]+alphaR[24]+alphaR[23]+alphaR[22]+alphaR[21]+alphaR[20]+alphaR[19]+alphaR[18]+alphaR[17]+alphaR[16])+0.1767766952966367*(alphaR[15]+alphaR[14]+alphaR[13]+alphaR[12]+alphaR[11]+alphaR[10]+alphaR[9]+alphaR[8]+alphaR[7]+alphaR[6])-0.1767766952966367*(alphaR[5]+alphaR[4]+alphaR[3]+alphaR[2]+alphaR[1])+0.1767766952966367*alphaR[0] > 0.) 
    sgn_alpha_surfR[0] = 1.0; 
  else  
    sgn_alpha_surfR[0] = -1.0; 
  
  if (0.1767766952966367*alphaR[31]-0.1767766952966367*(alphaR[30]+alphaR[29]+alphaR[28]+alphaR[27])+0.1767766952966367*(alphaR[26]+alphaR[25]+alphaR[24]+alphaR[23]+alphaR[22]+alphaR[21]+alphaR[20])-0.1767766952966367*(alphaR[19]+alphaR[18]+alphaR[17]+alphaR[16]+alphaR[15]+alphaR[14]+alphaR[13]+alphaR[12])+0.1767766952966367*(alphaR[11]+alphaR[10]+alphaR[9]+alphaR[8]+alphaR[7]+alphaR[6]+alphaR[5])-0.1767766952966367*(alphaR[4]+alphaR[3]+alphaR[2]+alphaR[1])+0.1767766952966367*alphaR[0] > 0.) 
    sgn_alpha_surfR[1] = 1.0; 
  else  
    sgn_alpha_surfR[1] = -1.0; 
  
  if (sgn_alpha_surfR[1] == sgn_alpha_surfR[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*alphaR[31]-0.1767766952966367*(alphaR[30]+alphaR[29]+alphaR[28])+0.1767766952966367*alphaR[27]-0.1767766952966367*alphaR[26]+0.1767766952966367*(alphaR[25]+alphaR[24]+alphaR[23])-0.1767766952966367*(alphaR[22]+alphaR[21]+alphaR[20])+0.1767766952966367*(alphaR[19]+alphaR[18]+alphaR[17])-0.1767766952966367*(alphaR[16]+alphaR[15])+0.1767766952966367*(alphaR[14]+alphaR[13]+alphaR[12])-0.1767766952966367*(alphaR[11]+alphaR[10]+alphaR[9])+0.1767766952966367*(alphaR[8]+alphaR[7]+alphaR[6])-0.1767766952966367*alphaR[5]+0.1767766952966367*alphaR[4]-0.1767766952966367*(alphaR[3]+alphaR[2]+alphaR[1])+0.1767766952966367*alphaR[0] > 0.) 
    sgn_alpha_surfR[2] = 1.0; 
  else  
    sgn_alpha_surfR[2] = -1.0; 
  
  if (sgn_alpha_surfR[2] == sgn_alpha_surfR[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*alphaR[31])+0.1767766952966367*(alphaR[30]+alphaR[29]+alphaR[28])-0.1767766952966367*(alphaR[27]+alphaR[26]+alphaR[25]+alphaR[24]+alphaR[23])+0.1767766952966367*(alphaR[22]+alphaR[21]+alphaR[20]+alphaR[19]+alphaR[18]+alphaR[17])-0.1767766952966367*alphaR[16]+0.1767766952966367*alphaR[15]-0.1767766952966367*(alphaR[14]+alphaR[13]+alphaR[12]+alphaR[11]+alphaR[10]+alphaR[9])+0.1767766952966367*(alphaR[8]+alphaR[7]+alphaR[6]+alphaR[5]+alphaR[4])-0.1767766952966367*(alphaR[3]+alphaR[2]+alphaR[1])+0.1767766952966367*alphaR[0] > 0.) 
    sgn_alpha_surfR[3] = 1.0; 
  else  
    sgn_alpha_surfR[3] = -1.0; 
  
  if (sgn_alpha_surfR[3] == sgn_alpha_surfR[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*alphaR[31]-0.1767766952966367*(alphaR[30]+alphaR[29])+0.1767766952966367*alphaR[28]-0.1767766952966367*(alphaR[27]+alphaR[26])+0.1767766952966367*alphaR[25]-0.1767766952966367*(alphaR[24]+alphaR[23])+0.1767766952966367*(alphaR[22]+alphaR[21])-0.1767766952966367*alphaR[20]+0.1767766952966367*(alphaR[19]+alphaR[18])-0.1767766952966367*alphaR[17]+0.1767766952966367*(alphaR[16]+alphaR[15])-0.1767766952966367*alphaR[14]+0.1767766952966367*(alphaR[13]+alphaR[12])-0.1767766952966367*alphaR[11]+0.1767766952966367*(alphaR[10]+alphaR[9])-0.1767766952966367*(alphaR[8]+alphaR[7])+0.1767766952966367*alphaR[6]-0.1767766952966367*(alphaR[5]+alphaR[4])+0.1767766952966367*alphaR[3]-0.1767766952966367*(alphaR[2]+alphaR[1])+0.1767766952966367*alphaR[0] > 0.) 
    sgn_alpha_surfR[4] = 1.0; 
  else  
    sgn_alpha_surfR[4] = -1.0; 
  
  if (sgn_alpha_surfR[4] == sgn_alpha_surfR[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*alphaR[31])+0.1767766952966367*(alphaR[30]+alphaR[29])-0.1767766952966367*alphaR[28]+0.1767766952966367*alphaR[27]-0.1767766952966367*(alphaR[26]+alphaR[25])+0.1767766952966367*(alphaR[24]+alphaR[23])-0.1767766952966367*(alphaR[22]+alphaR[21])+0.1767766952966367*(alphaR[20]+alphaR[19]+alphaR[18])-0.1767766952966367*alphaR[17]+0.1767766952966367*alphaR[16]-0.1767766952966367*alphaR[15]+0.1767766952966367*alphaR[14]-0.1767766952966367*(alphaR[13]+alphaR[12]+alphaR[11])+0.1767766952966367*(alphaR[10]+alphaR[9])-0.1767766952966367*(alphaR[8]+alphaR[7])+0.1767766952966367*(alphaR[6]+alphaR[5])-0.1767766952966367*alphaR[4]+0.1767766952966367*alphaR[3]-0.1767766952966367*(alphaR[2]+alphaR[1])+0.1767766952966367*alphaR[0] > 0.) 
    sgn_alpha_surfR[5] = 1.0; 
  else  
    sgn_alpha_surfR[5] = -1.0; 
  
  if (sgn_alpha_surfR[5] == sgn_alpha_surfR[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*alphaR[31])+0.1767766952966367*(alphaR[30]+alphaR[29])-0.1767766952966367*(alphaR[28]+alphaR[27])+0.1767766952966367*alphaR[26]-0.1767766952966367*alphaR[25]+0.1767766952966367*(alphaR[24]+alphaR[23]+alphaR[22]+alphaR[21])-0.1767766952966367*(alphaR[20]+alphaR[19]+alphaR[18])+0.1767766952966367*(alphaR[17]+alphaR[16])-0.1767766952966367*(alphaR[15]+alphaR[14])+0.1767766952966367*(alphaR[13]+alphaR[12]+alphaR[11])-0.1767766952966367*(alphaR[10]+alphaR[9]+alphaR[8]+alphaR[7])+0.1767766952966367*alphaR[6]-0.1767766952966367*alphaR[5]+0.1767766952966367*(alphaR[4]+alphaR[3])-0.1767766952966367*(alphaR[2]+alphaR[1])+0.1767766952966367*alphaR[0] > 0.) 
    sgn_alpha_surfR[6] = 1.0; 
  else  
    sgn_alpha_surfR[6] = -1.0; 
  
  if (sgn_alpha_surfR[6] == sgn_alpha_surfR[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*alphaR[31]-0.1767766952966367*(alphaR[30]+alphaR[29])+0.1767766952966367*(alphaR[28]+alphaR[27]+alphaR[26]+alphaR[25])-0.1767766952966367*(alphaR[24]+alphaR[23]+alphaR[22]+alphaR[21])+0.1767766952966367*alphaR[20]-0.1767766952966367*(alphaR[19]+alphaR[18])+0.1767766952966367*(alphaR[17]+alphaR[16]+alphaR[15]+alphaR[14])-0.1767766952966367*(alphaR[13]+alphaR[12])+0.1767766952966367*alphaR[11]-0.1767766952966367*(alphaR[10]+alphaR[9]+alphaR[8]+alphaR[7])+0.1767766952966367*(alphaR[6]+alphaR[5]+alphaR[4]+alphaR[3])-0.1767766952966367*(alphaR[2]+alphaR[1])+0.1767766952966367*alphaR[0] > 0.) 
    sgn_alpha_surfR[7] = 1.0; 
  else  
    sgn_alpha_surfR[7] = -1.0; 
  
  if (sgn_alpha_surfR[7] == sgn_alpha_surfR[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*alphaR[31]-0.1767766952966367*alphaR[30]+0.1767766952966367*alphaR[29]-0.1767766952966367*(alphaR[28]+alphaR[27]+alphaR[26]+alphaR[25])+0.1767766952966367*alphaR[24]-0.1767766952966367*alphaR[23]+0.1767766952966367*alphaR[22]-0.1767766952966367*alphaR[21]+0.1767766952966367*(alphaR[20]+alphaR[19])-0.1767766952966367*alphaR[18]+0.1767766952966367*(alphaR[17]+alphaR[16]+alphaR[15]+alphaR[14])-0.1767766952966367*alphaR[13]+0.1767766952966367*(alphaR[12]+alphaR[11])-0.1767766952966367*alphaR[10]+0.1767766952966367*alphaR[9]-0.1767766952966367*alphaR[8]+0.1767766952966367*alphaR[7]-0.1767766952966367*(alphaR[6]+alphaR[5]+alphaR[4]+alphaR[3])+0.1767766952966367*alphaR[2]-0.1767766952966367*alphaR[1]+0.1767766952966367*alphaR[0] > 0.) 
    sgn_alpha_surfR[8] = 1.0; 
  else  
    sgn_alpha_surfR[8] = -1.0; 
  
  if (sgn_alpha_surfR[8] == sgn_alpha_surfR[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*alphaR[31])+0.1767766952966367*alphaR[30]-0.1767766952966367*alphaR[29]+0.1767766952966367*(alphaR[28]+alphaR[27])-0.1767766952966367*alphaR[26]+0.1767766952966367*alphaR[25]-0.1767766952966367*alphaR[24]+0.1767766952966367*alphaR[23]-0.1767766952966367*alphaR[22]+0.1767766952966367*alphaR[21]-0.1767766952966367*alphaR[20]+0.1767766952966367*alphaR[19]-0.1767766952966367*alphaR[18]+0.1767766952966367*(alphaR[17]+alphaR[16])-0.1767766952966367*(alphaR[15]+alphaR[14])+0.1767766952966367*alphaR[13]-0.1767766952966367*alphaR[12]+0.1767766952966367*alphaR[11]-0.1767766952966367*alphaR[10]+0.1767766952966367*alphaR[9]-0.1767766952966367*alphaR[8]+0.1767766952966367*alphaR[7]-0.1767766952966367*alphaR[6]+0.1767766952966367*alphaR[5]-0.1767766952966367*(alphaR[4]+alphaR[3])+0.1767766952966367*alphaR[2]-0.1767766952966367*alphaR[1]+0.1767766952966367*alphaR[0] > 0.) 
    sgn_alpha_surfR[9] = 1.0; 
  else  
    sgn_alpha_surfR[9] = -1.0; 
  
  if (sgn_alpha_surfR[9] == sgn_alpha_surfR[8]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*alphaR[31])+0.1767766952966367*alphaR[30]-0.1767766952966367*alphaR[29]+0.1767766952966367*alphaR[28]-0.1767766952966367*alphaR[27]+0.1767766952966367*(alphaR[26]+alphaR[25])-0.1767766952966367*alphaR[24]+0.1767766952966367*(alphaR[23]+alphaR[22])-0.1767766952966367*alphaR[21]+0.1767766952966367*alphaR[20]-0.1767766952966367*alphaR[19]+0.1767766952966367*alphaR[18]-0.1767766952966367*alphaR[17]+0.1767766952966367*alphaR[16]-0.1767766952966367*alphaR[15]+0.1767766952966367*alphaR[14]-0.1767766952966367*alphaR[13]+0.1767766952966367*alphaR[12]-0.1767766952966367*alphaR[11]+0.1767766952966367*alphaR[10]-0.1767766952966367*(alphaR[9]+alphaR[8])+0.1767766952966367*alphaR[7]-0.1767766952966367*(alphaR[6]+alphaR[5])+0.1767766952966367*alphaR[4]-0.1767766952966367*alphaR[3]+0.1767766952966367*alphaR[2]-0.1767766952966367*alphaR[1]+0.1767766952966367*alphaR[0] > 0.) 
    sgn_alpha_surfR[10] = 1.0; 
  else  
    sgn_alpha_surfR[10] = -1.0; 
  
  if (sgn_alpha_surfR[10] == sgn_alpha_surfR[9]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*alphaR[31]-0.1767766952966367*alphaR[30]+0.1767766952966367*alphaR[29]-0.1767766952966367*alphaR[28]+0.1767766952966367*(alphaR[27]+alphaR[26])-0.1767766952966367*alphaR[25]+0.1767766952966367*alphaR[24]-0.1767766952966367*(alphaR[23]+alphaR[22])+0.1767766952966367*alphaR[21]-0.1767766952966367*(alphaR[20]+alphaR[19])+0.1767766952966367*alphaR[18]-0.1767766952966367*alphaR[17]+0.1767766952966367*(alphaR[16]+alphaR[15])-0.1767766952966367*alphaR[14]+0.1767766952966367*alphaR[13]-0.1767766952966367*(alphaR[12]+alphaR[11])+0.1767766952966367*alphaR[10]-0.1767766952966367*(alphaR[9]+alphaR[8])+0.1767766952966367*alphaR[7]-0.1767766952966367*alphaR[6]+0.1767766952966367*(alphaR[5]+alphaR[4])-0.1767766952966367*alphaR[3]+0.1767766952966367*alphaR[2]-0.1767766952966367*alphaR[1]+0.1767766952966367*alphaR[0] > 0.) 
    sgn_alpha_surfR[11] = 1.0; 
  else  
    sgn_alpha_surfR[11] = -1.0; 
  
  if (sgn_alpha_surfR[11] == sgn_alpha_surfR[10]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*alphaR[31])+0.1767766952966367*alphaR[30]-0.1767766952966367*(alphaR[29]+alphaR[28])+0.1767766952966367*(alphaR[27]+alphaR[26]+alphaR[25]+alphaR[24])-0.1767766952966367*(alphaR[23]+alphaR[22])+0.1767766952966367*(alphaR[21]+alphaR[20])-0.1767766952966367*alphaR[19]+0.1767766952966367*(alphaR[18]+alphaR[17])-0.1767766952966367*alphaR[16]+0.1767766952966367*alphaR[15]-0.1767766952966367*(alphaR[14]+alphaR[13])+0.1767766952966367*alphaR[12]-0.1767766952966367*(alphaR[11]+alphaR[10])+0.1767766952966367*(alphaR[9]+alphaR[8])-0.1767766952966367*(alphaR[7]+alphaR[6]+alphaR[5]+alphaR[4])+0.1767766952966367*(alphaR[3]+alphaR[2])-0.1767766952966367*alphaR[1]+0.1767766952966367*alphaR[0] > 0.) 
    sgn_alpha_surfR[12] = 1.0; 
  else  
    sgn_alpha_surfR[12] = -1.0; 
  
  if (sgn_alpha_surfR[12] == sgn_alpha_surfR[11]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*alphaR[31]-0.1767766952966367*alphaR[30]+0.1767766952966367*(alphaR[29]+alphaR[28])-0.1767766952966367*alphaR[27]+0.1767766952966367*alphaR[26]-0.1767766952966367*(alphaR[25]+alphaR[24])+0.1767766952966367*(alphaR[23]+alphaR[22])-0.1767766952966367*(alphaR[21]+alphaR[20]+alphaR[19])+0.1767766952966367*(alphaR[18]+alphaR[17])-0.1767766952966367*(alphaR[16]+alphaR[15])+0.1767766952966367*(alphaR[14]+alphaR[13])-0.1767766952966367*(alphaR[12]+alphaR[11]+alphaR[10])+0.1767766952966367*(alphaR[9]+alphaR[8])-0.1767766952966367*(alphaR[7]+alphaR[6])+0.1767766952966367*alphaR[5]-0.1767766952966367*alphaR[4]+0.1767766952966367*(alphaR[3]+alphaR[2])-0.1767766952966367*alphaR[1]+0.1767766952966367*alphaR[0] > 0.) 
    sgn_alpha_surfR[13] = 1.0; 
  else  
    sgn_alpha_surfR[13] = -1.0; 
  
  if (sgn_alpha_surfR[13] == sgn_alpha_surfR[12]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*alphaR[31]-0.1767766952966367*alphaR[30]+0.1767766952966367*(alphaR[29]+alphaR[28]+alphaR[27])-0.1767766952966367*(alphaR[26]+alphaR[25]+alphaR[24])+0.1767766952966367*alphaR[23]-0.1767766952966367*alphaR[22]+0.1767766952966367*(alphaR[21]+alphaR[20]+alphaR[19])-0.1767766952966367*(alphaR[18]+alphaR[17]+alphaR[16]+alphaR[15]+alphaR[14]+alphaR[13])+0.1767766952966367*(alphaR[12]+alphaR[11]+alphaR[10])-0.1767766952966367*alphaR[9]+0.1767766952966367*alphaR[8]-0.1767766952966367*(alphaR[7]+alphaR[6]+alphaR[5])+0.1767766952966367*(alphaR[4]+alphaR[3]+alphaR[2])-0.1767766952966367*alphaR[1]+0.1767766952966367*alphaR[0] > 0.) 
    sgn_alpha_surfR[14] = 1.0; 
  else  
    sgn_alpha_surfR[14] = -1.0; 
  
  if (sgn_alpha_surfR[14] == sgn_alpha_surfR[13]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*alphaR[31])+0.1767766952966367*alphaR[30]-0.1767766952966367*(alphaR[29]+alphaR[28]+alphaR[27]+alphaR[26])+0.1767766952966367*(alphaR[25]+alphaR[24])-0.1767766952966367*alphaR[23]+0.1767766952966367*alphaR[22]-0.1767766952966367*(alphaR[21]+alphaR[20])+0.1767766952966367*alphaR[19]-0.1767766952966367*(alphaR[18]+alphaR[17]+alphaR[16])+0.1767766952966367*(alphaR[15]+alphaR[14]+alphaR[13])-0.1767766952966367*alphaR[12]+0.1767766952966367*(alphaR[11]+alphaR[10])-0.1767766952966367*alphaR[9]+0.1767766952966367*alphaR[8]-0.1767766952966367*(alphaR[7]+alphaR[6])+0.1767766952966367*(alphaR[5]+alphaR[4]+alphaR[3]+alphaR[2])-0.1767766952966367*alphaR[1]+0.1767766952966367*alphaR[0] > 0.) 
    sgn_alpha_surfR[15] = 1.0; 
  else  
    sgn_alpha_surfR[15] = -1.0; 
  
  if (sgn_alpha_surfR[15] == sgn_alpha_surfR[14]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*(alphaR[31]+alphaR[30])-0.1767766952966367*(alphaR[29]+alphaR[28]+alphaR[27]+alphaR[26]+alphaR[25]+alphaR[24])+0.1767766952966367*alphaR[23]-0.1767766952966367*alphaR[22]+0.1767766952966367*(alphaR[21]+alphaR[20])-0.1767766952966367*alphaR[19]+0.1767766952966367*(alphaR[18]+alphaR[17]+alphaR[16]+alphaR[15]+alphaR[14]+alphaR[13])-0.1767766952966367*alphaR[12]+0.1767766952966367*(alphaR[11]+alphaR[10])-0.1767766952966367*alphaR[9]+0.1767766952966367*alphaR[8]-0.1767766952966367*(alphaR[7]+alphaR[6]+alphaR[5]+alphaR[4]+alphaR[3]+alphaR[2])+0.1767766952966367*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[16] = 1.0; 
  else  
    sgn_alpha_surfR[16] = -1.0; 
  
  if (sgn_alpha_surfR[16] == sgn_alpha_surfR[15]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*(alphaR[31]+alphaR[30]))+0.1767766952966367*(alphaR[29]+alphaR[28]+alphaR[27])-0.1767766952966367*alphaR[26]+0.1767766952966367*(alphaR[25]+alphaR[24])-0.1767766952966367*alphaR[23]+0.1767766952966367*alphaR[22]-0.1767766952966367*(alphaR[21]+alphaR[20]+alphaR[19])+0.1767766952966367*(alphaR[18]+alphaR[17]+alphaR[16])-0.1767766952966367*(alphaR[15]+alphaR[14]+alphaR[13])+0.1767766952966367*(alphaR[12]+alphaR[11]+alphaR[10])-0.1767766952966367*alphaR[9]+0.1767766952966367*alphaR[8]-0.1767766952966367*(alphaR[7]+alphaR[6])+0.1767766952966367*alphaR[5]-0.1767766952966367*(alphaR[4]+alphaR[3]+alphaR[2])+0.1767766952966367*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[17] = 1.0; 
  else  
    sgn_alpha_surfR[17] = -1.0; 
  
  if (sgn_alpha_surfR[17] == sgn_alpha_surfR[16]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*(alphaR[31]+alphaR[30]))+0.1767766952966367*(alphaR[29]+alphaR[28])-0.1767766952966367*alphaR[27]+0.1767766952966367*(alphaR[26]+alphaR[25]+alphaR[24])-0.1767766952966367*(alphaR[23]+alphaR[22])+0.1767766952966367*(alphaR[21]+alphaR[20]+alphaR[19])-0.1767766952966367*(alphaR[18]+alphaR[17])+0.1767766952966367*alphaR[16]-0.1767766952966367*alphaR[15]+0.1767766952966367*(alphaR[14]+alphaR[13])-0.1767766952966367*(alphaR[12]+alphaR[11]+alphaR[10])+0.1767766952966367*(alphaR[9]+alphaR[8])-0.1767766952966367*(alphaR[7]+alphaR[6]+alphaR[5])+0.1767766952966367*alphaR[4]-0.1767766952966367*(alphaR[3]+alphaR[2])+0.1767766952966367*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[18] = 1.0; 
  else  
    sgn_alpha_surfR[18] = -1.0; 
  
  if (sgn_alpha_surfR[18] == sgn_alpha_surfR[17]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*(alphaR[31]+alphaR[30])-0.1767766952966367*(alphaR[29]+alphaR[28])+0.1767766952966367*(alphaR[27]+alphaR[26])-0.1767766952966367*(alphaR[25]+alphaR[24])+0.1767766952966367*(alphaR[23]+alphaR[22])-0.1767766952966367*(alphaR[21]+alphaR[20])+0.1767766952966367*alphaR[19]-0.1767766952966367*(alphaR[18]+alphaR[17])+0.1767766952966367*(alphaR[16]+alphaR[15])-0.1767766952966367*(alphaR[14]+alphaR[13])+0.1767766952966367*alphaR[12]-0.1767766952966367*(alphaR[11]+alphaR[10])+0.1767766952966367*(alphaR[9]+alphaR[8])-0.1767766952966367*(alphaR[7]+alphaR[6])+0.1767766952966367*(alphaR[5]+alphaR[4])-0.1767766952966367*(alphaR[3]+alphaR[2])+0.1767766952966367*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[19] = 1.0; 
  else  
    sgn_alpha_surfR[19] = -1.0; 
  
  if (sgn_alpha_surfR[19] == sgn_alpha_surfR[18]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*(alphaR[31]+alphaR[30]))+0.1767766952966367*alphaR[29]-0.1767766952966367*alphaR[28]+0.1767766952966367*(alphaR[27]+alphaR[26]+alphaR[25])-0.1767766952966367*alphaR[24]+0.1767766952966367*(alphaR[23]+alphaR[22])-0.1767766952966367*alphaR[21]+0.1767766952966367*(alphaR[20]+alphaR[19])-0.1767766952966367*alphaR[18]+0.1767766952966367*alphaR[17]-0.1767766952966367*alphaR[16]+0.1767766952966367*alphaR[15]-0.1767766952966367*alphaR[14]+0.1767766952966367*alphaR[13]-0.1767766952966367*(alphaR[12]+alphaR[11])+0.1767766952966367*alphaR[10]-0.1767766952966367*(alphaR[9]+alphaR[8])+0.1767766952966367*alphaR[7]-0.1767766952966367*(alphaR[6]+alphaR[5]+alphaR[4])+0.1767766952966367*alphaR[3]-0.1767766952966367*alphaR[2]+0.1767766952966367*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[20] = 1.0; 
  else  
    sgn_alpha_surfR[20] = -1.0; 
  
  if (sgn_alpha_surfR[20] == sgn_alpha_surfR[19]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*(alphaR[31]+alphaR[30])-0.1767766952966367*alphaR[29]+0.1767766952966367*alphaR[28]-0.1767766952966367*alphaR[27]+0.1767766952966367*alphaR[26]-0.1767766952966367*alphaR[25]+0.1767766952966367*alphaR[24]-0.1767766952966367*(alphaR[23]+alphaR[22])+0.1767766952966367*alphaR[21]-0.1767766952966367*alphaR[20]+0.1767766952966367*alphaR[19]-0.1767766952966367*alphaR[18]+0.1767766952966367*alphaR[17]-0.1767766952966367*(alphaR[16]+alphaR[15])+0.1767766952966367*alphaR[14]-0.1767766952966367*alphaR[13]+0.1767766952966367*alphaR[12]-0.1767766952966367*alphaR[11]+0.1767766952966367*alphaR[10]-0.1767766952966367*(alphaR[9]+alphaR[8])+0.1767766952966367*alphaR[7]-0.1767766952966367*alphaR[6]+0.1767766952966367*alphaR[5]-0.1767766952966367*alphaR[4]+0.1767766952966367*alphaR[3]-0.1767766952966367*alphaR[2]+0.1767766952966367*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[21] = 1.0; 
  else  
    sgn_alpha_surfR[21] = -1.0; 
  
  if (sgn_alpha_surfR[21] == sgn_alpha_surfR[20]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*(alphaR[31]+alphaR[30])-0.1767766952966367*alphaR[29]+0.1767766952966367*(alphaR[28]+alphaR[27])-0.1767766952966367*(alphaR[26]+alphaR[25])+0.1767766952966367*alphaR[24]-0.1767766952966367*alphaR[23]+0.1767766952966367*alphaR[22]-0.1767766952966367*alphaR[21]+0.1767766952966367*alphaR[20]-0.1767766952966367*alphaR[19]+0.1767766952966367*alphaR[18]-0.1767766952966367*(alphaR[17]+alphaR[16]+alphaR[15]+alphaR[14])+0.1767766952966367*alphaR[13]-0.1767766952966367*alphaR[12]+0.1767766952966367*alphaR[11]-0.1767766952966367*alphaR[10]+0.1767766952966367*alphaR[9]-0.1767766952966367*alphaR[8]+0.1767766952966367*alphaR[7]-0.1767766952966367*(alphaR[6]+alphaR[5])+0.1767766952966367*(alphaR[4]+alphaR[3])-0.1767766952966367*alphaR[2]+0.1767766952966367*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[22] = 1.0; 
  else  
    sgn_alpha_surfR[22] = -1.0; 
  
  if (sgn_alpha_surfR[22] == sgn_alpha_surfR[21]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*(alphaR[31]+alphaR[30]))+0.1767766952966367*alphaR[29]-0.1767766952966367*(alphaR[28]+alphaR[27]+alphaR[26])+0.1767766952966367*alphaR[25]-0.1767766952966367*alphaR[24]+0.1767766952966367*alphaR[23]-0.1767766952966367*alphaR[22]+0.1767766952966367*alphaR[21]-0.1767766952966367*(alphaR[20]+alphaR[19])+0.1767766952966367*alphaR[18]-0.1767766952966367*(alphaR[17]+alphaR[16])+0.1767766952966367*(alphaR[15]+alphaR[14])-0.1767766952966367*alphaR[13]+0.1767766952966367*(alphaR[12]+alphaR[11])-0.1767766952966367*alphaR[10]+0.1767766952966367*alphaR[9]-0.1767766952966367*alphaR[8]+0.1767766952966367*alphaR[7]-0.1767766952966367*alphaR[6]+0.1767766952966367*(alphaR[5]+alphaR[4]+alphaR[3])-0.1767766952966367*alphaR[2]+0.1767766952966367*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[23] = 1.0; 
  else  
    sgn_alpha_surfR[23] = -1.0; 
  
  if (sgn_alpha_surfR[23] == sgn_alpha_surfR[22]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*(alphaR[31]+alphaR[30]+alphaR[29]))+0.1767766952966367*(alphaR[28]+alphaR[27]+alphaR[26])-0.1767766952966367*alphaR[25]+0.1767766952966367*(alphaR[24]+alphaR[23]+alphaR[22]+alphaR[21])-0.1767766952966367*alphaR[20]+0.1767766952966367*(alphaR[19]+alphaR[18])-0.1767766952966367*(alphaR[17]+alphaR[16])+0.1767766952966367*(alphaR[15]+alphaR[14])-0.1767766952966367*(alphaR[13]+alphaR[12])+0.1767766952966367*alphaR[11]-0.1767766952966367*(alphaR[10]+alphaR[9]+alphaR[8]+alphaR[7])+0.1767766952966367*alphaR[6]-0.1767766952966367*(alphaR[5]+alphaR[4]+alphaR[3])+0.1767766952966367*(alphaR[2]+alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[24] = 1.0; 
  else  
    sgn_alpha_surfR[24] = -1.0; 
  
  if (sgn_alpha_surfR[24] == sgn_alpha_surfR[23]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*(alphaR[31]+alphaR[30]+alphaR[29])-0.1767766952966367*(alphaR[28]+alphaR[27])+0.1767766952966367*(alphaR[26]+alphaR[25])-0.1767766952966367*(alphaR[24]+alphaR[23]+alphaR[22]+alphaR[21])+0.1767766952966367*(alphaR[20]+alphaR[19]+alphaR[18])-0.1767766952966367*(alphaR[17]+alphaR[16]+alphaR[15]+alphaR[14])+0.1767766952966367*(alphaR[13]+alphaR[12]+alphaR[11])-0.1767766952966367*(alphaR[10]+alphaR[9]+alphaR[8]+alphaR[7])+0.1767766952966367*(alphaR[6]+alphaR[5])-0.1767766952966367*(alphaR[4]+alphaR[3])+0.1767766952966367*(alphaR[2]+alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[25] = 1.0; 
  else  
    sgn_alpha_surfR[25] = -1.0; 
  
  if (sgn_alpha_surfR[25] == sgn_alpha_surfR[24]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*(alphaR[31]+alphaR[30]+alphaR[29])-0.1767766952966367*alphaR[28]+0.1767766952966367*alphaR[27]-0.1767766952966367*alphaR[26]+0.1767766952966367*alphaR[25]-0.1767766952966367*(alphaR[24]+alphaR[23])+0.1767766952966367*(alphaR[22]+alphaR[21])-0.1767766952966367*(alphaR[20]+alphaR[19]+alphaR[18])+0.1767766952966367*alphaR[17]-0.1767766952966367*(alphaR[16]+alphaR[15])+0.1767766952966367*alphaR[14]-0.1767766952966367*(alphaR[13]+alphaR[12]+alphaR[11])+0.1767766952966367*(alphaR[10]+alphaR[9])-0.1767766952966367*(alphaR[8]+alphaR[7])+0.1767766952966367*alphaR[6]-0.1767766952966367*alphaR[5]+0.1767766952966367*alphaR[4]-0.1767766952966367*alphaR[3]+0.1767766952966367*(alphaR[2]+alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[26] = 1.0; 
  else  
    sgn_alpha_surfR[26] = -1.0; 
  
  if (sgn_alpha_surfR[26] == sgn_alpha_surfR[25]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*(alphaR[31]+alphaR[30]+alphaR[29]))+0.1767766952966367*alphaR[28]-0.1767766952966367*(alphaR[27]+alphaR[26]+alphaR[25])+0.1767766952966367*(alphaR[24]+alphaR[23])-0.1767766952966367*(alphaR[22]+alphaR[21])+0.1767766952966367*alphaR[20]-0.1767766952966367*(alphaR[19]+alphaR[18])+0.1767766952966367*alphaR[17]-0.1767766952966367*alphaR[16]+0.1767766952966367*alphaR[15]-0.1767766952966367*alphaR[14]+0.1767766952966367*(alphaR[13]+alphaR[12])-0.1767766952966367*alphaR[11]+0.1767766952966367*(alphaR[10]+alphaR[9])-0.1767766952966367*(alphaR[8]+alphaR[7])+0.1767766952966367*(alphaR[6]+alphaR[5]+alphaR[4])-0.1767766952966367*alphaR[3]+0.1767766952966367*(alphaR[2]+alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[27] = 1.0; 
  else  
    sgn_alpha_surfR[27] = -1.0; 
  
  if (sgn_alpha_surfR[27] == sgn_alpha_surfR[26]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*(alphaR[31]+alphaR[30]+alphaR[29]+alphaR[28])-0.1767766952966367*(alphaR[27]+alphaR[26])+0.1767766952966367*(alphaR[25]+alphaR[24]+alphaR[23])-0.1767766952966367*(alphaR[22]+alphaR[21]+alphaR[20]+alphaR[19]+alphaR[18]+alphaR[17])+0.1767766952966367*(alphaR[16]+alphaR[15])-0.1767766952966367*(alphaR[14]+alphaR[13]+alphaR[12]+alphaR[11]+alphaR[10]+alphaR[9])+0.1767766952966367*(alphaR[8]+alphaR[7]+alphaR[6])-0.1767766952966367*(alphaR[5]+alphaR[4])+0.1767766952966367*(alphaR[3]+alphaR[2]+alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[28] = 1.0; 
  else  
    sgn_alpha_surfR[28] = -1.0; 
  
  if (sgn_alpha_surfR[28] == sgn_alpha_surfR[27]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*(alphaR[31]+alphaR[30]+alphaR[29]+alphaR[28]))+0.1767766952966367*alphaR[27]-0.1767766952966367*(alphaR[26]+alphaR[25]+alphaR[24]+alphaR[23])+0.1767766952966367*(alphaR[22]+alphaR[21]+alphaR[20])-0.1767766952966367*(alphaR[19]+alphaR[18]+alphaR[17])+0.1767766952966367*alphaR[16]-0.1767766952966367*alphaR[15]+0.1767766952966367*(alphaR[14]+alphaR[13]+alphaR[12])-0.1767766952966367*(alphaR[11]+alphaR[10]+alphaR[9])+0.1767766952966367*(alphaR[8]+alphaR[7]+alphaR[6]+alphaR[5])-0.1767766952966367*alphaR[4]+0.1767766952966367*(alphaR[3]+alphaR[2]+alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[29] = 1.0; 
  else  
    sgn_alpha_surfR[29] = -1.0; 
  
  if (sgn_alpha_surfR[29] == sgn_alpha_surfR[28]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*(alphaR[31]+alphaR[30]+alphaR[29]+alphaR[28]+alphaR[27]))+0.1767766952966367*alphaR[26]-0.1767766952966367*(alphaR[25]+alphaR[24]+alphaR[23]+alphaR[22]+alphaR[21]+alphaR[20])+0.1767766952966367*(alphaR[19]+alphaR[18]+alphaR[17]+alphaR[16])-0.1767766952966367*(alphaR[15]+alphaR[14]+alphaR[13]+alphaR[12])+0.1767766952966367*(alphaR[11]+alphaR[10]+alphaR[9]+alphaR[8]+alphaR[7]+alphaR[6])-0.1767766952966367*alphaR[5]+0.1767766952966367*(alphaR[4]+alphaR[3]+alphaR[2]+alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[30] = 1.0; 
  else  
    sgn_alpha_surfR[30] = -1.0; 
  
  if (sgn_alpha_surfR[30] == sgn_alpha_surfR[29]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*(alphaR[31]+alphaR[30]+alphaR[29]+alphaR[28]+alphaR[27]+alphaR[26]+alphaR[25]+alphaR[24]+alphaR[23]+alphaR[22]+alphaR[21]+alphaR[20]+alphaR[19]+alphaR[18]+alphaR[17]+alphaR[16]+alphaR[15]+alphaR[14]+alphaR[13]+alphaR[12]+alphaR[11]+alphaR[10]+alphaR[9]+alphaR[8]+alphaR[7]+alphaR[6]+alphaR[5]+alphaR[4]+alphaR[3]+alphaR[2]+alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[31] = 1.0; 
  else  
    sgn_alpha_surfR[31] = -1.0; 
  
  if (sgn_alpha_surfR[31] == sgn_alpha_surfR[30]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
