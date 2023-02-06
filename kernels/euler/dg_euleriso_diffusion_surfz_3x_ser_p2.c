#include <dg_euleriso_diffusion_kernels.h>

GKYL_CU_DH void
dg_euleriso_diffusion_surfz_3x_ser_p2(const double* w, const double* dx,
  const double* D_in,
  const double* uvarl, const double* uvarc, const double* uvarr,
  const double* statevecl, const double* statevecc, const double* statevecr,
  double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // D: Diffusion coefficient in the center cell
  // uvarl: Input velocity in the left cell
  // uvarc: Input velocity in the center cell
  // uvarr: Input velocity in the right cell
  // statevecl: Input field in the left cell
  // statevecc: Input field in the center cell
  // statevecr: Input field in the right cell
  // out: Incremented output

  const double J = 4/dx[2]/dx[2];
  const double *D = &D_in[20]; 
  double mu = D[0]; 
  const double *uvarxl = &uvarl[0]; 
  const double *uvaryl = &uvarl[20]; 
  const double *uvarzl = &uvarl[40]; 
  const double *uvarxc = &uvarc[0]; 
  const double *uvaryc = &uvarc[20]; 
  const double *uvarzc = &uvarc[40]; 
  const double *uvarxr = &uvarr[0]; 
  const double *uvaryr = &uvarr[20]; 
  const double *uvarzr = &uvarr[40]; 
  out[20] += J*(0.6708203932499369*uvarzr[9]*mu+0.6708203932499369*uvarzl[9]*mu-1.341640786499874*uvarzc[9]*mu+0.6708203932499369*uvaryr[8]*mu+0.6708203932499369*uvaryl[8]*mu-1.341640786499874*uvaryc[8]*mu+0.6708203932499369*uvarxr[7]*mu+0.6708203932499369*uvarxl[7]*mu-1.341640786499874*uvarxc[7]*mu-1.190784930203603*uvarzr[3]*mu+1.190784930203603*uvarzl[3]*mu-1.190784930203603*uvaryr[2]*mu+1.190784930203603*uvaryl[2]*mu-1.190784930203603*uvarxr[1]*mu+1.190784930203603*uvarxl[1]*mu+0.9375*uvarzr[0]*mu+0.9375*uvarzl[0]*mu-1.875*uvarzc[0]*mu+0.9375*uvaryr[0]*mu+0.9375*uvaryl[0]*mu-1.875*uvaryc[0]*mu+0.9375*uvarxr[0]*mu+0.9375*uvarxl[0]*mu-1.875*uvarxc[0]*mu);
  out[21] += J*(0.6708203932499369*uvarzr[15]*mu+0.6708203932499369*uvarzl[15]*mu-1.341640786499874*uvarzc[15]*mu+0.6708203932499369*uvaryr[12]*mu+0.6708203932499369*uvaryl[12]*mu-1.341640786499874*uvaryc[12]*mu+0.7382874503707888*uvarxr[7]*mu-0.7382874503707888*uvarxl[7]*mu-1.190784930203603*uvarzr[5]*mu+1.190784930203603*uvarzl[5]*mu-1.190784930203603*uvaryr[4]*mu+1.190784930203603*uvaryl[4]*mu+0.9375*uvarzr[1]*mu+0.9375*uvarzl[1]*mu-1.875*uvarzc[1]*mu+0.9375*uvaryr[1]*mu+0.9375*uvaryl[1]*mu-1.875*uvaryc[1]*mu-1.453125*uvarxr[1]*mu-1.453125*uvarxl[1]*mu-5.34375*uvarxc[1]*mu+1.190784930203603*uvarxr[0]*mu-1.190784930203603*uvarxl[0]*mu);
  out[22] += J*(0.6708203932499369*uvarzr[16]*mu+0.6708203932499369*uvarzl[16]*mu-1.341640786499874*uvarzc[16]*mu+0.6708203932499369*uvarxr[11]*mu+0.6708203932499369*uvarxl[11]*mu-1.341640786499874*uvarxc[11]*mu+0.7382874503707888*uvaryr[8]*mu-0.7382874503707888*uvaryl[8]*mu-1.190784930203603*uvarzr[6]*mu+1.190784930203603*uvarzl[6]*mu-1.190784930203603*uvarxr[4]*mu+1.190784930203603*uvarxl[4]*mu+0.9375*uvarzr[2]*mu+0.9375*uvarzl[2]*mu-1.875*uvarzc[2]*mu-1.453125*uvaryr[2]*mu-1.453125*uvaryl[2]*mu-5.34375*uvaryc[2]*mu+0.9375*uvarxr[2]*mu+0.9375*uvarxl[2]*mu-1.875*uvarxc[2]*mu+1.190784930203603*uvaryr[0]*mu-1.190784930203603*uvaryl[0]*mu);
  out[23] += J*(0.6708203932499369*uvaryr[14]*mu+0.6708203932499369*uvaryl[14]*mu-1.341640786499874*uvaryc[14]*mu+0.6708203932499369*uvarxr[13]*mu+0.6708203932499369*uvarxl[13]*mu-1.341640786499874*uvarxc[13]*mu+0.7382874503707888*uvarzr[9]*mu-0.7382874503707888*uvarzl[9]*mu-1.190784930203603*uvaryr[6]*mu+1.190784930203603*uvaryl[6]*mu-1.190784930203603*uvarxr[5]*mu+1.190784930203603*uvarxl[5]*mu-1.453125*uvarzr[3]*mu-1.453125*uvarzl[3]*mu-5.34375*uvarzc[3]*mu+0.9375*uvaryr[3]*mu+0.9375*uvaryl[3]*mu-1.875*uvaryc[3]*mu+0.9375*uvarxr[3]*mu+0.9375*uvarxl[3]*mu-1.875*uvarxc[3]*mu+1.190784930203603*uvarzr[0]*mu-1.190784930203603*uvarzl[0]*mu);
  out[24] += J*(0.6708203932499369*uvarzr[19]*mu+0.6708203932499369*uvarzl[19]*mu-1.341640786499874*uvarzc[19]*mu+0.7382874503707888*uvaryr[12]*mu-0.7382874503707888*uvaryl[12]*mu+0.7382874503707888*uvarxr[11]*mu-0.7382874503707888*uvarxl[11]*mu-1.190784930203603*uvarzr[10]*mu+1.190784930203603*uvarzl[10]*mu+0.9375*uvarzr[4]*mu+0.9375*uvarzl[4]*mu-1.875*uvarzc[4]*mu-1.453125*uvaryr[4]*mu-1.453125*uvaryl[4]*mu-5.34375*uvaryc[4]*mu-1.453125*uvarxr[4]*mu-1.453125*uvarxl[4]*mu-5.34375*uvarxc[4]*mu+1.190784930203603*uvarxr[2]*mu-1.190784930203603*uvarxl[2]*mu+1.190784930203603*uvaryr[1]*mu-1.190784930203603*uvaryl[1]*mu);
  out[25] += J*(0.6708203932499369*uvaryr[18]*mu+0.6708203932499369*uvaryl[18]*mu-1.341640786499874*uvaryc[18]*mu+0.7382874503707888*uvarzr[15]*mu-0.7382874503707888*uvarzl[15]*mu+0.7382874503707888*uvarxr[13]*mu-0.7382874503707888*uvarxl[13]*mu-1.190784930203603*uvaryr[10]*mu+1.190784930203603*uvaryl[10]*mu-1.453125*uvarzr[5]*mu-1.453125*uvarzl[5]*mu-5.34375*uvarzc[5]*mu+0.9375*uvaryr[5]*mu+0.9375*uvaryl[5]*mu-1.875*uvaryc[5]*mu-1.453125*uvarxr[5]*mu-1.453125*uvarxl[5]*mu-5.34375*uvarxc[5]*mu+1.190784930203603*uvarxr[3]*mu-1.190784930203603*uvarxl[3]*mu+1.190784930203603*uvarzr[1]*mu-1.190784930203603*uvarzl[1]*mu);
  out[26] += J*(0.6708203932499369*uvarxr[17]*mu+0.6708203932499369*uvarxl[17]*mu-1.341640786499874*uvarxc[17]*mu+0.7382874503707888*uvarzr[16]*mu-0.7382874503707888*uvarzl[16]*mu+0.7382874503707888*uvaryr[14]*mu-0.7382874503707888*uvaryl[14]*mu-1.190784930203603*uvarxr[10]*mu+1.190784930203603*uvarxl[10]*mu-1.453125*uvarzr[6]*mu-1.453125*uvarzl[6]*mu-5.34375*uvarzc[6]*mu-1.453125*uvaryr[6]*mu-1.453125*uvaryl[6]*mu-5.34375*uvaryc[6]*mu+0.9375*uvarxr[6]*mu+0.9375*uvarxl[6]*mu-1.875*uvarxc[6]*mu+1.190784930203603*uvaryr[3]*mu-1.190784930203603*uvaryl[3]*mu+1.190784930203603*uvarzr[2]*mu-1.190784930203603*uvarzl[2]*mu);
  out[27] += J*((-1.190784930203603*uvarzr[13]*mu)+1.190784930203603*uvarzl[13]*mu-1.190784930203603*uvaryr[11]*mu+1.190784930203603*uvaryl[11]*mu+0.9375*uvarzr[7]*mu+0.9375*uvarzl[7]*mu-1.875*uvarzc[7]*mu+0.9375*uvaryr[7]*mu+0.9375*uvaryl[7]*mu-1.875*uvaryc[7]*mu-0.140625*uvarxr[7]*mu-0.140625*uvarxl[7]*mu-6.28125*uvarxc[7]*mu-0.3025768239224545*uvarxr[1]*mu+0.3025768239224545*uvarxl[1]*mu+0.4192627457812106*uvarxr[0]*mu+0.4192627457812106*uvarxl[0]*mu-0.8385254915624212*uvarxc[0]*mu);
  out[28] += J*((-1.190784930203603*uvarzr[14]*mu)+1.190784930203603*uvarzl[14]*mu-1.190784930203603*uvarxr[12]*mu+1.190784930203603*uvarxl[12]*mu+0.9375*uvarzr[8]*mu+0.9375*uvarzl[8]*mu-1.875*uvarzc[8]*mu-0.140625*uvaryr[8]*mu-0.140625*uvaryl[8]*mu-6.28125*uvaryc[8]*mu+0.9375*uvarxr[8]*mu+0.9375*uvarxl[8]*mu-1.875*uvarxc[8]*mu-0.3025768239224545*uvaryr[2]*mu+0.3025768239224545*uvaryl[2]*mu+0.4192627457812106*uvaryr[0]*mu+0.4192627457812106*uvaryl[0]*mu-0.8385254915624212*uvaryc[0]*mu);
  out[29] += J*((-1.190784930203603*uvaryr[16]*mu)+1.190784930203603*uvaryl[16]*mu-1.190784930203603*uvarxr[15]*mu+1.190784930203603*uvarxl[15]*mu-0.140625*uvarzr[9]*mu-0.140625*uvarzl[9]*mu-6.28125*uvarzc[9]*mu+0.9375*uvaryr[9]*mu+0.9375*uvaryl[9]*mu-1.875*uvaryc[9]*mu+0.9375*uvarxr[9]*mu+0.9375*uvarxl[9]*mu-1.875*uvarxc[9]*mu-0.3025768239224545*uvarzr[3]*mu+0.3025768239224545*uvarzl[3]*mu+0.4192627457812106*uvarzr[0]*mu+0.4192627457812106*uvarzl[0]*mu-0.8385254915624212*uvarzc[0]*mu);
  out[30] += J*(0.7382874503707888*uvarzr[19]*mu-0.7382874503707888*uvarzl[19]*mu+0.7382874503707888*uvaryr[18]*mu-0.7382874503707888*uvaryl[18]*mu+0.7382874503707888*uvarxr[17]*mu-0.7382874503707888*uvarxl[17]*mu-1.453125*uvarzr[10]*mu-1.453125*uvarzl[10]*mu-5.34375*uvarzc[10]*mu-1.453125*uvaryr[10]*mu-1.453125*uvaryl[10]*mu-5.34375*uvaryc[10]*mu-1.453125*uvarxr[10]*mu-1.453125*uvarxl[10]*mu-5.34375*uvarxc[10]*mu+1.190784930203603*uvarxr[6]*mu-1.190784930203603*uvarxl[6]*mu+1.190784930203603*uvaryr[5]*mu-1.190784930203603*uvaryl[5]*mu+1.190784930203603*uvarzr[4]*mu-1.190784930203603*uvarzl[4]*mu);
  out[31] += J*((-1.190784930203603*uvarzr[17]*mu)+1.190784930203603*uvarzl[17]*mu+0.9375*uvarzr[11]*mu+0.9375*uvarzl[11]*mu-1.875*uvarzc[11]*mu-1.453125*uvaryr[11]*mu-1.453125*uvaryl[11]*mu-5.34375*uvaryc[11]*mu-0.140625*uvarxr[11]*mu-0.140625*uvarxl[11]*mu-6.28125*uvarxc[11]*mu+1.190784930203603*uvaryr[7]*mu-1.190784930203603*uvaryl[7]*mu-0.3025768239224544*uvarxr[4]*mu+0.3025768239224544*uvarxl[4]*mu+0.4192627457812105*uvarxr[2]*mu+0.4192627457812105*uvarxl[2]*mu-0.8385254915624211*uvarxc[2]*mu);
  out[32] += J*((-1.190784930203603*uvarzr[18]*mu)+1.190784930203603*uvarzl[18]*mu+0.9375*uvarzr[12]*mu+0.9375*uvarzl[12]*mu-1.875*uvarzc[12]*mu-0.140625*uvaryr[12]*mu-0.140625*uvaryl[12]*mu-6.28125*uvaryc[12]*mu-1.453125*uvarxr[12]*mu-1.453125*uvarxl[12]*mu-5.34375*uvarxc[12]*mu+1.190784930203603*uvarxr[8]*mu-1.190784930203603*uvarxl[8]*mu-0.3025768239224544*uvaryr[4]*mu+0.3025768239224544*uvaryl[4]*mu+0.4192627457812105*uvaryr[1]*mu+0.4192627457812105*uvaryl[1]*mu-0.8385254915624211*uvaryc[1]*mu);
  out[33] += J*((-1.190784930203603*uvaryr[17]*mu)+1.190784930203603*uvaryl[17]*mu-1.453125*uvarzr[13]*mu-1.453125*uvarzl[13]*mu-5.34375*uvarzc[13]*mu+0.9375*uvaryr[13]*mu+0.9375*uvaryl[13]*mu-1.875*uvaryc[13]*mu-0.140625*uvarxr[13]*mu-0.140625*uvarxl[13]*mu-6.28125*uvarxc[13]*mu+1.190784930203603*uvarzr[7]*mu-1.190784930203603*uvarzl[7]*mu-0.3025768239224544*uvarxr[5]*mu+0.3025768239224544*uvarxl[5]*mu+0.4192627457812105*uvarxr[3]*mu+0.4192627457812105*uvarxl[3]*mu-0.8385254915624211*uvarxc[3]*mu);
  out[34] += J*((-1.190784930203603*uvarxr[18]*mu)+1.190784930203603*uvarxl[18]*mu-1.453125*uvarzr[14]*mu-1.453125*uvarzl[14]*mu-5.34375*uvarzc[14]*mu-0.140625*uvaryr[14]*mu-0.140625*uvaryl[14]*mu-6.28125*uvaryc[14]*mu+0.9375*uvarxr[14]*mu+0.9375*uvarxl[14]*mu-1.875*uvarxc[14]*mu+1.190784930203603*uvarzr[8]*mu-1.190784930203603*uvarzl[8]*mu-0.3025768239224544*uvaryr[6]*mu+0.3025768239224544*uvaryl[6]*mu+0.4192627457812105*uvaryr[3]*mu+0.4192627457812105*uvaryl[3]*mu-0.8385254915624211*uvaryc[3]*mu);
  out[35] += J*((-1.190784930203603*uvaryr[19]*mu)+1.190784930203603*uvaryl[19]*mu-0.140625*uvarzr[15]*mu-0.140625*uvarzl[15]*mu-6.28125*uvarzc[15]*mu+0.9375*uvaryr[15]*mu+0.9375*uvaryl[15]*mu-1.875*uvaryc[15]*mu-1.453125*uvarxr[15]*mu-1.453125*uvarxl[15]*mu-5.34375*uvarxc[15]*mu+1.190784930203603*uvarxr[9]*mu-1.190784930203603*uvarxl[9]*mu-0.3025768239224544*uvarzr[5]*mu+0.3025768239224544*uvarzl[5]*mu+0.4192627457812105*uvarzr[1]*mu+0.4192627457812105*uvarzl[1]*mu-0.8385254915624211*uvarzc[1]*mu);
  out[36] += J*((-1.190784930203603*uvarxr[19]*mu)+1.190784930203603*uvarxl[19]*mu-0.140625*uvarzr[16]*mu-0.140625*uvarzl[16]*mu-6.28125*uvarzc[16]*mu-1.453125*uvaryr[16]*mu-1.453125*uvaryl[16]*mu-5.34375*uvaryc[16]*mu+0.9375*uvarxr[16]*mu+0.9375*uvarxl[16]*mu-1.875*uvarxc[16]*mu+1.190784930203603*uvaryr[9]*mu-1.190784930203603*uvaryl[9]*mu-0.3025768239224544*uvarzr[6]*mu+0.3025768239224544*uvarzl[6]*mu+0.4192627457812105*uvarzr[2]*mu+0.4192627457812105*uvarzl[2]*mu-0.8385254915624211*uvarzc[2]*mu);
  out[37] += J*((-1.453125*uvarzr[17]*mu)-1.453125*uvarzl[17]*mu-5.34375*uvarzc[17]*mu-1.453125*uvaryr[17]*mu-1.453125*uvaryl[17]*mu-5.34375*uvaryc[17]*mu-0.140625*uvarxr[17]*mu-0.140625*uvarxl[17]*mu-6.28125*uvarxc[17]*mu+1.190784930203603*uvaryr[13]*mu-1.190784930203603*uvaryl[13]*mu+1.190784930203603*uvarzr[11]*mu-1.190784930203603*uvarzl[11]*mu-0.3025768239224545*uvarxr[10]*mu+0.3025768239224545*uvarxl[10]*mu+0.4192627457812106*uvarxr[6]*mu+0.4192627457812106*uvarxl[6]*mu-0.8385254915624212*uvarxc[6]*mu);
  out[38] += J*((-1.453125*uvarzr[18]*mu)-1.453125*uvarzl[18]*mu-5.34375*uvarzc[18]*mu-0.140625*uvaryr[18]*mu-0.140625*uvaryl[18]*mu-6.28125*uvaryc[18]*mu-1.453125*uvarxr[18]*mu-1.453125*uvarxl[18]*mu-5.34375*uvarxc[18]*mu+1.190784930203603*uvarxr[14]*mu-1.190784930203603*uvarxl[14]*mu+1.190784930203603*uvarzr[12]*mu-1.190784930203603*uvarzl[12]*mu-0.3025768239224545*uvaryr[10]*mu+0.3025768239224545*uvaryl[10]*mu+0.4192627457812106*uvaryr[5]*mu+0.4192627457812106*uvaryl[5]*mu-0.8385254915624212*uvaryc[5]*mu);
  out[39] += J*((-0.140625*uvarzr[19]*mu)-0.140625*uvarzl[19]*mu-6.28125*uvarzc[19]*mu-1.453125*uvaryr[19]*mu-1.453125*uvaryl[19]*mu-5.34375*uvaryc[19]*mu-1.453125*uvarxr[19]*mu-1.453125*uvarxl[19]*mu-5.34375*uvarxc[19]*mu+1.190784930203603*uvarxr[16]*mu-1.190784930203603*uvarxl[16]*mu+1.190784930203603*uvaryr[15]*mu-1.190784930203603*uvaryl[15]*mu-0.3025768239224545*uvarzr[10]*mu+0.3025768239224545*uvarzl[10]*mu+0.4192627457812106*uvarzr[4]*mu+0.4192627457812106*uvarzl[4]*mu-0.8385254915624212*uvarzc[4]*mu);
  out[40] += J*(0.6708203932499369*uvarzr[9]*mu+0.6708203932499369*uvarzl[9]*mu-1.341640786499874*uvarzc[9]*mu+0.6708203932499369*uvaryr[8]*mu+0.6708203932499369*uvaryl[8]*mu-1.341640786499874*uvaryc[8]*mu+0.6708203932499369*uvarxr[7]*mu+0.6708203932499369*uvarxl[7]*mu-1.341640786499874*uvarxc[7]*mu-1.190784930203603*uvarzr[3]*mu+1.190784930203603*uvarzl[3]*mu-1.190784930203603*uvaryr[2]*mu+1.190784930203603*uvaryl[2]*mu-1.190784930203603*uvarxr[1]*mu+1.190784930203603*uvarxl[1]*mu+0.9375*uvarzr[0]*mu+0.9375*uvarzl[0]*mu-1.875*uvarzc[0]*mu+0.9375*uvaryr[0]*mu+0.9375*uvaryl[0]*mu-1.875*uvaryc[0]*mu+0.9375*uvarxr[0]*mu+0.9375*uvarxl[0]*mu-1.875*uvarxc[0]*mu);
  out[41] += J*(0.6708203932499369*uvarzr[15]*mu+0.6708203932499369*uvarzl[15]*mu-1.341640786499874*uvarzc[15]*mu+0.6708203932499369*uvaryr[12]*mu+0.6708203932499369*uvaryl[12]*mu-1.341640786499874*uvaryc[12]*mu+0.7382874503707888*uvarxr[7]*mu-0.7382874503707888*uvarxl[7]*mu-1.190784930203603*uvarzr[5]*mu+1.190784930203603*uvarzl[5]*mu-1.190784930203603*uvaryr[4]*mu+1.190784930203603*uvaryl[4]*mu+0.9375*uvarzr[1]*mu+0.9375*uvarzl[1]*mu-1.875*uvarzc[1]*mu+0.9375*uvaryr[1]*mu+0.9375*uvaryl[1]*mu-1.875*uvaryc[1]*mu-1.453125*uvarxr[1]*mu-1.453125*uvarxl[1]*mu-5.34375*uvarxc[1]*mu+1.190784930203603*uvarxr[0]*mu-1.190784930203603*uvarxl[0]*mu);
  out[42] += J*(0.6708203932499369*uvarzr[16]*mu+0.6708203932499369*uvarzl[16]*mu-1.341640786499874*uvarzc[16]*mu+0.6708203932499369*uvarxr[11]*mu+0.6708203932499369*uvarxl[11]*mu-1.341640786499874*uvarxc[11]*mu+0.7382874503707888*uvaryr[8]*mu-0.7382874503707888*uvaryl[8]*mu-1.190784930203603*uvarzr[6]*mu+1.190784930203603*uvarzl[6]*mu-1.190784930203603*uvarxr[4]*mu+1.190784930203603*uvarxl[4]*mu+0.9375*uvarzr[2]*mu+0.9375*uvarzl[2]*mu-1.875*uvarzc[2]*mu-1.453125*uvaryr[2]*mu-1.453125*uvaryl[2]*mu-5.34375*uvaryc[2]*mu+0.9375*uvarxr[2]*mu+0.9375*uvarxl[2]*mu-1.875*uvarxc[2]*mu+1.190784930203603*uvaryr[0]*mu-1.190784930203603*uvaryl[0]*mu);
  out[43] += J*(0.6708203932499369*uvaryr[14]*mu+0.6708203932499369*uvaryl[14]*mu-1.341640786499874*uvaryc[14]*mu+0.6708203932499369*uvarxr[13]*mu+0.6708203932499369*uvarxl[13]*mu-1.341640786499874*uvarxc[13]*mu+0.7382874503707888*uvarzr[9]*mu-0.7382874503707888*uvarzl[9]*mu-1.190784930203603*uvaryr[6]*mu+1.190784930203603*uvaryl[6]*mu-1.190784930203603*uvarxr[5]*mu+1.190784930203603*uvarxl[5]*mu-1.453125*uvarzr[3]*mu-1.453125*uvarzl[3]*mu-5.34375*uvarzc[3]*mu+0.9375*uvaryr[3]*mu+0.9375*uvaryl[3]*mu-1.875*uvaryc[3]*mu+0.9375*uvarxr[3]*mu+0.9375*uvarxl[3]*mu-1.875*uvarxc[3]*mu+1.190784930203603*uvarzr[0]*mu-1.190784930203603*uvarzl[0]*mu);
  out[44] += J*(0.6708203932499369*uvarzr[19]*mu+0.6708203932499369*uvarzl[19]*mu-1.341640786499874*uvarzc[19]*mu+0.7382874503707888*uvaryr[12]*mu-0.7382874503707888*uvaryl[12]*mu+0.7382874503707888*uvarxr[11]*mu-0.7382874503707888*uvarxl[11]*mu-1.190784930203603*uvarzr[10]*mu+1.190784930203603*uvarzl[10]*mu+0.9375*uvarzr[4]*mu+0.9375*uvarzl[4]*mu-1.875*uvarzc[4]*mu-1.453125*uvaryr[4]*mu-1.453125*uvaryl[4]*mu-5.34375*uvaryc[4]*mu-1.453125*uvarxr[4]*mu-1.453125*uvarxl[4]*mu-5.34375*uvarxc[4]*mu+1.190784930203603*uvarxr[2]*mu-1.190784930203603*uvarxl[2]*mu+1.190784930203603*uvaryr[1]*mu-1.190784930203603*uvaryl[1]*mu);
  out[45] += J*(0.6708203932499369*uvaryr[18]*mu+0.6708203932499369*uvaryl[18]*mu-1.341640786499874*uvaryc[18]*mu+0.7382874503707888*uvarzr[15]*mu-0.7382874503707888*uvarzl[15]*mu+0.7382874503707888*uvarxr[13]*mu-0.7382874503707888*uvarxl[13]*mu-1.190784930203603*uvaryr[10]*mu+1.190784930203603*uvaryl[10]*mu-1.453125*uvarzr[5]*mu-1.453125*uvarzl[5]*mu-5.34375*uvarzc[5]*mu+0.9375*uvaryr[5]*mu+0.9375*uvaryl[5]*mu-1.875*uvaryc[5]*mu-1.453125*uvarxr[5]*mu-1.453125*uvarxl[5]*mu-5.34375*uvarxc[5]*mu+1.190784930203603*uvarxr[3]*mu-1.190784930203603*uvarxl[3]*mu+1.190784930203603*uvarzr[1]*mu-1.190784930203603*uvarzl[1]*mu);
  out[46] += J*(0.6708203932499369*uvarxr[17]*mu+0.6708203932499369*uvarxl[17]*mu-1.341640786499874*uvarxc[17]*mu+0.7382874503707888*uvarzr[16]*mu-0.7382874503707888*uvarzl[16]*mu+0.7382874503707888*uvaryr[14]*mu-0.7382874503707888*uvaryl[14]*mu-1.190784930203603*uvarxr[10]*mu+1.190784930203603*uvarxl[10]*mu-1.453125*uvarzr[6]*mu-1.453125*uvarzl[6]*mu-5.34375*uvarzc[6]*mu-1.453125*uvaryr[6]*mu-1.453125*uvaryl[6]*mu-5.34375*uvaryc[6]*mu+0.9375*uvarxr[6]*mu+0.9375*uvarxl[6]*mu-1.875*uvarxc[6]*mu+1.190784930203603*uvaryr[3]*mu-1.190784930203603*uvaryl[3]*mu+1.190784930203603*uvarzr[2]*mu-1.190784930203603*uvarzl[2]*mu);
  out[47] += J*((-1.190784930203603*uvarzr[13]*mu)+1.190784930203603*uvarzl[13]*mu-1.190784930203603*uvaryr[11]*mu+1.190784930203603*uvaryl[11]*mu+0.9375*uvarzr[7]*mu+0.9375*uvarzl[7]*mu-1.875*uvarzc[7]*mu+0.9375*uvaryr[7]*mu+0.9375*uvaryl[7]*mu-1.875*uvaryc[7]*mu-0.140625*uvarxr[7]*mu-0.140625*uvarxl[7]*mu-6.28125*uvarxc[7]*mu-0.3025768239224545*uvarxr[1]*mu+0.3025768239224545*uvarxl[1]*mu+0.4192627457812106*uvarxr[0]*mu+0.4192627457812106*uvarxl[0]*mu-0.8385254915624212*uvarxc[0]*mu);
  out[48] += J*((-1.190784930203603*uvarzr[14]*mu)+1.190784930203603*uvarzl[14]*mu-1.190784930203603*uvarxr[12]*mu+1.190784930203603*uvarxl[12]*mu+0.9375*uvarzr[8]*mu+0.9375*uvarzl[8]*mu-1.875*uvarzc[8]*mu-0.140625*uvaryr[8]*mu-0.140625*uvaryl[8]*mu-6.28125*uvaryc[8]*mu+0.9375*uvarxr[8]*mu+0.9375*uvarxl[8]*mu-1.875*uvarxc[8]*mu-0.3025768239224545*uvaryr[2]*mu+0.3025768239224545*uvaryl[2]*mu+0.4192627457812106*uvaryr[0]*mu+0.4192627457812106*uvaryl[0]*mu-0.8385254915624212*uvaryc[0]*mu);
  out[49] += J*((-1.190784930203603*uvaryr[16]*mu)+1.190784930203603*uvaryl[16]*mu-1.190784930203603*uvarxr[15]*mu+1.190784930203603*uvarxl[15]*mu-0.140625*uvarzr[9]*mu-0.140625*uvarzl[9]*mu-6.28125*uvarzc[9]*mu+0.9375*uvaryr[9]*mu+0.9375*uvaryl[9]*mu-1.875*uvaryc[9]*mu+0.9375*uvarxr[9]*mu+0.9375*uvarxl[9]*mu-1.875*uvarxc[9]*mu-0.3025768239224545*uvarzr[3]*mu+0.3025768239224545*uvarzl[3]*mu+0.4192627457812106*uvarzr[0]*mu+0.4192627457812106*uvarzl[0]*mu-0.8385254915624212*uvarzc[0]*mu);
  out[50] += J*(0.7382874503707888*uvarzr[19]*mu-0.7382874503707888*uvarzl[19]*mu+0.7382874503707888*uvaryr[18]*mu-0.7382874503707888*uvaryl[18]*mu+0.7382874503707888*uvarxr[17]*mu-0.7382874503707888*uvarxl[17]*mu-1.453125*uvarzr[10]*mu-1.453125*uvarzl[10]*mu-5.34375*uvarzc[10]*mu-1.453125*uvaryr[10]*mu-1.453125*uvaryl[10]*mu-5.34375*uvaryc[10]*mu-1.453125*uvarxr[10]*mu-1.453125*uvarxl[10]*mu-5.34375*uvarxc[10]*mu+1.190784930203603*uvarxr[6]*mu-1.190784930203603*uvarxl[6]*mu+1.190784930203603*uvaryr[5]*mu-1.190784930203603*uvaryl[5]*mu+1.190784930203603*uvarzr[4]*mu-1.190784930203603*uvarzl[4]*mu);
  out[51] += J*((-1.190784930203603*uvarzr[17]*mu)+1.190784930203603*uvarzl[17]*mu+0.9375*uvarzr[11]*mu+0.9375*uvarzl[11]*mu-1.875*uvarzc[11]*mu-1.453125*uvaryr[11]*mu-1.453125*uvaryl[11]*mu-5.34375*uvaryc[11]*mu-0.140625*uvarxr[11]*mu-0.140625*uvarxl[11]*mu-6.28125*uvarxc[11]*mu+1.190784930203603*uvaryr[7]*mu-1.190784930203603*uvaryl[7]*mu-0.3025768239224544*uvarxr[4]*mu+0.3025768239224544*uvarxl[4]*mu+0.4192627457812105*uvarxr[2]*mu+0.4192627457812105*uvarxl[2]*mu-0.8385254915624211*uvarxc[2]*mu);
  out[52] += J*((-1.190784930203603*uvarzr[18]*mu)+1.190784930203603*uvarzl[18]*mu+0.9375*uvarzr[12]*mu+0.9375*uvarzl[12]*mu-1.875*uvarzc[12]*mu-0.140625*uvaryr[12]*mu-0.140625*uvaryl[12]*mu-6.28125*uvaryc[12]*mu-1.453125*uvarxr[12]*mu-1.453125*uvarxl[12]*mu-5.34375*uvarxc[12]*mu+1.190784930203603*uvarxr[8]*mu-1.190784930203603*uvarxl[8]*mu-0.3025768239224544*uvaryr[4]*mu+0.3025768239224544*uvaryl[4]*mu+0.4192627457812105*uvaryr[1]*mu+0.4192627457812105*uvaryl[1]*mu-0.8385254915624211*uvaryc[1]*mu);
  out[53] += J*((-1.190784930203603*uvaryr[17]*mu)+1.190784930203603*uvaryl[17]*mu-1.453125*uvarzr[13]*mu-1.453125*uvarzl[13]*mu-5.34375*uvarzc[13]*mu+0.9375*uvaryr[13]*mu+0.9375*uvaryl[13]*mu-1.875*uvaryc[13]*mu-0.140625*uvarxr[13]*mu-0.140625*uvarxl[13]*mu-6.28125*uvarxc[13]*mu+1.190784930203603*uvarzr[7]*mu-1.190784930203603*uvarzl[7]*mu-0.3025768239224544*uvarxr[5]*mu+0.3025768239224544*uvarxl[5]*mu+0.4192627457812105*uvarxr[3]*mu+0.4192627457812105*uvarxl[3]*mu-0.8385254915624211*uvarxc[3]*mu);
  out[54] += J*((-1.190784930203603*uvarxr[18]*mu)+1.190784930203603*uvarxl[18]*mu-1.453125*uvarzr[14]*mu-1.453125*uvarzl[14]*mu-5.34375*uvarzc[14]*mu-0.140625*uvaryr[14]*mu-0.140625*uvaryl[14]*mu-6.28125*uvaryc[14]*mu+0.9375*uvarxr[14]*mu+0.9375*uvarxl[14]*mu-1.875*uvarxc[14]*mu+1.190784930203603*uvarzr[8]*mu-1.190784930203603*uvarzl[8]*mu-0.3025768239224544*uvaryr[6]*mu+0.3025768239224544*uvaryl[6]*mu+0.4192627457812105*uvaryr[3]*mu+0.4192627457812105*uvaryl[3]*mu-0.8385254915624211*uvaryc[3]*mu);
  out[55] += J*((-1.190784930203603*uvaryr[19]*mu)+1.190784930203603*uvaryl[19]*mu-0.140625*uvarzr[15]*mu-0.140625*uvarzl[15]*mu-6.28125*uvarzc[15]*mu+0.9375*uvaryr[15]*mu+0.9375*uvaryl[15]*mu-1.875*uvaryc[15]*mu-1.453125*uvarxr[15]*mu-1.453125*uvarxl[15]*mu-5.34375*uvarxc[15]*mu+1.190784930203603*uvarxr[9]*mu-1.190784930203603*uvarxl[9]*mu-0.3025768239224544*uvarzr[5]*mu+0.3025768239224544*uvarzl[5]*mu+0.4192627457812105*uvarzr[1]*mu+0.4192627457812105*uvarzl[1]*mu-0.8385254915624211*uvarzc[1]*mu);
  out[56] += J*((-1.190784930203603*uvarxr[19]*mu)+1.190784930203603*uvarxl[19]*mu-0.140625*uvarzr[16]*mu-0.140625*uvarzl[16]*mu-6.28125*uvarzc[16]*mu-1.453125*uvaryr[16]*mu-1.453125*uvaryl[16]*mu-5.34375*uvaryc[16]*mu+0.9375*uvarxr[16]*mu+0.9375*uvarxl[16]*mu-1.875*uvarxc[16]*mu+1.190784930203603*uvaryr[9]*mu-1.190784930203603*uvaryl[9]*mu-0.3025768239224544*uvarzr[6]*mu+0.3025768239224544*uvarzl[6]*mu+0.4192627457812105*uvarzr[2]*mu+0.4192627457812105*uvarzl[2]*mu-0.8385254915624211*uvarzc[2]*mu);
  out[57] += J*((-1.453125*uvarzr[17]*mu)-1.453125*uvarzl[17]*mu-5.34375*uvarzc[17]*mu-1.453125*uvaryr[17]*mu-1.453125*uvaryl[17]*mu-5.34375*uvaryc[17]*mu-0.140625*uvarxr[17]*mu-0.140625*uvarxl[17]*mu-6.28125*uvarxc[17]*mu+1.190784930203603*uvaryr[13]*mu-1.190784930203603*uvaryl[13]*mu+1.190784930203603*uvarzr[11]*mu-1.190784930203603*uvarzl[11]*mu-0.3025768239224545*uvarxr[10]*mu+0.3025768239224545*uvarxl[10]*mu+0.4192627457812106*uvarxr[6]*mu+0.4192627457812106*uvarxl[6]*mu-0.8385254915624212*uvarxc[6]*mu);
  out[58] += J*((-1.453125*uvarzr[18]*mu)-1.453125*uvarzl[18]*mu-5.34375*uvarzc[18]*mu-0.140625*uvaryr[18]*mu-0.140625*uvaryl[18]*mu-6.28125*uvaryc[18]*mu-1.453125*uvarxr[18]*mu-1.453125*uvarxl[18]*mu-5.34375*uvarxc[18]*mu+1.190784930203603*uvarxr[14]*mu-1.190784930203603*uvarxl[14]*mu+1.190784930203603*uvarzr[12]*mu-1.190784930203603*uvarzl[12]*mu-0.3025768239224545*uvaryr[10]*mu+0.3025768239224545*uvaryl[10]*mu+0.4192627457812106*uvaryr[5]*mu+0.4192627457812106*uvaryl[5]*mu-0.8385254915624212*uvaryc[5]*mu);
  out[59] += J*((-0.140625*uvarzr[19]*mu)-0.140625*uvarzl[19]*mu-6.28125*uvarzc[19]*mu-1.453125*uvaryr[19]*mu-1.453125*uvaryl[19]*mu-5.34375*uvaryc[19]*mu-1.453125*uvarxr[19]*mu-1.453125*uvarxl[19]*mu-5.34375*uvarxc[19]*mu+1.190784930203603*uvarxr[16]*mu-1.190784930203603*uvarxl[16]*mu+1.190784930203603*uvaryr[15]*mu-1.190784930203603*uvaryl[15]*mu-0.3025768239224545*uvarzr[10]*mu+0.3025768239224545*uvarzl[10]*mu+0.4192627457812106*uvarzr[4]*mu+0.4192627457812106*uvarzl[4]*mu-0.8385254915624212*uvarzc[4]*mu);
  out[60] += J*(0.6708203932499369*uvarzr[9]*mu+0.6708203932499369*uvarzl[9]*mu-1.341640786499874*uvarzc[9]*mu+0.6708203932499369*uvaryr[8]*mu+0.6708203932499369*uvaryl[8]*mu-1.341640786499874*uvaryc[8]*mu+0.6708203932499369*uvarxr[7]*mu+0.6708203932499369*uvarxl[7]*mu-1.341640786499874*uvarxc[7]*mu-1.190784930203603*uvarzr[3]*mu+1.190784930203603*uvarzl[3]*mu-1.190784930203603*uvaryr[2]*mu+1.190784930203603*uvaryl[2]*mu-1.190784930203603*uvarxr[1]*mu+1.190784930203603*uvarxl[1]*mu+0.9375*uvarzr[0]*mu+0.9375*uvarzl[0]*mu-1.875*uvarzc[0]*mu+0.9375*uvaryr[0]*mu+0.9375*uvaryl[0]*mu-1.875*uvaryc[0]*mu+0.9375*uvarxr[0]*mu+0.9375*uvarxl[0]*mu-1.875*uvarxc[0]*mu);
  out[61] += J*(0.6708203932499369*uvarzr[15]*mu+0.6708203932499369*uvarzl[15]*mu-1.341640786499874*uvarzc[15]*mu+0.6708203932499369*uvaryr[12]*mu+0.6708203932499369*uvaryl[12]*mu-1.341640786499874*uvaryc[12]*mu+0.7382874503707888*uvarxr[7]*mu-0.7382874503707888*uvarxl[7]*mu-1.190784930203603*uvarzr[5]*mu+1.190784930203603*uvarzl[5]*mu-1.190784930203603*uvaryr[4]*mu+1.190784930203603*uvaryl[4]*mu+0.9375*uvarzr[1]*mu+0.9375*uvarzl[1]*mu-1.875*uvarzc[1]*mu+0.9375*uvaryr[1]*mu+0.9375*uvaryl[1]*mu-1.875*uvaryc[1]*mu-1.453125*uvarxr[1]*mu-1.453125*uvarxl[1]*mu-5.34375*uvarxc[1]*mu+1.190784930203603*uvarxr[0]*mu-1.190784930203603*uvarxl[0]*mu);
  out[62] += J*(0.6708203932499369*uvarzr[16]*mu+0.6708203932499369*uvarzl[16]*mu-1.341640786499874*uvarzc[16]*mu+0.6708203932499369*uvarxr[11]*mu+0.6708203932499369*uvarxl[11]*mu-1.341640786499874*uvarxc[11]*mu+0.7382874503707888*uvaryr[8]*mu-0.7382874503707888*uvaryl[8]*mu-1.190784930203603*uvarzr[6]*mu+1.190784930203603*uvarzl[6]*mu-1.190784930203603*uvarxr[4]*mu+1.190784930203603*uvarxl[4]*mu+0.9375*uvarzr[2]*mu+0.9375*uvarzl[2]*mu-1.875*uvarzc[2]*mu-1.453125*uvaryr[2]*mu-1.453125*uvaryl[2]*mu-5.34375*uvaryc[2]*mu+0.9375*uvarxr[2]*mu+0.9375*uvarxl[2]*mu-1.875*uvarxc[2]*mu+1.190784930203603*uvaryr[0]*mu-1.190784930203603*uvaryl[0]*mu);
  out[63] += J*(0.6708203932499369*uvaryr[14]*mu+0.6708203932499369*uvaryl[14]*mu-1.341640786499874*uvaryc[14]*mu+0.6708203932499369*uvarxr[13]*mu+0.6708203932499369*uvarxl[13]*mu-1.341640786499874*uvarxc[13]*mu+0.7382874503707888*uvarzr[9]*mu-0.7382874503707888*uvarzl[9]*mu-1.190784930203603*uvaryr[6]*mu+1.190784930203603*uvaryl[6]*mu-1.190784930203603*uvarxr[5]*mu+1.190784930203603*uvarxl[5]*mu-1.453125*uvarzr[3]*mu-1.453125*uvarzl[3]*mu-5.34375*uvarzc[3]*mu+0.9375*uvaryr[3]*mu+0.9375*uvaryl[3]*mu-1.875*uvaryc[3]*mu+0.9375*uvarxr[3]*mu+0.9375*uvarxl[3]*mu-1.875*uvarxc[3]*mu+1.190784930203603*uvarzr[0]*mu-1.190784930203603*uvarzl[0]*mu);
  out[64] += J*(0.6708203932499369*uvarzr[19]*mu+0.6708203932499369*uvarzl[19]*mu-1.341640786499874*uvarzc[19]*mu+0.7382874503707888*uvaryr[12]*mu-0.7382874503707888*uvaryl[12]*mu+0.7382874503707888*uvarxr[11]*mu-0.7382874503707888*uvarxl[11]*mu-1.190784930203603*uvarzr[10]*mu+1.190784930203603*uvarzl[10]*mu+0.9375*uvarzr[4]*mu+0.9375*uvarzl[4]*mu-1.875*uvarzc[4]*mu-1.453125*uvaryr[4]*mu-1.453125*uvaryl[4]*mu-5.34375*uvaryc[4]*mu-1.453125*uvarxr[4]*mu-1.453125*uvarxl[4]*mu-5.34375*uvarxc[4]*mu+1.190784930203603*uvarxr[2]*mu-1.190784930203603*uvarxl[2]*mu+1.190784930203603*uvaryr[1]*mu-1.190784930203603*uvaryl[1]*mu);
  out[65] += J*(0.6708203932499369*uvaryr[18]*mu+0.6708203932499369*uvaryl[18]*mu-1.341640786499874*uvaryc[18]*mu+0.7382874503707888*uvarzr[15]*mu-0.7382874503707888*uvarzl[15]*mu+0.7382874503707888*uvarxr[13]*mu-0.7382874503707888*uvarxl[13]*mu-1.190784930203603*uvaryr[10]*mu+1.190784930203603*uvaryl[10]*mu-1.453125*uvarzr[5]*mu-1.453125*uvarzl[5]*mu-5.34375*uvarzc[5]*mu+0.9375*uvaryr[5]*mu+0.9375*uvaryl[5]*mu-1.875*uvaryc[5]*mu-1.453125*uvarxr[5]*mu-1.453125*uvarxl[5]*mu-5.34375*uvarxc[5]*mu+1.190784930203603*uvarxr[3]*mu-1.190784930203603*uvarxl[3]*mu+1.190784930203603*uvarzr[1]*mu-1.190784930203603*uvarzl[1]*mu);
  out[66] += J*(0.6708203932499369*uvarxr[17]*mu+0.6708203932499369*uvarxl[17]*mu-1.341640786499874*uvarxc[17]*mu+0.7382874503707888*uvarzr[16]*mu-0.7382874503707888*uvarzl[16]*mu+0.7382874503707888*uvaryr[14]*mu-0.7382874503707888*uvaryl[14]*mu-1.190784930203603*uvarxr[10]*mu+1.190784930203603*uvarxl[10]*mu-1.453125*uvarzr[6]*mu-1.453125*uvarzl[6]*mu-5.34375*uvarzc[6]*mu-1.453125*uvaryr[6]*mu-1.453125*uvaryl[6]*mu-5.34375*uvaryc[6]*mu+0.9375*uvarxr[6]*mu+0.9375*uvarxl[6]*mu-1.875*uvarxc[6]*mu+1.190784930203603*uvaryr[3]*mu-1.190784930203603*uvaryl[3]*mu+1.190784930203603*uvarzr[2]*mu-1.190784930203603*uvarzl[2]*mu);
  out[67] += J*((-1.190784930203603*uvarzr[13]*mu)+1.190784930203603*uvarzl[13]*mu-1.190784930203603*uvaryr[11]*mu+1.190784930203603*uvaryl[11]*mu+0.9375*uvarzr[7]*mu+0.9375*uvarzl[7]*mu-1.875*uvarzc[7]*mu+0.9375*uvaryr[7]*mu+0.9375*uvaryl[7]*mu-1.875*uvaryc[7]*mu-0.140625*uvarxr[7]*mu-0.140625*uvarxl[7]*mu-6.28125*uvarxc[7]*mu-0.3025768239224545*uvarxr[1]*mu+0.3025768239224545*uvarxl[1]*mu+0.4192627457812106*uvarxr[0]*mu+0.4192627457812106*uvarxl[0]*mu-0.8385254915624212*uvarxc[0]*mu);
  out[68] += J*((-1.190784930203603*uvarzr[14]*mu)+1.190784930203603*uvarzl[14]*mu-1.190784930203603*uvarxr[12]*mu+1.190784930203603*uvarxl[12]*mu+0.9375*uvarzr[8]*mu+0.9375*uvarzl[8]*mu-1.875*uvarzc[8]*mu-0.140625*uvaryr[8]*mu-0.140625*uvaryl[8]*mu-6.28125*uvaryc[8]*mu+0.9375*uvarxr[8]*mu+0.9375*uvarxl[8]*mu-1.875*uvarxc[8]*mu-0.3025768239224545*uvaryr[2]*mu+0.3025768239224545*uvaryl[2]*mu+0.4192627457812106*uvaryr[0]*mu+0.4192627457812106*uvaryl[0]*mu-0.8385254915624212*uvaryc[0]*mu);
  out[69] += J*((-1.190784930203603*uvaryr[16]*mu)+1.190784930203603*uvaryl[16]*mu-1.190784930203603*uvarxr[15]*mu+1.190784930203603*uvarxl[15]*mu-0.140625*uvarzr[9]*mu-0.140625*uvarzl[9]*mu-6.28125*uvarzc[9]*mu+0.9375*uvaryr[9]*mu+0.9375*uvaryl[9]*mu-1.875*uvaryc[9]*mu+0.9375*uvarxr[9]*mu+0.9375*uvarxl[9]*mu-1.875*uvarxc[9]*mu-0.3025768239224545*uvarzr[3]*mu+0.3025768239224545*uvarzl[3]*mu+0.4192627457812106*uvarzr[0]*mu+0.4192627457812106*uvarzl[0]*mu-0.8385254915624212*uvarzc[0]*mu);
  out[70] += J*(0.7382874503707888*uvarzr[19]*mu-0.7382874503707888*uvarzl[19]*mu+0.7382874503707888*uvaryr[18]*mu-0.7382874503707888*uvaryl[18]*mu+0.7382874503707888*uvarxr[17]*mu-0.7382874503707888*uvarxl[17]*mu-1.453125*uvarzr[10]*mu-1.453125*uvarzl[10]*mu-5.34375*uvarzc[10]*mu-1.453125*uvaryr[10]*mu-1.453125*uvaryl[10]*mu-5.34375*uvaryc[10]*mu-1.453125*uvarxr[10]*mu-1.453125*uvarxl[10]*mu-5.34375*uvarxc[10]*mu+1.190784930203603*uvarxr[6]*mu-1.190784930203603*uvarxl[6]*mu+1.190784930203603*uvaryr[5]*mu-1.190784930203603*uvaryl[5]*mu+1.190784930203603*uvarzr[4]*mu-1.190784930203603*uvarzl[4]*mu);
  out[71] += J*((-1.190784930203603*uvarzr[17]*mu)+1.190784930203603*uvarzl[17]*mu+0.9375*uvarzr[11]*mu+0.9375*uvarzl[11]*mu-1.875*uvarzc[11]*mu-1.453125*uvaryr[11]*mu-1.453125*uvaryl[11]*mu-5.34375*uvaryc[11]*mu-0.140625*uvarxr[11]*mu-0.140625*uvarxl[11]*mu-6.28125*uvarxc[11]*mu+1.190784930203603*uvaryr[7]*mu-1.190784930203603*uvaryl[7]*mu-0.3025768239224544*uvarxr[4]*mu+0.3025768239224544*uvarxl[4]*mu+0.4192627457812105*uvarxr[2]*mu+0.4192627457812105*uvarxl[2]*mu-0.8385254915624211*uvarxc[2]*mu);
  out[72] += J*((-1.190784930203603*uvarzr[18]*mu)+1.190784930203603*uvarzl[18]*mu+0.9375*uvarzr[12]*mu+0.9375*uvarzl[12]*mu-1.875*uvarzc[12]*mu-0.140625*uvaryr[12]*mu-0.140625*uvaryl[12]*mu-6.28125*uvaryc[12]*mu-1.453125*uvarxr[12]*mu-1.453125*uvarxl[12]*mu-5.34375*uvarxc[12]*mu+1.190784930203603*uvarxr[8]*mu-1.190784930203603*uvarxl[8]*mu-0.3025768239224544*uvaryr[4]*mu+0.3025768239224544*uvaryl[4]*mu+0.4192627457812105*uvaryr[1]*mu+0.4192627457812105*uvaryl[1]*mu-0.8385254915624211*uvaryc[1]*mu);
  out[73] += J*((-1.190784930203603*uvaryr[17]*mu)+1.190784930203603*uvaryl[17]*mu-1.453125*uvarzr[13]*mu-1.453125*uvarzl[13]*mu-5.34375*uvarzc[13]*mu+0.9375*uvaryr[13]*mu+0.9375*uvaryl[13]*mu-1.875*uvaryc[13]*mu-0.140625*uvarxr[13]*mu-0.140625*uvarxl[13]*mu-6.28125*uvarxc[13]*mu+1.190784930203603*uvarzr[7]*mu-1.190784930203603*uvarzl[7]*mu-0.3025768239224544*uvarxr[5]*mu+0.3025768239224544*uvarxl[5]*mu+0.4192627457812105*uvarxr[3]*mu+0.4192627457812105*uvarxl[3]*mu-0.8385254915624211*uvarxc[3]*mu);
  out[74] += J*((-1.190784930203603*uvarxr[18]*mu)+1.190784930203603*uvarxl[18]*mu-1.453125*uvarzr[14]*mu-1.453125*uvarzl[14]*mu-5.34375*uvarzc[14]*mu-0.140625*uvaryr[14]*mu-0.140625*uvaryl[14]*mu-6.28125*uvaryc[14]*mu+0.9375*uvarxr[14]*mu+0.9375*uvarxl[14]*mu-1.875*uvarxc[14]*mu+1.190784930203603*uvarzr[8]*mu-1.190784930203603*uvarzl[8]*mu-0.3025768239224544*uvaryr[6]*mu+0.3025768239224544*uvaryl[6]*mu+0.4192627457812105*uvaryr[3]*mu+0.4192627457812105*uvaryl[3]*mu-0.8385254915624211*uvaryc[3]*mu);
  out[75] += J*((-1.190784930203603*uvaryr[19]*mu)+1.190784930203603*uvaryl[19]*mu-0.140625*uvarzr[15]*mu-0.140625*uvarzl[15]*mu-6.28125*uvarzc[15]*mu+0.9375*uvaryr[15]*mu+0.9375*uvaryl[15]*mu-1.875*uvaryc[15]*mu-1.453125*uvarxr[15]*mu-1.453125*uvarxl[15]*mu-5.34375*uvarxc[15]*mu+1.190784930203603*uvarxr[9]*mu-1.190784930203603*uvarxl[9]*mu-0.3025768239224544*uvarzr[5]*mu+0.3025768239224544*uvarzl[5]*mu+0.4192627457812105*uvarzr[1]*mu+0.4192627457812105*uvarzl[1]*mu-0.8385254915624211*uvarzc[1]*mu);
  out[76] += J*((-1.190784930203603*uvarxr[19]*mu)+1.190784930203603*uvarxl[19]*mu-0.140625*uvarzr[16]*mu-0.140625*uvarzl[16]*mu-6.28125*uvarzc[16]*mu-1.453125*uvaryr[16]*mu-1.453125*uvaryl[16]*mu-5.34375*uvaryc[16]*mu+0.9375*uvarxr[16]*mu+0.9375*uvarxl[16]*mu-1.875*uvarxc[16]*mu+1.190784930203603*uvaryr[9]*mu-1.190784930203603*uvaryl[9]*mu-0.3025768239224544*uvarzr[6]*mu+0.3025768239224544*uvarzl[6]*mu+0.4192627457812105*uvarzr[2]*mu+0.4192627457812105*uvarzl[2]*mu-0.8385254915624211*uvarzc[2]*mu);
  out[77] += J*((-1.453125*uvarzr[17]*mu)-1.453125*uvarzl[17]*mu-5.34375*uvarzc[17]*mu-1.453125*uvaryr[17]*mu-1.453125*uvaryl[17]*mu-5.34375*uvaryc[17]*mu-0.140625*uvarxr[17]*mu-0.140625*uvarxl[17]*mu-6.28125*uvarxc[17]*mu+1.190784930203603*uvaryr[13]*mu-1.190784930203603*uvaryl[13]*mu+1.190784930203603*uvarzr[11]*mu-1.190784930203603*uvarzl[11]*mu-0.3025768239224545*uvarxr[10]*mu+0.3025768239224545*uvarxl[10]*mu+0.4192627457812106*uvarxr[6]*mu+0.4192627457812106*uvarxl[6]*mu-0.8385254915624212*uvarxc[6]*mu);
  out[78] += J*((-1.453125*uvarzr[18]*mu)-1.453125*uvarzl[18]*mu-5.34375*uvarzc[18]*mu-0.140625*uvaryr[18]*mu-0.140625*uvaryl[18]*mu-6.28125*uvaryc[18]*mu-1.453125*uvarxr[18]*mu-1.453125*uvarxl[18]*mu-5.34375*uvarxc[18]*mu+1.190784930203603*uvarxr[14]*mu-1.190784930203603*uvarxl[14]*mu+1.190784930203603*uvarzr[12]*mu-1.190784930203603*uvarzl[12]*mu-0.3025768239224545*uvaryr[10]*mu+0.3025768239224545*uvaryl[10]*mu+0.4192627457812106*uvaryr[5]*mu+0.4192627457812106*uvaryl[5]*mu-0.8385254915624212*uvaryc[5]*mu);
  out[79] += J*((-0.140625*uvarzr[19]*mu)-0.140625*uvarzl[19]*mu-6.28125*uvarzc[19]*mu-1.453125*uvaryr[19]*mu-1.453125*uvaryl[19]*mu-5.34375*uvaryc[19]*mu-1.453125*uvarxr[19]*mu-1.453125*uvarxl[19]*mu-5.34375*uvarxc[19]*mu+1.190784930203603*uvarxr[16]*mu-1.190784930203603*uvarxl[16]*mu+1.190784930203603*uvaryr[15]*mu-1.190784930203603*uvaryl[15]*mu-0.3025768239224545*uvarzr[10]*mu+0.3025768239224545*uvarzl[10]*mu+0.4192627457812106*uvarzr[4]*mu+0.4192627457812106*uvarzl[4]*mu-0.8385254915624212*uvarzc[4]*mu);
}
