#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_stream_gen_geo_vol_3x3v_ser_p1(const double *w, const double *dxv, const double *basisVecComp, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // f:         Input distribution function.
  // out:       Incremented output.
  double w3Ddx0  = w[3]/dxv[0]; 
  double dv3Ddx0 = dxv[3]/dxv[0]; 
  double w4Ddx1  = w[4]/dxv[1]; 
  double dv4Ddx1 = dxv[4]/dxv[1]; 
  double w5Ddx2  = w[5]/dxv[2]; 
  double dv5Ddx2 = dxv[5]/dxv[2]; 
  double Gbar[8] = {0.0}; 

  Gbar[0] = 1.224744871391589*f[0]*w3Ddx0+0.3535533905932737*f[4]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[1]*w3Ddx0+0.3535533905932737*f[10]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[2]*w3Ddx0+0.3535533905932737*f[11]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[3]*w3Ddx0+0.3535533905932737*f[12]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[7]*w3Ddx0+0.3535533905932737*f[23]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[8]*w3Ddx0+0.3535533905932737*f[24]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[9]*w3Ddx0+0.3535533905932737*f[25]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[22]*w3Ddx0+0.3535533905932737*f[42]*dv3Ddx0; 
  out[1] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[2]*w3Ddx0+0.3535533905932737*f[11]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[7]*w3Ddx0+0.3535533905932737*f[23]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[0]*w3Ddx0+0.3535533905932737*f[4]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[9]*w3Ddx0+0.3535533905932737*f[25]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[1]*w3Ddx0+0.3535533905932737*f[10]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[22]*w3Ddx0+0.3535533905932737*f[42]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[3]*w3Ddx0+0.3535533905932737*f[12]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[8]*w3Ddx0+0.3535533905932737*f[24]*dv3Ddx0; 
  out[7] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[3]*w3Ddx0+0.3535533905932737*f[12]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[8]*w3Ddx0+0.3535533905932737*f[24]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[9]*w3Ddx0+0.3535533905932737*f[25]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[0]*w3Ddx0+0.3535533905932737*f[4]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[22]*w3Ddx0+0.3535533905932737*f[42]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[1]*w3Ddx0+0.3535533905932737*f[10]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[2]*w3Ddx0+0.3535533905932737*f[11]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[7]*w3Ddx0+0.3535533905932737*f[23]*dv3Ddx0; 
  out[8] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[4]*w3Ddx0+(0.3162277660168379*f[64]+0.3535533905932737*f[0])*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[10]*w3Ddx0+(0.3162277660168379*f[65]+0.3535533905932737*f[1])*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[11]*w3Ddx0+(0.3162277660168379*f[66]+0.3535533905932737*f[2])*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[12]*w3Ddx0+(0.3162277660168379*f[67]+0.3535533905932737*f[3])*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[23]*w3Ddx0+(0.3162277660168379*f[70]+0.3535533905932737*f[7])*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[24]*w3Ddx0+(0.3162277660168379*f[71]+0.3535533905932737*f[8])*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[25]*w3Ddx0+(0.3162277660168379*f[72]+0.3535533905932737*f[9])*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[42]*w3Ddx0+(0.3162277660168379*f[80]+0.3535533905932737*f[22])*dv3Ddx0; 
  out[10] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[5]*w3Ddx0+0.3535533905932737*f[16]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[13]*w3Ddx0+0.3535533905932737*f[29]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[14]*w3Ddx0+0.3535533905932737*f[30]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[15]*w3Ddx0+0.3535533905932737*f[31]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[26]*w3Ddx0+0.3535533905932737*f[44]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[27]*w3Ddx0+0.3535533905932737*f[45]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[28]*w3Ddx0+0.3535533905932737*f[46]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[43]*w3Ddx0+0.3535533905932737*f[57]*dv3Ddx0; 
  out[13] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[6]*w3Ddx0+0.3535533905932737*f[20]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[17]*w3Ddx0+0.3535533905932737*f[35]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[18]*w3Ddx0+0.3535533905932737*f[36]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[19]*w3Ddx0+0.3535533905932737*f[37]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[32]*w3Ddx0+0.3535533905932737*f[48]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[33]*w3Ddx0+0.3535533905932737*f[49]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[34]*w3Ddx0+0.3535533905932737*f[50]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[47]*w3Ddx0+0.3535533905932737*f[58]*dv3Ddx0; 
  out[17] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[9]*w3Ddx0+0.3535533905932737*f[25]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[22]*w3Ddx0+0.3535533905932737*f[42]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[3]*w3Ddx0+0.3535533905932737*f[12]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[2]*w3Ddx0+0.3535533905932737*f[11]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[8]*w3Ddx0+0.3535533905932737*f[24]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[7]*w3Ddx0+0.3535533905932737*f[23]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[0]*w3Ddx0+0.3535533905932737*f[4]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[1]*w3Ddx0+0.3535533905932737*f[10]*dv3Ddx0; 
  out[22] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[11]*w3Ddx0+(0.3162277660168379*f[66]+0.3535533905932737*f[2])*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[23]*w3Ddx0+(0.3162277660168379*f[70]+0.3535533905932737*f[7])*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[4]*w3Ddx0+(0.3162277660168379*f[64]+0.3535533905932737*f[0])*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[25]*w3Ddx0+(0.3162277660168379*f[72]+0.3535533905932737*f[9])*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[10]*w3Ddx0+(0.3162277660168379*f[65]+0.3535533905932737*f[1])*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[42]*w3Ddx0+(0.3162277660168379*f[80]+0.3535533905932737*f[22])*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[12]*w3Ddx0+(0.3162277660168379*f[67]+0.3535533905932737*f[3])*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[24]*w3Ddx0+(0.3162277660168379*f[71]+0.3535533905932737*f[8])*dv3Ddx0; 
  out[23] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[12]*w3Ddx0+(0.3162277660168379*f[67]+0.3535533905932737*f[3])*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[24]*w3Ddx0+(0.3162277660168379*f[71]+0.3535533905932737*f[8])*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[25]*w3Ddx0+(0.3162277660168379*f[72]+0.3535533905932737*f[9])*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[4]*w3Ddx0+(0.3162277660168379*f[64]+0.3535533905932737*f[0])*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[42]*w3Ddx0+(0.3162277660168379*f[80]+0.3535533905932737*f[22])*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[10]*w3Ddx0+(0.3162277660168379*f[65]+0.3535533905932737*f[1])*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[11]*w3Ddx0+(0.3162277660168379*f[66]+0.3535533905932737*f[2])*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[23]*w3Ddx0+(0.3162277660168379*f[70]+0.3535533905932737*f[7])*dv3Ddx0; 
  out[24] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[14]*w3Ddx0+0.3535533905932737*f[30]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[26]*w3Ddx0+0.3535533905932737*f[44]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[5]*w3Ddx0+0.3535533905932737*f[16]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[28]*w3Ddx0+0.3535533905932737*f[46]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[13]*w3Ddx0+0.3535533905932737*f[29]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[43]*w3Ddx0+0.3535533905932737*f[57]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[15]*w3Ddx0+0.3535533905932737*f[31]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[27]*w3Ddx0+0.3535533905932737*f[45]*dv3Ddx0; 
  out[26] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[15]*w3Ddx0+0.3535533905932737*f[31]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[27]*w3Ddx0+0.3535533905932737*f[45]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[28]*w3Ddx0+0.3535533905932737*f[46]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[5]*w3Ddx0+0.3535533905932737*f[16]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[43]*w3Ddx0+0.3535533905932737*f[57]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[13]*w3Ddx0+0.3535533905932737*f[29]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[14]*w3Ddx0+0.3535533905932737*f[30]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[26]*w3Ddx0+0.3535533905932737*f[44]*dv3Ddx0; 
  out[27] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[16]*w3Ddx0+(0.3162277660168379*f[68]+0.3535533905932737*f[5])*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[29]*w3Ddx0+(0.3162277660168379*f[73]+0.3535533905932737*f[13])*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[30]*w3Ddx0+(0.3162277660168379*f[74]+0.3535533905932737*f[14])*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[31]*w3Ddx0+(0.3162277660168379*f[75]+0.3535533905932737*f[15])*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[44]*w3Ddx0+(0.3162277660168379*f[81]+0.3535533905932737*f[26])*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[45]*w3Ddx0+(0.3162277660168379*f[82]+0.3535533905932737*f[27])*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[46]*w3Ddx0+(0.3162277660168379*f[83]+0.3535533905932737*f[28])*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[57]*w3Ddx0+(0.3162277660168379*f[90]+0.3535533905932737*f[43])*dv3Ddx0; 
  out[29] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[18]*w3Ddx0+0.3535533905932737*f[36]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[32]*w3Ddx0+0.3535533905932737*f[48]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[6]*w3Ddx0+0.3535533905932737*f[20]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[34]*w3Ddx0+0.3535533905932737*f[50]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[17]*w3Ddx0+0.3535533905932737*f[35]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[47]*w3Ddx0+0.3535533905932737*f[58]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[19]*w3Ddx0+0.3535533905932737*f[37]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[33]*w3Ddx0+0.3535533905932737*f[49]*dv3Ddx0; 
  out[32] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[19]*w3Ddx0+0.3535533905932737*f[37]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[33]*w3Ddx0+0.3535533905932737*f[49]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[34]*w3Ddx0+0.3535533905932737*f[50]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[6]*w3Ddx0+0.3535533905932737*f[20]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[47]*w3Ddx0+0.3535533905932737*f[58]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[17]*w3Ddx0+0.3535533905932737*f[35]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[18]*w3Ddx0+0.3535533905932737*f[36]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[32]*w3Ddx0+0.3535533905932737*f[48]*dv3Ddx0; 
  out[33] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[20]*w3Ddx0+(0.3162277660168379*f[69]+0.3535533905932737*f[6])*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[35]*w3Ddx0+(0.3162277660168379*f[76]+0.3535533905932737*f[17])*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[36]*w3Ddx0+(0.3162277660168379*f[77]+0.3535533905932737*f[18])*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[37]*w3Ddx0+(0.3162277660168379*f[78]+0.3535533905932737*f[19])*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[48]*w3Ddx0+(0.3162277660168379*f[84]+0.3535533905932737*f[32])*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[49]*w3Ddx0+(0.3162277660168379*f[85]+0.3535533905932737*f[33])*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[50]*w3Ddx0+(0.3162277660168379*f[86]+0.3535533905932737*f[34])*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[58]*w3Ddx0+(0.3162277660168379*f[91]+0.3535533905932737*f[47])*dv3Ddx0; 
  out[35] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[21]*w3Ddx0+0.3535533905932737*f[41]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[38]*w3Ddx0+0.3535533905932737*f[54]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[39]*w3Ddx0+0.3535533905932737*f[55]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[40]*w3Ddx0+0.3535533905932737*f[56]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[51]*w3Ddx0+0.3535533905932737*f[60]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[52]*w3Ddx0+0.3535533905932737*f[61]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[53]*w3Ddx0+0.3535533905932737*f[62]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[59]*w3Ddx0+0.3535533905932737*f[63]*dv3Ddx0; 
  out[38] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[25]*w3Ddx0+(0.3162277660168379*f[72]+0.3535533905932737*f[9])*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[42]*w3Ddx0+(0.3162277660168379*f[80]+0.3535533905932737*f[22])*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[12]*w3Ddx0+(0.3162277660168379*f[67]+0.3535533905932737*f[3])*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[11]*w3Ddx0+(0.3162277660168379*f[66]+0.3535533905932737*f[2])*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[24]*w3Ddx0+(0.3162277660168379*f[71]+0.3535533905932737*f[8])*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[23]*w3Ddx0+(0.3162277660168379*f[70]+0.3535533905932737*f[7])*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[4]*w3Ddx0+(0.3162277660168379*f[64]+0.3535533905932737*f[0])*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[10]*w3Ddx0+(0.3162277660168379*f[65]+0.3535533905932737*f[1])*dv3Ddx0; 
  out[42] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[28]*w3Ddx0+0.3535533905932737*f[46]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[43]*w3Ddx0+0.3535533905932737*f[57]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[15]*w3Ddx0+0.3535533905932737*f[31]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[14]*w3Ddx0+0.3535533905932737*f[30]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[27]*w3Ddx0+0.3535533905932737*f[45]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[26]*w3Ddx0+0.3535533905932737*f[44]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[5]*w3Ddx0+0.3535533905932737*f[16]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[13]*w3Ddx0+0.3535533905932737*f[29]*dv3Ddx0; 
  out[43] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[30]*w3Ddx0+(0.3162277660168379*f[74]+0.3535533905932737*f[14])*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[44]*w3Ddx0+(0.3162277660168379*f[81]+0.3535533905932737*f[26])*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[16]*w3Ddx0+(0.3162277660168379*f[68]+0.3535533905932737*f[5])*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[46]*w3Ddx0+(0.3162277660168379*f[83]+0.3535533905932737*f[28])*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[29]*w3Ddx0+(0.3162277660168379*f[73]+0.3535533905932737*f[13])*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[57]*w3Ddx0+(0.3162277660168379*f[90]+0.3535533905932737*f[43])*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[31]*w3Ddx0+(0.3162277660168379*f[75]+0.3535533905932737*f[15])*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[45]*w3Ddx0+(0.3162277660168379*f[82]+0.3535533905932737*f[27])*dv3Ddx0; 
  out[44] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[31]*w3Ddx0+(0.3162277660168379*f[75]+0.3535533905932737*f[15])*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[45]*w3Ddx0+(0.3162277660168379*f[82]+0.3535533905932737*f[27])*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[46]*w3Ddx0+(0.3162277660168379*f[83]+0.3535533905932737*f[28])*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[16]*w3Ddx0+(0.3162277660168379*f[68]+0.3535533905932737*f[5])*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[57]*w3Ddx0+(0.3162277660168379*f[90]+0.3535533905932737*f[43])*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[29]*w3Ddx0+(0.3162277660168379*f[73]+0.3535533905932737*f[13])*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[30]*w3Ddx0+(0.3162277660168379*f[74]+0.3535533905932737*f[14])*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[44]*w3Ddx0+(0.3162277660168379*f[81]+0.3535533905932737*f[26])*dv3Ddx0; 
  out[45] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[34]*w3Ddx0+0.3535533905932737*f[50]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[47]*w3Ddx0+0.3535533905932737*f[58]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[19]*w3Ddx0+0.3535533905932737*f[37]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[18]*w3Ddx0+0.3535533905932737*f[36]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[33]*w3Ddx0+0.3535533905932737*f[49]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[32]*w3Ddx0+0.3535533905932737*f[48]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[6]*w3Ddx0+0.3535533905932737*f[20]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[17]*w3Ddx0+0.3535533905932737*f[35]*dv3Ddx0; 
  out[47] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[36]*w3Ddx0+(0.3162277660168379*f[77]+0.3535533905932737*f[18])*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[48]*w3Ddx0+(0.3162277660168379*f[84]+0.3535533905932737*f[32])*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[20]*w3Ddx0+(0.3162277660168379*f[69]+0.3535533905932737*f[6])*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[50]*w3Ddx0+(0.3162277660168379*f[86]+0.3535533905932737*f[34])*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[35]*w3Ddx0+(0.3162277660168379*f[76]+0.3535533905932737*f[17])*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[58]*w3Ddx0+(0.3162277660168379*f[91]+0.3535533905932737*f[47])*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[37]*w3Ddx0+(0.3162277660168379*f[78]+0.3535533905932737*f[19])*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[49]*w3Ddx0+(0.3162277660168379*f[85]+0.3535533905932737*f[33])*dv3Ddx0; 
  out[48] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[37]*w3Ddx0+(0.3162277660168379*f[78]+0.3535533905932737*f[19])*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[49]*w3Ddx0+(0.3162277660168379*f[85]+0.3535533905932737*f[33])*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[50]*w3Ddx0+(0.3162277660168379*f[86]+0.3535533905932737*f[34])*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[20]*w3Ddx0+(0.3162277660168379*f[69]+0.3535533905932737*f[6])*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[58]*w3Ddx0+(0.3162277660168379*f[91]+0.3535533905932737*f[47])*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[35]*w3Ddx0+(0.3162277660168379*f[76]+0.3535533905932737*f[17])*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[36]*w3Ddx0+(0.3162277660168379*f[77]+0.3535533905932737*f[18])*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[48]*w3Ddx0+(0.3162277660168379*f[84]+0.3535533905932737*f[32])*dv3Ddx0; 
  out[49] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[39]*w3Ddx0+0.3535533905932737*f[55]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[51]*w3Ddx0+0.3535533905932737*f[60]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[21]*w3Ddx0+0.3535533905932737*f[41]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[53]*w3Ddx0+0.3535533905932737*f[62]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[38]*w3Ddx0+0.3535533905932737*f[54]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[59]*w3Ddx0+0.3535533905932737*f[63]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[40]*w3Ddx0+0.3535533905932737*f[56]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[52]*w3Ddx0+0.3535533905932737*f[61]*dv3Ddx0; 
  out[51] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[40]*w3Ddx0+0.3535533905932737*f[56]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[52]*w3Ddx0+0.3535533905932737*f[61]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[53]*w3Ddx0+0.3535533905932737*f[62]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[21]*w3Ddx0+0.3535533905932737*f[41]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[59]*w3Ddx0+0.3535533905932737*f[63]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[38]*w3Ddx0+0.3535533905932737*f[54]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[39]*w3Ddx0+0.3535533905932737*f[55]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[51]*w3Ddx0+0.3535533905932737*f[60]*dv3Ddx0; 
  out[52] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[41]*w3Ddx0+(0.3162277660168379*f[79]+0.3535533905932737*f[21])*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[54]*w3Ddx0+(0.3162277660168379*f[87]+0.3535533905932737*f[38])*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[55]*w3Ddx0+(0.3162277660168379*f[88]+0.3535533905932737*f[39])*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[56]*w3Ddx0+(0.3162277660168379*f[89]+0.3535533905932737*f[40])*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[60]*w3Ddx0+(0.3162277660168379*f[92]+0.3535533905932737*f[51])*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[61]*w3Ddx0+(0.3162277660168379*f[93]+0.3535533905932737*f[52])*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[62]*w3Ddx0+(0.3162277660168379*f[94]+0.3535533905932737*f[53])*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[63]*w3Ddx0+(0.3162277660168379*f[95]+0.3535533905932737*f[59])*dv3Ddx0; 
  out[54] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[46]*w3Ddx0+(0.3162277660168379*f[83]+0.3535533905932737*f[28])*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[57]*w3Ddx0+(0.3162277660168379*f[90]+0.3535533905932737*f[43])*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[31]*w3Ddx0+(0.3162277660168379*f[75]+0.3535533905932737*f[15])*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[30]*w3Ddx0+(0.3162277660168379*f[74]+0.3535533905932737*f[14])*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[45]*w3Ddx0+(0.3162277660168379*f[82]+0.3535533905932737*f[27])*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[44]*w3Ddx0+(0.3162277660168379*f[81]+0.3535533905932737*f[26])*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[16]*w3Ddx0+(0.3162277660168379*f[68]+0.3535533905932737*f[5])*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[29]*w3Ddx0+(0.3162277660168379*f[73]+0.3535533905932737*f[13])*dv3Ddx0; 
  out[57] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[50]*w3Ddx0+(0.3162277660168379*f[86]+0.3535533905932737*f[34])*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[58]*w3Ddx0+(0.3162277660168379*f[91]+0.3535533905932737*f[47])*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[37]*w3Ddx0+(0.3162277660168379*f[78]+0.3535533905932737*f[19])*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[36]*w3Ddx0+(0.3162277660168379*f[77]+0.3535533905932737*f[18])*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[49]*w3Ddx0+(0.3162277660168379*f[85]+0.3535533905932737*f[33])*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[48]*w3Ddx0+(0.3162277660168379*f[84]+0.3535533905932737*f[32])*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[20]*w3Ddx0+(0.3162277660168379*f[69]+0.3535533905932737*f[6])*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[35]*w3Ddx0+(0.3162277660168379*f[76]+0.3535533905932737*f[17])*dv3Ddx0; 
  out[58] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[53]*w3Ddx0+0.3535533905932737*f[62]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[59]*w3Ddx0+0.3535533905932737*f[63]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[40]*w3Ddx0+0.3535533905932737*f[56]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[39]*w3Ddx0+0.3535533905932737*f[55]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[52]*w3Ddx0+0.3535533905932737*f[61]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[51]*w3Ddx0+0.3535533905932737*f[60]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[21]*w3Ddx0+0.3535533905932737*f[41]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[38]*w3Ddx0+0.3535533905932737*f[54]*dv3Ddx0; 
  out[59] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[55]*w3Ddx0+(0.3162277660168379*f[88]+0.3535533905932737*f[39])*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[60]*w3Ddx0+(0.3162277660168379*f[92]+0.3535533905932737*f[51])*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[41]*w3Ddx0+(0.3162277660168379*f[79]+0.3535533905932737*f[21])*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[62]*w3Ddx0+(0.3162277660168379*f[94]+0.3535533905932737*f[53])*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[54]*w3Ddx0+(0.3162277660168379*f[87]+0.3535533905932737*f[38])*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[63]*w3Ddx0+(0.3162277660168379*f[95]+0.3535533905932737*f[59])*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[56]*w3Ddx0+(0.3162277660168379*f[89]+0.3535533905932737*f[40])*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[61]*w3Ddx0+(0.3162277660168379*f[93]+0.3535533905932737*f[52])*dv3Ddx0; 
  out[60] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[56]*w3Ddx0+(0.3162277660168379*f[89]+0.3535533905932737*f[40])*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[61]*w3Ddx0+(0.3162277660168379*f[93]+0.3535533905932737*f[52])*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[62]*w3Ddx0+(0.3162277660168379*f[94]+0.3535533905932737*f[53])*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[41]*w3Ddx0+(0.3162277660168379*f[79]+0.3535533905932737*f[21])*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[63]*w3Ddx0+(0.3162277660168379*f[95]+0.3535533905932737*f[59])*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[54]*w3Ddx0+(0.3162277660168379*f[87]+0.3535533905932737*f[38])*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[55]*w3Ddx0+(0.3162277660168379*f[88]+0.3535533905932737*f[39])*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[60]*w3Ddx0+(0.3162277660168379*f[92]+0.3535533905932737*f[51])*dv3Ddx0; 
  out[61] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[62]*w3Ddx0+(0.3162277660168379*f[94]+0.3535533905932737*f[53])*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[63]*w3Ddx0+(0.3162277660168379*f[95]+0.3535533905932737*f[59])*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[56]*w3Ddx0+(0.3162277660168379*f[89]+0.3535533905932737*f[40])*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[55]*w3Ddx0+(0.3162277660168379*f[88]+0.3535533905932737*f[39])*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[61]*w3Ddx0+(0.3162277660168379*f[93]+0.3535533905932737*f[52])*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[60]*w3Ddx0+(0.3162277660168379*f[92]+0.3535533905932737*f[51])*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[41]*w3Ddx0+(0.3162277660168379*f[79]+0.3535533905932737*f[21])*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[54]*w3Ddx0+(0.3162277660168379*f[87]+0.3535533905932737*f[38])*dv3Ddx0; 
  out[63] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[64]*w3Ddx0+0.3162277660168379*f[4]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[65]*w3Ddx0+0.3162277660168379*f[10]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[66]*w3Ddx0+0.3162277660168379*f[11]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[67]*w3Ddx0+0.3162277660168379*f[12]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[70]*w3Ddx0+0.3162277660168379*f[23]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[71]*w3Ddx0+0.3162277660168379*f[24]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[72]*w3Ddx0+0.3162277660168379*f[25]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[80]*w3Ddx0+0.3162277660168379*f[42]*dv3Ddx0; 
  out[65] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[66]*w3Ddx0+0.3162277660168379*f[11]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[70]*w3Ddx0+0.3162277660168379*f[23]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[64]*w3Ddx0+0.3162277660168379*f[4]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[72]*w3Ddx0+0.3162277660168379*f[25]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[65]*w3Ddx0+0.3162277660168379*f[10]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[80]*w3Ddx0+0.3162277660168379*f[42]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[67]*w3Ddx0+0.3162277660168379*f[12]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[71]*w3Ddx0+0.3162277660168379*f[24]*dv3Ddx0; 
  out[70] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[67]*w3Ddx0+0.3162277660168379*f[12]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[71]*w3Ddx0+0.3162277660168379*f[24]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[72]*w3Ddx0+0.3162277660168379*f[25]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[64]*w3Ddx0+0.3162277660168379*f[4]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[80]*w3Ddx0+0.3162277660168379*f[42]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[65]*w3Ddx0+0.3162277660168379*f[10]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[66]*w3Ddx0+0.3162277660168379*f[11]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[70]*w3Ddx0+0.3162277660168379*f[23]*dv3Ddx0; 
  out[71] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[68]*w3Ddx0+0.3162277660168379*f[16]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[73]*w3Ddx0+0.3162277660168379*f[29]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[74]*w3Ddx0+0.3162277660168379*f[30]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[75]*w3Ddx0+0.3162277660168379*f[31]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[81]*w3Ddx0+0.3162277660168379*f[44]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[82]*w3Ddx0+0.3162277660168379*f[45]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[83]*w3Ddx0+0.3162277660168379*f[46]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[90]*w3Ddx0+0.3162277660168379*f[57]*dv3Ddx0; 
  out[73] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[69]*w3Ddx0+0.3162277660168379*f[20]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[76]*w3Ddx0+0.3162277660168379*f[35]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[77]*w3Ddx0+0.3162277660168379*f[36]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[78]*w3Ddx0+0.3162277660168379*f[37]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[84]*w3Ddx0+0.3162277660168379*f[48]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[85]*w3Ddx0+0.3162277660168379*f[49]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[86]*w3Ddx0+0.3162277660168379*f[50]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[91]*w3Ddx0+0.3162277660168379*f[58]*dv3Ddx0; 
  out[76] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[72]*w3Ddx0+0.3162277660168379*f[25]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[80]*w3Ddx0+0.3162277660168379*f[42]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[67]*w3Ddx0+0.3162277660168379*f[12]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[66]*w3Ddx0+0.3162277660168379*f[11]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[71]*w3Ddx0+0.3162277660168379*f[24]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[70]*w3Ddx0+0.3162277660168379*f[23]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[64]*w3Ddx0+0.3162277660168379*f[4]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[65]*w3Ddx0+0.3162277660168379*f[10]*dv3Ddx0; 
  out[80] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[74]*w3Ddx0+0.3162277660168379*f[30]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[81]*w3Ddx0+0.3162277660168379*f[44]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[68]*w3Ddx0+0.3162277660168379*f[16]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[83]*w3Ddx0+0.3162277660168379*f[46]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[73]*w3Ddx0+0.3162277660168379*f[29]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[90]*w3Ddx0+0.3162277660168379*f[57]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[75]*w3Ddx0+0.3162277660168379*f[31]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[82]*w3Ddx0+0.3162277660168379*f[45]*dv3Ddx0; 
  out[81] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[75]*w3Ddx0+0.3162277660168379*f[31]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[82]*w3Ddx0+0.3162277660168379*f[45]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[83]*w3Ddx0+0.3162277660168379*f[46]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[68]*w3Ddx0+0.3162277660168379*f[16]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[90]*w3Ddx0+0.3162277660168379*f[57]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[73]*w3Ddx0+0.3162277660168379*f[29]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[74]*w3Ddx0+0.3162277660168379*f[30]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[81]*w3Ddx0+0.3162277660168379*f[44]*dv3Ddx0; 
  out[82] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[77]*w3Ddx0+0.3162277660168379*f[36]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[84]*w3Ddx0+0.3162277660168379*f[48]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[69]*w3Ddx0+0.3162277660168379*f[20]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[86]*w3Ddx0+0.3162277660168379*f[50]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[76]*w3Ddx0+0.3162277660168379*f[35]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[91]*w3Ddx0+0.3162277660168379*f[58]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[78]*w3Ddx0+0.3162277660168379*f[37]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[85]*w3Ddx0+0.3162277660168379*f[49]*dv3Ddx0; 
  out[84] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[78]*w3Ddx0+0.3162277660168379*f[37]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[85]*w3Ddx0+0.3162277660168379*f[49]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[86]*w3Ddx0+0.3162277660168379*f[50]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[69]*w3Ddx0+0.3162277660168379*f[20]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[91]*w3Ddx0+0.3162277660168379*f[58]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[76]*w3Ddx0+0.3162277660168379*f[35]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[77]*w3Ddx0+0.3162277660168379*f[36]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[84]*w3Ddx0+0.3162277660168379*f[48]*dv3Ddx0; 
  out[85] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[79]*w3Ddx0+0.3162277660168379*f[41]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[87]*w3Ddx0+0.3162277660168379*f[54]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[88]*w3Ddx0+0.3162277660168379*f[55]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[89]*w3Ddx0+0.3162277660168379*f[56]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[92]*w3Ddx0+0.3162277660168379*f[60]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[93]*w3Ddx0+0.3162277660168379*f[61]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[94]*w3Ddx0+0.3162277660168379*f[62]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[95]*w3Ddx0+0.3162277660168379*f[63]*dv3Ddx0; 
  out[87] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[83]*w3Ddx0+0.3162277660168379*f[46]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[90]*w3Ddx0+0.3162277660168379*f[57]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[75]*w3Ddx0+0.3162277660168379*f[31]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[74]*w3Ddx0+0.3162277660168379*f[30]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[82]*w3Ddx0+0.3162277660168379*f[45]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[81]*w3Ddx0+0.3162277660168379*f[44]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[68]*w3Ddx0+0.3162277660168379*f[16]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[73]*w3Ddx0+0.3162277660168379*f[29]*dv3Ddx0; 
  out[90] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[86]*w3Ddx0+0.3162277660168379*f[50]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[91]*w3Ddx0+0.3162277660168379*f[58]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[78]*w3Ddx0+0.3162277660168379*f[37]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[77]*w3Ddx0+0.3162277660168379*f[36]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[85]*w3Ddx0+0.3162277660168379*f[49]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[84]*w3Ddx0+0.3162277660168379*f[48]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[69]*w3Ddx0+0.3162277660168379*f[20]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[76]*w3Ddx0+0.3162277660168379*f[35]*dv3Ddx0; 
  out[91] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[88]*w3Ddx0+0.3162277660168379*f[55]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[92]*w3Ddx0+0.3162277660168379*f[60]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[79]*w3Ddx0+0.3162277660168379*f[41]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[94]*w3Ddx0+0.3162277660168379*f[62]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[87]*w3Ddx0+0.3162277660168379*f[54]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[95]*w3Ddx0+0.3162277660168379*f[63]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[89]*w3Ddx0+0.3162277660168379*f[56]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[93]*w3Ddx0+0.3162277660168379*f[61]*dv3Ddx0; 
  out[92] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[89]*w3Ddx0+0.3162277660168379*f[56]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[93]*w3Ddx0+0.3162277660168379*f[61]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[94]*w3Ddx0+0.3162277660168379*f[62]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[79]*w3Ddx0+0.3162277660168379*f[41]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[95]*w3Ddx0+0.3162277660168379*f[63]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[87]*w3Ddx0+0.3162277660168379*f[54]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[88]*w3Ddx0+0.3162277660168379*f[55]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[92]*w3Ddx0+0.3162277660168379*f[60]*dv3Ddx0; 
  out[93] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[94]*w3Ddx0+0.3162277660168379*f[62]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[95]*w3Ddx0+0.3162277660168379*f[63]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[89]*w3Ddx0+0.3162277660168379*f[56]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[88]*w3Ddx0+0.3162277660168379*f[55]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[93]*w3Ddx0+0.3162277660168379*f[61]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[92]*w3Ddx0+0.3162277660168379*f[60]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[79]*w3Ddx0+0.3162277660168379*f[41]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[87]*w3Ddx0+0.3162277660168379*f[54]*dv3Ddx0; 
  out[95] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[96]*w3Ddx0+0.3535533905932737*f[100]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[97]*w3Ddx0+0.3535533905932737*f[105]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[98]*w3Ddx0+0.3535533905932737*f[106]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[99]*w3Ddx0+0.3535533905932737*f[107]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[102]*w3Ddx0+0.3535533905932737*f[113]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[103]*w3Ddx0+0.3535533905932737*f[114]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[104]*w3Ddx0+0.3535533905932737*f[115]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[112]*w3Ddx0+0.3535533905932737*f[122]*dv3Ddx0; 
  out[97] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[98]*w3Ddx0+0.3535533905932737*f[106]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[102]*w3Ddx0+0.3535533905932737*f[113]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[96]*w3Ddx0+0.3535533905932737*f[100]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[104]*w3Ddx0+0.3535533905932737*f[115]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[97]*w3Ddx0+0.3535533905932737*f[105]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[112]*w3Ddx0+0.3535533905932737*f[122]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[99]*w3Ddx0+0.3535533905932737*f[107]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[103]*w3Ddx0+0.3535533905932737*f[114]*dv3Ddx0; 
  out[102] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[99]*w3Ddx0+0.3535533905932737*f[107]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[103]*w3Ddx0+0.3535533905932737*f[114]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[104]*w3Ddx0+0.3535533905932737*f[115]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[96]*w3Ddx0+0.3535533905932737*f[100]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[112]*w3Ddx0+0.3535533905932737*f[122]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[97]*w3Ddx0+0.3535533905932737*f[105]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[98]*w3Ddx0+0.3535533905932737*f[106]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[102]*w3Ddx0+0.3535533905932737*f[113]*dv3Ddx0; 
  out[103] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[100]*w3Ddx0+0.3535533905932737*f[96]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[105]*w3Ddx0+0.3535533905932737*f[97]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[106]*w3Ddx0+0.3535533905932737*f[98]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[107]*w3Ddx0+0.3535533905932737*f[99]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[113]*w3Ddx0+0.3535533905932737*f[102]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[114]*w3Ddx0+0.3535533905932737*f[103]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[115]*w3Ddx0+0.3535533905932737*f[104]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[122]*w3Ddx0+0.3535533905932737*f[112]*dv3Ddx0; 
  out[105] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[101]*w3Ddx0+0.3535533905932737*f[111]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[108]*w3Ddx0+0.3535533905932737*f[119]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[109]*w3Ddx0+0.3535533905932737*f[120]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[110]*w3Ddx0+0.3535533905932737*f[121]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[116]*w3Ddx0+0.3535533905932737*f[124]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[117]*w3Ddx0+0.3535533905932737*f[125]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[118]*w3Ddx0+0.3535533905932737*f[126]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[123]*w3Ddx0+0.3535533905932737*f[127]*dv3Ddx0; 
  out[108] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[104]*w3Ddx0+0.3535533905932737*f[115]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[112]*w3Ddx0+0.3535533905932737*f[122]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[99]*w3Ddx0+0.3535533905932737*f[107]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[98]*w3Ddx0+0.3535533905932737*f[106]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[103]*w3Ddx0+0.3535533905932737*f[114]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[102]*w3Ddx0+0.3535533905932737*f[113]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[96]*w3Ddx0+0.3535533905932737*f[100]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[97]*w3Ddx0+0.3535533905932737*f[105]*dv3Ddx0; 
  out[112] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[106]*w3Ddx0+0.3535533905932737*f[98]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[113]*w3Ddx0+0.3535533905932737*f[102]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[100]*w3Ddx0+0.3535533905932737*f[96]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[115]*w3Ddx0+0.3535533905932737*f[104]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[105]*w3Ddx0+0.3535533905932737*f[97]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[122]*w3Ddx0+0.3535533905932737*f[112]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[107]*w3Ddx0+0.3535533905932737*f[99]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[114]*w3Ddx0+0.3535533905932737*f[103]*dv3Ddx0; 
  out[113] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[107]*w3Ddx0+0.3535533905932737*f[99]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[114]*w3Ddx0+0.3535533905932737*f[103]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[115]*w3Ddx0+0.3535533905932737*f[104]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[100]*w3Ddx0+0.3535533905932737*f[96]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[122]*w3Ddx0+0.3535533905932737*f[112]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[105]*w3Ddx0+0.3535533905932737*f[97]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[106]*w3Ddx0+0.3535533905932737*f[98]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[113]*w3Ddx0+0.3535533905932737*f[102]*dv3Ddx0; 
  out[114] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[109]*w3Ddx0+0.3535533905932737*f[120]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[116]*w3Ddx0+0.3535533905932737*f[124]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[101]*w3Ddx0+0.3535533905932737*f[111]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[118]*w3Ddx0+0.3535533905932737*f[126]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[108]*w3Ddx0+0.3535533905932737*f[119]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[123]*w3Ddx0+0.3535533905932737*f[127]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[110]*w3Ddx0+0.3535533905932737*f[121]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[117]*w3Ddx0+0.3535533905932737*f[125]*dv3Ddx0; 
  out[116] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[110]*w3Ddx0+0.3535533905932737*f[121]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[117]*w3Ddx0+0.3535533905932737*f[125]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[118]*w3Ddx0+0.3535533905932737*f[126]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[101]*w3Ddx0+0.3535533905932737*f[111]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[123]*w3Ddx0+0.3535533905932737*f[127]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[108]*w3Ddx0+0.3535533905932737*f[119]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[109]*w3Ddx0+0.3535533905932737*f[120]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[116]*w3Ddx0+0.3535533905932737*f[124]*dv3Ddx0; 
  out[117] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[111]*w3Ddx0+0.3535533905932737*f[101]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[119]*w3Ddx0+0.3535533905932737*f[108]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[120]*w3Ddx0+0.3535533905932737*f[109]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[121]*w3Ddx0+0.3535533905932737*f[110]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[124]*w3Ddx0+0.3535533905932737*f[116]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[125]*w3Ddx0+0.3535533905932737*f[117]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[126]*w3Ddx0+0.3535533905932737*f[118]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[127]*w3Ddx0+0.3535533905932737*f[123]*dv3Ddx0; 
  out[119] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[115]*w3Ddx0+0.3535533905932737*f[104]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[122]*w3Ddx0+0.3535533905932737*f[112]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[107]*w3Ddx0+0.3535533905932737*f[99]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[106]*w3Ddx0+0.3535533905932737*f[98]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[114]*w3Ddx0+0.3535533905932737*f[103]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[113]*w3Ddx0+0.3535533905932737*f[102]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[100]*w3Ddx0+0.3535533905932737*f[96]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[105]*w3Ddx0+0.3535533905932737*f[97]*dv3Ddx0; 
  out[122] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[118]*w3Ddx0+0.3535533905932737*f[126]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[123]*w3Ddx0+0.3535533905932737*f[127]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[110]*w3Ddx0+0.3535533905932737*f[121]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[109]*w3Ddx0+0.3535533905932737*f[120]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[117]*w3Ddx0+0.3535533905932737*f[125]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[116]*w3Ddx0+0.3535533905932737*f[124]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[101]*w3Ddx0+0.3535533905932737*f[111]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[108]*w3Ddx0+0.3535533905932737*f[119]*dv3Ddx0; 
  out[123] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[120]*w3Ddx0+0.3535533905932737*f[109]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[124]*w3Ddx0+0.3535533905932737*f[116]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[111]*w3Ddx0+0.3535533905932737*f[101]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[126]*w3Ddx0+0.3535533905932737*f[118]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[119]*w3Ddx0+0.3535533905932737*f[108]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[127]*w3Ddx0+0.3535533905932737*f[123]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[121]*w3Ddx0+0.3535533905932737*f[110]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[125]*w3Ddx0+0.3535533905932737*f[117]*dv3Ddx0; 
  out[124] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[121]*w3Ddx0+0.3535533905932737*f[110]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[125]*w3Ddx0+0.3535533905932737*f[117]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[126]*w3Ddx0+0.3535533905932737*f[118]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[111]*w3Ddx0+0.3535533905932737*f[101]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[127]*w3Ddx0+0.3535533905932737*f[123]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[119]*w3Ddx0+0.3535533905932737*f[108]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[120]*w3Ddx0+0.3535533905932737*f[109]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[124]*w3Ddx0+0.3535533905932737*f[116]*dv3Ddx0; 
  out[125] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[126]*w3Ddx0+0.3535533905932737*f[118]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[127]*w3Ddx0+0.3535533905932737*f[123]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[121]*w3Ddx0+0.3535533905932737*f[110]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[120]*w3Ddx0+0.3535533905932737*f[109]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[125]*w3Ddx0+0.3535533905932737*f[117]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[124]*w3Ddx0+0.3535533905932737*f[116]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[111]*w3Ddx0+0.3535533905932737*f[101]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[119]*w3Ddx0+0.3535533905932737*f[108]*dv3Ddx0; 
  out[127] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[128]*w3Ddx0+0.3535533905932737*f[132]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[129]*w3Ddx0+0.3535533905932737*f[137]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[130]*w3Ddx0+0.3535533905932737*f[138]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[131]*w3Ddx0+0.3535533905932737*f[139]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[134]*w3Ddx0+0.3535533905932737*f[145]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[135]*w3Ddx0+0.3535533905932737*f[146]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[136]*w3Ddx0+0.3535533905932737*f[147]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[144]*w3Ddx0+0.3535533905932737*f[154]*dv3Ddx0; 
  out[129] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[130]*w3Ddx0+0.3535533905932737*f[138]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[134]*w3Ddx0+0.3535533905932737*f[145]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[128]*w3Ddx0+0.3535533905932737*f[132]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[136]*w3Ddx0+0.3535533905932737*f[147]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[129]*w3Ddx0+0.3535533905932737*f[137]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[144]*w3Ddx0+0.3535533905932737*f[154]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[131]*w3Ddx0+0.3535533905932737*f[139]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[135]*w3Ddx0+0.3535533905932737*f[146]*dv3Ddx0; 
  out[134] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[131]*w3Ddx0+0.3535533905932737*f[139]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[135]*w3Ddx0+0.3535533905932737*f[146]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[136]*w3Ddx0+0.3535533905932737*f[147]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[128]*w3Ddx0+0.3535533905932737*f[132]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[144]*w3Ddx0+0.3535533905932737*f[154]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[129]*w3Ddx0+0.3535533905932737*f[137]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[130]*w3Ddx0+0.3535533905932737*f[138]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[134]*w3Ddx0+0.3535533905932737*f[145]*dv3Ddx0; 
  out[135] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[132]*w3Ddx0+0.3535533905932737*f[128]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[137]*w3Ddx0+0.3535533905932737*f[129]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[138]*w3Ddx0+0.3535533905932737*f[130]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[139]*w3Ddx0+0.3535533905932737*f[131]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[145]*w3Ddx0+0.3535533905932737*f[134]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[146]*w3Ddx0+0.3535533905932737*f[135]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[147]*w3Ddx0+0.3535533905932737*f[136]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[154]*w3Ddx0+0.3535533905932737*f[144]*dv3Ddx0; 
  out[137] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[133]*w3Ddx0+0.3535533905932737*f[143]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[140]*w3Ddx0+0.3535533905932737*f[151]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[141]*w3Ddx0+0.3535533905932737*f[152]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[142]*w3Ddx0+0.3535533905932737*f[153]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[148]*w3Ddx0+0.3535533905932737*f[156]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[149]*w3Ddx0+0.3535533905932737*f[157]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[150]*w3Ddx0+0.3535533905932737*f[158]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[155]*w3Ddx0+0.3535533905932737*f[159]*dv3Ddx0; 
  out[140] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[136]*w3Ddx0+0.3535533905932737*f[147]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[144]*w3Ddx0+0.3535533905932737*f[154]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[131]*w3Ddx0+0.3535533905932737*f[139]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[130]*w3Ddx0+0.3535533905932737*f[138]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[135]*w3Ddx0+0.3535533905932737*f[146]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[134]*w3Ddx0+0.3535533905932737*f[145]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[128]*w3Ddx0+0.3535533905932737*f[132]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[129]*w3Ddx0+0.3535533905932737*f[137]*dv3Ddx0; 
  out[144] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[138]*w3Ddx0+0.3535533905932737*f[130]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[145]*w3Ddx0+0.3535533905932737*f[134]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[132]*w3Ddx0+0.3535533905932737*f[128]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[147]*w3Ddx0+0.3535533905932737*f[136]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[137]*w3Ddx0+0.3535533905932737*f[129]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[154]*w3Ddx0+0.3535533905932737*f[144]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[139]*w3Ddx0+0.3535533905932737*f[131]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[146]*w3Ddx0+0.3535533905932737*f[135]*dv3Ddx0; 
  out[145] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[139]*w3Ddx0+0.3535533905932737*f[131]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[146]*w3Ddx0+0.3535533905932737*f[135]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[147]*w3Ddx0+0.3535533905932737*f[136]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[132]*w3Ddx0+0.3535533905932737*f[128]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[154]*w3Ddx0+0.3535533905932737*f[144]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[137]*w3Ddx0+0.3535533905932737*f[129]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[138]*w3Ddx0+0.3535533905932737*f[130]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[145]*w3Ddx0+0.3535533905932737*f[134]*dv3Ddx0; 
  out[146] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[141]*w3Ddx0+0.3535533905932737*f[152]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[148]*w3Ddx0+0.3535533905932737*f[156]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[133]*w3Ddx0+0.3535533905932737*f[143]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[150]*w3Ddx0+0.3535533905932737*f[158]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[140]*w3Ddx0+0.3535533905932737*f[151]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[155]*w3Ddx0+0.3535533905932737*f[159]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[142]*w3Ddx0+0.3535533905932737*f[153]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[149]*w3Ddx0+0.3535533905932737*f[157]*dv3Ddx0; 
  out[148] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[142]*w3Ddx0+0.3535533905932737*f[153]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[149]*w3Ddx0+0.3535533905932737*f[157]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[150]*w3Ddx0+0.3535533905932737*f[158]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[133]*w3Ddx0+0.3535533905932737*f[143]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[155]*w3Ddx0+0.3535533905932737*f[159]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[140]*w3Ddx0+0.3535533905932737*f[151]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[141]*w3Ddx0+0.3535533905932737*f[152]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[148]*w3Ddx0+0.3535533905932737*f[156]*dv3Ddx0; 
  out[149] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[143]*w3Ddx0+0.3535533905932737*f[133]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[151]*w3Ddx0+0.3535533905932737*f[140]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[152]*w3Ddx0+0.3535533905932737*f[141]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[153]*w3Ddx0+0.3535533905932737*f[142]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[156]*w3Ddx0+0.3535533905932737*f[148]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[157]*w3Ddx0+0.3535533905932737*f[149]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[158]*w3Ddx0+0.3535533905932737*f[150]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[159]*w3Ddx0+0.3535533905932737*f[155]*dv3Ddx0; 
  out[151] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[147]*w3Ddx0+0.3535533905932737*f[136]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[154]*w3Ddx0+0.3535533905932737*f[144]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[139]*w3Ddx0+0.3535533905932737*f[131]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[138]*w3Ddx0+0.3535533905932737*f[130]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[146]*w3Ddx0+0.3535533905932737*f[135]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[145]*w3Ddx0+0.3535533905932737*f[134]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[132]*w3Ddx0+0.3535533905932737*f[128]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[137]*w3Ddx0+0.3535533905932737*f[129]*dv3Ddx0; 
  out[154] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[150]*w3Ddx0+0.3535533905932737*f[158]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[155]*w3Ddx0+0.3535533905932737*f[159]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[142]*w3Ddx0+0.3535533905932737*f[153]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[141]*w3Ddx0+0.3535533905932737*f[152]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[149]*w3Ddx0+0.3535533905932737*f[157]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[148]*w3Ddx0+0.3535533905932737*f[156]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[133]*w3Ddx0+0.3535533905932737*f[143]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[140]*w3Ddx0+0.3535533905932737*f[151]*dv3Ddx0; 
  out[155] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[152]*w3Ddx0+0.3535533905932737*f[141]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[156]*w3Ddx0+0.3535533905932737*f[148]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[143]*w3Ddx0+0.3535533905932737*f[133]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[158]*w3Ddx0+0.3535533905932737*f[150]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[151]*w3Ddx0+0.3535533905932737*f[140]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[159]*w3Ddx0+0.3535533905932737*f[155]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[153]*w3Ddx0+0.3535533905932737*f[142]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[157]*w3Ddx0+0.3535533905932737*f[149]*dv3Ddx0; 
  out[156] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[153]*w3Ddx0+0.3535533905932737*f[142]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[157]*w3Ddx0+0.3535533905932737*f[149]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[158]*w3Ddx0+0.3535533905932737*f[150]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[143]*w3Ddx0+0.3535533905932737*f[133]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[159]*w3Ddx0+0.3535533905932737*f[155]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[151]*w3Ddx0+0.3535533905932737*f[140]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[152]*w3Ddx0+0.3535533905932737*f[141]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[156]*w3Ddx0+0.3535533905932737*f[148]*dv3Ddx0; 
  out[157] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[158]*w3Ddx0+0.3535533905932737*f[150]*dv3Ddx0; 
  Gbar[1] = 1.224744871391589*f[159]*w3Ddx0+0.3535533905932737*f[155]*dv3Ddx0; 
  Gbar[2] = 1.224744871391589*f[153]*w3Ddx0+0.3535533905932737*f[142]*dv3Ddx0; 
  Gbar[3] = 1.224744871391589*f[152]*w3Ddx0+0.3535533905932737*f[141]*dv3Ddx0; 
  Gbar[4] = 1.224744871391589*f[157]*w3Ddx0+0.3535533905932737*f[149]*dv3Ddx0; 
  Gbar[5] = 1.224744871391589*f[156]*w3Ddx0+0.3535533905932737*f[148]*dv3Ddx0; 
  Gbar[6] = 1.224744871391589*f[143]*w3Ddx0+0.3535533905932737*f[133]*dv3Ddx0; 
  Gbar[7] = 1.224744871391589*f[151]*w3Ddx0+0.3535533905932737*f[140]*dv3Ddx0; 
  out[159] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[0]*w4Ddx1+0.3535533905932737*f[5]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[1]*w4Ddx1+0.3535533905932737*f[13]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[2]*w4Ddx1+0.3535533905932737*f[14]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[3]*w4Ddx1+0.3535533905932737*f[15]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[7]*w4Ddx1+0.3535533905932737*f[26]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[8]*w4Ddx1+0.3535533905932737*f[27]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[9]*w4Ddx1+0.3535533905932737*f[28]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[22]*w4Ddx1+0.3535533905932737*f[43]*dv4Ddx1; 
  out[2] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[1]*w4Ddx1+0.3535533905932737*f[13]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[0]*w4Ddx1+0.3535533905932737*f[5]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[7]*w4Ddx1+0.3535533905932737*f[26]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[8]*w4Ddx1+0.3535533905932737*f[27]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[2]*w4Ddx1+0.3535533905932737*f[14]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[3]*w4Ddx1+0.3535533905932737*f[15]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[22]*w4Ddx1+0.3535533905932737*f[43]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[9]*w4Ddx1+0.3535533905932737*f[28]*dv4Ddx1; 
  out[7] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[3]*w4Ddx1+0.3535533905932737*f[15]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[8]*w4Ddx1+0.3535533905932737*f[27]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[9]*w4Ddx1+0.3535533905932737*f[28]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[0]*w4Ddx1+0.3535533905932737*f[5]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[22]*w4Ddx1+0.3535533905932737*f[43]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[1]*w4Ddx1+0.3535533905932737*f[13]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[2]*w4Ddx1+0.3535533905932737*f[14]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[7]*w4Ddx1+0.3535533905932737*f[26]*dv4Ddx1; 
  out[9] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[4]*w4Ddx1+0.3535533905932737*f[16]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[10]*w4Ddx1+0.3535533905932737*f[29]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[11]*w4Ddx1+0.3535533905932737*f[30]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[12]*w4Ddx1+0.3535533905932737*f[31]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[23]*w4Ddx1+0.3535533905932737*f[44]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[24]*w4Ddx1+0.3535533905932737*f[45]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[25]*w4Ddx1+0.3535533905932737*f[46]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[42]*w4Ddx1+0.3535533905932737*f[57]*dv4Ddx1; 
  out[11] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[5]*w4Ddx1+(0.3162277660168379*f[96]+0.3535533905932737*f[0])*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[13]*w4Ddx1+(0.3162277660168379*f[97]+0.3535533905932737*f[1])*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[14]*w4Ddx1+(0.3162277660168379*f[98]+0.3535533905932737*f[2])*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[15]*w4Ddx1+(0.3162277660168379*f[99]+0.3535533905932737*f[3])*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[26]*w4Ddx1+(0.3162277660168379*f[102]+0.3535533905932737*f[7])*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[27]*w4Ddx1+(0.3162277660168379*f[103]+0.3535533905932737*f[8])*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[28]*w4Ddx1+(0.3162277660168379*f[104]+0.3535533905932737*f[9])*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[43]*w4Ddx1+(0.3162277660168379*f[112]+0.3535533905932737*f[22])*dv4Ddx1; 
  out[14] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[6]*w4Ddx1+0.3535533905932737*f[21]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[17]*w4Ddx1+0.3535533905932737*f[38]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[18]*w4Ddx1+0.3535533905932737*f[39]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[19]*w4Ddx1+0.3535533905932737*f[40]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[32]*w4Ddx1+0.3535533905932737*f[51]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[33]*w4Ddx1+0.3535533905932737*f[52]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[34]*w4Ddx1+0.3535533905932737*f[53]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[47]*w4Ddx1+0.3535533905932737*f[59]*dv4Ddx1; 
  out[18] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[8]*w4Ddx1+0.3535533905932737*f[27]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[3]*w4Ddx1+0.3535533905932737*f[15]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[22]*w4Ddx1+0.3535533905932737*f[43]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[1]*w4Ddx1+0.3535533905932737*f[13]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[9]*w4Ddx1+0.3535533905932737*f[28]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[0]*w4Ddx1+0.3535533905932737*f[5]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[7]*w4Ddx1+0.3535533905932737*f[26]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[2]*w4Ddx1+0.3535533905932737*f[14]*dv4Ddx1; 
  out[22] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[10]*w4Ddx1+0.3535533905932737*f[29]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[4]*w4Ddx1+0.3535533905932737*f[16]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[23]*w4Ddx1+0.3535533905932737*f[44]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[24]*w4Ddx1+0.3535533905932737*f[45]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[11]*w4Ddx1+0.3535533905932737*f[30]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[12]*w4Ddx1+0.3535533905932737*f[31]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[42]*w4Ddx1+0.3535533905932737*f[57]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[25]*w4Ddx1+0.3535533905932737*f[46]*dv4Ddx1; 
  out[23] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[12]*w4Ddx1+0.3535533905932737*f[31]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[24]*w4Ddx1+0.3535533905932737*f[45]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[25]*w4Ddx1+0.3535533905932737*f[46]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[4]*w4Ddx1+0.3535533905932737*f[16]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[42]*w4Ddx1+0.3535533905932737*f[57]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[10]*w4Ddx1+0.3535533905932737*f[29]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[11]*w4Ddx1+0.3535533905932737*f[30]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[23]*w4Ddx1+0.3535533905932737*f[44]*dv4Ddx1; 
  out[25] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[13]*w4Ddx1+(0.3162277660168379*f[97]+0.3535533905932737*f[1])*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[5]*w4Ddx1+(0.3162277660168379*f[96]+0.3535533905932737*f[0])*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[26]*w4Ddx1+(0.3162277660168379*f[102]+0.3535533905932737*f[7])*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[27]*w4Ddx1+(0.3162277660168379*f[103]+0.3535533905932737*f[8])*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[14]*w4Ddx1+(0.3162277660168379*f[98]+0.3535533905932737*f[2])*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[15]*w4Ddx1+(0.3162277660168379*f[99]+0.3535533905932737*f[3])*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[43]*w4Ddx1+(0.3162277660168379*f[112]+0.3535533905932737*f[22])*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[28]*w4Ddx1+(0.3162277660168379*f[104]+0.3535533905932737*f[9])*dv4Ddx1; 
  out[26] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[15]*w4Ddx1+(0.3162277660168379*f[99]+0.3535533905932737*f[3])*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[27]*w4Ddx1+(0.3162277660168379*f[103]+0.3535533905932737*f[8])*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[28]*w4Ddx1+(0.3162277660168379*f[104]+0.3535533905932737*f[9])*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[5]*w4Ddx1+(0.3162277660168379*f[96]+0.3535533905932737*f[0])*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[43]*w4Ddx1+(0.3162277660168379*f[112]+0.3535533905932737*f[22])*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[13]*w4Ddx1+(0.3162277660168379*f[97]+0.3535533905932737*f[1])*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[14]*w4Ddx1+(0.3162277660168379*f[98]+0.3535533905932737*f[2])*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[26]*w4Ddx1+(0.3162277660168379*f[102]+0.3535533905932737*f[7])*dv4Ddx1; 
  out[28] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[16]*w4Ddx1+(0.3162277660168379*f[100]+0.3535533905932737*f[4])*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[29]*w4Ddx1+(0.3162277660168379*f[105]+0.3535533905932737*f[10])*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[30]*w4Ddx1+(0.3162277660168379*f[106]+0.3535533905932737*f[11])*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[31]*w4Ddx1+(0.3162277660168379*f[107]+0.3535533905932737*f[12])*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[44]*w4Ddx1+(0.3162277660168379*f[113]+0.3535533905932737*f[23])*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[45]*w4Ddx1+(0.3162277660168379*f[114]+0.3535533905932737*f[24])*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[46]*w4Ddx1+(0.3162277660168379*f[115]+0.3535533905932737*f[25])*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[57]*w4Ddx1+(0.3162277660168379*f[122]+0.3535533905932737*f[42])*dv4Ddx1; 
  out[30] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[17]*w4Ddx1+0.3535533905932737*f[38]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[6]*w4Ddx1+0.3535533905932737*f[21]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[32]*w4Ddx1+0.3535533905932737*f[51]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[33]*w4Ddx1+0.3535533905932737*f[52]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[18]*w4Ddx1+0.3535533905932737*f[39]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[19]*w4Ddx1+0.3535533905932737*f[40]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[47]*w4Ddx1+0.3535533905932737*f[59]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[34]*w4Ddx1+0.3535533905932737*f[53]*dv4Ddx1; 
  out[32] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[19]*w4Ddx1+0.3535533905932737*f[40]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[33]*w4Ddx1+0.3535533905932737*f[52]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[34]*w4Ddx1+0.3535533905932737*f[53]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[6]*w4Ddx1+0.3535533905932737*f[21]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[47]*w4Ddx1+0.3535533905932737*f[59]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[17]*w4Ddx1+0.3535533905932737*f[38]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[18]*w4Ddx1+0.3535533905932737*f[39]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[32]*w4Ddx1+0.3535533905932737*f[51]*dv4Ddx1; 
  out[34] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[20]*w4Ddx1+0.3535533905932737*f[41]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[35]*w4Ddx1+0.3535533905932737*f[54]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[36]*w4Ddx1+0.3535533905932737*f[55]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[37]*w4Ddx1+0.3535533905932737*f[56]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[48]*w4Ddx1+0.3535533905932737*f[60]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[49]*w4Ddx1+0.3535533905932737*f[61]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[50]*w4Ddx1+0.3535533905932737*f[62]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[58]*w4Ddx1+0.3535533905932737*f[63]*dv4Ddx1; 
  out[36] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[21]*w4Ddx1+(0.3162277660168379*f[101]+0.3535533905932737*f[6])*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[38]*w4Ddx1+(0.3162277660168379*f[108]+0.3535533905932737*f[17])*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[39]*w4Ddx1+(0.3162277660168379*f[109]+0.3535533905932737*f[18])*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[40]*w4Ddx1+(0.3162277660168379*f[110]+0.3535533905932737*f[19])*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[51]*w4Ddx1+(0.3162277660168379*f[116]+0.3535533905932737*f[32])*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[52]*w4Ddx1+(0.3162277660168379*f[117]+0.3535533905932737*f[33])*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[53]*w4Ddx1+(0.3162277660168379*f[118]+0.3535533905932737*f[34])*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[59]*w4Ddx1+(0.3162277660168379*f[123]+0.3535533905932737*f[47])*dv4Ddx1; 
  out[39] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[24]*w4Ddx1+0.3535533905932737*f[45]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[12]*w4Ddx1+0.3535533905932737*f[31]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[42]*w4Ddx1+0.3535533905932737*f[57]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[10]*w4Ddx1+0.3535533905932737*f[29]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[25]*w4Ddx1+0.3535533905932737*f[46]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[4]*w4Ddx1+0.3535533905932737*f[16]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[23]*w4Ddx1+0.3535533905932737*f[44]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[11]*w4Ddx1+0.3535533905932737*f[30]*dv4Ddx1; 
  out[42] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[27]*w4Ddx1+(0.3162277660168379*f[103]+0.3535533905932737*f[8])*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[15]*w4Ddx1+(0.3162277660168379*f[99]+0.3535533905932737*f[3])*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[43]*w4Ddx1+(0.3162277660168379*f[112]+0.3535533905932737*f[22])*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[13]*w4Ddx1+(0.3162277660168379*f[97]+0.3535533905932737*f[1])*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[28]*w4Ddx1+(0.3162277660168379*f[104]+0.3535533905932737*f[9])*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[5]*w4Ddx1+(0.3162277660168379*f[96]+0.3535533905932737*f[0])*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[26]*w4Ddx1+(0.3162277660168379*f[102]+0.3535533905932737*f[7])*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[14]*w4Ddx1+(0.3162277660168379*f[98]+0.3535533905932737*f[2])*dv4Ddx1; 
  out[43] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[29]*w4Ddx1+(0.3162277660168379*f[105]+0.3535533905932737*f[10])*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[16]*w4Ddx1+(0.3162277660168379*f[100]+0.3535533905932737*f[4])*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[44]*w4Ddx1+(0.3162277660168379*f[113]+0.3535533905932737*f[23])*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[45]*w4Ddx1+(0.3162277660168379*f[114]+0.3535533905932737*f[24])*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[30]*w4Ddx1+(0.3162277660168379*f[106]+0.3535533905932737*f[11])*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[31]*w4Ddx1+(0.3162277660168379*f[107]+0.3535533905932737*f[12])*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[57]*w4Ddx1+(0.3162277660168379*f[122]+0.3535533905932737*f[42])*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[46]*w4Ddx1+(0.3162277660168379*f[115]+0.3535533905932737*f[25])*dv4Ddx1; 
  out[44] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[31]*w4Ddx1+(0.3162277660168379*f[107]+0.3535533905932737*f[12])*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[45]*w4Ddx1+(0.3162277660168379*f[114]+0.3535533905932737*f[24])*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[46]*w4Ddx1+(0.3162277660168379*f[115]+0.3535533905932737*f[25])*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[16]*w4Ddx1+(0.3162277660168379*f[100]+0.3535533905932737*f[4])*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[57]*w4Ddx1+(0.3162277660168379*f[122]+0.3535533905932737*f[42])*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[29]*w4Ddx1+(0.3162277660168379*f[105]+0.3535533905932737*f[10])*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[30]*w4Ddx1+(0.3162277660168379*f[106]+0.3535533905932737*f[11])*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[44]*w4Ddx1+(0.3162277660168379*f[113]+0.3535533905932737*f[23])*dv4Ddx1; 
  out[46] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[33]*w4Ddx1+0.3535533905932737*f[52]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[19]*w4Ddx1+0.3535533905932737*f[40]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[47]*w4Ddx1+0.3535533905932737*f[59]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[17]*w4Ddx1+0.3535533905932737*f[38]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[34]*w4Ddx1+0.3535533905932737*f[53]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[6]*w4Ddx1+0.3535533905932737*f[21]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[32]*w4Ddx1+0.3535533905932737*f[51]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[18]*w4Ddx1+0.3535533905932737*f[39]*dv4Ddx1; 
  out[47] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[35]*w4Ddx1+0.3535533905932737*f[54]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[20]*w4Ddx1+0.3535533905932737*f[41]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[48]*w4Ddx1+0.3535533905932737*f[60]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[49]*w4Ddx1+0.3535533905932737*f[61]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[36]*w4Ddx1+0.3535533905932737*f[55]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[37]*w4Ddx1+0.3535533905932737*f[56]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[58]*w4Ddx1+0.3535533905932737*f[63]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[50]*w4Ddx1+0.3535533905932737*f[62]*dv4Ddx1; 
  out[48] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[37]*w4Ddx1+0.3535533905932737*f[56]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[49]*w4Ddx1+0.3535533905932737*f[61]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[50]*w4Ddx1+0.3535533905932737*f[62]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[20]*w4Ddx1+0.3535533905932737*f[41]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[58]*w4Ddx1+0.3535533905932737*f[63]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[35]*w4Ddx1+0.3535533905932737*f[54]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[36]*w4Ddx1+0.3535533905932737*f[55]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[48]*w4Ddx1+0.3535533905932737*f[60]*dv4Ddx1; 
  out[50] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[38]*w4Ddx1+(0.3162277660168379*f[108]+0.3535533905932737*f[17])*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[21]*w4Ddx1+(0.3162277660168379*f[101]+0.3535533905932737*f[6])*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[51]*w4Ddx1+(0.3162277660168379*f[116]+0.3535533905932737*f[32])*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[52]*w4Ddx1+(0.3162277660168379*f[117]+0.3535533905932737*f[33])*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[39]*w4Ddx1+(0.3162277660168379*f[109]+0.3535533905932737*f[18])*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[40]*w4Ddx1+(0.3162277660168379*f[110]+0.3535533905932737*f[19])*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[59]*w4Ddx1+(0.3162277660168379*f[123]+0.3535533905932737*f[47])*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[53]*w4Ddx1+(0.3162277660168379*f[118]+0.3535533905932737*f[34])*dv4Ddx1; 
  out[51] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[40]*w4Ddx1+(0.3162277660168379*f[110]+0.3535533905932737*f[19])*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[52]*w4Ddx1+(0.3162277660168379*f[117]+0.3535533905932737*f[33])*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[53]*w4Ddx1+(0.3162277660168379*f[118]+0.3535533905932737*f[34])*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[21]*w4Ddx1+(0.3162277660168379*f[101]+0.3535533905932737*f[6])*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[59]*w4Ddx1+(0.3162277660168379*f[123]+0.3535533905932737*f[47])*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[38]*w4Ddx1+(0.3162277660168379*f[108]+0.3535533905932737*f[17])*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[39]*w4Ddx1+(0.3162277660168379*f[109]+0.3535533905932737*f[18])*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[51]*w4Ddx1+(0.3162277660168379*f[116]+0.3535533905932737*f[32])*dv4Ddx1; 
  out[53] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[41]*w4Ddx1+(0.3162277660168379*f[111]+0.3535533905932737*f[20])*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[54]*w4Ddx1+(0.3162277660168379*f[119]+0.3535533905932737*f[35])*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[55]*w4Ddx1+(0.3162277660168379*f[120]+0.3535533905932737*f[36])*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[56]*w4Ddx1+(0.3162277660168379*f[121]+0.3535533905932737*f[37])*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[60]*w4Ddx1+(0.3162277660168379*f[124]+0.3535533905932737*f[48])*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[61]*w4Ddx1+(0.3162277660168379*f[125]+0.3535533905932737*f[49])*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[62]*w4Ddx1+(0.3162277660168379*f[126]+0.3535533905932737*f[50])*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[63]*w4Ddx1+(0.3162277660168379*f[127]+0.3535533905932737*f[58])*dv4Ddx1; 
  out[55] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[45]*w4Ddx1+(0.3162277660168379*f[114]+0.3535533905932737*f[24])*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[31]*w4Ddx1+(0.3162277660168379*f[107]+0.3535533905932737*f[12])*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[57]*w4Ddx1+(0.3162277660168379*f[122]+0.3535533905932737*f[42])*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[29]*w4Ddx1+(0.3162277660168379*f[105]+0.3535533905932737*f[10])*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[46]*w4Ddx1+(0.3162277660168379*f[115]+0.3535533905932737*f[25])*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[16]*w4Ddx1+(0.3162277660168379*f[100]+0.3535533905932737*f[4])*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[44]*w4Ddx1+(0.3162277660168379*f[113]+0.3535533905932737*f[23])*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[30]*w4Ddx1+(0.3162277660168379*f[106]+0.3535533905932737*f[11])*dv4Ddx1; 
  out[57] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[49]*w4Ddx1+0.3535533905932737*f[61]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[37]*w4Ddx1+0.3535533905932737*f[56]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[58]*w4Ddx1+0.3535533905932737*f[63]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[35]*w4Ddx1+0.3535533905932737*f[54]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[50]*w4Ddx1+0.3535533905932737*f[62]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[20]*w4Ddx1+0.3535533905932737*f[41]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[48]*w4Ddx1+0.3535533905932737*f[60]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[36]*w4Ddx1+0.3535533905932737*f[55]*dv4Ddx1; 
  out[58] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[52]*w4Ddx1+(0.3162277660168379*f[117]+0.3535533905932737*f[33])*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[40]*w4Ddx1+(0.3162277660168379*f[110]+0.3535533905932737*f[19])*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[59]*w4Ddx1+(0.3162277660168379*f[123]+0.3535533905932737*f[47])*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[38]*w4Ddx1+(0.3162277660168379*f[108]+0.3535533905932737*f[17])*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[53]*w4Ddx1+(0.3162277660168379*f[118]+0.3535533905932737*f[34])*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[21]*w4Ddx1+(0.3162277660168379*f[101]+0.3535533905932737*f[6])*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[51]*w4Ddx1+(0.3162277660168379*f[116]+0.3535533905932737*f[32])*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[39]*w4Ddx1+(0.3162277660168379*f[109]+0.3535533905932737*f[18])*dv4Ddx1; 
  out[59] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[54]*w4Ddx1+(0.3162277660168379*f[119]+0.3535533905932737*f[35])*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[41]*w4Ddx1+(0.3162277660168379*f[111]+0.3535533905932737*f[20])*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[60]*w4Ddx1+(0.3162277660168379*f[124]+0.3535533905932737*f[48])*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[61]*w4Ddx1+(0.3162277660168379*f[125]+0.3535533905932737*f[49])*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[55]*w4Ddx1+(0.3162277660168379*f[120]+0.3535533905932737*f[36])*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[56]*w4Ddx1+(0.3162277660168379*f[121]+0.3535533905932737*f[37])*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[63]*w4Ddx1+(0.3162277660168379*f[127]+0.3535533905932737*f[58])*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[62]*w4Ddx1+(0.3162277660168379*f[126]+0.3535533905932737*f[50])*dv4Ddx1; 
  out[60] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[56]*w4Ddx1+(0.3162277660168379*f[121]+0.3535533905932737*f[37])*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[61]*w4Ddx1+(0.3162277660168379*f[125]+0.3535533905932737*f[49])*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[62]*w4Ddx1+(0.3162277660168379*f[126]+0.3535533905932737*f[50])*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[41]*w4Ddx1+(0.3162277660168379*f[111]+0.3535533905932737*f[20])*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[63]*w4Ddx1+(0.3162277660168379*f[127]+0.3535533905932737*f[58])*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[54]*w4Ddx1+(0.3162277660168379*f[119]+0.3535533905932737*f[35])*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[55]*w4Ddx1+(0.3162277660168379*f[120]+0.3535533905932737*f[36])*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[60]*w4Ddx1+(0.3162277660168379*f[124]+0.3535533905932737*f[48])*dv4Ddx1; 
  out[62] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[61]*w4Ddx1+(0.3162277660168379*f[125]+0.3535533905932737*f[49])*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[56]*w4Ddx1+(0.3162277660168379*f[121]+0.3535533905932737*f[37])*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[63]*w4Ddx1+(0.3162277660168379*f[127]+0.3535533905932737*f[58])*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[54]*w4Ddx1+(0.3162277660168379*f[119]+0.3535533905932737*f[35])*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[62]*w4Ddx1+(0.3162277660168379*f[126]+0.3535533905932737*f[50])*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[41]*w4Ddx1+(0.3162277660168379*f[111]+0.3535533905932737*f[20])*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[60]*w4Ddx1+(0.3162277660168379*f[124]+0.3535533905932737*f[48])*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[55]*w4Ddx1+(0.3162277660168379*f[120]+0.3535533905932737*f[36])*dv4Ddx1; 
  out[63] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[64]*w4Ddx1+0.3535533905932737*f[68]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[65]*w4Ddx1+0.3535533905932737*f[73]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[66]*w4Ddx1+0.3535533905932737*f[74]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[67]*w4Ddx1+0.3535533905932737*f[75]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[70]*w4Ddx1+0.3535533905932737*f[81]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[71]*w4Ddx1+0.3535533905932737*f[82]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[72]*w4Ddx1+0.3535533905932737*f[83]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[80]*w4Ddx1+0.3535533905932737*f[90]*dv4Ddx1; 
  out[66] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[65]*w4Ddx1+0.3535533905932737*f[73]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[64]*w4Ddx1+0.3535533905932737*f[68]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[70]*w4Ddx1+0.3535533905932737*f[81]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[71]*w4Ddx1+0.3535533905932737*f[82]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[66]*w4Ddx1+0.3535533905932737*f[74]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[67]*w4Ddx1+0.3535533905932737*f[75]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[80]*w4Ddx1+0.3535533905932737*f[90]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[72]*w4Ddx1+0.3535533905932737*f[83]*dv4Ddx1; 
  out[70] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[67]*w4Ddx1+0.3535533905932737*f[75]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[71]*w4Ddx1+0.3535533905932737*f[82]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[72]*w4Ddx1+0.3535533905932737*f[83]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[64]*w4Ddx1+0.3535533905932737*f[68]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[80]*w4Ddx1+0.3535533905932737*f[90]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[65]*w4Ddx1+0.3535533905932737*f[73]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[66]*w4Ddx1+0.3535533905932737*f[74]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[70]*w4Ddx1+0.3535533905932737*f[81]*dv4Ddx1; 
  out[72] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[68]*w4Ddx1+0.3535533905932737*f[64]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[73]*w4Ddx1+0.3535533905932737*f[65]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[74]*w4Ddx1+0.3535533905932737*f[66]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[75]*w4Ddx1+0.3535533905932737*f[67]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[81]*w4Ddx1+0.3535533905932737*f[70]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[82]*w4Ddx1+0.3535533905932737*f[71]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[83]*w4Ddx1+0.3535533905932737*f[72]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[90]*w4Ddx1+0.3535533905932737*f[80]*dv4Ddx1; 
  out[74] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[69]*w4Ddx1+0.3535533905932737*f[79]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[76]*w4Ddx1+0.3535533905932737*f[87]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[77]*w4Ddx1+0.3535533905932737*f[88]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[78]*w4Ddx1+0.3535533905932737*f[89]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[84]*w4Ddx1+0.3535533905932737*f[92]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[85]*w4Ddx1+0.3535533905932737*f[93]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[86]*w4Ddx1+0.3535533905932737*f[94]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[91]*w4Ddx1+0.3535533905932737*f[95]*dv4Ddx1; 
  out[77] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[71]*w4Ddx1+0.3535533905932737*f[82]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[67]*w4Ddx1+0.3535533905932737*f[75]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[80]*w4Ddx1+0.3535533905932737*f[90]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[65]*w4Ddx1+0.3535533905932737*f[73]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[72]*w4Ddx1+0.3535533905932737*f[83]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[64]*w4Ddx1+0.3535533905932737*f[68]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[70]*w4Ddx1+0.3535533905932737*f[81]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[66]*w4Ddx1+0.3535533905932737*f[74]*dv4Ddx1; 
  out[80] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[73]*w4Ddx1+0.3535533905932737*f[65]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[68]*w4Ddx1+0.3535533905932737*f[64]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[81]*w4Ddx1+0.3535533905932737*f[70]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[82]*w4Ddx1+0.3535533905932737*f[71]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[74]*w4Ddx1+0.3535533905932737*f[66]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[75]*w4Ddx1+0.3535533905932737*f[67]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[90]*w4Ddx1+0.3535533905932737*f[80]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[83]*w4Ddx1+0.3535533905932737*f[72]*dv4Ddx1; 
  out[81] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[75]*w4Ddx1+0.3535533905932737*f[67]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[82]*w4Ddx1+0.3535533905932737*f[71]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[83]*w4Ddx1+0.3535533905932737*f[72]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[68]*w4Ddx1+0.3535533905932737*f[64]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[90]*w4Ddx1+0.3535533905932737*f[80]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[73]*w4Ddx1+0.3535533905932737*f[65]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[74]*w4Ddx1+0.3535533905932737*f[66]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[81]*w4Ddx1+0.3535533905932737*f[70]*dv4Ddx1; 
  out[83] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[76]*w4Ddx1+0.3535533905932737*f[87]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[69]*w4Ddx1+0.3535533905932737*f[79]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[84]*w4Ddx1+0.3535533905932737*f[92]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[85]*w4Ddx1+0.3535533905932737*f[93]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[77]*w4Ddx1+0.3535533905932737*f[88]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[78]*w4Ddx1+0.3535533905932737*f[89]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[91]*w4Ddx1+0.3535533905932737*f[95]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[86]*w4Ddx1+0.3535533905932737*f[94]*dv4Ddx1; 
  out[84] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[78]*w4Ddx1+0.3535533905932737*f[89]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[85]*w4Ddx1+0.3535533905932737*f[93]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[86]*w4Ddx1+0.3535533905932737*f[94]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[69]*w4Ddx1+0.3535533905932737*f[79]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[91]*w4Ddx1+0.3535533905932737*f[95]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[76]*w4Ddx1+0.3535533905932737*f[87]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[77]*w4Ddx1+0.3535533905932737*f[88]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[84]*w4Ddx1+0.3535533905932737*f[92]*dv4Ddx1; 
  out[86] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[79]*w4Ddx1+0.3535533905932737*f[69]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[87]*w4Ddx1+0.3535533905932737*f[76]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[88]*w4Ddx1+0.3535533905932737*f[77]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[89]*w4Ddx1+0.3535533905932737*f[78]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[92]*w4Ddx1+0.3535533905932737*f[84]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[93]*w4Ddx1+0.3535533905932737*f[85]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[94]*w4Ddx1+0.3535533905932737*f[86]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[95]*w4Ddx1+0.3535533905932737*f[91]*dv4Ddx1; 
  out[88] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[82]*w4Ddx1+0.3535533905932737*f[71]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[75]*w4Ddx1+0.3535533905932737*f[67]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[90]*w4Ddx1+0.3535533905932737*f[80]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[73]*w4Ddx1+0.3535533905932737*f[65]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[83]*w4Ddx1+0.3535533905932737*f[72]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[68]*w4Ddx1+0.3535533905932737*f[64]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[81]*w4Ddx1+0.3535533905932737*f[70]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[74]*w4Ddx1+0.3535533905932737*f[66]*dv4Ddx1; 
  out[90] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[85]*w4Ddx1+0.3535533905932737*f[93]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[78]*w4Ddx1+0.3535533905932737*f[89]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[91]*w4Ddx1+0.3535533905932737*f[95]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[76]*w4Ddx1+0.3535533905932737*f[87]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[86]*w4Ddx1+0.3535533905932737*f[94]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[69]*w4Ddx1+0.3535533905932737*f[79]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[84]*w4Ddx1+0.3535533905932737*f[92]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[77]*w4Ddx1+0.3535533905932737*f[88]*dv4Ddx1; 
  out[91] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[87]*w4Ddx1+0.3535533905932737*f[76]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[79]*w4Ddx1+0.3535533905932737*f[69]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[92]*w4Ddx1+0.3535533905932737*f[84]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[93]*w4Ddx1+0.3535533905932737*f[85]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[88]*w4Ddx1+0.3535533905932737*f[77]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[89]*w4Ddx1+0.3535533905932737*f[78]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[95]*w4Ddx1+0.3535533905932737*f[91]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[94]*w4Ddx1+0.3535533905932737*f[86]*dv4Ddx1; 
  out[92] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[89]*w4Ddx1+0.3535533905932737*f[78]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[93]*w4Ddx1+0.3535533905932737*f[85]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[94]*w4Ddx1+0.3535533905932737*f[86]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[79]*w4Ddx1+0.3535533905932737*f[69]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[95]*w4Ddx1+0.3535533905932737*f[91]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[87]*w4Ddx1+0.3535533905932737*f[76]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[88]*w4Ddx1+0.3535533905932737*f[77]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[92]*w4Ddx1+0.3535533905932737*f[84]*dv4Ddx1; 
  out[94] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[93]*w4Ddx1+0.3535533905932737*f[85]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[89]*w4Ddx1+0.3535533905932737*f[78]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[95]*w4Ddx1+0.3535533905932737*f[91]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[87]*w4Ddx1+0.3535533905932737*f[76]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[94]*w4Ddx1+0.3535533905932737*f[86]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[79]*w4Ddx1+0.3535533905932737*f[69]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[92]*w4Ddx1+0.3535533905932737*f[84]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[88]*w4Ddx1+0.3535533905932737*f[77]*dv4Ddx1; 
  out[95] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[96]*w4Ddx1+0.3162277660168379*f[5]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[97]*w4Ddx1+0.3162277660168379*f[13]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[98]*w4Ddx1+0.3162277660168379*f[14]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[99]*w4Ddx1+0.3162277660168379*f[15]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[102]*w4Ddx1+0.3162277660168379*f[26]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[103]*w4Ddx1+0.3162277660168379*f[27]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[104]*w4Ddx1+0.3162277660168379*f[28]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[112]*w4Ddx1+0.3162277660168379*f[43]*dv4Ddx1; 
  out[98] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[97]*w4Ddx1+0.3162277660168379*f[13]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[96]*w4Ddx1+0.3162277660168379*f[5]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[102]*w4Ddx1+0.3162277660168379*f[26]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[103]*w4Ddx1+0.3162277660168379*f[27]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[98]*w4Ddx1+0.3162277660168379*f[14]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[99]*w4Ddx1+0.3162277660168379*f[15]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[112]*w4Ddx1+0.3162277660168379*f[43]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[104]*w4Ddx1+0.3162277660168379*f[28]*dv4Ddx1; 
  out[102] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[99]*w4Ddx1+0.3162277660168379*f[15]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[103]*w4Ddx1+0.3162277660168379*f[27]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[104]*w4Ddx1+0.3162277660168379*f[28]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[96]*w4Ddx1+0.3162277660168379*f[5]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[112]*w4Ddx1+0.3162277660168379*f[43]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[97]*w4Ddx1+0.3162277660168379*f[13]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[98]*w4Ddx1+0.3162277660168379*f[14]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[102]*w4Ddx1+0.3162277660168379*f[26]*dv4Ddx1; 
  out[104] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[100]*w4Ddx1+0.3162277660168379*f[16]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[105]*w4Ddx1+0.3162277660168379*f[29]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[106]*w4Ddx1+0.3162277660168379*f[30]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[107]*w4Ddx1+0.3162277660168379*f[31]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[113]*w4Ddx1+0.3162277660168379*f[44]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[114]*w4Ddx1+0.3162277660168379*f[45]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[115]*w4Ddx1+0.3162277660168379*f[46]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[122]*w4Ddx1+0.3162277660168379*f[57]*dv4Ddx1; 
  out[106] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[101]*w4Ddx1+0.3162277660168379*f[21]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[108]*w4Ddx1+0.3162277660168379*f[38]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[109]*w4Ddx1+0.3162277660168379*f[39]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[110]*w4Ddx1+0.3162277660168379*f[40]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[116]*w4Ddx1+0.3162277660168379*f[51]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[117]*w4Ddx1+0.3162277660168379*f[52]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[118]*w4Ddx1+0.3162277660168379*f[53]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[123]*w4Ddx1+0.3162277660168379*f[59]*dv4Ddx1; 
  out[109] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[103]*w4Ddx1+0.3162277660168379*f[27]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[99]*w4Ddx1+0.3162277660168379*f[15]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[112]*w4Ddx1+0.3162277660168379*f[43]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[97]*w4Ddx1+0.3162277660168379*f[13]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[104]*w4Ddx1+0.3162277660168379*f[28]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[96]*w4Ddx1+0.3162277660168379*f[5]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[102]*w4Ddx1+0.3162277660168379*f[26]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[98]*w4Ddx1+0.3162277660168379*f[14]*dv4Ddx1; 
  out[112] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[105]*w4Ddx1+0.3162277660168379*f[29]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[100]*w4Ddx1+0.3162277660168379*f[16]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[113]*w4Ddx1+0.3162277660168379*f[44]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[114]*w4Ddx1+0.3162277660168379*f[45]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[106]*w4Ddx1+0.3162277660168379*f[30]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[107]*w4Ddx1+0.3162277660168379*f[31]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[122]*w4Ddx1+0.3162277660168379*f[57]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[115]*w4Ddx1+0.3162277660168379*f[46]*dv4Ddx1; 
  out[113] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[107]*w4Ddx1+0.3162277660168379*f[31]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[114]*w4Ddx1+0.3162277660168379*f[45]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[115]*w4Ddx1+0.3162277660168379*f[46]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[100]*w4Ddx1+0.3162277660168379*f[16]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[122]*w4Ddx1+0.3162277660168379*f[57]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[105]*w4Ddx1+0.3162277660168379*f[29]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[106]*w4Ddx1+0.3162277660168379*f[30]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[113]*w4Ddx1+0.3162277660168379*f[44]*dv4Ddx1; 
  out[115] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[108]*w4Ddx1+0.3162277660168379*f[38]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[101]*w4Ddx1+0.3162277660168379*f[21]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[116]*w4Ddx1+0.3162277660168379*f[51]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[117]*w4Ddx1+0.3162277660168379*f[52]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[109]*w4Ddx1+0.3162277660168379*f[39]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[110]*w4Ddx1+0.3162277660168379*f[40]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[123]*w4Ddx1+0.3162277660168379*f[59]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[118]*w4Ddx1+0.3162277660168379*f[53]*dv4Ddx1; 
  out[116] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[110]*w4Ddx1+0.3162277660168379*f[40]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[117]*w4Ddx1+0.3162277660168379*f[52]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[118]*w4Ddx1+0.3162277660168379*f[53]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[101]*w4Ddx1+0.3162277660168379*f[21]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[123]*w4Ddx1+0.3162277660168379*f[59]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[108]*w4Ddx1+0.3162277660168379*f[38]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[109]*w4Ddx1+0.3162277660168379*f[39]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[116]*w4Ddx1+0.3162277660168379*f[51]*dv4Ddx1; 
  out[118] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[111]*w4Ddx1+0.3162277660168379*f[41]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[119]*w4Ddx1+0.3162277660168379*f[54]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[120]*w4Ddx1+0.3162277660168379*f[55]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[121]*w4Ddx1+0.3162277660168379*f[56]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[124]*w4Ddx1+0.3162277660168379*f[60]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[125]*w4Ddx1+0.3162277660168379*f[61]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[126]*w4Ddx1+0.3162277660168379*f[62]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[127]*w4Ddx1+0.3162277660168379*f[63]*dv4Ddx1; 
  out[120] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[114]*w4Ddx1+0.3162277660168379*f[45]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[107]*w4Ddx1+0.3162277660168379*f[31]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[122]*w4Ddx1+0.3162277660168379*f[57]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[105]*w4Ddx1+0.3162277660168379*f[29]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[115]*w4Ddx1+0.3162277660168379*f[46]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[100]*w4Ddx1+0.3162277660168379*f[16]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[113]*w4Ddx1+0.3162277660168379*f[44]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[106]*w4Ddx1+0.3162277660168379*f[30]*dv4Ddx1; 
  out[122] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[117]*w4Ddx1+0.3162277660168379*f[52]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[110]*w4Ddx1+0.3162277660168379*f[40]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[123]*w4Ddx1+0.3162277660168379*f[59]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[108]*w4Ddx1+0.3162277660168379*f[38]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[118]*w4Ddx1+0.3162277660168379*f[53]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[101]*w4Ddx1+0.3162277660168379*f[21]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[116]*w4Ddx1+0.3162277660168379*f[51]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[109]*w4Ddx1+0.3162277660168379*f[39]*dv4Ddx1; 
  out[123] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[119]*w4Ddx1+0.3162277660168379*f[54]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[111]*w4Ddx1+0.3162277660168379*f[41]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[124]*w4Ddx1+0.3162277660168379*f[60]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[125]*w4Ddx1+0.3162277660168379*f[61]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[120]*w4Ddx1+0.3162277660168379*f[55]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[121]*w4Ddx1+0.3162277660168379*f[56]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[127]*w4Ddx1+0.3162277660168379*f[63]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[126]*w4Ddx1+0.3162277660168379*f[62]*dv4Ddx1; 
  out[124] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[121]*w4Ddx1+0.3162277660168379*f[56]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[125]*w4Ddx1+0.3162277660168379*f[61]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[126]*w4Ddx1+0.3162277660168379*f[62]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[111]*w4Ddx1+0.3162277660168379*f[41]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[127]*w4Ddx1+0.3162277660168379*f[63]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[119]*w4Ddx1+0.3162277660168379*f[54]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[120]*w4Ddx1+0.3162277660168379*f[55]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[124]*w4Ddx1+0.3162277660168379*f[60]*dv4Ddx1; 
  out[126] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[125]*w4Ddx1+0.3162277660168379*f[61]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[121]*w4Ddx1+0.3162277660168379*f[56]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[127]*w4Ddx1+0.3162277660168379*f[63]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[119]*w4Ddx1+0.3162277660168379*f[54]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[126]*w4Ddx1+0.3162277660168379*f[62]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[111]*w4Ddx1+0.3162277660168379*f[41]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[124]*w4Ddx1+0.3162277660168379*f[60]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[120]*w4Ddx1+0.3162277660168379*f[55]*dv4Ddx1; 
  out[127] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[128]*w4Ddx1+0.3535533905932737*f[133]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[129]*w4Ddx1+0.3535533905932737*f[140]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[130]*w4Ddx1+0.3535533905932737*f[141]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[131]*w4Ddx1+0.3535533905932737*f[142]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[134]*w4Ddx1+0.3535533905932737*f[148]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[135]*w4Ddx1+0.3535533905932737*f[149]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[136]*w4Ddx1+0.3535533905932737*f[150]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[144]*w4Ddx1+0.3535533905932737*f[155]*dv4Ddx1; 
  out[130] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[129]*w4Ddx1+0.3535533905932737*f[140]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[128]*w4Ddx1+0.3535533905932737*f[133]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[134]*w4Ddx1+0.3535533905932737*f[148]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[135]*w4Ddx1+0.3535533905932737*f[149]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[130]*w4Ddx1+0.3535533905932737*f[141]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[131]*w4Ddx1+0.3535533905932737*f[142]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[144]*w4Ddx1+0.3535533905932737*f[155]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[136]*w4Ddx1+0.3535533905932737*f[150]*dv4Ddx1; 
  out[134] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[131]*w4Ddx1+0.3535533905932737*f[142]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[135]*w4Ddx1+0.3535533905932737*f[149]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[136]*w4Ddx1+0.3535533905932737*f[150]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[128]*w4Ddx1+0.3535533905932737*f[133]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[144]*w4Ddx1+0.3535533905932737*f[155]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[129]*w4Ddx1+0.3535533905932737*f[140]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[130]*w4Ddx1+0.3535533905932737*f[141]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[134]*w4Ddx1+0.3535533905932737*f[148]*dv4Ddx1; 
  out[136] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[132]*w4Ddx1+0.3535533905932737*f[143]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[137]*w4Ddx1+0.3535533905932737*f[151]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[138]*w4Ddx1+0.3535533905932737*f[152]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[139]*w4Ddx1+0.3535533905932737*f[153]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[145]*w4Ddx1+0.3535533905932737*f[156]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[146]*w4Ddx1+0.3535533905932737*f[157]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[147]*w4Ddx1+0.3535533905932737*f[158]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[154]*w4Ddx1+0.3535533905932737*f[159]*dv4Ddx1; 
  out[138] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[133]*w4Ddx1+0.3535533905932737*f[128]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[140]*w4Ddx1+0.3535533905932737*f[129]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[141]*w4Ddx1+0.3535533905932737*f[130]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[142]*w4Ddx1+0.3535533905932737*f[131]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[148]*w4Ddx1+0.3535533905932737*f[134]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[149]*w4Ddx1+0.3535533905932737*f[135]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[150]*w4Ddx1+0.3535533905932737*f[136]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[155]*w4Ddx1+0.3535533905932737*f[144]*dv4Ddx1; 
  out[141] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[135]*w4Ddx1+0.3535533905932737*f[149]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[131]*w4Ddx1+0.3535533905932737*f[142]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[144]*w4Ddx1+0.3535533905932737*f[155]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[129]*w4Ddx1+0.3535533905932737*f[140]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[136]*w4Ddx1+0.3535533905932737*f[150]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[128]*w4Ddx1+0.3535533905932737*f[133]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[134]*w4Ddx1+0.3535533905932737*f[148]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[130]*w4Ddx1+0.3535533905932737*f[141]*dv4Ddx1; 
  out[144] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[137]*w4Ddx1+0.3535533905932737*f[151]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[132]*w4Ddx1+0.3535533905932737*f[143]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[145]*w4Ddx1+0.3535533905932737*f[156]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[146]*w4Ddx1+0.3535533905932737*f[157]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[138]*w4Ddx1+0.3535533905932737*f[152]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[139]*w4Ddx1+0.3535533905932737*f[153]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[154]*w4Ddx1+0.3535533905932737*f[159]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[147]*w4Ddx1+0.3535533905932737*f[158]*dv4Ddx1; 
  out[145] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[139]*w4Ddx1+0.3535533905932737*f[153]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[146]*w4Ddx1+0.3535533905932737*f[157]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[147]*w4Ddx1+0.3535533905932737*f[158]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[132]*w4Ddx1+0.3535533905932737*f[143]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[154]*w4Ddx1+0.3535533905932737*f[159]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[137]*w4Ddx1+0.3535533905932737*f[151]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[138]*w4Ddx1+0.3535533905932737*f[152]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[145]*w4Ddx1+0.3535533905932737*f[156]*dv4Ddx1; 
  out[147] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[140]*w4Ddx1+0.3535533905932737*f[129]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[133]*w4Ddx1+0.3535533905932737*f[128]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[148]*w4Ddx1+0.3535533905932737*f[134]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[149]*w4Ddx1+0.3535533905932737*f[135]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[141]*w4Ddx1+0.3535533905932737*f[130]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[142]*w4Ddx1+0.3535533905932737*f[131]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[155]*w4Ddx1+0.3535533905932737*f[144]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[150]*w4Ddx1+0.3535533905932737*f[136]*dv4Ddx1; 
  out[148] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[142]*w4Ddx1+0.3535533905932737*f[131]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[149]*w4Ddx1+0.3535533905932737*f[135]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[150]*w4Ddx1+0.3535533905932737*f[136]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[133]*w4Ddx1+0.3535533905932737*f[128]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[155]*w4Ddx1+0.3535533905932737*f[144]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[140]*w4Ddx1+0.3535533905932737*f[129]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[141]*w4Ddx1+0.3535533905932737*f[130]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[148]*w4Ddx1+0.3535533905932737*f[134]*dv4Ddx1; 
  out[150] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[143]*w4Ddx1+0.3535533905932737*f[132]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[151]*w4Ddx1+0.3535533905932737*f[137]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[152]*w4Ddx1+0.3535533905932737*f[138]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[153]*w4Ddx1+0.3535533905932737*f[139]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[156]*w4Ddx1+0.3535533905932737*f[145]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[157]*w4Ddx1+0.3535533905932737*f[146]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[158]*w4Ddx1+0.3535533905932737*f[147]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[159]*w4Ddx1+0.3535533905932737*f[154]*dv4Ddx1; 
  out[152] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[146]*w4Ddx1+0.3535533905932737*f[157]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[139]*w4Ddx1+0.3535533905932737*f[153]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[154]*w4Ddx1+0.3535533905932737*f[159]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[137]*w4Ddx1+0.3535533905932737*f[151]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[147]*w4Ddx1+0.3535533905932737*f[158]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[132]*w4Ddx1+0.3535533905932737*f[143]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[145]*w4Ddx1+0.3535533905932737*f[156]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[138]*w4Ddx1+0.3535533905932737*f[152]*dv4Ddx1; 
  out[154] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[149]*w4Ddx1+0.3535533905932737*f[135]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[142]*w4Ddx1+0.3535533905932737*f[131]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[155]*w4Ddx1+0.3535533905932737*f[144]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[140]*w4Ddx1+0.3535533905932737*f[129]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[150]*w4Ddx1+0.3535533905932737*f[136]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[133]*w4Ddx1+0.3535533905932737*f[128]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[148]*w4Ddx1+0.3535533905932737*f[134]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[141]*w4Ddx1+0.3535533905932737*f[130]*dv4Ddx1; 
  out[155] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[151]*w4Ddx1+0.3535533905932737*f[137]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[143]*w4Ddx1+0.3535533905932737*f[132]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[156]*w4Ddx1+0.3535533905932737*f[145]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[157]*w4Ddx1+0.3535533905932737*f[146]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[152]*w4Ddx1+0.3535533905932737*f[138]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[153]*w4Ddx1+0.3535533905932737*f[139]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[159]*w4Ddx1+0.3535533905932737*f[154]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[158]*w4Ddx1+0.3535533905932737*f[147]*dv4Ddx1; 
  out[156] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[153]*w4Ddx1+0.3535533905932737*f[139]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[157]*w4Ddx1+0.3535533905932737*f[146]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[158]*w4Ddx1+0.3535533905932737*f[147]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[143]*w4Ddx1+0.3535533905932737*f[132]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[159]*w4Ddx1+0.3535533905932737*f[154]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[151]*w4Ddx1+0.3535533905932737*f[137]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[152]*w4Ddx1+0.3535533905932737*f[138]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[156]*w4Ddx1+0.3535533905932737*f[145]*dv4Ddx1; 
  out[158] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[157]*w4Ddx1+0.3535533905932737*f[146]*dv4Ddx1; 
  Gbar[1] = 1.224744871391589*f[153]*w4Ddx1+0.3535533905932737*f[139]*dv4Ddx1; 
  Gbar[2] = 1.224744871391589*f[159]*w4Ddx1+0.3535533905932737*f[154]*dv4Ddx1; 
  Gbar[3] = 1.224744871391589*f[151]*w4Ddx1+0.3535533905932737*f[137]*dv4Ddx1; 
  Gbar[4] = 1.224744871391589*f[158]*w4Ddx1+0.3535533905932737*f[147]*dv4Ddx1; 
  Gbar[5] = 1.224744871391589*f[143]*w4Ddx1+0.3535533905932737*f[132]*dv4Ddx1; 
  Gbar[6] = 1.224744871391589*f[156]*w4Ddx1+0.3535533905932737*f[145]*dv4Ddx1; 
  Gbar[7] = 1.224744871391589*f[152]*w4Ddx1+0.3535533905932737*f[138]*dv4Ddx1; 
  out[159] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[0]*w5Ddx2+0.3535533905932737*f[6]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[1]*w5Ddx2+0.3535533905932737*f[17]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[2]*w5Ddx2+0.3535533905932737*f[18]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[3]*w5Ddx2+0.3535533905932737*f[19]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[7]*w5Ddx2+0.3535533905932737*f[32]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[8]*w5Ddx2+0.3535533905932737*f[33]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[9]*w5Ddx2+0.3535533905932737*f[34]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[22]*w5Ddx2+0.3535533905932737*f[47]*dv5Ddx2; 
  out[3] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[1]*w5Ddx2+0.3535533905932737*f[17]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[0]*w5Ddx2+0.3535533905932737*f[6]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[7]*w5Ddx2+0.3535533905932737*f[32]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[8]*w5Ddx2+0.3535533905932737*f[33]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[2]*w5Ddx2+0.3535533905932737*f[18]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[3]*w5Ddx2+0.3535533905932737*f[19]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[22]*w5Ddx2+0.3535533905932737*f[47]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[9]*w5Ddx2+0.3535533905932737*f[34]*dv5Ddx2; 
  out[8] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[2]*w5Ddx2+0.3535533905932737*f[18]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[7]*w5Ddx2+0.3535533905932737*f[32]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[0]*w5Ddx2+0.3535533905932737*f[6]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[9]*w5Ddx2+0.3535533905932737*f[34]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[1]*w5Ddx2+0.3535533905932737*f[17]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[22]*w5Ddx2+0.3535533905932737*f[47]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[3]*w5Ddx2+0.3535533905932737*f[19]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[8]*w5Ddx2+0.3535533905932737*f[33]*dv5Ddx2; 
  out[9] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[4]*w5Ddx2+0.3535533905932737*f[20]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[10]*w5Ddx2+0.3535533905932737*f[35]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[11]*w5Ddx2+0.3535533905932737*f[36]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[12]*w5Ddx2+0.3535533905932737*f[37]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[23]*w5Ddx2+0.3535533905932737*f[48]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[24]*w5Ddx2+0.3535533905932737*f[49]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[25]*w5Ddx2+0.3535533905932737*f[50]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[42]*w5Ddx2+0.3535533905932737*f[58]*dv5Ddx2; 
  out[12] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[5]*w5Ddx2+0.3535533905932737*f[21]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[13]*w5Ddx2+0.3535533905932737*f[38]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[14]*w5Ddx2+0.3535533905932737*f[39]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[15]*w5Ddx2+0.3535533905932737*f[40]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[26]*w5Ddx2+0.3535533905932737*f[51]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[27]*w5Ddx2+0.3535533905932737*f[52]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[28]*w5Ddx2+0.3535533905932737*f[53]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[43]*w5Ddx2+0.3535533905932737*f[59]*dv5Ddx2; 
  out[15] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[6]*w5Ddx2+(0.3162277660168379*f[128]+0.3535533905932737*f[0])*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[17]*w5Ddx2+(0.3162277660168379*f[129]+0.3535533905932737*f[1])*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[18]*w5Ddx2+(0.3162277660168379*f[130]+0.3535533905932737*f[2])*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[19]*w5Ddx2+(0.3162277660168379*f[131]+0.3535533905932737*f[3])*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[32]*w5Ddx2+(0.3162277660168379*f[134]+0.3535533905932737*f[7])*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[33]*w5Ddx2+(0.3162277660168379*f[135]+0.3535533905932737*f[8])*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[34]*w5Ddx2+(0.3162277660168379*f[136]+0.3535533905932737*f[9])*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[47]*w5Ddx2+(0.3162277660168379*f[144]+0.3535533905932737*f[22])*dv5Ddx2; 
  out[19] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[7]*w5Ddx2+0.3535533905932737*f[32]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[2]*w5Ddx2+0.3535533905932737*f[18]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[1]*w5Ddx2+0.3535533905932737*f[17]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[22]*w5Ddx2+0.3535533905932737*f[47]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[0]*w5Ddx2+0.3535533905932737*f[6]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[9]*w5Ddx2+0.3535533905932737*f[34]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[8]*w5Ddx2+0.3535533905932737*f[33]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[3]*w5Ddx2+0.3535533905932737*f[19]*dv5Ddx2; 
  out[22] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[10]*w5Ddx2+0.3535533905932737*f[35]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[4]*w5Ddx2+0.3535533905932737*f[20]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[23]*w5Ddx2+0.3535533905932737*f[48]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[24]*w5Ddx2+0.3535533905932737*f[49]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[11]*w5Ddx2+0.3535533905932737*f[36]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[12]*w5Ddx2+0.3535533905932737*f[37]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[42]*w5Ddx2+0.3535533905932737*f[58]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[25]*w5Ddx2+0.3535533905932737*f[50]*dv5Ddx2; 
  out[24] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[11]*w5Ddx2+0.3535533905932737*f[36]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[23]*w5Ddx2+0.3535533905932737*f[48]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[4]*w5Ddx2+0.3535533905932737*f[20]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[25]*w5Ddx2+0.3535533905932737*f[50]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[10]*w5Ddx2+0.3535533905932737*f[35]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[42]*w5Ddx2+0.3535533905932737*f[58]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[12]*w5Ddx2+0.3535533905932737*f[37]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[24]*w5Ddx2+0.3535533905932737*f[49]*dv5Ddx2; 
  out[25] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[13]*w5Ddx2+0.3535533905932737*f[38]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[5]*w5Ddx2+0.3535533905932737*f[21]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[26]*w5Ddx2+0.3535533905932737*f[51]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[27]*w5Ddx2+0.3535533905932737*f[52]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[14]*w5Ddx2+0.3535533905932737*f[39]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[15]*w5Ddx2+0.3535533905932737*f[40]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[43]*w5Ddx2+0.3535533905932737*f[59]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[28]*w5Ddx2+0.3535533905932737*f[53]*dv5Ddx2; 
  out[27] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[14]*w5Ddx2+0.3535533905932737*f[39]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[26]*w5Ddx2+0.3535533905932737*f[51]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[5]*w5Ddx2+0.3535533905932737*f[21]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[28]*w5Ddx2+0.3535533905932737*f[53]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[13]*w5Ddx2+0.3535533905932737*f[38]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[43]*w5Ddx2+0.3535533905932737*f[59]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[15]*w5Ddx2+0.3535533905932737*f[40]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[27]*w5Ddx2+0.3535533905932737*f[52]*dv5Ddx2; 
  out[28] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[16]*w5Ddx2+0.3535533905932737*f[41]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[29]*w5Ddx2+0.3535533905932737*f[54]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[30]*w5Ddx2+0.3535533905932737*f[55]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[31]*w5Ddx2+0.3535533905932737*f[56]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[44]*w5Ddx2+0.3535533905932737*f[60]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[45]*w5Ddx2+0.3535533905932737*f[61]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[46]*w5Ddx2+0.3535533905932737*f[62]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[57]*w5Ddx2+0.3535533905932737*f[63]*dv5Ddx2; 
  out[31] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[17]*w5Ddx2+(0.3162277660168379*f[129]+0.3535533905932737*f[1])*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[6]*w5Ddx2+(0.3162277660168379*f[128]+0.3535533905932737*f[0])*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[32]*w5Ddx2+(0.3162277660168379*f[134]+0.3535533905932737*f[7])*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[33]*w5Ddx2+(0.3162277660168379*f[135]+0.3535533905932737*f[8])*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[18]*w5Ddx2+(0.3162277660168379*f[130]+0.3535533905932737*f[2])*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[19]*w5Ddx2+(0.3162277660168379*f[131]+0.3535533905932737*f[3])*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[47]*w5Ddx2+(0.3162277660168379*f[144]+0.3535533905932737*f[22])*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[34]*w5Ddx2+(0.3162277660168379*f[136]+0.3535533905932737*f[9])*dv5Ddx2; 
  out[33] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[18]*w5Ddx2+(0.3162277660168379*f[130]+0.3535533905932737*f[2])*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[32]*w5Ddx2+(0.3162277660168379*f[134]+0.3535533905932737*f[7])*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[6]*w5Ddx2+(0.3162277660168379*f[128]+0.3535533905932737*f[0])*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[34]*w5Ddx2+(0.3162277660168379*f[136]+0.3535533905932737*f[9])*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[17]*w5Ddx2+(0.3162277660168379*f[129]+0.3535533905932737*f[1])*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[47]*w5Ddx2+(0.3162277660168379*f[144]+0.3535533905932737*f[22])*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[19]*w5Ddx2+(0.3162277660168379*f[131]+0.3535533905932737*f[3])*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[33]*w5Ddx2+(0.3162277660168379*f[135]+0.3535533905932737*f[8])*dv5Ddx2; 
  out[34] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[20]*w5Ddx2+(0.3162277660168379*f[132]+0.3535533905932737*f[4])*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[35]*w5Ddx2+(0.3162277660168379*f[137]+0.3535533905932737*f[10])*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[36]*w5Ddx2+(0.3162277660168379*f[138]+0.3535533905932737*f[11])*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[37]*w5Ddx2+(0.3162277660168379*f[139]+0.3535533905932737*f[12])*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[48]*w5Ddx2+(0.3162277660168379*f[145]+0.3535533905932737*f[23])*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[49]*w5Ddx2+(0.3162277660168379*f[146]+0.3535533905932737*f[24])*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[50]*w5Ddx2+(0.3162277660168379*f[147]+0.3535533905932737*f[25])*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[58]*w5Ddx2+(0.3162277660168379*f[154]+0.3535533905932737*f[42])*dv5Ddx2; 
  out[37] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[21]*w5Ddx2+(0.3162277660168379*f[133]+0.3535533905932737*f[5])*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[38]*w5Ddx2+(0.3162277660168379*f[140]+0.3535533905932737*f[13])*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[39]*w5Ddx2+(0.3162277660168379*f[141]+0.3535533905932737*f[14])*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[40]*w5Ddx2+(0.3162277660168379*f[142]+0.3535533905932737*f[15])*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[51]*w5Ddx2+(0.3162277660168379*f[148]+0.3535533905932737*f[26])*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[52]*w5Ddx2+(0.3162277660168379*f[149]+0.3535533905932737*f[27])*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[53]*w5Ddx2+(0.3162277660168379*f[150]+0.3535533905932737*f[28])*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[59]*w5Ddx2+(0.3162277660168379*f[155]+0.3535533905932737*f[43])*dv5Ddx2; 
  out[40] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[23]*w5Ddx2+0.3535533905932737*f[48]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[11]*w5Ddx2+0.3535533905932737*f[36]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[10]*w5Ddx2+0.3535533905932737*f[35]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[42]*w5Ddx2+0.3535533905932737*f[58]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[4]*w5Ddx2+0.3535533905932737*f[20]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[25]*w5Ddx2+0.3535533905932737*f[50]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[24]*w5Ddx2+0.3535533905932737*f[49]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[12]*w5Ddx2+0.3535533905932737*f[37]*dv5Ddx2; 
  out[42] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[26]*w5Ddx2+0.3535533905932737*f[51]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[14]*w5Ddx2+0.3535533905932737*f[39]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[13]*w5Ddx2+0.3535533905932737*f[38]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[43]*w5Ddx2+0.3535533905932737*f[59]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[5]*w5Ddx2+0.3535533905932737*f[21]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[28]*w5Ddx2+0.3535533905932737*f[53]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[27]*w5Ddx2+0.3535533905932737*f[52]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[15]*w5Ddx2+0.3535533905932737*f[40]*dv5Ddx2; 
  out[43] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[29]*w5Ddx2+0.3535533905932737*f[54]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[16]*w5Ddx2+0.3535533905932737*f[41]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[44]*w5Ddx2+0.3535533905932737*f[60]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[45]*w5Ddx2+0.3535533905932737*f[61]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[30]*w5Ddx2+0.3535533905932737*f[55]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[31]*w5Ddx2+0.3535533905932737*f[56]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[57]*w5Ddx2+0.3535533905932737*f[63]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[46]*w5Ddx2+0.3535533905932737*f[62]*dv5Ddx2; 
  out[45] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[30]*w5Ddx2+0.3535533905932737*f[55]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[44]*w5Ddx2+0.3535533905932737*f[60]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[16]*w5Ddx2+0.3535533905932737*f[41]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[46]*w5Ddx2+0.3535533905932737*f[62]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[29]*w5Ddx2+0.3535533905932737*f[54]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[57]*w5Ddx2+0.3535533905932737*f[63]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[31]*w5Ddx2+0.3535533905932737*f[56]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[45]*w5Ddx2+0.3535533905932737*f[61]*dv5Ddx2; 
  out[46] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[32]*w5Ddx2+(0.3162277660168379*f[134]+0.3535533905932737*f[7])*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[18]*w5Ddx2+(0.3162277660168379*f[130]+0.3535533905932737*f[2])*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[17]*w5Ddx2+(0.3162277660168379*f[129]+0.3535533905932737*f[1])*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[47]*w5Ddx2+(0.3162277660168379*f[144]+0.3535533905932737*f[22])*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[6]*w5Ddx2+(0.3162277660168379*f[128]+0.3535533905932737*f[0])*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[34]*w5Ddx2+(0.3162277660168379*f[136]+0.3535533905932737*f[9])*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[33]*w5Ddx2+(0.3162277660168379*f[135]+0.3535533905932737*f[8])*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[19]*w5Ddx2+(0.3162277660168379*f[131]+0.3535533905932737*f[3])*dv5Ddx2; 
  out[47] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[35]*w5Ddx2+(0.3162277660168379*f[137]+0.3535533905932737*f[10])*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[20]*w5Ddx2+(0.3162277660168379*f[132]+0.3535533905932737*f[4])*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[48]*w5Ddx2+(0.3162277660168379*f[145]+0.3535533905932737*f[23])*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[49]*w5Ddx2+(0.3162277660168379*f[146]+0.3535533905932737*f[24])*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[36]*w5Ddx2+(0.3162277660168379*f[138]+0.3535533905932737*f[11])*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[37]*w5Ddx2+(0.3162277660168379*f[139]+0.3535533905932737*f[12])*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[58]*w5Ddx2+(0.3162277660168379*f[154]+0.3535533905932737*f[42])*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[50]*w5Ddx2+(0.3162277660168379*f[147]+0.3535533905932737*f[25])*dv5Ddx2; 
  out[49] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[36]*w5Ddx2+(0.3162277660168379*f[138]+0.3535533905932737*f[11])*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[48]*w5Ddx2+(0.3162277660168379*f[145]+0.3535533905932737*f[23])*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[20]*w5Ddx2+(0.3162277660168379*f[132]+0.3535533905932737*f[4])*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[50]*w5Ddx2+(0.3162277660168379*f[147]+0.3535533905932737*f[25])*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[35]*w5Ddx2+(0.3162277660168379*f[137]+0.3535533905932737*f[10])*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[58]*w5Ddx2+(0.3162277660168379*f[154]+0.3535533905932737*f[42])*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[37]*w5Ddx2+(0.3162277660168379*f[139]+0.3535533905932737*f[12])*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[49]*w5Ddx2+(0.3162277660168379*f[146]+0.3535533905932737*f[24])*dv5Ddx2; 
  out[50] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[38]*w5Ddx2+(0.3162277660168379*f[140]+0.3535533905932737*f[13])*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[21]*w5Ddx2+(0.3162277660168379*f[133]+0.3535533905932737*f[5])*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[51]*w5Ddx2+(0.3162277660168379*f[148]+0.3535533905932737*f[26])*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[52]*w5Ddx2+(0.3162277660168379*f[149]+0.3535533905932737*f[27])*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[39]*w5Ddx2+(0.3162277660168379*f[141]+0.3535533905932737*f[14])*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[40]*w5Ddx2+(0.3162277660168379*f[142]+0.3535533905932737*f[15])*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[59]*w5Ddx2+(0.3162277660168379*f[155]+0.3535533905932737*f[43])*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[53]*w5Ddx2+(0.3162277660168379*f[150]+0.3535533905932737*f[28])*dv5Ddx2; 
  out[52] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[39]*w5Ddx2+(0.3162277660168379*f[141]+0.3535533905932737*f[14])*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[51]*w5Ddx2+(0.3162277660168379*f[148]+0.3535533905932737*f[26])*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[21]*w5Ddx2+(0.3162277660168379*f[133]+0.3535533905932737*f[5])*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[53]*w5Ddx2+(0.3162277660168379*f[150]+0.3535533905932737*f[28])*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[38]*w5Ddx2+(0.3162277660168379*f[140]+0.3535533905932737*f[13])*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[59]*w5Ddx2+(0.3162277660168379*f[155]+0.3535533905932737*f[43])*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[40]*w5Ddx2+(0.3162277660168379*f[142]+0.3535533905932737*f[15])*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[52]*w5Ddx2+(0.3162277660168379*f[149]+0.3535533905932737*f[27])*dv5Ddx2; 
  out[53] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[41]*w5Ddx2+(0.3162277660168379*f[143]+0.3535533905932737*f[16])*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[54]*w5Ddx2+(0.3162277660168379*f[151]+0.3535533905932737*f[29])*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[55]*w5Ddx2+(0.3162277660168379*f[152]+0.3535533905932737*f[30])*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[56]*w5Ddx2+(0.3162277660168379*f[153]+0.3535533905932737*f[31])*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[60]*w5Ddx2+(0.3162277660168379*f[156]+0.3535533905932737*f[44])*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[61]*w5Ddx2+(0.3162277660168379*f[157]+0.3535533905932737*f[45])*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[62]*w5Ddx2+(0.3162277660168379*f[158]+0.3535533905932737*f[46])*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[63]*w5Ddx2+(0.3162277660168379*f[159]+0.3535533905932737*f[57])*dv5Ddx2; 
  out[56] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[44]*w5Ddx2+0.3535533905932737*f[60]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[30]*w5Ddx2+0.3535533905932737*f[55]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[29]*w5Ddx2+0.3535533905932737*f[54]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[57]*w5Ddx2+0.3535533905932737*f[63]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[16]*w5Ddx2+0.3535533905932737*f[41]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[46]*w5Ddx2+0.3535533905932737*f[62]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[45]*w5Ddx2+0.3535533905932737*f[61]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[31]*w5Ddx2+0.3535533905932737*f[56]*dv5Ddx2; 
  out[57] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[48]*w5Ddx2+(0.3162277660168379*f[145]+0.3535533905932737*f[23])*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[36]*w5Ddx2+(0.3162277660168379*f[138]+0.3535533905932737*f[11])*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[35]*w5Ddx2+(0.3162277660168379*f[137]+0.3535533905932737*f[10])*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[58]*w5Ddx2+(0.3162277660168379*f[154]+0.3535533905932737*f[42])*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[20]*w5Ddx2+(0.3162277660168379*f[132]+0.3535533905932737*f[4])*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[50]*w5Ddx2+(0.3162277660168379*f[147]+0.3535533905932737*f[25])*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[49]*w5Ddx2+(0.3162277660168379*f[146]+0.3535533905932737*f[24])*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[37]*w5Ddx2+(0.3162277660168379*f[139]+0.3535533905932737*f[12])*dv5Ddx2; 
  out[58] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[51]*w5Ddx2+(0.3162277660168379*f[148]+0.3535533905932737*f[26])*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[39]*w5Ddx2+(0.3162277660168379*f[141]+0.3535533905932737*f[14])*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[38]*w5Ddx2+(0.3162277660168379*f[140]+0.3535533905932737*f[13])*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[59]*w5Ddx2+(0.3162277660168379*f[155]+0.3535533905932737*f[43])*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[21]*w5Ddx2+(0.3162277660168379*f[133]+0.3535533905932737*f[5])*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[53]*w5Ddx2+(0.3162277660168379*f[150]+0.3535533905932737*f[28])*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[52]*w5Ddx2+(0.3162277660168379*f[149]+0.3535533905932737*f[27])*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[40]*w5Ddx2+(0.3162277660168379*f[142]+0.3535533905932737*f[15])*dv5Ddx2; 
  out[59] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[54]*w5Ddx2+(0.3162277660168379*f[151]+0.3535533905932737*f[29])*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[41]*w5Ddx2+(0.3162277660168379*f[143]+0.3535533905932737*f[16])*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[60]*w5Ddx2+(0.3162277660168379*f[156]+0.3535533905932737*f[44])*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[61]*w5Ddx2+(0.3162277660168379*f[157]+0.3535533905932737*f[45])*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[55]*w5Ddx2+(0.3162277660168379*f[152]+0.3535533905932737*f[30])*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[56]*w5Ddx2+(0.3162277660168379*f[153]+0.3535533905932737*f[31])*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[63]*w5Ddx2+(0.3162277660168379*f[159]+0.3535533905932737*f[57])*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[62]*w5Ddx2+(0.3162277660168379*f[158]+0.3535533905932737*f[46])*dv5Ddx2; 
  out[61] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[55]*w5Ddx2+(0.3162277660168379*f[152]+0.3535533905932737*f[30])*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[60]*w5Ddx2+(0.3162277660168379*f[156]+0.3535533905932737*f[44])*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[41]*w5Ddx2+(0.3162277660168379*f[143]+0.3535533905932737*f[16])*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[62]*w5Ddx2+(0.3162277660168379*f[158]+0.3535533905932737*f[46])*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[54]*w5Ddx2+(0.3162277660168379*f[151]+0.3535533905932737*f[29])*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[63]*w5Ddx2+(0.3162277660168379*f[159]+0.3535533905932737*f[57])*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[56]*w5Ddx2+(0.3162277660168379*f[153]+0.3535533905932737*f[31])*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[61]*w5Ddx2+(0.3162277660168379*f[157]+0.3535533905932737*f[45])*dv5Ddx2; 
  out[62] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[60]*w5Ddx2+(0.3162277660168379*f[156]+0.3535533905932737*f[44])*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[55]*w5Ddx2+(0.3162277660168379*f[152]+0.3535533905932737*f[30])*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[54]*w5Ddx2+(0.3162277660168379*f[151]+0.3535533905932737*f[29])*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[63]*w5Ddx2+(0.3162277660168379*f[159]+0.3535533905932737*f[57])*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[41]*w5Ddx2+(0.3162277660168379*f[143]+0.3535533905932737*f[16])*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[62]*w5Ddx2+(0.3162277660168379*f[158]+0.3535533905932737*f[46])*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[61]*w5Ddx2+(0.3162277660168379*f[157]+0.3535533905932737*f[45])*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[56]*w5Ddx2+(0.3162277660168379*f[153]+0.3535533905932737*f[31])*dv5Ddx2; 
  out[63] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[64]*w5Ddx2+0.3535533905932737*f[69]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[65]*w5Ddx2+0.3535533905932737*f[76]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[66]*w5Ddx2+0.3535533905932737*f[77]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[67]*w5Ddx2+0.3535533905932737*f[78]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[70]*w5Ddx2+0.3535533905932737*f[84]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[71]*w5Ddx2+0.3535533905932737*f[85]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[72]*w5Ddx2+0.3535533905932737*f[86]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[80]*w5Ddx2+0.3535533905932737*f[91]*dv5Ddx2; 
  out[67] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[65]*w5Ddx2+0.3535533905932737*f[76]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[64]*w5Ddx2+0.3535533905932737*f[69]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[70]*w5Ddx2+0.3535533905932737*f[84]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[71]*w5Ddx2+0.3535533905932737*f[85]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[66]*w5Ddx2+0.3535533905932737*f[77]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[67]*w5Ddx2+0.3535533905932737*f[78]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[80]*w5Ddx2+0.3535533905932737*f[91]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[72]*w5Ddx2+0.3535533905932737*f[86]*dv5Ddx2; 
  out[71] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[66]*w5Ddx2+0.3535533905932737*f[77]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[70]*w5Ddx2+0.3535533905932737*f[84]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[64]*w5Ddx2+0.3535533905932737*f[69]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[72]*w5Ddx2+0.3535533905932737*f[86]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[65]*w5Ddx2+0.3535533905932737*f[76]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[80]*w5Ddx2+0.3535533905932737*f[91]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[67]*w5Ddx2+0.3535533905932737*f[78]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[71]*w5Ddx2+0.3535533905932737*f[85]*dv5Ddx2; 
  out[72] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[68]*w5Ddx2+0.3535533905932737*f[79]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[73]*w5Ddx2+0.3535533905932737*f[87]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[74]*w5Ddx2+0.3535533905932737*f[88]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[75]*w5Ddx2+0.3535533905932737*f[89]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[81]*w5Ddx2+0.3535533905932737*f[92]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[82]*w5Ddx2+0.3535533905932737*f[93]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[83]*w5Ddx2+0.3535533905932737*f[94]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[90]*w5Ddx2+0.3535533905932737*f[95]*dv5Ddx2; 
  out[75] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[69]*w5Ddx2+0.3535533905932737*f[64]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[76]*w5Ddx2+0.3535533905932737*f[65]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[77]*w5Ddx2+0.3535533905932737*f[66]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[78]*w5Ddx2+0.3535533905932737*f[67]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[84]*w5Ddx2+0.3535533905932737*f[70]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[85]*w5Ddx2+0.3535533905932737*f[71]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[86]*w5Ddx2+0.3535533905932737*f[72]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[91]*w5Ddx2+0.3535533905932737*f[80]*dv5Ddx2; 
  out[78] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[70]*w5Ddx2+0.3535533905932737*f[84]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[66]*w5Ddx2+0.3535533905932737*f[77]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[65]*w5Ddx2+0.3535533905932737*f[76]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[80]*w5Ddx2+0.3535533905932737*f[91]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[64]*w5Ddx2+0.3535533905932737*f[69]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[72]*w5Ddx2+0.3535533905932737*f[86]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[71]*w5Ddx2+0.3535533905932737*f[85]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[67]*w5Ddx2+0.3535533905932737*f[78]*dv5Ddx2; 
  out[80] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[73]*w5Ddx2+0.3535533905932737*f[87]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[68]*w5Ddx2+0.3535533905932737*f[79]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[81]*w5Ddx2+0.3535533905932737*f[92]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[82]*w5Ddx2+0.3535533905932737*f[93]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[74]*w5Ddx2+0.3535533905932737*f[88]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[75]*w5Ddx2+0.3535533905932737*f[89]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[90]*w5Ddx2+0.3535533905932737*f[95]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[83]*w5Ddx2+0.3535533905932737*f[94]*dv5Ddx2; 
  out[82] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[74]*w5Ddx2+0.3535533905932737*f[88]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[81]*w5Ddx2+0.3535533905932737*f[92]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[68]*w5Ddx2+0.3535533905932737*f[79]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[83]*w5Ddx2+0.3535533905932737*f[94]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[73]*w5Ddx2+0.3535533905932737*f[87]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[90]*w5Ddx2+0.3535533905932737*f[95]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[75]*w5Ddx2+0.3535533905932737*f[89]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[82]*w5Ddx2+0.3535533905932737*f[93]*dv5Ddx2; 
  out[83] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[76]*w5Ddx2+0.3535533905932737*f[65]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[69]*w5Ddx2+0.3535533905932737*f[64]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[84]*w5Ddx2+0.3535533905932737*f[70]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[85]*w5Ddx2+0.3535533905932737*f[71]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[77]*w5Ddx2+0.3535533905932737*f[66]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[78]*w5Ddx2+0.3535533905932737*f[67]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[91]*w5Ddx2+0.3535533905932737*f[80]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[86]*w5Ddx2+0.3535533905932737*f[72]*dv5Ddx2; 
  out[85] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[77]*w5Ddx2+0.3535533905932737*f[66]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[84]*w5Ddx2+0.3535533905932737*f[70]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[69]*w5Ddx2+0.3535533905932737*f[64]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[86]*w5Ddx2+0.3535533905932737*f[72]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[76]*w5Ddx2+0.3535533905932737*f[65]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[91]*w5Ddx2+0.3535533905932737*f[80]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[78]*w5Ddx2+0.3535533905932737*f[67]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[85]*w5Ddx2+0.3535533905932737*f[71]*dv5Ddx2; 
  out[86] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[79]*w5Ddx2+0.3535533905932737*f[68]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[87]*w5Ddx2+0.3535533905932737*f[73]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[88]*w5Ddx2+0.3535533905932737*f[74]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[89]*w5Ddx2+0.3535533905932737*f[75]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[92]*w5Ddx2+0.3535533905932737*f[81]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[93]*w5Ddx2+0.3535533905932737*f[82]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[94]*w5Ddx2+0.3535533905932737*f[83]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[95]*w5Ddx2+0.3535533905932737*f[90]*dv5Ddx2; 
  out[89] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[81]*w5Ddx2+0.3535533905932737*f[92]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[74]*w5Ddx2+0.3535533905932737*f[88]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[73]*w5Ddx2+0.3535533905932737*f[87]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[90]*w5Ddx2+0.3535533905932737*f[95]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[68]*w5Ddx2+0.3535533905932737*f[79]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[83]*w5Ddx2+0.3535533905932737*f[94]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[82]*w5Ddx2+0.3535533905932737*f[93]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[75]*w5Ddx2+0.3535533905932737*f[89]*dv5Ddx2; 
  out[90] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[84]*w5Ddx2+0.3535533905932737*f[70]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[77]*w5Ddx2+0.3535533905932737*f[66]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[76]*w5Ddx2+0.3535533905932737*f[65]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[91]*w5Ddx2+0.3535533905932737*f[80]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[69]*w5Ddx2+0.3535533905932737*f[64]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[86]*w5Ddx2+0.3535533905932737*f[72]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[85]*w5Ddx2+0.3535533905932737*f[71]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[78]*w5Ddx2+0.3535533905932737*f[67]*dv5Ddx2; 
  out[91] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[87]*w5Ddx2+0.3535533905932737*f[73]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[79]*w5Ddx2+0.3535533905932737*f[68]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[92]*w5Ddx2+0.3535533905932737*f[81]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[93]*w5Ddx2+0.3535533905932737*f[82]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[88]*w5Ddx2+0.3535533905932737*f[74]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[89]*w5Ddx2+0.3535533905932737*f[75]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[95]*w5Ddx2+0.3535533905932737*f[90]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[94]*w5Ddx2+0.3535533905932737*f[83]*dv5Ddx2; 
  out[93] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[88]*w5Ddx2+0.3535533905932737*f[74]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[92]*w5Ddx2+0.3535533905932737*f[81]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[79]*w5Ddx2+0.3535533905932737*f[68]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[94]*w5Ddx2+0.3535533905932737*f[83]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[87]*w5Ddx2+0.3535533905932737*f[73]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[95]*w5Ddx2+0.3535533905932737*f[90]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[89]*w5Ddx2+0.3535533905932737*f[75]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[93]*w5Ddx2+0.3535533905932737*f[82]*dv5Ddx2; 
  out[94] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[92]*w5Ddx2+0.3535533905932737*f[81]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[88]*w5Ddx2+0.3535533905932737*f[74]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[87]*w5Ddx2+0.3535533905932737*f[73]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[95]*w5Ddx2+0.3535533905932737*f[90]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[79]*w5Ddx2+0.3535533905932737*f[68]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[94]*w5Ddx2+0.3535533905932737*f[83]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[93]*w5Ddx2+0.3535533905932737*f[82]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[89]*w5Ddx2+0.3535533905932737*f[75]*dv5Ddx2; 
  out[95] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[96]*w5Ddx2+0.3535533905932737*f[101]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[97]*w5Ddx2+0.3535533905932737*f[108]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[98]*w5Ddx2+0.3535533905932737*f[109]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[99]*w5Ddx2+0.3535533905932737*f[110]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[102]*w5Ddx2+0.3535533905932737*f[116]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[103]*w5Ddx2+0.3535533905932737*f[117]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[104]*w5Ddx2+0.3535533905932737*f[118]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[112]*w5Ddx2+0.3535533905932737*f[123]*dv5Ddx2; 
  out[99] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[97]*w5Ddx2+0.3535533905932737*f[108]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[96]*w5Ddx2+0.3535533905932737*f[101]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[102]*w5Ddx2+0.3535533905932737*f[116]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[103]*w5Ddx2+0.3535533905932737*f[117]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[98]*w5Ddx2+0.3535533905932737*f[109]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[99]*w5Ddx2+0.3535533905932737*f[110]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[112]*w5Ddx2+0.3535533905932737*f[123]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[104]*w5Ddx2+0.3535533905932737*f[118]*dv5Ddx2; 
  out[103] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[98]*w5Ddx2+0.3535533905932737*f[109]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[102]*w5Ddx2+0.3535533905932737*f[116]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[96]*w5Ddx2+0.3535533905932737*f[101]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[104]*w5Ddx2+0.3535533905932737*f[118]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[97]*w5Ddx2+0.3535533905932737*f[108]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[112]*w5Ddx2+0.3535533905932737*f[123]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[99]*w5Ddx2+0.3535533905932737*f[110]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[103]*w5Ddx2+0.3535533905932737*f[117]*dv5Ddx2; 
  out[104] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[100]*w5Ddx2+0.3535533905932737*f[111]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[105]*w5Ddx2+0.3535533905932737*f[119]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[106]*w5Ddx2+0.3535533905932737*f[120]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[107]*w5Ddx2+0.3535533905932737*f[121]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[113]*w5Ddx2+0.3535533905932737*f[124]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[114]*w5Ddx2+0.3535533905932737*f[125]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[115]*w5Ddx2+0.3535533905932737*f[126]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[122]*w5Ddx2+0.3535533905932737*f[127]*dv5Ddx2; 
  out[107] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[101]*w5Ddx2+0.3535533905932737*f[96]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[108]*w5Ddx2+0.3535533905932737*f[97]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[109]*w5Ddx2+0.3535533905932737*f[98]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[110]*w5Ddx2+0.3535533905932737*f[99]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[116]*w5Ddx2+0.3535533905932737*f[102]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[117]*w5Ddx2+0.3535533905932737*f[103]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[118]*w5Ddx2+0.3535533905932737*f[104]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[123]*w5Ddx2+0.3535533905932737*f[112]*dv5Ddx2; 
  out[110] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[102]*w5Ddx2+0.3535533905932737*f[116]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[98]*w5Ddx2+0.3535533905932737*f[109]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[97]*w5Ddx2+0.3535533905932737*f[108]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[112]*w5Ddx2+0.3535533905932737*f[123]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[96]*w5Ddx2+0.3535533905932737*f[101]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[104]*w5Ddx2+0.3535533905932737*f[118]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[103]*w5Ddx2+0.3535533905932737*f[117]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[99]*w5Ddx2+0.3535533905932737*f[110]*dv5Ddx2; 
  out[112] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[105]*w5Ddx2+0.3535533905932737*f[119]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[100]*w5Ddx2+0.3535533905932737*f[111]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[113]*w5Ddx2+0.3535533905932737*f[124]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[114]*w5Ddx2+0.3535533905932737*f[125]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[106]*w5Ddx2+0.3535533905932737*f[120]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[107]*w5Ddx2+0.3535533905932737*f[121]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[122]*w5Ddx2+0.3535533905932737*f[127]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[115]*w5Ddx2+0.3535533905932737*f[126]*dv5Ddx2; 
  out[114] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[106]*w5Ddx2+0.3535533905932737*f[120]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[113]*w5Ddx2+0.3535533905932737*f[124]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[100]*w5Ddx2+0.3535533905932737*f[111]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[115]*w5Ddx2+0.3535533905932737*f[126]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[105]*w5Ddx2+0.3535533905932737*f[119]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[122]*w5Ddx2+0.3535533905932737*f[127]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[107]*w5Ddx2+0.3535533905932737*f[121]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[114]*w5Ddx2+0.3535533905932737*f[125]*dv5Ddx2; 
  out[115] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[108]*w5Ddx2+0.3535533905932737*f[97]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[101]*w5Ddx2+0.3535533905932737*f[96]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[116]*w5Ddx2+0.3535533905932737*f[102]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[117]*w5Ddx2+0.3535533905932737*f[103]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[109]*w5Ddx2+0.3535533905932737*f[98]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[110]*w5Ddx2+0.3535533905932737*f[99]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[123]*w5Ddx2+0.3535533905932737*f[112]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[118]*w5Ddx2+0.3535533905932737*f[104]*dv5Ddx2; 
  out[117] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[109]*w5Ddx2+0.3535533905932737*f[98]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[116]*w5Ddx2+0.3535533905932737*f[102]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[101]*w5Ddx2+0.3535533905932737*f[96]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[118]*w5Ddx2+0.3535533905932737*f[104]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[108]*w5Ddx2+0.3535533905932737*f[97]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[123]*w5Ddx2+0.3535533905932737*f[112]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[110]*w5Ddx2+0.3535533905932737*f[99]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[117]*w5Ddx2+0.3535533905932737*f[103]*dv5Ddx2; 
  out[118] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[111]*w5Ddx2+0.3535533905932737*f[100]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[119]*w5Ddx2+0.3535533905932737*f[105]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[120]*w5Ddx2+0.3535533905932737*f[106]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[121]*w5Ddx2+0.3535533905932737*f[107]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[124]*w5Ddx2+0.3535533905932737*f[113]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[125]*w5Ddx2+0.3535533905932737*f[114]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[126]*w5Ddx2+0.3535533905932737*f[115]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[127]*w5Ddx2+0.3535533905932737*f[122]*dv5Ddx2; 
  out[121] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[113]*w5Ddx2+0.3535533905932737*f[124]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[106]*w5Ddx2+0.3535533905932737*f[120]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[105]*w5Ddx2+0.3535533905932737*f[119]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[122]*w5Ddx2+0.3535533905932737*f[127]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[100]*w5Ddx2+0.3535533905932737*f[111]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[115]*w5Ddx2+0.3535533905932737*f[126]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[114]*w5Ddx2+0.3535533905932737*f[125]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[107]*w5Ddx2+0.3535533905932737*f[121]*dv5Ddx2; 
  out[122] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[116]*w5Ddx2+0.3535533905932737*f[102]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[109]*w5Ddx2+0.3535533905932737*f[98]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[108]*w5Ddx2+0.3535533905932737*f[97]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[123]*w5Ddx2+0.3535533905932737*f[112]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[101]*w5Ddx2+0.3535533905932737*f[96]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[118]*w5Ddx2+0.3535533905932737*f[104]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[117]*w5Ddx2+0.3535533905932737*f[103]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[110]*w5Ddx2+0.3535533905932737*f[99]*dv5Ddx2; 
  out[123] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[119]*w5Ddx2+0.3535533905932737*f[105]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[111]*w5Ddx2+0.3535533905932737*f[100]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[124]*w5Ddx2+0.3535533905932737*f[113]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[125]*w5Ddx2+0.3535533905932737*f[114]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[120]*w5Ddx2+0.3535533905932737*f[106]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[121]*w5Ddx2+0.3535533905932737*f[107]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[127]*w5Ddx2+0.3535533905932737*f[122]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[126]*w5Ddx2+0.3535533905932737*f[115]*dv5Ddx2; 
  out[125] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[120]*w5Ddx2+0.3535533905932737*f[106]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[124]*w5Ddx2+0.3535533905932737*f[113]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[111]*w5Ddx2+0.3535533905932737*f[100]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[126]*w5Ddx2+0.3535533905932737*f[115]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[119]*w5Ddx2+0.3535533905932737*f[105]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[127]*w5Ddx2+0.3535533905932737*f[122]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[121]*w5Ddx2+0.3535533905932737*f[107]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[125]*w5Ddx2+0.3535533905932737*f[114]*dv5Ddx2; 
  out[126] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[124]*w5Ddx2+0.3535533905932737*f[113]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[120]*w5Ddx2+0.3535533905932737*f[106]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[119]*w5Ddx2+0.3535533905932737*f[105]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[127]*w5Ddx2+0.3535533905932737*f[122]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[111]*w5Ddx2+0.3535533905932737*f[100]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[126]*w5Ddx2+0.3535533905932737*f[115]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[125]*w5Ddx2+0.3535533905932737*f[114]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[121]*w5Ddx2+0.3535533905932737*f[107]*dv5Ddx2; 
  out[127] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[128]*w5Ddx2+0.3162277660168379*f[6]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[129]*w5Ddx2+0.3162277660168379*f[17]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[130]*w5Ddx2+0.3162277660168379*f[18]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[131]*w5Ddx2+0.3162277660168379*f[19]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[134]*w5Ddx2+0.3162277660168379*f[32]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[135]*w5Ddx2+0.3162277660168379*f[33]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[136]*w5Ddx2+0.3162277660168379*f[34]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[144]*w5Ddx2+0.3162277660168379*f[47]*dv5Ddx2; 
  out[131] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[129]*w5Ddx2+0.3162277660168379*f[17]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[128]*w5Ddx2+0.3162277660168379*f[6]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[134]*w5Ddx2+0.3162277660168379*f[32]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[135]*w5Ddx2+0.3162277660168379*f[33]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[130]*w5Ddx2+0.3162277660168379*f[18]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[131]*w5Ddx2+0.3162277660168379*f[19]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[144]*w5Ddx2+0.3162277660168379*f[47]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[136]*w5Ddx2+0.3162277660168379*f[34]*dv5Ddx2; 
  out[135] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[130]*w5Ddx2+0.3162277660168379*f[18]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[134]*w5Ddx2+0.3162277660168379*f[32]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[128]*w5Ddx2+0.3162277660168379*f[6]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[136]*w5Ddx2+0.3162277660168379*f[34]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[129]*w5Ddx2+0.3162277660168379*f[17]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[144]*w5Ddx2+0.3162277660168379*f[47]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[131]*w5Ddx2+0.3162277660168379*f[19]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[135]*w5Ddx2+0.3162277660168379*f[33]*dv5Ddx2; 
  out[136] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[132]*w5Ddx2+0.3162277660168379*f[20]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[137]*w5Ddx2+0.3162277660168379*f[35]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[138]*w5Ddx2+0.3162277660168379*f[36]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[139]*w5Ddx2+0.3162277660168379*f[37]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[145]*w5Ddx2+0.3162277660168379*f[48]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[146]*w5Ddx2+0.3162277660168379*f[49]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[147]*w5Ddx2+0.3162277660168379*f[50]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[154]*w5Ddx2+0.3162277660168379*f[58]*dv5Ddx2; 
  out[139] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[133]*w5Ddx2+0.3162277660168379*f[21]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[140]*w5Ddx2+0.3162277660168379*f[38]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[141]*w5Ddx2+0.3162277660168379*f[39]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[142]*w5Ddx2+0.3162277660168379*f[40]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[148]*w5Ddx2+0.3162277660168379*f[51]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[149]*w5Ddx2+0.3162277660168379*f[52]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[150]*w5Ddx2+0.3162277660168379*f[53]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[155]*w5Ddx2+0.3162277660168379*f[59]*dv5Ddx2; 
  out[142] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[134]*w5Ddx2+0.3162277660168379*f[32]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[130]*w5Ddx2+0.3162277660168379*f[18]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[129]*w5Ddx2+0.3162277660168379*f[17]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[144]*w5Ddx2+0.3162277660168379*f[47]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[128]*w5Ddx2+0.3162277660168379*f[6]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[136]*w5Ddx2+0.3162277660168379*f[34]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[135]*w5Ddx2+0.3162277660168379*f[33]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[131]*w5Ddx2+0.3162277660168379*f[19]*dv5Ddx2; 
  out[144] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[137]*w5Ddx2+0.3162277660168379*f[35]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[132]*w5Ddx2+0.3162277660168379*f[20]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[145]*w5Ddx2+0.3162277660168379*f[48]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[146]*w5Ddx2+0.3162277660168379*f[49]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[138]*w5Ddx2+0.3162277660168379*f[36]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[139]*w5Ddx2+0.3162277660168379*f[37]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[154]*w5Ddx2+0.3162277660168379*f[58]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[147]*w5Ddx2+0.3162277660168379*f[50]*dv5Ddx2; 
  out[146] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[138]*w5Ddx2+0.3162277660168379*f[36]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[145]*w5Ddx2+0.3162277660168379*f[48]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[132]*w5Ddx2+0.3162277660168379*f[20]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[147]*w5Ddx2+0.3162277660168379*f[50]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[137]*w5Ddx2+0.3162277660168379*f[35]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[154]*w5Ddx2+0.3162277660168379*f[58]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[139]*w5Ddx2+0.3162277660168379*f[37]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[146]*w5Ddx2+0.3162277660168379*f[49]*dv5Ddx2; 
  out[147] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[140]*w5Ddx2+0.3162277660168379*f[38]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[133]*w5Ddx2+0.3162277660168379*f[21]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[148]*w5Ddx2+0.3162277660168379*f[51]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[149]*w5Ddx2+0.3162277660168379*f[52]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[141]*w5Ddx2+0.3162277660168379*f[39]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[142]*w5Ddx2+0.3162277660168379*f[40]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[155]*w5Ddx2+0.3162277660168379*f[59]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[150]*w5Ddx2+0.3162277660168379*f[53]*dv5Ddx2; 
  out[149] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[141]*w5Ddx2+0.3162277660168379*f[39]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[148]*w5Ddx2+0.3162277660168379*f[51]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[133]*w5Ddx2+0.3162277660168379*f[21]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[150]*w5Ddx2+0.3162277660168379*f[53]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[140]*w5Ddx2+0.3162277660168379*f[38]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[155]*w5Ddx2+0.3162277660168379*f[59]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[142]*w5Ddx2+0.3162277660168379*f[40]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[149]*w5Ddx2+0.3162277660168379*f[52]*dv5Ddx2; 
  out[150] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[143]*w5Ddx2+0.3162277660168379*f[41]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[151]*w5Ddx2+0.3162277660168379*f[54]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[152]*w5Ddx2+0.3162277660168379*f[55]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[153]*w5Ddx2+0.3162277660168379*f[56]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[156]*w5Ddx2+0.3162277660168379*f[60]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[157]*w5Ddx2+0.3162277660168379*f[61]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[158]*w5Ddx2+0.3162277660168379*f[62]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[159]*w5Ddx2+0.3162277660168379*f[63]*dv5Ddx2; 
  out[153] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[145]*w5Ddx2+0.3162277660168379*f[48]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[138]*w5Ddx2+0.3162277660168379*f[36]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[137]*w5Ddx2+0.3162277660168379*f[35]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[154]*w5Ddx2+0.3162277660168379*f[58]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[132]*w5Ddx2+0.3162277660168379*f[20]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[147]*w5Ddx2+0.3162277660168379*f[50]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[146]*w5Ddx2+0.3162277660168379*f[49]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[139]*w5Ddx2+0.3162277660168379*f[37]*dv5Ddx2; 
  out[154] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[148]*w5Ddx2+0.3162277660168379*f[51]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[141]*w5Ddx2+0.3162277660168379*f[39]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[140]*w5Ddx2+0.3162277660168379*f[38]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[155]*w5Ddx2+0.3162277660168379*f[59]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[133]*w5Ddx2+0.3162277660168379*f[21]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[150]*w5Ddx2+0.3162277660168379*f[53]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[149]*w5Ddx2+0.3162277660168379*f[52]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[142]*w5Ddx2+0.3162277660168379*f[40]*dv5Ddx2; 
  out[155] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[151]*w5Ddx2+0.3162277660168379*f[54]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[143]*w5Ddx2+0.3162277660168379*f[41]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[156]*w5Ddx2+0.3162277660168379*f[60]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[157]*w5Ddx2+0.3162277660168379*f[61]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[152]*w5Ddx2+0.3162277660168379*f[55]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[153]*w5Ddx2+0.3162277660168379*f[56]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[159]*w5Ddx2+0.3162277660168379*f[63]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[158]*w5Ddx2+0.3162277660168379*f[62]*dv5Ddx2; 
  out[157] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[152]*w5Ddx2+0.3162277660168379*f[55]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[156]*w5Ddx2+0.3162277660168379*f[60]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[143]*w5Ddx2+0.3162277660168379*f[41]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[158]*w5Ddx2+0.3162277660168379*f[62]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[151]*w5Ddx2+0.3162277660168379*f[54]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[159]*w5Ddx2+0.3162277660168379*f[63]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[153]*w5Ddx2+0.3162277660168379*f[56]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[157]*w5Ddx2+0.3162277660168379*f[61]*dv5Ddx2; 
  out[158] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 

  Gbar[0] = 1.224744871391589*f[156]*w5Ddx2+0.3162277660168379*f[60]*dv5Ddx2; 
  Gbar[1] = 1.224744871391589*f[152]*w5Ddx2+0.3162277660168379*f[55]*dv5Ddx2; 
  Gbar[2] = 1.224744871391589*f[151]*w5Ddx2+0.3162277660168379*f[54]*dv5Ddx2; 
  Gbar[3] = 1.224744871391589*f[159]*w5Ddx2+0.3162277660168379*f[63]*dv5Ddx2; 
  Gbar[4] = 1.224744871391589*f[143]*w5Ddx2+0.3162277660168379*f[41]*dv5Ddx2; 
  Gbar[5] = 1.224744871391589*f[158]*w5Ddx2+0.3162277660168379*f[62]*dv5Ddx2; 
  Gbar[6] = 1.224744871391589*f[157]*w5Ddx2+0.3162277660168379*f[61]*dv5Ddx2; 
  Gbar[7] = 1.224744871391589*f[153]*w5Ddx2+0.3162277660168379*f[56]*dv5Ddx2; 
  out[159] += 0.25*((2.0*Gbar[7]+1.732050807568877*Gbar[4])*basisVecComp[7]+1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]+1.732050807568877*Gbar[2])*basisVecComp[6]+1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]+1.732050807568877*Gbar[1])*basisVecComp[5]+1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]+1.732050807568877*Gbar[0])*basisVecComp[3]+1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0])+0.25*((2.0*Gbar[7]-1.732050807568877*Gbar[4])*basisVecComp[7]-1.732050807568877*basisVecComp[4]*Gbar[7]+(2.0*Gbar[6]-1.732050807568877*Gbar[2])*basisVecComp[6]-1.732050807568877*basisVecComp[2]*Gbar[6]+(2.0*Gbar[5]-1.732050807568877*Gbar[1])*basisVecComp[5]-1.732050807568877*basisVecComp[1]*Gbar[5]+2.0*Gbar[4]*basisVecComp[4]+(2.0*Gbar[3]-1.732050807568877*Gbar[0])*basisVecComp[3]-1.732050807568877*basisVecComp[0]*Gbar[3]+2.0*Gbar[2]*basisVecComp[2]+2.0*Gbar[1]*basisVecComp[1]+2.0*Gbar[0]*basisVecComp[0]); 


  return 3.0*(fabs(w3Ddx0)+0.5*dv3Ddx0+fabs(w4Ddx1)+0.5*dv4Ddx1+fabs(w5Ddx2)+0.5*dv5Ddx2);
} 
