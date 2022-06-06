#include <gkyl_vlasov_bflux_kernels.h> 
GKYL_CU_DH void vlasov_bflux_1x3v_ser_p1(const int *idx, const double *w, enum gkyl_vel_edge edge, const double *dxv, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fIn:  Input Distribution function in cell.
  // out:       Incremented distribution function in center cell.
 if (edge == GKYL_VX_UPPER) {

  const double dx10 = 2/dxv[0]; 
  const double dv = dxv[1], wv = w[1]; 
  double Ghat[8]; 
  Ghat[0] = (1.224744871391589*fIn[1]+0.7071067811865475*fIn[0])*wv+(0.3535533905932737*fIn[5]+0.2041241452319314*fIn[2])*dv; 
  Ghat[1] = (1.224744871391589*fIn[5]+0.7071067811865475*fIn[2])*wv+(0.3535533905932737*fIn[1]+0.2041241452319314*fIn[0])*dv; 
  Ghat[2] = (1.224744871391589*fIn[6]+0.7071067811865475*fIn[3])*wv+(0.3535533905932737*fIn[11]+0.2041241452319314*fIn[7])*dv; 
  Ghat[3] = (1.224744871391589*fIn[8]+0.7071067811865475*fIn[4])*wv+(0.3535533905932737*fIn[12]+0.2041241452319314*fIn[9])*dv; 
  Ghat[4] = (1.224744871391589*fIn[11]+0.7071067811865475*fIn[7])*wv+(0.3535533905932737*fIn[6]+0.2041241452319314*fIn[3])*dv; 
  Ghat[5] = (1.224744871391589*fIn[12]+0.7071067811865475*fIn[9])*wv+(0.3535533905932737*fIn[8]+0.2041241452319314*fIn[4])*dv; 
  Ghat[6] = (1.224744871391589*fIn[13]+0.7071067811865475*fIn[10])*wv+(0.3535533905932737*fIn[15]+0.2041241452319314*fIn[14])*dv; 
  Ghat[7] = (1.224744871391589*fIn[15]+0.7071067811865475*fIn[14])*wv+(0.3535533905932737*fIn[13]+0.2041241452319314*fIn[10])*dv; 

  out[0] += -0.7071067811865475*Ghat[0]*dx10; 
  out[1] += -1.224744871391589*Ghat[0]*dx10; 
  out[2] += -0.7071067811865475*Ghat[1]*dx10; 
  out[3] += -0.7071067811865475*Ghat[2]*dx10; 
  out[4] += -0.7071067811865475*Ghat[3]*dx10; 
  out[5] += -1.224744871391589*Ghat[1]*dx10; 
  out[6] += -1.224744871391589*Ghat[2]*dx10; 
  out[7] += -0.7071067811865475*Ghat[4]*dx10; 
  out[8] += -1.224744871391589*Ghat[3]*dx10; 
  out[9] += -0.7071067811865475*Ghat[5]*dx10; 
  out[10] += -0.7071067811865475*Ghat[6]*dx10; 
  out[11] += -1.224744871391589*Ghat[4]*dx10; 
  out[12] += -1.224744871391589*Ghat[5]*dx10; 
  out[13] += -1.224744871391589*Ghat[6]*dx10; 
  out[14] += -0.7071067811865475*Ghat[7]*dx10; 
  out[15] += -1.224744871391589*Ghat[7]*dx10; 

 } else if (edge == GKYL_VX_LOWER) {

  const double dx10 = 2/dxv[0]; 
  const double dv = dxv[1], wv = w[1]; 
  double Ghat[8]; 
  Ghat[0] = -0.1178511301977578*((10.39230484541326*fIn[1]-6.0*fIn[0])*wv+(3.0*fIn[5]-1.732050807568877*fIn[2])*dv); 
  Ghat[1] = -0.1178511301977578*((10.39230484541326*fIn[5]-6.0*fIn[2])*wv+(3.0*fIn[1]-1.732050807568877*fIn[0])*dv); 
  Ghat[2] = -0.1178511301977578*((10.39230484541326*fIn[6]-6.0*fIn[3])*wv+(3.0*fIn[11]-1.732050807568877*fIn[7])*dv); 
  Ghat[3] = -0.1178511301977578*((10.39230484541326*fIn[8]-6.0*fIn[4])*wv+(3.0*fIn[12]-1.732050807568877*fIn[9])*dv); 
  Ghat[4] = -0.1178511301977578*((10.39230484541326*fIn[11]-6.0*fIn[7])*wv+(3.0*fIn[6]-1.732050807568877*fIn[3])*dv); 
  Ghat[5] = -0.1178511301977578*((10.39230484541326*fIn[12]-6.0*fIn[9])*wv+(3.0*fIn[8]-1.732050807568877*fIn[4])*dv); 
  Ghat[6] = -0.1178511301977578*((10.39230484541326*fIn[13]-6.0*fIn[10])*wv+(3.0*fIn[15]-1.732050807568877*fIn[14])*dv); 
  Ghat[7] = -0.1178511301977578*((10.39230484541326*fIn[15]-6.0*fIn[14])*wv+(3.0*fIn[13]-1.732050807568877*fIn[10])*dv); 

  out[0] += 0.7071067811865475*Ghat[0]*dx10; 
  out[1] += -1.224744871391589*Ghat[0]*dx10; 
  out[2] += 0.7071067811865475*Ghat[1]*dx10; 
  out[3] += 0.7071067811865475*Ghat[2]*dx10; 
  out[4] += 0.7071067811865475*Ghat[3]*dx10; 
  out[5] += -1.224744871391589*Ghat[1]*dx10; 
  out[6] += -1.224744871391589*Ghat[2]*dx10; 
  out[7] += 0.7071067811865475*Ghat[4]*dx10; 
  out[8] += -1.224744871391589*Ghat[3]*dx10; 
  out[9] += 0.7071067811865475*Ghat[5]*dx10; 
  out[10] += 0.7071067811865475*Ghat[6]*dx10; 
  out[11] += -1.224744871391589*Ghat[4]*dx10; 
  out[12] += -1.224744871391589*Ghat[5]*dx10; 
  out[13] += -1.224744871391589*Ghat[6]*dx10; 
  out[14] += 0.7071067811865475*Ghat[7]*dx10; 
  out[15] += -1.224744871391589*Ghat[7]*dx10; 
  } 

} 

