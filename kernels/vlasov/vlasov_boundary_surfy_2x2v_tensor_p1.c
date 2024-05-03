#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_boundary_surfy_2x2v_tensor_p1(const double *w, const double *dxv, const double *alpha_geo, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:     Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // alpha_geo:   Fields used only for general geometry.
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Incremented distribution function in center cell.
  const double dx11 = 2/dxv[1]; 
  const double dv = dxv[3], wv = w[3]; 
  double Ghat[16]; 

  if (edge == -1) { 

  if (wv>0) { 

  Ghat[0] = (1.224744871391589*fSkin[2]+0.7071067811865475*fSkin[0])*wv+(0.3535533905932737*fSkin[9]+0.2041241452319315*fSkin[4])*dv; 
  Ghat[1] = (1.224744871391589*fSkin[5]+0.7071067811865475*fSkin[1])*wv+(0.3535533905932737*fSkin[12]+0.2041241452319315*fSkin[8])*dv; 
  Ghat[2] = (1.224744871391589*fSkin[7]+0.7071067811865475*fSkin[3])*wv+(0.3535533905932737*fSkin[14]+0.2041241452319315*fSkin[10])*dv; 
  Ghat[3] = (1.224744871391589*fSkin[9]+0.7071067811865475*fSkin[4])*wv+(0.3535533905932737*fSkin[2]+0.2041241452319315*fSkin[0])*dv; 
  Ghat[4] = (1.224744871391589*fSkin[11]+0.7071067811865475*fSkin[6])*wv+(0.3535533905932737*fSkin[15]+0.2041241452319315*fSkin[13])*dv; 
  Ghat[5] = (1.224744871391589*fSkin[12]+0.7071067811865475*fSkin[8])*wv+(0.3535533905932737*fSkin[5]+0.2041241452319315*fSkin[1])*dv; 
  Ghat[6] = (1.224744871391589*fSkin[14]+0.7071067811865475*fSkin[10])*wv+(0.3535533905932737*fSkin[7]+0.2041241452319315*fSkin[3])*dv; 
  Ghat[7] = (1.224744871391589*fSkin[15]+0.7071067811865475*fSkin[13])*wv+(0.3535533905932737*fSkin[11]+0.2041241452319315*fSkin[6])*dv; 
  Ghat[12] = (0.3162277660168379*fSkin[9]+0.1825741858350554*fSkin[4])*dv; 
  Ghat[13] = (0.3162277660168379*fSkin[12]+0.1825741858350554*fSkin[8])*dv; 
  Ghat[14] = (0.3162277660168379*fSkin[14]+0.1825741858350554*fSkin[10])*dv; 
  Ghat[15] = (0.3162277660168379*fSkin[15]+0.1825741858350554*fSkin[13])*dv; 

  } else { 

  Ghat[0] = -0.1178511301977579*((10.39230484541326*fEdge[2]-6.0*fEdge[0])*wv+(3.0*fEdge[9]-1.732050807568877*fEdge[4])*dv); 
  Ghat[1] = -0.1178511301977579*((10.39230484541326*fEdge[5]-6.0*fEdge[1])*wv+(3.0*fEdge[12]-1.732050807568877*fEdge[8])*dv); 
  Ghat[2] = -0.1178511301977579*((10.39230484541326*fEdge[7]-6.0*fEdge[3])*wv+(3.0*fEdge[14]-1.732050807568877*fEdge[10])*dv); 
  Ghat[3] = -0.1178511301977579*((10.39230484541326*fEdge[9]-6.0*fEdge[4])*wv+(3.0*fEdge[2]-1.732050807568877*fEdge[0])*dv); 
  Ghat[4] = -0.1178511301977579*((10.39230484541326*fEdge[11]-6.0*fEdge[6])*wv+(3.0*fEdge[15]-1.732050807568877*fEdge[13])*dv); 
  Ghat[5] = -0.1178511301977579*((10.39230484541326*fEdge[12]-6.0*fEdge[8])*wv+(3.0*fEdge[5]-1.732050807568877*fEdge[1])*dv); 
  Ghat[6] = -0.1178511301977579*((10.39230484541326*fEdge[14]-6.0*fEdge[10])*wv+(3.0*fEdge[7]-1.732050807568877*fEdge[3])*dv); 
  Ghat[7] = -0.1178511301977579*((10.39230484541326*fEdge[15]-6.0*fEdge[13])*wv+(3.0*fEdge[11]-1.732050807568877*fEdge[6])*dv); 
  Ghat[12] = -0.04714045207910316*(6.708203932499369*fEdge[9]-3.872983346207417*fEdge[4])*dv; 
  Ghat[13] = -0.04714045207910316*(6.708203932499369*fEdge[12]-3.872983346207417*fEdge[8])*dv; 
  Ghat[14] = -0.04714045207910316*(6.708203932499369*fEdge[14]-3.872983346207417*fEdge[10])*dv; 
  Ghat[15] = -0.04714045207910316*(6.708203932499369*fEdge[15]-3.872983346207417*fEdge[13])*dv; 

  } 

  out[0] += -0.7071067811865475*Ghat[0]*dx11; 
  out[1] += -0.7071067811865475*Ghat[1]*dx11; 
  out[2] += -1.224744871391589*Ghat[0]*dx11; 
  out[3] += -0.7071067811865475*Ghat[2]*dx11; 
  out[4] += -0.7071067811865475*Ghat[3]*dx11; 
  out[5] += -1.224744871391589*Ghat[1]*dx11; 
  out[6] += -0.7071067811865475*Ghat[4]*dx11; 
  out[7] += -1.224744871391589*Ghat[2]*dx11; 
  out[8] += -0.7071067811865475*Ghat[5]*dx11; 
  out[9] += -1.224744871391589*Ghat[3]*dx11; 
  out[10] += -0.7071067811865475*Ghat[6]*dx11; 
  out[11] += -1.224744871391589*Ghat[4]*dx11; 
  out[12] += -1.224744871391589*Ghat[5]*dx11; 
  out[13] += -0.7071067811865475*Ghat[7]*dx11; 
  out[14] += -1.224744871391589*Ghat[6]*dx11; 
  out[15] += -1.224744871391589*Ghat[7]*dx11; 

  } else { 

  if (wv>0) { 

  Ghat[0] = (1.224744871391589*fEdge[2]+0.7071067811865475*fEdge[0])*wv+(0.3535533905932737*fEdge[9]+0.2041241452319315*fEdge[4])*dv; 
  Ghat[1] = (1.224744871391589*fEdge[5]+0.7071067811865475*fEdge[1])*wv+(0.3535533905932737*fEdge[12]+0.2041241452319315*fEdge[8])*dv; 
  Ghat[2] = (1.224744871391589*fEdge[7]+0.7071067811865475*fEdge[3])*wv+(0.3535533905932737*fEdge[14]+0.2041241452319315*fEdge[10])*dv; 
  Ghat[3] = (1.224744871391589*fEdge[9]+0.7071067811865475*fEdge[4])*wv+(0.3535533905932737*fEdge[2]+0.2041241452319315*fEdge[0])*dv; 
  Ghat[4] = (1.224744871391589*fEdge[11]+0.7071067811865475*fEdge[6])*wv+(0.3535533905932737*fEdge[15]+0.2041241452319315*fEdge[13])*dv; 
  Ghat[5] = (1.224744871391589*fEdge[12]+0.7071067811865475*fEdge[8])*wv+(0.3535533905932737*fEdge[5]+0.2041241452319315*fEdge[1])*dv; 
  Ghat[6] = (1.224744871391589*fEdge[14]+0.7071067811865475*fEdge[10])*wv+(0.3535533905932737*fEdge[7]+0.2041241452319315*fEdge[3])*dv; 
  Ghat[7] = (1.224744871391589*fEdge[15]+0.7071067811865475*fEdge[13])*wv+(0.3535533905932737*fEdge[11]+0.2041241452319315*fEdge[6])*dv; 
  Ghat[12] = (0.3162277660168379*fEdge[9]+0.1825741858350554*fEdge[4])*dv; 
  Ghat[13] = (0.3162277660168379*fEdge[12]+0.1825741858350554*fEdge[8])*dv; 
  Ghat[14] = (0.3162277660168379*fEdge[14]+0.1825741858350554*fEdge[10])*dv; 
  Ghat[15] = (0.3162277660168379*fEdge[15]+0.1825741858350554*fEdge[13])*dv; 

  } else { 

  Ghat[0] = -0.1178511301977579*((10.39230484541326*fSkin[2]-6.0*fSkin[0])*wv+(3.0*fSkin[9]-1.732050807568877*fSkin[4])*dv); 
  Ghat[1] = -0.1178511301977579*((10.39230484541326*fSkin[5]-6.0*fSkin[1])*wv+(3.0*fSkin[12]-1.732050807568877*fSkin[8])*dv); 
  Ghat[2] = -0.1178511301977579*((10.39230484541326*fSkin[7]-6.0*fSkin[3])*wv+(3.0*fSkin[14]-1.732050807568877*fSkin[10])*dv); 
  Ghat[3] = -0.1178511301977579*((10.39230484541326*fSkin[9]-6.0*fSkin[4])*wv+(3.0*fSkin[2]-1.732050807568877*fSkin[0])*dv); 
  Ghat[4] = -0.1178511301977579*((10.39230484541326*fSkin[11]-6.0*fSkin[6])*wv+(3.0*fSkin[15]-1.732050807568877*fSkin[13])*dv); 
  Ghat[5] = -0.1178511301977579*((10.39230484541326*fSkin[12]-6.0*fSkin[8])*wv+(3.0*fSkin[5]-1.732050807568877*fSkin[1])*dv); 
  Ghat[6] = -0.1178511301977579*((10.39230484541326*fSkin[14]-6.0*fSkin[10])*wv+(3.0*fSkin[7]-1.732050807568877*fSkin[3])*dv); 
  Ghat[7] = -0.1178511301977579*((10.39230484541326*fSkin[15]-6.0*fSkin[13])*wv+(3.0*fSkin[11]-1.732050807568877*fSkin[6])*dv); 
  Ghat[12] = -0.04714045207910316*(6.708203932499369*fSkin[9]-3.872983346207417*fSkin[4])*dv; 
  Ghat[13] = -0.04714045207910316*(6.708203932499369*fSkin[12]-3.872983346207417*fSkin[8])*dv; 
  Ghat[14] = -0.04714045207910316*(6.708203932499369*fSkin[14]-3.872983346207417*fSkin[10])*dv; 
  Ghat[15] = -0.04714045207910316*(6.708203932499369*fSkin[15]-3.872983346207417*fSkin[13])*dv; 

  } 

  out[0] += 0.7071067811865475*Ghat[0]*dx11; 
  out[1] += 0.7071067811865475*Ghat[1]*dx11; 
  out[2] += -1.224744871391589*Ghat[0]*dx11; 
  out[3] += 0.7071067811865475*Ghat[2]*dx11; 
  out[4] += 0.7071067811865475*Ghat[3]*dx11; 
  out[5] += -1.224744871391589*Ghat[1]*dx11; 
  out[6] += 0.7071067811865475*Ghat[4]*dx11; 
  out[7] += -1.224744871391589*Ghat[2]*dx11; 
  out[8] += 0.7071067811865475*Ghat[5]*dx11; 
  out[9] += -1.224744871391589*Ghat[3]*dx11; 
  out[10] += 0.7071067811865475*Ghat[6]*dx11; 
  out[11] += -1.224744871391589*Ghat[4]*dx11; 
  out[12] += -1.224744871391589*Ghat[5]*dx11; 
  out[13] += 0.7071067811865475*Ghat[7]*dx11; 
  out[14] += -1.224744871391589*Ghat[6]*dx11; 
  out[15] += -1.224744871391589*Ghat[7]*dx11; 

  } 
  return 0.;

} 