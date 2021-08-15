#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_boundary_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *qmem, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // qmem:        q/m*EM fields.
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv10 = 2/dxv[1]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double *E0 = &qmem[0]; 

  double Ghat[2]; 
  double alpha[2]; 

  alpha[0] = E0[0]; 
  alpha[1] = E0[1]; 

  double fUpwindQuad[2];
  double fUpwind[2];

  if (edge == -1) { 

  fUpwindQuad[0] = copysignf(1.0,alpha[0]-alpha[1])*((-0.8660254037844386*(fSkin[3]+fEdge[3]))+0.8660254037844386*(fSkin[2]+fEdge[2])-0.5*fSkin[1]+0.5*(fEdge[1]+fSkin[0])-0.5*fEdge[0])-0.8660254037844386*fSkin[3]+0.8660254037844386*(fEdge[3]+fSkin[2])-0.8660254037844386*fEdge[2]-0.5*(fSkin[1]+fEdge[1])+0.5*(fSkin[0]+fEdge[0]); 
  fUpwindQuad[1] = copysignf(1.0,alpha[1]+alpha[0])*(0.8660254037844386*(fSkin[3]+fEdge[3]+fSkin[2]+fEdge[2])+0.5*(fSkin[1]+fSkin[0])-0.5*(fEdge[1]+fEdge[0]))+0.8660254037844386*(fSkin[3]+fSkin[2])-0.8660254037844386*(fEdge[3]+fEdge[2])+0.5*(fSkin[1]+fEdge[1]+fSkin[0]+fEdge[0]); 

  fUpwind[0] = 0.3535533905932737*(fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.3535533905932737*(fUpwindQuad[1]-1.0*fUpwindQuad[0]); 

  Ghat[0] = 0.7071067811865475*(alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] = 0.7071067811865475*(alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 

  out[0] += -0.7071067811865475*Ghat[0]*dv10; 
  out[1] += -0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -1.224744871391589*Ghat[0]*dv10; 
  out[3] += -1.224744871391589*Ghat[1]*dv10; 

  } else { 

  fUpwindQuad[0] = copysignf(1.0,alpha[0]-alpha[1])*((-0.8660254037844386*(fSkin[3]+fEdge[3]))+0.8660254037844386*(fSkin[2]+fEdge[2])+0.5*fSkin[1]-0.5*(fEdge[1]+fSkin[0])+0.5*fEdge[0])+0.8660254037844386*fSkin[3]-0.8660254037844386*(fEdge[3]+fSkin[2])+0.8660254037844386*fEdge[2]-0.5*(fSkin[1]+fEdge[1])+0.5*(fSkin[0]+fEdge[0]); 
  fUpwindQuad[1] = copysignf(1.0,alpha[1]+alpha[0])*(0.8660254037844386*(fSkin[3]+fEdge[3]+fSkin[2]+fEdge[2])-0.5*(fSkin[1]+fSkin[0])+0.5*(fEdge[1]+fEdge[0]))-0.8660254037844386*(fSkin[3]+fSkin[2])+0.8660254037844386*(fEdge[3]+fEdge[2])+0.5*(fSkin[1]+fEdge[1]+fSkin[0]+fEdge[0]); 

  fUpwind[0] = 0.3535533905932737*(fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.3535533905932737*(fUpwindQuad[1]-1.0*fUpwindQuad[0]); 

  Ghat[0] = 0.7071067811865475*(alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] = 0.7071067811865475*(alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 

  out[0] += 0.7071067811865475*Ghat[0]*dv10; 
  out[1] += 0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -1.224744871391589*Ghat[0]*dv10; 
  out[3] += -1.224744871391589*Ghat[1]*dv10; 

  } 
} 
