#include <gkyl_vlasov_lbo_kernels.h> 
#include <gkyl_basis_ser_1x2v_p2_surfvy_quad.h> 
GKYL_CU_DH void vlasov_lbo_drag_boundary_surfvy_1x2v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[3]:         Cell-center coordinates. 
  // dxv[3]:       Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[6]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]: sum of thermal speeds squared time their respective collisionalities. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[2]; 

  const double *sumNuUy = &nuUSum[3]; 

  double alphaDrSurf[8] = {0.0}; 
  double fUpwindQuad[9] = {0.0};
  double fUpwind[8] = {0.0};;
  double drag_incr[8] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*(2.828427124746191*w[2]+1.414213562373095*dxv[2])-2.828427124746191*sumNuUy[0]); 
  alphaDrSurf[1] = 0.5*(nuSum[1]*(2.828427124746191*w[2]+1.414213562373095*dxv[2])-2.828427124746191*sumNuUy[1]); 
  alphaDrSurf[4] = 0.5*(2.828427124746191*nuSum[2]*w[2]-2.828427124746191*sumNuUy[2]+1.414213562373095*dxv[2]*nuSum[2]); 

  if (0.4472135954999579*alphaDrSurf[4]-0.6708203932499369*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = ser_1x2v_p2_surfvy_quad_0(1, fSkin); 
  } else { 
    fUpwindQuad[0] = ser_1x2v_p2_surfvy_quad_0(-1, fEdge); 
  } 
  if (0.5*alphaDrSurf[0]-0.5590169943749475*alphaDrSurf[4] < 0) { 
    fUpwindQuad[1] = ser_1x2v_p2_surfvy_quad_1(1, fSkin); 
  } else { 
    fUpwindQuad[1] = ser_1x2v_p2_surfvy_quad_1(-1, fEdge); 
  } 
  if (0.4472135954999579*alphaDrSurf[4]+0.6708203932499369*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[2] = ser_1x2v_p2_surfvy_quad_2(1, fSkin); 
  } else { 
    fUpwindQuad[2] = ser_1x2v_p2_surfvy_quad_2(-1, fEdge); 
  } 
  if (0.4472135954999579*alphaDrSurf[4]-0.6708203932499369*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = ser_1x2v_p2_surfvy_quad_3(1, fSkin); 
  } else { 
    fUpwindQuad[3] = ser_1x2v_p2_surfvy_quad_3(-1, fEdge); 
  } 
  if (0.5*alphaDrSurf[0]-0.5590169943749475*alphaDrSurf[4] < 0) { 
    fUpwindQuad[4] = ser_1x2v_p2_surfvy_quad_4(1, fSkin); 
  } else { 
    fUpwindQuad[4] = ser_1x2v_p2_surfvy_quad_4(-1, fEdge); 
  } 
  if (0.4472135954999579*alphaDrSurf[4]+0.6708203932499369*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[5] = ser_1x2v_p2_surfvy_quad_5(1, fSkin); 
  } else { 
    fUpwindQuad[5] = ser_1x2v_p2_surfvy_quad_5(-1, fEdge); 
  } 
  if (0.4472135954999579*alphaDrSurf[4]-0.6708203932499369*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[6] = ser_1x2v_p2_surfvy_quad_6(1, fSkin); 
  } else { 
    fUpwindQuad[6] = ser_1x2v_p2_surfvy_quad_6(-1, fEdge); 
  } 
  if (0.5*alphaDrSurf[0]-0.5590169943749475*alphaDrSurf[4] < 0) { 
    fUpwindQuad[7] = ser_1x2v_p2_surfvy_quad_7(1, fSkin); 
  } else { 
    fUpwindQuad[7] = ser_1x2v_p2_surfvy_quad_7(-1, fEdge); 
  } 
  if (0.4472135954999579*alphaDrSurf[4]+0.6708203932499369*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[8] = ser_1x2v_p2_surfvy_quad_8(1, fSkin); 
  } else { 
    fUpwindQuad[8] = ser_1x2v_p2_surfvy_quad_8(-1, fEdge); 
  } 

  fUpwind[0] = 0.154320987654321*fUpwindQuad[8]+0.2469135802469136*fUpwindQuad[7]+0.154320987654321*fUpwindQuad[6]+0.2469135802469136*fUpwindQuad[5]+0.3950617283950617*fUpwindQuad[4]+0.2469135802469136*fUpwindQuad[3]+0.154320987654321*fUpwindQuad[2]+0.2469135802469136*fUpwindQuad[1]+0.154320987654321*fUpwindQuad[0]; 
  fUpwind[1] = 0.2070433312499806*fUpwindQuad[8]-0.2070433312499806*fUpwindQuad[6]+0.3312693299999688*fUpwindQuad[5]-0.3312693299999688*fUpwindQuad[3]+0.2070433312499806*fUpwindQuad[2]-0.2070433312499806*fUpwindQuad[0]; 
  fUpwind[2] = 0.2070433312499806*fUpwindQuad[8]+0.3312693299999688*fUpwindQuad[7]+0.2070433312499806*fUpwindQuad[6]-0.2070433312499806*fUpwindQuad[2]-0.3312693299999688*fUpwindQuad[1]-0.2070433312499806*fUpwindQuad[0]; 
  fUpwind[3] = 0.2777777777777778*fUpwindQuad[8]-0.2777777777777778*fUpwindQuad[6]-0.2777777777777778*fUpwindQuad[2]+0.2777777777777778*fUpwindQuad[0]; 
  fUpwind[4] = 0.138028887499987*fUpwindQuad[8]-0.2760577749999741*fUpwindQuad[7]+0.138028887499987*fUpwindQuad[6]+0.2208462199999792*fUpwindQuad[5]-0.4416924399999584*fUpwindQuad[4]+0.2208462199999792*fUpwindQuad[3]+0.138028887499987*fUpwindQuad[2]-0.2760577749999741*fUpwindQuad[1]+0.138028887499987*fUpwindQuad[0]; 
  fUpwind[5] = 0.138028887499987*fUpwindQuad[8]+0.2208462199999792*fUpwindQuad[7]+0.138028887499987*fUpwindQuad[6]-0.2760577749999741*fUpwindQuad[5]-0.4416924399999584*fUpwindQuad[4]-0.2760577749999741*fUpwindQuad[3]+0.138028887499987*fUpwindQuad[2]+0.2208462199999792*fUpwindQuad[1]+0.138028887499987*fUpwindQuad[0]; 
  fUpwind[6] = 0.1851851851851853*fUpwindQuad[8]-0.3703703703703705*fUpwindQuad[7]+0.1851851851851853*fUpwindQuad[6]-0.1851851851851853*fUpwindQuad[2]+0.3703703703703705*fUpwindQuad[1]-0.1851851851851853*fUpwindQuad[0]; 
  fUpwind[7] = 0.1851851851851853*fUpwindQuad[8]-0.1851851851851853*fUpwindQuad[6]-0.3703703703703705*fUpwindQuad[5]+0.3703703703703705*fUpwindQuad[3]+0.1851851851851853*fUpwindQuad[2]-0.1851851851851853*fUpwindQuad[0]; 

  drag_incr[0] = 0.5*alphaDrSurf[4]*fUpwind[4]+0.5*alphaDrSurf[1]*fUpwind[1]+0.5*alphaDrSurf[0]*fUpwind[0]; 
  drag_incr[1] = 0.4472135954999579*alphaDrSurf[1]*fUpwind[4]+0.4472135954999579*fUpwind[1]*alphaDrSurf[4]+0.5*alphaDrSurf[0]*fUpwind[1]+0.5*fUpwind[0]*alphaDrSurf[1]; 
  drag_incr[2] = 0.5000000000000001*alphaDrSurf[4]*fUpwind[6]+0.5*alphaDrSurf[1]*fUpwind[3]+0.5*alphaDrSurf[0]*fUpwind[2]; 
  drag_incr[3] = 0.447213595499958*alphaDrSurf[1]*fUpwind[6]+0.4472135954999579*fUpwind[3]*alphaDrSurf[4]+0.5*alphaDrSurf[0]*fUpwind[3]+0.5*alphaDrSurf[1]*fUpwind[2]; 
  drag_incr[4] = 0.31943828249997*alphaDrSurf[4]*fUpwind[4]+0.5*alphaDrSurf[0]*fUpwind[4]+0.5*fUpwind[0]*alphaDrSurf[4]+0.4472135954999579*alphaDrSurf[1]*fUpwind[1]; 
  drag_incr[5] = 0.5000000000000001*alphaDrSurf[1]*fUpwind[7]+0.5*alphaDrSurf[0]*fUpwind[5]; 
  drag_incr[6] = 0.31943828249997*alphaDrSurf[4]*fUpwind[6]+0.5*alphaDrSurf[0]*fUpwind[6]+0.5000000000000001*fUpwind[2]*alphaDrSurf[4]+0.447213595499958*alphaDrSurf[1]*fUpwind[3]; 
  drag_incr[7] = 0.4472135954999579*alphaDrSurf[4]*fUpwind[7]+0.5*alphaDrSurf[0]*fUpwind[7]+0.5000000000000001*alphaDrSurf[1]*fUpwind[5]; 

  out[0] += 0.7071067811865475*drag_incr[0]*rdv2; 
  out[1] += 0.7071067811865475*drag_incr[1]*rdv2; 
  out[2] += 0.7071067811865475*drag_incr[2]*rdv2; 
  out[3] += 1.224744871391589*drag_incr[0]*rdv2; 
  out[4] += 0.7071067811865475*drag_incr[3]*rdv2; 
  out[5] += 1.224744871391589*drag_incr[1]*rdv2; 
  out[6] += 1.224744871391589*drag_incr[2]*rdv2; 
  out[7] += 0.7071067811865475*drag_incr[4]*rdv2; 
  out[8] += 0.7071067811865475*drag_incr[5]*rdv2; 
  out[9] += 1.58113883008419*drag_incr[0]*rdv2; 
  out[10] += 1.224744871391589*drag_incr[3]*rdv2; 
  out[11] += 0.7071067811865475*drag_incr[6]*rdv2; 
  out[12] += 0.7071067811865475*drag_incr[7]*rdv2; 
  out[13] += 1.224744871391589*drag_incr[4]*rdv2; 
  out[14] += 1.224744871391589*drag_incr[5]*rdv2; 
  out[15] += 1.58113883008419*drag_incr[1]*rdv2; 
  out[16] += 1.58113883008419*drag_incr[2]*rdv2; 
  out[17] += 1.224744871391589*drag_incr[6]*rdv2; 
  out[18] += 1.224744871391589*drag_incr[7]*rdv2; 
  out[19] += 1.58113883008419*drag_incr[3]*rdv2; 

  } else { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*(2.828427124746191*w[2]-1.414213562373095*dxv[2])-2.828427124746191*sumNuUy[0]); 
  alphaDrSurf[1] = 0.5*(nuSum[1]*(2.828427124746191*w[2]-1.414213562373095*dxv[2])-2.828427124746191*sumNuUy[1]); 
  alphaDrSurf[4] = 0.5*(2.828427124746191*nuSum[2]*w[2]-2.828427124746191*sumNuUy[2]-1.414213562373095*dxv[2]*nuSum[2]); 

  if (0.4472135954999579*alphaDrSurf[4]-0.6708203932499369*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = ser_1x2v_p2_surfvy_quad_0(1, fEdge); 
  } else { 
    fUpwindQuad[0] = ser_1x2v_p2_surfvy_quad_0(-1, fSkin); 
  } 
  if (0.5*alphaDrSurf[0]-0.5590169943749475*alphaDrSurf[4] < 0) { 
    fUpwindQuad[1] = ser_1x2v_p2_surfvy_quad_1(1, fEdge); 
  } else { 
    fUpwindQuad[1] = ser_1x2v_p2_surfvy_quad_1(-1, fSkin); 
  } 
  if (0.4472135954999579*alphaDrSurf[4]+0.6708203932499369*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[2] = ser_1x2v_p2_surfvy_quad_2(1, fEdge); 
  } else { 
    fUpwindQuad[2] = ser_1x2v_p2_surfvy_quad_2(-1, fSkin); 
  } 
  if (0.4472135954999579*alphaDrSurf[4]-0.6708203932499369*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = ser_1x2v_p2_surfvy_quad_3(1, fEdge); 
  } else { 
    fUpwindQuad[3] = ser_1x2v_p2_surfvy_quad_3(-1, fSkin); 
  } 
  if (0.5*alphaDrSurf[0]-0.5590169943749475*alphaDrSurf[4] < 0) { 
    fUpwindQuad[4] = ser_1x2v_p2_surfvy_quad_4(1, fEdge); 
  } else { 
    fUpwindQuad[4] = ser_1x2v_p2_surfvy_quad_4(-1, fSkin); 
  } 
  if (0.4472135954999579*alphaDrSurf[4]+0.6708203932499369*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[5] = ser_1x2v_p2_surfvy_quad_5(1, fEdge); 
  } else { 
    fUpwindQuad[5] = ser_1x2v_p2_surfvy_quad_5(-1, fSkin); 
  } 
  if (0.4472135954999579*alphaDrSurf[4]-0.6708203932499369*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[6] = ser_1x2v_p2_surfvy_quad_6(1, fEdge); 
  } else { 
    fUpwindQuad[6] = ser_1x2v_p2_surfvy_quad_6(-1, fSkin); 
  } 
  if (0.5*alphaDrSurf[0]-0.5590169943749475*alphaDrSurf[4] < 0) { 
    fUpwindQuad[7] = ser_1x2v_p2_surfvy_quad_7(1, fEdge); 
  } else { 
    fUpwindQuad[7] = ser_1x2v_p2_surfvy_quad_7(-1, fSkin); 
  } 
  if (0.4472135954999579*alphaDrSurf[4]+0.6708203932499369*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[8] = ser_1x2v_p2_surfvy_quad_8(1, fEdge); 
  } else { 
    fUpwindQuad[8] = ser_1x2v_p2_surfvy_quad_8(-1, fSkin); 
  } 

  fUpwind[0] = 0.154320987654321*fUpwindQuad[8]+0.2469135802469136*fUpwindQuad[7]+0.154320987654321*fUpwindQuad[6]+0.2469135802469136*fUpwindQuad[5]+0.3950617283950617*fUpwindQuad[4]+0.2469135802469136*fUpwindQuad[3]+0.154320987654321*fUpwindQuad[2]+0.2469135802469136*fUpwindQuad[1]+0.154320987654321*fUpwindQuad[0]; 
  fUpwind[1] = 0.2070433312499806*fUpwindQuad[8]-0.2070433312499806*fUpwindQuad[6]+0.3312693299999688*fUpwindQuad[5]-0.3312693299999688*fUpwindQuad[3]+0.2070433312499806*fUpwindQuad[2]-0.2070433312499806*fUpwindQuad[0]; 
  fUpwind[2] = 0.2070433312499806*fUpwindQuad[8]+0.3312693299999688*fUpwindQuad[7]+0.2070433312499806*fUpwindQuad[6]-0.2070433312499806*fUpwindQuad[2]-0.3312693299999688*fUpwindQuad[1]-0.2070433312499806*fUpwindQuad[0]; 
  fUpwind[3] = 0.2777777777777778*fUpwindQuad[8]-0.2777777777777778*fUpwindQuad[6]-0.2777777777777778*fUpwindQuad[2]+0.2777777777777778*fUpwindQuad[0]; 
  fUpwind[4] = 0.138028887499987*fUpwindQuad[8]-0.2760577749999741*fUpwindQuad[7]+0.138028887499987*fUpwindQuad[6]+0.2208462199999792*fUpwindQuad[5]-0.4416924399999584*fUpwindQuad[4]+0.2208462199999792*fUpwindQuad[3]+0.138028887499987*fUpwindQuad[2]-0.2760577749999741*fUpwindQuad[1]+0.138028887499987*fUpwindQuad[0]; 
  fUpwind[5] = 0.138028887499987*fUpwindQuad[8]+0.2208462199999792*fUpwindQuad[7]+0.138028887499987*fUpwindQuad[6]-0.2760577749999741*fUpwindQuad[5]-0.4416924399999584*fUpwindQuad[4]-0.2760577749999741*fUpwindQuad[3]+0.138028887499987*fUpwindQuad[2]+0.2208462199999792*fUpwindQuad[1]+0.138028887499987*fUpwindQuad[0]; 
  fUpwind[6] = 0.1851851851851853*fUpwindQuad[8]-0.3703703703703705*fUpwindQuad[7]+0.1851851851851853*fUpwindQuad[6]-0.1851851851851853*fUpwindQuad[2]+0.3703703703703705*fUpwindQuad[1]-0.1851851851851853*fUpwindQuad[0]; 
  fUpwind[7] = 0.1851851851851853*fUpwindQuad[8]-0.1851851851851853*fUpwindQuad[6]-0.3703703703703705*fUpwindQuad[5]+0.3703703703703705*fUpwindQuad[3]+0.1851851851851853*fUpwindQuad[2]-0.1851851851851853*fUpwindQuad[0]; 

  drag_incr[0] = 0.5*alphaDrSurf[4]*fUpwind[4]+0.5*alphaDrSurf[1]*fUpwind[1]+0.5*alphaDrSurf[0]*fUpwind[0]; 
  drag_incr[1] = 0.4472135954999579*alphaDrSurf[1]*fUpwind[4]+0.4472135954999579*fUpwind[1]*alphaDrSurf[4]+0.5*alphaDrSurf[0]*fUpwind[1]+0.5*fUpwind[0]*alphaDrSurf[1]; 
  drag_incr[2] = 0.5000000000000001*alphaDrSurf[4]*fUpwind[6]+0.5*alphaDrSurf[1]*fUpwind[3]+0.5*alphaDrSurf[0]*fUpwind[2]; 
  drag_incr[3] = 0.447213595499958*alphaDrSurf[1]*fUpwind[6]+0.4472135954999579*fUpwind[3]*alphaDrSurf[4]+0.5*alphaDrSurf[0]*fUpwind[3]+0.5*alphaDrSurf[1]*fUpwind[2]; 
  drag_incr[4] = 0.31943828249997*alphaDrSurf[4]*fUpwind[4]+0.5*alphaDrSurf[0]*fUpwind[4]+0.5*fUpwind[0]*alphaDrSurf[4]+0.4472135954999579*alphaDrSurf[1]*fUpwind[1]; 
  drag_incr[5] = 0.5000000000000001*alphaDrSurf[1]*fUpwind[7]+0.5*alphaDrSurf[0]*fUpwind[5]; 
  drag_incr[6] = 0.31943828249997*alphaDrSurf[4]*fUpwind[6]+0.5*alphaDrSurf[0]*fUpwind[6]+0.5000000000000001*fUpwind[2]*alphaDrSurf[4]+0.447213595499958*alphaDrSurf[1]*fUpwind[3]; 
  drag_incr[7] = 0.4472135954999579*alphaDrSurf[4]*fUpwind[7]+0.5*alphaDrSurf[0]*fUpwind[7]+0.5000000000000001*alphaDrSurf[1]*fUpwind[5]; 

  out[0] += -0.7071067811865475*drag_incr[0]*rdv2; 
  out[1] += -0.7071067811865475*drag_incr[1]*rdv2; 
  out[2] += -0.7071067811865475*drag_incr[2]*rdv2; 
  out[3] += 1.224744871391589*drag_incr[0]*rdv2; 
  out[4] += -0.7071067811865475*drag_incr[3]*rdv2; 
  out[5] += 1.224744871391589*drag_incr[1]*rdv2; 
  out[6] += 1.224744871391589*drag_incr[2]*rdv2; 
  out[7] += -0.7071067811865475*drag_incr[4]*rdv2; 
  out[8] += -0.7071067811865475*drag_incr[5]*rdv2; 
  out[9] += -1.58113883008419*drag_incr[0]*rdv2; 
  out[10] += 1.224744871391589*drag_incr[3]*rdv2; 
  out[11] += -0.7071067811865475*drag_incr[6]*rdv2; 
  out[12] += -0.7071067811865475*drag_incr[7]*rdv2; 
  out[13] += 1.224744871391589*drag_incr[4]*rdv2; 
  out[14] += 1.224744871391589*drag_incr[5]*rdv2; 
  out[15] += -1.58113883008419*drag_incr[1]*rdv2; 
  out[16] += -1.58113883008419*drag_incr[2]*rdv2; 
  out[17] += 1.224744871391589*drag_incr[6]*rdv2; 
  out[18] += 1.224744871391589*drag_incr[7]*rdv2; 
  out[19] += -1.58113883008419*drag_incr[3]*rdv2; 

  } 
} 