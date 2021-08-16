#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_boundary_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *qmem, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // qmem:        q/m*EM fields.
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv10 = 2/dxv[1]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double *E0 = &qmem[0]; 
  const double *B2 = &qmem[15]; 

  double Ghat[8]; 
  double alpha[8]; 

  alpha[0] = 1.414213562373095*(B2[0]*wv2+E0[0]); 
  alpha[1] = 1.414213562373095*(B2[1]*wv2+E0[1]); 
  alpha[2] = 0.408248290463863*B2[0]*dv2; 
  alpha[3] = 0.408248290463863*B2[1]*dv2; 
  alpha[4] = 1.414213562373095*(B2[2]*wv2+E0[2]); 
  alpha[6] = 0.408248290463863*B2[2]*dv2; 

  double fUpwindQuad[8];
  double fUpwind[8];

  if (edge == -1) { 

  if ((-0.5999999999999999*alpha[6])+0.4472135954999579*alpha[4]+0.9*alpha[3]-0.6708203932499369*(alpha[2]+alpha[1])+0.5*alpha[0] > 0) { 

    fUpwindQuad[0] = (-0.7348469228349525*fSkin[19])+1.42302494707577*fSkin[18]-0.7348469228349533*fSkin[17]+0.5477225575051661*fSkin[16]-0.4242640687119281*fSkin[15]-1.060660171779821*fSkin[14]-0.4242640687119285*fSkin[13]-1.060660171779821*fSkin[12]+0.5477225575051661*fSkin[11]+1.10227038425243*fSkin[10]+0.3162277660168379*fSkin[9]+0.7905694150420947*fSkin[8]+0.3162277660168379*fSkin[7]-0.8215838362577489*fSkin[6]+0.6363961030678926*fSkin[5]-0.8215838362577489*fSkin[4]-0.4743416490252568*fSkin[3]+0.6123724356957944*fSkin[2]-0.4743416490252568*fSkin[1]+0.3535533905932737*fSkin[0]; 
  } else { 

    fUpwindQuad[0] = 0.7348469228349525*fEdge[19]+1.42302494707577*fEdge[18]+0.7348469228349533*fEdge[17]-0.5477225575051661*fEdge[16]-0.4242640687119281*fEdge[15]-1.060660171779821*fEdge[14]-0.4242640687119285*fEdge[13]-1.060660171779821*fEdge[12]-0.5477225575051661*fEdge[11]-1.10227038425243*fEdge[10]+0.3162277660168379*fEdge[9]+0.7905694150420947*fEdge[8]+0.3162277660168379*fEdge[7]+0.8215838362577489*fEdge[6]+0.6363961030678926*fEdge[5]+0.8215838362577489*fEdge[4]-0.4743416490252568*fEdge[3]-0.6123724356957944*fEdge[2]-0.4743416490252568*fEdge[1]+0.3535533905932737*fEdge[0]; 
  } 
  if (0.75*alpha[6]-0.5590169943749475*alpha[4]-0.6708203932499369*alpha[2]+0.5*alpha[0] > 0) { 

    fUpwindQuad[1] = 0.9185586535436913*fSkin[17]+0.5477225575051661*fSkin[16]-1.060660171779821*fSkin[14]+0.5303300858899104*fSkin[13]-0.6846531968814574*fSkin[11]+0.3162277660168379*fSkin[9]+0.7905694150420947*fSkin[8]-0.3952847075210473*fSkin[7]-0.8215838362577489*fSkin[6]-0.4743416490252568*fSkin[3]+0.6123724356957944*fSkin[2]+0.3535533905932737*fSkin[0]; 
  } else { 

    fUpwindQuad[1] = (-0.9185586535436913*fEdge[17])-0.5477225575051661*fEdge[16]-1.060660171779821*fEdge[14]+0.5303300858899104*fEdge[13]+0.6846531968814574*fEdge[11]+0.3162277660168379*fEdge[9]+0.7905694150420947*fEdge[8]-0.3952847075210473*fEdge[7]+0.8215838362577489*fEdge[6]-0.4743416490252568*fEdge[3]-0.6123724356957944*fEdge[2]+0.3535533905932737*fEdge[0]; 
  } 
  if ((-0.5999999999999999*alpha[6])+0.4472135954999579*alpha[4]-0.9*alpha[3]-0.6708203932499369*alpha[2]+0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 

    fUpwindQuad[2] = 0.7348469228349525*fSkin[19]-1.42302494707577*fSkin[18]-0.7348469228349533*fSkin[17]+0.5477225575051661*fSkin[16]+0.4242640687119281*fSkin[15]-1.060660171779821*fSkin[14]-0.4242640687119285*fSkin[13]+1.060660171779821*fSkin[12]+0.5477225575051661*fSkin[11]-1.10227038425243*fSkin[10]+0.3162277660168379*fSkin[9]+0.7905694150420947*fSkin[8]+0.3162277660168379*fSkin[7]-0.8215838362577489*fSkin[6]-0.6363961030678926*fSkin[5]+0.8215838362577489*fSkin[4]-0.4743416490252568*fSkin[3]+0.6123724356957944*fSkin[2]+0.4743416490252568*fSkin[1]+0.3535533905932737*fSkin[0]; 
  } else { 

    fUpwindQuad[2] = (-0.7348469228349525*fEdge[19])-1.42302494707577*fEdge[18]+0.7348469228349533*fEdge[17]-0.5477225575051661*fEdge[16]+0.4242640687119281*fEdge[15]-1.060660171779821*fEdge[14]-0.4242640687119285*fEdge[13]+1.060660171779821*fEdge[12]-0.5477225575051661*fEdge[11]+1.10227038425243*fEdge[10]+0.3162277660168379*fEdge[9]+0.7905694150420947*fEdge[8]+0.3162277660168379*fEdge[7]+0.8215838362577489*fEdge[6]-0.6363961030678926*fEdge[5]-0.8215838362577489*fEdge[4]-0.4743416490252568*fEdge[3]-0.6123724356957944*fEdge[2]+0.4743416490252568*fEdge[1]+0.3535533905932737*fEdge[0]; 
  } 
  if (0.4472135954999579*alpha[4]-0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 

    fUpwindQuad[3] = 0.9185586535436913*fSkin[19]-0.6846531968814574*fSkin[16]+0.5303300858899104*fSkin[15]-1.060660171779821*fSkin[12]+0.5477225575051661*fSkin[11]-0.3952847075210473*fSkin[9]+0.7905694150420947*fSkin[8]+0.3162277660168379*fSkin[7]-0.8215838362577489*fSkin[4]+0.6123724356957944*fSkin[2]-0.4743416490252568*fSkin[1]+0.3535533905932737*fSkin[0]; 
  } else { 

    fUpwindQuad[3] = (-0.9185586535436913*fEdge[19])+0.6846531968814574*fEdge[16]+0.5303300858899104*fEdge[15]-1.060660171779821*fEdge[12]-0.5477225575051661*fEdge[11]-0.3952847075210473*fEdge[9]+0.7905694150420947*fEdge[8]+0.3162277660168379*fEdge[7]+0.8215838362577489*fEdge[4]-0.6123724356957944*fEdge[2]-0.4743416490252568*fEdge[1]+0.3535533905932737*fEdge[0]; 
  } 
  if (0.4472135954999579*alpha[4]+0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 

    fUpwindQuad[4] = (-0.9185586535436913*fSkin[19])-0.6846531968814574*fSkin[16]-0.5303300858899104*fSkin[15]+1.060660171779821*fSkin[12]+0.5477225575051661*fSkin[11]-0.3952847075210473*fSkin[9]+0.7905694150420947*fSkin[8]+0.3162277660168379*fSkin[7]+0.8215838362577489*fSkin[4]+0.6123724356957944*fSkin[2]+0.4743416490252568*fSkin[1]+0.3535533905932737*fSkin[0]; 
  } else { 

    fUpwindQuad[4] = 0.9185586535436913*fEdge[19]+0.6846531968814574*fEdge[16]-0.5303300858899104*fEdge[15]+1.060660171779821*fEdge[12]-0.5477225575051661*fEdge[11]-0.3952847075210473*fEdge[9]+0.7905694150420947*fEdge[8]+0.3162277660168379*fEdge[7]-0.8215838362577489*fEdge[4]-0.6123724356957944*fEdge[2]+0.4743416490252568*fEdge[1]+0.3535533905932737*fEdge[0]; 
  } 
  if (0.5999999999999999*alpha[6]+0.4472135954999579*alpha[4]-0.9*alpha[3]+0.6708203932499369*alpha[2]-0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 

    fUpwindQuad[5] = (-0.7348469228349525*fSkin[19])-1.42302494707577*fSkin[18]+0.7348469228349533*fSkin[17]+0.5477225575051661*fSkin[16]-0.4242640687119281*fSkin[15]+1.060660171779821*fSkin[14]+0.4242640687119285*fSkin[13]-1.060660171779821*fSkin[12]+0.5477225575051661*fSkin[11]-1.10227038425243*fSkin[10]+0.3162277660168379*fSkin[9]+0.7905694150420947*fSkin[8]+0.3162277660168379*fSkin[7]+0.8215838362577489*fSkin[6]-0.6363961030678926*fSkin[5]-0.8215838362577489*fSkin[4]+0.4743416490252568*fSkin[3]+0.6123724356957944*fSkin[2]-0.4743416490252568*fSkin[1]+0.3535533905932737*fSkin[0]; 
  } else { 

    fUpwindQuad[5] = 0.7348469228349525*fEdge[19]-1.42302494707577*fEdge[18]-0.7348469228349533*fEdge[17]-0.5477225575051661*fEdge[16]-0.4242640687119281*fEdge[15]+1.060660171779821*fEdge[14]+0.4242640687119285*fEdge[13]-1.060660171779821*fEdge[12]-0.5477225575051661*fEdge[11]+1.10227038425243*fEdge[10]+0.3162277660168379*fEdge[9]+0.7905694150420947*fEdge[8]+0.3162277660168379*fEdge[7]-0.8215838362577489*fEdge[6]-0.6363961030678926*fEdge[5]+0.8215838362577489*fEdge[4]+0.4743416490252568*fEdge[3]-0.6123724356957944*fEdge[2]-0.4743416490252568*fEdge[1]+0.3535533905932737*fEdge[0]; 
  } 
  if ((-0.75*alpha[6])-0.5590169943749475*alpha[4]+0.6708203932499369*alpha[2]+0.5*alpha[0] > 0) { 

    fUpwindQuad[6] = (-0.9185586535436913*fSkin[17])+0.5477225575051661*fSkin[16]+1.060660171779821*fSkin[14]-0.5303300858899104*fSkin[13]-0.6846531968814574*fSkin[11]+0.3162277660168379*fSkin[9]+0.7905694150420947*fSkin[8]-0.3952847075210473*fSkin[7]+0.8215838362577489*fSkin[6]+0.4743416490252568*fSkin[3]+0.6123724356957944*fSkin[2]+0.3535533905932737*fSkin[0]; 
  } else { 

    fUpwindQuad[6] = 0.9185586535436913*fEdge[17]-0.5477225575051661*fEdge[16]+1.060660171779821*fEdge[14]-0.5303300858899104*fEdge[13]+0.6846531968814574*fEdge[11]+0.3162277660168379*fEdge[9]+0.7905694150420947*fEdge[8]-0.3952847075210473*fEdge[7]-0.8215838362577489*fEdge[6]+0.4743416490252568*fEdge[3]-0.6123724356957944*fEdge[2]+0.3535533905932737*fEdge[0]; 
  } 
  if (0.5999999999999999*alpha[6]+0.4472135954999579*alpha[4]+0.9*alpha[3]+0.6708203932499369*(alpha[2]+alpha[1])+0.5*alpha[0] > 0) { 

    fUpwindQuad[7] = 0.7348469228349525*fSkin[19]+1.42302494707577*fSkin[18]+0.7348469228349533*fSkin[17]+0.5477225575051661*fSkin[16]+0.4242640687119281*fSkin[15]+1.060660171779821*fSkin[14]+0.4242640687119285*fSkin[13]+1.060660171779821*fSkin[12]+0.5477225575051661*fSkin[11]+1.10227038425243*fSkin[10]+0.3162277660168379*fSkin[9]+0.7905694150420947*fSkin[8]+0.3162277660168379*fSkin[7]+0.8215838362577489*fSkin[6]+0.6363961030678926*fSkin[5]+0.8215838362577489*fSkin[4]+0.4743416490252568*fSkin[3]+0.6123724356957944*fSkin[2]+0.4743416490252568*fSkin[1]+0.3535533905932737*fSkin[0]; 
  } else { 

    fUpwindQuad[7] = (-0.7348469228349525*fEdge[19])+1.42302494707577*fEdge[18]-0.7348469228349533*fEdge[17]-0.5477225575051661*fEdge[16]+0.4242640687119281*fEdge[15]+1.060660171779821*fEdge[14]+0.4242640687119285*fEdge[13]+1.060660171779821*fEdge[12]-0.5477225575051661*fEdge[11]-1.10227038425243*fEdge[10]+0.3162277660168379*fEdge[9]+0.7905694150420947*fEdge[8]+0.3162277660168379*fEdge[7]-0.8215838362577489*fEdge[6]+0.6363961030678926*fEdge[5]-0.8215838362577489*fEdge[4]+0.4743416490252568*fEdge[3]-0.6123724356957944*fEdge[2]+0.4743416490252568*fEdge[1]+0.3535533905932737*fEdge[0]; 
  } 

  fUpwind[0] = 0.05555555555555555*(fUpwindQuad[7]+8.0*fUpwindQuad[6]+fUpwindQuad[5]+8.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]+8.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.0414086662499961*(5.0*fUpwindQuad[7]+8.0*(fUpwindQuad[4]-1.0*fUpwindQuad[3])+5.0*fUpwindQuad[2]); 
  fUpwind[2] = 0.0414086662499961*(5.0*fUpwindQuad[7]+8.0*fUpwindQuad[6]+5.0*fUpwindQuad[5]-1.0*(5.0*fUpwindQuad[2]+8.0*fUpwindQuad[1]+5.0*fUpwindQuad[0])); 
  fUpwind[3] = 0.2777777777777778*(fUpwindQuad[7]-1.0*(fUpwindQuad[5]+fUpwindQuad[2])+fUpwindQuad[0]); 
  fUpwind[4] = 0.2484519974999766*(fUpwindQuad[7]-2.0*fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[2]-2.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[5] = 0.2484519974999766*(fUpwindQuad[7]+fUpwindQuad[5]-2.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]+fUpwindQuad[0]); 
  fUpwind[6] = 0.1851851851851853*(fUpwindQuad[7]-2.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[2]+2.0*fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[7] = 0.1851851851851853*(fUpwindQuad[7]-1.0*fUpwindQuad[5]+2.0*(fUpwindQuad[3]-1.0*fUpwindQuad[4])+fUpwindQuad[2]-1.0*fUpwindQuad[0]); 

  Ghat[0] = 0.5*(alpha[6]*fUpwind[6]+alpha[4]*fUpwind[4]+alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] = 0.1*(4.47213595499958*(alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6]+alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4])+5.0*(alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1])); 
  Ghat[2] = 0.1*(4.47213595499958*alpha[3]*fUpwind[7]+5.0*(alpha[4]*fUpwind[6]+fUpwind[4]*alpha[6])+4.47213595499958*alpha[2]*fUpwind[5]+5.0*(alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2])); 
  Ghat[3] = 0.1*((4.0*alpha[6]+4.47213595499958*alpha[2])*fUpwind[7]+4.47213595499958*(alpha[1]*fUpwind[6]+fUpwind[1]*alpha[6]+alpha[3]*(fUpwind[5]+fUpwind[4])+fUpwind[3]*alpha[4])+5.0*(alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2])); 
  Ghat[4] = 0.01428571428571429*((22.3606797749979*alpha[6]+35.0*alpha[2])*fUpwind[6]+35.0*fUpwind[2]*alpha[6]+(22.3606797749979*alpha[4]+35.0*alpha[0])*fUpwind[4]+35.0*fUpwind[0]*alpha[4]+31.30495168499706*(alpha[3]*fUpwind[3]+alpha[1]*fUpwind[1])); 
  Ghat[5] = 0.1*(5.0*alpha[1]*fUpwind[7]+4.47213595499958*alpha[6]*fUpwind[6]+5.0*alpha[0]*fUpwind[5]+4.47213595499958*(alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2])); 
  Ghat[6] = 0.01428571428571429*(28.0*alpha[3]*fUpwind[7]+(22.3606797749979*alpha[4]+35.0*alpha[0])*fUpwind[6]+(31.30495168499706*fUpwind[5]+22.3606797749979*fUpwind[4]+35.0*fUpwind[0])*alpha[6]+35.0*(alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4])+31.30495168499706*(alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3])); 
  Ghat[7] = 0.1*((4.47213595499958*alpha[4]+5.0*alpha[0])*fUpwind[7]+4.0*(alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6])+5.0*alpha[1]*fUpwind[5]+4.47213595499958*(alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3])); 

  out[0] += -0.7071067811865475*Ghat[0]*dv10; 
  out[1] += -0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -1.224744871391589*Ghat[0]*dv10; 
  out[3] += -0.7071067811865475*Ghat[2]*dv10; 
  out[4] += -1.224744871391589*Ghat[1]*dv10; 
  out[5] += -0.7071067811865475*Ghat[3]*dv10; 
  out[6] += -1.224744871391589*Ghat[2]*dv10; 
  out[7] += -0.7071067811865475*Ghat[4]*dv10; 
  out[8] += -1.58113883008419*Ghat[0]*dv10; 
  out[9] += -0.7071067811865475*Ghat[5]*dv10; 
  out[10] += -1.224744871391589*Ghat[3]*dv10; 
  out[11] += -1.224744871391589*Ghat[4]*dv10; 
  out[12] += -1.58113883008419*Ghat[1]*dv10; 
  out[13] += -0.7071067811865475*Ghat[6]*dv10; 
  out[14] += -1.58113883008419*Ghat[2]*dv10; 
  out[15] += -0.7071067811865475*Ghat[7]*dv10; 
  out[16] += -1.224744871391589*Ghat[5]*dv10; 
  out[17] += -1.224744871391589*Ghat[6]*dv10; 
  out[18] += -1.58113883008419*Ghat[3]*dv10; 
  out[19] += -1.224744871391589*Ghat[7]*dv10; 

  } else { 

  if ((-0.5999999999999999*alpha[6])+0.4472135954999579*alpha[4]+0.9*alpha[3]-0.6708203932499369*(alpha[2]+alpha[1])+0.5*alpha[0] > 0) { 

    fUpwindQuad[0] = (-0.7348469228349525*fEdge[19])+1.42302494707577*fEdge[18]-0.7348469228349533*fEdge[17]+0.5477225575051661*fEdge[16]-0.4242640687119281*fEdge[15]-1.060660171779821*fEdge[14]-0.4242640687119285*fEdge[13]-1.060660171779821*fEdge[12]+0.5477225575051661*fEdge[11]+1.10227038425243*fEdge[10]+0.3162277660168379*fEdge[9]+0.7905694150420947*fEdge[8]+0.3162277660168379*fEdge[7]-0.8215838362577489*fEdge[6]+0.6363961030678926*fEdge[5]-0.8215838362577489*fEdge[4]-0.4743416490252568*fEdge[3]+0.6123724356957944*fEdge[2]-0.4743416490252568*fEdge[1]+0.3535533905932737*fEdge[0]; 
  } else { 

    fUpwindQuad[0] = 0.7348469228349525*fSkin[19]+1.42302494707577*fSkin[18]+0.7348469228349533*fSkin[17]-0.5477225575051661*fSkin[16]-0.4242640687119281*fSkin[15]-1.060660171779821*fSkin[14]-0.4242640687119285*fSkin[13]-1.060660171779821*fSkin[12]-0.5477225575051661*fSkin[11]-1.10227038425243*fSkin[10]+0.3162277660168379*fSkin[9]+0.7905694150420947*fSkin[8]+0.3162277660168379*fSkin[7]+0.8215838362577489*fSkin[6]+0.6363961030678926*fSkin[5]+0.8215838362577489*fSkin[4]-0.4743416490252568*fSkin[3]-0.6123724356957944*fSkin[2]-0.4743416490252568*fSkin[1]+0.3535533905932737*fSkin[0]; 
  } 
  if (0.75*alpha[6]-0.5590169943749475*alpha[4]-0.6708203932499369*alpha[2]+0.5*alpha[0] > 0) { 

    fUpwindQuad[1] = 0.9185586535436913*fEdge[17]+0.5477225575051661*fEdge[16]-1.060660171779821*fEdge[14]+0.5303300858899104*fEdge[13]-0.6846531968814574*fEdge[11]+0.3162277660168379*fEdge[9]+0.7905694150420947*fEdge[8]-0.3952847075210473*fEdge[7]-0.8215838362577489*fEdge[6]-0.4743416490252568*fEdge[3]+0.6123724356957944*fEdge[2]+0.3535533905932737*fEdge[0]; 
  } else { 

    fUpwindQuad[1] = (-0.9185586535436913*fSkin[17])-0.5477225575051661*fSkin[16]-1.060660171779821*fSkin[14]+0.5303300858899104*fSkin[13]+0.6846531968814574*fSkin[11]+0.3162277660168379*fSkin[9]+0.7905694150420947*fSkin[8]-0.3952847075210473*fSkin[7]+0.8215838362577489*fSkin[6]-0.4743416490252568*fSkin[3]-0.6123724356957944*fSkin[2]+0.3535533905932737*fSkin[0]; 
  } 
  if ((-0.5999999999999999*alpha[6])+0.4472135954999579*alpha[4]-0.9*alpha[3]-0.6708203932499369*alpha[2]+0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 

    fUpwindQuad[2] = 0.7348469228349525*fEdge[19]-1.42302494707577*fEdge[18]-0.7348469228349533*fEdge[17]+0.5477225575051661*fEdge[16]+0.4242640687119281*fEdge[15]-1.060660171779821*fEdge[14]-0.4242640687119285*fEdge[13]+1.060660171779821*fEdge[12]+0.5477225575051661*fEdge[11]-1.10227038425243*fEdge[10]+0.3162277660168379*fEdge[9]+0.7905694150420947*fEdge[8]+0.3162277660168379*fEdge[7]-0.8215838362577489*fEdge[6]-0.6363961030678926*fEdge[5]+0.8215838362577489*fEdge[4]-0.4743416490252568*fEdge[3]+0.6123724356957944*fEdge[2]+0.4743416490252568*fEdge[1]+0.3535533905932737*fEdge[0]; 
  } else { 

    fUpwindQuad[2] = (-0.7348469228349525*fSkin[19])-1.42302494707577*fSkin[18]+0.7348469228349533*fSkin[17]-0.5477225575051661*fSkin[16]+0.4242640687119281*fSkin[15]-1.060660171779821*fSkin[14]-0.4242640687119285*fSkin[13]+1.060660171779821*fSkin[12]-0.5477225575051661*fSkin[11]+1.10227038425243*fSkin[10]+0.3162277660168379*fSkin[9]+0.7905694150420947*fSkin[8]+0.3162277660168379*fSkin[7]+0.8215838362577489*fSkin[6]-0.6363961030678926*fSkin[5]-0.8215838362577489*fSkin[4]-0.4743416490252568*fSkin[3]-0.6123724356957944*fSkin[2]+0.4743416490252568*fSkin[1]+0.3535533905932737*fSkin[0]; 
  } 
  if (0.4472135954999579*alpha[4]-0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 

    fUpwindQuad[3] = 0.9185586535436913*fEdge[19]-0.6846531968814574*fEdge[16]+0.5303300858899104*fEdge[15]-1.060660171779821*fEdge[12]+0.5477225575051661*fEdge[11]-0.3952847075210473*fEdge[9]+0.7905694150420947*fEdge[8]+0.3162277660168379*fEdge[7]-0.8215838362577489*fEdge[4]+0.6123724356957944*fEdge[2]-0.4743416490252568*fEdge[1]+0.3535533905932737*fEdge[0]; 
  } else { 

    fUpwindQuad[3] = (-0.9185586535436913*fSkin[19])+0.6846531968814574*fSkin[16]+0.5303300858899104*fSkin[15]-1.060660171779821*fSkin[12]-0.5477225575051661*fSkin[11]-0.3952847075210473*fSkin[9]+0.7905694150420947*fSkin[8]+0.3162277660168379*fSkin[7]+0.8215838362577489*fSkin[4]-0.6123724356957944*fSkin[2]-0.4743416490252568*fSkin[1]+0.3535533905932737*fSkin[0]; 
  } 
  if (0.4472135954999579*alpha[4]+0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 

    fUpwindQuad[4] = (-0.9185586535436913*fEdge[19])-0.6846531968814574*fEdge[16]-0.5303300858899104*fEdge[15]+1.060660171779821*fEdge[12]+0.5477225575051661*fEdge[11]-0.3952847075210473*fEdge[9]+0.7905694150420947*fEdge[8]+0.3162277660168379*fEdge[7]+0.8215838362577489*fEdge[4]+0.6123724356957944*fEdge[2]+0.4743416490252568*fEdge[1]+0.3535533905932737*fEdge[0]; 
  } else { 

    fUpwindQuad[4] = 0.9185586535436913*fSkin[19]+0.6846531968814574*fSkin[16]-0.5303300858899104*fSkin[15]+1.060660171779821*fSkin[12]-0.5477225575051661*fSkin[11]-0.3952847075210473*fSkin[9]+0.7905694150420947*fSkin[8]+0.3162277660168379*fSkin[7]-0.8215838362577489*fSkin[4]-0.6123724356957944*fSkin[2]+0.4743416490252568*fSkin[1]+0.3535533905932737*fSkin[0]; 
  } 
  if (0.5999999999999999*alpha[6]+0.4472135954999579*alpha[4]-0.9*alpha[3]+0.6708203932499369*alpha[2]-0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 

    fUpwindQuad[5] = (-0.7348469228349525*fEdge[19])-1.42302494707577*fEdge[18]+0.7348469228349533*fEdge[17]+0.5477225575051661*fEdge[16]-0.4242640687119281*fEdge[15]+1.060660171779821*fEdge[14]+0.4242640687119285*fEdge[13]-1.060660171779821*fEdge[12]+0.5477225575051661*fEdge[11]-1.10227038425243*fEdge[10]+0.3162277660168379*fEdge[9]+0.7905694150420947*fEdge[8]+0.3162277660168379*fEdge[7]+0.8215838362577489*fEdge[6]-0.6363961030678926*fEdge[5]-0.8215838362577489*fEdge[4]+0.4743416490252568*fEdge[3]+0.6123724356957944*fEdge[2]-0.4743416490252568*fEdge[1]+0.3535533905932737*fEdge[0]; 
  } else { 

    fUpwindQuad[5] = 0.7348469228349525*fSkin[19]-1.42302494707577*fSkin[18]-0.7348469228349533*fSkin[17]-0.5477225575051661*fSkin[16]-0.4242640687119281*fSkin[15]+1.060660171779821*fSkin[14]+0.4242640687119285*fSkin[13]-1.060660171779821*fSkin[12]-0.5477225575051661*fSkin[11]+1.10227038425243*fSkin[10]+0.3162277660168379*fSkin[9]+0.7905694150420947*fSkin[8]+0.3162277660168379*fSkin[7]-0.8215838362577489*fSkin[6]-0.6363961030678926*fSkin[5]+0.8215838362577489*fSkin[4]+0.4743416490252568*fSkin[3]-0.6123724356957944*fSkin[2]-0.4743416490252568*fSkin[1]+0.3535533905932737*fSkin[0]; 
  } 
  if ((-0.75*alpha[6])-0.5590169943749475*alpha[4]+0.6708203932499369*alpha[2]+0.5*alpha[0] > 0) { 

    fUpwindQuad[6] = (-0.9185586535436913*fEdge[17])+0.5477225575051661*fEdge[16]+1.060660171779821*fEdge[14]-0.5303300858899104*fEdge[13]-0.6846531968814574*fEdge[11]+0.3162277660168379*fEdge[9]+0.7905694150420947*fEdge[8]-0.3952847075210473*fEdge[7]+0.8215838362577489*fEdge[6]+0.4743416490252568*fEdge[3]+0.6123724356957944*fEdge[2]+0.3535533905932737*fEdge[0]; 
  } else { 

    fUpwindQuad[6] = 0.9185586535436913*fSkin[17]-0.5477225575051661*fSkin[16]+1.060660171779821*fSkin[14]-0.5303300858899104*fSkin[13]+0.6846531968814574*fSkin[11]+0.3162277660168379*fSkin[9]+0.7905694150420947*fSkin[8]-0.3952847075210473*fSkin[7]-0.8215838362577489*fSkin[6]+0.4743416490252568*fSkin[3]-0.6123724356957944*fSkin[2]+0.3535533905932737*fSkin[0]; 
  } 
  if (0.5999999999999999*alpha[6]+0.4472135954999579*alpha[4]+0.9*alpha[3]+0.6708203932499369*(alpha[2]+alpha[1])+0.5*alpha[0] > 0) { 

    fUpwindQuad[7] = 0.7348469228349525*fEdge[19]+1.42302494707577*fEdge[18]+0.7348469228349533*fEdge[17]+0.5477225575051661*fEdge[16]+0.4242640687119281*fEdge[15]+1.060660171779821*fEdge[14]+0.4242640687119285*fEdge[13]+1.060660171779821*fEdge[12]+0.5477225575051661*fEdge[11]+1.10227038425243*fEdge[10]+0.3162277660168379*fEdge[9]+0.7905694150420947*fEdge[8]+0.3162277660168379*fEdge[7]+0.8215838362577489*fEdge[6]+0.6363961030678926*fEdge[5]+0.8215838362577489*fEdge[4]+0.4743416490252568*fEdge[3]+0.6123724356957944*fEdge[2]+0.4743416490252568*fEdge[1]+0.3535533905932737*fEdge[0]; 
  } else { 

    fUpwindQuad[7] = (-0.7348469228349525*fSkin[19])+1.42302494707577*fSkin[18]-0.7348469228349533*fSkin[17]-0.5477225575051661*fSkin[16]+0.4242640687119281*fSkin[15]+1.060660171779821*fSkin[14]+0.4242640687119285*fSkin[13]+1.060660171779821*fSkin[12]-0.5477225575051661*fSkin[11]-1.10227038425243*fSkin[10]+0.3162277660168379*fSkin[9]+0.7905694150420947*fSkin[8]+0.3162277660168379*fSkin[7]-0.8215838362577489*fSkin[6]+0.6363961030678926*fSkin[5]-0.8215838362577489*fSkin[4]+0.4743416490252568*fSkin[3]-0.6123724356957944*fSkin[2]+0.4743416490252568*fSkin[1]+0.3535533905932737*fSkin[0]; 
  } 

  fUpwind[0] = 0.05555555555555555*(fUpwindQuad[7]+8.0*fUpwindQuad[6]+fUpwindQuad[5]+8.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]+8.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.0414086662499961*(5.0*fUpwindQuad[7]+8.0*(fUpwindQuad[4]-1.0*fUpwindQuad[3])+5.0*fUpwindQuad[2]); 
  fUpwind[2] = 0.0414086662499961*(5.0*fUpwindQuad[7]+8.0*fUpwindQuad[6]+5.0*fUpwindQuad[5]-1.0*(5.0*fUpwindQuad[2]+8.0*fUpwindQuad[1]+5.0*fUpwindQuad[0])); 
  fUpwind[3] = 0.2777777777777778*(fUpwindQuad[7]-1.0*(fUpwindQuad[5]+fUpwindQuad[2])+fUpwindQuad[0]); 
  fUpwind[4] = 0.2484519974999766*(fUpwindQuad[7]-2.0*fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[2]-2.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[5] = 0.2484519974999766*(fUpwindQuad[7]+fUpwindQuad[5]-2.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]+fUpwindQuad[0]); 
  fUpwind[6] = 0.1851851851851853*(fUpwindQuad[7]-2.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[2]+2.0*fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[7] = 0.1851851851851853*(fUpwindQuad[7]-1.0*fUpwindQuad[5]+2.0*(fUpwindQuad[3]-1.0*fUpwindQuad[4])+fUpwindQuad[2]-1.0*fUpwindQuad[0]); 

  Ghat[0] = 0.5*(alpha[6]*fUpwind[6]+alpha[4]*fUpwind[4]+alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] = 0.1*(4.47213595499958*(alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6]+alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4])+5.0*(alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1])); 
  Ghat[2] = 0.1*(4.47213595499958*alpha[3]*fUpwind[7]+5.0*(alpha[4]*fUpwind[6]+fUpwind[4]*alpha[6])+4.47213595499958*alpha[2]*fUpwind[5]+5.0*(alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2])); 
  Ghat[3] = 0.1*((4.0*alpha[6]+4.47213595499958*alpha[2])*fUpwind[7]+4.47213595499958*(alpha[1]*fUpwind[6]+fUpwind[1]*alpha[6]+alpha[3]*(fUpwind[5]+fUpwind[4])+fUpwind[3]*alpha[4])+5.0*(alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2])); 
  Ghat[4] = 0.01428571428571429*((22.3606797749979*alpha[6]+35.0*alpha[2])*fUpwind[6]+35.0*fUpwind[2]*alpha[6]+(22.3606797749979*alpha[4]+35.0*alpha[0])*fUpwind[4]+35.0*fUpwind[0]*alpha[4]+31.30495168499706*(alpha[3]*fUpwind[3]+alpha[1]*fUpwind[1])); 
  Ghat[5] = 0.1*(5.0*alpha[1]*fUpwind[7]+4.47213595499958*alpha[6]*fUpwind[6]+5.0*alpha[0]*fUpwind[5]+4.47213595499958*(alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2])); 
  Ghat[6] = 0.01428571428571429*(28.0*alpha[3]*fUpwind[7]+(22.3606797749979*alpha[4]+35.0*alpha[0])*fUpwind[6]+(31.30495168499706*fUpwind[5]+22.3606797749979*fUpwind[4]+35.0*fUpwind[0])*alpha[6]+35.0*(alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4])+31.30495168499706*(alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3])); 
  Ghat[7] = 0.1*((4.47213595499958*alpha[4]+5.0*alpha[0])*fUpwind[7]+4.0*(alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6])+5.0*alpha[1]*fUpwind[5]+4.47213595499958*(alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3])); 

  out[0] += 0.7071067811865475*Ghat[0]*dv10; 
  out[1] += 0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -1.224744871391589*Ghat[0]*dv10; 
  out[3] += 0.7071067811865475*Ghat[2]*dv10; 
  out[4] += -1.224744871391589*Ghat[1]*dv10; 
  out[5] += 0.7071067811865475*Ghat[3]*dv10; 
  out[6] += -1.224744871391589*Ghat[2]*dv10; 
  out[7] += 0.7071067811865475*Ghat[4]*dv10; 
  out[8] += 1.58113883008419*Ghat[0]*dv10; 
  out[9] += 0.7071067811865475*Ghat[5]*dv10; 
  out[10] += -1.224744871391589*Ghat[3]*dv10; 
  out[11] += -1.224744871391589*Ghat[4]*dv10; 
  out[12] += 1.58113883008419*Ghat[1]*dv10; 
  out[13] += 0.7071067811865475*Ghat[6]*dv10; 
  out[14] += 1.58113883008419*Ghat[2]*dv10; 
  out[15] += 0.7071067811865475*Ghat[7]*dv10; 
  out[16] += -1.224744871391589*Ghat[5]*dv10; 
  out[17] += -1.224744871391589*Ghat[6]*dv10; 
  out[18] += 1.58113883008419*Ghat[3]*dv10; 
  out[19] += -1.224744871391589*Ghat[7]*dv10; 

  } 
} 
