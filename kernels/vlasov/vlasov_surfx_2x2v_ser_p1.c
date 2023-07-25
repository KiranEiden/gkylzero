#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_surfx_2x2v_ser_p1(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // alpha_geo: Fields used only for general geometry.
  // fl/fc/fr:  Input Distribution function in left/center/right cells.
  // out:       Incremented distribution function in center cell.
  const double dx10 = 2/dxv[0]; 
  const double dv = dxv[2], wv = w[2]; 
  double Ghat_r[16]; 
  double Ghat_l[16]; 
  if (wv>0) { 

  Ghat_r[0] = (1.224744871391589*fc[1]+0.7071067811865475*fc[0])*wv+(0.3535533905932737*fc[6]+0.2041241452319315*fc[3])*dv; 
  Ghat_r[1] = (1.224744871391589*fc[5]+0.7071067811865475*fc[2])*wv+(0.3535533905932737*fc[11]+0.2041241452319315*fc[7])*dv; 
  Ghat_r[2] = (1.224744871391589*fc[6]+0.7071067811865475*fc[3])*wv+(0.3162277660168379*fc[17]+0.1825741858350554*fc[16]+0.3535533905932737*fc[1]+0.2041241452319315*fc[0])*dv; 
  Ghat_r[3] = (1.224744871391589*fc[8]+0.7071067811865475*fc[4])*wv+(0.3535533905932737*fc[13]+0.2041241452319315*fc[10])*dv; 
  Ghat_r[4] = (1.224744871391589*fc[11]+0.7071067811865475*fc[7])*wv+(0.3162277660168379*fc[20]+0.1825741858350554*fc[18]+0.3535533905932737*fc[5]+0.2041241452319315*fc[2])*dv; 
  Ghat_r[5] = (1.224744871391589*fc[12]+0.7071067811865475*fc[9])*wv+(0.3535533905932737*fc[15]+0.2041241452319315*fc[14])*dv; 
  Ghat_r[6] = (1.224744871391589*fc[13]+0.7071067811865475*fc[10])*wv+(0.3162277660168379*fc[21]+0.1825741858350554*fc[19]+0.3535533905932737*fc[8]+0.2041241452319315*fc[4])*dv; 
  Ghat_r[7] = (1.224744871391589*fc[15]+0.7071067811865475*fc[14])*wv+(0.3162277660168379*fc[23]+0.1825741858350554*fc[22]+0.3535533905932737*fc[12]+0.2041241452319315*fc[9])*dv; 
  Ghat_r[8] = (1.224744871391589*fc[17]+0.7071067811865475*fc[16])*wv+(0.3162277660168379*fc[6]+0.1825741858350554*fc[3])*dv; 
  Ghat_r[9] = (1.224744871391589*fc[20]+0.7071067811865475*fc[18])*wv+(0.3162277660168379*fc[11]+0.1825741858350554*fc[7])*dv; 
  Ghat_r[10] = (1.224744871391589*fc[21]+0.7071067811865475*fc[19])*wv+(0.3162277660168379*fc[13]+0.1825741858350554*fc[10])*dv; 
  Ghat_r[11] = (1.224744871391589*fc[23]+0.7071067811865475*fc[22])*wv+(0.3162277660168379*fc[15]+0.1825741858350554*fc[14])*dv; 
  Ghat_r[12] = (1.224744871391589*fc[25]+0.7071067811865475*fc[24])*wv+(0.3535533905932737*fc[29]+0.2041241452319315*fc[27])*dv; 
  Ghat_r[13] = (1.224744871391589*fc[28]+0.7071067811865475*fc[26])*wv+(0.3535533905932737*fc[31]+0.2041241452319315*fc[30])*dv; 
  Ghat_r[14] = (1.224744871391589*fc[29]+0.7071067811865475*fc[27])*wv+(0.3535533905932737*fc[25]+0.2041241452319315*fc[24])*dv; 
  Ghat_r[15] = (1.224744871391589*fc[31]+0.7071067811865475*fc[30])*wv+(0.3535533905932737*fc[28]+0.2041241452319315*fc[26])*dv; 

  Ghat_l[0] = (1.224744871391589*fl[1]+0.7071067811865475*fl[0])*wv+(0.3535533905932737*fl[6]+0.2041241452319315*fl[3])*dv; 
  Ghat_l[1] = (1.224744871391589*fl[5]+0.7071067811865475*fl[2])*wv+(0.3535533905932737*fl[11]+0.2041241452319315*fl[7])*dv; 
  Ghat_l[2] = (1.224744871391589*fl[6]+0.7071067811865475*fl[3])*wv+(0.3162277660168379*fl[17]+0.1825741858350554*fl[16]+0.3535533905932737*fl[1]+0.2041241452319315*fl[0])*dv; 
  Ghat_l[3] = (1.224744871391589*fl[8]+0.7071067811865475*fl[4])*wv+(0.3535533905932737*fl[13]+0.2041241452319315*fl[10])*dv; 
  Ghat_l[4] = (1.224744871391589*fl[11]+0.7071067811865475*fl[7])*wv+(0.3162277660168379*fl[20]+0.1825741858350554*fl[18]+0.3535533905932737*fl[5]+0.2041241452319315*fl[2])*dv; 
  Ghat_l[5] = (1.224744871391589*fl[12]+0.7071067811865475*fl[9])*wv+(0.3535533905932737*fl[15]+0.2041241452319315*fl[14])*dv; 
  Ghat_l[6] = (1.224744871391589*fl[13]+0.7071067811865475*fl[10])*wv+(0.3162277660168379*fl[21]+0.1825741858350554*fl[19]+0.3535533905932737*fl[8]+0.2041241452319315*fl[4])*dv; 
  Ghat_l[7] = (1.224744871391589*fl[15]+0.7071067811865475*fl[14])*wv+(0.3162277660168379*fl[23]+0.1825741858350554*fl[22]+0.3535533905932737*fl[12]+0.2041241452319315*fl[9])*dv; 
  Ghat_l[8] = (1.224744871391589*fl[17]+0.7071067811865475*fl[16])*wv+(0.3162277660168379*fl[6]+0.1825741858350554*fl[3])*dv; 
  Ghat_l[9] = (1.224744871391589*fl[20]+0.7071067811865475*fl[18])*wv+(0.3162277660168379*fl[11]+0.1825741858350554*fl[7])*dv; 
  Ghat_l[10] = (1.224744871391589*fl[21]+0.7071067811865475*fl[19])*wv+(0.3162277660168379*fl[13]+0.1825741858350554*fl[10])*dv; 
  Ghat_l[11] = (1.224744871391589*fl[23]+0.7071067811865475*fl[22])*wv+(0.3162277660168379*fl[15]+0.1825741858350554*fl[14])*dv; 
  Ghat_l[12] = (1.224744871391589*fl[25]+0.7071067811865475*fl[24])*wv+(0.3535533905932737*fl[29]+0.2041241452319315*fl[27])*dv; 
  Ghat_l[13] = (1.224744871391589*fl[28]+0.7071067811865475*fl[26])*wv+(0.3535533905932737*fl[31]+0.2041241452319315*fl[30])*dv; 
  Ghat_l[14] = (1.224744871391589*fl[29]+0.7071067811865475*fl[27])*wv+(0.3535533905932737*fl[25]+0.2041241452319315*fl[24])*dv; 
  Ghat_l[15] = (1.224744871391589*fl[31]+0.7071067811865475*fl[30])*wv+(0.3535533905932737*fl[28]+0.2041241452319315*fl[26])*dv; 

  } else { 

  Ghat_r[0] = -0.1178511301977579*((10.39230484541326*fr[1]-6.0*fr[0])*wv+(3.0*fr[6]-1.732050807568877*fr[3])*dv); 
  Ghat_r[1] = -0.1178511301977579*((10.39230484541326*fr[5]-6.0*fr[2])*wv+(3.0*fr[11]-1.732050807568877*fr[7])*dv); 
  Ghat_r[2] = -0.02357022603955158*((51.96152422706631*fr[6]-30.0*fr[3])*wv+(13.41640786499874*fr[17]-7.745966692414834*fr[16]+15.0*fr[1]-8.660254037844386*fr[0])*dv); 
  Ghat_r[3] = -0.1178511301977579*((10.39230484541326*fr[8]-6.0*fr[4])*wv+(3.0*fr[13]-1.732050807568877*fr[10])*dv); 
  Ghat_r[4] = -0.02357022603955158*((51.96152422706631*fr[11]-30.0*fr[7])*wv+(13.41640786499874*fr[20]-7.745966692414834*fr[18]+15.0*fr[5]-8.660254037844386*fr[2])*dv); 
  Ghat_r[5] = -0.1178511301977579*((10.39230484541326*fr[12]-6.0*fr[9])*wv+(3.0*fr[15]-1.732050807568877*fr[14])*dv); 
  Ghat_r[6] = -0.02357022603955158*((51.96152422706631*fr[13]-30.0*fr[10])*wv+(13.41640786499874*fr[21]-7.745966692414834*fr[19]+15.0*fr[8]-8.660254037844386*fr[4])*dv); 
  Ghat_r[7] = -0.02357022603955158*((51.96152422706631*fr[15]-30.0*fr[14])*wv+(13.41640786499874*fr[23]-7.745966692414834*fr[22]+15.0*fr[12]-8.660254037844386*fr[9])*dv); 
  Ghat_r[8] = -0.04714045207910316*((25.98076211353316*fr[17]-15.0*fr[16])*wv+(6.708203932499369*fr[6]-3.872983346207417*fr[3])*dv); 
  Ghat_r[9] = -0.04714045207910316*((25.98076211353316*fr[20]-15.0*fr[18])*wv+(6.708203932499369*fr[11]-3.872983346207417*fr[7])*dv); 
  Ghat_r[10] = -0.04714045207910316*((25.98076211353316*fr[21]-15.0*fr[19])*wv+(6.708203932499369*fr[13]-3.872983346207417*fr[10])*dv); 
  Ghat_r[11] = -0.04714045207910316*((25.98076211353316*fr[23]-15.0*fr[22])*wv+(6.708203932499369*fr[15]-3.872983346207417*fr[14])*dv); 
  Ghat_r[12] = -0.02357022603955158*((51.96152422706632*fr[25]-30.0*fr[24])*wv+(15.0*fr[29]-8.660254037844387*fr[27])*dv); 
  Ghat_r[13] = -0.02357022603955158*((51.96152422706632*fr[28]-30.0*fr[26])*wv+(15.0*fr[31]-8.660254037844387*fr[30])*dv); 
  Ghat_r[14] = -0.02357022603955158*((51.96152422706632*fr[29]-30.0*fr[27])*wv+(15.0*fr[25]-8.660254037844387*fr[24])*dv); 
  Ghat_r[15] = -0.02357022603955158*((51.96152422706632*fr[31]-30.0*fr[30])*wv+(15.0*fr[28]-8.660254037844387*fr[26])*dv); 

  Ghat_l[0] = -0.1178511301977579*((10.39230484541326*fc[1]-6.0*fc[0])*wv+(3.0*fc[6]-1.732050807568877*fc[3])*dv); 
  Ghat_l[1] = -0.1178511301977579*((10.39230484541326*fc[5]-6.0*fc[2])*wv+(3.0*fc[11]-1.732050807568877*fc[7])*dv); 
  Ghat_l[2] = -0.02357022603955158*((51.96152422706631*fc[6]-30.0*fc[3])*wv+(13.41640786499874*fc[17]-7.745966692414834*fc[16]+15.0*fc[1]-8.660254037844386*fc[0])*dv); 
  Ghat_l[3] = -0.1178511301977579*((10.39230484541326*fc[8]-6.0*fc[4])*wv+(3.0*fc[13]-1.732050807568877*fc[10])*dv); 
  Ghat_l[4] = -0.02357022603955158*((51.96152422706631*fc[11]-30.0*fc[7])*wv+(13.41640786499874*fc[20]-7.745966692414834*fc[18]+15.0*fc[5]-8.660254037844386*fc[2])*dv); 
  Ghat_l[5] = -0.1178511301977579*((10.39230484541326*fc[12]-6.0*fc[9])*wv+(3.0*fc[15]-1.732050807568877*fc[14])*dv); 
  Ghat_l[6] = -0.02357022603955158*((51.96152422706631*fc[13]-30.0*fc[10])*wv+(13.41640786499874*fc[21]-7.745966692414834*fc[19]+15.0*fc[8]-8.660254037844386*fc[4])*dv); 
  Ghat_l[7] = -0.02357022603955158*((51.96152422706631*fc[15]-30.0*fc[14])*wv+(13.41640786499874*fc[23]-7.745966692414834*fc[22]+15.0*fc[12]-8.660254037844386*fc[9])*dv); 
  Ghat_l[8] = -0.04714045207910316*((25.98076211353316*fc[17]-15.0*fc[16])*wv+(6.708203932499369*fc[6]-3.872983346207417*fc[3])*dv); 
  Ghat_l[9] = -0.04714045207910316*((25.98076211353316*fc[20]-15.0*fc[18])*wv+(6.708203932499369*fc[11]-3.872983346207417*fc[7])*dv); 
  Ghat_l[10] = -0.04714045207910316*((25.98076211353316*fc[21]-15.0*fc[19])*wv+(6.708203932499369*fc[13]-3.872983346207417*fc[10])*dv); 
  Ghat_l[11] = -0.04714045207910316*((25.98076211353316*fc[23]-15.0*fc[22])*wv+(6.708203932499369*fc[15]-3.872983346207417*fc[14])*dv); 
  Ghat_l[12] = -0.02357022603955158*((51.96152422706632*fc[25]-30.0*fc[24])*wv+(15.0*fc[29]-8.660254037844387*fc[27])*dv); 
  Ghat_l[13] = -0.02357022603955158*((51.96152422706632*fc[28]-30.0*fc[26])*wv+(15.0*fc[31]-8.660254037844387*fc[30])*dv); 
  Ghat_l[14] = -0.02357022603955158*((51.96152422706632*fc[29]-30.0*fc[27])*wv+(15.0*fc[25]-8.660254037844387*fc[24])*dv); 
  Ghat_l[15] = -0.02357022603955158*((51.96152422706632*fc[31]-30.0*fc[30])*wv+(15.0*fc[28]-8.660254037844387*fc[26])*dv); 

  } 
  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dx10; 
  out[1] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dx10; 
  out[2] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dx10; 
  out[3] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dx10; 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dx10; 
  out[5] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dx10; 
  out[6] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dx10; 
  out[7] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dx10; 
  out[8] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dx10; 
  out[9] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dx10; 
  out[10] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dx10; 
  out[11] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dx10; 
  out[12] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dx10; 
  out[13] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dx10; 
  out[14] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dx10; 
  out[15] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dx10; 
  out[16] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*dx10; 
  out[17] += -1.224744871391589*(Ghat_r[8]+Ghat_l[8])*dx10; 
  out[18] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*dx10; 
  out[19] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*dx10; 
  out[20] += -1.224744871391589*(Ghat_r[9]+Ghat_l[9])*dx10; 
  out[21] += -1.224744871391589*(Ghat_r[10]+Ghat_l[10])*dx10; 
  out[22] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*dx10; 
  out[23] += -1.224744871391589*(Ghat_r[11]+Ghat_l[11])*dx10; 
  out[24] += (0.7071067811865475*Ghat_l[12]-0.7071067811865475*Ghat_r[12])*dx10; 
  out[25] += -1.224744871391589*(Ghat_r[12]+Ghat_l[12])*dx10; 
  out[26] += (0.7071067811865475*Ghat_l[13]-0.7071067811865475*Ghat_r[13])*dx10; 
  out[27] += (0.7071067811865475*Ghat_l[14]-0.7071067811865475*Ghat_r[14])*dx10; 
  out[28] += -1.224744871391589*(Ghat_r[13]+Ghat_l[13])*dx10; 
  out[29] += -1.224744871391589*(Ghat_r[14]+Ghat_l[14])*dx10; 
  out[30] += (0.7071067811865475*Ghat_l[15]-0.7071067811865475*Ghat_r[15])*dx10; 
  out[31] += -1.224744871391589*(Ghat_r[15]+Ghat_l[15])*dx10; 

  return 0.;

} 
