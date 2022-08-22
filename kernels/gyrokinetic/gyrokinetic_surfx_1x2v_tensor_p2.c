#include <gkyl_gyrokinetic_kernels.h>
GKYL_CU_DH void gyrokinetic_surfx_1x2v_tensor_p2(const double *w, const double *dxv, const double q_, const double m_, const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi, const double *apar, const double *apardot, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // q_,m_: species charge and mass.
  // bmag: magnetic field amplitude.
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // cmag: coefficient multiplying parallel gradient.
  // b_i: covariant components of the field aligned unit vector.
  // phi: electrostatic potential .
  // apar: parallel component of magnetic vector potential.
  // apardot: time derivative of Apar.
  // fl,fc,fr: distribution function in left, center and right cells.
  // out: output increment in center cell.

  double wx = w[0];
  double rdx2 = 2.0/dxv[0];
  double wvpar = w[1];
  double rdvpar2 = 2.0/dxv[1];
  double wmu = w[2];
  double rdmu2 = 2.0/dxv[2];

  double wxSq = w[0]*w[0];
  double rdx2Sq = rdx2*rdx2;
  double wvparSq = w[1]*w[1];
  double rdvpar2Sq = rdvpar2*rdvpar2;
  double wmuSq = w[2]*w[2];
  double rdmu2Sq = rdmu2*rdmu2;

  const double *b_x = &b_i[0];
  const double *b_y = &b_i[3];
  const double *b_z = &b_i[6];

  double hamil[27] = {0.}; 
  hamil[0] = (0.2357022603955158*(3.0*rdvpar2Sq*(2.0*m_*wvparSq+2.828427124746191*(bmag[0]*wmu+phi[0]*q_))+2.0*m_))/rdvpar2Sq; 
  hamil[1] = 2.0*(bmag[1]*wmu+phi[1]*q_); 
  hamil[2] = (1.632993161855453*m_*wvpar)/rdvpar2; 
  hamil[3] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[5] = (1.154700538379252*bmag[1])/rdmu2; 
  hamil[7] = 2.0*(bmag[2]*wmu+phi[2]*q_); 
  hamil[8] = (0.421637021355784*m_)/rdvpar2Sq; 
  hamil[13] = (1.154700538379251*bmag[2])/rdmu2; 

  double BstarZdBmag[27] = {0.}; 
  BstarZdBmag[0] = (1.414213562373095*(1.732050807568877*(2.23606797749979*jacobtot_inv[1]*b_y[2]+jacobtot_inv[0]*b_y[1])*m_*rdx2*wvpar+(cmag[2]*jacobtot_inv[2]+cmag[1]*jacobtot_inv[1]+cmag[0]*jacobtot_inv[0])*q_))/q_; 
  BstarZdBmag[1] = (0.2828427124746191*(5.0*(1.732050807568877*(2.0*b_y[2]*jacobtot_inv[2]+b_y[1]*jacobtot_inv[1])*m_*rdx2*wvpar+(cmag[0]*jacobtot_inv[1]+jacobtot_inv[0]*cmag[1])*q_)+2.23606797749979*(8.660254037844386*jacobtot_inv[0]*b_y[2]*m_*rdx2*wvpar+2.0*(cmag[1]*jacobtot_inv[2]+jacobtot_inv[1]*cmag[2])*q_)))/q_; 
  BstarZdBmag[2] = (1.414213562373095*(2.23606797749979*jacobtot_inv[1]*b_y[2]+jacobtot_inv[0]*b_y[1])*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[4] = (1.414213562373095*(b_y[2]*(2.0*jacobtot_inv[2]+2.23606797749979*jacobtot_inv[0])+b_y[1]*jacobtot_inv[1])*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[7] = (0.04040610178208843*(60.6217782649107*(b_y[1]*jacobtot_inv[2]+2.0*jacobtot_inv[1]*b_y[2])*m_*rdx2*wvpar+(4.47213595499958*(5.0*cmag[2]*jacobtot_inv[2]+7.0*cmag[1]*jacobtot_inv[1])+35.0*(cmag[0]*jacobtot_inv[2]+jacobtot_inv[0]*cmag[2]))*q_))/q_; 
  BstarZdBmag[11] = (1.414213562373095*(b_y[1]*jacobtot_inv[2]+2.0*jacobtot_inv[1]*b_y[2])*m_*rdx2)/(q_*rdvpar2); 

  double alphaL[9] = {0.}; 
  alphaL[0] = (0.25*(8.660254037844387*hamil[8]*BstarZdBmag[11]+2.23606797749979*(1.732050807568877*(BstarZdBmag[2]*hamil[8]+hamil[2]*BstarZdBmag[7])-3.0*BstarZdBmag[4]*hamil[8])+(1.732050807568877*BstarZdBmag[0]-3.0*BstarZdBmag[1])*hamil[2])*rdvpar2)/m_; 
  alphaL[1] = (0.25*(3.872983346207417*hamil[2]*BstarZdBmag[11]+1.732050807568877*(5.0*BstarZdBmag[7]*hamil[8]+BstarZdBmag[2]*hamil[2])+2.23606797749979*(1.732050807568877*BstarZdBmag[0]-3.0*BstarZdBmag[1])*hamil[8]-3.0*hamil[2]*BstarZdBmag[4])*rdvpar2)/m_; 
  alphaL[4] = (0.5*hamil[8]*(3.872983346207417*BstarZdBmag[11]-3.0*BstarZdBmag[4]+1.732050807568877*BstarZdBmag[2])*rdvpar2)/m_; 

  double alphaR[9] = {0.}; 
  alphaR[0] = (0.25*(8.660254037844387*hamil[8]*BstarZdBmag[11]+2.23606797749979*(1.732050807568877*(BstarZdBmag[2]*hamil[8]+hamil[2]*BstarZdBmag[7])+3.0*BstarZdBmag[4]*hamil[8])+(3.0*BstarZdBmag[1]+1.732050807568877*BstarZdBmag[0])*hamil[2])*rdvpar2)/m_; 
  alphaR[1] = (0.25*(3.872983346207417*hamil[2]*BstarZdBmag[11]+1.732050807568877*(5.0*BstarZdBmag[7]*hamil[8]+BstarZdBmag[2]*hamil[2])+2.23606797749979*(3.0*BstarZdBmag[1]+1.732050807568877*BstarZdBmag[0])*hamil[8]+3.0*hamil[2]*BstarZdBmag[4])*rdvpar2)/m_; 
  alphaR[4] = (0.5*hamil[8]*(3.872983346207417*BstarZdBmag[11]+3.0*BstarZdBmag[4]+1.732050807568877*BstarZdBmag[2])*rdvpar2)/m_; 

  double fUpOrdL[9] = {0.};
  if (alphaL[4]-1.5*alphaL[1]+1.118033988749897*alphaL[0] > 0.) {
    fUpOrdL[0] = 0.6324555320336759*fl[26]+0.4898979485566357*fl[25]-0.9486832980505122*fl[24]-0.9486832980505137*fl[23]+0.2828427124746191*fl[22]+0.7071067811865475*(fl[21]+fl[20])-0.7348469228349525*fl[19]-0.7348469228349533*fl[18]+1.42302494707577*fl[17]-0.4242640687119281*fl[16]+0.5477225575051661*fl[15]-0.4242640687119285*fl[14]-1.060660171779821*fl[13]+0.5477225575051661*fl[12]-1.060660171779821*fl[11]+1.10227038425243*fl[10]+0.3162277660168379*(fl[9]+fl[8])+0.7905694150420947*fl[7]+0.6363961030678926*fl[6]-0.8215838362577489*(fl[5]+fl[4])-0.4743416490252568*(fl[3]+fl[2])+0.6123724356957944*fl[1]+0.3535533905932737*fl[0]; 
  } else {
    fUpOrdL[0] = 0.6324555320336759*fc[26]-0.4898979485566357*fc[25]-0.9486832980505122*fc[24]-0.9486832980505137*fc[23]+0.2828427124746191*fc[22]+0.7071067811865475*(fc[21]+fc[20])+0.7348469228349525*fc[19]+0.7348469228349533*fc[18]+1.42302494707577*fc[17]-0.4242640687119281*fc[16]-0.5477225575051661*fc[15]-0.4242640687119285*fc[14]-1.060660171779821*fc[13]-0.5477225575051661*fc[12]-1.060660171779821*fc[11]-1.10227038425243*fc[10]+0.3162277660168379*(fc[9]+fc[8])+0.7905694150420947*fc[7]+0.6363961030678926*fc[6]+0.8215838362577489*(fc[5]+fc[4])-0.4743416490252568*(fc[3]+fc[2])-0.6123724356957944*fc[1]+0.3535533905932737*fc[0]; 
  }
  if (alphaL[4]-1.5*alphaL[1]+1.118033988749897*alphaL[0] > 0.) {
    fUpOrdL[1] = (-0.7905694150420947*fl[26])-0.6123724356957944*fl[25]+1.185854122563141*fl[24]-0.3535533905932737*fl[22]-0.883883476483184*fl[21]+0.7071067811865475*fl[20]+0.9185586535436913*fl[19]+0.5303300858899104*fl[16]-0.6846531968814574*fl[15]+0.5477225575051661*fl[12]-1.060660171779821*fl[11]-0.3952847075210473*fl[9]+0.3162277660168379*fl[8]+0.7905694150420947*fl[7]-0.8215838362577489*fl[4]-0.4743416490252568*fl[2]+0.6123724356957944*fl[1]+0.3535533905932737*fl[0]; 
  } else {
    fUpOrdL[1] = (-0.7905694150420947*fc[26])+0.6123724356957944*fc[25]+1.185854122563141*fc[24]-0.3535533905932737*fc[22]-0.883883476483184*fc[21]+0.7071067811865475*fc[20]-0.9185586535436913*fc[19]+0.5303300858899104*fc[16]+0.6846531968814574*fc[15]-0.5477225575051661*fc[12]-1.060660171779821*fc[11]-0.3952847075210473*fc[9]+0.3162277660168379*fc[8]+0.7905694150420947*fc[7]+0.8215838362577489*fc[4]-0.4743416490252568*fc[2]-0.6123724356957944*fc[1]+0.3535533905932737*fc[0]; 
  }
  if (alphaL[4]-1.5*alphaL[1]+1.118033988749897*alphaL[0] > 0.) {
    fUpOrdL[2] = 0.6324555320336759*fl[26]+0.4898979485566357*fl[25]-0.9486832980505122*fl[24]+0.9486832980505137*fl[23]+0.2828427124746191*fl[22]+0.7071067811865475*(fl[21]+fl[20])-0.7348469228349525*fl[19]+0.7348469228349533*fl[18]-1.42302494707577*fl[17]-0.4242640687119281*fl[16]+0.5477225575051661*fl[15]+0.4242640687119285*fl[14]+1.060660171779821*fl[13]+0.5477225575051661*fl[12]-1.060660171779821*fl[11]-1.10227038425243*fl[10]+0.3162277660168379*(fl[9]+fl[8])+0.7905694150420947*fl[7]-0.6363961030678926*fl[6]+0.8215838362577489*fl[5]-0.8215838362577489*fl[4]+0.4743416490252568*fl[3]-0.4743416490252568*fl[2]+0.6123724356957944*fl[1]+0.3535533905932737*fl[0]; 
  } else {
    fUpOrdL[2] = 0.6324555320336759*fc[26]-0.4898979485566357*fc[25]-0.9486832980505122*fc[24]+0.9486832980505137*fc[23]+0.2828427124746191*fc[22]+0.7071067811865475*(fc[21]+fc[20])+0.7348469228349525*fc[19]-0.7348469228349533*fc[18]-1.42302494707577*fc[17]-0.4242640687119281*fc[16]-0.5477225575051661*fc[15]+0.4242640687119285*fc[14]+1.060660171779821*fc[13]-0.5477225575051661*fc[12]-1.060660171779821*fc[11]+1.10227038425243*fc[10]+0.3162277660168379*(fc[9]+fc[8])+0.7905694150420947*fc[7]-0.6363961030678926*fc[6]-0.8215838362577489*fc[5]+0.8215838362577489*fc[4]+0.4743416490252568*fc[3]-0.4743416490252568*fc[2]-0.6123724356957944*fc[1]+0.3535533905932737*fc[0]; 
  }
  if (alphaL[0]-1.118033988749893*alphaL[4] > 0.) {
    fUpOrdL[3] = (-0.7905694150420947*fl[26])-0.6123724356957944*fl[25]+1.185854122563142*fl[23]-0.3535533905932737*fl[22]+0.7071067811865475*fl[21]-0.883883476483184*fl[20]+0.9185586535436913*fl[18]+0.5477225575051661*fl[15]+0.5303300858899104*fl[14]-1.060660171779821*fl[13]-0.6846531968814574*fl[12]+0.3162277660168379*fl[9]-0.3952847075210473*fl[8]+0.7905694150420947*fl[7]-0.8215838362577489*fl[5]-0.4743416490252568*fl[3]+0.6123724356957944*fl[1]+0.3535533905932737*fl[0]; 
  } else {
    fUpOrdL[3] = (-0.7905694150420947*fc[26])+0.6123724356957944*fc[25]+1.185854122563142*fc[23]-0.3535533905932737*fc[22]+0.7071067811865475*fc[21]-0.883883476483184*fc[20]-0.9185586535436913*fc[18]-0.5477225575051661*fc[15]+0.5303300858899104*fc[14]-1.060660171779821*fc[13]+0.6846531968814574*fc[12]+0.3162277660168379*fc[9]-0.3952847075210473*fc[8]+0.7905694150420947*fc[7]+0.8215838362577489*fc[5]-0.4743416490252568*fc[3]-0.6123724356957944*fc[1]+0.3535533905932737*fc[0]; 
  }
  if (alphaL[0]-1.118033988749893*alphaL[4] > 0.) {
    fUpOrdL[4] = 0.9882117688026183*fl[26]+0.7654655446197428*fl[25]+0.441941738241592*fl[22]-0.883883476483184*(fl[21]+fl[20])-0.6846531968814574*(fl[15]+fl[12])-0.3952847075210473*(fl[9]+fl[8])+0.7905694150420947*fl[7]+0.6123724356957944*fl[1]+0.3535533905932737*fl[0]; 
  } else {
    fUpOrdL[4] = 0.9882117688026183*fc[26]-0.7654655446197428*fc[25]+0.441941738241592*fc[22]-0.883883476483184*(fc[21]+fc[20])+0.6846531968814574*(fc[15]+fc[12])-0.3952847075210473*(fc[9]+fc[8])+0.7905694150420947*fc[7]-0.6123724356957944*fc[1]+0.3535533905932737*fc[0]; 
  }
  if (alphaL[0]-1.118033988749893*alphaL[4] > 0.) {
    fUpOrdL[5] = (-0.7905694150420947*fl[26])-0.6123724356957944*fl[25]-1.185854122563142*fl[23]-0.3535533905932737*fl[22]+0.7071067811865475*fl[21]-0.883883476483184*fl[20]-0.9185586535436913*fl[18]+0.5477225575051661*fl[15]-0.5303300858899104*fl[14]+1.060660171779821*fl[13]-0.6846531968814574*fl[12]+0.3162277660168379*fl[9]-0.3952847075210473*fl[8]+0.7905694150420947*fl[7]+0.8215838362577489*fl[5]+0.4743416490252568*fl[3]+0.6123724356957944*fl[1]+0.3535533905932737*fl[0]; 
  } else {
    fUpOrdL[5] = (-0.7905694150420947*fc[26])+0.6123724356957944*fc[25]-1.185854122563142*fc[23]-0.3535533905932737*fc[22]+0.7071067811865475*fc[21]-0.883883476483184*fc[20]+0.9185586535436913*fc[18]-0.5477225575051661*fc[15]-0.5303300858899104*fc[14]+1.060660171779821*fc[13]+0.6846531968814574*fc[12]+0.3162277660168379*fc[9]-0.3952847075210473*fc[8]+0.7905694150420947*fc[7]-0.8215838362577489*fc[5]+0.4743416490252568*fc[3]-0.6123724356957944*fc[1]+0.3535533905932737*fc[0]; 
  }
  if (alphaL[4]+1.5*alphaL[1]+1.118033988749897*alphaL[0] > 0.) {
    fUpOrdL[6] = 0.6324555320336759*fl[26]+0.4898979485566357*fl[25]+0.9486832980505122*fl[24]-0.9486832980505137*fl[23]+0.2828427124746191*fl[22]+0.7071067811865475*(fl[21]+fl[20])+0.7348469228349525*fl[19]-0.7348469228349533*fl[18]-1.42302494707577*fl[17]+0.4242640687119281*fl[16]+0.5477225575051661*fl[15]-0.4242640687119285*fl[14]-1.060660171779821*fl[13]+0.5477225575051661*fl[12]+1.060660171779821*fl[11]-1.10227038425243*fl[10]+0.3162277660168379*(fl[9]+fl[8])+0.7905694150420947*fl[7]-0.6363961030678926*fl[6]-0.8215838362577489*fl[5]+0.8215838362577489*fl[4]-0.4743416490252568*fl[3]+0.4743416490252568*fl[2]+0.6123724356957944*fl[1]+0.3535533905932737*fl[0]; 
  } else {
    fUpOrdL[6] = 0.6324555320336759*fc[26]-0.4898979485566357*fc[25]+0.9486832980505122*fc[24]-0.9486832980505137*fc[23]+0.2828427124746191*fc[22]+0.7071067811865475*(fc[21]+fc[20])-0.7348469228349525*fc[19]+0.7348469228349533*fc[18]-1.42302494707577*fc[17]+0.4242640687119281*fc[16]-0.5477225575051661*fc[15]-0.4242640687119285*fc[14]-1.060660171779821*fc[13]-0.5477225575051661*fc[12]+1.060660171779821*fc[11]+1.10227038425243*fc[10]+0.3162277660168379*(fc[9]+fc[8])+0.7905694150420947*fc[7]-0.6363961030678926*fc[6]+0.8215838362577489*fc[5]-0.8215838362577489*fc[4]-0.4743416490252568*fc[3]+0.4743416490252568*fc[2]-0.6123724356957944*fc[1]+0.3535533905932737*fc[0]; 
  }
  if (alphaL[4]+1.5*alphaL[1]+1.118033988749897*alphaL[0] > 0.) {
    fUpOrdL[7] = (-0.7905694150420947*fl[26])-0.6123724356957944*fl[25]-1.185854122563141*fl[24]-0.3535533905932737*fl[22]-0.883883476483184*fl[21]+0.7071067811865475*fl[20]-0.9185586535436913*fl[19]-0.5303300858899104*fl[16]-0.6846531968814574*fl[15]+0.5477225575051661*fl[12]+1.060660171779821*fl[11]-0.3952847075210473*fl[9]+0.3162277660168379*fl[8]+0.7905694150420947*fl[7]+0.8215838362577489*fl[4]+0.4743416490252568*fl[2]+0.6123724356957944*fl[1]+0.3535533905932737*fl[0]; 
  } else {
    fUpOrdL[7] = (-0.7905694150420947*fc[26])+0.6123724356957944*fc[25]-1.185854122563141*fc[24]-0.3535533905932737*fc[22]-0.883883476483184*fc[21]+0.7071067811865475*fc[20]+0.9185586535436913*fc[19]-0.5303300858899104*fc[16]+0.6846531968814574*fc[15]-0.5477225575051661*fc[12]+1.060660171779821*fc[11]-0.3952847075210473*fc[9]+0.3162277660168379*fc[8]+0.7905694150420947*fc[7]-0.8215838362577489*fc[4]+0.4743416490252568*fc[2]-0.6123724356957944*fc[1]+0.3535533905932737*fc[0]; 
  }
  if (alphaL[4]+1.5*alphaL[1]+1.118033988749897*alphaL[0] > 0.) {
    fUpOrdL[8] = 0.6324555320336759*fl[26]+0.4898979485566357*fl[25]+0.9486832980505122*fl[24]+0.9486832980505137*fl[23]+0.2828427124746191*fl[22]+0.7071067811865475*(fl[21]+fl[20])+0.7348469228349525*fl[19]+0.7348469228349533*fl[18]+1.42302494707577*fl[17]+0.4242640687119281*fl[16]+0.5477225575051661*fl[15]+0.4242640687119285*fl[14]+1.060660171779821*fl[13]+0.5477225575051661*fl[12]+1.060660171779821*fl[11]+1.10227038425243*fl[10]+0.3162277660168379*(fl[9]+fl[8])+0.7905694150420947*fl[7]+0.6363961030678926*fl[6]+0.8215838362577489*(fl[5]+fl[4])+0.4743416490252568*(fl[3]+fl[2])+0.6123724356957944*fl[1]+0.3535533905932737*fl[0]; 
  } else {
    fUpOrdL[8] = 0.6324555320336759*fc[26]-0.4898979485566357*fc[25]+0.9486832980505122*fc[24]+0.9486832980505137*fc[23]+0.2828427124746191*fc[22]+0.7071067811865475*(fc[21]+fc[20])-0.7348469228349525*fc[19]-0.7348469228349533*fc[18]+1.42302494707577*fc[17]+0.4242640687119281*fc[16]-0.5477225575051661*fc[15]+0.4242640687119285*fc[14]+1.060660171779821*fc[13]-0.5477225575051661*fc[12]+1.060660171779821*fc[11]-1.10227038425243*fc[10]+0.3162277660168379*(fc[9]+fc[8])+0.7905694150420947*fc[7]+0.6363961030678926*fc[6]-0.8215838362577489*(fc[5]+fc[4])+0.4743416490252568*(fc[3]+fc[2])-0.6123724356957944*fc[1]+0.3535533905932737*fc[0]; 
  }

  double fUpL[9] = {0.};
  fUpL[0] = 0.154320987654321*fUpOrdL[8]+0.2469135802469136*fUpOrdL[7]+0.154320987654321*fUpOrdL[6]+0.2469135802469136*fUpOrdL[5]+0.3950617283950617*fUpOrdL[4]+0.2469135802469136*fUpOrdL[3]+0.154320987654321*fUpOrdL[2]+0.2469135802469136*fUpOrdL[1]+0.154320987654321*fUpOrdL[0]; 
  fUpL[1] = 0.2070433312499806*fUpOrdL[8]+0.3312693299999688*fUpOrdL[7]+0.2070433312499806*fUpOrdL[6]-0.2070433312499806*fUpOrdL[2]-0.3312693299999688*fUpOrdL[1]-0.2070433312499806*fUpOrdL[0]; 
  fUpL[2] = 0.2070433312499806*fUpOrdL[8]-0.2070433312499806*fUpOrdL[6]+0.3312693299999688*fUpOrdL[5]-0.3312693299999688*fUpOrdL[3]+0.2070433312499806*fUpOrdL[2]-0.2070433312499806*fUpOrdL[0]; 
  fUpL[3] = 0.2777777777777778*fUpOrdL[8]-0.2777777777777778*fUpOrdL[6]-0.2777777777777778*fUpOrdL[2]+0.2777777777777778*fUpOrdL[0]; 
  fUpL[4] = 0.138028887499987*fUpOrdL[8]+0.2208462199999792*fUpOrdL[7]+0.138028887499987*fUpOrdL[6]-0.2760577749999741*fUpOrdL[5]-0.4416924399999584*fUpOrdL[4]-0.2760577749999741*fUpOrdL[3]+0.138028887499987*fUpOrdL[2]+0.2208462199999792*fUpOrdL[1]+0.138028887499987*fUpOrdL[0]; 
  fUpL[5] = 0.138028887499987*fUpOrdL[8]-0.2760577749999741*fUpOrdL[7]+0.138028887499987*fUpOrdL[6]+0.2208462199999792*fUpOrdL[5]-0.4416924399999584*fUpOrdL[4]+0.2208462199999792*fUpOrdL[3]+0.138028887499987*fUpOrdL[2]-0.2760577749999741*fUpOrdL[1]+0.138028887499987*fUpOrdL[0]; 
  fUpL[6] = 0.1851851851851853*fUpOrdL[8]-0.1851851851851853*fUpOrdL[6]-0.3703703703703705*fUpOrdL[5]+0.3703703703703705*fUpOrdL[3]+0.1851851851851853*fUpOrdL[2]-0.1851851851851853*fUpOrdL[0]; 
  fUpL[7] = 0.1851851851851853*fUpOrdL[8]-0.3703703703703705*fUpOrdL[7]+0.1851851851851853*fUpOrdL[6]-0.1851851851851853*fUpOrdL[2]+0.3703703703703705*fUpOrdL[1]-0.1851851851851853*fUpOrdL[0]; 
  fUpL[8] = 0.1234567901234568*fUpOrdL[8]-0.2469135802469136*fUpOrdL[7]+0.1234567901234568*fUpOrdL[6]-0.2469135802469136*fUpOrdL[5]+0.4938271604938271*fUpOrdL[4]-0.2469135802469136*fUpOrdL[3]+0.1234567901234568*fUpOrdL[2]-0.2469135802469136*fUpOrdL[1]+0.1234567901234568*fUpOrdL[0]; 

  double GhatL[27] = {0.}; 
  GhatL[0] = 0.5*alphaL[4]*fUpL[4]+0.5*alphaL[1]*fUpL[1]+0.5*alphaL[0]*fUpL[0]; 
  GhatL[1] = 0.4472135954999579*alphaL[1]*fUpL[4]+0.4472135954999579*fUpL[1]*alphaL[4]+0.5*alphaL[0]*fUpL[1]+0.5*fUpL[0]*alphaL[1]; 
  GhatL[2] = 0.5000000000000001*alphaL[4]*fUpL[6]+0.5*alphaL[1]*fUpL[3]+0.5*alphaL[0]*fUpL[2]; 
  GhatL[3] = 0.447213595499958*alphaL[1]*fUpL[6]+0.4472135954999579*fUpL[3]*alphaL[4]+0.5*alphaL[0]*fUpL[3]+0.5*alphaL[1]*fUpL[2]; 
  GhatL[4] = 0.31943828249997*alphaL[4]*fUpL[4]+0.5*alphaL[0]*fUpL[4]+0.5*fUpL[0]*alphaL[4]+0.4472135954999579*alphaL[1]*fUpL[1]; 
  GhatL[5] = 0.5*alphaL[4]*fUpL[8]+0.5000000000000001*alphaL[1]*fUpL[7]+0.5*alphaL[0]*fUpL[5]; 
  GhatL[6] = 0.31943828249997*alphaL[4]*fUpL[6]+0.5*alphaL[0]*fUpL[6]+0.5000000000000001*fUpL[2]*alphaL[4]+0.447213595499958*alphaL[1]*fUpL[3]; 
  GhatL[7] = 0.447213595499958*alphaL[1]*fUpL[8]+0.4472135954999579*alphaL[4]*fUpL[7]+0.5*alphaL[0]*fUpL[7]+0.5000000000000001*alphaL[1]*fUpL[5]; 
  GhatL[8] = 0.31943828249997*alphaL[4]*fUpL[8]+0.5*alphaL[0]*fUpL[8]+0.447213595499958*alphaL[1]*fUpL[7]+0.5*alphaL[4]*fUpL[5]; 

  double fUpOrdR[9] = {0.};
  if (alphaR[4]-1.5*alphaR[1]+1.118033988749897*alphaR[0] > 0.) {
    fUpOrdR[0] = 0.6324555320336759*fc[26]+0.4898979485566357*fc[25]-0.9486832980505122*fc[24]-0.9486832980505137*fc[23]+0.2828427124746191*fc[22]+0.7071067811865475*(fc[21]+fc[20])-0.7348469228349525*fc[19]-0.7348469228349533*fc[18]+1.42302494707577*fc[17]-0.4242640687119281*fc[16]+0.5477225575051661*fc[15]-0.4242640687119285*fc[14]-1.060660171779821*fc[13]+0.5477225575051661*fc[12]-1.060660171779821*fc[11]+1.10227038425243*fc[10]+0.3162277660168379*(fc[9]+fc[8])+0.7905694150420947*fc[7]+0.6363961030678926*fc[6]-0.8215838362577489*(fc[5]+fc[4])-0.4743416490252568*(fc[3]+fc[2])+0.6123724356957944*fc[1]+0.3535533905932737*fc[0]; 
  } else {
    fUpOrdR[0] = 0.6324555320336759*fr[26]-0.4898979485566357*fr[25]-0.9486832980505122*fr[24]-0.9486832980505137*fr[23]+0.2828427124746191*fr[22]+0.7071067811865475*(fr[21]+fr[20])+0.7348469228349525*fr[19]+0.7348469228349533*fr[18]+1.42302494707577*fr[17]-0.4242640687119281*fr[16]-0.5477225575051661*fr[15]-0.4242640687119285*fr[14]-1.060660171779821*fr[13]-0.5477225575051661*fr[12]-1.060660171779821*fr[11]-1.10227038425243*fr[10]+0.3162277660168379*(fr[9]+fr[8])+0.7905694150420947*fr[7]+0.6363961030678926*fr[6]+0.8215838362577489*(fr[5]+fr[4])-0.4743416490252568*(fr[3]+fr[2])-0.6123724356957944*fr[1]+0.3535533905932737*fr[0]; 
  }
  if (alphaR[4]-1.5*alphaR[1]+1.118033988749897*alphaR[0] > 0.) {
    fUpOrdR[1] = (-0.7905694150420947*fc[26])-0.6123724356957944*fc[25]+1.185854122563141*fc[24]-0.3535533905932737*fc[22]-0.883883476483184*fc[21]+0.7071067811865475*fc[20]+0.9185586535436913*fc[19]+0.5303300858899104*fc[16]-0.6846531968814574*fc[15]+0.5477225575051661*fc[12]-1.060660171779821*fc[11]-0.3952847075210473*fc[9]+0.3162277660168379*fc[8]+0.7905694150420947*fc[7]-0.8215838362577489*fc[4]-0.4743416490252568*fc[2]+0.6123724356957944*fc[1]+0.3535533905932737*fc[0]; 
  } else {
    fUpOrdR[1] = (-0.7905694150420947*fr[26])+0.6123724356957944*fr[25]+1.185854122563141*fr[24]-0.3535533905932737*fr[22]-0.883883476483184*fr[21]+0.7071067811865475*fr[20]-0.9185586535436913*fr[19]+0.5303300858899104*fr[16]+0.6846531968814574*fr[15]-0.5477225575051661*fr[12]-1.060660171779821*fr[11]-0.3952847075210473*fr[9]+0.3162277660168379*fr[8]+0.7905694150420947*fr[7]+0.8215838362577489*fr[4]-0.4743416490252568*fr[2]-0.6123724356957944*fr[1]+0.3535533905932737*fr[0]; 
  }
  if (alphaR[4]-1.5*alphaR[1]+1.118033988749897*alphaR[0] > 0.) {
    fUpOrdR[2] = 0.6324555320336759*fc[26]+0.4898979485566357*fc[25]-0.9486832980505122*fc[24]+0.9486832980505137*fc[23]+0.2828427124746191*fc[22]+0.7071067811865475*(fc[21]+fc[20])-0.7348469228349525*fc[19]+0.7348469228349533*fc[18]-1.42302494707577*fc[17]-0.4242640687119281*fc[16]+0.5477225575051661*fc[15]+0.4242640687119285*fc[14]+1.060660171779821*fc[13]+0.5477225575051661*fc[12]-1.060660171779821*fc[11]-1.10227038425243*fc[10]+0.3162277660168379*(fc[9]+fc[8])+0.7905694150420947*fc[7]-0.6363961030678926*fc[6]+0.8215838362577489*fc[5]-0.8215838362577489*fc[4]+0.4743416490252568*fc[3]-0.4743416490252568*fc[2]+0.6123724356957944*fc[1]+0.3535533905932737*fc[0]; 
  } else {
    fUpOrdR[2] = 0.6324555320336759*fr[26]-0.4898979485566357*fr[25]-0.9486832980505122*fr[24]+0.9486832980505137*fr[23]+0.2828427124746191*fr[22]+0.7071067811865475*(fr[21]+fr[20])+0.7348469228349525*fr[19]-0.7348469228349533*fr[18]-1.42302494707577*fr[17]-0.4242640687119281*fr[16]-0.5477225575051661*fr[15]+0.4242640687119285*fr[14]+1.060660171779821*fr[13]-0.5477225575051661*fr[12]-1.060660171779821*fr[11]+1.10227038425243*fr[10]+0.3162277660168379*(fr[9]+fr[8])+0.7905694150420947*fr[7]-0.6363961030678926*fr[6]-0.8215838362577489*fr[5]+0.8215838362577489*fr[4]+0.4743416490252568*fr[3]-0.4743416490252568*fr[2]-0.6123724356957944*fr[1]+0.3535533905932737*fr[0]; 
  }
  if (alphaR[0]-1.118033988749893*alphaR[4] > 0.) {
    fUpOrdR[3] = (-0.7905694150420947*fc[26])-0.6123724356957944*fc[25]+1.185854122563142*fc[23]-0.3535533905932737*fc[22]+0.7071067811865475*fc[21]-0.883883476483184*fc[20]+0.9185586535436913*fc[18]+0.5477225575051661*fc[15]+0.5303300858899104*fc[14]-1.060660171779821*fc[13]-0.6846531968814574*fc[12]+0.3162277660168379*fc[9]-0.3952847075210473*fc[8]+0.7905694150420947*fc[7]-0.8215838362577489*fc[5]-0.4743416490252568*fc[3]+0.6123724356957944*fc[1]+0.3535533905932737*fc[0]; 
  } else {
    fUpOrdR[3] = (-0.7905694150420947*fr[26])+0.6123724356957944*fr[25]+1.185854122563142*fr[23]-0.3535533905932737*fr[22]+0.7071067811865475*fr[21]-0.883883476483184*fr[20]-0.9185586535436913*fr[18]-0.5477225575051661*fr[15]+0.5303300858899104*fr[14]-1.060660171779821*fr[13]+0.6846531968814574*fr[12]+0.3162277660168379*fr[9]-0.3952847075210473*fr[8]+0.7905694150420947*fr[7]+0.8215838362577489*fr[5]-0.4743416490252568*fr[3]-0.6123724356957944*fr[1]+0.3535533905932737*fr[0]; 
  }
  if (alphaR[0]-1.118033988749893*alphaR[4] > 0.) {
    fUpOrdR[4] = 0.9882117688026183*fc[26]+0.7654655446197428*fc[25]+0.441941738241592*fc[22]-0.883883476483184*(fc[21]+fc[20])-0.6846531968814574*(fc[15]+fc[12])-0.3952847075210473*(fc[9]+fc[8])+0.7905694150420947*fc[7]+0.6123724356957944*fc[1]+0.3535533905932737*fc[0]; 
  } else {
    fUpOrdR[4] = 0.9882117688026183*fr[26]-0.7654655446197428*fr[25]+0.441941738241592*fr[22]-0.883883476483184*(fr[21]+fr[20])+0.6846531968814574*(fr[15]+fr[12])-0.3952847075210473*(fr[9]+fr[8])+0.7905694150420947*fr[7]-0.6123724356957944*fr[1]+0.3535533905932737*fr[0]; 
  }
  if (alphaR[0]-1.118033988749893*alphaR[4] > 0.) {
    fUpOrdR[5] = (-0.7905694150420947*fc[26])-0.6123724356957944*fc[25]-1.185854122563142*fc[23]-0.3535533905932737*fc[22]+0.7071067811865475*fc[21]-0.883883476483184*fc[20]-0.9185586535436913*fc[18]+0.5477225575051661*fc[15]-0.5303300858899104*fc[14]+1.060660171779821*fc[13]-0.6846531968814574*fc[12]+0.3162277660168379*fc[9]-0.3952847075210473*fc[8]+0.7905694150420947*fc[7]+0.8215838362577489*fc[5]+0.4743416490252568*fc[3]+0.6123724356957944*fc[1]+0.3535533905932737*fc[0]; 
  } else {
    fUpOrdR[5] = (-0.7905694150420947*fr[26])+0.6123724356957944*fr[25]-1.185854122563142*fr[23]-0.3535533905932737*fr[22]+0.7071067811865475*fr[21]-0.883883476483184*fr[20]+0.9185586535436913*fr[18]-0.5477225575051661*fr[15]-0.5303300858899104*fr[14]+1.060660171779821*fr[13]+0.6846531968814574*fr[12]+0.3162277660168379*fr[9]-0.3952847075210473*fr[8]+0.7905694150420947*fr[7]-0.8215838362577489*fr[5]+0.4743416490252568*fr[3]-0.6123724356957944*fr[1]+0.3535533905932737*fr[0]; 
  }
  if (alphaR[4]+1.5*alphaR[1]+1.118033988749897*alphaR[0] > 0.) {
    fUpOrdR[6] = 0.6324555320336759*fc[26]+0.4898979485566357*fc[25]+0.9486832980505122*fc[24]-0.9486832980505137*fc[23]+0.2828427124746191*fc[22]+0.7071067811865475*(fc[21]+fc[20])+0.7348469228349525*fc[19]-0.7348469228349533*fc[18]-1.42302494707577*fc[17]+0.4242640687119281*fc[16]+0.5477225575051661*fc[15]-0.4242640687119285*fc[14]-1.060660171779821*fc[13]+0.5477225575051661*fc[12]+1.060660171779821*fc[11]-1.10227038425243*fc[10]+0.3162277660168379*(fc[9]+fc[8])+0.7905694150420947*fc[7]-0.6363961030678926*fc[6]-0.8215838362577489*fc[5]+0.8215838362577489*fc[4]-0.4743416490252568*fc[3]+0.4743416490252568*fc[2]+0.6123724356957944*fc[1]+0.3535533905932737*fc[0]; 
  } else {
    fUpOrdR[6] = 0.6324555320336759*fr[26]-0.4898979485566357*fr[25]+0.9486832980505122*fr[24]-0.9486832980505137*fr[23]+0.2828427124746191*fr[22]+0.7071067811865475*(fr[21]+fr[20])-0.7348469228349525*fr[19]+0.7348469228349533*fr[18]-1.42302494707577*fr[17]+0.4242640687119281*fr[16]-0.5477225575051661*fr[15]-0.4242640687119285*fr[14]-1.060660171779821*fr[13]-0.5477225575051661*fr[12]+1.060660171779821*fr[11]+1.10227038425243*fr[10]+0.3162277660168379*(fr[9]+fr[8])+0.7905694150420947*fr[7]-0.6363961030678926*fr[6]+0.8215838362577489*fr[5]-0.8215838362577489*fr[4]-0.4743416490252568*fr[3]+0.4743416490252568*fr[2]-0.6123724356957944*fr[1]+0.3535533905932737*fr[0]; 
  }
  if (alphaR[4]+1.5*alphaR[1]+1.118033988749897*alphaR[0] > 0.) {
    fUpOrdR[7] = (-0.7905694150420947*fc[26])-0.6123724356957944*fc[25]-1.185854122563141*fc[24]-0.3535533905932737*fc[22]-0.883883476483184*fc[21]+0.7071067811865475*fc[20]-0.9185586535436913*fc[19]-0.5303300858899104*fc[16]-0.6846531968814574*fc[15]+0.5477225575051661*fc[12]+1.060660171779821*fc[11]-0.3952847075210473*fc[9]+0.3162277660168379*fc[8]+0.7905694150420947*fc[7]+0.8215838362577489*fc[4]+0.4743416490252568*fc[2]+0.6123724356957944*fc[1]+0.3535533905932737*fc[0]; 
  } else {
    fUpOrdR[7] = (-0.7905694150420947*fr[26])+0.6123724356957944*fr[25]-1.185854122563141*fr[24]-0.3535533905932737*fr[22]-0.883883476483184*fr[21]+0.7071067811865475*fr[20]+0.9185586535436913*fr[19]-0.5303300858899104*fr[16]+0.6846531968814574*fr[15]-0.5477225575051661*fr[12]+1.060660171779821*fr[11]-0.3952847075210473*fr[9]+0.3162277660168379*fr[8]+0.7905694150420947*fr[7]-0.8215838362577489*fr[4]+0.4743416490252568*fr[2]-0.6123724356957944*fr[1]+0.3535533905932737*fr[0]; 
  }
  if (alphaR[4]+1.5*alphaR[1]+1.118033988749897*alphaR[0] > 0.) {
    fUpOrdR[8] = 0.6324555320336759*fc[26]+0.4898979485566357*fc[25]+0.9486832980505122*fc[24]+0.9486832980505137*fc[23]+0.2828427124746191*fc[22]+0.7071067811865475*(fc[21]+fc[20])+0.7348469228349525*fc[19]+0.7348469228349533*fc[18]+1.42302494707577*fc[17]+0.4242640687119281*fc[16]+0.5477225575051661*fc[15]+0.4242640687119285*fc[14]+1.060660171779821*fc[13]+0.5477225575051661*fc[12]+1.060660171779821*fc[11]+1.10227038425243*fc[10]+0.3162277660168379*(fc[9]+fc[8])+0.7905694150420947*fc[7]+0.6363961030678926*fc[6]+0.8215838362577489*(fc[5]+fc[4])+0.4743416490252568*(fc[3]+fc[2])+0.6123724356957944*fc[1]+0.3535533905932737*fc[0]; 
  } else {
    fUpOrdR[8] = 0.6324555320336759*fr[26]-0.4898979485566357*fr[25]+0.9486832980505122*fr[24]+0.9486832980505137*fr[23]+0.2828427124746191*fr[22]+0.7071067811865475*(fr[21]+fr[20])-0.7348469228349525*fr[19]-0.7348469228349533*fr[18]+1.42302494707577*fr[17]+0.4242640687119281*fr[16]-0.5477225575051661*fr[15]+0.4242640687119285*fr[14]+1.060660171779821*fr[13]-0.5477225575051661*fr[12]+1.060660171779821*fr[11]-1.10227038425243*fr[10]+0.3162277660168379*(fr[9]+fr[8])+0.7905694150420947*fr[7]+0.6363961030678926*fr[6]-0.8215838362577489*(fr[5]+fr[4])+0.4743416490252568*(fr[3]+fr[2])-0.6123724356957944*fr[1]+0.3535533905932737*fr[0]; 
  }

  double fUpR[9] = {0.};
  fUpR[0] = 0.154320987654321*fUpOrdR[8]+0.2469135802469136*fUpOrdR[7]+0.154320987654321*fUpOrdR[6]+0.2469135802469136*fUpOrdR[5]+0.3950617283950617*fUpOrdR[4]+0.2469135802469136*fUpOrdR[3]+0.154320987654321*fUpOrdR[2]+0.2469135802469136*fUpOrdR[1]+0.154320987654321*fUpOrdR[0]; 
  fUpR[1] = 0.2070433312499806*fUpOrdR[8]+0.3312693299999688*fUpOrdR[7]+0.2070433312499806*fUpOrdR[6]-0.2070433312499806*fUpOrdR[2]-0.3312693299999688*fUpOrdR[1]-0.2070433312499806*fUpOrdR[0]; 
  fUpR[2] = 0.2070433312499806*fUpOrdR[8]-0.2070433312499806*fUpOrdR[6]+0.3312693299999688*fUpOrdR[5]-0.3312693299999688*fUpOrdR[3]+0.2070433312499806*fUpOrdR[2]-0.2070433312499806*fUpOrdR[0]; 
  fUpR[3] = 0.2777777777777778*fUpOrdR[8]-0.2777777777777778*fUpOrdR[6]-0.2777777777777778*fUpOrdR[2]+0.2777777777777778*fUpOrdR[0]; 
  fUpR[4] = 0.138028887499987*fUpOrdR[8]+0.2208462199999792*fUpOrdR[7]+0.138028887499987*fUpOrdR[6]-0.2760577749999741*fUpOrdR[5]-0.4416924399999584*fUpOrdR[4]-0.2760577749999741*fUpOrdR[3]+0.138028887499987*fUpOrdR[2]+0.2208462199999792*fUpOrdR[1]+0.138028887499987*fUpOrdR[0]; 
  fUpR[5] = 0.138028887499987*fUpOrdR[8]-0.2760577749999741*fUpOrdR[7]+0.138028887499987*fUpOrdR[6]+0.2208462199999792*fUpOrdR[5]-0.4416924399999584*fUpOrdR[4]+0.2208462199999792*fUpOrdR[3]+0.138028887499987*fUpOrdR[2]-0.2760577749999741*fUpOrdR[1]+0.138028887499987*fUpOrdR[0]; 
  fUpR[6] = 0.1851851851851853*fUpOrdR[8]-0.1851851851851853*fUpOrdR[6]-0.3703703703703705*fUpOrdR[5]+0.3703703703703705*fUpOrdR[3]+0.1851851851851853*fUpOrdR[2]-0.1851851851851853*fUpOrdR[0]; 
  fUpR[7] = 0.1851851851851853*fUpOrdR[8]-0.3703703703703705*fUpOrdR[7]+0.1851851851851853*fUpOrdR[6]-0.1851851851851853*fUpOrdR[2]+0.3703703703703705*fUpOrdR[1]-0.1851851851851853*fUpOrdR[0]; 
  fUpR[8] = 0.1234567901234568*fUpOrdR[8]-0.2469135802469136*fUpOrdR[7]+0.1234567901234568*fUpOrdR[6]-0.2469135802469136*fUpOrdR[5]+0.4938271604938271*fUpOrdR[4]-0.2469135802469136*fUpOrdR[3]+0.1234567901234568*fUpOrdR[2]-0.2469135802469136*fUpOrdR[1]+0.1234567901234568*fUpOrdR[0]; 

  double GhatR[27] = {0.}; 
  GhatR[0] = 0.5*alphaR[4]*fUpR[4]+0.5*alphaR[1]*fUpR[1]+0.5*alphaR[0]*fUpR[0]; 
  GhatR[1] = 0.4472135954999579*alphaR[1]*fUpR[4]+0.4472135954999579*fUpR[1]*alphaR[4]+0.5*alphaR[0]*fUpR[1]+0.5*fUpR[0]*alphaR[1]; 
  GhatR[2] = 0.5000000000000001*alphaR[4]*fUpR[6]+0.5*alphaR[1]*fUpR[3]+0.5*alphaR[0]*fUpR[2]; 
  GhatR[3] = 0.447213595499958*alphaR[1]*fUpR[6]+0.4472135954999579*fUpR[3]*alphaR[4]+0.5*alphaR[0]*fUpR[3]+0.5*alphaR[1]*fUpR[2]; 
  GhatR[4] = 0.31943828249997*alphaR[4]*fUpR[4]+0.5*alphaR[0]*fUpR[4]+0.5*fUpR[0]*alphaR[4]+0.4472135954999579*alphaR[1]*fUpR[1]; 
  GhatR[5] = 0.5*alphaR[4]*fUpR[8]+0.5000000000000001*alphaR[1]*fUpR[7]+0.5*alphaR[0]*fUpR[5]; 
  GhatR[6] = 0.31943828249997*alphaR[4]*fUpR[6]+0.5*alphaR[0]*fUpR[6]+0.5000000000000001*fUpR[2]*alphaR[4]+0.447213595499958*alphaR[1]*fUpR[3]; 
  GhatR[7] = 0.447213595499958*alphaR[1]*fUpR[8]+0.4472135954999579*alphaR[4]*fUpR[7]+0.5*alphaR[0]*fUpR[7]+0.5000000000000001*alphaR[1]*fUpR[5]; 
  GhatR[8] = 0.31943828249997*alphaR[4]*fUpR[8]+0.5*alphaR[0]*fUpR[8]+0.447213595499958*alphaR[1]*fUpR[7]+0.5*alphaR[4]*fUpR[5]; 

  out[0] += (0.7071067811865475*GhatL[0]-0.7071067811865475*GhatR[0])*rdx2; 
  out[1] += ((-1.224744871391589*GhatR[0])-1.224744871391589*GhatL[0])*rdx2; 
  out[2] += (0.7071067811865475*GhatL[1]-0.7071067811865475*GhatR[1])*rdx2; 
  out[3] += (0.7071067811865475*GhatL[2]-0.7071067811865475*GhatR[2])*rdx2; 
  out[4] += ((-1.224744871391589*GhatR[1])-1.224744871391589*GhatL[1])*rdx2; 
  out[5] += ((-1.224744871391589*GhatR[2])-1.224744871391589*GhatL[2])*rdx2; 
  out[6] += (0.7071067811865475*GhatL[3]-0.7071067811865475*GhatR[3])*rdx2; 
  out[7] += (1.58113883008419*GhatL[0]-1.58113883008419*GhatR[0])*rdx2; 
  out[8] += (0.7071067811865475*GhatL[4]-0.7071067811865475*GhatR[4])*rdx2; 
  out[9] += (0.7071067811865475*GhatL[5]-0.7071067811865475*GhatR[5])*rdx2; 
  out[10] += ((-1.224744871391589*GhatR[3])-1.224744871391589*GhatL[3])*rdx2; 
  out[11] += (1.58113883008419*GhatL[1]-1.58113883008419*GhatR[1])*rdx2; 
  out[12] += ((-1.224744871391589*GhatR[4])-1.224744871391589*GhatL[4])*rdx2; 
  out[13] += (1.58113883008419*GhatL[2]-1.58113883008419*GhatR[2])*rdx2; 
  out[14] += (0.7071067811865475*GhatL[6]-0.7071067811865475*GhatR[6])*rdx2; 
  out[15] += ((-1.224744871391589*GhatR[5])-1.224744871391589*GhatL[5])*rdx2; 
  out[16] += (0.7071067811865475*GhatL[7]-0.7071067811865475*GhatR[7])*rdx2; 
  out[17] += (1.58113883008419*GhatL[3]-1.58113883008419*GhatR[3])*rdx2; 
  out[18] += ((-1.224744871391589*GhatR[6])-1.224744871391589*GhatL[6])*rdx2; 
  out[19] += ((-1.224744871391589*GhatR[7])-1.224744871391589*GhatL[7])*rdx2; 
  out[20] += (1.58113883008419*GhatL[4]-1.58113883008419*GhatR[4])*rdx2; 
  out[21] += (1.58113883008419*GhatL[5]-1.58113883008419*GhatR[5])*rdx2; 
  out[22] += (0.7071067811865475*GhatL[8]-0.7071067811865475*GhatR[8])*rdx2; 
  out[23] += (1.58113883008419*GhatL[6]-1.58113883008419*GhatR[6])*rdx2; 
  out[24] += (1.58113883008419*GhatL[7]-1.58113883008419*GhatR[7])*rdx2; 
  out[25] += ((-1.224744871391589*GhatR[8])-1.224744871391589*GhatL[8])*rdx2; 
  out[26] += (1.58113883008419*GhatL[8]-1.58113883008419*GhatR[8])*rdx2; 

} 