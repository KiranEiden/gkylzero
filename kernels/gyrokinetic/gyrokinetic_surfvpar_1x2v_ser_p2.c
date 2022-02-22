#include <gkyl_gyrokinetic_kernels.h>
GKYL_CU_DH void gyrokinetic_surfx_1x2v_ser_2(const double *w, const double *dxv, const double q_, const double m_, const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi, const double *apar, const double *apardot, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
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

  double hamil[20]; 
  hamil[0] = (0.2357022603955158*(3.0*rdvpar2Sq*(2.0*m_*wvparSq+2.828427124746191*(bmag[0]*wmu+phi[0]*q_))+2.0*m_))/rdvpar2Sq; 
  hamil[1] = 2.0*(bmag[1]*wmu+phi[1]*q_); 
  hamil[2] = (1.632993161855453*m_*wvpar)/rdvpar2; 
  hamil[3] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[5] = (1.154700538379252*bmag[1])/rdmu2; 
  hamil[7] = 2.0*(bmag[2]*wmu+phi[2]*q_); 
  hamil[8] = (0.421637021355784*m_)/rdvpar2Sq; 
  hamil[13] = (1.154700538379251*bmag[2])/rdmu2; 

  double BstarZdBmag[20]; 
  BstarZdBmag[0] = (1.414213562373095*(1.732050807568877*(2.23606797749979*jacobtot_inv[1]*b_y[2]+jacobtot_inv[0]*b_y[1])*m_*rdx2*wvpar+(cmag[2]*jacobtot_inv[2]+cmag[1]*jacobtot_inv[1]+cmag[0]*jacobtot_inv[0])*q_))/q_; 
  BstarZdBmag[1] = (0.2828427124746191*(5.0*(1.732050807568877*(2.0*b_y[2]*jacobtot_inv[2]+b_y[1]*jacobtot_inv[1])*m_*rdx2*wvpar+(cmag[0]*jacobtot_inv[1]+jacobtot_inv[0]*cmag[1])*q_)+2.23606797749979*(8.660254037844386*jacobtot_inv[0]*b_y[2]*m_*rdx2*wvpar+2.0*(cmag[1]*jacobtot_inv[2]+jacobtot_inv[1]*cmag[2])*q_)))/q_; 
  BstarZdBmag[2] = (1.414213562373095*(2.23606797749979*jacobtot_inv[1]*b_y[2]+jacobtot_inv[0]*b_y[1])*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[4] = (1.414213562373095*(b_y[2]*(2.0*jacobtot_inv[2]+2.23606797749979*jacobtot_inv[0])+b_y[1]*jacobtot_inv[1])*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[7] = (0.04040610178208843*(60.6217782649107*(b_y[1]*jacobtot_inv[2]+2.0*jacobtot_inv[1]*b_y[2])*m_*rdx2*wvpar+(4.47213595499958*(5.0*cmag[2]*jacobtot_inv[2]+7.0*cmag[1]*jacobtot_inv[1])+35.0*(cmag[0]*jacobtot_inv[2]+jacobtot_inv[0]*cmag[2]))*q_))/q_; 
  BstarZdBmag[11] = (1.414213562373095*(b_y[1]*jacobtot_inv[2]+2.0*jacobtot_inv[1]*b_y[2])*m_*rdx2)/(q_*rdvpar2); 

  double alphaL[8]; 
  alphaL[0] = -(0.25*((2.23606797749979*(1.732050807568877*BstarZdBmag[1]-3.0*BstarZdBmag[4])*hamil[7]+hamil[1]*(1.732050807568877*BstarZdBmag[0]-3.0*BstarZdBmag[2]))*rdx2+5.656854249492382*apardot[0]*q_))/m_; 
  alphaL[1] = (0.05*((30.0*hamil[7]*BstarZdBmag[11]+1.732050807568877*((-10.0*BstarZdBmag[7]*hamil[7])-5.0*BstarZdBmag[1]*hamil[1])+2.23606797749979*(15.0*BstarZdBmag[2]-8.660254037844386*BstarZdBmag[0])*hamil[7]+15.0*hamil[1]*BstarZdBmag[4])*rdx2-28.28427124746191*apardot[1]*q_))/m_; 
  alphaL[2] = (0.25*(3.872983346207417*(1.732050807568877*BstarZdBmag[4]-1.0*BstarZdBmag[1])*hamil[13]+(3.0*BstarZdBmag[2]-1.732050807568877*BstarZdBmag[0])*hamil[5])*rdx2)/m_; 
  alphaL[3] = -(0.05*((3.872983346207417*(4.47213595499958*BstarZdBmag[7]-8.660254037844386*BstarZdBmag[2]+5.0*BstarZdBmag[0])-30.0*BstarZdBmag[11])*hamil[13]+(8.660254037844386*BstarZdBmag[1]-15.0*BstarZdBmag[4])*hamil[5])*rdx2)/m_; 
  alphaL[4] = (0.05*((15.0*hamil[1]*BstarZdBmag[11]+1.732050807568877*((-10.0*BstarZdBmag[1]*hamil[7])-5.0*hamil[1]*BstarZdBmag[7])+30.0*BstarZdBmag[4]*hamil[7])*rdx2-28.28427124746191*apardot[2]*q_))/m_; 
  alphaL[6] = -(0.05*((17.32050807568877*BstarZdBmag[1]-30.0*BstarZdBmag[4])*hamil[13]+hamil[5]*(8.660254037844387*BstarZdBmag[7]-15.0*BstarZdBmag[11]))*rdx2)/m_; 

  double alphaR[8]; 
  alphaR[0] = -(0.25*((2.23606797749979*(3.0*BstarZdBmag[4]+1.732050807568877*BstarZdBmag[1])*hamil[7]+hamil[1]*(3.0*BstarZdBmag[2]+1.732050807568877*BstarZdBmag[0]))*rdx2+5.656854249492382*apardot[0]*q_))/m_; 
  alphaR[1] = -(0.05*((30.0*hamil[7]*BstarZdBmag[11]+1.732050807568877*(10.0*BstarZdBmag[7]*hamil[7]+5.0*BstarZdBmag[1]*hamil[1])+2.23606797749979*(15.0*BstarZdBmag[2]+8.660254037844386*BstarZdBmag[0])*hamil[7]+15.0*hamil[1]*BstarZdBmag[4])*rdx2+28.28427124746191*apardot[1]*q_))/m_; 
  alphaR[2] = -(0.25*(3.872983346207417*(1.732050807568877*BstarZdBmag[4]+BstarZdBmag[1])*hamil[13]+(3.0*BstarZdBmag[2]+1.732050807568877*BstarZdBmag[0])*hamil[5])*rdx2)/m_; 
  alphaR[3] = -(0.05*((30.0*BstarZdBmag[11]+3.872983346207417*(4.47213595499958*BstarZdBmag[7]+8.660254037844386*BstarZdBmag[2]+5.0*BstarZdBmag[0]))*hamil[13]+(15.0*BstarZdBmag[4]+8.660254037844386*BstarZdBmag[1])*hamil[5])*rdx2)/m_; 
  alphaR[4] = -(0.05*((15.0*hamil[1]*BstarZdBmag[11]+1.732050807568877*(10.0*BstarZdBmag[1]*hamil[7]+5.0*hamil[1]*BstarZdBmag[7])+30.0*BstarZdBmag[4]*hamil[7])*rdx2+28.28427124746191*apardot[2]*q_))/m_; 
  alphaR[6] = -(0.05*((30.0*BstarZdBmag[4]+17.32050807568877*BstarZdBmag[1])*hamil[13]+hamil[5]*(15.0*BstarZdBmag[11]+8.660254037844387*BstarZdBmag[7]))*rdx2)/m_; 

  double fUpOrdL[9];
  if (alphaL[0]-1.118033988749893*alphaL[4] > 0.) {
    fUpOrdL[0] = (-0.6846531968814574*(fl[16]+fl[11]))-0.3952847075210473*fl[9]+0.7905694150420947*fl[8]-0.3952847075210473*fl[7]+0.6123724356957944*fl[2]+0.3535533905932737*fl[0]; 
  } else {
    fUpOrdL[0] = 0.6846531968814574*(fc[16]+fc[11])-0.3952847075210473*fc[9]+0.7905694150420947*fc[8]-0.3952847075210473*fc[7]-0.6123724356957944*fc[2]+0.3535533905932737*fc[0]; 
  }
  if (alphaL[6]-0.7453559924999286*alphaL[4]-0.8944271909999143*alphaL[2]+0.6666666666666666*alphaL[0] > 0.) {
    fUpOrdL[1] = 0.9185586535436915*fl[17]+0.547722557505166*fl[16]-1.060660171779821*fl[14]+0.5303300858899105*fl[13]-0.6846531968814574*fl[11]+0.3162277660168378*fl[9]+0.7905694150420947*fl[8]-0.3952847075210473*fl[7]-0.821583836257749*fl[6]-0.4743416490252567*fl[3]+0.6123724356957944*fl[2]+0.3535533905932737*fl[0]; 
  } else {
    fUpOrdL[1] = (-0.9185586535436915*fc[17])-0.547722557505166*fc[16]-1.060660171779821*fc[14]+0.5303300858899105*fc[13]+0.6846531968814574*fc[11]+0.3162277660168378*fc[9]+0.7905694150420947*fc[8]-0.3952847075210473*fc[7]+0.821583836257749*fc[6]-0.4743416490252567*fc[3]-0.6123724356957944*fc[2]+0.3535533905932737*fc[0]; 
  }
  if ((-1.0*alphaL[6])-0.7453559924999286*alphaL[4]+0.8944271909999143*alphaL[2]+0.6666666666666666*alphaL[0] > 0.) {
    fUpOrdL[2] = (-0.9185586535436915*fl[17])+0.547722557505166*fl[16]+1.060660171779821*fl[14]-0.5303300858899105*fl[13]-0.6846531968814574*fl[11]+0.3162277660168378*fl[9]+0.7905694150420947*fl[8]-0.3952847075210473*fl[7]+0.821583836257749*fl[6]+0.4743416490252567*fl[3]+0.6123724356957944*fl[2]+0.3535533905932737*fl[0]; 
  } else {
    fUpOrdL[2] = 0.9185586535436915*fc[17]-0.547722557505166*fc[16]+1.060660171779821*fc[14]-0.5303300858899105*fc[13]+0.6846531968814574*fc[11]+0.3162277660168378*fc[9]+0.7905694150420947*fc[8]-0.3952847075210473*fc[7]-0.821583836257749*fc[6]+0.4743416490252567*fc[3]-0.6123724356957944*fc[2]+0.3535533905932737*fc[0]; 
  }
  if (alphaL[4]-1.5*alphaL[1]+1.118033988749897*alphaL[0] > 0.) {
    fUpOrdL[3] = 0.9185586535436915*fl[19]-0.6846531968814574*fl[16]+0.5303300858899105*fl[15]-1.060660171779821*fl[12]+0.547722557505166*fl[11]-0.3952847075210473*fl[9]+0.7905694150420947*fl[8]+0.3162277660168378*fl[7]-0.821583836257749*fl[4]+0.6123724356957944*fl[2]-0.4743416490252567*fl[1]+0.3535533905932737*fl[0]; 
  } else {
    fUpOrdL[3] = (-0.9185586535436915*fc[19])+0.6846531968814574*fc[16]+0.5303300858899105*fc[15]-1.060660171779821*fc[12]-0.547722557505166*fc[11]-0.3952847075210473*fc[9]+0.7905694150420947*fc[8]+0.3162277660168378*fc[7]+0.821583836257749*fc[4]-0.6123724356957944*fc[2]-0.4743416490252567*fc[1]+0.3535533905932737*fc[0]; 
  }
  if ((-1.0*alphaL[6])+0.7453559924999286*alphaL[4]+1.5*alphaL[3]-1.118033988749893*alphaL[2]-1.118033988749893*alphaL[1]+0.8333333333333334*alphaL[0] > 0.) {
    fUpOrdL[4] = (-0.7348469228349532*fl[19])+1.423024947075771*fl[18]-0.7348469228349532*fl[17]+0.547722557505166*fl[16]-0.4242640687119284*fl[15]-1.060660171779821*fl[14]-0.4242640687119284*fl[13]-1.060660171779821*fl[12]+0.547722557505166*fl[11]+1.10227038425243*fl[10]+0.3162277660168378*fl[9]+0.7905694150420947*fl[8]+0.3162277660168378*fl[7]-0.821583836257749*fl[6]+0.6363961030678926*fl[5]-0.821583836257749*fl[4]-0.4743416490252567*fl[3]+0.6123724356957944*fl[2]-0.4743416490252567*fl[1]+0.3535533905932737*fl[0]; 
  } else {
    fUpOrdL[4] = 0.7348469228349532*fc[19]+1.423024947075771*fc[18]+0.7348469228349532*fc[17]-0.547722557505166*fc[16]-0.4242640687119284*fc[15]-1.060660171779821*fc[14]-0.4242640687119284*fc[13]-1.060660171779821*fc[12]-0.547722557505166*fc[11]-1.10227038425243*fc[10]+0.3162277660168378*fc[9]+0.7905694150420947*fc[8]+0.3162277660168378*fc[7]+0.821583836257749*fc[6]+0.6363961030678926*fc[5]+0.821583836257749*fc[4]-0.4743416490252567*fc[3]-0.6123724356957944*fc[2]-0.4743416490252567*fc[1]+0.3535533905932737*fc[0]; 
  }
  if (alphaL[6]+0.7453559924999286*alphaL[4]-1.5*alphaL[3]+1.118033988749893*alphaL[2]-1.118033988749893*alphaL[1]+0.8333333333333334*alphaL[0] > 0.) {
    fUpOrdL[5] = (-0.7348469228349532*fl[19])-1.423024947075771*fl[18]+0.7348469228349532*fl[17]+0.547722557505166*fl[16]-0.4242640687119284*fl[15]+1.060660171779821*fl[14]+0.4242640687119284*fl[13]-1.060660171779821*fl[12]+0.547722557505166*fl[11]-1.10227038425243*fl[10]+0.3162277660168378*fl[9]+0.7905694150420947*fl[8]+0.3162277660168378*fl[7]+0.821583836257749*fl[6]-0.6363961030678926*fl[5]-0.821583836257749*fl[4]+0.4743416490252567*fl[3]+0.6123724356957944*fl[2]-0.4743416490252567*fl[1]+0.3535533905932737*fl[0]; 
  } else {
    fUpOrdL[5] = 0.7348469228349532*fc[19]-1.423024947075771*fc[18]-0.7348469228349532*fc[17]-0.547722557505166*fc[16]-0.4242640687119284*fc[15]+1.060660171779821*fc[14]+0.4242640687119284*fc[13]-1.060660171779821*fc[12]-0.547722557505166*fc[11]+1.10227038425243*fc[10]+0.3162277660168378*fc[9]+0.7905694150420947*fc[8]+0.3162277660168378*fc[7]-0.821583836257749*fc[6]-0.6363961030678926*fc[5]+0.821583836257749*fc[4]+0.4743416490252567*fc[3]-0.6123724356957944*fc[2]-0.4743416490252567*fc[1]+0.3535533905932737*fc[0]; 
  }
  if (alphaL[4]+1.5*alphaL[1]+1.118033988749897*alphaL[0] > 0.) {
    fUpOrdL[6] = (-0.9185586535436915*fl[19])-0.6846531968814574*fl[16]-0.5303300858899105*fl[15]+1.060660171779821*fl[12]+0.547722557505166*fl[11]-0.3952847075210473*fl[9]+0.7905694150420947*fl[8]+0.3162277660168378*fl[7]+0.821583836257749*fl[4]+0.6123724356957944*fl[2]+0.4743416490252567*fl[1]+0.3535533905932737*fl[0]; 
  } else {
    fUpOrdL[6] = 0.9185586535436915*fc[19]+0.6846531968814574*fc[16]-0.5303300858899105*fc[15]+1.060660171779821*fc[12]-0.547722557505166*fc[11]-0.3952847075210473*fc[9]+0.7905694150420947*fc[8]+0.3162277660168378*fc[7]-0.821583836257749*fc[4]-0.6123724356957944*fc[2]+0.4743416490252567*fc[1]+0.3535533905932737*fc[0]; 
  }
  if ((-1.0*alphaL[6])+0.7453559924999286*alphaL[4]-1.5*alphaL[3]-1.118033988749893*alphaL[2]+1.118033988749893*alphaL[1]+0.8333333333333334*alphaL[0] > 0.) {
    fUpOrdL[7] = 0.7348469228349532*fl[19]-1.423024947075771*fl[18]-0.7348469228349532*fl[17]+0.547722557505166*fl[16]+0.4242640687119284*fl[15]-1.060660171779821*fl[14]-0.4242640687119284*fl[13]+1.060660171779821*fl[12]+0.547722557505166*fl[11]-1.10227038425243*fl[10]+0.3162277660168378*fl[9]+0.7905694150420947*fl[8]+0.3162277660168378*fl[7]-0.821583836257749*fl[6]-0.6363961030678926*fl[5]+0.821583836257749*fl[4]-0.4743416490252567*fl[3]+0.6123724356957944*fl[2]+0.4743416490252567*fl[1]+0.3535533905932737*fl[0]; 
  } else {
    fUpOrdL[7] = (-0.7348469228349532*fc[19])-1.423024947075771*fc[18]+0.7348469228349532*fc[17]-0.547722557505166*fc[16]+0.4242640687119284*fc[15]-1.060660171779821*fc[14]-0.4242640687119284*fc[13]+1.060660171779821*fc[12]-0.547722557505166*fc[11]+1.10227038425243*fc[10]+0.3162277660168378*fc[9]+0.7905694150420947*fc[8]+0.3162277660168378*fc[7]+0.821583836257749*fc[6]-0.6363961030678926*fc[5]-0.821583836257749*fc[4]-0.4743416490252567*fc[3]-0.6123724356957944*fc[2]+0.4743416490252567*fc[1]+0.3535533905932737*fc[0]; 
  }
  if (alphaL[6]+0.7453559924999286*alphaL[4]+1.5*alphaL[3]+1.118033988749893*alphaL[2]+1.118033988749893*alphaL[1]+0.8333333333333334*alphaL[0] > 0.) {
    fUpOrdL[8] = 0.7348469228349532*fl[19]+1.423024947075771*fl[18]+0.7348469228349532*fl[17]+0.547722557505166*fl[16]+0.4242640687119284*fl[15]+1.060660171779821*fl[14]+0.4242640687119284*fl[13]+1.060660171779821*fl[12]+0.547722557505166*fl[11]+1.10227038425243*fl[10]+0.3162277660168378*fl[9]+0.7905694150420947*fl[8]+0.3162277660168378*fl[7]+0.821583836257749*fl[6]+0.6363961030678926*fl[5]+0.821583836257749*fl[4]+0.4743416490252567*fl[3]+0.6123724356957944*fl[2]+0.4743416490252567*fl[1]+0.3535533905932737*fl[0]; 
  } else {
    fUpOrdL[8] = (-0.7348469228349532*fc[19])+1.423024947075771*fc[18]-0.7348469228349532*fc[17]-0.547722557505166*fc[16]+0.4242640687119284*fc[15]+1.060660171779821*fc[14]+0.4242640687119284*fc[13]+1.060660171779821*fc[12]-0.547722557505166*fc[11]-1.10227038425243*fc[10]+0.3162277660168378*fc[9]+0.7905694150420947*fc[8]+0.3162277660168378*fc[7]-0.821583836257749*fc[6]+0.6363961030678926*fc[5]-0.821583836257749*fc[4]+0.4743416490252567*fc[3]-0.6123724356957944*fc[2]+0.4743416490252567*fc[1]+0.3535533905932737*fc[0]; 
  }

  double fUpL[8] = {0.};
  fUpL[0] = 0.154320987654321*fUpOrdL[8]+0.154320987654321*fUpOrdL[7]+0.2469135802469136*fUpOrdL[6]+0.154320987654321*fUpOrdL[5]+0.154320987654321*fUpOrdL[4]+0.2469135802469136*fUpOrdL[3]+0.2469135802469136*fUpOrdL[2]+0.2469135802469136*fUpOrdL[1]+0.3950617283950617*fUpOrdL[0]; 
  fUpL[1] = 0.2070433312499804*fUpOrdL[8]+0.2070433312499804*fUpOrdL[7]+0.3312693299999682*fUpOrdL[6]-0.2070433312499804*fUpOrdL[5]-0.2070433312499804*fUpOrdL[4]-0.3312693299999682*fUpOrdL[3]; 
  fUpL[2] = 0.2070433312499804*fUpOrdL[8]-0.2070433312499804*fUpOrdL[7]+0.2070433312499804*fUpOrdL[5]-0.2070433312499804*fUpOrdL[4]+0.3312693299999682*fUpOrdL[2]-0.3312693299999682*fUpOrdL[1]; 
  fUpL[3] = 0.2777777777777778*fUpOrdL[8]-0.2777777777777778*fUpOrdL[7]-0.2777777777777778*fUpOrdL[5]+0.2777777777777778*fUpOrdL[4]; 
  fUpL[4] = 0.138028887499987*fUpOrdL[8]+0.138028887499987*fUpOrdL[7]+0.2208462199999792*fUpOrdL[6]+0.138028887499987*fUpOrdL[5]+0.138028887499987*fUpOrdL[4]+0.2208462199999792*fUpOrdL[3]-0.2760577749999741*fUpOrdL[2]-0.2760577749999741*fUpOrdL[1]-0.4416924399999584*fUpOrdL[0]; 
  fUpL[5] = 0.138028887499987*fUpOrdL[8]+0.138028887499987*fUpOrdL[7]-0.2760577749999741*fUpOrdL[6]+0.138028887499987*fUpOrdL[5]+0.138028887499987*fUpOrdL[4]-0.2760577749999741*fUpOrdL[3]+0.2208462199999792*fUpOrdL[2]+0.2208462199999792*fUpOrdL[1]-0.4416924399999584*fUpOrdL[0]; 
  fUpL[6] = 0.1851851851851852*fUpOrdL[8]-0.1851851851851852*fUpOrdL[7]+0.1851851851851852*fUpOrdL[5]-0.1851851851851852*fUpOrdL[4]-0.3703703703703705*fUpOrdL[2]+0.3703703703703705*fUpOrdL[1]; 
  fUpL[7] = 0.1851851851851852*fUpOrdL[8]+0.1851851851851852*fUpOrdL[7]-0.3703703703703705*fUpOrdL[6]-0.1851851851851852*fUpOrdL[5]-0.1851851851851852*fUpOrdL[4]+0.3703703703703705*fUpOrdL[3]; 

  double GhatL[20] = {0.}; 
  GhatL[0] += 0.5*alphaL[6]*fUpL[6]+0.5*alphaL[4]*fUpL[4]+0.5*alphaL[3]*fUpL[3]+0.5*alphaL[2]*fUpL[2]+0.5*alphaL[1]*fUpL[1]+0.5*alphaL[0]*fUpL[0]; 
  GhatL[1] += 0.447213595499958*alphaL[3]*fUpL[6]+0.447213595499958*fUpL[3]*alphaL[6]+0.4472135954999579*alphaL[1]*fUpL[4]+0.4472135954999579*fUpL[1]*alphaL[4]+0.5*alphaL[2]*fUpL[3]+0.5*fUpL[2]*alphaL[3]+0.5*alphaL[0]*fUpL[1]+0.5*fUpL[0]*alphaL[1]; 
  GhatL[2] += 0.447213595499958*alphaL[3]*fUpL[7]+0.5000000000000001*alphaL[4]*fUpL[6]+0.5000000000000001*fUpL[4]*alphaL[6]+0.4472135954999579*alphaL[2]*fUpL[5]+0.5*alphaL[1]*fUpL[3]+0.5*fUpL[1]*alphaL[3]+0.5*alphaL[0]*fUpL[2]+0.5*fUpL[0]*alphaL[2]; 
  GhatL[3] += 0.4*alphaL[6]*fUpL[7]+0.447213595499958*alphaL[2]*fUpL[7]+0.447213595499958*alphaL[1]*fUpL[6]+0.447213595499958*fUpL[1]*alphaL[6]+0.4472135954999579*alphaL[3]*fUpL[5]+0.4472135954999579*alphaL[3]*fUpL[4]+0.4472135954999579*fUpL[3]*alphaL[4]+0.5*alphaL[0]*fUpL[3]+0.5*fUpL[0]*alphaL[3]+0.5*alphaL[1]*fUpL[2]+0.5*fUpL[1]*alphaL[2]; 
  GhatL[4] += 0.31943828249997*alphaL[6]*fUpL[6]+0.5000000000000001*alphaL[2]*fUpL[6]+0.5000000000000001*fUpL[2]*alphaL[6]+0.31943828249997*alphaL[4]*fUpL[4]+0.5*alphaL[0]*fUpL[4]+0.5*fUpL[0]*alphaL[4]+0.4472135954999579*alphaL[3]*fUpL[3]+0.4472135954999579*alphaL[1]*fUpL[1]; 
  GhatL[5] += 0.5000000000000001*alphaL[1]*fUpL[7]+0.4472135954999579*alphaL[6]*fUpL[6]+0.5*alphaL[0]*fUpL[5]+0.4472135954999579*alphaL[3]*fUpL[3]+0.4472135954999579*alphaL[2]*fUpL[2]; 
  GhatL[6] += 0.4*alphaL[3]*fUpL[7]+0.31943828249997*alphaL[4]*fUpL[6]+0.5*alphaL[0]*fUpL[6]+0.4472135954999579*fUpL[5]*alphaL[6]+0.31943828249997*fUpL[4]*alphaL[6]+0.5*fUpL[0]*alphaL[6]+0.5000000000000001*alphaL[2]*fUpL[4]+0.5000000000000001*fUpL[2]*alphaL[4]+0.447213595499958*alphaL[1]*fUpL[3]+0.447213595499958*fUpL[1]*alphaL[3]; 
  GhatL[7] += 0.4472135954999579*alphaL[4]*fUpL[7]+0.5*alphaL[0]*fUpL[7]+0.4*alphaL[3]*fUpL[6]+0.4*fUpL[3]*alphaL[6]+0.5000000000000001*alphaL[1]*fUpL[5]+0.447213595499958*alphaL[2]*fUpL[3]+0.447213595499958*fUpL[2]*alphaL[3]; 

  double fUpOrdR[9];
  if (alphaR[0]-1.118033988749893*alphaR[4] > 0.) {
    fUpOrdR[0] = (-0.6846531968814574*(fc[16]+fc[11]))-0.3952847075210473*fc[9]+0.7905694150420947*fc[8]-0.3952847075210473*fc[7]+0.6123724356957944*fc[2]+0.3535533905932737*fc[0]; 
  } else {
    fUpOrdR[0] = 0.6846531968814574*(fr[16]+fr[11])-0.3952847075210473*fr[9]+0.7905694150420947*fr[8]-0.3952847075210473*fr[7]-0.6123724356957944*fr[2]+0.3535533905932737*fr[0]; 
  }
  if (alphaR[6]-0.7453559924999286*alphaR[4]-0.8944271909999143*alphaR[2]+0.6666666666666666*alphaR[0] > 0.) {
    fUpOrdR[1] = 0.9185586535436915*fc[17]+0.547722557505166*fc[16]-1.060660171779821*fc[14]+0.5303300858899105*fc[13]-0.6846531968814574*fc[11]+0.3162277660168378*fc[9]+0.7905694150420947*fc[8]-0.3952847075210473*fc[7]-0.821583836257749*fc[6]-0.4743416490252567*fc[3]+0.6123724356957944*fc[2]+0.3535533905932737*fc[0]; 
  } else {
    fUpOrdR[1] = (-0.9185586535436915*fr[17])-0.547722557505166*fr[16]-1.060660171779821*fr[14]+0.5303300858899105*fr[13]+0.6846531968814574*fr[11]+0.3162277660168378*fr[9]+0.7905694150420947*fr[8]-0.3952847075210473*fr[7]+0.821583836257749*fr[6]-0.4743416490252567*fr[3]-0.6123724356957944*fr[2]+0.3535533905932737*fr[0]; 
  }
  if ((-1.0*alphaR[6])-0.7453559924999286*alphaR[4]+0.8944271909999143*alphaR[2]+0.6666666666666666*alphaR[0] > 0.) {
    fUpOrdR[2] = (-0.9185586535436915*fc[17])+0.547722557505166*fc[16]+1.060660171779821*fc[14]-0.5303300858899105*fc[13]-0.6846531968814574*fc[11]+0.3162277660168378*fc[9]+0.7905694150420947*fc[8]-0.3952847075210473*fc[7]+0.821583836257749*fc[6]+0.4743416490252567*fc[3]+0.6123724356957944*fc[2]+0.3535533905932737*fc[0]; 
  } else {
    fUpOrdR[2] = 0.9185586535436915*fr[17]-0.547722557505166*fr[16]+1.060660171779821*fr[14]-0.5303300858899105*fr[13]+0.6846531968814574*fr[11]+0.3162277660168378*fr[9]+0.7905694150420947*fr[8]-0.3952847075210473*fr[7]-0.821583836257749*fr[6]+0.4743416490252567*fr[3]-0.6123724356957944*fr[2]+0.3535533905932737*fr[0]; 
  }
  if (alphaR[4]-1.5*alphaR[1]+1.118033988749897*alphaR[0] > 0.) {
    fUpOrdR[3] = 0.9185586535436915*fc[19]-0.6846531968814574*fc[16]+0.5303300858899105*fc[15]-1.060660171779821*fc[12]+0.547722557505166*fc[11]-0.3952847075210473*fc[9]+0.7905694150420947*fc[8]+0.3162277660168378*fc[7]-0.821583836257749*fc[4]+0.6123724356957944*fc[2]-0.4743416490252567*fc[1]+0.3535533905932737*fc[0]; 
  } else {
    fUpOrdR[3] = (-0.9185586535436915*fr[19])+0.6846531968814574*fr[16]+0.5303300858899105*fr[15]-1.060660171779821*fr[12]-0.547722557505166*fr[11]-0.3952847075210473*fr[9]+0.7905694150420947*fr[8]+0.3162277660168378*fr[7]+0.821583836257749*fr[4]-0.6123724356957944*fr[2]-0.4743416490252567*fr[1]+0.3535533905932737*fr[0]; 
  }
  if ((-1.0*alphaR[6])+0.7453559924999286*alphaR[4]+1.5*alphaR[3]-1.118033988749893*alphaR[2]-1.118033988749893*alphaR[1]+0.8333333333333334*alphaR[0] > 0.) {
    fUpOrdR[4] = (-0.7348469228349532*fc[19])+1.423024947075771*fc[18]-0.7348469228349532*fc[17]+0.547722557505166*fc[16]-0.4242640687119284*fc[15]-1.060660171779821*fc[14]-0.4242640687119284*fc[13]-1.060660171779821*fc[12]+0.547722557505166*fc[11]+1.10227038425243*fc[10]+0.3162277660168378*fc[9]+0.7905694150420947*fc[8]+0.3162277660168378*fc[7]-0.821583836257749*fc[6]+0.6363961030678926*fc[5]-0.821583836257749*fc[4]-0.4743416490252567*fc[3]+0.6123724356957944*fc[2]-0.4743416490252567*fc[1]+0.3535533905932737*fc[0]; 
  } else {
    fUpOrdR[4] = 0.7348469228349532*fr[19]+1.423024947075771*fr[18]+0.7348469228349532*fr[17]-0.547722557505166*fr[16]-0.4242640687119284*fr[15]-1.060660171779821*fr[14]-0.4242640687119284*fr[13]-1.060660171779821*fr[12]-0.547722557505166*fr[11]-1.10227038425243*fr[10]+0.3162277660168378*fr[9]+0.7905694150420947*fr[8]+0.3162277660168378*fr[7]+0.821583836257749*fr[6]+0.6363961030678926*fr[5]+0.821583836257749*fr[4]-0.4743416490252567*fr[3]-0.6123724356957944*fr[2]-0.4743416490252567*fr[1]+0.3535533905932737*fr[0]; 
  }
  if (alphaR[6]+0.7453559924999286*alphaR[4]-1.5*alphaR[3]+1.118033988749893*alphaR[2]-1.118033988749893*alphaR[1]+0.8333333333333334*alphaR[0] > 0.) {
    fUpOrdR[5] = (-0.7348469228349532*fc[19])-1.423024947075771*fc[18]+0.7348469228349532*fc[17]+0.547722557505166*fc[16]-0.4242640687119284*fc[15]+1.060660171779821*fc[14]+0.4242640687119284*fc[13]-1.060660171779821*fc[12]+0.547722557505166*fc[11]-1.10227038425243*fc[10]+0.3162277660168378*fc[9]+0.7905694150420947*fc[8]+0.3162277660168378*fc[7]+0.821583836257749*fc[6]-0.6363961030678926*fc[5]-0.821583836257749*fc[4]+0.4743416490252567*fc[3]+0.6123724356957944*fc[2]-0.4743416490252567*fc[1]+0.3535533905932737*fc[0]; 
  } else {
    fUpOrdR[5] = 0.7348469228349532*fr[19]-1.423024947075771*fr[18]-0.7348469228349532*fr[17]-0.547722557505166*fr[16]-0.4242640687119284*fr[15]+1.060660171779821*fr[14]+0.4242640687119284*fr[13]-1.060660171779821*fr[12]-0.547722557505166*fr[11]+1.10227038425243*fr[10]+0.3162277660168378*fr[9]+0.7905694150420947*fr[8]+0.3162277660168378*fr[7]-0.821583836257749*fr[6]-0.6363961030678926*fr[5]+0.821583836257749*fr[4]+0.4743416490252567*fr[3]-0.6123724356957944*fr[2]-0.4743416490252567*fr[1]+0.3535533905932737*fr[0]; 
  }
  if (alphaR[4]+1.5*alphaR[1]+1.118033988749897*alphaR[0] > 0.) {
    fUpOrdR[6] = (-0.9185586535436915*fc[19])-0.6846531968814574*fc[16]-0.5303300858899105*fc[15]+1.060660171779821*fc[12]+0.547722557505166*fc[11]-0.3952847075210473*fc[9]+0.7905694150420947*fc[8]+0.3162277660168378*fc[7]+0.821583836257749*fc[4]+0.6123724356957944*fc[2]+0.4743416490252567*fc[1]+0.3535533905932737*fc[0]; 
  } else {
    fUpOrdR[6] = 0.9185586535436915*fr[19]+0.6846531968814574*fr[16]-0.5303300858899105*fr[15]+1.060660171779821*fr[12]-0.547722557505166*fr[11]-0.3952847075210473*fr[9]+0.7905694150420947*fr[8]+0.3162277660168378*fr[7]-0.821583836257749*fr[4]-0.6123724356957944*fr[2]+0.4743416490252567*fr[1]+0.3535533905932737*fr[0]; 
  }
  if ((-1.0*alphaR[6])+0.7453559924999286*alphaR[4]-1.5*alphaR[3]-1.118033988749893*alphaR[2]+1.118033988749893*alphaR[1]+0.8333333333333334*alphaR[0] > 0.) {
    fUpOrdR[7] = 0.7348469228349532*fc[19]-1.423024947075771*fc[18]-0.7348469228349532*fc[17]+0.547722557505166*fc[16]+0.4242640687119284*fc[15]-1.060660171779821*fc[14]-0.4242640687119284*fc[13]+1.060660171779821*fc[12]+0.547722557505166*fc[11]-1.10227038425243*fc[10]+0.3162277660168378*fc[9]+0.7905694150420947*fc[8]+0.3162277660168378*fc[7]-0.821583836257749*fc[6]-0.6363961030678926*fc[5]+0.821583836257749*fc[4]-0.4743416490252567*fc[3]+0.6123724356957944*fc[2]+0.4743416490252567*fc[1]+0.3535533905932737*fc[0]; 
  } else {
    fUpOrdR[7] = (-0.7348469228349532*fr[19])-1.423024947075771*fr[18]+0.7348469228349532*fr[17]-0.547722557505166*fr[16]+0.4242640687119284*fr[15]-1.060660171779821*fr[14]-0.4242640687119284*fr[13]+1.060660171779821*fr[12]-0.547722557505166*fr[11]+1.10227038425243*fr[10]+0.3162277660168378*fr[9]+0.7905694150420947*fr[8]+0.3162277660168378*fr[7]+0.821583836257749*fr[6]-0.6363961030678926*fr[5]-0.821583836257749*fr[4]-0.4743416490252567*fr[3]-0.6123724356957944*fr[2]+0.4743416490252567*fr[1]+0.3535533905932737*fr[0]; 
  }
  if (alphaR[6]+0.7453559924999286*alphaR[4]+1.5*alphaR[3]+1.118033988749893*alphaR[2]+1.118033988749893*alphaR[1]+0.8333333333333334*alphaR[0] > 0.) {
    fUpOrdR[8] = 0.7348469228349532*fc[19]+1.423024947075771*fc[18]+0.7348469228349532*fc[17]+0.547722557505166*fc[16]+0.4242640687119284*fc[15]+1.060660171779821*fc[14]+0.4242640687119284*fc[13]+1.060660171779821*fc[12]+0.547722557505166*fc[11]+1.10227038425243*fc[10]+0.3162277660168378*fc[9]+0.7905694150420947*fc[8]+0.3162277660168378*fc[7]+0.821583836257749*fc[6]+0.6363961030678926*fc[5]+0.821583836257749*fc[4]+0.4743416490252567*fc[3]+0.6123724356957944*fc[2]+0.4743416490252567*fc[1]+0.3535533905932737*fc[0]; 
  } else {
    fUpOrdR[8] = (-0.7348469228349532*fr[19])+1.423024947075771*fr[18]-0.7348469228349532*fr[17]-0.547722557505166*fr[16]+0.4242640687119284*fr[15]+1.060660171779821*fr[14]+0.4242640687119284*fr[13]+1.060660171779821*fr[12]-0.547722557505166*fr[11]-1.10227038425243*fr[10]+0.3162277660168378*fr[9]+0.7905694150420947*fr[8]+0.3162277660168378*fr[7]-0.821583836257749*fr[6]+0.6363961030678926*fr[5]-0.821583836257749*fr[4]+0.4743416490252567*fr[3]-0.6123724356957944*fr[2]+0.4743416490252567*fr[1]+0.3535533905932737*fr[0]; 
  }

  double fUpR[8] = {0.};
  fUpR[0] = 0.154320987654321*fUpOrdR[8]+0.154320987654321*fUpOrdR[7]+0.2469135802469136*fUpOrdR[6]+0.154320987654321*fUpOrdR[5]+0.154320987654321*fUpOrdR[4]+0.2469135802469136*fUpOrdR[3]+0.2469135802469136*fUpOrdR[2]+0.2469135802469136*fUpOrdR[1]+0.3950617283950617*fUpOrdR[0]; 
  fUpR[1] = 0.2070433312499804*fUpOrdR[8]+0.2070433312499804*fUpOrdR[7]+0.3312693299999682*fUpOrdR[6]-0.2070433312499804*fUpOrdR[5]-0.2070433312499804*fUpOrdR[4]-0.3312693299999682*fUpOrdR[3]; 
  fUpR[2] = 0.2070433312499804*fUpOrdR[8]-0.2070433312499804*fUpOrdR[7]+0.2070433312499804*fUpOrdR[5]-0.2070433312499804*fUpOrdR[4]+0.3312693299999682*fUpOrdR[2]-0.3312693299999682*fUpOrdR[1]; 
  fUpR[3] = 0.2777777777777778*fUpOrdR[8]-0.2777777777777778*fUpOrdR[7]-0.2777777777777778*fUpOrdR[5]+0.2777777777777778*fUpOrdR[4]; 
  fUpR[4] = 0.138028887499987*fUpOrdR[8]+0.138028887499987*fUpOrdR[7]+0.2208462199999792*fUpOrdR[6]+0.138028887499987*fUpOrdR[5]+0.138028887499987*fUpOrdR[4]+0.2208462199999792*fUpOrdR[3]-0.2760577749999741*fUpOrdR[2]-0.2760577749999741*fUpOrdR[1]-0.4416924399999584*fUpOrdR[0]; 
  fUpR[5] = 0.138028887499987*fUpOrdR[8]+0.138028887499987*fUpOrdR[7]-0.2760577749999741*fUpOrdR[6]+0.138028887499987*fUpOrdR[5]+0.138028887499987*fUpOrdR[4]-0.2760577749999741*fUpOrdR[3]+0.2208462199999792*fUpOrdR[2]+0.2208462199999792*fUpOrdR[1]-0.4416924399999584*fUpOrdR[0]; 
  fUpR[6] = 0.1851851851851852*fUpOrdR[8]-0.1851851851851852*fUpOrdR[7]+0.1851851851851852*fUpOrdR[5]-0.1851851851851852*fUpOrdR[4]-0.3703703703703705*fUpOrdR[2]+0.3703703703703705*fUpOrdR[1]; 
  fUpR[7] = 0.1851851851851852*fUpOrdR[8]+0.1851851851851852*fUpOrdR[7]-0.3703703703703705*fUpOrdR[6]-0.1851851851851852*fUpOrdR[5]-0.1851851851851852*fUpOrdR[4]+0.3703703703703705*fUpOrdR[3]; 

  double GhatR[20] = {0.}; 
  GhatR[0] += 0.5*alphaR[6]*fUpR[6]+0.5*alphaR[4]*fUpR[4]+0.5*alphaR[3]*fUpR[3]+0.5*alphaR[2]*fUpR[2]+0.5*alphaR[1]*fUpR[1]+0.5*alphaR[0]*fUpR[0]; 
  GhatR[1] += 0.447213595499958*alphaR[3]*fUpR[6]+0.447213595499958*fUpR[3]*alphaR[6]+0.4472135954999579*alphaR[1]*fUpR[4]+0.4472135954999579*fUpR[1]*alphaR[4]+0.5*alphaR[2]*fUpR[3]+0.5*fUpR[2]*alphaR[3]+0.5*alphaR[0]*fUpR[1]+0.5*fUpR[0]*alphaR[1]; 
  GhatR[2] += 0.447213595499958*alphaR[3]*fUpR[7]+0.5000000000000001*alphaR[4]*fUpR[6]+0.5000000000000001*fUpR[4]*alphaR[6]+0.4472135954999579*alphaR[2]*fUpR[5]+0.5*alphaR[1]*fUpR[3]+0.5*fUpR[1]*alphaR[3]+0.5*alphaR[0]*fUpR[2]+0.5*fUpR[0]*alphaR[2]; 
  GhatR[3] += 0.4*alphaR[6]*fUpR[7]+0.447213595499958*alphaR[2]*fUpR[7]+0.447213595499958*alphaR[1]*fUpR[6]+0.447213595499958*fUpR[1]*alphaR[6]+0.4472135954999579*alphaR[3]*fUpR[5]+0.4472135954999579*alphaR[3]*fUpR[4]+0.4472135954999579*fUpR[3]*alphaR[4]+0.5*alphaR[0]*fUpR[3]+0.5*fUpR[0]*alphaR[3]+0.5*alphaR[1]*fUpR[2]+0.5*fUpR[1]*alphaR[2]; 
  GhatR[4] += 0.31943828249997*alphaR[6]*fUpR[6]+0.5000000000000001*alphaR[2]*fUpR[6]+0.5000000000000001*fUpR[2]*alphaR[6]+0.31943828249997*alphaR[4]*fUpR[4]+0.5*alphaR[0]*fUpR[4]+0.5*fUpR[0]*alphaR[4]+0.4472135954999579*alphaR[3]*fUpR[3]+0.4472135954999579*alphaR[1]*fUpR[1]; 
  GhatR[5] += 0.5000000000000001*alphaR[1]*fUpR[7]+0.4472135954999579*alphaR[6]*fUpR[6]+0.5*alphaR[0]*fUpR[5]+0.4472135954999579*alphaR[3]*fUpR[3]+0.4472135954999579*alphaR[2]*fUpR[2]; 
  GhatR[6] += 0.4*alphaR[3]*fUpR[7]+0.31943828249997*alphaR[4]*fUpR[6]+0.5*alphaR[0]*fUpR[6]+0.4472135954999579*fUpR[5]*alphaR[6]+0.31943828249997*fUpR[4]*alphaR[6]+0.5*fUpR[0]*alphaR[6]+0.5000000000000001*alphaR[2]*fUpR[4]+0.5000000000000001*fUpR[2]*alphaR[4]+0.447213595499958*alphaR[1]*fUpR[3]+0.447213595499958*fUpR[1]*alphaR[3]; 
  GhatR[7] += 0.4472135954999579*alphaR[4]*fUpR[7]+0.5*alphaR[0]*fUpR[7]+0.4*alphaR[3]*fUpR[6]+0.4*fUpR[3]*alphaR[6]+0.5000000000000001*alphaR[1]*fUpR[5]+0.447213595499958*alphaR[2]*fUpR[3]+0.447213595499958*fUpR[2]*alphaR[3]; 

  out[0] += (0.7071067811865475*GhatL[0]-0.7071067811865475*GhatR[0])*rdvpar2; 
  out[1] += (0.7071067811865475*GhatL[1]-0.7071067811865475*GhatR[1])*rdvpar2; 
  out[2] += ((-1.224744871391589*GhatR[0])-1.224744871391589*GhatL[0])*rdvpar2; 
  out[3] += (0.7071067811865475*GhatL[2]-0.7071067811865475*GhatR[2])*rdvpar2; 
  out[4] += ((-1.224744871391589*GhatR[1])-1.224744871391589*GhatL[1])*rdvpar2; 
  out[5] += (0.7071067811865475*GhatL[3]-0.7071067811865475*GhatR[3])*rdvpar2; 
  out[6] += ((-1.224744871391589*GhatR[2])-1.224744871391589*GhatL[2])*rdvpar2; 
  out[7] += (0.7071067811865475*GhatL[4]-0.7071067811865475*GhatR[4])*rdvpar2; 
  out[8] += (1.58113883008419*GhatL[0]-1.58113883008419*GhatR[0])*rdvpar2; 
  out[9] += (0.7071067811865475*GhatL[5]-0.7071067811865475*GhatR[5])*rdvpar2; 
  out[10] += ((-1.224744871391589*GhatR[3])-1.224744871391589*GhatL[3])*rdvpar2; 
  out[11] += ((-1.224744871391589*GhatR[4])-1.224744871391589*GhatL[4])*rdvpar2; 
  out[12] += (1.58113883008419*GhatL[1]-1.58113883008419*GhatR[1])*rdvpar2; 
  out[13] += (0.7071067811865475*GhatL[6]-0.7071067811865475*GhatR[6])*rdvpar2; 
  out[14] += (1.58113883008419*GhatL[2]-1.58113883008419*GhatR[2])*rdvpar2; 
  out[15] += (0.7071067811865475*GhatL[7]-0.7071067811865475*GhatR[7])*rdvpar2; 
  out[16] += ((-1.224744871391589*GhatR[5])-1.224744871391589*GhatL[5])*rdvpar2; 
  out[17] += ((-1.224744871391589*GhatR[6])-1.224744871391589*GhatL[6])*rdvpar2; 
  out[18] += (1.58113883008419*GhatL[3]-1.58113883008419*GhatR[3])*rdvpar2; 
  out[19] += ((-1.224744871391589*GhatR[7])-1.224744871391589*GhatL[7])*rdvpar2; 

} 
