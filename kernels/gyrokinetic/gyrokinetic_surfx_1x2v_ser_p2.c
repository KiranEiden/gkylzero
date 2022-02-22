#include <gkyl_gyrokinetic_kernels.h>
GKYL_CU_DH void gyrokinetic_surfx_1x2v_ser(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out) 
{ 
  // Apar: parallel component of magnetic vector potential.
  // dApardt: time derivative of Apar.
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_i: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // phi: electrostatic potential .
  // fL,fC,fR: distribution function in left, center and right cells.
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
  BstarZdBmag[0] = (1.414213562373095*(1.732050807568877*(2.23606797749979*jacobTotInv[1]*b_y[2]+jacobTotInv[0]*b_y[1])*m_*rdx2*wvpar+(cmag[2]*jacobTotInv[2]+cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0])*q_))/q_; 
  BstarZdBmag[1] = (0.2828427124746191*(1.732050807568877*(b_y[2]*(10.0*jacobTotInv[2]+11.18033988749895*jacobTotInv[0])+5.0*b_y[1]*jacobTotInv[1])*m_*rdx2*wvpar+(4.47213595499958*(cmag[1]*jacobTotInv[2]+jacobTotInv[1]*cmag[2])+5.0*(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1]))*q_))/q_; 
  BstarZdBmag[2] = (1.414213562373095*(2.23606797749979*jacobTotInv[1]*b_y[2]+jacobTotInv[0]*b_y[1])*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[4] = (1.414213562373095*(b_y[2]*(2.0*jacobTotInv[2]+2.23606797749979*jacobTotInv[0])+b_y[1]*jacobTotInv[1])*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[7] = (0.04040610178208843*(60.6217782649107*(b_y[1]*jacobTotInv[2]+2.0*jacobTotInv[1]*b_y[2])*m_*rdx2*wvpar+((22.3606797749979*cmag[2]+35.0*cmag[0])*jacobTotInv[2]+7.0*(5.0*jacobTotInv[0]*cmag[2]+4.47213595499958*cmag[1]*jacobTotInv[1]))*q_))/q_; 
  BstarZdBmag[11] = (1.414213562373095*(b_y[1]*jacobTotInv[2]+2.0*jacobTotInv[1]*b_y[2])*m_*rdx2)/(q_*rdvpar2); 

  double alphaL[8]; 
  alphaL[0] = (0.25*(hamil[8]*(8.660254037844387*BstarZdBmag[11]-6.708203932499369*BstarZdBmag[4]+3.872983346207417*BstarZdBmag[2])+hamil[2]*(3.872983346207417*BstarZdBmag[7]-3.0*BstarZdBmag[1]+1.732050807568877*BstarZdBmag[0]))*rdvpar2)/m_; 
  alphaL[1] = (0.25*(3.872983346207417*hamil[2]*BstarZdBmag[11]+(8.660254037844386*BstarZdBmag[7]-6.708203932499369*BstarZdBmag[1]+3.872983346207417*BstarZdBmag[0])*hamil[8]+hamil[2]*(1.732050807568877*BstarZdBmag[2]-3.0*BstarZdBmag[4]))*rdvpar2)/m_; 
  alphaL[4] = (0.5*hamil[8]*(3.872983346207417*BstarZdBmag[11]-3.0*BstarZdBmag[4]+1.732050807568877*BstarZdBmag[2])*rdvpar2)/m_; 

  double alphaR[8]; 
  alphaR[0] = (0.25*(hamil[8]*(8.660254037844387*BstarZdBmag[11]+6.708203932499369*BstarZdBmag[4]+3.872983346207417*BstarZdBmag[2])+hamil[2]*(3.872983346207417*BstarZdBmag[7]+3.0*BstarZdBmag[1]+1.732050807568877*BstarZdBmag[0]))*rdvpar2)/m_; 
  alphaR[1] = (0.25*(3.872983346207417*hamil[2]*BstarZdBmag[11]+(8.660254037844386*BstarZdBmag[7]+6.708203932499369*BstarZdBmag[1]+3.872983346207417*BstarZdBmag[0])*hamil[8]+hamil[2]*(3.0*BstarZdBmag[4]+1.732050807568877*BstarZdBmag[2]))*rdvpar2)/m_; 
  alphaR[4] = (0.5*hamil[8]*(3.872983346207417*BstarZdBmag[11]+3.0*BstarZdBmag[4]+1.732050807568877*BstarZdBmag[2])*rdvpar2)/m_; 

  double fUpOrdL[9];
  if (alphaL[0]-1.118033988749893*alphaL[4] > 0.) {
    fUpOrdL[0] = (-0.6846531968814574*(fL[15]+fL[12]))-0.3952847075210473*(fL[9]+fL[8])+0.7905694150420947*fL[7]+0.6123724356957944*fL[1]+0.3535533905932737*fL[0]; 
  } else {
    fUpOrdL[0] = 0.6846531968814574*(fC[15]+fC[12])-0.3952847075210473*(fC[9]+fC[8])+0.7905694150420947*fC[7]-0.6123724356957944*fC[1]+0.3535533905932737*fC[0]; 
  }
  if (alphaL[0]-1.118033988749893*alphaL[4] > 0.) {
    fUpOrdL[1] = 0.9185586535436915*fL[18]+0.547722557505166*fL[15]+0.5303300858899105*fL[14]-1.060660171779821*fL[13]-0.6846531968814574*fL[12]+0.3162277660168378*fL[9]-0.3952847075210473*fL[8]+0.7905694150420947*fL[7]-0.821583836257749*fL[5]-0.4743416490252567*fL[3]+0.6123724356957944*fL[1]+0.3535533905932737*fL[0]; 
  } else {
    fUpOrdL[1] = (-0.9185586535436915*fC[18])-0.547722557505166*fC[15]+0.5303300858899105*fC[14]-1.060660171779821*fC[13]+0.6846531968814574*fC[12]+0.3162277660168378*fC[9]-0.3952847075210473*fC[8]+0.7905694150420947*fC[7]+0.821583836257749*fC[5]-0.4743416490252567*fC[3]-0.6123724356957944*fC[1]+0.3535533905932737*fC[0]; 
  }
  if (alphaL[0]-1.118033988749893*alphaL[4] > 0.) {
    fUpOrdL[2] = (-0.9185586535436915*fL[18])+0.547722557505166*fL[15]-0.5303300858899105*fL[14]+1.060660171779821*fL[13]-0.6846531968814574*fL[12]+0.3162277660168378*fL[9]-0.3952847075210473*fL[8]+0.7905694150420947*fL[7]+0.821583836257749*fL[5]+0.4743416490252567*fL[3]+0.6123724356957944*fL[1]+0.3535533905932737*fL[0]; 
  } else {
    fUpOrdL[2] = 0.9185586535436915*fC[18]-0.547722557505166*fC[15]-0.5303300858899105*fC[14]+1.060660171779821*fC[13]+0.6846531968814574*fC[12]+0.3162277660168378*fC[9]-0.3952847075210473*fC[8]+0.7905694150420947*fC[7]-0.821583836257749*fC[5]+0.4743416490252567*fC[3]-0.6123724356957944*fC[1]+0.3535533905932737*fC[0]; 
  }
  if (alphaL[4]-1.5*alphaL[1]+1.118033988749897*alphaL[0] > 0.) {
    fUpOrdL[3] = 0.9185586535436915*fL[19]+0.5303300858899105*fL[16]-0.6846531968814574*fL[15]+0.547722557505166*fL[12]-1.060660171779821*fL[11]-0.3952847075210473*fL[9]+0.3162277660168378*fL[8]+0.7905694150420947*fL[7]-0.821583836257749*fL[4]-0.4743416490252567*fL[2]+0.6123724356957944*fL[1]+0.3535533905932737*fL[0]; 
  } else {
    fUpOrdL[3] = (-0.9185586535436915*fC[19])+0.5303300858899105*fC[16]+0.6846531968814574*fC[15]-0.547722557505166*fC[12]-1.060660171779821*fC[11]-0.3952847075210473*fC[9]+0.3162277660168378*fC[8]+0.7905694150420947*fC[7]+0.821583836257749*fC[4]-0.4743416490252567*fC[2]-0.6123724356957944*fC[1]+0.3535533905932737*fC[0]; 
  }
  if (alphaL[4]-1.5*alphaL[1]+1.118033988749897*alphaL[0] > 0.) {
    fUpOrdL[4] = (-0.7348469228349532*(fL[19]+fL[18]))+1.423024947075771*fL[17]-0.4242640687119284*fL[16]+0.547722557505166*fL[15]-0.4242640687119284*fL[14]-1.060660171779821*fL[13]+0.547722557505166*fL[12]-1.060660171779821*fL[11]+1.10227038425243*fL[10]+0.3162277660168378*(fL[9]+fL[8])+0.7905694150420947*fL[7]+0.6363961030678926*fL[6]-0.821583836257749*(fL[5]+fL[4])-0.4743416490252567*(fL[3]+fL[2])+0.6123724356957944*fL[1]+0.3535533905932737*fL[0]; 
  } else {
    fUpOrdL[4] = 0.7348469228349532*(fC[19]+fC[18])+1.423024947075771*fC[17]-0.4242640687119284*fC[16]-0.547722557505166*fC[15]-0.4242640687119284*fC[14]-1.060660171779821*fC[13]-0.547722557505166*fC[12]-1.060660171779821*fC[11]-1.10227038425243*fC[10]+0.3162277660168378*(fC[9]+fC[8])+0.7905694150420947*fC[7]+0.6363961030678926*fC[6]+0.821583836257749*(fC[5]+fC[4])-0.4743416490252567*(fC[3]+fC[2])-0.6123724356957944*fC[1]+0.3535533905932737*fC[0]; 
  }
  if (alphaL[4]-1.5*alphaL[1]+1.118033988749897*alphaL[0] > 0.) {
    fUpOrdL[5] = (-0.7348469228349532*fL[19])+0.7348469228349532*fL[18]-1.423024947075771*fL[17]-0.4242640687119284*fL[16]+0.547722557505166*fL[15]+0.4242640687119284*fL[14]+1.060660171779821*fL[13]+0.547722557505166*fL[12]-1.060660171779821*fL[11]-1.10227038425243*fL[10]+0.3162277660168378*(fL[9]+fL[8])+0.7905694150420947*fL[7]-0.6363961030678926*fL[6]+0.821583836257749*fL[5]-0.821583836257749*fL[4]+0.4743416490252567*fL[3]-0.4743416490252567*fL[2]+0.6123724356957944*fL[1]+0.3535533905932737*fL[0]; 
  } else {
    fUpOrdL[5] = 0.7348469228349532*fC[19]-0.7348469228349532*fC[18]-1.423024947075771*fC[17]-0.4242640687119284*fC[16]-0.547722557505166*fC[15]+0.4242640687119284*fC[14]+1.060660171779821*fC[13]-0.547722557505166*fC[12]-1.060660171779821*fC[11]+1.10227038425243*fC[10]+0.3162277660168378*(fC[9]+fC[8])+0.7905694150420947*fC[7]-0.6363961030678926*fC[6]-0.821583836257749*fC[5]+0.821583836257749*fC[4]+0.4743416490252567*fC[3]-0.4743416490252567*fC[2]-0.6123724356957944*fC[1]+0.3535533905932737*fC[0]; 
  }
  if (alphaL[4]+1.5*alphaL[1]+1.118033988749897*alphaL[0] > 0.) {
    fUpOrdL[6] = (-0.9185586535436915*fL[19])-0.5303300858899105*fL[16]-0.6846531968814574*fL[15]+0.547722557505166*fL[12]+1.060660171779821*fL[11]-0.3952847075210473*fL[9]+0.3162277660168378*fL[8]+0.7905694150420947*fL[7]+0.821583836257749*fL[4]+0.4743416490252567*fL[2]+0.6123724356957944*fL[1]+0.3535533905932737*fL[0]; 
  } else {
    fUpOrdL[6] = 0.9185586535436915*fC[19]-0.5303300858899105*fC[16]+0.6846531968814574*fC[15]-0.547722557505166*fC[12]+1.060660171779821*fC[11]-0.3952847075210473*fC[9]+0.3162277660168378*fC[8]+0.7905694150420947*fC[7]-0.821583836257749*fC[4]+0.4743416490252567*fC[2]-0.6123724356957944*fC[1]+0.3535533905932737*fC[0]; 
  }
  if (alphaL[4]+1.5*alphaL[1]+1.118033988749897*alphaL[0] > 0.) {
    fUpOrdL[7] = 0.7348469228349532*fL[19]-0.7348469228349532*fL[18]-1.423024947075771*fL[17]+0.4242640687119284*fL[16]+0.547722557505166*fL[15]-0.4242640687119284*fL[14]-1.060660171779821*fL[13]+0.547722557505166*fL[12]+1.060660171779821*fL[11]-1.10227038425243*fL[10]+0.3162277660168378*(fL[9]+fL[8])+0.7905694150420947*fL[7]-0.6363961030678926*fL[6]-0.821583836257749*fL[5]+0.821583836257749*fL[4]-0.4743416490252567*fL[3]+0.4743416490252567*fL[2]+0.6123724356957944*fL[1]+0.3535533905932737*fL[0]; 
  } else {
    fUpOrdL[7] = (-0.7348469228349532*fC[19])+0.7348469228349532*fC[18]-1.423024947075771*fC[17]+0.4242640687119284*fC[16]-0.547722557505166*fC[15]-0.4242640687119284*fC[14]-1.060660171779821*fC[13]-0.547722557505166*fC[12]+1.060660171779821*fC[11]+1.10227038425243*fC[10]+0.3162277660168378*(fC[9]+fC[8])+0.7905694150420947*fC[7]-0.6363961030678926*fC[6]+0.821583836257749*fC[5]-0.821583836257749*fC[4]-0.4743416490252567*fC[3]+0.4743416490252567*fC[2]-0.6123724356957944*fC[1]+0.3535533905932737*fC[0]; 
  }
  if (alphaL[4]+1.5*alphaL[1]+1.118033988749897*alphaL[0] > 0.) {
    fUpOrdL[8] = 0.7348469228349532*(fL[19]+fL[18])+1.423024947075771*fL[17]+0.4242640687119284*fL[16]+0.547722557505166*fL[15]+0.4242640687119284*fL[14]+1.060660171779821*fL[13]+0.547722557505166*fL[12]+1.060660171779821*fL[11]+1.10227038425243*fL[10]+0.3162277660168378*(fL[9]+fL[8])+0.7905694150420947*fL[7]+0.6363961030678926*fL[6]+0.821583836257749*(fL[5]+fL[4])+0.4743416490252567*(fL[3]+fL[2])+0.6123724356957944*fL[1]+0.3535533905932737*fL[0]; 
  } else {
    fUpOrdL[8] = (-0.7348469228349532*(fC[19]+fC[18]))+1.423024947075771*fC[17]+0.4242640687119284*fC[16]-0.547722557505166*fC[15]+0.4242640687119284*fC[14]+1.060660171779821*fC[13]-0.547722557505166*fC[12]+1.060660171779821*fC[11]-1.10227038425243*fC[10]+0.3162277660168378*(fC[9]+fC[8])+0.7905694150420947*fC[7]+0.6363961030678926*fC[6]-0.821583836257749*(fC[5]+fC[4])+0.4743416490252567*(fC[3]+fC[2])-0.6123724356957944*fC[1]+0.3535533905932737*fC[0]; 
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
  GhatL[0] += 0.5*alphaL[4]*fUpL[4]+0.5*alphaL[1]*fUpL[1]+0.5*alphaL[0]*fUpL[0]; 
  GhatL[1] += 0.4472135954999579*alphaL[1]*fUpL[4]+0.4472135954999579*fUpL[1]*alphaL[4]+0.5*alphaL[0]*fUpL[1]+0.5*fUpL[0]*alphaL[1]; 
  GhatL[2] += 0.5000000000000001*alphaL[4]*fUpL[6]+0.5*alphaL[1]*fUpL[3]+0.5*alphaL[0]*fUpL[2]; 
  GhatL[3] += 0.447213595499958*alphaL[1]*fUpL[6]+0.4472135954999579*fUpL[3]*alphaL[4]+0.5*alphaL[0]*fUpL[3]+0.5*alphaL[1]*fUpL[2]; 
  GhatL[4] += 0.31943828249997*alphaL[4]*fUpL[4]+0.5*alphaL[0]*fUpL[4]+0.5*fUpL[0]*alphaL[4]+0.4472135954999579*alphaL[1]*fUpL[1]; 
  GhatL[5] += 0.5000000000000001*alphaL[1]*fUpL[7]+0.5*alphaL[0]*fUpL[5]; 
  GhatL[6] += 0.31943828249997*alphaL[4]*fUpL[6]+0.5*alphaL[0]*fUpL[6]+0.5000000000000001*fUpL[2]*alphaL[4]+0.447213595499958*alphaL[1]*fUpL[3]; 
  GhatL[7] += 0.4472135954999579*alphaL[4]*fUpL[7]+0.5*alphaL[0]*fUpL[7]+0.5000000000000001*alphaL[1]*fUpL[5]; 

  double fUpOrdR[9];
  if (alphaR[0]-1.118033988749893*alphaR[4] > 0.) {
    fUpOrdR[0] = (-0.6846531968814574*(fC[15]+fC[12]))-0.3952847075210473*(fC[9]+fC[8])+0.7905694150420947*fC[7]+0.6123724356957944*fC[1]+0.3535533905932737*fC[0]; 
  } else {
    fUpOrdR[0] = 0.6846531968814574*(fR[15]+fR[12])-0.3952847075210473*(fR[9]+fR[8])+0.7905694150420947*fR[7]-0.6123724356957944*fR[1]+0.3535533905932737*fR[0]; 
  }
  if (alphaR[0]-1.118033988749893*alphaR[4] > 0.) {
    fUpOrdR[1] = 0.9185586535436915*fC[18]+0.547722557505166*fC[15]+0.5303300858899105*fC[14]-1.060660171779821*fC[13]-0.6846531968814574*fC[12]+0.3162277660168378*fC[9]-0.3952847075210473*fC[8]+0.7905694150420947*fC[7]-0.821583836257749*fC[5]-0.4743416490252567*fC[3]+0.6123724356957944*fC[1]+0.3535533905932737*fC[0]; 
  } else {
    fUpOrdR[1] = (-0.9185586535436915*fR[18])-0.547722557505166*fR[15]+0.5303300858899105*fR[14]-1.060660171779821*fR[13]+0.6846531968814574*fR[12]+0.3162277660168378*fR[9]-0.3952847075210473*fR[8]+0.7905694150420947*fR[7]+0.821583836257749*fR[5]-0.4743416490252567*fR[3]-0.6123724356957944*fR[1]+0.3535533905932737*fR[0]; 
  }
  if (alphaR[0]-1.118033988749893*alphaR[4] > 0.) {
    fUpOrdR[2] = (-0.9185586535436915*fC[18])+0.547722557505166*fC[15]-0.5303300858899105*fC[14]+1.060660171779821*fC[13]-0.6846531968814574*fC[12]+0.3162277660168378*fC[9]-0.3952847075210473*fC[8]+0.7905694150420947*fC[7]+0.821583836257749*fC[5]+0.4743416490252567*fC[3]+0.6123724356957944*fC[1]+0.3535533905932737*fC[0]; 
  } else {
    fUpOrdR[2] = 0.9185586535436915*fR[18]-0.547722557505166*fR[15]-0.5303300858899105*fR[14]+1.060660171779821*fR[13]+0.6846531968814574*fR[12]+0.3162277660168378*fR[9]-0.3952847075210473*fR[8]+0.7905694150420947*fR[7]-0.821583836257749*fR[5]+0.4743416490252567*fR[3]-0.6123724356957944*fR[1]+0.3535533905932737*fR[0]; 
  }
  if (alphaR[4]-1.5*alphaR[1]+1.118033988749897*alphaR[0] > 0.) {
    fUpOrdR[3] = 0.9185586535436915*fC[19]+0.5303300858899105*fC[16]-0.6846531968814574*fC[15]+0.547722557505166*fC[12]-1.060660171779821*fC[11]-0.3952847075210473*fC[9]+0.3162277660168378*fC[8]+0.7905694150420947*fC[7]-0.821583836257749*fC[4]-0.4743416490252567*fC[2]+0.6123724356957944*fC[1]+0.3535533905932737*fC[0]; 
  } else {
    fUpOrdR[3] = (-0.9185586535436915*fR[19])+0.5303300858899105*fR[16]+0.6846531968814574*fR[15]-0.547722557505166*fR[12]-1.060660171779821*fR[11]-0.3952847075210473*fR[9]+0.3162277660168378*fR[8]+0.7905694150420947*fR[7]+0.821583836257749*fR[4]-0.4743416490252567*fR[2]-0.6123724356957944*fR[1]+0.3535533905932737*fR[0]; 
  }
  if (alphaR[4]-1.5*alphaR[1]+1.118033988749897*alphaR[0] > 0.) {
    fUpOrdR[4] = (-0.7348469228349532*(fC[19]+fC[18]))+1.423024947075771*fC[17]-0.4242640687119284*fC[16]+0.547722557505166*fC[15]-0.4242640687119284*fC[14]-1.060660171779821*fC[13]+0.547722557505166*fC[12]-1.060660171779821*fC[11]+1.10227038425243*fC[10]+0.3162277660168378*(fC[9]+fC[8])+0.7905694150420947*fC[7]+0.6363961030678926*fC[6]-0.821583836257749*(fC[5]+fC[4])-0.4743416490252567*(fC[3]+fC[2])+0.6123724356957944*fC[1]+0.3535533905932737*fC[0]; 
  } else {
    fUpOrdR[4] = 0.7348469228349532*(fR[19]+fR[18])+1.423024947075771*fR[17]-0.4242640687119284*fR[16]-0.547722557505166*fR[15]-0.4242640687119284*fR[14]-1.060660171779821*fR[13]-0.547722557505166*fR[12]-1.060660171779821*fR[11]-1.10227038425243*fR[10]+0.3162277660168378*(fR[9]+fR[8])+0.7905694150420947*fR[7]+0.6363961030678926*fR[6]+0.821583836257749*(fR[5]+fR[4])-0.4743416490252567*(fR[3]+fR[2])-0.6123724356957944*fR[1]+0.3535533905932737*fR[0]; 
  }
  if (alphaR[4]-1.5*alphaR[1]+1.118033988749897*alphaR[0] > 0.) {
    fUpOrdR[5] = (-0.7348469228349532*fC[19])+0.7348469228349532*fC[18]-1.423024947075771*fC[17]-0.4242640687119284*fC[16]+0.547722557505166*fC[15]+0.4242640687119284*fC[14]+1.060660171779821*fC[13]+0.547722557505166*fC[12]-1.060660171779821*fC[11]-1.10227038425243*fC[10]+0.3162277660168378*(fC[9]+fC[8])+0.7905694150420947*fC[7]-0.6363961030678926*fC[6]+0.821583836257749*fC[5]-0.821583836257749*fC[4]+0.4743416490252567*fC[3]-0.4743416490252567*fC[2]+0.6123724356957944*fC[1]+0.3535533905932737*fC[0]; 
  } else {
    fUpOrdR[5] = 0.7348469228349532*fR[19]-0.7348469228349532*fR[18]-1.423024947075771*fR[17]-0.4242640687119284*fR[16]-0.547722557505166*fR[15]+0.4242640687119284*fR[14]+1.060660171779821*fR[13]-0.547722557505166*fR[12]-1.060660171779821*fR[11]+1.10227038425243*fR[10]+0.3162277660168378*(fR[9]+fR[8])+0.7905694150420947*fR[7]-0.6363961030678926*fR[6]-0.821583836257749*fR[5]+0.821583836257749*fR[4]+0.4743416490252567*fR[3]-0.4743416490252567*fR[2]-0.6123724356957944*fR[1]+0.3535533905932737*fR[0]; 
  }
  if (alphaR[4]+1.5*alphaR[1]+1.118033988749897*alphaR[0] > 0.) {
    fUpOrdR[6] = (-0.9185586535436915*fC[19])-0.5303300858899105*fC[16]-0.6846531968814574*fC[15]+0.547722557505166*fC[12]+1.060660171779821*fC[11]-0.3952847075210473*fC[9]+0.3162277660168378*fC[8]+0.7905694150420947*fC[7]+0.821583836257749*fC[4]+0.4743416490252567*fC[2]+0.6123724356957944*fC[1]+0.3535533905932737*fC[0]; 
  } else {
    fUpOrdR[6] = 0.9185586535436915*fR[19]-0.5303300858899105*fR[16]+0.6846531968814574*fR[15]-0.547722557505166*fR[12]+1.060660171779821*fR[11]-0.3952847075210473*fR[9]+0.3162277660168378*fR[8]+0.7905694150420947*fR[7]-0.821583836257749*fR[4]+0.4743416490252567*fR[2]-0.6123724356957944*fR[1]+0.3535533905932737*fR[0]; 
  }
  if (alphaR[4]+1.5*alphaR[1]+1.118033988749897*alphaR[0] > 0.) {
    fUpOrdR[7] = 0.7348469228349532*fC[19]-0.7348469228349532*fC[18]-1.423024947075771*fC[17]+0.4242640687119284*fC[16]+0.547722557505166*fC[15]-0.4242640687119284*fC[14]-1.060660171779821*fC[13]+0.547722557505166*fC[12]+1.060660171779821*fC[11]-1.10227038425243*fC[10]+0.3162277660168378*(fC[9]+fC[8])+0.7905694150420947*fC[7]-0.6363961030678926*fC[6]-0.821583836257749*fC[5]+0.821583836257749*fC[4]-0.4743416490252567*fC[3]+0.4743416490252567*fC[2]+0.6123724356957944*fC[1]+0.3535533905932737*fC[0]; 
  } else {
    fUpOrdR[7] = (-0.7348469228349532*fR[19])+0.7348469228349532*fR[18]-1.423024947075771*fR[17]+0.4242640687119284*fR[16]-0.547722557505166*fR[15]-0.4242640687119284*fR[14]-1.060660171779821*fR[13]-0.547722557505166*fR[12]+1.060660171779821*fR[11]+1.10227038425243*fR[10]+0.3162277660168378*(fR[9]+fR[8])+0.7905694150420947*fR[7]-0.6363961030678926*fR[6]+0.821583836257749*fR[5]-0.821583836257749*fR[4]-0.4743416490252567*fR[3]+0.4743416490252567*fR[2]-0.6123724356957944*fR[1]+0.3535533905932737*fR[0]; 
  }
  if (alphaR[4]+1.5*alphaR[1]+1.118033988749897*alphaR[0] > 0.) {
    fUpOrdR[8] = 0.7348469228349532*(fC[19]+fC[18])+1.423024947075771*fC[17]+0.4242640687119284*fC[16]+0.547722557505166*fC[15]+0.4242640687119284*fC[14]+1.060660171779821*fC[13]+0.547722557505166*fC[12]+1.060660171779821*fC[11]+1.10227038425243*fC[10]+0.3162277660168378*(fC[9]+fC[8])+0.7905694150420947*fC[7]+0.6363961030678926*fC[6]+0.821583836257749*(fC[5]+fC[4])+0.4743416490252567*(fC[3]+fC[2])+0.6123724356957944*fC[1]+0.3535533905932737*fC[0]; 
  } else {
    fUpOrdR[8] = (-0.7348469228349532*(fR[19]+fR[18]))+1.423024947075771*fR[17]+0.4242640687119284*fR[16]-0.547722557505166*fR[15]+0.4242640687119284*fR[14]+1.060660171779821*fR[13]-0.547722557505166*fR[12]+1.060660171779821*fR[11]-1.10227038425243*fR[10]+0.3162277660168378*(fR[9]+fR[8])+0.7905694150420947*fR[7]+0.6363961030678926*fR[6]-0.821583836257749*(fR[5]+fR[4])+0.4743416490252567*(fR[3]+fR[2])-0.6123724356957944*fR[1]+0.3535533905932737*fR[0]; 
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
  GhatR[0] += 0.5*alphaR[4]*fUpR[4]+0.5*alphaR[1]*fUpR[1]+0.5*alphaR[0]*fUpR[0]; 
  GhatR[1] += 0.4472135954999579*alphaR[1]*fUpR[4]+0.4472135954999579*fUpR[1]*alphaR[4]+0.5*alphaR[0]*fUpR[1]+0.5*fUpR[0]*alphaR[1]; 
  GhatR[2] += 0.5000000000000001*alphaR[4]*fUpR[6]+0.5*alphaR[1]*fUpR[3]+0.5*alphaR[0]*fUpR[2]; 
  GhatR[3] += 0.447213595499958*alphaR[1]*fUpR[6]+0.4472135954999579*fUpR[3]*alphaR[4]+0.5*alphaR[0]*fUpR[3]+0.5*alphaR[1]*fUpR[2]; 
  GhatR[4] += 0.31943828249997*alphaR[4]*fUpR[4]+0.5*alphaR[0]*fUpR[4]+0.5*fUpR[0]*alphaR[4]+0.4472135954999579*alphaR[1]*fUpR[1]; 
  GhatR[5] += 0.5000000000000001*alphaR[1]*fUpR[7]+0.5*alphaR[0]*fUpR[5]; 
  GhatR[6] += 0.31943828249997*alphaR[4]*fUpR[6]+0.5*alphaR[0]*fUpR[6]+0.5000000000000001*fUpR[2]*alphaR[4]+0.447213595499958*alphaR[1]*fUpR[3]; 
  GhatR[7] += 0.4472135954999579*alphaR[4]*fUpR[7]+0.5*alphaR[0]*fUpR[7]+0.5000000000000001*alphaR[1]*fUpR[5]; 

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

} 
