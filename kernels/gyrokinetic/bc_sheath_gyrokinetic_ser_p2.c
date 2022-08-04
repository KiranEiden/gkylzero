#include <gkyl_bc_sheath_gyrokinetic.h> 


void bc_sheath_gyrokinetic_reflectedf_lower_1x2v_ser_p2(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double q2Dm, const double *phi, const double *phiWall, const double *f, double *fRefl) 
{ 
  double vcutSq; long double xc, b, xbarVal, fac; 
  double fReflZMuQuad[8][8]; 
  

  vcutSq = 0.5*(3.16227766016838*phiWall[2]-3.16227766016838*phi[2]-2.449489742783178*phiWall[1]+2.449489742783178*phi[1]+1.414213562373095*phiWall[0]-1.414213562373095*phi[0])*q2Dm; 
  if(vcutSq <= vlowerSq) { // absorb (no reflection) 
  fRefl[0] = 0.0; 
  fRefl[1] = 0.0; 
  fRefl[2] = 0.0; 
  fRefl[3] = 0.0; 
  fRefl[4] = 0.0; 
  fRefl[5] = 0.0; 
  fRefl[6] = 0.0; 
  fRefl[7] = 0.0; 
  fRefl[8] = 0.0; 
  fRefl[9] = 0.0; 
  fRefl[10] = 0.0; 
  fRefl[11] = 0.0; 
  fRefl[12] = 0.0; 
  fRefl[13] = 0.0; 
  fRefl[14] = 0.0; 
  fRefl[15] = 0.0; 
  fRefl[16] = 0.0; 
  fRefl[17] = 0.0; 
  fRefl[18] = 0.0; 
  fRefl[19] = 0.0; 
  } else if (vcutSq > vupperSq) { // full reflection 
  fRefl[0] = f[0]; 
  fRefl[1] = f[1]; 
  fRefl[2] = f[2]; 
  fRefl[3] = f[3]; 
  fRefl[4] = f[4]; 
  fRefl[5] = f[5]; 
  fRefl[6] = f[6]; 
  fRefl[7] = f[7]; 
  fRefl[8] = f[8]; 
  fRefl[9] = f[9]; 
  fRefl[10] = f[10]; 
  fRefl[11] = f[11]; 
  fRefl[12] = f[12]; 
  fRefl[13] = f[13]; 
  fRefl[14] = f[14]; 
  fRefl[15] = f[15]; 
  fRefl[16] = f[16]; 
  fRefl[17] = f[17]; 
  fRefl[18] = f[18]; 
  fRefl[19] = f[19]; 
  } else { // partial reflection 
  xbarVal = (0.9622504486493765*(2.0*(9.0*(f[19]+f[17])-6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[6]+f[4])-3.0*f[10])-5.0*f[2])))/(4.47213595499958*(6.708203932499369*(f[15]+f[13])-5.0*(f[9]+f[7]))+3.0*(11.18033988749895*(f[3]+f[1])-15.0*f[5])-25.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if (-0.02*(4.47213595499958*(6.708203932499369*(f[15]+f[13])-5.0*(f[9]+f[7]))+3.0*(11.18033988749895*(f[3]+f[1])-15.0*f[5])-25.0*f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[0][0] = 0.0; 
  fReflZMuQuad[0][1] = 0.0; 
  fReflZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][0] = (-0.02*(4.47213595499958*(6.708203932499369*(f[15]+f[13])-5.0*(f[9]+f[7]))+3.0*(11.18033988749895*(f[3]+f[1])-15.0*f[5])-25.0*f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][1] = (-0.03333333333333333*(2.0*(9.0*(f[19]+f[17])-6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[6]+f[4])-3.0*f[10])-5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][2] = (0.1*(9.0*f[18]-6.708203932499369*(f[14]+f[12])+5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][0] = (-0.02*(4.47213595499958*(6.708203932499369*(f[15]+f[13])-5.0*(f[9]+f[7]))+3.0*(11.18033988749895*(f[3]+f[1])-15.0*f[5])-25.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][1] = (-0.03333333333333333*(2.0*(9.0*(f[19]+f[17])-6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[6]+f[4])-3.0*f[10])-5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][2] = (0.1*(9.0*f[18]-6.708203932499369*(f[14]+f[12])+5.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(45.0*f[17]+6.708203932499369*(4.0*f[16]-5.0*f[11])+6.0*(5.0*f[2]-6.708203932499369*f[6])))/(2.23606797749979*(6.708203932499369*f[13]+4.0*f[9]-5.0*f[7])+2.0*(5.0*f[0]-6.708203932499369*f[3])); 
  // if f is not realizable, no reflection from this node 
  if (0.05*(2.23606797749979*(6.708203932499369*f[13]+4.0*f[9]-5.0*f[7])+2.0*(5.0*f[0]-6.708203932499369*f[3])) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[1][0] = 0.0; 
  fReflZMuQuad[1][1] = 0.0; 
  fReflZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][0] = (0.05*(2.23606797749979*(6.708203932499369*f[13]+4.0*f[9]-5.0*f[7])+2.0*(5.0*f[0]-6.708203932499369*f[3])))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][1] = (0.01666666666666667*(45.0*f[17]+6.708203932499369*(4.0*f[16]-5.0*f[11])+6.0*(5.0*f[2]-6.708203932499369*f[6])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][2] = (-0.1*(6.708203932499369*f[14]-5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][0] = (0.05*(2.23606797749979*(6.708203932499369*f[13]+4.0*f[9]-5.0*f[7])+2.0*(5.0*f[0]-6.708203932499369*f[3])))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][1] = (0.01666666666666667*(45.0*f[17]+6.708203932499369*(4.0*f[16]-5.0*f[11])+6.0*(5.0*f[2]-6.708203932499369*f[6])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][2] = (-0.1*(6.708203932499369*f[14]-5.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(9.0*(f[19]-1.0*f[17])+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[4]-1.0*f[6])-3.0*f[10])+5.0*f[2])))/(4.47213595499958*(6.708203932499369*(f[15]-1.0*f[13])+5.0*(f[9]+f[7]))+3.0*(11.18033988749895*(f[1]-1.0*f[3])-15.0*f[5])+25.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if (0.02*(4.47213595499958*(6.708203932499369*(f[15]-1.0*f[13])+5.0*(f[9]+f[7]))+3.0*(11.18033988749895*(f[1]-1.0*f[3])-15.0*f[5])+25.0*f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[2][0] = 0.0; 
  fReflZMuQuad[2][1] = 0.0; 
  fReflZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][0] = (0.02*(4.47213595499958*(6.708203932499369*(f[15]-1.0*f[13])+5.0*(f[9]+f[7]))+3.0*(11.18033988749895*(f[1]-1.0*f[3])-15.0*f[5])+25.0*f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][1] = (0.03333333333333333*(2.0*(9.0*(f[19]-1.0*f[17])+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[4]-1.0*f[6])-3.0*f[10])+5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][2] = (-0.1*(9.0*f[18]+6.708203932499369*f[14]-1.0*(6.708203932499369*f[12]+5.0*f[8])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][0] = (0.02*(4.47213595499958*(6.708203932499369*(f[15]-1.0*f[13])+5.0*(f[9]+f[7]))+3.0*(11.18033988749895*(f[1]-1.0*f[3])-15.0*f[5])+25.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][1] = (0.03333333333333333*(2.0*(9.0*(f[19]-1.0*f[17])+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[4]-1.0*f[6])-3.0*f[10])+5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][2] = (-0.1*(9.0*f[18]+6.708203932499369*f[14]-1.0*(6.708203932499369*f[12]+5.0*f[8])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(5.0*(9.0*f[19]-6.708203932499369*f[16])+2.0*(13.41640786499874*f[11]+3.0*(5.0*f[2]-6.708203932499369*f[4]))))/(2.23606797749979*(6.708203932499369*f[15]-5.0*f[9])+2.0*(2.23606797749979*(2.0*f[7]-3.0*f[1])+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if (0.05*(2.23606797749979*(6.708203932499369*f[15]-5.0*f[9])+2.0*(2.23606797749979*(2.0*f[7]-3.0*f[1])+5.0*f[0])) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[3][0] = 0.0; 
  fReflZMuQuad[3][1] = 0.0; 
  fReflZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][0] = (0.05*(2.23606797749979*(6.708203932499369*f[15]-5.0*f[9])+2.0*(2.23606797749979*(2.0*f[7]-3.0*f[1])+5.0*f[0])))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][1] = (0.01666666666666667*(5.0*(9.0*f[19]-6.708203932499369*f[16])+2.0*(13.41640786499874*f[11]+3.0*(5.0*f[2]-6.708203932499369*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][2] = (-0.1*(6.708203932499369*f[12]-5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][0] = (0.05*(2.23606797749979*(6.708203932499369*f[15]-5.0*f[9])+2.0*(2.23606797749979*(2.0*f[7]-3.0*f[1])+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][1] = (0.01666666666666667*(5.0*(9.0*f[19]-6.708203932499369*f[16])+2.0*(13.41640786499874*f[11]+3.0*(5.0*f[2]-6.708203932499369*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][2] = (-0.1*(6.708203932499369*f[12]-5.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(5.0*(9.0*f[19]+6.708203932499369*f[16])-2.0*(13.41640786499874*f[11]+3.0*(6.708203932499369*f[4]+5.0*f[2]))))/(2.23606797749979*(6.708203932499369*f[15]+5.0*f[9])-2.0*(2.23606797749979*(2.0*f[7]+3.0*f[1])+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if (-0.05*(2.23606797749979*(6.708203932499369*f[15]+5.0*f[9])-2.0*(2.23606797749979*(2.0*f[7]+3.0*f[1])+5.0*f[0])) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[4][0] = 0.0; 
  fReflZMuQuad[4][1] = 0.0; 
  fReflZMuQuad[4][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][0] = (-0.05*(2.23606797749979*(6.708203932499369*f[15]+5.0*f[9])-2.0*(2.23606797749979*(2.0*f[7]+3.0*f[1])+5.0*f[0])))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][1] = (-0.01666666666666667*(5.0*(9.0*f[19]+6.708203932499369*f[16])-2.0*(13.41640786499874*f[11]+3.0*(6.708203932499369*f[4]+5.0*f[2]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][2] = (0.1*(6.708203932499369*f[12]+5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][0] = (-0.05*(2.23606797749979*(6.708203932499369*f[15]+5.0*f[9])-2.0*(2.23606797749979*(2.0*f[7]+3.0*f[1])+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][1] = (-0.01666666666666667*(5.0*(9.0*f[19]+6.708203932499369*f[16])-2.0*(13.41640786499874*f[11]+3.0*(6.708203932499369*f[4]+5.0*f[2]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][2] = (0.1*(6.708203932499369*f[12]+5.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(9.0*f[19]-1.0*(9.0*f[17]+6.708203932499369*(f[16]+f[11])))+3.0*(3.0*(3.0*f[10]+2.23606797749979*(f[4]-1.0*f[6]))-5.0*f[2])))/(4.47213595499958*(6.708203932499369*f[15]-1.0*(6.708203932499369*f[13]+5.0*(f[9]+f[7])))+3.0*(15.0*f[5]+11.18033988749895*(f[1]-1.0*f[3]))-25.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if (-0.02*(4.47213595499958*(6.708203932499369*f[15]-1.0*(6.708203932499369*f[13]+5.0*(f[9]+f[7])))+3.0*(15.0*f[5]+11.18033988749895*(f[1]-1.0*f[3]))-25.0*f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[5][0] = 0.0; 
  fReflZMuQuad[5][1] = 0.0; 
  fReflZMuQuad[5][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][0] = (-0.02*(4.47213595499958*(6.708203932499369*f[15]-1.0*(6.708203932499369*f[13]+5.0*(f[9]+f[7])))+3.0*(15.0*f[5]+11.18033988749895*(f[1]-1.0*f[3]))-25.0*f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][1] = (-0.03333333333333333*(2.0*(9.0*f[19]-1.0*(9.0*f[17]+6.708203932499369*(f[16]+f[11])))+3.0*(3.0*(3.0*f[10]+2.23606797749979*(f[4]-1.0*f[6]))-5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][2] = (-0.1*(9.0*f[18]+6.708203932499369*(f[12]-1.0*f[14])-5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][0] = (-0.02*(4.47213595499958*(6.708203932499369*f[15]-1.0*(6.708203932499369*f[13]+5.0*(f[9]+f[7])))+3.0*(15.0*f[5]+11.18033988749895*(f[1]-1.0*f[3]))-25.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][1] = (-0.03333333333333333*(2.0*(9.0*f[19]-1.0*(9.0*f[17]+6.708203932499369*(f[16]+f[11])))+3.0*(3.0*(3.0*f[10]+2.23606797749979*(f[4]-1.0*f[6]))-5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][2] = (-0.1*(9.0*f[18]+6.708203932499369*(f[12]-1.0*f[14])-5.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(45.0*f[17]+6.708203932499369*(5.0*f[11]-4.0*f[16])-6.0*(6.708203932499369*f[6]+5.0*f[2])))/(2.23606797749979*(6.708203932499369*f[13]-4.0*f[9]+5.0*f[7])-2.0*(6.708203932499369*f[3]+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if (-0.05*(2.23606797749979*(6.708203932499369*f[13]-4.0*f[9]+5.0*f[7])-2.0*(6.708203932499369*f[3]+5.0*f[0])) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[6][0] = 0.0; 
  fReflZMuQuad[6][1] = 0.0; 
  fReflZMuQuad[6][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][0] = (-0.05*(2.23606797749979*(6.708203932499369*f[13]-4.0*f[9]+5.0*f[7])-2.0*(6.708203932499369*f[3]+5.0*f[0])))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][1] = (-0.01666666666666667*(45.0*f[17]+6.708203932499369*(5.0*f[11]-4.0*f[16])-6.0*(6.708203932499369*f[6]+5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][2] = (0.1*(6.708203932499369*f[14]+5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][0] = (-0.05*(2.23606797749979*(6.708203932499369*f[13]-4.0*f[9]+5.0*f[7])-2.0*(6.708203932499369*f[3]+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][1] = (-0.01666666666666667*(45.0*f[17]+6.708203932499369*(5.0*f[11]-4.0*f[16])-6.0*(6.708203932499369*f[6]+5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][2] = (0.1*(6.708203932499369*f[14]+5.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(9.0*(f[19]+f[17])+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(3.0*f[10]+2.23606797749979*(f[6]+f[4]))+5.0*f[2])))/(4.47213595499958*(6.708203932499369*(f[15]+f[13])+5.0*(f[9]+f[7]))+3.0*(15.0*f[5]+11.18033988749895*(f[3]+f[1]))+25.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if (0.02*(4.47213595499958*(6.708203932499369*(f[15]+f[13])+5.0*(f[9]+f[7]))+3.0*(15.0*f[5]+11.18033988749895*(f[3]+f[1]))+25.0*f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[7][0] = 0.0; 
  fReflZMuQuad[7][1] = 0.0; 
  fReflZMuQuad[7][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][0] = (0.02*(4.47213595499958*(6.708203932499369*(f[15]+f[13])+5.0*(f[9]+f[7]))+3.0*(15.0*f[5]+11.18033988749895*(f[3]+f[1]))+25.0*f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][1] = (0.03333333333333333*(2.0*(9.0*(f[19]+f[17])+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(3.0*f[10]+2.23606797749979*(f[6]+f[4]))+5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][2] = (0.1*(9.0*f[18]+6.708203932499369*(f[14]+f[12])+5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][0] = (0.02*(4.47213595499958*(6.708203932499369*(f[15]+f[13])+5.0*(f[9]+f[7]))+3.0*(15.0*f[5]+11.18033988749895*(f[3]+f[1]))+25.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][1] = (0.03333333333333333*(2.0*(9.0*(f[19]+f[17])+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(3.0*f[10]+2.23606797749979*(f[6]+f[4]))+5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][2] = (0.1*(9.0*f[18]+6.708203932499369*(f[14]+f[12])+5.0*f[8]))*fac; 
   } 
  } 
  fRefl[0] = 0.05555555555555555*(fReflZMuQuad[7][0]+8.0*fReflZMuQuad[6][0]+fReflZMuQuad[5][0]+8.0*(fReflZMuQuad[4][0]+fReflZMuQuad[3][0])+fReflZMuQuad[2][0]+8.0*fReflZMuQuad[1][0]+fReflZMuQuad[0][0]); 
  fRefl[1] = 4.469298760204439e-16*(4.63256860547201e+14*(fReflZMuQuad[7][0]-1.0*fReflZMuQuad[5][0])+7.4121097687552e+14*(fReflZMuQuad[4][0]-1.0*fReflZMuQuad[3][0])+4.63256860547201e+14*(fReflZMuQuad[2][0]-1.0*fReflZMuQuad[0][0])); 
  fRefl[2] = 0.05555555555555555*(fReflZMuQuad[7][1]+8.0*fReflZMuQuad[6][1]+fReflZMuQuad[5][1]+8.0*(fReflZMuQuad[4][1]+fReflZMuQuad[3][1])+fReflZMuQuad[2][1]+8.0*fReflZMuQuad[1][1]+fReflZMuQuad[0][1]); 
  fRefl[3] = 4.469298760204439e-16*(4.63256860547201e+14*fReflZMuQuad[7][0]+7.4121097687552e+14*fReflZMuQuad[6][0]+4.63256860547201e+14*fReflZMuQuad[5][0]-1.0*(4.63256860547201e+14*fReflZMuQuad[2][0]+7.4121097687552e+14*fReflZMuQuad[1][0]+4.63256860547201e+14*fReflZMuQuad[0][0])); 
  fRefl[4] = 4.469298760204439e-16*(4.63256860547201e+14*(fReflZMuQuad[7][1]-1.0*fReflZMuQuad[5][1])+7.4121097687552e+14*(fReflZMuQuad[4][1]-1.0*fReflZMuQuad[3][1])+4.63256860547201e+14*(fReflZMuQuad[2][1]-1.0*fReflZMuQuad[0][1])); 
  fRefl[5] = 0.2777777777777778*(fReflZMuQuad[7][0]-1.0*(fReflZMuQuad[5][0]+fReflZMuQuad[2][0])+fReflZMuQuad[0][0]); 
  fRefl[6] = 4.469298760204439e-16*(4.63256860547201e+14*fReflZMuQuad[7][1]+7.4121097687552e+14*fReflZMuQuad[6][1]+4.63256860547201e+14*fReflZMuQuad[5][1]-1.0*(4.63256860547201e+14*fReflZMuQuad[2][1]+7.4121097687552e+14*fReflZMuQuad[1][1]+4.63256860547201e+14*fReflZMuQuad[0][1])); 
  fRefl[7] = 0.2484519974999762*(fReflZMuQuad[7][0]-2.0*fReflZMuQuad[6][0]+fReflZMuQuad[5][0]+fReflZMuQuad[2][0]-2.0*fReflZMuQuad[1][0]+fReflZMuQuad[0][0]); 
  fRefl[8] = 0.05555555555555555*(fReflZMuQuad[7][2]+8.0*fReflZMuQuad[6][2]+fReflZMuQuad[5][2]+8.0*(fReflZMuQuad[4][2]+fReflZMuQuad[3][2])+fReflZMuQuad[2][2]+8.0*fReflZMuQuad[1][2]+fReflZMuQuad[0][2]); 
  fRefl[9] = 0.2484519974999762*(fReflZMuQuad[7][0]+fReflZMuQuad[5][0]-2.0*(fReflZMuQuad[4][0]+fReflZMuQuad[3][0])+fReflZMuQuad[2][0]+fReflZMuQuad[0][0]); 
  fRefl[10] = 0.2777777777777778*(fReflZMuQuad[7][1]-1.0*(fReflZMuQuad[5][1]+fReflZMuQuad[2][1])+fReflZMuQuad[0][1]); 
  fRefl[11] = 0.2484519974999762*(fReflZMuQuad[7][1]-2.0*fReflZMuQuad[6][1]+fReflZMuQuad[5][1]+fReflZMuQuad[2][1]-2.0*fReflZMuQuad[1][1]+fReflZMuQuad[0][1]); 
  fRefl[12] = 4.46929876020444e-16*(4.63256860547201e+14*(fReflZMuQuad[7][2]-1.0*fReflZMuQuad[5][2])+7.4121097687552e+14*fReflZMuQuad[4][2]+4.63256860547201e+14*(fReflZMuQuad[2][2]-1.0*fReflZMuQuad[0][2])); 
  fRefl[13] = 0.1851851851851852*(fReflZMuQuad[7][0]-2.0*fReflZMuQuad[6][0]+fReflZMuQuad[5][0]-1.0*fReflZMuQuad[2][0]+2.0*fReflZMuQuad[1][0]-1.0*fReflZMuQuad[0][0]); 
  fRefl[14] = 4.46929876020444e-16*(4.63256860547201e+14*fReflZMuQuad[7][2]+7.4121097687552e+14*fReflZMuQuad[6][2]+4.63256860547201e+14*fReflZMuQuad[5][2]-1.0*(4.63256860547201e+14*fReflZMuQuad[2][2]+7.4121097687552e+14*fReflZMuQuad[1][2]+4.63256860547201e+14*fReflZMuQuad[0][2])); 
  fRefl[15] = 0.1851851851851852*(fReflZMuQuad[7][0]-1.0*fReflZMuQuad[5][0]+2.0*(fReflZMuQuad[3][0]-1.0*fReflZMuQuad[4][0])+fReflZMuQuad[2][0]-1.0*fReflZMuQuad[0][0]); 
  fRefl[16] = 0.2484519974999762*(fReflZMuQuad[7][1]+fReflZMuQuad[5][1]-2.0*(fReflZMuQuad[4][1]+fReflZMuQuad[3][1])+fReflZMuQuad[2][1]+fReflZMuQuad[0][1]); 
  fRefl[17] = 0.1851851851851853*(fReflZMuQuad[7][1]-2.0*fReflZMuQuad[6][1]+fReflZMuQuad[5][1]-1.0*fReflZMuQuad[2][1]+2.0*fReflZMuQuad[1][1]-1.0*fReflZMuQuad[0][1]); 
  fRefl[18] = 0.2777777777777778*(fReflZMuQuad[7][2]-1.0*(fReflZMuQuad[5][2]+fReflZMuQuad[2][2])+fReflZMuQuad[0][2]); 
  fRefl[19] = 0.1851851851851853*(fReflZMuQuad[7][1]-1.0*fReflZMuQuad[5][1]+2.0*(fReflZMuQuad[3][1]-1.0*fReflZMuQuad[4][1])+fReflZMuQuad[2][1]-1.0*fReflZMuQuad[0][1]); 
  } 

 
}

void bc_sheath_gyrokinetic_reflectedf_upper_1x2v_ser_p2(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double q2Dm, const double *phi, const double *phiWall, const double *f, double *fRefl) 
{ 
  double vcutSq; long double xc, b, xbarVal, fac; 
  double fReflZMuQuad[8][8]; 
  

  vcutSq = 0.5*(3.16227766016838*phiWall[2]-3.16227766016838*phi[2]+2.449489742783178*phiWall[1]-2.449489742783178*phi[1]+1.414213562373095*phiWall[0]-1.414213562373095*phi[0])*q2Dm; 
  if(vcutSq <= vlowerSq) { // absorb (no reflection) 
  fRefl[0] = 0.0; 
  fRefl[1] = 0.0; 
  fRefl[2] = 0.0; 
  fRefl[3] = 0.0; 
  fRefl[4] = 0.0; 
  fRefl[5] = 0.0; 
  fRefl[6] = 0.0; 
  fRefl[7] = 0.0; 
  fRefl[8] = 0.0; 
  fRefl[9] = 0.0; 
  fRefl[10] = 0.0; 
  fRefl[11] = 0.0; 
  fRefl[12] = 0.0; 
  fRefl[13] = 0.0; 
  fRefl[14] = 0.0; 
  fRefl[15] = 0.0; 
  fRefl[16] = 0.0; 
  fRefl[17] = 0.0; 
  fRefl[18] = 0.0; 
  fRefl[19] = 0.0; 
  } else if (vcutSq > vupperSq) { // full reflection 
  fRefl[0] = f[0]; 
  fRefl[1] = f[1]; 
  fRefl[2] = f[2]; 
  fRefl[3] = f[3]; 
  fRefl[4] = f[4]; 
  fRefl[5] = f[5]; 
  fRefl[6] = f[6]; 
  fRefl[7] = f[7]; 
  fRefl[8] = f[8]; 
  fRefl[9] = f[9]; 
  fRefl[10] = f[10]; 
  fRefl[11] = f[11]; 
  fRefl[12] = f[12]; 
  fRefl[13] = f[13]; 
  fRefl[14] = f[14]; 
  fRefl[15] = f[15]; 
  fRefl[16] = f[16]; 
  fRefl[17] = f[17]; 
  fRefl[18] = f[18]; 
  fRefl[19] = f[19]; 
  } else { // partial reflection 
  xbarVal = (0.9622504486493765*(2.0*(9.0*(f[19]+f[17])-6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[6]+f[4])-3.0*f[10])-5.0*f[2])))/(4.47213595499958*(6.708203932499369*(f[15]+f[13])-5.0*(f[9]+f[7]))+3.0*(11.18033988749895*(f[3]+f[1])-15.0*f[5])-25.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if (-0.02*(4.47213595499958*(6.708203932499369*(f[15]+f[13])-5.0*(f[9]+f[7]))+3.0*(11.18033988749895*(f[3]+f[1])-15.0*f[5])-25.0*f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[0][0] = 0.0; 
  fReflZMuQuad[0][1] = 0.0; 
  fReflZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][0] = (-0.02*(4.47213595499958*(6.708203932499369*(f[15]+f[13])-5.0*(f[9]+f[7]))+3.0*(11.18033988749895*(f[3]+f[1])-15.0*f[5])-25.0*f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][1] = (-0.03333333333333333*(2.0*(9.0*(f[19]+f[17])-6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[6]+f[4])-3.0*f[10])-5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][2] = (0.1*(9.0*f[18]-6.708203932499369*(f[14]+f[12])+5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][0] = (-0.02*(4.47213595499958*(6.708203932499369*(f[15]+f[13])-5.0*(f[9]+f[7]))+3.0*(11.18033988749895*(f[3]+f[1])-15.0*f[5])-25.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][1] = (-0.03333333333333333*(2.0*(9.0*(f[19]+f[17])-6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[6]+f[4])-3.0*f[10])-5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][2] = (0.1*(9.0*f[18]-6.708203932499369*(f[14]+f[12])+5.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(45.0*f[17]+6.708203932499369*(4.0*f[16]-5.0*f[11])+6.0*(5.0*f[2]-6.708203932499369*f[6])))/(2.23606797749979*(6.708203932499369*f[13]+4.0*f[9]-5.0*f[7])+2.0*(5.0*f[0]-6.708203932499369*f[3])); 
  // if f is not realizable, no reflection from this node 
  if (0.05*(2.23606797749979*(6.708203932499369*f[13]+4.0*f[9]-5.0*f[7])+2.0*(5.0*f[0]-6.708203932499369*f[3])) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[1][0] = 0.0; 
  fReflZMuQuad[1][1] = 0.0; 
  fReflZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][0] = (0.05*(2.23606797749979*(6.708203932499369*f[13]+4.0*f[9]-5.0*f[7])+2.0*(5.0*f[0]-6.708203932499369*f[3])))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][1] = (0.01666666666666667*(45.0*f[17]+6.708203932499369*(4.0*f[16]-5.0*f[11])+6.0*(5.0*f[2]-6.708203932499369*f[6])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][2] = (-0.1*(6.708203932499369*f[14]-5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][0] = (0.05*(2.23606797749979*(6.708203932499369*f[13]+4.0*f[9]-5.0*f[7])+2.0*(5.0*f[0]-6.708203932499369*f[3])))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][1] = (0.01666666666666667*(45.0*f[17]+6.708203932499369*(4.0*f[16]-5.0*f[11])+6.0*(5.0*f[2]-6.708203932499369*f[6])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][2] = (-0.1*(6.708203932499369*f[14]-5.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(9.0*(f[19]-1.0*f[17])+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[4]-1.0*f[6])-3.0*f[10])+5.0*f[2])))/(4.47213595499958*(6.708203932499369*(f[15]-1.0*f[13])+5.0*(f[9]+f[7]))+3.0*(11.18033988749895*(f[1]-1.0*f[3])-15.0*f[5])+25.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if (0.02*(4.47213595499958*(6.708203932499369*(f[15]-1.0*f[13])+5.0*(f[9]+f[7]))+3.0*(11.18033988749895*(f[1]-1.0*f[3])-15.0*f[5])+25.0*f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[2][0] = 0.0; 
  fReflZMuQuad[2][1] = 0.0; 
  fReflZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][0] = (0.02*(4.47213595499958*(6.708203932499369*(f[15]-1.0*f[13])+5.0*(f[9]+f[7]))+3.0*(11.18033988749895*(f[1]-1.0*f[3])-15.0*f[5])+25.0*f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][1] = (0.03333333333333333*(2.0*(9.0*(f[19]-1.0*f[17])+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[4]-1.0*f[6])-3.0*f[10])+5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][2] = (-0.1*(9.0*f[18]+6.708203932499369*f[14]-1.0*(6.708203932499369*f[12]+5.0*f[8])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][0] = (0.02*(4.47213595499958*(6.708203932499369*(f[15]-1.0*f[13])+5.0*(f[9]+f[7]))+3.0*(11.18033988749895*(f[1]-1.0*f[3])-15.0*f[5])+25.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][1] = (0.03333333333333333*(2.0*(9.0*(f[19]-1.0*f[17])+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[4]-1.0*f[6])-3.0*f[10])+5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][2] = (-0.1*(9.0*f[18]+6.708203932499369*f[14]-1.0*(6.708203932499369*f[12]+5.0*f[8])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(5.0*(9.0*f[19]-6.708203932499369*f[16])+2.0*(13.41640786499874*f[11]+3.0*(5.0*f[2]-6.708203932499369*f[4]))))/(2.23606797749979*(6.708203932499369*f[15]-5.0*f[9])+2.0*(2.23606797749979*(2.0*f[7]-3.0*f[1])+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if (0.05*(2.23606797749979*(6.708203932499369*f[15]-5.0*f[9])+2.0*(2.23606797749979*(2.0*f[7]-3.0*f[1])+5.0*f[0])) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[3][0] = 0.0; 
  fReflZMuQuad[3][1] = 0.0; 
  fReflZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][0] = (0.05*(2.23606797749979*(6.708203932499369*f[15]-5.0*f[9])+2.0*(2.23606797749979*(2.0*f[7]-3.0*f[1])+5.0*f[0])))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][1] = (0.01666666666666667*(5.0*(9.0*f[19]-6.708203932499369*f[16])+2.0*(13.41640786499874*f[11]+3.0*(5.0*f[2]-6.708203932499369*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][2] = (-0.1*(6.708203932499369*f[12]-5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][0] = (0.05*(2.23606797749979*(6.708203932499369*f[15]-5.0*f[9])+2.0*(2.23606797749979*(2.0*f[7]-3.0*f[1])+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][1] = (0.01666666666666667*(5.0*(9.0*f[19]-6.708203932499369*f[16])+2.0*(13.41640786499874*f[11]+3.0*(5.0*f[2]-6.708203932499369*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][2] = (-0.1*(6.708203932499369*f[12]-5.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(5.0*(9.0*f[19]+6.708203932499369*f[16])-2.0*(13.41640786499874*f[11]+3.0*(6.708203932499369*f[4]+5.0*f[2]))))/(2.23606797749979*(6.708203932499369*f[15]+5.0*f[9])-2.0*(2.23606797749979*(2.0*f[7]+3.0*f[1])+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if (-0.05*(2.23606797749979*(6.708203932499369*f[15]+5.0*f[9])-2.0*(2.23606797749979*(2.0*f[7]+3.0*f[1])+5.0*f[0])) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[4][0] = 0.0; 
  fReflZMuQuad[4][1] = 0.0; 
  fReflZMuQuad[4][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][0] = (-0.05*(2.23606797749979*(6.708203932499369*f[15]+5.0*f[9])-2.0*(2.23606797749979*(2.0*f[7]+3.0*f[1])+5.0*f[0])))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][1] = (-0.01666666666666667*(5.0*(9.0*f[19]+6.708203932499369*f[16])-2.0*(13.41640786499874*f[11]+3.0*(6.708203932499369*f[4]+5.0*f[2]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][2] = (0.1*(6.708203932499369*f[12]+5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][0] = (-0.05*(2.23606797749979*(6.708203932499369*f[15]+5.0*f[9])-2.0*(2.23606797749979*(2.0*f[7]+3.0*f[1])+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][1] = (-0.01666666666666667*(5.0*(9.0*f[19]+6.708203932499369*f[16])-2.0*(13.41640786499874*f[11]+3.0*(6.708203932499369*f[4]+5.0*f[2]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][2] = (0.1*(6.708203932499369*f[12]+5.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(9.0*f[19]-1.0*(9.0*f[17]+6.708203932499369*(f[16]+f[11])))+3.0*(3.0*(3.0*f[10]+2.23606797749979*(f[4]-1.0*f[6]))-5.0*f[2])))/(4.47213595499958*(6.708203932499369*f[15]-1.0*(6.708203932499369*f[13]+5.0*(f[9]+f[7])))+3.0*(15.0*f[5]+11.18033988749895*(f[1]-1.0*f[3]))-25.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if (-0.02*(4.47213595499958*(6.708203932499369*f[15]-1.0*(6.708203932499369*f[13]+5.0*(f[9]+f[7])))+3.0*(15.0*f[5]+11.18033988749895*(f[1]-1.0*f[3]))-25.0*f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[5][0] = 0.0; 
  fReflZMuQuad[5][1] = 0.0; 
  fReflZMuQuad[5][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][0] = (-0.02*(4.47213595499958*(6.708203932499369*f[15]-1.0*(6.708203932499369*f[13]+5.0*(f[9]+f[7])))+3.0*(15.0*f[5]+11.18033988749895*(f[1]-1.0*f[3]))-25.0*f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][1] = (-0.03333333333333333*(2.0*(9.0*f[19]-1.0*(9.0*f[17]+6.708203932499369*(f[16]+f[11])))+3.0*(3.0*(3.0*f[10]+2.23606797749979*(f[4]-1.0*f[6]))-5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][2] = (-0.1*(9.0*f[18]+6.708203932499369*(f[12]-1.0*f[14])-5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][0] = (-0.02*(4.47213595499958*(6.708203932499369*f[15]-1.0*(6.708203932499369*f[13]+5.0*(f[9]+f[7])))+3.0*(15.0*f[5]+11.18033988749895*(f[1]-1.0*f[3]))-25.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][1] = (-0.03333333333333333*(2.0*(9.0*f[19]-1.0*(9.0*f[17]+6.708203932499369*(f[16]+f[11])))+3.0*(3.0*(3.0*f[10]+2.23606797749979*(f[4]-1.0*f[6]))-5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][2] = (-0.1*(9.0*f[18]+6.708203932499369*(f[12]-1.0*f[14])-5.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(45.0*f[17]+6.708203932499369*(5.0*f[11]-4.0*f[16])-6.0*(6.708203932499369*f[6]+5.0*f[2])))/(2.23606797749979*(6.708203932499369*f[13]-4.0*f[9]+5.0*f[7])-2.0*(6.708203932499369*f[3]+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if (-0.05*(2.23606797749979*(6.708203932499369*f[13]-4.0*f[9]+5.0*f[7])-2.0*(6.708203932499369*f[3]+5.0*f[0])) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[6][0] = 0.0; 
  fReflZMuQuad[6][1] = 0.0; 
  fReflZMuQuad[6][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][0] = (-0.05*(2.23606797749979*(6.708203932499369*f[13]-4.0*f[9]+5.0*f[7])-2.0*(6.708203932499369*f[3]+5.0*f[0])))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][1] = (-0.01666666666666667*(45.0*f[17]+6.708203932499369*(5.0*f[11]-4.0*f[16])-6.0*(6.708203932499369*f[6]+5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][2] = (0.1*(6.708203932499369*f[14]+5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][0] = (-0.05*(2.23606797749979*(6.708203932499369*f[13]-4.0*f[9]+5.0*f[7])-2.0*(6.708203932499369*f[3]+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][1] = (-0.01666666666666667*(45.0*f[17]+6.708203932499369*(5.0*f[11]-4.0*f[16])-6.0*(6.708203932499369*f[6]+5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][2] = (0.1*(6.708203932499369*f[14]+5.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(9.0*(f[19]+f[17])+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(3.0*f[10]+2.23606797749979*(f[6]+f[4]))+5.0*f[2])))/(4.47213595499958*(6.708203932499369*(f[15]+f[13])+5.0*(f[9]+f[7]))+3.0*(15.0*f[5]+11.18033988749895*(f[3]+f[1]))+25.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if (0.02*(4.47213595499958*(6.708203932499369*(f[15]+f[13])+5.0*(f[9]+f[7]))+3.0*(15.0*f[5]+11.18033988749895*(f[3]+f[1]))+25.0*f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[7][0] = 0.0; 
  fReflZMuQuad[7][1] = 0.0; 
  fReflZMuQuad[7][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][0] = (0.02*(4.47213595499958*(6.708203932499369*(f[15]+f[13])+5.0*(f[9]+f[7]))+3.0*(15.0*f[5]+11.18033988749895*(f[3]+f[1]))+25.0*f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][1] = (0.03333333333333333*(2.0*(9.0*(f[19]+f[17])+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(3.0*f[10]+2.23606797749979*(f[6]+f[4]))+5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][2] = (0.1*(9.0*f[18]+6.708203932499369*(f[14]+f[12])+5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][0] = (0.02*(4.47213595499958*(6.708203932499369*(f[15]+f[13])+5.0*(f[9]+f[7]))+3.0*(15.0*f[5]+11.18033988749895*(f[3]+f[1]))+25.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][1] = (0.03333333333333333*(2.0*(9.0*(f[19]+f[17])+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(3.0*f[10]+2.23606797749979*(f[6]+f[4]))+5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][2] = (0.1*(9.0*f[18]+6.708203932499369*(f[14]+f[12])+5.0*f[8]))*fac; 
   } 
  } 
  fRefl[0] = 0.05555555555555555*(fReflZMuQuad[7][0]+8.0*fReflZMuQuad[6][0]+fReflZMuQuad[5][0]+8.0*(fReflZMuQuad[4][0]+fReflZMuQuad[3][0])+fReflZMuQuad[2][0]+8.0*fReflZMuQuad[1][0]+fReflZMuQuad[0][0]); 
  fRefl[1] = 4.469298760204439e-16*(4.63256860547201e+14*(fReflZMuQuad[7][0]-1.0*fReflZMuQuad[5][0])+7.4121097687552e+14*(fReflZMuQuad[4][0]-1.0*fReflZMuQuad[3][0])+4.63256860547201e+14*(fReflZMuQuad[2][0]-1.0*fReflZMuQuad[0][0])); 
  fRefl[2] = 0.05555555555555555*(fReflZMuQuad[7][1]+8.0*fReflZMuQuad[6][1]+fReflZMuQuad[5][1]+8.0*(fReflZMuQuad[4][1]+fReflZMuQuad[3][1])+fReflZMuQuad[2][1]+8.0*fReflZMuQuad[1][1]+fReflZMuQuad[0][1]); 
  fRefl[3] = 4.469298760204439e-16*(4.63256860547201e+14*fReflZMuQuad[7][0]+7.4121097687552e+14*fReflZMuQuad[6][0]+4.63256860547201e+14*fReflZMuQuad[5][0]-1.0*(4.63256860547201e+14*fReflZMuQuad[2][0]+7.4121097687552e+14*fReflZMuQuad[1][0]+4.63256860547201e+14*fReflZMuQuad[0][0])); 
  fRefl[4] = 4.469298760204439e-16*(4.63256860547201e+14*(fReflZMuQuad[7][1]-1.0*fReflZMuQuad[5][1])+7.4121097687552e+14*(fReflZMuQuad[4][1]-1.0*fReflZMuQuad[3][1])+4.63256860547201e+14*(fReflZMuQuad[2][1]-1.0*fReflZMuQuad[0][1])); 
  fRefl[5] = 0.2777777777777778*(fReflZMuQuad[7][0]-1.0*(fReflZMuQuad[5][0]+fReflZMuQuad[2][0])+fReflZMuQuad[0][0]); 
  fRefl[6] = 4.469298760204439e-16*(4.63256860547201e+14*fReflZMuQuad[7][1]+7.4121097687552e+14*fReflZMuQuad[6][1]+4.63256860547201e+14*fReflZMuQuad[5][1]-1.0*(4.63256860547201e+14*fReflZMuQuad[2][1]+7.4121097687552e+14*fReflZMuQuad[1][1]+4.63256860547201e+14*fReflZMuQuad[0][1])); 
  fRefl[7] = 0.2484519974999762*(fReflZMuQuad[7][0]-2.0*fReflZMuQuad[6][0]+fReflZMuQuad[5][0]+fReflZMuQuad[2][0]-2.0*fReflZMuQuad[1][0]+fReflZMuQuad[0][0]); 
  fRefl[8] = 0.05555555555555555*(fReflZMuQuad[7][2]+8.0*fReflZMuQuad[6][2]+fReflZMuQuad[5][2]+8.0*(fReflZMuQuad[4][2]+fReflZMuQuad[3][2])+fReflZMuQuad[2][2]+8.0*fReflZMuQuad[1][2]+fReflZMuQuad[0][2]); 
  fRefl[9] = 0.2484519974999762*(fReflZMuQuad[7][0]+fReflZMuQuad[5][0]-2.0*(fReflZMuQuad[4][0]+fReflZMuQuad[3][0])+fReflZMuQuad[2][0]+fReflZMuQuad[0][0]); 
  fRefl[10] = 0.2777777777777778*(fReflZMuQuad[7][1]-1.0*(fReflZMuQuad[5][1]+fReflZMuQuad[2][1])+fReflZMuQuad[0][1]); 
  fRefl[11] = 0.2484519974999762*(fReflZMuQuad[7][1]-2.0*fReflZMuQuad[6][1]+fReflZMuQuad[5][1]+fReflZMuQuad[2][1]-2.0*fReflZMuQuad[1][1]+fReflZMuQuad[0][1]); 
  fRefl[12] = 4.46929876020444e-16*(4.63256860547201e+14*(fReflZMuQuad[7][2]-1.0*fReflZMuQuad[5][2])+7.4121097687552e+14*fReflZMuQuad[4][2]+4.63256860547201e+14*(fReflZMuQuad[2][2]-1.0*fReflZMuQuad[0][2])); 
  fRefl[13] = 0.1851851851851852*(fReflZMuQuad[7][0]-2.0*fReflZMuQuad[6][0]+fReflZMuQuad[5][0]-1.0*fReflZMuQuad[2][0]+2.0*fReflZMuQuad[1][0]-1.0*fReflZMuQuad[0][0]); 
  fRefl[14] = 4.46929876020444e-16*(4.63256860547201e+14*fReflZMuQuad[7][2]+7.4121097687552e+14*fReflZMuQuad[6][2]+4.63256860547201e+14*fReflZMuQuad[5][2]-1.0*(4.63256860547201e+14*fReflZMuQuad[2][2]+7.4121097687552e+14*fReflZMuQuad[1][2]+4.63256860547201e+14*fReflZMuQuad[0][2])); 
  fRefl[15] = 0.1851851851851852*(fReflZMuQuad[7][0]-1.0*fReflZMuQuad[5][0]+2.0*(fReflZMuQuad[3][0]-1.0*fReflZMuQuad[4][0])+fReflZMuQuad[2][0]-1.0*fReflZMuQuad[0][0]); 
  fRefl[16] = 0.2484519974999762*(fReflZMuQuad[7][1]+fReflZMuQuad[5][1]-2.0*(fReflZMuQuad[4][1]+fReflZMuQuad[3][1])+fReflZMuQuad[2][1]+fReflZMuQuad[0][1]); 
  fRefl[17] = 0.1851851851851853*(fReflZMuQuad[7][1]-2.0*fReflZMuQuad[6][1]+fReflZMuQuad[5][1]-1.0*fReflZMuQuad[2][1]+2.0*fReflZMuQuad[1][1]-1.0*fReflZMuQuad[0][1]); 
  fRefl[18] = 0.2777777777777778*(fReflZMuQuad[7][2]-1.0*(fReflZMuQuad[5][2]+fReflZMuQuad[2][2])+fReflZMuQuad[0][2]); 
  fRefl[19] = 0.1851851851851853*(fReflZMuQuad[7][1]-1.0*fReflZMuQuad[5][1]+2.0*(fReflZMuQuad[3][1]-1.0*fReflZMuQuad[4][1])+fReflZMuQuad[2][1]-1.0*fReflZMuQuad[0][1]); 
  } 

 
}
