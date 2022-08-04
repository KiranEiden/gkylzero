#include <gkyl_bc_sheath_gyrokinetic.h> 


void bc_sheath_gyrokinetic_reflectedf_lower_1x2v_ser_p1(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double q2Dm, const double *phi, const double *phiWall, const double *f, double *fRefl) 
{ 
  double vcutSq; long double xc, b, xbarVal, fac; 
  double fReflZMuQuad[4][6]; 
  

  vcutSq = -0.5*(2.449489742783178*phiWall[1]-2.449489742783178*phi[1]-1.414213562373095*phiWall[0]+1.414213562373095*phi[0])*q2Dm; 
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
  } else { // partial reflection 
  xbarVal = (0.5773502691896258*(f[7]-1.732050807568877*(f[6]+f[4])+3.0*f[2]))/(f[5]-1.732050807568877*(f[3]+f[1])+3.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if (0.1666666666666667*(f[5]-1.732050807568877*(f[3]+f[1])+3.0*f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[0][0] = 0.0; 
  fReflZMuQuad[0][1] = 0.0; 
  fReflZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][0] = (0.1666666666666667*(f[5]-1.732050807568877*(f[3]+f[1])+3.0*f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][1] = (0.1666666666666667*(f[7]-1.732050807568877*(f[6]+f[4])+3.0*f[2]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][2] = (0.03333333333333333*(5.0*f[11]-8.660254037844387*(f[10]+f[9])+15.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][0] = (0.1666666666666667*(f[5]-1.732050807568877*(f[3]+f[1])+3.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][1] = (0.1666666666666667*(f[7]-1.732050807568877*(f[6]+f[4])+3.0*f[2]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][2] = (0.03333333333333333*(5.0*f[11]-8.660254037844387*(f[10]+f[9])+15.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[7]+1.732050807568877*f[6]-1.0*(1.732050807568877*f[4]+3.0*f[2])))/(f[5]+1.732050807568877*f[3]-1.0*(1.732050807568877*f[1]+3.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if (-0.1666666666666667*(f[5]+1.732050807568877*f[3]-1.0*(1.732050807568877*f[1]+3.0*f[0])) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[1][0] = 0.0; 
  fReflZMuQuad[1][1] = 0.0; 
  fReflZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][0] = (-0.1666666666666667*(f[5]+1.732050807568877*f[3]-1.0*(1.732050807568877*f[1]+3.0*f[0])))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][1] = (-0.1666666666666667*(f[7]+1.732050807568877*f[6]-1.0*(1.732050807568877*f[4]+3.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][2] = (-0.03333333333333333*(5.0*f[11]+8.660254037844387*f[10]-1.0*(8.660254037844387*f[9]+15.0*f[8])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][0] = (-0.1666666666666667*(f[5]+1.732050807568877*f[3]-1.0*(1.732050807568877*f[1]+3.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][1] = (-0.1666666666666667*(f[7]+1.732050807568877*f[6]-1.0*(1.732050807568877*f[4]+3.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][2] = (-0.03333333333333333*(5.0*f[11]+8.660254037844387*f[10]-1.0*(8.660254037844387*f[9]+15.0*f[8])))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[7]+1.732050807568877*(f[4]-1.0*f[6])-3.0*f[2]))/(f[5]+1.732050807568877*(f[1]-1.0*f[3])-3.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if (-0.1666666666666667*(f[5]+1.732050807568877*(f[1]-1.0*f[3])-3.0*f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[2][0] = 0.0; 
  fReflZMuQuad[2][1] = 0.0; 
  fReflZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][0] = (-0.1666666666666667*(f[5]+1.732050807568877*(f[1]-1.0*f[3])-3.0*f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][1] = (-0.1666666666666667*(f[7]+1.732050807568877*(f[4]-1.0*f[6])-3.0*f[2]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][2] = (-0.03333333333333333*(5.0*f[11]+8.660254037844387*(f[9]-1.0*f[10])-15.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][0] = (-0.1666666666666667*(f[5]+1.732050807568877*(f[1]-1.0*f[3])-3.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][1] = (-0.1666666666666667*(f[7]+1.732050807568877*(f[4]-1.0*f[6])-3.0*f[2]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][2] = (-0.03333333333333333*(5.0*f[11]+8.660254037844387*(f[9]-1.0*f[10])-15.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[7]+1.732050807568877*(f[6]+f[4])+3.0*f[2]))/(f[5]+1.732050807568877*(f[3]+f[1])+3.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if (0.1666666666666667*(f[5]+1.732050807568877*(f[3]+f[1])+3.0*f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[3][0] = 0.0; 
  fReflZMuQuad[3][1] = 0.0; 
  fReflZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][0] = (0.1666666666666667*(f[5]+1.732050807568877*(f[3]+f[1])+3.0*f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][1] = (0.1666666666666667*(f[7]+1.732050807568877*(f[6]+f[4])+3.0*f[2]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][2] = (0.03333333333333333*(5.0*f[11]+8.660254037844387*(f[10]+f[9])+15.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][0] = (0.1666666666666667*(f[5]+1.732050807568877*(f[3]+f[1])+3.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][1] = (0.1666666666666667*(f[7]+1.732050807568877*(f[6]+f[4])+3.0*f[2]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][2] = (0.03333333333333333*(5.0*f[11]+8.660254037844387*(f[10]+f[9])+15.0*f[8]))*fac; 
   } 
  } 
  fRefl[0] = 0.5*(fReflZMuQuad[3][0]+fReflZMuQuad[2][0]+fReflZMuQuad[1][0]+fReflZMuQuad[0][0]); 
  fRefl[1] = 0.8660254037844385*(fReflZMuQuad[3][0]-1.0*fReflZMuQuad[2][0]+fReflZMuQuad[1][0]-1.0*fReflZMuQuad[0][0]); 
  fRefl[2] = 0.5*(fReflZMuQuad[3][1]+fReflZMuQuad[2][1]+fReflZMuQuad[1][1]+fReflZMuQuad[0][1]); 
  fRefl[3] = 0.8660254037844385*(fReflZMuQuad[3][0]+fReflZMuQuad[2][0]-1.0*(fReflZMuQuad[1][0]+fReflZMuQuad[0][0])); 
  fRefl[4] = 0.8660254037844385*(fReflZMuQuad[3][1]-1.0*fReflZMuQuad[2][1]+fReflZMuQuad[1][1]-1.0*fReflZMuQuad[0][1]); 
  fRefl[5] = 1.5*(fReflZMuQuad[3][0]-1.0*(fReflZMuQuad[2][0]+fReflZMuQuad[1][0])+fReflZMuQuad[0][0]); 
  fRefl[6] = 0.8660254037844385*(fReflZMuQuad[3][1]+fReflZMuQuad[2][1]-1.0*(fReflZMuQuad[1][1]+fReflZMuQuad[0][1])); 
  fRefl[7] = 1.5*(fReflZMuQuad[3][1]-1.0*(fReflZMuQuad[2][1]+fReflZMuQuad[1][1])+fReflZMuQuad[0][1]); 
  fRefl[8] = 0.5*(fReflZMuQuad[3][2]+fReflZMuQuad[2][2]+fReflZMuQuad[1][2]+fReflZMuQuad[0][2]); 
  fRefl[9] = 0.8660254037844384*(fReflZMuQuad[3][2]-1.0*fReflZMuQuad[2][2]+fReflZMuQuad[1][2]-1.0*fReflZMuQuad[0][2]); 
  fRefl[10] = 0.8660254037844384*(fReflZMuQuad[3][2]+fReflZMuQuad[2][2]-1.0*(fReflZMuQuad[1][2]+fReflZMuQuad[0][2])); 
  fRefl[11] = 1.5*(fReflZMuQuad[3][2]-1.0*(fReflZMuQuad[2][2]+fReflZMuQuad[1][2])+fReflZMuQuad[0][2]); 
  } 

 
}

void bc_sheath_gyrokinetic_reflectedf_upper_1x2v_ser_p1(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double q2Dm, const double *phi, const double *phiWall, const double *f, double *fRefl) 
{ 
  double vcutSq; long double xc, b, xbarVal, fac; 
  double fReflZMuQuad[4][6]; 
  

  vcutSq = 0.5*(2.449489742783178*phiWall[1]-2.449489742783178*phi[1]+1.414213562373095*phiWall[0]-1.414213562373095*phi[0])*q2Dm; 
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
  } else { // partial reflection 
  xbarVal = (0.5773502691896258*(f[7]-1.732050807568877*(f[6]+f[4])+3.0*f[2]))/(f[5]-1.732050807568877*(f[3]+f[1])+3.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if (0.1666666666666667*(f[5]-1.732050807568877*(f[3]+f[1])+3.0*f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[0][0] = 0.0; 
  fReflZMuQuad[0][1] = 0.0; 
  fReflZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][0] = (0.1666666666666667*(f[5]-1.732050807568877*(f[3]+f[1])+3.0*f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][1] = (0.1666666666666667*(f[7]-1.732050807568877*(f[6]+f[4])+3.0*f[2]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][2] = (0.03333333333333333*(5.0*f[11]-8.660254037844387*(f[10]+f[9])+15.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][0] = (0.1666666666666667*(f[5]-1.732050807568877*(f[3]+f[1])+3.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][1] = (0.1666666666666667*(f[7]-1.732050807568877*(f[6]+f[4])+3.0*f[2]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][2] = (0.03333333333333333*(5.0*f[11]-8.660254037844387*(f[10]+f[9])+15.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[7]+1.732050807568877*f[6]-1.0*(1.732050807568877*f[4]+3.0*f[2])))/(f[5]+1.732050807568877*f[3]-1.0*(1.732050807568877*f[1]+3.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if (-0.1666666666666667*(f[5]+1.732050807568877*f[3]-1.0*(1.732050807568877*f[1]+3.0*f[0])) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[1][0] = 0.0; 
  fReflZMuQuad[1][1] = 0.0; 
  fReflZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][0] = (-0.1666666666666667*(f[5]+1.732050807568877*f[3]-1.0*(1.732050807568877*f[1]+3.0*f[0])))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][1] = (-0.1666666666666667*(f[7]+1.732050807568877*f[6]-1.0*(1.732050807568877*f[4]+3.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][2] = (-0.03333333333333333*(5.0*f[11]+8.660254037844387*f[10]-1.0*(8.660254037844387*f[9]+15.0*f[8])))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][0] = (-0.1666666666666667*(f[5]+1.732050807568877*f[3]-1.0*(1.732050807568877*f[1]+3.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][1] = (-0.1666666666666667*(f[7]+1.732050807568877*f[6]-1.0*(1.732050807568877*f[4]+3.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][2] = (-0.03333333333333333*(5.0*f[11]+8.660254037844387*f[10]-1.0*(8.660254037844387*f[9]+15.0*f[8])))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[7]+1.732050807568877*(f[4]-1.0*f[6])-3.0*f[2]))/(f[5]+1.732050807568877*(f[1]-1.0*f[3])-3.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if (-0.1666666666666667*(f[5]+1.732050807568877*(f[1]-1.0*f[3])-3.0*f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[2][0] = 0.0; 
  fReflZMuQuad[2][1] = 0.0; 
  fReflZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][0] = (-0.1666666666666667*(f[5]+1.732050807568877*(f[1]-1.0*f[3])-3.0*f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][1] = (-0.1666666666666667*(f[7]+1.732050807568877*(f[4]-1.0*f[6])-3.0*f[2]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][2] = (-0.03333333333333333*(5.0*f[11]+8.660254037844387*(f[9]-1.0*f[10])-15.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][0] = (-0.1666666666666667*(f[5]+1.732050807568877*(f[1]-1.0*f[3])-3.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][1] = (-0.1666666666666667*(f[7]+1.732050807568877*(f[4]-1.0*f[6])-3.0*f[2]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][2] = (-0.03333333333333333*(5.0*f[11]+8.660254037844387*(f[9]-1.0*f[10])-15.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.5773502691896258*(f[7]+1.732050807568877*(f[6]+f[4])+3.0*f[2]))/(f[5]+1.732050807568877*(f[3]+f[1])+3.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if (0.1666666666666667*(f[5]+1.732050807568877*(f[3]+f[1])+3.0*f[0]) <= 0. || fabsl(xbarVal)>=.95) { 
  fReflZMuQuad[3][0] = 0.0; 
  fReflZMuQuad[3][1] = 0.0; 
  fReflZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : fabsl(b)<1e-10? (1.+xc)/2. : (expl(b*xc)-expl(-b))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][0] = (0.1666666666666667*(f[5]+1.732050807568877*(f[3]+f[1])+3.0*f[0]))*fac; 
    fac = (b>500 || fabsl(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*expl(b*xc)+(b+1)*expl(-b))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][1] = (0.1666666666666667*(f[7]+1.732050807568877*(f[6]+f[4])+3.0*f[2]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3-(2*(b*b+3*(b+1))*expl(-b))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][2] = (0.03333333333333333*(5.0*f[11]+8.660254037844387*(f[10]+f[9])+15.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : fabsl(b)<1e-10? (1.-xc)/2. : (expl(b)-expl(b*xc))/(2.*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][0] = (0.1666666666666667*(f[5]+1.732050807568877*(f[3]+f[1])+3.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || fabsl(b)<1e-8)? 0. : ((b-1)*expl(b)-(b*xc-1)*expl(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][1] = (0.1666666666666667*(f[7]+1.732050807568877*(f[6]+f[4])+3.0*f[2]))*fac; 
    fac = ((2*(b*b+3*(1-b))*expl(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*expl(b*xc))/3)/(-4*b*cosh(b) + 4/3*(3 + b*b)*sinh(b)); 
    if(isnan(fac) || isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][2] = (0.03333333333333333*(5.0*f[11]+8.660254037844387*(f[10]+f[9])+15.0*f[8]))*fac; 
   } 
  } 
  fRefl[0] = 0.5*(fReflZMuQuad[3][0]+fReflZMuQuad[2][0]+fReflZMuQuad[1][0]+fReflZMuQuad[0][0]); 
  fRefl[1] = 0.8660254037844385*(fReflZMuQuad[3][0]-1.0*fReflZMuQuad[2][0]+fReflZMuQuad[1][0]-1.0*fReflZMuQuad[0][0]); 
  fRefl[2] = 0.5*(fReflZMuQuad[3][1]+fReflZMuQuad[2][1]+fReflZMuQuad[1][1]+fReflZMuQuad[0][1]); 
  fRefl[3] = 0.8660254037844385*(fReflZMuQuad[3][0]+fReflZMuQuad[2][0]-1.0*(fReflZMuQuad[1][0]+fReflZMuQuad[0][0])); 
  fRefl[4] = 0.8660254037844385*(fReflZMuQuad[3][1]-1.0*fReflZMuQuad[2][1]+fReflZMuQuad[1][1]-1.0*fReflZMuQuad[0][1]); 
  fRefl[5] = 1.5*(fReflZMuQuad[3][0]-1.0*(fReflZMuQuad[2][0]+fReflZMuQuad[1][0])+fReflZMuQuad[0][0]); 
  fRefl[6] = 0.8660254037844385*(fReflZMuQuad[3][1]+fReflZMuQuad[2][1]-1.0*(fReflZMuQuad[1][1]+fReflZMuQuad[0][1])); 
  fRefl[7] = 1.5*(fReflZMuQuad[3][1]-1.0*(fReflZMuQuad[2][1]+fReflZMuQuad[1][1])+fReflZMuQuad[0][1]); 
  fRefl[8] = 0.5*(fReflZMuQuad[3][2]+fReflZMuQuad[2][2]+fReflZMuQuad[1][2]+fReflZMuQuad[0][2]); 
  fRefl[9] = 0.8660254037844384*(fReflZMuQuad[3][2]-1.0*fReflZMuQuad[2][2]+fReflZMuQuad[1][2]-1.0*fReflZMuQuad[0][2]); 
  fRefl[10] = 0.8660254037844384*(fReflZMuQuad[3][2]+fReflZMuQuad[2][2]-1.0*(fReflZMuQuad[1][2]+fReflZMuQuad[0][2])); 
  fRefl[11] = 1.5*(fReflZMuQuad[3][2]-1.0*(fReflZMuQuad[2][2]+fReflZMuQuad[1][2])+fReflZMuQuad[0][2]); 
  } 

 
}
