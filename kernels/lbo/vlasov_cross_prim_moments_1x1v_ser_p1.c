#include <gkyl_prim_lbo_vlasov_kernels.h> 
 
GKYL_CU_DH void vlasov_cross_prim_moments_1x1v_ser_p1(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *greene, const double m_self, const double *u_self, const double *vtsq_self, const double m_other, const double *u_other, const double *vtsq_other, const double *moms, const double *boundary_corrections) 
{ 
  // greene:               Greene's factor. 
  // m_:                   mass. 
  // moms:                 moments of the distribution function. 
  // u,vtSq:               self primitive moments: mean flow velocity and thermal speed squared. 
  // boundary_corrections: corrections to momentum and energy conservation due to finite velocity space. 
  // uCross,vtSqCross:     cross primitive moments: mean flow velocity and thermal speed squared. 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if (-0.5*(2.449489742783178*moms[1]-1.414213562373095*moms[0]) < 0) cellAvg = true; 
  if (0.5*(2.449489742783178*moms[1]+1.414213562373095*moms[0]) < 0) cellAvg = true; 
 
  double m0r[2] = {0.0}; 
  double m1r[2] = {0.0}; 
  double m2r[2] = {0.0}; 
  double cMr[2] = {0.0}; 
  double cEr[2] = {0.0}; 
  if (cellAvg) { 
    m0r[0] = moms[0]; 
    m0r[1] = 0.0; 
    m1r[0] = moms[2]; 
    m1r[1] = 0.0; 
    cMr[0] = boundary_corrections[0]; 
    cMr[1] = 0.0; 
    m2r[0] = moms[4]; 
    m2r[1] = 0.0; 
    cEr[0] = boundary_corrections[2]; 
    cEr[1] = 0.0; 
  } else { 
    m0r[0] = moms[0]; 
    m0r[1] = moms[1]; 
    m1r[0] = moms[2]; 
    m1r[1] = moms[3]; 
    m2r[0] = moms[4]; 
    m2r[1] = moms[5]; 
    cMr[0] = boundary_corrections[0]; 
    cMr[1] = boundary_corrections[1]; 
    cEr[0] = boundary_corrections[2]; 
    cEr[1] = boundary_corrections[3]; 
  } 
 
  double momRHS[2] = {0.0}; 
 
  // ... Block from weak multiply of m_self, nu, M0 and uCrossX ... // 
  gkyl_mat_set(A,0,0,1.414213562373095*m0r[0]); 
  gkyl_mat_set(A,0,1,1.414213562373095*m0r[1]); 
  gkyl_mat_set(A,1,0,1.414213562373095*m0r[1]); 
  gkyl_mat_set(A,1,1,1.414213562373095*m0r[0]); 
 
  // ... Block from correction to momentum conservation (self) ... // 
  gkyl_mat_set(A,0,2,-1.414213562373095*cMr[0]); 
  gkyl_mat_set(A,0,3,-1.414213562373095*cMr[1]); 
  gkyl_mat_set(A,1,2,-1.414213562373095*cMr[1]); 
  gkyl_mat_set(A,1,3,-1.414213562373095*cMr[0]); 
 
  // ... Block from weak multiply of m_self, nu, m1X and uCrossX ... // 
  gkyl_mat_set(A,2,0,(-0.5*m0r[1]*u_self[1])-0.5*m0r[1]*u_other[1]-0.5*m0r[0]*u_self[0]-0.5*m0r[0]*u_other[0]+1.414213562373095*m1r[0]); 
  gkyl_mat_set(A,2,1,(-0.5*m0r[0]*u_self[1])-0.5*m0r[0]*u_other[1]+1.414213562373095*m1r[1]-0.5*u_self[0]*m0r[1]-0.5*u_other[0]*m0r[1]); 
  gkyl_mat_set(A,3,0,(-0.5*m0r[0]*u_self[1])-0.5*m0r[0]*u_other[1]+1.414213562373095*m1r[1]-0.5*u_self[0]*m0r[1]-0.5*u_other[0]*m0r[1]); 
  gkyl_mat_set(A,3,1,(-0.9*m0r[1]*u_self[1])-0.9*m0r[1]*u_other[1]-0.5*m0r[0]*u_self[0]-0.5*m0r[0]*u_other[0]+1.414213562373095*m1r[0]); 
 
  momRHS[0] += -0.7071067811865475*(greene[1]*u_self[1]-1.0*greene[1]*u_other[1]+greene[0]*u_self[0]-1.0*greene[0]*u_other[0]-2.828427124746191*m1r[0]); 
  momRHS[1] += -0.7071067811865475*(greene[0]*u_self[1]-1.0*greene[0]*u_other[1]-2.828427124746191*m1r[1]+(u_self[0]-1.0*u_other[0])*greene[1]); 
 
  double ucMSelf[2] = {0.0}; 
  double ucMOther[2] = {0.0}; 
  for (int vd=0; vd<1; vd++) 
  { 
    int a0 = 2*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    ucMSelf[0] += 0.7071067811865475*cMr[a0+1]*u_self[a0+1]+0.7071067811865475*cMr[a0]*u_self[a0]; 
    ucMSelf[1] += 0.7071067811865475*cMr[a0]*u_self[a0+1]+0.7071067811865475*u_self[a0]*cMr[a0+1]; 
    ucMOther[0] += 0.7071067811865475*cMr[a0+1]*u_other[a0+1]+0.7071067811865475*cMr[a0]*u_other[a0]; 
    ucMOther[1] += 0.7071067811865475*cMr[a0]*u_other[a0+1]+0.7071067811865475*u_other[a0]*cMr[a0+1]; 
  } 
 
  // ... Block from correction to (self) 2nd moment of collision operator ... // 
  gkyl_mat_set(A,2,2,0.7071067811865475*ucMSelf[0]+0.7071067811865475*ucMOther[0]+1.414213562373095*m0r[0]-1.414213562373095*cEr[0]); 
  gkyl_mat_set(A,2,3,0.7071067811865475*ucMSelf[1]+0.7071067811865475*ucMOther[1]+1.414213562373095*m0r[1]-1.414213562373095*cEr[1]); 
  gkyl_mat_set(A,3,2,0.7071067811865475*ucMSelf[1]+0.7071067811865475*ucMOther[1]+1.414213562373095*m0r[1]-1.414213562373095*cEr[1]); 
  gkyl_mat_set(A,3,3,0.7071067811865475*ucMSelf[0]+0.7071067811865475*ucMOther[0]+1.414213562373095*m0r[0]-1.414213562373095*cEr[0]); 
 
  double uM1Self[2] = {0.0}; 
  double uM1Other[2] = {0.0}; 
  double uSumSq[2] = {0.0}; 
  for (int vd=0; vd<1; vd++) 
  { 
    int a0 = 2*vd; 
    // Dot product terms in energy equation RHS. 
    uM1Self[0] += 0.7071067811865475*m1r[1]*u_self[a0+1]+0.7071067811865475*m1r[0]*u_self[a0]; 
    uM1Self[1] += 0.7071067811865475*m1r[0]*u_self[a0+1]+0.7071067811865475*m1r[1]*u_self[a0]; 
    uM1Other[0] += 0.7071067811865475*m1r[1]*u_other[a0+1]+0.7071067811865475*m1r[0]*u_other[a0]; 
    uM1Other[1] += 0.7071067811865475*m1r[0]*u_other[a0+1]+0.7071067811865475*m1r[1]*u_other[a0]; 
  const double u_self0R2 = pow(u_self[a0],2);
  const double u_self1R2 = pow(u_self[a0+1],2);
  const double u_other0R2 = pow(u_other[a0],2);
  const double u_other1R2 = pow(u_other[a0+1],2);

  uSumSq[0] += 0.7071067811865475*u_self1R2-1.414213562373095*u_other[a0+1]*u_self[a0+1]+0.7071067811865475*u_other1R2+0.7071067811865475*u_self0R2-1.414213562373095*u_other[a0]*u_self[a0]+0.7071067811865475*u_other0R2; 
  uSumSq[1] += 1.414213562373095*u_self[a0]*u_self[a0+1]-1.414213562373095*u_other[a0]*u_self[a0+1]-1.414213562373095*u_self[a0]*u_other[a0+1]+1.414213562373095*u_other[a0]*u_other[a0+1]; 
  } 
 
  double enRHS[2] = {0.0}; 
  enRHS[0] = (-(2.0*greene[1]*vtsq_self[1]*m_self)/(2.828427124746191*m_self+2.828427124746191*m_other))-(1.0*greene[1]*uSumSq[1]*m_self)/(2.828427124746191*m_self+2.828427124746191*m_other)-(2.0*greene[0]*vtsq_self[0]*m_self)/(2.828427124746191*m_self+2.828427124746191*m_other)-(1.0*greene[0]*uSumSq[0]*m_self)/(2.828427124746191*m_self+2.828427124746191*m_other)+(2.0*greene[1]*vtsq_other[1]*m_other)/(2.828427124746191*m_self+2.828427124746191*m_other)+(greene[1]*uSumSq[1]*m_other)/(2.828427124746191*m_self+2.828427124746191*m_other)+(2.0*greene[0]*vtsq_other[0]*m_other)/(2.828427124746191*m_self+2.828427124746191*m_other)+(greene[0]*uSumSq[0]*m_other)/(2.828427124746191*m_self+2.828427124746191*m_other)-1.0*uM1Self[0]-1.0*uM1Other[0]+2.0*m2r[0]; 
  enRHS[1] = (-(2.0*greene[0]*vtsq_self[1]*m_self)/(2.828427124746191*m_self+2.828427124746191*m_other))-(1.0*greene[0]*uSumSq[1]*m_self)/(2.828427124746191*m_self+2.828427124746191*m_other)-(2.0*vtsq_self[0]*greene[1]*m_self)/(2.828427124746191*m_self+2.828427124746191*m_other)-(1.0*uSumSq[0]*greene[1]*m_self)/(2.828427124746191*m_self+2.828427124746191*m_other)+(2.0*greene[0]*vtsq_other[1]*m_other)/(2.828427124746191*m_self+2.828427124746191*m_other)+(greene[0]*uSumSq[1]*m_other)/(2.828427124746191*m_self+2.828427124746191*m_other)+(2.0*vtsq_other[0]*greene[1]*m_other)/(2.828427124746191*m_self+2.828427124746191*m_other)+(uSumSq[0]*greene[1]*m_other)/(2.828427124746191*m_self+2.828427124746191*m_other)-1.0*uM1Self[1]-1.0*uM1Other[1]+2.0*m2r[1]; 
 
  gkyl_mat_set(rhs,0,0,momRHS[0]); 
  gkyl_mat_set(rhs,1,0,momRHS[1]); 
  gkyl_mat_set(rhs,2,0,enRHS[0]); 
  gkyl_mat_set(rhs,3,0,enRHS[1]); 
} 
 
