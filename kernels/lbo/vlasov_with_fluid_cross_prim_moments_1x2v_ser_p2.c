#include <gkyl_prim_lbo_vlasov_kernels.h> 
 
GKYL_CU_DH void vlasov_with_fluid_cross_prim_moments_1x2v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *greene, const double m_self, const double *u_self, const double *vtsq_self, const double m_other, const double *u_other, const double *vtsq_other, const double *moms, const double *fluid, const double *boundary_corrections) 
{ 
  // greene:               Greene's factor. 
  // m:                    mass. 
  // moms:                 moments of the distribution function. 
  // u,vtSq:               self primitive moments: mean flow velocity and thermal speed squared. 
  // boundary_corrections: corrections to momentum and energy conservation due to finite velocity space. 
  // uCross,vtSqCross:     cross primitive moments: mean flow velocity and thermal speed squared. 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if (0.7071067811865475*(2.23606797749979*moms[2]-1.732050807568877*moms[1]+moms[0]) < 0) cellAvg = true; 
  if (0.7071067811865475*(2.23606797749979*moms[2]+1.732050807568877*moms[1]+moms[0]) < 0) cellAvg = true; 
 
  double m0r[3] = {0.0}; 
  double m1r[6] = {0.0}; 
  double m2r[3] = {0.0}; 
  double cMr[6] = {0.0}; 
  double cEr[3] = {0.0}; 
  if (cellAvg) { 
    m0r[0] = moms[0]; 
    m0r[1] = 0.0; 
    m0r[2] = 0.0; 
    m1r[0] = moms[3]; 
    m1r[1] = 0.0; 
    m1r[2] = 0.0; 
    cMr[0] = boundary_corrections[0]; 
    cMr[1] = 0.0; 
    cMr[2] = 0.0; 
    m1r[3] = moms[6]; 
    m1r[4] = 0.0; 
    m1r[5] = 0.0; 
    cMr[3] = boundary_corrections[3]; 
    cMr[4] = 0.0; 
    cMr[5] = 0.0; 
    m2r[0] = moms[9]+fluid[0]; 
    m2r[1] = 0.0; 
    m2r[2] = 0.0; 
    cEr[0] = boundary_corrections[6]; 
    cEr[1] = 0.0; 
    cEr[2] = 0.0; 
  } else { 
    m0r[0] = moms[0]; 
    m0r[1] = moms[1]; 
    m0r[2] = moms[2]; 
    m1r[0] = moms[3]; 
    m1r[1] = moms[4]; 
    m1r[2] = moms[5]; 
    m1r[3] = moms[6]; 
    m1r[4] = moms[7]; 
    m1r[5] = moms[8]; 
    m2r[0] = moms[9]+fluid[0]; 
    m2r[1] = moms[10]+fluid[1]; 
    m2r[2] = moms[11]+fluid[2]; 
    cMr[0] = boundary_corrections[0]; 
    cMr[1] = boundary_corrections[1]; 
    cMr[2] = boundary_corrections[2]; 
    cMr[3] = boundary_corrections[3]; 
    cMr[4] = boundary_corrections[4]; 
    cMr[5] = boundary_corrections[5]; 
    cEr[0] = boundary_corrections[6]; 
    cEr[1] = boundary_corrections[7]; 
    cEr[2] = boundary_corrections[8]; 
  } 
 
  double momRHS[6] = {0.0}; 
 
  // ... Block from weak multiply of m_self, nu, M0 and uCrossX ... // 
  gkyl_mat_set(A,0,0,1.414213562373095*m0r[0]); 
  gkyl_mat_set(A,0,1,1.414213562373095*m0r[1]); 
  gkyl_mat_set(A,0,2,1.414213562373095*m0r[2]); 
  gkyl_mat_set(A,1,0,1.414213562373095*m0r[1]); 
  gkyl_mat_set(A,1,1,1.264911064067352*m0r[2]+1.414213562373095*m0r[0]); 
  gkyl_mat_set(A,1,2,1.264911064067352*m0r[1]); 
  gkyl_mat_set(A,2,0,1.414213562373095*m0r[2]); 
  gkyl_mat_set(A,2,1,1.264911064067352*m0r[1]); 
  gkyl_mat_set(A,2,2,0.9035079029052515*m0r[2]+1.414213562373095*m0r[0]); 
 
  // ... Block from correction to momentum conservation (self) ... // 
  gkyl_mat_set(A,0,6,-1.414213562373095*cMr[0]); 
  gkyl_mat_set(A,0,7,-1.414213562373095*cMr[1]); 
  gkyl_mat_set(A,0,8,-1.414213562373095*cMr[2]); 
  gkyl_mat_set(A,1,6,-1.414213562373095*cMr[1]); 
  gkyl_mat_set(A,1,7,(-1.264911064067352*cMr[2])-1.414213562373095*cMr[0]); 
  gkyl_mat_set(A,1,8,-1.264911064067352*cMr[1]); 
  gkyl_mat_set(A,2,6,-1.414213562373095*cMr[2]); 
  gkyl_mat_set(A,2,7,-1.264911064067352*cMr[1]); 
  gkyl_mat_set(A,2,8,(-0.9035079029052515*cMr[2])-1.414213562373095*cMr[0]); 
 
  // ... Block from weak multiply of m_self, nu, m1X and uCrossX ... // 
  gkyl_mat_set(A,6,0,(-0.5*m0r[2]*u_self[2])-0.5*m0r[2]*u_other[2]-0.5*m0r[1]*u_self[1]-0.5*m0r[1]*u_other[1]-0.5*m0r[0]*u_self[0]-0.5*m0r[0]*u_other[0]+1.414213562373095*m1r[0]); 
  gkyl_mat_set(A,6,1,(-0.4472135954999579*m0r[1]*u_self[2])-0.4472135954999579*m0r[1]*u_other[2]-0.4472135954999579*u_self[1]*m0r[2]-0.4472135954999579*u_other[1]*m0r[2]-0.5*m0r[0]*u_self[1]-0.5*m0r[0]*u_other[1]+1.414213562373095*m1r[1]-0.5*u_self[0]*m0r[1]-0.5*u_other[0]*m0r[1]); 
  gkyl_mat_set(A,6,2,(-0.31943828249997*m0r[2]*u_self[2])-0.5*m0r[0]*u_self[2]-0.31943828249997*m0r[2]*u_other[2]-0.5*m0r[0]*u_other[2]+1.414213562373095*m1r[2]-0.5*u_self[0]*m0r[2]-0.5*u_other[0]*m0r[2]-0.4472135954999579*m0r[1]*u_self[1]-0.4472135954999579*m0r[1]*u_other[1]); 
  gkyl_mat_set(A,7,0,(-0.4472135954999579*m0r[1]*u_self[2])-0.4472135954999579*m0r[1]*u_other[2]-0.4472135954999579*u_self[1]*m0r[2]-0.4472135954999579*u_other[1]*m0r[2]-0.5*m0r[0]*u_self[1]-0.5*m0r[0]*u_other[1]+1.414213562373095*m1r[1]-0.5*u_self[0]*m0r[1]-0.5*u_other[0]*m0r[1]); 
  gkyl_mat_set(A,7,1,(-0.7857142857142857*m0r[2]*u_self[2])-0.4472135954999579*m0r[0]*u_self[2]-0.7857142857142857*m0r[2]*u_other[2]-0.4472135954999579*m0r[0]*u_other[2]+1.264911064067352*m1r[2]-0.4472135954999579*u_self[0]*m0r[2]-0.4472135954999579*u_other[0]*m0r[2]-0.9*m0r[1]*u_self[1]-0.9*m0r[1]*u_other[1]-0.5*m0r[0]*u_self[0]-0.5*m0r[0]*u_other[0]+1.414213562373095*m1r[0]); 
  gkyl_mat_set(A,7,2,(-0.7857142857142857*m0r[1]*u_self[2])-0.7857142857142857*m0r[1]*u_other[2]-0.7857142857142857*u_self[1]*m0r[2]-0.7857142857142857*u_other[1]*m0r[2]-0.4472135954999579*m0r[0]*u_self[1]-0.4472135954999579*m0r[0]*u_other[1]+1.264911064067352*m1r[1]-0.4472135954999579*u_self[0]*m0r[1]-0.4472135954999579*u_other[0]*m0r[1]); 
  gkyl_mat_set(A,8,0,(-0.31943828249997*m0r[2]*u_self[2])-0.5*m0r[0]*u_self[2]-0.31943828249997*m0r[2]*u_other[2]-0.5*m0r[0]*u_other[2]+1.414213562373095*m1r[2]-0.5*u_self[0]*m0r[2]-0.5*u_other[0]*m0r[2]-0.4472135954999579*m0r[1]*u_self[1]-0.4472135954999579*m0r[1]*u_other[1]); 
  gkyl_mat_set(A,8,1,(-0.7857142857142857*m0r[1]*u_self[2])-0.7857142857142857*m0r[1]*u_other[2]-0.7857142857142857*u_self[1]*m0r[2]-0.7857142857142857*u_other[1]*m0r[2]-0.4472135954999579*m0r[0]*u_self[1]-0.4472135954999579*m0r[0]*u_other[1]+1.264911064067352*m1r[1]-0.4472135954999579*u_self[0]*m0r[1]-0.4472135954999579*u_other[0]*m0r[1]); 
  gkyl_mat_set(A,8,2,(-1.071428571428571*m0r[2]*u_self[2])-0.31943828249997*m0r[0]*u_self[2]-1.071428571428571*m0r[2]*u_other[2]-0.31943828249997*m0r[0]*u_other[2]+0.9035079029052515*m1r[2]-0.31943828249997*u_self[0]*m0r[2]-0.31943828249997*u_other[0]*m0r[2]-0.7857142857142857*m0r[1]*u_self[1]-0.7857142857142857*m0r[1]*u_other[1]-0.5*m0r[0]*u_self[0]-0.5*m0r[0]*u_other[0]+1.414213562373095*m1r[0]); 
 
  momRHS[0] += -0.7071067811865475*(greene[2]*u_self[2]-1.0*greene[2]*u_other[2]+greene[1]*u_self[1]-1.0*greene[1]*u_other[1]+greene[0]*u_self[0]-1.0*greene[0]*u_other[0]-2.828427124746191*m1r[0]); 
  momRHS[1] += -0.1414213562373095*(4.47213595499958*greene[1]*u_self[2]-4.47213595499958*greene[1]*u_other[2]+(4.47213595499958*u_self[1]-4.47213595499958*u_other[1])*greene[2]+5.0*greene[0]*u_self[1]-5.0*greene[0]*u_other[1]-14.14213562373095*m1r[1]+(5.0*u_self[0]-5.0*u_other[0])*greene[1]); 
  momRHS[2] += -0.02020305089104421*((22.3606797749979*greene[2]+35.0*greene[0])*u_self[2]+((-22.3606797749979*greene[2])-35.0*greene[0])*u_other[2]-98.99494936611667*m1r[2]+(35.0*u_self[0]-35.0*u_other[0])*greene[2]+31.30495168499706*greene[1]*u_self[1]-31.30495168499706*greene[1]*u_other[1]); 
 
  // ... Block from weak multiply of m_self, nu, M0 and uCrossY ... // 
  gkyl_mat_set(A,3,3,1.414213562373095*m0r[0]); 
  gkyl_mat_set(A,3,4,1.414213562373095*m0r[1]); 
  gkyl_mat_set(A,3,5,1.414213562373095*m0r[2]); 
  gkyl_mat_set(A,4,3,1.414213562373095*m0r[1]); 
  gkyl_mat_set(A,4,4,1.264911064067352*m0r[2]+1.414213562373095*m0r[0]); 
  gkyl_mat_set(A,4,5,1.264911064067352*m0r[1]); 
  gkyl_mat_set(A,5,3,1.414213562373095*m0r[2]); 
  gkyl_mat_set(A,5,4,1.264911064067352*m0r[1]); 
  gkyl_mat_set(A,5,5,0.9035079029052515*m0r[2]+1.414213562373095*m0r[0]); 
 
  // ... Block from correction to momentum conservation (self) ... // 
  gkyl_mat_set(A,3,6,-1.414213562373095*cMr[3]); 
  gkyl_mat_set(A,3,7,-1.414213562373095*cMr[4]); 
  gkyl_mat_set(A,3,8,-1.414213562373095*cMr[5]); 
  gkyl_mat_set(A,4,6,-1.414213562373095*cMr[4]); 
  gkyl_mat_set(A,4,7,(-1.264911064067352*cMr[5])-1.414213562373095*cMr[3]); 
  gkyl_mat_set(A,4,8,-1.264911064067352*cMr[4]); 
  gkyl_mat_set(A,5,6,-1.414213562373095*cMr[5]); 
  gkyl_mat_set(A,5,7,-1.264911064067352*cMr[4]); 
  gkyl_mat_set(A,5,8,(-0.9035079029052515*cMr[5])-1.414213562373095*cMr[3]); 
 
  // ... Block from weak multiply of m_self, nu, m1Y and uCrossY ... // 
  gkyl_mat_set(A,6,3,(-0.5*m0r[2]*u_self[5])-0.5*m0r[2]*u_other[5]-0.5*m0r[1]*u_self[4]-0.5*m0r[1]*u_other[4]-0.5*m0r[0]*u_self[3]-0.5*m0r[0]*u_other[3]+1.414213562373095*m1r[3]); 
  gkyl_mat_set(A,6,4,(-0.4472135954999579*m0r[1]*u_self[5])-0.4472135954999579*m0r[1]*u_other[5]-0.4472135954999579*m0r[2]*u_self[4]-0.5*m0r[0]*u_self[4]-0.4472135954999579*m0r[2]*u_other[4]-0.5*m0r[0]*u_other[4]+1.414213562373095*m1r[4]-0.5*m0r[1]*u_self[3]-0.5*m0r[1]*u_other[3]); 
  gkyl_mat_set(A,6,5,(-0.31943828249997*m0r[2]*u_self[5])-0.5*m0r[0]*u_self[5]-0.31943828249997*m0r[2]*u_other[5]-0.5*m0r[0]*u_other[5]+1.414213562373095*m1r[5]-0.4472135954999579*m0r[1]*u_self[4]-0.4472135954999579*m0r[1]*u_other[4]-0.5*m0r[2]*u_self[3]-0.5*m0r[2]*u_other[3]); 
  gkyl_mat_set(A,7,3,(-0.4472135954999579*m0r[1]*u_self[5])-0.4472135954999579*m0r[1]*u_other[5]-0.4472135954999579*m0r[2]*u_self[4]-0.5*m0r[0]*u_self[4]-0.4472135954999579*m0r[2]*u_other[4]-0.5*m0r[0]*u_other[4]+1.414213562373095*m1r[4]-0.5*m0r[1]*u_self[3]-0.5*m0r[1]*u_other[3]); 
  gkyl_mat_set(A,7,4,(-0.7857142857142857*m0r[2]*u_self[5])-0.4472135954999579*m0r[0]*u_self[5]-0.7857142857142857*m0r[2]*u_other[5]-0.4472135954999579*m0r[0]*u_other[5]+1.264911064067352*m1r[5]-0.9*m0r[1]*u_self[4]-0.9*m0r[1]*u_other[4]-0.4472135954999579*m0r[2]*u_self[3]-0.5*m0r[0]*u_self[3]-0.4472135954999579*m0r[2]*u_other[3]-0.5*m0r[0]*u_other[3]+1.414213562373095*m1r[3]); 
  gkyl_mat_set(A,7,5,(-0.7857142857142857*m0r[1]*u_self[5])-0.7857142857142857*m0r[1]*u_other[5]-0.7857142857142857*m0r[2]*u_self[4]-0.4472135954999579*m0r[0]*u_self[4]-0.7857142857142857*m0r[2]*u_other[4]-0.4472135954999579*m0r[0]*u_other[4]+1.264911064067352*m1r[4]-0.4472135954999579*m0r[1]*u_self[3]-0.4472135954999579*m0r[1]*u_other[3]); 
  gkyl_mat_set(A,8,3,(-0.31943828249997*m0r[2]*u_self[5])-0.5*m0r[0]*u_self[5]-0.31943828249997*m0r[2]*u_other[5]-0.5*m0r[0]*u_other[5]+1.414213562373095*m1r[5]-0.4472135954999579*m0r[1]*u_self[4]-0.4472135954999579*m0r[1]*u_other[4]-0.5*m0r[2]*u_self[3]-0.5*m0r[2]*u_other[3]); 
  gkyl_mat_set(A,8,4,(-0.7857142857142857*m0r[1]*u_self[5])-0.7857142857142857*m0r[1]*u_other[5]-0.7857142857142857*m0r[2]*u_self[4]-0.4472135954999579*m0r[0]*u_self[4]-0.7857142857142857*m0r[2]*u_other[4]-0.4472135954999579*m0r[0]*u_other[4]+1.264911064067352*m1r[4]-0.4472135954999579*m0r[1]*u_self[3]-0.4472135954999579*m0r[1]*u_other[3]); 
  gkyl_mat_set(A,8,5,(-1.071428571428571*m0r[2]*u_self[5])-0.31943828249997*m0r[0]*u_self[5]-1.071428571428571*m0r[2]*u_other[5]-0.31943828249997*m0r[0]*u_other[5]+0.9035079029052515*m1r[5]-0.7857142857142857*m0r[1]*u_self[4]-0.7857142857142857*m0r[1]*u_other[4]-0.31943828249997*m0r[2]*u_self[3]-0.5*m0r[0]*u_self[3]-0.31943828249997*m0r[2]*u_other[3]-0.5*m0r[0]*u_other[3]+1.414213562373095*m1r[3]); 
 
  momRHS[3] += -0.7071067811865475*(greene[2]*u_self[5]-1.0*greene[2]*u_other[5]+greene[1]*u_self[4]-1.0*greene[1]*u_other[4]+greene[0]*u_self[3]-1.0*greene[0]*u_other[3]-2.828427124746191*m1r[3]); 
  momRHS[4] += -0.1414213562373095*(4.47213595499958*greene[1]*u_self[5]-4.47213595499958*greene[1]*u_other[5]+(4.47213595499958*greene[2]+5.0*greene[0])*u_self[4]+((-4.47213595499958*greene[2])-5.0*greene[0])*u_other[4]-14.14213562373095*m1r[4]+5.0*greene[1]*u_self[3]-5.0*greene[1]*u_other[3]); 
  momRHS[5] += -0.02020305089104421*((22.3606797749979*greene[2]+35.0*greene[0])*u_self[5]+((-22.3606797749979*greene[2])-35.0*greene[0])*u_other[5]-98.99494936611667*m1r[5]+31.30495168499706*greene[1]*u_self[4]-31.30495168499706*greene[1]*u_other[4]+35.0*greene[2]*u_self[3]-35.0*greene[2]*u_other[3]); 
 
  double ucMSelf[3] = {0.0}; 
  double ucMOther[3] = {0.0}; 
  for (int vd=0; vd<2; vd++) 
  { 
    int a0 = 3*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    ucMSelf[0] += 0.7071067811865475*cMr[a0+2]*u_self[a0+2]+0.7071067811865475*cMr[a0+1]*u_self[a0+1]+0.7071067811865475*cMr[a0]*u_self[a0]; 
    ucMSelf[1] += 0.6324555320336759*cMr[a0+1]*u_self[a0+2]+0.6324555320336759*u_self[a0+1]*cMr[a0+2]+0.7071067811865475*cMr[a0]*u_self[a0+1]+0.7071067811865475*u_self[a0]*cMr[a0+1]; 
    ucMSelf[2] += 0.4517539514526256*cMr[a0+2]*u_self[a0+2]+0.7071067811865475*cMr[a0]*u_self[a0+2]+0.7071067811865475*u_self[a0]*cMr[a0+2]+0.6324555320336759*cMr[a0+1]*u_self[a0+1]; 
    ucMOther[0] += 0.7071067811865475*cMr[a0+2]*u_other[a0+2]+0.7071067811865475*cMr[a0+1]*u_other[a0+1]+0.7071067811865475*cMr[a0]*u_other[a0]; 
    ucMOther[1] += 0.6324555320336759*cMr[a0+1]*u_other[a0+2]+0.6324555320336759*u_other[a0+1]*cMr[a0+2]+0.7071067811865475*cMr[a0]*u_other[a0+1]+0.7071067811865475*u_other[a0]*cMr[a0+1]; 
    ucMOther[2] += 0.4517539514526256*cMr[a0+2]*u_other[a0+2]+0.7071067811865475*cMr[a0]*u_other[a0+2]+0.7071067811865475*u_other[a0]*cMr[a0+2]+0.6324555320336759*cMr[a0+1]*u_other[a0+1]; 
  } 
 
  // ... Block from correction to (self) 2nd moment of collision operator ... // 
  gkyl_mat_set(A,6,6,0.7071067811865475*ucMSelf[0]+0.7071067811865475*ucMOther[0]+4.242640687119286*m0r[0]-1.414213562373095*cEr[0]); 
  gkyl_mat_set(A,6,7,0.7071067811865475*ucMSelf[1]+0.7071067811865475*ucMOther[1]+4.242640687119286*m0r[1]-1.414213562373095*cEr[1]); 
  gkyl_mat_set(A,6,8,0.7071067811865475*ucMSelf[2]+0.7071067811865475*ucMOther[2]+4.242640687119286*m0r[2]-1.414213562373095*cEr[2]); 
  gkyl_mat_set(A,7,6,0.7071067811865475*ucMSelf[1]+0.7071067811865475*ucMOther[1]+4.242640687119286*m0r[1]-1.414213562373095*cEr[1]); 
  gkyl_mat_set(A,7,7,0.6324555320336759*ucMSelf[2]+0.6324555320336759*ucMOther[2]+3.794733192202056*m0r[2]-1.264911064067352*cEr[2]+0.7071067811865475*ucMSelf[0]+0.7071067811865475*ucMOther[0]+4.242640687119286*m0r[0]-1.414213562373095*cEr[0]); 
  gkyl_mat_set(A,7,8,0.6324555320336759*ucMSelf[1]+0.6324555320336759*ucMOther[1]+3.794733192202056*m0r[1]-1.264911064067352*cEr[1]); 
  gkyl_mat_set(A,8,6,0.7071067811865475*ucMSelf[2]+0.7071067811865475*ucMOther[2]+4.242640687119286*m0r[2]-1.414213562373095*cEr[2]); 
  gkyl_mat_set(A,8,7,0.6324555320336759*ucMSelf[1]+0.6324555320336759*ucMOther[1]+3.794733192202056*m0r[1]-1.264911064067352*cEr[1]); 
  gkyl_mat_set(A,8,8,0.4517539514526256*ucMSelf[2]+0.4517539514526256*ucMOther[2]+2.710523708715754*m0r[2]-0.9035079029052515*cEr[2]+0.7071067811865475*ucMSelf[0]+0.7071067811865475*ucMOther[0]+4.242640687119286*m0r[0]-1.414213562373095*cEr[0]); 
 
  double uM1Self[3] = {0.0}; 
  double uM1Other[3] = {0.0}; 
  double uSumSq[3] = {0.0}; 
  for (int vd=0; vd<2; vd++) 
  { 
    int a0 = 3*vd; 
    // Dot product terms in energy equation RHS. 
    uM1Self[0] += 0.7071067811865475*m1r[5]*u_self[a0+2]+0.7071067811865475*m1r[4]*u_self[a0+1]+0.7071067811865475*m1r[3]*u_self[a0]; 
    uM1Self[1] += 0.6324555320336759*m1r[4]*u_self[a0+2]+0.6324555320336759*m1r[5]*u_self[a0+1]+0.7071067811865475*m1r[3]*u_self[a0+1]+0.7071067811865475*m1r[4]*u_self[a0]; 
    uM1Self[2] += 0.4517539514526256*m1r[5]*u_self[a0+2]+0.7071067811865475*m1r[3]*u_self[a0+2]+0.6324555320336759*m1r[4]*u_self[a0+1]+0.7071067811865475*m1r[5]*u_self[a0]; 
    uM1Other[0] += 0.7071067811865475*m1r[5]*u_other[a0+2]+0.7071067811865475*m1r[4]*u_other[a0+1]+0.7071067811865475*m1r[3]*u_other[a0]; 
    uM1Other[1] += 0.6324555320336759*m1r[4]*u_other[a0+2]+0.6324555320336759*m1r[5]*u_other[a0+1]+0.7071067811865475*m1r[3]*u_other[a0+1]+0.7071067811865475*m1r[4]*u_other[a0]; 
    uM1Other[2] += 0.4517539514526256*m1r[5]*u_other[a0+2]+0.7071067811865475*m1r[3]*u_other[a0+2]+0.6324555320336759*m1r[4]*u_other[a0+1]+0.7071067811865475*m1r[5]*u_other[a0]; 
  const double u_self0R2 = pow(u_self[a0],2);
  const double u_self1R2 = pow(u_self[a0+1],2);
  const double u_self2R2 = pow(u_self[a0+2],2);
  const double u_other0R2 = pow(u_other[a0],2);
  const double u_other1R2 = pow(u_other[a0+1],2);
  const double u_other2R2 = pow(u_other[a0+2],2);

  uSumSq[0] += 0.7071067811865475*u_self2R2-1.414213562373095*u_other[a0+2]*u_self[a0+2]+0.7071067811865475*u_other2R2+0.7071067811865475*u_self1R2-1.414213562373095*u_other[a0+1]*u_self[a0+1]+0.7071067811865475*u_other1R2+0.7071067811865475*u_self0R2-1.414213562373095*u_other[a0]*u_self[a0]+0.7071067811865475*u_other0R2; 
  uSumSq[1] += 1.264911064067352*u_self[a0+1]*u_self[a0+2]-1.264911064067352*u_other[a0+1]*u_self[a0+2]-1.264911064067352*u_self[a0+1]*u_other[a0+2]+1.264911064067352*u_other[a0+1]*u_other[a0+2]+1.414213562373095*u_self[a0]*u_self[a0+1]-1.414213562373095*u_other[a0]*u_self[a0+1]-1.414213562373095*u_self[a0]*u_other[a0+1]+1.414213562373095*u_other[a0]*u_other[a0+1]; 
  uSumSq[2] += 0.4517539514526256*u_self2R2-0.9035079029052515*u_other[a0+2]*u_self[a0+2]+1.414213562373095*u_self[a0]*u_self[a0+2]-1.414213562373095*u_other[a0]*u_self[a0+2]+0.4517539514526256*u_other2R2-1.414213562373095*u_self[a0]*u_other[a0+2]+1.414213562373095*u_other[a0]*u_other[a0+2]+0.6324555320336759*u_self1R2-1.264911064067352*u_other[a0+1]*u_self[a0+1]+0.6324555320336759*u_other1R2; 
  } 
 
  double enRHS[3] = {0.0}; 
 
  enRHS[0] = (-(6.0*greene[2]*vtsq_self[2]*m_self)/(2.828427124746191*m_self+2.828427124746191*m_other))-(1.0*greene[2]*uSumSq[2]*m_self)/(2.828427124746191*m_self+2.828427124746191*m_other)-(6.0*greene[1]*vtsq_self[1]*m_self)/(2.828427124746191*m_self+2.828427124746191*m_other)-(1.0*greene[1]*uSumSq[1]*m_self)/(2.828427124746191*m_self+2.828427124746191*m_other)-(6.0*greene[0]*vtsq_self[0]*m_self)/(2.828427124746191*m_self+2.828427124746191*m_other)-(1.0*greene[0]*uSumSq[0]*m_self)/(2.828427124746191*m_self+2.828427124746191*m_other)+(6.0*greene[2]*vtsq_other[2]*m_other)/(2.828427124746191*m_self+2.828427124746191*m_other)+(greene[2]*uSumSq[2]*m_other)/(2.828427124746191*m_self+2.828427124746191*m_other)+(6.0*greene[1]*vtsq_other[1]*m_other)/(2.828427124746191*m_self+2.828427124746191*m_other)+(greene[1]*uSumSq[1]*m_other)/(2.828427124746191*m_self+2.828427124746191*m_other)+(6.0*greene[0]*vtsq_other[0]*m_other)/(2.828427124746191*m_self+2.828427124746191*m_other)+(greene[0]*uSumSq[0]*m_other)/(2.828427124746191*m_self+2.828427124746191*m_other)-1.0*uM1Self[0]-1.0*uM1Other[0]+2.0*m2r[0]; 
  enRHS[1] = (-(26.83281572999748*greene[1]*vtsq_self[2]*m_self)/(14.14213562373095*m_self+14.14213562373095*m_other))-(4.47213595499958*greene[1]*uSumSq[2]*m_self)/(14.14213562373095*m_self+14.14213562373095*m_other)-(26.83281572999748*vtsq_self[1]*greene[2]*m_self)/(14.14213562373095*m_self+14.14213562373095*m_other)-(4.47213595499958*uSumSq[1]*greene[2]*m_self)/(14.14213562373095*m_self+14.14213562373095*m_other)-(30.0*greene[0]*vtsq_self[1]*m_self)/(14.14213562373095*m_self+14.14213562373095*m_other)-(5.0*greene[0]*uSumSq[1]*m_self)/(14.14213562373095*m_self+14.14213562373095*m_other)-(30.0*vtsq_self[0]*greene[1]*m_self)/(14.14213562373095*m_self+14.14213562373095*m_other)-(5.0*uSumSq[0]*greene[1]*m_self)/(14.14213562373095*m_self+14.14213562373095*m_other)+(26.83281572999748*greene[1]*vtsq_other[2]*m_other)/(14.14213562373095*m_self+14.14213562373095*m_other)+(4.47213595499958*greene[1]*uSumSq[2]*m_other)/(14.14213562373095*m_self+14.14213562373095*m_other)+(26.83281572999748*vtsq_other[1]*greene[2]*m_other)/(14.14213562373095*m_self+14.14213562373095*m_other)+(4.47213595499958*uSumSq[1]*greene[2]*m_other)/(14.14213562373095*m_self+14.14213562373095*m_other)+(30.0*greene[0]*vtsq_other[1]*m_other)/(14.14213562373095*m_self+14.14213562373095*m_other)+(5.0*greene[0]*uSumSq[1]*m_other)/(14.14213562373095*m_self+14.14213562373095*m_other)+(30.0*vtsq_other[0]*greene[1]*m_other)/(14.14213562373095*m_self+14.14213562373095*m_other)+(5.0*uSumSq[0]*greene[1]*m_other)/(14.14213562373095*m_self+14.14213562373095*m_other)-1.0*uM1Self[1]-1.0*uM1Other[1]+2.0*m2r[1]; 
  enRHS[2] = (-(134.1640786499874*greene[2]*vtsq_self[2]*m_self)/(98.99494936611667*m_self+98.99494936611667*m_other))-(210.0*greene[0]*vtsq_self[2]*m_self)/(98.99494936611667*m_self+98.99494936611667*m_other)-(22.3606797749979*greene[2]*uSumSq[2]*m_self)/(98.99494936611667*m_self+98.99494936611667*m_other)-(35.0*greene[0]*uSumSq[2]*m_self)/(98.99494936611667*m_self+98.99494936611667*m_other)-(210.0*vtsq_self[0]*greene[2]*m_self)/(98.99494936611667*m_self+98.99494936611667*m_other)-(35.0*uSumSq[0]*greene[2]*m_self)/(98.99494936611667*m_self+98.99494936611667*m_other)-(187.8297101099823*greene[1]*vtsq_self[1]*m_self)/(98.99494936611667*m_self+98.99494936611667*m_other)-(31.30495168499706*greene[1]*uSumSq[1]*m_self)/(98.99494936611667*m_self+98.99494936611667*m_other)+(134.1640786499874*greene[2]*vtsq_other[2]*m_other)/(98.99494936611667*m_self+98.99494936611667*m_other)+(210.0*greene[0]*vtsq_other[2]*m_other)/(98.99494936611667*m_self+98.99494936611667*m_other)+(22.3606797749979*greene[2]*uSumSq[2]*m_other)/(98.99494936611667*m_self+98.99494936611667*m_other)+(35.0*greene[0]*uSumSq[2]*m_other)/(98.99494936611667*m_self+98.99494936611667*m_other)+(210.0*vtsq_other[0]*greene[2]*m_other)/(98.99494936611667*m_self+98.99494936611667*m_other)+(35.0*uSumSq[0]*greene[2]*m_other)/(98.99494936611667*m_self+98.99494936611667*m_other)+(187.8297101099823*greene[1]*vtsq_other[1]*m_other)/(98.99494936611667*m_self+98.99494936611667*m_other)+(31.30495168499706*greene[1]*uSumSq[1]*m_other)/(98.99494936611667*m_self+98.99494936611667*m_other)-1.0*uM1Self[2]-1.0*uM1Other[2]+2.0*m2r[2]; 
 
  gkyl_mat_set(rhs,0,0,momRHS[0]); 
  gkyl_mat_set(rhs,1,0,momRHS[1]); 
  gkyl_mat_set(rhs,2,0,momRHS[2]); 
  gkyl_mat_set(rhs,3,0,momRHS[3]); 
  gkyl_mat_set(rhs,4,0,momRHS[4]); 
  gkyl_mat_set(rhs,5,0,momRHS[5]); 
  gkyl_mat_set(rhs,6,0,enRHS[0]); 
  gkyl_mat_set(rhs,7,0,enRHS[1]); 
  gkyl_mat_set(rhs,8,0,enRHS[2]); 
} 
 