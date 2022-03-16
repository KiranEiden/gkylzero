#include <gkyl_prim_lbo_gyrokinetic_kernels.h> 
 
GKYL_CU_DH void gyrokinetic_self_prim_moments_2x2v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *moms, const double *boundary_corrections) 
{ 
  // A:                    Matrix to be inverted to solve Ax = rhs (set by this function). 
  // rhs:                  right-hand side of Ax = rhs (set by this function). 
  // moms:                 moments of the distribution function (Zeroth, First, and Second in single array). 
  // boundary_corrections: boundary corrections to u and vtSq. 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if (-0.5*(3.872983346207417*(moms[7]+moms[6])-2.23606797749979*(moms[5]+moms[4])-3.0*moms[3]+1.732050807568877*(moms[2]+moms[1])-1.0*moms[0]) < 0) cellAvg = true; 
  if (0.5*(3.872983346207417*moms[7]-3.872983346207417*moms[6]+2.23606797749979*(moms[5]+moms[4])-3.0*moms[3]-1.732050807568877*moms[2]+1.732050807568877*moms[1]+moms[0]) < 0) cellAvg = true; 
  if (-0.5*(3.872983346207417*moms[7]-3.872983346207417*moms[6]-2.23606797749979*(moms[5]+moms[4])+3.0*moms[3]-1.732050807568877*moms[2]+1.732050807568877*moms[1]-1.0*moms[0]) < 0) cellAvg = true; 
  if (0.5*(3.872983346207417*(moms[7]+moms[6])+2.23606797749979*(moms[5]+moms[4])+3.0*moms[3]+1.732050807568877*(moms[2]+moms[1])+moms[0]) < 0) cellAvg = true; 
 
  double m0r[8] = {0.0}; 
  double m1r[8] = {0.0}; 
  double cMr[8] = {0.0}; 
  double cEr[8] = {0.0}; 
  if (cellAvg) { 
    m0r[0] = moms[0]; 
    m0r[1] = 0.0; 
    m0r[2] = 0.0; 
    m0r[3] = 0.0; 
    m0r[4] = 0.0; 
    m0r[5] = 0.0; 
    m0r[6] = 0.0; 
    m0r[7] = 0.0; 
    m1r[0] = moms[8]; 
    m1r[1] = 0.0; 
    m1r[2] = 0.0; 
    m1r[3] = 0.0; 
    m1r[4] = 0.0; 
    m1r[5] = 0.0; 
    m1r[6] = 0.0; 
    m1r[7] = 0.0; 
    gkyl_mat_set(rhs,0,0,moms[8]); 
    gkyl_mat_set(rhs,1,0,0.0); 
    gkyl_mat_set(rhs,2,0,0.0); 
    gkyl_mat_set(rhs,3,0,0.0); 
    gkyl_mat_set(rhs,4,0,0.0); 
    gkyl_mat_set(rhs,5,0,0.0); 
    gkyl_mat_set(rhs,6,0,0.0); 
    gkyl_mat_set(rhs,7,0,0.0); 
    cMr[0] = boundary_corrections[0]; 
    cMr[1] = 0.0; 
    cMr[2] = 0.0; 
    cMr[3] = 0.0; 
    cMr[4] = 0.0; 
    cMr[5] = 0.0; 
    cMr[6] = 0.0; 
    cMr[7] = 0.0; 
    cEr[0] = boundary_corrections[8]; 
    cEr[1] = 0.0; 
    cEr[2] = 0.0; 
    cEr[3] = 0.0; 
    cEr[4] = 0.0; 
    cEr[5] = 0.0; 
    cEr[6] = 0.0; 
    cEr[7] = 0.0; 
    gkyl_mat_set(rhs,8,0,moms[16]); 
    gkyl_mat_set(rhs,9,0,0.0); 
    gkyl_mat_set(rhs,10,0,0.0); 
    gkyl_mat_set(rhs,11,0,0.0); 
    gkyl_mat_set(rhs,12,0,0.0); 
    gkyl_mat_set(rhs,13,0,0.0); 
    gkyl_mat_set(rhs,14,0,0.0); 
    gkyl_mat_set(rhs,15,0,0.0); 
  } else { 
    m0r[0] = moms[0]; 
    m0r[1] = moms[1]; 
    m0r[2] = moms[2]; 
    m0r[3] = moms[3]; 
    m0r[4] = moms[4]; 
    m0r[5] = moms[5]; 
    m0r[6] = moms[6]; 
    m0r[7] = moms[7]; 
    m1r[0] = moms[8]; 
    m1r[1] = moms[9]; 
    m1r[2] = moms[10]; 
    m1r[3] = moms[11]; 
    m1r[4] = moms[12]; 
    m1r[5] = moms[13]; 
    m1r[6] = moms[14]; 
    m1r[7] = moms[15]; 
    gkyl_mat_set(rhs,0,0,moms[8]); 
    gkyl_mat_set(rhs,1,0,moms[9]); 
    gkyl_mat_set(rhs,2,0,moms[10]); 
    gkyl_mat_set(rhs,3,0,moms[11]); 
    gkyl_mat_set(rhs,4,0,moms[12]); 
    gkyl_mat_set(rhs,5,0,moms[13]); 
    gkyl_mat_set(rhs,6,0,moms[14]); 
    gkyl_mat_set(rhs,7,0,moms[15]); 
    cMr[0] = boundary_corrections[0]; 
    cMr[1] = boundary_corrections[1]; 
    cMr[2] = boundary_corrections[2]; 
    cMr[3] = boundary_corrections[3]; 
    cMr[4] = boundary_corrections[4]; 
    cMr[5] = boundary_corrections[5]; 
    cMr[6] = boundary_corrections[6]; 
    cMr[7] = boundary_corrections[7]; 
    cEr[0] = boundary_corrections[8]; 
    cEr[1] = boundary_corrections[9]; 
    cEr[2] = boundary_corrections[10]; 
    cEr[3] = boundary_corrections[11]; 
    cEr[4] = boundary_corrections[12]; 
    cEr[5] = boundary_corrections[13]; 
    cEr[6] = boundary_corrections[14]; 
    cEr[7] = boundary_corrections[15]; 
    gkyl_mat_set(rhs,8,0,moms[16]); 
    gkyl_mat_set(rhs,9,0,moms[17]); 
    gkyl_mat_set(rhs,10,0,moms[18]); 
    gkyl_mat_set(rhs,11,0,moms[19]); 
    gkyl_mat_set(rhs,12,0,moms[20]); 
    gkyl_mat_set(rhs,13,0,moms[21]); 
    gkyl_mat_set(rhs,14,0,moms[22]); 
    gkyl_mat_set(rhs,15,0,moms[23]); 
  } 
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  gkyl_mat_set(A,0,0,0.5*m0r[0]); 
  gkyl_mat_set(A,0,1,0.5*m0r[1]); 
  gkyl_mat_set(A,0,2,0.5*m0r[2]); 
  gkyl_mat_set(A,0,3,0.5*m0r[3]); 
  gkyl_mat_set(A,0,4,0.5*m0r[4]); 
  gkyl_mat_set(A,0,5,0.5*m0r[5]); 
  gkyl_mat_set(A,0,6,0.5*m0r[6]); 
  gkyl_mat_set(A,0,7,0.5*m0r[7]); 
  gkyl_mat_set(A,1,0,0.5*m0r[1]); 
  gkyl_mat_set(A,1,1,0.4472135954999579*m0r[4]+0.5*m0r[0]); 
  gkyl_mat_set(A,1,2,0.5*m0r[3]); 
  gkyl_mat_set(A,1,3,0.447213595499958*m0r[6]+0.5*m0r[2]); 
  gkyl_mat_set(A,1,4,0.4472135954999579*m0r[1]); 
  gkyl_mat_set(A,1,5,0.5000000000000001*m0r[7]); 
  gkyl_mat_set(A,1,6,0.447213595499958*m0r[3]); 
  gkyl_mat_set(A,1,7,0.5000000000000001*m0r[5]); 
  gkyl_mat_set(A,2,0,0.5*m0r[2]); 
  gkyl_mat_set(A,2,1,0.5*m0r[3]); 
  gkyl_mat_set(A,2,2,0.4472135954999579*m0r[5]+0.5*m0r[0]); 
  gkyl_mat_set(A,2,3,0.447213595499958*m0r[7]+0.5*m0r[1]); 
  gkyl_mat_set(A,2,4,0.5000000000000001*m0r[6]); 
  gkyl_mat_set(A,2,5,0.4472135954999579*m0r[2]); 
  gkyl_mat_set(A,2,6,0.5000000000000001*m0r[4]); 
  gkyl_mat_set(A,2,7,0.447213595499958*m0r[3]); 
  gkyl_mat_set(A,3,0,0.5*m0r[3]); 
  gkyl_mat_set(A,3,1,0.447213595499958*m0r[6]+0.5*m0r[2]); 
  gkyl_mat_set(A,3,2,0.447213595499958*m0r[7]+0.5*m0r[1]); 
  gkyl_mat_set(A,3,3,0.4472135954999579*m0r[5]+0.4472135954999579*m0r[4]+0.5*m0r[0]); 
  gkyl_mat_set(A,3,4,0.4472135954999579*m0r[3]); 
  gkyl_mat_set(A,3,5,0.4472135954999579*m0r[3]); 
  gkyl_mat_set(A,3,6,0.4*m0r[7]+0.447213595499958*m0r[1]); 
  gkyl_mat_set(A,3,7,0.4*m0r[6]+0.447213595499958*m0r[2]); 
  gkyl_mat_set(A,4,0,0.5*m0r[4]); 
  gkyl_mat_set(A,4,1,0.4472135954999579*m0r[1]); 
  gkyl_mat_set(A,4,2,0.5000000000000001*m0r[6]); 
  gkyl_mat_set(A,4,3,0.4472135954999579*m0r[3]); 
  gkyl_mat_set(A,4,4,0.31943828249997*m0r[4]+0.5*m0r[0]); 
  gkyl_mat_set(A,4,6,0.31943828249997*m0r[6]+0.5000000000000001*m0r[2]); 
  gkyl_mat_set(A,4,7,0.4472135954999579*m0r[7]); 
  gkyl_mat_set(A,5,0,0.5*m0r[5]); 
  gkyl_mat_set(A,5,1,0.5000000000000001*m0r[7]); 
  gkyl_mat_set(A,5,2,0.4472135954999579*m0r[2]); 
  gkyl_mat_set(A,5,3,0.4472135954999579*m0r[3]); 
  gkyl_mat_set(A,5,5,0.31943828249997*m0r[5]+0.5*m0r[0]); 
  gkyl_mat_set(A,5,6,0.4472135954999579*m0r[6]); 
  gkyl_mat_set(A,5,7,0.31943828249997*m0r[7]+0.5000000000000001*m0r[1]); 
  gkyl_mat_set(A,6,0,0.5*m0r[6]); 
  gkyl_mat_set(A,6,1,0.447213595499958*m0r[3]); 
  gkyl_mat_set(A,6,2,0.5000000000000001*m0r[4]); 
  gkyl_mat_set(A,6,3,0.4*m0r[7]+0.447213595499958*m0r[1]); 
  gkyl_mat_set(A,6,4,0.31943828249997*m0r[6]+0.5000000000000001*m0r[2]); 
  gkyl_mat_set(A,6,5,0.4472135954999579*m0r[6]); 
  gkyl_mat_set(A,6,6,0.4472135954999579*m0r[5]+0.31943828249997*m0r[4]+0.5*m0r[0]); 
  gkyl_mat_set(A,6,7,0.4*m0r[3]); 
  gkyl_mat_set(A,7,0,0.5*m0r[7]); 
  gkyl_mat_set(A,7,1,0.5000000000000001*m0r[5]); 
  gkyl_mat_set(A,7,2,0.447213595499958*m0r[3]); 
  gkyl_mat_set(A,7,3,0.4*m0r[6]+0.447213595499958*m0r[2]); 
  gkyl_mat_set(A,7,4,0.4472135954999579*m0r[7]); 
  gkyl_mat_set(A,7,5,0.31943828249997*m0r[7]+0.5000000000000001*m0r[1]); 
  gkyl_mat_set(A,7,6,0.4*m0r[3]); 
  gkyl_mat_set(A,7,7,0.31943828249997*m0r[5]+0.4472135954999579*m0r[4]+0.5*m0r[0]); 
 
  // ....... Block from correction to uX .......... // 
  gkyl_mat_set(A,0,8,-0.5*cMr[0]); 
  gkyl_mat_set(A,0,9,-0.5*cMr[1]); 
  gkyl_mat_set(A,0,10,-0.5*cMr[2]); 
  gkyl_mat_set(A,0,11,-0.5*cMr[3]); 
  gkyl_mat_set(A,0,12,-0.5*cMr[4]); 
  gkyl_mat_set(A,0,13,-0.5*cMr[5]); 
  gkyl_mat_set(A,0,14,-0.5*cMr[6]); 
  gkyl_mat_set(A,0,15,-0.5*cMr[7]); 
  gkyl_mat_set(A,1,8,-0.5*cMr[1]); 
  gkyl_mat_set(A,1,9,(-0.4472135954999579*cMr[4])-0.5*cMr[0]); 
  gkyl_mat_set(A,1,10,-0.5*cMr[3]); 
  gkyl_mat_set(A,1,11,(-0.447213595499958*cMr[6])-0.5*cMr[2]); 
  gkyl_mat_set(A,1,12,-0.4472135954999579*cMr[1]); 
  gkyl_mat_set(A,1,13,-0.5000000000000001*cMr[7]); 
  gkyl_mat_set(A,1,14,-0.447213595499958*cMr[3]); 
  gkyl_mat_set(A,1,15,-0.5000000000000001*cMr[5]); 
  gkyl_mat_set(A,2,8,-0.5*cMr[2]); 
  gkyl_mat_set(A,2,9,-0.5*cMr[3]); 
  gkyl_mat_set(A,2,10,(-0.4472135954999579*cMr[5])-0.5*cMr[0]); 
  gkyl_mat_set(A,2,11,(-0.447213595499958*cMr[7])-0.5*cMr[1]); 
  gkyl_mat_set(A,2,12,-0.5000000000000001*cMr[6]); 
  gkyl_mat_set(A,2,13,-0.4472135954999579*cMr[2]); 
  gkyl_mat_set(A,2,14,-0.5000000000000001*cMr[4]); 
  gkyl_mat_set(A,2,15,-0.447213595499958*cMr[3]); 
  gkyl_mat_set(A,3,8,-0.5*cMr[3]); 
  gkyl_mat_set(A,3,9,(-0.447213595499958*cMr[6])-0.5*cMr[2]); 
  gkyl_mat_set(A,3,10,(-0.447213595499958*cMr[7])-0.5*cMr[1]); 
  gkyl_mat_set(A,3,11,(-0.4472135954999579*cMr[5])-0.4472135954999579*cMr[4]-0.5*cMr[0]); 
  gkyl_mat_set(A,3,12,-0.4472135954999579*cMr[3]); 
  gkyl_mat_set(A,3,13,-0.4472135954999579*cMr[3]); 
  gkyl_mat_set(A,3,14,(-0.4*cMr[7])-0.447213595499958*cMr[1]); 
  gkyl_mat_set(A,3,15,(-0.4*cMr[6])-0.447213595499958*cMr[2]); 
  gkyl_mat_set(A,4,8,-0.5*cMr[4]); 
  gkyl_mat_set(A,4,9,-0.4472135954999579*cMr[1]); 
  gkyl_mat_set(A,4,10,-0.5000000000000001*cMr[6]); 
  gkyl_mat_set(A,4,11,-0.4472135954999579*cMr[3]); 
  gkyl_mat_set(A,4,12,(-0.31943828249997*cMr[4])-0.5*cMr[0]); 
  gkyl_mat_set(A,4,14,(-0.31943828249997*cMr[6])-0.5000000000000001*cMr[2]); 
  gkyl_mat_set(A,4,15,-0.4472135954999579*cMr[7]); 
  gkyl_mat_set(A,5,8,-0.5*cMr[5]); 
  gkyl_mat_set(A,5,9,-0.5000000000000001*cMr[7]); 
  gkyl_mat_set(A,5,10,-0.4472135954999579*cMr[2]); 
  gkyl_mat_set(A,5,11,-0.4472135954999579*cMr[3]); 
  gkyl_mat_set(A,5,13,(-0.31943828249997*cMr[5])-0.5*cMr[0]); 
  gkyl_mat_set(A,5,14,-0.4472135954999579*cMr[6]); 
  gkyl_mat_set(A,5,15,(-0.31943828249997*cMr[7])-0.5000000000000001*cMr[1]); 
  gkyl_mat_set(A,6,8,-0.5*cMr[6]); 
  gkyl_mat_set(A,6,9,-0.447213595499958*cMr[3]); 
  gkyl_mat_set(A,6,10,-0.5000000000000001*cMr[4]); 
  gkyl_mat_set(A,6,11,(-0.4*cMr[7])-0.447213595499958*cMr[1]); 
  gkyl_mat_set(A,6,12,(-0.31943828249997*cMr[6])-0.5000000000000001*cMr[2]); 
  gkyl_mat_set(A,6,13,-0.4472135954999579*cMr[6]); 
  gkyl_mat_set(A,6,14,(-0.4472135954999579*cMr[5])-0.31943828249997*cMr[4]-0.5*cMr[0]); 
  gkyl_mat_set(A,6,15,-0.4*cMr[3]); 
  gkyl_mat_set(A,7,8,-0.5*cMr[7]); 
  gkyl_mat_set(A,7,9,-0.5000000000000001*cMr[5]); 
  gkyl_mat_set(A,7,10,-0.447213595499958*cMr[3]); 
  gkyl_mat_set(A,7,11,(-0.4*cMr[6])-0.447213595499958*cMr[2]); 
  gkyl_mat_set(A,7,12,-0.4472135954999579*cMr[7]); 
  gkyl_mat_set(A,7,13,(-0.31943828249997*cMr[7])-0.5000000000000001*cMr[1]); 
  gkyl_mat_set(A,7,14,-0.4*cMr[3]); 
  gkyl_mat_set(A,7,15,(-0.31943828249997*cMr[5])-0.4472135954999579*cMr[4]-0.5*cMr[0]); 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  gkyl_mat_set(A,8,0,0.5*m1r[0]); 
  gkyl_mat_set(A,8,1,0.5*m1r[1]); 
  gkyl_mat_set(A,8,2,0.5*m1r[2]); 
  gkyl_mat_set(A,8,3,0.5*m1r[3]); 
  gkyl_mat_set(A,8,4,0.5*m1r[4]); 
  gkyl_mat_set(A,8,5,0.5*m1r[5]); 
  gkyl_mat_set(A,8,6,0.5*m1r[6]); 
  gkyl_mat_set(A,8,7,0.5*m1r[7]); 
  gkyl_mat_set(A,9,0,0.5*m1r[1]); 
  gkyl_mat_set(A,9,1,0.4472135954999579*m1r[4]+0.5*m1r[0]); 
  gkyl_mat_set(A,9,2,0.5*m1r[3]); 
  gkyl_mat_set(A,9,3,0.447213595499958*m1r[6]+0.5*m1r[2]); 
  gkyl_mat_set(A,9,4,0.4472135954999579*m1r[1]); 
  gkyl_mat_set(A,9,5,0.5000000000000001*m1r[7]); 
  gkyl_mat_set(A,9,6,0.447213595499958*m1r[3]); 
  gkyl_mat_set(A,9,7,0.5000000000000001*m1r[5]); 
  gkyl_mat_set(A,10,0,0.5*m1r[2]); 
  gkyl_mat_set(A,10,1,0.5*m1r[3]); 
  gkyl_mat_set(A,10,2,0.4472135954999579*m1r[5]+0.5*m1r[0]); 
  gkyl_mat_set(A,10,3,0.447213595499958*m1r[7]+0.5*m1r[1]); 
  gkyl_mat_set(A,10,4,0.5000000000000001*m1r[6]); 
  gkyl_mat_set(A,10,5,0.4472135954999579*m1r[2]); 
  gkyl_mat_set(A,10,6,0.5000000000000001*m1r[4]); 
  gkyl_mat_set(A,10,7,0.447213595499958*m1r[3]); 
  gkyl_mat_set(A,11,0,0.5*m1r[3]); 
  gkyl_mat_set(A,11,1,0.447213595499958*m1r[6]+0.5*m1r[2]); 
  gkyl_mat_set(A,11,2,0.447213595499958*m1r[7]+0.5*m1r[1]); 
  gkyl_mat_set(A,11,3,0.4472135954999579*m1r[5]+0.4472135954999579*m1r[4]+0.5*m1r[0]); 
  gkyl_mat_set(A,11,4,0.4472135954999579*m1r[3]); 
  gkyl_mat_set(A,11,5,0.4472135954999579*m1r[3]); 
  gkyl_mat_set(A,11,6,0.4*m1r[7]+0.447213595499958*m1r[1]); 
  gkyl_mat_set(A,11,7,0.4*m1r[6]+0.447213595499958*m1r[2]); 
  gkyl_mat_set(A,12,0,0.5*m1r[4]); 
  gkyl_mat_set(A,12,1,0.4472135954999579*m1r[1]); 
  gkyl_mat_set(A,12,2,0.5000000000000001*m1r[6]); 
  gkyl_mat_set(A,12,3,0.4472135954999579*m1r[3]); 
  gkyl_mat_set(A,12,4,0.31943828249997*m1r[4]+0.5*m1r[0]); 
  gkyl_mat_set(A,12,6,0.31943828249997*m1r[6]+0.5000000000000001*m1r[2]); 
  gkyl_mat_set(A,12,7,0.4472135954999579*m1r[7]); 
  gkyl_mat_set(A,13,0,0.5*m1r[5]); 
  gkyl_mat_set(A,13,1,0.5000000000000001*m1r[7]); 
  gkyl_mat_set(A,13,2,0.4472135954999579*m1r[2]); 
  gkyl_mat_set(A,13,3,0.4472135954999579*m1r[3]); 
  gkyl_mat_set(A,13,5,0.31943828249997*m1r[5]+0.5*m1r[0]); 
  gkyl_mat_set(A,13,6,0.4472135954999579*m1r[6]); 
  gkyl_mat_set(A,13,7,0.31943828249997*m1r[7]+0.5000000000000001*m1r[1]); 
  gkyl_mat_set(A,14,0,0.5*m1r[6]); 
  gkyl_mat_set(A,14,1,0.447213595499958*m1r[3]); 
  gkyl_mat_set(A,14,2,0.5000000000000001*m1r[4]); 
  gkyl_mat_set(A,14,3,0.4*m1r[7]+0.447213595499958*m1r[1]); 
  gkyl_mat_set(A,14,4,0.31943828249997*m1r[6]+0.5000000000000001*m1r[2]); 
  gkyl_mat_set(A,14,5,0.4472135954999579*m1r[6]); 
  gkyl_mat_set(A,14,6,0.4472135954999579*m1r[5]+0.31943828249997*m1r[4]+0.5*m1r[0]); 
  gkyl_mat_set(A,14,7,0.4*m1r[3]); 
  gkyl_mat_set(A,15,0,0.5*m1r[7]); 
  gkyl_mat_set(A,15,1,0.5000000000000001*m1r[5]); 
  gkyl_mat_set(A,15,2,0.447213595499958*m1r[3]); 
  gkyl_mat_set(A,15,3,0.4*m1r[6]+0.447213595499958*m1r[2]); 
  gkyl_mat_set(A,15,4,0.4472135954999579*m1r[7]); 
  gkyl_mat_set(A,15,5,0.31943828249997*m1r[7]+0.5000000000000001*m1r[1]); 
  gkyl_mat_set(A,15,6,0.4*m1r[3]); 
  gkyl_mat_set(A,15,7,0.31943828249997*m1r[5]+0.4472135954999579*m1r[4]+0.5*m1r[0]); 
 
  // ....... Block from correction to vtSq .......... // 
  gkyl_mat_set(A,8,8,1.5*m0r[0]-0.5*cEr[0]); 
  gkyl_mat_set(A,8,9,1.5*m0r[1]-0.5*cEr[1]); 
  gkyl_mat_set(A,8,10,1.5*m0r[2]-0.5*cEr[2]); 
  gkyl_mat_set(A,8,11,1.5*m0r[3]-0.5*cEr[3]); 
  gkyl_mat_set(A,8,12,1.5*m0r[4]-0.5*cEr[4]); 
  gkyl_mat_set(A,8,13,1.5*m0r[5]-0.5*cEr[5]); 
  gkyl_mat_set(A,8,14,1.5*m0r[6]-0.5*cEr[6]); 
  gkyl_mat_set(A,8,15,1.5*m0r[7]-0.5*cEr[7]); 
  gkyl_mat_set(A,9,8,1.5*m0r[1]-0.5*cEr[1]); 
  gkyl_mat_set(A,9,9,1.341640786499874*m0r[4]-0.4472135954999579*cEr[4]+1.5*m0r[0]-0.5*cEr[0]); 
  gkyl_mat_set(A,9,10,1.5*m0r[3]-0.5*cEr[3]); 
  gkyl_mat_set(A,9,11,1.341640786499874*m0r[6]-0.447213595499958*cEr[6]+1.5*m0r[2]-0.5*cEr[2]); 
  gkyl_mat_set(A,9,12,1.341640786499874*m0r[1]-0.4472135954999579*cEr[1]); 
  gkyl_mat_set(A,9,13,1.5*m0r[7]-0.5000000000000001*cEr[7]); 
  gkyl_mat_set(A,9,14,1.341640786499874*m0r[3]-0.447213595499958*cEr[3]); 
  gkyl_mat_set(A,9,15,1.5*m0r[5]-0.5000000000000001*cEr[5]); 
  gkyl_mat_set(A,10,8,1.5*m0r[2]-0.5*cEr[2]); 
  gkyl_mat_set(A,10,9,1.5*m0r[3]-0.5*cEr[3]); 
  gkyl_mat_set(A,10,10,1.341640786499874*m0r[5]-0.4472135954999579*cEr[5]+1.5*m0r[0]-0.5*cEr[0]); 
  gkyl_mat_set(A,10,11,1.341640786499874*m0r[7]-0.447213595499958*cEr[7]+1.5*m0r[1]-0.5*cEr[1]); 
  gkyl_mat_set(A,10,12,1.5*m0r[6]-0.5000000000000001*cEr[6]); 
  gkyl_mat_set(A,10,13,1.341640786499874*m0r[2]-0.4472135954999579*cEr[2]); 
  gkyl_mat_set(A,10,14,1.5*m0r[4]-0.5000000000000001*cEr[4]); 
  gkyl_mat_set(A,10,15,1.341640786499874*m0r[3]-0.447213595499958*cEr[3]); 
  gkyl_mat_set(A,11,8,1.5*m0r[3]-0.5*cEr[3]); 
  gkyl_mat_set(A,11,9,1.341640786499874*m0r[6]-0.447213595499958*cEr[6]+1.5*m0r[2]-0.5*cEr[2]); 
  gkyl_mat_set(A,11,10,1.341640786499874*m0r[7]-0.447213595499958*cEr[7]+1.5*m0r[1]-0.5*cEr[1]); 
  gkyl_mat_set(A,11,11,1.341640786499874*m0r[5]-0.4472135954999579*cEr[5]+1.341640786499874*m0r[4]-0.4472135954999579*cEr[4]+1.5*m0r[0]-0.5*cEr[0]); 
  gkyl_mat_set(A,11,12,1.341640786499874*m0r[3]-0.4472135954999579*cEr[3]); 
  gkyl_mat_set(A,11,13,1.341640786499874*m0r[3]-0.4472135954999579*cEr[3]); 
  gkyl_mat_set(A,11,14,1.2*m0r[7]-0.4*cEr[7]+1.341640786499874*m0r[1]-0.447213595499958*cEr[1]); 
  gkyl_mat_set(A,11,15,1.2*m0r[6]-0.4*cEr[6]+1.341640786499874*m0r[2]-0.447213595499958*cEr[2]); 
  gkyl_mat_set(A,12,8,1.5*m0r[4]-0.5*cEr[4]); 
  gkyl_mat_set(A,12,9,1.341640786499874*m0r[1]-0.4472135954999579*cEr[1]); 
  gkyl_mat_set(A,12,10,1.5*m0r[6]-0.5000000000000001*cEr[6]); 
  gkyl_mat_set(A,12,11,1.341640786499874*m0r[3]-0.4472135954999579*cEr[3]); 
  gkyl_mat_set(A,12,12,0.9583148474999099*m0r[4]-0.31943828249997*cEr[4]+1.5*m0r[0]-0.5*cEr[0]); 
  gkyl_mat_set(A,12,14,0.9583148474999099*m0r[6]-0.31943828249997*cEr[6]+1.5*m0r[2]-0.5000000000000001*cEr[2]); 
  gkyl_mat_set(A,12,15,1.341640786499874*m0r[7]-0.4472135954999579*cEr[7]); 
  gkyl_mat_set(A,13,8,1.5*m0r[5]-0.5*cEr[5]); 
  gkyl_mat_set(A,13,9,1.5*m0r[7]-0.5000000000000001*cEr[7]); 
  gkyl_mat_set(A,13,10,1.341640786499874*m0r[2]-0.4472135954999579*cEr[2]); 
  gkyl_mat_set(A,13,11,1.341640786499874*m0r[3]-0.4472135954999579*cEr[3]); 
  gkyl_mat_set(A,13,13,0.9583148474999099*m0r[5]-0.31943828249997*cEr[5]+1.5*m0r[0]-0.5*cEr[0]); 
  gkyl_mat_set(A,13,14,1.341640786499874*m0r[6]-0.4472135954999579*cEr[6]); 
  gkyl_mat_set(A,13,15,0.9583148474999099*m0r[7]-0.31943828249997*cEr[7]+1.5*m0r[1]-0.5000000000000001*cEr[1]); 
  gkyl_mat_set(A,14,8,1.5*m0r[6]-0.5*cEr[6]); 
  gkyl_mat_set(A,14,9,1.341640786499874*m0r[3]-0.447213595499958*cEr[3]); 
  gkyl_mat_set(A,14,10,1.5*m0r[4]-0.5000000000000001*cEr[4]); 
  gkyl_mat_set(A,14,11,1.2*m0r[7]-0.4*cEr[7]+1.341640786499874*m0r[1]-0.447213595499958*cEr[1]); 
  gkyl_mat_set(A,14,12,0.9583148474999099*m0r[6]-0.31943828249997*cEr[6]+1.5*m0r[2]-0.5000000000000001*cEr[2]); 
  gkyl_mat_set(A,14,13,1.341640786499874*m0r[6]-0.4472135954999579*cEr[6]); 
  gkyl_mat_set(A,14,14,1.341640786499874*m0r[5]-0.4472135954999579*cEr[5]+0.9583148474999099*m0r[4]-0.31943828249997*cEr[4]+1.5*m0r[0]-0.5*cEr[0]); 
  gkyl_mat_set(A,14,15,1.2*m0r[3]-0.4*cEr[3]); 
  gkyl_mat_set(A,15,8,1.5*m0r[7]-0.5*cEr[7]); 
  gkyl_mat_set(A,15,9,1.5*m0r[5]-0.5000000000000001*cEr[5]); 
  gkyl_mat_set(A,15,10,1.341640786499874*m0r[3]-0.447213595499958*cEr[3]); 
  gkyl_mat_set(A,15,11,1.2*m0r[6]-0.4*cEr[6]+1.341640786499874*m0r[2]-0.447213595499958*cEr[2]); 
  gkyl_mat_set(A,15,12,1.341640786499874*m0r[7]-0.4472135954999579*cEr[7]); 
  gkyl_mat_set(A,15,13,0.9583148474999099*m0r[7]-0.31943828249997*cEr[7]+1.5*m0r[1]-0.5000000000000001*cEr[1]); 
  gkyl_mat_set(A,15,14,1.2*m0r[3]-0.4*cEr[3]); 
  gkyl_mat_set(A,15,15,0.9583148474999099*m0r[5]-0.31943828249997*cEr[5]+1.341640786499874*m0r[4]-0.4472135954999579*cEr[4]+1.5*m0r[0]-0.5*cEr[0]); 
 
} 
 