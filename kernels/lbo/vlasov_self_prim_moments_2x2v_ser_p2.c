#include <gkyl_prim_lbo_vlasov_kernels.h> 
 
GKYL_CU_DH void vlasov_self_prim_moments_2x2v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *moms, const double *boundary_corrections) 
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
  double m1r[16] = {0.0}; 
  double cMr[16] = {0.0}; 
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
    m1r[8] = moms[16]; 
    m1r[9] = 0.0; 
    m1r[10] = 0.0; 
    m1r[11] = 0.0; 
    m1r[12] = 0.0; 
    m1r[13] = 0.0; 
    m1r[14] = 0.0; 
    m1r[15] = 0.0; 
    gkyl_mat_set(rhs,0,0,moms[16]); 
    gkyl_mat_set(rhs,1,0,0.0); 
    gkyl_mat_set(rhs,2,0,0.0); 
    gkyl_mat_set(rhs,3,0,0.0); 
    gkyl_mat_set(rhs,4,0,0.0); 
    gkyl_mat_set(rhs,5,0,0.0); 
    gkyl_mat_set(rhs,6,0,0.0); 
    gkyl_mat_set(rhs,7,0,0.0); 
    cMr[8] = boundary_corrections[8]; 
    cMr[9] = 0.0; 
    cMr[10] = 0.0; 
    cMr[11] = 0.0; 
    cMr[12] = 0.0; 
    cMr[13] = 0.0; 
    cMr[14] = 0.0; 
    cMr[15] = 0.0; 
    cEr[0] = boundary_corrections[16]; 
    cEr[1] = 0.0; 
    cEr[2] = 0.0; 
    cEr[3] = 0.0; 
    cEr[4] = 0.0; 
    cEr[5] = 0.0; 
    cEr[6] = 0.0; 
    cEr[7] = 0.0; 
    gkyl_mat_set(rhs,16,0,moms[24]); 
    gkyl_mat_set(rhs,17,0,0.0); 
    gkyl_mat_set(rhs,18,0,0.0); 
    gkyl_mat_set(rhs,19,0,0.0); 
    gkyl_mat_set(rhs,20,0,0.0); 
    gkyl_mat_set(rhs,21,0,0.0); 
    gkyl_mat_set(rhs,22,0,0.0); 
    gkyl_mat_set(rhs,23,0,0.0); 
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
    m1r[8] = moms[16]; 
    m1r[9] = moms[17]; 
    m1r[10] = moms[18]; 
    m1r[11] = moms[19]; 
    m1r[12] = moms[20]; 
    m1r[13] = moms[21]; 
    m1r[14] = moms[22]; 
    m1r[15] = moms[23]; 
    gkyl_mat_set(rhs,0,0,moms[8]); 
    gkyl_mat_set(rhs,1,0,moms[9]); 
    gkyl_mat_set(rhs,2,0,moms[10]); 
    gkyl_mat_set(rhs,3,0,moms[11]); 
    gkyl_mat_set(rhs,4,0,moms[12]); 
    gkyl_mat_set(rhs,5,0,moms[13]); 
    gkyl_mat_set(rhs,6,0,moms[14]); 
    gkyl_mat_set(rhs,7,0,moms[15]); 
    gkyl_mat_set(rhs,8,0,moms[16]); 
    gkyl_mat_set(rhs,9,0,moms[17]); 
    gkyl_mat_set(rhs,10,0,moms[18]); 
    gkyl_mat_set(rhs,11,0,moms[19]); 
    gkyl_mat_set(rhs,12,0,moms[20]); 
    gkyl_mat_set(rhs,13,0,moms[21]); 
    gkyl_mat_set(rhs,14,0,moms[22]); 
    gkyl_mat_set(rhs,15,0,moms[23]); 
    cMr[0] = boundary_corrections[0]; 
    cMr[1] = boundary_corrections[1]; 
    cMr[2] = boundary_corrections[2]; 
    cMr[3] = boundary_corrections[3]; 
    cMr[4] = boundary_corrections[4]; 
    cMr[5] = boundary_corrections[5]; 
    cMr[6] = boundary_corrections[6]; 
    cMr[7] = boundary_corrections[7]; 
    cMr[8] = boundary_corrections[8]; 
    cMr[9] = boundary_corrections[9]; 
    cMr[10] = boundary_corrections[10]; 
    cMr[11] = boundary_corrections[11]; 
    cMr[12] = boundary_corrections[12]; 
    cMr[13] = boundary_corrections[13]; 
    cMr[14] = boundary_corrections[14]; 
    cMr[15] = boundary_corrections[15]; 
    cEr[0] = boundary_corrections[16]; 
    cEr[1] = boundary_corrections[17]; 
    cEr[2] = boundary_corrections[18]; 
    cEr[3] = boundary_corrections[19]; 
    cEr[4] = boundary_corrections[20]; 
    cEr[5] = boundary_corrections[21]; 
    cEr[6] = boundary_corrections[22]; 
    cEr[7] = boundary_corrections[23]; 
    gkyl_mat_set(rhs,16,0,moms[24]); 
    gkyl_mat_set(rhs,17,0,moms[25]); 
    gkyl_mat_set(rhs,18,0,moms[26]); 
    gkyl_mat_set(rhs,19,0,moms[27]); 
    gkyl_mat_set(rhs,20,0,moms[28]); 
    gkyl_mat_set(rhs,21,0,moms[29]); 
    gkyl_mat_set(rhs,22,0,moms[30]); 
    gkyl_mat_set(rhs,23,0,moms[31]); 
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
  gkyl_mat_set(A,0,16,-0.5*cMr[0]); 
  gkyl_mat_set(A,0,17,-0.5*cMr[1]); 
  gkyl_mat_set(A,0,18,-0.5*cMr[2]); 
  gkyl_mat_set(A,0,19,-0.5*cMr[3]); 
  gkyl_mat_set(A,0,20,-0.5*cMr[4]); 
  gkyl_mat_set(A,0,21,-0.5*cMr[5]); 
  gkyl_mat_set(A,0,22,-0.5*cMr[6]); 
  gkyl_mat_set(A,0,23,-0.5*cMr[7]); 
  gkyl_mat_set(A,1,16,-0.5*cMr[1]); 
  gkyl_mat_set(A,1,17,(-0.4472135954999579*cMr[4])-0.5*cMr[0]); 
  gkyl_mat_set(A,1,18,-0.5*cMr[3]); 
  gkyl_mat_set(A,1,19,(-0.447213595499958*cMr[6])-0.5*cMr[2]); 
  gkyl_mat_set(A,1,20,-0.4472135954999579*cMr[1]); 
  gkyl_mat_set(A,1,21,-0.5000000000000001*cMr[7]); 
  gkyl_mat_set(A,1,22,-0.447213595499958*cMr[3]); 
  gkyl_mat_set(A,1,23,-0.5000000000000001*cMr[5]); 
  gkyl_mat_set(A,2,16,-0.5*cMr[2]); 
  gkyl_mat_set(A,2,17,-0.5*cMr[3]); 
  gkyl_mat_set(A,2,18,(-0.4472135954999579*cMr[5])-0.5*cMr[0]); 
  gkyl_mat_set(A,2,19,(-0.447213595499958*cMr[7])-0.5*cMr[1]); 
  gkyl_mat_set(A,2,20,-0.5000000000000001*cMr[6]); 
  gkyl_mat_set(A,2,21,-0.4472135954999579*cMr[2]); 
  gkyl_mat_set(A,2,22,-0.5000000000000001*cMr[4]); 
  gkyl_mat_set(A,2,23,-0.447213595499958*cMr[3]); 
  gkyl_mat_set(A,3,16,-0.5*cMr[3]); 
  gkyl_mat_set(A,3,17,(-0.447213595499958*cMr[6])-0.5*cMr[2]); 
  gkyl_mat_set(A,3,18,(-0.447213595499958*cMr[7])-0.5*cMr[1]); 
  gkyl_mat_set(A,3,19,(-0.4472135954999579*cMr[5])-0.4472135954999579*cMr[4]-0.5*cMr[0]); 
  gkyl_mat_set(A,3,20,-0.4472135954999579*cMr[3]); 
  gkyl_mat_set(A,3,21,-0.4472135954999579*cMr[3]); 
  gkyl_mat_set(A,3,22,(-0.4*cMr[7])-0.447213595499958*cMr[1]); 
  gkyl_mat_set(A,3,23,(-0.4*cMr[6])-0.447213595499958*cMr[2]); 
  gkyl_mat_set(A,4,16,-0.5*cMr[4]); 
  gkyl_mat_set(A,4,17,-0.4472135954999579*cMr[1]); 
  gkyl_mat_set(A,4,18,-0.5000000000000001*cMr[6]); 
  gkyl_mat_set(A,4,19,-0.4472135954999579*cMr[3]); 
  gkyl_mat_set(A,4,20,(-0.31943828249997*cMr[4])-0.5*cMr[0]); 
  gkyl_mat_set(A,4,22,(-0.31943828249997*cMr[6])-0.5000000000000001*cMr[2]); 
  gkyl_mat_set(A,4,23,-0.4472135954999579*cMr[7]); 
  gkyl_mat_set(A,5,16,-0.5*cMr[5]); 
  gkyl_mat_set(A,5,17,-0.5000000000000001*cMr[7]); 
  gkyl_mat_set(A,5,18,-0.4472135954999579*cMr[2]); 
  gkyl_mat_set(A,5,19,-0.4472135954999579*cMr[3]); 
  gkyl_mat_set(A,5,21,(-0.31943828249997*cMr[5])-0.5*cMr[0]); 
  gkyl_mat_set(A,5,22,-0.4472135954999579*cMr[6]); 
  gkyl_mat_set(A,5,23,(-0.31943828249997*cMr[7])-0.5000000000000001*cMr[1]); 
  gkyl_mat_set(A,6,16,-0.5*cMr[6]); 
  gkyl_mat_set(A,6,17,-0.447213595499958*cMr[3]); 
  gkyl_mat_set(A,6,18,-0.5000000000000001*cMr[4]); 
  gkyl_mat_set(A,6,19,(-0.4*cMr[7])-0.447213595499958*cMr[1]); 
  gkyl_mat_set(A,6,20,(-0.31943828249997*cMr[6])-0.5000000000000001*cMr[2]); 
  gkyl_mat_set(A,6,21,-0.4472135954999579*cMr[6]); 
  gkyl_mat_set(A,6,22,(-0.4472135954999579*cMr[5])-0.31943828249997*cMr[4]-0.5*cMr[0]); 
  gkyl_mat_set(A,6,23,-0.4*cMr[3]); 
  gkyl_mat_set(A,7,16,-0.5*cMr[7]); 
  gkyl_mat_set(A,7,17,-0.5000000000000001*cMr[5]); 
  gkyl_mat_set(A,7,18,-0.447213595499958*cMr[3]); 
  gkyl_mat_set(A,7,19,(-0.4*cMr[6])-0.447213595499958*cMr[2]); 
  gkyl_mat_set(A,7,20,-0.4472135954999579*cMr[7]); 
  gkyl_mat_set(A,7,21,(-0.31943828249997*cMr[7])-0.5000000000000001*cMr[1]); 
  gkyl_mat_set(A,7,22,-0.4*cMr[3]); 
  gkyl_mat_set(A,7,23,(-0.31943828249997*cMr[5])-0.4472135954999579*cMr[4]-0.5*cMr[0]); 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  gkyl_mat_set(A,16,0,0.5*m1r[0]); 
  gkyl_mat_set(A,16,1,0.5*m1r[1]); 
  gkyl_mat_set(A,16,2,0.5*m1r[2]); 
  gkyl_mat_set(A,16,3,0.5*m1r[3]); 
  gkyl_mat_set(A,16,4,0.5*m1r[4]); 
  gkyl_mat_set(A,16,5,0.5*m1r[5]); 
  gkyl_mat_set(A,16,6,0.5*m1r[6]); 
  gkyl_mat_set(A,16,7,0.5*m1r[7]); 
  gkyl_mat_set(A,17,0,0.5*m1r[1]); 
  gkyl_mat_set(A,17,1,0.4472135954999579*m1r[4]+0.5*m1r[0]); 
  gkyl_mat_set(A,17,2,0.5*m1r[3]); 
  gkyl_mat_set(A,17,3,0.447213595499958*m1r[6]+0.5*m1r[2]); 
  gkyl_mat_set(A,17,4,0.4472135954999579*m1r[1]); 
  gkyl_mat_set(A,17,5,0.5000000000000001*m1r[7]); 
  gkyl_mat_set(A,17,6,0.447213595499958*m1r[3]); 
  gkyl_mat_set(A,17,7,0.5000000000000001*m1r[5]); 
  gkyl_mat_set(A,18,0,0.5*m1r[2]); 
  gkyl_mat_set(A,18,1,0.5*m1r[3]); 
  gkyl_mat_set(A,18,2,0.4472135954999579*m1r[5]+0.5*m1r[0]); 
  gkyl_mat_set(A,18,3,0.447213595499958*m1r[7]+0.5*m1r[1]); 
  gkyl_mat_set(A,18,4,0.5000000000000001*m1r[6]); 
  gkyl_mat_set(A,18,5,0.4472135954999579*m1r[2]); 
  gkyl_mat_set(A,18,6,0.5000000000000001*m1r[4]); 
  gkyl_mat_set(A,18,7,0.447213595499958*m1r[3]); 
  gkyl_mat_set(A,19,0,0.5*m1r[3]); 
  gkyl_mat_set(A,19,1,0.447213595499958*m1r[6]+0.5*m1r[2]); 
  gkyl_mat_set(A,19,2,0.447213595499958*m1r[7]+0.5*m1r[1]); 
  gkyl_mat_set(A,19,3,0.4472135954999579*m1r[5]+0.4472135954999579*m1r[4]+0.5*m1r[0]); 
  gkyl_mat_set(A,19,4,0.4472135954999579*m1r[3]); 
  gkyl_mat_set(A,19,5,0.4472135954999579*m1r[3]); 
  gkyl_mat_set(A,19,6,0.4*m1r[7]+0.447213595499958*m1r[1]); 
  gkyl_mat_set(A,19,7,0.4*m1r[6]+0.447213595499958*m1r[2]); 
  gkyl_mat_set(A,20,0,0.5*m1r[4]); 
  gkyl_mat_set(A,20,1,0.4472135954999579*m1r[1]); 
  gkyl_mat_set(A,20,2,0.5000000000000001*m1r[6]); 
  gkyl_mat_set(A,20,3,0.4472135954999579*m1r[3]); 
  gkyl_mat_set(A,20,4,0.31943828249997*m1r[4]+0.5*m1r[0]); 
  gkyl_mat_set(A,20,6,0.31943828249997*m1r[6]+0.5000000000000001*m1r[2]); 
  gkyl_mat_set(A,20,7,0.4472135954999579*m1r[7]); 
  gkyl_mat_set(A,21,0,0.5*m1r[5]); 
  gkyl_mat_set(A,21,1,0.5000000000000001*m1r[7]); 
  gkyl_mat_set(A,21,2,0.4472135954999579*m1r[2]); 
  gkyl_mat_set(A,21,3,0.4472135954999579*m1r[3]); 
  gkyl_mat_set(A,21,5,0.31943828249997*m1r[5]+0.5*m1r[0]); 
  gkyl_mat_set(A,21,6,0.4472135954999579*m1r[6]); 
  gkyl_mat_set(A,21,7,0.31943828249997*m1r[7]+0.5000000000000001*m1r[1]); 
  gkyl_mat_set(A,22,0,0.5*m1r[6]); 
  gkyl_mat_set(A,22,1,0.447213595499958*m1r[3]); 
  gkyl_mat_set(A,22,2,0.5000000000000001*m1r[4]); 
  gkyl_mat_set(A,22,3,0.4*m1r[7]+0.447213595499958*m1r[1]); 
  gkyl_mat_set(A,22,4,0.31943828249997*m1r[6]+0.5000000000000001*m1r[2]); 
  gkyl_mat_set(A,22,5,0.4472135954999579*m1r[6]); 
  gkyl_mat_set(A,22,6,0.4472135954999579*m1r[5]+0.31943828249997*m1r[4]+0.5*m1r[0]); 
  gkyl_mat_set(A,22,7,0.4*m1r[3]); 
  gkyl_mat_set(A,23,0,0.5*m1r[7]); 
  gkyl_mat_set(A,23,1,0.5000000000000001*m1r[5]); 
  gkyl_mat_set(A,23,2,0.447213595499958*m1r[3]); 
  gkyl_mat_set(A,23,3,0.4*m1r[6]+0.447213595499958*m1r[2]); 
  gkyl_mat_set(A,23,4,0.4472135954999579*m1r[7]); 
  gkyl_mat_set(A,23,5,0.31943828249997*m1r[7]+0.5000000000000001*m1r[1]); 
  gkyl_mat_set(A,23,6,0.4*m1r[3]); 
  gkyl_mat_set(A,23,7,0.31943828249997*m1r[5]+0.4472135954999579*m1r[4]+0.5*m1r[0]); 
 
  // ....... Block from weak multiply of uY and m0  .......... // 
  gkyl_mat_set(A,8,8,0.5*m0r[0]); 
  gkyl_mat_set(A,8,9,0.5*m0r[1]); 
  gkyl_mat_set(A,8,10,0.5*m0r[2]); 
  gkyl_mat_set(A,8,11,0.5*m0r[3]); 
  gkyl_mat_set(A,8,12,0.5*m0r[4]); 
  gkyl_mat_set(A,8,13,0.5*m0r[5]); 
  gkyl_mat_set(A,8,14,0.5*m0r[6]); 
  gkyl_mat_set(A,8,15,0.5*m0r[7]); 
  gkyl_mat_set(A,9,8,0.5*m0r[1]); 
  gkyl_mat_set(A,9,9,0.4472135954999579*m0r[4]+0.5*m0r[0]); 
  gkyl_mat_set(A,9,10,0.5*m0r[3]); 
  gkyl_mat_set(A,9,11,0.447213595499958*m0r[6]+0.5*m0r[2]); 
  gkyl_mat_set(A,9,12,0.4472135954999579*m0r[1]); 
  gkyl_mat_set(A,9,13,0.5000000000000001*m0r[7]); 
  gkyl_mat_set(A,9,14,0.447213595499958*m0r[3]); 
  gkyl_mat_set(A,9,15,0.5000000000000001*m0r[5]); 
  gkyl_mat_set(A,10,8,0.5*m0r[2]); 
  gkyl_mat_set(A,10,9,0.5*m0r[3]); 
  gkyl_mat_set(A,10,10,0.4472135954999579*m0r[5]+0.5*m0r[0]); 
  gkyl_mat_set(A,10,11,0.447213595499958*m0r[7]+0.5*m0r[1]); 
  gkyl_mat_set(A,10,12,0.5000000000000001*m0r[6]); 
  gkyl_mat_set(A,10,13,0.4472135954999579*m0r[2]); 
  gkyl_mat_set(A,10,14,0.5000000000000001*m0r[4]); 
  gkyl_mat_set(A,10,15,0.447213595499958*m0r[3]); 
  gkyl_mat_set(A,11,8,0.5*m0r[3]); 
  gkyl_mat_set(A,11,9,0.447213595499958*m0r[6]+0.5*m0r[2]); 
  gkyl_mat_set(A,11,10,0.447213595499958*m0r[7]+0.5*m0r[1]); 
  gkyl_mat_set(A,11,11,0.4472135954999579*m0r[5]+0.4472135954999579*m0r[4]+0.5*m0r[0]); 
  gkyl_mat_set(A,11,12,0.4472135954999579*m0r[3]); 
  gkyl_mat_set(A,11,13,0.4472135954999579*m0r[3]); 
  gkyl_mat_set(A,11,14,0.4*m0r[7]+0.447213595499958*m0r[1]); 
  gkyl_mat_set(A,11,15,0.4*m0r[6]+0.447213595499958*m0r[2]); 
  gkyl_mat_set(A,12,8,0.5*m0r[4]); 
  gkyl_mat_set(A,12,9,0.4472135954999579*m0r[1]); 
  gkyl_mat_set(A,12,10,0.5000000000000001*m0r[6]); 
  gkyl_mat_set(A,12,11,0.4472135954999579*m0r[3]); 
  gkyl_mat_set(A,12,12,0.31943828249997*m0r[4]+0.5*m0r[0]); 
  gkyl_mat_set(A,12,14,0.31943828249997*m0r[6]+0.5000000000000001*m0r[2]); 
  gkyl_mat_set(A,12,15,0.4472135954999579*m0r[7]); 
  gkyl_mat_set(A,13,8,0.5*m0r[5]); 
  gkyl_mat_set(A,13,9,0.5000000000000001*m0r[7]); 
  gkyl_mat_set(A,13,10,0.4472135954999579*m0r[2]); 
  gkyl_mat_set(A,13,11,0.4472135954999579*m0r[3]); 
  gkyl_mat_set(A,13,13,0.31943828249997*m0r[5]+0.5*m0r[0]); 
  gkyl_mat_set(A,13,14,0.4472135954999579*m0r[6]); 
  gkyl_mat_set(A,13,15,0.31943828249997*m0r[7]+0.5000000000000001*m0r[1]); 
  gkyl_mat_set(A,14,8,0.5*m0r[6]); 
  gkyl_mat_set(A,14,9,0.447213595499958*m0r[3]); 
  gkyl_mat_set(A,14,10,0.5000000000000001*m0r[4]); 
  gkyl_mat_set(A,14,11,0.4*m0r[7]+0.447213595499958*m0r[1]); 
  gkyl_mat_set(A,14,12,0.31943828249997*m0r[6]+0.5000000000000001*m0r[2]); 
  gkyl_mat_set(A,14,13,0.4472135954999579*m0r[6]); 
  gkyl_mat_set(A,14,14,0.4472135954999579*m0r[5]+0.31943828249997*m0r[4]+0.5*m0r[0]); 
  gkyl_mat_set(A,14,15,0.4*m0r[3]); 
  gkyl_mat_set(A,15,8,0.5*m0r[7]); 
  gkyl_mat_set(A,15,9,0.5000000000000001*m0r[5]); 
  gkyl_mat_set(A,15,10,0.447213595499958*m0r[3]); 
  gkyl_mat_set(A,15,11,0.4*m0r[6]+0.447213595499958*m0r[2]); 
  gkyl_mat_set(A,15,12,0.4472135954999579*m0r[7]); 
  gkyl_mat_set(A,15,13,0.31943828249997*m0r[7]+0.5000000000000001*m0r[1]); 
  gkyl_mat_set(A,15,14,0.4*m0r[3]); 
  gkyl_mat_set(A,15,15,0.31943828249997*m0r[5]+0.4472135954999579*m0r[4]+0.5*m0r[0]); 
 
  // ....... Block from correction to uY .......... // 
  gkyl_mat_set(A,8,16,-0.5*cMr[8]); 
  gkyl_mat_set(A,8,17,-0.5*cMr[9]); 
  gkyl_mat_set(A,8,18,-0.5*cMr[10]); 
  gkyl_mat_set(A,8,19,-0.5*cMr[11]); 
  gkyl_mat_set(A,8,20,-0.5*cMr[12]); 
  gkyl_mat_set(A,8,21,-0.5*cMr[13]); 
  gkyl_mat_set(A,8,22,-0.5*cMr[14]); 
  gkyl_mat_set(A,8,23,-0.5*cMr[15]); 
  gkyl_mat_set(A,9,16,-0.5*cMr[9]); 
  gkyl_mat_set(A,9,17,(-0.4472135954999579*cMr[12])-0.5*cMr[8]); 
  gkyl_mat_set(A,9,18,-0.5*cMr[11]); 
  gkyl_mat_set(A,9,19,(-0.447213595499958*cMr[14])-0.5*cMr[10]); 
  gkyl_mat_set(A,9,20,-0.4472135954999579*cMr[9]); 
  gkyl_mat_set(A,9,21,-0.5000000000000001*cMr[15]); 
  gkyl_mat_set(A,9,22,-0.447213595499958*cMr[11]); 
  gkyl_mat_set(A,9,23,-0.5000000000000001*cMr[13]); 
  gkyl_mat_set(A,10,16,-0.5*cMr[10]); 
  gkyl_mat_set(A,10,17,-0.5*cMr[11]); 
  gkyl_mat_set(A,10,18,(-0.4472135954999579*cMr[13])-0.5*cMr[8]); 
  gkyl_mat_set(A,10,19,(-0.447213595499958*cMr[15])-0.5*cMr[9]); 
  gkyl_mat_set(A,10,20,-0.5000000000000001*cMr[14]); 
  gkyl_mat_set(A,10,21,-0.4472135954999579*cMr[10]); 
  gkyl_mat_set(A,10,22,-0.5000000000000001*cMr[12]); 
  gkyl_mat_set(A,10,23,-0.447213595499958*cMr[11]); 
  gkyl_mat_set(A,11,16,-0.5*cMr[11]); 
  gkyl_mat_set(A,11,17,(-0.447213595499958*cMr[14])-0.5*cMr[10]); 
  gkyl_mat_set(A,11,18,(-0.447213595499958*cMr[15])-0.5*cMr[9]); 
  gkyl_mat_set(A,11,19,(-0.4472135954999579*cMr[13])-0.4472135954999579*cMr[12]-0.5*cMr[8]); 
  gkyl_mat_set(A,11,20,-0.4472135954999579*cMr[11]); 
  gkyl_mat_set(A,11,21,-0.4472135954999579*cMr[11]); 
  gkyl_mat_set(A,11,22,(-0.4*cMr[15])-0.447213595499958*cMr[9]); 
  gkyl_mat_set(A,11,23,(-0.4*cMr[14])-0.447213595499958*cMr[10]); 
  gkyl_mat_set(A,12,16,-0.5*cMr[12]); 
  gkyl_mat_set(A,12,17,-0.4472135954999579*cMr[9]); 
  gkyl_mat_set(A,12,18,-0.5000000000000001*cMr[14]); 
  gkyl_mat_set(A,12,19,-0.4472135954999579*cMr[11]); 
  gkyl_mat_set(A,12,20,(-0.31943828249997*cMr[12])-0.5*cMr[8]); 
  gkyl_mat_set(A,12,22,(-0.31943828249997*cMr[14])-0.5000000000000001*cMr[10]); 
  gkyl_mat_set(A,12,23,-0.4472135954999579*cMr[15]); 
  gkyl_mat_set(A,13,16,-0.5*cMr[13]); 
  gkyl_mat_set(A,13,17,-0.5000000000000001*cMr[15]); 
  gkyl_mat_set(A,13,18,-0.4472135954999579*cMr[10]); 
  gkyl_mat_set(A,13,19,-0.4472135954999579*cMr[11]); 
  gkyl_mat_set(A,13,21,(-0.31943828249997*cMr[13])-0.5*cMr[8]); 
  gkyl_mat_set(A,13,22,-0.4472135954999579*cMr[14]); 
  gkyl_mat_set(A,13,23,(-0.31943828249997*cMr[15])-0.5000000000000001*cMr[9]); 
  gkyl_mat_set(A,14,16,-0.5*cMr[14]); 
  gkyl_mat_set(A,14,17,-0.447213595499958*cMr[11]); 
  gkyl_mat_set(A,14,18,-0.5000000000000001*cMr[12]); 
  gkyl_mat_set(A,14,19,(-0.4*cMr[15])-0.447213595499958*cMr[9]); 
  gkyl_mat_set(A,14,20,(-0.31943828249997*cMr[14])-0.5000000000000001*cMr[10]); 
  gkyl_mat_set(A,14,21,-0.4472135954999579*cMr[14]); 
  gkyl_mat_set(A,14,22,(-0.4472135954999579*cMr[13])-0.31943828249997*cMr[12]-0.5*cMr[8]); 
  gkyl_mat_set(A,14,23,-0.4*cMr[11]); 
  gkyl_mat_set(A,15,16,-0.5*cMr[15]); 
  gkyl_mat_set(A,15,17,-0.5000000000000001*cMr[13]); 
  gkyl_mat_set(A,15,18,-0.447213595499958*cMr[11]); 
  gkyl_mat_set(A,15,19,(-0.4*cMr[14])-0.447213595499958*cMr[10]); 
  gkyl_mat_set(A,15,20,-0.4472135954999579*cMr[15]); 
  gkyl_mat_set(A,15,21,(-0.31943828249997*cMr[15])-0.5000000000000001*cMr[9]); 
  gkyl_mat_set(A,15,22,-0.4*cMr[11]); 
  gkyl_mat_set(A,15,23,(-0.31943828249997*cMr[13])-0.4472135954999579*cMr[12]-0.5*cMr[8]); 
 
  // ....... Block from weak multiply of uY and m1Y  .......... // 
  gkyl_mat_set(A,16,8,0.5*m1r[8]); 
  gkyl_mat_set(A,16,9,0.5*m1r[9]); 
  gkyl_mat_set(A,16,10,0.5*m1r[10]); 
  gkyl_mat_set(A,16,11,0.5*m1r[11]); 
  gkyl_mat_set(A,16,12,0.5*m1r[12]); 
  gkyl_mat_set(A,16,13,0.5*m1r[13]); 
  gkyl_mat_set(A,16,14,0.5*m1r[14]); 
  gkyl_mat_set(A,16,15,0.5*m1r[15]); 
  gkyl_mat_set(A,17,8,0.5*m1r[9]); 
  gkyl_mat_set(A,17,9,0.4472135954999579*m1r[12]+0.5*m1r[8]); 
  gkyl_mat_set(A,17,10,0.5*m1r[11]); 
  gkyl_mat_set(A,17,11,0.447213595499958*m1r[14]+0.5*m1r[10]); 
  gkyl_mat_set(A,17,12,0.4472135954999579*m1r[9]); 
  gkyl_mat_set(A,17,13,0.5000000000000001*m1r[15]); 
  gkyl_mat_set(A,17,14,0.447213595499958*m1r[11]); 
  gkyl_mat_set(A,17,15,0.5000000000000001*m1r[13]); 
  gkyl_mat_set(A,18,8,0.5*m1r[10]); 
  gkyl_mat_set(A,18,9,0.5*m1r[11]); 
  gkyl_mat_set(A,18,10,0.4472135954999579*m1r[13]+0.5*m1r[8]); 
  gkyl_mat_set(A,18,11,0.447213595499958*m1r[15]+0.5*m1r[9]); 
  gkyl_mat_set(A,18,12,0.5000000000000001*m1r[14]); 
  gkyl_mat_set(A,18,13,0.4472135954999579*m1r[10]); 
  gkyl_mat_set(A,18,14,0.5000000000000001*m1r[12]); 
  gkyl_mat_set(A,18,15,0.447213595499958*m1r[11]); 
  gkyl_mat_set(A,19,8,0.5*m1r[11]); 
  gkyl_mat_set(A,19,9,0.447213595499958*m1r[14]+0.5*m1r[10]); 
  gkyl_mat_set(A,19,10,0.447213595499958*m1r[15]+0.5*m1r[9]); 
  gkyl_mat_set(A,19,11,0.4472135954999579*m1r[13]+0.4472135954999579*m1r[12]+0.5*m1r[8]); 
  gkyl_mat_set(A,19,12,0.4472135954999579*m1r[11]); 
  gkyl_mat_set(A,19,13,0.4472135954999579*m1r[11]); 
  gkyl_mat_set(A,19,14,0.4*m1r[15]+0.447213595499958*m1r[9]); 
  gkyl_mat_set(A,19,15,0.4*m1r[14]+0.447213595499958*m1r[10]); 
  gkyl_mat_set(A,20,8,0.5*m1r[12]); 
  gkyl_mat_set(A,20,9,0.4472135954999579*m1r[9]); 
  gkyl_mat_set(A,20,10,0.5000000000000001*m1r[14]); 
  gkyl_mat_set(A,20,11,0.4472135954999579*m1r[11]); 
  gkyl_mat_set(A,20,12,0.31943828249997*m1r[12]+0.5*m1r[8]); 
  gkyl_mat_set(A,20,14,0.31943828249997*m1r[14]+0.5000000000000001*m1r[10]); 
  gkyl_mat_set(A,20,15,0.4472135954999579*m1r[15]); 
  gkyl_mat_set(A,21,8,0.5*m1r[13]); 
  gkyl_mat_set(A,21,9,0.5000000000000001*m1r[15]); 
  gkyl_mat_set(A,21,10,0.4472135954999579*m1r[10]); 
  gkyl_mat_set(A,21,11,0.4472135954999579*m1r[11]); 
  gkyl_mat_set(A,21,13,0.31943828249997*m1r[13]+0.5*m1r[8]); 
  gkyl_mat_set(A,21,14,0.4472135954999579*m1r[14]); 
  gkyl_mat_set(A,21,15,0.31943828249997*m1r[15]+0.5000000000000001*m1r[9]); 
  gkyl_mat_set(A,22,8,0.5*m1r[14]); 
  gkyl_mat_set(A,22,9,0.447213595499958*m1r[11]); 
  gkyl_mat_set(A,22,10,0.5000000000000001*m1r[12]); 
  gkyl_mat_set(A,22,11,0.4*m1r[15]+0.447213595499958*m1r[9]); 
  gkyl_mat_set(A,22,12,0.31943828249997*m1r[14]+0.5000000000000001*m1r[10]); 
  gkyl_mat_set(A,22,13,0.4472135954999579*m1r[14]); 
  gkyl_mat_set(A,22,14,0.4472135954999579*m1r[13]+0.31943828249997*m1r[12]+0.5*m1r[8]); 
  gkyl_mat_set(A,22,15,0.4*m1r[11]); 
  gkyl_mat_set(A,23,8,0.5*m1r[15]); 
  gkyl_mat_set(A,23,9,0.5000000000000001*m1r[13]); 
  gkyl_mat_set(A,23,10,0.447213595499958*m1r[11]); 
  gkyl_mat_set(A,23,11,0.4*m1r[14]+0.447213595499958*m1r[10]); 
  gkyl_mat_set(A,23,12,0.4472135954999579*m1r[15]); 
  gkyl_mat_set(A,23,13,0.31943828249997*m1r[15]+0.5000000000000001*m1r[9]); 
  gkyl_mat_set(A,23,14,0.4*m1r[11]); 
  gkyl_mat_set(A,23,15,0.31943828249997*m1r[13]+0.4472135954999579*m1r[12]+0.5*m1r[8]); 
 
  // ....... Block from correction to vtSq .......... // 
  gkyl_mat_set(A,16,16,m0r[0]-0.5*cEr[0]); 
  gkyl_mat_set(A,16,17,m0r[1]-0.5*cEr[1]); 
  gkyl_mat_set(A,16,18,m0r[2]-0.5*cEr[2]); 
  gkyl_mat_set(A,16,19,m0r[3]-0.5*cEr[3]); 
  gkyl_mat_set(A,16,20,m0r[4]-0.5*cEr[4]); 
  gkyl_mat_set(A,16,21,m0r[5]-0.5*cEr[5]); 
  gkyl_mat_set(A,16,22,m0r[6]-0.5*cEr[6]); 
  gkyl_mat_set(A,16,23,m0r[7]-0.5*cEr[7]); 
  gkyl_mat_set(A,17,16,m0r[1]-0.5*cEr[1]); 
  gkyl_mat_set(A,17,17,0.8944271909999159*m0r[4]-0.4472135954999579*cEr[4]+m0r[0]-0.5*cEr[0]); 
  gkyl_mat_set(A,17,18,m0r[3]-0.5*cEr[3]); 
  gkyl_mat_set(A,17,19,0.8944271909999161*m0r[6]-0.447213595499958*cEr[6]+m0r[2]-0.5*cEr[2]); 
  gkyl_mat_set(A,17,20,0.8944271909999159*m0r[1]-0.4472135954999579*cEr[1]); 
  gkyl_mat_set(A,17,21,1.0*m0r[7]-0.5000000000000001*cEr[7]); 
  gkyl_mat_set(A,17,22,0.8944271909999161*m0r[3]-0.447213595499958*cEr[3]); 
  gkyl_mat_set(A,17,23,1.0*m0r[5]-0.5000000000000001*cEr[5]); 
  gkyl_mat_set(A,18,16,m0r[2]-0.5*cEr[2]); 
  gkyl_mat_set(A,18,17,m0r[3]-0.5*cEr[3]); 
  gkyl_mat_set(A,18,18,0.8944271909999159*m0r[5]-0.4472135954999579*cEr[5]+m0r[0]-0.5*cEr[0]); 
  gkyl_mat_set(A,18,19,0.8944271909999161*m0r[7]-0.447213595499958*cEr[7]+m0r[1]-0.5*cEr[1]); 
  gkyl_mat_set(A,18,20,1.0*m0r[6]-0.5000000000000001*cEr[6]); 
  gkyl_mat_set(A,18,21,0.8944271909999159*m0r[2]-0.4472135954999579*cEr[2]); 
  gkyl_mat_set(A,18,22,1.0*m0r[4]-0.5000000000000001*cEr[4]); 
  gkyl_mat_set(A,18,23,0.8944271909999161*m0r[3]-0.447213595499958*cEr[3]); 
  gkyl_mat_set(A,19,16,m0r[3]-0.5*cEr[3]); 
  gkyl_mat_set(A,19,17,0.8944271909999161*m0r[6]-0.447213595499958*cEr[6]+m0r[2]-0.5*cEr[2]); 
  gkyl_mat_set(A,19,18,0.8944271909999161*m0r[7]-0.447213595499958*cEr[7]+m0r[1]-0.5*cEr[1]); 
  gkyl_mat_set(A,19,19,0.8944271909999159*m0r[5]-0.4472135954999579*cEr[5]+0.8944271909999159*m0r[4]-0.4472135954999579*cEr[4]+m0r[0]-0.5*cEr[0]); 
  gkyl_mat_set(A,19,20,0.8944271909999159*m0r[3]-0.4472135954999579*cEr[3]); 
  gkyl_mat_set(A,19,21,0.8944271909999159*m0r[3]-0.4472135954999579*cEr[3]); 
  gkyl_mat_set(A,19,22,0.8*m0r[7]-0.4*cEr[7]+0.8944271909999161*m0r[1]-0.447213595499958*cEr[1]); 
  gkyl_mat_set(A,19,23,0.8*m0r[6]-0.4*cEr[6]+0.8944271909999161*m0r[2]-0.447213595499958*cEr[2]); 
  gkyl_mat_set(A,20,16,m0r[4]-0.5*cEr[4]); 
  gkyl_mat_set(A,20,17,0.8944271909999159*m0r[1]-0.4472135954999579*cEr[1]); 
  gkyl_mat_set(A,20,18,1.0*m0r[6]-0.5000000000000001*cEr[6]); 
  gkyl_mat_set(A,20,19,0.8944271909999159*m0r[3]-0.4472135954999579*cEr[3]); 
  gkyl_mat_set(A,20,20,0.6388765649999399*m0r[4]-0.31943828249997*cEr[4]+m0r[0]-0.5*cEr[0]); 
  gkyl_mat_set(A,20,22,0.6388765649999399*m0r[6]-0.31943828249997*cEr[6]+1.0*m0r[2]-0.5000000000000001*cEr[2]); 
  gkyl_mat_set(A,20,23,0.8944271909999159*m0r[7]-0.4472135954999579*cEr[7]); 
  gkyl_mat_set(A,21,16,m0r[5]-0.5*cEr[5]); 
  gkyl_mat_set(A,21,17,1.0*m0r[7]-0.5000000000000001*cEr[7]); 
  gkyl_mat_set(A,21,18,0.8944271909999159*m0r[2]-0.4472135954999579*cEr[2]); 
  gkyl_mat_set(A,21,19,0.8944271909999159*m0r[3]-0.4472135954999579*cEr[3]); 
  gkyl_mat_set(A,21,21,0.6388765649999399*m0r[5]-0.31943828249997*cEr[5]+m0r[0]-0.5*cEr[0]); 
  gkyl_mat_set(A,21,22,0.8944271909999159*m0r[6]-0.4472135954999579*cEr[6]); 
  gkyl_mat_set(A,21,23,0.6388765649999399*m0r[7]-0.31943828249997*cEr[7]+1.0*m0r[1]-0.5000000000000001*cEr[1]); 
  gkyl_mat_set(A,22,16,m0r[6]-0.5*cEr[6]); 
  gkyl_mat_set(A,22,17,0.8944271909999161*m0r[3]-0.447213595499958*cEr[3]); 
  gkyl_mat_set(A,22,18,1.0*m0r[4]-0.5000000000000001*cEr[4]); 
  gkyl_mat_set(A,22,19,0.8*m0r[7]-0.4*cEr[7]+0.8944271909999161*m0r[1]-0.447213595499958*cEr[1]); 
  gkyl_mat_set(A,22,20,0.6388765649999399*m0r[6]-0.31943828249997*cEr[6]+1.0*m0r[2]-0.5000000000000001*cEr[2]); 
  gkyl_mat_set(A,22,21,0.8944271909999159*m0r[6]-0.4472135954999579*cEr[6]); 
  gkyl_mat_set(A,22,22,0.8944271909999159*m0r[5]-0.4472135954999579*cEr[5]+0.6388765649999399*m0r[4]-0.31943828249997*cEr[4]+m0r[0]-0.5*cEr[0]); 
  gkyl_mat_set(A,22,23,0.8*m0r[3]-0.4*cEr[3]); 
  gkyl_mat_set(A,23,16,m0r[7]-0.5*cEr[7]); 
  gkyl_mat_set(A,23,17,1.0*m0r[5]-0.5000000000000001*cEr[5]); 
  gkyl_mat_set(A,23,18,0.8944271909999161*m0r[3]-0.447213595499958*cEr[3]); 
  gkyl_mat_set(A,23,19,0.8*m0r[6]-0.4*cEr[6]+0.8944271909999161*m0r[2]-0.447213595499958*cEr[2]); 
  gkyl_mat_set(A,23,20,0.8944271909999159*m0r[7]-0.4472135954999579*cEr[7]); 
  gkyl_mat_set(A,23,21,0.6388765649999399*m0r[7]-0.31943828249997*cEr[7]+1.0*m0r[1]-0.5000000000000001*cEr[1]); 
  gkyl_mat_set(A,23,22,0.8*m0r[3]-0.4*cEr[3]); 
  gkyl_mat_set(A,23,23,0.6388765649999399*m0r[5]-0.31943828249997*cEr[5]+0.8944271909999159*m0r[4]-0.4472135954999579*cEr[4]+m0r[0]-0.5*cEr[0]); 
 
} 
 
