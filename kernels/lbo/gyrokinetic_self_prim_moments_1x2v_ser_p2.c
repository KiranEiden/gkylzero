#include <gkyl_prim_lbo_gyrokinetic_kernels.h> 
 
GKYL_CU_DH void gyrokinetic_self_prim_moments_1x2v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *moms, const double *boundary_corrections) 
{ 
  // A:                    Matrix to be inverted to solve Ax = rhs (set by this function). 
  // rhs:                  right-hand side of Ax = rhs (set by this function). 
  // moms:                 moments of the distribution function (Zeroth, First, and Second in single array). 
  // boundary_corrections: boundary corrections to u and vtSq. 
 
  // If m0 or m2 is below zero at a corner, use cell averages.
  bool notCellAvg = true;
  if (notCellAvg && (0.7071067811865475*(2.23606797749979*moms[2]-1.732050807568877*moms[1]+moms[0]) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.7071067811865475*(2.23606797749979*moms[2]+1.732050807568877*moms[1]+moms[0]) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.7071067811865475*(2.23606797749979*moms[8]-1.732050807568877*moms[7]+moms[6]) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.7071067811865475*(2.23606797749979*moms[8]+1.732050807568877*moms[7]+moms[6]) < 0)) notCellAvg = false; 
 
  double m0r[3] = {0.0}; 
  double m1r[3] = {0.0}; 
  double cMr[3] = {0.0}; 
  double cEr[3] = {0.0}; 
  if (notCellAvg) { 
    m0r[0] = moms[0]; 
    m0r[1] = moms[1]; 
    m0r[2] = moms[2]; 
    m1r[0] = moms[3]; 
    m1r[1] = moms[4]; 
    m1r[2] = moms[5]; 
    gkyl_mat_set(rhs,0,0,moms[3]); 
    gkyl_mat_set(rhs,1,0,moms[4]); 
    gkyl_mat_set(rhs,2,0,moms[5]); 
    cMr[0] = boundary_corrections[0]; 
    cMr[1] = boundary_corrections[1]; 
    cMr[2] = boundary_corrections[2]; 
    cEr[0] = boundary_corrections[3]; 
    cEr[1] = boundary_corrections[4]; 
    cEr[2] = boundary_corrections[5]; 
    gkyl_mat_set(rhs,3,0,moms[6]); 
    gkyl_mat_set(rhs,4,0,moms[7]); 
    gkyl_mat_set(rhs,5,0,moms[8]); 
  } else { 
    m0r[0] = moms[0]; 
    m0r[1] = 0.0; 
    m0r[2] = 0.0; 
    m1r[0] = moms[3]; 
    m1r[1] = 0.0; 
    m1r[2] = 0.0; 
    gkyl_mat_set(rhs,0,0,moms[3]); 
    gkyl_mat_set(rhs,1,0,0.0); 
    gkyl_mat_set(rhs,2,0,0.0); 
    cMr[0] = boundary_corrections[0]; 
    cMr[1] = 0.0; 
    cMr[2] = 0.0; 
    cEr[0] = boundary_corrections[3]; 
    cEr[1] = 0.0; 
    cEr[2] = 0.0; 
    gkyl_mat_set(rhs,3,0,moms[6]); 
    gkyl_mat_set(rhs,4,0,0.0); 
    gkyl_mat_set(rhs,5,0,0.0); 
  } 
 
  // ....... Block from weak multiply of ux and m0  .......... // 
  gkyl_mat_set(A,0,0,0.7071067811865475*m0r[0]); 
  gkyl_mat_set(A,0,1,0.7071067811865475*m0r[1]); 
  gkyl_mat_set(A,0,2,0.7071067811865475*m0r[2]); 
  gkyl_mat_set(A,1,0,0.7071067811865475*m0r[1]); 
  gkyl_mat_set(A,1,1,0.6324555320336759*m0r[2]+0.7071067811865475*m0r[0]); 
  gkyl_mat_set(A,1,2,0.6324555320336759*m0r[1]); 
  gkyl_mat_set(A,2,0,0.7071067811865475*m0r[2]); 
  gkyl_mat_set(A,2,1,0.6324555320336759*m0r[1]); 
  gkyl_mat_set(A,2,2,0.4517539514526256*m0r[2]+0.7071067811865475*m0r[0]); 
 
  // ....... Block from correction to ux .......... // 
  gkyl_mat_set(A,0,3,-0.7071067811865475*cMr[0]); 
  gkyl_mat_set(A,0,4,-0.7071067811865475*cMr[1]); 
  gkyl_mat_set(A,0,5,-0.7071067811865475*cMr[2]); 
  gkyl_mat_set(A,1,3,-0.7071067811865475*cMr[1]); 
  gkyl_mat_set(A,1,4,(-0.6324555320336759*cMr[2])-0.7071067811865475*cMr[0]); 
  gkyl_mat_set(A,1,5,-0.6324555320336759*cMr[1]); 
  gkyl_mat_set(A,2,3,-0.7071067811865475*cMr[2]); 
  gkyl_mat_set(A,2,4,-0.6324555320336759*cMr[1]); 
  gkyl_mat_set(A,2,5,(-0.4517539514526256*cMr[2])-0.7071067811865475*cMr[0]); 
 
  // ....... Block from weak multiply of ux and m1x  .......... // 
  gkyl_mat_set(A,3,0,0.7071067811865475*m1r[0]); 
  gkyl_mat_set(A,3,1,0.7071067811865475*m1r[1]); 
  gkyl_mat_set(A,3,2,0.7071067811865475*m1r[2]); 
  gkyl_mat_set(A,4,0,0.7071067811865475*m1r[1]); 
  gkyl_mat_set(A,4,1,0.6324555320336759*m1r[2]+0.7071067811865475*m1r[0]); 
  gkyl_mat_set(A,4,2,0.6324555320336759*m1r[1]); 
  gkyl_mat_set(A,5,0,0.7071067811865475*m1r[2]); 
  gkyl_mat_set(A,5,1,0.6324555320336759*m1r[1]); 
  gkyl_mat_set(A,5,2,0.4517539514526256*m1r[2]+0.7071067811865475*m1r[0]); 
 
  // ....... Block from correction to vtSq .......... // 
  gkyl_mat_set(A,3,3,2.121320343559642*m0r[0]-0.7071067811865475*cEr[0]); 
  gkyl_mat_set(A,3,4,2.121320343559642*m0r[1]-0.7071067811865475*cEr[1]); 
  gkyl_mat_set(A,3,5,2.121320343559642*m0r[2]-0.7071067811865475*cEr[2]); 
  gkyl_mat_set(A,4,3,2.121320343559642*m0r[1]-0.7071067811865475*cEr[1]); 
  gkyl_mat_set(A,4,4,1.897366596101028*m0r[2]-0.6324555320336759*cEr[2]+2.121320343559642*m0r[0]-0.7071067811865475*cEr[0]); 
  gkyl_mat_set(A,4,5,1.897366596101028*m0r[1]-0.6324555320336759*cEr[1]); 
  gkyl_mat_set(A,5,3,2.121320343559642*m0r[2]-0.7071067811865475*cEr[2]); 
  gkyl_mat_set(A,5,4,1.897366596101028*m0r[1]-0.6324555320336759*cEr[1]); 
  gkyl_mat_set(A,5,5,1.355261854357877*m0r[2]-0.4517539514526256*cEr[2]+2.121320343559642*m0r[0]-0.7071067811865475*cEr[0]); 
 
} 
 
