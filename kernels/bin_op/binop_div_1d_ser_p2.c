#include <gkyl_mat.h> 
#include <gkyl_binop_div_ser.h> 
 
void binop_div_1d_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *f, const double *g, double *fdivg) 
{ 
  // A:       preallocated LHS matrix. 
  // rhs:     preallocated RHS vector. 
  // f:       numerator field (must be a scalar). 
  // g:       denominator field (must be a scalar). 
  // fdivg:   output field. 
 
  // If a corner value is below zero, use cell average g.
  bool avgg = false;
  if (1.581138830084189*g[2]-1.224744871391589*g[1]+0.7071067811865475*g[0] < 0.0) { 
    avgg = true;
  }
  if (1.581138830084189*g[2]+1.224744871391589*g[1]+0.7071067811865475*g[0] < 0.0) { 
    avgg = true;
  }
 
  double lhs[3]; 
  if (avgg) { 
    lhs[0] = g[0]; 
    lhs[1] = 0.0; 
    lhs[2] = 0.0; 
    gkyl_mat_set(rhs,1,0,f[0]); 
    gkyl_mat_set(rhs,2,0,0.0); 
    gkyl_mat_set(rhs,3,0,0.0); 
  } else { 
    lhs[0] = g[0]; 
    lhs[1] = g[1]; 
    lhs[2] = g[2]; 
    gkyl_mat_set(rhs,1,0,f[0]); 
    gkyl_mat_set(rhs,2,0,f[1]); 
    gkyl_mat_set(rhs,3,0,f[2]); 
  } 
 
  // Fill LHS matrix. 
  gkyl_mat_set(A,0,0,0.7071067811865475*lhs[0]); 
  gkyl_mat_set(A,0,1,0.7071067811865475*lhs[1]); 
  gkyl_mat_set(A,0,2,0.7071067811865475*lhs[2]); 
  gkyl_mat_set(A,1,0,0.7071067811865475*lhs[1]); 
  gkyl_mat_set(A,1,1,0.6324555320336759*lhs[2]+0.7071067811865475*lhs[0]); 
  gkyl_mat_set(A,1,2,0.6324555320336759*lhs[1]); 
  gkyl_mat_set(A,2,0,0.7071067811865475*lhs[2]); 
  gkyl_mat_set(A,2,1,0.6324555320336759*lhs[1]); 
  gkyl_mat_set(A,2,2,0.4517539514526256*lhs[2]+0.7071067811865475*lhs[0]); 

  // Solve the system of equations. 
  long ipiv[3]; 
  gkyl_mat_linsolve_lu(A,rhs,ipiv); 
  for(size_t i=0; i<3; i++) 
  { 
    fdivg[i] = gkyl_mat_get(rhs,i,0); 
  } 
} 
