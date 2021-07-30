#include <gkyl_mat.h> 
#include <gkyl_binop_div_ser.h> 
 
void binop_div_3d_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *f, const double *g, double *fdivg) 
{ 
  // A:       preallocated LHS matrix. 
  // rhs:     preallocated RHS vector. 
  // f:       numerator field (must be a scalar). 
  // g:       denominator field (must be a scalar). 
  // fdivg:   output field. 
 
  // If a corner value is below zero, use cell average g.
  bool avgg = false;
  if (2.371708245126284*(g[19]+g[18]+g[17])-1.369306393762915*(g[16]+g[15]+g[14]+g[13]+g[12]+g[11])-1.837117307087383*g[10]+0.7905694150420947*(g[9]+g[8]+g[7])+1.060660171779821*(g[6]+g[5]+g[4])-0.6123724356957944*(g[3]+g[2]+g[1])+0.3535533905932737*g[0] < 0.0) { 
    avgg = true;
  }
  if ((-2.371708245126284*(g[19]+g[18]))+2.371708245126284*g[17]-1.369306393762915*g[16]+1.369306393762915*g[15]-1.369306393762915*(g[14]+g[13])+1.369306393762915*g[12]-1.369306393762915*g[11]+1.837117307087383*g[10]+0.7905694150420947*(g[9]+g[8]+g[7])+1.060660171779821*g[6]-1.060660171779821*(g[5]+g[4])-0.6123724356957944*(g[3]+g[2])+0.6123724356957944*g[1]+0.3535533905932737*g[0] < 0.0) { 
    avgg = true;
  }
  if ((-2.371708245126284*g[19])+2.371708245126284*g[18]-2.371708245126284*g[17]+1.369306393762915*g[16]-1.369306393762915*(g[15]+g[14]+g[13]+g[12])+1.369306393762915*g[11]+1.837117307087383*g[10]+0.7905694150420947*(g[9]+g[8]+g[7])-1.060660171779821*g[6]+1.060660171779821*g[5]-1.060660171779821*g[4]-0.6123724356957944*g[3]+0.6123724356957944*g[2]-0.6123724356957944*g[1]+0.3535533905932737*g[0] < 0.0) { 
    avgg = true;
  }
  if (2.371708245126284*g[19]-2.371708245126284*(g[18]+g[17])+1.369306393762915*(g[16]+g[15])-1.369306393762915*(g[14]+g[13])+1.369306393762915*(g[12]+g[11])-1.837117307087383*g[10]+0.7905694150420947*(g[9]+g[8]+g[7])-1.060660171779821*(g[6]+g[5])+1.060660171779821*g[4]-0.6123724356957944*g[3]+0.6123724356957944*(g[2]+g[1])+0.3535533905932737*g[0] < 0.0) { 
    avgg = true;
  }
  if (2.371708245126284*g[19]-2.371708245126284*(g[18]+g[17])-1.369306393762915*(g[16]+g[15])+1.369306393762915*(g[14]+g[13])-1.369306393762915*(g[12]+g[11])+1.837117307087383*g[10]+0.7905694150420947*(g[9]+g[8]+g[7])-1.060660171779821*(g[6]+g[5])+1.060660171779821*g[4]+0.6123724356957944*g[3]-0.6123724356957944*(g[2]+g[1])+0.3535533905932737*g[0] < 0.0) { 
    avgg = true;
  }
  if ((-2.371708245126284*g[19])+2.371708245126284*g[18]-2.371708245126284*g[17]-1.369306393762915*g[16]+1.369306393762915*(g[15]+g[14]+g[13]+g[12])-1.369306393762915*g[11]-1.837117307087383*g[10]+0.7905694150420947*(g[9]+g[8]+g[7])-1.060660171779821*g[6]+1.060660171779821*g[5]-1.060660171779821*g[4]+0.6123724356957944*g[3]-0.6123724356957944*g[2]+0.6123724356957944*g[1]+0.3535533905932737*g[0] < 0.0) { 
    avgg = true;
  }
  if ((-2.371708245126284*(g[19]+g[18]))+2.371708245126284*g[17]+1.369306393762915*g[16]-1.369306393762915*g[15]+1.369306393762915*(g[14]+g[13])-1.369306393762915*g[12]+1.369306393762915*g[11]-1.837117307087383*g[10]+0.7905694150420947*(g[9]+g[8]+g[7])+1.060660171779821*g[6]-1.060660171779821*(g[5]+g[4])+0.6123724356957944*(g[3]+g[2])-0.6123724356957944*g[1]+0.3535533905932737*g[0] < 0.0) { 
    avgg = true;
  }
  if (2.371708245126284*(g[19]+g[18]+g[17])+1.369306393762915*(g[16]+g[15]+g[14]+g[13]+g[12]+g[11])+1.837117307087383*g[10]+0.7905694150420947*(g[9]+g[8]+g[7])+1.060660171779821*(g[6]+g[5]+g[4])+0.6123724356957944*(g[3]+g[2]+g[1])+0.3535533905932737*g[0] < 0.0) { 
    avgg = true;
  }
 
  double lhs[20]; 
  if (avgg) { 
    lhs[0] = g[0]; 
    lhs[1] = 0.0; 
    lhs[2] = 0.0; 
    lhs[3] = 0.0; 
    lhs[4] = 0.0; 
    lhs[5] = 0.0; 
    lhs[6] = 0.0; 
    lhs[7] = 0.0; 
    lhs[8] = 0.0; 
    lhs[9] = 0.0; 
    lhs[10] = 0.0; 
    lhs[11] = 0.0; 
    lhs[12] = 0.0; 
    lhs[13] = 0.0; 
    lhs[14] = 0.0; 
    lhs[15] = 0.0; 
    lhs[16] = 0.0; 
    lhs[17] = 0.0; 
    lhs[18] = 0.0; 
    lhs[19] = 0.0; 
    gkyl_mat_set(rhs,1,0,f[0]); 
    gkyl_mat_set(rhs,2,0,0.0); 
    gkyl_mat_set(rhs,3,0,0.0); 
    gkyl_mat_set(rhs,4,0,0.0); 
    gkyl_mat_set(rhs,5,0,0.0); 
    gkyl_mat_set(rhs,6,0,0.0); 
    gkyl_mat_set(rhs,7,0,0.0); 
    gkyl_mat_set(rhs,8,0,0.0); 
    gkyl_mat_set(rhs,9,0,0.0); 
    gkyl_mat_set(rhs,10,0,0.0); 
    gkyl_mat_set(rhs,11,0,0.0); 
    gkyl_mat_set(rhs,12,0,0.0); 
    gkyl_mat_set(rhs,13,0,0.0); 
    gkyl_mat_set(rhs,14,0,0.0); 
    gkyl_mat_set(rhs,15,0,0.0); 
    gkyl_mat_set(rhs,16,0,0.0); 
    gkyl_mat_set(rhs,17,0,0.0); 
    gkyl_mat_set(rhs,18,0,0.0); 
    gkyl_mat_set(rhs,19,0,0.0); 
    gkyl_mat_set(rhs,20,0,0.0); 
  } else { 
    lhs[0] = g[0]; 
    lhs[1] = g[1]; 
    lhs[2] = g[2]; 
    lhs[3] = g[3]; 
    lhs[4] = g[4]; 
    lhs[5] = g[5]; 
    lhs[6] = g[6]; 
    lhs[7] = g[7]; 
    lhs[8] = g[8]; 
    lhs[9] = g[9]; 
    lhs[10] = g[10]; 
    lhs[11] = g[11]; 
    lhs[12] = g[12]; 
    lhs[13] = g[13]; 
    lhs[14] = g[14]; 
    lhs[15] = g[15]; 
    lhs[16] = g[16]; 
    lhs[17] = g[17]; 
    lhs[18] = g[18]; 
    lhs[19] = g[19]; 
    gkyl_mat_set(rhs,1,0,f[0]); 
    gkyl_mat_set(rhs,2,0,f[1]); 
    gkyl_mat_set(rhs,3,0,f[2]); 
    gkyl_mat_set(rhs,4,0,f[3]); 
    gkyl_mat_set(rhs,5,0,f[4]); 
    gkyl_mat_set(rhs,6,0,f[5]); 
    gkyl_mat_set(rhs,7,0,f[6]); 
    gkyl_mat_set(rhs,8,0,f[7]); 
    gkyl_mat_set(rhs,9,0,f[8]); 
    gkyl_mat_set(rhs,10,0,f[9]); 
    gkyl_mat_set(rhs,11,0,f[10]); 
    gkyl_mat_set(rhs,12,0,f[11]); 
    gkyl_mat_set(rhs,13,0,f[12]); 
    gkyl_mat_set(rhs,14,0,f[13]); 
    gkyl_mat_set(rhs,15,0,f[14]); 
    gkyl_mat_set(rhs,16,0,f[15]); 
    gkyl_mat_set(rhs,17,0,f[16]); 
    gkyl_mat_set(rhs,18,0,f[17]); 
    gkyl_mat_set(rhs,19,0,f[18]); 
    gkyl_mat_set(rhs,20,0,f[19]); 
  } 
 
  // Fill LHS matrix. 
  gkyl_mat_set(A,0,0,0.3535533905932737*lhs[0]); 
  gkyl_mat_set(A,0,1,0.3535533905932737*lhs[1]); 
  gkyl_mat_set(A,0,2,0.3535533905932737*lhs[2]); 
  gkyl_mat_set(A,0,3,0.3535533905932737*lhs[3]); 
  gkyl_mat_set(A,0,4,0.3535533905932737*lhs[4]); 
  gkyl_mat_set(A,0,5,0.3535533905932737*lhs[5]); 
  gkyl_mat_set(A,0,6,0.3535533905932737*lhs[6]); 
  gkyl_mat_set(A,0,7,0.3535533905932737*lhs[7]); 
  gkyl_mat_set(A,0,8,0.3535533905932737*lhs[8]); 
  gkyl_mat_set(A,0,9,0.3535533905932737*lhs[9]); 
  gkyl_mat_set(A,0,10,0.3535533905932737*lhs[10]); 
  gkyl_mat_set(A,0,11,0.3535533905932737*lhs[11]); 
  gkyl_mat_set(A,0,12,0.3535533905932737*lhs[12]); 
  gkyl_mat_set(A,0,13,0.3535533905932737*lhs[13]); 
  gkyl_mat_set(A,0,14,0.3535533905932737*lhs[14]); 
  gkyl_mat_set(A,0,15,0.3535533905932737*lhs[15]); 
  gkyl_mat_set(A,0,16,0.3535533905932737*lhs[16]); 
  gkyl_mat_set(A,0,17,0.3535533905932737*lhs[17]); 
  gkyl_mat_set(A,0,18,0.3535533905932737*lhs[18]); 
  gkyl_mat_set(A,0,19,0.3535533905932737*lhs[19]); 
  gkyl_mat_set(A,1,0,0.3535533905932737*lhs[1]); 
  gkyl_mat_set(A,1,1,0.3162277660168379*lhs[7]+0.3535533905932737*lhs[0]); 
  gkyl_mat_set(A,1,2,0.3535533905932737*lhs[4]); 
  gkyl_mat_set(A,1,3,0.3535533905932737*lhs[5]); 
  gkyl_mat_set(A,1,4,0.3162277660168379*lhs[11]+0.3535533905932737*lhs[2]); 
  gkyl_mat_set(A,1,5,0.3162277660168379*lhs[13]+0.3535533905932737*lhs[3]); 
  gkyl_mat_set(A,1,6,0.3535533905932737*lhs[10]); 
  gkyl_mat_set(A,1,7,0.3162277660168379*lhs[1]); 
  gkyl_mat_set(A,1,8,0.3535533905932737*lhs[12]); 
  gkyl_mat_set(A,1,9,0.3535533905932737*lhs[15]); 
  gkyl_mat_set(A,1,10,0.3162277660168379*lhs[17]+0.3535533905932737*lhs[6]); 
  gkyl_mat_set(A,1,11,0.3162277660168379*lhs[4]); 
  gkyl_mat_set(A,1,12,0.3535533905932737*lhs[8]); 
  gkyl_mat_set(A,1,13,0.3162277660168379*lhs[5]); 
  gkyl_mat_set(A,1,14,0.3535533905932737*lhs[18]); 
  gkyl_mat_set(A,1,15,0.3535533905932737*lhs[9]); 
  gkyl_mat_set(A,1,16,0.3535533905932737*lhs[19]); 
  gkyl_mat_set(A,1,17,0.3162277660168379*lhs[10]); 
  gkyl_mat_set(A,1,18,0.3535533905932737*lhs[14]); 
  gkyl_mat_set(A,1,19,0.3535533905932737*lhs[16]); 
  gkyl_mat_set(A,2,0,0.3535533905932737*lhs[2]); 
  gkyl_mat_set(A,2,1,0.3535533905932737*lhs[4]); 
  gkyl_mat_set(A,2,2,0.3162277660168379*lhs[8]+0.3535533905932737*lhs[0]); 
  gkyl_mat_set(A,2,3,0.3535533905932737*lhs[6]); 
  gkyl_mat_set(A,2,4,0.3162277660168379*lhs[12]+0.3535533905932737*lhs[1]); 
  gkyl_mat_set(A,2,5,0.3535533905932737*lhs[10]); 
  gkyl_mat_set(A,2,6,0.3162277660168379*lhs[14]+0.3535533905932737*lhs[3]); 
  gkyl_mat_set(A,2,7,0.3535533905932737*lhs[11]); 
  gkyl_mat_set(A,2,8,0.3162277660168379*lhs[2]); 
  gkyl_mat_set(A,2,9,0.3535533905932737*lhs[16]); 
  gkyl_mat_set(A,2,10,0.3162277660168379*lhs[18]+0.3535533905932737*lhs[5]); 
  gkyl_mat_set(A,2,11,0.3535533905932737*lhs[7]); 
  gkyl_mat_set(A,2,12,0.3162277660168379*lhs[4]); 
  gkyl_mat_set(A,2,13,0.3535533905932737*lhs[17]); 
  gkyl_mat_set(A,2,14,0.3162277660168379*lhs[6]); 
  gkyl_mat_set(A,2,15,0.3535533905932737*lhs[19]); 
  gkyl_mat_set(A,2,16,0.3535533905932737*lhs[9]); 
  gkyl_mat_set(A,2,17,0.3535533905932737*lhs[13]); 
  gkyl_mat_set(A,2,18,0.3162277660168379*lhs[10]); 
  gkyl_mat_set(A,2,19,0.3535533905932737*lhs[15]); 
  gkyl_mat_set(A,3,0,0.3535533905932737*lhs[3]); 
  gkyl_mat_set(A,3,1,0.3535533905932737*lhs[5]); 
  gkyl_mat_set(A,3,2,0.3535533905932737*lhs[6]); 
  gkyl_mat_set(A,3,3,0.3162277660168379*lhs[9]+0.3535533905932737*lhs[0]); 
  gkyl_mat_set(A,3,4,0.3535533905932737*lhs[10]); 
  gkyl_mat_set(A,3,5,0.3162277660168379*lhs[15]+0.3535533905932737*lhs[1]); 
  gkyl_mat_set(A,3,6,0.3162277660168379*lhs[16]+0.3535533905932737*lhs[2]); 
  gkyl_mat_set(A,3,7,0.3535533905932737*lhs[13]); 
  gkyl_mat_set(A,3,8,0.3535533905932737*lhs[14]); 
  gkyl_mat_set(A,3,9,0.3162277660168379*lhs[3]); 
  gkyl_mat_set(A,3,10,0.3162277660168379*lhs[19]+0.3535533905932737*lhs[4]); 
  gkyl_mat_set(A,3,11,0.3535533905932737*lhs[17]); 
  gkyl_mat_set(A,3,12,0.3535533905932737*lhs[18]); 
  gkyl_mat_set(A,3,13,0.3535533905932737*lhs[7]); 
  gkyl_mat_set(A,3,14,0.3535533905932737*lhs[8]); 
  gkyl_mat_set(A,3,15,0.3162277660168379*lhs[5]); 
  gkyl_mat_set(A,3,16,0.3162277660168379*lhs[6]); 
  gkyl_mat_set(A,3,17,0.3535533905932737*lhs[11]); 
  gkyl_mat_set(A,3,18,0.3535533905932737*lhs[12]); 
  gkyl_mat_set(A,3,19,0.3162277660168379*lhs[10]); 
  gkyl_mat_set(A,4,0,0.3535533905932737*lhs[4]); 
  gkyl_mat_set(A,4,1,0.3162277660168379*lhs[11]+0.3535533905932737*lhs[2]); 
  gkyl_mat_set(A,4,2,0.3162277660168379*lhs[12]+0.3535533905932737*lhs[1]); 
  gkyl_mat_set(A,4,3,0.3535533905932737*lhs[10]); 
  gkyl_mat_set(A,4,4,0.3162277660168379*lhs[8]+0.3162277660168379*lhs[7]+0.3535533905932737*lhs[0]); 
  gkyl_mat_set(A,4,5,0.3162277660168379*lhs[17]+0.3535533905932737*lhs[6]); 
  gkyl_mat_set(A,4,6,0.3162277660168379*lhs[18]+0.3535533905932737*lhs[5]); 
  gkyl_mat_set(A,4,7,0.3162277660168379*lhs[4]); 
  gkyl_mat_set(A,4,8,0.3162277660168379*lhs[4]); 
  gkyl_mat_set(A,4,9,0.3535533905932737*lhs[19]); 
  gkyl_mat_set(A,4,10,0.3162277660168379*lhs[14]+0.3162277660168379*lhs[13]+0.3535533905932737*lhs[3]); 
  gkyl_mat_set(A,4,11,0.282842712474619*lhs[12]+0.3162277660168379*lhs[1]); 
  gkyl_mat_set(A,4,12,0.282842712474619*lhs[11]+0.3162277660168379*lhs[2]); 
  gkyl_mat_set(A,4,13,0.3162277660168379*lhs[10]); 
  gkyl_mat_set(A,4,14,0.3162277660168379*lhs[10]); 
  gkyl_mat_set(A,4,15,0.3535533905932737*lhs[16]); 
  gkyl_mat_set(A,4,16,0.3535533905932737*lhs[15]); 
  gkyl_mat_set(A,4,17,0.282842712474619*lhs[18]+0.3162277660168379*lhs[5]); 
  gkyl_mat_set(A,4,18,0.282842712474619*lhs[17]+0.3162277660168379*lhs[6]); 
  gkyl_mat_set(A,4,19,0.3535533905932737*lhs[9]); 
  gkyl_mat_set(A,5,0,0.3535533905932737*lhs[5]); 
  gkyl_mat_set(A,5,1,0.3162277660168379*lhs[13]+0.3535533905932737*lhs[3]); 
  gkyl_mat_set(A,5,2,0.3535533905932737*lhs[10]); 
  gkyl_mat_set(A,5,3,0.3162277660168379*lhs[15]+0.3535533905932737*lhs[1]); 
  gkyl_mat_set(A,5,4,0.3162277660168379*lhs[17]+0.3535533905932737*lhs[6]); 
  gkyl_mat_set(A,5,5,0.3162277660168379*lhs[9]+0.3162277660168379*lhs[7]+0.3535533905932737*lhs[0]); 
  gkyl_mat_set(A,5,6,0.3162277660168379*lhs[19]+0.3535533905932737*lhs[4]); 
  gkyl_mat_set(A,5,7,0.3162277660168379*lhs[5]); 
  gkyl_mat_set(A,5,8,0.3535533905932737*lhs[18]); 
  gkyl_mat_set(A,5,9,0.3162277660168379*lhs[5]); 
  gkyl_mat_set(A,5,10,0.3162277660168379*lhs[16]+0.3162277660168379*lhs[11]+0.3535533905932737*lhs[2]); 
  gkyl_mat_set(A,5,11,0.3162277660168379*lhs[10]); 
  gkyl_mat_set(A,5,12,0.3535533905932737*lhs[14]); 
  gkyl_mat_set(A,5,13,0.282842712474619*lhs[15]+0.3162277660168379*lhs[1]); 
  gkyl_mat_set(A,5,14,0.3535533905932737*lhs[12]); 
  gkyl_mat_set(A,5,15,0.282842712474619*lhs[13]+0.3162277660168379*lhs[3]); 
  gkyl_mat_set(A,5,16,0.3162277660168379*lhs[10]); 
  gkyl_mat_set(A,5,17,0.282842712474619*lhs[19]+0.3162277660168379*lhs[4]); 
  gkyl_mat_set(A,5,18,0.3535533905932737*lhs[8]); 
  gkyl_mat_set(A,5,19,0.282842712474619*lhs[17]+0.3162277660168379*lhs[6]); 
  gkyl_mat_set(A,6,0,0.3535533905932737*lhs[6]); 
  gkyl_mat_set(A,6,1,0.3535533905932737*lhs[10]); 
  gkyl_mat_set(A,6,2,0.3162277660168379*lhs[14]+0.3535533905932737*lhs[3]); 
  gkyl_mat_set(A,6,3,0.3162277660168379*lhs[16]+0.3535533905932737*lhs[2]); 
  gkyl_mat_set(A,6,4,0.3162277660168379*lhs[18]+0.3535533905932737*lhs[5]); 
  gkyl_mat_set(A,6,5,0.3162277660168379*lhs[19]+0.3535533905932737*lhs[4]); 
  gkyl_mat_set(A,6,6,0.3162277660168379*lhs[9]+0.3162277660168379*lhs[8]+0.3535533905932737*lhs[0]); 
  gkyl_mat_set(A,6,7,0.3535533905932737*lhs[17]); 
  gkyl_mat_set(A,6,8,0.3162277660168379*lhs[6]); 
  gkyl_mat_set(A,6,9,0.3162277660168379*lhs[6]); 
  gkyl_mat_set(A,6,10,0.3162277660168379*lhs[15]+0.3162277660168379*lhs[12]+0.3535533905932737*lhs[1]); 
  gkyl_mat_set(A,6,11,0.3535533905932737*lhs[13]); 
  gkyl_mat_set(A,6,12,0.3162277660168379*lhs[10]); 
  gkyl_mat_set(A,6,13,0.3535533905932737*lhs[11]); 
  gkyl_mat_set(A,6,14,0.282842712474619*lhs[16]+0.3162277660168379*lhs[2]); 
  gkyl_mat_set(A,6,15,0.3162277660168379*lhs[10]); 
  gkyl_mat_set(A,6,16,0.282842712474619*lhs[14]+0.3162277660168379*lhs[3]); 
  gkyl_mat_set(A,6,17,0.3535533905932737*lhs[7]); 
  gkyl_mat_set(A,6,18,0.282842712474619*lhs[19]+0.3162277660168379*lhs[4]); 
  gkyl_mat_set(A,6,19,0.282842712474619*lhs[18]+0.3162277660168379*lhs[5]); 
  gkyl_mat_set(A,7,0,0.3535533905932737*lhs[7]); 
  gkyl_mat_set(A,7,1,0.3162277660168379*lhs[1]); 
  gkyl_mat_set(A,7,2,0.3535533905932737*lhs[11]); 
  gkyl_mat_set(A,7,3,0.3535533905932737*lhs[13]); 
  gkyl_mat_set(A,7,4,0.3162277660168379*lhs[4]); 
  gkyl_mat_set(A,7,5,0.3162277660168379*lhs[5]); 
  gkyl_mat_set(A,7,6,0.3535533905932737*lhs[17]); 
  gkyl_mat_set(A,7,7,0.2258769757263127*lhs[7]+0.3535533905932737*lhs[0]); 
  gkyl_mat_set(A,7,10,0.3162277660168379*lhs[10]); 
  gkyl_mat_set(A,7,11,0.2258769757263127*lhs[11]+0.3535533905932737*lhs[2]); 
  gkyl_mat_set(A,7,12,0.3162277660168379*lhs[12]); 
  gkyl_mat_set(A,7,13,0.2258769757263127*lhs[13]+0.3535533905932737*lhs[3]); 
  gkyl_mat_set(A,7,15,0.3162277660168379*lhs[15]); 
  gkyl_mat_set(A,7,17,0.2258769757263127*lhs[17]+0.3535533905932737*lhs[6]); 
  gkyl_mat_set(A,7,18,0.3162277660168379*lhs[18]); 
  gkyl_mat_set(A,7,19,0.3162277660168379*lhs[19]); 
  gkyl_mat_set(A,8,0,0.3535533905932737*lhs[8]); 
  gkyl_mat_set(A,8,1,0.3535533905932737*lhs[12]); 
  gkyl_mat_set(A,8,2,0.3162277660168379*lhs[2]); 
  gkyl_mat_set(A,8,3,0.3535533905932737*lhs[14]); 
  gkyl_mat_set(A,8,4,0.3162277660168379*lhs[4]); 
  gkyl_mat_set(A,8,5,0.3535533905932737*lhs[18]); 
  gkyl_mat_set(A,8,6,0.3162277660168379*lhs[6]); 
  gkyl_mat_set(A,8,8,0.2258769757263127*lhs[8]+0.3535533905932737*lhs[0]); 
  gkyl_mat_set(A,8,10,0.3162277660168379*lhs[10]); 
  gkyl_mat_set(A,8,11,0.3162277660168379*lhs[11]); 
  gkyl_mat_set(A,8,12,0.2258769757263127*lhs[12]+0.3535533905932737*lhs[1]); 
  gkyl_mat_set(A,8,14,0.2258769757263127*lhs[14]+0.3535533905932737*lhs[3]); 
  gkyl_mat_set(A,8,16,0.3162277660168379*lhs[16]); 
  gkyl_mat_set(A,8,17,0.3162277660168379*lhs[17]); 
  gkyl_mat_set(A,8,18,0.2258769757263127*lhs[18]+0.3535533905932737*lhs[5]); 
  gkyl_mat_set(A,8,19,0.3162277660168379*lhs[19]); 
  gkyl_mat_set(A,9,0,0.3535533905932737*lhs[9]); 
  gkyl_mat_set(A,9,1,0.3535533905932737*lhs[15]); 
  gkyl_mat_set(A,9,2,0.3535533905932737*lhs[16]); 
  gkyl_mat_set(A,9,3,0.3162277660168379*lhs[3]); 
  gkyl_mat_set(A,9,4,0.3535533905932737*lhs[19]); 
  gkyl_mat_set(A,9,5,0.3162277660168379*lhs[5]); 
  gkyl_mat_set(A,9,6,0.3162277660168379*lhs[6]); 
  gkyl_mat_set(A,9,9,0.2258769757263127*lhs[9]+0.3535533905932737*lhs[0]); 
  gkyl_mat_set(A,9,10,0.3162277660168379*lhs[10]); 
  gkyl_mat_set(A,9,13,0.3162277660168379*lhs[13]); 
  gkyl_mat_set(A,9,14,0.3162277660168379*lhs[14]); 
  gkyl_mat_set(A,9,15,0.2258769757263127*lhs[15]+0.3535533905932737*lhs[1]); 
  gkyl_mat_set(A,9,16,0.2258769757263127*lhs[16]+0.3535533905932737*lhs[2]); 
  gkyl_mat_set(A,9,17,0.3162277660168379*lhs[17]); 
  gkyl_mat_set(A,9,18,0.3162277660168379*lhs[18]); 
  gkyl_mat_set(A,9,19,0.2258769757263127*lhs[19]+0.3535533905932737*lhs[4]); 
  gkyl_mat_set(A,10,0,0.3535533905932737*lhs[10]); 
  gkyl_mat_set(A,10,1,0.3162277660168379*lhs[17]+0.3535533905932737*lhs[6]); 
  gkyl_mat_set(A,10,2,0.3162277660168379*lhs[18]+0.3535533905932737*lhs[5]); 
  gkyl_mat_set(A,10,3,0.3162277660168379*lhs[19]+0.3535533905932737*lhs[4]); 
  gkyl_mat_set(A,10,4,0.3162277660168379*lhs[14]+0.3162277660168379*lhs[13]+0.3535533905932737*lhs[3]); 
  gkyl_mat_set(A,10,5,0.3162277660168379*lhs[16]+0.3162277660168379*lhs[11]+0.3535533905932737*lhs[2]); 
  gkyl_mat_set(A,10,6,0.3162277660168379*lhs[15]+0.3162277660168379*lhs[12]+0.3535533905932737*lhs[1]); 
  gkyl_mat_set(A,10,7,0.3162277660168379*lhs[10]); 
  gkyl_mat_set(A,10,8,0.3162277660168379*lhs[10]); 
  gkyl_mat_set(A,10,9,0.3162277660168379*lhs[10]); 
  gkyl_mat_set(A,10,10,0.3162277660168379*lhs[9]+0.3162277660168379*lhs[8]+0.3162277660168379*lhs[7]+0.3535533905932737*lhs[0]); 
  gkyl_mat_set(A,10,11,0.282842712474619*lhs[18]+0.3162277660168379*lhs[5]); 
  gkyl_mat_set(A,10,12,0.282842712474619*lhs[17]+0.3162277660168379*lhs[6]); 
  gkyl_mat_set(A,10,13,0.282842712474619*lhs[19]+0.3162277660168379*lhs[4]); 
  gkyl_mat_set(A,10,14,0.282842712474619*lhs[19]+0.3162277660168379*lhs[4]); 
  gkyl_mat_set(A,10,15,0.282842712474619*lhs[17]+0.3162277660168379*lhs[6]); 
  gkyl_mat_set(A,10,16,0.282842712474619*lhs[18]+0.3162277660168379*lhs[5]); 
  gkyl_mat_set(A,10,17,0.282842712474619*lhs[15]+0.282842712474619*lhs[12]+0.3162277660168379*lhs[1]); 
  gkyl_mat_set(A,10,18,0.282842712474619*lhs[16]+0.282842712474619*lhs[11]+0.3162277660168379*lhs[2]); 
  gkyl_mat_set(A,10,19,0.282842712474619*lhs[14]+0.282842712474619*lhs[13]+0.3162277660168379*lhs[3]); 
  gkyl_mat_set(A,11,0,0.3535533905932737*lhs[11]); 
  gkyl_mat_set(A,11,1,0.3162277660168379*lhs[4]); 
  gkyl_mat_set(A,11,2,0.3535533905932737*lhs[7]); 
  gkyl_mat_set(A,11,3,0.3535533905932737*lhs[17]); 
  gkyl_mat_set(A,11,4,0.282842712474619*lhs[12]+0.3162277660168379*lhs[1]); 
  gkyl_mat_set(A,11,5,0.3162277660168379*lhs[10]); 
  gkyl_mat_set(A,11,6,0.3535533905932737*lhs[13]); 
  gkyl_mat_set(A,11,7,0.2258769757263127*lhs[11]+0.3535533905932737*lhs[2]); 
  gkyl_mat_set(A,11,8,0.3162277660168379*lhs[11]); 
  gkyl_mat_set(A,11,10,0.282842712474619*lhs[18]+0.3162277660168379*lhs[5]); 
  gkyl_mat_set(A,11,11,0.3162277660168379*lhs[8]+0.2258769757263127*lhs[7]+0.3535533905932737*lhs[0]); 
  gkyl_mat_set(A,11,12,0.282842712474619*lhs[4]); 
  gkyl_mat_set(A,11,13,0.2258769757263127*lhs[17]+0.3535533905932737*lhs[6]); 
  gkyl_mat_set(A,11,14,0.3162277660168379*lhs[17]); 
  gkyl_mat_set(A,11,15,0.3162277660168379*lhs[19]); 
  gkyl_mat_set(A,11,17,0.3162277660168379*lhs[14]+0.2258769757263127*lhs[13]+0.3535533905932737*lhs[3]); 
  gkyl_mat_set(A,11,18,0.282842712474619*lhs[10]); 
  gkyl_mat_set(A,11,19,0.3162277660168379*lhs[15]); 
  gkyl_mat_set(A,12,0,0.3535533905932737*lhs[12]); 
  gkyl_mat_set(A,12,1,0.3535533905932737*lhs[8]); 
  gkyl_mat_set(A,12,2,0.3162277660168379*lhs[4]); 
  gkyl_mat_set(A,12,3,0.3535533905932737*lhs[18]); 
  gkyl_mat_set(A,12,4,0.282842712474619*lhs[11]+0.3162277660168379*lhs[2]); 
  gkyl_mat_set(A,12,5,0.3535533905932737*lhs[14]); 
  gkyl_mat_set(A,12,6,0.3162277660168379*lhs[10]); 
  gkyl_mat_set(A,12,7,0.3162277660168379*lhs[12]); 
  gkyl_mat_set(A,12,8,0.2258769757263127*lhs[12]+0.3535533905932737*lhs[1]); 
  gkyl_mat_set(A,12,10,0.282842712474619*lhs[17]+0.3162277660168379*lhs[6]); 
  gkyl_mat_set(A,12,11,0.282842712474619*lhs[4]); 
  gkyl_mat_set(A,12,12,0.2258769757263127*lhs[8]+0.3162277660168379*lhs[7]+0.3535533905932737*lhs[0]); 
  gkyl_mat_set(A,12,13,0.3162277660168379*lhs[18]); 
  gkyl_mat_set(A,12,14,0.2258769757263127*lhs[18]+0.3535533905932737*lhs[5]); 
  gkyl_mat_set(A,12,16,0.3162277660168379*lhs[19]); 
  gkyl_mat_set(A,12,17,0.282842712474619*lhs[10]); 
  gkyl_mat_set(A,12,18,0.2258769757263127*lhs[14]+0.3162277660168379*lhs[13]+0.3535533905932737*lhs[3]); 
  gkyl_mat_set(A,12,19,0.3162277660168379*lhs[16]); 
  gkyl_mat_set(A,13,0,0.3535533905932737*lhs[13]); 
  gkyl_mat_set(A,13,1,0.3162277660168379*lhs[5]); 
  gkyl_mat_set(A,13,2,0.3535533905932737*lhs[17]); 
  gkyl_mat_set(A,13,3,0.3535533905932737*lhs[7]); 
  gkyl_mat_set(A,13,4,0.3162277660168379*lhs[10]); 
  gkyl_mat_set(A,13,5,0.282842712474619*lhs[15]+0.3162277660168379*lhs[1]); 
  gkyl_mat_set(A,13,6,0.3535533905932737*lhs[11]); 
  gkyl_mat_set(A,13,7,0.2258769757263127*lhs[13]+0.3535533905932737*lhs[3]); 
  gkyl_mat_set(A,13,9,0.3162277660168379*lhs[13]); 
  gkyl_mat_set(A,13,10,0.282842712474619*lhs[19]+0.3162277660168379*lhs[4]); 
  gkyl_mat_set(A,13,11,0.2258769757263127*lhs[17]+0.3535533905932737*lhs[6]); 
  gkyl_mat_set(A,13,12,0.3162277660168379*lhs[18]); 
  gkyl_mat_set(A,13,13,0.3162277660168379*lhs[9]+0.2258769757263127*lhs[7]+0.3535533905932737*lhs[0]); 
  gkyl_mat_set(A,13,15,0.282842712474619*lhs[5]); 
  gkyl_mat_set(A,13,16,0.3162277660168379*lhs[17]); 
  gkyl_mat_set(A,13,17,0.3162277660168379*lhs[16]+0.2258769757263127*lhs[11]+0.3535533905932737*lhs[2]); 
  gkyl_mat_set(A,13,18,0.3162277660168379*lhs[12]); 
  gkyl_mat_set(A,13,19,0.282842712474619*lhs[10]); 
  gkyl_mat_set(A,14,0,0.3535533905932737*lhs[14]); 
  gkyl_mat_set(A,14,1,0.3535533905932737*lhs[18]); 
  gkyl_mat_set(A,14,2,0.3162277660168379*lhs[6]); 
  gkyl_mat_set(A,14,3,0.3535533905932737*lhs[8]); 
  gkyl_mat_set(A,14,4,0.3162277660168379*lhs[10]); 
  gkyl_mat_set(A,14,5,0.3535533905932737*lhs[12]); 
  gkyl_mat_set(A,14,6,0.282842712474619*lhs[16]+0.3162277660168379*lhs[2]); 
  gkyl_mat_set(A,14,8,0.2258769757263127*lhs[14]+0.3535533905932737*lhs[3]); 
  gkyl_mat_set(A,14,9,0.3162277660168379*lhs[14]); 
  gkyl_mat_set(A,14,10,0.282842712474619*lhs[19]+0.3162277660168379*lhs[4]); 
  gkyl_mat_set(A,14,11,0.3162277660168379*lhs[17]); 
  gkyl_mat_set(A,14,12,0.2258769757263127*lhs[18]+0.3535533905932737*lhs[5]); 
  gkyl_mat_set(A,14,14,0.3162277660168379*lhs[9]+0.2258769757263127*lhs[8]+0.3535533905932737*lhs[0]); 
  gkyl_mat_set(A,14,15,0.3162277660168379*lhs[18]); 
  gkyl_mat_set(A,14,16,0.282842712474619*lhs[6]); 
  gkyl_mat_set(A,14,17,0.3162277660168379*lhs[11]); 
  gkyl_mat_set(A,14,18,0.3162277660168379*lhs[15]+0.2258769757263127*lhs[12]+0.3535533905932737*lhs[1]); 
  gkyl_mat_set(A,14,19,0.282842712474619*lhs[10]); 
  gkyl_mat_set(A,15,0,0.3535533905932737*lhs[15]); 
  gkyl_mat_set(A,15,1,0.3535533905932737*lhs[9]); 
  gkyl_mat_set(A,15,2,0.3535533905932737*lhs[19]); 
  gkyl_mat_set(A,15,3,0.3162277660168379*lhs[5]); 
  gkyl_mat_set(A,15,4,0.3535533905932737*lhs[16]); 
  gkyl_mat_set(A,15,5,0.282842712474619*lhs[13]+0.3162277660168379*lhs[3]); 
  gkyl_mat_set(A,15,6,0.3162277660168379*lhs[10]); 
  gkyl_mat_set(A,15,7,0.3162277660168379*lhs[15]); 
  gkyl_mat_set(A,15,9,0.2258769757263127*lhs[15]+0.3535533905932737*lhs[1]); 
  gkyl_mat_set(A,15,10,0.282842712474619*lhs[17]+0.3162277660168379*lhs[6]); 
  gkyl_mat_set(A,15,11,0.3162277660168379*lhs[19]); 
  gkyl_mat_set(A,15,13,0.282842712474619*lhs[5]); 
  gkyl_mat_set(A,15,14,0.3162277660168379*lhs[18]); 
  gkyl_mat_set(A,15,15,0.2258769757263127*lhs[9]+0.3162277660168379*lhs[7]+0.3535533905932737*lhs[0]); 
  gkyl_mat_set(A,15,16,0.2258769757263127*lhs[19]+0.3535533905932737*lhs[4]); 
  gkyl_mat_set(A,15,17,0.282842712474619*lhs[10]); 
  gkyl_mat_set(A,15,18,0.3162277660168379*lhs[14]); 
  gkyl_mat_set(A,15,19,0.2258769757263127*lhs[16]+0.3162277660168379*lhs[11]+0.3535533905932737*lhs[2]); 
  gkyl_mat_set(A,16,0,0.3535533905932737*lhs[16]); 
  gkyl_mat_set(A,16,1,0.3535533905932737*lhs[19]); 
  gkyl_mat_set(A,16,2,0.3535533905932737*lhs[9]); 
  gkyl_mat_set(A,16,3,0.3162277660168379*lhs[6]); 
  gkyl_mat_set(A,16,4,0.3535533905932737*lhs[15]); 
  gkyl_mat_set(A,16,5,0.3162277660168379*lhs[10]); 
  gkyl_mat_set(A,16,6,0.282842712474619*lhs[14]+0.3162277660168379*lhs[3]); 
  gkyl_mat_set(A,16,8,0.3162277660168379*lhs[16]); 
  gkyl_mat_set(A,16,9,0.2258769757263127*lhs[16]+0.3535533905932737*lhs[2]); 
  gkyl_mat_set(A,16,10,0.282842712474619*lhs[18]+0.3162277660168379*lhs[5]); 
  gkyl_mat_set(A,16,12,0.3162277660168379*lhs[19]); 
  gkyl_mat_set(A,16,13,0.3162277660168379*lhs[17]); 
  gkyl_mat_set(A,16,14,0.282842712474619*lhs[6]); 
  gkyl_mat_set(A,16,15,0.2258769757263127*lhs[19]+0.3535533905932737*lhs[4]); 
  gkyl_mat_set(A,16,16,0.2258769757263127*lhs[9]+0.3162277660168379*lhs[8]+0.3535533905932737*lhs[0]); 
  gkyl_mat_set(A,16,17,0.3162277660168379*lhs[13]); 
  gkyl_mat_set(A,16,18,0.282842712474619*lhs[10]); 
  gkyl_mat_set(A,16,19,0.2258769757263127*lhs[15]+0.3162277660168379*lhs[12]+0.3535533905932737*lhs[1]); 
  gkyl_mat_set(A,17,0,0.3535533905932737*lhs[17]); 
  gkyl_mat_set(A,17,1,0.3162277660168379*lhs[10]); 
  gkyl_mat_set(A,17,2,0.3535533905932737*lhs[13]); 
  gkyl_mat_set(A,17,3,0.3535533905932737*lhs[11]); 
  gkyl_mat_set(A,17,4,0.282842712474619*lhs[18]+0.3162277660168379*lhs[5]); 
  gkyl_mat_set(A,17,5,0.282842712474619*lhs[19]+0.3162277660168379*lhs[4]); 
  gkyl_mat_set(A,17,6,0.3535533905932737*lhs[7]); 
  gkyl_mat_set(A,17,7,0.2258769757263127*lhs[17]+0.3535533905932737*lhs[6]); 
  gkyl_mat_set(A,17,8,0.3162277660168379*lhs[17]); 
  gkyl_mat_set(A,17,9,0.3162277660168379*lhs[17]); 
  gkyl_mat_set(A,17,10,0.282842712474619*lhs[15]+0.282842712474619*lhs[12]+0.3162277660168379*lhs[1]); 
  gkyl_mat_set(A,17,11,0.3162277660168379*lhs[14]+0.2258769757263127*lhs[13]+0.3535533905932737*lhs[3]); 
  gkyl_mat_set(A,17,12,0.282842712474619*lhs[10]); 
  gkyl_mat_set(A,17,13,0.3162277660168379*lhs[16]+0.2258769757263127*lhs[11]+0.3535533905932737*lhs[2]); 
  gkyl_mat_set(A,17,14,0.3162277660168379*lhs[11]); 
  gkyl_mat_set(A,17,15,0.282842712474619*lhs[10]); 
  gkyl_mat_set(A,17,16,0.3162277660168379*lhs[13]); 
  gkyl_mat_set(A,17,17,0.3162277660168379*lhs[9]+0.3162277660168379*lhs[8]+0.2258769757263127*lhs[7]+0.3535533905932737*lhs[0]); 
  gkyl_mat_set(A,17,18,0.2529822128134704*lhs[19]+0.282842712474619*lhs[4]); 
  gkyl_mat_set(A,17,19,0.2529822128134704*lhs[18]+0.282842712474619*lhs[5]); 
  gkyl_mat_set(A,18,0,0.3535533905932737*lhs[18]); 
  gkyl_mat_set(A,18,1,0.3535533905932737*lhs[14]); 
  gkyl_mat_set(A,18,2,0.3162277660168379*lhs[10]); 
  gkyl_mat_set(A,18,3,0.3535533905932737*lhs[12]); 
  gkyl_mat_set(A,18,4,0.282842712474619*lhs[17]+0.3162277660168379*lhs[6]); 
  gkyl_mat_set(A,18,5,0.3535533905932737*lhs[8]); 
  gkyl_mat_set(A,18,6,0.282842712474619*lhs[19]+0.3162277660168379*lhs[4]); 
  gkyl_mat_set(A,18,7,0.3162277660168379*lhs[18]); 
  gkyl_mat_set(A,18,8,0.2258769757263127*lhs[18]+0.3535533905932737*lhs[5]); 
  gkyl_mat_set(A,18,9,0.3162277660168379*lhs[18]); 
  gkyl_mat_set(A,18,10,0.282842712474619*lhs[16]+0.282842712474619*lhs[11]+0.3162277660168379*lhs[2]); 
  gkyl_mat_set(A,18,11,0.282842712474619*lhs[10]); 
  gkyl_mat_set(A,18,12,0.2258769757263127*lhs[14]+0.3162277660168379*lhs[13]+0.3535533905932737*lhs[3]); 
  gkyl_mat_set(A,18,13,0.3162277660168379*lhs[12]); 
  gkyl_mat_set(A,18,14,0.3162277660168379*lhs[15]+0.2258769757263127*lhs[12]+0.3535533905932737*lhs[1]); 
  gkyl_mat_set(A,18,15,0.3162277660168379*lhs[14]); 
  gkyl_mat_set(A,18,16,0.282842712474619*lhs[10]); 
  gkyl_mat_set(A,18,17,0.2529822128134704*lhs[19]+0.282842712474619*lhs[4]); 
  gkyl_mat_set(A,18,18,0.3162277660168379*lhs[9]+0.2258769757263127*lhs[8]+0.3162277660168379*lhs[7]+0.3535533905932737*lhs[0]); 
  gkyl_mat_set(A,18,19,0.2529822128134704*lhs[17]+0.282842712474619*lhs[6]); 
  gkyl_mat_set(A,19,0,0.3535533905932737*lhs[19]); 
  gkyl_mat_set(A,19,1,0.3535533905932737*lhs[16]); 
  gkyl_mat_set(A,19,2,0.3535533905932737*lhs[15]); 
  gkyl_mat_set(A,19,3,0.3162277660168379*lhs[10]); 
  gkyl_mat_set(A,19,4,0.3535533905932737*lhs[9]); 
  gkyl_mat_set(A,19,5,0.282842712474619*lhs[17]+0.3162277660168379*lhs[6]); 
  gkyl_mat_set(A,19,6,0.282842712474619*lhs[18]+0.3162277660168379*lhs[5]); 
  gkyl_mat_set(A,19,7,0.3162277660168379*lhs[19]); 
  gkyl_mat_set(A,19,8,0.3162277660168379*lhs[19]); 
  gkyl_mat_set(A,19,9,0.2258769757263127*lhs[19]+0.3535533905932737*lhs[4]); 
  gkyl_mat_set(A,19,10,0.282842712474619*lhs[14]+0.282842712474619*lhs[13]+0.3162277660168379*lhs[3]); 
  gkyl_mat_set(A,19,11,0.3162277660168379*lhs[15]); 
  gkyl_mat_set(A,19,12,0.3162277660168379*lhs[16]); 
  gkyl_mat_set(A,19,13,0.282842712474619*lhs[10]); 
  gkyl_mat_set(A,19,14,0.282842712474619*lhs[10]); 
  gkyl_mat_set(A,19,15,0.2258769757263127*lhs[16]+0.3162277660168379*lhs[11]+0.3535533905932737*lhs[2]); 
  gkyl_mat_set(A,19,16,0.2258769757263127*lhs[15]+0.3162277660168379*lhs[12]+0.3535533905932737*lhs[1]); 
  gkyl_mat_set(A,19,17,0.2529822128134704*lhs[18]+0.282842712474619*lhs[5]); 
  gkyl_mat_set(A,19,18,0.2529822128134704*lhs[17]+0.282842712474619*lhs[6]); 
  gkyl_mat_set(A,19,19,0.2258769757263127*lhs[9]+0.3162277660168379*lhs[8]+0.3162277660168379*lhs[7]+0.3535533905932737*lhs[0]); 

  // Solve the system of equations. 
  long ipiv[20]; 
  gkyl_mat_linsolve_lu(A,rhs,ipiv); 
  for(size_t i=0; i<20; i++) 
  { 
    fdivg[i] = gkyl_mat_get(rhs,i,0); 
  } 
} 
