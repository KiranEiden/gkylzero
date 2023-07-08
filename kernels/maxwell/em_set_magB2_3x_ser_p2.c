#include <gkyl_mat.h> 
#include <gkyl_maxwell_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH int em_set_magB2_3x_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *em) 
{ 
  // A:     preallocated LHS matrix. 
  // rhs:   preallocated RHS vector. 
  // em:    Input electromagnetic fields. 
 
  const double *B_x = &em[60]; 
  const double *B_y = &em[80]; 
  const double *B_z = &em[100]; 
 
  // Calculate |B|^2 and set matrix to solve for 1/|B|^2. 
  double B_x_sq[20] = {0.0}; 
  binop_mul_3d_ser_p2(B_x, B_x, B_x_sq); 
 
  double B_y_sq[20] = {0.0}; 
  binop_mul_3d_ser_p2(B_y, B_y, B_y_sq); 
 
  double B_z_sq[20] = {0.0}; 
  binop_mul_3d_ser_p2(B_z, B_z, B_z_sq); 
 
  double magB2[20] = {0.0}; 

  magB2[0] = B_z_sq[0]+B_y_sq[0]+B_x_sq[0]; 
  magB2[1] = B_z_sq[1]+B_y_sq[1]+B_x_sq[1]; 
  magB2[2] = B_z_sq[2]+B_y_sq[2]+B_x_sq[2]; 
  magB2[3] = B_z_sq[3]+B_y_sq[3]+B_x_sq[3]; 
  magB2[4] = B_z_sq[4]+B_y_sq[4]+B_x_sq[4]; 
  magB2[5] = B_z_sq[5]+B_y_sq[5]+B_x_sq[5]; 
  magB2[6] = B_z_sq[6]+B_y_sq[6]+B_x_sq[6]; 
  magB2[7] = B_z_sq[7]+B_y_sq[7]+B_x_sq[7]; 
  magB2[8] = B_z_sq[8]+B_y_sq[8]+B_x_sq[8]; 
  magB2[9] = B_z_sq[9]+B_y_sq[9]+B_x_sq[9]; 
  magB2[10] = B_z_sq[10]+B_y_sq[10]+B_x_sq[10]; 
  magB2[11] = B_z_sq[11]+B_y_sq[11]+B_x_sq[11]; 
  magB2[12] = B_z_sq[12]+B_y_sq[12]+B_x_sq[12]; 
  magB2[13] = B_z_sq[13]+B_y_sq[13]+B_x_sq[13]; 
  magB2[14] = B_z_sq[14]+B_y_sq[14]+B_x_sq[14]; 
  magB2[15] = B_z_sq[15]+B_y_sq[15]+B_x_sq[15]; 
  magB2[16] = B_z_sq[16]+B_y_sq[16]+B_x_sq[16]; 
  magB2[17] = B_z_sq[17]+B_y_sq[17]+B_x_sq[17]; 
  magB2[18] = B_z_sq[18]+B_y_sq[18]+B_x_sq[18]; 
  magB2[19] = B_z_sq[19]+B_y_sq[19]+B_x_sq[19]; 

  int cell_avg = 0;
  // Check if |B|^2 < 0 at control points. 
  if (2.371708245126284*magB2[19]+2.371708245126284*magB2[18]+2.371708245126284*magB2[17]-1.369306393762915*magB2[16]-1.369306393762915*magB2[15]-1.369306393762915*magB2[14]-1.369306393762915*magB2[13]-1.369306393762915*magB2[12]-1.369306393762915*magB2[11]-1.837117307087383*magB2[10]+0.7905694150420947*magB2[9]+0.7905694150420947*magB2[8]+0.7905694150420947*magB2[7]+1.060660171779821*magB2[6]+1.060660171779821*magB2[5]+1.060660171779821*magB2[4]-0.6123724356957944*magB2[3]-0.6123724356957944*magB2[2]-0.6123724356957944*magB2[1]+0.3535533905932737*magB2[0] < 0.0) cell_avg = 1; 
  if ((-1.185854122563142*magB2[17])-1.369306393762915*magB2[16]-1.369306393762915*magB2[14]+0.6846531968814574*magB2[13]+0.6846531968814574*magB2[11]+0.7905694150420947*magB2[9]+0.7905694150420947*magB2[8]-0.3952847075210473*magB2[7]+1.060660171779821*magB2[6]-0.6123724356957944*magB2[3]-0.6123724356957944*magB2[2]+0.3535533905932737*magB2[0] < 0.0) cell_avg = 1; 
  if ((-2.371708245126284*magB2[19])-2.371708245126284*magB2[18]+2.371708245126284*magB2[17]-1.369306393762915*magB2[16]+1.369306393762915*magB2[15]-1.369306393762915*magB2[14]-1.369306393762915*magB2[13]+1.369306393762915*magB2[12]-1.369306393762915*magB2[11]+1.837117307087383*magB2[10]+0.7905694150420947*magB2[9]+0.7905694150420947*magB2[8]+0.7905694150420947*magB2[7]+1.060660171779821*magB2[6]-1.060660171779821*magB2[5]-1.060660171779821*magB2[4]-0.6123724356957944*magB2[3]-0.6123724356957944*magB2[2]+0.6123724356957944*magB2[1]+0.3535533905932737*magB2[0] < 0.0) cell_avg = 1; 
  if ((-1.185854122563142*magB2[18])-1.369306393762915*magB2[15]+0.6846531968814574*magB2[14]-1.369306393762915*magB2[13]+0.6846531968814574*magB2[12]+0.7905694150420947*magB2[9]-0.3952847075210473*magB2[8]+0.7905694150420947*magB2[7]+1.060660171779821*magB2[5]-0.6123724356957944*magB2[3]-0.6123724356957944*magB2[1]+0.3535533905932737*magB2[0] < 0.0) cell_avg = 1; 
  if (1.185854122563142*magB2[18]+1.369306393762915*magB2[15]+0.6846531968814574*magB2[14]-1.369306393762915*magB2[13]-0.6846531968814574*magB2[12]+0.7905694150420947*magB2[9]-0.3952847075210473*magB2[8]+0.7905694150420947*magB2[7]-1.060660171779821*magB2[5]-0.6123724356957944*magB2[3]+0.6123724356957944*magB2[1]+0.3535533905932737*magB2[0] < 0.0) cell_avg = 1; 
  if ((-2.371708245126284*magB2[19])+2.371708245126284*magB2[18]-2.371708245126284*magB2[17]+1.369306393762915*magB2[16]-1.369306393762915*magB2[15]-1.369306393762915*magB2[14]-1.369306393762915*magB2[13]-1.369306393762915*magB2[12]+1.369306393762915*magB2[11]+1.837117307087383*magB2[10]+0.7905694150420947*magB2[9]+0.7905694150420947*magB2[8]+0.7905694150420947*magB2[7]-1.060660171779821*magB2[6]+1.060660171779821*magB2[5]-1.060660171779821*magB2[4]-0.6123724356957944*magB2[3]+0.6123724356957944*magB2[2]-0.6123724356957944*magB2[1]+0.3535533905932737*magB2[0] < 0.0) cell_avg = 1; 
  if (1.185854122563142*magB2[17]+1.369306393762915*magB2[16]-1.369306393762915*magB2[14]+0.6846531968814574*magB2[13]-0.6846531968814574*magB2[11]+0.7905694150420947*magB2[9]+0.7905694150420947*magB2[8]-0.3952847075210473*magB2[7]-1.060660171779821*magB2[6]-0.6123724356957944*magB2[3]+0.6123724356957944*magB2[2]+0.3535533905932737*magB2[0] < 0.0) cell_avg = 1; 
  if (2.371708245126284*magB2[19]-2.371708245126284*magB2[18]-2.371708245126284*magB2[17]+1.369306393762915*magB2[16]+1.369306393762915*magB2[15]-1.369306393762915*magB2[14]-1.369306393762915*magB2[13]+1.369306393762915*magB2[12]+1.369306393762915*magB2[11]-1.837117307087383*magB2[10]+0.7905694150420947*magB2[9]+0.7905694150420947*magB2[8]+0.7905694150420947*magB2[7]-1.060660171779821*magB2[6]-1.060660171779821*magB2[5]+1.060660171779821*magB2[4]-0.6123724356957944*magB2[3]+0.6123724356957944*magB2[2]+0.6123724356957944*magB2[1]+0.3535533905932737*magB2[0] < 0.0) cell_avg = 1; 
  if ((-1.185854122563142*magB2[19])+0.6846531968814574*magB2[16]+0.6846531968814574*magB2[15]-1.369306393762915*magB2[12]-1.369306393762915*magB2[11]-0.3952847075210473*magB2[9]+0.7905694150420947*magB2[8]+0.7905694150420947*magB2[7]+1.060660171779821*magB2[4]-0.6123724356957944*magB2[2]-0.6123724356957944*magB2[1]+0.3535533905932737*magB2[0] < 0.0) cell_avg = 1; 
  if (1.185854122563142*magB2[19]+0.6846531968814574*magB2[16]-0.6846531968814574*magB2[15]+1.369306393762915*magB2[12]-1.369306393762915*magB2[11]-0.3952847075210473*magB2[9]+0.7905694150420947*magB2[8]+0.7905694150420947*magB2[7]-1.060660171779821*magB2[4]-0.6123724356957944*magB2[2]+0.6123724356957944*magB2[1]+0.3535533905932737*magB2[0] < 0.0) cell_avg = 1; 
  if (1.185854122563142*magB2[19]-0.6846531968814574*magB2[16]+0.6846531968814574*magB2[15]-1.369306393762915*magB2[12]+1.369306393762915*magB2[11]-0.3952847075210473*magB2[9]+0.7905694150420947*magB2[8]+0.7905694150420947*magB2[7]-1.060660171779821*magB2[4]+0.6123724356957944*magB2[2]-0.6123724356957944*magB2[1]+0.3535533905932737*magB2[0] < 0.0) cell_avg = 1; 
  if ((-1.185854122563142*magB2[19])-0.6846531968814574*magB2[16]-0.6846531968814574*magB2[15]+1.369306393762915*magB2[12]+1.369306393762915*magB2[11]-0.3952847075210473*magB2[9]+0.7905694150420947*magB2[8]+0.7905694150420947*magB2[7]+1.060660171779821*magB2[4]+0.6123724356957944*magB2[2]+0.6123724356957944*magB2[1]+0.3535533905932737*magB2[0] < 0.0) cell_avg = 1; 
  if (2.371708245126284*magB2[19]-2.371708245126284*magB2[18]-2.371708245126284*magB2[17]-1.369306393762915*magB2[16]-1.369306393762915*magB2[15]+1.369306393762915*magB2[14]+1.369306393762915*magB2[13]-1.369306393762915*magB2[12]-1.369306393762915*magB2[11]+1.837117307087383*magB2[10]+0.7905694150420947*magB2[9]+0.7905694150420947*magB2[8]+0.7905694150420947*magB2[7]-1.060660171779821*magB2[6]-1.060660171779821*magB2[5]+1.060660171779821*magB2[4]+0.6123724356957944*magB2[3]-0.6123724356957944*magB2[2]-0.6123724356957944*magB2[1]+0.3535533905932737*magB2[0] < 0.0) cell_avg = 1; 
  if (1.185854122563142*magB2[17]-1.369306393762915*magB2[16]+1.369306393762915*magB2[14]-0.6846531968814574*magB2[13]+0.6846531968814574*magB2[11]+0.7905694150420947*magB2[9]+0.7905694150420947*magB2[8]-0.3952847075210473*magB2[7]-1.060660171779821*magB2[6]+0.6123724356957944*magB2[3]-0.6123724356957944*magB2[2]+0.3535533905932737*magB2[0] < 0.0) cell_avg = 1; 
  if ((-2.371708245126284*magB2[19])+2.371708245126284*magB2[18]-2.371708245126284*magB2[17]-1.369306393762915*magB2[16]+1.369306393762915*magB2[15]+1.369306393762915*magB2[14]+1.369306393762915*magB2[13]+1.369306393762915*magB2[12]-1.369306393762915*magB2[11]-1.837117307087383*magB2[10]+0.7905694150420947*magB2[9]+0.7905694150420947*magB2[8]+0.7905694150420947*magB2[7]-1.060660171779821*magB2[6]+1.060660171779821*magB2[5]-1.060660171779821*magB2[4]+0.6123724356957944*magB2[3]-0.6123724356957944*magB2[2]+0.6123724356957944*magB2[1]+0.3535533905932737*magB2[0] < 0.0) cell_avg = 1; 
  if (1.185854122563142*magB2[18]-1.369306393762915*magB2[15]-0.6846531968814574*magB2[14]+1.369306393762915*magB2[13]+0.6846531968814574*magB2[12]+0.7905694150420947*magB2[9]-0.3952847075210473*magB2[8]+0.7905694150420947*magB2[7]-1.060660171779821*magB2[5]+0.6123724356957944*magB2[3]-0.6123724356957944*magB2[1]+0.3535533905932737*magB2[0] < 0.0) cell_avg = 1; 
  if ((-1.185854122563142*magB2[18])+1.369306393762915*magB2[15]-0.6846531968814574*magB2[14]+1.369306393762915*magB2[13]-0.6846531968814574*magB2[12]+0.7905694150420947*magB2[9]-0.3952847075210473*magB2[8]+0.7905694150420947*magB2[7]+1.060660171779821*magB2[5]+0.6123724356957944*magB2[3]+0.6123724356957944*magB2[1]+0.3535533905932737*magB2[0] < 0.0) cell_avg = 1; 
  if ((-2.371708245126284*magB2[19])-2.371708245126284*magB2[18]+2.371708245126284*magB2[17]+1.369306393762915*magB2[16]-1.369306393762915*magB2[15]+1.369306393762915*magB2[14]+1.369306393762915*magB2[13]-1.369306393762915*magB2[12]+1.369306393762915*magB2[11]-1.837117307087383*magB2[10]+0.7905694150420947*magB2[9]+0.7905694150420947*magB2[8]+0.7905694150420947*magB2[7]+1.060660171779821*magB2[6]-1.060660171779821*magB2[5]-1.060660171779821*magB2[4]+0.6123724356957944*magB2[3]+0.6123724356957944*magB2[2]-0.6123724356957944*magB2[1]+0.3535533905932737*magB2[0] < 0.0) cell_avg = 1; 
  if ((-1.185854122563142*magB2[17])+1.369306393762915*magB2[16]+1.369306393762915*magB2[14]-0.6846531968814574*magB2[13]-0.6846531968814574*magB2[11]+0.7905694150420947*magB2[9]+0.7905694150420947*magB2[8]-0.3952847075210473*magB2[7]+1.060660171779821*magB2[6]+0.6123724356957944*magB2[3]+0.6123724356957944*magB2[2]+0.3535533905932737*magB2[0] < 0.0) cell_avg = 1; 
  if (2.371708245126284*magB2[19]+2.371708245126284*magB2[18]+2.371708245126284*magB2[17]+1.369306393762915*magB2[16]+1.369306393762915*magB2[15]+1.369306393762915*magB2[14]+1.369306393762915*magB2[13]+1.369306393762915*magB2[12]+1.369306393762915*magB2[11]+1.837117307087383*magB2[10]+0.7905694150420947*magB2[9]+0.7905694150420947*magB2[8]+0.7905694150420947*magB2[7]+1.060660171779821*magB2[6]+1.060660171779821*magB2[5]+1.060660171779821*magB2[4]+0.6123724356957944*magB2[3]+0.6123724356957944*magB2[2]+0.6123724356957944*magB2[1]+0.3535533905932737*magB2[0] < 0.0) cell_avg = 1; 
  if (cell_avg) { 
    magB2[1] = 0.0; 
    magB2[2] = 0.0; 
    magB2[3] = 0.0; 
    magB2[4] = 0.0; 
    magB2[5] = 0.0; 
    magB2[6] = 0.0; 
    magB2[7] = 0.0; 
    magB2[8] = 0.0; 
    magB2[9] = 0.0; 
    magB2[10] = 0.0; 
    magB2[11] = 0.0; 
    magB2[12] = 0.0; 
    magB2[13] = 0.0; 
    magB2[14] = 0.0; 
    magB2[15] = 0.0; 
    magB2[16] = 0.0; 
    magB2[17] = 0.0; 
    magB2[18] = 0.0; 
    magB2[19] = 0.0; 
  } 
  gkyl_mat_set(rhs,0,0,2.828427124746191); 
  gkyl_mat_set(A,0,0,0.3535533905932737*magB2[0]); 
  gkyl_mat_set(A,0,1,0.3535533905932737*magB2[1]); 
  gkyl_mat_set(A,0,2,0.3535533905932737*magB2[2]); 
  gkyl_mat_set(A,0,3,0.3535533905932737*magB2[3]); 
  gkyl_mat_set(A,0,4,0.3535533905932737*magB2[4]); 
  gkyl_mat_set(A,0,5,0.3535533905932737*magB2[5]); 
  gkyl_mat_set(A,0,6,0.3535533905932737*magB2[6]); 
  gkyl_mat_set(A,0,7,0.3535533905932737*magB2[7]); 
  gkyl_mat_set(A,0,8,0.3535533905932737*magB2[8]); 
  gkyl_mat_set(A,0,9,0.3535533905932737*magB2[9]); 
  gkyl_mat_set(A,0,10,0.3535533905932737*magB2[10]); 
  gkyl_mat_set(A,0,11,0.3535533905932737*magB2[11]); 
  gkyl_mat_set(A,0,12,0.3535533905932737*magB2[12]); 
  gkyl_mat_set(A,0,13,0.3535533905932737*magB2[13]); 
  gkyl_mat_set(A,0,14,0.3535533905932737*magB2[14]); 
  gkyl_mat_set(A,0,15,0.3535533905932737*magB2[15]); 
  gkyl_mat_set(A,0,16,0.3535533905932737*magB2[16]); 
  gkyl_mat_set(A,0,17,0.3535533905932737*magB2[17]); 
  gkyl_mat_set(A,0,18,0.3535533905932737*magB2[18]); 
  gkyl_mat_set(A,0,19,0.3535533905932737*magB2[19]); 
  gkyl_mat_set(A,1,0,0.3535533905932737*magB2[1]); 
  gkyl_mat_set(A,1,1,0.3162277660168379*magB2[7]+0.3535533905932737*magB2[0]); 
  gkyl_mat_set(A,1,2,0.3535533905932737*magB2[4]); 
  gkyl_mat_set(A,1,3,0.3535533905932737*magB2[5]); 
  gkyl_mat_set(A,1,4,0.3162277660168379*magB2[11]+0.3535533905932737*magB2[2]); 
  gkyl_mat_set(A,1,5,0.3162277660168379*magB2[13]+0.3535533905932737*magB2[3]); 
  gkyl_mat_set(A,1,6,0.3535533905932737*magB2[10]); 
  gkyl_mat_set(A,1,7,0.3162277660168379*magB2[1]); 
  gkyl_mat_set(A,1,8,0.3535533905932737*magB2[12]); 
  gkyl_mat_set(A,1,9,0.3535533905932737*magB2[15]); 
  gkyl_mat_set(A,1,10,0.3162277660168379*magB2[17]+0.3535533905932737*magB2[6]); 
  gkyl_mat_set(A,1,11,0.3162277660168379*magB2[4]); 
  gkyl_mat_set(A,1,12,0.3535533905932737*magB2[8]); 
  gkyl_mat_set(A,1,13,0.3162277660168379*magB2[5]); 
  gkyl_mat_set(A,1,14,0.3535533905932737*magB2[18]); 
  gkyl_mat_set(A,1,15,0.3535533905932737*magB2[9]); 
  gkyl_mat_set(A,1,16,0.3535533905932737*magB2[19]); 
  gkyl_mat_set(A,1,17,0.3162277660168379*magB2[10]); 
  gkyl_mat_set(A,1,18,0.3535533905932737*magB2[14]); 
  gkyl_mat_set(A,1,19,0.3535533905932737*magB2[16]); 
  gkyl_mat_set(A,2,0,0.3535533905932737*magB2[2]); 
  gkyl_mat_set(A,2,1,0.3535533905932737*magB2[4]); 
  gkyl_mat_set(A,2,2,0.3162277660168379*magB2[8]+0.3535533905932737*magB2[0]); 
  gkyl_mat_set(A,2,3,0.3535533905932737*magB2[6]); 
  gkyl_mat_set(A,2,4,0.3162277660168379*magB2[12]+0.3535533905932737*magB2[1]); 
  gkyl_mat_set(A,2,5,0.3535533905932737*magB2[10]); 
  gkyl_mat_set(A,2,6,0.3162277660168379*magB2[14]+0.3535533905932737*magB2[3]); 
  gkyl_mat_set(A,2,7,0.3535533905932737*magB2[11]); 
  gkyl_mat_set(A,2,8,0.3162277660168379*magB2[2]); 
  gkyl_mat_set(A,2,9,0.3535533905932737*magB2[16]); 
  gkyl_mat_set(A,2,10,0.3162277660168379*magB2[18]+0.3535533905932737*magB2[5]); 
  gkyl_mat_set(A,2,11,0.3535533905932737*magB2[7]); 
  gkyl_mat_set(A,2,12,0.3162277660168379*magB2[4]); 
  gkyl_mat_set(A,2,13,0.3535533905932737*magB2[17]); 
  gkyl_mat_set(A,2,14,0.3162277660168379*magB2[6]); 
  gkyl_mat_set(A,2,15,0.3535533905932737*magB2[19]); 
  gkyl_mat_set(A,2,16,0.3535533905932737*magB2[9]); 
  gkyl_mat_set(A,2,17,0.3535533905932737*magB2[13]); 
  gkyl_mat_set(A,2,18,0.3162277660168379*magB2[10]); 
  gkyl_mat_set(A,2,19,0.3535533905932737*magB2[15]); 
  gkyl_mat_set(A,3,0,0.3535533905932737*magB2[3]); 
  gkyl_mat_set(A,3,1,0.3535533905932737*magB2[5]); 
  gkyl_mat_set(A,3,2,0.3535533905932737*magB2[6]); 
  gkyl_mat_set(A,3,3,0.3162277660168379*magB2[9]+0.3535533905932737*magB2[0]); 
  gkyl_mat_set(A,3,4,0.3535533905932737*magB2[10]); 
  gkyl_mat_set(A,3,5,0.3162277660168379*magB2[15]+0.3535533905932737*magB2[1]); 
  gkyl_mat_set(A,3,6,0.3162277660168379*magB2[16]+0.3535533905932737*magB2[2]); 
  gkyl_mat_set(A,3,7,0.3535533905932737*magB2[13]); 
  gkyl_mat_set(A,3,8,0.3535533905932737*magB2[14]); 
  gkyl_mat_set(A,3,9,0.3162277660168379*magB2[3]); 
  gkyl_mat_set(A,3,10,0.3162277660168379*magB2[19]+0.3535533905932737*magB2[4]); 
  gkyl_mat_set(A,3,11,0.3535533905932737*magB2[17]); 
  gkyl_mat_set(A,3,12,0.3535533905932737*magB2[18]); 
  gkyl_mat_set(A,3,13,0.3535533905932737*magB2[7]); 
  gkyl_mat_set(A,3,14,0.3535533905932737*magB2[8]); 
  gkyl_mat_set(A,3,15,0.3162277660168379*magB2[5]); 
  gkyl_mat_set(A,3,16,0.3162277660168379*magB2[6]); 
  gkyl_mat_set(A,3,17,0.3535533905932737*magB2[11]); 
  gkyl_mat_set(A,3,18,0.3535533905932737*magB2[12]); 
  gkyl_mat_set(A,3,19,0.3162277660168379*magB2[10]); 
  gkyl_mat_set(A,4,0,0.3535533905932737*magB2[4]); 
  gkyl_mat_set(A,4,1,0.3162277660168379*magB2[11]+0.3535533905932737*magB2[2]); 
  gkyl_mat_set(A,4,2,0.3162277660168379*magB2[12]+0.3535533905932737*magB2[1]); 
  gkyl_mat_set(A,4,3,0.3535533905932737*magB2[10]); 
  gkyl_mat_set(A,4,4,0.3162277660168379*magB2[8]+0.3162277660168379*magB2[7]+0.3535533905932737*magB2[0]); 
  gkyl_mat_set(A,4,5,0.3162277660168379*magB2[17]+0.3535533905932737*magB2[6]); 
  gkyl_mat_set(A,4,6,0.3162277660168379*magB2[18]+0.3535533905932737*magB2[5]); 
  gkyl_mat_set(A,4,7,0.3162277660168379*magB2[4]); 
  gkyl_mat_set(A,4,8,0.3162277660168379*magB2[4]); 
  gkyl_mat_set(A,4,9,0.3535533905932737*magB2[19]); 
  gkyl_mat_set(A,4,10,0.3162277660168379*magB2[14]+0.3162277660168379*magB2[13]+0.3535533905932737*magB2[3]); 
  gkyl_mat_set(A,4,11,0.2828427124746191*magB2[12]+0.3162277660168379*magB2[1]); 
  gkyl_mat_set(A,4,12,0.2828427124746191*magB2[11]+0.3162277660168379*magB2[2]); 
  gkyl_mat_set(A,4,13,0.3162277660168379*magB2[10]); 
  gkyl_mat_set(A,4,14,0.3162277660168379*magB2[10]); 
  gkyl_mat_set(A,4,15,0.3535533905932737*magB2[16]); 
  gkyl_mat_set(A,4,16,0.3535533905932737*magB2[15]); 
  gkyl_mat_set(A,4,17,0.2828427124746191*magB2[18]+0.3162277660168379*magB2[5]); 
  gkyl_mat_set(A,4,18,0.2828427124746191*magB2[17]+0.3162277660168379*magB2[6]); 
  gkyl_mat_set(A,4,19,0.3535533905932737*magB2[9]); 
  gkyl_mat_set(A,5,0,0.3535533905932737*magB2[5]); 
  gkyl_mat_set(A,5,1,0.3162277660168379*magB2[13]+0.3535533905932737*magB2[3]); 
  gkyl_mat_set(A,5,2,0.3535533905932737*magB2[10]); 
  gkyl_mat_set(A,5,3,0.3162277660168379*magB2[15]+0.3535533905932737*magB2[1]); 
  gkyl_mat_set(A,5,4,0.3162277660168379*magB2[17]+0.3535533905932737*magB2[6]); 
  gkyl_mat_set(A,5,5,0.3162277660168379*magB2[9]+0.3162277660168379*magB2[7]+0.3535533905932737*magB2[0]); 
  gkyl_mat_set(A,5,6,0.3162277660168379*magB2[19]+0.3535533905932737*magB2[4]); 
  gkyl_mat_set(A,5,7,0.3162277660168379*magB2[5]); 
  gkyl_mat_set(A,5,8,0.3535533905932737*magB2[18]); 
  gkyl_mat_set(A,5,9,0.3162277660168379*magB2[5]); 
  gkyl_mat_set(A,5,10,0.3162277660168379*magB2[16]+0.3162277660168379*magB2[11]+0.3535533905932737*magB2[2]); 
  gkyl_mat_set(A,5,11,0.3162277660168379*magB2[10]); 
  gkyl_mat_set(A,5,12,0.3535533905932737*magB2[14]); 
  gkyl_mat_set(A,5,13,0.2828427124746191*magB2[15]+0.3162277660168379*magB2[1]); 
  gkyl_mat_set(A,5,14,0.3535533905932737*magB2[12]); 
  gkyl_mat_set(A,5,15,0.2828427124746191*magB2[13]+0.3162277660168379*magB2[3]); 
  gkyl_mat_set(A,5,16,0.3162277660168379*magB2[10]); 
  gkyl_mat_set(A,5,17,0.2828427124746191*magB2[19]+0.3162277660168379*magB2[4]); 
  gkyl_mat_set(A,5,18,0.3535533905932737*magB2[8]); 
  gkyl_mat_set(A,5,19,0.2828427124746191*magB2[17]+0.3162277660168379*magB2[6]); 
  gkyl_mat_set(A,6,0,0.3535533905932737*magB2[6]); 
  gkyl_mat_set(A,6,1,0.3535533905932737*magB2[10]); 
  gkyl_mat_set(A,6,2,0.3162277660168379*magB2[14]+0.3535533905932737*magB2[3]); 
  gkyl_mat_set(A,6,3,0.3162277660168379*magB2[16]+0.3535533905932737*magB2[2]); 
  gkyl_mat_set(A,6,4,0.3162277660168379*magB2[18]+0.3535533905932737*magB2[5]); 
  gkyl_mat_set(A,6,5,0.3162277660168379*magB2[19]+0.3535533905932737*magB2[4]); 
  gkyl_mat_set(A,6,6,0.3162277660168379*magB2[9]+0.3162277660168379*magB2[8]+0.3535533905932737*magB2[0]); 
  gkyl_mat_set(A,6,7,0.3535533905932737*magB2[17]); 
  gkyl_mat_set(A,6,8,0.3162277660168379*magB2[6]); 
  gkyl_mat_set(A,6,9,0.3162277660168379*magB2[6]); 
  gkyl_mat_set(A,6,10,0.3162277660168379*magB2[15]+0.3162277660168379*magB2[12]+0.3535533905932737*magB2[1]); 
  gkyl_mat_set(A,6,11,0.3535533905932737*magB2[13]); 
  gkyl_mat_set(A,6,12,0.3162277660168379*magB2[10]); 
  gkyl_mat_set(A,6,13,0.3535533905932737*magB2[11]); 
  gkyl_mat_set(A,6,14,0.2828427124746191*magB2[16]+0.3162277660168379*magB2[2]); 
  gkyl_mat_set(A,6,15,0.3162277660168379*magB2[10]); 
  gkyl_mat_set(A,6,16,0.2828427124746191*magB2[14]+0.3162277660168379*magB2[3]); 
  gkyl_mat_set(A,6,17,0.3535533905932737*magB2[7]); 
  gkyl_mat_set(A,6,18,0.2828427124746191*magB2[19]+0.3162277660168379*magB2[4]); 
  gkyl_mat_set(A,6,19,0.2828427124746191*magB2[18]+0.3162277660168379*magB2[5]); 
  gkyl_mat_set(A,7,0,0.3535533905932737*magB2[7]); 
  gkyl_mat_set(A,7,1,0.3162277660168379*magB2[1]); 
  gkyl_mat_set(A,7,2,0.3535533905932737*magB2[11]); 
  gkyl_mat_set(A,7,3,0.3535533905932737*magB2[13]); 
  gkyl_mat_set(A,7,4,0.3162277660168379*magB2[4]); 
  gkyl_mat_set(A,7,5,0.3162277660168379*magB2[5]); 
  gkyl_mat_set(A,7,6,0.3535533905932737*magB2[17]); 
  gkyl_mat_set(A,7,7,0.2258769757263128*magB2[7]+0.3535533905932737*magB2[0]); 
  gkyl_mat_set(A,7,8,0.0); 
  gkyl_mat_set(A,7,9,0.0); 
  gkyl_mat_set(A,7,10,0.3162277660168379*magB2[10]); 
  gkyl_mat_set(A,7,11,0.2258769757263128*magB2[11]+0.3535533905932737*magB2[2]); 
  gkyl_mat_set(A,7,12,0.3162277660168379*magB2[12]); 
  gkyl_mat_set(A,7,13,0.2258769757263128*magB2[13]+0.3535533905932737*magB2[3]); 
  gkyl_mat_set(A,7,14,0.0); 
  gkyl_mat_set(A,7,15,0.3162277660168379*magB2[15]); 
  gkyl_mat_set(A,7,16,0.0); 
  gkyl_mat_set(A,7,17,0.2258769757263128*magB2[17]+0.3535533905932737*magB2[6]); 
  gkyl_mat_set(A,7,18,0.3162277660168379*magB2[18]); 
  gkyl_mat_set(A,7,19,0.3162277660168379*magB2[19]); 
  gkyl_mat_set(A,8,0,0.3535533905932737*magB2[8]); 
  gkyl_mat_set(A,8,1,0.3535533905932737*magB2[12]); 
  gkyl_mat_set(A,8,2,0.3162277660168379*magB2[2]); 
  gkyl_mat_set(A,8,3,0.3535533905932737*magB2[14]); 
  gkyl_mat_set(A,8,4,0.3162277660168379*magB2[4]); 
  gkyl_mat_set(A,8,5,0.3535533905932737*magB2[18]); 
  gkyl_mat_set(A,8,6,0.3162277660168379*magB2[6]); 
  gkyl_mat_set(A,8,7,0.0); 
  gkyl_mat_set(A,8,8,0.2258769757263128*magB2[8]+0.3535533905932737*magB2[0]); 
  gkyl_mat_set(A,8,9,0.0); 
  gkyl_mat_set(A,8,10,0.3162277660168379*magB2[10]); 
  gkyl_mat_set(A,8,11,0.3162277660168379*magB2[11]); 
  gkyl_mat_set(A,8,12,0.2258769757263128*magB2[12]+0.3535533905932737*magB2[1]); 
  gkyl_mat_set(A,8,13,0.0); 
  gkyl_mat_set(A,8,14,0.2258769757263128*magB2[14]+0.3535533905932737*magB2[3]); 
  gkyl_mat_set(A,8,15,0.0); 
  gkyl_mat_set(A,8,16,0.3162277660168379*magB2[16]); 
  gkyl_mat_set(A,8,17,0.3162277660168379*magB2[17]); 
  gkyl_mat_set(A,8,18,0.2258769757263128*magB2[18]+0.3535533905932737*magB2[5]); 
  gkyl_mat_set(A,8,19,0.3162277660168379*magB2[19]); 
  gkyl_mat_set(A,9,0,0.3535533905932737*magB2[9]); 
  gkyl_mat_set(A,9,1,0.3535533905932737*magB2[15]); 
  gkyl_mat_set(A,9,2,0.3535533905932737*magB2[16]); 
  gkyl_mat_set(A,9,3,0.3162277660168379*magB2[3]); 
  gkyl_mat_set(A,9,4,0.3535533905932737*magB2[19]); 
  gkyl_mat_set(A,9,5,0.3162277660168379*magB2[5]); 
  gkyl_mat_set(A,9,6,0.3162277660168379*magB2[6]); 
  gkyl_mat_set(A,9,7,0.0); 
  gkyl_mat_set(A,9,8,0.0); 
  gkyl_mat_set(A,9,9,0.2258769757263128*magB2[9]+0.3535533905932737*magB2[0]); 
  gkyl_mat_set(A,9,10,0.3162277660168379*magB2[10]); 
  gkyl_mat_set(A,9,11,0.0); 
  gkyl_mat_set(A,9,12,0.0); 
  gkyl_mat_set(A,9,13,0.3162277660168379*magB2[13]); 
  gkyl_mat_set(A,9,14,0.3162277660168379*magB2[14]); 
  gkyl_mat_set(A,9,15,0.2258769757263128*magB2[15]+0.3535533905932737*magB2[1]); 
  gkyl_mat_set(A,9,16,0.2258769757263128*magB2[16]+0.3535533905932737*magB2[2]); 
  gkyl_mat_set(A,9,17,0.3162277660168379*magB2[17]); 
  gkyl_mat_set(A,9,18,0.3162277660168379*magB2[18]); 
  gkyl_mat_set(A,9,19,0.2258769757263128*magB2[19]+0.3535533905932737*magB2[4]); 
  gkyl_mat_set(A,10,0,0.3535533905932737*magB2[10]); 
  gkyl_mat_set(A,10,1,0.3162277660168379*magB2[17]+0.3535533905932737*magB2[6]); 
  gkyl_mat_set(A,10,2,0.3162277660168379*magB2[18]+0.3535533905932737*magB2[5]); 
  gkyl_mat_set(A,10,3,0.3162277660168379*magB2[19]+0.3535533905932737*magB2[4]); 
  gkyl_mat_set(A,10,4,0.3162277660168379*magB2[14]+0.3162277660168379*magB2[13]+0.3535533905932737*magB2[3]); 
  gkyl_mat_set(A,10,5,0.3162277660168379*magB2[16]+0.3162277660168379*magB2[11]+0.3535533905932737*magB2[2]); 
  gkyl_mat_set(A,10,6,0.3162277660168379*magB2[15]+0.3162277660168379*magB2[12]+0.3535533905932737*magB2[1]); 
  gkyl_mat_set(A,10,7,0.3162277660168379*magB2[10]); 
  gkyl_mat_set(A,10,8,0.3162277660168379*magB2[10]); 
  gkyl_mat_set(A,10,9,0.3162277660168379*magB2[10]); 
  gkyl_mat_set(A,10,10,0.3162277660168379*magB2[9]+0.3162277660168379*magB2[8]+0.3162277660168379*magB2[7]+0.3535533905932737*magB2[0]); 
  gkyl_mat_set(A,10,11,0.282842712474619*magB2[18]+0.3162277660168379*magB2[5]); 
  gkyl_mat_set(A,10,12,0.282842712474619*magB2[17]+0.3162277660168379*magB2[6]); 
  gkyl_mat_set(A,10,13,0.282842712474619*magB2[19]+0.3162277660168379*magB2[4]); 
  gkyl_mat_set(A,10,14,0.282842712474619*magB2[19]+0.3162277660168379*magB2[4]); 
  gkyl_mat_set(A,10,15,0.282842712474619*magB2[17]+0.3162277660168379*magB2[6]); 
  gkyl_mat_set(A,10,16,0.282842712474619*magB2[18]+0.3162277660168379*magB2[5]); 
  gkyl_mat_set(A,10,17,0.282842712474619*magB2[15]+0.282842712474619*magB2[12]+0.3162277660168379*magB2[1]); 
  gkyl_mat_set(A,10,18,0.282842712474619*magB2[16]+0.282842712474619*magB2[11]+0.3162277660168379*magB2[2]); 
  gkyl_mat_set(A,10,19,0.282842712474619*magB2[14]+0.282842712474619*magB2[13]+0.3162277660168379*magB2[3]); 
  gkyl_mat_set(A,11,0,0.3535533905932737*magB2[11]); 
  gkyl_mat_set(A,11,1,0.3162277660168379*magB2[4]); 
  gkyl_mat_set(A,11,2,0.3535533905932737*magB2[7]); 
  gkyl_mat_set(A,11,3,0.3535533905932737*magB2[17]); 
  gkyl_mat_set(A,11,4,0.2828427124746191*magB2[12]+0.3162277660168379*magB2[1]); 
  gkyl_mat_set(A,11,5,0.3162277660168379*magB2[10]); 
  gkyl_mat_set(A,11,6,0.3535533905932737*magB2[13]); 
  gkyl_mat_set(A,11,7,0.2258769757263128*magB2[11]+0.3535533905932737*magB2[2]); 
  gkyl_mat_set(A,11,8,0.3162277660168379*magB2[11]); 
  gkyl_mat_set(A,11,9,0.0); 
  gkyl_mat_set(A,11,10,0.282842712474619*magB2[18]+0.3162277660168379*magB2[5]); 
  gkyl_mat_set(A,11,11,0.3162277660168379*magB2[8]+0.2258769757263128*magB2[7]+0.3535533905932737*magB2[0]); 
  gkyl_mat_set(A,11,12,0.2828427124746191*magB2[4]); 
  gkyl_mat_set(A,11,13,0.2258769757263128*magB2[17]+0.3535533905932737*magB2[6]); 
  gkyl_mat_set(A,11,14,0.3162277660168379*magB2[17]); 
  gkyl_mat_set(A,11,15,0.3162277660168379*magB2[19]); 
  gkyl_mat_set(A,11,16,0.0); 
  gkyl_mat_set(A,11,17,0.3162277660168379*magB2[14]+0.2258769757263128*magB2[13]+0.3535533905932737*magB2[3]); 
  gkyl_mat_set(A,11,18,0.282842712474619*magB2[10]); 
  gkyl_mat_set(A,11,19,0.3162277660168379*magB2[15]); 
  gkyl_mat_set(A,12,0,0.3535533905932737*magB2[12]); 
  gkyl_mat_set(A,12,1,0.3535533905932737*magB2[8]); 
  gkyl_mat_set(A,12,2,0.3162277660168379*magB2[4]); 
  gkyl_mat_set(A,12,3,0.3535533905932737*magB2[18]); 
  gkyl_mat_set(A,12,4,0.2828427124746191*magB2[11]+0.3162277660168379*magB2[2]); 
  gkyl_mat_set(A,12,5,0.3535533905932737*magB2[14]); 
  gkyl_mat_set(A,12,6,0.3162277660168379*magB2[10]); 
  gkyl_mat_set(A,12,7,0.3162277660168379*magB2[12]); 
  gkyl_mat_set(A,12,8,0.2258769757263128*magB2[12]+0.3535533905932737*magB2[1]); 
  gkyl_mat_set(A,12,9,0.0); 
  gkyl_mat_set(A,12,10,0.282842712474619*magB2[17]+0.3162277660168379*magB2[6]); 
  gkyl_mat_set(A,12,11,0.2828427124746191*magB2[4]); 
  gkyl_mat_set(A,12,12,0.2258769757263128*magB2[8]+0.3162277660168379*magB2[7]+0.3535533905932737*magB2[0]); 
  gkyl_mat_set(A,12,13,0.3162277660168379*magB2[18]); 
  gkyl_mat_set(A,12,14,0.2258769757263128*magB2[18]+0.3535533905932737*magB2[5]); 
  gkyl_mat_set(A,12,15,0.0); 
  gkyl_mat_set(A,12,16,0.3162277660168379*magB2[19]); 
  gkyl_mat_set(A,12,17,0.282842712474619*magB2[10]); 
  gkyl_mat_set(A,12,18,0.2258769757263128*magB2[14]+0.3162277660168379*magB2[13]+0.3535533905932737*magB2[3]); 
  gkyl_mat_set(A,12,19,0.3162277660168379*magB2[16]); 
  gkyl_mat_set(A,13,0,0.3535533905932737*magB2[13]); 
  gkyl_mat_set(A,13,1,0.3162277660168379*magB2[5]); 
  gkyl_mat_set(A,13,2,0.3535533905932737*magB2[17]); 
  gkyl_mat_set(A,13,3,0.3535533905932737*magB2[7]); 
  gkyl_mat_set(A,13,4,0.3162277660168379*magB2[10]); 
  gkyl_mat_set(A,13,5,0.2828427124746191*magB2[15]+0.3162277660168379*magB2[1]); 
  gkyl_mat_set(A,13,6,0.3535533905932737*magB2[11]); 
  gkyl_mat_set(A,13,7,0.2258769757263128*magB2[13]+0.3535533905932737*magB2[3]); 
  gkyl_mat_set(A,13,8,0.0); 
  gkyl_mat_set(A,13,9,0.3162277660168379*magB2[13]); 
  gkyl_mat_set(A,13,10,0.282842712474619*magB2[19]+0.3162277660168379*magB2[4]); 
  gkyl_mat_set(A,13,11,0.2258769757263128*magB2[17]+0.3535533905932737*magB2[6]); 
  gkyl_mat_set(A,13,12,0.3162277660168379*magB2[18]); 
  gkyl_mat_set(A,13,13,0.3162277660168379*magB2[9]+0.2258769757263128*magB2[7]+0.3535533905932737*magB2[0]); 
  gkyl_mat_set(A,13,14,0.0); 
  gkyl_mat_set(A,13,15,0.2828427124746191*magB2[5]); 
  gkyl_mat_set(A,13,16,0.3162277660168379*magB2[17]); 
  gkyl_mat_set(A,13,17,0.3162277660168379*magB2[16]+0.2258769757263128*magB2[11]+0.3535533905932737*magB2[2]); 
  gkyl_mat_set(A,13,18,0.3162277660168379*magB2[12]); 
  gkyl_mat_set(A,13,19,0.282842712474619*magB2[10]); 
  gkyl_mat_set(A,14,0,0.3535533905932737*magB2[14]); 
  gkyl_mat_set(A,14,1,0.3535533905932737*magB2[18]); 
  gkyl_mat_set(A,14,2,0.3162277660168379*magB2[6]); 
  gkyl_mat_set(A,14,3,0.3535533905932737*magB2[8]); 
  gkyl_mat_set(A,14,4,0.3162277660168379*magB2[10]); 
  gkyl_mat_set(A,14,5,0.3535533905932737*magB2[12]); 
  gkyl_mat_set(A,14,6,0.2828427124746191*magB2[16]+0.3162277660168379*magB2[2]); 
  gkyl_mat_set(A,14,7,0.0); 
  gkyl_mat_set(A,14,8,0.2258769757263128*magB2[14]+0.3535533905932737*magB2[3]); 
  gkyl_mat_set(A,14,9,0.3162277660168379*magB2[14]); 
  gkyl_mat_set(A,14,10,0.282842712474619*magB2[19]+0.3162277660168379*magB2[4]); 
  gkyl_mat_set(A,14,11,0.3162277660168379*magB2[17]); 
  gkyl_mat_set(A,14,12,0.2258769757263128*magB2[18]+0.3535533905932737*magB2[5]); 
  gkyl_mat_set(A,14,13,0.0); 
  gkyl_mat_set(A,14,14,0.3162277660168379*magB2[9]+0.2258769757263128*magB2[8]+0.3535533905932737*magB2[0]); 
  gkyl_mat_set(A,14,15,0.3162277660168379*magB2[18]); 
  gkyl_mat_set(A,14,16,0.2828427124746191*magB2[6]); 
  gkyl_mat_set(A,14,17,0.3162277660168379*magB2[11]); 
  gkyl_mat_set(A,14,18,0.3162277660168379*magB2[15]+0.2258769757263128*magB2[12]+0.3535533905932737*magB2[1]); 
  gkyl_mat_set(A,14,19,0.282842712474619*magB2[10]); 
  gkyl_mat_set(A,15,0,0.3535533905932737*magB2[15]); 
  gkyl_mat_set(A,15,1,0.3535533905932737*magB2[9]); 
  gkyl_mat_set(A,15,2,0.3535533905932737*magB2[19]); 
  gkyl_mat_set(A,15,3,0.3162277660168379*magB2[5]); 
  gkyl_mat_set(A,15,4,0.3535533905932737*magB2[16]); 
  gkyl_mat_set(A,15,5,0.2828427124746191*magB2[13]+0.3162277660168379*magB2[3]); 
  gkyl_mat_set(A,15,6,0.3162277660168379*magB2[10]); 
  gkyl_mat_set(A,15,7,0.3162277660168379*magB2[15]); 
  gkyl_mat_set(A,15,8,0.0); 
  gkyl_mat_set(A,15,9,0.2258769757263128*magB2[15]+0.3535533905932737*magB2[1]); 
  gkyl_mat_set(A,15,10,0.282842712474619*magB2[17]+0.3162277660168379*magB2[6]); 
  gkyl_mat_set(A,15,11,0.3162277660168379*magB2[19]); 
  gkyl_mat_set(A,15,12,0.0); 
  gkyl_mat_set(A,15,13,0.2828427124746191*magB2[5]); 
  gkyl_mat_set(A,15,14,0.3162277660168379*magB2[18]); 
  gkyl_mat_set(A,15,15,0.2258769757263128*magB2[9]+0.3162277660168379*magB2[7]+0.3535533905932737*magB2[0]); 
  gkyl_mat_set(A,15,16,0.2258769757263128*magB2[19]+0.3535533905932737*magB2[4]); 
  gkyl_mat_set(A,15,17,0.282842712474619*magB2[10]); 
  gkyl_mat_set(A,15,18,0.3162277660168379*magB2[14]); 
  gkyl_mat_set(A,15,19,0.2258769757263128*magB2[16]+0.3162277660168379*magB2[11]+0.3535533905932737*magB2[2]); 
  gkyl_mat_set(A,16,0,0.3535533905932737*magB2[16]); 
  gkyl_mat_set(A,16,1,0.3535533905932737*magB2[19]); 
  gkyl_mat_set(A,16,2,0.3535533905932737*magB2[9]); 
  gkyl_mat_set(A,16,3,0.3162277660168379*magB2[6]); 
  gkyl_mat_set(A,16,4,0.3535533905932737*magB2[15]); 
  gkyl_mat_set(A,16,5,0.3162277660168379*magB2[10]); 
  gkyl_mat_set(A,16,6,0.2828427124746191*magB2[14]+0.3162277660168379*magB2[3]); 
  gkyl_mat_set(A,16,7,0.0); 
  gkyl_mat_set(A,16,8,0.3162277660168379*magB2[16]); 
  gkyl_mat_set(A,16,9,0.2258769757263128*magB2[16]+0.3535533905932737*magB2[2]); 
  gkyl_mat_set(A,16,10,0.282842712474619*magB2[18]+0.3162277660168379*magB2[5]); 
  gkyl_mat_set(A,16,11,0.0); 
  gkyl_mat_set(A,16,12,0.3162277660168379*magB2[19]); 
  gkyl_mat_set(A,16,13,0.3162277660168379*magB2[17]); 
  gkyl_mat_set(A,16,14,0.2828427124746191*magB2[6]); 
  gkyl_mat_set(A,16,15,0.2258769757263128*magB2[19]+0.3535533905932737*magB2[4]); 
  gkyl_mat_set(A,16,16,0.2258769757263128*magB2[9]+0.3162277660168379*magB2[8]+0.3535533905932737*magB2[0]); 
  gkyl_mat_set(A,16,17,0.3162277660168379*magB2[13]); 
  gkyl_mat_set(A,16,18,0.282842712474619*magB2[10]); 
  gkyl_mat_set(A,16,19,0.2258769757263128*magB2[15]+0.3162277660168379*magB2[12]+0.3535533905932737*magB2[1]); 
  gkyl_mat_set(A,17,0,0.3535533905932737*magB2[17]); 
  gkyl_mat_set(A,17,1,0.3162277660168379*magB2[10]); 
  gkyl_mat_set(A,17,2,0.3535533905932737*magB2[13]); 
  gkyl_mat_set(A,17,3,0.3535533905932737*magB2[11]); 
  gkyl_mat_set(A,17,4,0.2828427124746191*magB2[18]+0.3162277660168379*magB2[5]); 
  gkyl_mat_set(A,17,5,0.2828427124746191*magB2[19]+0.3162277660168379*magB2[4]); 
  gkyl_mat_set(A,17,6,0.3535533905932737*magB2[7]); 
  gkyl_mat_set(A,17,7,0.2258769757263128*magB2[17]+0.3535533905932737*magB2[6]); 
  gkyl_mat_set(A,17,8,0.3162277660168379*magB2[17]); 
  gkyl_mat_set(A,17,9,0.3162277660168379*magB2[17]); 
  gkyl_mat_set(A,17,10,0.282842712474619*magB2[15]+0.282842712474619*magB2[12]+0.3162277660168379*magB2[1]); 
  gkyl_mat_set(A,17,11,0.3162277660168379*magB2[14]+0.2258769757263128*magB2[13]+0.3535533905932737*magB2[3]); 
  gkyl_mat_set(A,17,12,0.282842712474619*magB2[10]); 
  gkyl_mat_set(A,17,13,0.3162277660168379*magB2[16]+0.2258769757263128*magB2[11]+0.3535533905932737*magB2[2]); 
  gkyl_mat_set(A,17,14,0.3162277660168379*magB2[11]); 
  gkyl_mat_set(A,17,15,0.282842712474619*magB2[10]); 
  gkyl_mat_set(A,17,16,0.3162277660168379*magB2[13]); 
  gkyl_mat_set(A,17,17,0.3162277660168379*magB2[9]+0.3162277660168379*magB2[8]+0.2258769757263128*magB2[7]+0.3535533905932737*magB2[0]); 
  gkyl_mat_set(A,17,18,0.2529822128134704*magB2[19]+0.2828427124746191*magB2[4]); 
  gkyl_mat_set(A,17,19,0.2529822128134704*magB2[18]+0.2828427124746191*magB2[5]); 
  gkyl_mat_set(A,18,0,0.3535533905932737*magB2[18]); 
  gkyl_mat_set(A,18,1,0.3535533905932737*magB2[14]); 
  gkyl_mat_set(A,18,2,0.3162277660168379*magB2[10]); 
  gkyl_mat_set(A,18,3,0.3535533905932737*magB2[12]); 
  gkyl_mat_set(A,18,4,0.2828427124746191*magB2[17]+0.3162277660168379*magB2[6]); 
  gkyl_mat_set(A,18,5,0.3535533905932737*magB2[8]); 
  gkyl_mat_set(A,18,6,0.2828427124746191*magB2[19]+0.3162277660168379*magB2[4]); 
  gkyl_mat_set(A,18,7,0.3162277660168379*magB2[18]); 
  gkyl_mat_set(A,18,8,0.2258769757263128*magB2[18]+0.3535533905932737*magB2[5]); 
  gkyl_mat_set(A,18,9,0.3162277660168379*magB2[18]); 
  gkyl_mat_set(A,18,10,0.282842712474619*magB2[16]+0.282842712474619*magB2[11]+0.3162277660168379*magB2[2]); 
  gkyl_mat_set(A,18,11,0.282842712474619*magB2[10]); 
  gkyl_mat_set(A,18,12,0.2258769757263128*magB2[14]+0.3162277660168379*magB2[13]+0.3535533905932737*magB2[3]); 
  gkyl_mat_set(A,18,13,0.3162277660168379*magB2[12]); 
  gkyl_mat_set(A,18,14,0.3162277660168379*magB2[15]+0.2258769757263128*magB2[12]+0.3535533905932737*magB2[1]); 
  gkyl_mat_set(A,18,15,0.3162277660168379*magB2[14]); 
  gkyl_mat_set(A,18,16,0.282842712474619*magB2[10]); 
  gkyl_mat_set(A,18,17,0.2529822128134704*magB2[19]+0.2828427124746191*magB2[4]); 
  gkyl_mat_set(A,18,18,0.3162277660168379*magB2[9]+0.2258769757263128*magB2[8]+0.3162277660168379*magB2[7]+0.3535533905932737*magB2[0]); 
  gkyl_mat_set(A,18,19,0.2529822128134704*magB2[17]+0.2828427124746191*magB2[6]); 
  gkyl_mat_set(A,19,0,0.3535533905932737*magB2[19]); 
  gkyl_mat_set(A,19,1,0.3535533905932737*magB2[16]); 
  gkyl_mat_set(A,19,2,0.3535533905932737*magB2[15]); 
  gkyl_mat_set(A,19,3,0.3162277660168379*magB2[10]); 
  gkyl_mat_set(A,19,4,0.3535533905932737*magB2[9]); 
  gkyl_mat_set(A,19,5,0.2828427124746191*magB2[17]+0.3162277660168379*magB2[6]); 
  gkyl_mat_set(A,19,6,0.2828427124746191*magB2[18]+0.3162277660168379*magB2[5]); 
  gkyl_mat_set(A,19,7,0.3162277660168379*magB2[19]); 
  gkyl_mat_set(A,19,8,0.3162277660168379*magB2[19]); 
  gkyl_mat_set(A,19,9,0.2258769757263128*magB2[19]+0.3535533905932737*magB2[4]); 
  gkyl_mat_set(A,19,10,0.282842712474619*magB2[14]+0.282842712474619*magB2[13]+0.3162277660168379*magB2[3]); 
  gkyl_mat_set(A,19,11,0.3162277660168379*magB2[15]); 
  gkyl_mat_set(A,19,12,0.3162277660168379*magB2[16]); 
  gkyl_mat_set(A,19,13,0.282842712474619*magB2[10]); 
  gkyl_mat_set(A,19,14,0.282842712474619*magB2[10]); 
  gkyl_mat_set(A,19,15,0.2258769757263128*magB2[16]+0.3162277660168379*magB2[11]+0.3535533905932737*magB2[2]); 
  gkyl_mat_set(A,19,16,0.2258769757263128*magB2[15]+0.3162277660168379*magB2[12]+0.3535533905932737*magB2[1]); 
  gkyl_mat_set(A,19,17,0.2529822128134704*magB2[18]+0.2828427124746191*magB2[5]); 
  gkyl_mat_set(A,19,18,0.2529822128134704*magB2[17]+0.2828427124746191*magB2[6]); 
  gkyl_mat_set(A,19,19,0.2258769757263128*magB2[9]+0.3162277660168379*magB2[8]+0.3162277660168379*magB2[7]+0.3535533905932737*magB2[0]); 
  return cell_avg;
} 
