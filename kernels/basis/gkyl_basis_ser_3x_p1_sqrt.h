GKYL_CU_DH static inline void 
ser_3x_p1_sqrt(const double *A, double *ASqrt) 
{ 
  // A:     Input DG field. 
  // ASqrt: Output DG field (expansion of sqrt(A)). 
 
  double AOrd[8] = {0.0}; 

  AOrd[0] = sqrt((-0.3535533905932737*A[7])+0.3535533905932737*(A[6]+A[5]+A[4])-0.3535533905932737*(A[3]+A[2]+A[1])+0.3535533905932737*A[0]); 
  AOrd[1] = sqrt(0.3535533905932737*A[7]-0.3535533905932737*(A[6]+A[5])+0.3535533905932737*(A[4]+A[3])-0.3535533905932737*(A[2]+A[1])+0.3535533905932737*A[0]); 
  AOrd[2] = sqrt(0.3535533905932737*A[7]-0.3535533905932737*A[6]+0.3535533905932737*A[5]-0.3535533905932737*(A[4]+A[3])+0.3535533905932737*A[2]-0.3535533905932737*A[1]+0.3535533905932737*A[0]); 
  AOrd[3] = sqrt((-0.3535533905932737*A[7])+0.3535533905932737*A[6]-0.3535533905932737*(A[5]+A[4])+0.3535533905932737*(A[3]+A[2])-0.3535533905932737*A[1]+0.3535533905932737*A[0]); 
  AOrd[4] = sqrt(0.3535533905932737*(A[7]+A[6])-0.3535533905932737*(A[5]+A[4]+A[3]+A[2])+0.3535533905932737*(A[1]+A[0])); 
  AOrd[5] = sqrt((-0.3535533905932737*(A[7]+A[6]))+0.3535533905932737*A[5]-0.3535533905932737*A[4]+0.3535533905932737*A[3]-0.3535533905932737*A[2]+0.3535533905932737*(A[1]+A[0])); 
  AOrd[6] = sqrt((-0.3535533905932737*(A[7]+A[6]+A[5]))+0.3535533905932737*A[4]-0.3535533905932737*A[3]+0.3535533905932737*(A[2]+A[1]+A[0])); 
  AOrd[7] = sqrt(0.3535533905932737*(A[7]+A[6]+A[5]+A[4]+A[3]+A[2]+A[1]+A[0])); 
  ASqrt[0] = 0.3535533905932737*AOrd[7]+0.3535533905932737*AOrd[6]+0.3535533905932737*AOrd[5]+0.3535533905932737*AOrd[4]+0.3535533905932737*AOrd[3]+0.3535533905932737*AOrd[2]+0.3535533905932737*AOrd[1]+0.3535533905932737*AOrd[0]; 
  ASqrt[1] = 0.3535533905932737*AOrd[7]+0.3535533905932737*AOrd[6]+0.3535533905932737*AOrd[5]+0.3535533905932737*AOrd[4]-0.3535533905932737*AOrd[3]-0.3535533905932737*AOrd[2]-0.3535533905932737*AOrd[1]-0.3535533905932737*AOrd[0]; 
  ASqrt[2] = 0.3535533905932737*AOrd[7]+0.3535533905932737*AOrd[6]-0.3535533905932737*AOrd[5]-0.3535533905932737*AOrd[4]+0.3535533905932737*AOrd[3]+0.3535533905932737*AOrd[2]-0.3535533905932737*AOrd[1]-0.3535533905932737*AOrd[0]; 
  ASqrt[3] = 0.3535533905932737*AOrd[7]-0.3535533905932737*AOrd[6]+0.3535533905932737*AOrd[5]-0.3535533905932737*AOrd[4]+0.3535533905932737*AOrd[3]-0.3535533905932737*AOrd[2]+0.3535533905932737*AOrd[1]-0.3535533905932737*AOrd[0]; 
  ASqrt[4] = 0.3535533905932737*AOrd[7]+0.3535533905932737*AOrd[6]-0.3535533905932737*AOrd[5]-0.3535533905932737*AOrd[4]-0.3535533905932737*AOrd[3]-0.3535533905932737*AOrd[2]+0.3535533905932737*AOrd[1]+0.3535533905932737*AOrd[0]; 
  ASqrt[5] = 0.3535533905932737*AOrd[7]-0.3535533905932737*AOrd[6]+0.3535533905932737*AOrd[5]-0.3535533905932737*AOrd[4]-0.3535533905932737*AOrd[3]+0.3535533905932737*AOrd[2]-0.3535533905932737*AOrd[1]+0.3535533905932737*AOrd[0]; 
  ASqrt[6] = 0.3535533905932737*AOrd[7]-0.3535533905932737*AOrd[6]-0.3535533905932737*AOrd[5]+0.3535533905932737*AOrd[4]+0.3535533905932737*AOrd[3]-0.3535533905932737*AOrd[2]-0.3535533905932737*AOrd[1]+0.3535533905932737*AOrd[0]; 
  ASqrt[7] = 0.3535533905932737*AOrd[7]-0.3535533905932737*AOrd[6]-0.3535533905932737*AOrd[5]+0.3535533905932737*AOrd[4]-0.3535533905932737*AOrd[3]+0.3535533905932737*AOrd[2]+0.3535533905932737*AOrd[1]-0.3535533905932737*AOrd[0]; 

} 
 