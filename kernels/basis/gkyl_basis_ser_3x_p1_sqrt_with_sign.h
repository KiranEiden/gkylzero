GKYL_CU_DH static inline void 
ser_3x_p1_sqrt_with_sign(const double *ASign, const double *A, double *ASqrt) 
{ 
  // ASign: Input DG field, used to get correct sign of Asqrt. 
  // A:     Input DG field. 
  // ASqrt: Output DG field (expansion of sqrt(A), with sign determined by Asign). 
 
  double AOrd[8] = {0.0}; 

  double temp = 0.0; 
  double temp_sign = 0.0; 
  temp = (-0.3535533905932737*A[7])+0.3535533905932737*(A[6]+A[5]+A[4])-0.3535533905932737*(A[3]+A[2]+A[1])+0.3535533905932737*A[0]; 
  temp_sign = (-0.3535533905932737*ASign[7])+0.3535533905932737*(ASign[6]+ASign[5]+ASign[4])-0.3535533905932737*(ASign[3]+ASign[2]+ASign[1])+0.3535533905932737*ASign[0]; 
  if (temp < 0.0) { 
  AOrd[0] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[0] = -sqrt(temp); 
  } else { 
  AOrd[0] = sqrt(temp); 
  } 
  temp = 0.3535533905932737*A[7]-0.3535533905932737*(A[6]+A[5])+0.3535533905932737*(A[4]+A[3])-0.3535533905932737*(A[2]+A[1])+0.3535533905932737*A[0]; 
  temp_sign = 0.3535533905932737*ASign[7]-0.3535533905932737*(ASign[6]+ASign[5])+0.3535533905932737*(ASign[4]+ASign[3])-0.3535533905932737*(ASign[2]+ASign[1])+0.3535533905932737*ASign[0]; 
  if (temp < 0.0) { 
  AOrd[1] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[1] = -sqrt(temp); 
  } else { 
  AOrd[1] = sqrt(temp); 
  } 
  temp = 0.3535533905932737*A[7]-0.3535533905932737*A[6]+0.3535533905932737*A[5]-0.3535533905932737*(A[4]+A[3])+0.3535533905932737*A[2]-0.3535533905932737*A[1]+0.3535533905932737*A[0]; 
  temp_sign = 0.3535533905932737*ASign[7]-0.3535533905932737*ASign[6]+0.3535533905932737*ASign[5]-0.3535533905932737*(ASign[4]+ASign[3])+0.3535533905932737*ASign[2]-0.3535533905932737*ASign[1]+0.3535533905932737*ASign[0]; 
  if (temp < 0.0) { 
  AOrd[2] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[2] = -sqrt(temp); 
  } else { 
  AOrd[2] = sqrt(temp); 
  } 
  temp = (-0.3535533905932737*A[7])+0.3535533905932737*A[6]-0.3535533905932737*(A[5]+A[4])+0.3535533905932737*(A[3]+A[2])-0.3535533905932737*A[1]+0.3535533905932737*A[0]; 
  temp_sign = (-0.3535533905932737*ASign[7])+0.3535533905932737*ASign[6]-0.3535533905932737*(ASign[5]+ASign[4])+0.3535533905932737*(ASign[3]+ASign[2])-0.3535533905932737*ASign[1]+0.3535533905932737*ASign[0]; 
  if (temp < 0.0) { 
  AOrd[3] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[3] = -sqrt(temp); 
  } else { 
  AOrd[3] = sqrt(temp); 
  } 
  temp = 0.3535533905932737*(A[7]+A[6])-0.3535533905932737*(A[5]+A[4]+A[3]+A[2])+0.3535533905932737*(A[1]+A[0]); 
  temp_sign = 0.3535533905932737*(ASign[7]+ASign[6])-0.3535533905932737*(ASign[5]+ASign[4]+ASign[3]+ASign[2])+0.3535533905932737*(ASign[1]+ASign[0]); 
  if (temp < 0.0) { 
  AOrd[4] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[4] = -sqrt(temp); 
  } else { 
  AOrd[4] = sqrt(temp); 
  } 
  temp = (-0.3535533905932737*(A[7]+A[6]))+0.3535533905932737*A[5]-0.3535533905932737*A[4]+0.3535533905932737*A[3]-0.3535533905932737*A[2]+0.3535533905932737*(A[1]+A[0]); 
  temp_sign = (-0.3535533905932737*(ASign[7]+ASign[6]))+0.3535533905932737*ASign[5]-0.3535533905932737*ASign[4]+0.3535533905932737*ASign[3]-0.3535533905932737*ASign[2]+0.3535533905932737*(ASign[1]+ASign[0]); 
  if (temp < 0.0) { 
  AOrd[5] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[5] = -sqrt(temp); 
  } else { 
  AOrd[5] = sqrt(temp); 
  } 
  temp = (-0.3535533905932737*(A[7]+A[6]+A[5]))+0.3535533905932737*A[4]-0.3535533905932737*A[3]+0.3535533905932737*(A[2]+A[1]+A[0]); 
  temp_sign = (-0.3535533905932737*(ASign[7]+ASign[6]+ASign[5]))+0.3535533905932737*ASign[4]-0.3535533905932737*ASign[3]+0.3535533905932737*(ASign[2]+ASign[1]+ASign[0]); 
  if (temp < 0.0) { 
  AOrd[6] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[6] = -sqrt(temp); 
  } else { 
  AOrd[6] = sqrt(temp); 
  } 
  temp = 0.3535533905932737*(A[7]+A[6]+A[5]+A[4]+A[3]+A[2]+A[1]+A[0]); 
  temp_sign = 0.3535533905932737*(ASign[7]+ASign[6]+ASign[5]+ASign[4]+ASign[3]+ASign[2]+ASign[1]+ASign[0]); 
  if (temp < 0.0) { 
  AOrd[7] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[7] = -sqrt(temp); 
  } else { 
  AOrd[7] = sqrt(temp); 
  } 
  ASqrt[0] = 0.3535533905932737*AOrd[7]+0.3535533905932737*AOrd[6]+0.3535533905932737*AOrd[5]+0.3535533905932737*AOrd[4]+0.3535533905932737*AOrd[3]+0.3535533905932737*AOrd[2]+0.3535533905932737*AOrd[1]+0.3535533905932737*AOrd[0]; 
  ASqrt[1] = 0.3535533905932737*AOrd[7]+0.3535533905932737*AOrd[6]+0.3535533905932737*AOrd[5]+0.3535533905932737*AOrd[4]-0.3535533905932737*AOrd[3]-0.3535533905932737*AOrd[2]-0.3535533905932737*AOrd[1]-0.3535533905932737*AOrd[0]; 
  ASqrt[2] = 0.3535533905932737*AOrd[7]+0.3535533905932737*AOrd[6]-0.3535533905932737*AOrd[5]-0.3535533905932737*AOrd[4]+0.3535533905932737*AOrd[3]+0.3535533905932737*AOrd[2]-0.3535533905932737*AOrd[1]-0.3535533905932737*AOrd[0]; 
  ASqrt[3] = 0.3535533905932737*AOrd[7]-0.3535533905932737*AOrd[6]+0.3535533905932737*AOrd[5]-0.3535533905932737*AOrd[4]+0.3535533905932737*AOrd[3]-0.3535533905932737*AOrd[2]+0.3535533905932737*AOrd[1]-0.3535533905932737*AOrd[0]; 
  ASqrt[4] = 0.3535533905932737*AOrd[7]+0.3535533905932737*AOrd[6]-0.3535533905932737*AOrd[5]-0.3535533905932737*AOrd[4]-0.3535533905932737*AOrd[3]-0.3535533905932737*AOrd[2]+0.3535533905932737*AOrd[1]+0.3535533905932737*AOrd[0]; 
  ASqrt[5] = 0.3535533905932737*AOrd[7]-0.3535533905932737*AOrd[6]+0.3535533905932737*AOrd[5]-0.3535533905932737*AOrd[4]-0.3535533905932737*AOrd[3]+0.3535533905932737*AOrd[2]-0.3535533905932737*AOrd[1]+0.3535533905932737*AOrd[0]; 
  ASqrt[6] = 0.3535533905932737*AOrd[7]-0.3535533905932737*AOrd[6]-0.3535533905932737*AOrd[5]+0.3535533905932737*AOrd[4]+0.3535533905932737*AOrd[3]-0.3535533905932737*AOrd[2]-0.3535533905932737*AOrd[1]+0.3535533905932737*AOrd[0]; 
  ASqrt[7] = 0.3535533905932737*AOrd[7]-0.3535533905932737*AOrd[6]-0.3535533905932737*AOrd[5]+0.3535533905932737*AOrd[4]-0.3535533905932737*AOrd[3]+0.3535533905932737*AOrd[2]+0.3535533905932737*AOrd[1]-0.3535533905932737*AOrd[0]; 

} 
 
