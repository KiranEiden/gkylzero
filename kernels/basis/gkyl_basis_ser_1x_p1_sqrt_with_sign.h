GKYL_CU_DH static inline void 
ser_1x_p1_sqrt_with_sign(const double *ASign, const double *A, double *ASqrt) 
{ 
  // ASign: Input DG field, used to get correct sign of Asqrt. 
  // A:     Input DG field. 
  // ASqrt: Output DG field (expansion of sqrt(A), with sign determined by Asign). 
 
  double AOrd[3] = {0.0}; 

  double temp = 0.0; 
  double temp_sign = 0.0; 
  temp = 0.7071067811865475*A[0]-0.9486832980505137*A[1]; 
  temp_sign = 0.7071067811865475*ASign[0]-0.9486832980505137*ASign[1]; 
  if (temp_sign < 0.0) { 
  AOrd[0] = -sqrt(temp); 
  } else { 
  AOrd[0] = sqrt(temp); 
  } 
  temp = 0.7071067811865475*A[0]; 
  temp_sign = 0.7071067811865475*ASign[0]; 
  if (temp_sign < 0.0) { 
  AOrd[1] = -sqrt(temp); 
  } else { 
  AOrd[1] = sqrt(temp); 
  } 
  temp = 0.9486832980505137*A[1]+0.7071067811865475*A[0]; 
  temp_sign = 0.9486832980505137*ASign[1]+0.7071067811865475*ASign[0]; 
  if (temp_sign < 0.0) { 
  AOrd[2] = -sqrt(temp); 
  } else { 
  AOrd[2] = sqrt(temp); 
  } 
  ASqrt[0] = 0.392837100659193*AOrd[2]+0.6285393610547091*AOrd[1]+0.392837100659193*AOrd[0]; 
  ASqrt[1] = 0.5270462766947298*AOrd[2]-0.5270462766947298*AOrd[0]; 

} 
 
