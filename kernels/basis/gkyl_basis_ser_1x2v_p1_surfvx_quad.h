static inline double 
ser_1x2v_p1_surfvx_quad_0(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return 0.6123724356957944*f[7]-0.6123724356957944*f[6]+0.3535533905932737*f[5]-0.6123724356957944*f[4]-0.3535533905932737*f[3]+0.6123724356957944*f[2]-0.3535533905932737*f[1]+0.3535533905932737*f[0]; 
  else 
    return (-0.6123724356957944*f[7])+0.6123724356957944*f[6]+0.3535533905932737*f[5]+0.6123724356957944*f[4]-0.3535533905932737*f[3]-0.6123724356957944*f[2]-0.3535533905932737*f[1]+0.3535533905932737*f[0]; 
} 
static inline double 
ser_1x2v_p1_surfvx_quad_1(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return (-0.6123724356957944*(f[7]+f[6]))-0.3535533905932737*f[5]+0.6123724356957944*f[4]-0.3535533905932737*f[3]+0.6123724356957944*f[2]+0.3535533905932737*(f[1]+f[0]); 
  else 
    return 0.6123724356957944*(f[7]+f[6])-0.3535533905932737*f[5]-0.6123724356957944*f[4]-0.3535533905932737*f[3]-0.6123724356957944*f[2]+0.3535533905932737*(f[1]+f[0]); 
} 
static inline double 
ser_1x2v_p1_surfvx_quad_2(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return (-0.6123724356957944*f[7])+0.6123724356957944*f[6]-0.3535533905932737*f[5]-0.6123724356957944*f[4]+0.3535533905932737*f[3]+0.6123724356957944*f[2]-0.3535533905932737*f[1]+0.3535533905932737*f[0]; 
  else 
    return 0.6123724356957944*f[7]-0.6123724356957944*f[6]-0.3535533905932737*f[5]+0.6123724356957944*f[4]+0.3535533905932737*f[3]-0.6123724356957944*f[2]-0.3535533905932737*f[1]+0.3535533905932737*f[0]; 
} 
static inline double 
ser_1x2v_p1_surfvx_quad_3(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return 0.6123724356957944*(f[7]+f[6])+0.3535533905932737*f[5]+0.6123724356957944*f[4]+0.3535533905932737*f[3]+0.6123724356957944*f[2]+0.3535533905932737*(f[1]+f[0]); 
  else 
    return (-0.6123724356957944*(f[7]+f[6]))+0.3535533905932737*f[5]-0.6123724356957944*f[4]+0.3535533905932737*f[3]-0.6123724356957944*f[2]+0.3535533905932737*(f[1]+f[0]); 
} 
