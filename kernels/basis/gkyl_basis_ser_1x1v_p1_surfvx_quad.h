static inline double 
ser_1x1v_p1_surfvx_quad_0(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return (-0.8660254037844386*f[3])+0.8660254037844386*f[2]-0.5*f[1]+0.5*f[0]; 
  else 
    return 0.8660254037844386*f[3]-0.8660254037844386*f[2]-0.5*f[1]+0.5*f[0]; 
} 
static inline double 
ser_1x1v_p1_surfvx_quad_1(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return 0.8660254037844386*(f[3]+f[2])+0.5*(f[1]+f[0]); 
  else 
    return 0.5*(f[1]+f[0])-0.8660254037844386*(f[3]+f[2]); 
} 
