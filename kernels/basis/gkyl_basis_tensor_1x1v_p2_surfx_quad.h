GKYL_CU_DH static inline double 
tensor_1x1v_p2_surfx_quad_0(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return f[8]+0.7745966692414833*f[7]-1.5*f[6]+0.4472135954999579*f[5]+1.118033988749895*f[4]-1.161895003862225*f[3]-0.6708203932499369*f[2]+0.8660254037844386*f[1]+0.5*f[0]; 
  else 
    return f[8]-0.7745966692414833*f[7]-1.5*f[6]+0.4472135954999579*f[5]+1.118033988749895*f[4]+1.161895003862225*f[3]-0.6708203932499369*f[2]-0.8660254037844386*f[1]+0.5*f[0]; 
} 
GKYL_CU_DH static inline double 
tensor_1x1v_p2_surfx_quad_1(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return (-1.25*f[8])-0.9682458365518543*f[7]-0.5590169943749475*f[5]+1.118033988749895*f[4]+0.8660254037844386*f[1]+0.5*f[0]; 
  else 
    return (-1.25*f[8])+0.9682458365518543*f[7]-0.5590169943749475*f[5]+1.118033988749895*f[4]-0.8660254037844386*f[1]+0.5*f[0]; 
} 
GKYL_CU_DH static inline double 
tensor_1x1v_p2_surfx_quad_2(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return f[8]+0.7745966692414833*f[7]+1.5*f[6]+0.4472135954999579*f[5]+1.118033988749895*f[4]+1.161895003862225*f[3]+0.6708203932499369*f[2]+0.8660254037844386*f[1]+0.5*f[0]; 
  else 
    return f[8]-0.7745966692414833*f[7]+1.5*f[6]+0.4472135954999579*f[5]+1.118033988749895*f[4]-1.161895003862225*f[3]+0.6708203932499369*f[2]-0.8660254037844386*f[1]+0.5*f[0]; 
} 
