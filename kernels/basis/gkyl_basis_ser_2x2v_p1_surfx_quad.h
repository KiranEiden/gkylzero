GKYL_CU_DH static inline double 
ser_2x2v_p1_surfx_quad_0(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return (-0.4330127018922193*f[15])-0.25*f[14]+0.4330127018922193*(f[13]+f[12]+f[11])+0.25*(f[10]+f[9])-0.4330127018922193*f[8]+0.25*f[7]-0.4330127018922193*(f[6]+f[5])-0.25*(f[4]+f[3]+f[2])+0.4330127018922193*f[1]+0.25*f[0]; 
  else 
    return 0.4330127018922193*f[15]-0.25*f[14]-0.4330127018922193*(f[13]+f[12]+f[11])+0.25*(f[10]+f[9])+0.4330127018922193*f[8]+0.25*f[7]+0.4330127018922193*(f[6]+f[5])-0.25*(f[4]+f[3]+f[2])-0.4330127018922193*f[1]+0.25*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_2x2v_p1_surfx_quad_1(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return 0.4330127018922193*f[15]+0.25*f[14]+0.4330127018922193*f[13]-0.4330127018922193*(f[12]+f[11])+0.25*f[10]-0.25*f[9]-0.4330127018922193*f[8]-0.25*f[7]-0.4330127018922193*f[6]+0.4330127018922193*f[5]-0.25*(f[4]+f[3])+0.25*f[2]+0.4330127018922193*f[1]+0.25*f[0]; 
  else 
    return (-0.4330127018922193*f[15])+0.25*f[14]-0.4330127018922193*f[13]+0.4330127018922193*(f[12]+f[11])+0.25*f[10]-0.25*f[9]+0.4330127018922193*f[8]-0.25*f[7]+0.4330127018922193*f[6]-0.4330127018922193*f[5]-0.25*(f[4]+f[3])+0.25*f[2]-0.4330127018922193*f[1]+0.25*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_2x2v_p1_surfx_quad_2(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return 0.4330127018922193*f[15]+0.25*f[14]-0.4330127018922193*f[13]+0.4330127018922193*f[12]-0.4330127018922193*f[11]-0.25*f[10]+0.25*f[9]-0.4330127018922193*f[8]-0.25*f[7]+0.4330127018922193*f[6]-0.4330127018922193*f[5]-0.25*f[4]+0.25*f[3]-0.25*f[2]+0.4330127018922193*f[1]+0.25*f[0]; 
  else 
    return (-0.4330127018922193*f[15])+0.25*f[14]+0.4330127018922193*f[13]-0.4330127018922193*f[12]+0.4330127018922193*f[11]-0.25*f[10]+0.25*f[9]+0.4330127018922193*f[8]-0.25*f[7]-0.4330127018922193*f[6]+0.4330127018922193*f[5]-0.25*f[4]+0.25*f[3]-0.25*f[2]-0.4330127018922193*f[1]+0.25*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_2x2v_p1_surfx_quad_3(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return (-0.4330127018922193*f[15])-0.25*f[14]-0.4330127018922193*(f[13]+f[12])+0.4330127018922193*f[11]-0.25*(f[10]+f[9])-0.4330127018922193*f[8]+0.25*f[7]+0.4330127018922193*(f[6]+f[5])-0.25*f[4]+0.25*(f[3]+f[2])+0.4330127018922193*f[1]+0.25*f[0]; 
  else 
    return 0.4330127018922193*f[15]-0.25*f[14]+0.4330127018922193*(f[13]+f[12])-0.4330127018922193*f[11]-0.25*(f[10]+f[9])+0.4330127018922193*f[8]+0.25*f[7]-0.4330127018922193*(f[6]+f[5])-0.25*f[4]+0.25*(f[3]+f[2])-0.4330127018922193*f[1]+0.25*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_2x2v_p1_surfx_quad_4(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return 0.4330127018922193*f[15]+0.25*f[14]-0.4330127018922193*(f[13]+f[12])+0.4330127018922193*f[11]-0.25*(f[10]+f[9])+0.4330127018922193*f[8]+0.25*f[7]-0.4330127018922193*(f[6]+f[5])+0.25*f[4]-0.25*(f[3]+f[2])+0.4330127018922193*f[1]+0.25*f[0]; 
  else 
    return (-0.4330127018922193*f[15])+0.25*f[14]+0.4330127018922193*(f[13]+f[12])-0.4330127018922193*f[11]-0.25*(f[10]+f[9])-0.4330127018922193*f[8]+0.25*f[7]+0.4330127018922193*(f[6]+f[5])+0.25*f[4]-0.25*(f[3]+f[2])-0.4330127018922193*f[1]+0.25*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_2x2v_p1_surfx_quad_5(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return (-0.4330127018922193*f[15])-0.25*f[14]-0.4330127018922193*f[13]+0.4330127018922193*f[12]-0.4330127018922193*f[11]-0.25*f[10]+0.25*f[9]+0.4330127018922193*f[8]-0.25*f[7]-0.4330127018922193*f[6]+0.4330127018922193*f[5]+0.25*f[4]-0.25*f[3]+0.25*f[2]+0.4330127018922193*f[1]+0.25*f[0]; 
  else 
    return 0.4330127018922193*f[15]-0.25*f[14]+0.4330127018922193*f[13]-0.4330127018922193*f[12]+0.4330127018922193*f[11]-0.25*f[10]+0.25*f[9]-0.4330127018922193*f[8]-0.25*f[7]+0.4330127018922193*f[6]-0.4330127018922193*f[5]+0.25*f[4]-0.25*f[3]+0.25*f[2]-0.4330127018922193*f[1]+0.25*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_2x2v_p1_surfx_quad_6(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return (-0.4330127018922193*f[15])-0.25*f[14]+0.4330127018922193*f[13]-0.4330127018922193*(f[12]+f[11])+0.25*f[10]-0.25*f[9]+0.4330127018922193*f[8]-0.25*f[7]+0.4330127018922193*f[6]-0.4330127018922193*f[5]+0.25*(f[4]+f[3])-0.25*f[2]+0.4330127018922193*f[1]+0.25*f[0]; 
  else 
    return 0.4330127018922193*f[15]-0.25*f[14]-0.4330127018922193*f[13]+0.4330127018922193*(f[12]+f[11])+0.25*f[10]-0.25*f[9]-0.4330127018922193*f[8]-0.25*f[7]-0.4330127018922193*f[6]+0.4330127018922193*f[5]+0.25*(f[4]+f[3])-0.25*f[2]-0.4330127018922193*f[1]+0.25*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_2x2v_p1_surfx_quad_7(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return 0.4330127018922193*f[15]+0.25*f[14]+0.4330127018922193*(f[13]+f[12]+f[11])+0.25*(f[10]+f[9])+0.4330127018922193*f[8]+0.25*f[7]+0.4330127018922193*(f[6]+f[5])+0.25*(f[4]+f[3]+f[2])+0.4330127018922193*f[1]+0.25*f[0]; 
  else 
    return (-0.4330127018922193*f[15])+0.25*f[14]-0.4330127018922193*(f[13]+f[12]+f[11])+0.25*(f[10]+f[9])-0.4330127018922193*f[8]+0.25*f[7]-0.4330127018922193*(f[6]+f[5])+0.25*(f[4]+f[3]+f[2])-0.4330127018922193*f[1]+0.25*f[0]; 
} 
