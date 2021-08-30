GKYL_CU_DH static inline double 
ser_1x3v_p2_surfvz_quad_0(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return (-1.349999999999999*f[47])+0.697137002317335*(f[46]+f[45]+f[44])+1.006230589874905*(f[43]+f[42]+f[41])-0.5196152422706628*(f[40]+f[39])-0.5196152422706631*(f[38]+f[37])-0.5196152422706628*f[36]-0.5196152422706631*f[35]+0.4024922359499621*(f[34]+f[33]+f[32])-1.045705503476002*f[31]-0.75*(f[30]+f[29]+f[28])+0.3872983346207416*(f[27]+f[26]+f[25])-0.2999999999999998*(f[24]+f[23])-0.2999999999999999*(f[22]+f[21])-0.2999999999999998*f[20]-0.2999999999999999*f[19]+0.7794228634059945*(f[18]+f[17]+f[16])-0.6037383539249431*f[15]+0.5590169943749475*f[14]+0.223606797749979*(f[13]+f[12]+f[11])-0.5809475019311124*(f[10]+f[9]+f[8])+0.45*(f[7]+f[6]+f[5])+0.4330127018922193*f[4]-0.3354101966249685*(f[3]+f[2]+f[1])+0.25*f[0]; 
  else 
    return (-1.349999999999999*f[47])-0.697137002317335*(f[46]+f[45]+f[44])+1.006230589874905*(f[43]+f[42]+f[41])+0.5196152422706628*(f[40]+f[39])+0.5196152422706631*(f[38]+f[37])+0.5196152422706628*f[36]+0.5196152422706631*f[35]+0.4024922359499621*(f[34]+f[33]+f[32])+1.045705503476002*f[31]-0.75*(f[30]+f[29]+f[28])-0.3872983346207416*(f[27]+f[26]+f[25])-0.2999999999999998*(f[24]+f[23])-0.2999999999999999*(f[22]+f[21])-0.2999999999999998*f[20]-0.2999999999999999*f[19]-0.7794228634059945*(f[18]+f[17]+f[16])-0.6037383539249431*f[15]+0.5590169943749475*f[14]+0.223606797749979*(f[13]+f[12]+f[11])+0.5809475019311124*(f[10]+f[9]+f[8])+0.45*(f[7]+f[6]+f[5])-0.4330127018922193*f[4]-0.3354101966249685*(f[3]+f[2]+f[1])+0.25*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_1x3v_p2_surfvz_quad_1(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return (-0.8714212528966688*f[44])+1.006230589874905*f[43]-0.5196152422706628*f[40]-0.5196152422706631*f[38]+0.6495190528383289*(f[37]+f[35])-0.5031152949374527*f[32]-0.75*(f[30]+f[29])+0.3872983346207416*(f[27]+f[26])-0.4841229182759271*f[25]-0.2999999999999998*f[24]-0.2999999999999999*f[22]+0.375*(f[21]+f[19])+0.7794228634059945*f[18]+0.5590169943749475*f[14]+0.223606797749979*(f[13]+f[12])-0.2795084971874737*f[11]-0.5809475019311124*(f[10]+f[9])+0.45*f[7]+0.4330127018922193*f[4]-0.3354101966249685*(f[3]+f[2])+0.25*f[0]; 
  else 
    return 0.8714212528966688*f[44]+1.006230589874905*f[43]+0.5196152422706628*f[40]+0.5196152422706631*f[38]-0.6495190528383289*(f[37]+f[35])-0.5031152949374527*f[32]-0.75*(f[30]+f[29])-0.3872983346207416*(f[27]+f[26])+0.4841229182759271*f[25]-0.2999999999999998*f[24]-0.2999999999999999*f[22]+0.375*(f[21]+f[19])-0.7794228634059945*f[18]+0.5590169943749475*f[14]+0.223606797749979*(f[13]+f[12])-0.2795084971874737*f[11]+0.5809475019311124*(f[10]+f[9])+0.45*f[7]-0.4330127018922193*f[4]-0.3354101966249685*(f[3]+f[2])+0.25*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_1x3v_p2_surfvz_quad_2(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return 1.349999999999999*f[47]-0.697137002317335*(f[46]+f[45])+0.697137002317335*f[44]+1.006230589874905*f[43]-1.006230589874905*(f[42]+f[41])-0.5196152422706628*f[40]+0.5196152422706628*f[39]-0.5196152422706631*(f[38]+f[37])+0.5196152422706628*f[36]-0.5196152422706631*f[35]-0.4024922359499621*(f[34]+f[33])+0.4024922359499621*f[32]+1.045705503476002*f[31]-0.75*(f[30]+f[29])+0.75*f[28]+0.3872983346207416*(f[27]+f[26]+f[25])-0.2999999999999998*f[24]+0.2999999999999998*f[23]-0.2999999999999999*(f[22]+f[21])+0.2999999999999998*f[20]-0.2999999999999999*f[19]+0.7794228634059945*f[18]-0.7794228634059945*(f[17]+f[16])+0.6037383539249431*f[15]+0.5590169943749475*f[14]+0.223606797749979*(f[13]+f[12]+f[11])-0.5809475019311124*(f[10]+f[9])+0.5809475019311124*f[8]+0.45*f[7]-0.45*(f[6]+f[5])+0.4330127018922193*f[4]-0.3354101966249685*(f[3]+f[2])+0.3354101966249685*f[1]+0.25*f[0]; 
  else 
    return 1.349999999999999*f[47]+0.697137002317335*(f[46]+f[45])-0.697137002317335*f[44]+1.006230589874905*f[43]-1.006230589874905*(f[42]+f[41])+0.5196152422706628*f[40]-0.5196152422706628*f[39]+0.5196152422706631*(f[38]+f[37])-0.5196152422706628*f[36]+0.5196152422706631*f[35]-0.4024922359499621*(f[34]+f[33])+0.4024922359499621*f[32]-1.045705503476002*f[31]-0.75*(f[30]+f[29])+0.75*f[28]-0.3872983346207416*(f[27]+f[26]+f[25])-0.2999999999999998*f[24]+0.2999999999999998*f[23]-0.2999999999999999*(f[22]+f[21])+0.2999999999999998*f[20]-0.2999999999999999*f[19]-0.7794228634059945*f[18]+0.7794228634059945*(f[17]+f[16])+0.6037383539249431*f[15]+0.5590169943749475*f[14]+0.223606797749979*(f[13]+f[12]+f[11])+0.5809475019311124*(f[10]+f[9])-0.5809475019311124*f[8]+0.45*f[7]-0.45*(f[6]+f[5])-0.4330127018922193*f[4]-0.3354101966249685*(f[3]+f[2])+0.3354101966249685*f[1]+0.25*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_1x3v_p2_surfvz_quad_3(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return (-0.8714212528966688*f[45])+1.006230589874905*f[42]-0.5196152422706628*f[39]+0.6495190528383289*f[38]-0.5196152422706631*f[37]+0.6495190528383289*f[36]-0.5031152949374527*f[33]-0.75*(f[30]+f[28])+0.3872983346207416*f[27]-0.4841229182759271*f[26]+0.3872983346207416*f[25]-0.2999999999999998*f[23]+0.375*f[22]-0.2999999999999999*f[21]+0.375*f[20]+0.7794228634059945*f[17]+0.5590169943749475*f[14]+0.223606797749979*f[13]-0.2795084971874737*f[12]+0.223606797749979*f[11]-0.5809475019311124*(f[10]+f[8])+0.45*f[6]+0.4330127018922193*f[4]-0.3354101966249685*(f[3]+f[1])+0.25*f[0]; 
  else 
    return 0.8714212528966688*f[45]+1.006230589874905*f[42]+0.5196152422706628*f[39]-0.6495190528383289*f[38]+0.5196152422706631*f[37]-0.6495190528383289*f[36]-0.5031152949374527*f[33]-0.75*(f[30]+f[28])-0.3872983346207416*f[27]+0.4841229182759271*f[26]-0.3872983346207416*f[25]-0.2999999999999998*f[23]+0.375*f[22]-0.2999999999999999*f[21]+0.375*f[20]-0.7794228634059945*f[17]+0.5590169943749475*f[14]+0.223606797749979*f[13]-0.2795084971874737*f[12]+0.223606797749979*f[11]+0.5809475019311124*(f[10]+f[8])+0.45*f[6]-0.4330127018922193*f[4]-0.3354101966249685*(f[3]+f[1])+0.25*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_1x3v_p2_surfvz_quad_4(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return 0.6495190528383289*(f[38]+f[37])-0.75*f[30]+0.3872983346207416*f[27]-0.4841229182759271*(f[26]+f[25])+0.375*(f[22]+f[21])+0.5590169943749475*f[14]+0.223606797749979*f[13]-0.2795084971874737*(f[12]+f[11])-0.5809475019311124*f[10]+0.4330127018922193*f[4]-0.3354101966249685*f[3]+0.25*f[0]; 
  else 
    return (-0.6495190528383289*(f[38]+f[37]))-0.75*f[30]-0.3872983346207416*f[27]+0.4841229182759271*(f[26]+f[25])+0.375*(f[22]+f[21])+0.5590169943749475*f[14]+0.223606797749979*f[13]-0.2795084971874737*(f[12]+f[11])+0.5809475019311124*f[10]-0.4330127018922193*f[4]-0.3354101966249685*f[3]+0.25*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_1x3v_p2_surfvz_quad_5(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return 0.8714212528966688*f[45]-1.006230589874905*f[42]+0.5196152422706628*f[39]+0.6495190528383289*f[38]-0.5196152422706631*f[37]-0.6495190528383289*f[36]+0.5031152949374527*f[33]-0.75*f[30]+0.75*f[28]+0.3872983346207416*f[27]-0.4841229182759271*f[26]+0.3872983346207416*f[25]+0.2999999999999998*f[23]+0.375*f[22]-0.2999999999999999*f[21]-0.375*f[20]-0.7794228634059945*f[17]+0.5590169943749475*f[14]+0.223606797749979*f[13]-0.2795084971874737*f[12]+0.223606797749979*f[11]-0.5809475019311124*f[10]+0.5809475019311124*f[8]-0.45*f[6]+0.4330127018922193*f[4]-0.3354101966249685*f[3]+0.3354101966249685*f[1]+0.25*f[0]; 
  else 
    return (-0.8714212528966688*f[45])-1.006230589874905*f[42]-0.5196152422706628*f[39]-0.6495190528383289*f[38]+0.5196152422706631*f[37]+0.6495190528383289*f[36]+0.5031152949374527*f[33]-0.75*f[30]+0.75*f[28]-0.3872983346207416*f[27]+0.4841229182759271*f[26]-0.3872983346207416*f[25]+0.2999999999999998*f[23]+0.375*f[22]-0.2999999999999999*f[21]-0.375*f[20]+0.7794228634059945*f[17]+0.5590169943749475*f[14]+0.223606797749979*f[13]-0.2795084971874737*f[12]+0.223606797749979*f[11]+0.5809475019311124*f[10]-0.5809475019311124*f[8]-0.45*f[6]-0.4330127018922193*f[4]-0.3354101966249685*f[3]+0.3354101966249685*f[1]+0.25*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_1x3v_p2_surfvz_quad_6(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return 1.349999999999999*f[47]-0.697137002317335*f[46]+0.697137002317335*f[45]-0.697137002317335*f[44]-1.006230589874905*f[43]+1.006230589874905*f[42]-1.006230589874905*f[41]+0.5196152422706628*f[40]-0.5196152422706628*f[39]-0.5196152422706631*(f[38]+f[37])-0.5196152422706628*f[36]+0.5196152422706631*f[35]-0.4024922359499621*f[34]+0.4024922359499621*f[33]-0.4024922359499621*f[32]+1.045705503476002*f[31]-0.75*f[30]+0.75*f[29]-0.75*f[28]+0.3872983346207416*(f[27]+f[26]+f[25])+0.2999999999999998*f[24]-0.2999999999999998*f[23]-0.2999999999999999*(f[22]+f[21])-0.2999999999999998*f[20]+0.2999999999999999*f[19]-0.7794228634059945*f[18]+0.7794228634059945*f[17]-0.7794228634059945*f[16]+0.6037383539249431*f[15]+0.5590169943749475*f[14]+0.223606797749979*(f[13]+f[12]+f[11])-0.5809475019311124*f[10]+0.5809475019311124*f[9]-0.5809475019311124*f[8]-0.45*f[7]+0.45*f[6]-0.45*f[5]+0.4330127018922193*f[4]-0.3354101966249685*f[3]+0.3354101966249685*f[2]-0.3354101966249685*f[1]+0.25*f[0]; 
  else 
    return 1.349999999999999*f[47]+0.697137002317335*f[46]-0.697137002317335*f[45]+0.697137002317335*f[44]-1.006230589874905*f[43]+1.006230589874905*f[42]-1.006230589874905*f[41]-0.5196152422706628*f[40]+0.5196152422706628*f[39]+0.5196152422706631*(f[38]+f[37])+0.5196152422706628*f[36]-0.5196152422706631*f[35]-0.4024922359499621*f[34]+0.4024922359499621*f[33]-0.4024922359499621*f[32]-1.045705503476002*f[31]-0.75*f[30]+0.75*f[29]-0.75*f[28]-0.3872983346207416*(f[27]+f[26]+f[25])+0.2999999999999998*f[24]-0.2999999999999998*f[23]-0.2999999999999999*(f[22]+f[21])-0.2999999999999998*f[20]+0.2999999999999999*f[19]+0.7794228634059945*f[18]-0.7794228634059945*f[17]+0.7794228634059945*f[16]+0.6037383539249431*f[15]+0.5590169943749475*f[14]+0.223606797749979*(f[13]+f[12]+f[11])+0.5809475019311124*f[10]-0.5809475019311124*f[9]+0.5809475019311124*f[8]-0.45*f[7]+0.45*f[6]-0.45*f[5]-0.4330127018922193*f[4]-0.3354101966249685*f[3]+0.3354101966249685*f[2]-0.3354101966249685*f[1]+0.25*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_1x3v_p2_surfvz_quad_7(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return 0.8714212528966688*f[44]-1.006230589874905*f[43]+0.5196152422706628*f[40]-0.5196152422706631*f[38]+0.6495190528383289*f[37]-0.6495190528383289*f[35]+0.5031152949374527*f[32]-0.75*f[30]+0.75*f[29]+0.3872983346207416*(f[27]+f[26])-0.4841229182759271*f[25]+0.2999999999999998*f[24]-0.2999999999999999*f[22]+0.375*f[21]-0.375*f[19]-0.7794228634059945*f[18]+0.5590169943749475*f[14]+0.223606797749979*(f[13]+f[12])-0.2795084971874737*f[11]-0.5809475019311124*f[10]+0.5809475019311124*f[9]-0.45*f[7]+0.4330127018922193*f[4]-0.3354101966249685*f[3]+0.3354101966249685*f[2]+0.25*f[0]; 
  else 
    return (-0.8714212528966688*f[44])-1.006230589874905*f[43]-0.5196152422706628*f[40]+0.5196152422706631*f[38]-0.6495190528383289*f[37]+0.6495190528383289*f[35]+0.5031152949374527*f[32]-0.75*f[30]+0.75*f[29]-0.3872983346207416*(f[27]+f[26])+0.4841229182759271*f[25]+0.2999999999999998*f[24]-0.2999999999999999*f[22]+0.375*f[21]-0.375*f[19]+0.7794228634059945*f[18]+0.5590169943749475*f[14]+0.223606797749979*(f[13]+f[12])-0.2795084971874737*f[11]+0.5809475019311124*f[10]-0.5809475019311124*f[9]-0.45*f[7]-0.4330127018922193*f[4]-0.3354101966249685*f[3]+0.3354101966249685*f[2]+0.25*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_1x3v_p2_surfvz_quad_8(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return (-1.349999999999999*f[47])+0.697137002317335*f[46]-0.697137002317335*(f[45]+f[44])-1.006230589874905*(f[43]+f[42])+1.006230589874905*f[41]+0.5196152422706628*(f[40]+f[39])-0.5196152422706631*(f[38]+f[37])+0.5196152422706628*f[36]+0.5196152422706631*f[35]+0.4024922359499621*f[34]-0.4024922359499621*(f[33]+f[32])-1.045705503476002*f[31]-0.75*f[30]+0.75*(f[29]+f[28])+0.3872983346207416*(f[27]+f[26]+f[25])+0.2999999999999998*(f[24]+f[23])-0.2999999999999999*(f[22]+f[21])+0.2999999999999998*f[20]+0.2999999999999999*f[19]-0.7794228634059945*(f[18]+f[17])+0.7794228634059945*f[16]-0.6037383539249431*f[15]+0.5590169943749475*f[14]+0.223606797749979*(f[13]+f[12]+f[11])-0.5809475019311124*f[10]+0.5809475019311124*(f[9]+f[8])-0.45*(f[7]+f[6])+0.45*f[5]+0.4330127018922193*f[4]-0.3354101966249685*f[3]+0.3354101966249685*(f[2]+f[1])+0.25*f[0]; 
  else 
    return (-1.349999999999999*f[47])-0.697137002317335*f[46]+0.697137002317335*(f[45]+f[44])-1.006230589874905*(f[43]+f[42])+1.006230589874905*f[41]-0.5196152422706628*(f[40]+f[39])+0.5196152422706631*(f[38]+f[37])-0.5196152422706628*f[36]-0.5196152422706631*f[35]+0.4024922359499621*f[34]-0.4024922359499621*(f[33]+f[32])+1.045705503476002*f[31]-0.75*f[30]+0.75*(f[29]+f[28])-0.3872983346207416*(f[27]+f[26]+f[25])+0.2999999999999998*(f[24]+f[23])-0.2999999999999999*(f[22]+f[21])+0.2999999999999998*f[20]+0.2999999999999999*f[19]+0.7794228634059945*(f[18]+f[17])-0.7794228634059945*f[16]-0.6037383539249431*f[15]+0.5590169943749475*f[14]+0.223606797749979*(f[13]+f[12]+f[11])+0.5809475019311124*f[10]-0.5809475019311124*(f[9]+f[8])-0.45*(f[7]+f[6])+0.45*f[5]-0.4330127018922193*f[4]-0.3354101966249685*f[3]+0.3354101966249685*(f[2]+f[1])+0.25*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_1x3v_p2_surfvz_quad_9(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return (-0.8714212528966688*f[46])+1.006230589874905*f[41]+0.6495190528383289*(f[40]+f[39])-0.5196152422706628*f[36]-0.5196152422706631*f[35]-0.5031152949374527*f[34]-0.75*(f[29]+f[28])-0.4841229182759271*f[27]+0.3872983346207416*(f[26]+f[25])+0.375*(f[24]+f[23])-0.2999999999999998*f[20]-0.2999999999999999*f[19]+0.7794228634059945*f[16]+0.5590169943749475*f[14]-0.2795084971874737*f[13]+0.223606797749979*(f[12]+f[11])-0.5809475019311124*(f[9]+f[8])+0.45*f[5]+0.4330127018922193*f[4]-0.3354101966249685*(f[2]+f[1])+0.25*f[0]; 
  else 
    return 0.8714212528966688*f[46]+1.006230589874905*f[41]-0.6495190528383289*(f[40]+f[39])+0.5196152422706628*f[36]+0.5196152422706631*f[35]-0.5031152949374527*f[34]-0.75*(f[29]+f[28])+0.4841229182759271*f[27]-0.3872983346207416*(f[26]+f[25])+0.375*(f[24]+f[23])-0.2999999999999998*f[20]-0.2999999999999999*f[19]-0.7794228634059945*f[16]+0.5590169943749475*f[14]-0.2795084971874737*f[13]+0.223606797749979*(f[12]+f[11])+0.5809475019311124*(f[9]+f[8])+0.45*f[5]-0.4330127018922193*f[4]-0.3354101966249685*(f[2]+f[1])+0.25*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_1x3v_p2_surfvz_quad_10(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return 0.6495190528383289*(f[40]+f[35])-0.75*f[29]-0.4841229182759271*f[27]+0.3872983346207416*f[26]-0.4841229182759271*f[25]+0.375*(f[24]+f[19])+0.5590169943749475*f[14]-0.2795084971874737*f[13]+0.223606797749979*f[12]-0.2795084971874737*f[11]-0.5809475019311124*f[9]+0.4330127018922193*f[4]-0.3354101966249685*f[2]+0.25*f[0]; 
  else 
    return (-0.6495190528383289*(f[40]+f[35]))-0.75*f[29]+0.4841229182759271*f[27]-0.3872983346207416*f[26]+0.4841229182759271*f[25]+0.375*(f[24]+f[19])+0.5590169943749475*f[14]-0.2795084971874737*f[13]+0.223606797749979*f[12]-0.2795084971874737*f[11]+0.5809475019311124*f[9]-0.4330127018922193*f[4]-0.3354101966249685*f[2]+0.25*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_1x3v_p2_surfvz_quad_11(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return 0.8714212528966688*f[46]-1.006230589874905*f[41]+0.6495190528383289*f[40]-0.6495190528383289*f[39]+0.5196152422706628*f[36]-0.5196152422706631*f[35]+0.5031152949374527*f[34]-0.75*f[29]+0.75*f[28]-0.4841229182759271*f[27]+0.3872983346207416*(f[26]+f[25])+0.375*f[24]-0.375*f[23]+0.2999999999999998*f[20]-0.2999999999999999*f[19]-0.7794228634059945*f[16]+0.5590169943749475*f[14]-0.2795084971874737*f[13]+0.223606797749979*(f[12]+f[11])-0.5809475019311124*f[9]+0.5809475019311124*f[8]-0.45*f[5]+0.4330127018922193*f[4]-0.3354101966249685*f[2]+0.3354101966249685*f[1]+0.25*f[0]; 
  else 
    return (-0.8714212528966688*f[46])-1.006230589874905*f[41]-0.6495190528383289*f[40]+0.6495190528383289*f[39]-0.5196152422706628*f[36]+0.5196152422706631*f[35]+0.5031152949374527*f[34]-0.75*f[29]+0.75*f[28]+0.4841229182759271*f[27]-0.3872983346207416*(f[26]+f[25])+0.375*f[24]-0.375*f[23]+0.2999999999999998*f[20]-0.2999999999999999*f[19]+0.7794228634059945*f[16]+0.5590169943749475*f[14]-0.2795084971874737*f[13]+0.223606797749979*(f[12]+f[11])+0.5809475019311124*f[9]-0.5809475019311124*f[8]-0.45*f[5]-0.4330127018922193*f[4]-0.3354101966249685*f[2]+0.3354101966249685*f[1]+0.25*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_1x3v_p2_surfvz_quad_12(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return 0.6495190528383289*(f[39]+f[36])-0.75*f[28]-0.4841229182759271*(f[27]+f[26])+0.3872983346207416*f[25]+0.375*(f[23]+f[20])+0.5590169943749475*f[14]-0.2795084971874737*(f[13]+f[12])+0.223606797749979*f[11]-0.5809475019311124*f[8]+0.4330127018922193*f[4]-0.3354101966249685*f[1]+0.25*f[0]; 
  else 
    return (-0.6495190528383289*(f[39]+f[36]))-0.75*f[28]+0.4841229182759271*(f[27]+f[26])-0.3872983346207416*f[25]+0.375*(f[23]+f[20])+0.5590169943749475*f[14]-0.2795084971874737*(f[13]+f[12])+0.223606797749979*f[11]+0.5809475019311124*f[8]-0.4330127018922193*f[4]-0.3354101966249685*f[1]+0.25*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_1x3v_p2_surfvz_quad_13(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return (-0.4841229182759271*(f[27]+f[26]+f[25]))+0.5590169943749475*f[14]-0.2795084971874737*(f[13]+f[12]+f[11])+0.4330127018922193*f[4]+0.25*f[0]; 
  else 
    return 0.4841229182759271*(f[27]+f[26]+f[25])+0.5590169943749475*f[14]-0.2795084971874737*(f[13]+f[12]+f[11])-0.4330127018922193*f[4]+0.25*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_1x3v_p2_surfvz_quad_14(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return (-0.6495190528383289*(f[39]+f[36]))+0.75*f[28]-0.4841229182759271*(f[27]+f[26])+0.3872983346207416*f[25]-0.375*(f[23]+f[20])+0.5590169943749475*f[14]-0.2795084971874737*(f[13]+f[12])+0.223606797749979*f[11]+0.5809475019311124*f[8]+0.4330127018922193*f[4]+0.3354101966249685*f[1]+0.25*f[0]; 
  else 
    return 0.6495190528383289*(f[39]+f[36])+0.75*f[28]+0.4841229182759271*(f[27]+f[26])-0.3872983346207416*f[25]-0.375*(f[23]+f[20])+0.5590169943749475*f[14]-0.2795084971874737*(f[13]+f[12])+0.223606797749979*f[11]-0.5809475019311124*f[8]-0.4330127018922193*f[4]+0.3354101966249685*f[1]+0.25*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_1x3v_p2_surfvz_quad_15(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return 0.8714212528966688*f[46]-1.006230589874905*f[41]-0.6495190528383289*f[40]+0.6495190528383289*f[39]-0.5196152422706628*f[36]+0.5196152422706631*f[35]+0.5031152949374527*f[34]+0.75*f[29]-0.75*f[28]-0.4841229182759271*f[27]+0.3872983346207416*(f[26]+f[25])-0.375*f[24]+0.375*f[23]-0.2999999999999998*f[20]+0.2999999999999999*f[19]-0.7794228634059945*f[16]+0.5590169943749475*f[14]-0.2795084971874737*f[13]+0.223606797749979*(f[12]+f[11])+0.5809475019311124*f[9]-0.5809475019311124*f[8]-0.45*f[5]+0.4330127018922193*f[4]+0.3354101966249685*f[2]-0.3354101966249685*f[1]+0.25*f[0]; 
  else 
    return (-0.8714212528966688*f[46])-1.006230589874905*f[41]+0.6495190528383289*f[40]-0.6495190528383289*f[39]+0.5196152422706628*f[36]-0.5196152422706631*f[35]+0.5031152949374527*f[34]+0.75*f[29]-0.75*f[28]+0.4841229182759271*f[27]-0.3872983346207416*(f[26]+f[25])-0.375*f[24]+0.375*f[23]-0.2999999999999998*f[20]+0.2999999999999999*f[19]+0.7794228634059945*f[16]+0.5590169943749475*f[14]-0.2795084971874737*f[13]+0.223606797749979*(f[12]+f[11])-0.5809475019311124*f[9]+0.5809475019311124*f[8]-0.45*f[5]-0.4330127018922193*f[4]+0.3354101966249685*f[2]-0.3354101966249685*f[1]+0.25*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_1x3v_p2_surfvz_quad_16(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return (-0.6495190528383289*(f[40]+f[35]))+0.75*f[29]-0.4841229182759271*f[27]+0.3872983346207416*f[26]-0.4841229182759271*f[25]-0.375*(f[24]+f[19])+0.5590169943749475*f[14]-0.2795084971874737*f[13]+0.223606797749979*f[12]-0.2795084971874737*f[11]+0.5809475019311124*f[9]+0.4330127018922193*f[4]+0.3354101966249685*f[2]+0.25*f[0]; 
  else 
    return 0.6495190528383289*(f[40]+f[35])+0.75*f[29]+0.4841229182759271*f[27]-0.3872983346207416*f[26]+0.4841229182759271*f[25]-0.375*(f[24]+f[19])+0.5590169943749475*f[14]-0.2795084971874737*f[13]+0.223606797749979*f[12]-0.2795084971874737*f[11]-0.5809475019311124*f[9]-0.4330127018922193*f[4]+0.3354101966249685*f[2]+0.25*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_1x3v_p2_surfvz_quad_17(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return (-0.8714212528966688*f[46])+1.006230589874905*f[41]-0.6495190528383289*(f[40]+f[39])+0.5196152422706628*f[36]+0.5196152422706631*f[35]-0.5031152949374527*f[34]+0.75*(f[29]+f[28])-0.4841229182759271*f[27]+0.3872983346207416*(f[26]+f[25])-0.375*(f[24]+f[23])+0.2999999999999998*f[20]+0.2999999999999999*f[19]+0.7794228634059945*f[16]+0.5590169943749475*f[14]-0.2795084971874737*f[13]+0.223606797749979*(f[12]+f[11])+0.5809475019311124*(f[9]+f[8])+0.45*f[5]+0.4330127018922193*f[4]+0.3354101966249685*(f[2]+f[1])+0.25*f[0]; 
  else 
    return 0.8714212528966688*f[46]+1.006230589874905*f[41]+0.6495190528383289*(f[40]+f[39])-0.5196152422706628*f[36]-0.5196152422706631*f[35]-0.5031152949374527*f[34]+0.75*(f[29]+f[28])+0.4841229182759271*f[27]-0.3872983346207416*(f[26]+f[25])-0.375*(f[24]+f[23])+0.2999999999999998*f[20]+0.2999999999999999*f[19]-0.7794228634059945*f[16]+0.5590169943749475*f[14]-0.2795084971874737*f[13]+0.223606797749979*(f[12]+f[11])-0.5809475019311124*(f[9]+f[8])+0.45*f[5]-0.4330127018922193*f[4]+0.3354101966249685*(f[2]+f[1])+0.25*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_1x3v_p2_surfvz_quad_18(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return 1.349999999999999*f[47]+0.697137002317335*f[46]-0.697137002317335*(f[45]+f[44])-1.006230589874905*(f[43]+f[42])+1.006230589874905*f[41]-0.5196152422706628*(f[40]+f[39])+0.5196152422706631*(f[38]+f[37])-0.5196152422706628*f[36]-0.5196152422706631*f[35]+0.4024922359499621*f[34]-0.4024922359499621*(f[33]+f[32])+1.045705503476002*f[31]+0.75*f[30]-0.75*(f[29]+f[28])+0.3872983346207416*(f[27]+f[26]+f[25])-0.2999999999999998*(f[24]+f[23])+0.2999999999999999*(f[22]+f[21])-0.2999999999999998*f[20]-0.2999999999999999*f[19]-0.7794228634059945*(f[18]+f[17])+0.7794228634059945*f[16]+0.6037383539249431*f[15]+0.5590169943749475*f[14]+0.223606797749979*(f[13]+f[12]+f[11])+0.5809475019311124*f[10]-0.5809475019311124*(f[9]+f[8])-0.45*(f[7]+f[6])+0.45*f[5]+0.4330127018922193*f[4]+0.3354101966249685*f[3]-0.3354101966249685*(f[2]+f[1])+0.25*f[0]; 
  else 
    return 1.349999999999999*f[47]-0.697137002317335*f[46]+0.697137002317335*(f[45]+f[44])-1.006230589874905*(f[43]+f[42])+1.006230589874905*f[41]+0.5196152422706628*(f[40]+f[39])-0.5196152422706631*(f[38]+f[37])+0.5196152422706628*f[36]+0.5196152422706631*f[35]+0.4024922359499621*f[34]-0.4024922359499621*(f[33]+f[32])-1.045705503476002*f[31]+0.75*f[30]-0.75*(f[29]+f[28])-0.3872983346207416*(f[27]+f[26]+f[25])-0.2999999999999998*(f[24]+f[23])+0.2999999999999999*(f[22]+f[21])-0.2999999999999998*f[20]-0.2999999999999999*f[19]+0.7794228634059945*(f[18]+f[17])-0.7794228634059945*f[16]+0.6037383539249431*f[15]+0.5590169943749475*f[14]+0.223606797749979*(f[13]+f[12]+f[11])-0.5809475019311124*f[10]+0.5809475019311124*(f[9]+f[8])-0.45*(f[7]+f[6])+0.45*f[5]-0.4330127018922193*f[4]+0.3354101966249685*f[3]-0.3354101966249685*(f[2]+f[1])+0.25*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_1x3v_p2_surfvz_quad_19(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return 0.8714212528966688*f[44]-1.006230589874905*f[43]-0.5196152422706628*f[40]+0.5196152422706631*f[38]-0.6495190528383289*f[37]+0.6495190528383289*f[35]+0.5031152949374527*f[32]+0.75*f[30]-0.75*f[29]+0.3872983346207416*(f[27]+f[26])-0.4841229182759271*f[25]-0.2999999999999998*f[24]+0.2999999999999999*f[22]-0.375*f[21]+0.375*f[19]-0.7794228634059945*f[18]+0.5590169943749475*f[14]+0.223606797749979*(f[13]+f[12])-0.2795084971874737*f[11]+0.5809475019311124*f[10]-0.5809475019311124*f[9]-0.45*f[7]+0.4330127018922193*f[4]+0.3354101966249685*f[3]-0.3354101966249685*f[2]+0.25*f[0]; 
  else 
    return (-0.8714212528966688*f[44])-1.006230589874905*f[43]+0.5196152422706628*f[40]-0.5196152422706631*f[38]+0.6495190528383289*f[37]-0.6495190528383289*f[35]+0.5031152949374527*f[32]+0.75*f[30]-0.75*f[29]-0.3872983346207416*(f[27]+f[26])+0.4841229182759271*f[25]-0.2999999999999998*f[24]+0.2999999999999999*f[22]-0.375*f[21]+0.375*f[19]+0.7794228634059945*f[18]+0.5590169943749475*f[14]+0.223606797749979*(f[13]+f[12])-0.2795084971874737*f[11]-0.5809475019311124*f[10]+0.5809475019311124*f[9]-0.45*f[7]-0.4330127018922193*f[4]+0.3354101966249685*f[3]-0.3354101966249685*f[2]+0.25*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_1x3v_p2_surfvz_quad_20(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return (-1.349999999999999*f[47])-0.697137002317335*f[46]+0.697137002317335*f[45]-0.697137002317335*f[44]-1.006230589874905*f[43]+1.006230589874905*f[42]-1.006230589874905*f[41]-0.5196152422706628*f[40]+0.5196152422706628*f[39]+0.5196152422706631*(f[38]+f[37])+0.5196152422706628*f[36]-0.5196152422706631*f[35]-0.4024922359499621*f[34]+0.4024922359499621*f[33]-0.4024922359499621*f[32]-1.045705503476002*f[31]+0.75*f[30]-0.75*f[29]+0.75*f[28]+0.3872983346207416*(f[27]+f[26]+f[25])-0.2999999999999998*f[24]+0.2999999999999998*f[23]+0.2999999999999999*(f[22]+f[21])+0.2999999999999998*f[20]-0.2999999999999999*f[19]-0.7794228634059945*f[18]+0.7794228634059945*f[17]-0.7794228634059945*f[16]-0.6037383539249431*f[15]+0.5590169943749475*f[14]+0.223606797749979*(f[13]+f[12]+f[11])+0.5809475019311124*f[10]-0.5809475019311124*f[9]+0.5809475019311124*f[8]-0.45*f[7]+0.45*f[6]-0.45*f[5]+0.4330127018922193*f[4]+0.3354101966249685*f[3]-0.3354101966249685*f[2]+0.3354101966249685*f[1]+0.25*f[0]; 
  else 
    return (-1.349999999999999*f[47])+0.697137002317335*f[46]-0.697137002317335*f[45]+0.697137002317335*f[44]-1.006230589874905*f[43]+1.006230589874905*f[42]-1.006230589874905*f[41]+0.5196152422706628*f[40]-0.5196152422706628*f[39]-0.5196152422706631*(f[38]+f[37])-0.5196152422706628*f[36]+0.5196152422706631*f[35]-0.4024922359499621*f[34]+0.4024922359499621*f[33]-0.4024922359499621*f[32]+1.045705503476002*f[31]+0.75*f[30]-0.75*f[29]+0.75*f[28]-0.3872983346207416*(f[27]+f[26]+f[25])-0.2999999999999998*f[24]+0.2999999999999998*f[23]+0.2999999999999999*(f[22]+f[21])+0.2999999999999998*f[20]-0.2999999999999999*f[19]+0.7794228634059945*f[18]-0.7794228634059945*f[17]+0.7794228634059945*f[16]-0.6037383539249431*f[15]+0.5590169943749475*f[14]+0.223606797749979*(f[13]+f[12]+f[11])-0.5809475019311124*f[10]+0.5809475019311124*f[9]-0.5809475019311124*f[8]-0.45*f[7]+0.45*f[6]-0.45*f[5]-0.4330127018922193*f[4]+0.3354101966249685*f[3]-0.3354101966249685*f[2]+0.3354101966249685*f[1]+0.25*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_1x3v_p2_surfvz_quad_21(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return 0.8714212528966688*f[45]-1.006230589874905*f[42]-0.5196152422706628*f[39]-0.6495190528383289*f[38]+0.5196152422706631*f[37]+0.6495190528383289*f[36]+0.5031152949374527*f[33]+0.75*f[30]-0.75*f[28]+0.3872983346207416*f[27]-0.4841229182759271*f[26]+0.3872983346207416*f[25]-0.2999999999999998*f[23]-0.375*f[22]+0.2999999999999999*f[21]+0.375*f[20]-0.7794228634059945*f[17]+0.5590169943749475*f[14]+0.223606797749979*f[13]-0.2795084971874737*f[12]+0.223606797749979*f[11]+0.5809475019311124*f[10]-0.5809475019311124*f[8]-0.45*f[6]+0.4330127018922193*f[4]+0.3354101966249685*f[3]-0.3354101966249685*f[1]+0.25*f[0]; 
  else 
    return (-0.8714212528966688*f[45])-1.006230589874905*f[42]+0.5196152422706628*f[39]+0.6495190528383289*f[38]-0.5196152422706631*f[37]-0.6495190528383289*f[36]+0.5031152949374527*f[33]+0.75*f[30]-0.75*f[28]-0.3872983346207416*f[27]+0.4841229182759271*f[26]-0.3872983346207416*f[25]-0.2999999999999998*f[23]-0.375*f[22]+0.2999999999999999*f[21]+0.375*f[20]+0.7794228634059945*f[17]+0.5590169943749475*f[14]+0.223606797749979*f[13]-0.2795084971874737*f[12]+0.223606797749979*f[11]-0.5809475019311124*f[10]+0.5809475019311124*f[8]-0.45*f[6]-0.4330127018922193*f[4]+0.3354101966249685*f[3]-0.3354101966249685*f[1]+0.25*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_1x3v_p2_surfvz_quad_22(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return (-0.6495190528383289*(f[38]+f[37]))+0.75*f[30]+0.3872983346207416*f[27]-0.4841229182759271*(f[26]+f[25])-0.375*(f[22]+f[21])+0.5590169943749475*f[14]+0.223606797749979*f[13]-0.2795084971874737*(f[12]+f[11])+0.5809475019311124*f[10]+0.4330127018922193*f[4]+0.3354101966249685*f[3]+0.25*f[0]; 
  else 
    return 0.6495190528383289*(f[38]+f[37])+0.75*f[30]-0.3872983346207416*f[27]+0.4841229182759271*(f[26]+f[25])-0.375*(f[22]+f[21])+0.5590169943749475*f[14]+0.223606797749979*f[13]-0.2795084971874737*(f[12]+f[11])-0.5809475019311124*f[10]-0.4330127018922193*f[4]+0.3354101966249685*f[3]+0.25*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_1x3v_p2_surfvz_quad_23(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return (-0.8714212528966688*f[45])+1.006230589874905*f[42]+0.5196152422706628*f[39]-0.6495190528383289*f[38]+0.5196152422706631*f[37]-0.6495190528383289*f[36]-0.5031152949374527*f[33]+0.75*(f[30]+f[28])+0.3872983346207416*f[27]-0.4841229182759271*f[26]+0.3872983346207416*f[25]+0.2999999999999998*f[23]-0.375*f[22]+0.2999999999999999*f[21]-0.375*f[20]+0.7794228634059945*f[17]+0.5590169943749475*f[14]+0.223606797749979*f[13]-0.2795084971874737*f[12]+0.223606797749979*f[11]+0.5809475019311124*(f[10]+f[8])+0.45*f[6]+0.4330127018922193*f[4]+0.3354101966249685*(f[3]+f[1])+0.25*f[0]; 
  else 
    return 0.8714212528966688*f[45]+1.006230589874905*f[42]-0.5196152422706628*f[39]+0.6495190528383289*f[38]-0.5196152422706631*f[37]+0.6495190528383289*f[36]-0.5031152949374527*f[33]+0.75*(f[30]+f[28])-0.3872983346207416*f[27]+0.4841229182759271*f[26]-0.3872983346207416*f[25]+0.2999999999999998*f[23]-0.375*f[22]+0.2999999999999999*f[21]-0.375*f[20]-0.7794228634059945*f[17]+0.5590169943749475*f[14]+0.223606797749979*f[13]-0.2795084971874737*f[12]+0.223606797749979*f[11]-0.5809475019311124*(f[10]+f[8])+0.45*f[6]-0.4330127018922193*f[4]+0.3354101966249685*(f[3]+f[1])+0.25*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_1x3v_p2_surfvz_quad_24(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return (-1.349999999999999*f[47])-0.697137002317335*(f[46]+f[45])+0.697137002317335*f[44]+1.006230589874905*f[43]-1.006230589874905*(f[42]+f[41])+0.5196152422706628*f[40]-0.5196152422706628*f[39]+0.5196152422706631*(f[38]+f[37])-0.5196152422706628*f[36]+0.5196152422706631*f[35]-0.4024922359499621*(f[34]+f[33])+0.4024922359499621*f[32]-1.045705503476002*f[31]+0.75*(f[30]+f[29])-0.75*f[28]+0.3872983346207416*(f[27]+f[26]+f[25])+0.2999999999999998*f[24]-0.2999999999999998*f[23]+0.2999999999999999*(f[22]+f[21])-0.2999999999999998*f[20]+0.2999999999999999*f[19]+0.7794228634059945*f[18]-0.7794228634059945*(f[17]+f[16])-0.6037383539249431*f[15]+0.5590169943749475*f[14]+0.223606797749979*(f[13]+f[12]+f[11])+0.5809475019311124*(f[10]+f[9])-0.5809475019311124*f[8]+0.45*f[7]-0.45*(f[6]+f[5])+0.4330127018922193*f[4]+0.3354101966249685*(f[3]+f[2])-0.3354101966249685*f[1]+0.25*f[0]; 
  else 
    return (-1.349999999999999*f[47])+0.697137002317335*(f[46]+f[45])-0.697137002317335*f[44]+1.006230589874905*f[43]-1.006230589874905*(f[42]+f[41])-0.5196152422706628*f[40]+0.5196152422706628*f[39]-0.5196152422706631*(f[38]+f[37])+0.5196152422706628*f[36]-0.5196152422706631*f[35]-0.4024922359499621*(f[34]+f[33])+0.4024922359499621*f[32]+1.045705503476002*f[31]+0.75*(f[30]+f[29])-0.75*f[28]-0.3872983346207416*(f[27]+f[26]+f[25])+0.2999999999999998*f[24]-0.2999999999999998*f[23]+0.2999999999999999*(f[22]+f[21])-0.2999999999999998*f[20]+0.2999999999999999*f[19]-0.7794228634059945*f[18]+0.7794228634059945*(f[17]+f[16])-0.6037383539249431*f[15]+0.5590169943749475*f[14]+0.223606797749979*(f[13]+f[12]+f[11])-0.5809475019311124*(f[10]+f[9])+0.5809475019311124*f[8]+0.45*f[7]-0.45*(f[6]+f[5])-0.4330127018922193*f[4]+0.3354101966249685*(f[3]+f[2])-0.3354101966249685*f[1]+0.25*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_1x3v_p2_surfvz_quad_25(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return (-0.8714212528966688*f[44])+1.006230589874905*f[43]+0.5196152422706628*f[40]+0.5196152422706631*f[38]-0.6495190528383289*(f[37]+f[35])-0.5031152949374527*f[32]+0.75*(f[30]+f[29])+0.3872983346207416*(f[27]+f[26])-0.4841229182759271*f[25]+0.2999999999999998*f[24]+0.2999999999999999*f[22]-0.375*(f[21]+f[19])+0.7794228634059945*f[18]+0.5590169943749475*f[14]+0.223606797749979*(f[13]+f[12])-0.2795084971874737*f[11]+0.5809475019311124*(f[10]+f[9])+0.45*f[7]+0.4330127018922193*f[4]+0.3354101966249685*(f[3]+f[2])+0.25*f[0]; 
  else 
    return 0.8714212528966688*f[44]+1.006230589874905*f[43]-0.5196152422706628*f[40]-0.5196152422706631*f[38]+0.6495190528383289*(f[37]+f[35])-0.5031152949374527*f[32]+0.75*(f[30]+f[29])-0.3872983346207416*(f[27]+f[26])+0.4841229182759271*f[25]+0.2999999999999998*f[24]+0.2999999999999999*f[22]-0.375*(f[21]+f[19])-0.7794228634059945*f[18]+0.5590169943749475*f[14]+0.223606797749979*(f[13]+f[12])-0.2795084971874737*f[11]-0.5809475019311124*(f[10]+f[9])+0.45*f[7]-0.4330127018922193*f[4]+0.3354101966249685*(f[3]+f[2])+0.25*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_1x3v_p2_surfvz_quad_26(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return 1.349999999999999*f[47]+0.697137002317335*(f[46]+f[45]+f[44])+1.006230589874905*(f[43]+f[42]+f[41])+0.5196152422706628*(f[40]+f[39])+0.5196152422706631*(f[38]+f[37])+0.5196152422706628*f[36]+0.5196152422706631*f[35]+0.4024922359499621*(f[34]+f[33]+f[32])+1.045705503476002*f[31]+0.75*(f[30]+f[29]+f[28])+0.3872983346207416*(f[27]+f[26]+f[25])+0.2999999999999998*(f[24]+f[23])+0.2999999999999999*(f[22]+f[21])+0.2999999999999998*f[20]+0.2999999999999999*f[19]+0.7794228634059945*(f[18]+f[17]+f[16])+0.6037383539249431*f[15]+0.5590169943749475*f[14]+0.223606797749979*(f[13]+f[12]+f[11])+0.5809475019311124*(f[10]+f[9]+f[8])+0.45*(f[7]+f[6]+f[5])+0.4330127018922193*f[4]+0.3354101966249685*(f[3]+f[2]+f[1])+0.25*f[0]; 
  else 
    return 1.349999999999999*f[47]-0.697137002317335*(f[46]+f[45]+f[44])+1.006230589874905*(f[43]+f[42]+f[41])-0.5196152422706628*(f[40]+f[39])-0.5196152422706631*(f[38]+f[37])-0.5196152422706628*f[36]-0.5196152422706631*f[35]+0.4024922359499621*(f[34]+f[33]+f[32])-1.045705503476002*f[31]+0.75*(f[30]+f[29]+f[28])-0.3872983346207416*(f[27]+f[26]+f[25])+0.2999999999999998*(f[24]+f[23])+0.2999999999999999*(f[22]+f[21])+0.2999999999999998*f[20]+0.2999999999999999*f[19]-0.7794228634059945*(f[18]+f[17]+f[16])+0.6037383539249431*f[15]+0.5590169943749475*f[14]+0.223606797749979*(f[13]+f[12]+f[11])-0.5809475019311124*(f[10]+f[9]+f[8])+0.45*(f[7]+f[6]+f[5])-0.4330127018922193*f[4]+0.3354101966249685*(f[3]+f[2]+f[1])+0.25*f[0]; 
} 
