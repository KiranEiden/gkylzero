GKYL_CU_DH static inline void 
tensor_3x_p2_sqrt(const double *A, double *ASqrt) 
{ 
  // A:     Input DG field. 
  // ASqrt: Output DG field (expansion of sqrt(A)). 
 
  double AOrd[27] = {0.0}; 

  AOrd[0] = sqrt(0.2529822128134704*A[26]-0.3794733192202044*(A[25]+A[24])-0.3794733192202055*A[23]+0.2828427124746191*(A[22]+A[21]+A[20])+0.5692099788303082*(A[19]+A[18]+A[17])-0.4242640687119281*(A[16]+A[15])-0.4242640687119285*(A[14]+A[13])-0.4242640687119281*A[12]-0.4242640687119285*A[11]-0.853814968245462*A[10]+0.3162277660168379*(A[9]+A[8]+A[7])+0.6363961030678926*(A[6]+A[5]+A[4])-0.4743416490252568*(A[3]+A[2]+A[1])+0.3535533905932737*A[0]); 
  AOrd[1] = sqrt((-0.3162277660168379*A[26])+0.4743416490252562*A[25]+0.4743416490252562*A[24]-0.3535533905932737*(A[22]+A[21])+0.2828427124746191*A[20]-0.711512473537885*A[19]+0.5303300858899104*(A[16]+A[15])-0.4242640687119281*A[12]-0.4242640687119285*A[11]-0.3952847075210473*A[9]+0.3162277660168379*(A[8]+A[7])+0.6363961030678926*A[4]-0.4743416490252568*(A[2]+A[1])+0.3535533905932737*A[0]); 
  AOrd[2] = sqrt(0.2529822128134704*A[26]-0.3794733192202044*(A[25]+A[24])+0.3794733192202055*A[23]+0.2828427124746191*(A[22]+A[21]+A[20])+0.5692099788303082*A[19]-0.5692099788303082*(A[18]+A[17])-0.4242640687119281*(A[16]+A[15])+0.4242640687119285*(A[14]+A[13])-0.4242640687119281*A[12]-0.4242640687119285*A[11]+0.853814968245462*A[10]+0.3162277660168379*(A[9]+A[8]+A[7])-0.6363961030678926*(A[6]+A[5])+0.6363961030678926*A[4]+0.4743416490252568*A[3]-0.4743416490252568*(A[2]+A[1])+0.3535533905932737*A[0]); 
  AOrd[3] = sqrt((-0.3162277660168379*A[26])+0.4743416490252562*A[25]+0.4743416490252568*A[23]-0.3535533905932737*A[22]+0.2828427124746191*A[21]-0.3535533905932737*A[20]-0.711512473537885*A[18]-0.4242640687119281*A[15]+0.5303300858899104*A[14]-0.4242640687119285*A[13]+0.5303300858899104*A[12]+0.3162277660168379*A[9]-0.3952847075210473*A[8]+0.3162277660168379*A[7]+0.6363961030678926*A[5]-0.4743416490252568*(A[3]+A[1])+0.3535533905932737*A[0]); 
  AOrd[4] = sqrt(0.3952847075210473*A[26]-0.592927061281571*A[25]+0.441941738241592*A[22]-0.3535533905932737*(A[21]+A[20])+0.5303300858899104*(A[15]+A[12])-0.3952847075210473*(A[9]+A[8])+0.3162277660168379*A[7]-0.4743416490252568*A[1]+0.3535533905932737*A[0]); 
  AOrd[5] = sqrt((-0.3162277660168379*A[26])+0.4743416490252562*A[25]-0.4743416490252568*A[23]-0.3535533905932737*A[22]+0.2828427124746191*A[21]-0.3535533905932737*A[20]+0.711512473537885*A[18]-0.4242640687119281*A[15]-0.5303300858899104*A[14]+0.4242640687119285*A[13]+0.5303300858899104*A[12]+0.3162277660168379*A[9]-0.3952847075210473*A[8]+0.3162277660168379*A[7]-0.6363961030678926*A[5]+0.4743416490252568*A[3]-0.4743416490252568*A[1]+0.3535533905932737*A[0]); 
  AOrd[6] = sqrt(0.2529822128134704*A[26]-0.3794733192202044*A[25]+0.3794733192202044*A[24]-0.3794733192202055*A[23]+0.2828427124746191*(A[22]+A[21]+A[20])-0.5692099788303082*A[19]+0.5692099788303082*A[18]-0.5692099788303082*A[17]+0.4242640687119281*A[16]-0.4242640687119281*A[15]-0.4242640687119285*(A[14]+A[13])-0.4242640687119281*A[12]+0.4242640687119285*A[11]+0.853814968245462*A[10]+0.3162277660168379*(A[9]+A[8]+A[7])-0.6363961030678926*A[6]+0.6363961030678926*A[5]-0.6363961030678926*A[4]-0.4743416490252568*A[3]+0.4743416490252568*A[2]-0.4743416490252568*A[1]+0.3535533905932737*A[0]); 
  AOrd[7] = sqrt((-0.3162277660168379*A[26])+0.4743416490252562*A[25]-0.4743416490252562*A[24]-0.3535533905932737*(A[22]+A[21])+0.2828427124746191*A[20]+0.711512473537885*A[19]-0.5303300858899104*A[16]+0.5303300858899104*A[15]-0.4242640687119281*A[12]+0.4242640687119285*A[11]-0.3952847075210473*A[9]+0.3162277660168379*(A[8]+A[7])-0.6363961030678926*A[4]+0.4743416490252568*A[2]-0.4743416490252568*A[1]+0.3535533905932737*A[0]); 
  AOrd[8] = sqrt(0.2529822128134704*A[26]-0.3794733192202044*A[25]+0.3794733192202044*A[24]+0.3794733192202055*A[23]+0.2828427124746191*(A[22]+A[21]+A[20])-0.5692099788303082*(A[19]+A[18])+0.5692099788303082*A[17]+0.4242640687119281*A[16]-0.4242640687119281*A[15]+0.4242640687119285*(A[14]+A[13])-0.4242640687119281*A[12]+0.4242640687119285*A[11]-0.853814968245462*A[10]+0.3162277660168379*(A[9]+A[8]+A[7])+0.6363961030678926*A[6]-0.6363961030678926*(A[5]+A[4])+0.4743416490252568*(A[3]+A[2])-0.4743416490252568*A[1]+0.3535533905932737*A[0]); 
  AOrd[9] = sqrt((-0.3162277660168379*A[26])+0.4743416490252562*A[24]+0.4743416490252568*A[23]+0.2828427124746191*A[22]-0.3535533905932737*(A[21]+A[20])-0.711512473537885*A[17]-0.4242640687119281*A[16]-0.4242640687119285*A[14]+0.5303300858899104*(A[13]+A[11])+0.3162277660168379*(A[9]+A[8])-0.3952847075210473*A[7]+0.6363961030678926*A[6]-0.4743416490252568*(A[3]+A[2])+0.3535533905932737*A[0]); 
  AOrd[10] = sqrt(0.3952847075210473*A[26]-0.592927061281571*A[24]-0.3535533905932737*A[22]+0.441941738241592*A[21]-0.3535533905932737*A[20]+0.5303300858899104*(A[16]+A[11])-0.3952847075210473*A[9]+0.3162277660168379*A[8]-0.3952847075210473*A[7]-0.4743416490252568*A[2]+0.3535533905932737*A[0]); 
  AOrd[11] = sqrt((-0.3162277660168379*A[26])+0.4743416490252562*A[24]-0.4743416490252568*A[23]+0.2828427124746191*A[22]-0.3535533905932737*(A[21]+A[20])+0.711512473537885*A[17]-0.4242640687119281*A[16]+0.4242640687119285*A[14]-0.5303300858899104*A[13]+0.5303300858899104*A[11]+0.3162277660168379*(A[9]+A[8])-0.3952847075210473*A[7]-0.6363961030678926*A[6]+0.4743416490252568*A[3]-0.4743416490252568*A[2]+0.3535533905932737*A[0]); 
  AOrd[12] = sqrt(0.3952847075210473*A[26]-0.592927061281571*A[23]-0.3535533905932737*(A[22]+A[21])+0.441941738241592*A[20]+0.5303300858899104*(A[14]+A[13])+0.3162277660168379*A[9]-0.3952847075210473*(A[8]+A[7])-0.4743416490252568*A[3]+0.3535533905932737*A[0]); 
  AOrd[13] = sqrt((-0.4941058844013091*A[26])+0.441941738241592*(A[22]+A[21]+A[20])-0.3952847075210473*(A[9]+A[8]+A[7])+0.3535533905932737*A[0]); 
  AOrd[14] = sqrt(0.3952847075210473*A[26]+0.592927061281571*A[23]-0.3535533905932737*(A[22]+A[21])+0.441941738241592*A[20]-0.5303300858899104*(A[14]+A[13])+0.3162277660168379*A[9]-0.3952847075210473*(A[8]+A[7])+0.4743416490252568*A[3]+0.3535533905932737*A[0]); 
  AOrd[15] = sqrt((-0.3162277660168379*A[26])-0.4743416490252562*A[24]+0.4743416490252568*A[23]+0.2828427124746191*A[22]-0.3535533905932737*(A[21]+A[20])+0.711512473537885*A[17]+0.4242640687119281*A[16]-0.4242640687119285*A[14]+0.5303300858899104*A[13]-0.5303300858899104*A[11]+0.3162277660168379*(A[9]+A[8])-0.3952847075210473*A[7]-0.6363961030678926*A[6]-0.4743416490252568*A[3]+0.4743416490252568*A[2]+0.3535533905932737*A[0]); 
  AOrd[16] = sqrt(0.3952847075210473*A[26]+0.592927061281571*A[24]-0.3535533905932737*A[22]+0.441941738241592*A[21]-0.3535533905932737*A[20]-0.5303300858899104*(A[16]+A[11])-0.3952847075210473*A[9]+0.3162277660168379*A[8]-0.3952847075210473*A[7]+0.4743416490252568*A[2]+0.3535533905932737*A[0]); 
  AOrd[17] = sqrt((-0.3162277660168379*A[26])-0.4743416490252562*A[24]-0.4743416490252568*A[23]+0.2828427124746191*A[22]-0.3535533905932737*(A[21]+A[20])-0.711512473537885*A[17]+0.4242640687119281*A[16]+0.4242640687119285*A[14]-0.5303300858899104*(A[13]+A[11])+0.3162277660168379*(A[9]+A[8])-0.3952847075210473*A[7]+0.6363961030678926*A[6]+0.4743416490252568*(A[3]+A[2])+0.3535533905932737*A[0]); 
  AOrd[18] = sqrt(0.2529822128134704*A[26]+0.3794733192202044*A[25]-0.3794733192202044*A[24]-0.3794733192202055*A[23]+0.2828427124746191*(A[22]+A[21]+A[20])-0.5692099788303082*(A[19]+A[18])+0.5692099788303082*A[17]-0.4242640687119281*A[16]+0.4242640687119281*A[15]-0.4242640687119285*(A[14]+A[13])+0.4242640687119281*A[12]-0.4242640687119285*A[11]+0.853814968245462*A[10]+0.3162277660168379*(A[9]+A[8]+A[7])+0.6363961030678926*A[6]-0.6363961030678926*(A[5]+A[4])-0.4743416490252568*(A[3]+A[2])+0.4743416490252568*A[1]+0.3535533905932737*A[0]); 
  AOrd[19] = sqrt((-0.3162277660168379*A[26])-0.4743416490252562*A[25]+0.4743416490252562*A[24]-0.3535533905932737*(A[22]+A[21])+0.2828427124746191*A[20]+0.711512473537885*A[19]+0.5303300858899104*A[16]-0.5303300858899104*A[15]+0.4242640687119281*A[12]-0.4242640687119285*A[11]-0.3952847075210473*A[9]+0.3162277660168379*(A[8]+A[7])-0.6363961030678926*A[4]-0.4743416490252568*A[2]+0.4743416490252568*A[1]+0.3535533905932737*A[0]); 
  AOrd[20] = sqrt(0.2529822128134704*A[26]+0.3794733192202044*A[25]-0.3794733192202044*A[24]+0.3794733192202055*A[23]+0.2828427124746191*(A[22]+A[21]+A[20])-0.5692099788303082*A[19]+0.5692099788303082*A[18]-0.5692099788303082*A[17]-0.4242640687119281*A[16]+0.4242640687119281*A[15]+0.4242640687119285*(A[14]+A[13])+0.4242640687119281*A[12]-0.4242640687119285*A[11]-0.853814968245462*A[10]+0.3162277660168379*(A[9]+A[8]+A[7])-0.6363961030678926*A[6]+0.6363961030678926*A[5]-0.6363961030678926*A[4]+0.4743416490252568*A[3]-0.4743416490252568*A[2]+0.4743416490252568*A[1]+0.3535533905932737*A[0]); 
  AOrd[21] = sqrt((-0.3162277660168379*A[26])-0.4743416490252562*A[25]+0.4743416490252568*A[23]-0.3535533905932737*A[22]+0.2828427124746191*A[21]-0.3535533905932737*A[20]+0.711512473537885*A[18]+0.4242640687119281*A[15]+0.5303300858899104*A[14]-0.4242640687119285*A[13]-0.5303300858899104*A[12]+0.3162277660168379*A[9]-0.3952847075210473*A[8]+0.3162277660168379*A[7]-0.6363961030678926*A[5]-0.4743416490252568*A[3]+0.4743416490252568*A[1]+0.3535533905932737*A[0]); 
  AOrd[22] = sqrt(0.3952847075210473*A[26]+0.592927061281571*A[25]+0.441941738241592*A[22]-0.3535533905932737*(A[21]+A[20])-0.5303300858899104*(A[15]+A[12])-0.3952847075210473*(A[9]+A[8])+0.3162277660168379*A[7]+0.4743416490252568*A[1]+0.3535533905932737*A[0]); 
  AOrd[23] = sqrt((-0.3162277660168379*A[26])-0.4743416490252562*A[25]-0.4743416490252568*A[23]-0.3535533905932737*A[22]+0.2828427124746191*A[21]-0.3535533905932737*A[20]-0.711512473537885*A[18]+0.4242640687119281*A[15]-0.5303300858899104*A[14]+0.4242640687119285*A[13]-0.5303300858899104*A[12]+0.3162277660168379*A[9]-0.3952847075210473*A[8]+0.3162277660168379*A[7]+0.6363961030678926*A[5]+0.4743416490252568*(A[3]+A[1])+0.3535533905932737*A[0]); 
  AOrd[24] = sqrt(0.2529822128134704*A[26]+0.3794733192202044*(A[25]+A[24])-0.3794733192202055*A[23]+0.2828427124746191*(A[22]+A[21]+A[20])+0.5692099788303082*A[19]-0.5692099788303082*(A[18]+A[17])+0.4242640687119281*(A[16]+A[15])-0.4242640687119285*(A[14]+A[13])+0.4242640687119281*A[12]+0.4242640687119285*A[11]-0.853814968245462*A[10]+0.3162277660168379*(A[9]+A[8]+A[7])-0.6363961030678926*(A[6]+A[5])+0.6363961030678926*A[4]-0.4743416490252568*A[3]+0.4743416490252568*(A[2]+A[1])+0.3535533905932737*A[0]); 
  AOrd[25] = sqrt((-0.3162277660168379*A[26])-0.4743416490252562*A[25]-0.4743416490252562*A[24]-0.3535533905932737*(A[22]+A[21])+0.2828427124746191*A[20]-0.711512473537885*A[19]-0.5303300858899104*(A[16]+A[15])+0.4242640687119281*A[12]+0.4242640687119285*A[11]-0.3952847075210473*A[9]+0.3162277660168379*(A[8]+A[7])+0.6363961030678926*A[4]+0.4743416490252568*(A[2]+A[1])+0.3535533905932737*A[0]); 
  AOrd[26] = sqrt(0.2529822128134704*A[26]+0.3794733192202044*(A[25]+A[24])+0.3794733192202055*A[23]+0.2828427124746191*(A[22]+A[21]+A[20])+0.5692099788303082*(A[19]+A[18]+A[17])+0.4242640687119281*(A[16]+A[15])+0.4242640687119285*(A[14]+A[13])+0.4242640687119281*A[12]+0.4242640687119285*A[11]+0.853814968245462*A[10]+0.3162277660168379*(A[9]+A[8]+A[7])+0.6363961030678926*(A[6]+A[5]+A[4])+0.4743416490252568*(A[3]+A[2]+A[1])+0.3535533905932737*A[0]); 
  ASqrt[0] = 0.06062300936098657*AOrd[26]+0.09699681497757856*AOrd[25]+0.06062300936098657*AOrd[24]+0.09699681497757856*AOrd[23]+0.1551949039641257*AOrd[22]+0.09699681497757856*AOrd[21]+0.06062300936098657*AOrd[20]+0.09699681497757856*AOrd[19]+0.06062300936098657*AOrd[18]+0.09699681497757856*AOrd[17]+0.1551949039641257*AOrd[16]+0.09699681497757856*AOrd[15]+0.1551949039641257*AOrd[14]+0.2483118463426013*AOrd[13]+0.1551949039641257*AOrd[12]+0.09699681497757856*AOrd[11]+0.1551949039641257*AOrd[10]+0.09699681497757856*AOrd[9]+0.06062300936098657*AOrd[8]+0.09699681497757856*AOrd[7]+0.06062300936098657*AOrd[6]+0.09699681497757856*AOrd[5]+0.1551949039641257*AOrd[4]+0.09699681497757856*AOrd[3]+0.06062300936098657*AOrd[2]+0.09699681497757856*AOrd[1]+0.06062300936098657*AOrd[0]; 
  ASqrt[1] = 0.08133430195906327*AOrd[26]+0.1301348831345013*AOrd[25]+0.08133430195906327*AOrd[24]+0.1301348831345013*AOrd[23]+0.2082158130152021*AOrd[22]+0.1301348831345013*AOrd[21]+0.08133430195906327*AOrd[20]+0.1301348831345013*AOrd[19]+0.08133430195906327*AOrd[18]-0.08133430195906327*AOrd[8]-0.1301348831345013*AOrd[7]-0.08133430195906327*AOrd[6]-0.1301348831345013*AOrd[5]-0.2082158130152021*AOrd[4]-0.1301348831345013*AOrd[3]-0.08133430195906327*AOrd[2]-0.1301348831345013*AOrd[1]-0.08133430195906327*AOrd[0]; 
  ASqrt[2] = 0.08133430195906327*AOrd[26]+0.1301348831345013*AOrd[25]+0.08133430195906327*AOrd[24]-0.08133430195906327*AOrd[20]-0.1301348831345013*AOrd[19]-0.08133430195906327*AOrd[18]+0.1301348831345013*AOrd[17]+0.2082158130152021*AOrd[16]+0.1301348831345013*AOrd[15]-0.1301348831345013*AOrd[11]-0.2082158130152021*AOrd[10]-0.1301348831345013*AOrd[9]+0.08133430195906327*AOrd[8]+0.1301348831345013*AOrd[7]+0.08133430195906327*AOrd[6]-0.08133430195906327*AOrd[2]-0.1301348831345013*AOrd[1]-0.08133430195906327*AOrd[0]; 
  ASqrt[3] = 0.08133430195906327*AOrd[26]-0.08133430195906327*AOrd[24]+0.1301348831345013*AOrd[23]-0.1301348831345013*AOrd[21]+0.08133430195906327*AOrd[20]-0.08133430195906327*AOrd[18]+0.1301348831345013*AOrd[17]-0.1301348831345013*AOrd[15]+0.2082158130152021*AOrd[14]-0.2082158130152021*AOrd[12]+0.1301348831345013*AOrd[11]-0.1301348831345013*AOrd[9]+0.08133430195906327*AOrd[8]-0.08133430195906327*AOrd[6]+0.1301348831345013*AOrd[5]-0.1301348831345013*AOrd[3]+0.08133430195906327*AOrd[2]-0.08133430195906327*AOrd[0]; 
  ASqrt[4] = 0.1091214168497758*AOrd[26]+0.1745942669596414*AOrd[25]+0.1091214168497758*AOrd[24]-0.1091214168497758*AOrd[20]-0.1745942669596414*AOrd[19]-0.1091214168497758*AOrd[18]-0.1091214168497758*AOrd[8]-0.1745942669596414*AOrd[7]-0.1091214168497758*AOrd[6]+0.1091214168497758*AOrd[2]+0.1745942669596414*AOrd[1]+0.1091214168497758*AOrd[0]; 
  ASqrt[5] = 0.1091214168497758*AOrd[26]-0.1091214168497758*AOrd[24]+0.1745942669596414*AOrd[23]-0.1745942669596414*AOrd[21]+0.1091214168497758*AOrd[20]-0.1091214168497758*AOrd[18]-0.1091214168497758*AOrd[8]+0.1091214168497758*AOrd[6]-0.1745942669596414*AOrd[5]+0.1745942669596414*AOrd[3]-0.1091214168497758*AOrd[2]+0.1091214168497758*AOrd[0]; 
  ASqrt[6] = 0.1091214168497758*AOrd[26]-0.1091214168497758*AOrd[24]-0.1091214168497758*AOrd[20]+0.1091214168497758*AOrd[18]+0.1745942669596414*AOrd[17]-0.1745942669596414*AOrd[15]-0.1745942669596414*AOrd[11]+0.1745942669596414*AOrd[9]+0.1091214168497758*AOrd[8]-0.1091214168497758*AOrd[6]-0.1091214168497758*AOrd[2]+0.1091214168497758*AOrd[0]; 
  ASqrt[7] = 0.05422286797270884*AOrd[26]+0.08675658875633419*AOrd[25]+0.05422286797270884*AOrd[24]+0.08675658875633419*AOrd[23]+0.1388105420101347*AOrd[22]+0.08675658875633419*AOrd[21]+0.05422286797270884*AOrd[20]+0.08675658875633419*AOrd[19]+0.05422286797270884*AOrd[18]-0.1084457359454177*AOrd[17]-0.1735131775126684*AOrd[16]-0.1084457359454177*AOrd[15]-0.1735131775126684*AOrd[14]-0.2776210840202695*AOrd[13]-0.1735131775126684*AOrd[12]-0.1084457359454177*AOrd[11]-0.1735131775126684*AOrd[10]-0.1084457359454177*AOrd[9]+0.05422286797270884*AOrd[8]+0.08675658875633419*AOrd[7]+0.05422286797270884*AOrd[6]+0.08675658875633419*AOrd[5]+0.1388105420101347*AOrd[4]+0.08675658875633419*AOrd[3]+0.05422286797270884*AOrd[2]+0.08675658875633419*AOrd[1]+0.05422286797270884*AOrd[0]; 
  ASqrt[8] = 0.05422286797270884*AOrd[26]+0.08675658875633419*AOrd[25]+0.05422286797270884*AOrd[24]-0.1084457359454177*AOrd[23]-0.1735131775126684*AOrd[22]-0.1084457359454177*AOrd[21]+0.05422286797270884*AOrd[20]+0.08675658875633419*AOrd[19]+0.05422286797270884*AOrd[18]+0.08675658875633419*AOrd[17]+0.1388105420101347*AOrd[16]+0.08675658875633419*AOrd[15]-0.1735131775126684*AOrd[14]-0.2776210840202695*AOrd[13]-0.1735131775126684*AOrd[12]+0.08675658875633419*AOrd[11]+0.1388105420101347*AOrd[10]+0.08675658875633419*AOrd[9]+0.05422286797270884*AOrd[8]+0.08675658875633419*AOrd[7]+0.05422286797270884*AOrd[6]-0.1084457359454177*AOrd[5]-0.1735131775126684*AOrd[4]-0.1084457359454177*AOrd[3]+0.05422286797270884*AOrd[2]+0.08675658875633419*AOrd[1]+0.05422286797270884*AOrd[0]; 
  ASqrt[9] = 0.05422286797270884*AOrd[26]-0.1084457359454177*AOrd[25]+0.05422286797270884*AOrd[24]+0.08675658875633419*AOrd[23]-0.1735131775126684*AOrd[22]+0.08675658875633419*AOrd[21]+0.05422286797270884*AOrd[20]-0.1084457359454177*AOrd[19]+0.05422286797270884*AOrd[18]+0.08675658875633419*AOrd[17]-0.1735131775126684*AOrd[16]+0.08675658875633419*AOrd[15]+0.1388105420101347*AOrd[14]-0.2776210840202695*AOrd[13]+0.1388105420101347*AOrd[12]+0.08675658875633419*AOrd[11]-0.1735131775126684*AOrd[10]+0.08675658875633419*AOrd[9]+0.05422286797270884*AOrd[8]-0.1084457359454177*AOrd[7]+0.05422286797270884*AOrd[6]+0.08675658875633419*AOrd[5]-0.1735131775126684*AOrd[4]+0.08675658875633419*AOrd[3]+0.05422286797270884*AOrd[2]-0.1084457359454177*AOrd[1]+0.05422286797270884*AOrd[0]; 
  ASqrt[10] = 0.1464017435263139*AOrd[26]-0.1464017435263139*AOrd[24]-0.1464017435263139*AOrd[20]+0.1464017435263139*AOrd[18]-0.1464017435263139*AOrd[8]+0.1464017435263139*AOrd[6]+0.1464017435263139*AOrd[2]-0.1464017435263139*AOrd[0]; 
  ASqrt[11] = 0.07274761123318395*AOrd[26]+0.1163961779730944*AOrd[25]+0.07274761123318395*AOrd[24]-0.07274761123318395*AOrd[20]-0.1163961779730944*AOrd[19]-0.07274761123318395*AOrd[18]-0.1454952224663679*AOrd[17]-0.2327923559461888*AOrd[16]-0.1454952224663679*AOrd[15]+0.1454952224663679*AOrd[11]+0.2327923559461888*AOrd[10]+0.1454952224663679*AOrd[9]+0.07274761123318395*AOrd[8]+0.1163961779730944*AOrd[7]+0.07274761123318395*AOrd[6]-0.07274761123318395*AOrd[2]-0.1163961779730944*AOrd[1]-0.07274761123318395*AOrd[0]; 
  ASqrt[12] = 0.07274761123318395*AOrd[26]+0.1163961779730944*AOrd[25]+0.07274761123318395*AOrd[24]-0.1454952224663679*AOrd[23]-0.2327923559461888*AOrd[22]-0.1454952224663679*AOrd[21]+0.07274761123318395*AOrd[20]+0.1163961779730944*AOrd[19]+0.07274761123318395*AOrd[18]-0.07274761123318395*AOrd[8]-0.1163961779730944*AOrd[7]-0.07274761123318395*AOrd[6]+0.1454952224663679*AOrd[5]+0.2327923559461888*AOrd[4]+0.1454952224663679*AOrd[3]-0.07274761123318395*AOrd[2]-0.1163961779730944*AOrd[1]-0.07274761123318395*AOrd[0]; 
  ASqrt[13] = 0.07274761123318395*AOrd[26]-0.07274761123318395*AOrd[24]+0.1163961779730944*AOrd[23]-0.1163961779730944*AOrd[21]+0.07274761123318395*AOrd[20]-0.07274761123318395*AOrd[18]-0.1454952224663679*AOrd[17]+0.1454952224663679*AOrd[15]-0.2327923559461888*AOrd[14]+0.2327923559461888*AOrd[12]-0.1454952224663679*AOrd[11]+0.1454952224663679*AOrd[9]+0.07274761123318395*AOrd[8]-0.07274761123318395*AOrd[6]+0.1163961779730944*AOrd[5]-0.1163961779730944*AOrd[3]+0.07274761123318395*AOrd[2]-0.07274761123318395*AOrd[0]; 
  ASqrt[14] = 0.07274761123318395*AOrd[26]-0.07274761123318395*AOrd[24]-0.1454952224663679*AOrd[23]+0.1454952224663679*AOrd[21]+0.07274761123318395*AOrd[20]-0.07274761123318395*AOrd[18]+0.1163961779730944*AOrd[17]-0.1163961779730944*AOrd[15]-0.2327923559461888*AOrd[14]+0.2327923559461888*AOrd[12]+0.1163961779730944*AOrd[11]-0.1163961779730944*AOrd[9]+0.07274761123318395*AOrd[8]-0.07274761123318395*AOrd[6]-0.1454952224663679*AOrd[5]+0.1454952224663679*AOrd[3]+0.07274761123318395*AOrd[2]-0.07274761123318395*AOrd[0]; 
  ASqrt[15] = 0.07274761123318395*AOrd[26]-0.1454952224663679*AOrd[25]+0.07274761123318395*AOrd[24]+0.1163961779730944*AOrd[23]-0.2327923559461888*AOrd[22]+0.1163961779730944*AOrd[21]+0.07274761123318395*AOrd[20]-0.1454952224663679*AOrd[19]+0.07274761123318395*AOrd[18]-0.07274761123318395*AOrd[8]+0.1454952224663679*AOrd[7]-0.07274761123318395*AOrd[6]-0.1163961779730944*AOrd[5]+0.2327923559461888*AOrd[4]-0.1163961779730944*AOrd[3]-0.07274761123318395*AOrd[2]+0.1454952224663679*AOrd[1]-0.07274761123318395*AOrd[0]; 
  ASqrt[16] = 0.07274761123318395*AOrd[26]-0.1454952224663679*AOrd[25]+0.07274761123318395*AOrd[24]-0.07274761123318395*AOrd[20]+0.1454952224663679*AOrd[19]-0.07274761123318395*AOrd[18]+0.1163961779730944*AOrd[17]-0.2327923559461888*AOrd[16]+0.1163961779730944*AOrd[15]-0.1163961779730944*AOrd[11]+0.2327923559461888*AOrd[10]-0.1163961779730944*AOrd[9]+0.07274761123318395*AOrd[8]-0.1454952224663679*AOrd[7]+0.07274761123318395*AOrd[6]-0.07274761123318395*AOrd[2]+0.1454952224663679*AOrd[1]-0.07274761123318395*AOrd[0]; 
  ASqrt[17] = 0.09760116235087592*AOrd[26]-0.09760116235087592*AOrd[24]-0.09760116235087592*AOrd[20]+0.09760116235087592*AOrd[18]-0.1952023247017519*AOrd[17]+0.1952023247017519*AOrd[15]+0.1952023247017519*AOrd[11]-0.1952023247017519*AOrd[9]+0.09760116235087592*AOrd[8]-0.09760116235087592*AOrd[6]-0.09760116235087592*AOrd[2]+0.09760116235087592*AOrd[0]; 
  ASqrt[18] = 0.09760116235087592*AOrd[26]-0.09760116235087592*AOrd[24]-0.1952023247017519*AOrd[23]+0.1952023247017519*AOrd[21]+0.09760116235087592*AOrd[20]-0.09760116235087592*AOrd[18]-0.09760116235087592*AOrd[8]+0.09760116235087592*AOrd[6]+0.1952023247017519*AOrd[5]-0.1952023247017519*AOrd[3]-0.09760116235087592*AOrd[2]+0.09760116235087592*AOrd[0]; 
  ASqrt[19] = 0.09760116235087592*AOrd[26]-0.1952023247017519*AOrd[25]+0.09760116235087592*AOrd[24]-0.09760116235087592*AOrd[20]+0.1952023247017519*AOrd[19]-0.09760116235087592*AOrd[18]-0.09760116235087592*AOrd[8]+0.1952023247017519*AOrd[7]-0.09760116235087592*AOrd[6]+0.09760116235087592*AOrd[2]-0.1952023247017519*AOrd[1]+0.09760116235087592*AOrd[0]; 
  ASqrt[20] = 0.04849840748878927*AOrd[26]+0.07759745198206287*AOrd[25]+0.04849840748878927*AOrd[24]-0.09699681497757856*AOrd[23]-0.1551949039641257*AOrd[22]-0.09699681497757856*AOrd[21]+0.04849840748878927*AOrd[20]+0.07759745198206287*AOrd[19]+0.04849840748878927*AOrd[18]-0.09699681497757856*AOrd[17]-0.1551949039641257*AOrd[16]-0.09699681497757856*AOrd[15]+0.1939936299551571*AOrd[14]+0.3103898079282515*AOrd[13]+0.1939936299551571*AOrd[12]-0.09699681497757856*AOrd[11]-0.1551949039641257*AOrd[10]-0.09699681497757856*AOrd[9]+0.04849840748878927*AOrd[8]+0.07759745198206287*AOrd[7]+0.04849840748878927*AOrd[6]-0.09699681497757856*AOrd[5]-0.1551949039641257*AOrd[4]-0.09699681497757856*AOrd[3]+0.04849840748878927*AOrd[2]+0.07759745198206287*AOrd[1]+0.04849840748878927*AOrd[0]; 
  ASqrt[21] = 0.04849840748878927*AOrd[26]-0.09699681497757856*AOrd[25]+0.04849840748878927*AOrd[24]+0.07759745198206287*AOrd[23]-0.1551949039641257*AOrd[22]+0.07759745198206287*AOrd[21]+0.04849840748878927*AOrd[20]-0.09699681497757856*AOrd[19]+0.04849840748878927*AOrd[18]-0.09699681497757856*AOrd[17]+0.1939936299551571*AOrd[16]-0.09699681497757856*AOrd[15]-0.1551949039641257*AOrd[14]+0.3103898079282515*AOrd[13]-0.1551949039641257*AOrd[12]-0.09699681497757856*AOrd[11]+0.1939936299551571*AOrd[10]-0.09699681497757856*AOrd[9]+0.04849840748878927*AOrd[8]-0.09699681497757856*AOrd[7]+0.04849840748878927*AOrd[6]+0.07759745198206287*AOrd[5]-0.1551949039641257*AOrd[4]+0.07759745198206287*AOrd[3]+0.04849840748878927*AOrd[2]-0.09699681497757856*AOrd[1]+0.04849840748878927*AOrd[0]; 
  ASqrt[22] = 0.04849840748878927*AOrd[26]-0.09699681497757856*AOrd[25]+0.04849840748878927*AOrd[24]-0.09699681497757856*AOrd[23]+0.1939936299551571*AOrd[22]-0.09699681497757856*AOrd[21]+0.04849840748878927*AOrd[20]-0.09699681497757856*AOrd[19]+0.04849840748878927*AOrd[18]+0.07759745198206287*AOrd[17]-0.1551949039641257*AOrd[16]+0.07759745198206287*AOrd[15]-0.1551949039641257*AOrd[14]+0.3103898079282515*AOrd[13]-0.1551949039641257*AOrd[12]+0.07759745198206287*AOrd[11]-0.1551949039641257*AOrd[10]+0.07759745198206287*AOrd[9]+0.04849840748878927*AOrd[8]-0.09699681497757856*AOrd[7]+0.04849840748878927*AOrd[6]-0.09699681497757856*AOrd[5]+0.1939936299551571*AOrd[4]-0.09699681497757856*AOrd[3]+0.04849840748878927*AOrd[2]-0.09699681497757856*AOrd[1]+0.04849840748878927*AOrd[0]; 
  ASqrt[23] = 0.06506744156725063*AOrd[26]-0.06506744156725063*AOrd[24]-0.1301348831345013*AOrd[23]+0.1301348831345013*AOrd[21]+0.06506744156725063*AOrd[20]-0.06506744156725063*AOrd[18]-0.1301348831345013*AOrd[17]+0.1301348831345013*AOrd[15]+0.2602697662690026*AOrd[14]-0.2602697662690026*AOrd[12]-0.1301348831345013*AOrd[11]+0.1301348831345013*AOrd[9]+0.06506744156725063*AOrd[8]-0.06506744156725063*AOrd[6]-0.1301348831345013*AOrd[5]+0.1301348831345013*AOrd[3]+0.06506744156725063*AOrd[2]-0.06506744156725063*AOrd[0]; 
  ASqrt[24] = 0.06506744156725063*AOrd[26]-0.1301348831345013*AOrd[25]+0.06506744156725063*AOrd[24]-0.06506744156725063*AOrd[20]+0.1301348831345013*AOrd[19]-0.06506744156725063*AOrd[18]-0.1301348831345013*AOrd[17]+0.2602697662690026*AOrd[16]-0.1301348831345013*AOrd[15]+0.1301348831345013*AOrd[11]-0.2602697662690026*AOrd[10]+0.1301348831345013*AOrd[9]+0.06506744156725063*AOrd[8]-0.1301348831345013*AOrd[7]+0.06506744156725063*AOrd[6]-0.06506744156725063*AOrd[2]+0.1301348831345013*AOrd[1]-0.06506744156725063*AOrd[0]; 
  ASqrt[25] = 0.06506744156725063*AOrd[26]-0.1301348831345013*AOrd[25]+0.06506744156725063*AOrd[24]-0.1301348831345013*AOrd[23]+0.2602697662690026*AOrd[22]-0.1301348831345013*AOrd[21]+0.06506744156725063*AOrd[20]-0.1301348831345013*AOrd[19]+0.06506744156725063*AOrd[18]-0.06506744156725063*AOrd[8]+0.1301348831345013*AOrd[7]-0.06506744156725063*AOrd[6]+0.1301348831345013*AOrd[5]-0.2602697662690026*AOrd[4]+0.1301348831345013*AOrd[3]-0.06506744156725063*AOrd[2]+0.1301348831345013*AOrd[1]-0.06506744156725063*AOrd[0]; 
  ASqrt[26] = 0.04337829437816709*AOrd[26]-0.08675658875633419*AOrd[25]+0.04337829437816709*AOrd[24]-0.08675658875633419*AOrd[23]+0.1735131775126684*AOrd[22]-0.08675658875633419*AOrd[21]+0.04337829437816709*AOrd[20]-0.08675658875633419*AOrd[19]+0.04337829437816709*AOrd[18]-0.08675658875633419*AOrd[17]+0.1735131775126684*AOrd[16]-0.08675658875633419*AOrd[15]+0.1735131775126684*AOrd[14]-0.3470263550253368*AOrd[13]+0.1735131775126684*AOrd[12]-0.08675658875633419*AOrd[11]+0.1735131775126684*AOrd[10]-0.08675658875633419*AOrd[9]+0.04337829437816709*AOrd[8]-0.08675658875633419*AOrd[7]+0.04337829437816709*AOrd[6]-0.08675658875633419*AOrd[5]+0.1735131775126684*AOrd[4]-0.08675658875633419*AOrd[3]+0.04337829437816709*AOrd[2]-0.08675658875633419*AOrd[1]+0.04337829437816709*AOrd[0]; 

} 
 
