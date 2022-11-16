#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH void lbo_vlasov_diff_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[3]: cell-center coordinates. 
  // dxv[3]: cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[9]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fl/fc/fr: distribution function in cells 
  // out: incremented distribution function in cell 

  const double *nuVtSqSum = &nuPrimMomsSum[6];

  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 
  double incr[20]; 

  incr[0] = 0.4743416490252568*nuVtSqSum[1]*fr[12]+0.4743416490252568*nuVtSqSum[1]*fl[12]-0.9486832980505138*nuVtSqSum[1]*fc[12]-0.842012099081717*nuVtSqSum[2]*fr[11]+0.842012099081717*nuVtSqSum[2]*fl[11]+0.4743416490252568*nuVtSqSum[0]*fr[8]+0.4743416490252568*nuVtSqSum[0]*fl[8]-0.9486832980505137*nuVtSqSum[0]*fc[8]+0.6629126073623879*nuVtSqSum[2]*fr[7]+0.6629126073623879*nuVtSqSum[2]*fl[7]-1.325825214724776*nuVtSqSum[2]*fc[7]-0.8420120990817169*nuVtSqSum[1]*fr[4]+0.8420120990817169*nuVtSqSum[1]*fl[4]-0.8420120990817169*nuVtSqSum[0]*fr[2]+0.8420120990817169*nuVtSqSum[0]*fl[2]+0.6629126073623879*fr[1]*nuVtSqSum[1]+0.6629126073623879*fl[1]*nuVtSqSum[1]-1.325825214724776*fc[1]*nuVtSqSum[1]+0.6629126073623879*fr[0]*nuVtSqSum[0]+0.6629126073623879*fl[0]*nuVtSqSum[0]-1.325825214724776*fc[0]*nuVtSqSum[0]; 
  incr[1] = 0.4242640687119285*nuVtSqSum[2]*fr[12]+0.4743416490252568*nuVtSqSum[0]*fr[12]+0.4242640687119285*nuVtSqSum[2]*fl[12]+0.4743416490252568*nuVtSqSum[0]*fl[12]-0.848528137423857*nuVtSqSum[2]*fc[12]-0.9486832980505138*nuVtSqSum[0]*fc[12]-0.7531185165696033*nuVtSqSum[1]*fr[11]+0.7531185165696033*nuVtSqSum[1]*fl[11]+0.4743416490252568*nuVtSqSum[1]*fr[8]+0.4743416490252568*nuVtSqSum[1]*fl[8]-0.9486832980505137*nuVtSqSum[1]*fc[8]+0.592927061281571*nuVtSqSum[1]*fr[7]+0.592927061281571*nuVtSqSum[1]*fl[7]-1.185854122563142*nuVtSqSum[1]*fc[7]-0.753118516569603*nuVtSqSum[2]*fr[4]-0.8420120990817169*nuVtSqSum[0]*fr[4]+0.753118516569603*nuVtSqSum[2]*fl[4]+0.8420120990817169*nuVtSqSum[0]*fl[4]+0.592927061281571*fr[1]*nuVtSqSum[2]+0.592927061281571*fl[1]*nuVtSqSum[2]-1.185854122563142*fc[1]*nuVtSqSum[2]-0.8420120990817169*nuVtSqSum[1]*fr[2]+0.8420120990817169*nuVtSqSum[1]*fl[2]+0.6629126073623879*fr[0]*nuVtSqSum[1]+0.6629126073623879*fl[0]*nuVtSqSum[1]-1.325825214724776*fc[0]*nuVtSqSum[1]+0.6629126073623879*nuVtSqSum[0]*fr[1]+0.6629126073623879*nuVtSqSum[0]*fl[1]-1.325825214724776*nuVtSqSum[0]*fc[1]; 
  incr[2] = 0.522048062622111*nuVtSqSum[1]*fr[12]-0.522048062622111*nuVtSqSum[1]*fl[12]-1.027514541411701*nuVtSqSum[2]*fr[11]-1.027514541411701*nuVtSqSum[2]*fl[11]-3.778601861965611*nuVtSqSum[2]*fc[11]+0.522048062622111*nuVtSqSum[0]*fr[8]-0.522048062622111*nuVtSqSum[0]*fl[8]+0.8420120990817169*nuVtSqSum[2]*fr[7]-0.8420120990817169*nuVtSqSum[2]*fl[7]-1.027514541411701*nuVtSqSum[1]*fr[4]-1.027514541411701*nuVtSqSum[1]*fl[4]-3.778601861965611*nuVtSqSum[1]*fc[4]-1.027514541411701*nuVtSqSum[0]*fr[2]-1.027514541411701*nuVtSqSum[0]*fl[2]-3.778601861965611*nuVtSqSum[0]*fc[2]+0.8420120990817169*fr[1]*nuVtSqSum[1]-0.8420120990817169*fl[1]*nuVtSqSum[1]+0.8420120990817169*fr[0]*nuVtSqSum[0]-0.8420120990817169*fl[0]*nuVtSqSum[0]; 
  incr[3] = 0.4743416490252568*nuVtSqSum[1]*fr[18]+0.4743416490252568*nuVtSqSum[1]*fl[18]-0.9486832980505137*nuVtSqSum[1]*fc[18]-0.8420120990817169*nuVtSqSum[2]*fr[17]+0.8420120990817169*nuVtSqSum[2]*fl[17]+0.4743416490252568*nuVtSqSum[0]*fr[14]+0.4743416490252568*nuVtSqSum[0]*fl[14]-0.9486832980505138*nuVtSqSum[0]*fc[14]+0.6629126073623879*nuVtSqSum[2]*fr[13]+0.6629126073623879*nuVtSqSum[2]*fl[13]-1.325825214724776*nuVtSqSum[2]*fc[13]-0.8420120990817169*nuVtSqSum[1]*fr[10]+0.8420120990817169*nuVtSqSum[1]*fl[10]-0.8420120990817169*nuVtSqSum[0]*fr[6]+0.8420120990817169*nuVtSqSum[0]*fl[6]+0.6629126073623879*nuVtSqSum[1]*fr[5]+0.6629126073623879*nuVtSqSum[1]*fl[5]-1.325825214724776*nuVtSqSum[1]*fc[5]+0.6629126073623879*nuVtSqSum[0]*fr[3]+0.6629126073623879*nuVtSqSum[0]*fl[3]-1.325825214724776*nuVtSqSum[0]*fc[3]; 
  incr[4] = 0.466933982218043*nuVtSqSum[2]*fr[12]+0.522048062622111*nuVtSqSum[0]*fr[12]-0.466933982218043*nuVtSqSum[2]*fl[12]-0.522048062622111*nuVtSqSum[0]*fl[12]-0.9190369449864345*nuVtSqSum[1]*fr[11]-0.9190369449864345*nuVtSqSum[1]*fl[11]-3.379684249304953*nuVtSqSum[1]*fc[11]+0.522048062622111*nuVtSqSum[1]*fr[8]-0.522048062622111*nuVtSqSum[1]*fl[8]+0.753118516569603*nuVtSqSum[1]*fr[7]-0.753118516569603*nuVtSqSum[1]*fl[7]-0.9190369449864345*nuVtSqSum[2]*fr[4]-1.027514541411701*nuVtSqSum[0]*fr[4]-0.9190369449864345*nuVtSqSum[2]*fl[4]-1.027514541411701*nuVtSqSum[0]*fl[4]-3.379684249304954*nuVtSqSum[2]*fc[4]-3.778601861965611*nuVtSqSum[0]*fc[4]+0.753118516569603*fr[1]*nuVtSqSum[2]-0.753118516569603*fl[1]*nuVtSqSum[2]-1.027514541411701*nuVtSqSum[1]*fr[2]-1.027514541411701*nuVtSqSum[1]*fl[2]-3.778601861965611*nuVtSqSum[1]*fc[2]+0.8420120990817169*fr[0]*nuVtSqSum[1]-0.8420120990817169*fl[0]*nuVtSqSum[1]+0.8420120990817169*nuVtSqSum[0]*fr[1]-0.8420120990817169*nuVtSqSum[0]*fl[1]; 
  incr[5] = 0.4242640687119285*nuVtSqSum[2]*fr[18]+0.4743416490252568*nuVtSqSum[0]*fr[18]+0.4242640687119285*nuVtSqSum[2]*fl[18]+0.4743416490252568*nuVtSqSum[0]*fl[18]-0.848528137423857*nuVtSqSum[2]*fc[18]-0.9486832980505137*nuVtSqSum[0]*fc[18]-0.753118516569603*nuVtSqSum[1]*fr[17]+0.753118516569603*nuVtSqSum[1]*fl[17]+0.4743416490252568*nuVtSqSum[1]*fr[14]+0.4743416490252568*nuVtSqSum[1]*fl[14]-0.9486832980505138*nuVtSqSum[1]*fc[14]+0.5929270612815709*nuVtSqSum[1]*fr[13]+0.5929270612815709*nuVtSqSum[1]*fl[13]-1.185854122563142*nuVtSqSum[1]*fc[13]-0.753118516569603*nuVtSqSum[2]*fr[10]-0.8420120990817169*nuVtSqSum[0]*fr[10]+0.753118516569603*nuVtSqSum[2]*fl[10]+0.8420120990817169*nuVtSqSum[0]*fl[10]-0.8420120990817169*nuVtSqSum[1]*fr[6]+0.8420120990817169*nuVtSqSum[1]*fl[6]+0.592927061281571*nuVtSqSum[2]*fr[5]+0.6629126073623879*nuVtSqSum[0]*fr[5]+0.592927061281571*nuVtSqSum[2]*fl[5]+0.6629126073623879*nuVtSqSum[0]*fl[5]-1.185854122563142*nuVtSqSum[2]*fc[5]-1.325825214724776*nuVtSqSum[0]*fc[5]+0.6629126073623879*nuVtSqSum[1]*fr[3]+0.6629126073623879*nuVtSqSum[1]*fl[3]-1.325825214724776*nuVtSqSum[1]*fc[3]; 
  incr[6] = 0.522048062622111*nuVtSqSum[1]*fr[18]-0.522048062622111*nuVtSqSum[1]*fl[18]-1.027514541411701*nuVtSqSum[2]*fr[17]-1.027514541411701*nuVtSqSum[2]*fl[17]-3.778601861965611*nuVtSqSum[2]*fc[17]+0.522048062622111*nuVtSqSum[0]*fr[14]-0.522048062622111*nuVtSqSum[0]*fl[14]+0.842012099081717*nuVtSqSum[2]*fr[13]-0.842012099081717*nuVtSqSum[2]*fl[13]-1.027514541411701*nuVtSqSum[1]*fr[10]-1.027514541411701*nuVtSqSum[1]*fl[10]-3.778601861965611*nuVtSqSum[1]*fc[10]-1.027514541411701*nuVtSqSum[0]*fr[6]-1.027514541411701*nuVtSqSum[0]*fl[6]-3.778601861965611*nuVtSqSum[0]*fc[6]+0.8420120990817169*nuVtSqSum[1]*fr[5]-0.8420120990817169*nuVtSqSum[1]*fl[5]+0.8420120990817169*nuVtSqSum[0]*fr[3]-0.8420120990817169*nuVtSqSum[0]*fl[3]; 
  incr[7] = 0.4242640687119285*nuVtSqSum[1]*fr[12]+0.4242640687119285*nuVtSqSum[1]*fl[12]-0.848528137423857*nuVtSqSum[1]*fc[12]-0.5379417975497165*nuVtSqSum[2]*fr[11]-0.842012099081717*nuVtSqSum[0]*fr[11]+0.5379417975497165*nuVtSqSum[2]*fl[11]+0.842012099081717*nuVtSqSum[0]*fl[11]+0.4743416490252568*nuVtSqSum[2]*fr[8]+0.4743416490252568*nuVtSqSum[2]*fl[8]-0.9486832980505137*nuVtSqSum[2]*fc[8]+0.4235193294868364*nuVtSqSum[2]*fr[7]+0.6629126073623879*nuVtSqSum[0]*fr[7]+0.4235193294868364*nuVtSqSum[2]*fl[7]+0.6629126073623879*nuVtSqSum[0]*fl[7]-0.8470386589736728*nuVtSqSum[2]*fc[7]-1.325825214724776*nuVtSqSum[0]*fc[7]-0.753118516569603*nuVtSqSum[1]*fr[4]+0.753118516569603*nuVtSqSum[1]*fl[4]-0.8420120990817169*fr[2]*nuVtSqSum[2]+0.8420120990817169*fl[2]*nuVtSqSum[2]+0.6629126073623879*fr[0]*nuVtSqSum[2]+0.6629126073623879*fl[0]*nuVtSqSum[2]-1.325825214724776*fc[0]*nuVtSqSum[2]+0.592927061281571*fr[1]*nuVtSqSum[1]+0.592927061281571*fl[1]*nuVtSqSum[1]-1.185854122563142*fc[1]*nuVtSqSum[1]; 
  incr[8] = (-0.09943689110435815*nuVtSqSum[1]*fr[12])-0.09943689110435815*nuVtSqSum[1]*fl[12]-4.441514469327998*nuVtSqSum[1]*fc[12]-0.2139541240254553*nuVtSqSum[2]*fr[11]+0.2139541240254553*nuVtSqSum[2]*fl[11]-0.09943689110435816*nuVtSqSum[0]*fr[8]-0.09943689110435816*nuVtSqSum[0]*fl[8]-4.441514469327998*nuVtSqSum[0]*fc[8]+0.2964635306407854*nuVtSqSum[2]*fr[7]+0.2964635306407854*nuVtSqSum[2]*fl[7]-0.592927061281571*nuVtSqSum[2]*fc[7]-0.2139541240254554*nuVtSqSum[1]*fr[4]+0.2139541240254554*nuVtSqSum[1]*fl[4]-0.2139541240254554*nuVtSqSum[0]*fr[2]+0.2139541240254554*nuVtSqSum[0]*fl[2]+0.2964635306407854*fr[1]*nuVtSqSum[1]+0.2964635306407854*fl[1]*nuVtSqSum[1]-0.592927061281571*fc[1]*nuVtSqSum[1]+0.2964635306407854*fr[0]*nuVtSqSum[0]+0.2964635306407854*fl[0]*nuVtSqSum[0]-0.592927061281571*fc[0]*nuVtSqSum[0]; 
  incr[9] = (-0.8420120990817169*nuVtSqSum[1]*fr[19])+0.8420120990817169*nuVtSqSum[1]*fl[19]-0.842012099081717*nuVtSqSum[0]*fr[16]+0.842012099081717*nuVtSqSum[0]*fl[16]+0.6629126073623879*nuVtSqSum[1]*fr[15]+0.6629126073623879*nuVtSqSum[1]*fl[15]-1.325825214724776*nuVtSqSum[1]*fc[15]+0.6629126073623879*nuVtSqSum[0]*fr[9]+0.6629126073623879*nuVtSqSum[0]*fl[9]-1.325825214724776*nuVtSqSum[0]*fc[9]; 
  incr[10] = 0.4669339822180429*nuVtSqSum[2]*fr[18]+0.522048062622111*nuVtSqSum[0]*fr[18]-0.4669339822180429*nuVtSqSum[2]*fl[18]-0.522048062622111*nuVtSqSum[0]*fl[18]-0.9190369449864345*nuVtSqSum[1]*fr[17]-0.9190369449864345*nuVtSqSum[1]*fl[17]-3.379684249304954*nuVtSqSum[1]*fc[17]+0.522048062622111*nuVtSqSum[1]*fr[14]-0.522048062622111*nuVtSqSum[1]*fl[14]+0.7531185165696033*nuVtSqSum[1]*fr[13]-0.7531185165696033*nuVtSqSum[1]*fl[13]-0.9190369449864345*nuVtSqSum[2]*fr[10]-1.027514541411701*nuVtSqSum[0]*fr[10]-0.9190369449864345*nuVtSqSum[2]*fl[10]-1.027514541411701*nuVtSqSum[0]*fl[10]-3.379684249304954*nuVtSqSum[2]*fc[10]-3.778601861965611*nuVtSqSum[0]*fc[10]-1.027514541411701*nuVtSqSum[1]*fr[6]-1.027514541411701*nuVtSqSum[1]*fl[6]-3.778601861965611*nuVtSqSum[1]*fc[6]+0.753118516569603*nuVtSqSum[2]*fr[5]+0.8420120990817169*nuVtSqSum[0]*fr[5]-0.753118516569603*nuVtSqSum[2]*fl[5]-0.8420120990817169*nuVtSqSum[0]*fl[5]+0.8420120990817169*nuVtSqSum[1]*fr[3]-0.8420120990817169*nuVtSqSum[1]*fl[3]; 
  incr[11] = 0.4669339822180429*nuVtSqSum[1]*fr[12]-0.4669339822180429*nuVtSqSum[1]*fl[12]-0.6564549607045962*nuVtSqSum[2]*fr[11]-1.027514541411701*nuVtSqSum[0]*fr[11]-0.6564549607045962*nuVtSqSum[2]*fl[11]-1.027514541411701*nuVtSqSum[0]*fl[11]-2.414060178074966*nuVtSqSum[2]*fc[11]-3.778601861965611*nuVtSqSum[0]*fc[11]+0.522048062622111*nuVtSqSum[2]*fr[8]-0.522048062622111*nuVtSqSum[2]*fl[8]+0.5379417975497165*nuVtSqSum[2]*fr[7]+0.842012099081717*nuVtSqSum[0]*fr[7]-0.5379417975497165*nuVtSqSum[2]*fl[7]-0.842012099081717*nuVtSqSum[0]*fl[7]-0.9190369449864345*nuVtSqSum[1]*fr[4]-0.9190369449864345*nuVtSqSum[1]*fl[4]-3.379684249304953*nuVtSqSum[1]*fc[4]-1.027514541411701*fr[2]*nuVtSqSum[2]-1.027514541411701*fl[2]*nuVtSqSum[2]-3.778601861965611*fc[2]*nuVtSqSum[2]+0.842012099081717*fr[0]*nuVtSqSum[2]-0.842012099081717*fl[0]*nuVtSqSum[2]+0.7531185165696033*fr[1]*nuVtSqSum[1]-0.7531185165696033*fl[1]*nuVtSqSum[1]; 
  incr[12] = (-0.0889390591922356*nuVtSqSum[2]*fr[12])-0.09943689110435816*nuVtSqSum[0]*fr[12]-0.0889390591922356*nuVtSqSum[2]*fl[12]-0.09943689110435816*nuVtSqSum[0]*fl[12]-3.972611310586524*nuVtSqSum[2]*fc[12]-4.441514469327998*nuVtSqSum[0]*fc[12]-0.1913663861549356*nuVtSqSum[1]*fr[11]+0.1913663861549356*nuVtSqSum[1]*fl[11]-0.09943689110435815*nuVtSqSum[1]*fr[8]-0.09943689110435815*nuVtSqSum[1]*fl[8]-4.441514469327998*nuVtSqSum[1]*fc[8]+0.2651650429449552*nuVtSqSum[1]*fr[7]+0.2651650429449552*nuVtSqSum[1]*fl[7]-0.5303300858899104*nuVtSqSum[1]*fc[7]-0.1913663861549357*nuVtSqSum[2]*fr[4]-0.2139541240254553*nuVtSqSum[0]*fr[4]+0.1913663861549357*nuVtSqSum[2]*fl[4]+0.2139541240254553*nuVtSqSum[0]*fl[4]+0.2651650429449552*fr[1]*nuVtSqSum[2]+0.2651650429449552*fl[1]*nuVtSqSum[2]-0.5303300858899104*fc[1]*nuVtSqSum[2]-0.2139541240254553*nuVtSqSum[1]*fr[2]+0.2139541240254553*nuVtSqSum[1]*fl[2]+0.2964635306407854*fr[0]*nuVtSqSum[1]+0.2964635306407854*fl[0]*nuVtSqSum[1]-0.5929270612815709*fc[0]*nuVtSqSum[1]+0.2964635306407854*nuVtSqSum[0]*fr[1]+0.2964635306407854*nuVtSqSum[0]*fl[1]-0.5929270612815709*nuVtSqSum[0]*fc[1]; 
  incr[13] = 0.4242640687119285*nuVtSqSum[1]*fr[18]+0.4242640687119285*nuVtSqSum[1]*fl[18]-0.848528137423857*nuVtSqSum[1]*fc[18]-0.5379417975497165*nuVtSqSum[2]*fr[17]-0.842012099081717*nuVtSqSum[0]*fr[17]+0.5379417975497165*nuVtSqSum[2]*fl[17]+0.842012099081717*nuVtSqSum[0]*fl[17]+0.4743416490252568*nuVtSqSum[2]*fr[14]+0.4743416490252568*nuVtSqSum[2]*fl[14]-0.9486832980505137*nuVtSqSum[2]*fc[14]+0.4235193294868364*nuVtSqSum[2]*fr[13]+0.6629126073623879*nuVtSqSum[0]*fr[13]+0.4235193294868364*nuVtSqSum[2]*fl[13]+0.6629126073623879*nuVtSqSum[0]*fl[13]-0.8470386589736728*nuVtSqSum[2]*fc[13]-1.325825214724776*nuVtSqSum[0]*fc[13]-0.7531185165696033*nuVtSqSum[1]*fr[10]+0.7531185165696033*nuVtSqSum[1]*fl[10]-0.842012099081717*nuVtSqSum[2]*fr[6]+0.842012099081717*nuVtSqSum[2]*fl[6]+0.5929270612815709*nuVtSqSum[1]*fr[5]+0.5929270612815709*nuVtSqSum[1]*fl[5]-1.185854122563142*nuVtSqSum[1]*fc[5]+0.6629126073623879*nuVtSqSum[2]*fr[3]+0.6629126073623879*nuVtSqSum[2]*fl[3]-1.325825214724776*nuVtSqSum[2]*fc[3]; 
  incr[14] = (-0.09943689110435815*nuVtSqSum[1]*fr[18])-0.09943689110435815*nuVtSqSum[1]*fl[18]-4.441514469327998*nuVtSqSum[1]*fc[18]-0.2139541240254553*nuVtSqSum[2]*fr[17]+0.2139541240254553*nuVtSqSum[2]*fl[17]-0.09943689110435816*nuVtSqSum[0]*fr[14]-0.09943689110435816*nuVtSqSum[0]*fl[14]-4.441514469327998*nuVtSqSum[0]*fc[14]+0.2964635306407854*nuVtSqSum[2]*fr[13]+0.2964635306407854*nuVtSqSum[2]*fl[13]-0.592927061281571*nuVtSqSum[2]*fc[13]-0.2139541240254553*nuVtSqSum[1]*fr[10]+0.2139541240254553*nuVtSqSum[1]*fl[10]-0.2139541240254553*nuVtSqSum[0]*fr[6]+0.2139541240254553*nuVtSqSum[0]*fl[6]+0.2964635306407854*nuVtSqSum[1]*fr[5]+0.2964635306407854*nuVtSqSum[1]*fl[5]-0.5929270612815709*nuVtSqSum[1]*fc[5]+0.2964635306407854*nuVtSqSum[0]*fr[3]+0.2964635306407854*nuVtSqSum[0]*fl[3]-0.5929270612815709*nuVtSqSum[0]*fc[3]; 
  incr[15] = (-0.7531185165696033*nuVtSqSum[2]*fr[19])-0.842012099081717*nuVtSqSum[0]*fr[19]+0.7531185165696033*nuVtSqSum[2]*fl[19]+0.842012099081717*nuVtSqSum[0]*fl[19]-0.8420120990817169*nuVtSqSum[1]*fr[16]+0.8420120990817169*nuVtSqSum[1]*fl[16]+0.592927061281571*nuVtSqSum[2]*fr[15]+0.6629126073623879*nuVtSqSum[0]*fr[15]+0.592927061281571*nuVtSqSum[2]*fl[15]+0.6629126073623879*nuVtSqSum[0]*fl[15]-1.185854122563142*nuVtSqSum[2]*fc[15]-1.325825214724776*nuVtSqSum[0]*fc[15]+0.6629126073623879*nuVtSqSum[1]*fr[9]+0.6629126073623879*nuVtSqSum[1]*fl[9]-1.325825214724776*nuVtSqSum[1]*fc[9]; 
  incr[16] = (-1.027514541411701*nuVtSqSum[1]*fr[19])-1.027514541411701*nuVtSqSum[1]*fl[19]-3.778601861965611*nuVtSqSum[1]*fc[19]-1.027514541411701*nuVtSqSum[0]*fr[16]-1.027514541411701*nuVtSqSum[0]*fl[16]-3.778601861965611*nuVtSqSum[0]*fc[16]+0.8420120990817169*nuVtSqSum[1]*fr[15]-0.8420120990817169*nuVtSqSum[1]*fl[15]+0.842012099081717*nuVtSqSum[0]*fr[9]-0.842012099081717*nuVtSqSum[0]*fl[9]; 
  incr[17] = 0.4669339822180429*nuVtSqSum[1]*fr[18]-0.4669339822180429*nuVtSqSum[1]*fl[18]-0.6564549607045962*nuVtSqSum[2]*fr[17]-1.027514541411701*nuVtSqSum[0]*fr[17]-0.6564549607045962*nuVtSqSum[2]*fl[17]-1.027514541411701*nuVtSqSum[0]*fl[17]-2.414060178074966*nuVtSqSum[2]*fc[17]-3.778601861965611*nuVtSqSum[0]*fc[17]+0.522048062622111*nuVtSqSum[2]*fr[14]-0.522048062622111*nuVtSqSum[2]*fl[14]+0.5379417975497165*nuVtSqSum[2]*fr[13]+0.842012099081717*nuVtSqSum[0]*fr[13]-0.5379417975497165*nuVtSqSum[2]*fl[13]-0.842012099081717*nuVtSqSum[0]*fl[13]-0.9190369449864345*nuVtSqSum[1]*fr[10]-0.9190369449864345*nuVtSqSum[1]*fl[10]-3.379684249304954*nuVtSqSum[1]*fc[10]-1.027514541411701*nuVtSqSum[2]*fr[6]-1.027514541411701*nuVtSqSum[2]*fl[6]-3.778601861965611*nuVtSqSum[2]*fc[6]+0.753118516569603*nuVtSqSum[1]*fr[5]-0.753118516569603*nuVtSqSum[1]*fl[5]+0.8420120990817169*nuVtSqSum[2]*fr[3]-0.8420120990817169*nuVtSqSum[2]*fl[3]; 
  incr[18] = (-0.0889390591922356*nuVtSqSum[2]*fr[18])-0.09943689110435816*nuVtSqSum[0]*fr[18]-0.0889390591922356*nuVtSqSum[2]*fl[18]-0.09943689110435816*nuVtSqSum[0]*fl[18]-3.972611310586524*nuVtSqSum[2]*fc[18]-4.441514469327998*nuVtSqSum[0]*fc[18]-0.1913663861549356*nuVtSqSum[1]*fr[17]+0.1913663861549356*nuVtSqSum[1]*fl[17]-0.09943689110435815*nuVtSqSum[1]*fr[14]-0.09943689110435815*nuVtSqSum[1]*fl[14]-4.441514469327998*nuVtSqSum[1]*fc[14]+0.2651650429449552*nuVtSqSum[1]*fr[13]+0.2651650429449552*nuVtSqSum[1]*fl[13]-0.5303300858899104*nuVtSqSum[1]*fc[13]-0.1913663861549356*nuVtSqSum[2]*fr[10]-0.2139541240254554*nuVtSqSum[0]*fr[10]+0.1913663861549356*nuVtSqSum[2]*fl[10]+0.2139541240254554*nuVtSqSum[0]*fl[10]-0.2139541240254554*nuVtSqSum[1]*fr[6]+0.2139541240254554*nuVtSqSum[1]*fl[6]+0.2651650429449552*nuVtSqSum[2]*fr[5]+0.2964635306407854*nuVtSqSum[0]*fr[5]+0.2651650429449552*nuVtSqSum[2]*fl[5]+0.2964635306407854*nuVtSqSum[0]*fl[5]-0.5303300858899105*nuVtSqSum[2]*fc[5]-0.592927061281571*nuVtSqSum[0]*fc[5]+0.2964635306407854*nuVtSqSum[1]*fr[3]+0.2964635306407854*nuVtSqSum[1]*fl[3]-0.592927061281571*nuVtSqSum[1]*fc[3]; 
  incr[19] = (-0.9190369449864345*nuVtSqSum[2]*fr[19])-1.027514541411701*nuVtSqSum[0]*fr[19]-0.9190369449864345*nuVtSqSum[2]*fl[19]-1.027514541411701*nuVtSqSum[0]*fl[19]-3.379684249304954*nuVtSqSum[2]*fc[19]-3.778601861965611*nuVtSqSum[0]*fc[19]-1.027514541411701*nuVtSqSum[1]*fr[16]-1.027514541411701*nuVtSqSum[1]*fl[16]-3.778601861965611*nuVtSqSum[1]*fc[16]+0.7531185165696033*nuVtSqSum[2]*fr[15]+0.842012099081717*nuVtSqSum[0]*fr[15]-0.7531185165696033*nuVtSqSum[2]*fl[15]-0.842012099081717*nuVtSqSum[0]*fl[15]+0.8420120990817169*nuVtSqSum[1]*fr[9]-0.8420120990817169*nuVtSqSum[1]*fl[9]; 

  out[0] += incr[0]*rdvSq4; 
  out[1] += incr[1]*rdvSq4; 
  out[2] += incr[2]*rdvSq4; 
  out[3] += incr[3]*rdvSq4; 
  out[4] += incr[4]*rdvSq4; 
  out[5] += incr[5]*rdvSq4; 
  out[6] += incr[6]*rdvSq4; 
  out[7] += incr[7]*rdvSq4; 
  out[8] += incr[8]*rdvSq4; 
  out[9] += incr[9]*rdvSq4; 
  out[10] += incr[10]*rdvSq4; 
  out[11] += incr[11]*rdvSq4; 
  out[12] += incr[12]*rdvSq4; 
  out[13] += incr[13]*rdvSq4; 
  out[14] += incr[14]*rdvSq4; 
  out[15] += incr[15]*rdvSq4; 
  out[16] += incr[16]*rdvSq4; 
  out[17] += incr[17]*rdvSq4; 
  out[18] += incr[18]*rdvSq4; 
  out[19] += incr[19]*rdvSq4; 
} 
