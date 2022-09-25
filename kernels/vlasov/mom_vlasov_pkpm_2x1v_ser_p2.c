#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void mom_vlasov_pkpm_2x1v_ser_p2(const double *w, const double *dxv, const int *idx, double mass, const double *bvar, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]/2.0; 
  const double wvpar = w[2], dvpar = dxv[2]; 
  const double wvpar_sq = wvpar*wvpar, dvpar_sq = dvpar*dvpar; 
  const double wvpar_cu = wvpar*wvpar*wvpar, dvpar_cu = dvpar*dvpar*dvpar; 

  const double *bx = &bvar[0]; 
  const double *by = &bvar[8]; 
  const double *bz = &bvar[16]; 

  out[0] += 1.414213562373095*f[0]*mass*volFact; 
  out[1] += 1.414213562373095*f[1]*mass*volFact; 
  out[2] += 1.414213562373095*f[2]*mass*volFact; 
  out[3] += 1.414213562373095*f[4]*mass*volFact; 
  out[4] += 1.414213562373095*f[7]*mass*volFact; 
  out[5] += 1.414213562373095*f[8]*mass*volFact; 
  out[6] += 1.414213562373095*f[11]*mass*volFact; 
  out[7] += 1.414213562373095*f[12]*mass*volFact; 
  out[8] += mass*volFact*(1.414213562373095*f[0]*wvpar_sq+0.8164965809277261*f[3]*dvpar*wvpar+0.105409255338946*f[9]*dvpar_sq+0.1178511301977579*f[0]*dvpar_sq); 
  out[9] += mass*volFact*(1.414213562373095*f[1]*wvpar_sq+0.8164965809277261*f[5]*dvpar*wvpar+0.105409255338946*f[15]*dvpar_sq+0.1178511301977579*f[1]*dvpar_sq); 
  out[10] += mass*volFact*(1.414213562373095*f[2]*wvpar_sq+0.8164965809277261*f[6]*dvpar*wvpar+0.105409255338946*f[16]*dvpar_sq+0.1178511301977579*f[2]*dvpar_sq); 
  out[11] += mass*volFact*(1.414213562373095*f[4]*wvpar_sq+0.8164965809277261*f[10]*dvpar*wvpar+0.105409255338946*f[19]*dvpar_sq+0.1178511301977579*f[4]*dvpar_sq); 
  out[12] += mass*volFact*(1.414213562373095*f[7]*wvpar_sq+0.816496580927726*f[13]*dvpar*wvpar+0.1178511301977579*f[7]*dvpar_sq); 
  out[13] += mass*volFact*(1.414213562373095*f[8]*wvpar_sq+0.816496580927726*f[14]*dvpar*wvpar+0.1178511301977579*f[8]*dvpar_sq); 
  out[14] += mass*volFact*(1.414213562373095*f[11]*wvpar_sq+0.816496580927726*f[17]*dvpar*wvpar+0.1178511301977579*f[11]*dvpar_sq); 
  out[15] += mass*volFact*(1.414213562373095*f[12]*wvpar_sq+0.816496580927726*f[18]*dvpar*wvpar+0.1178511301977579*f[12]*dvpar_sq); 
  out[16] += mass*volFact*(0.6123724356957944*bx[7]*f[18]*dvpar*wvpar_sq+0.6123724356957944*bx[6]*f[17]*dvpar*wvpar_sq+0.6123724356957944*bx[5]*f[14]*dvpar*wvpar_sq+0.6123724356957944*bx[4]*f[13]*dvpar*wvpar_sq+0.6123724356957944*bx[3]*f[10]*dvpar*wvpar_sq+0.6123724356957944*bx[2]*f[6]*dvpar*wvpar_sq+0.6123724356957944*bx[1]*f[5]*dvpar*wvpar_sq+0.6123724356957944*bx[0]*f[3]*dvpar*wvpar_sq+0.7071067811865475*bx[7]*f[12]*wvpar_cu+0.7071067811865475*bx[6]*f[11]*wvpar_cu+0.7071067811865475*bx[5]*f[8]*wvpar_cu+0.7071067811865475*bx[4]*f[7]*wvpar_cu+0.7071067811865475*bx[3]*f[4]*wvpar_cu+0.7071067811865475*bx[2]*f[2]*wvpar_cu+0.7071067811865475*bx[1]*f[1]*wvpar_cu+0.7071067811865475*bx[0]*f[0]*wvpar_cu+0.1581138830084189*bx[3]*f[19]*dvpar_sq*wvpar+0.1581138830084189*bx[2]*f[16]*dvpar_sq*wvpar+0.1581138830084189*bx[1]*f[15]*dvpar_sq*wvpar+0.1767766952966368*bx[7]*f[12]*dvpar_sq*wvpar+0.1767766952966368*bx[6]*f[11]*dvpar_sq*wvpar+0.1581138830084189*bx[0]*f[9]*dvpar_sq*wvpar+0.1767766952966368*bx[5]*f[8]*dvpar_sq*wvpar+0.1767766952966368*bx[4]*f[7]*dvpar_sq*wvpar+0.1767766952966368*bx[3]*f[4]*dvpar_sq*wvpar+0.1767766952966368*bx[2]*f[2]*dvpar_sq*wvpar+0.1767766952966368*bx[1]*f[1]*dvpar_sq*wvpar+0.1767766952966368*bx[0]*f[0]*dvpar_sq*wvpar+0.03061862178478971*bx[7]*f[18]*dvpar_cu+0.03061862178478971*bx[6]*f[17]*dvpar_cu+0.03061862178478971*bx[5]*f[14]*dvpar_cu+0.03061862178478971*bx[4]*f[13]*dvpar_cu+0.03061862178478971*bx[3]*f[10]*dvpar_cu+0.03061862178478971*bx[2]*f[6]*dvpar_cu+0.03061862178478971*bx[1]*f[5]*dvpar_cu+0.03061862178478971*bx[0]*f[3]*dvpar_cu); 
  out[17] += mass*volFact*(0.6123724356957944*bx[5]*f[18]*dvpar*wvpar_sq+0.5477225575051661*bx[3]*f[17]*dvpar*wvpar_sq+0.6123724356957944*bx[7]*f[14]*dvpar*wvpar_sq+0.5477225575051661*bx[1]*f[13]*dvpar*wvpar_sq+0.5477225575051661*bx[6]*f[10]*dvpar*wvpar_sq+0.6123724356957944*bx[2]*f[10]*dvpar*wvpar_sq+0.6123724356957944*bx[3]*f[6]*dvpar*wvpar_sq+0.5477225575051661*bx[4]*f[5]*dvpar*wvpar_sq+0.6123724356957944*bx[0]*f[5]*dvpar*wvpar_sq+0.6123724356957944*bx[1]*f[3]*dvpar*wvpar_sq+0.7071067811865475*bx[5]*f[12]*wvpar_cu+0.632455532033676*bx[3]*f[11]*wvpar_cu+0.7071067811865475*bx[7]*f[8]*wvpar_cu+0.6324555320336759*bx[1]*f[7]*wvpar_cu+0.632455532033676*f[4]*bx[6]*wvpar_cu+0.7071067811865475*bx[2]*f[4]*wvpar_cu+0.6324555320336759*f[1]*bx[4]*wvpar_cu+0.7071067811865475*f[2]*bx[3]*wvpar_cu+0.7071067811865475*bx[0]*f[1]*wvpar_cu+0.7071067811865475*f[0]*bx[1]*wvpar_cu+0.1414213562373095*bx[6]*f[19]*dvpar_sq*wvpar+0.1581138830084189*bx[2]*f[19]*dvpar_sq*wvpar+0.1581138830084189*bx[3]*f[16]*dvpar_sq*wvpar+0.1414213562373095*bx[4]*f[15]*dvpar_sq*wvpar+0.1581138830084189*bx[0]*f[15]*dvpar_sq*wvpar+0.1767766952966368*bx[5]*f[12]*dvpar_sq*wvpar+0.1581138830084189*bx[3]*f[11]*dvpar_sq*wvpar+0.1581138830084189*bx[1]*f[9]*dvpar_sq*wvpar+0.1767766952966368*bx[7]*f[8]*dvpar_sq*wvpar+0.1581138830084189*bx[1]*f[7]*dvpar_sq*wvpar+0.1581138830084189*f[4]*bx[6]*dvpar_sq*wvpar+0.1767766952966368*bx[2]*f[4]*dvpar_sq*wvpar+0.1581138830084189*f[1]*bx[4]*dvpar_sq*wvpar+0.1767766952966368*f[2]*bx[3]*dvpar_sq*wvpar+0.1767766952966368*bx[0]*f[1]*dvpar_sq*wvpar+0.1767766952966368*f[0]*bx[1]*dvpar_sq*wvpar+0.03061862178478971*bx[5]*f[18]*dvpar_cu+0.02738612787525829*bx[3]*f[17]*dvpar_cu+0.03061862178478971*bx[7]*f[14]*dvpar_cu+0.0273861278752583*bx[1]*f[13]*dvpar_cu+0.0273861278752583*bx[6]*f[10]*dvpar_cu+0.03061862178478971*bx[2]*f[10]*dvpar_cu+0.03061862178478971*bx[3]*f[6]*dvpar_cu+0.02738612787525829*bx[4]*f[5]*dvpar_cu+0.03061862178478971*bx[0]*f[5]*dvpar_cu+0.03061862178478971*bx[1]*f[3]*dvpar_cu); 
  out[18] += mass*volFact*(0.5477225575051661*bx[3]*f[18]*dvpar*wvpar_sq+0.6123724356957944*bx[4]*f[17]*dvpar*wvpar_sq+0.5477225575051661*bx[2]*f[14]*dvpar*wvpar_sq+0.6123724356957944*bx[6]*f[13]*dvpar*wvpar_sq+0.5477225575051661*bx[7]*f[10]*dvpar*wvpar_sq+0.6123724356957944*bx[1]*f[10]*dvpar*wvpar_sq+0.5477225575051661*bx[5]*f[6]*dvpar*wvpar_sq+0.6123724356957944*bx[0]*f[6]*dvpar*wvpar_sq+0.6123724356957944*bx[3]*f[5]*dvpar*wvpar_sq+0.6123724356957944*bx[2]*f[3]*dvpar*wvpar_sq+0.632455532033676*bx[3]*f[12]*wvpar_cu+0.7071067811865475*bx[4]*f[11]*wvpar_cu+0.6324555320336759*bx[2]*f[8]*wvpar_cu+0.7071067811865475*bx[6]*f[7]*wvpar_cu+0.632455532033676*f[4]*bx[7]*wvpar_cu+0.6324555320336759*f[2]*bx[5]*wvpar_cu+0.7071067811865475*bx[1]*f[4]*wvpar_cu+0.7071067811865475*f[1]*bx[3]*wvpar_cu+0.7071067811865475*bx[0]*f[2]*wvpar_cu+0.7071067811865475*f[0]*bx[2]*wvpar_cu+0.1414213562373095*bx[7]*f[19]*dvpar_sq*wvpar+0.1581138830084189*bx[1]*f[19]*dvpar_sq*wvpar+0.1414213562373095*bx[5]*f[16]*dvpar_sq*wvpar+0.1581138830084189*bx[0]*f[16]*dvpar_sq*wvpar+0.1581138830084189*bx[3]*f[15]*dvpar_sq*wvpar+0.1581138830084189*bx[3]*f[12]*dvpar_sq*wvpar+0.1767766952966368*bx[4]*f[11]*dvpar_sq*wvpar+0.1581138830084189*bx[2]*f[9]*dvpar_sq*wvpar+0.1581138830084189*bx[2]*f[8]*dvpar_sq*wvpar+0.1767766952966368*bx[6]*f[7]*dvpar_sq*wvpar+0.1581138830084189*f[4]*bx[7]*dvpar_sq*wvpar+0.1581138830084189*f[2]*bx[5]*dvpar_sq*wvpar+0.1767766952966368*bx[1]*f[4]*dvpar_sq*wvpar+0.1767766952966368*f[1]*bx[3]*dvpar_sq*wvpar+0.1767766952966368*bx[0]*f[2]*dvpar_sq*wvpar+0.1767766952966368*f[0]*bx[2]*dvpar_sq*wvpar+0.02738612787525829*bx[3]*f[18]*dvpar_cu+0.03061862178478971*bx[4]*f[17]*dvpar_cu+0.0273861278752583*bx[2]*f[14]*dvpar_cu+0.03061862178478971*bx[6]*f[13]*dvpar_cu+0.0273861278752583*bx[7]*f[10]*dvpar_cu+0.03061862178478971*bx[1]*f[10]*dvpar_cu+0.02738612787525829*bx[5]*f[6]*dvpar_cu+0.03061862178478971*bx[0]*f[6]*dvpar_cu+0.03061862178478971*bx[3]*f[5]*dvpar_cu+0.03061862178478971*bx[2]*f[3]*dvpar_cu); 
  out[19] += mass*volFact*(0.4898979485566356*bx[6]*f[18]*dvpar*wvpar_sq+0.5477225575051661*bx[2]*f[18]*dvpar*wvpar_sq+0.4898979485566356*bx[7]*f[17]*dvpar*wvpar_sq+0.5477225575051661*bx[1]*f[17]*dvpar*wvpar_sq+0.5477225575051661*bx[3]*f[14]*dvpar*wvpar_sq+0.5477225575051661*bx[3]*f[13]*dvpar*wvpar_sq+0.5477225575051661*bx[5]*f[10]*dvpar*wvpar_sq+0.5477225575051661*bx[4]*f[10]*dvpar*wvpar_sq+0.6123724356957944*bx[0]*f[10]*dvpar*wvpar_sq+0.5477225575051661*f[6]*bx[7]*dvpar*wvpar_sq+0.6123724356957944*bx[1]*f[6]*dvpar*wvpar_sq+0.5477225575051661*f[5]*bx[6]*dvpar*wvpar_sq+0.6123724356957944*bx[2]*f[5]*dvpar*wvpar_sq+0.6123724356957944*bx[3]*f[3]*dvpar*wvpar_sq+0.5656854249492381*bx[6]*f[12]*wvpar_cu+0.632455532033676*bx[2]*f[12]*wvpar_cu+0.5656854249492381*bx[7]*f[11]*wvpar_cu+0.632455532033676*bx[1]*f[11]*wvpar_cu+0.6324555320336759*bx[3]*f[8]*wvpar_cu+0.6324555320336759*bx[3]*f[7]*wvpar_cu+0.632455532033676*f[2]*bx[7]*wvpar_cu+0.632455532033676*f[1]*bx[6]*wvpar_cu+0.6324555320336759*f[4]*bx[5]*wvpar_cu+0.6324555320336759*bx[4]*f[4]*wvpar_cu+0.7071067811865475*bx[0]*f[4]*wvpar_cu+0.7071067811865475*f[0]*bx[3]*wvpar_cu+0.7071067811865475*bx[1]*f[2]*wvpar_cu+0.7071067811865475*f[1]*bx[2]*wvpar_cu+0.1414213562373095*bx[5]*f[19]*dvpar_sq*wvpar+0.1414213562373095*bx[4]*f[19]*dvpar_sq*wvpar+0.1581138830084189*bx[0]*f[19]*dvpar_sq*wvpar+0.1414213562373095*bx[7]*f[16]*dvpar_sq*wvpar+0.1581138830084189*bx[1]*f[16]*dvpar_sq*wvpar+0.1414213562373095*bx[6]*f[15]*dvpar_sq*wvpar+0.1581138830084189*bx[2]*f[15]*dvpar_sq*wvpar+0.1414213562373095*bx[6]*f[12]*dvpar_sq*wvpar+0.1581138830084189*bx[2]*f[12]*dvpar_sq*wvpar+0.1414213562373095*bx[7]*f[11]*dvpar_sq*wvpar+0.1581138830084189*bx[1]*f[11]*dvpar_sq*wvpar+0.1581138830084189*bx[3]*f[9]*dvpar_sq*wvpar+0.1581138830084189*bx[3]*f[8]*dvpar_sq*wvpar+0.1581138830084189*bx[3]*f[7]*dvpar_sq*wvpar+0.1581138830084189*f[2]*bx[7]*dvpar_sq*wvpar+0.1581138830084189*f[1]*bx[6]*dvpar_sq*wvpar+0.1581138830084189*f[4]*bx[5]*dvpar_sq*wvpar+0.1581138830084189*bx[4]*f[4]*dvpar_sq*wvpar+0.1767766952966368*bx[0]*f[4]*dvpar_sq*wvpar+0.1767766952966368*f[0]*bx[3]*dvpar_sq*wvpar+0.1767766952966368*bx[1]*f[2]*dvpar_sq*wvpar+0.1767766952966368*f[1]*bx[2]*dvpar_sq*wvpar+0.02449489742783177*bx[6]*f[18]*dvpar_cu+0.02738612787525829*bx[2]*f[18]*dvpar_cu+0.02449489742783177*bx[7]*f[17]*dvpar_cu+0.02738612787525829*bx[1]*f[17]*dvpar_cu+0.0273861278752583*bx[3]*f[14]*dvpar_cu+0.0273861278752583*bx[3]*f[13]*dvpar_cu+0.02738612787525829*bx[5]*f[10]*dvpar_cu+0.02738612787525829*bx[4]*f[10]*dvpar_cu+0.03061862178478971*bx[0]*f[10]*dvpar_cu+0.0273861278752583*f[6]*bx[7]*dvpar_cu+0.03061862178478971*bx[1]*f[6]*dvpar_cu+0.0273861278752583*f[5]*bx[6]*dvpar_cu+0.03061862178478971*bx[2]*f[5]*dvpar_cu+0.03061862178478971*bx[3]*f[3]*dvpar_cu); 
  out[20] += mass*volFact*(0.5477225575051661*bx[7]*f[18]*dvpar*wvpar_sq+0.3912303982179757*bx[6]*f[17]*dvpar*wvpar_sq+0.6123724356957944*bx[2]*f[17]*dvpar*wvpar_sq+0.3912303982179757*bx[4]*f[13]*dvpar*wvpar_sq+0.6123724356957944*bx[0]*f[13]*dvpar*wvpar_sq+0.5477225575051661*bx[3]*f[10]*dvpar*wvpar_sq+0.6123724356957944*bx[6]*f[6]*dvpar*wvpar_sq+0.5477225575051661*bx[1]*f[5]*dvpar*wvpar_sq+0.6123724356957944*f[3]*bx[4]*dvpar*wvpar_sq+0.6324555320336759*bx[7]*f[12]*wvpar_cu+0.4517539514526256*bx[6]*f[11]*wvpar_cu+0.7071067811865475*bx[2]*f[11]*wvpar_cu+0.4517539514526256*bx[4]*f[7]*wvpar_cu+0.7071067811865475*bx[0]*f[7]*wvpar_cu+0.7071067811865475*f[2]*bx[6]*wvpar_cu+0.6324555320336759*bx[3]*f[4]*wvpar_cu+0.7071067811865475*f[0]*bx[4]*wvpar_cu+0.6324555320336759*bx[1]*f[1]*wvpar_cu+0.1414213562373095*bx[3]*f[19]*dvpar_sq*wvpar+0.1581138830084189*bx[6]*f[16]*dvpar_sq*wvpar+0.1414213562373095*bx[1]*f[15]*dvpar_sq*wvpar+0.1581138830084189*bx[7]*f[12]*dvpar_sq*wvpar+0.1129384878631564*bx[6]*f[11]*dvpar_sq*wvpar+0.1767766952966368*bx[2]*f[11]*dvpar_sq*wvpar+0.1581138830084189*bx[4]*f[9]*dvpar_sq*wvpar+0.1129384878631564*bx[4]*f[7]*dvpar_sq*wvpar+0.1767766952966368*bx[0]*f[7]*dvpar_sq*wvpar+0.1767766952966368*f[2]*bx[6]*dvpar_sq*wvpar+0.1581138830084189*bx[3]*f[4]*dvpar_sq*wvpar+0.1767766952966368*f[0]*bx[4]*dvpar_sq*wvpar+0.1581138830084189*bx[1]*f[1]*dvpar_sq*wvpar+0.0273861278752583*bx[7]*f[18]*dvpar_cu+0.01956151991089878*bx[6]*f[17]*dvpar_cu+0.03061862178478971*bx[2]*f[17]*dvpar_cu+0.01956151991089878*bx[4]*f[13]*dvpar_cu+0.03061862178478971*bx[0]*f[13]*dvpar_cu+0.02738612787525829*bx[3]*f[10]*dvpar_cu+0.03061862178478971*bx[6]*f[6]*dvpar_cu+0.02738612787525829*bx[1]*f[5]*dvpar_cu+0.03061862178478971*f[3]*bx[4]*dvpar_cu); 
  out[21] += mass*volFact*(0.3912303982179757*bx[7]*f[18]*dvpar*wvpar_sq+0.6123724356957944*bx[1]*f[18]*dvpar*wvpar_sq+0.5477225575051661*bx[6]*f[17]*dvpar*wvpar_sq+0.3912303982179757*bx[5]*f[14]*dvpar*wvpar_sq+0.6123724356957944*bx[0]*f[14]*dvpar*wvpar_sq+0.5477225575051661*bx[3]*f[10]*dvpar*wvpar_sq+0.6123724356957944*f[5]*bx[7]*dvpar*wvpar_sq+0.5477225575051661*bx[2]*f[6]*dvpar*wvpar_sq+0.6123724356957944*f[3]*bx[5]*dvpar*wvpar_sq+0.4517539514526256*bx[7]*f[12]*wvpar_cu+0.7071067811865475*bx[1]*f[12]*wvpar_cu+0.6324555320336759*bx[6]*f[11]*wvpar_cu+0.4517539514526256*bx[5]*f[8]*wvpar_cu+0.7071067811865475*bx[0]*f[8]*wvpar_cu+0.7071067811865475*f[1]*bx[7]*wvpar_cu+0.7071067811865475*f[0]*bx[5]*wvpar_cu+0.6324555320336759*bx[3]*f[4]*wvpar_cu+0.6324555320336759*bx[2]*f[2]*wvpar_cu+0.1414213562373095*bx[3]*f[19]*dvpar_sq*wvpar+0.1414213562373095*bx[2]*f[16]*dvpar_sq*wvpar+0.1581138830084189*bx[7]*f[15]*dvpar_sq*wvpar+0.1129384878631564*bx[7]*f[12]*dvpar_sq*wvpar+0.1767766952966368*bx[1]*f[12]*dvpar_sq*wvpar+0.1581138830084189*bx[6]*f[11]*dvpar_sq*wvpar+0.1581138830084189*bx[5]*f[9]*dvpar_sq*wvpar+0.1129384878631564*bx[5]*f[8]*dvpar_sq*wvpar+0.1767766952966368*bx[0]*f[8]*dvpar_sq*wvpar+0.1767766952966368*f[1]*bx[7]*dvpar_sq*wvpar+0.1767766952966368*f[0]*bx[5]*dvpar_sq*wvpar+0.1581138830084189*bx[3]*f[4]*dvpar_sq*wvpar+0.1581138830084189*bx[2]*f[2]*dvpar_sq*wvpar+0.01956151991089878*bx[7]*f[18]*dvpar_cu+0.03061862178478971*bx[1]*f[18]*dvpar_cu+0.0273861278752583*bx[6]*f[17]*dvpar_cu+0.01956151991089878*bx[5]*f[14]*dvpar_cu+0.03061862178478971*bx[0]*f[14]*dvpar_cu+0.02738612787525829*bx[3]*f[10]*dvpar_cu+0.03061862178478971*f[5]*bx[7]*dvpar_cu+0.02738612787525829*bx[2]*f[6]*dvpar_cu+0.03061862178478971*f[3]*bx[5]*dvpar_cu); 
  out[22] += mass*volFact*(0.4898979485566356*bx[3]*f[18]*dvpar*wvpar_sq+0.5477225575051661*bx[5]*f[17]*dvpar*wvpar_sq+0.3912303982179757*bx[4]*f[17]*dvpar*wvpar_sq+0.6123724356957944*bx[0]*f[17]*dvpar*wvpar_sq+0.5477225575051661*bx[6]*f[14]*dvpar*wvpar_sq+0.3912303982179757*bx[6]*f[13]*dvpar*wvpar_sq+0.6123724356957944*bx[2]*f[13]*dvpar*wvpar_sq+0.4898979485566357*bx[7]*f[10]*dvpar*wvpar_sq+0.5477225575051661*bx[1]*f[10]*dvpar*wvpar_sq+0.6123724356957944*bx[4]*f[6]*dvpar*wvpar_sq+0.6123724356957944*f[3]*bx[6]*dvpar*wvpar_sq+0.5477225575051661*bx[3]*f[5]*dvpar*wvpar_sq+0.5656854249492381*bx[3]*f[12]*wvpar_cu+0.6324555320336759*bx[5]*f[11]*wvpar_cu+0.4517539514526256*bx[4]*f[11]*wvpar_cu+0.7071067811865475*bx[0]*f[11]*wvpar_cu+0.6324555320336759*bx[6]*f[8]*wvpar_cu+0.4517539514526256*bx[6]*f[7]*wvpar_cu+0.7071067811865475*bx[2]*f[7]*wvpar_cu+0.5656854249492381*f[4]*bx[7]*wvpar_cu+0.7071067811865475*f[0]*bx[6]*wvpar_cu+0.632455532033676*bx[1]*f[4]*wvpar_cu+0.7071067811865475*f[2]*bx[4]*wvpar_cu+0.632455532033676*f[1]*bx[3]*wvpar_cu+0.1264911064067352*bx[7]*f[19]*dvpar_sq*wvpar+0.1414213562373095*bx[1]*f[19]*dvpar_sq*wvpar+0.1581138830084189*bx[4]*f[16]*dvpar_sq*wvpar+0.1414213562373095*bx[3]*f[15]*dvpar_sq*wvpar+0.1414213562373095*bx[3]*f[12]*dvpar_sq*wvpar+0.1581138830084189*bx[5]*f[11]*dvpar_sq*wvpar+0.1129384878631564*bx[4]*f[11]*dvpar_sq*wvpar+0.1767766952966368*bx[0]*f[11]*dvpar_sq*wvpar+0.1581138830084189*bx[6]*f[9]*dvpar_sq*wvpar+0.1581138830084189*bx[6]*f[8]*dvpar_sq*wvpar+0.1129384878631564*bx[6]*f[7]*dvpar_sq*wvpar+0.1767766952966368*bx[2]*f[7]*dvpar_sq*wvpar+0.1414213562373095*f[4]*bx[7]*dvpar_sq*wvpar+0.1767766952966368*f[0]*bx[6]*dvpar_sq*wvpar+0.1581138830084189*bx[1]*f[4]*dvpar_sq*wvpar+0.1767766952966368*f[2]*bx[4]*dvpar_sq*wvpar+0.1581138830084189*f[1]*bx[3]*dvpar_sq*wvpar+0.02449489742783177*bx[3]*f[18]*dvpar_cu+0.0273861278752583*bx[5]*f[17]*dvpar_cu+0.01956151991089878*bx[4]*f[17]*dvpar_cu+0.03061862178478971*bx[0]*f[17]*dvpar_cu+0.0273861278752583*bx[6]*f[14]*dvpar_cu+0.01956151991089878*bx[6]*f[13]*dvpar_cu+0.03061862178478971*bx[2]*f[13]*dvpar_cu+0.02449489742783178*bx[7]*f[10]*dvpar_cu+0.0273861278752583*bx[1]*f[10]*dvpar_cu+0.03061862178478971*bx[4]*f[6]*dvpar_cu+0.03061862178478971*f[3]*bx[6]*dvpar_cu+0.0273861278752583*bx[3]*f[5]*dvpar_cu); 
  out[23] += mass*volFact*(0.3912303982179757*bx[5]*f[18]*dvpar*wvpar_sq+0.5477225575051661*bx[4]*f[18]*dvpar*wvpar_sq+0.6123724356957944*bx[0]*f[18]*dvpar*wvpar_sq+0.4898979485566356*bx[3]*f[17]*dvpar*wvpar_sq+0.3912303982179757*bx[7]*f[14]*dvpar*wvpar_sq+0.6123724356957944*bx[1]*f[14]*dvpar*wvpar_sq+0.5477225575051661*bx[7]*f[13]*dvpar*wvpar_sq+0.4898979485566357*bx[6]*f[10]*dvpar*wvpar_sq+0.5477225575051661*bx[2]*f[10]*dvpar*wvpar_sq+0.6123724356957944*f[3]*bx[7]*dvpar*wvpar_sq+0.5477225575051661*bx[3]*f[6]*dvpar*wvpar_sq+0.6123724356957944*bx[5]*f[5]*dvpar*wvpar_sq+0.4517539514526256*bx[5]*f[12]*wvpar_cu+0.6324555320336759*bx[4]*f[12]*wvpar_cu+0.7071067811865475*bx[0]*f[12]*wvpar_cu+0.5656854249492381*bx[3]*f[11]*wvpar_cu+0.4517539514526256*bx[7]*f[8]*wvpar_cu+0.7071067811865475*bx[1]*f[8]*wvpar_cu+0.6324555320336759*bx[7]*f[7]*wvpar_cu+0.7071067811865475*f[0]*bx[7]*wvpar_cu+0.5656854249492381*f[4]*bx[6]*wvpar_cu+0.7071067811865475*f[1]*bx[5]*wvpar_cu+0.632455532033676*bx[2]*f[4]*wvpar_cu+0.632455532033676*f[2]*bx[3]*wvpar_cu+0.1264911064067352*bx[6]*f[19]*dvpar_sq*wvpar+0.1414213562373095*bx[2]*f[19]*dvpar_sq*wvpar+0.1414213562373095*bx[3]*f[16]*dvpar_sq*wvpar+0.1581138830084189*bx[5]*f[15]*dvpar_sq*wvpar+0.1129384878631564*bx[5]*f[12]*dvpar_sq*wvpar+0.1581138830084189*bx[4]*f[12]*dvpar_sq*wvpar+0.1767766952966368*bx[0]*f[12]*dvpar_sq*wvpar+0.1414213562373095*bx[3]*f[11]*dvpar_sq*wvpar+0.1581138830084189*bx[7]*f[9]*dvpar_sq*wvpar+0.1129384878631564*bx[7]*f[8]*dvpar_sq*wvpar+0.1767766952966368*bx[1]*f[8]*dvpar_sq*wvpar+0.1581138830084189*bx[7]*f[7]*dvpar_sq*wvpar+0.1767766952966368*f[0]*bx[7]*dvpar_sq*wvpar+0.1414213562373095*f[4]*bx[6]*dvpar_sq*wvpar+0.1767766952966368*f[1]*bx[5]*dvpar_sq*wvpar+0.1581138830084189*bx[2]*f[4]*dvpar_sq*wvpar+0.1581138830084189*f[2]*bx[3]*dvpar_sq*wvpar+0.01956151991089878*bx[5]*f[18]*dvpar_cu+0.0273861278752583*bx[4]*f[18]*dvpar_cu+0.03061862178478971*bx[0]*f[18]*dvpar_cu+0.02449489742783177*bx[3]*f[17]*dvpar_cu+0.01956151991089878*bx[7]*f[14]*dvpar_cu+0.03061862178478971*bx[1]*f[14]*dvpar_cu+0.0273861278752583*bx[7]*f[13]*dvpar_cu+0.02449489742783178*bx[6]*f[10]*dvpar_cu+0.0273861278752583*bx[2]*f[10]*dvpar_cu+0.03061862178478971*f[3]*bx[7]*dvpar_cu+0.0273861278752583*bx[3]*f[6]*dvpar_cu+0.03061862178478971*bx[5]*f[5]*dvpar_cu); 
  out[24] += mass*volFact*(0.6123724356957944*by[7]*f[18]*dvpar*wvpar_sq+0.6123724356957944*by[6]*f[17]*dvpar*wvpar_sq+0.6123724356957944*by[5]*f[14]*dvpar*wvpar_sq+0.6123724356957944*by[4]*f[13]*dvpar*wvpar_sq+0.6123724356957944*by[3]*f[10]*dvpar*wvpar_sq+0.6123724356957944*by[2]*f[6]*dvpar*wvpar_sq+0.6123724356957944*by[1]*f[5]*dvpar*wvpar_sq+0.6123724356957944*by[0]*f[3]*dvpar*wvpar_sq+0.7071067811865475*by[7]*f[12]*wvpar_cu+0.7071067811865475*by[6]*f[11]*wvpar_cu+0.7071067811865475*by[5]*f[8]*wvpar_cu+0.7071067811865475*by[4]*f[7]*wvpar_cu+0.7071067811865475*by[3]*f[4]*wvpar_cu+0.7071067811865475*by[2]*f[2]*wvpar_cu+0.7071067811865475*by[1]*f[1]*wvpar_cu+0.7071067811865475*by[0]*f[0]*wvpar_cu+0.1581138830084189*by[3]*f[19]*dvpar_sq*wvpar+0.1581138830084189*by[2]*f[16]*dvpar_sq*wvpar+0.1581138830084189*by[1]*f[15]*dvpar_sq*wvpar+0.1767766952966368*by[7]*f[12]*dvpar_sq*wvpar+0.1767766952966368*by[6]*f[11]*dvpar_sq*wvpar+0.1581138830084189*by[0]*f[9]*dvpar_sq*wvpar+0.1767766952966368*by[5]*f[8]*dvpar_sq*wvpar+0.1767766952966368*by[4]*f[7]*dvpar_sq*wvpar+0.1767766952966368*by[3]*f[4]*dvpar_sq*wvpar+0.1767766952966368*by[2]*f[2]*dvpar_sq*wvpar+0.1767766952966368*by[1]*f[1]*dvpar_sq*wvpar+0.1767766952966368*by[0]*f[0]*dvpar_sq*wvpar+0.03061862178478971*by[7]*f[18]*dvpar_cu+0.03061862178478971*by[6]*f[17]*dvpar_cu+0.03061862178478971*by[5]*f[14]*dvpar_cu+0.03061862178478971*by[4]*f[13]*dvpar_cu+0.03061862178478971*by[3]*f[10]*dvpar_cu+0.03061862178478971*by[2]*f[6]*dvpar_cu+0.03061862178478971*by[1]*f[5]*dvpar_cu+0.03061862178478971*by[0]*f[3]*dvpar_cu); 
  out[25] += mass*volFact*(0.6123724356957944*by[5]*f[18]*dvpar*wvpar_sq+0.5477225575051661*by[3]*f[17]*dvpar*wvpar_sq+0.6123724356957944*by[7]*f[14]*dvpar*wvpar_sq+0.5477225575051661*by[1]*f[13]*dvpar*wvpar_sq+0.5477225575051661*by[6]*f[10]*dvpar*wvpar_sq+0.6123724356957944*by[2]*f[10]*dvpar*wvpar_sq+0.6123724356957944*by[3]*f[6]*dvpar*wvpar_sq+0.5477225575051661*by[4]*f[5]*dvpar*wvpar_sq+0.6123724356957944*by[0]*f[5]*dvpar*wvpar_sq+0.6123724356957944*by[1]*f[3]*dvpar*wvpar_sq+0.7071067811865475*by[5]*f[12]*wvpar_cu+0.632455532033676*by[3]*f[11]*wvpar_cu+0.7071067811865475*by[7]*f[8]*wvpar_cu+0.6324555320336759*by[1]*f[7]*wvpar_cu+0.632455532033676*f[4]*by[6]*wvpar_cu+0.7071067811865475*by[2]*f[4]*wvpar_cu+0.6324555320336759*f[1]*by[4]*wvpar_cu+0.7071067811865475*f[2]*by[3]*wvpar_cu+0.7071067811865475*by[0]*f[1]*wvpar_cu+0.7071067811865475*f[0]*by[1]*wvpar_cu+0.1414213562373095*by[6]*f[19]*dvpar_sq*wvpar+0.1581138830084189*by[2]*f[19]*dvpar_sq*wvpar+0.1581138830084189*by[3]*f[16]*dvpar_sq*wvpar+0.1414213562373095*by[4]*f[15]*dvpar_sq*wvpar+0.1581138830084189*by[0]*f[15]*dvpar_sq*wvpar+0.1767766952966368*by[5]*f[12]*dvpar_sq*wvpar+0.1581138830084189*by[3]*f[11]*dvpar_sq*wvpar+0.1581138830084189*by[1]*f[9]*dvpar_sq*wvpar+0.1767766952966368*by[7]*f[8]*dvpar_sq*wvpar+0.1581138830084189*by[1]*f[7]*dvpar_sq*wvpar+0.1581138830084189*f[4]*by[6]*dvpar_sq*wvpar+0.1767766952966368*by[2]*f[4]*dvpar_sq*wvpar+0.1581138830084189*f[1]*by[4]*dvpar_sq*wvpar+0.1767766952966368*f[2]*by[3]*dvpar_sq*wvpar+0.1767766952966368*by[0]*f[1]*dvpar_sq*wvpar+0.1767766952966368*f[0]*by[1]*dvpar_sq*wvpar+0.03061862178478971*by[5]*f[18]*dvpar_cu+0.02738612787525829*by[3]*f[17]*dvpar_cu+0.03061862178478971*by[7]*f[14]*dvpar_cu+0.0273861278752583*by[1]*f[13]*dvpar_cu+0.0273861278752583*by[6]*f[10]*dvpar_cu+0.03061862178478971*by[2]*f[10]*dvpar_cu+0.03061862178478971*by[3]*f[6]*dvpar_cu+0.02738612787525829*by[4]*f[5]*dvpar_cu+0.03061862178478971*by[0]*f[5]*dvpar_cu+0.03061862178478971*by[1]*f[3]*dvpar_cu); 
  out[26] += mass*volFact*(0.5477225575051661*by[3]*f[18]*dvpar*wvpar_sq+0.6123724356957944*by[4]*f[17]*dvpar*wvpar_sq+0.5477225575051661*by[2]*f[14]*dvpar*wvpar_sq+0.6123724356957944*by[6]*f[13]*dvpar*wvpar_sq+0.5477225575051661*by[7]*f[10]*dvpar*wvpar_sq+0.6123724356957944*by[1]*f[10]*dvpar*wvpar_sq+0.5477225575051661*by[5]*f[6]*dvpar*wvpar_sq+0.6123724356957944*by[0]*f[6]*dvpar*wvpar_sq+0.6123724356957944*by[3]*f[5]*dvpar*wvpar_sq+0.6123724356957944*by[2]*f[3]*dvpar*wvpar_sq+0.632455532033676*by[3]*f[12]*wvpar_cu+0.7071067811865475*by[4]*f[11]*wvpar_cu+0.6324555320336759*by[2]*f[8]*wvpar_cu+0.7071067811865475*by[6]*f[7]*wvpar_cu+0.632455532033676*f[4]*by[7]*wvpar_cu+0.6324555320336759*f[2]*by[5]*wvpar_cu+0.7071067811865475*by[1]*f[4]*wvpar_cu+0.7071067811865475*f[1]*by[3]*wvpar_cu+0.7071067811865475*by[0]*f[2]*wvpar_cu+0.7071067811865475*f[0]*by[2]*wvpar_cu+0.1414213562373095*by[7]*f[19]*dvpar_sq*wvpar+0.1581138830084189*by[1]*f[19]*dvpar_sq*wvpar+0.1414213562373095*by[5]*f[16]*dvpar_sq*wvpar+0.1581138830084189*by[0]*f[16]*dvpar_sq*wvpar+0.1581138830084189*by[3]*f[15]*dvpar_sq*wvpar+0.1581138830084189*by[3]*f[12]*dvpar_sq*wvpar+0.1767766952966368*by[4]*f[11]*dvpar_sq*wvpar+0.1581138830084189*by[2]*f[9]*dvpar_sq*wvpar+0.1581138830084189*by[2]*f[8]*dvpar_sq*wvpar+0.1767766952966368*by[6]*f[7]*dvpar_sq*wvpar+0.1581138830084189*f[4]*by[7]*dvpar_sq*wvpar+0.1581138830084189*f[2]*by[5]*dvpar_sq*wvpar+0.1767766952966368*by[1]*f[4]*dvpar_sq*wvpar+0.1767766952966368*f[1]*by[3]*dvpar_sq*wvpar+0.1767766952966368*by[0]*f[2]*dvpar_sq*wvpar+0.1767766952966368*f[0]*by[2]*dvpar_sq*wvpar+0.02738612787525829*by[3]*f[18]*dvpar_cu+0.03061862178478971*by[4]*f[17]*dvpar_cu+0.0273861278752583*by[2]*f[14]*dvpar_cu+0.03061862178478971*by[6]*f[13]*dvpar_cu+0.0273861278752583*by[7]*f[10]*dvpar_cu+0.03061862178478971*by[1]*f[10]*dvpar_cu+0.02738612787525829*by[5]*f[6]*dvpar_cu+0.03061862178478971*by[0]*f[6]*dvpar_cu+0.03061862178478971*by[3]*f[5]*dvpar_cu+0.03061862178478971*by[2]*f[3]*dvpar_cu); 
  out[27] += mass*volFact*(0.4898979485566356*by[6]*f[18]*dvpar*wvpar_sq+0.5477225575051661*by[2]*f[18]*dvpar*wvpar_sq+0.4898979485566356*by[7]*f[17]*dvpar*wvpar_sq+0.5477225575051661*by[1]*f[17]*dvpar*wvpar_sq+0.5477225575051661*by[3]*f[14]*dvpar*wvpar_sq+0.5477225575051661*by[3]*f[13]*dvpar*wvpar_sq+0.5477225575051661*by[5]*f[10]*dvpar*wvpar_sq+0.5477225575051661*by[4]*f[10]*dvpar*wvpar_sq+0.6123724356957944*by[0]*f[10]*dvpar*wvpar_sq+0.5477225575051661*f[6]*by[7]*dvpar*wvpar_sq+0.6123724356957944*by[1]*f[6]*dvpar*wvpar_sq+0.5477225575051661*f[5]*by[6]*dvpar*wvpar_sq+0.6123724356957944*by[2]*f[5]*dvpar*wvpar_sq+0.6123724356957944*by[3]*f[3]*dvpar*wvpar_sq+0.5656854249492381*by[6]*f[12]*wvpar_cu+0.632455532033676*by[2]*f[12]*wvpar_cu+0.5656854249492381*by[7]*f[11]*wvpar_cu+0.632455532033676*by[1]*f[11]*wvpar_cu+0.6324555320336759*by[3]*f[8]*wvpar_cu+0.6324555320336759*by[3]*f[7]*wvpar_cu+0.632455532033676*f[2]*by[7]*wvpar_cu+0.632455532033676*f[1]*by[6]*wvpar_cu+0.6324555320336759*f[4]*by[5]*wvpar_cu+0.6324555320336759*by[4]*f[4]*wvpar_cu+0.7071067811865475*by[0]*f[4]*wvpar_cu+0.7071067811865475*f[0]*by[3]*wvpar_cu+0.7071067811865475*by[1]*f[2]*wvpar_cu+0.7071067811865475*f[1]*by[2]*wvpar_cu+0.1414213562373095*by[5]*f[19]*dvpar_sq*wvpar+0.1414213562373095*by[4]*f[19]*dvpar_sq*wvpar+0.1581138830084189*by[0]*f[19]*dvpar_sq*wvpar+0.1414213562373095*by[7]*f[16]*dvpar_sq*wvpar+0.1581138830084189*by[1]*f[16]*dvpar_sq*wvpar+0.1414213562373095*by[6]*f[15]*dvpar_sq*wvpar+0.1581138830084189*by[2]*f[15]*dvpar_sq*wvpar+0.1414213562373095*by[6]*f[12]*dvpar_sq*wvpar+0.1581138830084189*by[2]*f[12]*dvpar_sq*wvpar+0.1414213562373095*by[7]*f[11]*dvpar_sq*wvpar+0.1581138830084189*by[1]*f[11]*dvpar_sq*wvpar+0.1581138830084189*by[3]*f[9]*dvpar_sq*wvpar+0.1581138830084189*by[3]*f[8]*dvpar_sq*wvpar+0.1581138830084189*by[3]*f[7]*dvpar_sq*wvpar+0.1581138830084189*f[2]*by[7]*dvpar_sq*wvpar+0.1581138830084189*f[1]*by[6]*dvpar_sq*wvpar+0.1581138830084189*f[4]*by[5]*dvpar_sq*wvpar+0.1581138830084189*by[4]*f[4]*dvpar_sq*wvpar+0.1767766952966368*by[0]*f[4]*dvpar_sq*wvpar+0.1767766952966368*f[0]*by[3]*dvpar_sq*wvpar+0.1767766952966368*by[1]*f[2]*dvpar_sq*wvpar+0.1767766952966368*f[1]*by[2]*dvpar_sq*wvpar+0.02449489742783177*by[6]*f[18]*dvpar_cu+0.02738612787525829*by[2]*f[18]*dvpar_cu+0.02449489742783177*by[7]*f[17]*dvpar_cu+0.02738612787525829*by[1]*f[17]*dvpar_cu+0.0273861278752583*by[3]*f[14]*dvpar_cu+0.0273861278752583*by[3]*f[13]*dvpar_cu+0.02738612787525829*by[5]*f[10]*dvpar_cu+0.02738612787525829*by[4]*f[10]*dvpar_cu+0.03061862178478971*by[0]*f[10]*dvpar_cu+0.0273861278752583*f[6]*by[7]*dvpar_cu+0.03061862178478971*by[1]*f[6]*dvpar_cu+0.0273861278752583*f[5]*by[6]*dvpar_cu+0.03061862178478971*by[2]*f[5]*dvpar_cu+0.03061862178478971*by[3]*f[3]*dvpar_cu); 
  out[28] += mass*volFact*(0.5477225575051661*by[7]*f[18]*dvpar*wvpar_sq+0.3912303982179757*by[6]*f[17]*dvpar*wvpar_sq+0.6123724356957944*by[2]*f[17]*dvpar*wvpar_sq+0.3912303982179757*by[4]*f[13]*dvpar*wvpar_sq+0.6123724356957944*by[0]*f[13]*dvpar*wvpar_sq+0.5477225575051661*by[3]*f[10]*dvpar*wvpar_sq+0.6123724356957944*by[6]*f[6]*dvpar*wvpar_sq+0.5477225575051661*by[1]*f[5]*dvpar*wvpar_sq+0.6123724356957944*f[3]*by[4]*dvpar*wvpar_sq+0.6324555320336759*by[7]*f[12]*wvpar_cu+0.4517539514526256*by[6]*f[11]*wvpar_cu+0.7071067811865475*by[2]*f[11]*wvpar_cu+0.4517539514526256*by[4]*f[7]*wvpar_cu+0.7071067811865475*by[0]*f[7]*wvpar_cu+0.7071067811865475*f[2]*by[6]*wvpar_cu+0.6324555320336759*by[3]*f[4]*wvpar_cu+0.7071067811865475*f[0]*by[4]*wvpar_cu+0.6324555320336759*by[1]*f[1]*wvpar_cu+0.1414213562373095*by[3]*f[19]*dvpar_sq*wvpar+0.1581138830084189*by[6]*f[16]*dvpar_sq*wvpar+0.1414213562373095*by[1]*f[15]*dvpar_sq*wvpar+0.1581138830084189*by[7]*f[12]*dvpar_sq*wvpar+0.1129384878631564*by[6]*f[11]*dvpar_sq*wvpar+0.1767766952966368*by[2]*f[11]*dvpar_sq*wvpar+0.1581138830084189*by[4]*f[9]*dvpar_sq*wvpar+0.1129384878631564*by[4]*f[7]*dvpar_sq*wvpar+0.1767766952966368*by[0]*f[7]*dvpar_sq*wvpar+0.1767766952966368*f[2]*by[6]*dvpar_sq*wvpar+0.1581138830084189*by[3]*f[4]*dvpar_sq*wvpar+0.1767766952966368*f[0]*by[4]*dvpar_sq*wvpar+0.1581138830084189*by[1]*f[1]*dvpar_sq*wvpar+0.0273861278752583*by[7]*f[18]*dvpar_cu+0.01956151991089878*by[6]*f[17]*dvpar_cu+0.03061862178478971*by[2]*f[17]*dvpar_cu+0.01956151991089878*by[4]*f[13]*dvpar_cu+0.03061862178478971*by[0]*f[13]*dvpar_cu+0.02738612787525829*by[3]*f[10]*dvpar_cu+0.03061862178478971*by[6]*f[6]*dvpar_cu+0.02738612787525829*by[1]*f[5]*dvpar_cu+0.03061862178478971*f[3]*by[4]*dvpar_cu); 
  out[29] += mass*volFact*(0.3912303982179757*by[7]*f[18]*dvpar*wvpar_sq+0.6123724356957944*by[1]*f[18]*dvpar*wvpar_sq+0.5477225575051661*by[6]*f[17]*dvpar*wvpar_sq+0.3912303982179757*by[5]*f[14]*dvpar*wvpar_sq+0.6123724356957944*by[0]*f[14]*dvpar*wvpar_sq+0.5477225575051661*by[3]*f[10]*dvpar*wvpar_sq+0.6123724356957944*f[5]*by[7]*dvpar*wvpar_sq+0.5477225575051661*by[2]*f[6]*dvpar*wvpar_sq+0.6123724356957944*f[3]*by[5]*dvpar*wvpar_sq+0.4517539514526256*by[7]*f[12]*wvpar_cu+0.7071067811865475*by[1]*f[12]*wvpar_cu+0.6324555320336759*by[6]*f[11]*wvpar_cu+0.4517539514526256*by[5]*f[8]*wvpar_cu+0.7071067811865475*by[0]*f[8]*wvpar_cu+0.7071067811865475*f[1]*by[7]*wvpar_cu+0.7071067811865475*f[0]*by[5]*wvpar_cu+0.6324555320336759*by[3]*f[4]*wvpar_cu+0.6324555320336759*by[2]*f[2]*wvpar_cu+0.1414213562373095*by[3]*f[19]*dvpar_sq*wvpar+0.1414213562373095*by[2]*f[16]*dvpar_sq*wvpar+0.1581138830084189*by[7]*f[15]*dvpar_sq*wvpar+0.1129384878631564*by[7]*f[12]*dvpar_sq*wvpar+0.1767766952966368*by[1]*f[12]*dvpar_sq*wvpar+0.1581138830084189*by[6]*f[11]*dvpar_sq*wvpar+0.1581138830084189*by[5]*f[9]*dvpar_sq*wvpar+0.1129384878631564*by[5]*f[8]*dvpar_sq*wvpar+0.1767766952966368*by[0]*f[8]*dvpar_sq*wvpar+0.1767766952966368*f[1]*by[7]*dvpar_sq*wvpar+0.1767766952966368*f[0]*by[5]*dvpar_sq*wvpar+0.1581138830084189*by[3]*f[4]*dvpar_sq*wvpar+0.1581138830084189*by[2]*f[2]*dvpar_sq*wvpar+0.01956151991089878*by[7]*f[18]*dvpar_cu+0.03061862178478971*by[1]*f[18]*dvpar_cu+0.0273861278752583*by[6]*f[17]*dvpar_cu+0.01956151991089878*by[5]*f[14]*dvpar_cu+0.03061862178478971*by[0]*f[14]*dvpar_cu+0.02738612787525829*by[3]*f[10]*dvpar_cu+0.03061862178478971*f[5]*by[7]*dvpar_cu+0.02738612787525829*by[2]*f[6]*dvpar_cu+0.03061862178478971*f[3]*by[5]*dvpar_cu); 
  out[30] += mass*volFact*(0.4898979485566356*by[3]*f[18]*dvpar*wvpar_sq+0.5477225575051661*by[5]*f[17]*dvpar*wvpar_sq+0.3912303982179757*by[4]*f[17]*dvpar*wvpar_sq+0.6123724356957944*by[0]*f[17]*dvpar*wvpar_sq+0.5477225575051661*by[6]*f[14]*dvpar*wvpar_sq+0.3912303982179757*by[6]*f[13]*dvpar*wvpar_sq+0.6123724356957944*by[2]*f[13]*dvpar*wvpar_sq+0.4898979485566357*by[7]*f[10]*dvpar*wvpar_sq+0.5477225575051661*by[1]*f[10]*dvpar*wvpar_sq+0.6123724356957944*by[4]*f[6]*dvpar*wvpar_sq+0.6123724356957944*f[3]*by[6]*dvpar*wvpar_sq+0.5477225575051661*by[3]*f[5]*dvpar*wvpar_sq+0.5656854249492381*by[3]*f[12]*wvpar_cu+0.6324555320336759*by[5]*f[11]*wvpar_cu+0.4517539514526256*by[4]*f[11]*wvpar_cu+0.7071067811865475*by[0]*f[11]*wvpar_cu+0.6324555320336759*by[6]*f[8]*wvpar_cu+0.4517539514526256*by[6]*f[7]*wvpar_cu+0.7071067811865475*by[2]*f[7]*wvpar_cu+0.5656854249492381*f[4]*by[7]*wvpar_cu+0.7071067811865475*f[0]*by[6]*wvpar_cu+0.632455532033676*by[1]*f[4]*wvpar_cu+0.7071067811865475*f[2]*by[4]*wvpar_cu+0.632455532033676*f[1]*by[3]*wvpar_cu+0.1264911064067352*by[7]*f[19]*dvpar_sq*wvpar+0.1414213562373095*by[1]*f[19]*dvpar_sq*wvpar+0.1581138830084189*by[4]*f[16]*dvpar_sq*wvpar+0.1414213562373095*by[3]*f[15]*dvpar_sq*wvpar+0.1414213562373095*by[3]*f[12]*dvpar_sq*wvpar+0.1581138830084189*by[5]*f[11]*dvpar_sq*wvpar+0.1129384878631564*by[4]*f[11]*dvpar_sq*wvpar+0.1767766952966368*by[0]*f[11]*dvpar_sq*wvpar+0.1581138830084189*by[6]*f[9]*dvpar_sq*wvpar+0.1581138830084189*by[6]*f[8]*dvpar_sq*wvpar+0.1129384878631564*by[6]*f[7]*dvpar_sq*wvpar+0.1767766952966368*by[2]*f[7]*dvpar_sq*wvpar+0.1414213562373095*f[4]*by[7]*dvpar_sq*wvpar+0.1767766952966368*f[0]*by[6]*dvpar_sq*wvpar+0.1581138830084189*by[1]*f[4]*dvpar_sq*wvpar+0.1767766952966368*f[2]*by[4]*dvpar_sq*wvpar+0.1581138830084189*f[1]*by[3]*dvpar_sq*wvpar+0.02449489742783177*by[3]*f[18]*dvpar_cu+0.0273861278752583*by[5]*f[17]*dvpar_cu+0.01956151991089878*by[4]*f[17]*dvpar_cu+0.03061862178478971*by[0]*f[17]*dvpar_cu+0.0273861278752583*by[6]*f[14]*dvpar_cu+0.01956151991089878*by[6]*f[13]*dvpar_cu+0.03061862178478971*by[2]*f[13]*dvpar_cu+0.02449489742783178*by[7]*f[10]*dvpar_cu+0.0273861278752583*by[1]*f[10]*dvpar_cu+0.03061862178478971*by[4]*f[6]*dvpar_cu+0.03061862178478971*f[3]*by[6]*dvpar_cu+0.0273861278752583*by[3]*f[5]*dvpar_cu); 
  out[31] += mass*volFact*(0.3912303982179757*by[5]*f[18]*dvpar*wvpar_sq+0.5477225575051661*by[4]*f[18]*dvpar*wvpar_sq+0.6123724356957944*by[0]*f[18]*dvpar*wvpar_sq+0.4898979485566356*by[3]*f[17]*dvpar*wvpar_sq+0.3912303982179757*by[7]*f[14]*dvpar*wvpar_sq+0.6123724356957944*by[1]*f[14]*dvpar*wvpar_sq+0.5477225575051661*by[7]*f[13]*dvpar*wvpar_sq+0.4898979485566357*by[6]*f[10]*dvpar*wvpar_sq+0.5477225575051661*by[2]*f[10]*dvpar*wvpar_sq+0.6123724356957944*f[3]*by[7]*dvpar*wvpar_sq+0.5477225575051661*by[3]*f[6]*dvpar*wvpar_sq+0.6123724356957944*by[5]*f[5]*dvpar*wvpar_sq+0.4517539514526256*by[5]*f[12]*wvpar_cu+0.6324555320336759*by[4]*f[12]*wvpar_cu+0.7071067811865475*by[0]*f[12]*wvpar_cu+0.5656854249492381*by[3]*f[11]*wvpar_cu+0.4517539514526256*by[7]*f[8]*wvpar_cu+0.7071067811865475*by[1]*f[8]*wvpar_cu+0.6324555320336759*by[7]*f[7]*wvpar_cu+0.7071067811865475*f[0]*by[7]*wvpar_cu+0.5656854249492381*f[4]*by[6]*wvpar_cu+0.7071067811865475*f[1]*by[5]*wvpar_cu+0.632455532033676*by[2]*f[4]*wvpar_cu+0.632455532033676*f[2]*by[3]*wvpar_cu+0.1264911064067352*by[6]*f[19]*dvpar_sq*wvpar+0.1414213562373095*by[2]*f[19]*dvpar_sq*wvpar+0.1414213562373095*by[3]*f[16]*dvpar_sq*wvpar+0.1581138830084189*by[5]*f[15]*dvpar_sq*wvpar+0.1129384878631564*by[5]*f[12]*dvpar_sq*wvpar+0.1581138830084189*by[4]*f[12]*dvpar_sq*wvpar+0.1767766952966368*by[0]*f[12]*dvpar_sq*wvpar+0.1414213562373095*by[3]*f[11]*dvpar_sq*wvpar+0.1581138830084189*by[7]*f[9]*dvpar_sq*wvpar+0.1129384878631564*by[7]*f[8]*dvpar_sq*wvpar+0.1767766952966368*by[1]*f[8]*dvpar_sq*wvpar+0.1581138830084189*by[7]*f[7]*dvpar_sq*wvpar+0.1767766952966368*f[0]*by[7]*dvpar_sq*wvpar+0.1414213562373095*f[4]*by[6]*dvpar_sq*wvpar+0.1767766952966368*f[1]*by[5]*dvpar_sq*wvpar+0.1581138830084189*by[2]*f[4]*dvpar_sq*wvpar+0.1581138830084189*f[2]*by[3]*dvpar_sq*wvpar+0.01956151991089878*by[5]*f[18]*dvpar_cu+0.0273861278752583*by[4]*f[18]*dvpar_cu+0.03061862178478971*by[0]*f[18]*dvpar_cu+0.02449489742783177*by[3]*f[17]*dvpar_cu+0.01956151991089878*by[7]*f[14]*dvpar_cu+0.03061862178478971*by[1]*f[14]*dvpar_cu+0.0273861278752583*by[7]*f[13]*dvpar_cu+0.02449489742783178*by[6]*f[10]*dvpar_cu+0.0273861278752583*by[2]*f[10]*dvpar_cu+0.03061862178478971*f[3]*by[7]*dvpar_cu+0.0273861278752583*by[3]*f[6]*dvpar_cu+0.03061862178478971*by[5]*f[5]*dvpar_cu); 
} 
