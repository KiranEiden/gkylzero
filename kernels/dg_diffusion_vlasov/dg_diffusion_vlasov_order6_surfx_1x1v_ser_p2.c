#include <gkyl_dg_diffusion_vlasov_kernels.h>

GKYL_CU_DH double dg_diffusion_vlasov_order6_surfx_1x1v_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],6.);

  out[0] += 0.0625*(563.489130329947*coeff[0]*qr[4]+563.489130329947*coeff[0]*ql[4]-1126.978260659894*coeff[0]*qc[4]-545.5960043841961*coeff[0]*qr[1]+545.5960043841961*coeff[0]*ql[1]+315.0*coeff[0]*qr[0]+315.0*coeff[0]*ql[0]-630.0*coeff[0]*qc[0])*Jfac; 
  out[1] += 0.0078125*(6587.944671898812*coeff[0]*qr[4]-6587.944671898812*coeff[0]*ql[4]-7245.0*coeff[0]*qr[1]-7245.0*coeff[0]*ql[1]-15750.0*coeff[0]*qc[1]+4364.768035073569*coeff[0]*qr[0]-4364.768035073569*coeff[0]*ql[0])*Jfac; 
  out[2] += 0.0625*(563.4891303299469*coeff[0]*qr[6]+563.4891303299469*coeff[0]*ql[6]-1126.978260659894*coeff[0]*qc[6]-545.5960043841961*coeff[0]*qr[3]+545.5960043841961*coeff[0]*ql[3]+315.0*coeff[0]*qr[2]+315.0*coeff[0]*ql[2]-630.0*coeff[0]*qc[2])*Jfac; 
  out[3] += 0.0078125*(6587.944671898817*coeff[0]*qr[6]-6587.944671898817*coeff[0]*ql[6]-7245.0*coeff[0]*qr[3]-7245.0*coeff[0]*ql[3]-15750.0*coeff[0]*qc[3]+4364.768035073569*coeff[0]*qr[2]-4364.768035073569*coeff[0]*ql[2])*Jfac; 
  out[4] += -0.0078125*(405.0*coeff[0]*qr[4]+405.0*coeff[0]*ql[4]+18090.0*coeff[0]*qc[4]+1568.558255214003*coeff[0]*qr[1]-1568.558255214003*coeff[0]*ql[1]-1609.968943799849*coeff[0]*qr[0]-1609.968943799849*coeff[0]*ql[0]+3219.937887599698*coeff[0]*qc[0])*Jfac; 
  out[5] += -0.0625*(545.5960043841964*coeff[0]*qr[7]-545.5960043841964*coeff[0]*ql[7]-315.0*coeff[0]*qr[5]-315.0*coeff[0]*ql[5]+630.0*coeff[0]*qc[5])*Jfac; 
  out[6] += -0.0078125*(405.0*coeff[0]*qr[6]+405.0*coeff[0]*ql[6]+18090.0*coeff[0]*qc[6]+1568.558255214004*coeff[0]*qr[3]-1568.558255214004*coeff[0]*ql[3]-1609.968943799848*coeff[0]*qr[2]-1609.968943799848*coeff[0]*ql[2]+3219.937887599697*coeff[0]*qc[2])*Jfac; 
  out[7] += -0.0078125*(7245.0*coeff[0]*qr[7]+7245.0*coeff[0]*ql[7]+15750.0*coeff[0]*qc[7]-4364.768035073571*coeff[0]*qr[5]+4364.768035073571*coeff[0]*ql[5])*Jfac; 

  return 0.;

}

GKYL_CU_DH double dg_diffusion_vlasov_order6_surfx_1x1v_ser_p2_varcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],6.);

  out[0] += 0.0078125*((9736.860376938761*coeff[2]+7246.369435793345*coeff[1]+3187.575881449729*coeff[0])*qr[4]+(9736.860376938761*coeff[2]-7246.369435793345*coeff[1]+3187.575881449729*coeff[0])*ql[4]+(7254.915574973979*coeff[2]-6375.151762899458*coeff[0])*qc[4]+((-3697.127263159871*qr[1])+3697.127263159871*ql[1]+1138.419957660617*qr[0]+1138.419957660617*ql[0]-2276.839915321235*qc[0])*coeff[2]+((-5791.204537917824*coeff[1])-3086.357075906804*coeff[0])*qr[1]+(3086.357075906804*coeff[0]-5791.204537917824*coeff[1])*ql[1]-9800.49998724555*coeff[1]*qc[1]+(3086.357075906804*qr[0]-3086.357075906804*ql[0])*coeff[1]+1781.909088590101*coeff[0]*qr[0]+1781.909088590101*coeff[0]*ql[0]-3563.818177180201*coeff[0]*qc[0])*Jfac; 
  out[1] += 0.00390625*((43980.58833167195*coeff[2]+25756.75154207146*coeff[1]+9316.760703162869*coeff[0])*qr[4]+((-43980.58833167195*coeff[2])+25756.75154207146*coeff[1]-9316.760703162869*coeff[0])*ql[4]-39560.09352870644*coeff[1]*qc[4]+((-28887.40642563815*qr[1])-28887.40642563815*ql[1]-37852.46359221552*qc[1]+13802.60844913019*qr[0]-13802.60844913019*ql[0])*coeff[2]+((-23698.81326142724*coeff[1])-10245.97725939307*coeff[0])*qr[1]+(23698.81326142724*coeff[1]-10245.97725939307*coeff[0])*ql[1]-22273.86360737625*coeff[0]*qc[1]+(13237.03894381218*qr[0]+13237.03894381218*ql[0]-26474.07788762436*qc[0])*coeff[1]+6172.71415181361*coeff[0]*qr[0]-6172.71415181361*coeff[0]*ql[0])*Jfac; 
  out[2] += 0.0078125*((9736.86037693876*coeff[2]+7246.369435793348*coeff[1]+3187.575881449728*coeff[0])*qr[6]+(9736.86037693876*coeff[2]-7246.369435793348*coeff[1]+3187.575881449728*coeff[0])*ql[6]+(7254.915574973978*coeff[2]-6375.151762899456*coeff[0])*qc[6]+((-3697.127263159871*coeff[2])-5791.204537917824*coeff[1]-3086.357075906804*coeff[0])*qr[3]+(3697.127263159871*coeff[2]-5791.204537917824*coeff[1]+3086.357075906804*coeff[0])*ql[3]-9800.49998724555*coeff[1]*qc[3]+(1138.419957660617*coeff[2]+3086.357075906804*coeff[1]+1781.909088590101*coeff[0])*qr[2]+(1138.419957660617*coeff[2]-3086.357075906804*coeff[1]+1781.909088590101*coeff[0])*ql[2]+((-2276.839915321235*coeff[2])-3563.818177180201*coeff[0])*qc[2])*Jfac; 
  out[3] += 0.00390625*((43980.58833167197*coeff[2]+25756.75154207145*coeff[1]+9316.760703162876*coeff[0])*qr[6]+((-43980.58833167197*coeff[2])+25756.75154207145*coeff[1]-9316.760703162876*coeff[0])*ql[6]-39560.09352870643*coeff[1]*qc[6]+((-28887.40642563815*coeff[2])-23698.81326142724*coeff[1]-10245.97725939307*coeff[0])*qr[3]+((-28887.40642563815*coeff[2])+23698.81326142724*coeff[1]-10245.97725939307*coeff[0])*ql[3]+((-37852.46359221552*coeff[2])-22273.86360737625*coeff[0])*qc[3]+(13802.60844913019*coeff[2]+13237.03894381218*coeff[1]+6172.71415181361*coeff[0])*qr[2]+((-13802.60844913019*coeff[2])+13237.03894381218*coeff[1]-6172.71415181361*coeff[0])*ql[2]-26474.07788762436*coeff[1]*qc[2])*Jfac; 
  out[4] += 0.00390625*((65032.24008136275*coeff[2]+20832.91026237092*coeff[1]-572.7564927611036*coeff[0])*qr[4]+(65032.24008136275*coeff[2]-20832.91026237092*coeff[1]-572.7564927611036*coeff[0])*ql[4]+((-70297.4323855431*coeff[2])-25583.12334332929*coeff[0])*qc[4]+((-56766.92478900014*qr[1])+56766.92478900014*ql[1]+30547.01294725887*qr[0]+30547.01294725887*ql[0]-61094.02589451776*qc[0])*coeff[2]+((-22910.70164791991*coeff[1])-2218.276357895922*coeff[0])*qr[1]+(2218.276357895922*coeff[0]-22910.70164791991*coeff[1])*ql[1]-49805.873147652*coeff[1]*qc[1]+(13802.60844913019*qr[0]-13802.60844913019*ql[0])*coeff[1]+2276.839915321235*coeff[0]*qr[0]+2276.839915321235*coeff[0]*ql[0]-4553.67983064247*coeff[0]*qc[0])*Jfac; 
  out[5] += -0.0078125*((3697.127263159873*coeff[2]+5791.204537917824*coeff[1]+3086.357075906806*coeff[0])*qr[7]+((-3697.127263159873*coeff[2])+5791.204537917824*coeff[1]-3086.357075906806*coeff[0])*ql[7]+9800.499987245548*coeff[1]*qc[7]+((-1138.419957660617*coeff[2])-3086.357075906804*coeff[1]-1781.909088590101*coeff[0])*qr[5]+((-1138.419957660617*coeff[2])+3086.357075906804*coeff[1]-1781.909088590101*coeff[0])*ql[5]+(2276.839915321235*coeff[2]+3563.818177180201*coeff[0])*qc[5])*Jfac; 
  out[6] += 0.00390625*((65032.24008136275*coeff[2]+20832.91026237092*coeff[1]-572.7564927611036*coeff[0])*qr[6]+(65032.24008136275*coeff[2]-20832.91026237092*coeff[1]-572.7564927611036*coeff[0])*ql[6]+((-70297.4323855431*coeff[2])-25583.12334332929*coeff[0])*qc[6]+((-56766.92478900017*coeff[2])-22910.70164791991*coeff[1]-2218.276357895923*coeff[0])*qr[3]+(56766.92478900017*coeff[2]-22910.70164791991*coeff[1]+2218.276357895923*coeff[0])*ql[3]-49805.87314765198*coeff[1]*qc[3]+(30547.01294725887*coeff[2]+13802.6084491302*coeff[1]+2276.839915321234*coeff[0])*qr[2]+(30547.01294725887*coeff[2]-13802.6084491302*coeff[1]+2276.839915321234*coeff[0])*ql[2]+((-61094.02589451776*coeff[2])-4553.679830642469*coeff[0])*qc[2])*Jfac; 
  out[7] += -0.00390625*((28887.40642563815*coeff[2]+23698.81326142724*coeff[1]+10245.97725939307*coeff[0])*qr[7]+(28887.40642563815*coeff[2]-23698.81326142724*coeff[1]+10245.97725939307*coeff[0])*ql[7]+(37852.46359221552*coeff[2]+22273.86360737625*coeff[0])*qc[7]+((-13802.6084491302*coeff[2])-13237.03894381217*coeff[1]-6172.714151813613*coeff[0])*qr[5]+(13802.6084491302*coeff[2]-13237.03894381217*coeff[1]+6172.714151813613*coeff[0])*ql[5]+26474.07788762436*coeff[1]*qc[5])*Jfac; 

  return 0.;

}

