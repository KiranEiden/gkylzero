#include <gkyl_dg_diffusion_kernels.h>

GKYL_CU_DH double dg_diffusion_order6_vlasov_surfx_1x2v_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],6.);

  out[0] += 0.0625*(563.489130329947*coeff[0]*qr[7]+563.489130329947*coeff[0]*ql[7]-1126.978260659894*coeff[0]*qc[7]-545.5960043841961*coeff[0]*qr[1]+545.5960043841961*coeff[0]*ql[1]+315.0*coeff[0]*qr[0]+315.0*coeff[0]*ql[0]-630.0*coeff[0]*qc[0])*Jfac; 
  out[1] += 0.0078125*(6587.944671898812*coeff[0]*qr[7]-6587.944671898812*coeff[0]*ql[7]-7245.0*coeff[0]*qr[1]-7245.0*coeff[0]*ql[1]-15750.0*coeff[0]*qc[1]+4364.768035073569*coeff[0]*qr[0]-4364.768035073569*coeff[0]*ql[0])*Jfac; 
  out[2] += 0.0625*(563.4891303299469*coeff[0]*qr[11]+563.4891303299469*coeff[0]*ql[11]-1126.978260659894*coeff[0]*qc[11]-545.5960043841961*coeff[0]*qr[4]+545.5960043841961*coeff[0]*ql[4]+315.0*coeff[0]*qr[2]+315.0*coeff[0]*ql[2]-630.0*coeff[0]*qc[2])*Jfac; 
  out[3] += 0.0625*(563.4891303299469*coeff[0]*qr[13]+563.4891303299469*coeff[0]*ql[13]-1126.978260659894*coeff[0]*qc[13]-545.5960043841961*coeff[0]*qr[5]+545.5960043841961*coeff[0]*ql[5]+315.0*coeff[0]*qr[3]+315.0*coeff[0]*ql[3]-630.0*coeff[0]*qc[3])*Jfac; 
  out[4] += 0.0078125*(6587.944671898817*coeff[0]*qr[11]-6587.944671898817*coeff[0]*ql[11]-7245.0*coeff[0]*qr[4]-7245.0*coeff[0]*ql[4]-15750.0*coeff[0]*qc[4]+4364.768035073569*coeff[0]*qr[2]-4364.768035073569*coeff[0]*ql[2])*Jfac; 
  out[5] += 0.0078125*(6587.944671898817*coeff[0]*qr[13]-6587.944671898817*coeff[0]*ql[13]-7245.0*coeff[0]*qr[5]-7245.0*coeff[0]*ql[5]-15750.0*coeff[0]*qc[5]+4364.768035073569*coeff[0]*qr[3]-4364.768035073569*coeff[0]*ql[3])*Jfac; 
  out[6] += 0.0625*(563.489130329947*coeff[0]*qr[17]+563.489130329947*coeff[0]*ql[17]-1126.978260659894*coeff[0]*qc[17]-545.5960043841961*coeff[0]*qr[10]+545.5960043841961*coeff[0]*ql[10]+315.0*coeff[0]*qr[6]+315.0*coeff[0]*ql[6]-630.0*coeff[0]*qc[6])*Jfac; 
  out[7] += -0.0078125*(405.0*coeff[0]*qr[7]+405.0*coeff[0]*ql[7]+18090.0*coeff[0]*qc[7]+1568.558255214003*coeff[0]*qr[1]-1568.558255214003*coeff[0]*ql[1]-1609.968943799849*coeff[0]*qr[0]-1609.968943799849*coeff[0]*ql[0]+3219.937887599698*coeff[0]*qc[0])*Jfac; 
  out[8] += -0.0625*(545.5960043841964*coeff[0]*qr[12]-545.5960043841964*coeff[0]*ql[12]-315.0*coeff[0]*qr[8]-315.0*coeff[0]*ql[8]+630.0*coeff[0]*qc[8])*Jfac; 
  out[9] += -0.0625*(545.5960043841964*coeff[0]*qr[15]-545.5960043841964*coeff[0]*ql[15]-315.0*coeff[0]*qr[9]-315.0*coeff[0]*ql[9]+630.0*coeff[0]*qc[9])*Jfac; 
  out[10] += 0.0078125*(6587.944671898812*coeff[0]*qr[17]-6587.944671898812*coeff[0]*ql[17]-7245.0*coeff[0]*qr[10]-7245.0*coeff[0]*ql[10]-15750.0*coeff[0]*qc[10]+4364.768035073569*coeff[0]*qr[6]-4364.768035073569*coeff[0]*ql[6])*Jfac; 
  out[11] += -0.0078125*(405.0*coeff[0]*qr[11]+405.0*coeff[0]*ql[11]+18090.0*coeff[0]*qc[11]+1568.558255214004*coeff[0]*qr[4]-1568.558255214004*coeff[0]*ql[4]-1609.968943799848*coeff[0]*qr[2]-1609.968943799848*coeff[0]*ql[2]+3219.937887599697*coeff[0]*qc[2])*Jfac; 
  out[12] += -0.0078125*(7245.0*coeff[0]*qr[12]+7245.0*coeff[0]*ql[12]+15750.0*coeff[0]*qc[12]-4364.768035073571*coeff[0]*qr[8]+4364.768035073571*coeff[0]*ql[8])*Jfac; 
  out[13] += -0.0078125*(405.0*coeff[0]*qr[13]+405.0*coeff[0]*ql[13]+18090.0*coeff[0]*qc[13]+1568.558255214004*coeff[0]*qr[5]-1568.558255214004*coeff[0]*ql[5]-1609.968943799848*coeff[0]*qr[3]-1609.968943799848*coeff[0]*ql[3]+3219.937887599697*coeff[0]*qc[3])*Jfac; 
  out[14] += -0.0625*(545.5960043841964*coeff[0]*qr[18]-545.5960043841964*coeff[0]*ql[18]-315.0*coeff[0]*qr[14]-315.0*coeff[0]*ql[14]+630.0*coeff[0]*qc[14])*Jfac; 
  out[15] += -0.0078125*(7245.0*coeff[0]*qr[15]+7245.0*coeff[0]*ql[15]+15750.0*coeff[0]*qc[15]-4364.768035073571*coeff[0]*qr[9]+4364.768035073571*coeff[0]*ql[9])*Jfac; 
  out[16] += -0.0625*(545.5960043841964*coeff[0]*qr[19]-545.5960043841964*coeff[0]*ql[19]-315.0*coeff[0]*qr[16]-315.0*coeff[0]*ql[16]+630.0*coeff[0]*qc[16])*Jfac; 
  out[17] += -0.0078125*(405.0*coeff[0]*qr[17]+405.0*coeff[0]*ql[17]+18090.0*coeff[0]*qc[17]+1568.558255214003*coeff[0]*qr[10]-1568.558255214003*coeff[0]*ql[10]-1609.968943799849*coeff[0]*qr[6]-1609.968943799849*coeff[0]*ql[6]+3219.937887599698*coeff[0]*qc[6])*Jfac; 
  out[18] += -0.0078125*(7245.0*coeff[0]*qr[18]+7245.0*coeff[0]*ql[18]+15750.0*coeff[0]*qc[18]-4364.768035073571*coeff[0]*qr[14]+4364.768035073571*coeff[0]*ql[14])*Jfac; 
  out[19] += -0.0078125*(7245.0*coeff[0]*qr[19]+7245.0*coeff[0]*ql[19]+15750.0*coeff[0]*qc[19]-4364.768035073571*coeff[0]*qr[16]+4364.768035073571*coeff[0]*ql[16])*Jfac; 

  return 0.;

}

GKYL_CU_DH double dg_diffusion_order6_vlasov_surfx_1x2v_ser_p2_varcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],6.);

  out[0] += 0.0110485434560398*((6885.0*coeff[2]+5123.95696703241*coeff[1]+2253.956521319788*coeff[0])*qr[7]+(6885.0*coeff[2]-5123.95696703241*coeff[1]+2253.956521319788*coeff[0])*ql[7]+(5130.0*coeff[2]-4507.913042639576*coeff[0])*qc[7]+((-2614.263758690006*qr[1])+2614.263758690006*ql[1]+804.9844718999244*qr[0]+804.9844718999244*ql[0]-1609.968943799849*qc[0])*coeff[2]+((-4095.0*coeff[1])-2182.384017536785*coeff[0])*qr[1]+(2182.384017536785*coeff[0]-4095.0*coeff[1])*ql[1]-6930.0*coeff[1]*qc[1]+(2182.384017536785*qr[0]-2182.384017536785*ql[0])*coeff[1]+1260.0*coeff[0]*qr[0]+1260.0*coeff[0]*ql[0]-2520.0*coeff[0]*qc[0])*Jfac; 
  out[1] += 0.005524271728019897*((31098.97224989918*coeff[2]+18212.77367673579*coeff[1]+6587.944671898812*coeff[0])*qr[7]+((-31098.97224989918*coeff[2])+18212.77367673579*coeff[1]-6587.944671898812*coeff[0])*ql[7]-27973.21039852238*coeff[1]*qc[7]+((-20426.48097446058*qr[1])-20426.48097446058*ql[1]-26765.73369067249*qc[1]+9759.918032442689*qr[0]-9759.918032442689*ql[0])*coeff[2]+((-16757.59156322888*coeff[1])-7245.0*coeff[0])*qr[1]+(16757.59156322888*coeff[1]-7245.0*coeff[0])*ql[1]-15750.0*coeff[0]*qc[1]+(9360.0*qr[0]+9360.0*ql[0]-18720.0*qc[0])*coeff[1]+4364.768035073569*coeff[0]*qr[0]-4364.768035073569*coeff[0]*ql[0])*Jfac; 
  out[2] += 0.0110485434560398*((6884.999999999997*coeff[2]+5123.956967032413*coeff[1]+2253.956521319787*coeff[0])*qr[11]+(6884.999999999997*coeff[2]-5123.956967032413*coeff[1]+2253.956521319787*coeff[0])*ql[11]+(5129.999999999999*coeff[2]-4507.913042639575*coeff[0])*qc[11]+((-2614.263758690006*coeff[2])-4095.0*coeff[1]-2182.384017536785*coeff[0])*qr[4]+(2614.263758690006*coeff[2]-4095.0*coeff[1]+2182.384017536785*coeff[0])*ql[4]-6930.0*coeff[1]*qc[4]+(804.9844718999244*coeff[2]+2182.384017536785*coeff[1]+1260.0*coeff[0])*qr[2]+(804.9844718999244*coeff[2]-2182.384017536785*coeff[1]+1260.0*coeff[0])*ql[2]+((-1609.968943799849*coeff[2])-2520.0*coeff[0])*qc[2])*Jfac; 
  out[3] += 0.0110485434560398*((6884.999999999997*coeff[2]+5123.956967032413*coeff[1]+2253.956521319787*coeff[0])*qr[13]+(6884.999999999997*coeff[2]-5123.956967032413*coeff[1]+2253.956521319787*coeff[0])*ql[13]+(5129.999999999999*coeff[2]-4507.913042639575*coeff[0])*qc[13]+((-2614.263758690006*coeff[2])-4095.0*coeff[1]-2182.384017536785*coeff[0])*qr[5]+(2614.263758690006*coeff[2]-4095.0*coeff[1]+2182.384017536785*coeff[0])*ql[5]-6930.0*coeff[1]*qc[5]+(804.9844718999244*coeff[2]+2182.384017536785*coeff[1]+1260.0*coeff[0])*qr[3]+(804.9844718999244*coeff[2]-2182.384017536785*coeff[1]+1260.0*coeff[0])*ql[3]+((-1609.968943799849*coeff[2])-2520.0*coeff[0])*qc[3])*Jfac; 
  out[4] += 0.005524271728019897*((31098.97224989919*coeff[2]+18212.77367673578*coeff[1]+6587.944671898817*coeff[0])*qr[11]+((-31098.97224989919*coeff[2])+18212.77367673578*coeff[1]-6587.944671898817*coeff[0])*ql[11]-27973.21039852237*coeff[1]*qc[11]+((-20426.48097446058*coeff[2])-16757.59156322888*coeff[1]-7245.0*coeff[0])*qr[4]+((-20426.48097446058*coeff[2])+16757.59156322888*coeff[1]-7245.0*coeff[0])*ql[4]+((-26765.73369067249*coeff[2])-15750.0*coeff[0])*qc[4]+(9759.918032442689*coeff[2]+9360.0*coeff[1]+4364.768035073569*coeff[0])*qr[2]+((-9759.918032442689*coeff[2])+9360.0*coeff[1]-4364.768035073569*coeff[0])*ql[2]-18720.0*coeff[1]*qc[2])*Jfac; 
  out[5] += 0.005524271728019897*((31098.97224989919*coeff[2]+18212.77367673578*coeff[1]+6587.944671898817*coeff[0])*qr[13]+((-31098.97224989919*coeff[2])+18212.77367673578*coeff[1]-6587.944671898817*coeff[0])*ql[13]-27973.21039852237*coeff[1]*qc[13]+((-20426.48097446058*coeff[2])-16757.59156322888*coeff[1]-7245.0*coeff[0])*qr[5]+((-20426.48097446058*coeff[2])+16757.59156322888*coeff[1]-7245.0*coeff[0])*ql[5]+((-26765.73369067249*coeff[2])-15750.0*coeff[0])*qc[5]+(9759.918032442689*coeff[2]+9360.0*coeff[1]+4364.768035073569*coeff[0])*qr[3]+((-9759.918032442689*coeff[2])+9360.0*coeff[1]-4364.768035073569*coeff[0])*ql[3]-18720.0*coeff[1]*qc[3])*Jfac; 
  out[6] += 0.0110485434560398*((6885.0*coeff[2]+5123.95696703241*coeff[1]+2253.956521319788*coeff[0])*qr[17]+(6885.0*coeff[2]-5123.95696703241*coeff[1]+2253.956521319788*coeff[0])*ql[17]+(5130.0*coeff[2]-4507.913042639576*coeff[0])*qc[17]+((-2614.263758690006*coeff[2])-4095.0*coeff[1]-2182.384017536785*coeff[0])*qr[10]+(2614.263758690006*coeff[2]-4095.0*coeff[1]+2182.384017536785*coeff[0])*ql[10]-6930.0*coeff[1]*qc[10]+(804.9844718999244*coeff[2]+2182.384017536785*coeff[1]+1260.0*coeff[0])*qr[6]+(804.9844718999244*coeff[2]-2182.384017536785*coeff[1]+1260.0*coeff[0])*ql[6]+((-1609.968943799849*coeff[2])-2520.0*coeff[0])*qc[6])*Jfac; 
  out[7] += 0.005524271728019897*((45984.73795728319*coeff[2]+14731.09211837329*coeff[1]-405.0*coeff[0])*qr[7]+(45984.73795728319*coeff[2]-14731.09211837329*coeff[1]-405.0*coeff[0])*ql[7]+((-49707.79113982034*coeff[2])-18090.0*coeff[0])*qc[7]+((-40140.27746540872*qr[1])+40140.27746540872*ql[1]+21600.0*qr[0]+21600.0*ql[0]-43200.0*qc[0])*coeff[2]+((-16200.31249698598*coeff[1])-1568.558255214003*coeff[0])*qr[1]+(1568.558255214003*coeff[0]-16200.31249698598*coeff[1])*ql[1]-35218.07064562169*coeff[1]*qc[1]+(9759.918032442689*qr[0]-9759.918032442689*ql[0])*coeff[1]+1609.968943799849*coeff[0]*qr[0]+1609.968943799849*coeff[0]*ql[0]-3219.937887599698*coeff[0]*qc[0])*Jfac; 
  out[8] += -0.0110485434560398*((2614.263758690007*coeff[2]+4095.0*coeff[1]+2182.384017536785*coeff[0])*qr[12]+((-2614.263758690007*coeff[2])+4095.0*coeff[1]-2182.384017536785*coeff[0])*ql[12]+6930.0*coeff[1]*qc[12]+((-804.9844718999244*coeff[2])-2182.384017536785*coeff[1]-1260.0*coeff[0])*qr[8]+((-804.9844718999244*coeff[2])+2182.384017536785*coeff[1]-1260.0*coeff[0])*ql[8]+(1609.968943799849*coeff[2]+2520.0*coeff[0])*qc[8])*Jfac; 
  out[9] += -0.0110485434560398*((2614.263758690007*coeff[2]+4095.0*coeff[1]+2182.384017536785*coeff[0])*qr[15]+((-2614.263758690007*coeff[2])+4095.0*coeff[1]-2182.384017536785*coeff[0])*ql[15]+6930.0*coeff[1]*qc[15]+((-804.9844718999244*coeff[2])-2182.384017536785*coeff[1]-1260.0*coeff[0])*qr[9]+((-804.9844718999244*coeff[2])+2182.384017536785*coeff[1]-1260.0*coeff[0])*ql[9]+(1609.968943799849*coeff[2]+2520.0*coeff[0])*qc[9])*Jfac; 
  out[10] += 0.005524271728019897*((31098.97224989918*coeff[2]+18212.77367673579*coeff[1]+6587.944671898812*coeff[0])*qr[17]+((-31098.97224989918*coeff[2])+18212.77367673579*coeff[1]-6587.944671898812*coeff[0])*ql[17]-27973.21039852238*coeff[1]*qc[17]+((-20426.48097446058*coeff[2])-16757.59156322888*coeff[1]-7245.0*coeff[0])*qr[10]+((-20426.48097446058*coeff[2])+16757.59156322888*coeff[1]-7245.0*coeff[0])*ql[10]+((-26765.73369067249*coeff[2])-15750.0*coeff[0])*qc[10]+(9759.918032442689*coeff[2]+9360.0*coeff[1]+4364.768035073569*coeff[0])*qr[6]+((-9759.918032442689*coeff[2])+9360.0*coeff[1]-4364.768035073569*coeff[0])*ql[6]-18720.0*coeff[1]*qc[6])*Jfac; 
  out[11] += 0.005524271728019897*((45984.73795728319*coeff[2]+14731.09211837329*coeff[1]-405.0*coeff[0])*qr[11]+(45984.73795728319*coeff[2]-14731.09211837329*coeff[1]-405.0*coeff[0])*ql[11]+((-49707.79113982034*coeff[2])-18090.0*coeff[0])*qc[11]+((-40140.27746540874*coeff[2])-16200.31249698597*coeff[1]-1568.558255214004*coeff[0])*qr[4]+(40140.27746540874*coeff[2]-16200.31249698597*coeff[1]+1568.558255214004*coeff[0])*ql[4]-35218.07064562168*coeff[1]*qc[4]+(21600.0*coeff[2]+9759.918032442692*coeff[1]+1609.968943799848*coeff[0])*qr[2]+(21600.0*coeff[2]-9759.918032442692*coeff[1]+1609.968943799848*coeff[0])*ql[2]+((-43199.99999999999*coeff[2])-3219.937887599697*coeff[0])*qc[2])*Jfac; 
  out[12] += -0.005524271728019897*((20426.48097446058*coeff[2]+16757.59156322888*coeff[1]+7245.0*coeff[0])*qr[12]+(20426.48097446058*coeff[2]-16757.59156322888*coeff[1]+7245.0*coeff[0])*ql[12]+(26765.73369067249*coeff[2]+15750.0*coeff[0])*qc[12]+((-9759.918032442692*coeff[2])-9360.0*coeff[1]-4364.768035073571*coeff[0])*qr[8]+(9759.918032442692*coeff[2]-9360.0*coeff[1]+4364.768035073571*coeff[0])*ql[8]+18720.0*coeff[1]*qc[8])*Jfac; 
  out[13] += 0.005524271728019897*((45984.73795728319*coeff[2]+14731.09211837329*coeff[1]-405.0*coeff[0])*qr[13]+(45984.73795728319*coeff[2]-14731.09211837329*coeff[1]-405.0*coeff[0])*ql[13]+((-49707.79113982034*coeff[2])-18090.0*coeff[0])*qc[13]+((-40140.27746540874*coeff[2])-16200.31249698597*coeff[1]-1568.558255214004*coeff[0])*qr[5]+(40140.27746540874*coeff[2]-16200.31249698597*coeff[1]+1568.558255214004*coeff[0])*ql[5]-35218.07064562168*coeff[1]*qc[5]+(21600.0*coeff[2]+9759.918032442692*coeff[1]+1609.968943799848*coeff[0])*qr[3]+(21600.0*coeff[2]-9759.918032442692*coeff[1]+1609.968943799848*coeff[0])*ql[3]+((-43199.99999999999*coeff[2])-3219.937887599697*coeff[0])*qc[3])*Jfac; 
  out[14] += -0.0110485434560398*((2614.263758690007*coeff[2]+4095.0*coeff[1]+2182.384017536785*coeff[0])*qr[18]+((-2614.263758690007*coeff[2])+4095.0*coeff[1]-2182.384017536785*coeff[0])*ql[18]+6930.0*coeff[1]*qc[18]+((-804.9844718999244*coeff[2])-2182.384017536785*coeff[1]-1260.0*coeff[0])*qr[14]+((-804.9844718999244*coeff[2])+2182.384017536785*coeff[1]-1260.0*coeff[0])*ql[14]+(1609.968943799849*coeff[2]+2520.0*coeff[0])*qc[14])*Jfac; 
  out[15] += -0.005524271728019897*((20426.48097446058*coeff[2]+16757.59156322888*coeff[1]+7245.0*coeff[0])*qr[15]+(20426.48097446058*coeff[2]-16757.59156322888*coeff[1]+7245.0*coeff[0])*ql[15]+(26765.73369067249*coeff[2]+15750.0*coeff[0])*qc[15]+((-9759.918032442692*coeff[2])-9360.0*coeff[1]-4364.768035073571*coeff[0])*qr[9]+(9759.918032442692*coeff[2]-9360.0*coeff[1]+4364.768035073571*coeff[0])*ql[9]+18720.0*coeff[1]*qc[9])*Jfac; 
  out[16] += -0.0110485434560398*((2614.263758690007*coeff[2]+4095.0*coeff[1]+2182.384017536785*coeff[0])*qr[19]+((-2614.263758690007*coeff[2])+4095.0*coeff[1]-2182.384017536785*coeff[0])*ql[19]+6930.0*coeff[1]*qc[19]+((-804.9844718999244*coeff[2])-2182.384017536785*coeff[1]-1260.0*coeff[0])*qr[16]+((-804.9844718999244*coeff[2])+2182.384017536785*coeff[1]-1260.0*coeff[0])*ql[16]+(1609.968943799849*coeff[2]+2520.0*coeff[0])*qc[16])*Jfac; 
  out[17] += 0.005524271728019897*((45984.73795728319*coeff[2]+14731.09211837329*coeff[1]-405.0*coeff[0])*qr[17]+(45984.73795728319*coeff[2]-14731.09211837329*coeff[1]-405.0*coeff[0])*ql[17]+((-49707.79113982034*coeff[2])-18090.0*coeff[0])*qc[17]+((-40140.27746540872*coeff[2])-16200.31249698598*coeff[1]-1568.558255214003*coeff[0])*qr[10]+(40140.27746540872*coeff[2]-16200.31249698598*coeff[1]+1568.558255214003*coeff[0])*ql[10]-35218.07064562169*coeff[1]*qc[10]+(21600.0*coeff[2]+9759.918032442689*coeff[1]+1609.968943799849*coeff[0])*qr[6]+(21600.0*coeff[2]-9759.918032442689*coeff[1]+1609.968943799849*coeff[0])*ql[6]+((-43200.0*coeff[2])-3219.937887599698*coeff[0])*qc[6])*Jfac; 
  out[18] += -0.005524271728019897*((20426.48097446058*coeff[2]+16757.59156322888*coeff[1]+7245.0*coeff[0])*qr[18]+(20426.48097446058*coeff[2]-16757.59156322888*coeff[1]+7245.0*coeff[0])*ql[18]+(26765.73369067249*coeff[2]+15750.0*coeff[0])*qc[18]+((-9759.918032442692*coeff[2])-9360.0*coeff[1]-4364.768035073571*coeff[0])*qr[14]+(9759.918032442692*coeff[2]-9360.0*coeff[1]+4364.768035073571*coeff[0])*ql[14]+18720.0*coeff[1]*qc[14])*Jfac; 
  out[19] += -0.005524271728019897*((20426.48097446058*coeff[2]+16757.59156322888*coeff[1]+7245.0*coeff[0])*qr[19]+(20426.48097446058*coeff[2]-16757.59156322888*coeff[1]+7245.0*coeff[0])*ql[19]+(26765.73369067249*coeff[2]+15750.0*coeff[0])*qc[19]+((-9759.918032442692*coeff[2])-9360.0*coeff[1]-4364.768035073571*coeff[0])*qr[16]+(9759.918032442692*coeff[2]-9360.0*coeff[1]+4364.768035073571*coeff[0])*ql[16]+18720.0*coeff[1]*qc[16])*Jfac; 

  return 0.;

}

