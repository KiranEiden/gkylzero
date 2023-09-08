#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order4_surfx_2x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],4.);

  out[0] += 0.0625*(107.3312629199899*coeff[0]*qr[4]+107.3312629199899*coeff[0]*ql[4]-214.6625258399798*coeff[0]*qc[4]-129.9038105676658*coeff[0]*qr[1]+129.9038105676658*coeff[0]*ql[1]+75.0*coeff[0]*qr[0]+75.0*coeff[0]*ql[0]-150.0*coeff[0]*qc[0])*Jfac; 
  out[1] += 0.03125*(290.4737509655563*coeff[0]*qr[4]-290.4737509655563*coeff[0]*ql[4]-405.0*coeff[0]*qr[1]-405.0*coeff[0]*ql[1]-990.0*coeff[0]*qc[1]+259.8076211353315*coeff[0]*qr[0]-259.8076211353315*coeff[0]*ql[0])*Jfac; 
  out[2] += 0.0625*(107.3312629199899*coeff[0]*qr[6]+107.3312629199899*coeff[0]*ql[6]-214.6625258399798*coeff[0]*qc[6]-129.9038105676658*coeff[0]*qr[3]+129.9038105676658*coeff[0]*ql[3]+75.0*coeff[0]*qr[2]+75.0*coeff[0]*ql[2]-150.0*coeff[0]*qc[2])*Jfac; 
  out[3] += 0.03125*(290.4737509655563*coeff[0]*qr[6]-290.4737509655563*coeff[0]*ql[6]-405.0*coeff[0]*qr[3]-405.0*coeff[0]*ql[3]-990.0*coeff[0]*qc[3]+259.8076211353315*coeff[0]*qr[2]-259.8076211353315*coeff[0]*ql[2])*Jfac; 
  out[4] += 0.03125*(21.0*coeff[0]*qr[4]+21.0*coeff[0]*ql[4]-1302.0*coeff[0]*qc[4]-151.0463505020892*coeff[0]*qr[1]+151.0463505020892*coeff[0]*ql[1]+134.1640786499874*coeff[0]*qr[0]+134.1640786499874*coeff[0]*ql[0]-268.3281572999748*coeff[0]*qc[0])*Jfac; 
  out[5] += -0.0625*(129.9038105676658*coeff[0]*qr[7]-129.9038105676658*coeff[0]*ql[7]-75.0*coeff[0]*qr[5]-75.0*coeff[0]*ql[5]+150.0*coeff[0]*qc[5])*Jfac; 
  out[6] += 0.03125*(21.0*coeff[0]*qr[6]+21.0*coeff[0]*ql[6]-1302.0*coeff[0]*qc[6]-151.0463505020893*coeff[0]*qr[3]+151.0463505020893*coeff[0]*ql[3]+134.1640786499874*coeff[0]*qr[2]+134.1640786499874*coeff[0]*ql[2]-268.3281572999747*coeff[0]*qc[2])*Jfac; 
  out[7] += -0.03125*(405.0*coeff[0]*qr[7]+405.0*coeff[0]*ql[7]+990.0*coeff[0]*qc[7]-259.8076211353317*coeff[0]*qr[5]+259.8076211353317*coeff[0]*ql[5])*Jfac; 

  return 0.;

}

GKYL_CU_DH double dg_diffusion_fluid_order4_surfx_2x_ser_p2_varcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],4.);

  out[0] += -0.015625*((495.0*coeff[7]+259.8076211353317*coeff[5])*qr[7]+(495.0*coeff[7]-259.8076211353317*coeff[5])*ql[7]+810.0*coeff[7]*qc[7]+(259.8076211353317*ql[5]-259.8076211353317*qr[5])*coeff[7]+((-795.0*coeff[6])-453.1390515062678*coeff[3]-214.6625258399798*coeff[2])*qr[6]+((-795.0*coeff[6])+453.1390515062678*coeff[3]-214.6625258399798*coeff[2])*ql[6]+(330.0*coeff[6]+429.3250516799596*coeff[2])*qc[6]+(755.2317525104464*qr[3]-755.2317525104464*ql[3]-335.4101966249685*qr[2]-335.4101966249685*ql[2]+670.820393249937*qc[2])*coeff[6]-150.0*coeff[5]*qr[5]-150.0*coeff[5]*ql[5]+300.0*coeff[5]*qc[5]+((-795.0*coeff[4])-453.1390515062677*coeff[1]-214.6625258399798*coeff[0])*qr[4]+((-795.0*coeff[4])+453.1390515062677*coeff[1]-214.6625258399798*coeff[0])*ql[4]+(330.0*coeff[4]+429.3250516799596*coeff[0])*qc[4]+(755.2317525104463*qr[1]-755.2317525104463*ql[1]-335.4101966249685*qr[0]-335.4101966249685*ql[0]+670.8203932499371*qc[0])*coeff[4]+(495.0*coeff[3]+259.8076211353315*coeff[2])*qr[3]+(495.0*coeff[3]-259.8076211353315*coeff[2])*ql[3]+810.0*coeff[3]*qc[3]+(259.8076211353315*ql[2]-259.8076211353315*qr[2])*coeff[3]-150.0*coeff[2]*qr[2]-150.0*coeff[2]*ql[2]+300.0*coeff[2]*qc[2]+(495.0*coeff[1]+259.8076211353315*coeff[0])*qr[1]+(495.0*coeff[1]-259.8076211353315*coeff[0])*ql[1]+810.0*coeff[1]*qc[1]+(259.8076211353315*ql[0]-259.8076211353315*qr[0])*coeff[1]-150.0*coeff[0]*qr[0]-150.0*coeff[0]*ql[0]+300.0*coeff[0]*qc[0])*Jfac; 
  out[1] += -0.015625*((779.4228634059946*coeff[7]+404.9999999999999*coeff[5])*qr[7]+(404.9999999999999*coeff[5]-779.4228634059946*coeff[7])*ql[7]+990.0*coeff[5]*qc[7]+((-450.0000000000001*qr[5])-450.0000000000001*ql[5]+900.0000000000001*qc[5])*coeff[7]+((-1195.115057222525*coeff[6])-643.9875775199394*coeff[3]-290.4737509655563*coeff[2])*qr[6]+(1195.115057222525*coeff[6]-643.9875775199394*coeff[3]+290.4737509655563*coeff[2])*ql[6]+1287.975155039879*coeff[3]*qc[6]+(1207.476707849886*qr[3]+1207.476707849886*ql[3]+1609.968943799848*qc[3]-580.9475019311126*qr[2]+580.9475019311126*ql[2])*coeff[6]-259.8076211353315*coeff[5]*qr[5]+259.8076211353315*coeff[5]*ql[5]+((-1195.115057222525*coeff[4])-643.9875775199395*coeff[1]-290.4737509655563*coeff[0])*qr[4]+(1195.115057222525*coeff[4]-643.9875775199395*coeff[1]+290.4737509655563*coeff[0])*ql[4]+1287.975155039879*coeff[1]*qc[4]+(1207.476707849887*qr[1]+1207.476707849887*ql[1]+1609.968943799849*qc[1]-580.9475019311126*qr[0]+580.9475019311126*ql[0])*coeff[4]+(779.4228634059946*coeff[3]+405.0*coeff[2])*qr[3]+(405.0*coeff[2]-779.4228634059946*coeff[3])*ql[3]+990.0*coeff[2]*qc[3]+((-450.0*qr[2])-450.0*ql[2]+900.0*qc[2])*coeff[3]-259.8076211353315*coeff[2]*qr[2]+259.8076211353315*coeff[2]*ql[2]+(779.4228634059946*coeff[1]+405.0*coeff[0])*qr[1]+(405.0*coeff[0]-779.4228634059946*coeff[1])*ql[1]+990.0*coeff[0]*qc[1]+((-450.0*qr[0])-450.0*ql[0]+900.0*qc[0])*coeff[1]-259.8076211353315*coeff[0]*qr[0]+259.8076211353315*coeff[0]*ql[0])*Jfac; 
  out[2] += -0.003125*((3377.49907475931*coeff[6]+2213.707297724792*coeff[3]+1161.895003862225*coeff[2])*qr[7]+((-3377.49907475931*coeff[6])+2213.707297724792*coeff[3]-1161.895003862225*coeff[2])*ql[7]+3622.430123549658*coeff[3]*qc[7]+((-2026.499444855586*qr[6])+2026.499444855586*ql[6]+2213.707297724792*qr[3]+2213.707297724792*ql[3]+3622.430123549658*qc[3]-1161.895003862225*qr[2]+1161.895003862225*ql[2])*coeff[7]+((-960.0000000000001*coeff[5])-3975.000000000001*coeff[4]-2265.695257531339*coeff[1]-1073.312629199899*coeff[0])*qr[6]+((-960.0000000000001*coeff[5])-3975.000000000001*coeff[4]+2265.695257531339*coeff[1]-1073.312629199899*coeff[0])*ql[6]+(1920.0*coeff[5]+1650.0*coeff[4]+2146.625258399798*coeff[0])*qc[6]+((-1500.0*qr[5])-1500.0*ql[5]+3000.0*qc[5]-3975.000000000001*qr[4]-3975.000000000001*ql[4]+1650.0*qc[4]+3776.158762552232*qr[1]-3776.158762552232*ql[1]-1677.050983124842*qr[0]-1677.050983124842*ql[0]+3354.101966249685*qc[0])*coeff[6]+((-1161.895003862225*coeff[3])-670.8203932499371*coeff[2])*qr[5]+(1161.895003862225*coeff[3]-670.8203932499371*coeff[2])*ql[5]+1341.640786499874*coeff[2]*qc[5]+(1161.895003862225*qr[3]-1161.895003862225*ql[3]-670.8203932499371*qr[2]-670.8203932499371*ql[2]+1341.640786499874*qc[2])*coeff[5]+((-2265.695257531339*coeff[3])-1073.312629199899*coeff[2])*qr[4]+(2265.695257531339*coeff[3]-1073.312629199899*coeff[2])*ql[4]+2146.625258399798*coeff[2]*qc[4]+(3776.158762552232*qr[3]-3776.158762552232*ql[3]-1677.050983124843*qr[2]-1677.050983124843*ql[2]+3354.101966249686*qc[2])*coeff[4]+(2475.0*coeff[1]+1299.038105676658*coeff[0])*qr[3]+(2475.0*coeff[1]-1299.038105676658*coeff[0])*ql[3]+4050.0*coeff[1]*qc[3]+(2475.0*qr[1]+2475.0*ql[1]+4050.0*qc[1]-1299.038105676658*qr[0]+1299.038105676658*ql[0])*coeff[3]+((-1299.038105676658*coeff[1])-750.0*coeff[0])*qr[2]+(1299.038105676658*coeff[1]-750.0*coeff[0])*ql[2]+1500.0*coeff[0]*qc[2]+(1299.038105676658*qr[1]-1299.038105676658*ql[1]-750.0*qr[0]-750.0*ql[0]+1500.0*qc[0])*coeff[2])*Jfac; 
  out[3] += -0.015625*((1080.0*coeff[6]+697.1370023173351*coeff[3]+362.2430123549658*coeff[2])*qr[7]+(1080.0*coeff[6]-697.1370023173351*coeff[3]+362.2430123549658*coeff[2])*ql[7]+(1440.0*coeff[6]+885.4829190899167*coeff[2])*qc[7]+((-576.0*qr[6])-576.0*ql[6]+1152.0*qc[6]+697.1370023173351*qr[3]-697.1370023173351*ql[3]-402.4922359499621*qr[2]-402.4922359499621*ql[2]+804.9844718999242*qc[2])*coeff[7]+((-259.8076211353317*coeff[5])-1195.115057222525*coeff[4]-643.9875775199394*coeff[1]-290.4737509655563*coeff[0])*qr[6]+(259.8076211353317*coeff[5]+1195.115057222525*coeff[4]-643.9875775199394*coeff[1]+290.4737509655563*coeff[0])*ql[6]+1287.975155039879*coeff[1]*qc[6]+((-519.6152422706633*qr[5])+519.6152422706633*ql[5]-1195.115057222525*qr[4]+1195.115057222525*ql[4]+1207.476707849886*qr[1]+1207.476707849886*ql[1]+1609.968943799848*qc[1]-580.9475019311126*qr[0]+580.9475019311126*ql[0])*coeff[6]+((-402.4922359499622*coeff[3])-232.379000772445*coeff[2])*qr[5]+(232.379000772445*coeff[2]-402.4922359499622*coeff[3])*ql[5]+804.9844718999244*coeff[3]*qc[5]+(362.243012354966*qr[3]+362.243012354966*ql[3]+885.4829190899168*qc[3]-232.379000772445*qr[2]+232.379000772445*ql[2])*coeff[5]+((-643.9875775199395*coeff[3])-290.4737509655563*coeff[2])*qr[4]+(290.4737509655563*coeff[2]-643.9875775199395*coeff[3])*ql[4]+1287.975155039879*coeff[3]*qc[4]+(1207.476707849887*qr[3]+1207.476707849887*ql[3]+1609.968943799849*qc[3]-580.9475019311126*qr[2]+580.9475019311126*ql[2])*coeff[4]+(779.4228634059946*coeff[1]+405.0*coeff[0])*qr[3]+(405.0*coeff[0]-779.4228634059946*coeff[1])*ql[3]+990.0*coeff[0]*qc[3]+(779.4228634059946*qr[1]-779.4228634059946*ql[1]-450.0*qr[0]-450.0*ql[0]+900.0*qc[0])*coeff[3]+((-450.0*coeff[1])-259.8076211353315*coeff[0])*qr[2]+(259.8076211353315*coeff[0]-450.0*coeff[1])*ql[2]+900.0*coeff[1]*qc[2]+(405.0*qr[1]+405.0*ql[1]+990.0*qc[1]-259.8076211353315*qr[0]+259.8076211353315*ql[0])*coeff[2])*Jfac; 
  out[4] += -0.0015625*((4930.529890387037*coeff[7]+1510.463505020893*coeff[5])*qr[7]+(4930.529890387037*coeff[7]-1510.463505020893*coeff[5])*ql[7]+12678.50543242381*coeff[7]*qc[7]+(3253.30601081423*ql[5]-3253.30601081423*qr[5])*coeff[7]+((-11034.99546896146*coeff[6])-3091.710691510446*coeff[3]-210.0*coeff[2])*qr[6]+((-11034.99546896146*coeff[6])+3091.710691510446*coeff[3]-210.0*coeff[2])*ql[6]+(7982.762679674251*coeff[6]+13020.0*coeff[2])*qc[6]+(12340.86200392825*qr[3]-12340.86200392825*ql[3]-6600.000000000001*qr[2]-6600.000000000001*ql[2]+13200.0*qc[2])*coeff[6]-1341.640786499874*coeff[5]*qr[5]-1341.640786499874*coeff[5]*ql[5]+2683.281572999748*coeff[5]*qc[5]+((-11034.99546896146*coeff[4])-3091.710691510445*coeff[1]-210.0*coeff[0])*qr[4]+((-11034.99546896146*coeff[4])+3091.710691510445*coeff[1]-210.0*coeff[0])*ql[4]+(7982.762679674251*coeff[4]+13020.0*coeff[0])*qc[4]+(12340.86200392825*qr[1]-12340.86200392825*ql[1]-6600.0*qr[0]-6600.0*ql[0]+13200.0*qc[0])*coeff[4]+(4930.529890387037*coeff[3]+1510.463505020893*coeff[2])*qr[3]+(4930.529890387037*coeff[3]-1510.463505020893*coeff[2])*ql[3]+12678.50543242381*coeff[3]*qc[3]+(3253.30601081423*ql[2]-3253.30601081423*qr[2])*coeff[3]-1341.640786499874*coeff[2]*qr[2]-1341.640786499874*coeff[2]*ql[2]+2683.281572999748*coeff[2]*qc[2]+(4930.529890387037*coeff[1]+1510.463505020893*coeff[0])*qr[1]+(4930.529890387037*coeff[1]-1510.463505020893*coeff[0])*ql[1]+12678.50543242381*coeff[1]*qc[1]+(3253.30601081423*ql[0]-3253.30601081423*qr[0])*coeff[1]-1341.640786499874*coeff[0]*qr[0]-1341.640786499874*coeff[0]*ql[0]+2683.281572999748*coeff[0]*qc[0])*Jfac; 
  out[5] += -4.464285714285714e-4*((11068.53648862396*coeff[7]+5809.475019311126*coeff[5]+26433.11133786562*coeff[4]+17325.0*coeff[1]+9093.266739736608*coeff[0])*qr[7]+(11068.53648862396*coeff[7]-5809.475019311126*coeff[5]-26433.11133786562*coeff[4]+17325.0*coeff[1]-9093.266739736608*coeff[0])*ql[7]+(18112.1506177483*coeff[7]+28349.99999999999*coeff[1])*qc[7]+((-5809.475019311126*qr[5])+5809.475019311126*ql[5]-15859.86680271937*qr[4]+15859.86680271937*ql[4]+17325.0*qr[1]+17325.0*ql[1]+28349.99999999999*qc[1]-9093.266739736608*qr[0]+9093.266739736608*ql[0])*coeff[7]+((-24887.43658957267*coeff[6])-14185.49611398911*coeff[3]-6720.000000000001*coeff[2])*qr[6]+((-24887.43658957267*coeff[6])+14185.49611398911*coeff[3]-6720.000000000001*coeff[2])*ql[6]+(10330.63405604903*coeff[6]+13440.0*coeff[2])*qc[6]+(23642.49352331518*qr[3]-23642.49352331518*ql[3]-10500.0*qr[2]-10500.0*ql[2]+21000.0*qc[2])*coeff[6]+((-3354.101966249686*coeff[5])-11739.3568818739*coeff[4]-9093.266739736604*coeff[1]-5250.0*coeff[0])*qr[5]+((-3354.101966249686*coeff[5])-11739.3568818739*coeff[4]+9093.266739736604*coeff[1]-5250.0*coeff[0])*ql[5]+(6708.203932499371*coeff[5]+23478.7137637478*coeff[4]+10500.0*coeff[0])*qc[5]+((-7513.188404399295*qr[4])-7513.188404399295*ql[4]+15026.37680879859*qc[4]+9093.266739736604*qr[1]-9093.266739736604*ql[1]-5250.0*qr[0]-5250.0*ql[0]+10500.0*qc[0])*coeff[5]+(15495.95108407355*coeff[3]+8133.265027035576*coeff[2])*qr[3]+(15495.95108407355*coeff[3]-8133.265027035576*coeff[2])*ql[3]+25357.01086484762*coeff[3]*qc[3]+(8133.265027035576*ql[2]-8133.265027035576*qr[2])*coeff[3]-4695.742752749559*coeff[2]*qr[2]-4695.742752749559*coeff[2]*ql[2]+9391.485505499119*coeff[2]*qc[2])*Jfac; 
  out[6] += -0.0015625*((11038.00253669114*coeff[6]+4410.0*coeff[3]+1350.999629903724*coeff[2])*qr[7]+((-11038.00253669114*coeff[6])+4410.0*coeff[3]-1350.999629903724*coeff[2])*ql[7]+11340.0*coeff[3]*qc[7]+((-2765.310109192096*qr[6])+2765.310109192096*ql[6]+4410.0*qr[3]+4410.0*ql[3]+11340.0*qc[3]-2909.845356715713*qr[2]+2909.845356715713*ql[2])*coeff[7]+((-187.8297101099823*coeff[5])-11034.99546896146*coeff[4]-3091.710691510445*coeff[1]-210.0*coeff[0])*qr[6]+((-187.8297101099823*coeff[5])-11034.99546896146*coeff[4]+3091.710691510445*coeff[1]-210.0*coeff[0])*ql[6]+(11645.44202681891*coeff[5]+7982.762679674251*coeff[4]+13020.0*coeff[0])*qc[6]+((-5903.219460599446*qr[5])-5903.219460599446*ql[5]+11806.43892119889*qc[5]-11034.99546896146*qr[4]-11034.99546896146*ql[4]+7982.762679674251*qc[4]+12340.86200392825*qr[1]-12340.86200392825*ql[1]-6600.0*qr[0]-6600.0*ql[0]+13200.0*qc[0])*coeff[6]+((-2909.845356715714*coeff[3])-1200.0*coeff[2])*qr[5]+(2909.845356715714*coeff[3]-1200.0*coeff[2])*ql[5]+2400.0*coeff[2]*qc[5]+(1350.999629903725*qr[3]-1350.999629903725*ql[3]-1200.0*qr[2]-1200.0*ql[2]+2400.0*qc[2])*coeff[5]+((-3091.710691510446*coeff[3])-210.0*coeff[2])*qr[4]+(3091.710691510446*coeff[3]-210.0*coeff[2])*ql[4]+13020.0*coeff[2]*qc[4]+(12340.86200392825*qr[3]-12340.86200392825*ql[3]-6600.000000000001*qr[2]-6600.000000000001*ql[2]+13200.0*qc[2])*coeff[4]+(4930.529890387036*coeff[1]+1510.463505020893*coeff[0])*qr[3]+(4930.529890387036*coeff[1]-1510.463505020893*coeff[0])*ql[3]+12678.5054324238*coeff[1]*qc[3]+(4930.529890387036*qr[1]+4930.529890387036*ql[1]+12678.5054324238*qc[1]-3253.30601081423*qr[0]+3253.30601081423*ql[0])*coeff[3]+((-3253.30601081423*coeff[1])-1341.640786499874*coeff[0])*qr[2]+(3253.30601081423*coeff[1]-1341.640786499874*coeff[0])*ql[2]+2683.281572999748*coeff[0]*qc[2]+(1510.463505020893*qr[1]-1510.463505020893*ql[1]-1341.640786499874*qr[0]-1341.640786499874*ql[0]+2683.281572999748*qc[0])*coeff[2])*Jfac; 
  out[7] += -0.002232142857142857*((3485.685011586676*coeff[7]+1811.21506177483*coeff[5]+8452.336954949207*coeff[4]+5455.960043841962*coeff[1]+2835.0*coeff[0])*qr[7]+((-3485.685011586676*coeff[7])+1811.21506177483*coeff[5]+8452.336954949207*coeff[4]-5455.960043841962*coeff[1]+2835.0*coeff[0])*ql[7]+(4427.414595449584*coeff[5]+11269.78260659894*coeff[4]+6930.0*coeff[0])*qc[7]+((-2012.461179749811*qr[5])-2012.461179749811*ql[5]+4024.922359499622*qc[5]-4507.913042639576*qr[4]-4507.913042639576*ql[4]+9015.826085279152*qc[4]+5455.960043841962*qr[1]-5455.960043841962*ql[1]-3150.0*qr[0]-3150.0*ql[0]+6300.0*qc[0])*coeff[7]+((-7482.60382487273*coeff[6])-4032.0*coeff[3]-1818.653347947321*coeff[2])*qr[6]+(7482.60382487273*coeff[6]-4032.0*coeff[3]+1818.653347947321*coeff[2])*ql[6]+8064.0*coeff[3]*qc[6]+(7560.0*qr[3]+7560.0*ql[3]+10080.0*qc[3]-3637.306695894642*qr[2]+3637.306695894642*ql[2])*coeff[6]+((-1161.895003862225*coeff[5])-4066.632513517788*coeff[4]-3150.0*coeff[1]-1818.653347947322*coeff[0])*qr[5]+(1161.895003862225*coeff[5]+4066.632513517788*coeff[4]-3150.0*coeff[1]+1818.653347947322*coeff[0])*ql[5]+6300.000000000001*coeff[1]*qc[5]+((-2033.316256758894*qr[4])+2033.316256758894*ql[4]+2834.999999999999*qr[1]+2834.999999999999*ql[1]+6930.0*qc[1]-1818.653347947322*qr[0]+1818.653347947322*ql[0])*coeff[5]+(4879.959016221346*coeff[3]+2535.70108648476*coeff[2])*qr[3]+(2535.70108648476*coeff[2]-4879.959016221346*coeff[3])*ql[3]+6198.380433629417*coeff[2]*qc[3]+((-2817.445651649735*qr[2])-2817.445651649735*ql[2]+5634.89130329947*qc[2])*coeff[3]-1626.653005407115*coeff[2]*qr[2]+1626.653005407115*coeff[2]*ql[2])*Jfac; 

  return 0.;

}

