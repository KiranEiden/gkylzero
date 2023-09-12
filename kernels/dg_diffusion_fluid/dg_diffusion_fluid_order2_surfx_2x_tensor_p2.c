#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order2_surfx_2x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],2.);

  out[0] += 0.0125*(53.66563145999496*coeff[0]*qr[4]+53.66563145999496*coeff[0]*ql[4]-107.3312629199899*coeff[0]*qc[4]-95.26279441628824*coeff[0]*qr[1]+95.26279441628824*coeff[0]*ql[1]+75.0*coeff[0]*qr[0]+75.0*coeff[0]*ql[0]-150.0*coeff[0]*qc[0])*Jfac; 
  out[1] += 0.003125*(236.2519841186524*coeff[0]*qr[4]-236.2519841186524*coeff[0]*ql[4]-465.0*coeff[0]*qr[1]-465.0*coeff[0]*ql[1]-1710.0*coeff[0]*qc[1]+381.051177665153*coeff[0]*qr[0]-381.051177665153*coeff[0]*ql[0])*Jfac; 
  out[2] += 0.0125*(53.66563145999495*coeff[0]*qr[6]+53.66563145999495*coeff[0]*ql[6]-107.3312629199899*coeff[0]*qc[6]-95.26279441628824*coeff[0]*qr[3]+95.26279441628824*coeff[0]*ql[3]+75.0*coeff[0]*qr[2]+75.0*coeff[0]*ql[2]-150.0*coeff[0]*qc[2])*Jfac; 
  out[3] += 0.003125*(236.2519841186524*coeff[0]*qr[6]-236.2519841186524*coeff[0]*ql[6]-465.0*coeff[0]*qr[3]-465.0*coeff[0]*ql[3]-1710.0*coeff[0]*qc[3]+381.051177665153*coeff[0]*qr[2]-381.051177665153*coeff[0]*ql[2])*Jfac; 
  out[4] += -0.015625*(9.0*coeff[0]*qr[4]+9.0*coeff[0]*ql[4]+402.0*coeff[0]*qc[4]+19.36491673103709*coeff[0]*qr[1]-19.36491673103709*coeff[0]*ql[1]-26.83281572999748*coeff[0]*qr[0]-26.83281572999748*coeff[0]*ql[0]+53.66563145999496*coeff[0]*qc[0])*Jfac; 
  out[5] += 0.0125*(53.66563145999496*coeff[0]*qr[8]+53.66563145999496*coeff[0]*ql[8]-107.3312629199899*coeff[0]*qc[8]-95.26279441628826*coeff[0]*qr[7]+95.26279441628826*coeff[0]*ql[7]+75.0*coeff[0]*qr[5]+75.0*coeff[0]*ql[5]-150.0*coeff[0]*qc[5])*Jfac; 
  out[6] += -0.015625*(9.0*coeff[0]*qr[6]+9.0*coeff[0]*ql[6]+402.0*coeff[0]*qc[6]+19.36491673103708*coeff[0]*qr[3]-19.36491673103708*coeff[0]*ql[3]-26.83281572999747*coeff[0]*qr[2]-26.83281572999747*coeff[0]*ql[2]+53.66563145999495*coeff[0]*qc[2])*Jfac; 
  out[7] += 0.003125*(236.2519841186524*coeff[0]*qr[8]-236.2519841186524*coeff[0]*ql[8]-465.0*coeff[0]*qr[7]-465.0*coeff[0]*ql[7]-1710.0*coeff[0]*qc[7]+381.051177665153*coeff[0]*qr[5]-381.051177665153*coeff[0]*ql[5])*Jfac; 
  out[8] += -0.015625*(9.0*coeff[0]*qr[8]+9.0*coeff[0]*ql[8]+402.0*coeff[0]*qc[8]+19.36491673103708*coeff[0]*qr[7]-19.36491673103708*coeff[0]*ql[7]-26.83281572999748*coeff[0]*qr[5]-26.83281572999748*coeff[0]*ql[5]+53.66563145999496*coeff[0]*qc[5])*Jfac; 

  return 0.;

}

GKYL_CU_DH double dg_diffusion_fluid_order2_surfx_2x_tensor_p2_varcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],2.);

  out[0] += 0.00625*((120.0*coeff[8]+92.95160030897802*coeff[7]+53.66563145999496*coeff[5])*qr[8]+(120.0*coeff[8]-92.95160030897802*coeff[7]+53.66563145999496*coeff[5])*ql[8]+((-240.0*coeff[8])-107.3312629199899*coeff[5])*qc[8]+((-213.0140840414079*qr[7])+213.0140840414079*ql[7]+167.7050983124843*qr[5]+167.7050983124843*ql[5]-335.4101966249685*qc[5])*coeff[8]+((-165.0*coeff[7])-95.26279441628826*coeff[5])*qr[7]+(95.26279441628826*coeff[5]-165.0*coeff[7])*ql[7]-330.0*coeff[7]*qc[7]+(129.9038105676658*qr[5]-129.9038105676658*ql[5])*coeff[7]+(120.0*coeff[6]+92.95160030897802*coeff[3]+53.66563145999495*coeff[2])*qr[6]+(120.0*coeff[6]-92.95160030897802*coeff[3]+53.66563145999495*coeff[2])*ql[6]+((-240.0*coeff[6])-107.3312629199899*coeff[2])*qc[6]+((-213.0140840414079*qr[3])+213.0140840414079*ql[3]+167.7050983124842*qr[2]+167.7050983124842*ql[2]-335.4101966249685*qc[2])*coeff[6]+75.0*coeff[5]*qr[5]+75.0*coeff[5]*ql[5]-150.0*coeff[5]*qc[5]+(120.0*coeff[4]+92.951600308978*coeff[1]+53.66563145999496*coeff[0])*qr[4]+(120.0*coeff[4]-92.951600308978*coeff[1]+53.66563145999496*coeff[0])*ql[4]+((-240.0*coeff[4])-107.3312629199899*coeff[0])*qc[4]+((-213.014084041408*qr[1])+213.014084041408*ql[1]+167.7050983124843*qr[0]+167.7050983124843*ql[0]-335.4101966249685*qc[0])*coeff[4]+((-165.0*coeff[3])-95.26279441628824*coeff[2])*qr[3]+(95.26279441628824*coeff[2]-165.0*coeff[3])*ql[3]-330.0*coeff[3]*qc[3]+(129.9038105676658*qr[2]-129.9038105676658*ql[2])*coeff[3]+75.0*coeff[2]*qr[2]+75.0*coeff[2]*ql[2]-150.0*coeff[2]*qc[2]+((-165.0*coeff[1])-95.26279441628824*coeff[0])*qr[1]+(95.26279441628824*coeff[0]-165.0*coeff[1])*ql[1]-330.0*coeff[1]*qc[1]+(129.9038105676658*qr[0]-129.9038105676658*ql[0])*coeff[1]+75.0*coeff[0]*qr[0]+75.0*coeff[0]*ql[0]-150.0*coeff[0]*qc[0])*Jfac; 
  out[1] += 0.0015625*((528.2754963085075*coeff[8]+409.2004398824615*coeff[7]+236.2519841186524*coeff[5])*qr[8]+((-528.2754963085075*coeff[8])+409.2004398824615*coeff[7]-236.2519841186524*coeff[5])*ql[8]-1757.549430314835*coeff[7]*qc[8]+((-1039.771609537402*qr[7])-1039.771609537402*ql[7]-1677.050983124842*qc[7]+852.0563361656318*qr[5]-852.0563361656318*ql[5])*coeff[8]+((-805.4036255195278*coeff[7])-465.0*coeff[5])*qr[7]+(805.4036255195278*coeff[7]-465.0*coeff[5])*ql[7]-1710.0*coeff[5]*qc[7]+(660.0000000000001*qr[5]+660.0000000000001*ql[5]-1320.0*qc[5])*coeff[7]+(528.2754963085075*coeff[6]+409.2004398824615*coeff[3]+236.2519841186524*coeff[2])*qr[6]+((-528.2754963085075*coeff[6])+409.2004398824615*coeff[3]-236.2519841186524*coeff[2])*ql[6]-1757.549430314835*coeff[3]*qc[6]+((-1039.771609537402*qr[3])-1039.771609537402*ql[3]-1677.050983124842*qc[3]+852.0563361656317*qr[2]-852.0563361656317*ql[2])*coeff[6]+381.051177665153*coeff[5]*qr[5]-381.051177665153*coeff[5]*ql[5]+(528.2754963085075*coeff[4]+409.2004398824615*coeff[1]+236.2519841186524*coeff[0])*qr[4]+((-528.2754963085075*coeff[4])+409.2004398824615*coeff[1]-236.2519841186524*coeff[0])*ql[4]-1757.549430314835*coeff[1]*qc[4]+((-1039.771609537402*qr[1])-1039.771609537402*ql[1]-1677.050983124843*qc[1]+852.0563361656318*qr[0]-852.0563361656318*ql[0])*coeff[4]+((-805.4036255195278*coeff[3])-465.0*coeff[2])*qr[3]+(805.4036255195278*coeff[3]-465.0*coeff[2])*ql[3]-1710.0*coeff[2]*qc[3]+(660.0*qr[2]+660.0*ql[2]-1320.0*qc[2])*coeff[3]+381.051177665153*coeff[2]*qr[2]-381.051177665153*coeff[2]*ql[2]+((-805.4036255195278*coeff[1])-465.0*coeff[0])*qr[1]+(805.4036255195278*coeff[1]-465.0*coeff[0])*ql[1]-1710.0*coeff[0]*qc[1]+(660.0*qr[0]+660.0*ql[0]-1320.0*qc[0])*coeff[1]+381.051177665153*coeff[0]*qr[0]-381.051177665153*coeff[0]*ql[0])*Jfac; 
  out[2] += 0.00125*((536.6563145999495*coeff[6]+415.6921938165305*coeff[3]+240.0*coeff[2])*qr[8]+(536.6563145999495*coeff[6]-415.6921938165305*coeff[3]+240.0*coeff[2])*ql[8]+((-1073.312629199899*coeff[6])-480.0*coeff[2])*qc[8]+(536.6563145999495*qr[6]+536.6563145999495*ql[6]-1073.312629199899*qc[6]-952.6279441628825*qr[3]+952.6279441628825*ql[3]+750.0*qr[2]+750.0*ql[2]-1500.0*qc[2])*coeff[8]+((-952.6279441628825*coeff[6])-737.9024325749306*coeff[3]-426.0281680828159*coeff[2])*qr[7]+(952.6279441628825*coeff[6]-737.9024325749306*coeff[3]+426.0281680828159*coeff[2])*ql[7]-1475.804865149861*coeff[3]*qc[7]+(415.6921938165305*qr[6]-415.6921938165305*ql[6]-737.9024325749306*qr[3]-737.9024325749306*ql[3]-1475.804865149861*qc[3]+580.9475019311126*qr[2]-580.9475019311126*ql[2])*coeff[7]+(240.0*coeff[5]+600.0000000000001*coeff[4]+464.7580015448901*coeff[1]+268.3281572999747*coeff[0])*qr[6]+(240.0*coeff[5]+600.0000000000001*coeff[4]-464.7580015448901*coeff[1]+268.3281572999747*coeff[0])*ql[6]+((-480.0000000000001*coeff[5])-1200.0*coeff[4]-536.6563145999495*coeff[0])*qc[6]+(750.0000000000001*qr[5]+750.0000000000001*ql[5]-1500.0*qc[5]+600.0000000000001*qr[4]+600.0000000000001*ql[4]-1200.0*qc[4]-1065.07042020704*qr[1]+1065.07042020704*ql[1]+838.5254915624212*qr[0]+838.5254915624212*ql[0]-1677.050983124842*qc[0])*coeff[6]+(580.9475019311126*coeff[3]+335.4101966249685*coeff[2])*qr[5]+(335.4101966249685*coeff[2]-580.9475019311126*coeff[3])*ql[5]-670.8203932499371*coeff[2]*qc[5]+((-426.0281680828159*qr[3])+426.0281680828159*ql[3]+335.4101966249685*qr[2]+335.4101966249685*ql[2]-670.8203932499371*qc[2])*coeff[5]+(464.7580015448901*coeff[3]+268.3281572999748*coeff[2])*qr[4]+(268.3281572999748*coeff[2]-464.7580015448901*coeff[3])*ql[4]-536.6563145999496*coeff[2]*qc[4]+((-1065.07042020704*qr[3])+1065.07042020704*ql[3]+838.5254915624214*qr[2]+838.5254915624214*ql[2]-1677.050983124843*qc[2])*coeff[4]+((-825.0*coeff[1])-476.3139720814412*coeff[0])*qr[3]+(476.3139720814412*coeff[0]-825.0*coeff[1])*ql[3]-1650.0*coeff[1]*qc[3]+((-825.0*qr[1])-825.0*ql[1]-1650.0*qc[1]+649.5190528383289*qr[0]-649.5190528383289*ql[0])*coeff[3]+(649.5190528383289*coeff[1]+375.0*coeff[0])*qr[2]+(375.0*coeff[0]-649.5190528383289*coeff[1])*ql[2]-750.0*coeff[0]*qc[2]+((-476.3139720814412*qr[1])+476.3139720814412*ql[1]+375.0*qr[0]+375.0*ql[0]-750.0*qc[0])*coeff[2])*Jfac; 
  out[3] += 3.125e-4*((2362.519841186524*coeff[6]+1830.0*coeff[3]+1056.550992617015*coeff[2])*qr[8]+((-2362.519841186524*coeff[6])+1830.0*coeff[3]-1056.550992617015*coeff[2])*ql[8]-7860.0*coeff[3]*qc[8]+(2362.519841186524*qr[6]-2362.519841186524*ql[6]-4650.0*qr[3]-4650.0*ql[3]-7500.0*qc[3]+3810.51177665153*qr[2]-3810.51177665153*ql[2])*coeff[8]+((-4650.0*coeff[6])-3601.874511972898*coeff[3]-2079.543219074804*coeff[2])*qr[7]+((-4650.0*coeff[6])+3601.874511972898*coeff[3]-2079.543219074804*coeff[2])*ql[7]+((-7500.0*coeff[6])-7647.352483049281*coeff[2])*qc[7]+(1830.0*qr[6]+1830.0*ql[6]-7860.0*qc[6]-3601.874511972898*qr[3]+3601.874511972898*ql[3]+2951.609730299722*qr[2]+2951.609730299722*ql[2]-5903.219460599445*qc[2])*coeff[7]+(1056.550992617015*coeff[5]+2641.377481542539*coeff[4]+2046.002199412307*coeff[1]+1181.259920593262*coeff[0])*qr[6]+((-1056.550992617015*coeff[5])-2641.377481542539*coeff[4]+2046.002199412307*coeff[1]-1181.259920593262*coeff[0])*ql[6]-8787.747151574174*coeff[1]*qc[6]+(3810.511776651531*qr[5]-3810.511776651531*ql[5]+2641.377481542539*qr[4]-2641.377481542539*ql[4]-5198.858047687011*qr[1]-5198.858047687011*ql[1]-8385.254915624211*qc[1]+4260.281680828159*qr[0]-4260.281680828159*ql[0])*coeff[6]+(2951.609730299723*coeff[3]+1704.112672331264*coeff[2])*qr[5]+(2951.609730299723*coeff[3]-1704.112672331264*coeff[2])*ql[5]-5903.219460599446*coeff[3]*qc[5]+((-2079.543219074805*qr[3])-2079.543219074805*ql[3]-7647.352483049282*qc[3]+1704.112672331264*qr[2]-1704.112672331264*ql[2])*coeff[5]+(2046.002199412308*coeff[3]+1181.259920593262*coeff[2])*qr[4]+(2046.002199412308*coeff[3]-1181.259920593262*coeff[2])*ql[4]-8787.747151574174*coeff[3]*qc[4]+((-5198.858047687012*qr[3])-5198.858047687012*ql[3]-8385.254915624215*qc[3]+4260.28168082816*qr[2]-4260.28168082816*ql[2])*coeff[4]+((-4027.018127597639*coeff[1])-2325.0*coeff[0])*qr[3]+(4027.018127597639*coeff[1]-2325.0*coeff[0])*ql[3]-8550.0*coeff[0]*qc[3]+((-4027.018127597639*qr[1])+4027.018127597639*ql[1]+3300.0*qr[0]+3300.0*ql[0]-6600.0*qc[0])*coeff[3]+(3300.0*coeff[1]+1905.255888325765*coeff[0])*qr[2]+(3300.0*coeff[1]-1905.255888325765*coeff[0])*ql[2]-6600.0*coeff[1]*qc[2]+((-2325.0*qr[1])-2325.0*ql[1]-8550.0*qc[1]+1905.255888325765*qr[0]-1905.255888325765*ql[0])*coeff[2])*Jfac; 
  out[4] += -0.0015625*((100.6230589874906*coeff[8]+77.94228634059948*coeff[7]+45.0*coeff[5])*qr[8]+(100.6230589874906*coeff[8]-77.94228634059948*coeff[7]+45.0*coeff[5])*ql[8]+(2010.0*coeff[5]-1945.379140424817*coeff[8])*qc[8]+(216.5063509461097*qr[7]-216.5063509461097*ql[7]-300.0*qr[5]-300.0*ql[5]+600.0*qc[5])*coeff[8]+(167.7050983124843*coeff[7]+96.82458365518542*coeff[5])*qr[7]+(167.7050983124843*coeff[7]-96.82458365518542*coeff[5])*ql[7]+1274.55874717488*coeff[7]*qc[7]+(232.379000772445*ql[5]-232.379000772445*qr[5])*coeff[7]+(100.6230589874906*coeff[6]+77.94228634059948*coeff[3]+45.0*coeff[2])*qr[6]+(100.6230589874906*coeff[6]-77.94228634059948*coeff[3]+45.0*coeff[2])*ql[6]+(2010.0*coeff[2]-1945.379140424817*coeff[6])*qc[6]+(216.5063509461097*qr[3]-216.5063509461097*ql[3]-300.0000000000001*qr[2]-300.0000000000001*ql[2]+600.0000000000001*qc[2])*coeff[6]-134.1640786499874*coeff[5]*qr[5]-134.1640786499874*coeff[5]*ql[5]+268.3281572999748*coeff[5]*qc[5]+(100.6230589874906*coeff[4]+77.94228634059945*coeff[1]+45.0*coeff[0])*qr[4]+(100.6230589874906*coeff[4]-77.94228634059945*coeff[1]+45.0*coeff[0])*ql[4]+(2010.0*coeff[0]-1945.379140424817*coeff[4])*qc[4]+(216.5063509461096*qr[1]-216.5063509461096*ql[1]-300.0*qr[0]-300.0*ql[0]+600.0*qc[0])*coeff[4]+(167.7050983124843*coeff[3]+96.82458365518544*coeff[2])*qr[3]+(167.7050983124843*coeff[3]-96.82458365518544*coeff[2])*ql[3]+1274.55874717488*coeff[3]*qc[3]+(232.379000772445*ql[2]-232.379000772445*qr[2])*coeff[3]-134.1640786499874*coeff[2]*qr[2]-134.1640786499874*coeff[2]*ql[2]+268.3281572999748*coeff[2]*qc[2]+(167.7050983124843*coeff[1]+96.82458365518544*coeff[0])*qr[1]+(167.7050983124843*coeff[1]-96.82458365518544*coeff[0])*ql[1]+1274.55874717488*coeff[1]*qc[1]+(232.379000772445*ql[0]-232.379000772445*qr[0])*coeff[1]-134.1640786499874*coeff[0]*qr[0]-134.1640786499874*coeff[0]*ql[0]+268.3281572999748*coeff[0]*qc[0])*Jfac; 
  out[5] += 1.785714285714285e-4*((2683.281572999748*coeff[8]+2078.460969082653*coeff[7]+1200.0*coeff[5]+4200.0*coeff[4]+3253.30601081423*coeff[1]+1878.297101099824*coeff[0])*qr[8]+(2683.281572999748*coeff[8]-2078.460969082653*coeff[7]+1200.0*coeff[5]+4200.0*coeff[4]-3253.30601081423*coeff[1]+1878.297101099824*coeff[0])*ql[8]+((-5366.563145999497*coeff[8])-2400.0*coeff[5]-8400.0*coeff[4]-3756.594202199647*coeff[0])*qc[8]+((-4763.139720814414*qr[7])+4763.139720814414*ql[7]+3750.0*qr[5]+3750.0*ql[5]-7500.0*qc[5]+4200.0*qr[4]+4200.0*ql[4]-8400.0*qc[4]-7455.492941449278*qr[1]+7455.492941449278*ql[1]+5869.678440936949*qr[0]+5869.678440936949*ql[0]-11739.3568818739*qc[0])*coeff[8]+((-3689.512162874654*coeff[7])-2130.140840414079*coeff[5]-7455.492941449278*coeff[4]-5775.000000000001*coeff[1]-3334.19780457009*coeff[0])*qr[7]+((-3689.512162874654*coeff[7])+2130.140840414079*coeff[5]+7455.492941449278*coeff[4]-5775.000000000001*coeff[1]+3334.19780457009*coeff[0])*ql[7]+((-7379.024325749308*coeff[7])-11550.0*coeff[1])*qc[7]+(2904.737509655563*qr[5]-2904.737509655563*ql[5]+3253.30601081423*qr[4]-3253.30601081423*ql[4]-5775.000000000001*qr[1]-5775.000000000001*ql[1]-11550.0*qc[1]+4546.633369868304*qr[0]-4546.633369868304*ql[0])*coeff[7]+(3756.594202199647*coeff[6]+2909.845356715714*coeff[3]+1680.0*coeff[2])*qr[6]+(3756.594202199647*coeff[6]-2909.845356715714*coeff[3]+1680.0*coeff[2])*ql[6]+((-7513.188404399295*coeff[6])-3360.0*coeff[2])*qc[6]+((-6668.39560914018*qr[3])+6668.39560914018*ql[3]+5250.000000000001*qr[2]+5250.000000000001*ql[2]-10500.0*qc[2])*coeff[6]+(1677.050983124843*coeff[5]+5869.678440936949*coeff[4]+4546.633369868302*coeff[1]+2625.0*coeff[0])*qr[5]+(1677.050983124843*coeff[5]+5869.678440936949*coeff[4]-4546.633369868302*coeff[1]+2625.0*coeff[0])*ql[5]+((-3354.101966249686*coeff[5])-11739.3568818739*coeff[4]-5250.0*coeff[0])*qc[5]+(1878.297101099824*qr[4]+1878.297101099824*ql[4]-3756.594202199647*qc[4]-3334.197804570088*qr[1]+3334.197804570088*ql[1]+2625.0*qr[0]+2625.0*ql[0]-5250.0*qc[0])*coeff[5]+((-5165.317028024515*coeff[3])-2982.197176579711*coeff[2])*qr[3]+(2982.197176579711*coeff[2]-5165.317028024515*coeff[3])*ql[3]-10330.63405604903*coeff[3]*qc[3]+(4066.632513517788*qr[2]-4066.632513517788*ql[2])*coeff[3]+2347.87137637478*coeff[2]*qr[2]+2347.87137637478*coeff[2]*ql[2]-4695.742752749559*coeff[2]*qc[2])*Jfac; 
  out[6] += -0.0015625*((90.0*coeff[6]+69.71370023173351*coeff[3]+40.24922359499621*coeff[2])*qr[8]+(90.0*coeff[6]-69.71370023173351*coeff[3]+40.24922359499621*coeff[2])*ql[8]+(1797.798653909831*coeff[2]-1740.0*coeff[6])*qc[8]+(90.0*qr[6]+90.0*ql[6]-1740.0*qc[6]+193.6491673103708*qr[3]-193.6491673103708*ql[3]-268.3281572999747*qr[2]-268.3281572999747*ql[2]+536.6563145999495*qc[2])*coeff[8]+(193.6491673103708*coeff[6]+150.0*coeff[3]+86.60254037844386*coeff[2])*qr[7]+((-193.6491673103708*coeff[6])+150.0*coeff[3]-86.60254037844386*coeff[2])*ql[7]+1140.0*coeff[3]*qc[7]+(69.71370023173351*qr[6]-69.71370023173351*ql[6]+150.0*qr[3]+150.0*ql[3]+1140.0*qc[3]-207.8460969082653*qr[2]+207.8460969082653*ql[2])*coeff[7]+(40.24922359499622*coeff[5]+100.6230589874906*coeff[4]+77.94228634059945*coeff[1]+45.0*coeff[0])*qr[6]+(40.24922359499622*coeff[5]+100.6230589874906*coeff[4]-77.94228634059945*coeff[1]+45.0*coeff[0])*ql[6]+(1797.798653909831*coeff[5]-1945.379140424817*coeff[4]+2010.0*coeff[0])*qc[6]+((-268.3281572999748*qr[5])-268.3281572999748*ql[5]+536.6563145999496*qc[5]+100.6230589874906*qr[4]+100.6230589874906*ql[4]-1945.379140424817*qc[4]+216.5063509461096*qr[1]-216.5063509461096*ql[1]-300.0*qr[0]-300.0*ql[0]+600.0*qc[0])*coeff[6]+((-207.8460969082653*coeff[3])-120.0*coeff[2])*qr[5]+(207.8460969082653*coeff[3]-120.0*coeff[2])*ql[5]+240.0*coeff[2]*qc[5]+(86.60254037844389*qr[3]-86.60254037844389*ql[3]-120.0*qr[2]-120.0*ql[2]+240.0*qc[2])*coeff[5]+(77.94228634059948*coeff[3]+45.0*coeff[2])*qr[4]+(45.0*coeff[2]-77.94228634059948*coeff[3])*ql[4]+2010.0*coeff[2]*qc[4]+(216.5063509461097*qr[3]-216.5063509461097*ql[3]-300.0000000000001*qr[2]-300.0000000000001*ql[2]+600.0000000000001*qc[2])*coeff[4]+(167.7050983124842*coeff[1]+96.82458365518542*coeff[0])*qr[3]+(167.7050983124842*coeff[1]-96.82458365518542*coeff[0])*ql[3]+1274.55874717488*coeff[1]*qc[3]+(167.7050983124842*qr[1]+167.7050983124842*ql[1]+1274.55874717488*qc[1]-232.379000772445*qr[0]+232.379000772445*ql[0])*coeff[3]+((-232.379000772445*coeff[1])-134.1640786499874*coeff[0])*qr[2]+(232.379000772445*coeff[1]-134.1640786499874*coeff[0])*ql[2]+268.3281572999747*coeff[0]*qc[2]+(96.82458365518542*qr[1]-96.82458365518542*ql[1]-134.1640786499874*qr[0]-134.1640786499874*ql[0]+268.3281572999747*qc[0])*coeff[2])*Jfac; 
  out[7] += 2.232142857142857e-4*((2362.519841186524*coeff[8]+1830.0*coeff[7]+1056.550992617015*coeff[5]+3697.928474159553*coeff[4]+2864.403079177231*coeff[1]+1653.763888830567*coeff[0])*qr[8]+((-2362.519841186524*coeff[8])+1830.0*coeff[7]-1056.550992617015*coeff[5]-3697.928474159553*coeff[4]+2864.403079177231*coeff[1]-1653.763888830567*coeff[0])*ql[8]+((-7860.0*coeff[7])-12302.84601220384*coeff[1])*qc[8]+((-4650.0*qr[7])-4650.0*ql[7]-7500.0*qc[7]+3810.511776651531*qr[5]-3810.511776651531*ql[5]+3697.928474159553*qr[4]-3697.928474159553*ql[4]-7278.401266761815*qr[1]-7278.401266761815*ql[1]-11739.35688187389*qc[1]+5964.394353159422*qr[0]-5964.394353159422*ql[0])*coeff[8]+((-3601.874511972898*coeff[7])-2079.543219074805*coeff[5]-7278.401266761817*coeff[4]-5637.825378636695*coeff[1]-3255.0*coeff[0])*qr[7]+(3601.874511972898*coeff[7]-2079.543219074805*coeff[5]-7278.401266761817*coeff[4]+5637.825378636695*coeff[1]-3255.0*coeff[0])*ql[7]+((-7647.352483049282*coeff[5])-11739.3568818739*coeff[4]-11970.0*coeff[0])*qc[7]+(2951.609730299723*qr[5]+2951.609730299723*ql[5]-5903.219460599446*qc[5]+2864.403079177231*qr[4]+2864.403079177231*ql[4]-12302.84601220384*qc[4]-5637.825378636695*qr[1]+5637.825378636695*ql[1]+4620.0*qr[0]+4620.0*ql[0]-9240.0*qc[0])*coeff[7]+(3307.527777661134*coeff[6]+2562.0*coeff[3]+1479.171389663821*coeff[2])*qr[6]+((-3307.527777661134*coeff[6])+2562.0*coeff[3]-1479.171389663821*coeff[2])*ql[6]-11004.0*coeff[3]*qc[6]+((-6510.0*qr[3])-6510.0*ql[3]-10500.0*qc[3]+5334.716487312142*qr[2]-5334.716487312142*ql[2])*coeff[6]+(1704.112672331263*coeff[5]+5964.394353159422*coeff[4]+4620.0*coeff[1]+2667.358243656071*coeff[0])*qr[5]+((-1704.112672331263*coeff[5])-5964.394353159422*coeff[4]+4620.0*coeff[1]-2667.358243656071*coeff[0])*ql[5]-9240.0*coeff[1]*qc[5]+(1653.763888830567*qr[4]-1653.763888830567*ql[4]-3255.0*qr[1]-3255.0*ql[1]-11970.0*qc[1]+2667.358243656071*qr[0]-2667.358243656071*ql[0])*coeff[5]+((-5042.624316762057*coeff[3])-2911.360506704726*coeff[2])*qr[3]+(5042.624316762057*coeff[3]-2911.360506704726*coeff[2])*ql[3]-10706.29347626899*coeff[2]*qc[3]+(4132.253622419611*qr[2]+4132.253622419611*ql[2]-8264.507244839222*qc[2])*coeff[3]+2385.757741263769*coeff[2]*qr[2]-2385.757741263769*coeff[2]*ql[2])*Jfac; 
  out[8] += -2.232142857142857e-4*((450.0*coeff[8]+348.5685011586676*coeff[7]+201.2461179749811*coeff[5]+704.3614129124339*coeff[4]+545.5960043841961*coeff[1]+315.0*coeff[0])*qr[8]+(450.0*coeff[8]-348.5685011586676*coeff[7]+201.2461179749811*coeff[5]+704.3614129124339*coeff[4]-545.5960043841961*coeff[1]+315.0*coeff[0])*ql[8]+((-8700.0*coeff[8])+8988.993269549157*coeff[5]-13617.65398297372*coeff[4]+14070.0*coeff[0])*qc[8]+(968.2458365518543*qr[7]-968.2458365518543*ql[7]-1341.640786499874*qr[5]-1341.640786499874*ql[5]+2683.281572999748*qc[5]+704.3614129124339*qr[4]+704.3614129124339*ql[4]-13617.65398297372*qc[4]+1515.544456622767*qr[1]-1515.544456622767*ql[1]-2100.0*qr[0]-2100.0*ql[0]+4200.0*qc[0])*coeff[8]+(750.0*coeff[7]+433.0127018922195*coeff[5]+1515.544456622768*coeff[4]+1173.93568818739*coeff[1]+677.772085586298*coeff[0])*qr[7]+(750.0*coeff[7]-433.0127018922195*coeff[5]-1515.544456622768*coeff[4]+1173.93568818739*coeff[1]-677.772085586298*coeff[0])*ql[7]+(5700.0*coeff[7]+8921.91123022416*coeff[1])*qc[7]+((-1039.230484541327*qr[5])+1039.230484541327*ql[5]+545.5960043841964*qr[4]-545.5960043841964*ql[4]+1173.93568818739*qr[1]+1173.93568818739*ql[1]+8921.91123022416*qc[1]-1626.653005407115*qr[0]+1626.653005407115*ql[0])*coeff[7]+(630.0*coeff[6]+487.9959016221346*coeff[3]+281.7445651649734*coeff[2])*qr[6]+(630.0*coeff[6]-487.9959016221346*coeff[3]+281.7445651649734*coeff[2])*ql[6]+(12584.59057736882*coeff[2]-12180.0*coeff[6])*qc[6]+(1355.544171172596*qr[3]-1355.544171172596*ql[3]-1878.297101099823*qr[2]-1878.297101099823*ql[2]+3756.594202199647*qc[2])*coeff[6]+((-600.0*coeff[5])-2100.0*coeff[4]-1626.653005407115*coeff[1]-939.1485505499119*coeff[0])*qr[5]+((-600.0*coeff[5])-2100.0*coeff[4]+1626.653005407115*coeff[1]-939.1485505499119*coeff[0])*ql[5]+(1200.0*coeff[5]+4200.0*coeff[4]+1878.297101099824*coeff[0])*qc[5]+(315.0*qr[4]+315.0*ql[4]+14070.0*qc[4]+677.7720855862981*qr[1]-677.7720855862981*ql[1]-939.1485505499119*qr[0]-939.1485505499119*ql[0]+1878.297101099824*qc[0])*coeff[5]+(1050.0*coeff[3]+606.217782649107*coeff[2])*qr[3]+(1050.0*coeff[3]-606.217782649107*coeff[2])*ql[3]+7980.0*coeff[3]*qc[3]+(1454.922678357857*ql[2]-1454.922678357857*qr[2])*coeff[3]-840.0*coeff[2]*qr[2]-840.0*coeff[2]*ql[2]+1680.0*coeff[2]*qc[2])*Jfac; 

  return 0.;

}
