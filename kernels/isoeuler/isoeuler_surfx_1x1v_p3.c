#include <gkyl_isoeuler_kernels.h> 
GKYL_CU_DH void isoeuler_surfx_1x1v_ser_p3(const double *w, const double *dxv, double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out)  
{ 
  // w: Cell-center coordinates. dxv[NDIM]: Cell spacing. statevec(l/c/r): [rho, rho ux, rho uy, rho uz] in (left/center/right) cell, uvar(l/c/r): [ux, uv, uz]  in (left/center/right) cell
  // out: output in left/right cells 
  const double dxl1 = 2.0/dxv[0]; 
  const double dxr1 = 2.0/dxv[0]; 
  const double *rhol = &statevecl[0]; 
  const double *rhou0l = &statevecl[4]; 
  const double *rhor = &statevecr[0]; 
  const double *rhou0r = &statevecr[4]; 
  const double *uvar0l = &uvarl[0]; 
  const double *uvar0r = &uvarr[0]; 
  double *outlrho = &out[0]; 
  double *outlrhoux = &out[4]; 
  double *outrrho = &out[8]; 
  double *outrrhoux = &out[12]; 
  double incr[4]; 

  double vthsq = vth*vth; 
  //FluxRho; 
  incr[0] = (-0.6614378277661477*rhou0r[3])+0.6614378277661477*rhou0l[3]+0.5590169943749475*rhou0r[2]+0.5590169943749475*rhou0l[2]-0.4330127018922193*rhou0r[1]+0.4330127018922193*rhou0l[1]+0.25*rhou0r[0]+0.25*rhou0l[0]; 
  incr[1] = 1.14564392373896*rhou0r[3]-1.14564392373896*rhou0l[3]-0.9682458365518543*rhou0r[2]-0.9682458365518543*rhou0l[2]+0.75*rhou0r[1]-0.75*rhou0l[1]-0.4330127018922193*rhou0r[0]-0.4330127018922193*rhou0l[0]; 
  incr[2] = (-1.479019945774904*rhou0r[3])+1.479019945774904*rhou0l[3]+1.25*rhou0r[2]+1.25*rhou0l[2]-0.9682458365518543*rhou0r[1]+0.9682458365518543*rhou0l[1]+0.5590169943749475*rhou0r[0]+0.5590169943749475*rhou0l[0]; 
  incr[3] = 1.75*rhou0r[3]-1.75*rhou0l[3]-1.479019945774904*rhou0r[2]-1.479019945774904*rhou0l[2]+1.14564392373896*rhou0r[1]-1.14564392373896*rhou0l[1]-0.6614378277661477*rhou0r[0]-0.6614378277661477*rhou0l[0]; 

  outrrho[0] += incr[0]*dxr1; 
  outrrho[1] += incr[1]*dxr1; 
  outrrho[2] += incr[2]*dxr1; 
  outrrho[3] += incr[3]*dxr1; 

  outlrho[0] += -1.0*incr[0]*dxl1; 
  outlrho[1] += incr[1]*dxl1; 
  outlrho[2] += -1.0*incr[2]*dxl1; 
  outlrho[3] += incr[3]*dxl1; 

 
  //FluxRhoUx; 
  incr[0] = 2.806243040080456*rhor[3]*vthsq+2.806243040080456*rhol[3]*vthsq-1.185854122563142*rhor[2]*vthsq+1.185854122563142*rhol[2]*vthsq+0.3061862178478972*rhor[1]*vthsq+0.3061862178478972*rhol[1]*vthsq+1.237436867076458*rhor[3]*uvar0r[3]-1.045825033167594*rhor[2]*uvar0r[3]+0.8100925873009825*rhor[1]*uvar0r[3]-0.4677071733467427*rhor[0]*uvar0r[3]+1.237436867076458*rhol[3]*uvar0l[3]+1.045825033167594*rhol[2]*uvar0l[3]+0.8100925873009825*rhol[1]*uvar0l[3]+0.4677071733467427*rhol[0]*uvar0l[3]-1.045825033167594*uvar0r[2]*rhor[3]+0.8100925873009825*uvar0r[1]*rhor[3]-0.4677071733467427*uvar0r[0]*rhor[3]+1.045825033167594*uvar0l[2]*rhol[3]+0.8100925873009825*uvar0l[1]*rhol[3]+0.4677071733467427*uvar0l[0]*rhol[3]+0.8838834764831843*rhor[2]*uvar0r[2]-0.6846531968814576*rhor[1]*uvar0r[2]+0.3952847075210474*rhor[0]*uvar0r[2]+0.8838834764831843*rhol[2]*uvar0l[2]+0.6846531968814576*rhol[1]*uvar0l[2]+0.3952847075210474*rhol[0]*uvar0l[2]-0.6846531968814576*uvar0r[1]*rhor[2]+0.3952847075210474*uvar0r[0]*rhor[2]+0.6846531968814576*uvar0l[1]*rhol[2]+0.3952847075210474*uvar0l[0]*rhol[2]+0.5303300858899106*rhor[1]*uvar0r[1]-0.3061862178478972*rhor[0]*uvar0r[1]+0.5303300858899106*rhol[1]*uvar0l[1]+0.3061862178478972*rhol[0]*uvar0l[1]-0.3061862178478972*uvar0r[0]*rhor[1]+0.3061862178478972*uvar0l[0]*rhol[1]+0.1767766952966369*rhor[0]*uvar0r[0]+0.1767766952966369*rhol[0]*uvar0l[0]; 
  incr[1] = 8.418729120241368*rhor[3]*vthsq-8.418729120241368*rhol[3]*vthsq-3.557562367689427*rhor[2]*vthsq-3.557562367689427*rhol[2]*vthsq+0.9185586535436916*rhor[1]*vthsq-0.9185586535436916*rhol[1]*vthsq-2.143303524935281*rhor[3]*uvar0r[3]+1.81142209327368*rhor[2]*uvar0r[3]-1.403121520040228*rhor[1]*uvar0r[3]+0.8100925873009825*rhor[0]*uvar0r[3]-2.143303524935281*rhol[3]*uvar0l[3]-1.81142209327368*rhol[2]*uvar0l[3]-1.403121520040228*rhol[1]*uvar0l[3]-0.8100925873009825*rhol[0]*uvar0l[3]+1.81142209327368*uvar0r[2]*rhor[3]-1.403121520040228*uvar0r[1]*rhor[3]+0.8100925873009825*uvar0r[0]*rhor[3]-1.81142209327368*uvar0l[2]*rhol[3]-1.403121520040228*uvar0l[1]*rhol[3]-0.8100925873009825*uvar0l[0]*rhol[3]-1.530931089239486*rhor[2]*uvar0r[2]+1.185854122563142*rhor[1]*uvar0r[2]-0.6846531968814576*rhor[0]*uvar0r[2]-1.530931089239486*rhol[2]*uvar0l[2]-1.185854122563142*rhol[1]*uvar0l[2]-0.6846531968814576*rhol[0]*uvar0l[2]+1.185854122563142*uvar0r[1]*rhor[2]-0.6846531968814576*uvar0r[0]*rhor[2]-1.185854122563142*uvar0l[1]*rhol[2]-0.6846531968814576*uvar0l[0]*rhol[2]-0.9185586535436916*rhor[1]*uvar0r[1]+0.5303300858899106*rhor[0]*uvar0r[1]-0.9185586535436916*rhol[1]*uvar0l[1]-0.5303300858899106*rhol[0]*uvar0l[1]+0.5303300858899106*uvar0r[0]*rhor[1]-0.5303300858899106*uvar0l[0]*rhol[1]-0.3061862178478972*rhor[0]*uvar0r[0]-0.3061862178478972*rhol[0]*uvar0l[0]; 
  incr[2] = 14.03121520040228*rhor[3]*vthsq+14.03121520040228*rhol[3]*vthsq-5.929270612815712*rhor[2]*vthsq+5.929270612815712*rhol[2]*vthsq+1.530931089239486*rhor[1]*vthsq+1.530931089239486*rhol[1]*vthsq+2.766992952647331*rhor[3]*uvar0r[3]-2.338535866733713*rhor[2]*uvar0r[3]+1.81142209327368*rhor[1]*uvar0r[3]-1.045825033167594*rhor[0]*uvar0r[3]+2.766992952647331*rhol[3]*uvar0l[3]+2.338535866733713*rhol[2]*uvar0l[3]+1.81142209327368*rhol[1]*uvar0l[3]+1.045825033167594*rhol[0]*uvar0l[3]-2.338535866733713*uvar0r[2]*rhor[3]+1.81142209327368*uvar0r[1]*rhor[3]-1.045825033167594*uvar0r[0]*rhor[3]+2.338535866733713*uvar0l[2]*rhol[3]+1.81142209327368*uvar0l[1]*rhol[3]+1.045825033167594*uvar0l[0]*rhol[3]+1.976423537605237*rhor[2]*uvar0r[2]-1.530931089239486*rhor[1]*uvar0r[2]+0.8838834764831843*rhor[0]*uvar0r[2]+1.976423537605237*rhol[2]*uvar0l[2]+1.530931089239486*rhol[1]*uvar0l[2]+0.8838834764831843*rhol[0]*uvar0l[2]-1.530931089239486*uvar0r[1]*rhor[2]+0.8838834764831843*uvar0r[0]*rhor[2]+1.530931089239486*uvar0l[1]*rhol[2]+0.8838834764831843*uvar0l[0]*rhol[2]+1.185854122563142*rhor[1]*uvar0r[1]-0.6846531968814576*rhor[0]*uvar0r[1]+1.185854122563142*rhol[1]*uvar0l[1]+0.6846531968814576*rhol[0]*uvar0l[1]-0.6846531968814576*uvar0r[0]*rhor[1]+0.6846531968814576*uvar0l[0]*rhol[1]+0.3952847075210474*rhor[0]*uvar0r[0]+0.3952847075210474*rhol[0]*uvar0l[0]; 
  incr[3] = 19.64370128056319*rhor[3]*vthsq-19.64370128056319*rhol[3]*vthsq-8.300978857941995*rhor[2]*vthsq-8.300978857941995*rhol[2]*vthsq+2.143303524935281*rhor[1]*vthsq-2.143303524935281*rhol[1]*vthsq-3.273950213427199*rhor[3]*uvar0r[3]+2.766992952647331*rhor[2]*uvar0r[3]-2.143303524935281*rhor[1]*uvar0r[3]+1.237436867076458*rhor[0]*uvar0r[3]-3.273950213427199*rhol[3]*uvar0l[3]-2.766992952647331*rhol[2]*uvar0l[3]-2.143303524935281*rhol[1]*uvar0l[3]-1.237436867076458*rhol[0]*uvar0l[3]+2.766992952647331*uvar0r[2]*rhor[3]-2.143303524935281*uvar0r[1]*rhor[3]+1.237436867076458*uvar0r[0]*rhor[3]-2.766992952647331*uvar0l[2]*rhol[3]-2.143303524935281*uvar0l[1]*rhol[3]-1.237436867076458*uvar0l[0]*rhol[3]-2.338535866733713*rhor[2]*uvar0r[2]+1.81142209327368*rhor[1]*uvar0r[2]-1.045825033167594*rhor[0]*uvar0r[2]-2.338535866733713*rhol[2]*uvar0l[2]-1.81142209327368*rhol[1]*uvar0l[2]-1.045825033167594*rhol[0]*uvar0l[2]+1.81142209327368*uvar0r[1]*rhor[2]-1.045825033167594*uvar0r[0]*rhor[2]-1.81142209327368*uvar0l[1]*rhol[2]-1.045825033167594*uvar0l[0]*rhol[2]-1.403121520040228*rhor[1]*uvar0r[1]+0.8100925873009825*rhor[0]*uvar0r[1]-1.403121520040228*rhol[1]*uvar0l[1]-0.8100925873009825*rhol[0]*uvar0l[1]+0.8100925873009825*uvar0r[0]*rhor[1]-0.8100925873009825*uvar0l[0]*rhol[1]-0.4677071733467427*rhor[0]*uvar0r[0]-0.4677071733467427*rhol[0]*uvar0l[0]; 

  outrrhoux[0] += incr[0]*dxr1; 
  outrrhoux[1] += incr[1]*dxr1; 
  outrrhoux[2] += incr[2]*dxr1; 
  outrrhoux[3] += incr[3]*dxr1; 

  outlrhoux[0] += -1.0*incr[0]*dxl1; 
  outlrhoux[1] += incr[1]*dxl1; 
  outlrhoux[2] += -1.0*incr[2]*dxl1; 
  outlrhoux[3] += incr[3]*dxl1; 

 
} 
