#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void euler_pkpm_pressure_2x_ser_p1(const double *u_i, const double *bvar, const double *vlasov_pkpm_moms, const double *statevec, double* GKYL_RESTRICT p_ij) 
{ 
  // u_i: [ux, uy, uz], Fluid flow.
  // bvar: Magnetic field unit vector and tensor.
  // vlasov_pkpm_moms: [rho, p_parallel, q_parallel], Moments computed from kinetic equation in pkpm model.
  // statevec: [rho ux, rho uy, rho uz, energy], Fluid input state vector.
  // p_ij: Output pressure tensor, p_ij = (p_parallel - p_perp)bb + p_perp I.

  const double *rhoux = &statevec[0]; 
  const double *rhouy = &statevec[4]; 
  const double *rhouz = &statevec[8]; 
  const double *energy = &statevec[12]; 
  const double *ux = &u_i[0]; 
  const double *uy = &u_i[4]; 
  const double *uz = &u_i[8]; 

  // Parallel pressure is first component of pkpm moment array and unit tensor are last six components of bvar array.
  const double *p_parallel = &vlasov_pkpm_moms[4]; 
  const double *bxbx = &bvar[12]; 
  const double *bxby = &bvar[16]; 
  const double *bxbz = &bvar[20]; 
  const double *byby = &bvar[24]; 
  const double *bybz = &bvar[28]; 
  const double *bzbz = &bvar[32]; 

  double *Pxx = &p_ij[0]; 
  double *Pxy = &p_ij[4]; 
  double *Pxz = &p_ij[8]; 
  double *Pyy = &p_ij[12]; 
  double *Pyz = &p_ij[16]; 
  double *Pzz = &p_ij[20]; 

  double p_perp[4] = {0.0}; 
  p_perp[0] = (-0.25*rhouz[3]*uz[3])-0.25*rhouy[3]*uy[3]-0.25*rhoux[3]*ux[3]-0.25*rhouz[2]*uz[2]-0.25*rhouy[2]*uy[2]-0.25*rhoux[2]*ux[2]-0.25*rhouz[1]*uz[1]-0.25*rhouy[1]*uy[1]-0.25*rhoux[1]*ux[1]-0.25*rhouz[0]*uz[0]-0.25*rhouy[0]*uy[0]-0.25*rhoux[0]*ux[0]-0.5*p_parallel[0]+energy[0]; 
  p_perp[1] = (-0.25*rhouz[2]*uz[3])-0.25*rhouy[2]*uy[3]-0.25*rhoux[2]*ux[3]-0.25*uz[2]*rhouz[3]-0.25*uy[2]*rhouy[3]-0.25*ux[2]*rhoux[3]-0.25*rhouz[0]*uz[1]-0.25*rhouy[0]*uy[1]-0.25*rhoux[0]*ux[1]-0.25*uz[0]*rhouz[1]-0.25*uy[0]*rhouy[1]-0.25*ux[0]*rhoux[1]-0.5*p_parallel[1]+energy[1]; 
  p_perp[2] = (-0.25*rhouz[1]*uz[3])-0.25*rhouy[1]*uy[3]-0.25*rhoux[1]*ux[3]-0.25*uz[1]*rhouz[3]-0.25*uy[1]*rhouy[3]-0.25*ux[1]*rhoux[3]-0.25*rhouz[0]*uz[2]-0.25*rhouy[0]*uy[2]-0.25*rhoux[0]*ux[2]-0.25*uz[0]*rhouz[2]-0.25*uy[0]*rhouy[2]-0.25*ux[0]*rhoux[2]-0.5*p_parallel[2]+energy[2]; 
  p_perp[3] = (-0.25*rhouz[0]*uz[3])-0.25*rhouy[0]*uy[3]-0.25*rhoux[0]*ux[3]-0.25*uz[0]*rhouz[3]-0.25*uy[0]*rhouy[3]-0.25*ux[0]*rhoux[3]-0.5*p_parallel[3]+energy[3]-0.25*rhouz[1]*uz[2]-0.25*rhouy[1]*uy[2]-0.25*rhoux[1]*ux[2]-0.25*uz[1]*rhouz[2]-0.25*uy[1]*rhouy[2]-0.25*ux[1]*rhoux[2]; 

  Pxx[0] = (-0.5*bxbx[3]*p_perp[3])+0.5*bxbx[3]*p_parallel[3]-0.5*bxbx[2]*p_perp[2]+0.5*bxbx[2]*p_parallel[2]-0.5*bxbx[1]*p_perp[1]+0.5*bxbx[1]*p_parallel[1]-0.5*bxbx[0]*p_perp[0]+p_perp[0]+0.5*bxbx[0]*p_parallel[0]; 
  Pxx[1] = (-0.5*bxbx[2]*p_perp[3])+0.5*bxbx[2]*p_parallel[3]-0.5*p_perp[2]*bxbx[3]+0.5*p_parallel[2]*bxbx[3]-0.5*bxbx[0]*p_perp[1]+p_perp[1]+0.5*bxbx[0]*p_parallel[1]-0.5*p_perp[0]*bxbx[1]+0.5*p_parallel[0]*bxbx[1]; 
  Pxx[2] = (-0.5*bxbx[1]*p_perp[3])+0.5*bxbx[1]*p_parallel[3]-0.5*p_perp[1]*bxbx[3]+0.5*p_parallel[1]*bxbx[3]-0.5*bxbx[0]*p_perp[2]+p_perp[2]+0.5*bxbx[0]*p_parallel[2]-0.5*p_perp[0]*bxbx[2]+0.5*p_parallel[0]*bxbx[2]; 
  Pxx[3] = (-0.5*bxbx[0]*p_perp[3])+p_perp[3]+0.5*bxbx[0]*p_parallel[3]-0.5*p_perp[0]*bxbx[3]+0.5*p_parallel[0]*bxbx[3]-0.5*bxbx[1]*p_perp[2]+0.5*bxbx[1]*p_parallel[2]-0.5*p_perp[1]*bxbx[2]+0.5*p_parallel[1]*bxbx[2]; 

  Pxy[0] = (-0.5*bxby[3]*p_perp[3])+0.5*bxby[3]*p_parallel[3]-0.5*bxby[2]*p_perp[2]+0.5*bxby[2]*p_parallel[2]-0.5*bxby[1]*p_perp[1]+0.5*bxby[1]*p_parallel[1]-0.5*bxby[0]*p_perp[0]+0.5*bxby[0]*p_parallel[0]; 
  Pxy[1] = (-0.5*bxby[2]*p_perp[3])+0.5*bxby[2]*p_parallel[3]-0.5*p_perp[2]*bxby[3]+0.5*p_parallel[2]*bxby[3]-0.5*bxby[0]*p_perp[1]+0.5*bxby[0]*p_parallel[1]-0.5*p_perp[0]*bxby[1]+0.5*p_parallel[0]*bxby[1]; 
  Pxy[2] = (-0.5*bxby[1]*p_perp[3])+0.5*bxby[1]*p_parallel[3]-0.5*p_perp[1]*bxby[3]+0.5*p_parallel[1]*bxby[3]-0.5*bxby[0]*p_perp[2]+0.5*bxby[0]*p_parallel[2]-0.5*p_perp[0]*bxby[2]+0.5*p_parallel[0]*bxby[2]; 
  Pxy[3] = (-0.5*bxby[0]*p_perp[3])+0.5*bxby[0]*p_parallel[3]-0.5*p_perp[0]*bxby[3]+0.5*p_parallel[0]*bxby[3]-0.5*bxby[1]*p_perp[2]+0.5*bxby[1]*p_parallel[2]-0.5*p_perp[1]*bxby[2]+0.5*p_parallel[1]*bxby[2]; 

  Pxz[0] = (-0.5*bxbz[3]*p_perp[3])+0.5*bxbz[3]*p_parallel[3]-0.5*bxbz[2]*p_perp[2]+0.5*bxbz[2]*p_parallel[2]-0.5*bxbz[1]*p_perp[1]+0.5*bxbz[1]*p_parallel[1]-0.5*bxbz[0]*p_perp[0]+0.5*bxbz[0]*p_parallel[0]; 
  Pxz[1] = (-0.5*bxbz[2]*p_perp[3])+0.5*bxbz[2]*p_parallel[3]-0.5*p_perp[2]*bxbz[3]+0.5*p_parallel[2]*bxbz[3]-0.5*bxbz[0]*p_perp[1]+0.5*bxbz[0]*p_parallel[1]-0.5*p_perp[0]*bxbz[1]+0.5*p_parallel[0]*bxbz[1]; 
  Pxz[2] = (-0.5*bxbz[1]*p_perp[3])+0.5*bxbz[1]*p_parallel[3]-0.5*p_perp[1]*bxbz[3]+0.5*p_parallel[1]*bxbz[3]-0.5*bxbz[0]*p_perp[2]+0.5*bxbz[0]*p_parallel[2]-0.5*p_perp[0]*bxbz[2]+0.5*p_parallel[0]*bxbz[2]; 
  Pxz[3] = (-0.5*bxbz[0]*p_perp[3])+0.5*bxbz[0]*p_parallel[3]-0.5*p_perp[0]*bxbz[3]+0.5*p_parallel[0]*bxbz[3]-0.5*bxbz[1]*p_perp[2]+0.5*bxbz[1]*p_parallel[2]-0.5*p_perp[1]*bxbz[2]+0.5*p_parallel[1]*bxbz[2]; 

  Pyy[0] = (-0.5*byby[3]*p_perp[3])+0.5*byby[3]*p_parallel[3]-0.5*byby[2]*p_perp[2]+0.5*byby[2]*p_parallel[2]-0.5*byby[1]*p_perp[1]+0.5*byby[1]*p_parallel[1]-0.5*byby[0]*p_perp[0]+p_perp[0]+0.5*byby[0]*p_parallel[0]; 
  Pyy[1] = (-0.5*byby[2]*p_perp[3])+0.5*byby[2]*p_parallel[3]-0.5*p_perp[2]*byby[3]+0.5*p_parallel[2]*byby[3]-0.5*byby[0]*p_perp[1]+p_perp[1]+0.5*byby[0]*p_parallel[1]-0.5*p_perp[0]*byby[1]+0.5*p_parallel[0]*byby[1]; 
  Pyy[2] = (-0.5*byby[1]*p_perp[3])+0.5*byby[1]*p_parallel[3]-0.5*p_perp[1]*byby[3]+0.5*p_parallel[1]*byby[3]-0.5*byby[0]*p_perp[2]+p_perp[2]+0.5*byby[0]*p_parallel[2]-0.5*p_perp[0]*byby[2]+0.5*p_parallel[0]*byby[2]; 
  Pyy[3] = (-0.5*byby[0]*p_perp[3])+p_perp[3]+0.5*byby[0]*p_parallel[3]-0.5*p_perp[0]*byby[3]+0.5*p_parallel[0]*byby[3]-0.5*byby[1]*p_perp[2]+0.5*byby[1]*p_parallel[2]-0.5*p_perp[1]*byby[2]+0.5*p_parallel[1]*byby[2]; 

  Pyz[0] = (-0.5*bybz[3]*p_perp[3])+0.5*bybz[3]*p_parallel[3]-0.5*bybz[2]*p_perp[2]+0.5*bybz[2]*p_parallel[2]-0.5*bybz[1]*p_perp[1]+0.5*bybz[1]*p_parallel[1]-0.5*bybz[0]*p_perp[0]+0.5*bybz[0]*p_parallel[0]; 
  Pyz[1] = (-0.5*bybz[2]*p_perp[3])+0.5*bybz[2]*p_parallel[3]-0.5*p_perp[2]*bybz[3]+0.5*p_parallel[2]*bybz[3]-0.5*bybz[0]*p_perp[1]+0.5*bybz[0]*p_parallel[1]-0.5*p_perp[0]*bybz[1]+0.5*p_parallel[0]*bybz[1]; 
  Pyz[2] = (-0.5*bybz[1]*p_perp[3])+0.5*bybz[1]*p_parallel[3]-0.5*p_perp[1]*bybz[3]+0.5*p_parallel[1]*bybz[3]-0.5*bybz[0]*p_perp[2]+0.5*bybz[0]*p_parallel[2]-0.5*p_perp[0]*bybz[2]+0.5*p_parallel[0]*bybz[2]; 
  Pyz[3] = (-0.5*bybz[0]*p_perp[3])+0.5*bybz[0]*p_parallel[3]-0.5*p_perp[0]*bybz[3]+0.5*p_parallel[0]*bybz[3]-0.5*bybz[1]*p_perp[2]+0.5*bybz[1]*p_parallel[2]-0.5*p_perp[1]*bybz[2]+0.5*p_parallel[1]*bybz[2]; 

  Pzz[0] = (-0.5*bzbz[3]*p_perp[3])+0.5*bzbz[3]*p_parallel[3]-0.5*bzbz[2]*p_perp[2]+0.5*bzbz[2]*p_parallel[2]-0.5*bzbz[1]*p_perp[1]+0.5*bzbz[1]*p_parallel[1]-0.5*bzbz[0]*p_perp[0]+p_perp[0]+0.5*bzbz[0]*p_parallel[0]; 
  Pzz[1] = (-0.5*bzbz[2]*p_perp[3])+0.5*bzbz[2]*p_parallel[3]-0.5*p_perp[2]*bzbz[3]+0.5*p_parallel[2]*bzbz[3]-0.5*bzbz[0]*p_perp[1]+p_perp[1]+0.5*bzbz[0]*p_parallel[1]-0.5*p_perp[0]*bzbz[1]+0.5*p_parallel[0]*bzbz[1]; 
  Pzz[2] = (-0.5*bzbz[1]*p_perp[3])+0.5*bzbz[1]*p_parallel[3]-0.5*p_perp[1]*bzbz[3]+0.5*p_parallel[1]*bzbz[3]-0.5*bzbz[0]*p_perp[2]+p_perp[2]+0.5*bzbz[0]*p_parallel[2]-0.5*p_perp[0]*bzbz[2]+0.5*p_parallel[0]*bzbz[2]; 
  Pzz[3] = (-0.5*bzbz[0]*p_perp[3])+p_perp[3]+0.5*bzbz[0]*p_parallel[3]-0.5*p_perp[0]*bzbz[3]+0.5*p_parallel[0]*bzbz[3]-0.5*bzbz[1]*p_perp[2]+0.5*bzbz[1]*p_parallel[2]-0.5*p_perp[1]*bzbz[2]+0.5*p_parallel[1]*bzbz[2]; 
} 
