#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void euler_pkpm_pressure_3x_ser_p1(const double *u_i, const double *bvar, const double *vlasov_pkpm_moms, const double *statevec, double* GKYL_RESTRICT p_ij) 
{ 
  // u_i: [ux, uy, uz], Fluid flow.
  // bvar: Magnetic field unit vector and tensor.
  // vlasov_pkpm_moms: [rho, p_parallel, q_parallel], Moments computed from kinetic equation in pkpm model.
  // statevec: [rho ux, rho uy, rho uz, energy], Fluid input state vector.
  // p_ij: Output pressure tensor, p_ij = (p_parallel - p_perp)bb + p_perp I.

  const double *rhoux = &statevec[0]; 
  const double *rhouy = &statevec[8]; 
  const double *rhouz = &statevec[16]; 
  const double *energy = &statevec[24]; 
  const double *ux = &u_i[0]; 
  const double *uy = &u_i[8]; 
  const double *uz = &u_i[16]; 

  // Parallel pressure is first component of pkpm moment array and unit tensor are last six components of bvar array.
  const double *p_parallel = &vlasov_pkpm_moms[8]; 
  const double *bxbx = &bvar[24]; 
  const double *bxby = &bvar[32]; 
  const double *bxbz = &bvar[40]; 
  const double *byby = &bvar[48]; 
  const double *bybz = &bvar[56]; 
  const double *bzbz = &bvar[64]; 

  double *Pxx = &p_ij[0]; 
  double *Pxy = &p_ij[8]; 
  double *Pxz = &p_ij[16]; 
  double *Pyy = &p_ij[24]; 
  double *Pyz = &p_ij[32]; 
  double *Pzz = &p_ij[40]; 

  double p_perp[8] = {0.0}; 
  p_perp[0] = (-0.1767766952966368*rhouz[7]*uz[7])-0.1767766952966368*rhouy[7]*uy[7]-0.1767766952966368*rhoux[7]*ux[7]-0.1767766952966368*rhouz[6]*uz[6]-0.1767766952966368*rhouy[6]*uy[6]-0.1767766952966368*rhoux[6]*ux[6]-0.1767766952966368*rhouz[5]*uz[5]-0.1767766952966368*rhouy[5]*uy[5]-0.1767766952966368*rhoux[5]*ux[5]-0.1767766952966368*rhouz[4]*uz[4]-0.1767766952966368*rhouy[4]*uy[4]-0.1767766952966368*rhoux[4]*ux[4]-0.1767766952966368*rhouz[3]*uz[3]-0.1767766952966368*rhouy[3]*uy[3]-0.1767766952966368*rhoux[3]*ux[3]-0.1767766952966368*rhouz[2]*uz[2]-0.1767766952966368*rhouy[2]*uy[2]-0.1767766952966368*rhoux[2]*ux[2]-0.1767766952966368*rhouz[1]*uz[1]-0.1767766952966368*rhouy[1]*uy[1]-0.1767766952966368*rhoux[1]*ux[1]-0.1767766952966368*rhouz[0]*uz[0]-0.1767766952966368*rhouy[0]*uy[0]-0.1767766952966368*rhoux[0]*ux[0]-0.5*p_parallel[0]+energy[0]; 
  p_perp[1] = (-0.1767766952966368*rhouz[6]*uz[7])-0.1767766952966368*rhouy[6]*uy[7]-0.1767766952966368*rhoux[6]*ux[7]-0.1767766952966368*uz[6]*rhouz[7]-0.1767766952966368*uy[6]*rhouy[7]-0.1767766952966368*ux[6]*rhoux[7]-0.1767766952966368*rhouz[3]*uz[5]-0.1767766952966368*rhouy[3]*uy[5]-0.1767766952966368*rhoux[3]*ux[5]-0.1767766952966368*uz[3]*rhouz[5]-0.1767766952966368*uy[3]*rhouy[5]-0.1767766952966368*ux[3]*rhoux[5]-0.1767766952966368*rhouz[2]*uz[4]-0.1767766952966368*rhouy[2]*uy[4]-0.1767766952966368*rhoux[2]*ux[4]-0.1767766952966368*uz[2]*rhouz[4]-0.1767766952966368*uy[2]*rhouy[4]-0.1767766952966368*ux[2]*rhoux[4]-0.1767766952966368*rhouz[0]*uz[1]-0.1767766952966368*rhouy[0]*uy[1]-0.1767766952966368*rhoux[0]*ux[1]-0.1767766952966368*uz[0]*rhouz[1]-0.1767766952966368*uy[0]*rhouy[1]-0.1767766952966368*ux[0]*rhoux[1]-0.5*p_parallel[1]+energy[1]; 
  p_perp[2] = (-0.1767766952966368*rhouz[5]*uz[7])-0.1767766952966368*rhouy[5]*uy[7]-0.1767766952966368*rhoux[5]*ux[7]-0.1767766952966368*uz[5]*rhouz[7]-0.1767766952966368*uy[5]*rhouy[7]-0.1767766952966368*ux[5]*rhoux[7]-0.1767766952966368*rhouz[3]*uz[6]-0.1767766952966368*rhouy[3]*uy[6]-0.1767766952966368*rhoux[3]*ux[6]-0.1767766952966368*uz[3]*rhouz[6]-0.1767766952966368*uy[3]*rhouy[6]-0.1767766952966368*ux[3]*rhoux[6]-0.1767766952966368*rhouz[1]*uz[4]-0.1767766952966368*rhouy[1]*uy[4]-0.1767766952966368*rhoux[1]*ux[4]-0.1767766952966368*uz[1]*rhouz[4]-0.1767766952966368*uy[1]*rhouy[4]-0.1767766952966368*ux[1]*rhoux[4]-0.1767766952966368*rhouz[0]*uz[2]-0.1767766952966368*rhouy[0]*uy[2]-0.1767766952966368*rhoux[0]*ux[2]-0.1767766952966368*uz[0]*rhouz[2]-0.1767766952966368*uy[0]*rhouy[2]-0.1767766952966368*ux[0]*rhoux[2]-0.5*p_parallel[2]+energy[2]; 
  p_perp[3] = (-0.1767766952966368*rhouz[4]*uz[7])-0.1767766952966368*rhouy[4]*uy[7]-0.1767766952966368*rhoux[4]*ux[7]-0.1767766952966368*uz[4]*rhouz[7]-0.1767766952966368*uy[4]*rhouy[7]-0.1767766952966368*ux[4]*rhoux[7]-0.1767766952966368*rhouz[2]*uz[6]-0.1767766952966368*rhouy[2]*uy[6]-0.1767766952966368*rhoux[2]*ux[6]-0.1767766952966368*uz[2]*rhouz[6]-0.1767766952966368*uy[2]*rhouy[6]-0.1767766952966368*ux[2]*rhoux[6]-0.1767766952966368*rhouz[1]*uz[5]-0.1767766952966368*rhouy[1]*uy[5]-0.1767766952966368*rhoux[1]*ux[5]-0.1767766952966368*uz[1]*rhouz[5]-0.1767766952966368*uy[1]*rhouy[5]-0.1767766952966368*ux[1]*rhoux[5]-0.1767766952966368*rhouz[0]*uz[3]-0.1767766952966368*rhouy[0]*uy[3]-0.1767766952966368*rhoux[0]*ux[3]-0.1767766952966368*uz[0]*rhouz[3]-0.1767766952966368*uy[0]*rhouy[3]-0.1767766952966368*ux[0]*rhoux[3]-0.5*p_parallel[3]+energy[3]; 
  p_perp[4] = (-0.1767766952966368*rhouz[3]*uz[7])-0.1767766952966368*rhouy[3]*uy[7]-0.1767766952966368*rhoux[3]*ux[7]-0.1767766952966368*uz[3]*rhouz[7]-0.1767766952966368*uy[3]*rhouy[7]-0.1767766952966368*ux[3]*rhoux[7]-0.1767766952966368*rhouz[5]*uz[6]-0.1767766952966368*rhouy[5]*uy[6]-0.1767766952966368*rhoux[5]*ux[6]-0.1767766952966368*uz[5]*rhouz[6]-0.1767766952966368*uy[5]*rhouy[6]-0.1767766952966368*ux[5]*rhoux[6]-0.1767766952966368*rhouz[0]*uz[4]-0.1767766952966368*rhouy[0]*uy[4]-0.1767766952966368*rhoux[0]*ux[4]-0.1767766952966368*uz[0]*rhouz[4]-0.1767766952966368*uy[0]*rhouy[4]-0.1767766952966368*ux[0]*rhoux[4]-0.5*p_parallel[4]+energy[4]-0.1767766952966368*rhouz[1]*uz[2]-0.1767766952966368*rhouy[1]*uy[2]-0.1767766952966368*rhoux[1]*ux[2]-0.1767766952966368*uz[1]*rhouz[2]-0.1767766952966368*uy[1]*rhouy[2]-0.1767766952966368*ux[1]*rhoux[2]; 
  p_perp[5] = (-0.1767766952966368*rhouz[2]*uz[7])-0.1767766952966368*rhouy[2]*uy[7]-0.1767766952966368*rhoux[2]*ux[7]-0.1767766952966368*uz[2]*rhouz[7]-0.1767766952966368*uy[2]*rhouy[7]-0.1767766952966368*ux[2]*rhoux[7]-0.1767766952966368*rhouz[4]*uz[6]-0.1767766952966368*rhouy[4]*uy[6]-0.1767766952966368*rhoux[4]*ux[6]-0.1767766952966368*uz[4]*rhouz[6]-0.1767766952966368*uy[4]*rhouy[6]-0.1767766952966368*ux[4]*rhoux[6]-0.1767766952966368*rhouz[0]*uz[5]-0.1767766952966368*rhouy[0]*uy[5]-0.1767766952966368*rhoux[0]*ux[5]-0.1767766952966368*uz[0]*rhouz[5]-0.1767766952966368*uy[0]*rhouy[5]-0.1767766952966368*ux[0]*rhoux[5]-0.5*p_parallel[5]+energy[5]-0.1767766952966368*rhouz[1]*uz[3]-0.1767766952966368*rhouy[1]*uy[3]-0.1767766952966368*rhoux[1]*ux[3]-0.1767766952966368*uz[1]*rhouz[3]-0.1767766952966368*uy[1]*rhouy[3]-0.1767766952966368*ux[1]*rhoux[3]; 
  p_perp[6] = (-0.1767766952966368*rhouz[1]*uz[7])-0.1767766952966368*rhouy[1]*uy[7]-0.1767766952966368*rhoux[1]*ux[7]-0.1767766952966368*uz[1]*rhouz[7]-0.1767766952966368*uy[1]*rhouy[7]-0.1767766952966368*ux[1]*rhoux[7]-0.1767766952966368*rhouz[0]*uz[6]-0.1767766952966368*rhouy[0]*uy[6]-0.1767766952966368*rhoux[0]*ux[6]-0.1767766952966368*uz[0]*rhouz[6]-0.1767766952966368*uy[0]*rhouy[6]-0.1767766952966368*ux[0]*rhoux[6]-0.5*p_parallel[6]+energy[6]-0.1767766952966368*rhouz[4]*uz[5]-0.1767766952966368*rhouy[4]*uy[5]-0.1767766952966368*rhoux[4]*ux[5]-0.1767766952966368*uz[4]*rhouz[5]-0.1767766952966368*uy[4]*rhouy[5]-0.1767766952966368*ux[4]*rhoux[5]-0.1767766952966368*rhouz[2]*uz[3]-0.1767766952966368*rhouy[2]*uy[3]-0.1767766952966368*rhoux[2]*ux[3]-0.1767766952966368*uz[2]*rhouz[3]-0.1767766952966368*uy[2]*rhouy[3]-0.1767766952966368*ux[2]*rhoux[3]; 
  p_perp[7] = (-0.1767766952966368*rhouz[0]*uz[7])-0.1767766952966368*rhouy[0]*uy[7]-0.1767766952966368*rhoux[0]*ux[7]-0.1767766952966368*uz[0]*rhouz[7]-0.1767766952966368*uy[0]*rhouy[7]-0.1767766952966368*ux[0]*rhoux[7]-0.5*p_parallel[7]+energy[7]-0.1767766952966368*rhouz[1]*uz[6]-0.1767766952966368*rhouy[1]*uy[6]-0.1767766952966368*rhoux[1]*ux[6]-0.1767766952966368*uz[1]*rhouz[6]-0.1767766952966368*uy[1]*rhouy[6]-0.1767766952966368*ux[1]*rhoux[6]-0.1767766952966368*rhouz[2]*uz[5]-0.1767766952966368*rhouy[2]*uy[5]-0.1767766952966368*rhoux[2]*ux[5]-0.1767766952966368*uz[2]*rhouz[5]-0.1767766952966368*uy[2]*rhouy[5]-0.1767766952966368*ux[2]*rhoux[5]-0.1767766952966368*rhouz[3]*uz[4]-0.1767766952966368*rhouy[3]*uy[4]-0.1767766952966368*rhoux[3]*ux[4]-0.1767766952966368*uz[3]*rhouz[4]-0.1767766952966368*uy[3]*rhouy[4]-0.1767766952966368*ux[3]*rhoux[4]; 

  Pxx[0] = (-0.3535533905932737*bxbx[7]*p_perp[7])+0.3535533905932737*bxbx[7]*p_parallel[7]-0.3535533905932737*bxbx[6]*p_perp[6]+0.3535533905932737*bxbx[6]*p_parallel[6]-0.3535533905932737*bxbx[5]*p_perp[5]+0.3535533905932737*bxbx[5]*p_parallel[5]-0.3535533905932737*bxbx[4]*p_perp[4]+0.3535533905932737*bxbx[4]*p_parallel[4]-0.3535533905932737*bxbx[3]*p_perp[3]+0.3535533905932737*bxbx[3]*p_parallel[3]-0.3535533905932737*bxbx[2]*p_perp[2]+0.3535533905932737*bxbx[2]*p_parallel[2]-0.3535533905932737*bxbx[1]*p_perp[1]+0.3535533905932737*bxbx[1]*p_parallel[1]-0.3535533905932737*bxbx[0]*p_perp[0]+p_perp[0]+0.3535533905932737*bxbx[0]*p_parallel[0]; 
  Pxx[1] = (-0.3535533905932737*bxbx[6]*p_perp[7])+0.3535533905932737*bxbx[6]*p_parallel[7]-0.3535533905932737*p_perp[6]*bxbx[7]+0.3535533905932737*p_parallel[6]*bxbx[7]-0.3535533905932737*bxbx[3]*p_perp[5]+0.3535533905932737*bxbx[3]*p_parallel[5]-0.3535533905932737*p_perp[3]*bxbx[5]+0.3535533905932737*p_parallel[3]*bxbx[5]-0.3535533905932737*bxbx[2]*p_perp[4]+0.3535533905932737*bxbx[2]*p_parallel[4]-0.3535533905932737*p_perp[2]*bxbx[4]+0.3535533905932737*p_parallel[2]*bxbx[4]-0.3535533905932737*bxbx[0]*p_perp[1]+p_perp[1]+0.3535533905932737*bxbx[0]*p_parallel[1]-0.3535533905932737*p_perp[0]*bxbx[1]+0.3535533905932737*p_parallel[0]*bxbx[1]; 
  Pxx[2] = (-0.3535533905932737*bxbx[5]*p_perp[7])+0.3535533905932737*bxbx[5]*p_parallel[7]-0.3535533905932737*p_perp[5]*bxbx[7]+0.3535533905932737*p_parallel[5]*bxbx[7]-0.3535533905932737*bxbx[3]*p_perp[6]+0.3535533905932737*bxbx[3]*p_parallel[6]-0.3535533905932737*p_perp[3]*bxbx[6]+0.3535533905932737*p_parallel[3]*bxbx[6]-0.3535533905932737*bxbx[1]*p_perp[4]+0.3535533905932737*bxbx[1]*p_parallel[4]-0.3535533905932737*p_perp[1]*bxbx[4]+0.3535533905932737*p_parallel[1]*bxbx[4]-0.3535533905932737*bxbx[0]*p_perp[2]+p_perp[2]+0.3535533905932737*bxbx[0]*p_parallel[2]-0.3535533905932737*p_perp[0]*bxbx[2]+0.3535533905932737*p_parallel[0]*bxbx[2]; 
  Pxx[3] = (-0.3535533905932737*bxbx[4]*p_perp[7])+0.3535533905932737*bxbx[4]*p_parallel[7]-0.3535533905932737*p_perp[4]*bxbx[7]+0.3535533905932737*p_parallel[4]*bxbx[7]-0.3535533905932737*bxbx[2]*p_perp[6]+0.3535533905932737*bxbx[2]*p_parallel[6]-0.3535533905932737*p_perp[2]*bxbx[6]+0.3535533905932737*p_parallel[2]*bxbx[6]-0.3535533905932737*bxbx[1]*p_perp[5]+0.3535533905932737*bxbx[1]*p_parallel[5]-0.3535533905932737*p_perp[1]*bxbx[5]+0.3535533905932737*p_parallel[1]*bxbx[5]-0.3535533905932737*bxbx[0]*p_perp[3]+p_perp[3]+0.3535533905932737*bxbx[0]*p_parallel[3]-0.3535533905932737*p_perp[0]*bxbx[3]+0.3535533905932737*p_parallel[0]*bxbx[3]; 
  Pxx[4] = (-0.3535533905932737*bxbx[3]*p_perp[7])+0.3535533905932737*bxbx[3]*p_parallel[7]-0.3535533905932737*p_perp[3]*bxbx[7]+0.3535533905932737*p_parallel[3]*bxbx[7]-0.3535533905932737*bxbx[5]*p_perp[6]+0.3535533905932737*bxbx[5]*p_parallel[6]-0.3535533905932737*p_perp[5]*bxbx[6]+0.3535533905932737*p_parallel[5]*bxbx[6]-0.3535533905932737*bxbx[0]*p_perp[4]+p_perp[4]+0.3535533905932737*bxbx[0]*p_parallel[4]-0.3535533905932737*p_perp[0]*bxbx[4]+0.3535533905932737*p_parallel[0]*bxbx[4]-0.3535533905932737*bxbx[1]*p_perp[2]+0.3535533905932737*bxbx[1]*p_parallel[2]-0.3535533905932737*p_perp[1]*bxbx[2]+0.3535533905932737*p_parallel[1]*bxbx[2]; 
  Pxx[5] = (-0.3535533905932737*bxbx[2]*p_perp[7])+0.3535533905932737*bxbx[2]*p_parallel[7]-0.3535533905932737*p_perp[2]*bxbx[7]+0.3535533905932737*p_parallel[2]*bxbx[7]-0.3535533905932737*bxbx[4]*p_perp[6]+0.3535533905932737*bxbx[4]*p_parallel[6]-0.3535533905932737*p_perp[4]*bxbx[6]+0.3535533905932737*p_parallel[4]*bxbx[6]-0.3535533905932737*bxbx[0]*p_perp[5]+p_perp[5]+0.3535533905932737*bxbx[0]*p_parallel[5]-0.3535533905932737*p_perp[0]*bxbx[5]+0.3535533905932737*p_parallel[0]*bxbx[5]-0.3535533905932737*bxbx[1]*p_perp[3]+0.3535533905932737*bxbx[1]*p_parallel[3]-0.3535533905932737*p_perp[1]*bxbx[3]+0.3535533905932737*p_parallel[1]*bxbx[3]; 
  Pxx[6] = (-0.3535533905932737*bxbx[1]*p_perp[7])+0.3535533905932737*bxbx[1]*p_parallel[7]-0.3535533905932737*p_perp[1]*bxbx[7]+0.3535533905932737*p_parallel[1]*bxbx[7]-0.3535533905932737*bxbx[0]*p_perp[6]+p_perp[6]+0.3535533905932737*bxbx[0]*p_parallel[6]-0.3535533905932737*p_perp[0]*bxbx[6]+0.3535533905932737*p_parallel[0]*bxbx[6]-0.3535533905932737*bxbx[4]*p_perp[5]+0.3535533905932737*bxbx[4]*p_parallel[5]-0.3535533905932737*p_perp[4]*bxbx[5]+0.3535533905932737*p_parallel[4]*bxbx[5]-0.3535533905932737*bxbx[2]*p_perp[3]+0.3535533905932737*bxbx[2]*p_parallel[3]-0.3535533905932737*p_perp[2]*bxbx[3]+0.3535533905932737*p_parallel[2]*bxbx[3]; 
  Pxx[7] = (-0.3535533905932737*bxbx[0]*p_perp[7])+p_perp[7]+0.3535533905932737*bxbx[0]*p_parallel[7]-0.3535533905932737*p_perp[0]*bxbx[7]+0.3535533905932737*p_parallel[0]*bxbx[7]-0.3535533905932737*bxbx[1]*p_perp[6]+0.3535533905932737*bxbx[1]*p_parallel[6]-0.3535533905932737*p_perp[1]*bxbx[6]+0.3535533905932737*p_parallel[1]*bxbx[6]-0.3535533905932737*bxbx[2]*p_perp[5]+0.3535533905932737*bxbx[2]*p_parallel[5]-0.3535533905932737*p_perp[2]*bxbx[5]+0.3535533905932737*p_parallel[2]*bxbx[5]-0.3535533905932737*bxbx[3]*p_perp[4]+0.3535533905932737*bxbx[3]*p_parallel[4]-0.3535533905932737*p_perp[3]*bxbx[4]+0.3535533905932737*p_parallel[3]*bxbx[4]; 

  Pxy[0] = (-0.3535533905932737*bxby[7]*p_perp[7])+0.3535533905932737*bxby[7]*p_parallel[7]-0.3535533905932737*bxby[6]*p_perp[6]+0.3535533905932737*bxby[6]*p_parallel[6]-0.3535533905932737*bxby[5]*p_perp[5]+0.3535533905932737*bxby[5]*p_parallel[5]-0.3535533905932737*bxby[4]*p_perp[4]+0.3535533905932737*bxby[4]*p_parallel[4]-0.3535533905932737*bxby[3]*p_perp[3]+0.3535533905932737*bxby[3]*p_parallel[3]-0.3535533905932737*bxby[2]*p_perp[2]+0.3535533905932737*bxby[2]*p_parallel[2]-0.3535533905932737*bxby[1]*p_perp[1]+0.3535533905932737*bxby[1]*p_parallel[1]-0.3535533905932737*bxby[0]*p_perp[0]+0.3535533905932737*bxby[0]*p_parallel[0]; 
  Pxy[1] = (-0.3535533905932737*bxby[6]*p_perp[7])+0.3535533905932737*bxby[6]*p_parallel[7]-0.3535533905932737*p_perp[6]*bxby[7]+0.3535533905932737*p_parallel[6]*bxby[7]-0.3535533905932737*bxby[3]*p_perp[5]+0.3535533905932737*bxby[3]*p_parallel[5]-0.3535533905932737*p_perp[3]*bxby[5]+0.3535533905932737*p_parallel[3]*bxby[5]-0.3535533905932737*bxby[2]*p_perp[4]+0.3535533905932737*bxby[2]*p_parallel[4]-0.3535533905932737*p_perp[2]*bxby[4]+0.3535533905932737*p_parallel[2]*bxby[4]-0.3535533905932737*bxby[0]*p_perp[1]+0.3535533905932737*bxby[0]*p_parallel[1]-0.3535533905932737*p_perp[0]*bxby[1]+0.3535533905932737*p_parallel[0]*bxby[1]; 
  Pxy[2] = (-0.3535533905932737*bxby[5]*p_perp[7])+0.3535533905932737*bxby[5]*p_parallel[7]-0.3535533905932737*p_perp[5]*bxby[7]+0.3535533905932737*p_parallel[5]*bxby[7]-0.3535533905932737*bxby[3]*p_perp[6]+0.3535533905932737*bxby[3]*p_parallel[6]-0.3535533905932737*p_perp[3]*bxby[6]+0.3535533905932737*p_parallel[3]*bxby[6]-0.3535533905932737*bxby[1]*p_perp[4]+0.3535533905932737*bxby[1]*p_parallel[4]-0.3535533905932737*p_perp[1]*bxby[4]+0.3535533905932737*p_parallel[1]*bxby[4]-0.3535533905932737*bxby[0]*p_perp[2]+0.3535533905932737*bxby[0]*p_parallel[2]-0.3535533905932737*p_perp[0]*bxby[2]+0.3535533905932737*p_parallel[0]*bxby[2]; 
  Pxy[3] = (-0.3535533905932737*bxby[4]*p_perp[7])+0.3535533905932737*bxby[4]*p_parallel[7]-0.3535533905932737*p_perp[4]*bxby[7]+0.3535533905932737*p_parallel[4]*bxby[7]-0.3535533905932737*bxby[2]*p_perp[6]+0.3535533905932737*bxby[2]*p_parallel[6]-0.3535533905932737*p_perp[2]*bxby[6]+0.3535533905932737*p_parallel[2]*bxby[6]-0.3535533905932737*bxby[1]*p_perp[5]+0.3535533905932737*bxby[1]*p_parallel[5]-0.3535533905932737*p_perp[1]*bxby[5]+0.3535533905932737*p_parallel[1]*bxby[5]-0.3535533905932737*bxby[0]*p_perp[3]+0.3535533905932737*bxby[0]*p_parallel[3]-0.3535533905932737*p_perp[0]*bxby[3]+0.3535533905932737*p_parallel[0]*bxby[3]; 
  Pxy[4] = (-0.3535533905932737*bxby[3]*p_perp[7])+0.3535533905932737*bxby[3]*p_parallel[7]-0.3535533905932737*p_perp[3]*bxby[7]+0.3535533905932737*p_parallel[3]*bxby[7]-0.3535533905932737*bxby[5]*p_perp[6]+0.3535533905932737*bxby[5]*p_parallel[6]-0.3535533905932737*p_perp[5]*bxby[6]+0.3535533905932737*p_parallel[5]*bxby[6]-0.3535533905932737*bxby[0]*p_perp[4]+0.3535533905932737*bxby[0]*p_parallel[4]-0.3535533905932737*p_perp[0]*bxby[4]+0.3535533905932737*p_parallel[0]*bxby[4]-0.3535533905932737*bxby[1]*p_perp[2]+0.3535533905932737*bxby[1]*p_parallel[2]-0.3535533905932737*p_perp[1]*bxby[2]+0.3535533905932737*p_parallel[1]*bxby[2]; 
  Pxy[5] = (-0.3535533905932737*bxby[2]*p_perp[7])+0.3535533905932737*bxby[2]*p_parallel[7]-0.3535533905932737*p_perp[2]*bxby[7]+0.3535533905932737*p_parallel[2]*bxby[7]-0.3535533905932737*bxby[4]*p_perp[6]+0.3535533905932737*bxby[4]*p_parallel[6]-0.3535533905932737*p_perp[4]*bxby[6]+0.3535533905932737*p_parallel[4]*bxby[6]-0.3535533905932737*bxby[0]*p_perp[5]+0.3535533905932737*bxby[0]*p_parallel[5]-0.3535533905932737*p_perp[0]*bxby[5]+0.3535533905932737*p_parallel[0]*bxby[5]-0.3535533905932737*bxby[1]*p_perp[3]+0.3535533905932737*bxby[1]*p_parallel[3]-0.3535533905932737*p_perp[1]*bxby[3]+0.3535533905932737*p_parallel[1]*bxby[3]; 
  Pxy[6] = (-0.3535533905932737*bxby[1]*p_perp[7])+0.3535533905932737*bxby[1]*p_parallel[7]-0.3535533905932737*p_perp[1]*bxby[7]+0.3535533905932737*p_parallel[1]*bxby[7]-0.3535533905932737*bxby[0]*p_perp[6]+0.3535533905932737*bxby[0]*p_parallel[6]-0.3535533905932737*p_perp[0]*bxby[6]+0.3535533905932737*p_parallel[0]*bxby[6]-0.3535533905932737*bxby[4]*p_perp[5]+0.3535533905932737*bxby[4]*p_parallel[5]-0.3535533905932737*p_perp[4]*bxby[5]+0.3535533905932737*p_parallel[4]*bxby[5]-0.3535533905932737*bxby[2]*p_perp[3]+0.3535533905932737*bxby[2]*p_parallel[3]-0.3535533905932737*p_perp[2]*bxby[3]+0.3535533905932737*p_parallel[2]*bxby[3]; 
  Pxy[7] = (-0.3535533905932737*bxby[0]*p_perp[7])+0.3535533905932737*bxby[0]*p_parallel[7]-0.3535533905932737*p_perp[0]*bxby[7]+0.3535533905932737*p_parallel[0]*bxby[7]-0.3535533905932737*bxby[1]*p_perp[6]+0.3535533905932737*bxby[1]*p_parallel[6]-0.3535533905932737*p_perp[1]*bxby[6]+0.3535533905932737*p_parallel[1]*bxby[6]-0.3535533905932737*bxby[2]*p_perp[5]+0.3535533905932737*bxby[2]*p_parallel[5]-0.3535533905932737*p_perp[2]*bxby[5]+0.3535533905932737*p_parallel[2]*bxby[5]-0.3535533905932737*bxby[3]*p_perp[4]+0.3535533905932737*bxby[3]*p_parallel[4]-0.3535533905932737*p_perp[3]*bxby[4]+0.3535533905932737*p_parallel[3]*bxby[4]; 

  Pxz[0] = (-0.3535533905932737*bxbz[7]*p_perp[7])+0.3535533905932737*bxbz[7]*p_parallel[7]-0.3535533905932737*bxbz[6]*p_perp[6]+0.3535533905932737*bxbz[6]*p_parallel[6]-0.3535533905932737*bxbz[5]*p_perp[5]+0.3535533905932737*bxbz[5]*p_parallel[5]-0.3535533905932737*bxbz[4]*p_perp[4]+0.3535533905932737*bxbz[4]*p_parallel[4]-0.3535533905932737*bxbz[3]*p_perp[3]+0.3535533905932737*bxbz[3]*p_parallel[3]-0.3535533905932737*bxbz[2]*p_perp[2]+0.3535533905932737*bxbz[2]*p_parallel[2]-0.3535533905932737*bxbz[1]*p_perp[1]+0.3535533905932737*bxbz[1]*p_parallel[1]-0.3535533905932737*bxbz[0]*p_perp[0]+0.3535533905932737*bxbz[0]*p_parallel[0]; 
  Pxz[1] = (-0.3535533905932737*bxbz[6]*p_perp[7])+0.3535533905932737*bxbz[6]*p_parallel[7]-0.3535533905932737*p_perp[6]*bxbz[7]+0.3535533905932737*p_parallel[6]*bxbz[7]-0.3535533905932737*bxbz[3]*p_perp[5]+0.3535533905932737*bxbz[3]*p_parallel[5]-0.3535533905932737*p_perp[3]*bxbz[5]+0.3535533905932737*p_parallel[3]*bxbz[5]-0.3535533905932737*bxbz[2]*p_perp[4]+0.3535533905932737*bxbz[2]*p_parallel[4]-0.3535533905932737*p_perp[2]*bxbz[4]+0.3535533905932737*p_parallel[2]*bxbz[4]-0.3535533905932737*bxbz[0]*p_perp[1]+0.3535533905932737*bxbz[0]*p_parallel[1]-0.3535533905932737*p_perp[0]*bxbz[1]+0.3535533905932737*p_parallel[0]*bxbz[1]; 
  Pxz[2] = (-0.3535533905932737*bxbz[5]*p_perp[7])+0.3535533905932737*bxbz[5]*p_parallel[7]-0.3535533905932737*p_perp[5]*bxbz[7]+0.3535533905932737*p_parallel[5]*bxbz[7]-0.3535533905932737*bxbz[3]*p_perp[6]+0.3535533905932737*bxbz[3]*p_parallel[6]-0.3535533905932737*p_perp[3]*bxbz[6]+0.3535533905932737*p_parallel[3]*bxbz[6]-0.3535533905932737*bxbz[1]*p_perp[4]+0.3535533905932737*bxbz[1]*p_parallel[4]-0.3535533905932737*p_perp[1]*bxbz[4]+0.3535533905932737*p_parallel[1]*bxbz[4]-0.3535533905932737*bxbz[0]*p_perp[2]+0.3535533905932737*bxbz[0]*p_parallel[2]-0.3535533905932737*p_perp[0]*bxbz[2]+0.3535533905932737*p_parallel[0]*bxbz[2]; 
  Pxz[3] = (-0.3535533905932737*bxbz[4]*p_perp[7])+0.3535533905932737*bxbz[4]*p_parallel[7]-0.3535533905932737*p_perp[4]*bxbz[7]+0.3535533905932737*p_parallel[4]*bxbz[7]-0.3535533905932737*bxbz[2]*p_perp[6]+0.3535533905932737*bxbz[2]*p_parallel[6]-0.3535533905932737*p_perp[2]*bxbz[6]+0.3535533905932737*p_parallel[2]*bxbz[6]-0.3535533905932737*bxbz[1]*p_perp[5]+0.3535533905932737*bxbz[1]*p_parallel[5]-0.3535533905932737*p_perp[1]*bxbz[5]+0.3535533905932737*p_parallel[1]*bxbz[5]-0.3535533905932737*bxbz[0]*p_perp[3]+0.3535533905932737*bxbz[0]*p_parallel[3]-0.3535533905932737*p_perp[0]*bxbz[3]+0.3535533905932737*p_parallel[0]*bxbz[3]; 
  Pxz[4] = (-0.3535533905932737*bxbz[3]*p_perp[7])+0.3535533905932737*bxbz[3]*p_parallel[7]-0.3535533905932737*p_perp[3]*bxbz[7]+0.3535533905932737*p_parallel[3]*bxbz[7]-0.3535533905932737*bxbz[5]*p_perp[6]+0.3535533905932737*bxbz[5]*p_parallel[6]-0.3535533905932737*p_perp[5]*bxbz[6]+0.3535533905932737*p_parallel[5]*bxbz[6]-0.3535533905932737*bxbz[0]*p_perp[4]+0.3535533905932737*bxbz[0]*p_parallel[4]-0.3535533905932737*p_perp[0]*bxbz[4]+0.3535533905932737*p_parallel[0]*bxbz[4]-0.3535533905932737*bxbz[1]*p_perp[2]+0.3535533905932737*bxbz[1]*p_parallel[2]-0.3535533905932737*p_perp[1]*bxbz[2]+0.3535533905932737*p_parallel[1]*bxbz[2]; 
  Pxz[5] = (-0.3535533905932737*bxbz[2]*p_perp[7])+0.3535533905932737*bxbz[2]*p_parallel[7]-0.3535533905932737*p_perp[2]*bxbz[7]+0.3535533905932737*p_parallel[2]*bxbz[7]-0.3535533905932737*bxbz[4]*p_perp[6]+0.3535533905932737*bxbz[4]*p_parallel[6]-0.3535533905932737*p_perp[4]*bxbz[6]+0.3535533905932737*p_parallel[4]*bxbz[6]-0.3535533905932737*bxbz[0]*p_perp[5]+0.3535533905932737*bxbz[0]*p_parallel[5]-0.3535533905932737*p_perp[0]*bxbz[5]+0.3535533905932737*p_parallel[0]*bxbz[5]-0.3535533905932737*bxbz[1]*p_perp[3]+0.3535533905932737*bxbz[1]*p_parallel[3]-0.3535533905932737*p_perp[1]*bxbz[3]+0.3535533905932737*p_parallel[1]*bxbz[3]; 
  Pxz[6] = (-0.3535533905932737*bxbz[1]*p_perp[7])+0.3535533905932737*bxbz[1]*p_parallel[7]-0.3535533905932737*p_perp[1]*bxbz[7]+0.3535533905932737*p_parallel[1]*bxbz[7]-0.3535533905932737*bxbz[0]*p_perp[6]+0.3535533905932737*bxbz[0]*p_parallel[6]-0.3535533905932737*p_perp[0]*bxbz[6]+0.3535533905932737*p_parallel[0]*bxbz[6]-0.3535533905932737*bxbz[4]*p_perp[5]+0.3535533905932737*bxbz[4]*p_parallel[5]-0.3535533905932737*p_perp[4]*bxbz[5]+0.3535533905932737*p_parallel[4]*bxbz[5]-0.3535533905932737*bxbz[2]*p_perp[3]+0.3535533905932737*bxbz[2]*p_parallel[3]-0.3535533905932737*p_perp[2]*bxbz[3]+0.3535533905932737*p_parallel[2]*bxbz[3]; 
  Pxz[7] = (-0.3535533905932737*bxbz[0]*p_perp[7])+0.3535533905932737*bxbz[0]*p_parallel[7]-0.3535533905932737*p_perp[0]*bxbz[7]+0.3535533905932737*p_parallel[0]*bxbz[7]-0.3535533905932737*bxbz[1]*p_perp[6]+0.3535533905932737*bxbz[1]*p_parallel[6]-0.3535533905932737*p_perp[1]*bxbz[6]+0.3535533905932737*p_parallel[1]*bxbz[6]-0.3535533905932737*bxbz[2]*p_perp[5]+0.3535533905932737*bxbz[2]*p_parallel[5]-0.3535533905932737*p_perp[2]*bxbz[5]+0.3535533905932737*p_parallel[2]*bxbz[5]-0.3535533905932737*bxbz[3]*p_perp[4]+0.3535533905932737*bxbz[3]*p_parallel[4]-0.3535533905932737*p_perp[3]*bxbz[4]+0.3535533905932737*p_parallel[3]*bxbz[4]; 

  Pyy[0] = (-0.3535533905932737*byby[7]*p_perp[7])+0.3535533905932737*byby[7]*p_parallel[7]-0.3535533905932737*byby[6]*p_perp[6]+0.3535533905932737*byby[6]*p_parallel[6]-0.3535533905932737*byby[5]*p_perp[5]+0.3535533905932737*byby[5]*p_parallel[5]-0.3535533905932737*byby[4]*p_perp[4]+0.3535533905932737*byby[4]*p_parallel[4]-0.3535533905932737*byby[3]*p_perp[3]+0.3535533905932737*byby[3]*p_parallel[3]-0.3535533905932737*byby[2]*p_perp[2]+0.3535533905932737*byby[2]*p_parallel[2]-0.3535533905932737*byby[1]*p_perp[1]+0.3535533905932737*byby[1]*p_parallel[1]-0.3535533905932737*byby[0]*p_perp[0]+p_perp[0]+0.3535533905932737*byby[0]*p_parallel[0]; 
  Pyy[1] = (-0.3535533905932737*byby[6]*p_perp[7])+0.3535533905932737*byby[6]*p_parallel[7]-0.3535533905932737*p_perp[6]*byby[7]+0.3535533905932737*p_parallel[6]*byby[7]-0.3535533905932737*byby[3]*p_perp[5]+0.3535533905932737*byby[3]*p_parallel[5]-0.3535533905932737*p_perp[3]*byby[5]+0.3535533905932737*p_parallel[3]*byby[5]-0.3535533905932737*byby[2]*p_perp[4]+0.3535533905932737*byby[2]*p_parallel[4]-0.3535533905932737*p_perp[2]*byby[4]+0.3535533905932737*p_parallel[2]*byby[4]-0.3535533905932737*byby[0]*p_perp[1]+p_perp[1]+0.3535533905932737*byby[0]*p_parallel[1]-0.3535533905932737*p_perp[0]*byby[1]+0.3535533905932737*p_parallel[0]*byby[1]; 
  Pyy[2] = (-0.3535533905932737*byby[5]*p_perp[7])+0.3535533905932737*byby[5]*p_parallel[7]-0.3535533905932737*p_perp[5]*byby[7]+0.3535533905932737*p_parallel[5]*byby[7]-0.3535533905932737*byby[3]*p_perp[6]+0.3535533905932737*byby[3]*p_parallel[6]-0.3535533905932737*p_perp[3]*byby[6]+0.3535533905932737*p_parallel[3]*byby[6]-0.3535533905932737*byby[1]*p_perp[4]+0.3535533905932737*byby[1]*p_parallel[4]-0.3535533905932737*p_perp[1]*byby[4]+0.3535533905932737*p_parallel[1]*byby[4]-0.3535533905932737*byby[0]*p_perp[2]+p_perp[2]+0.3535533905932737*byby[0]*p_parallel[2]-0.3535533905932737*p_perp[0]*byby[2]+0.3535533905932737*p_parallel[0]*byby[2]; 
  Pyy[3] = (-0.3535533905932737*byby[4]*p_perp[7])+0.3535533905932737*byby[4]*p_parallel[7]-0.3535533905932737*p_perp[4]*byby[7]+0.3535533905932737*p_parallel[4]*byby[7]-0.3535533905932737*byby[2]*p_perp[6]+0.3535533905932737*byby[2]*p_parallel[6]-0.3535533905932737*p_perp[2]*byby[6]+0.3535533905932737*p_parallel[2]*byby[6]-0.3535533905932737*byby[1]*p_perp[5]+0.3535533905932737*byby[1]*p_parallel[5]-0.3535533905932737*p_perp[1]*byby[5]+0.3535533905932737*p_parallel[1]*byby[5]-0.3535533905932737*byby[0]*p_perp[3]+p_perp[3]+0.3535533905932737*byby[0]*p_parallel[3]-0.3535533905932737*p_perp[0]*byby[3]+0.3535533905932737*p_parallel[0]*byby[3]; 
  Pyy[4] = (-0.3535533905932737*byby[3]*p_perp[7])+0.3535533905932737*byby[3]*p_parallel[7]-0.3535533905932737*p_perp[3]*byby[7]+0.3535533905932737*p_parallel[3]*byby[7]-0.3535533905932737*byby[5]*p_perp[6]+0.3535533905932737*byby[5]*p_parallel[6]-0.3535533905932737*p_perp[5]*byby[6]+0.3535533905932737*p_parallel[5]*byby[6]-0.3535533905932737*byby[0]*p_perp[4]+p_perp[4]+0.3535533905932737*byby[0]*p_parallel[4]-0.3535533905932737*p_perp[0]*byby[4]+0.3535533905932737*p_parallel[0]*byby[4]-0.3535533905932737*byby[1]*p_perp[2]+0.3535533905932737*byby[1]*p_parallel[2]-0.3535533905932737*p_perp[1]*byby[2]+0.3535533905932737*p_parallel[1]*byby[2]; 
  Pyy[5] = (-0.3535533905932737*byby[2]*p_perp[7])+0.3535533905932737*byby[2]*p_parallel[7]-0.3535533905932737*p_perp[2]*byby[7]+0.3535533905932737*p_parallel[2]*byby[7]-0.3535533905932737*byby[4]*p_perp[6]+0.3535533905932737*byby[4]*p_parallel[6]-0.3535533905932737*p_perp[4]*byby[6]+0.3535533905932737*p_parallel[4]*byby[6]-0.3535533905932737*byby[0]*p_perp[5]+p_perp[5]+0.3535533905932737*byby[0]*p_parallel[5]-0.3535533905932737*p_perp[0]*byby[5]+0.3535533905932737*p_parallel[0]*byby[5]-0.3535533905932737*byby[1]*p_perp[3]+0.3535533905932737*byby[1]*p_parallel[3]-0.3535533905932737*p_perp[1]*byby[3]+0.3535533905932737*p_parallel[1]*byby[3]; 
  Pyy[6] = (-0.3535533905932737*byby[1]*p_perp[7])+0.3535533905932737*byby[1]*p_parallel[7]-0.3535533905932737*p_perp[1]*byby[7]+0.3535533905932737*p_parallel[1]*byby[7]-0.3535533905932737*byby[0]*p_perp[6]+p_perp[6]+0.3535533905932737*byby[0]*p_parallel[6]-0.3535533905932737*p_perp[0]*byby[6]+0.3535533905932737*p_parallel[0]*byby[6]-0.3535533905932737*byby[4]*p_perp[5]+0.3535533905932737*byby[4]*p_parallel[5]-0.3535533905932737*p_perp[4]*byby[5]+0.3535533905932737*p_parallel[4]*byby[5]-0.3535533905932737*byby[2]*p_perp[3]+0.3535533905932737*byby[2]*p_parallel[3]-0.3535533905932737*p_perp[2]*byby[3]+0.3535533905932737*p_parallel[2]*byby[3]; 
  Pyy[7] = (-0.3535533905932737*byby[0]*p_perp[7])+p_perp[7]+0.3535533905932737*byby[0]*p_parallel[7]-0.3535533905932737*p_perp[0]*byby[7]+0.3535533905932737*p_parallel[0]*byby[7]-0.3535533905932737*byby[1]*p_perp[6]+0.3535533905932737*byby[1]*p_parallel[6]-0.3535533905932737*p_perp[1]*byby[6]+0.3535533905932737*p_parallel[1]*byby[6]-0.3535533905932737*byby[2]*p_perp[5]+0.3535533905932737*byby[2]*p_parallel[5]-0.3535533905932737*p_perp[2]*byby[5]+0.3535533905932737*p_parallel[2]*byby[5]-0.3535533905932737*byby[3]*p_perp[4]+0.3535533905932737*byby[3]*p_parallel[4]-0.3535533905932737*p_perp[3]*byby[4]+0.3535533905932737*p_parallel[3]*byby[4]; 

  Pyz[0] = (-0.3535533905932737*bybz[7]*p_perp[7])+0.3535533905932737*bybz[7]*p_parallel[7]-0.3535533905932737*bybz[6]*p_perp[6]+0.3535533905932737*bybz[6]*p_parallel[6]-0.3535533905932737*bybz[5]*p_perp[5]+0.3535533905932737*bybz[5]*p_parallel[5]-0.3535533905932737*bybz[4]*p_perp[4]+0.3535533905932737*bybz[4]*p_parallel[4]-0.3535533905932737*bybz[3]*p_perp[3]+0.3535533905932737*bybz[3]*p_parallel[3]-0.3535533905932737*bybz[2]*p_perp[2]+0.3535533905932737*bybz[2]*p_parallel[2]-0.3535533905932737*bybz[1]*p_perp[1]+0.3535533905932737*bybz[1]*p_parallel[1]-0.3535533905932737*bybz[0]*p_perp[0]+0.3535533905932737*bybz[0]*p_parallel[0]; 
  Pyz[1] = (-0.3535533905932737*bybz[6]*p_perp[7])+0.3535533905932737*bybz[6]*p_parallel[7]-0.3535533905932737*p_perp[6]*bybz[7]+0.3535533905932737*p_parallel[6]*bybz[7]-0.3535533905932737*bybz[3]*p_perp[5]+0.3535533905932737*bybz[3]*p_parallel[5]-0.3535533905932737*p_perp[3]*bybz[5]+0.3535533905932737*p_parallel[3]*bybz[5]-0.3535533905932737*bybz[2]*p_perp[4]+0.3535533905932737*bybz[2]*p_parallel[4]-0.3535533905932737*p_perp[2]*bybz[4]+0.3535533905932737*p_parallel[2]*bybz[4]-0.3535533905932737*bybz[0]*p_perp[1]+0.3535533905932737*bybz[0]*p_parallel[1]-0.3535533905932737*p_perp[0]*bybz[1]+0.3535533905932737*p_parallel[0]*bybz[1]; 
  Pyz[2] = (-0.3535533905932737*bybz[5]*p_perp[7])+0.3535533905932737*bybz[5]*p_parallel[7]-0.3535533905932737*p_perp[5]*bybz[7]+0.3535533905932737*p_parallel[5]*bybz[7]-0.3535533905932737*bybz[3]*p_perp[6]+0.3535533905932737*bybz[3]*p_parallel[6]-0.3535533905932737*p_perp[3]*bybz[6]+0.3535533905932737*p_parallel[3]*bybz[6]-0.3535533905932737*bybz[1]*p_perp[4]+0.3535533905932737*bybz[1]*p_parallel[4]-0.3535533905932737*p_perp[1]*bybz[4]+0.3535533905932737*p_parallel[1]*bybz[4]-0.3535533905932737*bybz[0]*p_perp[2]+0.3535533905932737*bybz[0]*p_parallel[2]-0.3535533905932737*p_perp[0]*bybz[2]+0.3535533905932737*p_parallel[0]*bybz[2]; 
  Pyz[3] = (-0.3535533905932737*bybz[4]*p_perp[7])+0.3535533905932737*bybz[4]*p_parallel[7]-0.3535533905932737*p_perp[4]*bybz[7]+0.3535533905932737*p_parallel[4]*bybz[7]-0.3535533905932737*bybz[2]*p_perp[6]+0.3535533905932737*bybz[2]*p_parallel[6]-0.3535533905932737*p_perp[2]*bybz[6]+0.3535533905932737*p_parallel[2]*bybz[6]-0.3535533905932737*bybz[1]*p_perp[5]+0.3535533905932737*bybz[1]*p_parallel[5]-0.3535533905932737*p_perp[1]*bybz[5]+0.3535533905932737*p_parallel[1]*bybz[5]-0.3535533905932737*bybz[0]*p_perp[3]+0.3535533905932737*bybz[0]*p_parallel[3]-0.3535533905932737*p_perp[0]*bybz[3]+0.3535533905932737*p_parallel[0]*bybz[3]; 
  Pyz[4] = (-0.3535533905932737*bybz[3]*p_perp[7])+0.3535533905932737*bybz[3]*p_parallel[7]-0.3535533905932737*p_perp[3]*bybz[7]+0.3535533905932737*p_parallel[3]*bybz[7]-0.3535533905932737*bybz[5]*p_perp[6]+0.3535533905932737*bybz[5]*p_parallel[6]-0.3535533905932737*p_perp[5]*bybz[6]+0.3535533905932737*p_parallel[5]*bybz[6]-0.3535533905932737*bybz[0]*p_perp[4]+0.3535533905932737*bybz[0]*p_parallel[4]-0.3535533905932737*p_perp[0]*bybz[4]+0.3535533905932737*p_parallel[0]*bybz[4]-0.3535533905932737*bybz[1]*p_perp[2]+0.3535533905932737*bybz[1]*p_parallel[2]-0.3535533905932737*p_perp[1]*bybz[2]+0.3535533905932737*p_parallel[1]*bybz[2]; 
  Pyz[5] = (-0.3535533905932737*bybz[2]*p_perp[7])+0.3535533905932737*bybz[2]*p_parallel[7]-0.3535533905932737*p_perp[2]*bybz[7]+0.3535533905932737*p_parallel[2]*bybz[7]-0.3535533905932737*bybz[4]*p_perp[6]+0.3535533905932737*bybz[4]*p_parallel[6]-0.3535533905932737*p_perp[4]*bybz[6]+0.3535533905932737*p_parallel[4]*bybz[6]-0.3535533905932737*bybz[0]*p_perp[5]+0.3535533905932737*bybz[0]*p_parallel[5]-0.3535533905932737*p_perp[0]*bybz[5]+0.3535533905932737*p_parallel[0]*bybz[5]-0.3535533905932737*bybz[1]*p_perp[3]+0.3535533905932737*bybz[1]*p_parallel[3]-0.3535533905932737*p_perp[1]*bybz[3]+0.3535533905932737*p_parallel[1]*bybz[3]; 
  Pyz[6] = (-0.3535533905932737*bybz[1]*p_perp[7])+0.3535533905932737*bybz[1]*p_parallel[7]-0.3535533905932737*p_perp[1]*bybz[7]+0.3535533905932737*p_parallel[1]*bybz[7]-0.3535533905932737*bybz[0]*p_perp[6]+0.3535533905932737*bybz[0]*p_parallel[6]-0.3535533905932737*p_perp[0]*bybz[6]+0.3535533905932737*p_parallel[0]*bybz[6]-0.3535533905932737*bybz[4]*p_perp[5]+0.3535533905932737*bybz[4]*p_parallel[5]-0.3535533905932737*p_perp[4]*bybz[5]+0.3535533905932737*p_parallel[4]*bybz[5]-0.3535533905932737*bybz[2]*p_perp[3]+0.3535533905932737*bybz[2]*p_parallel[3]-0.3535533905932737*p_perp[2]*bybz[3]+0.3535533905932737*p_parallel[2]*bybz[3]; 
  Pyz[7] = (-0.3535533905932737*bybz[0]*p_perp[7])+0.3535533905932737*bybz[0]*p_parallel[7]-0.3535533905932737*p_perp[0]*bybz[7]+0.3535533905932737*p_parallel[0]*bybz[7]-0.3535533905932737*bybz[1]*p_perp[6]+0.3535533905932737*bybz[1]*p_parallel[6]-0.3535533905932737*p_perp[1]*bybz[6]+0.3535533905932737*p_parallel[1]*bybz[6]-0.3535533905932737*bybz[2]*p_perp[5]+0.3535533905932737*bybz[2]*p_parallel[5]-0.3535533905932737*p_perp[2]*bybz[5]+0.3535533905932737*p_parallel[2]*bybz[5]-0.3535533905932737*bybz[3]*p_perp[4]+0.3535533905932737*bybz[3]*p_parallel[4]-0.3535533905932737*p_perp[3]*bybz[4]+0.3535533905932737*p_parallel[3]*bybz[4]; 

  Pzz[0] = (-0.3535533905932737*bzbz[7]*p_perp[7])+0.3535533905932737*bzbz[7]*p_parallel[7]-0.3535533905932737*bzbz[6]*p_perp[6]+0.3535533905932737*bzbz[6]*p_parallel[6]-0.3535533905932737*bzbz[5]*p_perp[5]+0.3535533905932737*bzbz[5]*p_parallel[5]-0.3535533905932737*bzbz[4]*p_perp[4]+0.3535533905932737*bzbz[4]*p_parallel[4]-0.3535533905932737*bzbz[3]*p_perp[3]+0.3535533905932737*bzbz[3]*p_parallel[3]-0.3535533905932737*bzbz[2]*p_perp[2]+0.3535533905932737*bzbz[2]*p_parallel[2]-0.3535533905932737*bzbz[1]*p_perp[1]+0.3535533905932737*bzbz[1]*p_parallel[1]-0.3535533905932737*bzbz[0]*p_perp[0]+p_perp[0]+0.3535533905932737*bzbz[0]*p_parallel[0]; 
  Pzz[1] = (-0.3535533905932737*bzbz[6]*p_perp[7])+0.3535533905932737*bzbz[6]*p_parallel[7]-0.3535533905932737*p_perp[6]*bzbz[7]+0.3535533905932737*p_parallel[6]*bzbz[7]-0.3535533905932737*bzbz[3]*p_perp[5]+0.3535533905932737*bzbz[3]*p_parallel[5]-0.3535533905932737*p_perp[3]*bzbz[5]+0.3535533905932737*p_parallel[3]*bzbz[5]-0.3535533905932737*bzbz[2]*p_perp[4]+0.3535533905932737*bzbz[2]*p_parallel[4]-0.3535533905932737*p_perp[2]*bzbz[4]+0.3535533905932737*p_parallel[2]*bzbz[4]-0.3535533905932737*bzbz[0]*p_perp[1]+p_perp[1]+0.3535533905932737*bzbz[0]*p_parallel[1]-0.3535533905932737*p_perp[0]*bzbz[1]+0.3535533905932737*p_parallel[0]*bzbz[1]; 
  Pzz[2] = (-0.3535533905932737*bzbz[5]*p_perp[7])+0.3535533905932737*bzbz[5]*p_parallel[7]-0.3535533905932737*p_perp[5]*bzbz[7]+0.3535533905932737*p_parallel[5]*bzbz[7]-0.3535533905932737*bzbz[3]*p_perp[6]+0.3535533905932737*bzbz[3]*p_parallel[6]-0.3535533905932737*p_perp[3]*bzbz[6]+0.3535533905932737*p_parallel[3]*bzbz[6]-0.3535533905932737*bzbz[1]*p_perp[4]+0.3535533905932737*bzbz[1]*p_parallel[4]-0.3535533905932737*p_perp[1]*bzbz[4]+0.3535533905932737*p_parallel[1]*bzbz[4]-0.3535533905932737*bzbz[0]*p_perp[2]+p_perp[2]+0.3535533905932737*bzbz[0]*p_parallel[2]-0.3535533905932737*p_perp[0]*bzbz[2]+0.3535533905932737*p_parallel[0]*bzbz[2]; 
  Pzz[3] = (-0.3535533905932737*bzbz[4]*p_perp[7])+0.3535533905932737*bzbz[4]*p_parallel[7]-0.3535533905932737*p_perp[4]*bzbz[7]+0.3535533905932737*p_parallel[4]*bzbz[7]-0.3535533905932737*bzbz[2]*p_perp[6]+0.3535533905932737*bzbz[2]*p_parallel[6]-0.3535533905932737*p_perp[2]*bzbz[6]+0.3535533905932737*p_parallel[2]*bzbz[6]-0.3535533905932737*bzbz[1]*p_perp[5]+0.3535533905932737*bzbz[1]*p_parallel[5]-0.3535533905932737*p_perp[1]*bzbz[5]+0.3535533905932737*p_parallel[1]*bzbz[5]-0.3535533905932737*bzbz[0]*p_perp[3]+p_perp[3]+0.3535533905932737*bzbz[0]*p_parallel[3]-0.3535533905932737*p_perp[0]*bzbz[3]+0.3535533905932737*p_parallel[0]*bzbz[3]; 
  Pzz[4] = (-0.3535533905932737*bzbz[3]*p_perp[7])+0.3535533905932737*bzbz[3]*p_parallel[7]-0.3535533905932737*p_perp[3]*bzbz[7]+0.3535533905932737*p_parallel[3]*bzbz[7]-0.3535533905932737*bzbz[5]*p_perp[6]+0.3535533905932737*bzbz[5]*p_parallel[6]-0.3535533905932737*p_perp[5]*bzbz[6]+0.3535533905932737*p_parallel[5]*bzbz[6]-0.3535533905932737*bzbz[0]*p_perp[4]+p_perp[4]+0.3535533905932737*bzbz[0]*p_parallel[4]-0.3535533905932737*p_perp[0]*bzbz[4]+0.3535533905932737*p_parallel[0]*bzbz[4]-0.3535533905932737*bzbz[1]*p_perp[2]+0.3535533905932737*bzbz[1]*p_parallel[2]-0.3535533905932737*p_perp[1]*bzbz[2]+0.3535533905932737*p_parallel[1]*bzbz[2]; 
  Pzz[5] = (-0.3535533905932737*bzbz[2]*p_perp[7])+0.3535533905932737*bzbz[2]*p_parallel[7]-0.3535533905932737*p_perp[2]*bzbz[7]+0.3535533905932737*p_parallel[2]*bzbz[7]-0.3535533905932737*bzbz[4]*p_perp[6]+0.3535533905932737*bzbz[4]*p_parallel[6]-0.3535533905932737*p_perp[4]*bzbz[6]+0.3535533905932737*p_parallel[4]*bzbz[6]-0.3535533905932737*bzbz[0]*p_perp[5]+p_perp[5]+0.3535533905932737*bzbz[0]*p_parallel[5]-0.3535533905932737*p_perp[0]*bzbz[5]+0.3535533905932737*p_parallel[0]*bzbz[5]-0.3535533905932737*bzbz[1]*p_perp[3]+0.3535533905932737*bzbz[1]*p_parallel[3]-0.3535533905932737*p_perp[1]*bzbz[3]+0.3535533905932737*p_parallel[1]*bzbz[3]; 
  Pzz[6] = (-0.3535533905932737*bzbz[1]*p_perp[7])+0.3535533905932737*bzbz[1]*p_parallel[7]-0.3535533905932737*p_perp[1]*bzbz[7]+0.3535533905932737*p_parallel[1]*bzbz[7]-0.3535533905932737*bzbz[0]*p_perp[6]+p_perp[6]+0.3535533905932737*bzbz[0]*p_parallel[6]-0.3535533905932737*p_perp[0]*bzbz[6]+0.3535533905932737*p_parallel[0]*bzbz[6]-0.3535533905932737*bzbz[4]*p_perp[5]+0.3535533905932737*bzbz[4]*p_parallel[5]-0.3535533905932737*p_perp[4]*bzbz[5]+0.3535533905932737*p_parallel[4]*bzbz[5]-0.3535533905932737*bzbz[2]*p_perp[3]+0.3535533905932737*bzbz[2]*p_parallel[3]-0.3535533905932737*p_perp[2]*bzbz[3]+0.3535533905932737*p_parallel[2]*bzbz[3]; 
  Pzz[7] = (-0.3535533905932737*bzbz[0]*p_perp[7])+p_perp[7]+0.3535533905932737*bzbz[0]*p_parallel[7]-0.3535533905932737*p_perp[0]*bzbz[7]+0.3535533905932737*p_parallel[0]*bzbz[7]-0.3535533905932737*bzbz[1]*p_perp[6]+0.3535533905932737*bzbz[1]*p_parallel[6]-0.3535533905932737*p_perp[1]*bzbz[6]+0.3535533905932737*p_parallel[1]*bzbz[6]-0.3535533905932737*bzbz[2]*p_perp[5]+0.3535533905932737*bzbz[2]*p_parallel[5]-0.3535533905932737*p_perp[2]*bzbz[5]+0.3535533905932737*p_parallel[2]*bzbz[5]-0.3535533905932737*bzbz[3]*p_perp[4]+0.3535533905932737*bzbz[3]*p_parallel[4]-0.3535533905932737*p_perp[3]*bzbz[4]+0.3535533905932737*p_parallel[3]*bzbz[4]; 
} 
