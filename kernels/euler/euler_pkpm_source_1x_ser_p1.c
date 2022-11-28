#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void euler_pkpm_source_1x_ser_p1(const double *qmem, const double *nu, const double *nu_vth_sq, const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *p_perp_source, double* out) 
{ 
  // qmem:             q/m*EM fields.
  // nu:               Collisionality.
  // nu_vth_sq:        nu*vth^2 = nu*T/m.
  // vlasov_pkpm_moms: [rho, p_parallel, q_parallel], Moments computed from kinetic equation in pkpm model.
  // euler_pkpm:       Input fluid variables.
  // p_perp_source:    Input perpendicular pressure compression source (p_perp div(u) - p_perp bb : grad(u)).
  // out:              Output increment
  const double *rho = &vlasov_pkpm_moms[0]; 

  const double *Ex = &qmem[0]; 
  const double *Ey = &qmem[2]; 
  const double *Ez = &qmem[4]; 
  const double *Bx = &qmem[6]; 
  const double *By = &qmem[8]; 
  const double *Bz = &qmem[10]; 

  const double *rhoux = &euler_pkpm[0]; 
  const double *rhouy = &euler_pkpm[2]; 
  const double *rhouz = &euler_pkpm[4]; 
  const double *p_perp = &euler_pkpm[6]; 

  double *outrhoux = &out[0]; 
  double *outrhouy = &out[2]; 
  double *outrhouz = &out[4]; 
  double *outp_perp = &out[6]; 

  outrhoux[0] += (-0.7071067811865475*By[1]*rhouz[1])+0.7071067811865475*Bz[1]*rhouy[1]+0.7071067811865475*Ex[1]*rho[1]-0.7071067811865475*By[0]*rhouz[0]+0.7071067811865475*Bz[0]*rhouy[0]+0.7071067811865475*Ex[0]*rho[0]; 
  outrhoux[1] += (-0.7071067811865475*By[0]*rhouz[1])+0.7071067811865475*Bz[0]*rhouy[1]+0.7071067811865475*Ex[0]*rho[1]+0.7071067811865475*rho[0]*Ex[1]+0.7071067811865475*rhouy[0]*Bz[1]-0.7071067811865475*rhouz[0]*By[1]; 

  outrhouy[0] += 0.7071067811865475*Bx[1]*rhouz[1]-0.7071067811865475*Bz[1]*rhoux[1]+0.7071067811865475*Ey[1]*rho[1]+0.7071067811865475*Bx[0]*rhouz[0]-0.7071067811865475*Bz[0]*rhoux[0]+0.7071067811865475*Ey[0]*rho[0]; 
  outrhouy[1] += 0.7071067811865475*Bx[0]*rhouz[1]-0.7071067811865475*Bz[0]*rhoux[1]+0.7071067811865475*Ey[0]*rho[1]+0.7071067811865475*rho[0]*Ey[1]-0.7071067811865475*rhoux[0]*Bz[1]+0.7071067811865475*rhouz[0]*Bx[1]; 

  outrhouz[0] += (-0.7071067811865475*Bx[1]*rhouy[1])+0.7071067811865475*By[1]*rhoux[1]+0.7071067811865475*Ez[1]*rho[1]-0.7071067811865475*Bx[0]*rhouy[0]+0.7071067811865475*By[0]*rhoux[0]+0.7071067811865475*Ez[0]*rho[0]; 
  outrhouz[1] += (-0.7071067811865475*Bx[0]*rhouy[1])+0.7071067811865475*By[0]*rhoux[1]+0.7071067811865475*Ez[0]*rho[1]+0.7071067811865475*rho[0]*Ez[1]+0.7071067811865475*rhoux[0]*By[1]-0.7071067811865475*rhouy[0]*Bx[1]; 

  outp_perp[0] += 0.7071067811865475*nu_vth_sq[1]*rho[1]-0.7071067811865475*nu[1]*p_perp[1]+0.7071067811865475*nu_vth_sq[0]*rho[0]+p_perp_source[0]-0.7071067811865475*nu[0]*p_perp[0]; 
  outp_perp[1] += 0.7071067811865475*nu_vth_sq[0]*rho[1]+p_perp_source[1]-0.7071067811865475*nu[0]*p_perp[1]+0.7071067811865475*rho[0]*nu_vth_sq[1]-0.7071067811865475*p_perp[0]*nu[1]; 

} 
