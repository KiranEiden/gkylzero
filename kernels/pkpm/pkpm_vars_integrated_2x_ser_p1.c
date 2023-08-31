#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void pkpm_vars_integrated_2x_ser_p1(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* prim, double* GKYL_RESTRICT int_pkpm_vars) 
{ 
  // vlasov_pkpm_moms: [rho, p_parallel, p_perp], Moments computed from kinetic equation in pkpm model.
  // euler_pkpm:       [rho ux, rho uy, rho uz], Fluid input state vector.
  // prim:             [ux, uy, uz, 1/rho div(p_par b), T_perp/m, m/T_perp].
  // int_pkpm_vars:    Output integrated variables.

  const double *rho = &vlasov_pkpm_moms[0]; 
  const double *p_parallel = &vlasov_pkpm_moms[4]; 
  const double *p_perp = &vlasov_pkpm_moms[8]; 

  const double *rhoux = &euler_pkpm[0]; 
  const double *rhouy = &euler_pkpm[4]; 
  const double *rhouz = &euler_pkpm[8]; 

  const double *ux = &prim[0]; 
  const double *uy = &prim[4]; 
  const double *uz = &prim[8]; 

  // Order of integrated variables is (rho, rhoux, rhouy, rhouz, rho ux^2, rho uy^2, rho uz^2, p_parallel, p_perp) 
  int_pkpm_vars[0] += 2.0*rho[0]; 
  int_pkpm_vars[1] += 2.0*rhoux[0]; 
  int_pkpm_vars[2] += 2.0*rhouy[0]; 
  int_pkpm_vars[3] += 2.0*rhouz[0]; 
  int_pkpm_vars[4] += rhoux[3]*ux[3]+rhoux[2]*ux[2]+rhoux[1]*ux[1]+rhoux[0]*ux[0]; 
  int_pkpm_vars[5] += rhouy[3]*uy[3]+rhouy[2]*uy[2]+rhouy[1]*uy[1]+rhouy[0]*uy[0]; 
  int_pkpm_vars[6] += rhouz[3]*uz[3]+rhouz[2]*uz[2]+rhouz[1]*uz[1]+rhouz[0]*uz[0]; 
  int_pkpm_vars[7] += 2.0*p_parallel[0]; 
  int_pkpm_vars[8] += 2.0*p_perp[0]; 
} 
