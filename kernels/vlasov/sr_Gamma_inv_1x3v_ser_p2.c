#include <gkyl_sr_Gamma_kernels.h> 
#include <gkyl_basis_ser_1x_p2_exp_sq.h> 
#include <gkyl_basis_ser_1x_p2_sqrt.h> 
GKYL_CU_DH void sr_Gamma_inv_1x3v_ser_p2(const double *V, double* GKYL_RESTRICT Gamma_inv) 
{ 
  // V:     Input velocity. 
  // Gamma: Gamma = 1/sqrt(1 - V^2/c^2). 
 
  const double *V_0 = &V[0]; 
  double V_0_sq[3] = {0.0}; 
  ser_1x_p2_exp_sq(V_0, V_0_sq); 
 
  const double *V_1 = &V[3]; 
  double V_1_sq[3] = {0.0}; 
  ser_1x_p2_exp_sq(V_1, V_1_sq); 
 
  const double *V_2 = &V[6]; 
  double V_2_sq[3] = {0.0}; 
  ser_1x_p2_exp_sq(V_2, V_2_sq); 
 
  double Gamma2_inv[3] = {0.0}; 
 
  Gamma2_inv[0] = (-1.0*V_2_sq[0])-1.0*V_1_sq[0]-1.0*V_0_sq[0]+1.414213562373095; 
  Gamma2_inv[1] = (-1.0*V_2_sq[1])-1.0*V_1_sq[1]-1.0*V_0_sq[1]; 
  Gamma2_inv[2] = (-1.0*V_2_sq[2])-1.0*V_1_sq[2]-1.0*V_0_sq[2]; 

  ser_1x_p2_sqrt(Gamma2_inv, Gamma_inv); 
} 
 
