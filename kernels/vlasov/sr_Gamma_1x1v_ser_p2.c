#include <gkyl_sr_Gamma_kernels.h> 
#include <gkyl_basis_ser_1x_p2_exp_sq.h> 
#include <gkyl_basis_ser_1x_p2_inv.h> 
#include <gkyl_basis_ser_1x_p2_sqrt.h> 
GKYL_CU_DH void sr_Gamma_1x1v_ser_p2(const double *V, double* GKYL_RESTRICT Gamma) 
{ 
  // V:     Input velocity. 
  // Gamma: Gamma = 1/sqrt(1 - V^2/c^2). 
 
  const double *V_0 = &V[0]; 
  double V_0_sq[3] = {0.0}; 
  ser_1x_p2_exp_sq(V_0, V_0_sq); 
 
  double Gamma2_inv[3] = {0.0}; 
  double Gamma2[3] = {0.0}; 
 
  Gamma2_inv[0] = 1.414213562373095-1.0*V_0_sq[0]; 
  Gamma2_inv[1] = -1.0*V_0_sq[1]; 
  Gamma2_inv[2] = -1.0*V_0_sq[2]; 

  bool notCellAvg = true;
  if (notCellAvg && (1.58113883008419*Gamma2_inv[2]-1.224744871391589*Gamma2_inv[1]+0.7071067811865475*Gamma2_inv[0] < 0)) notCellAvg = false; 
  if (notCellAvg && (1.58113883008419*Gamma2_inv[2]+1.224744871391589*Gamma2_inv[1]+0.7071067811865475*Gamma2_inv[0] < 0)) notCellAvg = false; 
 
  if (notCellAvg) { 
  ser_1x_p2_inv(Gamma2_inv, Gamma2); 
  } else { 
  Gamma2[0] = 2.0/Gamma2_inv[0]; 
  Gamma2[1] = 0.0; 
  Gamma2[2] = 0.0; 
  } 
  ser_1x_p2_sqrt(Gamma2, Gamma); 
} 
 
