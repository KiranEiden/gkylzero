#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void fluid_vars_source_3x_ser_p1(const double* app_accel, const double* fluid, double* GKYL_RESTRICT out) 
{ 
  // app_accel: External applied acceleration (external forces).
  // fluid:     [rho, rho ux, rho uy, rho uz...], Fluid input state vector.
  // out:       Output increment

  const double *app_accel_x = &app_accel[0]; 
  const double *app_accel_y = &app_accel[8]; 
  const double *app_accel_z = &app_accel[16]; 

  const double *rho = &fluid[0]; 
  const double *rhoux = &fluid[8]; 
  const double *rhouy = &fluid[16]; 
  const double *rhouz = &fluid[24]; 

  double *out_rhoux = &out[8]; 
  double *out_rhouy = &out[16]; 
  double *out_rhouz = &out[24]; 
  double *out_energy = &out[32]; 

  out_rhoux[0] += 0.3535533905932737*app_accel_x[7]*rho[7]+0.3535533905932737*app_accel_x[6]*rho[6]+0.3535533905932737*app_accel_x[5]*rho[5]+0.3535533905932737*app_accel_x[4]*rho[4]+0.3535533905932737*app_accel_x[3]*rho[3]+0.3535533905932737*app_accel_x[2]*rho[2]+0.3535533905932737*app_accel_x[1]*rho[1]+0.3535533905932737*app_accel_x[0]*rho[0]; 
  out_rhoux[1] += 0.3535533905932737*app_accel_x[6]*rho[7]+0.3535533905932737*rho[6]*app_accel_x[7]+0.3535533905932737*app_accel_x[3]*rho[5]+0.3535533905932737*rho[3]*app_accel_x[5]+0.3535533905932737*app_accel_x[2]*rho[4]+0.3535533905932737*rho[2]*app_accel_x[4]+0.3535533905932737*app_accel_x[0]*rho[1]+0.3535533905932737*rho[0]*app_accel_x[1]; 
  out_rhoux[2] += 0.3535533905932737*app_accel_x[5]*rho[7]+0.3535533905932737*rho[5]*app_accel_x[7]+0.3535533905932737*app_accel_x[3]*rho[6]+0.3535533905932737*rho[3]*app_accel_x[6]+0.3535533905932737*app_accel_x[1]*rho[4]+0.3535533905932737*rho[1]*app_accel_x[4]+0.3535533905932737*app_accel_x[0]*rho[2]+0.3535533905932737*rho[0]*app_accel_x[2]; 
  out_rhoux[3] += 0.3535533905932737*app_accel_x[4]*rho[7]+0.3535533905932737*rho[4]*app_accel_x[7]+0.3535533905932737*app_accel_x[2]*rho[6]+0.3535533905932737*rho[2]*app_accel_x[6]+0.3535533905932737*app_accel_x[1]*rho[5]+0.3535533905932737*rho[1]*app_accel_x[5]+0.3535533905932737*app_accel_x[0]*rho[3]+0.3535533905932737*rho[0]*app_accel_x[3]; 
  out_rhoux[4] += 0.3535533905932737*app_accel_x[3]*rho[7]+0.3535533905932737*rho[3]*app_accel_x[7]+0.3535533905932737*app_accel_x[5]*rho[6]+0.3535533905932737*rho[5]*app_accel_x[6]+0.3535533905932737*app_accel_x[0]*rho[4]+0.3535533905932737*rho[0]*app_accel_x[4]+0.3535533905932737*app_accel_x[1]*rho[2]+0.3535533905932737*rho[1]*app_accel_x[2]; 
  out_rhoux[5] += 0.3535533905932737*app_accel_x[2]*rho[7]+0.3535533905932737*rho[2]*app_accel_x[7]+0.3535533905932737*app_accel_x[4]*rho[6]+0.3535533905932737*rho[4]*app_accel_x[6]+0.3535533905932737*app_accel_x[0]*rho[5]+0.3535533905932737*rho[0]*app_accel_x[5]+0.3535533905932737*app_accel_x[1]*rho[3]+0.3535533905932737*rho[1]*app_accel_x[3]; 
  out_rhoux[6] += 0.3535533905932737*app_accel_x[1]*rho[7]+0.3535533905932737*rho[1]*app_accel_x[7]+0.3535533905932737*app_accel_x[0]*rho[6]+0.3535533905932737*rho[0]*app_accel_x[6]+0.3535533905932737*app_accel_x[4]*rho[5]+0.3535533905932737*rho[4]*app_accel_x[5]+0.3535533905932737*app_accel_x[2]*rho[3]+0.3535533905932737*rho[2]*app_accel_x[3]; 
  out_rhoux[7] += 0.3535533905932737*app_accel_x[0]*rho[7]+0.3535533905932737*rho[0]*app_accel_x[7]+0.3535533905932737*app_accel_x[1]*rho[6]+0.3535533905932737*rho[1]*app_accel_x[6]+0.3535533905932737*app_accel_x[2]*rho[5]+0.3535533905932737*rho[2]*app_accel_x[5]+0.3535533905932737*app_accel_x[3]*rho[4]+0.3535533905932737*rho[3]*app_accel_x[4]; 

  out_rhouy[0] += 0.3535533905932737*app_accel_y[7]*rho[7]+0.3535533905932737*app_accel_y[6]*rho[6]+0.3535533905932737*app_accel_y[5]*rho[5]+0.3535533905932737*app_accel_y[4]*rho[4]+0.3535533905932737*app_accel_y[3]*rho[3]+0.3535533905932737*app_accel_y[2]*rho[2]+0.3535533905932737*app_accel_y[1]*rho[1]+0.3535533905932737*app_accel_y[0]*rho[0]; 
  out_rhouy[1] += 0.3535533905932737*app_accel_y[6]*rho[7]+0.3535533905932737*rho[6]*app_accel_y[7]+0.3535533905932737*app_accel_y[3]*rho[5]+0.3535533905932737*rho[3]*app_accel_y[5]+0.3535533905932737*app_accel_y[2]*rho[4]+0.3535533905932737*rho[2]*app_accel_y[4]+0.3535533905932737*app_accel_y[0]*rho[1]+0.3535533905932737*rho[0]*app_accel_y[1]; 
  out_rhouy[2] += 0.3535533905932737*app_accel_y[5]*rho[7]+0.3535533905932737*rho[5]*app_accel_y[7]+0.3535533905932737*app_accel_y[3]*rho[6]+0.3535533905932737*rho[3]*app_accel_y[6]+0.3535533905932737*app_accel_y[1]*rho[4]+0.3535533905932737*rho[1]*app_accel_y[4]+0.3535533905932737*app_accel_y[0]*rho[2]+0.3535533905932737*rho[0]*app_accel_y[2]; 
  out_rhouy[3] += 0.3535533905932737*app_accel_y[4]*rho[7]+0.3535533905932737*rho[4]*app_accel_y[7]+0.3535533905932737*app_accel_y[2]*rho[6]+0.3535533905932737*rho[2]*app_accel_y[6]+0.3535533905932737*app_accel_y[1]*rho[5]+0.3535533905932737*rho[1]*app_accel_y[5]+0.3535533905932737*app_accel_y[0]*rho[3]+0.3535533905932737*rho[0]*app_accel_y[3]; 
  out_rhouy[4] += 0.3535533905932737*app_accel_y[3]*rho[7]+0.3535533905932737*rho[3]*app_accel_y[7]+0.3535533905932737*app_accel_y[5]*rho[6]+0.3535533905932737*rho[5]*app_accel_y[6]+0.3535533905932737*app_accel_y[0]*rho[4]+0.3535533905932737*rho[0]*app_accel_y[4]+0.3535533905932737*app_accel_y[1]*rho[2]+0.3535533905932737*rho[1]*app_accel_y[2]; 
  out_rhouy[5] += 0.3535533905932737*app_accel_y[2]*rho[7]+0.3535533905932737*rho[2]*app_accel_y[7]+0.3535533905932737*app_accel_y[4]*rho[6]+0.3535533905932737*rho[4]*app_accel_y[6]+0.3535533905932737*app_accel_y[0]*rho[5]+0.3535533905932737*rho[0]*app_accel_y[5]+0.3535533905932737*app_accel_y[1]*rho[3]+0.3535533905932737*rho[1]*app_accel_y[3]; 
  out_rhouy[6] += 0.3535533905932737*app_accel_y[1]*rho[7]+0.3535533905932737*rho[1]*app_accel_y[7]+0.3535533905932737*app_accel_y[0]*rho[6]+0.3535533905932737*rho[0]*app_accel_y[6]+0.3535533905932737*app_accel_y[4]*rho[5]+0.3535533905932737*rho[4]*app_accel_y[5]+0.3535533905932737*app_accel_y[2]*rho[3]+0.3535533905932737*rho[2]*app_accel_y[3]; 
  out_rhouy[7] += 0.3535533905932737*app_accel_y[0]*rho[7]+0.3535533905932737*rho[0]*app_accel_y[7]+0.3535533905932737*app_accel_y[1]*rho[6]+0.3535533905932737*rho[1]*app_accel_y[6]+0.3535533905932737*app_accel_y[2]*rho[5]+0.3535533905932737*rho[2]*app_accel_y[5]+0.3535533905932737*app_accel_y[3]*rho[4]+0.3535533905932737*rho[3]*app_accel_y[4]; 

  out_rhouz[0] += 0.3535533905932737*app_accel_z[7]*rho[7]+0.3535533905932737*app_accel_z[6]*rho[6]+0.3535533905932737*app_accel_z[5]*rho[5]+0.3535533905932737*app_accel_z[4]*rho[4]+0.3535533905932737*app_accel_z[3]*rho[3]+0.3535533905932737*app_accel_z[2]*rho[2]+0.3535533905932737*app_accel_z[1]*rho[1]+0.3535533905932737*app_accel_z[0]*rho[0]; 
  out_rhouz[1] += 0.3535533905932737*app_accel_z[6]*rho[7]+0.3535533905932737*rho[6]*app_accel_z[7]+0.3535533905932737*app_accel_z[3]*rho[5]+0.3535533905932737*rho[3]*app_accel_z[5]+0.3535533905932737*app_accel_z[2]*rho[4]+0.3535533905932737*rho[2]*app_accel_z[4]+0.3535533905932737*app_accel_z[0]*rho[1]+0.3535533905932737*rho[0]*app_accel_z[1]; 
  out_rhouz[2] += 0.3535533905932737*app_accel_z[5]*rho[7]+0.3535533905932737*rho[5]*app_accel_z[7]+0.3535533905932737*app_accel_z[3]*rho[6]+0.3535533905932737*rho[3]*app_accel_z[6]+0.3535533905932737*app_accel_z[1]*rho[4]+0.3535533905932737*rho[1]*app_accel_z[4]+0.3535533905932737*app_accel_z[0]*rho[2]+0.3535533905932737*rho[0]*app_accel_z[2]; 
  out_rhouz[3] += 0.3535533905932737*app_accel_z[4]*rho[7]+0.3535533905932737*rho[4]*app_accel_z[7]+0.3535533905932737*app_accel_z[2]*rho[6]+0.3535533905932737*rho[2]*app_accel_z[6]+0.3535533905932737*app_accel_z[1]*rho[5]+0.3535533905932737*rho[1]*app_accel_z[5]+0.3535533905932737*app_accel_z[0]*rho[3]+0.3535533905932737*rho[0]*app_accel_z[3]; 
  out_rhouz[4] += 0.3535533905932737*app_accel_z[3]*rho[7]+0.3535533905932737*rho[3]*app_accel_z[7]+0.3535533905932737*app_accel_z[5]*rho[6]+0.3535533905932737*rho[5]*app_accel_z[6]+0.3535533905932737*app_accel_z[0]*rho[4]+0.3535533905932737*rho[0]*app_accel_z[4]+0.3535533905932737*app_accel_z[1]*rho[2]+0.3535533905932737*rho[1]*app_accel_z[2]; 
  out_rhouz[5] += 0.3535533905932737*app_accel_z[2]*rho[7]+0.3535533905932737*rho[2]*app_accel_z[7]+0.3535533905932737*app_accel_z[4]*rho[6]+0.3535533905932737*rho[4]*app_accel_z[6]+0.3535533905932737*app_accel_z[0]*rho[5]+0.3535533905932737*rho[0]*app_accel_z[5]+0.3535533905932737*app_accel_z[1]*rho[3]+0.3535533905932737*rho[1]*app_accel_z[3]; 
  out_rhouz[6] += 0.3535533905932737*app_accel_z[1]*rho[7]+0.3535533905932737*rho[1]*app_accel_z[7]+0.3535533905932737*app_accel_z[0]*rho[6]+0.3535533905932737*rho[0]*app_accel_z[6]+0.3535533905932737*app_accel_z[4]*rho[5]+0.3535533905932737*rho[4]*app_accel_z[5]+0.3535533905932737*app_accel_z[2]*rho[3]+0.3535533905932737*rho[2]*app_accel_z[3]; 
  out_rhouz[7] += 0.3535533905932737*app_accel_z[0]*rho[7]+0.3535533905932737*rho[0]*app_accel_z[7]+0.3535533905932737*app_accel_z[1]*rho[6]+0.3535533905932737*rho[1]*app_accel_z[6]+0.3535533905932737*app_accel_z[2]*rho[5]+0.3535533905932737*rho[2]*app_accel_z[5]+0.3535533905932737*app_accel_z[3]*rho[4]+0.3535533905932737*rho[3]*app_accel_z[4]; 

  out_energy[0] += 0.3535533905932737*app_accel_z[7]*rhouz[7]+0.3535533905932737*app_accel_y[7]*rhouy[7]+0.3535533905932737*app_accel_x[7]*rhoux[7]+0.3535533905932737*app_accel_z[6]*rhouz[6]+0.3535533905932737*app_accel_y[6]*rhouy[6]+0.3535533905932737*app_accel_x[6]*rhoux[6]+0.3535533905932737*app_accel_z[5]*rhouz[5]+0.3535533905932737*app_accel_y[5]*rhouy[5]+0.3535533905932737*app_accel_x[5]*rhoux[5]+0.3535533905932737*app_accel_z[4]*rhouz[4]+0.3535533905932737*app_accel_y[4]*rhouy[4]+0.3535533905932737*app_accel_x[4]*rhoux[4]+0.3535533905932737*app_accel_z[3]*rhouz[3]+0.3535533905932737*app_accel_y[3]*rhouy[3]+0.3535533905932737*app_accel_x[3]*rhoux[3]+0.3535533905932737*app_accel_z[2]*rhouz[2]+0.3535533905932737*app_accel_y[2]*rhouy[2]+0.3535533905932737*app_accel_x[2]*rhoux[2]+0.3535533905932737*app_accel_z[1]*rhouz[1]+0.3535533905932737*app_accel_y[1]*rhouy[1]+0.3535533905932737*app_accel_x[1]*rhoux[1]+0.3535533905932737*app_accel_z[0]*rhouz[0]+0.3535533905932737*app_accel_y[0]*rhouy[0]+0.3535533905932737*app_accel_x[0]*rhoux[0]; 
  out_energy[1] += 0.3535533905932737*app_accel_z[6]*rhouz[7]+0.3535533905932737*app_accel_y[6]*rhouy[7]+0.3535533905932737*app_accel_x[6]*rhoux[7]+0.3535533905932737*rhouz[6]*app_accel_z[7]+0.3535533905932737*rhouy[6]*app_accel_y[7]+0.3535533905932737*rhoux[6]*app_accel_x[7]+0.3535533905932737*app_accel_z[3]*rhouz[5]+0.3535533905932737*app_accel_y[3]*rhouy[5]+0.3535533905932737*app_accel_x[3]*rhoux[5]+0.3535533905932737*rhouz[3]*app_accel_z[5]+0.3535533905932737*rhouy[3]*app_accel_y[5]+0.3535533905932737*rhoux[3]*app_accel_x[5]+0.3535533905932737*app_accel_z[2]*rhouz[4]+0.3535533905932737*app_accel_y[2]*rhouy[4]+0.3535533905932737*app_accel_x[2]*rhoux[4]+0.3535533905932737*rhouz[2]*app_accel_z[4]+0.3535533905932737*rhouy[2]*app_accel_y[4]+0.3535533905932737*rhoux[2]*app_accel_x[4]+0.3535533905932737*app_accel_z[0]*rhouz[1]+0.3535533905932737*app_accel_y[0]*rhouy[1]+0.3535533905932737*app_accel_x[0]*rhoux[1]+0.3535533905932737*rhouz[0]*app_accel_z[1]+0.3535533905932737*rhouy[0]*app_accel_y[1]+0.3535533905932737*rhoux[0]*app_accel_x[1]; 
  out_energy[2] += 0.3535533905932737*app_accel_z[5]*rhouz[7]+0.3535533905932737*app_accel_y[5]*rhouy[7]+0.3535533905932737*app_accel_x[5]*rhoux[7]+0.3535533905932737*rhouz[5]*app_accel_z[7]+0.3535533905932737*rhouy[5]*app_accel_y[7]+0.3535533905932737*rhoux[5]*app_accel_x[7]+0.3535533905932737*app_accel_z[3]*rhouz[6]+0.3535533905932737*app_accel_y[3]*rhouy[6]+0.3535533905932737*app_accel_x[3]*rhoux[6]+0.3535533905932737*rhouz[3]*app_accel_z[6]+0.3535533905932737*rhouy[3]*app_accel_y[6]+0.3535533905932737*rhoux[3]*app_accel_x[6]+0.3535533905932737*app_accel_z[1]*rhouz[4]+0.3535533905932737*app_accel_y[1]*rhouy[4]+0.3535533905932737*app_accel_x[1]*rhoux[4]+0.3535533905932737*rhouz[1]*app_accel_z[4]+0.3535533905932737*rhouy[1]*app_accel_y[4]+0.3535533905932737*rhoux[1]*app_accel_x[4]+0.3535533905932737*app_accel_z[0]*rhouz[2]+0.3535533905932737*app_accel_y[0]*rhouy[2]+0.3535533905932737*app_accel_x[0]*rhoux[2]+0.3535533905932737*rhouz[0]*app_accel_z[2]+0.3535533905932737*rhouy[0]*app_accel_y[2]+0.3535533905932737*rhoux[0]*app_accel_x[2]; 
  out_energy[3] += 0.3535533905932737*app_accel_z[4]*rhouz[7]+0.3535533905932737*app_accel_y[4]*rhouy[7]+0.3535533905932737*app_accel_x[4]*rhoux[7]+0.3535533905932737*rhouz[4]*app_accel_z[7]+0.3535533905932737*rhouy[4]*app_accel_y[7]+0.3535533905932737*rhoux[4]*app_accel_x[7]+0.3535533905932737*app_accel_z[2]*rhouz[6]+0.3535533905932737*app_accel_y[2]*rhouy[6]+0.3535533905932737*app_accel_x[2]*rhoux[6]+0.3535533905932737*rhouz[2]*app_accel_z[6]+0.3535533905932737*rhouy[2]*app_accel_y[6]+0.3535533905932737*rhoux[2]*app_accel_x[6]+0.3535533905932737*app_accel_z[1]*rhouz[5]+0.3535533905932737*app_accel_y[1]*rhouy[5]+0.3535533905932737*app_accel_x[1]*rhoux[5]+0.3535533905932737*rhouz[1]*app_accel_z[5]+0.3535533905932737*rhouy[1]*app_accel_y[5]+0.3535533905932737*rhoux[1]*app_accel_x[5]+0.3535533905932737*app_accel_z[0]*rhouz[3]+0.3535533905932737*app_accel_y[0]*rhouy[3]+0.3535533905932737*app_accel_x[0]*rhoux[3]+0.3535533905932737*rhouz[0]*app_accel_z[3]+0.3535533905932737*rhouy[0]*app_accel_y[3]+0.3535533905932737*rhoux[0]*app_accel_x[3]; 
  out_energy[4] += 0.3535533905932737*app_accel_z[3]*rhouz[7]+0.3535533905932737*app_accel_y[3]*rhouy[7]+0.3535533905932737*app_accel_x[3]*rhoux[7]+0.3535533905932737*rhouz[3]*app_accel_z[7]+0.3535533905932737*rhouy[3]*app_accel_y[7]+0.3535533905932737*rhoux[3]*app_accel_x[7]+0.3535533905932737*app_accel_z[5]*rhouz[6]+0.3535533905932737*app_accel_y[5]*rhouy[6]+0.3535533905932737*app_accel_x[5]*rhoux[6]+0.3535533905932737*rhouz[5]*app_accel_z[6]+0.3535533905932737*rhouy[5]*app_accel_y[6]+0.3535533905932737*rhoux[5]*app_accel_x[6]+0.3535533905932737*app_accel_z[0]*rhouz[4]+0.3535533905932737*app_accel_y[0]*rhouy[4]+0.3535533905932737*app_accel_x[0]*rhoux[4]+0.3535533905932737*rhouz[0]*app_accel_z[4]+0.3535533905932737*rhouy[0]*app_accel_y[4]+0.3535533905932737*rhoux[0]*app_accel_x[4]+0.3535533905932737*app_accel_z[1]*rhouz[2]+0.3535533905932737*app_accel_y[1]*rhouy[2]+0.3535533905932737*app_accel_x[1]*rhoux[2]+0.3535533905932737*rhouz[1]*app_accel_z[2]+0.3535533905932737*rhouy[1]*app_accel_y[2]+0.3535533905932737*rhoux[1]*app_accel_x[2]; 
  out_energy[5] += 0.3535533905932737*app_accel_z[2]*rhouz[7]+0.3535533905932737*app_accel_y[2]*rhouy[7]+0.3535533905932737*app_accel_x[2]*rhoux[7]+0.3535533905932737*rhouz[2]*app_accel_z[7]+0.3535533905932737*rhouy[2]*app_accel_y[7]+0.3535533905932737*rhoux[2]*app_accel_x[7]+0.3535533905932737*app_accel_z[4]*rhouz[6]+0.3535533905932737*app_accel_y[4]*rhouy[6]+0.3535533905932737*app_accel_x[4]*rhoux[6]+0.3535533905932737*rhouz[4]*app_accel_z[6]+0.3535533905932737*rhouy[4]*app_accel_y[6]+0.3535533905932737*rhoux[4]*app_accel_x[6]+0.3535533905932737*app_accel_z[0]*rhouz[5]+0.3535533905932737*app_accel_y[0]*rhouy[5]+0.3535533905932737*app_accel_x[0]*rhoux[5]+0.3535533905932737*rhouz[0]*app_accel_z[5]+0.3535533905932737*rhouy[0]*app_accel_y[5]+0.3535533905932737*rhoux[0]*app_accel_x[5]+0.3535533905932737*app_accel_z[1]*rhouz[3]+0.3535533905932737*app_accel_y[1]*rhouy[3]+0.3535533905932737*app_accel_x[1]*rhoux[3]+0.3535533905932737*rhouz[1]*app_accel_z[3]+0.3535533905932737*rhouy[1]*app_accel_y[3]+0.3535533905932737*rhoux[1]*app_accel_x[3]; 
  out_energy[6] += 0.3535533905932737*app_accel_z[1]*rhouz[7]+0.3535533905932737*app_accel_y[1]*rhouy[7]+0.3535533905932737*app_accel_x[1]*rhoux[7]+0.3535533905932737*rhouz[1]*app_accel_z[7]+0.3535533905932737*rhouy[1]*app_accel_y[7]+0.3535533905932737*rhoux[1]*app_accel_x[7]+0.3535533905932737*app_accel_z[0]*rhouz[6]+0.3535533905932737*app_accel_y[0]*rhouy[6]+0.3535533905932737*app_accel_x[0]*rhoux[6]+0.3535533905932737*rhouz[0]*app_accel_z[6]+0.3535533905932737*rhouy[0]*app_accel_y[6]+0.3535533905932737*rhoux[0]*app_accel_x[6]+0.3535533905932737*app_accel_z[4]*rhouz[5]+0.3535533905932737*app_accel_y[4]*rhouy[5]+0.3535533905932737*app_accel_x[4]*rhoux[5]+0.3535533905932737*rhouz[4]*app_accel_z[5]+0.3535533905932737*rhouy[4]*app_accel_y[5]+0.3535533905932737*rhoux[4]*app_accel_x[5]+0.3535533905932737*app_accel_z[2]*rhouz[3]+0.3535533905932737*app_accel_y[2]*rhouy[3]+0.3535533905932737*app_accel_x[2]*rhoux[3]+0.3535533905932737*rhouz[2]*app_accel_z[3]+0.3535533905932737*rhouy[2]*app_accel_y[3]+0.3535533905932737*rhoux[2]*app_accel_x[3]; 
  out_energy[7] += 0.3535533905932737*app_accel_z[0]*rhouz[7]+0.3535533905932737*app_accel_y[0]*rhouy[7]+0.3535533905932737*app_accel_x[0]*rhoux[7]+0.3535533905932737*rhouz[0]*app_accel_z[7]+0.3535533905932737*rhouy[0]*app_accel_y[7]+0.3535533905932737*rhoux[0]*app_accel_x[7]+0.3535533905932737*app_accel_z[1]*rhouz[6]+0.3535533905932737*app_accel_y[1]*rhouy[6]+0.3535533905932737*app_accel_x[1]*rhoux[6]+0.3535533905932737*rhouz[1]*app_accel_z[6]+0.3535533905932737*rhouy[1]*app_accel_y[6]+0.3535533905932737*rhoux[1]*app_accel_x[6]+0.3535533905932737*app_accel_z[2]*rhouz[5]+0.3535533905932737*app_accel_y[2]*rhouy[5]+0.3535533905932737*app_accel_x[2]*rhoux[5]+0.3535533905932737*rhouz[2]*app_accel_z[5]+0.3535533905932737*rhouy[2]*app_accel_y[5]+0.3535533905932737*rhoux[2]*app_accel_x[5]+0.3535533905932737*app_accel_z[3]*rhouz[4]+0.3535533905932737*app_accel_y[3]*rhouy[4]+0.3535533905932737*app_accel_x[3]*rhoux[4]+0.3535533905932737*rhouz[3]*app_accel_z[4]+0.3535533905932737*rhouy[3]*app_accel_y[4]+0.3535533905932737*rhoux[3]*app_accel_x[4]; 

} 
